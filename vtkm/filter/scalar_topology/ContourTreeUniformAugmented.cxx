//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
// Copyright (c) 2018, The Regents of the University of California, through
// Lawrence Berkeley National Laboratory (subject to receipt of any required approvals
// from the U.S. Dept. of Energy).  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National
//     Laboratory, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without
//     specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.
//
//=============================================================================
//
//  This code is an extension of the algorithm presented in the paper:
//  Parallel Peak Pruning for Scalable SMP Contour Tree Computation.
//  Hamish Carr, Gunther Weber, Christopher Sewell, and James Ahrens.
//  Proceedings of the IEEE Symposium on Large Data Analysis and Visualization
//  (LDAV), October 2016, Baltimore, Maryland.
//
//  The PPP2 algorithm and software were jointly developed by
//  Hamish Carr (University of Leeds), Gunther H. Weber (LBNL), and
//  Oliver Ruebel (LBNL)
//==============================================================================

#include <vtkm/filter/scalar_topology/ContourTreeUniformAugmented.h>
#include <vtkm/filter/scalar_topology/internal/ComputeBlockIndices.h>
#include <vtkm/filter/scalar_topology/worklet/ContourTreeUniformAugmented.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/ContourTreeMesh.h>
// Adding the new TopologyGraph Class
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/TopologyGraph.h>

// clang-format off
VTKM_THIRDPARTY_PRE_INCLUDE
#include <vtkm/thirdparty/diy/Configure.h>
#include <vtkm/thirdparty/diy/diy.h>
VTKM_THIRDPARTY_POST_INCLUDE
// clang-format on

#include <vtkm/filter/scalar_topology/worklet/contourtree_distributed/ContourTreeBlockData.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_distributed/MergeBlockFunctor.h>

#include <memory>


// for memory usage
#include <sys/resource.h>
#include <unistd.h>

#define WRITE_FILES 0


// New includes for writing the VTK-m parallel adjacency builder
// (add here)
#include <vtkm/worklet/WorkletMapTopology.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/Algorithm.h>
//#include <vtkm/worklet/FieldMap.h>

// Output type: vector of sets or similar (proxy for per-point neighbor sets)
using AdjacencyList = vtkm::cont::ArrayHandle<std::set<vtkm::Id>>;

class BuildAdjacencyWorklet : public vtkm::worklet::WorkletVisitPointsWithCells
{
public:
  using ControlSignature = void(CellSetIn cellSet,
                                FieldOut adjacencyList);
  using ExecutionSignature = void(WorkIndex, CellIndices, _2);
  using InputDomain = _1;

    template<typename CellIndicesVecType, typename OutSetType>
    VTKM_EXEC void operator()(vtkm::Id pointId,
                              const CellIndicesVecType& cellIndices,
                              OutSetType& outputSet) const
    {
      std::set<vtkm::Id> neighbors;

      for (vtkm::IdComponent c = 0; c < cellIndices.GetNumberOfComponents(); ++c)
      {
        auto cell = cellIndices[c];
        for (vtkm::IdComponent j = 0; j < cell.GetNumberOfComponents(); ++j)
        {
          vtkm::Id neighborId = cell[j];
          if (neighborId != pointId) // ← now correctly using point ID
          {
            neighbors.insert(neighborId);
          }
        }
      }

      outputSet = neighbors;
    }
};



class BuildFlatAdjacencyWorklet : public vtkm::worklet::WorkletVisitPointsWithCells
{
public:
  using ControlSignature = void(CellSetIn cellSet,
                                WholeArrayOut neighbors,
                                WholeArrayOut counts);
  using ExecutionSignature = void(CellIndices, _2, _3, WorkIndex);
  using InputDomain = _1;

  template<typename CellIndicesVecType, typename NeighborsArrayType, typename CountArrayType>
  VTKM_EXEC void operator()(const CellIndicesVecType& cellIndices,
                            NeighborsArrayType& neighbors,
                            CountArrayType& counts,
                            /* WorkIndex */ vtkm::Id pointId) const
  {
    vtkm::Id neighborCount = 0;
    const vtkm::Id maxNeighbors = 32; // Set an upper bound for stack buffer
    vtkm::Id temp[maxNeighbors];
    vtkm::Id tempCount = 0;

    for (vtkm::IdComponent c = 0; c < cellIndices.GetNumberOfComponents(); ++c)
    {
      auto cell = cellIndices[c];
      for (vtkm::IdComponent j = 0; j < cell.GetNumberOfComponents(); ++j)
      {
        vtkm::Id neighborId = cell[j];
        if (neighborId != pointId)
        {
          bool duplicate = false;
          for (vtkm::IdComponent k = 0; k < tempCount; ++k)
          {
            if (temp[k] == neighborId)
            {
              duplicate = true;
              break;
            }
          }
          if (!duplicate && tempCount < maxNeighbors)
          {
            temp[tempCount++] = neighborId;
          }
        }
      }
    }

    counts.Set(pointId, tempCount);
    for (vtkm::IdComponent i = 0; i < tempCount; ++i)
    {
      neighbors.Set(pointId * maxNeighbors + i, temp[i]);
    }
  }
};

//class CountNeighbors : public vtkm::worklet::WorkletVisitCellsWithPoints
//{
//public:
//  using ControlSignature = void(CellSetIn, AtomicArrayInOut counts);
//  using ExecutionSignature = void(PointCount, PointIndices, _2);
//  using InputDomain = _1;

//  template <typename PointIndexVec, typename CountArray>
//  VTKM_EXEC void operator()(vtkm::IdComponent numPointsInCell,
//                            const PointIndexVec& pointIndices,
//                            CountArray& counts) const
//  {
//    for (vtkm::IdComponent i = 0; i < 4; ++i)
////        for (vtkm::IdComponent i = 0; i < numPointsInCell; ++i)
//    {
//      // Each point is connected to all others in the same cell, except itself
//      vtkm::Id pointId = pointIndices[i];
////    vtkm::AtomicAdd(&counts.GetPortal().Get(pointId), numPointsInCell - 1);
////      vtkm::AtomicAdd(&counts, pointId, static_cast<vtkm::Id>(numPointsInCell - 1));
////      counts.Set(pointId, numPointsInCell - 1);
////      counts.Add(pointId, numPointsInCell - 1);
//      counts.Add(pointId, 1);
//    }
//  }
//};

class CountNeighborsWorklet : public vtkm::worklet::WorkletVisitPointsWithCells
{
public:
  using ControlSignature = void(CellSetIn cellSet, WholeArrayOut count);
  using ExecutionSignature = void(CellIndices, _2, WorkIndex);
  using InputDomain = _1;

  template <typename CellIndicesVecType, typename CountPortalType>
  VTKM_EXEC void operator()(const CellIndicesVecType& cellIndices,
                            CountPortalType& count,
                            vtkm::Id pointId) const
  {
    // Use a local fixed-size buffer to deduplicate small numbers of neighbors
    constexpr vtkm::IdComponent maxLocal = 64;
    vtkm::Id local[maxLocal];
    vtkm::IdComponent localCount = 0;

    for (vtkm::IdComponent i = 0; i < cellIndices.GetNumberOfComponents(); ++i)
    {
      auto cell = cellIndices[i];
      for (vtkm::IdComponent j = 0; j < cell.GetNumberOfComponents(); ++j)
      {
        vtkm::Id n = cell[j];
        if (n != pointId)
        {
          // Check for duplicates
          bool duplicate = false;
          for (vtkm::IdComponent k = 0; k < localCount; ++k)
          {
            if (local[k] == n)
            {
              duplicate = true;
              break;
            }
          }
          if (!duplicate && localCount < maxLocal)
          {
            local[localCount++] = n;
          }
        }
      }
    }

    count.Set(pointId, localCount);
  }
};


class FillNeighborsWorklet : public vtkm::worklet::WorkletVisitPointsWithCells
{
public:
  using ControlSignature = void(CellSetIn cellSet,
                                WholeArrayIn offsets,
                                WholeArrayOut neighbors);
  using ExecutionSignature = void(CellIndices, _2, _3, WorkIndex);
  using InputDomain = _1;

  template <typename CellIndicesVecType,
            typename OffsetsPortalType,
            typename NeighborPortalType>
  VTKM_EXEC void operator()(const CellIndicesVecType& cellIndices,
                            const OffsetsPortalType& offsets,
                            NeighborPortalType& flatNeighbors,
                            vtkm::Id pointId) const
  {
    constexpr vtkm::IdComponent maxLocal = 64;
    vtkm::Id local[maxLocal];
    vtkm::IdComponent localCount = 0;

    for (vtkm::IdComponent i = 0; i < cellIndices.GetNumberOfComponents(); ++i)
    {
      auto cell = cellIndices[i];
      for (vtkm::IdComponent j = 0; j < cell.GetNumberOfComponents(); ++j)
      {
        vtkm::Id n = cell[j];
        if (n != pointId)
        {
          bool duplicate = false;
          for (vtkm::IdComponent k = 0; k < localCount; ++k)
          {
            if (local[k] == n)
            {
              duplicate = true;
              break;
            }
          }
          if (!duplicate && localCount < maxLocal)
          {
            local[localCount++] = n;
          }
        }
      }
    }

    vtkm::Id offset = offsets.Get(pointId);
    for (vtkm::IdComponent i = 0; i < localCount; ++i)
    {
      flatNeighbors.Set(offset + i, local[i]);
    }
  }
};



class EmitEdges : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn cellId,
                                WholeArrayIn conn,
                                WholeArrayIn offsets,
                                WholeArrayOut edgeStarts,
                                WholeArrayOut edgeEnds);
  using ExecutionSignature = void(_1, _2, _3, _4, _5);

  template <typename ConnPortal, typename OffsetPortal, typename OutPortal>
  VTKM_EXEC void operator()(vtkm::Id cellId,
                            const ConnPortal& conn,
                            const OffsetPortal& offsets,
                            OutPortal& edgeStarts,
                            OutPortal& edgeEnds) const
  {
    const vtkm::Id start = offsets.Get(cellId);
    const vtkm::Id v[4] = {
      conn.Get(start),
      conn.Get(start + 1),
      conn.Get(start + 2),
      conn.Get(start + 3)
    };

    // Emit all 6 undirected edges
    const vtkm::Id edgePairs[6][2] = {
      {0, 1}, {0, 2}, {0, 3},
      {1, 2}, {1, 3}, {2, 3}
    };

    for (vtkm::Id i = 0; i < 6; ++i)
    {
      vtkm::Id u = v[edgePairs[i][0]];
      vtkm::Id w = v[edgePairs[i][1]];
      if (u > w)
      {
          vtkm::Id tmp = u;
          u = w;
          w = tmp;
      }

      vtkm::Id outIdx = cellId * 6 + i;
      edgeStarts.Set(outIdx, u);
      edgeEnds.Set(outIdx, w);
    }
  }
};



namespace vtkm
{
namespace filter
{
namespace scalar_topology
{

//-----------------------------------------------------------------------------
ContourTreeAugmented::ContourTreeAugmented(bool useMarchingCubes,
                                           unsigned int computeRegularStructure)
  : UseMarchingCubes(useMarchingCubes)
  , ComputeRegularStructure(computeRegularStructure)
  , MultiBlockTreeHelper(nullptr)
{
  this->SetOutputFieldName("resultData");
}

void ContourTreeAugmented::SetBlockIndices(
  vtkm::Id3 blocksPerDim,
  const vtkm::cont::ArrayHandle<vtkm::Id3>& localBlockIndices)
{
  if (this->MultiBlockTreeHelper)
  {
    this->MultiBlockTreeHelper.reset();
  }
  this->MultiBlockTreeHelper =
    std::make_unique<vtkm::worklet::contourtree_distributed::MultiBlockContourTreeHelper>(
      blocksPerDim, localBlockIndices);
}

const vtkm::worklet::contourtree_augmented::ContourTree& ContourTreeAugmented::GetContourTree()
  const
{
  return this->ContourTreeData;
}

const vtkm::worklet::contourtree_augmented::IdArrayType& ContourTreeAugmented::GetSortOrder() const
{
  return this->MeshSortOrder;
}

vtkm::Id ContourTreeAugmented::GetNumIterations() const
{
  return this->NumIterations;
}

// Returns current memory usage in KB
size_t getCurrentRSS_linux()
{
    std::ifstream status_file("/proc/self/status");
    std::string line;

    while(std::getline(status_file, line))
    {
        if(line.find("VmRSS:") == 0)
        {
            std::istringstream iss(line);
            std::string key;
            size_t memory; // memory value in kB
            std::string unit;

            iss >> key >> memory >> unit;
            return memory; // return in KB
        }
    }
    return 0; // Not found
}


void static printMemoryUsage(const std::string& message)
{
    // Red text formatting for highlighting some console output:
    const std::string ORANGE = "\033[38;2;255;165;0m";  // Start red text
    const std::string LIGHT_BLUE = "\033[38;5;117m";  // Light blue in 256-color
    const std::string RESET = "\033[0m"; // End red text

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

//    std::cout << ORANGE << message << LIGHT_BLUE << " - Memory usage (peak): " << usage.ru_maxrss
//              << " KB | (current) " << getCurrentRSS_linux() << " KB" << RESET << std::endl;
    std::cout << LIGHT_BLUE << message << " - Memory usage (peak): " << usage.ru_maxrss
              << " KB | (current) " << getCurrentRSS_linux() << " KB" << RESET << std::endl;
}

// Define type aliases for clarity:
//using AdjacencyList = std::unordered_map<vtkm::Id, std::unordered_set<vtkm::Id>>;
using AdjacencyList = std::unordered_map<vtkm::Id, std::set<vtkm::Id>>;

// ---------------------------------------------------------------------------
// Assumes pure tetrahedral mesh, no offset array, always 4 indices per cell:
// connectivity: [p0,p1,p2,p3, p4,p5,p6,p7, ...] (groups of 4 points per cell)
// numCells: total number of tetrahedral cells
// ---------------------------------------------------------------------------
AdjacencyList MakeAdjacencyTetrahedron(
    const vtkm::cont::ArrayHandle<vtkm::Id> &connectivity)
{
    AdjacencyList adjacency;

    auto connPortal = connectivity.ReadPortal();
    vtkm::Id numValues = connectivity.GetNumberOfValues();

    for (vtkm::Id idx = 0; idx < numValues; idx += 4)
    {
        vtkm::Id vertices[4] = {
            connPortal.Get(idx),
            connPortal.Get(idx+1),
            connPortal.Get(idx+2),
            connPortal.Get(idx+3)
        };

        // Loop through pairs of points in tetrahedron
        for (int i = 0; i < 4; ++i)
        {
            vtkm::Id vi = vertices[i];
            auto &neighbors = adjacency[vi]; // auto create if not exist

            for (int j = 0; j < 4; ++j)
            {
                if (i == j) continue;
                neighbors.insert(vertices[j]);
            }
        }
    }

    return adjacency;
}



// -----------------------------------------------------------------------------------
// General version using offsets array. Can handle arbitrary cell sizes/ connectivity.
// connectivity: Flat list of vertices
// offsets: Defines start of each cell's connectivity indices (offsets[i+1]-offsets[i])
// -----------------------------------------------------------------------------------
AdjacencyList MakeAdjacencyWithOffsets(
    const vtkm::cont::ArrayHandle<vtkm::Id> &connectivity,
    const vtkm::cont::ArrayHandle<vtkm::Id, vtkm::cont::StorageTagCounting> &offsets)
{
    AdjacencyList adjacency;

    auto connPortal = connectivity.ReadPortal();
    auto offsetPortal = offsets.ReadPortal();

    vtkm::Id numCells = offsets.GetNumberOfValues() - 1;

    for (vtkm::Id cellId = 0; cellId < numCells; ++cellId)
    {
        vtkm::Id start = offsetPortal.Get(cellId);
        vtkm::Id end = offsetPortal.Get(cellId+1);

        // Iterate through vertices in the cell:
        for (vtkm::Id i = start; i < end; ++i)
        {
            vtkm::Id vi = connPortal.Get(i);
            auto &neighbors = adjacency[vi];

            for (vtkm::Id j = start; j < end; ++j)
            {
                vtkm::Id vj = connPortal.Get(j);
                if(vi == vj) continue;
                neighbors.insert(vj);
            }
        }
    }

    return adjacency;
}


//-----------------------------------------------------------------------------
vtkm::cont::DataSet ContourTreeAugmented::DoExecute(const vtkm::cont::DataSet& input)
{
  std::cout << "{sc_tp/ContourTreeUniformAugmented.cxx - ContourTreeAugmented::DoExecute() Version for Point Clouds / Irregular Meshes\n";
  // Global benchmark timer for the whole ContourTreeAugmented::DoExecute() function
  vtkm::cont::Timer timer;
  timer.Start();

  // Local profiling timers to find the bottlenecks:
  vtkm::cont::Timer profiling1;
  vtkm::cont::Timer profiling2;
  vtkm::cont::Timer profiling3;
  vtkm::cont::Timer profiling4;
  vtkm::cont::Timer profiling5;
  vtkm::cont::Timer profiling6;
  profiling1.Start();

  // Input dataset is expected as a tetrahedral cell set (single cell type), ...
  // ... and irregular (no dimensions)

  // TODO blockIndex needs to change if we have multiple blocks per MPI rank and DoExecute is called for multiple blocks
  std::size_t blockIndex = 0;

  // Determine if and what augmentation we need to do
  unsigned int compRegularStruct = this->ComputeRegularStructure;
  // When running in parallel we need to at least augment with the boundary vertices
  if (compRegularStruct == 0)
  {
    if (this->MultiBlockTreeHelper)
    {
      if (this->MultiBlockTreeHelper->GetGlobalNumberOfBlocks() > 1)
      {
        compRegularStruct = 2; // Compute boundary augmentation
      }
    }
  }

  // Create the result object
  vtkm::cont::DataSet result;

  // FIXME: reduce the size of lambda.
  auto resolveType = [&](const auto& concrete)
  {// lambda start "auto resolveType = [&](const auto& concrete)"
    using T = typename std::decay_t<decltype(concrete)>::ValueType;

    // PACTBD-EDIT-FIXED
    // Tetrahedral connections used to be read from a file here:
    //      const std::string filename_vtk = "/home/sc17dd/modules/HCTC2024/VTK-m-topology-refactor/VTK-m_TopologyGraph-Refactor/examples/contour-visualiser/delaunay-parcels/200k-from-2M-sampled-excel-sorted.1-withvalues-manual.vtk";
    // ... the .vtk input is now read from the main applet (such as examples/ContourTreeApp.cxx) ...
    // ... and accessed from the parameter 'input' in DoExecute(const vtkm::cont::DataSet& input)

    int num_datapoints;
    vtkm::cont::ArrayHandle<vtkm::Id> nbor_connectivity_auto;
    vtkm::cont::ArrayHandle<vtkm::Id> nbor_offsets_auto;

    vtkm::cont::ArrayHandle<vtkm::Id> flatNeighbors;
    vtkm::cont::ArrayHandle<vtkm::Id> neighborOffsets;
    // Obtain portals to write data directly:

    // (the scoping deletes the reader right after populating the cont::DataSet)
    {
      num_datapoints = input.GetPointField("var").GetNumberOfValues();

      // Explicitly interpret as tetrahedral cell set
      // (already checking in the main applet)
      using TetCellSet = vtkm::cont::CellSetSingleType<>;
      printMemoryUsage("[ContourTreeUniformAugmented.cxx] BEFORE Creating Adjacency List");
      std::cout << "\t\t (1) Time-to-here: " << profiling1.GetElapsedTime() << " seconds" << std::endl;
      profiling2.Start();

////      /// START: Testing the new parallel worklet:

      vtkm::cont::CellSetSingleType<> cellSet =
          input.GetCellSet().AsCellSet<vtkm::cont::CellSetSingleType<>>();

//      vtkm::cont::Invoker invoke;
//      vtkm::cont::ArrayHandle<std::set<vtkm::Id>> adjacencyList;
//      invoke(BuildAdjacencyWorklet{}, cellSet, adjacencyList);

      vtkm::Id numPoints = num_datapoints;

      vtkm::Id numCells = cellSet.GetNumberOfCells();
      vtkm::Id numEdges = numCells * 6;

      // Step 1: Get connectivity and offsets
      auto conn = cellSet.GetConnectivityArray(
                    vtkm::TopologyElementTagCell{},
                    vtkm::TopologyElementTagPoint{});
      auto offsets = cellSet.GetOffsetsArray(
                       vtkm::TopologyElementTagCell{},
                       vtkm::TopologyElementTagPoint{});

      // Step 2: Prepare edge output arrays
      vtkm::cont::ArrayHandle<vtkm::Id> edgeStarts, edgeEnds;
      edgeStarts.Allocate(numEdges);
      edgeEnds.Allocate(numEdges);

      // Step 3: Generate cell IDs [0, 1, ..., numCells-1]
      vtkm::cont::ArrayHandleCounting<vtkm::Id> cellIds(0, 1, numCells);

      // Step 4: Dispatch EmitEdges
      vtkm::worklet::DispatcherMapField<EmitEdges> emitDispatcher;
      emitDispatcher.Invoke(cellIds, conn, offsets, edgeStarts, edgeEnds);

      // Step 5: Zip, sort, and deduplicate
      auto edgePairs = vtkm::cont::make_ArrayHandleZip(edgeStarts, edgeEnds);
      vtkm::cont::Algorithm::Sort(edgePairs);
      vtkm::cont::Algorithm::Unique(edgePairs);


      using EdgePairArray = vtkm::cont::ArrayHandleZip<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                       vtkm::cont::ArrayHandle<vtkm::Id>>;

      // Step 1: For each edge (u, v), emit two entries: (u → v) and (v → u)
      vtkm::cont::ArrayHandle<vtkm::Id> srcs, dsts;
      {
        auto numEdges = edgePairs.GetNumberOfValues();
        srcs.Allocate(numEdges * 2);
        dsts.Allocate(numEdges * 2);

        auto edgePortal = edgePairs.ReadPortal();
        auto srcPortal = srcs.WritePortal();
        auto dstPortal = dsts.WritePortal();

        for (vtkm::Id i = 0; i < numEdges; ++i)
        {
          auto edge = edgePortal.Get(i);
          vtkm::Id u = edge.first;
          vtkm::Id v = edge.second;

          srcPortal.Set(2 * i,     u);  dstPortal.Set(2 * i,     v);
          srcPortal.Set(2 * i + 1, v);  dstPortal.Set(2 * i + 1, u);
        }
      }

      // Step 2: Sort by source vertex (to group neighbors together)
      auto neighborPairs = vtkm::cont::make_ArrayHandleZip(srcs, dsts);
      vtkm::cont::Algorithm::SortByKey(srcs, dsts); // now dsts are grouped by srcs

//      // Step 3: Build CSR offsets
      vtkm::Id numVertices = numPoints + 1; // maxVertexId from your mesh
//      vtkm::cont::ArrayHandle<vtkm::Id> neighborOffsets;
//      vtkm::cont::Algorithm::UpperBounds(srcs, vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, numVertices), neighborOffsets);
//      // Step 1: Count neighbors per vertex
//      vtkm::cont::ArrayHandle<vtkm::Id> neighborCounts;
//      vtkm::cont::Algorithm::UpperBounds(
//          srcs,
//          vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, numVertices),
//          neighborCounts);

//      // Step 2: Convert counts to offsets with exclusive scan
//      vtkm::cont::ArrayHandle<vtkm::Id> neighborOffsets;
//      vtkm::cont::Algorithm::ScanExclusive(neighborCounts, neighborOffsets);

//      vtkm::cont::ArrayHandle<vtkm::Id> neighborOffsets;
      vtkm::cont::Algorithm::UpperBounds(
          srcs,
          vtkm::cont::ArrayHandleCounting<vtkm::Id>(-1, 1, numVertices),
//                  vtkm::cont::ArrayHandleCounting<vtkm::Id>(-1, 1, numVertices + 1), old, didn't start from 0
          neighborOffsets);

      // Step 4: The dsts array is your flat neighbor list
//      vtkm::cont::ArrayHandle<vtkm::Id>
              flatNeighbors = dsts;

      std::cout << "flatNeighbors&neighborOffsets IN PARALLEL! :" << std::endl;
      vtkm::cont::printSummary_ArrayHandle(flatNeighbors, std::cout);
      vtkm::cont::printSummary_ArrayHandle(neighborOffsets, std::cout);


      // Now edgePairs contains all unique undirected edges (u < v)


////      /// END: Testing the new parallel worklet:

//      // From the input file, which we get from TetGen, we get a tet-mesh ...
//      // ... and the .vtk format gives us the connectivity and offset arrays.
//      // For each tet, the four vertex indices are given that make up the tet
//      // However, for the topology graph mesh abstraction ...
//      // ... we later require an adjacency list, namely:
//      // vertex_id[i] = vertex_neighbourhood_list
//      // The following takes connectivity and offset arrays and computes an adjacency list


////      std::unordered_map<vtkm::Id, std::set<vtkm::Id>>adjacency_list = MakeAdjacencyWithOffsets(/* connectivity */
////                         input.GetCellSet().AsCellSet<TetCellSet>().GetConnectivityArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}),
////                                                                                                /* offsets */
////                         input.GetCellSet().AsCellSet<TetCellSet>().GetOffsetsArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}));



//      std::unordered_map<vtkm::Id, std::set<vtkm::Id>> adjacency_list;

//      auto connPortal = input.GetCellSet().AsCellSet<TetCellSet>().GetConnectivityArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}).ReadPortal();
//      auto offsetPortal = input.GetCellSet().AsCellSet<TetCellSet>().GetOffsetsArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}).ReadPortal();

////      vtkm::Id numCells = input.GetCellSet().AsCellSet<TetCellSet>().GetOffsetsArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}).GetNumberOfValues() - 1;

//      for (vtkm::Id cellId = 0; cellId < numCells; ++cellId)
//      {
//          vtkm::Id start = offsetPortal.Get(cellId);
//          vtkm::Id end = offsetPortal.Get(cellId+1);

//          // Iterate through vertices in the cell:
//          for (vtkm::Id i = start; i < end; ++i)
//          {
//              vtkm::Id vi = connPortal.Get(i);
//              auto &neighbors = adjacency_list[vi];

//              for (vtkm::Id j = start; j < end; ++j)
//              {
//                  vtkm::Id vj = connPortal.Get(j);
//                  if(vi == vj) continue;
//                  neighbors.insert(vj);
//              }
//          }
//      }




//      std::cout << "From data nbor{_connectivity/_offsets_}auto:" << std::endl;
//      vtkm::cont::printSummary_ArrayHandle(input.GetCellSet().AsCellSet<TetCellSet>().GetConnectivityArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}),
//                                           std::cout);
//      vtkm::cont::printSummary_ArrayHandle(input.GetCellSet().AsCellSet<TetCellSet>().GetOffsetsArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}),
//                                           std::cout);









//      // std::unordered_map<vtkm::Id, std::set<vtkm::Id>>adjacency_list = MakeAdjacencyTetrahedron(input.GetCellSet().AsCellSet<TetCellSet>().GetConnectivityArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}));

//      printMemoryUsage("[ContourTreeUniformAugmented.cxx] AFTER MakeAdjacencyTetrahedron");
      std::cout << "\t\t (2) Time-to-here: " << profiling2.GetElapsedTime() << " seconds" << std::endl;
      profiling3.Start();

//      // First, compute total sizes:
//      vtkm::Id total_nbors = 0;
//      for (vtkm::Id i = 0; i < num_datapoints; ++i)
//      {
//          total_nbors += adjacency_list[i].size();
////          total_nbors += adjacencyList[i].size();
//      }

//      nbor_connectivity_auto.Allocate(total_nbors);
//      // offsets array has size num_points + 1 because ...
//      // ... it includes the offset after the last position
//      nbor_offsets_auto.Allocate(num_datapoints + 1);

//      auto connectivityPortal = nbor_connectivity_auto.WritePortal();
//      auto offsetsPortal = nbor_offsets_auto.WritePortal();

//      // Populate the data
//      vtkm::Id offset_counter = 0;
//      offsetsPortal.Set(0, offset_counter); // initial offset is always 0

//      for (vtkm::Id i = 0; i < num_datapoints; ++i)
//      {
//          for (auto elem : adjacency_list[i])
////          for (auto elem : adjacencyList[i])
//          {
//              connectivityPortal.Set(offset_counter++, elem);
//          }
//          // offset for next datapoint starts here
//          offsetsPortal.Set(i + 1, offset_counter);
//      }

//      std::cout << "Real nbor{_connectivity/_offsets_}auto:" << std::endl;
//      vtkm::cont::printSummary_ArrayHandle(nbor_connectivity_auto, std::cout);
//      vtkm::cont::printSummary_ArrayHandle(nbor_offsets_auto, std::cout);

    }





    printMemoryUsage("[ContourTreeUniformAugmented.cxx] AFTER 'connectivity & offsets'");
    std::cout << "\t\t (3) Time-to-here: " << profiling3.GetElapsedTime() << " seconds" << std::endl;
    profiling4.Start();

      // populate 'freeby' arrays:
      // nodes_sorted - just a list of the nodes going from 0-to-N, ...
      // ... where N = num_datapoints
    // ------------------------------------------------------------------------------- //
      std::vector<vtkm::Id> std_nodes_sorted;
      for(vtkm::Id i = 0; i < num_datapoints; i++)
      {
        std_nodes_sorted.push_back(i);
      }

      vtkm::cont::ArrayHandle<vtkm::Id> nodes_sorted =
        vtkm::cont::make_ArrayHandle(std_nodes_sorted, vtkm::CopyFlag::Off);
    // ------------------------------------------------------------------------------- //
      std::vector<int> std_actual_values;
      for(int i = 0; i < num_datapoints; i++)
      {
        std_actual_values.push_back(i);
      }
     // note: this array does not actually matter for our case ...
      vtkm::cont::ArrayHandle<int> actual_values =
        vtkm::cont::make_ArrayHandle(std_actual_values, vtkm::CopyFlag::Off);
      // ------------------------------------------------------------------------------- //

      std::vector<vtkm::Id> std_global_inds = {0};
      vtkm::cont::ArrayHandle<vtkm::Id> global_inds =
        vtkm::cont::make_ArrayHandle(std_global_inds, vtkm::CopyFlag::Off);
      // ------------------------------------------------------------------------------- //

      // VTK FILE

#if WRITE_FILES
      std::ofstream file("ContourTreeGraph--recreate-CONNECTIONS.txt");

      vtkm::Id globalNbrID = 0;

//      std::cout << "nbor_connectivity.GetNumberOfValues() " << nbor_connectivity.GetNumberOfValues() << std::endl;
//      std::cout << "nbor_offsets.GetNumberOfValues() " << nbor_offsets.GetNumberOfValues() << std::endl;

      for (vtkm::Id vertexID = 1; vertexID <= num_datapoints; vertexID++) // skip 0 because offsets start at 0
      {
//        for (globalNbrID; globalNbrID < nbor_offsets.ReadPortal().Get(vertexID); globalNbrID++)
          for (globalNbrID; globalNbrID < nbor_offsets_auto.ReadPortal().Get(vertexID); globalNbrID++)
          {
//            file << nbor_connectivity.ReadPortal().Get(globalNbrID);
              file << nbor_connectivity_auto.ReadPortal().Get(globalNbrID);

//            if(globalNbrID+1 != nbor_offsets.ReadPortal().Get(vertexID))
              if(globalNbrID+1 != nbor_offsets_auto.ReadPortal().Get(vertexID))
              {
                  file << " ";
              }
          }
          file << std::endl;
      }
      file.close();
#endif

//       Using the 'TopologyGraph' constructor here ...
//       ... (to be moved to a separate class eventually)
      std::cout << "[ContourTreeUniformAugmented.cxx] Constructing a 'ContourTreeMesh' (future 'TopologyGraph')\n";
      vtkm::worklet::contourtree_augmented::ContourTreeMesh<int> mesh(nodes_sorted,
                              //arcs_list,
                                flatNeighbors, //nbor_connectivity_auto, //nbor_connectivity,
                                neighborOffsets, //nbor_offsets_auto,//nbor_offsets,
                                nodes_sorted,
                                // doesnt work out of the box:
                                // fieldArray, // testing fieldArray instead of manual actual_value
                                actual_values,
                                //nodes_sorted,
                                global_inds);
      std::cout << "[ContourTreeUniformAugmented.cxx] 'ContourTreeMesh' Constructed\n";

      std::cout << "\t\t (4) Time-to-here: " << profiling4.GetElapsedTime() << " seconds" << std::endl;
      profiling5.Start();

    // NOTE: The TopologyGraph ContourTree worklet no longer uses meshSize and this->useMarchingCubes
    // (because the mesh is irregular and does not have inherent dimensions/connectivity)
    vtkm::worklet::ContourTreeAugmented worklet;
    // Using the new TopologyGraph ContourTree worklet:
    worklet.Run(concrete,                                                                   // fieldArray?
                mesh,
                MultiBlockTreeHelper ? MultiBlockTreeHelper->LocalContourTrees[blockIndex]  //
                                     : this->ContourTreeData,                               // contourTree
                MultiBlockTreeHelper ? MultiBlockTreeHelper->LocalSortOrders[blockIndex]    //
                                     : this->MeshSortOrder,                                 // sortOrder
                this->NumIterations,                                                        // nIterations
                compRegularStruct);                                                         // computeRegularStructure

    std::cout << "[ContourTreeUniformAugmented.cxx] 'ContourTreeAugmented worklet' Finished\n";
    std::cout << "\t\t (5) Time-to-here: " << profiling5.GetElapsedTime() << " seconds" << std::endl;
    profiling6.Start();

    // If we run in parallel but with only one global block, then we need set our outputs correctly
    // here to match the expected behavior in parallel
    if (this->MultiBlockTreeHelper)
    {
      if (this->MultiBlockTreeHelper->GetGlobalNumberOfBlocks() == 1)
      {
        // Copy the contour tree and mesh sort order to the output
        this->ContourTreeData = this->MultiBlockTreeHelper->LocalContourTrees[0];
        this->MeshSortOrder = this->MultiBlockTreeHelper->LocalSortOrders[0];
        // In parallel we need the sorted values as output resulti
        // Construct the sorted values by permutting the input field
        auto fieldPermutted =
          vtkm::cont::make_ArrayHandlePermutation(this->MeshSortOrder, concrete);
        // FIXME: can sortedValues be ArrayHandleUnknown?
        vtkm::cont::ArrayHandle<T> sortedValues;
        vtkm::cont::Algorithm::Copy(fieldPermutted, sortedValues);

        // FIXME: is this the right way to create the DataSet? The original code creates an empty
        //  DataSet without any coordinate system etc.
        result = this->CreateResultField(input,
                                         this->GetOutputFieldName(),
                                         vtkm::cont::Field::Association::WholeDataSet,
                                         sortedValues);
        //        vtkm::cont::Field rfield(
        //          this->GetOutputFieldName(), vtkm::cont::Field::Association::WholeDataSet, sortedValues);
        //        result.AddField(rfield);
        //        return result;
      }
    }
    else
    {
      // Construct the expected result for serial execution. Note, in serial the result currently
      // not actually being used, but in parallel we need the sorted mesh values as output
      // This part is being hit when we run in serial or parallel with more than one rank.
      result =
        this->CreateResultFieldPoint(input, this->GetOutputFieldName(), ContourTreeData.Arcs);
      //  return CreateResultFieldPoint(input, ContourTreeData.Arcs, this->GetOutputFieldName());
    }
  }; // lambda end "auto resolveType = [&](const auto& concrete)"


  const auto& field = this->GetFieldFromDataSet(input); // Added new 2025-04-16
  this->CastAndCallScalarField(field, resolveType); // call the above lambda function

  printMemoryUsage("[ContourTreeUniformAugmented.cxx] AFTER resolveType lambda scope");
  std::cout << "\t\t (6) Time-to-here: " << profiling6.GetElapsedTime() << " seconds" << std::endl;

  VTKM_LOG_S(vtkm::cont::LogLevel::Warn,//vtkm::cont::LogLevel::Perf,
             std::endl
               << "    " << std::setw(38) << std::left << "Contour Tree Filter DoExecute"
               << ": " << timer.GetElapsedTime() << " seconds");
  /// DEBUG PRINT std::cout << "sc_tp/ContourTreeUniformAugmented.cxx : ContourTreeAugmented::DoExecute() finished}\n";
  return result;
} // ContourTreeAugmented::DoExecute



























////-------------------------- ORIGINAL(-ISH) ------------------------------------//
//vtkm::cont::DataSet ContourTreeAugmented::DoExecute(const vtkm::cont::DataSet& input)
//{
//  std::cout << "{sc_tp/ContourTreeUniformAugmented.cxx - ContourTreeAugmented::DoExecute()\n";
//  vtkm::cont::Timer timer;
//  timer.Start();

//  // Check that the field is Ok
//  const auto& field = this->GetFieldFromDataSet(input);
//  if (!field.IsPointField())
//  {
//    throw vtkm::cont::ErrorFilterExecution("Point field expected.");
//  }

//  // Use the GetPointDimensions struct defined in the header to collect the meshSize information
//  vtkm::Id3 meshSize;
//  const auto& cells = input.GetCellSet();
//  cells.CastAndCallForTypes<VTKM_DEFAULT_CELL_SET_LIST_STRUCTURED>(
//    vtkm::worklet::contourtree_augmented::GetPointDimensions(), meshSize);

//  std::cout << "mesh size: " << meshSize[0] << "x" << meshSize[1] << "x" << meshSize[2] << std::endl;

//  // TODO blockIndex needs to change if we have multiple blocks per MPI rank and DoExecute is called for multiple blocks
//  std::size_t blockIndex = 0;

//  // Determine if and what augmentation we need to do
//  unsigned int compRegularStruct = this->ComputeRegularStructure;
//  // When running in parallel we need to at least augment with the boundary vertices
//  if (compRegularStruct == 0)
//  {
//    if (this->MultiBlockTreeHelper)
//    {
//      if (this->MultiBlockTreeHelper->GetGlobalNumberOfBlocks() > 1)
//      {
//        compRegularStruct = 2; // Compute boundary augmentation
//      }
//    }
//  }

//  // Create the result object
//  vtkm::cont::DataSet result;

//  // FIXME: reduce the size of lambda.
//  auto resolveType = [&](const auto& concrete) {
//    using T = typename std::decay_t<decltype(concrete)>::ValueType;

//    vtkm::worklet::ContourTreeAugmented worklet;
//    // Run the worklet
//    worklet.Run(concrete,                                                                   // fieldArray?
//                MultiBlockTreeHelper ? MultiBlockTreeHelper->LocalContourTrees[blockIndex]  //
//                                     : this->ContourTreeData,                               // contourTree
//                MultiBlockTreeHelper ? MultiBlockTreeHelper->LocalSortOrders[blockIndex]    //
//                                     : this->MeshSortOrder,                                 // sortOrder
//                this->NumIterations,                                                        // nIterations
//                meshSize,                                                                   // meshSize
//                this->UseMarchingCubes,                                                     // useMarchingCubes
//                compRegularStruct);                                                         // computeRegularStructure

//    // If we run in parallel but with only one global block, then we need set our outputs correctly
//    // here to match the expected behavior in parallel
//    if (this->MultiBlockTreeHelper)
//    {
//      if (this->MultiBlockTreeHelper->GetGlobalNumberOfBlocks() == 1)
//      {
//        // Copy the contour tree and mesh sort order to the output
//        this->ContourTreeData = this->MultiBlockTreeHelper->LocalContourTrees[0];
//        this->MeshSortOrder = this->MultiBlockTreeHelper->LocalSortOrders[0];
//        // In parallel we need the sorted values as output resulti
//        // Construct the sorted values by permutting the input field
//        auto fieldPermutted =
//          vtkm::cont::make_ArrayHandlePermutation(this->MeshSortOrder, concrete);
//        // FIXME: can sortedValues be ArrayHandleUnknown?
//        vtkm::cont::ArrayHandle<T> sortedValues;
//        vtkm::cont::Algorithm::Copy(fieldPermutted, sortedValues);

//        // FIXME: is this the right way to create the DataSet? The original code creates an empty
//        //  DataSet without any coordinate system etc.
//        result = this->CreateResultField(input,
//                                         this->GetOutputFieldName(),
//                                         vtkm::cont::Field::Association::WholeDataSet,
//                                         sortedValues);
//        //        vtkm::cont::Field rfield(
//        //          this->GetOutputFieldName(), vtkm::cont::Field::Association::WholeDataSet, sortedValues);
//        //        result.AddField(rfield);
//        //        return result;
//      }
//    }
//    else
//    {
//      // Construct the expected result for serial execution. Note, in serial the result currently
//      // not actually being used, but in parallel we need the sorted mesh values as output
//      // This part is being hit when we run in serial or parallel with more than one rank.
//      result =
//        this->CreateResultFieldPoint(input, this->GetOutputFieldName(), ContourTreeData.Arcs);
//      //  return CreateResultFieldPoint(input, ContourTreeData.Arcs, this->GetOutputFieldName());
//    }
//  };
//  this->CastAndCallScalarField(field, resolveType);

//  VTKM_LOG_S(vtkm::cont::LogLevel::Perf,
//             std::endl
//               << "    " << std::setw(38) << std::left << "Contour Tree Filter DoExecute"
//               << ": " << timer.GetElapsedTime() << " seconds");
//  /// DEBUG PRINT std::cout << "sc_tp/ContourTreeUniformAugmented.cxx : ContourTreeAugmented::DoExecute() finished}\n";
//  return result;
//} // ContourTreeAugmented::DoExecute





























// TODO: is multiblock case ever tested?
VTKM_CONT vtkm::cont::PartitionedDataSet ContourTreeAugmented::DoExecutePartitions(
  const vtkm::cont::PartitionedDataSet& input)
{
  this->PreExecute(input);
  auto result = this->Filter::DoExecutePartitions(input);
  this->PostExecute(input, result);
  return result;
}

//-----------------------------------------------------------------------------
VTKM_CONT void ContourTreeAugmented::PreExecute(const vtkm::cont::PartitionedDataSet& input)
{
  if (this->MultiBlockTreeHelper)
  {
    if (input.GetGlobalNumberOfPartitions() !=
        this->MultiBlockTreeHelper->GetGlobalNumberOfBlocks())
    {
      throw vtkm::cont::ErrorFilterExecution(
        "Global number of block in MultiBlock dataset does not match the SpatialDecomposition");
    }
    if (this->MultiBlockTreeHelper->GetLocalNumberOfBlocks() != input.GetNumberOfPartitions())
    {
      throw vtkm::cont::ErrorFilterExecution(
        "Global number of block in MultiBlock dataset does not match the SpatialDecomposition");
    }
  }
  else
  {
    // No block indices set -> compute information automatically later
    this->MultiBlockTreeHelper =
      std::make_unique<vtkm::worklet::contourtree_distributed::MultiBlockContourTreeHelper>(input);
  }
}

//-----------------------------------------------------------------------------
template <typename T>
VTKM_CONT void ContourTreeAugmented::DoPostExecute(const vtkm::cont::PartitionedDataSet& input,
                                                   vtkm::cont::PartitionedDataSet& output)
{
  std::cout << "DoPostExecute()\n";
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  vtkm::Id size = comm.size();
  vtkm::Id rank = comm.rank();

  std::vector<vtkm::worklet::contourtree_augmented::ContourTreeMesh<T>*> localContourTreeMeshes;
  localContourTreeMeshes.resize(static_cast<std::size_t>(input.GetNumberOfPartitions()));
  // TODO need to allocate and free these ourselves. May need to update detail::MultiBlockContourTreeHelper::ComputeLocalContourTreeMesh
  std::vector<vtkm::worklet::contourtree_distributed::ContourTreeBlockData<T>*> localDataBlocks;
  localDataBlocks.resize(static_cast<size_t>(input.GetNumberOfPartitions()));
  std::vector<vtkmdiy::Link*> localLinks; // dummy links needed to make DIY happy
  localLinks.resize(static_cast<size_t>(input.GetNumberOfPartitions()));
  // We need to augment at least with the boundary vertices when running in parallel, even if the user requested at the end only the unaugmented contour tree
  unsigned int compRegularStruct =
    (this->ComputeRegularStructure > 0) ? this->ComputeRegularStructure : 2;

  for (std::size_t bi = 0; bi < static_cast<std::size_t>(input.GetNumberOfPartitions()); bi++)
  {
    // create the local contour tree mesh
    localLinks[bi] = new vtkmdiy::Link;
    auto currBlock = input.GetPartition(static_cast<vtkm::Id>(bi));
    auto currField =
      currBlock.GetField(this->GetActiveFieldName(), this->GetActiveFieldAssociation());

    vtkm::Id3 pointDimensions, globalPointDimensions, globalPointIndexStart;
    currBlock.GetCellSet().CastAndCallForTypes<VTKM_DEFAULT_CELL_SET_LIST_STRUCTURED>(
      vtkm::worklet::contourtree_augmented::GetLocalAndGlobalPointDimensions(),
      pointDimensions,
      globalPointDimensions,
      globalPointIndexStart);

    //const vtkm::cont::ArrayHandle<T,StorageType> &fieldData = currField.GetData().Cast<vtkm::cont::ArrayHandle<T,StorageType> >();
    vtkm::cont::ArrayHandle<T> fieldData;
    vtkm::cont::ArrayCopy(currField.GetData(), fieldData);
    auto currContourTreeMesh = vtkm::worklet::contourtree_distributed::MultiBlockContourTreeHelper::
      ComputeLocalContourTreeMesh<T>(globalPointIndexStart,
                                     pointDimensions,
                                     globalPointDimensions,
                                     fieldData,
                                     MultiBlockTreeHelper->LocalContourTrees[bi],
                                     MultiBlockTreeHelper->LocalSortOrders[bi],
                                     compRegularStruct);
    localContourTreeMeshes[bi] = currContourTreeMesh;
    // create the local data block structure
    localDataBlocks[bi] = new vtkm::worklet::contourtree_distributed::ContourTreeBlockData<T>();
    localDataBlocks[bi]->NumVertices = currContourTreeMesh->NumVertices;
    // localDataBlocks[bi]->SortOrder = currContourTreeMesh->SortOrder;
    localDataBlocks[bi]->SortedValue = currContourTreeMesh->SortedValues;
    localDataBlocks[bi]->GlobalMeshIndex = currContourTreeMesh->GlobalMeshIndex;
    localDataBlocks[bi]->NeighborConnectivity = currContourTreeMesh->NeighborConnectivity;
    localDataBlocks[bi]->NeighborOffsets = currContourTreeMesh->NeighborOffsets;
    localDataBlocks[bi]->MaxNeighbors = currContourTreeMesh->MaxNeighbors;
    localDataBlocks[bi]->BlockOrigin = globalPointIndexStart;
    localDataBlocks[bi]->BlockSize = pointDimensions;
    localDataBlocks[bi]->GlobalSize = globalPointDimensions;
    // We need to augment at least with the boundary vertices when running in parallel
    localDataBlocks[bi]->ComputeRegularStructure = compRegularStruct;
  }
  // Setup vtkmdiy to do global binary reduction of neighbouring blocks. See also RecuctionOperation struct for example

  // Create the vtkmdiy master
  vtkmdiy::Master master(comm,
                         1, // Use 1 thread, VTK-M will do the treading
                         -1 // All block in memory
  );

  // Compute the gids for our local blocks
  using RegularDecomposer = vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds>;

  RegularDecomposer::DivisionsVector diyDivisions;
  std::vector<int> vtkmdiyLocalBlockGids;
  vtkmdiy::DiscreteBounds diyBounds(0);
  if (this->MultiBlockTreeHelper->BlocksPerDimension[0] == -1)
  {
    VTKM_LOG_S(vtkm::cont::LogLevel::Info,
               "BlocksPerDimension not set. Computing block indices "
               "from information in CellSetStructured.");
    diyBounds = vtkm::filter::scalar_topology::internal::ComputeBlockIndices(
      input, diyDivisions, vtkmdiyLocalBlockGids);
  }
  else
  {
    VTKM_LOG_S(vtkm::cont::LogLevel::Info,
               "BlocksPerDimension set. Using information provided by caller.");
    diyBounds = vtkm::filter::scalar_topology::internal::ComputeBlockIndices(
      input,
      this->MultiBlockTreeHelper->BlocksPerDimension,
      this->MultiBlockTreeHelper->LocalBlockIndices,
      diyDivisions,
      vtkmdiyLocalBlockGids);
  }
  int numDims = diyBounds.min.dimension();
  int globalNumberOfBlocks =
    std::accumulate(diyDivisions.cbegin(), diyDivisions.cend(), 1, std::multiplies<int>{});

  // Add my local blocks to the vtkmdiy master.
  for (std::size_t bi = 0; bi < static_cast<std::size_t>(input.GetNumberOfPartitions()); bi++)
  {
    master.add(static_cast<int>(vtkmdiyLocalBlockGids[bi]), // block id
               localDataBlocks[bi],
               localLinks[bi]);
  }

  // Define the decomposition of the domain into regular blocks
  RegularDecomposer::BoolVector shareFace(3, true);
  RegularDecomposer::BoolVector wrap(3, false);
  RegularDecomposer::CoordinateVector ghosts(3, 1);
  RegularDecomposer decomposer(static_cast<int>(numDims),
                               diyBounds,
                               globalNumberOfBlocks,
                               shareFace,
                               wrap,
                               ghosts,
                               diyDivisions);

  // Define which blocks live on which rank so that vtkmdiy can manage them
  vtkmdiy::DynamicAssigner assigner(comm, static_cast<int>(size), globalNumberOfBlocks);
  for (vtkm::Id bi = 0; bi < input.GetNumberOfPartitions(); bi++)
  {
    assigner.set_rank(static_cast<int>(rank),
                      static_cast<int>(vtkmdiyLocalBlockGids[static_cast<size_t>(bi)]));
  }

  // Fix the vtkmdiy links. (NOTE: includes an MPI barrier)
  vtkmdiy::fix_links(master, assigner);

  // partners for merge over regular block grid
  vtkmdiy::RegularMergePartners partners(
    decomposer, // domain decomposition
    2,          // raix of k-ary reduction.
    true        // contiguous: true=distance doubling , false=distnace halving
  );
  // reduction
  vtkmdiy::reduce(
    master, assigner, partners, &vtkm::worklet::contourtree_distributed::MergeBlockFunctor<T>);

  comm.barrier(); // Be safe!

  if (rank == 0)
  {
    vtkm::Id3 dummy1, globalPointDimensions, dummy2;
    vtkm::cont::DataSet firstDS = input.GetPartition(0);
    firstDS.GetCellSet().CastAndCallForTypes<VTKM_DEFAULT_CELL_SET_LIST_STRUCTURED>(
      vtkm::worklet::contourtree_augmented::GetLocalAndGlobalPointDimensions(),
      dummy1,
      globalPointDimensions,
      dummy2);
    // Now run the contour tree algorithm on the last block to compute the final tree
    vtkm::Id currNumIterations;
    vtkm::worklet::contourtree_augmented::ContourTree currContourTree;
    vtkm::worklet::contourtree_augmented::IdArrayType currSortOrder;
    vtkm::worklet::ContourTreeAugmented worklet;
    vtkm::cont::ArrayHandle<T> currField;
    // Construct the contour tree mesh from the last block
    vtkm::worklet::contourtree_augmented::ContourTreeMesh<T> contourTreeMeshOut;
    contourTreeMeshOut.NumVertices = localDataBlocks[0]->NumVertices;
    contourTreeMeshOut.SortOrder = vtkm::cont::ArrayHandleIndex(contourTreeMeshOut.NumVertices);
    contourTreeMeshOut.SortIndices = vtkm::cont::ArrayHandleIndex(contourTreeMeshOut.NumVertices);
    contourTreeMeshOut.SortedValues = localDataBlocks[0]->SortedValue;
    contourTreeMeshOut.GlobalMeshIndex = localDataBlocks[0]->GlobalMeshIndex;
    contourTreeMeshOut.NeighborConnectivity = localDataBlocks[0]->NeighborConnectivity;
    contourTreeMeshOut.NeighborOffsets = localDataBlocks[0]->NeighborOffsets;
    contourTreeMeshOut.MaxNeighbors = localDataBlocks[0]->MaxNeighbors;
    // Construct the mesh boundary exectuion object needed for boundary augmentation
    vtkm::Id3 minIdx(0, 0, 0);
    vtkm::Id3 maxIdx = globalPointDimensions;
    maxIdx[0] = maxIdx[0] - 1;
    maxIdx[1] = maxIdx[1] - 1;
    maxIdx[2] = maxIdx[2] > 0 ? (maxIdx[2] - 1) : 0;
    auto meshBoundaryExecObj =
      contourTreeMeshOut.GetMeshBoundaryExecutionObject(globalPointDimensions, minIdx, maxIdx);
    // Run the worklet to compute the final contour tree
    std::cout << "using meshBoundaryExecObj inside worklet.Run()\n";
    std::cout << "... alternative suggested was copying this file and removing all meshBoundaryExecObj dependencies}\n";
    // worklet.Run function is defined in:
    worklet.Run(
      contourTreeMeshOut.SortedValues, // Unused param. Provide something to keep API happy
      contourTreeMeshOut,
      this->ContourTreeData,
      this->MeshSortOrder,
      currNumIterations,
      this->ComputeRegularStructure,
      meshBoundaryExecObj);

    // Set the final mesh sort order we need to use
    this->MeshSortOrder = contourTreeMeshOut.GlobalMeshIndex;
    // Remeber the number of iterations for the output
    this->NumIterations = currNumIterations;

    // Return the sorted values of the contour tree as the result
    // TODO the result we return for the parallel and serial case are different right now. This should be made consistent. However, only in the parallel case are we useing the result output
    vtkm::cont::DataSet temp;
    vtkm::cont::Field rfield(this->GetOutputFieldName(),
                             vtkm::cont::Field::Association::WholeDataSet,
                             contourTreeMeshOut.SortedValues);
    temp.AddField(rfield);
    output = vtkm::cont::PartitionedDataSet(temp);
  }
  else
  {
    this->ContourTreeData = MultiBlockTreeHelper->LocalContourTrees[0];
    this->MeshSortOrder = MultiBlockTreeHelper->LocalSortOrders[0];

    // Free allocated temporary pointers
    for (std::size_t bi = 0; bi < static_cast<std::size_t>(input.GetNumberOfPartitions()); bi++)
    {
      delete localContourTreeMeshes[bi];
      delete localDataBlocks[bi];
      // delete localLinks[bi];
    }
  }
  localContourTreeMeshes.clear();
  localDataBlocks.clear();
  localLinks.clear();
}

//-----------------------------------------------------------------------------
VTKM_CONT void ContourTreeAugmented::PostExecute(const vtkm::cont::PartitionedDataSet& input,
                                                 vtkm::cont::PartitionedDataSet& result)
{
  if (this->MultiBlockTreeHelper)
  {
    vtkm::cont::Timer timer;
    timer.Start();

    // We are running in parallel and need to merge the contour tree in PostExecute
    if (MultiBlockTreeHelper->GetGlobalNumberOfBlocks() == 1)
    {
      return;
    }

    auto field =
      input.GetPartition(0).GetField(this->GetActiveFieldName(), this->GetActiveFieldAssociation());

    // To infer and pass on the ValueType of the field.
    auto PostExecuteCaller = [&](const auto& concrete) {
      using T = typename std::decay_t<decltype(concrete)>::ValueType;
      this->DoPostExecute<T>(input, result);
    };
    this->CastAndCallScalarField(field, PostExecuteCaller);

    this->MultiBlockTreeHelper.reset();
    VTKM_LOG_S(vtkm::cont::LogLevel::Perf,
               std::endl
                 << "    " << std::setw(38) << std::left << "Contour Tree Filter PostExecute"
                 << ": " << timer.GetElapsedTime() << " seconds");
  }
}

} // namespace scalar_topology
} // namespace filter
} // namespace vtkm
