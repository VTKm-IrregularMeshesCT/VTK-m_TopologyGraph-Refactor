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

#ifndef vtk_m_worklet_ContourTreeUniformAugmented_h
#define vtk_m_worklet_ContourTreeUniformAugmented_h


#include <sstream>
#include <utility>

// VTKM includes
#include <vtkm/Math.h>
#include <vtkm/Types.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

// Contour tree worklet includes
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ActiveGraph.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ContourTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ContourTreeMaker.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/DataSetMesh.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/MergeTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/MeshExtrema.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/ContourTreeMesh.h>
// Adding the new TopologyGraph Class
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/TopologyGraph.h>

#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/MeshBoundary2D.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/MeshBoundary3D.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/MeshBoundaryContourTreeMesh.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/PrintVectors.h>

#include <chrono>
#include <thread>

#define PACT_DEBUG 0

using vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT;

namespace vtkm
{
namespace worklet
{

/// Compute the contour tree for 2d and 3d uniform grids and arbitrary topology graphs
class ContourTreeAugmented
{
public:
  /*!
  * Log level to be used for outputting timing information. Default is vtkm::cont::LogLevel::Perf
  * Use vtkm::cont::LogLevel::Off to disable outputing the results via vtkm logging here. The
  * results are saved in the TimingsLogString variable so we can use it to do our own logging
  */
  vtkm::cont::LogLevel TimingsLogLevel = vtkm::cont::LogLevel::Perf;

  /// Remember the results from our time-keeping so we can customize our logging
  std::string TimingsLogString;


  /*!
  * Run the contour tree to merge an existing set of contour trees
  *
  *  fieldArray   : Needed only as a pass-through value but not used in this case
  *  mesh : The ContourTreeMesh for which the contour tree should be computed
  *  contourTree  : The output contour tree to be computed (output)
  *  sortOrder    : The sort order for the mesh vertices (output)
  *  nIterations  : The number of iterations used to compute the contour tree (output)
  *  computeRegularStructure : 0=Off, 1=full augmentation with all vertices
  *                            2=boundary augmentation using meshBoundary
  *  meshBoundary : This parameter is generated by calling mesh.GetMeshBoundaryExecutionObject
  *                 For regular 2D/3D meshes this required no extra parameters, however, for a
  *                 ContourTreeMesh additional information about the block must be given. Rather
  *                 than generating the MeshBoundary descriptor here, we therefore, require it
  *                 as an input. The MeshBoundary is used to augment the contour tree with the
  *                 mesh boundary vertices. It is needed only if we want to augement by the
  *                 mesh boundary and computeRegularStructure is False (i.e., if we compute
  *                 the full regular strucuture this is not needed because all vertices
  *                 (including the boundary) will be addded to the tree anyways.
  */
  template <typename FieldType,
            typename StorageType,
            typename MeshType,
            typename MeshBoundaryMeshExecType>
  void Run(const vtkm::cont::ArrayHandle<FieldType, StorageType> fieldArray,
           MeshType& mesh,
           contourtree_augmented::ContourTree& contourTree,
           contourtree_augmented::IdArrayType& sortOrder,
           vtkm::Id& nIterations,
           unsigned int computeRegularStructure,
           const MeshBoundaryMeshExecType& meshBoundary)
  {
    RunContourTree(
      fieldArray, // Just a place-holder to fill the required field. Used when calling SortData on the contour tree which is a no-op
      contourTree,
      sortOrder,
      nIterations,
      mesh,
      computeRegularStructure,
      meshBoundary);
    return;
  }


//  template <typename FieldType, typename StorageType>
//  std::vector<vtkm::Id> sortInput(const vtkm::cont::ArrayHandle<FieldType, StorageType> fieldArray)
//  {
//      std::vector<vtkm::Id> sorted_vector;
////    optimise to resize before doing work:
////      sorted_vector.reserve(fieldArray);
////      sorted_vector.resize();





//      return sorted_vector;
//  }


//  VTKM_EXEC
//  inline vtkm::Pair<vtkm::Id, vtkm::Id> GetNeighbourComponentsMaskAndDegree(
//    vtkm::Id sortIndex, // note - this is used like the mesh index
//    bool getMaxComponents) const
//  { // GetNeighbourComponentsMaskAndDegree()
//    using namespace m2d_freudenthal;
//    // get data portals
//    // convert to a mesh index
//    vtkm::Id meshIndex = SortOrderPortal.Get(sortIndex);

//    // get the row and column
//    vtkm::Id2 pos = this->VertexPos(meshIndex);
//    vtkm::Int8 boundaryConfig = ((pos[0] == 0) ? LeftBit : 0) |
//      ((pos[0] == this->MeshSize[0] - 1) ? RightBit : 0) | ((pos[1] == 0) ? TopBit : 0) |
//      ((pos[1] == this->MeshSize[1] - 1) ? BottomBit : 0);

//    // and initialise the mask
//    vtkm::Id neighbourhoodMask = 0;
//    // in what follows, the boundary conditions always reset wasAscent
//    for (vtkm::Id edgeNo = 0; edgeNo < N_INCIDENT_EDGES; edgeNo++)
//    { // per edge
//      // ignore if at edge of data
//      if (!(boundaryConfig & EdgeBoundaryDetectionMasksPortal.Get(edgeNo)))
//      {
//        // calculate neighbour's ID and sort order
//        vtkm::Id nbrSortIndex = GetNeighbourIndex(sortIndex, edgeNo);

//        // if it's not a valid destination, ignore it
//        if (getMaxComponents ? (nbrSortIndex > sortIndex) : (nbrSortIndex < sortIndex))
//        {
//          // now set the flag in the neighbourhoodMask
//          neighbourhoodMask |= vtkm::Id{ 1 } << edgeNo;
//        }
//      }
//    } // per edge

//    // we now know which edges are outbound, so we count to get the outdegree
//    vtkm::Id outDegree = 0;
//    vtkm::Id neighbourComponentMask = 0;
//    // special case for local minimum
//    if (neighbourhoodMask == 0x3F)
//    {
//      outDegree = 1;
//    }
//    else
//    { // not a local minimum
//      if ((neighbourhoodMask & 0x30) == 0x20)
//      {
//        ++outDegree;
//        neighbourComponentMask |= vtkm::Id{ 1 } << 5;
//      }
//      if ((neighbourhoodMask & 0x18) == 0x10)
//      {
//        ++outDegree;
//        neighbourComponentMask |= vtkm::Id{ 1 } << 4;
//      }
//      if ((neighbourhoodMask & 0x0C) == 0x08)
//      {
//        ++outDegree;
//        neighbourComponentMask |= vtkm::Id{ 1 } << 3;
//      }
//      if ((neighbourhoodMask & 0x06) == 0x04)
//      {
//        ++outDegree;
//        neighbourComponentMask |= vtkm::Id{ 1 } << 2;
//      }
//      if ((neighbourhoodMask & 0x03) == 0x02)
//      {
//        ++outDegree;
//        neighbourComponentMask |= vtkm::Id{ 1 } << 1;
//      }
//      if ((neighbourhoodMask & 0x21) == 0x01)
//      {
//        ++outDegree;
//        neighbourComponentMask |= vtkm::Id{ 1 } << 0;
//      }
//    } // not a local minimum
//    vtkm::Pair<vtkm::Id, vtkm::Id> re(neighbourComponentMask, outDegree);
//    return re;
//  } // GetNeighbourComponentsMaskAndDegree()

int to1DIndex(int x, int y, int z, int width, int height) {
  return (z * width * height) + (y * width) + x;
}

// Function that returns a list of neighbor IDs for a given node ID.
std::vector<int> getNeighborIDs(int id, int width, int height, int depth) {
  std::vector<int> neighbors;
  int x, y, z;

  z = id / (width * height);
  id %= (width * height);
  y = id / width;
  x = id % width;

  // Check and add neighbor to the left
  if (x > 0)
      neighbors.push_back(to1DIndex(x - 1, y, z, width, height));
  // Check and add neighbor to the right
  if (x < width - 1)
      neighbors.push_back(to1DIndex(x + 1, y, z, width, height));
  // Check and add neighbor above
  if (y > 0)
      neighbors.push_back(to1DIndex(x, y - 1, z, width, height));
  // Check and add neighbor below
  if (y < height - 1)
      neighbors.push_back(to1DIndex(x, y + 1, z, width, height));
  // Check and add neighbor in front
  if (z > 0)
      neighbors.push_back(to1DIndex(x, y, z - 1, width, height));
  // Check and add neighbor behind
  if (z < depth - 1)
      neighbors.push_back(to1DIndex(x, y, z + 1, width, height));

  return neighbors;
}



// Function that returns an array with positions sorted from lowest to highest values
std::vector<int> getSortedPositions(const std::vector<double>& arr) {
    // Create a vector of pairs: value and original index
    std::vector<std::pair<int, int>> valueWithIndex;
    for (int i = 0; i < arr.size(); ++i) {
        valueWithIndex.emplace_back(arr[i], i);
    }

    // Sort the vector of pairs by the value, maintaining the original index
    std::sort(valueWithIndex.begin(), valueWithIndex.end(),
        [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a.first < b.first;
        });

    // Extract the sorted positions
    std::vector<int> sortedPositions;
    for (auto& vi : valueWithIndex) {
        sortedPositions.push_back(vi.second);
    }

    return sortedPositions;
}

// get "[sortID] = meshID" pairs
template <typename FieldType, typename StorageType>
std::vector<int> sortData(const vtkm::cont::ArrayHandle<FieldType, StorageType> fieldArray)
{
    // used for sorting, with epsilon*i
    std::vector<double> fieldArrayForSort;
    std::vector<int> resultSort;
    double epsilon = 0.00000001;

    // convert the fieldArray to doubles:
    int total_length = fieldArray.GetNumberOfValues();
    for(int i = 0; i < total_length; i++)
    {
        fieldArrayForSort.push_back((double)fieldArray.Get(i)+epsilon*i);
    }

    return getSortedPositions(fieldArrayForSort);
}

struct DelaunayMesh
{
    std::vector<vtkm::Id> std_nbor_connectivity;
    std::vector<vtkm::Id> std_nbor_offsets;
};


DelaunayMesh parseDelaunayASCII(const std::string& filePath)
{
    std::ifstream inputFile(filePath);
    std::string line;
    DelaunayMesh graph;
    size_t currentOffset = 0;

    if(!inputFile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return graph; // This will return empty vectors if file cannot be opened
    }

    while(std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        vtkm::Id id;

        // Push the currentOffset before processing the line
        graph.std_nbor_offsets.push_back(currentOffset);

        // Process each id in the current line
        while (iss >> id) {
            graph.std_nbor_connectivity.push_back(id);
            ++currentOffset; // Increment the offset for each id found
        }
    }

    // After processing all lines, the last offset should be equal to the total number of ids
    // This implies the end of the last vertex's neighborhood
    graph.std_nbor_offsets.push_back(currentOffset);

    return graph;

}



DelaunayMesh parseDelaunayDoubleASCII(const std::string& filePathUp,
                                      const std::string& filePathDown)
{
    std::ifstream inputFileUp(filePathUp);
    std::string lineUp;
    std::ifstream inputFileDown(filePathDown);
    std::string lineDown;

    DelaunayMesh graphUp;
    DelaunayMesh graphDown;
    DelaunayMesh graph;
    size_t currentOffsetUp = 0;
    size_t currentOffsetDown = 0;

    std::cout << "variables OK\n";

    if(!inputFileUp)
    {
        std::cerr << "Error opening file: " << filePathUp << std::endl;
        return graph; // This will return empty vectors if file cannot be opened
    }
    if(!inputFileDown)
    {
        std::cerr << "Error opening file: " << filePathDown << std::endl;
        return graph; // This will return empty vectors if file cannot be opened
    }

    std::cout << "files OK\n";

    // parse the up file:
    vtkm::Id processingID = 0;
    while(std::getline(inputFileUp, lineUp))
    {
        std::istringstream iss(lineUp);
        vtkm::Id id;
        vtkm::Id lineID;

        // Push the currentOffsetUp before processing the line
        graphUp.std_nbor_offsets.push_back(currentOffsetUp);

        // throw away the first element, since it is just the node ID:
        iss >> lineID;
        // check if equal to expected node (don't want to skip IDs):
        if(processingID != lineID)
        {
           std::cout << "skipped ID! " << processingID  << " " << lineID << "\n";
        }
        processingID++;

        // Process each id in the current line
        while (iss >> id)
        {
            graphUp.std_nbor_connectivity.push_back(id);
            ++currentOffsetUp; // Increment the offset for each id found
        }
    }

    // After processing all lines, the last offset should be equal to the total number of ids
    // This implies the end of the last vertex's neighborhood
    graphUp.std_nbor_offsets.push_back(currentOffsetUp);

    std::cout << "read up OK\n";

    // parse the down file:
    processingID = 0;
    while(std::getline(inputFileDown, lineDown))
    {
        std::istringstream iss(lineDown);
        vtkm::Id id;
        vtkm::Id lineID;

        // Push the currentOffsetDown before processing the line
        graphDown.std_nbor_offsets.push_back(currentOffsetDown);

        // throw away the first element, since it is just the node ID:
        iss >> lineID;
        // check if equal to expected node (don't want to skip IDs):
        if(processingID != lineID)
        {
           std::cout << "DOWN: skipped ID! " << processingID  << " " << lineID << "\n";
        }
        processingID++;
//        std::cout << lineID << " ";

        // Process each id in the current line
        while (iss >> id)
        {
            graphDown.std_nbor_connectivity.push_back(id);
//            std::cout << id << " ";
            ++currentOffsetDown; // Increment the offset for each id found
        }
//        std::cout << "\n";
    }

    // After processing all lines, the last offset should be equal to the total number of ids
    // This implies the end of the last vertex's neighborhood
    graphDown.std_nbor_offsets.push_back(currentOffsetDown);

    std::cout << "read down OK\n";

    std::cout << "processingID should match 125: " << processingID << "\n";

    // now merge both graphs:
    graph.std_nbor_offsets.push_back(graphUp.std_nbor_offsets[0] + graphDown.std_nbor_offsets[0]);
    for(vtkm::Id i = 0; i < processingID; i++)
    {
//        std::cout << i << " " << graphUp.std_nbor_offsets[i+1] << ", " << graphDown.std_nbor_offsets[i+1] << "\n";
        graph.std_nbor_offsets.push_back(graphUp.std_nbor_offsets[i+1] + graphDown.std_nbor_offsets[i+1]);
        // try without sorting at first ...
        for(vtkm::Id j = graphUp.std_nbor_offsets[i]; j < graphUp.std_nbor_offsets[i+1]; j++)
        {
            graph.std_nbor_connectivity.push_back(graphUp.std_nbor_connectivity[j]);
        }
        for(vtkm::Id j = graphDown.std_nbor_offsets[i]; j < graphDown.std_nbor_offsets[i+1]; j++)
        {
            graph.std_nbor_connectivity.push_back(graphDown.std_nbor_connectivity[j]);
        }
    }

    std::cout << "combined OK\n";


    return graph;

}






  /*!
   * Run the contour tree analysis. This helper function is used to
   * allow one to run the contour tree in a consistent fashion independent
   * of whether the data is 2D, 3D, or 3D_MC. This function initalizes
   * the approbritate mesh class from the contourtree_augmented worklet
   * and constructs ths mesh boundary exectuion object to be used. It the
   * subsequently calls RunContourTree method to compute the actual contour tree.
   *
   *  fieldArray   : Needed only as a pass-through value but not used in this case
   *  mesh : The ContourTreeMesh for which the contour tree should be computed
   *  contourTree  : The output contour tree to be computed (output)
   *  sortOrder    : The sort order for the mesh vertices (output)
   *  nIterations  : The number of iterations used to compute the contour tree (output)
   *  nRows        : Number of rows (i.e, x values) in the input mesh
   *  nCols        : Number of columns (i.e, y values) in the input mesh
   *  nSlices      : Number of slicex (i.e, z values) in the input mesh. Default is 1
   *                 to avoid having to set the nSlices for 2D input meshes
   *  useMarchingCubes : Boolean indicating whether marching cubes (true) or freudenthal (false)
   *                     connectivity should be used. Valid only for 3D input data. Default is false.
   *  computeRegularStructure : 0=Off, 1=full augmentation with all vertices
   *                            2=boundary augmentation using meshBoundary.
   */
  template <typename FieldType, typename StorageType>
  void Run(const vtkm::cont::ArrayHandle<FieldType, StorageType> fieldArray,
           contourtree_augmented::ContourTree& contourTree,
           contourtree_augmented::IdArrayType& sortOrder,
           vtkm::Id& nIterations,
           const vtkm::Id3 meshSize,
           bool useMarchingCubes = false,
           unsigned int computeRegularStructure = 1)
  {
    // use marching cubes for debugging the topology graph:
    useMarchingCubes = false;
    /// DEBUG PRINT std::cout << "{sc_tp/worklet/ContourTreeUniformAugmented.h : Running augmented CT ...\n";
    using namespace vtkm::worklet::contourtree_augmented;

    // 2D Contour Tree
    if (meshSize[2] == 1)
    {
      std::cout << "Creating 2D Mesh with dim: " << meshSize[0] << "x" << meshSize[1] << "\n";
      // Build the mesh and fill in the values
      // DataSetMeshTriangulation2DFreudenthal mesh(vtkm::Id2{ meshSize[0], meshSize[1] });
      //DataSetMeshTriangulation2DFreudenthal meshref(vtkm::Id2{ meshSize[0], meshSize[1] });

      std::cout << "Pairs:\n";
      for (int i = 0; i < 25; i++)
      {
          //vtkm::Pair<vtkm::Id, vtkm::Id> pair = meshref.GetNeighbourComponentsMaskAndDegree(i, false);

//          vtkm::Id comp   = pair.first();
//          vtkm::Id degree = pair.second();

//          std::cout << pair << "\n";
      }

      std::cout << "Contour Tree Mesh TEST from file ... \n";
      //vtkm::worklet::contourtree_augmented::PrintArrayHandle(
      PrintValues(
        "Field Array",
        fieldArray,
        -1,
        std::cout);

      //ContourTreeMesh mesh("ctinput.txt");
      // runtime error when reading CT from existing output file ...
      // ... (created by PPP-2.0 save function) ...
      //ContourTreeMesh<int> mesh("5x5.all.txt");
      // ... so trying the default constructor ...
      // ... which in the VTK-m version does nothing (no-op) ...
      // ... so ported values from PPP-2.0 to this one to try with a valid CT:
      //ContourTreeMesh<int> mesh;
      // 5x5 dataset values:
      // std::vector<int> nodes_sorted_values = {100,78,49,17,1,94,71,47,33,6,52,44,50,45,48,8,12,46,91,43,0,5,51,76,83};
      // std::vector<int> arcs =
      // Take this simple graph for example:
      //
      /* 4
          \
           \> 3 -> 1 <- 0
           /
          /
         2
        (Use this comment style to avoid warnings about multi-line comments triggered by '\' at
         the end of the line).
      */
      // make all nodes/sort/values the same for easier debugging:
      //std::vector<vtkm::Id> std_nodes_sorted = {0, 1, 2, 3, 4};
      //std::vector<vtkm::Id> std_nodes_sorted = {0, 1, 2, 3, 4, 5, 6, 7};
      // below is for 5x5:
//      std::vector<vtkm::Id> std_nodes_sorted = {24, 20, 14, 6, 1,
//                                                23, 18, 12, 7, 3,
//                                                17, 9,  15, 10, 13,
//                                                4,  5,  11, 22, 8,
//                                                0,  2,  16, 19, 21};

      std::vector<vtkm::Id> std_nodes_sorted = {0, 1, 2, 3, 4,
                                                5, 6, 7, 8, 9,
                                                10, 11, 12, 13, 14,
                                                15, 16, 17, 18, 19,
                                                20, 21, 22, 23, 24};

      // convert regular std::vector to an array that VTK likes ... picky
      vtkm::cont::ArrayHandle<vtkm::Id> nodes_sorted =
        vtkm::cont::make_ArrayHandle(std_nodes_sorted, vtkm::CopyFlag::Off);

      //IdArrayType nodes_sorted(0, 1, 2, 3, 4); doesnt work
      // take 2
      //std::vector<int> std_actual_values = {0, 1, 2, 3, 4};
      //std::vector<int> std_actual_values = {0, 1, 2, 3, 4, 5, 6, 7};
      // below for 5x5:
//      std::vector<int> std_actual_values = {24, 20, 14, 6, 1,
//                                            23, 18, 12, 7, 3,
//                                            17, 9,  15, 10, 13,
//                                            4,  5,  11, 22, 8,
//                                            0,  2,  16, 19, 21};

//      std::vector<int> std_actual_values = {0, 1, 2, 3, 4,
//                                            5, 6, 7, 8, 9,
//                                            10, 11, 12, 13, 14,
//                                            15, 16, 17, 18, 19,
//                                            20, 21, 22, 23, 24};

      std::vector<int> std_actual_values = {100, 78, 49, 17, 1,
                                            94, 71, 47, 33, 6,
                                            52, 44, 50, 45, 48,
                                            8, 12, 46, 91, 43,
                                            0, 5, 51, 76, 83};

      // convert regular std::vector to an array that VTK likes ... picky
      vtkm::cont::ArrayHandle<int> actual_values =
        vtkm::cont::make_ArrayHandle(std_actual_values, vtkm::CopyFlag::Off);
//          vtkm::cont::make_ArrayHandle(fieldArray, vtkm::CopyFlag::Off);

      // The contour tree algorithm stores this in an arcs array:
      //
      // idx:  0 1 2 3 4
      // arcs: 1 - 3 1 3 (- = NO_SUCH_ELEMENT, meaning no arc originating from this node)
      //std::vector<vtkm::Id> std_arcs_list = {1, NO_SUCH_ELEMENT,  3, 1, 3};
//      std::vector<vtkm::Id> std_arcs_list = {2, 2, 3, 4, 5, 7, 5, NO_SUCH_ELEMENT};
      // just part of my standard conversion routine at this point:
//      vtkm::cont::ArrayHandle<vtkm::Id> arcs_list =
//        vtkm::cont::make_ArrayHandle(std_arcs_list, vtkm::CopyFlag::Off);

      // after a lengthy process ... finally the last one:
      std::vector<vtkm::Id> std_global_inds = {0};
      vtkm::cont::ArrayHandle<vtkm::Id> global_inds =
        vtkm::cont::make_ArrayHandle(std_global_inds, vtkm::CopyFlag::Off);

      // uncomment if testing the ContourTreeMesh
//      ContourTreeMesh<int> mesh(nodes_sorted,
//                                arcs_list,
//                                nodes_sorted,
//                                actual_values,
//                                global_inds);

      // WARNING WARNING WARNING!
      // The following is to test the Topology Graph which is still in an unfinished state!
//      std::vector<vtkm::Id> std_nbor_connectivity = {2, 4,
//                                                     2, 3,
//                                                     0, 1, 3, 4,
//                                                     1, 2, 4, 5, 7,
//                                                     0, 2, 3, 5, 6,
//                                                     3, 4, 6, 7,
//                                                     4, 5,
//                                                     3, 5};
      // below for 5x5:
//      std::vector<vtkm::Id> std_nbor_connectivity = {18, 20, 23,                // [0] 24
//                                                     12, 14, 18, 24,            // [1] 20
//                                                     6,  7,  12, 20,            // [2] 14
//                                                     1,  3,  7,  14,            // [3] 6
//                                                     3,  6,                     // [4] 1
//                                                     9,  17, 18, 24,            // [5] 23
//                                                     9,  12, 15, 20, 23, 24,    // [6] 18
//                                                     7,  10, 14, 15, 18, 20,    // [7] 12
//                                                     3,   6, 10, 12, 13, 14,    // [8] 7
//                                                     1,  6,  7,  13,            // [9] 3
//                                                     4,  5,  9,  23,            // [10] 17
//                                                     5, 11, 15, 17, 18, 23,     // [11] 9
//                                                     9, 10, 11, 12, 18, 22,     // [12] 15
//                                                     7, 8, 12, 13, 15, 22,      // [13] 10
//                                                     3, 7, 8, 10,               // [14] 13
//                                                     0, 2, 5, 17,               // [15] 4
//                                                     2, 4, 9, 11, 16, 17,       // [16] 5
//                                                     5, 9, 15, 16, 19, 22,      // [17] 11
//                                                     8, 10, 11, 15, 19, 21,     // [18] 22
//                                                     10, 13, 21, 22,            // [19] 8
//                                                     2, 4,                      // [20] 0
//                                                     0, 4, 5, 16,               // [21] 2
//                                                     2, 5, 11, 19,              // [22] 16
//                                                     11, 16, 21, 22,            // [23] 19
//                                                     8, 19, 22};                // [24] 21

      std::vector<vtkm::Id> std_nbor_connectivity = {2,4,
                                                     3,6,
                                                     0,4,5,16,
                                                     1,6,7,13,
                                                     0,2,5,17,
                                                     2,4,9,11,16,17,
                                                     1,3,7,14,
                                                     3,6,10,12,13,14,
                                                     10,13,21,22,
                                                     5,11,15,17,18,23,
                                                     7,8,12,13,15,22,
                                                     5,9,15,16,19,22,
                                                     7,10,14,15,18,20,
                                                     3,7,8,10,
                                                     6,7,12,20,
                                                     9,10,11,12,18,22,
                                                     2,5,11,19,
                                                     4,5,9,23,
                                                     9,12,15,20,23,24,
                                                     11,16,21,22,
                                                     12,14,18,24,
                                                     8,19,22,
                                                     8,10,11,15,19,21,
                                                     9,17,18,24,
                                                     18,20,23,};                // [24] 21


//      std::vector<vtkm::Id> std_nbor_connectivity = {1,5,6,
//                                                     0,2,6,7,
//                                                     1,3,7,8,
//                                                     2,4,8,9,
//                                                     3,9,
//                                                     0,6,10,11,
//                                                     0,1,5,7,11,12,
//                                                     1,2,6,8,12,13,
//                                                     2,3,7,9,13,14,
//                                                     3,4,8,14,
//                                                     5,11,15,16,
//                                                     5,6,10,12,16,17,
//                                                     6,7,11,13,17,18,
//                                                     7,8,12,14,18,19,
//                                                     8,9,13,19,
//                                                     10,16,20,21,
//                                                     10,11,15,17,21,22,
//                                                     11,12,16,18,22,23,
//                                                     12,13,17,19,23,24,
//                                                     13,14,18,24,
//                                                     15,21,
//                                                     15,16,20,22,
//                                                     16,17,21,23,
//                                                     17,18,22,24,
//                                                     18,19,23,};                // [24] 21



      vtkm::cont::ArrayHandle<vtkm::Id> nbor_connectivity =
        vtkm::cont::make_ArrayHandle(std_nbor_connectivity, vtkm::CopyFlag::Off);

      //std::vector<vtkm::Id> std_nbor_offsets = {0, 2, 4, 8, 13, 18, 22, 24, 26};
      // below for 5x5:
    //std::vector<vtkm::Id> std_nbor_offsets = {0, 3, 7, 11, 15, 17, 21, 27, 33, 39, 43, 47, 53, 59, 65, 69, 73, 79, 85, 91, 95, 97, 101, 105, 109, 112};
      std::vector<vtkm::Id> std_nbor_offsets = {0,2,4,8,12,16,22,26,32,36,42,48,54,60,64,68,74,78,82,88,92,96,99,105,109,112};
      vtkm::cont::ArrayHandle<vtkm::Id> nbor_offsets =
        vtkm::cont::make_ArrayHandle(std_nbor_offsets, vtkm::CopyFlag::Off);

      // Using the 'TopologyGraph' constructor here ...
      // ... (to be moved to a separate class eventually)
      ContourTreeMesh<int> mesh(nodes_sorted,
                              //arcs_list,
                                nbor_connectivity,
                                nbor_offsets,
                                nodes_sorted,
                                // doesnt work out of the box:
                                // fieldArray, // testing fieldArray instead of manual actual_value
                                actual_values,
                                global_inds);

      //ContourTreeMesh<int> mesh("5x5.all.txt");
      std::cout << "Finished reading, printing content:\n";
      mesh.PrintContent(std::cout);
      std::cout << "Finished printing content\n";

      // Run the contour tree on the mesh
      RunContourTree(fieldArray,
                     contourTree,
                     sortOrder,
                     nIterations,
                     mesh,
                     computeRegularStructure,
                     //nullptr);
                     mesh.GetMeshBoundaryExecutionObject());
      std::cout << "sc_tp/worklet/ContourTreeUniformAugmented.h : Change to NULL ptr here}\n";
      return;
    }
    // 3D Contour Tree using marching cubes
    else if (useMarchingCubes)
    {
      std::cout << "USING MARCHING CUBES FOR 3D ...\n";

      std::cout << "Testing File Reading ...\n";
      // vtkmread
//      vtkm::io::VTKDataSetReader reader("data.vtk");

      //std::vector<int> sortID_to_meshID = sortData(fieldArray);
      //vtkm::cont::ArrayHandle<vtkm::Id> sortOrder_ =
      //  vtkm::cont::make_ArrayHandle(sortID_to_meshID, vtkm::CopyFlag::Off);

      std::vector<vtkm::Id> std_global_inds = {0};
      vtkm::cont::ArrayHandle<vtkm::Id> global_inds =
        vtkm::cont::make_ArrayHandle(std_global_inds, vtkm::CopyFlag::Off);

      // 5b-split
//      std::vector<vtkm::Id> std_nodes_sorted = {0, 1, 2,
//                                               3, 4, 5,
//                                               6, 7, 8,

//                                               9, 10, 11,
//                                               12, 13, 14,
//                                               15, 16, 17,

//                                               18, 19, 20,
//                                               21, 22, 23,
//                                               24, 25, 26};

//        5b-full
        std::vector<vtkm::Id> std_nodes_sorted = {  0,1,2,3,4,
                                                    5,6,7,8,9,
                                                    10,11,12,13,14,
                                                    15,16,17,18,19,
                                                    20,21,22,23,24,
                                                    25,26,27,28,29,
                                                    30,31,32,33,34,
                                                    35,36,37,38,39,
                                                    40,41,42,43,44,
                                                    45,46,47,48,49,
                                                    50,51,52,53,54,
                                                    55,56,57,58,59,
                                                    60,61,62,63,64,
                                                    65,66,67,68,69,
                                                    70,71,72,73,74,
                                                    75,76,77,78,79,
                                                    80,81,82,83,84,
                                                    85,86,87,88,89,
                                                    90,91,92,93,94,
                                                    95,96,97,98,99,
                                                    100,101,102,103,104,
                                                    105,106,107,108,109,
                                                    110,111,112,113,114,
                                                    115,116,117,118,119,
                                                    120,121,122,123,124};

      vtkm::cont::ArrayHandle<vtkm::Id> nodes_sorted =
        vtkm::cont::make_ArrayHandle(std_nodes_sorted, vtkm::CopyFlag::Off);

//      // 5b split
//      std::vector<int> std_actual_values = {0, 1, 2,
//                                                 3, 4, 5,
//                                                 6, 7, 8,

//                                                 9, 10, 11,
//                                                 12, 13, 14,
//                                                 15, 16, 17,

//                                                 18, 19, 20,
//                                                 21, 22, 23,
//                                                 24, 25, 26};

//      5b-full
      std::vector<int> std_actual_values = {0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,252,229,216,0,
                                            0,242,204,242,0,
                                            0,216,229,252,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,191,127,165,0,
                                            0,140,38,114,0,
                                            0,153,102,178,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,247,221,209,0,
                                            0,234,196,234,0,
                                            0,209,221,247,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0,
                                            0,0,0,0,0};

      // convert regular std::vector to an array that VTK likes ... picky
      vtkm::cont::ArrayHandle<int> actual_values =
        vtkm::cont::make_ArrayHandle(std_actual_values, vtkm::CopyFlag::Off);

      //PrintValues(sortOrder_);
      // 5b-split
//      std::vector<vtkm::Id> std_nbor_offsets = {0,3,9,14,20,30,38,43,51,57,64,74,81,91,98,104,112,117,125,130,133,139,145,155,162,172,182,196};
      // 5b-split
              //{0,3,7,10,14,19,23,26,30,33,37,42,46,51,55,58,62,65,69,72,75,79,83,88,92,97,102, 108}; //{0,2,4,8,12,16,22,26,32,36,42,48,54,60,64,68,74,78,82,88,92,96,99,105,109,112};

//      5b-full on 3D delaunay:
      std::vector<vtkm::Id> std_nbor_offsets = {0,5,13,21,29,36,43,53,64,74,83,91,102,114,126,135,144,154,166,177,187,194,202,210,217,222,229,239,249,259,268,279,291,301,313,322,334,341,353,366,377,384,391,401,411,421,430,442,454,464,476,488,500,507,519,531,543,551,558,569,579,589,598,609,620,631,642,655,665,674,686,698,709,718,723,732,738,744,748,754,766,777,788,794,800,812,824,835,842,848,860,873,885,894,899,906,914,924,931,945,958,973,987,1001,1016,1030,1046,1061,1075,1089,1104,1119,1133,1148,1162,1178,1192,1206,1219,1235,1249,1263,1277,1293,1307,1322};
//      5b full on marching cubes:
      //{0,3,7,11,15,18,22,27,32,37,41,45,50,55,60,64,68,73,78,83,87,90,94,98,102,105,109,114,119,124,128,133,138,143,148,153,158,162,167,172,177,181,185,190,195,200,204,209,214,219,224,229,234,238,243,248,253,257,261,266,271,276,280,285,290,295,300,305,310,314,319,324,329,333,336,340,344,348,351,355,360,365,370,374,378,383,388,393,397,401,406,411,416,420,423,427,431,435,438,444,450,456,462,468,474,480,486,492,498,504,510,516,522,528,534,540,546,552,558,564,570,576,582,588,594,600};
      vtkm::cont::ArrayHandle<vtkm::Id> nbor_offsets =
        vtkm::cont::make_ArrayHandle(std_nbor_offsets, vtkm::CopyFlag::Off);

      // 5b-split
//      std::vector<vtkm::Id> std_nbor_connectivity = {1,3,9,
//          0,2,3,4,9,10,
//          1,4,5,10,11,
//          0,1,4,6,9,12,
//          1,2,3,5,6,7,9,10,12,26,
//          2,4,7,8,10,11,24,26,
//          3,4,7,12,13,
//          4,5,6,8,12,13,25,26,
//          5,7,23,24,25,26,
//          0,1,3,4,10,12,14,
//          1,2,4,5,9,11,12,14,15,26,
//          2,5,10,15,16,24,26,
//          3,4,6,7,9,10,13,14,17,26,
//          6,7,12,17,18,25,26,
//          9,10,12,15,17,26,
//          10,11,14,16,17,22,24,26,
//          11,15,20,22,24,
//          12,13,14,15,18,22,25,26,
//          13,17,21,22,25,
//          20,21,23,
//          16,19,21,22,23,24,
//          18,19,20,22,23,25,
//          15,16,17,18,20,21,23,24,25,26,
//          8,19,20,21,22,24,25,
//          5,8,11,15,16,20,22,23,25,26,
//          7,8,13,17,18,21,22,23,24,26,
//          4,5,7,8,10,11,12,13,14,15,17,22,24,25,
//      };
//      {1,3,9,
//                                                            0,2,4,10,
//                                                            1,5,11,
//                                                            0,4,6,12,
//                                                            1,3,5,7,26,
//                                                            2,4,8,24,
//                                                            3,7,13,
//                                                            4,6,8,25,
//                                                            5,7,23,
//                                                            0,10,12,14,
//                                                            1,9,11,15,26,
//                                                            2,10,16,24,
//                                                            3,9,13,17,26,
//                                                            6,12,18,25,
//                                                            9,15,17,
//                                                            10,14,16,22,
//                                                            11,15,20,
//                                                            12,14,18,22,
//                                                            13,17,21,
//                                                            20,21,23,
//                                                            16,19,22,24,
//                                                            18,19,22,25,
//                                                            15,17,20,21,26,
//                                                            8,19,24,25,
//                                                            5,11,20,23,26,
//                                                            7,13,21,23,26,
//                                                            4,10,12,22,24,25,};         // [24] 21

//      5b full
      std::vector<vtkm::Id> std_nbor_connectivity = {1,5,6,25,30,
          0,2,5,6,25,26,30,124,
          1,3,6,7,26,27,115,124,
          2,4,7,8,27,28,111,115,
          3,8,9,28,29,31,111,
          0,1,6,10,11,25,30,
          0,1,2,5,7,10,11,12,30,124,
          2,3,6,8,11,12,13,108,115,119,124,
          3,4,7,9,12,13,108,111,115,120,
          4,8,13,14,29,31,33,111,120,
          5,6,11,15,16,30,32,124,
          5,6,7,10,12,15,16,17,32,119,124,
          6,7,8,11,13,16,17,18,108,112,116,119,
          7,8,9,12,14,17,18,19,108,116,120,123,
          9,13,18,19,31,33,35,120,123,
          10,11,16,20,21,32,34,112,119,
          10,11,12,15,17,20,21,22,112,119,
          11,12,13,16,18,21,22,23,37,38,112,116,
          12,13,14,17,19,22,23,38,39,116,123,
          13,14,18,23,24,33,35,39,40,123,
          15,16,21,34,36,37,112,
          15,16,17,20,22,36,37,112,
          16,17,18,21,23,37,38,116,
          17,18,19,22,24,38,39,
          19,23,35,39,40,
          0,1,5,26,30,41,46,
          1,2,25,27,30,41,42,46,106,124,
          2,3,26,28,42,43,101,106,115,124,
          3,4,27,29,43,44,101,104,111,115,
          4,9,28,31,44,45,47,104,111,
          0,1,5,6,10,25,26,32,41,46,124,
          4,9,14,29,33,45,47,49,100,104,111,120,
          10,11,15,30,34,46,48,106,119,124,
          9,14,19,31,35,47,49,51,100,105,120,123,
          15,20,32,36,48,50,102,112,119,
          14,19,24,33,39,40,49,51,55,56,105,123,
          20,21,34,37,50,52,112,
          17,20,21,22,36,38,50,52,53,103,112,116,
          17,18,22,23,37,39,53,54,55,99,103,116,123,
          18,19,23,24,35,38,40,54,55,56,123,
          19,24,35,39,51,55,56,
          25,26,30,42,46,57,58,
          26,27,41,43,46,57,58,106,121,124,
          27,28,42,44,58,59,101,106,113,121,
          28,29,43,45,59,60,101,104,109,113,
          29,31,44,47,60,61,63,104,109,
          25,26,30,32,41,42,48,57,58,62,106,124,
          29,31,33,45,49,61,63,65,100,104,109,118,
          32,34,46,50,62,64,102,106,119,121,
          31,33,35,47,51,63,65,67,100,105,118,122,
          34,36,37,48,52,53,64,66,102,103,112,117,
          33,35,40,49,55,56,65,67,72,100,105,122,
          36,37,50,53,66,68,103,
          37,38,50,52,54,66,68,69,70,99,103,110,
          38,39,53,55,69,70,71,99,105,110,114,123,
          35,38,39,40,51,54,56,70,71,72,105,123,
          35,39,40,51,55,71,72,105,
          41,42,46,58,62,73,74,
          41,42,43,46,57,59,62,73,74,106,121,
          43,44,58,60,74,75,79,80,113,121,
          44,45,59,61,75,76,80,81,109,113,
          45,47,60,63,76,77,81,82,109,
          46,48,57,58,64,73,74,78,79,106,121,
          45,47,49,61,65,81,82,86,87,109,118,
          48,50,62,66,78,79,83,84,102,117,121,
          47,49,51,63,67,86,87,92,109,118,122,
          50,52,53,64,68,69,83,84,88,89,103,110,117,
          49,51,65,72,91,92,97,105,118,122,
          52,53,66,69,88,89,93,94,110,
          53,54,66,68,70,89,90,93,94,95,110,114,
          53,54,55,69,71,90,94,95,96,110,114,122,
          54,55,56,70,72,95,96,97,105,114,122,
          51,55,56,67,71,96,97,105,122,
          57,58,62,74,78,
          57,58,59,62,73,75,78,79,121,
          59,60,74,76,79,80,
          60,61,75,77,80,81,
          61,76,81,82,
          62,64,73,74,79,83,
          59,62,64,74,75,78,80,83,84,113,117,121,
          59,60,75,76,79,81,84,85,107,109,113,
          60,61,63,76,77,80,82,85,86,87,109,
          61,63,77,81,86,87,
          64,66,78,79,84,88,
          64,66,79,80,83,85,88,89,107,110,113,117,
          80,81,84,86,89,90,91,107,109,110,114,118,
          63,65,81,82,85,87,90,91,92,109,118,
          63,65,81,82,86,91,92,
          66,68,83,84,89,93,
          66,68,69,84,85,88,90,93,94,95,110,114,
          69,70,85,86,89,91,94,95,96,110,114,118,122,
          67,85,86,87,90,92,95,96,97,114,118,122,
          65,67,86,87,91,96,97,118,122,
          68,69,88,89,94,
          68,69,70,89,90,93,95,
          69,70,71,89,90,91,94,96,
          70,71,72,90,91,92,95,97,114,122,
          67,71,72,91,92,96,122,
          99,100,101,102,103,104,107,108,111,113,115,116,117,120,
          38,53,54,98,100,103,105,107,110,114,116,120,123,
          31,33,47,49,51,98,99,104,105,107,111,114,118,120,122,
          27,28,43,44,98,102,104,106,108,111,113,115,117,121,
          34,48,50,64,98,101,103,106,108,112,115,117,119,121,
          37,38,50,52,53,66,98,99,102,107,108,110,112,116,117,
          28,29,31,44,45,47,98,100,101,107,109,111,113,118,
          33,35,49,51,54,55,56,67,71,72,99,100,114,120,122,123,
          26,27,32,42,43,46,48,58,62,101,102,115,119,121,124,
          80,84,85,98,99,100,103,104,109,110,113,114,117,118,
          7,8,12,13,98,101,102,103,111,112,115,116,119,120,
          44,45,47,60,61,63,65,80,81,85,86,104,107,113,118,
          53,54,66,68,69,70,84,85,89,90,99,103,107,114,117,
          3,4,8,9,28,29,31,98,100,101,104,108,115,120,
          12,15,16,17,20,21,34,36,37,50,102,103,108,116,119,
          43,44,59,60,79,80,84,98,101,104,107,109,117,121,
          54,69,70,71,85,89,90,91,96,99,100,105,107,110,118,122,
          2,3,7,8,27,28,98,101,102,106,108,111,119,124,
          12,13,17,18,22,37,38,98,99,103,108,112,120,123,
          50,64,66,79,84,98,101,102,103,107,110,113,121,
          47,49,63,65,67,85,86,90,91,92,100,104,107,109,114,122,
          7,11,12,15,16,32,34,48,102,106,108,112,115,124,
          8,9,13,14,31,33,98,99,100,105,108,111,116,123,
          42,43,48,58,59,62,64,74,79,101,102,106,113,117,
          49,51,65,67,70,71,72,90,91,92,96,97,100,105,114,118,
          13,14,18,19,33,35,38,39,54,55,99,105,116,120,
          1,2,6,7,10,11,26,27,30,32,42,46,106,115,119
      };

// //5b full on marching cubes manually
//        {1,5,25,
//                                                       0,2,6,26,
//                                                       1,3,7,27,
//                                                       2,4,8,28,
//                                                       3,9,29,
//                                                       0,6,10,30,
//                                                       1,5,7,11,124,
//                                                       2,6,8,12,115,
//                                                       3,7,9,13,111,
//                                                       4,8,14,31,
//                                                       5,11,15,32,
//                                                       6,10,12,16,119,
//                                                       7,11,13,17,108,
//                                                       8,12,14,18,120,
//                                                       9,13,19,33,
//                                                       10,16,20,34,
//                                                       11,15,17,21,112,
//                                                       12,16,18,22,116,
//                                                       13,17,19,23,123,
//                                                       14,18,24,35,
//                                                       15,21,36,
//                                                       16,20,22,37,
//                                                       17,21,23,38,
//                                                       18,22,24,39,
//                                                       19,23,40,
//                                                       0,26,30,41,
//                                                       1,25,27,42,124,
//                                                       2,26,28,43,115,
//                                                       3,27,29,44,111,
//                                                       4,28,31,45,
//                                                       5,25,32,46,124,
//                                                       9,29,33,47,111,
//                                                       10,30,34,48,119,
//                                                       14,31,35,49,120,
//                                                       15,32,36,50,112,
//                                                       19,33,40,51,123,
//                                                       20,34,37,52,
//                                                       21,36,38,53,112,
//                                                       22,37,39,54,116,
//                                                       23,38,40,55,123,
//                                                       24,35,39,56,
//                                                       25,42,46,57,
//                                                       26,41,43,58,106,
//                                                       27,42,44,59,101,
//                                                       28,43,45,60,104,
//                                                       29,44,47,61,
//                                                       30,41,48,62,106,
//                                                       31,45,49,63,104,
//                                                       32,46,50,64,102,
//                                                       33,47,51,65,100,
//                                                       34,48,52,66,103,
//                                                       35,49,56,67,105,
//                                                       36,50,53,68,
//                                                       37,52,54,69,103,
//                                                       38,53,55,70,99,
//                                                       39,54,56,71,105,
//                                                       40,51,55,72,
//                                                       41,58,62,73,
//                                                       42,57,59,74,121,
//                                                       43,58,60,75,113,
//                                                       44,59,61,76,109,
//                                                       45,60,63,77,
//                                                       46,57,64,78,121,
//                                                       47,61,65,82,109,
//                                                       48,62,66,83,117,
//                                                       49,63,67,87,118,
//                                                       50,64,68,88,110,
//                                                       51,65,72,92,122,
//                                                       52,66,69,93,
//                                                       53,68,70,94,110,
//                                                       54,69,71,95,114,
//                                                       55,70,72,96,122,
//                                                       56,67,71,97,
//                                                       57,74,78,
//                                                       58,73,75,79,
//                                                       59,74,76,80,
//                                                       60,75,77,81,
//                                                       61,76,82,
//                                                       62,73,79,83,
//                                                       74,78,80,84,121,
//                                                       75,79,81,85,113,
//                                                       76,80,82,86,109,
//                                                       63,77,81,87,
//                                                       64,78,84,88,
//                                                       79,83,85,89,117,
//                                                       80,84,86,90,107,
//                                                       81,85,87,91,118,
//                                                       65,82,86,92,
//                                                       66,83,89,93,
//                                                       84,88,90,94,110,
//                                                       85,89,91,95,114,
//                                                       86,90,92,96,122,
//                                                       67,87,91,97,
//                                                       68,88,94,
//                                                       69,89,93,95,
//                                                       70,90,94,96,
//                                                       71,91,95,97,
//                                                       72,92,96,
//                                                       99,100,101,102,107,108,
//                                                       54,98,103,105,114,116,
//                                                       49,98,104,105,118,120,
//                                                       43,98,104,106,113,115,
//                                                       48,98,103,106,117,119,
//                                                       50,53,99,102,110,112,
//                                                       44,47,100,101,109,111,
//                                                       51,55,99,100,122,123,
//                                                       42,46,101,102,121,124,
//                                                       85,98,113,114,117,118,
//                                                       12,98,115,116,119,120,
//                                                       60,63,81,104,113,118,
//                                                       66,69,89,103,114,117,
//                                                       8,28,31,104,115,120,
//                                                       16,34,37,103,116,119,
//                                                       59,80,101,107,109,121,
//                                                       70,90,99,107,110,122,
//                                                       7,27,101,108,111,124,
//                                                       17,38,99,108,112,123,
//                                                       64,84,102,107,110,121,
//                                                       65,86,100,107,109,122,
//                                                       11,32,102,108,112,124,
//                                                       13,33,100,108,111,123,
//                                                       58,62,79,106,113,117,
//                                                       67,71,91,105,114,118,
//                                                       18,35,39,105,116,120,
//                                                       6,26,30,106,115,119};

      vtkm::cont::ArrayHandle<vtkm::Id> nbor_connectivity =
        vtkm::cont::make_ArrayHandle(std_nbor_connectivity, vtkm::CopyFlag::Off);

      // Build the mesh and fill in the values
      //DataSetMeshTriangulation3DMarchingCubes mesh(meshSize);

      // Using the 'TopologyGraph' constructor here ...
      // ... (to be moved to a separate class eventually)
      ContourTreeMesh<int> mesh(nodes_sorted,
                              //arcs_list,
                                nbor_connectivity,
                                nbor_offsets,
                                nodes_sorted,
                                // doesnt work out of the box:
                                // fieldArray, // testing fieldArray instead of manual actual_value
                                actual_values,
                                //nodes_sorted,
                                global_inds);



      // Run the contour tree on the mesh
      RunContourTree(fieldArray,
                     contourTree,
                     sortOrder,
                     nIterations,
                     mesh,
                     computeRegularStructure,
                     mesh.GetMeshBoundaryExecutionObject());

//      // Run the contour tree on the mesh
//      RunContourTree(fieldArray,
//                     contourTree,
//                     sortOrder,
//                     nIterations,
//                     mesh,
//                     computeRegularStructure,
//                     //nullptr);
//                     mesh.GetMeshBoundaryExecutionObject());
      return;
    }
    // 3D Contour Tree with Freudenthal
    else
    {      
//        // Uncomment below for Freudenthal tests:
//        // NOTE: THIS EXPLICITLY CALLS THE CODE IN INITIALIZEACTIVEEDGES.H
//        // ...   BE CAREFUL TO REMOVE THE MESH SAVING CODE IF CRASHING
//        /// DEBUG PRINT std::cout << "USING 3D Freudenthal...\n";
//        // Build the mesh and fill in the values
//        DataSetMeshTriangulation3DFreudenthal mesh(meshSize);
//        /// DEBUG PRINT std::cout << "Freudenthal Mesh:\n";
//        // Run the contour tree on the mesh
//        RunContourTree(fieldArray,
//                         contourTree,
//                         sortOrder,
//                         nIterations,
//                         mesh,
//                         computeRegularStructure,
//                         mesh.GetMeshBoundaryExecutionObject());

//        return;


      // delaunay-mesh
      // read from: output-connections2.txt
//      int num_datapoints = 270985;
//      int num_datapoints = 125;
//      int num_datapoints = 25;
//        int num_datapoints = 27;
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/output-connections3.txt";
//      int num_datapoints = 65;
//      int num_datapoints = 4984;
//      int num_datapoints = 2456;
//      int num_datapoints = 1251;
//      int num_datapoints = 626;
//      int num_datapoints = 917;
//      int num_datapoints = 1026;
//      int num_datapoints = 128613;
//      int num_datapoints = 64995;
//      int num_datapoints = 94182;
//      int num_datapoints = 113646;
//      int num_datapoints = 104418;
//      int num_datapoints = 109688;
//      int num_datapoints = 111981;
//      int num_datapoints = 111155;
//      int num_datapoints = 111472;
//        int num_datapoints = 2160930;
//        int num_datapoints = 2384337;

//      int num_datapoints = 1781876;

//        int num_datapoints = 29791;
//        int num_datapoints = 56770560;
//        int num_datapoints = 13824;
//        int num_datapoints = 110592;
//        int num_datapoints = 884736;
//        int num_datapoints = 100001;
//        int num_datapoints = 100001;
//          int num_datapoints = 70;
//          int num_datapoints = 90318;

//       int num_datapoints = 16;

//       int num_datapoints = 9; // 3x3 2D Branch tet Volume dataset
//       int num_datapoints = 8; // 2x2x2 3D Branch tet Volume dataset

//       int num_datapoints = 270001; // 270k-jittered-cheat-sorted-sequential-values

//       const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build-4/BPECT-NastyW-16-connections.txt";
//       const std::string filename = "/home/user/HCTC/VTK-m-topology/vtkm-build-4/BPECT-NastyW-16-connections.txt";

//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-triang.txt";

    // 3x3 2D Branch tet Volume dataset
//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-3x3.txt";
    // 2x2x2 3D Branch tet Volume dataset
//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Cube-8-2x2x2.txt";

    //    270k-jittered-cheat-sorted-sequential-values
//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/270k-jittered-cheat-sorted-sequential-values-CONNECTIVITY.txt";
    // The above didn't work (Process 'Killed', probably because of excessive memory usage
    // Trying to scale up from 8 data points upwards, see where the memory limit it and measure peak memory usage etc.
    // ... to see if possible to run it on a supercomputer

    // PACTBD-EDIT
    int num_datapoints = 1001;
//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/8-from-2M-sampled-excel-sorted.1-CONNECTIVITY.txt";
    // PACTBD-EDIT
//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/101-from-2M-sampled-excel-sorted.1-CONNECTIVITY.txt";
    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/1k-from-2M-sampled-excel-sorted.1-CONNECTIVITY.txt";


//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/output-connections3.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/24x24x24-gridded-13k-output-connections.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/48x48x48-gridded-110k-output-connections.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/96x96x96-gridded-884k-output-connections.txt";

//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-freu3d.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-freud3d-duped.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-freud3d-duped2-bugfix.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/31x31x31-freud3d.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-delaunay3D.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/31x31x31-gridded-29k-output-connections.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/31x31x31-gridded-29k-output-connections-duped.txt";
//        const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/31x31x31-gridded-29k-output-connections-dupegd.txt";

      // LONG TERM RUNNING EXAMPLE 5x5 DATASET MANUAL
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5x5.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/64x64x64-parcels-130k-samplesoutput-connections2.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/64x64x64-parcels-270k-full-output-connections.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/128x128x128-parcels-2M-full-output-connections.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/128x128x128-parcels-100k-samples-output-connections.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/64x64x64-parcels-112k-62-samplesoutput-connections2.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/128x128x128-tetgen-2M-output-connections.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/2M-parcels-output-connections.txt";

//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/128x128x128-parcels-clean-2M-output-connections.txt";
//      const std::string filename = "/work/e710/e710/ddilys/PACT/VTK-m-topology/vtkm-build/128x128x128-parcels-clean-2M-output-connections.txt";



//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/VTK-m_TopologyGraph/examples/contour_tree_augmented/build/PACT_nbors_final_sort.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/VTK-m_TopologyGraph/examples/contour_tree_augmented/build/PACT_nbors_90k_final_sort.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/VTK-m_TopologyGraph/examples/contour_tree_augmented/build/PACT_nbors_90k_final_sort_duplicateCT.txt";
//      const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/VTK-m_TopologyGraph/examples/contour_tree_augmented/build/90k-dup.txt";

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


      DelaunayMesh delmesh = parseDelaunayASCII(filename);


//      // NOTE: putting DOWN as filename1 fixed a problem when running it on hh24-noE!
//      // combine up and down freudenthals into a single graph:
////        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-31x-down.txt";
////        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput384-down.txt";
////        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-384-up.txt";
////        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-hh24-noE-down.txt";
//        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-hh96-E-down.txt";

////            const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-31x-up.txt";
////        const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput384-up.txt";
////        const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-384-down.txt";
////        const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-hh24-noE-up.txt";
//        const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-hh96-E-up.txt";
////      const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/MeshOutput-up.txt";
//      std::cout << "Double ASCII reading: ... \n";
//      DelaunayMesh delmesh = parseDelaunayDoubleASCII(filename1, filename2);



//      // debug reader:
//      std::cout << "std_nbor_connectivity\n";
//      for(unsigned i = 0; i < delmesh.std_nbor_connectivity.size(); i++)
//      {
//          std::cout << i << " " << delmesh.std_nbor_connectivity[i] << "\n";
//      }
//      std::cout << "std_nbor_offsets\n";
//      for(unsigned i = 0; i < delmesh.std_nbor_offsets.size(); i++)
//      {
//          std::cout << i << " " << delmesh.std_nbor_offsets[i] << "\n";
//      }


//      std::vector<vtkm::Id> std_nbor_connectivity;
      vtkm::cont::ArrayHandle<vtkm::Id> nbor_connectivity =
        vtkm::cont::make_ArrayHandle(delmesh.std_nbor_connectivity, vtkm::CopyFlag::Off);

//      delmesh.std_nbor_connectivity.DebugPrint("std::cout");
//      vtkm::worklet::contourtree_augmented::PrintIndices(
//        "std_nbor_connectivity", nbor_connectivity);
//      nbor_connectivity.DebugPrint("std::cout");

//      std::cout << "Debug: delmesh.std_nbor_connectivity\n";
//      vtkm::worklet::contourtree_augmented::PrintValues(
//        "std_nbor_connectivity", nbor_connectivity);

      std::cout << "nbor_connectivity num vals: " << nbor_connectivity.GetNumberOfValues() << "\n";

      // ------------------------------------------------------------------------------- //

//      std::vector<vtkm::Id> std_nbor_offsets;
      vtkm::cont::ArrayHandle<vtkm::Id> nbor_offsets =
        vtkm::cont::make_ArrayHandle(delmesh.std_nbor_offsets, vtkm::CopyFlag::Off);

//      std::cout << "Debug: delmesh.std_nbor_offsets\n";
//      vtkm::worklet::contourtree_augmented::PrintValues(
//        "delmesh.std_nbor_offsets", nbor_offsets);

      std::cout << "nbor_offsets num vals: " << nbor_offsets.GetNumberOfValues() << "\n";

      // 112 nbor sleep
      std::this_thread::sleep_for(std::chrono::seconds(1));

//       Using the 'TopologyGraph' constructor here ...
//       ... (to be moved to a separate class eventually)


//// ----------------------------- MESHING VALIDATION ----------------------------- //
//      Uncomment below for PACT tests:
//      PACT:
      std::cout << "USING PACT (unoptimized)...\n";
      ContourTreeMesh<int> mesh(nodes_sorted,
                              //arcs_list,
                                nbor_connectivity,
                                nbor_offsets,
                                nodes_sorted,
                                // doesnt work out of the box:
                                // fieldArray, // testing fieldArray instead of manual actual_value
                                actual_values,
                                //nodes_sorted,
                                global_inds);

      // Run the contour tree on the mesh
      // PACT:
      RunContourTree(fieldArray,
                     contourTree,
                     sortOrder,
                     nIterations,
                     mesh,
                     computeRegularStructure,
                     mesh.GetMeshBoundaryExecutionObject());

//// Uncomment below for Freudenthal tests:
//// NOTE: THIS EXPLICITLY CALLS THE CODE IN INITIALIZEACTIVEEDGES.H
//// ...   BE CAREFUL TO REMOVE THE MESH SAVING CODE IF CRASHING
//      std::cout << "USING 3D Freudenthal...\n";
////       Build the mesh and fill in the values
//      DataSetMeshTriangulation3DFreudenthal mesh(meshSize);
//      std::cout << "Freudenthal Mesh:\n";
////       Run the contour tree on the mesh
//      RunContourTree(fieldArray,
//                     contourTree,
//                     sortOrder,
//                     nIterations,
//                     mesh,
//                     computeRegularStructure,
//                     mesh.GetMeshBoundaryExecutionObject());


//// Uncomment below for Marching Cubes tests:
//      std::cout << "USING 3D Marching Cubes...\n";
//      // Build the mesh and fill in the values
//      DataSetMeshTriangulation3DMarchingCubes mesh(meshSize);
//      // Run the contour tree on the mesh
//      RunContourTree(fieldArray,
//                     contourTree,
//                     sortOrder,
//                     nIterations,
//                     mesh,
//                     computeRegularStructure,
//                     mesh.GetMeshBoundaryExecutionObject());


      // Writing to a file
      std::ofstream file1("doubleASCII-conn.txt");
      std::ofstream file2("doubleASCII-offs.txt");

      std::cout << "std_nbor_connectivity\n";
      for(unsigned i = 0; i < delmesh.std_nbor_connectivity.size(); i++)
      {
          if (file1.is_open()) file1 << i << " " << delmesh.std_nbor_connectivity[i] << "\n";
//          std::cout << i << " " << delmesh.std_nbor_connectivity[i] << "\n";
      }
      std::cout << "std_nbor_offsets\n";
      for(unsigned i = 0; i < delmesh.std_nbor_offsets.size(); i++)
      {
          if (file2.is_open()) file2 << i << " " << delmesh.std_nbor_offsets[i] << "\n";
//          std::cout << i << " " << delmesh.std_nbor_offsets[i] << "\n";
      }


      return;
    }

//    // 3D Contour Tree with Freudenthal
//    else
//    {
//      std::cout << "USING 3D Freudenthal...\n";
//      // Build the mesh and fill in the values
//      DataSetMeshTriangulation3DFreudenthal mesh(meshSize);
//      // Run the contour tree on the mesh
//      RunContourTree(fieldArray,
//                     contourTree,
//                     sortOrder,
//                     nIterations,
//                     mesh,
//                     computeRegularStructure,
//                     mesh.GetMeshBoundaryExecutionObject());
//      return;
//    }

  }


private:
  /*!
  *  Run the contour tree for the given mesh. This function implements the main steps for
  *  computing the contour tree after the mesh has been constructed using the approbrite
  *  contour tree mesh class.
  *
  *  fieldArray   : The values of the mesh
  *  contourTree  : The output contour tree to be computed (output)
  *  sortOrder    : The sort order for the mesh vertices (output)
  *  nIterations  : The number of iterations used to compute the contour tree (output)
  *  mesh : The specific mesh (see vtkm/worklet/contourtree_augmented/mesh_dem_meshtypes
  *  computeRegularStructure : 0=Off, 1=full augmentation with all vertices
  *                            2=boundary augmentation using meshBoundary
  *  meshBoundary : This parameter is generated by calling mesh.GetMeshBoundaryExecutionObject
  *                 For regular 2D/3D meshes this required no extra parameters, however, for a
  *                 ContourTreeMesh additional information about the block must be given. Rather
  *                 than generating the MeshBoundary descriptor here, we therefore, require it
  *                 as an input. The MeshBoundary is used to augment the contour tree with the
  *                 mesh boundary vertices. It is needed only if we want to augement by the
  *                 mesh boundary and computeRegularStructure is False (i.e., if we compute
  *                 the full regular strucuture this is not needed because all vertices
  *                 (including the boundary) will be addded to the tree anyways.
  */
  template <typename FieldType,
            typename StorageType,
            typename MeshClass,
            typename MeshBoundaryClass>
  void RunContourTree(const vtkm::cont::ArrayHandle<FieldType, StorageType> fieldArray,
                      contourtree_augmented::ContourTree& contourTree,
                      contourtree_augmented::IdArrayType& sortOrder,
                      vtkm::Id& nIterations,
                      MeshClass& mesh,
                      unsigned int computeRegularStructure,
                      const MeshBoundaryClass& meshBoundary)
  {
    using namespace vtkm::worklet::contourtree_augmented;
    // Stage 1: Load the data into the mesh. This is done in the Run() method above and accessible
    //          here via the mesh parameter. The actual data load is performed outside of the
    //          worklet in the example contour tree app (or whoever uses the worklet)
    /// DEBUG PRINT std::cout << "S1. {worklet/ContourTreeUniformAugmented.h : RunContourTree}\n";
    // Stage 2 : Sort the data on the mesh to initialize sortIndex & indexReverse on the mesh
    // Start the timer for the mesh sort
    /// DEBUG PRINT std::cout << "S2\n";
    vtkm::cont::Timer timer;
    timer.Start();
    std::stringstream timingsStream; // Use a string stream to log in one message

    // Sort the mesh data
    mesh.SortData(fieldArray);
    timingsStream << "    " << std::setw(38) << std::left << "Sort Data"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
    timer.Start();

    // Stage 3: Assign every mesh vertex to a peak
    /// DEBUG PRINT std::cout << "S3\n";
    MeshExtrema extrema(mesh.NumVertices);
    extrema.SetStarts(mesh, true);
    extrema.BuildRegularChains(true);
    timingsStream << "    " << std::setw(38) << std::left << "Join Tree Regular Chains"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
    timer.Start();

    /// DEBUG PRINT std::cout << "Join Tree Regular Chains | Printing the extremums per vertex ...\n";
//    for(unsigned i = 0; i < mesh.NumVertices; i++)
//    {
////        std::cout << "i=" << i << MaskedIndex(extrema.Peaks[i]) << ", " <<  extrema.Peaks[i];
//        extrema.DebugPrint();
//    }
    const char* message = "debug";
    extrema.DebugPrint(message, message, 100);

    // Stage 4: Identify join saddles & construct Active Join Graph
    /// DEBUG PRINT std::cout << "S4\n";
    MergeTree joinTree(mesh.NumVertices, true);
    /// DEBUG PRINT std::cout << "join tree made with: " << mesh.NumVertices << "\n";
    ActiveGraph joinGraph(true);
    /// DEBUG PRINT std::cout << "ActiveGraph joinGraph(true) made\n";
    joinGraph.Initialise(mesh, extrema);
    /// DEBUG PRINT std::cout << "joinGraph initialised\n";
    timingsStream << "    " << std::setw(38) << std::left << "Join Tree Initialize Active Graph"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;

#ifdef DEBUG_PRINT
    joinGraph.DebugPrint("Active Graph Instantiated", __FILE__, __LINE__);
#endif
    timer.Start();

    // Stage 5: Compute Join Tree Hyperarcs from Active Join Graph
    /// DEBUG PRINT std::cout << "S5\n";
    joinGraph.MakeMergeTree(joinTree, extrema);
    timingsStream << "    " << std::setw(38) << std::left << "Join Tree Compute"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
#ifdef DEBUG_PRINT
    joinTree.DebugPrint("Join tree Computed", __FILE__, __LINE__);
    joinTree.DebugPrintTree("Join tree", __FILE__, __LINE__, mesh);
#endif
    timer.Start();

    // Stage 6: Assign every mesh vertex to a pit
    /// DEBUG PRINT std::cout << "S6\n";
    extrema.SetStarts(mesh, false);
    extrema.BuildRegularChains(false);
    timingsStream << "    " << std::setw(38) << std::left << "Split Tree Regular Chains"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
    timer.Start();

    #if PACT_DEBUG
        std::cout << "Split Tree Regular Chains | Printing the extremums per vertex ...\n";
    #endif

//    for(unsigned i = 0; i < mesh.NumVertices; i++)
//    {
////        std::cout << "i=" << i << MaskedIndex(extrema.Peaks[i]) << ", " <<  extrema.Peaks[i];
//        extrema.DebugPrint();
//    }
//    const char* message = "debug";
    #if PACT_DEBUG
        extrema.DebugPrint(message, message, 100);
    #endif

    // Stage 7:     Identify split saddles & construct Active Split Graph
    /// DEBUG PRINT std::cout << "S7\n";
    MergeTree splitTree(mesh.NumVertices, false);
    ActiveGraph splitGraph(false);
    splitGraph.Initialise(mesh, extrema);
    timingsStream << "    " << std::setw(38) << std::left << "Split Tree Initialize Active Graph"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
#ifdef DEBUG_PRINT
    splitGraph.DebugPrint("Active Graph Instantiated", __FILE__, __LINE__);
#endif
    timer.Start();

    // Stage 8: Compute Split Tree Hyperarcs from Active Split Graph
    /// DEBUG PRINT std::cout << "S8\n";
    splitGraph.MakeMergeTree(splitTree, extrema);
    timingsStream << "    " << std::setw(38) << std::left << "Split Tree Compute"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
#ifdef DEBUG_PRINT
    splitTree.DebugPrint("Split tree Computed", __FILE__, __LINE__);
    // Debug split and join tree
    joinTree.DebugPrintTree("Join tree", __FILE__, __LINE__, mesh);
    splitTree.DebugPrintTree("Split tree", __FILE__, __LINE__, mesh);
#endif
    timer.Start();

    // Stage 9: Join & Split Tree are Augmented, then combined to construct Contour Tree
    /// DEBUG PRINT std::cout << "S9\n";
    contourTree.Init(mesh.NumVertices);
    ContourTreeMaker treeMaker(contourTree, joinTree, splitTree);
    // 9.1 First we compute the hyper- and super- structure
    treeMaker.ComputeHyperAndSuperStructure();
    timingsStream << "    " << std::setw(38) << std::left
                  << "Contour Tree Hyper and Super Structure"
                  << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
    timer.Start();

    // 9.2 Then we compute the regular structure
    if (computeRegularStructure == 1) // augment with all vertices
    {
      treeMaker.ComputeRegularStructure(extrema);
      timingsStream << "    " << std::setw(38) << std::left << "Contour Tree Regular Structure"
                    << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
    }
    else if (computeRegularStructure == 2) // augment by the mesh boundary
    {
	  std::cout << "computeRegularStructure with dummy meshBoundary ...\n";
      treeMaker.ComputeBoundaryRegularStructure(extrema, mesh, meshBoundary);
      timingsStream << "    " << std::setw(38) << std::left
                    << "Contour Tree Boundary Regular Structure"
                    << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
    }
    timer.Start();

    // Collect the output data
    nIterations = treeMaker.ContourTreeResult.NumIterations;
    //  Need to make a copy of sortOrder since ContourTreeMesh uses a smart array handle
    // TODO: Check if we can just make sortOrder a return array with variable type or if we can make the SortOrder return optional
    // TODO/FIXME: According to Ken Moreland the short answer is no. We may need to go back and refactor this when we
    // improve the contour tree API. https://gitlab.kitware.com/vtk/vtk-m/-/merge_requests/2263#note_831128 for more details.
    vtkm::cont::Algorithm::Copy(mesh.SortOrder, sortOrder);
    // ProcessContourTree::CollectSortedSuperarcs<DeviceAdapter>(contourTree, mesh.SortOrder, saddlePeak);
    // contourTree.SortedArcPrint(mesh.SortOrder);
    // contourTree.PrintDotSuperStructure();

    // Log the collected timing results in one coherent log entry
    this->TimingsLogString = timingsStream.str();
    if (this->TimingsLogLevel != vtkm::cont::LogLevel::Off)
    {
      VTKM_LOG_S(this->TimingsLogLevel,
                 std::endl
                   << "    ------------------- Contour Tree Worklet Timings ----------------------"
                   << std::endl
                   << this->TimingsLogString);
    }
  }
};

} // namespace vtkm
} // namespace vtkm::worklet

#endif // vtk_m_worklet_ContourTreeUniformAugmented_h
