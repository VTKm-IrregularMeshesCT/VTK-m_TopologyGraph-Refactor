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

#include <vtkm/filter/scalar_topology/worklet/contourtree_distributed/PrintGraph.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

//#include <string>     // std::string, std::stof

// for sleeping
#include <chrono>
#include <thread>

// for memory usage
#include <sys/resource.h>
#include <unistd.h>

#define PACT_DEBUG 0
#define WRITE_FILES 1

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


void static printMemoryUsage(const std::string& message)
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << message << " - Memory usage: " << usage.ru_maxrss << " KB" << std::endl;
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
      std::cout << "CTUA in another Run() before RunContourTree()" << std::endl;

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
             vtkm::worklet::contourtree_augmented::ContourTreeMesh<int>& mesh,
             contourtree_augmented::ContourTree& contourTree,
             contourtree_augmented::IdArrayType& sortOrder,
             vtkm::Id& nIterations,
             unsigned int computeRegularStructure = 1)
    {
        std::cout << "worklet/ContourTreeUniformAugmented::RunContourTree() Run() [PACT For Irregular Data]" << std::endl;
        std::cout << "else: worklet/ContourTreeUniformAugmented::RunContourTree() ..." << std::endl;


          // Run the contour tree on the mesh
          // PACT:
          RunContourTree(fieldArray, //fakeFieldArray, //fieldArray,
                         contourTree,
                         sortOrder,
                         nIterations,
                         mesh,
                         computeRegularStructure,
                         mesh.GetMeshBoundaryExecutionObject());

#if WRITE_FILES
          std::ofstream outFile("CT-full-superdot.gv");

          vtkm::Id detailedMask =   vtkm::worklet::contourtree_distributed::SHOW_SUPER_STRUCTURE \
                                  | vtkm::worklet::contourtree_distributed::SHOW_SUPERNODE_ID \
                                  | vtkm::worklet::contourtree_distributed::SHOW_SUPERARC_ID \
                                  | vtkm::worklet::contourtree_distributed::SHOW_MESH_SORT_ID;
    //                              | vtkm::worklet::contourtree_distributed::SHOW_SUPERPARENT \
    //                              | vtkm::worklet::contourtree_distributed::SHOW_ITERATION \
    //                              | vtkm::worklet::contourtree_distributed::SHOW_DATA_VALUE \
    //                              | vtkm::worklet::contourtree_distributed::SHOW_HYPER_STRUCTURE \
    //                              | vtkm::worklet::contourtree_distributed::SHOW_ALL_IDS \
    //                              | vtkm::worklet::contourtree_distributed::SHOW_ALL_HYPERIDS;


          // Call the function after you've computed ContourTree and your associated data structures (`mesh` and `field`):
          outFile << vtkm::worklet::contourtree_distributed::ContourTreeDotGraphPrintSerial(
              "Contour Tree Super Dot",         // label/title
              mesh,                             // mesh (re)constructed above
              fieldArray, //fakeFieldArray, //fieldArray,     // scalar data array handle
              contourTree,                      // computed contour tree structure
              detailedMask,                     // detailed output with all info
              vtkm::cont::ArrayHandle<vtkm::Id>()); // global ids


          outFile.close();
#endif

          return;
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
      using namespace vtkm::worklet::contourtree_augmented;
      // 2D Contour Tree
      if (meshSize[2] == 1)
      {
        // Build the mesh and fill in the values
        DataSetMeshTriangulation2DFreudenthal mesh(vtkm::Id2{ meshSize[0], meshSize[1] });
        // Run the contour tree on the mesh
        RunContourTree(fieldArray,
                       contourTree,
                       sortOrder,
                       nIterations,
                       mesh,
                       computeRegularStructure,
                       mesh.GetMeshBoundaryExecutionObject());
        return;
      }
      // 3D Contour Tree using marching cubes
      else if (useMarchingCubes)
      {
        // Build the mesh and fill in the values
        DataSetMeshTriangulation3DMarchingCubes mesh(meshSize);
        // Run the contour tree on the mesh
        RunContourTree(fieldArray,
                       contourTree,
                       sortOrder,
                       nIterations,
                       mesh,
                       computeRegularStructure,
                       mesh.GetMeshBoundaryExecutionObject());
        return;
      }
      // 3D Contour Tree with Freudenthal
      else
      {
        // Build the mesh and fill in the values
        DataSetMeshTriangulation3DFreudenthal mesh(meshSize);
        // Run the contour tree on the mesh
        RunContourTree(fieldArray,
                       contourTree,
                       sortOrder,
                       nIterations,
                       mesh,
                       computeRegularStructure,
                       mesh.GetMeshBoundaryExecutionObject());
        return;
      }
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
