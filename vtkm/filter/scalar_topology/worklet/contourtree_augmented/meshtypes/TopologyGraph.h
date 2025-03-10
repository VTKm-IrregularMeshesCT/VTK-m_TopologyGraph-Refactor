//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
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

#ifndef vtk_m_worklet_contourtree_augmented_mesh_dem_topology_graph_h
#define vtk_m_worklet_contourtree_augmented_mesh_dem_topology_graph_h

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayHandlePermutation.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/ArrayRangeCompute.h>
#include <vtkm/cont/ArrayRangeComputeTemplate.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ArrayTransforms.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/data_set_mesh/IdRelabeler.h> // This is needed only as an unused default argument.
//#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/MeshStructureContourTreeMesh.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/MeshStructureTopologyGraph.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/AddToArrayElementsWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/ApplyLookupTableDecorator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/ArcComparator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/ArcValidDecorator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/CombinedSimulatedSimplicityIndexComparator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/CombinedVectorDifferentFromNext.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/CopyIntoCombinedArrayWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/CopyIntoCombinedNeighborsWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/CopyNeighborsToPackedArray.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/FindDuplicateInOtherWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/GetArcFromDecorator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/InitToCombinedSortOrderArraysWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/MergeSortedListsWithoutDuplicatesWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/contourtreemesh/ReplaceArcNumWithToVertexWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/ComputeMeshBoundaryContourTreeMesh.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/MeshBoundaryContourTreeMesh.h>
#include <vtkm/io/ErrorIO.h>
#include <vtkm/worklet/DispatcherMapField.h>

#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/PrintVectors.h> // TODO remove should not be needed

#include <vtkm/cont/ExecutionObjectBase.h>

#include <vtkm/cont/Invoker.h>

#include <algorithm>

VTKM_THIRDPARTY_PRE_INCLUDE
#include <vtkm/thirdparty/diy/Configure.h>
#include <vtkm/thirdparty/diy/diy.h>
VTKM_THIRDPARTY_POST_INCLUDE
// clang-format on

namespace topology_graph_inc_ns =
        vtkm::worklet::contourtree_augmented::mesh_dem_topology_graph_inc;
//  vtkm::worklet::contourtree_augmented::mesh_dem_topology_graph_inc;

// #define DEBUG_PRINT

namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{
template <typename FieldType>
class TopologyGraph : public vtkm::cont::ExecutionObjectBase
{ // class TopologyGraph
public:
  //Mesh dependent helper functions
  void SetPrepareForExecutionBehavior(bool getMax);

  topology_graph_inc_ns::MeshStructureTopologyGraph PrepareForExecution(
    vtkm::cont::DeviceAdapterId,
    vtkm::cont::Token& token) const;

  // ported the default contour tree mesh constructor ...
  // ... to actually construct a valid CT, ...
  // ... rather than just 'ContourTreeMesh() {}'
  TopologyGraph() {}

  // Constructor 0 - The TOPOLOGY GRAPH:
  TopologyGraph(const IdArrayType& nodes,
                //const IdArrayType& arcs, forget the arcs, ...
                // ... we do them manually for the topology graph:
                const IdArrayType& inNborConnectivity,
                const IdArrayType& inNborOffsets,
                // (arcs were used to compute nbor connectivity and offsets ...
                //  ... but because topology graphs no longer have ...
                //  ... the outdegree strictly capped at 1, ...
                //  ... we enter the information manually - for now anyway)
                const IdArrayType& inSortOrder,
                const vtkm::cont::ArrayHandle<FieldType>& values,
                const IdArrayType& inGlobalMeshIndex);

  vtkm::Id GetNumberOfVertices() const { return this->NumVertices; }

  // Combine two ContourTreeMeshes
  void MergeWith(TopologyGraph<FieldType>& other,
                 vtkm::cont::LogLevel TreeLogLevel = vtkm::cont::LogLevel::Perf,
                 std::string timingsMessage = "");

  // Save/Load the mesh helpers
  void Save(const char* filename) const;
  void Load(const char* filename);

  // Empty placeholder function to ensure compliance of this class with the interface
  // the other mesh classes. This is a no-op here since this class is initalized
  // from a known contour tree so sort is already done
  template <typename T, typename StorageType>
  void SortData(const vtkm::cont::ArrayHandle<T, StorageType>& values) const
  {
    (void)values; // Do nothing but avoid unsused param warning
  }

  // Public fields
  // WARNING! WARNING! WARNING!
  // MAX_OUTDEGREE WAS SET TO 20
  // ... NO LONGER THE CASE FOR TOPOLOGY GRAPHS!!!
  static const int MAX_OUTDEGREE = 150;// 20;

  vtkm::Id NumVertices;
  vtkm::cont::ArrayHandleIndex SortOrder;
  vtkm::cont::ArrayHandleIndex SortIndices;
  vtkm::cont::ArrayHandle<FieldType> SortedValues;
  IdArrayType GlobalMeshIndex;
  // NeighborConnectivity stores for each vertex the indices of its neighbors. For each vertex
  // the indices are sorted by value, i.e, the first neighbour has the lowest and
  // the last neighbour the highest value for the vertex. In the array we just
  // concatinate the list of neighbors from all vertices, i.e., we first
  // have the list of neighbors of the first vertex, then the second vertex and so on, i.e.:
  // [ n_1_1, n_1_2, n_2_1, n_2_2, n_2_3, etc.]
  IdArrayType NeighborConnectivity;
  // NeighborOffsets gives us for each vertex an index into the neighours array indicating
  // the index where the list of neighbors for the vertex begins
  IdArrayType NeighborOffsets;
  // the maximum number of neighbors of a vertex
  vtkm::Id MaxNeighbors;

  // Print Contents
  void PrintContent(std::ostream& outStream = std::cout) const;

  // Debug print routine
  void DebugPrint(const char* message, const char* fileName, long lineNum) const;

  // Get boundary execution object
  MeshBoundaryContourTreeMeshExec GetMeshBoundaryExecutionObject(vtkm::Id3 globalSize,
                                                                 vtkm::Id3 minIdx,
                                                                 vtkm::Id3 maxIdx) const;
  MeshBoundaryContourTreeMeshExec GetMeshBoundaryExecutionObject() const;

  void GetBoundaryVertices(IdArrayType& boundaryVertexArray,                    // output
                           IdArrayType& boundarySortIndexArray,                 // output
                           MeshBoundaryContourTreeMeshExec* meshBoundaryExecObj //input
  ) const;

  /// copies the global IDs for a set of sort IDs
  /// notice that the sort ID is the same as the mesh ID for the ContourTreeMesh class.
  /// To reduce memory usage we here use a fancy array handle rather than copy data
  /// as is needed for the DataSetMesh types.
  /// We here return a fancy array handle to convert values on-the-fly without requiring additional memory
  /// @param[in] sortIds Array with sort Ids to be converted from local to global Ids
  /// @param[in] localToGlobalIdRelabeler This parameter is here only for
  ///            consistency with the DataSetMesh types but is not
  ///            used here and as such can simply be set to nullptr
  inline vtkm::cont::ArrayHandlePermutation<IdArrayType, IdArrayType> GetGlobalIdsFromSortIndices(
    const IdArrayType& sortIds,
    const vtkm::worklet::contourtree_augmented::mesh_dem::IdRelabeler* localToGlobalIdRelabeler =
      nullptr) const
  {                                 // GetGlobalIDsFromSortIndices()
    (void)localToGlobalIdRelabeler; // avoid compiler warning
    return vtkm::cont::make_ArrayHandlePermutation(sortIds, this->GlobalMeshIndex);
  } // GetGlobalIDsFromSortIndices()

  /// copies the global IDs for a set of mesh IDs
  /// notice that the sort ID is the same as the mesh ID for the ContourTreeMesh class.
  /// To reduce memory usage we here use a fancy array handle rather than copy data
  /// as is needed for the DataSetMesh types.
  /// MeshIdArrayType must be an array if Ids. Usually this is a vtkm::worklet::contourtree_augmented::IdArrayType
  /// but in some cases it may also be a fancy array to avoid memory allocation
  /// We here return a fancy array handle to convert values on-the-fly without requiring additional memory
  /// @param[in] meshIds Array with mesh Ids to be converted from local to global Ids
  /// @param[in] localToGlobalIdRelabeler This parameter is here only for
  ///            consistency with the DataSetMesh types but is not
  ///            used here and as such can simply be set to nullptr
  template <typename MeshIdArrayType>
  inline vtkm::cont::ArrayHandlePermutation<MeshIdArrayType, IdArrayType>
  GetGlobalIdsFromMeshIndices(const MeshIdArrayType& meshIds,
                              const vtkm::worklet::contourtree_augmented::mesh_dem::IdRelabeler*
                                localToGlobalIdRelabeler = nullptr) const
  {                                 // GetGlobalIDsFromMeshIndices()
    (void)localToGlobalIdRelabeler; // avoid compiler warning
    return vtkm::cont::make_ArrayHandlePermutation(meshIds, this->GlobalMeshIndex);
  } // GetGlobalIDsFromMeshIndices()

private:
  vtkm::cont::Invoker Invoke;

  bool mGetMax; // Define the behavior for the PrepareForExecution function

  // Private init and helper functions
  void InitializeNeighborConnectivityFromArcs(const IdArrayType& arcs);
  void ComputeMaxNeighbors();

  // Private helper functions for saving data vectors
  // Internal helper function to save 1D index array to file
  template <typename ValueType>
  void SaveVector(std::ostream& os, const vtkm::cont::ArrayHandle<ValueType>& vec) const;

  // Internal helper function to Load 1D index array from file
  template <typename ValueType>
  void LoadVector(std::istream& is, vtkm::cont::ArrayHandle<ValueType>& vec);

}; // ContourTreeMesh

// print content
template <typename FieldType>
inline void TopologyGraph<FieldType>::PrintContent(std::ostream& outStream /*= std::cout*/) const
{ // PrintContent()
  PrintHeader(this->NumVertices, outStream);
  //PrintIndices("SortOrder", SortOrder, outStream);
  PrintValues("SortedValues", this->SortedValues, -1, outStream);
  PrintIndices("GlobalMeshIndex", this->GlobalMeshIndex, -1, outStream);
  PrintIndices("NeighborConnectivity", this->NeighborConnectivity, -1, outStream);
  PrintIndices("NeighborOffsets", this->NeighborOffsets, -1, outStream);
  outStream << "MaxNeighbors=" << this->MaxNeighbors << std::endl;
  outStream << "mGetMax=" << this->mGetMax << std::endl;
} // PrintContent()

// debug routine
template <typename FieldType>
inline void TopologyGraph<FieldType>::DebugPrint(const char* message,
                                                   const char* fileName,
                                                   long lineNum) const
{ // DebugPrint()
  std::cout << "---------------------------" << std::endl;
  std::cout << std::setw(30) << std::left << fileName << ":" << std::right << std::setw(4)
            << lineNum << std::endl;
  std::cout << std::left << std::string(message) << std::endl;
  std::cout << "Contour Tree Mesh Contains:     " << std::endl;
  std::cout << "---------------------------" << std::endl;
  std::cout << std::endl;

  PrintContent(std::cout);
} // DebugPrint()

// create the contour tree mesh from contour tree data
//template <typename FieldType>
//TopologyGraph<FieldType>::TopologyGraph(const IdArrayType& arcs,
//                                            const IdArrayType& inSortOrder,
//                                            const vtkm::cont::ArrayHandle<FieldType>& values,
//                                            const IdArrayType& inGlobalMeshIndex)
//  : SortOrder()
//  , SortedValues()
//  , GlobalMeshIndex(inGlobalMeshIndex)
//  , NeighborConnectivity()
//  , NeighborOffsets()
//{
//  std::cout << "{ContourTreeMesh.h : Constructor A\n";
//  std::cout << " ... arcs, inSortOrder, values, inGlobalMeshIndex}\n";
//  this->NumVertices = inSortOrder.GetNumberOfValues();
//  // Initalize the SortedIndices as a smart array handle
//  this->SortIndices = vtkm::cont::ArrayHandleIndex(this->NumVertices);
//  this->SortOrder = vtkm::cont::ArrayHandleIndex(this->NumVertices);
//  // values permuted by SortOrder to sort the values
//  auto permutedValues = vtkm::cont::make_ArrayHandlePermutation(inSortOrder, values);
//  // TODO check if we actually need to make this copy here. we could just store the permutedValues array to save memory
//  vtkm::cont::Algorithm::Copy(permutedValues, this->SortedValues);
//  this->InitializeNeighborConnectivityFromArcs(arcs);
//#ifdef DEBUG_PRINT
//  // Print the contents fo this for debugging
//  DebugPrint("ContourTreeMesh Initialized", __FILE__, __LINE__);
//#endif
//}


template <typename FieldType>
inline TopologyGraph<FieldType>::TopologyGraph(const IdArrayType& nodes,
                                                   //const IdArrayType& arcs, forget the arcs, ...
                                                   // ... we do them manually for the topology graph:
                                                   const IdArrayType& inNborConnectivity,
                                                   const IdArrayType& inNborOffsets,
                                                   // (arcs were used to compute nbor connectivity and offsets ...
                                                   //  ... but because topology graphs no longer have ...
                                                   //  ... the outdegree strictly capped at 1, ...
                                                   //  ... we enter the information manually - for now anyway)
                                                   const IdArrayType& inSortOrder,
                                                   const vtkm::cont::ArrayHandle<FieldType>& values,
                                                   const IdArrayType& inGlobalMeshIndex)
  : GlobalMeshIndex(inGlobalMeshIndex)
  , NeighborConnectivity(inNborConnectivity)
  , NeighborOffsets(inNborOffsets)
{
  std::cout << "{ContourTreeMesh.h : Constructor 0\n";
  std::cout << " ... +nodes, arcs, inSortOrder, values, inGlobalMeshIndex}\n";

  PrintIndices("inputNeighborConnectivity", inNborConnectivity, -1, std::cout);
  PrintIndices("NeighborConnectivity", this->NeighborConnectivity, -1, std::cout);
  PrintIndices("inputNeighborOffsets", inNborOffsets, -1, std::cout);
  PrintIndices("NeighborOffsets", this->NeighborOffsets, -1, std::cout);

  // Initialize the SortedValues array with values permuted by the SortOrder permuted by the nodes, i.e.,
  // this->SortedValues[v] = values[inSortOrder[nodes[v]]];
  vtkm::cont::ArrayHandlePermutation<IdArrayType, IdArrayType> permutedSortOrder(nodes,
                                                                                 inSortOrder);
  auto permutedValues = vtkm::cont::make_ArrayHandlePermutation(permutedSortOrder, values);
  vtkm::cont::Algorithm::Copy(permutedValues, this->SortedValues);
  // Initalize the SortedIndices as a smart array handle
  this->NumVertices = this->SortedValues.GetNumberOfValues();
  this->SortIndices = vtkm::cont::ArrayHandleIndex(this->NumVertices);
  this->SortOrder = vtkm::cont::ArrayHandleIndex(this->NumVertices);
  //this->InitializeNeighborConnectivityFromArcs(arcs); we do not do this here ...
  // nbor connectivity/offsets are now set in the member initializer list above

  // Still ... since we no longer call InitializeNeighborConnectivityFromArcs ...
  // ... we must not forget to compute maximum number of neighbors, ...
  // ... which we can here - since we have the nbor connection/offset information:
  this->ComputeMaxNeighbors();


#ifdef DEBUG_PRINT
  // Print the contents fo this for debugging
  DebugPrint("ContourTreeMesh Initialized", __FILE__, __LINE__);
#endif
}



//// Helper function to copy select set of indices of source array into
//// select set of inidces of destination array. Important: srcIndices
//// and  dstIndices must contain the same number of values.
//template <typename PT1, typename PT2, typename PT3, typename PT4>
//inline void CopyArrayByIndices(const PT1& srcArray,
//                               const PT2& srcIndices,
//                               PT3& dstArray,
//                               const PT4& dstIndices)
//{
//  VTKM_ASSERT(srcIndices.GetNumberOfValues() == dstIndices.GetNumberOfValues());
//  auto srcPermutation = make_ArrayHandlePermutation(srcIndices, srcArray);
//  auto dstPermuation = make_ArrayHandlePermutation(dstIndices, dstArray);
//  vtkm::cont::Algorithm::Copy(srcPermutation, dstPermuation);
//}

//// Helper function doing the same as previous function, but for
//// arrays of vector. This is necessary since we use an array
//// created with ArrayHandleGroupVecVariable as destination,
//// which breaks some conventions of ArrayHandle and does not work
//// with vtkm::cont::Algorithm::Copy.
//template <typename PT1, typename PT2, typename PT3, typename PT4>
//inline void CopyVecArrayByIndices(const PT1& srcArray,
//                                  const PT2& srcIndices,
//                                  PT3& dstArray,
//                                  const PT4& dstIndices)
//{
//  VTKM_ASSERT(srcIndices.GetNumberOfValues() == dstIndices.GetNumberOfValues());
//  auto srcPermutation = make_ArrayHandlePermutation(srcIndices, srcArray);
//  auto dstPermuation = make_ArrayHandlePermutation(dstIndices, dstArray);
//  // Use a worklet for copying data since ArrayHandleGroupVecVariable does
//  // not work as destination for vtkm::cont::Algorithm::Copy.
//  vtkm::cont::Invoker invoke;
//  invoke(
//    contourtree_mesh_inc_ns::CopyIntoCombinedNeighborsWorklet{}, srcPermutation, dstPermuation);
//}

// Initalize the contour tree from the arcs array
template <typename FieldType>
inline void TopologyGraph<FieldType>::InitializeNeighborConnectivityFromArcs(
  const IdArrayType& arcs)
{
  // This function computes the neighbor connecitvity (NeighborConnectivity, NeighborOffsets) from
  // an arc array. An arc array consists of undirected arcs. arc[i] connects contour tree nodes
  // i and arc[i]. For the neighbor connectiviy in the contour tree mesh, we first convert these
  // into two directed arcs that are then used to compute a list of neighbors in the mesh for
  // each node.
  //
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
  // (This is a graph with nodes 0 through 4 and edges 0 -> 1, 2 -> 3, 3 -> 1, 4 -> 3).
  // The way the graph is structured, each nodes has at most one outgoing edge.
  // The contour tree algorithm stores this in an arcs array:
  //
  // idx:  0 1 2 3 4
  // arcs: 1 - 3 1 3 (- = NO_SUCH_ELEMENT, meaning no arc originating from this node)
  //
  // This function translates this into the internal contour tree mesh represnetation,
  // which is "regular" vtk-m connectivity format, i.e., the connectity array is a
  // flat list of neighbor vertices and offsets give the start index of the
  // neighbor list for each vertex:
  //
  // connectivity: 1 0 3 3 1 2 4 3
  // counts: 1 2 1 3 1
  // offset: 0 1 3 4 7 8

  // Step 1: Implictely view arc array as directed arcs and add arcs in the opposite
  // direction. In the resulting arc list, arc 2*idx is the arc idx->arcs[idx] and arc
  // 2*idx+1 is the arc arcs[idx]->idx, i.e., in our example,
  // idx:  0 1 2 3 4 5 6 7 8 9
  // from: 0 1 1 - 2 3 3 1 4 3
  // to:   1 0 - 1 3 2 1 3 3 4
  vtkm::Id nArcsTotal = 2 * arcs.GetNumberOfValues();
  vtkm::cont::ArrayHandleIndex indexArray(nArcsTotal);
  auto arcIsValidArray = make_ArrayHandleDecorator(
    nArcsTotal, mesh_dem_contourtree_mesh_inc::ArcValidDecoratorImpl{}, arcs);
  // We first generate a list of "valid" arcs in this->NeighborConnectivity, in our
  // example:
  // connectivity: 0 1 4 5 6 7 8 9
  vtkm::cont::Algorithm::CopyIf(indexArray, arcIsValidArray, this->NeighborConnectivity);
  vtkm::Id nValidArcs = this->NeighborConnectivity.GetNumberOfValues();

  // Step 2: Sort arcs---by permuting their indices in the connectiviy array---so
  // that all arcs originating at the same vertex (same `from`) are adjacent.
  // All arcs are in neighbors array based on  sort index of their 'from' vertex
  // (and then within a run sorted by sort index of their 'to' vertex).
  // In our example this results in:
  // connectivity: 0 1 7 4 6 5 9 8
  // corresponding to an arc order of
  // from: 0 1 1 2 3 3 3 4
  // to:   1 0 3 3 1 2 4 3
  vtkm::cont::Algorithm::Sort(this->NeighborConnectivity,
                              contourtree_mesh_inc_ns::ArcComparator(arcs));

  // We can now obtain counts of the connectivity array by counting the number
  // of consecutive `from` entries with the same value. In our example:
  // counts: 1 2 1 3 1
  auto arcFrom = make_ArrayHandleDecorator(nValidArcs,
                                           mesh_dem_contourtree_mesh_inc::GetArcFromDecoratorImpl{},
                                           this->NeighborConnectivity,
                                           arcs);
  auto constOne = vtkm::cont::make_ArrayHandleConstant(vtkm::Id{ 1 }, nValidArcs);
  vtkm::cont::ArrayHandle<vtkm::Id> uniqueKeys;
  vtkm::cont::ArrayHandle<vtkm::Id> counts;
  vtkm::cont::Algorithm::ReduceByKey(arcFrom, constOne, uniqueKeys, counts, vtkm::Add());
  VTKM_ASSERT(uniqueKeys.GetNumberOfValues() == this->NumVertices);

  // Convert counts into offsts for the connectivity array
  vtkm::Id neighborOffsetsSize;
  vtkm::cont::ConvertNumComponentsToOffsets(counts, this->NeighborOffsets, neighborOffsetsSize);

  // Finally, the correct connectivity array correspons to the `to` array,
  // so replace arc indices with its `to`vertex. In our example, this results in:
  // connectivity: 1 0 3 3 1 2 4 3
  // which is exactly the array we needed to compute
  contourtree_mesh_inc_ns::ReplaceArcNumWithToVertexWorklet replaceArcNumWithToVertexWorklet;
  this->Invoke(replaceArcNumWithToVertexWorklet,
               this->NeighborConnectivity, // input/output
               arcs                        // input
  );

  // Compute maximum number of neighbors
  this->ComputeMaxNeighbors();

#ifdef DEBUG_PRINT
  std::cout << std::setw(30) << std::left << __FILE__ << ":" << std::right << std::setw(4)
            << __LINE__ << std::endl;
  auto neighborOffsetPortal = this->NeighborOffsets.ReadPortal();
  auto neighborConnectivityPortal = this->NeighborConnectivity.ReadPortal();
  for (vtkm::Id vtx = 0; vtx < NeighborOffsets.GetNumberOfValues(); ++vtx)
  {
    std::cout << vtx << ": ";
    vtkm::Id neighboursBeginIndex = neighborOffsetPortal.Get(vtx);
    vtkm::Id neighboursEndIndex = (vtx < this->NumVertices - 1)
      ? neighborOffsetPortal.Get(vtx + 1)
      : NeighborConnectivity.GetNumberOfValues();

    for (vtkm::Id ni = neighboursBeginIndex; ni < neighboursEndIndex; ++ni)
    {
      std::cout << neighborConnectivityPortal.Get(ni) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Max neighbours: " << this->MaxNeighbors << std::endl;
#endif
}

template <typename FieldType>
inline void TopologyGraph<FieldType>::ComputeMaxNeighbors()
{
  auto neighborCounts = make_ArrayHandleOffsetsToNumComponents(this->NeighborOffsets);

  vtkm::cont::ArrayHandle<vtkm::Range> rangeArray =
    vtkm::cont::ArrayRangeComputeTemplate(neighborCounts);
  this->MaxNeighbors = static_cast<vtkm::Id>(rangeArray.ReadPortal().Get(0).Max);
}

// Define the behavior for the execution object generate by the PrepareForExecution function
template <typename FieldType>
inline void TopologyGraph<FieldType>::SetPrepareForExecutionBehavior(bool getMax)
{
  this->mGetMax = getMax;
}

// Get VTKM execution object that represents the structure of the mesh and provides the mesh helper functions on the device
template <typename FieldType>
topology_graph_inc_ns::MeshStructureTopologyGraph inline TopologyGraph<
  FieldType>::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                  vtkm::cont::Token& token) const
{
  return topology_graph_inc_ns::MeshStructureTopologyGraph(this->NeighborConnectivity,
                                                               this->NeighborOffsets,
                                                               this->MaxNeighbors,
                                                               this->mGetMax,
                                                               device,
                                                               token);
}

//// Helper functor, basically negates criterion for CopyIf
//struct IsUnique
//{
//  VTKM_EXEC_CONT bool operator()(vtkm::IdComponent isInOther) const { return isInOther == 0; }
//};

// Combine two TopologyGraphs
template <typename FieldType>
inline void TopologyGraph<FieldType>::MergeWith(TopologyGraph<FieldType>& other,
                                                  vtkm::cont::LogLevel timingsLogLevel,
                                                  std::string timingsMessage)
{ // Merge With
#ifdef DEBUG_PRINT
  this->DebugPrint("THIS TopologyGraph", __FILE__, __LINE__);
  other.DebugPrint("OTHER TopologyGraph", __FILE__, __LINE__);
#endif
  // Track timing of main steps
  vtkm::cont::Timer totalTimer; // Total time for each call
  totalTimer.Start();
  vtkm::cont::Timer timer; // Time individual steps
  timer.Start();
  std::stringstream timingsStream;

  // Create combined sort order
  // TODO This vector could potentially be implemented purely as a smart array handle to reduce memory usage
  IdArrayType overallSortOrder;
  overallSortOrder.Allocate(this->NumVertices + other.NumVertices);

  { // Create a new scope so that the following two vectors get deleted when leaving the scope
    auto thisIndices = vtkm::cont::ArrayHandleIndex(this->NumVertices); // A regular index array
    MarkOther markOtherFunctor;
    auto otherIndices = vtkm::cont::make_ArrayHandleTransform(
      vtkm::cont::ArrayHandleIndex(other.NumVertices), markOtherFunctor);
    contourtree_mesh_inc_ns::CombinedSimulatedSimplicityIndexComparator<FieldType>
      cssicFunctorExecObj(
        this->GlobalMeshIndex, other.GlobalMeshIndex, this->SortedValues, other.SortedValues);
    contourtree_mesh_inc_ns::CopyIntoCombinedArrayWorklet<true>
      copyIntoCombinedArrayWorkletLowerBound;
    this->Invoke(copyIntoCombinedArrayWorkletLowerBound,
                 thisIndices,
                 otherIndices,
                 cssicFunctorExecObj,
                 overallSortOrder);
    contourtree_mesh_inc_ns::CopyIntoCombinedArrayWorklet<false>
      copyIntoCombinedArrayWorkletUpperBound;
    this->Invoke(copyIntoCombinedArrayWorkletUpperBound,
                 otherIndices,
                 thisIndices,
                 cssicFunctorExecObj,
                 overallSortOrder);
  }
  timingsStream << "    " << std::setw(38) << std::left << "Create OverallSortOrder"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

#ifdef DEBUG_PRINT
  std::cout << "OverallSortOrder.size  " << overallSortOrder.GetNumberOfValues() << std::endl;
  PrintIndices("overallSortOrder", overallSortOrder);
  std::cout << std::endl;
#endif

  IdArrayType overallSortIndex;
  overallSortIndex.Allocate(overallSortOrder.GetNumberOfValues());
  {
    // Array decorator with functor returning 0, 1 for each element depending
    // on whethern the current value is different from the next.
    auto differentFromNextArr = vtkm::cont::make_ArrayHandleDecorator(
      overallSortOrder.GetNumberOfValues() - 1,
      mesh_dem_contourtree_mesh_inc::CombinedVectorDifferentFromNextDecoratorImpl{},
      overallSortOrder,
      this->GlobalMeshIndex,
      other.GlobalMeshIndex);

    // Compute the extended scan of our transformed combined vector
    vtkm::cont::Algorithm::ScanExtended(differentFromNextArr, overallSortIndex);
  }
  vtkm::Id numVerticesCombined =
    ArrayGetValue(overallSortIndex.GetNumberOfValues() - 1, overallSortIndex) + 1;

#ifdef DEBUG_PRINT
  std::cout << "OverallSortIndex.size  " << overallSortIndex.GetNumberOfValues() << std::endl;
  PrintIndices("overallSortIndex", overallSortIndex);
  std::cout << "numVerticesCombined: " << numVerticesCombined << std::endl;
  std::cout << std::endl;
#endif
  timingsStream << "    " << std::setw(38) << std::left << "Create OverallSortIndex"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

  // thisToCombinedSortOrder and otherToCombinedSortOrder
  IdArrayType thisToCombinedSortOrder;
  thisToCombinedSortOrder.Allocate(this->NumVertices);
  IdArrayType otherToCombinedSortOrder;
  otherToCombinedSortOrder.Allocate(other.NumVertices);
  contourtree_mesh_inc_ns::InitToCombinedSortOrderArraysWorklet
    initToCombinedSortOrderArraysWorklet;
  this->Invoke(initToCombinedSortOrderArraysWorklet,
               overallSortIndex,
               overallSortOrder,
               thisToCombinedSortOrder,
               otherToCombinedSortOrder);

#ifdef DEBUG_PRINT
  PrintIndices("thisToCombinedSortOrder", thisToCombinedSortOrder);
  PrintIndices("otherToCombinedSortOrder", otherToCombinedSortOrder);
#endif
  timingsStream << "    " << std::setw(38) << std::left << "Create This/OtherCombinedSortOrder"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

  // Map neighbor IDs to global ID (ID in the combined) and group them
  auto neighborConnectivityGlobalThis =
    make_ArrayHandleDecorator(this->NeighborConnectivity.GetNumberOfValues(),
                              mesh_dem_contourtree_mesh_inc::ApplyLookupTableDecoratorImpl{},
                              this->NeighborConnectivity,
                              thisToCombinedSortOrder);
  auto neighborConnectivityGlobalGroupsThis = vtkm::cont::make_ArrayHandleGroupVecVariable(
    neighborConnectivityGlobalThis, this->NeighborOffsets);

  auto neighborConnectivityGlobalOther =
    make_ArrayHandleDecorator(other.NeighborConnectivity.GetNumberOfValues(),
                              mesh_dem_contourtree_mesh_inc::ApplyLookupTableDecoratorImpl{},
                              other.NeighborConnectivity,
                              otherToCombinedSortOrder);
  auto neighborConnectivityGlobalGroupsOther = vtkm::cont::make_ArrayHandleGroupVecVariable(
    neighborConnectivityGlobalOther, other.NeighborOffsets);

  // Merge the two neighborhood connecitivy lists. First, we split neighbor connectivity
  // into three groups (i) vertices only in this, (ii) vertices only in other, (iii)
  // vertices in both meshes. We then compute cobmined neighbor connectivity for vertices
  // in both meshes. Finally, we copy them into the combined array.

  // Split vertices into groups (i) uniuqe this, (ii) unique other, (iii) in both
  // ... compute arrays containing 1 if the element is in the other respective array
  vtkm::cont::ArrayHandle<vtkm::IdComponent> thisToCombinedSortOrderIsDuplicate;
  thisToCombinedSortOrderIsDuplicate.Allocate(thisToCombinedSortOrder.GetNumberOfValues());
  vtkm::cont::ArrayHandle<vtkm::IdComponent> otherToCombinedSortOrderIsDuplicate;
  otherToCombinedSortOrderIsDuplicate.AllocateAndFill(otherToCombinedSortOrder.GetNumberOfValues(),
                                                      vtkm::IdComponent{ 0 });
  this->Invoke(contourtree_mesh_inc_ns::FindDuplicateInOtherWorklet{},
               thisToCombinedSortOrder,
               otherToCombinedSortOrder,
               thisToCombinedSortOrderIsDuplicate,
               otherToCombinedSortOrderIsDuplicate);

#ifdef DEBUG_PRINT
  PrintIndices("thisToCombinedSortOrderIsDuplicate", thisToCombinedSortOrderIsDuplicate);
  PrintIndices("otherToCombinedSortOrderIsDuplicate", otherToCombinedSortOrderIsDuplicate);
#endif
  // ... create lists for all groups to be used to restrict operations to them
  vtkm::cont::ArrayHandleIndex indicesThis(thisToCombinedSortOrder.GetNumberOfValues());
  vtkm::cont::ArrayHandleIndex indicesOther(otherToCombinedSortOrder.GetNumberOfValues());

  IdArrayType indicesThisUnique, indicesThisDuplicate;
  vtkm::cont::Algorithm::CopyIf(
    indicesThis, thisToCombinedSortOrderIsDuplicate, indicesThisUnique, IsUnique{});
  vtkm::cont::Algorithm::CopyIf(
    indicesThis, thisToCombinedSortOrderIsDuplicate, indicesThisDuplicate);

#ifdef DEBUG_PRINT
  PrintIndices("indicesThisUnique", indicesThisUnique);
  PrintIndices("indicesThisDuplicate", indicesThisDuplicate);
#endif

  IdArrayType indicesOtherUnique, indicesOtherDuplicate;
  vtkm::cont::Algorithm::CopyIf(
    indicesOther, otherToCombinedSortOrderIsDuplicate, indicesOtherUnique, IsUnique{});
  vtkm::cont::Algorithm::CopyIf(
    indicesOther, otherToCombinedSortOrderIsDuplicate, indicesOtherDuplicate);

#ifdef DEBUG_PRINT
  PrintIndices("indicesOtherUnique", indicesOtherUnique);
  PrintIndices("indicesOtherDuplicate", indicesOtherDuplicate);
#endif

  VTKM_ASSERT(indicesThisDuplicate.GetNumberOfValues() ==
              indicesOtherDuplicate.GetNumberOfValues());

  // Merge the neighbor groups for vertices that occur in both meshes
  // ... compute combined counts (with duplicates)
  auto neighborCountsThis = make_ArrayHandleOffsetsToNumComponents(this->NeighborOffsets);
  auto permutedNeighborCountsThis =
    vtkm::cont::make_ArrayHandlePermutation(indicesThisDuplicate, neighborCountsThis);
  auto neighborCountsOther = make_ArrayHandleOffsetsToNumComponents(other.NeighborOffsets);
  auto permutedNeighborCountsOther =
    vtkm::cont::make_ArrayHandlePermutation(indicesOtherDuplicate, neighborCountsOther);
  vtkm::cont::ArrayHandle<vtkm::IdComponent> combinedCommonNeighborCountSums;
  vtkm::cont::Algorithm::Transform(permutedNeighborCountsThis,
                                   permutedNeighborCountsOther,
                                   combinedCommonNeighborCountSums,
                                   vtkm::Sum());

  // ... merge sorted lists
  // ...... create output arrays/groups
  vtkm::Id unpackedCombinedCommonNeighborConnectivitySize;
  IdArrayType unpackedCombinedCommonNeighborOffsets;
  vtkm::cont::ConvertNumComponentsToOffsets(combinedCommonNeighborCountSums,
                                            unpackedCombinedCommonNeighborOffsets,
                                            unpackedCombinedCommonNeighborConnectivitySize);
  IdArrayType unpackedCombinedCommonNeighborConnectivity;
  unpackedCombinedCommonNeighborConnectivity.Allocate(
    unpackedCombinedCommonNeighborConnectivitySize);
  auto unpackedCombinedCommonNeighborConnectivityGroups = make_ArrayHandleGroupVecVariable(
    unpackedCombinedCommonNeighborConnectivity, unpackedCombinedCommonNeighborOffsets);

  // ....... create permuted input arrays/groups
  auto permutedNeighborConnectivityGlobalGroupsThis =
    make_ArrayHandlePermutation(indicesThisDuplicate, neighborConnectivityGlobalGroupsThis);
  auto permutedNeighborConnectivityGlobalGroupsOther =
    make_ArrayHandlePermutation(indicesOtherDuplicate, neighborConnectivityGlobalGroupsOther);

  // ........ create array for actual counts of unique neighbors
  vtkm::cont::ArrayHandle<vtkm::IdComponent> packedCombinedCommonNeighborCounts;
  packedCombinedCommonNeighborCounts.Allocate(combinedCommonNeighborCountSums.GetNumberOfValues());

  // ........ perform merge
  this->Invoke(contourtree_mesh_inc_ns::MergeSortedListsWithoutDuplicatesWorklet{},
               permutedNeighborConnectivityGlobalGroupsThis,
               permutedNeighborConnectivityGlobalGroupsOther,
               unpackedCombinedCommonNeighborConnectivityGroups,
               packedCombinedCommonNeighborCounts);

  // ... pack sorted lists
  // ...... create the new offsets array for the merged groups (without duplicates).
  vtkm::Id packedCombinedCommonNeighborConnectivitySize;
  vtkm::cont::ArrayHandle<vtkm::Id> packedCombinedCommonNeighborOffsets;
  vtkm::cont::ConvertNumComponentsToOffsets(packedCombinedCommonNeighborCounts,
                                            packedCombinedCommonNeighborOffsets,
                                            packedCombinedCommonNeighborConnectivitySize);

  // ...... create a new grouped array for the packed connectivity
  vtkm::cont::ArrayHandle<vtkm::Id> packedCombinedCommonNeighborConnectivity;
  packedCombinedCommonNeighborConnectivity.Allocate(packedCombinedCommonNeighborConnectivitySize);
  auto packedCommonNeighborConnectivityGroups = vtkm::cont::make_ArrayHandleGroupVecVariable(
    packedCombinedCommonNeighborConnectivity, packedCombinedCommonNeighborOffsets);

  // ...... copy data to the packed array.
  this->Invoke(contourtree_mesh_inc_ns::CopyNeighborsToPackedArray{},
               unpackedCombinedCommonNeighborConnectivityGroups,
               packedCommonNeighborConnectivityGroups);

  // Create array for all three groups
  // ... create combined counts array
  IdArrayType combinedNeighborCounts;
  combinedNeighborCounts.Allocate(numVerticesCombined);

  auto thisOnlyToCombinedSortOrder =
    make_ArrayHandlePermutation(indicesThisUnique, thisToCombinedSortOrder);
  auto otherOnlyToCombinedSortOrder =
    make_ArrayHandlePermutation(indicesOtherUnique, otherToCombinedSortOrder);
  auto commonToCombinedSortOrder =
    make_ArrayHandlePermutation(indicesThisDuplicate, thisToCombinedSortOrder);

  CopyArrayByIndices(
    neighborCountsThis, indicesThisUnique, combinedNeighborCounts, thisOnlyToCombinedSortOrder);
  CopyArrayByIndices(
    neighborCountsOther, indicesOtherUnique, combinedNeighborCounts, otherOnlyToCombinedSortOrder);
  auto commonCombinedNeighborCounts =
    make_ArrayHandlePermutation(commonToCombinedSortOrder, combinedNeighborCounts);
  vtkm::cont::Algorithm::Copy(packedCombinedCommonNeighborCounts, commonCombinedNeighborCounts);

  // ... create offsets and allocate combinedNeighborConnectivity array
  vtkm::Id combinedNeighborConnectivitySize;
  vtkm::cont::ArrayHandle<vtkm::Id> combinedNeighborOffsets;
  vtkm::cont::ConvertNumComponentsToOffsets(
    combinedNeighborCounts, combinedNeighborOffsets, combinedNeighborConnectivitySize);
  IdArrayType combinedNeighborConnectivity;
  combinedNeighborConnectivity.Allocate(combinedNeighborConnectivitySize);
  auto combinedNeighborConnectivityGroups =
    make_ArrayHandleGroupVecVariable(combinedNeighborConnectivity, combinedNeighborOffsets);

  // ... copy the connectivity data including  previously merged lists
  CopyVecArrayByIndices(neighborConnectivityGlobalGroupsThis,
                        indicesThisUnique,
                        combinedNeighborConnectivityGroups,
                        thisOnlyToCombinedSortOrder);
  CopyVecArrayByIndices(neighborConnectivityGlobalGroupsOther,
                        indicesOtherUnique,
                        combinedNeighborConnectivityGroups,
                        otherOnlyToCombinedSortOrder);
  auto commonCombinedNeighborConnectivityGroups =
    make_ArrayHandlePermutation(commonToCombinedSortOrder, combinedNeighborConnectivityGroups);
  this->Invoke(contourtree_mesh_inc_ns::CopyIntoCombinedNeighborsWorklet{},
               packedCommonNeighborConnectivityGroups,
               commonCombinedNeighborConnectivityGroups);
  // Why doesn't the following copy work instead?
  // vtkm::cont::Algorithm::Copy(packedCommonNeighborConnectivityGroups, commonCombinedNeighborConnectivityGroups);

  timingsStream << "    " << std::setw(38) << std::left << "Compute CombinedNeighborConnectivity"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

  // Compute combined global mesh index arrays
  IdArrayType combinedGlobalMeshIndex;
  combinedGlobalMeshIndex.Allocate(numVerticesCombined);
  { // make sure arrays used for copy go out of scope
    auto permutedCombinedGlobalMeshIndex =
      vtkm::cont::make_ArrayHandlePermutation(thisToCombinedSortOrder, combinedGlobalMeshIndex);
    vtkm::cont::Algorithm::Copy(GlobalMeshIndex, permutedCombinedGlobalMeshIndex);
  }
  { // make sure arrays used for copy go out of scope
    auto permutedCombinedGlobalMeshIndex =
      vtkm::cont::make_ArrayHandlePermutation(otherToCombinedSortOrder, combinedGlobalMeshIndex);
    vtkm::cont::Algorithm::Copy(other.GlobalMeshIndex, permutedCombinedGlobalMeshIndex);
  }

  timingsStream << "    " << std::setw(38) << std::left << "Create CombinedGlobalMeshIndex"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

  // Compute combined sorted values
  vtkm::cont::ArrayHandle<FieldType> combinedSortedValues;
  combinedSortedValues.Allocate(numVerticesCombined);
  { // make sure arrays used for copy go out of scope
    auto permutedCombinedSortedValues =
      vtkm::cont::make_ArrayHandlePermutation(thisToCombinedSortOrder, combinedSortedValues);
    vtkm::cont::Algorithm::Copy(SortedValues, permutedCombinedSortedValues);
  }
  { // make sure arrays used for copy go out of scope
    auto permutedCombinedSortedValues =
      vtkm::cont::make_ArrayHandlePermutation(otherToCombinedSortOrder, combinedSortedValues);
    vtkm::cont::Algorithm::Copy(other.SortedValues, permutedCombinedSortedValues);
  }

  timingsStream << "    " << std::setw(38) << std::left << "Create CombinedSortedValues"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

  // Swap in combined version. VTKM ArrayHandles are smart so we can just swap in the new for the old
  this->SortedValues = combinedSortedValues;
  this->GlobalMeshIndex = combinedGlobalMeshIndex;
  this->NeighborConnectivity = combinedNeighborConnectivity;
  this->NeighborOffsets = combinedNeighborOffsets;
  this->NumVertices = SortedValues.GetNumberOfValues();
  this->SortIndices = vtkm::cont::ArrayHandleIndex(this->NumVertices);
  this->SortOrder = vtkm::cont::ArrayHandleIndex(this->NumVertices);

  timingsStream << "    " << std::setw(38) << std::left << "Swap in new arrays"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();

  // Re-compute maximum number of neigbours
  ComputeMaxNeighbors();

  timingsStream << "    " << std::setw(38) << std::left << "Compute MaxNeighbors"
                << ": " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();
  timingsStream << "    " << std::setw(38) << std::left << "Total time MergeWith"
                << ": " << totalTimer.GetElapsedTime() << " seconds" << std::endl;
  timer.Start();
  // Record the times we logged
  VTKM_LOG_S(timingsLogLevel,
             std::endl
               << "    ---------------- ContourTreeMesh MergeWith ---------------------"
               << std::endl
               << timingsMessage << timingsStream.str());
  // Prevent unused parameter warning when compiled without logging
  (void)timingsLogLevel;
  (void)timingsMessage;

#ifdef DEBUG_PRINT
  // Print the contents fo this for debugging
  DebugPrint("ContourTreeMeshes merged", __FILE__, __LINE__);
#endif
} // Merge With


template <typename FieldType>
inline void TopologyGraph<FieldType>::Save(const char* filename) const
{
  std::ofstream os(filename);
  SaveVector(os, this->SortedValues);
  SaveVector(os, this->GlobalMeshIndex);
  SaveVector(os, this->NeighborConnectivity);
  SaveVector(os, this->NeighborOffsets);
}

template <typename FieldType>
inline void TopologyGraph<FieldType>::Load(const char* filename)
{
  std::cout << "{ContourTreeMesh.h : Load File\n";
  std::cout << " ... filename = " << filename << "}\n";

  std::ifstream is(filename);
  if (!is.is_open())
  {
    throw vtkm::io::ErrorIO(std::string("Unable to open file: ") + std::string(filename));
  }
  LoadVector(is, this->SortedValues);
  LoadVector(is, this->GlobalMeshIndex);
  LoadVector(is, this->NeighborConnectivity);
  LoadVector(is, this->NeighborOffsets);
  this->ComputeMaxNeighbors();
  this->NumVertices = this->SortedValues.GetNumberOfValues();
  this->SortOrder = vtkm::cont::ArrayHandleIndex(this->NumVertices);
  this->SortIndices = vtkm::cont::ArrayHandleIndex(this->NumVertices);
}

template <typename FieldType>
template <typename ValueType>
inline void TopologyGraph<FieldType>::SaveVector(
  std::ostream& os,
  const vtkm::cont::ArrayHandle<ValueType>& vec) const
{
  vtkm::Id numVals = vec.GetNumberOfValues();
  //os.write(rXeinterpret_cast<const char*>(&numVals), sizeof(ValueType));
  os << numVals << ": ";
  auto vecPortal = vec.ReadPortal();
  for (vtkm::Id i = 0; i < numVals; ++i)
    os << vecPortal.Get(i) << " ";
  //os.write(reinterpret_cast<const char*>(vecPortal.Get(i)), sizeof(ValueType));
  os << std::endl;
}

template <typename FieldType>
template <typename ValueType>
inline void TopologyGraph<FieldType>::LoadVector(std::istream& is,
                                                   vtkm::cont::ArrayHandle<ValueType>& vec)
{
  vtkm::Id numVals;
  is >> numVals;
  char colon = is.get();
  if (colon != ':')
  {
    throw vtkm::io::ErrorIO("Error parsing file");
  }

  vec.Allocate(numVals);
  auto vecPortal = vec.WritePortal();
  ValueType val;
  for (vtkm::Id i = 0; i < numVals; ++i)
  {
    is >> val;
    vecPortal.Set(i, val);
  }
}

template <typename FieldType>
//inline MeshBoundaryTopologyGraphExec TopologyGraph<FieldType>::GetMeshBoundaryExecutionObject(
//        MeshBoundaryContourTreeMesh
inline MeshBoundaryContourTreeMeshExec TopologyGraph<FieldType>::GetMeshBoundaryExecutionObject(
  vtkm::Id3 globalSize,
  vtkm::Id3 minIdx,
  vtkm::Id3 maxIdx) const
{
  return MeshBoundaryContourTreeMeshExec(this->GlobalMeshIndex, globalSize, minIdx, maxIdx);
}

template <typename FieldType>
inline MeshBoundaryContourTreeMeshExec TopologyGraph<FieldType>::GetMeshBoundaryExecutionObject() const
{
  return MeshBoundaryContourTreeMeshExec(this->GlobalMeshIndex,
                                         vtkm::Id3{ this->NumVertices, this->NumVertices, this->NumVertices },
                                         vtkm::Id3{ 0, 0, 0 },
                                         vtkm::Id3{ this->NumVertices, this->NumVertices, this->NumVertices });
}


template <typename FieldType>
inline void TopologyGraph<FieldType>::GetBoundaryVertices(
  IdArrayType& boundaryVertexArray,                    // output
  IdArrayType& boundarySortIndexArray,                 // output
  MeshBoundaryContourTreeMeshExec* meshBoundaryExecObj //input
) const
{
  // start by generating a temporary array of indices
  auto indexArray = vtkm::cont::ArrayHandleIndex(this->GlobalMeshIndex.GetNumberOfValues());
  // compute the boolean array indicating which values lie on the boundary
  vtkm::cont::ArrayHandle<bool> isOnBoundary;

//  ComputeMeshBoundaryTopologyGraph computeMeshBoundaryTopologyGraphWorklet;
  ComputeMeshBoundaryContourTreeMesh computeMeshBoundaryTopologyGraphWorklet;
  this->Invoke(computeMeshBoundaryTopologyGraphWorklet,
               indexArray,           // input
               *meshBoundaryExecObj, // input
               isOnBoundary          // outut
  );

  // we will conditionally copy the boundary vertices' indices, capturing the end iterator to compute the # of boundary vertices
  vtkm::cont::Algorithm::CopyIf(indexArray, isOnBoundary, boundaryVertexArray);
  // duplicate these into the index array, since the BRACT uses indices into the underlying mesh anyway
  vtkm::cont::Algorithm::Copy(boundaryVertexArray, boundarySortIndexArray);
}

} // namespace contourtree_augmented
} // worklet
} // vtkm

//#undef DEBUG_PRINT

#endif
