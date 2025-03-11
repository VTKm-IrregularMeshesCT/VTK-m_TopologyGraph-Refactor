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

#ifndef vtk_m_worklet_contourtree_augmented_activegraph_h
#define vtk_m_worklet_contourtree_augmented_activegraph_h

#include <iomanip>
#include <numeric>

// local includes
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ArrayTransforms.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/MergeTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/MeshExtrema.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/PrintVectors.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/BuildChainsWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/BuildTrunkWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/CompactActiveEdgesComputeNewVertexOutdegree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/CompactActiveEdgesTransferActiveEdges.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/EdgePeakComparator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/FindGoverningSaddlesWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/FindSuperAndHyperNodesWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/HyperArcSuperNodeComparator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/InitializeActiveEdges.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/InitializeActiveGraphVertices.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/InitializeEdgeFarFromActiveIndices.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/InitializeHyperarcsFromActiveIndices.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/InitializeNeighbourhoodMasksAndOutDegrees.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SetArcsConnectNodes.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SetArcsSetSuperAndHypernodeArcs.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SetArcsSlideVertices.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SetHyperArcsWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SetSuperArcsSetTreeHyperparents.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SetSuperArcsSetTreeSuperarcs.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/SuperArcNodeComparator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/TransferRegularPointsWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/TransferSaddleStartsResetEdgeFar.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/TransferSaddleStartsSetNewOutdegreeForSaddles.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/activegraph/TransferSaddleStartsUpdateEdgeSorter.h>


//VTKM includes
#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayGetValues.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayHandlePermutation.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/Error.h>
#include <vtkm/cont/Invoker.h>

//#include <vtkm/filter/scalar_topology/worklet/contourtree/PrintVectors.h>

#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/ContourTreeMesh.h>
// Adding the new TopologyGraph Class
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/TopologyGraph.h>

#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/MeshBoundaryContourTreeMesh.h>

#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/DataSetMeshTriangulation3DFreudenthal.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/DataSetMeshTriangulation2DFreudenthal.h>

//#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/mesh_boundary/MeshBoundaryContourTreeMesh.h>

namespace active_graph_inc_ns = vtkm::worklet::contourtree_augmented::active_graph_inc;

#define DEBUG_PRINT_PACTBD 0


namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{


class ActiveGraph
{ // class ActiveGraph

  template <typename T, typename S>
  static T GetLastValue(const vtkm::cont::ArrayHandle<T, S>& ah)
  {
    return vtkm::cont::ArrayGetValue(ah.GetNumberOfValues() - 1, ah);
  }

public:
  vtkm::cont::Invoker Invoke;

  // we also need the orientation of the edges (i.e. is it join or split)
  bool IsJoinGraph;

  // we will store the number of iterations the computation took here
  vtkm::Id NumIterations;

  // ARRAYS FOR NODES IN THE TOPOLOGY GRAPH

  // for each vertex, we need to know where it is in global sort order / mesh
  IdArrayType GlobalIndex;

  // the hyperarcs - i.e. the pseudoextremum defining the hyperarc the vertex is on
  IdArrayType Hyperarcs;

  // the first edge for each vertex
  IdArrayType FirstEdge;

  // the outdegree for each vertex
  IdArrayType Outdegree;

  // ARRAYS FOR EDGES IN THE TOPOLOGY GRAPH

  // we will also need to keep track of both near and far ends of each edge
  IdArrayType EdgeFar;
  IdArrayType EdgeNear;

  // these now track the active nodes, edges, &c.:
  IdArrayType ActiveVertices;
  IdArrayType ActiveEdges;

  // and an array for sorting edges
  IdArrayType EdgeSorter;

  // temporary arrays for super/hyper ID numbers
  IdArrayType SuperID;
  IdArrayType HyperID;

  // variables tracking size of super/hyper tree
  vtkm::Id NumSupernodes;
  vtkm::Id NumHypernodes;

  // BASIC ROUTINES: CONSTRUCTOR, PRINT, &c.

  // constructor takes necessary references
  ActiveGraph(bool IsJoinGraph);

  // initialises the active graph
  template <class Mesh>
  void Initialise(Mesh& mesh, const MeshExtrema& meshExtrema);

// won't work like this ...
  // ... wishful thinking
//  template <>
//  void Initialise<ContourTreeMesh>(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema);
  // instead do overload for ContourTreeMesh (later TopologyGraph):
  // UNCOMMENT
  void Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<int>& mesh,  const MeshExtrema& meshExtrema);

//                vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation2DFreudenthal
//  void Initialise(vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation3DFreudenthal& mesh, const MeshExtrema& meshExtrema);
////     Initialise(vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation2DFreudenthal&, vtkm::worklet::contourtree_augmented::MeshExtrema const&)
//  void Initialise(vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation2DFreudenthal& mesh, const MeshExtrema& meshExtrema);

//  void Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<double>& mesh,  const MeshExtrema& meshExtrema);
//  void Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<float>& mesh,  const MeshExtrema& meshExtrema);
////  void vtkm::worklet::contourtree_augmented::ActiveGraph::Initialise<vtkm::worklet::contourtree_augmented::ContourTreeMesh<double> >(vtkm::worklet::contourtree_augmented::ContourTreeMesh<double>&, vtkm::worklet::contourtree_augmented::MeshExtrema const&)


  // routine that computes the merge tree from the active graph
  // was previously Compute()
  void MakeMergeTree(MergeTree& tree, MeshExtrema& meshExtrema);

  // sorts saddle starts to find governing saddles
  void FindGoverningSaddles();

  // marks now regular points for removal
  void TransferRegularPoints();

  // compacts the active vertex list
  void CompactActiveVertices();

  // compacts the active edge list
  void CompactActiveEdges();

  // builds the chains for the new active vertices
  void BuildChains();

  // suppresses non-saddles for the governing saddles pass
  void TransferSaddleStarts();

  // sets all remaining active vertices
  void BuildTrunk();

  // finds all super and hyper nodes, numbers them & sets up arrays for lookup
  void FindSuperAndHyperNodes(MergeTree& tree);

  // uses active graph to set superarcs & hyperparents in merge tree
  void SetSuperArcs(MergeTree& tree);

  // uses active graph to set hypernodes in merge tree
  void SetHyperArcs(MergeTree& tree);

  // uses active graph to set arcs in merge tree
  void SetArcs(MergeTree& tree, MeshExtrema& meshExtrema);

  // Allocate the vertex array
  void AllocateVertexArrays(vtkm::Id nElems);

  // Allocate the edge array
  void AllocateEdgeArrays(vtkm::Id nElems);

  // releases temporary arrays
  void ReleaseTemporaryArrays();

  // prints the contents of the active graph in a standard format
  void DebugPrint(const char* message, const char* fileName, long lineNum);

}; // class ActiveGraph


// constructor takes necessary references
inline ActiveGraph::ActiveGraph(bool isJoinGraph)
  : Invoke()
  , IsJoinGraph(isJoinGraph)
{ // constructor
  this->NumIterations = 0;
  this->NumSupernodes = 0;
  this->NumHypernodes = 0;
} // constructor



// TODO 2024-03-12
// ... only keep just bare minimum in what's needed in the ActiveGraph Initialise














































// UNCOMMMENT for fixed 64-nbor fix

// started refactoring 2024-03-15
// ... first by removing everything that we don't need
// ... just doing manual active graph computation at this point!

// 2024-03-21 refactoring, removing messy debug code


// C++ Partial Template Specialization

// initialises the active graph
//template </*using C++ Partial Template Specialization, leave empty */ >
//template <>
//inline void ActiveGraph::Initialise<ContourTreeMesh>(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema)
//inline void ActiveGraph::Initialise(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema)

//inline void ActiveGraph::Initialise(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema)
inline void ActiveGraph::Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<int>& mesh,
                                    const MeshExtrema& meshExtrema)
{ // InitialiseActiveGraph()
  // reference to the correct array in the extrema
  /// DEBUG PRINT std::cout << "--------------- ActiveGraph::Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh --------------- \n";
  /// DEBUG PRINT std::cout << "C++ Partial Template Specialization\n";
  #if DEBUG_PRINT_PACTBD
      PrintHeader(meshExtrema.Peaks.GetNumberOfValues());
      PrintIndices("MeshExtremaPeaks", meshExtrema.Peaks);
      PrintIndices("MeshExtremaPits", meshExtrema.Pits);
  #endif

  /// DEBUG PRINT std::cout << "A1\n";
  // if computing the Join Tree, we need Peaks, ...
  // ... and Pits for Merge Tree
  // ... extrema is the array which holds which peak/pit a vertex leads up/down to after pointer doubling
  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

  // Boolean flag to determine whether the individual vertices outDegrees count ...
  // ... the number of edge end-points higher (up) or lower (down) than itself
  bool up_or_down = this->IsJoinGraph ? true : false;

  #if DEBUG_PRINT_PACTBD
      for(unsigned i = 0; i < mesh.NumVertices; i++)
      {
          PrintIndices("mesh.NeighborConnectivity", mesh.NeighborConnectivity);
          PrintIndices("mesh.NeighborOffsets", mesh.NeighborOffsets);
          PrintIndices("MeshExtremaPeaksAFTER", extrema);
      }
      std::cout << "------------------------------- fine until here -------------------------------\n";
  #endif

  /// DEBUG PRINT std::cout << "A2\n";
  // No longer rely on 'neighbourhoodMasks' ...
  // Neighbourhood masks (one bit set per connected component in neighbourhood ...
  // ... as arbitrary graphs can have nbhoods greater than 64 ...
  // ... which is the bitcount of vtkm::Id's
  /// DEBUG PRINT std::cout << "A3\n";
  /// DEBUG PRINT std::cout << "A4\n";
  /// DEBUG PRINT std::cout << "A5\n";

  // Inverse index converts between sortID and meshID ...
  // ... here we just assume they're the same ...
  // ... which means each vertex ID is its rank number ...
  // ... 0-lowest N-highest
  IdArrayType inverseIndex;
  inverseIndex.Allocate(mesh.NumVertices);
  using IdPortalType = IdArrayType::WritePortalType;
  IdPortalType inversePortal = inverseIndex.WritePortal();
  for(vtkm::Id i = 0; i < mesh.NumVertices; i++)
  {
      inversePortal.Set(i, i);
  }

  /// DEBUG PRINT std::cout << "A6\n";
  // No optimisation - keep all mesh vertices in the active graph:
  vtkm::Id nCriticalPoints = mesh.NumVertices;
  /// DEBUG PRINT std::cout << "A7\n";
  // we need to keep track of what the index of each vertex is in the active graph
  // for most vertices, this should have the NO_SUCH_VERTEX flag set
  AllocateVertexArrays(
    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices
  /// DEBUG PRINT std::cout << "A8\n";
  // our processing now depends on the degree of the vertex
  // but basically, we want to set up the arrays for this vertex:
  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
  // GlobalIndex stores the index in the join tree for later access
  IdArrayType activeIndices;
  activeIndices.Allocate(mesh.NumVertices);
  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);
  /// DEBUG PRINT std::cout << "A9\n";


//  PrintIndices("FOR_INIT:sortIndexArray",    sortIndexArray);
//  PrintIndices("FOR_INIT:outDegrees",        outDegrees);
  #if DEBUG_PRINT_PACTBD
      PrintIndices("FOR_INIT:inverseIndex",      inverseIndex);
      PrintIndices("FOR_INIT:extrema",           extrema);
      PrintIndices("FOR_INIT:activeIndices",     activeIndices);
  #endif

  /// DEBUG PRINT std::cout << "A10\n";
  #if DEBUG_PRINT_PACTBD
      PrintIndices("INIT:activeIndices",     activeIndices);
      PrintIndices("INIT:GlobalIndex",    this->GlobalIndex);
      PrintIndices("INIT:Outdegree",      this->Outdegree);
      PrintIndices("INIT:Hyperarcs",      this->Hyperarcs);
      PrintIndices("INIT:ActiveVertices", this->ActiveVertices);
  #endif


  // ---------------------------------- Compute All Necessary Arrays ---------------------------------- //

  IdPortalType aixPortal = activeIndices.WritePortal();         // 1
  IdPortalType avPortal  = this->ActiveVertices.WritePortal();  // 2
  IdPortalType haPortal  = this->Hyperarcs.WritePortal();       // 3
  IdPortalType odPortal  = this->Outdegree.WritePortal();       // 4
  IdPortalType giPortal  = this->GlobalIndex.WritePortal();     // 5

  const auto exPortal  = extrema.ReadPortal();
  const auto outReadPortal = this->Outdegree.ReadPortal();

  const auto meshNbhoodConPortal = mesh.NeighborConnectivity.ReadPortal();
  const auto meshNbhoodOffPortal = mesh.NeighborOffsets.ReadPortal();

  std::vector<vtkm::Id> farEndsUP;
  std::vector<vtkm::Id> farEndsDOWN;

  // compute up/down outdegrees degrees
  for(vtkm::Id i = 0; i < nCriticalPoints; i++)
  {
      avPortal.Set(i, i);
      giPortal.Set(i, i);
      haPortal.Set(i, MaskedIndex(exPortal.Get(i)) );
      // For every vertex, work out whether it is critical
      // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
      // All vertices of outdegree 0 must be extrema
      // Saddle points must be at least outdegree 2, so this is a correct test
      // BUT it is possible to overestimate the degree of a non-extremum,
      // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph
      if(outReadPortal.Get(i) != 1)
      {
        aixPortal.Set(i, i);
      }
      else
      {
        aixPortal.Set(i, 0);
      }


      vtkm::Id start_interval = meshNbhoodOffPortal.Get(i);
      vtkm::Id end_interval   = meshNbhoodOffPortal.Get(i+1);

      int conditional_count_up   = 0;
      int conditional_count_down = 0;

      for(vtkm::Id j = start_interval; j < end_interval; j++)
      {
        if(meshNbhoodConPortal.Get(j) > i)
        {
            conditional_count_up++;
            farEndsUP.push_back(meshNbhoodConPortal.Get(j));
        }
        else
        {
            conditional_count_down++;
            farEndsDOWN.push_back(meshNbhoodConPortal.Get(j));
        }
      }

      if(up_or_down)
      {
          odPortal.Set(i, conditional_count_up);
      }
      else
      {
          odPortal.Set(i, conditional_count_down);
      }

  }

  // Compute up/down outdegrees degrees
  for(vtkm::Id i = 0; i < nCriticalPoints; i++)
  {
      if(outReadPortal.Get(i) != 1)
      {
        aixPortal.Set(i, i);
      }
      else
      {
        aixPortal.Set(i, 0);
      }
  }


  #if DEBUG_PRINT_PACTBD
      PrintIndices("VX:ActiveIndices",  activeIndices);
      PrintIndices("VX:GlobalIndex",    this->GlobalIndex);
      PrintIndices("VX:Outdegree",      this->Outdegree);
      PrintIndices("VX:Hyperarcs",      this->Hyperarcs);
      PrintIndices("VX:ActiveVertices", this->ActiveVertices);
  #endif
  // now we need to compute the FirstEdge array from the outDegrees
  this->FirstEdge.Allocate(nCriticalPoints);
  // STD Version of the prefix sum
  //this->FirstEdge.WritePortal().Set(0, 0);
  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
  // VTKM Version of the prefix sum
  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
  // Compute the number of critical edges

  vtkm::Id nCriticalEdges =
    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);
    // used to set nCriticalEdges = 56 manually, but the above works

  AllocateEdgeArrays(nCriticalEdges);
  /// DEBUG PRINT std::cout << "A11\n";
//  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
  active_graph_inc_ns::InitializeActiveEdges<vtkm::worklet::contourtree_augmented::ContourTreeMesh<int>> initActiveEdgesWorklet;

  #if DEBUG_PRINT_PACTBD
      std::cout << "vvv ALL ACTIVE EDGES INITIALIZED! vvv\n";
      PrintIndices("ACTIVE EDGES", this->ActiveEdges);
      PrintIndices("NEAR", this->EdgeNear);
      PrintIndices("FAR", this->EdgeFar);
      std::cout << "^^^ ALL ACTIVE EDGES INITIALIZED! ^^^\n";
      std::cout << "A11.1\n";
      // try doing this manually:
      std::cout << "initActiveEdgesWorkled params: Outdegree=\n";
    //  << this->Outdegree << ", GlobalIndex=" << this->GlobalIndex << "\n";
    //  PrintIndices("outDegrees?: ", outDegrees);
      PrintIndices("OUTDEGREE?: ", this->Outdegree);
      PrintIndices("GlobalIDX(sort)?: ", this->GlobalIndex);
      std::cout << "nCriticalEdges?: " << nCriticalEdges << "\n";
  #endif


  std::vector<vtkm::Id> indices;
  std::vector<vtkm::Id> low_arr;
  std::vector<vtkm::Id> up_arr;

  IdPortalType aePortal = this->ActiveEdges.WritePortal();
  IdPortalType enPortal = this->EdgeNear.WritePortal();
  IdPortalType efPortal = this->EdgeFar.WritePortal();

  auto readOutDegs = this->Outdegree.ReadPortal();

  unsigned total = 0;
  if(up_or_down)
  {
      this->ActiveEdges.Allocate(farEndsUP.size());
      this->EdgeNear.Allocate(farEndsUP.size());
      this->EdgeFar.Allocate(farEndsUP.size());

      for(unsigned i = 0; i < farEndsUP.size(); i++)
      {
        aePortal.Set(i, i);
        efPortal.Set(i, MaskedIndex(exPortal.Get(farEndsUP[i])) ); // replaced 2024-03-20
      }

      total = 0;
      for(unsigned i = 0; i < mesh.NumVertices; i++)
      {
          for (unsigned j = 0; j < readOutDegs.Get(i); j++)
          {
              enPortal.Set(total, i);
              total++;
          }
      }
  }
  else
  {
      this->ActiveEdges.Allocate(farEndsDOWN.size());
      this->EdgeNear.Allocate(farEndsDOWN.size());
      this->EdgeFar.Allocate(farEndsDOWN.size());

      for(unsigned i = 0; i < farEndsDOWN.size(); i++)
      {
          aePortal.Set(i, i);
          efPortal.Set(i, MaskedIndex(exPortal.Get(farEndsDOWN[i])) ); // replaced 2024-03-20
      }

      total = 0;
      for(unsigned i = 0; i < mesh.NumVertices; i++)
      {
          for (unsigned j = 0; j < readOutDegs.Get(i); j++)
          {
              enPortal.Set(total, i);
              total++;
          }
      }
  }

  /// DEBUG PRINT std::cout << "A11.1 FINISHED!\n";
  if(up_or_down)
  {
  /// DEBUG PRINT   std::cout << "UP\n";
  }
  else
  {
  /// DEBUG PRINT std::cout << "DOWN\n";
  }

  #if DEBUG_PRINT_PACTBD
        printf("INVOKER() on: %s\n", __PRETTY_FUNCTION__);
        std::cout << "Invoker3\n";
      std::cout << "vvv ALL ACTIVE EDGES GATHERED! vvv\n";
      PrintIndices("ACTIVE EDGES", this->ActiveEdges);          // 6
      PrintIndices("NEAR", this->EdgeNear);                     // 7
      PrintIndices("FAR", this->EdgeFar);                       // 8
      std::cout << "^^^ ALL ACTIVE EDGES GATHERED! ^^^\n";
  #endif

  /// DEBUG PRINT std::cout << "A11.2\n";
  // now we have to go through and set the far ends of the new edges using the
  // inverse index array
  /// DEBUG PRINT std::cout << "A11.3\n";
  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
  /// DEBUG PRINT std::cout << "A11.4\n";
  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

  DebugPrint("Active Graph Started", __FILE__, __LINE__);
  /// DEBUG PRINT std::cout << "A12\n";
  // then we loop through the active vertices to convert their indices to active graph indices
  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);
  /// DEBUG PRINT std::cout << "A13\n";
  // finally, allocate and initialise the edgeSorter array
  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
  /// DEBUG PRINT std::cout << "A14\n";
  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);

  #if DEBUG_PRINT_PACTBD
        printf("INVOKER() on: %s\n", __PRETTY_FUNCTION__);
        std::cout << "Invoker3\n";
  std::cout << "ActiveGraph - finished\n";
  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);

  std::cout << "C++ Partial Template Specialization\n";
  #endif

} // InitialiseActiveGraph()

// UNCOMMENT AUTOMATED

































//inline void Initialise(vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation3DFreudenthal& mesh, const MeshExtrema& meshExtrema)
//{

//}

//inline void Initialise(vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation2DFreudenthal& mesh, const MeshExtrema& meshExtrema)
//{

//}

//inline void Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<double>& mesh,  const MeshExtrema& meshExtrema)
//{

//}

//inline void Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<float>& mesh,  const MeshExtrema& meshExtrema)
//{

//}












//// initialises the active graph - FROM THE OFFICIAL VTK-M VERSION!
//// ... (WITH A MODIFICATION TO PRINT ALL THE EDGES OF THE MESH!!!!!!!!!!!!!!!!!
//template <class Mesh>
//inline void ActiveGraph::Initialise(Mesh& mesh, const MeshExtrema& meshExtrema)
////inline void ActiveGraph::Initialise(vtkm::worklet::contourtree_augmented::DataSetMeshTriangulation3DFreudenthal& mesh, const MeshExtrema& meshExtrema)
//{ // InitialiseActiveGraph()
//  // reference to the correct array in the extrema
//  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

//  // For every vertex, work out whether it is critical
//  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
//  // All vertices of outdegree 0 must be extrema
//  // Saddle points must be at least outdegree 2, so this is a correct test
//  // BUT it is possible to overestimate the degree of a non-extremum,
//  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph

//  // Neighbourhood mask (one bit set per connected component in neighbourhood
//  IdArrayType neighbourhoodMasks;
//  neighbourhoodMasks.Allocate(mesh.NumVertices);
//  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
//  outDegrees.Allocate(mesh.NumVertices);

//  // Initialize the nerighborhoodMasks and outDegrees arrays
//  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
//  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
//  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
//    this->IsJoinGraph);

//  this->Invoke(initNeighMasksAndOutDegWorklet,
//               sortIndexArray,
//               mesh,
//               neighbourhoodMasks, // output
//               outDegrees);        // output

//  std::cout << "Try to get all Freudenthal3D edges ...\n";
//  using IdPortalType = IdArrayType::ReadPortalType;
//  IdPortalType outDegPortal = outDegrees.ReadPortal();
//  IdPortalType nbhdmaPortal = neighbourhoodMasks.ReadPortal();

//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
//      std::cout << i << " " << outDegPortal.Get(i) << " - ";// << "\n";

//      for(unsigned j = 0; j < outDegPortal.Get(i); j++)
//      {
//        std::cout << MaskedIndex(nbhdmaPortal.Get(i) & (static_cast<vtkm::Id>(1) << j)) << " ";
//      }
//      std::cout << "\n";
//  }

//  // next, we compute where each vertex lands in the new array
//  // it needs to be one place offset, hence the +/- 1
//  // this should automatically parallelise
//  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
//  /*auto oneIfCritical = [](unsigned x) { return x!= 1 ? 1 : 0; };

//    // we need a temporary inverse index to change vertex IDs
//    IdArrayType inverseIndex;
//    inverseIndex.Allocate(mesh.NumVertices);
//    inverseIndex.WritePortal().Set(0,0);

//    std::partial_sum(
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(outDegrees.WritePortal()), oneIfCritical),
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(outDegrees.WritePortal())-1, oneIfCritical),
//            vtkm::cont::ArrayPortalToIteratorBegin(inverseIndex.WritePortal()) + 1);
//    */
//  IdArrayType inverseIndex;
//  OneIfCritical oneIfCriticalFunctor;
//  auto oneIfCriticalArrayHandle =
//    vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfCritical>(outDegrees, oneIfCriticalFunctor);
//  vtkm::cont::Algorithm::ScanExclusive(oneIfCriticalArrayHandle, inverseIndex);

//  // now we can compute how many critical points we carry forward
//  vtkm::Id nCriticalPoints =
//    this->GetLastValue(inverseIndex) + oneIfCriticalFunctor(this->GetLastValue(outDegrees));

//  // we need to keep track of what the index of each vertex is in the active graph
//  // for most vertices, this should have the NO_SUCH_VERTEX flag set
//  AllocateVertexArrays(
//    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices

//  // our processing now depends on the degree of the vertex
//  // but basically, we want to set up the arrays for this vertex:
//  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
//  // GlobalIndex stores the index in the join tree for later access
//  IdArrayType activeIndices;
//  activeIndices.Allocate(mesh.NumVertices);
//  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
//    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
//  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);

//  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
//  this->Invoke(initActiveGraphVerticesWorklet,
//               sortIndexArray,
//               outDegrees,
//               inverseIndex,
//               extrema,
//               activeIndices,
//               this->GlobalIndex,
//               this->Outdegree,
//               this->Hyperarcs,
//               this->ActiveVertices);

//  // now we need to compute the FirstEdge array from the outDegrees
//  this->FirstEdge.Allocate(nCriticalPoints);
//  // STD Version of the prefix sum
//  //this->FirstEdge.WritePortal().Set(0, 0);
//  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
//  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
//  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
//  // VTKM Version of the prefix sum
//  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
//  // Compute the number of critical edges

//  vtkm::Id nCriticalEdges =
//    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

//  AllocateEdgeArrays(nCriticalEdges);

//  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
//  this->Invoke(initActiveEdgesWorklet,
//               this->Outdegree,
//               mesh,
//               this->FirstEdge,
//               this->GlobalIndex,
//               extrema,
//               neighbourhoodMasks,
//               this->EdgeNear,
//               this->EdgeFar,
//               this->ActiveEdges);

//  // now we have to go through and set the far ends of the new edges using the
//  // inverse index array
//  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
//  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

//  DebugPrint("Active Graph Started", __FILE__, __LINE__);

//  // then we loop through the active vertices to convert their indices to active graph indices
//  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
//  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);

//  // finally, allocate and initialise the edgeSorter array
//  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
//  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);

//  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);
//} // InitialiseActiveGraph()















// initialises the active graph - FROM THE OFFICIAL VTK-M VERSION!
// ... MODIFICATION ABOVE TO GET THE MESH WITH ALL VERTICES MARKED AS CRITICAL !!!
template <class Mesh>
inline void ActiveGraph::Initialise(Mesh& mesh, const MeshExtrema& meshExtrema)
{ // InitialiseActiveGraph()
  // reference to the correct array in the extrema
  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

  // For every vertex, work out whether it is critical
  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
  // All vertices of outdegree 0 must be extrema
  // Saddle points must be at least outdegree 2, so this is a correct test
  // BUT it is possible to overestimate the degree of a non-extremum,
  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph

  // Neighbourhood mask (one bit set per connected component in neighbourhood
  IdArrayType neighbourhoodMasks;
  neighbourhoodMasks.Allocate(mesh.NumVertices);
  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
  outDegrees.Allocate(mesh.NumVertices);

  // Initialize the nerighborhoodMasks and outDegrees arrays
  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
    this->IsJoinGraph);

  this->Invoke(initNeighMasksAndOutDegWorklet,
               sortIndexArray,
               mesh,
               neighbourhoodMasks, // output
               outDegrees);        // output

  // next, we compute where each vertex lands in the new array
  // it needs to be one place offset, hence the +/- 1
  // this should automatically parallelise
  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
  /*auto oneIfCritical = [](unsigned x) { return x!= 1 ? 1 : 0; };

    // we need a temporary inverse index to change vertex IDs
    IdArrayType inverseIndex;
    inverseIndex.Allocate(mesh.NumVertices);
    inverseIndex.WritePortal().Set(0,0);

    std::partial_sum(
            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(outDegrees.WritePortal()), oneIfCritical),
            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(outDegrees.WritePortal())-1, oneIfCritical),
            vtkm::cont::ArrayPortalToIteratorBegin(inverseIndex.WritePortal()) + 1);
    */
  IdArrayType inverseIndex;
  OneIfCritical oneIfCriticalFunctor;
  auto oneIfCriticalArrayHandle =
    vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfCritical>(outDegrees, oneIfCriticalFunctor);
  vtkm::cont::Algorithm::ScanExclusive(oneIfCriticalArrayHandle, inverseIndex);

  // now we can compute how many critical points we carry forward
  vtkm::Id nCriticalPoints =
    this->GetLastValue(inverseIndex) + oneIfCriticalFunctor(this->GetLastValue(outDegrees));

  // we need to keep track of what the index of each vertex is in the active graph
  // for most vertices, this should have the NO_SUCH_VERTEX flag set
  AllocateVertexArrays(
    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices

  // our processing now depends on the degree of the vertex
  // but basically, we want to set up the arrays for this vertex:
  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
  // GlobalIndex stores the index in the join tree for later access
  IdArrayType activeIndices;
  activeIndices.Allocate(mesh.NumVertices);
  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);

  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
  this->Invoke(initActiveGraphVerticesWorklet,
               sortIndexArray,
               outDegrees,
               inverseIndex,
               extrema,
               activeIndices,
               this->GlobalIndex,
               this->Outdegree,
               this->Hyperarcs,
               this->ActiveVertices);

  // now we need to compute the FirstEdge array from the outDegrees
  this->FirstEdge.Allocate(nCriticalPoints);
  // STD Version of the prefix sum
  //this->FirstEdge.WritePortal().Set(0, 0);
  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
  // VTKM Version of the prefix sum
  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
  // Compute the number of critical edges

  vtkm::Id nCriticalEdges =
    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

  AllocateEdgeArrays(nCriticalEdges);

  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
  this->Invoke(initActiveEdgesWorklet,
               this->Outdegree,
               mesh,
               this->FirstEdge,
               this->GlobalIndex,
               extrema,
               neighbourhoodMasks,
               this->EdgeNear,
               this->EdgeFar,
               this->ActiveEdges);

  // WAIT FOR MESH SLEEP
  /// DEBUG PRINT std::cout << "Check the MeshOutput file ... \n";
//  std::this_thread::sleep_for(std::chrono::seconds(3));

  // now we have to go through and set the far ends of the new edges using the
  // inverse index array
  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

  DebugPrint("Active Graph Started", __FILE__, __LINE__);

  // then we loop through the active vertices to convert their indices to active graph indices
  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);

  // finally, allocate and initialise the edgeSorter array
  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);

  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);
} // InitialiseActiveGraph()





//// started refactoring 2024-03-15
//// ... first by removing everything that we don't need
//// ... just doing manual active graph computation at this point!

//// initialises the active graph
//template <class Mesh>
//inline void ActiveGraph::Initialise(Mesh& mesh, const MeshExtrema& meshExtrema)
//{ // InitialiseActiveGraph()
//  // reference to the correct array in the extrema

//  const char* message = "debug";
//  long workaround = 100;
//  std::cout << "--------------- extrema.DebugPrint(message, message, 100); --------------- \n";
//  printf("ActiveGraph::Initialise() on: %s\n", __PRETTY_FUNCTION__);
////  meshExtrema.DebugPrint(message, message, workaround);
//  PrintHeader(meshExtrema.Peaks.GetNumberOfValues());
//  PrintIndices("MeshExtremaPeaks", meshExtrema.Peaks);
//  PrintIndices("MeshExtremaPits", meshExtrema.Pits);
//  std::cout << "A1\n";
//  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;
//  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;

//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
//      PrintIndices("MeshExtremaPeaksAFTER", extrema);
//  }
//  std::cout << "------------------------------- fine until here -------------------------------\n";


////  / *
//  // For every vertex, work out whether it is critical
//  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
//  // All vertices of outdegree 0 must be extrema
//  // Saddle points must be at least outdegree 2, so this is a correct test
//  // BUT it is possible to overestimate the degree of a non-extremum,
//  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph
//  std::cout << "A2\n";
//  // Neighbourhood mask (one bit set per connected component in neighbourhood
////  IdArrayType neighbourhoodMasks;
////  neighbourhoodMasks.Allocate(mesh.NumVertices);
////  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
////  outDegrees.Allocate(mesh.NumVertices);
//  std::cout << "A3\n";
//  // Initialize the nerighborhoodMasks and outDegrees arrays
////  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
////  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
////  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
////    this->IsJoinGraph);
//  std::cout << "A4\n";
//  //  invoker-disabled-2024-03-18:
////  this->Invoke(initNeighMasksAndOutDegWorklet,
////               sortIndexArray,
////               mesh,
////               neighbourhoodMasks, // output
////               outDegrees);        // output

//  std::cout << "------------------------------- fine until here -------------------------------\n";
//  std::cout << "(Just initialised sort indices and outdegs up there...)\n";
//  std::cout << "------------------------------- fine until here -------------------------------\n";
//  std::cout << "A5\n";
//  IdArrayType inverseIndex;
//  inverseIndex.Allocate(mesh.NumVertices);
//  using IdPortalType = IdArrayType::WritePortalType;
//  IdPortalType inversePortal = inverseIndex.WritePortal();
//  for(vtkm::Id i = 0; i < mesh.NumVertices; i++)
//  {
//      inversePortal.Set(i, i);
//  }

//  std::cout << "A6\n";
//  vtkm::Id nCriticalPoints = mesh.NumVertices;
//  std::cout << "A7\n";
//  // we need to keep track of what the index of each vertex is in the active graph
//  // for most vertices, this should have the NO_SUCH_VERTEX flag set
//  AllocateVertexArrays(
//    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices
//  std::cout << "A8\n";
//  // our processing now depends on the degree of the vertex
//  // but basically, we want to set up the arrays for this vertex:
//  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
//  // GlobalIndex stores the index in the join tree for later access
//  IdArrayType activeIndices;
//  activeIndices.Allocate(mesh.NumVertices);
//  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
//    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
//  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);
//  std::cout << "A9\n";

//  // print removed 2024-03-18
////  auto readSortIdx =  sortIndexArray.ReadPortal();
////  for(unsigned i = 0; i < mesh.NumVertices; i++)
////  {
////      std::cout << readSortIdx.Get(i) << " ";
////  }
////  std::cout << "\n";

////  auto readNhbOffsets =  mesh.NeighborOffsets.ReadPortal();
////  PrintIndices(mesh.NeighborOffsets);
//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
////      std::cout << readNhbOffsets.Get(i) << " ";
//  }


////  PrintIndices("FOR_INIT:sortIndexArray",    sortIndexArray);
////  PrintIndices("FOR_INIT:outDegrees",        outDegrees);
//  PrintIndices("FOR_INIT:inverseIndex",      inverseIndex);
//  PrintIndices("FOR_INIT:extrema",           extrema);
//  PrintIndices("FOR_INIT:activeIndices",     activeIndices);

////  invoker-disabled-2024-03-18:
////  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
////  this->Invoke(initActiveGraphVerticesWorklet,
////               sortIndexArray,
////               outDegrees,
////               inverseIndex,
////               extrema,
////               activeIndices,
////               this->GlobalIndex,
////               this->Outdegree,
////               this->Hyperarcs,
////               this->ActiveVertices);
//  std::cout << "A10\n";
//  PrintIndices("INIT:activeIndices",     activeIndices);
//  PrintIndices("INIT:GlobalIndex",    this->GlobalIndex);
//  PrintIndices("INIT:Outdegree",      this->Outdegree);
//  PrintIndices("INIT:Hyperarcs",      this->Hyperarcs);
//  PrintIndices("INIT:ActiveVertices", this->ActiveVertices);
////  * /

//  IdPortalType aixPortal = activeIndices.WritePortal();
//  IdPortalType avPortal  = this->ActiveVertices.WritePortal();
//  IdPortalType haPortal  = this->Hyperarcs.WritePortal();
//  IdPortalType odPortal  = this->Outdegree.WritePortal();
//  IdPortalType giPortal  = this->GlobalIndex.WritePortal();

//  auto exPortal = extrema.ReadPortal();

//  vtkm::Id testactiveidx[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,15,0,0,18,19,0,0,22,0,24};

//  vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,1,2,1,1,3,2,1,1,0,1,0}; // bare
////vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,2,3,2,0,0}; // optimised

//  vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,24,22,22,24,24,22,24,22,22,24,24}; // bare
////vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,  22,   24,22,22,24,24,22,24,22,22,24,24}; // optimised

//  vtkm::Id testoutdeg2[] = {0,0,1,1,2,2,2,2,0,1,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // bare
////vtkm::Id testoutdeg2[] = {0,0,2,2,2,2,0,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // optimised

//  vtkm::Id testextrem2[] = {0,1,0,1,0,0,1,1,8,0,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // bare
////vtkm::Id testextrem2[] = {0,1,0,0,1,1,8,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // optimised

//  // hacking around:
//  if(up_or_down)
//  {
//      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//      {
//          avPortal.Set(i, i);
//    //      haPortal.Set(i, (vtkm::Id)exPortal.Get(i));
//          haPortal.Set(i, testextrem1[i]);
//    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
//          odPortal.Set(i, testoutdeg1[i]);
//          giPortal.Set(i, i);

//          // 2024-03-18:
//          aixPortal.Set(i, testactiveidx[i]);
//      }
//  }
//  else
//  {
//      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//      {
//          avPortal.Set(i, i);
//    //      haPortal.Set(i, (vtkm::Id)exPortal.Get(i));
//          haPortal.Set(i, testextrem2[i]);
//    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
//          odPortal.Set(i, testoutdeg2[i]);
//          giPortal.Set(i, i);


//          // 2024-03-18:
//          aixPortal.Set(i, testactiveidx[i]);
//      }
//  }

//  PrintIndices("VX:ActiveIndices",  activeIndices);
//  PrintIndices("VX:GlobalIndex",    this->GlobalIndex);
//  PrintIndices("VX:Outdegree",      this->Outdegree);
//  PrintIndices("VX:Hyperarcs",      this->Hyperarcs);
//  PrintIndices("VX:ActiveVertices", this->ActiveVertices);
//  // now we need to compute the FirstEdge array from the outDegrees
//  this->FirstEdge.Allocate(nCriticalPoints);
//  // STD Version of the prefix sum
//  //this->FirstEdge.WritePortal().Set(0, 0);
//  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
//  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
//  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
//  // VTKM Version of the prefix sum
//  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
//  // Compute the number of critical edges

//  vtkm::Id nCriticalEdges =
//    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

////  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;
//  if(up_or_down)
//  {
//      nCriticalEdges = 56;
//  }
//  else
//  {
////      nCriticalEdges = 53;
//      nCriticalEdges = 56;
//  }

//  AllocateEdgeArrays(nCriticalEdges);
//  std::cout << "A11\n";
//  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
//  std::cout << "vvv ALL ACTIVE EDGES INITIALIZED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES INITIALIZED! ^^^\n";
//  std::cout << "A11.1\n";
//  // try doing this manually:
//  std::cout << "initActiveEdgesWorkled params: Outdegree=\n";
////  << this->Outdegree << ", GlobalIndex=" << this->GlobalIndex << "\n";
////  PrintIndices("outDegrees?: ", outDegrees);
//  PrintIndices("OUTDEGREE?: ", this->Outdegree);
//  PrintIndices("GlobalIDX(sort)?: ", this->GlobalIndex);
//  std::cout << "nCriticalEdges?: " << nCriticalEdges << "\n";

//  // currently using REAL IDs (trying to make ALL vertices active)
//  vtkm::Id low_arr1[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,15,15,16,17,18,18,18,19,19,20,21,23};
//  // ACTIVE IDS:
//                      //{0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,14,14,14,14,15,15,15,16,16,16,16,17};
//                      //{0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
//  vtkm::Id up_arr1[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,24,22,22,24,24,24,24,22,22,24,22,24};
////  {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
////  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
//  // unoptimised:
//  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
//  int num1 = 56; //50;

//  // currently using REAL IDs (trying to make ALL vertices active)
//  vtkm::Id low_arr2[] = {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  // shortened 53:
////  vtkm::Id low_arr2[] = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  // peaks as pits: {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,15,15,15,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,18,18,18};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//          //{2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  vtkm::Id up_arr2[] = {0,1,0,0,0,0,1,1,1,1,0,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  // shortened 53
////  vtkm::Id up_arr2[] = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //{0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //shortened 53
//  //  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
//  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52, 53,54,55};
//  int num2 = 56;

//  int num=0;

//  std::vector<vtkm::Id> indices;
//  std::vector<vtkm::Id> low_arr;
//  std::vector<vtkm::Id> up_arr;

//  IdPortalType aePortal = this->ActiveEdges.WritePortal();
//  IdPortalType enPortal = this->EdgeNear.WritePortal();
//  IdPortalType efPortal = this->EdgeFar.WritePortal();

//  if(up_or_down)
//  {
//      this->ActiveEdges.Allocate(num1);
//      this->EdgeNear.Allocate(num1);
//      this->EdgeFar.Allocate(num1);

//      for(vtkm::Id i = 0; i < num1; i++)
//      {
//        indices.push_back(indices1[i]);
//        low_arr.push_back(low_arr1[i]);
//        up_arr.push_back(up_arr1[i]);

//        aePortal.Set(i, indices1[i]);
//        enPortal.Set(i, low_arr1[i]);
//        efPortal.Set(i, up_arr1[i]);
//      }
//  }
//  else
//  {
//      this->ActiveEdges.Allocate(num2);
//      this->EdgeNear.Allocate(num2);
//      this->EdgeFar.Allocate(num2);

//      for(vtkm::Id i = 0; i < num2; i++)
//      {
//        indices.push_back(indices2[i]);
//        low_arr.push_back(low_arr2[i]);
//        up_arr.push_back(up_arr2[i]);

//        aePortal.Set(i, indices2[i]);
//        enPortal.Set(i, low_arr2[i]);
//        efPortal.Set(i, up_arr2[i]);

//      }
//  }

//  std::cout << "A11.1 FINISHED!\n";
//  if(up_or_down)
//  {
//    std::cout << "UP\n";
//  }
//  else
//  {
//    std::cout << "DOWN\n";
//  }

//  std::cout << "vvv ALL ACTIVE EDGES GATHERED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES GATHERED! ^^^\n";

//  std::cout << "A11.2\n";
//  // now we have to go through and set the far ends of the new edges using the
//  // inverse index array
//  std::cout << "A11.3\n";
//  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
//  std::cout << "A11.4\n";
//  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

//  DebugPrint("Active Graph Started", __FILE__, __LINE__);
//  std::cout << "A12\n";
//  // then we loop through the active vertices to convert their indices to active graph indices
//  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
//  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);
//  std::cout << "A13\n";
//  // finally, allocate and initialise the edgeSorter array
//  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
//  std::cout << "A14\n";
//  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);
//  std::cout << "ActiveGraph - finished\n";
//  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);
//} // InitialiseActiveGraph()
// UNCOMMENT MANUAL































































































































//// UNCOMMMENT for fixed 64-nbor fix - MESSY, UNCLEANED VERSION AS OF 2024-03-21

//// started refactoring 2024-03-15
//// ... first by removing everything that we don't need
//// ... just doing manual active graph computation at this point!


//// C++ Partial Template Specialization

//// initialises the active graph
////template </*using C++ Partial Template Specialization, leave empty */ >
////template <>
////inline void ActiveGraph::Initialise<ContourTreeMesh>(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema)
////inline void ActiveGraph::Initialise(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema)

////inline void ActiveGraph::Initialise(ContourTreeMesh& mesh, const MeshExtrema& meshExtrema)
//inline void ActiveGraph::Initialise(vtkm::worklet::contourtree_augmented::ContourTreeMesh<int>& mesh,
//                                    const MeshExtrema& meshExtrema)
//{ // InitialiseActiveGraph()
//  // reference to the correct array in the extrema

//  const char* message = "debug";
//  long workaround = 100;
//  std::cout << "--------------- extrema.DebugPrint(message, message, 100); --------------- \n";
//  std::cout << "C++ Partial Template Specialization\n";
////  meshExtrema.DebugPrint(message, message, workaround);
//  #if DEBUG_PRINT_PACTBD
//      PrintHeader(meshExtrema.Peaks.GetNumberOfValues());
//      PrintIndices("MeshExtremaPeaks", meshExtrema.Peaks);
//      PrintIndices("MeshExtremaPits", meshExtrema.Pits);
//  #endif



//  std::cout << "A1\n";
//  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;
//  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;

//  #if DEBUG_PRINT_PACTBD
//      for(unsigned i = 0; i < mesh.NumVertices; i++)
//      {
//          PrintIndices("mesh.NeighborConnectivity", mesh.NeighborConnectivity);
//          PrintIndices("mesh.NeighborOffsets", mesh.NeighborOffsets);
//          PrintIndices("MeshExtremaPeaksAFTER", extrema);
//      }
//      std::cout << "------------------------------- fine until here -------------------------------\n";
//  #endif


////  / *
//  // For every vertex, work out whether it is critical
//  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
//  // All vertices of outdegree 0 must be extrema
//  // Saddle points must be at least outdegree 2, so this is a correct test
//  // BUT it is possible to overestimate the degree of a non-extremum,
//  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph
//  std::cout << "A2\n";
//  // Neighbourhood mask (one bit set per connected component in neighbourhood
////  IdArrayType neighbourhoodMasks;
////  neighbourhoodMasks.Allocate(mesh.NumVertices);
////  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
////  outDegrees.Allocate(mesh.NumVertices);
//  std::cout << "A3\n";
//  // Initialize the nerighborhoodMasks and outDegrees arrays
////  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
////  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
////  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
////    this->IsJoinGraph);
//  std::cout << "A4\n";
//  //  invoker-disabled-2024-03-18:
////  this->Invoke(initNeighMasksAndOutDegWorklet,
////               sortIndexArray,
////               mesh,
////               neighbourhoodMasks, // output
////               outDegrees);        // output

//  #if DEBUG_PRINT_PACTBD
//      std::cout << "------------------------------- fine until here -------------------------------\n";
//      std::cout << "(Just initialised sort indices and outdegs up there...)\n";
//      std::cout << "------------------------------- fine until here -------------------------------\n";
//  #endif
//  std::cout << "A5\n";
//  IdArrayType inverseIndex;
//  inverseIndex.Allocate(mesh.NumVertices);
//  using IdPortalType = IdArrayType::WritePortalType;
//  IdPortalType inversePortal = inverseIndex.WritePortal();
//  for(vtkm::Id i = 0; i < mesh.NumVertices; i++)
//  {
//      inversePortal.Set(i, i);
//  }

//  std::cout << "A6\n";
//  vtkm::Id nCriticalPoints = mesh.NumVertices;
//  std::cout << "A7\n";
//  // we need to keep track of what the index of each vertex is in the active graph
//  // for most vertices, this should have the NO_SUCH_VERTEX flag set
//  AllocateVertexArrays(
//    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices
//  std::cout << "A8\n";
//  // our processing now depends on the degree of the vertex
//  // but basically, we want to set up the arrays for this vertex:
//  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
//  // GlobalIndex stores the index in the join tree for later access
//  IdArrayType activeIndices;
//  activeIndices.Allocate(mesh.NumVertices);
//  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
//    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
//  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);
//  std::cout << "A9\n";

//  // print removed 2024-03-18
////  auto readSortIdx =  sortIndexArray.ReadPortal();
////  for(unsigned i = 0; i < mesh.NumVertices; i++)
////  {
////      std::cout << readSortIdx.Get(i) << " ";
////  }
////  std::cout << "\n";

////  auto readNhbOffsets =  mesh.NeighborOffsets.ReadPortal();
////  PrintIndices(mesh.NeighborOffsets);
//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
////      std::cout << readNhbOffsets.Get(i) << " ";
//  }


////  PrintIndices("FOR_INIT:sortIndexArray",    sortIndexArray);
////  PrintIndices("FOR_INIT:outDegrees",        outDegrees);
//  #if DEBUG_PRINT_PACTBD
//      PrintIndices("FOR_INIT:inverseIndex",      inverseIndex);
//      PrintIndices("FOR_INIT:extrema",           extrema);
//      PrintIndices("FOR_INIT:activeIndices",     activeIndices);
//  #endif

////  invoker-disabled-2024-03-18:
////  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
////  this->Invoke(initActiveGraphVerticesWorklet,
////               sortIndexArray,
////               outDegrees,
////               inverseIndex,
////               extrema,
////               activeIndices,
////               this->GlobalIndex,
////               this->Outdegree,
////               this->Hyperarcs,
////               this->ActiveVertices);
//  std::cout << "A10\n";
//  #if DEBUG_PRINT_PACTBD
//      PrintIndices("INIT:activeIndices",     activeIndices);
//      PrintIndices("INIT:GlobalIndex",    this->GlobalIndex);
//      PrintIndices("INIT:Outdegree",      this->Outdegree);
//      PrintIndices("INIT:Hyperarcs",      this->Hyperarcs);
//      PrintIndices("INIT:ActiveVertices", this->ActiveVertices);
//  #endif
////  * /

//  IdPortalType aixPortal = activeIndices.WritePortal();         // 1
//  IdPortalType avPortal  = this->ActiveVertices.WritePortal();  // 2
//  IdPortalType haPortal  = this->Hyperarcs.WritePortal();       // 3
//  IdPortalType odPortal  = this->Outdegree.WritePortal();       // 4
//  IdPortalType giPortal  = this->GlobalIndex.WritePortal();     // 5

////  auto exPortal = extrema.ReadPortal();
//  const auto exPortal  = extrema.ReadPortal();
//  const auto outReadPortal = this->Outdegree.ReadPortal();

////  vtkm::Id testactiveidx[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,15,0,0,18,19,0,0,22,0,24};

//  // compute from 'outDegrees' and 'nbhoodConnectivity'
////  vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,1,2,1,1,3,2,1,1,0,1,0}; // bare
////vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,2,3,2,0,0}; // optimised

//  // take from 'extrema'
////  vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,24,22,22,24,24,22,24,22,22,24,24}; // bare
////vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,  22,   24,22,22,24,24,22,24,22,22,24,24}; // optimised

//  // compute from 'outDegrees' and 'nbhoodConnectivity'
////  vtkm::Id testoutdeg2[] = {0,0,1,1,2,2,2,2,0,1,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // bare
////vtkm::Id testoutdeg2[] = {0,0,2,2,2,2,0,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // optimised

//  // take from 'extrema'
////  vtkm::Id testextrem2[] = {0,1,0,1,0,0,1,1,8,0,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // bare
////vtkm::Id testextrem2[] = {0,1,0,0,1,1,8,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // optimised


//  const auto meshNbhoodConPortal = mesh.NeighborConnectivity.ReadPortal();
//  const auto meshNbhoodOffPortal = mesh.NeighborOffsets.ReadPortal();

//  std::vector<vtkm::Id> farEndsUP;
//  std::vector<vtkm::Id> farEndsDOWN;

//  // compute up/down outdegrees degrees
//  for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//  {
//      avPortal.Set(i, i);
//      giPortal.Set(i, i);
//      haPortal.Set(i, MaskedIndex(exPortal.Get(i)) );
//      // 2024-03-18 original aixPortal:
//      // aixPortal.Set(i, testactiveidx[i]);

//      // 2024-03-20 update:
////      if(MaskedIndex(exPortal.Get(i)) != 1)


////      std::cout << outReadPortal.Get(i) << "\n";

//      if(outReadPortal.Get(i) != 1)
//      {
//        aixPortal.Set(i, i);
//      }
//      else
//      {
//        aixPortal.Set(i, 0);
//      }


//      vtkm::Id start_interval = meshNbhoodOffPortal.Get(i);
//      vtkm::Id end_interval   = meshNbhoodOffPortal.Get(i+1);

//      int conditional_count_up   = 0;
//      int conditional_count_down = 0;

//      for(vtkm::Id j = start_interval; j < end_interval; j++)
//      {
//        if(meshNbhoodConPortal.Get(j) > i)
//        {
//            conditional_count_up++;
//            farEndsUP.push_back(meshNbhoodConPortal.Get(j));
//        }
//        else
//        {
//            conditional_count_down++;
//            farEndsDOWN.push_back(meshNbhoodConPortal.Get(j));
//        }
//      }

//      if(up_or_down)
//      {
////          std::cout << conditional_count_up << "\n";
//          odPortal.Set(i, conditional_count_up);
//      }
//      else
//      {
////          std::cout << conditional_count_down << "\n";
//          odPortal.Set(i, conditional_count_down);
//      }

//  }


//  // sort out active indices after out degrees are computed!
////  std::cout << "!!!!!!!!!compute up/down outdegrees degrees\n";
//  for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//  {
////      std::cout << outReadPortal.Get(i) << "\n";

//      if(outReadPortal.Get(i) != 1)
//      {
//        aixPortal.Set(i, i);
//      }
//      else
//      {
//        aixPortal.Set(i, 0);
//      }
//  }

////  // hacking around:
////  if(up_or_down)
////  {
////      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
////      {

////          std::cout << "exPortal.GetUp(i):" << MaskedIndex(exPortal.Get(i)) << "\n";

////          // NOTE TO SELF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////          // ... IF PRINTINDEX ADDS .t.... THEN IT IS ACTUALLY A MASKEDINDEX NUMBER!!!!!!!!!


////          avPortal.Set(i, i);
////          haPortal.Set(i, MaskedIndex(exPortal.Get(i)) ); // !!! Use MaskedIndex !!!
//////          haPortal.Set(i, testextrem1[i]); // should be got from 'exPortal' - same as extrema!
////    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
////          odPortal.Set(i, testoutdeg1[i]);
////          giPortal.Set(i, i);

////          // 2024-03-18:
////          aixPortal.Set(i, testactiveidx[i]);
////      }
//////      std::cout << "\n";
////  }
////  else
////  {
////      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
////      {
////          std::cout << "exPortal.GetDown(i):" << MaskedIndex(exPortal.Get(i)) << "\n";

////          avPortal.Set(i, i);
////          haPortal.Set(i, MaskedIndex(exPortal.Get(i)) );

//////          haPortal.Set(i, testextrem2[i]);
////    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
////          odPortal.Set(i, testoutdeg2[i]);
////          giPortal.Set(i, i);


////          // 2024-03-18:
////          aixPortal.Set(i, testactiveidx[i]);
////      }
////  }

//  #if DEBUG_PRINT_PACTBD
//      PrintIndices("VX:ActiveIndices",  activeIndices);
//      PrintIndices("VX:GlobalIndex",    this->GlobalIndex);
//      PrintIndices("VX:Outdegree",      this->Outdegree);
//      PrintIndices("VX:Hyperarcs",      this->Hyperarcs);
//      PrintIndices("VX:ActiveVertices", this->ActiveVertices);
//  #endif
//  // now we need to compute the FirstEdge array from the outDegrees
//  this->FirstEdge.Allocate(nCriticalPoints);
//  // STD Version of the prefix sum
//  //this->FirstEdge.WritePortal().Set(0, 0);
//  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
//  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
//  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
//  // VTKM Version of the prefix sum
//  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
//  // Compute the number of critical edges

//  vtkm::Id nCriticalEdges =
//    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

////  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;
////  if(up_or_down)
////  {
////      nCriticalEdges = 56;
////  }
////  else
////  {
//////      nCriticalEdges = 53;
////      nCriticalEdges = 56;
////  }

//  AllocateEdgeArrays(nCriticalEdges);
//  std::cout << "A11\n";
////  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
//  active_graph_inc_ns::InitializeActiveEdges<vtkm::worklet::contourtree_augmented::ContourTreeMesh<int>> initActiveEdgesWorklet;

//  #if DEBUG_PRINT_PACTBD
//      std::cout << "vvv ALL ACTIVE EDGES INITIALIZED! vvv\n";
//      PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//      PrintIndices("NEAR", this->EdgeNear);
//      PrintIndices("FAR", this->EdgeFar);
//      std::cout << "^^^ ALL ACTIVE EDGES INITIALIZED! ^^^\n";
//      std::cout << "A11.1\n";
//      // try doing this manually:
//      std::cout << "initActiveEdgesWorkled params: Outdegree=\n";
//    //  << this->Outdegree << ", GlobalIndex=" << this->GlobalIndex << "\n";
//    //  PrintIndices("outDegrees?: ", outDegrees);
//      PrintIndices("OUTDEGREE?: ", this->Outdegree);
//      PrintIndices("GlobalIDX(sort)?: ", this->GlobalIndex);
//      std::cout << "nCriticalEdges?: " << nCriticalEdges << "\n";
//  #endif

//  // currently using REAL IDs (trying to make ALL vertices active)
////  vtkm::Id low_arr1[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,15,15,16,17,18,18,18,19,19,20,21,23};
//  // ACTIVE IDS:
//                      //{0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,14,14,14,14,15,15,15,16,16,16,16,17};
//                      //{0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////  vtkm::Id up_arr1[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,24,22,22,24,24,24,24,22,22,24,22,24};
////  {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
////  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
//  // unoptimised:
////  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
////  int num1 = 56; //50;

//  // currently using REAL IDs (trying to make ALL vertices active)
////  vtkm::Id low_arr2[] = {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  // shortened 53:
////  vtkm::Id low_arr2[] = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  // peaks as pits: {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,15,15,15,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,18,18,18};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//          //{2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
////  vtkm::Id up_arr2[] = {0,1,0,0,0,0,1,1,1,1,0,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  // shortened 53
////  vtkm::Id up_arr2[] = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //{0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //shortened 53
//  //  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
////  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52, 53,54,55};
////  int num2 = 56;

////  int num=0;

//  std::vector<vtkm::Id> indices;
//  std::vector<vtkm::Id> low_arr;
//  std::vector<vtkm::Id> up_arr;

//  IdPortalType aePortal = this->ActiveEdges.WritePortal();
//  IdPortalType enPortal = this->EdgeNear.WritePortal();
//  IdPortalType efPortal = this->EdgeFar.WritePortal();

//  auto readOutDegs = this->Outdegree.ReadPortal();


////  std::cout << "Far Ends test:\n";
//  unsigned total = 0;
//  if(up_or_down)
//  {
//      this->ActiveEdges.Allocate(farEndsUP.size());
//      this->EdgeNear.Allocate(farEndsUP.size());
//      this->EdgeFar.Allocate(farEndsUP.size());

//      for(unsigned i = 0; i < farEndsUP.size(); i++)
//      {
////        std::cout << i << ": " << farEndsUP[i] << " -> " << MaskedIndex(exPortal.Get(farEndsUP[i])) << "\n";
////        std::cout << MaskedIndex(exPortal.Get(farEndsUP[i])) << "\n";

//        aePortal.Set(i, i);
//        efPortal.Set(i, MaskedIndex(exPortal.Get(farEndsUP[i])) ); // replaced 2024-03-20
//      }

//      total = 0;
//      for(unsigned i = 0; i < mesh.NumVertices; i++)
//      {
//          for (unsigned j = 0; j < readOutDegs.Get(i); j++)
//          {
//              enPortal.Set(total, i);
//              total++;
//          }
//      }
//  }
//  else
//  {
//      this->ActiveEdges.Allocate(farEndsDOWN.size());
//      this->EdgeNear.Allocate(farEndsDOWN.size());
//      this->EdgeFar.Allocate(farEndsDOWN.size());

//      for(unsigned i = 0; i < farEndsDOWN.size(); i++)
//      {
////        std::cout << i << ": " << farEndsDOWN[i] << " -> " << MaskedIndex(exPortal.Get(farEndsDOWN[i])) << "\n";
////          std::cout << MaskedIndex(exPortal.Get(farEndsDOWN[i])) << "\n";

//          aePortal.Set(i, i);
//          efPortal.Set(i, MaskedIndex(exPortal.Get(farEndsDOWN[i])) ); // replaced 2024-03-20
//      }

//      total = 0;
//      for(unsigned i = 0; i < mesh.NumVertices; i++)
//      {
//          for (unsigned j = 0; j < readOutDegs.Get(i); j++)
//          {
//              enPortal.Set(total, i);
//              total++;
//          }
//      }
//  }


////  if(up_or_down)
////  {
//////      this->ActiveEdges.Allocate(num1);
//////      this->EdgeNear.Allocate(num1);
//////      this->EdgeFar.Allocate(num1);

////      for(vtkm::Id i = 0; i < num1; i++)
////      {
//////        indices.push_back(indices1[i]);
////        low_arr.push_back(low_arr1[i]);
////        up_arr.push_back(up_arr1[i]);

//////        aePortal.Set(i, indices1[i]);
//////        aePortal.Set(i, i);            // replaced 2024-03-20
//////        enPortal.Set(i, low_arr1[i]); // replaced later on 2024-03-20
//////        efPortal.Set(i, up_arr1[i]); // replaced 2024-03-20
////      }
////  }
////  else
////  {
//////      this->ActiveEdges.Allocate(num2);
//////      this->EdgeNear.Allocate(num2);
//////      this->EdgeFar.Allocate(num2);

////      for(vtkm::Id i = 0; i < num2; i++)
////      {
//////        indices.push_back(indices2[i]);
////        low_arr.push_back(low_arr2[i]);
////        up_arr.push_back(up_arr2[i]);

//////        aePortal.Set(i, indices2[i]);
//////        aePortal.Set(i, i);            // replaced 2024-03-20
//////        enPortal.Set(i, low_arr2[i]);  // replaced later 2024-03-20
//////        efPortal.Set(i, up_arr2[i]); // replaced 2024-03-20

////      }
////  }

//  std::cout << "A11.1 FINISHED!\n";
//  if(up_or_down)
//  {
//    std::cout << "UP\n";
//  }
//  else
//  {
//    std::cout << "DOWN\n";
//  }

//  #if DEBUG_PRINT_PACTBD
//        printf("INVOKER() on: %s\n", __PRETTY_FUNCTION__);
//        std::cout << "Invoker3\n";
//      std::cout << "vvv ALL ACTIVE EDGES GATHERED! vvv\n";
//      PrintIndices("ACTIVE EDGES", this->ActiveEdges);          // 6
//      PrintIndices("NEAR", this->EdgeNear);                     // 7
//      PrintIndices("FAR", this->EdgeFar);                       // 8
//      std::cout << "^^^ ALL ACTIVE EDGES GATHERED! ^^^\n";
//  #endif

//  std::cout << "A11.2\n";
//  // now we have to go through and set the far ends of the new edges using the
//  // inverse index array
//  std::cout << "A11.3\n";
//  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
//  std::cout << "A11.4\n";
//  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

//  DebugPrint("Active Graph Started", __FILE__, __LINE__);
//  std::cout << "A12\n";
//  // then we loop through the active vertices to convert their indices to active graph indices
//  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
//  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);
//  std::cout << "A13\n";
//  // finally, allocate and initialise the edgeSorter array
//  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
//  std::cout << "A14\n";
//  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);

//  #if DEBUG_PRINT_PACTBD
//        printf("INVOKER() on: %s\n", __PRETTY_FUNCTION__);
//        std::cout << "Invoker3\n";
//  std::cout << "ActiveGraph - finished\n";
//  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);

//  std::cout << "C++ Partial Template Specialization\n";
//  #endif

//} // InitialiseActiveGraph()

//// UNCOMMENT AUTOMATED























//  OLD commented below (before 2024-03-15)

//template <class Mesh>
//inline void ActiveGraph::Initialise(Mesh& mesh, const MeshExtrema& meshExtrema)
//{ // InitialiseActiveGraph()
//  // reference to the correct array in the extrema

//  // 2024-03-10: We will be writing portals
//  using IdPortalType = IdArrayType::WritePortalType;
////  using ReadPortalType = typename vtkm::cont::ArrayHandle<vtkm::Id,Storage>::ReadPortalType;


//  const char* message = "debug";
//  long workaround = 100;
//  std::cout << "--------------- extrema.DebugPrint(message, message, 100); --------------- \n";
////  meshExtrema.DebugPrint(message, message, workaround);
//  PrintHeader(meshExtrema.Peaks.GetNumberOfValues());
//  PrintIndices("MeshExtremaPeaks", meshExtrema.Peaks);
//  PrintIndices("MeshExtremaPits", meshExtrema.Pits);
//  std::cout << "A1\n";
//  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

//  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;

//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
////        std::cout << "i=" << i << MaskedIndex(extrema.Peaks[i]) << ", " <<  extrema.Peaks[i];
////      extrema.Get(i);
//      PrintIndices("MeshExtremaPeaksAFTER", extrema);
//  }
//  std::cout << "------------------------------- fine until here -------------------------------\n";


////  / *
//  // For every vertex, work out whether it is critical
//  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
//  // All vertices of outdegree 0 must be extrema
//  // Saddle points must be at least outdegree 2, so this is a correct test
//  // BUT it is possible to overestimate the degree of a non-extremum,
//  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph
//  std::cout << "A2\n";
//  // Neighbourhood mask (one bit set per connected component in neighbourhood
//  IdArrayType neighbourhoodMasks;
//  neighbourhoodMasks.Allocate(mesh.NumVertices);
//  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
//  outDegrees.Allocate(mesh.NumVertices);
//  std::cout << "A3\n";
//  // Initialize the nerighborhoodMasks and outDegrees arrays
//  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
//  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
//  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
//    this->IsJoinGraph);
//  std::cout << "A4\n";
//  this->Invoke(initNeighMasksAndOutDegWorklet,
//               sortIndexArray,
//               mesh,
//               neighbourhoodMasks, // output
//               outDegrees);        // output

//  // next, we compute where each vertex lands in the new array
//  // it needs to be one place offset, hence the +/- 1
//  // this should automatically parallelise
//  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
//  /*auto oneIfCritical = [](unsigned x) { return x!= 1 ? 1 : 0; };

//    // we need a temporary inverse index to change vertex IDs
//    IdArrayType inverseIndex;
//    inverseIndex.Allocate(mesh.NumVertices);
//    inverseIndex.WritePortal().Set(0,0);

//    std::partial_sum(
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(outDegrees.WritePortal()), oneIfCritical),
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(outDegrees.WritePortal())-1, oneIfCritical),
//            vtkm::cont::ArrayPortalToIteratorBegin(inverseIndex.WritePortal()) + 1);
//COMMENT    */
//  std::cout << "------------------------------- fine until here -------------------------------\n";
//  std::cout << "(Just initialised sort indices and outdegs up there...)\n";
//  std::cout << "------------------------------- fine until here -------------------------------\n";
//  std::cout << "A5\n";
//  IdArrayType inverseIndex;
//  inverseIndex.Allocate(mesh.NumVertices);
////  OneIfCritical oneIfCriticalFunctor;
////  auto oneIfCriticalArrayHandle =
////    vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfCritical>(outDegrees, oneIfCriticalFunctor);
////  vtkm::cont::Algorithm::ScanExclusive(oneIfCriticalArrayHandle, inverseIndex);

//  IdPortalType inversePortal = inverseIndex.WritePortal();
//  for(vtkm::Id i = 0; i < mesh.NumVertices; i++)
//  {
//      inversePortal.Set(i, i);
//  }


//  std::cout << "A6\n";
//  // now we can compute how many critical points we carry forward
////  vtkm::Id nCriticalPoints =
////    this->GetLastValue(inverseIndex) + oneIfCriticalFunctor(this->GetLastValue(outDegrees));


//  vtkm::Id nCriticalPoints = mesh.NumVertices;
////  if(up_or_down)
////  {

////  }


//  std::cout << "A7\n";
//  // we need to keep track of what the index of each vertex is in the active graph
//  // for most vertices, this should have the NO_SUCH_VERTEX flag set
//  AllocateVertexArrays(
//    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices
//  std::cout << "A8\n";
//  // our processing now depends on the degree of the vertex
//  // but basically, we want to set up the arrays for this vertex:
//  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
//  // GlobalIndex stores the index in the join tree for later access
//  IdArrayType activeIndices;
//  activeIndices.Allocate(mesh.NumVertices);
//  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
//    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
//  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);
//  std::cout << "A9\n";

//  auto readSortIdx =  sortIndexArray.ReadPortal();
//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
//      std::cout << readSortIdx.Get(i) << " ";
//  }
//  std::cout << "\n";

////  auto readNhbOffsets =  mesh.NeighborOffsets.ReadPortal();
//////  PrintIndices(mesh.NeighborOffsets);
////  for(unsigned i = 0; i < mesh.NumVertices; i++)
////  {
////      std::cout << readNhbOffsets.Get(i) << " ";
////  }


////  PrintIndices("FOR_INIT:sortIndexArray",    sortIndexArray);
//  PrintIndices("FOR_INIT:outDegrees",        outDegrees);
//  PrintIndices("FOR_INIT:inverseIndex",      inverseIndex);
//  PrintIndices("FOR_INIT:extrema",           extrema);
//  PrintIndices("FOR_INIT:activeIndices",     activeIndices);

//  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
//  this->Invoke(initActiveGraphVerticesWorklet,
//               sortIndexArray,
//               outDegrees,
//               inverseIndex,
//               extrema,
//               activeIndices,
//               this->GlobalIndex,
//               this->Outdegree,
//               this->Hyperarcs,
//               this->ActiveVertices);
//  std::cout << "A10\n";
//  PrintIndices("INIT:GlobalIndex",    this->GlobalIndex);
//  PrintIndices("INIT:Outdegree",      this->Outdegree);
//  PrintIndices("INIT:Hyperarcs",      this->Hyperarcs);
//  PrintIndices("INIT:ActiveVertices", this->ActiveVertices);
////  * /

////  // ------------------------- EXPERIMENTAL STUFF ------------------------ //
////  vtkm::Id nCriticalPoints = mesh.NumVertices;
////  AllocateVertexArrays(
////    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices

////  std::cout << "A8\n";
////  // our processing now depends on the degree of the vertex
////  // but basically, we want to set up the arrays for this vertex:
////  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
////  // GlobalIndex stores the index in the join tree for later access
////  IdArrayType activeIndices;
////  activeIndices.Allocate(mesh.NumVertices);
////  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
////    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
////  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);

////  // --------------------- END OF EXPERIMENTAL STUFF --------------------- //

//  IdPortalType avPortal = this->ActiveVertices.WritePortal();
//  IdPortalType haPortal = this->Hyperarcs.WritePortal();
//  IdPortalType odPortal = this->Outdegree.WritePortal();
//  IdPortalType giPortal = this->GlobalIndex.WritePortal();

//  auto exPortal = extrema.ReadPortal();

//  vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,1,2,1,1,3,2,1,1,0,1,0}; // bare
////vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,2,3,2,0,0}; // optimised

//  vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,24,22,22,24,24,22,24,22,22,24,24}; // bare
////vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,  22,   24,22,22,24,24,22,24,22,22,24,24}; // optimised

//  vtkm::Id testoutdeg2[] = {0,0,1,1,2,2,2,2,0,1,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // bare
////vtkm::Id testoutdeg2[] = {0,0,2,2,2,2,0,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // optimised

//  vtkm::Id testextrem2[] = {0,1,0,1,0,0,1,1,8,0,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // bare
////vtkm::Id testextrem2[] = {0,1,0,0,1,1,8,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // optimised

//  // hacking around:
//  if(up_or_down)
//  {
//      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//      {
//          avPortal.Set(i, i);
//    //      haPortal.Set(i, (vtkm::Id)exPortal.Get(i));
//          haPortal.Set(i, testextrem1[i]);
//    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
//          odPortal.Set(i, testoutdeg1[i]);
//          giPortal.Set(i, i);
//      }
//  }
//  else
//  {
//      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//      {
//          avPortal.Set(i, i);
//    //      haPortal.Set(i, (vtkm::Id)exPortal.Get(i));
//          haPortal.Set(i, testextrem2[i]);
//    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
//          odPortal.Set(i, testoutdeg2[i]);
//          giPortal.Set(i, i);
//      }
//  }

//  PrintIndices("VX:activeIndices",  activeIndices);
//  PrintIndices("VX:GlobalIndex",    this->GlobalIndex);
//  PrintIndices("VX:Outdegree",      this->Outdegree);
//  PrintIndices("VX:Hyperarcs",      this->Hyperarcs);
//  PrintIndices("VX:ActiveVertices", this->ActiveVertices);
//  // now we need to compute the FirstEdge array from the outDegrees
//  this->FirstEdge.Allocate(nCriticalPoints);
//  // STD Version of the prefix sum
//  //this->FirstEdge.WritePortal().Set(0, 0);
//  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
//  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
//  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
//  // VTKM Version of the prefix sum
//  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
//  // Compute the number of critical edges

//  vtkm::Id nCriticalEdges =
//    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

////  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;
//  if(up_or_down)
//  {
//      nCriticalEdges = 56;
//  }
//  else
//  {
////      nCriticalEdges = 53;
//      nCriticalEdges = 56;
//  }

//  AllocateEdgeArrays(nCriticalEdges);
//  std::cout << "A11\n";
//  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
//  std::cout << "vvv ALL ACTIVE EDGES INITIALIZED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES INITIALIZED! ^^^\n";
//  std::cout << "A11.1\n";
//  // try doing this manually:
//  std::cout << "initActiveEdgesWorkled params: Outdegree=\n";
////  << this->Outdegree << ", GlobalIndex=" << this->GlobalIndex << "\n";
////  PrintIndices("outDegrees?: ", outDegrees);
//  PrintIndices("OUTDEGREE?: ", this->Outdegree);
//  PrintIndices("GlobalIDX(sort)?: ", this->GlobalIndex);
//  std::cout << "nCriticalEdges?: " << nCriticalEdges << "\n";

////                              const MeshStructureType& meshStructure,
////                              const vtkm::Id& firstEdgeIndex,
////                              const vtkm::Id& sortIndex, // = GlobalIndex.Get(activeIndex)
////                              const InFieldPortalType& extrema,
////                              const InFieldPortalType& neighbourhoodMasks,
////                              const OutFieldPortalType& edgeNear,
////                              const OutFieldPortalType& edgeFar,
////                              const OutFieldPortalType& activeEdges) const

////  this->Invoke(initActiveEdgesWorklet,
////               /* _1 */ this->Outdegree,         // FieldIn outdegree            // const vtkm::Id& outdegree
////               /* InputIndex */          // const vtkm::Id activeIndex,
////               /* _2 */ mesh,                    // ExecObject meshStructure     // const vtkm::Id activeIndex,
////               /* _3 */ this->FirstEdge,         // FieldIn firstEdge
////               /* _4 */ this->GlobalIndex,       // FieldIn globalIndex
////               /* _5 */ extrema,                 // WholeArrayIn extrema
////               /* _6 */ neighbourhoodMasks,      // WholeArrayIn neighbourhoodMasks
////               /* _7 */ this->EdgeNear,          // WholeArrayOut edgeNear
////               /* _8 */ this->EdgeFar,           // WholeArrayOut edgeFar
////               /* _9 */ this->ActiveEdges);      // WholeArrayOut activeEdges


//  // currently using REAL IDs (trying to make ALL vertices active)
//  vtkm::Id low_arr1[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,15,15,16,17,18,18,18,19,19,20,21,23};
//  // ACTIVE IDS:
//                      //{0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,14,14,14,14,15,15,15,16,16,16,16,17};

//                    //  {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
//  vtkm::Id up_arr1[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,24,22,22,24,24,24,24,22,22,24,22,24};
////  {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
////  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
//  // unoptimised:
//  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
//  int num1 = 56; //50;

//  //      std::vector<vtkm::Id> indices =
//  //      std::vector<vtkm::Id> low_arr =
//  //      std::vector<vtkm::Id> up_arr =
//  // currently using REAL IDs (trying to make ALL vertices active)
//  vtkm::Id low_arr2[] = {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  // shortened 53:
////  vtkm::Id low_arr2[] = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  // peaks as pits: {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,15,15,15,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,18,18,18};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//          //{2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  vtkm::Id up_arr2[] = {0,1,0,0,0,0,1,1,1,1,0,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  // shortened 53
////  vtkm::Id up_arr2[] = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //{0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //shortened 53
//  //  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
//  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52, 53,54,55};
//  int num2 = 56;

//  int num=0;

//  std::vector<vtkm::Id> indices;
//  std::vector<vtkm::Id> low_arr;
//  std::vector<vtkm::Id> up_arr;

////  using IdPortalType = IdArrayType::WritePortalType;

//  IdPortalType aePortal = this->ActiveEdges.WritePortal();
//  IdPortalType enPortal = this->EdgeNear.WritePortal();
//  IdPortalType efPortal = this->EdgeFar.WritePortal();

//  if(up_or_down)
//  {
//      this->ActiveEdges.Allocate(num1);
//      this->EdgeNear.Allocate(num1);
//      this->EdgeFar.Allocate(num1);

//      //num = num1;
////      for(unsigned i = 0; i < num1; i++)
//      for(vtkm::Id i = 0; i < num1; i++)
//      {
//        indices.push_back(indices1[i]);
//        low_arr.push_back(low_arr1[i]);
//        up_arr.push_back(up_arr1[i]);

//        aePortal.Set(i, indices1[i]);
//        enPortal.Set(i, low_arr1[i]);
//        efPortal.Set(i, up_arr1[i]);
//      }
//  }
//  else
//  {
//      this->ActiveEdges.Allocate(num2);
//      this->EdgeNear.Allocate(num2);
//      this->EdgeFar.Allocate(num2);

////      for(unsigned i = 0; i < num2; i++)
//      for(vtkm::Id i = 0; i < num2; i++)
//      {
//        indices.push_back(indices2[i]);
//        low_arr.push_back(low_arr2[i]);
//        up_arr.push_back(up_arr2[i]);

//        aePortal.Set(i, indices2[i]);
//        enPortal.Set(i, low_arr2[i]);
//        efPortal.Set(i, up_arr2[i]);

//      }
//  }



////  std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
////  std::vector<vtkm::Id> low_arr = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////  std::vector<vtkm::Id> up_arr = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};


////  this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////  this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////  this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);


////  if(up_or_down)
////  {
////      std::cout << "UP\n";
////      std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
////      std::vector<vtkm::Id> low_arr = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////      std::vector<vtkm::Id> up_arr = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};

////      this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////      this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////      this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);
////  }
////  else
////  {
////      std::cout << "DOWN\n";
////      std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
////      std::vector<vtkm::Id> low_arr = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
////      std::vector<vtkm::Id> up_arr = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};

////      this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////      this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////      this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);
////  }



////  this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr, vtkm::CopyFlag::Off);
////  this->EdgeFar  = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);

////  for (vtkm::Id edge = 0; edge < 50; edge++)
////  { // per edge
//    //edgeNear.Set(edgeID, activeIndex);
////      this->EdgeNear.Set(edge, low_arr[edge]);
////      this->EdgeNear[edge] = low_arr[edge];
//    //edgeFar.Set(edgeID, MaskedIndex(extrema.Get(neigbourComponents[edge])) );
////      this->EdgeFar.Set(edge, up_arr[edge] );
////      this->EdgeFar[edge] = up_arr[edge];

//      // and save the edge itself
////      this->ActiveEdges.Set(edge, edge);
////      this->ActiveEdges[edge] = edge;
////  }
//  std::cout << "A11.1 FINISHED!\n";
//  if(up_or_down)
//  {
//    std::cout << "UP\n";
//  }
//  else
//  {
//    std::cout << "DOWN\n";
//  }

//  std::cout << "vvv ALL ACTIVE EDGES GATHERED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES GATHERED! ^^^\n";

//  std::cout << "A11.2\n";
//  // now we have to go through and set the far ends of the new edges using the
//  // inverse index array
//  std::cout << "A11.3\n";
//  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
//  std::cout << "A11.4\n";
//  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

//  DebugPrint("Active Graph Started", __FILE__, __LINE__);
//  std::cout << "A12\n";
//  // then we loop through the active vertices to convert their indices to active graph indices
//  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
//  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);
//  std::cout << "A13\n";
//  // finally, allocate and initialise the edgeSorter array
//  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
//  std::cout << "A14\n";
//  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);
//  std::cout << "ActiveGraph - finished\n";
//  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);
//} // InitialiseActiveGraph()








































































//// initialises the active graph
//template <class Mesh>
//inline void ActiveGraph::Initialise(Mesh& mesh, const MeshExtrema& meshExtrema)
//{ // InitialiseActiveGraph()
//  // reference to the correct array in the extrema

//  // 2024-03-10: We will be writing portals
//  using IdPortalType = IdArrayType::WritePortalType;
////  using ReadPortalType = typename vtkm::cont::ArrayHandle<vtkm::Id,Storage>::ReadPortalType;


//  const char* message = "debug";
//  long workaround = 100;
//  std::cout << "--------------- extrema.DebugPrint(message, message, 100); --------------- \n";
////  meshExtrema.DebugPrint(message, message, workaround);
//  PrintHeader(meshExtrema.Peaks.GetNumberOfValues());
//  PrintIndices("MeshExtremaPeaks", meshExtrema.Peaks);
//  PrintIndices("MeshExtremaPits", meshExtrema.Pits);
//  std::cout << "A1\n";
//  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

//  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;

//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
////        std::cout << "i=" << i << MaskedIndex(extrema.Peaks[i]) << ", " <<  extrema.Peaks[i];
////      extrema.Get(i);
//      PrintIndices("MeshExtremaPeaksAFTER", extrema);
//  }
//  std::cout << "------------------------------- fine until here -------------------------------\n";


////  / *
//  // For every vertex, work out whether it is critical
//  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
//  // All vertices of outdegree 0 must be extrema
//  // Saddle points must be at least outdegree 2, so this is a correct test
//  // BUT it is possible to overestimate the degree of a non-extremum,
//  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph
//  std::cout << "A2\n";
//  // Neighbourhood mask (one bit set per connected component in neighbourhood
//  IdArrayType neighbourhoodMasks;
//  neighbourhoodMasks.Allocate(mesh.NumVertices);
//  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
//  outDegrees.Allocate(mesh.NumVertices);
//  std::cout << "A3\n";
//  // Initialize the nerighborhoodMasks and outDegrees arrays
//  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
//  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
//  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
//    this->IsJoinGraph);
//  std::cout << "A4\n";
//  this->Invoke(initNeighMasksAndOutDegWorklet,
//               sortIndexArray,
//               mesh,
//               neighbourhoodMasks, // output
//               outDegrees);        // output

//  // next, we compute where each vertex lands in the new array
//  // it needs to be one place offset, hence the +/- 1
//  // this should automatically parallelise
//  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
//  /*auto oneIfCritical = [](unsigned x) { return x!= 1 ? 1 : 0; };

//    // we need a temporary inverse index to change vertex IDs
//    IdArrayType inverseIndex;
//    inverseIndex.Allocate(mesh.NumVertices);
//    inverseIndex.WritePortal().Set(0,0);

//    std::partial_sum(
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(outDegrees.WritePortal()), oneIfCritical),
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(outDegrees.WritePortal())-1, oneIfCritical),
//            vtkm::cont::ArrayPortalToIteratorBegin(inverseIndex.WritePortal()) + 1);
//COMMENT    */
//  std::cout << "------------------------------- fine until here -------------------------------\n";
//  std::cout << "(Just initialised sort indices and outdegs up there...)\n";
//  std::cout << "------------------------------- fine until here -------------------------------\n";
//  std::cout << "A5\n";
//  IdArrayType inverseIndex;
//  inverseIndex.Allocate(mesh.NumVertices);
////  OneIfCritical oneIfCriticalFunctor;
////  auto oneIfCriticalArrayHandle =
////    vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfCritical>(outDegrees, oneIfCriticalFunctor);
////  vtkm::cont::Algorithm::ScanExclusive(oneIfCriticalArrayHandle, inverseIndex);

//  IdPortalType inversePortal = inverseIndex.WritePortal();
//  for(vtkm::Id i = 0; i < mesh.NumVertices; i++)
//  {
//      inversePortal.Set(i, i);
//  }


//  std::cout << "A6\n";
//  // now we can compute how many critical points we carry forward
////  vtkm::Id nCriticalPoints =
////    this->GetLastValue(inverseIndex) + oneIfCriticalFunctor(this->GetLastValue(outDegrees));


//  vtkm::Id nCriticalPoints = mesh.NumVertices;
////  if(up_or_down)
////  {

////  }


//  std::cout << "A7\n";
//  // we need to keep track of what the index of each vertex is in the active graph
//  // for most vertices, this should have the NO_SUCH_VERTEX flag set
//  AllocateVertexArrays(
//    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices
//  std::cout << "A8\n";
//  // our processing now depends on the degree of the vertex
//  // but basically, we want to set up the arrays for this vertex:
//  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
//  // GlobalIndex stores the index in the join tree for later access
//  IdArrayType activeIndices;
//  activeIndices.Allocate(mesh.NumVertices);
//  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
//    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
//  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);
//  std::cout << "A9\n";

//  auto readSortIdx =  sortIndexArray.ReadPortal();
//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
//      std::cout << readSortIdx.Get(i) << " ";
//  }
//  std::cout << "\n";

////  auto readNhbOffsets =  mesh.NeighborOffsets.ReadPortal();
////  PrintIndices(mesh.NeighborOffsets);
//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
////      std::cout << readNhbOffsets.Get(i) << " ";
//  }


////  PrintIndices("FOR_INIT:sortIndexArray",    sortIndexArray);
//  PrintIndices("FOR_INIT:outDegrees",        outDegrees);
//  PrintIndices("FOR_INIT:inverseIndex",      inverseIndex);
//  PrintIndices("FOR_INIT:extrema",           extrema);
//  PrintIndices("FOR_INIT:activeIndices",     activeIndices);

//  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
//  this->Invoke(initActiveGraphVerticesWorklet,
//               sortIndexArray,
//               outDegrees,
//               inverseIndex,
//               extrema,
//               activeIndices,
//               this->GlobalIndex,
//               this->Outdegree,
//               this->Hyperarcs,
//               this->ActiveVertices);
//  std::cout << "A10\n";
//  PrintIndices("INIT:GlobalIndex",    this->GlobalIndex);
//  PrintIndices("INIT:Outdegree",      this->Outdegree);
//  PrintIndices("INIT:Hyperarcs",      this->Hyperarcs);
//  PrintIndices("INIT:ActiveVertices", this->ActiveVertices);
////  * /

////  // ------------------------- EXPERIMENTAL STUFF ------------------------ //
////  vtkm::Id nCriticalPoints = mesh.NumVertices;
////  AllocateVertexArrays(
////    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices

////  std::cout << "A8\n";
////  // our processing now depends on the degree of the vertex
////  // but basically, we want to set up the arrays for this vertex:
////  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
////  // GlobalIndex stores the index in the join tree for later access
////  IdArrayType activeIndices;
////  activeIndices.Allocate(mesh.NumVertices);
////  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
////    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
////  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);

////  // --------------------- END OF EXPERIMENTAL STUFF --------------------- //

//  IdPortalType avPortal = this->ActiveVertices.WritePortal();
//  IdPortalType haPortal = this->Hyperarcs.WritePortal();
//  IdPortalType odPortal = this->Outdegree.WritePortal();
//  IdPortalType giPortal = this->GlobalIndex.WritePortal();

//  auto exPortal = extrema.ReadPortal();

//  vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,1,2,1,1,3,2,1,1,0,1,0}; // bare
////vtkm::Id testoutdeg1[] = {2,2,3,3,2,4,2,4,4,5,4,4,4,0,2,3,2,0,0}; // optimised

//  vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,24,22,22,24,24,22,24,22,22,24,24}; // bare
////vtkm::Id testextrem1[] = {24,24,22,13,24,24,24,24,22,24,22,22,24,13,  22,   24,22,22,24,24,22,24,22,22,24,24}; // optimised

//  vtkm::Id testoutdeg2[] = {0,0,1,1,2,2,2,2,0,1,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // bare
////vtkm::Id testoutdeg2[] = {0,0,2,2,2,2,0,2,2,2,4,3,4,3,3,3,2,3,2,6,3,3}; // optimised

//  vtkm::Id testextrem2[] = {0,1,0,1,0,0,1,1,8,0,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // bare
////vtkm::Id testextrem2[] = {0,1,0,0,1,1,8,1,0,1,1,1,0,0,0,0,0,1,8,8,0,0}; // optimised

//  // hacking around:
//  if(up_or_down)
//  {
//      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//      {
//          avPortal.Set(i, i);
//    //      haPortal.Set(i, (vtkm::Id)exPortal.Get(i));
//          haPortal.Set(i, testextrem1[i]);
//    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
//          odPortal.Set(i, testoutdeg1[i]);
//          giPortal.Set(i, i);
//      }
//  }
//  else
//  {
//      for(vtkm::Id i = 0; i < nCriticalPoints; i++)
//      {
//          avPortal.Set(i, i);
//    //      haPortal.Set(i, (vtkm::Id)exPortal.Get(i));
//          haPortal.Set(i, testextrem2[i]);
//    //      haPortal.Set(i, extrema.ReadPortal().Get(i));
//          odPortal.Set(i, testoutdeg2[i]);
//          giPortal.Set(i, i);
//      }
//  }

//PrintIndices("VX:activeIndices",  activeIndices);
//  PrintIndices("VX:GlobalIndex",    this->GlobalIndex);
//  PrintIndices("VX:Outdegree",      this->Outdegree);
//  PrintIndices("VX:Hyperarcs",      this->Hyperarcs);
//  PrintIndices("VX:ActiveVertices", this->ActiveVertices);
//  // now we need to compute the FirstEdge array from the outDegrees
//  this->FirstEdge.Allocate(nCriticalPoints);
//  // STD Version of the prefix sum
//  //this->FirstEdge.WritePortal().Set(0, 0);
//  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
//  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
//  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
//  // VTKM Version of the prefix sum
//  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
//  // Compute the number of critical edges

//  vtkm::Id nCriticalEdges =
//    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

////  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;
//  if(up_or_down)
//  {
//      nCriticalEdges = 56;
//  }
//  else
//  {
////      nCriticalEdges = 53;
//      nCriticalEdges = 56;
//  }

//  AllocateEdgeArrays(nCriticalEdges);
//  std::cout << "A11\n";
//  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
//  std::cout << "vvv ALL ACTIVE EDGES INITIALIZED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES INITIALIZED! ^^^\n";
//  std::cout << "A11.1\n";
//  // try doing this manually:
//  std::cout << "initActiveEdgesWorkled params: Outdegree=\n";
////  << this->Outdegree << ", GlobalIndex=" << this->GlobalIndex << "\n";
////  PrintIndices("outDegrees?: ", outDegrees);
//  PrintIndices("OUTDEGREE?: ", this->Outdegree);
//  PrintIndices("GlobalIDX(sort)?: ", this->GlobalIndex);
//  std::cout << "nCriticalEdges?: " << nCriticalEdges << "\n";

////                              const MeshStructureType& meshStructure,
////                              const vtkm::Id& firstEdgeIndex,
////                              const vtkm::Id& sortIndex, // = GlobalIndex.Get(activeIndex)
////                              const InFieldPortalType& extrema,
////                              const InFieldPortalType& neighbourhoodMasks,
////                              const OutFieldPortalType& edgeNear,
////                              const OutFieldPortalType& edgeFar,
////                              const OutFieldPortalType& activeEdges) const

////  this->Invoke(initActiveEdgesWorklet,
////               /* _1 */ this->Outdegree,         // FieldIn outdegree            // const vtkm::Id& outdegree
////               /* InputIndex */          // const vtkm::Id activeIndex,
////               /* _2 */ mesh,                    // ExecObject meshStructure     // const vtkm::Id activeIndex,
////               /* _3 */ this->FirstEdge,         // FieldIn firstEdge
////               /* _4 */ this->GlobalIndex,       // FieldIn globalIndex
////               /* _5 */ extrema,                 // WholeArrayIn extrema
////               /* _6 */ neighbourhoodMasks,      // WholeArrayIn neighbourhoodMasks
////               /* _7 */ this->EdgeNear,          // WholeArrayOut edgeNear
////               /* _8 */ this->EdgeFar,           // WholeArrayOut edgeFar
////               /* _9 */ this->ActiveEdges);      // WholeArrayOut activeEdges


//  // currently using REAL IDs (trying to make ALL vertices active)
//  vtkm::Id low_arr1[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,15,15,16,17,18,18,18,19,19,20,21,23};
//  // ACTIVE IDS:
//                      //{0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,14,14,14,14,15,15,15,16,16,16,16,17};

//                    //  {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
//  vtkm::Id up_arr1[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,24,22,22,24,24,24,24,22,22,24,22,24};
////  {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
////  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
//  // unoptimised:
//  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
//  int num1 = 56; //50;

//  //      std::vector<vtkm::Id> indices =
//  //      std::vector<vtkm::Id> low_arr =
//  //      std::vector<vtkm::Id> up_arr =
//  // currently using REAL IDs (trying to make ALL vertices active)
//  vtkm::Id low_arr2[] = {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  // shortened 53:
////  vtkm::Id low_arr2[] = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  // peaks as pits: {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,15,15,15,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,18,18,18};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//          //{2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  vtkm::Id up_arr2[] = {0,1,0,0,0,0,1,1,1,1,0,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  // shortened 53
////  vtkm::Id up_arr2[] = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //{0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //shortened 53
//  //  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
//  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52, 53,54,55};
//  int num2 = 56;

//  int num=0;

//  std::vector<vtkm::Id> indices;
//  std::vector<vtkm::Id> low_arr;
//  std::vector<vtkm::Id> up_arr;

////  using IdPortalType = IdArrayType::WritePortalType;

//  IdPortalType aePortal = this->ActiveEdges.WritePortal();
//  IdPortalType enPortal = this->EdgeNear.WritePortal();
//  IdPortalType efPortal = this->EdgeFar.WritePortal();

//  if(up_or_down)
//  {
//      this->ActiveEdges.Allocate(num1);
//      this->EdgeNear.Allocate(num1);
//      this->EdgeFar.Allocate(num1);

//      //num = num1;
////      for(unsigned i = 0; i < num1; i++)
//      for(vtkm::Id i = 0; i < num1; i++)
//      {
//        indices.push_back(indices1[i]);
//        low_arr.push_back(low_arr1[i]);
//        up_arr.push_back(up_arr1[i]);

//        aePortal.Set(i, indices1[i]);
//        enPortal.Set(i, low_arr1[i]);
//        efPortal.Set(i, up_arr1[i]);
//      }
//  }
//  else
//  {
//      this->ActiveEdges.Allocate(num2);
//      this->EdgeNear.Allocate(num2);
//      this->EdgeFar.Allocate(num2);

////      for(unsigned i = 0; i < num2; i++)
//      for(vtkm::Id i = 0; i < num2; i++)
//      {
//        indices.push_back(indices2[i]);
//        low_arr.push_back(low_arr2[i]);
//        up_arr.push_back(up_arr2[i]);

//        aePortal.Set(i, indices2[i]);
//        enPortal.Set(i, low_arr2[i]);
//        efPortal.Set(i, up_arr2[i]);

//      }
//  }



////  std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
////  std::vector<vtkm::Id> low_arr = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////  std::vector<vtkm::Id> up_arr = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};


////  this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////  this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////  this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);


////  if(up_or_down)
////  {
////      std::cout << "UP\n";
////      std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
////      std::vector<vtkm::Id> low_arr = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////      std::vector<vtkm::Id> up_arr = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};

////      this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////      this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////      this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);
////  }
////  else
////  {
////      std::cout << "DOWN\n";
////      std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
////      std::vector<vtkm::Id> low_arr = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
////      std::vector<vtkm::Id> up_arr = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};

////      this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////      this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////      this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);
////  }



////  this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr, vtkm::CopyFlag::Off);
////  this->EdgeFar  = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);

////  for (vtkm::Id edge = 0; edge < 50; edge++)
////  { // per edge
//    //edgeNear.Set(edgeID, activeIndex);
////      this->EdgeNear.Set(edge, low_arr[edge]);
////      this->EdgeNear[edge] = low_arr[edge];
//    //edgeFar.Set(edgeID, MaskedIndex(extrema.Get(neigbourComponents[edge])) );
////      this->EdgeFar.Set(edge, up_arr[edge] );
////      this->EdgeFar[edge] = up_arr[edge];

//      // and save the edge itself
////      this->ActiveEdges.Set(edge, edge);
////      this->ActiveEdges[edge] = edge;
////  }
//  std::cout << "A11.1 FINISHED!\n";
//  if(up_or_down)
//  {
//    std::cout << "UP\n";
//  }
//  else
//  {
//    std::cout << "DOWN\n";
//  }

//  std::cout << "vvv ALL ACTIVE EDGES GATHERED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES GATHERED! ^^^\n";

//  std::cout << "A11.2\n";
//  // now we have to go through and set the far ends of the new edges using the
//  // inverse index array
//  std::cout << "A11.3\n";
//  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
//  std::cout << "A11.4\n";
//  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

//  DebugPrint("Active Graph Started", __FILE__, __LINE__);
//  std::cout << "A12\n";
//  // then we loop through the active vertices to convert their indices to active graph indices
//  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
//  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);
//  std::cout << "A13\n";
//  // finally, allocate and initialise the edgeSorter array
//  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
//  std::cout << "A14\n";
//  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);
//  std::cout << "ActiveGraph - finished\n";
//  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);
//} // InitialiseActiveGraph()
































//// initialises the active graph
//template <class Mesh>
//inline void ActiveGraph::Initialise(Mesh& mesh, const MeshExtrema& meshExtrema)
//{ // InitialiseActiveGraph()
//  // reference to the correct array in the extrema
//  const char* message = "debug";
//  long workaround = 100;
//  std::cout << "--------------- extrema.DebugPrint(message, message, 100); --------------- \n";
////  meshExtrema.DebugPrint(message, message, workaround);
//  PrintHeader(meshExtrema.Peaks.GetNumberOfValues());
//  PrintIndices("MeshExtremaPeaks", meshExtrema.Peaks);
//  PrintIndices("MeshExtremaPits", meshExtrema.Pits);
//  std::cout << "A1\n";
//  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

//  for(unsigned i = 0; i < mesh.NumVertices; i++)
//  {
////        std::cout << "i=" << i << MaskedIndex(extrema.Peaks[i]) << ", " <<  extrema.Peaks[i];
////      extrema.Get(i);
//      PrintIndices("MeshExtremaPeaksAFTER", extrema);
//  }
//  std::cout << "------------------------------- fine until here -------------------------------\n";

//  // For every vertex, work out whether it is critical
//  // We do so by computing outdegree in the mesh & suppressing the vertex if outdegree is 1
//  // All vertices of outdegree 0 must be extrema
//  // Saddle points must be at least outdegree 2, so this is a correct test
//  // BUT it is possible to overestimate the degree of a non-extremum,
//  // The test is therefore necessary but not sufficient, and extra vertices are put in the active graph
//  std::cout << "A2\n";
//  // Neighbourhood mask (one bit set per connected component in neighbourhood
//  IdArrayType neighbourhoodMasks;
//  neighbourhoodMasks.Allocate(mesh.NumVertices);
//  IdArrayType outDegrees; // TODO Should we change this to an unsigned type
//  outDegrees.Allocate(mesh.NumVertices);
//  std::cout << "A3\n";
//  // Initialize the nerighborhoodMasks and outDegrees arrays
//  mesh.SetPrepareForExecutionBehavior(this->IsJoinGraph);
//  vtkm::cont::ArrayHandleIndex sortIndexArray(mesh.NumVertices);
//  active_graph_inc_ns::InitializeNeighbourhoodMasksAndOutDegrees initNeighMasksAndOutDegWorklet(
//    this->IsJoinGraph);
//  std::cout << "A4\n";
//  this->Invoke(initNeighMasksAndOutDegWorklet,
//               sortIndexArray,
//               mesh,
//               neighbourhoodMasks, // output
//               outDegrees);        // output

//  // next, we compute where each vertex lands in the new array
//  // it needs to be one place offset, hence the +/- 1
//  // this should automatically parallelise
//  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
//  /*auto oneIfCritical = [](unsigned x) { return x!= 1 ? 1 : 0; };

//    // we need a temporary inverse index to change vertex IDs
//    IdArrayType inverseIndex;
//    inverseIndex.Allocate(mesh.NumVertices);
//    inverseIndex.WritePortal().Set(0,0);

//    std::partial_sum(
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(outDegrees.WritePortal()), oneIfCritical),
//            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(outDegrees.WritePortal())-1, oneIfCritical),
//            vtkm::cont::ArrayPortalToIteratorBegin(inverseIndex.WritePortal()) + 1);
//    */
//  std::cout << "A5\n";
//  IdArrayType inverseIndex;
//  OneIfCritical oneIfCriticalFunctor;
//  auto oneIfCriticalArrayHandle =
//    vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfCritical>(outDegrees, oneIfCriticalFunctor);
//  vtkm::cont::Algorithm::ScanExclusive(oneIfCriticalArrayHandle, inverseIndex);
//  std::cout << "A6\n";
//  // now we can compute how many critical points we carry forward
//  vtkm::Id nCriticalPoints =
//    this->GetLastValue(inverseIndex) + oneIfCriticalFunctor(this->GetLastValue(outDegrees));
//  std::cout << "A7\n";
//  // we need to keep track of what the index of each vertex is in the active graph
//  // for most vertices, this should have the NO_SUCH_VERTEX flag set
//  AllocateVertexArrays(
//    nCriticalPoints); // allocates outdegree, GlobalIndex, Hyperarcs, ActiveVertices
//  std::cout << "A8\n";
//  // our processing now depends on the degree of the vertex
//  // but basically, we want to set up the arrays for this vertex:
//  // activeIndex gets the next available ID in the active graph (was called nearIndex before)
//  // GlobalIndex stores the index in the join tree for later access
//  IdArrayType activeIndices;
//  activeIndices.Allocate(mesh.NumVertices);
//  vtkm::cont::ArrayHandleConstant<vtkm::Id> noSuchElementArray(
//    static_cast<vtkm::Id>(NO_SUCH_ELEMENT), mesh.NumVertices);
//  vtkm::cont::Algorithm::Copy(noSuchElementArray, activeIndices);
//  std::cout << "A9\n";
//  active_graph_inc_ns::InitializeActiveGraphVertices initActiveGraphVerticesWorklet;
//  this->Invoke(initActiveGraphVerticesWorklet,
//               sortIndexArray,
//               outDegrees,
//               inverseIndex,
//               extrema,
//               activeIndices,
//               this->GlobalIndex,
//               this->Outdegree,
//               this->Hyperarcs,
//               this->ActiveVertices);
//  std::cout << "A10\n";
//  // now we need to compute the FirstEdge array from the outDegrees
//  this->FirstEdge.Allocate(nCriticalPoints);
//  // STD Version of the prefix sum
//  //this->FirstEdge.WritePortal().Set(0, 0);
//  //std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(this->Outdegree.GetPortalControl()),
//  //                 vtkm::cont::ArrayPortalToIteratorEnd(this->Outdegree.GetPortalControl()) - 1,
//  //                 vtkm::cont::ArrayPortalToIteratorBegin(this->firstEdge.GetPortalControl()) + 1);
//  // VTKM Version of the prefix sum
//  vtkm::cont::Algorithm::ScanExclusive(this->Outdegree, this->FirstEdge);
//  // Compute the number of critical edges

//  vtkm::Id nCriticalEdges =
//    this->GetLastValue(this->FirstEdge) + this->GetLastValue(this->Outdegree);

//  bool up_or_down = this->IsJoinGraph ? true : false; //meshExtrema.Peaks : meshExtrema.Pits;
//  if(up_or_down)
//  {
//      nCriticalEdges = 56;
//  }
//  else
//  {
//      nCriticalEdges = 53;
//  }

//  AllocateEdgeArrays(nCriticalEdges);
//  std::cout << "A11\n";
//  active_graph_inc_ns::InitializeActiveEdges<Mesh> initActiveEdgesWorklet;
//  std::cout << "vvv ALL ACTIVE EDGES INITIALIZED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES INITIALIZED! ^^^\n";
//  std::cout << "A11.1\n";
//  // try doing this manually:
//  std::cout << "initActiveEdgesWorkled params: Outdegree=\n";
////  << this->Outdegree << ", GlobalIndex=" << this->GlobalIndex << "\n";
//  PrintIndices("outDegrees?: ", outDegrees);
//  PrintIndices("OUTDEGREE?: ", this->Outdegree);
//  PrintIndices("GlobalIDX(sort)?: ", this->GlobalIndex);
//  std::cout << "nCriticalEdges?: " << nCriticalEdges << "\n";

////                              const MeshStructureType& meshStructure,
////                              const vtkm::Id& firstEdgeIndex,
////                              const vtkm::Id& sortIndex, // = GlobalIndex.Get(activeIndex)
////                              const InFieldPortalType& extrema,
////                              const InFieldPortalType& neighbourhoodMasks,
////                              const OutFieldPortalType& edgeNear,
////                              const OutFieldPortalType& edgeFar,
////                              const OutFieldPortalType& activeEdges) const

////  this->Invoke(initActiveEdgesWorklet,
////               /* _1 */ this->Outdegree,         // FieldIn outdegree            // const vtkm::Id& outdegree
////               /* InputIndex */          // const vtkm::Id activeIndex,
////               /* _2 */ mesh,                    // ExecObject meshStructure     // const vtkm::Id activeIndex,
////               /* _3 */ this->FirstEdge,         // FieldIn firstEdge
////               /* _4 */ this->GlobalIndex,       // FieldIn globalIndex
////               /* _5 */ extrema,                 // WholeArrayIn extrema
////               /* _6 */ neighbourhoodMasks,      // WholeArrayIn neighbourhoodMasks
////               /* _7 */ this->EdgeNear,          // WholeArrayOut edgeNear
////               /* _8 */ this->EdgeFar,           // WholeArrayOut edgeFar
////               /* _9 */ this->ActiveEdges);      // WholeArrayOut activeEdges



//  vtkm::Id low_arr1[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,14,14,14,14,15,15,15,16,16,16,16,17};
////  {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
//  vtkm::Id up_arr1[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,24,22,22,24,24,24,24,22,22,24,22,24};
////  {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
////  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
//  // unoptimised:
//  vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
//  int num1 = 56; //50;

//  //      std::vector<vtkm::Id> indices =
//  //      std::vector<vtkm::Id> low_arr =
//  //      std::vector<vtkm::Id> up_arr =
//  vtkm::Id low_arr2[] = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  // peaks as pits: {2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,15,15,15,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,18,18,18};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//          //{2,3,4,4,5,5,6,6,7,7,9,10,10,11,11,12,12,13,13,13,13,14,14,14,15,15,15,15,16,16,16,17,17,17,18,18,18,19,19,20,20,20,21,21,22,22,22,22,22,22,23,23,23,24,24,24};
//  //{2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
//  vtkm::Id up_arr2[] = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  //{0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
//  vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
//  int num2 = 53;

//  int num=0;

//  std::vector<vtkm::Id> indices;
//  std::vector<vtkm::Id> low_arr;
//  std::vector<vtkm::Id> up_arr;

//  using IdPortalType = IdArrayType::WritePortalType;

//  IdPortalType aePortal = this->ActiveEdges.WritePortal();
//  IdPortalType enPortal = this->EdgeNear.WritePortal();
//  IdPortalType efPortal = this->EdgeFar.WritePortal();

//  if(up_or_down)
//  {
//      this->ActiveEdges.Allocate(num1);
//      this->EdgeNear.Allocate(num1);
//      this->EdgeFar.Allocate(num1);

//      //num = num1;
////      for(unsigned i = 0; i < num1; i++)
//      for(vtkm::Id i = 0; i < num1; i++)
//      {
//        indices.push_back(indices1[i]);
//        low_arr.push_back(low_arr1[i]);
//        up_arr.push_back(up_arr1[i]);

//        aePortal.Set(i, indices1[i]);
//        enPortal.Set(i, low_arr1[i]);
//        efPortal.Set(i, up_arr1[i]);
//      }
//  }
//  else
//  {
//      this->ActiveEdges.Allocate(num2);
//      this->EdgeNear.Allocate(num2);
//      this->EdgeFar.Allocate(num2);

////      for(unsigned i = 0; i < num2; i++)
//      for(vtkm::Id i = 0; i < num2; i++)
//      {
//        indices.push_back(indices2[i]);
//        low_arr.push_back(low_arr2[i]);
//        up_arr.push_back(up_arr2[i]);

//        aePortal.Set(i, indices2[i]);
//        enPortal.Set(i, low_arr2[i]);
//        efPortal.Set(i, up_arr2[i]);

//      }
//  }



////  std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
////  std::vector<vtkm::Id> low_arr = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////  std::vector<vtkm::Id> up_arr = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};


////  this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////  this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////  this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);


////  if(up_or_down)
////  {
////      std::cout << "UP\n";
////      std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
////      std::vector<vtkm::Id> low_arr = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
////      std::vector<vtkm::Id> up_arr = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};

////      this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////      this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////      this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);
////  }
////  else
////  {
////      std::cout << "DOWN\n";
////      std::vector<vtkm::Id> indices = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
////      std::vector<vtkm::Id> low_arr = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
////      std::vector<vtkm::Id> up_arr = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};

////      this->ActiveEdges = vtkm::cont::make_ArrayHandle(indices, vtkm::CopyFlag::Off);
////      this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr , vtkm::CopyFlag::Off);
////      this->EdgeFar = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);
////  }



////  this->EdgeNear = vtkm::cont::make_ArrayHandle(low_arr, vtkm::CopyFlag::Off);
////  this->EdgeFar  = vtkm::cont::make_ArrayHandle(up_arr, vtkm::CopyFlag::Off);

////  for (vtkm::Id edge = 0; edge < 50; edge++)
////  { // per edge
//    //edgeNear.Set(edgeID, activeIndex);
////      this->EdgeNear.Set(edge, low_arr[edge]);
////      this->EdgeNear[edge] = low_arr[edge];
//    //edgeFar.Set(edgeID, MaskedIndex(extrema.Get(neigbourComponents[edge])) );
////      this->EdgeFar.Set(edge, up_arr[edge] );
////      this->EdgeFar[edge] = up_arr[edge];

//      // and save the edge itself
////      this->ActiveEdges.Set(edge, edge);
////      this->ActiveEdges[edge] = edge;
////  }
//  std::cout << "A11.1 FINISHED!\n";
//  if(up_or_down)
//  {
//    std::cout << "UP\n";
//  }
//  else
//  {
//    std::cout << "DOWN\n";
//  }

//  std::cout << "vvv ALL ACTIVE EDGES GATHERED! vvv\n";
//  PrintIndices("ACTIVE EDGES", this->ActiveEdges);
//  PrintIndices("NEAR", this->EdgeNear);
//  PrintIndices("FAR", this->EdgeFar);
//  std::cout << "^^^ ALL ACTIVE EDGES GATHERED! ^^^\n";

//  std::cout << "A11.2\n";
//  // now we have to go through and set the far ends of the new edges using the
//  // inverse index array
//  std::cout << "A11.3\n";
//  active_graph_inc_ns::InitializeEdgeFarFromActiveIndices initEdgeFarWorklet;
//  std::cout << "A11.4\n";
//  this->Invoke(initEdgeFarWorklet, this->EdgeFar, extrema, activeIndices);

//  DebugPrint("Active Graph Started", __FILE__, __LINE__);
//  std::cout << "A12\n";
//  // then we loop through the active vertices to convert their indices to active graph indices
//  active_graph_inc_ns::InitializeHyperarcsFromActiveIndices initHyperarcsWorklet;
//  this->Invoke(initHyperarcsWorklet, this->Hyperarcs, activeIndices);
//  std::cout << "A13\n";
//  // finally, allocate and initialise the edgeSorter array
//  this->EdgeSorter.Allocate(this->ActiveEdges.GetNumberOfValues());
//  std::cout << "A14\n";
//  vtkm::cont::Algorithm::Copy(this->ActiveEdges, this->EdgeSorter);
//  std::cout << "ActiveGraph - finished\n";
//  //DebugPrint("Active Graph Initialised", __FILE__, __LINE__);
//} // InitialiseActiveGraph()


















// routine that computes the merge tree from the active graph
// was previously Compute()
inline void ActiveGraph::MakeMergeTree(MergeTree& tree, MeshExtrema& meshExtrema)
{ // MakeMergeTree()
  DebugPrint("Active Graph Computation Starting", __FILE__, __LINE__);
  /// DEBUG PRINT std::cout << "ActiveGraph::MakeMergeTree - 1\n";
  // loop until we run out of active edges
  vtkm::Id maxNumIterations = this->EdgeSorter.GetNumberOfValues();
  this->NumIterations = 0;
  while (true)
  { // main loop
    // choose the subset of edges for the governing saddles
    TransferSaddleStarts();

    // test whether there are any left (if not, we're on the trunk)
    if (this->EdgeSorter.GetNumberOfValues() <= 0)
    {
      break;
    }
    // test whether we are in a bad infinite loop due to bad input data. Usually
    // this is not an issue for the merge tree (only for the contour tree), but
    // we check just to make absolutely sure we won't get stuck in an infinite loop
    if (this->NumIterations >= maxNumIterations)
    {
      throw new vtkm::cont::ErrorInternal("Bad iteration. Merge tree unable to process all edges.");
    }

    // find & label the extrema with their governing saddles
    FindGoverningSaddles();

    // label the regular points
    TransferRegularPoints();

    // compact the active set of vertices & edges
    CompactActiveVertices();
    CompactActiveEdges();

    // rebuild the chains
    BuildChains();

    // increment the iteration count
    this->NumIterations++;
  } // main loop

  // final pass to label the trunk vertices
  BuildTrunk();

  // transfer results to merge tree
  FindSuperAndHyperNodes(tree);
  SetSuperArcs(tree);
  SetHyperArcs(tree);
  SetArcs(tree, meshExtrema);

  // we can now release many of the arrays to free up space
  ReleaseTemporaryArrays();

  /// DEBUG PRINT std::cout << "ActiveGraph::MakeMergeTree Merge Tree Computed - 5\n";
  DebugPrint("Merge Tree Computed", __FILE__, __LINE__);
} // MakeMergeTree()


// suppresses non-saddles for the governing saddles pass
inline void ActiveGraph::TransferSaddleStarts()
{ // TransferSaddleStarts()
  // update all of the edges so that the far end resets to the result of the ascent in the previous step

  active_graph_inc_ns::TransferSaddleStartsResetEdgeFar transferSaddleResetWorklet;
  this->Invoke(transferSaddleResetWorklet, this->ActiveEdges, this->Hyperarcs, this->EdgeFar);

  // in parallel, we need to create a vector to count the first edge for each vertex
  IdArrayType newOutdegree;
  newOutdegree.Allocate(this->ActiveVertices.GetNumberOfValues());

  // this will be a stream compaction later, but for now we'll do it the serial way
  active_graph_inc_ns::TransferSaddleStartsSetNewOutdegreeForSaddles transferOutDegree;
  this->Invoke(transferOutDegree,
               this->ActiveVertices,
               this->FirstEdge,
               this->Outdegree,
               this->ActiveEdges,
               this->Hyperarcs,
               this->EdgeFar,
               newOutdegree);

  // now do a parallel prefix sum using the offset partial sum trick.
  IdArrayType newFirstEdge;
  newFirstEdge.Allocate(this->ActiveVertices.GetNumberOfValues());
  // STD version of the prefix sum
  // newFirstEdge.WritePortal().Set(0, 0);
  // std::partial_sum(vtkm::cont::ArrayPortalToIteratorBegin(newOutdegree.WritePortal()),
  //                 vtkm::cont::ArrayPortalToIteratorEnd(newOutdegree.WritePortal()) - 1,
  //                 vtkm::cont::ArrayPortalToIteratorBegin(newFirstEdge.WritePortal()) + 1);
  // VTK:M version of the prefix sum
  vtkm::cont::Algorithm::ScanExclusive(newOutdegree, newFirstEdge);

  vtkm::Id nEdgesToSort = this->GetLastValue(newFirstEdge) + this->GetLastValue(newOutdegree);

  // now we write only the active saddle edges to the sorting array
  this->EdgeSorter
    .ReleaseResources(); // TODO is there a single way to resize an array handle without calling ReleaseResources followed by Allocate
  this->EdgeSorter.Allocate(nEdgesToSort);

  // this will be a stream compaction later, but for now we'll do it the serial way
  active_graph_inc_ns::TransferSaddleStartsUpdateEdgeSorter updateEdgeSorterWorklet;
  this->Invoke(updateEdgeSorterWorklet,
               this->ActiveVertices,
               this->ActiveEdges,
               this->FirstEdge,
               newFirstEdge,
               newOutdegree,
               this->EdgeSorter);

  DebugPrint("Saddle Starts Transferred", __FILE__, __LINE__);
} // TransferSaddleStarts()


// sorts saddle starts to find governing saddles
inline void ActiveGraph::FindGoverningSaddles()
{ // FindGoverningSaddles()
  // sort with the comparator
  vtkm::cont::Algorithm::Sort(
    this->EdgeSorter,
    active_graph_inc_ns::EdgePeakComparator(this->EdgeFar, this->EdgeNear, this->IsJoinGraph));

  // DebugPrint("After Sorting", __FILE__, __LINE__);

  // now loop through the edges to find the governing saddles
  active_graph_inc_ns::FindGoverningSaddlesWorklet findGovSaddlesWorklet;
  vtkm::cont::ArrayHandleIndex edgeIndexArray(this->EdgeSorter.GetNumberOfValues());

  this->Invoke(findGovSaddlesWorklet,
               edgeIndexArray,
               this->EdgeSorter,
               this->EdgeFar,
               this->EdgeNear,
               this->Hyperarcs,
               this->Outdegree);

  DebugPrint("Governing Saddles Set", __FILE__, __LINE__);
} // FindGoverningSaddles()


// marks now regular points for removal
inline void ActiveGraph::TransferRegularPoints()
{ // TransferRegularPointsWorklet
  // we need to label the regular points that have been identified
  active_graph_inc_ns::TransferRegularPointsWorklet transRegPtWorklet(this->IsJoinGraph);
  this->Invoke(transRegPtWorklet, this->ActiveVertices, this->Hyperarcs, this->Outdegree);

  DebugPrint("Regular Points Should Now Be Labelled", __FILE__, __LINE__);
} // TransferRegularPointsWorklet()


// compacts the active vertex list
inline void ActiveGraph::CompactActiveVertices()
{ // CompactActiveVertices()
  using PermuteIndexType = vtkm::cont::ArrayHandlePermutation<IdArrayType, IdArrayType>;

  // create a temporary array the same size
  vtkm::cont::ArrayHandle<vtkm::Id> newActiveVertices;

  // Use only the current this->ActiveVertices this->Outdegree to match size on CopyIf
  vtkm::cont::ArrayHandle<vtkm::Id> outdegreeLookup;
  vtkm::cont::Algorithm::Copy(PermuteIndexType(this->ActiveVertices, this->Outdegree),
                              outdegreeLookup);

  // compact the this->ActiveVertices array to keep only the ones of interest
  vtkm::cont::Algorithm::CopyIf(this->ActiveVertices, outdegreeLookup, newActiveVertices);

  this->ActiveVertices.ReleaseResources();
  vtkm::cont::Algorithm::Copy(newActiveVertices, this->ActiveVertices);

  DebugPrint("Active Vertex List Compacted", __FILE__, __LINE__);
} // CompactActiveVertices()


// compacts the active edge list
inline void ActiveGraph::CompactActiveEdges()
{ // CompactActiveEdges()
  // grab the size of the array for easier reference
  vtkm::Id nActiveVertices = this->ActiveVertices.GetNumberOfValues();

  // first, we have to work out the first edge for each active vertex
  // we start with a temporary new outdegree
  IdArrayType newOutdegree;
  newOutdegree.Allocate(nActiveVertices);

  // Run workflet to compute newOutdegree for each vertex
  active_graph_inc_ns::CompactActiveEdgesComputeNewVertexOutdegree computeNewOutdegreeWorklet;
  this->Invoke(computeNewOutdegreeWorklet,
               this->ActiveVertices, // (input)
               this->ActiveEdges,    // (input)
               this->EdgeFar,        // (input)
               this->FirstEdge,      // (input)
               this->Outdegree,      // (input)
               this->Hyperarcs,      // (input/output)
               newOutdegree          // (output)
  );

  // now we do a reduction to compute the offsets of each vertex
  vtkm::cont::ArrayHandle<vtkm::Id> newPosition;
  // newPosition.Allocate(nActiveVertices);   // Not necessary. ScanExclusive takes care of this.
  vtkm::cont::Algorithm::ScanExclusive(newOutdegree, newPosition);

  vtkm::Id nNewEdges = vtkm::cont::ArrayGetValue(nActiveVertices - 1, newPosition) +
    vtkm::cont::ArrayGetValue(nActiveVertices - 1, newOutdegree);

  // create a temporary vector for copying
  IdArrayType newActiveEdges;
  newActiveEdges.Allocate(nNewEdges);
  // overwriting Hyperarcs in parallel is safe, as the worst than can happen is
  // that another valid ascent is found; for comparison and validation purposes
  // however it makes sense to have a `canoical' computation. To achieve this
  // canonical computation, we need to write into a new array during computation
  // ensuring that we always use the same information. The following is left in
  // commented out for future debugging and validation

  //DebugPrint("Active Edges Counted", __FILE__, __LINE__);

  // now copy the relevant edges into the active edge array
  active_graph_inc_ns::CompactActiveEdgesTransferActiveEdges transferActiveEdgesWorklet;
  this->Invoke(transferActiveEdgesWorklet,
               this->ActiveVertices,
               newPosition,       // (input)
               newOutdegree,      // (input)
               this->ActiveEdges, // (input)
               newActiveEdges,    // (output)
               this->EdgeFar,     // (input/output)
               this->FirstEdge,   // (input/output)
               this->Outdegree,   // (input/output)
               this->Hyperarcs    // (input/output)
  );

  // resize the original array and recopy
  //vtkm::cont::Algorithm::::Copy(newActiveEdges, this-ActiveEdges);
  this->ActiveEdges.ReleaseResources();
  // vtkm ArrayHandles are smart, so we can just swap it in without having to copy
  this->ActiveEdges = newActiveEdges;

  // for canonical computation: swap in newly computed hyperarc array
  //      this->Hyperarcs.swap(newHyperarcs);

  DebugPrint("Active Edges Now Compacted", __FILE__, __LINE__);
} // CompactActiveEdges()



// builds the chains for the new active vertices
inline void ActiveGraph::BuildChains()
{ // BuildChains()
  // 1. compute the number of log steps required in this pass
  vtkm::Id numLogSteps = 1;
  for (vtkm::Id shifter = this->ActiveVertices.GetNumberOfValues(); shifter != 0; shifter >>= 1)
    numLogSteps++;

  // 2.   Use path compression / step doubling to collect vertices along chains
  //              until every vertex has been assigned to *an* extremum
  for (vtkm::Id logStep = 0; logStep < numLogSteps; logStep++)
  { // per log step
    active_graph_inc_ns::BuildChainsWorklet buildChainsWorklet;
    this->Invoke(buildChainsWorklet, this->ActiveVertices, this->Hyperarcs);
  } // per log step
  DebugPrint("Chains Built", __FILE__, __LINE__);
} // BuildChains()


// sets all remaining active vertices
inline void ActiveGraph::BuildTrunk()
{ //BuildTrunk
  // all remaining vertices belong to the trunk
  active_graph_inc_ns::BuildTrunkWorklet buildTrunkWorklet;
  this->Invoke(buildTrunkWorklet, this->ActiveVertices, this->Hyperarcs);

  DebugPrint("Trunk Built", __FILE__, __LINE__);
} //BuildTrunk


// finds all super and hyper nodes, numbers them & sets up arrays for lookup
inline void ActiveGraph::FindSuperAndHyperNodes(MergeTree& tree)
{ // FindSuperAndHyperNodes()
  // allocate memory for nodes
  this->HyperID.ReleaseResources();
  this->HyperID.Allocate(this->GlobalIndex.GetNumberOfValues());

  // compute new node positions
  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
  /*auto oneIfSupernode = [](vtkm::Id v) { return IsSupernode(v) ? 1 : 0; };
    IdArrayType newSupernodePosition;
    newSupernodePosition.Allocate(this->Hyperarcs.GetNumberOfValues());
    newSupernodePosition.WritePortal().Set(0, 0);

    std::partial_sum(
            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(this->Hyperarcs.WritePortal()), oneIfSupernode),
            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(this->Hyperarcs.WritePortal()) - 1, oneIfSupernode),
            vtkm::cont::ArrayPortalToIteratorBegin(newSupernodePosition.GetPortalControl()) + 1);*/

  IdArrayType newSupernodePosition;
  OneIfSupernode oneIfSupernodeFunctor;
  auto oneIfSupernodeArrayHandle = vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfSupernode>(
    this->Hyperarcs, oneIfSupernodeFunctor);
  vtkm::cont::Algorithm::ScanExclusive(oneIfSupernodeArrayHandle, newSupernodePosition);

  this->NumSupernodes = this->GetLastValue(newSupernodePosition) +
    oneIfSupernodeFunctor(this->GetLastValue(this->Hyperarcs));

  tree.Supernodes.ReleaseResources();
  tree.Supernodes.Allocate(this->NumSupernodes);

  // The following commented code block is variant ported directly from PPP2 using std::partial_sum. This has been replaced here with vtkm's ScanExclusive.
  /*
    auto oneIfHypernode = [](vtkm::Id v) { return IsHypernode(v) ? 1 : 0; };
    IdArrayType newHypernodePosition;
    newHypernodePosition.Allocate(this->Hyperarcs.GetNumberOfValues());
    newHypernodePosition.WritePortal().Set(0, 0);
    std::partial_sum(
            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorBegin(this->Hyperarcs.WritePortal()), oneIfHypernode),
            boost::make_transform_iterator(vtkm::cont::ArrayPortalToIteratorEnd(this->Hyperarcs.WritePortal()) - 1, oneIfHypernode),
             vtkm::cont::ArrayPortalToIteratorBegin(newHypernodePosition.GetPortalControl()) + 1);
    */
  IdArrayType newHypernodePosition;
  OneIfHypernode oneIfHypernodeFunctor;
  auto oneIfHypernodeArrayHandle = vtkm::cont::ArrayHandleTransform<IdArrayType, OneIfHypernode>(
    this->Hyperarcs, oneIfHypernodeFunctor);
  vtkm::cont::Algorithm::ScanExclusive(oneIfHypernodeArrayHandle, newHypernodePosition);

  this->NumHypernodes = this->GetLastValue(newHypernodePosition) +
    oneIfHypernodeFunctor(this->GetLastValue(this->Hyperarcs));

  tree.Hypernodes.ReleaseResources();
  tree.Hypernodes.Allocate(this->GlobalIndex.GetNumberOfValues());

  // perform stream compression
  active_graph_inc_ns::FindSuperAndHyperNodesWorklet findSuperAndHyperNodesWorklet;
  vtkm::cont::ArrayHandleIndex graphVertexIndex(this->GlobalIndex.GetNumberOfValues());
  this->Invoke(findSuperAndHyperNodesWorklet,
               graphVertexIndex,
               this->Hyperarcs,
               newHypernodePosition,
               newSupernodePosition,
               this->HyperID,
               tree.Hypernodes,
               tree.Supernodes);

  DebugPrint("Super/Hypernodes Found", __FILE__, __LINE__);
  tree.DebugPrint("Super/Hypernodes Found", __FILE__, __LINE__);
} // FindSuperAndHyperNodes()


// uses active graph to set superarcs & hyperparents in merge tree
inline void ActiveGraph::SetSuperArcs(MergeTree& tree)
{ // SetSuperArcs()
  using PermutedIdArrayType = vtkm::cont::ArrayHandlePermutation<IdArrayType, IdArrayType>;

  //      1.      set the hyperparents
  // allocate space for the hyperparents
  tree.Hyperparents.ReleaseResources();
  tree.Hyperparents.Allocate(this->NumSupernodes);

  // execute the worklet to set the hyperparents
  active_graph_inc_ns::SetSuperArcsSetTreeHyperparents setTreeHyperparentsWorklet;
  this->Invoke(setTreeHyperparentsWorklet, tree.Supernodes, this->Hyperarcs, tree.Hyperparents);

  tree.DebugPrint("Hyperparents Set", __FILE__, __LINE__);
  //      a.      And the super ID array needs setting up
  this->SuperID.ReleaseResources();
  vtkm::cont::Algorithm::Copy(
    vtkm::cont::make_ArrayHandleConstant(NO_SUCH_ELEMENT, this->GlobalIndex.GetNumberOfValues()),
    this->SuperID);
  vtkm::cont::ArrayHandleIndex supernodeIndex(this->NumSupernodes);
  PermutedIdArrayType permutedSuperID(tree.Supernodes, this->SuperID);
  vtkm::cont::Algorithm::Copy(supernodeIndex, permutedSuperID);

  //      2.      Sort the supernodes into segments according to hyperparent
  //              See comparator for details
  vtkm::cont::Algorithm::Sort(tree.Supernodes,
                              active_graph_inc_ns::HyperArcSuperNodeComparator(
                                tree.Hyperparents, this->SuperID, tree.IsJoinTree));

  //      3.      Now update the other arrays to match
  IdArrayType hyperParentsTemp;
  hyperParentsTemp.Allocate(this->NumSupernodes);
  auto permutedTreeHyperparents = vtkm::cont::make_ArrayHandlePermutation(
    vtkm::cont::make_ArrayHandlePermutation(tree.Supernodes, this->SuperID), tree.Hyperparents);

  vtkm::cont::Algorithm::Copy(permutedTreeHyperparents, hyperParentsTemp);
  vtkm::cont::Algorithm::Copy(hyperParentsTemp, tree.Hyperparents);
  hyperParentsTemp.ReleaseResources();
  //      a.      And the super ID array needs setting up // TODO Check if we really need this?
  vtkm::cont::Algorithm::Copy(supernodeIndex, permutedSuperID);

  DebugPrint("Supernodes Sorted", __FILE__, __LINE__);
  tree.DebugPrint("Supernodes Sorted", __FILE__, __LINE__);

  //      4.      Allocate memory for superarcs
  tree.Superarcs.ReleaseResources();
  tree.Superarcs.Allocate(this->NumSupernodes);
  tree.FirstSuperchild.ReleaseResources();
  tree.FirstSuperchild.Allocate(this->NumHypernodes);

  //      5.      Each supernode points to its neighbour in the list, except at the end of segments
  // execute the worklet to set the tree.Hyperparents and tree.FirstSuperchild
  active_graph_inc_ns::SetSuperArcsSetTreeSuperarcs setTreeSuperarcsWorklet;
  this->Invoke(setTreeSuperarcsWorklet,
               tree.Supernodes,     // (input)
               this->Hyperarcs,     // (input)
               tree.Hyperparents,   // (input)
               this->SuperID,       // (input)
               this->HyperID,       // (input)
               tree.Superarcs,      // (output)
               tree.FirstSuperchild // (output)
  );

  // 6.   Now we can reset the supernodes to mesh IDs
  PermutedIdArrayType permuteGlobalIndex(tree.Supernodes, this->GlobalIndex);
  vtkm::cont::Algorithm::Copy(permuteGlobalIndex, tree.Supernodes);

  // 7.   and the hyperparent to point to a hyperarc rather than a graph index
  PermutedIdArrayType permuteHyperID(tree.Hyperparents, this->HyperID);
  vtkm::cont::Algorithm::Copy(permuteHyperID, tree.Hyperparents);

  tree.DebugPrint("Superarcs Set", __FILE__, __LINE__);
} // SetSuperArcs()


// uses active graph to set hypernodes in merge tree
inline void ActiveGraph::SetHyperArcs(MergeTree& tree)
{ // SetHyperArcs()
  //      1.      Allocate memory for hypertree
  tree.Hypernodes.Allocate(
    this->NumHypernodes, // Has been allocated previously.
    vtkm::CopyFlag::On); // The values are needed but the size may be too large.
  tree.Hyperarcs.ReleaseResources();
  tree.Hyperarcs.Allocate(this->NumHypernodes); // Has not been allocated yet

  //      2.      Use the superIDs already set to fill in the Hyperarcs array
  active_graph_inc_ns::SetHyperArcsWorklet setHyperArcsWorklet;
  this->Invoke(
    setHyperArcsWorklet, tree.Hypernodes, tree.Hyperarcs, this->Hyperarcs, this->SuperID);

  // Debug output
  DebugPrint("Hyperarcs Set", __FILE__, __LINE__);
  tree.DebugPrint("Hyperarcs Set", __FILE__, __LINE__);
} // SetHyperArcs()


// uses active graph to set arcs in merge tree
inline void ActiveGraph::SetArcs(MergeTree& tree, MeshExtrema& meshExtrema)
{ // SetArcs()
  using PermuteIndexType = vtkm::cont::ArrayHandlePermutation<IdArrayType, IdArrayType>;

  // reference to the correct array in the extrema
  const IdArrayType& extrema = this->IsJoinGraph ? meshExtrema.Peaks : meshExtrema.Pits;

  // 1.   Set the arcs for the super/hypernodes based on where they prune to
  active_graph_inc_ns::SetArcsSetSuperAndHypernodeArcs setSuperAndHypernodeArcsWorklet;
  this->Invoke(setSuperAndHypernodeArcsWorklet,
               this->GlobalIndex,
               this->Hyperarcs,
               this->HyperID,
               tree.Arcs,
               tree.Superparents);

  DebugPrint("Sliding Arcs Set", __FILE__, __LINE__);
  tree.DebugPrint("Sliding Arcs Set", __FILE__, __LINE__);

  // 2.   Loop through all vertices to slide down Hyperarcs
  active_graph_inc_ns::SetArcsSlideVertices slideVerticesWorklet(
    this->IsJoinGraph, this->NumSupernodes, this->NumHypernodes);
  this->Invoke(slideVerticesWorklet,
               tree.Arcs,            // (input)
               extrema,              // (input)  i.e,. meshExtrema.Peaks or meshExtrema.Pits
               tree.FirstSuperchild, // (input)
               tree.Supernodes,      // (input)
               tree.Superparents);   // (input/output)

  tree.DebugPrint("Sliding Finished", __FILE__, __LINE__);

  // 3.   Now set the superparents correctly for the supernodes
  PermuteIndexType permuteTreeSuperparents(tree.Supernodes, tree.Superparents);
  vtkm::cont::ArrayHandleIndex supernodesIndex(this->NumSupernodes);
  vtkm::cont::Algorithm::Copy(supernodesIndex, permuteTreeSuperparents);

  tree.DebugPrint("Superparents Set", __FILE__, __LINE__);

  // 4.   Finally, sort all of the vertices onto their superarcs
  IdArrayType nodes;
  vtkm::cont::ArrayHandleIndex nodesIndex(tree.Arcs.GetNumberOfValues());
  vtkm::cont::Algorithm::Copy(nodesIndex, nodes);

  //  5.  Sort the nodes into segments according to superparent
  //      See comparator for details
  vtkm::cont::Algorithm::Sort(
    nodes, active_graph_inc_ns::SuperArcNodeComparator(tree.Superparents, tree.IsJoinTree));

  //  6. Connect the nodes to each other
  active_graph_inc_ns::SetArcsConnectNodes connectNodesWorklet;
  this->Invoke(connectNodesWorklet,
               tree.Arcs,         // (input/output)
               nodes,             // (input)
               tree.Superparents, // (input)
               tree.Superarcs,    // (input)
               tree.Supernodes);  // (input)

  tree.DebugPrint("Arcs Set", __FILE__, __LINE__);
} // SetArcs()


// Allocate the vertex array
inline void ActiveGraph::AllocateVertexArrays(vtkm::Id nElems)
{
  this->GlobalIndex.Allocate(nElems);
  this->Outdegree.Allocate(nElems);
  this->Hyperarcs.Allocate(nElems);
  this->ActiveVertices.Allocate(nElems);
}


// Allocate the edge array
inline void ActiveGraph::AllocateEdgeArrays(vtkm::Id nElems)
{
  this->ActiveEdges.Allocate(nElems);
  this->EdgeNear.Allocate(nElems);
  this->EdgeFar.Allocate(nElems);
}


// releases temporary arrays
inline void ActiveGraph::ReleaseTemporaryArrays()
{
  this->GlobalIndex.ReleaseResources();
  this->FirstEdge.ReleaseResources();
  this->Outdegree.ReleaseResources();
  this->EdgeNear.ReleaseResources();
  this->EdgeFar.ReleaseResources();
  this->ActiveEdges.ReleaseResources();
  this->ActiveVertices.ReleaseResources();
  this->EdgeSorter.ReleaseResources();
  this->Hyperarcs.ReleaseResources();
  this->HyperID.ReleaseResources();
  this->SuperID.ReleaseResources();
}


// prints the contents of the active graph in a standard format
inline void ActiveGraph::DebugPrint(const char* message, const char* fileName, long lineNum)
{ // DebugPrint()
#if DEBUG_PRINT_PACTBD // WE WANT DEBUG PRINT ACTIVE GRAPH
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::setw(30) << std::left << fileName << ":" << std::right << std::setw(4)
            << lineNum << std::endl;
  std::cout << std::left << std::string(message) << std::endl;
  std::cout << "Active Graph Contains:                                " << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

  std::cout << "Is Join Graph? " << (this->IsJoinGraph ? "T" : "F") << std::endl;
  std::cout << "NumIterations    " << this->NumIterations << std::endl;
  std::cout << "nSupernodes    " << this->NumSupernodes << std::endl;
  std::cout << "nHypernodes    " << this->NumHypernodes << std::endl;

  // Full Vertex Arrays
  std::cout << "Full Vertex Arrays - Size:  " << this->GlobalIndex.GetNumberOfValues() << std::endl;
  PrintHeader(this->GlobalIndex.GetNumberOfValues());
  PrintIndices("Global Index", this->GlobalIndex);
  PrintIndices("First Edge", this->FirstEdge);
  PrintIndices("Outdegree", this->Outdegree);
  PrintIndices("Hyperarc ID", this->Hyperarcs);
  PrintIndices("Hypernode ID", this->HyperID);
  PrintIndices("Supernode ID", this->SuperID);
  std::cout << std::endl;

  // Active Vertex Arrays
  IdArrayType activeIndices;
  PermuteArray<vtkm::Id>(this->GlobalIndex, this->ActiveVertices, activeIndices);
  IdArrayType activeFirst;
  PermuteArray<vtkm::Id>(this->FirstEdge, this->ActiveVertices, activeFirst);
  IdArrayType activeOutdegree;
  PermuteArray<vtkm::Id>(this->Outdegree, this->ActiveVertices, activeOutdegree);
  IdArrayType activeHyperarcs;
  PermuteArray<vtkm::Id>(this->Hyperarcs, this->ActiveVertices, activeHyperarcs);
  std::cout << "Active Vertex Arrays - Size: " << this->ActiveVertices.GetNumberOfValues()
            << std::endl;
  PrintHeader(this->ActiveVertices.GetNumberOfValues());
  PrintIndices("Active Vertices", this->ActiveVertices);
  PrintIndices("Active Indices", activeIndices);
  PrintIndices("Active First Edge", activeFirst);
  PrintIndices("Active Outdegree", activeOutdegree);
  PrintIndices("Active Hyperarc ID", activeHyperarcs);
  std::cout << std::endl;

  // Full Edge Arrays
  IdArrayType farIndices;
  PermuteArray<vtkm::Id>(this->GlobalIndex, this->EdgeFar, farIndices);
  IdArrayType nearIndices;
  PermuteArray<vtkm::Id>(this->GlobalIndex, this->EdgeNear, nearIndices);
  std::cout << "Full Edge Arrays - Size:     " << this->EdgeNear.GetNumberOfValues() << std::endl;
  PrintHeader(this->EdgeFar.GetNumberOfValues());
  PrintIndices("Near", this->EdgeNear);
  PrintIndices("Far", this->EdgeFar);
  PrintIndices("Near Index", nearIndices);
  PrintIndices("Far Index", farIndices);
  std::cout << std::endl;

  // Active Edge Arrays
  IdArrayType activeFarIndices;
  PermuteArray<vtkm::Id>(this->EdgeFar, this->ActiveEdges, activeFarIndices);
  IdArrayType activeNearIndices;
  PermuteArray<vtkm::Id>(this->EdgeNear, this->ActiveEdges, activeNearIndices);
  std::cout << "Active Edge Arrays - Size:   " << this->ActiveEdges.GetNumberOfValues()
            << std::endl;
  PrintHeader(this->ActiveEdges.GetNumberOfValues());
  PrintIndices("Active Edges", this->ActiveEdges);
  PrintIndices("Edge Near Index", activeNearIndices);
  PrintIndices("Edge Far Index", activeFarIndices);
  std::cout << std::endl;

  // Edge Sorter Array
  IdArrayType sortedFarIndices;
  PermuteArray<vtkm::Id>(this->EdgeFar, this->EdgeSorter, sortedFarIndices);
  IdArrayType sortedNearIndices;
  PermuteArray<vtkm::Id>(this->EdgeNear, this->EdgeSorter, sortedNearIndices);
  std::cout << "Edge Sorter - Size:          " << this->EdgeSorter.GetNumberOfValues() << std::endl;
  PrintHeader(this->EdgeSorter.GetNumberOfValues());
  PrintIndices("Edge Sorter", this->EdgeSorter);
  PrintIndices("Sorted Near Index", sortedNearIndices);
  PrintIndices("Sorted Far Index", sortedFarIndices);
  std::cout << std::endl;

  std::cout << "---------------------------" << std::endl;
  std::cout << std::endl;
#else
  // Prevent unused parameter warning
  (void)message;
  (void)fileName;
  (void)lineNum;
#endif
} // DebugPrint()



} // namespace contourtree_augmented
} // worklet
} // vtkm

#endif
