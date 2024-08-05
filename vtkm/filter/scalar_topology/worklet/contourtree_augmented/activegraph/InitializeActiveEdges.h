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

#ifndef vtk_m_worklet_contourtree_augmented_active_graph_initialize_active_edges_h
#define vtk_m_worklet_contourtree_augmented_active_graph_initialize_active_edges_h

#include <vtkm/exec/arg/BasicArg.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <chrono>
#include <thread>

#include <bitset>

namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{
namespace active_graph_inc
{


// Worklet for computing the sort indices from the sort order
template <class MeshClassType>
class InitializeActiveEdges : public vtkm::worklet::WorkletMapField
{
public:
  /// Additional basic execution argument tags
  //struct _10 : vtkm::exec::arg::BasicArg<10> {  };
  //struct _11 : vtkm::exec::arg::BasicArg<11> {  };

  typedef void ControlSignature(
    FieldIn outdegree,               // (input) outdegree
    ExecObject meshStructure,        // (input) execution object with the mesh structure
    FieldIn firstEdge,               // (input)
    FieldIn globalIndex,             // (input) ActiveGraph.GlobalIndex
    WholeArrayIn extrema,            // (input)
    WholeArrayIn neighbourhoodMasks, // (input)
    WholeArrayOut edgeNear,          // (output) edgeNear
    WholeArrayOut edgeFar,           // (output) edgeFar
    WholeArrayOut activeEdges);      // (output) activeEdges
  typedef void ExecutionSignature(_1, InputIndex, _2, _3, _4, _5, _6, _7, _8, _9);

  using InputDomain = _1;

  // Default Constructor
  VTKM_EXEC_CONT
  InitializeActiveEdges() {}

  template <typename MeshStructureType, typename InFieldPortalType, typename OutFieldPortalType>
  VTKM_EXEC void operator()(const vtkm::Id& outdegree,
                            const vtkm::Id activeIndex,
                            const MeshStructureType& meshStructure,
                            const vtkm::Id& firstEdgeIndex,
                            const vtkm::Id& sortIndex, // = GlobalIndex.Get(activeIndex)
                            const InFieldPortalType& extrema,
                            const InFieldPortalType& neighbourhoodMasks,
                            const OutFieldPortalType& edgeNear,
                            const OutFieldPortalType& edgeFar,
                            const OutFieldPortalType& activeEdges) const
  {

  // DEBUG_PRINT  std::cout << "InitializeActiveEdges 1\n";

    bool GET_MESH = false;

////    // UNCOMMENT FOR WRITING FOR FREUDENTHAL MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//    // Writing to a file
//    std::ofstream file("MeshOutput-hh96-E.txt");

////    for(vtkm::Id i = 0; i < 125; i++)
////    for(vtkm::Id i = 0; i < 29791; i++)
////    for(vtkm::Id i = 0; i < 56770560; i++)
////    for(vtkm::Id i = 0; i < 13824; i++)
//    for(vtkm::Id i = 0; i < 884736; i++)
//    {
//      if (file.is_open()) file << i << " "; //",";
//      else std::cerr << "Unable to open file for writing.\n";
////      std::cout << i << " ";
//      // loop through all POSSIBLE nbors ...
//      for (vtkm::Id nbrNo = 0; nbrNo < meshStructure.GetMaxNumberOfNeighbours(); ++nbrNo)
//      {
//        if ( (neighbourhoodMasks.Get(i) & (static_cast<vtkm::Id>(1) << nbrNo)) )
//        { // if the current vertex actually HAS a neighbour at number bween [0, max] ...
//          // ... and also implicitly if that nbor is higher (nbhood mask returns black/white)
//          if(GET_MESH)
//          {
//              if (file.is_open()) file << meshStructure.GetNeighbourIndex(i, nbrNo) << " ";
//              else std::cerr << "Unable to open file for writing.\n";
////            std::cout << meshStructure.GetNeighbourIndex(i, nbrNo) << " ";
//          }
//        }
//      }
//      if (file.is_open()) file << "\n";
//      else std::cerr << "Unable to open file for writing.\n";
////      std::cout << "\n";
//    }

//    file.close(); // Make sure to close the file when done






























//    bool GET_MESH = true;

//    if(GET_MESH) std::cout << sortIndex << " ";

    if (outdegree != 0)
    {
      // DEBUG_PRINT  if(!GET_MESH) std::cout << "InitializeActiveEdges 2 - outdegree: " << outdegree << "\n";
      // temporary array for storing edges
      vtkm::Id neigbourComponents[MeshClassType::MAX_OUTDEGREE];
      int currNbrNo = 0;
      if(!GET_MESH)
      {
          // DEBUG_PRINT std::cout << "InitializeActiveEdges 3 - meshStructure.GetMaxNumberOfNeighbours(): "
          // DEBUG_PRINT       << meshStructure.GetMaxNumberOfNeighbours()
          // DEBUG_PRINT        << " -vs- " << MeshClassType::MAX_OUTDEGREE << "\n";
      }

      // loop through all POSSIBLE nbors ...
      // DEBUG_PRINT for (vtkm::Id nbrNo = 0; nbrNo < meshStructure.GetMaxNumberOfNeighbours(); ++nbrNo)
      for (vtkm::Id nbrNo = 0; nbrNo < 63; ++nbrNo)
      {
//        std::cout << "v" << neighbourhoodMasks.GetNumberOfValues() << "\n";
        // DEBUG_PRINT if(!GET_MESH) std::cout << nbrNo << " " << "s" << sortIndex << " ";
//        std::this_thread::sleep_for(std::chrono::nanoseconds(100));
        // added (currNbrNo+1 <= outdegree)

        if ( (currNbrNo+1 <= outdegree) && (neighbourhoodMasks.Get(sortIndex) & (static_cast<vtkm::Id>(1) << nbrNo)) )
        { // if the current vertex actually HAS a neighbour at number bween [0, max] ...
          // ... and also implicitly if that nbor is higher (nbhood mask returns black/white)
          std::bitset<64> x(neighbourhoodMasks.Get(sortIndex));
          std::bitset<64> y((static_cast<vtkm::Id>(1) << nbrNo));

          // DEBUG_PRINT if(!GET_MESH) std::cout << "x=" << x << " y=" << y << " ";

          // DEBUG_PRINT if(!GET_MESH) std::cout << "[" << currNbrNo << "]->";
          // DEBUG_PRINT if(!GET_MESH) std::cout << neigbourComponents[currNbrNo] << "->";
//          neigbourComponents[currNbrNo++] = meshStructure.GetNeighbourIndex(sortIndex, nbrNo);
          neigbourComponents[currNbrNo] = meshStructure.GetNeighbourIndex(sortIndex, nbrNo);

//          if(GET_MESH)
//          {
//            std::cout << neigbourComponents[currNbrNo] << " ";
//          }

//          if(GET_MESH)
//          {
//            std::cout << "\nx " << sortIndex+1 << " - ";
//            std::cout << meshStructure.GetNeighbourIndex(sortIndex+1, 0) << " ";
//            std::cout << meshStructure.GetNeighbourIndex(sortIndex+1, 1) << " ";
//          }


          currNbrNo++;

          // DEBUG_PRINT if(!GET_MESH) std::cout << "[" << currNbrNo << "]->";
          // DEBUG_PRINT if(!GET_MESH) std::cout << neigbourComponents[currNbrNo] << "->";
          // DEBUG_PRINT if(!GET_MESH) std::cout << "[" << currNbrNo << "] ";

        }
//        std::cout << "If success\n";
      }
      // DEBUG_PRINT if(!GET_MESH) std::cout << "InitializeActiveEdges 4\n";

      if(outdegree > 64) // == meshStructure.GetMaxNumberOfNeighbours())
      {
        std::this_thread::sleep_for(std::chrono::seconds(1));
      }


      // arcs stores the ID from the join tree - i.e. the chain extremum
      //              we cannot store the correct ID yet, because it may not have been assigned yet
      //              in serial, we could hack around this by processing the vertices in a given order
      //              but in parallel we can't, so we have two stages:
      //                      in this stage, we store the join tree ID (after suppressing flags)
      //                      in a later stage, we convert it to an active graph ID
      // firstEdge / outdegree / edgeNear / edgeFar are straightforward
      // as with earlier versions, the parallel equivalent will need to use stream compression
      // but the serial version can be expressed more simply.

//      int low_arr[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
//      int up_arr[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
      // removed duplicates:
      // low1     {0,0,1,1,2,2,3,3,4,5,5,6,7,7,7,8,8,9,9,10,10,10,11,12,12,14,14,15,16};//
      // up1      {22,24,13,24,22,24,13,24,24,22,24,24,13,22,24,13,22,22,24,13,22,24,22,22,24,22,24,24,22};//
      // indices1 {0,1,2,3,6,4,9,7,10,13,12,16,20,18,19,23,22,26,28,32,33,31,35,40,39,44,43,45,48}; //
      vtkm::Id low_arr1[] = {0,0,1,1,2,2,2,3,3,3,4,4,5,5,5,5,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,14,14,15,15,15,16,16};
      vtkm::Id up_arr1[] = {22,24,13,24,24,24,22,24,24,13,24,24,24,22,22,24,24,24,22,24,13,24,22,13,22,22,22,22,24,24,24,24,13,22,22,22,22,22,22,24,22,24,24,24,22,24,24,24,22,22};
      vtkm::Id indices1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
      int num1 = 50; //29

      //      std::vector<vtkm::Id> indices =
      //      std::vector<vtkm::Id> low_arr =
      //      std::vector<vtkm::Id> up_arr =
      // removed duplicates:
      // low2     {2,3,4,5,7,7,8,9,10,10,11,12,12,13,14,15,15,16,17,17,18,18,19,19,19,20,21,21};//
      // up2      {0,0,1,1,1,8,0,1,1,8,1,0,1,0,0,0,1,0,0,1,0,8,0,1,8,0,0,1};//
      // indices2 {0,2,4,6,8,9,10,12,14,16,18,21,22,25,28,31,32,34,38,36,40,39,43,42,41,47,50,51}; //
      // original output:
      vtkm::Id low_arr2[] = {2,2,3,3,4,4,5,5,7,7,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,16,17,17,17,18,18,19,19,19,19,19,19,20,20,20,21,21,21};
      vtkm::Id up_arr2[] = {0,0,0,0,1,1,1,1,1,8,0,0,1,1,1,1,8,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,8,0,8,1,0,0,0,8,0,0,0,0,1,0};
      vtkm::Id indices2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
      int num2 = 53; //28

      std::vector<vtkm::Id> indices;
      std::vector<vtkm::Id> low_arr;
      std::vector<vtkm::Id> up_arr;

      if(MaskedIndex(extrema.Get(0)))
      {
          //num = num1;
          for(unsigned i = 0; i < num1; i++)
          {
//            indices.push_back(indices1[i]);
//            low_arr.push_back(low_arr1[i]);
//            up_arr.push_back(up_arr1[i]);
          }
      }
      else
      {
          for(unsigned i = 0; i < num2; i++)
          {
//            indices.push_back(indices2[i]);
//            low_arr.push_back(low_arr2[i]);
//            up_arr.push_back(up_arr2[i]);
          }
      }

      for (vtkm::Id edge = 0; edge < outdegree; edge++)
      { // per edge
        // compute the edge index in the edge arrays

        vtkm::Id edgeID = firstEdgeIndex + edge;
        // DEBUG_PRINT if(!GET_MESH) std::cout << "\nEXTREMA: " << MaskedIndex(extrema.Get(0)) << "\n";
        // DEBUG_PRINT if(!GET_MESH) std::cout << "\nActive: " << activeIndex << "\n";
        // DEBUG_PRINT if(!GET_MESH) std::cout << edge << "vs edgeID: " << edgeID << "\n";
        // DEBUG_PRINT if(!GET_MESH) std::cout << "\nFT: " << sortIndex << "->" << neigbourComponents[edge] << "\n";
        // DEBUG_PRINT if(!GET_MESH) std::cout << "LH: " << neigbourComponents[edge]  << "->" << MaskedIndex(extrema.Get(neigbourComponents[edge])) << "\n";

        // now set the low and high ends
        edgeNear.Set(edgeID, activeIndex);
//        edgeNear.Set(edgeID, sortIndex);
//        edgeNear.Set(edgeID, low_arr[edgeID]);
//        edgeNear.Set(edgeID, 1);
        edgeFar.Set(edgeID, MaskedIndex(extrema.Get(neigbourComponents[edge])) );
//        edgeFar.Set(edgeID, up_arr[edgeID] );
//        edgeFar.Set(edgeID, 1 );

        // and save the edge itself
        activeEdges.Set(edgeID, edgeID);
      } // per edge

    }

    // DEBUG_PRINT if(!GET_MESH) std::cout << "\nInitializeActiveEdges 5 - End\n";

    /* DEBUG_PRINT if(!GET_MESH)
    {
      std::cout << "\n";
    }
    */

    // This operator implements the following loop from the serial code
    /*for (indexType activeIndex = 0; activeIndex < nCriticalPoints; ++activeIndex)
      {
        indexType outdegree = outdegree[activeIndex];
        if (outdegree != 0)
          {
            indexType sortIndex = globalIndex[activeIndex];

            // temporary array for storing edges
            indexType neigbourComponents[Mesh::MAX_OUTDEGREE];
            int currNbrNo = 0;
            for (vtkm::Int32 nbrNo = 0; nbrNo < mesh.GetMaxNumberOfNeighbours(); ++nbrNo)
              if (neighbourhoodMasks[sortIndex] & 1 << nbrNo)
                {
                   neigbourComponents[currNbrNo++] = mesh.GetNeighbourIndex(sortIndex, nbrNo);
                }


            // arcs stores the ID from the join tree - i.e. the chain extremum
            //              we cannot store the correct ID yet, because it may not have been assigned yet
            //              in serial, we could hack around this by processing the vertices in a given order
            //              but in parallel we can't, so we have two stages:
            //                      in this stage, we store the join tree ID (after suppressing flags)
            //                      in a later stage, we convert it to an active graph ID
            // firstEdge / outdegree / edgeNear / edgeFar are straightforward
            // as with earlier versions, the parallel equivalent will need to use stream compression
            // but the serial version can be expressed more simply.

            for (indexType edge = 0; edge < outdegree; edge++)
              { // per edge
                // compute the edge index in the edge arrays
                indexType edgeID = firstEdge[activeIndex] + edge;

                // now set the low and high ends
                edgeNear[edgeID] = activeIndex;
                edgeFar[edgeID] = MaskedIndex(extrema[neigbourComponents[edge]]);

                // and save the edge itself
                activeEdges[edgeID] = edgeID;
              } // per edge
          }
        } // per activeIndex*/
  }
}; // Mesh2D_DEM_VertexStarter


} // namespace mesh_dem_triangulation_worklets
} // namespace contourtree_augmented
} // namespace worklet
} // namespace vtkm

#endif
