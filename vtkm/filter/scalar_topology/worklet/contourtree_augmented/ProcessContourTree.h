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


#ifndef vtk_m_worklet_contourtree_augmented_process_contourtree_h
#define vtk_m_worklet_contourtree_augmented_process_contourtree_h

// global includes
#include <algorithm>
#include <iomanip>
#include <iostream>

// Additional includes for coordinates:
#include <map>
#include <tuple>

// local includes
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/PrintVectors.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>

//VTKM includes
#include <vtkm/Pair.h>
#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayHandleView.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ArrayTransforms.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ContourTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/PrintVectors.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/Branch.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/SuperArcVolumetricComparator.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/SuperNodeBranchComparator.h>

#include <vtkm/cont/Invoker.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/HypersweepWorklets.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/PointerDoubling.h>

//#define DEBUG_PRINT


namespace process_contourtree_inc_ns =
  vtkm::worklet::contourtree_augmented::process_contourtree_inc;

using ValueType = vtkm::Float32;
using FloatArrayType = vtkm::cont::ArrayHandle<ValueType>;

namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{

struct Coordinates
{
    double x;
    double y;
    double z;
};

struct Triangle
{
    int p1;
    int p2;
    int p3;
};

// Adjusted function to read the file.
//std::map<vtkm::Id, Coordinates>


std::vector<Coordinates> ReadCoordinatesFromFile(const std::string& filename) {
//    std::map<vtkm::Id, Coordinates> coordinatesMap;
    std::vector<Coordinates> coordinatesMap;
    std::ifstream file(filename);
    std::string line;

    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return coordinatesMap; // Return an empty map if file opening fails.
    }

    std::cout << "READING ... \n";

    while (getline(file, line)) {
        std::istringstream iss(line);
        vtkm::Id id;
        Coordinates coords;

        // Extracting ID and coordinates from the current line, assuming the ID is a vtkm::Id.
        if (!(iss >> coords.x >> coords.y >> coords.z)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue; // Skip to the next line if parsing fails.
        }

        std::cout << coords.x << " " << coords.y << " " << coords.z << std::endl;

        // Store the extracted data in the map.
        coordinatesMap.push_back(coords);
    }

    return coordinatesMap;
}


std::vector<Triangle> ReadTrianglesFromFile(const std::string& filename) {
//    std::map<vtkm::Id, Coordinates> coordinatesMap;
    std::vector<Triangle> triangleMap;
    std::ifstream file(filename);
    std::string line;

    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return triangleMap; // Return an empty map if file opening fails.
    }

    std::cout << "READING ... \n";

    while (getline(file, line)) {
        std::istringstream iss(line);
        vtkm::Id id;
        Triangle triang;

        // Extracting ID and coordinates from the current line, assuming the ID is a vtkm::Id.
        if (!(iss >> id >> triang.p1 >> triang.p2 >> triang.p3)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue; // Skip to the next line if parsing fails.
        }

        std::cout << triang.p1 << " " << triang.p2 << " " << triang.p3 << std::endl;

        // Store the extracted data in the map.
        triangleMap.push_back(triang);
    }

    return triangleMap;
}




// TODO Many of the post processing routines still need to be parallelized
// Class with routines for post processing the contour tree
class ProcessContourTree
{ // class ProcessContourTree
public:
  // initialises contour tree arrays - rest is done by another class
  ProcessContourTree()
  { // ProcessContourTree()
  } // ProcessContourTree()

  // collect the sorted arcs
  void static CollectSortedArcs(const ContourTree& contourTree,
                                const IdArrayType& sortOrder,
                                EdgePairArray& sortedArcs)
  { // CollectSortedArcs
    // create an array for sorting the arcs
    std::vector<EdgePair> arcSorter;

    // fill it up
    auto arcsPortal = contourTree.Arcs.ReadPortal();
    auto sortOrderPortal = sortOrder.ReadPortal();

    for (vtkm::Id node = 0; node < contourTree.Arcs.GetNumberOfValues(); node++)
    { // per node
      // retrieve ID of target supernode
      vtkm::Id arcTo = arcsPortal.Get(node);

      // if this is true, it is the last pruned vertex & is omitted
      if (NoSuchElement(arcTo))
        continue;

      // otherwise, strip out the flags
      arcTo = MaskedIndex(arcTo);

      // now convert to mesh IDs from sort IDs
      // otherwise, we need to convert the IDs to regular mesh IDs
      vtkm::Id regularID = sortOrderPortal.Get(node);

      // retrieve the regular ID for it
      vtkm::Id regularTo = sortOrderPortal.Get(arcTo);

      // how we print depends on which end has lower ID
      if (regularID < regularTo)
        arcSorter.push_back(EdgePair(regularID, regularTo));
      else
        arcSorter.push_back(EdgePair(regularTo, regularID));
    } // per vertex

    // now sort it
    // Setting saddlePeak reference to the make_ArrayHandle directly does not work
    sortedArcs = vtkm::cont::make_ArrayHandle(arcSorter, vtkm::CopyFlag::On);
    vtkm::cont::Algorithm::Sort(sortedArcs, SaddlePeakSort());
  } // CollectSortedArcs

  // collect the sorted superarcs
  void static CollectSortedSuperarcs(const ContourTree& contourTree,
                                     const IdArrayType& sortOrder,
                                     EdgePairArray& saddlePeak)
  { // CollectSortedSuperarcs()
    // create an array for sorting the arcs
    std::vector<EdgePair> superarcSorter;

    // fill it up
    auto supernodesPortal = contourTree.Supernodes.ReadPortal();
    auto superarcsPortal = contourTree.Superarcs.ReadPortal();
    auto sortOrderPortal = sortOrder.ReadPortal();

    for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
         supernode++)
    { // per supernode
      // sort ID of the supernode
      vtkm::Id sortID = supernodesPortal.Get(supernode);

      // retrieve ID of target supernode
      vtkm::Id superTo = superarcsPortal.Get(supernode);

      // if this is true, it is the last pruned vertex & is omitted
      if (NoSuchElement(superTo))
        continue;

      // otherwise, strip out the flags
      superTo = MaskedIndex(superTo);

      // otherwise, we need to convert the IDs to regular mesh IDs
      vtkm::Id regularID = sortOrderPortal.Get(MaskedIndex(sortID));

      // retrieve the regular ID for it
      vtkm::Id regularTo = sortOrderPortal.Get(MaskedIndex(supernodesPortal.Get(superTo)));

      // how we print depends on which end has lower ID
      if (regularID < regularTo)
      { // from is lower
        // extra test to catch duplicate edge
        if (superarcsPortal.Get(superTo) != supernode)
        {
          superarcSorter.push_back(EdgePair(regularID, regularTo));
        }
      } // from is lower
      else
      {
        superarcSorter.push_back(EdgePair(regularTo, regularID));
      }
    } // per vertex

    // Setting saddlePeak reference to the make_ArrayHandle directly does not work
    saddlePeak = vtkm::cont::make_ArrayHandle(superarcSorter, vtkm::CopyFlag::On);

    // now sort it
    vtkm::cont::Algorithm::Sort(saddlePeak, SaddlePeakSort());
  } // CollectSortedSuperarcs()



















//// ORIGINAL
//  // routine to compute the volume for each hyperarc and superarc
//  void static ComputeVolumeWeightsSerial(const ContourTree& contourTree,
//                                         const vtkm::Id nIterations,
//                                         IdArrayType& superarcIntrinsicWeight,
//                                         IdArrayType& superarcDependentWeight,
//                                         IdArrayType& supernodeTransferWeight,
//                                         IdArrayType& hyperarcDependentWeight)
//  { // ContourTreeMaker::ComputeWeights()
//    // start by storing the first sorted vertex ID for each superarc
//    IdArrayType firstVertexForSuperparent;
//    firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
//    superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());
//    auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();
//    auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();
//    auto superparentsPortal = contourTree.Superparents.ReadPortal();
//    auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
//    auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
//    auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
//    // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
//    auto nodesPortal = contourTree.Nodes.ReadPortal();
//    // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();
//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    { // per node in sorted order
//      vtkm::Id sortID = nodesPortal.Get(sortedNode);
//      vtkm::Id superparent = superparentsPortal.Get(sortID);

//      std::cout << sortID << " - " << superparent << std::endl;

//      if (sortedNode == 0)
//        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
//      else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
//        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
//    } // per node in sorted order

//    std::cout << std::endl;

//    // now we use that to compute the intrinsic weights
//    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//      if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          contourTree.Arcs.GetNumberOfValues() -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      else
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          firstVertexForSuperparentPortal.Get(superarc + 1) -
//                                            firstVertexForSuperparentPortal.Get(superarc));

//    // now initialise the arrays for transfer & dependent weights
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Superarcs.GetNumberOfValues()),
//      superarcDependentWeight);
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Supernodes.GetNumberOfValues()),
//      supernodeTransferWeight);
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Hyperarcs.GetNumberOfValues()),
//      hyperarcDependentWeight);

//    // set up the array which tracks which supernodes to deal with on which iteration
//    auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
//    auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
//    auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
//    auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
//    auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

//    /*
//    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
//                          firstSupernodePerIteration);
//    auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
//    for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
//         supernode++)
//    { // per supernode
//      vtkm::Id when = MaskedIndex(whenTransferredPortal.Get(supernode));
//      if (supernode == 0)
//      { // zeroth supernode
//        firstSupernodePerIterationPortal.Set(when, supernode);
//      } // zeroth supernode
//      else if (when != MaskedIndex(whenTransferredPortal.Get(supernode - 1)))
//      { // non-matching supernode
//        firstSupernodePerIterationPortal.Set(when, supernode);
//      } // non-matching supernode
//    }   // per supernode
//    for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
//      if (firstSupernodePerIterationPortal.Get(iteration) == 0)
//        firstSupernodePerIterationPortal.Set(iteration,
//                                             firstSupernodePerIterationPortal.Get(iteration + 1));

//    // set the sentinel at the end of the array
//    firstSupernodePerIterationPortal.Set(nIterations, contourTree.Supernodes.GetNumberOfValues());

//    // now use that array to construct a similar array for hypernodes
//    IdArrayType firstHypernodePerIteration;
//    firstHypernodePerIteration.Allocate(nIterations + 1);
//    auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();
//    auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
//    auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
//    auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();
//    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
//      firstHypernodePerIterationPortal.Set(
//        iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
//    firstHypernodePerIterationPortal.Set(nIterations, contourTree.Hypernodes.GetNumberOfValues());
//    */

//    // now iterate, propagating weights inwards
//    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
//    { // per iteration

//      std::cout << "Iteration: " << iteration << std::endl;

//      // pull the array bounds into register
//      vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
//      vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
//      vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
//      vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

//      // Recall that the superarcs are sorted by (iteration, hyperarc), & that all superarcs for a given hyperarc are processed
//      // in the same iteration.  Assume therefore that:
//      //      i. we now have the intrinsic weight assigned for each superarc, and
//      // ii. we also have the transfer weight assigned for each supernode.
//      //
//      // Suppose we have a sequence of superarcs
//      //                      s11 s12 s13 s14 s21 s22 s23 s31
//      // with transfer weights at their origins and intrinsic weights along them
//      //      sArc                     s11 s12 s13 s14 s21 s22 s23 s31
//      //      transfer wt               0   1   2   1   2   3   1   0
//      //      intrinsic wt              1   2   1   5   2   6   1   1
//      //
//      //  now, if we do a prefix sum on each of these and add the two sums together, we get:
//      //      sArc                                  s11 s12 s13 s14 s21 s22 s23 s31
//      //      hyperparent sNode ID                  s11 s11 s11 s11 s21 s21 s21 s31
//      //      transfer weight                       0   1   2   1   2   3   1   0
//      //      intrinsic weight                      1   2   1   5   2   6   1   1
//      //      sum(xfer + intrinsic)                 1   3   3   6   4   9   2   1
//      //  prefix sum (xfer + int)                   1   4   7  13  17  26  28  29
//      //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  16


//      std::cout << "SUPERARCS: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << supernode << " ";
//      }
//      std::cout << std::endl;


//      std::cout << "transfer: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << supernodeTransferWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      std::cout << "intrinsic: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << superarcIntrinsicWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      std::cout << "step 1: ";
//      // so, step 1: add xfer + int & store in dependent weight
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//        superarcDependentWeightPortal.Set(supernode,
//                                          supernodeTransferWeightPortal.Get(supernode) +
//                                            superarcIntrinsicWeightPortal.Get(supernode));

//        std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << " - DEPENDENT = TRANSFER + INTRINSIC" << std::endl;

//      std::cout << "step 2: " << std::endl;
//      std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
//      // step 2: perform prefix sum on the dependent weight range
//      for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
//      {
//        superarcDependentWeightPortal.Set(supernode,
//                                          superarcDependentWeightPortal.Get(supernode) +
//                                            superarcDependentWeightPortal.Get(supernode - 1));
//        //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

//        std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

//      }
//      std::cout << std::endl;
////      std::cout << " - DEPENDENT = DEPENDENT[CURRENT] + DEPENDENT[PREVIOUS]" << std::endl;

//      // step 3: subtract out the dependent weight of the prefix to the entire hyperarc. This will be a transfer, but for now, it's easier
//      // to show it in serial. NB: Loops backwards so that computation uses the correct value
//      // As a bonus, note that we test > firstsupernode, not >=.  This is because we've got unsigned integers, & otherwise it will not terminate
//      // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
//      std::cout << "subtract:\n";
//      for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
//      { // per supernode
//        // retrieve the hyperparent & convert to a supernode ID
//        vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
//        vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

//        // if the hyperparent is the first in the sequence, dependent weight is already correct
//        if (hyperparent == firstHypernode)
//          continue;

//        // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
//        superarcDependentWeightPortal.Set(
//          supernode,
//          superarcDependentWeightPortal.Get(supernode) -
//            superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

//        std::cout << supernode << "(" << hyperparentSuperID << ")" << " - " << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << std::endl;

//      } // per supernode


//      std::cout << "target transfer weights:\n";
//      // step 4: transfer the dependent weight to the hyperarc's target supernode
//      for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
//      { // per hypernode
//        // last superarc for the hyperarc
//        vtkm::Id lastSuperarc;
//        // special case for the last hyperarc
//        if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
//          // take the last superarc in the array
//          lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
//        else
//          // otherwise, take the next hypernode's ID and subtract 1
//          lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

//        // now, given the last superarc for the hyperarc, transfer the dependent weight
//        hyperarcDependentWeightPortal.Set(hypernode,
//                                          superarcDependentWeightPortal.Get(lastSuperarc));

//        // note that in parallel, this will have to be split out as a sort & partial sum in another array
//        vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
//        supernodeTransferWeightPortal.Set(hyperarcTarget,
//                                          supernodeTransferWeightPortal.Get(hyperarcTarget) +
//                                            hyperarcDependentWeightPortal.Get(hypernode));

//        std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

//      } // per hypernode

//      std::cout << std::endl;
//      std::cout << "final:\n";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
////        superarcDependentWeightPortal.Set(supernode,
////                                          superarcDependentWeightPortal.Get(supernode) +
////                                            superarcDependentWeightPortal.Get(supernode - 1));
//        //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

//        std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

//      }
//      std::cout << std::endl;

//    }   // per iteration
//  }     // ContourTreeMaker::ComputeWeights()









    // compute the triangle area from the coordinate points of the triangle:
    double static ComputeTriangleArea(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double x3, double y3, double z3)
//    double static ComputeTriangleArea(Coordinates c1, Coordinates c2, Coordinates c3)
    {
        return 0.5 * abs( (x2-x1)*(y3-y1) - (x3 - x1) * (y2 -y1) );
    }






    // 2024-08-01 COMPUTE THE FLOAT VERSION OF THE WEIGHTS
    void static ComputeVolumeWeightsSerialFloat(const ContourTree& contourTree,
                                                const vtkm::Id nIterations,
                                                FloatArrayType & superarcIntrinsicWeight,
                                                FloatArrayType & superarcDependentWeight,
                                                FloatArrayType & supernodeTransferWeight,
                                                FloatArrayType & hyperarcDependentWeight)
    { // ContourTreeMaker::ComputeWeights()
        // start by storing the first sorted vertex ID for each superarc
        IdArrayType firstVertexForSuperparent;
        firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
        superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();
        auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();
        auto superparentsPortal = contourTree.Superparents.ReadPortal();
        auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
        auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
        auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
        // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto nodesPortal = contourTree.Nodes.ReadPortal();
        // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();


//            const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-coordinates.txt";
//            const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-triangles.txt";

            const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-coordinates.txt";
            const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-triang.txt";

        //    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-coordinates.txt";

        //    std::map<vtkm::Id, Coordinates>
            std::vector<Coordinates> coordlist    = ReadCoordinatesFromFile(filename1);
            std::vector<Triangle> trianglelist = ReadTrianglesFromFile(filename2);

            std::vector<double> weightList;

            std::cout << "PRINT THE ARRAYS OF COORDINATES: \n";

            // Print the read data for demonstration purposes.
        //    for (const auto& pair : coordlist)
            for (int i = 0; i < coordlist.size(); i++)
            {
                std::cout << i << ": " << coordlist[i].x << ", " << coordlist[i].y << ", " << coordlist[i].z << std::endl;
                weightList.push_back(0.0);
            }
            std::cout << "PRINT THE ARRAYS OF TRIANGLES: \n";
            for (int i = 0; i < trianglelist.size(); i++)
            {
                std::cout << i << ": " << trianglelist[i].p1 << ", " << trianglelist[i].p2 << ", " << trianglelist[i].p3; //<< std::endl;
                double area = ComputeTriangleArea(
                                coordlist[trianglelist[i].p1].x, coordlist[trianglelist[i].p1].y, coordlist[trianglelist[i].p1].z,
                                coordlist[trianglelist[i].p2].x, coordlist[trianglelist[i].p2].y, coordlist[trianglelist[i].p2].z,
                                coordlist[trianglelist[i].p3].x, coordlist[trianglelist[i].p3].y, coordlist[trianglelist[i].p3].z);

                //                double area = 0.5;

                double wt = area / 3.0;
                // for each vertex comprising the triangle ...
                // ... add 1/3rd of the triangle's area to the vertice's weight:
                weightList[trianglelist[i].p1] += wt;
                weightList[trianglelist[i].p2] += wt;
                weightList[trianglelist[i].p3] += wt;

                std::cout << " = " << area << std::endl;
            }

            for (int i = 0; i < weightList.size(); i++)
            {
                std::cout << i << ": " << weightList[i] << std::endl;
            }


        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);

            // initialise the intrinsic weight array counter:
//            superarcIntrinsicWeightPortal.Set(superparent, 0);
            superarcIntrinsicWeightPortal.Set(superparent, 0.f);
        }

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        { // per node in sorted order
          vtkm::Id sortID = nodesPortal.Get(sortedNode);
          vtkm::Id superparent = superparentsPortal.Get(sortID);

          std::cout << sortID << " - " << superparent << std::endl;

          if (sortedNode == 0)
            firstVertexForSuperparentPortal.Set(superparent, sortedNode);
          else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
            firstVertexForSuperparentPortal.Set(superparent, sortedNode);


          // CHANGES:
          // UPDATE AT REGULAR NODE: +1
          //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

          // UPDATE AT REGULAR NODE: +area/3
//          superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
//          superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);
          superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);


        } // per node in sorted order

        std::cout << std::endl;

  //      std::cout << "target transfer weights:\n";
  //      // step 4: transfer the dependent weight to the hyperarc's target supernode
  //      for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
  //      { // per hypernode
  //        // last superarc for the hyperarc
  //        vtkm::Id lastSuperarc;
  //        // special case for the last hyperarc
  //        if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
  //          // take the last superarc in the array
  //          lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
  //        else
  //          // otherwise, take the next hypernode's ID and subtract 1
  //          lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

  //        // now, given the last superarc for the hyperarc, transfer the dependent weight
  //        hyperarcDependentWeightPortal.Set(hypernode,
  //                                          superarcDependentWeightPortal.Get(lastSuperarc));

  //        // note that in parallel, this will have to be split out as a sort & partial sum in another array
  //        vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
  //        supernodeTransferWeightPortal.Set(hyperarcTarget,
  //                                          supernodeTransferWeightPortal.Get(hyperarcTarget) +
  //                                            hyperarcDependentWeightPortal.Get(hypernode));

  //        std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

  //      } // per hypernode

  //      // COMMS: old trick to compute the intrinsic wts of branches ...
  //      // COMMS: ... now we replace that with an array pass above
  //      // now we use that to compute the intrinsic weights
  //      for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
  //        if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
  //          superarcIntrinsicWeightPortal.Set(superarc,
  //                                            contourTree.Arcs.GetNumberOfValues() -
  //                                              firstVertexForSuperparentPortal.Get(superarc));
  //        else
  //          superarcIntrinsicWeightPortal.Set(superarc,
  //                                            firstVertexForSuperparentPortal.Get(superarc + 1) -
  //                                              firstVertexForSuperparentPortal.Get(superarc));

        // now initialise the arrays for transfer & dependent weights
        vtkm::cont::ArrayCopy(
          vtkm::cont::ArrayHandleConstant<ValueType>(0.f, contourTree.Superarcs.GetNumberOfValues()),
          superarcDependentWeight);
        vtkm::cont::ArrayCopy(
          vtkm::cont::ArrayHandleConstant<ValueType>(0.f, contourTree.Supernodes.GetNumberOfValues()),
          supernodeTransferWeight);
        vtkm::cont::ArrayCopy(
          vtkm::cont::ArrayHandleConstant<ValueType>(0.f, contourTree.Hyperarcs.GetNumberOfValues()),
          hyperarcDependentWeight);

        // set up the array which tracks which supernodes to deal with on which iteration
        auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
        auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
        auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
        auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
        auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

        /*
        vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
                              firstSupernodePerIteration);
        auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
        for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
             supernode++)
        { // per supernode
          vtkm::Id when = MaskedIndex(whenTransferredPortal.Get(supernode));
          if (supernode == 0)
          { // zeroth supernode
            firstSupernodePerIterationPortal.Set(when, supernode);
          } // zeroth supernode
          else if (when != MaskedIndex(whenTransferredPortal.Get(supernode - 1)))
          { // non-matching supernode
            firstSupernodePerIterationPortal.Set(when, supernode);
          } // non-matching supernode
        }   // per supernode
        for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
          if (firstSupernodePerIterationPortal.Get(iteration) == 0)
            firstSupernodePerIterationPortal.Set(iteration,
                                                 firstSupernodePerIterationPortal.Get(iteration + 1));

        // set the sentinel at the end of the array
        firstSupernodePerIterationPortal.Set(nIterations, contourTree.Supernodes.GetNumberOfValues());

        // now use that array to construct a similar array for hypernodes
        IdArrayType firstHypernodePerIteration;
        firstHypernodePerIteration.Allocate(nIterations + 1);
        auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();
        auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
        auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
        auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();
        for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
          firstHypernodePerIterationPortal.Set(
            iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
        firstHypernodePerIterationPortal.Set(nIterations, contourTree.Hypernodes.GetNumberOfValues());
        */

        // now iterate, propagating weights inwards
        for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
        { // per iteration

          std::cout << "Iteration: " << iteration << std::endl;

          // pull the array bounds into register
          vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
          vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
          vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
          vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

          // Recall that the superarcs are sorted by (iteration, hyperarc), & that all superarcs for a given hyperarc are processed
          // in the same iteration.  Assume therefore that:
          //      i. we now have the intrinsic weight assigned for each superarc, and
          // ii. we also have the transfer weight assigned for each supernode.
          //
          // Suppose we have a sequence of superarcs
          //                      s11 s12 s13 s14 s21 s22 s23 s31
          // with transfer weights at their origins and intrinsic weights along them
          //      sArc                     s11 s12 s13 s14 s21 s22 s23 s31
          //      transfer wt               0   1   2   1   2   3   1   0
          //      intrinsic wt              1   2   1   5   2   6   1   1
          //
          //  now, if we do a prefix sum on each of these and add the two sums together, we get:
          //      sArc                                  s11 s12 s13 s14 s21 s22 s23 s31
          //      hyperparent sNode ID                  s11 s11 s11 s11 s21 s21 s21 s31
          //      transfer weight                       0   1   2   1   2   3   1   0
          //      intrinsic weight                      1   2   1   5   2   6   1   1
          //      sum(xfer + intrinsic)                 1   3   3   6   4   9   2   1
          //  prefix sum (xfer + int)                   1   4   7  13  17  26  28  29
          //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  16


          std::cout << "SUPERARCS: ";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << supernode << " ";
          }
          std::cout << std::endl;


          std::cout << "transfer: ";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << supernodeTransferWeightPortal.Get(supernode) << " ";
          }
          std::cout << std::endl;

          std::cout << "intrinsic: ";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << superarcIntrinsicWeightPortal.Get(supernode) << " ";
          }
          std::cout << std::endl;

          std::cout << "step 1: ";
          // so, step 1: add xfer + int & store in dependent weight
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
            superarcDependentWeightPortal.Set(supernode,
                                              supernodeTransferWeightPortal.Get(supernode) +
                                                superarcIntrinsicWeightPortal.Get(supernode));

            std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
          }
          std::cout << " - DEPENDENT = TRANSFER + INTRINSIC" << std::endl;

          std::cout << "step 2: " << std::endl;
          std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
          // step 2: perform prefix sum on the dependent weight range
          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
          {
            superarcDependentWeightPortal.Set(supernode,
                                              superarcDependentWeightPortal.Get(supernode) +
                                                superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

            std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

          }
          std::cout << std::endl;
    //      std::cout << " - DEPENDENT = DEPENDENT[CURRENT] + DEPENDENT[PREVIOUS]" << std::endl;

          // step 3: subtract out the dependent weight of the prefix to the entire hyperarc. This will be a transfer, but for now, it's easier
          // to show it in serial. NB: Loops backwards so that computation uses the correct value
          // As a bonus, note that we test > firstsupernode, not >=.  This is because we've got unsigned integers, & otherwise it will not terminate
          // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
          std::cout << "subtract:\n";
          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
              continue;

            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
            superarcDependentWeightPortal.Set(
              supernode,
              superarcDependentWeightPortal.Get(supernode) -
                superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

            std::cout << supernode << "(" << hyperparentSuperID << ")" << " - " << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << std::endl;

          } // per supernode


          std::cout << "target transfer weights:\n";
          // step 4: transfer the dependent weight to the hyperarc's target supernode
          for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
          { // per hypernode
            // last superarc for the hyperarc
            vtkm::Id lastSuperarc;
            // special case for the last hyperarc
            if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
              // take the last superarc in the array
              lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
            else
              // otherwise, take the next hypernode's ID and subtract 1
              lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

            // now, given the last superarc for the hyperarc, transfer the dependent weight
            hyperarcDependentWeightPortal.Set(hypernode,
                                              superarcDependentWeightPortal.Get(lastSuperarc));

            // note that in parallel, this will have to be split out as a sort & partial sum in another array
            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
            supernodeTransferWeightPortal.Set(hyperarcTarget,
                                              supernodeTransferWeightPortal.Get(hyperarcTarget) +
                                                hyperarcDependentWeightPortal.Get(hypernode));

            std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

          } // per hypernode

          std::cout << std::endl;
          std::cout << "final:\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
    //        superarcDependentWeightPortal.Set(supernode,
    //                                          superarcDependentWeightPortal.Get(supernode) +
    //                                            superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

            std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

          }
          std::cout << std::endl;

        }   // per iteration

        std::cout << std::endl << "Superarc Intrinsic Weight Portal:" << std::endl;
        for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << superarcIntrinsicWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;

        std::cout << std::endl << "superarc Dependent Weight Portal:" << std::endl;
        for(int i = 0; i < superarcDependentWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << superarcDependentWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;


        std::cout << std::endl << "supernodeTransferWeight Portal:" << std::endl;
        for(int i = 0; i < supernodeTransferWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << supernodeTransferWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;

        std::cout << std::endl << "hyperarcDependentWeight Portal:" << std::endl;
        for(int i = 0; i < hyperarcDependentWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << hyperarcDependentWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;

    }













































    void static print2Darray(std::vector<std::vector<double>> vxtc)
    {
        std::cout << "\n   ";
        for(int j=0; j<vxtc[0].size(); j++)
        {
            std::cout << j << " ";
        }
        std::cout << std::endl;
        for(int i=0; i < vxtc.size(); i++)
        {
            std::cout << i << ") ";
            for(int j=0; j<vxtc[i].size(); j++)
            {
                std::cout<<vxtc[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
    }



    // 2024-08-23 COMPUTE THE FLOAT VERSION OF THE WEIGHTS WITH COEFFICIENTS
    void static ComputeVolumeWeightsSerialFloatCoefficients(const ContourTree& contourTree,
                                                const vtkm::Id nIterations,
                                                FloatArrayType & superarcIntrinsicWeight,
                                                FloatArrayType & superarcDependentWeight,
                                                FloatArrayType & supernodeTransferWeight,
                                                FloatArrayType & hyperarcDependentWeight)
    { // ContourTreeMaker::ComputeWeights()
      // START ComputeVolumeWeightsSerialFloatCoefficients
        // start by storing the first sorted vertex ID for each superarc
        IdArrayType firstVertexForSuperparent;
        firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
        superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();
        auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();
        auto superparentsPortal = contourTree.Superparents.ReadPortal();
        auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
        auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
        auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
        // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto nodesPortal = contourTree.Nodes.ReadPortal();
        // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();

        std::cout << "CALL FROM THE COEFFICIENT-BASED FLOAT FUNCTION" << std::endl;

    //            const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-coordinates.txt";
    //            const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-triangles.txt";

            const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-coordinates.txt";
            const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-triang.txt";

        //    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-coordinates.txt";

        //    std::map<vtkm::Id, Coordinates>
            std::vector<Coordinates> coordlist    = ReadCoordinatesFromFile(filename1);
            std::vector<Triangle> trianglelist = ReadTrianglesFromFile(filename2);

            std::vector<double> weightList;

            std::cout << "PRINT THE ARRAYS OF COORDINATES: \n";

            // Print the read data for demonstration purposes.
        //    for (const auto& pair : coordlist)
            for (int i = 0; i < coordlist.size(); i++)
            {
                std::cout << i << ": " << coordlist[i].x << ", " << coordlist[i].y << ", " << coordlist[i].z << std::endl;
                weightList.push_back(0.0);
            }
            std::cout << "PRINT THE ARRAYS OF TRIANGLES: \n";
            for (int i = 0; i < trianglelist.size(); i++)
            {
                std::cout << i << ": " << trianglelist[i].p1 << ", " << trianglelist[i].p2 << ", " << trianglelist[i].p3; //<< std::endl;
                double area = ComputeTriangleArea(
                                coordlist[trianglelist[i].p1].x, coordlist[trianglelist[i].p1].y, coordlist[trianglelist[i].p1].z,
                                coordlist[trianglelist[i].p2].x, coordlist[trianglelist[i].p2].y, coordlist[trianglelist[i].p2].z,
                                coordlist[trianglelist[i].p3].x, coordlist[trianglelist[i].p3].y, coordlist[trianglelist[i].p3].z);

                //                double area = 0.5;

                double wt = area / 3.0;
                // for each vertex comprising the triangle ...
                // ... add 1/3rd of the triangle's area to the vertice's weight:
                weightList[trianglelist[i].p1] += wt;
                weightList[trianglelist[i].p2] += wt;
                weightList[trianglelist[i].p3] += wt;

                std::cout << " = " << area << std::endl;
            }

            std::cout << "PRINT THE WEIGHTS SO FAR: \n";
            for (int i = 0; i < weightList.size(); i++)
            {
                std::cout << i << ": " << weightList[i] << std::endl;
            }


        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);

            // initialise the intrinsic weight array counter:
            // superarcIntrinsicWeightPortal.Set(superparent, 0);
            superarcIntrinsicWeightPortal.Set(superparent, 0.f);
        }


        // ----------------------------------- PRE-PROCESS ----------------------------------- //
        // Here we basically populate the table as such:
        // vertexID / triangleID
        //          triangle1 triangle2   row sum(del. h1/h2)   col (prefix) sum del. h1/h2
        // deltas:  h1 h2     h1 h2
        // v1       x  y      z  w

        // (below is actually the length of the hypotenuse projected onto the base from the 90-degree angle)
        double fake_area = sqrt(2.0) / 2.0;



        std::cout << "PREPROCESSING STEP:" << std::endl;

        // Initializing the 2-D vector
        std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
                                              std::vector<double> (trianglelist.size(), 0.0));

        std::cout << "\nInitialised vxtch1:" << std::endl;
        print2Darray(vxtch1);

        std::vector<std::vector<double>> vxtch2(contourTree.Arcs.GetNumberOfValues(),
                                              std::vector<double> (trianglelist.size(), 0.0));
        print2Darray(vxtch2);

        for (int i = 0; i < trianglelist.size(); i++)
        {
            // ASSUMPTION #1:
            // FOR EACH TRIANGLE, ITS VERTICES ARE GIVEN IN INCREASING ORDER OF VALUE:
            // HIGH, MID, LOW
            // lowvalue, midvalue, highvalue
            double Lv, Mv, Hv;
            Hv = (double)trianglelist[i].p1;
            Mv = (double)trianglelist[i].p2;
            Lv = (double)trianglelist[i].p3;

            // ASSUMPTION #2:
            // THE POINTS ARE THE SAME AS THEIR VALUES GIVEN IN THE FILE

            // m=h1=slope for a line equation:
            // h1 for equation starting at the LOW vertex:
            double Mml = fake_area / (Mv-Lv);
            // h1 for equation starting at the HIGH vertex:
            double Mmh = fake_area / (Mv-Hv);

            // c=h2=intercept for a line equation:
            // h2 for equation starting at the LOW vertex:
            double cL = -Mml * Lv;
            double cH = -Mmh * Hv;

            // deltas:
            // Low-vertex deltas (same as original since adding from 0)
            double ld1 = Mml;
            double ld2 = cL;

            vxtch1[trianglelist[i].p3][i] = ld1;
            vxtch2[trianglelist[i].p3][i] = ld2;

            // Middle-vertex deltas:
            double md1 = Mmh - ld1;
            double md2 = cH - ld2;

            vxtch1[trianglelist[i].p2][i] = md1;
            vxtch2[trianglelist[i].p2][i] = md2;

            // High-vertex deltas:
            double hd1 = -Mmh;
            double hd2 = -cH;

            vxtch1[trianglelist[i].p1][i] = hd1;
            vxtch2[trianglelist[i].p1][i] = hd2;

//            std::cout << " " << i << ") " << delta_h1 << ", " << delta_h2 << std::endl;
        }

        print2Darray(vxtch1);
        print2Darray(vxtch2);


        // now with deltas computed, also compute the prefix sums of the different components:

        std::vector<double> vx_delta_h1_sum;
        std::vector<double> vx_delta_h2_sum;

        std::vector<double> delta_h1_pfixsum; // = 0.0;
        std::vector<double> delta_h2_pfixsum; // = 0.0;

        delta_h1_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
        delta_h2_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());

        for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
        {
            vx_delta_h1_sum.push_back(0.0);
            vx_delta_h2_sum.push_back(0.0);

            for (int j = 0; j < trianglelist.size(); j++)
            {
                vx_delta_h1_sum[i] += vxtch1[i][j];
                vx_delta_h2_sum[i] += vxtch2[i][j];
            }

            std::cout << "del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " ~ "; //<< std::endl;

            if (i == 0)
            {
                delta_h1_pfixsum[0] = vx_delta_h1_sum[0];
//                delta_h1_pfixsum.push_back(vx_delta_h1_sum[0]);
                delta_h2_pfixsum[0] = vx_delta_h2_sum[0];
//                delta_h2_pfixsum.push_back(vx_delta_h2_sum[0]);
            }
            else
            {
                delta_h1_pfixsum[i] += delta_h1_pfixsum[i-1] + vx_delta_h1_sum[i];
                delta_h2_pfixsum[i] += delta_h2_pfixsum[i-1] + vx_delta_h2_sum[i];
            }

            std::cout << "pfix(h1)= (" << i << ")" << delta_h1_pfixsum[i] << " pfix(h2)=" << delta_h2_pfixsum[i] << std::endl;
        }

//        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//        {
//            vtkm::Id sortID = nodesPortal.Get(sortedNode);
//            vtkm::Id superparent = superparentsPortal.Get(sortID);

//            double delta_h1, delta_h2 = 0.0;

//            std::cout << sortID;

//            // each vertex will have a different coefficient on each triangle
//            // ... we index them with a 2D [vertexID][triangleID] array ...
//            // ... to find the corresponding coefficients:
//            // Initializing the 2-D vector
////            std::vector<std::vector<int>> vxtc(contourTree.Arcs.GetNumberOfValues(),
////                                     trianglelist.size());

//            // Initializing the 2-D vector
//            std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
//                                                  std::vector<double> (trianglelist.size(), 0.0));

//            std::cout << "\nInitialised vxtch1:" << std::endl;
//            (vxtch1);

//            std::vector<std::vector<double>> vxtch2(contourTree.Arcs.GetNumberOfValues(),
//                                                  std::vector<double> (trianglelist.size(), 0.0));
//            print2Darray(vxtch2);

//            // (vertex IDs are the same as vertex values for now which makes things easier)

//            // first initialise the array to all 0s:
////            for

//            for (int i = 0; i < trianglelist.size(); i++)
//            {
//                // ASSUMPTION #1:
//                // FOR EACH TRIANGLE, ITS VERTICES ARE GIVEN IN INCREASING ORDER OF VALUE:
//                // HIGH, MID, LOW
//                // lowvalue, midvalue, highvalue
//                double Lv, Mv, Hv;
//                Hv = (double)trianglelist[i].p1;
//                Mv = (double)trianglelist[i].p2;
//                Lv = (double)trianglelist[i].p3;

//                // ASSUMPTION #2:
//                // THE POINTS ARE THE SAME AS THEIR VALUES GIVEN IN THE FILE

//                // m=h1=slope for a line equation:
//                // h1 for equation starting at the LOW vertex:
//                double Mml = fake_area / (Mv-Lv);
//                // h1 for equation starting at the HIGH vertex:
//                double Mmh = fake_area / (Mv-Hv);


//                // c=h2=intercept for a line equation:
//                // h2 for equation starting at the LOW vertex:
//                double cL = -Mml * Lv;
//                double cH = -Mmh * Hv;

//                // deltas:
//                // Low-vertex deltas (same as original since adding from 0)
//                double ld1 = Mml;
//                double ld2 = cL;

//                // Middle-vertex deltas:
//                double md1 = Mmh - ld1;
//                double md2 = cH - ld2;

//                // High-vertex deltas:
//                double hd1 = -Mmh;
//                double hd2 = -cH;

//                std::cout << " " << i << ") " << delta_h1 << ", " << delta_h2 << std::endl;
//            }
//        }



        // ----------------------------------- PRE-PROCESS ----------------------------------- //



        // for the sweep, we will be using the pre-computed delta coefficients from the mesh

        // -------------------------------------- SWEEP  ------------------------------------- //

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        { // per node in sorted order
          vtkm::Id sortID = nodesPortal.Get(sortedNode);
          vtkm::Id superparent = superparentsPortal.Get(sortID);

          std::cout << sortID << " - " << superparent << std::endl;

          if (sortedNode == 0)
            firstVertexForSuperparentPortal.Set(superparent, sortedNode);
          else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
            firstVertexForSuperparentPortal.Set(superparent, sortedNode);


          // CHANGES:
          // UPDATE AT REGULAR NODE: +1
          //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

          // UPDATE AT REGULAR NODE: +area/3
    //          superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
    //          superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);
          superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);


        } // per node in sorted order

        // -------------------------------------- SWEEP  ------------------------------------- //

        std::cout << std::endl;

    //      std::cout << "target transfer weights:\n";
    //      // step 4: transfer the dependent weight to the hyperarc's target supernode
    //      for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
    //      { // per hypernode
    //        // last superarc for the hyperarc
    //        vtkm::Id lastSuperarc;
    //        // special case for the last hyperarc
    //        if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
    //          // take the last superarc in the array
    //          lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
    //        else
    //          // otherwise, take the next hypernode's ID and subtract 1
    //          lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

    //        // now, given the last superarc for the hyperarc, transfer the dependent weight
    //        hyperarcDependentWeightPortal.Set(hypernode,
    //                                          superarcDependentWeightPortal.Get(lastSuperarc));

    //        // note that in parallel, this will have to be split out as a sort & partial sum in another array
    //        vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
    //        supernodeTransferWeightPortal.Set(hyperarcTarget,
    //                                          supernodeTransferWeightPortal.Get(hyperarcTarget) +
    //                                            hyperarcDependentWeightPortal.Get(hypernode));

    //        std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

    //      } // per hypernode

    //      // COMMS: old trick to compute the intrinsic wts of branches ...
    //      // COMMS: ... now we replace that with an array pass above
    //      // now we use that to compute the intrinsic weights
    //      for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
    //        if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
    //          superarcIntrinsicWeightPortal.Set(superarc,
    //                                            contourTree.Arcs.GetNumberOfValues() -
    //                                              firstVertexForSuperparentPortal.Get(superarc));
    //        else
    //          superarcIntrinsicWeightPortal.Set(superarc,
    //                                            firstVertexForSuperparentPortal.Get(superarc + 1) -
    //                                              firstVertexForSuperparentPortal.Get(superarc));

        // now initialise the arrays for transfer & dependent weights
        vtkm::cont::ArrayCopy(
          vtkm::cont::ArrayHandleConstant<ValueType>(0.f, contourTree.Superarcs.GetNumberOfValues()),
          superarcDependentWeight);
        vtkm::cont::ArrayCopy(
          vtkm::cont::ArrayHandleConstant<ValueType>(0.f, contourTree.Supernodes.GetNumberOfValues()),
          supernodeTransferWeight);
        vtkm::cont::ArrayCopy(
          vtkm::cont::ArrayHandleConstant<ValueType>(0.f, contourTree.Hyperarcs.GetNumberOfValues()),
          hyperarcDependentWeight);

        // set up the array which tracks which supernodes to deal with on which iteration
        auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
        auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
        auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
        auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
        auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

        /*
        vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
                              firstSupernodePerIteration);
        auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
        for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
             supernode++)
        { // per supernode
          vtkm::Id when = MaskedIndex(whenTransferredPortal.Get(supernode));
          if (supernode == 0)
          { // zeroth supernode
            firstSupernodePerIterationPortal.Set(when, supernode);
          } // zeroth supernode
          else if (when != MaskedIndex(whenTransferredPortal.Get(supernode - 1)))
          { // non-matching supernode
            firstSupernodePerIterationPortal.Set(when, supernode);
          } // non-matching supernode
        }   // per supernode
        for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
          if (firstSupernodePerIterationPortal.Get(iteration) == 0)
            firstSupernodePerIterationPortal.Set(iteration,
                                                 firstSupernodePerIterationPortal.Get(iteration + 1));

        // set the sentinel at the end of the array
        firstSupernodePerIterationPortal.Set(nIterations, contourTree.Supernodes.GetNumberOfValues());

        // now use that array to construct a similar array for hypernodes
        IdArrayType firstHypernodePerIteration;
        firstHypernodePerIteration.Allocate(nIterations + 1);
        auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();
        auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
        auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
        auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();
        for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
          firstHypernodePerIterationPortal.Set(
            iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
        firstHypernodePerIterationPortal.Set(nIterations, contourTree.Hypernodes.GetNumberOfValues());
        */

        // now iterate, propagating weights inwards
        for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
        { // per iteration

          std::cout << "Iteration: " << iteration << std::endl;

          // pull the array bounds into register
          vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
          vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
          vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
          vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

          // Recall that the superarcs are sorted by (iteration, hyperarc), & that all superarcs for a given hyperarc are processed
          // in the same iteration.  Assume therefore that:
          //      i. we now have the intrinsic weight assigned for each superarc, and
          // ii. we also have the transfer weight assigned for each supernode.
          //
          // Suppose we have a sequence of superarcs
          //                      s11 s12 s13 s14 s21 s22 s23 s31
          // with transfer weights at their origins and intrinsic weights along them
          //      sArc                     s11 s12 s13 s14 s21 s22 s23 s31
          //      transfer wt               0   1   2   1   2   3   1   0
          //      intrinsic wt              1   2   1   5   2   6   1   1
          //
          //  now, if we do a prefix sum on each of these and add the two sums together, we get:
          //      sArc                                  s11 s12 s13 s14 s21 s22 s23 s31
          //      hyperparent sNode ID                  s11 s11 s11 s11 s21 s21 s21 s31
          //      transfer weight                       0   1   2   1   2   3   1   0
          //      intrinsic weight                      1   2   1   5   2   6   1   1
          //      sum(xfer + intrinsic)                 1   3   3   6   4   9   2   1
          //  prefix sum (xfer + int)                   1   4   7  13  17  26  28  29
          //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  16


          std::cout << "SUPERARCS: ";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << supernode << " ";
          }
          std::cout << std::endl;


          std::cout << "transfer: ";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << supernodeTransferWeightPortal.Get(supernode) << " ";
          }
          std::cout << std::endl;

          std::cout << "intrinsic: ";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << superarcIntrinsicWeightPortal.Get(supernode) << " ";
          }
          std::cout << std::endl;

          std::cout << "step 1: ";
          // so, step 1: add xfer + int & store in dependent weight
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
            superarcDependentWeightPortal.Set(supernode,
                                              supernodeTransferWeightPortal.Get(supernode) +
                                                superarcIntrinsicWeightPortal.Get(supernode));

            std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
          }
          std::cout << " - DEPENDENT = TRANSFER + INTRINSIC" << std::endl;

          std::cout << "step 2: " << std::endl;
          std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
          // step 2: perform prefix sum on the dependent weight range
          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
          {
            superarcDependentWeightPortal.Set(supernode,
                                              superarcDependentWeightPortal.Get(supernode) +
                                                superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

            std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

          }
          std::cout << std::endl;
    //      std::cout << " - DEPENDENT = DEPENDENT[CURRENT] + DEPENDENT[PREVIOUS]" << std::endl;

          // step 3: subtract out the dependent weight of the prefix to the entire hyperarc. This will be a transfer, but for now, it's easier
          // to show it in serial. NB: Loops backwards so that computation uses the correct value
          // As a bonus, note that we test > firstsupernode, not >=.  This is because we've got unsigned integers, & otherwise it will not terminate
          // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
          std::cout << "subtract:\n";
          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
              continue;

            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
            superarcDependentWeightPortal.Set(
              supernode,
              superarcDependentWeightPortal.Get(supernode) -
                superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

            std::cout << supernode << "(" << hyperparentSuperID << ")" << " - " << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << std::endl;

          } // per supernode


          std::cout << "target transfer weights:\n";
          // step 4: transfer the dependent weight to the hyperarc's target supernode
          for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
          { // per hypernode
            // last superarc for the hyperarc
            vtkm::Id lastSuperarc;
            // special case for the last hyperarc
            if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
              // take the last superarc in the array
              lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
            else
              // otherwise, take the next hypernode's ID and subtract 1
              lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

            // now, given the last superarc for the hyperarc, transfer the dependent weight
            hyperarcDependentWeightPortal.Set(hypernode,
                                              superarcDependentWeightPortal.Get(lastSuperarc));

            // note that in parallel, this will have to be split out as a sort & partial sum in another array
            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
            supernodeTransferWeightPortal.Set(hyperarcTarget,
                                              supernodeTransferWeightPortal.Get(hyperarcTarget) +
                                                hyperarcDependentWeightPortal.Get(hypernode));

            std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

          } // per hypernode

          std::cout << std::endl;
          std::cout << "final:\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
    //        superarcDependentWeightPortal.Set(supernode,
    //                                          superarcDependentWeightPortal.Get(supernode) +
    //                                            superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

            std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

          }
          std::cout << std::endl;

        }   // per iteration

        std::cout << std::endl << "Superarc Intrinsic Weight Portal:" << std::endl;
        for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << superarcIntrinsicWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;

        std::cout << std::endl << "superarc Dependent Weight Portal:" << std::endl;
        for(int i = 0; i < superarcDependentWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << superarcDependentWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;


        std::cout << std::endl << "supernodeTransferWeight Portal:" << std::endl;
        for(int i = 0; i < supernodeTransferWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << supernodeTransferWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;

        std::cout << std::endl << "hyperarcDependentWeight Portal:" << std::endl;
        for(int i = 0; i < hyperarcDependentWeightPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << " -> " << hyperarcDependentWeightPortal.Get(i) << std::endl;
        }
        std::cout << std::endl;

        std::cout << "END ComputeVolumeWeightsSerialFloatCoefficients" << std::endl;

    }  // END ComputeVolumeWeightsSerialFloatCoefficients



























    struct ContourLengthCoef
    {
        double v1h1; // slope     v1
        double v1h2; // intercept v2
        double v2h1; // slope     v1
        double v2h2; // intercept v2
    };


    struct ContourAreaCoef
    {
        double h1; // a
        double h2; // b
        double h3; // c
    };



    ContourLengthCoef static Compute1DCoeffs(Triangle singleT)
    {// Pass in a single triangle (singleT) and compute itsssss
     //       double area = ComputeTriangleArea(singleT);

        double area = (sqrt(2.0)/2.0);

        // m1 slope
        double v1h1m = area / (singleT.p2 - singleT.p3);
        double v2h1m = area / (singleT.p2 - singleT.p1);
        // m1 slope
        double v1h2c = -v1h1m * singleT.p3;
        double v2h2c = -v2h1m * singleT.p3;

        return {v1h1m, v2h1m, v1h2c, v2h2c};
    }




  // 2024-07-28 MODIFICATION
    // routine to compute the volume for each hyperarc and superarc
    void static ComputeVolumeWeightsSerial(const ContourTree& contourTree,
                                           const vtkm::Id nIterations,
                                           IdArrayType& superarcIntrinsicWeight,
                                           IdArrayType& superarcDependentWeight,
                                           IdArrayType& supernodeTransferWeight,
                                           IdArrayType& hyperarcDependentWeight)
    { // ContourTreeMaker::ComputeWeights()
      // start by storing the first sorted vertex ID for each superarc
      IdArrayType firstVertexForSuperparent;
      firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
      superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());
      auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();
      auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();
      auto superparentsPortal = contourTree.Superparents.ReadPortal();
      auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
      auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
      auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
      // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
      auto nodesPortal = contourTree.Nodes.ReadPortal();
      // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();


//          const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-coordinates.txt";
//          const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-triangles.txt";

          const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-coordinates.txt";
          const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-triang.txt";

      //    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-coordinates.txt";

      //    std::map<vtkm::Id, Coordinates>
          std::vector<Coordinates> coordlist    = ReadCoordinatesFromFile(filename1);
          std::vector<Triangle> trianglelist = ReadTrianglesFromFile(filename2);

          std::vector<double> weightList;

          std::cout << "PRINT THE ARRAYS OF COORDINATES: \n";

          // Print the read data for demonstration purposes.
      //    for (const auto& pair : coordlist)
          for (int i = 0; i < coordlist.size(); i++)
          {
              std::cout << i << ": " << coordlist[i].x << ", " << coordlist[i].y << ", " << coordlist[i].z << std::endl;
              weightList.push_back(0.0);
          }
          std::cout << "PRINT THE ARRAYS OF TRIANGLES: \n";
          for (int i = 0; i < trianglelist.size(); i++)
          {
              std::cout << std::endl << i << ": " << trianglelist[i].p1 << ", " << trianglelist[i].p2 << ", " << trianglelist[i].p3; //<< std::endl;
              double area = ComputeTriangleArea(coordlist[trianglelist[i].p1].x, coordlist[trianglelist[i].p1].y, coordlist[trianglelist[i].p1].z,
                                      coordlist[trianglelist[i].p2].x, coordlist[trianglelist[i].p2].y, coordlist[trianglelist[i].p2].z,
                                      coordlist[trianglelist[i].p3].x, coordlist[trianglelist[i].p3].y, coordlist[trianglelist[i].p3].z);

              double wt = area / 3.0;
              weightList[trianglelist[i].p1] += wt;
              weightList[trianglelist[i].p2] += wt;
              weightList[trianglelist[i].p3] += wt;

              std::cout << " = " << area; // << std::endl;

              ContourLengthCoef vxh1sh2s = Compute1DCoeffs(trianglelist[i]);

              std::cout << "\n1) h1 h2 " << vxh1sh2s.v1h1 << " " << vxh1sh2s.v1h2 << std::endl;
              std::cout << "2) h1 h2 " << vxh1sh2s.v2h1 << " " << vxh1sh2s.v2h2 << std::endl;


          }

          for (int i = 0; i < weightList.size(); i++)
          {
              std::cout << i << ": " << weightList[i] << std::endl;
          }


      for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
      {
          vtkm::Id sortID = nodesPortal.Get(sortedNode);
          vtkm::Id superparent = superparentsPortal.Get(sortID);

          // initialise the transfer weight array counter:
          superarcIntrinsicWeightPortal.Set(superparent, 0);
      }

      for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
      { // per node in sorted order
        vtkm::Id sortID = nodesPortal.Get(sortedNode);
        vtkm::Id superparent = superparentsPortal.Get(sortID);

        std::cout << sortID << " - " << superparent << std::endl;

        if (sortedNode == 0)
          firstVertexForSuperparentPortal.Set(superparent, sortedNode);
        else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
          firstVertexForSuperparentPortal.Set(superparent, sortedNode);


        // CHANGES:
        // UPDATE AT REGULAR NODE: +1
        //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

        // UPDATE AT REGULAR NODE: +area/3
        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);


      } // per node in sorted order

      std::cout << std::endl;

//      std::cout << "target transfer weights:\n";
//      // step 4: transfer the dependent weight to the hyperarc's target supernode
//      for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
//      { // per hypernode
//        // last superarc for the hyperarc
//        vtkm::Id lastSuperarc;
//        // special case for the last hyperarc
//        if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
//          // take the last superarc in the array
//          lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
//        else
//          // otherwise, take the next hypernode's ID and subtract 1
//          lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

//        // now, given the last superarc for the hyperarc, transfer the dependent weight
//        hyperarcDependentWeightPortal.Set(hypernode,
//                                          superarcDependentWeightPortal.Get(lastSuperarc));

//        // note that in parallel, this will have to be split out as a sort & partial sum in another array
//        vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
//        supernodeTransferWeightPortal.Set(hyperarcTarget,
//                                          supernodeTransferWeightPortal.Get(hyperarcTarget) +
//                                            hyperarcDependentWeightPortal.Get(hypernode));

//        std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

//      } // per hypernode

//      // COMMS: old trick to compute the intrinsic wts of branches ...
//      // COMMS: ... now we replace that with an array pass above
//      // now we use that to compute the intrinsic weights
//      for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//        if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
//          superarcIntrinsicWeightPortal.Set(superarc,
//                                            contourTree.Arcs.GetNumberOfValues() -
//                                              firstVertexForSuperparentPortal.Get(superarc));
//        else
//          superarcIntrinsicWeightPortal.Set(superarc,
//                                            firstVertexForSuperparentPortal.Get(superarc + 1) -
//                                              firstVertexForSuperparentPortal.Get(superarc));

      // now initialise the arrays for transfer & dependent weights
      vtkm::cont::ArrayCopy(
        vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Superarcs.GetNumberOfValues()),
        superarcDependentWeight);
      vtkm::cont::ArrayCopy(
        vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Supernodes.GetNumberOfValues()),
        supernodeTransferWeight);
      vtkm::cont::ArrayCopy(
        vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Hyperarcs.GetNumberOfValues()),
        hyperarcDependentWeight);

      // set up the array which tracks which supernodes to deal with on which iteration
      auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
      auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
      auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
      auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
      auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

      /*
      vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
                            firstSupernodePerIteration);
      auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
      for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
           supernode++)
      { // per supernode
        vtkm::Id when = MaskedIndex(whenTransferredPortal.Get(supernode));
        if (supernode == 0)
        { // zeroth supernode
          firstSupernodePerIterationPortal.Set(when, supernode);
        } // zeroth supernode
        else if (when != MaskedIndex(whenTransferredPortal.Get(supernode - 1)))
        { // non-matching supernode
          firstSupernodePerIterationPortal.Set(when, supernode);
        } // non-matching supernode
      }   // per supernode
      for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
        if (firstSupernodePerIterationPortal.Get(iteration) == 0)
          firstSupernodePerIterationPortal.Set(iteration,
                                               firstSupernodePerIterationPortal.Get(iteration + 1));

      // set the sentinel at the end of the array
      firstSupernodePerIterationPortal.Set(nIterations, contourTree.Supernodes.GetNumberOfValues());

      // now use that array to construct a similar array for hypernodes
      IdArrayType firstHypernodePerIteration;
      firstHypernodePerIteration.Allocate(nIterations + 1);
      auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();
      auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
      auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
      auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();
      for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
        firstHypernodePerIterationPortal.Set(
          iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
      firstHypernodePerIterationPortal.Set(nIterations, contourTree.Hypernodes.GetNumberOfValues());
      */

      // now iterate, propagating weights inwards
      for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
      { // per iteration

        std::cout << "Iteration: " << iteration << std::endl;

        // pull the array bounds into register
        vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
        vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
        vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
        vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

        // Recall that the superarcs are sorted by (iteration, hyperarc), & that all superarcs for a given hyperarc are processed
        // in the same iteration.  Assume therefore that:
        //      i. we now have the intrinsic weight assigned for each superarc, and
        // ii. we also have the transfer weight assigned for each supernode.
        //
        // Suppose we have a sequence of superarcs
        //                      s11 s12 s13 s14 s21 s22 s23 s31
        // with transfer weights at their origins and intrinsic weights along them
        //      sArc                     s11 s12 s13 s14 s21 s22 s23 s31
        //      transfer wt               0   1   2   1   2   3   1   0
        //      intrinsic wt              1   2   1   5   2   6   1   1
        //
        //  now, if we do a prefix sum on each of these and add the two sums together, we get:
        //      sArc                                  s11 s12 s13 s14 s21 s22 s23 s31
        //      hyperparent sNode ID                  s11 s11 s11 s11 s21 s21 s21 s31
        //      transfer weight                       0   1   2   1   2   3   1   0
        //      intrinsic weight                      1   2   1   5   2   6   1   1
        //      sum(xfer + intrinsic)                 1   3   3   6   4   9   2   1
        //  prefix sum (xfer + int)                   1   4   7  13  17  26  28  29
        //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  16


        std::cout << "SUPERARCS: ";
        for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
        {
            std::cout << supernode << " ";
        }
        std::cout << std::endl;


        std::cout << "transfer: ";
        for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
        {
            std::cout << supernodeTransferWeightPortal.Get(supernode) << " ";
        }
        std::cout << std::endl;

        std::cout << "intrinsic: ";
        for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
        {
            std::cout << superarcIntrinsicWeightPortal.Get(supernode) << " ";
        }
        std::cout << std::endl;

        std::cout << "step 1: ";
        // so, step 1: add xfer + int & store in dependent weight
        for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
        {
          superarcDependentWeightPortal.Set(supernode,
                                            supernodeTransferWeightPortal.Get(supernode) +
                                              superarcIntrinsicWeightPortal.Get(supernode));

          std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
        }
        std::cout << " - DEPENDENT = TRANSFER + INTRINSIC" << std::endl;

        std::cout << "step 2: " << std::endl;
        std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
        // step 2: perform prefix sum on the dependent weight range
        for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
        {
          superarcDependentWeightPortal.Set(supernode,
                                            superarcDependentWeightPortal.Get(supernode) +
                                              superarcDependentWeightPortal.Get(supernode - 1));
          //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

          std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

        }
        std::cout << std::endl;
  //      std::cout << " - DEPENDENT = DEPENDENT[CURRENT] + DEPENDENT[PREVIOUS]" << std::endl;

        // step 3: subtract out the dependent weight of the prefix to the entire hyperarc. This will be a transfer, but for now, it's easier
        // to show it in serial. NB: Loops backwards so that computation uses the correct value
        // As a bonus, note that we test > firstsupernode, not >=.  This is because we've got unsigned integers, & otherwise it will not terminate
        // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
        std::cout << "subtract:\n";
        for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
        { // per supernode
          // retrieve the hyperparent & convert to a supernode ID
          vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
          vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

          // if the hyperparent is the first in the sequence, dependent weight is already correct
          if (hyperparent == firstHypernode)
            continue;

          // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
          superarcDependentWeightPortal.Set(
            supernode,
            superarcDependentWeightPortal.Get(supernode) -
              superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

          std::cout << supernode << "(" << hyperparentSuperID << ")" << " - " << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << std::endl;

        } // per supernode


        std::cout << "target transfer weights:\n";
        // step 4: transfer the dependent weight to the hyperarc's target supernode
        for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
        { // per hypernode
          // last superarc for the hyperarc
          vtkm::Id lastSuperarc;
          // special case for the last hyperarc
          if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
            // take the last superarc in the array
            lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
          else
            // otherwise, take the next hypernode's ID and subtract 1
            lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

          // now, given the last superarc for the hyperarc, transfer the dependent weight
          hyperarcDependentWeightPortal.Set(hypernode,
                                            superarcDependentWeightPortal.Get(lastSuperarc));

          // note that in parallel, this will have to be split out as a sort & partial sum in another array
          vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
          supernodeTransferWeightPortal.Set(hyperarcTarget,
                                            supernodeTransferWeightPortal.Get(hyperarcTarget) +
                                              hyperarcDependentWeightPortal.Get(hypernode));

          std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

        } // per hypernode

        std::cout << std::endl;
        std::cout << "final:\n";
        for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
        {
  //        superarcDependentWeightPortal.Set(supernode,
  //                                          superarcDependentWeightPortal.Get(supernode) +
  //                                            superarcDependentWeightPortal.Get(supernode - 1));
          //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

          std::cout << supernode << " - " << superarcDependentWeightPortal.Get(supernode) << std::endl;

        }
        std::cout << std::endl;

      }   // per iteration

      std::cout << std::endl << "Superarc Intrinsic Weight Portal:" << std::endl;
      for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << superarcIntrinsicWeightPortal.Get(i) << std::endl;
      }
      std::cout << std::endl;

      std::cout << std::endl << "superarc Dependent Weight Portal:" << std::endl;
      for(int i = 0; i < superarcDependentWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << superarcDependentWeightPortal.Get(i) << std::endl;
      }
      std::cout << std::endl;


      std::cout << std::endl << "supernodeTransferWeight Portal:" << std::endl;
      for(int i = 0; i < supernodeTransferWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << supernodeTransferWeightPortal.Get(i) << std::endl;
      }
      std::cout << std::endl;

      std::cout << std::endl << "hyperarcDependentWeight Portal:" << std::endl;
      for(int i = 0; i < hyperarcDependentWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << hyperarcDependentWeightPortal.Get(i) << std::endl;
      }
      std::cout << std::endl;



    }     // ContourTreeMaker::ComputeWeights()







































































//  // 2024-07-25 MODIFICATION:
//  // routine to compute the volume for each hyperarc and superarc
//  void static ComputeVolumeWeightsSerial(const ContourTree& contourTree,
//                                         const vtkm::Id nIterations,
//                                         IdArrayType& superarcIntrinsicWeight,
//                                         IdArrayType& superarcDependentWeight,
//                                         IdArrayType& supernodeTransferWeight,
//                                         IdArrayType& hyperarcDependentWeight)
//  { // ContourTreeMaker::ComputeWeights()
//    std::cout << "ComputeVolumeWeightsSerial()\n";
//    std::cout << "NUMBER OF ITERATIONS: " << nIterations << std::endl;
//    // start by storing the first sorted vertex ID for each superarc
//    IdArrayType firstVertexForSuperparent;
//    IdArrayType prefixSumSuperarc;
//    // COMMS: One vertex per superarc, so allocate the size = number of superarcs
//    firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
//    prefixSumSuperarc.Allocate(contourTree.Superarcs.GetNumberOfValues());

//    // COMMS: Each superarc will have one intrinsic weight (sum of all the regular nodes (arcs) on it
//    superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());

//    // COMMS: set up all the arrays we will need to write to:
//    auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();
//    auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();

//    auto prefixSumSuperarcPortal    = prefixSumSuperarc.WritePortal();
//    auto prefixSumSuperarcGetPortal = prefixSumSuperarc.ReadPortal();

//    auto superparentsPortal = contourTree.Superparents.ReadPortal();
//    auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
//    auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
//    auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
//    // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
//    auto nodesPortal = contourTree.Nodes.ReadPortal();

//    // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();


//    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-coordinates.txt";

////    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/5b-coordinates.txt";

////    std::map<vtkm::Id, Coordinates>
//    std::vector<Coordinates> coordlist = ReadCoordinatesFromFile(filename);

//    std::cout << "PRINT THE ARRAYS OF COORDINATES: \n";

//    // Print the read data for demonstration purposes.
////    for (const auto& pair : coordlist)
//    for (int i = 0; i < coordlist.size(); i++)
//    {
//        std::cout << i << ": " << coordlist[i].x << ", " << coordlist[i].y << ", " << coordlist[i].z << std::endl;
//    }





//    // ======================================= v SWEEP HERE v ========================================= //

//    std::vector<float> vertex_weights(contourTree.Arcs.GetNumberOfValues()); // same as contourArcs.size()

//    std::cout << "Sorted nodes:" << std::endl;
//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    { // per node in sorted order
//      vtkm::Id sortID = nodesPortal.Get(sortedNode);

//      std::cout << sortID << " ";

//      // HERE WE WILL GET TRIANGLE AREAS ETC.
//      // NOW JUST SET TO 1.f TO BE SAME AS THE SEGMENTED ARRAY PREFIX SUM
//      vertex_weights[sortID] = 1.f;
//    }
//    std::cout << std::endl;

//    std::cout << "weight array:" << std::endl;
//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    {
//        std::cout << vertex_weights[sortedNode] << " ";
//    }

//    // ======================================= ^ SWEEP HERE ^ ========================================= //


//    // ---------------------------------- TASK 1: v REPLACE THIS v ---------------------------------- //
//    // COMMS: STEP 1: Create a segmented array of superparent IDs ...
//    // COMMS: ... for each node, by its sorted ID populate the SP Segmented Array
//    // COMMS: The Segmented Superparent ID array will look like this:
//    // COMMS: ... superparent ID                [ 11, 12, 13, 14, 21, 22, 23, 31 ]
//    // COMMS: ... first node sorted ID in arc:  [  0,  1,  3,  4,  9, 11, 17, 18 ]

//    std::cout << std::endl << "NUM OF ARCS IN CT: << " << contourTree.Arcs.GetNumberOfValues() << std::endl;

//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    { // per node in sorted order

//      // get the node's sort ID:
//      vtkm::Id sortID = nodesPortal.Get(sortedNode);
//      // ... and then get its superparent's ID:
//      vtkm::Id superparent = superparentsPortal.Get(sortID);

//      std::cout << sortedNode << "(" << sortID << ")" << " - > " << superparent << std::endl;

//      if (sortedNode == 0)
//      {
//        firstVertexForSuperparentPortal.Set(superparent, sortID); //sortedNode);
//        prefixSumSuperarcPortal.Set(superparent, sortID);
//      }
//      else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
//      { // COMMS: This gets executed when we're already on a new branch/superarc
//        // COMMS:
//        // if we have the last entry of the array that is not already the superparent in the array ...
//        // ... add a new superparent only then
//        firstVertexForSuperparentPortal.Set(superparent, sortID); //sortedNode);

//        // get data about the last superarc:
//        // previous node:
////        vtkm::Id sortedNodePrev = sortedNode-1;
//        // get the previous node's sort ID:
//        vtkm::Id sortIDPrev = nodesPortal.Get(sortedNode - 1);
//        // ... and then get the previous node's superparent's ID:
//        vtkm::Id superparentPrev = superparentsPortal.Get(sortIDPrev);

//        std::cout << firstVertexForSuperparentPortal.Get(superparentPrev) << " - " << sortIDPrev << " = "
//                  << firstVertexForSuperparentPortal.Get(superparentPrev) - sortIDPrev <<  std::endl;

//        prefixSumSuperarcPortal.Set(superparentPrev, abs(firstVertexForSuperparentPortal.Get(superparentPrev) - sortIDPrev));
//      }

//      if (sortedNode == contourTree.Arcs.GetNumberOfValues()-1)
//      {
//        prefixSumSuperarcPortal.Set(superparent, 0);
//      }

//    } // per node in sorted order

////    DOES NOT WORK HERE ...
////    PrintValues(
////      "firstVertexForSuperparentPortal print:",
////      firstVertexForSuperparent,
////      -1,
////      std::cout);

//    auto firstVertexForSuperparentReadPortal = firstVertexForSuperparent.ReadPortal();

//    std::cout << std::endl << "First Vertices of Each Superarc" << std::endl;
//    for(int i = 0; i < firstVertexForSuperparentReadPortal.GetNumberOfValues(); i++)
//    {
//        std::cout << firstVertexForSuperparentReadPortal.Get(i) << " ";
//    }
//    std::cout << std::endl << std::endl;

//    std::cout << std::endl << "Weights by value of each Superarc" << std::endl;
//    for(int i = 0; i < firstVertexForSuperparentReadPortal.GetNumberOfValues(); i++)
//    {
//        std::cout << prefixSumSuperarcGetPortal.Get(i) << " ";
//    }
//    std::cout << std::endl << std::endl;



//    // COMMS: STEP 2
//    // now we use that to compute the intrinsic weights
//    // COMMS: (The intrinsic weight = the number of regular vertices on the superarc)
//    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//    {
//      if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
//      { // COMMS: deal with the last case manually, ...
//        // COMMS: ... because there is no right nbor, use absolute last arc ID number
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          contourTree.Arcs.GetNumberOfValues() -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      }
//      else
//      { // COMMS: look at the right nbor, get the difference in sort IDs, know the count
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          firstVertexForSuperparentPortal.Get(superarc + 1) -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      }
//    }

//    std::cout << std::endl << "Superarc Intrinsic Weight Portal:" << std::endl;
//    for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
//    {
//        std::cout << superarcIntrinsicWeightPortal.Get(i) << " ";
//    }
//    std::cout << std::endl << std::endl;


//    // ---------------------------------- TASK 1: ^ REPLACE THIS ^ ---------------------------------- //














//    // =================================== v MAIN CHANGES v ==================================== //

//    // COMMS: For each superarc, we get its low and high bounds ...
//    // COMMS: ... we then loop only from lowID to highID to
////    IdArrayType lastVertexForSuperparent;
////    auto lastVertexForSuperparentPortal = lastVertexForSuperparent.WritePortal();

////    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//    std::vector<vtkm::Id> lastVertexSuperparent(contourTree.Superarcs.GetNumberOfValues(), 0);
//    std::vector<vtkm::Id> firstVertexSuperparent(contourTree.Superarcs.GetNumberOfValues(), contourTree.Superarcs.GetNumberOfValues());

//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    { // per node in sorted order

//      // get the node's sort ID:
//      vtkm::Id sortID = nodesPortal.Get(sortedNode);
//      // ... and then get its superparent's ID:
//      vtkm::Id superparent = superparentsPortal.Get(sortID);

//      auto CTSuperarcPortal = contourTree.Superarcs.ReadPortal();

//      std::cout << sortID << " -> " << superparent << "(" << CTSuperarcPortal.Get(superparent) << ")" << std::endl;

//      if (sortID > lastVertexSuperparent[superparent])
//      {
//          lastVertexSuperparent[sortedNode] = sortID;
//      }

//      if (sortID < firstVertexSuperparent[superparent])
//      {
//          firstVertexSuperparent[sortedNode] = sortID;
//      }

////      if (sortedNode == 0)
////        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
////      else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
////        // if we have the last entry of the array that is not already the superparent in the array ...
////        // ... add a new superparent only then
////        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
//    } // per node in sorted order

//    std::cout << std::endl;

//    std::cout << "firstVertexSuperparent:" << std::endl;
//    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//    {
//        std::cout << firstVertexSuperparent[superarc] << " ";
//    }

//    std::cout << std::endl;

//    std::cout << "lastVertexSuperparent:" << std::endl;
//    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//    {
//        std::cout << lastVertexSuperparent[superarc] << " ";
//    }
//    std::cout << std::endl;


//    // THE MAIN CHANGE FOR COMPUTING THE VOLUMES HERE ...
//    // ... WE CHANGE SO THAT WE SUM OVER THE REGULAR NODES OF EACH SUPERARC
//    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//    {

//      float prefix_sum_of_regular_nodes_within_superarc = 0.f;

////      // we want to run through each regular node within the superarc
////      for (vtkm::Id regularnode = 0; )
////      {
////        prefix_sum_of_regular_nodes_within_superarc += vertex_weights[regularnode];
////      }


//      if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
//      { // COMMS: deal with the last case manually, ...
//        // COMMS: ... because there is no right nbor, use absolute last arc ID number
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          contourTree.Arcs.GetNumberOfValues() -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      }
//      else
//      { // COMMS: look at the right nbor, get the difference in sort IDs, know the count
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          firstVertexForSuperparentPortal.Get(superarc + 1) -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      }

////      superarcIntrinsicWeightPortal.Set(superarc, prefix_sum)

//    }


//    // =================================== ^ MAIN CHANGES ^ ==================================== //













//    // COMMS: STEP 3 - Dependent weights

//    // now initialise the arrays for transfer & dependent weights
//    // COMMS: (target is the top-end of an edge/arc)
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Superarcs.GetNumberOfValues()),
//      superarcDependentWeight);
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Supernodes.GetNumberOfValues()),
//      supernodeTransferWeight);
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Hyperarcs.GetNumberOfValues()),
//      hyperarcDependentWeight);

//    // set up the array which tracks which supernodes to deal with on which iteration
//    // COMMS: segment the array based on the target
//    auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
//    auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
//    auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
//    auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
//    auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

//    /*
//    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
//                          firstSupernodePerIteration);
//    auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
//    for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
//         supernode++)
//    { // per supernode
//      vtkm::Id when = MaskedIndex(whenTransferredPortal.Get(supernode));
//      if (supernode == 0)
//      { // zeroth supernode
//        firstSupernodePerIterationPortal.Set(when, supernode);
//      } // zeroth supernode
//      else if (when != MaskedIndex(whenTransferredPortal.Get(supernode - 1)))
//      { // non-matching supernode
//        firstSupernodePerIterationPortal.Set(when, supernode);
//      } // non-matching supernode
//    }   // per supernode
//    for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
//      if (firstSupernodePerIterationPortal.Get(iteration) == 0)
//        firstSupernodePerIterationPortal.Set(iteration,
//                                             firstSupernodePerIterationPortal.Get(iteration + 1));

//    // set the sentinel at the end of the array
//    firstSupernodePerIterationPortal.Set(nIterations, contourTree.Supernodes.GetNumberOfValues());

//    // now use that array to construct a similar array for hypernodes
//    IdArrayType firstHypernodePerIteration;
//    firstHypernodePerIteration.Allocate(nIterations + 1);
//    auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();
//    auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
//    auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
//    auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();
//    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
//      firstHypernodePerIterationPortal.Set(
//        iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
//    firstHypernodePerIterationPortal.Set(nIterations, contourTree.Hypernodes.GetNumberOfValues());
//    */

//    // now iterate, propagating weights inwards
//    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
//    { // per iteration
//      // pull the array bounds into register
//      vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
//      vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
//      vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
//      vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

//      // Recall that the superarcs are sorted by (iteration, hyperarc), ...
//      // ... & that all superarcs for a given hyperarc are processed in the same iteration.
//      // Assume therefore that:
//      //      i.  we now have the intrinsic weight assigned for each superarc, and
//      //      ii. we also have the transfer weight assigned for each supernode.
//      //
//      // Suppose we have a sequence of superarcs
//      //                      s11 s12 s13 s14 s21 s22 s23 s31
//      // with transfer weights at their origins and intrinsic weights along them
//      //      sArc                     s11 s12 s13 s14 s21 s22 s23 s31
//      //      transfer wt               0   1   2   1   2   3   1   0
//      //      intrinsic wt              1   2   1   5   2   6   1   1
//      //
//      //  now, if we do a prefix sum on each of these and add the two sums together, we get:
//      //      sArc                                  s11 s12 s13 s14 s21 s22 s23 s31
//      //      hyperparent sNode ID                  s11 s11 s11 s11 s21 s21 s21 s31
//      //      transfer weight                       0   1   2   1   2   3   1   0
//      //      intrinsic weight                      1   2   1   5   2   6   1   1
//      //      sum(xfer + intrinsic)                 1   3   3   6   4   9   2   1
//      //  prefix sum (xfer + int)                   1   4   7  13  17  26  28  29
//      //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  16


//      std::cout << "transfer: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << supernodeTransferWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      std::cout << "intrinsic: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << superarcIntrinsicWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      // so, step 1: add xfer + int & store in dependent weight
//      std::cout << "step 1: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//        superarcDependentWeightPortal.Set(supernode,
//                                          supernodeTransferWeightPortal.Get(supernode) +
//                                            superarcIntrinsicWeightPortal.Get(supernode));

//        std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      // step 2: perform prefix sum on the dependent weight range
//      std::cout << "step 2: ";
//      for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
//      {
//        superarcDependentWeightPortal.Set(supernode,
//                                          superarcDependentWeightPortal.Get(supernode) +
//                                            superarcDependentWeightPortal.Get(supernode - 1));

//        std::cout << superarcDependentWeightPortal.Get(supernode) + superarcDependentWeightPortal.Get(supernode - 1) << " ";
//      }
//      std::cout << std::endl;

//      // step 3: subtract out the dependent weight of the prefix to the entire hyperarc.
//      // This will be a transfer, ...
//      // ... but for now, it's easier to show it in serial.
//      // NB: Loops backwards so that computation uses the correct value
//      // As a bonus, note that we test > firstsupernode, not >=.
//      // This is because we've got unsigned integers, & otherwise it will not terminate
//      // But the first is always correct anyway
//      // (same reason as the short-cut termination on hyperparent), so we're fine
//      std::cout << "step 3: ";
//      for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
//      { // per supernode
//        // retrieve the hyperparent & convert to a supernode ID
//        vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
//        vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

//        // if the hyperparent is the first in the sequence, dependent weight is already correct
//        if (hyperparent == firstHypernode)
//          continue;

//        // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
//        superarcDependentWeightPortal.Set(
//          supernode,
//          superarcDependentWeightPortal.Get(supernode) -
//            superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

//        std::cout << superarcDependentWeightPortal.Get(supernode) - superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << " ";

//      } // per supernode
//      std::cout << std::endl;

//      // step 4: transfer the dependent weight to the hyperarc's target supernode
//      std::cout << "step 4: ";
//      for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
//      { // per hypernode
//        // last superarc for the hyperarc
//        vtkm::Id lastSuperarc;
//        // special case for the last hyperarc
//        if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
//          // take the last superarc in the array
//          lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
//        else
//          // otherwise, take the next hypernode's ID and subtract 1
//          lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

//        // now, given the last superarc for the hyperarc, transfer the dependent weight
//        hyperarcDependentWeightPortal.Set(hypernode,
//                                          superarcDependentWeightPortal.Get(lastSuperarc));

//        // note that in parallel, this will have to be split out as a sort & partial sum in another array
//        vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
//        supernodeTransferWeightPortal.Set(hyperarcTarget,
//                                          supernodeTransferWeightPortal.Get(hyperarcTarget) +
//                                            hyperarcDependentWeightPortal.Get(hypernode));
//      } // per hypernode

//      std::cout << std::endl << std::endl;

//    }   // per iteration
//  }     // ContourTreeMaker::ComputeWeights()

























//  // routine to compute the volume for each hyperarc and superarc
//  void static ComputeVolumeWeightsSerial(const ContourTree& contourTree,
//                                         const vtkm::Id nIterations,
//                                         IdArrayType& superarcIntrinsicWeight,
//                                         IdArrayType& superarcDependentWeight,
//                                         IdArrayType& supernodeTransferWeight,
//                                         IdArrayType& hyperarcDependentWeight)
//  { // ContourTreeMaker::ComputeWeights()
//    std::cout << "ComputeVolumeWeightsSerial()\n";
//    std::cout << "NUMBER OF ITERATIONS: " << nIterations << std::endl;
//    // start by storing the first sorted vertex ID for each superarc
//    IdArrayType firstVertexForSuperparent;
//    // COMMS: One vertex per superarc, so allocate the size = number of superarcs
//    firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
//    // COMMS: Each superarc will have one intrinsic weight (sum of all the regular nodes (arcs) on it
//    superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());

//    // COMMS: set up all the arrays we will need to write to:
//    auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();
//    auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();
//    auto superparentsPortal = contourTree.Superparents.ReadPortal();
//    auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
//    auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
//    auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
//    // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
//    auto nodesPortal = contourTree.Nodes.ReadPortal();

//    // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();


//    // ======================================== SWEEP HERE ========================================== //

//    std::vector<float> vertex_weights(contourTree.Arcs.GetNumberOfValues());

//    std::cout << "Sorted nodes:" << std::endl;
//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    { // per node in sorted order
//      vtkm::Id sortID = nodesPortal.Get(sortedNode);

//      std::cout << sortID << " ";

//      // HERE WE WILL GET TRIANGLE AREAS ETC.
//      // NOW JUST SET TO 1.f TO BE SAME AS THE SEGMENTED ARRAY PREFIX SUM
//      vertex_weights[sortID] = 1.f;
//    }
//    std::cout << std::endl;

//    std::cout << "weight array:" << std::endl;
//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    {
//        std::cout << vertex_weights[sortedNode] << " ";
//    }

//    // ======================================== SWEEP HERE ========================================== //


//    // ---------------------------------- TASK 1: v REPLACE THIS v ---------------------------------- //
//    // COMMS: STEP 1: Create a segmented array of superparent IDs ...
//    // COMMS: ... for each node, by its sorted ID populate the SP Segmented Array
//    // COMMS: The Segmented Superparent ID array will look like this:
//    // COMMS: ... superparent ID                [ 11, 12, 13, 14, 21, 22, 23, 31 ]
//    // COMMS: ... first node sorted ID in arc:  [  0,  1,  3,  4,  9, 11, 17, 18 ]
//    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//    { // per node in sorted order
//      vtkm::Id sortID = nodesPortal.Get(sortedNode);
//      vtkm::Id superparent = superparentsPortal.Get(sortID);
//      if (sortedNode == 0)
//        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
//      else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
//        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
//    } // per node in sorted order

////    DOES NOT WORK HERE ...
////    PrintValues(
////      "firstVertexForSuperparentPortal print:",
////      firstVertexForSuperparent,
////      -1,
////      std::cout);

//    auto firstVertexForSuperparentReadPortal = firstVertexForSuperparent.ReadPortal();

//    std::cout << std::endl << "First Vertices of Each Superarc" << std::endl;
//    for(int i = 0; i < firstVertexForSuperparentReadPortal.GetNumberOfValues(); i++)
//    {
//        std::cout << firstVertexForSuperparentReadPortal.Get(i) << " ";
//    }
//    std::cout << std::endl << std::endl;


//    // COMMS: STEP 2
//    // now we use that to compute the intrinsic weights
//    // COMMS: (The intrinsic weight = the number of regular vertices on the superarc)
//    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
//    {
//      if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
//      { // COMMS: deal with the last case manually, ...
//        // COMMS: ... because there is no right nbor, use absolute last arc ID number
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          contourTree.Arcs.GetNumberOfValues() -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      }
//      else
//      { // COMMS: look at the right nbor, get the difference in sort IDs, know the count
//        superarcIntrinsicWeightPortal.Set(superarc,
//                                          firstVertexForSuperparentPortal.Get(superarc + 1) -
//                                            firstVertexForSuperparentPortal.Get(superarc));
//      }
//    }


//    // ---------------------------------- TASK 1: ^ REPLACE THIS ^ ---------------------------------- //

//    // COMMS: STEP 3 - Dependent weights

//    // now initialise the arrays for transfer & dependent weights
//    // COMMS: (target is the top-end of an edge/arc)
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Superarcs.GetNumberOfValues()),
//      superarcDependentWeight);
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Supernodes.GetNumberOfValues()),
//      supernodeTransferWeight);
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Hyperarcs.GetNumberOfValues()),
//      hyperarcDependentWeight);

//    // set up the array which tracks which supernodes to deal with on which iteration
//    // COMMS: segment the array based on the target
//    auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
//    auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
//    auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
//    auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
//    auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

//    /*
//    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
//                          firstSupernodePerIteration);
//    auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
//    for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues();
//         supernode++)
//    { // per supernode
//      vtkm::Id when = MaskedIndex(whenTransferredPortal.Get(supernode));
//      if (supernode == 0)
//      { // zeroth supernode
//        firstSupernodePerIterationPortal.Set(when, supernode);
//      } // zeroth supernode
//      else if (when != MaskedIndex(whenTransferredPortal.Get(supernode - 1)))
//      { // non-matching supernode
//        firstSupernodePerIterationPortal.Set(when, supernode);
//      } // non-matching supernode
//    }   // per supernode
//    for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
//      if (firstSupernodePerIterationPortal.Get(iteration) == 0)
//        firstSupernodePerIterationPortal.Set(iteration,
//                                             firstSupernodePerIterationPortal.Get(iteration + 1));

//    // set the sentinel at the end of the array
//    firstSupernodePerIterationPortal.Set(nIterations, contourTree.Supernodes.GetNumberOfValues());

//    // now use that array to construct a similar array for hypernodes
//    IdArrayType firstHypernodePerIteration;
//    firstHypernodePerIteration.Allocate(nIterations + 1);
//    auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();
//    auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
//    auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
//    auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();
//    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
//      firstHypernodePerIterationPortal.Set(
//        iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
//    firstHypernodePerIterationPortal.Set(nIterations, contourTree.Hypernodes.GetNumberOfValues());
//    */

//    // now iterate, propagating weights inwards
//    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
//    { // per iteration
//      // pull the array bounds into register
//      vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
//      vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
//      vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
//      vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

//      // Recall that the superarcs are sorted by (iteration, hyperarc), ...
//      // ... & that all superarcs for a given hyperarc are processed in the same iteration.
//      // Assume therefore that:
//      //      i.  we now have the intrinsic weight assigned for each superarc, and
//      //      ii. we also have the transfer weight assigned for each supernode.
//      //
//      // Suppose we have a sequence of superarcs
//      //                      s11 s12 s13 s14 s21 s22 s23 s31
//      // with transfer weights at their origins and intrinsic weights along them
//      //      sArc                     s11 s12 s13 s14 s21 s22 s23 s31
//      //      transfer wt               0   1   2   1   2   3   1   0
//      //      intrinsic wt              1   2   1   5   2   6   1   1
//      //
//      //  now, if we do a prefix sum on each of these and add the two sums together, we get:
//      //      sArc                                  s11 s12 s13 s14 s21 s22 s23 s31
//      //      hyperparent sNode ID                  s11 s11 s11 s11 s21 s21 s21 s31
//      //      transfer weight                       0   1   2   1   2   3   1   0
//      //      intrinsic weight                      1   2   1   5   2   6   1   1
//      //      sum(xfer + intrinsic)                 1   3   3   6   4   9   2   1
//      //  prefix sum (xfer + int)                   1   4   7  13  17  26  28  29
//      //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  16


//      std::cout << "transfer: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << supernodeTransferWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      std::cout << "intrinsic: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//          std::cout << superarcIntrinsicWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      // so, step 1: add xfer + int & store in dependent weight
//      std::cout << "step 1: ";
//      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
//      {
//        superarcDependentWeightPortal.Set(supernode,
//                                          supernodeTransferWeightPortal.Get(supernode) +
//                                            superarcIntrinsicWeightPortal.Get(supernode));

//        std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
//      }
//      std::cout << std::endl;

//      // step 2: perform prefix sum on the dependent weight range
//      std::cout << "step 2: ";
//      for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
//      {
//        superarcDependentWeightPortal.Set(supernode,
//                                          superarcDependentWeightPortal.Get(supernode) +
//                                            superarcDependentWeightPortal.Get(supernode - 1));

//        std::cout << superarcDependentWeightPortal.Get(supernode) + superarcDependentWeightPortal.Get(supernode - 1) << " ";
//      }
//      std::cout << std::endl;

//      // step 3: subtract out the dependent weight of the prefix to the entire hyperarc.
//      // This will be a transfer, ...
//      // ... but for now, it's easier to show it in serial.
//      // NB: Loops backwards so that computation uses the correct value
//      // As a bonus, note that we test > firstsupernode, not >=.
//      // This is because we've got unsigned integers, & otherwise it will not terminate
//      // But the first is always correct anyway
//      // (same reason as the short-cut termination on hyperparent), so we're fine
//      std::cout << "step 3: ";
//      for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
//      { // per supernode
//        // retrieve the hyperparent & convert to a supernode ID
//        vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
//        vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

//        // if the hyperparent is the first in the sequence, dependent weight is already correct
//        if (hyperparent == firstHypernode)
//          continue;

//        // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
//        superarcDependentWeightPortal.Set(
//          supernode,
//          superarcDependentWeightPortal.Get(supernode) -
//            superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

//        std::cout << superarcDependentWeightPortal.Get(supernode) - superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << " ";

//      } // per supernode
//      std::cout << std::endl;

//      // step 4: transfer the dependent weight to the hyperarc's target supernode
//      std::cout << "step 4: ";
//      for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
//      { // per hypernode
//        // last superarc for the hyperarc
//        vtkm::Id lastSuperarc;
//        // special case for the last hyperarc
//        if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
//          // take the last superarc in the array
//          lastSuperarc = contourTree.Supernodes.GetNumberOfValues() - 1;
//        else
//          // otherwise, take the next hypernode's ID and subtract 1
//          lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

//        // now, given the last superarc for the hyperarc, transfer the dependent weight
//        hyperarcDependentWeightPortal.Set(hypernode,
//                                          superarcDependentWeightPortal.Get(lastSuperarc));

//        // note that in parallel, this will have to be split out as a sort & partial sum in another array
//        vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
//        supernodeTransferWeightPortal.Set(hyperarcTarget,
//                                          supernodeTransferWeightPortal.Get(hyperarcTarget) +
//                                            hyperarcDependentWeightPortal.Get(hypernode));
//      } // per hypernode

//      std::cout << std::endl << std::endl;

//    }   // per iteration
//  }     // ContourTreeMaker::ComputeWeights()
































    // routine to compute the branch decomposition by volume
    void static ComputeVolumeBranchDecompositionSerialFloat(const ContourTree& contourTree,
                                                       const FloatArrayType& superarcDependentWeight,
                                                       const FloatArrayType& superarcIntrinsicWeight,
                                                       IdArrayType& whichBranch,
                                                       IdArrayType& branchMinimum,
                                                       IdArrayType& branchMaximum,
                                                       IdArrayType& branchSaddle,
                                                       IdArrayType& branchParent)
    { // ComputeVolumeBranchDecomposition()
      std::cout << "ComputeVolumeBranchDecompositionSerialFloat()\n";

      // COMMS: NOTE: both intrinsic and dependent weights come pre-computed!!!
      auto superarcDependentWeightPortal = superarcDependentWeight.ReadPortal();
      auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.ReadPortal();

      // cache the number of non-root supernodes & superarcs
      vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
      vtkm::Id nSuperarcs = nSupernodes - 1;

      // STAGE I:  Find the upward and downwards weight for each superarc, and set up arrays
      // COMMS: just allocate up weight, but do not set values
      IdArrayType upWeight;
      upWeight.Allocate(nSuperarcs);
      auto upWeightPortal = upWeight.WritePortal();
      // COMMS: just allocate down weight, but do not set values
      IdArrayType downWeight;
      downWeight.Allocate(nSuperarcs);
      auto downWeightPortal = downWeight.WritePortal();
      // set up
      // initialise to a known value, indicating that no best up/down is known (yet)
      IdArrayType bestUpward;
      auto noSuchElementArray =
        vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nSupernodes);
      vtkm::cont::ArrayCopy(noSuchElementArray, bestUpward);
      IdArrayType bestDownward;
      vtkm::cont::ArrayCopy(noSuchElementArray, bestDownward);
      vtkm::cont::ArrayCopy(noSuchElementArray, whichBranch);
      auto bestUpwardPortal = bestUpward.WritePortal();
      auto bestDownwardPortal = bestDownward.WritePortal();

      // STAGE II: Pick the best (largest volume) edge upwards and downwards
      // II A. Pick the best upwards weight by sorting on lower vertex then processing by segments
      // II A 1.  Sort the superarcs by lower vertex
      // II A 2.  Per segment, best superarc writes to the best upwards array
      //          We want all of the superarcs to be listed low-end -> high-end in that order
      //          Because some are ascending and some descending, we will need to copy them carefully
      //          (the (super)arcs are oriented inwards, this means careful processing)
      vtkm::cont::ArrayHandle<EdgePair> superarcList;
      // [ask Petar why EdgePair(-1, -1) instead of NO_SUCH_ELEMENT]
      vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<EdgePair>(EdgePair(-1, -1), nSuperarcs),
                            superarcList);
      auto superarcListWritePortal = superarcList.WritePortal();
      // TODO: take the sum over the - total volume of the mesh
      vtkm::Id totalVolume = contourTree.Nodes.GetNumberOfValues();

    #ifdef DEBUG_PRINT
      std::cout << "Total Volume: " << totalVolume << std::endl;
    #endif
      // superarcs array stores the destination (supernode ID) of each superarc
      // the origin (source) of the superarc is always the supernode ID as the superarc ID
      auto superarcsPortal = contourTree.Superarcs.ReadPortal();

      // NB: Last element in array is guaranteed to be root superarc to infinity,
      // WARNING WARNING WARNING: This changes in the distributed version!
      // so we can easily skip it by not indexing to the full size
      // i.e we loop to N superarcs, not to N supernodes since N superarcs = supernodes-1
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++) //  for (vtkm::Id supernode = 0; supernode < nsupernodes-1; supernode++) vtkm::Id supernode = superarc (temporary variable)
      { // per superarc
        if (IsAscending(superarcsPortal.Get(superarc))) // flag on the ID of superarc
        { // ascending superarc
            // put the lower-end first
          superarcListWritePortal.Set(superarc, // each superarc starts at the supernode at the same ID and goes to the supernode whose ID is stored in the superarc's array
                                      // pair the origin and the destination of that superarc. We store them in an edge pair and write it to the array
                                      EdgePair(superarc, MaskedIndex(superarcsPortal.Get(superarc))));
          // because this is an ascending superarc, the dependent weight refers to the weight at the upper end
          // so, we set the up weight based on the dependent weight
          upWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          // at the inner end, dependent weight is the total in the subtree.  Then there are vertices along the edge itself (intrinsic weight), including the supernode at the outer end
          // So, to get the "dependent" weight in the other direction, we start with totalVolume - dependent, then subtract (intrinsic - 1)
          // set the weight at the down end by using the invert operator:
          downWeightPortal.Set(superarc,
                               // below is the invert operator for node count!
                               (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
                                 (superarcIntrinsicWeightPortal.Get(superarc) - 1));
        } // ascending superarc
        else
        { // descending superarc
          // lower-end is also first, but in the reverse order compared to IsAscending
          superarcListWritePortal.Set(superarc,
                                      EdgePair(MaskedIndex(superarcsPortal.Get(superarc)), superarc));
          downWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          // at the inner end, dependent weight is the total in the subtree.  Then there are vertices along the edge itself (intrinsic weight), including the supernode at the outer end
          // So, to get the "dependent" weight in the other direction, we start with totalVolume - dependent, then subtract (intrinsic - 1)
          upWeightPortal.Set(superarc,
                             (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));
        } // descending superarc
      }   // per superarc

    #ifdef DEBUG_PRINT
      std::cout << "II A. Weights Computed" << std::endl;
      PrintHeader(upWeight.GetNumberOfValues());
      //PrintIndices("Intrinsic Weight", superarcIntrinsicWeight);
      //PrintIndices("Dependent Weight", superarcDependentWeight);
      PrintIndices("Upwards Weight", upWeight);
      PrintIndices("Downwards Weight", downWeight);
      std::cout << std::endl;
    #endif

      // II B. Pick the best downwards weight by sorting on upper vertex then processing by segments
      // II B 1.      Sort the superarcs by upper vertex
      IdArrayType superarcSorter;
      superarcSorter.Allocate(nSuperarcs);
      auto superarcSorterPortal = superarcSorter.WritePortal();
      // make the array of indices for indirect sorting
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
        superarcSorterPortal.Set(superarc, superarc);

      vtkm::cont::Algorithm::Sort(
        superarcSorter,
                  // false / true = either ascending/descending
                  // sort by up/down weight so that we have a segmented array
        process_contourtree_inc_ns::SuperArcVolumetricComparator(upWeight, superarcList, false));

      // Initialize after in-place sort algorithm. (Kokkos)
      auto superarcSorterReadPortal = superarcSorter.ReadPortal();

      // II B 2.  Per segment, best superarc writes to the best upward array
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
      { // per superarc
        vtkm::Id superarcID = superarcSorterReadPortal.Get(superarc);
        const EdgePair& edge = superarcListWritePortal.Get(superarcID);
        // if it's the last one
        if (superarc == nSuperarcs - 1)
          bestDownwardPortal.Set(edge.second, edge.first);
        else
        { // not the last one
          const EdgePair& nextEdge =
            superarcListWritePortal.Get(superarcSorterReadPortal.Get(superarc + 1));
          // if the next edge belongs to another, we're the highest
          if (nextEdge.second != edge.second)
            bestDownwardPortal.Set(edge.second, edge.first);
        } // not the last one
      }   // per superarc

      // II B 3.  Repeat for lower vertex
      vtkm::cont::Algorithm::Sort(
        superarcSorter,
        process_contourtree_inc_ns::SuperArcVolumetricComparator(downWeight, superarcList, true));

      // Re-initialize after in-place sort algorithm. (Kokkos)
      superarcSorterReadPortal = superarcSorter.ReadPortal();

      // II B 2.  Per segment, best superarc writes to the best upward array
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
      { // per superarc
        vtkm::Id superarcID = superarcSorterReadPortal.Get(superarc);
        const EdgePair& edge = superarcListWritePortal.Get(superarcID);
        // if it's the last one
        if (superarc == nSuperarcs - 1)
          bestUpwardPortal.Set(edge.first, edge.second);
        else
        { // not the last one
          const EdgePair& nextEdge =
            superarcListWritePortal.Get(superarcSorterReadPortal.Get(superarc + 1));
          // if the next edge belongs to another, we're the highest
          if (nextEdge.first != edge.first)
            bestUpwardPortal.Set(edge.first, edge.second);
        } // not the last one
      }   // per superarc

    #ifdef DEBUG_PRINT
      std::cout << "II. Best Edges Selected" << std::endl;
      PrintHeader(bestUpward.GetNumberOfValues());
      PrintIndices("Best Upwards", bestUpward);
      PrintIndices("Best Downwards", bestDownward);
      std::cout << std::endl;
    #endif

      ProcessContourTree::ComputeBranchData(contourTree,
                                            whichBranch,
                                            branchMinimum,
                                            branchMaximum,
                                            branchSaddle,
                                            branchParent,
                                            bestUpward,
                                            bestDownward);
    } // ComputeVolumeBranchDecomposition()






















































  // routine to compute the branch decomposition by volume
  void static ComputeVolumeBranchDecompositionSerial(const ContourTree& contourTree,
                                                     const IdArrayType& superarcDependentWeight,
                                                     const IdArrayType& superarcIntrinsicWeight,
                                                     IdArrayType& whichBranch,
                                                     IdArrayType& branchMinimum,
                                                     IdArrayType& branchMaximum,
                                                     IdArrayType& branchSaddle,
                                                     IdArrayType& branchParent)
  { // ComputeVolumeBranchDecomposition()
    std::cout << "ComputeVolumeBranchDecompositionSerial()\n";

    // COMMS: NOTE: both intrinsic and dependent weights come pre-computed!!!
    auto superarcDependentWeightPortal = superarcDependentWeight.ReadPortal();
    auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.ReadPortal();

    // cache the number of non-root supernodes & superarcs
    vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
    vtkm::Id nSuperarcs = nSupernodes - 1;

    // STAGE I:  Find the upward and downwards weight for each superarc, and set up arrays
    // COMMS: just allocate up weight, but do not set values
    IdArrayType upWeight;
    upWeight.Allocate(nSuperarcs);
    auto upWeightPortal = upWeight.WritePortal();
    // COMMS: just allocate down weight, but do not set values
    IdArrayType downWeight;
    downWeight.Allocate(nSuperarcs);
    auto downWeightPortal = downWeight.WritePortal();
    // set up
    // initialise to a known value, indicating that no best up/down is known (yet)
    IdArrayType bestUpward;
    auto noSuchElementArray =
      vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nSupernodes);
    vtkm::cont::ArrayCopy(noSuchElementArray, bestUpward);
    IdArrayType bestDownward;
    vtkm::cont::ArrayCopy(noSuchElementArray, bestDownward);
    vtkm::cont::ArrayCopy(noSuchElementArray, whichBranch);
    auto bestUpwardPortal = bestUpward.WritePortal();
    auto bestDownwardPortal = bestDownward.WritePortal();

    // STAGE II: Pick the best (largest volume) edge upwards and downwards
    // II A. Pick the best upwards weight by sorting on lower vertex then processing by segments
    // II A 1.  Sort the superarcs by lower vertex
    // II A 2.  Per segment, best superarc writes to the best upwards array
    //          We want all of the superarcs to be listed low-end -> high-end in that order
    //          Because some are ascending and some descending, we will need to copy them carefully
    //          (the (super)arcs are oriented inwards, this means careful processing)
    vtkm::cont::ArrayHandle<EdgePair> superarcList;
    // [ask Petar why EdgePair(-1, -1) instead of NO_SUCH_ELEMENT]
    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<EdgePair>(EdgePair(-1, -1), nSuperarcs),
                          superarcList);
    auto superarcListWritePortal = superarcList.WritePortal();
    // TODO: take the sum over the - total volume of the mesh
    vtkm::Id totalVolume = contourTree.Nodes.GetNumberOfValues();

#ifdef DEBUG_PRINT
    std::cout << "Total Volume: " << totalVolume << std::endl;
#endif
    // superarcs array stores the destination (supernode ID) of each superarc
    // the origin (source) of the superarc is always the supernode ID as the superarc ID
    auto superarcsPortal = contourTree.Superarcs.ReadPortal();

    // NB: Last element in array is guaranteed to be root superarc to infinity,
    // WARNING WARNING WARNING: This changes in the distributed version!
    // so we can easily skip it by not indexing to the full size
    // i.e we loop to N superarcs, not to N supernodes since N superarcs = supernodes-1
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++) //  for (vtkm::Id supernode = 0; supernode < nsupernodes-1; supernode++) vtkm::Id supernode = superarc (temporary variable)
    { // per superarc
      if (IsAscending(superarcsPortal.Get(superarc))) // flag on the ID of superarc
      { // ascending superarc
          // put the lower-end first
        superarcListWritePortal.Set(superarc, // each superarc starts at the supernode at the same ID and goes to the supernode whose ID is stored in the superarc's array
                                    // pair the origin and the destination of that superarc. We store them in an edge pair and write it to the array
                                    EdgePair(superarc, MaskedIndex(superarcsPortal.Get(superarc))));
        // because this is an ascending superarc, the dependent weight refers to the weight at the upper end
        // so, we set the up weight based on the dependent weight
        upWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
        // at the inner end, dependent weight is the total in the subtree.  Then there are vertices along the edge itself (intrinsic weight), including the supernode at the outer end
        // So, to get the "dependent" weight in the other direction, we start with totalVolume - dependent, then subtract (intrinsic - 1)
        // set the weight at the down end by using the invert operator:
        downWeightPortal.Set(superarc,
                             // below is the invert operator for node count!
                             (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));
      } // ascending superarc
      else
      { // descending superarc
        // lower-end is also first, but in the reverse order compared to IsAscending
        superarcListWritePortal.Set(superarc,
                                    EdgePair(MaskedIndex(superarcsPortal.Get(superarc)), superarc));
        downWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
        // at the inner end, dependent weight is the total in the subtree.  Then there are vertices along the edge itself (intrinsic weight), including the supernode at the outer end
        // So, to get the "dependent" weight in the other direction, we start with totalVolume - dependent, then subtract (intrinsic - 1)
        upWeightPortal.Set(superarc,
                           (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
                             (superarcIntrinsicWeightPortal.Get(superarc) - 1));
      } // descending superarc
    }   // per superarc

#ifdef DEBUG_PRINT
    std::cout << "II A. Weights Computed" << std::endl;
    PrintHeader(upWeight.GetNumberOfValues());
    //PrintIndices("Intrinsic Weight", superarcIntrinsicWeight);
    //PrintIndices("Dependent Weight", superarcDependentWeight);
    PrintIndices("Upwards Weight", upWeight);
    PrintIndices("Downwards Weight", downWeight);
    std::cout << std::endl;
#endif

    // II B. Pick the best downwards weight by sorting on upper vertex then processing by segments
    // II B 1.      Sort the superarcs by upper vertex
    IdArrayType superarcSorter;
    superarcSorter.Allocate(nSuperarcs);
    auto superarcSorterPortal = superarcSorter.WritePortal();
    // make the array of indices for indirect sorting
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
      superarcSorterPortal.Set(superarc, superarc);

    vtkm::cont::Algorithm::Sort(
      superarcSorter,
                // false / true = either ascending/descending
                // sort by up/down weight so that we have a segmented array
      process_contourtree_inc_ns::SuperArcVolumetricComparator(upWeight, superarcList, false));

    // Initialize after in-place sort algorithm. (Kokkos)
    auto superarcSorterReadPortal = superarcSorter.ReadPortal();

    // II B 2.  Per segment, best superarc writes to the best upward array
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    { // per superarc
      vtkm::Id superarcID = superarcSorterReadPortal.Get(superarc);
      const EdgePair& edge = superarcListWritePortal.Get(superarcID);
      // if it's the last one
      if (superarc == nSuperarcs - 1)
        bestDownwardPortal.Set(edge.second, edge.first);
      else
      { // not the last one
        const EdgePair& nextEdge =
          superarcListWritePortal.Get(superarcSorterReadPortal.Get(superarc + 1));
        // if the next edge belongs to another, we're the highest
        if (nextEdge.second != edge.second)
          bestDownwardPortal.Set(edge.second, edge.first);
      } // not the last one
    }   // per superarc

    // II B 3.  Repeat for lower vertex
    vtkm::cont::Algorithm::Sort(
      superarcSorter,
      process_contourtree_inc_ns::SuperArcVolumetricComparator(downWeight, superarcList, true));

    // Re-initialize after in-place sort algorithm. (Kokkos)
    superarcSorterReadPortal = superarcSorter.ReadPortal();

    // II B 2.  Per segment, best superarc writes to the best upward array
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    { // per superarc
      vtkm::Id superarcID = superarcSorterReadPortal.Get(superarc);
      const EdgePair& edge = superarcListWritePortal.Get(superarcID);
      // if it's the last one
      if (superarc == nSuperarcs - 1)
        bestUpwardPortal.Set(edge.first, edge.second);
      else
      { // not the last one
        const EdgePair& nextEdge =
          superarcListWritePortal.Get(superarcSorterReadPortal.Get(superarc + 1));
        // if the next edge belongs to another, we're the highest
        if (nextEdge.first != edge.first)
          bestUpwardPortal.Set(edge.first, edge.second);
      } // not the last one
    }   // per superarc

#ifdef DEBUG_PRINT
    std::cout << "II. Best Edges Selected" << std::endl;
    PrintHeader(bestUpward.GetNumberOfValues());
    PrintIndices("Best Upwards", bestUpward);
    PrintIndices("Best Downwards", bestDownward);
    std::cout << std::endl;
#endif

    ProcessContourTree::ComputeBranchData(contourTree,
                                          whichBranch,
                                          branchMinimum,
                                          branchMaximum,
                                          branchSaddle,
                                          branchParent,
                                          bestUpward,
                                          bestDownward);
  } // ComputeVolumeBranchDecomposition()

  // routine to compute the branch decomposition by volume
  void static ComputeBranchData(const ContourTree& contourTree,
                                IdArrayType& whichBranch,
                                IdArrayType& branchMinimum,
                                IdArrayType& branchMaximum,
                                IdArrayType& branchSaddle,
                                IdArrayType& branchParent,
                                IdArrayType& bestUpward,
                                IdArrayType& bestDownward)
  { // ComputeBranchData()

    // Each superarc has an up and a down supernode.
    // Each supernode has a best up and a best down.
    // Suppose we have supernode N, and its best up is superarc U.
    // Then, the down supernode of U must be N
    // However, for a superarc A incident to N from above - that is not the best up from N - this is not true
    // i.e any superarc A for which bestUp[down[A]] != A is at the end of the branch, as is any leaf (with no bestUP/bestDown)
    // (we have now got a test to where the branches stop and we loop those arcs to point to themselves or to the supernode at that end)
    //

    // Set up constants
    vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
    // initially, none of them belong to a branch:
    auto noSuchElementArray =
      vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nSupernodes);
    vtkm::cont::ArrayCopy(noSuchElementArray, whichBranch);

    // STAGE III: For each vertex, identify which neighbours are on same branch
    // Let v = BestUp(u). Then if u = BestDown(v), copy BestUp(u) to whichBranch(u)
    // Otherwise, let whichBranch(u) = BestUp(u) | TERMINAL to mark the end of the side branch
    // NB 1: Leaves already have the flag set, but it's redundant so its safe
    // NB 2: We don't need to do it downwards because it's symmetric
    vtkm::cont::Invoker invoke;
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::PropagateBestUpDown
      propagateBestUpDownWorklet;
    invoke(propagateBestUpDownWorklet, bestUpward, bestDownward, whichBranch);

#ifdef DEBUG_PRINT
    std::cout << "III. Branch Neighbours Identified" << std::endl;
    PrintHeader(whichBranch.GetNumberOfValues());
    PrintIndices("Which Branch", whichBranch);
    std::cout << std::endl;
#endif

    // STAGE IV: Use pointer-doubling on whichBranch to propagate branches
    // Compute the number of log steps required in this pass
    vtkm::Id numLogSteps = 1;
    for (vtkm::Id shifter = nSupernodes; shifter != 0; shifter >>= 1)
      numLogSteps++;

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::PointerDoubling pointerDoubling(
      nSupernodes);

    // use pointer-doubling to build the branches
    for (vtkm::Id iteration = 0; iteration < numLogSteps; iteration++)
    { // per iteration
      invoke(pointerDoubling, whichBranch);
    } // per iteration


#ifdef DEBUG_PRINT
    std::cout << "IV. Branch Chains Propagated" << std::endl;
    PrintHeader(whichBranch.GetNumberOfValues());
    PrintIndices("Which Branch", whichBranch);
    std::cout << std::endl;
#endif

    // Initialise
    IdArrayType chainToBranch;
    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nSupernodes), chainToBranch);

    // Set 1 to every relevant
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::PrepareChainToBranch
      prepareChainToBranchWorklet;
    invoke(prepareChainToBranchWorklet, whichBranch, chainToBranch);

    // Prefix scanto get IDs
    vtkm::Id nBranches = vtkm::cont::Algorithm::ScanInclusive(chainToBranch, chainToBranch);

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::FinaliseChainToBranch
      finaliseChainToBranchWorklet;
    invoke(finaliseChainToBranchWorklet, whichBranch, chainToBranch);

    // V B.  Create the arrays for the branches
    auto noSuchElementArrayNBranches =
      vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nBranches);
    vtkm::cont::ArrayCopy(noSuchElementArrayNBranches, branchMinimum);
    vtkm::cont::ArrayCopy(noSuchElementArrayNBranches, branchMaximum);
    vtkm::cont::ArrayCopy(noSuchElementArrayNBranches, branchSaddle);
    vtkm::cont::ArrayCopy(noSuchElementArrayNBranches, branchParent);

#ifdef DEBUG_PRINT
    std::cout << "V. Branch Arrays Created" << std::endl;
    PrintHeader(chainToBranch.GetNumberOfValues());
    PrintIndices("Chain To Branch", chainToBranch);
    PrintHeader(nBranches);
    PrintIndices("Branch Minimum", branchMinimum);
    PrintIndices("Branch Maximum", branchMaximum);
    PrintIndices("Branch Saddle", branchSaddle);
    PrintIndices("Branch Parent", branchParent);
#endif

    IdArrayType supernodeSorter;
    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleIndex(nSupernodes), supernodeSorter);

    vtkm::cont::Algorithm::Sort(
      supernodeSorter,
      process_contourtree_inc_ns::SuperNodeBranchComparator(whichBranch, contourTree.Supernodes));

    IdArrayType permutedBranches;
    permutedBranches.Allocate(nSupernodes);
    PermuteArray<vtkm::Id>(whichBranch, supernodeSorter, permutedBranches);

    IdArrayType permutedRegularID;
    permutedRegularID.Allocate(nSupernodes);
    PermuteArray<vtkm::Id>(contourTree.Supernodes, supernodeSorter, permutedRegularID);

#ifdef DEBUG_PRINT
    std::cout << "VI A. Sorted into Branches" << std::endl;
    PrintHeader(nSupernodes);
    PrintIndices("Supernode IDs", supernodeSorter);
    PrintIndices("Branch", permutedBranches);
    PrintIndices("Regular ID", permutedRegularID);
#endif

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::WhichBranchNewId
      whichBranchNewIdWorklet;
    invoke(whichBranchNewIdWorklet, chainToBranch, whichBranch);

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::BranchMinMaxSet
      branchMinMaxSetWorklet(nSupernodes);
    invoke(branchMinMaxSetWorklet, supernodeSorter, whichBranch, branchMinimum, branchMaximum);

#ifdef DEBUG_PRINT
    std::cout << "VI. Branches Set" << std::endl;
    PrintHeader(nBranches);
    PrintIndices("Branch Maximum", branchMaximum);
    PrintIndices("Branch Minimum", branchMinimum);
    PrintIndices("Branch Saddle", branchSaddle);
    PrintIndices("Branch Parent", branchParent);
#endif

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::BranchSaddleParentSet
      branchSaddleParentSetWorklet;
    invoke(branchSaddleParentSetWorklet,
           whichBranch,
           branchMinimum,
           branchMaximum,
           bestDownward,
           bestUpward,
           branchSaddle,
           branchParent);

#ifdef DEBUG_PRINT
    std::cout << "VII. Branches Constructed" << std::endl;
    PrintHeader(nBranches);
    PrintIndices("Branch Maximum", branchMaximum);
    PrintIndices("Branch Minimum", branchMinimum);
    PrintIndices("Branch Saddle", branchSaddle);
    PrintIndices("Branch Parent", branchParent);
#endif

  } // ComputeBranchData()

  // Create branch decomposition from contour tree
  template <typename T, typename StorageType>
  static process_contourtree_inc_ns::Branch<T>* ComputeBranchDecomposition(
    const IdArrayType& contourTreeSuperparents,
    const IdArrayType& contourTreeSupernodes,
    const IdArrayType& whichBranch,
    const IdArrayType& branchMinimum,
    const IdArrayType& branchMaximum,
    const IdArrayType& branchSaddle,
    const IdArrayType& branchParent,
    const IdArrayType& sortOrder,
    const vtkm::cont::ArrayHandle<T, StorageType>& dataField,
    bool dataFieldIsSorted)
  {
    return process_contourtree_inc_ns::Branch<T>::ComputeBranchDecomposition(
      contourTreeSuperparents,
      contourTreeSupernodes,
      whichBranch,
      branchMinimum,
      branchMaximum,
      branchSaddle,
      branchParent,
      sortOrder,
      dataField,
      dataFieldIsSorted);
  }



  void static ComputeVolumeBranchDecomposition(const ContourTree& contourTree,
                                               const vtkm::Id nIterations,
                                               IdArrayType& whichBranch,
                                               IdArrayType& branchMinimum,
                                               IdArrayType& branchMaximum,
                                               IdArrayType& branchSaddle,
                                               IdArrayType& branchParent)
  // modified version of the volume computation ...
  { // ComputeHeightBranchDecomposition()

    std::cout << "ComputeVolumeBranchDecomposition() - PLL???\n";

    vtkm::cont::Invoker Invoke;

    // STEP 1. Compute the number of nodes in every superarc, that's the intrinsic weight
    // MODIFICATION FOR DELAUNAY MESHES:#
    // ... remove segmented array for computing intrinsic weights ...
    // ... replace with prefix sum over the regular vertices
    IdArrayType superarcIntrinsicWeight;
    superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());

    // for vertex v:
    //  for each triangle T:
    //      compute A(T) / 3
    //      add to weight wt(V)



    IdArrayType firstVertexForSuperparent;
    firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());

    // Compute the number of regular nodes on every superarcs (the intrinsic weight)
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetFirstVertexForSuperparent
      setFirstVertexForSuperparent;
    Invoke(setFirstVertexForSuperparent,
           contourTree.Nodes,
           contourTree.Superparents,
           firstVertexForSuperparent);

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::ComputeIntrinsicWeight
      computeIntrinsicWeight;
    Invoke(computeIntrinsicWeight,
           contourTree.Arcs,
           contourTree.Superarcs,
           firstVertexForSuperparent,
           superarcIntrinsicWeight);


    // Cache the number of non-root supernodes & superarcs
    vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
    auto noSuchElementArray =
      vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nSupernodes);

    // Set up bestUpward and bestDownward array, these are the things we want to compute in this routine.
    IdArrayType bestUpward, bestDownward;
    vtkm::cont::ArrayCopy(noSuchElementArray, bestUpward);
    vtkm::cont::ArrayCopy(noSuchElementArray, bestDownward);

    // We initiale with the weight of the superarcs, once we sum those up we'll get the hypersweep weight
    IdArrayType sumValues;
    vtkm::cont::ArrayCopy(superarcIntrinsicWeight, sumValues);

    // This should be 0 here, because we're not changing the root
    vtkm::cont::ArrayHandle<vtkm::Id> howManyUsed;
    vtkm::cont::ArrayCopy(
      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Hyperarcs.GetNumberOfValues()),
      howManyUsed);

    // Perform a sum hypersweep
    hyperarcScan<decltype(vtkm::Sum())>(contourTree.Supernodes,
                                        contourTree.Hypernodes,
                                        contourTree.Hyperarcs,
                                        contourTree.Hyperparents,
                                        contourTree.Hyperparents,
                                        contourTree.WhenTransferred,
                                        howManyUsed,
                                        nIterations,
                                        vtkm::Sum(),
                                        sumValues);

    // For every directed arc store the volume of it's associate subtree
    vtkm::cont::ArrayHandle<vtkm::worklet::contourtree_augmented::EdgeDataVolume> arcs;
    arcs.Allocate(contourTree.Superarcs.GetNumberOfValues() * 2 - 2);

    vtkm::Id totalVolume = contourTree.Nodes.GetNumberOfValues();
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::InitialiseArcsVolume initArcs(
      totalVolume);
    Invoke(initArcs, sumValues, superarcIntrinsicWeight, contourTree.Superarcs, arcs);

    // Sort arcs to obtain the best up and down
    vtkm::cont::Algorithm::Sort(arcs, vtkm::SortLess());

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetBestUpDown setBestUpDown;
    Invoke(setBestUpDown, bestUpward, bestDownward, arcs);

    ProcessContourTree::ComputeBranchData(contourTree,
                                          whichBranch,
                                          branchMinimum,
                                          branchMaximum,
                                          branchSaddle,
                                          branchParent,
                                          bestUpward,
                                          bestDownward);

  } // ComputeHeightBranchDecomposition()




  // routine to compute the branch decomposition by volume
//  void static ComputeVolumeBranchDecomposition(const ContourTree& contourTree,
//                                               const vtkm::Id nIterations,
//                                               IdArrayType& whichBranch,
//                                               IdArrayType& branchMinimum,
//                                               IdArrayType& branchMaximum,
//                                               IdArrayType& branchSaddle,
//                                               IdArrayType& branchParent)
//  { // ComputeHeightBranchDecomposition()

//    vtkm::cont::Invoker Invoke;

//    // STEP 1. Compute the number of nodes in every superarc, that's the intrinsic weight
//    IdArrayType superarcIntrinsicWeight;
//    superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());

//    IdArrayType firstVertexForSuperparent;
//    firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());

//    // Compute the number of regular nodes on every superarcs (the intrinsic weight)
//    vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetFirstVertexForSuperparent
//      setFirstVertexForSuperparent;
//    Invoke(setFirstVertexForSuperparent,
//           contourTree.Nodes,
//           contourTree.Superparents,
//           firstVertexForSuperparent);

//    vtkm::worklet::contourtree_augmented::process_contourtree_inc::ComputeIntrinsicWeight
//      computeIntrinsicWeight;
//    Invoke(computeIntrinsicWeight,
//           contourTree.Arcs,
//           contourTree.Superarcs,
//           firstVertexForSuperparent,
//           superarcIntrinsicWeight);


//    // Cache the number of non-root supernodes & superarcs
//    vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
//    auto noSuchElementArray =
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nSupernodes);

//    // Set up bestUpward and bestDownward array, these are the things we want to compute in this routine.
//    IdArrayType bestUpward, bestDownward;
//    vtkm::cont::ArrayCopy(noSuchElementArray, bestUpward);
//    vtkm::cont::ArrayCopy(noSuchElementArray, bestDownward);

//    // We initiale with the weight of the superarcs, once we sum those up we'll get the hypersweep weight
//    IdArrayType sumValues;
//    vtkm::cont::ArrayCopy(superarcIntrinsicWeight, sumValues);

//    // This should be 0 here, because we're not changing the root
//    vtkm::cont::ArrayHandle<vtkm::Id> howManyUsed;
//    vtkm::cont::ArrayCopy(
//      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, contourTree.Hyperarcs.GetNumberOfValues()),
//      howManyUsed);

//    // Perform a sum hypersweep
//    hyperarcScan<decltype(vtkm::Sum())>(contourTree.Supernodes,
//                                        contourTree.Hypernodes,
//                                        contourTree.Hyperarcs,
//                                        contourTree.Hyperparents,
//                                        contourTree.Hyperparents,
//                                        contourTree.WhenTransferred,
//                                        howManyUsed,
//                                        nIterations,
//                                        vtkm::Sum(),
//                                        sumValues);

//    // For every directed arc store the volume of it's associate subtree
//    vtkm::cont::ArrayHandle<vtkm::worklet::contourtree_augmented::EdgeDataVolume> arcs;
//    arcs.Allocate(contourTree.Superarcs.GetNumberOfValues() * 2 - 2);

//    vtkm::Id totalVolume = contourTree.Nodes.GetNumberOfValues();
//    vtkm::worklet::contourtree_augmented::process_contourtree_inc::InitialiseArcsVolume initArcs(
//      totalVolume);
//    Invoke(initArcs, sumValues, superarcIntrinsicWeight, contourTree.Superarcs, arcs);

//    // Sort arcs to obtain the best up and down
//    vtkm::cont::Algorithm::Sort(arcs, vtkm::SortLess());

//    vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetBestUpDown setBestUpDown;
//    Invoke(setBestUpDown, bestUpward, bestDownward, arcs);

//    ProcessContourTree::ComputeBranchData(contourTree,
//                                          whichBranch,
//                                          branchMinimum,
//                                          branchMaximum,
//                                          branchSaddle,
//                                          branchParent,
//                                          bestUpward,
//                                          bestDownward);

//  } // ComputeHeightBranchDecomposition()

  // routine to compute the branch decomposition by volume
  void static ComputeHeightBranchDecomposition(const ContourTree& contourTree,
                                               const cont::ArrayHandle<Float64> fieldValues,
                                               const IdArrayType& ctSortOrder,
                                               const vtkm::Id nIterations,
                                               IdArrayType& whichBranch,
                                               IdArrayType& branchMinimum,
                                               IdArrayType& branchMaximum,
                                               IdArrayType& branchSaddle,
                                               IdArrayType& branchParent)
  { // ComputeHeightBranchDecomposition()

    std::cout << "ComputeHeightBranchDecomposition() - PLL???\n";

    // Cache the number of non-root supernodes & superarcs
    vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
    auto noSuchElementArray =
      vtkm::cont::ArrayHandleConstant<vtkm::Id>((vtkm::Id)NO_SUCH_ELEMENT, nSupernodes);

    // Set up bestUpward and bestDownward array, these are the things we want to compute in this routine.
    IdArrayType bestUpward, bestDownward;
    vtkm::cont::ArrayCopy(noSuchElementArray, bestUpward);
    vtkm::cont::ArrayCopy(noSuchElementArray, bestDownward);

    // maxValues and minValues store the values from the max and min hypersweep respectively.
    IdArrayType minValues, maxValues;
    vtkm::cont::ArrayCopy(contourTree.Supernodes, maxValues);
    vtkm::cont::ArrayCopy(contourTree.Supernodes, minValues);

    // Store the direction of the superarcs in the min and max hypersweep (find a way to get rid of these, the only differing direction is on the path from the root to the min/max).
    IdArrayType minParents, maxParents;
    vtkm::cont::ArrayCopy(contourTree.Superarcs, minParents);
    vtkm::cont::ArrayCopy(contourTree.Superarcs, maxParents);

    auto minParentsPortal = minParents.WritePortal();
    auto maxParentsPortal = maxParents.WritePortal();

    // Cache the glonal minimum and global maximum (these will be the roots in the min and max hypersweep)
    Id minSuperNode = MaskedIndex(contourTree.Superparents.ReadPortal().Get(0));
    Id maxSuperNode = MaskedIndex(
      contourTree.Superparents.ReadPortal().Get(contourTree.Nodes.GetNumberOfValues() - 1));

    // Find the path from the global minimum to the root, not parallelisable (but it's fast, no need to parallelise)
    auto minPath = findSuperPathToRoot(contourTree.Superarcs.ReadPortal(), minSuperNode);

    // Find the path from the global minimum to the root, not parallelisable (but it's fast, no need to parallelise)
    auto maxPath = findSuperPathToRoot(contourTree.Superarcs.ReadPortal(), maxSuperNode);

    // Reserve the direction of the superarcs on the min path.
    for (std::size_t i = 1; i < minPath.size(); i++)
    {
      minParentsPortal.Set(minPath[i], minPath[i - 1]);
    }
    minParentsPortal.Set(minPath[0], 0);

    // Reserve the direction of the superarcs on the max path.
    for (std::size_t i = 1; i < maxPath.size(); i++)
    {
      maxParentsPortal.Set(maxPath[i], maxPath[i - 1]);
    }
    maxParentsPortal.Set(maxPath[0], 0);

    vtkm::cont::Invoker Invoke;
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::UnmaskArray unmaskArrayWorklet;
    Invoke(unmaskArrayWorklet, minValues);
    Invoke(unmaskArrayWorklet, maxValues);

    // Thse arrays hold the changes hyperarcs in the min and max hypersweep respectively
    vtkm::cont::ArrayHandle<vtkm::Id> minHyperarcs, maxHyperarcs;
    vtkm::cont::ArrayCopy(contourTree.Hyperarcs, minHyperarcs);
    vtkm::cont::ArrayCopy(contourTree.Hyperarcs, maxHyperarcs);

    // These arrays hold the changed hyperarcs for the min and max hypersweep
    vtkm::cont::ArrayHandle<vtkm::Id> minHyperparents, maxHyperparents;
    vtkm::cont::ArrayCopy(contourTree.Hyperparents, minHyperparents);
    vtkm::cont::ArrayCopy(contourTree.Hyperparents, maxHyperparents);

    auto minHyperparentsPortal = minHyperparents.WritePortal();
    auto maxHyperparentsPortal = maxHyperparents.WritePortal();

    for (std::size_t i = 0; i < minPath.size(); i++)
    {
      // Set a unique dummy Id (something that the prefix scan by key will leave alone)
      minHyperparentsPortal.Set(minPath[i],
                                contourTree.Hypernodes.GetNumberOfValues() + minPath[i]);
    }

    for (std::size_t i = 0; i < maxPath.size(); i++)
    {
      // Set a unique dummy Id (something that the prefix scan by key will leave alone)
      maxHyperparentsPortal.Set(maxPath[i],
                                contourTree.Hypernodes.GetNumberOfValues() + maxPath[i]);
    }

    // These arrays hold the number of nodes in each hypearcs that are on the min or max path for the min and max hypersweep respectively.
    vtkm::cont::ArrayHandle<vtkm::Id> minHowManyUsed, maxHowManyUsed;
    vtkm::cont::ArrayCopy(
      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, maxHyperarcs.GetNumberOfValues()),
      minHowManyUsed);
    vtkm::cont::ArrayCopy(
      vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, maxHyperarcs.GetNumberOfValues()),
      maxHowManyUsed);

    // Min Hypersweep
    const auto minOperator = vtkm::Minimum();

    // Cut hyperarcs at the first node on the path from the max to the root
    editHyperarcs(contourTree.Hyperparents.ReadPortal(),
                  minPath,
                  minHyperarcs.WritePortal(),
                  minHowManyUsed.WritePortal());

    // Perform an ordinary hypersweep on those new hyperarcs
    hyperarcScan<decltype(vtkm::Minimum())>(contourTree.Supernodes,
                                            contourTree.Hypernodes,
                                            minHyperarcs,
                                            contourTree.Hyperparents,
                                            minHyperparents,
                                            contourTree.WhenTransferred,
                                            minHowManyUsed,
                                            nIterations,
                                            vtkm::Minimum(),
                                            minValues);

    // Prefix sum along the path from the min to the root
    fixPath(vtkm::Minimum(), minPath, minValues.WritePortal());

    // Max Hypersweep
    const auto maxOperator = vtkm::Maximum();

    // Cut hyperarcs at the first node on the path from the max to the root
    editHyperarcs(contourTree.Hyperparents.ReadPortal(),
                  maxPath,
                  maxHyperarcs.WritePortal(),
                  maxHowManyUsed.WritePortal());

    // Perform an ordinary hypersweep on those new hyperarcs
    hyperarcScan<decltype(vtkm::Maximum())>(contourTree.Supernodes,
                                            contourTree.Hypernodes,
                                            maxHyperarcs,
                                            contourTree.Hyperparents,
                                            maxHyperparents,
                                            contourTree.WhenTransferred,
                                            maxHowManyUsed,
                                            nIterations,
                                            vtkm::Maximum(),
                                            maxValues);

    // Prefix sum along the path from the max to the root
    fixPath(vtkm::Maximum(), maxPath, maxValues.WritePortal());

    // For every directed edge (a, b) consider that subtree who's root is b and does not contain a.
    // We have so far found the min and max in all sub subtrees, now we compare those to a and incorporate a into that.
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::IncorporateParent<decltype(
      vtkm::Minimum())>
      incorporateParentMinimumWorklet(minOperator);
    Invoke(incorporateParentMinimumWorklet, minParents, contourTree.Supernodes, minValues);

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::IncorporateParent<decltype(
      vtkm::Maximum())>
      incorporateParentMaximumWorklet(maxOperator);
    Invoke(incorporateParentMaximumWorklet, maxParents, contourTree.Supernodes, maxValues);

    // Initialise all directed superarcs in the contour tree. Those will correspond to subtrees whos height we need for the branch decomposition.
    vtkm::cont::ArrayHandle<vtkm::worklet::contourtree_augmented::EdgeDataHeight> arcs;
    arcs.Allocate(contourTree.Superarcs.GetNumberOfValues() * 2 - 2);

    vtkm::worklet::contourtree_augmented::process_contourtree_inc::InitialiseArcs initArcs(
      0, contourTree.Arcs.GetNumberOfValues() - 1, minPath[minPath.size() - 1]);

    Invoke(initArcs, minParents, maxParents, minValues, maxValues, contourTree.Superarcs, arcs);

    // Use the min & max to compute the height of all subtrees
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::ComputeSubtreeHeight
      computeSubtreeHeight;
    Invoke(computeSubtreeHeight, fieldValues, ctSortOrder, contourTree.Supernodes, arcs);

    // Sort all directed edges based on the height of their subtree
    vtkm::cont::Algorithm::Sort(arcs, vtkm::SortLess());

    // Select a best up and best down neighbour for every vertex in the contour tree using heights of all subtrees
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetBestUpDown setBestUpDown;
    Invoke(setBestUpDown, bestUpward, bestDownward, arcs);

    // Having computed the bestUp/Down we can propagte those to obtain the branches of the branch decomposition
    ProcessContourTree::ComputeBranchData(contourTree,
                                          whichBranch,
                                          branchMinimum,
                                          branchMaximum,
                                          branchSaddle,
                                          branchParent,
                                          bestUpward,
                                          bestDownward);

  } // ComputeHeightBranchDecomposition()

  std::vector<Id> static findSuperPathToRoot(
    vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType parentsPortal,
    vtkm::Id vertex)
  {
    // Initialise the empty path and starting vertex
    std::vector<vtkm::Id> path;
    vtkm::Id current = vertex;

    // Go up the parent list until we reach the root
    while (MaskedIndex(parentsPortal.Get(current)) != 0)
    {
      path.push_back(current);
      current = MaskedIndex(parentsPortal.Get(current));
    }
    path.push_back(current);

    return path;
  }

  // Given a path from a leaf (the global min/max) to the root of the contour tree and a hypersweep (where all hyperarcs are cut at the path)
  // This function performs a prefix scan along that path to obtain the correct hypersweep values (as in the global min/max is the root of the hypersweep)
  void static fixPath(const std::function<vtkm::Id(vtkm::Id, vtkm::Id)> operation,
                      const std::vector<vtkm::Id> path,
                      vtkm::cont::ArrayHandle<vtkm::Id>::WritePortalType minMaxIndex)
  {
    using vtkm::worklet::contourtree_augmented::MaskedIndex;

    // Fix path from the old root to the new root. Parallelisble with a prefix scan, but sufficiently fast for now.
    for (auto i = path.size() - 2; i > 0; i--)
    {
      const auto vertex = path[i + 1];
      const auto parent = path[i];

      const auto vertexValue = minMaxIndex.Get(vertex);
      const auto parentValue = minMaxIndex.Get(parent);

      minMaxIndex.Set(parent, operation(vertexValue, parentValue));
    }
  }

  // This function edits all the hyperarcs which contain vertices which are on the supplied path. This path is usually the path between the global min/max to the root of the tree.
  // This function effectively cuts hyperarcs at the first node they encounter along that path.
  // In addition to this it computed the number of supernodes every hyperarc has on that path. This helps in the function hyperarcScan for choosing the new target of the cut hyperarcs.
  // NOTE: It is assumed that the supplied path starts at a leaf and ends at the root of the contour tree.
  void static editHyperarcs(
    const vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType hyperparentsPortal,
    const std::vector<vtkm::Id> path,
    vtkm::cont::ArrayHandle<vtkm::Id>::WritePortalType hyperarcsPortal,
    vtkm::cont::ArrayHandle<vtkm::Id>::WritePortalType howManyUsedPortal)
  {
    using vtkm::worklet::contourtree_augmented::MaskedIndex;

    std::size_t i = 0;
    while (i < path.size())
    {
      // Cut the hyperacs at the first point
      hyperarcsPortal.Set(MaskedIndex(hyperparentsPortal.Get(path[i])), path[i]);

      Id currentHyperparent = MaskedIndex(hyperparentsPortal.Get(path[i]));

      // Skip the rest of the supernodes which are on the same hyperarc
      while (i < path.size() && MaskedIndex(hyperparentsPortal.Get(path[i])) == currentHyperparent)
      {
        const auto value = howManyUsedPortal.Get(MaskedIndex(hyperparentsPortal.Get(path[i])));
        howManyUsedPortal.Set(MaskedIndex(hyperparentsPortal.Get(path[i])), value + 1);
        i++;
      }
    }
  }

  template <class BinaryFunctor>
  void static hyperarcScan(const vtkm::cont::ArrayHandle<vtkm::Id> supernodes,
                           const vtkm::cont::ArrayHandle<vtkm::Id> hypernodes,
                           const vtkm::cont::ArrayHandle<vtkm::Id> hyperarcs,
                           const vtkm::cont::ArrayHandle<vtkm::Id> hyperparents,
                           const vtkm::cont::ArrayHandle<vtkm::Id> hyperparentKeys,
                           const vtkm::cont::ArrayHandle<vtkm::Id> whenTransferred,
                           const vtkm::cont::ArrayHandle<vtkm::Id> howManyUsed,
                           const vtkm::Id nIterations,
                           const BinaryFunctor operation,
                           vtkm::cont::ArrayHandle<vtkm::Id> minMaxIndex)
  {
    using vtkm::worklet::contourtree_augmented::MaskedIndex;

    vtkm::cont::Invoker invoke;

    auto supernodesPortal = supernodes.ReadPortal();
    auto hypernodesPortal = hypernodes.ReadPortal();
    auto hyperparentsPortal = hyperparents.ReadPortal();

    // Set the first supernode per iteration
    vtkm::cont::ArrayHandle<vtkm::Id> firstSupernodePerIteration;
    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
                          firstSupernodePerIteration);

    // The first different from the previous is the first in the iteration
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetFirstSupernodePerIteration
      setFirstSupernodePerIteration;
    invoke(setFirstSupernodePerIteration, whenTransferred, firstSupernodePerIteration);

    auto firstSupernodePerIterationPortal = firstSupernodePerIteration.WritePortal();
    for (vtkm::Id iteration = 1; iteration < nIterations; ++iteration)
    {
      if (firstSupernodePerIterationPortal.Get(iteration) == 0)
      {
        firstSupernodePerIterationPortal.Set(iteration,
                                             firstSupernodePerIterationPortal.Get(iteration + 1));
      }
    }

    // set the sentinel at the end of the array
    firstSupernodePerIterationPortal.Set(nIterations, supernodesPortal.GetNumberOfValues());

    // Set the first hypernode per iteration
    vtkm::cont::ArrayHandle<vtkm::Id> firstHypernodePerIteration;
    vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleConstant<vtkm::Id>(0, nIterations + 1),
                          firstHypernodePerIteration);
    auto firstHypernodePerIterationPortal = firstHypernodePerIteration.WritePortal();

    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
    {
      firstHypernodePerIterationPortal.Set(
        iteration, hyperparentsPortal.Get(firstSupernodePerIterationPortal.Get(iteration)));
    }

    // Set the sentinel at the end of the array
    firstHypernodePerIterationPortal.Set(nIterations, hypernodesPortal.GetNumberOfValues());

    // This workled is used in every iteration of the following loop, so it's initialised outside.
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::AddDependentWeightHypersweep<
      BinaryFunctor>
      addDependentWeightHypersweepWorklet(operation);

    // For every iteration do a prefix scan on all hyperarcs in that iteration and then transfer the scanned value to every hyperarc's target supernode
    for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
    {
      // Determine the first and last hypernode in the current iteration (all hypernodes between them are also in the current iteration)
      vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
      vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);
      lastHypernode = vtkm::Minimum()(lastHypernode, hypernodes.GetNumberOfValues() - 1);

      // Determine the first and last supernode in the current iteration (all supernode between them are also in the current iteration)
      vtkm::Id firstSupernode = MaskedIndex(hypernodesPortal.Get(firstHypernode));
      vtkm::Id lastSupernode = MaskedIndex(hypernodesPortal.Get(lastHypernode));
      lastSupernode = vtkm::Minimum()(lastSupernode, hyperparents.GetNumberOfValues() - 1);

      // Prefix scan along all hyperarcs in the current iteration
      auto subarrayValues = vtkm::cont::make_ArrayHandleView(
        minMaxIndex, firstSupernode, lastSupernode - firstSupernode);
      auto subarrayKeys = vtkm::cont::make_ArrayHandleView(
        hyperparentKeys, firstSupernode, lastSupernode - firstSupernode);
      vtkm::cont::Algorithm::ScanInclusiveByKey(
        subarrayKeys, subarrayValues, subarrayValues, operation);

      // Array containing the Ids of the hyperarcs in the current iteration
      vtkm::cont::ArrayHandleCounting<vtkm::Id> iterationHyperarcs(
        firstHypernode, 1, lastHypernode - firstHypernode);

      // Transfer the value accumulated in the last entry of the prefix scan to the hypernode's targe supernode
      invoke(addDependentWeightHypersweepWorklet,
             iterationHyperarcs,
             hypernodes,
             hyperarcs,
             howManyUsed,
             minMaxIndex);
    }
  }
}; // class ProcessContourTree
} // namespace contourtree_augmented
} // worklet
} // vtkm

#endif
