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
#include <vtkm/Math.h>
#include <vtkm/VectorAnalysis.h>

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


struct Tetrahedron
{
    int p1;
    int p2;
    int p3;
    int p4;
};


struct Coefficients
{
    double h1;
    double h2;
    double h3;
    double h4;
};

// since in our code vectors are associated with points in space ...
// ... we wrap vtkm vectors, which are directional, to a position vector ...
// ... by keeping track of start and end points
class PositionVector
{
public:
    vtkm::Vec3f_32 start;
    vtkm::Vec3f_32 end;
    // we define difference as end-start
    vtkm::Vec3f_32 difference;
    // direction is just the normalised (unit) difference between start and end
    vtkm::Vec3f_32 direction;


    PositionVector(vtkm::Vec3f_32 aStart, vtkm::Vec3f_32 aEnd)
    {
        start = aStart;
        end = aEnd;
        difference = aEnd - aStart;
        direction = vtkm::Normal(difference);
    }

    void lerp(double interpolant)
    {
        // start does not change, only update the end and difference
        difference = difference * interpolant;
        end = start + difference;

        direction = vtkm::Normal(difference);

    }

    vtkm::Vec3f_32 lerp2point(double interpolant)
    {
        return start + (difference * interpolant);
    }

    double mag()
    {
        return vtkm::Magnitude(difference);
    }

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



std::vector<Tetrahedron> ReadTetsFromFile(const std::string& filename) {
//    std::map<vtkm::Id, Coordinates> coordinatesMap;
    std::vector<Tetrahedron> tetMap;
    std::ifstream file(filename);
    std::string line;

    if (!file)
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return tetMap; // Return an empty map if file opening fails.
    }

    std::cout << "READING ... \n";

    while (getline(file, line))
    {
        std::istringstream iss(line);
        vtkm::Id id;
        Tetrahedron tet;

        // Extracting ID and coordinates from the current line, assuming the ID is a vtkm::Id.
        if (!(iss >> id >> tet.p1 >> tet.p2 >> tet.p3 >> tet.p4))
        {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue; // Skip to the next line if parsing fails.
        }

        std::cout << tet.p1 << " " << tet.p2 << " " << tet.p3 << " " << tet.p4 << std::endl;

        // Store the extracted data in the map.
        tetMap.push_back(tet);
    }

    return tetMap;
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
          //  prefix sum (xfer + int - previous hArc)   1   4   7  13  4   13  15  1 (29-28=1)

          // intrinsic - weight of the superarc itself
          // dependent - weight of itself + dependent subtree (this is unambiguous because the superarc is oriented)
          // transfer weight - sum of the dependent weights of the subtrees connected to the supernode that have already been processed
          // once all of the dependent weights have been transfered to the superarc, the supernode can then be processed in another iteration
          // of the hypersweep, at which point the transfer weight is added to the intrinsic weight and becomes the dependent weight of the superarc.




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


    void static print2DarrayInt(std::vector<std::vector<int>> vxtc)
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




    vtkm::Vec3f_32 static fromTo(vtkm::Vec3f_32 a, vtkm::Vec3f_32 b)
    {
        vtkm::Vec3f_32 start = a;
        vtkm::Vec3f_32 result;

        result = start + (a - b);

        return result;
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

        // Files for basic 1D experiments of Contour Tree Branch node-count weight computations
        // const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-coordinates.txt";
        // const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-triangles.txt";

        // Files for 2D experiments of Contour Tree Branch length/area-based weight computations
        // (For debugging, I am currently keeping both 2D and 3D files, as to quickly flick between them to compare correctness)
        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-coordinates.txt";
        const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-triang.txt";

        // Files for 3D experiments of Contour Tree Branch volume-based weight computations
        // (For debugging, I am currently keeping both 2D and 3D files, as to quickly flick between them to compare correctness)
        const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Cube-8-coordinates.txt";
        const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Cube-8-tets.txt";


        // get the total number of values to be swept ...
        // ... this will be one more than the total number of datapoints ...
        // ... because we include the region beyond the last isovalue N, as a range [N, +inf)
        // (also, the number of values corresponds to the total number of different data values in a dataset ...
        //  ... because of the simulation of simplicity, which ensures ever data point has a unique value)
        int num_sweep_values = contourTree.Arcs.GetNumberOfValues() + 1;

        // Coordinate lists and connectivities of 2D triangle-based and 3D tetrahedral datasets.
        // The coordinates are needed for computing the areas/volumes.
        // Note that the coordinate information is not previously needed for the base contour tree computation ...
        // ... and is only introduced as a requirement at this stage.
        // (For debugging, I am currently keeping both 2D and 3D files, as to quickly flick between them to compare correctness)
        std::vector<Coordinates> coordlist    = ReadCoordinatesFromFile(filename1);
        // triangle list has pairs of 3 vertices that make up the triangle
        std::vector<Triangle> trianglelist    = ReadTrianglesFromFile(filename2);

        std::vector<Coordinates> coordlist3D  = ReadCoordinatesFromFile(filename3D1);
        // tetrahedron list has pairs of 4 vertices that make up the tetrahedron
        std::vector<Tetrahedron> tetlist      = ReadTetsFromFile(filename3D2);

        // Weight list is used for the basic implementation that was used for my learning.
        // It is only used in the naive area weight implementation, ...
        // ... where the weight of each vertex of the triangle is area/3 (area div by 3)
        std::vector<double> weightList;

        std::cout << "PRINT THE ARRAYS OF COORDINATES: \n";

        // Print the coordinates data to check if it was read correctly.
        for (int i = 0; i < coordlist.size(); i++)
        {
            std::cout << i << ": " << coordlist[i].x << ", " << coordlist[i].y << ", " << coordlist[i].z << std::endl;
            // initialise the weight list array while at it
            weightList.push_back(0.0);
        }

        std::cout << "PRINT THE ARRAYS OF TRIANGLES: \n";
        std::cout << "num. of triangles: " << trianglelist.size() << std::endl;
        for (int i = 0; i < trianglelist.size(); i++)
        {
            std::cout << i << ": " << trianglelist[i].p1 << ", " << trianglelist[i].p2 << ", " << trianglelist[i].p3;
            double area = ComputeTriangleArea(
                            coordlist[trianglelist[i].p1].x, coordlist[trianglelist[i].p1].y, coordlist[trianglelist[i].p1].z,
                            coordlist[trianglelist[i].p2].x, coordlist[trianglelist[i].p2].y, coordlist[trianglelist[i].p2].z,
                            coordlist[trianglelist[i].p3].x, coordlist[trianglelist[i].p3].y, coordlist[trianglelist[i].p3].z);
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

        // Now move on to set up 3D tetrahedral variables

        std::cout << "PRINT THE ARRAYS OF TETS: \n";
        std::cout << "num. of tets: " << tetlist.size() << std::endl;

        // keep track of all tets in a sorted list
        // (the 'sort' refers to sorting by the vertex ID.
        // The vertex ID corresponds to a data value, ...
        //  ... for example a tet X, Y, Z, W might have vertices with IDs X=6, Y=5, Z=3, W=7)
        //  ... we sort them as 3, 5, 6, 7 ...
        //  ... and refer to vertices as A, B, C, D
        //  ... where we then assign A = 3, B = 5, C = 6, D = 7)
        //  (we do not keep track of the original ordering X, Y, Z, W) ...
        std::vector<std::vector<int>> tetlistSorted(tetlist.size(),
                                                    std::vector<int> (4, 0));

        // Vertices that define the Tetrahedron ABCD (Entire tetrahedron) ...
        // ... with their corresponding isovalues
        std::vector<vtkm::Vec3f_32> verticesA; // vertex A contains the lowest isovalue h1
        std::vector<int> teth1s;
        std::vector<vtkm::Vec3f_32> verticesB; // vertex B - h2
        std::vector<int> teth2s;
        std::vector<vtkm::Vec3f_32> verticesC; // vertex C - h3
        std::vector<int> teth3s;
        std::vector<vtkm::Vec3f_32> verticesD; // vertex D - h4
        std::vector<int> teth4s;

        // Deriving middle slab triangle vertices E, F, G, H
        // Plane Points at isovalue h=h2 (4) (for interval h1->h2)              - FIRST TET
        std::vector<vtkm::Vec3f_32> verticesE;
        std::vector<vtkm::Vec3f_32> verticesF;

        // Plane Points at isovalue h=h3 (5) (for interval h4->h3)              - LAST TET
        std::vector<vtkm::Vec3f_32> verticesG;
        std::vector<vtkm::Vec3f_32> verticesH;

        // Plane Points between isovalues h2 and h3  [Vertices P, Q, R, S]      - MIDDLE QUAD SLAB
        std::vector<vtkm::Vec3f_32> verticesP;
        std::vector<vtkm::Vec3f_32> verticesQ;
        std::vector<vtkm::Vec3f_32> verticesR;
        std::vector<vtkm::Vec3f_32> verticesS;

//        // Initializing the 2-D vector
//        std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
//                                              std::vector<double> (trianglelist.size(), 0.0));


        for (int i = 0; i < tetlist.size(); i++)
        {
            std::cout << i << ": " << tetlist[i].p1 << ", " << tetlist[i].p2 << ", " << tetlist[i].p3 << ", " << tetlist[i].p4; //<< std::endl;
            double volume = 0.1666666666667; // HARDCODED TODO CHANGE
            std::cout << " = " << volume << std::endl;

            tetlistSorted[i][0] = tetlist[i].p1;
            tetlistSorted[i][1] = tetlist[i].p2;
            tetlistSorted[i][2] = tetlist[i].p3;
            tetlistSorted[i][3] = tetlist[i].p4;
        }

//        std::sort()
        std::cout << "PRINT TETS UNSORTED:" << std::endl;
        print2DarrayInt(tetlistSorted);

        for (int i = 0; i < tetlist.size(); i++)
        {
            std::sort(tetlistSorted[i].begin(), tetlistSorted[i].end());
        }

        std::cout << "PRINT TETS SORTED:" << std::endl;
        print2DarrayInt(tetlistSorted);

        std::cout << "============================================================================================" << std::endl;

        //        std::sort(tetlist.begin(), tetlist.end(), std::greater<int>());

        bool dim1 = false;
        bool dim2 = false;
        bool dim3 = true;

        // ----------------------------------- PRE-PROCESS ----------------------------------- //
        // Here we basically populate the table as such:
        // vertexID / triangleID
        //          triangle1 triangle2   row sum(del. h1/h2)   col (prefix) sum del. h1/h2
        // deltas:  h1 h2     h1 h2
        // v1       x  y      z  w
        std::vector<double> vx_delta_h1;
        std::vector<double> vx_delta_h1_sum;
        std::vector<double> vx_down_delta_h1_sum;

        std::vector<double> vx_delta_h2;     // 1D coefficient deltas
        std::vector<double> vx_delta_h2_sum; // 1D coefficients
        std::vector<double> vx_down_delta_h2_sum; // 1D coefficients

        std::vector<double> vx_delta_h3;     // 2D coefficient deltas
        std::vector<double> vx_delta_h3_sum; // 2D coefficients
        std::vector<double> vx_down_delta_h3_sum; // 2D coefficients

        std::vector<double> vx_delta_h4;     // 3D coefficient deltas
        std::vector<double> vx_delta_h4_sum; // 3D coefficients
        std::vector<double> vx_down_delta_h4_sum; // 3D coefficients



        std::vector<double> delta_h1_pfixsum;
        std::vector<double> delta_h2_pfixsum; // 1D coefficients
        std::vector<double> delta_h3_pfixsum; // 2D coefficients
        std::vector<double> delta_h4_pfixsum; // 3D coefficients

        if(dim1)
        {
            // (below is actually the length of the hypotenuse projected onto the base from the 90-degree angle)
            // TODO: replace this with actual area of each triangle!
            double fake_area = sqrt(2.0) / 2.0;


            std::cout << "PREPROCESSING STEP:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));

            std::cout << "\nInitialised vxtch1:" << std::endl;
            print2Darray(vxtch1);

            std::vector<std::vector<double>> vxtch2(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));
            std::cout << "\nInitialised vxtch2:" << std::endl;
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
            std::cout << "vxtch1" << std::endl;
            print2Darray(vxtch1);
            std::cout << "vxtch2" << std::endl;
            print2Darray(vxtch2);


            // now with deltas computed, also compute the prefix sums of the different components:

            delta_h1_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
            delta_h2_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());

            std::cout << "deltas:" << std::endl;

            for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            {
                vx_delta_h1_sum.push_back(0.0);
                vx_delta_h2_sum.push_back(0.0);

                for (int j = 0; j < trianglelist.size(); j++)
                {
                    vx_delta_h1_sum[i] += vxtch1[i][j];
                    vx_delta_h2_sum[i] += vxtch2[i][j];
                }

                std::cout << i << ") del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " ~ "; //<< std::endl;

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

                std::cout << i << ") pfix(h1)= (" << i << ")" << delta_h1_pfixsum[i] << " pfix(h2)=" << delta_h2_pfixsum[i] << std::endl;
            }
        }
        else if(dim2)
        {
            std::cout << "Computing 2D coefficients:" << std::endl;

            // (below is the area of a full 1x1 90-degree triangle)
            // TODO: replace this with actual area of each triangle!
            double fake_area = 0.5; //sqrt(2.0) / 2.0;

            double a_mid[] = {0.2148571, 0.0625, 0.16667, 0.25, 0.4375, 0.14285, 0.16667, 0.375};

            std::cout << "Printing a_mid areas:" << std::endl;
            for(int i = 0; i < 8; i++)
            {
                std::cout << a_mid[i] << " ";
            }
            std::cout << std::endl;

            std::cout << "PREPROCESSING STEP:" << std::endl;

            std::cout << "// ----------------------------------- PRE-PROCESS ----------------------------------- //" << std::endl;


            // Initializing the 2-D vector
            std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));

            std::cout << "\nInitialised vxtch1:" << std::endl;
            print2Darray(vxtch1);

            std::vector<std::vector<double>> vxtch2(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));
            std::cout << "\nInitialised vxtch2:" << std::endl;
            print2Darray(vxtch2);

            std::vector<std::vector<double>> vxtch3(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));
            print2Darray(vxtch3);
            std::cout << "\nInitialised vxtch3:" << std::endl;

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

                double Dn = Mv - Lv;
                double Dm = Hv - Mv;

                double Ba1 = 1.0/Dn;
                double Ba2 = Lv/Dn;

                double Bb1 = Hv / Dm;
                double Bb2 = 1.0/ Dm;


                // deltas:
                // Low-vertex deltas (same as original since adding from 0)
                double H12D = Ba1*Ba1 * a_mid[i];
                double H22D = 2 * Ba1 * Ba2 * a_mid[i];
                double H32D = Ba2*Ba2 * a_mid[i];

//                vxtch1[trianglelist[i].p3][i] = H12D;
//                vxtch2[trianglelist[i].p3][i] = H22D;
//                vxtch3[trianglelist[i].p3][i] = H32D;

                double H1b2D = -1.0 * (Bb2*Bb2) * (fake_area - a_mid[i]);
                double H2b2D = -2.0 * Bb2 * Bb1 * (fake_area - a_mid[i]);
                double H3b2D = fake_area - (Bb1 * Bb1) * (fake_area - a_mid[i]);

//                // Middle-vertex deltas:
//                vxtch1[trianglelist[i].p2][i] = H1b2D - H12D;
//                vxtch2[trianglelist[i].p2][i] = H2b2D - H22D;
//                vxtch3[trianglelist[i].p2][i] = H3b2D - H32D;

//                // High-vertex deltas:
//                vxtch1[trianglelist[i].p1][i] = -H1b2D;
//                vxtch2[trianglelist[i].p1][i] = -H2b2D;
//                vxtch3[trianglelist[i].p1][i] = -H3b2D;

                // invert operation:
//                // low:
//                vxtch1[trianglelist[i].p3][i] = -((H1b2D - H12D) -H1b2D); // = H12D
//                vxtch2[trianglelist[i].p3][i] = -((H2b2D - H22D) -H2b2D); // = H22D
//                vxtch3[trianglelist[i].p3][i] = -((H3b2D - H32D) -H3b2D); // = H32D

                // middle:
                vxtch1[trianglelist[i].p2][i] = -(-H1b2D + H12D);
                vxtch2[trianglelist[i].p2][i] = -(-H2b2D + H22D);
                vxtch3[trianglelist[i].p2][i] = -(-H3b2D + H32D);

                // high:
                vxtch1[trianglelist[i].p1][i] = -(H12D + (H1b2D - H12D));
                vxtch2[trianglelist[i].p1][i] = -(H22D + (H2b2D - H22D));
                vxtch3[trianglelist[i].p1][i] = fake_area-(H32D + (H3b2D - H32D));

                // low:
                vxtch1[trianglelist[i].p3][i] = -(vxtch1[trianglelist[i].p2][i] + vxtch1[trianglelist[i].p1][i]);  // -((H1b2D - H12D) -H1b2D); // = H12D
                vxtch2[trianglelist[i].p3][i] = -(vxtch2[trianglelist[i].p2][i] + vxtch2[trianglelist[i].p1][i]);  // -((H2b2D - H22D) -H2b2D); // = H22D
                vxtch3[trianglelist[i].p3][i] = fake_area-(vxtch3[trianglelist[i].p2][i] + vxtch3[trianglelist[i].p1][i]);  // fake_area-((H3b2D - H32D) -H3b2D); // = H32D

    //            std::cout << " " << i << ") " << delta_h1 << ", " << delta_h2 << std::endl;
            }


//            else if (dim4)
//            {

//            }

            std::cout << "vxtch1" << std::endl;
            print2Darray(vxtch1);
            std::cout << "vxtch2" << std::endl;
            print2Darray(vxtch2);
            std::cout << "vxtch3" << std::endl;
            print2Darray(vxtch3);


            // now with deltas computed, also compute the prefix sums of the different components:

//            std::vector<double> vx_delta_h1_sum;
//            std::vector<double> vx_delta_h2_sum;

//            std::vector<double> delta_h1_pfixsum; // = 0.0;
//            std::vector<double> delta_h2_pfixsum; // = 0.0;

            delta_h1_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
            delta_h2_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
            delta_h3_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());

            std::cout << "2D deltas:" << std::endl;

            for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            {
                vx_delta_h1_sum.push_back(0.0);
                vx_delta_h2_sum.push_back(0.0);
                vx_delta_h3_sum.push_back(0.0);

                for (int j = 0; j < trianglelist.size(); j++)
                {
                    vx_delta_h1_sum[i] += vxtch1[i][j];
                    vx_delta_h2_sum[i] += vxtch2[i][j];
                    vx_delta_h3_sum[i] += vxtch3[i][j];
                }

                std::cout << i << ") del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " del(h3)=" << vx_delta_h3_sum[i] << " ~ "; //<< std::endl;

                if (i == 0)
                {
                    delta_h1_pfixsum[0] = vx_delta_h1_sum[0];
    //                delta_h1_pfixsum.push_back(vx_delta_h1_sum[0]);
                    delta_h2_pfixsum[0] = vx_delta_h2_sum[0];
    //                delta_h2_pfixsum.push_back(vx_delta_h2_sum[0]);
                    delta_h3_pfixsum[0] = vx_delta_h3_sum[0];
                }
                else
                {
                    delta_h1_pfixsum[i] += delta_h1_pfixsum[i-1] + vx_delta_h1_sum[i];
                    delta_h2_pfixsum[i] += delta_h2_pfixsum[i-1] + vx_delta_h2_sum[i];
                    delta_h3_pfixsum[i] += delta_h3_pfixsum[i-1] + vx_delta_h3_sum[i];
                }

                std::cout << i << ") pfix(h1)= (" << i << ")" << delta_h1_pfixsum[i] << " pfix(h2)=" << delta_h2_pfixsum[i] << " pfix(h3)=" << delta_h3_pfixsum[i] << std::endl;
            }
        }
        else if(dim3)
        {
            std::cout << "Computing 3D Coefficients ..." << std::endl;

            // process one tet at a time ...
            // ... starting with they local coordinates and individual volumes
            for(int i = 0; i < tetlist.size(); i++)
            {
                std::cout << tetlist[i].p1 << " " << tetlist[i].p2 << " " << tetlist[i].p3 << " " << tetlist[i].p4 << std::endl;
                std::cout << "(" << tetlistSorted[i][0] << " " << tetlistSorted[i][1] << " " << tetlistSorted[i][2] << " " << tetlistSorted[i][3] << ")" << std::endl;

                std::cout << "\t " << tetlistSorted[i][0] << " = <" << coordlist3D[tetlistSorted[i][0]].x << " " << coordlist3D[tetlistSorted[i][0]].y << " " << coordlist3D[tetlistSorted[i][0]].z << ">" << std::endl;
                std::cout << "\t " << tetlistSorted[i][1] << " = <" << coordlist3D[tetlistSorted[i][1]].x << " " << coordlist3D[tetlistSorted[i][1]].y << " " << coordlist3D[tetlistSorted[i][1]].z << ">" << std::endl;
                std::cout << "\t " << tetlistSorted[i][2] << " = <" << coordlist3D[tetlistSorted[i][2]].x << " " << coordlist3D[tetlistSorted[i][2]].y << " " << coordlist3D[tetlistSorted[i][2]].z << ">" << std::endl;
                std::cout << "\t " << tetlistSorted[i][3] << " = <" << coordlist3D[tetlistSorted[i][3]].x << " " << coordlist3D[tetlistSorted[i][3]].y << " " << coordlist3D[tetlistSorted[i][3]].z << ">" << std::endl;
            }


            // ==================== \/ Step 1: Name the vertices that define each tetrahedron A, B, C, D \/ ==================== //
            // ... together with h1, h2, h3, h4
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                verticesA.emplace_back( coordlist3D[tetlistSorted[i][0]].x, coordlist3D[tetlistSorted[i][0]].y, coordlist3D[tetlistSorted[i][0]].z );
                verticesB.emplace_back( coordlist3D[tetlistSorted[i][1]].x, coordlist3D[tetlistSorted[i][1]].y, coordlist3D[tetlistSorted[i][1]].z );
                verticesC.emplace_back( coordlist3D[tetlistSorted[i][2]].x, coordlist3D[tetlistSorted[i][2]].y, coordlist3D[tetlistSorted[i][2]].z );
                verticesD.emplace_back( coordlist3D[tetlistSorted[i][3]].x, coordlist3D[tetlistSorted[i][3]].y, coordlist3D[tetlistSorted[i][3]].z );

                // keep track of tetrahedron boundary values h1, h2, h3, and h4
                teth1s.emplace_back(tetlistSorted[i][0]);
                teth2s.emplace_back(tetlistSorted[i][1]);
                teth3s.emplace_back(tetlistSorted[i][2]);
                teth4s.emplace_back(tetlistSorted[i][3]);
            }

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << "===" << std::endl;
                std::cout << "\t" << teth1s[i] << " = " << verticesA[i] << std::endl;
                std::cout << "\t" << teth2s[i] << " = " << verticesB[i] << std::endl;
                std::cout << "\t" << teth3s[i] << " = " << verticesC[i] << std::endl;
                std::cout << "\t" << teth4s[i] << " = " << verticesD[i] << std::endl;
                std::cout << "===" << std::endl;
            }

            // =================== \/ Step 1.5: from points A, B, C, D define vectors between them \/ ================== //

//            std::vector<vtkm::Vec3f_32> vectorsAB;
//            std::vector<vtkm::Vec3f_32> vectorsAC;

            std::vector<PositionVector> vectorsAB;
            std::vector<PositionVector> vectorsAC;
            std::vector<PositionVector> vectorsAD;
            std::vector<PositionVector> vectorsBD;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsAB.emplace_back(verticesA[i], verticesB[i]);
                vectorsAC.emplace_back(verticesA[i], verticesC[i]);
                vectorsAD.emplace_back(verticesA[i], verticesD[i]);
                vectorsBD.emplace_back(verticesB[i], verticesD[i]);

            }

            std::cout << "Vectors AB for each tet:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << i << " " << vectorsAD[i].start << " " << vectorsAD[i].end << " - " << vectorsAD[i].difference << std::endl;
                std::cout << i << " " << vectorsAC[i].start << " " << vectorsAC[i].end << " + " << vectorsAC[i].difference << std::endl;
//                vectorsAB[i].lerp(0.5);
//                std::cout << i << " " << vectorsAB[i].start << " " << vectorsAB[i].end << " + " << vectorsAB[i].difference << std::endl;
            }


            // =================== /\ Step 1.5: from points A, B, C, D define vectors between them /\ ================== //

            // ==================== \/ Step 2: Deriving middle slab triangle vertices E, F, G, H \/ ==================== //

            // Plane Points at isovalue h=h2 (4) (for interval h1->h2)              - FIRST TET

            // We will now be computing:
            //            verticesE;
            //            verticesF;
            // ... from vectors AD and AC respectively.


            for(int i = 0; i < tetlistSorted.size(); i++)
            {
                // lerp at h=2 on AD between 1 and 4 (a=1, d=4)
                double lerpADh2_h1_h4 = double(teth2s[i] - teth1s[i]) / double(teth4s[i] - teth1s[i]);
                verticesE.push_back(vectorsAD[i].lerp2point(lerpADh2_h1_h4));

                double lerpACh2_h1_h3 = double(teth2s[i] - teth1s[i]) / double(teth3s[i] - teth1s[i]);
                verticesF.push_back(vectorsAC[i].lerp2point(lerpACh2_h1_h3));

                std::cout << "lerpADh2_h1_h4: " << lerpADh2_h1_h4 << " E = " << verticesE[i] << std::endl;
                std::cout << "lerpACh2_h1_h3: " << lerpACh2_h1_h3 << " F = " << verticesF[i] << std::endl;
            }


            // Plane Points at isovalue h=h3 (5) (for interval h4->h3)              - LAST TET SLAB
            // We will now be computing:
            //            verticesG;
            //            verticesH;
            // ... from vectors AD and BD respectively.

            for(int i = 0; i < tetlistSorted.size(); i++)
            {
                // lerp at h=3 on AD between h2 and h3 (a=h1, d=h4)
                double lerpADh3_h1_h4 = double(teth3s[i] - teth1s[i]) / double(teth4s[i] - teth1s[i]);
                verticesG.push_back(vectorsAD[i].lerp2point(lerpADh3_h1_h4));

                // lerp at h3 on BD between h2 and h4 (b=h2, d=h4)

                double lerpBDh2_h2_h4 = double(teth3s[i] - teth2s[i]) / double(teth4s[i] - teth2s[i]);
                verticesH.push_back(vectorsBD[i].lerp2point(lerpBDh2_h2_h4));

                std::cout << "lerpADh3_h1_h4: " << lerpADh3_h1_h4 << " G = " << verticesG[i] << std::endl;
                std::cout << "lerpBDh2_h2_h4: " << lerpBDh2_h2_h4 << " H = " << verticesH[i] << std::endl;
            }


            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << "===" << std::endl;
                std::cout << "\t" << teth1s[i] << " = " << verticesA[i] << std::endl;
                std::cout << "\t" << teth2s[i] << " = " << verticesB[i] << std::endl;
                std::cout << "\t" << teth3s[i] << " = " << verticesC[i] << std::endl;
                std::cout << "\t" << teth4s[i] << " = " << verticesD[i] << std::endl;
                std::cout << "===" << std::endl;
            }
            // ==================== /\ Step 2: Deriving middle slab triangle vertices E, F, G, H /\ ==================== //




            // =========================== \/ Step 3: Compute Entire (full) Tet Volumes \/ ============================ //

            std::vector<double> full_tet_volumes;

            std::cout << "TET VOLUMES:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
//                vtkm::Vec3f_32 a_vol = verticesA[i];
//                vtkm::Vec3f_32 b_vol = verticesB[i];
//                vtkm::Vec3f_32 c_vol = verticesC[i];

                PositionVector a_vol(verticesA[i], verticesC[i]);
                PositionVector b_vol(verticesA[i], verticesD[i]);
                PositionVector c_vol(verticesA[i], verticesB[i]);

                full_tet_volumes.push_back((1.0/6.0) * abs(vtkm::Dot(vtkm::Cross(a_vol.difference, b_vol.difference), c_vol.difference) ) );
                std::cout << i << " = " << full_tet_volumes[i] << std::endl;
            }



            // =========================== /\ Step 3: Compute Entire (full) Tet Volumes  /\ ============================ //



            // ---------------------------------------------- FIRST SLAB ----------------------------------------------- //



            // ============== \/ Step 4: Compute the first slab volume (defined from isovalues h1-h2) \/ =============== //

            std::vector<double> slab1_h1h2_tet_volumes;

            std::cout << "SLAB1 VOLUMES:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
//                vtkm::Vec3f_32 a_vol = verticesA[i];
//                vtkm::Vec3f_32 b_vol = verticesB[i];
//                vtkm::Vec3f_32 c_vol = verticesC[i];
                double alpha_h2= std::max(0.0, std::min(1.0,
                                          double(teth2s[i]-teth1s[i])/double(teth3s[i]-teth1s[i])));

                double beta_h2 = std::max(0.0, std::min(1.0,
                                          double(teth2s[i]-teth1s[i])/double(teth4s[i]-teth1s[i])));

                // NOTE: gamma can be optimised away, since we reach point B at h2 (because B holds h2) ...
                // ... gamma will always be 1.0
                double gamma_h2= std::max(0.0, std::min(1.0,
                                          double(teth2s[i]-teth1s[i])/double(teth2s[i]-teth1s[i])));


                PositionVector a_h1h2_vol(verticesA[i], verticesC[i]);
                PositionVector b_h1h2_vol(verticesA[i], verticesD[i]);
                PositionVector c_h1h2_vol(verticesA[i], verticesB[i]);

                a_h1h2_vol.lerp(alpha_h2);
                b_h1h2_vol.lerp(beta_h2);
                c_h1h2_vol.lerp(gamma_h2);

                slab1_h1h2_tet_volumes.push_back((1.0/6.0) * abs(vtkm::Dot(vtkm::Cross(a_h1h2_vol.difference, b_h1h2_vol.difference),
                                                                           c_h1h2_vol.difference) ) );
                std::cout << i << " = " << slab1_h1h2_tet_volumes[i] << "(" << alpha_h2 << ", " << beta_h2 << ", " << gamma_h2 << ")" << std::endl;
            }



            // ==========  \/ Step 5: Compute the first slab coefficients (defined from isovalues h1-h2) \/ ============ //

            //
            std::vector<double> a_h1h2;
            std::vector<double> b_h1h2;
            std::vector<double> c_h1h2;
            std::vector<double> d_h1h2;


            std::cout << "SLAB1 h1h2 coefficients:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                a_h1h2.push_back(slab1_h1h2_tet_volumes[i]                                  / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );
                b_h1h2.push_back(-(3.0 * slab1_h1h2_tet_volumes[i] * teth1s[i])              / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );
                c_h1h2.push_back((3.0 * slab1_h1h2_tet_volumes[i] * std::pow(teth1s[i], 2)) / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );
                d_h1h2.push_back(-(slab1_h1h2_tet_volumes[i] * std::pow(teth1s[i], 3))      / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );


                std::cout << i << " = " << "a = " << a_h1h2[i]<< ", b = " << b_h1h2[i]<< ", c = " << c_h1h2[i] << ", d = " << d_h1h2[i] << std::endl;
            }

            // ========== /\ Step 5: Compute the first slab coefficients (defined from isovalues h1-h2) /\ ============ //



            // ---------------------------------------------- LAST SLAB ----------------------------------------------- //




            // ============== \/ Step 6: Compute the last slab volume (defined from isovalues h3-h4) \/ =============== //


            std::vector<double> slab3_h3h4_tet_volumes; // only the slab volume
            std::vector<double> slab3_h3h4_tet_volumes_sweeping_up; // (total volume) - (slab)

            std::cout << "SLAB3 VOLUMES:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
//                vtkm::Vec3f_32 a_vol = verticesA[i];
//                vtkm::Vec3f_32 b_vol = verticesB[i];
//                vtkm::Vec3f_32 c_vol = verticesC[i];
                double alpha_h3= std::max(0.0, std::min(1.0,
                                          double(teth3s[i]-teth4s[i])/double(teth2s[i]-teth4s[i])));

                double beta_h3 = std::max(0.0, std::min(1.0,
                                          double(teth3s[i]-teth4s[i])/double(teth1s[i]-teth4s[i])));

                double gamma_h3= std::max(0.0, std::min(1.0,
                                          double(teth3s[i]-teth4s[i])/double(teth3s[i]-teth4s[i])));


                PositionVector a_h3h4_vol(verticesD[i], verticesB[i]);
                PositionVector b_h3h4_vol(verticesD[i], verticesA[i]);
                PositionVector c_h3h4_vol(verticesD[i], verticesC[i]);

                a_h3h4_vol.lerp(alpha_h3);
                b_h3h4_vol.lerp(beta_h3);
                c_h3h4_vol.lerp(gamma_h3);

                slab3_h3h4_tet_volumes.push_back((1.0/6.0) * abs(vtkm::Dot(vtkm::Cross(a_h3h4_vol.difference, b_h3h4_vol.difference),
                                                                           c_h3h4_vol.difference) ) );

                slab3_h3h4_tet_volumes_sweeping_up.push_back(full_tet_volumes[i] - slab3_h3h4_tet_volumes[i]);

                std::cout << i << " = " << slab3_h3h4_tet_volumes[i] << "(" << alpha_h3 << ", " << beta_h3 << ", " << gamma_h3 << ")" << std::endl;
            }

            // ============== /\ Step 6: Compute the last slab volume (defined from isovalues h3-h4) /\ =============== //

            // ==========  \/ Step 7: Compute the last slab coefficients (defined from isovalues h3-h4) \/ ============ //

            //
            std::vector<double> a_h3h4;
            std::vector<double> b_h3h4;
            std::vector<double> c_h3h4;
            std::vector<double> d_h3h4;
            std::vector<double> d_h3h4_down;

            // compute the coefficients using 'full_tet_volumes' as v_volumeh4moveup ...
            // ... and 'slab3_h3h4_tet_volumes' as v_volumeh3moveup


            std::cout << "SLAB3 h3h4 coefficients:" << std::endl;
//            double d_h3h4; // dealing with the fourth coefficient separately, as we need to move the last segment up to start from slab volume at h3
            double d_h3h3_to0;              //       updated 'd' coefficient that moves the slab function to start at y=0
            std::vector<double> d_h3h4_up;  // final updated 'd' coefficient that moves the slab function to start at y=volume_h3
            double vol_h3 = 0.0;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << slab3_h3h4_tet_volumes_sweeping_up[i] << " " << full_tet_volumes[i] << std::endl;
                a_h3h4.push_back((slab3_h3h4_tet_volumes_sweeping_up[i] - full_tet_volumes[i])                                                              / double( std::pow((teth3s[i] - teth4s[i]), 3) ) );
                b_h3h4.push_back((-3.0*teth4s[i] * slab3_h3h4_tet_volumes_sweeping_up[i] + 3.0*teth4s[i] * full_tet_volumes[i])                             / double( std::pow((teth3s[i] - teth4s[i]), 3) ) );
                c_h3h4.push_back(( 3.0*std::pow(teth4s[i],2) * slab3_h3h4_tet_volumes_sweeping_up[i] - 3.0*std::pow(teth4s[i],2) * full_tet_volumes[i])     / double( std::pow((teth3s[i] - teth4s[i]), 3) ) );

                // compute the base d coefficient, update it later to lift the function up:
                d_h3h4.push_back((-std::pow(teth4s[i],3) * slab3_h3h4_tet_volumes_sweeping_up[i] + std::pow(teth4s[i],3) * full_tet_volumes[i])                     / double( std::pow((teth3s[i] - teth4s[i]), 3)));


                vol_h3 = a_h3h4[i]*std::pow(teth3s[i], 3) + b_h3h4[i]*std::pow(teth3s[i], 2) + c_h3h4[i]*teth3s[i] + d_h3h4[i];

                d_h3h3_to0 = d_h3h4[i] - vol_h3;

                // deal with the fourth coefficient (d - the constant) separately, as it helps to move the function up/down
                d_h3h4_up.push_back( d_h3h3_to0 + slab3_h3h4_tet_volumes_sweeping_up[i] );


                std::cout << i << " = " << "a = " << a_h3h4[i]<< ", b = " << b_h3h4[i]<< ", c = " << c_h3h4[i] << ", d = " << d_h3h4_up[i] << std::endl;
            }

            // ========== /\ Step 7: Compute the last slab coefficients (defined from isovalues h3-h4) /\ ============ //



            // ---------------------------------------------- MID SLAB ----------------------------------------------- //


            // =========================  \/ Step 8: Compute sin(theta1) and sin(theta2) \/ ========================== //
            // sin theta 1 computation:

            std::vector<double> areas_CGH;
            std::vector<double> sin_theta_1s;

            std::vector<PositionVector> vectorsGH;
            std::vector<PositionVector> vectorsCH;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsGH.emplace_back(verticesG[i], verticesH[i]);
                vectorsCH.emplace_back(verticesC[i], verticesH[i]);
            }

            std::cout << "Areas CGH:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                areas_CGH.push_back(1.0/2.0 * vtkm::Magnitude(vtkm::Cross( vectorsGH[i].difference, vectorsCH[i].difference )) );
                sin_theta_1s.push_back(2.0 * areas_CGH[i]/ (vectorsGH[i].mag() * vectorsCH[i].mag()) );

                std::cout << "Area of " << i << " = " << areas_CGH[i] << std::endl;
                std::cout << "sin(theta1) of " << i << " = " << sin_theta_1s[i] << std::endl;
            }



            // sin theta 2 computation:

            std::vector<double> areas_BEF;
            std::vector<double> sin_theta_2s;

            std::vector<PositionVector> vectorsFB;
            std::vector<PositionVector> vectorsFE;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsFB.emplace_back(verticesF[i], verticesB[i]);
                vectorsFE.emplace_back(verticesF[i], verticesE[i]);
            }

            std::cout << "Areas BEF:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                areas_BEF.push_back(1.0/2.0 * vtkm::Magnitude(vtkm::Cross( vectorsFB[i].difference, vectorsFE[i].difference )) );
                sin_theta_2s.push_back(2.0 * areas_BEF[i]/ (vectorsFB[i].mag() * vectorsFE[i].mag()) );

                std::cout << "Area of " << i << " = " << areas_BEF[i] << std::endl;
                std::cout << "sin(theta1) of " << i << " = " << sin_theta_2s[i] << std::endl;
            }


            // =========================  /\ Step 8: Compute sin(theta1) and sin(theta2) /\ ========================== //


            // ===========================  \/ Step 9: Compute Middle Slab Coefficients \/ =========================== //

            // ----------------------------------------------- BATCH 1 ----------------------------------------------- //
            // a_s1, b_s1, c_s1
            std::vector<PositionVector> vectorsHG;
            std::vector<PositionVector> vectorsBE;
            std::vector<PositionVector> vectorsHC;
            std::vector<PositionVector> vectorsCG;
            // FE - already defined
            // FB - already defined


            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsHG.emplace_back(verticesH[i], verticesG[i]);
                vectorsBE.emplace_back(verticesB[i], verticesE[i]);
                vectorsHC.emplace_back(verticesH[i], verticesC[i]);
                vectorsCG.emplace_back(verticesC[i], verticesG[i]);

            }


            std::vector<double> tetk2s;
            std::vector<double> a_s1;
            std::vector<double> b_s1;
            std::vector<double> c_s1;

            // local variables for simplifying notation:
            double n1;
            double n2;
            double n3;
            double n4;

            std::cout << "Mid slab pre-coefficients" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // noting down repeating terms as I am writing the code for the first time:
                //                (tetk2s[i] * (vectorsHG[i].mag() - vectorsBE[i].mag() )

                tetk2s.push_back( 1.0 / double(teth3s[i] - teth2s[i]) );

                n1 = (tetk2s[i] * (vectorsHG[i].mag() - vectorsBE[i].mag() ));
                n2 = (vectorsHC[i].mag() * tetk2s[i]);

                n3 = ( teth2s[i] * tetk2s[i] * vectorsHG[i].mag() * vectorsHC[i].mag() * tetk2s[i] );
                n4 = ( tetk2s[i] * teth3s[i] * vectorsBE[i].mag() * vectorsHC[i].mag() * tetk2s[i] );


                a_s1.push_back( n1 * n2 );
                b_s1.push_back( -( (n1 * n2 * teth2s[i]) + n3 - n4 ) );
                c_s1.push_back( n3 * teth2s[i] - n4 * teth2s[i] );

                std::cout << "[a_s1] " << i << " == " << a_s1[i] << std::endl;
                std::cout << "[b_s1] " << i << " == " << b_s1[i] << std::endl;
                std::cout << "[c_s1] " << i << " == " << c_s1[i] << std::endl;

            }

            // ----------------------------------------------- BATCH 2 ----------------------------------------------- //
            // a_s2, b_s2, c_s2

            // HG - already defined
            // BE - already defined
            // HC - already defined
            // CG - already defined
            // FE - already defined
            // FB - already defined

            std::vector<double> a_s2;
            std::vector<double> b_s2;
            std::vector<double> c_s2;

            // local variables for simplifying notation:
            double m1;
            double m2;
            double m3;
            double m4;

            std::cout << "Mid slab pre-coefficients" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                m1 = (tetk2s[i] * (vectorsCG[i].mag() - vectorsFE[i].mag() ));
                m2 = (-vectorsFB[i].mag() * tetk2s[i]);

                m3 = (-teth2s[i] * tetk2s[i] * vectorsCG[i].mag() * vectorsFB[i].mag() * tetk2s[i] );
                m4 = ( tetk2s[i] * teth3s[i] * vectorsFE[i].mag() * vectorsFB[i].mag() * tetk2s[i] );


                a_s2.push_back( m1 * m2 );
                b_s2.push_back( m1 * -m2 * teth3s[i] - m3 -m4);
                c_s2.push_back( m3 * teth3s[i] + m4 * teth3s[i] );

                std::cout << "[a_s2] " << i << " == " << a_s2[i] << std::endl;
                std::cout << "[b_s2] " << i << " == " << b_s2[i] << std::endl;
                std::cout << "[c_s2] " << i << " == " << c_s2[i] << std::endl;

            }

            // ------------------------------------------ COMBINE BATCH 1+2 ------------------------------------------- //
            std::vector<double> a_mid;
            std::vector<double> b_mid;
            std::vector<double> c_mid;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                a_mid.push_back(sin_theta_1s[i] / 2.0 * a_s1[i] + sin_theta_2s[i] / 2.0 * a_s2[i]);
                b_mid.push_back(sin_theta_1s[i] / 2.0 * b_s1[i] + sin_theta_2s[i] / 2.0 * b_s2[i]);
                c_mid.push_back(sin_theta_1s[i] / 2.0 * c_s1[i] + sin_theta_2s[i] / 2.0 * c_s2[i]);

                std::cout << "comb " << i << " " << a_mid[i] << " " << b_mid[i] << " " << c_mid[i] << std::endl;
            }


            // ---------------------------- Compute the Integration correction coefficient ---------------------------- //
            std::vector<vtkm::Vec3f_32> plane_normals;
            std::vector<double> plane_distances;
            vtkm::Vec3f_32 FExFB_cross_product;

            std::vector<double> correction_factor_nominators;


            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                FExFB_cross_product = vtkm::Cross(vectorsFE[i].difference, vectorsFB[i].difference);

                plane_normals.push_back( FExFB_cross_product / (vtkm::Magnitude(FExFB_cross_product)) );
                plane_distances.push_back( vtkm::Magnitude(vtkm::Dot(plane_normals[i], verticesB[i]) - vtkm::Dot(plane_normals[i], verticesH[i]) ) / (vtkm::Magnitude(plane_normals[i]) ) );
//                plane_distances.push_back( (vtkm::Dot(plane_normals[i], verticesB[i]) - vtkm::Dot(plane_normals[i], verticesH[i]) ) / (vtkm::Magnitude(plane_normals[i]) ) );


                correction_factor_nominators.push_back(plane_distances[i] * tetk2s[i]);

                std::cout << "corr: " << i << " " << correction_factor_nominators[i] << std::endl;

            }


            // ---------------------------- Compute the Integration correction coefficient ---------------------------- //

            double d_h2h3_to0;
//            double d_h2h3_up;
//            std::vector<double> d_mid;

            std::vector<double> a_h2h3;
            std::vector<double> b_h2h3;
            std::vector<double> c_h2h3;
            std::vector<double> d_h2h3;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                a_h2h3.push_back(correction_factor_nominators[i]/3.0 * a_mid[i]);
                b_h2h3.push_back(correction_factor_nominators[i]/2.0 * b_mid[i]);
                c_h2h3.push_back(correction_factor_nominators[i]     * c_mid[i]);

                d_h2h3_to0 = a_h2h3[i]  * std::pow(teth2s[i], 3) +\
                             b_h2h3[i]  * std::pow(teth2s[i], 2) +\
                             c_h2h3[i]  *          teth2s[i];

//                d_mid.push_back();


                d_h2h3.push_back(-d_h2h3_to0 + slab1_h1h2_tet_volumes[i]);


            }

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << "mid coeffs " << i << " : " << a_h2h3[i] << " " << b_h2h3[i] << " " << c_h2h3[i] << " " << d_h2h3[i] << std::endl;
            }


            // ===========================  /\ Step 9: Compute Middle Slab Coefficients /\ ===========================  //





            // ============================  \/ Step 10: Compute UP coefficients table  \/ ============================  //

            std::cout << "Initialised coefficient tables:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_coeffs_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_coeffs_h1:" << std::endl;
            print2Darray(tet_coeffs_h1);

            std::vector<std::vector<double>> tet_coeffs_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_coeffs_h2:" << std::endl;
            print2Darray(tet_coeffs_h2);

            std::vector<std::vector<double>> tet_coeffs_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_coeffs_h3:" << std::endl;
            print2Darray(tet_coeffs_h3);

            std::vector<std::vector<double>> tet_coeffs_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_coeffs_h4:" << std::endl;
            print2Darray(tet_coeffs_h4);



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in just the h1 coeffs:
                for(int h = teth1s[i]; h < teth2s[i]; h++)
                {
                    tet_coeffs_h1[h][i] = a_h1h2[i];
                    tet_coeffs_h2[h][i] = b_h1h2[i];
                    tet_coeffs_h3[h][i] = c_h1h2[i];
                    tet_coeffs_h4[h][i] = d_h1h2[i];
                }
                for(int h = teth2s[i]; h < teth3s[i]; h++)
                {
                    tet_coeffs_h1[h][i] = a_h2h3[i];
                    tet_coeffs_h2[h][i] = b_h2h3[i];
                    tet_coeffs_h3[h][i] = c_h2h3[i];
                    tet_coeffs_h4[h][i] = d_h2h3[i];
                }
                for(int h = teth3s[i]; h < teth4s[i]; h++)
                {
                    tet_coeffs_h1[h][i] = a_h3h4[i];
                    tet_coeffs_h2[h][i] = b_h3h4[i];
                    tet_coeffs_h3[h][i] = c_h3h4[i];
                    tet_coeffs_h4[h][i] = d_h3h4_up[i];
                }
                for(int h = teth4s[i]; h < 8; h++)
                {
                    tet_coeffs_h1[h][i] = 0.0;
                    tet_coeffs_h2[h][i] = 0.0;
                    tet_coeffs_h3[h][i] = 0.0;
                    tet_coeffs_h4[h][i] = full_tet_volumes[i];
                }
            }

            std::cout << "\nh1 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h1);
            std::cout << "\nh2 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h2);
            std::cout << "\nh3 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h3);
            std::cout << "\nh4 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h4);



            // ============================  /\ Step 10: Compute UP coefficients table  /\ ============================  //





            // =========================  \/ Step 11: Compute UP coefficient deltas table  \/ =========================  //

            std::cout << "Initialised coefficient deltas tables:" << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "  << std::endl;
            std::cout << "Contour Tree Number of Arcs: " << contourTree.Arcs.GetNumberOfValues() << std::endl;
            std::cout << "Total Sweep Values then +1 : " << num_sweep_values << std::endl; // 2024-11-18 updated sweep value
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "  << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_deltas_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_deltas_h1:" << std::endl;
            print2Darray(tet_deltas_h1);

            std::vector<std::vector<double>> tet_deltas_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_deltas_h2:" << std::endl;
            print2Darray(tet_deltas_h2);

            std::vector<std::vector<double>> tet_deltas_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_deltas_h3:" << std::endl;
            print2Darray(tet_deltas_h3);

            std::vector<std::vector<double>> tet_deltas_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_deltas_h4:" << std::endl;
            print2Darray(tet_deltas_h4 );



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in just the h1 coeffs:
//                for(int h = teth1s[i]; h < teth2s[i]; h++)
//                {
                    tet_deltas_h1[teth1s[i]+1][i] = a_h1h2[i] - 0.0;
                    tet_deltas_h2[teth1s[i]+1][i] = b_h1h2[i] - 0.0;
                    tet_deltas_h3[teth1s[i]+1][i] = c_h1h2[i] - 0.0;
                    tet_deltas_h4[teth1s[i]+1][i] = d_h1h2[i] - 0.0;
//                }
//                for(int h = teth2s[i]; h < teth3s[i]; h++)
//                {
                    tet_deltas_h1[teth2s[i]+1][i] = a_h2h3[i] - a_h1h2[i];
                    tet_deltas_h2[teth2s[i]+1][i] = b_h2h3[i] - b_h1h2[i];
                    tet_deltas_h3[teth2s[i]+1][i] = c_h2h3[i] - c_h1h2[i];
                    tet_deltas_h4[teth2s[i]+1][i] = d_h2h3[i] - d_h1h2[i];
//                }
//                for(int h = teth3s[i]; h < teth4s[i]; h++)
//                {
                    tet_deltas_h1[teth3s[i]+1][i] = a_h3h4[i] - a_h2h3[i];
                    tet_deltas_h2[teth3s[i]+1][i] = b_h3h4[i] - b_h2h3[i];
                    tet_deltas_h3[teth3s[i]+1][i] = c_h3h4[i] - c_h2h3[i];
                    tet_deltas_h4[teth3s[i]+1][i] = d_h3h4_up[i] - d_h2h3[i];
//                }
//                for(int h = teth4s[i]; h < 8; h++)
//                {
                    tet_deltas_h1[teth4s[i]+1][i] = - a_h3h4[i];
                    tet_deltas_h2[teth4s[i]+1][i] = - b_h3h4[i];
                    tet_deltas_h3[teth4s[i]+1][i] = - c_h3h4[i];
                    tet_deltas_h4[teth4s[i]+1][i] = full_tet_volumes[i] - d_h3h4_up[i];
//                }
            }

            std::cout << "\nh1 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h1);
            std::cout << "\nh2 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h2);
            std::cout << "\nh3 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h3);
            std::cout << "\nh4 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h4);



            // =========================  /\ Step 11: Compute UP coefficient deltas table  /\ =========================  //






            // ============================ \/ Step 12: Compute DOWN coefficient tables  \/ ==========================  //

            std::vector<double> d_h1h2_down;
            std::vector<double> d_h2h3_down;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                d_h2h3_down.push_back( (  -a_h2h3[i] * std::pow(teth3s[i], 3) - b_h2h3[i] * std::pow(teth3s[i], 2) - c_h2h3[i] * teth3s[i] ) \
                                       - (-a_h3h4[i] * std::pow(teth3s[i], 3) - b_h3h4[i] * std::pow(teth3s[i], 2) - c_h3h4[i] * teth3s[i] - d_h3h4[i] ) );

                d_h1h2_down.push_back( (  -a_h1h2[i] * std::pow(teth2s[i], 3) - b_h1h2[i] * std::pow(teth2s[i], 2) - c_h1h2[i] * teth2s[i] ) \
                                       - (-a_h2h3[i] * std::pow(teth2s[i], 3) - b_h2h3[i] * std::pow(teth2s[i], 2) - c_h2h3[i] * teth2s[i] - d_h2h3_down[i] ) );

            }







            std::cout << "Initialised DOWN coefficient tables:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_down_coeffs_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_down_coeffs_h1:" << std::endl;
            print2Darray(tet_down_coeffs_h1);

            std::vector<std::vector<double>> tet_down_coeffs_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_coeffs_h2:" << std::endl;
            print2Darray(tet_down_coeffs_h2);

            std::vector<std::vector<double>> tet_down_coeffs_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_coeffs_h3:" << std::endl;
            print2Darray(tet_down_coeffs_h3);

            std::vector<std::vector<double>> tet_down_coeffs_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_coeffs_h4:" << std::endl;
            print2Darray(tet_down_coeffs_h4);



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in before the h1 coeffs:
                for(int h = 0; h < teth1s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = 0.0;
                    tet_down_coeffs_h2[h][i] = 0.0;
                    tet_down_coeffs_h3[h][i] = 0.0;
                    tet_down_coeffs_h4[h][i] = full_tet_volumes[i];
                }

                // fill in just the h1 coeffs:
                for(int h = teth1s[i]; h < teth2s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = -a_h1h2[i];
                    tet_down_coeffs_h2[h][i] = -b_h1h2[i];
                    tet_down_coeffs_h3[h][i] = -c_h1h2[i];
                    tet_down_coeffs_h4[h][i] = -d_h1h2_down[i];
                }
                for(int h = teth2s[i]; h < teth3s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = -a_h2h3[i];
                    tet_down_coeffs_h2[h][i] = -b_h2h3[i];
                    tet_down_coeffs_h3[h][i] = -c_h2h3[i];
                    tet_down_coeffs_h4[h][i] = -d_h2h3_down[i];
                }
                for(int h = teth3s[i]; h < teth4s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = -a_h3h4[i];
                    tet_down_coeffs_h2[h][i] = -b_h3h4[i];
                    tet_down_coeffs_h3[h][i] = -c_h3h4[i];
                    tet_down_coeffs_h4[h][i] = -d_h3h4[i];
                }
                for(int h = teth4s[i]; h < 8; h++)
                {
                    tet_down_coeffs_h1[h][i] = 0.0;
                    tet_down_coeffs_h2[h][i] = 0.0;
                    tet_down_coeffs_h3[h][i] = 0.0;
                    tet_down_coeffs_h4[h][i] = 0.0;
                }
            }

            std::cout << "\nh1 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h1);
            std::cout << "\nh2 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h2);
            std::cout << "\nh3 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h3);
            std::cout << "\nh4 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h4);



            // ============================ /\ Step 12: Compute DOWN coefficient tables  /\ ==========================  //








            // ======================== \/ Step 13: Compute DOWN coefficient delta tables  \/ =======================  //

            std::cout << "Initialised DOWN coefficient deltas tables:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_down_deltas_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_down_deltas_h1:" << std::endl;
            print2Darray(tet_down_deltas_h1);

            std::vector<std::vector<double>> tet_down_deltas_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_deltas_h2:" << std::endl;
            print2Darray(tet_down_deltas_h2);

            std::vector<std::vector<double>> tet_down_deltas_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_deltas_h3:" << std::endl;
            print2Darray(tet_down_deltas_h3);

            std::vector<std::vector<double>> tet_down_deltas_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_deltas_h4:" << std::endl;
            print2Darray(tet_down_deltas_h4 );



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in just the h1 coeffs:
//                for(int h = teth1s[i]; h < teth2s[i]; h++)
//                {
                    tet_down_deltas_h1[teth1s[i]][i] = 0.0 + a_h1h2[i];
                    tet_down_deltas_h2[teth1s[i]][i] = 0.0 + b_h1h2[i];
                    tet_down_deltas_h3[teth1s[i]][i] = 0.0 + c_h1h2[i];
                    tet_down_deltas_h4[teth1s[i]][i] = full_tet_volumes[i] + d_h1h2_down[i];
//                }
//                for(int h = teth2s[i]; h < teth3s[i]; h++)
//                {
                    tet_down_deltas_h1[teth2s[i]][i] = -a_h1h2[i] + a_h2h3[i];
                    tet_down_deltas_h2[teth2s[i]][i] = -b_h1h2[i] + b_h2h3[i];
                    tet_down_deltas_h3[teth2s[i]][i] = -c_h1h2[i] + c_h2h3[i];
                    tet_down_deltas_h4[teth2s[i]][i] = -d_h1h2_down[i] + d_h2h3_down[i];
//                }
//                for(int h = teth3s[i]; h < teth4s[i]; h++)
//                {
                    tet_down_deltas_h1[teth3s[i]][i] = -a_h2h3[i] + a_h3h4[i];
                    tet_down_deltas_h2[teth3s[i]][i] = -b_h2h3[i] + b_h3h4[i];
                    tet_down_deltas_h3[teth3s[i]][i] = -c_h2h3[i] + c_h3h4[i];
                    tet_down_deltas_h4[teth3s[i]][i] = -d_h2h3_down[i] + d_h3h4[i];
//                }
//                for(int h = teth4s[i]; h < 8; h++)
//                {
                    tet_down_deltas_h1[teth4s[i]][i] = -a_h3h4[i];
                    tet_down_deltas_h2[teth4s[i]][i] = -b_h3h4[i];
                    tet_down_deltas_h3[teth4s[i]][i] = -c_h3h4[i];
                    tet_down_deltas_h4[teth4s[i]][i] = -d_h3h4[i];
//                }
            }

            std::cout << "\nh1 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h1);
            std::cout << "\nh2 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h2);
            std::cout << "\nh3 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h3);
            std::cout << "\nh4 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h4);


            // ======================== /\ Step 13: Compute DOWN coefficient delta tables  /\ =======================  //



            // ======================= \/ Step 14: Compute delta table up/down prefix sums \/ =======================  //

            std::vector<std::vector<double>> tet_up_deltas_pfix(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (4, 0.0));


            std::vector<std::vector<double>> tet_down_deltas_pfix(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (4, 0.0));

            for (int i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            // num_sweep_values, // 2024-11-18 updated sweep value
            for (int i = 0; i < num_sweep_values; i++)
            {
                double a_coeff_sum = 0.0;
                double b_coeff_sum = 0.0;
                double c_coeff_sum = 0.0;
                double d_coeff_sum = 0.0;

                // up deltas
                for(int k = 0; k < tetlistSorted.size(); k++)
                {
                    a_coeff_sum += tet_deltas_h1[i][k];
                    b_coeff_sum += tet_deltas_h2[i][k];
                    c_coeff_sum += tet_deltas_h3[i][k];
                    d_coeff_sum += tet_deltas_h4[i][k];
                }
                tet_up_deltas_pfix[i][0] =  a_coeff_sum;
                tet_up_deltas_pfix[i][1] =  b_coeff_sum;
                tet_up_deltas_pfix[i][2] =  c_coeff_sum;
                tet_up_deltas_pfix[i][3] =  d_coeff_sum;


                a_coeff_sum = 0.0;
                b_coeff_sum = 0.0;
                c_coeff_sum = 0.0;
                d_coeff_sum = 0.0;

                // down deltas
                for(int k = 0; k < tetlistSorted.size(); k++)
                {
                    a_coeff_sum += tet_down_deltas_h1[i][k];
                    b_coeff_sum += tet_down_deltas_h2[i][k];
                    c_coeff_sum += tet_down_deltas_h3[i][k];
                    d_coeff_sum += tet_down_deltas_h4[i][k];
                }
                tet_down_deltas_pfix[i][0] =  a_coeff_sum;
                tet_down_deltas_pfix[i][1] =  b_coeff_sum;
                tet_down_deltas_pfix[i][2] =  c_coeff_sum;
                tet_down_deltas_pfix[i][3] =  d_coeff_sum;

            }

            std::cout << "tet_up_deltas_pfix:" << std::endl;
            print2Darray(tet_up_deltas_pfix);

            std::cout << "tet_down_deltas_pfix:" << std::endl;
            print2Darray(tet_down_deltas_pfix);

            // ======================= /\ Step 14: Compute delta table up/down prefix sums /\ =======================  //

            std::cout << std::endl;

//            for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            for (vtkm::Id i = 0; i < num_sweep_values; i++)
            {
//                vx_delta_h1_sum.push_back(tet_up_deltas_pfix[i][0]);
//                vx_delta_h2_sum.push_back(tet_up_deltas_pfix[i][1]);
//                vx_delta_h3_sum.push_back(tet_up_deltas_pfix[i][2]);
//                vx_delta_h4_sum.push_back(tet_up_deltas_pfix[i][3]);

                vx_delta_h1_sum.push_back(tet_down_deltas_pfix[i][0]);
                vx_delta_h2_sum.push_back(tet_down_deltas_pfix[i][1]);
                vx_delta_h3_sum.push_back(tet_down_deltas_pfix[i][2]);
                vx_delta_h4_sum.push_back(tet_down_deltas_pfix[i][3]);

                vx_down_delta_h1_sum.push_back(tet_down_deltas_pfix[i][0]);
                vx_down_delta_h2_sum.push_back(tet_down_deltas_pfix[i][1]);
                vx_down_delta_h3_sum.push_back(tet_down_deltas_pfix[i][2]);
                vx_down_delta_h4_sum.push_back(tet_down_deltas_pfix[i][3]);
            }

            std::cout << "up delta totals:" << std::endl;

            for (vtkm::Id i = 0; i < num_sweep_values; i++)
            {
                std::cout << i << ") del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " del(h3)=" << vx_delta_h3_sum[i] << " del(h4)=" << vx_delta_h4_sum[i] << std::endl;
            }

            std::cout << "down delta totals:" << std::endl;

            for (vtkm::Id i = 0; i < num_sweep_values; i++)
            {
                std::cout << i << ") del(h1)=" << vx_down_delta_h1_sum[i] << " del(h2)=" << vx_down_delta_h2_sum[i] << " del(h3)=" << vx_down_delta_h3_sum[i] << " del(h4)=" << vx_down_delta_h4_sum[i] << std::endl;
            }



/* THE FOLLOWING REQUIRES SWEEP ISOVALUE H */
            // ==================== \/ Step 3: Deriving middle slab quad vertices P, Q, R, S \/ ==================== //

//            // First we define control vectors for the middle quad within the middle slab ...
//            // ... from our new middle slab points E F G H from Step 2:
//            std::vector<PositionVector> vectorsBH;
//            std::vector<PositionVector> vectorsEG;
//            std::vector<PositionVector> vectorsFC;
//            std::vector<PositionVector> vectorsBC;

//            for (int i = 0; i < tetlistSorted.size(); i++)
//            {
//                vectorsBH.emplace_back(verticesB[i], verticesH[i]);
//                vectorsEG.emplace_back(verticesE[i], verticesG[i]);
//                vectorsFC.emplace_back(verticesF[i], verticesC[i]);
//                vectorsBC.emplace_back(verticesB[i], verticesC[i]);

//                // interpolation value will be same for all points E F G H:
//                double lerpEFGH = double(h - teth2s[i]) / double(teth3s[i] - teth2s[i]);

//            }

            // ==================== /\ Step 3: Deriving middle slab quad vertices P, Q, R, S /\ ==================== //
/* THE PAST REQUIRES SWEEP ISOVALUE H */
        }



        // ----------------------------------- PRE-PROCESS ----------------------------------- //


        std::cout << "// ----------------------------------- PRE-PROCESS ----------------------------------- //" << std::endl;


        // for the sweep, we will be using the pre-computed delta coefficients from the mesh

        // -------------------------------------- SWEEP  ------------------------------------- //


        std::cout << "// -------------------------------------- SWEEP  ------------------------------------- //" << std::endl;

        std::vector<double> coefficientweightList;
        coefficientweightList.resize(superparentsPortal.GetNumberOfValues());


        std::vector<double> delta_h1_partial_pfixsum; // = 0.0;
        std::vector<double> delta_h2_partial_pfixsum; // = 0.0;
        std::vector<double> delta_h3_partial_pfixsum; // = 0.0;
        std::vector<double> delta_h4_partial_pfixsum; // = 0.0;


        delta_h1_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());
        delta_h2_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());
        delta_h3_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());
        delta_h4_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());


        for(int i = 0; i < coefficientweightList.size(); i++)
        {
            coefficientweightList[i] = 0.0;
            delta_h1_partial_pfixsum[i] = 0.0;
            delta_h2_partial_pfixsum[i] = 0.0;
            delta_h3_partial_pfixsum[i] = 0.0;
            delta_h3_partial_pfixsum[i] = 0.0;
        }

        auto arcsPortal = contourTree.Arcs.ReadPortal();
        auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto supernodesPortal = contourTree.Supernodes.ReadPortal();

        vtkm::Id prevIndex = -1;
        vtkm::Id superNodeID = prevIndex;

        std::cout << "Num of supernodes: " << contourTree.Supernodes.GetNumberOfValues() << std::endl;

        std::map<vtkm::Id, vtkm::Id> tailends;

        for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues(); supernode++)
        {
            vtkm::Id superNode = supernodesPortal.Get(supernode);

            std::cout << supernode << " - " << superNode << "->" << supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode))) << std::endl;

            tailends.insert(std::make_pair(superNode, supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode)))));
        }

        std::cout << "-----------------------------" << std::endl;

        if(dim1)
        {// test 1D coefficients
            for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
            { // per node in sorted order
              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              if (prevIndex != superparent)
              {
                  prevIndex = superparent;
                  // tail-end of a branch
                  superNodeID = sortID;
              }

              std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;


              if (sortedNode == 0)
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);
              else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);


              // CHANGES:
              // UPDATE AT REGULAR NODE: +1
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

              // UPDATE AT REGULAR NODE: +area/3
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);

              // weights before 2024-08-27:
              // superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);

              // weights after 2024-08-27:
              delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[sortID];
              delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[sortID];

              coefficientweightList[superparent] += delta_h1_partial_pfixsum[superparent] * sortID + delta_h2_partial_pfixsum[superparent];

              std::cout << "\t\t" << sortID << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                        << sortID << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << std::endl;


              if(sortedNode != contourTree.Arcs.GetNumberOfValues()-1)
              {
                  vtkm::Id nextSortID = nodesPortal.Get(sortedNode+1);
                  vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

                  if(nextSuperparent != superparent)
                  {
                      // std::cout << nextSuperparent << " -vs- " << superparent <<  " -- TRIGGER\n" << std::endl;

                      delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[tailends[superNodeID]];
                      delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[tailends[superNodeID]];

                      std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;

                      coefficientweightList[superparent] += delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                      std::cout << "\t\t" << tailends[superNodeID] << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                                << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;

                  }

              }
              else
              {
                // std::cout << "COMPUTE " << superparent << "->" << tailends[superNodeID] << std::endl;
                std::cout << tailends[superNodeID] << " - " << superparent << "->" << tailends[superNodeID] << "\n";
                coefficientweightList[superparent] += delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                std::cout << "\t\t" << sortID << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                          << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;


              }

              superarcIntrinsicWeightPortal.Set(superparent,
                                                superarcIntrinsicWeightPortal.Get(superparent)+coefficientweightList[superparent]);

            } // per node in sorted order
        }

        else if (dim2)
        {// test 2D coefficients

            double a_mid[] = {0.2148571, 0.0625, 0.16667, 0.25, 0.4375, 0.14285, 0.16667, 0.375};

            for(int i = 0; i < 8; i++)
            {
                std::cout << a_mid[i] << " ";
            }
            std::cout << std::endl;

            for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
            { // per node in sorted order
              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              if (prevIndex != superparent)
              {
                  prevIndex = superparent;
                  // tail-end of a branch
                  superNodeID = sortID;
              }

              std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;


              if (sortedNode == 0)
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);
              else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);


              // CHANGES:
              // UPDATE AT REGULAR NODE: +1
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

              // UPDATE AT REGULAR NODE: +area/3
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);

              // weights before 2024-08-27:
              // superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);

              // weights after 2024-08-27:
              delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[sortID];
              delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[sortID];
              delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[sortID];

              // 1/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
              coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (sortID*sortID) - delta_h2_partial_pfixsum[superparent] * sortID + delta_h3_partial_pfixsum[superparent];

              std::cout << "\t\t" << sortID << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                        << sortID << "^2 - " << delta_h2_partial_pfixsum[superparent] << " * " << sortID << " + " << delta_h3_partial_pfixsum[superparent]  << " = " << coefficientweightList[superparent] << std::endl;


              if(sortedNode != contourTree.Arcs.GetNumberOfValues()-1)
              {
                  vtkm::Id nextSortID = nodesPortal.Get(sortedNode+1);
                  vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

                  if(nextSuperparent != superparent)
                  {
                      // std::cout << nextSuperparent << " -vs- " << superparent <<  " -- TRIGGER\n" << std::endl;

                      delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[tailends[superNodeID]];
                      delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[tailends[superNodeID]];
                      delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[tailends[superNodeID]];

                      std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;

                      // 2/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
                      coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]) - delta_h2_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h3_partial_pfixsum[superparent];
                              //delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                      std::cout << "\t\t" << tailends[superNodeID]  << " -h " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                                << tailends[superNodeID]  << "^2 - " << delta_h2_partial_pfixsum[superparent] << " * " << tailends[superNodeID]  << " + " << delta_h3_partial_pfixsum[superparent]  << " = " << coefficientweightList[superparent] << std::endl;

                  }

              }
              else
              {
                // std::cout << "COMPUTE " << superparent << "->" << tailends[superNodeID] << std::endl;
                std::cout << tailends[superNodeID] << " - " << superparent << "->" << tailends[superNodeID] << "\n";
                // 3/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
                coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                std::cout << "\t\t" << sortID << " -g " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                          << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;


              }

              superarcIntrinsicWeightPortal.Set(superparent,
                                                superarcIntrinsicWeightPortal.Get(superparent)+coefficientweightList[superparent]);

            } // per node in sorted order
        }

        else if (dim3)
        {// test 3D coefficients
            std::cout << "3D Coefficients Sweep" << std::endl;
            for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
            { // per node in sorted order
              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              if (prevIndex != superparent)
              {
                  prevIndex = superparent;
                  // tail-end of a branch
                  superNodeID = sortID;
              }

              std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;


              if (sortedNode == 0)
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);
              else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);


              // CHANGES:
              // UPDATE AT REGULAR NODE: +1
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

              // UPDATE AT REGULAR NODE: +area/3
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);

              // weights before 2024-08-27:
              // superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);

              // weights after 2024-08-27:
              delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[sortID];
              delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[sortID];
              delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[sortID];
              delta_h4_partial_pfixsum[superparent] += vx_delta_h4_sum[sortID];

              // 1/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
              coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (sortID*sortID*sortID) +\
                                                   delta_h2_partial_pfixsum[superparent] * (sortID*sortID) +\
                                                   delta_h3_partial_pfixsum[superparent] * sortID + \
                                                   delta_h4_partial_pfixsum[superparent];

              std::cout << "\t\t" << sortID << " - " << superparent << "\t"
                        << delta_h1_partial_pfixsum[superparent] << " * " << sortID << "^3 + "
                        << delta_h2_partial_pfixsum[superparent] << " * " << sortID << "^2 + "
                        << delta_h3_partial_pfixsum[superparent] << " * " << sortID << " + "
                        << delta_h4_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << std::endl;


              if(sortedNode != contourTree.Arcs.GetNumberOfValues()-1)
              {
                  vtkm::Id nextSortID = nodesPortal.Get(sortedNode+1);
                  vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

                  if(nextSuperparent != superparent)
                  {
                      // std::cout << nextSuperparent << " -vs- " << superparent <<  " -- TRIGGER\n" << std::endl;

                      delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[tailends[superNodeID]];
                      delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[tailends[superNodeID]];
                      delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[tailends[superNodeID]];
                      delta_h4_partial_pfixsum[superparent] += vx_delta_h4_sum[tailends[superNodeID]];
                      std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;

                      // 2/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
//                      coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]) - delta_h2_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h3_partial_pfixsum[superparent];
                              //delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];
                      coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]*tailends[superNodeID]) +\
                                                           delta_h2_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]) +\
                                                           delta_h3_partial_pfixsum[superparent] * tailends[superNodeID] + \
                                                           delta_h4_partial_pfixsum[superparent];

//                      std::cout << "\t\t" << tailends[superNodeID]  << " -h " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
//                                << tailends[superNodeID]  << "^2 - " << delta_h2_partial_pfixsum[superparent] << " * " << tailends[superNodeID]  << " + " << delta_h3_partial_pfixsum[superparent]  << " = " << coefficientweightList[superparent] << std::endl;

                      std::cout << "\t\t" << sortID << " -h " << superparent << "\t"
                                << delta_h1_partial_pfixsum[superparent] << " * " << tailends[superNodeID] << "^3 + "
                                << delta_h2_partial_pfixsum[superparent] << " * " << tailends[superNodeID] << "^2 + "
                                << delta_h3_partial_pfixsum[superparent] << " * " << tailends[superNodeID] << " + "
                                << delta_h4_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << std::endl;


                  }

              }
              else
              {
                // std::cout << "COMPUTE " << superparent << "->" << tailends[superNodeID] << std::endl;
                std::cout << tailends[superNodeID] << " - " << superparent << "->" << tailends[superNodeID] << "\n";
                // 3/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
                coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                std::cout << "\t\t" << sortID << " -g " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                          << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;


              }

              superarcIntrinsicWeightPortal.Set(superparent,
                                                superarcIntrinsicWeightPortal.Get(superparent)+coefficientweightList[superparent]);

            } // per node in sorted order
        }


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

































    // 2024-11-25 COMPUTE THE STRUCT COEFFICIENTS VERSION OF THE WEIGHTS WITH COEFFICIENTS
    void static ComputeVolumeWeightsSerialStructCoefficients(const ContourTree& contourTree,
                                                const vtkm::Id nIterations,
                                                vtkm::cont::ArrayHandle<Coefficients>& superarcIntrinsicWeightCoeff,
                                                vtkm::cont::ArrayHandle<Coefficients>& superarcDependentWeightCoeff,
                                                vtkm::cont::ArrayHandle<Coefficients>& supernodeTransferWeightCoeff,
                                                vtkm::cont::ArrayHandle<Coefficients>& hyperarcDependentWeightCoeff)
    // 2) COEFFICIENTS:
    { // ContourTreeMaker::ComputeWeights()
      // START ComputeVolumeWeightsSerialStructCoefficients

        // Add old arrays that hold the actual volume for compatability:
        FloatArrayType superarcIntrinsicWeight;
        FloatArrayType superarcDependentWeight;
        FloatArrayType supernodeTransferWeight;
        FloatArrayType hyperarcDependentWeight;

        // start by storing the first sorted vertex ID for each superarc
        IdArrayType firstVertexForSuperparent;
        firstVertexForSuperparent.Allocate(contourTree.Superarcs.GetNumberOfValues());
        superarcIntrinsicWeight.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.WritePortal();

        // 2) COEFFICIENTS:
        superarcIntrinsicWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcIntrinsicWeightCoeffPortal = superarcIntrinsicWeightCoeff.WritePortal();

        auto firstVertexForSuperparentPortal = firstVertexForSuperparent.WritePortal();
        auto superparentsPortal = contourTree.Superparents.ReadPortal();
        auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
        auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
        auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
        // auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto nodesPortal = contourTree.Nodes.ReadPortal();
        // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();

        std::cout << "CALL FROM THE COEFFICIENT-BASED STRUCT FUNCTION" << std::endl;

        // Files for basic 1D experiments of Contour Tree Branch node-count weight computations
        // const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-coordinates.txt";
        // const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/BPECT-WW-16-triangles.txt";

        // Files for 2D experiments of Contour Tree Branch length/area-based weight computations
        // (For debugging, I am currently keeping both 2D and 3D files, as to quickly flick between them to compare correctness)
        const std::string filename1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-coordinates.txt";
        const std::string filename2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Square-9-triang.txt";

        // Files for 3D experiments of Contour Tree Branch volume-based weight computations
        // (For debugging, I am currently keeping both 2D and 3D files, as to quickly flick between them to compare correctness)
        const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Cube-8-coordinates.txt";
        const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/Cube-8-tets.txt";


        // get the total number of values to be swept ...
        // ... this will be one more than the total number of datapoints ...
        // ... because we include the region beyond the last isovalue N, as a range [N, +inf)
        // (also, the number of values corresponds to the total number of different data values in a dataset ...
        //  ... because of the simulation of simplicity, which ensures ever data point has a unique value)
        int num_sweep_values = contourTree.Arcs.GetNumberOfValues() + 1;

        // Coordinate lists and connectivities of 2D triangle-based and 3D tetrahedral datasets.
        // The coordinates are needed for computing the areas/volumes.
        // Note that the coordinate information is not previously needed for the base contour tree computation ...
        // ... and is only introduced as a requirement at this stage.
        // (For debugging, I am currently keeping both 2D and 3D files, as to quickly flick between them to compare correctness)
        std::vector<Coordinates> coordlist    = ReadCoordinatesFromFile(filename1);
        // triangle list has pairs of 3 vertices that make up the triangle
        std::vector<Triangle> trianglelist    = ReadTrianglesFromFile(filename2);

        std::vector<Coordinates> coordlist3D  = ReadCoordinatesFromFile(filename3D1);
        // tetrahedron list has pairs of 4 vertices that make up the tetrahedron
        std::vector<Tetrahedron> tetlist      = ReadTetsFromFile(filename3D2);

        // Weight list is used for the basic implementation that was used for my learning.
        // It is only used in the naive area weight implementation, ...
        // ... where the weight of each vertex of the triangle is area/3 (area div by 3)
        std::vector<double> weightList;

        std::cout << "PRINT THE ARRAYS OF COORDINATES: \n";

        // Print the coordinates data to check if it was read correctly.
        for (int i = 0; i < coordlist.size(); i++)
        {
            std::cout << i << ": " << coordlist[i].x << ", " << coordlist[i].y << ", " << coordlist[i].z << std::endl;
            // initialise the weight list array while at it
            weightList.push_back(0.0);
        }

        std::cout << "PRINT THE ARRAYS OF TRIANGLES: \n";
        std::cout << "num. of triangles: " << trianglelist.size() << std::endl;
        for (int i = 0; i < trianglelist.size(); i++)
        {
            std::cout << i << ": " << trianglelist[i].p1 << ", " << trianglelist[i].p2 << ", " << trianglelist[i].p3;
            double area = ComputeTriangleArea(
                            coordlist[trianglelist[i].p1].x, coordlist[trianglelist[i].p1].y, coordlist[trianglelist[i].p1].z,
                            coordlist[trianglelist[i].p2].x, coordlist[trianglelist[i].p2].y, coordlist[trianglelist[i].p2].z,
                            coordlist[trianglelist[i].p3].x, coordlist[trianglelist[i].p3].y, coordlist[trianglelist[i].p3].z);
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

        // Now move on to set up 3D tetrahedral variables

        std::cout << "PRINT THE ARRAYS OF TETS: \n";
        std::cout << "num. of tets: " << tetlist.size() << std::endl;

        // keep track of all tets in a sorted list
        // (the 'sort' refers to sorting by the vertex ID.
        // The vertex ID corresponds to a data value, ...
        //  ... for example a tet X, Y, Z, W might have vertices with IDs X=6, Y=5, Z=3, W=7)
        //  ... we sort them as 3, 5, 6, 7 ...
        //  ... and refer to vertices as A, B, C, D
        //  ... where we then assign A = 3, B = 5, C = 6, D = 7)
        //  (we do not keep track of the original ordering X, Y, Z, W) ...
        std::vector<std::vector<int>> tetlistSorted(tetlist.size(),
                                                    std::vector<int> (4, 0));

        // Vertices that define the Tetrahedron ABCD (Entire tetrahedron) ...
        // ... with their corresponding isovalues
        std::vector<vtkm::Vec3f_32> verticesA; // vertex A contains the lowest isovalue h1
        std::vector<int> teth1s;
        std::vector<vtkm::Vec3f_32> verticesB; // vertex B - h2
        std::vector<int> teth2s;
        std::vector<vtkm::Vec3f_32> verticesC; // vertex C - h3
        std::vector<int> teth3s;
        std::vector<vtkm::Vec3f_32> verticesD; // vertex D - h4
        std::vector<int> teth4s;

        // Deriving middle slab triangle vertices E, F, G, H
        // Plane Points at isovalue h=h2 (4) (for interval h1->h2)              - FIRST TET
        std::vector<vtkm::Vec3f_32> verticesE;
        std::vector<vtkm::Vec3f_32> verticesF;

        // Plane Points at isovalue h=h3 (5) (for interval h4->h3)              - LAST TET
        std::vector<vtkm::Vec3f_32> verticesG;
        std::vector<vtkm::Vec3f_32> verticesH;

        // Plane Points between isovalues h2 and h3  [Vertices P, Q, R, S]      - MIDDLE QUAD SLAB
        std::vector<vtkm::Vec3f_32> verticesP;
        std::vector<vtkm::Vec3f_32> verticesQ;
        std::vector<vtkm::Vec3f_32> verticesR;
        std::vector<vtkm::Vec3f_32> verticesS;

//        // Initializing the 2-D vector
//        std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
//                                              std::vector<double> (trianglelist.size(), 0.0));


        for (int i = 0; i < tetlist.size(); i++)
        {
            std::cout << i << ": " << tetlist[i].p1 << ", " << tetlist[i].p2 << ", " << tetlist[i].p3 << ", " << tetlist[i].p4; //<< std::endl;
            double volume = 0.1666666666667; // HARDCODED TODO CHANGE
            std::cout << " = " << volume << std::endl;

            tetlistSorted[i][0] = tetlist[i].p1;
            tetlistSorted[i][1] = tetlist[i].p2;
            tetlistSorted[i][2] = tetlist[i].p3;
            tetlistSorted[i][3] = tetlist[i].p4;
        }

//        std::sort()
        std::cout << "PRINT TETS UNSORTED:" << std::endl;
        print2DarrayInt(tetlistSorted);

        for (int i = 0; i < tetlist.size(); i++)
        {
            std::sort(tetlistSorted[i].begin(), tetlistSorted[i].end());
        }

        std::cout << "PRINT TETS SORTED:" << std::endl;
        print2DarrayInt(tetlistSorted);

        std::cout << "============================================================================================" << std::endl;

        //        std::sort(tetlist.begin(), tetlist.end(), std::greater<int>());

        bool dim1 = false;
        bool dim2 = false;
        bool dim3 = true;

        // ----------------------------------- PRE-PROCESS ----------------------------------- //
        // Here we basically populate the table as such:
        // vertexID / triangleID
        //          triangle1 triangle2   row sum(del. h1/h2)   col (prefix) sum del. h1/h2
        // deltas:  h1 h2     h1 h2
        // v1       x  y      z  w
        std::vector<double> vx_delta_h1;
        std::vector<double> vx_delta_h1_sum;
        std::vector<double> vx_down_delta_h1_sum;

        std::vector<double> vx_delta_h2;     // 1D coefficient deltas
        std::vector<double> vx_delta_h2_sum; // 1D coefficients
        std::vector<double> vx_down_delta_h2_sum; // 1D coefficients

        std::vector<double> vx_delta_h3;     // 2D coefficient deltas
        std::vector<double> vx_delta_h3_sum; // 2D coefficients
        std::vector<double> vx_down_delta_h3_sum; // 2D coefficients

        std::vector<double> vx_delta_h4;     // 3D coefficient deltas
        std::vector<double> vx_delta_h4_sum; // 3D coefficients
        std::vector<double> vx_down_delta_h4_sum; // 3D coefficients



        std::vector<double> delta_h1_pfixsum;
        std::vector<double> delta_h2_pfixsum; // 1D coefficients
        std::vector<double> delta_h3_pfixsum; // 2D coefficients
        std::vector<double> delta_h4_pfixsum; // 3D coefficients

        if(dim1)
        {
            // (below is actually the length of the hypotenuse projected onto the base from the 90-degree angle)
            // TODO: replace this with actual area of each triangle!
            double fake_area = sqrt(2.0) / 2.0;


            std::cout << "PREPROCESSING STEP:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));

            std::cout << "\nInitialised vxtch1:" << std::endl;
            print2Darray(vxtch1);

            std::vector<std::vector<double>> vxtch2(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));
            std::cout << "\nInitialised vxtch2:" << std::endl;
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
            std::cout << "vxtch1" << std::endl;
            print2Darray(vxtch1);
            std::cout << "vxtch2" << std::endl;
            print2Darray(vxtch2);


            // now with deltas computed, also compute the prefix sums of the different components:

            delta_h1_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
            delta_h2_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());

            std::cout << "deltas:" << std::endl;

            for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            {
                vx_delta_h1_sum.push_back(0.0);
                vx_delta_h2_sum.push_back(0.0);

                for (int j = 0; j < trianglelist.size(); j++)
                {
                    vx_delta_h1_sum[i] += vxtch1[i][j];
                    vx_delta_h2_sum[i] += vxtch2[i][j];
                }

                std::cout << i << ") del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " ~ "; //<< std::endl;

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

                std::cout << i << ") pfix(h1)= (" << i << ")" << delta_h1_pfixsum[i] << " pfix(h2)=" << delta_h2_pfixsum[i] << std::endl;
            }
        }
        else if(dim2)
        {
            std::cout << "Computing 2D coefficients:" << std::endl;

            // (below is the area of a full 1x1 90-degree triangle)
            // TODO: replace this with actual area of each triangle!
            double fake_area = 0.5; //sqrt(2.0) / 2.0;

            double a_mid[] = {0.2148571, 0.0625, 0.16667, 0.25, 0.4375, 0.14285, 0.16667, 0.375};

            std::cout << "Printing a_mid areas:" << std::endl;
            for(int i = 0; i < 8; i++)
            {
                std::cout << a_mid[i] << " ";
            }
            std::cout << std::endl;

            std::cout << "PREPROCESSING STEP:" << std::endl;

            std::cout << "// ----------------------------------- PRE-PROCESS ----------------------------------- //" << std::endl;


            // Initializing the 2-D vector
            std::vector<std::vector<double>> vxtch1(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));

            std::cout << "\nInitialised vxtch1:" << std::endl;
            print2Darray(vxtch1);

            std::vector<std::vector<double>> vxtch2(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));
            std::cout << "\nInitialised vxtch2:" << std::endl;
            print2Darray(vxtch2);

            std::vector<std::vector<double>> vxtch3(contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (trianglelist.size(), 0.0));
            print2Darray(vxtch3);
            std::cout << "\nInitialised vxtch3:" << std::endl;

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

                double Dn = Mv - Lv;
                double Dm = Hv - Mv;

                double Ba1 = 1.0/Dn;
                double Ba2 = Lv/Dn;

                double Bb1 = Hv / Dm;
                double Bb2 = 1.0/ Dm;


                // deltas:
                // Low-vertex deltas (same as original since adding from 0)
                double H12D = Ba1*Ba1 * a_mid[i];
                double H22D = 2 * Ba1 * Ba2 * a_mid[i];
                double H32D = Ba2*Ba2 * a_mid[i];

//                vxtch1[trianglelist[i].p3][i] = H12D;
//                vxtch2[trianglelist[i].p3][i] = H22D;
//                vxtch3[trianglelist[i].p3][i] = H32D;

                double H1b2D = -1.0 * (Bb2*Bb2) * (fake_area - a_mid[i]);
                double H2b2D = -2.0 * Bb2 * Bb1 * (fake_area - a_mid[i]);
                double H3b2D = fake_area - (Bb1 * Bb1) * (fake_area - a_mid[i]);

//                // Middle-vertex deltas:
//                vxtch1[trianglelist[i].p2][i] = H1b2D - H12D;
//                vxtch2[trianglelist[i].p2][i] = H2b2D - H22D;
//                vxtch3[trianglelist[i].p2][i] = H3b2D - H32D;

//                // High-vertex deltas:
//                vxtch1[trianglelist[i].p1][i] = -H1b2D;
//                vxtch2[trianglelist[i].p1][i] = -H2b2D;
//                vxtch3[trianglelist[i].p1][i] = -H3b2D;

                // invert operation:
//                // low:
//                vxtch1[trianglelist[i].p3][i] = -((H1b2D - H12D) -H1b2D); // = H12D
//                vxtch2[trianglelist[i].p3][i] = -((H2b2D - H22D) -H2b2D); // = H22D
//                vxtch3[trianglelist[i].p3][i] = -((H3b2D - H32D) -H3b2D); // = H32D

                // middle:
                vxtch1[trianglelist[i].p2][i] = -(-H1b2D + H12D);
                vxtch2[trianglelist[i].p2][i] = -(-H2b2D + H22D);
                vxtch3[trianglelist[i].p2][i] = -(-H3b2D + H32D);

                // high:
                vxtch1[trianglelist[i].p1][i] = -(H12D + (H1b2D - H12D));
                vxtch2[trianglelist[i].p1][i] = -(H22D + (H2b2D - H22D));
                vxtch3[trianglelist[i].p1][i] = fake_area-(H32D + (H3b2D - H32D));

                // low:
                vxtch1[trianglelist[i].p3][i] = -(vxtch1[trianglelist[i].p2][i] + vxtch1[trianglelist[i].p1][i]);  // -((H1b2D - H12D) -H1b2D); // = H12D
                vxtch2[trianglelist[i].p3][i] = -(vxtch2[trianglelist[i].p2][i] + vxtch2[trianglelist[i].p1][i]);  // -((H2b2D - H22D) -H2b2D); // = H22D
                vxtch3[trianglelist[i].p3][i] = fake_area-(vxtch3[trianglelist[i].p2][i] + vxtch3[trianglelist[i].p1][i]);  // fake_area-((H3b2D - H32D) -H3b2D); // = H32D

    //            std::cout << " " << i << ") " << delta_h1 << ", " << delta_h2 << std::endl;
            }


//            else if (dim4)
//            {

//            }

            std::cout << "vxtch1" << std::endl;
            print2Darray(vxtch1);
            std::cout << "vxtch2" << std::endl;
            print2Darray(vxtch2);
            std::cout << "vxtch3" << std::endl;
            print2Darray(vxtch3);


            // now with deltas computed, also compute the prefix sums of the different components:

//            std::vector<double> vx_delta_h1_sum;
//            std::vector<double> vx_delta_h2_sum;

//            std::vector<double> delta_h1_pfixsum; // = 0.0;
//            std::vector<double> delta_h2_pfixsum; // = 0.0;

            delta_h1_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
            delta_h2_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());
            delta_h3_pfixsum.resize(contourTree.Arcs.GetNumberOfValues());

            std::cout << "2D deltas:" << std::endl;

            for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            {
                vx_delta_h1_sum.push_back(0.0);
                vx_delta_h2_sum.push_back(0.0);
                vx_delta_h3_sum.push_back(0.0);

                for (int j = 0; j < trianglelist.size(); j++)
                {
                    vx_delta_h1_sum[i] += vxtch1[i][j];
                    vx_delta_h2_sum[i] += vxtch2[i][j];
                    vx_delta_h3_sum[i] += vxtch3[i][j];
                }

                std::cout << i << ") del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " del(h3)=" << vx_delta_h3_sum[i] << " ~ "; //<< std::endl;

                if (i == 0)
                {
                    delta_h1_pfixsum[0] = vx_delta_h1_sum[0];
    //                delta_h1_pfixsum.push_back(vx_delta_h1_sum[0]);
                    delta_h2_pfixsum[0] = vx_delta_h2_sum[0];
    //                delta_h2_pfixsum.push_back(vx_delta_h2_sum[0]);
                    delta_h3_pfixsum[0] = vx_delta_h3_sum[0];
                }
                else
                {
                    delta_h1_pfixsum[i] += delta_h1_pfixsum[i-1] + vx_delta_h1_sum[i];
                    delta_h2_pfixsum[i] += delta_h2_pfixsum[i-1] + vx_delta_h2_sum[i];
                    delta_h3_pfixsum[i] += delta_h3_pfixsum[i-1] + vx_delta_h3_sum[i];
                }

                std::cout << i << ") pfix(h1)= (" << i << ")" << delta_h1_pfixsum[i] << " pfix(h2)=" << delta_h2_pfixsum[i] << " pfix(h3)=" << delta_h3_pfixsum[i] << std::endl;
            }
        }
        else if(dim3)
        {
            std::cout << "Computing 3D Coefficients ..." << std::endl;

            // process one tet at a time ...
            // ... starting with they local coordinates and individual volumes
            for(int i = 0; i < tetlist.size(); i++)
            {
                std::cout << tetlist[i].p1 << " " << tetlist[i].p2 << " " << tetlist[i].p3 << " " << tetlist[i].p4 << std::endl;
                std::cout << "(" << tetlistSorted[i][0] << " " << tetlistSorted[i][1] << " " << tetlistSorted[i][2] << " " << tetlistSorted[i][3] << ")" << std::endl;

                std::cout << "\t " << tetlistSorted[i][0] << " = <" << coordlist3D[tetlistSorted[i][0]].x << " " << coordlist3D[tetlistSorted[i][0]].y << " " << coordlist3D[tetlistSorted[i][0]].z << ">" << std::endl;
                std::cout << "\t " << tetlistSorted[i][1] << " = <" << coordlist3D[tetlistSorted[i][1]].x << " " << coordlist3D[tetlistSorted[i][1]].y << " " << coordlist3D[tetlistSorted[i][1]].z << ">" << std::endl;
                std::cout << "\t " << tetlistSorted[i][2] << " = <" << coordlist3D[tetlistSorted[i][2]].x << " " << coordlist3D[tetlistSorted[i][2]].y << " " << coordlist3D[tetlistSorted[i][2]].z << ">" << std::endl;
                std::cout << "\t " << tetlistSorted[i][3] << " = <" << coordlist3D[tetlistSorted[i][3]].x << " " << coordlist3D[tetlistSorted[i][3]].y << " " << coordlist3D[tetlistSorted[i][3]].z << ">" << std::endl;
            }


            // ==================== \/ Step 1: Name the vertices that define each tetrahedron A, B, C, D \/ ==================== //
            // ... together with h1, h2, h3, h4
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                verticesA.emplace_back( coordlist3D[tetlistSorted[i][0]].x, coordlist3D[tetlistSorted[i][0]].y, coordlist3D[tetlistSorted[i][0]].z );
                verticesB.emplace_back( coordlist3D[tetlistSorted[i][1]].x, coordlist3D[tetlistSorted[i][1]].y, coordlist3D[tetlistSorted[i][1]].z );
                verticesC.emplace_back( coordlist3D[tetlistSorted[i][2]].x, coordlist3D[tetlistSorted[i][2]].y, coordlist3D[tetlistSorted[i][2]].z );
                verticesD.emplace_back( coordlist3D[tetlistSorted[i][3]].x, coordlist3D[tetlistSorted[i][3]].y, coordlist3D[tetlistSorted[i][3]].z );

                // keep track of tetrahedron boundary values h1, h2, h3, and h4
                teth1s.emplace_back(tetlistSorted[i][0]);
                teth2s.emplace_back(tetlistSorted[i][1]);
                teth3s.emplace_back(tetlistSorted[i][2]);
                teth4s.emplace_back(tetlistSorted[i][3]);
            }

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << "===" << std::endl;
                std::cout << "\t" << teth1s[i] << " = " << verticesA[i] << std::endl;
                std::cout << "\t" << teth2s[i] << " = " << verticesB[i] << std::endl;
                std::cout << "\t" << teth3s[i] << " = " << verticesC[i] << std::endl;
                std::cout << "\t" << teth4s[i] << " = " << verticesD[i] << std::endl;
                std::cout << "===" << std::endl;
            }

            // =================== \/ Step 1.5: from points A, B, C, D define vectors between them \/ ================== //

//            std::vector<vtkm::Vec3f_32> vectorsAB;
//            std::vector<vtkm::Vec3f_32> vectorsAC;

            std::vector<PositionVector> vectorsAB;
            std::vector<PositionVector> vectorsAC;
            std::vector<PositionVector> vectorsAD;
            std::vector<PositionVector> vectorsBD;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsAB.emplace_back(verticesA[i], verticesB[i]);
                vectorsAC.emplace_back(verticesA[i], verticesC[i]);
                vectorsAD.emplace_back(verticesA[i], verticesD[i]);
                vectorsBD.emplace_back(verticesB[i], verticesD[i]);

            }

            std::cout << "Vectors AB for each tet:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << i << " " << vectorsAD[i].start << " " << vectorsAD[i].end << " - " << vectorsAD[i].difference << std::endl;
                std::cout << i << " " << vectorsAC[i].start << " " << vectorsAC[i].end << " + " << vectorsAC[i].difference << std::endl;
//                vectorsAB[i].lerp(0.5);
//                std::cout << i << " " << vectorsAB[i].start << " " << vectorsAB[i].end << " + " << vectorsAB[i].difference << std::endl;
            }


            // =================== /\ Step 1.5: from points A, B, C, D define vectors between them /\ ================== //

            // ==================== \/ Step 2: Deriving middle slab triangle vertices E, F, G, H \/ ==================== //

            // Plane Points at isovalue h=h2 (4) (for interval h1->h2)              - FIRST TET

            // We will now be computing:
            //            verticesE;
            //            verticesF;
            // ... from vectors AD and AC respectively.


            for(int i = 0; i < tetlistSorted.size(); i++)
            {
                // lerp at h=2 on AD between 1 and 4 (a=1, d=4)
                double lerpADh2_h1_h4 = double(teth2s[i] - teth1s[i]) / double(teth4s[i] - teth1s[i]);
                verticesE.push_back(vectorsAD[i].lerp2point(lerpADh2_h1_h4));

                double lerpACh2_h1_h3 = double(teth2s[i] - teth1s[i]) / double(teth3s[i] - teth1s[i]);
                verticesF.push_back(vectorsAC[i].lerp2point(lerpACh2_h1_h3));

                std::cout << "lerpADh2_h1_h4: " << lerpADh2_h1_h4 << " E = " << verticesE[i] << std::endl;
                std::cout << "lerpACh2_h1_h3: " << lerpACh2_h1_h3 << " F = " << verticesF[i] << std::endl;
            }


            // Plane Points at isovalue h=h3 (5) (for interval h4->h3)              - LAST TET SLAB
            // We will now be computing:
            //            verticesG;
            //            verticesH;
            // ... from vectors AD and BD respectively.

            for(int i = 0; i < tetlistSorted.size(); i++)
            {
                // lerp at h=3 on AD between h2 and h3 (a=h1, d=h4)
                double lerpADh3_h1_h4 = double(teth3s[i] - teth1s[i]) / double(teth4s[i] - teth1s[i]);
                verticesG.push_back(vectorsAD[i].lerp2point(lerpADh3_h1_h4));

                // lerp at h3 on BD between h2 and h4 (b=h2, d=h4)

                double lerpBDh2_h2_h4 = double(teth3s[i] - teth2s[i]) / double(teth4s[i] - teth2s[i]);
                verticesH.push_back(vectorsBD[i].lerp2point(lerpBDh2_h2_h4));

                std::cout << "lerpADh3_h1_h4: " << lerpADh3_h1_h4 << " G = " << verticesG[i] << std::endl;
                std::cout << "lerpBDh2_h2_h4: " << lerpBDh2_h2_h4 << " H = " << verticesH[i] << std::endl;
            }


            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << "===" << std::endl;
                std::cout << "\t" << teth1s[i] << " = " << verticesA[i] << std::endl;
                std::cout << "\t" << teth2s[i] << " = " << verticesB[i] << std::endl;
                std::cout << "\t" << teth3s[i] << " = " << verticesC[i] << std::endl;
                std::cout << "\t" << teth4s[i] << " = " << verticesD[i] << std::endl;
                std::cout << "===" << std::endl;
            }
            // ==================== /\ Step 2: Deriving middle slab triangle vertices E, F, G, H /\ ==================== //




            // =========================== \/ Step 3: Compute Entire (full) Tet Volumes \/ ============================ //

            std::vector<double> full_tet_volumes;

            std::cout << "TET VOLUMES:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
//                vtkm::Vec3f_32 a_vol = verticesA[i];
//                vtkm::Vec3f_32 b_vol = verticesB[i];
//                vtkm::Vec3f_32 c_vol = verticesC[i];

                PositionVector a_vol(verticesA[i], verticesC[i]);
                PositionVector b_vol(verticesA[i], verticesD[i]);
                PositionVector c_vol(verticesA[i], verticesB[i]);

                full_tet_volumes.push_back((1.0/6.0) * abs(vtkm::Dot(vtkm::Cross(a_vol.difference, b_vol.difference), c_vol.difference) ) );
                std::cout << i << " = " << full_tet_volumes[i] << std::endl;
            }



            // =========================== /\ Step 3: Compute Entire (full) Tet Volumes  /\ ============================ //



            // ---------------------------------------------- FIRST SLAB ----------------------------------------------- //



            // ============== \/ Step 4: Compute the first slab volume (defined from isovalues h1-h2) \/ =============== //

            std::vector<double> slab1_h1h2_tet_volumes;

            std::cout << "SLAB1 VOLUMES:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
//                vtkm::Vec3f_32 a_vol = verticesA[i];
//                vtkm::Vec3f_32 b_vol = verticesB[i];
//                vtkm::Vec3f_32 c_vol = verticesC[i];
                double alpha_h2= std::max(0.0, std::min(1.0,
                                          double(teth2s[i]-teth1s[i])/double(teth3s[i]-teth1s[i])));

                double beta_h2 = std::max(0.0, std::min(1.0,
                                          double(teth2s[i]-teth1s[i])/double(teth4s[i]-teth1s[i])));

                // NOTE: gamma can be optimised away, since we reach point B at h2 (because B holds h2) ...
                // ... gamma will always be 1.0
                double gamma_h2= std::max(0.0, std::min(1.0,
                                          double(teth2s[i]-teth1s[i])/double(teth2s[i]-teth1s[i])));


                PositionVector a_h1h2_vol(verticesA[i], verticesC[i]);
                PositionVector b_h1h2_vol(verticesA[i], verticesD[i]);
                PositionVector c_h1h2_vol(verticesA[i], verticesB[i]);

                a_h1h2_vol.lerp(alpha_h2);
                b_h1h2_vol.lerp(beta_h2);
                c_h1h2_vol.lerp(gamma_h2);

                slab1_h1h2_tet_volumes.push_back((1.0/6.0) * abs(vtkm::Dot(vtkm::Cross(a_h1h2_vol.difference, b_h1h2_vol.difference),
                                                                           c_h1h2_vol.difference) ) );
                std::cout << i << " = " << slab1_h1h2_tet_volumes[i] << "(" << alpha_h2 << ", " << beta_h2 << ", " << gamma_h2 << ")" << std::endl;
            }



            // ==========  \/ Step 5: Compute the first slab coefficients (defined from isovalues h1-h2) \/ ============ //

            //
            std::vector<double> a_h1h2;
            std::vector<double> b_h1h2;
            std::vector<double> c_h1h2;
            std::vector<double> d_h1h2;


            std::cout << "SLAB1 h1h2 coefficients:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                a_h1h2.push_back(slab1_h1h2_tet_volumes[i]                                  / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );
                b_h1h2.push_back(-(3.0 * slab1_h1h2_tet_volumes[i] * teth1s[i])              / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );
                c_h1h2.push_back((3.0 * slab1_h1h2_tet_volumes[i] * std::pow(teth1s[i], 2)) / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );
                d_h1h2.push_back(-(slab1_h1h2_tet_volumes[i] * std::pow(teth1s[i], 3))      / double( std::pow((-teth1s[i] + teth2s[i]), 3) ) );


                std::cout << i << " = " << "a = " << a_h1h2[i]<< ", b = " << b_h1h2[i]<< ", c = " << c_h1h2[i] << ", d = " << d_h1h2[i] << std::endl;
            }

            // ========== /\ Step 5: Compute the first slab coefficients (defined from isovalues h1-h2) /\ ============ //



            // ---------------------------------------------- LAST SLAB ----------------------------------------------- //




            // ============== \/ Step 6: Compute the last slab volume (defined from isovalues h3-h4) \/ =============== //


            std::vector<double> slab3_h3h4_tet_volumes; // only the slab volume
            std::vector<double> slab3_h3h4_tet_volumes_sweeping_up; // (total volume) - (slab)

            std::cout << "SLAB3 VOLUMES:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
//                vtkm::Vec3f_32 a_vol = verticesA[i];
//                vtkm::Vec3f_32 b_vol = verticesB[i];
//                vtkm::Vec3f_32 c_vol = verticesC[i];
                double alpha_h3= std::max(0.0, std::min(1.0,
                                          double(teth3s[i]-teth4s[i])/double(teth2s[i]-teth4s[i])));

                double beta_h3 = std::max(0.0, std::min(1.0,
                                          double(teth3s[i]-teth4s[i])/double(teth1s[i]-teth4s[i])));

                double gamma_h3= std::max(0.0, std::min(1.0,
                                          double(teth3s[i]-teth4s[i])/double(teth3s[i]-teth4s[i])));


                PositionVector a_h3h4_vol(verticesD[i], verticesB[i]);
                PositionVector b_h3h4_vol(verticesD[i], verticesA[i]);
                PositionVector c_h3h4_vol(verticesD[i], verticesC[i]);

                a_h3h4_vol.lerp(alpha_h3);
                b_h3h4_vol.lerp(beta_h3);
                c_h3h4_vol.lerp(gamma_h3);

                slab3_h3h4_tet_volumes.push_back((1.0/6.0) * abs(vtkm::Dot(vtkm::Cross(a_h3h4_vol.difference, b_h3h4_vol.difference),
                                                                           c_h3h4_vol.difference) ) );

                slab3_h3h4_tet_volumes_sweeping_up.push_back(full_tet_volumes[i] - slab3_h3h4_tet_volumes[i]);

                std::cout << i << " = " << slab3_h3h4_tet_volumes[i] << "(" << alpha_h3 << ", " << beta_h3 << ", " << gamma_h3 << ")" << std::endl;
            }

            // ============== /\ Step 6: Compute the last slab volume (defined from isovalues h3-h4) /\ =============== //

            // ==========  \/ Step 7: Compute the last slab coefficients (defined from isovalues h3-h4) \/ ============ //

            //
            std::vector<double> a_h3h4;
            std::vector<double> b_h3h4;
            std::vector<double> c_h3h4;
            std::vector<double> d_h3h4;
            std::vector<double> d_h3h4_down;

            // compute the coefficients using 'full_tet_volumes' as v_volumeh4moveup ...
            // ... and 'slab3_h3h4_tet_volumes' as v_volumeh3moveup


            std::cout << "SLAB3 h3h4 coefficients:" << std::endl;
//            double d_h3h4; // dealing with the fourth coefficient separately, as we need to move the last segment up to start from slab volume at h3
            double d_h3h3_to0;              //       updated 'd' coefficient that moves the slab function to start at y=0
            std::vector<double> d_h3h4_up;  // final updated 'd' coefficient that moves the slab function to start at y=volume_h3
            double vol_h3 = 0.0;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << slab3_h3h4_tet_volumes_sweeping_up[i] << " " << full_tet_volumes[i] << std::endl;
                a_h3h4.push_back((slab3_h3h4_tet_volumes_sweeping_up[i] - full_tet_volumes[i])                                                              / double( std::pow((teth3s[i] - teth4s[i]), 3) ) );
                b_h3h4.push_back((-3.0*teth4s[i] * slab3_h3h4_tet_volumes_sweeping_up[i] + 3.0*teth4s[i] * full_tet_volumes[i])                             / double( std::pow((teth3s[i] - teth4s[i]), 3) ) );
                c_h3h4.push_back(( 3.0*std::pow(teth4s[i],2) * slab3_h3h4_tet_volumes_sweeping_up[i] - 3.0*std::pow(teth4s[i],2) * full_tet_volumes[i])     / double( std::pow((teth3s[i] - teth4s[i]), 3) ) );

                // compute the base d coefficient, update it later to lift the function up:
                d_h3h4.push_back((-std::pow(teth4s[i],3) * slab3_h3h4_tet_volumes_sweeping_up[i] + std::pow(teth4s[i],3) * full_tet_volumes[i])                     / double( std::pow((teth3s[i] - teth4s[i]), 3)));


                vol_h3 = a_h3h4[i]*std::pow(teth3s[i], 3) + b_h3h4[i]*std::pow(teth3s[i], 2) + c_h3h4[i]*teth3s[i] + d_h3h4[i];

                d_h3h3_to0 = d_h3h4[i] - vol_h3;

                // deal with the fourth coefficient (d - the constant) separately, as it helps to move the function up/down
                d_h3h4_up.push_back( d_h3h3_to0 + slab3_h3h4_tet_volumes_sweeping_up[i] );


                std::cout << i << " = " << "a = " << a_h3h4[i]<< ", b = " << b_h3h4[i]<< ", c = " << c_h3h4[i] << ", d = " << d_h3h4_up[i] << std::endl;
            }

            // ========== /\ Step 7: Compute the last slab coefficients (defined from isovalues h3-h4) /\ ============ //



            // ---------------------------------------------- MID SLAB ----------------------------------------------- //


            // =========================  \/ Step 8: Compute sin(theta1) and sin(theta2) \/ ========================== //
            // sin theta 1 computation:

            std::vector<double> areas_CGH;
            std::vector<double> sin_theta_1s;

            std::vector<PositionVector> vectorsGH;
            std::vector<PositionVector> vectorsCH;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsGH.emplace_back(verticesG[i], verticesH[i]);
                vectorsCH.emplace_back(verticesC[i], verticesH[i]);
            }

            std::cout << "Areas CGH:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                areas_CGH.push_back(1.0/2.0 * vtkm::Magnitude(vtkm::Cross( vectorsGH[i].difference, vectorsCH[i].difference )) );
                sin_theta_1s.push_back(2.0 * areas_CGH[i]/ (vectorsGH[i].mag() * vectorsCH[i].mag()) );

                std::cout << "Area of " << i << " = " << areas_CGH[i] << std::endl;
                std::cout << "sin(theta1) of " << i << " = " << sin_theta_1s[i] << std::endl;
            }



            // sin theta 2 computation:

            std::vector<double> areas_BEF;
            std::vector<double> sin_theta_2s;

            std::vector<PositionVector> vectorsFB;
            std::vector<PositionVector> vectorsFE;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsFB.emplace_back(verticesF[i], verticesB[i]);
                vectorsFE.emplace_back(verticesF[i], verticesE[i]);
            }

            std::cout << "Areas BEF:" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                areas_BEF.push_back(1.0/2.0 * vtkm::Magnitude(vtkm::Cross( vectorsFB[i].difference, vectorsFE[i].difference )) );
                sin_theta_2s.push_back(2.0 * areas_BEF[i]/ (vectorsFB[i].mag() * vectorsFE[i].mag()) );

                std::cout << "Area of " << i << " = " << areas_BEF[i] << std::endl;
                std::cout << "sin(theta1) of " << i << " = " << sin_theta_2s[i] << std::endl;
            }


            // =========================  /\ Step 8: Compute sin(theta1) and sin(theta2) /\ ========================== //


            // ===========================  \/ Step 9: Compute Middle Slab Coefficients \/ =========================== //

            // ----------------------------------------------- BATCH 1 ----------------------------------------------- //
            // a_s1, b_s1, c_s1
            std::vector<PositionVector> vectorsHG;
            std::vector<PositionVector> vectorsBE;
            std::vector<PositionVector> vectorsHC;
            std::vector<PositionVector> vectorsCG;
            // FE - already defined
            // FB - already defined


            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                vectorsHG.emplace_back(verticesH[i], verticesG[i]);
                vectorsBE.emplace_back(verticesB[i], verticesE[i]);
                vectorsHC.emplace_back(verticesH[i], verticesC[i]);
                vectorsCG.emplace_back(verticesC[i], verticesG[i]);

            }


            std::vector<double> tetk2s;
            std::vector<double> a_s1;
            std::vector<double> b_s1;
            std::vector<double> c_s1;

            // local variables for simplifying notation:
            double n1;
            double n2;
            double n3;
            double n4;

            std::cout << "Mid slab pre-coefficients" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // noting down repeating terms as I am writing the code for the first time:
                //                (tetk2s[i] * (vectorsHG[i].mag() - vectorsBE[i].mag() )

                tetk2s.push_back( 1.0 / double(teth3s[i] - teth2s[i]) );

                n1 = (tetk2s[i] * (vectorsHG[i].mag() - vectorsBE[i].mag() ));
                n2 = (vectorsHC[i].mag() * tetk2s[i]);

                n3 = ( teth2s[i] * tetk2s[i] * vectorsHG[i].mag() * vectorsHC[i].mag() * tetk2s[i] );
                n4 = ( tetk2s[i] * teth3s[i] * vectorsBE[i].mag() * vectorsHC[i].mag() * tetk2s[i] );


                a_s1.push_back( n1 * n2 );
                b_s1.push_back( -( (n1 * n2 * teth2s[i]) + n3 - n4 ) );
                c_s1.push_back( n3 * teth2s[i] - n4 * teth2s[i] );

                std::cout << "[a_s1] " << i << " == " << a_s1[i] << std::endl;
                std::cout << "[b_s1] " << i << " == " << b_s1[i] << std::endl;
                std::cout << "[c_s1] " << i << " == " << c_s1[i] << std::endl;

            }

            // ----------------------------------------------- BATCH 2 ----------------------------------------------- //
            // a_s2, b_s2, c_s2

            // HG - already defined
            // BE - already defined
            // HC - already defined
            // CG - already defined
            // FE - already defined
            // FB - already defined

            std::vector<double> a_s2;
            std::vector<double> b_s2;
            std::vector<double> c_s2;

            // local variables for simplifying notation:
            double m1;
            double m2;
            double m3;
            double m4;

            std::cout << "Mid slab pre-coefficients" << std::endl;
            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                m1 = (tetk2s[i] * (vectorsCG[i].mag() - vectorsFE[i].mag() ));
                m2 = (-vectorsFB[i].mag() * tetk2s[i]);

                m3 = (-teth2s[i] * tetk2s[i] * vectorsCG[i].mag() * vectorsFB[i].mag() * tetk2s[i] );
                m4 = ( tetk2s[i] * teth3s[i] * vectorsFE[i].mag() * vectorsFB[i].mag() * tetk2s[i] );


                a_s2.push_back( m1 * m2 );
                b_s2.push_back( m1 * -m2 * teth3s[i] - m3 -m4);
                c_s2.push_back( m3 * teth3s[i] + m4 * teth3s[i] );

                std::cout << "[a_s2] " << i << " == " << a_s2[i] << std::endl;
                std::cout << "[b_s2] " << i << " == " << b_s2[i] << std::endl;
                std::cout << "[c_s2] " << i << " == " << c_s2[i] << std::endl;

            }

            // ------------------------------------------ COMBINE BATCH 1+2 ------------------------------------------- //
            std::vector<double> a_mid;
            std::vector<double> b_mid;
            std::vector<double> c_mid;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                a_mid.push_back(sin_theta_1s[i] / 2.0 * a_s1[i] + sin_theta_2s[i] / 2.0 * a_s2[i]);
                b_mid.push_back(sin_theta_1s[i] / 2.0 * b_s1[i] + sin_theta_2s[i] / 2.0 * b_s2[i]);
                c_mid.push_back(sin_theta_1s[i] / 2.0 * c_s1[i] + sin_theta_2s[i] / 2.0 * c_s2[i]);

                std::cout << "comb " << i << " " << a_mid[i] << " " << b_mid[i] << " " << c_mid[i] << std::endl;
            }


            // ---------------------------- Compute the Integration correction coefficient ---------------------------- //
            std::vector<vtkm::Vec3f_32> plane_normals;
            std::vector<double> plane_distances;
            vtkm::Vec3f_32 FExFB_cross_product;

            std::vector<double> correction_factor_nominators;


            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                FExFB_cross_product = vtkm::Cross(vectorsFE[i].difference, vectorsFB[i].difference);

                plane_normals.push_back( FExFB_cross_product / (vtkm::Magnitude(FExFB_cross_product)) );
                plane_distances.push_back( vtkm::Magnitude(vtkm::Dot(plane_normals[i], verticesB[i]) - vtkm::Dot(plane_normals[i], verticesH[i]) ) / (vtkm::Magnitude(plane_normals[i]) ) );
//                plane_distances.push_back( (vtkm::Dot(plane_normals[i], verticesB[i]) - vtkm::Dot(plane_normals[i], verticesH[i]) ) / (vtkm::Magnitude(plane_normals[i]) ) );


                correction_factor_nominators.push_back(plane_distances[i] * tetk2s[i]);

                std::cout << "corr: " << i << " " << correction_factor_nominators[i] << std::endl;

            }


            // ---------------------------- Compute the Integration correction coefficient ---------------------------- //

            double d_h2h3_to0;
//            double d_h2h3_up;
//            std::vector<double> d_mid;

            std::vector<double> a_h2h3;
            std::vector<double> b_h2h3;
            std::vector<double> c_h2h3;
            std::vector<double> d_h2h3;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                a_h2h3.push_back(correction_factor_nominators[i]/3.0 * a_mid[i]);
                b_h2h3.push_back(correction_factor_nominators[i]/2.0 * b_mid[i]);
                c_h2h3.push_back(correction_factor_nominators[i]     * c_mid[i]);

                d_h2h3_to0 = a_h2h3[i]  * std::pow(teth2s[i], 3) +\
                             b_h2h3[i]  * std::pow(teth2s[i], 2) +\
                             c_h2h3[i]  *          teth2s[i];

//                d_mid.push_back();


                d_h2h3.push_back(-d_h2h3_to0 + slab1_h1h2_tet_volumes[i]);


            }

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                std::cout << "mid coeffs " << i << " : " << a_h2h3[i] << " " << b_h2h3[i] << " " << c_h2h3[i] << " " << d_h2h3[i] << std::endl;
            }


            // ===========================  /\ Step 9: Compute Middle Slab Coefficients /\ ===========================  //





            // ============================  \/ Step 10: Compute UP coefficients table  \/ ============================  //

            std::cout << "Initialised coefficient tables:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_coeffs_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_coeffs_h1:" << std::endl;
            print2Darray(tet_coeffs_h1);

            std::vector<std::vector<double>> tet_coeffs_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_coeffs_h2:" << std::endl;
            print2Darray(tet_coeffs_h2);

            std::vector<std::vector<double>> tet_coeffs_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_coeffs_h3:" << std::endl;
            print2Darray(tet_coeffs_h3);

            std::vector<std::vector<double>> tet_coeffs_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_coeffs_h4:" << std::endl;
            print2Darray(tet_coeffs_h4);



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in just the h1 coeffs:
                for(int h = teth1s[i]; h < teth2s[i]; h++)
                {
                    tet_coeffs_h1[h][i] = a_h1h2[i];
                    tet_coeffs_h2[h][i] = b_h1h2[i];
                    tet_coeffs_h3[h][i] = c_h1h2[i];
                    tet_coeffs_h4[h][i] = d_h1h2[i];
                }
                for(int h = teth2s[i]; h < teth3s[i]; h++)
                {
                    tet_coeffs_h1[h][i] = a_h2h3[i];
                    tet_coeffs_h2[h][i] = b_h2h3[i];
                    tet_coeffs_h3[h][i] = c_h2h3[i];
                    tet_coeffs_h4[h][i] = d_h2h3[i];
                }
                for(int h = teth3s[i]; h < teth4s[i]; h++)
                {
                    tet_coeffs_h1[h][i] = a_h3h4[i];
                    tet_coeffs_h2[h][i] = b_h3h4[i];
                    tet_coeffs_h3[h][i] = c_h3h4[i];
                    tet_coeffs_h4[h][i] = d_h3h4_up[i];
                }
                for(int h = teth4s[i]; h < 8; h++)
                {
                    tet_coeffs_h1[h][i] = 0.0;
                    tet_coeffs_h2[h][i] = 0.0;
                    tet_coeffs_h3[h][i] = 0.0;
                    tet_coeffs_h4[h][i] = full_tet_volumes[i];
                }
            }

            std::cout << "\nh1 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h1);
            std::cout << "\nh2 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h2);
            std::cout << "\nh3 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h3);
            std::cout << "\nh4 UP coefficients written to table:" << std::endl;
            print2Darray(tet_coeffs_h4);



            // ============================  /\ Step 10: Compute UP coefficients table  /\ ============================  //





            // =========================  \/ Step 11: Compute UP coefficient deltas table  \/ =========================  //

            std::cout << "Initialised coefficient deltas tables:" << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "  << std::endl;
            std::cout << "Contour Tree Number of Arcs: " << contourTree.Arcs.GetNumberOfValues() << std::endl;
            std::cout << "Total Sweep Values then +1 : " << num_sweep_values << std::endl; // 2024-11-18 updated sweep value
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "  << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_deltas_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_deltas_h1:" << std::endl;
            print2Darray(tet_deltas_h1);

            std::vector<std::vector<double>> tet_deltas_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_deltas_h2:" << std::endl;
            print2Darray(tet_deltas_h2);

            std::vector<std::vector<double>> tet_deltas_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_deltas_h3:" << std::endl;
            print2Darray(tet_deltas_h3);

            std::vector<std::vector<double>> tet_deltas_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_deltas_h4:" << std::endl;
            print2Darray(tet_deltas_h4 );



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in just the h1 coeffs:
//                for(int h = teth1s[i]; h < teth2s[i]; h++)
//                {
                    tet_deltas_h1[teth1s[i]+1][i] = a_h1h2[i] - 0.0;
                    tet_deltas_h2[teth1s[i]+1][i] = b_h1h2[i] - 0.0;
                    tet_deltas_h3[teth1s[i]+1][i] = c_h1h2[i] - 0.0;
                    tet_deltas_h4[teth1s[i]+1][i] = d_h1h2[i] - 0.0;
//                }
//                for(int h = teth2s[i]; h < teth3s[i]; h++)
//                {
                    tet_deltas_h1[teth2s[i]+1][i] = a_h2h3[i] - a_h1h2[i];
                    tet_deltas_h2[teth2s[i]+1][i] = b_h2h3[i] - b_h1h2[i];
                    tet_deltas_h3[teth2s[i]+1][i] = c_h2h3[i] - c_h1h2[i];
                    tet_deltas_h4[teth2s[i]+1][i] = d_h2h3[i] - d_h1h2[i];
//                }
//                for(int h = teth3s[i]; h < teth4s[i]; h++)
//                {
                    tet_deltas_h1[teth3s[i]+1][i] = a_h3h4[i] - a_h2h3[i];
                    tet_deltas_h2[teth3s[i]+1][i] = b_h3h4[i] - b_h2h3[i];
                    tet_deltas_h3[teth3s[i]+1][i] = c_h3h4[i] - c_h2h3[i];
                    tet_deltas_h4[teth3s[i]+1][i] = d_h3h4_up[i] - d_h2h3[i];
//                }
//                for(int h = teth4s[i]; h < 8; h++)
//                {
                    tet_deltas_h1[teth4s[i]+1][i] = - a_h3h4[i];
                    tet_deltas_h2[teth4s[i]+1][i] = - b_h3h4[i];
                    tet_deltas_h3[teth4s[i]+1][i] = - c_h3h4[i];
                    tet_deltas_h4[teth4s[i]+1][i] = full_tet_volumes[i] - d_h3h4_up[i];
//                }
            }

            std::cout << "\nh1 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h1);
            std::cout << "\nh2 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h2);
            std::cout << "\nh3 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h3);
            std::cout << "\nh4 UP coefficient deltas written to table:" << std::endl;
            print2Darray(tet_deltas_h4);



            // =========================  /\ Step 11: Compute UP coefficient deltas table  /\ =========================  //






            // ============================ \/ Step 12: Compute DOWN coefficient tables  \/ ==========================  //

            std::vector<double> d_h1h2_down;
            std::vector<double> d_h2h3_down;

            for (int i = 0; i < tetlistSorted.size(); i++)
            {
                d_h2h3_down.push_back( (  -a_h2h3[i] * std::pow(teth3s[i], 3) - b_h2h3[i] * std::pow(teth3s[i], 2) - c_h2h3[i] * teth3s[i] ) \
                                       - (-a_h3h4[i] * std::pow(teth3s[i], 3) - b_h3h4[i] * std::pow(teth3s[i], 2) - c_h3h4[i] * teth3s[i] - d_h3h4[i] ) );

                d_h1h2_down.push_back( (  -a_h1h2[i] * std::pow(teth2s[i], 3) - b_h1h2[i] * std::pow(teth2s[i], 2) - c_h1h2[i] * teth2s[i] ) \
                                       - (-a_h2h3[i] * std::pow(teth2s[i], 3) - b_h2h3[i] * std::pow(teth2s[i], 2) - c_h2h3[i] * teth2s[i] - d_h2h3_down[i] ) );

            }







            std::cout << "Initialised DOWN coefficient tables:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_down_coeffs_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_down_coeffs_h1:" << std::endl;
            print2Darray(tet_down_coeffs_h1);

            std::vector<std::vector<double>> tet_down_coeffs_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_coeffs_h2:" << std::endl;
            print2Darray(tet_down_coeffs_h2);

            std::vector<std::vector<double>> tet_down_coeffs_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_coeffs_h3:" << std::endl;
            print2Darray(tet_down_coeffs_h3);

            std::vector<std::vector<double>> tet_down_coeffs_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_coeffs_h4:" << std::endl;
            print2Darray(tet_down_coeffs_h4);



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in before the h1 coeffs:
                for(int h = 0; h < teth1s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = 0.0;
                    tet_down_coeffs_h2[h][i] = 0.0;
                    tet_down_coeffs_h3[h][i] = 0.0;
                    tet_down_coeffs_h4[h][i] = full_tet_volumes[i];
                }

                // fill in just the h1 coeffs:
                for(int h = teth1s[i]; h < teth2s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = -a_h1h2[i];
                    tet_down_coeffs_h2[h][i] = -b_h1h2[i];
                    tet_down_coeffs_h3[h][i] = -c_h1h2[i];
                    tet_down_coeffs_h4[h][i] = -d_h1h2_down[i];
                }
                for(int h = teth2s[i]; h < teth3s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = -a_h2h3[i];
                    tet_down_coeffs_h2[h][i] = -b_h2h3[i];
                    tet_down_coeffs_h3[h][i] = -c_h2h3[i];
                    tet_down_coeffs_h4[h][i] = -d_h2h3_down[i];
                }
                for(int h = teth3s[i]; h < teth4s[i]; h++)
                {
                    tet_down_coeffs_h1[h][i] = -a_h3h4[i];
                    tet_down_coeffs_h2[h][i] = -b_h3h4[i];
                    tet_down_coeffs_h3[h][i] = -c_h3h4[i];
                    tet_down_coeffs_h4[h][i] = -d_h3h4[i];
                }
                for(int h = teth4s[i]; h < 8; h++)
                {
                    tet_down_coeffs_h1[h][i] = 0.0;
                    tet_down_coeffs_h2[h][i] = 0.0;
                    tet_down_coeffs_h3[h][i] = 0.0;
                    tet_down_coeffs_h4[h][i] = 0.0;
                }
            }

            std::cout << "\nh1 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h1);
            std::cout << "\nh2 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h2);
            std::cout << "\nh3 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h3);
            std::cout << "\nh4 DOWN coefficients written to table:" << std::endl;
            print2Darray(tet_down_coeffs_h4);



            // ============================ /\ Step 12: Compute DOWN coefficient tables  /\ ==========================  //








            // ======================== \/ Step 13: Compute DOWN coefficient delta tables  \/ =======================  //

            std::cout << "Initialised DOWN coefficient deltas tables:" << std::endl;

            // Initializing the 2-D vector
            std::vector<std::vector<double>> tet_down_deltas_h1(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));

            std::cout << "\nInitialised tet_down_deltas_h1:" << std::endl;
            print2Darray(tet_down_deltas_h1);

            std::vector<std::vector<double>> tet_down_deltas_h2(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_deltas_h2:" << std::endl;
            print2Darray(tet_down_deltas_h2);

            std::vector<std::vector<double>> tet_down_deltas_h3(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_deltas_h3:" << std::endl;
            print2Darray(tet_down_deltas_h3);

            std::vector<std::vector<double>> tet_down_deltas_h4(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (tetlistSorted.size(), 0.0));
            std::cout << "\nInitialised tet_down_deltas_h4:" << std::endl;
            print2Darray(tet_down_deltas_h4 );



            for (int i = 0; i < tetlistSorted.size(); i++)
            {

                // fill in just the h1 coeffs:
//                for(int h = teth1s[i]; h < teth2s[i]; h++)
//                {
                    tet_down_deltas_h1[teth1s[i]][i] = 0.0 + a_h1h2[i];
                    tet_down_deltas_h2[teth1s[i]][i] = 0.0 + b_h1h2[i];
                    tet_down_deltas_h3[teth1s[i]][i] = 0.0 + c_h1h2[i];
                    tet_down_deltas_h4[teth1s[i]][i] = full_tet_volumes[i] + d_h1h2_down[i];
//                }
//                for(int h = teth2s[i]; h < teth3s[i]; h++)
//                {
                    tet_down_deltas_h1[teth2s[i]][i] = -a_h1h2[i] + a_h2h3[i];
                    tet_down_deltas_h2[teth2s[i]][i] = -b_h1h2[i] + b_h2h3[i];
                    tet_down_deltas_h3[teth2s[i]][i] = -c_h1h2[i] + c_h2h3[i];
                    tet_down_deltas_h4[teth2s[i]][i] = -d_h1h2_down[i] + d_h2h3_down[i];
//                }
//                for(int h = teth3s[i]; h < teth4s[i]; h++)
//                {
                    tet_down_deltas_h1[teth3s[i]][i] = -a_h2h3[i] + a_h3h4[i];
                    tet_down_deltas_h2[teth3s[i]][i] = -b_h2h3[i] + b_h3h4[i];
                    tet_down_deltas_h3[teth3s[i]][i] = -c_h2h3[i] + c_h3h4[i];
                    tet_down_deltas_h4[teth3s[i]][i] = -d_h2h3_down[i] + d_h3h4[i];
//                }
//                for(int h = teth4s[i]; h < 8; h++)
//                {
                    tet_down_deltas_h1[teth4s[i]][i] = -a_h3h4[i];
                    tet_down_deltas_h2[teth4s[i]][i] = -b_h3h4[i];
                    tet_down_deltas_h3[teth4s[i]][i] = -c_h3h4[i];
                    tet_down_deltas_h4[teth4s[i]][i] = -d_h3h4[i];
//                }
            }

            std::cout << "\nh1 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h1);
            std::cout << "\nh2 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h2);
            std::cout << "\nh3 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h3);
            std::cout << "\nh4 DOWN coefficient deltas written to table:" << std::endl;
            print2Darray(tet_down_deltas_h4);


            // ======================== /\ Step 13: Compute DOWN coefficient delta tables  /\ =======================  //



            // ======================= \/ Step 14: Compute delta table up/down prefix sums \/ =======================  //

            std::vector<std::vector<double>> tet_up_deltas_pfix(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (4, 0.0));


            std::vector<std::vector<double>> tet_down_deltas_pfix(num_sweep_values, // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues(),
                                                  std::vector<double> (4, 0.0));

            for (int i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            // num_sweep_values, // 2024-11-18 updated sweep value
            for (int i = 0; i < num_sweep_values; i++)
            {
                double a_coeff_sum = 0.0;
                double b_coeff_sum = 0.0;
                double c_coeff_sum = 0.0;
                double d_coeff_sum = 0.0;

                // up deltas
                for(int k = 0; k < tetlistSorted.size(); k++)
                {
                    a_coeff_sum += tet_deltas_h1[i][k];
                    b_coeff_sum += tet_deltas_h2[i][k];
                    c_coeff_sum += tet_deltas_h3[i][k];
                    d_coeff_sum += tet_deltas_h4[i][k];
                }
                tet_up_deltas_pfix[i][0] =  a_coeff_sum;
                tet_up_deltas_pfix[i][1] =  b_coeff_sum;
                tet_up_deltas_pfix[i][2] =  c_coeff_sum;
                tet_up_deltas_pfix[i][3] =  d_coeff_sum;


                a_coeff_sum = 0.0;
                b_coeff_sum = 0.0;
                c_coeff_sum = 0.0;
                d_coeff_sum = 0.0;

                // down deltas
                for(int k = 0; k < tetlistSorted.size(); k++)
                {
                    a_coeff_sum += tet_down_deltas_h1[i][k];
                    b_coeff_sum += tet_down_deltas_h2[i][k];
                    c_coeff_sum += tet_down_deltas_h3[i][k];
                    d_coeff_sum += tet_down_deltas_h4[i][k];
                }
                tet_down_deltas_pfix[i][0] =  a_coeff_sum;
                tet_down_deltas_pfix[i][1] =  b_coeff_sum;
                tet_down_deltas_pfix[i][2] =  c_coeff_sum;
                tet_down_deltas_pfix[i][3] =  d_coeff_sum;

            }

            std::cout << "tet_up_deltas_pfix:" << std::endl;
            print2Darray(tet_up_deltas_pfix);

            std::cout << "tet_down_deltas_pfix:" << std::endl;
            print2Darray(tet_down_deltas_pfix);

            // ======================= /\ Step 14: Compute delta table up/down prefix sums /\ =======================  //

            std::cout << std::endl;

//            for (vtkm::Id i = 0; i < contourTree.Arcs.GetNumberOfValues(); i++)
            for (vtkm::Id i = 0; i < num_sweep_values; i++)
            {
//                vx_delta_h1_sum.push_back(tet_up_deltas_pfix[i][0]);
//                vx_delta_h2_sum.push_back(tet_up_deltas_pfix[i][1]);
//                vx_delta_h3_sum.push_back(tet_up_deltas_pfix[i][2]);
//                vx_delta_h4_sum.push_back(tet_up_deltas_pfix[i][3]);

                vx_delta_h1_sum.push_back(tet_down_deltas_pfix[i][0]);
                vx_delta_h2_sum.push_back(tet_down_deltas_pfix[i][1]);
                vx_delta_h3_sum.push_back(tet_down_deltas_pfix[i][2]);
                vx_delta_h4_sum.push_back(tet_down_deltas_pfix[i][3]);

                vx_down_delta_h1_sum.push_back(tet_down_deltas_pfix[i][0]);
                vx_down_delta_h2_sum.push_back(tet_down_deltas_pfix[i][1]);
                vx_down_delta_h3_sum.push_back(tet_down_deltas_pfix[i][2]);
                vx_down_delta_h4_sum.push_back(tet_down_deltas_pfix[i][3]);
            }

            std::cout << "up delta totals:" << std::endl;

            for (vtkm::Id i = 0; i < num_sweep_values; i++)
            {
                std::cout << i << ") del(h1)=" << vx_delta_h1_sum[i] << " del(h2)=" << vx_delta_h2_sum[i] << " del(h3)=" << vx_delta_h3_sum[i] << " del(h4)=" << vx_delta_h4_sum[i] << std::endl;
            }

            std::cout << "down delta totals:" << std::endl;

            for (vtkm::Id i = 0; i < num_sweep_values; i++)
            {
                std::cout << i << ") del(h1)=" << vx_down_delta_h1_sum[i] << " del(h2)=" << vx_down_delta_h2_sum[i] << " del(h3)=" << vx_down_delta_h3_sum[i] << " del(h4)=" << vx_down_delta_h4_sum[i] << std::endl;
            }



/* THE FOLLOWING REQUIRES SWEEP ISOVALUE H */
            // ==================== \/ Step 3: Deriving middle slab quad vertices P, Q, R, S \/ ==================== //

//            // First we define control vectors for the middle quad within the middle slab ...
//            // ... from our new middle slab points E F G H from Step 2:
//            std::vector<PositionVector> vectorsBH;
//            std::vector<PositionVector> vectorsEG;
//            std::vector<PositionVector> vectorsFC;
//            std::vector<PositionVector> vectorsBC;

//            for (int i = 0; i < tetlistSorted.size(); i++)
//            {
//                vectorsBH.emplace_back(verticesB[i], verticesH[i]);
//                vectorsEG.emplace_back(verticesE[i], verticesG[i]);
//                vectorsFC.emplace_back(verticesF[i], verticesC[i]);
//                vectorsBC.emplace_back(verticesB[i], verticesC[i]);

//                // interpolation value will be same for all points E F G H:
//                double lerpEFGH = double(h - teth2s[i]) / double(teth3s[i] - teth2s[i]);

//            }

            // ==================== /\ Step 3: Deriving middle slab quad vertices P, Q, R, S /\ ==================== //
/* THE PAST REQUIRES SWEEP ISOVALUE H */
        }



        // ----------------------------------- PRE-PROCESS ----------------------------------- //


        std::cout << "// ----------------------------------- PRE-PROCESS ----------------------------------- //" << std::endl;


        // for the sweep, we will be using the pre-computed delta coefficients from the mesh

        // -------------------------------------- SWEEP  ------------------------------------- //


        std::cout << "// -------------------------------------- SWEEP  ------------------------------------- //" << std::endl;

        std::vector<double> coefficientweightList;
        coefficientweightList.resize(superparentsPortal.GetNumberOfValues());


        std::vector<double> delta_h1_partial_pfixsum; // = 0.0;
        std::vector<double> delta_h2_partial_pfixsum; // = 0.0;
        std::vector<double> delta_h3_partial_pfixsum; // = 0.0;
        std::vector<double> delta_h4_partial_pfixsum; // = 0.0;


        delta_h1_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());
        delta_h2_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());
        delta_h3_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());
        delta_h4_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());


        for(int i = 0; i < coefficientweightList.size(); i++)
        {
            coefficientweightList[i] = 0.0;
            delta_h1_partial_pfixsum[i] = 0.0;
            delta_h2_partial_pfixsum[i] = 0.0;
            delta_h3_partial_pfixsum[i] = 0.0;
            delta_h3_partial_pfixsum[i] = 0.0;
        }

        auto arcsPortal = contourTree.Arcs.ReadPortal();
        auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto supernodesPortal = contourTree.Supernodes.ReadPortal();

        vtkm::Id prevIndex = -1;
        vtkm::Id superNodeID = prevIndex;

        std::cout << "Num of supernodes: " << contourTree.Supernodes.GetNumberOfValues() << std::endl;

        std::map<vtkm::Id, vtkm::Id> tailends;

        for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues(); supernode++)
        {
            vtkm::Id superNode = supernodesPortal.Get(supernode);

            std::cout << supernode << " - " << superNode << "->" << supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode))) << std::endl;

            tailends.insert(std::make_pair(superNode, supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode)))));
        }

        std::cout << "-----------------------------" << std::endl;

        if(dim1)
        {// test 1D coefficients
            for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
            { // per node in sorted order
              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              if (prevIndex != superparent)
              {
                  prevIndex = superparent;
                  // tail-end of a branch
                  superNodeID = sortID;
              }

              std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;


              if (sortedNode == 0)
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);
              else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);


              // CHANGES:
              // UPDATE AT REGULAR NODE: +1
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

              // UPDATE AT REGULAR NODE: +area/3
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);

              // weights before 2024-08-27:
              // superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);

              // weights after 2024-08-27:
              delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[sortID];
              delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[sortID];

              coefficientweightList[superparent] += delta_h1_partial_pfixsum[superparent] * sortID + delta_h2_partial_pfixsum[superparent];

              std::cout << "\t\t" << sortID << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                        << sortID << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << std::endl;


              if(sortedNode != contourTree.Arcs.GetNumberOfValues()-1)
              {
                  vtkm::Id nextSortID = nodesPortal.Get(sortedNode+1);
                  vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

                  if(nextSuperparent != superparent)
                  {
                      // std::cout << nextSuperparent << " -vs- " << superparent <<  " -- TRIGGER\n" << std::endl;

                      delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[tailends[superNodeID]];
                      delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[tailends[superNodeID]];

                      std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;

                      coefficientweightList[superparent] += delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                      std::cout << "\t\t" << tailends[superNodeID] << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                                << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;

                  }

              }
              else
              {
                // std::cout << "COMPUTE " << superparent << "->" << tailends[superNodeID] << std::endl;
                std::cout << tailends[superNodeID] << " - " << superparent << "->" << tailends[superNodeID] << "\n";
                coefficientweightList[superparent] += delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                std::cout << "\t\t" << sortID << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                          << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;


              }

              superarcIntrinsicWeightPortal.Set(superparent,
                                                superarcIntrinsicWeightPortal.Get(superparent)+coefficientweightList[superparent]);

            } // per node in sorted order
        }

        else if (dim2)
        {// test 2D coefficients

            double a_mid[] = {0.2148571, 0.0625, 0.16667, 0.25, 0.4375, 0.14285, 0.16667, 0.375};

            for(int i = 0; i < 8; i++)
            {
                std::cout << a_mid[i] << " ";
            }
            std::cout << std::endl;

            for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
            { // per node in sorted order
              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              if (prevIndex != superparent)
              {
                  prevIndex = superparent;
                  // tail-end of a branch
                  superNodeID = sortID;
              }

              std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;


              if (sortedNode == 0)
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);
              else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);


              // CHANGES:
              // UPDATE AT REGULAR NODE: +1
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

              // UPDATE AT REGULAR NODE: +area/3
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);

              // weights before 2024-08-27:
              // superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);

              // weights after 2024-08-27:
              delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[sortID];
              delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[sortID];
              delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[sortID];

              // 1/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
              coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (sortID*sortID) - delta_h2_partial_pfixsum[superparent] * sortID + delta_h3_partial_pfixsum[superparent];

              std::cout << "\t\t" << sortID << " - " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                        << sortID << "^2 - " << delta_h2_partial_pfixsum[superparent] << " * " << sortID << " + " << delta_h3_partial_pfixsum[superparent]  << " = " << coefficientweightList[superparent] << std::endl;


              if(sortedNode != contourTree.Arcs.GetNumberOfValues()-1)
              {
                  vtkm::Id nextSortID = nodesPortal.Get(sortedNode+1);
                  vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

                  if(nextSuperparent != superparent)
                  {
                      // std::cout << nextSuperparent << " -vs- " << superparent <<  " -- TRIGGER\n" << std::endl;

                      delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[tailends[superNodeID]];
                      delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[tailends[superNodeID]];
                      delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[tailends[superNodeID]];

                      std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;

                      // 2/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
                      coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]) - delta_h2_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h3_partial_pfixsum[superparent];
                              //delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                      std::cout << "\t\t" << tailends[superNodeID]  << " -h " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                                << tailends[superNodeID]  << "^2 - " << delta_h2_partial_pfixsum[superparent] << " * " << tailends[superNodeID]  << " + " << delta_h3_partial_pfixsum[superparent]  << " = " << coefficientweightList[superparent] << std::endl;

                  }

              }
              else
              {
                // std::cout << "COMPUTE " << superparent << "->" << tailends[superNodeID] << std::endl;
                std::cout << tailends[superNodeID] << " - " << superparent << "->" << tailends[superNodeID] << "\n";
                // 3/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
                coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                std::cout << "\t\t" << sortID << " -g " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                          << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;


              }

              superarcIntrinsicWeightPortal.Set(superparent,
                                                superarcIntrinsicWeightPortal.Get(superparent)+coefficientweightList[superparent]);

            } // per node in sorted order
        }

        else if (dim3)
        {// test 3D coefficients
            std::cout << "3D Coefficients Sweep" << std::endl;
            for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
            { // per node in sorted order
              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              if (prevIndex != superparent)
              {
                  prevIndex = superparent;
                  // tail-end of a branch
                  superNodeID = sortID;
              }

              std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;


              if (sortedNode == 0)
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);
              else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
                firstVertexForSuperparentPortal.Set(superparent, sortedNode);


              // CHANGES:
              // UPDATE AT REGULAR NODE: +1
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+1);

              // UPDATE AT REGULAR NODE: +area/3
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[superparent]);
              //        superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortedNode]);

              // weights before 2024-08-27:
              // superarcIntrinsicWeightPortal.Set(superparent, superarcIntrinsicWeightPortal.Get(superparent)+weightList[sortID]);

              // weights after 2024-08-27:
              delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[sortID];
              delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[sortID];
              delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[sortID];
              delta_h4_partial_pfixsum[superparent] += vx_delta_h4_sum[sortID];

              // 1/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
              coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (sortID*sortID*sortID) +\
                                                   delta_h2_partial_pfixsum[superparent] * (sortID*sortID) +\
                                                   delta_h3_partial_pfixsum[superparent] * sortID + \
                                                   delta_h4_partial_pfixsum[superparent];

              std::cout << "\t\t" << sortID << " - " << superparent << "\t"
                        << delta_h1_partial_pfixsum[superparent] << " * " << sortID << "^3 + "
                        << delta_h2_partial_pfixsum[superparent] << " * " << sortID << "^2 + "
                        << delta_h3_partial_pfixsum[superparent] << " * " << sortID << " + "
                        << delta_h4_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << std::endl;


              if(sortedNode != contourTree.Arcs.GetNumberOfValues()-1)
              {
                  vtkm::Id nextSortID = nodesPortal.Get(sortedNode+1);
                  vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

                  if(nextSuperparent != superparent)
                  {
                      // std::cout << nextSuperparent << " -vs- " << superparent <<  " -- TRIGGER\n" << std::endl;

                      delta_h1_partial_pfixsum[superparent] += vx_delta_h1_sum[tailends[superNodeID]];
                      delta_h2_partial_pfixsum[superparent] += vx_delta_h2_sum[tailends[superNodeID]];
                      delta_h3_partial_pfixsum[superparent] += vx_delta_h3_sum[tailends[superNodeID]];
                      delta_h4_partial_pfixsum[superparent] += vx_delta_h4_sum[tailends[superNodeID]];
                      std::cout << sortID << " - " << superparent << "->" << tailends[superNodeID] << "\n";// << hypernode << " " << hyperarct << std::endl;

                      // 2/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
//                      coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]) - delta_h2_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h3_partial_pfixsum[superparent];
                              //delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];
                      coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]*tailends[superNodeID]) +\
                                                           delta_h2_partial_pfixsum[superparent] * (tailends[superNodeID]*tailends[superNodeID]) +\
                                                           delta_h3_partial_pfixsum[superparent] * tailends[superNodeID] + \
                                                           delta_h4_partial_pfixsum[superparent];

//                      std::cout << "\t\t" << tailends[superNodeID]  << " -h " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
//                                << tailends[superNodeID]  << "^2 - " << delta_h2_partial_pfixsum[superparent] << " * " << tailends[superNodeID]  << " + " << delta_h3_partial_pfixsum[superparent]  << " = " << coefficientweightList[superparent] << std::endl;

                      std::cout << "\t\t" << sortID << " -h " << superparent << "\t"
                                << delta_h1_partial_pfixsum[superparent] << " * " << tailends[superNodeID] << "^3 + "
                                << delta_h2_partial_pfixsum[superparent] << " * " << tailends[superNodeID] << "^2 + "
                                << delta_h3_partial_pfixsum[superparent] << " * " << tailends[superNodeID] << " + "
                                << delta_h4_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << std::endl;


                  }

              }
              else
              {
                // std::cout << "COMPUTE " << superparent << "->" << tailends[superNodeID] << std::endl;
                std::cout << tailends[superNodeID] << " - " << superparent << "->" << tailends[superNodeID] << "\n";
                // 3/3 2024-11-17 FIX: replace '+=' with '=' as the delta_h1_partial_pfixsum already prefixes the volumes
                coefficientweightList[superparent] = delta_h1_partial_pfixsum[superparent] * tailends[superNodeID] + delta_h2_partial_pfixsum[superparent];

                std::cout << "\t\t" << sortID << " -g " << superparent << " (" << delta_h1_partial_pfixsum[superparent] << " * "
                          << tailends[superNodeID] << " + " << delta_h2_partial_pfixsum[superparent] << " = " << coefficientweightList[superparent] << "\n" << std::endl;


              }

              superarcIntrinsicWeightPortal.Set(superparent,
                                                superarcIntrinsicWeightPortal.Get(superparent)+coefficientweightList[superparent]);

            } // per node in sorted order
        }





        // -------------------------------------- SWEEP  ------------------------------------- //























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

        // set up the array which tracks which supernodes to deal with on which iteration:
        // 1) VOLUMES:
        auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
        auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
        auto supernodeTransferWeightPortal = supernodeTransferWeight.WritePortal();
        auto superarcDependentWeightPortal = superarcDependentWeight.WritePortal();
        auto hyperarcDependentWeightPortal = hyperarcDependentWeight.WritePortal();

        auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();

        // 2) COEFFICIENTS:
//        auto firstSupernodePerIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
//        auto firstHypernodePerIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();
//        auto supernodeTransferWeightCoeffPortal    = supernodeTransferWeightCoeff.WritePortal();
//        auto superarcDependentWeightCoeffPortal    = superarcDependentWeightCoeff.WritePortal();
//        auto hyperarcDependentWeightCoeffPortal    = hyperarcDependentWeightCoeff.WritePortal();

//        superarcIntrinsicWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
//        auto superarcIntrinsicWeightCoeffPortal = superarcIntrinsicWeightCoeff.WritePortal();

        superarcDependentWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcDependentWeightCoeffPortal = superarcDependentWeightCoeff.WritePortal();

        supernodeTransferWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto supernodeTransferWeightCoeffPortal = supernodeTransferWeightCoeff.WritePortal();

        hyperarcDependentWeightCoeff.Allocate(contourTree.Hyperarcs.GetNumberOfValues());
        auto hyperarcDependentWeightCoeffPortal = hyperarcDependentWeightCoeff.WritePortal();

        Coefficients SAlocalIntrinsic;
        Coefficients SNlocalTransfer;
        Coefficients SAlocalDependent;
        Coefficients HAlocalDependent;

        // initialise the coefficient arrays:
        // intrinsic to

        std::cout << "number of arcs:"      << contourTree.Arcs.GetNumberOfValues() << std::endl;
        std::cout << "number of superarcs:" << contourTree.Superarcs.GetNumberOfValues() << std::endl;
        std::cout << "number of hyperarcs:" << contourTree.Hyperarcs.GetNumberOfValues() << std::endl;


//        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
//            contourTree.Supernodes.GetNumberOfValues()

//        contourTree.

        // Initialise the super{node/arc} weight arrays with 0s
        // (data below depends on the super{nodes/arcs}
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Supernodes.GetNumberOfValues(); sortedNode++)
        { // per node in sorted order
            SAlocalIntrinsic.h1 = 0.0;
            SAlocalIntrinsic.h2 = 0.0;
            SAlocalIntrinsic.h3 = 0.0;
            SAlocalIntrinsic.h4 = 0.0;

            SNlocalTransfer.h1 = 0.0;
            SNlocalTransfer.h2 = 0.0;
            SNlocalTransfer.h3 = 0.0;
            SNlocalTransfer.h4 = 0.0;

            SAlocalDependent.h1 = 0.0;
            SAlocalDependent.h2 = 0.0;
            SAlocalDependent.h3 = 0.0;
            SAlocalDependent.h4 = 0.0;

            superarcIntrinsicWeightCoeffPortal.Set(sortedNode, SAlocalIntrinsic);
            supernodeTransferWeightCoeffPortal.Set(sortedNode, SNlocalTransfer);
            superarcDependentWeightCoeffPortal.Set(sortedNode, SAlocalDependent);
        }

        // Initialise the hyper{node/arc} weight arrays with 0s
        // (data below depends on the hyper{nodes/arcs}
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Hypernodes.GetNumberOfValues(); sortedNode++)
        { // per node in sorted order
            HAlocalDependent.h1 = 0.0;
            HAlocalDependent.h2 = 0.0;
            HAlocalDependent.h3 = 0.0;
            HAlocalDependent.h4 = 0.0;

            hyperarcDependentWeightCoeffPortal.Set(sortedNode, HAlocalDependent);
        }


        std::cout << "SuperARC check: " << std::endl;
//        for(vtkm::Id nodeID = 0; nodeID < contourTree.Hyperarcs.GetNumberOfValues(); nodeID++)
        for (vtkm::Id hypernode = 0; hypernode < contourTree.Hypernodes.GetNumberOfValues(); hypernode++)
        {
            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
//            std::cout << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

            std::cout << hypernode << ") " << hyperarcTarget << std::endl;
        }


        vtkm::Id previousParent = -1;
        vtkm::Id nextParent = -1;

        std::cout << "----------------------------------------------------------------------" << std::endl;

        std::cout << "Arcs count: " << contourTree.Arcs.GetNumberOfValues() << std::endl;

        std::cout << "Supernodes   count: " << contourTree.Supernodes.GetNumberOfValues() << std::endl;
        std::cout << "Superarcs    count: " << contourTree.Superarcs.GetNumberOfValues() << std::endl;
        std::cout << "Superparents count: " << contourTree.Superparents.GetNumberOfValues() << std::endl;

        std::cout << "Hypernodes   count: " << contourTree.Hypernodes.GetNumberOfValues() << std::endl;
        std::cout << "Hyperarcs    count: "  << contourTree.Hyperarcs.GetNumberOfValues() << std::endl;
        std::cout << "Hyperparents count: "  << contourTree.Hyperparents.GetNumberOfValues() << std::endl;

        std::cout << "----------------------------------------------------------------------" << std::endl;

        std::cout << "'Node -> Superarc' relationships" << std::endl;

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);
//            vtkm::Id hyperparent = hyperparentsPortal.Get(sortID);

            vtkm::Id nextSortID = (sortedNode+1 == contourTree.Arcs.GetNumberOfValues()) ? 0 : nodesPortal.Get(sortedNode+1);
//            vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

            std::cout << sortID << ") superparent = " << superparent << ", next = " << nextSortID << std::endl;
        }

        std::cout << "'Supernode -> Hyperarc (Hyperparent)' relationships" << std::endl;


//        for (vtkm::Id sortedSupernode = 0; sortedSupernode < contourTree.Supernodes.GetNumberOfValues(); sortedSupernode++)
//        {
//            std::cout << "sortedSupernode = " << sortedSupernode << std::endl;
//            vtkm::Id sortID = supernodesPortal.Get(sortedSupernode );
//            std::cout << "sortID = " << sortID << std::endl;
//            vtkm::Id hyperarc = hyperarcsPortal.Get(sortID);
//            vtkm::Id hyperparent = hyperparentsPortal.Get(sortID);


//            vtkm::Id nextSortID = (sortedSupernode+1 == contourTree.Superarcs.GetNumberOfValues()) ? 0 : nodesPortal.Get(sortedSupernode+1);
////            vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);

//            std::cout << sortID << ") hyperarc = " << hyperarc << ", hyperparent = " << hyperparent << ", next = " << nextSortID << std::endl;
//        }

        std::cout << "----------------------------------------------------------------------" << std::endl;

        // TODO: Assumption - for now, treat hyperarcs same as superarcs


        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);
//            vtkm::Id hyperparent = hyperparentsPortal.Get(sortID);

            vtkm::Id nextSortID = (sortedNode+1 == contourTree.Arcs.GetNumberOfValues()) ? 0 : nodesPortal.Get(sortedNode+1);
            vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);
//            vtkm::Id nextHyperparent = hyperparentsPortal.Get(nextSortID);


            if(superparent != previousParent )
            {
                std::cout << "----------------------------------------------------------------------" << std::endl;
                previousParent = superparent;
            }

            SAlocalIntrinsic.h1 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h1 + vx_delta_h1_sum[sortID];
            SAlocalIntrinsic.h2 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h2 + vx_delta_h2_sum[sortID];
            SAlocalIntrinsic.h3 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h3 + vx_delta_h3_sum[sortID];
            SAlocalIntrinsic.h4 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h4 + vx_delta_h4_sum[sortID];

            superarcIntrinsicWeightCoeffPortal.Set(superparent, SAlocalIntrinsic);

            vtkm::Id supertarget = vtkm::cont::ArrayGetValue(superparent, contourTree.Superarcs);

            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(superparent)); //sortedNode));

            std::cout << "SP: " << superparent << " - " << sortID << " (Target: " << hyperarcTarget << "), next=" << nextSuperparent << "\n";

            std::cout << "Intrinsic: " << vx_delta_h1_sum[sortID] << " " << vx_delta_h2_sum[sortID] << " " << vx_delta_h3_sum[sortID] << " " << vx_delta_h4_sum[sortID] << std::endl;
            std::cout << "\tIntrinsic Total: " << superarcIntrinsicWeightCoeffPortal.Get(superparent).h1 << " " << superarcIntrinsicWeightCoeffPortal.Get(superparent).h2 << " " << superarcIntrinsicWeightCoeffPortal.Get(superparent).h3 << " " << superarcIntrinsicWeightCoeffPortal.Get(superparent).h4 << std::endl << std::endl;

//            if (!vtkm::worklet::contourtree_augmented::NoSuchElement(supertarget))
//            {
//                SNlocalTransfer.h1 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h1 + vx_delta_h1_sum[sortID];
//                SNlocalTransfer.h2 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h2 + vx_delta_h2_sum[sortID];
//                SNlocalTransfer.h3 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h3 + vx_delta_h3_sum[sortID];
//                SNlocalTransfer.h4 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h4 + vx_delta_h4_sum[sortID];

//                supernodeTransferWeightCoeffPortal.Set(hyperarcTarget, SNlocalTransfer);

//                std::cout << "Transfer: " << vx_delta_h1_sum[sortID] << " " << vx_delta_h2_sum[sortID] << " " << vx_delta_h3_sum[sortID] << " " << vx_delta_h4_sum[sortID] << " to target: " << hyperarcTarget << std::endl;
//                std::cout << "\tTarget Total: " << SNlocalTransfer.h1 << " " << SNlocalTransfer.h2 << " " << SNlocalTransfer.h3 << " " << SNlocalTransfer.h4 << " to target: " << hyperarcTarget << std::endl;
//            }


            // TOGGLEWEIGHT Enable Dependent Weight computation

            SAlocalDependent.h1 = 0.0; //SAlocalIntrinsic.h1 + supernodeTransferWeightCoeffPortal.Get(superparent).h1;
            SAlocalDependent.h2 = 0.0; //SAlocalIntrinsic.h2 + supernodeTransferWeightCoeffPortal.Get(superparent).h2;
            SAlocalDependent.h3 = 0.0; //SAlocalIntrinsic.h3 + supernodeTransferWeightCoeffPortal.Get(superparent).h3;
            SAlocalDependent.h4 = 0.0; //SAlocalIntrinsic.h4 + supernodeTransferWeightCoeffPortal.Get(superparent).h4;

            superarcDependentWeightCoeffPortal.Set(superparent, SAlocalDependent);

            // TODO: Assumption - for now, treat hyperarcs same as superarcs
            hyperarcDependentWeightCoeffPortal.Set(superparent, SAlocalDependent);

            std::cout << "Dependent: " << supernodeTransferWeightCoeffPortal.Get(superparent).h1 << " " << supernodeTransferWeightCoeffPortal.Get(superparent).h2 << " " << supernodeTransferWeightCoeffPortal.Get(superparent).h3 << " " << supernodeTransferWeightCoeffPortal.Get(superparent).h4 << " to target: " << superparent << std::endl;
            std::cout << "\tDependent Total: " << SAlocalDependent.h1 << " " << SAlocalDependent.h2 << " " << SAlocalDependent.h3 << " " << SAlocalDependent.h4 << " of target: " << superparent << std::endl << std::endl;


//            if(superparent != previousParent )
//            {
//                std::cout << "----------------------------------------------------------------------" << std::endl;
//                previousParent = superparent;
//            }



            // TOGGLEWEIGHT Enable Transfer Weight computation

            if ((!vtkm::worklet::contourtree_augmented::NoSuchElement(supertarget)) && (superparent != nextSuperparent ))
            {
                SNlocalTransfer.h1 = 0.0; // supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h1 + SAlocalDependent.h1;
                SNlocalTransfer.h2 = 0.0; // supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h2 + SAlocalDependent.h2;
                SNlocalTransfer.h3 = 0.0; // supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h3 + SAlocalDependent.h3;
                SNlocalTransfer.h4 = 0.0; // supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h4 + SAlocalDependent.h4;

                supernodeTransferWeightCoeffPortal.Set(hyperarcTarget, SNlocalTransfer);

                std::cout << "Transfer: " << SAlocalDependent.h1 << " " << SAlocalDependent.h2 << " " << SAlocalDependent.h3 << " " << SAlocalDependent.h4 << " to target: " << hyperarcTarget << std::endl;
                std::cout << "\tTarget (" << hyperarcTarget << ") Total: " << SNlocalTransfer.h1 << " " << SNlocalTransfer.h2 << " " << SNlocalTransfer.h3 << " " << SNlocalTransfer.h4 << " to target: " << hyperarcTarget << std::endl << std::endl;
            }


//             std::cout << "SP: " << superparentsPortal.Get(sortedNode) << " - " << nodesPortal.Get(sortedNode) << "\n";
             std::cout << std::endl;
//             std::cout << "SP: " << superparent << " - " << sortID << "(" << sortedNode << ")\n";
        }

        std::cout << "------------ setting SAlocalIntrinsic ------------" << std::endl;

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Supernodes.GetNumberOfValues(); sortedNode++)
        { // per node in sorted order

            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);

//            SAlocalIntrinsic.h1 = vx_delta_h1_sum[sortedNode];
//            SAlocalIntrinsic.h2 = vx_delta_h2_sum[sortedNode];
//            SAlocalIntrinsic.h3 = vx_delta_h3_sum[sortedNode];
//            SAlocalIntrinsic.h4 = vx_delta_h4_sum[sortedNode];

//            superarcIntrinsicWeightCoeffPortal.Set(sortedNode, SAlocalIntrinsic);
            std::cout << "SP: " << superparent << " - "; //superparentsPortal.Get(sortedNode) << " - ";
            std::cout << "SA: " << sortedNode << std::setw(10) << " " << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h1 << "\t";
            std::cout << std::setw(10) << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h2 << "\t";
            std::cout << std::setw(10) << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h3 << "\t";
            std::cout << std::setw(10) << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h4 << "\t"; //std::endl;

            std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h1 << "\t";
            std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h2 << "\t";
            std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h3 << "\t";
            std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h4 << "\t";


            std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h1 << "\t";
            std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h2 << "\t";
            std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h3 << "\t";
            std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h4 << std::endl;

//            std::cout << "SA: " << sortedNode << " " << SAlocalIntrinsic.h1 << " ";
//            std::cout << SAlocalIntrinsic.h2 << " ";
//            std::cout << SAlocalIntrinsic.h3 << " ";
//            std::cout << SAlocalIntrinsic.h4 << std::endl;

        }








        // ===================== ITERATIVE WEIGHT PROPAGATION INWARDS ===================== //







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
        //        for (vtkm::Id iteration = 0; iteration < nIterations; iteration++)
        // try to run for one more iteration to capture the whole tree
        for (vtkm::Id iteration = 0; iteration < nIterations+1; iteration++)
        { // per iteration

          std::string indent = "\t";

          std::cout << "Iteration: " << iteration << "\n{\n" << std::endl;

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


          std::cout << indent << "SUPERARCS:\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << indent << indent << supernode << " ";
          }
          std::cout << std::endl << indent << "}\n" << std::endl;


          std::cout << indent << "transfer (simple):\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << indent << indent << supernodeTransferWeightPortal.Get(supernode) << " ";
          }
          std::cout << std::endl << indent << "}\n" << std::endl;

          std::cout << indent << "transfer (coeff):\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
               std::cout << indent << indent << supernodeTransferWeightCoeffPortal.Get(supernode).h1 << " " \
                         << supernodeTransferWeightCoeffPortal.Get(supernode).h2 << " " \
                         << supernodeTransferWeightCoeffPortal.Get(supernode).h3 << " " \
                         << supernodeTransferWeightCoeffPortal.Get(supernode).h4 << std::endl;
          }
          std::cout << indent << "}\n" << std::endl;


          std::cout << indent << "intrinsic (simple):\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout  << indent << indent << superarcIntrinsicWeightPortal.Get(supernode) << " ";
          }
          std::cout << std::endl << indent << "}\n" << std::endl;


          std::cout << indent << "intrinsic (coeff):\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
              std::cout << indent << indent << superarcIntrinsicWeightCoeffPortal.Get(supernode).h1 << " " \
                        << superarcIntrinsicWeightCoeffPortal.Get(supernode).h2 << " " \
                        << superarcIntrinsicWeightCoeffPortal.Get(supernode).h3 << " " \
                        << superarcIntrinsicWeightCoeffPortal.Get(supernode).h4 << std::endl;
          }
          std::cout << indent << "}\n" << std::endl;


          Coefficients step1Dependent;
          Coefficients step2Dependent;
          Coefficients step3Dependent;
          Coefficients step4HyperarcDependent;
          Coefficients step4SupernodeTransfer;

          std::map<vtkm::Id, vtkm::Id> tailends;

          std::cout << indent << "tailends:\n" << indent << "{\n";
          for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues(); supernode++)
          {
              vtkm::Id superNode = supernodesPortal.Get(supernode);

              std::cout << indent << indent << supernode << " - " << superNode << "->" << supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode))) << std::endl;

              tailends.insert(std::make_pair(superNode, supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode)))));
          }
          std::cout << indent << "}\n" << std::endl;

          std::cout << indent <<  "// step 1: Calculate DEPENDENT weight, which is INTRINSIC + TRANSFER" << std::endl;
          std::cout << indent << "step 1 (simple):\n" << indent << "{\n";
          // so, step 1: add xfer + int & store in dependent weight
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
            std::cout << indent << indent << "DependentSA[" << supernode << "] = TransferSN[" << supernode << "] + " << "IntrinsicSA[" << supernode << "]\n" << std::endl;

            superarcDependentWeightPortal.Set(supernode,
                                              supernodeTransferWeightPortal.Get(supernode) +
                                                superarcIntrinsicWeightPortal.Get(supernode));

            std::cout << indent << indent << supernode << " = " << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << "\n" << std::endl;
          }
          std::cout << std::endl << indent << indent << "(SIMPL.) - DEPENDENT = TRANSFER + INTRINSIC\n" << indent << "}\n" << std::endl;




          std::cout << indent << "step 1 (coeff):\n" << indent << "{\n";
          // so, step 1: add xfer + int & store in dependent weight
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
            std::cout << indent << indent << "DependentSA[" << supernode << "] = TransferSN[" << supernode << "] + " << "IntrinsicSA[" << supernode << "]\n" << std::endl;

            step1Dependent.h1 = supernodeTransferWeightCoeffPortal.Get(supernode).h1 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h1;
            step1Dependent.h2 = supernodeTransferWeightCoeffPortal.Get(supernode).h2 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h2;
            step1Dependent.h3 = supernodeTransferWeightCoeffPortal.Get(supernode).h3 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h3;
            step1Dependent.h4 = supernodeTransferWeightCoeffPortal.Get(supernode).h4 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h4;

            superarcDependentWeightCoeffPortal.Set(supernode, step1Dependent);

//            std::cout << supernodeTransferWeightPortal.Get(supernode) + superarcIntrinsicWeightPortal.Get(supernode) << " ";
            std::cout << indent << indent << supernode << " = " << step1Dependent.h1 << " " \
                              << step1Dependent.h2 << " " \
                              << step1Dependent.h3 << " " \
                              << step1Dependent.h4 << "\n" << std::endl;
          }
          std::cout << std::endl << indent << indent << "(COEFF.) - DEPENDENT = TRANSFER + INTRINSIC\n" << indent << "}\n" << std::endl;




          std::cout << indent << "// step 2: Calculate DEPENDENT weight, which is INTRINSIC + TRANSFER" << std::endl;

          std::cout << indent << "step 2 (simple):\n" << indent << "{\n" << indent << indent;
          std::cout << "from: " << firstSupernode+1 << " to " << lastSupernode << std::endl;
          //          std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
          // step 2: perform prefix sum on the dependent weight range
          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
          {
            std::cout << indent << indent << "DependentSA[" << supernode << "] = DependentSA[" << supernode << "] + " << "DependentSA[" << supernode - 1 << "]\n" << std::endl;

            superarcDependentWeightPortal.Set(supernode,
                                              superarcDependentWeightPortal.Get(supernode) +
                                                superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

            std::cout << indent << indent << supernode << " = " << superarcDependentWeightPortal.Get(supernode) << std::endl;

          }

          std::cout << std::endl << indent << indent << "(SIMPL.) - DEPENDENT = DEPENDENT[CURRENT] + DEPENDENT[PREVIOUS]\n" << indent << "}\n" << std::endl;

          std::cout << std::endl;
          std::cout << " " << std::endl;


          std::cout << indent << "step 2 (coeff):\n" << indent << "{\n" << indent << indent;
          //          std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
          std::cout << "from: " << firstSupernode+1 << " to " << lastSupernode << std::endl;
          // step 2: perform prefix sum on the dependent weight range
          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
          {
              std::cout << indent << indent << "DependentSA[" << supernode << "] = DependentSA[" << supernode << "] + " << "DependentSA[" << supernode - 1 << "]\n" << std::endl;

              step2Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 + superarcDependentWeightCoeffPortal.Get(supernode-1).h1;
              step2Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 + superarcDependentWeightCoeffPortal.Get(supernode-1).h2;
              step2Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 + superarcDependentWeightCoeffPortal.Get(supernode-1).h3;
              step2Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 + superarcDependentWeightCoeffPortal.Get(supernode-1).h4;

              superarcDependentWeightCoeffPortal.Set(supernode, step2Dependent);

              std::cout << indent << indent << supernode << " = " << step2Dependent.h1 << " " \
                                << step2Dependent.h2 << " " \
                                << step2Dependent.h3 << " " \
                                << step2Dependent.h4 << std::endl;
          }

          std::cout << std::endl << indent << indent << "(COEFF.) - DEPENDENT = DEPENDENT[CURRENT] + DEPENDENT[PREVIOUS]\n" << indent << "}\n" << std::endl;










          // step 3: subtract out the dependent weight of the prefix to the entire hyperarc. This will be a transfer, but for now, it's easier
          // to show it in serial. NB: Loops backwards so that computation uses the correct value
          // As a bonus, note that we test > firstsupernode, not >=.  This is because we've got unsigned integers, & otherwise it will not terminate
          // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
          std::cout << indent << "// step 3: Subtract the last SA dependent weight of the iteration from the total" << std::endl;
          std::cout << indent << "step 3 (simple) - subtract:\n" << indent << "{\n";
          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
              continue;

            std::cout << indent << indent << "DependentSA[" << supernode << "] = DependentSA[" << supernode << "] - " << "DependentSA[" << hyperparentSuperID - 1 << "]\n";
            std::cout << indent << indent << "[" << superarcDependentWeightPortal.Get(supernode) << "] = [" << superarcDependentWeightPortal.Get(supernode) << "] - " << "[" << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << "]\n" << std::endl;

            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
            superarcDependentWeightPortal.Set(
              supernode,
              superarcDependentWeightPortal.Get(supernode) -
                superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

            //std::cout << indent << indent << supernode << "(" << hyperparentSuperID << ")" << " - " << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << std::endl;
            std::cout << indent << indent << supernode << "(" << hyperparentSuperID << ")" << " - " << superarcDependentWeightPortal.Get(hyperparentSuperID - 1) << std::endl;

          } // per supernode

          std::cout << std::endl << indent << "}\n" << std::endl;


          std::cout << indent << "step 3 (coeff) - subtract:\n" << indent << "{\n";
          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
              continue;

            std::cout << indent << indent << "DependentSA[" << supernode << "] = DependentSA[" << supernode << "] + " << "DependentSA[" << supernode - 1 << "]\n" << std::endl;

            step3Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h1;
            step3Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h2;
            step3Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h3;
            step3Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h4;

            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
            superarcDependentWeightCoeffPortal.Set(supernode, step3Dependent);

            std::cout << indent << indent << supernode << " = " << step3Dependent.h1 << " " \
                              << step3Dependent.h2 << " " \
                              << step3Dependent.h3 << " " \
                              << step3Dependent.h4 << std::endl;

          } // per supernode

          std::cout << std::endl << indent << "}\n" << std::endl;


          std::cout << indent << "step 4 (simple) - target transfer weights:\n" << indent << "{\n";
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
            std::cout << indent << indent << "DependentHA[" << hypernode << "] = DependentSA[" << lastSuperarc << "]\n";
            std::cout << indent << indent << "[" << hyperarcDependentWeightPortal.Get(hypernode) << "] = [" << superarcDependentWeightPortal.Get(lastSuperarc) << "]\n";
            hyperarcDependentWeightPortal.Set(hypernode,
                                              superarcDependentWeightPortal.Get(lastSuperarc));

            // note that in parallel, this will have to be split out as a sort & partial sum in another array
            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
            std::cout << indent << indent << "TransferSN[" << hyperarcTarget << "] = TransferSN[" << hyperarcTarget << "] + DependentHA[" << hypernode << "]\n";
            std::cout << indent << indent << "[" << supernodeTransferWeightPortal.Get(hyperarcTarget) << "] = [" << supernodeTransferWeightPortal.Get(hyperarcTarget)  << "] + [" << hyperarcDependentWeightPortal.Get(hypernode) << "]\n" << std::endl;
            supernodeTransferWeightPortal.Set(hyperarcTarget,
                                              supernodeTransferWeightPortal.Get(hyperarcTarget) +
                                                hyperarcDependentWeightPortal.Get(hypernode));


            std::cout << indent << indent << "hyperarcDependentWeightPortal:" << std::endl << indent << indent;
            std::cout << indent << indent << hypernode << " - " << hyperarcDependentWeightPortal.Get(hypernode) << std::endl;

            std::cout << indent << indent << "supernodeTransferWeightPortal:" << std::endl << indent << indent;
            std::cout << indent << indent << hyperarcTarget << " - " << supernodeTransferWeightPortal.Get(hyperarcTarget) << std::endl;

          } // per hypernode

          std::cout << std::endl << indent << "}\n" << std::endl;



          std::cout << indent << "step 4 (coeff) - target transfer weights:\n" << indent << "{\n";
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

            std::cout << indent << indent << "DependentHA[" << hypernode << "] = DependentSA[" << lastSuperarc << "]\n";
            step4HyperarcDependent.h1 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h1;
            step4HyperarcDependent.h2 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h2;
            step4HyperarcDependent.h3 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h3;
            step4HyperarcDependent.h4 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h4;

            // now, given the last superarc for the hyperarc, transfer the dependent weight
            hyperarcDependentWeightCoeffPortal.Set(hypernode, step4HyperarcDependent);

            // note that in parallel, this will have to be split out as a sort & partial sum in another array
            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));
            std::cout << indent << indent << "TransferSN[" << hyperarcTarget << "] = TransferSN[" << hyperarcTarget << "] + DependentHA[" << hypernode << "]\n" << std::endl;
            step4SupernodeTransfer.h1 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h1 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h1;
            step4SupernodeTransfer.h2 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h2 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h2;
            step4SupernodeTransfer.h3 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h3 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h3;
            step4SupernodeTransfer.h4 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h4 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h4;

            // now, given the last superarc for the hyperarc, transfer the dependent weight
//            hyperarcDependentWeightCoeffPortal.Set(hypernode, step4HyperarcDependent);

            supernodeTransferWeightCoeffPortal.Set(hyperarcTarget, step4SupernodeTransfer);

            std::cout << indent << indent << "hyperarcDependentWeightCoeffPortal:" << std::endl << indent << indent;
            std::cout << hypernode << " = " << step4HyperarcDependent.h1 << " " \
                              << step4HyperarcDependent.h2 << " " \
                              << step4HyperarcDependent.h3 << " " \
                              << step4HyperarcDependent.h4 << std::endl;

            std::cout << indent << indent << "supernodeTransferWeightCoeffPortal:" << std::endl << indent << indent;
            std::cout << hyperarcTarget << " = " << step4SupernodeTransfer.h1 << " " \
                              << step4SupernodeTransfer.h2 << " " \
                              << step4SupernodeTransfer.h3 << " " \
                              << step4SupernodeTransfer.h4 << std::endl << std::endl;

          } // per hypernode

          std::cout << std::endl << indent << "}\n" << std::endl;


          std::cout << std::endl;
          std::cout << indent << "final (simple):\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
    //        superarcDependentWeightPortal.Set(supernode,
    //                                          superarcDependentWeightPortal.Get(supernode) +
    //                                            superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";

            std::cout << indent << indent << supernode << " = " << superarcDependentWeightPortal.Get(supernode) << std::endl;

          }
          std::cout << std::endl << indent << "}\n" << std::endl;

//          for (vtkm::Id nodeID = firstSupernode; nodeID < lastSupernode; nodeID++)

          //JUMP
          std::cout << indent << "sorted nodes:\n" << indent << "{\n";
          for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues(); supernode++)
          {
              vtkm::Id superNode = supernodesPortal.Get(supernode);

              std::cout << indent << indent << supernode << " - " << superNode << "->" << tailends[superNode] << std::endl;
          }

          std::cout << std::endl << indent << "}\n" << std::endl;

          std::cout << indent << "final (coeff):\n" << indent << "{\n";
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {

            vtkm::Id superNode = supernodesPortal.Get(supernode);

            std::cout << indent << indent << supernode << "(" << superNode << "->" << tailends[superNode] << ") = " << superarcDependentWeightCoeffPortal.Get(supernode).h1 << " " \
                              << superarcDependentWeightCoeffPortal.Get(supernode).h2 << " " \
                              << superarcDependentWeightCoeffPortal.Get(supernode).h3 << " " \
                              << superarcDependentWeightCoeffPortal.Get(supernode).h4 << " = "; // std::endl;

            double a_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h1 * std::pow(tailends[superNode], 3);
            double b_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h2 * std::pow(tailends[superNode], 2);
            double c_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h3 * tailends[superNode];
            double d_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h4;

            std::cout << a_coeff + b_coeff + c_coeff + d_coeff << std::endl;
          }

          std::cout << std::endl << indent << "}\n" << std::endl;

          std::cout << std::endl << "}\n" << std::endl;



          for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Supernodes.GetNumberOfValues(); sortedNode++)
          { // per node in sorted

              vtkm::Id sortID = nodesPortal.Get(sortedNode);
              vtkm::Id superparent = superparentsPortal.Get(sortID);

              std::cout << "SP: " << superparent << " - "; //superparentsPortal.Get(sortedNode) << " - ";
              std::cout << "SA: " << sortedNode << std::setw(10) << " " << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h1 << "\t";
              std::cout << std::setw(10) << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h2 << "\t";
              std::cout << std::setw(10) << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h3 << "\t";
              std::cout << std::setw(10) << superarcIntrinsicWeightCoeffPortal.Get(sortedNode).h4 << "\t"; //std::endl;

              std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h1 << "\t";
              std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h2 << "\t";
              std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h3 << "\t";
              std::cout << std::setw(10) << superarcDependentWeightCoeffPortal.Get(sortedNode).h4 << "\t";

              std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h1 << "\t";
              std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h2 << "\t";
              std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h3 << "\t";
              std::cout << std::setw(10) << supernodeTransferWeightCoeffPortal.Get(sortedNode).h4 << std::endl;

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

        std::cout << "END ComputeVolumeWeightsSerialStructCoefficients" << std::endl;

    }  // END ComputeVolumeWeightsSerialStructCoefficients



























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

      // COMMS: Both 'intrinsic' and 'dependent' weights come precomputed from 'ComputeVolumeWeights()'
      // ... the following just sets up the read portals for both
      auto superarcDependentWeightPortal = superarcDependentWeight.ReadPortal();
      auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.ReadPortal();

      // CHANGE: because at the moment the dependent and intrinsic weights have some errors ...
      // ... I write the arrays with the correct values for each to see how the BT should work.
      FloatArrayType superarcDependentWeightCorrect; //= superarcDependentWeight.ReadPortal();
      FloatArrayType superarcIntrinsicWeightCorrect; //= superarcIntrinsicWeight.ReadPortal();

      superarcDependentWeightCorrect.Allocate(superarcDependentWeightPortal.GetNumberOfValues());
      superarcIntrinsicWeightCorrect.Allocate(superarcIntrinsicWeightPortal.GetNumberOfValues());


      std::string indent = "\t";
      std::array<double, 6> realIntrinsic = {0.0208333, 0.14127, 0.178175, 0.0236112,  0.636111,                   0.0};
      std::array<double, 6> realDependent = {0.0208333, 0.14127, 0.178175, 0.0236112,  0.636111+0.0208333+0.14127, 1.0};

      auto superarcDependentWeightCorrectReadPortal = superarcDependentWeightCorrect.ReadPortal();
      auto superarcIntrinsicWeightCorrectReadPortal = superarcIntrinsicWeightCorrect.ReadPortal();

      auto superarcDependentWeightCorrectWritePortal = superarcDependentWeightCorrect.WritePortal();
      auto superarcIntrinsicWeightCorrectWritePortal = superarcIntrinsicWeightCorrect.WritePortal();


      std::cout << std::endl << "Superarc Intrinsic Weight Portal:" << std::endl;
      for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << superarcIntrinsicWeightPortal.Get(i) << std::endl;

          superarcIntrinsicWeightCorrectWritePortal.Set(i, realIntrinsic[i]);

          std::cout << indent << i << " -> " << superarcIntrinsicWeightCorrectReadPortal.Get(i) << std::endl;


      }
      std::cout << std::endl;

      std::cout << std::endl << "superarc Dependent Weight Portal:" << std::endl;
      for(int i = 0; i < superarcDependentWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << superarcDependentWeightPortal.Get(i) << std::endl;

          superarcDependentWeightCorrectWritePortal.Set(i, realDependent[i]);

          std::cout << indent << i << " -> " << superarcDependentWeightCorrectReadPortal.Get(i) << std::endl;


      }
      std::cout << std::endl;

      std::cout << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "============================ START ===============================" << std::endl;
      std::cout << "===================== BRANCH DECOMPOSITION =======================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << std::endl;

      std::cout << "ANALYSIS ComputeVolumeBranchDecompositionSerialFloat ANALYSIS" << std::endl;

      // cache the number of non-root supernodes & superarcs
      vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
      vtkm::Id nSuperarcs = nSupernodes - 1;

      std::cout << "nSupernodes = " << nSupernodes << std::endl;
      std::cout << "nSuperarcs  = " << nSuperarcs << std::endl;

      // STAGE I:  Find the upward and downwards weight for each superarc, and set up arrays

      // Allocation/initialization
      // COMMS: just allocate up-weight arrays, do not set values yet
      IdArrayType upWeight;
      upWeight.Allocate(nSuperarcs);
      auto upWeightPortal = upWeight.WritePortal();

          // FLOAT versions
          FloatArrayType upWeightFloat;
          upWeightFloat.Allocate(nSuperarcs);
          auto upWeightFloatPortal = upWeightFloat.WritePortal();

          // FLOAT CORRECT versions
          FloatArrayType upWeightFloatCorrect;
          upWeightFloatCorrect.Allocate(nSuperarcs);
          auto upWeightFloatCorrectPortal = upWeightFloatCorrect.WritePortal();

      // COMMS: just allocate down-weight arrays, do not set values yet
      IdArrayType downWeight;
      downWeight.Allocate(nSuperarcs);
      auto downWeightPortal = downWeight.WritePortal();

          // FLOAT versions
          FloatArrayType downWeightFloat;
          downWeightFloat.Allocate(nSuperarcs);
          auto downWeightFloatPortal = downWeightFloat.WritePortal();

          // FLOAT versions
          FloatArrayType downWeightFloatCorrect;
          downWeightFloatCorrect.Allocate(nSuperarcs);
          auto downWeightFloatCorrectPortal = downWeightFloatCorrect.WritePortal();

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
      double   totalVolumeFloat = 1.0; // HACK - for the small example we know total volume is 1.0
      // ... later we will get the volume from the last superarc! (dependent wght of the root holds the total volume)

      std::cout << "totalVolume = " << totalVolume << std::endl;

      std::cout << "\n" << "----------------- Computing Up/Down Weights -----------------" << std::endl;

    #ifdef DEBUG_PRINT
      std::cout << "Total Volume: " << totalVolume << std::endl;
    #endif
      // superarcs array stores the destination (supernode ID) of each superarc
      // the origin (source) of the superarc is always the supernode ID as the superarc ID
      auto superarcsPortal = contourTree.Superarcs.ReadPortal();

      // CHANGE: added for getting 'real' supernode values:
      auto supernodesPortal = contourTree.Supernodes.ReadPortal();

      // NB: Last element in array is guaranteed to be root superarc to infinity,
      // WARNING WARNING WARNING: This changes in the distributed version!
      // so we can easily skip it by not indexing to the full size
      // i.e we loop to N superarcs, not to N supernodes since N superarcs = supernodes-1
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++) //  for (vtkm::Id supernode = 0; supernode < nsupernodes-1; supernode++) vtkm::Id supernode = superarc (temporary variable)
      { // per superarc

        std::cout << "Processing superarc: " << superarc << std::endl << "{" << std::endl;

        if (IsAscending(superarcsPortal.Get(superarc))) // flag on the ID of superarc
        { // ascending superarc

          std::cout << indent << "ASCENDING\n";

            // put the lower-end first
          superarcListWritePortal.Set(superarc, // each superarc starts at the supernode at the same ID and goes to the supernode whose ID is stored in the superarc's array
                                      // pair the origin and the destination of that superarc. We store them in an edge pair and write it to the array
                                      EdgePair(superarc, MaskedIndex(superarcsPortal.Get(superarc))));

          vtkm::Id superNode = supernodesPortal.Get(superarc);

          std::cout << indent << superarc << " = " << superNode << " -> " << MaskedIndex(superarcsPortal.Get(superarc)) << std::endl << std::endl;

          // because this is an ascending superarc, the dependent weight refers to the weight at the upper end
          // so, we set the up weight based on the dependent weight
          upWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          upWeightFloatPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          upWeightFloatCorrectPortal.Set(superarc, superarcDependentWeightCorrectReadPortal.Get(superarc));


          // at the inner end, dependent weight is the total in the subtree.
          // Then there are vertices along the edge itself (intrinsic weight), ...
          // ... including the supernode at the outer end
          // So, to get the "dependent" weight in the other direction, ...
          // ... we start with totalVolume - dependent, then subtract (intrinsic - 1)
          // set the weight at the down end by using the invert operator:
          downWeightPortal.Set(superarc,
                               // below is the invert operator for node count!
                               (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
                                 (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          downWeightFloatPortal.Set(superarc,
                                   // below is the invert operator for node count!
                                   (totalVolumeFloat - superarcDependentWeightPortal.Get(superarc)) +
                                     (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          downWeightFloatCorrectPortal.Set(superarc,
                                   // below is the invert operator for node count!
                                   (totalVolumeFloat - superarcDependentWeightCorrectReadPortal.Get(superarc)) +
                                     (superarcIntrinsicWeightCorrectReadPortal.Get(superarc) - 1));

          std::cout << indent << "upWeightPortal             = " << upWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatPortal        = " << upWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatCorrectPortal = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
          std::cout << std::endl;
          std::cout << indent << "downWeightPortal             = " << downWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatPortal        = " << downWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatCorrectPortal = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;

          std::cout << "\n}\n";

        } // ascending superarc
        else
        { // descending superarc

          std::cout << indent << "DESCENDING\n";

          // lower-end is also first, but in the reverse order compared to IsAscending
          superarcListWritePortal.Set(superarc,
                                      EdgePair(MaskedIndex(superarcsPortal.Get(superarc)), superarc));

          vtkm::Id superNode = supernodesPortal.Get(superarc);

          std::cout << indent << superarc << " = " << superNode << " -> " << MaskedIndex(superarcsPortal.Get(superarc)) << std::endl << std::endl;

          downWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          downWeightFloatPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          downWeightFloatCorrectPortal.Set(superarc, superarcDependentWeightCorrectReadPortal.Get(superarc));

          // at the inner end, dependent weight is the total in the subtree.
          // Then there are vertices along the edge itself (intrinsic weight), ...
          // ... including the supernode at the outer end
          // So, to get the "dependent" weight in the other direction, ...
          // ... we start with totalVolume - dependent, then subtract (intrinsic - 1)
          upWeightPortal.Set(superarc,
                             (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          upWeightFloatPortal.Set(superarc,
                             (totalVolumeFloat - superarcDependentWeightPortal.Get(superarc)) +
                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          upWeightFloatCorrectPortal.Set(superarc,
                                   // below is the invert operator for node count!
                                   (totalVolumeFloat - superarcDependentWeightCorrectReadPortal.Get(superarc)) +
                                     (superarcIntrinsicWeightCorrectReadPortal.Get(superarc) - 1));

          std::cout << indent << "upWeightPortal      = " << upWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatPortal = " << upWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatCorrectPortal = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
          std::cout << std::endl;
          std::cout << indent << "downWeightPortal      = " << downWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatPortal = " << downWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatCorrectPortal = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;

          std::cout << "\n}\n";


        } // descending superarc
      }   // per superarc

    std::cout << std::endl;

    //JUMP

    std::cout << "Up Weights" << std::endl;
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    {
        vtkm::Id superNode = supernodesPortal.Get(superarc);

        std::cout << superarc << "(" << superNode  << ") = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
    }
    std::cout << "Down Weights" << std::endl;
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    {
        vtkm::Id superNode = supernodesPortal.Get(superarc);

        std::cout << superarc << "(" << superNode  << ") = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;
    }

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

      // OLD: Vertex count sort
//      vtkm::cont::Algorithm::Sort(
//        superarcSorter,
//                  // false / true = either ascending/descending
//                  // sort by up/down weight so that we have a segmented array
//        process_contourtree_inc_ns::SuperArcVolumetricComparator(upWeight, superarcList, false));

      // CHANGE: sort by the float up weights
      vtkm::cont::Algorithm::Sort(
        superarcSorter,
                  // false / true = either ascending/descending
                  // sort by up/down weight so that we have a segmented array
        process_contourtree_inc_ns::SuperArcVolumetricComparator(upWeightFloatCorrect, superarcList, false));



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

      // OLD: Vertex count sort
//      vtkm::cont::Algorithm::Sort(
//        superarcSorter,
//        process_contourtree_inc_ns::SuperArcVolumetricComparator(downWeightFloatCorrect, superarcList, true));


      // CHANGE: sort by the float down weights
      vtkm::cont::Algorithm::Sort(
        superarcSorter,
        process_contourtree_inc_ns::SuperArcVolumetricComparator(downWeightFloatCorrect, superarcList, true));




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

      std::cout << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "============================= START ==============================" << std::endl;
      std::cout << "========================== BRANCH DATA ===========================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << std::endl;


      //JUMP

      ProcessContourTree::ComputeBranchData(contourTree,
                                            whichBranch,    // (output)
                                            branchMinimum,  // (output)
                                            branchMaximum,  // (output)
                                            branchSaddle,   // (output)
                                            branchParent,   // (output)
                              /* (input) */ bestUpward,
                              /* (input) */ bestDownward);


      std::cout << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "============================ FINISH ==============================" << std::endl;
      std::cout << "========================== BRANCH DATA ===========================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << std::endl;

      std::cout << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "============================ FINISH ==============================" << std::endl;
      std::cout << "===================== BRANCH DECOMPOSITION =======================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << "==================================================================" << std::endl;
      std::cout << std::endl;

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
                                IdArrayType& whichBranch,       // (output)
                                IdArrayType& branchMinimum,     // (output)
                                IdArrayType& branchMaximum,     // (output)
                                IdArrayType& branchSaddle,      // (output)
                                IdArrayType& branchParent,      // (output)
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

    //JUMP

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
