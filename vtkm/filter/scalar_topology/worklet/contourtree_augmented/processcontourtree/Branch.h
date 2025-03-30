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

#ifndef vtk_m_worklet_contourtree_augmented_process_contourtree_inc_branch_h
#define vtk_m_worklet_contourtree_augmented_process_contourtree_inc_branch_h

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ContourTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/PiecewiseLinearFunction.h>

#include <cmath>
#include <algorithm>

#define DEBUG_PRINT_PACTBD 0
#define SLEEP_ON 0

namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{
namespace process_contourtree_inc
{
//    using ValueType = vtkm::Float32;
using ValueType = vtkm::Float64; //vtkm::FloatDefault;
using FloatArrayType = vtkm::cont::ArrayHandle<ValueType>;



// TODO The pointered list structure and use of std::vector don't seem to fit well with using Branch with VTKM
template <typename T>
class Branch
{
public:
  vtkm::Id OriginalId;              // Index of the extremum in the mesh
  vtkm::Id Extremum;                // Index of the extremum in the mesh
  T ExtremumVal;                    // Value at the extremum:w
  vtkm::Id Saddle;                  // Index of the saddle in the mesh (or minimum for root branch)
  T SaddleVal;                      // Corresponding value
  vtkm::Id Volume;                  // Volume
  ValueType VolumeFloat;
  Branch<T>* Parent;                // Pointer to parent, or nullptr if no parent
  std::vector<Branch<T>*> Children; // List of pointers to children

  // Create branch decomposition from contour tree
  template <typename StorageType>
  static Branch<T>* ComputeBranchDecomposition(
    const IdArrayType& contourTreeSuperparents,
    const IdArrayType& contourTreeSupernodes,
    const IdArrayType& whichBranch,
    const IdArrayType& branchMinimum,
    const IdArrayType& branchMaximum,
    const IdArrayType& branchSaddle,
    const IdArrayType& branchParent,
    const IdArrayType& sortOrder,
    const vtkm::cont::ArrayHandle<T, StorageType>& dataField,
    bool dataFieldIsSorted);


  template <typename StorageType>
  static Branch<T>* ComputeBranchDecomposition(
    const IdArrayType& contourTreeSuperparents,
    const IdArrayType& contourTreeSupernodes,
          const IdArrayType& contourTreeSuperarcs, // NEW: passed superarcs for surface component IDs vis
    const IdArrayType& whichBranch,
    const IdArrayType& branchMinimum,
    const IdArrayType& branchMaximum,
    const IdArrayType& branchSaddle,
    const IdArrayType& branchParent,
    const IdArrayType& sortOrder,
    const vtkm::cont::ArrayHandle<T, StorageType>& dataField,
    bool dataFieldIsSorted,
    const FloatArrayType& superarcDependentWeight,            // NEW: passed intrincid
    const FloatArrayType& superarcIntrinsicWeight);


  // Simplify branch composition down to target size (i.e., consisting of targetSize branches)
  void SimplifyToSize(vtkm::Id targetSize, bool usePersistenceSorter = true);

  // Print the branch decomposition
  void PrintBranchDecomposition(std::ostream& os, std::string::size_type indent = 0) const;
  //  void PrintDotBranchDecomposition(std::ostream& os, std::vector<vtkm::Id>& saddles, std::string::size_type indent = 0); // = std::vector<vtkm::Id>());
  void PrintDotBranchDecomposition(std::ostream& os, std::vector<vtkm::Id>& saddles,
                                   std::vector<vtkm::Id>& local_branches,
                                   std::vector<vtkm::Id>& depth,
                                   std::vector<vtkm::FloatDefault>& branch_weights,
                                   std::vector<vtkm::FloatDefault>& branch_weights_write,
                                   std::vector<bool>& main_branch_flags,
                                   std::vector<vtkm::Id>& depths_write,
                                   vtkm::Id parent_saddle,
                                   vtkm::Id parent_extremum,
                                   int iteration = 0,
                                   std::string::size_type indent = 0); // = std::vector<vtkm::Id>());

  // Persistence of branch
  T Persistence() { return std::fabs(ExtremumVal - SaddleVal); }

  // Destroy branch (deleting children and propagating Volume to parent)
  ~Branch();

  // Compute list of relevant/interesting isovalues
  void GetRelevantValues(int type, T eps, std::vector<T>& values) const;

  void AccumulateIntervals(int type, T eps, PiecewiseLinearFunction<T>& plf) const;

  void PrintBranchInformation(); // TODO

//  static void PrintBranchInformation(Branch<T>* root);// TODO
  static void PrintBranchInformation(Branch<T>* branch,
                                     std::vector<std::vector<vtkm::Id>>& branch_SP_map,
                                     vtkm::Id bid)
  {
//      std::cout << bid << ":" << std::endl;
      if(branch == nullptr)
      {
          std::cout << " non ";
          return;
      }

//      std::cout << branch->OriginalId << " -> ";

//      if(branch->Children.empty())
//      {
//          std::cout << " non ";
//      }
      for(int j = 0; j < branch_SP_map[bid].size(); j+=3) //j++)
      {
          if(branch_SP_map[bid][j] != bid)
          {
              std::cout << branch_SP_map[bid][j] << " -> ";
          }
      }

      if(!branch->Children.empty())
      {
          std::cout << std::endl << "Children of i=" << bid << std::endl;
      }

      for(auto child : branch->Children)
      {
          PrintBranchInformation(child, branch_SP_map, child->OriginalId);  // recursively print children
      }
  }


private:
  // Private default constructore to ensure that branch decomposition can only be created from a contour tree or loaded from storate (via static methods)
  Branch()
    : Extremum((vtkm::Id)vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT)
    , ExtremumVal(0)
    , Saddle((vtkm::Id)vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT)
    , SaddleVal(0)
    , Volume(0)
    , VolumeFloat(0.f) //0.023232323f) // VolumeFloat Placeholder - Initialisation
    , Parent(nullptr)
    , Children()
  {
  }

  // Remove symbolic perturbation, i.e., branches with zero persistence
  void removeSymbolicPerturbation();
}; // class Branch


template <typename T>
struct PersistenceSorter
{ // PersistenceSorter()
  inline bool operator()(Branch<T>* a, Branch<T>* b) { return a->Persistence() < b->Persistence(); }
}; // PersistenceSorter()


template <typename T>
struct VolumeSorter
{ // VolumeSorter()
  inline bool operator()(Branch<T>* a, Branch<T>* b)
  {
#if DEBUG_PRINT_PACTBD
      std::cout << "< VolumeSorter Branch Comparitor >: ";
      std::cout << "(" << a->ExtremumVal << ")->(" << a->SaddleVal << ") vs (" << a->ExtremumVal << ")->(" << a->SaddleVal << ") ";
      std::cout << a->Volume << " < " << b->Volume << std::endl;
#endif
      return a->Volume < b->Volume;
  }
}; // VolumeSorter()

template <typename T>
struct VolumeSorterFloat
{ // VolumeSorter()
  inline bool operator()(Branch<T>* a, Branch<T>* b)
  {
      // FIXME: No simulation of simplicity - WHY NOT?
      // (no fallback for them being identical)
#if DEBUG_PRINT_PACTBD
      std::cout << "< VolumeSorterFloat Branch Comparitor > | ";
      std::cout << "(" << a->ExtremumVal << ")->(" << a->SaddleVal << ") vs (" << b->ExtremumVal << ")->(" << b->SaddleVal << ") |";
      std::cout << a->VolumeFloat << " < " << b->VolumeFloat << " |" << std::endl << std::endl;
#endif
      return a->VolumeFloat < b->VolumeFloat;
  }
}; // VolumeSorter()



template <typename T>
template <typename StorageType>
Branch<T>* Branch<T>::ComputeBranchDecomposition(
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
{ // C)omputeBranchDecomposition()

  std::cout << "-Branch.h-> Calling: ComputeBranchDecomposition()" << std::endl;

  std::cout << "##################################################################" << std::endl;
  std::cout << "############################# START ##############################" << std::endl;
  std::cout << "##################### Branch.H Decomposition #####################" << std::endl;
  std::cout << "##################################################################" << std::endl;

  auto branchMinimumPortal = branchMinimum.ReadPortal();
  auto branchMaximumPortal = branchMaximum.ReadPortal();
  auto branchSaddlePortal = branchSaddle.ReadPortal();
  auto branchParentPortal = branchParent.ReadPortal();
  auto sortOrderPortal = sortOrder.ReadPortal();
  auto supernodesPortal = contourTreeSupernodes.ReadPortal();
  auto dataFieldPortal = dataField.ReadPortal();
  vtkm::Id nBranches = branchSaddle.GetNumberOfValues();
  std::vector<Branch<T>*> branches;
  Branch<T>* root = nullptr;
  branches.reserve(static_cast<std::size_t>(nBranches));

  for (int branchID = 0; branchID < nBranches; ++branchID)
    branches.push_back(new Branch<T>);

  // Reconstruct explicit branch decomposition from array representation
  for (std::size_t branchID = 0; branchID < static_cast<std::size_t>(nBranches); ++branchID)
  {
    branches[branchID]->OriginalId = static_cast<vtkm::Id>(branchID);
    if (!NoSuchElement(branchSaddlePortal.Get(static_cast<vtkm::Id>(branchID))))
    {
      branches[branchID]->Saddle = MaskedIndex(
        supernodesPortal.Get(MaskedIndex(branchSaddlePortal.Get(static_cast<vtkm::Id>(branchID)))));
      vtkm::Id branchMin = MaskedIndex(supernodesPortal.Get(
        MaskedIndex(branchMinimumPortal.Get(static_cast<vtkm::Id>(branchID)))));
      vtkm::Id branchMax = MaskedIndex(supernodesPortal.Get(
        MaskedIndex(branchMaximumPortal.Get(static_cast<vtkm::Id>(branchID)))));
      if (branchMin < branches[branchID]->Saddle)
        branches[branchID]->Extremum = branchMin;
      else if (branchMax > branches[branchID]->Saddle)
        branches[branchID]->Extremum = branchMax;
      else
      {
        std::cerr << "Internal error";
        return 0;
      }
    }
    else
    {
      branches[branchID]->Saddle =
        supernodesPortal.Get(MaskedIndex(branchMinimumPortal.Get(static_cast<vtkm::Id>(branchID))));
      branches[branchID]->Extremum =
        supernodesPortal.Get(MaskedIndex(branchMaximumPortal.Get(static_cast<vtkm::Id>(branchID))));
    }

    if (dataFieldIsSorted)
    {
      branches[branchID]->SaddleVal = dataFieldPortal.Get(branches[branchID]->Saddle);
      branches[branchID]->ExtremumVal = dataFieldPortal.Get(branches[branchID]->Extremum);
    }
    else
    {
      branches[branchID]->SaddleVal =
        dataFieldPortal.Get(sortOrderPortal.Get(branches[branchID]->Saddle));
      branches[branchID]->ExtremumVal =
        dataFieldPortal.Get(sortOrderPortal.Get(branches[branchID]->Extremum));
    }

    branches[branchID]->Saddle = sortOrderPortal.Get(branches[branchID]->Saddle);
    branches[branchID]->Extremum = sortOrderPortal.Get(branches[branchID]->Extremum);

    if (NoSuchElement(branchParentPortal.Get(static_cast<vtkm::Id>(branchID))))
    {
      root = branches[branchID]; // No parent -> this is the root branch
    }
    else
    {
      branches[branchID]->Parent = branches[static_cast<size_t>(
        MaskedIndex(branchParentPortal.Get(static_cast<vtkm::Id>(branchID))))];
      branches[branchID]->Parent->Children.push_back(branches[branchID]);
    }
  }

  // FIXME: This is a somewhat hackish way to compute the Volume, but it works
  // It would probably be better to compute this from the already computed Volume information
  auto whichBranchPortal = whichBranch.ReadPortal();
  auto superparentsPortal = contourTreeSuperparents.ReadPortal();
  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
  {
    branches[static_cast<size_t>(
               MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]
      ->Volume++; // Increment Volume
  }

  // 2025-01-13
#if DEBUG_PRINT_PACTBD
  std::cout << "------------- vvv Branch INTEGER weights vvv -------------" << std::endl;


  // loop through all the regular nodes, ...
  // ... then counting how many regular nodes are on each branch
  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
  {
    size_t branchID = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))));
    //      vtkm::Id sortID = supernodesPortal.Get(i);
    //    std::cout << i << " [sortID=" << sortID << "]" << " (" << MaskedIndex(superparentsPortal.Get(i)) << ") "


    std::cout << i << "[" << branchID << "] " << " (" << MaskedIndex(superparentsPortal.Get(i)) << ") "
              << branches[branchID]->Extremum << "->" << branches[branchID]->Saddle << " = "
              << branches[branchID]->Volume
              << std::endl;

  }

  std::cout << "------------- ^^^ Branch INTEGER weights ^^^ -------------" << std::endl;
#endif

#if DEBUG_PRINT_PACTBD
  std::cout << "------------- vvv Branch FLOATING weights vvv -------------" << std::endl;
  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
  {
    size_t branchID = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))));

    std::cout << MaskedIndex(superparentsPortal.Get(i)) << ") "
              << branches[branchID]->Extremum << "->" << branches[branchID]->Saddle << " = "
              << branches[branchID]->VolumeFloat
              << std::endl;
  }
  std::cout << "------------- ^^^ Branch FLOATING weights ^^^ -------------" << std::endl;
#endif

  // 2024-08-16 Get the existing volume information instead of just counting nodes

  std::cout << "Number of SuperParents: " << contourTreeSuperparents.GetNumberOfValues() << std::endl;

  // 2025-01-05 getting the right (floating point) volume information
  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
  {
    branches[static_cast<size_t>(
               MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]
      ->VolumeFloat = 0.0f; // Increment Volume

//    vtkm::Id sortID = supernodesPortal.Get(i);

//    // retrieve ID of target supernode
////    vtkm::Id superTo = superarcsPortal.Get(supernode);

//    // if this is true, it is the last pruned vertex & is omitted
//    if (NoSuchElement(superTo))
//      continue;

//    // otherwise, strip out the flags
//    superTo = MaskedIndex(superTo);

    // otherwise, we need to convert the IDs to regular mesh IDs
//    vtkm::Id regularID = sortOrderPortal.Get(MaskedIndex(sortID));

    std::cout << "'Incremented' Volume ... of [Branch " << MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i))))
              << "] of SP: (" << superparentsPortal.Get(i) << ") = " << std::endl; // << regularID << std::endl;

  }



  if (root)
  {
    root->removeSymbolicPerturbation();
  }

  std::cout << "##################################################################" << std::endl;
  std::cout << "############################# FINISH #############################" << std::endl;
  std::cout << "##################### Branch.H Decomposition #####################" << std::endl;
  std::cout << "##################################################################" << std::endl;


  return root;
} // ComputeBranchDecomposition()















































//template <typename T>
//template <typename StorageType>
//Branch<T>* Branch<T>::ComputeBranchDecomposition(
//  const IdArrayType& contourTreeSuperparents,
//  const IdArrayType& contourTreeSupernodes,
//  const IdArrayType& whichBranch,
//  const IdArrayType& branchMinimum,
//  const IdArrayType& branchMaximum,
//  const IdArrayType& branchSaddle,
//  const IdArrayType& branchParent,
//  const IdArrayType& sortOrder,
//  const vtkm::cont::ArrayHandle<T, StorageType>& dataField,
//  bool dataFieldIsSorted,
//  const FloatArrayType& superarcIntrinsicWeight,            // NEW: passed intrincid
//  const FloatArrayType& superarcDependentWeight)


//static recursiveBranching(bran



template <typename T>
template <typename StorageType>
Branch<T>* Branch<T>::ComputeBranchDecomposition(
  const IdArrayType& contourTreeSuperparents,
  const IdArrayType& contourTreeSupernodes,
  const IdArrayType& contourTreeSuperarcs,
  const IdArrayType& whichBranch,
  const IdArrayType& branchMinimum,
  const IdArrayType& branchMaximum,
  const IdArrayType& branchSaddle,
  const IdArrayType& branchParent,
  const IdArrayType& sortOrder,
  const vtkm::cont::ArrayHandle<T, StorageType>& dataField,
  bool dataFieldIsSorted,
  const FloatArrayType& superarcIntrinsicWeight,            // NEW: passed intrincid
  const FloatArrayType& superarcDependentWeight)
{ // C)omputeBranchDecomposition()

  std::cout << "ContourTreeApp->ProcessContourTree->(Branch.h->ComputeBranchDecomposition())" << std::endl;

  std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;
  std::cout << "##################################################################" << std::endl;
  std::cout << "############################# START ##############################" << std::endl;
  std::cout << "############# !MODIFIED! Branch.H Decomposition ##################" << std::endl;
  std::cout << "##################################################################" << std::endl;

  auto branchMinimumPortal = branchMinimum.ReadPortal();
  auto branchMaximumPortal = branchMaximum.ReadPortal();
  auto branchSaddlePortal = branchSaddle.ReadPortal();
  auto branchParentPortal = branchParent.ReadPortal();
  auto sortOrderPortal = sortOrder.ReadPortal();
  auto supernodesPortal = contourTreeSupernodes.ReadPortal();
  auto dataFieldPortal = dataField.ReadPortal();

  // NEW: add the read portals for branch intrinsic weights:
  auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.ReadPortal();
  // NEW: and the superarcs
  auto superarcsPortal = contourTreeSuperarcs.ReadPortal();







#if DEBUG_PRINT_PACTBD
  std::cout << std::endl << "(Branch.h->ComputeBranchDecomposition) Superarc Intrinsic Weight Portal (PASSED IN):" << std::endl;
  for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
  {
      std::cout << i << " -> " << superarcIntrinsicWeightPortal.Get(i) << std::endl;
  }
  std::cout << std::endl;
#endif

//  NOT USING DEPENDENT WEIGHTS YET
//  std::cout << std::endl << "(Branch.h->ComputeBranchDecomposition) Superarc Dependent Weight Portal:" << std::endl;
//  for(int i = 0; i < superarcDependentWeightNEWPortal.GetNumberOfValues(); i++)
//  {
//      std::cout << i << " -> " << superarcDependentWeightNEWPortal.Get(i) << std::endl;

//      superarcDependentWeightCorrectWritePortal.Set(i, realDependent[i]);

//      std::cout << indent << i << " -> " << superarcDependentWeightCorrectReadPortal.Get(i) << std::endl;
//  }







  vtkm::Id nBranches = branchSaddle.GetNumberOfValues();
  std::vector<Branch<T>*> branches;
  Branch<T>* root = nullptr;
  branches.reserve(static_cast<std::size_t>(nBranches));


  std::cout << "Number of Branches, for the Branch Decomposition:" << nBranches << std::endl;

  for (int branchID = 0; branchID < nBranches; ++branchID)
    branches.push_back(new Branch<T>);

  // Reconstruct explicit branch decomposition from array representation
  for (std::size_t branchID = 0; branchID < static_cast<std::size_t>(nBranches); ++branchID)
  {
    branches[branchID]->OriginalId = static_cast<vtkm::Id>(branchID);
    if (!NoSuchElement(branchSaddlePortal.Get(static_cast<vtkm::Id>(branchID))))
    {
      branches[branchID]->Saddle = MaskedIndex(
        supernodesPortal.Get(MaskedIndex(branchSaddlePortal.Get(static_cast<vtkm::Id>(branchID)))));
      vtkm::Id branchMin = MaskedIndex(supernodesPortal.Get(
        MaskedIndex(branchMinimumPortal.Get(static_cast<vtkm::Id>(branchID)))));
      vtkm::Id branchMax = MaskedIndex(supernodesPortal.Get(
        MaskedIndex(branchMaximumPortal.Get(static_cast<vtkm::Id>(branchID)))));
      if (branchMin < branches[branchID]->Saddle)
        branches[branchID]->Extremum = branchMin;
      else if (branchMax > branches[branchID]->Saddle)
        branches[branchID]->Extremum = branchMax;
      else
      {
        std::cerr << "Internal error";
        return 0;
      }
    }
    else
    {
      branches[branchID]->Saddle =
        supernodesPortal.Get(MaskedIndex(branchMinimumPortal.Get(static_cast<vtkm::Id>(branchID))));
      branches[branchID]->Extremum =
        supernodesPortal.Get(MaskedIndex(branchMaximumPortal.Get(static_cast<vtkm::Id>(branchID))));
    }

    if (dataFieldIsSorted)
    {
      branches[branchID]->SaddleVal = dataFieldPortal.Get(branches[branchID]->Saddle);
      branches[branchID]->ExtremumVal = dataFieldPortal.Get(branches[branchID]->Extremum);
    }
    else
    {
      branches[branchID]->SaddleVal =
        dataFieldPortal.Get(sortOrderPortal.Get(branches[branchID]->Saddle));
      branches[branchID]->ExtremumVal =
        dataFieldPortal.Get(sortOrderPortal.Get(branches[branchID]->Extremum));
    }

    branches[branchID]->Saddle = sortOrderPortal.Get(branches[branchID]->Saddle);
    branches[branchID]->Extremum = sortOrderPortal.Get(branches[branchID]->Extremum);

    if (NoSuchElement(branchParentPortal.Get(static_cast<vtkm::Id>(branchID))))
    {
      root = branches[branchID]; // No parent -> this is the root branch
    }
    else
    {
      branches[branchID]->Parent = branches[static_cast<size_t>(
        MaskedIndex(branchParentPortal.Get(static_cast<vtkm::Id>(branchID))))];
      branches[branchID]->Parent->Children.push_back(branches[branchID]);
    }
  }

//  // FIXME: This is a somewhat hackish way to compute the Volume, but it works
//  // It would probably be better to compute this from the already computed Volume information
//  // (already replaced, 2025-03-06 commented out until floats)
//  std::cout << "Computing Integer Volumes" << std::endl;
  auto whichBranchPortal = whichBranch.ReadPortal();
  auto superparentsPortal = contourTreeSuperparents.ReadPortal();
//  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
//  {
//    branches[static_cast<size_t>(
//               MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]
//      ->Volume++; // Increment Volume

//    std::cout << "branch[" << static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))
//              << "]" << branches[static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]->Volume
//              << std::endl;
//  }

//  std::cout << std::endl;

//  // 2025-01-13
//  std::cout << "------------- vvv Branch INTEGER weights vvv -------------" << std::endl;

//  // loop through all the regular nodes, ...
//  // ... then counting how many regular nodes are on each branch
//  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
//  {
//    size_t branchID = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))));
//    //      vtkm::Id sortID = supernodesPortal.Get(i);
//    //    std::cout << i << " [sortID=" << sortID << "]" << " (" << MaskedIndex(superparentsPortal.Get(i)) << ") "
//    std::cout << i << "[" << branchID << "] " << " (" << MaskedIndex(superparentsPortal.Get(i)) << ") "
//              << branches[branchID]->Extremum << "->" << branches[branchID]->Saddle << " = "
//              << branches[branchID]->Volume
//              << std::endl;
//  }
//  std::cout << "------------- ^^^ Branch INTEGER weights ^^^ -------------" << std::endl;

//  std::cout << std::endl;

  vtkm::Id sortID; // = supernodesPortal.Get(i);
  int current_superparent;
  int supernode_tailend;
  std::vector<std::vector<vtkm::Id>> branch_SP_map(nBranches);

// 2025-03-10 Commented out the branch initialisation, replaced with a branch-parent-array
#if DEBUG_PRINT_PACTBD
  std::cout << "------- vvv Branch VALUE TYPE INITIAL weights vvv --------" << std::endl;
#endif
  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
  {
    current_superparent = MaskedIndex(superparentsPortal.Get(i));
//    sortID = supernodesPortal.Get(i);
    supernode_tailend = supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(current_superparent)));

    size_t branchID = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(current_superparent)));

    vtkm::Id branchMin = MaskedIndex(supernodesPortal.Get(
      MaskedIndex(branchMinimumPortal.Get(static_cast<vtkm::Id>(branchID)))));
    vtkm::Id branchMax = MaskedIndex(supernodesPortal.Get(
      MaskedIndex(branchMaximumPortal.Get(static_cast<vtkm::Id>(branchID)))));

//    std::cout << MaskedIndex(superparentsPortal.Get(i)) << ") "
//              << branches[branchID]->Extremum << "->" << branches[branchID]->Saddle << " = "
//              << branches[branchID]->VolumeFloat
//              << std::endl;

#if DEBUG_PRINT_PACTBD
    std::cout << i << "[" << branchID << "] " << " (" << current_superparent << ") "
              << branches[branchID]->Extremum << "->" << branches[branchID]->Saddle << " = "
              << branches[branchID]->VolumeFloat
              << std::endl;
#endif
    if (std::find(branch_SP_map[branchID].begin(), branch_SP_map[branchID].end(), current_superparent) == branch_SP_map[branchID].end())
    {
#if DEBUG_PRINT_PACTBD
      // 2025-03-10
      std::cout << "-> ADDED BRANCH CHAIN SP: " << current_superparent << std::endl;
#endif
      branch_SP_map[branchID].push_back(current_superparent); // supernode_tailend
      branch_SP_map[branchID].push_back(i); // NEW 2025-03-14
      branch_SP_map[branchID].push_back(supernode_tailend); // NEW 2025-03-14
    }


  }
#if DEBUG_PRINT_PACTBD
  std::cout << "------- ^^^ Branch VALUE TYPE INITIAL weights ^^^ --------" << std::endl << std::endl;
  std::cout << "Printing the supernode/branch mappings" << std::endl;
#endif

//  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
//  {

//      std::cout << i << ")" << MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))) << std::endl;

//    branches[static_cast<size_t>(
//               MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]
//      ->Volume++; // Increment Volume
//  }
#if DEBUG_PRINT_PACTBD
std::cout << "Printing the supernode/branch mappings" << std::endl;
#endif
  for(int i = 0; i < nBranches; i++)
  {
#if DEBUG_PRINT_PACTBD
      std::cout << "branch[" << i << "] "  << "\t";
#endif

#if DEBUG_PRINT_PACTBD
    std::cout << i << " -> "
              << branches[i]->Extremum << "->" << branches[i]->Saddle << " = "
              << branches[i]->VolumeFloat << std::endl << "\t";
#endif

    for(int j = 0; j < branch_SP_map[i].size(); j+=3) //j++)
    {
#if DEBUG_PRINT_PACTBD
        std::cout << branch_SP_map[i][j] << " ";
#endif
//        if (!std::isnan(superarcIntrinsicWeightPortal.Get(branch_SP_map[i][j])))
        branches[i]->VolumeFloat += superarcIntrinsicWeightPortal.Get(branch_SP_map[i][j]);

#if DEBUG_PRINT_PACTBD
//        std::cout << "(" << branch_SP_map[i][j] << " = " << branch_SP_map[i][j+1] << " -> " << branch_SP_map[i][j+2] << ") ";
        std::cout << branch_SP_map[i][j] << " -> ";
#endif

#if DEBUG_PRINT_PACTBD
        size_t branchIDPrint = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(branch_SP_map[i][j])));

        if(!branches[branchIDPrint]->Children.empty())
        {
            std::cout << std::endl << "\t";
            PrintBranchInformation(branches[branchIDPrint]->Children[0], branch_SP_map, branches[branchIDPrint]->OriginalId);
        }

                for (Branch<T>* c : branches[branchIDPrint]->Children)
                {
                    PrintBranchInformation(c, branch_SP_map, c->OriginalId);
                    std::cout << c->OriginalId << " -> ";
                }

        if(!branches[i]->Children.empty())
        {
            std::cout << "Children of i=" << i << std::endl;
        }
        else
        {
            std::cout << std::endl;
        }

        for (Branch<T>* c : branches[i]->Children)
        {
//            PrintBranchInformation(c, branch_SP_map, c->OriginalId);
            std::cout << c->OriginalId << " -> ";
        }
#endif


    }

#if DEBUG_PRINT_PACTBD
    if(!branches[i]->Children.empty())
    {
        std::cout << std::endl << "Children of i=" << i << std::endl;
    }

    for (Branch<T>* c : branches[i]->Children)
    {
        std::cout << c->OriginalId << " -> ";
//                    std::cout << std::endl << "\t";
            PrintBranchInformation(c, branch_SP_map, c->OriginalId);
    }
#endif


//    size_t branchIDPrint = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(branch_SP_map[i][0])));

//    if(!branches[branchIDPrint]->Children.empty())
//    {
//        std::cout << std::endl << "\t";
//    }

//    for (Branch<T>* c : branches[branchIDPrint]->Children)
//    {
//        PrintBranchInformation(c, branch_SP_map, branchIDPrint);
////        std::cout << " (" << c->OriginalId << ") ";
//    }

//    if(!branches[i]->Children.empty())
//    {
//        std::cout << std::endl << "\t";
//    }

//// 2025-03-15 NEW
//    for (Branch<T>* c : branches[i]->Children)
//    {
//        PrintBranchInformation(c);
////        std::cout << c->OriginalId << " -> ";
//    }

#if DEBUG_PRINT_PACTBD
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << i << " -> "
              << branches[i]->Extremum << "->" << branches[i]->Saddle << " = "
              << branches[i]->VolumeFloat << std::endl;
#endif

  }

//  // 2025-03-10 NO LONGER DOING THIS:
//  int previous_superparent = MaskedIndex(superparentsPortal.Get(0));
//  //  int current_superparent;

//  int previous_branchID = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(0)))));
//  int current_branchID;
//  bool last_iteration = false;

//  std::cout << std::endl;
//  std::cout << "------------- vvv Weight computation by summing branch intrinsic VALUE TYPE weights vvv -------------" << std::endl;
//  // loop through all the regular nodes, ...
//  // ... then counting how many regular nodes are on each branch
//  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
//  {
//    if(i+1 == contourTreeSuperparents.GetNumberOfValues())
//    {
//        last_iteration = true;
//    }
//    current_superparent = MaskedIndex(superparentsPortal.Get(i));

//    size_t branchID = static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(current_superparent)));

//    std::cout << i << "[" << branchID << "] " << " (" << current_superparent << ") "
//              << branches[branchID]->Extremum << "->" << branches[branchID]->Saddle << " = "
//              << branches[branchID]->VolumeFloat
//              << std::endl;

//    if(current_superparent != previous_superparent)
//    {
//        std::cout << "Current superparent: " << current_superparent << " vs Previous: " << previous_superparent << std::endl;
//        std::cout << "\tbranch[" << previous_branchID << "] = " << branches[previous_branchID]->VolumeFloat << std::endl;

//        std::cout << "\t\tbranch[" << previous_branchID << "] += " << superarcIntrinsicWeightPortal.Get(previous_superparent) << std::endl;

//        branches[previous_branchID]->VolumeFloat += superarcIntrinsicWeightPortal.Get(previous_superparent);

//        std::cout << "\tbranch[" << previous_branchID << "] = " << branches[previous_branchID]->VolumeFloat << std::endl << std::endl;
//    }


//    previous_superparent = current_superparent;
//    previous_branchID = branchID;

//    if(last_iteration)
//    {
//        std::cout << "\tbranch[" << previous_branchID << "] = " << branches[previous_branchID]->VolumeFloat << std::endl;

//        std::cout << "\t\tbranch[" << previous_branchID << "] += " << superarcIntrinsicWeightPortal.Get(previous_superparent) << std::endl;

//        branches[previous_branchID]->VolumeFloat += superarcIntrinsicWeightPortal.Get(previous_superparent);

//        std::cout << "\tbranch[" << previous_branchID << "] = " << branches[previous_branchID]->VolumeFloat << std::endl << std::endl;
//    }


//  }
//  std::cout << "------------- ^^^ Weight computation by summing branch intrinsic VALUE TYPE weights ^^^ -------------" << std::endl;
//  std::cout << std::endl;


//  // 2025-03-10 REMOVE THE NODE-COUNTING CODE BELOW:
//  // 2024-08-16 Get the existing volume information instead of just counting nodes

//  std::cout << "Number of SPs: " << contourTreeSuperparents.GetNumberOfValues() << std::endl;

//  // 2025-01-05 getting the right (floating point) volume information
//  for (vtkm::Id i = 0; i < contourTreeSuperparents.GetNumberOfValues(); i++)
//  {
////    branches[static_cast<size_t>(
////               MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]
////      ->VolumeFloat = 0.5f; // Increment Volume

////    vtkm::Id sortID = supernodesPortal.Get(i);

////    // retrieve ID of target supernode
//////    vtkm::Id superTo = superarcsPortal.Get(supernode);

////    // if this is true, it is the last pruned vertex & is omitted
////    if (NoSuchElement(superTo))
////      continue;

////    // otherwise, strip out the flags
////    superTo = MaskedIndex(superTo);

//    // otherwise, we need to convert the IDs to regular mesh IDs
////    vtkm::Id regularID = sortOrderPortal.Get(MaskedIndex(sortID));

////    std::cout << "'Incremented' Volume ... of [Branch " << MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i))))
////              << "] of SP: (" << superparentsPortal.Get(i) << ") = " << std::endl; // << regularID << std::endl;

//    std::cout << "branch[" << static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i))))) << "] = "
//              << branches[static_cast<size_t>(MaskedIndex(whichBranchPortal.Get(MaskedIndex(superparentsPortal.Get(i)))))]->VolumeFloat
//              << std::endl;

//  }



  if (root)
  {
    root->removeSymbolicPerturbation();
  }

  std::cout << "##################################################################" << std::endl;
  std::cout << "############################ FINISH ##############################" << std::endl;
  std::cout << "############# !MODIFIED! Branch.H Decomposition ##################" << std::endl;
  std::cout << "##################################################################" << std::endl;
  std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;

  std::cout << "ContourTreeApp->ProcessContourTree->(Branch.h->ComputeBranchDecomposition())" << std::endl;

  return root;
} // ComputeBranchDecomposition()





template <typename T>
// print the graph in python dict format:
static void PrintBranchInformation(Branch<T>* root)
{
    std::cout<< root->OriginalId << " -> ";
}
























template <typename T>
void Branch<T>::SimplifyToSize(vtkm::Id targetSize, bool usePersistenceSorter)
{ // SimplifyToSize()

  std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~ Branch.h SimplifyToSize ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;



  std::cout << "-Branch.h-> Calling: SimplifyToSize() with targetSize = " << targetSize << std::endl;

  if (targetSize <= 1)
    return;

  // Top-down simplification:
  // starting from one branch (the main one) ...
  // ... and adding in the rest on a biggest-first basis
  std::vector<Branch<T>*> q;
  q.push_back(this);



  // Print the branch (root) information - using 'this' to get the (main) branch this function was called from:
  vtkm::Id OriginalId;              // Index of the extremum in the mesh
  vtkm::Id Extremum;                // Index of the extremum in the mesh
  T ExtremumVal;                    // Value at the extremum:w
  vtkm::Id Saddle;                  // Index of the saddle in the mesh (or minimum for root branch)
  T SaddleVal;                      // Corresponding value
  vtkm::Id Volume;                  // Volume
  ValueType VolumeFloat;
  Branch<T>* Parent;                // Pointer to parent, or nullptr if no parent
  std::vector<Branch<T>*> Children; // List of pointers to children
  std::cout << "------------------branchDecompostionRoot-----------------" << std::endl;
  std::cout << "branch ID: "            << this->OriginalId  << std::endl;      // Index of the extremum in the mesh
  std::cout << "branch Extremum: "      << this->Extremum    << std::endl;      // Index of the extremum in the mesh
  std::cout << "branch ExtremumVal: "   << this->ExtremumVal << std::endl;      // Value at the extremum:w
  std::cout << "branch Saddle: "        << this->Saddle      << std::endl;      // Index of the saddle in the mesh (or minimum for root branch)
  std::cout << "branch SaddleVal: "     << this->SaddleVal   << std::endl;      // Corresponding value
  std::cout << "branch Volume: "        << this->Volume      << std::endl;      // Volume
  std::cout << "branch VolumeFloat: "   << this->VolumeFloat << std::endl;
//  Branch<T>* Parent;                // Pointer to parent, or nullptr if no parent
//  std::vector<Branch<T>*> Children; // List of pointers to children
  std::cout << "------------------branchDecompostionRoot-----------------\n" << std::endl;


  int num_active_branches = 0;

  // local array - active definition:
  std::vector<Branch<T>*> active;
  while (active.size() < static_cast<std::size_t>(targetSize) && !q.empty())
  {
    std::cout << "====================== iteration " << num_active_branches << " ======================" << std::endl;
    std::cout << "-Branch.h-> active.size() = " << active.size() << std::endl;

    if (usePersistenceSorter)
    {
      std::cout << "-Branch.h-> usePersistenceSorter = " << usePersistenceSorter<< std::endl;
      std::pop_heap(
        q.begin(),
        q.end(),
        PersistenceSorter<
          T>()); // FIXME: This should be Volume, but we were doing this wrong for the demo, so let's start with doing this wrong here, too
    }
    else
    {
      std::cout << "-Branch.h-> VolumeSorter: usePersistenceSorter = " << usePersistenceSorter<< std::endl;
      std::pop_heap(
        q.begin(),
        q.end(),
        // VolumeSorter<
        VolumeSorterFloat<
          T>()); // FIXME: This should be Volume, but we were doing this wrong for the demo, so let's start with doing this wrong here, too

        //      std::push_heap(q.begin(), q.end(), VolumeSorterFloat<T>());

    }
    Branch<T>* b = q.back();
    q.pop_back();

    std::cout << "----------> (Active) Processing Branch: " << b->ExtremumVal << std::endl;

    active.push_back(b);

    std::cout << "-Branch.h-> active.size() = " << active.size() << std::endl;

    std::cout << "\n---------------------------------------------------------" << std::endl;
    std::cout << "branch ID: "            << b->OriginalId      << std::endl;      // Index of the extremum in the mesh
    std::cout << "branch Extremum: "      << b->Extremum        << std::endl;      // Index of the extremum in the mesh
    std::cout << "branch ExtremumVal: "   << b->ExtremumVal     << std::endl;      // Value at the extremum:w
    std::cout << "branch Saddle: "        << b->Saddle          << std::endl;      // Index of the saddle in the mesh (or minimum for root branch)
    std::cout << "branch SaddleVal: "     << b->SaddleVal       << std::endl;      // Corresponding value
    std::cout << "branch Volume: "        << b->Volume          << std::endl;      // Volume
    std::cout << "branch VolumeFloat: "   << b->VolumeFloat     << std::endl;
    std::cout << "branch Children: "      << b->Children.size() << std::endl;
  //  Branch<T>* Parent;                // Pointer to parent, or nullptr if no parent
  //  std::vector<Branch<T>*> Children; // List of pointers to children
    std::cout << "---------------------------------------------------------\n" << std::endl;

    int child_num = 0;

    for (Branch<T>* c : b->Children)
    {
#if DEBUG_PRINT_PACTBD
      std::cout << "---------> (Child " << child_num << ") Processing Branch: " << c->ExtremumVal << std::endl;
#endif

      q.push_back(c);

#if DEBUG_PRINT_PACTBD
      for (auto e : q)
      {
          std::cout << e->ExtremumVal << ' ';
          std::cout << '\n';
      }
#endif

      if (usePersistenceSorter)
      {
        std::push_heap(q.begin(), q.end(), PersistenceSorter<T>());
      }
      else
      {
#if DEBUG_PRINT_PACTBD
        std::cout << "\t-> pushing child" << std::endl;
#endif
        // std::push_heap(q.begin(), q.end(), VolumeSorter<T>());
        std::push_heap(q.begin(), q.end(), VolumeSorterFloat<T>());
      }

#if DEBUG_PRINT_PACTBD
      for (auto e : q)
      {
          std::cout << e->ExtremumVal << ' ';
          std::cout << '\n';
      }
#endif

      child_num++;

    }
#if DEBUG_PRINT_PACTBD
    std::cout << "====================== iteration " << num_active_branches << " ======================" << std::endl;
#endif

    num_active_branches++;

  }

  // Rest are inactive
  for (Branch<T>* b : q)
  {
    // Hackish, remove c from its parents child list
    if (b->Parent)
      b->Parent->Children.erase(
        std::remove(b->Parent->Children.begin(), b->Parent->Children.end(), b));

    delete b;
  }

  std::cout << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FINISH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~ Branch.h SimplifyToSize ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << std::endl << std::endl;

} // SimplifyToSize()



template <typename T>
// print the graph in python dict format:
void Branch<T>::PrintBranchDecomposition(std::ostream& os, std::string::size_type indent) const
{ // PrintBranchDecomposition()

  os << std::string(indent, ' ') << "{" << std::endl;
  os << std::string(indent, ' ') << "  'Saddle' : " << SaddleVal << ","
     << std::endl;
  os << std::string(indent, ' ') << "  'Extremum' : " << ExtremumVal << ","
     << std::endl;
  os << std::string(indent, ' ') << "  'Volume' : " << Volume << ","<< std::endl;
  os << std::string(indent, ' ') << "  'VolumeFloat' : " << VolumeFloat << ","<< std::endl;

  if (!Children.empty())
  {
    os << std::string(indent, ' ') << "  'Children' : [" << std::endl;
    for (Branch<T>* c : Children)
    {
      c->PrintBranchDecomposition(os, indent + 4);
    }
    os << std::string(indent, ' ') << std::string(indent, ' ') << "  ]," << std::endl;
  }
  os << std::string(indent, ' ') << "}," << std::endl;
} // PrintBranchDecomposition()








template <typename T>
// print the graph in dot (.gv) format
void Branch<T>::PrintDotBranchDecomposition(std::ostream& os,
                                            std::vector<vtkm::Id>& saddles,
                                            std::vector<vtkm::Id>& local_branches,
                                            std::vector<vtkm::Id>& depth,
                                            std::vector<vtkm::FloatDefault>& branch_weights,
                                            std::vector<vtkm::FloatDefault>& branch_weights_write,
                                            std::vector<bool>& main_branch_flags,
                                            std::vector<vtkm::Id>& depths_write,
                                            vtkm::Id parent_saddle,
                                            vtkm::Id parent_extremum,
                                            int iteration,
                                            std::string::size_type indent)
{ // PrintBranchDecomposition()







    // NEW: use the data values passed into the field here:
    // NEW PACTBD-EDIT
    //    int num_datapoints = 101;
        int num_datapoints = 1001;
    //    int num_datapoints = 10001;
//        int num_datapoints = 99972;
//        int num_datapoints = 200001;
//        int num_datapoints = 985181;
//        int num_datapoints = 2160930;
//      const std::string field_filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/101-field.txt";
//    const std::string field_filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/10k-field.txt";
//      const std::string field_filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/1M-field.txt";
//      const std::string field_filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/2M-parcels-20250225-field-sorted.txt";

    // ARCHER2
//    const std::string field_filename = "/work/e710/e710/ddilys/PACTBD/data/10k-field.txt";
    const std::string field_filename = "/work/e710/e710/ddilys/PACTBD/data/100k-field.txt";


    std::ifstream field_input(field_filename);
    vtkm::cont::ArrayHandle<vtkm::Float64> fakeFieldArray;
    fakeFieldArray.Allocate(num_datapoints);
    auto fakeFieldArrayWritePortal = fakeFieldArray.WritePortal();
    auto fakeFieldArrayReadPortal = fakeFieldArray.ReadPortal();



    std::vector<vtkm::Float64> std_field;
    if(field_input.is_open())
    {
        std::string line;
        int i = 0;
        while(getline(field_input, line))
        {
            std_field.push_back(static_cast<vtkm::Float64>(std::stof(line)));
        }

        for(vtkm::Id i = 0; i < num_datapoints; i++)
        {
          fakeFieldArrayWritePortal.Set(i, std_field[i]);
        }
    }
    else
    {
        std::cerr << "Unable to open file: " << field_filename << "\n";
    }
    field_input.close();










  std::string tab = "\t";
  bool write_triggerGV = false;

  if(iteration == 0)
  {// if we are in the root call of recursion ...
   // ... trigger the writing to file flag ...
   // ... and write the header to the file

    write_triggerGV = true;

    os << "digraph G" << std::endl;
    os << tab << "{" << std::endl;
    os << tab << "size=\"6.5, 9\"" << std::endl;
    os << tab << "ratio=\"fill\"" << std::endl;

  }

  std::vector<vtkm::Id> saddles_local; // saddle array (nodes on the main branch)
  std::vector<vtkm::Id> nodes;

  if (!Children.empty())
  {// HAS children (a or b cases)

    local_branches.push_back(SaddleVal);
    depth.push_back(iteration);
    branch_weights.push_back(VolumeFloat);


    for (Branch<T>* c : Children)
    {
      c->PrintDotBranchDecomposition(os, saddles, local_branches, depth, branch_weights, branch_weights_write, main_branch_flags, depths_write,
                                     SaddleVal, ExtremumVal, iteration+1, indent + 4);
    }
#if DEBUG_PRINT_PACTBD
    std::cout << "Finished branch segment: " << SaddleVal << " " << ExtremumVal << " " << local_branches.size() << " it: " << iteration+1 << std::endl;
#endif

    std::vector<vtkm::Id> local_iteration_branches;
    for(int i = 0; i < local_branches.size(); i++)
    {
#if DEBUG_PRINT_PACTBD
        std::cout << depth[i] << " ";
#endif
        if(depth[i] == iteration+1)
        { // only copy branches from the same depth
           local_iteration_branches.push_back(local_branches[i]);
        }
    }
#if DEBUG_PRINT_PACTBD
    std::cout << std::endl << "local: " << local_iteration_branches.size() << std::endl;
#endif

    std::sort(local_iteration_branches.begin(), local_iteration_branches.end());

    int for_iterator = 0;
    for (auto c : local_iteration_branches)
    {
        if(for_iterator)
        {
#if DEBUG_PRINT_PACTBD
            std::cout << c << "->" << parent_saddle << "[" << iteration+1 << "a" << for_iterator << "] (ORANGE)" << std::endl;
#endif
            saddles.push_back(c);
            saddles.push_back(parent_saddle);
            branch_weights_write.push_back(VolumeFloat);
            branch_weights_write.push_back(VolumeFloat);
            main_branch_flags.push_back(write_triggerGV);
            main_branch_flags.push_back(write_triggerGV);
            depths_write.push_back(iteration);
            depths_write.push_back(iteration);
        }
        else
        {// if it's the first iteration then the parent_saddle is the original
#if DEBUG_PRINT_PACTBD
            std::cout << c << "->" << SaddleVal << "[" << iteration+1 << "a" << for_iterator << "] (ORANGE)" << std::endl;
#endif
            saddles.push_back(c);
            saddles.push_back(SaddleVal);
            branch_weights_write.push_back(VolumeFloat);
            branch_weights_write.push_back(VolumeFloat);
            main_branch_flags.push_back(write_triggerGV);
            main_branch_flags.push_back(write_triggerGV);
            depths_write.push_back(iteration);
            depths_write.push_back(iteration);
        }

        parent_saddle = c;
        for_iterator++;
    }
#if DEBUG_PRINT_PACTBD
    std::cout << ExtremumVal << "->" << local_iteration_branches[for_iterator-1] << "[" << iteration+1 << "b] (RED)" << std::endl;
#endif
    saddles.push_back(ExtremumVal);
    saddles.push_back(local_iteration_branches[for_iterator-1]);
    branch_weights_write.push_back(VolumeFloat);
    branch_weights_write.push_back(VolumeFloat);
    main_branch_flags.push_back(write_triggerGV);
    main_branch_flags.push_back(write_triggerGV);
    depths_write.push_back(iteration);
    depths_write.push_back(iteration);

  }
  else
  {// HAS NO children (pure cases Extremum -> Saddle)
#if DEBUG_PRINT_PACTBD
      std::cout << ExtremumVal << "->" << SaddleVal << "[" << iteration << "] (GREEN)" << std::endl;
#endif

      saddles.push_back(ExtremumVal);
      saddles.push_back(SaddleVal);

      branch_weights_write.push_back(VolumeFloat);
      branch_weights_write.push_back(VolumeFloat);

      main_branch_flags.push_back(write_triggerGV);
      main_branch_flags.push_back(write_triggerGV);

      depths_write.push_back(iteration);
      depths_write.push_back(iteration);

      local_branches.push_back(SaddleVal);
      depth.push_back(iteration);

      // hack - push the weights for both nodes representing a branch
      branch_weights.push_back(VolumeFloat);
  }


  if (write_triggerGV)
  { // write from the root recursion call (main branch)
      std::cout << "Writing to graphviz file" << std::endl;

      std::vector<vtkm::Id> nodes;
      std::vector<vtkm::Id> depths;

//      for(vtkm::Id nodeID : saddles)
      for(int nodeID = 0; nodeID < saddles.size(); nodeID++) // : saddles)
      {
          if (std::find(nodes.begin(), nodes.end(), saddles[nodeID]) == nodes.end())
          {// if saddle value not already in saddles
            nodes.push_back(saddles[nodeID]);
            depths.push_back(depths_write[nodeID]);
          }
      }

      std::map<int, std::string> colourMap;
      std::string colour_label;

      std::cout << "... Nodes gathered" << std::endl;

      for(vtkm::Id iterator = 0; iterator < saddles.size()-1; iterator+=2)
      {
          if(main_branch_flags[iterator])
          {
              colour_label = "green";
          }
          else
          {
              colour_label = "red";
          }

          colourMap.insert({saddles[iterator],   colour_label});
          colourMap.insert({saddles[iterator+1], colour_label});
      }

      std::cout << "... Mappings made" << std::endl;

      int depth_iter = 0;
//      for(vtkm::Id nodeID : nodes)
      for(int nodeID = 0; nodeID < nodes.size(); nodeID++) // : saddles)
      {
#if DEBUG_PRINT_PACTBD
          std::cout << "s" << nodes[nodeID] << std::endl;
#endif

//          if(main_branch_flags[iterator])


          if(depths[nodeID])
          {
              os << tab << "s" << nodes[nodeID] << "[style=filled,fillcolor=" << colourMap[nodes[nodeID]] << ", label=\"" << fakeFieldArrayReadPortal.Get(nodes[nodeID]) << "\"]" << std::endl;
          }
          else
          {
//              os << tab << "s" << fakeFieldArrayReadPortal.Get(nodes[nodeID]) << "[style=filled,fillcolor=" << colourMap[nodes[nodeID]] << "]" << std::endl;
              os << tab << "s" << nodes[nodeID] << "[style=filled,fillcolor=" << colourMap[nodes[nodeID]] << ", label=\"" << fakeFieldArrayReadPortal.Get(nodes[nodeID]) << "\"]" << std::endl;
//              os << tab << "s" << nodes[nodeID] << "[style=filled,fillcolor=" << colourMap[nodes[nodeID]] << "]" << std::endl;
          }

          depth_iter++;
      }
#if DEBUG_PRINT_PACTBD
      std::cout << std::endl;
#endif
      int local_iterator = 0;
//      for(vtkm::Id nodeID : saddles)
      for(vtkm::Id iterator = 0; iterator < saddles.size()-1; iterator+=2)
      {
#if DEBUG_PRINT_PACTBD
          std::cout << "s" << saddles[iterator] << "->" << "s" << saddles[iterator+1] << std::endl;
#endif
          if(main_branch_flags[iterator])
          {
              os << tab << "s" << saddles[iterator] << " -> " << "s" << saddles[iterator+1] << "[label=\"(main) " << branch_weights_write[iterator] << "\"]" << std::endl;
          }
          else
          {
              os << tab << "s" << saddles[iterator] << " -> " << "s" << saddles[iterator+1] << "[label=\"" << branch_weights_write[iterator] << "\"]" << std::endl;
          }

          local_iterator++;
      }
#if DEBUG_PRINT_PACTBD
      std::cout << std::endl;
#endif
      os << tab << "}" << std::endl;
  }

} // PrintBranchDecomposition()





// OLD FORMAT
//void Branch<T>::PrintBranchDecomposition(std::ostream& os, std::string::size_type indent) const
//{ // PrintBranchDecomposition()
//  os << std::string(indent, ' ') << "{" << std::endl;
//  os << std::string(indent, ' ') << "  Saddle = " << SaddleVal << " (" << Saddle << ")"
//     << std::endl;
//  os << std::string(indent, ' ') << "  Extremum = " << ExtremumVal << " (" << Extremum << ")"
//     << std::endl;
//  os << std::string(indent, ' ') << "  Volume = " << Volume << std::endl;
//  if (!Children.empty())
//  {
//    os << std::string(indent, ' ') << "  Children = [" << std::endl;
//    for (Branch<T>* c : Children)
//      c->PrintBranchDecomposition(os, indent + 4);
//    os << std::string(indent, ' ') << std::string(indent, ' ') << "  ]" << std::endl;
//  }
//  os << std::string(indent, ' ') << "}" << std::endl;
//} // PrintBranchDecomposition()



template <typename T>
Branch<T>::~Branch()
{ // ~Branch()
  for (Branch<T>* c : Children)
    delete c;
  if (Parent)
  {
    Parent->Volume += Volume;
    Parent->VolumeFloat += VolumeFloat;
  }
} // ~Branch()


// TODO this recursive accumlation of values does not lend itself well to the use of VTKM data structures
template <typename T>
void Branch<T>::GetRelevantValues(int type, T eps, std::vector<T>& values) const
{ // GetRelevantValues()
  T val;

  bool isMax = false;
  if (ExtremumVal > SaddleVal)
    isMax = true;

  switch (type)
  {
    default:
    case 0:
      val = SaddleVal + (isMax ? +eps : -eps);
      break;
    case 1:
      val = T(0.5f) * (ExtremumVal + SaddleVal);
      break;
    case 2:
      val = ExtremumVal + (isMax ? -eps : +eps);
      break;
  }
  if (Parent)
    values.push_back({ val });
  for (Branch* c : Children)
    c->GetRelevantValues(type, eps, values);
} // GetRelevantValues()


template <typename T>
void Branch<T>::AccumulateIntervals(int type, T eps, PiecewiseLinearFunction<T>& plf) const
{ //AccumulateIntervals()
  bool isMax = (ExtremumVal > SaddleVal);
  T val;

  switch (type)
  {
    default:
    case 0:
      val = SaddleVal + (isMax ? +eps : -eps);
      break;
    case 1:
      val = T(0.5f) * (ExtremumVal + SaddleVal);
      break;
    case 2:
      val = ExtremumVal + (isMax ? -eps : +eps);
      break;
  }

  if (Parent)
  {
    PiecewiseLinearFunction<T> addPLF;
    addPLF.addSample(SaddleVal, 0.0);
    addPLF.addSample(ExtremumVal, 0.0);
    addPLF.addSample(val, 1.0);
    plf += addPLF;
  }
  for (Branch<T>* c : Children)
    c->AccumulateIntervals(type, eps, plf);
} // AccumulateIntervals()


template <typename T>
void Branch<T>::removeSymbolicPerturbation()
{                                      // removeSymbolicPerturbation()
  std::vector<Branch<T>*> newChildren; // Temporary list of children that are not flat

  for (Branch<T>* c : Children)
  {
    // First recursively remove symbolic perturbation (zero persistence branches) for  all children below the current child
    // Necessary to be able to detect whether we can remove the current child
    c->removeSymbolicPerturbation();

    // Does child have zero persistence (flat region)
    if (c->ExtremumVal == c->SaddleVal && c->Children.empty())
    {
      // If yes, then we get its associated Volume and delete it
      delete c; // Will add Volume to parent, i.e., us
    }
    else
    {
      // Otherwise, keep child
      newChildren.push_back(c);
    }
  }
  // Swap out new list of children
  Children.swap(newChildren);
} // removeSymbolicPerturbation()

} // process_contourtree_inc
} // namespace contourtree_augmented
} // namespace worklet
} // namespace vtkm

#endif // vtk_m_worklet_contourtree_augmented_process_contourtree_inc_branch_h
