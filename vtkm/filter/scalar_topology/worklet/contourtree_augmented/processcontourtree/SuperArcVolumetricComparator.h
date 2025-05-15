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

#ifndef vtk_m_worklet_contourtree_augmented_process_contourtree_inc_superarc_volumetric_comperator_h
#define vtk_m_worklet_contourtree_augmented_process_contourtree_inc_superarc_volumetric_comperator_h

#include <vtkm/Pair.h>
#include <vtkm/Types.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ExecutionObjectBase.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>

namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{
namespace process_contourtree_inc
{

class SuperArcVolumetricComparatorImpl
{ // SuperArcVolumetricComparatorImpl
public:
  using IdPortalType = vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType;
  using EdgePairArrayPortalType = EdgePairArray::ReadPortalType;

//    using ValueType = vtkm::Float32;
using ValueType = vtkm::Float64; //vtkm::FloatDefault;
    using FloatArrayType = vtkm::cont::ArrayHandle<ValueType>;
    using FloatPortalType = vtkm::cont::ArrayHandle<ValueType>::ReadPortalType;

  IdPortalType weightPortal;
  FloatPortalType weightFloatPortal;

  bool pairsAtLowEnd;
  EdgePairArrayPortalType superarcListPortal;

  // constructor
  SuperArcVolumetricComparatorImpl(const IdArrayType& Weight,
                                   const EdgePairArray& SuperarcList,
                                   bool PairsAtLowEnd,
                                   vtkm::cont::DeviceAdapterId device,
                                   vtkm::cont::Token& token)
    : pairsAtLowEnd(PairsAtLowEnd)
  { // constructor
    weightPortal = Weight.PrepareForInput(device, token);
    superarcListPortal = SuperarcList.PrepareForInput(device, token);

    std::cout << "IdArrayType Comparator Implementation" << std::endl;

  } // constructor

  // float constructor
  SuperArcVolumetricComparatorImpl(const FloatArrayType& WeightFloat,
                                   const EdgePairArray& SuperarcList,
                                   bool PairsAtLowEnd,
                                   vtkm::cont::DeviceAdapterId device,
                                   vtkm::cont::Token& token)
    : pairsAtLowEnd(PairsAtLowEnd)
  { // constructor
    weightFloatPortal = WeightFloat.PrepareForInput(device, token);
    superarcListPortal = SuperarcList.PrepareForInput(device, token);

//    std::cout << "FloatArrayType Comparator Implementation" << std::endl;

  } // constructor

  // () operator - gets called to do comparison (vtkm::Id)
  VTKM_EXEC
  bool operator()(const vtkm::Id& i1, const vtkm::Id& i2) const
  { // operator()

    int intLen = weightPortal.GetNumberOfValues();
    int floatLen = weightFloatPortal.GetNumberOfValues();

//    std::cout << "vtkm::Id Comparator ()" << std::endl;
//    std::cout << "\tid    len: " << weightPortal.GetNumberOfValues() << std::endl;
//    std::cout << "\tfloat len: " << weightFloatPortal.GetNumberOfValues() << std::endl;
//    std::cout << "id1 val: " << weightPortal.Get(i1) << "id2 val: " << weightPortal.Get(i2) << std::endl;


    // get local references to the edge details
    EdgePair e1 = superarcListPortal.Get(i1);
    EdgePair e2 = superarcListPortal.Get(i2);

    if (pairsAtLowEnd)
    { // pairs at low end
      // test by low end ID
      if (e1.first < e2.first)
        return true;
      if (e1.first > e2.first)
        return false;

      if (intLen > floatLen)
      {
          // test by volumetric measure
          if (weightPortal.Get(i1) < weightPortal.Get(i2))
            return true;
          if (weightPortal.Get(i1) > weightPortal.Get(i2))
            return false;
      }

      if (floatLen > intLen)
      {
          // test by volumetric measure
          if (weightFloatPortal.Get(i1) < weightFloatPortal.Get(i2))
            return true;
          if (weightFloatPortal.Get(i1) > weightFloatPortal.Get(i2))
            return false;
      }



      // test by ID (persistence)
      if (e1.second < e2.second)
        return true;
      if (e1.second > e2.second)
        return false;

      // fallback
      return false;
    } // pairs at low end
    else
    { // pairs at high end
      // test by high end ID
      if (e1.second < e2.second)
        return true;
      if (e1.second > e2.second)
        return false;

      if (intLen > floatLen)
      {
          // test by volumetric measure
          if (weightPortal.Get(i1) < weightPortal.Get(i2))
            return true;
          if (weightPortal.Get(i1) > weightPortal.Get(i2))
            return false;
      }

      if (floatLen > intLen)
      {
          // test by volumetric measure
          if (weightFloatPortal.Get(i1) < weightFloatPortal.Get(i2))
            return true;
          if (weightFloatPortal.Get(i1) > weightFloatPortal.Get(i2))
            return false;
      }

      // test by ID (persistence)
      // Note the reversal from above - we want the greatest difference, not
      // the greatest value
      if (e1.first > e2.first)
        return true;
      if (e1.first < e2.first)
        return false;

      // fallback
      return false;
    } // pairs at high end
  }   // operator()
};    // SuperArcVolumetricComparatorImpl

class SuperArcVolumetricComparator : public vtkm::cont::ExecutionObjectBase
{ // SuperArcVolumetricComparator


//    using ValueType = vtkm::Float32;
using ValueType = vtkm::Float64; //vtkm::FloatDefault;
    using FloatArrayType = vtkm::cont::ArrayHandle<ValueType>;

public:
  // constructor
  SuperArcVolumetricComparator(const IdArrayType& weight,
                               const EdgePairArray& superArcList,
                               bool pairsAtLowEnd)
    : Weight(weight)
    , SuperArcList(superArcList)
    , PairsAtLowEnd(pairsAtLowEnd)
  {
      std::cout << "IdArrayType Comparator Constructor" << std::endl;
  }

  // constructor for floating point data
  SuperArcVolumetricComparator(const FloatArrayType& weight,
                               const EdgePairArray& superArcList,
                               bool pairsAtLowEnd)
    : WeightFloat(weight)
    , SuperArcList(superArcList)
    , PairsAtLowEnd(pairsAtLowEnd)
  {
    //      std::cout << "FloatArrayType Comparator Constructor" << std::endl;
  }

  VTKM_CONT SuperArcVolumetricComparatorImpl PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                                                 vtkm::cont::Token& token)
  {

      int intLen = this->Weight.GetNumberOfValues();
      int floatLen = this->WeightFloat.GetNumberOfValues();

      //      std::cout << "SuperArcVolumetricComparatorImpl" << std::endl;
      //      std::cout << "\tid    len: " << this->Weight.GetNumberOfValues() << std::endl;
      //      std::cout << "\tfloat len: " << this->WeightFloat.GetNumberOfValues() << std::endl;

      if (intLen > floatLen)
      {
          return SuperArcVolumetricComparatorImpl(
            this->Weight, this->SuperArcList, this->PairsAtLowEnd, device, token);
      }

      if (floatLen > intLen)
      {
          return SuperArcVolumetricComparatorImpl(
            this->WeightFloat, this->SuperArcList, this->PairsAtLowEnd, device, token);
      }
  }

private:
  IdArrayType Weight;
  FloatArrayType WeightFloat;
  EdgePairArray SuperArcList;
  bool PairsAtLowEnd;
}; // SuperArcVolumetricComparator

} // namespace process_contourtree_inc
} // namespace contourtree_augmented
} // namespace worklet
} // namespace vtkm

#endif
