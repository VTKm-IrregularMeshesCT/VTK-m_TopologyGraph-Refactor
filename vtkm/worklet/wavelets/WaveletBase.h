//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 Sandia Corporation.
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#ifndef vtk_m_worklet_wavelets_waveletbase_h
#define vtk_m_worklet_wavelets_waveletbase_h


#include <vtkm/worklet/wavelets/WaveletFilter.h>
#include <vtkm/worklet/wavelets/WaveletTransforms.h>

#include <vtkm/Math.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>

namespace vtkm {
namespace worklet {

namespace wavelets {

// Functionalities are similar to MatWaveBase in VAPoR.
class WaveletBase
{
public:

  // Constructor
  WaveletBase( WaveletName name ) : wname ( name ),
                                    filter( name )
  {
    if( wname == CDF9_7 || wname == BIOR4_4 || 
        wname == CDF5_3 || wname == BIOR2_2 )
    {
      this->wmode = SYMW;   // Default extension mode, see MatWaveBase.cpp
    }
    else if( wname == HAAR   || wname == BIOR1_1 ||
             wname == CDF8_4 || wname == BIOR3_3 )
    {
      this->wmode = SYMH;
    }
  }

  // Returns length of approximation coefficients from a decompostition pass.
  vtkm::Id GetApproxLength( vtkm::Id sigInLen )
  {
    vtkm::Id filterLen = this->filter.GetFilterLength();

    if (this->filter.isSymmetric()) 
    {
      if ( (this->wmode == SYMW && (filterLen % 2 != 0)) ||
           (this->wmode == SYMH && (filterLen % 2 == 0)) )  
      {
        if (sigInLen % 2 != 0)
          return((sigInLen+1) / 2);
        else 
          return((sigInLen) / 2);
      }
    }

    return ((sigInLen + filterLen - 1) / 2);
  }

  // Returns length of detail coefficients from a decompostition pass
  vtkm::Id GetDetailLength( vtkm::Id sigInLen )
  {
    vtkm::Id filterLen = this->filter.GetFilterLength();

    if (this->filter.isSymmetric()) 
    {
      if ( (this->wmode == SYMW && (filterLen % 2 != 0)) ||
           (this->wmode == SYMH && (filterLen % 2 == 0)) )  
      {
        if (sigInLen % 2 != 0)
          return((sigInLen-1) / 2);
        else 
          return((sigInLen) / 2);
      }
    }

    return static_cast<vtkm::Id>( vtkm::Floor(
           static_cast<vtkm::Float64>(sigInLen + filterLen - 1) / 2.0 ) );
  }

  // Returns length of coefficients generated in a decompostition pass
  vtkm::Id GetCoeffLength( vtkm::Id sigInLen )
  {
    return( GetApproxLength( sigInLen ) + GetDetailLength( sigInLen ) );
  }
  vtkm::Id GetCoeffLength2( vtkm::Id sigInX, vtkm::Id sigInY )
  {
    return( GetCoeffLength( sigInX) * GetCoeffLength( sigInY ) );
  }
  vtkm::Id GetCoeffLength3( vtkm::Id sigInX, vtkm::Id sigInY, vtkm::Id sigInZ)
  {
    return( GetCoeffLength( sigInX) * GetCoeffLength( sigInY ) * GetCoeffLength( sigInZ ) );
  }

  // Returns maximum wavelet decompostion level
  vtkm::Id GetWaveletMaxLevel( vtkm::Id sigInLen )
  {
    vtkm::Id filterLen = this->filter.GetFilterLength(); 
    vtkm::Id level;
    this->WaveLengthValidate( sigInLen, filterLen, level );
    return level;
  }

  // perform a device copy. The whole 1st array to a certain start location of the 2nd array
  template< typename ArrayType1, typename ArrayType2, typename DeviceTag >
  void DeviceCopyStartX( const ArrayType1   &srcArray, 
                               ArrayType2   &dstArray,
                               vtkm::Id     startIdx,
                               DeviceTag                )
  {
      typedef vtkm::worklet::wavelets::CopyWorklet CopyType;
      CopyType cp( startIdx );
      vtkm::worklet::DispatcherMapField< CopyType, DeviceTag > dispatcher( cp  );
      dispatcher.Invoke( srcArray, dstArray );
  }

  // Assign zero value to a certain location of an array
  template< typename ArrayType, typename DeviceTag >
  void DeviceAssignZero( ArrayType &array, vtkm::Id index, DeviceTag )
  {
    typedef vtkm::worklet::wavelets::AssignZeroWorklet ZeroWorklet;
    ZeroWorklet worklet( index );
    vtkm::worklet::DispatcherMapField< ZeroWorklet, DeviceTag > dispatcher( worklet );
    dispatcher.Invoke( array );
  }

  // Sort by the absolute value on device
  struct SortLessAbsFunctor
  { 
    template< typename T >
    VTKM_EXEC_EXPORT 
    bool operator()(const T& x, const T& y) const 
    { 
      return vtkm::Abs(x) < vtkm::Abs(y); 
    } 
  }; 
  template< typename ArrayType, typename DeviceTag >
  void DeviceSort( ArrayType &array, DeviceTag )
  {
    vtkm::cont::DeviceAdapterAlgorithm< DeviceTag >::Sort
          ( array, SortLessAbsFunctor() );
  }
  
  // Reduce to the sum of all values on device
  template< typename ArrayType, typename DeviceTag >
  typename ArrayType::ValueType DeviceSum( const ArrayType &array, DeviceTag )
  {
    return vtkm::cont::DeviceAdapterAlgorithm< DeviceTag >::Reduce
              ( array, 0.0 );
  }

  // Find the max and min of an array
  struct minFunctor
  {
    template< typename FieldType >
    VTKM_EXEC_EXPORT
    FieldType operator()(const FieldType &x, const FieldType &y) const {
      return Min(x, y);
    }
  };
  struct maxFunctor
  {
    template< typename FieldType >
    VTKM_EXEC_EXPORT
    FieldType operator()(const FieldType& x, const FieldType& y) const {
      return vtkm::Max(x, y);
    }
  };
  template< typename ArrayType, typename DeviceTag >
  typename ArrayType::ValueType DeviceMax( const ArrayType &array, DeviceTag )
  {
    typename ArrayType::ValueType initVal = array.GetPortalConstControl().Get(0);
    return vtkm::cont::DeviceAdapterAlgorithm< DeviceTag >::Reduce
              ( array, initVal, maxFunctor() );
  }
  template< typename ArrayType, typename DeviceTag >
  typename ArrayType::ValueType DeviceMin( const ArrayType &array, DeviceTag )
  {
    typename ArrayType::ValueType initVal = array.GetPortalConstControl().Get(0);
    return vtkm::cont::DeviceAdapterAlgorithm< DeviceTag >::Reduce
              ( array, initVal, minFunctor() );
  }

  // Max absolute value of an array
  struct maxAbsFunctor
  {
    template< typename FieldType >
    VTKM_EXEC_EXPORT
    FieldType operator()(const FieldType& x, const FieldType& y) const {
      return vtkm::Max( vtkm::Abs(x), vtkm::Abs(y) );
    }
  };
  template< typename ArrayType, typename DeviceTag >
  typename ArrayType::ValueType DeviceMaxAbs( const ArrayType &array, DeviceTag )
  {
    typename ArrayType::ValueType initVal = array.GetPortalConstControl().Get(0);
    return vtkm::cont::DeviceAdapterAlgorithm< DeviceTag >::Reduce
              ( array, initVal, maxAbsFunctor() );
  }

  // Calculate variance of an array
  template< typename ArrayType, typename DeviceTag >
  vtkm::Float64 DeviceCalculateVariance( ArrayType &array, DeviceTag )
  {
    vtkm::Float64 mean = static_cast<vtkm::Float64>(this->DeviceSum( array, DeviceTag() )) / 
                         static_cast<vtkm::Float64>(array.GetNumberOfValues());
    
    vtkm::cont::ArrayHandle< vtkm::Float64 > squaredDeviation;
    
    // Use a worklet
    typedef vtkm::worklet::wavelets::SquaredDeviation SDWorklet;
    SDWorklet sdw( mean );
    vtkm::worklet::DispatcherMapField< SDWorklet, DeviceTag > dispatcher( sdw  );
    dispatcher.Invoke( array, squaredDeviation );

    vtkm::Float64 sdMean = this->DeviceSum( squaredDeviation, DeviceTag() ) / 
                           static_cast<vtkm::Float64>( squaredDeviation.GetNumberOfValues() );

    return sdMean;
  }

  // Transpose a matrix in an array
  template< typename InputArrayType, typename OutputArrayType, typename DeviceTag >
  void DeviceTranspose( const InputArrayType    &inputArray,
                             OutputArrayType   &outputArray,
                                    vtkm::Id   inputX,
                                    vtkm::Id   inputY,
                                   DeviceTag            )
  {
    // use a worklet
    typedef vtkm::worklet::wavelets::TransposeWorklet TransposeType;
    TransposeType tw ( inputX, inputY );
    vtkm::worklet::DispatcherMapField< TransposeType, DeviceTag > dispatcher( tw );
    dispatcher.Invoke( inputArray, outputArray );
  }

  // Copy a small rectangle to a big rectangle
  template< typename SmallArrayType, typename BigArrayType, typename DeviceTag>
  void DeviceRectangleCopyTo( const SmallArrayType    &smallRect,
                                    vtkm::Id          smallX,
                                    vtkm::Id          smallY,
                                    BigArrayType      &bigRect,
                                    vtkm::Id          bigX,
                                    vtkm::Id          bigY,
                                    vtkm::Id          startX,
                                    vtkm::Id          startY,
                                    DeviceTag                  )
  {
    typedef vtkm::worklet::wavelets::RectangleCopyTo  CopyToWorklet;
    CopyToWorklet cp( smallX, smallY, bigX, bigY, startX, startY );
    vtkm::worklet::DispatcherMapField< CopyToWorklet, DeviceTag > dispatcher( cp  );
    dispatcher.Invoke(smallRect, bigRect);
  }

  // Fill a small rectangle from a portion of a big rectangle
  template< typename SmallArrayType, typename BigArrayType, typename DeviceTag >
  void DeviceRectangleCopyFrom(       SmallArrayType    &smallRect,
                                      vtkm::Id          smallX,
                                      vtkm::Id          smallY,
                                const BigArrayType      &bigRect,
                                      vtkm::Id          bigX,
                                      vtkm::Id          bigY,
                                      vtkm::Id          startX,
                                      vtkm::Id          startY,
                                      DeviceTag                      )
  {
    smallRect.PrepareForOutput( smallX*smallY, DeviceTag() );
    typedef vtkm::worklet::wavelets::RectangleCopyFrom  CopyFromWorklet;
    CopyFromWorklet cpFrom( smallX, smallY, bigX, bigY, startX, startY );
    vtkm::worklet::DispatcherMapField< CopyFromWorklet, DeviceTag > dispatcherFrom( cpFrom );
    dispatcherFrom.Invoke( smallRect, bigRect );
  }



protected:
  WaveletName                               wname;
  DWTMode                                   wmode;
  WaveletFilter                             filter;

  void WaveLengthValidate( vtkm::Id sigInLen, vtkm::Id filterLength, vtkm::Id &level)
  {
    if( sigInLen < filterLength )
      level = 0;
    else
      level = static_cast<vtkm::Id>( vtkm::Floor( 
                  vtkm::Log2( static_cast<vtkm::Float64>(sigInLen) / 
                              static_cast<vtkm::Float64>(filterLength) ) + 1.0 ) );
  }
};    // class WaveletBase.


}     // namespace wavelets

}     // namespace worklet
}     // namespace vtkm

#endif 
