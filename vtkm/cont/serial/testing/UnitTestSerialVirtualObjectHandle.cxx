//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2017 UT-Battelle, LLC.
//  Copyright 2017 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#include <vtkm/cont/RuntimeDeviceTracker.h>
#include <vtkm/cont/testing/TestingVirtualObjectHandle.h>

namespace
{

void TestVirtualObjectHandle()
{
  using DeviceAdapterList = vtkm::ListTagBase<vtkm::cont::DeviceAdapterTagSerial>;
  vtkm::cont::testing::TestingVirtualObjectHandle<DeviceAdapterList>::Run();
}

} // anonymous namespace

int UnitTestSerialVirtualObjectHandle(int argc, char* argv[])
{
  auto tracker = vtkm::cont::GetRuntimeDeviceTracker();
  tracker.ForceDevice(vtkm::cont::DeviceAdapterTagSerial{});
  return vtkm::cont::testing::Testing::Run(TestVirtualObjectHandle, argc, argv);
}
