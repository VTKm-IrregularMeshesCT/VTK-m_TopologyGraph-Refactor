//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_cont_Invoker_h
#define vtk_m_cont_Invoker_h

#include <vtkm/worklet/internal/MaskBase.h>
#include <vtkm/worklet/internal/ScatterBase.h>

#include <vtkm/cont/TryExecute.h>

#define PACT_DEBUG 0

namespace vtkm
{
namespace cont
{

namespace detail
{
template <typename T>
using scatter_or_mask = std::integral_constant<bool,
                                               vtkm::worklet::internal::is_mask<T>::value ||
                                                 vtkm::worklet::internal::is_scatter<T>::value>;
}

/// \brief Allows launching any worklet without a dispatcher.
///
/// \c Invoker is a generalized \c Dispatcher that is able to automatically
/// determine how to properly launch/invoke any worklet that is passed to it.
/// When an \c Invoker is constructed it is provided the desired device adapter
/// that all worklets invoked by it should be launched on.
///
/// \c Invoker is designed to not only reduce the verbosity of constructing
/// multiple dispatchers inside a block of logic, but also makes it easier to
/// make sure all worklets execute on the same device.
struct Invoker
{

  /// Constructs an Invoker that will try to launch worklets on any device
  /// that is enabled.
  ///
  explicit Invoker()
    : DeviceId(vtkm::cont::DeviceAdapterTagAny{})
  {
        /// DEBUG PRINT std::cout << "ExplicitInvoker1\n";
  }

  /// Constructs an Invoker that will try to launch worklets only on the
  /// provided device adapter.
  ///
  explicit Invoker(vtkm::cont::DeviceAdapterId device)
    : DeviceId(device)
  {
        /// DEBUG PRINT std::cout << "ExplicitInvoker2\n";
  }

  /// Launch the worklet that is provided as the first parameter.
  /// Optional second parameter is either the scatter or mask type associated with the worklet.
  /// Any additional parameters are the ControlSignature arguments for the worklet.
  ///
  template <typename Worklet,
            typename T,
            typename... Args,
            typename std::enable_if<detail::scatter_or_mask<T>::value, int>::type* = nullptr>
  inline void operator()(Worklet&& worklet, T&& scatterOrMask, Args&&... args) const
  {
    using WorkletType = vtkm::internal::remove_cvref<Worklet>;
    using DispatcherType = typename WorkletType::template Dispatcher<WorkletType>;
    // Invoke() cont::Invoker Invoker.h Invoker_h
    //    std::cout << "Invoker1\n";
    DispatcherType dispatcher(worklet, scatterOrMask);
    dispatcher.SetDevice(this->DeviceId);
    dispatcher.Invoke(std::forward<Args>(args)...);
  }

  /// Launch the worklet that is provided as the first parameter.
  /// Optional second parameter is either the scatter or mask type associated with the worklet.
  /// Optional third parameter is either the scatter or mask type associated with the worklet.
  /// Any additional parameters are the ControlSignature arguments for the worklet.
  ///
  template <
    typename Worklet,
    typename T,
    typename U,
    typename... Args,
    typename std::enable_if<detail::scatter_or_mask<T>::value && detail::scatter_or_mask<U>::value,
                            int>::type* = nullptr>
  inline void operator()(Worklet&& worklet,
                         T&& scatterOrMaskA,
                         U&& scatterOrMaskB,
                         Args&&... args) const
  {
    using WorkletType = vtkm::internal::remove_cvref<Worklet>;
    using DispatcherType = typename WorkletType::template Dispatcher<WorkletType>;
    std::cout << "Invoker2\n";
    DispatcherType dispatcher(worklet, scatterOrMaskA, scatterOrMaskB);
    dispatcher.SetDevice(this->DeviceId);
    dispatcher.Invoke(std::forward<Args>(args)...);
  }

  /// Launch the worklet that is provided as the first parameter.
  /// Optional second parameter is either the scatter or mask type associated with the worklet.
  /// Any additional parameters are the ControlSignature arguments for the worklet.
  ///
  template <typename Worklet,
            typename T,
            typename... Args,
            typename std::enable_if<!detail::scatter_or_mask<T>::value, int>::type* = nullptr>
  inline void operator()(Worklet&& worklet, T&& t, Args&&... args) const
  {
    #if PACT_DEBUG
        printf("INVOKER() on: %s\n", __PRETTY_FUNCTION__);
        std::cout << "Invoker3\n";
    #endif
    using WorkletType = vtkm::internal::remove_cvref<Worklet>;
    #if PACT_DEBUG
    std::cout << "Invoker3.1\n";
    #endif
    using DispatcherType = typename WorkletType::template Dispatcher<WorkletType>;
    #if PACT_DEBUG
    std::cout << "Invoker3.2\n";
    #endif
    DispatcherType dispatcher(worklet);
    #if PACT_DEBUG
    std::cout << "Invoker3.3\n";
    #endif
    dispatcher.SetDevice(this->DeviceId);
    #if PACT_DEBUG
    std::cout << "Invoker3.4\n";
    #endif
    dispatcher.Invoke(std::forward<T>(t), std::forward<Args>(args)...);

    #if PACT_DEBUG
    std::cout << "Invoker3.5 - end\n";
    #endif
  }

  /// Get the device adapter that this Invoker is bound too
  ///
  vtkm::cont::DeviceAdapterId GetDevice() const { return DeviceId; }

private:
  vtkm::cont::DeviceAdapterId DeviceId;
};
}
}

#endif
