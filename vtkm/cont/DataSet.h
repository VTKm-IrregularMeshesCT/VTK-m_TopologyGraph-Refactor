#ifndef vtk_m_cont_DataSet_h
#define vtk_m_cont_DataSet_h

#include <vtkm/CellType.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/DynamicArrayHandle.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>
#include <vtkm/cont/ExplicitConnectivity.h>
#include <vtkm/cont/RegularConnectivity.h>

namespace vtkm {
namespace cont {

class DataSet
{
public:
  DataSet() {}

  template <typename T>
  void AddFieldViaCopy(T *ptr, int nvals)
  {
    vtkm::cont::ArrayHandle<T> tmp = vtkm::cont::make_ArrayHandle(ptr, nvals);
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> array;
    vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::
      Copy(tmp, array);
    Fields.resize(Fields.size()+1);
    Fields[Fields.size()-1].SetData(array);
    /*
    Fields.resize(Fields.size()+1);
    Fields[Fields.size()-1].CopyIntoData(tmp);
    */
  }
  vtkm::cont::Field &GetField(int index)
  {
    return Fields[index];
  }

  vtkm::Id x_idx, y_idx, z_idx;

  ExplicitConnectivity conn;
  RegularConnectivity3D reg;
    //TODO: Logical structure: vtkm::Extents?  Use EAVL logicalStructure?

  //traditional data-model
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault,3> > Points;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault,1> > Field;

private:
  std::vector<vtkm::cont::Field> Fields;


};

}
} // namespace vtkm::cont


#endif //vtk_m_cont_DataSet_h
