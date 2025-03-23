//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/io/VTKUnstructuredGridReader.h>

#include <vtkm/io/internal/VTKDataSetCells.h>

#include <vtkm/cont/ConvertNumComponentsToOffsets.h>

namespace vtkm
{
namespace io
{

VTKUnstructuredGridReader::VTKUnstructuredGridReader(const char* fileName)
  : VTKDataSetReaderBase(fileName)
{
}

VTKUnstructuredGridReader::VTKUnstructuredGridReader(const std::string& fileName)
  : VTKDataSetReaderBase(fileName)
{
}

void VTKUnstructuredGridReader::Read()
{
  std::cout << "VTKUnstructuredGridReader::Read() called" << std::endl;
  if (this->DataFile->Structure != vtkm::io::internal::DATASET_UNSTRUCTURED_GRID)
  {
    throw vtkm::io::ErrorIO("Incorrect DataSet type");
  }

  //We need to be able to handle VisIt files which dump Field data
  //at the top of a VTK file
  std::string tag;
  this->DataFile->Stream >> tag;
  if (tag == "FIELD")
  {
    std::cout << "VTKUnstructuredGridReader::Read() --> Parsing 'FIELD'" << std::endl;
    this->ReadGlobalFields();
    this->DataFile->Stream >> tag;
  }

  // Read the points
  internal::parseAssert(tag == "POINTS");
  std::cout << "VTKUnstructuredGridReader::Read() --> parseAssert 'POINTS'" << std::endl;
  std::cout << "VTKUnstructuredGridReader::Read() --> Call to VTKDataSetReaderBase::ReadPoints()" << std::endl;
  this->ReadPoints();

  vtkm::Id numPoints = this->DataSet.GetNumberOfPoints();
  std::cout << "VTKUnstructuredGridReader::Read() --> " << numPoints << " POINTS just read" << std::endl;

  // Read the cellset
  vtkm::cont::ArrayHandle<vtkm::Id> connectivity;
  vtkm::cont::ArrayHandle<vtkm::IdComponent> numIndices;
  vtkm::cont::ArrayHandle<vtkm::UInt8> shapes;

  std::cout << "VTKUnstructuredGridReader::Read() --> Reading the CellSet" << std::endl;

  this->DataFile->Stream >> tag;
  std::cout << "VTKUnstructuredGridReader::Read() --> this->DataFile->Stream >> tag=" << tag << std::endl;
  internal::parseAssert(tag == "CELLS");
  std::cout << "VTKUnstructuredGridReader::Read() --> parseAssert 'CELLS'" << std::endl;
  std::cout << "VTKUnstructuredGridReader::Read() --> Call to VTKDataSetReaderBase::ReadCells(connectivity, numIndices)" << std::endl;
  this->ReadCells(connectivity, numIndices);
  // connectivity stores N vertices per CELL (where N is the number of vertices in the cell)
  std::cout << "VTKUnstructuredGridReader::Read() --> " << connectivity.GetNumberOfValues()/4 << " 'CELLS' just read" << std::endl;

  std::cout << "VTKUnstructuredGridReader::Read() --> Reading the 'CELL_TYPES'" << std::endl;
  this->ReadShapes(shapes);
  std::cout << "VTKUnstructuredGridReader::Read() --> " << shapes.GetNumberOfValues() << "'CELL_TYPES' just read" << std::endl;

  vtkm::cont::ArrayHandle<vtkm::Id> permutation;
  std::cout << "VTKUnstructuredGridReader::Read() --> Call to vtkm::io::internal::FixupCellSet(...)" << std::endl;
  vtkm::io::internal::FixupCellSet(connectivity, numIndices, shapes, permutation);
  std::cout << "VTKUnstructuredGridReader::Read() --> Call to SetCellsPermutation(permutation);" << std::endl;
  this->SetCellsPermutation(permutation);

  if (vtkm::io::internal::IsSingleShape(shapes))
  {
    std::cout << "VTKUnstructuredGridReader::Read() --> CellSet is Single Type (all 'CELL_TYPES' are the same)" << std::endl;
    vtkm::cont::CellSetSingleType<> cellSet;
    cellSet.Fill(
      numPoints, shapes.ReadPortal().Get(0), numIndices.ReadPortal().Get(0), connectivity);
    std::cout << "VTKUnstructuredGridReader::Read() --> Call to this->DataSet.SetCellSet(cellSet)" << std::endl;
    this->DataSet.SetCellSet(cellSet);
  }
  else
  {
    auto offsets = vtkm::cont::ConvertNumComponentsToOffsets(numIndices);
    vtkm::cont::CellSetExplicit<> cellSet;
    cellSet.Fill(numPoints, shapes, connectivity, offsets);
    this->DataSet.SetCellSet(cellSet);
  }

  // Read points and cell attributes
  std::cout << "VTKUnstructuredGridReader::Read() --> Read points and cell attributes" << std::endl;
  this->ReadAttributes();
  std::cout << "VTKUnstructuredGridReader::Read() -->  -> Returning" << std::endl;;
}
}
}
