//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
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

#include <vtkm/Types.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DeviceAdapterTag.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/RuntimeDeviceTracker.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/io/BOVDataSetReader.h>

#include <vtkm/filter/MapFieldPermutation.h>
#include <vtkm/filter/scalar_topology/ContourTreeUniformAugmented.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/PrintVectors.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ProcessContourTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/Branch.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKPolyDataReader.h>
//#include <vtkNew.h>

// clang-format off
VTKM_THIRDPARTY_PRE_INCLUDE
#include <vtkm/thirdparty/diy/Configure.h>
#include <vtkm/thirdparty/diy/diy.h>
VTKM_THIRDPARTY_POST_INCLUDE
// clang-format on

//#define WITH_MPI

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include <chrono>
#include <thread>

//#include <valgrind/callgrind.h>

#define DEBUG_PRINT_PACTBD 0
#define SLEEP_ON 0
#define WRITE_FILES 1

//using vtkm::FloatDefault = vtkm::Float64;

//    using ValueType = vtkm::Float32;
using ValueType = vtkm::Float64; //vtkm::FloatDefault;
using FloatArrayType = vtkm::cont::ArrayHandle<ValueType>;
using BranchType = vtkm::worklet::contourtree_augmented::process_contourtree_inc::Branch<ValueType>;

namespace ctaug_ns = vtkm::worklet::contourtree_augmented;
using Coefficients = vtkm::worklet::contourtree_augmented::Coefficients;

// Simple helper class for parsing the command line options
class ParseCL
{
public:
  ParseCL() {}

  void parse(int& argc, char** argv)
  {
    mCLOptions.resize(static_cast<std::size_t>(argc));
    for (std::size_t i = 1; i < static_cast<std::size_t>(argc); ++i)
    {
      this->mCLOptions[i] = std::string(argv[i]);
    }
  }

  vtkm::Id findOption(const std::string& option) const
  {
    auto it =
      std::find_if(this->mCLOptions.begin(),
                   this->mCLOptions.end(),
                   [option](const std::string& val) -> bool { return val.find(option) == 0; });
    if (it == this->mCLOptions.end())
    {
      return -1;
    }
    else
    {
      return static_cast<vtkm::Id>(it - this->mCLOptions.begin());
    }
  }

  bool hasOption(const std::string& option) const { return this->findOption(option) >= 0; }

  std::string getOption(const std::string& option) const
  {
    std::size_t index = static_cast<std::size_t>(this->findOption(option));
    std::string val = this->mCLOptions[index];
    auto valPos = val.find("=");
    if (valPos)
    {
      return val.substr(valPos + 1);
    }
    return std::string("");
  }

  const std::vector<std::string>& getOptions() const { return this->mCLOptions; }

private:
  std::vector<std::string> mCLOptions;
};

inline vtkm::Id3 ComputeNumberOfBlocksPerAxis(vtkm::Id3 globalSize, vtkm::Id numberOfBlocks)
{
  vtkm::Id currNumberOfBlocks = numberOfBlocks;
  vtkm::Id3 blocksPerAxis{ 1, 1, 1 };
  while (currNumberOfBlocks > 1)
  {
    vtkm::IdComponent splitAxis = 0;
    for (vtkm::IdComponent d = 1; d < 3; ++d)
    {
      if (globalSize[d] > globalSize[splitAxis])
      {
        splitAxis = d;
      }
    }
    if (currNumberOfBlocks % 2 == 0)
    {
      blocksPerAxis[splitAxis] *= 2;
      globalSize[splitAxis] /= 2;
      currNumberOfBlocks /= 2;
    }
    else
    {
      blocksPerAxis[splitAxis] *= currNumberOfBlocks;
      break;
    }
  }
  return blocksPerAxis;
}

inline std::tuple<vtkm::Id3, vtkm::Id3, vtkm::Id3> ComputeBlockExtents(vtkm::Id3 globalSize,
                                                                       vtkm::Id3 blocksPerAxis,
                                                                       vtkm::Id blockNo)
{
  // DEBUG: std::cout << "ComputeBlockExtents("<<globalSize <<", " << blocksPerAxis << ", " << blockNo << ")" << std::endl;
  // DEBUG: std::cout << "Block " << blockNo;

  vtkm::Id3 blockIndex, blockOrigin, blockSize;
  for (vtkm::IdComponent d = 0; d < 3; ++d)
  {
    blockIndex[d] = blockNo % blocksPerAxis[d];
    blockNo /= blocksPerAxis[d];

    float dx = float(globalSize[d] - 1) / float(blocksPerAxis[d]);
    blockOrigin[d] = vtkm::Id(blockIndex[d] * dx);
    vtkm::Id maxIdx =
      blockIndex[d] < blocksPerAxis[d] - 1 ? vtkm::Id((blockIndex[d] + 1) * dx) : globalSize[d] - 1;
    blockSize[d] = maxIdx - blockOrigin[d] + 1;
    // DEBUG: std::cout << " " << blockIndex[d] <<  dx << " " << blockOrigin[d] << " " << maxIdx << " " << blockSize[d] << "; ";
  }
  // DEBUG: std::cout << " -> " << blockIndex << " "  << blockOrigin << " " << blockSize << std::endl;
  return std::make_tuple(blockIndex, blockOrigin, blockSize);
}

inline vtkm::cont::DataSet CreateSubDataSet(const vtkm::cont::DataSet& ds,
                                            vtkm::Id3 blockOrigin,
                                            vtkm::Id3 blockSize,
                                            const std::string& fieldName)
{
  vtkm::Id3 globalSize;
  ds.GetCellSet().CastAndCallForTypes<VTKM_DEFAULT_CELL_SET_LIST_STRUCTURED>(
    vtkm::worklet::contourtree_augmented::GetPointDimensions(), globalSize);
  const vtkm::Id nOutValues = blockSize[0] * blockSize[1] * blockSize[2];

  const auto inDataArrayHandle = ds.GetPointField(fieldName).GetData();

  vtkm::cont::ArrayHandle<vtkm::Id> copyIdsArray;
  copyIdsArray.Allocate(nOutValues);
  auto copyIdsPortal = copyIdsArray.WritePortal();

  vtkm::Id3 outArrIdx;
  for (outArrIdx[2] = 0; outArrIdx[2] < blockSize[2]; ++outArrIdx[2])
    for (outArrIdx[1] = 0; outArrIdx[1] < blockSize[1]; ++outArrIdx[1])
      for (outArrIdx[0] = 0; outArrIdx[0] < blockSize[0]; ++outArrIdx[0])
      {
        vtkm::Id3 inArrIdx = outArrIdx + blockOrigin;
        vtkm::Id inIdx = (inArrIdx[2] * globalSize[1] + inArrIdx[1]) * globalSize[0] + inArrIdx[0];
        vtkm::Id outIdx =
          (outArrIdx[2] * blockSize[1] + outArrIdx[1]) * blockSize[0] + outArrIdx[0];
        VTKM_ASSERT(inIdx >= 0 && inIdx < inDataArrayHandle.GetNumberOfValues());
        VTKM_ASSERT(outIdx >= 0 && outIdx < nOutValues);
        copyIdsPortal.Set(outIdx, inIdx);
      }
  // DEBUG: std::cout << copyIdsPortal.GetNumberOfValues() << std::endl;

  vtkm::cont::Field permutedField;
  bool success =
    vtkm::filter::MapFieldPermutation(ds.GetPointField(fieldName), copyIdsArray, permutedField);
  if (!success)
    throw vtkm::cont::ErrorBadType("Field copy failed (probably due to invalid type)");


  vtkm::cont::DataSetBuilderUniform dsb;
  if (globalSize[2] <= 1) // 2D Data Set
  {
    vtkm::Id2 dimensions{ blockSize[0], blockSize[1] };
    vtkm::cont::DataSet dataSet = dsb.Create(dimensions);
    vtkm::cont::CellSetStructured<2> cellSet;
    cellSet.SetPointDimensions(dimensions);
    cellSet.SetGlobalPointDimensions(vtkm::Id2{ globalSize[0], globalSize[1] });
    cellSet.SetGlobalPointIndexStart(vtkm::Id2{ blockOrigin[0], blockOrigin[1] });
    dataSet.SetCellSet(cellSet);
    dataSet.AddField(permutedField);
    return dataSet;
  }
  else
  {
    vtkm::cont::DataSet dataSet = dsb.Create(blockSize);
    vtkm::cont::CellSetStructured<3> cellSet;
    cellSet.SetPointDimensions(blockSize);
    cellSet.SetGlobalPointDimensions(globalSize);
    cellSet.SetGlobalPointIndexStart(blockOrigin);
    dataSet.SetCellSet(cellSet);
    dataSet.AddField(permutedField);
    return dataSet;
  }
}

void static printMemoryUsage(const std::string& message)
{
    // Red text formatting for highlighting some console output:
    const std::string ORANGE = "\033[38;2;255;165;0m";  // Start red text
    const std::string LIGHT_BLUE = "\033[38;5;117m";  // Light blue in 256-color
    const std::string RESET = "\033[0m"; // End red text

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

//    std::cout << ORANGE << message << LIGHT_BLUE << " - Memory usage: " << usage.ru_maxrss << " KB" << RESET << std::endl;
    std::cout << LIGHT_BLUE << message << " - Memory usage: " << usage.ru_maxrss << " KB" << RESET << std::endl;
}

// Compute and render an isosurface for a uniform grid example
int main(int argc, char* argv[])
{
    // Red text formatting for highlighting some console output:
    const std::string RED = "\033[31m";  // Start red text
    const std::string ORANGE = "\033[38;2;255;165;0m";  // Start red text
    const std::string YELLOW = "\033[38;2;240;240;13m";  // Warm, readable yellow
    const std::string RESET = "\033[0m"; // End red text
    /////////////////////////////////////////////////
    // START MAIN-0 Initialise VTK-m (totalTime) //
    /////////////////////////////////////////////////

    std::cout << ORANGE << std::endl;
    std::cout << "///////////////////////////////////////////////" << std::endl;
    std::cout << "// START MAIN-0 Initialise VTK-m (totalTime) //" << std::endl;
    std::cout << "///////////////////////////////////////////////" << std::endl;
    std::cout << RESET << std::endl;

#if WRITE_FILES
    int file_io_counter = 0;
#endif

#ifdef WITH_MPI
  // Setup the MPI environment.
  MPI_Init(&argc, &argv);
  auto comm = MPI_COMM_WORLD;

  // Tell VTK-m which communicator it should use.
  vtkm::cont::EnvironmentTracker::SetCommunicator(vtkmdiy::mpi::communicator());

  // get the rank and size
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  int numBlocks = size;
#endif

  // initialize vtkm-m (e.g., logging via -v and device via the -d option)
  vtkm::cont::InitializeOptions vtkm_initialize_options = vtkm::cont::InitializeOptions::RequireDevice;
  vtkm::cont::InitializeResult vtkm_config = vtkm::cont::Initialize(argc, argv, vtkm_initialize_options);
  auto device = vtkm_config.Device; // gets the device from --vtkm-device command line argument (mandatory)

#ifdef WITH_MPI
  VTKM_LOG_IF_S(vtkm::cont::LogLevel::Info, rank == 0, "Running with MPI. #ranks=" << size);
#else
  VTKM_LOG_S(vtkm::cont::LogLevel::Info, "Single node run");
  int rank = 0;
#endif

  // Setup timing for the whole program (for performance runs)
  ///TIMING///////////////////////////////////////////////////////
  vtkm::cont::Timer totalTime;
  totalTime.Start();
  ////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////
  // START MAIN-1 Parse the command line options (startUpTimeDisplay) //
  //////////////////////////////////////////////////////////////////////

  std::cout << ORANGE << std::endl;
  std::cout << "////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "// MAIN-1 Parse the command line options (startUpTimeDisplay) //" << std::endl;
  std::cout << "////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << RESET << std::endl;

  // Time how long parsing command line arguments takes:
  ///TIMING///////////////////////////////////////////////////////
  vtkm::Float64 currTime = 0; // TIMING: start time here from 0
  vtkm::Float64 prevTime = 0; // TIMING: no other previous yet
  // first timed category will be the Start-up (parsing arguments)
  vtkm::Float64 startUpTimeDisplay;// startUpTimeDisplay timer begins
  ////////////////////////////////////////////////////////////////


  // Create the command-line parser object and get the input filename
  ParseCL parser;
  parser.parse(argc, argv);
  std::string filename = parser.getOptions().back();
  // some flags will be true by default and the command-line will be just to override them
  unsigned int augmentComputeRegularStructure = 1; // 1=fully augmented, augment by default, ...
                                                   // ... required for the branch decomposition
  bool useMarchingCubes = false;                   // use Freudenthal by default (on regular data) ...
                                                   // ... TODO: deprecate for irregular grids
  bool computeBranchDecomposition = true;          // Requires --augmentTree (Default=True)
  bool printContourTree = false;                   // Print the contour tree. (Default=False)
  bool printDebug       = true;                    // Print debug information (Default=True)

  // start parsing arguments
  if (parser.hasOption("--augmentTree"))
  { //      --augmentTree     1 = compute the fully augmented contour tree (Default)
    //                        2 = compute the boundary augmented contour tree
    //                        0 = no augmentation. NOTE: When using MPI, local ranks use
    augmentComputeRegularStructure = static_cast<unsigned int>(std::stoi(parser.getOption("--augmentTree")));
  }
  if (parser.hasOption("--mc"))
  {// TOGGLE ONLY
   // Use marching cubes interpolation for contour tree calculation. (Default=False)
      useMarchingCubes = true;
  }
  if (parser.hasOption("--printCT"))
  {// TOGGLE ONLY
   // --printCT         Print the contour tree. (Default=False)
      printContourTree = true;
  }
  if (parser.hasOption("--branchDecomp"))
  {// BOOLEAN 1 or 0 only (Default=True)
   // Compute the volume branch decomposition for the contour tree.
   // (Requires --augmentTree, which is Default=True)
    computeBranchDecomposition = std::stoi(parser.getOption("--branchDecomp"));
  }
  // We need the fully augmented tree to compute the branch decomposition, check that here:
  if (computeBranchDecomposition && (augmentComputeRegularStructure != 1))
  {
    VTKM_LOG_S(vtkm::cont::LogLevel::Warn,
               "Regular structure is required for branch decomposition."
               " Disabling branch decomposition");
    computeBranchDecomposition = false;
  }

  // Iso value selection parameters
  // Approach to be used to select contours based on the tree
  vtkm::Id contourType = 0;
  // Error away from critical point
  ValueType eps = 0.00001f;
  // Number of iso levels to be selected. By default we disable the isovalue selection.
  vtkm::Id numLevels = 0;
  // Number of components the tree should be simplified to
  vtkm::Id numComp = numLevels + 1;
  // Method to be used to compute the relevant iso values
  vtkm::Id contourSelectMethod = 0;
  bool usePersistenceSorter = true;
  if (parser.hasOption("--levels"))
    numLevels = std::stoi(parser.getOption("--levels"));
  if (parser.hasOption("--type"))
    contourType = std::stoi(parser.getOption("--type"));
  if (parser.hasOption("--eps"))
    eps = std::stof(parser.getOption("--eps"));
  if (parser.hasOption("--method"))
    contourSelectMethod = std::stoi(parser.getOption("--method"));
  if (parser.hasOption("--comp"))
    numComp = std::stoi(parser.getOption("--comp"));
  if (contourSelectMethod == 0)
    numComp = numLevels + 1;
  if (parser.hasOption("--useVolumeSorter"))
    usePersistenceSorter = false;
  if ((numLevels > 0) && (!computeBranchDecomposition))
  {
    VTKM_LOG_S(vtkm::cont::LogLevel::Warn,
               "Iso level selection only available when branch decomposition is enabled. "
               "Disabling iso value selection");
    numLevels = 0;
  }

  if (rank == 0 && (argc < 2 || parser.hasOption("--help") || parser.hasOption("-h")))
  {
    std::cout << "ContourTreeAugmented <options> <fileName>" << std::endl;
    std::cout << std::endl;
    std::cout << "<fileName>       Name of the input data file." << std::endl;
    std::cout << "The file is expected to be ASCII with either: " << std::endl;
    std::cout << "  - xdim ydim integers for 2D or" << std::endl;
    std::cout << "  - xdim ydim zdim integers for 3D" << std::endl;
    std::cout << "followed by vector data last dimension varying fastest" << std::endl;
    std::cout << std::endl;
    std::cout << "----------------------------- VTKM Options -----------------------------"
              << std::endl;
    std::cout << vtkm_config.Usage << std::endl;
    std::cout << std::endl;
    std::cout << "------------------------- Contour Tree Options -------------------------"
              << std::endl;
    std::cout << "Options: (Bool options are give via int, i.e. =0 for False and =1 for True)"
              << std::endl;
    std::cout << "--mc              Use marching cubes interpolation for contour tree calculation. "
                 "(Default=False)"
              << std::endl;
    std::cout << "--augmentTree     1 = compute the fully augmented contour tree (Default)"
              << std::endl;
    std::cout << "                  2 = compute the boundary augmented contour tree " << std::endl;
    std::cout << "                  0 = no augmentation. NOTE: When using MPI, local ranks use"
              << std::endl;
    std::cout << "                      boundary augmentation to support parallel merge of blocks"
              << std::endl;
    std::cout << "--branchDecomp    Compute the volume branch decomposition for the contour tree. "
                 "Requires --augmentTree (Default=True)"
              << std::endl;
    std::cout << "--printCT         Print the contour tree. (Default=False)" << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------- Isovalue Selection Options ----------------------"
              << std::endl;
    std::cout << "Isovalue selection options: (require --branchDecomp=1 and augmentTree=1)"
              << std::endl;
    std::cout << "--levels=<int>  Number of iso-contour levels to be used (default=0, i.e., "
                 "disable isovalue computation)"
              << std::endl;
    std::cout << "--comp=<int>    Number of components the contour tree should be simplified to. "
                 "Only used if method==1. (default=0)"
              << std::endl;
    std::cout
      << "--eps=<float>   Floating point offset awary from the critical point. (default=0.00001)"
      << std::endl;
    std::cout << "--type=<int>    Approach to be used for selection of iso-values. 0=saddle+-eps; "
                 "1=mid point between saddle and extremum, 2=extremum+-eps. (default=0)"
              << std::endl;
    std::cout << "--method=<int>  Method used for selecting relevant iso-values. (default=0)"
              << std::endl;
    std::cout << std::endl;

#ifdef WITH_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
  }

  if (rank == 0)
  {
    std::stringstream logmessage;
    logmessage << "    ------------ Settings -----------" << std::endl
               << "    filename=" << filename << std::endl
               << "    device=" << device.GetName() << std::endl
               << "    mc=" << useMarchingCubes << std::endl
               << "    augmentTree=" << augmentComputeRegularStructure << std::endl
               << "    branchDecomp=" << computeBranchDecomposition << std::endl
               << "    printCT=" << printContourTree << std::endl
               <<
#ifdef WITH_MPI
      "    nblocks=" << numBlocks << std::endl
               <<
#endif
      "    computeIsovalues=" << (numLevels > 0);
    VTKM_LOG_S(vtkm::cont::LogLevel::Info, std::endl << logmessage.str());
    VTKM_LOG_IF_S(vtkm::cont::LogLevel::Info,
                  numLevels > 0,
                  std::endl
                    << "    ------------ Settings Isolevel Selection -----------" << std::endl
                    << "    levels=" << numLevels << std::endl
                    << "    eps=" << eps << std::endl
                    << "    comp" << numComp << std::endl
                    << "    type=" << contourType << std::endl
                    << "    method=" << contourSelectMethod << std::endl
                    << "    mc=" << useMarchingCubes << std::endl
                    << "    use" << (usePersistenceSorter ? "PersistenceSorter" : "VolumeSorter"));
  }



// Redirect stdout to file if we are using MPI with Debugging
#ifdef WITH_MPI
#ifdef DEBUG_PRINT
  // From https://www.unix.com/302983597-post2.html
  char cstr_filename[32];
  snprintf(cstr_filename, sizeof(cstr_filename), "cout_%d.log", rank);
  int out = open(cstr_filename, O_RDWR | O_CREAT | O_TRUNC | O_APPEND, 0600);
  if (-1 == out)
  {
    perror("opening cout.log");
    return 255;
  }

  snprintf(cstr_filename, sizeof(cstr_filename), "cerr_%d.log", rank);
  int err = open(cstr_filename, O_RDWR | O_CREAT | O_TRUNC | O_APPEND, 0600);
  if (-1 == err)
  {
    perror("opening cerr.log");
    return 255;
  }

  int save_out = dup(fileno(stdout));
  int save_err = dup(fileno(stderr));

  if (-1 == dup2(out, fileno(stdout)))
  {
    perror("cannot redirect stdout");
    return 255;
  }
  if (-1 == dup2(err, fileno(stderr)))
  {
    perror("cannot redirect stderr");
    return 255;
  }
#endif
#endif

  //////////////////////////////////////////////////////////////////////
  // END MAIN-1 Parse the command line options (startUpTimeDisplay) //
  //////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // START MAIN-2 Read the input data (buildDatasetTimeDisplay) //
  ////////////////////////////////////////////////////////////

  std::cout << ORANGE << std::endl;
  std::cout << "//////////////////////////////////////////////////////////" << std::endl;
  std::cout << "// MAIN-2 Read the input data (buildDatasetTimeDisplay) //" << std::endl;
  std::cout << "//////////////////////////////////////////////////////////" << std::endl;
  std::cout << RESET << std::endl;

  // Finish timing the 'start-up time', begin timing 'input data reading-in'
  ///TIMING/////////////////////////////////////////////////////////////////
  currTime = totalTime.GetElapsedTime();
  // TIMING: Start-up (finish) - timing saved to 'startUpTimeDisplay'
  startUpTimeDisplay = currTime - prevTime;
  // TIMING: Data Read (begin) - timing saved to 'dataReadTimeDisplay'
  //                             ... and 'buildDatasetTimeDisplay'
  // ('read in' from file and 'build'/parse to VTK-m format)
  vtkm::Float64 dataReadTimeDisplay = 0;
  vtkm::Float64 buildDatasetTimeDisplay = 0;
  // (for .vtk files the 'dataReadTimeDisplay' time will be same as 'buildDatasetTimeDisplay' ...
  // ... because the dataset gets read-in as an internal vtk data structure already, ...
  // ... unlike ASCII files where we have to store the raw values to a newly created vtk data object)
  prevTime = currTime;
  //////////////////////////////////////////////////////////////////////////

  std::vector<vtkm::Float32>::size_type nDims = 0;
  vtkm::cont::DataSet inDataSet;
  std::vector<ValueType> values;
  // TODO: for ContourTreeDelaunayApp.cxx only support irregular meshes, ...
  // ... do not check dimensions / disable marching cubes as a parameter
  std::vector<vtkm::Id> dims;

  printMemoryUsage("[ContourTreeApp.cxx] BEFORE READING IN VTK");

  if (filename.compare(filename.length() - 3, 3, "bov") == 0)
  {
    vtkm::io::BOVDataSetReader reader(filename);
    inDataSet = reader.ReadDataSet();
    nDims = 3;
    ///TIMING///////////////////////////////////////////////////////
    currTime = totalTime.GetElapsedTime();
    dataReadTimeDisplay = currTime - prevTime;
    prevTime = currTime; // 'buildDatasetTimeDisplay' starts now
    ////////////////////////////////////////////////////////////////
  }
  else if (filename.compare(filename.length() - 3, 3, "bin") == 0)
  {// Read binary data input
      std::cout << "BINARY INPUT\n";
      std::ifstream inFile(filename, std::ios::in | std::ios::binary);
      if (!inFile)
      {
          // Handle file open error
          std::cerr << "Error opening file: " << filename << std::endl;
          return 0;
      }

      // Read the dimensions of the mesh (ASCII)
      std::string header;
      std::getline(inFile, header);
      std::istringstream headerStream(header);
      vtkm::Id dimVertices;
      while (headerStream >> dimVertices)
      {
          dims.push_back(dimVertices);
      }

      // Compute the number of vertices, i.e., xdim * ydim * zdim
      nDims = static_cast<unsigned short>(dims.size());
      std::size_t numVertices = std::accumulate(dims.begin(), dims.end(),
                                                std::size_t(1), std::multiplies<std::size_t>());

      // Rest of the data checks...

      // Read binary data
      values.resize(numVertices);
      inFile.read(reinterpret_cast<char*>(values.data()), numVertices * sizeof(values[0]));

      // Make sure we successfully read the expected amount of data
      if (!inFile)
      {
          // Handle read error
          std::cerr << "Error reading binary data from file: " << filename << std::endl;
          inFile.close();
          return 0;
      }

      // file reading completed, close it to free up some memory:
      inFile.close();

      // ---------------------------------------------------------- //

      ///TIMING///////////////////////////////////////////////////////
      currTime = totalTime.GetElapsedTime();
      dataReadTimeDisplay = currTime - prevTime;
      prevTime = currTime; // 'buildDatasetTimeDisplay' starts now
      ////////////////////////////////////////////////////////////////

      // parse the raw data into VTK-m format (building the dataset):

      // swap dims order
      std::swap(dims[0], dims[1]);

      // build the input dataset
      vtkm::cont::DataSetBuilderUniform dsb;
      // 2D data
      if (nDims == 2)
      {
        vtkm::Id2 vdims;
        vdims[0] = static_cast<vtkm::Id>(dims[0]);
        vdims[1] = static_cast<vtkm::Id>(dims[1]);
        inDataSet = dsb.Create(vdims);
      }
      // 3D data
      else
      {
        vtkm::Id3 vdims;
        vdims[0] = static_cast<vtkm::Id>(dims[0]);
        vdims[1] = static_cast<vtkm::Id>(dims[1]);
        vdims[2] = static_cast<vtkm::Id>(dims[2]);
        inDataSet = dsb.Create(vdims);
      }
      inDataSet.AddPointField("values", values);

      /// DEBUG PRINT std::cout << "inDataSet ASCII summary\n";
      inDataSet.PrintSummary(std::cout);

      // DEBUG SLEEP std::this_thread::sleep_for(std::chrono::seconds(10));
  }
  else if (filename.compare(filename.length() - 3, 3, "vtk") == 0)
  {
      std::cout << "VTK file (a Delaunay output by TetGen expected): " << filename << std::endl;
      // const std::string filename_vtk = "/home/sc17dd/modules/HCTC2024/VTK-m-topology-refactor/VTK-m_TopologyGraph-Refactor/examples/contour-visualiser/delaunay-parcels/200k-from-2M-sampled-excel-sorted.1-withvalues-manual.vtk";
      vtkm::io::VTKDataSetReader reader(filename);

      // read the data from a VTK file:
      reader.PrintSummary(std::cout);
      inDataSet = reader.ReadDataSet();

      ///TIMING///////////////////////////////////////////////////////
      currTime = totalTime.GetElapsedTime();
      dataReadTimeDisplay = currTime - prevTime;
      prevTime = currTime; // 'buildDatasetTimeDisplay' starts now
      ////////////////////////////////////////////////////////////////

      // Explicitly interpret as tetrahedral cell set
      // Expecting only single type (tetrahedral) cells from the TetGen Delaunay tetrahedralisation ...
      // ... hence check the type early to match vtkm::cont::CellSetSingleType<>
      if ( !inDataSet.GetCellSet().IsType< vtkm::cont::CellSetSingleType<> >() )
      {
          std::cerr << "Dataset is NOT CellSetSingleType. Check input!" << std::endl;
          return 0;
      }

  }
  else if (filename.compare(filename.length() - 3, 3, "foo") == 0)
  {
    std::cout << "Foo file: " << filename << "\n";
    // fake file, no file reading:
    ///TIMING///////////////////////////////////////////////////////
    currTime = totalTime.GetElapsedTime();
    dataReadTimeDisplay = currTime - prevTime;
    prevTime = currTime; // 'buildDatasetTimeDisplay' starts now
    ////////////////////////////////////////////////////////////////

    // build the input dataset
    vtkm::cont::DataSetBuilderUniform dsb;

    int dimsx = 2;
    values.resize(dimsx*dimsx*dimsx);

    vtkm::Id3 vdims;
    vdims[0] = static_cast<vtkm::Id>(dimsx);
    vdims[1] = static_cast<vtkm::Id>(dimsx);
    vdims[2] = static_cast<vtkm::Id>(dimsx);

    inDataSet = dsb.Create(vdims);

    inDataSet.AddPointField("values", values);

    /// DEBUG PRINT std::cout << "inDataSet ASCII summary\n";
    inDataSet.PrintSummary(std::cout);
  }
  else // Read ASCII data input
  {
    std::cout << "Reading file: " << filename << "\n";
    std::ifstream inFile(filename);
    if (inFile.bad())
      return 0;

    // Read the dimensions of the mesh, i.e,. number of elementes in x, y, and z
    std::string line;
    getline(inFile, line);
    std::istringstream linestream(line);
    vtkm::Id dimVertices;
    while (linestream >> dimVertices)
    {
      dims.push_back(dimVertices);
    }

    // Compute the number of vertices, i.e., xdim * ydim * zdim
    nDims = static_cast<unsigned short>(dims.size());
    std::size_t numVertices = static_cast<std::size_t>(
      std::accumulate(dims.begin(), dims.end(), std::size_t(1), std::multiplies<std::size_t>()));

    // Check for fatal input errors
    // Check the the number of dimensiosn is either 2D or 3D
    bool invalidNumDimensions = (nDims < 2 || nDims > 3);
    // Log any errors if found on rank 0
    VTKM_LOG_IF_S(vtkm::cont::LogLevel::Error,
                  invalidNumDimensions && (rank == 0),
                  "The input mesh is " << nDims << "D. The input data must be either 2D or 3D.");
    // If we found any errors in the setttings than finalize MPI and exit the execution
    if (invalidNumDimensions)
    {
#ifdef WITH_MPI
      MPI_Finalize();
#endif
      return EXIT_SUCCESS;
    }

    // Read data
    values.resize(numVertices);
    for (std::size_t vertex = 0; vertex < numVertices; ++vertex)
    {
      inFile >> values[vertex];
    }

    // file reading completed, close it to free up some memory:
    inFile.close();


    ///TIMING///////////////////////////////////////////////////////
    currTime = totalTime.GetElapsedTime();
    dataReadTimeDisplay = currTime - prevTime;
    prevTime = currTime; // 'buildDatasetTimeDisplay' starts now
    ////////////////////////////////////////////////////////////////


    // parse the raw data into VTK-m format (building the dataset):

    // swap dims order
    std::swap(dims[0], dims[1]);

    // build the input dataset
    vtkm::cont::DataSetBuilderUniform dsb;
    // 2D data
    if (nDims == 2)
    {
      vtkm::Id2 vdims;
      vdims[0] = static_cast<vtkm::Id>(dims[0]);
      vdims[1] = static_cast<vtkm::Id>(dims[1]);
      inDataSet = dsb.Create(vdims);
    }
    // 3D data
    else
    {
      vtkm::Id3 vdims;
      vdims[0] = static_cast<vtkm::Id>(dims[0]);
      vdims[1] = static_cast<vtkm::Id>(dims[1]);
      vdims[2] = static_cast<vtkm::Id>(dims[2]);
      inDataSet = dsb.Create(vdims);
    }
    inDataSet.AddPointField("values", values);

    /// DEBUG PRINT std::cout << "inDataSet ASCII summary\n";
    inDataSet.PrintSummary(std::cout);

    // DEBUG_PRINT std::this_thread::sleep_for(std::chrono::seconds(10));

  } // END ASCII Read


#ifndef WITH_MPI                              // construct regular, single-block VTK-M input dataset
  vtkm::cont::DataSet useDataSet = inDataSet; // Single block dataset
#else  // Create a multi-block dataset for multi-block DIY-paralle processing
  // Determine split
  vtkm::Id3 globalSize = nDims == 3 ? vtkm::Id3(static_cast<vtkm::Id>(dims[0]),
                                                static_cast<vtkm::Id>(dims[1]),
                                                static_cast<vtkm::Id>(dims[2]))
                                    : vtkm::Id3(static_cast<vtkm::Id>(dims[0]),
                                                static_cast<vtkm::Id>(dims[1]),
                                                static_cast<vtkm::Id>(1));
  vtkm::Id3 blocksPerDim = ComputeNumberOfBlocksPerAxis(globalSize, numBlocks);
  vtkm::Id blocksPerRank = numBlocks / size;
  vtkm::Id numRanksWithExtraBlock = numBlocks % size;
  vtkm::Id blocksOnThisRank, startBlockNo;
  if (rank < numRanksWithExtraBlock)
  {
    blocksOnThisRank = blocksPerRank + 1;
    startBlockNo = (blocksPerRank + 1) * rank;
  }
  else
  {
    blocksOnThisRank = blocksPerRank;
    startBlockNo = numRanksWithExtraBlock * (blocksPerRank + 1) +
      (rank - numRanksWithExtraBlock) * blocksPerRank;
  }

  if (blocksOnThisRank != 1)
  {
    std::cerr << "Currently only one block per rank supported!";
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // Created partitioned (split) data set
  vtkm::cont::PartitionedDataSet useDataSet;
  vtkm::cont::ArrayHandle<vtkm::Id3> localBlockIndices;
  localBlockIndices.Allocate(blocksPerRank);
  auto localBlockIndicesPortal = localBlockIndices.WritePortal();

  for (vtkm::Id blockNo = 0; blockNo < blocksOnThisRank; ++blockNo)
  {
    vtkm::Id3 blockOrigin, blockSize, blockIndex;
    std::tie(blockIndex, blockOrigin, blockSize) =
      ComputeBlockExtents(globalSize, blocksPerDim, startBlockNo + blockNo);
    useDataSet.AppendPartition(CreateSubDataSet(inDataSet, blockOrigin, blockSize, "values"));
    localBlockIndicesPortal.Set(blockNo, blockIndex);
  }
#endif // WITH_MPI construct input dataset

  printMemoryUsage("[ContourTreeApp.cxx] AFTER READING IN VTK");

  ////////////////////////////////////////////////////////////
  // END MAIN-2 Read the input data (buildDatasetTimeDisplay) //
  ////////////////////////////////////////////////////////////

  // Wrap up the timing of building the dataset, ...
  // ... start the contour tree timer for the next stage
  ///TIMING///////////////////////////////////////////////////////
  currTime = totalTime.GetElapsedTime();
  // TIMING: Build VTKM Dataset (finish) - timing saved to 'buildDatasetTimeDisplay'
  buildDatasetTimeDisplay = currTime - prevTime;
  // TIMING: Compute Contour Tree (begin) - timing will be saved to 'computeContourTreeTimeDisplay'
  vtkm::Float64 computeContourTreeTimeDisplay = 0;
  prevTime = currTime;
  ////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////
  // START MAIN-3 Compute Contour Tree (computeContourTreeTimeDisplay) //
  /////////////////////////////////////////////////////////////////

  std::cout << ORANGE << std::endl;
  std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "// MAIN-3 Compute Contour Tree (computeContourTreeTimeDisplay) //" << std::endl;
  std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << RESET << std::endl;

  // Convert the mesh of values into contour tree, pairs of vertex ids
  vtkm::filter::scalar_topology::ContourTreeAugmented filter(useMarchingCubes,
                                                             augmentComputeRegularStructure);

#ifdef WITH_MPI
  filter.SetBlockIndices(blocksPerDim, localBlockIndices);
#endif
//  filter.SetActiveField("values");
  filter.SetActiveField("var"); // ContourTreeDelaunayApp.cxx

  // Execute the contour tree analysis. NOTE: If MPI is used the result  will be
  // a vtkm::cont::PartitionedDataSet instead of a vtkm::cont::DataSet
  std::cout << "FILTER.EXECUTE(useDataSet) ... \n";
  auto result = filter.Execute(useDataSet);
  std::cout << "... DONE: FILTER.EXECUTE(useDataSet) ... \n";

  ///TIMING///////////////////////////////////////////////////////
  currTime = totalTime.GetElapsedTime();
  // TIMING: Compute Contour Tree (finish) - timing saved to 'computeContourTreeTimeDisplay'
  computeContourTreeTimeDisplay = currTime - prevTime;
  // TIMING: Compute Branch Decomposition (begin) - timing will be saved to 'computeSimplifyBranchDecompTimeDisplay'
  vtkm::Float64 computeSimplifyBranchDecompTimeDisplay = 0;
  prevTime = currTime;
  ////////////////////////////////////////////////////////////////

#ifdef WITH_MPI
#ifdef DEBUG_PRINT
  std::cout << std::flush;
  close(out);
  std::cerr << std::flush;
  close(err);

  dup2(save_out, fileno(stdout));
  dup2(save_err, fileno(stderr));

  close(save_out);
  close(save_err);
#endif
#endif












  //////////////////////////////////////////////////////////////////////////////////////
  // MAIN-3.5 Compute the branch decomposition (computeSimplifyBranchDecompTimeDisplay) //
  //////////////////////////////////////////////////////////////////////////////////////

  std::cout << ORANGE << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "// MAIN-3.5 Compute the branch decomposition (computeSimplifyBranchDecompTimeDisplay) //" << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << RESET << std::endl;

  int compute_betti = 1;

  if (rank == 0 && compute_betti)
  {
      // Time branch decompostion: Volume Weight Computation:
      vtkm::cont::Timer computeVolumeWeightsTimeDisplay;
      computeVolumeWeightsTimeDisplay.Start();

      // ---------------------------- FLOAT WEIGHTS ---------------------------- //

      FloatArrayType superarcIntrinsicBetti;
      FloatArrayType superarcDependentBetti;
      FloatArrayType supernodeTransferBetti;
      FloatArrayType hyperarcDependentBetti;

      vtkm::cont::ArrayHandle<Coefficients> superarcIntrinsicWeightBetti;
      vtkm::cont::ArrayHandle<Coefficients> superarcDependentWeightBetti;
      vtkm::cont::ArrayHandle<Coefficients> supernodeTransferWeightBetti;
      vtkm::cont::ArrayHandle<Coefficients> hyperarcDependentWeightBetti;

  #if DEBUG_PRINT_PACTBD
      std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
  #endif
      std::cout << YELLOW << "\n[STAGE 3.5 Start - Coeff. Weights (IDThD)] ContourTreeApp.cxx:ComputeVolumeWeightsSerialStructCoefficients START ..." << RESET << std::endl;
  #if SLEEP_ON
      std::this_thread::sleep_for(std::chrono::seconds(3));
  #endif

      //    CALLGRIND_START_INSTRUMENTATION;
      ctaug_ns::ProcessContourTree::ComputeBettiNumbersForRegularArcs(useDataSet,
                                                                      filter.GetContourTree(),
                                                                      filter.GetNumIterations(),
                                                                      // The following four outputs are the coefficient tuples
                                                                      // (such as h1, h2, h3, h4 pairs)
                                                                      superarcIntrinsicWeightBetti,  // (output)
                                                                      superarcDependentWeightBetti,  // (output)
                                                                      supernodeTransferWeightBetti,  // (output)
                                                                      hyperarcDependentWeightBetti,
                                                                      // 2025-01-30 added additional output ...
                                                                      // ... to have access to "collapsed" TODO termdefine
                                                                      // ("collapsed" = computed single value weight, ...
                                                                      //  ... instead of N-length coefficient tuples)
                                                                      // These "collapsed" weights are used for ...
                                                                      // ... computing branch weights without relying on ...
                                                                      // ... the node count on the branches
                                                                      superarcIntrinsicBetti,  // (output)
                                                                      superarcDependentBetti,  // (output)
                                                                      supernodeTransferBetti,  // (output)
                                                                      hyperarcDependentBetti); // (output)

      //    CALLGRIND_STOP_INSTRUMENTATION;
      //    CALLGRIND_DUMP_STATS;

      std::cout << YELLOW << "[STAGE 3.5 End - Coeff. Weights (IDThD)] ContourTreeApp.cxx:ComputeVolumeWeightsSerialStructCoefficients ... END\n" << RESET << std::endl;

  }

























  //////////////////////////////////////////////////////////////////////////////////////
  // MAIN-4 Compute the branch decomposition (computeSimplifyBranchDecompTimeDisplay) //
  //////////////////////////////////////////////////////////////////////////////////////

  std::cout << ORANGE << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "// MAIN-4 Compute the branch decomposition (computeSimplifyBranchDecompTimeDisplay) //" << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << RESET << std::endl;

  if (rank == 0 && computeBranchDecomposition && augmentComputeRegularStructure)
  {
    // Time branch decompostion: Volume Weight Computation:
    vtkm::cont::Timer computeVolumeWeightsTimeDisplay;
    computeVolumeWeightsTimeDisplay.Start();

    // ---------------------------- FLOAT WEIGHTS ---------------------------- //

    FloatArrayType superarcIntrinsicWeightNEW;
    FloatArrayType superarcDependentWeightNEW;
    FloatArrayType supernodeTransferWeightNEW;
    FloatArrayType hyperarcDependentWeightNEW;

    vtkm::cont::ArrayHandle<Coefficients> superarcIntrinsicWeightCoeffs;
    vtkm::cont::ArrayHandle<Coefficients> superarcDependentWeightCoeffs;
    vtkm::cont::ArrayHandle<Coefficients> supernodeTransferWeightCoeffs;
    vtkm::cont::ArrayHandle<Coefficients> hyperarcDependentWeightCoeffs;

#if DEBUG_PRINT_PACTBD
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
#endif
    std::cout << YELLOW << "\n[STAGE 4.1 Start - Coeff. Weights (IDThD)] ContourTreeApp.cxx:ComputeVolumeWeightsSerialStructCoefficients START ..." << RESET << std::endl;
#if SLEEP_ON
    std::this_thread::sleep_for(std::chrono::seconds(3));
#endif

    //    CALLGRIND_START_INSTRUMENTATION;
    ctaug_ns::ProcessContourTree::ComputeVolumeWeightsSerialStructCoefficients(useDataSet,
                                                                              filter.GetContourTree(),
                                                                              filter.GetNumIterations(),
                                                                              // The following four outputs are the coefficient tuples
                                                                              // (such as h1, h2, h3, h4 pairs)
                                                                              superarcIntrinsicWeightCoeffs,  // (output)
                                                                              superarcDependentWeightCoeffs,  // (output)
                                                                              supernodeTransferWeightCoeffs,  // (output)
                                                                              hyperarcDependentWeightCoeffs,
                                                                              // 2025-01-30 added additional output ...
                                                                              // ... to have access to "collapsed" TODO termdefine
                                                                              // ("collapsed" = computed single value weight, ...
                                                                              //  ... instead of N-length coefficient tuples)
                                                                              // These "collapsed" weights are used for ...
                                                                              // ... computing branch weights without relying on ...
                                                                              // ... the node count on the branches
                                                                              superarcIntrinsicWeightNEW,  // (output)
                                                                              superarcDependentWeightNEW,  // (output)
                                                                              supernodeTransferWeightNEW,  // (output)
                                                                              hyperarcDependentWeightNEW); // (output)

    //    CALLGRIND_STOP_INSTRUMENTATION;
    //    CALLGRIND_DUMP_STATS;

    std::cout << YELLOW << "[STAGE 4.1 End - Coeff. Weights (IDThD)] ContourTreeApp.cxx:ComputeVolumeWeightsSerialStructCoefficients ... END\n" << RESET << std::endl;

#if DEBUG_PRINT_PACTBD
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
#endif

    // ---------------------------- FLOAT WEIGHTS ---------------------------- //

    std::stringstream timingsStream; // Use a string stream to log in one message
    timingsStream << std::endl;
    timingsStream << "    --------------- Branch Decomposition Timings " << rank
                  << " --------------" << std::endl;
    timingsStream << "    " << std::setw(50) << std::left << "Compute Volume Weights"
                  << ": " << RED << computeVolumeWeightsTimeDisplay.GetElapsedTime() << " seconds" << RESET << std::endl;

    VTKM_LOG_S(vtkm::cont::LogLevel::Warn, timingsStream.str());

    // Time branch decompostion 2/3: Branch Decomposition:
    vtkm::cont::Timer computeBranchDecompositionTimeDisplay;
    computeBranchDecompositionTimeDisplay.Start();

    // compute the branch decomposition by volume
    ctaug_ns::IdArrayType whichBranch;
    ctaug_ns::IdArrayType branchMinimum;
    ctaug_ns::IdArrayType branchMaximum;
    ctaug_ns::IdArrayType branchSaddle;
    ctaug_ns::IdArrayType branchParent;


#if DEBUG_PRINT_PACTBD
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
#endif

    std::cout << YELLOW << "\n[STAGE 4.2 Start - Branch Decomposition Arrays] ContourTreeApp.cxx:ComputeVolumeBranchDecompositionSerialFloat() START ..." << RESET << std::endl;

    ctaug_ns::ProcessContourTree::ComputeVolumeBranchDecompositionSerialFloat(filter.GetContourTree(),
                                                                              superarcDependentWeightNEW,
                                                                              superarcIntrinsicWeightNEW,
                                                                              whichBranch,   // (output)
                                                                              branchMinimum, // (output)
                                                                              branchMaximum, // (output)
                                                                              branchSaddle,  // (output)
                                                                              branchParent); // (output)

    std::cout << YELLOW << "[STAGE 4.2 End - Branch Decomposition Arrays] ContourTreeApp.cxx:ComputeVolumeBranchDecompositionSerialFloat() ... END\n" << RESET << std::endl;

#if DEBUG_PRINT_PACTBD
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
#endif

#if WRITE_FILES

    auto whichBranchPortal   = whichBranch.ReadPortal();
    auto branchMinimumPortal = branchMinimum.ReadPortal();
    auto branchMaximumPortal = branchMaximum.ReadPortal();
    auto branchSaddlePortal  = branchSaddle.ReadPortal();
    auto branchParentPortal  = branchParent.ReadPortal();

    std::cout << "(ContourTreeApp.cxx) Printing the arrays output from the Branch Decomposition to 'ContourTreeGraph--original-fullCT-BRANCH-COLLAPSED.txt'" << std::endl;
    std::cout << "(ContourTreeApp.cxx) whichBranch:";

    std::ofstream file("ContourTreeGraph--original-fullCT-SN-TO-BRANCH-COLLAPSED.txt");

    for (vtkm::Id branchID = 0; branchID < whichBranch.GetNumberOfValues(); branchID++)
    {
        // std::cout << branchID << " = " << whichBranchPortal.Get(branchID) << std::endl;
        file << branchID << "," << whichBranchPortal.Get(branchID) << std::endl;
    }
    file.close();
#endif


    // Record and log the branch decompostion timings
    timingsStream << "    " << std::setw(50) << std::left << "Compute Volume Branch Decomposition"
                  << ": " << computeBranchDecompositionTimeDisplay.GetElapsedTime() << " seconds" << std::endl;
//    VTKM_LOG_S(vtkm::cont::LogLevel::Info, timingsStream.str());
    VTKM_LOG_S(vtkm::cont::LogLevel::Warn, timingsStream.str());



    /// DEBUG PRINT std::cout << "... Computing the Branch Decomposition: LOGGING DONE\n";

    //----main branch decompostion end
    //----Isovalue seleciton start

    if (numLevels > 0) // if compute isovalues
    {
        // Time branch decompostion 2/3: Branch Decomposition:
        vtkm::cont::Timer computeIsovaluesFromBranchDecompositionTimeDisplay;
        computeIsovaluesFromBranchDecompositionTimeDisplay.Start();

      // Get the data values for computing the explicit branch decomposition
//      vtkm::cont::ArrayHandle<ValueType> dataField;
#ifdef WITH_MPI
      result.GetPartitions()[0].GetField("values").GetData().AsArrayHandle(dataField);
      bool dataFieldIsSorted = true;
#else
//      useDataSet.GetField("values").GetData().AsArrayHandle(dataField);
//      bool dataFieldIsSorted = false;

      bool dataFieldIsSorted = true;

      std::vector<vtkm::FloatDefault> std_nodes_sorted;
      // PACTBD-EDIT-FIXED
//      for(vtkm::FloatDefault i = 0.f; i < 29791.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 16.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 9.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 102.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 1002.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 10002.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 99973.f; i += 1.f)

      float num_datapoints_sort = (float)useDataSet.GetPointField("var").GetNumberOfValues();

//      for(vtkm::FloatDefault i = 0.f; i < 200002.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 985182.f; i += 1.f)
//      for(vtkm::FloatDefault i = 0.f; i < 2160931.f; i += 1.f)
      for(vtkm::FloatDefault i = 0.f; i < num_datapoints_sort+1.f; i += 1.f)
      {
        std_nodes_sorted.push_back(i);
      }

      vtkm::cont::ArrayHandle<vtkm::FloatDefault> dataField =
        vtkm::cont::make_ArrayHandle(std_nodes_sorted, vtkm::CopyFlag::Off);
#if DEBUG_PRINT_PACTBD
      for(unsigned i = 0; i < dataField.GetNumberOfValues(); i++)
      {
          std::cout << dataField.ReadPortal().Get(i) << " ";
//          dataField.WritePortal().Set(i, i);

      }

      std::cout << "\n";
#endif

#endif

#if DEBUG_PRINT_PACTBD
      std::cout << std::endl << "(ContourTreeApp.cxx) Superarc Intrinsic from decomposition:" << std::endl;
      for(int i = 0; i < superarcIntrinsicWeightNEWPortal.GetNumberOfValues(); i++)
      {
          std::cout << i << " -> " << superarcIntrinsicWeightNEWPortal.Get(i) << std::endl;
      }
#endif

      /// DEBUG PRINT std::cout << "... Computing the Branch Decomposition: create explicit representation of the branch decompostion from the array representation\n";

      std::cout << YELLOW << "\n[STAGE 4.3 Start - Explicit Branch Decomposition Aggregate Volume Weights] ctaug_ns::ProcessContourTree::ComputeBranchDecomposition<ValueType>() START ..." << RESET << std::endl;

      std::cout << "(ContourTreeApp.cxx) -Branch.h->ComputeBranchDecomposition " << std::endl;
      std::cout << "(ContourTreeApp)->ProcessContourTree->Branch.h->ComputeBranchDecomposition()" << std::endl;

      BranchType* branchDecompostionRoot =
          ctaug_ns::ProcessContourTree::ComputeBranchDecomposition<ValueType>(
            filter.GetContourTree().Superparents,
            filter.GetContourTree().Supernodes,
            filter.GetContourTree().Superarcs,
            whichBranch,
            branchMinimum,
            branchMaximum,
            branchSaddle,
            branchParent,
            filter.GetSortOrder(),
            dataField,
            dataFieldIsSorted,
            superarcIntrinsicWeightNEW,   // used to use manually set values for BD: superarcIntrinsicWeightCorrect,
            superarcDependentWeightNEW); // used to use manually set values for BD: superarcDependentWeightCorrect );

      // The preceding is taken from ProcessContourTree.h and hardcoded here for testing
      std::cout << "(ContourTreeApp)->ProcessContourTree->Branch.h->ComputeBranchDecomposition()" << std::endl;
      std::cout << YELLOW << "[STAGE 4.3 End - Explicit Branch Decomposition Aggregate Volume Weights] ctaug_ns::ProcessContourTree::ComputeBranchDecomposition<ValueType>() ... END\n" << RESET << std::endl;

      /// DEBUG PRINT
#if DEBUG_PRINT_PACTBD
      std::cout << "(ContourTreeApp) Computing the Branch Decomposition: PRINTING\n";
      branchDecompostionRoot->PrintBranchDecomposition(std::cout);
#endif




#if WRITE_FILES
std::cout << "(ContourTreeApp) PRINTING DOT FORMAT: The Branch Decomposition:\n";
// FILE IO START
      file_io_counter++;
      VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
                 std::endl
                   << file_io_counter << ") WRITING: (PrintDot) ContourTreeGraph--branch-decomposition-fullCT.gv" << std::endl);
      std::vector<vtkm::Id> saddle_rootingFullBD = std::vector<vtkm::Id>();
      std::vector<vtkm::Id> local_branchesFullBD = std::vector<vtkm::Id>();
      std::vector<vtkm::Id> depth_FullBD = std::vector<vtkm::Id>();
      std::vector<vtkm::FloatDefault>    branch_weightsFullBD = std::vector<vtkm::FloatDefault>();
      std::vector<vtkm::FloatDefault>    branch_weights_write_FullBD = std::vector<vtkm::FloatDefault>();
      std::vector<bool>     main_branch_flags_FullBD = std::vector<bool>();
      std::vector<vtkm::Id> depth_write_FullBD = std::vector<vtkm::Id>();


      std::ofstream filegvbdfullBD("ContourTreeGraph--branch-decomposition-fullCT.gv");
      branchDecompostionRoot->PrintDotBranchDecomposition(useDataSet, filegvbdfullBD, saddle_rootingFullBD, local_branchesFullBD,
                                                          depth_FullBD, branch_weightsFullBD, branch_weights_write_FullBD,
                                                          main_branch_flags_FullBD, depth_write_FullBD, 0, 0);
// FILE IO END
#endif


#if WRITE_FILES
// FILE IO START
      file_io_counter++;
            VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
                       std::endl
                         << file_io_counter << ") WRITING: (PrintBranch) ContourTreeGraph--branch-decomposition-fullCT.txt" << std::endl);
      std::ofstream filegvbdfull("ContourTreeGraph--branch-decomposition-fullCT.txt");
      branchDecompostionRoot->PrintBranchDecomposition(filegvbdfull);
// FILE IO END
#endif

#ifdef DEBUG_PRINT
      branchDecompostionRoot->PrintBranchDecomposition(std::cout);
#endif

      // Simplify the contour tree of the branch decompostion
      std::cout << std::endl;
      std::cout << "(ContourTreeApp) APPLYING BRANCH SIMPLIFICATION (BrS):\n";

      std::cout << YELLOW << "\n[STAGE 4.4 Start - Simplifying to N branches] ContourTreeApp.cxx:branchDecompostionRoot->SimplifyToSize() START ..." << RESET << std::endl;




      usePersistenceSorter = false;
//      branchDecompostionRoot->SimplifyToSize(2, usePersistenceSorter);
      branchDecompostionRoot->SimplifyToSize(10, usePersistenceSorter);
//       branchDecompostionRoot->SimplifyToSize(numComp, usePersistenceSorter);
      /// DEBUG PRINT
       std::cout << YELLOW << "[STAGE 4.4 End - Simplifying to N branches] ContourTreeApp.cxx:branchDecompostionRoot->SimplifyToSize() ... END\n" << RESET << std::endl;

#if WRITE_FILES
// FILE IO START
      std::cout << "(REFACTOR VERSION)\n";
      std::cout << "(ContourTreeApp) Computing the Branch Decomposition: PRINTING AFTER SIMPLIFICATION\n";

      std::vector<vtkm::Id> saddle_rooting = std::vector<vtkm::Id>();
      std::vector<vtkm::Id> local_branches = std::vector<vtkm::Id>();
      std::vector<vtkm::Id> depth = std::vector<vtkm::Id>();
      std::vector<vtkm::FloatDefault>    branch_weights = std::vector<vtkm::FloatDefault>();
      std::vector<vtkm::FloatDefault>    branch_weights_write = std::vector<vtkm::FloatDefault>();
      std::vector<bool>     main_branch_flags = std::vector<bool>();
      std::vector<vtkm::Id> depth_write = std::vector<vtkm::Id>();

      file_io_counter++;
      branchDecompostionRoot->PrintBranchDecomposition(std::cout);
            VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
                       std::endl
                         << file_io_counter << ") WRITING: ContourTreeGraph--branch-decomposition-simplifiedCT.gv" << std::endl);
      std::ofstream filegvbdsimplified("ContourTreeGraph--branch-decomposition-simplifiedCT.gv");
      branchDecompostionRoot->PrintDotBranchDecomposition(useDataSet, filegvbdsimplified, saddle_rooting,
                                                          local_branches, depth, branch_weights, branch_weights_write,
                                                          main_branch_flags, depth_write, 0, 0);
// FILE IO END
#endif

      // Compute the relevant iso-values
      std::vector<ValueType> isoValues;
      switch (contourSelectMethod)
      {
        default:
        case 0:
        {
          branchDecompostionRoot->GetRelevantValues(static_cast<int>(contourType), eps, isoValues);
        }
        break;
        case 1:
        {
          vtkm::worklet::contourtree_augmented::process_contourtree_inc::PiecewiseLinearFunction<
            ValueType>
            plf;
          branchDecompostionRoot->AccumulateIntervals(static_cast<int>(contourType), eps, plf);
          isoValues = plf.nLargest(static_cast<unsigned int>(numLevels));
        }
        break;
      }

      // Print the compute iso values
      std::stringstream isoStream; // Use a string stream to log in one message
      isoStream << std::endl;
      isoStream << "    ------------------- Isovalue Suggestions --------------------" << std::endl;
      std::sort(isoValues.begin(), isoValues.end());
      isoStream << "    Isovalues: ";
      for (ValueType val : isoValues)
      {
        isoStream << val << " ";
      }
      isoStream << std::endl;
      // Unique isovalues
      std::vector<ValueType>::iterator it = std::unique(isoValues.begin(), isoValues.end());
      isoValues.resize(static_cast<std::size_t>(std::distance(isoValues.begin(), it)));
      isoStream << "    Unique Isovalues (" << isoValues.size() << "):";
      for (ValueType val : isoValues)
      {
        isoStream << val << " ";
      }
//      VTKM_LOG_S(vtkm::cont::LogLevel::Info, isoStream.str());
      VTKM_LOG_S(vtkm::cont::LogLevel::Warn, isoStream.str());


      std::cout << isoStream.str() << std::endl << std::endl << std::endl << std::endl;

      // 2024-05-28
//      std::ofstream filebdgv("ContourTreeGraph-29K-branch-decomposition-CT.gv");
//      std::ofstream filebdgv("ContourTreeGraph-13k-branch-decomposition-prunedCT.txt");
//      std::ofstream filebdgv("ContourTreeGraph-56M-branch-decomposition-prunedCT.txt");
//      std::ofstream filebdgv("ContourTreeGraph--1024--branch-decomposition-prunedCT.txt");
//      std::ofstream filebdgv("ContourTreeGraph--NastyW-16--branch-decomposition-prunedCT.txt");

#if WRITE_FILES
// FILE IO START
      file_io_counter++;
            VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
                       std::endl
                         << file_io_counter << ") WRITING: ContourTreeGraph--branch-decomposition-simplifiedCT.txt" << std::endl);
      std::ofstream filebdgv("ContourTreeGraph--branch-decomposition-simplifiedCT.txt");
      branchDecompostionRoot->PrintBranchDecomposition(filebdgv);
// FILE IO END
#endif

      timingsStream << "    " << std::setw(50) << std::left << "Compute Isovalues from Branch Decomposition"
                    << ": " << computeIsovaluesFromBranchDecompositionTimeDisplay.GetElapsedTime() << " seconds" << std::endl;
  //    VTKM_LOG_S(vtkm::cont::LogLevel::Info, timingsStream.str());
      VTKM_LOG_S(vtkm::cont::LogLevel::Warn, timingsStream.str());

    } //end if compute isovalue

    ///TIMING///////////////////////////////////////////////////////
    currTime = totalTime.GetElapsedTime();
    // TIMING: Compute Branch Decomposition (finish) - timing saved to 'computeSimplifyBranchDecompTimeDisplay'
    computeSimplifyBranchDecompTimeDisplay = currTime - prevTime;
    prevTime = currTime;
    ////////////////////////////////////////////////////////////////

  }

  // previously getting branch time here, need to move it above before FILE IO
//  currTime = totalTime.GetElapsedTime();
//  // TIMING: Compute Branch Decomposition (finish)
//  vtkm::Float64 computeSimplifyBranchDecompTimeDisplay = currTime - prevTime;
//  prevTime = currTime;

  //vtkm::cont::Field resultField =  result.GetField();
  //vtkm::cont::ArrayHandle<vtkm::Pair<vtkm::Id, vtkm::Id> > saddlePeak;
  //resultField.GetData().AsArrayHandle(saddlePeak);

  // Dump out contour tree for comparison
  if (rank == 0 && printContourTree)
  {
    std::cout << "FINISHED Contour Tree" << std::endl;
    std::cout << "============" << std::endl;
    ctaug_ns::EdgePairArray saddlePeak;
    ctaug_ns::ProcessContourTree::CollectSortedSuperarcs(
      filter.GetContourTree(), filter.GetSortOrder(), saddlePeak);
    std::cout << "ctaug_ns::PrintEdgePairArrayColumnLayout(saddlePeak, std::cout);" << std::endl;
    //ctaug_ns::PrintEdgePairArrayColumnLayout(saddlePeak, std::cout);
    std::cout << "NOTE: ... skipped printing to cout, only printing to the file" << std::endl;

#if WRITE_FILES
// FILE IO START
    file_io_counter++;
    VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
               std::endl
                 << file_io_counter << ") WRITING: ContourTreeGraph--original-fullCT-ColumnFormat.txt" << std::endl);

    std::ofstream file("ContourTreeGraph--original-fullCT-ColumnFormat.txt");


    if (file.is_open()) {
        ctaug_ns::PrintEdgePairArrayColumnLayout(saddlePeak, file);
        file.close(); // Make sure to close the file when done
    } else {
        std::cerr << "Unable to open file for writing.\n";
    }
// FILE IO END
#endif

#if WRITE_FILES
// FILE IO START
    file_io_counter++;
          VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
                     std::endl
                       << file_io_counter << ") WRITING: ContourTreeGraph--original-fullCT.gv" << std::endl);

    std::cout << "Saving the Contour Tree As Dot GraphViz File .gv" << std::endl;
    std::ofstream filegv("ContourTreeGraph--original-fullCT.gv");
    filter.GetContourTree().PrintDotSuperStructure(filegv);
    std::cout << " Finished PrintDotSuperStructure " << std::endl;
// FILE IO END
#endif

//    std::ofstream filegv("ContourTreeGraph-branch-decomposition-LT2M-PACT.gv");
//    branchDecompostionRoot->PrintBranchDecomposition(std::cout);

  }

#ifdef WITH_MPI
  // Force a simple round-robin on the ranks for the summary prints. Its not perfect for MPI but
  // it works well enough to sort the summaries from the ranks for small-scale debugging.
  if (rank > 0)
  {
    int temp;
    MPI_Status status;
    MPI_Recv(&temp, 1, MPI_INT, (rank - 1), 0, comm, &status);
  }
#endif
  currTime = totalTime.GetElapsedTime();
  //VTKM_LOG_S(vtkm::cont::LogLevel::Info,
  VTKM_LOG_S(vtkm::cont::LogLevel::Warn,
         std::endl
           << "    -------------------------- Totals " << rank
           << " -----------------------------" << std::endl
           << std::setw(42) << std::left << "    Start-up"                      << ": " << std::setw(8) << std::left << startUpTimeDisplay             << std::setw(12) << std::left << " seconds" << "(startUpTimeDisplay)"  << std::endl
           << std::setw(42) << std::left << "    Data Read"                     << ": " << std::setw(8) << std::left << dataReadTimeDisplay            << std::setw(12) << std::left << " seconds" << "(dataReadTimeDisplay)" << std::endl
           << std::setw(42) << std::left << "    Build VTKM Dataset"            << ": " << std::setw(8) << std::left << buildDatasetTimeDisplay        << std::setw(12) << std::left << " seconds" << "(buildDatasetTimeDisplay)" << std::endl
           << std::setw(42) << std::left << "    Compute Contour Tree"          << ": " << std::setw(8) << std::left << computeContourTreeTimeDisplay  << std::setw(12) << std::left << " seconds" << "(computeContourTreeTimeDisplay)" << std::endl
           << std::setw(42) << std::left << "    Compute Branch Decomposition"  << ": " << std::setw(8) << std::left << computeSimplifyBranchDecompTimeDisplay << std::setw(12) << std::left << " seconds" << "(computeSimplifyBranchDecompTimeDisplay)" << std::endl
           << std::setw(42) << std::left << "    Total Time"                    << ": " << std::setw(8) << std::left << currTime                << std::setw(12) << std::left << " seconds" << "(currTime)");

  const ctaug_ns::ContourTree& ct = filter.GetContourTree();
  VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
             std::endl
               << "    ---------------- Contour Tree Array Sizes ---------------------" << std::endl
               << ct.PrintArraySizes());
  // Print hyperstructure statistics
  VTKM_LOG_S(vtkm::cont::LogLevel::Warn, //Info,
             std::endl
               << ct.PrintHyperStructureStatistics(false) << std::endl);

  std::cout << "When Transferred - local:" << std::endl;
  auto whenTransferredPortal = ct.WhenTransferred.ReadPortal();
  for(int i = 0; i < whenTransferredPortal.GetNumberOfValues(); i++)
  {
      std::cout << i << "\t" << whenTransferredPortal.Get(i) << "\t"
                << vtkm::worklet::contourtree_augmented::MaskedIndex(whenTransferredPortal.Get(i)) << std::endl;
  }

  // Flush ouput streams just to make sure everything has been logged (in particular when using MPI)
  std::cout << std::flush;
  std::cerr << std::flush;

#ifdef WITH_MPI
  // Let the next rank know that it is time to print their summary.
  if (rank < (size - 1))
  {
    int message = 1;
    MPI_Send(&message, 1, MPI_INT, (rank + 1), 0, comm);
  }
#endif

#ifdef WITH_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
