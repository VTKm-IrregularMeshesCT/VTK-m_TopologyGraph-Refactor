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

// for sleeping
#include <chrono>
#include <thread>

// for memory usage
#include <sys/resource.h>
#include <unistd.h>

// file IO
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#define DEBUG_PRINT_PACTBD 0
#define SLEEP_ON 0
#define PROFILING_PACTBD 1
#define WRITE_FILES 0


namespace process_contourtree_inc_ns =
  vtkm::worklet::contourtree_augmented::process_contourtree_inc;

//    using ValueType = vtkm::Float32;
using ValueType = vtkm::Float64; //vtkm::FloatDefault;
using FloatArrayType = vtkm::cont::ArrayHandle<ValueType>;



namespace vtkm
{
namespace worklet
{
namespace contourtree_augmented
{

struct Coordinates
{
    long double x;
    long double y;
    long double z;

    Coordinates()
    {

    }

    Coordinates(double _x, double _y, double _z)
    {
        x = (long double) _x;
        y = (long double) _y;
        z = (long double) _z;
    }
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
    long double h1;
    long double h2;
    long double h3;
    long double h4;
};

// since in our code vectors are associated with points in space ...
// ... we wrap vtkm vectors, which are directional, to a position vector ...
// ... by keeping track of start and end points
class PositionVector
{
public:
    vtkm::Vec3f_64 start;
    vtkm::Vec3f_64 end;
    // we define difference as end-start
    vtkm::Vec3f_64 difference;
    // direction is just the normalised (unit) difference between start and end
    vtkm::Vec3f_64 direction;


    PositionVector(vtkm::Vec3f_64 aStart, vtkm::Vec3f_64 aEnd)
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

    vtkm::Vec3f_64 lerp2point(double interpolant)
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

    std::cout << "READING Coordinates: " << filename << std::endl;

    while (getline(file, line)) {
        std::istringstream iss(line);
        vtkm::Id id;
        Coordinates coords;

        // Extracting ID and coordinates from the current line, assuming the ID is a vtkm::Id.
        if (!(iss >> coords.x >> coords.y >> coords.z)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue; // Skip to the next line if parsing fails.
        }
#if DEBUG_PRINT_PACTBD
        std::cout << coords.x << " " << coords.y << " " << coords.z << std::endl;
#endif

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

    std::cout << "READING Triangles:" << filename << std::endl;

    while (getline(file, line)) {
        std::istringstream iss(line);
        vtkm::Id id;
        Triangle triang;

        // Extracting ID and coordinates from the current line, assuming the ID is a vtkm::Id.
        if (!(iss >> id >> triang.p1 >> triang.p2 >> triang.p3)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue; // Skip to the next line if parsing fails.
        }
#if DEBUG_PRINT_PACTBD
        std::cout << triang.p1 << " " << triang.p2 << " " << triang.p3 << std::endl;
#endif
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

    std::cout << "READING Tets: " << filename << std::endl;

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
#if DEBUG_PRINT_PACTBD
        std::cout << tet.p1 << " " << tet.p2 << " " << tet.p3 << " " << tet.p4 << std::endl;
#endif

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









    // compute the triangle area from the coordinate points of the triangle:
    double static ComputeTriangleArea(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double x3, double y3, double z3)
//    double static ComputeTriangleArea(Coordinates c1, Coordinates c2, Coordinates c3)
    {
        return 0.5 * abs( (x2-x1)*(y3-y1) - (x3 - x1) * (y2 -y1) );
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

    void static print2Darray(std::vector<std::vector<float>> vxtc)
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




    vtkm::Vec3f_64 static fromTo(vtkm::Vec3f_64 a, vtkm::Vec3f_64 b)
    {
        vtkm::Vec3f_64 start = a;
        vtkm::Vec3f_64 result;

        result = start + (a - b);

        return result;
    }


    void static printMemoryUsage(const std::string& message)
    {
        // Red text formatting for highlighting some console output:
        const std::string ORANGE = "\033[38;2;255;165;0m";  // Start red text
        const std::string LIGHT_BLUE = "\033[38;5;117m";  // Light blue in 256-color
        const std::string RESET = "\033[0m"; // End red text

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);

        std::ifstream status_file("/proc/self/status");
        std::string line;

        size_t current_usage = 0;

        while(std::getline(status_file, line))
        {
            if(line.find("VmRSS:") == 0)
            {
                std::istringstream iss(line);
                std::string key;
                size_t memory; // memory value in kB
                std::string unit;

                iss >> key >> memory >> unit;
                current_usage = memory; // return in KB
            }
        }

        std::cout << ORANGE << message << LIGHT_BLUE << " - Memory usage (peak): " << usage.ru_maxrss
                  << " KB | (current) " << current_usage << " KB" << RESET << std::endl;
    }





    // 2024-11-25 COMPUTE THE STRUCT COEFFICIENTS VERSION OF THE WEIGHTS WITH COEFFICIENTS
    void static ComputeVolumeWeightsSerialStructCoefficients(const vtkm::cont::DataSet& input, // the coefficient-based version additionally requires tetrahedral connections and vertex coordinates
                                                const ContourTree& contourTree,
                                                const vtkm::Id nIterations,
                                                vtkm::cont::ArrayHandle<Coefficients>& superarcIntrinsicWeightCoeff,
                                                vtkm::cont::ArrayHandle<Coefficients>& superarcDependentWeightCoeff,
                                                vtkm::cont::ArrayHandle<Coefficients>& supernodeTransferWeightCoeff,
                                                vtkm::cont::ArrayHandle<Coefficients>& hyperarcDependentWeightCoeff,
                                                // Added 2025-01-30
                                                // We use simple weights for the branch decomposition
                                                FloatArrayType& superarcIntrinsicWeight,
                                                FloatArrayType& superarcDependentWeight,
                                                FloatArrayType& supernodeTransferWeight,
                                                FloatArrayType& hyperarcDependentWeight)
    // 2) COEFFICIENTS:
    { // ContourTreeMaker::ComputeWeights()
      // START ComputeVolumeWeightsSerialStructCoefficients

        std::cout << "[ProcessContourTree.h::ComputeVolumeWeightsSerialStructCoefficients] Compute h1,h2,h3,h4 Volume Weights using Coefficients" << std::endl;
        printMemoryUsage("[ProcessContourTree.h::ComputeVolumeWeightsSerialStructCoefficients] Checkpoint 1/4 - START");

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

        // Files for 3D experiments of Contour Tree Branch volume-based weight computations
        // PACTBD-EDIT-FIXED
//        const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/101-from-2M-sampled-excel-sorted.1-COORDINATES.txt";
//        const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/200k-from-2M-sampled-excel-sorted.1-COORDINATES.txt";
//        const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/1M-from-2M-sampled-excel-sorted.1-COORDINATES.txt";
//        const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/2M-parcels-20250225-sorted.1-valued-COORDINATES.txt";
        // PACTBD-EDIT-FIXED
//        const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/101-from-2M-sampled-excel-sorted.1-TETS.txt";
//        const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/200k-from-2M-sampled-excel-sorted.1-TETS.txt";
//        const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/1M-from-2M-sampled-excel-sorted.1-TETS.txt";
//        const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/2M-parcels-20250225-sorted.1-valued-TETS.txt";

        // get the total number of values to be swept ...
        // ... this will be one more than the total number of datapoints ...
        // ... because we include the region beyond the last isovalue N, as a range [N, +inf)
        // (also, the number of values corresponds to the total number of different data values in a dataset ...
        //  ... because of the simulation of simplicity, which ensures ever data point has a unique value)
        int num_sweep_values = contourTree.Arcs.GetNumberOfValues() + 1;

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {// for each sortedNode
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);

            // initialise the intrinsic weight array counter:
            // superarcIntrinsicWeightPortal.Set(superparent, 0);
            superarcIntrinsicWeightPortal.Set(superparent, 0.f);
        }// for each sortedNode

        // Now move on to set up 3D tetrahedral variables
        // \/\/\/ Timing and performance profiling \/\/\/ //
        // For printing text as red in the console:
        const std::string RED = "\033[31m";
         const std::string RESET = "\033[0m";
        // the same timer will be used thoghout, resetting with GetTimeElapsed ...
        /// ... and restarted for separate sections
        vtkm::cont::Timer timer;
        timer.Start();
        // /\/\/\ Timing and performance profiling /\/\/\ //

//        // Vertices of the tet are given as sort IDs and there are 4 of them:
//        std::vector<std::vector<int>> tetlistSorted(tetlist.size(),
//                                                    std::vector<int> (4, 0));

//        std::cout << "    " << RED << std::setw(38) << std::left << "tetlistSorted allocation"
//                      << ": " << timer.GetElapsedTime() << " seconds" << RESET << std::endl;

//        timer.Start();


        // Keep track of all tetrahedra vertices in a sorted list:
        //  ... for example a tet X, Y, Z, W might have vertices with sort IDs:
        //      X=6, Y=5, Z=3, W=7
        //  ... we sort them in increasing order as 3, 5, 6, 7 ...
        //  ... and refer to vertices as A, B, C, D
        //  ... where we then assign A = 3, B = 5, C = 6, D = 7)
        //  (we do not keep track of the original ordering X, Y, Z, W) ...

        // PACTBD-EDIT-FIXED
//        const std::string filename_vtk = "/home/sc17dd/modules/HCTC2024/VTK-m-topology-refactor/VTK-m_TopologyGraph-Refactor/examples/contour-visualiser/delaunay-parcels/200k-from-2M-sampled-excel-sorted.1-withvalues-manual.vtk";

        cont::ArrayHandle<vtkm::Vec3f> coordinatesVTK;
        coordinatesVTK = input.GetPointField("coordinates").GetData().AsArrayHandle<cont::ArrayHandle<vtkm::Vec3f>>();

        // Explicitly interpret as tetrahedral cell set
        using TetCellSet = vtkm::cont::CellSetSingleType<>;
        // Now safely access connectivity data (moving them up to avoid memory duplication)
        // Get the tetrahedral connectivity array
        // NOTE: At this point it is assumed that the file is 'CellSetSingleType' ...
        //       ... with 'input.GetCellSet().IsType<TetCellSet>()'
        vtkm::cont::ArrayHandle<vtkm::Id> tet_connectivity = input.GetCellSet().AsCellSet<TetCellSet>().GetConnectivityArray(vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{});


#if WRITE_FILES
        std::ofstream fileCoords("ContourTreeGraph--recreate-COORDINATES.txt");
#endif
        auto coordinatesAutoVTK = coordinatesVTK.ReadPortal();
        // Coordinate lists and connectivities of 2D triangle-based and 3D tetrahedral datasets.
        // The coordinates are needed for computing the areas/volumes.
        // Note that the coordinate information is not previously needed for the base contour tree computation ...
        // ... and is only introduced as a requirement at this stage.
        std::vector<Coordinates> coordlist3D;//  = ReadCoordinatesFromFile(filename3D1);
        coordlist3D.reserve(coordinatesVTK.GetNumberOfValues());

        for(int i = 0; i < coordinatesVTK.GetNumberOfValues(); i++)
        {
            coordlist3D.emplace_back(coordinatesAutoVTK.Get(i)[0],
                                     coordinatesAutoVTK.Get(i)[1],
                                     coordinatesAutoVTK.Get(i)[2]);
#if WRITE_FILES
            fileCoords << std::fixed << std::setprecision(13) << coordinatesAutoVTK.Get(i)[0] << " "
                       << std::fixed << std::setprecision(13) << coordinatesAutoVTK.Get(i)[1] << " "
                       << std::fixed << std::setprecision(13) << coordinatesAutoVTK.Get(i)[2] << std::endl;
#endif
        }
#if WRITE_FILES
        fileCoords.close();
#endif

#if WRITE_FILES
        std::ofstream fileTets("ContourTreeGraph--recreate-TETS.txt");
#endif

        // Vertices of the tet are given as sort IDs and there are 4 of them:
        std::vector<std::vector<int>> tetlistSorted(tet_connectivity.GetNumberOfValues()/4,
                                                    std::vector<int> (4, 0));

        std::cout << "    " << RED << std::setw(38) << std::left << "tetlistSorted allocation"
                      << ": " << timer.GetElapsedTime() << " seconds" << RESET << std::endl;

        timer.Start();

        auto tet_connectivityRead = tet_connectivity.ReadPortal();
        for (int i = 0; i < tet_connectivity.GetNumberOfValues(); i+=4)
        {
#if WRITE_FILES
            fileTets << tet_connectivityRead.Get(i)    << " " << tet_connectivityRead.Get(i+1)  << " "
                 << tet_connectivityRead.Get(i+2)  << " " << tet_connectivityRead.Get(i+3)  << std::endl;
#endif
            tetlistSorted[i/4][0] = tet_connectivityRead.Get(i);
            tetlistSorted[i/4][1] = tet_connectivityRead.Get(i+1);
            tetlistSorted[i/4][2] = tet_connectivityRead.Get(i+2);
            tetlistSorted[i/4][3] = tet_connectivityRead.Get(i+3);
        }
#if WRITE_FILES
        fileTets.close();
#endif

        for (int i = 0; i < tet_connectivity.GetNumberOfValues()/4; i++)
        {// for each tet
            std::sort(tetlistSorted[i].begin(), tetlistSorted[i].end());
        }// for each tet

        // Time profiling: node .GetElapsedTime() resets the timer so we can reuse it for the next section with .Start()
        printMemoryUsage("[ProcessContourTree.h::ComputeVolumeWeightsSerialStructCoefficients] Checkpoint 2/4 - BEFORE populating coefficient tables");
        std::cout << "    " << RED << std::setw(38) << std::left << "BEFORE PRE-PROCESS"
                      << ": " << timer.GetElapsedTime() << " seconds" << RESET << std::endl;
        timer.Start();

        vtkm::cont::Timer preprocess_total_timer;
        preprocess_total_timer.Start();

        std::cout << "==========================================PRE-PROCESS=====================================" << std::endl;

        // ----------------------------------- PRE-PROCESS ----------------------------------- //
        // Here we basically populate the table as such:
        // vertexID / triangleID
        //          triangle1 triangle2   row sum(del. h1/h2)   col (prefix) sum del. h1/h2
        // deltas:  h1 h2     h1 h2
        // v1       x  y      z  w
        std::vector<long double> vx_delta_h1_sum;
        vx_delta_h1_sum.reserve(tetlistSorted.size());
        std::vector<long double> vx_delta_h2_sum; // 1D coefficients
        vx_delta_h2_sum.reserve(tetlistSorted.size());
        std::vector<long double> vx_delta_h3_sum; // 2D coefficients
        vx_delta_h3_sum.reserve(tetlistSorted.size());
        std::vector<long double> vx_delta_h4_sum; // 3D coefficients
        vx_delta_h4_sum.reserve(tetlistSorted.size());

        std::cout << "Computing 3D Coefficients ..." << std::endl;

        long double min_volume = 38358000.0l;
        long double max_volume = 0.0l;

        // Preallocation
        // not using the 2D array, writing directly to vx_delta arrays:
//        std::vector<std::vector<long double>> tet_down_deltas_pfix(num_sweep_values, std::vector<long double>(4, 0.0l));

        // before your loop ensure these are initialized correctly:
        vx_delta_h1_sum.assign(num_sweep_values, 0.0l);
        vx_delta_h2_sum.assign(num_sweep_values, 0.0l);
        vx_delta_h3_sum.assign(num_sweep_values, 0.0l);
        vx_delta_h4_sum.assign(num_sweep_values, 0.0l);

        for (vtkm::Id i = 0; i < tetlistSorted.size(); ++i)
        {
            // Step 1
            // Vertices that define the Tetrahedron ABCD (The entire tetrahedron) ...
            // ... with their corresponding isovalues (A holds h1, B - h2, C - h3 and D - h4)
            vtkm::Vec3f_64 verticesA = {coordlist3D[tetlistSorted[i][0]].x, coordlist3D[tetlistSorted[i][0]].y, coordlist3D[tetlistSorted[i][0]].z};
            vtkm::Vec3f_64 verticesB = {coordlist3D[tetlistSorted[i][1]].x, coordlist3D[tetlistSorted[i][1]].y, coordlist3D[tetlistSorted[i][1]].z};
            vtkm::Vec3f_64 verticesC = {coordlist3D[tetlistSorted[i][2]].x, coordlist3D[tetlistSorted[i][2]].y, coordlist3D[tetlistSorted[i][2]].z};
            vtkm::Vec3f_64 verticesD = {coordlist3D[tetlistSorted[i][3]].x, coordlist3D[tetlistSorted[i][3]].y, coordlist3D[tetlistSorted[i][3]].z};

            // Step 2
            // Deriving middle slab triangle vertices E, F, G, H
            // Plane Points at isovalue h=h2 (4) (for interval h1->h2)              - FIRST TET
            PositionVector vAC(verticesA, verticesC), vAD(verticesA, verticesD), vBD(verticesB, verticesD);
            long double lerpADh2 = (long double)(tetlistSorted[i][1] - tetlistSorted[i][0])/(long double)(tetlistSorted[i][3] - tetlistSorted[i][0]);
            vtkm::Vec3f_64 verticesE = (vAD.lerp2point(lerpADh2));
            long double lerpACh2 = (long double)(tetlistSorted[i][1] - tetlistSorted[i][0])/(long double)(tetlistSorted[i][2] - tetlistSorted[i][0]);
            vtkm::Vec3f_64 verticesF = (vAC.lerp2point(lerpACh2));

            // Plane Points at isovalue h=h3 (5) (for interval h4->h3)              - LAST TET
            long double lerpADh3 = (long double)(tetlistSorted[i][2] - tetlistSorted[i][0])/(long double)(tetlistSorted[i][3] - tetlistSorted[i][0]);
            vtkm::Vec3f_64 verticesG = (vAD.lerp2point(lerpADh3));
            long double lerpBDh2 = (long double)(tetlistSorted[i][2] - tetlistSorted[i][1])/(long double)(tetlistSorted[i][3] - tetlistSorted[i][1]);
            vtkm::Vec3f_64 verticesH = (vBD.lerp2point(lerpBDh2));

            // Step 3
            PositionVector a_vol(verticesA, verticesC),
                           b_vol(verticesA, verticesD),
                           c_vol(verticesA, verticesB);
            long double full_tet_vol = (1.0l / 6.0l) * abs(vtkm::Dot(vtkm::Cross(a_vol.difference, b_vol.difference), c_vol.difference));
            min_volume = std::min(min_volume, full_tet_vol);
            max_volume = std::max(max_volume, full_tet_vol);

            // Step 4: Slab 1 Volumes
//                long double alpha_h2 = std::max(0.0l, std::min(1.0l, (long double)(teth2s[i]-teth1s[i])/(teth3s[i]-teth1s[i])));
//                long double beta_h2  = std::max(0.0l, std::min(1.0l, (long double)(teth2s[i]-teth1s[i])/(teth4s[i]-teth1s[i])));
            long double alpha_h2= std::max(0.0l, std::min(1.0l,
                                      (long double)(tetlistSorted[i][1]-tetlistSorted[i][0])/(long double)(tetlistSorted[i][2]-tetlistSorted[i][0])));

            long double beta_h2 = std::max(0.0l, std::min(1.0l,
                                      (long double)(tetlistSorted[i][1]-tetlistSorted[i][0])/(long double)(tetlistSorted[i][3]-tetlistSorted[i][0])));
            PositionVector a_h1h2_vol = a_vol, b_h1h2_vol = b_vol, c_h1h2_vol = c_vol;
            a_h1h2_vol.lerp(alpha_h2);
            b_h1h2_vol.lerp(beta_h2);
            long double slab1_h1h2_vol = (1.0l/6.0l)*abs(vtkm::Dot(vtkm::Cross(a_h1h2_vol.difference,b_h1h2_vol.difference),c_h1h2_vol.difference));

            // Step 5: Slab 1 Coefficients
            long double denom_h1h2 = (long double)( std::pow((-tetlistSorted[i][0] + tetlistSorted[i][1]), 3) ); //pow(-tetlistSorted[i][0] + tetlistSorted[i][1], 3);
            long double a_h1h2 = slab1_h1h2_vol / denom_h1h2;
            long double b_h1h2 = -(3.0l * slab1_h1h2_vol * tetlistSorted[i][0]) / denom_h1h2;
            long double c_h1h2 = (3.0l * slab1_h1h2_vol * std::pow(tetlistSorted[i][0], 2)) / denom_h1h2;
            long double d_h1h2 = -(slab1_h1h2_vol * std::pow(tetlistSorted[i][0], 3)) / denom_h1h2;

            // Step 6,7: Slab 3 Volumes & Coefficients
//                long double alpha_h3 = std::max(0.0l,std::min(1.0l,(long double)(teth3s[i]-teth4s[i])/(teth2s[i]-teth4s[i])));
//                long double beta_h3 = std::max(0.0l,std::min(1.0l,(long double)(teth3s[i]-teth4s[i])/(teth1s[i]-teth4s[i])));
            long double alpha_h3= std::max(0.0l, std::min(1.0l,
                                      (long double)(tetlistSorted[i][2]-tetlistSorted[i][3])/(long double)(tetlistSorted[i][1]-tetlistSorted[i][3])));

            long double beta_h3 = std::max(0.0l, std::min(1.0l,
                                      (long double)(tetlistSorted[i][2]-tetlistSorted[i][3])/(long double)(tetlistSorted[i][0]-tetlistSorted[i][3])));
            PositionVector a_h3h4_vol(verticesD,verticesB), b_h3h4_vol(verticesD,verticesA), c_h3h4_vol(verticesD,verticesC);
            a_h3h4_vol.lerp(alpha_h3);
            b_h3h4_vol.lerp(beta_h3);
            long double slab3_h3h4_vol = (1.0l/6.0l)*abs(vtkm::Dot(vtkm::Cross(a_h3h4_vol.difference,b_h3h4_vol.difference),c_h3h4_vol.difference));
            long double slab3_up_vol = full_tet_vol - slab3_h3h4_vol;
            long double denom_h3h4 = (long double)( std::pow((tetlistSorted[i][2] - tetlistSorted[i][3]), 3) ); //pow(tetlistSorted[i][2]-tetlistSorted[i][3],3);
            long double a_h3h4 = (slab3_up_vol - full_tet_vol)/denom_h3h4;
            long double b_h3h4 = (-3.0l*tetlistSorted[i][3]*slab3_up_vol+3.0l*tetlistSorted[i][3]*full_tet_vol)/denom_h3h4;
            long double c_h3h4 = (3.0l*tetlistSorted[i][3]*tetlistSorted[i][3]*slab3_up_vol-3.0l*tetlistSorted[i][3]*tetlistSorted[i][3]*full_tet_vol)/denom_h3h4;
            long double d_h3h4 = (-pow(tetlistSorted[i][3],3)*slab3_up_vol+pow(tetlistSorted[i][3],3)*full_tet_vol)/denom_h3h4;


            // Step 9 calculations (exactly as in original, unchanged form)
            PositionVector vectorsHG(verticesH, verticesG);
            PositionVector vectorsBE(verticesB, verticesE);
            PositionVector vectorsHC(verticesH, verticesC);
            PositionVector vectorsCG(verticesC, verticesG);
            long double tetk2s = ( 1.0l / (long double)(tetlistSorted[i][2] - tetlistSorted[i][1]) );

            long double n1 = tetk2s * (vectorsHG.mag() - vectorsBE.mag());
            long double n2 = vectorsHC.mag() * tetk2s;
            long double n3 = tetlistSorted[i][1] * tetk2s * vectorsHG.mag() * vectorsHC.mag() * tetk2s;
            long double n4 = tetk2s * tetlistSorted[i][2] * vectorsBE.mag() * vectorsHC.mag() * tetk2s;
            long double a_s1 = n1 * n2;
            long double b_s1 = -( (n1 * n2 * tetlistSorted[i][1]) + n3 - n4 );
            long double c_s1 = n3 * tetlistSorted[i][1] - n4 * tetlistSorted[i][1];

            PositionVector vectorsGH(verticesG, verticesH);
            PositionVector vectorsCH(verticesC, verticesH);
            long double areas_CGH = (1.0l/2.0l) * vtkm::Magnitude(vtkm::Cross(vectorsGH.difference, vectorsCH.difference));
            long double sin_theta_1 = 2.0l * areas_CGH / (vectorsGH.mag() * vectorsCH.mag());

            PositionVector vectorsFB(verticesF, verticesB);
            PositionVector vectorsFE(verticesF, verticesE);
            long double areas_BEF = (1.0l/2.0l) * vtkm::Magnitude(vtkm::Cross(vectorsFB.difference, vectorsFE.difference));
            long double sin_theta_2 = 2.0l * areas_BEF / (vectorsFB.mag() * vectorsFE.mag());

            long double m1 = tetk2s * (vectorsCG.mag() - vectorsFE.mag());
            long double m2 = -vectorsFB.mag() * tetk2s;
            long double m3 = -tetlistSorted[i][1] * tetk2s * vectorsCG.mag() * vectorsFB.mag() * tetk2s;
            long double m4 = tetk2s * tetlistSorted[i][2] * vectorsFE.mag() * vectorsFB.mag() * tetk2s;

            long double a_s2 = m1 * m2;
            long double b_s2 = ( m1 * -m2 * tetlistSorted[i][2] - m3 - m4 );
            long double c_s2 = m3 * tetlistSorted[i][2] + m4 * tetlistSorted[i][2];

            long double a_mid = (sin_theta_1 / 2.0l * a_s1 + sin_theta_2 / 2.0l * a_s2);
            long double b_mid = (sin_theta_1 / 2.0l * b_s1 + sin_theta_2 / 2.0l * b_s2);
            long double c_mid = (sin_theta_1 / 2.0l * c_s1 + sin_theta_2 / 2.0l * c_s2);

            vtkm::Vec3f_64 FExFB_cross_product = vtkm::Cross(vectorsFE.difference, vectorsFB.difference);
            vtkm::Vec3f_64 plane_normal( FExFB_cross_product / (vtkm::Magnitude(FExFB_cross_product)) );
            long double plane_distance = ( vtkm::Magnitude(vtkm::Dot(plane_normal, verticesB) - vtkm::Dot(plane_normal, verticesH)) / vtkm::Magnitude(plane_normal) );
            long double correction_factor = (plane_distance * tetk2s);

            long double a_h2h3 = correction_factor / 3.0l * a_mid;
            long double b_h2h3 = correction_factor / 2.0l * b_mid;
            long double c_h2h3 = correction_factor * c_mid;

            long double d_h2h3 = slab1_h1h2_vol - a_h2h3 * pow(tetlistSorted[i][1],3) - b_h2h3 * pow(tetlistSorted[i][1],2) - c_h2h3 * tetlistSorted[i][1];

            long double d_h2h3_down = ( ( -a_h2h3 * pow(tetlistSorted[i][2], 3) - b_h2h3 * pow(tetlistSorted[i][2], 2) - c_h2h3 * tetlistSorted[i][2] ) -
                                       (-a_h3h4 * pow(tetlistSorted[i][2], 3) - b_h3h4 * pow(tetlistSorted[i][2], 2) - c_h3h4 * tetlistSorted[i][2] - d_h3h4) );

            long double d_h1h2_down = ( ( -a_h1h2 * pow(tetlistSorted[i][1], 3) - b_h1h2 * pow(tetlistSorted[i][1], 2) - c_h1h2 * tetlistSorted[i][1] ) -
                                       (-a_h2h3 * pow(tetlistSorted[i][1], 3) - b_h2h3 * pow(tetlistSorted[i][1], 2) - c_h2h3 * tetlistSorted[i][1] - d_h2h3_down) );

//            // Then directly add to accumulation
//            tet_down_deltas_pfix[tetlistSorted[i][0]][0] += a_h1h2;
//            tet_down_deltas_pfix[tetlistSorted[i][0]][1] += b_h1h2;
//            tet_down_deltas_pfix[tetlistSorted[i][0]][2] += c_h1h2;
//            tet_down_deltas_pfix[tetlistSorted[i][0]][3] += full_tet_vol+d_h1h2_down;

//            tet_down_deltas_pfix[tetlistSorted[i][1]][0] += -a_h1h2+a_h2h3;
//            tet_down_deltas_pfix[tetlistSorted[i][1]][1] += -b_h1h2+b_h2h3;
//            tet_down_deltas_pfix[tetlistSorted[i][1]][2] += -c_h1h2+c_h2h3;
//            tet_down_deltas_pfix[tetlistSorted[i][1]][3] += -d_h1h2_down+d_h2h3_down;

//            tet_down_deltas_pfix[tetlistSorted[i][2]][0] += -a_h2h3+a_h3h4;
//            tet_down_deltas_pfix[tetlistSorted[i][2]][1] += -b_h2h3+b_h3h4;
//            tet_down_deltas_pfix[tetlistSorted[i][2]][2] += -c_h2h3+c_h3h4;
//            tet_down_deltas_pfix[tetlistSorted[i][2]][3] += -d_h2h3_down+d_h3h4;

//            tet_down_deltas_pfix[tetlistSorted[i][3]][0] += -a_h3h4;
//            tet_down_deltas_pfix[tetlistSorted[i][3]][1] += -b_h3h4;
//            tet_down_deltas_pfix[tetlistSorted[i][3]][2] += -c_h3h4;
//            tet_down_deltas_pfix[tetlistSorted[i][3]][3] += -d_h3h4;

            // Directly accumulate, NO extra memory needed
            vx_delta_h1_sum[tetlistSorted[i][0]] += a_h1h2;
            vx_delta_h2_sum[tetlistSorted[i][0]] += b_h1h2;
            vx_delta_h3_sum[tetlistSorted[i][0]] += c_h1h2;
            vx_delta_h4_sum[tetlistSorted[i][0]] += (full_tet_vol + d_h1h2_down);

            vx_delta_h1_sum[tetlistSorted[i][1]] += (-a_h1h2 + a_h2h3);
            vx_delta_h2_sum[tetlistSorted[i][1]] += (-b_h1h2 + b_h2h3);
            vx_delta_h3_sum[tetlistSorted[i][1]] += (-c_h1h2 + c_h2h3);
            vx_delta_h4_sum[tetlistSorted[i][1]] += (-d_h1h2_down + d_h2h3_down);

            vx_delta_h1_sum[tetlistSorted[i][2]] += (-a_h2h3 + a_h3h4);
            vx_delta_h2_sum[tetlistSorted[i][2]] += (-b_h2h3 + b_h3h4);
            vx_delta_h3_sum[tetlistSorted[i][2]] += (-c_h2h3 + c_h3h4);
            vx_delta_h4_sum[tetlistSorted[i][2]] += (-d_h2h3_down + d_h3h4);

            vx_delta_h1_sum[tetlistSorted[i][3]] += (-a_h3h4);
            vx_delta_h2_sum[tetlistSorted[i][3]] += (-b_h3h4);
            vx_delta_h3_sum[tetlistSorted[i][3]] += (-c_h3h4);
            vx_delta_h4_sum[tetlistSorted[i][3]] += (-d_h3h4);

        }

//        for (vtkm::Id i = 0; i < num_sweep_values; i++)
//        {
//            vx_delta_h1_sum.push_back(tet_down_deltas_pfix[i][0]);
//            vx_delta_h2_sum.push_back(tet_down_deltas_pfix[i][1]);
//            vx_delta_h3_sum.push_back(tet_down_deltas_pfix[i][2]);
//            vx_delta_h4_sum.push_back(tet_down_deltas_pfix[i][3]);
//        }










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
//                long double lerpEFGH = (long double)(h - teth2s[i]) / (long double)(teth3s[i] - teth2s[i]);

//            }

            // ==================== /\ Step 3: Deriving middle slab quad vertices P, Q, R, S /\ ==================== //
/* THE PAST REQUIRES SWEEP ISOVALUE H */




        // ----------------------------------- PRE-PROCESS ----------------------------------- //

        printMemoryUsage("[ProcessContourTree.h::ComputeVolumeWeightsSerialStructCoefficients] Checkpoint 3/4 - AFTER populating coefficient tables");
        std::cout << "// ----------------------------------- PRE-PROCESS ----------------------------------- //" << std::endl;
        std::cout << "    " << RED << std::setw(38) << std::left << "TOTAL PRE-PROCESS"
                      << ": " << preprocess_total_timer.GetElapsedTime() << " seconds" << RESET << std::endl;

        // for the sweep, we will be using the pre-computed delta coefficients from the mesh

        // -------------------------------------- SWEEP  ------------------------------------- //


        std::cout << "// -------------------------------------- SWEEP  ------------------------------------- //" << std::endl;

        vtkm::cont::Timer sweep_total_timer;
        sweep_total_timer.Start();

        std::vector<long double> coefficientweightList;
        coefficientweightList.resize(superparentsPortal.GetNumberOfValues());


        std::vector<long double> delta_h1_partial_pfixsum; // = 0.0l;
        std::vector<long double> delta_h2_partial_pfixsum; // = 0.0l;
        std::vector<long double> delta_h3_partial_pfixsum; // = 0.0l;
        std::vector<long double> delta_h4_partial_pfixsum; // = 0.0l;


        delta_h1_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());
        delta_h2_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());
        delta_h3_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());
        delta_h4_partial_pfixsum.resize(num_sweep_values); // 2024-11-18 updated sweep value contourTree.Arcs.GetNumberOfValues());contourTree.Arcs.GetNumberOfValues());


        for(int i = 0; i < coefficientweightList.size(); i++)
        {// for each superparent, initialise four coefficients
            coefficientweightList[i] = 0.0l;
            delta_h1_partial_pfixsum[i] = 0.0l;
            delta_h2_partial_pfixsum[i] = 0.0l;
            delta_h3_partial_pfixsum[i] = 0.0l;
            delta_h3_partial_pfixsum[i] = 0.0l;
        }// for each superparent, initialised four coefficients

        auto arcsPortal = contourTree.Arcs.ReadPortal();
        auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto supernodesPortal = contourTree.Supernodes.ReadPortal();

        vtkm::Id prevIndex = -1;
        vtkm::Id superNodeID = prevIndex;

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


        superarcDependentWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcDependentWeightCoeffPortal = superarcDependentWeightCoeff.WritePortal();

        supernodeTransferWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto supernodeTransferWeightCoeffPortal = supernodeTransferWeightCoeff.WritePortal();

        // TODO: Assumption - for now, treat hyperarcs same as superarcs
        // NOTE: 2025-03-03 ran into the above assumption that breaks the code
//        hyperarcDependentWeightCoeff.Allocate(contourTree.Hyperarcs.GetNumberOfValues());
        hyperarcDependentWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto hyperarcDependentWeightCoeffPortal = hyperarcDependentWeightCoeff.WritePortal();

        Coefficients SAlocalIntrinsic;
        Coefficients SNlocalTransfer;
        Coefficients SAlocalDependent;
        Coefficients HAlocalDependent;

        // initialise the coefficient arrays:
        // Initialise the super{node/arc} weight arrays with 0s
        // (data below depends on the super{nodes/arcs}
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Supernodes.GetNumberOfValues(); sortedNode++)
        { // per node in sorted order
            SAlocalIntrinsic.h1 = 0.0l;
            SAlocalIntrinsic.h2 = 0.0l;
            SAlocalIntrinsic.h3 = 0.0l;
            SAlocalIntrinsic.h4 = 0.0l;

            SNlocalTransfer.h1 = 0.0l;
            SNlocalTransfer.h2 = 0.0l;
            SNlocalTransfer.h3 = 0.0l;
            SNlocalTransfer.h4 = 0.0l;

            SAlocalDependent.h1 = 0.0l;
            SAlocalDependent.h2 = 0.0l;
            SAlocalDependent.h3 = 0.0l;
            SAlocalDependent.h4 = 0.0l;

            superarcIntrinsicWeightCoeffPortal.Set(sortedNode, SAlocalIntrinsic);
            supernodeTransferWeightCoeffPortal.Set(sortedNode, SNlocalTransfer);
            superarcDependentWeightCoeffPortal.Set(sortedNode, SAlocalDependent);
        }//per node in sorted order

        // Initialise the hyper{node/arc} weight arrays with 0s
        // (data below depends on the hyper{nodes/arcs}
//        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Hypernodes.GetNumberOfValues(); sortedNode++)
        // TODO: Assumption - for now, treat hyperarcs same as superarcs
        // NOTE: 2025-03-03 ran into the above assumption that breaks the code
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Supernodes.GetNumberOfValues(); sortedNode++)
        {// per node in sorted order
            HAlocalDependent.h1 = 0.0l;
            HAlocalDependent.h2 = 0.0l;
            HAlocalDependent.h3 = 0.0l;
            HAlocalDependent.h4 = 0.0l;

            hyperarcDependentWeightCoeffPortal.Set(sortedNode, HAlocalDependent);
        }// per node in sorted order


        // TODO: Assumption - for now, treat hyperarcs same as superarcs
        // NOTE: 2025-03-03 ran into the above assumption that breaks the code

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {// per node in sorted order, add its weight to its superparent
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);

            SAlocalIntrinsic.h1 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h1 + vx_delta_h1_sum[sortID];
            SAlocalIntrinsic.h2 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h2 + vx_delta_h2_sum[sortID];
            SAlocalIntrinsic.h3 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h3 + vx_delta_h3_sum[sortID];
            SAlocalIntrinsic.h4 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h4 + vx_delta_h4_sum[sortID];

            superarcIntrinsicWeightCoeffPortal.Set(superparent, SAlocalIntrinsic);

        }//per node in sorted order, added its weight to its superparent

        std::cout << "    " << RED << std::setw(38) << std::left << "TOTAL SWEEP TIME"
                      << ": " << sweep_total_timer.GetElapsedTime() << " seconds" << RESET << std::endl;


std::cout << "// ================================= ITERATIONS =================================== //" << std::endl;

        vtkm::cont::Timer iterations_total_timer;
        iterations_total_timer.Start();

        // now iterate, propagating weights inwards
        // try to run for one more iteration to capture the whole tree
        std::cout << "Iteration: ";
        for (vtkm::Id iteration = 0; iteration < nIterations+1; iteration++)
        {// per iteration
          std::cout << iteration << " ";

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

          Coefficients step1Dependent;
          Coefficients step2Dependent;
          Coefficients step3Dependent;
          Coefficients step4HyperarcDependent;
          Coefficients step4SupernodeTransfer;

          std::map<vtkm::Id, vtkm::Id> tailends;

          for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues(); supernode++)
          {
              vtkm::Id superNode = supernodesPortal.Get(supernode);
              tailends.insert(std::make_pair(superNode, supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode)))));
          }

          // so, step 1: add xfer + int & store in dependent weight
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
            vtkm::Id superNode = supernodesPortal.Get(supernode);
            superarcDependentWeightPortal.Set(supernode,
                                              supernodeTransferWeightPortal.Get(supernode) +
                                                superarcIntrinsicWeightPortal.Get(supernode));
          }

          // so, step 1: add xfer + int & store in dependent weight
          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {
            step1Dependent.h1 = supernodeTransferWeightCoeffPortal.Get(supernode).h1 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h1;
            step1Dependent.h2 = supernodeTransferWeightCoeffPortal.Get(supernode).h2 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h2;
            step1Dependent.h3 = supernodeTransferWeightCoeffPortal.Get(supernode).h3 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h3;
            step1Dependent.h4 = supernodeTransferWeightCoeffPortal.Get(supernode).h4 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h4;

            superarcDependentWeightCoeffPortal.Set(supernode, step1Dependent);
          }

          //          std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
          // step 2: perform prefix sum on the dependent weight range
          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
          {
            superarcDependentWeightPortal.Set(supernode,
                                              superarcDependentWeightPortal.Get(supernode) +
                                                superarcDependentWeightPortal.Get(supernode - 1));
            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";
          }

          // step 2: perform prefix sum on the dependent weight range
          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
          {
              step2Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 + superarcDependentWeightCoeffPortal.Get(supernode-1).h1;
              step2Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 + superarcDependentWeightCoeffPortal.Get(supernode-1).h2;
              step2Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 + superarcDependentWeightCoeffPortal.Get(supernode-1).h3;
              step2Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 + superarcDependentWeightCoeffPortal.Get(supernode-1).h4;

              superarcDependentWeightCoeffPortal.Set(supernode, step2Dependent);
          }




// !!! CONVERTING FROM COEFFICIENT BASED TO SIMPLE(VALUE-TYPE) FOR DEPENDENT WEIGHTS !!! //

for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
{
    vtkm::Id superNode = supernodesPortal.Get(supernode);

    long double a_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h1 * std::pow(tailends[superNode], 3);
    long double b_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h2 * std::pow(tailends[superNode], 2);
    long double c_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h3 * tailends[superNode];
    long double d_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h4;

    // NEW: 2025-01-24 actually save the value to the dependent array SADWP:
    superarcDependentWeightPortal.Set(supernode, a_coeff + b_coeff + c_coeff + d_coeff);

}






          // step 3: subtract out the dependent weight of the prefix to the entire hyperarc. This will be a transfer, but for now, it's easier
          // to show it in serial. NB: Loops backwards so that computation uses the correct value
          // As a bonus, note that we test > firstsupernode, not >=.  This is because we've got unsigned integers, & otherwise it will not terminate
          // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
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

          } // per supernode

          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--)
          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
              continue;

            step3Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h1;
            step3Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h2;
            step3Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h3;
            step3Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h4;

            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
            superarcDependentWeightCoeffPortal.Set(supernode, step3Dependent);

          } // per supernode


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

          } // per hypernode

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

            step4HyperarcDependent.h1 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h1;
            step4HyperarcDependent.h2 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h2;
            step4HyperarcDependent.h3 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h3;
            step4HyperarcDependent.h4 = superarcDependentWeightCoeffPortal.Get(lastSuperarc).h4;

            // now, given the last superarc for the hyperarc, transfer the dependent weight
            hyperarcDependentWeightCoeffPortal.Set(hypernode, step4HyperarcDependent);

            // note that in parallel, this will have to be split out as a sort & partial sum in another array
            vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hypernode));

            step4SupernodeTransfer.h1 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h1 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h1;
            step4SupernodeTransfer.h2 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h2 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h2;
            step4SupernodeTransfer.h3 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h3 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h3;
            step4SupernodeTransfer.h4 = supernodeTransferWeightCoeffPortal.Get(hyperarcTarget).h4 + hyperarcDependentWeightCoeffPortal.Get(hypernode).h4;

            // now, given the last superarc for the hyperarc, transfer the dependent weight
//            hyperarcDependentWeightCoeffPortal.Set(hypernode, step4HyperarcDependent);

            supernodeTransferWeightCoeffPortal.Set(hyperarcTarget, step4SupernodeTransfer);

          } // per hypernode


          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
          {

            vtkm::Id superNode = supernodesPortal.Get(supernode);

            long double a_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h1 * std::pow(tailends[superNode], 3);
            long double b_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h2 * std::pow(tailends[superNode], 2);
            long double c_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h3 * tailends[superNode];
            long double d_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h4;

            // NEW: 2025-01-24 actually save the value to the dependent array SADWP:
            superarcDependentWeightPortal.Set(supernode, a_coeff + b_coeff + c_coeff + d_coeff);

            a_coeff = supernodeTransferWeightCoeffPortal.Get(supernode).h1 * std::pow(tailends[superNode], 3);
            b_coeff = supernodeTransferWeightCoeffPortal.Get(supernode).h2 * std::pow(tailends[superNode], 2);
            c_coeff = supernodeTransferWeightCoeffPortal.Get(supernode).h3 * tailends[superNode];
            d_coeff = supernodeTransferWeightCoeffPortal.Get(supernode).h4;

            // NEW: 2025-03-08 actually save the value to the dependent array SATWP:
            supernodeTransferWeightPortal.Set(supernode, a_coeff + b_coeff + c_coeff + d_coeff);
          }


        }// per iteration

        std::cout << std::endl;

std::cout << std::endl << "REINITILISING TRANSFER" << std::endl;
for(int i = 0; i < supernodeTransferWeightPortal.GetNumberOfValues(); i++)
{
    supernodeTransferWeightPortal.Set(i, 0.0l);
}

vtkm::Id previousParent = -1;

// 2025-03-09 RECOMPUTE THE INTRINSIC (FLOATING POINT) FROM THE CORRECT DEPENDENT AND TRANSFER CONNECTIONS:
for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
{
    vtkm::Id sortID = nodesPortal.Get(sortedNode);
    vtkm::Id superparent = superparentsPortal.Get(sortID);
    vtkm::Id hyperparent = hyperparentsPortal.Get(superparent); // uncommented 2025-03-03

    vtkm::Id nextSortID = (sortedNode+1 == contourTree.Arcs.GetNumberOfValues()) ? 0 : nodesPortal.Get(sortedNode+1);
    vtkm::Id nextSuperparent = superparentsPortal.Get(nextSortID);
//            vtkm::Id nextHyperparent = hyperparentsPortal.Get(nextSuperparent); // uncommented 2025-03-03

    vtkm::Id supertarget = MaskedIndex(superarcsPortal.Get(MaskedIndex(superparent))); //vtkm::cont::ArrayGetValue(superparent, contourTree.Superarcs);
    vtkm::Id hyperarcTarget = MaskedIndex(hyperarcsPortal.Get(hyperparent));

    if(superparent != previousParent )
    {
        if(supertarget)
        {
            supernodeTransferWeightPortal.Set(supertarget,
                                              supernodeTransferWeightPortal.Get(supertarget)+superarcDependentWeightPortal.Get(superparent));
        }

        previousParent = superparent;
    }


}

// 2025-03-09 RECOMPUTE THE INTRINSIC (FLOATING POINT) FROM THE CORRECT DEPENDENT AND TRANSFER CONNECTIONS:
for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
{
    superarcIntrinsicWeightPortal.Set(i, superarcDependentWeightPortal.Get(i) - supernodeTransferWeightPortal.Get(i));
}


        std::cout << "END ComputeVolumeWeightsSerialStructCoefficients" << std::endl;

#if PROFILING_PACTBD
            printMemoryUsage("[ProcessContourTree.h::ComputeVolumeWeightsSerialStructCoefficients] Checkpoint 4/4 - END");

            std::cout << "    " << RED << std::setw(38) << std::left << "Iterations Total Time"
                          << ": " << iterations_total_timer.GetElapsedTime() << " seconds" << RESET << std::endl;
#endif

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
      std::cout << "[ProcessContourTree.h::ComputeVolumeBranchDecompositionSerialFloat()] START" << std::endl;

      // COMMS: Both 'intrinsic' and 'dependent' weights come precomputed from 'ComputeVolumeWeights()'
      // ... the following just sets up the read portals for both
      auto superarcDependentWeightPortal = superarcDependentWeight.ReadPortal();
      auto superarcIntrinsicWeightPortal = superarcIntrinsicWeight.ReadPortal();

      // CHANGE: because at the moment the dependent and intrinsic weights have some errors ...
      // ... I write the arrays with the correct values for each to see how the BT should work.
//      FloatArrayType superarcDependentWeightCorrect= superarcDependentWeight.ReadPortal();
//      FloatArrayType superarcIntrinsicWeightCorrect= superarcIntrinsicWeight.ReadPortal();

//      superarcDependentWeightCorrect.Allocate(superarcDependentWeightPortal.GetNumberOfValues());
//      superarcIntrinsicWeightCorrect.Allocate(superarcIntrinsicWeightPortal.GetNumberOfValues());


      std::string indent = "\t";
//      std::array<double, 6> realIntrinsic = {0.0208333, 0.14127, 0.178175, 0.0236112,
//                                             0.636111,                   0.0};
//      std::array<double, 6> realDependent = {0.0208333, 0.14127, 0.178175, 0.0236112,  0.636111+0.0208333+0.14127, 1.0};

      auto superarcDependentWeightCorrectReadPortal = superarcDependentWeight.ReadPortal();
      auto superarcIntrinsicWeightCorrectReadPortal = superarcIntrinsicWeight.ReadPortal();

      auto superarcDependentWeightCorrectWritePortal = superarcDependentWeight.WritePortal();
      auto superarcIntrinsicWeightCorrectWritePortal = superarcIntrinsicWeight.WritePortal();

#if DEBUG_PRINT_PACTBD
      std::cout << std::endl << "(ComputeVolumeBranchDecompositionSerialFloat) Superarc Intrinsic Weight Portal (vs Correct):" << std::endl;
      for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << indent << i << " -> " << superarcIntrinsicWeightCorrectReadPortal.Get(i) << std::endl;
      }
      std::cout << std::endl;

      std::cout << std::endl << "(ComputeVolumeBranchDecompositionSerialFloat) superarc Dependent Weight Portal:" << std::endl;
      for(int i = 0; i < superarcDependentWeightPortal.GetNumberOfValues(); i++)
      {
          std::cout << indent << i << " -> " << superarcDependentWeightCorrectReadPortal.Get(i) << std::endl;
      }
      std::cout << std::endl;
#endif



      // cache the number of non-root supernodes & superarcs
      vtkm::Id nSupernodes = contourTree.Supernodes.GetNumberOfValues();
      vtkm::Id nSuperarcs = nSupernodes - 1;

      // STAGE I:  Find the upward and downwards weight for each superarc, and set up arrays

      // Allocation/initialization
      // COMMS: just allocate up-weight arrays, do not set values yet
//      IdArrayType upWeight;
//      upWeight.Allocate(nSuperarcs);
//      auto upWeightPortal = upWeight.WritePortal();

          // FLOAT versions
          FloatArrayType upWeightFloat;
          upWeightFloat.Allocate(nSuperarcs);
          auto upWeightFloatPortal = upWeightFloat.WritePortal();

          // FLOAT CORRECT versions
          FloatArrayType upWeightFloatCorrect;
          upWeightFloatCorrect.Allocate(nSuperarcs);
          auto upWeightFloatCorrectPortal = upWeightFloatCorrect.WritePortal();

      // COMMS: just allocate down-weight arrays, do not set values yet
//      IdArrayType downWeight;
//      downWeight.Allocate(nSuperarcs);
//      auto downWeightPortal = downWeight.WritePortal();

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

      totalVolumeFloat = superarcDependentWeightPortal.Get(superarcDependentWeightPortal.GetNumberOfValues()-1);

      std::cout << std::setw(55) << "Total Volume Int (Nodes)" << "= " << totalVolume << std::endl;
      std::cout << std::setw(55) << "Total Volume Float" << "= " << totalVolumeFloat << std::endl;

#if DEBUG_PRINT_PACTBD
      std::cout << "\n" << "----------------- Computing Up/Down Weights -----------------" << std::endl;

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
#if DEBUG_PRINT_PACTBD
        std::cout << "Processing superarc: " << superarc << std::endl << "{" << std::endl;
#endif
        if (IsAscending(superarcsPortal.Get(superarc))) // flag on the ID of superarc
        { // ascending superarc
#if DEBUG_PRINT_PACTBD
          std::cout << indent << "ASCENDING\n";
#endif
            // put the lower-end first
          superarcListWritePortal.Set(superarc, // each superarc starts at the supernode at the same ID and goes to the supernode whose ID is stored in the superarc's array
                                      // pair the origin and the destination of that superarc. We store them in an edge pair and write it to the array
                                      EdgePair(superarc, MaskedIndex(superarcsPortal.Get(superarc))));

          vtkm::Id superNode = supernodesPortal.Get(superarc);
#if DEBUG_PRINT_PACTBD
          std::cout << indent << superarc << " = " << superNode << " -> " << MaskedIndex(superarcsPortal.Get(superarc)) << std::endl << std::endl;
#endif
          // because this is an ascending superarc, the dependent weight refers to the weight at the upper end
          // so, we set the up weight based on the dependent weight
//          upWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          upWeightFloatPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          upWeightFloatCorrectPortal.Set(superarc, superarcDependentWeightCorrectReadPortal.Get(superarc));


          // at the inner end, dependent weight is the total in the subtree.
          // Then there are vertices along the edge itself (intrinsic weight), ...
          // ... including the supernode at the outer end
          // So, to get the "dependent" weight in the other direction, ...
          // ... we start with totalVolume - dependent, then subtract (intrinsic - 1)
          // set the weight at the down end by using the invert operator:
//          downWeightPortal.Set(superarc,
//                               // below is the invert operator for node count!
//                               (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
//                                 (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          downWeightFloatPortal.Set(superarc,
                                   // below is the invert operator for node count!
                                   (totalVolumeFloat - superarcDependentWeightPortal.Get(superarc)) +
                                     (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          downWeightFloatCorrectPortal.Set(superarc,
                                   // below is the invert operator for node count!
                                   (totalVolumeFloat - superarcDependentWeightCorrectReadPortal.Get(superarc)) +
                                     (superarcIntrinsicWeightCorrectReadPortal.Get(superarc) - 1));
#if DEBUG_PRINT_PACTBD
//          std::cout << indent << "upWeightPortal             = " << upWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatPortal        = " << upWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatCorrectPortal = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
          std::cout << std::endl;
//          std::cout << indent << "downWeightPortal             = " << downWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatPortal        = " << downWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatCorrectPortal = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;

          std::cout << "\n}\n";
#endif
        } // ascending superarc
        else
        { // descending superarc
#if DEBUG_PRINT_PACTBD
          std::cout << indent << "DESCENDING\n";
#endif
          // lower-end is also first, but in the reverse order compared to IsAscending
          superarcListWritePortal.Set(superarc,
                                      EdgePair(MaskedIndex(superarcsPortal.Get(superarc)), superarc));

          vtkm::Id superNode = supernodesPortal.Get(superarc);
#if DEBUG_PRINT_PACTBD
          std::cout << indent << superarc << " = " << superNode << " -> " << MaskedIndex(superarcsPortal.Get(superarc)) << std::endl << std::endl;
#endif

//          downWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          downWeightFloatPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
          downWeightFloatCorrectPortal.Set(superarc, superarcDependentWeightCorrectReadPortal.Get(superarc));

          // at the inner end, dependent weight is the total in the subtree.
          // Then there are vertices along the edge itself (intrinsic weight), ...
          // ... including the supernode at the outer end
          // So, to get the "dependent" weight in the other direction, ...
          // ... we start with totalVolume - dependent, then subtract (intrinsic - 1)
//          upWeightPortal.Set(superarc,
//                             (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
//                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          upWeightFloatPortal.Set(superarc,
                             (totalVolumeFloat - superarcDependentWeightPortal.Get(superarc)) +
                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));

          upWeightFloatCorrectPortal.Set(superarc,
                                   // below is the invert operator for node count!
                                   (totalVolumeFloat - superarcDependentWeightCorrectReadPortal.Get(superarc)) +
                                     (superarcIntrinsicWeightCorrectReadPortal.Get(superarc) - 1));
#if DEBUG_PRINT_PACTBD
//          std::cout << indent << "upWeightPortal      = " << upWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatPortal = " << upWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "upWeightFloatCorrectPortal = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
          std::cout << std::endl;
//          std::cout << indent << "downWeightPortal      = " << downWeightPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatPortal = " << downWeightFloatPortal.Get(superarc) << std::endl;
          std::cout << indent << "downWeightFloatCorrectPortal = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;

          std::cout << "\n}\n";
#endif

        } // descending superarc
      }   // per superarc

#if DEBUG_PRINT_PACTBD
    std::cout << "Up Weights (float)" << std::endl;
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    {
        vtkm::Id superNode = supernodesPortal.Get(superarc);

        std::cout << superarc << "(" << superNode  << ") = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
    }
    std::cout << "Down Weights (float)" << std::endl;
    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    {
        vtkm::Id superNode = supernodesPortal.Get(superarc);

        std::cout << superarc << "(" << superNode  << ") = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;
    }

    std::cout << std::endl;
//    std::cout << "Up Weights (int)" << std::endl;
//    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
//    {
//        vtkm::Id superNode = supernodesPortal.Get(superarc);

//        std::cout << superarc << "(" << superNode  << ") = " << upWeightPortal.Get(superarc) << std::endl;
//    }
//    std::cout << "Down Weights (int)" << std::endl;
//    for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
//    {
//        vtkm::Id superNode = supernodesPortal.Get(superarc);

//        std::cout << superarc << "(" << superNode  << ") = " << downWeightPortal.Get(superarc) << std::endl;
//    }
#endif

#if DEBUG_PRINT_PACTBD
      std::cout << "II A. Weights Computed" << std::endl;
      PrintHeader(upWeightFloatCorrect.GetNumberOfValues());
      PrintIndices("Intrinsic Weight", superarcIntrinsicWeight);
      PrintIndices("Dependent Weight", superarcDependentWeight);
//      PrintIndices("Upwards Weight",           upWeight);
      PrintValues("Upwards Weight (float)",   upWeightFloatCorrect);
//      PrintIndices("Downwards Weight",         downWeight);
      PrintValues("Downwards Weight (float)", downWeightFloatCorrect);
      std::cout << std::endl;
#endif

      // II B. Pick the best downwards weight by sorting on upper vertex then processing by segments
      // II B 1.      Sort the superarcs by upper vertex
      IdArrayType superarcSorter;
      superarcSorter.Allocate(nSuperarcs);
      auto superarcSorterPortal = superarcSorter.WritePortal();
      // make the array of indices for indirect sorting

#if DEBUG_PRINT_PACTBD
      std::cout << "Unsorted arcs:" << std::endl;
#endif
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
      {
        superarcSorterPortal.Set(superarc, superarc);
#if DEBUG_PRINT_PACTBD
        std::cout << superarc << ") " << superarcSorterPortal.Get(superarc) << " - " << upWeightFloatCorrectPortal.Get(superarcSorterPortal.Get(superarc)) << " = " << upWeightFloatCorrectPortal.Get(superarcSorterPortal.Get(superarc))  << std::endl;
#endif
      }

//      // OLD: Vertex count sort
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

#if DEBUG_PRINT_PACTBD
      std::cout << "Sorted arcs (by float upweight):" << std::endl;
      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
      {
        std::cout << superarc << ") " << superarcSorterPortal.Get(superarc) << " - " << upWeightFloatCorrectPortal.Get(superarcSorterPortal.Get(superarc)) << " = " << upWeightFloatCorrectPortal.Get(superarcSorterPortal.Get(superarc))  << std::endl;
      }
#endif


//      std::cout << "Sorted arcs (by int supernode count):" << std::endl;
//      for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
//      {
//        std::cout << superarc << ") " << superarcSorterPortal.Get(superarc) << " - " << upWeightPortal.Get(superarcSorterPortal.Get(superarc)) << " = " << upWeightFloatCorrectPortal.Get(superarcSorterPortal.Get(superarc))  << std::endl;
//      }


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

//      // OLD: Vertex count sort
//      vtkm::cont::Algorithm::Sort(
//        superarcSorter,
//        process_contourtree_inc_ns::SuperArcVolumetricComparator(downWeight, superarcList, true));


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

#if DEBUG_PRINT_PACTBD
      std::cout << "II. Best Edges Selected" << std::endl;
      PrintHeader(bestUpward.GetNumberOfValues());
      PrintIndices("Best Upwards", bestUpward);
      PrintIndices("Best Downwards", bestDownward);
      std::cout << std::endl;
#endif

      ProcessContourTree::ComputeBranchData(contourTree,
                                            whichBranch,    // (output)
                                            branchMinimum,  // (output)
                                            branchMaximum,  // (output)
                                            branchSaddle,   // (output)
                                            branchParent,   // (output)
                              /* (input) */ bestUpward,
                              /* (input) */ bestDownward);

      std::cout << std::setw(55) << "Num. of Branches in the Branch Decomposition" << "= " << branchSaddle.GetNumberOfValues() << std::endl;


      std::cout << "[ProcessContourTree.h::ComputeVolumeBranchDecompositionSerialFloat()] END" << std::endl;

    } // ComputeVolumeBranchDecompositionSerialFloat()














  // routine to compute the branch decomposition by volume (the original)
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


    std::cout << "Total Volume: " << totalVolume << std::endl;

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

#if DEBUG_PRINT_PACTBD
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

#if DEBUG_PRINT_PACTBD
    std::cout << "II. Best Edges Selected" << std::endl;
    PrintHeader(bestUpward.GetNumberOfValues());
    PrintIndices("Best Upwards", bestUpward);
    PrintIndices("Best Downwards", bestDownward);
    std::cout << std::endl;
#endif

    ProcessContourTree::ComputeBranchData(contourTree,
                                          whichBranch,      // (output)
                                          branchMinimum,    // (output)
                                          branchMaximum,    // (output)
                                          branchSaddle,     // (output)
                                          branchParent,     // (output)
                          /* input */     bestUpward,
                          /* input */     bestDownward);

  } // ComputeVolumeBranchDecomposition()

  // routine to compute the branch decomposition by volume
  void static ComputeBranchData(const ContourTree& contourTree,
                                IdArrayType& whichBranch,       // (output)
                                IdArrayType& branchMinimum,     // (output)
                                IdArrayType& branchMaximum,     // (output)
                                IdArrayType& branchSaddle,      // (output)
                                IdArrayType& branchParent,      // (output)
                /* (input) */   IdArrayType& bestUpward,
                /* (input) */   IdArrayType& bestDownward)
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
    // Let v = BestUp(u).
    // Then:
    //      if u = BestDown(v), ...
    //      ... copy v (BestUp(u) = v) to whichBranch[u]
    // Otherwise, ...
    //      let whichBranch[u] = u | TERMINAL to mark the end of the side branch
    //      ( | TERMINAL is a bitmask which looks as '.t...' when printed in terminal)
    //      [this comment used to be incorrect - it had whichBranch[u] = BestUp(u)] ...
    //
    // NB 1: Leaves already have the flag set, but it's redundant so its safe
    // NB 2: We don't need to do it downwards because it's symmetric
    vtkm::cont::Invoker invoke;
    vtkm::worklet::contourtree_augmented::process_contourtree_inc::PropagateBestUpDown
      propagateBestUpDownWorklet;
    invoke(propagateBestUpDownWorklet, bestUpward, bestDownward, whichBranch);

#if DEBUG_PRINT_PACTBD
    std::cout <<  std::endl << std::endl << "III. Branch Neighbours Identified" << std::endl;
    PrintHeader(whichBranch.GetNumberOfValues());
    PrintIndices("Which Branch", whichBranch);
    std::cout << std::endl;
#endif

    // -------------------------------------------------- POINTER DOUBLING START --------------------------------------------------- //

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


    // --------------------------------------------------- POINTER DOUBLING END --------------------------------------------------- //

#if DEBUG_PRINT_PACTBD
    std::cout <<  std::endl << std::endl << "IV. Branch Chains Propagated" << std::endl;
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

#if DEBUG_PRINT_PACTBD
    std::cout <<  std::endl << std::endl << "V. Branch Arrays Created" << std::endl;
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

#if DEBUG_PRINT_PACTBD
    std::cout <<  std::endl << std::endl << "VI A. Sorted into Branches" << std::endl;
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



//    std::cout << "VI. Branches Set" << std::endl;
//    PrintHeader(nBranches);
//    PrintIndices("Branch Maximum", branchMaximum);
//    PrintIndices("Branch Minimum", branchMinimum);
//    PrintIndices("Branch Saddle", branchSaddle);
//    PrintIndices("Branch Parent", branchParent);


#if DEBUG_PRINT_PACTBD
    std::cout <<  std::endl << std::endl << "VI. Branches Set" << std::endl;
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

#if DEBUG_PRINT_PACTBD
    std::cout <<  std::endl << std::endl << "VII. Branches Constructed" << std::endl;
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



  // Create branch decomposition from contour tree
  template <typename T, typename StorageType>
  static process_contourtree_inc_ns::Branch<T>* ComputeBranchDecomposition(
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
    const FloatArrayType& superarcDependentWeight,            // NEW: passed intrincid
    const FloatArrayType& superarcIntrinsicWeight)
  {
    std::cout << "ContourTreeApp->(ProcessContourTree)->Branch.h->ComputeBranchDecomposition()" << std::endl;

    return process_contourtree_inc_ns::Branch<T>::ComputeBranchDecomposition(
      contourTreeSuperparents,
      contourTreeSupernodes,
                contourTreeSuperarcs,
      whichBranch,
      branchMinimum,
      branchMaximum,
      branchSaddle,
      branchParent,
      sortOrder,
      dataField,
      dataFieldIsSorted,
      superarcDependentWeight,
      superarcIntrinsicWeight);
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
