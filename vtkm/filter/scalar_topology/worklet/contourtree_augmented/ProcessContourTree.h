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


// For 2025-10 2025 October Betti Insertion:
// wrong #include <vtkm/filter/scalar_topology/worklet/contourtree_distributed/hierarchical_augmenter/ResizeArraysBuildNewSupernodeIdsWorklet.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h> // actually defined the ResizeVector

// for sleeping
#include <chrono>
#include <thread>

// for memory usage
#include <sys/resource.h>
#include <unistd.h>

#include <float.h>

// file IO
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#define DEBUG_PRINT_PACTBD 1
#define SLEEP_ON 0
#define PROFILING_PACTBD 1
#define WRITE_FILES 0
#define ENABLE_DEBUG_MINMAX 1

#define UPDATE_MINMAX(val, min_var, max_var) \
    do {                                     \
        if ((val) < (min_var)) (min_var) = (val); \
        if ((val) > (max_var)) (max_var) = (val); \
    } while(0)

#define UPDATE_MINMAX_ABS(val, min_var, max_var)     \
    do {                                             \
        long double abs_val = fabsl(val);            \
        if (abs_val < (min_var)) (min_var) = abs_val; \
        if (abs_val > (max_var)) (max_var) = abs_val; \
    } while(0)

#define UPDATE_MINMAX_ABS_NONZERO(val, min_var, max_var)       \
    do {                                                       \
        long double abs_val = fabsl(val);                      \
        if (abs_val > 0.0L) {                                  \
            if (abs_val < (min_var)) (min_var) = abs_val;      \
            if (abs_val > (max_var)) (max_var) = abs_val;      \
        }                                                      \
    } while(0)

#ifdef ENABLE_DEBUG_MINMAX
#define TRACK_MINMAX(val, minv, maxv) UPDATE_MINMAX_ABS_NONZERO(val, minv, maxv) //UPDATE_MINMAX(val, minv, maxv)
#else
#define TRACK_MINMAX(val, minv, maxv) // nothing
#endif


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


struct BettiCoefficients
{
    long num_vtx;
    long num_edg;
    long num_fac;
    long num_tet;

    long betti0;
    long betti1;
    long betti3;
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

//        std::cout << ORANGE << message << LIGHT_BLUE << " - Memory usage (peak): " << usage.ru_maxrss
//                  << " KB | (current) " << current_usage << " KB" << RESET << std::endl;

        std::cout << LIGHT_BLUE << message << " - Memory usage (peak): " << usage.ru_maxrss
                  << " KB | (current) " << current_usage << " KB" << RESET << std::endl;

    }












    // Hash and equality for 3-vertex faces
    struct FaceHash {
        std::size_t operator()(const std::array<vtkm::Id, 3>& f) const noexcept {
            return std::hash<vtkm::Id>()(f[0]) ^ (std::hash<vtkm::Id>()(f[1]) << 1) ^ (std::hash<vtkm::Id>()(f[2]) << 2);
        }
    };
    struct FaceEq {
        bool operator()(const std::array<vtkm::Id, 3>& a, const std::array<vtkm::Id, 3>& b) const noexcept {
            return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
        }
    };

    struct TriangleFace {
        vtkm::Id v0, v1, v2;
        bool boundary;
    };

    std::vector<TriangleFace> static GetTrianglesWithBoundary(const vtkm::cont::DataSet& input,
                                                       bool onlyBoundary = false)
    {
        using TetCellSet = vtkm::cont::CellSetSingleType<>;

        const auto& unknown = input.GetCellSet();
        if (!unknown.IsType<TetCellSet>()) {
            throw std::runtime_error("Error: Dataset does not contain CellSetSingleType<> (tets).");
        }

        const auto& cellSet = unknown.AsCellSet<TetCellSet>();
        vtkm::Id numCells = cellSet.GetNumberOfCells();

        std::unordered_map<std::array<vtkm::Id, 3>, int, FaceHash, FaceEq> faceCounter;

        // Loop over all tets and count their faces
        for (vtkm::Id cellId = 0; cellId < numCells; ++cellId) {
            vtkm::Id ptIds[4];
            cellSet.GetCellPointIds(cellId, ptIds);

            // Each tetrahedron has 4 triangular faces
            std::array<std::array<vtkm::Id, 3>, 4> faces = {{
                {ptIds[0], ptIds[1], ptIds[2]},
                {ptIds[0], ptIds[1], ptIds[3]},
                {ptIds[0], ptIds[2], ptIds[3]},
                {ptIds[1], ptIds[2], ptIds[3]}
            }};

            for (auto& f : faces) {
                std::sort(f.begin(), f.end()); // canonical order
                faceCounter[f]++;
            }
        }

        // Build result
        std::vector<TriangleFace> result;
        result.reserve(faceCounter.size());

        for (const auto& [f, count] : faceCounter) {
            bool isBoundary = (count == 1);
            if (onlyBoundary && !isBoundary)
                continue;

            result.push_back({f[0], f[1], f[2], isBoundary});
        }

        // Sort by vertex indices for deterministic output
        std::sort(result.begin(), result.end(),
                  [](const TriangleFace& a, const TriangleFace& b) {
                      if (a.v0 != b.v0) return a.v0 < b.v0;
                      if (a.v1 != b.v1) return a.v1 < b.v1;
                      return a.v2 < b.v2;
                  });

        return result;
    }







    // Hash and equality for 2-vertex edges
    struct EdgeHash {
        std::size_t operator()(const std::array<vtkm::Id, 2>& e) const noexcept {
            return std::hash<vtkm::Id>()(e[0]) ^ (std::hash<vtkm::Id>()(e[1]) << 1);
        }
    };
    struct EdgeEq {
        bool operator()(const std::array<vtkm::Id, 2>& a, const std::array<vtkm::Id, 2>& b) const noexcept {
            return a[0] == b[0] && a[1] == b[1];
        }
    };

    // A simple edge struct for clarity
    struct Edge {
        vtkm::Id v0, v1;
    };

    std::vector<Edge> static GetEdgesFromVTK(const vtkm::cont::DataSet& input)
    {
        using TetCellSet = vtkm::cont::CellSetSingleType<>;

        const auto& unknown = input.GetCellSet();
        if (!unknown.IsType<TetCellSet>()) {
            throw std::runtime_error("Error: Dataset does not contain CellSetSingleType<> (expected tetrahedra).");
        }

        const auto& cellSet = unknown.AsCellSet<TetCellSet>();
        vtkm::Id numCells = cellSet.GetNumberOfCells();

        std::unordered_set<std::array<vtkm::Id, 2>, EdgeHash, EdgeEq> edgeSet;

        for (vtkm::Id cellId = 0; cellId < numCells; ++cellId) {
            vtkm::Id ptIds[4];
            cellSet.GetCellPointIds(cellId, ptIds);

            // Ensure each tetrahedron's vertex order is canonical
            std::sort(ptIds, ptIds + 4);

            // Add all 6 unique edges per tetrahedron
            std::array<std::array<vtkm::Id, 2>, 6> edges = {{
                {ptIds[0], ptIds[1]},
                {ptIds[0], ptIds[2]},
                {ptIds[0], ptIds[3]},
                {ptIds[1], ptIds[2]},
                {ptIds[1], ptIds[3]},
                {ptIds[2], ptIds[3]}
            }};

            for (auto& e : edges) {
                std::sort(e.begin(), e.end());
                edgeSet.insert(e);
            }
        }

        // Move to a sorted vector for deterministic output
        std::vector<Edge> result;
        result.reserve(edgeSet.size());
        for (const auto& e : edgeSet)
            result.push_back({e[0], e[1]});

        std::sort(result.begin(), result.end(),
                  [](const Edge& a, const Edge& b) {
                      return (a.v0 < b.v0) || (a.v0 == b.v0 && a.v1 < b.v1);
                  });

        return result;
    }




    void static LUstars(// INPUTS
                        int numVertices, // in sort order, enough to have the total number, since we start from 0 incrementing by 1 up to N
                        std::vector<Edge>&              edges,
                        std::vector<TriangleFace>&      triangles,
                        std::vector<std::vector<int>>&  tetrahedra,
                        // OUTPUTS
                        std::vector<int>& lowerStars,
                        std::vector<int>& upperStars,
                        std::vector<int>& deltaBoundary)
    {
        // 1. initialise LU, US, dB:
        lowerStars.resize(numVertices, 1);
        upperStars.resize(numVertices, 1);
        deltaBoundary.resize(numVertices, 0);
        std::cout << "LU\tUS\tdB:" << std::endl;
        for(int i = 0; i < numVertices; i++)
        {
            std::cout << i << "\t" << lowerStars[i] << "\t" << upperStars[i] << "\t" << deltaBoundary[i] << std::endl;
        }

        // 2. for each edge:
        int i,j;
        for(int it = 0; it < edges.size(); it++)
        {
            i = edges[it].v0;
            j = edges[it].v1;
            if(i < j)
            {
                lowerStars[j]--;
                upperStars[i]--;
            }
        }

        std::cout << "(Edge)LU\tUS\tdB:" << std::endl;
        for(int i = 0; i < numVertices; i++)
        {
            std::cout << i << "\t" << lowerStars[i] << "\t" << upperStars[i] << "\t" << deltaBoundary[i] << std::endl;
        }

        // 3. for each edge:
        i=0;
        j=0;
        int k;
        bool b;
        for(int it = 0; it < triangles.size(); it++)
        {
            i = triangles[it].v0;
            j = triangles[it].v1;
            k = triangles[it].v2;
            b = triangles[it].boundary;
            if ((i < j) && (j < k))
            {
                lowerStars[k]++;
                upperStars[i]++;
                if (b)
                {
                    deltaBoundary[k]--;
                    deltaBoundary[i]++;
                }
            }
        }

        std::cout << "(Triangle)LU\tUS\tdB:" << std::endl;
        for(int i = 0; i < numVertices; i++)
        {
            std::cout << i << "\t" << lowerStars[i] << "\t" << upperStars[i] << "\t" << deltaBoundary[i] << std::endl;
        }

        // 3. for each tetrahedron:
        i=0;
        j=0;
        k=0;
        int l;
        for(int it = 0; it < tetrahedra.size(); it++)
        {
            i = tetrahedra[it][0];
            j = tetrahedra[it][1];
            k = tetrahedra[it][2];
            l = tetrahedra[it][3];

            if ((i < j) && (j < k) && (k < l))
            {
                lowerStars[l]--;
                upperStars[i]--;
            }
        }

        std::cout << "(Tet)LU\tUS\tdB:" << std::endl;
        for(int i = 0; i < numVertices; i++)
        {
            std::cout << i << "\t" << lowerStars[i] << "\t" << upperStars[i] << "\t" << deltaBoundary[i] << std::endl;
        }

    }




   // 2025-10-11 COMPUTE BETTI NUMBERS FOR EACH REGULAR BRANCH
   // BASED ON 2004 Pascucci Parallel Computation of the Topology of Level Sets
    void static ComputeBettiNumbersForRegularArcs(const vtkm::cont::DataSet& input, // the coefficient-based version additionally requires tetrahedral connections and vertex coordinates
                                                  const ContourTree& contourTree,
                                                  const vtkm::Id nIterations,
                                                  vtkm::cont::ArrayHandle<Coefficients>& superarcIntrinsicWeightCoeff, // (output)
                                                  vtkm::cont::ArrayHandle<Coefficients>& superarcDependentWeightCoeff, // (output)
                                                  vtkm::cont::ArrayHandle<Coefficients>& supernodeTransferWeightCoeff, // (output)
                                                  vtkm::cont::ArrayHandle<Coefficients>& hyperarcDependentWeightCoeff, // (output)
                                                  // Added 2025-01-30
                                                  // We use simple weights for the branch decomposition
                                                  FloatArrayType& superarcIntrinsicWeight, // (output)
                                                  FloatArrayType& superarcDependentWeight, // (output)
                                                  FloatArrayType& supernodeTransferWeight, // (output)
                                                  FloatArrayType& hyperarcDependentWeight) // (output))
    {
        std::cout << "[ProcessContourTree.h::ComputeBettiNumbersForRegularArcs] Compute Betti Numbers for each Regular Arc" << std::endl;
        printMemoryUsage("[ProcessContourTree.h::ComputeVolumeWeightsSerialStructCoefficients] Checkpoint 1/4 - START");

        using TetCellSet = vtkm::cont::CellSetSingleType<>;
        const auto& unknown = input.GetCellSet();


        std::vector<int> lowerStars;
        std::vector<int> upperStars;
        std::vector<int> deltaBoundary;
        if (unknown.IsType<TetCellSet>())
        {
            const auto& cellSet = unknown.AsCellSet<TetCellSet>();
            vtkm::Id numCells = cellSet.GetNumberOfCells();

            std::vector<std::vector<int>> tetlistSorted(numCells,                 // size
                                                        std::vector<int> (4, 0)); // 4 ints (vertices) per tet initialised to 0

            std::cout << "Number of tets: " << numCells << std::endl;

            for (vtkm::Id cellId = 0; cellId < numCells; ++cellId)
            {
                vtkm::IdComponent npts = cellSet.GetNumberOfPointsInCell(cellId);
                vtkm::Id ptIds[4];
                cellSet.GetCellPointIds(cellId, ptIds);
                // process tetrahedron here
                std::cout << cellId << "\t" << ptIds[0] << "\t" << ptIds[1] << "\t" << ptIds[2] << "\t" << ptIds[3] << std::endl;
                tetlistSorted[cellId][0] = ptIds[0];
                tetlistSorted[cellId][1] = ptIds[1];
                tetlistSorted[cellId][2] = ptIds[2];
                tetlistSorted[cellId][3] = ptIds[3];
            }

            std::cout << "Tetrahedra | Expecting Last:\n383\t0\t1\t6\t124" << std::endl;

            //  Sort connected tetrahedral vertices in increasing order
            // (then we can assume the first vertex of the tet holds the lowest value)
            for (int i = 0; i < numCells; i++)
            {// for each tet
                std::sort(tetlistSorted[i].begin(), tetlistSorted[i].end());
            }// for each tet

            // Print the sorted tets for checking against a manually computed list:
            for (vtkm::Id i = 0; i < tetlistSorted.size(); ++i)
            {
                std::cout << i << "\t" << tetlistSorted[i][0]
                               << "\t" << tetlistSorted[i][1]
                               << "\t" << tetlistSorted[i][2]
                               << "\t" << tetlistSorted[i][3] << std::endl;
            }

            std::cout << "Tetrahedra | Expecting Last:\n383\t0\t1\t6\t124" << std::endl;

            auto triangles = GetTrianglesWithBoundary(input, /*onlyBoundary=*/false);

            int triangle_id = 0;
            for (const auto& tri : triangles)
            {
                std::cout << triangle_id
                          << "\t" << tri.v0 << "\t" << tri.v1 << "\t" << tri.v2
                          << "\t" << tri.boundary << "\n";
                triangle_id++;
            }

            std::cout << "Triangles | Expecting Last:\n863	112	116	119	0" << std::endl;

            auto edges = GetEdgesFromVTK(input);
            std::cout << "Number of unique edges: " << edges.size() << "\n";

            int edge_id = 0;
            for (const auto& e : edges) {
                std::cout << edge_id
                          << "\t" << e.v0 << "\t" << e.v1 << "\n";
            }

            std::cout << "Edges | Expecting Last:\n603	120	124" << std::endl;

            std::cout << "Running LUstars ..." << std::endl;
            LUstars(contourTree.Arcs.GetNumberOfValues(),
                    edges,
                    triangles,
                    tetlistSorted,
                    lowerStars,
                    upperStars,
                    deltaBoundary);

        }

        std::cout << "Betti Number Regular Instrinsic Pre-Processing:" << std::endl;

        auto arcsPortal = contourTree.Arcs.ReadPortal();
        auto nodesPortal = contourTree.Nodes.ReadPortal();
        auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto supernodesPortal = contourTree.Supernodes.ReadPortal();
        auto superparentsPortal = contourTree.Superparents.ReadPortal();

        std::vector<vtkm::Id> chi_xij, chi_x;
        std::vector<vtkm::Id> be_ij,   bei;
        chi_xij.resize( contourTree.Arcs.GetNumberOfValues(), 0);
        chi_x.resize(   contourTree.Arcs.GetNumberOfValues(), 0);
        bei.resize(     contourTree.Arcs.GetNumberOfValues(), 0);
        be_ij.resize(   contourTree.Arcs.GetNumberOfValues(), 0);

        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {// for each sortedNode
            vtkm::Id i_sortID = nodesPortal.Get(sortedNode);
            vtkm::Id i_superparent = superparentsPortal.Get(i_sortID);

            vtkm::Id j_sortID = nodesPortal.Get(sortedNode+1);
            vtkm::Id j_superparent = superparentsPortal.Get(j_sortID );

            vtkm::Id tailend = supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(i_superparent)));

            vtkm::Id delta = 1;

            if (i_superparent == j_superparent)
            {
                tailend = j_sortID;
            }

            std::cout << sortedNode << "\t" << i_sortID << "\t" << tailend << std::endl;

            if(i_sortID > tailend)
            {
                delta = -1;
            }

            chi_xij[sortedNode] = delta * (chi_x[i_sortID] - upperStars[i_sortID] + lowerStars[i_sortID]);
            be_ij[sortedNode]   = delta * (bei[i_sortID] + deltaBoundary[i_sortID]);

            chi_x[tailend] += delta * chi_xij[sortedNode];
            bei[tailend] += delta * be_ij[sortedNode];
        }

        std::vector<vtkm::Id> regular_nodes_to_insert;
        std::vector<vtkm::Id> node_ascend;

        std::vector<vtkm::Id> nodes_to_relabel_superparent; // 2026-01-03 addition
        std::vector<vtkm::Id> nodes_to_relabel_hyperparent;
        std::vector<double> nodes_to_relabel_dataflip;
        std::vector<vtkm::Id> nodes_to_relabel_regularID;
        int previous_betti1 = 0;


        std::cout << "VTK-m FIELDS:" << std::endl;
        // Loop over all fields
        const vtkm::cont::Field& field = input.GetPointField("var");

        // Get the UnknownArrayHandle
        vtkm::cont::UnknownArrayHandle ua = field.GetData();

        // Cast to the correct array type:
        // (replace float with the actual value type)
        vtkm::cont::ArrayHandle<double> array;
        ua.AsArrayHandle(array);

        // Get read-only access
        auto dataPortal = array.ReadPortal();
        vtkm::Id n = dataPortal.GetNumberOfValues();

        for (vtkm::Id i = 0; i < n; ++i)
        {
            std::cout << i << "\t" << dataPortal.Get(i) << std::endl;
        }


//        auto superparentsPortal = contourTree.Superparents.ReadPortal();
        auto hyperparentsPortal = contourTree.Hyperparents.ReadPortal();
        auto hypernodesPortal = contourTree.Hypernodes.ReadPortal();
        auto hyperarcsPortal = contourTree.Hyperarcs.ReadPortal();
//        auto superarcsPortal = contourTree.Superarcs.ReadPortal();
//        auto nodesPortal = contourTree.Nodes.ReadPortal();


        std::cout << "Euler chi & delta border edges" << std::endl;
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {
            vtkm::Id i_sortID = nodesPortal.Get(sortedNode);
            vtkm::Id i_superparent = superparentsPortal.Get(i_sortID);

            vtkm::Id j_sortID = nodesPortal.Get(sortedNode+1);
            vtkm::Id j_superparent = superparentsPortal.Get(j_sortID );

            vtkm::Id tailend = supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(i_superparent)));

            vtkm::Id delta = 1;

            if (i_superparent == j_superparent)
            {
                tailend = j_sortID;
            }

            vtkm::Id betti0 = 1;
            vtkm::Id betti1 = 1;
            vtkm::Id betti2 = 1;

            if(be_ij[sortedNode] > 0)
            {
                betti2 = 0; // no void, border detected
            }

            betti1 = betti0 + betti2 - chi_xij[sortedNode];

            std::cout << std::setw(10) << i_sortID << "->" << tailend << "\t" << chi_xij[sortedNode] << "\t" << be_ij[sortedNode]
                                  << "\t" << betti0 << "\t" << betti1 << "\t" << betti2;//<< std::endl;

            std::cout << "\t\t" << i_sortID << "\t" << supernodesPortal.Get(i_superparent);
            std::cout << "\t\t" << i_sortID << "\t" << supernodesPortal.Get(i_superparent) << std::endl;


            if(betti1 != previous_betti1)
            {
                regular_nodes_to_insert.push_back(i_sortID);

                nodes_to_relabel_superparent.push_back(i_superparent); // 2026-01-03 addition (original SPs of nodes to be upgraded to supernodes)
                nodes_to_relabel_hyperparent.push_back(hyperparentsPortal.Get(i_superparent)); // 2026-01-03 hyperparents failing when not matched
                nodes_to_relabel_regularID.push_back(i_sortID);

                if(i_sortID > tailend)
                {
                    node_ascend.push_back(-1);
//                    nodes_to_relabel_dataflip.push_back(dataPortal.Get(i_sortID) * -1.0); //data value can be the same, regular ID won't
                        nodes_to_relabel_dataflip.push_back((double)i_sortID * -1.0);
                }
                else
                {
                    node_ascend.push_back(1);
//                    nodes_to_relabel_dataflip.push_back(dataPortal.Get(i_sortID) * 1.0);  //data value can be the same, regular ID won't
                        nodes_to_relabel_dataflip.push_back((double)i_sortID * 1.0);
                }

            }
            else
            {// if betti number didn't change, but dealing with a supernode ...
             // ... still need to relabel it
                if(i_sortID == supernodesPortal.Get(i_superparent))
                {
                    nodes_to_relabel_superparent.push_back(i_superparent);
                    nodes_to_relabel_hyperparent.push_back(hyperparentsPortal.Get(i_superparent)); // 2026-01-03 hyperparents failing when not matched
                    nodes_to_relabel_regularID.push_back(i_sortID);

                    if(i_sortID > tailend)
                    {
                        node_ascend.push_back(-1);
//                        nodes_to_relabel_dataflip.push_back(dataPortal.Get(i_sortID) * -1.0); //data value can be the same, regular ID won't
                        nodes_to_relabel_dataflip.push_back((double)i_sortID * -1.0);
                    }
                    else
                    {
                        node_ascend.push_back(1);
//                        nodes_to_relabel_dataflip.push_back(dataPortal.Get(i_sortID) * 1.0);  //data value can be the same, regular ID won't
                        nodes_to_relabel_dataflip.push_back((double)i_sortID * 1.0);
                    }

                }
            }

            previous_betti1 = betti1;
        }


        std::cout << "Augment the tree with Betti Numbers ..." << std::endl;

//        for(int i = 0; i < ; i++)
//        {
//            vtkm::Id regularId = nodesPortal.Get(sortID);
//            vtkm::Id superparentId = superparentsPortal.Get(regularId);
//            vtkm::Id hyperparentId = hyperparentsPortal.Get(superparentId);
//        }

        for(int i = 0; i < regular_nodes_to_insert.size(); i++)
        {
            std::cout << regular_nodes_to_insert[i] << std::endl;
        }


        // Augment the tree with Betti Numbers
//        HierarchicalAugmenter<MESH_TRIANGULATION_T> betti_augmenter;
//        betti_augmenter.Initialize(block, &hierarchicalTrees[block], &augmentedTrees[block], &meshes[block]);
//        betti_augmenter.BuildAugmentedTree();

        auto augNodesPortal = contourTree.Augmentnodes.ReadPortal();
        auto augArcsPortal = contourTree.Augmentarcs.ReadPortal();
        auto transferPortal = contourTree.WhenTransferred.ReadPortal();

        std::cout << "Augmented Nodes (" << augNodesPortal.GetNumberOfValues() << ")" << std::endl;
        for(int i = 0; i < augNodesPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id augnode = augNodesPortal.Get(i);
            std::cout << i << "\t" << augnode << std::endl;
        }
        std::cout << "Augmented Arcs (" << augArcsPortal.GetNumberOfValues() << ")" << std::endl;
        for(int i = 0; i < augArcsPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id augarc = augArcsPortal.Get(i);
            std::cout << i << "\t" << augarc << std::endl;
        }
        std::cout << "Transferred When (" << transferPortal.GetNumberOfValues() << ")" << std::endl;
        for(int i = 0; i < transferPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << vtkm::worklet::contourtree_augmented::MaskedIndex(transferPortal.Get(i)) << std::endl;
        }

        std::cout << "Nodes" << std::endl;
        for(int sortID = 0; sortID < nodesPortal.GetNumberOfValues(); sortID++)
        {
            vtkm::Id regularId = nodesPortal.Get(sortID);
            std::cout << sortID << "\t" << regularId << std::endl;
        }

        std::cout << "Arcs" << std::endl;
        for(int sortID = 0; sortID < arcsPortal.GetNumberOfValues(); sortID++)
        {
            vtkm::Id regularId = arcsPortal.Get(sortID);
            vtkm::Id maskedArc = vtkm::worklet::contourtree_augmented::MaskedIndex(regularId);
            std::cout << sortID << "\t" << regularId << "\t" << maskedArc << std::endl;
        }


        std::cout << "Supernodes:" << std::endl;
//         auto supernodesPortal = contourTree.Supernodes.ReadPortal();
        for(int i = 0; i < supernodesPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << supernodesPortal.Get(i) << std::endl;
        }
        std::cout << "Superparents:" << std::endl;
        for(int i = 0; i < superparentsPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << superparentsPortal.Get(i) << std::endl;
        }
        std::cout << "Superarcs:" << std::endl;
        for(int i = 0; i < superarcsPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id maskedSuperarc = vtkm::worklet::contourtree_augmented::MaskedIndex(superarcsPortal.Get(i));
            std::cout << i << "\t" << superarcsPortal.Get(i) << "\t" << maskedSuperarc << std::endl;
        }

        auto firstSupernodeIterationPortal = contourTree.FirstSupernodePerIteration.ReadPortal();
        auto firstHypernodeIterationPortal = contourTree.FirstHypernodePerIteration.ReadPortal();

        std::cout << "firstSupernodeIterationPortal:" << std::endl;
        for(int i = 0; i < firstSupernodeIterationPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << firstSupernodeIterationPortal.Get(i) << std::endl;
        }

        std::cout << "Hypernodes:" << std::endl;
        for(int i = 0; i < hypernodesPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id hypernodeID = hypernodesPortal.Get(i);
            vtkm::Id hyperparentID = hyperparentsPortal.Get(hypernodeID);
            std::cout << i << "\t" << hypernodeID << std::endl; //<< "\t" << hyperparentID << std::endl;
        }

        std::cout << "Hyperparents:" << std::endl;
        for(int i = 0; i < hyperparentsPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << ") " << hyperparentsPortal.Get(i) << std::endl;
        }

        std::cout << "Hyperarcs:" << std::endl;
        for(int i = 0; i < hyperarcsPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id maskedHyperarc = vtkm::worklet::contourtree_augmented::MaskedIndex(hyperarcsPortal.Get(i));
            std::cout << i << "\t" << hyperarcsPortal.Get(i) << "\t" << maskedHyperarc << std::endl;
        }

        std::cout << "firstHypernodeIterationPortal:" << std::endl;
        for(int i = 0; i < firstHypernodeIterationPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << firstHypernodeIterationPortal.Get(i) << std::endl;
        }

        std::cout << "node->supernode->superarc(superparent)->hypernode->hyperarc(hyperparent) mappings" << std::endl;
        for(int sortID = 0; sortID < nodesPortal.GetNumberOfValues(); sortID++)
        {
            vtkm::Id regularId = nodesPortal.Get(sortID);
            vtkm::Id superparentId = superparentsPortal.Get(regularId);
            vtkm::Id hyperparentId = hyperparentsPortal.Get(superparentId);

//            std::cout << sortID << ")" << regularID << "->" << superparentID << "(" << hyperparentID << ")" << std::endl;
//            std::cout << sortID << "\t" << regularId << "\t" << superparentId << "\t" << hyperparentId << std::endl;

//            std::cout << "Probing HyperPath\n";
//            std::cout << "Node:        " << sortID << std::endl;
//            std::cout << "Regular ID: ";
//            vtkm::worklet::contourtree_augmented::PrintIndexType(regularId, std::cout);
//            resultStream << "  Value: " << vtkm::cont::ArrayGetValue(regularId, this->DataValues);
//            resultStream << " Global ID: ";
//            vtkm::worklet::contourtree_augmented::PrintIndexType(
//              vtkm::cont::ArrayGetValue(regularId, this->RegularNodeGlobalIds), std::cout);
//            resultStream << " Regular ID: ";
//            vtkm::worklet::contourtree_augmented::PrintIndexType(regularId, std::cout);
//            resultStream << " SNode ID: ";
//            vtkm::worklet::contourtree_augmented::PrintIndexType(
//              vtkm::cont::ArrayGetValue(regularId, this->Regular2Supernode), std::cout);
//            std::cout << "Superparents: ";
//            vtkm::worklet::contourtree_augmented::PrintIndexType(
//              vtkm::cont::ArrayGetValue(regularId, contourTree.Superparents), std::cout << "\t");

//            vtkm::Id hypertarget = vtkm::cont::ArrayGetValue(hyperparentId, contourTree.Hyperarcs);
//            std::cout << "Hypertarget: " << vtkm::cont::ArrayGetValue(hypertarget, contourTree.Hyperarcs) << std::endl;
//            std::cout << "Hypertargets: ";
            vtkm::Id hypertarget = vtkm::cont::ArrayGetValue(hyperparentId, contourTree.Hyperarcs);
            vtkm::Id maskedHypertarget = vtkm::worklet::contourtree_augmented::MaskedIndex(hypertarget);
//            vtkm::worklet::contourtree_augmented::PrintIndexType(
//              vtkm::cont::ArrayGetValue(hypertarget, contourTree.Superparents));
//            std::cout << hypertarget << " - " /*<< vtkm::cont::ArrayGetValue(hypertarget, contourTree.Superparents)*/ << std::endl;

            vtkm::Id supertarget = vtkm::cont::ArrayGetValue(superparentId, contourTree.Superarcs);
            vtkm::Id maskedSupertarget = vtkm::worklet::contourtree_augmented::MaskedIndex(supertarget);
//            std::cout << supertarget << " - " /*<< vtkm::cont::ArrayGetValue(supertarget, contourTree.Superparents)*/ << std::endl;

             std::cout << sortID << "\t" << regularId << "\t" << superparentId << "(" << maskedSupertarget << ")"
                       << "\t" << hyperparentId << "(" << maskedHypertarget << ")" << std::endl;
        }



        std::cout << "The BID array for rehooking up the super{hyper}structure:" << std::endl;
//        int num_betti_change_nodes = 10;
//        int bid_size = supernodesPortal.GetNumberOfValues() + num_betti_change_nodes;

        using vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT;

        // following the template to resize arrays (from HierarchicalAugmenter.h ResizeArrays(vtkm::Id roundNumber):959)
//        vtkm::worklet::contourtree_augmented::ResizeVector(
//          &contourTree.Supernodes,
//          bid_size,
//          vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT);

        // first    : sort on HP
        // secondary: on data value (with its comparitor) with the flip ascending/descending
        // tertiary : sort on regular ID (out of paranoia)
//        std::cout << "id)\thpID\tval_filp\tregularID" << std::endl;
//        for(int i = 0; i < supernodesPortal.GetNumberOfValues(); i++)
//        {
//            vtkm::Id regularId = supernodesPortal.Get(i);
//            vtkm::Id superparentId = superparentsPortal.Get(regularId);
//            vtkm::Id hyperparentId = hyperparentsPortal.Get(superparentId);

////            vtkm::Id i_sortID = supernodesPortal.Get(i+1);
////            vtkm::Id i_superparent = superparentsPortal.Get(i_sortID);

//            vtkm::Id regularId_j = supernodesPortal.Get(i+1);
////            vtkm::Id j_superparent = superparentsPortal.Get(j_sortID );

//            std::cout << i << ")  " << supernodesPortal.Get(i)
//                      << "\t" << hyperparentId
//                      << "\t" << regularId
//                      << "\t(" << supernodesPortal.Get(superparentId) << ")"
//                      << std::endl;
//        }

//        for(int i = 0; i < regular_nodes_to_insert.size(); i++)
//        {
//            vtkm::Id regularId = regular_nodes_to_insert[i];
//            vtkm::Id superparentId = superparentsPortal.Get(regularId);
//            vtkm::Id hyperparentId = hyperparentsPortal.Get(superparentId);

//            std::cout << i << ")  " << regular_nodes_to_insert[i]
//                      << "\t" << hyperparentId
//                      << "\t" << regularId
//                      << "\t(" << supernodesPortal.Get(superparentId) << ")"
//                      << std::endl;
//        }
        std::cout << "id)\thpID\tval_filp\tregularID\tSP" << std::endl;
        for(int i = 0; i < nodes_to_relabel_regularID.size(); i++)
        {
//            vtkm::Id regularId = nodes_to_relabel_regularID[i];
            vtkm::Id superparentId = superparentsPortal.Get(nodes_to_relabel_regularID[i]);
            vtkm::Id hyperparentId = hyperparentsPortal.Get(superparentId);

            std::cout << i << ")  " << nodes_to_relabel_hyperparent[i]
                           << "\t"  << nodes_to_relabel_dataflip[i]
                           << "\t"  << nodes_to_relabel_regularID[i]
                           << "\t(" << supernodesPortal.Get(superparentId) << ")"
                           << "\t"  << nodes_to_relabel_superparent[i]
                           << std::endl;
        }

        vtkm::cont::ArrayHandle<vtkm::Id> ah_super;
        ah_super.Allocate(nodes_to_relabel_superparent.size());
        auto portal_sp = ah_super.WritePortal();
        for (vtkm::Id i = 0; i < vtkm::Id(nodes_to_relabel_superparent.size()); ++i)
            portal_sp.Set(i, nodes_to_relabel_superparent[i]);


        vtkm::cont::ArrayHandle<vtkm::Id> ah_hyper;
        ah_hyper.Allocate(nodes_to_relabel_hyperparent.size());
        auto portal = ah_hyper.WritePortal();
        for (vtkm::Id i = 0; i < vtkm::Id(nodes_to_relabel_hyperparent.size()); ++i)
            portal.Set(i, nodes_to_relabel_hyperparent[i]);

        vtkm::cont::ArrayHandle<double> ah_data;
        ah_data.Allocate(nodes_to_relabel_dataflip.size());
        auto portal_data = ah_data.WritePortal();
        for (vtkm::Id i = 0; i < vtkm::Id(nodes_to_relabel_dataflip.size()); ++i)
            portal_data.Set(i, nodes_to_relabel_dataflip[i]);

        vtkm::cont::ArrayHandle<vtkm::Id> ah_regular;
        ah_regular.Allocate(nodes_to_relabel_regularID.size());
        auto portal_regular = ah_regular.WritePortal();
        for (vtkm::Id i = 0; i < vtkm::Id(nodes_to_relabel_regularID.size()); ++i)
            portal_regular.Set(i, nodes_to_relabel_regularID[i]);

        // first zip 2 arrays
//        auto zipped12 = vtkm::cont::make_ArrayHandleZip(ah_hyper, ah_data);
//        auto zipped34 = vtkm::cont::make_ArrayHandleZip(ah_regular, ah_super); // 2026-01-03

        auto zipped12 = vtkm::cont::make_ArrayHandleZip(ah_super, ah_data);
        auto zipped34 = vtkm::cont::make_ArrayHandleZip(ah_regular, ah_hyper); // 2026-01-03

        // then zip the previous zip with the final array
//        auto zipped123 = vtkm::cont::make_ArrayHandleZip(zipped12, ah_regular); // 2026-01-03
//        auto zipped123 = vtkm::cont::make_ArrayHandleZip(zipped12, ah_regular); // 2026-01-03
        auto zipped123 = vtkm::cont::make_ArrayHandleZip(zipped12, zipped34); // 2026-01-03

        // Sort lexicographically
        vtkm::cont::Algorithm::Sort(zipped123);

        std::cout << "!!!!!!!!!!!!!!!! PREVIOUS SUPERNODE LIST MAX ID (LEN) !!!!!!!!!!!!!!!!" << std::endl;
        int num_original_supernodes = supernodesPortal.GetNumberOfValues();
        std::cout << num_original_supernodes << std::endl;

        std::cout << "!!!!!!!!!!!!!!!! SORTED !!!!!!!!!!!!!!!!" << std::endl;

        std::cout << "i"
                  << "\t" << "HP"
                  << "\t" << "valflip"
                  << "\t" << "regID"
                  << "\t" << "SP"
//                  << "\t" << "isNew"
                  << "\t" << "+1==HP"
                  << "\t" << "HT"
//                  << "\t" << "pfixsum"
                  << "\t" << "relabel"
                  << "\t" << "new ST"
                  << std::endl;

        auto zipPortal = zipped123.ReadPortal();
        int new_nodes_pfix_sum = 0;

        std::vector<vtkm::Id> newSuperIDsRelabelled; // contain either old or new super ID names in single array
        std::vector<vtkm::Id> newSupernodes;

        // array holding the REGULAR IDs of the to-become new supernodes thanks to betti augmentation
        // using NEWSUPERID to index these regular ids
        newSupernodes.resize(nodes_to_relabel_regularID.size()); //20); // hack-resolved

        std::cout << "nodes_to_relabel_regularID size = " << nodes_to_relabel_regularID.size() << std::endl;

        for (vtkm::Id i = 0; i < zipPortal.GetNumberOfValues(); ++i)
        {
            auto triple = zipPortal.Get(i);
            // Because of nested pairs:
            vtkm::Id hyperparent = triple.second.second; //triple.first.first;    // (first of outer pair) -> first of inner pair
            double  dataflip     = triple.first.second;   // (second of inner pair)
            vtkm::Id regularID   = triple.second.first;         // (second of outer pair) triple.second;
            vtkm::Id superparent   = triple.first.first; //  triple.second.second         // (second of outer pair) // 2026-01-03 was triple.second before
            vtkm::Id superID     = superparentsPortal.Get(regularID);

            bool isNew = false;
            vtkm::Id new_superID_relabel = superparent; //hyperparent; // 2026-01-03 actually superparent here

            if(regularID != supernodesPortal.Get(superID))
            {
                isNew = true;
                new_nodes_pfix_sum++;
                new_superID_relabel = (num_original_supernodes-1)+new_nodes_pfix_sum;
            }

            newSuperIDsRelabelled.push_back(new_superID_relabel);

            // superID to regularID NEW mapping:
            newSupernodes[new_superID_relabel] = regularID;

        }

        std::cout << "newSupernodes array: " << newSupernodes.size() << std::endl;
        for(int i = 0; i < newSupernodes.size(); i++)
        {
            std::cout << i << "\t" << newSupernodes[i] << std::endl;
        }

        std::this_thread::sleep_for(std::chrono::seconds(3));

        vtkm::Id num_added_supernodes = newSupernodes.size() - contourTree.Supernodes.GetNumberOfValues();

        bool plus1test = false;
        vtkm::Id newSuperTarget = -1;

        std::vector<vtkm::Id> newSuperTargets;

        std::cout << "i"
                  << "\thyperparent"
                  << "\tdataflip"
                  << "\tregularID"
                  << "\tsupernodesPortal.Get(superID)"
                  << "\tnborSuperparent"
                  << "\tplus1test"
                  << "\tsuperparent"
                  << "\tsupertarget"
                  << "\tnewSuperIDsRelabelled[i]"
                  << "\tnewSuperTargets[i]"
                  << std::endl;

        for (vtkm::Id i = 0; i < zipPortal.GetNumberOfValues(); ++i)
        {
            plus1test = false;
            auto triple = zipPortal.Get(i);
            // Because of nested pairs:
            vtkm::Id hyperparent = triple.second.second; //triple.first.first;    // (first of outer pair) -> first of inner pair
            double  dataflip     = triple.first.second;   // (second of inner pair)
            vtkm::Id regularID   = triple.second.first;         // (second of outer pair) triple.second;
            vtkm::Id superparent   = triple.first.first; //  triple.second.second         // (second of outer pair) // 2026-01-03 was triple.second before

            auto nextTriple = zipPortal.Get(i+1);
//                vtkm::Id nborHyperparent = nextTriple.first.first;    // (first of outer pair) -> first of inner pair
//            vtkm::Id nborSuperparent = nextTriple.second.second;    // 2026-01-09 bug found - this is now HP! (first of outer pair) -> first of inner pair 2026-01-03 use SPs
            vtkm::Id nborSuperparent = nextTriple.first.first;    // (first of outer pair) -> first of inner pair 2026-01-03 use SPs


            if(i+1 >= zipPortal.GetNumberOfValues())
            {
                // if the last element, do nothing
            }
            else
            {
//                auto nextTriple = zipPortal.Get(i+1);
////                vtkm::Id nborHyperparent = nextTriple.first.first;    // (first of outer pair) -> first of inner pair
//                vtkm::Id nborSuperparent = nextTriple.second.second;    // (first of outer pair) -> first of inner pair 2026-01-03 use SPs

                if(superparent == nborSuperparent)
//                    if(hyperparent == nborHyperparent)
                {
                    // segmented test:
                    // if the next row in the sort is on the same superarc, ...
                    // ... set the segmented flag
                    plus1test = true;
                }
            }

            if(plus1test)
            { // if not end of segment yet ...
              // ... set the supertarget as the next new super ID
//                newSuperTarget = newSuperIDsRelabelled[i+1];
                newSuperTargets.push_back(newSuperIDsRelabelled[i+1]);
            }
            else
            {// if it's the end of the segment, ...
             // ... set the supertarget to be the previous hypertarget (not the next supernode)
                // (reached end of the last superarc)
//                newSuperTarget = vtkm::worklet::contourtree_augmented::MaskedIndex(vtkm::cont::ArrayGetValue(hyperparent, contourTree.Hyperarcs));
//                newSuperTargets.push_back(vtkm::worklet::contourtree_augmented::MaskedIndex(vtkm::cont::ArrayGetValue(hyperparent, contourTree.Hyperarcs)));
//                newSuperTargets.push_back(vtkm::worklet::contourtree_augmented::MaskedIndex(vtkm::cont::ArrayGetValue(hyperparent, contourTree.Hyperarcs)));
                newSuperTargets.push_back(vtkm::worklet::contourtree_augmented::MaskedIndex(vtkm::cont::ArrayGetValue(superparent, contourTree.Superarcs)));
            }

            vtkm::Id superID     = superparentsPortal.Get(regularID);

//            bool isNew = false;
//            vtkm::Id new_superID_relabel = hyperparent;

//            if(regularID != supernodesPortal.Get(superID))
//            {
//                isNew = true;
//                new_nodes_pfix_sum++;

//                new_superID_relabel = (num_original_supernodes-1)+new_nodes_pfix_sum;
//            }

            std::cout << i
                      << "\t" << hyperparent
                      << "\t" << dataflip
                      << "\t" << regularID
                      << "\t" << supernodesPortal.Get(superID)
                      << "\t" << nborSuperparent
//                      << "\t" << isNew
                      << "\t" << plus1test
//                      << "\t" << vtkm::worklet::contourtree_augmented::MaskedIndex(vtkm::cont::ArrayGetValue(hyperparent, contourTree.Hyperarcs))
                      << "\t" << superparent
                      << "\t" << vtkm::worklet::contourtree_augmented::MaskedIndex(vtkm::cont::ArrayGetValue(superparent, contourTree.Superarcs))
//                      << "\t" << new_nodes_pfix_sum
                      << "\t" << newSuperIDsRelabelled[i]
                         << "\t" << newSuperTargets[i]
                      << std::endl;
        }

        std::this_thread::sleep_for(std::chrono::seconds(3));


        std::cout << "sortID\tregularID\tSP" << std::endl;
        for(int sortID = 0; sortID < nodesPortal.GetNumberOfValues(); sortID++)
        {
            vtkm::Id regularId = nodesPortal.Get(sortID);
            std::cout << sortID << "\t" << regularId << "\t" << superparentsPortal.Get(regularId) << std::endl;
        }

        auto arcsWritePortal = contourTree.Arcs.WritePortal();
        auto nodesWritePortal = contourTree.Nodes.WritePortal();
//        auto supernodesWritePortal = contourTree.Supernodes.WritePortal();


        vtkm::Id originalSize = superarcsPortal.GetNumberOfValues();

        std::cout << "ORIGINAL Superarcs:" << std::endl;
        for(int i = 0; i < superarcsPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id maskedSuperarc = vtkm::worklet::contourtree_augmented::MaskedIndex(superarcsPortal.Get(i));
            std::cout << i << "\t" << superarcsPortal.Get(i) << "\t" << maskedSuperarc << std::endl;
        }
        std::cout << "Increasing array size to: " << originalSize + num_added_supernodes << std::endl; // hack-resolved
        contourTree.Superarcs.Allocate(contourTree.Superarcs.GetNumberOfValues() + num_added_supernodes, vtkm::CopyFlag::On); // hack-resolved



        auto superarcsWritePortal = contourTree.Superarcs.WritePortal();


        // RELABEL SUPERARCS
        std::cout << "RELABEL SUPERARCS: " << std::endl;
        for(int i = 0; i < newSuperIDsRelabelled.size(); i++)
        {
            // check superarc direction, if ascending, add the IS_ASCENDING flag
            if (newSupernodes[newSuperIDsRelabelled[i]] < newSupernodes[newSuperTargets[i]])
            {
                // RELABEL--
                superarcsWritePortal.Set(newSuperIDsRelabelled[i], newSuperTargets[i] | vtkm::worklet::contourtree_augmented::IS_ASCENDING);
                std::cout << newSuperIDsRelabelled[i] << "\t" << newSupernodes[newSuperIDsRelabelled[i]] << " < "
                          << newSupernodes[newSuperTargets[i]] << "\tIS_ASCENDING" << std::endl;
            }
            else
            {   // id descending, don't need a flag
                // RELABEL--
                superarcsWritePortal.Set(newSuperIDsRelabelled[i], newSuperTargets[i]);
                std::cout << newSuperIDsRelabelled[i] << "\t" << newSupernodes[newSuperIDsRelabelled[i]] << " > "
                          << newSupernodes[newSuperTargets[i]] << std::endl;
            }

            if(newSuperTargets[i] == 0)
            {// if root node, the flag is NO_SUCH_ELEMENT by convention
                // RELABEL--
                superarcsWritePortal.Set(newSuperIDsRelabelled[i], newSuperTargets[i] | vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT);
                std::cout << newSuperIDsRelabelled[i] << "\t" << "\tROOT" << std::endl;
            }
        }

        auto superarcsReinvPortal = contourTree.Superarcs.ReadPortal();
        std::cout << "RELABELLED SUPERARCS: " << superarcsReinvPortal.GetNumberOfValues() << std::endl;
        for(int i = 0; i < superarcsReinvPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id maskedSuperarc = vtkm::worklet::contourtree_augmented::MaskedIndex(superarcsReinvPortal.Get(i));
            std::cout << i << "\t" << superarcsReinvPortal.Get(i) << "\t" << maskedSuperarc << "\t"
                      << vtkm::worklet::contourtree_augmented::IsAscending(superarcsReinvPortal.Get(i)) << std::endl;

//            vtkm::worklet::contourtree_augmented::NoSuchElement for finding NAN (such as root) nodes
//            std::cout << newSupernodes[i] << "\t" << newSupernodes[maskedSuperarc] << std::endl;
        }


        std::vector<vtkm::Id> newSupernodeHyperparents;

        std::cout << "NEW HYPERPARENTS:" << std::endl;
        for(int i = 0; i < newSupernodes.size(); i++)
        {
            std::cout << i << "\t" << newSupernodes[i] << "\t"
                      << superparentsPortal.Get(newSupernodes[i]) << "\t"
                      << hyperparentsPortal.Get(superparentsPortal.Get(newSupernodes[i])) << std::endl;

            newSupernodeHyperparents.push_back(hyperparentsPortal.Get(superparentsPortal.Get(newSupernodes[i])));
        }

        contourTree.Hyperparents.Allocate(contourTree.Hyperparents.GetNumberOfValues() + num_added_supernodes, vtkm::CopyFlag::On);

        auto hyperparentsWritePortal = contourTree.Hyperparents.WritePortal();

        for(int i = 0; i < newSupernodeHyperparents.size(); i++)
        {
            hyperparentsWritePortal.Set(i, newSupernodeHyperparents[i]);
        }

//        auto hyperparentsRelabelledPortal = contourTree.Hyperparents.ReadPortal();

//        std::cout << "RELABELLED HYPERPARENTS" << std::endl;
//        for(int i = 0; i < hyperparentsRelabelledPortal.GetNumberOfValues(); i++)
//        {
//            std::cout << i << "\t " << hyperparentsRelabelledPortal.Get(i) << std::endl;
//        }



        std::cout << "do a regular walk from superarcs to relabel superparents:" << std::endl;
        std::cout << "from\tto\tset\tgetseg" << std::endl;

        std::vector<vtkm::Id> replaceSuperparentsWith;
        std::vector<vtkm::Id> targetSegments;
        std::vector<vtkm::Id> fromRegID;
        std::vector<vtkm::Id> toRegID;

        std::vector<vtkm::Id> oldSuperparents; // used for telling whenTransferred
        oldSuperparents.resize(superarcsReinvPortal.GetNumberOfValues());

        for(int i = 0; i < superarcsReinvPortal.GetNumberOfValues(); i++)
        {
            vtkm::Id maskedSuperarc = vtkm::worklet::contourtree_augmented::MaskedIndex(superarcsReinvPortal.Get(i));
            std::cout << newSupernodes[i] << "\t" << newSupernodes[maskedSuperarc] << "\t" << i <<  "\t"
                      << superparentsPortal.Get(newSupernodes[i]);// << std::endl;

            oldSuperparents[i] = superparentsPortal.Get(newSupernodes[i]);

            if(i != superparentsPortal.Get(newSupernodes[i]))
            {
                targetSegments.push_back(superparentsPortal.Get(newSupernodes[i]));
                std::cout << "\tadded";
                fromRegID.push_back(newSupernodes[i]);
                toRegID.push_back(newSupernodes[maskedSuperarc]);

                replaceSuperparentsWith.push_back(i);
            }

            std::cout << std::endl;
        }

        std::vector<vtkm::Id> segmentA, segmentB;

        std::cout << "RELABEL SUPERPARENTS ... " << std::endl;
        std::cout << "sortID\tregularID\tSP" << std::endl;
        for(int sortID = 0; sortID < nodesPortal.GetNumberOfValues(); sortID++)
        {
            vtkm::Id regularId = nodesPortal.Get(sortID);
            std::cout << sortID << "\t"
                      << regularId << "\t"
                      << superparentsPortal.Get(regularId)
                      << std::endl;

            segmentA.push_back(regularId);
            segmentB.push_back(superparentsPortal.Get(regularId));
        }


//        std::this_thread::sleep_for(std::chrono::seconds(3));

        vtkm::Id targetSegment = 6;

        auto superparentsWritePortal = contourTree.Superparents.WritePortal();

        for(int i = 0; i < targetSegments.size(); i++)
        {
            // Find the range of segment
            auto beginIt = std::lower_bound(
                segmentB.begin(), segmentB.end(), targetSegments[i]);

            auto endIt = std::upper_bound(
                segmentB.begin(), segmentB.end(), targetSegments[i]);

            std::size_t begin = std::distance(segmentB.begin(), beginIt);
            std::size_t end   = std::distance(segmentB.begin(), endIt);

            std::cout << "Segment " << targetSegments[i]
                      << " reg range: [" << fromRegID[i] << ", " << toRegID[i] << ")\n"
                      << " idx range: [" << begin << ", " << end << ")\n";

//            oldSuperparents[replaceSuperparentsWith[i]] =

            // Extract corresponding A values
            for (std::size_t j = begin; j < end; j++)
            {
                std::cout << "A[" << j << "] = " << segmentA[j];// << "\n";
                if(
                    ( (segmentA[j] <= fromRegID[i]) && (segmentA[j] > toRegID[i]))
                    ||
                    ( (segmentA[j] >= fromRegID[i]) && (segmentA[j] < toRegID[i]))
                  )
                {
//                    std::cout << "\t" << superparentsPortal.Get(j) << "\treplace to: " << replaceSuperparentsWith[i];
                    std::cout << "\t" << superparentsPortal.Get(segmentA[j]) << "\treplace to: " << replaceSuperparentsWith[i];

                    // RELABEL SUPERPARENTS:
                    // RELABEL--
                    superparentsWritePortal.Set(segmentA[j], replaceSuperparentsWith[i]);

                }
                std::cout << std::endl;
            }

        }

        auto superparentsRewritePortal = contourTree.Superparents.ReadPortal();

        std::cout << "REPLACED SUPERPARENTS ..." << std::endl;
        for(int i = 0; i < superparentsRewritePortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << superparentsRewritePortal.Get(i) << std::endl;
        }

        std::this_thread::sleep_for(std::chrono::seconds(3));

        std::cout << "NEW SUPERNODES:" << std::endl;
        for(int i = 0; i < newSupernodes.size(); i++)
        {
            std::cout << i << "\t" << newSupernodes[i] << "\t"
                      << superparentsPortal.Get(newSupernodes[i]) << "\t"
                      << hyperparentsPortal.Get(superparentsPortal.Get(newSupernodes[i])) << std::endl;
        }

        std::cout << "Increasing array size to: " << originalSize + num_added_supernodes << std::endl; // hack-resolved
//        vtkm::Id num_original_supernodes already defined ...
//        vtkm::Id num_added_supernodes = newSupernodes.size() - contourTree.Supernodes.GetNumberOfValues();
        contourTree.Supernodes.Allocate(contourTree.Supernodes.GetNumberOfValues() + num_added_supernodes, vtkm::CopyFlag::On);


        std::cout << "RELABEL SUPERNODES:" << std::endl;
        auto supernodesWritePortal = contourTree.Supernodes.WritePortal();
        for(int i = num_original_supernodes;
                i < num_original_supernodes + num_added_supernodes;
                i++)
        {
            std::cout << i << "\t" << newSupernodes[i] << std::endl;
            // RELABEL--
            supernodesWritePortal.Set(i,newSupernodes[i]);
        }

//        std::cout << "RELABELLED SUPERNODES:" << std::endl;
//        auto supernodesRelabelledPortal = contourTree.Supernodes.ReadPortal();
//        for(int i = 0; i < supernodesRelabelledPortal.GetNumberOfValues(); i++)
//        {
//            std::cout << i << "\t" << supernodesRelabelledPortal.Get(i) << std::endl;
//        }

        std::cout << "WhenTransferred" << std::endl;
//                auto whenTransferredWritePortal = contourTree.WhenTransferred.WritePortal();
        auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();

        std::vector<vtkm::Id> whenTransferredVec;

        for(int i = 0; i < oldSuperparents.size(); i++)
        {

            std::cout << i << "\t" << oldSuperparents[i] << "\t"
                      << whenTransferredPortal.Get(oldSuperparents[i]) << "\t"
                      << vtkm::worklet::contourtree_augmented::MaskedIndex(whenTransferredPortal.Get(oldSuperparents[i]))
                      << std::endl;
            whenTransferredVec.push_back(whenTransferredPortal.Get(oldSuperparents[i]));
        }

        contourTree.WhenTransferred.Allocate(contourTree.WhenTransferred.GetNumberOfValues() + num_added_supernodes, vtkm::CopyFlag::On);

        auto whenTransferredWritePortal = contourTree.WhenTransferred.WritePortal();
        std::cout << "RELABEL WhenTransferred (" << contourTree.WhenTransferred.GetNumberOfValues() << ")\n";
        for(int i = 0; i < contourTree.WhenTransferred.GetNumberOfValues(); i++)
        {
            // RELABEL--
            whenTransferredWritePortal.Set(i, whenTransferredVec[i]);
        }

//        std::cout << "RELABELLED WhenTransferred (" << contourTree.WhenTransferred.GetNumberOfValues() << ")\n";
//        auto whenTransferredRelabelledPortal = contourTree.WhenTransferred.ReadPortal();
//        for(int i = 0; i < whenTransferredRelabelledPortal.GetNumberOfValues(); i++)
//        {
//            std::cout << i << "\t" << whenTransferredRelabelledPortal.Get(i) << std::endl;
//        }


//        vtkm::cont::ArrayHandle<vtkm::Id> Ahandle =
//            vtkm::cont::make_ArrayHandle(segmentA, vtkm::CopyFlag::Off);

//        vtkm::cont::ArrayHandle<vtkm::Id> Bhandle =
//            vtkm::cont::make_ArrayHandle(segmentB, vtkm::CopyFlag::Off);

//        // Segment ID(s) to extract
//        vtkm::cont::ArrayHandle<vtkm::Int32> segments2get =
//            vtkm::cont::make_ArrayHandle<vtkm::Int32>({ 6 });

//        vtkm::cont::ArrayHandle<vtkm::Id> outputLow =
//            vtkm::cont::Algorithm::LowerBounds(Bhandle, segments2get);

//        vtkm::cont::ArrayHandle<vtkm::Id> outputHigh =
//            vtkm::cont::Algorithm::UpperBounds(Bhandle, segments2get);

//        // Portals
//        auto aPortal  = Ahandle.ReadPortal();
//        auto olPortal = outputLow.ReadPortal();
//        auto ohPortal = outputHigh.ReadPortal();

//        // ✅ Correct range extraction
//        vtkm::Id begin = olPortal.Get(0);
//        vtkm::Id end   = ohPortal.Get(0);

//        for (vtkm::Id i = begin; i < end; i++)
//        {
//            std::cout << aPortal.Get(i) << std::endl;
//        }

    }


    // 2024-11-25 COMPUTE THE STRUCT COEFFICIENTS VERSION OF THE WEIGHTS WITH COEFFICIENTS
    void static ComputeVolumeWeightsSerialStructCoefficients(const vtkm::cont::DataSet& input, // the coefficient-based version additionally requires tetrahedral connections and vertex coordinates
                                                const ContourTree& contourTree,
                                                const vtkm::Id nIterations,
                                                vtkm::cont::ArrayHandle<Coefficients>& superarcIntrinsicWeightCoeff, // (output)
                                                vtkm::cont::ArrayHandle<Coefficients>& superarcDependentWeightCoeff, // (output)
                                                vtkm::cont::ArrayHandle<Coefficients>& supernodeTransferWeightCoeff, // (output)
                                                vtkm::cont::ArrayHandle<Coefficients>& hyperarcDependentWeightCoeff, // (output)
                                                // Added 2025-01-30
                                                // We use simple weights for the branch decomposition
                                                FloatArrayType& superarcIntrinsicWeight, // (output)
                                                FloatArrayType& superarcDependentWeight, // (output)
                                                FloatArrayType& supernodeTransferWeight, // (output)
                                                FloatArrayType& hyperarcDependentWeight) // (output)
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
//         auto superarcsPortal = contourTree.Superarcs.ReadPortal();
        auto nodesPortal = contourTree.Nodes.ReadPortal();
        // auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal();

        // Files for 3D experiments of Contour Tree Branch volume-based weight computations
        // PACTBD-EDIT-FIXED
        // const std::string filename3D1 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/101-from-2M-sampled-excel-sorted.1-COORDINATES.txt";
        // PACTBD-EDIT-FIXED
        // const std::string filename3D2 = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/2M-parcels-20250225-sorted.1-valued-TETS.txt";

        // get the total number of values to be swept ...
        // ... this will be one more than the total number of datapoints ...
        // ... because we include the region beyond the last isovalue N, as a range [N, +inf)
        // (also, the number of values corresponds to the total number of different data values in a dataset ...
        //  ... because of the simulation of simplicity, which ensures ever data point has a unique value)
        int num_sweep_values = contourTree.Arcs.GetNumberOfValues() + 1;

        // initialise the intrinsic weight array:
        std::cout << "All REGULAR Nodes by sort ID:" << std::endl;
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {// for each sortedNode
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);
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

        // Keep track of all tetrahedra vertices in a sorted list:
        //  ... for example a tet X, Y, Z, W might have vertices with sort IDs:
        //      X=6, Y=5, Z=3, W=7
        //  ... we sort them in increasing order as 3, 5, 6, 7 ...
        //  ... and refer to vertices as A, B, C, D
        //  ... where we then assign A = 3, B = 5, C = 6, D = 7)
        //  (we do not keep track of the original ordering X, Y, Z, W) ...

        // PACTBD-EDIT-FIXED
        // const std::string filename_vtk = "/home/sc17dd/modules/HCTC2024/VTK-m-topology-refactor/VTK-m_TopologyGraph-Refactor/examples/contour-visualiser/delaunay-parcels/200k-from-2M-sampled-excel-sorted.1-withvalues-manual.vtk";

        cont::ArrayHandle<vtkm::Vec3f> coordinatesVTK;
        coordinatesVTK = input.GetPointField("coordinates").GetData().AsArrayHandle<cont::ArrayHandle<vtkm::Vec3f>>();

        // Explicitly interpret as a tetrahedral cell set
        using TetCellSet = vtkm::cont::CellSetSingleType<>;
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

        //  Sort connected tetrahedral vertices in increasing order
        // (then we can assume the first vertex of the tet holds the lowest value)
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

        long double max_a, max_b, max_c, max_d;
        max_a = max_b = max_c = max_d = 0.0l; // -LDBL_MAX;
        long double min_a, min_b, min_c, min_d;
        min_a = min_b = min_c = min_d = LDBL_MAX;

        // Preallocation
        // not using the 2D array, writing directly to vx_delta arrays:
//        std::vector<std::vector<long double>> tet_down_deltas_pfix(num_sweep_values, std::vector<long double>(4, 0.0l));

        // before your loop ensure these are initialized correctly:
        vx_delta_h1_sum.assign(num_sweep_values, 0.0l);
        vx_delta_h2_sum.assign(num_sweep_values, 0.0l);
        vx_delta_h3_sum.assign(num_sweep_values, 0.0l);
        vx_delta_h4_sum.assign(num_sweep_values, 0.0l);


        #pragma omp parallel for
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
            PositionVector a_vol(verticesA, verticesC), // AC already defined, potentially extranneous
                           b_vol(verticesA, verticesD), // AD already defined, potentially extranneous
                           c_vol(verticesA, verticesB); // AB already defined, potentially extranneous
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
            PositionVector a_h1h2_vol = a_vol, b_h1h2_vol = b_vol, c_h1h2_vol = c_vol; // makes copies of AC, AD, AB for interpolation
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

//            for (vtkm::Id i = 0; i < tetlistSorted.size(); ++i)
//            {
//                 // first slab is where the sliced tet produces a triangle face
//                  vx_delta_sum_slab1[tetlistSorted[i][0]] += vx_delta_slab1;
//                  fc_delta_sum_slab1[tetlistSorted[i][0]] += fc_delta_slab1;
//                  ed_delta_sum_slab1[tetlistSorted[i][0]] += ed_delta_slab1;
//                  // second slab is where the sliced tet produces a quad (two triangle faces)
//                  // repeat ...
//                  // third (last) slab is where the sliced tet produces a triangle face
//                  // repeat ...
//            }

//            long double tmp;

//            // Vertex 0
//            tmp = a_h1h2;
//            TRACK_MINMAX(tmp, min_a, max_a);
//            vx_delta_h1_sum[tetlistSorted[i][0]] += tmp;

//            tmp = b_h1h2;
//            TRACK_MINMAX(tmp, min_b, max_b);
//            vx_delta_h2_sum[tetlistSorted[i][0]] += tmp;

//            tmp = c_h1h2;
//            TRACK_MINMAX(tmp, min_c, max_c);
//            vx_delta_h3_sum[tetlistSorted[i][0]] += tmp;

//            tmp = full_tet_vol + d_h1h2_down;
//            TRACK_MINMAX(tmp, min_d, max_d);
//            vx_delta_h4_sum[tetlistSorted[i][0]] += tmp;

//            // Vertex 1
//            tmp = -a_h1h2 + a_h2h3;
//            TRACK_MINMAX(tmp, min_a, max_a);
//            vx_delta_h1_sum[tetlistSorted[i][1]] += tmp;

//            tmp = -b_h1h2 + b_h2h3;
//            TRACK_MINMAX(tmp, min_b, max_b);
//            vx_delta_h2_sum[tetlistSorted[i][1]] += tmp;

//            tmp = -c_h1h2 + c_h2h3;
//            TRACK_MINMAX(tmp, min_c, max_c);
//            vx_delta_h3_sum[tetlistSorted[i][1]] += tmp;

//            tmp = -d_h1h2_down + d_h2h3_down;
//            TRACK_MINMAX(tmp, min_d, max_d);
//            vx_delta_h4_sum[tetlistSorted[i][1]] += tmp;

//            // Vertex 2
//            tmp = -a_h2h3 + a_h3h4;
//            TRACK_MINMAX(tmp, min_a, max_a);
//            vx_delta_h1_sum[tetlistSorted[i][2]] += tmp;

//            tmp = -b_h2h3 + b_h3h4;
//            TRACK_MINMAX(tmp, min_b, max_b);
//            vx_delta_h2_sum[tetlistSorted[i][2]] += tmp;

//            tmp = -c_h2h3 + c_h3h4;
//            TRACK_MINMAX(tmp, min_c, max_c);
//            vx_delta_h3_sum[tetlistSorted[i][2]] += tmp;

//            tmp = -d_h2h3_down + d_h3h4;
//            TRACK_MINMAX(tmp, min_d, max_d);
//            vx_delta_h4_sum[tetlistSorted[i][2]] += tmp;

//            // Vertex 3
//            tmp = -a_h3h4;
//            TRACK_MINMAX(tmp, min_a, max_a);
//            vx_delta_h1_sum[tetlistSorted[i][3]] += tmp;

//            tmp = -b_h3h4;
//            TRACK_MINMAX(tmp, min_b, max_b);
//            vx_delta_h2_sum[tetlistSorted[i][3]] += tmp;

//            tmp = -c_h3h4;
//            TRACK_MINMAX(tmp, min_c, max_c);
//            vx_delta_h3_sum[tetlistSorted[i][3]] += tmp;

//            tmp = -d_h3h4;
//            TRACK_MINMAX(tmp, min_d, max_d);
//            vx_delta_h4_sum[tetlistSorted[i][3]] += tmp;

        }

//        printf("Dynamic range A: %Le - %Le\n", max_a, min_a);
//        printf("Dynamic range B: %Le - %Le\n", max_b, min_b);
//        printf("Dynamic range C: %Le - %Le\n", max_c, min_c);
//        printf("Dynamic range D: %Le - %Le\n", max_d, min_d);

//        printf("Dynamic range A: %Le\n", fabsl(max_a) / fabsl(min_a));
//        printf("Dynamic range B: %Le\n", fabsl(max_b) / fabsl(min_b));
//        printf("Dynamic range C: %Le\n", fabsl(max_c) / fabsl(min_c));
//        printf("Dynamic range D: %Le\n", fabsl(max_d) / fabsl(min_d));




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

        auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal(); // only using after the betti augmented node update


        superarcDependentWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto superarcDependentWeightCoeffPortal = superarcDependentWeightCoeff.WritePortal();

        supernodeTransferWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues());
        auto supernodeTransferWeightCoeffPortal = supernodeTransferWeightCoeff.WritePortal();

        // TODO: Assumption - for now, treat hyperarcs same as superarcs
        // NOTE: 2025-03-03 ran into the above assumption that breaks the code
        hyperarcDependentWeightCoeff.Allocate(contourTree.Hyperarcs.GetNumberOfValues());
//        hyperarcDependentWeightCoeff.Allocate(contourTree.Superarcs.GetNumberOfValues()); // EXPLICIT-SNs
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
        // NOTE: 2025-12-07 ran into the above assumption when simplifying on relabelled superstructure ...
        //       ... as now it treats the hyperstructure as the superstructure
//      removed 2025-12-07:  for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Supernodes.GetNumberOfValues(); sortedNode++)
        for (vtkm::Id hyperNode = 0; hyperNode < contourTree.Hypernodes.GetNumberOfValues(); hyperNode++)
        {// per node in sorted order
            HAlocalDependent.h1 = 0.0l;
            HAlocalDependent.h2 = 0.0l;
            HAlocalDependent.h3 = 0.0l;
            HAlocalDependent.h4 = 0.0l;

//            hyperarcDependentWeightCoeffPortal.Set(sortedNode, HAlocalDependent);
            hyperarcDependentWeightCoeffPortal.Set(hyperNode, HAlocalDependent);
        }// per node in sorted order


        // TODO: Assumption - for now, treat hyperarcs same as superarcs
        // NOTE: 2025-03-03 ran into the above assumption that breaks the code
        // NOTE: 2025-12-07 ran into the above assumption when simplifying on relabelled superstructure ...
        //       ... as now it treats the hyperstructure as the superstructure

        std::cout << "[SWEEP] - Compute Intrinsic, (from deltas, to sums per arc)" << std::endl;

        std::vector<std::vector<vtkm::Id>> iterationSupernodes; // 2025-12-10 hack-resolved
//        std::vector<vtkm::Id> translateSupernodes;

        iterationSupernodes.resize(nIterations+1);
//        iterationSupernodes.resize(contourTree.Rootnode+1);
        vtkm::Id loopIteration = 0;
        vtkm::Id loopSupernode = -1;


        std::cout << "iterationSNs:" << std::endl;
        for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
        {// per node in sorted order, add its weight to its superparent
            vtkm::Id sortID = nodesPortal.Get(sortedNode);
            vtkm::Id superparent = superparentsPortal.Get(sortID);

            SAlocalIntrinsic.h1 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h1 + vx_delta_h1_sum[sortID];
            SAlocalIntrinsic.h2 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h2 + vx_delta_h2_sum[sortID];
            SAlocalIntrinsic.h3 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h3 + vx_delta_h3_sum[sortID];
            SAlocalIntrinsic.h4 = superarcIntrinsicWeightCoeffPortal.Get(superparent).h4 + vx_delta_h4_sum[sortID];

            superarcIntrinsicWeightCoeffPortal.Set(superparent, SAlocalIntrinsic);

            // NEW: add explicit order for processing supernodes ...
            // ... after the Betti augmentation, supernodes are not processed in a continuous order
            // we use the first index as the iteration, following by the array of supernodes to be processed:

            std::cout << sortID << "\t" << superparent << "\t" << loopIteration << "\t" << MaskedIndex(whenTransferredPortal.Get(superparent)) << std::endl;

            if(loopIteration != MaskedIndex(whenTransferredPortal.Get(superparent)))
            {
                loopIteration = MaskedIndex(whenTransferredPortal.Get(superparent));
            }

            if(loopSupernode != superparent) // new segment
            {
                loopSupernode = superparent;
                iterationSupernodes[loopIteration].push_back(loopSupernode);

//                contourTree.translateSupernodes.push_back(loopSupernode);
            }


        }//per node in sorted order, added its weight to its superparent

        std::cout << "iterationSupernodes.size() = " << iterationSupernodes.size() << std::endl;
        std::cout << "iterationSupernodes[0].size() = " << iterationSupernodes[0].size() << std::endl;

        for(int i = 0; i < iterationSupernodes.size(); i++)
        {
            std::cout << i << ")\t";
            for(int j = 0; j < iterationSupernodes[i].size(); j++)
            {
                std::cout << iterationSupernodes[i][j] << "\t";
            }
            std::cout << std::endl;
        }

        std::cout << "contourTree.translateSupernodes.size() = " << contourTree.translateSupernodes.size() << std::endl;

        for(int i = 0; i < contourTree.translateSupernodes.size(); i++)
        {
            std::cout << i << ")\t" << contourTree.translateSupernodes[i] << "\t";
        }

        std::cout << std::endl;


        std::cout << "    " << RED << std::setw(38) << std::left << "TOTAL SWEEP TIME"
                      << ": " << sweep_total_timer.GetElapsedTime() << " seconds" << RESET << std::endl;

        std::cout << "Superparent Local Intrinsic Coeffs:" << std::endl;
        for (vtkm::Id sp = 0; sp < contourTree.Supernodes.GetNumberOfValues(); sp++)
        {
            std::cout << sp << "\t"
                      << superarcIntrinsicWeightCoeffPortal.Get(sp).h1 << "\t"
                      << superarcIntrinsicWeightCoeffPortal.Get(sp).h2 << "\t"
                      << superarcIntrinsicWeightCoeffPortal.Get(sp).h3 << "\t"
                      << superarcIntrinsicWeightCoeffPortal.Get(sp).h4 << std::endl;
        }


std::cout << "// ================================= ITERATIONS =================================== //" << std::endl;

        std::cout << "[HYPERSWEEP] - Compute Transfers, Dependents" << std::endl;

        vtkm::cont::Timer iterations_total_timer;
        iterations_total_timer.Start();

        // now iterate, propagating weights inwards
        // try to run for one more iteration to capture the whole tree
        std::cout << "Iteration: ";


//        // original iteration sequence: 2025-12-10
//        iterationSupernodes.push_back({0, 1, 2, 3});
//        iterationSupernodes.push_back({4, 5});
//        iterationSupernodes.push_back({6, 7});
//        iterationSupernodes.push_back({8});
//        iterationSupernodes.push_back({9});

        // EXPLICIT-SNs

        // relabelled iteration sequence:
//        iterationSupernodes.push_back({0, 1, 2, 3});                // same
//        iterationSupernodes.push_back({4, 5});                      // same
//        iterationSupernodes.push_back({6, 10, 11, 7, 12, 13});      // insert 10 11 (6) and 12 13 (7)
////        iterationSupernodes.push_back({6, 7, 10, 11, 12, 13});      // insert 10 11 (6) and 12 13 (7)
//        iterationSupernodes.push_back({8, 14, 15, 16, 17, 18, 19}); // insert 14 15 16 17 18 19
//        iterationSupernodes.push_back({9});                         // same


        std::cout << "For deriving the iteration vector from the whenTransferredPortal" << std::endl;
        for(int i = 0; i < whenTransferredPortal.GetNumberOfValues(); i++)
        {
            std::cout << i << "\t" << MaskedIndex(whenTransferredPortal.Get(i)) << std::endl;
        }

//        for(int i = 0; i < nIterations; i++)
//        {

//        }

        for (vtkm::Id iteration = 0; iteration < nIterations+1; iteration++)
        {// per iteration
          std::cout << iteration << " ";

          // pull the array bounds into register
          vtkm::Id firstSupernode = firstSupernodePerIterationPortal.Get(iteration);
          vtkm::Id lastSupernode = firstSupernodePerIterationPortal.Get(iteration + 1);
          vtkm::Id firstHypernode = firstHypernodePerIterationPortal.Get(iteration);
          vtkm::Id lastHypernode = firstHypernodePerIterationPortal.Get(iteration + 1);

          // 2025-12-08 Betti number hyperstructure relabel with simplification:
          // The first-last supernode is no longer guaranteed to be in a continuous sequence ...
          // ... therefore, we replace the {first/last}Supernode id's with an explicit array ...
          // ... for each iteration




          std::cout << "first -> last (SN)\t" << firstSupernode << "\t" << lastSupernode << std::endl;
          std::cout << "first -> last (HN)\t" << firstHypernode << "\t" << lastHypernode << std::endl;

          std::cout << "Explicit iteration-supernode pairs:" << std::endl;

          for(int i = 0; i < iterationSupernodes[iteration].size(); i++)
          {
              std::cout << iterationSupernodes[iteration][i] << "\t";
          }
          std::cout << std::endl;

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

          std::cout << "Tailends:" << std::endl;
          for (vtkm::Id supernode = 0; supernode < contourTree.Supernodes.GetNumberOfValues(); supernode++)
          {
              vtkm::Id superNode = supernodesPortal.Get(supernode);
              tailends.insert(std::make_pair(superNode, supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode)))));

              std::cout << superNode << "\t" << supernodesPortal.Get(MaskedIndex(superarcsPortal.Get(supernode))) << std::endl;
          }

          // so, step 1: add xfer + int & store in dependent weight
//          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++) EXPLICIT-SNs
          for(int i = 0; i < iterationSupernodes[iteration].size(); i++)
          {
            vtkm::Id supernode = iterationSupernodes[iteration][i];

            vtkm::Id superNode = supernodesPortal.Get(supernode);
            superarcDependentWeightPortal.Set(supernode,
                                              supernodeTransferWeightPortal.Get(supernode) +
                                                superarcIntrinsicWeightPortal.Get(supernode));
          }

          // so, step 1: add xfer + int & store in dependent weight
//          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++) EXPLICIT-SNs
          for(int i = 0; i < iterationSupernodes[iteration].size(); i++)
          {
            vtkm::Id supernode = iterationSupernodes[iteration][i];

//            std::cout << "(@" << supernode << ")" << std::endl;

            step1Dependent.h1 = supernodeTransferWeightCoeffPortal.Get(supernode).h1 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h1;
            step1Dependent.h2 = supernodeTransferWeightCoeffPortal.Get(supernode).h2 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h2;
            step1Dependent.h3 = supernodeTransferWeightCoeffPortal.Get(supernode).h3 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h3;
            step1Dependent.h4 = supernodeTransferWeightCoeffPortal.Get(supernode).h4 + superarcIntrinsicWeightCoeffPortal.Get(supernode).h4;

            superarcDependentWeightCoeffPortal.Set(supernode, step1Dependent);
          }

          //          std::cout << firstSupernode << " - " << superarcDependentWeightPortal.Get(firstSupernode) << std::endl;
          // step 2: perform prefix sum on the dependent weight range
//          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++) EXPLICIT-SNs
          for(int i = 1; i < iterationSupernodes[iteration].size(); i++)
          {
            vtkm::Id supernode         = iterationSupernodes[iteration][i];
            vtkm::Id supernodePrevious = iterationSupernodes[iteration][i-1]; // 2025-12-07 fix - SAs no longer continuous

//            superarcDependentWeightPortal.Set(supernode,
//                                              superarcDependentWeightPortal.Get(supernode) +
//                                                superarcDependentWeightPortal.Get(supernode - 1));

            superarcDependentWeightPortal.Set(supernode,
                                              superarcDependentWeightPortal.Get(supernode) +
                                                superarcDependentWeightPortal.Get(supernodePrevious));

            //std::cout << superarcDependentWeightPortal.Get(supernode) << " "; // + superarcDependentWeightPortal.Get(supernode - 1) << " ";
          }

          // step 2: perform prefix sum on the dependent weight range
//          for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++) EXPLICIT-SNs
          for(int i = 1; i < iterationSupernodes[iteration].size(); i++)
          {
            vtkm::Id supernode = iterationSupernodes[iteration][i];
            vtkm::Id supernodePrevious = iterationSupernodes[iteration][i-1]; // 2025-12-07 fix - SAs no longer continuous

              step2Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 + superarcDependentWeightCoeffPortal.Get(supernodePrevious).h1;
              step2Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 + superarcDependentWeightCoeffPortal.Get(supernodePrevious).h2;
              step2Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 + superarcDependentWeightCoeffPortal.Get(supernodePrevious).h3;
              step2Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 + superarcDependentWeightCoeffPortal.Get(supernodePrevious).h4;

//              2025-12-07 removed:
//              step2Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 + superarcDependentWeightCoeffPortal.Get(supernode-1).h1;
//              step2Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 + superarcDependentWeightCoeffPortal.Get(supernode-1).h2;
//              step2Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 + superarcDependentWeightCoeffPortal.Get(supernode-1).h3;
//              step2Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 + superarcDependentWeightCoeffPortal.Get(supernode-1).h4;

              superarcDependentWeightCoeffPortal.Set(supernode, step2Dependent);
          }




// !!! CONVERTING FROM COEFFICIENT BASED TO SIMPLE(VALUE-TYPE) FOR DEPENDENT WEIGHTS !!! //
std::cout << "!!! CONVERTING FROM COEFFICIENT BASED TO SIMPLE(VALUE-TYPE) FOR DEPENDENT WEIGHTS !!!" << std::endl;
//for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++) // EXPLICIT-SNs
for(int i = 0; i < iterationSupernodes[iteration].size(); i++)
{
    vtkm::Id supernode = iterationSupernodes[iteration][i];

//    std::cout << "(" << supernode << ")" << std::endl;

    vtkm::Id superNode = supernodesPortal.Get(supernode);

    long double a_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h1 * std::pow(tailends[superNode], 3);
    long double b_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h2 * std::pow(tailends[superNode], 2);
    long double c_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h3 * tailends[superNode];
    long double d_coeff = superarcDependentWeightCoeffPortal.Get(supernode).h4;

    // NEW: 2025-01-24 actually save the value to the dependent array SADWP:
    superarcDependentWeightPortal.Set(supernode, a_coeff + b_coeff + c_coeff + d_coeff);

    std::cout << supernode << "\t" << a_coeff + b_coeff + c_coeff + d_coeff << std::endl;

}


//        std::vector<vtkm::Id> hyperparents = {0,1,2,3,4,5,6,7,8,9,// }; //,
//                                              6,6,7,7,8,8,8,8,8,8}; // 2025-12-10 with betti hack-resolved

//        std::cout << "Manual vs real hyperparents:" << std::endl;
//        for(int i = 0; i < hyperparents.size(); i++)
//        {
//            std::cout << i << "\t" << hyperparents[i] << "\t" << hyperparentsPortal.Get(i) << std::endl;
//        }



          // step 3: subtract out the dependent weight of the prefix to the entire hyperarc.
          // This will be a transfer, but for now, it's easier to show it in serial.
          // NB: Loops backwards so that computation uses the correct value
          // As a bonus, note that we test > firstsupernode, not >=.
          // This is because we've got unsigned integers, & otherwise it will not terminate
          // But the first is always correct anyway (same reason as the short-cut termination on hyperparent), so we're fine
//          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--) EXPLICIT-SNs
          std::cout << "UPDATED SN DEPENDENT:" << std::endl;
          for(int i = iterationSupernodes[iteration].size()-1; i > 0; i--)
          {
            vtkm::Id supernode = iterationSupernodes[iteration][i];

//          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            std::cout << "sn(" << supernode << ") hp(" << hyperparent << ") hpsid(" << hyperparentSuperID << ") ";

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
            {
              std::cout << std::endl;
              continue;
            }


            // 2025-12-16 - Rewrite of 2025-12-07 to use array handles directly:
            vtkm::Id lastSuperarc;
            for (vtkm::Id i = hyperparentsPortal.GetNumberOfValues() - 1; i >= 0; --i)
            {
              if (hyperparentsPortal.Get(i) == hyperparentSuperID - 1)
              {
                lastSuperarc = i;
                break;
              }
            }

            if (lastSuperarc >= 0)
            {
//                vtkm::Id lastHPSuperarc = hyperparents.size() - 1 - std::distance(hyperparents.rbegin(), it);

                superarcDependentWeightPortal.Set(
                  supernode,
                  superarcDependentWeightPortal.Get(supernode) -
                    superarcDependentWeightPortal.Get(lastSuperarc));

                std::cout << supernode << "\t" << superarcDependentWeightPortal.Get(supernode) << std::endl;

            }
            else
            {
              std::cout << "Value not found\n";
            }


//            // 2025-12-07:
//            auto it = std::find(hyperparents.rbegin(), hyperparents.rend(), hyperparentSuperID - 1);
//            if (it != hyperparents.rend())
//            {
//                vtkm::Id lastHPSuperarc = hyperparents.size() - 1 - std::distance(hyperparents.rbegin(), it);
//                std::cout << "lasthp(" << lastHPSuperarc << ")\n";

//                superarcDependentWeightPortal.Set(
//                  supernode,
//                  superarcDependentWeightPortal.Get(supernode) -
//                    superarcDependentWeightPortal.Get(lastHPSuperarc));

//                std::cout << supernode << "\t" << superarcDependentWeightPortal.Get(supernode) << std::endl;

//            }
//            else {
//                std::cout << "Value not found\n";
//            }
            // end of 2025-12-07

//            // old
//            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
//            superarcDependentWeightPortal.Set(
//              supernode,
//              superarcDependentWeightPortal.Get(supernode) -
//                superarcDependentWeightPortal.Get(hyperparentSuperID - 1));

//            std::cout << supernode << "\t" << superarcDependentWeightPortal.Get(supernode) << std::endl;


          } // per supernode

//          for (vtkm::Id supernode = lastSupernode - 1; supernode > firstSupernode; supernode--) EXPLICIT-SNs
        for(int i = iterationSupernodes[iteration].size()-1; i > 0; i--)
        {
          vtkm::Id supernode = iterationSupernodes[iteration][i];
//          { // per supernode
            // retrieve the hyperparent & convert to a supernode ID
            vtkm::Id hyperparent = hyperparentsPortal.Get(supernode);
            vtkm::Id hyperparentSuperID = hypernodesPortal.Get(hyperparent);

            // if the hyperparent is the first in the sequence, dependent weight is already correct
            if (hyperparent == firstHypernode)
              continue;

            // 2025-12-16 - Rewrite of 2025-12-07 to use array handles directly:
            vtkm::Id lastSuperarc;
            for (vtkm::Id i = hyperparentsPortal.GetNumberOfValues() - 1; i >= 0; --i)
            {
              if (hyperparentsPortal.Get(i) == hyperparentSuperID - 1)
              {
                lastSuperarc = i;
                break;
              }
            }

            if (lastSuperarc >= 0)
            {
//                vtkm::Id lastHPSuperarc = hyperparents.size() - 1 - std::distance(hyperparents.rbegin(), it);
                std::cout << "lasthp(" << lastSuperarc << ")\n";

                step3Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 - superarcDependentWeightCoeffPortal.Get(lastSuperarc).h1;
                step3Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 - superarcDependentWeightCoeffPortal.Get(lastSuperarc).h2;
                step3Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 - superarcDependentWeightCoeffPortal.Get(lastSuperarc).h3;
                step3Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 - superarcDependentWeightCoeffPortal.Get(lastSuperarc).h4;

                std::cout << supernode << "\t" << superarcDependentWeightPortal.Get(supernode) << std::endl;
            }
            else
            {
              std::cout << "Value not found\n";
            }


//            // 2025-12-07:
//            auto it = std::find(hyperparents.rbegin(), hyperparents.rend(), hyperparentSuperID - 1);
//            if (it != hyperparents.rend())
//            {
//                vtkm::Id lastHPSuperarc = hyperparents.size() - 1 - std::distance(hyperparents.rbegin(), it);
//                std::cout << "lasthp(" << lastHPSuperarc << ")\n";

//                step3Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 - superarcDependentWeightCoeffPortal.Get(lastHPSuperarc).h1;
//                step3Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 - superarcDependentWeightCoeffPortal.Get(lastHPSuperarc).h2;
//                step3Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 - superarcDependentWeightCoeffPortal.Get(lastHPSuperarc).h3;
//                step3Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 - superarcDependentWeightCoeffPortal.Get(lastHPSuperarc).h4;

//                std::cout << supernode << "\t" << superarcDependentWeightPortal.Get(supernode) << std::endl;

//            }
//            else {
//                std::cout << "Value not found\n";
//            }


//            // old
//            step3Dependent.h1 = superarcDependentWeightCoeffPortal.Get(supernode).h1 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h1;
//            step3Dependent.h2 = superarcDependentWeightCoeffPortal.Get(supernode).h2 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h2;
//            step3Dependent.h3 = superarcDependentWeightCoeffPortal.Get(supernode).h3 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h3;
//            step3Dependent.h4 = superarcDependentWeightCoeffPortal.Get(supernode).h4 - superarcDependentWeightCoeffPortal.Get(hyperparentSuperID - 1).h4;

            // otherwise, subtract out the dependent weight *immediately* before the hyperparent's supernode
            superarcDependentWeightCoeffPortal.Set(supernode, step3Dependent);

          } // per supernode

// 2025-12-07 moved one loop up
//        std::vector<vtkm::Id> hyperparents = {0,1,2,3,4,5,6,7,8,9,
//                                              6,6,7,7,8,8,8,8,8,8};

          // step 4: transfer the dependent weight to the hyperarc's target supernode
          std::cout << "HYPERNODE DEPENDENTS:" << std::endl;
          for (vtkm::Id hypernode = firstHypernode; hypernode < lastHypernode; hypernode++)
          { // per hypernode
            // last superarc for the hyperarc
            vtkm::Id lastSuperarc;
            // special case for the last hyperarc
            if (hypernode == contourTree.Hypernodes.GetNumberOfValues() - 1)
              // take the last superarc in the array
              lastSuperarc = contourTree.Rootnode; //9; // contourTree. Supernodes.GetNumberOfValues() - 1;
            else
            {
                // otherwise, take the next hypernode's ID and subtract 1
                // old lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1; // 2025-12-07 9;

                // 2025-12-16 - Rewrite of 2025-12-07 to use array handles directly:
                for (vtkm::Id i = hyperparentsPortal.GetNumberOfValues() - 1; i >= 0; --i)
                {
                  if (hyperparentsPortal.Get(i) == hypernode)
                  {
                    lastSuperarc = i;
                    break;
                  }
                }

                if (lastSuperarc >= 0)
                {
                  std::cout << "Last index: " << lastSuperarc << "\n";
                }
                else
                {
                  std::cout << "Value not found\n";
                }


//                // 2025-12-07 with std::vector:
//                auto it = std::find(hyperparents.rbegin(), hyperparents.rend(), hypernode);

//                if (it != hyperparents.rend())
//                {
//                    lastSuperarc = hyperparents.size() - 1 - std::distance(hyperparents.rbegin(), it);
//                    std::cout << "Last index: " << lastSuperarc << "\n";
//                }
//                else {
//                    std::cout << "Value not found\n";
//                }
            }

            // now, given the last superarc for the hyperarc, transfer the dependent weight

            hyperarcDependentWeightPortal.Set(hypernode,
                                              superarcDependentWeightPortal.Get(lastSuperarc));

            std::cout << "hn-" << hypernode << "\t" << hyperarcDependentWeightPortal.Get(hypernode) << std::endl;

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
              lastSuperarc = contourTree.Rootnode; // 9; //contourTree.Supernodes.GetNumberOfValues() - 1;
            else
            {
              // otherwise, take the next hypernode's ID and subtract 1
              // old lastSuperarc = hypernodesPortal.Get(hypernode + 1) - 1;

              // 2025-12-16 - Rewrite of 2025-12-07 to use array handles directly:
              for (vtkm::Id i = hyperparentsPortal.GetNumberOfValues() - 1; i >= 0; --i)
              {
                if (hyperparentsPortal.Get(i) == hypernode)
                {
                  lastSuperarc = i;
                  break;
                }
              }

              if (lastSuperarc >= 0)
              {
                std::cout << "Last index: " << lastSuperarc << "\n";
              }
              else
              {
                std::cout << "Value not found\n";
              }
              // 2025-12-16 end

//              // 2025-12-07:
//                auto it = std::find(hyperparents.rbegin(), hyperparents.rend(), hypernode);

//                if (it != hyperparents.rend())
//                {
//                    lastSuperarc = hyperparents.size() - 1 - std::distance(hyperparents.rbegin(), it);
//                    std::cout << "Last index: " << lastSuperarc << "\n";
//                }
//                else {
//                    std::cout << "Value not found\n";
//                }
//              // 2025-12-07 end

            }

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


//          for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++) EXPLICIT-SNs
          for(int i = 0; i < iterationSupernodes[iteration].size(); i++)
          {
            vtkm::Id supernode = iterationSupernodes[iteration][i];

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
std::cout << "RECOMPUTE THE INTRINSIC (FLOATING POINT) FROM THE CORRECT DEPENDENT AND TRANSFER CONNECTIONS" << std::endl;
for(int i = 0; i < superarcIntrinsicWeightPortal.GetNumberOfValues(); i++)
{
    superarcIntrinsicWeightPortal.Set(i, superarcDependentWeightPortal.Get(i) - supernodeTransferWeightPortal.Get(i));

    std::cout << i << ")\t" << superarcDependentWeightPortal.Get(i) << " - "
              << supernodeTransferWeightPortal.Get(i) << " = "
              << superarcIntrinsicWeightPortal.Get(i) << std::endl;

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
      // crash function

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

      // 2025-12-07 prior to betti number implementation, ...
      // ... the total volume used to be the dependent weight of the last supernode
      // Now with betti number augmentation there are superarcs labelled beyond the previous last supernode ...
      // ... so we set the total volume manually:
      // 2025-12-12 add the ROOT node explicitly, ...
      // ... since it is not guaranteed to be the last one after betti augmentation
//      vtkm::Id rootNode = 9; // hack-resolved
      totalVolumeFloat = superarcDependentWeightPortal.Get(contourTree.Rootnode); // 2025-12-07 hack-resolved
//      totalVolumeFloat = superarcDependentWeightPortal.Get(superarcDependentWeightPortal.GetNumberOfValues()-1); old

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

      // 2025-12-12 Adding a translation array for dealing with supernodes not in a continuous sequence
//      std::vector<vtkm::Id> translateSupernodes = {0,1,2,3,4,5,6,10,11,7,12,13,8,14,15,16,17,18,19,9}; // hack-resolved

      std::vector<vtkm::Id> translateSupernodes;

      auto superparentsPortal = contourTree.Superparents.ReadPortal();
      auto nodesPortal = contourTree.Nodes.ReadPortal();
      auto whenTransferredPortal = contourTree.WhenTransferred.ReadPortal(); // only using after the betti augmented node update

      vtkm::Id loopIteration = 0;
      vtkm::Id loopSupernode = -1;

      std::cout << "translateSupernodes:" << std::endl;
      for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
      {// per node in sorted order, add its weight to its superparent
          vtkm::Id sortID = nodesPortal.Get(sortedNode);
          vtkm::Id superparent = superparentsPortal.Get(sortID);

          // NEW: add explicit order for processing supernodes ...
          // ... after the Betti augmentation, supernodes are not processed in a continuous order
          // we use the first index as the iteration, following by the array of supernodes to be processed:

          std::cout << sortID << "\t" << superparent << "\t" << loopIteration << "\t" << MaskedIndex(whenTransferredPortal.Get(superparent)) << std::endl;

          if(loopIteration != MaskedIndex(whenTransferredPortal.Get(superparent)))
          {
              loopIteration = MaskedIndex(whenTransferredPortal.Get(superparent));
          }

          if(loopSupernode != superparent) // new segment
          {
              loopSupernode = superparent;
              translateSupernodes.push_back(loopSupernode);
          }


      }//per node in sorted order, added its weight to its superparent

      std::cout << "translateSupernodes.size() = " << translateSupernodes.size() << std::endl;

      for(int i = 0; i < translateSupernodes.size(); i++)
      {
//          std::cout << i << ")\t" << translateSupernodes[i] << "\t";
          std::cout << translateSupernodes[i] << std::endl; //<< "\t";
      }
      std::cout << std::endl;

//      // original iteration sequence: 2025-12-10
//      translateSupernodes.push_back({0, 1, 2, 3});
//      translateSupernodes.push_back({4, 5});
//      translateSupernodes.push_back({6, 7});
//      translateSupernodes.push_back({8});
//      translateSupernodes.push_back({9});

//      translateSupernodes.push_back(0); //, 1, 2, 3});                // same
//      translateSupernodes.push_back(1); //, 2, 3});                // same
//      translateSupernodes.push_back(2); //, 3});                // same
//      translateSupernodes.push_back(3);                 // same
//      translateSupernodes.push_back(4); // 5                 // same
//      translateSupernodes.push_back(5); // 5                 // same
//      translateSupernodes.push_back(6); // 10, 11, 7, 12, 13                // same
//      translateSupernodes.push_back(10); // 11, 7, 12, 13                // same
//      translateSupernodes.push_back(11); // 7, 12, 13                // same
//      translateSupernodes.push_back(7); // 12, 13                // same
//      translateSupernodes.push_back(12); // 13                // same
//      translateSupernodes.push_back(13);                 // same
//      translateSupernodes.push_back(8); // 14, 15, 16, 17, 18, 19                 // same
//      translateSupernodes.push_back(14); // 15, 16, 17, 18, 19                 // same
//      translateSupernodes.push_back(15); // 15, 16, 17, 18, 19                 // same
//      translateSupernodes.push_back({8, 14, 15, 16, 17, 18, 19}); // insert 14 15 16 17 18 19
//      translateSupernodes.push_back({9});                         // same

      // NB: Last element in array is guaranteed to be root superarc to infinity,
      // WARNING WARNING WARNING: This changes in the distributed version!
      // so we can easily skip it by not indexing to the full size
      // i.e we loop to N superarcs, not to N supernodes since N superarcs = supernodes-1
//      for (vtkm::Id superarc = 0; superarc <= nSuperarcs; superarc++) //  for (vtkm::Id supernode = 0; supernode < nsupernodes-1; supernode++) vtkm::Id supernode = superarc (temporary variable)
//      { // per superarc

      std::cout << "translated vs sequence: " << translateSupernodes.size() << "\t" << nSuperarcs << "(" << contourTree.Rootnode << ")" << std::endl;
      for(int i = 0; i < translateSupernodes.size(); i++)
      {// per superarc (translated by iteration)
          vtkm::Id superarc = translateSupernodes[i];
          vtkm::Id superarcSeqID = i;

          std::cout << superarc << "\t" << superarcSeqID << std::endl;

          //jump

          // after implementing betti number augmentation, the last superarc in the list is no longer guaranteed to be the root ...
          // ... so we check here explicitly and make the loop <= nSuperarcs instead of < nSuperarcs
          if(superarc != contourTree.Rootnode) // 2025-12-12
          { // if not root supernode


    #if DEBUG_PRINT_PACTBD
            std::cout << "Processing superarc: " << superarc << std::endl << "{" << std::endl;
    #endif
            if (IsAscending(superarcsPortal.Get(superarc))) // flag on the ID of superarc
            { // ascending superarc
    #if DEBUG_PRINT_PACTBD
              std::cout << indent << "ASCENDING\n";
    #endif
                // put the lower-end first
              superarcListWritePortal.Set(superarcSeqID, // each superarc starts at the supernode at the same ID and goes to the supernode whose ID is stored in the superarc's array
                                          // pair the origin and the destination of that superarc. We store them in an edge pair and write it to the array
                                          EdgePair(superarc, MaskedIndex(superarcsPortal.Get(superarc))));

              vtkm::Id superNode = supernodesPortal.Get(superarc);
    #if DEBUG_PRINT_PACTBD
              std::cout << indent << superarc << " = " << superNode << " -> " << MaskedIndex(superarcsPortal.Get(superarc)) << std::endl << std::endl;
    #endif
              // because this is an ascending superarc, the dependent weight refers to the weight at the upper end
              // so, we set the up weight based on the dependent weight
    //          upWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
              upWeightFloatPortal.Set(superarcSeqID, superarcDependentWeightPortal.Get(superarc));
              upWeightFloatCorrectPortal.Set(superarcSeqID, superarcDependentWeightCorrectReadPortal.Get(superarc));


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

              downWeightFloatPortal.Set(superarcSeqID,
                                       // below is the invert operator for node count!
                                       (totalVolumeFloat - superarcDependentWeightPortal.Get(superarc)) +
                                         (superarcIntrinsicWeightPortal.Get(superarc) - 1));

              downWeightFloatCorrectPortal.Set(superarcSeqID,
                                       // below is the invert operator for node count!
                                       (totalVolumeFloat - superarcDependentWeightCorrectReadPortal.Get(superarc)) +
                                         (superarcIntrinsicWeightCorrectReadPortal.Get(superarc) - 1));
    #if DEBUG_PRINT_PACTBD
    //          std::cout << indent << "upWeightPortal             = " << upWeightPortal.Get(superarc) << std::endl;
              std::cout << indent << "upWeightFloatPortal        = " << upWeightFloatPortal.Get(superarcSeqID) << std::endl;
              std::cout << indent << "upWeightFloatCorrectPortal = " << upWeightFloatCorrectPortal.Get(superarcSeqID) << std::endl;
              std::cout << std::endl;
    //          std::cout << indent << "downWeightPortal             = " << downWeightPortal.Get(superarc) << std::endl;
              std::cout << indent << "downWeightFloatPortal        = " << downWeightFloatPortal.Get(superarcSeqID) << std::endl;
              std::cout << indent << "downWeightFloatCorrectPortal = " << downWeightFloatCorrectPortal.Get(superarcSeqID) << std::endl;

              std::cout << "\n}\n";
    #endif
            } // ascending superarc
            else
            { // descending superarc
    #if DEBUG_PRINT_PACTBD
              std::cout << indent << "DESCENDING\n";
    #endif
              // lower-end is also first, but in the reverse order compared to IsAscending
              superarcListWritePortal.Set(superarcSeqID,
                                          EdgePair(MaskedIndex(superarcsPortal.Get(superarc)), superarc));

              vtkm::Id superNode = supernodesPortal.Get(superarc);
    #if DEBUG_PRINT_PACTBD
              std::cout << indent << superarc << " = " << superNode << " -> " << MaskedIndex(superarcsPortal.Get(superarc)) << std::endl << std::endl;
    #endif

    //          downWeightPortal.Set(superarc, superarcDependentWeightPortal.Get(superarc));
              downWeightFloatPortal.Set(superarcSeqID, superarcDependentWeightPortal.Get(superarc));
              downWeightFloatCorrectPortal.Set(superarcSeqID, superarcDependentWeightCorrectReadPortal.Get(superarc));

              // at the inner end, dependent weight is the total in the subtree.
              // Then there are vertices along the edge itself (intrinsic weight), ...
              // ... including the supernode at the outer end
              // So, to get the "dependent" weight in the other direction, ...
              // ... we start with totalVolume - dependent, then subtract (intrinsic - 1)
    //          upWeightPortal.Set(superarc,
    //                             (totalVolume - superarcDependentWeightPortal.Get(superarc)) +
    //                               (superarcIntrinsicWeightPortal.Get(superarc) - 1));

              upWeightFloatPortal.Set(superarcSeqID,
                                 (totalVolumeFloat - superarcDependentWeightPortal.Get(superarc)) +
                                   (superarcIntrinsicWeightPortal.Get(superarc) - 1));

              upWeightFloatCorrectPortal.Set(superarcSeqID,
                                       // below is the invert operator for node count!
                                       (totalVolumeFloat - superarcDependentWeightCorrectReadPortal.Get(superarc)) +
                                         (superarcIntrinsicWeightCorrectReadPortal.Get(superarc) - 1));
    #if DEBUG_PRINT_PACTBD
    //          std::cout << indent << "upWeightPortal      = " << upWeightPortal.Get(superarc) << std::endl;
              std::cout << indent << "upWeightFloatPortal = " << upWeightFloatPortal.Get(superarcSeqID) << std::endl;
              std::cout << indent << "upWeightFloatCorrectPortal = " << upWeightFloatCorrectPortal.Get(superarcSeqID) << std::endl;
              std::cout << std::endl;
    //          std::cout << indent << "downWeightPortal      = " << downWeightPortal.Get(superarc) << std::endl;
              std::cout << indent << "downWeightFloatPortal = " << downWeightFloatPortal.Get(superarcSeqID) << std::endl;
              std::cout << indent << "downWeightFloatCorrectPortal = " << downWeightFloatCorrectPortal.Get(superarcSeqID) << std::endl;

              std::cout << "\n}\n";
    #endif

            } // descending superarc


        } // end else if not root supernode

      }   // per superarc

#if DEBUG_PRINT_PACTBD
    std::cout << "Up Weights (float)" << std::endl;
    for (vtkm::Id superarc = 0; superarc <= nSuperarcs; superarc++)
//        for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    {
        if(superarc != contourTree.Rootnode) // 2025-12-12
        {
            vtkm::Id superNode = supernodesPortal.Get(superarc);

            std::cout << superarc << "(" << superNode  << ") = " << upWeightFloatCorrectPortal.Get(superarc) << std::endl;
        }
    }
    std::cout << "Down Weights (float)" << std::endl;
    for (vtkm::Id superarc = 0; superarc <= nSuperarcs; superarc++)
//        for (vtkm::Id superarc = 0; superarc < nSuperarcs; superarc++)
    {
        if(superarc != contourTree.Rootnode) // 2025-12-12
        {
            vtkm::Id superNode = supernodesPortal.Get(superarc);

            std::cout << superarc << "(" << superNode  << ") = " << downWeightFloatCorrectPortal.Get(superarc) << std::endl;
        }
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
      std::cout << "II A Weights Computed" << std::endl;
      std::cout << "(Rootnode) = " << contourTree.Rootnode << std::endl;
      PrintHeader(upWeightFloatCorrect.GetNumberOfValues());
      std::cout << "(Intrinsic Weight - Size) = " << superarcIntrinsicWeightPortal.GetNumberOfValues() << std::endl;
//      PrintIndices("Intrinsic Weight", superarcIntrinsicWeight);
      PrintValues("Intrinsic Weight", superarcIntrinsicWeight);
      std::cout << "(Dependent Weight - Size) = " << superarcDependentWeight.GetNumberOfValues() << std::endl;
//      PrintIndices("Dependent Weight", superarcDependentWeight);
      PrintValues("Dependent Weight", superarcDependentWeight);
//      PrintIndices("Upwards Weight",           upWeight);
      std::cout << "(Upwards Weight - Size) = " << upWeightFloatCorrect.GetNumberOfValues() << std::endl;
      PrintValues("Upwards Weight (float)",   upWeightFloat); // might be a corrupt 'correct' array ...
//      PrintIndices("Downwards Weight",         downWeight);
      std::cout << "(Downwards Weight - Size) = " << downWeightFloatCorrect.GetNumberOfValues() << std::endl;
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
      std::cout << "Unsorted arcs: (" << nSuperarcs << ")" << std::endl;
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
//        process_contourtree_inc_ns::SuperArcVolumetricComparator(upWeightFloatCorrect, superarcList, false)); // 2025-12-28 crash
        process_contourtree_inc_ns::SuperArcVolumetricComparator(upWeightFloat, superarcList, false)); // using upWeightFloat

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

      // crash

    } // ComputeVolumeBranchDecompositionSerialFloat()

    //jump














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
    std::cout << "II A) Weights Computed" << std::endl;
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
    for (vtkm::Id sortedNode = 0; sortedNode < contourTree.Arcs.GetNumberOfValues(); sortedNode++)
    { // per node in sorted order
      vtkm::Id sortID = nodesPortal.Get(sortedNode);
      vtkm::Id superparent = superparentsPortal.Get(sortID);
      if (sortedNode == 0)
        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
      else if (superparent != superparentsPortal.Get(nodesPortal.Get(sortedNode - 1)))
        firstVertexForSuperparentPortal.Set(superparent, sortedNode);
    } // per node in sorted order
    // now we use that to compute the intrinsic weights
    for (vtkm::Id superarc = 0; superarc < contourTree.Superarcs.GetNumberOfValues(); superarc++)
      if (superarc == contourTree.Superarcs.GetNumberOfValues() - 1)
        superarcIntrinsicWeightPortal.Set(superarc,
                                          contourTree.Arcs.GetNumberOfValues() -
                                            firstVertexForSuperparentPortal.Get(superarc));
      else
        superarcIntrinsicWeightPortal.Set(superarc,
                                          firstVertexForSuperparentPortal.Get(superarc + 1) -
                                            firstVertexForSuperparentPortal.Get(superarc));

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

      // so, step 1: add xfer + int & store in dependent weight
      for (vtkm::Id supernode = firstSupernode; supernode < lastSupernode; supernode++)
      {
        superarcDependentWeightPortal.Set(supernode,
                                          supernodeTransferWeightPortal.Get(supernode) +
                                            superarcIntrinsicWeightPortal.Get(supernode));
      }

      // step 2: perform prefix sum on the dependent weight range
      for (vtkm::Id supernode = firstSupernode + 1; supernode < lastSupernode; supernode++)
        superarcDependentWeightPortal.Set(supernode,
                                          superarcDependentWeightPortal.Get(supernode) +
                                            superarcDependentWeightPortal.Get(supernode - 1));

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
    }   // per iteration
  }     // ContourTreeMaker::ComputeWeights()


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
          const vtkm::cont::ArrayHandle<T, StorageType>& valueField,
    const vtkm::cont::ArrayHandle<T, StorageType>& dataField,
    bool dataFieldIsSorted,
    const FloatArrayType& superarcDependentWeight,            // NEW: passed intrincid
    const FloatArrayType& superarcIntrinsicWeight,
          const vtkm::Id& contourTreeRootnode)                // NEW: used to get the augmented betti nodes (which are past the root node in index)
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
                valueField,
      dataField,
      dataFieldIsSorted,
      superarcDependentWeight,
      superarcIntrinsicWeight,
                contourTreeRootnode);
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
