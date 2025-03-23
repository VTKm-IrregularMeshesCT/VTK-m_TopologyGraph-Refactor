//============================================================================
// Copyright (c) 2019, The Regents of the University of California, through
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
//============================================================================

#include "interface.h"

#include "./ct.h"
#include "./mc.h"
#include "./filter.h"

#include <string>
#include <cassert>
#include <algorithm>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKUnstructuredGridReader.h>


#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleConstant.h>
//#include <vtkm/cont/DataSetFieldAdd.h> deprecated since VTK 1.6 https://gitlab.kitware.com/vtk/vtk-m/-/releases/v1.6.0
#include <vtkm/cont/DataSetBuilderExplicit.h>

#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/Types.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/ProcessContourTree.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/Branch.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/processcontourtree/SetTriangleSuperarcId.h>


//#include <vtkm/io/VTKUnstructuredGridReader.h>
//#include <vtkm/filter/scalar_topology/worklet/ContourTreeUniformAugmented.h>
#include <vtkm/filter/scalar_topology/worklet/contourtree_augmented/meshtypes/ContourTreeMesh.h>



struct DelaunayMesh
{
    std::vector<vtkm::Id> std_nbor_connectivity;
    std::vector<vtkm::Id> std_nbor_offsets;
};


DelaunayMesh parseDelaunayASCII(const std::string& filePath)
{
    std::ifstream inputFile(filePath);
    std::string line;
    DelaunayMesh graph;
    size_t currentOffset = 0;

    if(!inputFile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return graph; // This will return empty vectors if file cannot be opened
    }

    while(std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        vtkm::Id id;

        // Push the currentOffset before processing the line
        graph.std_nbor_offsets.push_back(currentOffset);

        // Process each id in the current line
        while (iss >> id) {
            graph.std_nbor_connectivity.push_back(id);
            ++currentOffset; // Increment the offset for each id found
        }
    }

    // After processing all lines, the last offset should be equal to the total number of ids
    // This implies the end of the last vertex's neighborhood
    graph.std_nbor_offsets.push_back(currentOffset);

    return graph;

}


using namespace std;
using namespace vtkm;

using Coefficients = vtkm::worklet::contourtree_augmented::Coefficients;
using FloatArrayType = vtkm::cont::ArrayHandle<vtkm::Float64>;
namespace ctaug_ns = vtkm::worklet::contourtree_augmented;

vtkm::cont::PartitionedDataSet cv1k::interface::computeMostSignificantContours(vtkm::cont::DataSet inputData,
                                                                               std::string fieldName,
                                                                               std::string inputCTFilename,
                                                                               int simplificationThreshold,
                                                                               std::string decompositionType,
                                                                               string selectionType,
                                                                               std::string branchIsovalueFlag,
                                                                               bool performanceRun)
{//computeMostSignificantContours
    // debug >= 1, print what is being currently computed with timings
    // debug >= 2, print contour tree arrays and isosurfac triangles

    std::cout << "computeMostSignificantContours() " << std::endl;

    int debugLevel = 1;

    //
    // Compute Contour Tree
    //
    if (debugLevel >= 1)
    {
        //cout << "Computing Contour Tree ..." << endl;
    }

    vtkm::cont::Timer timer;
    timer.Start();

    worklet::contourtree_augmented::ContourTree ct;
    vtkm::cont::ArrayHandle<vtkm::Id> ctSortOrder;
    vtkm::cont::ArrayHandle<vtkm::Id> ctSortIndices;
    Id ctNumIterations;

    if (true == inputCTFilename.empty())
    {
        std::cout << "No ContourTree input, generating it from scratch ..." << std::endl;
        tie(ct, ctSortOrder, ctSortIndices, ctNumIterations) = cv1k::ct::getContourTree(inputData, fieldName);
    }
    else
    {
        std::cout << "Reading in the ContourTree from file: " << inputCTFilename << std::endl;
        tie(ct, ctSortOrder, ctNumIterations) = cv1k::ct::readContourTree(inputCTFilename);
    }

    timer.Stop();

    if (debugLevel >= 1)
    {
        //std::cout << "Computing CT took " << timer.GetElapsedTime() << " seconds." << std::endl;
    }

    //
    // Compute Branch Decomposition
    //

    // ============================================== Traditional CT ============================================== //

    cont::ArrayHandle<vtkm::Id> whichBranch;
    cont::ArrayHandle<vtkm::Id> branchMinimum;
    cont::ArrayHandle<vtkm::Id> branchMaximum;
    cont::ArrayHandle<vtkm::Id> branchSaddle;
    cont::ArrayHandle<vtkm::Id> branchParent;

    cont::ArrayHandle<Id> superarcIntrinsicWeight;
    cont::ArrayHandle<Id> superarcDependentWeight;
    cont::ArrayHandle<Id> supernodeTransferWeight;
    cont::ArrayHandle<Id> hyperarcDependentWeigh;

    // ================================================= PACT-BD ================================================= //

    // Floating point type weights are required for PACT-BD
    FloatArrayType superarcIntrinsicWeightNEW;
    FloatArrayType superarcDependentWeightNEW;
    FloatArrayType supernodeTransferWeightNEW;
    FloatArrayType hyperarcDependentWeightNEW;

    // compute the branch decomposition by volume (these remain the same for PACT-BD)
    // The following arrays are already defined
    //        ctaug_ns::IdArrayType whichBranch;
    //        ctaug_ns::IdArrayType branchMinimum;
    //        ctaug_ns::IdArrayType branchMaximum;
    //        ctaug_ns::IdArrayType branchSaddle;
    //        ctaug_ns::IdArrayType branchParent;

    if ("volume" == decompositionType)
    {
        std::cout << "Branch Weights Chosen: VOLUME" << std::endl;
        worklet::contourtree_augmented::ProcessContourTree::ComputeVolumeWeightsSerial(
                ct,
                ctNumIterations,
                superarcIntrinsicWeight,
                superarcDependentWeight,
                supernodeTransferWeight,
                hyperarcDependentWeigh
                );

        worklet::contourtree_augmented::ProcessContourTree::ComputeVolumeBranchDecompositionSerial(
                ct,
                superarcDependentWeight,
                superarcIntrinsicWeight,
                whichBranch,
                branchMinimum,
                branchMaximum,
                branchSaddle,
                branchParent
                );

    }
    else if ("pactbd" == decompositionType)
    {
        std::cout << "Branch Weights Chosen: PACTBD" << std::endl;

        vtkm::cont::ArrayHandle<Coefficients> superarcIntrinsicWeightCoeffs;
        vtkm::cont::ArrayHandle<Coefficients> superarcDependentWeightCoeffs;
        vtkm::cont::ArrayHandle<Coefficients> supernodeTransferWeightCoeffs;
        vtkm::cont::ArrayHandle<Coefficients> hyperarcDependentWeightCoeffs;

        ctaug_ns::ProcessContourTree::ComputeVolumeWeightsSerialStructCoefficients(ct,
                                                                                  ctNumIterations,
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


        ctaug_ns::ProcessContourTree::ComputeVolumeBranchDecompositionSerialFloat(ct,
                                                                                  superarcDependentWeightNEW,
                                                                                  superarcIntrinsicWeightNEW,
                                                                                  whichBranch,   // (output)
                                                                                  branchMinimum, // (output)
                                                                                  branchMaximum, // (output)
                                                                                  branchSaddle,  // (output)
                                                                                  branchParent); // (output)
    }

    else if("height" == decompositionType)
    {
        std::cout << "Branch Weights Chosen: HEIGHT" << std::endl;
        vtkm::cont::Timer timer;
        timer.Start();

        worklet::contourtree_augmented::ProcessContourTree::ComputeHeightBranchDecomposition(
                ct,
                inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>(),
                ctSortOrder,
                ctNumIterations,
                whichBranch,
                branchMinimum,
                branchMaximum,
                branchSaddle,
                branchParent
                );
        timer.Stop();
        //std::cout << "Height branch decomposition took " << timer.GetElapsedTime() << " seconds." << std::endl;
        if (true == performanceRun)
        {
            std::cout << timer.GetElapsedTime() << std::endl;
        }
    }
    else
    {
        std::cout << "Branch Weights Chosen: NOT IMPLEMENTED" << std::endl;
    }

    // Bit hacky, but no need to do anything else if this is a performance run
    if (true == performanceRun)
    {
        return vtkm::cont::PartitionedDataSet();
    }

    //
    // Compute Additioanl Data for the Branches
    //

    // Which branch do regular vertices belong to
    vector<Id> whichBranchRegular(ct.Nodes.GetNumberOfValues());

    // Endpoints of path in the contour tree
    cont::ArrayHandle<Vec<Id, 2>> branchEndpoints;
    branchEndpoints.Allocate(branchMinimum.GetNumberOfValues());

    cont::ArrayHandle<Vec<Id, 2>> branchEndpointsRegular;
    branchEndpointsRegular.Allocate(branchMinimum.GetNumberOfValues());

    // Endpoints of path in the contour tree
    cont::ArrayHandle<Float64> branchIsovalueHeightArray;
    branchIsovalueHeightArray.Allocate(branchMinimum.GetNumberOfValues());

    // Each point has it's own isovalue
    cont::ArrayHandle<Float64> branchIsovalueArray;
    branchIsovalueArray.Allocate(branchMinimum.GetNumberOfValues());

    // Each point has it's own isovalue
    cont::ArrayHandle<Float64> branchPointVolumeArray;
    branchPointVolumeArray.Allocate(branchMinimum.GetNumberOfValues());

    // Endpoints of path in the contour tree
    cont::ArrayHandle<Id> branchHeightArray;
    branchHeightArray.Allocate(branchMinimum.GetNumberOfValues());


    // Populate the previously innitialsed arrays
    if ("volume" == decompositionType)
    {
        computeAdditionalBranchData(
                inputData,
                fieldName,
                ct,
                ctSortOrder,
                whichBranch,
                branchMinimum,
                branchSaddle,
                branchMaximum,
                whichBranchRegular,
                branchEndpointsRegular,
                branchEndpoints,
                branchHeightArray,
                branchIsovalueHeightArray,
                branchIsovalueArray,
                branchPointVolumeArray,
                superarcIntrinsicWeight,
                superarcDependentWeight,
                branchIsovalueFlag
                );
    }
    else if ("pactbd" == decompositionType)
    {
        // PACTBD-EDIT
        int num_datapoints = 10001;

        const std::string field_filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/10k-field.txt";
        std::ifstream field_input(field_filename);
        vtkm::cont::ArrayHandle<Float64> fakeFieldArray;
        fakeFieldArray.Allocate(num_datapoints);
        auto fakeFieldArrayWritePortal = fakeFieldArray.WritePortal();

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
        computeAdditionalBranchDataFloat(fakeFieldArray,
                        //                inputData,
                        //                fieldName,
                                        ct,
                                        ctSortOrder,
                                        whichBranch,
                                        branchMinimum,
                                        branchSaddle,
                                        branchMaximum,
                                        whichBranchRegular,
                                        branchEndpointsRegular,
                                        branchEndpoints,
                                        branchHeightArray,
                                        branchIsovalueHeightArray,
                                        branchIsovalueArray,
                                        branchPointVolumeArray,
                                        superarcIntrinsicWeightNEW,
                                        superarcDependentWeightNEW,
                                        branchIsovalueFlag);
    }

    // Prevent getting more branches than we have available
    int numberOfBranches = min(static_cast<unsigned long long>(simplificationThreshold), static_cast<unsigned long long>(branchMaximum.GetNumberOfValues()));

    std::cout << "Extracting " << numberOfBranches << " most significant contours" << std::endl;

    // 
    // We also need the mesh and extrema for the extraction of IDs (for marching cubes). This only works for 3D data.
    //
    if (debugLevel >= 1)
    {
        cout << "Mesh and extrema ..." << endl;
    }
    timer.Reset();
    timer.Start();


    // Now we have to recreate some components that are used in the Contour Tree computation ...
    // ... that are deleted afterwards for saving memory such as the mesh, and the extrema chains

    cout << "Reating with VTKUnstructuredGridReader" << endl;
    // Manually sending the VTK file from here since it is only needed here
//    vtkm::io::VTKUnstructuredGridReader reader("../delaunay-parcels/10k-from-2M-sampled-excel-sorted.1.vtk");
    //vtkm::io::VTKDataSetReader reader("../delaunay-parcels/10k-from-2M-sampled-excel-sorted.1.vtk");
    vtkm::io::VTKDataSetReader reader("../delaunay-parcels/10k-from-2M-sampled-excel-sorted-withvalues.vtk");

    reader.PrintSummary(std::cout);
    cont::DataSet inputDataVTK = reader.ReadDataSet();

//    vtkm::io::VTKUnstructuredGridReader reader("/home/sc17dd/modules/HCTC2024/VTK-m-topology-refactor/VTK-m_TopologyGraph-Refactor/examples/contour-visualiser/build/10k-from-2M-sampled-excel-sorted.1.vtk");
////    vtkm::io::VTKUnstructuredGridReader reader("../delaunay-parcels/10k-from-2M-sampled-excel-sorted-field.vtk");
//    reader.PrintSummary(std::cout);
//    cont::DataSet inputDataVTK = reader.ReadDataSet();

    cout << "Done!" << endl;

    reader.PrintSummary(std::cout);

    cout << "Summary done!" << endl;


//    // Build the mesh (Regular)
//        vtkm::Id3 pointDimensions = inputData.GetCellSet().AsCellSet<vtkm::cont::CellSetStructured<3>>().GetPointDimensions();
//        worklet::contourtree_augmented::DataSetMeshTriangulation3DMarchingCubes mesh(vtkm::Id3(pointDimensions[0], pointDimensions[1], pointDimensions[2]));
//        mesh.SortData(inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>());

    // Build the mesh (PACT-BD)
    // PACTBD-EDIT
    int num_datapoints = 10001;
    const std::string filename = "/home/sc17dd/modules/HCTC2024/VTK-m-topology/vtkm-build/10k-from-2M-sampled-excel-sorted.1-CONNECTIVITY.txt";
    DelaunayMesh delmesh = parseDelaunayASCII(filename);
    // Get CONNECTIVITY
    vtkm::cont::ArrayHandle<vtkm::Id> nbor_connectivity =
      vtkm::cont::make_ArrayHandle(delmesh.std_nbor_connectivity, vtkm::CopyFlag::Off);
    std::cout << "nbor_connectivity num vals: " << nbor_connectivity.GetNumberOfValues() << "\n";
    // Get OFFSETS
    vtkm::cont::ArrayHandle<vtkm::Id> nbor_offsets =
      vtkm::cont::make_ArrayHandle(delmesh.std_nbor_offsets, vtkm::CopyFlag::Off);
    std::cout << "nbor_offsets num vals: " << nbor_offsets.GetNumberOfValues() << "\n";
    // Get VALUES
    //    std::vector<int> std_actual_values;
    std::vector<vtkm::Float64> std_actual_values;
    for(int i = 0; i < num_datapoints; i++)
    {
      std_actual_values.push_back((vtkm::Float64)i);
    }
   // note: this array does not actually matter for our case ...
//    vtkm::cont::ArrayHandle<int> actual_values =
    vtkm::cont::ArrayHandle<vtkm::Float64> actual_values =
      vtkm::cont::make_ArrayHandle(std_actual_values, vtkm::CopyFlag::Off);
    // Get GLOBAL IDs
    std::vector<vtkm::Id> std_global_inds = {0};
    vtkm::cont::ArrayHandle<vtkm::Id> global_inds =
      vtkm::cont::make_ArrayHandle(std_global_inds, vtkm::CopyFlag::Off);

    std::cout << "USING PACT (unoptimized)...\n";
//    vtkm::worklet::contourtree_augmented::ContourTreeMesh<int> mesh(ctSortOrder,
    vtkm::worklet::contourtree_augmented::ContourTreeMesh<vtkm::Float64> mesh(ctSortOrder, // const IdArrayType& nodes,
                            //arcs_list,
                              nbor_connectivity,                                           // const IdArrayType& inNborConnectivity
                              nbor_offsets,                                                // const IdArrayType& inNborOffsets
                              ctSortOrder,                                                 // const IdArrayType& inSortOrder
                              // doesnt work out of the box:
                              // fieldArray, // testing fieldArray instead of manual actual_value
                              actual_values,                                               // const vtkm::cont::ArrayHandle<FieldType>& values
                              //nodes_sorted,
                              global_inds);                                                // const IdArrayType& inGlobalMeshIndex


    // Set up extremal chains (unchanged)
    worklet::contourtree_augmented::MeshExtrema extrema(mesh.NumVertices);
    extrema.SetStarts(mesh, true);
    extrema.BuildRegularChains(true);
    extrema.SetStarts(mesh, false);
    extrema.BuildRegularChains(false);

    timer.Stop();

    if (debugLevel >= 1)
    {
        std::cout << "Computing Mesh and Extrema took " << timer.GetElapsedTime() << " seconds." << std::endl;
        cout << "Branch Decomposition ..." << endl;
    }


    vtkm::cont::ArrayHandle<Float64> regularValueSecondaryImportance;
    regularValueSecondaryImportance.Allocate(branchMinimum.GetNumberOfValues());

    for (int i = 0 ; i < regularValueSecondaryImportance.GetNumberOfValues() ; i++)
    {
        // We've set B to be the bigger one in computeAdditionalBranchData
        //vtkm::Id endpointBRegularID = static_cast<vtkm::Float64>(ctSortIndices.ReadPortal().Get(branchEndpoints.ReadPortal().Get(i)[1]));
        vtkm::Id endpointBRegularID = 11;
        regularValueSecondaryImportance.WritePortal().Set(i, endpointBRegularID);
    }


    // 
    // Pick the importance metric. That is either height or volume (which is independent from the simplification type)
    //
    vtkm::cont::ArrayHandle<Float64> branchImportance;
    vtkm::cont::ArrayHandle<Float64> branchImportanceSecondary;
    if ("volume" == decompositionType)
    {
        branchImportance = branchPointVolumeArray;
        branchImportanceSecondary = branchIsovalueHeightArray;
    }
    else if ("pactbd" == decompositionType)
    {
        branchImportance = branchPointVolumeArray;
        branchImportanceSecondary = branchIsovalueHeightArray;
    }
    else if("height" == decompositionType)
    {
        branchImportance = branchIsovalueHeightArray;
        branchImportanceSecondary = regularValueSecondaryImportance;
    }
    else
    {
        branchImportance = branchIsovalueHeightArray;
        branchImportanceSecondary = regularValueSecondaryImportance;
    }



    // Get branch decomposition & simplify and get relevant isovalues.
    vector<std::pair<vtkm::Float64, vtkm::Id>> vals;

    if ("sort" == selectionType)
    {
        vector<Id> branchOrder = getBranchesSortedOrder(branchImportance, branchImportanceSecondary);

        for (int i = 0 ; i  < branchOrder.size() ; i++)
        {
            vals.push_back({branchIsovalueArray.ReadPortal().Get(branchOrder[i]), branchOrder[i]});
        }
    }
    else
    {
        assert(false);
    }

    //
    // Extract one contour per branch
    //

    // We could have less branches than we specified.
    numberOfBranches = std::min(static_cast<unsigned long long>(numberOfBranches), static_cast<unsigned long long>(vals.size()));

    // Prepare output

    vtkm::cont::PartitionedDataSet outputContours;


    // ==================================================================================================================================== //
    // ==================================================================================================================================== //
    // ==================================================================================================================================== //

    // SET UP (creating the CT, branch decomposition, simplification, and choosing N most important) DONE. Now actualy extracting contours  //

    // ==================================================================================================================================== //
    // ==================================================================================================================================== //
    // ==================================================================================================================================== //

    for (int k = 0 ; k < numberOfBranches; k++)
    {
        //
        // Compute MC Triangles for each isovalue and filter out the triangles that do not belong to the current branch
        //
        const Id branchId = vals[k].second;
        const Float64 branchIsovalue = vals[k].first;

        // 
        // Print some debug info
        //
        if (debugLevel >= 1)
        {
            cout << endl << endl << "=======================================================" << std::endl;
            printf("Getting the triangles for isovalue %f, index %i and branch ID %llu and importance %f and height %f volume %f.\n",
                                           branchIsovalue,        k,          branchId,
                                                            branchImportance.ReadPortal().Get(branchId),
                                                                  branchIsovalueHeightArray.ReadPortal().Get(branchId),
                                                                              branchPointVolumeArray.ReadPortal().Get(branchId));
        }

        // Compute an isosurface for the whole data set
//        cont::ArrayHandle<cv1k::Triangle> mcTriangles = cv1k::mc::getMarchingCubeTriangles(inputData, {branchIsovalue}, fieldName);
        cont::ArrayHandle<cv1k::Triangle> mcTriangles = cv1k::mc::getMarchingCubeTriangles(inputDataVTK, {branchIsovalue}, fieldName);

        std::cout << "mcTriangles extracted total:" << mcTriangles.GetNumberOfValues() << std::endl;

        // Compute the superarc ID of all the triangles
//         cv1k::filter::computeTriangleIds(ct, mesh, extrema, inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>(), mcTriangles, branchIsovalue);
        cv1k::filter::computeTriangleIds(ct, mesh, extrema,
                                         inputDataVTK.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<vtkm::Float64>>(),
                                         mcTriangles, branchIsovalue);
        std::cout << "Triangle IDs computed" << std::endl;
        std::cout << "mcTriangles remaining total:" << mcTriangles.GetNumberOfValues() << std::endl;

        //
        // Compute the superarc from the current branch sits at that isovalue 
        //
        cont::ArrayHandle<Vec<Id, 2>> endpoints;
        cont::ArrayHandle<Id> superarcIds;
        cont::ArrayHandle<Float64> isovalueArray;
        endpoints.Allocate(1);
        superarcIds.Allocate(1);
        isovalueArray.Allocate(1);


        isovalueArray.WritePortal().Set(0, branchIsovalue);
        //endpoints.WritePortal().Set(0, branchEndpoints.ReadPortal().Get(branchId));
        endpoints.WritePortal().Set(0, {
                ctSortOrder.ReadPortal().Get(branchEndpointsRegular.ReadPortal().Get(branchId)[0]),
                ctSortOrder.ReadPortal().Get(branchEndpointsRegular.ReadPortal().Get(branchId)[1]),
                });

        // Set up the worklet
        vtkm::worklet::contourtree_augmented::process_contourtree_inc::SetTriangleSuperarcId setTrianglesId(ct.Hypernodes.GetNumberOfValues(),
                                                                                                            ct.Supernodes.GetNumberOfValues());
        cont::Invoker Invoke;

        //// Run the worklet
        Invoke(
                setTrianglesId,
                endpoints,
                    inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>(),
//                    inputDataVTK.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<vtkm::Float64>>(),
                isovalueArray,
                mesh.SortOrder, // (input)
                mesh.SortIndices, // (input)
                ct.Superparents, // (input)
                ct.WhenTransferred, // (input)
                ct.Hyperparents, // (input)
                ct.Hyperarcs, // (input)
                ct.Hypernodes, // (input)
                ct.Supernodes, // (input)
                extrema.Peaks, // (input)
                extrema.Pits,
                superarcIds
              ); // (input)

        //const vtkm::Id branchSuperarcID = superarcIds.ReadPortal().Get(0);
        // assign the branch name to the last superarc, not the first one for PACTBD somehow:
        const vtkm::Id branchSuperarcID = superarcIds.ReadPortal().Get(superarcIds.GetNumberOfValues()-1);

        std::cout << "branchSuperarcID: " << branchSuperarcID << std::endl;


        //
        // Filter ouf the triangles that do not belong to the isosurface of the current branch
        //
        vector<cv1k::Triangle> branchTriangles;
        for (int j = 0 ; j < mcTriangles.GetNumberOfValues() ; j++)
        {
            Triangle currentTriangle = mcTriangles.ReadPortal().Get(j);
//            branchTriangles.push_back(currentTriangle);

            if (currentTriangle.superarcId == branchSuperarcID)
            {
                currentTriangle.superarcId = k;
                currentTriangle.isovalue = branchIsovalue;
                currentTriangle.importance = branchImportance.ReadPortal().Get(branchId);
                branchTriangles.push_back(currentTriangle);
            }
        }

        std::cout << "branchTriangles remaining total:" << branchTriangles.size() << std::endl;


        if (debugLevel >= 1)
        {
            cout << "The size of the contour is " << (branchTriangles.size()) << " from " << mcTriangles.GetNumberOfValues()<< " that is " <<  100.0 * (static_cast<double>(branchTriangles.size()) / static_cast<double>(mcTriangles.GetNumberOfValues())) << "%" << endl;
        }

        std::vector<vtkm::Vec3f_32> pointCoordinates;
        std::vector<vtkm::UInt8> shapes;
        std::vector<vtkm::IdComponent> numIndices;
        std::vector<vtkm::Id> connectivity;

        for (size_t i = 0 ; i < branchTriangles.size() ; i++)
        {
            for (int j = 0 ; j < 3 ; j++)
            {
                pointCoordinates.push_back({branchTriangles[i].points[j][0], branchTriangles[i].points[j][1], branchTriangles[i].points[j][2]});
            }
            shapes.push_back(vtkm::CELL_SHAPE_TRIANGLE);
            numIndices.push_back(3);
            connectivity.push_back(i * 3);
            connectivity.push_back(i * 3 + 1);
            connectivity.push_back(i * 3 + 2);
        }



        vtkm::cont::DataSetBuilderExplicit  dataSetBuilder;
        vtkm::cont::DataSet contourDataSet = dataSetBuilder.Create(pointCoordinates, shapes, numIndices, connectivity);


        // ----------------- Write Cell Data ----------------- //

        cont::ArrayHandle<int> branchIDCellField;
        branchIDCellField.Allocate(contourDataSet.GetNumberOfCells());
        auto branchIDCellFieldWritePortal = branchIDCellField.WritePortal();
        for (int i = 0 ; i < branchIDCellField.GetNumberOfValues() ; i++)
        {
            branchIDCellFieldWritePortal.Set(i, branchId);
        }
        contourDataSet.AddCellField("branchId", branchIDCellField);

        cont::ArrayHandle<int> importanceCellField;
        importanceCellField.Allocate(contourDataSet.GetNumberOfCells());
        auto importanceCellFieldWritePortal = importanceCellField.WritePortal();
        for (int i = 0 ; i < importanceCellField.GetNumberOfValues() ; i++)
        {
            importanceCellFieldWritePortal.Set(i, k);
        }
        contourDataSet.AddCellField("importance", importanceCellField);

        std::cout << "4 branchIsovalue:" << branchIsovalue << std::endl;
        std::cout << "4 branchId      :" << branchId << std::endl;

        cont::ArrayHandle<vtkm::Float64> isovalueCellField;
        isovalueCellField.Allocate(contourDataSet.GetNumberOfCells());
        auto isovalueCellFieldWritePortal = isovalueCellField.WritePortal();
        for (int i = 0 ; i < isovalueCellField.GetNumberOfValues() ; i++)
        {
            isovalueCellFieldWritePortal.Set(i, branchIsovalue);
        }
        contourDataSet.AddCellField("branchIsovalue",    isovalueCellField);
        std::cout << "5 branchIsovalue:" << branchIsovalue << std::endl;

//        if ("volume" == decompositionType)
//        {
//            cont::ArrayHandle<int> volumeCellField;
//            volumeCellField.Allocate(contourDataSet.GetNumberOfCells());
//            auto volumeCellFieldWritePortal = volumeCellField.WritePortal();
//            for (int i = 0 ; i < volumeCellField.GetNumberOfValues() ; i++)
//            {
//                volumeCellFieldWritePortal.Set(i, branchPointVolumeArray.ReadPortal().Get(branchId));
//            }
//            contourDataSet.AddCellField("branchVolume",volumeCellField);
//        }

//        else if ("pactbd" == decompositionType)
//        {
            cont::ArrayHandle<vtkm::Float64> volumeCellField;
            volumeCellField.Allocate(contourDataSet.GetNumberOfCells());
            auto volumeCellFieldWritePortal = volumeCellField.WritePortal();
            for (int i = 0 ; i < volumeCellField.GetNumberOfValues() ; i++)
            {
                volumeCellFieldWritePortal.Set(i, branchPointVolumeArray.ReadPortal().Get(branchId));
            }
            contourDataSet.AddCellField("branchVolume",volumeCellField);
//        }


        cont::ArrayHandle<int> superarcIdCellField;
        superarcIdCellField.Allocate(mcTriangles.GetNumberOfValues());
        auto superarcIdCellFieldWritePortal = superarcIdCellField.WritePortal();
        for (int j = 0 ; j < mcTriangles.GetNumberOfValues() ; j++)
        {
            Triangle currentTriangle = mcTriangles.ReadPortal().Get(j);
            superarcIdCellFieldWritePortal.Set(j, currentTriangle.superarcId);
        }
        contourDataSet.AddCellField("superarcID",superarcIdCellField);




        // ----------------- Write Point Data (shown first) ----------------- //

        cont::ArrayHandle<vtkm::Float64> branchIDPointField;
        branchIDPointField.Allocate(contourDataSet.GetNumberOfPoints());
        auto branchIDPointFieldWritePortal = branchIDPointField.WritePortal();
        for (int i = 0 ; i < branchIDPointField.GetNumberOfValues() ; i++)
        {
            branchIDPointFieldWritePortal.Set(i, branchId);//0.5);
        }

        contourDataSet.AddPointField("branchIDpt", branchIDPointField);

        cont::ArrayHandle<vtkm::Float64> isovaluePointField;
        isovaluePointField.Allocate(contourDataSet.GetNumberOfPoints());
        auto isovaluePointFieldWritePortal = isovaluePointField.WritePortal();
        for (int i = 0 ; i < isovaluePointField.GetNumberOfValues() ; i++)
        {
            isovaluePointFieldWritePortal.Set(i, branchIsovalue);//0.5);
        }

        contourDataSet.AddPointField("isovaluePoints", isovaluePointField);


        //        contourDataSet.AddCellField("branchImportance",  branchImportance.ReadPortal().Get(branchId));
        //        contourDataSet.AddCellField("branchHeight",      branchIsovalueHeightArray.ReadPortal().Get(branchId));
        outputContours.AppendPartition(contourDataSet);
    }

    return outputContours;
}//computeMostSignificantContours


vector<Id> cv1k::interface::getBranchesSortedOrder(cont::ArrayHandle<Float64> branchImportance,
                                                   cont::ArrayHandle<Float64> secondaryBranchImportance)
{
    // Add only the branches which have non zero height
    vector<tuple<Id, Float64, Float64>> sortedBranches;
    for (int i = 0 ; i < branchImportance.GetNumberOfValues() ; i++)
    {
        if (branchImportance.ReadPortal().Get(i) > 0 && secondaryBranchImportance.ReadPortal().Get(i) > 0)
        {
            sortedBranches.push_back({i, branchImportance.ReadPortal().Get(i), secondaryBranchImportance.ReadPortal().Get(i)});
        }
        else
        {
        }
    }

    // Sort by first criteria, than volume, when the first criteria is volume this is redundant.
    std::sort(sortedBranches.begin(), sortedBranches.end(), [](const tuple<Id, Float64, Float64> a, const tuple<Id, Float64, Float64> b){ 
            if (std::get<1>(a) == std::get<1>(b))
            {
                return std::get<2>(a) > std::get<2>(b);
            }
            return std::get<1>(a) > std::get<1>(b);
            });

    vector<Id> sortedBranchesMap(sortedBranches.size());
    for (int i = 0 ; i < sortedBranchesMap.size() ; i++)
    {
        sortedBranchesMap[i] = std::get<0>(sortedBranches[i]);
    }

    return sortedBranchesMap;
}


void cv1k::interface::computeAdditionalBranchData(
        const cont::DataSet inputData,
        const string fieldName,
        const worklet::contourtree_augmented::ContourTree ct, 
        const vtkm::cont::ArrayHandle<vtkm::Id> ctSortOrder,
        const cont::ArrayHandle<Id> whichBranch, 
        const cont::ArrayHandle<Id> branchMinimum, 
        const cont::ArrayHandle<Id> branchSaddle, 
        const cont::ArrayHandle<Id> branchMaximum, 
        vector<Id> &whichBranchRegular,
        cont::ArrayHandle<Vec<Id, 2>> &branchEndpointsRegular,
        cont::ArrayHandle<Vec<Id, 2>> &branchEndpoints,
        cont::ArrayHandle<Id> &branchHeightArray,
        cont::ArrayHandle<Float64> &branchIsovalueHeightArray,
        cont::ArrayHandle<Float64> &branchIsovalueArray,
        cont::ArrayHandle<Float64> &branchPointVolumeArray,
        const cont::ArrayHandle<Id> superarcIntrinsicWeight, 
        const cont::ArrayHandle<Id> superarcDependentWeight,
        const std::string branchIsovalueFlag
        )
{
    // Compute the volume of every branch.
    for (int i = 0 ; i < branchPointVolumeArray.GetNumberOfValues() ; i++)
    {
        branchPointVolumeArray.WritePortal().Set(i, 1);
    }

    using vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT;
    // @TODO Should I use intrinsic of dependent weight?
    std::cout << "SuperArc to Branch ID mappings (num. mappings: "
              << superarcIntrinsicWeight.GetNumberOfValues() << "=num. of superarcs):" << std::endl;
    for (int i = 0 ; i < superarcIntrinsicWeight.GetNumberOfValues() ; i++)
    {
        Id branchId = whichBranch.ReadPortal().Get(i);

#ifdef DEBUG_PRINT
        std::cout << "SA " << i << " -> BID: " << branchId << std::endl;
#endif

        if (branchId != NO_SUCH_ELEMENT)
        {
            branchPointVolumeArray.WritePortal().Set(branchId, 1);
            Id currentWeight = branchPointVolumeArray.ReadPortal().Get(branchId);
            //Id weight = superarcIntrinsicWeight.ReadPortal().Get(i);
            Id weight = superarcDependentWeight.ReadPortal().Get(i);

            branchPointVolumeArray.WritePortal().Set(branchId, currentWeight + weight);
        }
    }

    std::vector<std::tuple<vtkm::Id, vtkm::Id, vtkm::Id, vtkm::Float64>> branchEndpointsStd;

    // Converts the branches to mesh vertex endpoints and compute isovalue
    std::cout << "Computing isovalues for each branch (total num. isovalues: "
              << branchMinimum.GetNumberOfValues() << "=num. of branches):" << std::endl;

    for (int i = 0 ; i < branchMinimum.GetNumberOfValues() ; i++)
    {
        using vtkm::worklet::contourtree_augmented::MaskedIndex;

        Id branchHeight = 0;
        Vec<Id, 2> endpoints;
        Vec<Id, 2> regularEndpoints;

        Id max = MaskedIndex(ct.Supernodes.ReadPortal().Get(MaskedIndex(branchMaximum.ReadPortal().Get(i))));
        Id min = MaskedIndex(ct.Supernodes.ReadPortal().Get(MaskedIndex(branchMinimum.ReadPortal().Get(i))));
        Id saddle = MaskedIndex(ct.Supernodes.ReadPortal().Get(MaskedIndex(branchSaddle.ReadPortal().Get(i))));

        // 0 for descending, 1 for ascending, 2 for Master,
        int branchType = -1;

        // Determine the type of the branch and set its endpoints
        using vtkm::worklet::contourtree_augmented::NoSuchElement;
        if (false == NoSuchElement(branchSaddle.ReadPortal().Get(static_cast<vtkm::Id>(i))))
        {
            if (min < saddle)
            {
                endpoints = {ctSortOrder.ReadPortal().Get(min), ctSortOrder.ReadPortal().Get(saddle)};
                regularEndpoints = {min, saddle};
                branchHeight = saddle - min;
                branchType = 0;
            }
            else if (saddle < max)
            {
                endpoints = {ctSortOrder.ReadPortal().Get(saddle), ctSortOrder.ReadPortal().Get(max)};
                regularEndpoints = {saddle, max};
                branchHeight = max - saddle;
                branchType = 1;

            }
            else
            {
                assert(false);
            }
        }
        else
        {
            endpoints = {ctSortOrder.ReadPortal().Get(min), ctSortOrder.ReadPortal().Get(max)};
            regularEndpoints = {min, max};
            branchHeight = max - min;
            branchType = 2;
        }

        Float64 a = inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>().ReadPortal().Get(endpoints[0]);
        Float64 b = inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>().ReadPortal().Get(endpoints[1]);

        Float64 isovalue = 0.0;

        if ("half" == branchIsovalueFlag)
        {
            isovalue = (a + b) / 2.0;
#ifdef DEBUG_PRINT
            std::cout << "Extracting at half the branch isovalue: " << isovalue
                      << ", where its end values are: " << a << " and " << b << std::endl;
#endif
        }
        else
        {
            vtkm::Float64 epsilon = 0.00000001f;
            if (0 == branchType)
            {
                isovalue = b - epsilon;
#ifdef DEBUG_PRINT
                std::cout << "Extracting at branch isovalue - epsilon: " << isovalue
                          << ", where its top end value is: " << b << std::endl;
#endif
            }
            else if (1 == branchType)
            {
                isovalue = a + epsilon;
#ifdef DEBUG_PRINT
                std::cout << "Extracting at branch isovalue + epsilon: " << isovalue
                          << ", where its bottom end value is: " << a << std::endl;
#endif
            }
            else if (2 == branchType)
            {
                isovalue = a + (b - a) / 2.0;
#ifdef DEBUG_PRINT
                std::cout << "Extracting at branch isovalue + (b-a)/2: " << isovalue
                          << ", where its end values are: " << a << " and " << b << std::endl;
#endif
            }
            else
            {   // Print even without debug mode because something weird is going on
                std::cout << "Unexpected branch type encountered" << std::endl;
                assert(false);
            }

        }

#ifdef DEBUG_PRINT
        std::cout << "BID " << i << "\t range = [" << std::setw(8) << std::fixed << std::setprecision(8) << a
                  << " -> " << std::setw(8) << std::fixed <<std::setprecision(8) << b
                  << "], isovalue = " << isovalue << " (volume = " << branchPointVolumeArray.ReadPortal().Get(i)
                  << ")" << std::endl;
#endif

        branchIsovalueArray.WritePortal().Set(i, isovalue);
        branchHeightArray.WritePortal().Set(i, branchHeight);
        branchIsovalueHeightArray.WritePortal().Set(i, abs(a - b));

        branchEndpoints.WritePortal().Set(i, {endpoints[0], endpoints[1]});
        branchEndpointsRegular.WritePortal().Set(i, {regularEndpoints[0], regularEndpoints[1]});
    }
}


void cv1k::interface::computeAdditionalBranchDataFloat(
        const vtkm::cont::ArrayHandle<Float64> fakeFieldArray,
//        const cont::DataSet inputData,    // passing a fakeFieldArray instead
//        const string fieldName,           // passing a fakeFieldArray instead
        const worklet::contourtree_augmented::ContourTree ct,
        const vtkm::cont::ArrayHandle<vtkm::Id> ctSortOrder,
        const cont::ArrayHandle<Id> whichBranch,
        const cont::ArrayHandle<Id> branchMinimum,
        const cont::ArrayHandle<Id> branchSaddle,
        const cont::ArrayHandle<Id> branchMaximum,
        vector<Id> &whichBranchRegular,
        cont::ArrayHandle<Vec<Id, 2>> &branchEndpointsRegular,
        cont::ArrayHandle<Vec<Id, 2>> &branchEndpoints,
        cont::ArrayHandle<Id> &branchHeightArray,
        cont::ArrayHandle<Float64> &branchIsovalueHeightArray,
        cont::ArrayHandle<Float64> &branchIsovalueArray,
        cont::ArrayHandle<Float64> &branchPointVolumeArray,
        const cont::ArrayHandle<Float64> superarcIntrinsicWeight,
        const cont::ArrayHandle<Float64> superarcDependentWeight,
        const std::string branchIsovalueFlag
        )
{
    // Compute the volume of every branch.
    for (int i = 0 ; i < branchPointVolumeArray.GetNumberOfValues() ; i++)
    {
        branchPointVolumeArray.WritePortal().Set(i, 1);
    }

    using vtkm::worklet::contourtree_augmented::NO_SUCH_ELEMENT;
    // @TODO Should I use intrinsic of dependent weight?
    std::cout << "SuperArc to Branch ID mappings (num. mappings: "
              << superarcIntrinsicWeight.GetNumberOfValues() << "=num. of superarcs):" << std::endl;
    for (int i = 0 ; i < superarcIntrinsicWeight.GetNumberOfValues() ; i++)
    {
        Id branchId = whichBranch.ReadPortal().Get(i);

#ifdef DEBUG_PRINT
        std::cout << "SA " << i << " -> BID: " << branchId << std::endl;
#endif

        if (branchId != NO_SUCH_ELEMENT)
        {
            branchPointVolumeArray.WritePortal().Set(branchId, 1);
            Id currentWeight = branchPointVolumeArray.ReadPortal().Get(branchId);
            //Id weight = superarcIntrinsicWeight.ReadPortal().Get(i);
            Id weight = superarcDependentWeight.ReadPortal().Get(i);

            branchPointVolumeArray.WritePortal().Set(branchId, currentWeight + weight);
        }
    }

    std::vector<std::tuple<vtkm::Id, vtkm::Id, vtkm::Id, vtkm::Float64>> branchEndpointsStd;

    // Converts the branches to mesh vertex endpoints and compute isovalue
    std::cout << "Computing isovalues for each branch (total num. isovalues: "
              << branchMinimum.GetNumberOfValues() << "=num. of branches):" << std::endl;

    for (int i = 0 ; i < branchMinimum.GetNumberOfValues() ; i++)
    {
        using vtkm::worklet::contourtree_augmented::MaskedIndex;

        Id branchHeight = 0;
        Vec<Id, 2> endpoints;
        Vec<Id, 2> regularEndpoints;

        Id max = MaskedIndex(ct.Supernodes.ReadPortal().Get(MaskedIndex(branchMaximum.ReadPortal().Get(i))));
        Id min = MaskedIndex(ct.Supernodes.ReadPortal().Get(MaskedIndex(branchMinimum.ReadPortal().Get(i))));
        Id saddle = MaskedIndex(ct.Supernodes.ReadPortal().Get(MaskedIndex(branchSaddle.ReadPortal().Get(i))));

        // 0 for descending, 1 for ascending, 2 for Master,
        int branchType = -1;

        // Determine the type of the branch and set its endpoints
        using vtkm::worklet::contourtree_augmented::NoSuchElement;
        if (false == NoSuchElement(branchSaddle.ReadPortal().Get(static_cast<vtkm::Id>(i))))
        {
            if (min < saddle)
            {
                endpoints = {ctSortOrder.ReadPortal().Get(min), ctSortOrder.ReadPortal().Get(saddle)};
                regularEndpoints = {min, saddle};
                branchHeight = saddle - min;
                branchType = 0;
            }
            else if (saddle < max)
            {
                endpoints = {ctSortOrder.ReadPortal().Get(saddle), ctSortOrder.ReadPortal().Get(max)};
                regularEndpoints = {saddle, max};
                branchHeight = max - saddle;
                branchType = 1;

            }
            else
            {
                assert(false);
            }
        }
        else
        {
            endpoints = {ctSortOrder.ReadPortal().Get(min), ctSortOrder.ReadPortal().Get(max)};
            regularEndpoints = {min, max};
            branchHeight = max - min;
            branchType = 2;
        }


        std::cout << "Branch ENDPOINTS    : " << endpoints[0] << " and " << endpoints[1] << std::endl;

        Float64 a = fakeFieldArray.ReadPortal().Get(endpoints[0]);
                //inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>().ReadPortal().Get(endpoints[0]);
        Float64 b = fakeFieldArray.ReadPortal().Get(endpoints[1]);
                //inputData.GetField(fieldName).GetData().AsArrayHandle<cont::ArrayHandle<Float64>>().ReadPortal().Get(endpoints[1]);

        std::cout << "Branch ENDPOINT VALS: " << a << " and " << b << std::endl;

        Float64 isovalue = 0.0;

        if ("half" == branchIsovalueFlag)
        {
            isovalue = (a + b) / 2.0;
#ifdef DEBUG_PRINT
            std::cout << "Extracting at half the branch isovalue: " << isovalue
                      << ", where its end values are: " << a << " and " << b << std::endl;
#endif
        }
        else
        {
            vtkm::Float64 epsilon = 0.00000001f;
            if (0 == branchType)
            {
                isovalue = b - epsilon;
#ifdef DEBUG_PRINT
                std::cout << "Extracting at branch isovalue - epsilon: " << isovalue
                          << ", where its top end value is: " << b << std::endl;
#endif
            }
            else if (1 == branchType)
            {
                isovalue = a + epsilon;
#ifdef DEBUG_PRINT
                std::cout << "Extracting at branch isovalue + epsilon: " << isovalue
                          << ", where its bottom end value is: " << a << std::endl;
#endif
            }
            else if (2 == branchType)
            {
                isovalue = a + (b - a) / 2.0;
#ifdef DEBUG_PRINT
                std::cout << "Extracting at branch isovalue + (b-a)/2: " << isovalue
                          << ", where its end values are: " << a << " and " << b << std::endl;
#endif
            }
            else
            {   // Print even without debug mode because something weird is going on
                std::cout << "Unexpected branch type encountered" << std::endl;
                assert(false);
            }

        }

#ifdef DEBUG_PRINT
        std::cout << "BID " << i << "\t range = [" << std::setw(8) << std::fixed << std::setprecision(8) << a
                  << " -> " << std::setw(8) << std::fixed <<std::setprecision(8) << b
                  << "], isovalue = " << isovalue << " (volume = " << branchPointVolumeArray.ReadPortal().Get(i)
                  << ")" << std::endl;
#endif

        branchIsovalueArray.WritePortal().Set(i, isovalue);
        branchHeightArray.WritePortal().Set(i, branchHeight);
        branchIsovalueHeightArray.WritePortal().Set(i, abs(a - b));

        branchEndpoints.WritePortal().Set(i, {endpoints[0], endpoints[1]});
        branchEndpointsRegular.WritePortal().Set(i, {regularEndpoints[0], regularEndpoints[1]});

        std::cout << "Branch ISOVALUE: " << isovalue << std::endl;
    }
}
