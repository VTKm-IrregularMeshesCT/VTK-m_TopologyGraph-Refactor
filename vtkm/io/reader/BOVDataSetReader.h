//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2015 Sandia Corporation.
//  Copyright 2015 UT-Battelle, LLC.
//  Copyright 2015 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef vtk_m_io_reader_BOVDataSetReader_h
#define vtk_m_io_reader_BOVDataSetReader_h

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/ErrorIO.h>
#include <fstream>

namespace vtkm {
namespace io {
namespace reader {

class BOVDataSetReader
{
public:
    BOVDataSetReader(const char *fileName) : FileName(fileName), Loaded(false), DataSet()
    {
    }

    const vtkm::cont::DataSet &
    ReadDataSet()
    {
        try
        {
            LoadFile();
        }
        catch (std::ifstream::failure &e)
        {
            std::string message("IO Error: ");
            throw vtkm::io::ErrorIO(message + e.what());
        }
        
        return this->DataSet;
    }

private:
    typedef enum {ByteData, ShortData, IntegerData, FloatData, DoubleData} DataFormat;

    void parseAssert(bool v)
    {
        if (!v)
            throw vtkm::io::ErrorIO("ParseError");
    }

    void LoadFile()
    {
        if (this->Loaded)
            return;

        std::ifstream stream(this->FileName);
        if (stream.fail())
            throw vtkm::io::ErrorIO("Failed to open file: " + this->FileName);

        DataFormat dataFormat;
        std::string bovFile, line, token, options, variableName;
        vtkm::Id numComponents = 1;
        vtkm::Id3 dim;
        vtkm::Vec<vtkm::FloatDefault,3> origin(0,0,0);
        vtkm::Vec<vtkm::FloatDefault,3> spacing(1,1,1);
        bool spacingSet = false;
        
        while (stream.good())
        {
            std::getline(stream, line);
            if (line.size() == 0 || line[0] == '#')
                continue;
            //std::cout<<"::"<<line<<"::"<<std::endl;
            std::size_t pos = line.find(":");
            if (pos == std::string::npos)
                throw vtkm::io::ErrorIO("Unsupported option: "+line);
            token = line.substr(0,pos);
            options = line.substr(pos+1, line.size()-1);
            //std::cout<<token<<"::"<<options<<std::endl;
            
            std::stringstream strStream(options);

            //Format supports both space and "_" seperated tokens...
            if (token.find("DATA") != std::string::npos &&
                token.find("FILE") != std::string::npos)
            {
                strStream >> bovFile >> std::ws;
            }
            else if (token.find("DATA") != std::string::npos &&
                token.find("SIZE") != std::string::npos)
            {
                strStream >> dim[0] >> dim[1] >> dim[2] >> std::ws;
            }
            else if (token.find("BRICK") != std::string::npos &&
                     token.find("ORIGIN") != std::string::npos)
            {
                strStream >> origin[0] >> origin[1] >> origin[2] >> std::ws;
            }            
            else if (token.find("BRICK") != std::string::npos &&
                     token.find("SIZE") != std::string::npos)
            {
                strStream >> spacing[0] >> spacing[1] >> spacing[2] >> std::ws;
                spacingSet = true;
            }
            else if (token.find("DATA") != std::string::npos &&
                token.find("FORMAT") != std::string::npos)
            {
                std::string opt;
                strStream >> opt >> std::ws;
                if (opt.find("FLOAT") != std::string::npos ||
                    opt.find("REAL") != std::string::npos)
                    dataFormat = FloatData;
                else if (opt.find("DOUBLE") != std::string::npos)
                    dataFormat = DoubleData;
                else
                    throw vtkm::io::ErrorIO("Unsupported data type: "+token);
            }
            else if (token.find("DATA") != std::string::npos ||
                     token.find("COMPONENTS") != std::string::npos)
            {
                strStream >> numComponents >> std::ws;
                if (numComponents != 1 && numComponents != 3)
                    throw vtkm::io::ErrorIO("Unsupported number of components");
            }
            else if (token.find("VARIABLE") != std::string::npos &&
                     token.find("PALETTE") == std::string::npos)
            {
                strStream >> variableName >> std::ws;
                if (variableName[0] == '"')
                    variableName = variableName.substr(1, variableName.size()-2);
            }
            else
                std::cerr<<"Unsupported BOV option: "<<token<<std::endl;
        }

        if (spacingSet)
        {
            spacing[0] = (spacing[0]) / static_cast<vtkm::FloatDefault>(dim[0]-1);
            spacing[1] = (spacing[1]) / static_cast<vtkm::FloatDefault>(dim[1]-1);
            spacing[2] = (spacing[2]) / static_cast<vtkm::FloatDefault>(dim[2]-1);
        }

        //Get whole path for data file.
        std::string fullPathDataFile;
        if (bovFile[0] == '/')
            fullPathDataFile = bovFile;
        else
        {
            //Get base dir.
            std::string baseDir, baseFile;
            std::size_t pos = FileName.rfind("/");
            if (pos != std::string::npos)
                baseDir = this->FileName.substr(0, pos);
            
            if (bovFile.substr(0,2) == "./")
                baseFile = bovFile.substr(2, bovFile.size()-2);
            fullPathDataFile = baseDir + "/" + baseFile;
        }

        vtkm::cont::DataSetBuilderUniform dataSetBuilder;
        vtkm::cont::DataSetFieldAdd dsf;        
        this->DataSet = dataSetBuilder.Create(dim, origin, spacing);

        vtkm::Id numTuples = dim[0]*dim[1]*dim[2];
        if (numComponents == 1)
        {
            if (dataFormat == FloatData)
            {
                vtkm::cont::ArrayHandle<vtkm::Float32> var;
                ReadScalar(fullPathDataFile, numTuples, var);
                dsf.AddPointField(this->DataSet, variableName, var);
            }
            else if (dataFormat == DoubleData)
            {
                vtkm::cont::ArrayHandle<vtkm::Float64> var;
                ReadScalar(fullPathDataFile, numTuples, var);
                dsf.AddPointField(this->DataSet, variableName, var);                
            }            
        }
        else if (numComponents == 3)
        {
            if (dataFormat == FloatData)
            {
                vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3> > var;
                ReadVector(fullPathDataFile, numTuples, var);
                dsf.AddPointField(this->DataSet, variableName, var);
            }
            else if (dataFormat == DoubleData)
            {
                vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3> > var;                
                ReadVector(fullPathDataFile, numTuples, var);
                dsf.AddPointField(this->DataSet, variableName, var);
            }                        
        }
    }

    template<typename T>
    void
    ReadBuffer(const std::string &fName, const vtkm::Id &sz,
               std::vector<T> &buff)
    {
        FILE *fp = fopen(fName.c_str(), "rb");
        if (fp == NULL)
            throw vtkm::io::ErrorIO("Unable to open data file: "+fName);
        buff.resize(sz);
        std::size_t nread = fread(&buff[0], sizeof(T), sz, fp);
        if (nread != sz)
            throw vtkm::io::ErrorIO("Data file read failed: "+fName);
        fclose(fp);
    }
    
    template<typename T>
    void
    ReadScalar(const std::string &fName, const vtkm::Id &nTuples,
               vtkm::cont::ArrayHandle<T> &var)
    {
        std::vector<T> buff;
        ReadBuffer(fName, nTuples, buff);

        var.Allocate(nTuples);
        for (vtkm::Id i = 0; i < nTuples; i++)
            var.GetPortalControl().Set(i, buff[i]);
    }

    template<typename T>
    void
    ReadVector(const std::string &fName, const vtkm::Id &nTuples,
               vtkm::cont::ArrayHandle<vtkm::Vec<T,3> > &var)
    {
        std::vector<T> buff;
        ReadBuffer(fName, nTuples*3, buff);

        var.Allocate(nTuples);
        vtkm::Vec<T,3> v;
        for (vtkm::Id i = 0; i < nTuples; i++)
        {
            v[0] = buff[i*3+0];
            v[1] = buff[i*3+1];
            v[2] = buff[i*3+2];            
            var.GetPortalControl().Set(i, v);
        }
    }    

    bool Loaded;
    std::string FileName;
    vtkm::cont::DataSet DataSet;
};

}
}
} // vtkm::io::reader

#endif // vtk_m_io_reader_BOVReader_h    
