// -*- c++ -*-
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCGNSReaderInternal.h

  Copyright (c) 2013-2014 Mickael Philit
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
#include "vtkCGNSReaderInternal.h"

#include <cgns_io.h> // Low level IO for fast parsing
#include <iostream>
#include <algorithm>

namespace CGNSRead
{

//------------------------------------------------------------------------------
bool testValidVector(const CGNSVector& item)
{
  // apply some logic and return true or false
  return (item.numComp == 0);
}

//------------------------------------------------------------------------------
void fillVectorsFromVars(std::vector< CGNSRead::CGNSVariable >&  vars,
                            std::vector< CGNSRead::CGNSVector >& vectors,
                            const int physicalDim)
{
  // get number of scalars and vectors
  const size_t nvar = vars.size();
  size_t len;
  char_33 name;

  for (size_t n = 0; n < nvar; ++n)
    {
    vars[n].isComponent = false;
    vars[n].xyzIndex = 0;
    }

  for (size_t n = 0; n < nvar; ++n)
    {
    len = strlen(vars[n].name) - 1;
    switch (vars[n].name[len])
      {
      case 'X':
        vars[n].xyzIndex = 1;
        vars[n].isComponent = true;
        break;
      case 'Y':
        vars[n].xyzIndex = 2;
        vars[n].isComponent = true;
        break;
      case 'Z':
        vars[n].xyzIndex = 3;
        vars[n].isComponent = true;
        break;
      }
    if (vars[n].isComponent == true)
      {
      strcpy(name, vars[n].name);
      name[len] = '\0';
      std::vector< CGNSRead::CGNSVector >::iterator iter =
        CGNSRead::getVectorFromName(vectors, name);
      if (iter != vectors.end())
        {
        iter->numComp += vars[n].xyzIndex ;
        iter->xyzIndex[vars[n].xyzIndex-1] = (int) n ;
        }
      else
        {
        CGNSRead::CGNSVector newVector;
        newVector.xyzIndex[0] = -1;
        newVector.xyzIndex[1] = -1;
        newVector.xyzIndex[2] = -1;
        newVector.numComp = vars[n].xyzIndex ;
        newVector.xyzIndex[vars[n].xyzIndex-1] = (int) n ;
        strcpy(newVector.name , name);
        vectors.push_back(newVector);
        }
      }
    }

  // Detect and tag invalid vector :
  bool invalid = false;
  for (std::vector<CGNSRead::CGNSVector>::iterator iter = vectors.begin();
       iter != vectors.end(); ++iter)
    {
    // Check if number of component agrees with phys_dim
    if (((physicalDim == 3) && (iter->numComp != 6)) ||
        ((physicalDim == 2) && (iter->numComp != 3)))
      {
      for (int index = 0; index < physicalDim; index++)
        {
        int nv = iter->xyzIndex[index];
        if (nv >= 0)
          {
          vars[nv].isComponent = false;
          }
        }
      iter->numComp = 0;
      invalid = true;
      }
    // Check if a variable is present with a similar
    // name as the vector being built
    if (CGNSRead::isACGNSVariable(vars, iter->name) == true)
      {
      //vtkWarningMacro ( "Warning, vector " << iter->name
      //                  << " can't be assembled." << std::endl );
      for (int index = 0; index < physicalDim; index++)
        {
        int n = iter->xyzIndex[index];
        if (n >= 0)
          {
          vars[n].isComponent = false;
          }
        }
      iter->numComp = 0;
      invalid = true;
      }
    if (iter->numComp > 0)
      {
      // Check if DataType_t are identical for all components
      if ((vars[iter->xyzIndex[0]].dt != vars[iter->xyzIndex[1]].dt) ||
          (vars[iter->xyzIndex[0]].dt != vars[iter->xyzIndex[physicalDim-1]].dt))
        {
        for (int index = 0; index < physicalDim; index++)
          {
          int n = iter->xyzIndex[index];
          if (n >= 0)
            {
            vars[n].isComponent = false;
            }
          }
        iter->numComp = 0;
        invalid = true;
        }
      }
    }
  // Remove invalid vectors
  if (invalid == true)
    {
    vectors.erase(std::remove_if(vectors.begin(), vectors.end(),
                                 CGNSRead::testValidVector),
                  vectors.end());
    }
}

//------------------------------------------------------------------------------
template <typename T>
int readNodeData(int cgioNum, double nodeId, std::vector<T>& data)
{
  int n;
  cgsize_t size = 1;
  cgsize_t dimVals[12];
  int ndim;
  T* tmpData;

  if (cgio_get_dimensions(cgioNum, nodeId, &ndim, dimVals) != CG_OK)
    {
    cgio_error_exit("cgio_get_dimensions");
    return 1;
    }

  // allocate data
  for (n = 0; n < ndim; n++)
    {
    size *= dimVals[n];
    }
  if (size <= 0)
    {
    return 1;
    }

  data.resize(size);
  tmpData = new T[size];

  // read data
  if (cgio_read_all_data(cgioNum, nodeId, tmpData) != CG_OK)
    {
    delete [] tmpData;
    return 1;
    }

  for (size_t ii = 0; ii < data.size(); ++ii)
    {
    data[ii] = tmpData[ii];
    }

  delete [] tmpData;
  return 0;
}

//------------------------------------------------------------------------------
// Specialize char array
template <>
int readNodeData<char>(int cgioNum, double nodeId, std::vector<char>& data)
{
  int n;
  cgsize_t size=1;
  cgsize_t dimVals[12];
  int ndim;
  char* tmpData;

  if (cgio_get_dimensions(cgioNum, nodeId, &ndim, dimVals) != CG_OK)
    {
    cgio_error_exit("cgio_get_dimensions");
    return 1;
    }

  // allocate data
  for (n = 0; n < ndim; n++)
    {
    size *= dimVals[n];
    }
  if (size <= 0)
    {
    return 1;
    }

  data.resize(size+1);
  tmpData = new char[size];

  // read data
  if (cgio_read_all_data(cgioNum, nodeId, tmpData) != CG_OK)
    {
    delete [] tmpData;
    return 1;
    }

  for (size_t ii=0; ii < data.size(); ++ii)
    {
    data[ii] = tmpData[ii];
    }
  data[size] = '\0';

  delete [] tmpData;
  return 0;
}

//------------------------------------------------------------------------------
int getNodeChildrenId(int cgioNum, double fatherId,
                          std::vector<double>& childrenIds)
{
  int nchildren;
  int len;

  cgio_number_children(cgioNum, fatherId, &nchildren);

  childrenIds.resize(nchildren);
  double *idList = new double[nchildren];

  cgio_children_ids(cgioNum, fatherId, 1, nchildren, &len, idList);

  if (len != nchildren)
    {
    delete[] idList;
    std::cerr << "Mismatch in number of children and child IDs read"
              << std::endl;
    return 1;
    }

  for (int child = 0; child < nchildren; child++)
    {
    childrenIds[child] = idList[child];
    }

  delete [] idList;
  return 0;
}

//------------------------------------------------------------------------------
int readBaseIds(int cgioNum, double rootId,
                    std::vector<double>& baseIds)
{
  CGNSRead::char_33 nodeLabel;
  size_t nbases = 0;
  size_t nc;

  baseIds.clear();
  getNodeChildrenId(cgioNum, rootId, baseIds);
  if (baseIds.size()  < 1)
    {
    std::cerr << "Error: Not enough nodes under the root description file."
              << std::endl;
    return 1;
    }

  for (nbases = 0, nc = 0; nc < baseIds.size(); nc++)
    {
    if (cgio_get_label(cgioNum, baseIds[nc], nodeLabel) != CG_OK)
      {
      return 1;
      }
    if (strcmp(nodeLabel, "CGNSBase_t") == 0)
      {
      if (nbases < nc)
        {
        baseIds[nbases] = baseIds[nc];
        }
      nbases++;
      }
    else
      {
      cgio_release_id(cgioNum, baseIds[nc]);
      }
    }
  baseIds.resize(nbases);

  if (nbases < 1)
    {
    std::cerr << "Error: Not enough bases in the file." << std::endl;
    return 1;
    }

  return 0;
}

//------------------------------------------------------------------------------
int readBaseCoreInfo(int cgioNum, double baseId,
                          CGNSRead::BaseInformation &baseInfo)
{
  CGNSRead::char_33 dataType;
  std::vector<int> mdata;

  if (cgio_get_name(cgioNum, baseId, baseInfo.name) != CG_OK)
    {
    std::cerr << "cgio_get_name" << std::endl;
    return 1;
    }

  // read node data type
  if (cgio_get_data_type(cgioNum , baseId, dataType) != CG_OK)
    {
    return 1;
    }

  if (strcmp(dataType, "I4") != 0)
    {
    std::cerr << "Unexpected data type for dimension data of base"
              << std::endl;
    return 1;
    }

  CGNSRead::readNodeData<int>(cgioNum, baseId, mdata);

  baseInfo.cellDim = mdata[0];
  baseInfo.physicalDim = mdata[1];

  return 0;
}

//------------------------------------------------------------------------------
int readBaseIteration(int cgioNum, double nodeId,
                          CGNSRead::BaseInformation& baseInfo)
{
  CGNSRead::char_33 nodeLabel;
  CGNSRead::char_33 nodeName;
  CGNSRead::char_33 dataType;

  bool createTimeStates = true;
  bool createIterStates = true;

  std::vector<int> ndata;
  // read node data type
  if (cgio_get_data_type(cgioNum , nodeId , dataType) != CG_OK)
    {
    return 1;
    }

  if (strcmp(dataType, "I4") != 0)
    {
    std::cerr << "Unexpected data type for iteration number of steps"
              << std::endl;
    return 1;
    }

  CGNSRead::readNodeData<int>(cgioNum, nodeId, ndata);

  int nstates = ndata[0];
  std::vector<double> childrenIterative;

  getNodeChildrenId(cgioNum, nodeId, childrenIterative);

  for(size_t nc = 0; nc < childrenIterative.size(); ++nc)
    {
    if (cgio_get_label(cgioNum, childrenIterative[nc], nodeLabel) != CG_OK)
      {
      return 1;
      }

    if (cgio_get_name(cgioNum, childrenIterative[nc], nodeName) != CG_OK)
      {
      return 1;
      }

    bool isDataArray = (strcmp(nodeLabel, "DataArray_t") == 0);
    if (isDataArray && (strcmp(nodeName, "TimeValues") == 0))
      {
      // Read time values
      // read node data type
      if (cgio_get_data_type(cgioNum,
                             childrenIterative[nc],
                             dataType) != CG_OK)
        {
        return 1;
        }
      baseInfo.times.clear();
      if (strcmp(dataType, "R8") == 0)
        {
        CGNSRead::readNodeData<double>(cgioNum,
                                       childrenIterative[nc],
                                       baseInfo.times);
        }
      else if (strcmp(dataType, "R4") == 0)
        {
        std::vector<float> iteData;
        CGNSRead::readNodeData<float>(cgioNum,
                                      childrenIterative[nc],
                                      iteData);
        baseInfo.times.resize(iteData.size());
        for (size_t ii = 0; ii < iteData.size(); ii++)
          {
          baseInfo.times[ii] = (double) iteData[ii];
          }
        }
      else
        {
        std::cerr << "Unexpected data type for iterative data"
                  << std::endl;
        return 1;
        }

      if (static_cast<int>(baseInfo.times.size()) != nstates)
        {
        std::cerr << "Error reading times node";
        return 1;
        }

      createTimeStates = false;
      }
    else if (isDataArray && (strcmp(nodeName, "IterationValues") == 0))
      {
      // Read time steps
      // read node data type
      if (cgio_get_data_type(cgioNum, childrenIterative[nc],
                              dataType) != CG_OK)
        {
        return 1;
        }

      if (strcmp(dataType, "I4") != 0)
        {
        std::cerr << "Unexpected data type for iterative data"
                  << std::endl;
        return 1;
        }

      baseInfo.steps.clear();
      CGNSRead::readNodeData<int> ( cgioNum, childrenIterative[nc],
                                      baseInfo.steps );
      if (static_cast<int>(baseInfo.steps.size()) != nstates)
        {
        std::cerr << "Error reading steps node";
        return 1;
        }
      createIterStates = false;
      }
    else
      {
      cgio_release_id(cgioNum, childrenIterative[nc]);
      }
    }

  if (createIterStates == true)
    {
    for (int i = 0; i < nstates; ++i)
      {
      baseInfo.steps.push_back(i);
      }
    }
  if (createTimeStates == true)
    {
    for (int i = 0; i < nstates; ++i)
      {
      baseInfo.times.push_back((double) baseInfo.steps[i]);
      }
    }
  return 0;
}

//------------------------------------------------------------------------------
int readZoneIterInfo(int cgioNum, double nodeId,
                        CGNSRead::BaseInformation& baseInfo)
{
  CGNSRead::char_33 nodeLabel;
  CGNSRead::char_33 nodeName;
  std::vector<double> iterChildId;

  getNodeChildrenId(cgioNum, nodeId, iterChildId);

  for (size_t nn = 0; nn < iterChildId.size(); nn++)
    {

    if (cgio_get_name(cgioNum, iterChildId[nn], nodeName) != CG_OK)
      {
      return 1;
      }
    if (cgio_get_label(cgioNum, iterChildId[nn], nodeLabel) != CG_OK)
      {
      return 1;
      }
    bool isDataArray = (strcmp(nodeLabel, "DataArray_t") == 0);
    if ( isDataArray &&
         (strcmp(nodeName, "GridCoordinatesPointers") == 0))
      {
      baseInfo.useGridPointers = true ;
      }
    else if ( isDataArray &&
              ( strcmp ( nodeName, "FlowSolutionPointers" ) == 0 ) )
      {
      baseInfo.useFlowPointers = true ;
      // Maybe load FlowSolutionPointers once and for all
      }
    cgio_release_id ( cgioNum, iterChildId[nn] );
    }
  return 0;
}

//------------------------------------------------------------------------------
int readSolInfo(int cgioNum, double nodeId,
                  CGNSRead::BaseInformation& baseInfo)
{
  CGNS_ENUMT(GridLocation_t) varCentering ;
  varCentering = CGNS_ENUMV(Vertex);

  CGNSRead::char_33 nodeLabel;

  std::vector<double> solChildId;

  getNodeChildrenId(cgioNum, nodeId, solChildId);

  std::vector< CGNSRead::CGNSVariable > cgnsVars;
  std::vector< CGNSRead::CGNSVector > cgnsVectors;

  size_t nn;
  size_t nvars = 0;

  for (nvars = 0, nn = 0; nn < solChildId.size(); nn++)
    {
    if (cgio_get_label(cgioNum, solChildId[nn], nodeLabel) != CG_OK)
      {
      std::cerr << "Error while reading nodelabel" << std::endl ;
      return 1;
      }

    if ( strcmp(nodeLabel, "DataArray_t") == 0)
      {
      CGNSRead::CGNSVariable curVar;

      if ( cgio_get_name ( cgioNum, solChildId[nn], curVar.name ) != CG_OK )
        {
        return 1;
        }
      curVar.isComponent = false;
      curVar.xyzIndex = 0;

      // read node data type
      CGNSRead::char_33 dataType;
      if (cgio_get_data_type(cgioNum , solChildId[nn], dataType))
        {
        continue;
        }
      if (strcmp(dataType, "R8") == 0)
        {
        curVar.dt = CGNS_ENUMV(RealDouble);
        }
      else if (strcmp(dataType, "R4") == 0)
        {
        curVar.dt = CGNS_ENUMV(RealSingle);
        }
      else if (strcmp(dataType, "I4") == 0)
        {
        curVar.dt = CGNS_ENUMV(Integer);
        }
      else if (strcmp(dataType, "I8") == 0)
        {
        curVar.dt = CGNS_ENUMV(LongInteger);
        }
      else
        {
        continue;
        }

      cgnsVars.push_back(curVar);
      cgio_release_id(cgioNum, solChildId[nn]);
      nvars++;
      }
    else if (strcmp(nodeLabel, "GridLocation_t") == 0)
      {
      CGNSRead::char_33 dataType;

      if (cgio_get_data_type(cgioNum, solChildId[nn], dataType) != CG_OK)
        {
        return 1;
        }

      if (strcmp(dataType, "C1") != 0)
        {
        std::cerr << "Unexpected data type for GridLocation_t node"
                  << std::endl;
        return 1;
        }

      std::vector<char> location;
      CGNSRead::readNodeData<char>(cgioNum, solChildId[nn], location);

      if (strcmp(location.data(), "Vertex") == 0)
        {
        varCentering = CGNS_ENUMV(Vertex);
        }
      else if (strcmp(location.data(), "CellCenter") == 0)
        {
        varCentering = CGNS_ENUMV(CellCenter);
        }
      else
        {
        varCentering = CGNS_ENUMV(GridLocationNull);
        }
      }
    else
      {
      cgio_release_id(cgioNum, solChildId[nn]);
      }
    }
  //
  if (varCentering != CGNS_ENUMV(Vertex) &&
      varCentering != CGNS_ENUMV(CellCenter))
    {
    std::cerr << "Bad Centering encountered !"
              "Only Vertex and CellCenter are supported"
              << std::endl;
    return 1;
    }

  CGNSRead::fillVectorsFromVars(cgnsVars, cgnsVectors, baseInfo.physicalDim);

  for (size_t ii=0; ii < cgnsVars.size(); ++ii)
    {
    if (cgnsVars[ii].isComponent == true)
      {
      continue ;
      }
    switch (varCentering)
      {
      case CGNS_ENUMV(Vertex):
        if (!baseInfo.PointDataArraySelection.HasArray(cgnsVars[ii].name))
          {
          baseInfo.PointDataArraySelection.AddArray(cgnsVars[ii].name, false);
          }

        break;
      case CGNS_ENUMV(CellCenter):
        if (!baseInfo.CellDataArraySelection.HasArray(cgnsVars[ii].name))
          {
          baseInfo.CellDataArraySelection.AddArray(cgnsVars[ii].name, false);
          }
        break;
      }
    }
  for ( size_t jj = 0; jj <cgnsVectors.size(); ++jj )
    {
    switch (varCentering)
      {
      case CGNS_ENUMV(Vertex):
        if (!baseInfo.PointDataArraySelection.HasArray(cgnsVectors[jj].name))
          {
          baseInfo.PointDataArraySelection.AddArray(cgnsVectors[jj].name, false);
          }

        break;
      case CGNS_ENUMV(CellCenter):
        if (!baseInfo.CellDataArraySelection.HasArray(cgnsVectors[jj].name))
          {
          baseInfo.CellDataArraySelection.AddArray(cgnsVectors[jj].name, false);
          }
        break;
      }
    }
  return 0;
}

//------------------------------------------------------------------------------
int readBaseFamily(int cgioNum, double nodeId,
                     CGNSRead::BaseInformation& baseInfo)
{
  CGNSRead::FamilyInformation curFamily;
  CGNSRead::char_33 nodeLabel;
  std::vector<double> famChildId;

  if (cgio_get_name(cgioNum, nodeId, curFamily.name) != CG_OK)
    {
    return 1;
    }
  curFamily.isBC = false ;

  getNodeChildrenId(cgioNum, nodeId, famChildId);

  for (size_t nn = 0; nn < famChildId.size(); nn++)
    {
    if (cgio_get_label(cgioNum, famChildId[nn], nodeLabel) != CG_OK)
      {
      return 1;
      }
    if (strcmp(nodeLabel, "FamilyBC_t") == 0)
      {
      curFamily.isBC = true ;
      }
    cgio_release_id(cgioNum, famChildId[nn]);
    }
  baseInfo.family.push_back(curFamily);

  return 0;
}

//------------------------------------------------------------------------------
int readBaseReferenceState(int cgioNum, double nodeId,
                           CGNSRead::BaseInformation& baseInfo)
{
  CGNSRead::char_33 nodeLabel;
  CGNSRead::char_33 curName;

  std::vector<double> children;
  getNodeChildrenId(cgioNum, nodeId, children);

  size_t nn;
  for (nn = 0; nn < children.size(); nn++)
    {
    if (cgio_get_label(cgioNum, children[nn], nodeLabel) != CG_OK)
      {
      std::cerr << "Error while reading nodelabel" << std::endl;
      return 1;
      }
    if (cgio_get_name(cgioNum, children[nn], curName) != CG_OK)
      {
      return 1;
      }
    bool isDataArray = strcmp(nodeLabel, "DataArray_t") == 0;
    if (isDataArray && ((strcmp(curName, "Mach") == 0 ) ||
                        (strcmp(curName, "SpecificHeatRatio") == 0) ||
                        (strcmp(curName, "IdealGasConstant") == 0) ||
                        (strcmp(curName, "SpecificHeatVolume") == 0) ||
                        (strcmp(curName, "SpecificHeatPressure") == 0)))
      {
      // read node data type
      CGNSRead::char_33 dataType;
      if (cgio_get_data_type(cgioNum , children[nn], dataType))
        {
        return 1;
        }

      if (strcmp(dataType, "R8") == 0)
        {
        std::vector<double> bdata;
        CGNSRead::readNodeData<double>(cgioNum, children[nn], bdata);
        baseInfo.referenceState[curName] = (double) bdata[0];
        }
      else if (strcmp(dataType, "R4") == 0)
        {
        std::vector<float> bdata;
        CGNSRead::readNodeData<float> ( cgioNum, children[nn], bdata );
        baseInfo.referenceState[curName] = (double) bdata[0];
        }
      else
        {
        std::cerr << "Unexpected data in ReferenceState_t"
                  << std::endl;
        return 1;
        }
      }
    cgio_release_id(cgioNum, children[nn]);
    }
  return 0;
}

//------------------------------------------------------------------------------
int readZoneInfo(int cgioNum, double nodeId,
                 CGNSRead::BaseInformation& baseInfo)
{
  CGNSRead::char_33 nodeLabel;
  std::vector<double> zoneChildId;

  getNodeChildrenId(cgioNum, nodeId, zoneChildId);

  int nflows = 0;
  size_t nn;
  for (nflows = 0, nn = 0; nn < zoneChildId.size(); nn++)
    {

    if (cgio_get_label(cgioNum, zoneChildId[nn], nodeLabel) != CG_OK)
      {
      std::cerr << "Error while reading nodelabel" << std::endl;
      return 1;
      }

    if ((nflows < 3) && (strcmp(nodeLabel, "FlowSolution_t") == 0))
      {
      // Read only 3 Flow Solution to have a chance
      // to get Cell and Vertex variables
      // C=Cell V=Vertex
      // Layout sample :
      //    1. C init state (this one may be not processed due
      //                     to FlowSolutionPointers but we still
      //                     want some information about the two next node)
      //    2. C time 1s
      //    3. V time 1s

      readSolInfo(cgioNum, zoneChildId[nn], baseInfo);
      nflows++;
      }
    else if (strcmp(nodeLabel, "ZoneIterativeData_t") == 0)
      {
      // get time information
      readZoneIterInfo(cgioNum, zoneChildId[nn], baseInfo);
      }
    cgio_release_id(cgioNum, zoneChildId[nn]);
    }
  return 0;
}

//------------------------------------------------------------------------------
bool vtkCGNSMetaData::Parse(const char* cgnsFileName)
{

  if (!cgnsFileName)
    {
    return false;
    }

  if (this->LastReadFilename == cgnsFileName)
    {
    return true;
    }

  int cgioNum;
  int ier;
  double rootId;
  char nodeLabel[CGIO_MAX_NAME_LENGTH+1];

  // use cgio routine to open the file
  if (cgio_open_file(cgnsFileName, CGIO_MODE_READ, 0, &cgioNum ) != CG_OK)
    {
    cgio_error_exit("cgio_file_open");
    }
  if ( cgio_get_root_id(cgioNum, &rootId) != CG_OK)
    {
    cgio_error_exit("cgio_get_root_id");
    }

  // Get base id list :
  std::vector<double> baseIds;
  ier = readBaseIds(cgioNum, rootId, baseIds);
  if ( ier != 0 )
    {
    return false;
    }

  if ( this->baseList.size() > 0 )
    {
    this->baseList.clear();
    }
  this->baseList.resize(baseIds.size());
  // Read base list
  for (size_t numBase=0; numBase < baseIds.size(); numBase++)
    {
    // base names for later selection
    readBaseCoreInfo(cgioNum, baseIds[numBase], this->baseList[numBase]);

    std::vector<double> baseChildId;

    getNodeChildrenId(cgioNum, baseIds[numBase], baseChildId);

    size_t nzones = 0;
    size_t nn;
    for (nzones = 0, nn = 0; nn < baseChildId.size(); ++nn)
      {
      if (cgio_get_label(cgioNum, baseChildId[nn], nodeLabel) != CG_OK)
        {
        return false;
        }

      if (strcmp(nodeLabel, "Zone_t") == 0)
        {
        if (nzones < nn)
          {
          baseChildId[nzones] = baseChildId[nn];
          }
        nzones++;
        }
      else if (strcmp(nodeLabel, "Family_t") == 0)
        {
        readBaseFamily(cgioNum, baseChildId[nn],
                           this->baseList[numBase]);
        }
      else if (strcmp(nodeLabel, "BaseIterativeData_t") == 0)
        {
        readBaseIteration(cgioNum, baseChildId[nn],
                              this->baseList[numBase]);
        }
      else if (strcmp(nodeLabel, "ReferenceState_t") == 0)
        {
        readBaseReferenceState(cgioNum, baseChildId[nn],
                                    this->baseList[numBase]);
        }
      else
        {
        cgio_release_id(cgioNum, baseChildId[nn]);
        }
      }
    this->baseList[numBase].nzones = static_cast<int>(nzones);

    if (this->baseList[numBase].times.size() < 1)
      {
      // If no time information were found
      // just put default values
      this->baseList[numBase].steps.clear();
      this->baseList[numBase].times.clear();
      this->baseList[numBase].steps.push_back(0);
      this->baseList[numBase].times.push_back(0.0);
      }

    if ( nzones > 0 )
      {
      // variable name and more, based on first zone only
      readZoneInfo(cgioNum, baseChildId[0], this->baseList[numBase]);
      }

    }

  // Same Timesteps in all root nodes
  // or separated time range by root nodes
  // timesteps need to be sorted for each root node
  this->GlobalTime.clear();
  for (size_t numBase=0; numBase < baseList.size(); numBase++)
    {
    if (numBase == 0)
      {
      this->GlobalTime = this->baseList[numBase].times;
      continue;
      }
    const std::vector<double>& times = this->baseList[numBase].times;
    if (times.front() > this->GlobalTime.back())
      {
      this->GlobalTime.insert( this->GlobalTime.end(),
                               times.begin(), times.end());
      }

    if (times.back() < this->GlobalTime.front())
      {
      this->GlobalTime.insert( this->GlobalTime.begin(),
                               times.begin(), times.end());
      }
    }

  this->LastReadFilename = cgnsFileName;
  cgio_close_file ( cgioNum );
  return true;
}

//------------------------------------------------------------------------------
vtkCGNSMetaData::vtkCGNSMetaData()
{
}
//------------------------------------------------------------------------------
vtkCGNSMetaData::~vtkCGNSMetaData()
{
}
//------------------------------------------------------------------------------
void vtkCGNSMetaData::PrintSelf(std::ostream& os)
{
  os << "--> vtkCGNSMetaData"  << std::endl;
  os << "LastReadFileName: " << this->LastReadFilename << std::endl;
  os << "Base information:"  << std::endl;
  for (size_t b=0; b < this->baseList.size(); b++)
    {
    os << "  Base name: "  << this->baseList[b].name << std::endl ;
    os << "    number of zones: " << this->baseList[b].nzones << std::endl;
    os << "    number of time steps: "<< this->baseList[b].times.size()
       << std::endl;
    os << "    use unsteady grid: "<< this->baseList[b].useGridPointers
       << std::endl;
    os << "    use unsteady flow: "<< this->baseList[b].useFlowPointers
       << std::endl;

    for (int i = 0;
         i < this->baseList[b].PointDataArraySelection.GetNumberOfArrays();
         ++i)
      {
        os << "      Vertex :: ";
        os << this->baseList[b].PointDataArraySelection.GetArrayName(i)
           << std::endl;
      }
    for (int i = 0;
         i < this->baseList[b].CellDataArraySelection.GetNumberOfArrays();
         ++i)
      {
        os << "      Cell :: ";
        os << this->baseList[b].CellDataArraySelection.GetArrayName(i)
           << std::endl;
      }

    os << "    Family Number: "<< this->baseList[b].family.size() << std::endl;
    for (size_t fam=0; fam< this->baseList[b].family.size(); fam++)
      {
      os << "      Family: " << this->baseList[b].family[fam].name << " is BC: "
         << this->baseList[b].family[fam].isBC << std::endl;
      }

    os << "    Reference State:"<< std::endl;
    std::map<std::string, double>::iterator iter;
    for (iter = this->baseList[b].referenceState.begin();
         iter !=this->baseList[b].referenceState.end(); iter++)
      {
      os << "  Variable: " << iter->first;
      os << "  Value: " << iter->second << std::endl;
      }
    }
}

#ifdef PARAVIEW_USE_MPI
//------------------------------------------------------------------------------
static void BroadcastCGNSString(vtkMultiProcessController* ctrl,
                                CGNSRead::char_33 & str)
{
  int len = 33;
  if ( str )
    {
    ctrl->Broadcast(&len, 1, 0);
    ctrl->Broadcast(&str[0], len, 0);
    }
  else
    {
    len = 0;
    ctrl->Broadcast(&len, 1, 0);
    }
}

//------------------------------------------------------------------------------
static void BroadcastString(vtkMultiProcessController* controller,
                            std::string& str, int rank)
{
  unsigned long len = static_cast<unsigned long>(str.size()) + 1;
  controller->Broadcast(&len, 1, 0);
  if (len)
    {
    if (rank)
      {
      std::vector<char> tmp;
      tmp.resize(len);
      controller->Broadcast(&(tmp[0]), len, 0);
      str = &tmp[0];
      }
    else
      {
      const char* start = str.c_str();
      std::vector<char> tmp(start, start + len);
      controller->Broadcast(&tmp[0], len, 0);
      }
    }
}
//------------------------------------------------------------------------------
static void BroadcastDoubleVector(vtkMultiProcessController* controller,
                                  std::vector<double>& dvec, int rank)
{
  unsigned long len = static_cast<unsigned long>(dvec.size());
  controller->Broadcast(&len, 1, 0);
  if (rank)
    {
    dvec.resize(len);
    }
  if (len)
    {
    controller->Broadcast(&dvec[0], len, 0);
    }
}
//------------------------------------------------------------------------------
static void BroadcastIntVector(vtkMultiProcessController* controller,
                               std::vector<int>& ivec, int rank)
{
  unsigned long len = static_cast<unsigned long>(ivec.size());
  controller->Broadcast(&len, 1, 0);
  if (rank)
    {
    ivec.resize(len);
    }
  if (len)
    {
    controller->Broadcast(&ivec[0], len, 0);
    }
}
//------------------------------------------------------------------------------
static void BroadcastSelection(vtkMultiProcessController* controller,
                               CGNSRead::vtkCGNSArraySelection & selection,
                               int rank)
{
  unsigned long len = static_cast<unsigned long>(selection.size());
  controller->Broadcast(&len, 1, 0);
  if ( rank == 0 )
    {
    std::map<std::string, bool>::iterator ite;
    int tmp;
    for (ite = selection.begin(); ite != selection.end(); ++ ite)
      {
      unsigned long len = static_cast<unsigned long>(ite->first.size()) + 1;
      controller->Broadcast(&len, 1, 0);
      if (len)
        {
          const char* start = ite->first.c_str();
          std::vector<char> tmp(start, start + len);
          controller->Broadcast(&tmp[0], len, 0);
        }
      tmp = (int) ite->second;
      controller->Broadcast(&tmp, 1, 0);
      }
    }
  else
    {
    unsigned long i;
    for (i = 0; i < len; ++ i)
      {
      std::string key;
      int tmp;
      CGNSRead::BroadcastString(controller, key, rank);
      selection[key] = false;
      controller->Broadcast(&tmp, 1, 0);
      selection[key] = (bool) tmp;
      }
    }
}

//------------------------------------------------------------------------------
static void BroadcastRefState(vtkMultiProcessController* controller,
                              std::map<std::string,double> & refInfo, int rank)
{
  unsigned long len = static_cast<unsigned long>(refInfo.size());
  controller->Broadcast(&len, 1, 0);
  if (rank == 0)
    {
    std::map<std::string, double>::iterator ite;
    for (ite = refInfo.begin(); ite != refInfo.end(); ++ ite)
      {
      unsigned long len = static_cast<unsigned long>(ite->first.size()) + 1;
      controller->Broadcast(&len, 1, 0);
      if (len)
        {
          const char* start = ite->first.c_str();
          std::vector<char> tmp(start, start + len);
          controller->Broadcast(&tmp[0], len, 0);
        }
      controller->Broadcast(&ite->second, 1, 0);
      }
    }
  else
    {
    for (unsigned long i = 0; i < len; ++i)
      {
      std::string key;
      CGNSRead::BroadcastString(controller, key, rank);
      refInfo[key] = 0.0;
      controller->Broadcast(&refInfo[key], 1, 0);
      }
    }
}

//------------------------------------------------------------------------------
static void BroadcastFamilies(vtkMultiProcessController* controller,
                              std::vector<CGNSRead::FamilyInformation>& famInfo,
                              int rank)
{
  unsigned long len = static_cast<unsigned long>(famInfo.size());
  controller->Broadcast(&len, 1, 0);
  if (rank != 0)
    {
    famInfo.resize(len);
    }
  std::vector<CGNSRead::FamilyInformation>::iterator ite;
  for (ite = famInfo.begin(); ite != famInfo.end(); ++ite)
    {
    BroadcastCGNSString(controller, ite->name);
    int flags = 0;
    if (rank == 0)
      {
      if (ite->isBC == true)
        {
        flags = 1;
        }
      controller->Broadcast(&flags, 1, 0);
      }
    else
      {
      controller->Broadcast(&flags, 1, 0);
      if ((flags & 1) != 0)
        {
        ite->isBC = true;
        }
      }
    }
}

//------------------------------------------------------------------------------
void vtkCGNSMetaData::Broadcast(vtkMultiProcessController* controller,
                                int rank)
{
  unsigned long len = static_cast<unsigned long>(this->baseList.size());
  controller->Broadcast(&len, 1, 0);
  if (rank != 0)
    {
    this->baseList.resize(len);
    }
  std::vector<CGNSRead::BaseInformation>::iterator ite;
  for (ite = this->baseList.begin(); ite != baseList.end(); ++ite)
    {
    CGNSRead::BroadcastCGNSString(controller, ite->name);
    controller->Broadcast(&ite->cellDim, 1, 0);
    controller->Broadcast(&ite->physicalDim, 1, 0);
    controller->Broadcast(&ite->baseNumber, 1, 0);
    controller->Broadcast(&ite->nzones, 1, 0);

    int flags = 0;
    if (rank == 0)
      {
      if (ite->useGridPointers == true)
        {
        flags = 1;
        }
      if (ite->useFlowPointers == true)
        {
        flags = (flags | 2);
        }
      controller->Broadcast(&flags, 1, 0);
      }
    else
      {
      controller->Broadcast(&flags, 1, 0);
      if ((flags & 1) != 0)
        {
        ite->useGridPointers = true;
        }
      if ((flags & 2) != 0)
        {
        ite->useFlowPointers = true;
        }
      }

    CGNSRead::BroadcastRefState(controller, ite->referenceState, rank);
    CGNSRead::BroadcastFamilies(controller, ite->family, rank);

    CGNSRead::BroadcastSelection(controller, ite->PointDataArraySelection,
                                 rank);
    CGNSRead::BroadcastSelection(controller, ite->CellDataArraySelection,
                                 rank);

    BroadcastIntVector(controller, ite->steps, rank);
    BroadcastDoubleVector(controller, ite->times, rank);
    }
  CGNSRead::BroadcastString(controller, this->LastReadFilename, rank);
  BroadcastDoubleVector(controller, this->GlobalTime, rank);
}
#endif

}
