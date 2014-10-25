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

#include "cgio_helpers.h"

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
