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
// .NAME vtkCGNSReaderInternal -- class(es) used by vtkCGNSReader
// .SECTION Description
//     parse a file in "CGNS" format
//
// .SECTION Caveats
//     Only Cell/Vertex data are supported.
//
// .SECTION Thanks
// Thanks to .

#ifndef __CGNS_READER_INTERNAL__
#define __CGNS_READER_INTERNAL__

#include "vtkPVConfig.h"     // For PARAVIEW_USE_MPI

#include <cgnslib.h> // DataType, and other definition

#include <vector>
#include <map>
#include <string>
#include <string.h> // for inline strcmp

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
#endif

namespace CGNSRead
{

typedef char char_33[33];

//------------------------------------------------------------------------------
class vtkCGNSArraySelection : public std::map<std::string, bool>
{
public:
  void Merge(const vtkCGNSArraySelection& other)
    {
    vtkCGNSArraySelection::const_iterator iter = other.begin();
    for (; iter != other.end(); ++iter)
      {
      (*this)[iter->first] = iter->second;
      }
    }

  void AddArray(const char* name, bool status=true)
    {
    (*this)[name] = status;
    }

  bool ArrayIsEnabled(const char* name)
    {
    vtkCGNSArraySelection::iterator iter = this->find(name);
    if (iter != this->end())
      {
      return iter->second;
      }

    // don't know anything about this array, enable it by default.
    return true;
    }

  bool HasArray(const char* name)
    {
    vtkCGNSArraySelection::iterator iter = this->find(name);
    return (iter != this->end());
    }

  int GetArraySetting(const char* name)
    {
    return this->ArrayIsEnabled(name)? 1 : 0;
    }

  void SetArrayStatus(const char* name, bool status)
    {
    this->AddArray(name, status);
    }

  const char* GetArrayName(int index)
  {
    int cc=0;
    for (vtkCGNSArraySelection::iterator iter = this->begin();
         iter != this->end(); ++iter)
      {

      if (cc==index)
        {
        return iter->first.c_str();
        }
      cc++;
      }
    return NULL;
  }

  int GetNumberOfArrays()
  {
    return static_cast<int>(this->size());
  }
};

//------------------------------------------------------------------------------
typedef struct
{
  int cnt; // 0  1 or 3
  int pos; // variable position in zone
  int xyzIndex;
  int owner_pos ;
  CGNS_ENUMT(DataType_t) dt;
  char_33 name;
} Variable;

//------------------------------------------------------------------------------
typedef struct
{
  int xyzIndex;
  bool isComponent;
  CGNS_ENUMT(DataType_t) dt;
  char_33 name;
} CGNSVariable ;

//------------------------------------------------------------------------------
typedef struct
{
  int numComp;
  char_33 name;
  int xyzIndex[3];
} CGNSVector;

//------------------------------------------------------------------------------
typedef struct
{
  bool isVector;
  int xyzIndex;
  char_33 name;
} VTKVariable;

//------------------------------------------------------------------------------
class ZoneInformation
{
public :
  char_33 name;
};

//------------------------------------------------------------------------------
class FamilyInformation
{
public:
  char_33 name;
  bool isBC;
};

//------------------------------------------------------------------------------
class BaseInformation
{
public:
  char_33 name;

  int cellDim;
  int physicalDim;
  //
  int baseNumber;

  std::vector<int>    steps;
  std::vector<double> times;

  // For unsteady meshes :
  // if useGridPointers == True:
  //    loadGridPointers for first zone
  //    and assume every zone use the same
  //    notation
  // else :
  //    assume only one grid is stored
  //    only first grid is read
  //
  // For unsteady flow
  // if useFlowPointers == True :
  //    same behavior as GridPointers
  // else if ( nstates > 1 ) :
  //    assume flow_solution are sorted
  //    to keep VisIt like behavior
  // else :
  //    only first solution is read
  //

  bool useGridPointers ; // for unsteady mesh
  bool useFlowPointers ; // for unsteady flow

  std::vector<CGNSRead::FamilyInformation> family;
  std::map< std::string, double> referenceState;

  int nzones;

  //std::vector<CGNSRead::zone> zone;
  vtkCGNSArraySelection PointDataArraySelection;
  vtkCGNSArraySelection CellDataArraySelection;
};

//------------------------------------------------------------------------------
class vtkCGNSMetaData
{
public:
  // Description:
  // quick parsing of cgns file to get interesting information
  // from a VTK point of view
  bool Parse(const char* cgnsFileName);

  // Description:
  // return number of base nodes
  int GetNumberOfBaseNodes()
    {
    return static_cast<int>(this->baseList.size());
    }

  // Description:
  // return const reference to a base information
  const CGNSRead::BaseInformation & GetBase(int numBase)
    {
    return this->baseList[numBase];
    }

  // Description
  // return reference to GlobalTime
  std::vector<double> & GetTimes()
    {
    return this->GlobalTime;
    }

  // Description:
  // print object debugging purpose
  void PrintSelf(std::ostream& os);
  
#ifdef PARAVIEW_USE_MPI
  void Broadcast ( vtkMultiProcessController* controller, int rank );
#endif  

  // Description
  // Constructor/Destructor
  vtkCGNSMetaData();
  ~vtkCGNSMetaData();

private:
  vtkCGNSMetaData(const vtkCGNSMetaData&);  // Not implemented.
  void operator=(const vtkCGNSMetaData&);  // Not implemented.

  std::vector<CGNSRead::BaseInformation> baseList;
  std::string LastReadFilename;
  // Not very elegant :
  std::vector<double> GlobalTime;
};

//------------------------------------------------------------------------------
// sort variables by name helper function
static int sortVariablesByName(const void* vOne, const void* vTwo)
{
  Variable* varOne = (Variable*) vOne;
  Variable* varTwo = (Variable*) vTwo;

  return (strcmp(varOne->name, varTwo->name));
}

//------------------------------------------------------------------------------
// compare name return true if name1 == name2
inline bool compareName(const char_33 nameOne ,const char_33 nameTwo)
{
  return (strncmp(nameOne , nameTwo , 32) == 0);
}

//------------------------------------------------------------------------------
// get vector from name
inline std::vector<CGNSVector>::iterator getVectorFromName(
  std::vector< CGNSVector >& vectorList, const char_33 name)
{
  for (std::vector<CGNSVector>::iterator iter = vectorList.begin();
       iter != vectorList.end(); ++iter)
    {
    if (strncmp (iter->name, name , 31) == 0)
      {
      return iter;
      }
    }
  return vectorList.end();
}

//------------------------------------------------------------------------------
inline bool isACGNSVariable(const std::vector<CGNSVariable>& varList,
                     const char_33 name)
{
  for (std::vector<CGNSVariable>::const_iterator iter = varList.begin();
       iter != varList.end(); ++iter)
    {
    if (strncmp(iter->name, name, 32) == 0)
      {
      return true;
      }
    }
  return false;
}

//------------------------------------------------------------------------------
void fillVectorsFromVars(std::vector< CGNSRead::CGNSVariable >&  vars,
                            std::vector< CGNSRead::CGNSVector >& vectors,
                            const int physicalDim);

}

#endif //__CGNS_READER_INTERNAL__
