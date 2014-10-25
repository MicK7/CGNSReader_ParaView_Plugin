// -*- c++ -*-
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCGNSReader.h

  Copyright (c) 2013-2014 Mickael Philit
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
// .NAME vtkCGNSReader -- reads a dataset in "CGNS" format
// .SECTION Description
// vtkCGNSReader creates a multi-block dataset and reads unstructured grids,
// and structured meshes from binary files stored in CGNS file format,
// with data stored at the nodes or at the cells.
//
// vtkCGNSReader is inspired by the VisIt CGNS reader originally written by
// B. Whitlock. vtkCGNSReader relies on the low level CGNS API to load DataSet
// and reduce memory footprint.
//
// .SECTION Caveats
//   ...
//
// .SECTION Thanks
// Thanks to .

#ifndef __vtkCGNSReader_h
#define __vtkCGNSReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkPVConfig.h"     // For PARAVIEW_USE_MPI

#include "vtkCGNSReaderInternal.h" // For parsing information request

class vtkDataArraySelection;
class vtkCallbackCommand;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class vtkCGNSReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkCGNSReader *New();
  vtkTypeMacro(vtkCGNSReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of CGNS datafile to read
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Is the given file name a CGNS file?
  int CanReadFile(const char* filename);
  

  // The following methods allow selective reading of solutions fields.
  int GetBaseArrayStatus(const char* name);
  void SetBaseArrayStatus(const char* name, int status);
  void DisableAllBases();
  void EnableAllBases();

  int GetNumberOfBaseArrays(); 
  int GetNumberOfPointArrays();
  int GetNumberOfCellArrays();

  const char* GetBaseArrayName(int index);
  const char* GetPointArrayName(int index);
  const char* GetCellArrayName(int index);
  
  int GetPointArrayStatus(const char* name);
  int GetCellArrayStatus(const char* name);
  
  void SetPointArrayStatus(const char* name, int status);
  void SetCellArrayStatus(const char* name, int status);

  void DisableAllPointArrays();
  void EnableAllPointArrays();

  void DisableAllCellArrays();
  void EnableAllCellArrays();

  vtkSetMacro(DoublePrecisionMesh,int);
  vtkGetMacro(DoublePrecisionMesh,int);
  vtkBooleanMacro(DoublePrecisionMesh,int);

  vtkSetMacro(LoadBndPatch,int);
  vtkGetMacro(LoadBndPatch,int);
  vtkBooleanMacro(LoadBndPatch,int);  
  
  vtkSetMacro(CreateEachSolutionAsBlock,int);
  vtkGetMacro(CreateEachSolutionAsBlock,int);
  vtkBooleanMacro(CreateEachSolutionAsBlock,int);
  
#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set/get the communication object used to relay a list of files
  // from the rank 0 process to all others. This is the only interprocess
  // communication required by vtkPExodusIIReader.
  void SetController(vtkMultiProcessController* c);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  // Description:
  // Sends metadata (that read from the input file, not settings modified
  // through this API) from the rank 0 node to all other processes in a job.
  void Broadcast( vtkMultiProcessController* ctrl );
#endif


protected:
  vtkCGNSReader();
  ~vtkCGNSReader();

  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
  virtual int RequestInformation(vtkInformation*,
                                 vtkInformationVector**,
                                 vtkInformationVector*);

  
  vtkDataArraySelection* BaseSelection;
  vtkDataArraySelection* PointDataArraySelection;
  vtkDataArraySelection* CellDataArraySelection;

  // The observer to modify this object when the array selections are
  // modified.
  vtkCallbackCommand* SelectionObserver;

  // Callback registered with the SelectionObserver.
  static void SelectionModifiedCallback(vtkObject* caller, unsigned long eid,
                                        void* clientdata, void* calldata);

  int GetCurvilinearZone(int  base, int zone,
                         int cell_dim, int phys_dim, cgsize_t *zsize,
                         vtkMultiBlockDataSet *mbase);

  int GetUnstructuredZone(int  base, int zone,
                          int cell_dim, int phys_dim, cgsize_t *zsize,
                          vtkMultiBlockDataSet *mbase);
#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  vtkIdType ProcRank;
  vtkIdType ProcSize;
#endif

  //BTX
  bool IsVarEnabled(CGNS_ENUMT(GridLocation_t) varcentering,
                    const CGNSRead::char_33 name);
  //ETX

  int getGridAndSolutionName(int base,
                             CGNSRead::char_33 GridCoordName, CGNSRead::char_33 SolutionName,
                             bool& readGridCoordName, bool& readSolutionName);
  
  int getCoordsIdAndFillRind(const CGNSRead::char_33 GridCoordName,
                             const int physicalDim, size_t& nCoordsArray,
                             std::vector<double>& gridChildId, int* rind);
  
  //BTX
  int getVarsIdAndFillRind(const double cgioSolId,
                           size_t& nVarArray, CGNS_ENUMT(GridLocation_t)& varCentering,
                           std::vector<double>& solChildId, int* rind);
  //ETX
  
  int fillArrayInformation(const std::vector<double>& solChildId,
                           const int physicalDim,
                           std::vector< CGNSRead::CGNSVariable >& cgnsVars,
                           std::vector< CGNSRead::CGNSVector >& cgnsVectors);
  //BTX
  int AllocateVtkArray(const int physicalDim, const vtkIdType nVals,
                       const CGNS_ENUMT ( GridLocation_t ) varCentering,
                       const std::vector< CGNSRead::CGNSVariable >& cgnsVars,
                       const std::vector< CGNSRead::CGNSVector >& cgnsVectors,
                       std::vector<vtkDataArray *>& vtkVars);
  //ETX
  int AttachReferenceValue(const int base, vtkDataSet* ds);
  
private:
  vtkCGNSReader(const vtkCGNSReader&);  // Not implemented.
  void operator=(const vtkCGNSReader&);  // Not implemented.

  CGNSRead::vtkCGNSMetaData Internal;  // Metadata

  char *FileName; // cgns file name
  int LoadBndPatch; // option to set section loading for unstructured grid
  int DoublePrecisionMesh; // option to set mesh loading to double precision
  int CreateEachSolutionAsBlock; // debug option to create
  
  // For internal cgio calls (low level IO)
  int cgioNum; // cgio file reference
  double rootId; // id of root node
  double currentId; // id of node currently being read (zone)
  //
  unsigned int NumberOfBases;
  int ActualTimeStep;
};

#endif // __vtkCGNSReader_h
