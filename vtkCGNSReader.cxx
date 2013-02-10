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

#include "vtkCGNSReader.h"

#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCallbackCommand.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkErrorCode.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTypeInt32Array.h"
#include "vtkTypeInt64Array.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVertex.h"
#include "vtkPolyhedron.h"
#include "vtkCellArray.h"
#include "vtkErrorCode.h"
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkInformationStringKey.h>

#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <map>
#include <string>

#include <tr1/array>

#include <vtksys/SystemTools.hxx>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
#endif

#include <cgns_io.h> // Low level IO
//#include <cgnslib.h> // CGNS_VERSION

vtkStandardNewMacro ( vtkCGNSReader );

//------------------------------------------------------------------------------
namespace CGNSRead
{
int GetVTKElemType (CGNS_ENUMT(ElementType_t) elemtype, bool &higherOrderWarning ,
                    bool &cgnsOrderFlag )
{
  int cell_type ;
  higherOrderWarning = false;
  cgnsOrderFlag = false;
  // Bring it into a clean GetVTKElemType implementation
  // see Xdmf
  switch (elemtype)
    {
    case CGNS_ENUMV(NODE):
      cell_type = VTK_VERTEX ;
      break;
    case CGNS_ENUMV(BAR_2):
      cell_type = VTK_LINE;
      break;
    case CGNS_ENUMV(BAR_3):
      cell_type = VTK_QUADRATIC_EDGE;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(BAR_4):
      cell_type = VTK_CUBIC_LINE;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(TRI_3):
      cell_type = VTK_TRIANGLE;
      break;
    case CGNS_ENUMV(TRI_6):
      cell_type = VTK_QUADRATIC_TRIANGLE;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(QUAD_4):
      cell_type = VTK_QUAD;
      break;
    case CGNS_ENUMV(QUAD_8):
      cell_type = VTK_QUADRATIC_QUAD;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(QUAD_9):
      cell_type = VTK_BIQUADRATIC_QUAD;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(TETRA_4):
      cell_type = VTK_TETRA;
      break;
    case CGNS_ENUMV(TETRA_10):
      cell_type = VTK_QUADRATIC_TETRA;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(PYRA_5):
      cell_type = VTK_PYRAMID;
      break;
    case CGNS_ENUMV(PYRA_14):
      cell_type = VTK_QUADRATIC_PYRAMID;
      higherOrderWarning = true;
      break;
    case CGNS_ENUMV(PENTA_6):
      cell_type = VTK_WEDGE;
      break;
    case CGNS_ENUMV(PENTA_15):
      cell_type = VTK_QUADRATIC_WEDGE;
      higherOrderWarning = true;
      cgnsOrderFlag = true;
      break;
    case CGNS_ENUMV(PENTA_18):
      cell_type = VTK_BIQUADRATIC_QUADRATIC_WEDGE;
      higherOrderWarning = true;
      cgnsOrderFlag = true;
      break;
    case CGNS_ENUMV(HEXA_8):
      cell_type = VTK_HEXAHEDRON;
      break;
    case CGNS_ENUMV(HEXA_20):
      cell_type = VTK_QUADRATIC_HEXAHEDRON;
      higherOrderWarning = true;
      cgnsOrderFlag = true;
      break;
    case CGNS_ENUMV(HEXA_27):
      cell_type = VTK_TRIQUADRATIC_HEXAHEDRON;
      higherOrderWarning = true;
      cgnsOrderFlag = true;
      break;
    default:
      cell_type = VTK_EMPTY_CELL;
      break;
    }
  return cell_type;
}

//------------------------------------------------------------------------------
void CGNS2VTKorder(const vtkIdType size, const int *cells_types,
                   vtkIdType *elements)
{
  const int maxPointsPerCells = 27;
  //static const int NULL_translate = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
  //                                  16,17,18,19,20,21,22,23,24,25,26};

  //CGNS --> VTK ordering of Elements
  static const int NODE_ToVTK[1]  = {0};

  static const int BAR_2_ToVTK[2] = {0,1};

  static const int BAR_3_ToVTK[3] = {0,1,2};

  static const int BAR_4_ToVTK[4] = {0,1,2,3};

  static const int TRI_3_ToVTK[3] = {0,1,2};

  static const int TRI_6_ToVTK[6] = {0,1,2,3,4,5};

  static const int QUAD_4_ToVTK[4] = {0,1,2,3};

  static const int QUAD_8_ToVTK[8] = {0,1,2,3,4,5,6,7};

  static const int QUAD_9_ToVTK[9] = {0,1,2,3,4,5,6,7,8};

  static const int TETRA_4_ToVTK[4] = {0,1,2,3};

  static const int TETRA_10_ToVTK[10] = {0,1,2,3,4,5,6,7,8,9};

  static const int PYRA_5_ToVTK[5] = {0,1,2,3,4};

  static const int PYRA_14_ToVTK[14] = {0,1,2,3,4,
                                        5,6,7,8,9,
                                        10,11,12,13};

  static const int PENTA_6_ToVTK[6] = {0,1,2,3,4,5};

  static const int PENTA_15_ToVTK[15] = {0,1,2,3,4,5,6,7,8,
                                         12,13,14,
                                         9,10,11};

  static const int PENTA_18_ToVTK[18] = {0,1,2,3,4,5,6,7,8,
                                         12,13,14,
                                         9,10,11,
                                         15,16,17};

  static const int HEXA_8_ToVTK[8] = {0,1,2,3,4,5,6,7};

  static const int HEXA_20_ToVTK[20] = {0,1,2,3,4,5,6,7,
                                        8,9,10,11,
                                        16,17,18,19,
                                        12,13,14,15};

  static const int HEXA_27_ToVTK[27] = {0,1,2,3,4,5,6,7,
                                        8,9,10,11,
                                        16,17,18,19,
                                        12,13,14,15,
                                        24,22,21,23,
                                        20,25,26};


  int tmp[maxPointsPerCells];
  const int *translator;
  vtkIdType pos = 0;
  for (vtkIdType icell = 0; icell < size; ++icell)
    {
    switch (cells_types[icell])
      {
      case VTK_VERTEX:
      case VTK_LINE:
      case VTK_QUADRATIC_EDGE:
      case VTK_CUBIC_LINE:
      case VTK_TRIANGLE:
      case VTK_QUADRATIC_TRIANGLE:
      case VTK_QUAD:
      case VTK_QUADRATIC_QUAD:
      case VTK_BIQUADRATIC_QUAD:
      case VTK_TETRA:
      case VTK_QUADRATIC_TETRA:
      case VTK_PYRAMID:
      case VTK_QUADRATIC_PYRAMID:
      case VTK_WEDGE:
        translator = NULL;
        break;
      case VTK_QUADRATIC_WEDGE:
        translator = PENTA_15_ToVTK;
        break;
      case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
        translator = PENTA_18_ToVTK;
        break;
      case VTK_HEXAHEDRON:
        translator = NULL;
        break;
      case VTK_QUADRATIC_HEXAHEDRON:
        translator = HEXA_20_ToVTK;
        break;
      case VTK_TRIQUADRATIC_HEXAHEDRON:
        translator = HEXA_27_ToVTK;
        break;
      default:
        translator = NULL;
        break;
      }
    vtkIdType numPointsPerCell = elements[pos];
    pos++;
    if (translator != NULL)
      {
      for (vtkIdType ip = 0; ip < numPointsPerCell; ++ip)
        {
        tmp[ip] = elements[translator[ip]+pos];
        }
      for (vtkIdType ip = 0; ip < numPointsPerCell; ++ip)
        {
        elements[pos+ip] = tmp[ip];
        }
      }
    pos += numPointsPerCell;
    }
}
//------------------------------------------------------------------------------
void CGNS2VTKorderMonoElem(const vtkIdType size, const int cell_type,
                           vtkIdType *elements)
{
  const int maxPointsPerCells = 27;
  //static const int NULL_translate = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
  //                                  16,17,18,19,20,21,22,23,24,25,26};

  //CGNS --> VTK ordering of Elements
  static const int NODE_ToVTK[1]  = {0};

  static const int BAR_2_ToVTK[2] = {0,1};

  static const int BAR_3_ToVTK[3] = {0,1,2};

  static const int BAR_4_ToVTK[4] = {0,1,2,3};

  static const int TRI_3_ToVTK[3] = {0,1,2};

  static const int TRI_6_ToVTK[6] = {0,1,2,3,4,5};

  static const int QUAD_4_ToVTK[4] = {0,1,2,3};

  static const int QUAD_8_ToVTK[8] = {0,1,2,3,4,5,6,7};

  static const int QUAD_9_ToVTK[9] = {0,1,2,3,4,5,6,7,8};

  static const int TETRA_4_ToVTK[4] = {0,1,2,3};

  static const int TETRA_10_ToVTK[10] = {0,1,2,3,4,5,6,7,8,9};

  static const int PYRA_5_ToVTK[5] = {0,1,2,3,4};

  static const int PYRA_14_ToVTK[14] = {0,1,2,3,4,
                                        5,6,7,8,9,
                                        10,11,12,13};

  static const int PENTA_6_ToVTK[6] = {0,1,2,3,4,5};

  static const int PENTA_15_ToVTK[15] = {0,1,2,3,4,5,6,7,8,
                                         12,13,14,
                                         9,10,11};

  static const int PENTA_18_ToVTK[18] = {0,1,2,3,4,5,6,7,8,
                                         12,13,14,
                                         9,10,11,
                                         15,16,17};

  static const int HEXA_8_ToVTK[8] = {0,1,2,3,4,5,6,7};

  static const int HEXA_20_ToVTK[20] = {0,1,2,3,4,5,6,7,
                                        8,9,10,11,
                                        16,17,18,19,
                                        12,13,14,15};

  static const int HEXA_27_ToVTK[27] = {0,1,2,3,4,5,6,7,
                                        8,9,10,11,
                                        16,17,18,19,
                                        12,13,14,15,
                                        24,22,21,23,
                                        20,25,26};

  int tmp[maxPointsPerCells];
  const int *translator;
  switch (cell_type)
    {
    case VTK_VERTEX:
    case VTK_LINE:
    case VTK_QUADRATIC_EDGE:
    case VTK_CUBIC_LINE:
    case VTK_TRIANGLE:
    case VTK_QUADRATIC_TRIANGLE:
    case VTK_QUAD:
    case VTK_QUADRATIC_QUAD:
    case VTK_BIQUADRATIC_QUAD:
    case VTK_TETRA:
    case VTK_QUADRATIC_TETRA:
    case VTK_PYRAMID:
    case VTK_QUADRATIC_PYRAMID:
    case VTK_WEDGE:
      translator = NULL;
      break;
    case VTK_QUADRATIC_WEDGE:
      translator = PENTA_15_ToVTK;
      break;
    case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
      translator = PENTA_18_ToVTK;
      break;
    case VTK_HEXAHEDRON:
      translator = NULL;
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      translator = HEXA_20_ToVTK;
      break;
    case VTK_TRIQUADRATIC_HEXAHEDRON:
      translator = HEXA_27_ToVTK;
      break;
    default:
      translator = NULL;
      break;
    }
  if (translator == NULL)
    {
    return;
    }

  vtkIdType pos = 0;
  for (vtkIdType icell = 0; icell < size; ++icell)
    {
    vtkIdType numPointsPerCell = elements[pos];
    pos++;
    for (vtkIdType ip = 0; ip < numPointsPerCell; ++ip)
      {
      tmp[ip] = elements[translator[ip]+pos];
      }
    for (vtkIdType ip = 0; ip < numPointsPerCell; ++ip)
      {
      elements[pos+ip] = tmp[ip];
      }
    pos += numPointsPerCell;
    }
}

}

//----------------------------------------------------------------------------
vtkCGNSReader::vtkCGNSReader()
{
  this->FileName  = NULL;

  this->LoadBndPatch = 0;
  this->NumberOfBases = 0;
  this->ActualTimeStep = 0;

  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->BaseSelection = vtkDataArraySelection::New();

  // Setup the selection callback to modify this object when an array
  // selection is changed.
  this->SelectionObserver = vtkCallbackCommand::New();
  this->SelectionObserver->SetCallback ( &vtkCGNSReader::SelectionModifiedCallback );
  this->SelectionObserver->SetClientData ( this );
  this->PointDataArraySelection->AddObserver ( vtkCommand::ModifiedEvent,
                                               this->SelectionObserver );
  this->CellDataArraySelection->AddObserver ( vtkCommand::ModifiedEvent,
                                              this->SelectionObserver );

  this->BaseSelection->AddObserver ( vtkCommand::ModifiedEvent,
                                     this->SelectionObserver );

  this->SetNumberOfInputPorts ( 0 );
  this->SetNumberOfOutputPorts ( 1 );

#ifdef PARAVIEW_USE_MPI
  this->ProcRank = 0;
  this->ProcSize = 1;
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
vtkCGNSReader::~vtkCGNSReader()
{
  this->SetFileName ( 0 );

  this->PointDataArraySelection->RemoveObserver ( this->SelectionObserver );
  this->PointDataArraySelection->Delete();
  this->CellDataArraySelection->RemoveObserver ( this->SelectionObserver );
  this->CellDataArraySelection->Delete();
  this->BaseSelection->RemoveObserver ( this->SelectionObserver );
  this->BaseSelection->Delete();
  this->SelectionObserver->Delete();

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif

}

#ifdef PARAVIEW_USE_MPI
//----------------------------------------------------------------------------
void vtkCGNSReader::SetController( vtkMultiProcessController* c )
{
  if (this->Controller == c)
    {
    return;
    }

  this->Modified();

  if (this->Controller)
    {
    this->Controller->UnRegister(this);
    }

  this->Controller = c;

  if (this->Controller)
    {
    this->Controller->Register(this);
    this->ProcRank = this->Controller->GetLocalProcessId();
    this->ProcSize = this->Controller->GetNumberOfProcesses();
    }

  if (! this->Controller || this->ProcSize <= 0)
    {
    this->ProcRank = 0;
    this->ProcSize = 1;
    }
}
#endif

//------------------------------------------------------------------------------
bool vtkCGNSReader::IsVarEnabled(CGNS_ENUMT(GridLocation_t) varcentering,
                                 const CGNSRead::char_33 name)
{
  vtkDataArraySelection *DataSelection = 0;
  if (varcentering == CGNS_ENUMV(Vertex))
    {
    DataSelection = this->PointDataArraySelection;
    }
  else
    {
    DataSelection = this->CellDataArraySelection;
    }

  return (DataSelection->ArrayIsEnabled( name ) != 0);
}

//------------------------------------------------------------------------------
int vtkCGNSReader::GetCurvilinearZone ( int fn, int  base, int zone,
                                        int cell_dim, int phys_dim,
                                        cgsize_t *zsize,
                                        vtkMultiBlockDataSet *mbase )
{
  bool err;
  int ngrids = 0;
  int rind[6];
  int i;

  int n;
  int ier;
  int ncoords = 0;

  // Check the number of grids stored in zone
  vtkDebugMacro( << "get number of grids in zone "
                 << zone << "\n" );

  if ( cg_ngrids ( fn, base, zone, &ngrids ) != CG_OK )
    {
    vtkErrorMacro( << "Could not get number of grids in zone "
                   << zone << "\n"
                   << cg_get_error() );
    return 1;
    }

  if ( ngrids < 1 ){
    vtkErrorMacro ( << "No Grid found in zone " << zone << "\n" );
    return 1;
    }

  CGNS_ENUMT(DataType_t) MeshType = CGNS_ENUMV(RealDouble) ;//TODO modify to get user/pipe preference

  // Get the number of Coordinates in GridCoordinates node
  if ( cg_ncoords ( fn, base, zone, &ncoords ) != CG_OK )
    {
    vtkErrorMacro( << "\t\tCould not get the number of coords\n"
                   <<  cg_get_error() << "\n" );
    return 1;
    }

  if ( cg_goto ( fn, base, "Zone_t", zone,
                 "GridCoordinates_t", 1, "end" ) )
    {
    vtkErrorMacro( << "No GridCoordinates_t\n");
    return 1;
    }
  // check for rind
  if ( ( ier = cg_rind_read ( rind ) ) != CG_OK )
    {
    if ( ier != CG_NODE_NOT_FOUND )
      {
      vtkErrorMacro(<< "Unexpected problem while reading rind information\n");
      return 1;
      }
    for ( n = 0; n < 6; n++ )
      {
      rind[n] = 0;
      }
    }

  // source layout
  cgsize_t s_start[3]  = {1,1,1};
  cgsize_t s_stride[3] = {1,1,1};
  cgsize_t s_end[3];
  // memory destination layout
  cgsize_t m_start[3]  = {1,1,1} ;
  cgsize_t m_stride[3] = {3,1,1};
  cgsize_t m_end[3] = {1,1,1};
  cgsize_t m_dims[3] = {1,1,1};

  vtkIdType nPts = 0;

  // m_num_dims = cell_dim;

  // get grid coordinate range
  for ( n = 0; n < cell_dim; n++ )
    {
    s_start[n] = rind[2*n] + 1;
    s_end[n]   = rind[2*n] + zsize[n];
    m_end[n]   = zsize[n];
    m_dims[n]  = zsize[n];
    }
  // Compute number of points
  nPts = (vtkIdType) ( m_end[0] * m_end[1] * m_end[2] );

  // Set up points
  vtkPoints *points  = vtkPoints::New();

  // Populate the points array
  int extent[6] = {0,0,0,0,0,0};
  extent[1] = m_end[0]-1;
  extent[3] = m_end[1]-1;
  extent[5] = m_end[2]-1;
  
  // wacky hack ...
  m_end[0] *= 3; //for aliasing
  
  //
  // vtkPoints assumes float data type
  //
  if ( MeshType == CGNS_ENUMV(RealDouble) ){
    points->SetDataTypeToDouble();
    }
  //
  // Resize vtkPoints to fit data
  //
  points->SetNumberOfPoints ( nPts );


  //----------------------------------------------------------------------------
  // Get Coordinates and FlowSolution node names
  bool readGridCoordName = true;
  bool readSolutionName = true;
  CGNSRead::char_33 GridCoordName;
  CGNSRead::char_33 SolutionName;

  if ((this->Internal.GetBase(base-1).useGridPointers == true) ||
      (this->Internal.GetBase(base-1).useFlowPointers == true))
    {
    CGNSRead::char_33 nameZoneIter;
    size_t ptSize = 32*this->Internal.GetBase(base-1).steps.size()+1;
    char *pointers = new char[ptSize];

    cg_ziter_read(fn, base, zone, nameZoneIter);
    if ( cg_goto(fn,base, "Zone_t",zone,nameZoneIter,0,NULL) == CG_OK )
      {
      int arraycount,arrayno;
      cg_narrays( &arraycount );
      for( arrayno = 1; arrayno<=arraycount; arrayno++ )
        {
        CGNSRead::char_33 arrayname;
        CGNS_ENUMT(DataType_t) arraydt;
        int dimension;
        cgsize_t elements[ 12 ];
        cg_array_info(arrayno,arrayname,&arraydt,&dimension,elements);
        if( strcmp( arrayname,"GridCoordinatesPointers" ) == 0 )
          {
          cg_array_read( arrayno, pointers ) ;
          strncpy(GridCoordName, &pointers[this->ActualTimeStep*32], 32);
          GridCoordName[32]='\0';
          readGridCoordName = false;
          }

        if( strcmp( arrayname,"FlowSolutionPointers" ) == 0 )
          {
          cg_array_read( arrayno, pointers ) ;
          strncpy(SolutionName, &pointers[this->ActualTimeStep*32], 32);
          SolutionName[32]='\0';
          readSolutionName = false;
          }

        }
      }
    else
      {
      strcpy(GridCoordName, "GridCoordinates");
      strcpy(SolutionName, "FlowSolution");
      }
    delete [] pointers;
    }
  //----------------------------------------------------------------------------

  if ( readGridCoordName )
    {
    int requiredgrid = 1;
    cg_grid_read ( fn, base, zone, requiredgrid, GridCoordName );
    }

  int narrays=0;
  // std::cout << "Reading mesh node " << GridCoordName << std::endl;

  if ( cg_goto ( fn, base, "Zone_t", zone, GridCoordName, 0, "end" ) != CG_OK )
    {
    vtkErrorMacro( << cg_get_error() );
    return 1;
    }

  cg_narrays ( &narrays );
  if ( narrays < ncoords )
    {
    vtkErrorMacro( << "Not enough coordinates in node "
                   << GridCoordName << "\n" );
    return 1;
    }

  // Get GridCoordinate node ID for low level access
  double gridid, coordid ;
  cgio_get_node_id ( this->cgio_desc, this->currentid, GridCoordName, &gridid );

  //
  // Populate the coordinates.  Put in 3D points with z=0 if the mesh is 2D.
  //
  if ( MeshType == CGNS_ENUMV(RealDouble) ) // DOUBLE PRECISION MESHPOINTS
    {
    double *coords = (double *) points->GetVoidPointer ( 0 ) ;
    double *currentCoord = (double *) &(coords[0]);

    CGNSRead::char_33 coordname;
    size_t len;
    CGNS_ENUMT(DataType_t) ct;

    memset ( coords, 0, 3*nPts*sizeof ( double ) );

    for ( int c = 1; c <= ncoords; ++c )
      {
      if ( cg_coord_info ( fn, base, zone, c, &ct,
                           coordname ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        break;
        }
      // Determine direction X,Y,Z
      len = strlen ( coordname ) - 1;
      switch ( coordname[len] )
        {
        case 'X':
          currentCoord = (double *) &(coords[0]) ;
          break;
        case 'Y':
          currentCoord = (double *) &(coords[1]) ;
          break;
        case 'Z':
          currentCoord = (double *) &(coords[2]) ;
          break;
        }

      cgio_get_node_id ( this->cgio_desc, gridid, coordname, &coordid );

      // quick transfer of data if same data types
      if ( ct == CGNS_ENUMV(RealDouble) )
        {
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, cell_dim , m_end,
                              m_start, m_end, m_stride, (void *) currentCoord ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( << "cgio_read_data :" << message );
          }
        }
      else
        {
        float *array_data = 0;
        const cgsize_t m_nostride[3] = {1,1,1};

        if ( ct != CGNS_ENUMV(RealSingle) )
          {
          vtkErrorMacro( << "Invalid datatype for GridCoordinates\n" );
          break;
          }

        // need to read into temp array to convert data
        array_data = new float[nPts];
        if ( array_data == 0 )
          {
          vtkErrorMacro ( "Error allocating buffer array" );
          break;
          }
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, cell_dim, m_dims,
                              m_start, m_dims, m_nostride, ( void* ) array_data ) )
          {
          delete [] array_data;
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "Buffer array cgio_read_data :" << message );
          break;
          }

        for (vtkIdType ii = 0; ii < nPts; ii++ )
          {
          currentCoord[m_stride[0]*ii] = ( double ) array_data[ii];
          }

        delete [] array_data;
        }
      }
    }
  else  // SINGLE PRECISION MESHPOINTS
    {
    float *coords = (float *) points->GetVoidPointer ( 0 ) ;
    float *currentCoord = (float *) &(coords[0]);

    CGNSRead::char_33 coordname;
    size_t len;
    CGNS_ENUMT(DataType_t) ct;

    memset ( coords, 0, 3*nPts*sizeof ( float ) );

    for ( int c = 1; c <= ncoords; ++c )
      {
      if ( cg_coord_info ( fn, base, zone, c, &ct,
                           coordname ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        break;
        }
      // Determine direction X,Y,Z
      len = strlen ( coordname ) - 1;
      switch ( coordname[len] )
        {
        case 'X':
          currentCoord = (float *) &(coords[0]) ;
          break;
        case 'Y':
          currentCoord = (float *) &(coords[1]) ;
          break;
        case 'Z':
          currentCoord = (float *) &(coords[2]) ;
          break;
        }

      cgio_get_node_id ( this->cgio_desc, gridid, coordname, &coordid );

      // quick transfer of data if same data types
      if ( ct == CGNS_ENUMV(RealSingle) )
        {
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, cell_dim , m_end,
                              m_start, m_end, m_stride, (void *) currentCoord ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "cgio_read_data :" << message );
          }
        }
      else
        {
        double *array_data = 0;
        const cgsize_t m_nostride[3] = {1,1,1};

        if ( ct != CGNS_ENUMV(RealDouble) )
          {
          vtkErrorMacro( << "Invalid datatype for GridCoordinates\n" );
          break;
          }

        // need to read into temp array to convert data
        array_data = new double[nPts];
        if ( array_data == 0 )
          {
          vtkErrorMacro ( << "Error allocating buffer array" );
          break;
          }
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, cell_dim, m_dims,
                              m_start, m_dims, m_nostride, (void *) array_data ) )
          {
          delete [] array_data;
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( << "Buffer array cgio_read_data :" << message );
          break;
          }

        for (vtkIdType ii = 0; ii < nPts; ii++ )
          {
          currentCoord[m_stride[0]*ii] = (float) array_data[ii];
          }

        delete [] array_data;
        }
      }
    }
  cgio_release_id ( this->cgio_desc, gridid );

  //----------------------------------------------------------------------------
  // Handle solutions
  //----------------------------------------------------------------------------
  int nsols=0;
  ier = cg_nsols ( fn, base, zone, &nsols );

  vtkMultiBlockDataSet* mzone = vtkMultiBlockDataSet::New();

  bool CreateEachSolutionAsBlock = false ; // Debugging mode !
  
  if ((nsols > 0) && (readSolutionName != true))
    {
    int sol=1;
    CGNS_ENUMT(GridLocation_t) varcentering= CGNS_ENUMV(Vertex);
    double sol_cgio_id, var_cgio_id ;

    if ( cg_goto ( fn, base, "Zone_t", zone, SolutionName, 0, "end" ) != CG_OK )
      {
      vtkErrorMacro( << cg_get_error() );
      return 1;
      }
    cg_gridlocation_read( &varcentering );
    {
    int fileNumber;
    int curBase;
    int curDepth = 2;

    CGNSRead::char_33 curLabel[2];
    int numList[2];
    char ** llabel = new char*[2];
    llabel[0] = (char *) &curLabel[0];
    llabel[1] = (char *) &curLabel[1];

    cg_where(&fileNumber, &curBase, &curDepth, llabel, numList);
    sol = numList[1];
    delete [] llabel;
    }
    cgio_get_node_id ( this->cgio_desc, this->currentid, SolutionName, &sol_cgio_id );

    vtkStructuredGrid *sgrid   = vtkStructuredGrid::New();
    sgrid->SetExtent ( extent );
    sgrid->SetPoints ( points );

    bool skip = false ;
    while( skip != true )
      {

      if (( varcentering != CGNS_ENUMV(Vertex)) &&
          ( varcentering != CGNS_ENUMV(CellCenter)))
        {
        vtkWarningMacro( << "Solution " << SolutionName <<
                         " centering is not supported");

        skip =true;
        break;
        }
      vtkDebugMacro("Reading solution :" << SolutionName << "\n" );
      int nfields ;
      CGNSRead::char_33 fieldname ;

      ier = cg_nfields ( fn, base, zone, sol, &nfields );

      std::vector< CGNSRead::CGNSVariable > vars ( nfields ) ;
      std::vector< CGNSRead::CGNSVector > listOfVector ;

      // Read variable names
      for ( int f = 0; f < nfields; ++f )
        {
        cg_field_info(fn, base, zone, sol, f+1,
                      &(vars[f].dt), vars[f].name);
        vars[f].isComponent = false;
        vars[f].xyzIndex = 0;
        }
      // Create vector name from available variable
      // when VarX, VarY, VarZ is detected
      CGNSRead::fillVectorsFromVars(vars, listOfVector, phys_dim);

      //------------------------------------------------------------
      if ( cg_goto ( fn, base, "Zone_t", zone,"FlowSolution_t", sol, "end" ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        return 1;
        }

      // check for rind
      if ( ( ier = cg_rind_read ( rind ) ) != CG_OK )
        {
        if ( ier != CG_NODE_NOT_FOUND ){
          vtkErrorMacro( << "Unexpected error while reading rind"
                         "information in solution node\n" );
          return 1;
          }
        for ( n = 0; n < 6; n++ ){
          rind[n] = 0;
          }
        }

      // source
      cgsize_t fs_start[3]  = {1,1,1};
      cgsize_t fs_stride[3] = {1,1,1};
      cgsize_t fs_end[3];
      //
      int m_num_dims;
      // destination
      cgsize_t fm_start[3]  = {1,1,1} ;
      cgsize_t fm_stride[3] = {1,1,1};
      cgsize_t fm_end[3] = {1,1,1};
      cgsize_t fm_dims[3] = {1,1,1};

      vtkIdType nVals = 0;

      m_num_dims = cell_dim;

      // get solution data range
      int nsc = varcentering == CGNS_ENUMV ( Vertex ) ? 0 : cell_dim;

      for ( n = 0; n < cell_dim; ++n )
        {
        fs_start[n] = rind[2*n] + 1;
        fs_end[n]   = rind[2*n] + zsize[n+nsc];
        fm_end[n]   = zsize[n+nsc];
        fm_dims[n]  = zsize[n+nsc];
        }

      // compute number of field values
      nVals = (vtkIdType) ( fm_end[0] * fm_end[1] * fm_end[2] );
      //---------------------------------------------------------
      // VECTORS aliasing ...
      // destination
      cgsize_t fvm_start[3]  = {1,1,1} ;
      cgsize_t fvm_stride[3] = {3,1,1};
      cgsize_t fvm_end[3] = {1,1,1};
      cgsize_t fvm_dims[3] = {1,1,1};

      fvm_stride[0] = ( cgsize_t ) phys_dim;

      // Sioux ...
      fvm_dims[0] = fm_dims[0]*fvm_stride[0]; //*phys_dim to support 2D aliasing
      fvm_dims[1] = fm_dims[1];
      fvm_dims[2] = fm_dims[2];
      fvm_end[0] = fm_end[0]*fvm_stride[0]; //*phys_dim for aliasing
      fvm_end[1] = fm_end[1];
      fvm_end[2] = fm_end[2];

      //----------------------------------------------------------
      //
      // Count number of vars and vectors
      // Assign vars and vectors to a vtkvars array
      std::vector<vtkDataArray *> vtkvars ( nfields );

      for ( int f = 0; f < nfields; f++ )
        {
        vtkvars[f] = 0 ;

        if ( vars[f].isComponent == false )
          {
          if (IsVarEnabled( varcentering, vars[f].name ) == false)
            {
            continue;
            }

          switch ( vars[f].dt )
            {
            // Other case to handle
            case CGNS_ENUMV(Integer):
              vtkvars[f] = vtkIntArray::New();
              break;
            case CGNS_ENUMV(RealSingle):
              vtkvars[f] = vtkFloatArray::New();
              break;
            case CGNS_ENUMV(RealDouble):
              vtkvars[f] = vtkDoubleArray::New();
              break;
            case CGNS_ENUMV(Character):
              vtkvars[f] = vtkCharArray::New();
              break;
            }
          vtkvars[f]->SetName ( vars[f].name );
          vtkvars[f]->SetNumberOfComponents ( 1 );
          vtkvars[f]->SetNumberOfTuples ( nVals );
          }

        }

      for ( std::vector<CGNSRead::CGNSVector>::iterator it = listOfVector.begin();
            it != listOfVector.end(); ++it )
        {
        vtkDataArray *arr = 0;

        if (IsVarEnabled( varcentering, it->name ) == false){
          continue;
          }

        int nv = it->xyzIndex[0];
        switch ( vars[ nv ].dt )
          {
          // TODO Other cases to be done
          case CGNS_ENUMV(Integer):
            arr = vtkIntArray::New();
            break;
          case CGNS_ENUMV(RealSingle):
            arr = vtkFloatArray::New();
            break;
          case CGNS_ENUMV(RealDouble):
            arr = vtkDoubleArray::New();
            break;
          case CGNS_ENUMV(Character):
            arr = vtkCharArray::New();
            break;
          }

        arr->SetName ( it->name );
        arr->SetNumberOfComponents ( phys_dim );
        arr->SetNumberOfTuples ( nVals );

        for ( int dim=0; dim< phys_dim; dim++ )
          {
          arr->SetComponentName ( ( vtkIdType ) ( dim ) , vars[ it->xyzIndex[dim] ].name );
          vtkvars[ it->xyzIndex[dim] ] = arr ;
          }
        }

      // Load Data :
      for ( int f = 0; f < nfields; ++f )
        {
        // only read allocated fields
        if ( vtkvars[f] == 0 )
          {
          continue;
          }

        cgio_get_node_id ( this->cgio_desc, sol_cgio_id, vars[f].name, &var_cgio_id );

        // quick transfer of data because data types is given by cgns database
        if ( vars[f].isComponent == false )
          {
          if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                                fs_start, fs_end, fs_stride, cell_dim, fm_dims,
                                fm_start, fm_end, fm_stride,
                                ( void* ) vtkvars[f]->GetVoidPointer ( 0 ) ) != 0 )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          }
        else
          {
          // quick transfer of data because data types is given by cgns database
          if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                                fs_start, fs_end, fs_stride, cell_dim, fvm_dims,
                                fvm_start, fvm_end, fvm_stride,
                                ( void* ) vtkvars[f]->GetVoidPointer ( vars[f].xyzIndex-1 ) ) )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }

          }
        cgio_release_id ( this->cgio_desc, var_cgio_id );
        }
      cgio_release_id ( this->cgio_desc, sol_cgio_id );

      // Append data to StructuredGrid
      vtkDataSetAttributes* dsa = 0;
      if ( varcentering == CGNS_ENUMV ( Vertex ) ) //ON_NODES
        {
        dsa = sgrid->GetPointData();
        }
      if ( varcentering == CGNS_ENUMV ( CellCenter ) ) //ON_CELL
        {
        dsa = sgrid->GetCellData();
        }
      // SetData in vtk Structured Zone + Clean Pointers
      for ( int nv = 0; nv < nfields; nv++ )
        {
        // only transfer allocated fields
        if ( vtkvars[nv] == 0 ){
          continue;
          }

        if ( ( vars[nv].isComponent == false ) ||
             ( vars[nv].xyzIndex == 1 ) )
          {
          dsa->AddArray ( vtkvars[nv] );
          vtkvars[nv]->Delete();
          }
        vtkvars[nv] = 0 ;
        }
      skip = true;
      }
    // Handle Ref Values
    const std::map< std::string, double>& arrState = this->Internal.GetBase(base-1).referenceState ;
    std::map< std::string, double>::const_iterator iteRef = arrState.begin();
    for ( iteRef = arrState.begin(); iteRef != arrState.end(); iteRef++)
      {
      vtkDoubleArray* refValArray = vtkDoubleArray::New();
      refValArray->SetNumberOfComponents(1);
      refValArray->SetName( iteRef->first.c_str() );
      refValArray->InsertNextValue( iteRef->second );
      sgrid->GetFieldData()->AddArray( refValArray );
      refValArray->Delete();
      }
    //
    mbase->SetBlock ( ( zone-1 ), sgrid );
    sgrid->Delete();

    }
  else if ( (nsols > 0)  && (! CreateEachSolutionAsBlock) )
    {
    vtkStructuredGrid *sgrid   = vtkStructuredGrid::New();
    sgrid->SetExtent ( extent );
    sgrid->SetPoints ( points );
    
    int requiredsol = 1;
    int cell_solution = 0;
    int point_solution = 0;

    for (int sol=1; sol<=nsols; sol++ )
      {
      CGNSRead::char_33 solname;
      CGNS_ENUMT(GridLocation_t) varcentering;
      double sol_cgio_id, var_cgio_id ;

      cg_sol_info ( fn, base, zone, sol, solname, &varcentering );
      cgio_get_node_id ( this->cgio_desc, this->currentid, solname, &sol_cgio_id );

      bool skip = false ;

      if ( varcentering != CGNS_ENUMV ( Vertex ) )
        {
        point_solution++;
        skip = ( point_solution != requiredsol ) ;
        }
      else if ( varcentering != CGNS_ENUMV ( CellCenter ) )
        {
        cell_solution++;
        skip = ( cell_solution != requiredsol ) ;
        }
      else
        {
        vtkWarningMacro( << "Solution " << solname <<
                         " centering is not supported");
        skip = true;
        }

      if (skip)
        {
        continue;
        }

      int nfields ;
      CGNSRead::char_33 fieldname ;

      ier = cg_nfields ( fn, base, zone, sol, &nfields );

      std::vector< CGNSRead::CGNSVariable > vars ( nfields ) ;
      std::vector< CGNSRead::CGNSVector > listOfVector ;

      // Read variable names
      for ( int f = 0; f < nfields; ++f )
        {
        cg_field_info(fn, base, zone, sol, f+1,
                      &(vars[f].dt), vars[f].name);
        vars[f].isComponent = false;
        vars[f].xyzIndex = 0;
        }
      // Create vector name from available variable
      // when VarX, VarY, VarZ is detected
      CGNSRead::fillVectorsFromVars(vars, listOfVector, phys_dim);

      //------------------------------------------------------------
      if ( cg_goto ( fn, base, "Zone_t", zone,"FlowSolution_t", sol, "end" ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        return 1;
        }

      // check for rind
      if ( ( ier = cg_rind_read ( rind ) ) != CG_OK )
        {
        if ( ier != CG_NODE_NOT_FOUND ){
          vtkErrorMacro( << "Unexpected error while reading rind"
                         "information in solution node\n" );
          return 1;
          }
        for ( n = 0; n < 6; n++ ){
          rind[n] = 0;
          }
        }

      // source
      cgsize_t fs_start[3]  = {1,1,1};
      cgsize_t fs_stride[3] = {1,1,1};
      cgsize_t fs_end[3];
      //
      int m_num_dims;
      // destination
      cgsize_t fm_start[3]  = {1,1,1} ;
      cgsize_t fm_stride[3] = {1,1,1};
      cgsize_t fm_end[3] = {1,1,1};
      cgsize_t fm_dims[3] = {1,1,1};

      vtkIdType nVals = 0;

      m_num_dims = cell_dim;

      // get solution data range
      int nsc = varcentering == CGNS_ENUMV ( Vertex ) ? 0 : cell_dim;

      for ( n = 0; n < cell_dim; ++n )
        {
        fs_start[n] = rind[2*n] + 1;
        fs_end[n]   = rind[2*n] + zsize[n+nsc];
        fm_end[n]   = zsize[n+nsc];
        fm_dims[n]  = zsize[n+nsc];
        }

      // compute number of field values
      nVals = (vtkIdType) ( fm_end[0] * fm_end[1] * fm_end[2] );
      //---------------------------------------------------------
      // VECTORS aliasing ...
      // destination
      cgsize_t fvm_start[3]  = {1,1,1} ;
      cgsize_t fvm_stride[3] = {3,1,1};
      cgsize_t fvm_end[3] = {1,1,1};
      cgsize_t fvm_dims[3] = {1,1,1};

      fvm_stride[0] = ( cgsize_t ) phys_dim;

      // Sioux ...
      fvm_dims[0] = fm_dims[0]*fvm_stride[0]; //*phys_dim to support 2D aliasing
      fvm_dims[1] = fm_dims[1];
      fvm_dims[2] = fm_dims[2];
      fvm_end[0] = fm_end[0]*fvm_stride[0]; //*phys_dim for aliasing
      fvm_end[1] = fm_end[1];
      fvm_end[2] = fm_end[2];

      //----------------------------------------------------------
      //
      // Count number of vars and vectors
      // Assign vars and vectors to a vtkvars array
      std::vector<vtkDataArray *> vtkvars ( nfields );

      for ( int f = 0; f < nfields; f++ )
        {
        vtkvars[f] = 0 ;

        if ( vars[f].isComponent == false )
          {
          if (IsVarEnabled( varcentering, vars[f].name ) == false)
            {
            continue;
            }

          switch ( vars[f].dt )
            {
            // Other case to handle
            case CGNS_ENUMV(Integer):
              vtkvars[f] = vtkIntArray::New();
              break;
            case CGNS_ENUMV(RealSingle):
              vtkvars[f] = vtkFloatArray::New();
              break;
            case CGNS_ENUMV(RealDouble):
              vtkvars[f] = vtkDoubleArray::New();
              break;
            case CGNS_ENUMV(Character):
              vtkvars[f] = vtkCharArray::New();
              break;
            }
          vtkvars[f]->SetName ( vars[f].name );
          vtkvars[f]->SetNumberOfComponents ( 1 );
          vtkvars[f]->SetNumberOfTuples ( nVals );
          }

        }

      for ( std::vector<CGNSRead::CGNSVector>::iterator it = listOfVector.begin();
            it != listOfVector.end(); ++it )
        {
        vtkDataArray *arr = 0;

        if (IsVarEnabled( varcentering, it->name ) == false){
          continue;
          }

        int nv = it->xyzIndex[0];
        switch ( vars[ nv ].dt )
          {
          // TODO Other cases to be done
          case CGNS_ENUMV(Integer):
            arr = vtkIntArray::New();
            break;
          case CGNS_ENUMV(RealSingle):
            arr = vtkFloatArray::New();
            break;
          case CGNS_ENUMV(RealDouble):
            arr = vtkDoubleArray::New();
            break;
          case CGNS_ENUMV(Character):
            arr = vtkCharArray::New();
            break;
          }

        arr->SetName ( it->name );
        arr->SetNumberOfComponents ( phys_dim );
        arr->SetNumberOfTuples ( nVals );

        for ( int dim=0; dim< phys_dim; dim++ )
          {
          arr->SetComponentName ( ( vtkIdType ) ( dim ) , vars[ it->xyzIndex[dim] ].name );
          vtkvars[ it->xyzIndex[dim] ] = arr ;
          }
        }

      // Load Data :
      for ( int f = 0; f < nfields; ++f )
        {
        // only read allocated fields
        if ( vtkvars[f] == 0 ){
          continue;
          }

        cgio_get_node_id ( this->cgio_desc, sol_cgio_id, vars[f].name, &var_cgio_id );

        // quick transfer of data because data types is given by cgns database
        if ( vars[f].isComponent == false )
          {
          if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                                fs_start, fs_end, fs_stride, cell_dim, fm_dims,
                                fm_start, fm_end, fm_stride,
                                ( void* ) vtkvars[f]->GetVoidPointer ( 0 ) ) != 0 )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          }
        else
          {
          // quick transfer of data because data types is given by cgns database
          if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                                fs_start, fs_end, fs_stride, cell_dim, fvm_dims,
                                fvm_start, fvm_end, fvm_stride,
                                ( void* ) vtkvars[f]->GetVoidPointer ( vars[f].xyzIndex-1 ) ) )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }

          }
        cgio_release_id ( this->cgio_desc, var_cgio_id );
        }
      cgio_release_id ( this->cgio_desc, sol_cgio_id );

      // Append data to StructuredGrid
      vtkDataSetAttributes* dsa = 0;
      if ( varcentering == CGNS_ENUMV ( Vertex ) ) //ON_NODES
        {
        dsa = sgrid->GetPointData();
        }
      if ( varcentering == CGNS_ENUMV ( CellCenter ) ) //ON_CELL
        {
        dsa = sgrid->GetCellData();
        }
      // SetData in vtk Structured Zone + Clean Pointers
      for ( int nv = 0; nv < nfields; nv++ )
        {
        // only transfer allocated fields
        if ( vtkvars[nv] == 0 ){
          continue;
          }

        if ( ( vars[nv].isComponent == false ) ||
             ( vars[nv].xyzIndex == 1 ) )
          {
          dsa->AddArray ( vtkvars[nv] );
          vtkvars[nv]->Delete();
          }
        vtkvars[nv] = 0 ;
        }
      }
    // Handle Ref Values
    const std::map< std::string, double>& arrState = this->Internal.GetBase(base-1).referenceState ;
    std::map< std::string, double>::const_iterator iteRef = arrState.begin();
    for ( iteRef = arrState.begin(); iteRef != arrState.end(); iteRef++)
      {
      vtkDoubleArray* refValArray = vtkDoubleArray::New();
      refValArray->SetNumberOfComponents(1);
      refValArray->SetName( iteRef->first.c_str() );
      refValArray->InsertNextValue( iteRef->second );
      sgrid->GetFieldData()->AddArray( refValArray );
      refValArray->Delete();
      }
    //
    mbase->SetBlock ( ( zone-1 ), sgrid );
    sgrid->Delete();
    }
  else if ( (nsols > 0 ) && (CreateEachSolutionAsBlock) )
    {
    mzone->SetNumberOfBlocks ( nsols );
    for (int sol=1; sol<=nsols; sol++ )
      {
      CGNSRead::char_33 solname;
      CGNS_ENUMT(GridLocation_t) varcentering;
      double sol_cgio_id, var_cgio_id ;

      cg_sol_info ( fn, base, zone, sol, solname, &varcentering );
      cgio_get_node_id ( this->cgio_desc, this->currentid, solname, &sol_cgio_id );

      mzone->GetMetaData((unsigned int)sol - 1)->Set(vtkCompositeDataSet::NAME(), solname);

      if ( varcentering != CGNS_ENUMV ( Vertex ) &&
           varcentering != CGNS_ENUMV ( CellCenter ) )
        {
        vtkWarningMacro( << "Solution " << solname << " centering is not supported");
        }

      vtkStructuredGrid *sgrid   = vtkStructuredGrid::New();
      sgrid->SetExtent ( extent );
      sgrid->SetPoints ( points );
      //mzone->SetBlock ( ( sol-1 ), sgrid );
      //sgrid->Delete();

      int nfields ;
      CGNSRead::char_33 fieldname ;

      ier = cg_nfields ( fn, base, zone, sol, &nfields );

      std::vector< CGNSRead::CGNSVariable > vars ( nfields ) ;
      std::vector< CGNSRead::CGNSVector > listOfVector ;

      // Read variable names
      for ( int f = 0; f < nfields; ++f )
        {
        cg_field_info ( fn, base, zone, sol, f+1, & ( vars[f].dt ), vars[f].name );

        vars[f].isComponent = false;
        vars[f].xyzIndex = 0;
        }
      // Create vector name from available variable
      // when VarX, VarY, VarZ is detected
      CGNSRead::fillVectorsFromVars(vars, listOfVector, phys_dim);

      //-----------------------------------------
      if ( cg_goto ( fn, base, "Zone_t", zone,"FlowSolution_t", sol, "end" ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        return 1;
        }

      /* check for rind */
      if ( ( ier = cg_rind_read ( rind ) ) != CG_OK )
        {
        if ( ier != CG_NODE_NOT_FOUND ){
          vtkErrorMacro( << "Unexpected error while reading rind information in solution node\n" );
          return 1;
          }
        for ( n = 0; n < 6; n++ ){
          rind[n] = 0;
          }
        }

      // source
      cgsize_t fs_start[3]  = {1,1,1};
      cgsize_t fs_stride[3] = {1,1,1};
      cgsize_t fs_end[3];
      //
      int m_num_dims;
      // destination
      cgsize_t fm_start[3]  = {1,1,1} ;
      cgsize_t fm_stride[3] = {1,1,1};
      cgsize_t fm_end[3] = {1,1,1};
      cgsize_t fm_dims[3] = {1,1,1};

      vtkIdType nVals = 0;

      m_num_dims = cell_dim;

      // get solution data range
      int nsc = varcentering == CGNS_ENUMV ( Vertex ) ? 0 : cell_dim;

      for ( n = 0; n < cell_dim; ++n )
        {
        fs_start[n] = rind[2*n] + 1;
        fs_end[n]   = rind[2*n] + zsize[n+nsc];
        fm_end[n]   = zsize[n+nsc];
        fm_dims[n]  = zsize[n+nsc];
        }

      // compute number of field values
      nVals = (vtkIdType) ( fm_end[0] * fm_end[1] * fm_end[2] );
      //----------------------------------------------------
      // VECTORS aliasing ...
      // destination
      cgsize_t fvm_start[3]  = {1,1,1} ;
      cgsize_t fvm_stride[3] = {3,1,1};
      cgsize_t fvm_end[3] = {1,1,1};
      cgsize_t fvm_dims[3] = {1,1,1};

      fvm_stride[0] = ( cgsize_t ) phys_dim;

      // Sioux ...
      fvm_dims[0] = fm_dims[0]*fvm_stride[0]; //*phys_dim to support 2D aliasing
      fvm_dims[1] = fm_dims[1];
      fvm_dims[2] = fm_dims[2];
      fvm_end[0] = fm_end[0]*fvm_stride[0]; //*phys_dim for aliasing
      fvm_end[1] = fm_end[1];
      fvm_end[2] = fm_end[2];

      //==========================================================
      //
      // Count number of vars and vectors
      // Assign vars and vectors to a vtkvars array
      std::vector<vtkDataArray *> vtkvars ( nfields );

      for ( int f = 0; f < nfields; f++ )
        {
        vtkvars[f] = 0 ;

        if ( vars[f].isComponent == false )
          {
          if (IsVarEnabled( varcentering, vars[f].name ) == false){
            continue;
            }

          switch ( vars[f].dt )
            {
            // Autres cas a faire
            case CGNS_ENUMV(Integer):
              vtkvars[f] = vtkIntArray::New();
              break;
            case CGNS_ENUMV(RealSingle):
              vtkvars[f] = vtkFloatArray::New();
              break;
            case CGNS_ENUMV(RealDouble):
              vtkvars[f] = vtkDoubleArray::New();
              break;
            case CGNS_ENUMV(Character):
              vtkvars[f] = vtkCharArray::New();
              break;
            }
          vtkvars[f]->SetName ( vars[f].name );
          vtkvars[f]->SetNumberOfComponents ( 1 );
          vtkvars[f]->SetNumberOfTuples ( nVals );
          }

        }

      for ( std::vector<CGNSRead::CGNSVector>::iterator it = listOfVector.begin(); it != listOfVector.end(); ++it )
        {
        vtkDataArray *arr = 0;

        if (IsVarEnabled( varcentering, it->name ) == false){
          continue;
          }

        int nv = it->xyzIndex[0];
        switch ( vars[ nv ].dt )
          {
          // TODO Other cases to be done
          case CGNS_ENUMV(Integer):
            arr = vtkIntArray::New();
            break;
          case CGNS_ENUMV(RealSingle):
            arr = vtkFloatArray::New();
            break;
          case CGNS_ENUMV(RealDouble):
            arr = vtkDoubleArray::New();
            break;
          case CGNS_ENUMV(Character):
            arr = vtkCharArray::New();
            break;
          }

        arr->SetName ( it->name );
        arr->SetNumberOfComponents ( phys_dim );
        arr->SetNumberOfTuples ( nVals );

        for ( int dim=0; dim< phys_dim; dim++ )
          {
          arr->SetComponentName ( ( vtkIdType ) ( dim ) , vars[ it->xyzIndex[dim] ].name );
          vtkvars[ it->xyzIndex[dim] ] = arr ;
          }
        }

      // Load Data :
      for ( int f = 0; f < nfields; ++f )
        {
        // only read allocated fields
        if ( vtkvars[f] == 0 ){
          continue;
          }

        cgio_get_node_id ( this->cgio_desc, sol_cgio_id, vars[f].name, &var_cgio_id );

        // quick transfer of data because data types is given by cgns database
        if ( vars[f].isComponent == false )
          {
          if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                                fs_start, fs_end, fs_stride, cell_dim, fm_dims,
                                fm_start, fm_end, fm_stride,
                                ( void* ) vtkvars[f]->GetVoidPointer ( 0 ) ) != 0 )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          }
        else
          {
          // quick transfer of data because data types is given by cgns database
          if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                                fs_start, fs_end, fs_stride, cell_dim, fvm_dims,
                                fvm_start, fvm_end, fvm_stride,
                                ( void* ) vtkvars[f]->GetVoidPointer ( vars[f].xyzIndex-1 ) ) )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }

          }
        cgio_release_id ( this->cgio_desc, var_cgio_id );
        }
      cgio_release_id ( this->cgio_desc, sol_cgio_id );

      // Append data to StructuredGrid
      vtkDataSetAttributes* dsa = 0;
      if ( varcentering == CGNS_ENUMV ( Vertex ) ) //ON_NODES
        {
        dsa = sgrid->GetPointData();
        //dsa = mzone->GetBlock(sol-1)->GetPointData();
        }
      if ( varcentering == CGNS_ENUMV ( CellCenter ) ) //ON_CELL
        {
        dsa = sgrid->GetCellData();
        //dsa = mzone->GetBlock(sol-1)->GetCellData();
        }
      // SetData in vtk Structured Zone + Clean Pointers
      for ( int nv = 0; nv < nfields; nv++ )
        {
        // only transfer allocated fields
        if ( vtkvars[nv] == 0 ){
          continue;
          }

        if ( ( vars[nv].isComponent == false ) ||
             ( vars[nv].xyzIndex == 1 ) )
          {
          dsa->AddArray ( vtkvars[nv] );
          vtkvars[nv]->Delete();
          }
        vtkvars[nv] = 0 ;
        }

      //if ( nfields != 0 )
      // Handle Ref Values
      const std::map< std::string, double>& arrState = this->Internal.GetBase ( base-1 ).referenceState ;
      std::map< std::string, double>::const_iterator iteRef = arrState.begin();
      for ( iteRef = arrState.begin(); iteRef != arrState.end(); iteRef++ )
        {
        vtkDoubleArray* refValArray = vtkDoubleArray::New();
        refValArray->SetNumberOfComponents ( 1 );
        refValArray->SetName ( iteRef->first.c_str() );
        refValArray->InsertNextValue ( iteRef->second );
        sgrid->GetFieldData()->AddArray ( refValArray );
        refValArray->Delete();
        }
      //
      
      mzone->SetBlock ( ( sol-1 ), sgrid );
      sgrid->Delete();
      }
    mbase->SetBlock ( ( zone-1 ), mzone );
    }
  else
    {
    vtkStructuredGrid *sgrid   = vtkStructuredGrid::New();
    sgrid->SetExtent ( extent );
    sgrid->SetPoints ( points );
    mbase->SetBlock ( ( zone-1 ), sgrid );
    sgrid->Delete();
    }
  points->Delete();
  mzone->Delete();
  return 0;
}


//------------------------------------------------------------------------------
int vtkCGNSReader::GetUnstructuredZone ( int fn, int  base, int zone,
                                         int cell_dim, int phys_dim,
                                         cgsize_t *zsize,
                                         vtkMultiBlockDataSet *mbase )
{
  //
  bool err;
  int ngrids = 0;
  int rind[6];
  int i;

  int n;
  int ier;
  int ncoords = 0;

  // Check the number of grids stored in zone
  vtkDebugMacro( << "get number of grids in zone "
                 << zone << "\n" );

  if ( cg_ngrids ( fn, base, zone, &ngrids ) != CG_OK )
    {
    vtkErrorMacro( << "Could not get number of grids in zone "
                   << zone << "\n"
                   << cg_get_error() );
    return 1;
    }

  if ( ngrids < 1 ){
    vtkErrorMacro ( << "No Grid found in zone " << zone << "\n" );
    return 1;
    }

  CGNS_ENUMT(DataType_t) MeshType = CGNS_ENUMV(RealDouble) ;// TODO modify this line

  // Get the number of Coordinates in GridCoordinates node
  if ( cg_ncoords ( fn, base, zone, &ncoords ) != CG_OK )
    {
    vtkErrorMacro( << "\t\tCould not get the number of coords\n"
                   <<  cg_get_error() << "\n" );
    return 1;
    }

  if ( cg_goto ( fn, base, "Zone_t", zone,
                 "GridCoordinates_t", 1, "end" ) )
    {
    vtkErrorMacro( << "No GridCoordinates_t\n");
    return 1;
    }
  // check for rind
  if ( ( ier = cg_rind_read ( rind ) ) != CG_OK )
    {
    if ( ier != CG_NODE_NOT_FOUND ){
      vtkErrorMacro( << "Unexpected problem while reading rind information\n" );
      return 1;
      }
    for ( n = 0; n < 6; n++ ){
      rind[n] = 0;
      }
    }
  // source layout
  cgsize_t s_start[3]  = {1,1,1};
  cgsize_t s_stride[3] = {1,1,1};
  cgsize_t s_end[3];
  // Number of Dimension
  const int m_num_dims = 1;
  // memory destination layout
  cgsize_t m_start[3]  = {1,1,1} ;
  cgsize_t m_stride[3] = {3,1,1};
  cgsize_t m_end[3] = {1,1,1};
  cgsize_t m_dims[3] = {1,1,1};

  vtkIdType nPts = 0;

  /* get grid coordinate range */
  s_start[0] = rind[0] + 1;
  s_end[0]   = rind[0] + zsize[0];
  m_end[0]   = zsize[0];
  m_dims[0]  = zsize[0];

  // Compute number of points
  nPts = (vtkIdType) zsize[0];

  // Set up points
  vtkPoints *points  = vtkPoints::New();

  //
  // wacky hack ...
  m_end[0] *= 3; //for memory aliasing
  //
  // vtkPoints assumes float data type
  //
  if ( MeshType == CGNS_ENUMV(RealDouble) ){
    points->SetDataTypeToDouble();
    }
  //
  // Resize vtkPoints to fit data
  //
  points->SetNumberOfPoints ( nPts );

  //----------------------------------------------------------------------------
  // Get Coordinates and FlowSolution node names
  bool readGridCoordName = true;
  bool readSolutionName = true;
  CGNSRead::char_33 GridCoordName;
  CGNSRead::char_33 SolutionName;

  if ((this->Internal.GetBase(base-1).useGridPointers == true) ||
      (this->Internal.GetBase(base-1).useFlowPointers == true))
    {
    CGNSRead::char_33 nameZoneIter;
    size_t ptSize = 32*this->Internal.GetBase(base-1).steps.size()+1;
    char *pointers = new char[ptSize];

    cg_ziter_read(fn, base, zone, nameZoneIter);
    if ( cg_goto(fn,base,"Zone_t",zone,nameZoneIter,0,NULL) == CG_OK )
      {
      int arraycount,arrayno;
      cg_narrays( &arraycount );
      for( arrayno = 1; arrayno<=arraycount; arrayno++ )
        {
        CGNSRead::char_33 arrayname;
        CGNS_ENUMT(DataType_t) arraydt;
        int dimension;
        cgsize_t elements[ 12 ];
        cg_array_info(arrayno,arrayname,&arraydt,&dimension,elements);
        if( strcmp( arrayname,"GridCoordinatesPointers" ) == 0 )
          {
          cg_array_read( arrayno, pointers ) ;
          strncpy(GridCoordName, &pointers[this->ActualTimeStep*32], 32);
          GridCoordName[32]='\0';
          readGridCoordName = false;
          }

        if( strcmp( arrayname,"FlowSolutionPointers" ) == 0 )
          {
          cg_array_read( arrayno, pointers ) ;
          strncpy(SolutionName, &pointers[this->ActualTimeStep*32], 32);
          SolutionName[32]='\0';
          readSolutionName = false;
          }

        }
      }
    else
      {
      strcpy(GridCoordName, "GridCoordinates");
      strcpy(SolutionName, "FlowSolution");
      }
    delete [] pointers;
    }
  //----------------------------------------------------------------------------
  if ( readGridCoordName )
    {
    int requiredgrid = 1;
    cg_grid_read ( fn, base, zone, requiredgrid, GridCoordName );
    }

  int narrays=0;

  //std::cout << "Reading mesh node " << GridCoordName << std::endl;

  if ( cg_goto ( fn, base, "Zone_t", zone, GridCoordName, 0, "end" ) != CG_OK )
    {
    vtkErrorMacro( << cg_get_error() );
    return 1;
    }

  cg_narrays ( &narrays );
  if ( narrays < ncoords )
    {
    vtkErrorMacro( << "Not enough coordinates in node "
                   << GridCoordName << "\n" );
    return 1;
    }

  double gridid, coordid ;
  cgio_get_node_id ( this->cgio_desc, this->currentid, GridCoordName, &gridid );

  if ( MeshType == CGNS_ENUMV(RealDouble) )
    {
    double *coords = (double *) points->GetVoidPointer ( 0 ) ;
    double *currentCoord = (double *) &(coords[0]);

    CGNSRead::char_33 coordname;
    size_t len;
    CGNS_ENUMT(DataType_t) ct;

    memset ( coords, 0, 3*nPts*sizeof ( double ) );

    for ( int c = 1; c <= ncoords; ++c )
      {
      if ( cg_coord_info ( fn, base, zone, c, &ct,
                           coordname ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        break;
        }
      // Determine direction X,Y,Z
      len = strlen ( coordname ) - 1;
      switch ( coordname[len] )
        {
        case 'X':
          currentCoord = (double *) &(coords[0]) ;
          break;
        case 'Y':
          currentCoord = (double *) &(coords[1]) ;
          break;
        case 'Z':
          currentCoord = (double *) &(coords[2]) ;
          break;
        }

      cgio_get_node_id ( this->cgio_desc, gridid, coordname, &coordid );

      /* quick transfer of data if same data types */
      if ( ct == CGNS_ENUMV(RealDouble) )
        {
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, m_num_dims, m_end,
                              m_start, m_end, m_stride, ( void* ) currentCoord ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( << "cgio_read_data :" << message );
          }
        }
      else
        {
        float *array_data = 0;
        const cgsize_t m_nostride[3] = {1,1,1};

        if ( ct != CGNS_ENUMV(RealSingle) )
          {
          vtkErrorMacro( << "Invalid datatype for GridCoordinates\n" );
          break;
          }

        /* need to read into temp array to convert data */
        array_data = new float[nPts];
        if ( array_data == 0 )
          {
          vtkErrorMacro ( "Error allocating buffer array" );
          break;
          }
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, m_num_dims, m_dims,
                              m_start, m_dims, m_nostride, ( void* ) array_data ) )
          {
          delete [] array_data;
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "Buffer array cgio_read_data :" << message );
          break;
          }

        for (vtkIdType ii = 0; ii < nPts; ii++ )
          {
          currentCoord[m_stride[0]*ii] = ( double ) array_data[ii];
          }

        delete [] array_data;
        }
      }
    }
  else
    {
    float *coords = (float *) points->GetVoidPointer ( 0 ) ;
    float *currentCoord = (float *) &(coords[0]);

    CGNSRead::char_33 coordname;
    size_t len;
    CGNS_ENUMT(DataType_t) ct;

    for ( int c = 1; c <= ncoords; ++c )
      {
      if ( cg_coord_info ( fn, base, zone, c, &ct,
                           coordname ) != CG_OK )
        {
        vtkErrorMacro( << cg_get_error() );
        break;
        }
      // Determine direction X,Y,Z
      len = strlen ( coordname ) - 1;
      switch ( coordname[len] )
        {
        case 'X':
          currentCoord = (float *) &(coords[0]) ;
          break;
        case 'Y':
          currentCoord = (float *) &(coords[1]) ;
          break;
        case 'Z':
          currentCoord = (float *) &(coords[2]) ;
          break;
        }

      cgio_get_node_id ( this->cgio_desc, gridid, coordname, &coordid );

      /* quick transfer of data if same data types */
      if ( ct == CGNS_ENUMV(RealSingle) )
        {
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, m_num_dims , m_end,
                              m_start, m_end, m_stride, (void *) currentCoord ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "cgio_read_data :" << message );
          }
        }
      else
        {
        double *array_data = 0;
        const cgsize_t m_nostride[3] = {1,1,1};

        if ( ct != CGNS_ENUMV(RealDouble) )
          {
          vtkErrorMacro( << "Invalid datatype for GridCoordinates\n" );
          break;
          }

        /* need to read into temp array to convert data */
        array_data = new double[nPts];
        if ( array_data == 0 )
          {
          vtkErrorMacro ( << "Error allocating buffer array" );
          break;
          }
        if ( cgio_read_data ( this->cgio_desc, coordid,
                              s_start, s_end, s_stride, m_num_dims, m_dims,
                              m_start, m_dims, m_nostride, (void *) array_data ) )
          {
          delete [] array_data;
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( << "Buffer array cgio_read_data :" << message );
          break;
          }

        for (vtkIdType ii = 0; ii < nPts; ii++ )
          {
          currentCoord[m_stride[0]*ii] = (float) array_data[ii];
          }

        delete [] array_data;

        }
      }
    }
  cgio_release_id ( this->cgio_desc, gridid );
  // points are now loaded
  //----------------------------------------------------------------------------
  // Read the number of sections, for the zone.
  int nsections = 0;
  if ( cg_nsections ( fn, base, zone, &nsections ) != CG_OK )
    {
    vtkErrorMacro( << cg_get_error() );
    return 1;
    }

  cgsize_t is, ie, nelems, maxsize, maxelems;
  int ns, nn, ip;
  cgsize_t size;
  char name[33];

  // Find section layout
  // Section is composed of => 1 Volume + bnd surfaces
  //                        => multi-Volume + Bnd surfaces
  // determine dim to allocate for connectivity reading
  cgsize_t elementCoreSize = 0;
  vtkIdType numCoreCells = 0;
  //
  std::vector<int> coreSec;
  std::vector<int> bndSec;
  std::vector<int> sizeSec;
  std::vector<int> startSec;
  //
  numCoreCells = 0; // force initialize
  for ( int sec = 1; sec <= nsections; ++sec )
    {
    char sectionname[33];
    CGNS_ENUMT(ElementType_t) elemtype = CGNS_ENUMV(ElementTypeNull);
    cgsize_t start = 1, end = 1;
    cgsize_t elementSize = 0;
    int bound = 0, parent_flag = 0;
    if ( cg_section_read ( fn, base, zone, sec, sectionname, &elemtype,
                           &start, &end, &bound, &parent_flag ) != CG_OK )
      {
      std::cout << cg_get_error() << std::endl;
      continue;
      }
    elementSize = end-start+1; // Interior Volume + Bnd

    if ( elemtype != CGNS_ENUMV ( MIXED ) )
      {
      // all cells are of the same type.
      int numPointsPerCell = 0;
      int cell_type ;
      bool higherOrderWarning;
      bool reOrderElements;

      if ( cg_npe ( elemtype, &numPointsPerCell ) || numPointsPerCell == 0 )
        {
        vtkErrorMacro( << "Element type error\n" );
        }
      // make a clean GetVTKElemType
      // see Xdmf
      cell_type = CGNSRead::GetVTKElemType ( elemtype, higherOrderWarning,
                                             reOrderElements );

      cgsize_t eDataSize = 0;
      cgsize_t EltsEnd = elementSize + start -1;
      if ( cg_ElementPartialSize ( fn, base, zone, sec,
                                   start, EltsEnd, &eDataSize ) != CG_OK )
        {
        vtkErrorMacro( << "Could not determine ElementDataSize\n" );
        continue;
        }

      if ( eDataSize != numPointsPerCell*elementSize )
        {
        vtkErrorMacro( << "FATAL wrong elements dimensions\n");
        }
      if ( start > zsize[1] )
        {

        vtkDebugMacro( << "@@1: Boundary Section not accounted" << "\n" );
        bndSec.push_back ( sec );
        }
      else
        {
        sizeSec.push_back ( ( numPointsPerCell+1 ) *elementSize );
        startSec.push_back ( start-1 );
        elementCoreSize += ( numPointsPerCell+1 ) *elementSize ;
        numCoreCells += elementSize;
        coreSec.push_back ( sec );
        }
      }
    else if ( elemtype == CGNS_ENUMV ( MIXED ) )
      {
      if ( start > zsize[1] )
        {
        vtkDebugMacro( << "@@ Boundary Section not accounted" << "\n");
        bndSec.push_back ( sec );
        }
      else
        {
        cgsize_t eDataSize = 0;
        if ( cg_ElementPartialSize ( fn, base, zone, sec, start, end,
                                     &eDataSize ) != CG_OK )
          {
          vtkErrorMacro( << "Could not determine ElementDataSize\n");
          continue;
          }
        sizeSec.push_back ( eDataSize );  //+ elementSize
        startSec.push_back ( start-1 );
        elementCoreSize += ( eDataSize ) ; // + elementSize
        numCoreCells += elementSize;
        coreSec.push_back ( sec );
        }
      }
    }
  vtkIdType* startArraySec =  new vtkIdType[coreSec.size()];
  for ( int sec = 0 ; sec<coreSec.size(); sec++ )
    {
    int curstart = startSec[sec];
    vtkIdType curarraystart = 0;
    for ( int lse = 0 ; lse<coreSec.size(); lse++ )
      {
      if ( startSec[lse] < curstart )
        {
        curarraystart += sizeSec[lse] ;
        }
      }
    startArraySec[sec] = curarraystart;
    }

  // Create Cell Array
  vtkCellArray* cells = vtkCellArray::New();
  // Modification for memory reliability
  vtkIdTypeArray *cellLocations = vtkIdTypeArray::New();
  cellLocations->SetNumberOfValues ( elementCoreSize );
  vtkIdType* elements = cellLocations->GetPointer ( 0 );
  //vtkIdType* elements = new vtkIdType[elementCoreSize];

  if ( elements == 0 )
    {
    vtkErrorMacro( << "Could not allocate memory for connectivity\n" );
    return 1;
    }

  int *cells_types = new int[numCoreCells];
  if ( cells_types == 0 )
    {
    vtkErrorMacro( << "Could not allocate memory for connectivity\n" );
    return 1;
    }

  // Iterate over core sections.
  for ( std::vector<int>::iterator it = coreSec.begin() ; it != coreSec.end(); ++it )
    {
    int sec = *it;
    char sectionname[33];
    CGNS_ENUMT(ElementType_t) elemtype = CGNS_ENUMV(ElementTypeNull);
    cgsize_t start = 1, end = 1;
    cgsize_t elementSize = 0;
    int bound = 0, parent_flag = 0;
    if ( cg_section_read ( fn, base, zone, sec, sectionname, &elemtype,
                           &start, &end, &bound, &parent_flag ) != CG_OK )
      {
      std::cout << cg_get_error() << std::endl;
      continue;
      }
    elementSize = end-start+1; // Interior Volume + Bnd
    if ( start > zsize[1] )
      {
      vtkErrorMacro( << "ERROR:: Boundary Section " << end );
      }

    double section_cgio_id;
    double elem_cgio_connect;

    cgio_get_node_id ( this->cgio_desc, this->currentid, sectionname, &section_cgio_id );
    
    if ( elemtype != CGNS_ENUMV ( MIXED ) )
      {
      // all cells are of the same type.
      int numPointsPerCell = 0;
      int cell_type ;
      bool higherOrderWarning;
      bool reOrderElements;
      //
      if ( cg_npe ( elemtype, &numPointsPerCell ) || numPointsPerCell == 0 )
        {
        vtkErrorMacro(<< "Invalid numPointsPerCell\n");
        }
      // TODO : create cleaner GetVTKElemType, see Xdmf reader
      cell_type = CGNSRead::GetVTKElemType ( elemtype, higherOrderWarning,
                                             reOrderElements );
      
      //
      for ( vtkIdType i=start-1; i<end; i++ )
        {
        cells_types[i]= cell_type;
        }
      //
      cgsize_t eDataSize = 0;
      cgsize_t EltsEnd = elementSize + start -1;
      if ( cg_ElementPartialSize ( fn, base, zone, sec,
                                   start, EltsEnd, &eDataSize ) != CG_OK )
        {
        vtkErrorMacro(<< "Could not determine ElementDataSize\n");
        continue;
        }

      vtkDebugMacro( << "Element data size for sec " << sec <<
                     " is: " << eDataSize << "\n" );

      // Test at compilation time with static assert ... to be done
      // In case  cgsize_t < vtkIdType one could try to start from the array end
      if ( sizeof ( cgsize_t ) > sizeof ( vtkIdType ) )
        {
        vtkErrorMacro(<< " Impossible to load data with sizeof cgsize_t bigger than sizeof vtkIdType\n");
        return 1;
        }

      // TODO Warning at compilation time ??
      if ( sizeof ( cgsize_t ) != sizeof ( vtkIdType ) )
        {
        vtkWarningMacro(<< "Warning cgsize_t do not have same size as vtkIdType\n");
        vtkWarningMacro(<< "sizeof vtkIdType = " << sizeof ( vtkIdType ) << "\n");
        vtkWarningMacro(<< "sizeof cgsize_t = " << sizeof ( cgsize_t ) << "\n") ;
        }

      if ( eDataSize != numPointsPerCell*elementSize )
        {
        vtkErrorMacro(<< "FATAL wrong elements dimensions\n");
        }

      // strided_read not fully efficient !!
      // because it cannot skip points ... for the time being
      // pointer on start !!
      vtkIdType* locelements = & ( elements[ startArraySec[sec-1] ] );

      if ( cg_goto ( fn, base, "Zone_t", zone, "Elements_t" , sec ,"end" ) != CG_OK )
        {
        vtkErrorMacro(<< "tonk" << cg_get_error() << "\n");
        }
      // WARNING !!!
      int narrays=0;
      cg_narrays ( &narrays );

      cgsize_t m_dim[2];

      cgsize_t npe = numPointsPerCell;

      s_start[0]  = 1 ;// + npe* ( wantedstartperprocess-start ) ; startoffset
      s_start[1]  = 1;

      s_end[0] = ( EltsEnd-start+ 1 ) *npe;
      s_stride[0] = 1;

      m_start[0]  = 2;
      m_start[1]  = 1;
      m_end[0]    = npe+1;
      m_end[1]    = EltsEnd-start+1;
      m_stride[0] = 1;
      m_stride[1] = 1;
      m_dim[0]    = npe+1;
      m_dim[1]    = EltsEnd-start+1;

      const char *connectpath = "ElementConnectivity";
      char data_type[3];
      size_t sizeof_cnt;

      memset (locelements, 1, sizeof(vtkIdType)*(npe+1)*(EltsEnd-start+1));
      cgio_get_node_id ( this->cgio_desc, section_cgio_id, connectpath, &elem_cgio_connect );
      cgio_get_data_type ( this->cgio_desc, elem_cgio_connect, data_type );

      if ( strcmp ( data_type, "I4" ) == 0 )
        {
        sizeof_cnt = sizeof ( int );
        }
      else if ( strcmp ( data_type, "I8" ) == 0 )
        {
        sizeof_cnt = sizeof ( cglong_t );
        }
      else
        {
        vtkErrorMacro ( "ElementConnectivity data_type unknown" );
        }

      if ( sizeof_cnt == sizeof ( vtkIdType ) )
        {
        if ( cgio_read_data ( this->cgio_desc, elem_cgio_connect,
                              s_start, s_end, s_stride, 2, m_dim,
                              m_start, m_end, m_stride, ( void* ) locelements ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "cgio_read_data :" << message );
          }

        }
      else
        {
        // need to read into temp array to convert data
        cgsize_t nn = ( m_dim[0]*m_dim[1] );
        if (sizeof_cnt == sizeof (int))
          {
          int *data = new int[nn];
          if ( data == 0 )
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if ( cgio_read_data ( this->cgio_desc, elem_cgio_connect,
                                s_start, s_end, s_stride, 2, m_dim,
                                m_start, m_end, m_stride,
                                (void *)data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < nn; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        else if (sizeof_cnt == sizeof(cglong_t))
          {
          cglong_t* data = new cglong_t[nn];
          if (data == 0)
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if ( cgio_read_data ( this->cgio_desc, elem_cgio_connect,
                                s_start, s_end, s_stride, 2, m_dim,
                                m_start, m_end, m_stride,
                                (void *)data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < nn; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        }

      cgio_release_id ( this->cgio_desc, elem_cgio_connect );

      //-----------------------------------------
      // Add numptspercell and do -1 on indexes
      for ( vtkIdType icell = 0; icell < elementSize; ++icell )
        {
        locelements[icell* ( numPointsPerCell+1 )] = ( vtkIdType ) numPointsPerCell;
        for ( vtkIdType ip = 0; ip < numPointsPerCell; ip++ )
          {
          locelements[ip+icell* ( numPointsPerCell+1 ) +1] = locelements[ip+icell* ( numPointsPerCell+1 ) +1] - 1;
          }
        }

      if (reOrderElements == true)
        {
        CGNSRead::CGNS2VTKorderMonoElem(elementSize, cell_type, locelements);
        }
      }
    else if ( elemtype == CGNS_ENUMV ( MIXED ) )
      {
      //
      // all cells are of the same type.
      int numPointsPerCell = 0;
      int cell_type ;
      bool higherOrderWarning;
      bool reOrderElements;
      //
      // strided_read not fully efficient !!
      // because it cannot skip points ...
      //pointer on start !!     
      vtkIdType* locelements = & ( elements[ startArraySec[sec-1] ] );
      if ( cg_goto ( fn, base, "Zone_t", zone, "Elements_t" , sec ,"end" ) != CG_OK )
        {
        vtkErrorMacro( "message: "<< cg_get_error() << "\n");
        }
      // WARNING !!!
      int narrays=0;
      cg_narrays ( &narrays );

      //ElementConnectivity node !! --> which position ??
      cgsize_t EltsEnd = elementSize + start -1;

      cgsize_t eDataSize = 0;
      if ( cg_ElementDataSize ( fn, base, zone, sec,&eDataSize ) != CG_OK )
        {
        vtkErrorMacro("Could not determine ElementDataSize\n");
        continue;
        }

      const char *connectpath = "ElementConnectivity";
      char data_type[3];
      size_t sizeof_cnt;

      cgio_get_node_id ( this->cgio_desc, section_cgio_id, connectpath, &elem_cgio_connect );
      cgio_get_data_type ( this->cgio_desc, elem_cgio_connect, data_type );

      if ( strcmp ( data_type, "I4" ) == 0 )
        {
        sizeof_cnt = sizeof ( int );
        }
      else if ( strcmp ( data_type, "I8" ) == 0 )
        {
        sizeof_cnt = sizeof ( cglong_t );
        }
      else
        {
        vtkErrorMacro ( "ElementConnectivity data_type unknown\n" );
        }

      if ( sizeof_cnt == sizeof ( vtkIdType ) )
        {

        if ( cgio_read_all_data ( this->cgio_desc, elem_cgio_connect, ( void* ) locelements ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "cgio_read_data :" << message );
          }
        }
      else
        {
        // need to read into temp array to convert data
        if (sizeof_cnt == sizeof (int))
          {
          int *data = new int[eDataSize];
          if ( data == 0 )
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if (  cgio_read_all_data ( this->cgio_desc, elem_cgio_connect,
                                     ( void* ) data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < eDataSize; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        else if (sizeof_cnt == sizeof(cglong_t))
          {
          cglong_t* data = new cglong_t[eDataSize];
          if (data == 0)
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if ( cgio_read_all_data ( this->cgio_desc, elem_cgio_connect,
                                    ( void* ) data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < eDataSize; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        }
      cgio_release_id ( this->cgio_desc, elem_cgio_connect );


      vtkIdType pos = 0;
      reOrderElements = false;
      for ( vtkIdType icell = 0, i=start-1; icell < elementSize; ++icell, ++i )
        {
        bool orderFlag;
        elemtype = ( CGNS_ENUMT(ElementType_t) ) locelements[pos];
        cg_npe ( elemtype, &numPointsPerCell );
        cell_type = CGNSRead::GetVTKElemType ( elemtype ,higherOrderWarning,
                                               orderFlag );
        reOrderElements = reOrderElements|orderFlag;
        cells_types[i]= cell_type;
        locelements[pos] = ( vtkIdType ) numPointsPerCell;
        pos++;
        for ( vtkIdType ip = 0; ip < numPointsPerCell; ip++ )
          {
          locelements[ip+pos] = locelements[ip+pos] - 1;
          }
        pos += numPointsPerCell;
        }

        if (reOrderElements == true)
          {
          CGNSRead::CGNS2VTKorder(elementSize, &cells_types[start-1], locelements);
          }
      }
    
    cgio_release_id ( this->cgio_desc,section_cgio_id );
    }
  delete[]  startArraySec;

  cells->SetCells ( numCoreCells , cellLocations );
  cellLocations->Delete();
  //
  bool requiredPatch = ( this->LoadBndPatch != 0 );
  // SetUp zone Blocks
  vtkMultiBlockDataSet* mzone = vtkMultiBlockDataSet::New();
  if ( bndSec.size() > 0 && requiredPatch )
    {
    mzone->SetNumberOfBlocks ( 2 );
    }
  else
    {
    mzone->SetNumberOfBlocks ( 1 );
    }
  mzone->GetMetaData ( ( unsigned int ) 0 )->Set(vtkCompositeDataSet::NAME(),"Internal");

  // Set up ugrid
  // Create an unstructured grid to contain the points.
  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
  ugrid->SetPoints ( points );

  ugrid->SetCells ( cells_types, cells );

  cells->Delete();
  delete [] cells_types;

  //----------------------------------------------------------------------------
  // Handle solutions
  //----------------------------------------------------------------------------
  int nsols=0;

  ier = cg_nsols ( fn, base, zone, &nsols );

  if ( nsols > 0 )
    {
    int sol=1;
    CGNS_ENUMT(GridLocation_t) varcentering=CGNS_ENUMV(Vertex);
    double sol_cgio_id, var_cgio_id ;

    if ( readSolutionName == true )
      {
      int requiredsol = 1;
      cg_sol_info ( fn, base, zone, requiredsol, SolutionName, &varcentering );
      sol = requiredsol;
      }

    if ( cg_goto ( fn, base, "Zone_t", zone, SolutionName, 0, "end" ) != CG_OK )
      {
      vtkErrorMacro("message: " << cg_get_error());
      return 1;
      }

    if ( readSolutionName != true )
      {
      cg_gridlocation_read( &varcentering );
      int fileNumber;
      int curBase;
      int curDepth = 2;

      CGNSRead::char_33 curLabel[2];
      int numList[2];
      char ** llabel = new char*[2];
      llabel[0] = (char *) &curLabel[0];
      llabel[1] = (char *) &curLabel[1];

      cg_where(&fileNumber, &curBase, &curDepth, llabel, numList);
      sol = numList[1];
      delete [] llabel;
      }

    cgio_get_node_id ( this->cgio_desc, this->currentid, SolutionName, &sol_cgio_id );

    {
    bool fieldNotFound = true;

    int nfields ;
    ier = cg_nfields ( fn, base, zone, sol, &nfields );

    std::vector< CGNSRead::CGNSVariable > vars ( nfields ) ;
    std::vector< CGNSRead::CGNSVector > listOfVector ;
    // Read variable names
    for ( int f = 0; f < nfields; ++f )
      {
      cg_field_info ( fn, base, zone, sol, f+1, & ( vars[f].dt ), vars[f].name );

      vars[f].isComponent = false;
      vars[f].xyzIndex = 0;
      }
    // Create vector name from available variable
    // when VarX, VarY, VarZ is detected
    CGNSRead::fillVectorsFromVars(vars, listOfVector, phys_dim);

    // check for rind
    if ( ( ier = cg_rind_read ( rind ) ) != CG_OK )
      {
      if ( ier != CG_NODE_NOT_FOUND ){
        vtkErrorMacro( "Unexpected error while reading rind information in solution node " << SolutionName );
        return 1;
        }
      for ( n = 0; n < 6; n++ ){
        rind[n] = 0;
        }
      }

    // source
    cgsize_t fs_start[3]  = {1,1,1};
    cgsize_t fs_stride[3] = {1,1,1};
    cgsize_t fs_end[3];
    //
    const int m_num_dims = 1;
    // destination
    cgsize_t fm_start[3]  = {1,1,1} ;
    cgsize_t fm_stride[3] = {1,1,1};
    cgsize_t fm_end[3] = {1,1,1};
    cgsize_t fm_dims[3] = {1,1,1};

    vtkIdType nVals = 0;

    // get solution data range
    int nsc = varcentering == CGNS_ENUMV ( Vertex ) ? 0 : 1;
    fs_start[0] = rind[0] + 1;
    fs_end[0]   = rind[0] + zsize[nsc];
    fm_end[0]   = zsize[nsc];
    fm_dims[0]  = zsize[nsc];
    // compute number of elements
    nVals = (vtkIdType) fm_end[0] ;
    //==========================================================
    // VECTORS aliasing ...
    // destination
    cgsize_t fvm_start[3]  = {1,1,1} ;
    cgsize_t fvm_stride[3] = {3,1,1}; //
    cgsize_t fvm_end[3] = {1,1,1};
    cgsize_t fvm_dims[3] = {1,1,1};

    fvm_stride[0] = ( cgsize_t ) phys_dim;

    // Sioux ...
    fvm_dims[0] = fm_dims[0]*fvm_stride[0]; //*phys_dim to support 2D aliasing
    fvm_dims[1] = fm_dims[1];
    fvm_dims[2] = fm_dims[2];
    fvm_end[0] = fm_end[0]*fvm_stride[0]; //*phys_dim for aliasing
    fvm_end[1] = fm_end[1];
    fvm_end[2] = fm_end[2];

    //==========================================================

    //
    // Count number of vars and vectors
    // Assign vars and vectors to a vtkvars array
    std::vector<vtkDataArray *> vtkvars ( nfields );

    for ( int f = 0; f < nfields; f++ )
      {
      vtkvars[f] = 0 ;

      if ( vars[f].isComponent == false )
        {
        if (IsVarEnabled( varcentering, vars[f].name ) == false)
          {
          continue;
          }

        switch ( vars[f].dt )
          {
          // Autres cas a faire
          case CGNS_ENUMV(Integer):
            vtkvars[f] = vtkIntArray::New();
            break;
          case CGNS_ENUMV(RealSingle):
            vtkvars[f] = vtkFloatArray::New();
            break;
          case CGNS_ENUMV(RealDouble):
            vtkvars[f] = vtkDoubleArray::New();
            break;
          case CGNS_ENUMV(Character):
            vtkvars[f] = vtkCharArray::New();
            break;
          }
        vtkvars[f]->SetName ( vars[f].name );
        vtkvars[f]->SetNumberOfComponents ( 1 );
        vtkvars[f]->SetNumberOfTuples ( nVals );
        }

      }

    for (std::vector<CGNSRead::CGNSVector>::iterator it = listOfVector.begin();
         it != listOfVector.end(); ++it)
      {
      vtkDataArray *arr = 0;

      if (IsVarEnabled( varcentering, it->name ) == false)
        {
        continue;
        }

      int nv = it->xyzIndex[0];
      switch ( vars[ nv ].dt )
        {
        // TODO: other cases
        case CGNS_ENUMV(Integer):
          arr = vtkIntArray::New();
          break;
        case CGNS_ENUMV(RealSingle):
          arr = vtkFloatArray::New();
          break;
        case CGNS_ENUMV(RealDouble):
          arr = vtkDoubleArray::New();
          break;
        case CGNS_ENUMV(Character):
          arr = vtkCharArray::New();
          break;
        }

      arr->SetName ( it->name );
      arr->SetNumberOfComponents ( phys_dim );
      arr->SetNumberOfTuples ( nVals );

      for ( int dim=0; dim< phys_dim; dim++ )
        {
        arr->SetComponentName ( ( vtkIdType ) ( dim ) , vars[ it->xyzIndex[dim] ].name );
        vtkvars[ it->xyzIndex[dim] ] = arr ;
        }
      }

    // Load Data :
    for ( int f = 0; f < nfields; ++f )
      {
      // only read allocated fields
      if ( vtkvars[f] == 0 ){
        continue;
        }

      cgio_get_node_id ( this->cgio_desc, sol_cgio_id, vars[f].name, &var_cgio_id );

      // quick transfer of data because data types come from cgns database
      if ( vars[f].isComponent == false )
        {
        if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                              fs_start, fs_end, fs_stride, m_num_dims, fm_dims,
                              fm_start, fm_end, fm_stride,
                              ( void* ) vtkvars[f]->GetVoidPointer ( 0 ) ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "cgio_read_data :" << message );
          }

        }
      else
        {
        // quick transfer of data because data types come from cgns database
        if ( cgio_read_data ( this->cgio_desc, var_cgio_id,
                              fs_start, fs_end, fs_stride, m_num_dims, fvm_dims,
                              fvm_start, fvm_end, fvm_stride,
                              ( void* ) vtkvars[f]->GetVoidPointer ( vars[f].xyzIndex-1 ) ) )
          {
          char message[81];
          cgio_error_message ( message );
          vtkErrorMacro ( "cgio_read_data :" << message );
          }

        }
      cgio_release_id ( this->cgio_desc, var_cgio_id );
      }
    cgio_release_id ( this->cgio_desc, sol_cgio_id );

    // Append data to UnstructuredGrid
    vtkDataSetAttributes* dsa;
    if ( varcentering == CGNS_ENUMV ( Vertex ) ) //ON_NODES
      {
      dsa = ugrid->GetPointData();
      }
    if ( varcentering == CGNS_ENUMV ( CellCenter ) ) //ON_CELL
      {
      dsa = ugrid->GetCellData();
      }
    // SetData in vtk Unstructured Zone + Clean Pointers
    for ( int nv = 0; nv < nfields; nv++ )
      {

      // only transfer allocated fields
      if ( vtkvars[nv] == 0 ){
        continue;
        }

      if ( ( vars[nv].isComponent == false ) ||
           ( vars[nv].xyzIndex == 1 ) )
        {
        dsa->AddArray ( vtkvars[nv] );
        vtkvars[nv]->Delete();
        }
      vtkvars[nv] = 0 ;
      }
    }
    }
  
    // Handle Ref Values
    const std::map< std::string, double>& arrState = this->Internal.GetBase(base-1).referenceState ;
    std::map< std::string, double>::const_iterator iteRef = arrState.begin();
    for ( iteRef = arrState.begin(); iteRef != arrState.end(); iteRef++)
    {
      vtkDoubleArray* refValArray = vtkDoubleArray::New();
      refValArray->SetNumberOfComponents(1);
      refValArray->SetName( iteRef->first.c_str());
      refValArray->InsertNextValue( iteRef->second );      
      ugrid->GetFieldData()->AddArray( refValArray );
      refValArray->Delete();
    }
  //--------------------------------------------------
  // Read patch boundary Sections
  //--------------------------------------------------
  // Iterate over bnd sections.
  vtkIntArray *ugrid_id_arr = vtkIntArray::New();
  ugrid_id_arr->SetNumberOfTuples ( 1 );
  ugrid_id_arr->SetValue ( 0, 0 );
  ugrid_id_arr->SetName ( "ispatch" );
  ugrid->GetFieldData()->AddArray ( ugrid_id_arr );
  ugrid_id_arr->Delete();
  
  if ( bndSec.size() > 0  && requiredPatch )
    {
    // mzone Set Blocks
    mzone->SetBlock ( 0, ugrid );
    ugrid->Delete();
    vtkMultiBlockDataSet* mpatch = vtkMultiBlockDataSet::New();
    mpatch->SetNumberOfBlocks ( bndSec.size() );

    int bndnum = 0;
    for ( std::vector<int>::iterator it = bndSec.begin() ; it != bndSec.end(); ++it )
      {
      int sec = *it;
      char sectionname[33];
      CGNS_ENUMT ( ElementType_t ) elemtype = CGNS_ENUMV(ElementTypeNull);
      cgsize_t start = 1, end = 1;
      cgsize_t elementSize = 0;
      int bound = 0, parent_flag = 0;

      if ( cg_section_read ( fn, base, zone, sec, sectionname, &elemtype,
                             &start, &end, &bound, &parent_flag ) != CG_OK )
        {
        vtkErrorMacro( "message" << cg_get_error() << "\n");
        continue;
        }
      mpatch->GetMetaData ( ( unsigned int ) bndnum )->Set ( vtkCompositeDataSet::NAME(), sectionname );
      elementSize = end-start+1; // Bnd Volume + Bnd
      if ( start < zsize[1] )
        {
        vtkErrorMacro( "ERROR:: Internal Section " << end );
        }

      int *bnd_cells_types = new int[elementSize];
      if ( bnd_cells_types == 0 )
        {
        vtkErrorMacro( "Could not allocate memory for connectivity\n");
        return 1;
        }

      cgsize_t eDataSize = 0;
      cgsize_t EltsEnd = elementSize + start -1;
      if ( cg_ElementPartialSize ( fn, base, zone, sec,
                                   start, EltsEnd, &eDataSize ) != CG_OK )
        {
        vtkErrorMacro( "Could not determine ElementDataSize\n");
        continue;
        }
      vtkDebugMacro( "Element data size for sec " << sec << " is: " << eDataSize << "\n");
      // Allocate
      cgsize_t elementBndSize = 0;
      elementBndSize = eDataSize;
      //
      // Modification ...
      vtkIdTypeArray* IdBndArray_ptr = vtkIdTypeArray::New();
      vtkIdType* bnd_elements = NULL;

      double section_cgio_id;
      double elem_cgio_connect;

      cgio_get_node_id ( this->cgio_desc, this->currentid, sectionname, &section_cgio_id );


      if ( elemtype != CGNS_ENUMV ( MIXED ) )
        {
        // all cells are of the same type.
        int numPointsPerCell = 0;
        int cell_type ;
        bool higherOrderWarning;
        bool reOrderElements;
        //
        //
        if ( cg_npe ( elemtype, &numPointsPerCell ) || numPointsPerCell == 0 )
          {
          vtkErrorMacro(<< "Bug\n");
          }
        // A mettre dans une fonction GetVTKElemType propre
        // voir Xdmf
        cell_type = CGNSRead::GetVTKElemType ( elemtype ,higherOrderWarning,
                                               reOrderElements );
        //
        //
        for ( vtkIdType i=0; i<elementSize; i++ )
          {
          bnd_cells_types[i]= cell_type;
          }
        //
        elementBndSize = ( numPointsPerCell+1 ) *elementSize;
        IdBndArray_ptr->SetNumberOfValues ( elementBndSize );
        bnd_elements = IdBndArray_ptr->GetPointer ( 0 );
        //bnd_elements = new vtkIdType[elementBndSize];
        if ( bnd_elements == 0 )
          {
          std::cout << "Could not allocate memory for bnd connectivity\n";
          return 1;
          }
        // In particular case cgsize_t < vtkIdType it may be possible to start
        // at the array ending to relocate cgsize_t
        if ( sizeof ( cgsize_t ) > sizeof ( vtkIdType ) )
          {
          std::cout << " Impossible to load data with sizeof cgsize_t bigger than sizeof vtkIdType\n";
          return 1;
          }
        //
        if ( sizeof ( cgsize_t ) != sizeof ( vtkIdType ) )
          {
          std::cout << "Warning cgsize_t do not have same size as vtkIdType\n" ;
          std::cout << "sizeof vtkIdType = " << sizeof ( vtkIdType ) << "\n" ;
          std::cout << "sizeof cgsize_t = " << sizeof ( cgsize_t ) << "\n" ;
          }

        if ( eDataSize != numPointsPerCell*elementSize )
          {
          std::cout << "FATAL wrong elements dimensions\n";
          }

        // strided_read not fully efficient !!
        // because it cannot skip points ...
        //pointer on start !!
        vtkIdType* locelements = & ( bnd_elements[0] );
        if ( cg_goto ( fn, base, "Zone_t", zone, "Elements_t" , sec ,"end" ) != CG_OK )
          {
          std::cout << "tonk" << cg_get_error() << "\n";
          }
        // WARNING !!!
        int narrays=0;
        cg_narrays ( &narrays );

        cgsize_t m_dim[2];
        cgsize_t npe = numPointsPerCell;

        s_start[0]  = 1 ;// + npe* ( wantedstartperprocess-start ) ; startoffset
        s_start[1]  = 1;

        s_end[0] = ( EltsEnd-start+ 1 ) *npe;
        s_stride[0] = 1;

        m_start[0]  = 2;
        m_start[1]  = 1;
        m_end[0]    = npe+1;
        m_end[1]    = EltsEnd-start+1;
        m_stride[0] = 1;
        m_stride[1] = 1;
        m_dim[0]    = npe+1;
        m_dim[1]    = EltsEnd-start+1;

        const char *connectpath = "ElementConnectivity";
        char data_type[3];
        size_t sizeof_cnt;

        cgio_get_node_id ( this->cgio_desc, section_cgio_id, connectpath, &elem_cgio_connect );
        cgio_get_data_type ( this->cgio_desc, elem_cgio_connect, data_type );

        if ( strcmp ( data_type, "I4" ) == 0 )
          {
          sizeof_cnt = sizeof ( int );
          }
        else if ( strcmp ( data_type, "I8" ) == 0 )
          {
          sizeof_cnt = sizeof ( cglong_t );
          }
        else
          {
          vtkErrorMacro ( "ElementConnectivity data_type unknown" );
          }

        if ( sizeof_cnt == sizeof ( vtkIdType ) )
          {

          if ( cgio_read_data ( this->cgio_desc, elem_cgio_connect,
                                s_start, s_end, s_stride, 2, m_dim,
                                m_start, m_end, m_stride, ( void* ) locelements ) )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }

          }
        else
          {
        // need to read into temp array to convert data
        cgsize_t nn = ( m_dim[0]*m_dim[1] );
        if (sizeof_cnt == sizeof (int))
          {
          int *data = new int[nn];
          if ( data == 0 )
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if ( cgio_read_data ( this->cgio_desc, elem_cgio_connect,
                                s_start, s_end, s_stride, 2, m_dim,
                                m_start, m_end, m_stride,
                                (void *)data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < nn; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        else if (sizeof_cnt == sizeof(cglong_t))
          {
          cglong_t* data = new cglong_t[nn];
          if (data == 0)
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if ( cgio_read_data ( this->cgio_desc, elem_cgio_connect,
                                s_start, s_end, s_stride, 2, m_dim,
                                m_start, m_end, m_stride,
                                (void *)data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < nn; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        }
        
        cgio_release_id ( this->cgio_desc, elem_cgio_connect );

        // Add numptspercell and do -1 on indexes
        for ( vtkIdType icell = 0; icell < elementSize; ++icell )
          {
          locelements[icell* ( numPointsPerCell+1 )] = ( vtkIdType ) numPointsPerCell;
          for ( vtkIdType ip = 0; ip < numPointsPerCell; ip++ )
            {
            locelements[ip+icell* ( numPointsPerCell+1 ) +1] = locelements[ip+icell* ( numPointsPerCell+1 ) +1] - 1;
            }
          }

        }
      else if ( elemtype == CGNS_ENUMV ( MIXED ) )
        {
        //
        // all cells are of the same type.
        int numPointsPerCell = 0;
        int cell_type ;
        bool higherOrderWarning;
        bool reOrderElements;
        //
        // strided_read not fully efficient !!
        // because it cannot skip points ...
        //pointer on start !!

        //bnd_elements = new vtkIdType[elementBndSize];
        IdBndArray_ptr->SetNumberOfValues ( elementBndSize );
        bnd_elements = IdBndArray_ptr->GetPointer ( 0 );

        if ( bnd_elements == 0 )
          {
          std::cout << "Could not allocate memory for bnd connectivity\n";
          return 1;
          }
        //
        vtkIdType* locelements = & ( bnd_elements[0] );
        if ( cg_goto ( fn, base, "Zone_t", zone, "Elements_t" , sec ,"end" ) != CG_OK )
          {
          std::cout << "tonk" << cg_get_error() << "\n";
          }
        // WARNING !!!
        int narrays=0;
        cg_narrays ( &narrays );
        //ElementConnectivity node !! --> which position ??

        const char *connectpath = "ElementConnectivity";
        char data_type[3];
        size_t sizeof_cnt;

        cgio_get_node_id ( this->cgio_desc, section_cgio_id, connectpath, &elem_cgio_connect );
        cgio_get_data_type ( this->cgio_desc, elem_cgio_connect, data_type );

        if ( strcmp ( data_type, "I4" ) == 0 )
          {
          sizeof_cnt = sizeof ( int );
          }
        else if ( strcmp ( data_type, "I8" ) == 0 )
          {
          sizeof_cnt = sizeof ( cglong_t );
          }
        else
          {
          vtkErrorMacro ( "ElementConnectivity data_type unknown" );
          }

        if ( sizeof_cnt == sizeof ( vtkIdType ) )
          {

          if ( cgio_read_all_data ( this->cgio_desc, elem_cgio_connect, ( void* ) locelements ) )
            {
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          }
        else
          {
        // need to read into temp array to convert data
        if (sizeof_cnt == sizeof (int))
          {
          int *data = new int[eDataSize];
          if ( data == 0 )
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if (  cgio_read_all_data ( this->cgio_desc, elem_cgio_connect,
                                     ( void* ) data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < eDataSize; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        else if (sizeof_cnt == sizeof(cglong_t))
          {
          cglong_t* data = new cglong_t[eDataSize];
          if (data == 0)
            {
            vtkErrorMacro ( "malloc failed for tempory array" );
            }
          if ( cgio_read_all_data ( this->cgio_desc, elem_cgio_connect,
                                    ( void* ) data ) != CG_OK )
            {
            delete[] data;
            char message[81];
            cgio_error_message ( message );
            vtkErrorMacro ( "cgio_read_data :" << message );
            }
          for ( cgsize_t n = 0; n < eDataSize; n++ )
            {
            locelements[n] = static_cast<vtkIdType>(data[n]);
            }
          delete[] data;
          }
        }
        
        cgio_release_id ( this->cgio_desc, elem_cgio_connect );


        vtkIdType pos = 0;
        for ( vtkIdType icell = 0; icell < elementSize; ++icell )
          {
          elemtype = ( CGNS_ENUMT(ElementType_t) ) locelements[pos];
          cg_npe ( elemtype, &numPointsPerCell );
          cell_type = CGNSRead::GetVTKElemType ( elemtype ,higherOrderWarning,
                                                 reOrderElements );
          bnd_cells_types[icell]= cell_type;
          locelements[pos] = ( vtkIdType ) numPointsPerCell;
          pos++;
          for ( vtkIdType ip = 0; ip < numPointsPerCell; ip++ )
            {
            locelements[ip+pos] = locelements[ip+pos] - 1;
            }
          pos += numPointsPerCell;
          }
        }

      // Create Cell Array
      vtkCellArray* bnd_cells = vtkCellArray::New();
      //vtkIdTypeArray* IdBndArray_ptr = vtkIdTypeArray::New();
      //IdBndArray_ptr->SetArray ( bnd_elements, elementBndSize , 1 );
      bnd_cells->SetCells ( elementSize , IdBndArray_ptr );
      IdBndArray_ptr->Delete();
      // Set up ugrid
      // Create an unstructured grid to contain the points.
      vtkUnstructuredGrid *bndugrid = vtkUnstructuredGrid::New();
      bndugrid->SetPoints ( points );
      bndugrid->SetCells ( bnd_cells_types, bnd_cells );
      bnd_cells->Delete();
      delete [] bnd_cells_types ;

      //
      // Add ispatch 0=false/1=true as field data
      //
      vtkIntArray *bnd_id_arr = vtkIntArray::New();
      bnd_id_arr->SetNumberOfTuples ( 1 );
      bnd_id_arr->SetValue ( 0, 1 );
      bnd_id_arr->SetName ( "ispatch" );
      bndugrid->GetFieldData()->AddArray ( bnd_id_arr );
      bnd_id_arr->Delete();

      // Handle Ref Values
      const std::map< std::string, double>& arrState = this->Internal.GetBase(base-1).referenceState ;
      std::map< std::string, double>::const_iterator iteRef = arrState.begin();
      for ( iteRef = arrState.begin(); iteRef != arrState.end(); iteRef++)
      {
      vtkDoubleArray* refValArray = vtkDoubleArray::New();
      refValArray->SetNumberOfComponents(1);
      refValArray->SetName( iteRef->first.c_str() );
      refValArray->InsertNextValue( iteRef->second );      
      bndugrid->GetFieldData()->AddArray( refValArray );
      refValArray->Delete();
      }
      
      // Copy PointData if exists
      vtkPointData* temp = ugrid->GetPointData();
      if ( temp != NULL )
        {
        int NumArray = temp->GetNumberOfArrays();
        for ( int i = 0 ; i< NumArray; ++i )
          {
          vtkDataArray* datatmp = temp->GetArray ( i );
          bndugrid->GetPointData()->AddArray ( datatmp );
          }
        }
      mpatch->SetBlock ( ( bndnum ), bndugrid );
      bndugrid->Delete();
      bndnum++;
      }
    mzone->SetBlock ( 1 , mpatch );
    mpatch->Delete();
    mzone->GetMetaData ( ( unsigned int ) 1 )->Set ( vtkCompositeDataSet::NAME(), "Patches" );
    }
  //
  points->Delete();
  vtkDebugMacro( << "Points released\n") ;
  
  if ( bndSec.size() > 0  && requiredPatch )
    {
    mbase->SetBlock ( ( zone-1 ), mzone );
    }
  else
    {
    mbase->SetBlock ( ( zone-1 ), ugrid );
    ugrid->Delete();
    }
  mzone->Delete();
  return 0;
}


//----------------------------------------------------------------------------
class WithinTolerance: public std::binary_function<double, double, bool>
{
public:
  result_type operator()(first_argument_type a, second_argument_type b) const
  {
    bool result = (fabs(a-b)<=(a*1E-6));
    return (result_type)result;
  }
};


//----------------------------------------------------------------------------
int vtkCGNSReader::RequestData ( vtkInformation *vtkNotUsed ( request ),
                                 vtkInformationVector **vtkNotUsed ( inputVector ),
                                 vtkInformationVector *outputVector )
{
  int ier;
  int nbases;
  int nzones;
  int fn;
  int nSelectedBases = 0;
  unsigned int blockIndex = 0 ;

#ifdef PARAVIEW_USE_MPI
  int processNumber;
  int numProcessors;
  int startRange, endRange;

  vtkInformation* outInfo = outputVector->GetInformationObject( 0 );
  // get the output
  vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get( vtkDataObject::DATA_OBJECT() ) );

  // The whole notion of pieces for this reader is really
  // just a division of zones between processors
  processNumber =
    outInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER() );
  numProcessors =
    outInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES() );

  int numBases = this->Internal.GetNumberOfBaseNodes();
  int numZones = 0;
  for ( int bb=0; bb< numBases; bb++)
    {
    numZones += this->Internal.GetBase(bb).nzones ;
    }

  // Divide the files evenly between processors
  int num_zones_per_process = numZones / numProcessors;

  // This if/else logic is for when you don't have a nice even division of files
  // Each process computes which sequence of files it needs to read in
  int left_over_zones = numZones - (num_zones_per_process*numProcessors);
  // base --> startZone,endZone
  typedef std::tr1::array<int, 2> duo_t;
  std::map<int, duo_t> baseToZoneRange;

  // REDO this part !!!!
  if ( processNumber < left_over_zones )
    {
    int accumulated = 0;
    startRange = (num_zones_per_process+1) * processNumber;
    endRange = startRange + (num_zones_per_process+1);
    for ( int bb=0; bb< numBases; bb++)
      {
      duo_t zoneRange = {0,0};
      startRange = startRange - accumulated;
      endRange = endRange  - accumulated;
      int startInterZone = std::max(startRange, 0);
      int endInterZone = std::min(endRange, this->Internal.GetBase(bb).nzones);

    if ((endInterZone - startInterZone) > 0)
        {
        zoneRange[0] = startInterZone;
        zoneRange[1] = endInterZone;
        }
      accumulated = this->Internal.GetBase(bb).nzones;
      baseToZoneRange[bb] = zoneRange;
      }
    }
  else
    {
    int accumulated = 0;
    startRange = num_zones_per_process * processNumber + left_over_zones;
    endRange = startRange + num_zones_per_process;
    for ( int bb=0; bb< numBases; bb++)
      {
      duo_t zoneRange = {0,0};
      startRange = startRange - accumulated;
      endRange = endRange  - accumulated;
      int startInterZone = std::max(startRange, 0);
      int endInterZone = std::min(endRange, this->Internal.GetBase(bb).nzones);
      if ((endInterZone - startInterZone) > 0)
        {
        zoneRange[0] = startInterZone;
        zoneRange[1] = endInterZone;
        }
      accumulated = this->Internal.GetBase(bb).nzones;
      baseToZoneRange[bb] = zoneRange;
      }
    }
    
    //Bnd Sections Not implemented yet for parallel
    if ( numProcessors > 1 )
      {
      this->LoadBndPatch = 0;
      }
#endif

  if (!this->Internal.Parse(this->FileName))
    {
    return 0;
    }

#ifndef PARAVIEW_USE_MPI
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject ( 0 );

  // get the output
  vtkMultiBlockDataSet *output =
      vtkMultiBlockDataSet::SafeDownCast (
        outInfo->Get ( vtkMultiBlockDataSet::DATA_OBJECT() ) );
#endif

  // this->CurrentOutput = output;
  vtkMultiBlockDataSet* rootNode = output;

  vtkDebugMacro ( << "Start Loading CGNS data" );

  this->UpdateProgress ( 0.0 );
  // Setup Global Time Information
  if ( outInfo->Has ( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() ) )
    {

    // Get the requested time step. We only support requests of a single time
    // step in this reader right now
    double requestedTimeValue =
        outInfo->Get ( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );

    // Snapping and clamping of the requested time step to be in bounds
    // can be done by FileSeriesReader class or not ?
    vtkDebugMacro ( << "RequestData: requested time value: "
                    << requestedTimeValue );

    // Clamp requestedTimeValue to available time range.
    if ( requestedTimeValue < this->Internal.GetTimes().front() )
      {
      requestedTimeValue = this->Internal.GetTimes().front() ;
      }

    if ( requestedTimeValue > this->Internal.GetTimes().back() )
      {
      requestedTimeValue = this->Internal.GetTimes().back() ;
      }

    std::vector<double>::iterator timeIte = std::find_if(
          this->Internal.GetTimes().begin(), this->Internal.GetTimes().end(),
          vtkstd::bind2nd( WithinTolerance( ), requestedTimeValue ));
    //
    if ( timeIte == this->Internal.GetTimes().end() )
      {
      return 0;
      }
    requestedTimeValue = *timeIte;
    output->GetInformation()->Set ( vtkDataObject::DATA_TIME_STEP(),
                                    requestedTimeValue );
    }

  vtkDebugMacro ( << "CGNSReader::RequestData: Reading from file <"
                  << this->FileName << ">..." );

  // Openning with cgio layer
  ier = cgio_open_file ( this->FileName, CGIO_MODE_READ, 0, & ( this->cgio_desc ) );
  if ( ier != CG_OK )
    {
    vtkErrorMacro ( "" << "Error Reading file with cgio" );
    return 0;
    }
  cgio_get_root_id ( this->cgio_desc, & ( this->rootid ) );

  // Opening with mid-level
  ier = cg_open ( this->FileName, CG_MODE_READ, &fn );
  
  if ( ier != CG_OK )
    {
    vtkErrorMacro ( "" << "Error Reading file" );
    return 0;
    }

  ier = cg_nbases ( fn, &nbases );
  if ( ier != CG_OK )
    {
    vtkErrorMacro ( "" << "Error Reading file" );
    goto errorData;
    }
  
  // if only 1 base return base else return root --> TODO
  //rootNode->SetNumberOfBlocks ( nbases );
  nSelectedBases = this->BaseSelection->GetNumberOfArraysEnabled();
  rootNode->SetNumberOfBlocks ( nSelectedBases );
  blockIndex = 0 ;
  for (int base = 1; base <= nbases ; ++base)
    {
    CGNSRead::char_33 baseName;
    int cell_dim = 0;
    int phys_dim = 0;

    ier = cg_base_read(fn, base, baseName, &cell_dim, &phys_dim);
    if ( ier  != CG_OK )
      {
      vtkErrorMacro ("" << "Error reading base number " << base);
      }

    // skip unselected base
    if (this->BaseSelection->ArrayIsEnabled(baseName) == 0)
      {
      continue;
      }

    const CGNSRead::BaseInformation & curBaseInfo = this->Internal.GetBase(base-1);
    // Get timesteps here !!
    // Operate on Global time scale :
    // clamp requestedTimeValue to available time range
    // if < timemin --> timemin
    // if > timemax --> timemax
    // Then for each base get Index for TimeStep
    // if useFlowSolution read flowSolution and take name with index
    // same for use
    // Setup Global Time Information
    this->ActualTimeStep = 0;
    bool skipBase = false;

    if ( outInfo->Has ( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() ) )
      {

      // Get the requested time step. We only support requests of a single time
      // step in this reader right now
      double requestedTimeValue =
          outInfo->Get ( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );

      vtkDebugMacro ( << "RequestData: requested time value: "
                      << requestedTimeValue );

      // Clamp requestedTimeValue to available time range.
      if ( ( requestedTimeValue < this->Internal.GetTimes().front() ) ||
           ( requestedTimeValue > this->Internal.GetTimes().back() ))
        {
        skipBase = true;
        requestedTimeValue = this->Internal.GetTimes().front();
        }

      std::vector<double>::const_iterator iter ;
      iter = std::upper_bound(curBaseInfo.times.begin(),
                              curBaseInfo.times.end(), requestedTimeValue);

      if (iter == curBaseInfo.times.begin() )
        {
        // The requested time step is before any time
        this->ActualTimeStep = 0;
        }
      else
        {
        iter--;
        this->ActualTimeStep = static_cast<int>(iter - curBaseInfo.times.begin());
        }
      }
    if ( skipBase == true )
      {
      continue;
      }
    vtkMultiBlockDataSet* mbase = vtkMultiBlockDataSet::New();
    rootNode->GetMetaData ( blockIndex )->Set ( vtkCompositeDataSet::NAME(),
                                                baseName );

    cg_nzones ( fn, base, &nzones );
    if ( nzones == 0 )
      {
      vtkWarningMacro ( "" << "No zones in base " << baseName );
      }
    else
      {
      mbase->SetNumberOfBlocks ( nzones );
      }

#ifdef PARAVIEW_USE_MPI
    int zonemin = baseToZoneRange[base-1][0]+1;
    int zonemax = baseToZoneRange[base-1][1];
    for ( int zone=zonemin; zone <= zonemax; ++zone )
      {
#else
    for ( int zone=1; zone <= nzones; ++zone )
      {
#endif
      CGNSRead::char_33 zonename;
      cgsize_t zsize[9];
      CGNS_ENUMT(ZoneType_t) zt = CGNS_ENUMV(ZoneTypeNull);
      memset ( zonename, 0, 33 );
      memset ( zsize, 0, 9*sizeof ( cgsize_t ) );

      ier = cg_zone_read ( fn, base, zone , zonename, zsize );
      if ( ier != CG_OK )
        {
        vtkErrorMacro ( << "Problem while reading zone number " << zone );
        }

      mbase->GetMetaData(zone - 1)->Set(vtkCompositeDataSet::NAME(), zonename);

      char * FamilyName = NULL;
      if ( cg_famname_read(FamilyName) == CG_OK )
        {
        vtkInformationStringKey* zonefamily = new vtkInformationStringKey("FAMILY","vtkCompositeDataSet");
        mbase->GetMetaData(zone-1)->Set(zonefamily, FamilyName);
        }

      //setUp pathTozone = "/" namebase "/"+ zonename ;
      char pathTozone[67];
      pathTozone[0] = '/';
      pathTozone[1] = '\0';
      strcat ( pathTozone,baseName );
      strcat ( pathTozone,"/" );
      strcat ( pathTozone,zonename );
      cgio_get_node_id(this->cgio_desc, this->rootid, pathTozone,
                       &(this->currentid));

      ier = cg_zone_type ( fn, base, zone, &zt );
      switch ( zt )
        {
        case CGNS_ENUMV(ZoneTypeNull):
          break;
        case CGNS_ENUMV(ZoneTypeUserDefined):
          break;
        case CGNS_ENUMV(Structured):
          {
          GetCurvilinearZone(fn, base,zone, cell_dim, phys_dim, zsize, mbase);
          break;
          }
        case CGNS_ENUMV(Unstructured):
          GetUnstructuredZone(fn, base,zone, cell_dim, phys_dim, zsize, mbase);
          break;
        }
      }
    rootNode->SetBlock ( blockIndex, mbase );
    mbase->Delete();
    blockIndex++;
    }

errorData:
  cgio_close_file ( this->cgio_desc );
  cg_close ( fn );
  fn = 0;

  this->UpdateProgress ( 1.0 );
  return 1;
}

//------------------------------------------------------------------------------
int vtkCGNSReader::RequestInformation ( vtkInformation * request,
                                        vtkInformationVector **vtkNotUsed (inputVector),
                                        vtkInformationVector *outputVector )
{

#ifdef PARAVIEW_USE_MPI
  // Setting maximum number of pieces to -1 indicates to the
  // upstream consumer that I can provide the same number of pieces
  // as there are number of processors
  // get the info object
  {
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  }

  if (this->ProcRank == 0)
    {
#endif
    if (!this->FileName )
      {
      vtkErrorMacro("File name not set");
      return 0;
      }

    // First make sure the file exists.  This prevents an empty file
    // from being created on older compilers.
    if (!vtksys::SystemTools::FileExists(this->FileName))
      {
      vtkErrorMacro("Error opening file " << this->FileName);
      return false;
      }

    vtkDebugMacro(<< "CGNSReader::RequestInformation: Parsing file "
                  << this->FileName << " for fields and time steps");

    // Parse the file...
    if (!this->Internal.Parse(this->FileName))
      {
      vtkErrorMacro("Failed to parse cgns file: " << this->FileName);
      return false;
      }
#ifdef PARAVIEW_USE_MPI
    } // End_ProcRank_0

  if (this->ProcSize>1)
    {
    this->Broadcast( this->Controller );
    }
#endif

  this->NumberOfBases = this->Internal.GetNumberOfBaseNodes() ;

  // Set up time information
  if ( this->Internal.GetTimes().size() != 0 )
    {
    std::vector<double> timeSteps(this->Internal.GetTimes().begin(),
                                  this->Internal.GetTimes().end());

    vtkInformation* outInfo = outputVector->GetInformationObject ( 0 );
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                   timeSteps.data(),
                   static_cast<int>(timeSteps.size()));
    double timeRange[2];
    timeRange[0] = timeSteps.front();
    timeRange[1] = timeSteps.back();
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                 timeRange, 2);
    }

  for (int base = 0; base < this->Internal.GetNumberOfBaseNodes() ; ++base)
    {
    const CGNSRead::BaseInformation& curBase = this->Internal.GetBase(base);
    // Fill base names
    if (base == 0 && (!this->BaseSelection->ArrayExists(curBase.name)))
      {
      this->BaseSelection->EnableArray(curBase.name);
      }
    else if ( !this->BaseSelection->ArrayExists(curBase.name))
      {
      this->BaseSelection->DisableArray(curBase.name);
      }

    // Fill Variable Vertex/Cell names ... perhaps should be improved
    CGNSRead::vtkCGNSArraySelection::const_iterator iter;
    for (iter = curBase.PointDataArraySelection.begin();
         iter != curBase.PointDataArraySelection.end(); ++iter)
      {
      if (!this->PointDataArraySelection->ArrayExists(iter->first.c_str()))
        {
        this->PointDataArraySelection->DisableArray(iter->first.c_str());
        }
      }
    for (iter = curBase.CellDataArraySelection.begin();
         iter != curBase.CellDataArraySelection.end(); ++iter)
      {
      if (!this->CellDataArraySelection->ArrayExists(iter->first.c_str()))
        {
        this->CellDataArraySelection->DisableArray(iter->first.c_str());
        }
      }

    }
  return 1;
}


//------------------------------------------------------------------------------
void vtkCGNSReader::PrintSelf ( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf ( os,indent );

  os << indent << "File Name: "
     << ( this->FileName ? this->FileName : "(none)" ) << "\n";

}

//------------------------------------------------------------------------------
int vtkCGNSReader::CanReadFile ( const char *name )
{
  // return value 0: can not read
  // return value 1: can read
  int cgioNum, ierr = 1;
  double rootid, childid;
  float FileVersion = 0.0;
  char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
  char errmsg[CGIO_MAX_ERROR_LENGTH+1];
  int ndim = 0;
  cgsize_t dim_vals[12];
  int file_type = CG_FILE_NONE;

  if (cgio_open_file( name, CG_MODE_READ, CG_FILE_NONE, &cgioNum) != CG_OK)
    {
    cgio_error_message(errmsg);
    vtkErrorMacro("vtkCGNSReader::CanReadFile : "<< errmsg);
    return 0;
    }

  cgio_get_root_id( cgioNum, &rootid);
  cgio_get_file_type ( cgioNum, &file_type );

  if ( cgio_get_node_id ( cgioNum, rootid, "CGNSLibraryVersion", &childid ) )
    {
    cgio_error_message ( errmsg );
    vtkErrorMacro ( "vtkCGNSReader::CanReadFile : "<< errmsg );
    ierr = 0;
    goto CanReadError;
    }

  if ( cgio_get_data_type ( cgioNum, childid, data_type ) )
    {
    vtkErrorMacro ( "CGNS Version data type" );
    ierr = 0;
    goto CanReadError;
    }

  if ( cgio_get_dimensions ( cgioNum, childid, &ndim, dim_vals ) )
    {
    vtkErrorMacro ( "cgio_get_dimensions" );
    ierr = 0;
    goto CanReadError;
    }

  // check data type
  if ( strcmp ( data_type,"R4" ) !=0 )
    {
    vtkErrorMacro ( "Unexpected data type for CGNS-Library-Version="
                    << data_type );
    ierr = 0;
    goto CanReadError;
    }

  // check data dim
  if ( ndim != 1 || ( dim_vals[0]!=1 ) )
    {
    vtkDebugMacro ( "Wrong data dimension for CGNS-Library-Version" );
    ierr = 0;
    goto CanReadError;
    }

  // read data
  if ( cgio_read_all_data ( cgioNum, childid, &FileVersion ) )
    {
    vtkErrorMacro ( "read CGNS version number" );
    ierr = 0;
    goto CanReadError;
    }

  // Check that the library version is at least as recent as the one used
  //   to create the file being read

  if ( ( FileVersion*1000 ) > CGNS_VERSION )
    {
    // This code allows reading version newer than the lib,
    // as long as the 1st digit of the versions are equal
    if ( ( FileVersion ) > ( CGNS_VERSION / 1000 ) )
      {
      vtkErrorMacro("The file " << name <<
                    " was written with a more recent version"
                    "of the CGNS library.  You must update your CGNS"
                    "library before trying to read this file.");
      ierr = 0;
      }
    // warn only if different in second digit
    if ( ( FileVersion*10 ) > ( CGNS_VERSION / 100 ) )
      {
      vtkWarningMacro("The file being read is more recent"
                      "than the CGNS library used" );
      }
    }
  if ( ( FileVersion*100 ) < 255 )
    {
    vtkWarningMacro ( "The file being read was written with an old version"
                      "of the CGNS library. Please update your file"
                      "to a more recent version." );
    }
  vtkDebugMacro ( "FileVersion=" <<  FileVersion << "\n" );

CanReadError:
  cgio_close_file(cgioNum);
  return ierr ? 1 : 0;
}

//------------------------------------------------------------------------------
int vtkCGNSReader::FillOutputPortInformation(int vtkNotUsed (port),
                                             vtkInformation *info )
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

//----------------------------------------------------------------------------
void vtkCGNSReader::DisableAllBases()
{
  this->BaseSelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
void vtkCGNSReader::EnableAllBases()
{
  this->BaseSelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
int vtkCGNSReader::GetNumberOfBaseArrays()
{
  return this->BaseSelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
int vtkCGNSReader::GetBaseArrayStatus(const char* name)
{
  return this->BaseSelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkCGNSReader::SetBaseArrayStatus ( const char* name, int status )
{
  if (status)
    this->BaseSelection->EnableArray(name);
  else
    this->BaseSelection->DisableArray(name);
}

//----------------------------------------------------------------------------
const char* vtkCGNSReader::GetBaseArrayName ( int index )
{
  if (index >= (int) this->NumberOfBases || index < 0)
    return NULL;
  else
    return this->BaseSelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
void vtkCGNSReader::DisableAllPointArrays()
{
  this->PointDataArraySelection->DisableAllArrays();
}


//----------------------------------------------------------------------------
void vtkCGNSReader::EnableAllPointArrays()
{
  this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
int vtkCGNSReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}


//----------------------------------------------------------------------------
const char* vtkCGNSReader::GetPointArrayName ( int index )
{
  if ( index >= ( int ) this->GetNumberOfPointArrays() || index < 0 )
    return NULL;
  else
    return this->PointDataArraySelection->GetArrayName ( index );
}

//----------------------------------------------------------------------------
int vtkCGNSReader::GetPointArrayStatus ( const char* name )
{
  return this->PointDataArraySelection->ArrayIsEnabled ( name );
}


//----------------------------------------------------------------------------
void vtkCGNSReader::SetPointArrayStatus ( const char* name, int status )
{
  if ( status )
    this->PointDataArraySelection->EnableArray ( name );
  else
    this->PointDataArraySelection->DisableArray ( name );
}

//----------------------------------------------------------------------------
void vtkCGNSReader::DisableAllCellArrays()
{
  this->CellDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
void vtkCGNSReader::EnableAllCellArrays()
{
  this->CellDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
int vtkCGNSReader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkCGNSReader::GetCellArrayName ( int index )
{
  if ( index >= ( int ) this->GetNumberOfCellArrays() || index < 0 )
    return NULL;
  else
    return this->CellDataArraySelection->GetArrayName ( index );
}

//----------------------------------------------------------------------------
int vtkCGNSReader::GetCellArrayStatus ( const char* name )
{
  return this->CellDataArraySelection->ArrayIsEnabled ( name );
}

//----------------------------------------------------------------------------
void vtkCGNSReader::SetCellArrayStatus ( const char* name, int status )
{
  if ( status )
    this->CellDataArraySelection->EnableArray ( name );
  else
    this->CellDataArraySelection->DisableArray ( name );
}

//----------------------------------------------------------------------------
void vtkCGNSReader::SelectionModifiedCallback ( vtkObject*, unsigned long,
                                                void* clientdata, void* )
{
  static_cast<vtkCGNSReader*>(clientdata)->Modified();
}

#ifdef PARAVIEW_USE_MPI
//------------------------------------------------------------------------------
void vtkCGNSReader::Broadcast( vtkMultiProcessController* ctrl )
{
  if ( ctrl )
    {
    int rank = ctrl->GetLocalProcessId();
    this->Internal.Broadcast(ctrl, rank);
    }
}
#endif
