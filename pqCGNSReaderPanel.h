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
// .NAME pqCGNSReaderPanel -- custom panel for CGNS reader
// .SECTION Description
// This panel will display properties of the "CGNS" database being reader
//  like "families" available and more
//  Right now it does nothing
////
// .SECTION Caveats
//   does not work 
//
// .SECTION Thanks
// Thanks to .

#ifndef __pqCGNSReaderPanel_h
#define __pqCGNSReaderPanel_h

// ParaView Server Manager includes
#include <vtkSMIntVectorProperty.h>
#include <vtkSMProxy.h>

// Qt Includes.
#include <QCheckBox>
#include <QtDebug>

// ParaView Includes.
#include <pqPropertyWidget.h>
#include <pqTreeWidget.h>
#include <pqTreeWidgetItemObject.h>
#include <vtkPVConfig.h>

#if defined(_WIN32) || defined(WIN32)
#  if !defined(__CYGWIN__)
#    include <winsock.h> // WSADATA, include before sys/types.h
#  endif
#endif


//#include <pqNamedObjectPanel.h>
#include <pqComponentsModule.h>


class pqOutputPort;
class pqTreeWidgetItemObject;
class QPixmap;
class QTreeWidget;
class QTreeWidgetItem;
class vtkPVArrayInformation;
class vtkSMProperty;



/// Custom panel for CGNSReader source.
/// This panel is only provided to add capabilities to hide/gray out
/// GUI elements in case they do not apply and to link the checkbox
/// ImportTracers to PointArrayStatus.
class pqCGNSReaderPanel : public pqPropertyWidget
{
  Q_OBJECT;
  typedef pqPropertyWidget Superclass;

public:
  pqCGNSReaderPanel(vtkSMProxy *smproxy, vtkSMPropertyGroup *smproperty, QWidget *parentObject=0);
 
private:
  Q_DISABLE_COPY(pqCGNSReaderPanel)

};

#endif //__pqCGNSReaderPanel_h
