set(PLUGIN_CGNS_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
set(PLUGIN_CGNS_CMAKE_PATH ${PLUGIN_CGNS_DIR}/CMake)

#include the CGNS cmake directory on the cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PLUGIN_CGNS_CMAKE_PATH})

#
find_package(CGNS REQUIRED)
IF(NOT CGNS_FOUND)
  message(FATAL_ERROR "Unable to locate CGNS library.")
ENDIF(NOT CGNS_FOUND)

#
include_directories(${CGNS_INCLUDE_DIR})

SET_SOURCE_FILES_PROPERTIES(
vtkCGNSReaderInternal
cgio_helpers

WRAP_EXCLUDE
)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
IF(PARAVIEW_BUILD_QT_GUI)
  INCLUDE(${QT_USE_FILE})
ENDIF(PARAVIEW_BUILD_QT_GUI)

# -----------------------------------------------------------------------------
# Disable some warnings
# -----------------------------------------------------------------------------
IF (WIN32)
  IF (MSVC_VERSION GREATER 1400)
    # -------------------------------------------------------------------------
    # Disable deprecation warnings for standard C and STL functions in VS2005
    # and newer
    # -------------------------------------------------------------------------
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
  ENDIF (MSVC_VERSION GREATER 1400)
ELSEIF (UNIX)
  IF (CMAKE_COMPILER_IS_GNUCXX)
    # -------------------------------------------------------------------------
    # Disable "ignoring return value of ..." warnings from included gmvread.c
    # -------------------------------------------------------------------------
    INCLUDE(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-Wno-unused-result Wno-unused-result)
    IF (Wno-unused-result)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
    ENDIF (Wno-unused-result)
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
ENDIF (WIN32)


SET(SOURCES
cgio_helpers.cxx
cgio_helpers.h
vtkCGNSReaderInternal.cxx
vtkCGNSReaderInternal.h
vtkCGNSReader.cxx
vtkCGNSReader.h
)

# ---------------------------------------------------------------------------
# Combined plugin
# -----------------------------------------------------------------------------
IF (PARAVIEW_BUILD_QT_GUI)
  
  QT4_WRAP_CPP(MOC_SRCS pqCGNSReaderPanel.h)
  ADD_PARAVIEW_PROPERTY_GROUP_WIDGET(IFACES IFACE_SRCS
                               TYPE "cgns_widget_type"
                               CLASS_NAME pqCGNSReaderPanel
                              )

  ADD_PARAVIEW_PLUGIN(CGNSReader "1.0"
    REQUIRED_ON_SERVER
    REQUIRED_ON_CLIENT
    SERVER_MANAGER_XML CGNSReaderSM.xml
    SERVER_MANAGER_SOURCES ${SOURCES}
    GUI_INTERFACES ${IFACES}
    GUI_SOURCES ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS} pqCGNSReaderPanel.cxx
    )
ELSE (PARAVIEW_BUILD_QT_GUI)
  ADD_PARAVIEW_PLUGIN(CGNSReader "1.0"
    REQUIRED_ON_SERVER
    REQUIRED_ON_CLIENT
    SERVER_MANAGER_XML CGNSReaderSM.xml
    SERVER_MANAGER_SOURCES ${SOURCES})
ENDIF (PARAVIEW_BUILD_QT_GUI)


SET( CGNS_LINK_TO_HDF5 ON CACHE BOOL "use -lhdf5 for linking cgns to hdf5" )
MARK_AS_ADVANCED( CGNS_LINK_TO_HDF5 )

IF ( CGNS_LINK_TO_HDF5 )
  target_link_libraries(CGNSReader LINK_PRIVATE ${CGNS_LIBRARIES} -lhdf5)
ELSE (CGNS_LINK_TO_HDF5 )
  target_link_libraries(CGNSReader LINK_PRIVATE ${CGNS_LIBRARIES})
ENDIF ( CGNS_LINK_TO_HDF5 )

# Add testing if necessary
# IF (BUILD_TESTING)
#   ADD_SUBDIRECTORY(Testing)
# ENDIF (BUILD_TESTING)
