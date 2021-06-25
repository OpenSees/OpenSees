cmake_minimum_required(VERSION 3.18)
#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
# Major sections:
# - General Setup
#   Project wide includes, properties, libraries, etc
#
# - Target definitions
#   Define component libraries
#
# - Target Configuration
#   Specify target-specific compiler options, linking, etc
#
# - Option Application
#==============================================================================
#project(
#                              OpenSees
#)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/ETC/cmake)
set(OPS_EXTERN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/OTHER/)
set(OPS_EXTERNALS_DIR ${PROJECT_SOURCE_DIR}/OTHER/)

set(OPS_BUNDLED_DIR ${PROJECT_SOURCE_DIR}/OTHER/)
set(OPS_SRC_DIR ${PROJECT_SOURCE_DIR}/SRC/)

include(${PROJECT_SOURCE_DIR}/Conf.cmake)
#message("OPS >>> BLAS:    ${BLAS_LIBRARIES}")
#message("OPS >>> CBLAS:   ${CBLAS_LIBRARY}")
#message("OPS >>> LAPACK:  ${LAPACK_LIBRARIES}")
#message("OPS >>> ARPACK:  ${ARPACK_LIBRARIES}")
#message("OPS >>> SUPERLU: ${SUPERLU_LIBRARIES}")
#message("OPS >>> TCL:     ${TCL_LIBRARIES}")
 

set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS 
    G3 OpenSeesTcl OpenSeesMP OpenSeesSP OpenSeesPy
)

#==============================================================================
#                            General Setup
#
#==============================================================================
##                            Properties
#
#define_property(TARGET
#    PROPERTY   OPS_INTERPRETER_GLOBAL #TODO
#    BRIEF_DOCS "Include functionality for using global interpreter"
#    FULL_DOCS  "..."
#)


#----------------------------------------------------------------
# Compilers
#----------------------------------------------------------------

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------
# Warnings
opensees_add_cxx_flag(
    GNU  -Wall
    MSVC /W4
)

# Floating-point
opensees_add_cxx_flag(
    GNU  -ffloat-store
    MSVC /fp:precise
)

#----------------------------------------------------------------
# OS Configuration
#----------------------------------------------------------------

if(APPLE)
 message(STATUS ">>> MacOS")
endif()

if(UNIX AND NOT APPLE)
   message(STATUS ">>> LINUX")
   add_compile_definitions(_LINUX _UNIX)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ffloat-store")
endif()

if(WIN32)
 message(STATUS ">>> WIN32")
endif()

#----------------------------------------------------------------
# Global Includes
#----------------------------------------------------------------
#
# include paths to main abstract classes
#
# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(src
#    src/base
#    src/matrix
#    src/handler
#    src/database
#    #src/element
    src/element/forceBeamColumn
    src/element/nonlinearBeamColumn/matrixutil
#    src/coordTransformation
#    src/tagged
#    src/tagged/storage
#    src/recorder
#    src/renderer
#    src/damage
#    src/recorder/response
#    src/material
#    src/material/section
#    src/material/uniaxial
#    #src/material/nD
#    src/graph/graph
#    src/graph/numberer
#    src/graph/partitioner
#    src/domain/component
#    src/domain/domain
#    src/domain/subdomain
#    src/domain/load
#    src/domain/loadBalancer
#    src/domain/pattern
#    src/domain/groundMotion
#    src/domain/node
#    src/domain/constraints
#    src/domain/region
#    src/analysis/algorithm
#    src/analysis/dof_grp
#    src/analysis/fe_ele
#    src/analysis/algorithm/equiSolnAlgo
#    src/analysis/algorithm/eigenAlgo
#    src/analysis/algorithm/domainDecompAlgo
#    src/analysis/analysis
#    src/analysis/integrator
#    src/analysis/handler
#    src/analysis/numberer
#    src/analysis/model
#    src/convergenceTest
#    src/modelbuilder
#    src/system_of_eqn
#    src/system_of_eqn/linearSOE
#    src/system_of_eqn/eigenSOE
#    src/actor/actor
#    src/actor/channel
#    src/actor/objectBroker
#    src/actor/message
)
#


#==============================================================================
#                            Define Targets
#
#==============================================================================
add_library(G3_Matrix              OBJECT)
add_library(G3_Actor               OBJECT)
add_library(G3_Domain              OBJECT)
add_library(G3_SysOfEqn            OBJECT)
add_library(G3_Handler             OBJECT)
add_library(G3_Analysis            OBJECT)
add_library(G3_Recorder            OBJECT)
add_library(G3_Tagged              OBJECT)
add_library(G3_Utilities           OBJECT)
add_library(G3_Element             OBJECT)
add_library(G3_Material            OBJECT)
add_library(G3_ModelBuilder        OBJECT)
add_library(G3_Numerics            OBJECT 
    ${ARPACK_LIBRARIES}
    ${SUPERLU_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

# Core OpenSees
add_library(OPS_Tcl                STATIC)
add_library(OPS_Domain             OBJECT)
add_library(OPS_SysOfEqn           OBJECT)
add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)

# Optional Extensions
add_library(OPS_DRM                OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_PFEM               OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_ASDEA              OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Reliability        OBJECT EXCLUDE_FROM_ALL)


add_library(G3 )
add_library(OpenSees         EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/interpreter/OpenSeesCommands.cpp)
add_executable(OpenSeesTcl   EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclAppInit.cpp)
add_executable(OpenSeesSP    EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclAppInit.cpp)
add_executable(OpenSeesMP    EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclAppInit.cpp)
add_executable(OpenSeesPy    EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/interpreter/pythonMain.cpp)

# Build selected frontend by default
set_target_properties(${OPS_FINAL_TARGET} PROPERTIES EXCLUDE_FROM_ALL OFF)
# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})

#==============================================================================
#                            Configure targets
#
#==============================================================================
#----------------------------
# G3
#----------------------------
target_link_libraries(G3 
    G3_Actor 
    G3_ModelBuilder
    G3_Analysis 
    G3_Element
    G3_Domain
    G3_Matrix
    G3_Recorder
    G3_Handler
    G3_SysOfEqn
    G3_Tagged 
    G3_Utilities
    graph 
    OPS_ConvergenceTest
    G3_Numerics 
)
#----------------------------
# OpenSees
#----------------------------
target_link_libraries(OpenSees
    #OPS_Paraview
    OPS_Element
    OPS_Domain
    OPS_SysOfEqn
    OPS_Analysis
    G3 
)

target_include_directories(OpenSees PUBLIC ${TCL_INCLUDE_DIRS})

#----------------------------
# OPS_Tcl
#----------------------------
file(GLOB_RECURSE OPS_TCL_SOURCES CONFIGURE_DEPENDS 
	${OPS_SRC_DIR}/Tcl*.cpp 
)
# TODO: Properly list out Tcl files
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX "${OPS_SRC_DIR}/database/.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*PFEM.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*/reliability/.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*/drm/.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*tclMain.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*tclAppInit.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*Interpreter.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*/petsc/.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclFeViewer.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclVideoP.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*nonlinearBeamColumn.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclNewUniaxialMaterial.*")
# END

# Add sources to OPS_Tcl target
target_sources(OPS_Tcl PRIVATE 
	 ${OPS_TCL_SOURCES}
	"${OPS_SRC_DIR}/api/elementAPI_TCL.cpp"
	"${OPS_SRC_DIR}/tcl/commands.cpp"
)
message("TCL INCLUDE: ${TCL_INCLUDE_DIRS}")
target_include_directories(OPS_Tcl PUBLIC ${TCL_INCLUDE_DIRS})
target_link_libraries(OPS_Tcl PRIVATE ${TCL_LIBRARIES})

opensees_add_cxx_flag(TARGETS OPS_Tcl
    GNU -fpermissive -fvisibility=hidden
)

#----------------------------
# OpenSeesTcl
#----------------------------
target_link_libraries(OpenSeesTcl
	PRIVATE OpenSees OPS_Tcl #${TCL_LIBRARY}
)

#----------------------------
# OpenSeesMP
#----------------------------
target_compile_definitions(OpenSeesMP
    PUBLIC _OPS_PARALLEL_PROCESSING 
)
target_link_libraries(OpenSeesMP OpenSees)


#----------------------------
# OpenSeesSP
#----------------------------
target_compile_definitions(OpenSeesSP 
    PUBLIC _OPS_PARALLEL_INTERPRETERS
)
#----------------------------
# OpenSeesPy
#----------------------------
opensees_add_cxx_flag(TARGETS OpenSeesPy
    GNU  -fPIC
)

#==============================================================================
#                            Apply Options
#
#==============================================================================
#----------------------------
# FMK
#----------------------------
if(FMK)
    add_compile_definitions(
        _HAVE_Damage2p    
        _HAVE_PSUMAT
        _HAVE_PML
        _FILIP_LHNMYS
    )    
endif()

#----------------------------
# Renderer
#----------------------------
if (OPS_Use_Renderer)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Renderer)
endif()

#----------------------------
# Reliability
#----------------------------
if (OPS_Use_RELIABILITY)
    add_compile_definitions(_RELIABILITY)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Reliability)
endif()

#----------------------------
# ASDEA
#----------------------------
if (OPS_Use_ASDEA)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_ASDEA)
endif()

#----------------------------
# PFEM
#----------------------------
if (OPS_Use_PFEM)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_PFEM)
endif()

#----------------------------
# HDF5
#----------------------------
if(OPS_Use_HDF5)
   find_package(HDF5)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(_hdf5_libs hdf5 hdf5_cpp)
    add_compile_definitions(-D_H5DRM)
    else()
        message(STATUS "OPS >>> Could not find HDF5")
    endif()
endif()

#----------------------------
# Thermal
#----------------------------
if (OPS_Use_Thermal)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Thermal)
endif()

#----------------------------
# DRM
#----------------------------
if (OPS_Use_DRM)
    target_sources(OPS_Tcl 
	PRIVATE
	   ${OPS_SRC_DIR}/domain/pattern/drm/TclPatternCommand.cpp
    )
    target_link_libraries(${OPS_FINAL_TARGET} OPS_DRM)
endif()

