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
project(
                              OpenSees
)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/ETC/cmake)
set(OPS_EXTERN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/OTHER/)
set(OPS_EXTERNALS_DIR ${PROJECT_SOURCE_DIR}/OTHER/)

set(OPS_BUNDLED_DIR ${PROJECT_SOURCE_DIR}/OTHER/)
set(OPS_SRC_DIR ${PROJECT_SOURCE_DIR}/SRC/)
include(OpenSeesFunctions)
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS 
    G3 OpenSeesTcl OpenSeesMP OpenSeesSP OpenSeesPy
)

#==============================================================================
#                            OS Configuration
#
#==============================================================================
if (UNIX)
    include(OpenSeesDependenciesUnix)

    if(APPLE)
     message(STATUS ">>> MacOS")
    
    else()
       message(STATUS ">>> LINUX")
       add_compile_definitions(_LINUX _UNIX)
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ffloat-store")
    endif()
endif()

if(WIN32)
    include(OpenSeesDependenciesBundled)
    message(STATUS ">>> WIN32")
endif()

message("OPS >>> BLAS:    ${BLAS_LIBRARIES}")
message("OPS >>> CBLAS:   ${CBLAS_LIBRARY}\n")

message("OPS >>> LAPACK:  ${LAPACK_LIBRARIES}")

message("OPS >>> ARPACK:  ${ARPACK_LIBRARIES}")

message("OPS >>> SUPERLU: ${SUPERLU_LIBRARIES}")

message("OPS >>> TCL:     ${TCL_LIBRARY}")
set(TCL_LIBRARIES ${TCL_LIBRARY})

message("OPS >>> AMD:     ${AMD_LIBRARIES}")


#==============================================================================
#                            General Setup
#
#==============================================================================
#                             Properties

#define_property(TARGET
#    PROPERTY   OPS_INTERPRETER_GLOBAL #TODO
#    BRIEF_DOCS "Include functionality for using global interpreter"
#    FULL_DOCS  "..."
#)


# Fortran
#--------------------------------------
enable_language(Fortran)


# C++
#--------------------------------------

set(CMAKE_CXX_STANDARD 11)

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
# Global Includes
#----------------------------------------------------------------
 
# include paths to main abstract classes
 
# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(
    ${OPS_SRC_DIR}
    ${OPS_SRC_DIR}/matrix
    ${OPS_SRC_DIR}/handler
    ${OPS_SRC_DIR}/database
    ${OPS_SRC_DIR}/element
    ${OPS_SRC_DIR}/element/forceBeamColumn
    ${OPS_SRC_DIR}/element/nonlinearBeamColumn/matrixutil
    ${OPS_SRC_DIR}/coordTransformation
    ${OPS_SRC_DIR}/tagged
    ${OPS_SRC_DIR}/tagged/storage
    ${OPS_SRC_DIR}/recorder
    ${OPS_SRC_DIR}/renderer
    ${OPS_SRC_DIR}/damage
    ${OPS_SRC_DIR}/recorder/response
    ${OPS_SRC_DIR}/material
    ${OPS_SRC_DIR}/material/section
    ${OPS_SRC_DIR}/material/uniaxial
    ${OPS_SRC_DIR}/material/nD
    ${OPS_SRC_DIR}/graph/graph
    ${OPS_SRC_DIR}/graph/numberer
    ${OPS_SRC_DIR}/graph/partitioner
    ${OPS_SRC_DIR}/domain/component
    ${OPS_SRC_DIR}/domain/domain
    ${OPS_SRC_DIR}/domain/subdomain
    ${OPS_SRC_DIR}/domain/load
    ${OPS_SRC_DIR}/domain/loadBalancer
    ${OPS_SRC_DIR}/domain/pattern
    ${OPS_SRC_DIR}/domain/groundMotion
    ${OPS_SRC_DIR}/domain/node
    ${OPS_SRC_DIR}/domain/constraints
    ${OPS_SRC_DIR}/domain/region
    ${OPS_SRC_DIR}/analysis/algorithm
    ${OPS_SRC_DIR}/analysis/dof_grp
    ${OPS_SRC_DIR}/analysis/fe_ele
    ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo
    ${OPS_SRC_DIR}/analysis/algorithm/eigenAlgo
    ${OPS_SRC_DIR}/analysis/algorithm/domainDecompAlgo
    ${OPS_SRC_DIR}/analysis/analysis
    ${OPS_SRC_DIR}/analysis/integrator
    ${OPS_SRC_DIR}/analysis/handler
    ${OPS_SRC_DIR}/analysis/numberer
    ${OPS_SRC_DIR}/analysis/model
    ${OPS_SRC_DIR}/convergenceTest
    ${OPS_SRC_DIR}/modelbuilder
    ${OPS_SRC_DIR}/system_of_eqn
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE
    ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
    ${OPS_SRC_DIR}/actor/actor
    ${OPS_SRC_DIR}/actor/channel
    ${OPS_SRC_DIR}/actor/objectBroker
    ${OPS_SRC_DIR}/actor/message
)
# Temporary fix
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/incl)
include_directories(${TCL_INCLUDE_PATH})
include_directories(${MYSQL_INCLUDE_DIR})

file(GLOB_RECURSE ops_include_files CONFIGURE_DEPENDS 
	${OPS_SRC_DIR}/*.h
	${OPS_BUNDLED_DIR}/*.h
)
set(OPS_INCLUDE_DIRS "")
foreach(filepath ${ops_include_files})
    get_filename_component(dir_path ${filepath} PATH)
    set(OPS_INCLUDE_DIRS ${OPS_INCLUDE_DIRS} ${dir_path})
endforeach()
list(REMOVE_DUPLICATES OPS_INCLUDE_DIRS)

include_directories(${OPS_INCLUDE_DIRS})
#message(${OPS_INCLUDE_DIRS})


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
add_library(G3_Numerics            INTERFACE)
target_link_libraries(G3_Numerics INTERFACE
    ${ARPACK_LIBRARIES}
    ${SUPERLU_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    ${AMD_LIBRARIES}
)

# Core OpenSees
add_library(OPS_Tcl                STATIC)
add_library(OPS_Domain             OBJECT)
add_library(OPS_SysOfEqn           OBJECT)
add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_Damage             OBJECT)
add_library(OPS_Database           OBJECT)



# Optional Extensions
add_library(OPS_DRM                OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_PFEM               OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_ASDEA              OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
#add_library(OPS_Reliability        OBJECT EXCLUDE_FROM_ALL)


# Packaged libraries
add_library(G3) 
add_library(OpenSees         EXCLUDE_FROM_ALL )

# Executables
#------------------------------------------------------------------------------
add_executable(OpenSeesTcl   EXCLUDE_FROM_ALL 
    ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
    ${OPS_SRC_DIR}/tcl/tclMain.cpp
)
add_executable(OpenSeesSP    EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclAppInit.cpp)
add_executable(OpenSeesMP    EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclAppInit.cpp)
add_executable(OpenSeesPy    EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/interpreter/pythonMain.cpp)


# Set selected frontend to be built by default
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
    G3_Analysis 
    G3_ModelBuilder
    G3_Element
    G3_Actor 
    G3_Material
    G3_Domain
    G3_Matrix
    G3_Recorder
    G3_Handler
    G3_SysOfEqn
    G3_Tagged 
    G3_Utilities
    graph
    coordTransformation
    OPS_ConvergenceTest
    G3_Numerics 
)

#----------------------------
# OpenSees
#----------------------------
target_sources(OpenSees
    PRIVATE
	${OPS_SRC_DIR}/interpreter/OpenSeesBeamIntegrationCommands.cpp
	${OPS_SRC_DIR}/interpreter/OpenSeesCrdTransfCommands.cpp
	${OPS_SRC_DIR}/interpreter/OpenSeesFrictionModelCommands.cpp
	${OPS_SRC_DIR}/interpreter/OpenSeesMiscCommands.cpp
	${OPS_SRC_DIR}/interpreter/OpenSeesNDMaterialCommands.cpp
	${OPS_SRC_DIR}/interpreter/OpenSeesOutputCommands.cpp
        ${OPS_SRC_DIR}/interpreter/OpenSeesCommands.cpp
        ${OPS_SRC_DIR}/interpreter/OpenSeesElementCommands.cpp
	${OPS_SRC_DIR}/interpreter/OpenSeesSectionCommands.cpp
        ${OPS_SRC_DIR}/interpreter/OpenSeesUniaxialMaterialCommands.cpp
)
target_link_libraries(OpenSees
    #OPS_Paraview
    OPS_Renderer
    OPS_Database
    OPS_ConvergenceTest
    OPS_Element
    OPS_Material
    OPS_Material_YieldSurface
    OPS_Section_Repres
    OPS_Section_YieldSurface
    OPS_Damage
    OPS_SysOfEqn
    OPS_Analysis
    OPS_Domain
    OPS_Api
    G3
    ${MYSQL_LIBRARIES}
)

target_include_directories(OpenSees PUBLIC ${TCL_INCLUDE_PATH})

#----------------------------
# OPS_Tcl
#----------------------------
file(GLOB_RECURSE OPS_TCL_SOURCES CONFIGURE_DEPENDS 
	${OPS_SRC_DIR}/Tcl*.cpp 
)
# TODO: Properly list out Tcl files
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX "${OPS_SRC_DIR}/database/TclBerkeley.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*PFEM.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*/reliability/.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*/drm/.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*tclMain.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*tclAppInit.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*Interpreter.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*/petsc/.*")

#list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclFeViewer.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclVideoP.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclFeap.*")

list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*nonlinearBeamColumn.*")
list(FILTER OPS_TCL_SOURCES EXCLUDE REGEX ".*TclNewUniaxialMaterial.*")
# END

# Add sources to OPS_Tcl target
target_sources(OPS_Tcl PRIVATE 
	 ${OPS_TCL_SOURCES}
	"${OPS_SRC_DIR}/api/elementAPI_TCL.cpp"
	"${OPS_SRC_DIR}/tcl/commands.cpp"
)

target_include_directories(OPS_Tcl PUBLIC ${TCL_INCLUDE_PATH})
target_link_libraries(OPS_Tcl PRIVATE ${TCL_LIBRARIES})

target_compile_definitions(OPS_Tcl PUBLIC _TCL85)

opensees_add_cxx_flag(TARGETS OPS_Tcl
    GNU -fpermissive -fvisibility=hidden
)

#----------------------------
# OpenSeesTcl
#----------------------------

target_link_libraries(OpenSeesTcl 	
	OPS_Tcl ${OPS_Element_List} OPS_Thermal OpenSees ${CMAKE_DL_LIBS}
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
if (OPS_Use_Graphics)
    message("OPS >>> Including OPS_Graphics option")
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Renderer)
else()
    add_compile_definitions(_NOGRAPHICS)
endif()

#----------------------------
# Reliability
#----------------------------
if (OPS_Use_Reliability)
    add_compile_definitions(_RELIABILITY)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Reliability)
endif()

#----------------------------
# ASDEA
#----------------------------
if (OPS_Use_ASDEA)
    target_link_libraries(${OPS_FINAL_TARGET}  OPS_ASDEA)
endif()

#----------------------------
# PFEM
#----------------------------
if (OPS_Use_PFEM)
    message("OPS >>> Including OPS_PFEM option")
    #target_link_libraries(${OPS_FINAL_TARGET}  OPS_PFEM)
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
    message("OPS >>> Including OPS_Thermal option")
    #target_link_libraries(${OPS_FINAL_TARGET} PUBLIC OPS_Thermal)
    target_link_libraries(OpenSees OPS_Thermal)
endif()

#----------------------------
# DRM
#----------------------------
if (OPS_Use_DRM)
    message("OPS >>> Including OPS_DRM option")
    target_sources(OPS_Tcl 
	PRIVATE
	   ${OPS_SRC_DIR}/domain/pattern/drm/TclPatternCommand.cpp
    )
    target_link_libraries(${OPS_FINAL_TARGET} OPS_DRM)
endif()

