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
#
# Major sections:
# - General Setup
#   Project wide includes, properties, libraries, etc
#
# - Target definitions
#   Declarations/definitions of component libraries
#
# - Target Configuration
#   Specify target-specific compiler options, linking, etc
#
# - Option Application
#
#==============================================================================
cmake_minimum_required(VERSION 3.16)
project(OpenSees C CXX Fortran)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(OPS_EXTERN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/OTHER/")

set(OPS_BUNDLED_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_SRC_DIR "${PROJECT_SOURCE_DIR}/SRC/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/ETC/cmake")
include(OpenSeesFunctions)

# pseudo target modeling all exteral packages
add_library(OPS_External_packages INTERFACE)
#add_library(OPS_Element_libraries INTERFACE)

# include user config
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

message("OPS >>> BLAS:    ${BLAS_LIBRARIES}")
message("OPS >>> CBLAS:   ${CBLAS_LIBRARY}\n")
message("OPS >>> LAPACK:  ${LAPACK_LIBRARIES}")
message("OPS >>> ARPACK:  ${ARPACK_LIBRARIES}")
message("OPS >>> SUPERLU: ${SUPERLU_LIBRARIES}")
message("OPS >>> TCL:     ${TCL_LIBRARY}")
set(TCL_LIBRARIES ${TCL_LIBRARY})
message("OPS >>> AMD:     ${AMD_LIBRARIES}")

# define user-selectable options for end target
set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS 
    G3 OpenSeesTcl OpenSeesMP OpenSeesSP OpenSeesPy)

#==============================================================================
#                            OS Configuration
#
#==============================================================================
add_library(OPS_OS_Specific_libs INTERFACE)

if (UNIX)
    include(OpenSeesDependenciesUnix)
    if(APPLE)
     message(STATUS ">>> MacOS")
    else()
       message(STATUS ">>> LINUX")
       add_compile_definitions(_LINUX _UNIX _TCL85)
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ffloat-store")
    endif()
endif()

if(WIN32) # NOTE: this will execute for both 32-bit and 64-bit builds.
    include(OpenSeesDependenciesWin)
    target_link_libraries(OPS_OS_Specific_libs INTERFACE wsock32 ws2_32)
    add_compile_definitions(_WIN32)
    message(STATUS ">>> WIN32")
endif()



#==============================================================================
#                            General Setup
#
#==============================================================================
#----------------------------------------------------------------
# Compilers
#---------------------------------------------------------------- 

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------
#add_link_options(-rdynamic)
set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cqls <TARGET> <LINK_FLAGS> <OBJECTS>")
#add_compile_options(-g)

# Warnings
opensees_add_cxx_flag(
    GNU  -Wall
    MSVC /W0
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

include_directories(${TCL_INCLUDE_PATH})
include_directories(${MYSQL_INCLUDE_DIR})

# TODO: Temporary fix; include all directories with .h files through glob
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
# end TODO

#==============================================================================
#                            Define Targets
#
#==============================================================================
add_library(OPS_Numerics           INTERFACE)
target_link_libraries(OPS_Numerics INTERFACE
  ${ARPACK_LIBRARIES} 
  ${SUPERLU_LIBRARIES}
  ${AMD_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} 
)
# Core OpenSees
add_library(OPS_Matrix             OBJECT)
add_library(OPS_Actor              OBJECT)
add_library(OPS_ObjectBroker       OBJECT)
add_library(OPS_Handler            OBJECT)
add_library(OPS_Recorder           OBJECT)
add_library(OPS_Tagged             OBJECT)
add_library(OPS_Utilities          OBJECT)
add_library(OPS_ModelBuilder       OBJECT)
add_library(OPS_Domain             OBJECT)

add_library(OPS_SysOfEqn           OBJECT)

add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_Damage             OBJECT)
add_library(OPS_Database           OBJECT)
add_library(OPS_Interp_Tcl         STATIC)

# Optional Extensions
add_library(OPS_DRM                OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_ASDEA              OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
#add_library(OPS_Reliability        OBJECT EXCLUDE_FROM_ALL)

opensees_library(OPS_SysOfEqn_UMF LINK ${UMFPACK_LIBRARIES})
opensees_library(OPS_PFEM         LINK ${CSPARSE_LIBRARIES} ${UMFPACK_LIBRARIES})


# Packaged libraries
add_library(G3) 
add_library(OpenSees         EXCLUDE_FROM_ALL)

# Executables and final targets
#------------------------------------------------------------------------------
# Adding EXCLUDE_FROM_ALL prevents these from being made by default.

add_executable(OpenSeesTcl   EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
)

add_executable(OpenSeesSP    EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
)

add_executable(OpenSeesMP    EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
)

add_executable(OpenSeesPy    EXCLUDE_FROM_ALL
  ${OPS_SRC_DIR}/interpreter/pythonMain.cpp
)

# Set selected executable to be built by default
set_target_properties(${OPS_FINAL_TARGET} PROPERTIES EXCLUDE_FROM_ALL OFF)

# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})

#==============================================================================
#                            Configure targets
#
#==============================================================================
target_link_libraries(OPS_Numerics     INTERFACE ${OPS_Numlib_List})
message("Numlibs: ${OPS_Numlib_List}")
target_link_libraries(OPS_Element      PRIVATE ${OPS_Element_List})
target_link_libraries(OPS_ObjectBroker PRIVATE ${OPS_Element_List})
target_link_libraries(OPS_Interp_Tcl   PRIVATE ${OPS_Element_List})

#----------------------------
# G3
#----------------------------
target_link_libraries(G3 
    OPS_Matrix
    OPS_Analysis 
    OPS_ModelBuilder
    OPS_Domain
    OPS_ConvergenceTest
    OPS_Element
    OPS_Material
    OPS_Recorder
    OPS_Handler
    OPS_SysOfEqn
    ${OPS_SysOfEqn_List}
    OPS_Tagged 
    OPS_Utilities
    graph
    coordTransformation
    OPS_Actor 
    OPS_ObjectBroker
    OPS_Numerics 
)

#----------------------------
# OpenSees
#----------------------------
#<<<<<<< Updated upstream
#=======
#target_sources(OpenSees
#    PRIVATE
#    "${OPS_SRC_DIR}/interpreter/OpenSeesBeamIntegrationCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesCrdTransfCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesFrictionModelCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesMiscCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesNDMaterialCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesOutputCommands.cpp"
#  #"${OPS_SRC_DIR}/interpreter/OpenSeesCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesElementCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesSectionCommands.cpp"
#    "${OPS_SRC_DIR}/interpreter/OpenSeesUniaxialMaterialCommands.cpp"
#)
#>>>>>>> Stashed changes

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
    OPS_OS_Specific_libs
    G3
    ${MYSQL_LIBRARIES}
)

target_include_directories(OpenSees PUBLIC ${TCL_INCLUDE_PATH})

#----------------------------
# OPS_Interp_Tcl
#----------------------------
# Add sources to OPS_Interp_Tcl target

target_include_directories(OPS_Interp_Tcl PUBLIC ${TCL_INCLUDE_PATH})
target_link_libraries(OPS_Interp_Tcl PRIVATE ${TCL_LIBRARIES})

target_compile_definitions(OPS_Interp_Tcl PUBLIC _TCL85)

opensees_add_cxx_flag(TARGETS OPS_Interp_Tcl
    GNU -fpermissive -fvisibility=hidden
)

#----------------------------
# OpenSeesTcl
#----------------------------
target_link_libraries(OpenSeesTcl 	
	OPS_Interp_Tcl #${OPS_Element_List} 
  OPS_Thermal OpenSees ${CMAKE_DL_LIBS}
)

#----------------------------
# OpenSeesMP
#----------------------------
target_compile_definitions(OpenSeesMP
    PUBLIC _PARALLEL_PROCESSING 
)
target_link_libraries(OpenSeesMP OpenSees)


#----------------------------
# OpenSeesSP
#----------------------------
target_compile_definitions(OpenSeesSP 
    PUBLIC _PARALLEL_INTERPRETERS
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
# Extensions
#----------------------------
message("OPS >>> Configuring OpenSees extensions")
foreach(extension IN LISTS OPS_SysOfEqn_List OPS_Element_List OPS_Extension_List)
    string(TOUPPER "${extension}" ext_flag) 
    string(REGEX REPLACE "^OPS_" "OPSDEF_" ext_flag "${ext_flag}")
    # message("    Adding macro definition '${ext_flag}'")
    add_compile_definitions(${ext_flag})
endforeach()

# PFEM
#----------------------------
if (OPS_Use_PFEM)
  add_subdirectory("${OPS_EXTERNALS_DIR}/Tetgen")
  #add_subdirectory("${OPS_EXTERNALS_DIR}/UMFPACK")
  #add_subdirectory("${OPS_EXTERNALS_DIR}/CSPARSE")
  target_link_libraries(${OPS_FINAL_TARGET}
    OPS_Element_PFEMElement OPS_PFEM  #tet
  )
endif()

# Renderer
#----------------------------
if (OPS_Use_Graphics)
    message("OPS >>> Including OPS_Graphics option")
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Renderer)
else()
    add_compile_definitions(_NOGRAPHICS)
endif()


# Reliability
#----------------------------
if ("OPS_Reliability" IN_LIST OPS_Extension_List)
    add_compile_definitions(_RELIABILITY)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Reliability)
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
# DRM
#----------------------------
if (OPS_Use_DRM)
    message("OPS >>> Including OPS_DRM option")
    target_sources(OPS_Interp_Tcl 
      PRIVATE
	      ${OPS_SRC_DIR}/domain/pattern/drm/TclPatternCommand.cpp
    )
    target_link_libraries(${OPS_FINAL_TARGET} OPS_DRM)
endif()

