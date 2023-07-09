#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#------------------------------------------------------------------------------
#
#  NOTE: THIS BUILD SOLUTION HAS BEEN DEVELOPED FOR USE WITH THE OPENSEES 
#     FRAMEWORK AND MAY DIFFER CONSIDERABLY FROM THE UPSTREAM LIBRARY.
# 
#           AMD Version 2.4, Copyright (c) 1996-2013 by Timothy A. Davis
#                  people.engr.tamu.edu/davis/suitesparse.html
#==============================================================================
# Claudio Perez

project(AMD)
#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------
set(AMD_SOURCES
  amd_aat.c
  amd_1.c
  amd_2.c
  amd_dump.c
  amd_postorder.c
  amd_post_tree.c
  amd_defaults.c
  amd_order.c
  amd_control.c
  amd_info.c
  amd_valid.c
  amd_preprocess.c
)

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

add_library(AMD_DINT ${AMD_SOURCES})
target_compile_definitions(AMD_DINT PRIVATE DINT)
target_compile_options(AMD_DINT PRIVATE -fexceptions -fPIC -I.)

add_library(AMD_DLNG ${AMD_SOURCES})
target_compile_definitions(AMD_DLNG PRIVATE DLONG)
target_compile_options(AMD_DLNG PRIVATE -fexceptions -fPIC -I.)

#-------------------------------------------------------------------------------

add_library(AMD)
target_sources(AMD PRIVATE
  SuiteSparse_config.c
  $<TARGET_OBJECTS:AMD_DINT>
  $<TARGET_OBJECTS:AMD_DLNG>
)
#target_link_libraries(AMD PUBLIC AMD_DINT AMD_DLNG)


