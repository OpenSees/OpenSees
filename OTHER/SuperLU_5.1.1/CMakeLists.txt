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
##-Extracted from accompanying Makefile-###############################
#  The files are organized as follows:
#
#   ALLAUX -- Auxiliary routines called from all precisions of SuperLU
#   SLUSRC -- Single precision real SuperLU routines
#   DLUSRC -- Double precision real SuperLU routines
#   CLUSRC -- Single precision complex SuperLU routines
#   ZLUSRC -- Double precision complex SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single double
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################
# Claudio Perez

set(SUPERLU_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/SRC")

if (WIN32)
  #target_compile_definitions(SUPERLU PUBLIC UpCase)
  add_definitions(-DUpCase)
endif()

add_library(SUPERLU)

add_library(DZLAUX  OBJECT
    ${SUPERLU_SRC_DIR}/dmach.c
)

add_library(ALLAUX OBJECT
    ${SUPERLU_SRC_DIR}/superlu_timer.c
    ${SUPERLU_SRC_DIR}/util.c
    ${SUPERLU_SRC_DIR}/memory.c
    ${SUPERLU_SRC_DIR}/get_perm_c.c
    ${SUPERLU_SRC_DIR}/mmd.c
    ${SUPERLU_SRC_DIR}/sp_coletree.c
    ${SUPERLU_SRC_DIR}/sp_preorder.c
    ${SUPERLU_SRC_DIR}/sp_ienv.c
    ${SUPERLU_SRC_DIR}/relax_snode.c
    ${SUPERLU_SRC_DIR}/heap_relax_snode.c
    ${SUPERLU_SRC_DIR}/colamd.c
    ${SUPERLU_SRC_DIR}/ilu_relax_snode.c
    ${SUPERLU_SRC_DIR}/ilu_heap_relax_snode.c
    ${SUPERLU_SRC_DIR}/mark_relax.c
    ${SUPERLU_SRC_DIR}/mc64ad.c
    ${SUPERLU_SRC_DIR}/qselect.c
    ${SUPERLU_SRC_DIR}/input_error.c
    ${SUPERLU_SRC_DIR}/smach.c
    #${SUPERLU_SRC_DIR}/dmach.c
)

add_library(DLUSRC OBJECT
    ${SUPERLU_SRC_DIR}/dgssv.c
    ${SUPERLU_SRC_DIR}/dgssvx.c
    ${SUPERLU_SRC_DIR}/dsp_blas2.c
    ${SUPERLU_SRC_DIR}/dsp_blas3.c
    ${SUPERLU_SRC_DIR}/dgscon.c
    ${SUPERLU_SRC_DIR}/dlangs.c
    ${SUPERLU_SRC_DIR}/dgsequ.c
    ${SUPERLU_SRC_DIR}/dlaqgs.c
    ${SUPERLU_SRC_DIR}/dpivotgrowth.c
    ${SUPERLU_SRC_DIR}/dgsrfs.c
    ${SUPERLU_SRC_DIR}/dgstrf.c
    ${SUPERLU_SRC_DIR}/dgstrs.c
    ${SUPERLU_SRC_DIR}/dcopy_to_ucol.c
    ${SUPERLU_SRC_DIR}/dsnode_dfs.c
    ${SUPERLU_SRC_DIR}/dsnode_bmod.c
    ${SUPERLU_SRC_DIR}/dpanel_dfs.c
    ${SUPERLU_SRC_DIR}/dpanel_bmod.c
    ${SUPERLU_SRC_DIR}/dreadhb.c
    ${SUPERLU_SRC_DIR}/dreadrb.c
    ${SUPERLU_SRC_DIR}/dreadtriple.c
    ${SUPERLU_SRC_DIR}/dcolumn_dfs.c
    ${SUPERLU_SRC_DIR}/dcolumn_bmod.c
    ${SUPERLU_SRC_DIR}/dpivotL.c
    ${SUPERLU_SRC_DIR}/dpruneL.c
    ${SUPERLU_SRC_DIR}/dmemory.c
    ${SUPERLU_SRC_DIR}/dutil.c
    ${SUPERLU_SRC_DIR}/dmyblas2.c
    ${SUPERLU_SRC_DIR}/dgsisx.c
    ${SUPERLU_SRC_DIR}/dgsitrf.c
    ${SUPERLU_SRC_DIR}/dldperm.c
    ${SUPERLU_SRC_DIR}/ilu_ddrop_row.c
    ${SUPERLU_SRC_DIR}/ilu_dsnode_dfs.c
    ${SUPERLU_SRC_DIR}/ilu_dcolumn_dfs.c
    ${SUPERLU_SRC_DIR}/ilu_dpanel_dfs.c
    ${SUPERLU_SRC_DIR}/ilu_dcopy_to_ucol.c
    ${SUPERLU_SRC_DIR}/ilu_dpivotL.c
    ${SUPERLU_SRC_DIR}/ddiagonal.c
    ${SUPERLU_SRC_DIR}/dlacon2.c
     ## dgstrsL.c dgstrsU.c
)


include_directories(${SUPERLU_SRC_DIR})
target_link_libraries(SUPERLU DLUSRC ALLAUX DZLAUX ${LAPACK_LIBRARIES})

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

export(TARGETS ALLAUX DZLAUX DLUSRC SUPERLU
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/SUPERLU.cmake"
)

