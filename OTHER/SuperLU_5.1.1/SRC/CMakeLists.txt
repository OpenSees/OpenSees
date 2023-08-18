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

if (WIN32)
  #target_compile_definitions(SUPERLU PUBLIC UpCase)
  add_definitions(-DUpCase)
endif()

add_library(SUPERLU)

add_library(DZLAUX  OBJECT
     dmach.c
)

add_library(ALLAUX OBJECT
    superlu_timer.c
    util.c
    memory.c
    get_perm_c.c
    mmd.c
    sp_coletree.c
    sp_preorder.c
    sp_ienv.c
    relax_snode.c
    heap_relax_snode.c
    colamd.c
    ilu_relax_snode.c
    ilu_heap_relax_snode.c
    mark_relax.c
    mc64ad.c
    qselect.c
    input_error.c
    smach.c
    #dmach.c
)

add_library(DLUSRC OBJECT
    dgssv.c
    dgssvx.c
    dsp_blas2.c
    dsp_blas3.c
    dgscon.c
    dlangs.c
    dgsequ.c
    dlaqgs.c
    dpivotgrowth.c
    dgsrfs.c
    dgstrf.c
    dgstrs.c
    dcopy_to_ucol.c
    dsnode_dfs.c
    dsnode_bmod.c
    dpanel_dfs.c
    dpanel_bmod.c
    dreadhb.c
    dreadrb.c
    dreadtriple.c
    dcolumn_dfs.c
    dcolumn_bmod.c
    dpivotL.c
    dpruneL.c
    dmemory.c
    dutil.c
    dmyblas2.c
    dgsisx.c
    dgsitrf.c
    dldperm.c
    ilu_ddrop_row.c
    ilu_dsnode_dfs.c
    ilu_dcolumn_dfs.c
    ilu_dpanel_dfs.c
    ilu_dcopy_to_ucol.c
    ilu_dpivotL.c
    ddiagonal.c
    dlacon2.c
     ## dgstrsL.c dgstrsU.c
)

target_link_libraries(SUPERLU DLUSRC ALLAUX DZLAUX ${LAPACK_LIBRARIES})

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

export(TARGETS ALLAUX DZLAUX DLUSRC SUPERLU
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/SUPERLU.cmake"
)

