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

#project(SUPERLU C)
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

#add_library(SLUSRC OBJECT
#    ${SUPERLU_SRC_DIR}/sgssv.c
#    ${SUPERLU_SRC_DIR}/sgssvx.c
#    ${SUPERLU_SRC_DIR}/ssp_blas2.c
#    ${SUPERLU_SRC_DIR}/ssp_blas3.c
#    ${SUPERLU_SRC_DIR}/sgscon.c
#    ${SUPERLU_SRC_DIR}/slangs.c
#    ${SUPERLU_SRC_DIR}/sgsequ.c
#    ${SUPERLU_SRC_DIR}/slaqgs.c
#    ${SUPERLU_SRC_DIR}/spivotgrowth.c
#    ${SUPERLU_SRC_DIR}/sgsrfs.c
#    ${SUPERLU_SRC_DIR}/sgstrf.c
#    ${SUPERLU_SRC_DIR}/sgstrs.c
#    ${SUPERLU_SRC_DIR}/scopy_to_ucol.c
#    ${SUPERLU_SRC_DIR}/ssnode_dfs.c
#    ${SUPERLU_SRC_DIR}/ssnode_bmod.c
#    ${SUPERLU_SRC_DIR}/spanel_dfs.c
#    ${SUPERLU_SRC_DIR}/spanel_bmod.c
#    ${SUPERLU_SRC_DIR}/sreadhb.c
#    ${SUPERLU_SRC_DIR}/sreadrb.c
#    ${SUPERLU_SRC_DIR}/sreadtriple.c
#    ${SUPERLU_SRC_DIR}/scolumn_dfs.c
#    ${SUPERLU_SRC_DIR}/scolumn_bmod.c
#    ${SUPERLU_SRC_DIR}/spivotL.c
#    ${SUPERLU_SRC_DIR}/spruneL.c
#    ${SUPERLU_SRC_DIR}/smemory.c
#    ${SUPERLU_SRC_DIR}/sutil.c
#    ${SUPERLU_SRC_DIR}/smyblas2.c
#    ${SUPERLU_SRC_DIR}/sgsisx.c
#    ${SUPERLU_SRC_DIR}/sgsitrf.c
#    ${SUPERLU_SRC_DIR}/sldperm.c
#    ${SUPERLU_SRC_DIR}/ilu_sdrop_row.c
#    ${SUPERLU_SRC_DIR}/ilu_ssnode_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_scolumn_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_spanel_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_scopy_to_ucol.c
#    ${SUPERLU_SRC_DIR}/ilu_spivotL.c
#    ${SUPERLU_SRC_DIR}/sdiagonal.c
#    ${SUPERLU_SRC_DIR}/slacon2.c
#)

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

#add_library(CLUSRC OBJECT
#    ${SUPERLU_SRC_DIR}/scomplex.c
#    ${SUPERLU_SRC_DIR}/cgssv.c
#    ${SUPERLU_SRC_DIR}/cgssvx.c
#    ${SUPERLU_SRC_DIR}/csp_blas2.c
#    ${SUPERLU_SRC_DIR}/csp_blas3.c
#    ${SUPERLU_SRC_DIR}/cgscon.c
#    ${SUPERLU_SRC_DIR}/clangs.c
#    ${SUPERLU_SRC_DIR}/cgsequ.c
#    ${SUPERLU_SRC_DIR}/claqgs.c
#    ${SUPERLU_SRC_DIR}/cpivotgrowth.c
#    ${SUPERLU_SRC_DIR}/cgsrfs.c
#    ${SUPERLU_SRC_DIR}/cgstrf.c
#    ${SUPERLU_SRC_DIR}/cgstrs.c
#    ${SUPERLU_SRC_DIR}/ccopy_to_ucol.c
#    ${SUPERLU_SRC_DIR}/csnode_dfs.c
#    ${SUPERLU_SRC_DIR}/csnode_bmod.c
#    ${SUPERLU_SRC_DIR}/cpanel_dfs.c
#    ${SUPERLU_SRC_DIR}/cpanel_bmod.c
#    ${SUPERLU_SRC_DIR}/creadhb.c
#    ${SUPERLU_SRC_DIR}/creadrb.c
#    ${SUPERLU_SRC_DIR}/creadtriple.c
#    ${SUPERLU_SRC_DIR}/ccolumn_dfs.c
#    ${SUPERLU_SRC_DIR}/ccolumn_bmod.c
#    ${SUPERLU_SRC_DIR}/cpivotL.c
#    ${SUPERLU_SRC_DIR}/cpruneL.c
#    ${SUPERLU_SRC_DIR}/cmemory.c
#    ${SUPERLU_SRC_DIR}/cutil.c
#    ${SUPERLU_SRC_DIR}/cmyblas2.c
#    ${SUPERLU_SRC_DIR}/cgsisx.c
#    ${SUPERLU_SRC_DIR}/cgsitrf.c
#    ${SUPERLU_SRC_DIR}/cldperm.c
#    ${SUPERLU_SRC_DIR}/ilu_cdrop_row.c
#    ${SUPERLU_SRC_DIR}/ilu_csnode_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_ccolumn_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_cpanel_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_ccopy_to_ucol.c
#    ${SUPERLU_SRC_DIR}/ilu_cpivotL.c
#    ${SUPERLU_SRC_DIR}/cdiagonal.c
#    ${SUPERLU_SRC_DIR}/clacon2.c
#    ${SUPERLU_SRC_DIR}/scsum1.c
#    ${SUPERLU_SRC_DIR}/icmax1.c
#)

#add_library(ZLUSRC OBJECT 
#    ${SUPERLU_SRC_DIR}/dcomplex.c
#    ${SUPERLU_SRC_DIR}/zgssv.c
#    ${SUPERLU_SRC_DIR}/zgssvx.c
#    ${SUPERLU_SRC_DIR}/zsp_blas2.c
#    ${SUPERLU_SRC_DIR}/zsp_blas3.c
#    ${SUPERLU_SRC_DIR}/zgscon.c
#    ${SUPERLU_SRC_DIR}/zlangs.c
#    ${SUPERLU_SRC_DIR}/zgsequ.c
#    ${SUPERLU_SRC_DIR}/zlaqgs.c
#    ${SUPERLU_SRC_DIR}/zpivotgrowth.c
#    ${SUPERLU_SRC_DIR}/zgsrfs.c
#    ${SUPERLU_SRC_DIR}/zgstrf.c
#    ${SUPERLU_SRC_DIR}/zgstrs.c
#    ${SUPERLU_SRC_DIR}/zcopy_to_ucol.c
#    ${SUPERLU_SRC_DIR}/zsnode_dfs.c
#    ${SUPERLU_SRC_DIR}/zsnode_bmod.c
#    ${SUPERLU_SRC_DIR}/zpanel_dfs.c
#    ${SUPERLU_SRC_DIR}/zpanel_bmod.c
#    ${SUPERLU_SRC_DIR}/zreadhb.c
#    ${SUPERLU_SRC_DIR}/zreadrb.c
#    ${SUPERLU_SRC_DIR}/zreadtriple.c
#    ${SUPERLU_SRC_DIR}/zcolumn_dfs.c
#    ${SUPERLU_SRC_DIR}/zcolumn_bmod.c
#    ${SUPERLU_SRC_DIR}/zpivotL.c
#    ${SUPERLU_SRC_DIR}/zpruneL.c
#    ${SUPERLU_SRC_DIR}/zmemory.c
#    ${SUPERLU_SRC_DIR}/zutil.c
#    ${SUPERLU_SRC_DIR}/zmyblas2.c
#    ${SUPERLU_SRC_DIR}/zgsisx.c
#    ${SUPERLU_SRC_DIR}/zgsitrf.c
#    ${SUPERLU_SRC_DIR}/zldperm.c
#    ${SUPERLU_SRC_DIR}/ilu_zdrop_row.c
#    ${SUPERLU_SRC_DIR}/ilu_zsnode_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_zcolumn_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_zpanel_dfs.c
#    ${SUPERLU_SRC_DIR}/ilu_zcopy_to_ucol.c
#    ${SUPERLU_SRC_DIR}/ilu_zpivotL.c
#    ${SUPERLU_SRC_DIR}/zdiagonal.c
#    ${SUPERLU_SRC_DIR}/zlacon2.c
#    ${SUPERLU_SRC_DIR}/dzsum1.c
#    ${SUPERLU_SRC_DIR}/izmax1.c
#)

include_directories(${SUPERLU_SRC_DIR})
target_link_libraries(SUPERLU DLUSRC ALLAUX DZLAUX ${LAPACK_LIBRARIES})

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

export(TARGETS ALLAUX DZLAUX DLUSRC SUPERLU
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/SUPERLU.cmake"
)

