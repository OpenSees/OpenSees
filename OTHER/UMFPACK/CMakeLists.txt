#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
#        (Licensing @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
#
# NOTE: THIS BUILD SOLUTION HAS BEEN DEVELOPED FOR USE WITH THE OPENSEES 
# FRAMEWORK AND MAY DIFFER CONSIDERABLY FROM THE UPSTREAM LIBRARY.
# 
#                         UMFPACK V5.7.1 (Oct 10, 2014)
#                  people.engr.tamu.edu/davis/suitesparse.html
#==============================================================================
# Claudio Perez
cmake_minimum_required(VERSION 3.16)


include_directories(. ../AMD)
add_compile_options(-fexceptions -fPIC)
add_compile_definitions(NCHOLMOD)

function (add_umflib libname libsource)
  #list(APPEND UMF_CREATED ${libname})
  set(UMF_CREATED ${UMF_CREATED} ${libname} PARENT_SCOPE)
  add_library(${libname} OBJECT)
  target_sources(${libname} PRIVATE ${libsource})
  foreach(def ${ARGN})
    target_compile_definitions(${libname} PRIVATE ${def})
  endforeach()
endfunction()

#----------------------------------------------------------------------
# Files with different variants
#----------------------------------------------------------------------
# non-user-callable, created from umf_ltsolve.c, umf_utsolve.c,
# umf_triplet.c, and umf_assemble.c , with int/SuiteSparse_long
# and real/complex versions:
#          library name           source file      definitions...
add_umflib(umf_triplet_nomap_nox  umf_triplet.c)
add_umflib(umf_triplet_map_x      umf_triplet.c    DO_MAP  DO_VALUES)
add_umflib(umf_triplet_map_nox    umf_triplet.c    DO_MAP)
add_umflib(umf_triplet_nomap_x    umf_triplet.c            DO_VALUES)
add_umflib(umf_ltsolve            umf_ltsolve.c)
add_umflib(umf_lhsolve            umf_ltsolve.c    CONJUGATE_SOLVE)
add_umflib(umf_utsolve            umf_utsolve.c)
add_umflib(umf_uhsolve            umf_utsolve.c    CONJUGATE_SOLVE)
add_umflib(umf_assemble_di        umf_assemble.c)
add_umflib(umf_assemble_di_fixq   umf_assemble.c   FIXQ)
add_umflib(umf_store_lu           umf_store_lu.c)
add_umflib(umf_store_lu_drop      umf_store_lu.c   DROP)
add_umflib(umfpack_solve          umfpack_solve.c)
add_umflib(umfpack_di_wsolve      umfpack_solve.c  WSOLVE)

#----------------------------------------------------------
# source files
#----------------------------------------------------------
# non-user-callable umf_*.[ch] files:
set(UMFCH
  umf_assemble.c
  umf_blas3_update.c
  umf_build_tuples.c
  umf_create_element.c
  umf_extend_front.c
  umf_garbage_collection.c
  umf_get_memory.c
  umf_init_front.c
  umf_kernel.c
  umf_kernel_init.c
  umf_kernel_wrapup.c
  umf_local_search.c
  umf_lsolve.c
  umf_ltsolve.c
  umf_mem_alloc_element.c
  umf_mem_alloc_head_block.c
  umf_mem_alloc_tail_block.c
  umf_mem_free_tail_block.c
  umf_mem_init_memoryspace.c
  umf_report_vector.c
  umf_row_search.c
  umf_scale_column.c
  umf_set_stats.c
  umf_solve.c
  umf_symbolic_usage.c
  umf_transpose.c
  umf_tuple_lengths.c
  umf_usolve.c
  umf_utsolve.c
  umf_valid_numeric.c
  umf_valid_symbolic.c
  umf_grow_front.c
  umf_start_front.c  
  umf_store_lu.c
  umf_scale.c
)

# non-user-callable umf_*.[ch] files, int/SuiteSparse_long versions only
# (no real/complex):
set(UMFINT
  umf_analyze.c
  umf_apply_order.c
  umf_colamd.c
  umf_free.c
  umf_fsize.c
  umf_is_permutation.c
  umf_malloc.c
  umf_realloc.c
  umf_report_perm.c
  umf_singletons.c
  umf_cholmod.c
)

# user-callable umfpack_*.[ch] files (int/SuiteSparse_long and real/complex):
set(UMFPACK_USER
  umfpack_col_to_triplet.c
  umfpack_defaults.c
  umfpack_free_numeric.c
  umfpack_free_symbolic.c
  umfpack_get_numeric.c
  umfpack_get_lunz.c
  umfpack_get_symbolic.c
  umfpack_get_determinant.c
  umfpack_numeric.c
  umfpack_qsymbolic.c
  umfpack_report_control.c
  umfpack_report_info.c
  umfpack_report_matrix.c
  umfpack_report_numeric.c
  umfpack_report_perm.c
  umfpack_report_status.c
  umfpack_report_symbolic.c
  umfpack_report_triplet.c
  umfpack_report_vector.c
  umfpack_symbolic.c
  umfpack_transpose.c
  umfpack_triplet_to_col.c
  umfpack_scale.c
  umfpack_load_numeric.c
  umfpack_save_numeric.c
  umfpack_load_symbolic.c
  umfpack_save_symbolic.c
)  

# user-callable, created from umfpack_solve.c (umfpack_wsolve.h exists, though):
# with int/SuiteSparse_long and real/complex versions:
#add_library(UMFPACKW PUBLIC umfpack_wsolve)


# user-callable, only one version for int/SuiteSparse_long,
# real/complex, *.[ch] files:
set(GENERIC umfpack_timer.c umfpack_tictoc.c)

#----------------------------------------------------------------------
# Component libraries
#----------------------------------------------------------------------
# non-user-callable, int/SuiteSparse_long and real/complex versions:
add_library(UMF OBJECT ${UMFCH})

add_library(UMFPACK_GN OBJECT ${GENERIC})
add_library(UMFPACK_II OBJECT ${UMFINT})
add_library(UMFPACK_DI OBJECT ${UMF} ${UMFPACK_USER})
target_compile_definitions(UMF         PRIVATE DINT)
target_compile_definitions(UMFPACK_II  PRIVATE DINT)
target_compile_definitions(UMFPACK_DI  PRIVATE DINT)

#----------------------------------------------------------------------
# Create the libumfpack.a library
#----------------------------------------------------------------------
add_library(UMFPACK #INTERFACE)
#target_link_libraries(UMFPACK INTERFACE
#target_sources(umfpack
  #"-Wl,--whole-archive"
  $<TARGET_OBJECTS:UMFPACK_DI>
  $<TARGET_OBJECTS:UMFPACK_II>
  $<TARGET_OBJECTS:UMFPACK_GN>
  $<TARGET_OBJECTS:UMF>
  #${UMF_CREATED}
  #"-Wl,--no-whole-archive"
  #${AMD_LIBRARIES}
)
target_link_libraries(UMFPACK PUBLIC ${UMF_CREATED})
set_target_properties(UMFPACK PROPERTIES LINKER_LANGUAGE C)

