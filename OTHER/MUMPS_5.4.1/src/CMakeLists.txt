
target_sources(MUMPS_C PUBLIC
		       mumps_io_basic.c
		       mumps_numa.c
		       mumps_scotch_int.c
		       dmumps_gpu.c
		       mumps_io_err.c
		       mumps_pord.c
		       mumps_size.c
		       mumps_c.c
		       mumps_io_thread.c
		       mumps_register_thread.c
		       mumps_thread.c
		       mumps_common.c
		       mumps_metis.c
		       mumps_save_restore_C.c
		       mumps_thread_affinity.c
		       mumps_config_file_C.c
		       mumps_metis64.c
		       mumps_scotch.c
		       mumps_io.c
		       mumps_metis_int.c
		       mumps_scotch64.c
)

add_library(MUMPS_F)

target_sources(MUMPS_F PUBLIC
		       ana_AMDMF.F
		       ana_blk.F		
		       ana_blk_m.F	       
		       ana_omp_m.F
		       ana_orderings.F
		       ana_orderings_wrappers_m.F
		       ana_set_ordering.F
		       bcast_errors.F
		       dana_LDLT_preprocess.F
		       dana_aux.F
		       dana_driver.F
		       dana_aux_ELT.F
		       dana_aux_par.F
		       dana_dist_m.F
		       dana_lr.F
		       dana_mtrans.F
		       dana_reordertree.F
		       darrowheads.F
		       dbcast_int.F
		       dend_driver.F
		       dfac_asm.F
		       dfac_asm_ELT.F
		       dfac_asm_master_ELT_m.F
		       dfac_asm_master_m.F
		       dfac_b.F
		       dfac_determinant.F
		       dfac_distrib_ELT.F
		       dfac_distrib_distentry.F
		       dfac_driver.F
		       dfac_front_LDLT_type1.F
		       dfac_front_LDLT_type2.F
		       dfac_front_LU_type1.F
		       dfac_front_LU_type2.F
		       dfac_front_aux.F
		       dfac_front_type2_aux.F
		       dfac_lastrtnelind.F
		       dfac_lr.F
		       dfac_mem_alloc_cb.F
		       dfac_mem_compress_cb.F
		       dfac_mem_dynamic.F
		       dfac_mem_free_block_cb.F
		       dfac_mem_stack.F
		       dfac_mem_stack_aux.F
		       dfac_omp_m.F
		       dfac_par_m.F
		       dfac_process_band.F
		       dfac_process_bf.F
		       dfac_process_blfac_slave.F
		       dfac_process_blocfacto.F
		       dfac_process_blocfacto_LDLT.F
		       dfac_process_contrib_type1.F
		       dfac_process_contrib_type2.F
		       dfac_process_contrib_type3.F
		       dfac_process_end_facto_slave.F
		       dfac_process_maprow.F
		       dfac_process_master2.F
		       dfac_process_message.F
		       dfac_process_root2slave.F
		       dfac_process_root2son.F
		       dfac_process_rtnelind.F
		       dfac_root_parallel.F
		       dfac_scalings.F
		       dfac_scalings_simScaleAbs.F
		       dfac_scalings_simScale_util.F
		       dfac_sispointers_m.F
		       dfac_sol_l0omp_m.F
		       dfac_sol_pool.F
		       dfac_type3_symmetrize.F
		       dini_defaults.F
		       dini_driver.F
		       dlr_core.F
		       dlr_stats.F
		       dlr_type.F
		       dmumps_comm_buffer.F
		       dmumps_config_file.F
		       dmumps_driver.F
		       dmumps_f77.F
		       dmumps_iXamax.F
		       dmumps_load.F
		       dmumps_lr_data_m.F
		       dmumps_ooc.F
		       dmumps_ooc_buffer.F
		       dmumps_save_restore.F
		       dmumps_save_restore_files.F
		       dmumps_sol_es.F
		       dmumps_struc_def.F
		       domp_tps_m.F
		       dooc_panel_piv.F
		       double_linked_list.F
		       drank_revealing.F
		       dsol_aux.F
		       dsol_bwd.F
		       dsol_bwd_aux.F
		       dsol_c.F
		       dsol_distrhs.F
		       dsol_driver.F
		       dsol_fwd.F
		       dsol_fwd_aux.F
		       dsol_lr.F
		       dsol_matvec.F
		       dsol_omp_m.F
		       dsol_root_parallel.F
		       dstatic_ptr_m.F
		       dtools.F
		       dtype3_root.F
		       estim_flops.F
		       fac_asm_build_sort_index_ELT_m.F
		       fac_asm_build_sort_index_m.F
		       fac_descband_data_m.F
		       fac_future_niv2_mod.F
		       fac_ibct_data_m.F
		       fac_maprow_data_m.F
		       front_data_mgt_m.F
		       lr_common.F
		       mumps_comm_ibcast.F
		       mumps_l0_omp_m.F
		       mumps_memory_mod.F
		       mumps_mpitoomp_m.F
		       mumps_ooc_common.F
		       mumps_print_defined.F
		       mumps_static_mapping.F
		       mumps_type2_blocking.F
		       mumps_type_size.F
		       mumps_version.F
		       omp_tps_common_m.F
		       sol_common.F
		       tools_common.F
)

target_include_directories(MUMPS_F PRIVATE ${MPI_CXX_INCLUDE_PATH})  

set_target_properties(MUMPS_F PROPERTIES LINKER_LANGUAGE CXX)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
 MESSAGE("Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
 target_compile_options(MUMPS_C PUBLIC -DAdd_ -Dpord)
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel") 
 MESSAGE("Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
 MESSAGE("MSVC")
else()
 MESSAFE("UKNOWN COMPILER ${CMAKE_CXX_COMPILER_ID}")
endif()


include_directories(${CMAKE_CURRENT_LIST_DIR})


