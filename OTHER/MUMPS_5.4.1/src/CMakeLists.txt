
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
		       ana_AMDMF.f
		       ana_blk.f		
		       ana_blk_m.f	       
		       ana_omp_m.f
		       ana_orderings.f
		       ana_orderings_wrappers_m.f
		       ana_set_ordering.f
		       bcast_errors.f
		       dana_LDLT_preprocess.f
		       dana_aux.f
		       dana_driver.f
		       dana_aux_ELT.f
		       dana_aux_par.f
		       dana_dist_m.f
		       dana_lr.f
		       dana_mtrans.f
		       dana_reordertree.f
		       darrowheads.f
		       dbcast_int.f
		       dend_driver.f
		       dfac_asm.f
		       dfac_asm_ELT.f
		       dfac_asm_master_ELT_m.f
		       dfac_asm_master_m.f
		       dfac_b.f
		       dfac_determinant.f
		       dfac_distrib_ELT.f
		       dfac_distrib_distentry.f
		       dfac_driver.f
		       dfac_front_LDLT_type1.f
		       dfac_front_LDLT_type2.f
		       dfac_front_LU_type1.f
		       dfac_front_LU_type2.f
		       dfac_front_aux.f
		       dfac_front_type2_aux.f
		       dfac_lastrtnelind.f
		       dfac_lr.f
		       dfac_mem_alloc_cb.f
		       dfac_mem_compress_cb.f
		       dfac_mem_dynamic.f
		       dfac_mem_free_block_cb.f
		       dfac_mem_stack.f
		       dfac_mem_stack_aux.f
		       dfac_omp_m.f
		       dfac_par_m.f
		       dfac_process_band.f
		       dfac_process_bf.f
		       dfac_process_blfac_slave.f
		       dfac_process_blocfacto.f
		       dfac_process_blocfacto_LDLT.f
		       dfac_process_contrib_type1.f
		       dfac_process_contrib_type2.f
		       dfac_process_contrib_type3.f
		       dfac_process_end_facto_slave.f
		       dfac_process_maprow.f
		       dfac_process_master2.f
		       dfac_process_message.f
		       dfac_process_root2slave.f
		       dfac_process_root2son.f
		       dfac_process_rtnelind.f
		       dfac_root_parallel.f
		       dfac_scalings.f
		       dfac_scalings_simScaleAbs.f
		       dfac_scalings_simScale_util.f
		       dfac_sispointers_m.f
		       dfac_sol_l0omp_m.f
		       dfac_sol_pool.f
		       dfac_type3_symmetrize.f
		       dini_defaults.f
		       dini_driver.f
		       dlr_core.f
		       dlr_stats.f
		       dlr_type.f
		       dmumps_comm_buffer.f
		       dmumps_config_file.f
		       dmumps_driver.f
		       dmumps_f77.f
		       dmumps_iXamax.f
		       dmumps_load.f
		       dmumps_lr_data_m.f
		       dmumps_ooc.f
		       dmumps_ooc_buffer.f
		       dmumps_save_restore.f
		       dmumps_save_restore_files.f
		       dmumps_sol_es.f
		       dmumps_struc_def.f
		       domp_tps_m.f
		       dooc_panel_piv.f
		       double_linked_list.f
		       drank_revealing.f
		       dsol_aux.f
		       dsol_bwd.f
		       dsol_bwd_aux.f
		       dsol_c.f
		       dsol_distrhs.f
		       dsol_driver.f
		       dsol_fwd.f
		       dsol_fwd_aux.f
		       dsol_lr.f
		       dsol_matvec.f
		       dsol_omp_m.f
		       dsol_root_parallel.f
		       dstatic_ptr_m.f
		       dtools.f
		       dtype3_root.f
		       estim_flops.f
		       fac_asm_build_sort_index_ELT_m.f
		       fac_asm_build_sort_index_m.f
		       fac_descband_data_m.f
		       fac_future_niv2_mod.f
		       fac_ibct_data_m.f
		       fac_maprow_data_m.f
		       front_data_mgt_m.f
		       lr_common.f
		       mumps_comm_ibcast.f
		       mumps_l0_omp_m.f
		       mumps_memory_mod.f
		       mumps_mpitoomp_m.f
		       mumps_ooc_common.f
		       mumps_print_defined.f
		       mumps_static_mapping.f
		       mumps_type2_blocking.f
		       mumps_type_size.f
		       mumps_version.f
		       omp_tps_common_m.f
		       sol_common.f
		       tools_common.f
)

target_include_directories(MUMPS_F PRIVATE ${MPI_CXX_INCLUDE_PATH})  

set_target_properties(MUMPS_F PROPERTIES LINKER_LANGUAGE CXX)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
 MESSAGE("Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
 target_compile_options(MUMPS_C PUBLIC -DAdd_ -Dpord)
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel") 
 MESSAGE("Intel")
 target_compile_options(MUMPS_F PUBLIC -Dnounderscore)
 target_compile_options(MUMPS_C PUBLIC -DAdd_ -Dpord)
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
 MESSAGE("MSVC")
else()
 MESSAFE("UKNOWN COMPILER ${CMAKE_CXX_COMPILER_ID}")
endif()


include_directories(${CMAKE_CURRENT_LIST_DIR})


