include ../../Makefile.def

I = -I. -I../AMD
C = $(CC) $(OPT_FLAG) -fexceptions -fPIC -DNCHOLMOD $(I)

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

# non-user-callable umf_*.[ch] files:
UMFCH = umf_assemble umf_blas3_update umf_build_tuples umf_create_element \
	umf_extend_front umf_garbage_collection umf_get_memory \
	umf_init_front umf_kernel umf_kernel_init umf_kernel_wrapup \
	umf_local_search umf_lsolve umf_ltsolve umf_mem_alloc_element \
	umf_mem_alloc_head_block umf_mem_alloc_tail_block \
	umf_mem_free_tail_block umf_mem_init_memoryspace \
	umf_report_vector umf_row_search umf_scale_column \
	umf_set_stats umf_solve umf_symbolic_usage umf_transpose \
	umf_tuple_lengths umf_usolve umf_utsolve umf_valid_numeric \
	umf_valid_symbolic umf_grow_front umf_start_front \
	umf_store_lu umf_scale

# non-user-callable umf_*.[ch] files, int/SuiteSparse_long versions only
# (no real/complex):
UMFINT = umf_analyze umf_apply_order umf_colamd umf_free umf_fsize \
	umf_is_permutation umf_malloc umf_realloc umf_report_perm \
	umf_singletons umf_cholmod

# non-user-callable, created from umf_ltsolve.c, umf_utsolve.c,
# umf_triplet.c, and umf_assemble.c , with int/SuiteSparse_long
# and real/complex versions:
UMF_CREATED = umf_lhsolve umf_uhsolve umf_triplet_map_nox \
	umf_triplet_nomap_x umf_triplet_nomap_nox umf_triplet_map_x \
	umf_assemble_fixq umf_store_lu_drop

# non-user-callable, int/SuiteSparse_long and real/complex versions:
UMF = $(UMF_CREATED) $(UMFCH)

# user-callable umfpack_*.[ch] files (int/SuiteSparse_long and real/complex):
UMFPACK = umfpack_col_to_triplet umfpack_defaults umfpack_free_numeric \
	umfpack_free_symbolic umfpack_get_numeric umfpack_get_lunz \
	umfpack_get_symbolic umfpack_get_determinant umfpack_numeric \
	umfpack_qsymbolic umfpack_report_control umfpack_report_info \
	umfpack_report_matrix umfpack_report_numeric umfpack_report_perm \
	umfpack_report_status umfpack_report_symbolic umfpack_report_triplet \
	umfpack_report_vector umfpack_solve umfpack_symbolic \
	umfpack_transpose umfpack_triplet_to_col umfpack_scale \
	umfpack_load_numeric umfpack_save_numeric \
	umfpack_load_symbolic umfpack_save_symbolic

# user-callable, created from umfpack_solve.c (umfpack_wsolve.h exists, though):
# with int/SuiteSparse_long and real/complex versions:
UMFPACKW = umfpack_wsolve

USER = $(UMFPACKW) $(UMFPACK)

# user-callable, only one version for int/SuiteSparse_long,
# real/complex, *.[ch] files:
GENERIC = umfpack_timer umfpack_tictoc

#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

INC = umfpack.h ../AMD/SuiteSparse_config.h \
	umf_config.h umf_version.h \
	umf_internal.h umf_triplet.h \
	$(addprefix ./,  $(addsuffix .h,$(UMFCH))) \
	$(addprefix ./,  $(addsuffix .h,$(UMFINT))) \
	$(addprefix ./, $(addsuffix .h,$(USER))) \
	$(addprefix ./, $(addsuffix .h,$(GENERIC))) \
	../AMD/amd_internal.h ../AMD/amd.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------
UMFPACK_DI = $(addsuffix .o, $(subst umf_,umf_di_,$(UMF)) $(subst umfpack_,umfpack_di_,$(USER)))
UMFPACK_II = $(addsuffix .o, $(subst umf_,umf_i_,$(UMFINT)))
UMFPACK_GN = $(addsuffix .o, $(subst umfpack_,umfpack_gn_,$(GENERIC)))


# obj
UMFPACK_OBJ = $(UMFPACK_DI) $(UMFPACK_II) $(UMFPACK_GN)



all: umfpack

#-------------------------------------------------------------------------------
# compile each int and SuiteSparse_long routine (with no real/complex version)
#-------------------------------------------------------------------------------

umf_i_%.o: umf_%.c $(INC)
	$(C) -DDINT -c $< -o $@

#-------------------------------------------------------------------------------
# compile each routine in the DI version
#-------------------------------------------------------------------------------

umf_di_%.o: umf_%.c $(INC)
	$(C) -DDINT -c $< -o $@

umf_di_%hsolve.o: umf_%tsolve.c $(INC)
	$(C) -DDINT -DCONJUGATE_SOLVE -c $< -o $@

umf_di_triplet_map_x.o: umf_triplet.c $(INC)
	$(C) -DDINT -DDO_MAP -DDO_VALUES -c $< -o $@

umf_di_triplet_map_nox.o: umf_triplet.c $(INC)
	$(C) -DDINT -DDO_MAP -c $< -o $@

umf_di_triplet_nomap_x.o: umf_triplet.c $(INC)
	$(C) -DDINT -DDO_VALUES -c $< -o $@

umf_di_triplet_nomap_nox.o: umf_triplet.c $(INC)
	$(C) -DDINT -c $< -o $@

umf_di_assemble_fixq.o: umf_assemble.c $(INC)
	$(C) -DDINT -DFIXQ -c $< -o $@

umf_di_store_lu_drop.o: umf_store_lu.c $(INC)
	$(C) -DDINT -DDROP -c $< -o $@

umfpack_di_wsolve.o: umfpack_solve.c $(INC)
	$(C) -DDINT -DWSOLVE -c $< -o $@

umfpack_di_%.o: umfpack_%.c $(INC)
	$(C) -DDINT -c $< -o $@

#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

umfpack_gn_%.o: umfpack_%.c $(INC)
	$(C) -c $< -o $@

#-------------------------------------------------------------------------------
# Create the libumfpack.a library
#-------------------------------------------------------------------------------

umfpack: $(UMFPACK_OBJ)
	@$(RM) $(RMFLAGS) $(UMFPACK_LIBRARY)
	@$(AR) $(ARFLAGS) $(UMFPACK_LIBRARY) $(UMFPACK_OBJ)
	@$(RANLIB) $(UMFPACK_LIBRARY)

clean:
	@$(RM) $(RMFLAGS) *.o

wipe: clean

