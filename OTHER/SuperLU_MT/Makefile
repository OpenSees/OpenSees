# makefile for sparse supernodal LU, implemented in ANSI C, targeted for 
# shared memory machines (SMPs).
#

include ../../Makefile.def

#######################################################################
#  This is the makefile to create a library for SuperLU.
#  The files are organized as follows:
#
#       ALLAUX -- Auxiliary routines called from all precisions
#       DLUSRC -- Double precision real SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
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

ALLAUX = superlu_timer.o dclock.o sp_ienv.o lsame.o xerbla.o \
	util.o pmemory.o qrnzcnt.o dmyblas2.o await.o \
	get_perm_c.o mmd.o colamd.o sp_coletree.o sp_colorder.o \
	pxgstrf_relax_snode.o \
	pxgstrf_synch.o \
	pxgstrf_scheduler.o \
	pxgstrf_mark_busy_descends.o \
	pxgstrf_pruneL.o \
	pxgstrf_super_bnd_dfs.o

DZAUX  = dlamch.o

DLUSRC = dreadhb.o dmatgen.o \
	pdgssv.o pdgssvx.o dgstrs.o dgsrfs.o \
	dgscon.o dlacon.o dlangs.o dgsequ.o dlaqgs.o dpivotgrowth.o \
	pdmemory.o pdutil.o \
	pdgstrf.o pdgstrf_init.o pdgstrf_finalize.o \
	pdgstrf_thread.o pdgstrf_thread_init.o pdgstrf_thread_finalize.o \
	pdgstrf_factor_snode.o pdgstrf_snode_dfs.o pdgstrf_snode_bmod.o \
	pdgstrf_panel_dfs.o pdgstrf_panel_bmod.o \
	pdgstrf_copy_to_ucol.o pdgstrf_pivotL.o \
	pdgstrf_column_dfs.o pdgstrf_column_bmod.o \
	pdgstrf_bmod1D.o pdgstrf_bmod2D.o \
	pdgstrf_bmod1D_mv2.o pdgstrf_bmod2D_mv2.o \
	dsp_blas2.o dsp_blas3.o

double: $(DLUSRC) $(ALLAUX) $(DZAUX)
	@$(RM) $(RMFLAGS) $(SUPERLU_LIBRARY)
	@$(AR) $(ARFLAGS) $(SUPERLU_LIBRARY) $(ALLAUX) $(DLUSRC) $(DZAUX)
	@$(RANLIB) $(SUPERLU_LIBRARY)


##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c  $<
slamch.o: slamch.c ; $(CC) -c  $<
timer.o:  timer.c ;  $(CC) -c  $<

# Do not optimize this on Cray; related to 'volatile' variable.
await.o: await.c
	$(CC) -c $(NOOPTS) $< $(VERBOSE)

pdgstrf_bmod1D.o: pdgstrf_bmod1D.c
	$(CC) $(CFLAGS) $(BLASOPTS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

pdgstrf_bmod2D.o: pdgstrf_bmod2D.c
	$(CC) $(CFLAGS) $(BLASOPTS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

dmyblas2.o: dmyblas2.c
	$(CC) $(CFLAGS) $(BLASOPTS) $(BLASDEF) -c $< $(VERBOSE)

dgemv2.o: dgemv2.f
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

#.c.o:
#	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

clean:	
	rm -f *.o core ../$(SUPERLULIB)
