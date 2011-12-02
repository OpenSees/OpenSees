# makefile for sparse supernodal LU, implemented in ANSI C
include ../../Makefile.def

ALLAUX = superlu_timer.o lsame.o util.o memory.o get_perm_c.o mmd.o \
	 sp_coletree.o sp_preorder.o sp_ienv.o relax_snode.o \
	 xerbla.o

DZLAUX = dlamch.o

DLUSRC = \
	dgssv.o dgssvx.o \
	dsp_blas2.o dsp_blas3.o dgscon.o dlacon.o \
	dlangs.o dgsequ.o dlaqgs.o dpivotgrowth.o  \
	dgsrfs.o dgstrf.o dgstrs.o dcopy_to_ucol.o \
	dsnode_dfs.o dsnode_bmod.o \
	dpanel_dfs.o dpanel_bmod.o dreadhb.o \
	dcolumn_dfs.o dcolumn_bmod.o dpivotL.o dpruneL.o \
	dmemory.o dutil.o


all:    double

double: $(DLUSRC) $(ALLAUX) $(DZLAUX)
	@$(RM) $(RMFLAGS) $(SUPERLU_LIBRARY)
	@$(AR) $(ARFLAGS) $(SUPERLU_LIBRARY) $(DLUSRC) $(ALLAUX) $(DZLAUX)
	@$(RANLIB) $(SUPERLU_LIBRARY)

##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c $<

superlu_timer.o:  superlu_timer.c ; $(CC) -c $<

clean:	
	@$(RM) $(RMFLAGS) *.o

wipe: clean

