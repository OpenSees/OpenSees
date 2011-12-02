############################################################################
#
#  Program:         SuperLU
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   October 2, 1995
#
#  Modified:        February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
# makefile for sparse supernodal LU, implemented in ANSI C

include ../../Makefile.def

ALLAUX = superlu_timer.o lsame.o util.o memory.o get_perm_c.o mmd.o \
	 sp_coletree.o sp_preorder.o sp_ienv.o relax_snode.o \
	 xerbla.o colamd.o

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
	@$(AR) $(ARFLAGS) $(SUPERLU_LIBRARY) $(ALLAUX) $(DLUSRC) $(DZLAUX)
	@$(RANLIB) $(SUPERLU_LIBRARY)

##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c $<

clean:	
	@$(RM) $(RMFLAGS) *.o

wipe: clean

