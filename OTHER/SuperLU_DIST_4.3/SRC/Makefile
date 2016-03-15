
include ../../../Makefile.def
DSUPERLULIB = $(DISTRIBUTED_SUPERLU_LIBRARY)
ARCH = $(AR)
ARCHFLAGS = $(ARFLAGS)

#######################################################################
#
#  This makefile creates a library for distributed SuperLU.
#  The files are organized as follows:
#
#       ALLAUX  -- Auxiliary routines called from all precisions
#       DSLUSRC -- Double precision real serial SuperLU routines
#       DPLUSRC -- Double precision real parallel SuperLU routines
#       ZSLUSRC -- Double precision complex serial SuperLU routines
#       ZPLUSRC -- Double precision complex parallel SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the two precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make double
#       make double complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all two precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################
# include ../make.inc
#
# Precision independent routines
#
ALLAUX 	= sp_ienv.o etree.o sp_colorder.o get_perm_c.o \
	  mmd.o comm.o memory.o util.o superlu_grid.o \
	  pxerr_dist.o superlu_timer.o GetDiagU.o symbfact.o \
	  psymbfact.o psymbfact_util.o mc64ad_dist.o \
	  static_schedule.o xerr_dist.o smach_dist.o dmach_dist.o

ifeq "${ACC}" "GPU"
ALLAUX += cublas_utils.o
endif

#
# Routines literally taken from SuperLU, but renamed with suffix _dist
#
DSLUSRC	= dlangs_dist.o dgsequ_dist.o dlaqgs_dist.o dutil_dist.o \
	  dmemory_dist.o dmyblas2_dist.o dsp_blas2_dist.o dsp_blas3_dist.o
ZSLUSRC	= dcomplex_dist.o zlangs_dist.o zgsequ_dist.o zlaqgs_dist.o \
	  zutil_dist.o zmemory_dist.o zmyblas2_dist.o \
	  zsp_blas2_dist.o zsp_blas3_dist.o

#
# Routines for double precision parallel SuperLU
DPLUSRC = pdgssvx.o pdgssvx_ABglobal.o \
	  dreadhb.o dreadrb.o dreadtriple.o dreadMM.o \
	  pdgsequ.o pdlaqgs.o dldperm_dist.o pdlangs.o pdutil.o \
	  pdsymbfact_distdata.o ddistribute.o pddistribute.o \
	  pdgstrf.o pdgstrf2.o \
	  pdgstrs.o pdgstrs1.o pdgstrs_lsum.o pdgstrs_Bglobal.o \
	  pdgsrfs.o pdgsmv.o pdgsrfs_ABXglobal.o pdgsmv_AXglobal.o 

#
# Routines for double complex parallel SuperLU
ZPLUSRC = pzgssvx.o pzgssvx_ABglobal.o \
	  zreadhb.o zreadrb.o zreadtriple.o zreadMM.o \
	  pzgsequ.o pzlaqgs.o zldperm_dist.o pzlangs.o pzutil.o \
	  pzsymbfact_distdata.o zdistribute.o pzdistribute.o \
	  pzgstrf.o pzgstrf2.o \
	  pzgstrs.o pzgstrs1.o pzgstrs_lsum.o pzgstrs_Bglobal.o \
	  pzgsrfs.o pzgsmv.o pzgsrfs_ABXglobal.o pzgsmv_AXglobal.o 

all:  double 

double: $(DSLUSRC) $(DPLUSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) $(DSUPERLULIB) \
		$(DSLUSRC) $(DPLUSRC) $(ALLAUX)
	$(RANLIB) $(DSUPERLULIB)

complex16: $(ZSLUSRC) $(ZPLUSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) $(DSUPERLULIB) \
		$(ZSLUSRC) $(ZPLUSRC) $(ALLAUX)
	$(RANLIB) $(DSUPERLULIB)


.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o $(DSUPERLULIB)
