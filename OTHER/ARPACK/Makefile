############################################################################
#
#  Program:         ARPACK
#
#  Module:          Makefile
#
#  Purpose:         Sources Makefile
#
#  Creation date:   February 22, 1996
#
#  Modified:        September 9, 1996
#
#  Send bug reports, comments or suggestions to arpack.caam.rice.edu
#
############################################################################

include ../../Makefile.def

DOBJ  = dgetv0.o dlaqrb.o dstqrb.o dsortc.o dsortr.o dstatn.o dstats.o \
	dnaitr.o dnapps.o dnaup2.o dnaupd.o dnconv.o dneigh.o dngets.o \
	dsaitr.o dsapps.o dsaup2.o dsaupd.o dsconv.o dseigt.o dsgets.o \
	dneupd.o dseupd.o dsesrt.o dvout.o  dmout.o \
	icnteq.o icopy.o iset.o iswap.o ivout.o second.o

all: double

double: $(DOBJ)
	@$(RM) $(RMFLAGS) $(ARPACK_LIBRARY)
	@$(AR) $(ARFLAGS) $(ARPACK_LIBRARY) $(DOBJ)
	@$(RANLIB) $(ARPACK_LIBRARY)


clean:
	@$(RM) $(RMFLAGS) *.o a.out core

wipe: clean
