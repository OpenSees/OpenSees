 
include ../../Makefile.def

DOBJ  = dgeqr2.o dlabad.o dlacon.o dlacpy.o dladiv.o dlae2.o dlaev2.o\
        dlaexc.o dlagtm.o dlahqr.o dlaln2.o dlamch.o dlange.o dlanhs.o\
        dlanst.o dlanv2.o dlaptm.o dlapy2.o dlapy3.o dlaran.o dlarf.o\
        dlarfg.o dlarfx.o dlarnd.o dlarnv.o dlartg.o dlaruv.o dlascl.o\
        dlaset.o dlasr.o  dlasrt.o dlassq.o dlasy2.o dorm2r.o dsteqr.o\
        dtrevc.o dtrexc.o dtrsen.o dtrsyl.o dgbtf2.o dgbtrf.o dgbtrs.o\
	dgttrf.o dgttrs.o dpttrf.o dpttrs.o dlaswp.o dgbsv.o dpbtrs.o\
	dgesv.o dgetrs.o dpbsv.o dgetrf.o dgetf2.o dpbtrf.o dpotf2.o\
	dpbtf2.o dgels.o dgeqrf.o dgelqf.o dormqr.o dormlq.o\
	dgelq2.o dlarft.o dlarfb.o dorml2.o\
	dlaebz.o dlagtf.o dlagts.o dlansb.o dlar2v.o dlargv.o\
	dlartv.o dsbevx.o dsbtrd.o dstebz.o dstein.o dsterf.o\
	dgetri.o dtrtri.o dtrti2.o


IOBJ  = ilaenv.o lsame.o lsamen.o xerbla.o xlaenv.o ieeeck.o


all: double 
 
double: $(DOBJ) $(IOBJ)
	@$(RM) $(RMFLAGS) $(LAPACK_LIBRARY)
	@$(AR) $(ARFLAGS) $(LAPACK_LIBRARY) $(DOBJ) $(IOBJ)
	@$(RANLIB) $(LAPACK_LIBRARY)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core go

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o test go

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
