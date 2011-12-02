include ../../Makefile.def

DOBJ  = idamax.o dasum.o daxpy.o dcopy.o ddot.o  dnrm2.o drot.o  drotg.o\
        dscal.o  dswap.o dgemv.o dger.o  dsymv.o dsyr2.o dtrmm.o dgbmv.o\
	dgemm.o dtrsm.o dtbsv.o dsyr.o dtrsv.o dsyrk.o dtrmv.o


all: double 
 
double: $(DOBJ) $(IOBJ)
	@$(RM) $(RMFLAGS) $(BLAS_LIBRARY)
	@$(AR) $(ARFLAGS) $(BLAS_LIBRARY) $(DOBJ) $(IOBJ)
	@$(RANLIB) $(BLAS_LIBRARY)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*#  core

clean: tidy
	@$(RM) $(RMFLAGS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
