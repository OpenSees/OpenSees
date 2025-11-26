include ../../../Makefile.def

OBJS       = SparsePythonCompressedEigenSOE.o \
	SparsePythonCompressedEigenSolver.o \
	SparsePythonCOOEigenSOE.o \
	SparsePythonCOOEigenSolver.o \
	SparsePythonEigenCommon.o \
	SparsePythonEigenFactory.o

all:    $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS)

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

