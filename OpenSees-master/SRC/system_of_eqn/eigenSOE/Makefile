include ../../../Makefile.def

OBJS       = EigenSOE.o \
	EigenSolver.o \
	ArpackSOE.o \
	ArpackSolver.o \
	SymBandEigenSOE.o \
	SymBandEigenSolver.o \
	FullGenEigenSOE.o \
	FullGenEigenSolver.o

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
