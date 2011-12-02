include ../../../../Makefile.def

OBJS       = BandSPDLinSOE.o \
	BandSPDLinSolver.o \
	BandSPDLinLapackSolver.o \
	DistributedBandSPDLinSOE.o

PROGRAM = go

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS)

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
