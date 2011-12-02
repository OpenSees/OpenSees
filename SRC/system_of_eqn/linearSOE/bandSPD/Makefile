include ../../../../Makefile.def

OBJS       = BandSPDLinSOE.o BandSPDLinSolver.o BandSPDLinLapackSolver.o 

PROGRAM = go

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM)

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
