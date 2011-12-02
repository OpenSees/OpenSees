include ../../../Makefile.def

OBJS       = EigenSOE.o EigenSolver.o BandArpackSOE.o  BandArpackSolver.o \
	     SymArpackSOE.o  SymArpackSolver.o

all:    $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
