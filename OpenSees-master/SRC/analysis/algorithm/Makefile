include ../../../Makefile.def

OBJS       = SolutionAlgorithm.o

# Compilation control

all:         $(OBJS)
	     @$(CD) $(FE)/analysis/algorithm/domainDecompAlgo; $(MAKE);
	     @$(CD) $(FE)/analysis/algorithm/equiSolnAlgo; $(MAKE);
	     @$(CD) $(FE)/analysis/algorithm/eigenAlgo; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	     @$(CD) $(FE)/analysis/algorithm/domainDecompAlgo; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/algorithm/equiSolnAlgo; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/algorithm/eigenAlgo; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.






