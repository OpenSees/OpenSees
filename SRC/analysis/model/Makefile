include ../../../Makefile.def

OBJS       = AnalysisModel.o  FE_EleIter.o DOF_GrpIter.o

# Compilation control

all:         $(OBJS)
	@$(CD) simple; $(MAKE); 

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) simple; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
