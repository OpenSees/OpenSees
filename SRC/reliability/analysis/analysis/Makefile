include ../../../../Makefile.def

OBJS       = 	FORMAnalysis.o \
	FOSMAnalysis.o \
	FragilityAnalysis.o \
	GFunVisualizationAnalysis.o \
	OutCrossingAnalysis.o \
	ParametricReliabilityAnalysis.o \
	SamplingAnalysis.o \
	ReliabilityAnalysis.o \
	SORMAnalysis.o \
	SystemAnalysis.o


# Compilation control
all:    $(OBJS)
	@$(CD) $(FE)/reliability/analysis/analysis/system; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/reliability/analysis/analysis/system; $(MAKE) wipe;

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
