include ../../../../Makefile.def

OBJS       = 	FORMAnalysis.o\
	FragilityAnalysis.o \
	GFunVisualizationAnalysis.o \
	MVFOSMAnalysis.o \
	OutCrossingAnalysis.o \
	SamplingAnalysis.o \
	ReliabilityAnalysis.o \
	SORMAnalysis.o \
	SystemAnalysis.o


# Compilation control
all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
