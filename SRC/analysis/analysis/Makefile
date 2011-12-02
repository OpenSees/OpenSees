include ../../../Makefile.def

OBJS       = DomainUser.o Analysis.o StaticAnalysis.o TransientAnalysis.o \
	     DirectIntegrationAnalysis.o DomainDecompositionAnalysis.o \
	     SubstructuringAnalysis.o EigenAnalysis.o \
	     VariableTimeStepDirectIntegrationAnalysis.o \
	     StaticDomainDecompositionAnalysis.o \
	     TransientDomainDecompositionAnalysis.o

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
