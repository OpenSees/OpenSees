include ../../../../Makefile.def

OBJS       = 	CorrelationCoefficient.o CorrelationCoefficientIter.o \
		Cutset.o CutsetIter.o \
		RandomVariable.o RandomVariableIter.o \
		ReliabilityDomainComponent.o \
		ReliabilityDomain.o \
	PrincipalPlane.o

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
