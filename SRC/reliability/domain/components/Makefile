include ../../../../Makefile.def

OBJS       = 	CorrelationCoefficient.o \
		LimitStateFunction.o \
		RandomVariable.o \
		RandomVariablePositioner.o \
		ReliabilityDomainComponent.o \
		ReliabilityDomain.o \
		ParameterPositioner.o

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
