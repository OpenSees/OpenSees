include ../../../Makefile.def

OBJS       = ConstraintFunction.o \
	DesignVariable.o \
	DesignVariablePositioner.o \
	ObjectiveFunction.o \
	SnoptProblem.o \
	ToyOptFunction.o \
	ToyFunction.o \
	SnoptWrapper.o \
	SnoptAnalysis.o 


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
