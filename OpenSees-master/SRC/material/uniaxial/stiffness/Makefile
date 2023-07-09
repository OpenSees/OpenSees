include ../../../../Makefile.def

OBJS       = 	StiffnessDegradation.o \
	ConstantStiffnessDegradation.o DuctilityStiffnessDegradation.o \
	EnergyStiffnessDegradation.o PincheiraStiffnessDegradation.o \
	TclModelBuilderStiffnessDegradationCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
