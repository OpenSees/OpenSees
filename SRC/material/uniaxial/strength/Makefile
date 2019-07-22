include ../../../../Makefile.def

OBJS       = 	StrengthDegradation.o \
	ACIStrengthDegradation.o ConstantStrengthDegradation.o \
	DuctilityStrengthDegradation.o EnergyStrengthDegradation.o \
	PetrangeliStrengthDegradation.o SectionStrengthDegradation.o \
	TclModelBuilderStrengthDegradationCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
