include ../../../Makefile.def

OBJS       = ConstraintHandler.o PlainHandler.o \
	PenaltyConstraintHandler.o LagrangeConstraintHandler.o \
	TransformationConstraintHandler.o PenaltyHandlerNoHomoSPMultipliers.o

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
