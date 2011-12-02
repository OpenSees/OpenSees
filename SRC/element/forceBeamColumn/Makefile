include ../../../Makefile.def

OBJS       = ForceBeamColumn2d.o ForceBeamColumn3d.o \
	TclForceBeamColumnCommand.o \
	BeamIntegration.o \
	LobattoBeamIntegration.o UserDefinedBeamIntegration.o \
	HingeMidpointBeamIntegration2d.o HingeMidpointBeamIntegration3d.o \
	HingeRadauBeamIntegration2d.o HingeRadauBeamIntegration3d.o \
	HingeRadauTwoBeamIntegration2d.o HingeRadauTwoBeamIntegration3d.o \
	UserDefinedHingeIntegration2d.o UserDefinedHingeIntegration3d.o

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









