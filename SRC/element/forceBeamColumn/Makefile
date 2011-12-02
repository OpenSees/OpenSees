include ../../../Makefile.def

OBJS       = ForceBeamColumn2d.o ForceBeamColumn3d.o \
	TclForceBeamColumnCommand.o \
	BeamIntegration.o \
	LobattoBeamIntegration.o LegendreBeamIntegration.o \
	NewtonCotesBeamIntegration.o \
	HingeMidpointBeamIntegration2d.o HingeMidpointBeamIntegration3d.o \
	HingeMidpointBeamIntegration.o \
	HingeEndpointBeamIntegration2d.o HingeEndpointBeamIntegration3d.o \
	HingeEndpointBeamIntegration.o \
	HingeRadauBeamIntegration2d.o HingeRadauBeamIntegration3d.o \
	HingeRadauBeamIntegration.o \
	HingeRadauTwoBeamIntegration2d.o HingeRadauTwoBeamIntegration3d.o \
	HingeRadauTwoBeamIntegration.o \
	HingeModRadauTwoBeamIntegration2d.o \
	UserDefinedBeamIntegration.o \
	UserDefinedHingeIntegration2d.o UserDefinedHingeIntegration3d.o \

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









