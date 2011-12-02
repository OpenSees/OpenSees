include ../../../Makefile.def

OBJS       = ForceBeamColumn2d.o ForceBeamColumn3d.o \
	ElasticForceBeamColumn2d.o ElasticForceBeamColumn3d.o \
	TclForceBeamColumnCommand.o \
	BeamIntegration.o \
	LobattoBeamIntegration.o \
	LegendreBeamIntegration.o \
	RadauBeamIntegration.o \
	NewtonCotesBeamIntegration.o \
	HingeMidpointBeamIntegration.o \
	HingeEndpointBeamIntegration.o \
	HingeRadauBeamIntegration.o \
	HingeRadauTwoBeamIntegration.o \
	UserDefinedBeamIntegration.o \
	UserDefinedHingeIntegration.o \
	DistHingeIntegration.o \
	TrapezoidalBeamIntegration.o \
	FixedLocationBeamIntegration.o \
	LowOrderBeamIntegration.o \
	MidDistanceBeamIntegration.o \
	ElasticForceBeamColumn2d.o \
	ElasticForceBeamColumn3d.o

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









