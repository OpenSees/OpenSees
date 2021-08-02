include ../../../Makefile.def

OBJS       = ForceBeamColumn2d.o ForceBeamColumn3d.o \
	ForceBeamColumn2dThermal.o \
	ElasticForceBeamColumn2d.o ElasticForceBeamColumn3d.o \
	ElasticForceBeamColumnWarping2d.o \
	ForceBeamColumnCBDI2d.o ForceBeamColumnCBDI3d.o \
	ForceBeamColumnWarping2d.o \
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
	SimpsonBeamIntegration.o \
	UserDefinedHingeIntegration.o \
	DistHingeIntegration.o \
	RegularizedHingeIntegration.o \
	TrapezoidalBeamIntegration.o \
	CompositeSimpsonBeamIntegration.o \
	FixedLocationBeamIntegration.o \
	LowOrderBeamIntegration.o \
	MidDistanceBeamIntegration.o \
	GaussQBeamIntegration.o \
	ChebyshevBeamIntegration.o \
	gaussq.o d1mach.o

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









