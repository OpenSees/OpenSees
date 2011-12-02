include ../../../Makefile.def

ifeq ($(RELIABILITY), YES_RELIABILITY)

RELIABILITY_OBJS = DiscretizedRandomProcessSeries.o \
	SimulatedRandomProcessSeries.o

else

RELIABILITY_OBJS =

endif

OBJS       = $(RELIABILITY_OBJS) \
	EarthquakePattern.o \
	PBowlLoading.o \
	LinearSeries.o \
	LoadPattern.o \
	LoadPatternIter.o \
	PathSeries.o \
	PathTimeSeries.o \
	RectangularSeries.o \
	TimeSeries.o \
	TclPatternCommand.o \
	TrigSeries.o \
	ConstantSeries.o \
	UniformExcitation.o \
	TimeSeriesIntegrator.o \
	TrapezoidalTimeSeriesIntegrator.o \
	MultiSupportPattern.o \
	TclSeriesCommand.o \
	PulseSeries.o \
	TriangleSeries.o \
	TclSeriesIntegratorCommand.o \
	PeerMotion.o \
	PeerNGAMotion.o

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
