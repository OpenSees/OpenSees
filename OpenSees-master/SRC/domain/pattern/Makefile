include ../../../Makefile.def

ifeq ($(RELIABILITY), YES_RELIABILITY)

RELIABILITY_OBJS = DiscretizedRandomProcessSeries.o \
	SimulatedRandomProcessSeries.o

else

RELIABILITY_OBJS =

endif

OBJS       = $(RELIABILITY_OBJS) \
	EarthquakePattern.o \
	LinearSeries.o \
	LoadPattern.o \
	FireLoadPattern.o \
	LoadPatternIter.o \
	PathSeries.o \
	PathTimeSeries.o \
	PathTimeSeriesThermal.o \
	RectangularSeries.o \
	TimeSeries.o \
	TclPatternCommand.o \
	TrigSeries.o \
	ConstantSeries.o \
	UniformExcitation.o \
	TimeSeriesIntegrator.o \
	TrapezoidalTimeSeriesIntegrator.o \
	SimpsonTimeSeriesIntegrator.o \
	MultiSupportPattern.o \
	TclSeriesCommand.o \
	PulseSeries.o \
	TriangleSeries.o \
	TclSeriesIntegratorCommand.o \
	PeerMotion.o \
	PeerNGAMotion.o

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/domain/pattern/drm; $(MAKE);


# Miscellaneous
tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/domain/pattern/drm; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.
