include ../../../Makefile.def

OBJS       = EarthquakePattern.o PBowlLoading.o LinearSeries.o \
	LoadPattern.o LoadPatternIter.o PathSeries.o PathTimeSeries.o \
	RectangularSeries.o TimeSeries.o TclPatternCommand.o \
	TrigSeries.o ConstantSeries.o UniformExcitation.o \
	TimeSeriesIntegrator.o TrapezoidalTimeSeriesIntegrator.o \
	DiscretizedRandomProcessSeries.o SimulatedRandomProcessSeries.o \
	MultiSupportPattern.o TclSeriesCommand.o TclSeriesIntegratorCommand.o \
	PulseSeries.o TriangleSeries.o

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
