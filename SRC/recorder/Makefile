include ../../Makefile.def

OBJS       = Recorder.o \
	DatastoreRecorder.o \
	ElementRecorder.o \
	NodeRecorder.o \
	EnvelopeElementRecorder.o \
	EnvelopeNodeRecorder.o \
	TclRecorderCommands.o \
	FilePlotter.o \
	AlgorithmIncrements.o \
	DriftRecorder.o \
	EnvelopeDriftRecorder.o \
	PatternRecorder.o \
	DamageRecorder.o

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/recorder/response; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/recorder/response; $(MAKE) wipe;

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
