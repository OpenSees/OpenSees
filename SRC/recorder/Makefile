include ../../Makefile.def

OBJS       = Recorder.o DatastoreRecorder.o MaxNodeDispRecorder.o \
	ElementRecorder.o TclRecorderCommands.o FilePlotter.o \
	NodeRecorder.o AlgorithmIncrements.o DriftRecorder.o \
	EnvelopeNodeRecorder.o  EnvelopeElementRecorder.o GSA_Recorder.o \
	YsVisual.o PatternRecorder.o DamageRecorder.o

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
