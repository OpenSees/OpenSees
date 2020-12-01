include ../../Makefile.def

ifeq ($(GRAPHICS), NONE)
GRAPHIC_OBJECTS = 
else
GRAPHIC_OBJECTS = FilePlotter.o AlgorithmIncrements.o 
endif

OBJS       = Recorder.o \
	DatastoreRecorder.o \
	ElementRecorder.o \
	ElementRecorderRMS.o \
	NodeRecorder.o \
	NodeRecorderRMS.o \
	EnvelopeElementRecorder.o \
	NormElementRecorder.o \
	NormEnvelopeElementRecorder.o \
	EnvelopeNodeRecorder.o \
	TclRecorderCommands.o \
	DriftRecorder.o \
	EnvelopeDriftRecorder.o \
	PatternRecorder.o \
	RemoveRecorder.o \
	DamageRecorder.o $(GRAPHIC_OBJECTS) \
	PVDRecorder.o MPCORecorder.o GmshRecorder.o \
	VTK_Recorder.o


# Compilation control

all:      $(OBJS) $(GRAPHIC_OBJECTS)
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
