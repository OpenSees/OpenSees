include ../../Makefile.def


OBJS       = commands.o TclFeViewer.o TclVideoPlayer.o


tcl: $(OBJS)  tclMain.o tclAppInit.o

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS       = commands.o 
tcl: $(OBJS) tclMain.o mpiMain.o

endif

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

tcl: $(OBJS) tclMain.o mpiParameterMain.o
OBJS       = commands.o 
endif

java: $(OBJS)  OpenSeesEvaluator.o OpenSeesEvaluatorImpl.o

all: $(OBJS)

tk: $(OBJS)  tkMain.o tkAppInit.o


FEAP_ELE_LIBS = $(FE)/element/feap/elmt01.o \
	$(FE)/element/feap/elmt02.o \
	$(FE)/element/feap/elmt03.o \
	$(FE)/element/feap/elmt04.o \
	$(FE)/element/feap/elmt05.o \
	$(FE)/element/feap/fillCommon.o \
	$(FE)/element/feap/getCommon.o \
	$(FE)/element/feap/dummyFeap.o \
	$(FE)/element/feap/fElement.o \
	$(FE)/element/feap/fElmt02.o 

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.





