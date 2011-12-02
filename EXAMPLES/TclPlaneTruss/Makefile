include ../../Makefile.def

PROGRAM         = myOpenSees

all:         $(PROGRAM)

OBJS = TclPlaneTruss.o MyTruss.o myCommands.o

$(PROGRAM):  $(OBJS)
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(LINKER) $(LINKFLAGS) $(OBJS) $(INTERPRETER_LIBS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(TCL_LIBRARY) \
	-o $(PROGRAM)

INTERPRETER_LIBS = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/tclAppInit.o \
	$(FE)/tcl/commands.o \
	$(FE)/tcl/TclFeViewer.o \
	$(FE)/tcl/TclVideoPlayer.o 

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(PROGRAM)

# DO NOT DELETE THIS LINE -- make depend depends on it.
