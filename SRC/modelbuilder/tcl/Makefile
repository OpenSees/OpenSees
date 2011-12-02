# Makefile for model objects

include ../../../Makefile.def

OBJS       = TclModelBuilder.o myCommands.o 

# Compilation control

all: $(OBJS)  

tcl: myCommands.o
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(LINKER) $(LINKFLAGS) myCommands.o $(INTERPRETER_LIBS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

INTERPRETER_LIBS = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/tclAppInit.o \
	$(FE)/tcl/commands.o \
	$(FE)/tcl/TclFeViewer.o \
	$(FE)/tcl/TclVideoPlayer.o 

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o g3

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
