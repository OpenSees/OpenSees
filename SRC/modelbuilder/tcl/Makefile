# Makefile for model objects

include ../../../Makefile.def

OBJS       = TclModelBuilder.o myCommands.o TclUniaxialMaterialTester.o \
	Block2D.o Block3D.o

# Compilation control

all: $(OBJS)  

tcl: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TCL) \
	$(GRAPHIC_LIBS_TCL) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

tk: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tk;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TK) \
	$(GRAPHIC_LIBS_TCL) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

ifeq ($(PROGRAMMING_MODE), PARALLEL)

tcl: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

tk: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tk;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TK) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

endif

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

tcl: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI_PARAMETERS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

tk: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tk;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TK) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

endif


INTERPRETER_LIBS_TCL = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/tclAppInit.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o 

INTERPRETER_LIBS_MPI = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/mpiMain.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE)/actor/machineBroker/MPI_MachineBroker.o \
	$(FE)/domain/subdomain/ActorSubdomain.o \
	$(FE)/domain/subdomain/ShadowSubdomain.o 

INTERPRETER_LIBS_MPI_PARAMETERS = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/mpiParameterMain.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE)/actor/machineBroker/MPI_MachineBroker.o \
	$(FE)/domain/subdomain/ActorSubdomain.o \
	$(FE)/domain/subdomain/ShadowSubdomain.o 


INTERPRETER_LIBS_TK = $(FE)/tcl/tkMain.o \
	$(FE)/tcl/tkAppInit.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o 

GRAPHIC_LIBS_TCL = $(FE)/tcl/TclFeViewer.o \
	$(FE)/tcl/TclVideoPlayer.o 


ifeq ($(GRAPHICS), NONE)

GRAPHIC_LIBS_TCL = 

endif

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o g3

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
