
include ../../Makefile.def

OBJS       = DL_Interpreter.o \
	TclInterpreter.o \
	PythonInterpreter.o \
	OpenSeesCommandsTcl.o \
	OpenSeesCommandsPython.o \
	OpenSeesUniaxialMaterialCommands.o \
	PythonModelBuilder.o PythonAnalysisBuilder.o

OBJS  = DL_Interpreter.o OpenSeesCommands.o OpenSeesUniaxialMaterialCommands.o OpenSeesElementCommands.o OpenSeesTimeSeriesCommands.o OpenSeesPatternCommands.o OpenSeesSectionCommands.o OpenSeesOutputCommands.o OpenSeesCrdTransfCommands.o OpenSeesBeamIntegrationCommands.o OpenSeesNDMaterialCommands.o OpenSeesMiscCommands.o OpenSeesParameterCommands.o OpenSeesFrictionModelCommands.o OpenSeesReliabilityCommands.o

TclOBJS =  tclMain.o TclInterpreter.o TclWrapper.o $(OBJS)
PythonOBJS = pythonMain.o PythonInterpreter.o PythonWrapper.o $(OBJS)
PythonModuleOBJS = PythonModule.o PythonWrapper.o $(OBJS)

# Compilation control

all:  pythonmodule

tcl: $(TclOBJS)
	$(LINKER) $(LINKFLAGS) $(TclOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	 -o tclInterpreter

#PYTHON_LIBRARY = -framework Python
PYTHON_LIBRARY =

python:  $(PythonOBJS)
	$(LINKER) $(LINKFLAGS) $(PythonOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS)  \
	 -o pythonInterpreter

#	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace pythonExample.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o fmkSum.dylib

OBJSm = OpenSeesCommandsPython.o OpenSeesUniaxialMaterialCommands.o PythonModelBuilder.o PythonAnalysisBuilder.o

LINKFLAGS = -Wl,-rpath,\$$ORIGIN/lib

pythonmodule: $(PythonModuleOBJS)
	$(RM) $(RMFLAGS) opensees.so;
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(MAKE); \
		$(CD) ..; \
	done );
	$(CD) $(FE)/interpreter;
	$(LINKER) $(LINKFLAGS) -shared $(PythonModuleOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(PARALLEL_LIB) \
	$(INTERPRETER_LIBS_MPI_PARAMETERS) \
	 -o opensees.so


INTERPRETER_LIBS_MPI_PARAMETERS =

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)
INTERPRETER_LIBS_MPI_PARAMETERS = \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE)/actor/machineBroker/MPI_MachineBroker.o \
	$(FE)/actor/machineBroker/MachineBroker.o \
	$(FE)/actor/actor/Actor.o \
	$(FE)/actor/shadow/Shadow.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.o \
	$(FE)/domain/subdomain/ActorSubdomain.o \
	$(FE)/domain/subdomain/ShadowSubdomain.o
endif


# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
