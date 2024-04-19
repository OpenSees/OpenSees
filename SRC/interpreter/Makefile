
include ../../Makefile.def

OBJS  = DL_Interpreter.o OpenSeesCommands.o OpenSeesUniaxialMaterialCommands.o OpenSeesElementCommands.o OpenSeesTimeSeriesCommands.o OpenSeesPatternCommands.o OpenSeesSectionCommands.o OpenSeesOutputCommands.o OpenSeesCrdTransfCommands.o OpenSeesDampingCommands.o OpenSeesBeamIntegrationCommands.o OpenSeesNDMaterialCommands.o OpenSeesMiscCommands.o OpenSeesParameterCommands.o OpenSeesFrictionModelCommands.o OpenSeesReliabilityCommands.o OpenSeesNDTestCommands.o OpenSeesIGACommands.o 

PythonOtherFiles = ../reliability/domain/functionEvaluator/PythonEvaluator.o

TclOBJS =  tclMain.o TclInterpreter.o TclWrapper.o $(OBJS)
PythonOBJS = pythonMain.o PythonModule.o PythonWrapper.o PythonStream.o $(OBJS)
PythonModuleOBJS = PythonModule.o PythonWrapper.o PythonStream.o $(PythonOtherFiles) $(OBJS)

# Compilation control

all:  pythonmodule

tcl: $(TclOBJS)
	$(LINKER) $(LINKFLAGS) $(TclOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(PARALLEL_LIB) \
	$(INTERPRETER_LIBS_MPI_PARAMETERS) \
	 -o openseestcl

#PYTHON_LIBRARY = -framework Python

python:  $(PythonOBJS)
	$(RM) $(RMFLAGS) openseespy;
	$(LINKER) $(LINKFLAGS) -Wl,-rpath,\$$ORIGIN/lib $(PythonOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(PARALLEL_LIB) \
	$(INTERPRETER_LIBS_MPI_PARAMETERS) \
	 -o openseespy

#	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace pythonExample.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o fmkSum.dylib

OBJSm = OpenSeesCommandsPython.o OpenSeesUniaxialMaterialCommands.o PythonModelBuilder.o PythonAnalysisBuilder.o

pythonmodule: $(PythonModuleOBJS)
	$(RM) $(RMFLAGS) opensees.so;
	$(LINKER) $(LINKFLAGS) -Wl,-rpath,\$$ORIGIN/lib -shared $(PythonModuleOBJS) \
	$(INTERPRETER_LIBS_MPI_PARAMETERS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(PARALLEL_LIB) \
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
