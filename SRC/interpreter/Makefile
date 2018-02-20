
-include ../../Makefile.def

OBJS       = DL_Interpreter.o \
	TclInterpreter.o \
	PythonInterpreter.o \
	OpenSeesCommandsTcl.o \
	OpenSeesCommandsPython.o \
	OpenSeesUniaxialMaterialCommands.o \
	PythonModelBuilder.o PythonAnalysisBuilder.o

OBJS  = DL_Interpreter.o OpenSeesCommands.o OpenSeesUniaxialMaterialCommands.o OpenSeesElementCommands.o OpenSeesTimeSeriesCommands.o OpenSeesPatternCommands.o OpenSeesSectionCommands.o OpenSeesOutputCommands.o OpenSeesCrdTransfCommands.o OpenSeesBeamIntegrationCommands.o OpenSeesNDMaterialCommands.o OpenSeesMiscCommands.o

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

PYTHON_LIBRARY = -framework Python

python:  $(PythonOBJS)
	$(LINKER) $(LINKFLAGS) $(PythonOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS)  \
	 -o pythonInterpreter

#	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace pythonExample.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o fmkSum.dylib

OBJSm = OpenSeesCommandsPython.o OpenSeesUniaxialMaterialCommands.o PythonModelBuilder.o PythonAnalysisBuilder.o

pythonmodule: $(PythonModuleOBJS)
	$(LINKER) $(LINKFLAGS) -shared $(PythonModuleOBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS)  \
	 -o opensees.so

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
