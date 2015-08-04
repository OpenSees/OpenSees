
-include ../../Makefile.def

OBJS       = DL_Interpreter.o \
	TclInterpreter.o \
	PythonInterpreter.o \
	OpenSeesCommandsTcl.o \
	OpenSeesCommandsPython.o \
	OpenSeesUniaxialMaterialCommands.o

# Compilation control

all:         python

tcl: tclMain.o $(OBJS)
	$(LINKER) $(LINKFLAGS) tclMain.o DL_Interpreter.o TclInterpreter.o OpenSeesCommandsTcl.o OpenSeesUniaxialMaterialCommands.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(FE)/api/elementAPI_Tcl.o \
	 -o tclInterpreter

PYTHON_LIBRARY = -framework python

python: pythonMain.o PythonInterpreter.o OpenSeesCommandsPython.o $(OBJS)
	$(LINKER) $(LINKFLAGS) pythonMain.o DL_Interpreter.o PythonInterpreter.o OpenSeesCommandsPython.o OpenSeesUniaxialMaterialCommands.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS)  \
	 -o pythonInterpreter

#	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace pythonExample.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o fmkSum.dylib

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
