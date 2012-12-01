include ../../Makefile.def

OBJS       = DL_Interpreter.o \
	TclInterpreter.o

# Compilation control

all:         $(OBJS)

tcl: tclMain.o $(OBJS)
	$(LINKER) $(LINKFLAGS) tclMain.o DL_Interpreter.o TclInterpreter.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
	 -o tclInterpreter

PYTHON_LIBRARY = -framework python

python: pythonMain.o PythonInterpreter.o $(OBJS)
	$(LINKER) $(LINKFLAGS) pythonMain.o DL_Interpreter.o PythonInterpreter.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PYTHON_LIBRARY) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
	 -o pythonInterpreter

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
