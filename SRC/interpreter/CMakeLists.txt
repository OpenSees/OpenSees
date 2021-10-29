#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_library(OPS_Commands)
target_sources(OPS_Commands
  PRIVATE
    DL_Interpreter.cpp
    OpenSeesCommands.cpp
    OpenSeesUniaxialMaterialCommands.cpp
    OpenSeesElementCommands.cpp
    OpenSeesTimeSeriesCommands.cpp
    OpenSeesPatternCommands.cpp
    OpenSeesSectionCommands.cpp
    OpenSeesOutputCommands.cpp
    OpenSeesCrdTransfCommands.cpp
    OpenSeesBeamIntegrationCommands.cpp
    OpenSeesNDMaterialCommands.cpp
    OpenSeesMiscCommands.cpp
    OpenSeesParameterCommands.cpp
    OpenSeesFrictionModelCommands.cpp
    OpenSeesReliabilityCommands.cpp
)

# Tcl interpreter
#-------------------------------------------------------------------
add_library(OPS_Interpreter_Tcl)
target_sources(OPS_Interpreter_Tcl
  PRIVATE
    tclMain.cpp
    TclInterpreter.cpp
    TclWrapper.cpp
)
target_link_libraries(OPS_Interpreter_Tcl PRIVATE OPS_Commands)
#tcl: $(TclOBJS)
#	$(LINKER) $(LINKFLAGS) $(TclOBJS)
#	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY)
#	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(PARALLEL_LIB)
#	$(INTERPRETER_LIBS_MPI_PARAMETERS)
#	 -o openseestcl

# Python interpreter
#-------------------------------------------------------------------
add_library(OPS_Interpreter_Python)
target_sources(OPS_Interpreter_Python
  PRIVATE
    pythonMain.cpp
    PythonModule.cpp
    PythonWrapper.cpp
)
target_link_libraries(OPS_Interpreter_Python PRIVATE OPS_Commands)
target_include_directories(OPS_Interpreter_Python PRIVATE ${Python_INCLUDE_DIRS})

# Python module
#-------------------------------------------------------------------
add_library(OPS_InterpreterModule_Python OBJECT)
target_sources(OPS_InterpreterModule_Python
  PRIVATE
    PythonModule.cpp
    PythonWrapper.cpp
)
target_link_libraries(OPS_InterpreterModule_Python PRIVATE OPS_Commands)
target_include_directories(OPS_InterpreterModule_Python PRIVATE ${Python_INCLUDE_DIRS})



#OBJSm = OpenSeesCommandsPython.o OpenSeesUniaxialMaterialCommands.o PythonModelBuilder.o PythonAnalysisBuilder.o


#INTERPRETER_LIBS_MPI_PARAMETERS =
#
#ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)
#INTERPRETER_LIBS_MPI_PARAMETERS
#	$(FE)/actor/channel/MPI_Channel.o
#	$(FE)/actor/machineBroker/MPI_MachineBroker.o
#	$(FE)/actor/machineBroker/MachineBroker.o
#	$(FE)/actor/actor/Actor.o
#	$(FE)/actor/shadow/Shadow.o
#	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.o
#	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.o
#	$(FE)/domain/subdomain/ActorSubdomain.o
#	$(FE)/domain/subdomain/ShadowSubdomain.o
#endif

