
target_sources(OPS_INTERPRETER
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
    OpenSeesDampingCommands.cpp
    OpenSeesBeamIntegrationCommands.cpp
    OpenSeesNDMaterialCommands.cpp
    OpenSeesMiscCommands.cpp
    OpenSeesParameterCommands.cpp
    OpenSeesFrictionModelCommands.cpp
    OpenSeesReliabilityCommands.cpp
    OpenSeesNDTestCommands.cpp
    OpenSeesIGACommands.cpp
)


target_sources(OpenSeesPy
   PRIVATE
     PythonModule.cpp
     PythonWrapper.cpp
     PythonStream.cpp
)
