
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
    OpenSeesBeamIntegrationCommands.cpp
    OpenSeesNDMaterialCommands.cpp
    OpenSeesMiscCommands.cpp
    OpenSeesParameterCommands.cpp
    OpenSeesFrictionModelCommands.cpp
    OpenSeesReliabilityCommands.cpp
)


target_sources(OpenSeesPy
   PRIVATE
     PythonModule.cpp
     PythonWrapper.cpp
)
