#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Domain
    PRIVATE
    Load.cpp
    #ElementalLoad.cpp
    ElementalLoadIter.cpp
    NodalLoadIter.cpp
    Beam2dUniformLoad.cpp
    Beam2dPointLoad.cpp
    Beam2dPartialUniformLoad.cpp
    Beam3dUniformLoad.cpp
    Beam3dPointLoad.cpp
    Beam3dPartialUniformLoad.cpp
    BrickSelfWeight.cpp
    SurfaceLoader.cpp
    SelfWeight.cpp
    LysmerVelocityLoader.cpp
    PUBLIC
    Load.h
    #ElementalLoad.h
    ElementalLoadIter.h
    NodalLoadIter.h
    Beam2dUniformLoad.h
    Beam2dPointLoad.h
    Beam2dPartialUniformLoad.h
    Beam3dUniformLoad.h
    Beam3dPointLoad.h
    Beam3dPartialUniformLoad.h
    BrickSelfWeight.h
    SurfaceLoader.h
    SelfWeight.h
    LysmerVelocityLoader.h
)

target_sources(OPS_Thermal
    PRIVATE
        Beam2dTempLoad.cpp
        Beam2dThermalAction.cpp
        NodalThermalAction.cpp
        ThermalActionWrapper.cpp
        Beam3dThermalAction.cpp
        ShellThermalAction.cpp
    PUBLIC
        Beam2dTempLoad.h
        Beam2dThermalAction.h
        NodalThermalAction.h
        ThermalActionWrapper.h
        Beam3dThermalAction.h
        ShellThermalAction.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
