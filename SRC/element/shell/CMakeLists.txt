#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        ASDShellQ4.cpp
        R3vectors.cpp
        ShellANDeS.cpp
        ShellDKGQ.cpp
        ShellDKGT.cpp
        ShellMITC4.cpp
        ShellMITC4Thermal.cpp
        ShellMITC9.cpp
        ShellNLDKGQ.cpp
        ShellNLDKGQThermal.cpp
        ShellNLDKGT.cpp
    PUBLIC
        ASDEICR.h
        ASDMath.h
        ASDShellQ4CorotationalTransformation.h
        ASDShellQ4.h
        ASDShellQ4LocalCoordinateSystem.h
        ASDShellQ4Transformation.h
        R3vectors.h
        ShellANDeS.h
        ShellDKGQ.h
        ShellDKGT.h
        ShellMITC4.h
        ShellMITC4Thermal.h
        ShellMITC9.h
        ShellNLDKGQ.h
        ShellNLDKGQThermal.h
        ShellNLDKGT.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

