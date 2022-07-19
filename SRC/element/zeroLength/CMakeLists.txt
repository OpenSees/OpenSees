#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
    ZeroLength.cpp
    CoupledZeroLength.cpp
    ZeroLengthContact2D.cpp
    ZeroLengthContact3D.cpp
    ZeroLengthContactASDimplex.cpp
    ZeroLengthContactNTS2D.cpp
    ZeroLengthImpact3D.cpp
    ZeroLengthInterface2D.cpp
    ZeroLengthND.cpp
    ZeroLengthRocking.cpp
    ZeroLengthSection.cpp
    ZeroLengthVG_HG.cpp
    PUBLIC
    ZeroLength.h
    CoupledZeroLength.h
    ZeroLengthContact2D.h
    ZeroLengthContact3D.h
    ZeroLengthContactASDimplex.h
    ZeroLengthContactNTS2D.h
    ZeroLengthImpact3D.h
    ZeroLengthInterface2D.h
    ZeroLengthND.h
    ZeroLengthRocking.h
    ZeroLengthSection.h
    ZeroLengthVG_HG.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
