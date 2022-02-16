#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_updatedLagrangianBeamColumn OBJECT)

target_sources(OPS_Element
    PRIVATE
        BilinearCyclic.cpp
        CyclicModel.cpp
        Elastic2DGNL.cpp
        Inelastic2DYS01.cpp
        Inelastic2DYS02.cpp
        Inelastic2DYS03.cpp
        InelasticYS2DGNL.cpp
        LinearCyclic.cpp
        QuadraticCyclic.cpp
        #TclCyclicModelCommands.cpp
        #TclElement2dGNL.cpp
        #TclElement2dYS.cpp
        UpdatedLagrangianBeam2D.cpp
    PUBLIC
        BilinearCyclic.h
        CyclicModel.h
        Elastic2DGNL.h
        Inelastic2DYS01.h
        Inelastic2DYS02.h
        Inelastic2DYS03.h
        InelasticYS2DGNL.h
        LinearCyclic.h
        QuadraticCyclic.h
        UpdatedLagrangianBeam2D.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

