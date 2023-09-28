#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
    SimpleContact2D.cpp
    SimpleContact3D.cpp
    BeamContact2D.cpp
    BeamContact2Dp.cpp
    BeamContact3D.cpp
    BeamContact3Dp.cpp
    BeamEndContact3D.cpp
    BeamEndContact3Dp.cpp
    SSPquad.cpp
    SSPquadUP.cpp
    SSPbrick.cpp
    SSPbrickUP.cpp
    Quad4FiberOverlay.cpp
    Brick8FiberOverlay.cpp
    QuadBeamEmbedContact.cpp
    PileToe3D.cpp
    EmbeddedBeamInterfaceL.cpp
    EmbeddedBeamInterfaceP.cpp
    EmbeddedEPBeamInterface.cpp
    Tcl_generateInterfacePoints.cpp
    PUBLIC
    SimpleContact2D.h
    SimpleContact3D.h
    BeamContact2D.h
    BeamContact2Dp.h
    BeamContact3D.h
    BeamContact3Dp.h
    BeamEndContact3D.h
    BeamEndContact3Dp.h
    SSPquad.h
    SSPquadUP.h
    SSPbrick.h
    SSPbrickUP.h
    Quad4FiberOverlay.h
    Brick8FiberOverlay.h
    QuadBeamEmbedContact.h
    PileToe3D.h
    EmbeddedBeamInterfaceL.h
    EmbeddedBeamInterfaceP.h
    EmbeddedEPBeamInterface.h
    Tcl_generateInterfacePoints.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
