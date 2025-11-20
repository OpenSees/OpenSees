#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    SectionIntegration.cpp 
    WideFlangeSectionIntegration.cpp 
    RCSectionIntegration.cpp
    RCTBeamSectionIntegration.cpp 
    RCCircularSectionIntegration.cpp
    TubeSectionIntegration.cpp 
    RCTunnelSectionIntegration.cpp
    MVLEMSectionIntegration.cpp
    SFIMVLEMSectionIntegration.cpp
    HSSSectionIntegration.cpp    
    PUBLIC
    SectionIntegration.h 
    WideFlangeSectionIntegration.h 
    RCSectionIntegration.h
    RCTBeamSectionIntegration.h 
    RCCircularSectionIntegration.h
    TubeSectionIntegration.h 
    RCTunnelSectionIntegration.h
    MVLEMSectionIntegration.h
    SFIMVLEMSectionIntegration.h
    HSSSectionIntegration.h    
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
