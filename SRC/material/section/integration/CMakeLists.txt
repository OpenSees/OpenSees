target_sources(material
    PRIVATE
	SectionIntegration.cpp 
	WideFlangeSectionIntegration.cpp 
	RCSectionIntegration.cpp
	RCTBeamSectionIntegration.cpp 
	RCCircularSectionIntegration.cpp
	TubeSectionIntegration.cpp 
	RCTunnelSectionIntegration.cpp
    PUBLIC
	SectionIntegration.h 
	WideFlangeSectionIntegration.h 
	RCSectionIntegration.h
	RCTBeamSectionIntegration.h 
	RCCircularSectionIntegration.h
	TubeSectionIntegration.h 
	RCTunnelSectionIntegration.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))
