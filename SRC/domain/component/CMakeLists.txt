
target_sources(domain
    PRIVATE
	DomainComponent.cpp
	Parameter.cpp
	ElementParameter.cpp
	NodeResponseParameter.cpp
	LoadFactorParameter.cpp
	MaterialStageParameter.cpp
	MatParameter.cpp
	InitialStateParameter.cpp
	ElementStateParameter.cpp
    PUBLIC
	DomainComponent.h
	Parameter.h
	ElementParameter.h
	NodeResponseParameter.h
	LoadFactorParameter.h
	MaterialStageParameter.h
	MatParameter.h
	InitialStateParameter.h
	ElementStateParameter.h
)

if(RELIABILITY)
target_sources(domain
    PRIVATE
	RVParameter.cpp
	DVParameter.cpp
    PUBLIC
	RVParameter.h
	DVParameter.h
)
endif()


target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
