
target_sources(material
    PRIVATE
	UnloadingRule.cpp
	ConstantUnloadingRule.cpp 
	EnergyUnloadingRule.cpp
	KarsanUnloadingRule.cpp 
	TakedaUnloadingRule.cpp
    PUBLIC
	UnloadingRule.h
	ConstantUnloadingRule.h 
	EnergyUnloadingRule.h
	KarsanUnloadingRule.h 
	TakedaUnloadingRule.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))


