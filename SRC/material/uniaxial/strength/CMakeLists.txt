
target_sources(material
    PRIVATE
	StrengthDegradation.cpp
	ACIStrengthDegradation.cpp 
	ConstantStrengthDegradation.cpp
	DuctilityStrengthDegradation.cpp 
	EnergyStrengthDegradation.cpp
	PetrangeliStrengthDegradation.cpp 
	SectionStrengthDegradation.cpp
    PUBLIC
	StrengthDegradation.h
	ACIStrengthDegradation.h 
	ConstantStrengthDegradation.h
	DuctilityStrengthDegradation.h 
	EnergyStrengthDegradation.h
	PetrangeliStrengthDegradation.h 
	SectionStrengthDegradation.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))


