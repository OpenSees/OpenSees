
target_sources(domain
    PRIVATE
	Domain.cpp
    PUBLIC
	Domain.h
	ElementIter.h
	LoadCaseIter.h
	MP_ConstraintIter.h
	NodeIter.h
	ParameterIter.h
	Pressure_ConstraintIter.h
	SP_ConstraintIter.h
	SubdomainIter.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(single)
add_subdirectory(partitioned)
