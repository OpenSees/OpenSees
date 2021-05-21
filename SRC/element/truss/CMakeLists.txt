target_sources(element
    PRIVATE
	Truss.cpp
	TrussSection.cpp
	CorotTruss.cpp
	CorotTrussSection.cpp
	Truss2.cpp
	N4BiaxialTruss.cpp
	CorotTruss2.cpp
    PUBLIC
	Truss.h
	TrussSection.h
	CorotTruss.h
	CorotTrussSection.h
	Truss2.h
	N4BiaxialTruss.h
	CorotTruss2.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
