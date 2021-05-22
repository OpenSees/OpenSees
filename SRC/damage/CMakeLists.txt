add_library(damage "")

target_sources(damage
    PRIVATE
	DamageModel.cpp
	HystereticEnergy.cpp
	Kratzig.cpp
	Mehanny.cpp
	ParkAng.cpp
	NormalizedPeak.cpp
	DamageResponse.cpp
    PUBLIC
	DamageModel.h
	HystereticEnergy.h
	Kratzig.h
	Mehanny.h
	ParkAng.h
	NormalizedPeak.h
	DamageResponse.h
)

target_include_directories(damage PUBLIC $(CMAKE_CURRENT_LIST_DIR))
