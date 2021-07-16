target_sources(OPS_Material
    PRIVATE
	stressDensity.cpp
	SDM-UC.f
	SDM3D.f90
	SDM2D.f90
	MODULE_ALL_INTERFACES_3D.f90
	MODULE_ALL_INTERFACES_2D.f90
	StressDensityModel.cpp
	StressDensityModel2D.cpp 
	StressDensityModel3D.cpp
    PUBLIC
	stressDensity.h
	StressDensityModel.h
	StressDensityModel2D.h
	StressDensityModel3D.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
