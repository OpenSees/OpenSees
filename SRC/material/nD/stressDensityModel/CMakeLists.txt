target_sources(material
    PRIVATE
	stressDensity.cpp
	SDM-UC.f
	sdm3d.f90
	sdm2d.f90
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

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))
