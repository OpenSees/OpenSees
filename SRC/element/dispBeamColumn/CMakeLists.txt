target_sources(element
    PRIVATE
	DispBeamColumn2d.cpp
	DispBeamColumn3d.cpp
	DispBeamColumn2dThermal.cpp
	DispBeamColumn3dThermal.cpp
	DispBeamColumnNL2d.cpp
	DispBeamColumn2dWithSensitivity.cpp
	DispBeamColumn3dWithSensitivity.cpp
	AxEqDispBeamColumn2d.cpp
	TimoshenkoBeamColumn2d.cpp
	DispBeamColumn3dID.cpp
    PUBLIC
	DispBeamColumn2d.h
	DispBeamColumn3d.h
	DispBeamColumn2dThermal.h
	DispBeamColumn3dThermal.h
	DispBeamColumnNL2d.h
	DispBeamColumn2dWithSensitivity.h
	DispBeamColumn3dWithSensitivity.h
	AxEqDispBeamColumn2d.h
	TimoshenkoBeamColumn2d.h
	DispBeamColumn3dID.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
