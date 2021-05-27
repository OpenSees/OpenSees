target_sources(element
    PRIVATE
	Brick.cpp
	Twenty_Node_Brick.cpp
	shp3d.cpp
	BbarBrick.cpp
	BbarBrickWithSensitivity.cpp
    PUBLIC
	Brick.h
	Twenty_Node_Brick.h
	shp3d.h
	BbarBrick.h
	BbarBrickWithSensitivity.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
