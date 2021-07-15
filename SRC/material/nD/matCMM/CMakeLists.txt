target_sources(material
    PRIVATE
	MaterialCMM.cpp
    PUBLIC
	MaterialCMM.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))
