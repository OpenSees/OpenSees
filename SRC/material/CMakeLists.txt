add_library(material OBJECT)

target_sources(G3_SysOfEqn
    PRIVATE
      Material.cpp
    PUBLIC
      Material.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(uniaxial)
add_subdirectory(nD)
add_subdirectory(section)
