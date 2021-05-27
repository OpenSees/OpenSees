add_library(material "")

target_sources(material
    PRIVATE
      Material.cpp
    PUBLIC
      Material.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(uniaxial)
add_subdirectory(nD)
add_subdirectory(section)
