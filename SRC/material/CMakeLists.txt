#add_library(G3_Material OBJECT)

target_sources(G3_Material
    PRIVATE
      Material.cpp
    PUBLIC
      Material.h
)

target_include_directories(OPS_Material PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(uniaxial)
add_subdirectory(nD)
add_subdirectory(section)
add_subdirectory(yieldSurface)

