#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Material OBJECT)

target_sources(OPS_Material
    PRIVATE
      Material.cpp
    PUBLIC
      Material.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(uniaxial)
add_subdirectory(nD)
add_subdirectory(section)
add_subdirectory(yieldSurface)

