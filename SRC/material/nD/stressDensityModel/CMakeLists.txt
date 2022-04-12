#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
  PRIVATE
    stressDensity.cpp
    #StressDensityModel.cpp
    StressDensityModel2D.cpp 
    StressDensityModel3D.cpp
  PUBLIC
    stressDensity.h
    #StressDensityModel.h
    StressDensityModel2D.h
    StressDensityModel3D.h
)

#add_library(OPS_MaterialFortran)
target_sources(OPS_MaterialFortran
  PRIVATE
    SDM-UC.f
    SDM3D.f90
    SDM2D.f90
    MODULE_ALL_INTERFACES_3D.f90
    MODULE_ALL_INTERFACES_2D.f90
)

#target_link_libraries(OPS_Material_StressDensity PRIVATE OPS_Material_nD_StressDensity_f)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

