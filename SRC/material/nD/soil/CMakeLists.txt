#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    FluidSolidPorousMaterial.cpp
    MultiYieldSurface.cpp
    MultiYieldSurfaceClay.cpp
    PressureDependMultiYield.cpp
    PressureIndependMultiYield.cpp
    T2Vector.cpp
    PressureDependMultiYield02.cpp
    PressureDependMultiYield03.cpp
    PUBLIC
    FluidSolidPorousMaterial.h
    MultiYieldSurface.h
    MultiYieldSurfaceClay.h
    PressureDependMultiYield.h
    PressureIndependMultiYield.h
    T2Vector.h
    PressureDependMultiYield02.h
    PressureDependMultiYield03.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
