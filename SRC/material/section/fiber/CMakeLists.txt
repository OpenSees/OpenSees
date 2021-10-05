#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    Fiber.cpp
    UniaxialFiber2d.cpp
    UniaxialFiber3d.cpp
    NDFiber2d.cpp
    NDFiber3d.cpp
    PUBLIC
    Fiber.h
    UniaxialFiber2d.h
    UniaxialFiber3d.h
    NDFiber2d.h
    NDFiber3d.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

