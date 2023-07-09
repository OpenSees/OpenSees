#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(damping "")

target_sources(damping
    PRIVATE
    Damping.cpp
    SecStifDamping.cpp
    UniformDamping.cpp
    URDDamping.cpp
    URDDampingbeta.cpp
    PUBLIC
    Damping.h
    SecStifDamping.h
    UniformDamping.h
    URDDamping.h
    URDDampingbeta.h
)

target_include_directories(damping PUBLIC ${CMAKE_CURRENT_LIST_DIR})
