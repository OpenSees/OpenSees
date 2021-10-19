#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(coordTransformation "")

target_sources(coordTransformation
    PRIVATE
    CrdTransf.cpp
    LinearCrdTransf2d.cpp
    LinearCrdTransf3d.cpp
    PDeltaCrdTransf2d.cpp
    PDeltaCrdTransf3d.cpp
    CorotCrdTransf2d.cpp
    CorotCrdTransf3d.cpp
    CorotCrdTransfWarping2d.cpp
    PUBLIC
    CrdTransf.cpp
    LinearCrdTransf2d.cpp
    LinearCrdTransf3d.cpp
    PDeltaCrdTransf2d.cpp
    PDeltaCrdTransf3d.cpp
    CorotCrdTransf2d.cpp
    CorotCrdTransf3d.cpp
    CorotCrdTransfWarping2d.cpp
)

target_include_directories(coordTransformation PUBLIC ${CMAKE_CURRENT_LIST_DIR})
