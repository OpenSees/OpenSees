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
    CorotCrdTransfWarping3d.cpp    
    PUBLIC
    CrdTransf.h
    LinearCrdTransf2d.h
    LinearCrdTransf3d.h
    PDeltaCrdTransf2d.h
    PDeltaCrdTransf3d.h
    CorotCrdTransf2d.h
    CorotCrdTransf3d.h
    CorotCrdTransfWarping2d.h
    CorotCrdTransfWarping3d.h        
)

target_include_directories(coordTransformation PUBLIC ${CMAKE_CURRENT_LIST_DIR})
