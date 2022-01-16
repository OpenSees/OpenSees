#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_mvlem OBJECT)

target_sources(OPS_Element
    PRIVATE
        MVLEM_3D.cpp
        MVLEM.cpp
        SFI_MVLEM_3D.cpp
        SFI_MVLEM.cpp
    PUBLIC
        MVLEM_3D.h
        MVLEM.h
        SFI_MVLEM_3D.h
        SFI_MVLEM.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

