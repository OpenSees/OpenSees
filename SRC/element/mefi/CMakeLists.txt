#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_mefi OBJECT)

target_sources(OPS_Element
    PRIVATE
        MEFI.cpp
		MEFI_3D.cpp
    PUBLIC
        MEFI.h
		MEFI_3D.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

