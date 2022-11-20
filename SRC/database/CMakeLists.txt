#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#add_library(OPS_Database OBJECT)

target_sources(OPS_Database
    PRIVATE
        FE_Datastore.cpp
        FileDatastore.cpp
    PUBLIC
        FE_Datastore.h
        FileDatastore.h
)
target_include_directories(OPS_Database PUBLIC ${CMAKE_CURRENT_LIST_DIR})

