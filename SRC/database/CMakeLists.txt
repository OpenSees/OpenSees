#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#add_library(OPS_Database OBJECT)

target_sources(OPS_Database
    PRIVATE
        #BerkeleyDbDatastore.cpp
        FE_Datastore.cpp
        FileDatastore.cpp
        MySqlDatastore.cpp
        OracleDatastore.cpp
    PUBLIC
        #BerkeleyDbDatastore.h
        FE_Datastore.h
        FileDatastore.h
        MySqlDatastore.h
        OracleDatastore.h
)
target_include_directories(OPS_Database PUBLIC ${CMAKE_CURRENT_LIST_DIR})

