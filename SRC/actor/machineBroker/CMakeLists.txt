#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Actor
    PRIVATE
      MachineBroker.cpp
    PUBLIC
      MachineBroker.h
)



target_sources(OpenSeesMP
    PRIVATE
      MPI_MachineBroker.cpp
    PUBLIC
      MPI_MachineBroker.h
)

target_sources(OpenSeesSP
     PRIVATE
       MPI_MachineBroker.cpp
     PUBLIC
       MPI_MachineBroker.h
)


target_include_directories(OPS_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})

