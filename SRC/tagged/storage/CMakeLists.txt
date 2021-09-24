#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Tagged
    PRIVATE
      ArrayOfTaggedObjects.cpp 
      ArrayOfTaggedObjectsIter.cpp
      MapOfTaggedObjectsIter.cpp 
      MapOfTaggedObjects.cpp
    PUBLIC
      ArrayOfTaggedObjects.h 
      ArrayOfTaggedObjectsIter.h
      MapOfTaggedObjectsIter.h 
      MapOfTaggedObjects.h
)

target_include_directories(OPS_Tagged PUBLIC ${CMAKE_CURRENT_LIST_DIR})

