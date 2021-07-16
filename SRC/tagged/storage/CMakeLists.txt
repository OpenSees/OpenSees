target_sources(G3_Tagged
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

target_include_directories(G3_Tagged PUBLIC ${CMAKE_CURRENT_LIST_DIR})

