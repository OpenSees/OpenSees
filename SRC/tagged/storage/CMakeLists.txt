target_sources(tagged
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

target_include_directories(tagged PUBLIC $(CMAKE_CURRENT_LIST_DIR))

