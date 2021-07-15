add_library(tagged "")

target_sources(tagged
    PRIVATE
      TaggedObject.cpp
    PUBLIC
      TaggedObject.h
)

target_include_directories(tagged PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(storage)
