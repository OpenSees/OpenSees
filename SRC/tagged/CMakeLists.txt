target_sources(G3_Tagged
    PRIVATE
      TaggedObject.cpp
    PUBLIC
      TaggedObject.h
)

target_include_directories(G3_Tagged PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(storage)
