target_sources(OPS_Tagged
    PRIVATE
      TaggedObject.cpp
    PUBLIC
      TaggedObject.h
)

target_include_directories(OPS_Tagged PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(storage)
