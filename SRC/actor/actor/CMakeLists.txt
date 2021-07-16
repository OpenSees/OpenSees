target_sources(G3_Actor
    PRIVATE
      Actor.cpp
      MovableObject.cpp
    PUBLIC
      Actor.h
      MovableObject.h
)

target_include_directories(G3_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
