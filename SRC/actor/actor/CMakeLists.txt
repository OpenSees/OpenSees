target_sources(actor
    PRIVATE
      Actor.cpp
      MovableObject.cpp
    PUBLIC
      Actor.h
      MovableObject.h
)

target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
