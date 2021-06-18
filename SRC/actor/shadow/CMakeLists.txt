target_sources(actor
    PRIVATE
      Shadow.cpp      
    PUBLIC
      Shadow.h
)

target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
