target_sources(actor
    PRIVATE
      Message.cpp      
    PUBLIC
      Message.h
)

target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
