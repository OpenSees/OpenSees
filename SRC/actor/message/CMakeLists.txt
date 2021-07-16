target_sources(G3_Actor
    PRIVATE
      Message.cpp      
    PUBLIC
      Message.h
)

target_include_directories(G3_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
