target_sources(G3_Actor
    PRIVATE
      Channel.cpp
      HTTP.cpp
      Socket.cpp
      TCP_Socket.cpp
      UDP_Socket.cpp      
    PUBLIC
      Channel.h
      Socket.h
      TCP_Socket.h
      UDP_Socket.h      
)


if(PARALLEL_INTERPRETERS)
target_sources(G3_Actor
    PRIVATE
      MPI_Channel.cpp
    PUBLIC
      MPI_Channel.h
)
elseif(PARALLEL)
target_sources(G3_Actor
    PRIVATE
      MPI_Channel.cpp
    PUBLIC
      MPI_Channel.h
)
endif()


target_include_directories(G3_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
