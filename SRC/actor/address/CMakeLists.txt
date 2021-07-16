target_sources(G3_Actor
    PRIVATE
      ChannelAddress.cpp
    PUBLIC
      ChannelAddress.h
)


if(PARALLEL_INTERPRETERS)
target_sources(G3_Actor
    PRIVATE
      MPI_Address.cpp
    PUBLIC
      MPI_Address.h
)
elseif(PARALLEL)
target_sources(G3_Actor
    PRIVATE
      MPI_Address.cpp
    PUBLIC
      MPI_Address.h
)
endif()


target_include_directories(G3_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})

