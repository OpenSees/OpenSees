target_sources(G3_Actor
    PRIVATE
      ChannelAddress.cpp
    PUBLIC
      ChannelAddress.h
)

target_sources(OpenSeesMP
  PRIVATE
    MPI_ChannelAddress.cpp
  PUBLIC
    MPI_ChannelAddress.h
)

target_sources(OpenSeesSP
  PRIVATE
    MPI_ChannelAddress.cpp
  PUBLIC
    MPI_ChannelAddress.h
)

target_include_directories(G3_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})

