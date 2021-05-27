target_sources(actor
    PRIVATE
      ChannelAddress.cpp
    PUBLIC
      ChannelAddress.h
)


if(PARALLEL_INTERPRETERS)
target_sources(actor
    PRIVATE
      MPI_Address.cpp
    PUBLIC
      MPI_Address.h
)
elseif(PARALLEL)
target_sources(actor
    PRIVATE
      MPI_Address.cpp
    PUBLIC
      MPI_Address.h
)
endif()


target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))

