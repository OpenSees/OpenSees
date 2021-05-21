target_sources(actor
    PRIVATE
      MachineBroker.cpp
    PUBLIC
      MachineBroker.h
)


if(PARALLEL_INTERPRETERS)
  target_sources(actor
    PRIVATE
      MPI_Machine.cpp
    PUBLIC
      MPI_MachineBroker.h
)
elseif(PARALLEL)
   target_sources(actor
     PRIVATE
       MPI_Machine.cpp
     PUBLIC
       MPI_MachineBroker.h
)
endif()


target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
