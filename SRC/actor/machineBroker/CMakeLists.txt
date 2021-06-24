target_sources(G3_Actor
    PRIVATE
      MachineBroker.cpp
    PUBLIC
      MachineBroker.h
)


if(PARALLEL_INTERPRETERS)
  target_sources(G3_Actor
    PRIVATE
      MPI_Machine.cpp
    PUBLIC
      MPI_MachineBroker.h
)
elseif(PARALLEL)
   target_sources(G3_Actor
     PRIVATE
       MPI_Machine.cpp
     PUBLIC
       MPI_MachineBroker.h
)
endif()


target_include_directories(G3_Actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
