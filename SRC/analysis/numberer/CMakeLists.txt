
target_sources(analysis
    PRIVATE
      DOF_Numberer.cpp
      PlainNumberer.cpp
      ParallelNumberer.cpp
    PUBLIC
      DOF_Numberer.h
      PlainNumberer.h
      ParallelNumberer.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
