
target_sources(G3_Analysis
    PRIVATE
      DOF_Numberer.cpp
      PlainNumberer.cpp
      ParallelNumberer.cpp
    PUBLIC
      DOF_Numberer.h
      PlainNumberer.h
      ParallelNumberer.h
)

#target_include_directories(G3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
