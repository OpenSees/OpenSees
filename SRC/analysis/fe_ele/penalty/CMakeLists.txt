target_sources(G3_Analysis
    PRIVATE
      PenaltySP_FE.cpp
      PenaltyMP_FE.cpp
    PUBLIC
      PenaltySP_FE.h
      PenaltyMP_FE.h
)

#target_include_directories(G3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
