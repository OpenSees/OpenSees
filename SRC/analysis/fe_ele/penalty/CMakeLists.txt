target_sources(analysis
    PRIVATE
      PenaltySP_FE.cpp
      PenaltyMP_FE.cpp
    PUBLIC
      PenaltySP_FE.h
      PenaltyMP_FE.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
