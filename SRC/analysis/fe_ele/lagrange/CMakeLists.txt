target_sources(analysis
    PRIVATE
      LagrangeSP_FE.cpp
      LagrangeMP_FE.cpp
    PUBLIC
      LagrangeSP_FE.h
      LagrangeMP_FE.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
