target_sources(G3_Analysis
    PRIVATE
      LagrangeSP_FE.cpp
      LagrangeMP_FE.cpp
    PUBLIC
      LagrangeSP_FE.h
      LagrangeMP_FE.h
)

target_include_directories(G3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
