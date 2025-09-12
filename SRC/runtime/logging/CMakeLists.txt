#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================

target_include_directories(OPS_Logging PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Logging 
  PRIVATE
    "logging.cpp"
  PUBLIC
    "Logging.h"
    "AnsiColors.h"
)

