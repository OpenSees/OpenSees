#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#------------------------------------------------------------------------------
# Claudio Perez
project(ITPACK LANGUAGES Fortran)
add_library(ITPACK)
target_sources(ITPACK
  PRIVATE
    dsrc2c.f
    jcg.f
    jsi.f
    rscg.f
    rssi.f
    sor.f
    ssorcg.f
    ssorsi.f
)

