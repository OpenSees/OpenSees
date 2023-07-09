#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#------------------------------------------------------------------------------
# Claudio Perez

add_library(ITPACK STATIC
    dsrc2c.f
    jcg.f
    jsi.f
    rscg.f
    rssi.f
    sor.f
    ssorcg.f
    ssorsi.f
)

