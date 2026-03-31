#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
      ElasticBeam2d.cpp
      ElasticBeam3d.cpp
      ElasticBeamWarping3d.cpp      
      ElasticTimoshenkoBeam2d.cpp
      ElasticTimoshenkoBeam3d.cpp
      ModElasticBeam2d.cpp
      ModElasticBeam3d.cpp
      WheelRail.cpp
    PUBLIC
      ElasticBeam2d.h
      ElasticBeam3d.h
      ElasticBeamWarping3d.h      
      ElasticTimoshenkoBeam2d.h
      ElasticTimoshenkoBeam3d.h
      ModElasticBeam2d.h
      ModElasticBeam3d.h
      WheelRail.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
