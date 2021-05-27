target_sources(element
    PRIVATE
      ElasticBeam2d.cpp
      ElasticBeam3d.cpp
      ElasticTimoshenkoBeam2d.cpp
      ElasticTimoshenkoBeam3d.cpp
      ModElasticBeam2d.cpp
      WheelRail.cpp
    PUBLIC
      ElasticBeam2d.h
      ElasticBeam3d.h
      ElasticTimoshenkoBeam2d.h
      ElasticTimoshenkoBeam3d.h
      ModElasticBeam2d.h
      WheelRail.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
