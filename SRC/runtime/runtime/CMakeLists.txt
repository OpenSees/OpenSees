#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================

target_sources(OPS_Runtime
    PRIVATE
      BasicAnalysisBuilder.cpp
      BasicModelBuilder.cpp
      TclPackageClassBroker.cpp

    PUBLIC
      BasicAnalysisBuilder.h
      BasicModelBuilder.h
      TclPackageClassBroker.h
)

add_subdirectory(SectionBuilder)
add_subdirectory(Renderer)


