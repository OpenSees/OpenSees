#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_ModelBuilder
    PRIVATE
        ModelBuilder.cpp
    #PartitionedModelBuilder.cpp
    PUBLIC
        ModelBuilder.h
    #PartitionedModelBuilder.h
)


add_library(OPS_ModelBuilder_PlaneFrame EXCLUDE_FROM_ALL)
target_sources(OPS_ModelBuilder_PlaneFrame PRIVATE PlaneFrame.cpp)
target_link_libraries(OPS_ModelBuilder_PlaneFrame PRIVATE OPS_ModelBuilder)


add_subdirectory(tcl)

