target_sources(OPS_ModelBuilder
    PRIVATE
        ModelBuilder.cpp
    #PartitionedModelBuilder.cpp
    PUBLIC
        ModelBuilder.h
    #PartitionedModelBuilder.h
)
add_subdirectory(tcl)

