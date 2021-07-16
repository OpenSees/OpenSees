
target_sources(OPS_DRM
    PRIVATE
    Mesh3DSubdomain.cpp
    GeometricBrickDecorator.cpp
    DRMBoundaryLayerDecorator.cpp
    DRMLoadPattern.cpp
    DRMLoadPatternWrapper.cpp
    DRMInputHandler.cpp
    PlaneDRMInputHandler.cpp
    PUBLIC
    Mesh3DSubdomain.h
    GeometricBrickDecorator.h
    DRMBoundaryLayerDecorator.h
    DRMLoadPattern.h
    DRMLoadPatternWrapper.h
    DRMInputHandler.h
    PlaneDRMInputHandler.h
)
if (OPS_Use_HDF5)
target_sources(OPS_DRM
    PRIVATE
        H5DRM.cpp
    PUBLIC
        H5DRM.h    
)
endif()

target_include_directories(OPS_DRM PUBLIC $(CMAKE_CURRENT_LIST_DIR))

