
target_sources(domain
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

if(DRM)
  target_sources(domain
    PRIVATE
	H5DRM.cpp
    PUBLIC
	H5DRM.h    
  )
endif()

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
