#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        BeamColumnJoint2d.cpp
        BeamColumnJoint3d.cpp
        ElasticTubularJoint.cpp
        Joint2D.cpp
        Joint3D.cpp
		Inno3DPnPJoint.cpp
        LehighJoint2d.cpp
        MP_Joint2D.cpp
        MP_Joint3D.cpp
        #TclBeamColumnJointCommand.cpp
        #TclJoint2dCommand.cpp
        #TclJoint3dCommand.cpp
    PUBLIC
        BeamColumnJoint2d.h
        BeamColumnJoint3d.h
        ElasticTubularJoint.h
        Joint2D.h
        Joint3D.h
		Inno3DPnPJoint.h
        LehighJoint2d.h
        MP_Joint2D.h
        MP_Joint3D.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

