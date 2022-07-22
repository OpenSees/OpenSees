#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# TODO: Create new dispBeamColumn target
target_sources(OPS_Element
  PRIVATE
    DispBeamColumn2d.cpp
    DispBeamColumn3d.cpp
    DispBeamColumnAsym3d.cpp
    DispBeamColumn2dThermal.cpp
    DispBeamColumn3dThermal.cpp
    DispBeamColumnNL2d.cpp
    DispBeamColumnNL3d.cpp
    DispBeamColumnWarping3d.cpp    
    DispBeamColumn2dWithSensitivity.cpp
    DispBeamColumn3dWithSensitivity.cpp
    AxEqDispBeamColumn2d.cpp
    TimoshenkoBeamColumn2d.cpp
    DispBeamColumn3dID.cpp
  PUBLIC
    DispBeamColumn2d.h
    DispBeamColumn3d.h
    DispBeamColumnAsym3d.h
    DispBeamColumn2dThermal.h
    DispBeamColumn3dThermal.h
    DispBeamColumnNL2d.h
    DispBeamColumnNL3d.h
    DispBeamColumnWarping3d.h    
    DispBeamColumn2dWithSensitivity.h
    DispBeamColumn3dWithSensitivity.h
    AxEqDispBeamColumn2d.h
    TimoshenkoBeamColumn2d.h
    DispBeamColumn3dID.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
