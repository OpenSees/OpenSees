#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    MP_Constraint.cpp 
    RigidBeam.cpp 
    RigidDiaphragm.cpp
    RigidRod.cpp 
    SP_Constraint.cpp 
    Pressure_Constraint.cpp
  PUBLIC
    MP_Constraint.h 
    RigidBeam.h 
    RigidDiaphragm.h
    RigidRod.h 
    SP_Constraint.h 
    Pressure_Constraint.h
)
target_sources(OPS_Domain
  PRIVATE
    ImposedMotionSP.cpp 
    ImposedMotionSP1.cpp
  PUBLIC
    ImposedMotionSP.h 
    ImposedMotionSP1.h
)
target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
