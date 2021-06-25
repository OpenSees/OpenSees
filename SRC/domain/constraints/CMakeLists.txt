
target_sources(G3_Domain
    PRIVATE
	MP_Constraint.cpp 
	RigidBeam.cpp 
	RigidDiaphragm.cpp
	RigidRod.cpp 
	SP_Constraint.cpp 
    PUBLIC
	MP_Constraint.h 
	RigidBeam.h 
	RigidDiaphragm.h
	RigidRod.h 
	SP_Constraint.h 
)
target_sources(OPS_Domain
    PRIVATE
	ImposedMotionSP.cpp 
	ImposedMotionSP1.cpp
	Pressure_Constraint.cpp
    PUBLIC
	ImposedMotionSP.h 
	ImposedMotionSP1.h
	Pressure_Constraint.h
)
target_include_directories(OPS_Domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
