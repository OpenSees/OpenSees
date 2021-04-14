
target_sources(domain
    PRIVATE
	SP_Constraint.cpp 
	MP_Constraint.cpp 
	RigidDiaphragm.cpp
	RigidRod.cpp 
	RigidBeam.cpp 
	ImposedMotionSP.cpp 
	ImposedMotionSP1.cpp
	Pressure_Constraint.cpp
    PUBLIC
	SP_Constraint.h 
	MP_Constraint.h 
	RigidDiaphragm.h
	RigidRod.h 
	RigidBeam.h 
	ImposedMotionSP.h 
	ImposedMotionSP1.h
	Pressure_Constraint.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
