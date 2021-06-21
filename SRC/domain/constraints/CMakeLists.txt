
target_sources(libg3
    PRIVATE
	ImposedMotionSP.cpp 
	ImposedMotionSP1.cpp
	SP_Constraint.cpp 
	MP_Constraint.cpp 
	RigidDiaphragm.cpp
	RigidRod.cpp 
	RigidBeam.cpp 
    PUBLIC
	ImposedMotionSP.h 
	ImposedMotionSP1.h
	SP_Constraint.h 
	MP_Constraint.h 
	RigidDiaphragm.h
	RigidRod.h 
	RigidBeam.h 
)
target_sources(domain
    PRIVATE
	Pressure_Constraint.cpp
    PUBLIC
	Pressure_Constraint.h
)
target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
