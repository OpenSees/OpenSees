target_sources(G3_SysOfEqn
    PRIVATE
    PUBLIC
	Solver.h
	SystemOfEqn.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(linearSOE)
add_subdirectory(eigenSOE)
