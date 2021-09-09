add_library(systemOfEqn "")

target_sources(systemOfEqn
    PRIVATE
    PUBLIC
	Solver.h
	SystemOfEqn.h
)

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(linearSOE)
add_subdirectory(eigenSOE)
