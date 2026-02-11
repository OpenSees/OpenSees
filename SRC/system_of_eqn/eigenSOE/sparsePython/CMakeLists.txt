target_sources(OPS_SysOfEqn
  PRIVATE
    SparsePythonCompressedEigenSOE.cpp
    SparsePythonCompressedEigenSolver.cpp
    SparsePythonCOOEigenSOE.cpp
    SparsePythonCOOEigenSolver.cpp
    SparsePythonEigenCommon.cpp
    SparsePythonEigenFactory.cpp
  PUBLIC
    SparsePythonCompressedEigenSOE.h
    SparsePythonCompressedEigenSolver.h
    SparsePythonCOOEigenSOE.h
    SparsePythonCOOEigenSolver.h
    SparsePythonEigenCommon.h
    SparsePythonEigenFactory.h
)

find_package(Python3 COMPONENTS Development REQUIRED)
target_link_libraries(OPS_SysOfEqn PRIVATE Python3::Python)

