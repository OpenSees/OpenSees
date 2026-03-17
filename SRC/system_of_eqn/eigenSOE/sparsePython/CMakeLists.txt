if(Python_FOUND)
  target_sources(OPS_SysOfEqnPython
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

  target_link_libraries(OPS_SysOfEqnPython PRIVATE Python::Python)
endif()

