if(Python_FOUND)
  target_sources(OPS_SysOfEqnPython
    PRIVATE
      SparsePythonCompressedLinSOE.cpp
      SparsePythonCompressedLinSolver.cpp
      SparsePythonCOOLinSOE.cpp
      SparsePythonCOOLinSolver.cpp
      SparsePythonFactory.cpp
      SparsePythonCommon.cpp
    PUBLIC
      SparsePythonCompressedLinSOE.h
      SparsePythonCompressedLinSolver.h
      SparsePythonCOOLinSOE.h
      SparsePythonCOOLinSolver.h
      SparsePythonCommon.h
      SparsePythonFactory.h
  )

  target_link_libraries(OPS_SysOfEqnPython PRIVATE Python::Python)
endif()

