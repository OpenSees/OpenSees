## Regression Tests

This directory contains deterministic verification tests using OpenSeesPy. The
same suite is registered with CTest when `OPS_ENABLE_TESTS=ON` so it can run
against a freshly built module in local builds and continuous integration.

### How to Run

Run the tests directly with `pytest` when an OpenSeesPy installation is
available:

```console
pytest -v
```

For a source build, configure and run the CTest entry point:

```console
python3 -m pip install pytest
cmake -S . -B build/Release -DOPS_ENABLE_TESTS=ON
cmake --build build/Release --target OpenSeesPy
ctest --test-dir build/Release --output-on-failure
```

Sanitizer builds can be tested the same way by adding the compiler sanitizer
flags to the CMake configure command.

### Import Statements

So that the tests will work with the latest source code on GitHub 
Actions, first try a local import, then use the standard pip install if 
the local library is not found.

```python
try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
   
```

### Contributing

Add focused, deterministic test scripts to this folder via a PR. Prefer
physical or analytical checks with tolerances over exact text snapshots.
