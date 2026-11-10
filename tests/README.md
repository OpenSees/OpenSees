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
cmake --build build/Release --target OpenSeesPy matrix-vector-otherfact-test
ctest --test-dir build/Release --output-on-failure
```

CTest also runs the native Matrix/Vector regression when `BUILD_TESTING=ON`,
independently of `OPS_ENABLE_TESTS`. Build both targets above before running
the combined suite. Set `BUILD_TESTING=OFF` to disable all CTest registration.

Sanitizer builds can be tested the same way by adding the compiler sanitizer
flags to the CMake configure command. When a GNU AddressSanitizer-instrumented
module is loaded into an ordinary Python executable, preload the matching
runtime before running CTest:

```console
LD_PRELOAD="$(gcc -print-file-name=libasan.so)" ctest --test-dir build/Sanitize --output-on-failure
```

The runner imports the exact binary passed to `--module`, without creating,
overwriting, or deleting an `opensees.so`/`opensees.pyd` alias. It does not fall
back to a pip installation. On Windows, dependent DLLs must be present beside
the built module before the tests run.

### Numerical Coverage

`test_numerical_benchmarks.py` checks cantilever displacement, reactions and
strain energy against closed-form results, lumped and consistent truss mass,
and single-DOF frequency and undamped total-energy conservation. The beam load
and cable mass fixes carry their own focused tests, which this same directory
discovery will execute when those fixes are merged; this CI PR does not depend
on either numerical fix being accepted first.

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
