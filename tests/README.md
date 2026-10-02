## Test Scripts

A directory of verification tests using OpenSeesPy.

### How to Run

Run locally with `pytest`

```console
pytest -v
```

To test a locally compiled extension, use the Python interpreter it was built
against. This loads the specified binary and fails if it cannot be loaded:

```console
python tests/run_pytest.py --module build/Release/OpenSeesPy.dylib --tests-dir tests/test_catenary_cable.py
```

Use the actual build path and extension suffix for your platform. The
CatenaryCable tests cover the documented verification model, 3-/6-DOF
equivalence, self-weight equilibrium, translational mass, mixed-DOF elements,
force recorder metadata, and coupling to a beam with free rotations.

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

Add test scripts to this folder via a PR.
