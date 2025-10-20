## Test Scripts

A directory of verification tests using OpenSeesPy.

### How to Run

Run locally with `pytest`

```console
pytest -v
```

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
