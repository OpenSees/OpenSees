Binary recorder floating-point storage
=====================================

Use `-float32` with `-binary` to store recorder values as 32-bit floats:

```tcl
recorder Node -binary displacement.bin -float32 -time -node 2 -dof 1 disp
recorder Element -binary forces.bin -float32 -time -ele 1 localForce
```

```python
ops.recorder('Node', '-binary', 'displacement.bin', '-float32',
             '-time', '-node', 2, '-dof', 1, 'disp')
```

The default remains 64-bit doubles. `-float64` explicitly selects the default.
Place storage flags before the response name; they may precede or follow
`-binary filename`. If repeated, the last storage flag wins. Storage flags
require binary output.

The option applies to node, element, envelope, RMS, and drift recorders that
support binary output. It changes output storage only; analysis calculations
remain double precision. All numeric columns, including `-time`, use the
selected width. Finite values outside the float32 range produce an error
message and stop writes to that stream.

Files retain the existing native byte order and one newline byte after each
row, with no header or type metadata. Readers must know the column count and
storage width and skip the newline byte after each row. The existing
`binaryToText` converter expects 64-bit files.

Run the regression tests with `pytest tests/test_binary_recorder.py` against a
build containing this feature. Set `OPENSEES_TCL` to the Tcl executable path
to include Tcl recorder tests.
