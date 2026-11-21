"""printA -file <fileName> -sparse <baseIndex> -binary

writes the matrix of the system of equations to a binary file, in the byte order of the machine:

    int64   number of rows, number of columns, number of entries (nnz), baseIndex
    int32   row index of each entry     (nnz of them)
    int32   column index of each entry  (nnz of them)
    double  value of each entry         (nnz of them)

The file must hold exactly what printA -sparse <baseIndex> -ret returns.
"""

import array
import struct

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _portal_frame(system):
    """A portal frame of three elastic beams, pushed sideways, analysed with the given system."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for tag, x, y in ((1, 0.0, 0.0), (2, 0.0, 3.0), (3, 4.0, 3.0), (4, 4.0, 0.0)):
        ops.node(tag, x, y)
    ops.fix(1, 1, 1, 1)
    ops.fix(4, 1, 1, 0)
    ops.geomTransf("Linear", 1)
    for tag, i, j in ((1, 1, 2), (2, 2, 3), (3, 3, 4)):
        ops.element("elasticBeamColumn", tag, i, j, 0.01, 2.0e11, 1.0e-4, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1000.0, 0.0, 0.0)
    ops.system(system)
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _read_binary(path):
    """The header and the three arrays of a file written by printA -binary."""
    data = path.read_bytes()
    rows, cols, nnz, base = struct.unpack_from("=4q", data, 0)
    offset = struct.calcsize("=4q")
    arrays = []
    for code in ("i", "i", "d"):
        a = array.array(code)
        assert a.itemsize == (4 if code == "i" else 8)
        a.frombytes(data[offset:offset + nnz * a.itemsize])
        arrays.append(list(a))
        offset += nnz * a.itemsize
    assert offset == len(data)
    return (rows, cols, nnz, base), arrays[0], arrays[1], arrays[2]


@pytest.mark.parametrize("system", ["SparseGeneral", "UmfPack", "FullGeneral"])
@pytest.mark.parametrize("base", [0, 1])
def test_printA_binary_holds_what_sparse_ret_returns(tmp_path, system, base):
    path = tmp_path / "A.bin"
    try:
        _portal_frame(system)
        assert ops.printA("-file", str(path), "-sparse", base, "-binary") == 0
        expected = ops.printA("-sparse", base, "-ret")
        size = ops.systemSize()
    finally:
        ops.wipe()

    header, rows, cols, values = _read_binary(path)
    assert header == (size, size, len(expected["values"]), base)
    assert rows == list(expected["rowIndices"])
    assert cols == list(expected["colIndices"])
    # exactly the same doubles: nothing is rounded on the way through the file
    assert values == list(expected["values"])


def test_printA_binary_holds_the_matrix(tmp_path):
    path = tmp_path / "A.bin"
    try:
        _portal_frame("FullGeneral")
        ops.printA("-file", str(path), "-sparse", 1, "-binary")
        dense = ops.printA("-ret")
        size = ops.systemSize()
    finally:
        ops.wipe()

    (n, _, _, base), rows, cols, values = _read_binary(path)
    assert n == size and size > 0
    full = [[0.0] * n for _ in range(n)]
    for i, j, v in zip(rows, cols, values):
        full[i - base][j - base] = v
    # printA -ret returns the dense matrix as it is stored, column by column
    for i in range(n):
        for j in range(n):
            assert full[i][j] == dense[j * n + i]


def test_printA_binary_matches_the_text_file(tmp_path):
    text, binary = tmp_path / "A.mtx", tmp_path / "A.bin"
    try:
        _portal_frame("SparseGeneral")
        ops.printA("-file", str(text), "-sparse", 1, "-precision", 16)
        ops.printA("-file", str(binary), "-sparse", 1, "-binary")
    finally:
        ops.wipe()

    lines = [line.split() for line in text.read_text().splitlines()
             if line.strip() and not line.startswith("%")]
    header, rows, cols, values = _read_binary(binary)
    assert header == (int(lines[0][0]), int(lines[0][1]), int(lines[0][2]), 1)
    entries = lines[1:]
    assert rows == [int(e[0]) for e in entries]
    assert cols == [int(e[1]) for e in entries]
    # the text holds 16 significant digits
    assert values == pytest.approx([float(e[2]) for e in entries], rel=1e-15, abs=0.0)


@pytest.mark.parametrize("args", [
    ("-binary",),
    ("-file", "A.bin", "-binary"),
    ("-file", "A.bin", "-sparse", 0, "-binary", "-ret"),
])
def test_printA_binary_needs_a_file_and_sparse_and_no_ret(tmp_path, args):
    path = tmp_path / "A.bin"
    args = tuple(str(path) if a == "A.bin" else a for a in args)
    try:
        _portal_frame("SparseGeneral")
        with pytest.raises(Exception):
            ops.printA(*args)
    finally:
        ops.wipe()
    assert not path.exists()


def test_printA_binary_needs_a_system_of_equations(tmp_path):
    path = tmp_path / "A.bin"
    ops.wipe()
    try:
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0)
        with pytest.raises(Exception):
            ops.printA("-file", str(path), "-sparse", 0, "-binary")
    finally:
        ops.wipe()
    assert not path.exists()
