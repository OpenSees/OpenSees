"""Binary recorder storage tests; set OPENSEES_TCL to also test the Tcl executable."""

import os
import struct
import subprocess

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


KINDS = ("Node", "EnvelopeNode", "NodeRMS", "Element", "EnvelopeElement", "ElementRMS")
PYTHON_KINDS = tuple(kind for kind in KINDS if not kind.endswith("RMS"))
MODES = {
    "default": (),
    "float64": ("-float64",),
    "float32": ("-float32",),
    "before": ("-float32",),
    "reset": ("-float32", "-float64"),
}


def recorder_args(kind):
    if kind in ("Drift", "EnvelopeDrift"):
        return ("-time", "-iNode", 1, "-jNode", 2, "-dof", 2, "-perpDirn", 2)
    time = () if kind.endswith("RMS") else ("-time",)
    if "Node" in kind:
        return time + ("-node", 2, "-dof", 2, "disp")
    return time + ("-ele", 1, "axialForce")


def output_args(directory, mode):
    binary = ("-binary", str(directory / (mode + ".bin")))
    return MODES[mode] + binary if mode == "before" else binary + MODES[mode]


def read_rows(path, columns, code):
    raw = path.read_bytes()
    row = struct.Struct("=" + code * columns)
    stride = row.size + 1
    assert raw and len(raw) % stride == 0
    result = []
    for offset in range(0, len(raw), stride):
        assert raw[offset + row.size:offset + stride] == b"\n"
        result.append(row.unpack_from(raw, offset))
    return result


def run_python(kind, directory):
    ops.wipe()
    try:
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, 1.0)
        ops.fix(1, 1, 1)
        ops.fix(2, 1, 0)
        ops.uniaxialMaterial("Elastic", 1, 10.0)
        ops.element("truss", 1, 1, 2, 1.0, 1)
        for mode in MODES:
            assert ops.recorder(kind, *output_args(directory, mode), *recorder_args(kind)) >= 0
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0.0, 1.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 0.125)
        ops.analysis("Static")
        assert ops.analyze(3) == 0
    finally:
        ops.wipe()


def run_tcl(kind, directory):
    executable = os.environ.get("OPENSEES_TCL")
    if not executable:
        pytest.skip("set OPENSEES_TCL to test the Tcl executable")
    commands = []
    for flag in ("-float32", "-float64"):
        args = (kind, "-file", str(directory / "invalid.txt"), flag, *recorder_args(kind))
        command = "recorder " + " ".join("{" + str(arg).replace("\\", "/") + "}" for arg in args)
        commands.append("if {![catch {" + command + "}]} {error {expected storage flag rejection}}")
    for mode in MODES:
        args = (kind, *output_args(directory, mode), *recorder_args(kind))
        commands.append("recorder " + " ".join("{" + str(arg).replace("\\", "/") + "}" for arg in args))
    script = directory / "record.tcl"
    script.write_text("""
wipe
model basic -ndm 2 -ndf 2
node 1 0 0
node 2 0 1
fix 1 1 1
fix 2 1 0
uniaxialMaterial Elastic 1 10
element truss 1 1 2 1 1
""" + "\n".join(commands) + """
timeSeries Linear 1
pattern Plain 1 1 {load 2 0 1}
constraints Plain
numberer Plain
system BandGeneral
algorithm Linear
integrator LoadControl 0.125
analysis Static
if {[analyze 3] != 0} {exit 1}
wipe
puts BINARY_TEST_OK
exit
""")
    result = subprocess.run([executable, str(script)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "BINARY_TEST_OK" in result.stdout + result.stderr, result.stdout + result.stderr


@pytest.mark.parametrize("backend,kind", [("python", k) for k in PYTHON_KINDS]
                         + [("tcl", k) for k in (*KINDS, "Drift", "EnvelopeDrift")])
def test_binary_storage(tmp_path, backend, kind):
    (run_python if backend == "python" else run_tcl)(kind, tmp_path)
    columns = 1 if kind.endswith("RMS") else 2
    reference = read_rows(tmp_path / "default.bin", columns, "d")
    assert (tmp_path / "float64.bin").read_bytes() == (tmp_path / "default.bin").read_bytes()
    assert (tmp_path / "reset.bin").read_bytes() == (tmp_path / "default.bin").read_bytes()
    expected = [tuple(struct.unpack("=f", struct.pack("=f", value))[0] for value in row)
                for row in reference]
    for mode in ("float32", "before"):
        assert read_rows(tmp_path / (mode + ".bin"), columns, "f") == expected
    if kind == "Node":
        assert [row[0] for row in reference] == [0.125, 0.25, 0.375]
        assert [row[1] for row in reference] == pytest.approx([0.0125, 0.025, 0.0375])


@pytest.mark.parametrize("kind", PYTHON_KINDS)
@pytest.mark.parametrize("flag", ("-float32", "-float64"))
def test_storage_flag_requires_binary(tmp_path, kind, flag):
    ops.wipe()
    try:
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        with pytest.raises(Exception):
            ops.recorder(kind, "-file", str(tmp_path / "out.txt"), flag, *recorder_args(kind))
    finally:
        ops.wipe()
