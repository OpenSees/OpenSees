"""Check exact-module loading without touching existing module aliases."""

import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest


def _runner():
    path = Path(__file__).with_name("run_pytest.py")
    spec = importlib.util.spec_from_file_location("regression_runner", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


@pytest.mark.parametrize("windows", [False, True])
def test_exact_binary_alias_preservation_and_exit_code(tmp_path, monkeypatch, windows):
    runner = _runner()
    binary = tmp_path / ("OpenSeesPy.dll" if windows else "OpenSeesPy.so")
    binary.write_bytes(b"mock binary")
    alias = tmp_path / ("opensees.pyd" if windows else "opensees.so")
    alias.write_bytes(b"existing alias must survive")
    calls = []
    extension = SimpleNamespace()
    loader = SimpleNamespace(exec_module=lambda module: calls.append(("exec", module)))
    dll_handle = SimpleNamespace(close=lambda: calls.append(("close",)))
    monkeypatch.setattr(runner, "os", SimpleNamespace(
        name="nt" if windows else "posix",
        add_dll_directory=lambda path: calls.append(("dll", path)) or dll_handle))
    monkeypatch.setattr(runner.importlib.machinery, "ExtensionFileLoader",
                        lambda name, path: calls.append(("load", name, path)) or loader)
    monkeypatch.setattr(runner.importlib.util, "spec_from_file_location", lambda *a, **kw: None)
    monkeypatch.setattr(runner.importlib.util, "module_from_spec", lambda spec: extension)
    monkeypatch.setattr(runner.pytest, "main", lambda args: calls.append(("pytest", args)) or 5)
    monkeypatch.setattr(sys, "argv", ["run_pytest.py", "--module", str(binary),
                                    "--tests-dir", "first", "--tests-dir", "second"])
    monkeypatch.setitem(sys.modules, "opensees", sys.modules.get("opensees"))
    assert runner.main() == 5
    assert ("load", "opensees", str(binary.resolve())) in calls
    assert ("exec", extension) in calls
    assert ("pytest", ["-q", "first", "second"]) in calls
    assert binary.read_bytes() == b"mock binary"
    assert alias.read_bytes() == b"existing alias must survive"
    assert (("close",) in calls) == windows


def test_missing_binary_is_an_error(tmp_path, monkeypatch):
    runner = _runner()
    monkeypatch.setattr(sys, "argv", ["run_pytest.py", "--module", str(tmp_path / "missing.so"),
                                    "--tests-dir", "."])
    with pytest.raises(SystemExit) as error:
        runner.main()
    assert error.value.code == 2
