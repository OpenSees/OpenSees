#!/usr/bin/env python3
"""Run pytest using a specific locally built OpenSeesPy extension."""

import argparse
import importlib.machinery
import importlib.util
from pathlib import Path
import sys

import pytest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", type=Path, required=True)
    parser.add_argument("--tests-dir", action="append", required=True)
    args = parser.parse_args()
    module = args.module.resolve()
    if not module.is_file():
        parser.error(f"OpenSeesPy module does not exist: {module}")

    # CMake calls the library OpenSeesPy, but its entry point is PyInit_opensees.
    loader = importlib.machinery.ExtensionFileLoader("opensees", str(module))
    spec = importlib.util.spec_from_file_location("opensees", module, loader=loader)
    extension = importlib.util.module_from_spec(spec)
    loader.exec_module(extension)
    sys.modules["opensees"] = extension
    print(f"Testing OpenSeesPy module: {module}", flush=True)
    return pytest.main(["-q", *args.tests_dir])


if __name__ == "__main__":
    raise SystemExit(main())
