#!/usr/bin/env python3
"""Run OpenSeesPy regression tests against a freshly built module."""

import argparse
import importlib.machinery
import importlib.util
import os
from pathlib import Path
import sys

import pytest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", required=True, help="Path to the built OpenSeesPy module")
    parser.add_argument(
        "--tests-dir",
        action="append",
        required=True,
        help="Directory containing pytest tests; may be supplied more than once",
    )
    args = parser.parse_args()

    module = Path(args.module).resolve()
    if not module.is_file():
        parser.error(f"OpenSeesPy module does not exist: {module}")

    # Load the requested binary by its exported PyInit_opensees entry point;
    # never replace an existing alias or silently test a pip installation.
    dll_directory = os.add_dll_directory(str(module.parent)) if os.name == "nt" else None
    try:
        loader = importlib.machinery.ExtensionFileLoader("opensees", str(module))
        spec = importlib.util.spec_from_file_location("opensees", module, loader=loader)
        extension = importlib.util.module_from_spec(spec)
        loader.exec_module(extension)
        sys.modules["opensees"] = extension
        print(f"Testing OpenSeesPy module: {module}", flush=True)
        return pytest.main(["-q", *args.tests_dir])
    finally:
        if dll_directory is not None:
            dll_directory.close()


if __name__ == "__main__":
    raise SystemExit(main())
