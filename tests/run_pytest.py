#!/usr/bin/env python3
"""Run OpenSeesPy regression tests against a freshly built module."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


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

    # The extension exports PyInit_opensees, while the CMake target is named
    # OpenSeesPy. Keep the alias beside the built module so its dependencies
    # resolve consistently on Linux, macOS, and Windows.
    alias_name = "opensees.pyd" if os.name == "nt" else "opensees.so"
    alias = module.parent / alias_name
    shutil.copy2(module, alias)

    env = os.environ.copy()
    python_path = str(module.parent)
    if env.get("PYTHONPATH"):
        python_path += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = python_path

    try:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", *args.tests_dir],
            cwd=module.parent,
            env=env,
            check=False,
        )
        return result.returncode
    finally:
        alias.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
