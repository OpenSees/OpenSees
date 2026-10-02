#!/usr/bin/env python3
"""
Setup script for OpenSeesPy.

Builds the OpenSeesPy native extension from source using CMake,
then packages it as a pip-installable Python package.

The native C extension is internally named "opensees" (PyInit_opensees)
and is installed directly as a top-level module (not inside a package).
"""

import os
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

import setuptools
from setuptools import setup
from setuptools.command.build_ext import build_ext

# ── Package metadata ──────────────────────────────────────────────────

HERE = Path(__file__).resolve().parent
BUILD_DIR = HERE / "build"

# Base version from the project (reported as 3.8.0 by the module)
VERSION = "3.8.0+cc6dof"


# ── Custom build command ──────────────────────────────────────────────

class CMakeBuild(build_ext):
    """Build the OpenSeesPy native library via CMake."""

    def run(self):
        self._ensure_cmake_configured()
        self._build_openseespy()
        # Copy the built .so to the setuptools build output path
        self._copy_so_to_build()

    # ─ helpers ────────────────────────────────────────────────────────

    def _cmake_args(self):
        """Return list of cmake configuration arguments."""
        args = [
            "-S", str(HERE),
            "-B", str(BUILD_DIR),
        ]
        # Use the same Python that is running pip
        args += [
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DPython_INCLUDE_DIRS={sysconfig.get_path('include')}",
            f"-DPython_LIBRARIES={sysconfig.get_config_var('LIBDIR') or ''}",
        ]
        # Set the final target to OpenSeesPy so we only build what we need
        args += ["-DOPS_FINAL_TARGET=OpenSeesPy"]
        return args

    def _ensure_cmake_configured(self):
        """Run cmake configure step if the build tree is not already configured."""
        cache_file = BUILD_DIR / "CMakeCache.txt"
        if cache_file.exists():
            self.announce(
                "CMake build directory already configured, skipping configure step.",
                level=3,
            )
            return

        self.announce("Configuring CMake build ...", level=3)
        cmake_args = self._cmake_args()
        subprocess.check_call(["cmake"] + cmake_args, cwd=HERE)

    def _build_openseespy(self):
        """Build the OpenSeesPy target."""
        self.announce("Building OpenSeesPy target ...", level=3)
        nproc = os.cpu_count() or 1
        subprocess.check_call(
            ["cmake", "--build", str(BUILD_DIR), "--target", "OpenSeesPy",
             "--", f"-j{nproc}"],
            cwd=HERE,
        )

    def _find_so(self):
        """Return the path to the built OpenSeesPy shared library."""
        # Try opensees.so (symlink created by CMake install target)
        src = BUILD_DIR / "opensees.so"
        if src.exists():
            return src.resolve()
        # Fall back to the real library name
        src = BUILD_DIR / "OpenSeesPy.so"
        if src.exists():
            return src.resolve()
        # Search in the build directory
        for p in BUILD_DIR.rglob("OpenSeesPy*.so"):
            return p.resolve()
        raise RuntimeError(
            "Could not find built OpenSeesPy shared library in "
            f"{BUILD_DIR}. Build may have failed."
        )

    def _copy_so_to_build(self):
        """Copy the built .so to the setuptools build output path."""
        ext = self.extensions[0]
        dest_path = self.get_ext_fullpath(ext.name)
        src_path = self._find_so()

        self.announce(f"Copying {src_path} -> {dest_path}", level=3)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        shutil.copy2(src_path, dest_path)
        os.chmod(dest_path, 0o644)


# ── Setup ─────────────────────────────────────────────────────────────

setup(
    name="openseespy",
    version=VERSION,
    description="OpenSeesPy - Python interface for OpenSees",
    long_description=(
        "OpenSees (Open System for Earthquake Engineering Simulation) is a "
        "software framework for developing applications to simulate structural "
        "and geotechnical systems primarily in the field of earthquake "
        "engineering. This package provides the Python interface (OpenSeesPy)."
    ),
    long_description_content_type="text/plain",
    author="Pacific Earthquake Engineering Research Center (PEER)",
    author_email="opensees@berkeley.edu",
    url="https://opensees.berkeley.edu/",
    license="BSD 3-Clause",
    # No packages — the native extension is installed directly as a top-level module
    py_modules=[],  
    zip_safe=False,
    cmdclass={"build_ext": CMakeBuild},
    ext_modules=[
        setuptools.Extension(
            "opensees",
            sources=[],
        ),
    ],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: C++",
        "Programming Language :: Fortran",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)