from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps
from conan.tools.files import copy, export_conandata_patches, get
import os

class OpenSeesDependencies(ConanFile):
    name = "OpenSeesDependencies"
    version = "1.0.0"
    description = "Provides Software Packages needed to build OpenSees"
    license = "BSD-3-Clause"
    author = "fmk <fmckenna@berkeley.edu>"
    
    # Standard Conan 2 settings
    settings = "os", "compiler", "build_type", "arch"

    options = {
        "shared": [True, False],
    }
    default_options = {
        "shared": False,
        "mkl-static:threaded": False,
        "ipp-static:simcenter_backend": True,
    }

    # Conan 2 requires explicit generators in layout()
    exports_sources = "CMakeLists.txt", "source_subfolder/*"

    # Dependencies (same as before)
    requires = (
        "hdf5/1.14.0",
        "tcl/8.6.11",
        "zlib/1.2.13",
        "eigen/3.4.0"
    )

    short_paths = True

    def layout(self):
        self.folders.source = "source_subfolder"
        self.folders.build = "build_subfolder"

    # Instead of SCM mode (removed in Conan 2), fetch explicitly if needed
    # def source(self):
    #     get(self, url="https://example.com/source.tar.gz", strip_root=True)

    def generate(self):
        """
        Conan 2 replacement for generators="cmake", "cmake_find_package"
        """
        tc = CMakeToolchain(self)
        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        # Copy license
        copy(self, "LICENSE",
             src=self.source_folder,
             dst=os.path.join(self.package_folder, "licenses"))

        cmake = CMake(self)
        cmake.install()

        # Copy applications/bin output
        apps_src = os.path.join(self.source_folder, "applications")
        apps_dst = os.path.join(self.package_folder, "bin")
        copy(self, "*", src=apps_src, dst=apps_dst)

    def package_info(self):
        """
        Conan 2 no longer uses env_info.
        Add binary folder to PATH via cpp_info.
        """
        self.cpp_info.bindirs.append("bin")
