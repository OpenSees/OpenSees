
from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
from conan.tools.files import copy
from conan.errors import ConanInvalidConfiguration
import os

class OpenSeesDependencies(ConanFile):
    name = "OpenSeesDependencies"
    version = "1.0.0"
    description = "Provides Software Packages needed to build OpenSees"
    license = "BSD 3-Clause"
    author = "fmk fmckenna@berkeley.edu"

    # 1. Settings: Explicitly define the settings needed for compilation
    settings = "os", "compiler", "build_type", "arch"

    options = {}
    default_options = {}

    exports_sources = "CMakeLists.txt", "src/*", "applications/*", "LICENSE"

    def layout(self):
        # 5. Layout: This standardizes where source/build files go. 
        # It replaces the manual '_source_subfolder' logic.
        cmake_layout(self)

    def requirements(self):
        self.requires("hdf5/1.14.0")
        self.requires("tcl/8.6.11")
        self.requires("zlib/1.3.1")
        self.requires("eigen/3.4.0")
        # Note: If you rely on 'mkl-static', you must add self.requires("mkl-static/...") here.

    def generate(self):
        # 6. Generators: Standard toolchain and dependency generator
        tc = CMakeToolchain(self)
        tc.generate()
        
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        # 7. Packaging: Use 'self.source_folder' which is auto-set by cmake_layout()
        copy(self, "LICENSE", src=self.source_folder, dst=os.path.join(self.package_folder, "licenses"))
        
        cmake = CMake(self)
        cmake.install()
        
        # Copy extra applications if CMake install didn't catch them
        copy(self, "*", 
             src=os.path.join(self.source_folder, "applications"), 
             dst=os.path.join(self.package_folder, "bin"))

    def package_info(self):
        # 8. RunEnv: Correctly add the bin folder to PATH for consumers
        self.runenv_info.append_path("PATH", os.path.join(self.package_folder, "bin"))
