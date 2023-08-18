import os
from conan import ConanFile
from conan.tools.files import save, load, copy, update_conandata
from conan.tools.cmake import CMake, CMakeToolchain,CMakeDeps, cmake_layout
from conan.tools.gnu import AutotoolsToolchain, AutotoolsDeps
from conan.tools.microsoft import unix_path, VCVars, is_msvc
from conan.errors import ConanInvalidConfiguration
from conan.errors import ConanException
from conan.tools.intel import IntelCC


class OpenSeesDependencies(ConanFile):
    name = "OpenSeesDependencies"
    version = "2.0.0"
    description = "Provides Software Packages needed to build OpenSees"
    license = "BSD 3-Clause"
    author = "fmk fmckenna@berkeley.edu"
    settings = {"os", "build_type", "compiler", "arch"}
    options = {"shared": [True, False]}
    default_options = {"mkl-static/*:threaded": False, "ipp-static/*:simcenter_backend": True}    
    generators = "CMakeToolchain","CMakeDeps"#,"IntelCC"
    user_presets_path = 'CMakeUserPresets.json'
    build_policy = "missing"
    requires = "hdf5/1.12.2","tcl/8.6.11"
    # Custom attributes for Bincrafters recipe conventions
    _source_subfolder = "source_subfolder"
    _build_subfolder = "build_subfolder"
    # Set short paths for Windows
    short_paths = True    
    
    def export(self):
        git = Git(self, self.recipe_folder)
        scm_url, scm_commit = git.get_url_and_commit()
        update_conandata(self, {"sources": {"commit": scm_commit, "url": scm_url}})

    def source(self):
        git = Git(self)
        sources = self.conan_data["sources"]
        git.clone(url=sources["url"], target=".")
        git.checkout(commit=sources["commit"])
    #scm = {
    #    "type": "git",  # Use "type": "svn", if local repo is managed using SVN
    #    "subfolder": _source_subfolder,
    #    "url": "auto",
    #    "revision": "auto"
    #}
    
    def layout(self):
        cmake_layout(self)
        
    def configure(self):
        self.options.shared = False

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder=self._source_subfolder)
        return cmake
    
    def build(self):
        cmake = self.configure_cmake()
        cmake.build()

    def package(self):
        #self.copy(pattern="LICENSE", dst="licenses", src=self._source_subfolder)
        copy(self, pattern="LICENSE", src=self._source_subfolder, dst="licenses")
        cmake = self.configure_cmake()
        cmake.install()
        #self.copy("*", dst="bin", src=self._source_subfolder + "/applications")
        copy(self, "*", self._source_subfolder + "/applications", "bin")

    def package_info(self):
        #self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.buildenv_info.PATH.append(os.path.join(self.package_folder, "bin"))
