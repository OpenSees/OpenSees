from conans import ConanFile, CMake, tools
import os

class OpenSeesDependencies(ConanFile):
    name = "OpenSeesDependencies"
    version = "1.0.0"
    description = "rovides Software Packages needed to build OpenSees"
    license = "BSD 3-Clause"
    author = "fmk fmckenna@berkeley.edu"
    settings = {"os": None, "build_type": None, "compiler": None, "arch": ["x86_64"]}
    options = {"shared": [True, False]}
    default_options = {"mkl-static:threaded": False, "ipp-static:simcenter_backend": True}    
    generators = "cmake"
    build_policy = "missing"
    requires = "libmysqlclient/8.0.25", \
        "tcl/8.6.10"
    # Custom attributes for Bincrafters recipe conventions
    _source_subfolder = "source_subfolder"
    _build_subfolder = "build_subfolder"
    # Set short paths for Windows
    short_paths = True    
    scm = {
        "type": "git",  # Use "type": "svn", if local repo is managed using SVN
        "subfolder": _source_subfolder,
        "url": "auto",
        "revision": "auto"
    }


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
        self.copy(pattern="LICENSE", dst="licenses", src=self._source_subfolder)
        cmake = self.configure_cmake()
        cmake.install()
        self.copy("*", dst="bin", src=self._source_subfolder + "/applications")

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
