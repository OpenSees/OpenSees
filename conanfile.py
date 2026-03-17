# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "conan",
# ]
# ///

if __name__ == "__main__":
    from conan.cli.cli import main
    from pathlib import Path
    from sys import argv

    this_dir = Path(__file__).parent.relative_to(Path.cwd())

    args = argv[1:]
    if args:
        main(args)
    else:
        try:
            main(["profile", "detect", "--exist-ok"])
        except SystemExit as e:
            if e.code != 0:
                raise

        main(["build", str(this_dir), "--build=missing"])


import json
import shutil
import subprocess
from pathlib import Path
from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
from conan.tools.env.environment import create_env_script
from conan.tools.microsoft import VCVars, is_msvc

_DEFAULT_GENERATOR = "Ninja"

class OpenSees(ConanFile):
    settings = "os", "arch", "compiler", "build_type"

    def requirements(self):
        self.requires("hdf5/1.14.6")
        self.requires("tcl/8.6.13")

    def build_requirements(self):

        def missing_tool(tool_name, major_version, minor_version=0):
            exe = shutil.which(tool_name)
            if not exe:
                return True

            try:
                # Capture stderr as well, as some tools might print version there
                output = subprocess.check_output(
                    [exe, "--version"], 
                    stderr=subprocess.STDOUT, 
                    text=True
                )
                
                import re
                match = re.search(r"(\d+)\.(\d+)", output)
                if match:
                    major, minor = match.groups()
                    if int(major) != int(major_version):
                        return True
                    if int(minor) < int(minor_version):
                        return True
                    return False
                else:
                    return True # Can't parse version, assume incompatible
            except (subprocess.CalledProcessError, FileNotFoundError):
                return True

        if missing_tool("cmake", 4):
            self.tool_requires("cmake/4.1.1")
        
        if self.settings.os == "Linux":
            if missing_tool("ninja", 1, 11):
                self.tool_requires("ninja/1.13.2")
            
            if missing_tool("patchelf", 0, 16):
                self.tool_requires("patchelf/0.18")

    def layout(self):
        cmake_layout(self, generator=_DEFAULT_GENERATOR)

    def generate(self):
        tc = CMakeToolchain(self, generator=_DEFAULT_GENERATOR)

        if self.settings.os == "Windows":
            if tc.generator.startswith("Visual Studio"):
                # Sets CMAKE_GENERATOR_TOOLSET
                tc.blocks["generic_system"].values["toolset"] += ",fortran=ifx"
            else:
                # For non-Visual Studio generators, we need to set the Fortran compiler explicitly
                tc.variables["CMAKE_Fortran_COMPILER"] = "ifx"
                tc.variables["CMAKE_NINJA_FORCE_RESPONSE_FILE"] = "TRUE"
                tc.variables["CMAKE_Fortran_USE_RESPONSE_FILE_FOR_OBJECTS"] = "TRUE"
                tc.variables["CMAKE_Fortran_USE_RESPONSE_FILE_FOR_LIBRARIES"] = "TRUE"
                tc.variables["CMAKE_Fortran_USE_RESPONSE_FILE_FOR_INCLUDES"] = "TRUE"

        # This stops extra configurations from being generated
        # See https://github.com/conan-io/docs/issues/3517
        if tc.is_multi_configuration:
            configuration_types = set([self.settings.get_safe("build_type")])

            preset_path = Path(self.generators_folder) / "CMakePresets.json"
            if preset_path.exists():
                preset = json.loads(preset_path.read_text())
                # Conan generated presets should have only 1 configurePreset, no more
                cache_variables = preset["configurePresets"][0]["cacheVariables"]
                if existing := cache_variables.get("CMAKE_CONFIGURATION_TYPES"):
                    configuration_types.update(existing.split(";"))

            configuration_types_str = ";".join(sorted(configuration_types))
            tc.cache_variables["CMAKE_CONFIGURATION_TYPES"] = configuration_types_str

        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()

        if is_msvc(self):
            # This is done before oneAPI env script
            # to ensure the correct settings are used,
            # e.g., windows SDK version
            # If not called beforehand,
            # oneapi-vars.bat will call it with default values
            VCVars(self).generate()

            oneapi_vars_bat = r"C:\Program Files (x86)\Intel\oneAPI\2025.0\oneapi-vars.bat"
            create_env_script(self, f'call "{oneapi_vars_bat}"', "oneapi_vars.bat")

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        config_build_folder = Path(self.build_folder)
        if cmake.is_multi_configuration:
            config_build_folder /= str(self.settings.build_type)
        install_dir = config_build_folder / "install"
        if install_dir.exists():
            shutil.rmtree(install_dir)
        cmake.install(cli_args=["--prefix", str(install_dir)])
