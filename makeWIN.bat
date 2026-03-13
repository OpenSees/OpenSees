#rd /s /q build

conan install . -s arch=x86_64 -s compiler.runtime=static --build=missing -c tools.cmake.cmaketoolchain:generator=Ninja  

"C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"  -S . -B build/Release -G "Ninja" -DCMAKE_TOOLCHAIN_FILE=build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER="C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/ifx.exe" -DBLA_STATIC=ON -DMKL_LINK=static -DMKL_INTERFACE_FULL=intel_lp64 -DMUMPS_DIR="..\..\mumps\build" -DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON -DCMAKE_C_USE_RESPONSE_FILE_FOR_OBJECTS=ON -DCMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS=ON -DCMAKE_INSTALL_PREFIX=%USERPROFILE%\bin\OpenSees3.8.0

cd build\Release
cmake --build . --target OpenSees --config Release --parallel 10
cmake --build . --target OpenSeesPy --config Release
cmake --install .