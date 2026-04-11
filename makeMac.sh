rm -fr build
conan install . --build missing
cmake -S . -B build/Release -DCMAKE_TOOLCHAIN_FILE=build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/bin/OpenSees3.8.0
cd build/Release
cmake --build . --parallel 10
cmake --build . --target OpenSeesPy
cmake --install .
