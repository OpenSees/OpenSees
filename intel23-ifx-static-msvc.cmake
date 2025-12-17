# ------------------------------------------------------------
# Intel oneAPI ifx + MSVC + static runtime toolchain
# ------------------------------------------------------------

# We are targeting Windows
set(CMAKE_SYSTEM_NAME Windows)

# --- COMPILERS ---
# IMPORTANT: use FULL PATHS
set(CMAKE_Fortran_COMPILER "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/ifx.exe" CACHE FILEPATH "")
set(CMAKE_C_COMPILER       "cl.exe" CACHE STRING "")
set(CMAKE_CXX_COMPILER     "cl.exe" CACHE STRING "")

# --- MSVC RUNTIME: STATIC ---
# This is the key line to avoid multithreaded DLLs
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>" CACHE STRING "")

# --- FORTRAN FLAGS ---
# Ensure static runtime for ifx
set(CMAKE_Fortran_FLAGS_INIT "/libs:static" CACHE STRING "")

# --- OPTIONAL: avoid auto-switching generators ---
set(CMAKE_GENERATOR_TOOLSET "" CACHE STRING "")

# ------------------------------------------------------------