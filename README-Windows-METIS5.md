# METIS 5 for OpenSeesMP

OpenSeesMP's Tcl `partition` command requires the METIS 5 headers and library.
The copy under `OTHER/METIS` is METIS 4 and is retained only for legacy graph
code. If METIS 5 is unavailable, OpenSeesMP can still build, but CMake warns
that `partition` will remain a stub.

## Ubuntu

Install the METIS 5 development package:

```bash
sudo apt update
sudo apt install libmetis-dev
```

The default `makeOpenSeesMP_Ubuntu.sh` configuration uses `/usr` as the METIS
prefix. For a nonstandard installation, set `METIS5_DIR`:

```bash
METIS5_DIR=/path/to/metis/install ./makeOpenSeesMP_Ubuntu.sh
```

The prefix must contain:

```text
include/metis.h
lib/libmetis.so (or the platform-specific library directory)
```

## Windows

There is no official Windows METIS binary. Build METIS 5.1.0 from source and
point OpenSees to its installation prefix. A convenient layout is:

```text
...\OpenSees
...\mumps\build          <- MUMPS_DIR
...\metis-5.1.0\install  <- METIS5_DIR
```

### Obtain METIS 5.1.0

Use the classic 5.1.0 release with bundled GKlib. This fork includes Windows
build support:

```bat
cd ..
git clone --depth 1 https://github.com/eric2003/METIS-5.1.0-Modified.git metis-5.1.0
cd metis-5.1.0
```

Do not substitute the KarypisLab `master` branch (METIS 5.2.x) for this build.
OpenSees's `partition` implementation targets the METIS 5.1 API, and
`METIS_PartMeshNodal` from 5.2.x has caused access violations with MSVC.

### Build with Visual Studio 2022

The OpenSees Windows script uses the static MSVC runtime (`/MT`). Build METIS
the same way. Keep 32-bit `idx_t` because OpenSeesMP communicates these arrays
using `MPI_INT`.

From the `metis-5.1.0` directory:

```bat
set CMAKE=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe
set PREFIX=%CD%\install

"%CMAKE%" -S . -B build -G "Visual Studio 17 2022" -A x64 ^
  -DMETIS_ENABLE_64BIT=OFF ^
  -DSHARED=FALSE ^
  -DCMAKE_INSTALL_PREFIX="%PREFIX%" ^
  -DCMAKE_POLICY_DEFAULT_CMP0091=NEW ^
  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded

"%CMAKE%" --build build --config Release --target metis --parallel 8
```

Only the `metis` library target is needed. Assemble the installation prefix:

```bat
mkdir install\include install\lib
copy include\metis.h install\include\
copy build\libmetis\Release\metis.lib install\lib\
```

The result should be:

```text
metis-5.1.0\install\include\metis.h
metis-5.1.0\install\lib\metis.lib
```

GKlib is included in `metis.lib`; it does not require a separate link step.

### Build OpenSeesMP

`makeOpenSeesMP_WIN.bat` defaults to the sibling directory shown above. Edit
`METIS5_DIR` near the top of the script if METIS is installed elsewhere, then
run:

```bat
makeOpenSeesMP_WIN.bat
```

The executable is written to:

```text
build-mp\Release\OpenSeesMP.exe
```

## Confirm METIS was detected

A successful CMake configuration includes messages similar to:

```text
-- METIS 5 (partition / OpenSeesMP): ... HAVE_METIS5_HEADER=TRUE
-- MPI executables: linking METIS 5 at .../libmetis...
```

If the configuration instead reports:

```text
METIS 5 headers not found; OpenSeesMP partition() will remain a stub
```

check `METIS5_DIR`, `include/metis.h`, and the library path. If the prefix was
changed after an earlier configuration, remove the stale CMake cache entry or
pass `-UOPENSEES_METIS5_LIBRARY`; the Windows build script already does this.

## Quick partition check

From `EXAMPLES/Partition`, use a smaller mesh for a fast check:

```bash
OpenSees Example.tcl serial 4 4 20
mpiexec -n 4 OpenSeesMP Example.tcl metis 4 4 20
mpiexec -n 4 OpenSeesMP Example.tcl custom 4 4 20
mpiexec -n 4 OpenSeesMP Example.tcl samePart 4 4 20
```

The gravity displacement, eigenvalues, and transient displacement should agree
between serial, METIS, custom, and `samePart` runs up to normal floating-point
roundoff.
