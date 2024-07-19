# OpenSeesRT

## Features

- Use of Tcl API encourages users to design well-structured commands
  with semantic meaning (e.g. using lists when meaningful as opposed to
  sequences of arguments)

- **Direct access to components** Developing OpenSees components
  is made easier by providing a direct interface to core object
  types, allowing calls constructs like `material.getStress(0.002, commit=False)`

- **Structured problem representation**

  - Allows introspective add-ons ([Torsion]())
  - Better visualization tools, independent of core source code.

- Promotes stability of OpenSees core.

- Supports more Python versions (3.6+)

Additional minor features:

- Verbosity control
- new `invoke` Tcl command and Python constructs
- new `progress` command
- new `export` command
- new `=` command, fixes vexing operator precedence in `expr`
- improved `print` command
  - new `-registry` option
  - more reliable JSON printing
  - includes `MP_Constraint` information
- consolidated `rigidLink` command
  - Fewer files
  - Better error handling (errors dont happen in a constructor, so they can
    be handled properly).


## User Changes

When `OpenSeesRT` is loaded as a Tcl library, there are a few minor
changes from the classic `OpenSees` interpreter:

- `puts` command prints to `stdout` by default, whereas classic OpenSees
  writes only to `stderr`. Errors (ie `opserr`) still writes to `stderr`.

- Removed `SimulationInformation` functionality

- Dropped section functions:

      // extern void *OPS_WFSection2d(G3_Runtime*);
      // extern void *OPS_RCCircularSection(G3_Runtime*);
      // extern void *OPS_RCSection2d(G3_Runtime*);
      // extern void *OPS_RCTBeamSection2d(G3_Runtime*);
      // extern void *OPS_RCTunnelSection(G3_Runtime*);
      // extern void *OPS_TubeSection(G3_Runtime*);




## Developer Changes

- No more `OPS_GetInt(void)`; use your host's API

- `ModelBuilder` namespacing functionality

  - Eliminates random code in important places like `Domain::Print`

- Load classes no longer derive from `DomainComponent`;

- Several stale classes removed:
  - `utility/SimulationInformation.*`
  - `utility/StringContainer.*`
  - removed all use of ad-hoc `bool.h`

## Proposed changes

- remove `sdfResponse` command; see [`sdof`](https://pypi.org/project/sdof) python package.
- remove `database` command.

## Codebase changes

- `socket.c` removed
- `utility/PeerNGA.cpp` removed
- Put `RigidDiaphragm` and friends into one file,
  change to functions.

## Cleaning & TODO

Remove dependence on

- Rendering commands to retire:
  ```
     "prp"
     "vrp"
     "vup"
     "vpn"
     "viewWindow"
     "plane"
     "port"
     "projection"
     "fill"
     "display"
  ```

- Remove TimeSeriesIntegrators from C++; handle in pre-processing?

Note to self:

  cmake -S /home/claudio/packages/opensees-pypi -B /home/claudio/packages/opensees-pypi/build/temp.linux-x86_64-cpython-39_rt -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=/home/claudio/packages/opensees-pypi/build/lib.linux-x86_64-cpython-39/opensees -G "Unix Makefiles" -DDependencies=Conda -DCMAKE_BUILD_TYPE=DEBUG -DPYTHON_EXECUTABLE:FILEPATH=/home/claudio/mambaforge/envs/py39/bin/python

  CC="clang" CXX="clang++" cmake -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE=include-what-you-use ..


