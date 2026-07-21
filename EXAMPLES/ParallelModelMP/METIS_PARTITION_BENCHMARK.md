# OpenSeesMP METIS partition benchmark

The benchmark constructs a 40-by-20 two-dimensional truss grid containing
1,660 elements. Every MPI rank constructs the same model before the `partition`
command distributes its elements.

## Environment

- macOS 26.5.2, Apple Silicon
- OpenSees 3.8.0 (`098e28c`)
- MPICH 5.0.1
- METIS 5.1.0
- Tcl 8.6.18
- Release build targeting macOS 26.0

## Results

### macOS / MPICH

| Ranks | Objective | METIS value | Elements per rank | Sum | Maximum imbalance |
|---:|---|---:|---|---:|---:|
| 2 | volume | 42 | 832, 828 | 1,660 | 0.24% |
| 4 | volume | 112 | 410, 417, 416, 417 | 1,660 | 0.48% |
| 8 | volume | 186 | 209, 211, 208, 211, 211, 202, 211, 197 | 1,660 | 1.69% |
| 4 | cut | 65 | 408, 411, 421, 420 | 1,660 | 1.45% |

The four-rank volume run was repeated with seed 42 and produced the same
objective value and per-rank element counts. An invalid objective was rejected
on every rank with a Tcl error.

### Ubuntu / Open MPI

A clean Ubuntu 24.04 ARM64 build linked against Open MPI 4.1.6 and system
METIS 5.1.0 passed all seven pytest cases. Its four-rank volume reference
run produced objective value 113 and element counts 420, 419, 410, and 411
(1.20% maximum imbalance).

### Ubuntu / MPICH

An independent Ubuntu build linked against MPICH 4.2.0 also passed all seven
cases. Ubuntu's ARM64 Hydra/PMIx package combination created singleton
communicators even for a minimal MPI program, so that binary was launched by
Open MPI's PMIx-aware launcher. The MPI implementation loaded by OpenSeesMP
remained MPICH.

## Partition-state regression

`metisPartitionState.tcl` adds a four-rank structural model that verifies
partitioned state beyond element counts. It checks that element loads follow
their owning elements, nodal mass remains on exactly one rank, floating fixed
nodes retain their constraints and nodal loads, and mesh-to-floating plus
floating-to-floating `equalDOF` constraints remain complete.

## Run

```sh
export TCL_LIBRARY=/opt/homebrew/opt/tcl-tk@8/lib/tcl8.6
mpiexec -n 4 OpenSeesMP metisPartitionBenchmark.tcl volume 42
mpiexec -n 4 OpenSeesMP metisPartitionBenchmark.tcl cut 42
```

OpenSeesMP now synchronizes ranks before `MPI_Finalize`, and all benchmark runs
exit cleanly after printing their final results.

The pytest launcher and executable can be selected without editing the test:

```sh
OPENSEESMP=/path/to/OpenSeesMP \
MPIEXEC=/path/to/mpiexec \
MPIEXEC_FLAGS="--oversubscribe" \
python -m pytest -q tests/test_metis_partition_benchmark.py
```
