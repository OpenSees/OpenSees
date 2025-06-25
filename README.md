# OpenSees - Open System for Earthquake Engineering Simulation

This repository contains the OpenSees source code with earthquake modeling capabilities, Python analysis tools, and a modernized build system.

## Table of Contents
- [Overview](#overview)
- [Project Structure](#project-structure)
- [Earthquake Modeling](#earthquake-modeling)
- [Python Analysis Tools](#python-analysis-tools)
- [Build Instructions](#build-instructions)
- [Usage Examples](#usage-examples)
- [Documentation](#documentation)

## Overview

OpenSees is an open-source software framework for simulating the response of structural and geotechnical systems subjected to earthquakes. This repository includes:

- **Core OpenSees**: Main simulation engine with 2D/3D finite element capabilities
- **Earthquake Models**: 2D elastic soil models with absorbing boundary conditions
- **Python Analysis Tools**: Visualization and post-processing scripts
- **Modern Build System**: CMake-based build with Conan dependency management

## Project Structure

### Core Directories
- `SRC/` - Main OpenSees source code
- `OTHER/` - Third-party libraries (cleaned up)
- `EXAMPLES/` - Example files and tests
- `SCRIPTS/` - Build and utility scripts
- `cmake/` - CMake modules and utilities
- `docker/` - Docker configuration
- `.github/` - GitHub Actions CI/CD

### Key Files
- `earthquake_model.tcl` - 2D elastic soil earthquake model
- `eq_plot_results.py` - Python visualization script for earthquake results
- `cantilever_q4.tcl` - Cantilever beam example
- `cantilever_q4.py` - Python version of cantilever analysis
- `CMakeLists.txt` - Main CMake configuration
- `conanfile.py` - Conan dependency management

## Earthquake Modeling

### 2D Elastic Soil Model (`earthquake_model.tcl`)

A sophisticated 2D earthquake model featuring:

- **Mesh**: Configurable Q4 element mesh (default: 10×5 elements)
- **Material**: Linear elastic isotropic soil
- **Boundary Conditions**: Absorbing boundary conditions (viscous dashpots)
- **Loading**: Earthquake motion applied as shear stress at base nodes
- **Analysis**: Transient dynamic analysis with Rayleigh damping

#### Model Parameters
```tcl
set E_soil 3.0e7         ;# Young's Modulus (Pa)
set nu_soil 0.3          ;# Poisson's Ratio
set rho_soil 2000.0      ;# Density (kg/m³)
set accel_factor 50.0    ;# Earthquake acceleration amplitude
```

#### Running the Earthquake Model
```bash
# Set up environment
source /opt/intel/oneapi/setvars.sh
module load mpi/openmpi-x86_64

# Run the analysis
./build/bin/OpenSees earthquake_model.tcl
```

#### Output Files
All results are saved in `RESULTS/` directory:
- `eq_displacements.out` - Node displacement history
- `eq_accelerations.out` - Node acceleration history
- `eq_element_stress.out` - Element stress history
- `eq_nodes.txt` - Node coordinates
- `eq_elements.txt` - Element connectivity

### Visualization (`eq_plot_results.py`)

Comprehensive Python visualization script that generates:

1. **Displacement History**: Top node displacement vs time
2. **Mesh Comparison**: Original vs deformed mesh
3. **Stress Distribution**: Contour plots of σxx, σyy, τxy stresses
4. **Deformation Animation**: MP4 animation of mesh deformation

#### Running Visualization
```bash
python eq_plot_results.py
```

## Python Analysis Tools

### Cantilever Beam Analysis (`cantilever_q4.py`)

A Python implementation of 2D Q4 element cantilever beam analysis.

#### Model Description
- **Geometry**: 100" × 10" × 1" (Length × Height × Thickness)
- **Material**: Linear elastic isotropic (E = 30,000 ksi, ν = 0.3)
- **Mesh**: 20 × 4 Q4 elements
- **Boundary Conditions**: Fixed at left end
- **Loading**: 100 kips distributed load at right end

#### Prerequisites
```bash
# Create conda environment
conda create -n openS python=3.8
conda activate openS

# Install OpenSees (method depends on your setup)
```

#### Running the Analysis
```bash
conda activate openS
python cantilever_q4.py
```

#### Output Files
- `displacements.out` - Node displacements
- `elementForces.out` - Element forces
- `elementStrains.out` - Element strains
- `reactions.out` - Support reactions
- `stress_distribution.png` - Stress visualization

## Build Instructions

### Prerequisites
- CMake 3.16 or higher
- Conan package manager
- C++ compiler with C++17 support
- Intel MKL libraries (for optimal performance)

### Building with CMake (Recommended)
```bash
# Create build directory
mkdir build
cd build

# Install dependencies with Conan
conan install .. --build missing

# Configure and build
cmake ..
make -j$(nproc)
```

### Environment Setup
```bash
# Source Intel environment (for MKL support)
source /opt/intel/oneapi/setvars.sh

# Load MPI module
module load mpi/openmpi-x86_64
```

### Available Targets
- `OpenSees` - Main Tcl interpreter
- `OpenSeesPy` - Python module

## Usage Examples

### Basic Earthquake Analysis
```bash
# Run earthquake model
./build/bin/OpenSees earthquake_model.tcl

# Visualize results
python eq_plot_results.py
```

### Cantilever Beam Analysis
```bash
# Run cantilever analysis
python cantilever_q4.py
```

### Custom Model Development
```bash
# Create your own Tcl model
# Use the existing models as templates
# Run with OpenSees interpreter
./build/bin/OpenSees your_model.tcl
```

## Documentation

### Official Documentation
- **Main Documentation**: https://opensees.github.io/OpenSeesDocumentation
- **Build Instructions**: https://opensees.github.io/OpenSeesDocumentation/developer/build.html
- **Source Repository**: https://github.com/OpenSees/OpenSeesDocumentation

### Community Resources
- **Message Board**: https://opensees.berkeley.edu/community
- **Facebook Group**: https://facebook.com/groups/opensees

### Modeling Questions
For modeling questions, please use:
- OpenSees message board
- OpenSees Facebook group
- GitHub issues (for bugs only)

## Project Cleanup Summary

This repository has been cleaned up to remove redundancy and modernize the build system:

### Removed Components
- Classic Make build system (replaced by CMake)
- Duplicated libraries (Tetgen, eigenAPI, ITPACK, AMGCL)
- Developer directory (experimental code)
- Obsolete CI configuration files

### Benefits
- **Reduced Complexity**: Single build system
- **Smaller Repository**: ~8MB reduction
- **Modern Dependencies**: Conan package management
- **Better Maintainability**: Simplified configuration
- **Preserved Functionality**: All essential features maintained

### Third-Party Libraries (Preserved)
- AMD, ARPACK, BLAS, CBLAS, CSPARSE
- LAPACK, METIS, MUMPS
- SuperLU, Triangle, UMFPACK
- tetgen1.4.3

## Contributing

If you plan to collaborate or use OpenSees as your base code:
1. **FORK** this repository to your own account
2. Work on your fork
3. Submit **PULL REQUESTS** for consideration
4. Direct writes to this repository are not allowed

For forking workflow guidance:
- https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow
- https://www.atlassian.com/git/tutorials/saving-changes

## License

This project contains all revisions to OpenSees source code since Version 2.3.2. Older revisions are available upon request.

---

**Note**: This repository contains all revisions to OpenSees source code since Version 2.3.2. For older revisions, please contact the maintainers.
