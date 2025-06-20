# Cantilever Beam Analysis - Python Version

This directory contains a Python implementation of the 2D Q4 element cantilever beam analysis using OpenSees.

## Files

- `cantilever_q4.py` - Main Python script for the cantilever beam analysis
- `cantilever_q4.tcl` - Original Tcl version for reference
- `RESULTS/` - Directory containing analysis outputs

## Prerequisites

1. **Conda Environment**: You need a conda environment named 'openS' with OpenSees installed
2. **OpenSees Python**: The OpenSees Python module must be available in the environment

## Setup

If you haven't set up the conda environment yet:

```bash
# Create the conda environment
conda create -n openS python=3.8

# Activate the environment
conda activate openS

# Install OpenSees (you may need to install from source or use a specific channel)
# This depends on your OpenSees installation method
```

## Running the Analysis

```bash
# Activate the conda environment
conda activate openS

# Run the Python script
python3 cantilever_q4.py
```

The script will automatically:
- Check if OpenSees is available in the environment
- Create the RESULTS directory if it doesn't exist
- Perform the analysis and save all outputs
- Print model information and results summary

## Model Description

The analysis models a cantilever beam with the following properties:

- **Geometry**: 100" × 10" × 1" (Length × Height × Thickness)
- **Material**: Linear elastic isotropic (E = 30,000 ksi, ν = 0.3)
- **Mesh**: 20 × 4 Q4 elements
- **Boundary Conditions**: Fixed at left end
- **Loading**: 100 kips distributed load at right end

## Output Files

All results are saved in the `RESULTS/` directory:

- `displacements.out` - Node displacements
- `elementForces.out` - Element forces
- `elementStrains.out` - Element strains
- `reactions.out` - Support reactions
- `model_info.txt` - Model information and statistics

## Expected Results

The analysis should produce:
- Tip displacement at the right end
- Complete displacement field
- Element forces and strains
- Support reactions at the fixed end

## Troubleshooting

1. **Import Error**: Make sure the 'openS' conda environment is activated and OpenSees is installed
2. **Analysis Failure**: Check the model parameters and ensure the mesh is valid
3. **Permission Issues**: Ensure you have write permissions in the current directory

## Comparison with Tcl Version

The Python version provides the same analysis as the original Tcl script but with:
- Better error handling
- More detailed output and model information
- Structured code organization
- Additional model statistics 