# Earthquake Model Project

This project simulates a 2D elastic soil model subjected to earthquake loading using OpenSees.

## Running the Simulation

To run the earthquake simulation, execute the Tcl script:

```bash
OpenSees earthquake_model.tcl
```

This will create a `RESULTS/` directory (if it doesn't exist) and populate it with output files from the simulation, including:
- `eq_displacements.out`: Nodal displacements over time.
- `eq_accelerations.out`: Nodal accelerations over time.
- `eq_element_stress.out`: Element stresses over time.
- `eq_nodes.txt`: Information about node coordinates.
- `eq_elements.txt`: Information about element connectivity.
- `eq_model_info.txt`: Summary of model parameters and analysis settings.

## Visualizing the Model

A Python script is provided to visualize the generated mesh, boundary conditions, and an example of the applied load.

To run the visualizer:

```bash
python eq_print_model.py
```

This will generate an image named `RESULTS/eq_model.png` showing the model setup.

## Input Files

- `earthquake_model.tcl`: The main OpenSees script that defines the model, applies loads, and runs the analysis.
- `earthquake_motion.txt`: Contains the time history data for the earthquake motion. (Note: `earthquake_model.tcl` currently uses a trigonometric function for the motion, but this file is included for potential future use or if the script is modified to read from it).
- `eq_print_model.py`: Python script to generate a visual representation of the model.

## Dependencies

- **OpenSees**: Required to run the `earthquake_model.tcl` script. Ensure it's installed and accessible in your PATH.
- **Python 3**: Required to run `eq_print_model.py`.
- **Matplotlib**: A Python library used by `eq_print_model.py` for plotting. Install using pip:
  ```bash
  pip install matplotlib numpy
  ```
