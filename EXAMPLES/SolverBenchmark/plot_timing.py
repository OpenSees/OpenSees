#!/usr/bin/env python3
"""
Plot timing results from benchmark CSV file.
Plots num_equations vs time_seconds for different solvers.
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter

def main():
    if len(sys.argv) != 3:
        print("Usage: python plot_timing.py <input.csv> <output.png>")
        sys.exit(1)
    
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    input_csv_arg = sys.argv[1]
    output_png_arg = sys.argv[2]
    
    # Resolve input CSV path
    if os.path.isabs(input_csv_arg):
        input_csv = input_csv_arg
    else:
        # Check in script directory
        input_csv = os.path.join(script_dir, input_csv_arg)
        if not os.path.exists(input_csv):
            # If not found, try the path as-is (might be relative to current working directory)
            input_csv = input_csv_arg
    
    if not os.path.exists(input_csv):
        print(f"Error: Input file not found: {input_csv}")
        sys.exit(1)
    
    # Get directory of input CSV file
    csv_dir = os.path.dirname(os.path.abspath(input_csv))
    
    # Resolve output PNG path
    if os.path.isabs(output_png_arg):
        output_png = output_png_arg
    elif os.path.dirname(output_png_arg):
        # Has a directory component, treat as relative to current working directory
        output_png = output_png_arg
    else:
        # Just a filename, store in same directory as input CSV
        output_png = os.path.join(csv_dir, output_png_arg)
    
    # Read CSV
    df = pd.read_csv(input_csv)
    
    # Filter out skipped runs (status == -999)
    df = df[df['status'] != -999]
    
    # Filter out rows with invalid num_equations or time_seconds
    df = df[(df['num_equations'] > 0) & (df['time_seconds'].notna())]
    
    # Get unique solvers
    solvers = df['solver'].unique()
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    # Plot each solver
    for solver in sorted(solvers):
        solver_data = df[df['solver'] == solver]
        solver_data = solver_data.sort_values('num_equations')
        plt.plot(solver_data['num_equations'], solver_data['time_seconds'], 
                marker='o', label=solver, linewidth=2, markersize=6)
    
    # Set labels and title
    plt.xlabel('Number of Equations', fontsize=12)
    plt.ylabel('Time (seconds)', fontsize=12)
    plt.title('Solver Performance: Number of Equations vs Time', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # Apply engineering notation to x-axis
    ax = plt.gca()
    ax.xaxis.set_major_formatter(EngFormatter(unit=''))
    
    # Use log scale for better visualization if needed
    # Uncomment the following lines if you want log scale:
    # plt.xscale('log')
    # plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {os.path.abspath(output_png)}")

if __name__ == '__main__':
    main()

