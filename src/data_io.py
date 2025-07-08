"""Data input/output functions for vibronic fitting."""

import os
import numpy as np
import json
import pandas as pd

# Default results directory
RESULTS_DIR = "results"

def setup_directories(results_dir=RESULTS_DIR):
    """Create necessary directories for storing results."""
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(os.path.join(results_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(results_dir, "data"), exist_ok=True)
    os.makedirs(os.path.join(results_dir, "params"), exist_ok=True)
    
    return {
        'main': results_dir,
        'plots': os.path.join(results_dir, "plots"),
        'data': os.path.join(results_dir, "data"),
        'parameters': os.path.join(results_dir, "params"),
    }

def read_data(data_dir, spec_data_filename="Data_buffer.txt"):
    """Read in raw spectral data and baselines."""
    spec_data = np.loadtxt(os.path.join(data_dir, spec_data_filename))
    baseline_0 = np.loadtxt(os.path.join(data_dir, 'Baseline_0.txt'))
    baseline_100 = np.loadtxt(os.path.join(data_dir, 'Baseline_100.txt'))
    
    return spec_data, baseline_0, baseline_100

def read_ca_data(data_dir, results_dir=RESULTS_DIR, ca_data_filename="CA.csv"):
    """Read the chronoamperometry (CA) data file."""
    ca_file_path = os.path.join(data_dir, ca_data_filename)
    
    if not os.path.exists(ca_file_path):
        raise FileNotFoundError(f"CA data file not found at {ca_file_path}. Cannot proceed without experimental CA data.")
    
    # Try various delimiters that might be used in the CSV file
    try:
        ca_data = np.loadtxt(ca_file_path, delimiter=';', skiprows=1)
    except Exception:
        try:
            ca_data = np.loadtxt(ca_file_path, delimiter=',', skiprows=1)
        except Exception:
            try:
                ca_data = np.loadtxt(ca_file_path, delimiter='\t', skiprows=1)
            except Exception:
                try:
                    ca_data = np.genfromtxt(ca_file_path, delimiter=';', skip_header=1)
                except Exception:
                    raise ValueError("Could not parse CA data file. Please check format.")
    
    if len(ca_data) == 0:
        raise ValueError("CA data file is empty. Cannot proceed without experimental data.")
    
    # Ensure we have at least 3 columns (time, current, voltage)
    if ca_data.shape[1] < 3:
        raise ValueError(f"CA data has insufficient columns: found {ca_data.shape[1]}, need at least 3 (time, current, voltage)")
    
    return ca_data

def save_base_parameters(parameters_dict, filename="params.json", results_dir=RESULTS_DIR):
    """Save base model parameters to a JSON file."""
    filepath = os.path.join(results_dir, 'params', filename)
    with open(filepath, 'w') as f:
        json.dump(parameters_dict, f, indent=4)

def save_fit_summary(param_headers, param_values, spectrum_idx, voltage, results_dir=RESULTS_DIR):
    """Save a summary of fitting results to CSV."""
    fit_summary = pd.DataFrame()
    fit_summary['Parameter'] = param_headers[::4]
    fit_summary['Value'] = param_values[::4]
    fit_summary['Error'] = param_values[1::4]
    fit_summary['95% CI Low'] = param_values[2::4]
    fit_summary['95% CI High'] = param_values[3::4]
    
    filename = os.path.join(results_dir, 'params', f'fit_v{voltage:.2f}.csv')
    
    fit_summary.to_csv(filename, index=False)
    
    return fit_summary
