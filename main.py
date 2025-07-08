#!/usr/bin/env python3
"""Main script for vibronic fitting of spectroelectrochemical data."""

import os
import sys
import argparse
import numpy as np
import json
import matplotlib.pyplot as plt

from src import data_io, processing, fitting, visualization, utils

def parse_arguments():
    parser = argparse.ArgumentParser(description='Vibronic fitting for spectroelectrochemical data')
    parser.add_argument('--config', '-c', type=str, default=None,
                       help='JSON configuration file with experiment parameters (optional)')
    return parser.parse_args()

def main():
    plt.rcParams['figure.max_open_warning'] = 50
    plt.ioff()
    
    args = parse_arguments()
    
    config = utils.VibFittingConfig()
    
    if args.config is not None and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            config_dict = json.load(f)
            config = utils.VibFittingConfig.from_dict(config_dict)
            print(f"Loaded configuration from {args.config}")

    utils.configure_matplotlib()
    data_io.setup_directories(config.results_dir)
    
    # Look directly in data/raw for our data files
    data_dir = os.path.join(os.getcwd(), 'data', 'raw')
    
    if not os.path.exists(data_dir):
        print(f"Error: Data directory {data_dir} not found!")
        sys.exit(1)
    
    required_files = ['Buffer.txt', 'Baseline_0.txt', 'Baseline_100.txt', 'CA.csv']
    missing_files = [f for f in required_files if not os.path.exists(os.path.join(data_dir, f))]

    if missing_files:
        print(f"Error: The following required files are missing: {', '.join(missing_files)}")
        print(f"Please ensure these files are present in {data_dir}")
        sys.exit(1)
    
    print(f"Using data from: {data_dir}")
    
    spec_data_filepath = os.path.join(data_dir, "Buffer.txt")
    ca_data_filepath = os.path.join(data_dir, "CA.csv")
    
    try:
        spec_data, baseline_0, baseline_100 = data_io.read_data(
            data_dir, 
            spec_data_filename=os.path.basename(spec_data_filepath)
        )
        ca_data = data_io.read_ca_data(
            data_dir, 
            ca_data_filename=os.path.basename(ca_data_filepath)
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during data loading: {e}")
        sys.exit(1)

    try:
        wavelength, energy, absorbance, absorbance_variance = processing.process_spectra(
            spec_data, baseline_100,
            config.detector_min_ev, 
            config.detector_max_ev, 
            config.int_time
        )
    except Exception as e:
        print(f"Error processing spectral data: {e}")
        sys.exit(1)
    
    try:
        if ca_data.size == 0:
            print("Error: CA data is empty after loading.")
            sys.exit(1)
            
        voltage_value = np.mean(ca_data[:, 2])
        voltage_array = np.array([voltage_value])
        
        print(f"Processing voltage step at {voltage_value:.3f} V")
    except Exception as e:
        print(f"Error processing voltage information: {e}")
        sys.exit(1)
    
    visualization.plot_current_time_single_step(
        ca_data,
        voltage_array[0],
        config.times_to_average,
        config.time_at_each_voltage,
        config.results_dir
    )

    avg_absorbance, avg_variance = processing.average_spectra(
        absorbance,
        absorbance_variance,
        num_voltages=1,
        time_at_each_voltage=config.time_at_each_voltage,
        spectrum_sample_rate=config.spectrum_sample_rate,
        times_to_average=config.times_to_average
    )
    
    if avg_absorbance.ndim > 1 and avg_absorbance.shape[1] == 1:
        spectrum_data = avg_absorbance[:, 0]
        spectrum_variance = avg_variance[:, 0]
    else:
        spectrum_data = avg_absorbance
        spectrum_variance = avg_variance
    
    # Set initial guess for Scale parameter
    config.initial_guess[0] = 0.48
    
    spectrum_index = 0
    current_voltage = voltage_array[0]
    
    print(f'Performing vibronic fit at {current_voltage:.2f} V')
    
    try:
        param_headers, param_table, fitted_params, plot_data, cov_diag, params_dict = fitting.fit_vibronic_peaks(
            energy,
            config.vib_fit_range,
            spectrum_data,
            spectrum_variance,
            config.ep,
            config.initial_guess,
            (1.0 / np.abs(config.upper_bounds)),
            spectrum_index,
            voltage_array,
            huang_rhys=config.huang_rhys
        )
        
    except Exception as e:
        print(f'Error during vibronic fitting: {e}')
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    data_io.save_fit_summary(param_headers, param_table, 0, current_voltage, config.results_dir)
    
    all_params = {"0": params_dict}
    data_io.save_base_parameters(all_params, results_dir=config.results_dir)

    print('Generating spectrum plot with fit parameters...')
    
    fit_params_for_plot = {
        "base_values": params_dict,
        "normout_array": param_table,
        "cov_diag": cov_diag
    }

    visualization.plot_spectrum_with_params(
        energy,
        spectrum_data.reshape(-1, 1),
        voltage_array,
        all_plot_data={"0": plot_data},
        fit_parameters_detailed=fit_params_for_plot,
        vib_fit_range=config.vib_fit_range, 
        ep_value=config.ep,
        huang_rhys_factor=config.huang_rhys,
        spectrum_index=0,
        results_dir=config.results_dir
    )

    # Save plot data for CSV export
    data_to_save = np.column_stack((
        plot_data[:, 0],  # Energy
        spectrum_data,    # Experimental data
        plot_data[:, 3],  # Total model
        plot_data[:, 4],  # 0-0 transition
        plot_data[:, 5],  # 0-1 transition
        plot_data[:, 6],  # 0-2 transition
        plot_data[:, 7],  # 0-3 transition
        plot_data[:, 8]   # 0-4 transition
    ))
    
    plot_data_header = "Energy_eV,Experimental_Abs,Total_Model,Vib_0-0,Vib_0-1,Vib_0-2,Vib_0-3,Vib_0-4"
    plot_data_dir = os.path.join(config.results_dir, "data", "plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    plot_data_file = os.path.join(plot_data_dir, f'fit_v{current_voltage:.2f}.csv')
    
    np.savetxt(plot_data_file, data_to_save, delimiter=',', header=plot_data_header, comments='')
    
    print("\nExporting averaged spectrum to CSV...")
    avg_export_dir = os.path.join(config.results_dir, "data", "avg_spectra")
    os.makedirs(avg_export_dir, exist_ok=True)
    
    avg_spec_filename = f"avg_spectrum_v{current_voltage:.2f}.csv"
    avg_spec_filepath = os.path.join(avg_export_dir, avg_spec_filename)
    
    avg_spec_export = np.column_stack((energy, spectrum_data))
    avg_spec_header = "Energy_eV,Averaged_Absorbance"
    
    np.savetxt(avg_spec_filepath, avg_spec_export, delimiter=',', header=avg_spec_header, comments='')
    
    print("\nFitting complete.")

if __name__ == '__main__':
    main()
