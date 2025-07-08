"""Visualization functions for vibronic fitting."""

import os
import numpy as np
import matplotlib.pyplot as plt
from src import processing

def plot_current_time_single_step(ca_data, voltage_value, times_to_average, time_at_each_voltage, results_dir="results"):
    """Plot current vs. time for a single voltage step."""
    fig, ax = plt.subplots(figsize=(8, 5))

    time_data = ca_data[:, 0]
    current_data = ca_data[:, 1]

    ax.plot(time_data, current_data, color='royalblue', linewidth=2)
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Current (A)', fontsize=12)
    ax.set_title(f'Current vs. Time (Voltage: {voltage_value:.3f} V)', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)

    # Calculate absolute time window
    avg_start = times_to_average[0]
    avg_end = times_to_average[1]

    # Ensure time window is within data range
    avg_start = max(avg_start, time_data.min())
    avg_end = min(avg_end, time_data.max())
    
    if avg_start < avg_end:
        ax.axvspan(avg_start, avg_end, color='lightcoral', alpha=0.4, 
                  label=f'Averaging Window ({times_to_average[0]}-{times_to_average[1]}s)')
        ax.legend(loc='best')
    
    plt.tight_layout()
    save_path = os.path.join(results_dir, 'plots', 'current_vs_time.png')
    plt.savefig(save_path, dpi=300)
    plt.show()
    
    return fig

def plot_spectrum_with_params(energy, absorbance_data, voltage_array, all_plot_data, 
                             fit_parameters_detailed, vib_fit_range, ep_value, huang_rhys_factor, 
                             spectrum_index=0, results_dir="results"):
    """Create a detailed plot of a spectrum with fit components."""
    idx = spectrum_index 

    # Check if we have data for the requested spectrum index
    if str(idx) not in all_plot_data or idx >= absorbance_data.shape[1]:
        print(f"Warning: Spectrum index {idx} is out of bounds or data not found.")
        error_fig, error_ax = plt.subplots(figsize=(10, 6))
        error_ax.text(0.5, 0.5, f"Invalid spectrum index: {idx} or data missing",
                      ha='center', va='center', fontsize=14, color='red')
        plt.tight_layout()
        error_save_path = os.path.join(results_dir, 'plots', f'error_idx_{idx}.png')
        plt.savefig(error_save_path, dpi=300)
        plt.show()
        return error_fig

    # Create main plot
    fig, ax = plt.subplots(figsize=(8, 6)) 

    # Extract data for easier referencing
    experimental_energy = energy 
    experimental_absorbance = absorbance_data[:, idx]

    plot_data = all_plot_data[str(idx)]
    fit_energy = plot_data[:, 0]      
    total_fit = plot_data[:, 3]
    
    # Extract individual vibronic components
    vibronic_peaks = [
        plot_data[:, 4],  # 0-0 transition
        plot_data[:, 5],  # 0-1 transition
        plot_data[:, 6],  # 0-2 transition
        plot_data[:, 7],  # 0-3 transition
        plot_data[:, 8]   # 0-4 transition
    ]

    # Plot experimental data and total fit
    ax.plot(experimental_energy, experimental_absorbance, color='darkblue', linewidth=3, label='Experimental Data')
    ax.plot(fit_energy, total_fit, color='darkgreen', linewidth=2, label='Total Fit')

    # Plot individual vibronic components
    vib_colors = ['gold', 'orange', 'darkorange', 'red', 'darkred'] 
    vib_labels = ['0-0', '0-1', '0-2', '0-3', '0-4']

    for i, component in enumerate(vibronic_peaks):
        ax.plot(fit_energy, component, color=vib_colors[i], linewidth=1.0) 
        ax.fill_between(fit_energy, 0, component, color=vib_colors[i], alpha=0.5, label=f'Vibronic {vib_labels[i]}')

    # Add amorphous contribution to the plot
    amorphous_energy_limit = 1.9
    amorphous_mask = fit_energy > amorphous_energy_limit 
    if np.any(amorphous_mask):
        high_energy = fit_energy[amorphous_mask]
        model_high_energy = total_fit[amorphous_mask]
        exp_high_energy = experimental_absorbance[amorphous_mask] 
        amorphous_signal = np.maximum(0, exp_high_energy - model_high_energy)
        ax.fill_between(
            high_energy, 
            model_high_energy,
            model_high_energy + amorphous_signal,
            hatch='//', alpha=0.3, color='navy', label='Amorphous'
        )

    # Add title and labels
    voltage_val = voltage_array[idx]
    ax.set_title(f'Voltage: {voltage_val:.2f}V', fontsize=14)
    ax.set_xlabel('Energy (eV)', fontsize=12)
    ax.set_ylabel('Absorbance (a.u.)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Set reasonable axis limits
    ax.set_xlim(np.min(experimental_energy) - 0.05, np.max(experimental_energy) + 0.05) 
    ax.set_ylim(bottom=0, top=np.max(experimental_absorbance) * 1.15)

    # Add parameter text box
    fit_params = fit_parameters_detailed['base_values'] 
    param_table = fit_parameters_detailed['normout_array']   
    cov_diag = fit_parameters_detailed.get('cov_diag') 
    
    param_names = ['Scale', 'E00', 'Width', 'W'] 
    param_text = "Fit Parameters:\n"
    param_text += "---------------------\n"
    for i, name in enumerate(param_names):
        if i*4 < len(param_table):
            value = param_table[i*4]      
            error = param_table[i*4 + 1]  
            param_text += f"{name}: {value:.4f} ± {error:.4f}\n"
    
    param_text += f"S: {huang_rhys_factor:.2f}\n"
    final_ep = fit_params.get('EP_optimized', ep_value) 
    param_text += f"Ep: {final_ep:.4f} eV" 
    if 'EP_optimized' in fit_params:
        param_text += " (optimized)"
    param_text += "\n"

    # Calculate aggregate fraction
    agg_frac, agg_frac_error = processing.calculate_aggregate_fraction(
        energy=fit_energy, 
        plot_data=plot_data, 
        cov_diag=cov_diag 
    )
    param_text += f"Agg. Frac.: {agg_frac:.3f} ± {agg_frac_error:.3f}\n"

    # Add parameter textbox at the top right
    ax.text(0.98, 0.98, param_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.5', fc='aliceblue', alpha=0.9))

    # Add legend
    ax.legend(loc='center right', fontsize=9) 
    plt.tight_layout()
    
    # Save the figure
    save_filename = f'spectrum_v{voltage_val:.2f}.png'
    save_path = os.path.join(results_dir, 'plots', save_filename)
    plt.savefig(save_path, dpi=300)
    plt.show()
    
    return fig
