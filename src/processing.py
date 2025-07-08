"""Data processing functions for vibronic fitting."""

import numpy as np

def process_spectra(spec_data, baseline_100, detector_min_ev=1.1, detector_max_ev=3.4, int_time=10e-3):
    """Process raw spectral data into wavelength and absorbance."""
    wavelength = spec_data[:, 0]
    energy = 1240 / wavelength  # Convert wavelength (nm) to energy (eV)
    
    absorbance = -np.log10(spec_data[:, 1:] / baseline_100[:, 1:2])
    absorbance = np.real(absorbance)
    
    # Calculate variance components for error propagation
    shot_noise = np.sqrt(spec_data[:, 1:] / (40/2.5))
    read_noise = 40
    dark_noise = 200 * int_time
    total_noise = np.sqrt(shot_noise**2 + read_noise**2 + dark_noise**2)
    
    baseline_shot = np.sqrt(baseline_100[:, 1:2])
    baseline_total = np.sqrt(baseline_shot**2 + read_noise**2 + dark_noise**2)
    
    absorbance_variance = np.sqrt((1/spec_data[:, 1:])**2 * (1/np.log(10))**2 * total_noise**2 + 
                                  (1/np.log(10))**2 * (1/baseline_100[:, 1:2])**2 * baseline_total**2)
    
    # Filter by detector range
    energy_mask = (energy >= detector_min_ev) & (energy <= detector_max_ev)
    wavelength_filtered = wavelength[energy_mask]
    energy_filtered = energy[energy_mask]
    absorbance_filtered = absorbance[energy_mask, :]
    variance_filtered = absorbance_variance[energy_mask, :]
    
    return wavelength_filtered, energy_filtered, absorbance_filtered, variance_filtered

def average_spectra(absorbance_data, variance_data, num_voltages, time_at_each_voltage=30, 
                    spectrum_sample_rate=0.1, times_to_average=(25, 29)):
    """Average spectra at each voltage over the specified time range."""
    num_wavelengths = absorbance_data.shape[0]
    avg_absorbance = np.zeros((num_wavelengths, num_voltages))
    avg_variance = np.zeros((num_wavelengths, num_voltages))
    
    spectra_per_voltage = int(time_at_each_voltage / spectrum_sample_rate)
    total_spectra = absorbance_data.shape[1]
    required_spectra = num_voltages * spectra_per_voltage
    
    if required_spectra > total_spectra:
        available_voltage_steps = total_spectra // spectra_per_voltage
        
        if available_voltage_steps < num_voltages:
            num_voltages = available_voltage_steps
            
            if num_voltages == 0:
                raise ValueError("Not enough spectral data to process even a single voltage step.")
    
    for voltage_idx in range(num_voltages):
        start_idx = voltage_idx * spectra_per_voltage
        
        if start_idx + spectra_per_voltage > total_spectra:
            end_idx = total_spectra
            
            if end_idx - start_idx < 5:  # Need minimum 5 spectra for meaningful average
                continue
        else:
            end_idx = (voltage_idx + 1) * spectra_per_voltage
        
        spec_slice = absorbance_data[:, start_idx:end_idx]
        var_slice = variance_data[:, start_idx:end_idx]
        
        # Convert time-based average window to index-based
        start_time_idx = min(int(times_to_average[0] / spectrum_sample_rate), spec_slice.shape[1] - 2)
        end_time_idx = min(int(times_to_average[1] / spectrum_sample_rate), spec_slice.shape[1] - 1)
        
        if start_time_idx >= end_time_idx:
            start_time_idx = 0
            end_time_idx = spec_slice.shape[1] - 1
        
        # Calculate average spectrum
        mean_spectrum = np.mean(spec_slice[:, start_time_idx:end_time_idx+1], axis=1)
        
        # Calculate variance of the average
        variance_sum = np.zeros(num_wavelengths)
        for time_idx in range(start_time_idx, end_time_idx + 1):
            if time_idx < var_slice.shape[1]:
                variance_sum += var_slice[:, time_idx]**2
            
        n_points = end_time_idx - start_time_idx + 1
        mean_variance = np.sqrt(variance_sum) / n_points
        
        avg_absorbance[:, voltage_idx] = mean_spectrum
        avg_variance[:, voltage_idx] = mean_variance
    
    return avg_absorbance, avg_variance

def calculate_aggregate_fraction(energy, plot_data, energy_limit=1.9, eta_agg=1.39, eta_amor=1, cov_diag=None):
    """Calculate the aggregate fraction from the fitted data."""
    # Extract arrays from plot data
    data = plot_data[:, 1]          # Experimental data
    variance = plot_data[:, 2] if plot_data.shape[1] > 2 else None
    fit_total = plot_data[:, 3]     # Total fit
    fit00 = plot_data[:, 4]         # 0-0 transition
    fit01 = plot_data[:, 5]         # 0-1 transition
    fit02 = plot_data[:, 6]         # 0-2 transition
    fit03 = plot_data[:, 7]         # 0-3 transition
    fit04 = plot_data[:, 8]         # 0-4 transition
    
    # Calculate the aggregate component (sum of all vibronic components)
    vibronic_components = [fit00, fit01, fit02, fit03, fit04]
    agg_sum = 0
    for component in vibronic_components:
        agg_sum += np.trapz(y=component[::-1], x=energy[::-1])
    
    # Filter data based on energy threshold
    mask = energy > energy_limit
    energy_filtered = energy[mask]
    data_filtered = data[mask]
    fit_filtered = fit_total[mask]
    
    # Calculate the amorphous contribution
    diff = np.maximum(0, data_filtered - fit_filtered)
    amor_int = np.trapz(y=diff[::-1], x=energy_filtered[::-1])
    
    # Calculate aggregate fraction
    numerator = eta_agg * agg_sum
    denominator = eta_agg * agg_sum + eta_amor * amor_int
    
    if denominator > 0:
        agg_fraction = numerator / denominator
    else:
        agg_fraction = 1.0 if amor_int <= 1e-10 else 0.0
    
    # Calculate uncertainty in aggregate fraction
    if cov_diag is not None and denominator > 0 and np.all(np.isfinite(cov_diag)):
        A = eta_agg * agg_sum
        B = eta_amor * amor_int
        
        # Calculate uncertainty in aggregate area
        agg_jacobian = np.zeros(len(cov_diag))
        
        if len(cov_diag) > 0:
            agg_jacobian[0] = agg_sum
            
        for i in range(1, len(cov_diag)):
            param_importance = np.sqrt(cov_diag[i]) / np.sum(np.sqrt(cov_diag[1:]))
            agg_jacobian[i] = agg_sum * param_importance * 0.5
        
        agg_variance = np.sum(agg_jacobian**2 * cov_diag)
        
        # Estimate uncertainty in amorphous contribution
        if np.any(mask) and variance is not None:
            data_variance_in_region = np.mean(variance[mask]) if np.any(variance[mask] > 0) else 0
            amor_rel_uncertainty = np.sqrt(data_variance_in_region / np.mean(diff)**2) if np.mean(diff) > 0 else 0
            amor_uncertainty = amor_int * amor_rel_uncertainty
        else:
            amor_uncertainty = 0
        
        # Propagate uncertainty through fraction formula
        dAF_dA = B / denominator**2
        dAF_dB = -A / denominator**2
        
        af_variance = (dAF_dA**2 * agg_variance) + (dAF_dB**2 * amor_uncertainty**2)
        
        return agg_fraction, np.sqrt(af_variance)
    else:
        return agg_fraction, 0.07
