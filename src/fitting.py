"""Vibronic peak fitting functions."""

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from math import factorial
import warnings
from scipy import stats

def exciton_fit_polaron(energy, params, phonon_energy, start_vib, end_vib, param_scaler, huang_rhys=None):
    """Fit polaron model to exciton data."""
    # Unscale parameters
    scale = params[0] / param_scaler[0]      # Amplitude scaling
    e00 = params[1] / param_scaler[1]        # Energy of 0-0 transition
    width = params[2] / param_scaler[2]      # Peak width
    bandwidth = params[3] / param_scaler[3]  # Exciton bandwidth
    
    if huang_rhys is None:
        raise ValueError("Huang-Rhys factor must be provided, not None")
    
    s = huang_rhys  # Huang-Rhys parameter
    peak_widths = np.ones(end_vib + 1) * width
    
    model_output = np.zeros_like(energy)
    
    for vib_index in range(start_vib, end_vib + 1):
        # Calculate interaction term Gm
        g_term = 0
        for n in range(6):  # Sum over first 6 terms
            if n != vib_index and n - vib_index != 0:
                g_term += (s**n) / (factorial(n) * (n - vib_index))
        
        # Calculate peak position with Holstein model
        peak_position = e00 + vib_index * phonon_energy + 0.5 * bandwidth * s**vib_index * np.exp(-s)
        current_width = peak_widths[vib_index]
        
        # Calculate individual vibronic peak intensity
        vibronic_term = ((np.exp(-s) * s**vib_index) / factorial(vib_index)) * \
                        (1 - (bandwidth * np.exp(-s) / (2 * phonon_energy)) * g_term)**2 * \
                        np.exp(-(energy - peak_position)**2 / (2 * current_width**2)) / current_width
        
        model_output += vibronic_term
    
    return model_output * scale

def fit_vibronic_peaks(energy, fit_range, abs_data, abs_variance, phonon_energy, initial_guess, 
                      param_scaler, spectrum_index, voltage_array, huang_rhys):
    """Fit vibronic peaks to the data."""
    # Parameter names
    param_names = ['Scale', 'E00', 'Width', 'W']
    
    # Prepare initial guess and scaling
    initial_params = initial_guess.copy()
    
    # Extract fitting region
    upper_idx = np.argmax(energy <= fit_range[0])
    lower_idx = np.argmax(energy <= fit_range[1])
    
    fit_energy = energy[lower_idx:upper_idx]
    fit_data = abs_data[lower_idx:upper_idx]
    fit_variance = abs_variance[lower_idx:upper_idx]
    
    # Normalize by maximum for numerical stability
    y_max = np.max(fit_data)
    if y_max <= 1e-9 or np.isnan(y_max):
        y_max = 1.0
    
    fit_data_norm = fit_data / y_max
    fit_variance_norm = fit_variance / (y_max**2) if fit_variance is not None else None
    
    # Adjust initial scale parameter for normalization
    norm_initial_params = initial_params.copy()
    norm_initial_params[0] /= y_max
    scaled_norm_params = norm_initial_params * param_scaler
    
    # Apply light smoothing for better fitting
    smoothed_data = savgol_filter(fit_data_norm, 5, 1)
    
    # Define residual function for least squares
    def residual_func(params):
        return exciton_fit_polaron(fit_energy, params, phonon_energy, 0, 4, param_scaler, huang_rhys) - smoothed_data
    
    # Perform optimization
    result = least_squares(residual_func, scaled_norm_params, method='lm', 
                          ftol=1e-14, xtol=1e-14, gtol=1e-14, max_nfev=5000)
    
    optimized_params = result.x
    residuals = result.fun
    jacobian = result.jac
    
    # Calculate full model and individual components
    full_model = exciton_fit_polaron(energy, optimized_params, phonon_energy, 0, 4, param_scaler, huang_rhys)
    peaks = [
        exciton_fit_polaron(energy, optimized_params, phonon_energy, i, i, param_scaler, huang_rhys)
        for i in range(5)  # 0 to 4 transitions
    ]
    
    # Scale back to original data magnitude
    full_model_scaled = full_model * y_max
    peaks_scaled = [peak * y_max for peak in peaks]
    
    # Calculate degrees of freedom
    dof = upper_idx - lower_idx - len(param_names)
    
    # Calculate covariance matrix and errors
    try:
        cov_matrix = np.linalg.inv(jacobian.T @ jacobian) * np.var(residuals)
        cov_diag = np.diag(cov_matrix)
    except np.linalg.LinAlgError:
        warnings.warn("Covariance matrix calculation failed - matrix nearly singular")
        cov_matrix = np.zeros((len(optimized_params), len(optimized_params)))
        cov_diag = np.zeros(len(optimized_params))
    
    t_crit = stats.t.ppf(0.975, dof)
    
    # Calculate standard errors and confidence intervals
    param_errors = np.zeros(len(optimized_params))
    confidence_intervals = np.zeros((len(optimized_params), 2))
    
    try:
        param_errors = np.sqrt(np.diag(cov_matrix))
        confidence_intervals[:, 0] = optimized_params - t_crit * param_errors
        confidence_intervals[:, 1] = optimized_params + t_crit * param_errors
    except:
        pass
    
    # Calculate chi-square goodness of fit
    chi_square = 0
    model_slice = full_model[lower_idx:upper_idx]
    
    for i in range(len(fit_data_norm)):
        if fit_variance_norm[i] > 0:
            chi_square += (fit_data_norm[i] - model_slice[i])**2 / fit_variance_norm[i]
    
    reduced_chi_square = chi_square / dof
    
    # Unscale parameters for reporting
    unscaled_params = optimized_params / param_scaler
    unscaled_params[0] *= y_max  # Re-scale amplitude by y_max
    
    unscaled_errors = param_errors / param_scaler
    unscaled_errors[0] *= y_max  # Also scale errors
    
    unscaled_ci = np.zeros_like(confidence_intervals)
    for i in range(confidence_intervals.shape[0]):
        unscaled_ci[i, 0] = confidence_intervals[i, 0] / param_scaler[i]
        unscaled_ci[i, 1] = confidence_intervals[i, 1] / param_scaler[i]
    
    unscaled_ci[0, :] *= y_max  # Scale confidence intervals for amplitude
    
    # Report fitted parameters
    voltage = voltage_array[spectrum_index]
    print(f'Voltage: {voltage:.2f} V - Reduced chi square: {reduced_chi_square:.6f}')
    print("\nFitted Parameters:")
    for i, name in enumerate(param_names):
        if i < len(unscaled_params):
            print(f"{name:15s}: {unscaled_params[i]:.6f} ± {unscaled_errors[i]:.6f}")
    print(f"Y-max normalization factor: {y_max:.6f}")
    
    # Prepare output for parameter table
    header = []
    for name in param_names:
        header.append(name)
        header.append(f'{name} error')
        header.append(f'{name} 95% low')
        header.append(f'{name} 95% high')
    
    # Fill output array with parameter values and statistics
    param_table = np.zeros(len(header))
    for i, name in enumerate(param_names):
        if i < len(unscaled_params):
            param_table[i*4] = unscaled_params[i]
            param_table[i*4 + 1] = unscaled_errors[i]
            if i < len(unscaled_ci):
                param_table[i*4 + 2] = unscaled_ci[i, 0]
                param_table[i*4 + 3] = unscaled_ci[i, 1]
    
    # Combine all output data for plotting
    plot_data = np.column_stack((
        energy, abs_data, abs_variance, full_model_scaled,
        peaks_scaled[0], peaks_scaled[1], peaks_scaled[2], peaks_scaled[3], peaks_scaled[4]
    ))
    
    # Create parameter dictionary for output
    param_dict = {name: unscaled_params[i] for i, name in enumerate(param_names) if i < len(unscaled_params)}
    param_dict["red_chi_square"] = reduced_chi_square
    
    # Alternative scale calculation for reference
    global_max = np.max(abs_data)
    param_dict["Scale_alternative"] = unscaled_params[0] * (global_max / y_max) if y_max > 0 else unscaled_params[0]
    
    return header, param_table, unscaled_params, plot_data, cov_diag, param_dict