"""Utility functions for vibronic fitting."""

import numpy as np
import matplotlib.pyplot as plt
import json

# Physical constants
PLANCK_CONST = 6.62607015e-34 
SPEED_OF_LIGHT = 2.99792458e8  
ELECTRON_CHARGE = 1.60217663e-19 
WL_CONVERSION = PLANCK_CONST*SPEED_OF_LIGHT/ELECTRON_CHARGE*1e9 

def configure_matplotlib(font_size=14, line_width=1.5):
    """Configure matplotlib settings for consistent plots."""
    plt.rcParams['font.size'] = font_size
    plt.rcParams['axes.linewidth'] = line_width
    plt.rcParams['axes.grid'] = True
    plt.rcParams['lines.linewidth'] = line_width
    plt.rcParams['lines.markersize'] = 4
    plt.rcParams['figure.figsize'] = [6, 4.5]
    plt.rcParams['figure.facecolor'] = 'white'

class VibFittingConfig:
    """Configuration parameters for vibronic fitting."""
    
    def __init__(self):
        # Default parameters for molar absorptivity ratio
        self.eta_amor = 1
        self.eta_agg = 1.39
        
        # Detector limits parameters
        self.detector_min_ev = 1.1
        self.detector_max_ev = 3.4
        
        # Directory and file parameters
        self.results_dir = "results"
        self.time_step = 0.1
        self.int_time = 10e-3
        self.time_at_each_voltage = 30
        self.spectrum_sample_rate = 0.1
        self.times_to_average = [24, 29]
        
        # Vibronic fitting parameters
        self.ep = 0.17
        self.e00 = 1.77
        self.w = 0.04
        self.sigma = 0.07
        self.huang_rhys = 0.95
        
        # Initial guesses and bounds
        self.initial_guess = [0.3, 1.80, 0.06, 0.0]
        self.lower_bounds = [0.001, 1.5, 0.005, 0]
        self.upper_bounds = [2.0, 2.2, 0.15, 0.05]
        
        # Fitting range
        self.vib_fit_range = [1.76, 1.97]

    def to_dict(self):
        """Convert configuration to dictionary."""
        return {
            'eta_amor': self.eta_amor,
            'eta_agg': self.eta_agg,
            'detector_min_ev': self.detector_min_ev,
            'detector_max_ev': self.detector_max_ev,
            'results_dir': self.results_dir,
            'time_step': self.time_step,
            'int_time': self.int_time,
            'time_at_each_voltage': self.time_at_each_voltage,
            'spectrum_sample_rate': self.spectrum_sample_rate,
            'times_to_average': self.times_to_average,
            'ep': self.ep,
            'e00': self.e00,
            'w': self.w,
            'sigma': self.sigma,
            'huang_rhys': self.huang_rhys,
            'initial_guess': self.initial_guess,
            'lower_bounds': self.lower_bounds,
            'upper_bounds': self.upper_bounds,
            'vib_fit_range': self.vib_fit_range
        }
    
    @classmethod
    def from_dict(cls, config_dict):
        """Create configuration from dictionary."""
        config = cls()
        for key, value in config_dict.items():
            if hasattr(config, key):
                setattr(config, key, value)
        return config

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super(NumpyEncoder, self).default(obj)