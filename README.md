![README_1](https://github.com/user-attachments/assets/d2c1e814-6852-4c31-84bf-428afcb143d1)


# Vibronic Fitting Tool
A Python package for analyzing and fitting vibronic transitions in spectroscopic data of polymer semiconductors.

## Physical Model and Theory

This tool implements the Holstein-Spano model for vibronic transitions in conjugated polymers. The model describes optical absorption spectra with a series of vibronic peaks corresponding to electronic transitions coupled with vibrational excitations.

The key parameters in the model include:

- **E00**: The 0-0 transition energy (energy gap between ground state and first excited state)
- **W**: Exciton bandwidth, which relates to intermolecular coupling
- **σ**: Gaussian width parameter for each vibronic transition
- **S**: Huang-Rhys factor, which describes electron-phonon coupling strength
- **Ep**: Phonon energy

The absorption spectrum is modeled as a sum of Gaussian peaks corresponding to different vibronic transitions (0-0, 0-1, 0-2, etc.). The relative intensities of these transitions are governed by the Huang-Rhys factor and the exciton bandwidth.

<div style="display: flex; gap: 10px; align-items: flex-start;">
  <img src="https://github.com/user-attachments/assets/aa4be63a-1722-421f-9df0-7f6ac54a32a9" height="250" />
  <img src="https://github.com/user-attachments/assets/8f79479c-b526-446b-b088-429700fb2c9c" height="250" />
</div>


### Mathematical formulation (Holstein–Spano aggregate)

$$
A(E)=A_0 \sum_{m=0}^{M}
\left[
e^{-S} \ \frac{S^{m}}{m!} \ \ 
\left(1-\frac{W e^{-S}}{2E_p} \ G_m(S)\right)^{2}
\ \ \frac{\exp\left(-\frac{(E-E_m)^2}{2\sigma^{2}}\right)}{\sigma}
\right]
$$

with shared Gaussian width $\sigma$ and vibronic orders $m=0,\dots,M$ (here $M=4$). $A_0$ is a global amplitude (prefactor) capturing overall intensity. Peak centers are

$$
E_m = E_{00} + m \ E_p + \frac{1}{2} \ W \ S^{m}e^{-S},
$$

and the interference term

$$
G_m(S) = \sum_{\substack{n=0, \\\\ n\neq m}}^{N}\frac{S^{n}}{n! \ (n-m)}
$$

### Aggregate fraction from fitted decomposition

Let $D(E)$ be the measured spectrum (experimental data) and $A_{\mathrm{fit}}(E)$ the fitted total. We define the aggregate area from fitted vibronic components $A_m(E)$ and the amorphous excess for $E>E_c$:

$$
A_{\mathrm{agg}}=\sum_{m=0}^{M}\int A_m(E) \ \mathrm{d}E,\qquad
I_{\mathrm{amor}}=\int_{E>E_c} \max \left\\{0, D(E)-A_{\mathrm{fit}}(E) \right\\} \mathrm{d}E.
$$

With radiative-efficiency weights (relative extinction coefficients) $\eta_{\mathrm{agg}}$ and $\eta_{\mathrm{amor}}$,

$$
\phi_{\mathrm{agg}}=\frac{\eta_{\mathrm{agg}} \ A_{\mathrm{agg}}}{\eta_{\mathrm{agg}} \ A_{\mathrm{agg}}+\eta_{\mathrm{amor}} \ I_{\mathrm{amor}}}.
$$

## Installation

1. Clone the repository:
```bash
git clone https://github.com/lukasbongartz/vibronic-fitting-tool.git
cd vibronic-fitting-tool
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Data Requirements

The tool requires the following data files in the `data/raw` directory:

- `Buffer.txt`: Raw spectral data
- `Baseline_0.txt`: 0% transmission baseline
- `Baseline_100.txt`: 100% transmission baseline
- `CA.csv`: Chronoamperometry data (time, current, voltage)


## Data Processing Workflow

The data processing follows these steps:

1. **Data Loading**: Raw spectral data and chronoamperometry data are loaded from the provided files.

2. **Spectral Processing**:
   - Conversion of wavelength (nm) to energy (eV)
   - Normalization using the baseline calibrations
   - Calculation of absorbance values and their statistical variance

3. **Time Correlation**:
   - Extraction of voltage values from chronoamperometry data
   - Synchronization of optical and electrical measurements based on time

4. **Spectral Averaging**:
   - Averaging of spectra within specific time windows for better signal-to-noise ratio
   - The `times_to_average` parameter controls which timepoints are used (typically steady-state regions)

5. **Vibronic Fitting**:
   - The averaged spectrum is fitted with the Holstein-Spano model
   - The fit is performed over a specific energy range (`vib_fit_range`)
   - Multiple vibronic transitions (0-0, 0-1, 0-2, etc.) are summed to model the experimental spectrum

## Data Outputs

The analysis produces several outputs:

1. **Fitted Parameters**: The Scale, E00, W, and offset parameters that best describe the experimental spectrum.

2. **Vibronic Components**: Individual contributions from each vibronic transition (0-0, 0-1, etc.).

3. **Visualizations**: Plots showing the experimental data, fitted model, and individual vibronic components.

4. **CSV Exports**: Processed data files containing both the experimental and fitted data for further analysis.

## Usage

Run the vibronic fitting analysis on the example dataset:

```bash
python main.py
```

To use a custom configuration:

```bash
python main.py --config your_config.json
```

## Configuration

You can customize the analysis by providing a JSON configuration file:

```json
{
  "detector_min_ev": 1.1,
  "detector_max_ev": 3.4,
  "time_at_each_voltage": 30,
  "spectrum_sample_rate": 0.1,
  "times_to_average": [24, 29],
  "ep": 0.17,
  "huang_rhys": 0.95,
  "vib_fit_range": [1.76, 1.97],
  "initial_guess": [0.3, 1.80, 0.06, 0.0]
}
```

## Output

Results are stored in the `results` directory with the following structure:

```
results/
├── data/
│   ├── avg_spectra/      # Averaged absorption spectra
│   └── plot_data/        # Fitted spectra data for each voltage
├── params/               # Fitting parameters
└── plots/                # Visualization outputs
    ├── ca_data/          # Chronoamperometry plots
    └── spectra/          # Spectrum plots with fitted peaks
```

## Key Components

- `main.py`: Main script to run the analysis
- `src/data_io.py`: Data input/output functions
- `src/processing.py`: Spectral data processing
- `src/fitting.py`: Implementation of the vibronic fitting algorithm
- `src/visualization.py`: Plot generation functions
- `src/utils.py`: Utility functions and configuration

## Physical Interpretation

The fitted parameters provide insights into:

1. **Electronic structure**: The E00 transition energy corresponds to the _optical_ band gap
2. **Molecular packing**: The exciton bandwidth W relates to intermolecular coupling strength
3. **Conformational disorder**: The Gaussian width σ reflects energetic disorder

## Citation

If you use this tool in your research, please cite:

```
@article{bongartz2025electron,
	author = {Bongartz, Lukas M. and LeCroy, Garrett and Quill, Tyler J and Siemons, Nicholas and Dijk, Gerwin and Marks, Adam and Cheng, Christina and Kleemann, Hans and Leo, Karl and Salleo, Alberto},
	journal = {Communications Materials},
	title = {Electron--ion coupling breaks energy symmetry in bistable organic electrochemical transistors},
  volume = {6},
  number = {1},
  pages = {241},
	year = {2025}}
```

## License

This project is licensed under the MIT License.
