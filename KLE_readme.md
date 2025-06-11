# 3D Karhunen-Loève Expansion Tool for Reservoir Permeability Field Generation

## Overview

This repository contains a set of Python scripts for generating three-dimensional (3D) permeability field realizations using the Karhunen-Loève (KL) expansion method. These tools are specifically designed to aid stochastic modeling of reservoir properties such as the permeability field. The python scripts are also useful for generating large training datasets  of realistic permeability fields, which can be used to train artificial intelligence (AI)-based surrogate reservoir models (SRMs).

## Theory

### Karhunen-Loève Expansion

The Karhunen-Loève expansion is a method for representing a stochastic process as an infinite linear combination of orthogonal functions. In the context of reservoir modeling, it is used to generate multiple realizations of spatially correlated permeability fields.

For a random field <img src="https://latex.codecogs.com/svg.image?K(x)" alt="K(x)" style="vertical-align: middle;"/> with known mean <img src="https://latex.codecogs.com/svg.image?\mu(x)" alt="mu(x)" style="vertical-align: middle;"/> and covariance function <img src="https://latex.codecogs.com/svg.image?C(x,y)" alt="C(x,y)" style="vertical-align: middle;"/>, the KL expansion is given by:

<!--
$$K(x) = \mu(x) + \sum_{i=1}^{\infty} \sqrt{\lambda_i} \xi_i \phi_i(x)$$
-->

<img src="https://latex.codecogs.com/svg.image?K(x)%20=%20\mu(x)%20+%20\sum_{i=1}^{\infty}%20\sqrt{\lambda_i}%20\xi_i%20\phi_i(x)" alt="KL Expansion Equation" align="center"/>


Where:
- <!--
$\lambda_i$ and $\phi_i(x)$ are eigenvalues and eigenfunctions of the covariance function.
-->
<img src="https://latex.codecogs.com/svg.image?\lambda_i" alt="lambda_i" style="vertical-align: middle;"/> and <img src="https://latex.codecogs.com/svg.image?\phi_i(x)" alt="phi_i(x)" style="vertical-align: middle;"/> are eigenvalues and eigenfunctions of the covariance function.
- <!--
$\xi_i$ are independent standard normal random variables.
-->
<img src="https://latex.codecogs.com/svg.image?\xi_i" alt="xi_i" style="vertical-align: middle;"/> are independent standard normal random variables.
- In practice, the infinite sum is truncated to a finite number of terms based on an energy threshold.

### Log-Normal Permeability Model

Since permeability is strictly positive, we model it as a log-normal field. If <!--
$Z(x)$
--> <img src="https://latex.codecogs.com/svg.image?Z(x)" alt="Z(x)" style="vertical-align: middle;"/> is a normal random field, then <!--
$K(x) = e^{Z(x)}$
-->
<img src="https://latex.codecogs.com/svg.image?K(x)%20=%20e^{Z(x)}" alt="K(x) = e^{Z(x)}" align="center"/> is log-normal. The real distribution parameters are related to the log-normal distribution parameters by:

<!--
$$\sigma_{log} = \sqrt{\ln(1 + (\sigma_{real}/\mu_{real})^2)}$$
$$\mu_{log} = \ln(\mu_{real}) - 0.5\sigma_{log}^2$$
-->
<img src="https://latex.codecogs.com/svg.image?\sigma_{log}%20=%20\sqrt{\ln(1%20+%20(\sigma_{real}/\mu_{real})^2)}" alt="sigma_log equation" align="center"/>
<img src="https://latex.codecogs.com/svg.image?\mu_{log}%20=%20\ln(\mu_{real})%20-%200.5\sigma_{log}^2" alt="mu_log equation" align="center"/>


Where <img src="https://latex.codecogs.com/svg.image?\mu_{real}" alt="mu_real" style="vertical-align: middle;"/> and <img src="https://latex.codecogs.com/svg.image?\sigma_{real}" alt="sigma_real" style="vertical-align: middle;"/> are the desired mean and standard deviation of the physical permeability field.

### Conditional Simulation

The implementation supports conditional simulation. This can be used to ensure that the generated realizations honor known permeability values at specific locations, e.g., at well locations where the property is known with high certainty from core, well log data, pressure transient tests, etc. This is achieved through a kriging adjustment:

<!--
$$K_{cond}(x) = K_{uncond}(x) + \mathbf{C}_{x,obs} \mathbf{C}_{obs,obs}^{-1} (\mathbf{z}_{obs} - \mathbf{z}_{uncond,obs})$$
-->
<p align="center"><img src="https://latex.codecogs.com/svg.image?K_{cond}(x)%20=%20K_{uncond}(x)%20+%20\mathbf{C}_{x,obs}%20\mathbf{C}_{obs,obs}^{-1}%20(\mathbf{z}_{obs}%20-%20\mathbf{z}_{uncond,obs})" alt="Kriging adjustment equation"/></p>


Where:
- <img src="https://latex.codecogs.com/svg.image?K_{cond}(x)" alt="K_cond(x)" style="vertical-align: middle;"/> is the conditioned field at location <img src="https://latex.codecogs.com/svg.image?x" alt="x" style="vertical-align: middle;"/>.
- <img src="https://latex.codecogs.com/svg.image?K_{uncond}(x)" alt="K_uncond(x)" style="vertical-align: middle;"/> is the unconditioned field at location <img src="https://latex.codecogs.com/svg.image?x" alt="x" style="vertical-align: middle;"/>.
- <img src="https://latex.codecogs.com/svg.image?\mathbf{C}_{x,obs}" alt="C_x,obs" style="vertical-align: middle;"/> is the covariance between location <img src="https://latex.codecogs.com/svg.image?x" alt="x" style="vertical-align: middle;"/> and observation locations.
- <img src="https://latex.codecogs.com/svg.image?\mathbf{C}_{obs,obs}" alt="C_obs,obs" style="vertical-align: middle;"/> is the covariance matrix of observation locations.
- <img src="https://latex.codecogs.com/svg.image?\mathbf{z}_{obs}" alt="z_obs" style="vertical-align: middle;"/> are the observed log-permeability values.
- <img src="https://latex.codecogs.com/svg.image?\mathbf{z}_{uncond,obs}" alt="z_uncond,obs" style="vertical-align: middle;"/> are the unconditioned field values at observation locations.

## Implementation

The implementation consists of two main Python modules:

### 1. KL_expansion.py

This module provides the core functionality for generating permeability realizations using the KL expansion. This is a standalone module that can be used to generate realizations for any reservoir model. Key functions:

- `generate_kl_log_normal_real_params_3D()`: Generates multiple 3D log-normal permeability field realizations.
- `plot_realizations_3D()`: Visualizes 2D slices of the 3D realizations.
- `plot_model_3d_grid()`: Creates 3D voxel visualizations of permeability fields.

### 2. kl_realizations_generator.py

This module provides a higher-level wrapper interface for generating, saving, and organizing KL realizations. Key features:

- `KLConfig` class: Configuration settings for the KL expansion.
- Functions for saving realizations in various formats (numpy arrays, .dat files).
- Data splitting functions for organizing realizations into training, validation, and test datasets.
- Directory structure creation with unique identifiers based on configuration hash.

The module is integrated with other modules in the processing pipeline for generating realizations, which are then used to train AI-based SRMs.

### 3. default_configurations.py

This module contains default configuration settings for various aspects of the AI-based SRM including the architecture of the AI-based SRM,optimizer settings and data processing pipeline. It also contains the characterization parameters for the reservoir, fluid, and well configurations, for which the AI-based SRM is developed. 

- Reservoir properties (dimensions, grid size, permeability statistics).
- General settings (working directory, data types, normalization methods).
- Well configurations.
- Data splitting parameters.
- AI-based SRM architecture.
- Optimizer settings.
- Data processing pipeline

## Usage

### Basic Usage

```python
from KL_expansion import generate_kl_log_normal_real_params_3D
import numpy as np

# Generate 10 realizations of a 3D permeability field
permeability_fields, num_modes, grid = generate_kl_log_normal_real_params_3D(
    n_realizations=10,
    Nx=39, Ny=39, Nz=5,
    Lx=2900.0, Ly=2900.0, Lz=80.0,
    real_mean=3.0, real_std=1.5,
    corr_length_fac=0.2,
    energy_threshold=0.95,
    seed=2000
)

# Access the first realization
first_realization = permeability_fields[0]
```

### Conditional Simulation

```python
# Define known permeability values at specific locations
cond_values = {
    (29, 29, 0): 2.0,  # (i,j,k) coordinates and permeability value
    (29, 9, 0): 1.5,
    (9, 9, 0): 1.0,
    (9, 29, 0): 0.5
}

# Generate conditional realizations
permeability_fields, num_modes, grid = generate_kl_log_normal_real_params_3D(
    n_realizations=10,
    Nx=39, Ny=39, Nz=1,
    Lx=2900.0, Ly=2900.0, Lz=80.0,
    real_mean=3.0, real_std=1.5,
    corr_length_fac=0.2,
    energy_threshold=0.95,
    seed=2000,
    cond_values=cond_values
)
```

### Using the Realization Generator

```python
from kl_realizations_generator import KLConfig, generate_and_save_realizations
from default_configurations import DEFAULT_GENERAL_CONFIG, DEFAULT_RESERVOIR_CONFIG, DEFAULT_WELLS_CONFIG

# Create a KL configuration object
kl_config = KLConfig(
    number_of_realizations=100,
    Nx=39, Ny=39, Nz=1,
    Lx=2900.0, Ly=2900.0, Lz=80.0,
    mean=3.0, std=1.5,
    correlation_length_factor=0.2,
    energy_threshold=0.95,
    seed=2000
)

# Generate and save realizations
output_dir = generate_and_save_realizations(
    kl_config,
    general_config=DEFAULT_GENERAL_CONFIG,
    reservoir_config=DEFAULT_RESERVOIR_CONFIG,
    wells_config=DEFAULT_WELLS_CONFIG,
    save_compressed=True,
    overwrite=True
)
```

### Using Command Line Interface

```bash
python kl_realizations_generator.py --num-realizations 100 --nx 39 --ny 39 --nz 1 --lx 2900 --ly 2900 --lz 80 --mean 3.0 --std 1.5 --corr-length 0.2 --energy-threshold 0.95 --seed 2000 --output-dir ./output
```

## Parameter Documentation

### KL Expansion Parameters

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `n_realizations` | Number of permeability field realizations to generate | 10 |
| `Nx`, `Ny`, `Nz` | Grid dimensions in x, y, and z directions | 30, 30, 10 |
| `Lx`, `Ly`, `Lz` | Physical dimensions along x, y, and z in feet/meters | 100.0, 50.0, 20.0 |
| `real_mean` | Mean permeability (physical value) | 3.0 |
| `real_std` | Standard deviation of permeability (physical value) | 1.0 |
| `corr_length_fac` | Correlation length factor (as fraction of max dimension) | 0.2 |
| `energy_threshold` | Fraction of total energy to capture (determines # of KL modes) | 0.95 |
| `seed` | Random seed for reproducibility | 2000 |
| `reverse_order` | Whether to output fields in (Nz,Ny,Nx) order (z,y,x) | False |
| `cond_values` | Dictionary of conditional values at specific locations | None |
| `dtype` | Data type for numerical calculations | np.float32 |

### Configuration and Output Options

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `output_keyword` | Keyword to include in .dat files (e.g., "PERMX") | "PERMX" |
| `comment_prefix` | Prefix for comment lines in .dat files | "--" |
| `add_comments` | Whether to include comments in output files | True |
| `save_compressed` | Whether to save compressed versions of the realizations | False |
| `split_ratio` | Ratio for train/val/test split | (0.7, 0.0, 0.3) |
| `split_sampling_method` | Method for splitting (random or sequential) | "random" |

## Outputs

The KL realization generator produces several outputs:

1. **Realization files**: 
   - Each realization as a separate .dat file in simulator-compatible format
   - Complete numpy arrays of all realizations (.npy and optional .npz compressed)

2. **Grid information**:
   - Grid coordinates as numpy arrays (.npy format)
   - Summary grid information in JSON format

3. **Configuration**:
   - Serialized configuration settings (JSON)
   - Unique directories based on configuration hash for tracking experiments

4. **Data splits**:
   - Training, validation, and test subdirectories (when splitting is enabled)
   - Indices for each split saved as numpy arrays

## Examples

Example of generating and visualizing 3D permeability fields:

```python
from KL_expansion import generate_kl_log_normal_real_params_3D, plot_realizations_3D, plot_model_3d_grid
import numpy as np

# Example usage: Generate 100 realizations of a 3D permeability field with conditional values
n_realizations = 100
permeability_fields_3D, num_modes, grid = generate_kl_log_normal_real_params_3D(
    n_realizations,
    Nx=39, Ny=39, Nz=1,
    Lx=2900.0, Ly=2900.0, Lz=80.0,
    real_mean=3.0, real_std=3,
    corr_length_fac=0.1, energy_threshold=0.95,
    seed=2000,
    reverse_order=False,
    cond_values={
        (29, 29, 0): 2.0, 
        (29, 9, 0): 1.5, 
        (9, 9, 0): 1.0, 
        (9, 29, 0): 0.5
    },
    dtype=np.float32
)

# Plot specified 2D slices from the 3D realizations (slicing along the z dimension)
plot_realizations_3D(
    permeability_fields_3D,
    realization_indices=(0, 100, 5),  # Plot every 5th realization
    z_slices=[0],  # Plot the first z-slice
    Lx=2900.0, Ly=2900.0, Lz=80.0,
    title="3D Permeability Realizations"
)

# Create a 3D grid (voxel) visualization of the first realization
plot_model_3d_grid(
    permeability_fields_3D, grid,
    Lx=2900.0, Ly=2900.0, Lz=80.0,
    sample_index=0, cmap='viridis'
)
```


## References

1. Karhunen, K. (1947). Über lineare Methoden in der Wahrscheinlichkeitsrechnung. Annales Academiae Scientiarum Fennicae, Series A1: Mathematica-Physica.
2. Loève, M. (1978). Probability Theory Vol. II (4th ed.). Springer-Verlag.
3. Zhang, D. (2002). Stochastic Methods for Flow in Porous Media: Coping with Uncertainties. Academic Press.
4. Victor C. Molokwu, Bonaventure C. Molokwu and Mahmoud Jamiolahmady. (2024). Application and effects of physics-based and non-physics-based regularizations in artificial intelligence-based surrogate modelling for highly compressible subsurface flow. Geoenergy Science and Engineering. https://doi.org/10.1016/j.geoen.2023.212474
