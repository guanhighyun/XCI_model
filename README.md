# X-Chromosome Inactivation (XCI) Computational Model

## Overview

This repository contains computational models and analysis code for studying X-chromosome inactivation (XCI) mechanisms. The project implements mathematical models to understand the molecular interactions between Xist, SPEN, and activator proteins during X-chromosome inactivation, using both deterministic ordinary differential equations and stochastic spatial simulations.

## Project Structure

### Main Directories

- **`Fig1/`** - Code for generating Figure 1 (basic model dynamics)
- **`Fig2/`** - Code for generating Figure 2 (parameter analysis and comparison)
- **`Fig3/`** - Code for generating Figure 3 (model variants)
- **`Fig4/`** - Code for generating Figure 4 (evolutionary algorithm results)
- **`Fig5/`** - Code for generating Figure 5 (spatial simulations)
- **`Supplemental figures/`** - Code for all supplemental figures
- **`EA_for_activator_inhibition_model/`** - Evolutionary algorithm for parameter optimization (activator inhibition model)
- **`EA_for_SPEN_depletion_model/`** - Evolutionary algorithm for parameter optimization (SPEN depletion model)
- **`Smoldyn_simulations_for_Fig5/`** - Spatial stochastic simulations using Smoldyn

### Key Files

- **`Parameters_for_activator_inhibition.csv`** - Optimized parameters for the activator inhibition model
- **`Parameters_for_SPEN_depletion.csv`** - Optimized parameters for the SPEN depletion model

## Models

The repository implements two main mechanistic models for X-chromosome inactivation:

### 1. Activator Inhibition Model
This model considers how SPEN binding to Xist affects activator-mediated gene expression regulation.

### 2. SPEN Depletion Model  
This model examines the effects of SPEN depletion on X-chromosome inactivation dynamics.

Both models are described by systems of ordinary differential equations that capture the molecular interactions between:
- **Activator proteins** - Regulatory factors affecting gene expression
- **Xist RNA** - Long non-coding RNA essential for X-inactivation
- **SPEN protein** - Silencing factor that binds to Xist

## Requirements

### Software Dependencies

- **Python 3.x** with packages:
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `pandas`
  - `seaborn`
  - `deap` (Distributed Evolutionary Algorithms)
  - `jupyter`

- **MATLAB** (for .m files)
  - Required for ODE integration and figure generation
  - Uses `ode15s` solver for stiff differential equations

- **Smoldyn** (for spatial simulations)
  - Stochastic spatial simulation software
  - Required for Figure 5 spatial simulations

### Installation

```bash
# Install Python dependencies
pip install numpy scipy matplotlib pandas seaborn deap jupyter

# MATLAB installation required separately
# Smoldyn installation required for spatial simulations
```

## Usage

### Running Figure Generation Code

Each figure directory contains self-contained code:

```bash
# For MATLAB figures
cd Fig1/
matlab -batch "run('Fig1C.m')"

# For Python/Jupyter analysis
cd Fig2/
jupyter notebook Fig2D.ipynb
```

### Parameter Optimization

To run evolutionary algorithm optimization:

```bash
# For activator inhibition model
cd EA_for_activator_inhibition_model/
python EA.py

# For SPEN depletion model  
cd EA_for_SPEN_depletion_model/
python EA.py
```

### Spatial Simulations

To run Smoldyn spatial simulations:

```bash
cd Smoldyn_simulations_for_Fig5/
python Execute_smoldyn.py
```

## Model Parameters

### Key Parameters (see parameter CSV files for full sets)

- **`a_act`** - Activator synthesis rate (molecules/min)
- **`d_act`** - Activator degradation rate (1/min)
- **`a_x`** - Xist synthesis rate (molecules/min)  
- **`d_x`** - Xist degradation rate (1/min)
- **`K_n`** - Half-saturation constant for SPEN inhibition
- **`n`** - Hill coefficient for SPEN inhibition
- **`m`** - Hill coefficient for Xist-SPEN binding
- **`K_S`** - SPEN binding affinity
- **`sT`** - Total SPEN concentration
- **`XbsT`** - Total bound Xist-SPEN complex

### Parameter Units

All parameters are in molecular quantities assuming a well-stirred nuclear compartment:
- Concentrations: molecules
- Rate constants: molecules/min or 1/min
- Binding constants: molecules

## Evolutionary Algorithm Details

The parameter optimization uses DEAP (Distributed Evolutionary Algorithms in Python) with:

- **Population size**: Typically 100-300 individuals
- **Generations**: 50-100 generations  
- **Selection**: Tournament selection
- **Crossover**: Uniform crossover
- **Mutation**: Gaussian mutation
- **Fitness function**: Least squares difference to experimental data

### Optimization Ranges

Parameters are optimized within biologically reasonable ranges (see `EA.py` files for specific bounds).

## File Formats

- **`.m`** - MATLAB scripts for ODE integration and plotting
- **`.py`** - Python scripts for evolutionary algorithms and analysis
- **`.ipynb`** - Jupyter notebooks for interactive analysis
- **`.csv`** - Parameter sets and optimization results
- **`.mat`** - MATLAB data files containing simulation results
- **`.cfg`** - Smoldyn configuration files for spatial simulations

## Results and Analysis

### Statistical Analysis

The project includes statistical comparison between models using:
- Kolmogorov-Smirnov tests for parameter distribution comparison
- Parameter sensitivity analysis
- Model selection criteria

### Visualization

All figures are generated with publication-quality formatting:
- High-resolution outputs suitable for manuscripts
- Consistent styling and color schemes  
- Statistical significance indicators

## Troubleshooting

### Common Issues

1. **MATLAB Path Issues**: Ensure all functions are in MATLAB path
2. **Python Dependencies**: Install all required packages listed above
3. **Smoldyn Installation**: Follow Smoldyn documentation for platform-specific installation
4. **Memory Issues**: Large parameter sweeps may require significant RAM

### Performance Tips

- Use parallel processing where available (MATLAB Parallel Computing Toolbox)
- Adjust ODE solver tolerances for speed vs accuracy tradeoffs
- Consider reducing parameter space for faster optimization

## Contributing

When contributing to this project:
1. Maintain consistent file naming conventions
2. Document all new parameters and their biological meaning
3. Include appropriate error handling in analysis scripts
4. Test code with provided example parameter sets

## Citation

Please cite the associated publication when using this code in your research.

## License

[Include appropriate license information]

## Contact

[Include contact information for questions and support]

---

*This README provides a comprehensive guide to using the XCI computational model. For specific implementation details, please refer to the comments within individual code files.*
