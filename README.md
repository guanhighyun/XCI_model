# XCI Model 

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

- **Smoldyn** (for spatial simulations)
  - Stochastic spatial simulation software
  - Required for Figure 5 spatial simulations

## Citation

Please cite the associated publication when using this code in your research.

---

*This README provides a comprehensive guide to using the XCI computational model. For specific implementation details, please refer to the comments within individual code files.*
