# README
This repository contains the code and instructions to reproduce the results from this [publication](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5208643).

If you use this code or data in your research, please cite:

```bibtex
Yael Frischholz, Océane Hames, Michael Lehning,
Optimizing snow distribution in alpine PV systems: CFD-based design guidelines for power plant layout,
Cold Regions Science and Technology,
2025,
104654,
ISSN 0165-232X,
https://doi.org/10.1016/j.coldregions.2025.104654.
```

To reproduce the results follow the instructions in Section 2.  
## 1. Repository content

### Setup
`requirements.yaml`: conda environment settings

### Parse simulations output
`computeSimuStats.py`: parses the simulations output, computes and writes out statistics.\
`ParserSnowBedFoam.py`: parser for snowBedFoam simulations output. 

### Plot simulations statistics
`hp_plot.py`: plots results.\
`HP_plotter.py`: object containing parameters for plots

## 2. Reproduce results:
The analysis code supports the following platforms : Linux, Darwin, Windows. Here conda-based python environment management is used. 
1. clone this repository:\
`git clone https://github.com/frischwood/snowbed-helio.git`
2. Download simulations folder from EnviDat (https://www.doi.org/10.16904/envidat.589) to your preferred location: `path_to_simu_folder`
3. Setup environment:\
`conda env create -f requirements.yaml`
4. Compute statistics (needed once only):\
`python computeSimuStats.py -sf "path_to_simu_folder"`
5. Plot results:\
`python hp_plot.py -sf "path_to_simu_folder" -v "var"`
var: "mean" or "median" 