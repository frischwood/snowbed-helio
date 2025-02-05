# README
This repository contains the code used to produce the results presented in the following publication:\
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
`git clone "my_repo`
2. Download simulations folder from EnviDat () to your preferred location: `path_to_simu_folder`
3. Setup environment:\
`conda env create -f requirements.yaml`
4. Compute statistics (needed once only):\
`python computeSimuStats.py -sf "path_to_simu_folder"`
5. Plot results:\
`python hp_plot.py -sf "path_to_simu_folder" -v "var"`
var: "mean" or "median" 