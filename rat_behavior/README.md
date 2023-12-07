# rat_behavior: analysis of rat behavior in a temporal wagering task

This is the code used in [Mah, A., Schiereck, S.S., Bossio, V., Constantinople, C.M. (2023). Distinct value computations support rapid sequential decisions. Nature Communications, 2023](https://www.biorxiv.org/content/10.1101/2023.03.14.532617v2](https://www.nature.com/articles/s41467-023-43250-x)

The analyses are broken to two primary sets:
1. Analysis of rat's behavior during the task
2. Analysis of behavioral models fit to the rat's behavior
	
## Data
The data associated with this work can be found on Zenodo at LINK TO COME

## Demos
Figures directory contains code to generate every figure in the paper

## Code dependencies
This project was developed with Matlab 2019a, and uses the following toolboxes
  Curve Fitter
	Optimization
	Signal Analyzer

## Directory structure
- analysis/ : code for analyzing rat beahvior
- behavioralModel : code for generating model estimates of wait time or trial initiation time
  - TrialInitiationTimeModels - models for generating trial initiation times data
  - WaitTimeModels - models for generating wait time data
  - utils : utility functions for the models
- Figures/: data visualization code
- utils/: random use fitting and visualization code
