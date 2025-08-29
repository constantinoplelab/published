# OFCpaper: analysis of behavior and neural activity early and late in learning of hidden states in a temporal wagering task

This is the code used in [Schiereck, S. S., Pérez-Rivera, D. T., Mah, A., DeMaegd, M. L.,Hocker, D., Ward, R. M., Savin, C., & Constantinople, C. M. Neural dynamics in the orbitofrontal cortex reveal cognitive strategies] (https://www.biorxiv.org/content/10.1101/2024.10.29.620879v2).

The analyses and figures are broken into 5 primary sets:
1. Analysis of expert rat's behavior during the task along with model predictions
2. Analysis of OFC inactivation data
3. Analysis of expert electrophysiology data
4. Analysis of behavior early in learning (block naive) along with model predictions
5. Analysis of block naive electrophysiology data
	
## Data
The data associated with this work can be found on Zenodo at: XXX

## Demos
vis directory contains code to generate every figure in the paper.
process directory contains code to perform every analysis in the paper and must be run before the associated visualization function

## Code dependencies
This project was developed with Matlab 2023a, and uses the following toolboxes
	Curve fitting
	XX
  XX
	

## Directory structure
- analysis/ : code for performing analyses
- process/ : code for creating the processed data files that are used in the visualization functions
- vis/: data visualization code
- utils/: miscellaneous use for processing and visualization
