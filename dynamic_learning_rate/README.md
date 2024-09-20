# rat_behavior: analysis of rat behavior in a temporal wagering task

This is the code used in [Mah, A., Golden, C. E., Constantinople, C. M. C. Dopamine transients encode reward prediction errors independent of learning rates.](https://www.biorxiv.org/content/10.1101/2024.04.18.590090v2)

The analyses are broken to three primary sets:
1. Analysis of rat's behavior during the task
2. Analysis of behavioral models fit to the rat's behavior
3. Analysis of photometry data
	
## Data
The data associated with this work can be found on Zenodo at: 10.5281/zenodo.13748709

## Demos
vis directory contains code to generate every figure in the paper.
process directory contains code to perform every analysis in the paper.

## Code dependencies
This project was developed with Matlab 2023a, and uses the following toolboxes
  Curve Fitter
	Optimization
	Signal Analyzer

## Directory structure
- process/ : code for analyzing rat beahvior and photometry signals
- model : code for generating model estimates of trial initiation time
- vis/: data visualization code
- utils/: miscellaneous use fitting and visualization code
