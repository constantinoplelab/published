# ofc_glm: generalized linear modeling and clustering of lOFC responses

This is the code used in Hocker D, Brody CD, Savin C, Constantinople CM. Subpopulations of neurons in lOFC encode previous and current rewards at time of choice. eLife 2021;10:e70129

The analysis in this work was based on two primary sets of analysis: 1) building and fitting a generalized linear model to lOFC responses and 2) performing clustering on lOFC neural data to understand the population-level organization of neural responses. 

1. GLM. The general workflow of the code is the following
	- to first load in behavioral and neural data and build a design matrix based on the behavioral data 
	- Choose a basis set for parametrizing the kernels
	- Fit the model and hyperparameters
	- Check model fit, behavior, and encoding features
	
2. Clustering. The general workflow for clustering is the following:
	- Create a feature space for clustering
	- Calculate the PAIRS statistic to check if clusters exist
	- Calculate the gap statistic for k-means clustering on that feature space
	
## Data
The data associated with this work can be found on Zenodo at https://doi.org/10.5281/zenodo.5592702

## Demos
Simple examples have been provided for the following:
* visualizing cluster-averaged responsed
* fitting a single neuron GLM model
* running the clustering pipeline
	- build a feature space
	- calculate the PAIRS statistic
	- calculate the gap statistic
* calculating CPD for a sample neuron
* calculating mutual information for a sample neuron


## Code dependencies
This project was developed with Matlab 2019a, and uses the following toolboxes
+ communication_toolbox
+ optimization_toolbox
+ signal_blocks
+ statistics_toolbox

## Directory structure
- process/ : 
	- preprocess/: download and clean up raw data
	- encoding/: fit a GLM to lOFC data
	- cluster/: parses GLM data to create feature spaces for clustering
- analysis/ : code that works on/analyses/does computation on processed data
	-encoding/: mutual information and coefficient of partial determination claculations
	-cluster/: k-means clustering and evaluation of gap statistic, PAIRS
- vis/: data visualization code
- demos/: Examples of code use
- utils/: random use math and visualization coce


