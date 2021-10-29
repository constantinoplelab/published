# Demos
 These are demos that show the core functionality and analyses for the glm and clustering. In order to run these, you will first need to download a copy of the raw data from zenodo at (https://doi.org/10.5281/zenodo.5592702)[https://doi.org/10.5281/zenodo.5592702].

You should only have to set three main path variables in these script in order to run on your machine:
+ `datadir`: this is where the raw data from zenodo lives on your machine
+ `codedir`: the directory on your local machine where this repo lives
+ `savedir`: where you want results of glm fits and anlyses to be stored

The mutual information and CPD analyses are model based, and will first require you to generate glm fits for data to peform them. Other than that, the other demos are standalone and can be run in any order

 ## Contents 

+ demo_clusteringPipeline.m: demo for running the clustering pipeline and the gap statistic
+ demo_condPSTH.m: plot conditional PSTHs of neurons
+ demo_cpd.m: calcualtes CPD for an example cell
+ demo_figure1.m: recreates clustering results of Figure 1, on PSTH-based feature space
+ demo_glm.m: fits a sample GLM 
+ demo_MI.m: mutual information between spikes and trial-level behavior/outcome stimuli
