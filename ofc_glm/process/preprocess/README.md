# Preprocessing functions
Functions to simplify ofc data, as well and extract basic trial contingencies
 
## Contents 

+ batch_concatdata_ofc.m: compress the start, choice, and poke data so that you aren't opening up tons of files
+ batch_ofc_scrubdata.m: script to scrub and simplify the ofc data. removed optOut trials and low firing rate
+ batch_resmooth_OFCdata.m: re-perform spike-train smoothing
+ binspikes.m: put spikes into equal time bins
+ nspikespertrials.m: counts number or spikes on each trial in a session
+ ofc_Ttests.m: a wrapper for selectiveNeuronCalc to tests a neuron's selectivity 
+ parse_choices.m: determines if trial was rewarded, and probability  and volume  of reward
+ selectiveNeuronCalc.m: decides via t test which neurons are selective at the spike count/trial level
+ spiketimes_pertrial.m: parses spike timings into different trials
