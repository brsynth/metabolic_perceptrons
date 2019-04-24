# 10/01/19:

The aim of this set of scripts is to reproduce the data and supplementary data from the in vivo modeling part of the Metabolic Perceptrons for Neural Computing in Biological Systems article. 

# General remarks:
All data is in the data folder, adder for the adder and transducers for the biosensor and transducers.
For this paper, run_number = 100.
helper_function contains useful functions for model fitting and visualisation and should not be modified.
All images and csv files have been removed for simplicity. They can be generated from scripts presented below.

# Biosensor fitting
3 scripts concern the initial biosensor fitting.
- fitting biosensor is not necessary for the pipeline development but allows visualisation of individual parameter fitting runs
- fitting_biosensor_generate_random is necessary for the rest of the pipeline. It allows fitting of run_number simulations, computes statistics on the parameter fitting results and outputs an image of the run_nuber simulation results as well as a file containing necessary parameter properties for the rest of the project.
- biosensor_visualise_from_pars_sampling: this file is not necessary for the rest of the pipeline, but allows the user to verify the model parameters are well constrained (cf Supplementary session or Methods)

# Transducers fitting
3 scripts allow the user to fit hippurate and benzaldehyde transducers, using parameters from the biosensor fitting.
- fitting transducers: same idea as fitting biosensor above
- fitting_transducers_generate_random: same as for biosensor, on hippurate and cocaine transducers from adder and single experiments. Visualisation on the full adder as well.
- transducers_visualise_from_pars_sampling: idem

# Cocaine fitting
3 scripts allow the user to toy with cocaine sensor as well.
- fitting_cocaine
- fitting_cocaine_generate_random
- cocaine_visualise_fromçpars_sampling

# Scoring
scoring_mean_pars generates scores from the mean of fitted parameters on all previous simulations. It just necessitates parameters obtained previously (present in the current repository).
