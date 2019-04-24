# 06/02/19:

The aim of this set of scripts is to reproduce the data and supplementary data from the cell free modeling part of the Metabolic Perceptrons for Neural Computing in Biological Systems article. 

Data is in the data folder.

# Biosensor:
- fitting_biosensor_linear: is not necessary for the pipeline development but allows visualisation of individual parameter fitting runs
- fitting_biosensor_linear_generate_random: is necessary for the rest of the pipeline
It allows fitting of run_number simulations, computes statistics on the 
parameter fitting results and outputs an image of the run_nuber simulation 
results as well as a file containing necessary parameter properties 
for the rest of the project.
- biosensor_linear_visualise_from_pars_sampling generate from pars sampling to see whether it's well constrained or not.
- scoring_mean_pars_linear_biosensor generates scores and images for the biosensor.

# Transducers and adders
- fitting_t_linear: is not necessary for the pipeline development but allows visualisation of individual parameter fitting runs
- fitting_t_linear_generate_random: is necessary for the rest of the pipeline. Remarks: for Km and slope, from 95%, while larger (d) for fold change and baseline.
- transducers_linear_visualise_from_pars_sampling: also for visualisation
- scoring_mean_pars_linear_transducers generates scores and images for the transducers.

# Clusters:
- for 2D clustering experiments presented in Supplementary Figure S7.
Each script corresponds to a subfigure.

# Perceptron prediction:

- perceptron_prediction_complex_classifier: predicts complex classifier (best model)
- perceptron_prediction_or_classifier: predicts or classifier (best model)
- perceptron_cluster_complex_classifier: robustness
- perceptron_cluster_or_classifier: robustness
