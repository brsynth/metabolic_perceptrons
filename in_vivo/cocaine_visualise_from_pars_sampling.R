# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: visualise the model results from sampling on parameters estimation

rm(list = ls())

# Importing sources

library(ggplot2)
library(deSolve)
library(reshape2)

# Adress of the helper functions file
source("helper_function_synthetic_metabolic_circuits_in_vivo.R")

# Define where your working directory is
current_folder = "/home/mkoch/Documents/articles/perceptron_Amir/code/in_vivo"

# Define where your data is
folder_for_data <- paste(current_folder, "data/", sep = '/')

# Define where you want to store your images
folder_for_images <- paste(current_folder, "cocaine_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

# collecting_all_data

means_data_total = read.csv(paste(folder_for_data, "transducers", "_", "means_data", '.csv', sep = ''))
sd_data_total = read.csv(paste(folder_for_data, "transducers", "_", "sd_data", '.csv', sep = ''))

means_data_coc = means_data_total[, c("concentrations", "coc")]
sd_data_coc =  sd_data_total[, c("concentrations", "coc")]

run_number = 100
set.seed(42)

parameters_from_trasnducers_fitting = read.csv("cocaine_full_run_1/mean_and_ci.csv")
rownames(parameters_from_trasnducers_fitting) = c("mean", "sd", "se" ,"95_ci")

# Parameters from the biosensor: taken from the biosensor fitting results.
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "hill_transfer"], 
                                    sd = parameters_from_trasnducers_fitting["se", "hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "Km"],
                         sd = parameters_from_trasnducers_fitting["se", "Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "fold_change"], 
                                  sd = parameters_from_trasnducers_fitting["se", "fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "baseline"], 
                               sd = parameters_from_trasnducers_fitting["se", "baseline"])

# Parameters for enzymes
range_CocE_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "range_CocE"], 
                                 sd = parameters_from_trasnducers_fitting["se", "range_CocE"])

# Parameters for resource competition:
total_enzyme_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "total_enzyme"], 
                                   sd = parameters_from_trasnducers_fitting["se", "total_enzyme"])
ratio_hip_benz_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "ratio_hip_benz"], 
                                     sd = parameters_from_trasnducers_fitting["se", "ratio_hip_benz"])
cooperativity_resource_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "cooperativity_resource"], 
                                             sd = parameters_from_trasnducers_fitting["se", "cooperativity_resource"])
range_resource_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "range_resource"], 
                                     sd = parameters_from_trasnducers_fitting["se", "range_resource"])

all_pars = cbind(hill_transfer_starting_pars, Km_starting_pars, fold_change_starting_pars, baseline_starting_pars,
                 range_CocE_starting_pars, 
                 total_enzyme_starting_pars, ratio_hip_benz_starting_pars, cooperativity_resource_starting_pars, range_resource_starting_pars)
colnames(all_pars) = c("hill_transfer", "Km", "fold_change", "baseline", 'range_CocE', 
                       'total_enzyme', "ratio_hip_benz", "cooperativity_resource", "range_resource")

concentrations_adder = c(0, 1, 10, 20, 100, 500, 1000)
concentrations_transducers = c(0, 1, 10, 20, 100, 200, 500, 1000)

# Simulating and visualising the hippurate transducer

random_fit_cocaine_transducer = data.frame(cbind(concentrations_transducers, apply(all_pars,1, calculate_coc_transducer, concentrations_coc = concentrations_transducers)))

unmelted_random_fit = random_fit_cocaine_transducer
colnames(random_fit_cocaine_transducer) = c("concentrations", 1:run_number)
random_fit_cocaine_transducer = melt(data= random_fit_cocaine_transducer, id.vars = "concentrations")
df_observed = melt(means_data_coc, id.vars = "concentrations")
sd_observed = melt(sd_data_coc, id.vars = "concentrations")

plotting_random_pars_and_experiments_local(df_plot_total = random_fit_cocaine_transducer, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                           legend_title = "Constructs", 
                                           x_axis = "concentrations", y_axis = 'value', x_lab = "Cocaine concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                           title = "Cocaine transducer", subtitle = "Comparing data and 100 best fits", 
                                           xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")

saving_graph(filename = paste("cocaine_from_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
cocaine_mean_fit = data.frame(cbind(concentrations_transducers, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(cocaine_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(cocaine_mean_fit, "cocaine_mean_pars_sampling.csv", row.names = FALSE)
