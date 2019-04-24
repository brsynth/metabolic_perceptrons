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
folder_for_images <- paste(current_folder, "biosensor_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

# collecting_all_data

means_data = read.csv(paste(folder_for_data, "transducers", "_", "means_data", '.csv', sep = ''))
sd_data = read.csv(paste(folder_for_data, "transducers", "_", "sd_data", '.csv', sep = ''))

means_data = means_data[, c("concentrations", "biosensor")]
sd_data =  sd_data[, c("concentrations", "biosensor")]

run_number = 100

# I copy pasted parameters from that file. Will necessitate one further round of simplification next time.
parameters_from_biosensor_fitting = read.csv("biosensor_full_run_1/mean_and_ci.csv")
rownames(parameters_from_biosensor_fitting) = c("mean", "sd", "se","95_ci")

# mean = mean of fit, sd
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","hill_transfer"], 
                                    sd = parameters_from_biosensor_fitting["se","hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","Km"], 
                         sd = parameters_from_biosensor_fitting["se","Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","fold_change"], 
                                  sd = parameters_from_biosensor_fitting["se","fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","baseline"],
                               sd = parameters_from_biosensor_fitting["se","baseline"])

all_pars = cbind(hill_transfer_starting_pars, Km_starting_pars, fold_change_starting_pars, baseline_starting_pars)
colnames(all_pars) = c("hill_transfer", "Km", "fold_change", "baseline")

# Either at each data poitn or smoothed.
concentrations = means_data[, c("concentrations")]
concentrations = seq(0, 1000)

random_fit_biosensor = data.frame(cbind(concentrations, apply(all_pars,1, calculate_biosensor, concentrations = concentrations)))
unmelted_random_fit = random_fit_biosensor
colnames(random_fit_biosensor) = c("concentrations", 1:run_number)
random_fit_biosensor = melt(data= random_fit_biosensor, id.vars = "concentrations")
df_observed = melt(means_data, id.vars = "concentrations")
sd_observed = melt(sd_data, id.vars = "concentrations")

plotting_random_pars_and_experiments_local(df_plot_total = random_fit_biosensor, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                           legend_title = "Constructs", 
                                           x_axis = "concentrations", y_axis = 'value', x_lab = "Benzoic acid concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                           title = "Benzoic acid biosensor", subtitle = "Comparing data and 100 best fits", 
                                           xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
folder_to_save = folder_for_images
saving_graph(filename = paste("biosensor_from_pars_sampling.jpeg", sep = '_'), path =folder_to_save)

biosensor_mean_fit = data.frame(cbind(concentrations, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(biosensor_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(biosensor_mean_fit, "biosenseor_mean_pars_sampling.csv", row.names = FALSE)
