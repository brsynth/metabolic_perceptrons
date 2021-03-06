# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: fit and visualise the model results on biosensor fitting for 100 parameters 
# to obtain parameters confidenc eintervals

rm(list = ls())

# Importing sources

library(ggplot2)
library(deSolve)
library(reshape2)

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

means_data_total = read.csv(paste(folder_for_data, "transducers", "_", "means_data", '.csv', sep = ''))
sd_data_total = read.csv(paste(folder_for_data, "transducers", "_", "sd_data", '.csv', sep = ''))

means_data = means_data_total[, c("concentrations", "biosensor")]
sd_data =  sd_data_total[, c("concentrations", "biosensor")]
benzoic_acid_range = means_data[, c("concentrations")]

# Defining function to optimise
F_to_optimise_biosensor <- function(pars_to_optimise, list_of_arguments_of_F) {
  # Obtain data 
  means_data <- list_of_arguments_of_F$means_data
  sd_data <- list_of_arguments_of_F$sd_data
  # Define the correct objective function
  objective_function_name = list_of_arguments_of_F$objective_function_name
  if (objective_function_name == "weighted_LSE") {
    objective_function = weighted_LSE
  } else {
    objective_function = LSE
  }
  
  # Calculate results with current par sets and objective (usually Least Squared Errors)
  biosensor_df = data.frame(cbind(benzoic_acid_range, calculate_biosensor(parameters = pars_to_optimise, concentrations = benzoic_acid_range)))
  colnames(biosensor_df) = c("concentrations", "biosensor")
  current_objective = objective_function(prediction = biosensor_df, observed = means_data, variable_name = 'biosensor', sd_observed = sd_data)
  
  # If problem, visualise the current_objective
  # print(current_objective)
  return(current_objective)
}

# Defining bounds for the optimisation

lower_bounds = c("hill_transfer" = 0.5, "Km" = 100, "fold_change" = 1, "baseline" = 100)
upper_bounds = c("hill_transfer" = 2.2, "Km" = 200, "fold_change" = 200, "baseline" = 200)

run_number = 100
# This ensures reproducibility of the parameter fitting process.
set.seed(42)

hill_transfer_starting_pars = runif(run_number, min = lower_bounds["hill_transfer"], max = upper_bounds["hill_transfer"])
Km_starting_pars = runif(run_number, min = lower_bounds["Km"], max = upper_bounds["Km"])
fold_change_starting_pars = runif(run_number, min = lower_bounds["fold_change"], max = upper_bounds["fold_change"])
baseline_starting_pars = runif(run_number, min = lower_bounds["baseline"], max = upper_bounds["baseline"])

# Fit run_number times and keep obtained parameters
all_pars = NULL
for (i in 1:run_number) {
  starting_pars = c("hill_transfer" = hill_transfer_starting_pars[i], 
                    "Km" = Km_starting_pars[i],  
                    "fold_change" = fold_change_starting_pars[i],
                    "baseline" = baseline_starting_pars[i])
  list_of_arguments_of_F <- list(
    "means_data" = means_data, 
    "sd_data" = sd_data, 
    "variable_name" = 'biosensor', 
    "objective_function_name" = "LSE"
  )
  exp_optim <- optim(par = starting_pars, hessian =TRUE, fn = F_to_optimise_biosensor, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, list_of_arguments_of_F = list_of_arguments_of_F)
  pars_optimised <- exp_optim$par
  all_pars = rbind(all_pars, pars_optimised)
}

# Saving properties of the resulting parameters
all_pars_mean = data.frame(t(apply(all_pars,2,mean)))  
all_pars_sd= data.frame(t(apply(all_pars,2,sd))) 
all_pars_sem = all_pars_sd/sqrt(run_number)
all_pars_95_CI = 1.96 * all_pars_sd/sqrt(run_number)
parameters_mean_and_ci = rbind(all_pars_mean, all_pars_sd, all_pars_sem, all_pars_95_CI)
write.csv(x = parameters_mean_and_ci, paste(folder_for_images, "mean_and_ci.csv", sep = '/'), row.names = FALSE)

# Visualisating the results of those parameters fitting
random_fit_biosensor = data.frame(cbind(benzoic_acid_range, apply(all_pars,1, calculate_biosensor, concentrations = benzoic_acid_range)))
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

saving_graph(filename = paste("biosensor_mean_simulation.jpeg", sep = '_'), path =folder_for_images)

biosensor_mean_fit = data.frame(cbind(benzoic_acid_range, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(biosensor_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(biosensor_mean_fit, "biosensor_mean_simulation_fit.csv", row.names = FALSE)
