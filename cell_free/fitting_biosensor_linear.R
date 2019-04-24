# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: generate and visualise individual fittings on the model on the biosensor for cell free

rm(list = ls())

# Importing sources

library(ggplot2)
library(deSolve)
library(reshape2)

# Adress of the helper functions file
source("helper_function_synthetic_metabolic_circuits_cell_free.R")

# Define where your working directory is
current_folder = "/home/mkoch/Documents/articles/perceptron_Amir/code/cell_free"

# Define where your data is
folder_for_data <- paste(current_folder, "data/", sep = '/')

# Define where you want to store your images
folder_for_images <- paste(current_folder, "biosensor_linear_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

# collecting_all_data

means_data_total = read.csv(paste(folder_for_data, "biosensor", "_", "means_data", '.csv', sep = ''))
sd_data_total = read.csv(paste(folder_for_data, "biosensor", "_", "sd_data", '.csv', sep = ''), sep = ";")

lower_bounds = c("hill_transfer" = 2, "Km" = 5, "fold_change" = 100, "baseline" = 0.03, "slower_slope" = 8)
upper_bounds = c("hill_transfer" = 2.2, "Km" = 10, "fold_change" = 180, "baseline" = 0.035,  "slower_slope" = 10)

# Function to optimise
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
  biosensor_df = data.frame(cbind(benzoic_acid_range, calculate_biosensor_linear(parameters = pars_to_optimise, concentrations = benzoic_acid_range)))
  colnames(biosensor_df) = c("concentrations", "biosensor")
  current_objective = objective_function(prediction = biosensor_df, observed = means_data, variable_name = 'biosensor', sd_observed = sd_data)
  
  # If problem, visualise the current_objective
  # print(current_objective)
  return(current_objective)
}

# Generating the initial starting parameters
run_number = 10
set.seed(42)

hill_transfer_starting_pars = runif(run_number, min = lower_bounds["hill_transfer"], max = upper_bounds["hill_transfer"])
Km_starting_pars = runif(run_number, min = lower_bounds["Km"], max = upper_bounds["Km"])
fold_change_starting_pars = runif(run_number, min = lower_bounds["fold_change"], max = upper_bounds["fold_change"])
baseline_starting_pars = runif(run_number, min = lower_bounds["baseline"], max = upper_bounds["baseline"])
slower_slope_starting_pars = runif(run_number, min = lower_bounds["slower_slope"], max = upper_bounds["slower_slope"])


all_pars = NULL
correlations = NULL

for (i in 1:run_number) {
  means_data = means_data_total[, c("concentrations", "biosensor")]
  sd_data =  sd_data_total[, c("concentrations", "biosensor")]
  benzoic_acid_range = means_data[, c("concentrations")]
  starting_pars = c("hill_transfer" = hill_transfer_starting_pars[i], 
                    "Km" = Km_starting_pars[i],  
                    "fold_change" = fold_change_starting_pars[i],
                    "baseline" = baseline_starting_pars[i],
                    "slower_slope" = slower_slope_starting_pars[i])
  list_of_arguments_of_F <- list(
    "means_data" = means_data, 
    "sd_data" = sd_data,
    "objective_function_name" = "LSE"
  )
  exp_optim <- optim(par = starting_pars, hessian =TRUE, fn = F_to_optimise_biosensor, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, list_of_arguments_of_F = list_of_arguments_of_F)
  pars_optimised <- exp_optim$par
  all_pars = rbind(all_pars, pars_optimised)
  
  # Visualise the fitting result
  simulated_data =data.frame(cbind(benzoic_acid_range, calculate_biosensor_linear(parameters = pars_optimised, concentrations = benzoic_acid_range)))
  colnames(simulated_data) = c("concentrations", "biosensor")
  
  means_data = cbind(means_data, unname(simulated_data["biosensor"]))
  sd_data = cbind(sd_data, 0)
  
  colnames(means_data) = c("concentrations", "Data", "Model")
  colnames(sd_data) = c("concentrations", "Data", "Model")

  correlation = cor(means_data[,c("Data")], means_data[,c("Model")])
  correlations = rbind(correlations, c("run" = i, "correlation" = correlation, "enzyme" = "fraction", "biosensor" = "fraction"))
  
  title= paste("Modeling the benzoic acid biosensor", correlation, sep =' ')
  
  df_observed <- melt(data= means_data, id.vars = "concentrations")
  sd_observed = melt(data= sd_data, id.vars = "concentrations")
  
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "Benzoic acid concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  folder_to_save = folder_for_images
  saving_graph(filename = paste(i, "biosensor_model_and_data.jpeg", sep = '_'), path =folder_to_save)
}

all_pars_mean = data.frame(t(apply(all_pars,2,mean)))  
all_pars_sd= data.frame(t(apply(all_pars,2,sd))) 
all_pars_sem = all_pars_sd/sqrt(run_number)
all_pars_95_CI = 1.96 * all_pars_sd/sqrt(run_number)
parameters_mean_and_ci = rbind(all_pars_mean, all_pars_sd, all_pars_sem, all_pars_95_CI)

write.csv(x = correlations, paste(folder_for_images, "correlation.csv", sep = '/'), row.names = FALSE)
write.csv(x = parameters_mean_and_ci, paste(folder_for_images, "mean_and_ci.csv", sep = '/'), row.names = FALSE)