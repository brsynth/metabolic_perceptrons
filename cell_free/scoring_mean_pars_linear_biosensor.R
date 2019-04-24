# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: Score and visualise the results from biosensor simulation with the mean sampled parameters

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
folder_for_images <- paste(current_folder, "biosensor_linear_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

parameters_from_fitting = read.csv("biosensor_linear_full_run_1/mean_and_ci.csv")

# collecting_all_data

biosensor_means_data = read.csv(paste(folder_for_data, "biosensor_means_data.csv", sep = '/'), sep = ",")
biosensor_sd_data = read.csv(paste(folder_for_data, "biosensor_sd_data.csv", sep = '/'), sep = ";")

benzoic_acid_range = biosensor_means_data[, c("concentrations")]
# For enzyme inducer adders
rownames(parameters_from_fitting) = c("mean", "sd", "se","95_ci")

mean_pars = unlist(parameters_from_fitting[c("mean"),])

# Visualise biosensor
visu_biosensor= TRUE
if (visu_biosensor) {
  benzoic_acid_range = biosensor_means_data[, c("concentrations")]
  biosensor_df = unlist(calculate_biosensor_linear(parameters = mean_pars, concentrations = benzoic_acid_range))
  means_data_visu = cbind(biosensor_means_data, biosensor_df)
  sd_data_visu = cbind(biosensor_sd_data, 0)
  data_and_model_for_export = cbind(biosensor_means_data, biosensor_sd_data[,c("biosensor")], biosensor_df)
  colnames(data_and_model_for_export) = c("Concentrations", "Data mean", "Data sd", "Model")
  print(data_and_model_for_export)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_visu, sd_data = sd_data_visu)
  
  RMSD_biosensor = all_scores["RMSD"]
  correlation_biosensor = all_scores["correlation"]
  r_squared_biosensor = all_scores["r_squared"]
  r_squared_unweighted_biosensor =all_scores["r_squared_unweighted"]
  
  biosensor_table_results = c(RMSD_biosensor, correlation_biosensor, r_squared_biosensor, r_squared_unweighted_biosensor)
  title= paste("Modeling the benzoic acid biosensor", correlation_biosensor, sep =' ')
  
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "Benzoic acid concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("biosensor_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  write.csv(means_data_visu, paste(folder_for_images, "biosensor_from_mean_of_pars.csv", sep = '/'))
  write.csv(data_and_model_for_export, paste(folder_for_images, "all_biosensor_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
}


table_results_dataframe = cbind(biosensor_table_results)

rownames(table_results_dataframe) = c("RMSD", "Correlation", "R_squared", "Unweighted R squared")
write.csv(table_results_dataframe, paste(folder_for_images, "table_results_correct_order.csv", sep = '/'), row.names = TRUE, quote = FALSE)
