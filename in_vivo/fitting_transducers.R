# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: fit and visualise the model results on individual fitting on the 
# hippurate and benzaldehyde transducers and adder

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
folder_for_images <- paste(current_folder, "transducers_single_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

# collecting_all_data

means_data_total = read.csv(paste(folder_for_data, "transducers", "_", "means_data", '.csv', sep = ''))
sd_data_total = read.csv(paste(folder_for_data, "transducers", "_", "sd_data", '.csv', sep = ''))

means_data_adder = read.csv(paste(folder_for_data, "adder", "_", "means_data", '.csv', sep = ''))
sd_data_adder = read.csv(paste(folder_for_data, "adder", "_", "sd_data", '.csv', sep = ''))

means_data_benzal = means_data_total[, c("concentrations", "benzal")]
sd_data_benzal =  sd_data_total[, c("concentrations", "benzal")]

means_data_hip = means_data_total[, c("concentrations", "hip")]
sd_data_hip =  sd_data_total[, c("concentrations", "hip")]

# Remark: biosensor parameters taken from the 1000 fitting dataset

parameters_from_biosensor_fitting = read.csv("biosensor_full_run_1/mean_and_ci.csv")
rownames(parameters_from_biosensor_fitting) = c("mean", "sd", "se","95_ci")

lower_bounds = c("hill_transfer" = parameters_from_biosensor_fitting["mean","hill_transfer"] - parameters_from_biosensor_fitting["95_ci","hill_transfer"], 
                 "Km" = parameters_from_biosensor_fitting["mean","Km"] - parameters_from_biosensor_fitting["95_ci","Km"], 
                 "fold_change" = parameters_from_biosensor_fitting["mean","fold_change"] - parameters_from_biosensor_fitting["95_ci","fold_change"], 
                 "baseline" = parameters_from_biosensor_fitting["mean","baseline"] - parameters_from_biosensor_fitting["95_ci","baseline"], 
                "range_HipO" = 0.7,
                "range_BenZ" = 0.7,
                 "total_enzyme" = 2, "ratio_hip_benz" = 0.5, "cooperativity_resource" = 0.5, "range_resource" = 0.8
                 )


upper_bounds = c("hill_transfer" = parameters_from_biosensor_fitting["mean","hill_transfer"] + parameters_from_biosensor_fitting["95_ci","hill_transfer"], 
                 "Km" = parameters_from_biosensor_fitting["mean","Km"] + parameters_from_biosensor_fitting["95_ci","Km"],
                 "fold_change" = parameters_from_biosensor_fitting["mean","fold_change"] + parameters_from_biosensor_fitting["95_ci","fold_change"],
                 "baseline" = parameters_from_biosensor_fitting["mean","baseline"] + parameters_from_biosensor_fitting["95_ci","baseline"],
                 "range_HipO" = 1.1,
                 "range_BenZ" = 1.1,
                "total_enzyme" = 10, "ratio_hip_benz" = 5, "cooperativity_resource" = 2.2, "range_resource" = 5
                 )


F_to_optimise_transducer <- function(pars_to_optimise, list_of_arguments_of_F) {
  # Obtaining the data necessary for fitting
  means_data_adder <- list_of_arguments_of_F$means_data_adder
  sd_data_adder <- list_of_arguments_of_F$sd_data_adder
  
  means_data_benzal <- list_of_arguments_of_F$means_data_benzal
  sd_data_benzal <- list_of_arguments_of_F$sd_data_benzal
  
  means_data_hip <- list_of_arguments_of_F$means_data_hip
  sd_data_hip <- list_of_arguments_of_F$sd_data_hip
  
  variable_name_adder = list_of_arguments_of_F$variable_name_adder
  variable_name_benzal = list_of_arguments_of_F$variable_name_benzal
  variable_name_hip = list_of_arguments_of_F$variable_name_hip
  
  means_data_benzal_adder = means_data_adder[which(means_data_adder[,c("hip")] == 0), c("benzal", "hip_benzal_adder")]
  sd_data_benzal_adder = sd_data_adder[which(sd_data_adder[,c("hip")] == 0), c("benzal", "hip_benzal_adder")]
  
  means_data_hip_adder = means_data_adder[which(means_data_adder[,c("benzal")] == 0), c("hip", "hip_benzal_adder")]
  sd_data_hip_adder = sd_data_adder[which(sd_data_adder[,c("benzal")] == 0), c("hip", "hip_benzal_adder")]
  
  # Defining the objective function
  objective_function_name = list_of_arguments_of_F$objective_function_name
  if (objective_function_name == "weighted_LSE") {
    objective_function = weighted_LSE
  } else if (objective_function_name == "relative_LSE") {
    objective_function = relative_LSE
  } else {
    objective_function = LSE
  }  
  
  # Calculate simulation results with those parameters
  
  benz_adder_df = data.frame(cbind(concentrations_adder, calculate_hip_benz_adder(parameters = pars_to_optimise, concentrations_benzal = concentrations_adder, concentrations_hip = 0)))
  colnames(benz_adder_df) = c("benzal", variable_name_adder)
  benzal_adder_objective = objective_function(prediction = benz_adder_df, observed = means_data_benzal_adder, variable_name = variable_name_adder, sd_observed = sd_data_benzal_adder)
  
  # Hip transducer from adder:
  hip_adder_df =  data.frame(cbind(concentrations_adder, calculate_hip_benz_adder(parameters = pars_to_optimise, concentrations_benzal = 0, concentrations_hip = concentrations_adder)))
  colnames(hip_adder_df) = c("hip", variable_name_adder)
  hip_adder_objective = objective_function(prediction = hip_adder_df, observed = means_data_hip_adder, variable_name = variable_name_adder, sd_observed = sd_data_hip_adder)
  
  # Benzal transducer objective:
  
  benzal_df =  data.frame(cbind(concentrations_transducers, calculate_benz_transducer(parameters = pars_to_optimise, concentrations_benzal = concentrations_transducers)))
  colnames(benzal_df) = c("concentrations", variable_name_benzal)
  benzal_objective = objective_function(prediction = benzal_df, observed = means_data_benzal, variable_name = variable_name_benzal, sd_observed = sd_data_benzal)
  
  # Hip trasnducer objective:
  
  hip_df = data.frame(cbind(concentrations_transducers, calculate_hip_transducer(parameters = pars_to_optimise, concentrations_hip = concentrations_transducers)))
  colnames(hip_df) = c("concentrations", variable_name_hip)
  hip_objective = objective_function(prediction = hip_df, observed = means_data_hip, variable_name = variable_name_hip, sd_observed = sd_data_hip)
  
  current_objective = hip_objective + benzal_objective + hip_adder_objective + benzal_adder_objective
  # print(paste( hip_objective, benzal_objective, hip_adder_objective, benzal_adder_objective))
  return(current_objective)
}

run_number = 10
set.seed(42)

# Parameters from the biosensor: taken from the biosensor fitting results.
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","hill_transfer"], 
                                    sd = parameters_from_biosensor_fitting["se","hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","Km"], 
                         sd = parameters_from_biosensor_fitting["se","Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","fold_change"], 
                                  sd = parameters_from_biosensor_fitting["se","fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean","baseline"],
                               sd = parameters_from_biosensor_fitting["se","baseline"])

# Parameters for enzymes
range_BenZ_starting_pars = runif(run_number, min = lower_bounds["range_BenZ"], max = upper_bounds["range_BenZ"])
range_HipO_starting_pars = runif(run_number, min = lower_bounds["range_HipO"], max = upper_bounds["range_HipO"])

# Parameters for resource competition: "total_enzyme" = 10, "ratio_hip_benz" = 5, "cooperativity_resource" = 2.2, "range_resource" = 5
total_enzyme_starting_pars = runif(run_number, min = lower_bounds["total_enzyme"], max = upper_bounds["total_enzyme"])
ratio_hip_benz_starting_pars = runif(run_number, min = lower_bounds["ratio_hip_benz"], max = upper_bounds["ratio_hip_benz"])
cooperativity_resource_starting_pars = runif(run_number, min = lower_bounds["cooperativity_resource"], max = upper_bounds["cooperativity_resource"])
range_resource_starting_pars = runif(run_number, min = lower_bounds["range_resource"], max = upper_bounds["range_resource"])

all_pars = NULL
correlations= NULL
all_starting_pars = NULL

concentrations_adder = c(0, 1, 10, 20, 100, 500, 1000)
concentrations_transducers = c(0, 1, 10, 20, 100, 200, 500, 1000)

weighted = FALSE

for (i in 1:run_number) {
  starting_pars = c(# Defined from the biosensor fitting
                    "hill_transfer" = hill_transfer_starting_pars[i],
                    "Km" = Km_starting_pars[i],
                    "fold_change" = fold_change_starting_pars[i],
                    "baseline" = baseline_starting_pars[i],
                    # enzyme characteristics
                    "range_BenZ" = range_BenZ_starting_pars[i],
                    "range_HipO" = range_HipO_starting_pars[i],
                    # resource competition characteristics
                    "total_enzyme" = total_enzyme_starting_pars[i],
                    "ratio_hip_benz" = ratio_hip_benz_starting_pars[i],
                    "cooperativity_resource"= cooperativity_resource_starting_pars[i],
                    "range_resource" = range_resource_starting_pars[i]
  )
  all_starting_pars = rbind(all_starting_pars, data.frame(t(starting_pars)))
  if (weighted) {
    objective = "weighted_LSE"
  } else {
    objective = "LSE"
  }
  list_of_arguments_of_F <- list(
    "means_data_adder" = means_data_adder, 
    "sd_data_adder" = sd_data_adder, 
    "means_data_hip" = means_data_hip, 
    "sd_data_hip" = sd_data_hip, 
    "means_data_benzal" = means_data_benzal, 
    "sd_data_benzal" = sd_data_benzal, 
    "variable_name_adder" = 'hip_benzal_adder', 
    "variable_name_hip" = 'hip', 
    "variable_name_benzal" = 'benzal', 
    "objective_function_name" = objective
  )
  exp_optim <- optim(par = starting_pars, hessian =TRUE, fn = F_to_optimise_transducer, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, list_of_arguments_of_F = list_of_arguments_of_F)
  pars_optimised <- exp_optim$par
  all_pars = rbind(all_pars, pars_optimised)
  
  # Visualise results:
  
  # Hip transducer:
  
  hip_df = calculate_hip_transducer(parameters = pars_optimised, concentrations_hip = concentrations_transducers)
  means_data_hip_visu = cbind(means_data_hip, hip_df)
  sd_data_hip_visu = cbind(sd_data_hip, 0)
  
  colnames(means_data_hip_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_hip_visu) = c("concentrations", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_hip_visu, sd_data = sd_data_hip_visu)
  
  RMSD_hip = all_scores["RMSD"]
  correlation_hip = all_scores["correlation"]
  r_squared_hip = all_scores["r_squared"]
  r_squared_unweighted_hip =all_scores["r_squared_unweighted"]
  
  title= paste("Modeling the hippurate transducer", correlation_hip, sep =' ')   
  subtitle = ""
  df_observed <- melt(data= means_data_hip_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_hip_visu, id.vars = "concentrations")
  
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", 
                                        x_lab = "Hippuric acid concentration (µM)" ,y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = subtitle, 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  folder_to_save = folder_for_images
  saving_graph(filename = paste(i, "hip_transducers_from_and_data.jpeg", sep = '_'), path =folder_to_save)  
  
  # Benzaldehyde transducer
  
  ben_df = calculate_benz_transducer(parameters = pars_optimised, concentrations_ben = concentrations_transducers)
  means_data_ben_visu = cbind(means_data_benzal, ben_df)
  sd_data_ben_visu = cbind(sd_data_benzal, 0)
  
  colnames(means_data_ben_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_ben_visu) = c("concentrations", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_ben_visu, sd_data = sd_data_ben_visu)
  
  RMSD_ben = all_scores["RMSD"]
  correlation_ben = all_scores["correlation"]
  r_squared_ben = all_scores["r_squared"]
  r_squared_unweighted_ben =all_scores["r_squared_unweighted"]
  
  title= paste("Modeling the benzaldehyde transducer", correlation_ben, sep =' ')  
  subtitle = ''
  df_observed <- melt(data= means_data_ben_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_ben_visu, id.vars = "concentrations")
  
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", 
                                        x_lab = "Benzaldehyde concentration (µM)" ,y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = subtitle, 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  folder_to_save = folder_for_images
  saving_graph(filename = paste(i, "ben_transducers_model_and_data.jpeg", sep = '_'), path =folder_to_save)  
  
  # Hippurate benzaldehyde adder
  
  adder_df =  cbind(means_data_adder[,c("hip", "benzal")], calculate_hip_benz_adder(parameters = pars_optimised, concentrations_benzal = concentrations_adder, concentrations_hip = concentrations_adder))
  colnames(adder_df) = c("Hippuric acid", "Benzaldehyde", "hip_benzal_adder")
  
  plot_2D_heatmap(df_means= adder_df, df_sd = NULL, title_base = "Visualising hippuric acid - benzaldehyde adder", 
                  saving = TRUE, x_lab = "(µM)", y_lab = "(µM)",
                  name_for_saving = paste(i, "adder_model_and_data.jpeg", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste(i, "adder_model.jpeg", sep = '_'), path =folder_to_save)
  colnames(adder_df) = c("hip", "benzal", "hip_benzal_adder")
  
  visu_detailed = TRUE
  # Visualise individual rows and columns from the 2D heatmap
  if (visu_detailed) {
    means_data_visu = NULL
    sd_data_visu = NULL
    adder_df_visu = NULL
    for (hip in concentrations_adder) {
      means_data_visu = means_data_adder[which(means_data_adder[,c("hip")] == hip),c("benzal", "hip_benzal_adder")]
      sd_data_visu = sd_data_adder[which(sd_data_adder[,c("hip")] == hip),c("benzal", "hip_benzal_adder")]
      adder_df_visu =  adder_df[which(adder_df[,c("hip")] == hip),c("hip_benzal_adder")]
      
      means_data_visu = cbind(means_data_visu, unname(adder_df_visu))
      sd_data_visu = cbind(sd_data_visu, 0)
      
      colnames(means_data_visu) = c("concentrations", "Data", "Model")
      colnames(sd_data_visu) = c("concentrations", "Data", "Model")
      
      title= paste("From adder with hippuric acid at",hip, "µM", sep =' ')
      subtitle = paste("Hippuric acid is constant at", sep = ' ')
      
      df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
      sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
      
      print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, 
                                            variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                            x_axis = "concentrations",
                                            y_axis = "value", x_lab = "Benzaldehyde concentration (µM)" ,
                                            y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), 
                                            title = title, subtitle = subtitle,
                                            names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
      folder_to_save = folder_for_images
      saving_graph(filename = paste(i, "hip", hip, "model_and_data.jpeg", sep = '_'), path =folder_to_save)
    }
    
    means_data_visu = NULL
    sd_data_visu = NULL
    adder_df_visu = NULL
    
    for (benzal in concentrations_adder) {
      means_data_visu = means_data_adder[which(means_data_adder[,c("benzal")] == benzal),c("hip", "hip_benzal_adder")]
      sd_data_visu = sd_data_adder[which(sd_data_adder[,c("benzal")] == benzal),c("hip", "hip_benzal_adder")]
      adder_df_visu =  adder_df[which(adder_df[,c("benzal")] == benzal),c("hip_benzal_adder")]
      
      means_data_visu = cbind(means_data_visu, unname(adder_df_visu))
      sd_data_visu = cbind(sd_data_visu, 0)
      
      colnames(means_data_visu) = c("concentrations", "Data", "Model")
      colnames(sd_data_visu) = c("concentrations", "Data", "Model")
      
      title= paste("From adder with benzaldehyde at",benzal, "µM", sep =' ')
      subtitle = paste("Benzaldehyde is constant at", benzal, "µM", sep = ' ')
      
      df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
      sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
      
      print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, 
                                            variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                            x_axis = "concentrations",
                                            y_axis = "value", x_lab = "Hippuric acid concentration (µM)" ,
                                            y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), 
                                            title = title, subtitle = subtitle,
                                            names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
      folder_to_save = folder_for_images
      saving_graph(filename = paste(i, "benzal", benzal, "model_and_data.jpeg", sep = '_'), path =folder_to_save)
    }
  }
  
  # Verifi order before melting them
  # adder_df
  # means_data_adder
  
  means_data_adder_score = cbind(means_data_adder, unname(adder_df[,c("hip_benzal_adder")]))
  sd_data_adder_score = cbind(sd_data_adder, 0)
  
  colnames(means_data_adder_score) = c("hip", "benzal", "Data", "Model")
  colnames(sd_data_adder_score) = c("hip", "benzal", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_adder_score, sd_data = sd_data_adder_score)
  
  RMSD_adder = all_scores["RMSD"]
  correlation_adder = all_scores["correlation"]
  r_squared_adder = all_scores["r_squared"]
  r_squared_unweighted_adder =all_scores["r_squared_unweighted"]
  
  # Saving all calculated scores 
  
  correlations = rbind(correlations, c("run" = i, 
                                       "RMSD_adder" = RMSD_adder, 
                                       "RMSD_hip" = RMSD_hip,
                                       "RMSD_benzal" = RMSD_ben,
                                       "r_squared_adder" = r_squared_adder, 
                                       "r_squared_hip" = r_squared_hip,
                                       "r_squared_benzal" = r_squared_ben,
                                       "r_squared_unweighted_adder" = r_squared_unweighted_adder, 
                                       "r_squared_unweighted_hip" = r_squared_unweighted_hip,
                                       "r_squared_unweighted_ben" = r_squared_unweighted_ben,
                                       "adder_correlation" = correlation_adder, 
                                       "hip_correlation" = correlation_hip,
                                       "benzal_correlation" = correlation_ben,
                                       "fucntion" = list_of_arguments_of_F$objective_function_name))
} 
    
all_pars_mean = data.frame(t(apply(all_pars,2,mean)))  
all_pars_sd= data.frame(t(apply(all_pars,2,sd))) 
all_pars_sem = all_pars_sd/sqrt(run_number)
all_pars_95_CI = 1.96 * all_pars_sd/sqrt(run_number)
parameters_mean_and_ci = rbind(all_pars_mean, all_pars_sd, all_pars_sem, all_pars_95_CI)

write.csv(x = correlations, paste(folder_for_images, "correlation.csv", sep = '/'), row.names = FALSE)
write.csv(x = parameters_mean_and_ci, paste(folder_for_images, "mean_and_ci.csv", sep = '/'), row.names = FALSE)

