# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: fit and visualise the model results on hippuric acid and benzaldehyde fitting for 
# 100 parameters to obtain parameters confidence intervals

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
folder_for_images <- paste(current_folder, "transducers_full_run_1", sep = '/')
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

# Remark: change this for biosensor, should be taken from the pars fited on biosensor

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

run_number = 100
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

# Parameters for resource competition:
total_enzyme_starting_pars = runif(run_number, min = lower_bounds["total_enzyme"], max = upper_bounds["total_enzyme"])
ratio_hip_benz_starting_pars = runif(run_number, min = lower_bounds["ratio_hip_benz"], max = upper_bounds["ratio_hip_benz"])
cooperativity_resource_starting_pars = runif(run_number, min = lower_bounds["cooperativity_resource"], max = upper_bounds["cooperativity_resource"])
range_resource_starting_pars = runif(run_number, min = lower_bounds["range_resource"], max = upper_bounds["range_resource"])

all_pars = NULL
correlations= NULL

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
}

# All simulations have been run and resulting parameters have been saved

# Saving properties of the resulting parameters
all_pars_mean = data.frame(t(apply(all_pars,2,mean)))  
all_pars_sd= data.frame(t(apply(all_pars,2,sd))) 
all_pars_sem = all_pars_sd/sqrt(run_number)
all_pars_95_CI = 1.96 * all_pars_sd/sqrt(run_number)
parameters_mean_and_ci = rbind(all_pars_mean, all_pars_sd, all_pars_sem, all_pars_95_CI)
write.csv(x = parameters_mean_and_ci, paste(folder_for_images, "mean_and_ci.csv", sep = '/'), row.names = FALSE)

# Simulating and visualising the hippurate transducer

random_fit_hippurate_transducer = data.frame(cbind(concentrations_transducers, apply(all_pars,1, calculate_hip_transducer, concentrations_hip = concentrations_transducers)))

unmelted_random_fit = random_fit_hippurate_transducer
colnames(random_fit_hippurate_transducer) = c("concentrations", 1:run_number)
random_fit_hippurate_transducer = melt(data= random_fit_hippurate_transducer, id.vars = "concentrations")
df_observed = melt(means_data_hip, id.vars = "concentrations")
sd_observed = melt(sd_data_hip, id.vars = "concentrations")

plotting_random_pars_and_experiments_local(df_plot_total = random_fit_hippurate_transducer, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                           legend_title = "Constructs", 
                                           x_axis = "concentrations", y_axis = 'value', x_lab = "Hippuric acid concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                           title = "Hippurate transducer", subtitle = "Comparing data and 100 best fits", 
                                           xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")


saving_graph(filename = paste("hippurate_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
hippurate_mean_fit = data.frame(cbind(concentrations_transducers, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(hippurate_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(hippurate_mean_fit, "hippurate_mean_simulation_fit.csv", row.names = FALSE)

# Visualise benzaldehyde transducer

random_fit_benzaldehyde_transducer = data.frame(cbind(concentrations_transducers, apply(all_pars,1, calculate_benz_transducer, concentrations_ben = concentrations_transducers)))

unmelted_random_fit = random_fit_benzaldehyde_transducer
colnames(random_fit_benzaldehyde_transducer) = c("concentrations", 1:run_number)
random_fit_benzaldehyde_transducer = melt(data= random_fit_benzaldehyde_transducer, id.vars = "concentrations")
df_observed = melt(means_data_hip, id.vars = "concentrations")
sd_observed = melt(sd_data_hip, id.vars = "concentrations")

plotting_random_pars_and_experiments_local(df_plot_total = random_fit_benzaldehyde_transducer, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                           legend_title = "Constructs", 
                                           x_axis = "concentrations", y_axis = 'value', x_lab = "Benzaldehyde concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                           title = "Benzaldehyde transducer", subtitle = "Comparing data and 100 best fits", 
                                           xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")


saving_graph(filename = paste("benzaldehyde_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
benzal_mean_fit = data.frame(cbind(concentrations_transducers, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(benzal_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(benzal_mean_fit, "benzal_mean_simulation_fit.csv", row.names = FALSE)

# Calculate and visualise adder

concentration_columns = means_data_adder[,c("hip", "benzal")]
results = apply(all_pars,1, calculate_hip_benz_adder, concentrations_benzal = concentrations_adder, concentrations_hip =concentrations_adder)

random_fit_adder = data.frame(cbind(concentration_columns, results))
adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
colnames(adder_mean_fit) = c("Hippurate concentrations", "Benzaldehyde concentrations","Fluorescence level")
write.csv(adder_mean_fit, "adder_mean_simulations_fit.csv", row.names = FALSE)

colnames(adder_mean_fit) = c("Hippuric acid", "Benzaldehyde","Fluorescence level")

plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising hippuric acid - benzaldehyde adder", saving = FALSE, 
                name_for_saving = NULL, 
                x_lab = "(µM)", y_lab = "(µM)",
                folder_for_concentration_images_strategy = folder_for_images)
saving_graph(filename = paste("adder_heatmap_mean_simulation.jpeg", sep = '_'), path =folder_for_images)


# Visualise individual rows and columns from the 2D heatmap
colnames(adder_mean_fit) = c("hip", "benzal","Fluorescence level")

means_data_visu = NULL
sd_data_visu = NULL
adder_df_visu = NULL
for (hip in concentrations_adder) {
  means_data_visu = means_data_adder[which(means_data_adder[,c("hip")] == hip),c("benzal", "hip_benzal_adder")]
  sd_data_visu = sd_data_adder[which(sd_data_adder[,c("hip")] == hip),c("benzal", "hip_benzal_adder")]
  random_fit_hip = random_fit_adder[which(random_fit_adder[,c("hip")] == hip),]
  random_fit_hip = data.frame(cbind("benzal" = random_fit_hip[,c("benzal")], random_fit_hip[,3:(run_number +2)]))
  
  random_fit_hip = melt(data= random_fit_hip, id.vars = "benzal")
  df_observed = melt(means_data_visu, id.vars = "benzal")
  sd_observed = melt(sd_data_visu, id.vars = "benzal")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_hip, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "benzal", y_axis = 'value', x_lab = "Benzaldehyde concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                             title = paste("Transducer of benzaldehyde with hip =", hip, "µM"), subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("hip", hip, "model_and_data.jpeg", sep = '_'), path =folder_for_images)
}

means_data_visu = NULL
sd_data_visu = NULL
adder_df_visu = NULL

for (benzal in concentrations_adder) {
  means_data_visu = means_data_adder[which(means_data_adder[,c("benzal")] == benzal),c("hip", "hip_benzal_adder")]
  sd_data_visu = sd_data_adder[which(sd_data_adder[,c("benzal")] == benzal),c("hip", "hip_benzal_adder")]
  
  random_fit_benzal = random_fit_adder[which(random_fit_adder[,c("benzal")] == benzal),]
  random_fit_benzal = data.frame(cbind("hip" = random_fit_benzal[,c("hip")], random_fit_benzal[,3:(run_number +2)]))
  
  random_fit_benzal = melt(data= random_fit_benzal, id.vars = "hip")
  df_observed = melt(means_data_visu, id.vars = "hip")
  sd_observed = melt(sd_data_visu, id.vars = "hip")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_benzal, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "hip", y_axis = 'value', x_lab = "Hippuric acid concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                             title = paste("Transducer of hippuric acid with benzaldehyde =", benzal, "µM"), subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("benzal", benzal, "model_and_data.jpeg", sep = '_'), path =folder_for_images)
}

