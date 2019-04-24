# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: fit and visualise the model results on individual fitting on all experiments except biosensor, for 100 runs

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
folder_for_images <- paste(current_folder, "t_linear_full_range_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

# Collecting all useful data

biosensor_means_data = read.csv(paste(folder_for_data, "biosensor_means_data.csv", sep = '/'), sep = ",")
biosensor_sd_data = read.csv(paste(folder_for_data, "biosensor_sd_data.csv", sep = '/'), sep = ";")

coce_means_data = read.csv(paste(folder_for_data, "coce_coc_means_data.csv", sep = '/'), sep = ",")
coce_sd_data = read.csv(paste(folder_for_data, "coce_coc_sd_data.csv", sep = '/'), sep = ",")

hipo_means_data = read.csv(paste(folder_for_data, "hipo_hip_means_data.csv", sep = '/'), sep = ",")
hipo_sd_data = read.csv(paste(folder_for_data, "hipo_hip_sd_data.csv", sep = '/'), sep = ",")

fixed_inducer_adder_means_data = read.csv(paste(folder_for_data, "fixed_inducer_adder_means_data.csv", sep = '/'), sep = ";")
fixed_inducer_adder_sd_data = read.csv(paste(folder_for_data, "fixed_inducer_adder_sd_data.csv", sep = '/'), sep = ";")

fixed_enzyme_adder_means_data = read.csv(paste(folder_for_data, "adder_1_3_means_data.csv", sep = '/'), sep = ",")
fixed_enzyme_adder_sd_data = read.csv(paste(folder_for_data, "adder_1_3_sd_data.csv", sep = '/'), sep = ",")

benzamid_means_data = read.csv(paste(folder_for_data, "benzamid_means_data.csv", sep = '/'), sep = ",")
benzamid_sd_data = read.csv(paste(folder_for_data, "benzamid_sd_data.csv", sep = '/'), sep = ",")

biphenyl_means_data = read.csv(paste(folder_for_data, "biphenyl_diol_means_data.csv", sep = '/'), sep = ",")
biphenyl_sd_data = read.csv(paste(folder_for_data, "biphenyl_diol_sd_data.csv", sep = '/'), sep = ",")

# For visualisation on a smaller subset of data
benzamid_means_data_small = benzamid_means_data[which(benzamid_means_data[,c("benzamid")] %in% c(0, 10, 100)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]
benzamid_sd_data_small = benzamid_sd_data[which(benzamid_sd_data[,c("benzamid")] %in% c(0, 10, 100)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]

benzamid_means_data_small = benzamid_means_data_small[which(benzamid_means_data_small[,c("benzamid_enz")] %in% c(0.1, 0.3, 1, 3)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]
benzamid_sd_data_small = benzamid_sd_data_small[which(benzamid_sd_data_small[,c("benzamid_enz")] %in%  c(0.1, 0.3, 1, 3)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]

biphenyl_means_data_small = biphenyl_means_data[which(biphenyl_means_data[,c("biphenyl")] %in% c(0, 10, 100)), c("biphenyl", "biphenyl_enz", "biphenyl_enz_biphenyl")]
biphenyl_sd_data_small = biphenyl_sd_data[which(biphenyl_sd_data[,c("biphenyl")] %in% c(0, 10, 100)), c("biphenyl", "biphenyl_enz", "biphenyl_enz_biphenyl")]

# Parameters from biosensor fitting will be taken from the 100 random transducers

parameters_from_biosensor_fitting = read.csv("biosensor_linear_full_run_1/mean_and_ci.csv")
rownames(parameters_from_biosensor_fitting) = c("mean", "sd", "se" ,"95_ci")


lower_bounds = c("hill_transfer" = 2.1, 
                 "Km" = parameters_from_biosensor_fitting["mean","Km"] - parameters_from_biosensor_fitting["95_ci","Km"], 
                 "fold_change" = parameters_from_biosensor_fitting["mean","fold_change"] - parameters_from_biosensor_fitting["sd","fold_change"], 
                 "baseline" = parameters_from_biosensor_fitting["mean","baseline"] - parameters_from_biosensor_fitting["sd","baseline"], 
                 "slower_slope"  = parameters_from_biosensor_fitting["mean","slower_slope"] - parameters_from_biosensor_fitting["95_ci","slower_slope"],
                 
                 "range_HipO" = 200, "HipO_constant" = 0.1, "hippurate_constant" = 20, "hill_HipO" = 0.5, "hill_hippurate" = 0.2,
                 "range_CocE" = 100, "CocE_constant" = 0.1, "cocaine_constant" = 10, "hill_CocE" = 0.5, "hill_cocaine" = 0.2,
                 
                 "range_benzamid_enz" = 50, "benzamid_enz_constant" = 0.1, "benzamid_constant" = 1, "hill_benzamid_enz" = 0.2, "hill_benzamid" = 0.2,
                 "range_biphenyl_enz" = 20, "biphenyl_enz_constant" = 0.1, "biphenyl_constant" = 10, "hill_biphenyl_enz" = 0.5, "hill_biphenyl" = 0.2)


upper_bounds = c("hill_transfer" = 2.2,
                 "Km" = parameters_from_biosensor_fitting["mean","Km"] + parameters_from_biosensor_fitting["95_ci","Km"], 
                 "fold_change" = parameters_from_biosensor_fitting["mean","fold_change"] + parameters_from_biosensor_fitting["sd","fold_change"], 
                 "baseline" = parameters_from_biosensor_fitting["mean","baseline"] + parameters_from_biosensor_fitting["sd","baseline"], 
                 "slower_slope"  = parameters_from_biosensor_fitting["mean","slower_slope"] + parameters_from_biosensor_fitting["95_ci","slower_slope"], 
                 
                 "range_HipO" = 800, "HipO_constant" = 0.6, "hippurate_constant" = 500, "hill_HipO" = 2, "hill_hippurate" = 1.5,
                 "range_CocE" = 600, "CocE_constant" = 0.8, "cocaine_constant" = 100, "hill_CocE" = 2, "hill_cocaine" = 1.5,
                 
                 "range_benzamid_enz" = 400, "benzamid_enz_constant" = 5, "benzamid_constant" = 100, "hill_benzamid_enz" = 2.5, "hill_benzamid" = 1.5,
                 "range_biphenyl_enz" = 100, "biphenyl_enz_constant" = 10, "biphenyl_constant" = 100, "hill_biphenyl_enz" = 2.5, "hill_biphenyl" = 4)

F_to_optimise_everything <- function(pars_to_optimise, list_of_arguments_of_F) {
  coce_means_data = list_of_arguments_of_F$coce_means_data
  coce_sd_data = list_of_arguments_of_F$coce_sd_data
  hipo_means_data = list_of_arguments_of_F$hipo_means_data
  hipo_sd_data = list_of_arguments_of_F$hipo_sd_data
  
  biphenyl_means_data = list_of_arguments_of_F$biphenyl_means_data
  biphenyl_sd_data = list_of_arguments_of_F$biphenyl_sd_data
  benzamid_means_data = list_of_arguments_of_F$benzamid_means_data
  benzamid_sd_data = list_of_arguments_of_F$benzamid_sd_data
  
  objective_function_name = list_of_arguments_of_F$objective_function_name
  if (objective_function_name == "weighted_LSE") {
    objective_function = weighted_LSE
  } else if (objective_function_name == "relative_LSE") {
    objective_function = relative_LSE
  } else {
    objective_function = LSE
  }  
  
  pars_to_optimise = c(pars_to_optimise)
  
  # CocE coc adder
  coce_coc_adder =  data.frame(cbind(coce_means_data[,c("coc", "cocE")], calculate_coce_coc_adder_linear(parameters = pars_to_optimise, concentrations_coc = inducer_range, concentrations_CocE = enzyme_range)))
  colnames(coce_coc_adder) = c("coc", "coce", "coce_coc")
  coce_coc_objective = objective_function(prediction = coce_coc_adder, observed = coce_means_data, variable_name = "coce_coc", sd_observed = coce_sd_data)
  
  # HipO hip adder
  hipo_hip_adder =  data.frame(cbind(hipo_means_data[,c("hip", "hipo")], calculate_hipo_hip_adder_linear(parameters = pars_to_optimise, concentrations_hip =inducer_range, concentrations_HipO = enzyme_range)))
  colnames(hipo_hip_adder) = c("hip", "hipo", "hipo_hip")
  hipo_hip_objective = objective_function(prediction = hipo_hip_adder, observed = hipo_means_data, variable_name = "hipo_hip", sd_observed = hipo_sd_data)
  
  # Benzamid enzyme benzamid adder
  benzamid_enzyme_benzamid_adder =  data.frame(cbind(benzamid_means_data[,c("benzamid", "benzamid_enz")], calculate_benzamid_enz_benzamid_adder_linear(parameters = pars_to_optimise, concentrations_benzamid = inducer_range, concentrations_benzamid_enz = enzyme_range)))
  colnames(benzamid_enzyme_benzamid_adder) = c("benzamid", "benzamid_enz", 'benzamid_enz_benzamid')
  benzamid_objective = objective_function(prediction = benzamid_enzyme_benzamid_adder, observed = benzamid_means_data, variable_name = "benzamid_enz_benzamid", sd_observed = benzamid_sd_data)
  
  # Biphenyl enzyme adder
  biphenyl_enzyme_biphenyl_adder =  data.frame(cbind(biphenyl_means_data[,c("biphenyl", "biphenyl_enz")], calculate_biphenyl_enz_biphenyl_adder_linear(parameters = pars_to_optimise, concentrations_biphenyl = inducer_range, concentrations_biphenyl_enz = enzyme_range)))
  colnames(biphenyl_enzyme_biphenyl_adder) = c("biphenyl", "biphenyl_enz", 'biphenyl_enz_biphenyl')
  biphenyl_objective = objective_function(prediction = biphenyl_enzyme_biphenyl_adder, observed = biphenyl_means_data, variable_name = "biphenyl_enz_biphenyl", sd_observed = biphenyl_sd_data)
  
  current_objective = coce_coc_objective + hipo_hip_objective + benzamid_objective + biphenyl_objective
  return(current_objective)
}

run_number = 100
set.seed(42)

# Parameters from the biosensor: taken from the biosensor fitting results.
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean", "hill_transfer"], 
                                    sd = parameters_from_biosensor_fitting["se", "hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean", "Km"],
                         sd = parameters_from_biosensor_fitting["se", "Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean", "fold_change"], 
                                  sd = parameters_from_biosensor_fitting["se", "fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean", "baseline"], 
                               sd = parameters_from_biosensor_fitting["se", "baseline"])
slower_slope_starting_pars = rnorm(run_number, mean = parameters_from_biosensor_fitting["mean", "slower_slope"], 
                               sd = parameters_from_biosensor_fitting["se", "slower_slope"])


# Parameters for the enzymes
range_HipO_starting_pars = runif(run_number, min = lower_bounds["range_HipO"], max = upper_bounds["range_HipO"])
HipO_constant_starting_pars = runif(run_number, min = lower_bounds["HipO_constant"], max = upper_bounds["HipO_constant"])
hippurate_constant_starting_pars = runif(run_number, min = lower_bounds["hippurate_constant"], max = upper_bounds["hippurate_constant"])
hill_HipO_starting_pars = runif(run_number, min = lower_bounds["hill_HipO"], max = upper_bounds["hill_HipO"])
hill_hippurate_starting_pars = runif(run_number, min = lower_bounds["hill_hippurate"], max = upper_bounds["hill_hippurate"])

range_CocE_starting_pars = runif(run_number, min = lower_bounds["range_CocE"], max = upper_bounds["range_CocE"])
CocE_constant_starting_pars = runif(run_number, min = lower_bounds["CocE_constant"], max = upper_bounds["CocE_constant"])
cocaine_constant_starting_pars = runif(run_number, min = lower_bounds["cocaine_constant"], max = upper_bounds["cocaine_constant"])
hill_CocE_starting_pars = runif(run_number, min = lower_bounds["hill_CocE"], max = upper_bounds["hill_CocE"])
hill_cocaine_starting_pars = runif(run_number, min = lower_bounds["hill_cocaine"], max = upper_bounds["hill_cocaine"])

range_benzamid_enz_starting_pars = runif(run_number, min = lower_bounds["range_benzamid_enz"], max = upper_bounds["range_benzamid_enz"])
benzamid_enz_constant_starting_pars = runif(run_number, min = lower_bounds["benzamid_enz_constant"], max = upper_bounds["benzamid_enz_constant"])
benzamid_constant_starting_pars = runif(run_number, min = lower_bounds["benzamid_constant"], max = upper_bounds["benzamid_constant"])
hill_benzamid_enz_starting_pars = runif(run_number, min = lower_bounds["hill_benzamid_enz"], max = upper_bounds["hill_benzamid_enz"])
hill_benzamid_starting_pars = runif(run_number, min = lower_bounds["hill_benzamid"], max = upper_bounds["hill_benzamid"])

range_biphenyl_enz_starting_pars = runif(run_number, min = lower_bounds["range_biphenyl_enz"], max = upper_bounds["range_biphenyl_enz"])
biphenyl_enz_constant_starting_pars = runif(run_number, min = lower_bounds["biphenyl_enz_constant"], max = upper_bounds["biphenyl_enz_constant"])
biphenyl_constant_starting_pars = runif(run_number, min = lower_bounds["biphenyl_constant"], max = upper_bounds["biphenyl_constant"])
hill_biphenyl_enz_starting_pars = runif(run_number, min = lower_bounds["hill_biphenyl_enz"], max = upper_bounds["hill_biphenyl_enz"])
hill_biphenyl_starting_pars = runif(run_number, min = lower_bounds["hill_biphenyl"], max = upper_bounds["hill_biphenyl"])

all_pars = NULL
correlations = NULL

# For enzyme inducer adders
enzyme_range = c(0.1, 0.3, 1, 3, 10)
enzyme_small_range = c(0.1, 0.3, 1, 3)
inducer_range =  c(0, 10, 100, 1000)
inducer_small_range = c(0, 10, 100)

# For weighted enzyme adder
CocE_range = c(0, 0.1, 0.3, 1, 3, 10)
HipO_range = c(0, 0.1, 0.3, 1, 3, 10)
# For adder
coc_adder_range = c(0, 1, 10, 20, 100, 500, 1000)
hip_adder_range = c(0, 1, 10, 20, 100, 500, 1000)

weighted = FALSE

for (i in 1:run_number) {
  starting_pars = c("hill_transfer" = hill_transfer_starting_pars[i], 
                    "Km" = Km_starting_pars[i],  
                    "fold_change" = fold_change_starting_pars[i], 
                    "baseline" = baseline_starting_pars[i],
                    "slower_slope" = slower_slope_starting_pars[i],
                    "range_HipO" = range_HipO_starting_pars[i], 
                    "HipO_constant" = HipO_constant_starting_pars[i], 
                    "hippurate_constant" = hippurate_constant_starting_pars[i], 
                    "hill_HipO" = hill_HipO_starting_pars[i], 
                    "hill_hippurate" = hill_hippurate_starting_pars[i], 
                    "range_CocE" = range_CocE_starting_pars[i],  
                    "CocE_constant" = CocE_constant_starting_pars[i], 
                    "cocaine_constant" = cocaine_constant_starting_pars[i],  
                    "hill_CocE" = hill_CocE_starting_pars[i], 
                    "hill_cocaine" = hill_cocaine_starting_pars[i],
                    "range_benzamid_enz" = range_benzamid_enz_starting_pars[i], 
                    "benzamid_enz_constant" = benzamid_enz_constant_starting_pars[i], 
                    "benzamid_constant" = benzamid_constant_starting_pars[i], 
                    "hill_benzamid_enz" = hill_benzamid_enz_starting_pars[i], 
                    "hill_benzamid" = hill_benzamid_starting_pars[i], 
                    "range_biphenyl_enz" = range_biphenyl_enz_starting_pars[i],  
                    "biphenyl_enz_constant" = biphenyl_enz_constant_starting_pars[i], 
                    "biphenyl_constant" = biphenyl_constant_starting_pars[i],  
                    "hill_biphenyl_enz" = hill_biphenyl_enz_starting_pars[i], 
                    "hill_biphenyl" = hill_biphenyl_starting_pars[i])
  
  # print(starting_pars)
  list_of_arguments_of_F <- list(
    "hipo_means_data" = hipo_means_data, 
    "hipo_sd_data" = hipo_sd_data, 
    "coce_means_data" = coce_means_data, 
    "coce_sd_data" = coce_sd_data,  
    "benzamid_means_data" = benzamid_means_data, 
    "benzamid_sd_data" = benzamid_sd_data, 
    "biphenyl_means_data" = biphenyl_means_data, 
    "biphenyl_sd_data" = biphenyl_sd_data, 
    "objective_function_name" = "LSE"
  )
  
  exp_optim <- optim(par = starting_pars, hessian =TRUE, fn = F_to_optimise_everything, 
                     method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, 
                     list_of_arguments_of_F = list_of_arguments_of_F)
  
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

visu_biosensor= TRUE
if (visu_biosensor) {
  benzoic_acid_range = biosensor_means_data[, c("concentrations")]
  # Visualisating the results of those parameters fitting
  random_fit_biosensor = data.frame(cbind(benzoic_acid_range, apply(all_pars,1, calculate_biosensor_linear, concentrations = benzoic_acid_range)))
  unmelted_random_fit = random_fit_biosensor
  colnames(random_fit_biosensor) = c("concentrations", 1:run_number)
  random_fit_biosensor = melt(data= random_fit_biosensor, id.vars = "concentrations")
  df_observed = melt(biosensor_means_data, id.vars = "concentrations")
  sd_observed = melt(biosensor_sd_data, id.vars = "concentrations")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_biosensor, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "concentrations", y_axis = 'value', x_lab = "Benzoic acid concentration (µM)" ,y_lab ="RFU", 
                                             title = "Benzoic acid biosensor", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("biosensor_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  biosensor_mean_fit = data.frame(cbind(benzoic_acid_range, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
  colnames(biosensor_mean_fit) = c("Concentrations", "Fluorescence level")
  write.csv(biosensor_mean_fit, paste(folder_for_images, "biosensor_mean_simulation_fit.csv", sep = '/'), row.names = FALSE)
}

visu_coce_coc_adder = TRUE
if (visu_coce_coc_adder) {
  concentration_columns = coce_means_data[,c("coc", "cocE")]
  results = apply(all_pars,1, calculate_coce_coc_adder_linear, concentrations_coc = inducer_range, concentrations_CocE = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("coc", "cocE", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Cocaine", "CocE", "Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "cocaine_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising cocaine adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("cocaine_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  # Visualise the transducer at 100
  
  # Visualisating the results of those parameters fitting
  random_fit_transducer_coc = random_fit_adder[which(random_fit_adder[,c("coc")] == 100), c("cocE", 1:run_number)]
  random_fit_transducer_coc = melt(data= random_fit_transducer_coc, id.vars = "cocE")
  transducer_from_coce = coce_means_data[which(coce_means_data[,c("coc")] == 100), c("cocE", "coce_coc")]
  transducer_from_coce_sd = coce_sd_data[which(coce_sd_data[,c("coc")] == 100), c("cocE", "coce_coc")]
  
  df_observed = melt(transducer_from_coce, id.vars = "cocE")
  sd_observed = melt(transducer_from_coce_sd, id.vars = "cocE")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_transducer_coc, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "cocE", y_axis = 'value', x_lab = "CocE (nM)" ,y_lab ="RFU", 
                                             title = "Cocaine transducer from adder at 100 µM", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_cocaine_adder.jpeg", sep = '_'), path =folder_for_images)
}

visu_hipo_hip_adder = TRUE
if (visu_hipo_hip_adder) {
  concentration_columns = hipo_means_data[,c("hip", "hipo")]
  results = apply(all_pars,1, calculate_hipo_hip_adder_linear, concentrations_hip = inducer_range, concentrations_HipO = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("hip", "hipo", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Hippuric acid", "HipO", "Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "hippurate_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising hippuric acid adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("hippurate_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  # Visualisating the results of those parameters fitting
  random_fit_transducer_hip = random_fit_adder[which(random_fit_adder[,c("hip")] == 100), c("hipo", 1:run_number)]
  random_fit_transducer_hip = melt(data= random_fit_transducer_hip, id.vars = "hipo")
  transducer_from_hipo = hipo_means_data[which(hipo_means_data[,c("hip")] == 100), c("hipo", "hipo_hip")]
  transducer_from_hipo_sd = hipo_sd_data[which(hipo_sd_data[,c("hip")] == 100), c("hipo", "hipo_hip")]
  
  df_observed = melt(transducer_from_hipo, id.vars = "hipo")
  sd_observed = melt(transducer_from_hipo_sd, id.vars = "hipo")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_transducer_hip, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "hipo", y_axis = 'value', x_lab = "HipO (nM)" ,y_lab ="RFU", 
                                             title = "Hippurate transducer from adder at 100 µM", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_hippurate_adder.jpeg", sep = '_'), path =folder_for_images)
}

visu_benzamide_full = TRUE
if (visu_benzamide_full) {
  concentration_columns = benzamid_means_data[,c("benzamid", "benzamid_enz")]
  results = apply(all_pars,1, calculate_benzamid_enz_benzamid_adder_linear, concentrations_benzamid = inducer_range, concentrations_benzamid_enz = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("benzamid", "benzamid_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Benzamid", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "benzamid_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Benzamid adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("benzamid_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  # Visualisating the results of those parameters fitting
  random_fit_transducer_benzamid = random_fit_adder[which(random_fit_adder[,c("benzamid")] == 100), c("benzamid_enz", 1:run_number)]
  random_fit_transducer_benzamid = melt(data= random_fit_transducer_benzamid, id.vars = "benzamid_enz")
  transducer_from_benzamid = benzamid_means_data[which(benzamid_means_data[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")]
  transducer_from_benzamid_sd = benzamid_sd_data[which(benzamid_sd_data[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")]
  
  df_observed = melt(transducer_from_benzamid, id.vars = "benzamid_enz")
  sd_observed = melt(transducer_from_benzamid_sd, id.vars = "benzamid_enz")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_transducer_benzamid, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "benzamid_enz", y_axis = 'value', x_lab = "Benzamide enzyme (nM)" ,y_lab ="RFU", 
                                             title = "Benzamide transducer from adder at 100 µM", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_benzamide_adder.jpeg", sep = '_'), path =folder_for_images)
}

visu_benzamide_small = TRUE
if (visu_benzamide_small) {
  concentration_columns = benzamid_means_data_small[,c("benzamid", "benzamid_enz")]
  results = apply(all_pars,1, calculate_benzamid_enz_benzamid_adder_linear, concentrations_benzamid = inducer_small_range, concentrations_benzamid_enz = enzyme_small_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("benzamid", "benzamid_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Benzamid", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "benzamid_small_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Benzamid adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("benzamid_small_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  # Visualisating the results of those parameters fitting
  random_fit_transducer_benzamid = random_fit_adder[which(random_fit_adder[,c("benzamid")] == 100), c("benzamid_enz", 1:run_number)]
  random_fit_transducer_benzamid = melt(data= random_fit_transducer_benzamid, id.vars = "benzamid_enz")
  transducer_from_benzamid = benzamid_means_data_small[which(benzamid_means_data_small[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")]
  transducer_from_benzamid_sd = benzamid_sd_data_small[which(benzamid_sd_data_small[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")]
  
  df_observed = melt(transducer_from_benzamid, id.vars = "benzamid_enz")
  sd_observed = melt(transducer_from_benzamid_sd, id.vars = "benzamid_enz")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_transducer_benzamid, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "benzamid_enz", y_axis = 'value', x_lab = "Benzamide enzyme (nM)" ,y_lab ="RFU", 
                                             title = "Benzamide transducer from adder at 100 µM", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_benzamide_small_adder.jpeg", sep = '_'), path =folder_for_images)
}

visu_biphenyle_full = TRUE
if (visu_biphenyle_full) {
  concentration_columns = biphenyl_means_data[,c("biphenyl", "biphenyl_enz")]
  results = apply(all_pars,1, calculate_biphenyl_enz_biphenyl_adder_linear, concentrations_biphenyl = inducer_range, concentrations_biphenyl_enz = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("biphenyl", "biphenyl_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Biphenyl", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "biphenyl_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Biphenyl adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("biphenyl_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  # Visualisating the results of those parameters fitting
  random_fit_transducer_biphenyl = random_fit_adder[which(random_fit_adder[,c("biphenyl")] == 100), c("biphenyl_enz", 1:run_number)]
  random_fit_transducer_biphenyl = melt(data= random_fit_transducer_biphenyl, id.vars = "biphenyl_enz")
  transducer_from_biphenyl = biphenyl_means_data[which(biphenyl_means_data[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")]
  transducer_from_biphenyl_sd = biphenyl_sd_data[which(biphenyl_sd_data[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")]
  
  df_observed = melt(transducer_from_biphenyl, id.vars = "biphenyl_enz")
  sd_observed = melt(transducer_from_biphenyl_sd, id.vars = "biphenyl_enz")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_transducer_biphenyl, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "biphenyl_enz", y_axis = 'value', x_lab = "Biphenyl enzyme (nM)" ,y_lab ="RFU", 
                                             title = "Biphenyl transducer from adder at 100 µM", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_biphenyl_adder.jpeg", sep = '_'), path =folder_for_images)
}

visu_biphenyle_small = TRUE
if (visu_biphenyle_small) {
  concentration_columns = biphenyl_means_data_small[,c("biphenyl", "biphenyl_enz")]
  results = apply(all_pars,1, calculate_biphenyl_enz_biphenyl_adder_linear, concentrations_biphenyl = inducer_small_range, concentrations_biphenyl_enz = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("biphenyl", "biphenyl_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Biphenyl", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "biphenyl_small_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Biphenyl adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("biphenyl_small_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
  
  # Visualisating the results of those parameters fitting
  random_fit_transducer_biphenyl = random_fit_adder[which(random_fit_adder[,c("biphenyl")] == 100), c("biphenyl_enz", 1:run_number)]
  random_fit_transducer_biphenyl = melt(data= random_fit_transducer_biphenyl, id.vars = "biphenyl_enz")
  transducer_from_biphenyl = biphenyl_means_data[which(biphenyl_means_data[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")]
  transducer_from_biphenyl_sd = biphenyl_sd_data[which(biphenyl_sd_data[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")]
  
  df_observed = melt(transducer_from_biphenyl, id.vars = "biphenyl_enz")
  sd_observed = melt(transducer_from_biphenyl_sd, id.vars = "biphenyl_enz")
  
  plotting_random_pars_and_experiments_local(df_plot_total = random_fit_transducer_biphenyl, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                             legend_title = "Constructs", 
                                             x_axis = "biphenyl_enz", y_axis = 'value', x_lab = "Biphenyl enzyme (nM)" ,y_lab ="RFU", 
                                             title = "Biphenyl transducer from adder at 100 µM", subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_biphenyl_small_adder.jpeg", sep = '_'), path =folder_for_images)
}

visu_fixed_enzyme_adder = TRUE
if (visu_fixed_enzyme_adder) {
  concentration_columns = fixed_enzyme_adder_means_data[,c("coc", "hip")]
  results = apply(all_pars,1, calculate_hip_coc_adder_linear, concentrations_hip = hip_adder_range, concentrations_coc =coc_adder_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Cocaine", "Hippuric acid","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "fixed_enzyme_adder_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising hippuric acid - cocaine adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(µM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("fixed_enzyme_adder_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
}

visu_fixed_inducer_adder = TRUE
if (visu_fixed_inducer_adder) {
  concentration_columns = fixed_inducer_adder_means_data[,c("hipo", "coce")]
  results = apply(all_pars,1, calculate_hip_coc_enzymatic_adder_linear, concentrations_hipo = HipO_range, concentrations_coce =CocE_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  enzyme_adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(enzyme_adder_mean_fit) = c("HipO", "CocE", "Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "fixed_inducer_adder_mean_simulations_fit.csv", sep = '/'), row.names = FALSE)
  
  plot_2D_heatmap(df_means= enzyme_adder_mean_fit, df_sd = NULL, title_base =  "Visualising HipO - CocE adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(nM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("fixed_inducer_adder_mean_simulation.jpeg", sep = '_'), path =folder_for_images)
}


