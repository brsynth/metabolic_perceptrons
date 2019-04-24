# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: fit and visualise the model results on individual fitting on all experiments except biosensor - previously fitted

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
folder_for_images <- paste(current_folder, "t_linear_full_range_run_1", sep = '/')
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
  
  # CocE coc adder
  coce_coc_adder =  data.frame(cbind(coce_means_data[,c("coc", "cocE")], calculate_coce_coc_adder_linear(parameters = pars_to_optimise, concentrations_coc = inducer_range, concentrations_CocE = enzyme_range)))
  colnames(coce_coc_adder) = c("coc", "coce", "coce_coc")
  coce_coc_objective = objective_function(prediction = coce_coc_adder, observed = coce_means_data, variable_name = "coce_coc", sd_observed = coce_sd_data)
  
  # HipO hip adder
  hipo_hip_adder =  data.frame(cbind(hipo_means_data[,c("hip", "hipo")], calculate_hipo_hip_adder_linear(parameters = pars_to_optimise, concentrations_hip = inducer_range, concentrations_HipO = enzyme_range)))
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

run_number = 10
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
  
  visu_biosensor= TRUE
  if (visu_biosensor) {
    benzoic_acid_range = biosensor_means_data[, c("concentrations")]
    simulated_data =data.frame(cbind(benzoic_acid_range, calculate_biosensor_linear(parameters = pars_optimised, concentrations = benzoic_acid_range)))
    colnames(simulated_data) = c("concentrations", "biosensor")
    
    means_data = cbind(biosensor_means_data, unname(simulated_data["biosensor"]))
    sd_data = cbind(biosensor_sd_data, 0)
    
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
    saving_graph(filename = paste(i, "biosensor_model_and_data.jpeg", sep = '_'), path =folder_for_images)
  }
  
  visu_coce_coc_adder = TRUE
  if (visu_coce_coc_adder) {
    coce_coc_adder =  data.frame(cbind(coce_means_data[,c("coc", "cocE")], calculate_coce_coc_adder_linear(parameters = pars_optimised, concentrations_coc = inducer_range, concentrations_CocE = enzyme_range)))
    colnames(coce_coc_adder) = c("Cocaine", "CocE", "coce_coc")
    plot_2D_heatmap(df_means= coce_coc_adder, df_sd = NULL, title_base = "Visualising cocaine CocE adder", 
                    x_lab = "(µM)", y_lab = "(nM)",
                    saving = TRUE, name_for_saving = paste(i, "coce_coc_adder", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  }
  
  visu_hipo_hip_adder = TRUE
  if (visu_hipo_hip_adder) {
    hipo_hip_adder =  data.frame(cbind(hipo_means_data[,c("hip", "hipo")], calculate_hipo_hip_adder_linear(parameters = pars_optimised, concentrations_hip = inducer_range, concentrations_HipO = enzyme_range)))
    colnames(hipo_hip_adder) = c("Hippuric acid", "HipO", "hipo_hip")
    plot_2D_heatmap(df_means= hipo_hip_adder, df_sd = NULL, title_base = "Visualising hippurate HipO adder", saving = TRUE, 
                    x_lab = "(µM)", y_lab = "(nM)",
                    name_for_saving = paste(i, "hipo_hip_adder", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  }
  
  visu_benzamide_full = TRUE
  if (visu_benzamide_full) {
    benzamid_enzyme_benzamid_adder =  data.frame(cbind(benzamid_means_data[,c("benzamid", "benzamid_enz")], calculate_benzamid_enz_benzamid_adder_linear(parameters = pars_optimised, concentrations_benzamid = inducer_range, concentrations_benzamid_enz = enzyme_range)))
    colnames(benzamid_enzyme_benzamid_adder) = c("Benzamide", "Enzyme", 'benzamid_enz_benzamid')
    plot_2D_heatmap(df_means= benzamid_enzyme_benzamid_adder, df_sd = NULL, title_base = "Visualising Benzamide Enzyme adder", saving = TRUE, 
                    x_lab = "(µM)", y_lab = "(nM)",
                    name_for_saving = paste(i, "benzamid_full", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
    
  }
  
  visu_benzamide_small = TRUE
  if (visu_benzamide_small) {
    benzamid_enzyme_benzamid_adder =  data.frame(cbind(benzamid_means_data_small[,c("benzamid", "benzamid_enz")], calculate_benzamid_enz_benzamid_adder_linear(parameters = pars_optimised, concentrations_benzamid = inducer_small_range, concentrations_benzamid_enz = enzyme_small_range)))
    colnames(benzamid_enzyme_benzamid_adder) = c("Benzamide", "Enzyme", 'benzamid_enz_benzamid')
    plot_2D_heatmap(df_means= benzamid_enzyme_benzamid_adder, df_sd = NULL, title_base = "Visualising Benzamide Enzyme adder", saving = TRUE, 
                    x_lab = "(µM)", y_lab = "(nM)",
                    name_for_saving = paste(i, "benzamid_small", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
    
  }
  
  visu_biphenyle_full = TRUE
  if (visu_biphenyle_full) {
    biphenyl_enzyme_biphenyl_adder =  data.frame(cbind(biphenyl_means_data[,c("biphenyl", "biphenyl_enz")], calculate_biphenyl_enz_biphenyl_adder_linear(parameters = pars_optimised, concentrations_biphenyl = inducer_range, concentrations_biphenyl_enz = enzyme_range)))
    colnames(biphenyl_enzyme_biphenyl_adder) = c("Biphenyl", "Enzyme", 'biphenyl_enz_biphenyl')
    plot_2D_heatmap(df_means= biphenyl_enzyme_biphenyl_adder, df_sd = NULL, title_base = "Visualising Biphenyl Enzyme adder", saving = TRUE, 
                    x_lab = "(µM)", y_lab = "(nM)",
                    name_for_saving = paste(i, "biphenyl_full", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
    
  }
  
  visu_biphenyle_small = TRUE
  if (visu_biphenyle_small) {
    biphenyl_enzyme_biphenyl_adder =  data.frame(cbind(biphenyl_means_data_small[,c("biphenyl", "biphenyl_enz")], calculate_biphenyl_enz_biphenyl_adder_linear(parameters = pars_optimised, concentrations_biphenyl = inducer_small_range, concentrations_biphenyl_enz = enzyme_range)))
    colnames(biphenyl_enzyme_biphenyl_adder) = c("Biphenyl", "Enzyme", 'biphenyl_enz_biphenyl')
    plot_2D_heatmap(df_means= biphenyl_enzyme_biphenyl_adder, df_sd = NULL, title_base = "Visualising Biphenyl Enzyme adder", saving = TRUE, 
                    x_lab = "(µM)", y_lab = "(nM)",
                    name_for_saving = paste(i, "biphenyl_small", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
    
  }
  
  visu_fixed_enzyme_adder = TRUE
  if (visu_fixed_enzyme_adder) {
    adder_df =  data.frame(cbind(fixed_enzyme_adder_means_data[,c("coc", "hip")], calculate_hip_coc_adder_linear(parameters = pars_optimised, concentrations_hip = hip_adder_range, concentrations_coc = coc_adder_range)))
    colnames(adder_df) = c("coc", "hip", "coc_hip")
    colnames(adder_df) = c("Cocaine", "Hippuric acid", "coc_hip")
    
    plot_2D_heatmap(df_means= adder_df, df_sd = NULL, title_base = "Visualising cocaine - hippuric acid adder", 
                    saving = TRUE, x_lab = "(µM)", y_lab = "(µM)",
                    name_for_saving = paste(i, "fixed_enzyme_adder_model.jpeg", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
    saving_graph(filename = paste(i, "fixed_enzyme_adder_model.jpeg", sep = '_'), path =folder_for_images)
  }
  
  visu_fixed_inducer_adder = TRUE
  if (visu_fixed_inducer_adder) {
    fixed_inducer_adder =  data.frame(cbind(fixed_inducer_adder_means_data[,c("hipo", "coce")], calculate_hip_coc_enzymatic_adder_linear(parameters = pars_optimised, concentrations_hipo = HipO_range, concentrations_coce = CocE_range)))
    colnames(fixed_inducer_adder) = c("HipO", "CocE", "hipo_coce")
    
    plot_2D_heatmap(df_means= fixed_inducer_adder, df_sd = NULL, title_base = "Visualising HipO - CocE adder", 
                    saving = TRUE, x_lab = "(nM)", y_lab = "(nM)",
                    name_for_saving = paste(i, "fixed_inducer_adder_model.jpeg", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
    saving_graph(filename = paste(i, "fixed_inducer_adder_model.jpeg", sep = '_'), path =folder_for_images)
  }
}

all_pars_mean = data.frame(t(apply(all_pars,2,mean)))  
all_pars_sd= data.frame(t(apply(all_pars,2,sd))) 
all_pars_sem = all_pars_sd/sqrt(run_number)
all_pars_95_CI = 1.96 * all_pars_sd/sqrt(run_number)
parameters_mean_and_ci = rbind(all_pars_mean, all_pars_sd, all_pars_sem, all_pars_95_CI)

write.csv(x = correlations, paste(folder_for_images, "correlation.csv", sep = '/'), row.names = FALSE)
write.csv(x = parameters_mean_and_ci, paste(folder_for_images, "mean_and_ci.csv", sep = '/'), row.names = FALSE)


