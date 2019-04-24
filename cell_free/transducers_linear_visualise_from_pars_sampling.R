# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: visualise the model results from sampling on parameters estimation for all experiments

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

# Define where you want to store your images:

folder_for_images <- paste(current_folder, "t_linear_full_range_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

parameters_from_fitting = read.csv("t_linear_full_range_full_run_1/mean_and_ci.csv")
rownames(parameters_from_fitting) = c("mean", "sd", "se","95_ci")

# collecting_all_data

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

run_number = 100
set.seed(42)

# I copy pasted parameters from that file. Will necessitate one further round of simplification next time.

# mean = mean of fit, sd
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_transfer"], 
                                    sd = parameters_from_fitting["se","hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","Km"], 
                         sd = parameters_from_fitting["se","Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","fold_change"], 
                                  sd = parameters_from_fitting["se","fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","baseline"],
                               sd = parameters_from_fitting["se","baseline"])
slower_slope_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean", "slower_slope"], 
                                   sd = parameters_from_fitting["se", "slower_slope"])

range_HipO_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","range_HipO"], 
                                 sd = parameters_from_fitting["se","range_HipO"])
HipO_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","HipO_constant"], 
                                    sd = parameters_from_fitting["se","HipO_constant"])
hippurate_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hippurate_constant"], 
                                         sd = parameters_from_fitting["se","hippurate_constant"])
hill_HipO_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_HipO"], 
                                sd = parameters_from_fitting["se","hill_HipO"])
hill_hippurate_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_hippurate"], 
                                     sd = parameters_from_fitting["se","hill_hippurate"])

range_CocE_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","range_CocE"], 
                                 sd = parameters_from_fitting["se","range_CocE"])
CocE_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","CocE_constant"], 
                                    sd = parameters_from_fitting["se","CocE_constant"])
cocaine_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","cocaine_constant"], 
                                         sd = parameters_from_fitting["se","cocaine_constant"])
hill_CocE_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_CocE"], 
                                sd = parameters_from_fitting["se","hill_CocE"])
hill_cocaine_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_cocaine"], 
                                     sd = parameters_from_fitting["se","hill_cocaine"])

range_benzamid_enz_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","range_benzamid_enz"], 
                                         sd = parameters_from_fitting["se","range_benzamid_enz"])
benzamid_enz_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","benzamid_enz_constant"], 
                                            sd = parameters_from_fitting["se","benzamid_enz_constant"])
benzamid_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","benzamid_constant"], 
                                        sd = parameters_from_fitting["se","benzamid_constant"])
hill_benzamid_enz_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_benzamid_enz"], 
                                        sd = parameters_from_fitting["se","hill_benzamid_enz"])
hill_benzamid_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_benzamid"], 
                                    sd = parameters_from_fitting["se","hill_benzamid"])

range_biphenyl_enz_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","range_biphenyl_enz"], 
                                         sd = parameters_from_fitting["se","range_biphenyl_enz"])
biphenyl_enz_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","biphenyl_enz_constant"], 
                                            sd = parameters_from_fitting["se","biphenyl_enz_constant"])
biphenyl_constant_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","biphenyl_constant"], 
                                        sd = parameters_from_fitting["se","biphenyl_constant"])
hill_biphenyl_enz_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_biphenyl_enz"], 
                                        sd = parameters_from_fitting["se","hill_biphenyl_enz"])
hill_biphenyl_starting_pars = rnorm(run_number, mean = parameters_from_fitting["mean","hill_biphenyl"], 
                                   sd = parameters_from_fitting["se","hill_biphenyl"])


all_pars = cbind(hill_transfer_starting_pars, Km_starting_pars, fold_change_starting_pars, baseline_starting_pars, slower_slope_starting_pars,
                 range_HipO_starting_pars, HipO_constant_starting_pars, hippurate_constant_starting_pars, hill_HipO_starting_pars, hill_hippurate_starting_pars,
                 range_CocE_starting_pars, CocE_constant_starting_pars, cocaine_constant_starting_pars, hill_CocE_starting_pars,hill_cocaine_starting_pars,
                 range_benzamid_enz_starting_pars, benzamid_enz_constant_starting_pars, benzamid_constant_starting_pars, hill_benzamid_enz_starting_pars, hill_benzamid_starting_pars,
                 range_biphenyl_enz_starting_pars, biphenyl_enz_constant_starting_pars, biphenyl_constant_starting_pars, hill_biphenyl_enz_starting_pars, hill_biphenyl_starting_pars
                 )
colnames(all_pars) = c("hill_transfer", "Km", "fold_change", "baseline", "slower_slope",
                       "range_HipO", "HipO_constant", "hippurate_constant", "hill_HipO", "hill_hippurate",
                       "range_CocE", "CocE_constant", "cocaine_constant", "hill_CocE", "hill_cocaine",
                       "range_benzamid_enz", "benzamid_enz_constant", "benzamid_constant", "hill_benzamid_enz", "hill_benzamid",
                       "range_biphenyl_enz", "biphenyl_enz_constant", "biphenyl_constant", "hill_biphenyl_enz", "hill_biphenyl")

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
                                             x_axis = "concentrations", y_axis = 'value', x_lab = "Benzoic acid concentration (µM)" , y_lab ="RFU", 
                                             title = "Benzoic acid biosensor", subtitle = "", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("biosensor_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
  biosensor_mean_fit = data.frame(cbind(benzoic_acid_range, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
  colnames(biosensor_mean_fit) = c("Concentrations", "Fluorescence level")
  write.csv(biosensor_mean_fit, paste(folder_for_images, "biosensor_mean_pars_sampling_fit.csv", sep = '/'), row.names = FALSE)
}

visu_coce_coc_adder = TRUE
if (visu_coce_coc_adder) {
  concentration_columns = coce_means_data[,c("coc", "cocE")]
  results = apply(all_pars,1, calculate_coce_coc_adder_linear, concentrations_coc = inducer_range, concentrations_CocE = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("coc", "cocE", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Cocaine", "CocE", "Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "cocaine_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising cocaine adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("cocaine_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
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
                                             title = "Cocaine transducer from adder", subtitle = "Cocaine concentration of 100 µM", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_cocaine_adder_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_hipo_hip_adder = TRUE
if (visu_hipo_hip_adder) {
  concentration_columns = hipo_means_data[,c("hip", "hipo")]
  results = apply(all_pars,1, calculate_hipo_hip_adder_linear, concentrations_hip = inducer_range, concentrations_HipO = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("hip", "hipo", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Hippuric acid", "HipO", "Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "hippurate_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising hippuric acid adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("hippurate_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
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
                                             title = "Hippurate transducer from adder", subtitle = "Hippurate concentration of 100 µM", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_hippurate_adder_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_benzamide_full = TRUE
if (visu_benzamide_full) {
  concentration_columns = benzamid_means_data[,c("benzamid", "benzamid_enz")]
  results = apply(all_pars,1, calculate_benzamid_enz_benzamid_adder_linear, concentrations_benzamid = inducer_range, concentrations_benzamid_enz = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("benzamid", "benzamid_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Benzamid", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "benzamid_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Benzamid adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("benzamid_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
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
                                             title = "Benzamide transducer from adder", subtitle = "Benzamide concentration of 100 µM",  
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_benzamide_adder_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_benzamide_small = TRUE
if (visu_benzamide_small) {
  concentration_columns = benzamid_means_data_small[,c("benzamid", "benzamid_enz")]
  results = apply(all_pars,1, calculate_benzamid_enz_benzamid_adder_linear, concentrations_benzamid = inducer_small_range, concentrations_benzamid_enz = enzyme_small_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("benzamid", "benzamid_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Benzamid", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "benzamid_small_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Benzamid adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("benzamid_small_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
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
  
  saving_graph(filename = paste("transducer_from_benzamide_small_adder_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_biphenyle_full = TRUE
if (visu_biphenyle_full) {
  concentration_columns = biphenyl_means_data[,c("biphenyl", "biphenyl_enz")]
  results = apply(all_pars,1, calculate_biphenyl_enz_biphenyl_adder_linear, concentrations_biphenyl = inducer_range, concentrations_biphenyl_enz = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("biphenyl", "biphenyl_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Biphenyl", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "biphenyl_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Biphenyl adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("biphenyl_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
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
                                             title = "Biphenyl transducer from adder", subtitle = "Biphenyl concentration of 100 µM",  
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("transducer_from_biphenyl_adder_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_biphenyle_small = TRUE
if (visu_biphenyle_small) {
  concentration_columns = biphenyl_means_data_small[,c("biphenyl", "biphenyl_enz")]
  results = apply(all_pars,1, calculate_biphenyl_enz_biphenyl_adder_linear, concentrations_biphenyl = inducer_small_range, concentrations_biphenyl_enz = enzyme_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  colnames(random_fit_adder) = c("biphenyl", "biphenyl_enz", 1:run_number)
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Biphenyl", "Enzyme","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "biphenyl_small_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising Biphenyl adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("biphenyl_small_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
  
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
  
  saving_graph(filename = paste("transducer_from_biphenyl_small_adder_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_fixed_enzyme_adder = TRUE
if (visu_fixed_enzyme_adder) {
  concentration_columns = fixed_enzyme_adder_means_data[,c("coc", "hip")]
  results = apply(all_pars,1, calculate_hip_coc_adder_linear, concentrations_hip = hip_adder_range, concentrations_coc =coc_adder_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(adder_mean_fit) = c("Cocaine", "Hippuric acid","Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "fixed_enzyme_adder_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, title_base =  "Visualising hippuric acid - cocaine adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(µM)", y_lab = "(µM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("fixed_enzyme_adder_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

visu_fixed_inducer_adder = TRUE
if (visu_fixed_inducer_adder) {
  concentration_columns = fixed_inducer_adder_means_data[,c("hipo", "coce")]
  results = apply(all_pars,1, calculate_hip_coc_enzymatic_adder_linear, concentrations_hipo = HipO_range, concentrations_coce =CocE_range)
  random_fit_adder = data.frame(cbind(concentration_columns, results))
  enzyme_adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
  colnames(enzyme_adder_mean_fit) = c("HipO", "CocE", "Fluorescence level")
  write.csv(adder_mean_fit, paste(folder_for_images, "fixed_inducer_adder_mean_pars_sampling.csv", sep = '/'), row.names = FALSE)
  
  plot_2D_heatmap(df_means= enzyme_adder_mean_fit, df_sd = NULL, title_base =  "Visualising HipO - CocE adder", saving = FALSE, 
                  name_for_saving = NULL, 
                  x_lab = "(nM)", y_lab = "(nM)",
                  folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("fixed_inducer_adder_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}



