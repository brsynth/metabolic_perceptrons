# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: Score and visualise the results from simulation with the mean sampled parameters

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

parameters_from_fitting = read.csv("t_linear_full_range_full_run_1/mean_and_ci.csv")

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

benzamid_means_data_small = benzamid_means_data[which(benzamid_means_data[,c("benzamid")] %in% c(0, 10, 100)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]
benzamid_sd_data_small = benzamid_sd_data[which(benzamid_sd_data[,c("benzamid")] %in% c(0, 10, 100)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]

benzamid_means_data_small = benzamid_means_data_small[which(benzamid_means_data_small[,c("benzamid_enz")] %in% c(0.1, 0.3, 1, 3)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]
benzamid_sd_data_small = benzamid_sd_data_small[which(benzamid_sd_data_small[,c("benzamid_enz")] %in%  c(0.1, 0.3, 1, 3)), c("benzamid", "benzamid_enz", "benzamid_enz_benzamid")]

biphenyl_means_data_small = biphenyl_means_data[which(biphenyl_means_data[,c("biphenyl")] %in% c(0, 10, 100)), c("biphenyl", "biphenyl_enz", "biphenyl_enz_biphenyl")]
biphenyl_sd_data_small = biphenyl_sd_data[which(biphenyl_sd_data[,c("biphenyl")] %in% c(0, 10, 100)), c("biphenyl", "biphenyl_enz", "biphenyl_enz_biphenyl")]

benzoic_acid_range = biosensor_means_data[, c("concentrations")]
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


rownames(parameters_from_fitting) = c("mean", "sd", "se","95_ci")

mean_pars = unlist(parameters_from_fitting[c("mean"),])

# Visualise biosensor
visu_biosensor= TRUE
if (visu_biosensor) {
  benzoic_acid_range = biosensor_means_data[, c("concentrations")]
  biosensor_df = unlist(calculate_biosensor_linear(parameters = mean_pars, concentrations = benzoic_acid_range))
  means_data_visu = cbind(biosensor_means_data, biosensor_df)
  sd_data_visu = cbind(biosensor_sd_data, 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  all_data_and_model = cbind(biosensor_means_data, biosensor_sd_data[,c("biosensor")], biosensor_df)
  colnames(all_data_and_model) = c("Concentrations", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_biosensor_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
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
}

# Visualise cocaine heatmap and 100 transducer

visu_coce_coc_adder = TRUE
if (visu_coce_coc_adder) {
  coce_coc_adder =  data.frame(cbind(coce_means_data[,c("coc", "cocE")], calculate_coce_coc_adder_linear(parameters = mean_pars, concentrations_coc = inducer_range, concentrations_CocE = enzyme_range)))
  colnames(coce_coc_adder) = c("Cocaine", "CocE", "coce_coc")
  plot_2D_heatmap(df_means= coce_coc_adder, df_sd = NULL, title_base = "Visualising cocaine CocE adder", 
                  x_lab = "(µM)", y_lab = "(nM)",
                  saving = TRUE, name_for_saving = paste("coce_coc_adder_mean_of_pars", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  # Transducer
  means_data_visu = cbind(coce_means_data[which(coce_means_data[,c("coc")] == 100), c("cocE", "coce_coc")], 
                          coce_coc_adder[which(coce_coc_adder[,c("Cocaine")] == 100), c("coce_coc")])
  sd_data_visu = cbind(coce_sd_data[which(coce_sd_data[,c("coc")] == 100), c("cocE", "coce_coc")], 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  all_data_and_model = cbind(coce_means_data, coce_sd_data[,c("coce_coc")], coce_coc_adder[,c("coce_coc")])
  colnames(all_data_and_model) = c("Cocaine", "CocE", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_cocaine_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
  title = 'Cocaine transducer'
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "Cocaine concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("coc_trasnducer_from_adder_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  # Save data and calculate scores
  means_data_coc_adder_score = cbind(coce_means_data, unname(coce_coc_adder[,c("coce_coc")]))
  sd_data_coc_adder_score = cbind(coce_sd_data, 0)
  
  colnames(means_data_coc_adder_score) = c("coc", "CocE", "Data", "Model")
  colnames(sd_data_coc_adder_score) = c("coc", "CocE", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_coc_adder_score, sd_data = sd_data_coc_adder_score)
  
  RMSD_coc_adder = all_scores["RMSD"]
  correlation_coc_adder = all_scores["correlation"]
  r_squared_coc_adder = all_scores["r_squared"]
  r_squared_coc_unweighted_adder =all_scores["r_squared_unweighted"]
  
  coc_adder_table_results = c(RMSD_coc_adder, correlation_coc_adder, r_squared_coc_adder, r_squared_coc_unweighted_adder)
  write.csv(coce_coc_adder, paste(folder_for_images, "coce_coc_adder_from_mean_of_pars.csv", sep = '/'))
}

# Visualise hippurate heatmap and 100 transducer

visu_hipo_hip_adder = TRUE
if (visu_hipo_hip_adder) {
  hipo_hip_adder =  data.frame(cbind(hipo_means_data[,c("hip", "hipo")], calculate_hipo_hip_adder_linear(parameters = mean_pars, concentrations_hip = inducer_range, concentrations_HipO = enzyme_range)))
  colnames(hipo_hip_adder) = c("Hippurate", "HipO", "hipo_hip")
  plot_2D_heatmap(df_means= hipo_hip_adder, df_sd = NULL, title_base = "Visualising hippurate HipO adder", 
                  x_lab = "(µM)", y_lab = "(nM)",
                  saving = TRUE, name_for_saving = paste("hipo_hip_adder_mean_of_pars", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  # Transducer
  means_data_visu = cbind(hipo_means_data[which(hipo_means_data[,c("hip")] == 100), c("hipo", "hipo_hip")], 
                          hipo_hip_adder[which(hipo_hip_adder[,c("Hippurate")] == 100), c("hipo_hip")])
  sd_data_visu = cbind(hipo_sd_data[which(hipo_sd_data[,c("hip")] == 100), c("hipo", "hipo_hip")], 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  all_data_and_model = cbind(hipo_means_data, hipo_sd_data[,c("hipo_hip")], hipo_hip_adder[,c("hipo_hip")])
  colnames(all_data_and_model) = c("Hippurate", "HipO", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_hippurate_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  title = 'Hippuric acid transducer'
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "Hippurate concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("hip_trasnducer_from_adder_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  # Save data and calculate scores
  means_data_hip_adder_score = cbind(hipo_means_data, unname(hipo_hip_adder[,c("hipo_hip")]))
  sd_data_hip_adder_score = cbind(hipo_sd_data, 0)
  
  colnames(means_data_hip_adder_score) = c("hip", "hipo", "Data", "Model")
  colnames(sd_data_hip_adder_score) = c("hip", "hipo", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_hip_adder_score, sd_data = sd_data_hip_adder_score)
  
  RMSD_hip_adder = all_scores["RMSD"]
  correlation_hip_adder = all_scores["correlation"]
  r_squared_hip_adder = all_scores["r_squared"]
  r_squared_hip_unweighted_adder =all_scores["r_squared_unweighted"]
  
  hip_adder_table_results = c(RMSD_hip_adder, correlation_hip_adder, r_squared_hip_adder, r_squared_hip_unweighted_adder)
  write.csv(hipo_hip_adder, paste(folder_for_images, "hipo_hip_adder_from_mean_of_pars.csv", sep = '/'))
}

# Visualise benzamid heatmap and 100 transducer

visu_benzamid_adder = TRUE
if (visu_benzamid_adder) {
  benzamid_adder =  data.frame(cbind(benzamid_means_data[,c("benzamid", "benzamid_enz")], calculate_benzamid_enz_benzamid_adder_linear(parameters = mean_pars, concentrations_benzamid = inducer_range, concentrations_benzamid_enz = enzyme_range)))
  colnames(benzamid_adder) = c("Benzamide", "Enzyme", 'benzamid_enz_benzamid')
  plot_2D_heatmap(df_means= benzamid_adder, df_sd = NULL, title_base = "Visualising benzamide enzyme adder", 
                  x_lab = "(µM)", y_lab = "(nM)",
                  saving = TRUE, name_for_saving = paste("benzamid_adder_mean_of_pars", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  # Transducer
  means_data_visu = cbind(benzamid_means_data[which(benzamid_means_data[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")], 
                          benzamid_adder[which(benzamid_adder[,c("Benzamide")] == 100), c("benzamid_enz_benzamid")])
  sd_data_visu = cbind(benzamid_sd_data[which(benzamid_sd_data[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")], 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  all_data_and_model = cbind(benzamid_means_data, benzamid_sd_data[,c("benzamid_enz_benzamid")], benzamid_adder[,c("benzamid_enz_benzamid")])
  colnames(all_data_and_model) = c("Benzamide", "Enzyme", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_benzamid_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  title = 'Hippuric acid transducer'
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "Benzamid concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("benzamid_transducer_from_adder_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  # Save data and calculate scores
  means_data_benzamid_adder_score = cbind(benzamid_means_data, unname(benzamid_adder[,c("benzamid_enz_benzamid")]))
  sd_data_benzamid_adder_score = cbind(benzamid_sd_data, 0)
  
  colnames(means_data_benzamid_adder_score) = c("benzamid", "benzamid_enz", "Data", "Model")
  colnames(sd_data_benzamid_adder_score) = c("benzamid", "benzamid_enz", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_benzamid_adder_score, sd_data = sd_data_benzamid_adder_score)
  
  RMSD_benzamid_adder = all_scores["RMSD"]
  correlation_benzamid_adder = all_scores["correlation"]
  r_squared_benzamid_adder = all_scores["r_squared"]
  r_squared_benzamid_unweighted_adder =all_scores["r_squared_unweighted"]
  
  benzamid_adder_table_results = c(RMSD_benzamid_adder, correlation_benzamid_adder, r_squared_benzamid_adder, r_squared_benzamid_unweighted_adder)
  write.csv(benzamid_adder, paste(folder_for_images, "benzamid_enz_benzamid_adder_from_mean_of_pars.csv", sep = '/'))
}

# Visualise benzamid small heatmap and 100 transducer

visu_benzamid_adder_small = TRUE
if (visu_benzamid_adder_small) {
  benzamid_adder =data.frame(cbind(benzamid_means_data_small[,c("benzamid", "benzamid_enz")], calculate_benzamid_enz_benzamid_adder_linear(parameters = mean_pars, concentrations_benzamid = inducer_small_range, concentrations_benzamid_enz = enzyme_small_range)))
  colnames(benzamid_adder) = c("Benzamide", "Enzyme", 'benzamid_enz_benzamid')
  plot_2D_heatmap(df_means= benzamid_adder, df_sd = NULL, title_base = "Visualising benzamide enzyme adder", 
                  x_lab = "(µM)", y_lab = "(nM)",
                  saving = TRUE, name_for_saving = paste("benzamid_small_adder_mean_of_pars", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  # Transducer
  means_data_visu = cbind(benzamid_means_data_small[which(benzamid_means_data_small[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")], 
                          benzamid_adder[which(benzamid_adder[,c("Benzamide")] == 100), c("benzamid_enz_benzamid")])
  sd_data_visu = cbind(benzamid_sd_data_small[which(benzamid_sd_data_small[,c("benzamid")] == 100), c("benzamid_enz", "benzamid_enz_benzamid")], 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  title = 'Hippuric acid transducer'
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "Benzamid concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("benzamid_transducer_from_small_adder_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  # Save data and calculate scores
  means_data_benzamid_adder_score = cbind(benzamid_means_data_small, unname(benzamid_adder[,c("benzamid_enz_benzamid")]))
  sd_data_benzamid_adder_score = cbind(benzamid_sd_data_small, 0)
  
  colnames(means_data_benzamid_adder_score) = c("benzamid", "benzamid_enz", "Data", "Model")
  colnames(sd_data_benzamid_adder_score) = c("benzamid", "benzamid_enz", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_benzamid_adder_score, sd_data = sd_data_benzamid_adder_score)
  
  RMSD_benzamid_small_adder = all_scores["RMSD"]
  correlation_benzamid_small_adder = all_scores["correlation"]
  r_squared_benzamid_small_adder = all_scores["r_squared"]
  r_squared_benzamid_small_unweighted_adder =all_scores["r_squared_unweighted"]
  
  benzamid_small_adder_table_results = c(RMSD_benzamid_small_adder, correlation_benzamid_small_adder, r_squared_benzamid_small_adder, r_squared_benzamid_small_unweighted_adder)
  # write.csv(benzamid_adder, paste(folder_for_images, "benzamid_enz_benzamid_adder_from_mean_of_pars.csv", sep = '/'))
}

# Visu biphenyl heatmap and 100 transducer
visu_biphenyl_adder = TRUE
if (visu_biphenyl_adder) {
  biphenyl_adder =  data.frame(cbind(biphenyl_means_data[,c("biphenyl", "biphenyl_enz")], calculate_biphenyl_enz_biphenyl_adder_linear(parameters = mean_pars, concentrations_biphenyl = inducer_range, concentrations_biphenyl_enz = enzyme_range)))
  colnames(biphenyl_adder) = c("biphenyl", "Enzyme", 'biphenyl_enz_biphenyl')
  plot_2D_heatmap(df_means= biphenyl_adder, df_sd = NULL, title_base = "Visualising biphenyl enzyme adder", 
                  x_lab = "(µM)", y_lab = "(nM)",
                  saving = TRUE, name_for_saving = paste("biphenyl_adder_mean_of_pars", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  # Transducer
  means_data_visu = cbind(biphenyl_means_data[which(biphenyl_means_data[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")], 
                          biphenyl_adder[which(biphenyl_adder[,c("biphenyl")] == 100), c("biphenyl_enz_biphenyl")])
  sd_data_visu = cbind(biphenyl_sd_data[which(biphenyl_sd_data[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")], 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  all_data_and_model = cbind(biphenyl_means_data, biphenyl_sd_data[,c("biphenyl_enz_biphenyl")], biphenyl_adder[,c("biphenyl_enz_biphenyl")])
  colnames(all_data_and_model) = c("Biphenyl", "Enzyme", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_biphenyl_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  title = 'Hippuric acid transducer'
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "biphenyl concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("biphenyl_transducer_from_adder_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  # Save data and calculate scores
  means_data_biphenyl_adder_score = cbind(biphenyl_means_data, unname(biphenyl_adder[,c("biphenyl_enz_biphenyl")]))
  sd_data_biphenyl_adder_score = cbind(biphenyl_sd_data, 0)
  
  colnames(means_data_biphenyl_adder_score) = c("biphenyl", "biphenyl_enz", "Data", "Model")
  colnames(sd_data_biphenyl_adder_score) = c("biphenyl", "biphenyl_enz", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_biphenyl_adder_score, sd_data = sd_data_biphenyl_adder_score)
  
  RMSD_biphenyl_adder = all_scores["RMSD"]
  correlation_biphenyl_adder = all_scores["correlation"]
  r_squared_biphenyl_adder = all_scores["r_squared"]
  r_squared_biphenyl_unweighted_adder =all_scores["r_squared_unweighted"]
  
  biphenyl_adder_table_results = c(RMSD_biphenyl_adder, correlation_biphenyl_adder, r_squared_biphenyl_adder, r_squared_biphenyl_unweighted_adder)
  write.csv(biphenyl_adder, paste(folder_for_images, "biphenyl_enz_biphenyl_adder_from_mean_of_pars.csv", sep = '/'))
} 

# Visualise biphenyl small heatmap and 100 transducer

visu_biphenyl_adder_small = TRUE
if (visu_biphenyl_adder_small) {
  biphenyl_adder =data.frame(cbind(biphenyl_means_data_small[,c("biphenyl", "biphenyl_enz")], calculate_biphenyl_enz_biphenyl_adder_linear(parameters = mean_pars, concentrations_biphenyl = inducer_small_range, concentrations_biphenyl_enz = enzyme_range)))
  colnames(biphenyl_adder) = c("biphenyl", "Enzyme", 'biphenyl_enz_biphenyl')
  plot_2D_heatmap(df_means= biphenyl_adder, df_sd = NULL, title_base = "Visualising biphenyl enzyme adder", 
                  x_lab = "(µM)", y_lab = "(nM)",
                  saving = TRUE, name_for_saving = paste("biphenyl_small_adder_mean_of_pars", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  # Transducer
  means_data_visu = cbind(biphenyl_means_data_small[which(biphenyl_means_data_small[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")], 
                          biphenyl_adder[which(biphenyl_adder[,c("biphenyl")] == 100), c("biphenyl_enz_biphenyl")])
  sd_data_visu = cbind(biphenyl_sd_data_small[which(biphenyl_sd_data_small[,c("biphenyl")] == 100), c("biphenyl_enz", "biphenyl_enz_biphenyl")], 0)
  colnames(means_data_visu) = c("concentrations", "Data", "Model")
  colnames(sd_data_visu) = c("concentrations", "Data", "Model")
  
  df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
  sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
  title = 'Hippuric acid transducer'
  print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                        x_axis = "concentrations",y_axis = "value", x_lab = "biphenyl concentration (µM)" ,y_lab = expression('RFU'), title = title, subtitle = "", 
                                        names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
  saving_graph(filename = paste("biphenyl_transducer_from_small_adder_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  # Save data and calculate scores
  means_data_biphenyl_adder_score = cbind(biphenyl_means_data_small, unname(biphenyl_adder[,c("biphenyl_enz_biphenyl")]))
  sd_data_biphenyl_adder_score = cbind(biphenyl_sd_data_small, 0)
  
  colnames(means_data_biphenyl_adder_score) = c("biphenyl", "biphenyl_enz", "Data", "Model")
  colnames(sd_data_biphenyl_adder_score) = c("biphenyl", "biphenyl_enz", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_biphenyl_adder_score, sd_data = sd_data_biphenyl_adder_score)
  
  RMSD_biphenyl_small_adder = all_scores["RMSD"]
  correlation_biphenyl_small_adder = all_scores["correlation"]
  r_squared_biphenyl_small_adder = all_scores["r_squared"]
  r_squared_biphenyl_small_unweighted_adder =all_scores["r_squared_unweighted"]
  
  biphenyl_small_adder_table_results = c(RMSD_biphenyl_small_adder, correlation_biphenyl_small_adder, r_squared_biphenyl_small_adder, r_squared_biphenyl_small_unweighted_adder)
  # write.csv(biphenyl_adder, paste(folder_for_images, "biphenyl_enz_biphenyl_adder_from_mean_of_pars.csv", sep = '/'))
}

# Visu fixed enzyme adder

visu_fixed_enzyme_adder = TRUE
if (visu_fixed_enzyme_adder) {
  adder_df =  data.frame(cbind(fixed_enzyme_adder_means_data[,c("coc", "hip")], calculate_hip_coc_adder_linear(parameters = mean_pars, concentrations_hip = hip_adder_range, concentrations_coc = coc_adder_range)))
  colnames(adder_df) = c("coc", "hip", "coc_hip")
  colnames(adder_df) = c("Cocaine", "Hippuric acid", "coc_hip")
  
  all_data_and_model = cbind(fixed_enzyme_adder_means_data, fixed_enzyme_adder_sd_data[,c("coc_hip")], adder_df[,c("coc_hip")])
  colnames(all_data_and_model) = c("Cocaine", "Hippurate", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_fixed_enzyme_adder_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
  write.csv(adder_df, paste(folder_for_images, "fixed_enzyme_adder_from_mean_of_pars.csv", sep = '/'))
  plot_2D_heatmap(df_means= adder_df, df_sd = NULL, title_base = "Visualising cocaine - hippuric acid adder", 
                  saving = TRUE, x_lab = "(µM)", y_lab = "(µM)",
                  name_for_saving = paste("fixed_enzyme_adder_model_from_mean_of_pars.jpeg", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("fixed_enzyme_adder_model_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  
  means_data_adder_score = cbind(fixed_enzyme_adder_means_data, unname(adder_df[,c("coc_hip")]))
  sd_data_adder_score = cbind(fixed_enzyme_adder_sd_data, 0)
  
  colnames(means_data_adder_score) = c("coc", "hip", "Data", "Model")
  colnames(sd_data_adder_score) = c("coc", "hip", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_adder_score, sd_data = sd_data_adder_score)
  
  RMSD_adder = all_scores["RMSD"]
  correlation_adder = all_scores["correlation"]
  r_squared_adder = all_scores["r_squared"]
  r_squared_unweighted_adder =all_scores["r_squared_unweighted"]
  within_SEM_adder =all_scores["within_SEM"]
  
  fixed_enzyme_adder_table_results = c(RMSD_adder, correlation_adder, r_squared_adder, r_squared_unweighted_adder)
}

# Visu weighted adder

visu_fixed_inducer_adder = TRUE
if (visu_fixed_inducer_adder) {
  fixed_inducer_adder =  data.frame(cbind(fixed_inducer_adder_means_data[,c("hipo", "coce")], calculate_hip_coc_enzymatic_adder_linear(parameters = mean_pars, concentrations_hipo = HipO_range, concentrations_coce = CocE_range)))
  colnames(fixed_inducer_adder) = c("HipO", "CocE", "hipo_coce")
  
  all_data_and_model = cbind(fixed_inducer_adder_means_data, fixed_inducer_adder_sd_data[,c("coce_hipo")], fixed_inducer_adder[,c("hipo_coce")])
  colnames(all_data_and_model) = c("HipO", "CocE", "Data Means", "Data Sd", "Model")
  write.csv(all_data_and_model,paste(folder_for_images, "all_fixed_inducer_adder_from_mean_of_pars.csv", sep = '/'), row.names = FALSE)
  
  write.csv(fixed_inducer_adder, paste(folder_for_images, "fixed_inducer_adder_from_mean_of_pars.csv", sep = '/'))
  plot_2D_heatmap(df_means= fixed_inducer_adder, df_sd = NULL, title_base = "Visualising CocE - HipO adder", 
                  saving = TRUE, x_lab = "(µM)", y_lab = "(µM)",
                  name_for_saving = paste("fixed_inducer_adder_model_from_mean_of_pars.jpeg", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)
  saving_graph(filename = paste("fixed_inducer_adder_model_from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)
  
  means_data_weighted_adder_score = cbind(fixed_inducer_adder_means_data, unname(fixed_inducer_adder[,c("hipo_coce")]))
  sd_data_weighted_adder_score = cbind(fixed_inducer_adder_sd_data, 0)
  
  colnames(means_data_weighted_adder_score) = c("hipo", "coce", "Data", "Model")
  colnames(sd_data_weighted_adder_score) = c("hipo", "coce", "Data", "Model")
  
  all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_weighted_adder_score, sd_data = sd_data_weighted_adder_score)
  
  RMSD_fixed_inducer_adder = all_scores["RMSD"]
  correlation_fixed_inducer_adder = all_scores["correlation"]
  r_squared_fixed_inducer_adder = all_scores["r_squared"]
  r_squared_unweighted_fixed_inducer_adder =all_scores["r_squared_unweighted"]
  
  fixed_inducer_adder_table_results = c(RMSD_fixed_inducer_adder, correlation_fixed_inducer_adder, 
                                        r_squared_fixed_inducer_adder, r_squared_unweighted_fixed_inducer_adder)
}


table_results_dataframe = cbind(biosensor_table_results, coc_adder_table_results, hip_adder_table_results, benzamid_adder_table_results,
                                benzamid_small_adder_table_results, biphenyl_adder_table_results, biphenyl_small_adder_table_results,
                                fixed_enzyme_adder_table_results, fixed_inducer_adder_table_results)


rownames(table_results_dataframe) = c("RMSD", "Correlation", "R_squared", "Unweighted R squared")

write.csv(table_results_dataframe, paste(folder_for_images, "table_results_correct_order.csv", sep = '/'), row.names = TRUE, quote = FALSE)
