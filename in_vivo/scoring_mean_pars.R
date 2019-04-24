# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: Score and visualise the results from simulation with the mean sampled parameters

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
folder_for_images <- paste(current_folder, "scoring_full_run_1", sep = '/')
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

means_data_coc = means_data_total[, c("concentrations", "coc")]
sd_data_coc =  sd_data_total[, c("concentrations", "coc")]

means_data_biosensor = means_data_total[, c("concentrations", "biosensor")]
sd_data_biosensor =  sd_data_total[, c("concentrations", "biosensor")]
benzoic_acid_range = means_data_biosensor[, c("concentrations")]

benzoic_acid_range = c(0, 1, 10, 20, 100, 200, 500, 1000)
concentrations_adder = c(0, 1, 10, 20, 100, 500, 1000)
concentrations_transducers = c(0, 1, 10, 20, 100, 200, 500, 1000)

parameters_from_trasnducers_fitting = read.csv("transducers_full_run_1/mean_and_ci.csv")
rownames(parameters_from_trasnducers_fitting) = c("mean", "sd", "se" ,"95_ci")

parameters_from_cocaine_fitting = read.csv("cocaine_full_run_1/mean_and_ci.csv")
rownames(parameters_from_cocaine_fitting) = c("mean", "sd", "se" ,"95_ci")

mean_pars = c(unlist(parameters_from_trasnducers_fitting["mean",]), 
              "range_CocE" = parameters_from_cocaine_fitting["mean", "range_CocE"])

# Visualise biosensor

biosensor_df = calculate_biosensor(parameters = mean_pars, concentrations = benzoic_acid_range)

means_data_visu = cbind(means_data_biosensor, biosensor_df)
sd_data_visu = cbind(sd_data_biosensor, 0)

colnames(means_data_visu) = c("concentrations", "Data", "Model")
all_data_and_model = cbind(means_data_biosensor, sd_data_biosensor[,c("biosensor")], biosensor_df)
colnames(all_data_and_model) = c("Concentrations", "Data Mean", "Data Sd", "Model")
write.csv(all_data_and_model, "all_biosensor_from_mean_of_pars.csv", row.names = FALSE)

colnames(sd_data_visu) = c("concentrations", "Data", "Model")
write.csv(means_data_visu, "biosensor_from_mean_of_pars.csv")
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
                                      x_axis = "concentrations",y_axis = "value", x_lab = "Benzoic acid concentration (µM)" ,y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = "", 
                                      names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)

saving_graph(filename = paste("biosensor", "from_mean_of_pars.jpeg", sep = '_'), path =folder_for_images)

# Visualise hippurate
  
hip_df = unlist(calculate_hip_transducer(parameters = mean_pars, concentrations_hip = concentrations_transducers))
means_data_hip_visu = cbind(means_data_hip, hip_df)
sd_data_hip_visu = cbind(sd_data_hip, 0)

colnames(means_data_hip_visu) = c("concentrations", "Data", "Model")
colnames(sd_data_hip_visu) = c("concentrations", "Data", "Model")
all_data_and_model = cbind(means_data_hip, sd_data_hip[,c("hip")], hip_df)
colnames(all_data_and_model) = c("Concentrations", "Data Means", "Data Sd", "Model")
write.csv(all_data_and_model, "all_hippurate_from_mean_of_pars.csv", row.names = FALSE)

write.csv(means_data_hip_visu, "hippurate_from_mean_of_pars.csv")

all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_hip_visu, sd_data = sd_data_hip_visu)

RMSD_hip = all_scores["RMSD"]
correlation_hip = all_scores["correlation"]
r_squared_hip = all_scores["r_squared"]
r_squared_unweighted_hip =all_scores["r_squared_unweighted"]

hip_table_results = c(RMSD_hip, correlation_hip, r_squared_hip, r_squared_unweighted_hip)

title= paste("Modeling the hippurate transducer", correlation_hip, sep =' ')    
df_observed <- melt(data= means_data_hip_visu, id.vars = "concentrations")
sd_observed = melt(data= sd_data_hip_visu, id.vars = "concentrations")

print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                      x_axis = "concentrations",y_axis = "value", x_lab = "Hippuric acid concentration (µM)" ,y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = "", 
                                      names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
folder_to_save = folder_for_images
saving_graph(filename = paste("hip_transducer_from_mean_of_pars.jpeg", sep = '_'), path =folder_to_save)  
  
# Benzaldehyde transducer

ben_df = unlist(calculate_benz_transducer(parameters = mean_pars, concentrations_ben = concentrations_transducers))
means_data_ben_visu = cbind(means_data_benzal, ben_df)
sd_data_ben_visu = cbind(sd_data_benzal, 0)

colnames(means_data_ben_visu) = c("concentrations", "Data", "Model")
colnames(sd_data_ben_visu) = c("concentrations", "Data", "Model")
all_data_and_model = cbind(means_data_benzal, sd_data_benzal[,c("benzal")], ben_df)
colnames(all_data_and_model) = c("Concentrations", "Data Means", "Data Sd", "Model")
write.csv(all_data_and_model, "all_benzal_from_mean_of_pars.csv", row.names = FALSE)

write.csv(means_data_ben_visu, "benzaldehyde_from_mean_of_pars.csv")

all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_ben_visu, sd_data = sd_data_ben_visu)

RMSD_ben = all_scores["RMSD"]
correlation_ben = all_scores["correlation"]
r_squared_ben = all_scores["r_squared"]
r_squared_unweighted_ben =all_scores["r_squared_unweighted"]

ben_table_results = c(RMSD_ben, correlation_ben, r_squared_ben, r_squared_unweighted_ben)

title= paste("Modeling the benzaldehyde transducer", correlation_ben, sep =' ')    
df_observed <- melt(data= means_data_ben_visu, id.vars = "concentrations")
sd_observed = melt(data= sd_data_ben_visu, id.vars = "concentrations")

print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                      x_axis = "concentrations",y_axis = "value", x_lab = "Benzaldehyde concentration (µM)" ,y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = "", 
                                      names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
folder_to_save = folder_for_images
saving_graph(filename = paste("ben_transducers_from_mean_pars.jpeg", sep = '_'), path =folder_to_save)  

# Visualise cocaine

coc_df = unlist(calculate_coc_transducer(parameters = mean_pars, concentrations_coc = concentrations_transducers))
means_data_coc_visu = cbind(means_data_coc, coc_df)
sd_data_coc_visu = cbind(sd_data_coc, 0)

colnames(means_data_coc_visu) = c("concentrations", "Data", "Model")
colnames(sd_data_coc_visu) = c("concentrations", "Data", "Model")
all_data_and_model = cbind(means_data_coc, sd_data_coc[,c("coc")], coc_df)
colnames(all_data_and_model) = c("Concentrations", "Data Means", "Data Sd", "Model")
write.csv(all_data_and_model, "all_cocaine_from_mean_of_pars.csv", row.names = FALSE)

write.csv(means_data_coc_visu, "cocaine_from_mean_of_pars.csv")
all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_coc_visu, sd_data = sd_data_coc_visu)

RMSD_coc = all_scores["RMSD"]
correlation_coc = all_scores["correlation"]
r_squared_coc = all_scores["r_squared"]
r_squared_unweighted_coc =all_scores["r_squared_unweighted"]

coc_table_results = c(RMSD_coc, correlation_coc, r_squared_coc, r_squared_unweighted_coc)

title= paste("Modeling the cocaine transducer", correlation_coc, sep =' ')    
df_observed <- melt(data= means_data_coc_visu, id.vars = "concentrations")
sd_observed = melt(data= sd_data_coc_visu, id.vars = "concentrations")

print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                      x_axis = "concentrations",y_axis = "value", x_lab = "Cocaine concentration (µM)" ,y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = "", 
                                      names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
folder_to_save = folder_for_images
saving_graph(filename = paste("coc_transducers_from_mean_pars.jpeg", sep = '_'), path =folder_to_save)  



# Adder

adder_df =  cbind(means_data_adder[,c("hip", "benzal")], unlist(calculate_hip_benz_adder(parameters = mean_pars, concentrations_benzal = concentrations_adder, concentrations_hip = concentrations_adder)))
colnames(adder_df) = c("Hippuric acid", "Benzaldehyde", "hip_benzal_adder")

plot_2D_heatmap(df_means= adder_df, df_sd = NULL, title_base = "Visualising adder 2D heatmap", saving = TRUE, 
                x_lab = "(µM)", y_lab = "(µM)",
                name_for_saving = paste("adder", "from_mean_of_pars.jpeg", sep = '_'), folder_for_concentration_images_strategy = folder_for_images)

colnames(adder_df) = c("hip", "benzal", "hip_benzal_adder")

all_data_and_model = cbind(means_data_adder, sd_data_adder[,c("hip_benzal_adder")], adder_df[,c("hip_benzal_adder")])
colnames(all_data_and_model) = c("Hippurate", "Benzaldehyde" , "Data Means", "Data Sd", "Model")
write.csv(all_data_and_model, "all_adder_from_mean_of_pars.csv", row.names = FALSE)

write.csv(adder_df, "adder_from_mean_of_pars.csv")

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
    
    title= paste("Modeling the adder with constant hippurate at", hip, sep =' ')
    df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
    sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
    
    print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, 
                                          variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                          x_axis = "concentrations",
                                          y_axis = "value", x_lab = "Benzaldehyde concentration (µM)" ,
                                          y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = "",
                                          names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
    folder_to_save = folder_for_images
    saving_graph(filename = paste("hip", hip, "from_mean_of_pars.jpeg", sep = '_'), path =folder_to_save)
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
    
    title= paste("Modeling the adder with constant benzaldehyde at", benzal, sep =' ')
    df_observed <- melt(data= means_data_visu, id.vars = "concentrations")
    sd_observed = melt(data= sd_data_visu, id.vars = "concentrations")
    
    print_same_graph_error_bars_biosensor(df = df_observed, df_errors = sd_observed, df_observed= NULL, df_sd_observed = NULL, 
                                          variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                          x_axis = "concentrations",
                                          y_axis = "value", x_lab = "Hippuric acid concentration (µM)" ,
                                          y_lab = expression('Fold change of GFP/OD'[600]*' (AU)'), title = title, subtitle = "",
                                          names_titles = c("Data", "Model"), legends = NULL, xlog =TRUE, print = TRUE, name_condition = "", legend = "bottom", label_names = NULL)
    folder_to_save = folder_for_images
    saving_graph(filename = paste("benzal", benzal, "from_mean_of_pars.jpeg", sep = '_'), path =folder_to_save)
  }
}


means_data_adder_score = cbind(means_data_adder, unname(adder_df[,c("hip_benzal_adder")]))
sd_data_adder_score = cbind(sd_data_adder, 0)

colnames(means_data_adder_score) = c("hip", "benzal", "Data", "Model")
colnames(sd_data_adder_score) = c("hip", "benzal", "Data", "Model")

all_scores = calculate_all_score_pred_true_true_sd(means_data = means_data_adder_score, sd_data = sd_data_adder_score)

RMSD_adder = all_scores["RMSD"]
correlation_adder = all_scores["correlation"]
r_squared_adder = all_scores["r_squared"]
r_squared_unweighted_adder =all_scores["r_squared_unweighted"]
within_SEM_adder =all_scores["within_SEM"]

adder_table_results = c(RMSD_adder, correlation_adder, r_squared_adder, r_squared_unweighted_adder)
# Saving all calculated scores 

table_results_dataframe = cbind(biosensor_table_results, ben_table_results, hip_table_results, 
                                coc_table_results,
                                adder_table_results)

# table_results_dataframe = cbind(biosensor_table_results, ben_table_results, coc_table_results, hip_table_results, 
#                                 benzal_from_trans_table_results, hip_from_trans_table_results, adder_table_results,
#                                 adder_table_results_unfit)

rownames(table_results_dataframe) = c("RMSD", "Correlation", "R_squared", "Unweighted R squared")

write.table(table_results_dataframe, "table_results_correct_order.csv", row.names = TRUE, quote = FALSE, sep = ";")