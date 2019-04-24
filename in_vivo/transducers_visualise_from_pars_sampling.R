# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: visualise the model results from sampling on parameters estimation

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

run_number = 100
set.seed(42)

parameters_from_trasnducers_fitting = read.csv("transducers_full_run_1/mean_and_ci.csv")
rownames(parameters_from_trasnducers_fitting) = c("mean", "sd", "se" ,"95_ci")

# Parameters from the biosensor: taken from the biosensor fitting results.
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "hill_transfer"], 
                                    sd = parameters_from_trasnducers_fitting["se", "hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "Km"],
                         sd = parameters_from_trasnducers_fitting["se", "Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "fold_change"], 
                                  sd = parameters_from_trasnducers_fitting["se", "fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "baseline"], 
                               sd = parameters_from_trasnducers_fitting["se", "baseline"])

# Parameters for enzymes
range_BenZ_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "range_BenZ"], 
                                 sd = parameters_from_trasnducers_fitting["se", "range_BenZ"])
range_HipO_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "range_HipO"], 
                                 sd = parameters_from_trasnducers_fitting["se", "range_HipO"])

# Parameters for resource competition: "total_enzyme" = 10, "ratio_hip_benz" = 5, "cooperativity_resource" = 2.2, "range_resource" = 5
total_enzyme_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "total_enzyme"], 
                                   sd = parameters_from_trasnducers_fitting["se", "total_enzyme"])
ratio_hip_benz_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "ratio_hip_benz"], 
                                     sd = parameters_from_trasnducers_fitting["se", "ratio_hip_benz"])
cooperativity_resource_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "cooperativity_resource"], 
                                             sd = parameters_from_trasnducers_fitting["se", "cooperativity_resource"])
range_resource_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "range_resource"], 
                                     sd = parameters_from_trasnducers_fitting["se", "range_resource"])

all_pars = cbind(hill_transfer_starting_pars, Km_starting_pars, fold_change_starting_pars, baseline_starting_pars,
                 range_BenZ_starting_pars, range_HipO_starting_pars, 
                 total_enzyme_starting_pars, ratio_hip_benz_starting_pars, cooperativity_resource_starting_pars, range_resource_starting_pars)
colnames(all_pars) = c("hill_transfer", "Km", "fold_change", "baseline", 'range_BenZ', 'range_HipO', 
                       'total_enzyme', "ratio_hip_benz", "cooperativity_resource", "range_resource")

concentrations_adder = c(0, 1, 10, 20, 100, 500, 1000)
concentrations_transducers = c(0, 1, 10, 20, 100, 200, 500, 1000)

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


saving_graph(filename = paste("hippurate_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
hippurate_mean_fit = data.frame(cbind(concentrations_transducers, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(hippurate_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(hippurate_mean_fit, "hippurate_mean_pars_sampling.csv", row.names = FALSE)

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


saving_graph(filename = paste("benzaldehyde_pars_sampling.jpeg", sep = '_'), path =folder_for_images)
benzal_mean_fit = data.frame(cbind(concentrations_transducers, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(benzal_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(benzal_mean_fit, "benzal_mean_pars_sampling.csv", row.names = FALSE)

# Calculate and visualise adder

concentration_columns = means_data_adder[,c("hip", "benzal")]
results = apply(all_pars,1, calculate_hip_benz_adder, concentrations_benzal = concentrations_adder, concentrations_hip =concentrations_adder)

random_fit_adder = data.frame(cbind(concentration_columns, results))

adder_mean_fit = data.frame(cbind(concentration_columns, means_random_fit = apply(random_fit_adder[,3:(run_number +2)],1, mean)))
colnames(adder_mean_fit) = c("Hippurate concentrations", "Benzaldehyde concentrations","Fluorescence level")
write.csv(adder_mean_fit, "adder_mean_pars_sampling.csv", row.names = FALSE)

colnames(adder_mean_fit) = c("Hippuric acid", "Benzaldehyde","Fluorescence level")
plot_2D_heatmap(df_means= adder_mean_fit, df_sd = NULL, 
                x_lab = "(µM)", y_lab = "(µM)", title_base = "Visualising hippuric acid - benzaldehyde adder", saving = FALSE, 
                name_for_saving = NULL, folder_for_concentration_images_strategy = folder_for_images)
saving_graph(filename = paste("adder_heatmap_mean_pars_sampling.jpeg", sep = '_'), path =folder_for_images)


# Visualise individual rows and columns from the 2D heatmap

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
                                             x_axis = "benzal", y_axis = 'value', x_lab = "Benzoic acid concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                             title = paste("Trasnducer of benzaldehyde with hippurate =", hip), subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("hip", hip, "pars_sampling.jpeg", sep = '_'), path =folder_for_images)
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
                                             x_axis = "hip", y_axis = 'value', x_lab = "Benzoic acid concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                             title = paste("Transducer of hippurate with benzaldehyde =", benzal), subtitle = "Comparing data and 100 best fits", 
                                             xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")
  
  saving_graph(filename = paste("benzal", benzal, "pars_sampling.jpeg", sep = '_'), path =folder_for_images)
}

