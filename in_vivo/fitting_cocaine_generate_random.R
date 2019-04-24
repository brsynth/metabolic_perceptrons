# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: fit and visualise the model results on cocaine transducer fitting for 100 parameters 
# to obtain parameters confidenc eintervals

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
folder_for_images <- paste(current_folder, "cocaine_full_run_1", sep = '/')
if (!file.exists(folder_for_images)) {
  dir.create(folder_for_images)
}

# Collecting all useful data

means_data_total = read.csv(paste(folder_for_data, "transducers", "_", "means_data", '.csv', sep = ''))
sd_data_total = read.csv(paste(folder_for_data, "transducers", "_", "sd_data", '.csv', sep = ''))

means_data_coc = means_data_total[, c("concentrations", "coc")]
sd_data_coc =  sd_data_total[, c("concentrations", "coc")]

# Parameters from biosensor and resource competition will be taken from the 100 random transducers

parameters_from_trasnducers_fitting = read.csv("transducers_full_run_1/mean_and_ci.csv")
rownames(parameters_from_trasnducers_fitting) = c("mean", "sd", "se" ,"95_ci")


lower_bounds = c("hill_transfer" = parameters_from_trasnducers_fitting["mean","hill_transfer"] - parameters_from_trasnducers_fitting["95_ci","hill_transfer"], 
                 "Km" = parameters_from_trasnducers_fitting["mean","Km"] - parameters_from_trasnducers_fitting["95_ci","Km"], 
                 "fold_change" = parameters_from_trasnducers_fitting["mean","fold_change"] - parameters_from_trasnducers_fitting["95_ci","fold_change"], 
                 "baseline" = parameters_from_trasnducers_fitting["mean","baseline"] - parameters_from_trasnducers_fitting["95_ci","baseline"], 
                 "range_CocE" = 0.1,
                 "total_enzyme" = parameters_from_trasnducers_fitting["mean","total_enzyme"] - parameters_from_trasnducers_fitting["95_ci","total_enzyme"],  
                 "ratio_hip_benz" = parameters_from_trasnducers_fitting["mean","ratio_hip_benz"] - parameters_from_trasnducers_fitting["95_ci","ratio_hip_benz"], 
                 "cooperativity_resource" = parameters_from_trasnducers_fitting["mean","cooperativity_resource"] - parameters_from_trasnducers_fitting["95_ci","cooperativity_resource"], 
                 "range_resource" = parameters_from_trasnducers_fitting["mean","range_resource"] - parameters_from_trasnducers_fitting["95_ci","range_resource"]
)


upper_bounds = c("hill_transfer" = parameters_from_trasnducers_fitting["mean","hill_transfer"] + parameters_from_trasnducers_fitting["95_ci","hill_transfer"], 
                 "Km" = parameters_from_trasnducers_fitting["mean","Km"] + parameters_from_trasnducers_fitting["95_ci","Km"], 
                 "fold_change" = parameters_from_trasnducers_fitting["mean","fold_change"] + parameters_from_trasnducers_fitting["95_ci","fold_change"], 
                 "baseline" = parameters_from_trasnducers_fitting["mean","baseline"] + parameters_from_trasnducers_fitting["95_ci","baseline"], 
                 "range_CocE" = 1.2,
                 "total_enzyme" = parameters_from_trasnducers_fitting["mean","total_enzyme"] + parameters_from_trasnducers_fitting["95_ci","total_enzyme"],  
                 "ratio_hip_benz" = parameters_from_trasnducers_fitting["mean","ratio_hip_benz"] + parameters_from_trasnducers_fitting["95_ci","ratio_hip_benz"], 
                 "cooperativity_resource" = parameters_from_trasnducers_fitting["mean","cooperativity_resource"] + parameters_from_trasnducers_fitting["95_ci","cooperativity_resource"], 
                 "range_resource" = parameters_from_trasnducers_fitting["mean","range_resource"] + parameters_from_trasnducers_fitting["95_ci","range_resource"]
)

coce_lower_bound = c("range_CocE" = 0.1)
coce_upper_bound = c("range_CocE" = 1.2)

F_to_optimise_cocaine <- function(pars_to_optimise, list_of_arguments_of_F) {
  means_data_coc <- list_of_arguments_of_F$means_data_coc
  sd_data_coc <- list_of_arguments_of_F$sd_data_coc
  variable_name_coc = list_of_arguments_of_F$variable_name_coc
  fixed_pars = list_of_arguments_of_F$fixed_pars
  
  objective_function_name = list_of_arguments_of_F$objective_function_name
  if (objective_function_name == "weighted_LSE") {
    objective_function = weighted_LSE
  } else if (objective_function_name == "relative_LSE") {
    objective_function = relative_LSE
  } else {
    objective_function = LSE
  }  
  
  coc_df = data.frame(cbind(concentrations_transducers, calculate_coc_transducer(parameters = c(fixed_pars, pars_to_optimise), concentrations_coc = concentrations_transducers)))
  colnames(coc_df) = c("concentrations", variable_name_coc)
  coc_objective = objective_function(prediction = coc_df, observed = means_data_coc, variable_name = variable_name_coc, sd_observed = sd_data_coc)
  
  current_objective = coc_objective
    
  return(current_objective)
}

run_number = 100

# Parameters from the biosensor: taken from the biosensor fitting results.
hill_transfer_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "hill_transfer"], 
                                    sd = parameters_from_trasnducers_fitting["se", "hill_transfer"])
Km_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "Km"],
                         sd = parameters_from_trasnducers_fitting["se", "Km"])
fold_change_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "fold_change"], 
                                  sd = parameters_from_trasnducers_fitting["se", "fold_change"])
baseline_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "baseline"], 
                               sd = parameters_from_trasnducers_fitting["se", "baseline"])

# Parameters for the enzyme
range_CocE_starting_pars = runif(run_number, min = lower_bounds["range_CocE"], max = upper_bounds["range_CocE"])

# Parameters for resource competition:
total_enzyme_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "total_enzyme"], 
                                   sd = parameters_from_trasnducers_fitting["se", "total_enzyme"])
ratio_hip_benz_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "ratio_hip_benz"], 
                                     sd = parameters_from_trasnducers_fitting["se", "ratio_hip_benz"])
cooperativity_resource_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "cooperativity_resource"], 
                                             sd = parameters_from_trasnducers_fitting["se", "cooperativity_resource"])
range_resource_starting_pars = rnorm(run_number, mean = parameters_from_trasnducers_fitting["mean", "range_resource"], 
                                     sd = parameters_from_trasnducers_fitting["se", "range_resource"])

all_pars = NULL
correlations = NULL

concentrations_transducers = c(0, 1, 10, 20, 100, 200, 500, 1000)
weighted = FALSE

for (i in 1:run_number) {
  starting_pars = c("range_CocE" = range_CocE_starting_pars[i]
                    )
  fixed_pars = c("hill_transfer" = hill_transfer_starting_pars[i], 
                 "Km" = Km_starting_pars[i],  
                 "fold_change" = fold_change_starting_pars[i],
                 "baseline" = baseline_starting_pars[i],  
                 "total_enzyme" = total_enzyme_starting_pars[i],
                 "ratio_hip_benz" = ratio_hip_benz_starting_pars[i],
                 "cooperativity_resource" = cooperativity_resource_starting_pars[i],
                 "range_resource" = range_resource_starting_pars[i]
  )
  # print(starting_pars)
  list_of_arguments_of_F <- list(
    "fixed_pars" = fixed_pars,
    "means_data_coc" = means_data_coc, 
    "sd_data_coc" = sd_data_coc, 
    "variable_name_coc" = 'coc', 
    "objective_function_name" = "LSE"
  )
  
  exp_optim <- optim(par = starting_pars, hessian =TRUE, fn = F_to_optimise_cocaine, 
                     method = "L-BFGS-B", lower = coce_lower_bound, upper = coce_upper_bound, 
                     list_of_arguments_of_F = list_of_arguments_of_F)
  
  pars_optimised <- exp_optim$par
  pars_optimised = c(pars_optimised, fixed_pars)
  all_pars = rbind(all_pars, pars_optimised)

}

all_pars_mean = data.frame(t(apply(all_pars,2,mean)))  
all_pars_sd= data.frame(t(apply(all_pars,2,sd))) 
all_pars_sem = all_pars_sd/sqrt(run_number)
all_pars_95_CI = 1.96 * all_pars_sd/sqrt(run_number)
parameters_mean_and_ci = rbind(all_pars_mean, all_pars_sd, all_pars_sem, all_pars_95_CI)
write.csv(x = parameters_mean_and_ci, paste(folder_for_images, "mean_and_ci.csv", sep = '/'), row.names = FALSE)

# Simulating and visualising the hippurate transducer

random_fit_cocaine_transducer = data.frame(cbind(concentrations_transducers, apply(all_pars,1, calculate_coc_transducer, concentrations_coc = concentrations_transducers)))

unmelted_random_fit = random_fit_cocaine_transducer
colnames(random_fit_cocaine_transducer) = c("concentrations", 1:run_number)
random_fit_cocaine_transducer = melt(data= random_fit_cocaine_transducer, id.vars = "concentrations")
df_observed = melt(means_data_coc, id.vars = "concentrations")
sd_observed = melt(sd_data_coc, id.vars = "concentrations")

plotting_random_pars_and_experiments_local(df_plot_total = random_fit_cocaine_transducer, df_observed = df_observed, sd_observed = sd_observed, repeats = 3, 
                                           legend_title = "Constructs", 
                                           x_axis = "concentrations", y_axis = 'value', x_lab = "Cocaine concentration (µM)" ,y_lab ="GFP by OD (AU)", 
                                           title = "Cocaine transducer", subtitle = "Comparing data and 100 best fits", 
                                           xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random")


saving_graph(filename = paste("cocaine_simulation_fit.jpeg", sep = '_'), path =folder_for_images)
cocaine_mean_fit = data.frame(cbind(concentrations_transducers, means_random_fit = apply(unmelted_random_fit[,2:(run_number +1)],1, mean)))
colnames(cocaine_mean_fit) = c("Concentrations", "Fluorescence level")
write.csv(cocaine_mean_fit, "cocaine_mean_simulation_fit.csv", row.names = FALSE)

