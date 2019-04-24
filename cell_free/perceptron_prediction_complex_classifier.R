# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: Predict the full complex classifier for best prediction and working experimental implementation

#↓cleaning
rm(list = ls())

#importing packages: 

library(rootSolve)
library(deSolve)
library(ggplot2)
library('reshape2')

#importing source

source("helper_function_synthetic_metabolic_circuits_cell_free.R")

all_classifiers = read.csv(file = "implemented_classifiers.csv")
parameters_from_fitting = read.csv("t_linear_full_range_full_run_1/mean_and_ci.csv")

rownames(parameters_from_fitting) = c("mean", "sd", "se","95_ci")
chosen_parameter_set = unlist(parameters_from_fitting[c("mean"),])

# Concentration for ON inputs (inducer concentrations)

highs = c("cocaine_high" = 100,
          "hippurate_high" = 100,
          "benzamid_high" = 100,
          "biphenyl_high" = 100)

# EXPERIMENTAL IMPLEMENTATION

weights = c("cocaine_weight" = 0.1,
            "hippurate_weight" = 0.03,
            "benzamid_weight" = 1,
            "biphenyl_weight" = 10)

results_dataframe = all_classifiers[c("coc", "hip", "benzamid", "biphenyl")]

names_of_combinations = c()
calculated_gate = c()
for (element in 1:dim(all_classifiers)[1]) {
  conditions = all_classifiers[element, c("coc", "hip", "benzamid", "biphenyl")]
  condition_result = calculating_point_value(conditions, chosen_parameter_set, highs, weights) 
  
  calculated_gate = c(calculated_gate, unname(unlist(condition_result)))
  names_of_combinations = c(names_of_combinations, paste(conditions, sep = "", collapse = "_"))
}

results_dataframe = cbind(results_dataframe, "names" = names_of_combinations , "expected" = all_classifiers[,c("coc_and_hip_or_biphenyl_or_benzamid")], "value" = calculated_gate)

testing_threshold = c(1, 1.2, 1.5,2, 2.5,2.7, 3, 3.1, 3.3, 3.5, 3.6, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.5)
colnames = colnames(results_dataframe)

for (threshold in testing_threshold) {
  colnames = c(colnames, threshold)
  results_dataframe = cbind(results_dataframe, calculation = as.numeric(results_dataframe[,c("value")] > threshold))
}
results_dataframe = data.frame(results_dataframe)
colnames(results_dataframe) = colnames
write.table(results_dataframe, paste("complex_", paste(weights, collapse = "_"), ".csv", sep = ''), sep = ";")
plot_logic_gates(results_dataframe, testing_threshold = testing_threshold, 
                 name_for_saving = paste("complex_", paste(weights, collapse = "_"), sep = ''),
                 title = "Predicting complex gate")


# Best prediction

weights = c("cocaine_weight" = 0.1,
            "hippurate_weight" = 0.1,
            "benzamid_weight" = 1,
            "biphenyl_weight" = 10)

results_dataframe = all_classifiers[c("coc", "hip", "benzamid", "biphenyl")]

names_of_combinations = c()
calculated_gate = c()
for (element in 1:dim(all_classifiers)[1]) {
  conditions = all_classifiers[element, c("coc", "hip", "benzamid", "biphenyl")]
  condition_result = calculating_point_value(conditions, chosen_parameter_set, highs, weights) 
  
  calculated_gate = c(calculated_gate, unname(unlist(condition_result)))
  names_of_combinations = c(names_of_combinations, paste(conditions, sep = "", collapse = "_"))
}

results_dataframe = cbind(results_dataframe, "names" = names_of_combinations , "expected" = all_classifiers[,c("coc_and_hip_or_biphenyl_or_benzamid")], "value" = calculated_gate)

testing_threshold = c(1, 1.2, 1.5,2, 2.5,2.7, 3, 3.1, 3.3, 3.5, 3.6, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.5)
colnames = colnames(results_dataframe)

for (threshold in testing_threshold) {
  colnames = c(colnames, threshold)
  results_dataframe = cbind(results_dataframe, calculation = as.numeric(results_dataframe[,c("value")] > threshold))
}
results_dataframe = data.frame(results_dataframe)
colnames(results_dataframe) = colnames
write.table(results_dataframe, paste("complex_", paste(weights, collapse = "_"), ".csv", sep = ''), sep = ";")
plot_logic_gates(results_dataframe, testing_threshold = testing_threshold, 
                 name_for_saving = paste("complex_", paste(weights, collapse = "_"), sep = ''),
                 title = "Predicting complex gate")
