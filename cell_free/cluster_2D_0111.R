# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: Predict the 2D cluster for Figure 2D

# cleaning
rm(list = ls())

#importing packages: 

library(rootSolve)
library(deSolve)
library(ggplot2)
library('reshape2')

#importing source

source("helper_function_synthetic_metabolic_circuits_cell_free.R")

parameters_from_fitting = read.csv("t_linear_full_range_full_run_1/mean_and_ci.csv")
rownames(parameters_from_fitting) = c("mean", "sd", "se","95_ci")
chosen_parameter_set = unlist(parameters_from_fitting[c("mean"),])

# Concentration for ON inputs (inducer concentrations)

highs = c("cocaine_high" = 100,
          "hippurate_high" = 100)

hill_transfer = parameters_from_fitting["mean","hill_transfer"]
Km = parameters_from_fitting["mean","Km"]
fold_change = parameters_from_fitting["mean","fold_change"]
baseline = parameters_from_fitting["mean","baseline"]
slower_slope = parameters_from_fitting["mean","slower_slope"]


for (coc_weight in c(0.1)) {
  for (hip_weight in c(0.3, 1)) {
    weights = c("cocaine_weight" = coc_weight,
                "hippurate_weight" = hip_weight)
    
    n_simulations = 50
    for (max_cluster in c(20)) {
      for (min_cluster in c(2)) {
        # Calculations for the logic gate predictions
        # x_low = runif(n_simulations, min = 0, max = 0.01  * min_cluster)
        x_low = runif(n_simulations, min = 1 - 0.01  * max_cluster, max = 1)
        y_low = runif(n_simulations, min = 0, max = 0.01  * min_cluster)
        
        x_high = runif(n_simulations, min = 1 - 0.01  * max_cluster, max = 1)
        y_high = runif(n_simulations, min = 1 - 0.01  * max_cluster, max = 1)
        
        expected_behavior = c(rep(0, times = n_simulations), rep(1, times = n_simulations))
        x = c(x_low, x_high)
        y = c(y_low, y_high)
        
        all_conditions = cbind(x, y)
        colnames(all_conditions) = c("coc", 'hip')
        
        all_results = c()
        for (i in seq(1, nrow(all_conditions))) {
          conditions = all_conditions[i, ]
          result = calculating_point_value_2D_cluster(conditions = conditions, 
                                                   chosen_parameter_set = chosen_parameter_set, 
                                                   highs = highs, weights = weights)
          all_results = c(all_results, result)
        }
        all_conditions = cbind(all_conditions, all_results)
        rownames(all_conditions) = NULL
        all_conditions = data.frame(all_conditions)

        calculating_threshold_dataframe = data.frame(cbind(all_conditions, expected_behavior))
        min_threshold = max(calculating_threshold_dataframe[which(calculating_threshold_dataframe[,c("expected_behavior")] == 0), c('all_results')])
        max_threshold = min(calculating_threshold_dataframe[which(calculating_threshold_dataframe[,c("expected_behavior")] == 1), c('all_results')])
        
        # boundary_min = reciprocate_transfer(min_threshold, hill_transfer, Km, fold_change, baseline, slower_slope)
        # result_min = data.frame(line_from_boundary(boundary_min))
        # 
        # boundary_max = reciprocate_transfer(max_threshold, hill_transfer, Km, fold_change, baseline, slower_slope)
        # result_max = data.frame(line_from_boundary(boundary_max))
        
        all_thresholds = data.frame()
        for (threshold in c(2.5, 3, 3.5)) {
          boundary_threshold = reciprocate_transfer(threshold, hill_transfer, Km, fold_change, baseline, slower_slope)
          result_threshold = data.frame(line_from_boundary_horizontal(boundary_threshold))
          all_thresholds = rbind(all_thresholds, cbind(result_threshold, rep(threshold, times = nrow(result_threshold))))
        }
        colnames(all_thresholds) = c("x", "y", "threshold")
        RFU = all_conditions[,c("all_results")]
        shapes = expected_behavior
        p <- ggplot(data = all_conditions, aes(x= all_conditions[,c("coc")], y= all_conditions[,c("hip")], colour = RFU))
        p <- p + geom_point()
        p <- p + geom_line(data = all_thresholds, aes(x= all_thresholds[,c("x")], y= all_thresholds[,c("y")],
                                                      group = all_thresholds[,c("threshold")], colour = all_thresholds[,c("threshold")]), size = 2)
        p <- p + scale_color_gradient2(mid="blue", high="red")
        p <- p + labs(x = "Cocaine", y = "Hippurate", title = "title", subtitle = "subtitle")
        p <- p + theme_bw() #theme black and white
        print(p)
        
        name= paste("2D_cluster_01_11_HipO", weights["hippurate_weight"], "CocE",weights["cocaine_weight"],  '_max_', (100-max_cluster), '_min_', min_cluster, sep = '_')
        saving_graph(filename = paste(name, '.jpeg', sep = ''),path =  "images/")
      }
    }
  }
}



