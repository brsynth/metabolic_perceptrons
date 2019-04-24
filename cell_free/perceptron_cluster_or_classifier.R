# Author: Mathilde Koch, INRA, date: 06/02/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: Predict the full or classifier with clustering of inputs- imperfect signal

#↓cleaning
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
          "hippurate_high" = 100,
          "benzamid_high" = 100,
          "biphenyl_high" = 100)

# Experimental and best prediction

weights = c("cocaine_weight" = 1,
            "hippurate_weight" = 1,
            "benzamid_weight" = 1,
            "biphenyl_weight" = 10)

for (max_cluster in c(20)) {
  for (min_cluster in c(2)) {
    
    all_classifiers_t = read.csv(file = "implemented_classifiers_transposed.csv")
    
    # names(conditions) = c("coc", "hip", "benzamid", "biphenyl")
    condition_bis = runif(10, min = 0, max = 0.05)
    condition_high= runif(10, min = 0.95, max = 1)
    # Calculations for the logic gate predictions
    expected_behavior = all_classifiers_t[which(all_classifiers_t[,c("input")] == "or"),]

    new_results_cluster = data.frame(all_classifiers_t[which(all_classifiers_t[,c("input")] != "complex"),])
    rownames(new_results_cluster) = NULL
    n_simulations = 50
    for (i in 1:n_simulations) {
      resulting_row = c()
      for (gate_name in colnames(all_classifiers_t)[2:length(colnames(all_classifiers_t))]) {
        conditions = all_classifiers_t[, c(gate_name)]
        names(conditions)= c("coc", "hip", "benzamid", "biphenyl","or", "complex")
        random_change = runif(4, min = 0, max = 1)
        altered_conditions = c()
        for (i in 1:4) {
          # Changing the high and low values taking clustering into account
          if (conditions[i] == 0) {
            new_condition = random_change[i] * min_cluster * 0.01
          } else {
            new_condition = 1 - random_change[i] * max_cluster * 0.01
          }
          altered_conditions = c(altered_conditions, new_condition)
        }
        names(altered_conditions) = c("coc", "hip", "benzamid", "biphenyl")
        result = calculating_point_value_cluster(conditions = altered_conditions, 
                                                 chosen_parameter_set = chosen_parameter_set, 
                                                 highs = highs, weights = weights)
        resulting_row = c(resulting_row, result)
      }
      resulting_row = c("input" = i, resulting_row)
      names(resulting_row) = colnames(new_results_cluster)
      new_results_cluster = rbind(new_results_cluster, resulting_row)
    }
    
    calculating_threshold_dataframe = new_results_cluster[5:nrow(new_results_cluster),]
    max_and_min = c()
    for (i in 2:ncol(calculating_threshold_dataframe)) {
      if (calculating_threshold_dataframe[1,i] == 0) {
        threshold = max(as.numeric(calculating_threshold_dataframe[2:nrow(calculating_threshold_dataframe),i]))
        max_and_min = c(max_and_min, threshold)
      } else {
        threshold = min(as.numeric(calculating_threshold_dataframe[2:nrow(calculating_threshold_dataframe),i]))
        max_and_min = c(max_and_min, threshold)
      }
    }
    calculating_threshold_dataframe = data.frame(t(data.frame(rbind(calculating_threshold_dataframe[1,2:ncol(calculating_threshold_dataframe)], max_and_min))))
    colnames(calculating_threshold_dataframe) = c("Expected", "Value")
    
    min_threshold = max(calculating_threshold_dataframe[which(calculating_threshold_dataframe[,c("Expected")] == 0), c('Value')])
    max_threshold = min(calculating_threshold_dataframe[which(calculating_threshold_dataframe[,c("Expected")] == 1), c('Value')])
    
    rownames(new_results_cluster) = NULL
    visualisation_dataframe = new_results_cluster[6:nrow(new_results_cluster),]
    visualisation_dataframe = rbind(colnames(visualisation_dataframe), visualisation_dataframe)
    rownames(visualisation_dataframe) = NULL
    colnames(visualisation_dataframe) = NULL
    df_for_visu = data.frame(t(visualisation_dataframe))
    rownames(df_for_visu) = NULL
    df_for_visu= df_for_visu[2:nrow(df_for_visu),]
    
    colnames(df_for_visu) = c("combi", 2:ncol(df_for_visu))
    colours= unlist(as.vector(expected_behavior)[2:length(as.vector(expected_behavior))])
    colours = rep(x = colours, times = n_simulations)
    # colours = factor(colours)
    
    
    df_for_visu[,c("combi")] <- as.factor(df_for_visu[,c("combi")])
    df_for_visu = melt(df_for_visu, id.vars = "combi",  measure.vars = 2:ncol(df_for_visu))
    df_for_visu[,c("value")] = as.numeric(df_for_visu[,c("value")])
    # jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    # geom_line(data = df_observed, aes(x = df_observed[,x_axis], y = df_observed[,variable_of_interest]), color = "black",lwd = 1 )
    digital = factor(colours)
    p <- ggplot(data = df_for_visu, aes(x= df_for_visu[,c("combi")], y= df_for_visu[,c("value")], colour = digital), group = "variable")
    # p <- ggplot(data = df_for_visu, aes(x= df_for_visu[,c("combi")], y= df_for_visu[,c("value")], colour = factor(colours)), group = "variable")
    p <- p + geom_point()
    p <- p + geom_hline(yintercept =min_threshold, colour = 2)
    p <- p + geom_hline(yintercept =max_threshold, colour = 5)
    # p <- p + scale_colour_manual(labels = c('YEs', "no"), aesthetics = "colour")
    p <- p + labs(x = "Combination", y = "Fluorescence level", title = "title", subtitle = "subtitle")
    p <- p + theme_bw() #theme black and white
    print(p)
    
    name= paste("cluster_or_gate_hipo_", weights["hippurate_weight"], '_max_', (100-max_cluster), '_min_', min_cluster, sep = '_')
    saving_graph(filename = paste(name, '.jpeg', sep = ''),path =  "clusters/")
    file_name = paste("clusters/", name, ".csv", sep = '')
    write.table(t(new_results_cluster), file_name, eol = "\n \n", sep = ";", row.names = TRUE, col.names = FALSE)
  }
}
