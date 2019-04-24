# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: the aim of this file is to store all functions necessary for generating and analysing data from the Metabolic Perceptrons for Neural Computing in Biological Systems

# Models

# Hill function characterising the biosensor
transfer_function <-function(total, hill_transfer, Km, fold_change, baseline) {
  fraction = total^hill_transfer/(total^hill_transfer + Km^hill_transfer) * baseline
  return(fraction * fold_change + baseline)
}

# Function characterising the transducer
transducer= function(inducer, range_enzyme) {
  total = range_enzyme * inducer
  return(total)
}

# Function to calculate the biosensor results on desired concentrations
calculate_biosensor = function(parameters, concentrations = 0) {
  calculated_biosensor = NULL
  
  for (benzoic_acid in concentrations) {
    value_here = transfer_function(total = benzoic_acid, hill_transfer = parameters["hill_transfer"], 
                                   Km = parameters["Km"], fold_change = parameters["fold_change"],
                                   baseline = parameters["baseline"])
    calculated_biosensor = rbind(calculated_biosensor, value_here)
  }
  calculated_biosensor = data.frame(calculated_biosensor)
  colnames(calculated_biosensor) = 'biosensor'
  return(calculated_biosensor[,c("biosensor")])
}

# Function to calculate the adder results on desired concentrations
calculate_hip_benz_adder = function(parameters, concentrations_benzal = 0, concentrations_hip = 0) {
  # Plasmid copy number
  benz_enz = 5
  hipo = 5
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Resource parameters
  range_resource = parameters["range_resource"] 
  total_enzyme = parameters["total_enzyme"] 
  cooperativity_resource = parameters["cooperativity_resource"] 
  ratio_hip_benz = parameters["ratio_hip_benz"] 
  
  adder_df = NULL
  for (hip in concentrations_hip) {
    for (benzal in concentrations_benzal) {
      benzal_total = transducer(inducer = benzal, range_enzyme = parameters["range_BenZ"])
      hip_total = transducer(inducer = hip, range_enzyme = parameters["range_HipO"])
      total = benzal_total + hip_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      value_here = range_resource * value_here * total_enzyme^cooperativity_resource/((benz_enz + ratio_hip_benz * hipo)^cooperativity_resource + total_enzyme^cooperativity_resource)
      adder_df = rbind(adder_df, c("hip" = hip, "benzal" = benzal, 'GFP' = value_here))
    }
  }
  adder_df = data.frame(adder_df)
  colnames(adder_df) = c("hip", "benzal", 'adder')
  return(adder_df[,c("adder")])
}

# Function to calculate the benzaldehyde transducer results on desired concentrations
calculate_benz_transducer = function(parameters, concentrations_benzal = 0) {
  # Plasmid copy number
  benz_enz = 5
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Resource parameters
  range_resource = parameters["range_resource"] 
  total_enzyme = parameters["total_enzyme"] 
  cooperativity_resource = parameters["cooperativity_resource"] 
  ratio_hip_benz = parameters["ratio_hip_benz"] 
  
  transducer_df = NULL
  for (benzal in concentrations_benzal) {
    benzal_total = transducer(inducer = benzal, range_enzyme = parameters["range_BenZ"])
    total = benzal_total
    value_here = transfer_function(total,
                                   hill_transfer = hill_transfer,
                                   Km = Km,
                                   fold_change = fold_change, 
                                   baseline = baseline)
    # Competition
    value_here = range_resource * value_here * total_enzyme^cooperativity_resource/((benz_enz + ratio_hip_benz * 0)^cooperativity_resource + total_enzyme^cooperativity_resource)
    transducer_df = rbind(transducer_df, c("benzal" = benzal, 'GFP' = value_here))
  }
  
  transducer_df = data.frame(transducer_df)
  colnames(transducer_df) = c("benzal", 'transducer')
  return(transducer_df[,c("transducer")])
}

# Function to calculate the benzaldehyde transducer results on desired concentrations
calculate_hip_transducer = function(parameters, concentrations_hip = 0) {
  # Plasmid copy number
  hipo = 5
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Resource parameters
  range_resource = parameters["range_resource"] 
  total_enzyme = parameters["total_enzyme"] 
  cooperativity_resource = parameters["cooperativity_resource"] 
  ratio_hip_benz = parameters["ratio_hip_benz"] 
  
  transducer_df = NULL
  for (hip in concentrations_hip) {
    hip_total = transducer(inducer = hip, range_enzyme = parameters["range_HipO"])
    total = hip_total
    value_here = transfer_function(total,
                                   hill_transfer = hill_transfer,
                                   Km = Km,
                                   fold_change = fold_change, 
                                   baseline = baseline)
    # Competition
    value_here = range_resource * value_here * total_enzyme^cooperativity_resource/((0 + ratio_hip_benz * hipo)^cooperativity_resource + total_enzyme^cooperativity_resource)
    transducer_df = rbind(transducer_df, c("hip" = hip, 'GFP' = value_here))
  }
  transducer_df = data.frame(transducer_df)
  colnames(transducer_df) = c("hip", 'transducer')
  return(transducer_df[,c("transducer")])
}

# Function to calculate the cocaine transducer results on desired concentrations
calculate_coc_transducer = function(parameters, concentrations_coc = 0) {
  # Plasmid copy number
  coce = 5
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Resource parameters
  range_resource = parameters["range_resource"] 
  total_enzyme = parameters["total_enzyme"] 
  cooperativity_resource = parameters["cooperativity_resource"] 
  ratio_hip_benz = parameters["ratio_hip_benz"] 
  
  transducer_df = NULL
  for (coc in concentrations_coc) {
    coc_total = transducer(inducer = coc, range_enzyme = parameters["range_CocE"])
    total = coc_total
    value_here = transfer_function(total,
                                   hill_transfer = hill_transfer,
                                   Km = Km,
                                   fold_change = fold_change, 
                                   baseline = baseline)
    # Competition
    value_here = range_resource * value_here * total_enzyme^cooperativity_resource/((coce)^cooperativity_resource + total_enzyme^cooperativity_resource)
    transducer_df = rbind(transducer_df, c("coc" = coc, 'GFP' = value_here))
  }
  transducer_df = data.frame(transducer_df)
  colnames(transducer_df) = c("coc", 'transducer')
  return(transducer_df[,c("transducer")])
}


# Objective functions

# Weighted least squared errors
weighted_LSE <- function(prediction, observed, sd_observed, variable_name = "coce_coc") {
  colnames(prediction) = colnames(observed)
  sum <- 0
  for (i in 1:nrow(observed)) {
    errort <- 0
    obs <- observed[i, variable_name]
    pred <- prediction[i, variable_name]
    sd <- sd_observed[i, variable_name]
    errort <- errort + (obs - pred)*(obs - pred)/(sd *sd)
    sum <- c(sum, unname(errort))
  }
  
  if (class(sum) == "list") {
    sum = unlist(sum, recursive = TRUE, use.names = FALSE)
  }
  return(sum(sum))
}

# Least squared errors
LSE <- function(prediction, observed, variable_name = "coce_coc", sd_observed = NULL) {
  colnames(prediction) = colnames(observed)
  #  I need to know the variance of my data, I divide each error by the variance of the data. Check that when I understand my data
  sum <- 0
  for (i in 1:nrow(observed)) {
    errort <- 0
    obs <- observed[i, variable_name]
    pred <- prediction[i, variable_name]
    errort <- errort + (obs - pred)*(obs - pred)
    sum <- c(sum, unname(errort))
  }
  
  if (class(sum) == "list") {
    sum = unlist(sum, recursive = TRUE, use.names = FALSE)
  }
  return(sum(sum))
}

# Scoring functions

r_squared = function(y_true, y_pred, y_std = NULL) {
  if (is.null(y_std)) {
    weights = rep(1, length(y_true))
  } else {
    weights = 1/(y_std * y_std)
  }
  y_true_mean = mean(y_true)
  numerator = 0
  denominator = 0
  for (i in 1:length(y_true)) {
    numerator = numerator + weights[i] * (y_true[i] - y_pred[i])  * (y_true[i] - y_pred[i]) 
    denominator = denominator + weights[i] * (y_true[i] - y_true_mean)  * (y_true[i] - y_true_mean) 
  }
  return(1 - numerator/denominator)
}

RMSD = function(y_true, y_pred) {
  sum = 0
  for (i in 1:length(y_pred)) {
    sum = sum + (y_true[i] - y_pred[i])  * (y_true[i] - y_pred[i])
  }
  RMSD = sqrt(sum/length(y_pred))
  return(RMSD)
}

calculate_all_score_pred_true_true_sd = function(means_data, sd_data) {
  RMSD = RMSD(means_data[,c("Data")], means_data[,c("Model")])
  correlation = cor(means_data[,c("Data")], means_data[,c("Model")])
  r_squared = r_squared(y_true = means_data[,c("Data")], 
                        y_pred = means_data[,c("Model")], 
                        y_std = sd_data[,c("Data")])
  r_squared_unweighted = r_squared(y_true = means_data[,c("Data")], 
                                   y_pred = means_data[,c("Model")])
  return(c("RMSD" = RMSD, "correlation"= correlation, "r_squared" = r_squared,
           "r_squared_unweighted" = r_squared_unweighted))
}

# Graph visualisation and saving

# Wrapper for 2D heatmap visualisation
plot_2D_heatmap = function(df_means, df_sd = NULL, title_base = "Visualising adder in 2D heatmap", name_for_saving, folder_for_concentration_images_strategy, saving = TRUE, x_lab = "(nM)", y_lab = "(µM)", subtitle_here = "GFP") {
  experiment_1 = colnames(df_means)[1]
  experiment_2 = colnames(df_means)[2]
  experiments = c(experiment_1, experiment_2)
  
  GFP_by_OD_factor <- melt(df_means, id.vars = experiments)
  GFP_by_OD_factor[,c(experiment_1)] <- as.factor(GFP_by_OD_factor[,c(experiment_1)])
  GFP_by_OD_factor[,c(experiment_2)] <- as.factor(GFP_by_OD_factor[,c(experiment_2)])
  
  xlab = paste(experiment_1, x_lab, sep = ' ')
  ylab = paste(experiment_2, y_lab, sep = ' ')
  
  print_heatmap(data_as_factor = GFP_by_OD_factor, x_lab = xlab, y_lab = ylab, x_axis = experiment_1, y_axis = experiment_2, title = title_base, subtitle = subtitle_here)
  if (saving) {
    folder_to_save = folder_for_concentration_images_strategy
    saving_graph(filename = paste("heat_map", 'GFP_by_OD',  name_for_saving, sep = '_'), path = folder_to_save)
  }
}

# Plotting the 2D heatmap
print_heatmap <- function(data_as_factor, x_axis = "coc", y_axis = "hip", x_lab = "Compound 1", y_lab ="Compound 2", title = "", subtitle = "", 
                          print = TRUE) {
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  p <- ggplot(data = data_as_factor, aes(x= data_as_factor[,x_axis], y= data_as_factor[,y_axis], fill=value))
  p <- p + geom_tile()
  p <- p + scale_fill_gradientn(colours = jet.colors(7))
  p <- p + labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  p <- p + theme_bw() #theme black and white
  if (print) {
    print(p)
  } else {
    return(p)
  }
}

# Plotting biosensor or transducers model and data on the same graph
print_same_graph_error_bars_biosensor <- function(df, df_errors = NULL, df_observed= NULL, df_sd_observed = NULL, variable_of_interest = NULL, repeats = 3, error_bars = TRUE, 
                                                  x_axis = "time",y_axis = 'value', x_lab = "Time (s)" ,y_lab ="Normalised signal", title = "", subtitle = "Subtitle", 
                                                  names_titles = NULL, legends = NULL, xlog =TRUE, print = TRUE, name_condition = "Condition", legend = "bottom", label_names = NULL, log_ticks = FALSE) {
  p <- ggplot(data = df) #plotting the dataframe
  p <- p + geom_line(aes(x = df[,x_axis], y = df[,y_axis], group= variable, colour = variable)) #plotting value againt inducer ratio
  if (!is.null(df_errors)){
    if (error_bars) {
      limits <- aes(x = df[,x_axis], ymax = df[,y_axis] + df_errors[,y_axis]/sqrt(repeats), ymin = df[,y_axis] - df_errors[,y_axis]/sqrt(repeats), group= variable, colour = variable)
      p <- p + geom_errorbar(limits, width = 0.2)
    }
  }
  
  if (!is.null(df_observed)) {
    p <- p + geom_line(data = df_observed, aes(x = df_observed[,x_axis], y = df_observed[,variable_of_interest]), color = "black",lwd = 1 )
    if (!is.null(df_sd_observed)){
      if (error_bars) {
        limits <- aes(x = df_observed[,x_axis], ymax = df_observed[,variable_of_interest] + df_sd_observed[,variable_of_interest]/sqrt(repeats), 
                      ymin = df_observed[,variable_of_interest] - df_sd_observed[,variable_of_interest]/sqrt(repeats))
        p <- p + geom_errorbar(limits, color = "black",lwd = 1, linetype = "dotted")
      } 
    }
  }
  if (xlog == TRUE) {
    p <- p + scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
    if (log_ticks) {    
      p <- p + annotation_logticks(sides = "b",scaled = TRUE) #log ticks
    } 
    p <- p + scale_y_continuous()
  } else {
    p <- p + scale_x_continuous()
    p <- p + scale_y_continuous()
  }
  p <- p + labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  #legend
  if (is.null(label_names)) {
    label_names = names_titles
  }
  p <- p + scale_colour_hue(name = name_condition,
                            breaks = names_titles,
                            labels = label_names)
  
  p <- p + theme_bw() #theme black and white
  p <- p + theme(legend.position=legend)
  if (print == TRUE) {
    print(p)
  } else {
    return(p)
  }
}

# Plotting biosensor or transducers data and run_number simulations
plotting_random_pars_and_experiments_local <- function(df_plot_total, df_observed, sd_observed, repeats = 1, legend_title = "Constructs", 
                                                       x_axis = "concentrations", y_axis = 'value', x_lab = "Concentrations (µM)" ,y_lab ="RFP by OD (AU)", title = "Comparing best fit and experimental data", subtitle = "Substitle", 
                                                       xlog =TRUE, legend = "bottom", folder_to_save = NULL, title_to_save = "random") {
  p <- ggplot(data = df_observed) #plotting the dataframe
  p <- p + geom_line(aes(x = df_observed[,x_axis], y = df_observed[,y_axis], group= 'variable'), colour = 'black') #plotting value againt inducer ratio
  limits <- aes(x = df_observed[,x_axis], ymax = df_observed[,y_axis] + sd_observed[,y_axis]/sqrt(repeats), 
                ymin = df_observed[,y_axis] - sd_observed[,y_axis]/sqrt(repeats))
  p <- p + geom_errorbar(limits, color = "black",lwd = 0.5)
  #  p <- p + geom_line(data = df_plot_total, aes(x = df_plot_total[,x_axis], y = df_plot_total[,y_axis], group = variable, alpha = 0.01, colour = "black", linetype = 'dotted'))
  
  p <- p + geom_line(data = df_plot_total, aes(x = df_plot_total[,x_axis], y = df_plot_total[,y_axis], group = variable, alpha = 0.01), colour = "blue")
  if (xlog == TRUE) {
    p <- p + scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
    p <- p + annotation_logticks(sides = "b",scaled = TRUE) #log ticks
    p <- p + scale_y_continuous()
  } else {
    p <- p + scale_x_continuous()
    p <- p + scale_y_continuous()
  }
  
  p <- p + guides(alpha=FALSE, linetype = FALSE)
  p <- p + scale_colour_discrete(name = legend_title,
                                 breaks = colours,
                                 labels = c("Data", "Experiment"))
  
  p <- p + labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  p <- p + theme_bw() #theme black and white
  p <- p + theme(legend.position=legend)
  
  p
  if (!is.null(folder_to_save)) {
    saving_graph(filename = title_to_save, path =folder_to_save)
  }
  print(p)
}

# Saves the laast graph to disk
saving_graph <- function(plot = last_plot(),filename, path = getwd(), device = "jpeg", 
                         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                         dpi = 300, limitsize = TRUE) {
  
  ggsave(filename = filename, plot = plot, device = device, path = path, scale = scale, width = width, height = height, units = units,
         dpi = dpi, limitsize = limitsize)
}

