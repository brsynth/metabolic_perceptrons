# Author: Mathilde Koch, date: 10/01/19
# Article: Metabolic Perceptrons for Neural Computing in Biological Systems

# Aim: the aim of this file is to store all functions necessary for 
# generating and analysing data from the synthetic metabolic circuits article


# Models

# Hill function characterising the biosensor
transfer_function <-function(total, hill_transfer, Km, fold_change, baseline) {
  fraction = total^hill_transfer/(total^hill_transfer + Km^hill_transfer) * baseline
  return(fraction * fold_change + baseline)
}


transfer_function_linear <-function(total, hill_transfer, Km, fold_change, baseline,
                                    slower_slope) {
  fraction = total^hill_transfer/(total^hill_transfer + Km^hill_transfer) * baseline
  return(fraction * fold_change + baseline + slower_slope * 0.0001* total)
}

# Calculating benzoic acid for a fluoresnce value of fluo
reciprocate_transfer <-function(fluo, hill_transfer, Km, fold_change, baseline, slower_slope) {
  best_of_x = 0
  best_of_fluo = 10000
  for (x in seq(0, 1000)) {
    calculated_fluo = transfer_function_linear(x, hill_transfer, Km, fold_change, baseline,
                                               slower_slope)
    if (abs(calculated_fluo - fluo) <= abs(best_of_fluo - fluo)) {
      best_of_fluo = calculated_fluo
      best_of_x = x
    }
    
  }
  return(best_of_x)
}
# Function characterising the transducer: to modify for cell-free, much more complex
# transducer= function(inducer, range_enzyme) {
#   total = range_enzyme * inducer
#   return(total)
# }

transducer <- function(enzyme, inducer, range_enzyme, enzyme_constant, inducer_constant, hill_enzyme, hill_inducer) {
  total = range_enzyme * enzyme^hill_enzyme/(enzyme^hill_enzyme + enzyme_constant^hill_enzyme) * inducer^hill_inducer/(inducer^hill_inducer + inducer_constant^hill_inducer)
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

# Function to calculate the linear biosensor results on desired concentrations
calculate_biosensor_linear = function(parameters, concentrations = 0) {
  calculated_biosensor = NULL
  
  for (benzoic_acid in concentrations) {
    value_here = transfer_function_linear(total = benzoic_acid, hill_transfer = parameters["hill_transfer"], 
                                   Km = parameters["Km"], fold_change = parameters["fold_change"],
                                   baseline = parameters["baseline"], slower_slope = parameters["slower_slope"])
    calculated_biosensor = rbind(calculated_biosensor, value_here)
  }
  calculated_biosensor = data.frame(calculated_biosensor)
  colnames(calculated_biosensor) = 'biosensor'
  return(calculated_biosensor[,c("biosensor")])
}

# Function to calculate the non-linear weighted enzyme adder on desired concentrations

calculate_hip_coc_enzymatic_adder = function(parameters, concentrations_hipo = 0, concentrations_coce = 0) {
  # It is at the moment written only for hippuric acid and cocaine
  # Inducer concentrations
  hip = 100
  coc = 100
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Enzyme parameters
  range_CocE = parameters["range_CocE"] 
  CocE_constant = parameters["CocE_constant"] 
  cocaine_constant = parameters["cocaine_constant"]
  hill_CocE = parameters["hill_CocE"] 
  hill_cocaine = parameters["hill_cocaine"]
  
  range_HipO = parameters["range_HipO"] 
  HipO_constant = parameters["HipO_constant"] 
  hippurate_constant = parameters["hippurate_constant"]
  hill_HipO = parameters["hill_HipO"] 
  hill_hippurate = parameters["hill_hippurate"]
  
  enzymatic_adder_df = NULL
  for (hipo in concentrations_hipo) {
    for (coce in concentrations_coce) {
      HipO = hipo
      CocE = coce
      
      hip_total = transducer(enzyme = HipO, inducer = hip, range_enzyme = range_HipO, enzyme_constant = HipO_constant, 
                             inducer_constant = hippurate_constant, hill_enzyme = hill_HipO, hill_inducer = hill_hippurate)
      
      coc_total = transducer(enzyme = CocE, inducer = coc, range_enzyme = range_CocE, enzyme_constant = CocE_constant, 
                             inducer_constant = cocaine_constant, hill_enzyme = hill_CocE, hill_inducer = hill_cocaine)
      
      total = hip_total + coc_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      enzymatic_adder_df = rbind(enzymatic_adder_df, c("HipO" = hipo, "CocE" = coce, 'GFP' = value_here))
    }
  }
  enzymatic_adder_df = data.frame(enzymatic_adder_df)
  colnames(enzymatic_adder_df) = c("hipo", "coce", 'coce_hipo')
  return(enzymatic_adder_df[,c("coce_hipo")])
}

calculate_hip_coc_adder = function(parameters, concentrations_hip = 0, concentrations_coc = 0) {
  # It is at the moment written only for hippuric acid and cocaine
  # Inducer concentrations
  hipo = 1
  coce = 3
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Enzyme parameters
  range_CocE = parameters["range_CocE"] 
  CocE_constant = parameters["CocE_constant"] 
  cocaine_constant = parameters["cocaine_constant"]
  hill_CocE = parameters["hill_CocE"] 
  hill_cocaine = parameters["hill_cocaine"]
  
  range_HipO = parameters["range_HipO"] 
  HipO_constant = parameters["HipO_constant"] 
  hippurate_constant = parameters["hippurate_constant"]
  hill_HipO = parameters["hill_HipO"] 
  hill_hippurate = parameters["hill_hippurate"]
  
  adder_df = NULL
  for (hip in concentrations_hip) {
    for (coc in concentrations_coc) {
      HipO = hipo
      CocE = coce
      
      hip_total = transducer(enzyme = HipO, inducer = hip, range_enzyme = range_HipO, enzyme_constant = HipO_constant, 
                             inducer_constant = hippurate_constant, hill_enzyme = hill_HipO, hill_inducer = hill_hippurate)
      
      coc_total = transducer(enzyme = CocE, inducer = coc, range_enzyme = range_CocE, enzyme_constant = CocE_constant, 
                             inducer_constant = cocaine_constant, hill_enzyme = hill_CocE, hill_inducer = hill_cocaine)
      
      total = hip_total + coc_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      adder_df = rbind(adder_df, c("cocaine" = coc,"hippurate" = hip, 'GFP' = value_here))
    }
  }
  adder_df = data.frame(adder_df)
  colnames(adder_df) = c("coc", "hip", 'coc_hip')
  return(adder_df[,c("coc_hip")])
}

calculate_hipo_hip_adder = function(parameters, concentrations_hip = 0, concentrations_HipO = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Enzyme parameters
  range_HipO = parameters["range_HipO"] 
  HipO_constant = parameters["HipO_constant"] 
  hippurate_constant = parameters["hippurate_constant"]
  hill_HipO = parameters["hill_HipO"] 
  hill_hippurate = parameters["hill_hippurate"]
  
  hip_HipO_adder_df = NULL
  for (hipo in concentrations_HipO) {
    for (hip in concentrations_hip) {
      HipO = hipo
      hip_total = transducer(enzyme = HipO, inducer = hip, range_enzyme = range_HipO, enzyme_constant = HipO_constant, 
                             inducer_constant = hippurate_constant, hill_enzyme = hill_HipO, hill_inducer = hill_hippurate)
      
      total = hip_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      hip_HipO_adder_df = rbind(hip_HipO_adder_df, c("hippurate" = hip, "HipO" = hipo, 'GFP' = value_here))
    }
  }
  hip_HipO_adder_df = data.frame(hip_HipO_adder_df)
  colnames(hip_HipO_adder_df) = c("hip", "hipo", 'hipo_hip')
  return(hip_HipO_adder_df[,c("hipo_hip")])
}

calculate_coce_coc_adder = function(parameters, concentrations_coc = 0, concentrations_CocE = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Enzyme parameters
  range_CocE = parameters["range_CocE"] 
  CocE_constant = parameters["CocE_constant"] 
  cocaine_constant = parameters["cocaine_constant"]
  hill_CocE = parameters["hill_CocE"] 
  hill_cocaine = parameters["hill_cocaine"]
  
  coc_CocE_adder_df = NULL
  for (coce in concentrations_CocE) {
    for (coc in concentrations_coc) {
      CocE = coce
      coc_total = transducer(enzyme = CocE, inducer = coc, range_enzyme = range_CocE, enzyme_constant = CocE_constant, 
                             inducer_constant = cocaine_constant, hill_enzyme = hill_CocE, hill_inducer = hill_cocaine)
      
      total = coc_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      coc_CocE_adder_df = rbind(coc_CocE_adder_df, c("cocaine" = coc, "CocE" = coce, 'GFP' = value_here))
    }
  }
  coc_CocE_adder_df = data.frame(coc_CocE_adder_df)
  colnames(coc_CocE_adder_df) = c("coc", "coce", 'coce_coc')
  return(coc_CocE_adder_df[,c("coce_coc")])
}

calculate_benzamid_enz_benzamid_adder = function(parameters, concentrations_benzamid = 0, concentrations_benzamid_enz = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Enzyme parameters
  range_benzamid_enz = parameters["range_benzamid_enz"] 
  benzamid_enz_constant = parameters["benzamid_enz_constant"] 
  benzamid_constant = parameters["benzamid_constant"]
  hill_benzamid_enz = parameters["hill_benzamid_enz"] 
  hill_benzamid = parameters["hill_benzamid"]
  
  benzamid_benzamid_enz_adder_df = NULL
  for (benzamid_enz in concentrations_benzamid_enz) {
    for (benzamid in concentrations_benzamid) {
      benzamid_total = transducer(enzyme = benzamid_enz, inducer = benzamid, range_enzyme = range_benzamid_enz, enzyme_constant = benzamid_enz_constant, 
                             inducer_constant = benzamid_constant, hill_enzyme = hill_benzamid_enz, hill_inducer = hill_benzamid)
      
      total = benzamid_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      benzamid_benzamid_enz_adder_df = rbind(benzamid_benzamid_enz_adder_df, c("benzamid" = benzamid, "benzamid_enz" = benzamid_enz, 'GFP' = value_here))
    }
  }
  benzamid_benzamid_enz_adder_df = data.frame(benzamid_benzamid_enz_adder_df)
  colnames(benzamid_benzamid_enz_adder_df) = c("benzamid", "benzamid_enz", 'benzamid_enz_benzamid')
  return(benzamid_benzamid_enz_adder_df[,c("benzamid_enz_benzamid")])
}

calculate_biphenyl_enz_biphenyl_adder = function(parameters, concentrations_biphenyl = 0, concentrations_biphenyl_enz = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  # Enzyme parameters
  range_biphenyl_enz = parameters["range_biphenyl_enz"] 
  biphenyl_enz_constant = parameters["biphenyl_enz_constant"] 
  biphenyl_constant = parameters["biphenyl_constant"]
  hill_biphenyl_enz = parameters["hill_biphenyl_enz"] 
  hill_biphenyl = parameters["hill_biphenyl"]
  
  biphenyl_biphenyl_enz_adder_df = NULL
  for (biphenyl_enz in concentrations_biphenyl_enz) {
    for (biphenyl in concentrations_biphenyl) {
      biphenyl_total = transducer(enzyme = biphenyl_enz, inducer = biphenyl, range_enzyme = range_biphenyl_enz, enzyme_constant = biphenyl_enz_constant, 
                                  inducer_constant = biphenyl_constant, hill_enzyme = hill_biphenyl_enz, hill_inducer = hill_biphenyl)
      
      total = biphenyl_total
      value_here = transfer_function(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline)
      # Competition
      biphenyl_biphenyl_enz_adder_df = rbind(biphenyl_biphenyl_enz_adder_df, c("biphenyl" = biphenyl, "biphenyl_enz" = biphenyl_enz, 'GFP' = value_here))
    }
  }
  biphenyl_biphenyl_enz_adder_df = data.frame(biphenyl_biphenyl_enz_adder_df)
  colnames(biphenyl_biphenyl_enz_adder_df) = c("biphenyl", "biphenyl_enz", 'biphenyl_enz_biphenyl')
  return(biphenyl_biphenyl_enz_adder_df[,c("biphenyl_enz_biphenyl")])
}


# Function to calculate the linear weighted enzyme adder on desired concentrations

calculate_hip_coc_enzymatic_adder_linear = function(parameters, concentrations_hipo = 0, concentrations_coce = 0) {
  # It is at the moment written only for hippuric acid and cocaine
  # Inducer concentrations
  hip = 100
  coc = 100
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  slower_slope = parameters["slower_slope"]
  # Enzyme parameters
  range_CocE = parameters["range_CocE"] 
  CocE_constant = parameters["CocE_constant"] 
  cocaine_constant = parameters["cocaine_constant"]
  hill_CocE = parameters["hill_CocE"] 
  hill_cocaine = parameters["hill_cocaine"]
  
  range_HipO = parameters["range_HipO"] 
  HipO_constant = parameters["HipO_constant"] 
  hippurate_constant = parameters["hippurate_constant"]
  hill_HipO = parameters["hill_HipO"] 
  hill_hippurate = parameters["hill_hippurate"]
  
  enzymatic_adder_df = NULL
  for (hipo in concentrations_hipo) {
    for (coce in concentrations_coce) {
      HipO = hipo
      CocE = coce
      
      hip_total = transducer(enzyme = HipO, inducer = hip, range_enzyme = range_HipO, enzyme_constant = HipO_constant, 
                             inducer_constant = hippurate_constant, hill_enzyme = hill_HipO, hill_inducer = hill_hippurate)
      
      coc_total = transducer(enzyme = CocE, inducer = coc, range_enzyme = range_CocE, enzyme_constant = CocE_constant, 
                             inducer_constant = cocaine_constant, hill_enzyme = hill_CocE, hill_inducer = hill_cocaine)
      
      total = hip_total + coc_total
      value_here = transfer_function_linear(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline,
                                     slower_slope = slower_slope)
      # Competition
      enzymatic_adder_df = rbind(enzymatic_adder_df, c("HipO" = hipo, "CocE" = coce, 'GFP' = value_here))
    }
  }
  enzymatic_adder_df = data.frame(enzymatic_adder_df)
  colnames(enzymatic_adder_df) = c("hipo", "coce", 'coce_hipo')
  return(enzymatic_adder_df[,c("coce_hipo")])
}

calculate_hip_coc_adder_linear = function(parameters, concentrations_hip = 0, concentrations_coc = 0) {
  # It is at the moment written only for hippuric acid and cocaine
  # Inducer concentrations
  hipo = 1
  coce = 3
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  slower_slope = parameters["slower_slope"]
  # Enzyme parameters
  range_CocE = parameters["range_CocE"] 
  CocE_constant = parameters["CocE_constant"] 
  cocaine_constant = parameters["cocaine_constant"]
  hill_CocE = parameters["hill_CocE"] 
  hill_cocaine = parameters["hill_cocaine"]
  
  range_HipO = parameters["range_HipO"] 
  HipO_constant = parameters["HipO_constant"] 
  hippurate_constant = parameters["hippurate_constant"]
  hill_HipO = parameters["hill_HipO"] 
  hill_hippurate = parameters["hill_hippurate"]
  
  adder_df = NULL
  for (hip in concentrations_hip) {
    for (coc in concentrations_coc) {
      HipO = hipo
      CocE = coce
      
      hip_total = transducer(enzyme = HipO, inducer = hip, range_enzyme = range_HipO, enzyme_constant = HipO_constant, 
                             inducer_constant = hippurate_constant, hill_enzyme = hill_HipO, hill_inducer = hill_hippurate)
      
      coc_total = transducer(enzyme = CocE, inducer = coc, range_enzyme = range_CocE, enzyme_constant = CocE_constant, 
                             inducer_constant = cocaine_constant, hill_enzyme = hill_CocE, hill_inducer = hill_cocaine)
      
      total = hip_total + coc_total
      value_here = transfer_function_linear(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline,
                                     slower_slope = slower_slope)
      # Competition
      adder_df = rbind(adder_df, c("cocaine" = coc,"hippurate" = hip, 'GFP' = value_here))
    }
  }
  adder_df = data.frame(adder_df)
  colnames(adder_df) = c("coc", "hip", 'coc_hip')
  return(adder_df[,c("coc_hip")])
}

calculate_hipo_hip_adder_linear = function(parameters, concentrations_hip = 0, concentrations_HipO = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  slower_slope = parameters["slower_slope"]
  # Enzyme parameters
  range_HipO = parameters["range_HipO"] 
  HipO_constant = parameters["HipO_constant"] 
  hippurate_constant = parameters["hippurate_constant"]
  hill_HipO = parameters["hill_HipO"] 
  hill_hippurate = parameters["hill_hippurate"]
  
  hip_HipO_adder_df = NULL
  for (hipo in concentrations_HipO) {
    for (hip in concentrations_hip) {
      HipO = hipo
      hip_total = transducer(enzyme = HipO, inducer = hip, range_enzyme = range_HipO, enzyme_constant = HipO_constant, 
                             inducer_constant = hippurate_constant, hill_enzyme = hill_HipO, hill_inducer = hill_hippurate)
      
      total = hip_total
      value_here = transfer_function_linear(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline,
                                     slower_slope = slower_slope)
      # Competition
      hip_HipO_adder_df = rbind(hip_HipO_adder_df, c("hippurate" = hip, "HipO" = hipo, 'GFP' = value_here))
    }
  }
  hip_HipO_adder_df = data.frame(hip_HipO_adder_df)
  colnames(hip_HipO_adder_df) = c("hip", "hipo", 'hipo_hip')
  return(hip_HipO_adder_df[,c("hipo_hip")])
}

calculate_coce_coc_adder_linear = function(parameters, concentrations_coc = 0, concentrations_CocE = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  slower_slope = parameters["slower_slope"]
  # Enzyme parameters
  range_CocE = parameters["range_CocE"] 
  CocE_constant = parameters["CocE_constant"] 
  cocaine_constant = parameters["cocaine_constant"]
  hill_CocE = parameters["hill_CocE"] 
  hill_cocaine = parameters["hill_cocaine"]
  
  coc_CocE_adder_df = NULL
  for (coce in concentrations_CocE) {
    for (coc in concentrations_coc) {
      CocE = coce
      coc_total = transducer(enzyme = CocE, inducer = coc, range_enzyme = range_CocE, enzyme_constant = CocE_constant, 
                             inducer_constant = cocaine_constant, hill_enzyme = hill_CocE, hill_inducer = hill_cocaine)
      
      total = coc_total
      value_here = transfer_function_linear(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline,
                                     slower_slope = slower_slope)
      # Competition
      coc_CocE_adder_df = rbind(coc_CocE_adder_df, c("cocaine" = coc, "CocE" = coce, 'GFP' = value_here))
    }
  }
  coc_CocE_adder_df = data.frame(coc_CocE_adder_df)
  colnames(coc_CocE_adder_df) = c("coc", "coce", 'coce_coc')
  return(coc_CocE_adder_df[,c("coce_coc")])
}

calculate_benzamid_enz_benzamid_adder_linear = function(parameters, concentrations_benzamid = 0, concentrations_benzamid_enz = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  slower_slope = parameters["slower_slope"]
  # Enzyme parameters
  range_benzamid_enz = parameters["range_benzamid_enz"] 
  benzamid_enz_constant = parameters["benzamid_enz_constant"] 
  benzamid_constant = parameters["benzamid_constant"]
  hill_benzamid_enz = parameters["hill_benzamid_enz"] 
  hill_benzamid = parameters["hill_benzamid"]
  
  benzamid_benzamid_enz_adder_df = NULL
  for (benzamid_enz in concentrations_benzamid_enz) {
    for (benzamid in concentrations_benzamid) {
      benzamid_total = transducer(enzyme = benzamid_enz, inducer = benzamid, range_enzyme = range_benzamid_enz, enzyme_constant = benzamid_enz_constant, 
                                  inducer_constant = benzamid_constant, hill_enzyme = hill_benzamid_enz, hill_inducer = hill_benzamid)
      
      total = benzamid_total
      value_here = transfer_function_linear(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline,
                                     slower_slope = slower_slope)
      # Competition
      benzamid_benzamid_enz_adder_df = rbind(benzamid_benzamid_enz_adder_df, c("benzamid" = benzamid, "benzamid_enz" = benzamid_enz, 'GFP' = value_here))
    }
  }
  benzamid_benzamid_enz_adder_df = data.frame(benzamid_benzamid_enz_adder_df)
  colnames(benzamid_benzamid_enz_adder_df) = c("benzamid", "benzamid_enz", 'benzamid_enz_benzamid')
  return(benzamid_benzamid_enz_adder_df[,c("benzamid_enz_benzamid")])
}

calculate_biphenyl_enz_biphenyl_adder_linear = function(parameters, concentrations_biphenyl = 0, concentrations_biphenyl_enz = 0) {
  # Biosensor parameters
  hill_transfer = parameters["hill_transfer"] 
  Km = parameters["Km"]
  fold_change = parameters["fold_change"]
  baseline = parameters["baseline"]
  slower_slope = parameters["slower_slope"]
  # Enzyme parameters
  range_biphenyl_enz = parameters["range_biphenyl_enz"] 
  biphenyl_enz_constant = parameters["biphenyl_enz_constant"] 
  biphenyl_constant = parameters["biphenyl_constant"]
  hill_biphenyl_enz = parameters["hill_biphenyl_enz"] 
  hill_biphenyl = parameters["hill_biphenyl"]
  
  biphenyl_biphenyl_enz_adder_df = NULL
  for (biphenyl_enz in concentrations_biphenyl_enz) {
    for (biphenyl in concentrations_biphenyl) {
      biphenyl_total = transducer(enzyme = biphenyl_enz, inducer = biphenyl, range_enzyme = range_biphenyl_enz, enzyme_constant = biphenyl_enz_constant, 
                                  inducer_constant = biphenyl_constant, hill_enzyme = hill_biphenyl_enz, hill_inducer = hill_biphenyl)
      
      total = biphenyl_total
      value_here = transfer_function_linear(total,
                                     hill_transfer = hill_transfer,
                                     Km = Km,
                                     fold_change = fold_change, 
                                     baseline = baseline,
                                     slower_slope = slower_slope)
      # Competition
      biphenyl_biphenyl_enz_adder_df = rbind(biphenyl_biphenyl_enz_adder_df, c("biphenyl" = biphenyl, "biphenyl_enz" = biphenyl_enz, 'GFP' = value_here))
    }
  }
  biphenyl_biphenyl_enz_adder_df = data.frame(biphenyl_biphenyl_enz_adder_df)
  colnames(biphenyl_biphenyl_enz_adder_df) = c("biphenyl", "biphenyl_enz", 'biphenyl_enz_biphenyl')
  return(biphenyl_biphenyl_enz_adder_df[,c("biphenyl_enz_biphenyl")])
}

# Calculations for the logic gate predictions

calculating_point_value <- function(conditions, chosen_parameter_set, highs, weights) {
  cocaine_quantity = conditions["coc"] * highs["cocaine_high"]
  cocaine_total = transducer(enzyme = weights["cocaine_weight"], inducer = cocaine_quantity, 
                             range_enzyme = chosen_parameter_set["range_CocE"], 
                             enzyme_constant = chosen_parameter_set["CocE_constant"], 
                             hill_enzyme = chosen_parameter_set["hill_CocE"],
                             inducer_constant = chosen_parameter_set["cocaine_constant"], 
                             hill_inducer = chosen_parameter_set["hill_cocaine"])
  
  
  hippurate_quantity = conditions["hip"] * highs["hippurate_high"]
  hippurate_total = transducer(enzyme = weights["hippurate_weight"], inducer = hippurate_quantity, 
                               range_enzyme = chosen_parameter_set["range_HipO"], 
                               enzyme_constant = chosen_parameter_set["HipO_constant"], 
                               hill_enzyme = chosen_parameter_set["hill_HipO"],
                               inducer_constant = chosen_parameter_set["hippurate_constant"], 
                               hill_inducer = chosen_parameter_set["hill_hippurate"])
  
  
  benzamid_quantity = conditions["benzamid"] * highs["benzamid_high"]
  benzamid_total = transducer(enzyme = weights["benzamid_weight"], inducer = benzamid_quantity, 
                              range_enzyme = chosen_parameter_set["range_benzamid_enz"], 
                              enzyme_constant = chosen_parameter_set["benzamid_enz_constant"], 
                              hill_enzyme = chosen_parameter_set["hill_benzamid_enz"],
                              inducer_constant = chosen_parameter_set["benzamid_constant"], 
                              hill_inducer = chosen_parameter_set["hill_benzamid"])
  
  
  biphenyl_quantity = conditions["biphenyl"] * highs["biphenyl_high"]
  biphenyl_total = transducer(enzyme = weights["biphenyl_weight"], inducer = biphenyl_quantity, 
                              range_enzyme = chosen_parameter_set["range_biphenyl_enz"], 
                              enzyme_constant = chosen_parameter_set["biphenyl_enz_constant"], 
                              hill_enzyme = chosen_parameter_set["hill_biphenyl_enz"],
                              inducer_constant = chosen_parameter_set["biphenyl_constant"], 
                              hill_inducer = chosen_parameter_set["hill_biphenyl"])
  total = cocaine_total + hippurate_total + benzamid_total + biphenyl_total
  value_here = transfer_function_linear(total, 
                                        hill_transfer = chosen_parameter_set["hill_transfer"],
                                        Km = chosen_parameter_set["Km"],
                                        fold_change = chosen_parameter_set["fold_change"],
                                        baseline = chosen_parameter_set["baseline"],
                                        slower_slope = chosen_parameter_set["slower_slope"])
  return(value_here)
}

calculating_point_value_cluster <- function(conditions, chosen_parameter_set, highs, weights) {
  cocaine_quantity = conditions["coc"] * highs["cocaine_high"]
  cocaine_total = transducer(enzyme = weights["cocaine_weight"], inducer = cocaine_quantity, 
                             range_enzyme = chosen_parameter_set["range_CocE"], 
                             enzyme_constant = chosen_parameter_set["CocE_constant"], 
                             hill_enzyme = chosen_parameter_set["hill_CocE"],
                             inducer_constant = chosen_parameter_set["cocaine_constant"], 
                             hill_inducer = chosen_parameter_set["hill_cocaine"])
  
  hippurate_quantity = conditions["hip"] * highs["hippurate_high"]
  hippurate_total = transducer(enzyme = weights["hippurate_weight"], inducer = hippurate_quantity, 
                               range_enzyme = chosen_parameter_set["range_HipO"], 
                               enzyme_constant = chosen_parameter_set["HipO_constant"], 
                               hill_enzyme = chosen_parameter_set["hill_HipO"],
                               inducer_constant = chosen_parameter_set["hippurate_constant"], 
                               hill_inducer = chosen_parameter_set["hill_hippurate"])
  benzamid_quantity = conditions["benzamid"] * highs["benzamid_high"]
  benzamid_total = transducer(enzyme = weights["benzamid_weight"], inducer = benzamid_quantity, 
                              range_enzyme = chosen_parameter_set["range_benzamid_enz"], 
                              enzyme_constant = chosen_parameter_set["benzamid_enz_constant"], 
                              hill_enzyme = chosen_parameter_set["hill_benzamid_enz"],
                              inducer_constant = chosen_parameter_set["benzamid_constant"], 
                              hill_inducer = chosen_parameter_set["hill_benzamid"])
  
  
  biphenyl_quantity = conditions["biphenyl"] * highs["biphenyl_high"]
  biphenyl_total = transducer(enzyme = weights["biphenyl_weight"], inducer = biphenyl_quantity, 
                              range_enzyme = chosen_parameter_set["range_biphenyl_enz"], 
                              enzyme_constant = chosen_parameter_set["biphenyl_enz_constant"], 
                              hill_enzyme = chosen_parameter_set["hill_biphenyl_enz"],
                              inducer_constant = chosen_parameter_set["biphenyl_constant"], 
                              hill_inducer = chosen_parameter_set["hill_biphenyl"])
  
  total = cocaine_total + hippurate_total + benzamid_total + biphenyl_total
  value_here = transfer_function_linear(total, 
                                        hill_transfer = chosen_parameter_set["hill_transfer"],
                                        Km = chosen_parameter_set["Km"],
                                        fold_change = chosen_parameter_set["fold_change"],
                                        baseline = chosen_parameter_set["baseline"],
                                        slower_slope = chosen_parameter_set["slower_slope"])
  return(value_here)
}

# Clusters for 2D gates

calculating_point_value_2D_cluster <- function(conditions, chosen_parameter_set, highs, weights) {
  cocaine_quantity = conditions["coc"] * highs["cocaine_high"]
  cocaine_total = transducer(enzyme = weights["cocaine_weight"], inducer = cocaine_quantity, 
                             range_enzyme = chosen_parameter_set["range_CocE"], 
                             enzyme_constant = chosen_parameter_set["CocE_constant"], 
                             hill_enzyme = chosen_parameter_set["hill_CocE"],
                             inducer_constant = chosen_parameter_set["cocaine_constant"], 
                             hill_inducer = chosen_parameter_set["hill_cocaine"])
  
  hippurate_quantity = conditions["hip"] * highs["hippurate_high"]
  hippurate_total = transducer(enzyme = weights["hippurate_weight"], inducer = hippurate_quantity, 
                               range_enzyme = chosen_parameter_set["range_HipO"], 
                               enzyme_constant = chosen_parameter_set["HipO_constant"], 
                               hill_enzyme = chosen_parameter_set["hill_HipO"],
                               inducer_constant = chosen_parameter_set["hippurate_constant"], 
                               hill_inducer = chosen_parameter_set["hill_hippurate"])
  
  total = cocaine_total + hippurate_total
  # total = cocaine_total + hippurate_total + benzamid_total + biphenyl_total
  value_here = transfer_function_linear(total, 
                                        hill_transfer = chosen_parameter_set["hill_transfer"],
                                        Km = chosen_parameter_set["Km"],
                                        fold_change = chosen_parameter_set["fold_change"],
                                        baseline = chosen_parameter_set["baseline"],
                                        slower_slope = chosen_parameter_set["slower_slope"])
  return(value_here)
}

line_from_boundary_vertical <-function(boundary) {
  x_vector = seq(0, 1, by = 0.01)
  y_vector =  seq(0, 1, by = 0.001)
  some_impossible_boundaries = FALSE
  y_result = c()
  x_result = c()
  # Solving by y , makes more sense given what it's going to look like in the end.
  for (y in y_vector) {
    hippurate_quantity = y * highs["hippurate_high"]
    hippurate_total = transducer(enzyme = weights["hippurate_weight"], inducer = hippurate_quantity, 
                                 range_enzyme = chosen_parameter_set["range_HipO"], 
                                 enzyme_constant = chosen_parameter_set["HipO_constant"], 
                                 hill_enzyme = chosen_parameter_set["hill_HipO"],
                                 inducer_constant = chosen_parameter_set["hippurate_constant"], 
                                 hill_inducer = chosen_parameter_set["hill_hippurate"])
    aim = boundary - hippurate_total
    if (aim <0) {
      some_impossible_boundaries = TRUE
      # print(paste("No value for x =", x, sep = ''))
    } else {
      best_x = 0
      best_x_value = 1000
      for (x in x_vector) {
        cocaine_quantity = x * highs["cocaine_high"]
        cocaine_total = transducer(enzyme = weights["cocaine_weight"], inducer = cocaine_quantity, 
                                   range_enzyme = chosen_parameter_set["range_CocE"], 
                                   enzyme_constant = chosen_parameter_set["CocE_constant"], 
                                   hill_enzyme = chosen_parameter_set["hill_CocE"],
                                   inducer_constant = chosen_parameter_set["cocaine_constant"], 
                                   hill_inducer = chosen_parameter_set["hill_cocaine"])
        
        
        if (abs(cocaine_total-aim) < abs(best_x_value - aim)) {
          best_x_value = cocaine_total
          best_x = x
          # print(paste("Aim is ", aim, "and hipp total is ", hippurate_total, sep = ' '))
        }
      }
      y_result = c(y_result, y)
      x_result = c(x_result, best_x)
    }
  }
  return(cbind(x_result, y_result))
}


line_from_boundary_horizontal <-function(boundary) {
  x_vector = seq(0, 1, by = 0.001)
  y_vector =  seq(0, 1, by = 0.001)
  some_impossible_boundaries = FALSE
  y_result = c()
  x_result = c()
  for (x in x_vector) {
    cocaine_quantity = x * highs["cocaine_high"]
    cocaine_total = transducer(enzyme = weights["cocaine_weight"], inducer = cocaine_quantity, 
                               range_enzyme = chosen_parameter_set["range_CocE"], 
                               enzyme_constant = chosen_parameter_set["CocE_constant"], 
                               hill_enzyme = chosen_parameter_set["hill_CocE"],
                               inducer_constant = chosen_parameter_set["cocaine_constant"], 
                               hill_inducer = chosen_parameter_set["hill_cocaine"])
    
    aim = boundary - cocaine_total
    if (aim <0) {
      some_impossible_boundaries = TRUE
      # print(paste("No value for x =", x, sep = ''))
    } else {
      best_y = 0
      best_y_value = 1000
      for (y in y_vector) {
        hippurate_quantity = y * highs["hippurate_high"]
        hippurate_total = transducer(enzyme = weights["hippurate_weight"], inducer = hippurate_quantity, 
                                     range_enzyme = chosen_parameter_set["range_HipO"], 
                                     enzyme_constant = chosen_parameter_set["HipO_constant"], 
                                     hill_enzyme = chosen_parameter_set["hill_HipO"],
                                     inducer_constant = chosen_parameter_set["hippurate_constant"], 
                                     hill_inducer = chosen_parameter_set["hill_hippurate"])
        
        if (abs(hippurate_total-aim) < abs(best_y_value - aim)) {
          best_y_value = hippurate_total
          best_y = y
          # print(paste("Aim is ", aim, "and hipp total is ", hippurate_total, sep = ' '))
        }
      }
      y_result = c(y_result, best_y)
      x_result = c(x_result, x)
    }
  }
  return(cbind(x_result, y_result))
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
    saving_graph(filename = name_for_saving, path = folder_to_save)
  }
}

# Plotting the 2D heatmap
print_heatmap <- function(data_as_factor, x_axis = "coc", y_axis = "hip", x_lab = "Compound 1", y_lab ="Compound 2", title = "", subtitle = "", 
                          print = TRUE, vertical = FALSE) {
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  p <- ggplot(data = data_as_factor, aes(x= data_as_factor[,x_axis], y= data_as_factor[,y_axis], fill=value))
  p <- p + geom_tile()
  p <- p + scale_fill_gradientn(colours = jet.colors(7))
  p <- p + labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  p <- p + theme_bw() #theme black and white
  if (vertical) {
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }

  # p <- p + theme_bw(legend.direction="vertical")
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
  p <- p + theme(legend.position=legend, plot.title = element_text(hjust = 0.5))
  
  p
  if (!is.null(folder_to_save)) {
    saving_graph(filename = title_to_save, path =folder_to_save)
  }
  print(p)
}


plot_logic_gates <- function(results_dataframe_for_plotting, testing_threshold = NULL, 
                             name_for_saving = "logic_gates", title = "Logic gates",
                             folder_to_save = NULL, saving = TRUE) {
  # For plotting logic gates predictions with various thresholds
  results_dataframe_for_plotting = results_dataframe[c("names", "expected", "value", testing_threshold)]
  
  heatmap_columns = as.character(results_dataframe_for_plotting[,c("names")])
  results_dataframe_for_plotting = melt(results_dataframe_for_plotting, id.vars = "names")
  heatmap_rownames = c("expected", "value", testing_threshold)
  
  results_dataframe_for_plotting[,"names"] <- as.factor(results_dataframe_for_plotting[,c("names")])
  results_dataframe_for_plotting[,"variable"] <- as.factor(results_dataframe_for_plotting[,c("variable")])
  
  xlab = "Thresholds"
  ylab = "Combinations"
  
  print_heatmap(data_as_factor = results_dataframe_for_plotting, x_lab = xlab, y_lab = ylab, x_axis = "variable", 
                y_axis = "names", title = title, subtitle = "", vertical = TRUE)
  
  if (saving) {  
    saving_graph(filename = paste(name_for_saving, '.jpeg', sep = ''), path = folder_to_save)
  }
  
}

# Saves the laast graph to disk
saving_graph <- function(plot = last_plot(),filename, path = getwd(), device = "jpeg", 
                         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                         dpi = 300, limitsize = TRUE) {
  
  ggsave(filename = filename, plot = plot, device = device, path = path, scale = scale, width = width, height = height, units = units,
         dpi = dpi, limitsize = limitsize)
}

