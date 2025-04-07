all_variables <- c(
  "ergm_balance", "ergm_triangle", "ergm_twopath", "igraph_avg_path_length",
  "igraph_avg_degree", "igraph_transitivity", "igraph_modularity"
)

all_variables_combined <- c(
  "ergm_balance", "I(ergm_balance^2)", "log(ergm_balance)",
  "ergm_triangle", "I(ergm_triangle^2)", "log(ergm_triangle)",
  "ergm_twopath", "I(ergm_twopath^2)", "log(ergm_twopath)",
  "igraph_avg_path_length", "I(igraph_avg_path_length^2)", "log(igraph_avg_path_length)",
  "igraph_avg_degree", "I(igraph_avg_degree^2)", "log(igraph_avg_degree)",
  "igraph_transitivity", "I(igraph_transitivity^2)", "log(igraph_transitivity)",
  "igraph_modularity", "I(igraph_modularity^2)", "log(igraph_modularity)"
)

# Function to calculate percentiles for a list of target models
find_models_percentiles <- function(target_formulas, all_models) {
  
  # Extract AIC, BIC, and formulas from all models
  aic_values <- unlist(lapply(all_models, function(x) x$aic))
  bic_values <- unlist(lapply(all_models, function(x) x$bic))
  formulas_all_models <- sapply(all_models, function(x) as.character(x$formula)[3])
  formulas_all_models <- na.omit(formulas_all_models)
  #formulas_all_models <- formulas_all_models[3,]
  
  formulas_all_models <- Filter(Negate(is.null), formulas_all_models)
  
  # Create a data frame for all models
  all_models_df <- data.frame(
    aic = aic_values,
    bic = bic_values,
    formula = formulas_all_models
  )
  
  # Normalize AIC and BIC values
  mean_AIC <- mean(aic_values)
  mean_BIC <- mean(bic_values)
  sd_AIC <- sd(aic_values)
  sd_BIC <- sd(bic_values)
  
  all_models_df$AIC_z <- (all_models_df$aic - mean_AIC) / sd_AIC
  all_models_df$BIC_z <- (all_models_df$bic - mean_BIC) / sd_BIC
  
  # Calculate the mean z-score for each model
  all_models_df$Mean_z <- rowMeans(all_models_df[, c("AIC_z", "BIC_z")])
  
  # Initialize a list to store percentiles for target models
  target_percentiles <- list()
  
  # Loop through each target formula and find its percentile
  for (target_formula in target_formulas) {
    target_index <- which(all_models_df$formula == target_formula)
    
    if (length(target_index) == 0) {
      warning(paste("The specified model formula was not found:", target_formula))
      next
    }
    
    # Get the mean z-score of the target model
    target_mean_z <- all_models_df$Mean_z[target_index]
    
    # Calculate the percentile based on the number of models with a higher or equal z-score
    percentile <- sum(all_models_df$Mean_z >= target_mean_z) / nrow(all_models_df) * 100
    
    # Store the percentile in the list with the formula as the name
    target_percentiles[[target_formula]] <- percentile
  }
  
  return(target_percentiles)
}

################################## Peak Preval Models ------------------------------

peak_preval_combined_file <- 'data/03-model-comparison-combined/03-combined_models_peak_preval.RData'

if (!file.exists(peak_preval_combined_file)) {
  stop('Model data file does not exist.')
} else {
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
}

# Define a list of target formulas
peak_preval_target_formulas <- c(
  "I(factor(nettype)) + ergm_twopath + igraph_avg_path_length + igraph_avg_degree",
  "I(factor(nettype)) + ergm_balance + ergm_twopath + igraph_avg_path_length",
  "I(factor(nettype)) + ergm_triangle + ergm_twopath + igraph_avg_path_length + igraph_avg_degree",
  "I(factor(nettype)) + ergm_balance + ergm_triangle + ergm_twopath + igraph_avg_path_length",
  "I(factor(nettype)) + ergm_twopath + igraph_avg_path_length + igraph_avg_degree + igraph_transitivity",
  "I(factor(nettype)) + ergm_balance + ergm_twopath + igraph_avg_path_length + igraph_transitivity"
)

# Get percentiles for the list of target formulas
peak_preval_formulas_percentiles <- find_models_percentiles(peak_preval_target_formulas, peak_preval_combined_models)

# Print percentiles
print(peak_preval_formulas_percentiles)

################################## Peak Time Models ------------------------------

peak_time_combined_file <- 'data/03-model-comparison-combined/03-combined_models_peak_time.RData'

if (!file.exists(peak_time_combined_file)) {
  stop('Model data file does not exist.')
} else {
  peak_time_combined_models <- readRDS(peak_time_combined_file)
}

# Define a list of target formulas
peak_time_target_formulas <- c(
  "I(factor(nettype)) + ergm_balance + igraph_avg_degree",
  "I(factor(nettype)) + ergm_balance + igraph_avg_path_length",
  "I(factor(nettype)) + igraph_avg_path_length + igraph_avg_degree"
)

# Get percentiles for the list of target formulas
peak_time_formulas_percentiles <- find_models_percentiles(peak_time_target_formulas, peak_time_combined_models)

# Print percentiles
print(peak_time_formulas_percentiles)

################################## Gentime Models ------------------------------

gentime_combined_file <- 'data/03-model-comparison-combined/03-combined_models_gentime.RData'

if (!file.exists(gentime_combined_file)) {
  stop('Model data file does not exist.')
} else {
  gentime_combined_models <- readRDS(gentime_combined_file)
}

# Define a list of target formulas
gentime_target_formulas <- c(
  "I(factor(nettype)) + ergm_balance + igraph_avg_degree",
  "I(factor(nettype)) + ergm_balance + igraph_avg_path_length",
  "I(factor(nettype)) + igraph_avg_path_length + igraph_avg_degree"
)

# Get percentiles for the list of target formulas
gentime_formulas_percentiles <- find_models_percentiles(gentime_target_formulas, gentime_combined_models)

# Print percentiles
print(gentime_formulas_percentiles)

################################## Rep Num Models ------------------------------

rep_num_combined_file <- 'data/03-model-comparison-combined/03-combined_models_rep_num.RData'

if (!file.exists(rep_num_combined_file)) {
  stop('Model data file does not exist.')
} else {
  rep_num_combined_models <- readRDS(rep_num_combined_file)
}

# Define a list of target formulas
rep_num_target_formulas <- c(
  "I(factor(nettype)) + ergm_balance + igraph_avg_degree",
  "I(factor(nettype)) + ergm_balance + igraph_avg_path_length",
  "I(factor(nettype)) + igraph_avg_path_length + igraph_avg_degree"
)

# Get percentiles for the list of target formulas
rep_num_formulas_percentiles <- find_models_percentiles(rep_num_target_formulas, rep_num_combined_models)

# Print percentiles
print(rep_num_formulas_percentiles)
