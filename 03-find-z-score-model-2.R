
all_variables <- c(
  "ergm_balance", "ergm_triangles", "igraph_avg_path_length",
  "igraph_avg_degree", "igraph_local_transitivity", "igraph_modularity",
  "igraph_degree_variance", "igraph_assortativity", "igraph_spectral_radius",
  "igraph_avg_betweenness"
)

all_variables_combined <- c(
  "ergm_balance", "I(ergm_balance^2)", "log(ergm_balance)",
  "ergm_triangles", "I(ergm_triangles^2)", "log(ergm_triangles)",
  "igraph_avg_path_length", "I(igraph_avg_path_length^2)", "log(igraph_avg_path_length)",
  "igraph_avg_degree", "I(igraph_avg_degree^2)", "log(igraph_avg_degree)",
  "igraph_local_transitivity", "I(igraph_local_transitivity^2)", "log(igraph_local_transitivity)",
  "igraph_modularity", "I(igraph_modularity^2)", "log(igraph_modularity)",
  "igraph_degree_variance", "I(igraph_degree_variance^2)", "log(igraph_degree_variance)",
  "igraph_assortativity", "I(igraph_assortativity^2)", "log(igraph_assortativity)",
  "igraph_spectral_radius", "I(igraph_spectral_radius^2)", "log(igraph_spectral_radius)",
  "igraph_avg_betweenness", "I(igraph_avg_betweenness^2)", "log(igraph_avg_betweenness)"
)

# Function to calculate percentiles for a list of target models
find_models_percentiles <- function(target_formulas, all_models) {
  
  # Extract AIC, BIC, and formulas from all models
  aic_values <- unlist(lapply(all_models, function(x) x$aic))
  bic_values <- unlist(lapply(all_models, function(x) x$bic))
  formulas_all_models <- sapply(all_models, function(x) as.character(x$formula)[3])
  formulas_all_models <- na.omit(formulas_all_models)
  
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
