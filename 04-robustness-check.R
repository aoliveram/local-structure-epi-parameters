library(randomForest)
library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)

# Start Global Timer
global_start_time <- Sys.time()
message(paste("Robustness Check (Random Forest) started at:", global_start_time))

# Create output directory
output_dir <- "04-robustness-check"
dir.create(output_dir, showWarnings = FALSE)

################################## Data Preparation ----------------------------
# (Copied from 03-model-comparison-combined-new-var.R to ensure consistency)

# Load the new dataset
simstats <- fread("00-data/02-dataprep-results-2.csv.gz")
simstats <- simstats[peak_preval > 1]

# Rescaling variables
n <- 534
twopath_complete <- (n * (n - 1) * (n - 2)) / 2
balance_complete <- ergm_triangle <- (n * (n - 1) * (n - 2)) / 6
triangle_complete <- balance_complete

# Rescale ERGM counts to percentages of max possible
simstats$ergm_twopath <- (simstats$ergm_twopath / twopath_complete) * 100
simstats$ergm_balance <- (simstats$ergm_balance / balance_complete) * 100
simstats$ergm_triangles <- (simstats$ergm_triangles / triangle_complete) * 100

# Rescale 0-1 metrics to 0-100 for better coefficient readability
simstats$igraph_local_transitivity <- simstats$igraph_local_transitivity * 100
simstats$igraph_transitivity <- simstats$igraph_transitivity * 100
simstats$igraph_modularity <- simstats$igraph_modularity * 100
# simstats$igraph_assortativity <- simstats$igraph_assortativity * 100 # Removed as per previous script

# Rescale Variance (Divide by 100 to make coefficients comparable)
simstats$igraph_degree_variance <- simstats$igraph_degree_variance / 100

simstats[, nettype := factor(
  nettype,
  levels = c("ergm", "sf", "swp01", "swp02", "degseq", "er"),
  labels = c("ERGM", "Scale-free", "Small-world (p=0.1)", "Small-world (p=0.2)", "Degree-sequence", "Erdos-Renyi")
)]

# Define Variables
dependent_vars <- c("peak_preval", "peak_time", "gentime", "rt_0")
dependent_vars_labels <- c("Peak Prevalence", "Peak Time", "Generation Time", "Reproductive Number")

predictors <- c(
  "igraph_avg_degree",
  "igraph_avg_path_length",
  "igraph_local_transitivity",
  "igraph_modularity",
  "igraph_assortativity",
  "igraph_degree_variance",
  "nettype"
)

predictors_labels <- c(
  "Avg. Degree",
  "Avg. Path Length",
  "Local Transitivity",
  "Modularity",
  "Assortativity",
  "Degree Variance",
  "Network Type"
)

names(predictors_labels) <- predictors

################################## Random Forest Analysis ----------------------

run_random_forest <- function(dep_var, label) {
  message(paste("Running Random Forest for:", label))
  
  # Prepare formula
  rf_formula <- as.formula(paste(dep_var, "~", paste(predictors, collapse = " + ")))
  
  # Run Random Forest
  set.seed(123) # For reproducibility
  rf_model <- randomForest(rf_formula, data = simstats, importance = TRUE, ntree = 500, na.action = na.omit)
  
  # Print summary
  print(rf_model)
  
  # Extract Importance
  imp_df <- as.data.frame(importance(rf_model))
  imp_df$Variable <- rownames(imp_df)
  
  # Use %IncMSE (Percentage Increase in MSE) as the importance metric
  # This shows how much the model accuracy decreases if we permute the variable.
  # Higher is better.
  
  # Clean variable names for plotting
  imp_df$Label <- predictors_labels[imp_df$Variable]
  
  # Plot Variable Importance
  p <- ggplot(imp_df, aes(x = reorder(Label, `%IncMSE`), y = `%IncMSE`)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("Variable Importance (Random Forest) -", label),
      subtitle = paste("Variance Explained:", round(rf_model$rsq[length(rf_model$rsq)] * 100, 2), "%"),
      x = "Variable",
      y = "% Increase in MSE (Importance)"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text = element_text(size = 10)
    )
  
  # Save Plot
  ggsave(file.path(output_dir, paste0("RF_Importance_", dep_var, ".pdf")), p, width = 8, height = 6)
  
  return(rf_model)
}

# Run for all dependent variables
rf_results <- list()

for (i in seq_along(dependent_vars)) {
  rf_results[[dependent_vars[i]]] <- run_random_forest(dependent_vars[i], dependent_vars_labels[i])
}

# Save RF models
saveRDS(rf_results, file.path(output_dir, "rf_models.rds"))

message("Random Forest analysis completed.")
