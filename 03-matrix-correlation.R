library(data.table)
library(corrplot)

# Load the dataset
message("Loading dataset...")
simstats <- fread("00-data/02-dataprep-results-2.csv.gz")

# Filter valid simulations (consistent with analysis scripts)
simstats <- simstats[peak_preval > 1]

# Select all numeric variables of interest
# We exclude IDs and filenames. We include outcomes and network metrics.
vars_of_interest <- c(
  # Outcomes
  "peak_time", "peak_preval", "rt", "rt_mean", "gentime", "final_preval",
  
  # ERGM metrics
  "ergm_balance", "ergm_triangles", "ergm_twopath",
  
  # igraph metrics
  "igraph_avg_degree", "igraph_avg_path_length", 
  "igraph_local_transitivity", "igraph_transitivity", 
  "igraph_modularity", "igraph_degree_moment_2", 
  "igraph_degree_variance", "igraph_assortativity", 
  "igraph_spectral_radius", "igraph_avg_betweenness", 
  "igraph_diameter", "igraph_avg_closeness", 
  "igraph_avg_eigenvector", "igraph_density"
)

# Check if all columns exist
missing_vars <- setdiff(vars_of_interest, names(simstats))
if(length(missing_vars) > 0) {
  warning(paste("Missing variables:", paste(missing_vars, collapse = ", ")))
  vars_of_interest <- intersect(vars_of_interest, names(simstats))
}

# Subset data
data_for_cor <- simstats[, ..vars_of_interest]

# Calculate correlation matrix
message("Calculating correlation matrix...")
cor_matrix <- cor(data_for_cor, use = "complete.obs")

# Save correlation matrix to CSV for detailed inspection
write.csv(cor_matrix, "03-matrix-correlation-all-values.csv")
message("Correlation matrix values saved to '03-matrix-correlation-all-values.csv'")

# Plotting
output_pdf <- "03-matrix-correlation-all.pdf"
pdf(output_pdf, width = 12, height = 12)

# Improve readability for the plot
# Shorten names for the plot
colnames(cor_matrix) <- gsub("igraph_", "", colnames(cor_matrix))
colnames(cor_matrix) <- gsub("ergm_", "", colnames(cor_matrix))
rownames(cor_matrix) <- colnames(cor_matrix)

corrplot(cor_matrix, 
         method = "color", 
         type = "lower", 
         diag = FALSE,
         tl.col = "black", 
         tl.srt = 45, 
         addCoef.col = "black", # Add numbers
         number.cex = 0.5,      # Smaller numbers to fit
         tl.cex = 0.8,          # Smaller text labels
         cl.ratio = 0.1)

dev.off()
message(paste("Correlation plot saved to", output_pdf))
