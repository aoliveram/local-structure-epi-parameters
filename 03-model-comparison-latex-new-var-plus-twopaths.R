library(parallel)
library(doParallel)
library(foreach)
library(data.table)
library(ggplot2)
library(lmtest)
library(texreg)
library(dplyr)
library(MASS)
library(grid)
library(stargazer)
library(corrplot)
library(car) # For GVIF
library(reshape2)

# Start Global Timer
global_start_time <- Sys.time()
message(paste("Analysis started at:", global_start_time))

################################## Data & Functions ----------------------------

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
simstats$igraph_transitivity <- simstats$igraph_transitivity * 100 # Added for consistency
simstats$igraph_modularity <- simstats$igraph_modularity * 100
simstats$igraph_assortativity <- simstats$igraph_assortativity * 100

# Renaming
vallabs <- c(
  peak_time = "Peak time",
  peak_preval = "Peak prevalence",
  rt = "Reproductive number",
  rt_0 = "Reproductive number",
  rt_mean = "Average reproductive number",
  dispersion = "Dispersion",
  gentime = "Generation time",
  final_preval = "Final prevalence"
)

simstats[, nettype := factor(
  nettype,
  levels = c("ergm", "sf", "swp01", "swp02", "degseq", "er"),
  labels = c("ERGM", "Scale-free", "Small-world (p=0.1)", "Small-world (p=0.2)", "Degree-sequence", "Erdos-Renyi")
)]

# The New Selected Variables + Two-Paths
all_variables <- c(
  "igraph_avg_degree",
  "igraph_avg_path_length",
  "igraph_local_transitivity",
  "igraph_modularity",
  "igraph_assortativity",
  "igraph_degree_variance",
  "ergm_twopath"
)

all_variables_combined <- c(
  "igraph_avg_degree", "I(igraph_avg_degree^2)", "log(igraph_avg_degree)",
  "igraph_avg_path_length", "I(igraph_avg_path_length^2)", "log(igraph_avg_path_length)",
  "igraph_local_transitivity", "I(igraph_local_transitivity^2)", "log(igraph_local_transitivity)",
  "igraph_modularity", "I(igraph_modularity^2)", # Skipped log because it can be negative
  "igraph_assortativity", "I(igraph_assortativity^2)", # Skipped log because it can be negative
  "igraph_degree_variance", "I(igraph_degree_variance^2)", "log(igraph_degree_variance)",
  "ergm_twopath", "I(ergm_twopath^2)", "log(ergm_twopath)"
)

# Function to select best models for table
Q_best_models_table <- function(Q_values, all_models, all_variables_combined, porcent) {
  
  aic_values <- unlist(lapply(all_models, function(x) x$aic))
  bic_values <- unlist(lapply(all_models, function(x) x$bic))
  
  formulas_all_models <- as.character(unlist(lapply(all_models, function(x) x$formula)))
  
  all_models_df <- as.data.frame(list(
    aic = aic_values,
    bic = bic_values,
    formula = formulas_all_models
  ))
  
  # Normalize AIC and BIC values
  mean_AIC <- mean(aic_values)
  mean_BIC <- mean(bic_values)
  sd_AIC <- sd(aic_values)
  sd_BIC <- sd(bic_values)
  
  all_models_df$AIC_z <- (all_models_df$aic - mean_AIC) / sd_AIC
  all_models_df$BIC_z <- (all_models_df$bic - mean_BIC) / sd_BIC
  
  # Mean of both normalized AIC and BIC
  all_models_df$Mean_z <- rowMeans(all_models_df[, c("AIC_z", "BIC_z")])
  all_models_df <- all_models_df[order(all_models_df$Mean_z), ]
  
  # Selecting Q% best models
  top_Q_models <- list()
  
  for (Q in Q_values) {
    if (porcent) {
      num_models <- floor((Q / 100) * nrow(all_models_df)) 
    } else {
      num_models <- Q
    }
    
    top_Q_models <- all_models_df$formula[1:num_models]
  }
  
  return(top_Q_models)
}

################################## Peak Preval Models ----------------------------

peak_preval_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_peak_preval.RData'

if (file.exists(peak_preval_combined_file)) {
  print('Loading peak_preval_combined_models')
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
  peak_preval_combined_models <- Filter(Negate(is.null), peak_preval_combined_models)
  
  # Table
  Q_value <- c(10); porcent <- FALSE
  Q_best_formulas_peak_preval <- Q_best_models_table(Q_value, peak_preval_combined_models, all_variables_combined, porcent)
  
  Q_best_models_peak_preval <- lapply(Q_best_formulas_peak_preval, function(formula) {
    glm(as.formula(formula), data = simstats)
  }) # Runing best models
  
  aic_values_peak_preval <- round(sapply(Q_best_models_peak_preval, AIC), 3)
  bic_values_peak_preval <- round(sapply(Q_best_models_peak_preval, BIC), 3)
  
  tex_file_peak_preval <- '03-model-comparison-combined-new-var-plus-twopaths/03-best_models_table_peak_preval.tex'
  
  stargazer(Q_best_models_peak_preval, type = "latex", out = tex_file_peak_preval,
            title = "Best 10 models for Peak Prevalence",
            label = "tab:best_models_peak_preval",
            dep.var.labels = "Peak Prevalence",
            # covariate.labels = best_variables_peak_preval,
            omit.stat = c("LL", "ser", "f", "aic"),
            omit.table.layout = "d", # not showing dependent var
            add.lines = list(c("AIC", aic_values_peak_preval),
                             c("BIC", bic_values_peak_preval)),
            no.space = TRUE,
            digits = 3)
} else {
  print(paste("File not found:", peak_preval_combined_file))
}

################################## Peak Time Models ----------------------------

peak_time_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_peak_time.RData'

if (file.exists(peak_time_combined_file)) {
  print('Loading peak_time_combined_models')
  peak_time_combined_models <- readRDS(peak_time_combined_file)
  peak_time_combined_models <- Filter(Negate(is.null), peak_time_combined_models)
  
  # Table
  Q_value <- c(10); porcent <- FALSE
  Q_best_formulas_peak_time <- Q_best_models_table(Q_value, peak_time_combined_models, all_variables_combined, porcent)
  
  Q_best_models_peak_time <- lapply(Q_best_formulas_peak_time, function(formula) {
    glm(as.formula(formula), data = simstats)
  }) # Runing best models
  
  aic_values_peak_time <- round(sapply(Q_best_models_peak_time, AIC), 3)
  bic_values_peak_time <- round(sapply(Q_best_models_peak_time, BIC), 3)
  
  tex_file_peak_time <- '03-model-comparison-combined-new-var-plus-twopaths/03-best_models_table_peak_time.tex'
  
  stargazer(Q_best_models_peak_time, type = "latex", out = tex_file_peak_time,
            title = "Best 10 models for Peak Time",
            label = "tab:best_models_peak_time",
            dep.var.labels = "Peak Time",
            # covariate.labels = best_variables_peak_time,
            omit.stat = c("LL", "ser", "f", "aic"),
            omit.table.layout = "d", # not showing dependent var
            add.lines = list(c("AIC", aic_values_peak_time),
                             c("BIC", bic_values_peak_time)),
            no.space = TRUE,
            digits = 3)
} else {
  print(paste("File not found:", peak_time_combined_file))
}

################################## Gentime Models ------------------------------

gentime_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_gentime.RData'

if (file.exists(gentime_combined_file)) {
  print('Loading gentime_combined_models')
  gentime_combined_models <- readRDS(gentime_combined_file)
  gentime_combined_models <- Filter(Negate(is.null), gentime_combined_models)
  
  # Table
  Q_value <- c(10); porcent <- FALSE
  Q_best_formulas_gentime <- Q_best_models_table(Q_value, gentime_combined_models, all_variables_combined, porcent)
  
  Q_best_models_gentime <- lapply(Q_best_formulas_gentime, function(formula) {
    glm(as.formula(formula), data = simstats)
  }) # Runing best models
  
  aic_values_gentime <- round(sapply(Q_best_models_gentime, AIC), 3)
  bic_values_gentime <- round(sapply(Q_best_models_gentime, BIC), 3)
  
  tex_file_gentime <- '03-model-comparison-combined-new-var-plus-twopaths/03-best_models_table_gentime.tex'
  
  stargazer(Q_best_models_gentime, type = "latex", out = tex_file_gentime,
            title = "Best 10 models for Generation Time",
            label = "tab:best_models_gentime",
            dep.var.labels = "Generation Time",
            # covariate.labels = best_variables_gentime,
            omit.stat = c("LL", "ser", "f", "aic"),
            omit.table.layout = "d", # not showing dependent var
            add.lines = list(c("AIC", aic_values_gentime),
                             c("BIC", bic_values_gentime)),
            no.space = TRUE,
            digits = 3)
} else {
  print(paste("File not found:", gentime_combined_file))
}

################################## Rep Num Models ------------------------------

rep_num_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_rep_num.RData'

if (file.exists(rep_num_combined_file)) {
  print('Loading rep_num_combined_models')
  rep_num_combined_models <- readRDS(rep_num_combined_file)
  rep_num_combined_models <- Filter(Negate(is.null), rep_num_combined_models)
  
  # Table
  Q_value <- c(10); porcent <- FALSE
  Q_best_formulas_rep_num <- Q_best_models_table(Q_value, rep_num_combined_models, all_variables_combined, porcent)
  
  Q_best_models_rep_num <- lapply(Q_best_formulas_rep_num, function(formula) {
    glm.nb(as.formula(formula), data = simstats)
  }) # Runing best models
  
  aic_values_rep_num <- round(sapply(Q_best_models_rep_num, AIC), 3)
  bic_values_rep_num <- round(sapply(Q_best_models_rep_num, BIC), 3)
  
  tex_file_rep_num <- '03-model-comparison-combined-new-var-plus-twopaths/03-best_models_table_rep_num.tex'
  
  stargazer(Q_best_models_rep_num, type = "latex", out = tex_file_rep_num,
            title = "Best 10 models for Reproductive Number",
            label = "tab:best_models_rep_num",
            dep.var.labels = "Reproductive Number",
            # covariate.labels = best_variables_rep_num,
            omit.stat = c("LL", "ser", "f", "aic"),
            omit.table.layout = "d", # not showing dependent var
            add.lines = list(c("AIC", aic_values_rep_num),
                             c("BIC", bic_values_rep_num)),
            no.space = TRUE,
            digits = 3)
} else {
  print(paste("File not found:", rep_num_combined_file))
}

# End Global Timer
global_end_time <- Sys.time()
message(paste("Analysis finished at:", global_end_time))
message(paste("Total duration:", round(difftime(global_end_time, global_start_time, units = "mins"), 2), "minutes"))
