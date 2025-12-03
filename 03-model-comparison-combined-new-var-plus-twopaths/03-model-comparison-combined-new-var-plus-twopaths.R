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
# simstats$igraph_assortativity <- simstats$igraph_assortativity * 100 # Removed to keep coefficients larger

# Rescale Variance (Divide by 100 to make coefficients comparable)
simstats$igraph_degree_variance <- simstats$igraph_degree_variance / 100

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

# Let's see the factors betwwen the minimum and maximum
variables_metadata <- data.frame(t(
  sapply(all_variables, function(var) {
    max_val <- max(simstats[[var]], na.rm = TRUE)
    min_val <- min(simstats[[var]], na.rm = TRUE)
    mean_val <- mean(simstats[[var]], na.rm = TRUE)
    sd_val <- sd(simstats[[var]], na.rm = TRUE)
    c(Mean = mean_val, SD = sd_val, Factor = max_val / min_val)})
))

# Create output directory
dir.create("03-model-comparison-combined-new-var-plus-twopaths", showWarnings = FALSE)

# correlation matrix
matrix_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-matrix_correlation.pdf'

if (!file.exists(matrix_file)) {
  correlation_matrix <- cor(simstats[, ..all_variables], use = "complete.obs")
  # Update labels for the new set
  all_variables_labels <- c('Avg. Degree', 'Avg. Path Length', 'Local Transitivity', 
                            'Modularity', 'Assortativity', 'Degree Variance', 'Two-Paths')
  rownames(correlation_matrix) <- all_variables_labels
  colnames(correlation_matrix) <- all_variables_labels
  
  pdf(matrix_file, width = 8, height = 8)
  corrplot(correlation_matrix, method = "square", type = "lower", diag = FALSE,
           tl.col = "black", tl.srt = 45, addCoef.col = "grey40", cl.ratio = 0.225)
  dev.off()
} else {
  print('matrix_file already exists')
}

# Functions

combined_formula_generator <- function(dependent_var, all_vars) {
  formulas <- list()
  
  # All combinations of variables
  comb <- unlist(lapply(1:length(all_vars), function(k) combn(all_vars, k, simplify = FALSE)), recursive = FALSE)
  
  for (vars in comb) {
    # formulas with linear variables
    linear_comb <- vars
    formula_str <- paste(dependent_var, "~ I(factor(nettype)) +", paste(linear_comb, collapse = " + "))
    formulas <- c(formulas, as.formula(formula_str))
    
    # formulas with quadratic variables
    for (i in seq_along(vars)) {
      quadratic_comb <- vars
      quadratic_comb[i] <- paste0("I(", vars[i], "^2)")
      formula_str <- paste(dependent_var, "~ I(factor(nettype)) +", paste(quadratic_comb, collapse = " + "))
      formulas <- c(formulas, as.formula(formula_str))
    }
    
    # formulas with logaritmic variables
    vars_for_log <- vars
    # Identify vars that can be negative
    negative_vars <- c("igraph_assortativity", "igraph_modularity")
    
    for (i in seq_along(vars)) {
      if (vars[i] %in% negative_vars) next # Skip log for negative vars
      
      logarithmic_comb <- vars
      logarithmic_comb[i] <- paste0("log(", vars[i], ")")
      formula_str <- paste(dependent_var, "~ I(factor(nettype)) +", paste(logarithmic_comb, collapse = " + "))
      formulas <- c(formulas, as.formula(formula_str))
    }
  }
  
  return(formulas)
}

run_combined_regression_models <- function(dependent_var, regression_method, output_file) {
  
  formulas_combined <- combined_formula_generator(dependent_var, all_variables)
  message(paste("Total models to run:", length(formulas_combined)))
  
  run_model <- function(formula, method) {
    tryCatch({
      model <- do.call(method, list(formula, data = simstats))
      
      if (any(is.na(coef(model)))) return(NULL)
      
      vif_val <- try(car::vif(model), silent = TRUE)
      if (!inherits(vif_val, "try-error")) {
        if (is.matrix(vif_val)) {
           # Use GVIF^(1/(2*Df)) which corrects for degrees of freedom
           if ("GVIF^(1/(2*Df))" %in% colnames(vif_val)) {
             max_vif <- max(vif_val[, "GVIF^(1/(2*Df))"])
           } else {
             max_vif <- max(vif_val[, ncol(vif_val)])
           }
        } else {
           max_vif <- max(vif_val)
        }
        
        if (max_vif > 10) return(NULL)
      }
      
      list(
        formula = formula,
        coeff = coef(model),
        aic = AIC(model),
        bic = BIC(model)
      )
    }, error = function(e) return(NULL))
  }
  
  start_time <- Sys.time()
  packages <- if (regression_method == "glm.nb") c("MASS", "car") else c("car")
  
  cluster <- makeForkCluster(8)
  registerDoParallel(cluster)
  
  clusterEvalQ(cluster, {
    library(car)
    library(MASS)
  })
  
  models <- clusterMap(cluster, run_model, formulas_combined, MoreArgs = list(method = regression_method), SIMPLIFY = FALSE)
  
  stopCluster(cluster)
  
  models <- Filter(Negate(is.null), models)
  message(paste("Models kept after filtering:", length(models)))
  
  saveRDS(models, file = output_file)
  
  end_time <- Sys.time()
  time_taken_parall <- end_time - start_time
  print(paste("Time taken for parallel execution: ", time_taken_parall))
  
  return(models)
}

Q_best_models <- function(Q_values, all_models, all_variables_combined, porcent) {
  
  aic_values <- unlist(lapply(all_models, function(x) x$aic))
  bic_values <- unlist(lapply(all_models, function(x) x$bic))
  
  formulas_all_models <- as.character(unlist(lapply(all_models, function(x) x$formula)))
  
  all_models_df <- as.data.frame(list(
    aic = aic_values,
    bic = bic_values,
    formula = formulas_all_models
  ))
  
  mean_AIC <- mean(aic_values)
  mean_BIC <- mean(bic_values)
  sd_AIC <- sd(aic_values)
  sd_BIC <- sd(bic_values)
  
  all_models_df$AIC_z <- (all_models_df$aic - mean_AIC) / sd_AIC
  all_models_df$BIC_z <- (all_models_df$bic - mean_BIC) / sd_BIC
  
  all_models_df$Mean_z <- rowMeans(all_models_df[, c("AIC_z", "BIC_z")])
  all_models_df <- all_models_df[order(all_models_df$Mean_z), ]
  
  variable_freq <- data.frame(Variable = all_variables_combined)
  top_Q_formulas <- list()
  variable_count_freq_list <- list()
  
  for (Q in Q_values) {
    if (porcent) {
      num_models <- floor((Q / 100) * nrow(all_models_df)) 
    } else {
      num_models <- Q
    }
    
    top_Q_formulas[[paste0("top_Q_formulas_Q", Q)]] <- all_models_df$formula[1:num_models]
    
    freq <- vector("numeric", length(all_variables_combined))
    names(freq) <- all_variables_combined
    
    for (i in seq_along(all_variables_combined)) {
      var <- all_variables_combined[i]
      count <- 0
      for (term in top_Q_formulas[[paste0("top_Q_formulas_Q", Q)]]) {
        if (var %in% labels(terms(as.formula(term)))) {count <- count + 1}
      }
      freq[i] <- count
    }
    variable_freq[[paste0("Q_", Q)]] <- freq / num_models
    
    formulas <- top_Q_formulas[[paste0("top_Q_formulas_Q", Q)]]
    
    variable_counts <- sapply(formulas, function(formula) {
      variables <- all.vars(as.formula(formula))[-1]
      variables <- variables[variables != "nettype"]
      return(length(variables))
    }) |> table()
    
    variable_count_freq_list[[paste0("Q_", Q)]] <- variable_counts / num_models
  }
  
  variable_count_freq_df <- do.call(rbind, lapply(names(variable_count_freq_list), function(name) {
    df <- as.data.frame(variable_count_freq_list[[name]])
    df$Q <- name
    df }))
  variable_count_freq_df <- variable_count_freq_df %>% rename(Num_var = Var1)
  
  return(list(variable_freq = variable_freq, variable_count_freq_df = variable_count_freq_df))
}

graphic_generator <- function(variable_count_freq_df, var_freq_long, all_variables_combined, main_title) {
  
  var_freq_sum <- aggregate(Frequency ~ Variable, data = var_freq_long, sum)
  max_main <- max(var_freq_sum$Frequency)
  
  ordered_vars <- var_freq_sum$Variable[order(-var_freq_sum$Frequency)]
  
  rename_map <- data.frame(
    original = all_variables_combined,
    aesthetic = all_variables_combined 
  )
  
  rename_map$aesthetic <- gsub("igraph_", "", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("ergm_", "", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("avg_", "Avg. ", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("local_", "Local ", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("degree_variance", "Variance", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("twopath", "Two-Paths", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("assortativity", "Assortativity", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("path_length", "Path Length", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("transitivity", "Transitivity", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("modularity", "Modularity", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("triangles", "Triangles", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("balance", "Balance", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("I\\((.*)\\^2\\)", "[\\1^2]", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("log\\((.*)\\)", "Log(\\1)", rename_map$aesthetic)
  
  var_freq_long$Variable <- factor(var_freq_long$Variable, levels = ordered_vars)
  levels(var_freq_long$Variable) <- rename_map$aesthetic[match(levels(var_freq_long$Variable), rename_map$original)]
  
  small_plot <- ggplot(variable_count_freq_df, aes(x = Num_var, y = Freq, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                      labels = c("Q_50" = "50%", "Q_25" = "25%", "Q_10" = "10%")) +
    scale_y_continuous(breaks = seq(0, 2.50, by = 0.50)) + 
    scale_x_discrete(limits = as.character(1:5)) +
    theme_minimal() +
    labs(x = "N° Variables", y = NULL) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank())
  
  main_plot <- ggplot(var_freq_long, aes(x = Variable, y = Frequency, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                      labels = c("Q_50" = "50%", "Q_25" = "25%", "Q_10" = "10%")) +
    scale_y_continuous(breaks = seq(0, 2.00, by = 0.25)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          legend.position = c(0.92, 0.75)) + 
    labs(title = main_title, x = NULL, y = "Sum of Frequencies", fill = "Q") + 
    annotation_custom(ggplotGrob(small_plot), xmin = 8.5, xmax = 14.0, ymin = max_main*0.35, ymax = max_main*1.02)
  
  print(main_plot)
  
  main_title_clean <- gsub(" ", "_", main_title)
  ggsave(file.path('03-model-comparison-combined-new-var-plus-twopaths', paste0(main_title_clean, ".pdf")), main_plot, width = 5, height = 3.4)
}

################################## Peak Preval Models ----------------------------

peak_preval_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_peak_preval.RData'

if (!file.exists(peak_preval_combined_file)) {
  print('Running peak_preval_combined_models')
  peak_preval_combined_models <- run_combined_regression_models('peak_preval', 'glm', peak_preval_combined_file)
} else {
  print('Loading peak_preval_combined_models')
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
}

peak_preval_combined_models <- Filter(Negate(is.null), peak_preval_combined_models)

Q_values <- c(50, 25, 10)
porcent <- TRUE
Q_best_peak_preval_combined <- Q_best_models(Q_values, peak_preval_combined_models, all_variables_combined, porcent)
var_freq_peak_preval_combined <- Q_best_peak_preval_combined$variable_freq
var_count_freq_df_peak_preval_combined <- Q_best_peak_preval_combined$variable_count_freq_df

var_freq_long_peak_preval_combined <- melt(var_freq_peak_preval_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_preval_combined, var_freq_long_peak_preval_combined, all_variables_combined, "Peak Prevalence")

################################## Peak Time Models ----------------------------

peak_time_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_peak_time.RData'

if (!file.exists(peak_time_combined_file)) {
  print('Running peak_time_combined_models')
  peak_time_combined_models <- run_combined_regression_models('peak_time', 'glm', peak_time_combined_file)
} else {
  print('Loading peak_time_combined_models')
  peak_time_combined_models <- readRDS(peak_time_combined_file)
}

peak_time_combined_models <- Filter(Negate(is.null), peak_time_combined_models)

Q_values <- c(50, 25, 10)
porcent <- TRUE
Q_best_peak_time_combined <- Q_best_models(Q_values, peak_time_combined_models, all_variables_combined, porcent)
var_freq_peak_time_combined <- Q_best_peak_time_combined$variable_freq
var_count_freq_df_peak_time_combined <- Q_best_peak_time_combined$variable_count_freq_df

var_freq_long_peak_time_combined <- melt(var_freq_peak_time_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_time_combined, var_freq_long_peak_time_combined, all_variables_combined, "Peak Time")

################################## Gentime Models ------------------------------

gentime_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_gentime.RData'

if (!file.exists(gentime_combined_file)) {
  print('Running gentime_combined_models')
  gentime_combined_models <- run_combined_regression_models('gentime', 'glm', gentime_combined_file)
} else {
  print('Loading gentime_combined_models')
  gentime_combined_models <- readRDS(gentime_combined_file)
}

gentime_combined_models <- Filter(Negate(is.null), gentime_combined_models)

Q_values <- c(50, 25, 10)
porcent <- TRUE
Q_best_gentime_combined <- Q_best_models(Q_values, gentime_combined_models, all_variables_combined, porcent)
var_freq_gentime_combined <- Q_best_gentime_combined$variable_freq
var_count_freq_df_gentime_combined <- Q_best_gentime_combined$variable_count_freq_df

var_freq_long_gentime_combined <- melt(var_freq_gentime_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_gentime_combined, var_freq_long_gentime_combined, all_variables_combined, 'Generation time')

################################## Rep Num Models ------------------------------

rep_num_combined_file <- '03-model-comparison-combined-new-var-plus-twopaths/03-combined_models_rep_num.RData'

if (!file.exists(rep_num_combined_file)) {
  print('Running rep_num_combined_models')
  rep_num_combined_models <- run_combined_regression_models('rt', 'glm.nb', rep_num_combined_file)
} else {
  print('Loading rep_num_combined_models')
  rep_num_combined_models <- readRDS(rep_num_combined_file)
}

rep_num_combined_models <- Filter(Negate(is.null), rep_num_combined_models)

Q_values <- c(50, 25, 10)
porcent <- TRUE
Q_best_rep_num_combined <- Q_best_models(Q_values, rep_num_combined_models, all_variables_combined, porcent)
var_freq_rep_num_combined <- Q_best_rep_num_combined$variable_freq
var_count_freq_df_rep_num_combined <- Q_best_rep_num_combined$variable_count_freq_df

var_freq_long_rep_num_combined <- melt(var_freq_rep_num_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_rep_num_combined, var_freq_long_rep_num_combined, all_variables_combined, 'Reproductive number')

# End Global Timer
global_end_time <- Sys.time()
message(paste("Analysis finished at:", global_end_time))
message(paste("Total duration:", round(difftime(global_end_time, global_start_time, units = "mins"), 2), "minutes"))
