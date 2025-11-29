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

# The 10 Selected Variables
all_variables <- c(
  "igraph_avg_degree",
  "igraph_avg_path_length",
  "igraph_local_transitivity",
  "igraph_modularity",
  "ergm_balance",
  "ergm_triangles",
  "igraph_degree_variance",
  "igraph_assortativity",
  "igraph_spectral_radius",
  "igraph_avg_betweenness"
)

all_variables_combined <- c(
  "igraph_avg_degree", "I(igraph_avg_degree^2)", "log(igraph_avg_degree)",
  "igraph_avg_path_length", "I(igraph_avg_path_length^2)", "log(igraph_avg_path_length)",
  "igraph_local_transitivity", "I(igraph_local_transitivity^2)", "log(igraph_local_transitivity)",
  "igraph_modularity", "I(igraph_modularity^2)", "log(igraph_modularity)",
  "ergm_balance", "I(ergm_balance^2)", "log(ergm_balance)",
  "ergm_triangles", "I(ergm_triangles^2)", "log(ergm_triangles)",
  "igraph_degree_variance", "I(igraph_degree_variance^2)", "log(igraph_degree_variance)",
  "igraph_assortativity", "I(igraph_assortativity^2)", "log(igraph_assortativity)",
  "igraph_spectral_radius", "I(igraph_spectral_radius^2)", "log(igraph_spectral_radius)",
  "igraph_avg_betweenness", "I(igraph_avg_betweenness^2)", "log(igraph_avg_betweenness)"
)

# Let's see the factors betwwen the minimum and maximum
variables_metadata <- data.frame(t(
  sapply(all_variables, function(var) {
    # Handle potential negative values for log
    # Shift if necessary for log? The formula generator just does log(x). 
    # If x <= 0, log(x) is NaN. 
    # Assortativity can be negative. Modularity can be negative.
    # We should probably handle this in the formula generator or data prep.
    # For now, let's just compute stats.
    max_val <- max(simstats[[var]], na.rm = TRUE)
    min_val <- min(simstats[[var]], na.rm = TRUE)
    mean_val <- mean(simstats[[var]], na.rm = TRUE)
    sd_val <- sd(simstats[[var]], na.rm = TRUE)
    c(Mean = mean_val, SD = sd_val, Factor = max_val / min_val)})
))

# Create output directory
dir.create("03-model-comparison-combined-2", showWarnings = FALSE)

# correlation matrix
matrix_file <- '03-model-comparison-combined-2/03-matrix_correlation.pdf'

if (!file.exists(matrix_file)) {
  correlation_matrix <- cor(simstats[, ..all_variables], use = "complete.obs")
  all_variables_labels <- c('Avg. Degree', 'Avg. Path Length', 'Local Transitivity', 'Modularity', 
                            'Balance', 'Triangles', 'Variance', 'Assortativity', 
                            'Spectral Radius', 'Avg. Betweenness')
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
  # Limit combination size to avoid explosion if needed, but user asked for ALL.
  # With 10 vars, 1023 combinations.
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
    # Check if variables are positive before adding log
    # We will add them, but the model might fail if data is <= 0.
    # We should probably filter out log models for variables that have <= 0 values.
    # Assortativity and Modularity can be negative.
    # We will skip log for them in the loop if they are in the combination.
    
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
      
      # Check for aliased coefficients (perfect collinearity)
      if (any(is.na(coef(model)))) return(NULL)
      
      # Check GVIF
      # If GVIF > 10, discard
      # Note: GVIF is for generalized linear models.
      # For simple models we can use VIF.
      # car::vif works for glm and lm.
      # However, calculating VIF for every model is expensive.
      # The user mentioned "if the value was above 10, the model was dropped".
      
      vif_val <- try(car::vif(model), silent = TRUE)
      if (!inherits(vif_val, "try-error")) {
        # Handle GVIF output (it can be a matrix or vector)
        if (is.matrix(vif_val)) {
           # For GVIF, we look at GVIF^(1/(2*Df))
           # If any is > sqrt(10) ~ 3.16? Or just GVIF > 10?
           # User said "if the value was above 10".
           # Usually for GVIF it's the adjusted one.
           # Let's assume standard VIF > 10 check.
           # If matrix (GVIF), take the max of GVIF column
           max_vif <- max(vif_val[,1])
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
  
  # Use 8 cores for M4 Pro
  cluster <- makeForkCluster(8)
  registerDoParallel(cluster)
  
  # Export necessary libraries to workers
  clusterEvalQ(cluster, {
    library(car)
    library(MASS)
  })
  
  models <- clusterMap(cluster, run_model, formulas_combined, MoreArgs = list(method = regression_method), SIMPLIFY = FALSE)
  
  stopCluster(cluster)
  
  # Filter NULLs
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
    
    # Counting the frequencies of variables #
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
    
    # Counting the number of variables in Q% best models #
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
  
  # Calcular la suma de 'Frequency' para cada variable
  var_freq_sum <- aggregate(Frequency ~ Variable, data = var_freq_long, sum)
  max_main <- max(var_freq_sum$Frequency)
  
  # Ordenar las variables por la suma de 'Frequency' en orden descendente
  ordered_vars <- var_freq_sum$Variable[order(-var_freq_sum$Frequency)]
  
  # Renombrar las variables para el gráfico
  # Need to update this map for the new variables
  # This is purely aesthetic mapping
  
  # Create a dynamic map based on what we have
  rename_map <- data.frame(
    original = all_variables_combined,
    aesthetic = all_variables_combined # Default to original
  )
  
  # Custom mappings
  rename_map$aesthetic <- gsub("igraph_", "", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("ergm_", "", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("avg_", "Avg. ", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("local_", "Local ", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("degree_variance", "Variance", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("spectral_radius", "Spectral Radius", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("assortativity", "Assortativity", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("betweenness", "Betweenness", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("path_length", "Path Length", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("transitivity", "Transitivity", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("modularity", "Modularity", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("balance", "Balance", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("triangles", "Triangles", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("I\\((.*)\\^2\\)", "[\\1^2]", rename_map$aesthetic)
  rename_map$aesthetic <- gsub("log\\((.*)\\)", "Log(\\1)", rename_map$aesthetic)
  
  # Renombramos, manteniendo el orden
  var_freq_long$Variable <- factor(var_freq_long$Variable, levels = ordered_vars)
  levels(var_freq_long$Variable) <- rename_map$aesthetic[match(levels(var_freq_long$Variable), rename_map$original)]
  
  # Gráfico pequeño
  small_plot <- ggplot(variable_count_freq_df, aes(x = Num_var, y = Freq, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                      labels = c("Q_50" = "50%", "Q_25" = "25%", "Q_10" = "10%")) +
    scale_y_continuous(breaks = seq(0, 2.50, by = 0.75)) +  # Marcas en el eje-y
    #scale_x_discrete(limits = as.character(1:10)) +  # Marcas en el eje-x
    theme_minimal() +
    labs(x = "N° Variables", y = NULL) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank())
  
  # Gráfico principal
  main_plot <- ggplot(var_freq_long, aes(x = Variable, y = Frequency, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                      labels = c("Q_50" = "50%", "Q_25" = "25%", "Q_10" = "10%")) +
    scale_y_continuous(breaks = seq(0, 2.00, by = 0.2)) +  # Marcas en el eje-y
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          legend.position = c(0.92, 0.75)) + # Posicionar la leyenda dentro del gráfico
    labs(title = main_title, x = NULL, y = "Sum of Frequencies", fill = "Q") + 
    annotation_custom(ggplotGrob(small_plot), xmin = length(ordered_vars)*0.6, xmax = length(ordered_vars), ymin = max_main*0.35, ymax = max_main*1.02)
  
  print(main_plot)
  
  main_title_clean <- gsub(" ", "_", main_title)
  ggsave(file.path('03-model-comparison-combined-2', paste0(main_title_clean, ".pdf")), main_plot, width = 10, height = 6)
}

################################## Peak Preval Models ----------------------------

peak_preval_combined_file <- '03-model-comparison-combined-2/03-combined_models_peak_preval.RData'

if (!file.exists(peak_preval_combined_file)) {
  print('Running peak_preval_combined_models')
  peak_preval_combined_models <- run_combined_regression_models('peak_preval', 'glm', peak_preval_combined_file)
} else {
  print('Loading peak_preval_combined_models')
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
}

peak_preval_combined_models <- Filter(Negate(is.null), peak_preval_combined_models)

# Figure
Q_values <- c(50, 25, 10)
porcent <- TRUE
Q_best_peak_preval_combined <- Q_best_models(Q_values, peak_preval_combined_models, all_variables_combined, porcent)
var_freq_peak_preval_combined <- Q_best_peak_preval_combined$variable_freq
var_count_freq_df_peak_preval_combined <- Q_best_peak_preval_combined$variable_count_freq_df

var_freq_long_peak_preval_combined <- melt(var_freq_peak_preval_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_preval_combined, var_freq_long_peak_preval_combined, all_variables_combined, "Peak Prevalence")

################################## Peak Time Models ----------------------------

peak_time_combined_file <- '03-model-comparison-combined-2/03-combined_models_peak_time.RData'

if (!file.exists(peak_time_combined_file)) {
  print('Running peak_time_combined_models')
  peak_time_combined_models <- run_combined_regression_models('peak_time', 'glm', peak_time_combined_file)
} else {
  print('Loading peak_time_combined_models')
  peak_time_combined_models <- readRDS(peak_time_combined_file)
}

peak_time_combined_models <- Filter(Negate(is.null), peak_time_combined_models)

# Figure
Q_values <- c(50, 25, 10)
porcent <- TRUE
Q_best_peak_time_combined <- Q_best_models(Q_values, peak_time_combined_models, all_variables_combined, porcent)
var_freq_peak_time_combined <- Q_best_peak_time_combined$variable_freq
var_count_freq_df_peak_time_combined <- Q_best_peak_time_combined$variable_count_freq_df

var_freq_long_peak_time_combined <- melt(var_freq_peak_time_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_time_combined, var_freq_long_peak_time_combined, all_variables_combined, "Peak Time")

################################## Gentime Models ------------------------------

gentime_combined_file <- '03-model-comparison-combined-2/03-combined_models_gentime.RData'

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

rep_num_combined_file <- '03-model-comparison-combined-2/03-combined_models_rep_num.RData'

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
