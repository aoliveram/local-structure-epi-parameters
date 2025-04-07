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


################################## Data & Functions ----------------------------

simstats <- fread("data/02-dataprep-network-stats-new-exposed.csv.gz")
simstats <- simstats[peak_preval > 1]

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

####

all_variables <- c("ergm_twopath", "ergm_balance", "ergm_triangle", "igraph_avg_degree", "igraph_avg_path_length", "igraph_density", "igraph_transitivity", "igraph_modularity")

all_variables_combined <- c(
  "ergm_twopath", "I(ergm_twopath^2)", "log(ergm_twopath)",
  "ergm_balance", "I(ergm_balance^2)", "log(ergm_balance)",
  "ergm_triangle", "I(ergm_triangle^2)", "log(ergm_triangle)",
  "igraph_avg_degree", "I(igraph_avg_degree^2)", "log(igraph_avg_degree)",
  "igraph_avg_path_length", "I(igraph_avg_path_length^2)", "log(igraph_avg_path_length)"
)

combined_formula_generator <- function(dependent_var, comb) {
  formulas <- list()
  
  # Generar todas las combinaciones de las variables
  comb <- unlist(lapply(1:length(all_vars), function(k) combn(all_vars, k, simplify = FALSE)), recursive = FALSE)
  
  for (vars in comb) {
    # Crear una fórmula con cada variable como cuadrática o logarítmica
    for (i in seq_along(vars)) {
      quadratic_comb <- vars
      quadratic_comb[i] <- paste0("I(", vars[i], "^2)")
      formula_str <- paste(dependent_var, "~ I(factor(nettype)) +", paste(quadratic_comb, collapse = " + "))
      formulas <- c(formulas, as.formula(formula_str))
      
      logarithmic_comb <- vars
      logarithmic_comb[i] <- paste0("log(", vars[i], ")")
      formula_str <- paste(dependent_var, "~ I(factor(nettype)) +", paste(logarithmic_comb, collapse = " + "))
      formulas <- c(formulas, as.formula(formula_str))
    }
  }
  
  return(formulas)
}

all_vars <- c("ergm_twopath", "ergm_balance", "ergm_triangle")

combinations <- unlist(lapply(1:length(all_vars), function(k) combn(all_vars, k, simplify = FALSE)), recursive = FALSE)
formulas_combined <- combined_formula_generator(dependent_var, combinations)

dependent_var <- 'peak_time'
regression_method <- 'glm'

#formulas_from_function <- list()
models <- list()
for (i in seq_along(formulas_combined)) {
  formula <- formulas_combined[[i]]
  model <- do.call(regression_method, list(formula, data = simstats))
  res_i <- list(
    formula = formula,
    coeff = summary(model)$coefficients,
    aic = AIC(model),
    bic = BIC(model)
  )
  models <- append(models, list(res_i))
}

run_combined_regression_models <- function(dependent_var, regression_method, output_file, all_vars) {
  
  combinations <- unlist(lapply(1:length(all_vars), function(k) combn(all_vars, k, simplify = FALSE)), recursive = FALSE)
  
  start_time <- Sys.time()
  packages <- if (regression_method == "glm.nb") c("MASS") else character(0)
  
  #cores <- detectCores() - 1
  cluster <- makeForkCluster(20)
  # registerDoParallel(cluster)
  
  # Usar clusterMap para paralelizar las tareas
  models <- list()
  
  models <- clusterMap(cluster, function(comb) {
    tryCatch({
      
      comb <- unlist(lapply(1:length(all_vars), function(k) combn(all_vars, k, simplify = FALSE)), recursive = FALSE)
      
      formulas_combined <- combined_formula_generator(dependent_var, comb)
      
      local_model_list <- list()
      #formulas_from_function <- list()
      
      for (i in seq_along(formulas_combined)) {
        formula <- formulas_combined[[i]]
        
        #formulas_from_function <- append(formulas_from_function, formula)
        
        # if (regression_method == "glm") {
        #   model <- glm(formula, data = simstats)
        # } else {
        #   model <- glm.nb(formula, data = simstats)
        # }
        # model <- glm(formula, data = simstats)
        model <- do.call(regression_method, list(formula, data = simstats))
        
        res_i <- list(
          formula = formula,
          coeff = summary(model)$coefficients,
          aic = AIC(model),
          bic = BIC(model)
        )
        local_model_list <- append(local_model_list, list(res_i)) # list(model))
        #list_fomulas <- append(list_fomulas, length(formula))
      }
      
      #return(formulas_from_function)
      return(local_model_list)
    }, error = function(e) e)
  }, combinations, .scheduling = "dynamic")
  
  
  #save(models, file = output_file)
  saveRDS(models, file = output_file)
  
  stopCluster(cluster)
  
  end_time <- Sys.time()
  time_taken_parall <- end_time - start_time
  print(paste("Time taken for parallel execution: ", time_taken_parall))
  
  return(NULL)
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
  
  for (Q in Q_values) {
    if (porcent) {
      num_models <- floor((Q / 100) * nrow(all_models_df)) 
    } else {
      num_models <- Q
    }
    
    top_Q_formulas[[paste0("top_Q_formulas_Q", Q)]] <- all_models_df$formula[1:num_models]
    
    # Counting the frequencies of variables
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
  }
  
  # Counting the number of variables in Q% best models
  variable_count_freq_list <- list()
  for (Q in Q_values) {
    formulas <- top_Q_formulas[[paste0("top_Q_formulas_Q", Q)]]
    
    variable_counts <- sapply(formulas, function(formula) {
      variables <- all.vars(as.formula(formula))[-1]
      variables <- variables[variables != "nettype"]
      return(length(variables))
    }) |> table()
    
    variable_count_freq_list[[paste0("Q_", Q)]] <- variable_counts
  }
  
  variable_count_freq_df <- do.call(rbind, lapply(names(variable_count_freq_list), function(name) {
    df <- as.data.frame(variable_count_freq_list[[name]])
    df$Q <- name
    df }))
  variable_count_freq_df <- variable_count_freq_df %>% rename(Num_var = Var1)
  
  return(list(variable_freq = variable_freq, variable_count_freq_df = variable_count_freq_df))
}

graphic_generator <- function(variable_count_freq_df, var_freq_long) {
  # Obtener el nombre del dataframe pasado como argumento
  df_name <- deparse(substitute(variable_count_freq_df))
  
  # Extraer el sufijo del nombre del dataframe
  suffix <- sub("var_count_freq_df_", "", df_name)
  
  max_y_small <- max(variable_count_freq_df$Freq)
  
  # Crear el gráfico pequeño
  small_plot <- ggplot(variable_count_freq_df, aes(x = Num_var, y = Freq, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    #scale_y_continuous(limits = c(0, max_y_small)) +
    labs(x = "N° Variables", y = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Determinar el título del gráfico principal
  title <- paste("Frecuencia de Variables", suffix)
  
  # Crear el gráfico principal
  main_plot <- ggplot(var_freq_long, aes(x = Variable, y = Frequency, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    #scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = title, x = NULL, y = "Frecuencia", fill = "Q") +
    annotation_custom(ggplotGrob(small_plot), xmin = 5, xmax = 11, ymin = 0.65, ymax = 1.05)
  
  # Mostrar el gráfico principal
  print(main_plot)
}

################################## Peak Preval Models ----------------------------

peak_preval_combined_file <- 'data/03-model-comparison-combined/03-peak_preval_combined_models.RData'

if (!file.exists(peak_preval_combined_file)) {
  print('Running peak_preval_combined_models')
  run_combined_regression_models('peak_preval', 'glm', peak_preval_combined_file)
} else {
  print('Loading peak_preval_combined_models')
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
}

peak_preval_combined_models <- unlist(peak_preval_combined_models, recursive = FALSE)

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_peak_preval_combined <- Q_best_models(Q_values, peak_preval_combined_models, all_variables_combined, porcent)
var_freq_peak_preval_combined <- Q_best_peak_preval_combined$variable_freq
var_count_freq_df_peak_preval_combined <- Q_best_peak_preval_combined$variable_count_freq_df

var_freq_long_peak_preval_combined <- melt(var_freq_peak_preval_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_preval_combined, var_freq_long_peak_preval_combined)

################################## Peak Time Models ----------------------------

all_vars <- c("ergm_twopath", "ergm_balance", "ergm_triangle")

all_variables_combined <- c(
  "ergm_twopath", "I(ergm_twopath^2)", "log(ergm_twopath)", 
  "ergm_balance", "I(ergm_balance^2)", "log(ergm_balance)",
  "ergm_triangle", "I(ergm_triangle^2)", "log(ergm_triangle)"
  )

peak_time_combined_file <- 'data/03-model-comparison-combined/03-peak_time_combined_models_test_2.RData'

if (!file.exists(peak_time_combined_file)) {
  print('Running peak_time_combined_models')
  run_combined_regression_models('peak_time', 'glm', peak_time_combined_file, all_vars)
} else {
  print('Loading peak_time_combined_models')
  peak_time_combined_models_2 <- readRDS(peak_time_combined_file)
}

peak_time_combined_models <- unlist(peak_time_combined_models, recursive = FALSE)

formulas_from_models <- list()
for (i in 1:24) {
  formulas_from_models_unique <- peak_time_combined_models[[i]]$formula
  formulas_from_models <- append(formulas_from_models, formulas_from_models_unique)
}



Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_peak_time_combined <- Q_best_models(Q_values, peak_time_combined_models, all_variables_combined, porcent)
var_freq_peak_time_combined <- Q_best_peak_time_combined$variable_freq
var_count_freq_df_peak_time_combined <- Q_best_peak_time_combined$variable_count_freq_df

var_freq_long_peak_time_combined <- melt(var_freq_peak_time_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_time_combined, var_freq_long_peak_time_combined)

################################## Gentime Models ------------------------------

gentime_combined_file <- 'data/03-model-comparison-combined/03-gentime_combined_models.RData'

if (!file.exists(gentime_combined_file)) {
  print('Running gentime_combined_models')
  run_combined_regression_models('gentime', 'glm', gentime_combined_file)
} else {
  print('Loading gentime_combined_models')
  gentime_combined_models <- readRDS(gentime_combined_file)
}

gentime_combined_models <- unlist(gentime_combined_models, recursive = FALSE)

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_gentime_combined <- Q_best_models(Q_values, gentime_combined_models, all_variables_combined, porcent)
var_freq_gentime_combined <- Q_best_gentime_combined$variable_freq
var_count_freq_df_gentime_combined <- Q_best_gentime_combined$variable_count_freq_df

var_freq_long_gentime_combined <- melt(var_freq_gentime_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_gentime_combined, var_freq_long_gentime_combined)

################################## Rep Num Models ------------------------------

rep_num_combined_file <- 'data/03-model-comparison-combined/03-rep_num_combined_models_6var.RData'

if (!file.exists(rep_num_combined_file)) {
  print('Running rep_num_combined_models')
  run_combined_regression_models('rt', 'glm.nb', rep_num_combined_file)
} else {
  print('Loading rep_num_combined_models')
  rep_num_combined_models_6var <- readRDS(rep_num_combined_file)
}

rep_num_combined_models_6var <- unlist(rep_num_combined_models_6var, recursive = FALSE)

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_rep_num_combined <- Q_best_models(Q_values, rep_num_combined_models, all_variables_combined, porcent)
var_freq_rep_num_combined <- Q_best_rep_num_combined$variable_freq
var_count_freq_df_rep_num_combined <- Q_best_rep_num_combined$variable_count_freq_df

var_freq_long_rep_num_combined <- melt(var_freq_rep_num_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_rep_num_combined, var_freq_long_rep_num_combined)
