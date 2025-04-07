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

################################## Data & Functions ----------------------------

simstats <- fread("data/02-dataprep-network-stats-new-exposed.csv.gz")
simstats <- simstats[peak_preval > 1]

# Rescaling variables
n <- 534
twopath_complete <- (n * (n - 1) * (n - 2)) / 2
balance_complete <- ergm_triangle <- (n * (n - 1) * (n - 2)) / 6
triangle_complete <- balance_complete

simstats$ergm_twopath <- (simstats$ergm_twopath / twopath_complete) * 100
simstats$ergm_balance <- (simstats$ergm_balance / balance_complete) * 100
simstats$ergm_triangle <- (simstats$ergm_triangle / triangle_complete) * 100
simstats$igraph_density <- simstats$igraph_density * 100

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

# Variables (without density)
all_variables <- c("ergm_balance", "ergm_triangle", "ergm_twopath", "igraph_avg_path_length",
                   "igraph_avg_degree", "igraph_transitivity", "igraph_modularity")

all_variables_combined <- c(
  "ergm_balance", "I(ergm_balance^2)", "log(ergm_balance)",
  "ergm_triangle", "I(ergm_triangle^2)", "log(ergm_triangle)",
  "ergm_twopath", "I(ergm_twopath^2)", "log(ergm_twopath)",
  "igraph_avg_path_length", "I(igraph_avg_path_length^2)", "log(igraph_avg_path_length)",
  "igraph_avg_degree", "I(igraph_avg_degree^2)", "log(igraph_avg_degree)",
  "igraph_transitivity", "I(igraph_transitivity^2)", "log(igraph_transitivity)",
  "igraph_modularity", "I(igraph_modularity^2)", "log(igraph_modularity)"
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

# correlation matrix

matrix_file <- 'data/03-model-comparison-combined/03-matrix_correlation.pdf'

if (!file.exists(matrix_file)) {
  correlation_matrix <- cor(simstats[, ..all_variables], use = "complete.obs")
  all_variables_labels <- c('Balance', 'Triangle', 'Two paths', 'Avg. path length', 'Avg. degree', 'Transitivity', 'Modularity')
  rownames(correlation_matrix) <- all_variables_labels
  colnames(correlation_matrix) <- all_variables_labels
  
  pdf(matrix_file, width = 5.5, height = 5.5)
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
    for (i in seq_along(vars)) {
      logarithmic_comb <- vars
      logarithmic_comb[i] <- paste0("log(", vars[i], ")")
      formula_str <- paste(dependent_var, "~ I(factor(nettype)) +", paste(logarithmic_comb, collapse = " + "))
      formulas <- c(formulas, as.formula(formula_str))
    }
  }
  
  return(formulas)
}

run_combined_regression_models <- function(dependent_var, regression_method, output_file, combinations) {
  
  formulas_combined <- combined_formula_generator(dependent_var, all_variables)
  
  run_model <- function(formula, method) {
    model <- do.call(method, list(formula, data = simstats))
    
    coefficients <- model$coefficients
    if (any(is.na(coefficients))) {
      return(NULL)
    }
    
    list(
      formula = formula,
      coeff = coefficients,
      aic = AIC(model),
      bic = BIC(model)
    )
  }
  
  start_time <- Sys.time()
  packages <- if (regression_method == "glm.nb") c("MASS") else character(0)
  
  cluster <- makeForkCluster(20)
  registerDoParallel(cluster)
  
  models <- clusterMap(cluster, run_model, formulas_combined, MoreArgs = list(method = regression_method), SIMPLIFY = FALSE)
  
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
  all_variables_aesthetic <- c(
    'Balance', '[Balance^2]', 'Log(Balance)',
    'Triangle', '[Triangle^2]', 'Log(Triangle)',
    'Two path', '[Two path^2]', 'Log(Two path)',
    'Avg. path length', '[Avg. path length^2]', 'Log(Avg. path length)',
    'Avg. degree', '[Avg. degree^2]', 'Log(Avg. degree)',
    'Transitivity', '[Transitivity^2]', 'Log(Transitivity)',
    'Modularity', '[Modularity^2]', 'Log(Modularity)'
  )
  
  rename_map <- data.frame(
    original = all_variables_combined,
    aesthetic = all_variables_aesthetic
  )
  
  # Renombramos, manteniendo el orden
  var_freq_long$Variable <- factor(var_freq_long$Variable, levels = ordered_vars)
  levels(var_freq_long$Variable) <- rename_map$aesthetic[match(levels(var_freq_long$Variable), rename_map$original)]
  
  # Gráfico pequeño
  small_plot <- ggplot(variable_count_freq_df, aes(x = Num_var, y = Freq, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    #scale_fill_manual(values = c("Q_30" = "deepskyblue3", "Q_20" = "chartreuse3", "Q_10" = "red"),
    #                  labels = c("Q_30" = "30%", "Q_20" = "20%", "Q_10" = "10%")) +
    scale_fill_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                      labels = c("Q_50" = "50%", "Q_25" = "25%", "Q_10" = "10%")) +
    scale_y_continuous(breaks = seq(0, 2.50, by = 0.75)) +  # Marcas en el eje-y
    scale_x_discrete(limits = as.character(1:5)) +  # Marcas en el eje-x
    theme_minimal() +
    labs(x = "N° Variables", y = NULL) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank())  # Eliminar la cuadrícula mayor
          #panel.grid.minor = element_blank())  # Eliminar la cuadrícula menor
  
  # Gráfico principal
  main_plot <- ggplot(var_freq_long, aes(x = Variable, y = Frequency, fill = Q)) +
    geom_bar(stat = "identity", position = "stack") +
    #scale_fill_manual(values = c("Q_30" = "deepskyblue3", "Q_20" = "chartreuse3", "Q_10" = "red"),
    #                  labels = c("Q_30" = "30%", "Q_20" = "20%", "Q_10" = "10%")) +
    scale_fill_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                      labels = c("Q_50" = "50%", "Q_25" = "25%", "Q_10" = "10%")) +
    scale_y_continuous(breaks = seq(0, 2.00, by = 0.2)) +  # Marcas en el eje-y
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          legend.position = c(0.92, 0.75)) + # Posicionar la leyenda dentro del gráfico
    labs(title = main_title, x = NULL, y = "Sum of Frequencies", fill = "Q") + 
    annotation_custom(ggplotGrob(small_plot), xmin = 8.5, xmax = 14.0, ymin = max_main*0.35, ymax = max_main*1.02)
  
  print(main_plot)
  
  main_title_clean <- gsub(" ", "_", main_title)
  ggsave(file.path('data/03-model-comparison-combined', paste0(main_title_clean, ".pdf")), main_plot, width = 6, height = 4)
  #ggsave(file.path('data/03-model-comparison-combined.pdf'), main_plot, width = 6, height = 4)
}

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

peak_preval_combined_file <- 'data/03-model-comparison-combined/03-combined_models_peak_preval.RData'

if (!file.exists(peak_preval_combined_file)) {
  print('Running peak_preval_combined_models')
  run_combined_regression_models('peak_preval', 'glm', peak_preval_combined_file, all_vars)
} else {
  print('Loading peak_preval_combined_models')
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
}

peak_preval_combined_models <- Filter(Negate(is.null), peak_preval_combined_models) # 1023 originally, 32 dropped

# Figure

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_peak_preval_combined <- Q_best_models(Q_values, peak_preval_combined_models, all_variables_combined, porcent)
var_freq_peak_preval_combined <- Q_best_peak_preval_combined$variable_freq
var_count_freq_df_peak_preval_combined <- Q_best_peak_preval_combined$variable_count_freq_df

var_freq_long_peak_preval_combined <- melt(var_freq_peak_preval_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_preval_combined, var_freq_long_peak_preval_combined, all_variables_combined, "Peak Prevalence")

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_peak_preval <- Q_best_models_table(Q_value, peak_preval_combined_models, all_variables_combined, porcent)

Q_best_models_peak_preval <- lapply(Q_best_formulas_peak_preval, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Runing best models

aic_values_peak_preval <- round(sapply(Q_best_models_peak_preval, AIC), 3)
bic_values_peak_preval <- round(sapply(Q_best_models_peak_preval, BIC), 3)
best_variables_peak_preval <- c("Net(Scale-free)", "Net(SW) (p=0.1)", "Net(SW) (p=0.1)", "Net(Degree sequence)",
                                "Net(Erdös–Rényi)", "Balance", "Triangle", "Two path", "Log(Avg. path length)",
                                "sq(Triangle)", "Avg. path length", "sq(Avg. path length)", "Avg. degree", "Log(Avg. Degree)")

tex_file_peak_preval <- 'data/03-model-comparison-combined/03-best_models_table_peak_preval.tex'

stargazer(Q_best_models_peak_preval, type = "latex", out = tex_file_peak_preval,
          title = "Best 10 models for Peak Prevalence",
          label = "tab:best_models_peak_preval",
          dep.var.labels = "Peak Prevalence",
          covariate.labels = best_variables_peak_preval,
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_peak_preval),
                           c("BIC", bic_values_peak_preval)),
          no.space = TRUE,
          digits = 3)

################################## Peak Time Models ----------------------------

peak_time_combined_file <- 'data/03-model-comparison-combined/03-combined_models_peak_time.RData'

if (!file.exists(peak_time_combined_file)) {
  print('Running peak_time_combined_models')
  run_combined_regression_models('peak_time', 'glm', peak_time_combined_file, all_vars)
} else {
  print('Loading peak_time_combined_models')
  peak_time_combined_models <- readRDS(peak_time_combined_file)
}

peak_time_combined_models <- Filter(Negate(is.null), peak_time_combined_models) # 1023 originally, 32 dropped also

# Figure

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_peak_time_combined <- Q_best_models(Q_values, peak_time_combined_models, all_variables_combined, porcent)
var_freq_peak_time_combined <- Q_best_peak_time_combined$variable_freq
var_count_freq_df_peak_time_combined <- Q_best_peak_time_combined$variable_count_freq_df

var_freq_long_peak_time_combined <- melt(var_freq_peak_time_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_time_combined, var_freq_long_peak_time_combined, all_variables_combined, "Peak Time")

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_peak_time <- Q_best_models_table(Q_value, peak_time_combined_models, all_variables_combined, porcent)

Q_best_models_peak_time <- lapply(Q_best_formulas_peak_time, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Runing best models

aic_values_peak_time <- round(sapply(Q_best_models_peak_time, AIC), 3)
bic_values_peak_time <- round(sapply(Q_best_models_peak_time, BIC), 3)
best_variables_peak_time <- c("Net(Scale-free)", "Net(SW) (p=0.1)", "Net(SW) (p=0.1)", "Net(Degree sequence)",
                                "Net(Erdös–Rényi)", "Log(Balance)", "sq(Balance)", "Avg. path length", "Avg. degree",
                                "Log(Modularity)", "Modularity", "sq(Avg. degree)", "Triangle", "Log(Two path)", "Transitivity")

tex_file_peak_time <- 'data/03-model-comparison-combined/03-best_models_table_peak_time.tex'

stargazer(Q_best_models_peak_time, type = "latex", out = tex_file_peak_time,
          title = "Best 10 models for Peak Time",
          label = "tab:best_models_peak_time",
          dep.var.labels = "Peak Time",
          covariate.labels = best_variables_peak_time,
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_peak_time),
                           c("BIC", bic_values_peak_time)),
          no.space = TRUE,
          digits = 3)

################################## Gentime Models ------------------------------

gentime_combined_file <- 'data/03-model-comparison-combined/03-combined_models_gentime.RData'

if (!file.exists(gentime_combined_file)) {
  print('Running gentime_combined_models')
  run_combined_regression_models('gentime', 'glm', gentime_combined_file)
} else {
  print('Loading gentime_combined_models')
  gentime_combined_models <- readRDS(gentime_combined_file)
}

gentime_combined_models <- Filter(Negate(is.null), gentime_combined_models) # 1023 originally, 32 dropped also!

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_gentime_combined <- Q_best_models(Q_values, gentime_combined_models, all_variables_combined, porcent)
var_freq_gentime_combined <- Q_best_gentime_combined$variable_freq
var_count_freq_df_gentime_combined <- Q_best_gentime_combined$variable_count_freq_df

var_freq_long_gentime_combined <- melt(var_freq_gentime_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_gentime_combined, var_freq_long_gentime_combined, all_variables_combined, 'Generation time')

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_gentime <- Q_best_models_table(Q_value, gentime_combined_models, all_variables_combined, porcent)

Q_best_models_gentime <- lapply(Q_best_formulas_gentime, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Runing best models

aic_values_gentime <- round(sapply(Q_best_models_gentime, AIC), 3)
bic_values_gentime <- round(sapply(Q_best_models_gentime, BIC), 3)
best_variables_gentime <- c("Net(Scale-free)", "Net(SW) (p=0.1)", "Net(SW) (p=0.1)", "Net(Degree sequence)",
                              "Net(Erdös–Rényi)", "Balance", "sq(Avg. degree)", "sq(Modularity)", "Two path",
                              "sq(Two path)", "Avg. degree", "Triangle", "Log(Avg. degree)", "Log(Modularity)")

tex_file_gentime <- 'data/03-model-comparison-combined/03-best_models_table_gentime.tex'

stargazer(Q_best_models_gentime, type = "latex", out = tex_file_gentime,
          title = "Best 10 models for Generation Time",
          label = "tab:best_models_gentime",
          dep.var.labels = "Peak Time",
          covariate.labels = best_variables_gentime,
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_gentime),
                           c("BIC", bic_values_gentime)),
          no.space = TRUE,
          digits = 3)

################################## Rep Num Models ------------------------------

rep_num_combined_file <- 'data/03-model-comparison-combined/03-combined_models_rep_num.RData'

if (!file.exists(rep_num_combined_file)) {
  print('Running rep_num_combined_models')
  run_combined_regression_models('rt', 'glm.nb', rep_num_combined_file)
} else {
  print('Loading rep_num_combined_models')
  rep_num_combined_models <- readRDS(rep_num_combined_file)
}

rep_num_combined_models <- Filter(Negate(is.null), rep_num_combined_models) # 1023 originally, 32 dropped also!!

Q_values <- c(30, 20, 10)
porcent <- TRUE
Q_best_rep_num_combined <- Q_best_models(Q_values, rep_num_combined_models, all_variables_combined, porcent)
var_freq_rep_num_combined <- Q_best_rep_num_combined$variable_freq
var_count_freq_df_rep_num_combined <- Q_best_rep_num_combined$variable_count_freq_df

var_freq_long_rep_num_combined <- melt(var_freq_rep_num_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_rep_num_combined, var_freq_long_rep_num_combined, all_variables_combined, 'Reproductive number')

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_rep_num <- Q_best_models_table(Q_value, rep_num_combined_models, all_variables_combined, porcent)

Q_best_models_rep_num <- lapply(Q_best_formulas_rep_num, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Runing best models

aic_values_rep_num <- round(sapply(Q_best_models_rep_num, AIC), 3)
bic_values_rep_num <- round(sapply(Q_best_models_rep_num, BIC), 3)
best_variables_rep_num <- c("Net(Scale-free)", "Net(SW) (p=0.1)", "Net(SW) (p=0.1)", "Net(Degree sequence)",
                            "Net(Erdös–Rényi)", "sq(Av. path length)", "Av. path length", "Log(Av. path length)", "sq(Balance)",
                            "sq(Avg. degree)", "Balance", "Avg. degree", "Log(Balance)", "Log(Avg. degree)", "Log(Two path)")

tex_file_rep_num <- 'data/03-model-comparison-combined/03-best_models_table_rep_num.tex'

stargazer(Q_best_models_rep_num, type = "latex", out = tex_file_rep_num,
          title = "Best 10 models for Generation Time",
          label = "tab:best_models_rep_num",
          dep.var.labels = "Peak Time",
          covariate.labels = best_variables_rep_num,
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_rep_num),
                           c("BIC", bic_values_rep_num)),
          no.space = TRUE,
          digits = 3)
