library(car)
library(ggplot2)

run_combined_regression_models <- function(dependent_var, regression_method, output_file, combinations) {
  
  formulas_combined <- combined_formula_generator(dependent_var, all_variables)
  
  run_model <- function(formula, method) {
    model <- do.call(method, list(formula, data = simstats))
    
    coefficients <- model$coefficients
    if (any(is.na(coefficients))) {
      return(list(
        formula = formula,
        coeff = NULL,
        aic = NULL,
        bic = NULL,
        vif = NULL
      ))
    }
    
    # Calculate VIF values
    vif_values <- tryCatch({
      vif(model)
    }, error = function(e) {
      return(NA)  # Return NA if VIF calculation fails
    })
    
    list(
      formula = formula,
      coeff = coefficients,
      aic = AIC(model),
      bic = BIC(model),
      vif = vif_values
    )
  }
  
  start_time <- Sys.time()
  packages <- if (regression_method == "glm.nb") c("MASS") else character(0)
  
  cluster <- makeForkCluster(16)
  registerDoParallel(cluster)
  
  models <- clusterMap(cluster, run_model, formulas_combined, MoreArgs = list(method = regression_method), SIMPLIFY = FALSE)
  
  saveRDS(models, file = output_file)
  
  stopCluster(cluster)
  
  end_time <- Sys.time()
  time_taken_parall <- end_time - start_time
  print(paste("Time taken for parallel execution: ", time_taken_parall))
  
  return(NULL)
}

models_classifier <- function(models) {
  # Initialize a vector to store GVIF categories
  gvif_categories <- character(length(models))
  
  # Classify each model based on GVIF
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    if (!is.null(model$vif)) {
      gvif_values <- model$vif[, "GVIF^(1/(2*Df))"]
      max_gvif <- max(gvif_values, na.rm = TRUE)
      
      gvif_categories[i] <- if (max_gvif > 10) {
        "bad"
      } else if (max_gvif > 4) {
        "moderate"
      } else {
        "good"
      }
    } else {gvif_categories[i] <- "unknown"}
  }
  return(gvif_categories)
}

models_classifier_with_gvif <- function(models) {
  gvif_categories <- character(length(models))
  gvif_values <- numeric(length(models))
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    if (!is.null(model$vif)) {
      gvif_val <- model$vif[, "GVIF^(1/(2*Df))"]
      max_gvif <- max(gvif_val, na.rm = TRUE)
      gvif_values[i] <- max_gvif
      
      gvif_categories[i] <- if (max_gvif > 10) {
        "bad"
      } else if (max_gvif > 4) {
        "moderate"
      } else {
        "good"
      }
    } else {
      gvif_values[i] <- NA
      gvif_categories[i] <- "unknown"
    }
  }
  
  return(data.frame(Index = seq_along(models), GVIF_Adj = gvif_values, Category = gvif_categories))
}


################################## Peak Preval Models ----------------------------

peak_preval_combined_file <- 'data/03-model-comparison-combined/03-vif_combined_models_peak_preval.RData'

if (!file.exists(peak_preval_combined_file)) {
  print('Running peak_preval_combined_models')
  run_combined_regression_models('peak_preval', 'glm', peak_preval_combined_file, all_vars)
} else {
  print('Loading peak_preval_combined_models')
  peak_preval_combined_models <- readRDS(peak_preval_combined_file)
}

# Plot 1

# Classify models and obtain GVIF data
gvif_data_peak_preval <- models_classifier_with_gvif(peak_preval_combined_models)

# Count the number of models in each category
good_count_pp <- sum(gvif_data$Category == "good")
moderate_count_pp <- sum(gvif_data$Category == "moderate")
bad_count_pp <- sum(gvif_data$Category == "bad")

# Create the plot
ggplot(gvif_data_peak_preval) +
  geom_point(data = subset(gvif_data_peak_preval, Category == "bad"), aes(x = Index, y = GVIF_Adj), color = "red", size = 0.7) +
  geom_point(data = subset(gvif_data_peak_preval, Category == "moderate"), aes(x = Index, y = GVIF_Adj), color = "orange", size = 0.8) +
  geom_point(data = subset(gvif_data_peak_preval, Category == "good"), aes(x = Index, y = GVIF_Adj), color = "green", size = 0.8) +
  scale_color_manual(values = c("good" = "green", "moderate" = "orange", "bad" = "red", "unknown" = "gray")) +
  labs(title = "GVIF^(1/(2*Df)) for Models",
       x = "Model Index",
       y = "GVIF^(1/(2*Df))",
       caption = paste("Good:", good_count_pp, "| Moderate:", moderate_count_pp, "| Bad:", bad_count_pp)) +
  ylim(0, 75) +
  theme_minimal()

# Find the highest moderate index and its details
highest_moderate_index_pp <- max(gvif_data_peak_preval$Index[gvif_data_peak_preval$Category == "moderate"], na.rm = TRUE)
highest_moderate_index_pp
peak_preval_combined_models[[highest_moderate_index_pp]]$formula
gvif_data_peak_preval[highest_moderate_index_pp, 'GVIF_Adj']

# Plot 2

filtered_peak_preval_models <- peak_preval_combined_models[models_classifier(peak_preval_combined_models) %in% c("good", "moderate")]

# Run Q_best_models with the filtered models
#Q_values <- c(30, 20, 10); porcent <- FALSE
Q_values <- c(50, 25, 10); porcent <- TRUE
Q_best_peak_preval_combined <- Q_best_models(Q_values, filtered_peak_preval_models, all_variables_combined, porcent)

var_freq_peak_preval_combined <- Q_best_peak_preval_combined$variable_freq
var_count_freq_df_peak_preval_combined <- Q_best_peak_preval_combined$variable_count_freq_df

var_freq_long_peak_preval_combined <- melt(var_freq_peak_preval_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_preval_combined, var_freq_long_peak_preval_combined, all_variables_combined, "Peak Prevalence")

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_peak_preval <- Q_best_models_table(Q_value, filtered_peak_preval_models, all_variables_combined, porcent)
Q_best_models_peak_preval <- lapply(Q_best_formulas_peak_preval, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Runing best models

aic_values_peak_preval <- round(sapply(Q_best_models_peak_preval, AIC), 3)
bic_values_peak_preval <- round(sapply(Q_best_models_peak_preval, BIC), 3)
#best_variables_peak_preval <- c("Net(Scale-free)", "Net(SW) (p=0.1)", "Net(SW) (p=0.1)", "Net(Degree sequence)",
#                                "Net(Erdös–Rényi)", "Balance", "Triangle", "Two path", "Log(Avg. path length)",
#                                "sq(Triangle)", "Avg. path length", "sq(Avg. path length)", "Avg. degree", "Log(Avg. Degree)")

tex_file_peak_preval <- 'data/03-model-comparison-combined/03-best_models_table_peak_preval.tex'

stargazer(Q_best_models_peak_preval, type = "latex", out = tex_file_peak_preval,
          title = "Best 10 models for Peak Prevalence",
          label = "tab:best_models_peak_preval",
          dep.var.labels = "Peak Prevalence",
          #covariate.labels = best_variables_peak_preval,
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_peak_preval),
                           c("BIC", bic_values_peak_preval)),
          no.space = TRUE,
          digits = 3)


################################## Peak Time Models ----------------------------

peak_time_combined_file <- 'data/03-model-comparison-combined/03-vif_combined_models_peak_time.RData'

if (!file.exists(peak_time_combined_file)) {
  print('Running peak_time_combined_models')
  run_combined_regression_models('peak_time', 'glm', peak_time_combined_file, all_vars)
} else {
  print('Loading peak_time_combined_models')
  peak_time_combined_models <- readRDS(peak_time_combined_file)
}

# Plot 1

# Classify models and obtain GVIF data
gvif_data_peak_time <- models_classifier_with_gvif(peak_time_combined_models)

# Count the number of models in each category
good_count_pt <- sum(gvif_data_peak_time$Category == "good")
moderate_count_pt <- sum(gvif_data_peak_time$Category == "moderate")
bad_count_pt <- sum(gvif_data_peak_time$Category == "bad")

# Create the plot
ggplot(gvif_data_peak_time) +
  geom_point(data = subset(gvif_data_peak_time, Category == "bad"), aes(x = Index, y = GVIF_Adj), color = "red", size = 0.7) +
  geom_point(data = subset(gvif_data_peak_time, Category == "moderate"), aes(x = Index, y = GVIF_Adj), color = "orange", size = 0.8) +
  geom_point(data = subset(gvif_data_peak_time, Category == "good"), aes(x = Index, y = GVIF_Adj), color = "green", size = 0.8) +
  scale_color_manual(values = c("good" = "green", "moderate" = "orange", "bad" = "red", "unknown" = "gray")) +
  labs(title = "GVIF^(1/(2*Df)) for Models",
       x = "Model Index",
       y = "GVIF^(1/(2*Df))",
       caption = paste("Good:", good_count_pt, "| Moderate:", moderate_count_pt, "| Bad:", bad_count_pt)) +
  ylim(0, 75) +
  theme_minimal()

# Find the highest moderate index and its details
highest_moderate_index_pt <- max(gvif_data_peak_time$Index[gvif_data_peak_time$Category == "moderate"], na.rm = TRUE)
highest_moderate_index_pt
peak_time_combined_models[[highest_moderate_index_pt]]$formula
gvif_data_peak_time[highest_moderate_index_pt, 'GVIF_Adj']

# Plot 2

filtered_peak_time_models <- peak_time_combined_models[models_classifier(peak_time_combined_models) %in% c("good", "moderate")]

# Run Q_best_models with the filtered models
Q_values <- c(30, 20, 10); porcent <- FALSE
Q_values <- c(50, 25, 10); porcent <- TRUE
Q_best_peak_time_combined <- Q_best_models(Q_values, filtered_peak_time_models, all_variables_combined, porcent)

var_freq_peak_time_combined <- Q_best_peak_time_combined$variable_freq
var_count_freq_df_peak_time_combined <- Q_best_peak_time_combined$variable_count_freq_df

var_freq_long_peak_time_combined <- melt(var_freq_peak_time_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_peak_time_combined, var_freq_long_peak_time_combined, all_variables_combined, "Peak Time")

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_peak_time <- Q_best_models_table(Q_value, filtered_peak_time_models, all_variables_combined, porcent)
Q_best_models_peak_time <- lapply(Q_best_formulas_peak_time, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Running best models

aic_values_peak_time <- round(sapply(Q_best_models_peak_time, AIC), 3)
bic_values_peak_time <- round(sapply(Q_best_models_peak_time, BIC), 3)

tex_file_peak_time <- 'data/03-model-comparison-combined/03-best_models_table_peak_time.tex'

stargazer(Q_best_models_peak_time, type = "latex", out = tex_file_peak_time,
          title = "Best 10 models for Peak Time",
          label = "tab:best_models_peak_time",
          dep.var.labels = "Peak Time",
          #covariate.labels = best_variables_peak_time, # Uncomment and define if needed
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_peak_time),
                           c("BIC", bic_values_peak_time)),
          no.space = TRUE,
          digits = 3)

################################## Gentime Models ----------------------------

gentime_combined_file <- 'data/03-model-comparison-combined/03-vif_combined_models_gentime.RData'

if (!file.exists(gentime_combined_file)) {
  print('Running gentime_combined_models')
  run_combined_regression_models('gentime', 'glm', gentime_combined_file, all_vars)
} else {
  print('Loading gentime_combined_models')
  gentime_combined_models <- readRDS(gentime_combined_file)
}

# Plot 1

# Classify models and obtain GVIF data
gvif_data_gentime <- models_classifier_with_gvif(gentime_combined_models)

# Count the number of models in each category
good_count_gt <- sum(gvif_data_gentime$Category == "good")
moderate_count_gt <- sum(gvif_data_gentime$Category == "moderate")
bad_count_gt <- sum(gvif_data_gentime$Category == "bad")

# Create the plot
ggplot(gvif_data_gentime) +
  geom_point(data = subset(gvif_data_gentime, Category == "bad"), aes(x = Index, y = GVIF_Adj), color = "red", size = 0.7) +
  geom_point(data = subset(gvif_data_gentime, Category == "moderate"), aes(x = Index, y = GVIF_Adj), color = "orange", size = 0.8) +
  geom_point(data = subset(gvif_data_gentime, Category == "good"), aes(x = Index, y = GVIF_Adj), color = "green", size = 0.8) +
  scale_color_manual(values = c("good" = "green", "moderate" = "orange", "bad" = "red", "unknown" = "gray")) +
  labs(title = "GVIF^(1/(2*Df)) for Models",
       x = "Model Index",
       y = "GVIF^(1/(2*Df))",
       caption = paste("Good:", good_count_gt, "| Moderate:", moderate_count_gt, "| Bad:", bad_count_gt)) +
  ylim(0, 75) +
  theme_minimal()

# Find the highest moderate index and its details
highest_moderate_index_gt <- max(gvif_data_gentime$Index[gvif_data_gentime$Category == "moderate"], na.rm = TRUE)
highest_moderate_index_gt
gentime_combined_models[[highest_moderate_index_gt]]$formula
gvif_data_gentime[highest_moderate_index_gt, 'GVIF_Adj']

# Plot 2

filtered_gentime_models <- gentime_combined_models[models_classifier(gentime_combined_models) %in% c("good", "moderate")]

# Run Q_best_models with the filtered models
Q_values <- c(30, 20, 10); porcent <- FALSE
Q_values <- c(50, 25, 10); porcent <- TRUE
Q_best_gentime_combined <- Q_best_models(Q_values, filtered_gentime_models, all_variables_combined, porcent)

var_freq_gentime_combined <- Q_best_gentime_combined$variable_freq
var_count_freq_df_gentime_combined <- Q_best_gentime_combined$variable_count_freq_df

var_freq_long_gentime_combined <- melt(var_freq_gentime_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_gentime_combined, var_freq_long_gentime_combined, all_variables_combined, "Generation Time")

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_gentime <- Q_best_models_table(Q_value, filtered_gentime_models, all_variables_combined, porcent)
Q_best_models_gentime <- lapply(Q_best_formulas_gentime, function(formula) {
  glm(as.formula(formula), data = simstats)
}) # Running best models

aic_values_gentime <- round(sapply(Q_best_models_gentime, AIC), 3)
bic_values_gentime <- round(sapply(Q_best_models_gentime, BIC), 3)

tex_file_gentime <- 'data/03-model-comparison-combined/03-best_models_table_gentime.tex'

stargazer(Q_best_models_gentime, type = "latex", out = tex_file_gentime,
          title = "Best 10 models for Generation Time",
          label = "tab:best_models_gentime",
          dep.var.labels = "Generation Time",
          #covariate.labels = best_variables_gentime, # Uncomment and define if needed
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_gentime),
                           c("BIC", bic_values_gentime)),
          no.space = TRUE,
          digits = 3)


################################## Rep Num Models ----------------------------

rep_num_combined_file <- 'data/03-model-comparison-combined/03-vif_combined_models_rep_num_0.RData'

if (!file.exists(rep_num_combined_file)) {
  print('Running rep_num_combined_models')
  run_combined_regression_models('rt_0', 'glm.nb', rep_num_combined_file, all_vars)
} else {
  print('Loading rep_num_combined_models')
  rep_num_combined_models <- readRDS(rep_num_combined_file)
}

# Plot 1

# Classify models and obtain GVIF data
gvif_data_rep_num <- models_classifier_with_gvif(rep_num_combined_models)

# Count the number of models in each category
good_count_rn <- sum(gvif_data_rep_num$Category == "good")
moderate_count_rn <- sum(gvif_data_rep_num$Category == "moderate")
bad_count_rn <- sum(gvif_data_rep_num$Category == "bad")
good_count_rn + moderate_count_rn

# Create the plot
ggplot(gvif_data_rep_num) +
  geom_point(data = subset(gvif_data_rep_num, Category == "bad"), aes(x = Index, y = GVIF_Adj), color = "red", size = 0.7) +
  geom_point(data = subset(gvif_data_rep_num, Category == "moderate"), aes(x = Index, y = GVIF_Adj), color = "orange", size = 0.8) +
  geom_point(data = subset(gvif_data_rep_num, Category == "good"), aes(x = Index, y = GVIF_Adj), color = "green", size = 0.8) +
  scale_color_manual(values = c("good" = "green", "moderate" = "orange", "bad" = "red", "unknown" = "gray")) +
  labs(title = "GVIF^(1/(2*Df)) for Models",
       x = "Model Index",
       y = "GVIF^(1/(2*Df))",
       caption = paste("Good:", good_count_rn, "| Moderate:", moderate_count_rn, "| Bad:", bad_count_rn)) +
  ylim(0, 75) +
  theme_minimal()

# Find the highest moderate index and its details
highest_moderate_index_rn <- max(gvif_data_rep_num$Index[gvif_data_rep_num$Category == "moderate"], na.rm = TRUE)
highest_moderate_index_rn
rep_num_combined_models[[highest_moderate_index_rn]]$formula
gvif_data_rep_num[highest_moderate_index_rn, 'GVIF_Adj']

# Plot 2

filtered_rep_num_models <- rep_num_combined_models[models_classifier(rep_num_combined_models) %in% c("good", "moderate")]

# Run Q_best_models with the filtered models
Q_values <- c(50, 25, 10); porcent <- TRUE
Q_best_rep_num_combined <- Q_best_models(Q_values, filtered_rep_num_models, all_variables_combined, porcent)

var_freq_rep_num_combined <- Q_best_rep_num_combined$variable_freq
var_count_freq_df_rep_num_combined <- Q_best_rep_num_combined$variable_count_freq_df

var_freq_long_rep_num_combined <- melt(var_freq_rep_num_combined, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")
graphic_generator(var_count_freq_df_rep_num_combined, var_freq_long_rep_num_combined, all_variables_combined, "Reproduction Number")

# Table

Q_value <- c(10); porcent <- FALSE
Q_best_formulas_rep_num <- Q_best_models_table(Q_value, filtered_rep_num_models, all_variables_combined, porcent)
Q_best_models_rep_num <- lapply(Q_best_formulas_rep_num, function(formula) {
  glm.nb(as.formula(formula), data = simstats)
}) # Running best models using glm.nb for negative binomial regression

aic_values_rep_num <- round(sapply(Q_best_models_rep_num, AIC), 3)
bic_values_rep_num <- round(sapply(Q_best_models_rep_num, BIC), 3)

tex_file_rep_num <- 'data/03-model-comparison-combined/03-best_models_table_rep_num_0.tex'

stargazer(Q_best_models_rep_num, type = "latex", out = tex_file_rep_num,
          title = "Best 10 models for Reproduction Number",
          label = "tab:best_models_rep_num",
          dep.var.labels = "Reproduction Number",
          #covariate.labels = best_variables_rep_num, # Uncomment and define if needed
          omit.stat = c("LL", "ser", "f", "aic"),
          omit.table.layout = "d", # not showing dependent var
          add.lines = list(c("AIC", aic_values_rep_num),
                           c("BIC", bic_values_rep_num)),
          no.space = TRUE,
          digits = 3)
