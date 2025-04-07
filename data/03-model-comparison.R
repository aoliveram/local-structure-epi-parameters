library(data.table)
library(ggplot2)
library(lmtest)
library(texreg)
library(dplyr)


################################## Data -------------------------------------

simstats <- fread("data/02-dataprep-network-stats-new-exposed.csv.gz")

simstats <- simstats[peak_preval > 1]

# Variable labels
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

# More legible names for nettype
simstats[, nettype := factor(
  nettype,
  levels = c("ergm", "sf", "swp01", "swp02", "degseq", "er"),
  labels = c("ERGM", "Scale-free", "Small-world (p=0.1)", "Small-world (p=0.2)", "Degree-sequence", "Erdos-Renyi")
)]

# Variables to iterate for
all_variables <- c("ergm_twopath", "ergm_balance", "ergm_triangle", 
                    "igraph_avg_degree", "igraph_avg_path_length", 
                    "igraph_density", "igraph_transitivity", "igraph_modularity")

# Extracting AIC and BIC
aic_bic_info <- function(model) {
  c(AIC = AIC(model), BIC = BIC(model))
}

################################## Peak Time ---------------------------------

peak_time_models <- list()

# Function for generate formulas in regressions
formula_generator_peak_time <- function(variables) {
  formula <- as.formula(paste("peak_time ~ I(factor(nettype))", 
                              paste(variables, collapse = " + "), sep = " + "))
  return(formula)}

# Add fixed effect only
peak_time_models[[1]] <- glm(peak_time ~ I(factor(nettype)), data = simstats)

# Run all others regressions
for (k in 1:length(all_variables)) {
  combinations <- combn(all_variables, k, simplify = FALSE)
  for (combination in combinations) {
    formula <- formula_generator_peak_time(combination)
    model <- glm(formula, data = simstats)
    peak_time_models <- append(peak_time_models, list(model))
  }}

# Obtaining AIC and BIC for all the models
inf_crit_peak_time <- as.data.frame(do.call(rbind, lapply(peak_time_models, aic_bic_info)))
inf_crit_peak_time$model <- paste0("model_", seq_along(peak_time_models))

# Normalize AIC and BIC values
mean_AIC <- mean(inf_crit_peak_time$AIC); sd_AIC <- sd(inf_crit_peak_time$AIC)
mean_BIC <- mean(inf_crit_peak_time$BIC); sd_BIC <- sd(inf_crit_peak_time$BIC)
inf_crit_peak_time$AIC_z <- (inf_crit_peak_time$AIC - mean_AIC) / sd_AIC
inf_crit_peak_time$BIC_z <- (inf_crit_peak_time$BIC - mean_BIC) / sd_BIC

# Mean of both normalized AIC and BIC
inf_crit_peak_time$Mean_z <- rowMeans(inf_crit_peak_time[, c("AIC_z", "BIC_z")])
Mean_z_peak_time_sort <- inf_crit_peak_time[order(inf_crit_peak_time$Mean_z), ]

# Selecting Q best models
Q_values <- c(10, 20, 30)
variable_freq_peak_time <- data.frame(Variable = all_variables)

for (Q in Q_values) {
  top_Q_models_peak_time <- head(Mean_z_peak_time_sort, Q)
  top_Q_formulas_peak_time <- sapply(top_Q_models_peak_time$model, function(x) peak_time_models[[as.numeric(sub("model_", "", x))]]$formula)
  
  # Frequency of variables in Q best formulas
  freq <- sapply(all_variables, function(var) {
    sum(sapply(top_Q_formulas_peak_time, function(formula) grepl(var, formula)))
  })
  
  variable_freq_peak_time[[paste0("Q_", Q)]] <- freq / Q
}

variable_freq_peak_time <- variable_freq_peak_time[order(-variable_freq_peak_time$Q_30), ]
variable_freq_peak_time

# Generate top_Q_formulas_peak_time for each Q value
top_Q_formulas_peak_time <- list()
for (Q in Q_values) {
  top_Q_models_peak_time <- head(Mean_z_peak_time_sort, Q)
  top_Q_formulas_peak_time[[paste0("top_Q_formulas_peak_time_Q", Q)]] <- sapply(top_Q_models_peak_time$model, function(x) peak_time_models[[as.numeric(sub("model_", "", x))]]$formula)
}
top_Q_formulas_peak_time

variable_count_freq_list_peak_time <- list()
for (i in seq_along(Q_values)) {
  Q <- Q_values[i]
  formulas <- top_Q_formulas_peak_time[[i]]
  
  variable_counts <- sapply(formulas, function(formula) {
    variables <- all.vars(formula)[-1]  # Excluir la variable dependiente
    variables <- variables[variables != "factor(nettype)"]
    return(length(variables))
  }) |> table()
  
  variable_count_freq_list_peak_time[[paste0("Q_", Q)]] <- variable_counts
}
variable_count_freq_list_peak_time

################################## Peak Preval ---------------------------------

peak_preval_models <- list()

# Function for generate formulas in regressions
formula_generator_peak_preval <- function(variables) {
  formula <- as.formula(paste("peak_preval ~ I(factor(nettype))", 
                              paste(variables, collapse = " + "), sep = " + "))
  return(formula)}

# Add fixed effect only
peak_preval_models[[1]] <- glm(peak_preval ~ I(factor(nettype)), data = simstats)

# Run all others regressions
for (k in 1:length(all_variables)) {
  combinations <- combn(all_variables, k, simplify = FALSE)
  for (combination in combinations) {
    formula <- formula_generator_peak_preval(combination)
    model <- glm(formula, data = simstats)
    peak_preval_models <- append(peak_preval_models, list(model))
  }}

# Obtaining AIC and BIC for all the models
inf_crit_peak_preval <- as.data.frame(do.call(rbind, lapply(peak_preval_models, aic_bic_info)))
inf_crit_peak_preval$model <- paste0("model_", seq_along(peak_preval_models))

# Normalize AIC and BIC values
mean_AIC <- mean(inf_crit_peak_preval$AIC); sd_AIC <- sd(inf_crit_peak_preval$AIC)
mean_BIC <- mean(inf_crit_peak_preval$BIC); sd_BIC <- sd(inf_crit_peak_preval$BIC)
inf_crit_peak_preval$AIC_z <- (inf_crit_peak_preval$AIC - mean_AIC) / sd_AIC
inf_crit_peak_preval$BIC_z <- (inf_crit_peak_preval$BIC - mean_BIC) / sd_BIC

# Mean of both normalized AIC and BIC
inf_crit_peak_preval$Mean_z <- rowMeans(inf_crit_peak_preval[, c("AIC_z", "BIC_z")])
Mean_z_peak_preval_sort <- inf_crit_peak_preval[order(inf_crit_peak_preval$Mean_z), ]

# Selecting Q best models
Q_values <- c(10, 20, 30)
variable_freq_peak_preval <- data.frame(Variable = all_variables)

for (Q in Q_values) {
  top_Q_models_peak_preval <- head(Mean_z_peak_preval_sort, Q)
  top_Q_formulas_peak_preval <- sapply(top_Q_models_peak_preval$model, function(x) peak_preval_models[[as.numeric(sub("model_", "", x))]]$formula)
  
  # Frequency of variables in Q best formulas
  freq <- sapply(all_variables, function(var) {
    sum(sapply(top_Q_formulas_peak_preval, function(formula) grepl(var, formula)))
  })
  
  variable_freq_peak_preval[[paste0("Q_", Q)]] <- freq / Q
}

variable_freq_peak_preval <- variable_freq_peak_preval[order(-variable_freq_peak_preval$Q_30), ]
variable_freq_peak_preval

# Generate top_Q_formulas_peak_preval for each Q value
top_Q_formulas_peak_preval <- list()
for (Q in Q_values) {
  top_Q_models_peak_preval <- head(Mean_z_peak_preval_sort, Q)
  top_Q_formulas_peak_preval[[paste0("top_Q_formulas_peak_preval_Q", Q)]] <- sapply(top_Q_models_peak_preval$model, function(x) peak_preval_models[[as.numeric(sub("model_", "", x))]]$formula)
}
top_Q_formulas_peak_preval


variable_count_freq_list_peak_preval <- list()
for (i in seq_along(Q_values)) {
  Q <- Q_values[i]
  formulas <- top_Q_formulas_peak_preval[[i]]
  
  variable_counts <- sapply(formulas, function(formula) {
    variables <- all.vars(formula)[-1]  # Excluir la variable dependiente
    variables <- variables[variables != "factor(nettype)"]
    return(length(variables))
  }) |> table()
  
  variable_count_freq_list_peak_preval[[paste0("Q_", Q)]] <- variable_counts
}
variable_count_freq_list_peak_preval


################################## Gentime ---------------------------------

gentime_models <- list()

# Function for generate formulas in regressions
formula_generator_gentime <- function(variables) {
  formula <- as.formula(paste("gentime ~ I(factor(nettype))", 
                              paste(variables, collapse = " + "), sep = " + "))
  return(formula)}

# Add fixed effect only
gentime_models[[1]] <- glm(gentime ~ I(factor(nettype)), data = simstats)

# Run all others regressions
for (k in 1:length(all_variables)) {
  combinations <- combn(all_variables, k, simplify = FALSE)
  for (combination in combinations) {
    formula <- formula_generator_gentime(combination)
    model <- glm(formula, data = simstats)
    gentime_models <- append(gentime_models, list(model))
  }}

# Obtaining AIC and BIC for all the models
inf_crit_gentime <- as.data.frame(do.call(rbind, lapply(gentime_models, aic_bic_info)))
inf_crit_gentime$model <- paste0("model_", seq_along(gentime_models))

# Normalize AIC and BIC values
mean_AIC <- mean(inf_crit_gentime$AIC); sd_AIC <- sd(inf_crit_gentime$AIC)
mean_BIC <- mean(inf_crit_gentime$BIC); sd_BIC <- sd(inf_crit_gentime$BIC)
inf_crit_gentime$AIC_z <- (inf_crit_gentime$AIC - mean_AIC) / sd_AIC
inf_crit_gentime$BIC_z <- (inf_crit_gentime$BIC - mean_BIC) / sd_BIC

# Mean of both normalized AIC and BIC
inf_crit_gentime$Mean_z <- rowMeans(inf_crit_gentime[, c("AIC_z", "BIC_z")])
Mean_z_gentime_sort <- inf_crit_gentime[order(inf_crit_gentime$Mean_z), ]

# Selecting Q best models
Q_values <- c(10, 20, 30)
variable_freq_gentime <- data.frame(Variable = all_variables)

for (Q in Q_values) {
  top_Q_models_gentime <- head(Mean_z_gentime_sort, Q)
  top_Q_formulas_gentime <- sapply(top_Q_models_gentime$model, function(x) peak_preval_models[[as.numeric(sub("model_", "", x))]]$formula)
  
  # Frequency of variables in Q best formulas
  freq <- sapply(all_variables, function(var) {
    sum(sapply(top_Q_formulas_gentime, function(formula) grepl(var, formula)))
  })
  
  variable_freq_gentime[[paste0("Q_", Q)]] <- freq / Q
}

variable_freq_gentime <- variable_freq_gentime[order(-variable_freq_gentime$Q_30), ]
variable_freq_gentime

# Generate top_Q_formulas_gentime for each Q value
top_Q_formulas_gentime <- list()
for (Q in Q_values) {
  top_Q_models_gentime <- head(Mean_z_gentime_sort, Q)
  top_Q_formulas_gentime[[paste0("top_Q_formulas_gentime_Q", Q)]] <- sapply(top_Q_models_gentime$model, function(x) peak_preval_models[[as.numeric(sub("model_", "", x))]]$formula)
}
top_Q_formulas_gentime

variable_count_freq_list_gentime <- list()
for (i in seq_along(Q_values)) {
  Q <- Q_values[i]
  formulas <- top_Q_formulas_gentime[[i]]
  
  variable_counts <- sapply(formulas, function(formula) {
    variables <- all.vars(formula)[-1]  # Excluir la variable dependiente
    variables <- variables[variables != "factor(nettype)"]
    return(length(variables))
  }) |> table()
  
  variable_count_freq_list_gentime[[paste0("Q_", Q)]] <- variable_counts
}
variable_count_freq_list_gentime


################################## Rep Num -----------------------------------

library(MASS)

rep_num_models <- list()

# Function for generate formulas in regressions
formula_generator_rep_num <- function(variables) {
  formula <- as.formula(paste("rt ~ I(factor(nettype))", 
                              paste(variables, collapse = " + "), sep = " + "))
  return(formula)}

# Add fixed effect only
rep_num_models[[1]] <- glm.nb(rt ~ I(factor(nettype)), data = simstats)

start_time <- Sys.time()
# Run all others regressions
for (k in 1:length(all_variables)) {
  combinations <- combn(all_variables, k, simplify = FALSE)
  for (combination in combinations) {
    formula <- formula_generator_rep_num(combination)
    model <- glm.nb(formula, data = simstats)
    rep_num_models <- append(rep_num_models, list(model))
  }}
end_time <- Sys.time()
time_taken_no_parall <- end_time - start_time # ~170 seg

# Obtaining AIC and BIC for all the models
inf_crit_rep_num <- as.data.frame(do.call(rbind, lapply(rep_num_models, aic_bic_info)))
inf_crit_rep_num$model <- paste0("model_", seq_along(rep_num_models))

# Normalize AIC and BIC values
mean_AIC <- mean(inf_crit_rep_num$AIC); sd_AIC <- sd(inf_crit_rep_num$AIC)
mean_BIC <- mean(inf_crit_rep_num$BIC); sd_BIC <- sd(inf_crit_rep_num$BIC)
inf_crit_rep_num$AIC_z <- (inf_crit_rep_num$AIC - mean_AIC) / sd_AIC
inf_crit_rep_num$BIC_z <- (inf_crit_rep_num$BIC - mean_BIC) / sd_BIC

# Mean of both normalized AIC and BIC
inf_crit_rep_num$Mean_z <- rowMeans(inf_crit_rep_num[, c("AIC_z", "BIC_z")])
Mean_z_rep_num_sort <- inf_crit_rep_num[order(inf_crit_rep_num$Mean_z), ]

# Selecting Q best models
Q_values <- c(10, 20, 30)
variable_freq_rep_num <- data.frame(Variable = all_variables)

for (Q in Q_values) {
  #Q <- 10
  top_Q_models_rep_num <- head(Mean_z_rep_num_sort, Q)
  top_Q_formulas_rep_num <- sapply(top_Q_models_rep_num$model, function(x) {
    model_index <- as.numeric(sub("model_", "", x))
    formula(rep_num_models[[model_index]])
  })
  
  # Frequency of variables in Q best formulas
  #var <- all_variables[1]
  freq <- sapply(all_variables, function(var) {
    sum(sapply(top_Q_formulas_rep_num, function(formula) grepl(var, formula)))
  })
  
  variable_freq_rep_num[[paste0("Q_", Q)]] <- freq / Q
}

variable_freq_rep_num <- variable_freq_rep_num[order(-variable_freq_rep_num$Q_30), ]
variable_freq_rep_num

# Generate top_Q_formulas_rep_num for each Q value
top_Q_formulas_rep_num <- list()
for (Q in Q_values) {
  top_Q_models_rep_num <- head(Mean_z_rep_num_sort, Q)
  top_Q_formulas_rep_num[[paste0("top_Q_formulas_rep_num_Q", Q)]] <- sapply(top_Q_models_rep_num$model, function(x) {
    model_index <- as.numeric(sub("model_", "", x))
    formula(rep_num_models[[model_index]])
  })
}
top_Q_formulas_rep_num


variable_count_freq_list_rep_num <- list()
for (i in seq_along(Q_values)) {
  Q <- Q_values[i]
  formulas <- top_Q_formulas_rep_num[[i]]
  
  variable_counts <- sapply(formulas, function(formula) {
    variables <- all.vars(formula)[-1]  # Excluir la variable dependiente
    variables <- variables[variables != "factor(nettype)"]
    return(length(variables))
  }) |> table()
  
  variable_count_freq_list_rep_num[[paste0("Q_", Q)]] <- variable_counts
}
variable_count_freq_list_rep_num


################################## Plot variables ------------------------------

library(ggplot2)
library(reshape2)

# Long format for ggplot2 - Peak Time 
var_freq_long_peak_time <- melt(variable_freq_peak_time, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")

ggplot(var_freq_long_peak_time, aes(x = Variable, y = Frequency, color = Q, group = Q)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Frecuencia de Variables para Peak Time",
       x = NULL,
       y = "Frecuencia",
       color = "Q")

# Long format for ggplot2 - Peak Preval
var_freq_long_peak_preval <- melt(variable_freq_peak_preval, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")

ggplot(var_freq_long_peak_preval, aes(x = Variable, y = Frequency, color = Q, group = Q)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Frecuencia de Variables para Peak Prevalence",
       x = NULL,
       y = "Frecuencia",
       color = "Q")

# Long format for ggplot2 - Gentime 
var_freq_long_gentime <- melt(variable_freq_gentime, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")

ggplot(var_freq_long_gentime, aes(x = Variable, y = Frequency, color = Q, group = Q)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Frecuencia de Variables para Gentime",
       x = NULL,
       y = "Frecuencia",
       color = "Q")

# Long format for ggplot2 - Rep Num 
var_freq_long_rep_num <- melt(variable_freq_rep_num, id.vars = "Variable", variable.name = "Q", value.name = "Frequency")

ggplot(var_freq_long_rep_num, aes(x = Variable, y = Frequency, color = Q, group = Q)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Frecuencia de Variables para Rep Num",
       x = NULL,
       y = "Frecuencia",
       color = "Q")

variable_count_freq_list_peak_time
variable_count_freq_list_peak_preval
variable_count_freq_list_gentime
variable_count_freq_list_rep_num

################################## Plots Regressions ---------------------------

# Better labels

coefmap <- list(
  "I(factor(nettype))Scale-free" = "Scale-free",
  "I(factor(nettype))Small-world (p=0.1)" = "Small-world (p=0.1)",
  "I(factor(nettype))Small-world (p=0.2)" = "Small-world (p=0.2)",
  "I(factor(nettype))Degree-sequence" = "Degree-sequence",
  "I(factor(nettype))Erdos-Renyi" = "Erdos-Renyi",
  igraph_avg_degree = "Average degree",
  ergm_twopath = "Two-path",
  igraph_transitivity = "Transitivity",
  ergm_triangle = "Triangles",
  ergm_balance = "Balance"
)

results <- list(
  `Peak preval` = model_peak_preval5,
  `Peak time`   = model_peak_time6,
  `Gen time`    = model_gentime5,
  Rt            = model_rt5
)


hreg <- htmlreg(
  results,
  custom.coef.map = coefmap,
  single.row = TRUE,
  custom.note = "%stars. In the case of Rt, we used a negative binomial regression model.",
  outer.rules = 0,
  inline.css = TRUE,
  groups = list(`Fixed effects` = 1:5, `Network structure` = 6:10),
  caption = NULL
)
print(hreg, file = "03-regression-comparison.html")

library(jtools)

coefmap2 <- setNames(names(coefmap), unlist(coefmap))

do.call(
  plot_coefs,
  list(
    results,
    coefs = coefmap2,
    groups = list(
      `Fixed effects` = names(coefmap2)[1:5],
      `Network structure` = names(coefmap2)[6:10]
    )
  )) +
  scale_colour_brewer(palette = "Set1")

ggsave(
  "figures/03-regression-comparison.png", width = 4, height = 2.5,
  scale = 1.5,
  dpi = 300
)

screenreg(
  list(
    `Peak preval` = model_peak_preval5,
    `Peak time`   = model_peak_time6,
    `Gen time`    = model_gentime5,
    Rt             = model_rt5
  ),
  custom.coef.map = coefmap,
  single.row = TRUE,
  custom.note = "%stars. In the case of Rt, we used a negative binomial regression model."
)
