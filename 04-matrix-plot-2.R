library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)
library(dplyr)
library(reshape2)
library(stringr)

# Load data and models
base_dir <- "03-model-comparison-combined-new-var"
output_dir <- "03-model-comparison-combined-new-var"
dependent_vars <- c("peak_preval", "peak_time", "gentime", "rt")
dependent_vars_labels <- c("Peak Prevalence", "Peak Time", "Generation Time", "Reproductive Number")

# Helper function to extract model data
extract_model_data <- function(models_list) {
  aic_values <- unlist(lapply(models_list, function(x) x$aic))
  bic_values <- unlist(lapply(models_list, function(x) x$bic))
  formulas <- as.character(unlist(lapply(models_list, function(x) x$formula)))
  
  df <- data.frame(
    id = 1:length(models_list),
    aic = aic_values,
    bic = bic_values,
    formula = formulas,
    stringsAsFactors = FALSE
  )
  
  # Calculate Z-scores
  mean_AIC <- mean(aic_values)
  mean_BIC <- mean(bic_values)
  sd_AIC <- sd(aic_values)
  sd_BIC <- sd(bic_values)
  
  df$AIC_z <- (df$aic - mean_AIC) / sd_AIC
  df$BIC_z <- (df$bic - mean_BIC) / sd_BIC
  df$Mean_z <- rowMeans(df[, c("AIC_z", "BIC_z")])
  
  # Rank models
  df <- df[order(df$Mean_z), ]
  df$rank <- 1:nrow(df)
  
  return(df)
}

# Helper function to check significance
check_significance <- function(model_obj, var_name) {
  if (!is.null(model_obj$p_vals) && var_name %in% names(model_obj$p_vals)) {
    return(model_obj$p_vals[var_name] < 0.05)
  } else if (is.null(model_obj$p_vals)) {
    return(var_name %in% labels(terms(as.formula(model_obj$formula))))
  }
  return(FALSE)
}

# Helper to parse variable name and functional form
parse_variable <- function(term) {
  # Shapes: Linear (Circle), Quadratic (Square), Logarithmic (Triangle)
  
  base_var <- term
  form <- "Linear"
  
  if (grepl("^I\\(.*\\^2\\)$", term)) {
    form <- "Quadratic"
    base_var <- gsub("^I\\((.*)\\^2\\)$", "\\1", term)
  } else if (grepl("^log\\(.*\\)$", term)) {
    form <- "Logarithmic"
    base_var <- gsub("^log\\((.*)\\)$", "\\1", term)
  }
  
  # Clean base variable name
  base_var <- gsub("igraph_", "", base_var)
  base_var <- gsub("degree_variance", "Variance", base_var)
  base_var <- gsub("avg_", "Avg. ", base_var)
  base_var <- gsub("local_", "Local ", base_var)
  base_var <- gsub("assortativity", "Assortativity", base_var)
  base_var <- gsub("path_length", "Path Length", base_var)
  base_var <- gsub("transitivity", "Transitivity", base_var)
  base_var <- gsub("modularity", "Modularity", base_var)
  base_var <- gsub("degree", "Degree", base_var)
  
  return(list(base = base_var, form = form))
}

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) > 0) legend <- tmp$grobs[[leg]]
  else legend <- NULL
  return(legend)
}

# Main plotting function
create_matrix_plot_v2 <- function(dep_var_key, dep_var_label) {
  # Try specific filename for reproductive number if needed
  rdata_file <- file.path(base_dir, paste0("03-combined_models_", dep_var_key, ".RData"))
  if (!file.exists(rdata_file) && dep_var_key == "rt") {
     rdata_file <- file.path(base_dir, "03-combined_models_rep_num.RData")
  }

  if (!file.exists(rdata_file)) {
    message(paste("File not found:", rdata_file))
    return(NULL)
  }
  
  models_list <- readRDS(rdata_file)
  models_list <- Filter(Negate(is.null), models_list)
  
  model_stats <- extract_model_data(models_list)
  
  # Select top 50% models
  n_total <- nrow(model_stats)
  n_top50 <- floor(0.50 * n_total)
  top_models_df <- model_stats[1:n_top50, ]
  
  # Define Q thresholds
  n_top10 <- floor(0.10 * n_total)
  n_top25 <- floor(0.25 * n_total)
  
  top_models_df$Q_group <- "Q_50"
  top_models_df$Q_group[top_models_df$rank <= n_top25] <- "Q_25"
  top_models_df$Q_group[top_models_df$rank <= n_top10] <- "Q_10"
  
  # Identify all unique terms
  all_terms <- unique(unlist(lapply(top_models_df$formula, function(f) {
    attr(terms(as.formula(f)), "term.labels")
  })))
  all_terms <- all_terms[all_terms != "I(factor(nettype))"] 
  
  # Build matrix data
  matrix_data <- list()
  
  for (i in 1:nrow(top_models_df)) {
    model_idx <- top_models_df$id[i]
    model_obj <- models_list[[model_idx]]
    rank <- top_models_df$rank[i]
    q_group <- top_models_df$Q_group[i]
    
    for (term in all_terms) {
      if (check_significance(model_obj, term)) {
        parsed <- parse_variable(term)
        matrix_data[[length(matrix_data) + 1]] <- data.frame(
          rank = rank,
          term = term,
          base_var = parsed$base,
          form = parsed$form,
          Q_group = q_group,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  matrix_df <- do.call(rbind, matrix_data)
  
  # Frequency Calculations (Sum of frequencies)
  # For each Q group, we count occurences of each base_var (across all forms)
  # Then normalize by N_models_in_Q
  
  # Filter data for each threshold
  data_Q10 <- matrix_df[matrix_df$rank <= n_top10, ]
  data_Q25 <- matrix_df[matrix_df$rank <= n_top25, ]
  data_Q50 <- matrix_df[matrix_df$rank <= n_top50, ]
  
  # Count occurences per base_var
  count_Q10 <- data_Q10 %>% group_by(base_var) %>% summarise(count = n())
  count_Q25 <- data_Q25 %>% group_by(base_var) %>% summarise(count = n())
  count_Q50 <- data_Q50 %>% group_by(base_var) %>% summarise(count = n())
  
  # Normalize to get frequency (avg appearances per model)
  # NOTE: A variable can appear multiple times in one model (e.g. Linear + Quad).
  # This metric reflects "intensity" of the variable.
  count_Q10$freq <- count_Q10$count / n_top10
  count_Q25$freq <- count_Q25$count / n_top25
  count_Q50$freq <- count_Q50$count / n_top50
  
  all_expected_vars <- c("Avg. Degree", "Avg. Path Length", "Local Transitivity", 
                       "Modularity", "Assortativity", "Variance")
  
  freq_df <- data.frame(base_var = all_expected_vars)
  
  freq_df <- merge(freq_df, count_Q10[, c("base_var", "freq")], by="base_var", all.x=TRUE); colnames(freq_df)[2] <- "freq_Q10"
  freq_df <- merge(freq_df, count_Q25[, c("base_var", "freq")], by="base_var", all.x=TRUE); colnames(freq_df)[3] <- "freq_Q25"
  freq_df <- merge(freq_df, count_Q50[, c("base_var", "freq")], by="base_var", all.x=TRUE); colnames(freq_df)[4] <- "freq_Q50"
  freq_df[is.na(freq_df)] <- 0
  
  # Calculate Total Sum of Frequencies for ordering
  freq_df$total_sum <- freq_df$freq_Q10 + freq_df$freq_Q25 + freq_df$freq_Q50
  
  # Organize order (Highest total frequency at top)
  ordered_vars <- freq_df$base_var[order(freq_df$total_sum, decreasing = FALSE)] # Ascending for Coord Flip logic?
  # Actually, for geom_point with categorical Y, the first level is at bottom.
  # So we want "Best" variable as the LAST level to appear at TOP.
  matrix_df$base_var <- factor(matrix_df$base_var, levels = ordered_vars)
  
  # Prepare Bar Data (Long format)
  freq_long <- melt(freq_df[, c("base_var", "freq_Q10", "freq_Q25", "freq_Q50")], 
                    id.vars = "base_var", variable.name = "Q_group", value.name = "frequency")
  
  # Stack order: Red (Q10), Green (Q25), Blue (Q50)
  # In ggplot stack, factor order determines order. usually reverse.
  freq_long$Q_group <- factor(freq_long$Q_group, levels = c("freq_Q50", "freq_Q25", "freq_Q10"))
  freq_long$base_var <- factor(freq_long$base_var, levels = ordered_vars)
  
  
  # --- Plot 1: The Matrix (Middle/Left) ---
  p_matrix <- ggplot(matrix_df, aes(x = rank, y = base_var)) +
    # Grid lines first
    geom_hline(yintercept = 1:length(ordered_vars), color = "grey90", linetype = "dashed") +
    # Points with Shapes
    geom_point(aes(shape = form, color = Q_group), size = 6, stroke = 1.5) +
    # Scales
    scale_shape_manual(values = c("Linear" = 16, "Quadratic" = 15, "Logarithmic" = 17),  # Circle, Square, Triangle
                       labels = c("Linear", "Quadratic", "Logarithmic")) +
    scale_color_manual(values = c("Q_50" = "deepskyblue3", "Q_25" = "chartreuse3", "Q_10" = "red"),
                       labels = c("Top 50%", "Top 25%", "Top 10%")) +
    scale_x_continuous(breaks = seq(1, n_top50, by = 5), expand = c(0.02, 0.02)) +
    # Ensure all variables are shown on Y axis
    scale_y_discrete(drop = FALSE) + 
    labs(x = "Model Rank (Best to Worst)", y = NULL, color = "Top Q%", shape = "Functional Form") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 24, face = "bold", color = "black"),
      axis.text.x = element_text(size = 20),
      axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 10)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 20)
    ) +
    guides(color = "none", shape = guide_legend(title="Functional Form", nrow=1, override.aes = list(size=10, color="black")))
    
  # --- Plot 2: The Bars (Right) ---
  # We want horizontal bars aligned with the matrix rows
  p_bars <- ggplot(freq_long, aes(x = base_var, y = frequency, fill = Q_group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("freq_Q10" = "red", "freq_Q25" = "chartreuse3", "freq_Q50" = "deepskyblue3"),
                      labels = c("10%", "25%", "50%")) +
    coord_flip() + # Flip to make bars horizontal
    scale_x_discrete(drop = FALSE) + # Ensure all vars shown
    labs(y = "Sum of Frequencies", x = NULL, fill = "Top Q%") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(), # Hide Y labels as they are in the matrix plot
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 10)),
      axis.text.x = element_text(size = 20),
      legend.position = "bottom",
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 20)
    ) +
    guides(fill = guide_legend(title = "Top Q%", nrow = 1))

  # Extract legends
  legend_shape <- get_legend(p_matrix)
  legend_fill <- get_legend(p_bars)
  
  # Remove legends from main plots
  p_matrix_nolegend <- p_matrix + theme(legend.position = "none")
  p_bars_nolegend <- p_bars + theme(legend.position = "none")
  
  # Final Assemble
  pdf_file <- file.path(output_dir, paste0("Matrix_Plot_v2_", dep_var_label, ".pdf"))
  pdf(pdf_file, width = 20, height = 8)
  
  # Layout:
  # Title
  # Matrix | Bars
  # Shape Leg | Fill Leg
  
  grid.arrange(
    arrangeGrob(p_matrix_nolegend, p_bars_nolegend, nrow = 1, widths = c(3, 1)),
    arrangeGrob(legend_shape, legend_fill, nrow = 1, widths = c(3, 1)),
    nrow = 2,
    heights = c(10, 1),
    top = textGrob(dep_var_label, gp = gpar(fontsize = 36, fontface = "bold"))
  )
  
  dev.off()
  message(paste("Saved plot for", dep_var_label))
}

# Run for all dependent variables
for (i in seq_along(dependent_vars)) {
  create_matrix_plot_v2(dependent_vars[i], dependent_vars_labels[i])
}
