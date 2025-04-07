###################################### -----------------------------------

library(dplyr)
library(purrr)
library(data.table)

rep_num_sims_file <- 'data/02-reproductive-number-mean/02-reproductive-number-sims.RData'

if (!file.exists(rep_num_sims_file)) {
  print('Collecting data')
  
  simfiles <- list.files("data/sims", pattern = "-sim-[0-9]+\\.rds$", full.names = TRUE)
  
  rep_num_avg_all <- list()
  
  for (i in seq_along(simfiles)) {
    x <- readRDS(simfiles[i])
    rep_num_avg <- x$repnum$avg
    
    # Verificar si rep_num_avg es NULL
    if (!is.null(rep_num_avg)) {
      rep_num_avg_all[[i]] <- list(rep_num_avg = rep_num_avg, simfile = simfiles[i])
    }
  }
  
  # Convertir la lista de listas en un data.frame, excluyendo elementos NULL
  rep_num_avg_df <- map_dfr(rep_num_avg_all, function(sublist) {
    if (!is.null(sublist$rep_num_avg)) {
      netid <- gsub(".+/([0-9]+-[a-z0-9]+)-sim.+", "\\1", sublist$simfile)
      nettype <- gsub(".+-([a-z0-9]+)-sim.+", "\\1", sublist$simfile)
      data.frame(matrix(unlist(sublist$rep_num_avg), nrow = 1, byrow = TRUE),
                 netid = netid,
                 nettype = nettype)
    } else {
      NULL
    }
  })
  
  # Guardar el data.frame en un archivo RDS
  saveRDS(rep_num_avg_df, file = rep_num_sims_file)
} else {
  print('Loading data')
  rep_num_avg_df <- readRDS(rep_num_sims_file)
}

setDT(rep_num_avg_df) # convertir data.frame

rep_num_avg_df[, nettype := factor(
  nettype,
  levels = c("ergm", "sf", "swp01", "swp02", "degseq", "er"),
  labels = c("ERGM", "Scale-free", "Small-world (p=0.1)", "Small-world (p=0.2)", "Degree-sequence", "Erdös–Rényi")
)]

## PLOT

library(tidyr)
library(ggplot2)

# Convertir data.frame a formato largo
long_df <- melt(rep_num_avg_df, id.vars = c("netid", "nettype"), variable.name = "day", value.name = "value")

long_df[value == 0.0, value := NA] # Tratar los valores 0.0 como NA
long_df <- long_df[!is.na(value)] # Filtrar los NA

nettypes <- unique(long_df$nettype)
plots <- list()

for (net in nettypes) {
  
  # Filtro por tipo de red
  net_df <- long_df[nettype == net]
  
  stats_df <- net_df %>%
    group_by(day) %>%
    summarize(
      avg_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE)
    )
  
  stats_df <- stats_df %>%
    mutate(sd_value = ifelse(is.na(sd_value), 0, sd_value))
  
  net_df_out <- net_df %>%
    inner_join(stats_df, by = "day") %>%
    filter(value < avg_value - sd_value | value > avg_value + sd_value)
  
  p <- ggplot() +
    geom_point(data = net_df, aes(x = as.numeric(day), y = value), alpha = 0.5, color = "darkgreen", size = 0.3) + 
    geom_point(data = net_df_out, aes(x = as.numeric(day), y = value), alpha = 0.5, color = "green", size = 0.3) + 
    geom_line(data = stats_df, aes(x = as.numeric(day), y = avg_value), color = "red", size = 1) +
    scale_y_log10(limits = c(0.05, 7.5)) +
    scale_x_continuous(breaks = seq(0, max(as.numeric(long_df$day)), by = 10), labels = seq(0, max(as.numeric(long_df$day)), by = 10)) + 
    labs(title = paste("Network type:", net), x = "Day", y = "Reproductive number (log10)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = 52.5, y = 1.5, label = "R=1", color = "black", size = 5, hjust = 0)
  
  ggsave(filename = paste0("data/02-reproductive-number-mean/02-reproductive-number-mean-", net, ".pdf"), plot = p, width = 2.5, height = 2)
  
  plots[[net]] <- p
}

for (net in nettypes) {
  print(plots[[net]])
}
