library(network)
library(ergm)
library(igraph)
library(intergraph)
library(data.table)
library(parallel)

# Setup
ncores <- 8
set.seed(1231)

# Files
networks <- list.files("00-data/graphs", pattern = "[0-9]+-(ergm|sf|swp[0-9]{2}|degseq|er)\\.rds$", full.names = TRUE)
fn_sim <- "00-data/02-dataprep-sim.csv.gz"
fn_out <- "00-data/02-dataprep-results-2.csv.gz"

# 1. Compute Network Statistics
message("Computing network statistics...")

compute_net_stats <- function(fn) {
  tryCatch({
    # Read as network object (for ERGM)
    net_ergm <- readRDS(fn)
    
    # ERGM Stats
    # Balance and Triangles
    ergm_stats <- summary(net_ergm ~ balance + triangle)
    names(ergm_stats) <- c("balance", "triangles")
    
    # Convert to igraph
    net_ig <- asIgraph(net_ergm)
    
    # IGRAPH Stats
    # 1. Avg Degree
    avg_degree <- mean(degree(net_ig), na.rm = TRUE)
    
    # 2. Avg Path Length
    avg_path_length <- mean_distance(net_ig)
    
    # 3. Local Clustering (Transitivity type="average") & Global
    local_transitivity <- transitivity(net_ig, type = "average")
    global_transitivity <- transitivity(net_ig, type = "global")
    
    # 4. Modularity
    # Need community detection first
    comm <- cluster_fast_greedy(net_ig)
    modularity_val <- modularity(comm)
    
    # 7. Variance <k^2> and Variance(k)
    deg <- degree(net_ig)
    degree_moment_2 <- mean(deg^2) # <k^2>
    degree_var_stat <- var(deg)    # Variance
    
    # 8. Assortativity
    assortativity_val <- assortativity_degree(net_ig)
    
    # 9. Spectral Radius
    spectral_radius <- eigen_centrality(net_ig)$value
    
    # 10. Avg Betweenness
    avg_betweenness <- mean(betweenness(net_ig), na.rm = TRUE)
    
    data.table(
      netfile = fn,
      ergm_balance = ergm_stats["balance"],
      ergm_triangles = ergm_stats["triangles"],
      igraph_avg_degree = avg_degree,
      igraph_avg_path_length = avg_path_length,
      igraph_local_transitivity = local_transitivity,
      igraph_transitivity = global_transitivity, # Keeping original name for global
      igraph_modularity = modularity_val,
      igraph_degree_moment_2 = degree_moment_2,
      igraph_degree_variance = degree_var_stat,
      igraph_assortativity = assortativity_val,
      igraph_spectral_radius = spectral_radius,
      igraph_avg_betweenness = avg_betweenness
    )
    
  }, error = function(e) {
    message(paste("Error in", fn, ":", conditionMessage(e)))
    return(NULL)
  })
}

# Run parallel
# Using mclapply for forking on macOS
net_stats_list <- mclapply(networks, compute_net_stats, mc.cores = ncores)
net_stats <- rbindlist(net_stats_list)

# Add nettype and netid
net_stats[, nettype := gsub(".+-([a-z0-9]+)\\.rds", "\\1", netfile)]
net_stats[, netid   := gsub(".+/([0-9]+-[a-z0-9]+)\\.rds", "\\1", netfile)]

# 2. Simulation Results
message("Processing simulation results...")

if (file.exists(fn_sim)) {
  message("Loading existing simulation results...")
  simres <- fread(fn_sim)
} else {
  message("Computing simulation results from scratch (this might take a while)...")
  simfiles <- list.files("00-data/sims", pattern = "-sim-[0-9]+\\.rds$", full.names = TRUE)
  
  process_sim <- function(fn) {
    tryCatch({
      x <- readRDS(fn)
      
      # Extract metrics
      peak_idx    <- which.max(x$incidence$Exposed)
      peak_preval <- x$incidence$Exposed[peak_idx]
      peak_time   <- as.integer(rownames(x$incidence)[peak_idx])
      
      rt_idx     <- with(x$repnum, which.min(abs(peak_time - date)))
      rt         <- x$repnum$avg[rt_idx]
      
      r_mean     <- with(x$repnum, sum(avg * n, na.rm = TRUE)/sum(n, na.rm = TRUE))
      
      gentime <- with(x$gentime, sum(avg * n, na.rm = TRUE)/sum(n, na.rm = TRUE))
      
      final_preval <- with(x$history, tail(counts[state == "Removed"], 1))
      
      data.table(
        simfile           = fn,
        simid             = x$simid,
        peak_time         = peak_time,
        peak_preval       = peak_preval,
        rt                = rt,
        rt_mean           = r_mean,
        gentime           = gentime,
        final_preval      = final_preval
      )
    }, error = function(e) return(NULL))
  }
  
  simres_list <- mclapply(simfiles, process_sim, mc.cores = ncores)
  simres <- rbindlist(simres_list)
  
  # Save sim results just in case
  fwrite(simres, fn_sim, compress = "auto")
}

simres[, netid := gsub(".+/([0-9]+-[a-z0-9]+)-sim.+", "\\1", simfile)]

# 3. Merge
message("Merging datasets...")
final_data <- merge(simres, net_stats, by = "netid", all.x = TRUE)

# Save
fwrite(final_data, fn_out, compress = "auto")
message("Done! Saved to ", fn_out)
