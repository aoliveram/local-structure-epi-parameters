#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --ntasks=31
#SBATCH --mem=128GB
#SBATCH --job-name=02-dataprep
#SBATCH --time=24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=george.vegayon@utah.edu

library(igraph)
library(network)
# library(netplot)
library(intergraph)
library(ergm)
library(data.table)

ncores <- 30

# Read the Simulated data rds file from data/
message("Loading the networks files")
networks <- list.files(
  "data/graphs", pattern = "[0-9]+-(ergm|sf|sw|degseq)\\.rds$",full.names = TRUE) 

# File name for each type
fn_ergm <- "data/02-dataprep-ergm.csv.gz"
fn_igraph <- "data/02-dataprep-igraph.csv.gz"
fn_sim <- "data/02-dataprep-sim.csv.gz"

# Computing statistics using ERGM
if (!file.exists(fn_ergm)) {
  message("Computing statistics using ERGM")
  S_ergm <- parallel::mclapply(seq_along(networks), \(i) {

    n <- networks[[i]]

    # Reading the network
    n <- readRDS(n)

    res <- summary_formula(
      n ~ edges + nodematch("gender") + nodematch("grade") + triangles + balance +
        twopath + gwdegree(decay = .25, fixed = TRUE) + isolates
    ) |> as.list() |> as.data.table()

    res[, netfile := networks[[i]]]

    res
  }, mc.cores = ncores) |> rbindlist()

  # Pre-appending `ergm_` to all column names, except netfile
  setnames(S_ergm, new = paste0("ergm_", names(S_ergm)))
  setnames(S_ergm, old = "ergm_netfile", new = "netfile")

  fwrite(S_ergm, fn_ergm, compress = "auto")

  message("Done computing statistics using ERGM")
} else {
  message("Loading statistics using ERGM")
  S_ergm <- fread(fn_ergm)
}

head(S_ergm)

# Computing statistics based on igraph ----------------------------------------

if (!file.exists(fn_igraph)) {

  message("Computing statistics based on igraph")
  S_igraph <- parallel::mclapply(networks, \(i) { 

    # Reading the network and converting it into igraph
    n <- readRDS(i)
    n <- intergraph::asIgraph(n)
    
    data.table(
      netfile         = i,
      modularity      = modularity(cluster_fast_greedy(n)),
      transitivity    = transitivity(n),
      density         = igraph::edge_density(n),
      diameter        = diameter(n),
      avg_path_length = igraph::mean_distance(n),
      avg_degree      = mean(degree(n), na.rm = TRUE),
      avg_betweenness = mean(betweenness(n), na.rm = TRUE),
      avg_closeness   = mean(closeness(n), na.rm = TRUE),
      avg_eigenvector = mean(eigen_centrality(n)$vector, na.rm = TRUE),
      components      = components(n)$no
    )
  }, mc.cores = ncores) |> rbindlist()

  head(S_igraph)

  setnames(S_igraph, new = paste0("igraph_", names(S_igraph)))
  setnames(S_igraph, old = "igraph_netfile", new = "netfile")

  fwrite(S_igraph, fn_igraph, compress = "auto")

  message("Done computing statistics based on igraph")
} else {
  message("Loading statistics based on igraph")
  S_igraph <- fread(fn_igraph)
}

# Combining the datasets
S <- merge(S_igraph, S_ergm, by = "netfile", all = TRUE)

S[, nettype := gsub(".+-([a-z]+)\\.rds", "\\1", netfile)]
S[, netid   := gsub(".+/([0-9]+-[a-z]+)\\.rds", "\\1", netfile)]


# Reading simulation results ---------------------------------------------------

if (!file.exists(fn_sim)) {

  message("Processing the simulation results")
  simfiles <- list.files("data/graphs", pattern = "-sim-[0-9]+\\.rds$", full.names = TRUE)

  simres <- parallel::mclapply(simfiles, \(fn) {

    x <- readRDS(fn)
    
    res <- tryCatch({
      history <- data.table(x$history)

      # Index of the max prevalence
      peak_idx    <- which.max(x$incidence$Exposed)
      peak_preval <- x$incidence$Exposed[peak_idx]
      peak_time   <- as.integer(rownames(x$incidence)[peak_idx])
      
      # Rt in the peak
      rt_idx     <- with(x$repnum, which.min(abs(peak_time - date)))
      rt         <- x$repnum$avg[rt_idx]
      dispersion <- 1/x$repnum$sd[rt_idx]^2

      # Mean Rt
      r_mean     <- with(x$repnum, sum(avg * n, na.rm = TRUE)/
        sum(n, na.rm = TRUE))

      #  Rt per day
      rt_0 <- with(x$repnum[x$repnum$date == 0,,drop=FALSE], avg)
      rt_1 <- with(x$repnum[x$repnum$date == 1,,drop=FALSE], avg)
      rt_2 <- with(x$repnum[x$repnum$date == 2,,drop=FALSE], avg)
      rt_3 <- with(x$repnum[x$repnum$date == 3,,drop=FALSE], avg)
      rt_4 <- with(x$repnum[x$repnum$date == 4,,drop=FALSE], avg)
      rt_5 <- with(x$repnum[x$repnum$date == 5,,drop=FALSE], avg)

      gentime <- with(x$gentime, sum(avg * n, na.rm = TRUE)/
        sum(n, na.rm = TRUE))

      # Final prevalence
      final_preval <- with(
        x$history, tail(counts[state == "Removed"], 1)
        )

      # Return a data.table
      data.table(
        simfile           = fn,
        simid             = x$simid,
        peak_time         = peak_time,
        peak_preval       = peak_preval,
        rt                = rt,
        rt_mean           = r_mean,
        rt_0              = rt_0,
        rt_1              = rt_1,
        rt_2              = rt_2,
        rt_3              = rt_3,
        rt_4              = rt_4,
        rt_5              = rt_5,
        dispersion        = dispersion,
        gentime           = gentime,
        final_preval      = final_preval,
        transmission_rate = x$param$transmission_rate,
        inc_days          = x$param$inc_days,
        recovery_rate     = x$param$recovery_rate
      )
    }, error=function(e) e)

    if (inherits(res, "error"))
      return(NULL)

    res

  }, mc.cores = ncores) |> rbindlist()

  fwrite(simres, fn_sim, compress = "auto")

  message("Done processing the simulation results")

} else {

  message("Loading the simulation results")
  simres <- fread(fn_sim)
  
}

simres[, netid := gsub(".+/([0-9]+-[a-z]+)-sim.+", "\\1", simfile)]

print(head(S))

print(head(simres))

# Merge the datasets
S <- merge(simres, S, by = "netid", all.x = TRUE, all.y = FALSE)

fwrite(
  S,
  file = "data/02-dataprep-network-stats.csv.gz"
)

message("Done!")
