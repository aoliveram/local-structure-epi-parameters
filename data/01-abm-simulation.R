#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=abm-simulation-main
#SBATCH --output=abm-simulation-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --mem=32G

library(epiworldR)
library(slurmR)
library(network)

# Set the seed for reproducibility
set.seed(1231)

# Load the simulated networks from the RDS file
networks <- readRDS("data/Simulated_1000_networks.rds")

# Combining with the sf, sw, and r networks
networks_sf <- readRDS("data/Simulated_1000_networks_sf.rds")
networks_sw <- readRDS("data/Simulated_1000_networks_sw.rds")
networks_r <- readRDS("data/Simulated_1000_networks_r.rds")

# Turning the networks to edgelists
elists <- lapply(networks, \(n) {
  network::as.edgelist(n)
})

# Combining with the sf, sw, and r networks
as_edgelist_int <- function(x) {
  matrix(as.integer(
    igraph::as_edgelist(x$net)
  ), ncol = 2)
}
elists <- c(
  elists,
  lapply(networks_sf, as_edgelist_int),
  lapply(networks_sw, as_edgelist_int),
  lapply(networks_r, as_edgelist_int)
  )

# Vector with sizes
sizes <- c(
  lapply(networks, \(x) {network::network.size(x)}),
  lapply(networks_sf, \(x) {x$n}),
  lapply(networks_sw, \(x) {x$n}),
  lapply(networks_r, \(x) {x$n})
)

# Vector with network types
net_types <- rep(c("ergm", "sf", "sw", "r"), each = length(networks))

# Set the parameters
nsims <- 20000
Njobs <- 20L

# Set the slurmR options
SB_OPTS <- list(
    account    = "vegayon-np",
    partition  = "vegayon-shared-np",
    mem        = "64G"
    )

# Sampling infectiousness from a beta distribution
# This has mean 0.3 and sd 0.05
alpha <- 20
beta <- 100
infectiousness <- rbeta(nsims, alpha, beta)

# Sampling incubation days from a Gamma distribution
# This has mean 7 and variance 7
alpha <- 7 * 2
beta <- 1 * 2
incubation_days <- ceiling(rgamma(nsims, alpha, beta))

# Sampling recovery rate from a beta distribution
# This has mean 0.3 and variance 0.02
alpha <- 20
beta <- 50
recovery_rate <- rbeta(nsims, alpha, beta)

# Set the temporary path of slurmR
# opts_slurmR$set_tmp_path("/tmp")

# Group infectiousness and recovery_rate into a list of length
# nsims featuring one element per simulation
net_ids <- sample.int(length(elists), nsims, replace = TRUE)
seeds   <- sample.int(.Machine$integer.max, nsims, replace = TRUE)

params <- Map(\(a, b, n, i, netid, simid) {
  list(
    netid = netid, infectiousness = a, recovery_rate = b,
    net   = n, inc_days = i, seed = seeds[simid],
    size  = sizes[[netid]], #network::network.size(networks[[netid]]),
    simid = simid
    )
  },
  a = infectiousness,
  b = recovery_rate,
  i = incubation_days,
  n = elists[net_ids],
  netid = net_ids,
  simid = 1:nsims
  )

# Sampling

# opts_slurmR$debug_on()

res <- Slurm_lapply(params, FUN = \(param) {

  message(
    "Initiating model for netid ", param$netid, " and simid ", param$simid,
    "...", appendLF = FALSE
    )

  ans <- tryCatch({
    # Creating the SEIR model
    model <- ModelSEIR(
      name = "SEIR",
      prevalence = 0.01,
      transmission_rate = param$infectiousness,
      incubation_days = param$inc_days,
      recovery_rate = param$recovery_rate
    )

    verbose_off(model)

    # Adding the network to model
    agents_from_edgelist(
      model,
      size   = param$size,
      source = param$net[, 1] - 1L,
      target = param$net[, 2] - 1L,
      directed = FALSE
      )

    # Running the simulation
    run(model, ndays = 100, seed = param$seed)

    # Get the results
    list(
      simid     = param$simid,
      netid     = param$netid,
      history   = get_hist_total(model),
      repnum    = plot_reproductive_number(model, plot = FALSE),
      incidence = plot_incidence(model, plot = FALSE),
      gentime   = plot_generation_time(model, plot = FALSE),
      params    = param
    )
  }, error = function(e) e)

  if (inherits(ans, "error")) {
    message(
      "Error for netid ", param$netid, " and simid ", param$simid,
      appendLF = FALSE
      )
    return(ans)
  } 

  message("done")

  ans

}, njobs = Njobs, sbatch_opt = SB_OPTS, job_name = "abm-simulation-lapply",
plan = "wait")

print(res)

# Saving the results under data/
saveRDS(
  Slurm_collect(res, any. = TRUE), file = "data/01-abm-simulation.rds", 
  compress = FALSE
  )

