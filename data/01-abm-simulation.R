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

# Set the seed for reproducibility
set.seed(1231)

# Load the simulated networks from the RDS file
networks <- readRDS("data/Simulated_1000_networks.rds")

# Turning the networks to edgelists
elists <- lapply(networks, \(n) {
  network::as.edgelist(n)
})

# Set the parameters
nsims <- 5000
Njobs <- 20L

# Set the slurmR options
SB_OPTS <- list(
    account    = "vegayon-np",
    partition  = "vegayon-shared-np",
    mem        = "32G"
    )

# Sampling infectiousness from a beta distribution
alpha <- 2
beta <- 5
infectiousness <- rbeta(nsims, alpha, beta)

# Sampling incubation days from a beta distribution
alpha <- 2
beta <- 5
incubation_days <- rbeta(nsims, alpha, beta)

# Sampling recovery rate from a beta distribution
alpha <- 5
beta <- 2
recovery_rate <- rbeta(nsims, alpha, beta)

# Set the temporary path of slurmR
# opts_slurmR$set_tmp_path("/tmp")

# Group infectiousness and recovery_rate into a list of length
# nsims featuring one element per simulation
net_ids <- sample.int(length(elists), nsims, replace = TRUE)
seeds   <- sample.int(.Machine$integer.max, nsims, replace = TRUE)

params <- Map(\(a, b, n, i, netid) {
  list(
    netid = netid, infectiousness = a, recovery_rate = b,
    net = n, inc_days = i, seed = seeds[netid],
    size = network::network.size(networks[[netid]])
    )
  },
  a = infectiousness,
  b = recovery_rate,
  i = incubation_days,
  n = elists[net_ids],
  netid = net_ids
  )

# Sampling

# opts_slurmR$debug_on()

res <- Slurm_lapply(params, FUN = \(param) {

  tryCatch({
    # Creating the SEIR model
    model <- ModelSEIR(
      name = "SEIR",
      prevalence = 0.01,
      infectiousness = param$infectiousness,
      incubation_days = param$inc_days,
      recovery = param$recovery_rate
    )

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
      netid     = param$netid,
      history   = get_hist_total(model),
      repnum    = plot_reproductive_number(model, plot = FALSE),
      incidence = plot_incidence(model, plot = FALSE),
      gentime   = plot_generation_time(model, plot = FALSE)
    )
  }, error = function(e) e)

}, njobs = Njobs, sbatch_opt = SB_OPTS, job_name = "abm-simulation-lapply")


# Saving the results under data/
saveRDS(res, file = "data/01-abm-simulation.rds")

