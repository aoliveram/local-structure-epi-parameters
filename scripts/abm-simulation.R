library(epiworldR)
library(slurmR)

# Set the seed for reproducibility
set.seed(1231)

# Load the simulated networks from the RDS file
networks <- readRDS("data/networks.rds")

# Set the parameters
nsims <- 1000L
Njobs <- 20L

# Set the slurmR options
SB_OPTS <- list(
    account = "vegayon-np-shared",
    partition = "vegayon-np-shared",
    job_name = "abm-simulation"
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
opts_slurmR$set_tmp_path("/tmp")

# Group infectiousness and recovery_rate into a list of length
# nsims featuring one element per simulation
params <- Map(\(a, b, n, i) {
    c(infectiousness = a, recovery_rate = b, net = n, inc_days = i)
}, a = infectiousness, b = recovery_rate, n = networks, i = incubation_days)

res <- Slurm_lapply(params, FUN = \(param) {

    # Creating the SEIR model
    model <- ModelSEIR(
        name = "SEIR",
        prevalence = 0.01,
        infectiousness = param$infectiousness,
        incubation_days = param$inc_days, recovery = param$recovery_rate
    )

    # Adding the network to model
    agents_from_edgelist(
        model,
        source = params$net$source,
        target = params$net$target
        )

    # Running the simulation
    run(model, ndays = 100)

    # Get the results
    list(
        history = get_hist_total(model),
        repnum = get_reproductive_number(model)
    )

}, njobs = Njobs, sbatch_opts = SB_OPTS)

