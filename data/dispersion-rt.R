library(epiworldR)
library(data.table)

r_rate   <- 1/3
c_rate   <- 5 # Average observed
nsims    <- 100
nthreads <- 20
n        <- 5e4

rts <- c(.75, 1.5, 2, 2.5, 3)
res <- NULL

set.seed(1233)

# Creating filenames
fnames <- sprintf("%s/%%05lu-simulation-%s.csv", tempdir(), rts)
fnames <- setNames(fnames, rts)

for (rt in rts) {

  message("Running simulations for Rt = ", rt)

  # Figuring out the transmission rate
  t_rate <- rt * r_rate/(c_rate - rt + rt * r_rate) # Analytical sol

  model <- ModelSIRCONN(
    name = "A", n = n, prevalence = 1/n, contact_rate = c_rate,
    transmission_rate = t_rate, recovery_rate = r_rate
    )

  run_multiple(
    model, 100, nsims = nsims,
    saver = make_saver("reproductive", fn = fnames[as.character(rt)]),
    nthreads = nthreads
    )

  repnums <- data.table(run_multiple_get_results(model)[[1]])

  # Making sure we clean all the files
  unlink(list.files(tempdir(), full.names = TRUE, pattern = "simulation"))

  repnums <- repnums[source_exposure_date == 0]

  repnums <- c(repnums$rt, rep(0, nsims - nrow(repnums)))

  res <- rbind(
    res, data.table(
    rt     = rt,
    rt_obs = repnums,
    t_rate = t_rate
  ))

  print(res[, mean(rt_obs), by = "rt"])

}

res[, mean(rt_obs), by = .(rt, t_rate)]
