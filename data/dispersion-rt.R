library(epiworldR)

model <- ModelSIRCONN(
  name = "A", n = 1e3, prevalence = 1/1e3, contact_rate = 2,
  transmission_rate = .3, recovery_rate = 1/3
  )

set.seed(1233)
run_multiple(
  model, 100, nsims = 500, saver = make_saver("reproductive"),
  nthreads = 6
  )

plot(model)

repnums <- run_multiple_get_results(model)[[1]]

subset(repnums, source_exposure_date == 0)$rt |>
  hist()
