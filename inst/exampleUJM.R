library(survival)
data(dataLong)
data(dataSurv)

\dontrun{

A=UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "linear", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "quadratic", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "linear", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "quadratic", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

}
