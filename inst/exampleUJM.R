library(survival)
data(dataLong)
data(dataSurv)

\dontrun{

A <- UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

A$PMean$alpha

A <- UJM(
  formFixed = Y2 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1,
  dataLong, dataSurv, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

A$PMean$alpha





data(dataLong_ind)
data(dataSurv_ind)


B <- UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1+ w2 + x2,
  dataLong_ind, dataSurv_ind, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

B$PMean$alpha



B <- UJM(
  formFixed = Y2 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, CR) ~ w1 + x1+ w2 + x2,
  dataLong_ind, dataSurv_ind, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)

B$PMean$alpha


}
