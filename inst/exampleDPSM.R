\dontrun{
rm(list = ls())
data(dataLong_ind)
data(dataSurv_ind)
set.seed(2)
INDTRAIN <- sample(dataSurv_ind$id, 0.8 * (dim(dataSurv_ind)[1]))
INDVALID <- dataSurv_ind$id[-INDTRAIN]
dataLong_t <- subset(
  dataLong_ind,
  dataLong_ind$id %in% INDTRAIN
)
dataSurv_t <- subset(
  dataSurv_ind,
  dataSurv_ind$id %in% INDTRAIN
)
names(dataSurv_t)

dataLong_v <- subset(
  dataLong_ind,
  dataLong_ind$id %in% INDVALID
)
dataSurv_v <- subset(
  dataSurv_ind,
  dataSurv_ind$id %in% INDVALID
)


formFixed <- list(
  Y1 ~ obstime + x1 + x2, Y2 ~ obstime + x1 + x2, Y3 ~ obstime + x1 + x2,
  Y4 ~ obstime + x1 + x2, Y5 ~ obstime + x1 + x2, Y6 ~ obstime + x1 + x2,
  Y7 ~ obstime + x1 + x2, Y8 ~ obstime + x1 + x2, Y9 ~ obstime + x1 + x2,
  Y10 ~ obstime + x1 + x2
)
formRandom <- list(
  ~obstime, ~obstime, ~obstime, ~obstime, ~obstime, ~obstime,
  ~obstime, ~obstime, ~obstime, ~obstime
)
formGroup <- list(~id, ~id, ~id, ~id, ~id, ~id, ~id, ~id, ~id, ~id)
formSurv <- survival::Surv(survtime, CR) ~ w1 + x1 + w2 + x2
model <- list(
  "linear", "linear", "linear", "linear", "linear", "linear",
  "linear", "linear", "quadratic", "quadratic"
)


VS <- VS(formFixed, formRandom, formGroup, formSurv,
  nmark = 2, K1 = 15, K2 = 15,
  model = model, n.chains1 = 2, n.iter1 = 2000, n.burnin1 = 1000,
  n.thin1 = 1, n.chains2 = 2, n.iter2 = 2000, n.burnin2 = 1000,
  n.thin2 = 1, simplify = TRUE, Obstime = "obstime", Method = "DS",
  ncl = 2,
  DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)

VS$sim_step1[[1]]$PMean$alpha
VS$sim_step1[[2]]$PMean$alpha

SVS(VS)


Step2 <- VS2(VS,
  Method = "LBFDR", n.chains = 2, n.iter = 1000, n.burnin = 500,
  n.thin = 1, dataLong = dataLong_t, dataSurv = dataSurv_t
)





D1 <- DPSM(VS, Step2,
  N_marker = c(1, 2), s = 0.1, t = 0.5, cause_main = 1, n.chains = 1,
  n.iter = 2000, n.burnin = 1000,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE, dataLong_v, dataSurv_v
)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = D1$DP$est, cause = 1
)$Cri


D2 <- DPOM(VS,
  N_marker = 1, s = 0.1, t = 0.5, cause_main = 1, n.chains = 1,
  n.iter = 1000, n.burnin = 500,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE, dataLong_v, dataSurv_v
)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = D2$DP$est, cause = 1
)$Cri

}
