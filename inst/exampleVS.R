rm(list=ls())
data(dataLong)
data(dataSurv)
set.seed(2)
INDTRAIN <- sample(dataSurv$id, 0.5 * (dim(dataSurv)[1]))
INDVALID <- dataSurv$id[-INDTRAIN]
dataLong_t <- subset(
  dataLong,
  dataLong$id %in% INDTRAIN
)
dataSurv_t <- subset(
  dataSurv,
  dataSurv$id %in% INDTRAIN
)
names(dataSurv_t)

dataLong_v <- subset(
  dataLong,
  dataLong$id %in% INDVALID
)
dataSurv_v <- subset(
  dataSurv,
  dataSurv$id %in% INDVALID
)



formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, CR) ~ w1 + x1
model <- list("intercept", "linear", "quadratic")



VS <- VS(formFixed, formRandom, formGroup, formSurv,
  nmark = 3, K1 = 15, K2 = 15,
  model = model, n.chains1 = 1, n.iter1 = 30, n.burnin1 = 10,
  n.thin1 = 1,  n.chains2 = 1, n.iter2 = 30, n.burnin2 = 10,
  n.thin2 = 1, simplify=TRUE, Obstime = "obstime", Method = "DS", ncl = 2,
  DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)

SVS <- SVS(VS)

Step2 <- VS2(VS,
  Method = "LBFDR", n.chains = 2, n.iter = 100, n.burnin = 50,
  n.thin = 1, dataLong = dataLong_t, dataSurv = dataSurv_t
)

DP <- DP(VS, Step2,
  Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 100, n.burnin = 50,
  n.thin = 1,cause_main=1,
  DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)



MCDP <- MCDP(VS, Step2,
             Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 20, n.burnin = 10,
             n.thin = 1,cause_main=1, mi=5,
             DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)
######################################################################
\dontrun{
library(DPCri);library(survival)

Criteria(s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
         CR = dataSurv_v$CR, P = DP$DP[,2], cause = 1)



Criteria(s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
         CR = dataSurv_v$CR, P = MCDP$DP[,2], cause = 1)

######################################################################
rm(list=ls())
data(dataLong_m)
data(dataSurv_m)
set.seed(2)
INDTRAIN <- sample(dataSurv_m$id, 0.5 * (dim(dataSurv_m)[1]))
INDVALID <- dataSurv_m$id[-INDTRAIN]
dataLong_t <- subset(
  dataLong_m,
  dataLong_m$id %in% INDTRAIN
)
dataSurv_t <- subset(
  dataSurv_m,
  dataSurv_m$id %in% INDTRAIN
)
names(dataSurv_t)

dataLong_v <- subset(
  dataLong_m,
  dataLong_m$id %in% INDVALID
)
dataSurv_v <- subset(
  dataSurv_m,
  dataSurv_m$id %in% INDVALID
)



formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, CR) ~ w1 + x1
model <- list("intercept", "linear", "quadratic")



VS <- VS(formFixed, formRandom, formGroup, formSurv,
         nmark = 2, K1 = 15, K2 = 15,
         model = model, n.chains1 = 2, n.iter1 = 30, n.burnin1 = 10,
         n.thin1 = 1,  n.chains2 = 2, n.iter2 = 30, n.burnin2 = 10,
         n.thin2 = 1, simplify=TRUE, Obstime = "obstime", Method = "DS", ncl = 2,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)

SVS <- SVS(VS)

Step2 <- VS2(VS,
             Method = "LBFDR", n.chains = 2, n.iter = 30, n.burnin = 20,
             n.thin = 1, dataLong = dataLong_t, dataSurv = dataSurv_t
)

DP <- DP(VS, Step2,
         Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 30, n.burnin = 20,
         n.thin = 1,cause_main=1,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)


Criteria(s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
         CR = dataSurv_v$CR, P =DP$DP$est, cause = 1)


MCDP <- MCDP(VS, Step2,
         Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 20, n.burnin = 10,
         n.thin = 1,cause_main=1, mi=5,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)

Criteria(s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
         CR = dataSurv_v$CR, P = MCDP$DP$est, cause = 1)

######################################################################

  rm(list=ls())
  set.seed(2)
  INDTRAIN=sample(dataSurv$id,.5*(dim(dataSurv)[1]))
  INDVALID <- dataSurv$id[-INDTRAIN]

  dataLong_t=subset(dataLong,
                    dataLong$id %in% INDTRAIN)

  dataSurv_t=subset(dataSurv,
                    dataSurv$id %in% INDTRAIN)
  names(dataSurv_t)
  dataLong_v <- subset(
    dataLong,
    dataLong$id %in% INDVALID
  )
  dataSurv_v <- subset(
    dataSurv,
    dataSurv$id %in% INDVALID
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
  formSurv <- survival::Surv(survtime, CR) ~ w1 + x1
  model <- list(
    "linear", "linear", "linear", "linear", "linear", "linear",
    "linear", "linear", "quadratic", "quadratic"
  )




  VS <- VS(formFixed, formRandom, formGroup, formSurv,
           nmark = 10, K1 = 15, K2 = 15,
           model = model, n.chains1 = 2, n.iter1 = 20, n.burnin1 = 10,
           n.thin1 = 1,  n.chains2 = 2, n.iter2 = 20, n.burnin2 = 10,
           n.thin2 = 1, simplify=TRUE, Obstime = "obstime", Method = "DS", ncl = 2,
           DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
  )

  SVS <- SVS(VS)

  Step2 <- VS2(VS,
               Method = "LBFDR", n.chains = 2, n.iter = 20, n.burnin = 10,
               n.thin = 1, dataLong = dataLong_t, dataSurv = dataSurv_t
  )

  DP <- DP(VS, Step2,
           Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 30, n.burnin = 20,
           n.thin = 1,cause_main=1,
           DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
  )


  Criteria(s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
           CR = dataSurv_v$CR, P =DP$DP$est, cause = 1)


  MCDP <- MCDP(VS, Step2,
               Method = "LBFDR", s = 0, t = 0.5, n.chains = 1, n.iter = 20, n.burnin = 10,
               n.thin = 1,cause_main=1, mi=5,
               DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
  )

  Criteria(s = 0, t = 0.5, Survt = dataSurv_v$survtime,
           CR = dataSurv_v$CR, P = MCDP$DP$est, cause = 1)





}




