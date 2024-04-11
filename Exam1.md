Getting Started
---------------

```
library(TSJM)
```
Generating summary statistics


Loading data of the package:

```
data(dataLong)
data(dataSurv)
```

Dividing data to 80% training data and 20% validation set


```
set.seed(2)
INDTRAIN <- sample(dataSurv$id, 0.8 * (dim(dataSurv)[1]))
INDVALID <- dataSurv$id[-INDTRAIN]
dataLong_t <- subset(
  dataLong,
  dataLong$id %in% INDTRAIN
)
dataSurv_t <- subset(
  dataSurv,
  dataSurv$id %in% INDTRAIN
)

dataLong_v <- subset(
  dataLong,
  dataLong$id %in% INDVALID
)
dataSurv_v <- subset(
  dataSurv,
  dataSurv$id %in% INDVALID
)
```

We are considering three markers; therefore, we require three lists as follows: one for the fixed effects model, one for the random effects model, and another for the survival model.

```
formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, death) ~ x1 + x2
```

We need to choose the model for the marker trend among "intercept," "linear," and "quadratic." For instance, if we consider a covariate $x_1$, the options are as follows:

intercept:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_1+b_{0ki}+\varepsilon_{ikt}$


linear:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_1+b_{0ki}+b_{1ki} t+\varepsilon_{ikt}$


quadratic:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}t^2+\beta_{3k}x_1+b_{0ki}+b_{1ki} t+b_{1ki} t^2+\varepsilon_{ikt}$


This has been done by considering the following command:

```
model <- list("intercept", "linear", "quadratic")
```
Finally, we have to use the TS function with the following arguments:


- formFixed a list of formulas for fixed part of longitudinal model
- formRandom a list of formulas for random part of longitudinal model
- formGroup a list of formulas specifying the cluster variable for Y (e.g. = list (~ subject, ~ subject,...))
- formSurv formula for survival model
- dataLong data set of observed longitudinal variables.
- dataSurv data set of observed survival variables.
- nmark the number of longitudinal markers
- K1 Number of nodes and weights for calculating Gaussian quadrature in the first stage.
- K2 Number of nodes and weights for calculating Gaussian quadrature in the second stage.
- model a list of the models for the longitudinal part which includes "linear" or "quadratic".
- Obstime the observed time in longitudinal data
- ncl the number of nodes to be forked for parallel computing
- n.chains1 the number of parallel chains for the model in the first stage; default is 1.
- n.iter1 integer specifying the total number of iterations in the first stage; default is 1000.
- n.burnin1 integer specifying how many of n.iter to discard as burn-in in the first stage; default is 5000.
- n.thin1 integer specifying the thinning of the chains in the first stage; default is 1.
- n.chains2 the number of parallel chains for the model in the second stage; default is 1.
- n.iter2 integer specifying the total number of iterations in the second stage; default is 1000.
- n.burnin2 integer specifying how many of n.iter to discard as burn-in in the second stage; default is 5000.
- n.thin2 integer specifying the thinning of the chains in the second stage; default is 1.
- simplify Logical; the option for simplifying the use of CS and DS; default is TRUE.
- DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
- quiet Logical, whether to suppress stdout in jags.model().
-----------------

As an example, consider the following command, where this implementation has been performed on training data:


```
TS0 <- TS(formFixed, formRandom, formGroup, formSurv,
         nmark = 3, K1 = 15, K2 = 15,
         model = model, n.chains1 = 1, n.iter1 = 2000, n.burnin1 = 1000,
         n.thin1 = 1,  n.chains2 = 1, n.iter2 = 3000, n.burnin2 = 1000,
         n.thin2 = 1, Obstime = "obstime", ncl = 3,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)
```

The outputs of this function is as follows: 

```
> TS0$Longitudinal
[[1]]
[[1]]$Longitudinal_model
                     Est          SD       L_CI      U_CI
(Intercept) -0.003167629 0.051649161 -0.1016439 0.0865190
obstime      0.651195083 0.024117669  0.6057510 0.6996882
sigma2_e     0.304434795 0.006494876  0.2912420 0.3178109
sigma2_b     1.318007024 0.087255212  1.1564923 1.5036663


[[2]]
[[2]]$Longitudinal_model
                  Est         SD       L_CI      U_CI
(Intercept) 0.1079177 0.02551374 0.05766939 0.1556610
obstime     0.7318347 0.02821918 0.68070271 0.7879759
sigma2      0.1988282 0.00442706 0.19091613 0.2079220

[[2]]$Sigma
          Intercept      Time
Intercept 1.0611350 0.5255265
Time      0.5255265 1.0764285


[[3]]
[[3]]$Longitudinal_model
                    Est          SD        L_CI        U_CI
(Intercept) -0.01248891 0.020887352 -0.05146358  0.02923079
obstime      0.53419759 0.073747348  0.39215866  0.64446072
obstime2    -0.63622137 0.052401701 -0.70784160 -0.53403394
sigma2       0.19828364 0.004728854  0.18956590  0.20868304

[[3]]$Sigma
           Intercept       Time      Time2
Intercept  1.0512538  0.5376876 -0.1235974
Time       0.5376876  1.0998180 -0.2271203
Time2     -0.1235974 -0.2271203  0.4990175


> TS0$TDsurvival
$S_model
                Est         SD        L_CI        U_CI
x1       0.27251737 0.14241424 -0.01891624  0.54445881
x2       0.02010240 0.13661497 -0.25578120  0.27991697
Marker1 -0.03986755 0.05585649 -0.14893049  0.07103536
Marker2 -0.27965352 0.05428684 -0.38340941 -0.17559974
Marker3  0.21850350 0.05160108  0.11619111  0.32225479
h1       1.05664630 0.19742373  0.71824690  1.49755694
h2       0.90545060 0.16341794  0.61461051  1.24709989
h3       1.00960526 0.18004706  0.69618933  1.39670388
h4       0.97249676 0.17158193  0.67793851  1.33887908
h5       0.95898814 0.18192923  0.64538832  1.35723983
```


If we consider n.chains1 or n.chains2 > 1, the values of the Gelman-Rubin criteria are also provided, which helps in checking the convergence of the MCMC.
