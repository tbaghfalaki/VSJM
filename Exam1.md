Getting Started
---------------

```
library(VSJM)
```
Loading the data from the package includes both longitudinal data in long format and survival data. It's essential to ensure that the same subject (ID) is present in both datasets.

```
data(dataLong)
data(dataSurv)
```

Dividing data to 50% training data and 50% validation set:

```
set.seed(2)
INDTRAIN <- sample(dataSurv$id, .5 * (dim(dataSurv)[1]))
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
```

We are considering ten markers; therefore, we require three lists as follows: one for the fixed effects model, one for the random effects model, and another for the survival model.

```
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

```

We need to choose the model for the marker trend among "intercept," "linear," and "quadratic." For instance, if we consider a covariate $x_1$, the options are as follows:

intercept:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_1+b_{0ki}+\varepsilon_{ikt}$


linear:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_1+b_{0ki}+b_{1ki} t+\varepsilon_{ikt}$


quadratic:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}t^2+\beta_{3k}x_1+b_{0ki}+b_{1ki} t+b_{1ki} t^2+\varepsilon_{ikt}$

For example, for three markers, we can consider:
```
model <- list("intercept", "linear", "quadratic")
```

In our example, we consider all of them to be linear, as follows:

```
model <- list(
  "linear", "linear", "linear", "linear", "linear", "linear",
  "linear", "linear", "linear", "linear"
)
```


Finally, we have to use the VS function with the following arguments:

-  formFixed a list of formulas for fixed part of longitudinal model
-  formRandom a list of formulas for random part of longitudinal model
-  formGroup a list of formulas specifying the cluster variable for Y (e.g. = list (~ subject, ~ subject,...))
-  formSurv formula for survival model
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  nmark the number of longitudinal markers
-  K1 Number of nodes and weights for calculating Gaussian quadrature in the first stage.
-  K2 Number of nodes and weights for calculating Gaussian quadrature in the second stage.
-  model a list of the models for the longitudinal part which includes "linear" or "quadratic".
-  Obstime the observed time in longitudinal data
-  ncl the number of nodes to be forked for parallel computing
-  Method the method for variable selection including "CS" for continues spike and "DS" for Dirac spike.
-  n.chains1 the number of parallel chains for the model in the first stage; default is 1.
-  n.iter1 integer specifying the total number of iterations in the first stage; default is 1000.
-  n.burnin1 integer specifying how many of n.iter to discard as burn-in in the first stage; default is 5000.
-  n.thin1 integer specifying the thinning of the chains in the first stage; default is 1.
-  n.chains2 the number of parallel chains for the model in the second stage; default is 1.
-  n.iter2 integer specifying the total number of iterations in the second stage; default is 1000.
-  n.burnin2 integer specifying how many of n.iter to discard as burn-in in the second stage; default is 5000.
-  n.thin2 integer specifying the thinning of the chains in the second stage; default is 1.
-  simplify Logical; the option for simplifying the use of CS and DS; default is TRUE.
-  DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
-  quiet Logical, whether to suppress stdout in jags.model().

-----------------


$\gamma_{lj} | \pi_{\gamma_{lj}},\sigma_{\gamma_{lj}}^2 \sim (1 - \pi_{\gamma_{lj}}) N(0,\sigma_{\gamma_{lj}}^2)+ \pi_{\gamma_{lj}}{\delta _0}(\gamma _{lj}),\sigma{\gamma{lj}}^2>0, l=1,...,L, j=1,...,p^\gamma_l,$  

$\pi_{\gamma_{lj}} |{a_{\gamma_{lj}}},{b_{\gamma_{lj}}} \sim Beta(a_{\gamma_{lj}},{b_{\gamma_{lj}}}),$

$\sigma^2_{\gamma _{jk}}$ | 

$a_{\sigma_{\gamma _{jk}}^2},b_{\sigma_{\gamma _{jk}}^2}$


$\sigma_{\gamma_{jk}}^2$




$\sim I\Gamma$ 


$(a_{\sigma_{\gamma _{jk}}^2},b_{\sigma_{\gamma _{jk}}^2})$




As an example, consider the following command, where this implementation has been performed on training data. We consider "DS" method.



```
VS <- VS(formFixed, formRandom, formGroup, formSurv,
  nmark = 10, K1 = 15, K2 = 15,
  model = model, n.chains1 = 2, n.iter1 = 2000, n.burnin1 = 1000,
  n.thin1 = 1, n.chains2 = 2, n.iter2 = 2000, n.burnin2 = 1000,
  n.thin2 = 1, simplify = TRUE, Obstime = "obstime", Method = "DS", ncl = 6,
  DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)
```

The outputs of this function is as follows: 

```

```


If we consider n.chains1 or n.chains2 > 1, the values of the Gelman-Rubin criteria are also provided, which helps in checking the convergence of the MCMC.
