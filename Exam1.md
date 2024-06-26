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
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_{1i}+b_{0ki}+\varepsilon_{ikt}$


linear:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_{1i}+b_{0ki}+b_{1ki} t+\varepsilon_{ikt}$


quadratic:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}t^2+\beta_{3k}x_{1i}+b_{0ki}+b_{1ki} t+b_{1ki} t^2+\varepsilon_{ikt}$

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
This approach considers CS or DS priors as follows:

DS:

$\gamma_{lj} | \pi_{\gamma_{lj}},\sigma_{\gamma_{lj}}^2 \sim (1 - \pi_{\gamma_{lj}}) N(0,\sigma_{\gamma_{lj}}^2)+ \pi_{\gamma_{lj}}{\delta _0}(\gamma _{lj}), l=1,...,L, j=1,...,p^\gamma_l,$  

CS:

$\gamma_{lj} | \pi_{\gamma_{lj}},\sigma_{\gamma_{lj}}^2 \sim (1 - \pi_{\gamma_{lj}}) N(0,\sigma_{\gamma_{lj}}^2)+ \pi_{\gamma_{lj}}N(0,10^{-3}), l=1,...,L, j=1,...,p^\gamma_l,$  


It is important to note that if you consider "simplify=TRUE", instead of hirarchical setup for the variance $\sigma_{\gamma_{lj}}^2$ an  inverse gamma prior with parameters (0.01,0.01) are considered for both CS and DS. 


As an example, consider the following command, where this implementation has been performed on training data using the "DS" method:

```
VS <- VS(formFixed, formRandom, formGroup, formSurv,
  nmark = 10, K1 = 15, K2 = 15,
  model = model, n.chains1 = 2, n.iter1 = 2000, n.burnin1 = 1000,
  n.thin1 = 1, n.chains2 = 2, n.iter2 = 2000, n.burnin2 = 1000,
  n.thin2 = 1, simplify = TRUE, Obstime = "obstime", Method = "DS", ncl = 6,
  DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)
```

For summarizing the outputs of this function, we utilize the following function:

```
SVS <- SVS(VS)
```

The output of the function is as follows:

```
> SVS
$gamma
$gamma$LBFDR
           w1    x1
Cause 1 0.000 0.996
Cause 2 0.453 0.967

$gamma$BF
              w1          x1
Cause 1      Inf 0.004016064
Cause 2 1.207506 0.034126163

$gamma$I_LBFDR
        w1 x1
Cause 1  1  0
Cause 2  0  0

$gamma$I_BF
        w1 x1
Cause 1  1  0
Cause 2  1  0


$alpha
$alpha$LBFDR
        Marker1 Marker2 Marker3 Marker4 Marker5 Marker6 Marker7 Marker8 Marker9 Marker10
Cause 1       0       0       1       1   0.143       1  0.1365       1       1        1
Cause 2       1       0       0       1   1.000       1  1.0000       1       1        1

$alpha$BF
        Marker1 Marker2 Marker3 Marker4  Marker5 Marker6  Marker7 Marker8 Marker9 Marker10
Cause 1     Inf     Inf       0       0 5.993007       0 6.326007       0       0        0
Cause 2       0     Inf     Inf       0 0.000000       0 0.000000       0       0        0

$alpha$I_LBFDR
        Marker1 Marker2 Marker3 Marker4 Marker5 Marker6 Marker7 Marker8 Marker9 Marker10
Cause 1       1       1       0       0       0       0       0       0       0        0
Cause 2       0       1       1       0       0       0       0       0       0        0

$alpha$I_BF
        Marker1 Marker2 Marker3 Marker4 Marker5 Marker6 Marker7 Marker8 Marker9 Marker10
Cause 1       1       1       0       0       1       0       1       0       0        0
Cause 2       0       1       1       0       0       0       0       0       0        0
```
At this stage variable selection has been done. The next stage is risk prediction. 

Dynamic prediction
---------------
To reduce estimation biases resulting from variable selection, we propose incorporating an additional stage to calculate dynamic predictions. After variable selection using CS or DS prior, we recommend re-estimating the proportional hazard model by substituting CS or Ds with non-informative normal priors for the association parameters of the selected markers and the regression coefficients of the selected covariates. This has been done by considering *VS2* function in the package. The main arguments in this function are:

- object an object inheriting from class VS function.
- Method the method for variable selection including "LBFDR" for LBFDR and "BF" for Bayes factor.

The following command is considered for this aim:

```
Step2 <- VS2(VS,
  Method = "LBFDR", n.chains = 2, n.iter = 2000, n.burnin = 1000,
  n.thin = 1, dataLong = dataLong_t, dataSurv = dataSurv_t
)
```
Finally, for dynamic prediction, we should utilize the *DP* function, specifying the following arguments:


- object an object inheriting from class VS
- object2 an object inheriting from class VS2
- Method the method for variable selection including "LBFDR" for LBFDR and "BF" for Bayes factor.
- s the landmark time for prediction
- t the window of prediction for prediction
- cause_main the main cause for prediction
- dataLong data set of observed longitudinal variables (validation set).
- dataSurv data set of observed survival variables (validation set).



```
DP <- DP(VS, Step2,
  Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 3000, n.burnin = 2000,
  n.thin = 1, cause_main = 1,
  DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)
```

The outputs of this function are as follows:

```
$DP
     id         est
1     2 0.710044463
2     5 0.121999798
3     7 0.155945232
4    10 0.059012615
5    11 0.031432560
6    15 0.059668201
7    16 0.029397494
8    18 0.104735078
9    19 0.134540067
10   21 0.223344760
11   22 0.405120681
12   24 0.009959658
.
.
.
243 488 0.054126771
244 492 0.051628731
245 493 0.089307181
246 494 0.007899571
247 495 0.058865561
248 497 0.125211111
249 498 0.045489998
250 500 0.054972112

$s
[1] 0.1

$t
[1] 0.5
```
Computing AUC and BS for the predictions
---------------
For this purpose, we use DPCri package <https://github.com/tbaghfalaki/DPCri>.

Computing the criteria using this package is straightforward, as demonstrated by the following commands:

- s the landmark time for prediction
- t the window of prediction for prediction
- Survt the survival time
- CR the indicator for competing risks or censoring
- P the risk predictions
- cause the main cause for prediction


Consider the following command: 

```
library(survival)
library(DPCri)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = DP$DP$est, cause = 1
)
```
with the following outputs:

```
$Cri
           est         sd
AUC 0.74090038 0.05008035
BS  0.09049356 0.01382442
```

Other beneficial functions
---------

Although all the requirements for computing risk prediction are completed at this stage, we propose some other beneficial functions.


- The first one is *DP0*, which computes dynamic prediction without reestimation and the use of *VS2* functions.

```
DP0 <- DP0(VS, s = 0.1, t = 0.5, n.chains = 1, n.iter = 3000, n.burnin = 2000,
         n.thin = 1, cause_main = 1,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)

$Cri
          est         sd
AUC 0.7344633 0.04313503
BS  0.1251009 0.01484688
```

- The second one is dynamic prediction for one-marker by *DPOM* as follows:

```
D1 <- DPOM(VS,
  N_marker = 1, s = 0.1, t = 0.5, cause_main = 1, n.chains = 1,
  n.iter = 2000, n.burnin = 1000,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE, dataLong_v, dataSurv_v
)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = D1$DP$est, cause = 1
)$Cri
```
with the following outputs:

```
$DP
     id        est
1     2 0.35079126
2     5 0.13189275
3     7 0.08357039
4    10 0.12484472
5    11 0.14788052
.
.
.
247 495 0.12207898
248 497 0.08566377
249 498 0.09357884
250 500 0.08983822

$s
[1] 0.1

$t
[1] 0.5


    est         sd
AUC 0.68654215 0.05618483
BS  0.09340079 0.01390426
```

- The third one is dynamic prediction for some markers JM by *DPSM* as follows which markers 1 and 2 are considered:
 
```
D1 <- DPSM(VS, Step2,
  N_markers = c(1, 2), s = 0.1, t = 0.5, cause_main = 1, n.chains = 1,
  n.iter = 2000, n.burnin = 1000,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE, dataLong_v, dataSurv_v
)
```
with the following outputs:

```
$DP
     id         est
1     2 0.636124878
2     5 0.108507783
3     7 0.133040704
4    10 0.062460364
5    11 0.030595321
6    15 0.054849917
7    16 0.029366528
8    18 0.102721955
9    19 0.126131630
10   21 0.216189878
11   22 0.370319961
12   24 0.010574606
.
.
.
247 495 0.054092313
248 497 0.121026926
249 498 0.047753484
250 500 0.053321847

     est         sd
AUC 0.74880268 0.04967864
BS  0.08925278 0.01368505
````


- The last one is the Monte Carlo approximation of dynamic prediction, which is the Monte Carlo version of the *DP* function with a new argument as follows:

- mi the number of multiple imputation for Monte-Carlo approximation; default is 10.


Using this function facilitates the computation of credible intervals for each prediction.

```
MCDP <- MCDP(VS, Step2,
  Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
  n.thin = 1, cause_main = 1, mi = 10,
  DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = MCDP$DP$est, cause = 1
)
```
with the following outputs:
```
$DP
     id         est       lower       upper
1     2 0.716386507 0.675167218 0.755475954
2     5 0.114635792 0.107778261 0.124362366
3     7 0.149360907 0.144127688 0.160760439
4    10 0.061768634 0.058670334 0.064408627
5    11 0.031963961 0.030191980 0.034225293
.
.
.
247 495 0.056800414 0.051971578 0.061358109
248 497 0.119515472 0.114530735 0.125109349
249 498 0.045570271 0.042719707 0.048772580
250 500 0.058485596 0.054425111 0.063311557

$Cri
           est         sd
AUC 0.74209770 0.05029602
BS  0.08990641 0.01377804
```
