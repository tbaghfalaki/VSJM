### Variable Selection and Dynamic Prediction in Joint Models of Multiple Longitudinal Measures and Competing Risks
Employing a two-stage strategy for variable selection in joint modeling of multiple longitudinal measures and competing risks outcomes. The initial stage involves the estimation of K one-marker joint models for both the competing risks and each marker. Subsequently, in the second stage, time-varying covariates are integrated into a proportional hazard model to estimate parameters. Both phases follow the Bayesian paradigm, with the latter stage using a spike-and-slab prior to improve variable selection within the competing risks model. Additionally, this approach facilitates the computation of dynamic predictions based on the estimated parameter values.


### Installation
To acquire the latest development version of VSJM, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/VSJM")
```
This will seamlessly fetch and install the most up-to-date version of VSJM for your use.

### Reference 
Baghfalaki, T., Hashemi, R. & Jacqmin-Gadda, H. (2024). A Two-stage Approach for Variable Selection in Joint Modeling of Multiple Longitudinal Markers and Competing Risk Outcomes. 
