#'  Variable Selection Joint Modeling 2
#'
#' @description
#' Running a time-dependent proportional hazard model on the results of VS
#'
#'
#' @details
#' After variable selection by VS, the second stage is re-estimated by replacing CS or Ds with non-informative normal distributions as priors for the association parameters of the selected markers and the regression coefficients of the selected covariates.
#'
#' @param object an object inheriting from class VSJM
#' @param Method the method for variable selection including "LBFDR" for LBFDR and "BF" for Bayes factor.
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#'
#'
#' @importFrom stats quantile rnorm model.frame model.matrix
#'
#' @return
#' - MCMC chains for the unknown parameters
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 25% list of posterior median for each parameter
#' - 50% list of posterior median for each parameter
#' - 75% list of posterior median for each parameter
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#'
#' @example inst/exampleVS.R
#'
#' @md
#' @export
VS2 <- function(object, Method = "LBFDR", n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                DIC = TRUE, quiet = FALSE, dataLong, dataSurv) {
  LBFDR2 <- object$LBFDRY
  BF2 <- object$BFY

  I_LBFDR_alpha <- 1 * (LBFDR2 < 0.05)
  I_BF_alpha <- 1 * (BF2 > 1)
  if (Method == "LBFDR") {
    I_alpha <- I_LBFDR_alpha
  } else {
    I_alpha <- I_BF_alpha
  }

  LBFDR1 <- object$LBFDRX
  BF1 <- object$BFX

  I_LBFDR_betaS <- 1 * (LBFDR1 < 0.05)
  I_BF_betaS <- 1 * (BF1 > 1)
  if (Method == "LBFDR") {
    I_betaS <- I_LBFDR_betaS
  } else {
    I_betaS <- I_BF_betaS
  }
  #######
  formFixed <- object$formFixed
  formRandom <- object$formRandom
  formGroup <- object$formGroup
  formSurv <- object$formSurv
  model <- object$model
  Obstime <- object$Obstime

  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  CR <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  suppressWarnings({
    mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
  })
  XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]

  C <- object$C
  nmark <- object$nmark


  n <- dim(dataSurv)[1]
  mu1 <- object$mu1
  #######

  Lp1 <- object$Lp1
  Lp2 <- object$Lp2
  Lp3 <- object$Lp3

  LP1 <- LP2 <- LP3 <- array(0, c(dim(Lp1), C))
  for (l in 1:C) {
    LP1[, , l] <- Lp1
    LP2[, , l] <- Lp2
    LP3[, , l] <- Lp3
  }

  ########## Step_2 without variable selection  ##############
  model_S2 <- "model{
  for(k in 1:n){
    # Scaling Gauss-Kronrod/Legendre quadrature
  for(j in 1:K){
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
    }
    for(l in 1:C){
      Alpha0[k,l]<- inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+inprod(alpha[l,],LP1[k,])
      Alpha1[k,l]<- inprod(alpha[l,],LP2[k,])
      Alpha2[k,l]<- inprod(alpha[l,],LP3[k,])

      haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                    (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                    (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                    (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                    (h[5,l]*step(Time[k]-s[4])))*
        exp(inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+inprod(alpha[l,],mu1[k,]))


      for(j in 1:K){
        chaz[k,j,l]<-  ((h[1,l]*step(s[1]-xk11[k,j]))+
                          (h[2,l]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                          (h[3,l]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                          (h[4,l]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                          (h[5,l]*step(xk11[k,j]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*xk11[k,j]+Alpha2[k,l]*pow(xk11[k,j],2))
      }
      logSurv[k,l]<- -inprod(wk11[k,],chaz[k,,l])
      phi1[k,l]<--equals(CR[k],l)*log(haz[k,l])-logSurv[k,l]
    }
    #Definition of the survival log-likelihood using zeros trick
    phi[k]<-KK+sum(phi1[k,])
    zeros[k]~dpois(phi[k])
  }
  #Prior distributions


for(k in 1:NbetasS){
    for(l in 1:C){
      betaS[l,k]<-I_betaS[l,k]*B11[l,k]
      B11[l,k]~dnorm(0,0.0001)
    }}



  for(j in 1:J){
    for(l in 1:C){
      h[j,l]~dgamma(0.1,0.1)
    }}


for(k in 1:NbetasS){
    for(l in 1:C){
      alpha[l,k]<-I_alpha[l,k]*a11[l,k]
      a11[l,k]~dnorm(0,0.0001)
    }}


}"

  glq <- statmod::gauss.quad(15, kind = "legendre")
  xk <- glq$nodes # Nodes
  wk <- glq$weights # Weights
  K <- length(xk) # K-points

  peice <- object$peice
  C <- length(unique(CR)) - 1
  NbetasS <- dim(XS)[2]
  d.jags <- list(
    n = n, mu1 = mu1, zeros = rep(0, n), NbetasS = dim(XS)[2], C = C, nmark = nmark, I_alpha = I_alpha, I_betaS = I_betaS,
    LP1 = Lp1, LP2 = Lp2, LP3 = Lp3, Time = Time, CR = CR, XS = XS, KK = 100000,
    s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
  )

  i.jags <- function() {
    list(alpha1 = matrix(1, C, nmark) * I_alpha, betaS1 = matrix(NA, C, NbetasS) * I_betaS)
  }


  parameters <- c("alpha", "betaS", "h")
  model.file <- textConnection(model_S2)

  sim1 <- jagsUI::jags(
    data = d.jags,
    inits = i.jags,
    parameters.to.save = parameters,
    model.file = model.file,
    n.chains = n.chains,
    parallel = FALSE,
    n.adapt = FALSE,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    DIC = TRUE
  )

  ############
  MCMC <- list(
    beta = sim1$sims.list$alpha, gamma = sim1$sims.list$betaS,
    h = sim1$sims.list$h
  )

  if (n.chains > 1) {
    colnames(sim1$mean$alpha) <-
      colnames(sim1$sd$alpha) <-
      colnames(sim1$q2.5$alpha) <-
      colnames(sim1$q97.5$alpha) <-
      colnames(sim1$Rhat$alpha) <- paste0("Marker", 1:nmark)
    rownames(sim1$mean$alpha) <-
      rownames(sim1$sd$alpha) <-
      rownames(sim1$q2.5$alpha) <-
      rownames(sim1$q97.5$alpha) <-
      rownames(sim1$Rhat$alpha) <- paste0("Cause", 1:C)



    colnames(sim1$mean$betaS) <-
      colnames(sim1$sd$betaS) <-
      colnames(sim1$q2.5$betaS) <-
      colnames(sim1$q97.5$betaS) <-
      colnames(sim1$Rhat$betaS) <- colnames(XS)

    rownames(sim1$mean$betaS) <-
      rownames(sim1$sd$betaS) <-
      rownames(sim1$q2.5$betaS) <-
      rownames(sim1$q97.5$betaS) <-
      rownames(sim1$Rhat$betaS) <- paste0("Cause", 1:C)



    rownames(sim1$mean$h) <-
      rownames(sim1$sd$h) <-
      rownames(sim1$q2.5$h) <-
      rownames(sim1$q97.5$h) <-
      rownames(sim1$Rhat$h) <- paste0("lambda", 1:(length(peice) + 1))


    colnames(sim1$mean$h) <-
      colnames(sim1$sd$h) <-
      colnames(sim1$q2.5$h) <-
      colnames(sim1$q97.5$h) <-
      colnames(sim1$Rhat$h) <- paste0("Cause", 1:C)


    LM <- list(
      gamma = list(
        Est = sim1$mean$betaS,
        SD = sim1$sd$betaS,
        L_CI = sim1$q2.5$betaS,
        U_CI = sim1$q97.5$betaS,
        R_hat = sim1$Rhat$betaS
      ),
      alpha = list(
        Est = sim1$mean$alpha,
        SD = sim1$sd$alpha,
        L_CI = sim1$q2.5$alpha,
        U_CI = sim1$q97.5$alpha,
        R_hat = sim1$Rhat$alpha
      ),
      lambda = list(
        Est = sim1$mean$h,
        SD = sim1$sd$h,
        L_CI = sim1$q2.5$h,
        U_CI = sim1$q97.5$h,
        R_hat = sim1$Rhat$h
      )
    )



    results <- list(Survival_model = LM)
  } else {
    colnames(sim1$mean$alpha) <-
      colnames(sim1$sd$alpha) <-
      colnames(sim1$q2.5$alpha) <-
      colnames(sim1$q97.5$alpha) <- paste0("Marker", 1:nmark)
    rownames(sim1$mean$alpha) <-
      rownames(sim1$sd$alpha) <-
      rownames(sim1$q2.5$alpha) <-
      rownames(sim1$q97.5$alpha) <- paste0("Cause", 1:C)



    colnames(sim1$mean$betaS) <-
      colnames(sim1$sd$betaS) <-
      colnames(sim1$q2.5$betaS) <-
      colnames(sim1$q97.5$betaS) <- colnames(XS)

    rownames(sim1$mean$betaS) <-
      rownames(sim1$sd$betaS) <-
      rownames(sim1$q2.5$betaS) <-
      rownames(sim1$q97.5$betaS) <- paste0("Cause", 1:C)



    rownames(sim1$mean$h) <-
      rownames(sim1$sd$h) <-
      rownames(sim1$q2.5$h) <-
      rownames(sim1$q97.5$h) <- paste0("lambda", 1:(length(peice) + 1))


    colnames(sim1$mean$h) <-
      colnames(sim1$sd$h) <-
      colnames(sim1$q2.5$h) <-
      colnames(sim1$q97.5$h) <- paste0("Cause", 1:C)


    LM <- list(
      gamma = list(
        Est = sim1$mean$betaS,
        SD = sim1$sd$betaS,
        L_CI = sim1$q2.5$betaS,
        U_CI = sim1$q97.5$betaS
      ),
      alpha = list(
        Est = sim1$mean$alpha,
        SD = sim1$sd$alpha,
        L_CI = sim1$q2.5$alpha,
        U_CI = sim1$q97.5$alpha
      ),
      lambda = list(
        Est = sim1$mean$h,
        SD = sim1$sd$h,
        L_CI = sim1$q2.5$h,
        U_CI = sim1$q97.5$h
      )
    )



    results <- list(Survival_model = LM)
  }

  KK <- 100000
  DIC <- sim1$DIC - 2 * KK * n


  list(
    Results = sim1,
    MCMC = MCMC,
    Estimation = results, DIC = DIC
  )
}

