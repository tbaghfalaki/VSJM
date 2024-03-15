#' Variable Selection Joint Modeling
#'
#' @description
#' Employ a two-stage approach to efficiently select variables for joint modeling of multiple longitudinal markers and competing risks outcomes within a Bayesian framework.
#'
#' @details
#' A two-stage variable selection approach within a Bayesian framework is designed to handle multiple longitudinal measurements and competing risks outcomes. In the first stage, we address the estimation of a one-marker joint model for the event along with each longitudinal marker. Leveraging these estimates, we derive predictions for the expected values or slopes of individual marker trajectories. In the second stage, we utilize a proportional hazard model that incorporates the expected current values and/or slopes of all markers as time-dependent covariates. Here, we adopt a spike-and-slab prior to effectively select crucial longitudinal markers and covariates, enhancing the model's predictive power and interpretability. All implementations have been done with JAGS.
#'
#' @param formFixed a list of formulas for fixed part of longitudinal model
#' @param formRandom a list of formulas for random part of longitudinal model
#' @param formGroup a list of formulas specifying the cluster variable for Y (e.g. = list (~ subject, ~ subject,...))
#' @param formSurv formula for survival model
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param nmark the number of longitudinal markers
#' @param K1 Number of nodes and weights for calculating Gaussian quadrature in the first stage.
#' @param K2 Number of nodes and weights for calculating Gaussian quadrature in the second stage.
#' @param model a list of the models for the longitudinal part which includes "linear" or "quadratic".
#' @param Obstime the observed time in longitudinal data
#' @param ncl the number of nodes to be forked for parallel computing
#' @param Method the method for variable selection including "CS" for continues spike and "DS" for Dirac spike.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#' @importFrom stats quantile rnorm model.frame model.matrix
#'
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


VS <- function(formFixed, formRandom, formGroup, formSurv, nmark, K1 = K1, K2 = K2,
               model = model, n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
               n.thin = max(1, floor((n.iter - n.burnin) / 1000)), Obstime = "obstime", Method = "DS", ncl = ncl,
               DIC = TRUE, quiet = FALSE, dataLong, dataSurv) {
  j <- 1:nmark
  boot_fx <- function(j) {
    A1 <- UJM(
      formFixed = formFixed[[j]], formRandom = formRandom[[j]],
      formGroup = formGroup[[j]], formSurv = formSurv, dataLong = dataLong,
      dataSurv = dataSurv, K = K1, model = model[[j]], Obstime = Obstime,
      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = DIC, quiet = quiet
    )


    list(sim = A1$MCMC, PMean = A1$PMean, Long = A1$Estimation)
  }

  results <- parallel::mclapply(j, boot_fx, mc.cores = ncl)

  print(results[[1]]$PMean$linearpred)
  print(results[[2]]$PMean$linearpred)
  print(results[[3]]$PMean$linearpred)
  print(length(results[[1]]$PMean$linearpred))
  print(length(results[[2]]$PMean$linearpred))
  print(length(results[[3]]$PMean$linearpred))
  #############
  gamma <- sigma <- c()
  K <- K2
  X <- Z <- Xv <- Zv <- Nb <- list()
  indB <- indtime <- list()
  for (j in 1:nmark) {
    if (model[[j]] == "intercept") {
      data_long <- dataLong[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      data_long <- data_long[is.na(y) == FALSE, ]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      mfX <- stats::model.frame(formFixed[[j]], data = data_long) # , na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX) # , na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long) # , na.action = NULL)
      id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
      n2 <- length(unique(id))

      M <- table(id)
      id2 <- rep(1:length(M), M)

      Obstime <- Obstime
      Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications

      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime] # index of time
    }
    ####
    if (model[[j]] == "linear") {
      data_long <- dataLong[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      data_long <- data_long[is.na(y) == FALSE, ]
      y <- data_long[all.vars(formFixed[[j]])][, 1]

      mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
      Z[[j]] <- stats::model.matrix(formRandom[[j]], mfU, na.action = NULL)
      id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])

      n2 <- length(unique(id))
      Nb[[j]] <- dim(Z[[j]])[2]
      Obstime <- Obstime
      Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1]

      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime] # index of time
    }
    #############
    if (model[[j]] == "quadratic") {
      data_long <- dataLong[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      data_long <- data_long[is.na(y) == FALSE, ]
      y <- data_long[all.vars(formFixed[[j]])][, 1]


      mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
      Z[[j]] <- stats::model.matrix(formRandom[[j]], mfU, na.action = NULL)


      colnamesmfX <- colnames(X[[j]])
      colnamesmfU <- colnames(Z[[j]])
      Obstime <- Obstime

      Xtime2 <- mfX[, Obstime]^2

      X[[j]] <- cbind(X[[j]], Xtime2)
      colnames(X[[j]]) <- c(colnamesmfX, "obstime2")
      Z[[j]] <- cbind(Z[[j]], Xtime2)
      colnames(Z[[j]]) <- c(colnamesmfU, "obstime2")
      Obstime2 <- "obstime2"

      id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])


      n2 <- length(unique(id))
      Nb[[j]] <- dim(Z[[j]])[2]

      Obstime2n <- c(Obstime, Obstime2)
      Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime2n)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications


      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime2n] # index of time
    }
  }
  ##########
  n <- length(id)

  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  CR <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  suppressWarnings({
    mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
  })
  XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]
  #########
  n2 <- dim(dataSurv)[1]
  gamma <- sigma <- c()
  mu1 <- matrix(0, n2, nmark)
  betaL <- b <- list()
  for (j in 1:nmark) {
    betaL[[j]] <- results[[j]]$PMean$beta
    gamma <- append(gamma, results[[j]]$PMean$alpha)
    sigma <- append(sigma, results[[j]]$PMean$sigma)
    mu1[, j] <- results[[j]]$PMean$linearpred
    b[[j]] <- results[[j]]$PMean$b
  }
  indtime <- nindtime <- list()
  for (j in 1:nmark) {
    indB <- 1:dim(X[[j]])[2]
    if (model[[j]] == "intercept") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
    }
    if (model[[j]] == "linear") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
    }
    if (model[[j]] == "quadratic") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% c(Obstime, Obstime2)]
    }
    nindtime[[j]] <- c(1:dim(X[[j]])[2])[-indtime[[j]]]
  }
  Lp1 <- Lp2 <- Lp3 <- matrix(0, n2, nmark)
  for (i in 1:n2) {
    for (j in 1:nmark) {
      if (is.matrix(Xv[[j]]) == TRUE) {
        if (model[[j]] != "intercept") {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i, 1]
        } else {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i]
        }
      } else {
        if (model[[j]] != "intercept") {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i, 1]
        } else {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i]
        }
      }
      if (model[[j]] != "intercept") {
        Lp2[i, j] <- betaL[[j]][indtime[[j]][1]] + b[[j]][i, 2]
      } else {
        Lp2[i, j] <- betaL[[j]][indtime[[j]][1]]
      }
      Lp3[i, j] <- 0


      if (model[[j]] == "quadratic") (Lp3[i, j] <- betaL[[j]][indtime[[j]][2]] + b[[j]][i, 3])
    }
  }

  ########  Gauss-Legendre quadrature (15 points)  ########

  glq <- statmod::gauss.quad(K, kind = "legendre")
  xk <- glq$nodes # Nodes
  wk <- glq$weights # Weights
  K <- length(xk) # K-points
  ####################### Univariate #######################
  peice <- quantile(Time, seq(.2, 0.8, length = 4))
  delta <- nnet::class.ind(arules::discretize(Time, method = "fixed", c(0, peice, max(Time))))

  model_S2CS <- "model{

  # Scaling Gauss-Kronrod/Legendre quadrature

  for(k in 1:n){
# Scaling Gauss-Kronrod/Legendre quadrature
  for(j in 1:K){
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
  }
    #
    for(l in 1:C){
      Alpha0[k,l]<-  inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+inprod(alpha[l,],LP1[k,])
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
    phi[k]<-100000+sum(phi1[k,])
    zeros[k]~dpois(phi[k])
  }
  #Prior distributions
  for(k in 1:NbetasS){
    for(l in 1:C){
      betaS[l,k]<-xi1[l,k]*B11[l,k]+ (1- xi1[l,k])*B21[l,k]
      B11[l,k]~dnorm(0,tau1[l,k])
      B21[l,k]~dnorm(0,1000)
      xi1[l,k] ~ dbern(kappa1[l,k])
      kappa1[l,k]~dbeta(.1,.1)
      lBFDR1[l,k]<- 1-xi1[l,k]
      tau1[l,k]~dgamma(0.01,0.01)

    }}

  for(j in 1:J){
    for(l in 1:C){
      h[j,l]~dgamma(0.1,0.1)
    }}


    for(j in 1:nmark){
    for(l in 1:C){
     alpha[l,j]<-xi2[l,j]*B12[l,j]+ (1- xi2[l,j])*B22[l,j]
      B12[l,j]~dnorm(0,tau2[l,j])
      B22[l,j]~dnorm(0,1000)
      xi2[l,j] ~ dbern(kappa2[l,j])
      kappa2[l,j]~dbeta(.1,.1)
      lBFDR2[l,j]<- 1-xi2[l,j]
      tau2[l,j]~dgamma(0.01,0.01)
  }}

}"
  ######
  model_S2DS <- "model{


  for(k in 1:n){
    # Scaling Gauss-Kronrod/Legendre quadrature
  for(j in 1:K){
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
  }
    #
    for(l in 1:C){
      Alpha0[k,l]<-  inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+inprod(alpha[l,],LP1[k,])
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
    phi[k]<-100000+sum(phi1[k,])
    zeros[k]~dpois(phi[k])
  }
  #Prior distributions
  for(k in 1:NbetasS){
    for(l in 1:C){
      betaS[l,k]<-ind1[l,k]*B21[l,k]
      ind1[l,k] <- equals(xi1[l,k], 0)
      xi1[l,k] ~ dbern(kappa1[l,k])
      kappa1[l,k]~dbeta(.01,.01)
      B21[l,k]~dnorm(0,tau1[l,k])
      lBFDR1[l,k]<- 1-ind1[l,k]
      tau1[l,k]~dgamma(0.01,0.01)

    }}

  for(j in 1:J){
    for(l in 1:C){
      h[j,l]~dgamma(0.1,0.1)
    }}


  for(j in 1:nmark){
    for(l in 1:C){
      alpha[l,j]<-ind2[l,j]*B22[l,j]
      ind2[l,j] <- equals(xi2[l,j], 0)
      xi2[l,j] ~ dbern(kappa2[l,j])
      kappa2[l,j]~dbeta(.01,.01)
      B22[l,j]~dnorm(0,tau2[l,j])
      lBFDR2[l,j]<- 1-ind2[l,j]
      tau2[l,j]~dgamma(0.01,0.01)
    }}

}"



  C <- length(unique(CR)) - 1
  NbetasS <- dim(XS)[2]
  d.jags <- list(
    n = n2, mu1 = mu1, zeros = rep(0, n2), NbetasS = dim(XS)[2], C = C, nmark = nmark,
    LP1 = Lp1, LP2 = Lp2, LP3 = Lp3, Time = Time, CR = CR, XS = XS,
    s = peice, J = dim(delta)[2], xk = xk, wk = wk, K = K
  )

  i.jagsCS <- function() {
    list(B12 = matrix(0, C, nmark), B22 = matrix(0, C, nmark), B11 = matrix(0, C, NbetasS), B21 = matrix(0, C, NbetasS))
  }

  i.jagsDS <- function() {
    list(B22 = matrix(0, C, nmark), B21 = matrix(0, C, NbetasS))
  }
  parameters <- c("alpha", "betaS", "h", "lBFDR1", "lBFDR2")
  if (Method == "CS") {
    model.file <- textConnection(model_S2CS)
    i.jags <- i.jagsCS
  } else {
    model.file <- textConnection(model_S2DS)
    i.jags <- i.jagsDS
  }


  step2_unijm <- jagsUI::jags(
    data = d.jags,
    inits=i.jags,
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


  LBFDR1 <- step2_unijm$mean$lBFDR1
  BF1 <- (1 - LBFDR1) / LBFDR1
  rownames(LBFDR1) <- rownames(BF1) <- paste0("Cause ", 1:C)
  colnames(LBFDR1) <- colnames(BF1) <- colnames(XS)


  LBFDR2 <- step2_unijm$mean$lBFDR2
  BF2 <- (1 - LBFDR2) / LBFDR2
  rownames(LBFDR2) <- rownames(BF2) <- paste0("Cause ", 1:C)
  colnames(LBFDR2) <- colnames(BF2) <- paste0("Marker", 1:nmark)

  Longitudinal <- list()
  for (j in 1:nmark) {
    Longitudinal[[j]] <- results[[j]]$Long
  }

  list(
    formFixed = formFixed, formRandom = formRandom, formGroup = formGroup, formSurv = formSurv,
    model = model, Obstime = Obstime,
    mu1 = mu1,
    C = C, nmark = nmark, XS = XS, Lp1 = Lp1, Lp2 = Lp2, Lp3 = Lp3,
    sim_step1 = results, Longitudinal = Longitudinal, TDsurvival = step2_unijm, peice = peice,
    LBFDRX = LBFDR1, BFX = BF1, LBFDRY = LBFDR2, BFY = BF2
  )
}
