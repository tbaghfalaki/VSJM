#' One-marker Joint Modeling
#'
#' @description
#' Fits Bayesian linear or quadratic one-marker joint modeling of longitudinal measurements and competing risks using MCMC
#'
#' @details
#' Use the 'JAGS' software to estimate Bayesian linear or quadratic one-marker joint modeling of longitudinal measurements and competing risks using MCMC
#'
#' @param formFixed formula for fixed part of longitudinal model
#' @param formRandom formula for random part of longitudinal model
#' @param formGroup formula specifying the cluster variable for Y (e.g. = ~ subject)
#' @param formSurv formula for survival model
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param K Number of nodes and weights for calculating Gaussian quadrature.
#' @param model the model for the longitudinal part which includes "intercept", "linear" or "quadratic".
#' @param Obstime the observed time in longitudinal data
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
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#' @example inst/exampleUJM.R
#'
#' @md
#' @export
UJM <- function(formFixed, formRandom, formGroup, formSurv, dataLong, dataSurv, K = K, model = "linear", Obstime = "Obstime",
                n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                DIC = TRUE, quiet = FALSE) {
  KK <- 1000000
  if (model == "intercept") {
    data_long <- dataLong[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
    y <- data_long[all.vars(formFixed)][, 1]
    data_long <- data_long[is.na(y) == FALSE, ]
    y <- data_long[all.vars(formFixed)][, 1]
    mfX <- stats::model.frame(formFixed, data = data_long) # , na.action = NULL)
    X <- stats::model.matrix(formFixed, mfX) # , na.action = NULL)
    mfU <- stats::model.frame(formRandom, data = data_long) # , na.action = NULL)
    id <- as.integer(data_long[all.vars(formGroup)][, 1])
    n2 <- length(unique(id))

    M <- table(id)
    id2 <- rep(1:length(M), M)

    Obstime <- Obstime
    Xvtime <- cbind(id, X[, colnames(X) %in% setdiff(colnames(X), Obstime)])
    Xv <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications

    indB <- 1:dim(X)[2]
    indtime <- indB[colnames(X) %in% Obstime] # index of time

    ## Number of patients and number of longitudinal observations per patient
    n <- length(id)
    tmp <- dataSurv[all.vars(formSurv)]
    Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
    CR <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
    nTime <- length(Time) # number of subject having Time
    # design matrice
    suppressWarnings({
      mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
    })
    # XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]
    XS <- stats::model.matrix(formSurv, mfZ)[, -1]

    ########  Gauss-Legendre quadrature (15 points)  ########

    glq <- statmod::gauss.quad(K, kind = "legendre")
    xk <- glq$nodes # Nodes
    wk <- glq$weights # Weights
    K <- length(xk) # K-points
    ################

    peice <- stats::quantile(Time, seq(.2, 0.8, length = 4))

    ########### univariate_jm
    model_fileI1 <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+b[id[i]]
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }}
  for(k in 1:n2){
    linearpred[k]<-betaL[1]+b[k]+(betaL[indtime])*Time[k]
    for(l in 1:C){
    Alpha0[k,l]<- inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+
                 gamma1[l]*(betaL[1]+b[k])


    Alpha1[k,l]<- gamma1[l]*(betaL[indtime])
    haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5,l]*step(Time[k]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j,l]<-  ((h[1,l]*step(s[1]-xk11[k,j]))+
                      (h[2,l]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3,l]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4,l]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5,l]*step(xk11[k,j]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*xk11[k,j])

    }


    logSurv[k,l]<- -inprod(wk11[k,],chaz[k,,l])
    phi1[k,l]<--equals(CR[k],l)*log(haz[k,l])-logSurv[k,l]

    }
    #Definition of the survival log-likelihood using zeros trick

    phi[k]<-KK+sum(phi1[k,])

    zeros[k]~dpois(phi[k])
    #Random effects
    b[k]~dnorm(0,Omega)

  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.001)
  }

  for(k in 1:NbetasS){
    for(l in 1:C){
    betaS[l,k]~dnorm(0,0.001)
    }}


  for(j in 1:J){
    for(l in 1:C){

    h[j,l]~dgamma(0.1,0.1)
  }}

  for(l in 1:C){
  gamma1[l]~dnorm(0,0.001)
  }
  Omega~dgamma(0.01,0.01)
  #Derive dquantity
  Sigma<-1/Omega

  tau1~dgamma(0.01,0.01)
  sigma1<-1/tau1
}"

    model_fileI <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+b[id[i]]
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }}
  for(k in 1:n2){
    linearpred[k]<-inprod(betaL[nindtime],Xv[k,])+b[k]+(betaL[indtime])*Time[k]
    for(l in 1:C){
    Alpha0[k,l]<- inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+
                 gamma1[l]*(inprod(betaL[nindtime],Xv[k,])+b[k])


    Alpha1[k,l]<- gamma1[l]*(betaL[indtime])
    haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5,l]*step(Time[k]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j,l]<-  ((h[1,l]*step(s[1]-xk11[k,j]))+
                      (h[2,l]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3,l]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4,l]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5,l]*step(xk11[k,j]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*xk11[k,j])

    }


    logSurv[k,l]<- -inprod(wk11[k,],chaz[k,,l])
    phi1[k,l]<--equals(CR[k],l)*log(haz[k,l])-logSurv[k,l]

    }
    #Definition of the survival log-likelihood using zeros trick

    phi[k]<-KK+sum(phi1[k,])

    zeros[k]~dpois(phi[k])
    #Random effects
    b[k]~dnorm(0,Omega)

  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.001)
  }

  for(k in 1:NbetasS){
    for(l in 1:C){
    betaS[l,k]~dnorm(0,0.001)
    }}


  for(j in 1:J){
    for(l in 1:C){

    h[j,l]~dgamma(0.1,0.1)
  }}

  for(l in 1:C){
  gamma1[l]~dnorm(0,0.001)
  }
  Omega~dgamma(0.01,0.01)
  #Derive dquantity
  Sigma<-1/Omega

  tau1~dgamma(0.01,0.01)
  sigma1<-1/tau1
}"



    C <- length(unique(CR)) - 1
    i.jags <- function() {
      list(
        gamma1 = stats::rnorm(C),
        betaL = stats::rnorm(dim(X)[2]), tau1 = 1, Omega = 1, b = rep(0, n2)
      )
    }

    parameters <- c(
      "betaL", "betaS", "linearpred", "b",
      "gamma1", "sigma1", "Sigma", "h"
    )
    ###############################
    d.jags1 <- list(
      n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
      X = X, id = id2, indtime = indtime,
      CR = CR, zeros = rep(0, n2),
      NbetasL = dim(X)[2], s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, KK = KK
    )

    d.jags <- list(
      n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
      X = X, id = id2, Xv = Xv, indtime = indtime, nindtime = c(1:dim(X)[2])[-indtime],
      CR = CR, zeros = rep(0, n2),
      NbetasL = dim(X)[2], s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, KK = KK
    )
    if (is.matrix(Xv) == FALSE) {
      model_fileL_last <- textConnection(model_fileI1)

      d.jags <- d.jags1
    } else {
      model_fileL_last <- textConnection(model_fileI)
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      inits = i.jags,
      parameters.to.save = parameters,
      model.file = model_fileL_last,
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
      beta = sim1$sims.list$betaL, gamma = sim1$sims.list$betaS,
      alpha = sim1$sims.list$gamma1,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma1,
      Sigma = sim1$sims.list$Sigma,
      linearpred = sim1$sims.list$linearpred,
      b = sim1$sims.list$b
    )
    PMean <- list(
      beta = sim1$mean$betaL, gamma = sim1$mean$betaS,
      alpha = sim1$mean$gamma1,
      h = sim1$mean$h,
      sigma = sim1$mean$sigma1,
      Sigma = sim1$mean$Sigma,
      linearpred = sim1$mean$linearpred,
      b = sim1$mean$b
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL) <-
        names(sim1$sd$betaL) <-
        names(sim1$q2.5$betaL) <-
        names(sim1$q97.5$betaL) <-
        names(sim1$Rhat$betaL) <- colnames(X)

      names(sim1$mean$sigma1) <-
        names(sim1$sd$sigma1) <-
        names(sim1$q2.5$sigma1) <-
        names(sim1$q97.5$sigma1) <-
        names(sim1$Rhat$sigma1) <- "sigma2_e"


      names(sim1$mean$Sigma) <-
        names(sim1$sd$Sigma) <-
        names(sim1$q2.5$Sigma) <-
        names(sim1$q97.5$Sigma) <-
        names(sim1$Rhat$Sigma) <- "sigma2_b"





      LM <- rbind(
        cbind(sim1$mean$betaL, sim1$sd$betaL, sim1$q2.5$betaL, sim1$q97.5$betaL, sim1$Rhat$betaL),
        cbind(sim1$mean$sigma1, sim1$sd$sigma1, sim1$q2.5$sigma1, sim1$q97.5$sigma1, sim1$Rhat$sigma1),
        cbind(sim1$mean$Sigma, sim1$sd$Sigma, sim1$q2.5$Sigma, sim1$q97.5$Sigma, sim1$Rhat$Sigma)
      )


      colnames(LM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Longitudinal_model = LM)
    } else {
      names(sim1$mean$betaL) <-
        names(sim1$sd$betaL) <-
        names(sim1$q2.5$betaL) <-
        names(sim1$q97.5$betaL) <- colnames(X)

      names(sim1$mean$sigma1) <-
        names(sim1$sd$sigma1) <-
        names(sim1$q2.5$sigma1) <-
        names(sim1$q97.5$sigma1) <- "sigma2_e"


      names(sim1$mean$Sigma) <-
        names(sim1$sd$Sigma) <-
        names(sim1$q2.5$Sigma) <-
        names(sim1$q97.5$Sigma) <- "sigma2_b"





      LM <- rbind(
        cbind(sim1$mean$betaL, sim1$sd$betaL, sim1$q2.5$betaL, sim1$q97.5$betaL),
        cbind(sim1$mean$sigma1, sim1$sd$sigma1, sim1$q2.5$sigma1, sim1$q97.5$sigma1),
        cbind(sim1$mean$Sigma, sim1$sd$Sigma, sim1$q2.5$Sigma, sim1$q97.5$Sigma)
      )


      colnames(LM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Longitudinal_model = LM)
    }
  }
  #####


  if (model == "linear") {
    data_long <- dataLong[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]

    y <- data_long[all.vars(formFixed)][, 1]
    data_long <- data_long[is.na(y) == FALSE, ]
    y <- data_long[all.vars(formFixed)][, 1]
    mfX <- stats::model.frame(formFixed, data = data_long) # , na.action = NULL)
    X <- stats::model.matrix(formFixed, mfX) # , na.action = NULL)
    mfU <- stats::model.frame(formRandom, data = data_long) # , na.action = NULL)
    Z <- stats::model.matrix(formRandom, mfU) # , na.action = NULL)
    id <- as.integer(data_long[all.vars(formGroup)][, 1])
    n2 <- length(unique(id))

    M <- table(id)
    id2 <- rep(1:length(M), M)


    Nb <- dim(Z)[2]
    Obstime <- Obstime
    Xvtime <- cbind(id, X[, colnames(X) %in% setdiff(colnames(X), Obstime)])
    Xv <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications

    indB <- 1:dim(X)[2]
    indtime <- indB[colnames(X) %in% Obstime] # index of time

    ## Number of patients and number of longitudinal observations per patient
    n <- length(id)

    tmp <- dataSurv[all.vars(formSurv)]
    Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
    CR <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
    nTime <- length(Time) # number of subject having Time
    # design matrice
    suppressWarnings({
      mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
    })
    # XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]
    XS <- stats::model.matrix(formSurv, mfZ)[, -1]

    ########  Gauss-Legendre quadrature (15 points)  ########

    glq <- statmod::gauss.quad(K, kind = "legendre")
    xk <- glq$nodes # Nodes
    wk <- glq$weights # Weights
    K <- length(xk) # K-points
    ################

    peice <- stats::quantile(Time, seq(.2, 0.8, length = 4))
    delta <- nnet::class.ind(arules::discretize(Time, method = "fixed", c(0, peice, max(Time))))

    ########### univariate_jm
    model_fileL1 <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }}
  for(k in 1:n2){
    linearpred[k]<-betaL[1]+b[k,1]+(betaL[indtime]+b[k,2])*Time[k]
    for(l in 1:C){
    Alpha0[k,l]<- inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+
                 gamma1[l]*(betaL[1]+b[k,1])


    Alpha1[k,l]<- gamma1[l]*(betaL[indtime]+b[k,2])
    haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5,l]*step(Time[k]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j,l]<-  ((h[1,l]*step(s[1]-xk11[k,j]))+
                      (h[2,l]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3,l]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4,l]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5,l]*step(xk11[k,j]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*xk11[k,j])

    }


    logSurv[k,l]<- -inprod(wk11[k,],chaz[k,,l])
    phi1[k,l]<--equals(CR[k],l)*log(haz[k,l])-logSurv[k,l]

    }
    #Definition of the survival log-likelihood using zeros trick

    phi[k]<-KK+sum(phi1[k,])

    zeros[k]~dpois(phi[k])
    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.001)
  }

  for(k in 1:NbetasS){
    for(l in 1:C){
    betaS[l,k]~dnorm(0,0.001)
    }}


  for(j in 1:J){
    for(l in 1:C){

    h[j,l]~dgamma(0.1,0.1)
  }}

  for(l in 1:C){
  gamma1[l]~dnorm(0,0.001)
  }
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Derive dquantity
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])

  tau1~dgamma(0.01,0.01)
  sigma1<-1/tau1
}"

    model_fileL <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }}
  for(k in 1:n2){
    linearpred[k]<-inprod(betaL[nindtime],Xv[k,])+b[k,1]+(betaL[indtime]+b[k,2])*Time[k]
    for(l in 1:C){
    Alpha0[k,l]<- inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+
                 gamma1[l]*(inprod(betaL[nindtime],Xv[k,])+b[k,1])


    Alpha1[k,l]<- gamma1[l]*(betaL[indtime]+b[k,2])
    haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5,l]*step(Time[k]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j,l]<-  ((h[1,l]*step(s[1]-xk11[k,j]))+
                      (h[2,l]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3,l]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4,l]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5,l]*step(xk11[k,j]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*xk11[k,j])

    }


    logSurv[k,l]<- -inprod(wk11[k,],chaz[k,,l])
    phi1[k,l]<--equals(CR[k],l)*log(haz[k,l])-logSurv[k,l]

    }
    #Definition of the survival log-likelihood using zeros trick

    phi[k]<-KK+sum(phi1[k,])

    zeros[k]~dpois(phi[k])
    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.001)
  }

  for(k in 1:NbetasS){
    for(l in 1:C){
    betaS[l,k]~dnorm(0,0.001)
    }}


  for(j in 1:J){
    for(l in 1:C){

    h[j,l]~dgamma(0.1,0.1)
  }}

  for(l in 1:C){
  gamma1[l]~dnorm(0,0.001)
  }
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Derive dquantity
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])

  tau1~dgamma(0.01,0.01)
  sigma1<-1/tau1
}"



    C <- length(unique(CR)) - 1
    i.jags <- function() {
      list(
        gamma1 = stats::rnorm(C),
        betaL = stats::rnorm(dim(X)[2]), tau1 = 1, Omega = diag(stats::runif(Nb)), b = matrix(0, n2, Nb)
      )
    }

    parameters <- c(
      "betaL", "betaS", "linearpred", "b",
      "gamma1", "sigma1", "Sigma", "h"
    )
    ###############################
    d.jags1 <- list(
      n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
      X = X, Z = Z, id = id2, indtime = indtime,
      CR = CR, mub = rep(0, Nb), V = diag(1, Nb), Nb = Nb, zeros = rep(0, n2),
      NbetasL = dim(X)[2], s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, KK = KK
    )

    d.jags <- list(
      n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
      X = X, Z = Z, id = id2, Xv = Xv, indtime = indtime, nindtime = c(1:dim(X)[2])[-indtime],
      CR = CR, mub = rep(0, Nb), V = diag(1, Nb), Nb = Nb, zeros = rep(0, n2),
      NbetasL = dim(X)[2], s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, KK = KK
    )
    if (is.matrix(Xv) == FALSE) {
      model_fileL_last <- textConnection(model_fileL1)
      d.jags <- d.jags1
    } else {
      model_fileL_last <- textConnection(model_fileL)
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      inits = i.jags,
      parameters.to.save = parameters,
      model.file = model_fileL_last,
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
      beta = sim1$sims.list$betaL, gamma = sim1$sims.list$betaS,
      alpha = sim1$sims.list$gamma1,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma1,
      Sigma = sim1$sims.list$Sigma,
      linearpred = sim1$sims.list$linearpred,
      b = sim1$sims.list$b
    )

    PMean <- list(
      beta = sim1$mean$betaL, gamma = sim1$mean$betaS,
      alpha = sim1$mean$gamma1,
      h = sim1$mean$h,
      sigma = sim1$mean$sigma1,
      Sigma = sim1$mean$Sigma,
      linearpred = sim1$mean$linearpred,
      b = sim1$mean$b
    )
    if (n.chains > 1) {
      names(sim1$mean$betaL) <-
        names(sim1$sd$betaL) <-
        names(sim1$q2.5$betaL) <-
        names(sim1$q97.5$betaL) <-
        names(sim1$Rhat$betaL) <- colnames(X)

      names(sim1$mean$sigma1) <-
        names(sim1$sd$sigma1) <-
        names(sim1$q2.5$sigma1) <-
        names(sim1$q97.5$sigma1) <-
        names(sim1$Rhat$sigma1) <- "sigma2"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c("Intercept", "Time")



      LM <- rbind(
        cbind(sim1$mean$betaL, sim1$sd$betaL, sim1$q2.5$betaL, sim1$q97.5$betaL, sim1$Rhat$betaL),
        cbind(sim1$mean$sigma1, sim1$sd$sigma1, sim1$q2.5$sigma1, sim1$q97.5$sigma1, sim1$Rhat$sigma1)
      )

      Sigma <- sim1$mean$Sigma

      colnames(LM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Longitudinal_model = LM, Sigma = Sigma)
    } else {
      names(sim1$mean$betaL) <-
        names(sim1$sd$betaL) <-
        names(sim1$q2.5$betaL) <-
        names(sim1$q97.5$betaL) <- colnames(X)

      names(sim1$mean$sigma1) <-
        names(sim1$sd$sigma1) <-
        names(sim1$q2.5$sigma1) <-
        names(sim1$q97.5$sigma1) <- "sigma2"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c("Intercept", "Time")



      LM <- rbind(
        cbind(sim1$mean$betaL, sim1$sd$betaL, sim1$q2.5$betaL, sim1$q97.5$betaL),
        cbind(sim1$mean$sigma1, sim1$sd$sigma1, sim1$q2.5$sigma1, sim1$q97.5$sigma1)
      )


      colnames(LM) <- c("Est", "SD", "L_CI", "U_CI")

      Sigma <- sim1$mean$Sigma


      results <- list(Longitudinal_model = LM, Sigma = Sigma)
    }
  }
  #####
  if (model == "quadratic") {
    data_long <- dataLong[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
    y <- data_long[all.vars(formFixed)][, 1]
    data_long <- data_long[is.na(y) == FALSE, ]
    y <- data_long[all.vars(formFixed)][, 1]
    mfX <- stats::model.frame(formFixed, data = data_long) # , na.action = NULL)
    X <- stats::model.matrix(formFixed, mfX) # , na.action = NULL)
    mfU <- stats::model.frame(formRandom, data = data_long) # , na.action = NULL)
    Z <- stats::model.matrix(formRandom, mfU) # , na.action = NULL)
    id <- as.integer(data_long[all.vars(formGroup)][, 1])

    M <- table(id)
    id2 <- rep(1:length(M), M)


    colnamesmfX <- colnames(X)
    colnamesmfU <- colnames(Z)
    Obstime <- Obstime

    Xtime2 <- mfX[, Obstime]^2

    X <- cbind(X, Xtime2)
    colnames(X) <- c(colnamesmfX, "obstime2")
    Z <- cbind(Z, Xtime2)
    colnames(Z) <- c(colnamesmfU, "obstime2")
    Obstime2 <- "obstime2"



    n2 <- length(unique(id))
    Nb <- dim(Z)[2]

    Obstime2n <- c(Obstime, Obstime2)
    Xvtime <- cbind(id, X[, colnames(X) %in% setdiff(colnames(X), Obstime2n)])
    Xv <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications


    indB <- 1:dim(X)[2]
    indtime <- indB[colnames(X) %in% Obstime2n] # index of time


    ## Number of patients and number of longitudinal observations per patient
    n <- length(id)

    tmp <- dataSurv[all.vars(formSurv)]
    Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
    CR <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
    nTime <- length(Time) # number of subject having Time
    # design matrice

    # design matrice
    suppressWarnings({
      mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
    })
    XS <- stats::model.matrix(formSurv, mfZ)[, -1]
    #    XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]

    ########  Gauss-Legendre quadrature (15 points)  ########
    glq <- statmod::gauss.quad(K, kind = "legendre")
    xk <- glq$nodes # Nodes
    wk <- glq$weights # Weights
    K <- length(xk) # K-points
    ################

    peice <- stats::quantile(Time, seq(.2, 0.8, length = 4))
    delta <- nnet::class.ind(arules::discretize(Time, method = "fixed", c(0, peice, max(Time))))

    ########### univariate_jm

    model_fileQ1 <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  for(k in 1:n2){
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
    }}
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
    linearpred[k]<-betaL[1]+b[k,1]+(betaL[indtime[1]]+b[k,2])*Time[k]+
      (betaL[indtime[2]]+b[k,3])*pow(Time[k],2)

    for(l in 1:C){

    Alpha0[k,l]<-  inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+
                 gamma1[l]*(betaL[1]+b[k,1])
    Alpha1[k,l]<- gamma1[l]*(betaL[indtime[1]]+b[k,2])
    Alpha2[k,l]<- gamma1[l]*(betaL[indtime[2]]+b[k,3])

    haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                  (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                  (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                  (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                  (h[5,l]*step(Time[k]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*Time[k]+Alpha2[k,l]*pow(Time[k],2))


    for(j in 1:K){
      #  Hazard function at Gauss-Kronrod/Legendre nodes
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
    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.001)
  }

  for(k in 1:NbetasS){
    for(l in 1:C){
      betaS[l,k]~dnorm(0,0.001)
    }}

  for(j in 1:J){
    for(l in 1:C){
      h[j,l]~dgamma(0.1,0.1)
    }}

  for(l in 1:C){
    gamma1[l]~dnorm(0,0.001)
  }
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Derive dquantity
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])



  tau1~dgamma(0.01,0.01)
  sigma1<-1/tau1
}"

    model_fileQ <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  for(k in 1:n2){
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
    }}
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
    linearpred[k]<-inprod(betaL[nindtime],Xv[k,])+b[k,1]+(betaL[indtime[1]]+b[k,2])*Time[k]+
      (betaL[indtime[2]]+b[k,3])*pow(Time[k],2)

    for(l in 1:C){

    Alpha0[k,l]<-  inprod(betaS[l,1:NbetasS],XS[k,1:NbetasS])+
                 gamma1[l]*(inprod(betaL[nindtime],Xv[k,])+b[k,1])
    Alpha1[k,l]<- gamma1[l]*(betaL[indtime[1]]+b[k,2])
    Alpha2[k,l]<- gamma1[l]*(betaL[indtime[2]]+b[k,3])

    haz[k,l]<- ((h[1,l]*step(s[1]-Time[k]))+
                  (h[2,l]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                  (h[3,l]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                  (h[4,l]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                  (h[5,l]*step(Time[k]-s[4])))*exp(Alpha0[k,l]+Alpha1[k,l]*Time[k]+Alpha2[k,l]*pow(Time[k],2))


    for(j in 1:K){
      #  Hazard function at Gauss-Kronrod/Legendre nodes
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
    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.001)
  }

  for(k in 1:NbetasS){
    for(l in 1:C){
      betaS[l,k]~dnorm(0,0.001)
    }}

  for(j in 1:J){
    for(l in 1:C){
      h[j,l]~dgamma(0.1,0.1)
    }}

  for(l in 1:C){
    gamma1[l]~dnorm(0,0.001)
  }
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Derive dquantity
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])



  tau1~dgamma(0.01,0.01)
  sigma1<-1/tau1
}"

    C <- length(unique(CR)) - 1
    set.seed(2)
    i.jags <- function() {
      list(
        gamma1 = stats::rnorm(C),
        betaL = stats::rnorm(dim(X)[2]), tau1 = 1, Omega = diag(stats::runif(Nb)),
        b = matrix(0, n2, Nb)
      )
    }

    parameters <- c(
      "betaL", "betaS", "linearpred", "b",
      "gamma1", "sigma1", "Sigma", "h"
    )
    ###############################

    d.jags1 <- list(
      n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
      X = X, Z = Z, id = id2, indtime = indtime,
      CR = CR, mub = rep(0, Nb), V = diag(1, Nb), Nb = Nb, zeros = rep(0, n2),
      NbetasL = dim(X)[2], s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, KK = KK
    )

    d.jags <- list(
      n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
      X = X, Z = Z, id = id2, Xv = Xv, indtime = indtime, nindtime = c(1:dim(X)[2])[-indtime],
      CR = CR, mub = rep(0, Nb), V = diag(1, Nb), Nb = Nb, zeros = rep(0, n2),
      NbetasL = dim(X)[2], s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, KK = KK
    )
    set.seed(2)
    if (is.matrix(Xv) == FALSE) {
      model_fileQ_final <- textConnection(model_fileQ1)
      d.jags <- d.jags1
    } else {
      model_fileQ_final <- textConnection(model_fileQ)
    }



    sim1 <- jagsUI::jags(
      data = d.jags,
      inits = i.jags,
      parameters.to.save = parameters,
      model.file = model_fileQ_final,
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
      beta = sim1$sims.list$betaL, gamma = sim1$sims.list$betaS,
      alpha = sim1$sims.list$gamma1,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma1,
      Sigma = sim1$sims.list$Sigma,
      linearpred = sim1$sims.list$linearpred,
      b = sim1$sims.list$b
    )
    PMean <- list(
      beta = sim1$mean$betaL, gamma = sim1$mean$betaS,
      alpha = sim1$mean$gamma1,
      h = sim1$mean$h,
      sigma = sim1$mean$sigma1,
      Sigma = sim1$mean$Sigma,
      linearpred = sim1$mean$linearpred,
      b = sim1$mean$b
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL) <-
        names(sim1$sd$betaL) <-
        names(sim1$q2.5$betaL) <-
        names(sim1$q97.5$betaL) <-
        names(sim1$Rhat$betaL) <- colnames(X)

      names(sim1$mean$sigma1) <-
        names(sim1$sd$sigma1) <-
        names(sim1$q2.5$sigma1) <-
        names(sim1$q97.5$sigma1) <-
        names(sim1$Rhat$sigma1) <- "sigma2"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c("Intercept", "Time", "Time2")



      LM <- rbind(
        cbind(sim1$mean$betaL, sim1$sd$betaL, sim1$q2.5$betaL, sim1$q97.5$betaL, sim1$Rhat$betaL),
        cbind(sim1$mean$sigma1, sim1$sd$sigma1, sim1$q2.5$sigma1, sim1$q97.5$sigma1, sim1$Rhat$sigma1)
      )

      Sigma <- sim1$mean$Sigma
      colnames(LM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Longitudinal_model = LM, Sigma = Sigma)
    } else {
      names(sim1$mean$betaL) <-
        names(sim1$sd$betaL) <-
        names(sim1$q2.5$betaL) <-
        names(sim1$q97.5$betaL) <- colnames(X)

      names(sim1$mean$sigma1) <-
        names(sim1$sd$sigma1) <-
        names(sim1$q2.5$sigma1) <-
        names(sim1$q97.5$sigma1) <- "sigma2"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c("Intercept", "Time", "Time2")



      LM <- rbind(
        cbind(sim1$mean$betaL, sim1$sd$betaL, sim1$q2.5$betaL, sim1$q97.5$betaL),
        cbind(sim1$mean$sigma1, sim1$sd$sigma1, sim1$q2.5$sigma1, sim1$q97.5$sigma1)
      )


      colnames(LM) <- c("Est", "SD", "L_CI", "U_CI")

      Sigma <- sim1$mean$Sigma


      results <- list(Longitudinal_model = LM, Sigma = Sigma)
    }
  }

  DIC <- sim1$DIC - 2 * KK * n2


  list(
    PMean = PMean,
    MCMC = MCMC,
    formFixed = formFixed, formRandom = formRandom, formGroup = formGroup, formSurv = formSurv,
    Estimation = results, DIC = DIC
  )
}
