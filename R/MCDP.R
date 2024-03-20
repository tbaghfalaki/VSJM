#'  Monte-Carlo Approximation of Dynamic prediction
#'
#' @description
#' Monte-Carlo approximation ofDynamic prediction for VSJM
#'
#'
#' @details
#' Estimate Monte-Carlo approximation of DP for joint modeling based on VS
#'
#' @param object an object inheriting from class VS
#' @param object2 an object inheriting from class VS2
#' @param Method the method for variable selection including "LBFDR" for LBFDR and "BF" for Bayes factor.
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param s the landmark time for prediction
#' @param t the window of prediction for prediction
#' @param cause_main the main cause for prediction
#' @param mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
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
#' @example inst/exampleVSJM.R
#'
#' @md
#' @export

MCDP <- function(object, object2, Method = "LBFDR", s = s, t = t, cause_main=cause_main, mi=10,
                 n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
               n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
               DIC = TRUE, quiet = FALSE, dataLong, dataSurv) {
  Dt <- t
  KK <- 1000000

  formFixed <- object$formFixed
  formRandom <- object$formRandom
  formGroup <- object$formGroup
  formSurv <- object$formSurv
  model <- object$model
  Obstime <- object$Obstime
  C <- object$C
  nmark <- object$nmark
  mu1 <- object$mu1
  #######
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

  ########### univariate_jm_random_effect_estimation

  model_fileI1b <- "model{
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

  Omega<-1/Sigma

  tau1<-1/sigma1
}"

  model_fileIb <- "model{
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

 Omega<-1/Sigma

  tau1<-1/sigma1
}"


  model_fileL1b <- "model{
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


  Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1

}"

  model_fileLb <- "model{
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

   Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1
}"



  model_fileQ1b <- "model{
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


    Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1
}"

  model_fileQb <- "model{
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

    Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1
}"


  ########
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  CR <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  suppressWarnings({
    mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
  })
  XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]

  glq <- statmod::gauss.quad(15, kind = "legendre")
  xk <- glq$nodes # Nodes
  wk <- glq$weights # Weights
  K <- length(xk) # K-points
  ################

  peice <- stats::quantile(Time, seq(.2, 0.8, length = 4))
  delta <- nnet::class.ind(arules::discretize(Time, method = "fixed", c(0, peice, max(Time))))


  data_Long_s <- dataLong[dataLong$obstime <= s, ]

  DP_tot=NULL
  for(ttt in 1:mi){
  X <- Z <- Xv <- Zv <- Nb <- list()
  indB <- indtime <- list()
  bhat_mean <- bhat_chain <- list()
  for (j in c(1:nmark)[apply(I_alpha, 2, max) > 0]) {

    ###############################   ###############################  ###############################
    betaLmc <- object$sim_step1[[j]]$sim$beta
    gamma1mc <- object$sim_step1[[j]]$sim$alpha
    sigma1mc <- object$sim_step1[[j]]$sim$sigma

    Sample=sample(1:length(sigma1mc),mi)

    betaSmc <- object$sim_step1[[j]]$sim$gamma
    Sigmamc <- object$sim_step1[[j]]$sim$Sigma
    hmc <- object$sim_step1[[j]]$sim$h
    ###############################  ###############################  ###############################
    if(is.matrix(betaLmc)==TRUE){
      betaLmc=betaLmc[Sample,]
    }else{
      betaLmc=betaLmc[Sample]
    }
    gamma1mc=gamma1mc[Sample,]
    sigma1mc=sigma1mc[Sample]
    if(is.array(betaSmc)==TRUE){
      betaSmc=betaSmc[Sample,,]
    }else{
      betaSmc=betaSmc[Sample,]
    }

    if(is.array(Sigmamc)==TRUE){
      Sigmamc=Sigmamc[Sample,,]
    }else{
      Sigmamc=Sigmamc[Sample]
    }

    hmc=hmc[Sample,,]




    if(is.matrix(betaLmc)==TRUE){
      betaL=betaLmc[ttt,]
    }else{
      betaL=betaLmc[ttt]
    }

    gamma1 <- gamma1mc[ttt,]
    sigma1 <- sigma1mc[ttt]
    h <- hmc[ttt,,]



    if(is.array(betaSmc)==TRUE){
      betaS=betaSmc[ttt,,]
    }else{
      betaS=betaSmc[ttt,]
    }

    if(is.array(Sigmamc)==TRUE){
      Sigma=Sigmamc[ttt,,]
    }else{
      Sigma=Sigmamc[ttt]
    }


    ###############################  ###############################  ###############################


    if (model[[j]] == "intercept") {
      data_long <- data_Long_s[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      #data_long <- data_long[is.na(y) == FALSE, ]
      #y <- data_long[all.vars(formFixed[[j]])][, 1]
      mfX <- stats::model.frame(formFixed[[j]], data = data_long , na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX , na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long , na.action = NULL)
      id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
      n2 <- length(unique(id))
      n <- length(id)

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
      data_long <- data_Long_s[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      #data_long <- data_long[is.na(y) == FALSE, ]
      #y <- data_long[all.vars(formFixed[[j]])][, 1]

      mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
      Z[[j]] <- stats::model.matrix(formRandom[[j]], mfU, na.action = NULL)
      id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
      n <- length(id)
      M <- table(id)
      id2 <- rep(1:length(M), M)
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
      data_long <- data_Long_s[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      # data_long <- data_long[is.na(y) == FALSE, ]
      # y <- data_long[all.vars(formFixed[[j]])][, 1]


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
      M <- table(id)
      id2 <- rep(1:length(M), M)

      n <- length(id)
      n2 <- length(unique(id))
      Nb[[j]] <- dim(Z[[j]])[2]

      Obstime2n <- c(Obstime, Obstime2)
      Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime2n)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications


      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime2n] # index of time
    }
    # }






    ########  BUGS code  ########
    M1 <- table(id)
    NbetasS <- dim(XS)[2]


    if(model[[j]]=="intercept"){
      i.jags <- function() {
        list(
          b = rep(0, n2)
        )
      }
    }else{
      i.jags <- function() {
        list(
          b = matrix(0, n2, Nb[[j]])
        )
      }
    }
    parameters <- c("b")

    #### Data


    NbetasS <- dim(XS)[2]

    if (is.matrix(Xv[[j]]) == FALSE) {

      if(model[[j]]=="intercept"){
        model_fileLb_last <- textConnection(model_fileI1b)
        d.jags <- list(
          betaL = betaL, betaS = betaS,
          gamma1 = gamma1, sigma1 = sigma1,
          Sigma = Sigma, h = h,
          n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
          X = X[[j]], id = id2, indtime = indtime[[j]],
          CR = CR, zeros = rep(0, n2),
          s = peice, xk = xk, wk = wk, K = K, KK = KK
        )
      }

      if(model[[j]]=="linear"){
        model_fileLb_last <- textConnection(model_fileL1b)
        d.jags <- list(
          betaL = betaL, betaS = betaS,
          gamma1 = gamma1, sigma1 = sigma1,
          Sigma = Sigma, h = h,
          n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
          X = X[[j]], Z = Z[[j]], id = id2, indtime = indtime[[j]],
          CR = CR, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
          s = peice,  xk = xk, wk = wk, K = K, KK = KK
        )
      }
      if(model[[j]]=="quadratic"){
        model_fileLb_last <- textConnection(model_fileQ1b)
        d.jags <- list(
          betaL = betaL, betaS = betaS,
          gamma1 = gamma1, sigma1 = sigma1,
          Sigma = Sigma, h = h,
          n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
          X = X[[j]], Z = Z[[j]], id = id2, indtime = indtime[[j]],
          CR = CR, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
          s = peice,  xk = xk, wk = wk, K = K, KK = KK
        )
      }
    }else{
      if(model[[j]]=="intercept"){
        model_fileLb_last <- textConnection(model_fileIb)
        d.jags <- list(
          betaL = betaL, betaS = betaS,
          gamma1 = gamma1, sigma1 = sigma1,
          Sigma = Sigma, h = h,
          n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
          X = X[[j]], id = id2, Xv = Xv[[j]], indtime = indtime[[j]], nindtime = c(1:dim(X[[j]])[2])[-indtime[[j]]],
          CR = CR, zeros = rep(0, n2),
          s = peice, xk = xk, wk = wk, K = K, KK = KK
        )
      }
      if(model[[j]]=="linear"){
        model_fileLb_last <- textConnection(model_fileLb)
        d.jags <- list(
          betaL = betaL, betaS = betaS,
          gamma1 = gamma1, sigma1 = sigma1,
          Sigma = Sigma, h = h,
          n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
          X = X[[j]], Z = Z[[j]], id = id2, Xv = Xv[[j]], indtime = indtime[[j]], nindtime = c(1:dim(X[[j]])[2])[-indtime[[j]]],
          CR = CR, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
          s = peice, xk = xk, wk = wk, K = K, KK = KK
        )



      }
      if(model[[j]]=="quadratic"){
        model_fileLb_last <- textConnection(model_fileQb)
        d.jags <- list(
          betaL = betaL, betaS = betaS,
          gamma1 = gamma1, sigma1 = sigma1,
          Sigma = Sigma, h = h,
          n = n, Time = Time, Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2], C = C,
          X = X[[j]], Z = Z[[j]], id = id2, Xv = Xv[[j]], indtime = indtime[[j]],
          nindtime = c(1:dim(X[[j]])[2])[-indtime[[j]]],
          CR = CR,  mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
          s = peice, xk = xk, wk = wk, K = K, KK = KK
        )
      }
    }
    sim1 <- jagsUI::jags(
      data = d.jags,
      inits=i.jags,
      parameters.to.save = parameters,
      model.file = model_fileLb_last,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    bhat_mean[[j]] <- sim1$mean$b
    bhat_chain[[j]] <- sim1$sims.list$b
  }

  ###################################
  n2 <- dim(dataSurv)[1]
  sigma <- c()
  # mu1 <- matrix(0, n2, nmark)
  betaL <- b <- list()

  for (j in c(1:nmark)[apply(I_alpha, 2, max) > 0]) {
    betaL[[j]] <- object$sim_step1[[j]]$PMean$beta
    sigma <- append(sigma, object$sim_step1[[j]]$PMean$sigma)
    # mu1[,j]<- object$sim_step1[[j]]$PMean$linearpred
    b[[j]] <- bhat_mean[[j]]
  }
  indtime <- nindtime <- list()
  for (j in c(1:nmark)[apply(I_alpha, 2, max) > 0]) {
    indB <- 1:dim(X[[j]])[2]
    if (model[[j]] != "quadratic") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
    }
    if (model[[j]] == "quadratic"){
      indtime[[j]] <- indB[colnames(X[[j]]) %in% c(Obstime, Obstime2)]
    }
    nindtime[[j]] <- c(1:dim(X[[j]])[2])[-indtime[[j]]]
  }
  LP1 <- LP2 <- LP3 <- matrix(0, n2, nmark)
  for (i in 1:n2) {
    for (j in c(1:nmark)[apply(I_alpha, 2, max) > 0]) {
      if (is.matrix(Xv[[j]]) == TRUE) {
        if (model[[j]] == "intercept") {
          LP1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i]
        } else {
          LP1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i, 1]
        }
      } else {
        if (model[[j]] == "intercept") {
          LP1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i]

        } else {
          LP1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i, 1]

        }
      }
      if (model[[j]] == "intercept") {
        LP2[i, j] <- betaL[[j]][indtime[[j]][1]]
      } else {
        LP2[i, j] <- betaL[[j]][indtime[[j]][1]] + b[[j]][i, 2]
      }
      LP3[i, j] <- 0


      if (model[[j]] == "quadratic") (LP3[i, j] <- betaL[[j]][indtime[[j]][2]] + b[[j]][i, 3])
    }
  }

  ########################### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  betaS <- object2$Estimation$Survival_model$gamma$Est
  alpha <- object2$Estimation$Survival_model$alpha$Est
  h <- object2$Estimation$Survival_model$lambda$Est
  ########################### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  xk11 <- wk11 <- matrix(0, n, K)
  for (k in 1:n) {
    for (j in 1:K) {
      xk11[k, j] <- (xk[j] + 1) / 2 * s
      wk11[k, j] <- wk[j] * s / 2
    }
  }

  step <- function(x) {
    z <- 0
    if (x >= 0) (z <- 1)
    z
  }

  Alpha0 <- Alpha1 <- Alpha2 <-  CH <- matrix(0, n2, C)
  chaz <- array(0, c(n2, K, C))
  for (k in 1:n2) {
    for (l in 1:C) {
      Alpha0[k, l] <- betaS[l, ] %*% XS[k, ] + alpha[l, ] %*% LP1[k, ]
      Alpha1[k, l] <- alpha[l, ] %*% LP2[k, ]
      Alpha2[k, l] <- alpha[l, ] %*% LP3[k, ]

      for (j in 1:K) {
        chaz[k, j, l] <- ((h[1, l] * step(peice[1] - xk11[k, j])) +
                            (h[2, l] * step(xk11[k, j] - peice[1]) * step(peice[2] - xk11[k, j])) +
                            (h[3, l] * step(xk11[k, j] - peice[2]) * step(peice[3] - xk11[k, j])) +
                            (h[4, l] * step(xk11[k, j] - peice[3]) * step(peice[4] - xk11[k, j])) +
                            (h[5, l] * step(xk11[k, j] - peice[4]))) * exp(Alpha0[k, l] + Alpha1[k, l] * xk11[k, j] + Alpha2[k, l] * (xk11[k, j]^2))
      }
      CH[k, l] <- wk11[k, ] %*% chaz[k, , l]
    }
  }
  DENOM <- exp(-apply(CH, 1, sum))
  ########
  haz <- matrix(0, n2, K)
  NUM <- c()
  NUM1 <- matrix(0, n2, K)
  xk1 <- wk1 <- c()

  for (j in 1:K) {
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk1[j] <- (xk[j] + 1) / 2 * Dt + s
    wk1[j] <- wk[j] * Dt / 2
  }
  chaz <- array(0, c(n2, K, C))
  for (k in 1:n2) {
    for (j in 1:K) {
      ### Hazard for the cause of interest
      haz[k, j] <- ((h[1, cause_main] * step(peice[1] - xk1[j])) +
                      (h[2, cause_main] * step(xk1[j] - peice[1]) * step(peice[2] - xk1[j])) +
                      (h[3, cause_main] * step(xk1[j] - peice[2]) * step(peice[3] - xk1[j])) +
                      (h[4, cause_main] * step(xk1[j] - peice[3]) * step(peice[4] - xk1[j])) +
                      (h[5, cause_main] * step(xk1[j] - peice[4]))) *
        exp(Alpha0[k, l] + Alpha1[k, l] * xk1[j] + Alpha2[k, l] * (xk1[j]^2))


      ################################################

      sumint <- function(v) {
        xkn1 <- wkn1 <- c()
        CH <- matrix(0, n2, C)
        for (l in 1:C) {
          for (jj in 1:K) {
            # Scaling Gauss-Kronrod/Legendre quadrature
            xkn1[jj] <- (xk[jj] + 1) / 2 * v
            wkn1[jj] <- wk[jj] * v / 2



            chaz[k, jj, l] <- ((h[1, l] * step(peice[1] - xkn1[jj])) +
                                 (h[2, l] * step(xkn1[jj] - peice[1]) * step(peice[2] - xkn1[jj])) +
                                 (h[3, l] * step(xkn1[jj] - peice[2]) * step(peice[3] - xkn1[jj])) +
                                 (h[4, l] * step(xkn1[jj] - peice[3]) * step(peice[4] - xkn1[jj])) +
                                 (h[5, l] * step(xkn1[jj] - peice[4]))) * exp(Alpha0[k, l] + Alpha1[k, l] * xkn1[jj] + Alpha2[k, l] * (xkn1[jj]^2))
          }
          CH[k, l] <- wkn1 %*% chaz[k, , l]
        }
        exp(-apply(CH, 1, sum))
      }




      NUM1[, j] <- sumint(xk1[j]) * haz[, j]
    }

    NUM[k] <- NUM1[k, ] %*% wk1
  }
  #####################
  DP <- NUM / DENOM
  #####################
  DP_tot=rbind(DP_tot,DP)
  }

  DP=apply(DP_tot,2,mean)
  DPQ=t(apply(DP_tot,2,quantile,c(0.025,0.975)))

  DP_last=cbind(unique(id), DP,DPQ)
  colnames(DP_last)=c("id","est","lower","upper")
  DP_last=data.frame(DP_last)
  list(DP =DP_last, s = s, t = Dt)
}

