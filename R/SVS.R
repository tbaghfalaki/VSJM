#' Summary function of VS
#'
#' @description
#' Computing dynamic prediction based on the output of summary of the outputs of the VS function, a two-stage approach to efficiently select variables for joint modeling of multiple longitudinal markers and competing risks outcomes within a Bayesian framework.
#'
#' @details
#' It presents a summary of the output of VS function, including LBFDR, BF, and parameter estimations.
#'
#'
#' @param object 	an object inheriting from class VS
#'
#'
#' @return
#' - gamma list of posterior mean for each gamma including estimation, standard deviation, LBFDR and BF
#' - alpha list of posterior mean for each alpha including estimation, standard deviation, LBFDR and BF
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#'
#' @example inst/exampleVS.R
#'
#' @md
#' @export

SVS <- function(object) {
  XS <- object$XS
  C <- object$C
  nmark <- object$nmark

  LBFDR1 <- object$LBFDRX
  BF1 <- object$BFX
  rownames(LBFDR1) <- rownames(BF1)  <- paste0("Cause ", 1:C)
  colnames(LBFDR1) <- colnames(BF1)  <- colnames(XS)


  LBFDR2 <- object$LBFDRY
  BF2 <- object$BFY
  rownames(LBFDR2) <- rownames(BF2) <- paste0("Cause ", 1:C)
  colnames(LBFDR2) <- colnames(BF2) <- paste0("Marker", 1:nmark)


  Gamma_r <- list(LBFDR = LBFDR1, BF = BF1, I_LBFDR = 1 * (LBFDR1 < 0.05), I_BF = 1 * (BF1 > 1))
  Alpha_r <- list(LBFDR = LBFDR2, BF = BF2, I_LBFDR = 1 * (LBFDR2 < 0.05), I_BF = 1 * (BF2 > 1))

  list(gamma = Gamma_r, alpha = Alpha_r)
}
