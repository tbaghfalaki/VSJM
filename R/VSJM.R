#' VSJM Package
#'
#' @description
#' Employing a two-stage strategy for variable selection in joint modeling of multiple longitudinal measures and competing risks outcomes. The initial stage involves estimating K one-marker joint models for both the competing risks and each marker. Subsequently, in the second stage, time-varying covariates are integrated into a proportional hazard model to estimate parameters. Both phases follow the Bayesian paradigm, with the latter stage using a spike-and-slab prior to enhance variable selection within the competing risks model. Additionally, this approach facilitates the computation of dynamic predictions based on the estimated parameter values.
#'
#' @docType package
#'
#' @keywords VSJM
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#' @references
#' \enumerate{
#' \item
#' Baghfalaki, T., Hashemi, R. & Jacqmin-Gadda, H. (2024). A two-stage approach for variable selection in joint modeling of multiple longitudinal markers and competing risk outcomes. Submitted.
#'
#' }
#' @name VSJM
NULL
