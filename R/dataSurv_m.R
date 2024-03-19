#' Simulated competing risks data
#'
#' Simulated survival data were generated in which "survtime" follows a survival time with a event indicator called "CR" and some covariates.
#'
#'
#' @name dataSurv_m
#' @format A data frame which contains id, survtime, CR, w1, x1, w2, x2.
#' \describe{
#'   \item{id}{patients identifier}
#'   \item{survtime}{survival time (the response variable)}
#'   \item{CR}{event indicator, 1=the first cause, 2=the second cause, 0=censored}
#'   \item{w1}{a continuous covariate}
#'   \item{x1}{a binary covariate}
#'   \item{w2}{a continuous covariate}
#'   \item{x2}{a binary covariate}
#' }
#' @seealso \code{\link{UJM}}
"dataSurv_m"
