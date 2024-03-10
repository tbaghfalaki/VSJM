#' Simulated multiple longitudinal data
#'
#' Simulated data were generated in which "Y1",...,"Y10" follows a Gaussian distribution, "obstime" represents the observed time with some covariates.
#'
#'
#' @name dataLong
#' @format A data frame which contains id, obstime, x1, x2, w1, w2, and Y1,...,Y10.
#' \describe{
#'   \item{id}{patients identifier}
#'   \item{obstime}{observed times for longitudinal measerements}
#'   \item{Y1}{the first longitudinal response}
#'   \item{Y2}{the second longitudinal response}
#'   \item{Y3}{the third longitudinal response}
#'   \item{Y4}{the fourth longitudinal response}
#'   \item{Y5}{the fifth longitudinal response}
#'   \item{Y6}{the sixth longitudinal response}
#'   \item{Y7}{the seventh longitudinal response}
#'   \item{Y8}{the eighth longitudinal response}
#'   \item{Y9}{the ninth longitudinal response}
#'   \item{Y10}{the tenth longitudinal response}
#'   \item{w1}{a continuous covariate}
#'   \item{x1}{a binary covariate}
#'   \item{w2}{a continuous covariate}
#'   \item{x2}{a binary covariate}
#' }
#' @seealso \code{\link{UJM}}
"dataLong"
