#' Create a matrix of temporal weights.
#'
#' The \code{createWeightT} function creates a symmetric matrix of temporal weights
#' based on temporal autocorrelation.   At present only AR(1) structure has been
#' implemented.
#'
#' Description of AR1 correlation.
#'
#' @param n integer.  This should specify the number of time points in the data.
#' @param rho scalar.  Rho must be a value between 0 and 1 reflecting the
#'   strength of temporal autocorrelation.
#' @param corr.form character string.  At present only "AR1" is implemented for
#'   autoregressive temporal correlation of lag 1.
#' @return An \eqn{n \times n} matrix of temporal weights.
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' weightT <- createWeightT(n = 5, rho = 0.5)
#' @export
createWeightT <- function(n, rho, corr.form = "AR1") {
  if(rho <= 0 | rho >=1) {stop("rho must be between 0 and 1")}

  if(corr.form=="AR1") {
  weightT <- diag(n)
  weightT <- rho^abs(row(weightT)-col(weightT))
  }
  return(weightT)
}


