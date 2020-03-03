#' Wrapper function for \code{missMDA} package to impute missing values
#'
#' A wrapper function to impute missing values using the \code{missMDA} package
#' (Josse and Husson (2012)).  The function estimates the number of principal
#' components used to impute missing values and carries out multiple imputation.
#' Diagnostic plots can also be produced if required.
#'
#' This function uses "K-fold" cross validation to estimate the number of
#' principal component to use to impute missing values.  Alternative cross
#' validation methods are available in the \code{\link[missMDA]{estim_ncpPCA}}
#' function in the \code{missMDA} package.  In order to assess uncertainty due
#' to missingness, this function also implements multiple imputation, based on
#' the \code{\link[missMDA]{MIPCA}} function in the \code{missMDA} package.  At
#' present, \code{method.mi} from \code{\link[missMDA]{MIPCA}} is set to
#' \code{"Boot"}.  Further information can be found by consulting the help
#' documentation in the \code{missMDA} package.
#'
#' \code{pca.mode} The \code{pca.mode} must be specified otherwise an error
#' message will be produced.  This can be specified as either \code{"Smode"} or
#' \code{"Tmode"}.  For S-mode, PCA is performed on a data matrix where the
#' columns are monitoring sites and the rows are regularly spaced time points.
#' For T-mode, PCA will be performed on a data matrix where the columns are
#' regularly spaced time points and the rows are monitoring sites.  See the
#' Common_Patterns vignette or Richman (1986) for more information.
#'
#' \code{ncp.min} If this equals 0 then missing values have been imputed using
#' column means.  This suggests that variables (columns) in the data are not
#' correlated.
#'
#' \code{estim.plot} This plots the number of principal components used to
#' impute missing values against the K-fold cross validation score.  The number
#' of principal components used to produce the completed data is the one
#' corresponding to the minimum cross validation score.
#'
#' \code{pNA} This is the proportion of values to remove from the observed data
#' that will be used for "Kfold" cross validation.  Further details about this
#' can be found in help file for \code{\link[missMDA]{estim_ncpPCA}} function in
#' the \code{missMDA} package.
#'
#' \code{nboot} is the number of simulated data sets to create during multiple
#' imputation, to assess variability of results due to missingness.  Further
#' information can be found in the \code{\link[missMDA]{MIPCA}} help file in the
#' \code{missMDA} package.
#'
#' \code{mi.plot} This produces the \code{"dim"} plot from the \code{choice}
#' argument in \code{\link[missMDA]{plot.MIPCA}} function in the \code{missMDA}
#' package. The plot can only be produced if the number of principal components
#' (\code{ncp}) used to calculate missing values is 2 or more.  If \code{TRUE}
#' then the plot will be produced for the first 2 principal components otherwise
#' if \code{ncp} < 2 a warning message is produced.  The function will run a bit
#' quicker if this is set to \code{FALSE} as the plot takes a few minutes to be
#' produced.
#'
#' @import missMDA
#'
#' @param x data.frame.  Missing values should be indicated by NA.  Columns MUST
#'   be monitoring sites while rows MUST be timepoints.
#' @param pca.mode character string.  MUST be specified.  \code{"Smode"} is PCA
#'   performed on data where columns are monitoring sites.  \code{"T-mode"} is
#'   for PCA performed on data where columns are time points.  See Details for
#'   further information.
#' @param ncp.min integer.  Minimum number of principal components used to
#'   estimate missing values.  Default is 0.  See Details for more information.
#' @param ncp.max integer.  Maximum number of principal components used to
#'   estimate missing values.  Default is 5.
#' @param scale logical.  TRUE will divide values in x by column standard
#'   deviation and should be used if PCA is to be performed on correlation
#'   matrix.  Default is FALSE and is equivalent to PCA performed on covariance
#'   matrix.
#' @param mean.plot logical.  If TRUE will produce a plot of the mean time
#'   series of incomplete data and mean time series of data completed using
#'   imputation.  Default is TRUE.
#' @param response.var character string.  The text for the y-axis label in
#'   \code{mean.plot}.  Default is "response (units)".  Only used if
#'   \code{mean.plot = TRUE}.
#' @param estim.plot logical.  If TRUE (default) then a plot is produced showing
#'   the cross validation score for \code{ncp.min} to \code{ncp.max} principal
#'   components.  See Details for further information.
#' @param nbsim integer.  Number of simulations for K-fold cross validation.
#' @param pNA numeric.  Must be in the range (0, 1).  Specifies the proportion
#'   of values to leave out during K-fold cross validation.  See Details for
#'   further information.
#' @param nboot integer.  Number of simulated data sets to create during
#'   multiple imputation.  See Details for further information.
#' @param mi.plot logical.  If TRUE will produce a diagnostic plot to assess
#'   variability from imputed data.  Default is FALSE.  See Details for further
#'   information.
#' @return \code{ncp} is the number of components to use to calculate missing
#'   values, estimated using K-fold cross validation.
#'
#' \code{criterion} are the cross validation values for ncp.min to ncp.max.
#'
#' \code{res.MIPCA} contains \code{res.MI = nboot} simulated data sets used to
#' assess variability from missing values and to produce \code{mi.plot},
#' \code{res.imputePCA} is the completed data set with no missing values and
#' calculated using \code{ncp} principal components, \code{call} is the call to
#' the function.
#'
#' The function will produce plots showing \code{criterion} against \code{ncp},
#' \code{mi.plot} if TRUE and \code{ncp} > 1, \code{mean.plot} if TRUE.
#'
#' @references Josse, J. and F. Husson (2012). Handling missing values in
#'   exploratory multivariate data analysis methods. Journal de la Societe
#'   Francaise de Statistique, 153(2), 79-99.
#'
#'   Richman, M. B. (1986), Rotation of principal components. J. Climatol., 6:
#'   293-335. doi:10.1002/joc.3370060305
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' library(stpca)
#'
#' ## load a demo dataset containing missing values
#' data(demoYmiss)
#'
#' ## Imputation in T-mode
#' Tmode.impute <- completeData(demoYmiss, pca.mode = "Tmode", ncp.min = 0,
#'                              ncp.max = 1, scale = FALSE, mean.plot = TRUE,
#'                              response.var ="response (units)",
#'                              estim.plot = FALSE, nbsim = 10, pNA = 0.05,
#'                              nboot = 10, mi.plot = TRUE)
#' ## Imputation in S-mode
#' Smode.impute <- completeData(x = demoYmiss, pca.mode = "Smode", ncp.min = 0,
#'                              ncp.max = 4, scale = FALSE, mean.plot = TRUE,
#'                              response.var ="response  (units)",
#'                              estim.plot = TRUE, nbsim = 10, pNA = 0.05,
#'                              nboot = 10, mi.plot = TRUE)
#'
#' @export
completeData <- function(x, pca.mode = NULL,
                         ncp.min = 0, ncp.max = 5, scale = FALSE,
                         mean.plot = TRUE, response.var="response (units)",
                         estim.plot=TRUE, nbsim = 100, pNA = 0.05,
                         nboot = 100, mi.plot = FALSE) {


  if(missing(pca.mode)) {stop("you must specify pca.mode")}

  results <- list()
  # read in data and scale as appropriate
  data <- x

  if(pca.mode=="Tmode") {
    data <- t(data)
  }

  # estimate the number of principal components to use for imputation
  num.comp <- missMDA::estim_ncpPCA(data, ncp.min=ncp.min, ncp.max=ncp.max,
                                    method = "Regularized", scale=scale,
                                    method.cv="Kfold", nbsim=nbsim, pNA=pNA)
  ncp <- num.comp$ncp

  if(estim.plot==TRUE) {

    if(length(ncp.min:ncp.max) < 2) {warning("cannot produce estim.plot when ncp.min=ncp.max")}

    if(length(ncp > 1)) {
      par(mfrow=c(1,1))
      plot(ncp.min:ncp.max, num.comp$criterion, type="b",
           xlab="number of principal components",
           ylab=paste("K-fold cross validation criterion"),
           cex.lab=1.5, cex.axis=1.5, cex=1.5, xaxt="n")
      abline(v=num.comp$ncp, col=2, lty=2)
      axis(1, at=c(ncp.min:ncp.max), labels=c(ncp.min:ncp.max), cex.axis=1.5)
      axis(1, at=num.comp$ncp, labels=num.comp$ncp, col=2, col.axis=2, cex.axis=1.5)
    }
  }

  results[["ncp"]] <- num.comp$ncp
  results[["criterion"]] <- num.comp$criterion

  # do multiple imputation of missing values
  mult.data <- missMDA::MIPCA(data, ncp=ncp, scale=scale, nboot=nboot,
                              method.mi="Boot")

  results[["res.MIPCA"]] <- mult.data

  # plot simulated data sets
  if(mi.plot==TRUE) {
    if(ncp < 2) {warning("cannot produce mi.plot when ncp < 2")}

    if(ncp >= 2) {
    par(mfrow=c(1,1))
    plot(mult.data, new.plot=FALSE, choice="dim")
    }
  }


  # plot mean incomplete time series and mean complete time series
  if(mean.plot==TRUE) {

    par(mfrow=c(1,1))
    if(pca.mode=="Smode") {
      plot(1:nrow(x), apply(x, 1, mean, na.rm=TRUE), type="l", lwd=2, col=2,
           xlab="Date", ylab=response.var,
           main="Mean time series of original and completed data",
           cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
      lines(1:nrow(x), apply(mult.data$res.imputePCA, 1, mean), lwd=2, col=1)
      legend("bottomright", legend=c("original data", "completed data"),
             lwd=c(2,2), col=c(2,1), bty="n")
    }

    if(pca.mode=="Tmode") {
      plot(1:nrow(x), apply(x, 1, mean, na.rm=TRUE), type="l", lwd=2, col=2,
           xlab="Date", ylab=response.var,
           main="Mean time series of original and completed data",
           cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
      lines(1:nrow(x), apply(mult.data$res.imputePCA, 2, mean), lwd=2, col=1)
      legend("bottomright", legend=c("original data", "completed data"),
             lwd=c(2,2), col=c(2,1), bty="n")
    }

  }

  return(results)

}



