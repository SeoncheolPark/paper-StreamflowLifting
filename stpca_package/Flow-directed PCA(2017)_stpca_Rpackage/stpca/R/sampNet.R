#' Create sampled networks.
#'
#' Investigate the effect of reducing the size of the monitoring network on
#' covariance function parameters and predictions based on reduced networks.
#' This function allows you to sample networks using random, stratified, and
#' weighted sampling schemes.  This function relies heavily on the
#' \code{\link[SSN]{SSN}}, geoR (\url{http://www.leg.ufpr.br/geoR}), and
#' sampling
#' (\url{https://cran.r-project.org/web/packages/sampling/sampling.pdf})
#' packages.
#'
#' This function creates \code{nsim} sampled networks and fits a model with a
#' spatial covariance function to each sampled network.  Various model outputs
#' and predicted values based on the spatial covariance function are stored.
#'
#' If data are from a \code{\link[SSN]{SpatialStreamNetwork-class}} object the
#' spatial covariance model with Epanechnikov Tail-up, Gaussian Euclidean, and
#' nugget components is fitted to each sampled network using
#' \code{\link[SSN]{glmssn}} from the \code{\link[SSN]{SSN}} package.  If the data are
#' from a \code{data.frame} then an exponential spatial covariance model is fitted
#' using \code{\link[geoR]{likfit}} in the geoR package.
#' For both types of data, the spatial covariance model is fitted using REML.
#'
#' If the data are from a \code{data.frame} then initialization values are
#' required to fit the spatial covariance model.  The \code{\link[geoR]{variog}}
#' and \code{\link[geoR]{variofit}} functions from the geoR package are used to
#' estimate the initialization values using all default options, and specifying
#' \code{cov.model = "exponential"} in \code{\link[geoR]{variofit}}.
#'
#' The spatial covariance model is fitted assuming no spatial trend and so for
#' the \code{\link[SSN]{glmssn}} model, this is equivalent of \code{formula =
#' response.col ~ 1}.  For the \code{\link[geoR]{likfit}} model, this is
#' equivalent to specifying \code{trend = "cte"}.
#'
#' The number of monitoring sites in a sampled network is equal to the sum of
#' sites to be included in each stratum if "stratified" is included in
#' samp.scheme.  Otherwise, the number of monitoring sites in a sampled network
#' is \code{ceiling(nrow(data) * prop.sites)} where \code{data} is the name of
#' the dataframe. The number of monitoring sites for each value in
#' \code{prop.sites} is calculated using some rounding and so there can be
#' differences between the total number of monitoring sites in "stratified" and
#' "random" for the same \code{prop.sites} value.  Because of this, the number
#' of sites to be included in a sampled network is set to be the total included
#' in the stratified samples so that the same number of sites are included for
#' all sampling schemes.  This allows comparisons to be made between sampling
#' schemes.
#'
#' \code{samp.scheme} is a character vector and can take any combination of
#' "random" for random sampling, "stratified" for stratified sampling, or
#' "weighted" for weighted sampling.  In "random" sampling, each monitoring site
#' has the same probability of being included in the sampled network.  For
#' "stratified" sampling, proportional stratification is used so that the
#' sampled network has the same proportion of sites within each strata as the
#' full network.  \code{strat.col} is the column in the data object containing
#' the variable used for stratification.  This column should be a factor
#' variable. For "weighted" sampling, a column should be included in the data
#' object containing weights in the range [0,1] reflecting the probability of a
#' monitoring site being included in the sampled network.  If weights are not in
#' the range [0,1] then the function will stop running and a warning message is
#' produced saying that \code{weight.col} does not contain appropriate values.
#' It is recommended that if any of the weights are 0 then a small constant
#' (0.001 say) is added, unless the user really does not want this monitoring
#' site included in the sampled network.  Likewise, if any of the weights are 1
#' then a small constant should be subtracted, unless the user wants this site
#' to be included in every sampled network.
#'
#' \code{sub.filepath} is a location to temporarily store sampled networks.  For
#' each \code{nsim}, a new \code{\link[SSN]{SpatialStreamNetwork-class}} object
#' is created containing the monitoring sites included in the sampled network. A
#' new SSN object is created and stored in \code{sub.filepath} for each
#' \code{nsim} and removed once the SSN model has been fitted and various
#' outputs calculated and stored.
#'
#' \code{print.sim} can be switched off if you do not require a message to be
#' printed after each sampled network is created.  The message will show you
#' which comination of `samp.scheme` 100*`prop.sites` `nsim` has just been
#' evaluated.  The other part of the message that is printed is a result of
#' importing the sampled SSN object using \code{\link[SSN]{importSSN}} within
#' the \code{sampNet} function (this part cannot be suppressed).
#'
#' @import SSN
#' @importFrom sampling strata
#' @importFrom geoR as.geodata
#' @importFrom geoR variog
#' @importFrom geoR variofit
#' @importFrom geoR likfit
#' @importFrom geoR krige.conv
#' @importFrom geoR krige.control
#' @importFrom geoR output.control
#' @importFrom geoR xvalid
#'
#' @param x	an object of class \code{"data.frame"} or
#'   \code{\link[SSN]{SpatialStreamNetwork-class}}
#' @param nsim	integer.  Specifies the number of sampled networks to be created.
#' @param prop.sites	numeric.  Specifies a vector containing the proportions of
#'   monitoring sites to retain.  Must only take values between 0 and 1.
#' @param samp.scheme	character vector.  Specifies the type of sampling scheme
#'   to use.  Can include "random", "stratified", or "weighted".  See Details
#'   for more information.
#' @param sub.filepath	character string.  Specifies the filepath to temporarily
#'   store sampled networks. Only used if \code{x} is a
#'   \code{\link[SSN]{SpatialStreamNetwork-class}} object.  See Details for
#'   further information.
#' @param predpts	character string.  Specifies the name of the dataframe in the
#'   SSN object containing prediction locations.  Only used if \code{x} is a
#'   \code{\link[SSN]{SpatialStreamNetwork-class}} object.
#' @param response.col	character string.  Specifies the name of the column in
#'   the dataframe of observed values in the SSN object that contains the
#'   response variable.  Only used if \code{x} is a
#'   \code{\link[SSN]{SpatialStreamNetwork-class}} object
#' @param addfunccol	character string.  Contains the name of the column in the
#'   dataframe of observed values in the SSN objcet containing the additive
#'   function values.  Only used if \code{x} is a
#'   \code{\link[SSN]{SpatialStreamNetwork-class}} object.
#' @param siteID	character string.  Specifies the name of the column in the
#'   data.frame or SSN object that contains monitoring site ID's.
#' @param strat.col	character string.  Specifies the name of the column in the
#'   data.frame or SSN object that contains the stratum information.  Only used
#'   when stratified sampling scheme is implemented.  See Details for further
#'   information.
#' @param weight.col	character string.  Specifies the name of the column in the
#'   data.frame or SSN object that contains the values used for weighted
#'   sampling.  Only used when weighted sampling scheme is implemented.  See
#'   Details for further information.
#' @param coords.col a vector with the column numbers corresponding to the
#'   spatial coordinates.  Only used if \code{x} is of class \code{"data.frame"}.
#' @param data.col a scalar with column number corresponding to the column in
#'   which the observed data values are stored.   Only used if \code{x} is of
#'   class \code{"data.frame"}
#' @param preds.coords an N x 2 matrix or data-frame with the 2-D coordinates of
#'   the N prediction location.   Only used if \code{x} is of class
#'   \code{"data.frame"}
#' @param print.sim logical.  If TRUE (default) then a message will be printed
#'   after each sampled network is created.  See Details for further information.
#'
#' @return Returns a list with results for each combination of
#'   \code{samp.scheme} and \code{prop.sites}.  The name of each list entry
#'   contains the sampling scheme and percentage of sites retained in the
#'   sampled networks.  For example "random.50" means that \code{samp.scheme =
#'   "random"} and \code{prop.sites = 0.5}.  Each of these will contain the
#'   following:
#'
#' \code{samp.site} is a matrix with \code{nsim} columns, each of which contains
#' the ID's of monitoring sites included in each sampled network.
#'
#' \code{cov.params} is a matrix with \code{nsim} rows and 5 columns if \code{x}
#' is an SSN object: Tail-up partial sill, Tail-up range, Euclidean partial
#' sill, Euclidean partial range, nugget.  If \code{x} is of class
#' \code{"data.frame"} then this will have 3 columns: range, partial sill, and
#' nugget.
#'
#' \code{crossval.pred} is a matrix with \code{nsim} columns and the number of
#' rows is the  number of monitoring sites in each sampled network.  Values are
#' predicted values from leave one out cross validation where the spatial
#' covariance model is fitted to all of the monitoring sites except one, and the
#' response at the omitted site is estimated from the fitted model.  See
#' \code{\link[SSN]{CrossValidationSSN}} if \code{x} is an SSN object or
#' \code{\link[geoR]{xvalid}} if \code{x} is of class \code{"data.frame"} for
#' further details.
#'
#' \code{crossval.se} is a matrix with \code{nsim} columns and number of rows is
#' equal to the number of monitoring sites in a sampled network.  Values are the
#' standard errors of \code{crossval.pred}.
#'
#' \code{cross.stats} a matrix containing root mean square prediction error
#' (RMSPE).  RMSPE is calculated as:
#'
#' \deqn{\sqrt{\frac{1}{N}\sum^N_{i=1}\left(\hat{Y}(s_i)-Y(s_i)\right)^2}}
#'
#' where \eqn{N} is the number of monitoring sites, \eqn{\hat{Y}(s_i)} is the
#' value at site \eqn{i} predicted from a model fitted to \eqn{\mathbf{s}_{-i}},
#' data from all monitoring sites except site \eqn{i}, and \eqn{Y(s_i)} is the
#' observed value at site \eqn{i}.  This is also known as Leave One Out Cross
#' Validation (LOOCV).  This is a measure of how well the model fits the data.
#' Low prediction error is desireable.
#'
#' \code{preds.value} is a matrix with \code{nsim} columns and number of rows
#' equal to the number of prediction locations.  Values are predicted response
#' at unobserved locations based on the spatial covariance function estimated
#' for sampled network and were obtained using \code{\link[SSN]{predict.glmssn}}
#' if \code{x} is an SSN object, or \code{\link[geoR]{xvalid}} if \code{x} is of
#' class \code{"data.frame"}.
#'
#' \code{preds.se} is a matrix with \code{nsim} columns and number of rows equal
#' to the number of prediction locations.  Values are standard error of
#' predicted values for \code{preds.value}.  See
#' \code{\link[SSN]{predict.glmssn}} if \code{x} is an SSN object, or
#' \code{\link[geoR]{xvalid}} if \code{x} is of class \code{"data.frame"} for
#' further details.  (Note that if \code{x} is of class \code{"data.frame"} then
#' the value here is the square root of \code{krige.var} from
#' \code{\link[geoR]{krige.conv}}).
#'
#' Other results, all of which are used in \code{\link{plot_sampNet}} include:
#'
#' \code{results.names} a character string containing all combinations of
#' \code{samp.scheme} and \code{prop.sites}.
#'
#' \code{prop.sites} is the vector of values entered for the \code{prop.sites}
#' argument.
#'
#' \code{samp.scheme} is the vector of values entered for the \code{samp.scheme}
#' argument.
#'
#' \code{all} contains various outputs related to fitting a spatial covariance
#' model to the full monitoring network.  Outputs are model details
#' (\code{all.model}), covariance function parameter estimates
#' (\code{all.covparams}), cross validation statistics (\code{all.cross.stats}),
#' and predicted values with corresponding standard errors (\code{all.preds}).
#'
#' \code{columns} is a vector containing \code{response.col} and the column
#' containing standard errors of predicted values.
#'
#' \code{call} is a list containing the argument values specified in the
#' function.
#'
#' @references
#' \code{geoR}:
#'
#' Paulo J. Ribeiro Jr and Peter J. Diggle (2015). geoR: Analysis of
#' Geostatistical Data. Rpackage version 1.7-5.1.
#' \url{https://CRAN.R-project.org/package=geoR}
#'
#' \code{sampling}:
#'
#' Yves Tille and Alina Matei (2015). sampling: Survey Sampling. R package
#' version 2.7.\url{https://CRAN.R-project.org/package=sampling}
#'
#' \code{SSN}:
#'
#'   Ver Hoef, J. M. and Peterson, E. E. (2010) A moving average approach for
#'   spatial statistical models of stream networks (with discussion). Journal of
#'   the American Statistical Association 105:6-18.
#'   \url{http://dx.doi.org/10.1198/jasa.2009.ap08248}
#'
#'   Jay M. Ver Hoef, Erin E. Peterson, David Clifford, Rohan Shah (2014). SSN:
#'   An R Package for Spatial Statistical Modeling on Stream Networks. Journal
#'   of Statistical Software, 56(3), 1-43. \url{http://www.jstatsoft.org/v56/i03/.}
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' library(stpca)
#' library(SSN)
#'
#' ## get the data
#' x <- importSSN(system.file("demoSSN/demoNet.ssn", package = "stpca"),
#'                predpts = "preds")
#'
#' ## x is an SSN object
#' test.samp <- sampNet(x = x, nsim = 3, prop.sites = c(0.5),
#'                      samp.scheme = c("random"),
#'                      sub.filepath = paste(tempdir(),"/subset1.ssn", sep = ""),
#'                      predpts = "preds", response.col = "Sim_Values",
#'                      addfunccol = "addfunccol", siteID = "locID")
#'
#'
#' ## x is of class "data.frame"
#' ## in order to demonstrate this we need to first extract
#' ## the data from the SSN object so that the data are stored
#' ## in a data.frame.  sampNet() is then applied to data in a
#' ## data.frame rather than an SSN object as in the previous
#' ## example.
#'
#' # get the data.frame containing observed values
#' data(demoNet)
#' x <- getSSNdata.frame(demoNet, "Obs")
#'
#' # get the data.frame containing prediction locations
#' preds.coords <- getSSNdata.frame(demoNet, "preds")[,c("NEAR_X", "NEAR_Y")]
#'
#' test.samp <- sampNet(x, nsim = 5, prop.sites = c(0.8, 0.7),
#'                      samp.scheme = c("random", "stratified", "weighted"),
#'                      data.col = "Sim_Values", coords.col = c(9, 10),
#'                      preds.coords = preds.coords, siteID = "locID",
#'                      strat.col = "strata", weight.col = "weight")
#'
#' @export
sampNet <- function(x, nsim, prop.sites, samp.scheme = c("random"),
                    sub.filepath, predpts, response.col,
                    addfunccol, siteID, strat.col, weight.col,
                    coords.col, data.col, preds.coords, print.sim = TRUE) {

    if(class(x) == "SpatialStreamNetwork") {

    formula <- as.formula(paste(response.col, "~", 1))

    # fit SSN model to all data
    all <- list()
    all.model <- SSN::glmssn(formula , ssn.object=x,
                             CorModels=c("Epanech.tailup", "Gaussian.Euclid"),
                             addfunccol=addfunccol, EstMeth="REML")
    all.covparams <- SSN::covparms(all.model)
    all.varcomps <- SSN::varcomp(all.model)
    all.cross.stats <- SSN::CrossValidationStatsSSN(all.model)
    all.preds <- SSN::predict.glmssn(all.model, predpts)

    # store results
    all[["all.model"]] <- all.model
    all[["all.covparams"]] <- all.covparams
    all[["all.varcomps"]] <- all.varcomps
    all[["all.cross.stats"]] <- all.cross.stats
    all[["all.preds"]] <- all.preds
    all[["class"]] <- class(x)

    # data for subsets
    obs.data <- SSN::getSSNdata.frame(x, "Obs")


    results <- list()
    results.names <- paste(rep(samp.scheme, each=length(prop.sites)),
                           rep(prop.sites*100, times=length(samp.scheme)),sep=".")

    n.obs <- nrow(obs.data)
    n.obs1 <- rep(ceiling(prop.sites*n.obs), times=length(samp.scheme))




    if("stratified" %in% samp.scheme) {
      n.strat <- numeric(length(prop.sites))

      for(k in 1:length(prop.sites)) {
        stratdata <- obs.data[,c(siteID, strat.col)]
        stratdata <- stratdata[order(stratdata[,strat.col]),]
        size  <- summary(stratdata[,strat.col])
        prop <- round(size/n.obs,4)
        n <- ceiling(n.obs*as.numeric(strsplit(results.names,"[.]")[[k]][2])/100)
        n.strat[k] <- sum(round(n*prop))
      }
      n.obs1 <- rep(n.strat, times=length(samp.scheme))
    }

    # do  the simulations
    for(j in 1:length(results.names)) {

      n.obs2 <- n.obs1[j]

      prop.sites1 <- as.numeric(strsplit(results.names, "[.]")[[j]][2])/100

      samp.site <- matrix(0, nrow=n.obs2, ncol=nsim)
      colnames(samp.site) <- c(paste("sim", 1:nsim, sep=""))

      cov.params <- matrix(0, nrow=nsim, ncol=5)
      colnames(cov.params) <- c("TUparsill", "TUrange", "EucParsill", "EucRange", "Nugget")

      var.comps <- matrix(0, nrow=nsim, ncol=4)
      colnames(var.comps) <- c("Rsq", "TU", "Euc", "Nugget")

      crossval.pred <- matrix(0, nrow=n.obs2, ncol=nsim)
      colnames(crossval.pred) <- c(paste("sim", 1:nsim, sep=""))

      crossval.se <- matrix(0, nrow=n.obs2, ncol=nsim)
      colnames(crossval.se) <- c(paste("sim", 1:nsim, sep=""))

      cross.stats <- matrix(0, nrow=nsim, ncol=4)
      colnames(cross.stats) <- c("RMSPE", "cov80", "cov90", "cov.95")

      pred.n <- dim(getSSNdata.frame(x, predpts))[1]

      preds.value <- matrix(0, nrow=pred.n, ncol=nsim)
      colnames(preds.value) <- c(paste("sim", 1:nsim, sep=""))

      preds.se <- matrix(0, nrow=pred.n, ncol=nsim)
      colnames(preds.se) <- c(paste("sim", 1:nsim, sep=""))

      if(strsplit(results.names, "[.]")[[j]][1] == "random") {

        for (i in 1:nsim) {
          # copy ssn object so don't need to load data each time
          x2 <- x

          # create subset of data
          samp.sites <- sample(1:n.obs, size=n.obs2, replace=FALSE)
          obs.data$rows <- c(1:n.obs)
          obs.data$sample <- obs.data$rows %in% samp.sites

          x2 <- SSN::putSSNdata.frame(obs.data, x2, "Obs")
          sub <- SSN::subsetSSN(x2, filename=sub.filepath,
                                subset = sample == TRUE)

          samp.data <- SSN::importSSN(sub.filepath,
                                      predpts=predpts)
          samp.obs <- SSN::getSSNdata.frame(samp.data, "Obs")
          samp.obs$netID2 <- factor(samp.obs$netID,levels=sort(as.numeric(levels(samp.obs$netID))))

          drops <- "netID"
          samp.obs <- samp.obs[,!names(samp.obs) %in% drops]
          names(samp.obs)[names(samp.obs)=="netID2"] <- "netID"
          samp.data <- SSN::putSSNdata.frame(samp.obs, samp.data, "Obs")

          SSN::createDistMat(samp.data, predpts=predpts, amongpreds=TRUE, o.write=TRUE)

          # fit model to subset of data
          model <- SSN::glmssn(formula , ssn.object=samp.data,
                               CorModels=c("Epanech.tailup", "Gaussian.Euclid"),
                               addfunccol=addfunccol, EstMeth="REML")

          samp.site[,i] <- samp.sites

          cov.params[i,] <- SSN::covparms(model)[,3]

          var.comps[i,] <- SSN::varcomp(model)[,2]

          cv <- SSN::CrossValidationSSN(model)
          crossval.pred[,i] <- cv[,1]
          crossval.se[,i] <- cv[,2]

          cvstats <- SSN::CrossValidationStatsSSN(model)
          cross.stats[i,1] <- cvstats[1,3]
          cross.stats[i,2] <- cvstats[1,6]
          cross.stats[i,3] <- cvstats[1,7]
          cross.stats[i,4] <- cvstats[1,8]

          predictions <- SSN::predict.glmssn(model, predpointsID=predpts)
          predsdata <- SSN::getSSNdata.frame(predictions, predpts)
          preds.value[,i] <- predsdata[,response.col]
          predsSEcol <- paste(response.col, "predSE", sep=".")
          preds.se[,i] <- predsdata[,predsSEcol]

          unlink(sub.filepath, recursive=TRUE)

          if(print.sim==TRUE) {
            print(paste(strsplit(results.names, "[.]")[[j]][1],
                        strsplit(results.names, "[.]")[[j]][2], i))
          }


        }
      }


      if(strsplit(results.names, "[.]")[[j]][1] == "stratified") {

        for (i in 1:nsim) {
          # copy ssn object so don't need to load data each time
          x2 <- x

          # create subset of data
          stratdata <- obs.data[,c(siteID, strat.col)]
          stratdata <- stratdata[order(stratdata[,strat.col]),]
          stratdata[,strat.col] <- as.factor(stratdata[,strat.col])
          size  <- summary(stratdata[,strat.col])
          prop <- round(size/n.obs,4)
          n <- ceiling(n.obs*as.numeric(strsplit(results.names,"[.]")[[j]][2])/100)
          ni <- round(n*prop)

          samp.sites <- sampling::strata(stratdata, stratanames=strat.col,
                                         size=ni,
                                         method="srswor")
          samp.sites2 <- stratdata[samp.sites$ID_unit, siteID]
          obs.data$sample <- obs.data[, siteID] %in% samp.sites2

          x2 <- SSN::putSSNdata.frame(obs.data, x2, "Obs")
          sub <- SSN::subsetSSN(x2, filename=sub.filepath,
                                subset = sample == TRUE)

          samp.data <- SSN::importSSN(sub.filepath,
                                      predpts=predpts)
          samp.obs <- SSN::getSSNdata.frame(samp.data, "Obs")
          samp.obs$netID2 <- factor(samp.obs$netID,levels=sort(as.numeric(levels(samp.obs$netID))))

          drops <- "netID"
          samp.obs <- samp.obs[,!names(samp.obs) %in% drops]
          names(samp.obs)[names(samp.obs)=="netID2"] <- "netID"
          samp.data <- SSN::putSSNdata.frame(samp.obs, samp.data, "Obs")

          SSN::createDistMat(samp.data, predpts=predpts, amongpreds=TRUE, o.write=TRUE)

          # fit model to subset of data
          model <- SSN::glmssn(formula , ssn.object=samp.data,
                               CorModels=c("Epanech.tailup", "Gaussian.Euclid"),
                               addfunccol=addfunccol, EstMeth="REML")

          samp.site[,i] <- samp.sites2

          cov.params[i,] <- SSN::covparms(model)[,3]

          var.comps[i,] <- SSN::varcomp(model)[,2]

          cv <- SSN::CrossValidationSSN(model)
          crossval.pred[,i] <- cv[,1]
          crossval.se[,i] <- cv[,2]

          cvstats <- SSN::CrossValidationStatsSSN(model)
          cross.stats[i,1] <- cvstats[1,3]
          cross.stats[i,2] <- cvstats[1,6]
          cross.stats[i,3] <- cvstats[1,7]
          cross.stats[i,4] <- cvstats[1,8]

          predictions <- SSN::predict.glmssn(model, predpointsID=predpts)
          predsdata <- SSN::getSSNdata.frame(predictions, predpts)
          preds.value[,i] <- predsdata[,response.col]
          predsSEcol <- paste(response.col, "predSE", sep=".")
          preds.se[,i] <- predsdata[,predsSEcol]

          unlink(sub.filepath, recursive=TRUE)

          if(print.sim==TRUE) {
            print(paste(strsplit(results.names, "[.]")[[j]][1],
                        strsplit(results.names, "[.]")[[j]][2], i))
          }


        }
      }


      if(strsplit(results.names, "[.]")[[j]][1] == "weighted") {

        for (i in 1:nsim) {
          # copy ssn object so don't need to load data each time
          x2 <- x
          x2.obs <- getSSNdata.frame(x2, "Obs")

          if(sum(x2.obs[,weight.col] >= 0 & sum(x2.obs[,weight.col] <= 1)) < nrow(x2.obs)) {
            stop("weights must be in the range [0,1]")
          }

          # create subset of data
          samp.sites <- sample(x2.obs[, siteID],
                               size=ceiling(n.obs2),
                               replace=FALSE,
                               prob=(obs.data[, weight.col]))
          obs.data$sample <- obs.data[, siteID] %in% samp.sites

          x2 <- SSN::putSSNdata.frame(obs.data, x2, "Obs")
          sub <- SSN::subsetSSN(x2, filename=sub.filepath,
                                subset = sample == TRUE)

          samp.data <- SSN::importSSN(sub.filepath,
                                      predpts=predpts)
          samp.obs <- SSN::getSSNdata.frame(samp.data, "Obs")
          samp.obs$netID2 <- factor(samp.obs$netID,levels=sort(as.numeric(levels(samp.obs$netID))))

          drops <- "netID"
          samp.obs <- samp.obs[,!names(samp.obs) %in% drops]
          names(samp.obs)[names(samp.obs)=="netID2"] <- "netID"
          samp.data <- SSN::putSSNdata.frame(samp.obs, samp.data, "Obs")

          SSN::createDistMat(samp.data, predpts=predpts, amongpreds=TRUE, o.write=TRUE)

          # fit model to subset of data
          model <- SSN::glmssn(formula , ssn.object=samp.data,
                               CorModels=c("Epanech.tailup", "Gaussian.Euclid"),
                               addfunccol=addfunccol, EstMeth="REML")

          samp.site[,i] <- samp.sites

          cov.params[i,] <- SSN::covparms(model)[,3]

          var.comps[i,] <- SSN::varcomp(model)[,2]

          cv <- SSN::CrossValidationSSN(model)
          crossval.pred[,i] <- cv[,1]
          crossval.se[,i] <- cv[,2]

          cvstats <- SSN::CrossValidationStatsSSN(model)
          cross.stats[i,1] <- cvstats[1,3]
          cross.stats[i,2] <- cvstats[1,6]
          cross.stats[i,3] <- cvstats[1,7]
          cross.stats[i,4] <- cvstats[1,8]

          predictions <- SSN::predict.glmssn(model, predpointsID=predpts)
          predsdata <- SSN::getSSNdata.frame(predictions, predpts)
          preds.value[,i] <- predsdata[,response.col]
          predsSEcol <- paste(response.col, "predSE", sep=".")
          preds.se[,i] <- predsdata[,predsSEcol]

          unlink(sub.filepath, recursive=TRUE)

          if(print.sim==TRUE) {
            print(paste(strsplit(results.names, "[.]")[[j]][1],
                        strsplit(results.names, "[.]")[[j]][2], i))
          }


        }
      }


      results[[results.names[j]]] <- list(samp.site = samp.site,
                                          cov.params = cov.params,
                                          var.comps = var.comps,
                                          crossval.pred = crossval.pred,
                                          crossval.se = crossval.se,
                                          cross.stats = cross.stats,
                                          preds.value = preds.value,
                                          preds.se = preds.se)

    }
    results[["results.names"]] <- results.names
    results[["prop.sites"]] <- prop.sites
    results[["samp.scheme"]] <- samp.scheme
    results[["all"]] <- all
    results[["columns"]] <- c(response.col, predsSEcol)
    results[["call"]] <- as.list(match.call())
    return(results)

  }

  if(class(x) == "data.frame") {
    all <- list()

    all.geodata <- geoR::as.geodata(x, coords.col=coords.col, data.col=data.col)

    # get ini.cov.pars
    all.variog <- geoR::variog(all.geodata, messages=FALSE)
    all.variofit <- geoR::variofit(all.variog, cov.model="exponential", messages=FALSE)

    all.likfit <- geoR::likfit(all.geodata, cov.model="exponential",
                               lik.method="REML",
                               ini.cov.pars=all.variofit,
                               messages=FALSE)
    all.covparams <- matrix(nrow=1, ncol=3)
    all.covparams[1,1:2] <- all.likfit$cov.pars
    all.covparams[,3] <- all.likfit$nugget
    colnames(all.covparams) <- c("parsill", "range", "nugget")
    all.cross.stats <- geoR::xvalid(all.geodata, model=all.likfit, messages=FALSE)
    all.preds <- geoR::krige.conv(all.geodata, locations=preds.coords,
                                  krige=geoR::krige.control(obj.model=all.likfit),
                                  output=geoR::output.control(messages=FALSE))


    # store results
    all[["all.likfit"]] <- all.likfit
    all[["all.covparams"]] <- all.covparams
    all[["all.cross.stats"]] <- all.cross.stats
    all[["all.preds"]] <- all.preds
    all[["class"]] <- class(x)

    obs.data <- x

    # do simulations
    results <- list()
    results.names <- paste(rep(samp.scheme, each=length(prop.sites)),
                           rep(prop.sites*100, times=length(samp.scheme)),sep=".")

    n.obs <- nrow(obs.data)
    n.obs1 <- rep(ceiling(prop.sites*n.obs), times=length(samp.scheme))




    if("stratified" %in% samp.scheme) {
      n.strat <- numeric(length(prop.sites))

      for(k in 1:length(prop.sites)) {
        stratdata <- obs.data[,c(siteID, strat.col)]
        stratdata <- stratdata[order(stratdata[,strat.col]),]
        size  <- summary(stratdata[,strat.col])
        prop <- round(size/n.obs,4)
        n <- ceiling(n.obs*as.numeric(strsplit(results.names,"[.]")[[k]][2])/100)
        n.strat[k] <- sum(round(n*prop))
      }
      n.obs1 <- rep(n.strat, times=length(samp.scheme))
    }




    for(j in 1:length(results.names)) {

      n.obs2 <- n.obs1[j]

      prop.sites1 <- as.numeric(strsplit(results.names, "[.]")[[j]][2])/100

      samp.site <- matrix(0, nrow=n.obs2, ncol=nsim)
      colnames(samp.site) <- c(paste("sim", 1:nsim, sep=""))

      cov.params <- matrix(0, nrow=nsim, ncol=3)
      colnames(cov.params) <- c("parsill", "range", "Nugget")

      crossval.pred <- matrix(0, nrow=n.obs2, ncol=nsim)
      colnames(crossval.pred) <- c(paste("sim", 1:nsim, sep=""))

      crossval.se <- matrix(0, nrow=n.obs2, ncol=nsim)
      colnames(crossval.se) <- c(paste("sim", 1:nsim, sep=""))

      cross.stats <- matrix(0, nrow=nsim, ncol=1)
      colnames(cross.stats) <- "RMSPE"

      pred.n <- dim(preds.coords)[1]

      preds.value <- matrix(0, nrow=pred.n, ncol=nsim)
      colnames(preds.value) <- c(paste("sim", 1:nsim, sep=""))

      preds.se <- matrix(0, nrow=pred.n, ncol=nsim)
      colnames(preds.se) <- c(paste("sim", 1:nsim, sep=""))

      if(strsplit(results.names, "[.]")[[j]][1] == "random") {

        for (i in 1:nsim) {

          # create subset of data
          samp.sites <- sample(1:n.obs, size=n.obs2, replace=FALSE)

          samp.data <- obs.data[samp.sites,]

          # fit model to subset of data
          samp.geodata <- geoR::as.geodata(samp.data, coords.col=coords.col, data.col=data.col)

          # get ini.cov.pars
          samp.variog <- geoR::variog(samp.geodata, messages=FALSE)
          samp.variofit <- geoR::variofit(samp.variog, cov.model="exponential", messages=FALSE)

          samp.likfit <- geoR::likfit(samp.geodata, cov.model="exponential", lik.method="REML",
                                      ini.cov.pars=samp.variofit,
                                      messages=FALSE)

          samp.site[,i] <- samp.sites

          cov.params[i,] <- c(samp.likfit$cov.pars, samp.likfit$nugget)

          cv <- geoR::xvalid(samp.geodata, model=samp.likfit,
                             message=FALSE)
          crossval.pred[,i] <- cv$predicted
          crossval.se[,i] <- cv$std.error
          cross.stats[i,1] <- sqrt((sum(cv$error^2))/length(cv))

          predictions <- geoR::krige.conv(samp.geodata, locations=preds.coords,
                                          krige=geoR::krige.control(obj.model=samp.likfit),
                                          output=geoR::output.control(messages=FALSE))
          preds.value[,i] <- predictions$predict
          preds.se[,i] <- sqrt(predictions$krige.var)

          if(print.sim==TRUE) {
            print(paste(strsplit(results.names, "[.]")[[j]][1],
                        strsplit(results.names, "[.]")[[j]][2], i))
          }


        }
      }


      if(strsplit(results.names, "[.]")[[j]][1] == "stratified") {

        for (i in 1:nsim) {

          # create subset of data
          stratdata <- obs.data[,c(siteID, strat.col)]
          stratdata <- stratdata[order(stratdata[,strat.col]),]
          stratdata[,strat.col] <- as.factor(stratdata[,strat.col])
          size  <- summary(stratdata[,strat.col])
          prop <- round(size/n.obs,4)
          n <- ceiling(n.obs*as.numeric(strsplit(results.names,"[.]")[[j]][2])/100)
          ni <- round(n*prop)

          samp.sites <- sampling::strata(stratdata, stratanames=strat.col,
                                         size=ni,
                                         method="srswor")

          samp.data <- obs.data[samp.sites$ID_unit,]

          # fit model to subset of data
          samp.geodata <- geoR::as.geodata(samp.data, coords.col=coords.col, data.col=data.col)

          # get ini.cov.pars
          samp.variog <- geoR::variog(samp.geodata, messages=FALSE)
          samp.variofit <- geoR::variofit(samp.variog, cov.model="exponential", messages=FALSE)

          samp.likfit <- geoR::likfit(samp.geodata, cov.model="exponential", lik.method="REML",
                                      ini.cov.pars=samp.variofit,
                                      messages=FALSE)

          samp.site[,i] <- samp.sites$ID_unit

          cov.params[i,] <- c(samp.likfit$cov.pars, samp.likfit$nugget)

          cv <- geoR::xvalid(samp.geodata, model=samp.likfit,
                             message=FALSE)
          crossval.pred[,i] <- cv$predicted
          crossval.se[,i] <- cv$std.error
          cross.stats[i,1] <- sqrt((sum(cv$error^2))/length(cv))

          predictions <- geoR::krige.conv(samp.geodata, locations=preds.coords,
                                          krige=geoR::krige.control(obj.model=samp.likfit),
                                          output=geoR::output.control(messages=FALSE))
          preds.value[,i] <- predictions$predict
          preds.se[,i] <- sqrt(predictions$krige.var)

          if(print.sim==TRUE) {
            print(paste(strsplit(results.names, "[.]")[[j]][1],
                        strsplit(results.names, "[.]")[[j]][2], i))
          }


        }
      }


      if(strsplit(results.names, "[.]")[[j]][1] == "weighted") {

        for (i in 1:nsim) {

          # check weights are in correct range
          if(sum(obs.data[,weight.col] >= 0 & sum(obs.data[,weight.col] <= 1)) < nrow(obs.data)) {
            stop("weights must be in the range [0,1]")
          }

          # create subset of data
          samp.sites <- sample(1:nrow(obs.data),
                               size=ceiling(n.obs2),
                               replace=FALSE,
                               prob=(obs.data[, weight.col]))

          samp.data <- obs.data[samp.sites,]

          # fit model to subset of data
          samp.geodata <- geoR::as.geodata(samp.data, coords.col=coords.col, data.col=data.col)

          # get ini.cov.pars
          samp.variog <- geoR::variog(samp.geodata, messages=FALSE)
          samp.variofit <- geoR::variofit(samp.variog, cov.model="exponential", messages=FALSE)

          samp.likfit <- geoR::likfit(samp.geodata, cov.model="exponential", lik.method="REML",
                                      ini.cov.pars=samp.variofit,
                                      messages=FALSE)

          samp.site[,i] <- samp.sites

          cov.params[i,] <- c(samp.likfit$cov.pars, samp.likfit$nugget)

          cv <- geoR::xvalid(samp.geodata, model=samp.likfit,
                             message=FALSE)
          crossval.pred[,i] <- cv$predicted
          crossval.se[,i] <- cv$std.error
          cross.stats[i,1] <- sqrt((sum(cv$error^2))/length(cv))

          predictions <- geoR::krige.conv(samp.geodata, locations=preds.coords,
                                          krige=geoR::krige.control(obj.model=samp.likfit),
                                          output=geoR::output.control(messages=FALSE))
          preds.value[,i] <- predictions$predict
          preds.se[,i] <- sqrt(predictions$krige.var)

          if(print.sim==TRUE) {
            print(paste(strsplit(results.names, "[.]")[[j]][1],
                        strsplit(results.names, "[.]")[[j]][2], i))
          }



        }
      }


      results[[results.names[j]]] <- list(samp.site = samp.site,
                                          cov.params = cov.params,
                                          crossval.pred = crossval.pred,
                                          crossval.se = crossval.se,
                                          cross.stats = cross.stats,
                                          preds.value = preds.value,
                                          preds.se = preds.se)

    }
    results[["results.names"]] <- results.names
    results[["prop.sites"]] <- prop.sites
    results[["samp.scheme"]] <- samp.scheme
    results[["all"]] <- all
    results[["call"]] <- as.list(match.call())
    return(results)
  }

}
