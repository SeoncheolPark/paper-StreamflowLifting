#' Sampled networks from \code{\link{demoNet}}
#'
#' Results from applying \code{\link{sampNet}} to \code{\link{demoNet}} to
#' illustrate \code{\link{plot_sampNet}}.
#'
#' @format \code{\link{sampNet}} was applied to the data from
#'   \code{\link{demoNet}} using the following options:
#'
#' First, set the seed as shown below, and get the data:
#'
#' \code{set.seed(618999)} \cr
#' \code{obs.data <- getSSNdata.frame(demoNet, "Obs")} \cr
#' \code{preds.coords <- getSSNdata.frame(demoNet, "preds")}
#'
#' Next, run \code{sampDemo <- sampNet()} with the following arguments:
#'
#' \itemize{
#' \item{\code{x = }}{\code{obs.data}}
#' \item{\code{nsim = }}{\code{5}}
#' \item{\code{prop.sites = }}{\code{c(0.8, 0.7)}}
#' \item{\code{samp.scheme = }}{\code{c("random", "stratified", "weighted")}}
#' \item{\code{data.col = }}{\code{"Sim_Values"}}
#' \item{\code{coords.col = }}{\code{c(9, 10)}}
#' \item{\code{preds.coords = }}{\code{preds.coords[,c("NEAR_X", "NEAR_Y")]}}
#' \item{\code{siteID = }}{\code{"locID"}}
#' \item{\code{strat.col = }}{\code{"strata"}}
#' \item{\code{weight.col = }}{\code{"weight"}}
#' }
#'
#'  The \code{sampDemo} object contains results for random, weighted and
#'  stratified sampling schemes where 70\% and 80\% of monitoring sites were
#'  retained.  Five sampled networks were created for each \code{prop.sites} and
#'  \code{samp.scheme} combination.
#'
#' @source See vignette for Creating Demo Data
"sampDemo"
