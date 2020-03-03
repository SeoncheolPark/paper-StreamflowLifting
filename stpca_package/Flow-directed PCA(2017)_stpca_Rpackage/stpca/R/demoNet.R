#' Simulated observations of a single water quality parameter on a river
#' network.
#'
#' Simulated values of a single water quality parameter recorded at 30
#' monitoring sites on a simulated river network.
#'
#' @format An object of \code{\link[SSN]{SpatialStreamNetwork-class}} containing
#'   simulated observations at a single time point of a single water quality
#'   parameter.  The SSN object contains two data.frames.  "Obs" contains data
#'   and other variables related to 30 observed locations and "preds" contains
#'   the same information but for 30 unobserved prediction locations.  Full
#'   details of the structure of an SSN object can be found by looking at
#'   \code{\link[SSN]{SpatialStreamNetwork-class}}.  The "Obs" and "preds"
#'   data.frames for this simulated network contain the following variables:
#'
#'   \describe{
#'   \item{locID}{unique identifier for each monitoring site.}
#'   \item{upDist}{the upstream distance (unitless) from outlet to the monitoring site.}
#'   \item{pid}{unique identifier for each row of the "Obs" and "preds" data.frames.}
#'   \item{netID}{identifier for the river network.}
#'   \item{rid}{identifier for the stream segment.}
#'   \item{ratio}{location of monitoring site on a stream segment as a proportion of the segment length from the downstream end of the segment.}
#'   \item{shreve}{Shreve's number for the stream segment.}
#'   \item{addfunccol}{additive function value, based on Shreve's number.}
#'   \item{NEAR_X}{x-coordinate of monitoring site location (equivalent to Easting in real data).}
#'   \item{NEAR_Y}{y-coordinate of monitoring site location (equivalent to Northing in real data).}
#'   \item{Sim_Values}{simulated observed values of a single water quality parameter.}
#'   \item{strata}{a randomly allocated category used for stratified sampling.}
#'   \item{weight}{a value in the range (0,1) used for weighted sampling.}
#'   }
#'
#' @source See vignette for Creating Demo Data
"demoNet"
