#' Create an asymmetric matrix of spatial weights from an \code{\link[SSN]{SSN}}
#' object.
#'
#' The \code{createWeightS} function extracts the stream distance matrix for
#' observed data from a \code{\link[SSN]{SpatialStreamNetwork-class}} object and
#' constructs an asymmetric weight matrix based on the specified additive
#' function column.
#'
#' Peterson and ver Hoef (2010) Appendix A shows how weights can be calculated
#' for the tail-up model, where weights include information about the flow
#' connected structure in the river network, network distance and relative
#' influence of upstream sites on downstream sites.  Weights are based on the
#' variable used to calculate the additive function in the SSN object.
#' \code{createWeightS} follows the steps described in Appendix A but stops
#' before the matrix is forced to symmetry.
#'
#' @import SSN
#'
#' @importFrom Matrix bdiag
#'
#' @importFrom matrixcalc hadamard.prod
#'
#' @param ssndata filepath for an object of
#'   \code{\link[SSN]{SpatialStreamNetwork-class}}.  Include the .ssn folder in
#'   the path.
#' @param afvcol character string. This should be the column name in the SSN
#'   object containing additive function values.
#'
#' @return Produces an asymmetric \eqn{p \times p} matrix (p = number of
#'   monitoring sites) where columns = upstream sites, rows = downstream sites.
#'
#' @references Peterson, E. E. and J. M. ver Hoef (2010). A mixed-model moving
#'   average approach to geostatistical modeling in stream networks. Ecology 91,
#'   644-651.
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' library(stpca)
#'
#' ## get filepath for SSN object
#' ssndata <- system.file("demoSSN/demoNet.ssn", package = "stpca")
#'
#' ## create matrix of spatial weights
#' weightS <- createWeightS(ssndata = ssndata, afvcol = "addfunccol")
#'
#' @export
createWeightS <- function(ssndata, afvcol) {

  x <- SSN::importSSN(ssndata)
  if(class(x) != "SpatialStreamNetwork") {stop("x is not an SSN object")}

  path = paste0(x@path, "/distance/", "Obs")
  distlist = list.files(path)
  distmats = vector("list", length(distlist))
  pid.names <- list()
  for (i in 1:length(distlist)) {
    path1 = paste0(path, "/", distlist[i])
    file_handle <- file(path1, open = "rb")
    distmat <- unserialize(file_handle)
    close(file_handle)

    distmat1 <- matrix(0, nrow=nrow(distmat), ncol=ncol(distmat))
    for (j in 1:ncol(distmat)) {
      for (k in 1:nrow(distmat)) {
        distmat1[k,j] <- ifelse(distmat[k,j] != 0 & distmat[j,k] == 0, 1, 0)
      }
    }
    distmats[[i]] = distmat1
    pid.names[[i]] <- colnames(distmat)
  }


  # put all networks together in one big matrix
  connectedness <- as.matrix(Matrix::bdiag(distmats))
  colnames(connectedness) <- unlist(pid.names)
  rownames(connectedness) <- unlist(pid.names)

  # get data frame containing additive function values
  afvdata <- SSN::getSSNdata.frame(x)

  # order afvdata by pid
  afvdata <- afvdata[order(afvdata$pid),]

  # order connectedness by pid.names so in same order as ssndata
  ordrow <- order(as.numeric(rownames(connectedness)))
  ordcol <- order(as.numeric(colnames(connectedness)))
  connectedness <- connectedness[ordrow, ordcol, drop = F]

  # check colnames in connectedness is same order as in afvdata
  if(mean(afvdata$pid==colnames(connectedness)) != 1) {stop("data frame observations are not in same order as distance matrix")}

  # create weights based on afv values
  x <- afvdata[,afvcol] #i=upstream
  y <- afvdata[,afvcol] #j=downstream
  weight <- matrix(0, nrow=length(x), ncol=length(x))
  for(i in 1:length(x)) {
    weight[i,] <- sqrt(x/y[i])
  }

  weight.connected <- matrixcalc::hadamard.prod(connectedness, weight)

  # make the diagonal = 1 as in Peterson (2010) appendix A
  # colums are "from" and rows are "to"
  diag(weight.connected) <- rep(1, length(diag(weight.connected)))
  if(sum(range(weight.connected)%in%c(0,1)) != 2) {stop("weights are not in range [0,1]")}

  # put weight.connected rows/columns into pid order in SSN object
  x2 <- SSN::importSSN(ssndata)
  order.data <- SSN::getSSNdata.frame(x2)
  ordrow2 <- as.character(order.data$pid)
  ordcol2 <- as.character(order.data$pid)
  weight.connected <- weight.connected[ordrow2, ordcol2, drop = F]


  return(weight.connected)
}



