#' Perform PCA to identify common spatiotemporal patterns, adjusting for spatial
#' and/or temporal autocorrelation where appropriate.
#'
#' PCA is performed on spatiotemporal data in T-mode or S-mode using singular
#' value decomposition and can be adjusted for spatial and/or temporal
#' autocorrelation if appropriate.  If spatial or temporal weights are used then
#' the returned loadings and principal components are backtransformed using
#' transformation described in paper (being written).
#'
#' \code{x} must always have columns = monitoring sites and rows = time points,
#' regardless of the analysis the user intends to perform.  If \code{pca.mode =
#' "Tmode"} is specified, the function will automatically transpose the data.
#'
#' \code{x} If column names in \code{spatial.wt} are not the same or in a
#' different order from \code{rownames(x)} then an error message is produced.
#'
#' \code{scale} If TRUE then this is equivalent to performing PCA on the
#' correlation matrix.  If FALSE then this is equivalent to performing PCA on
#' the covariance matrix.
#'
#' \code{pca.mode} PCA can be applied to spatiotemporal data in either S-mode or
#' T-mode.  In S-mode, the columns of \code{x} are monitoring sites while the
#' rows are time points.  In T-mode, the columns of \code{x} are time points
#' while the rows are monitoring sites.
#'
#' \code{spatial.wt} The column names of the spatial weights matrix must be of
#' the same form and in the same order as \code{rownames(x)}.  An error is
#' produced if these do no match.  The user should order \code{x} and set
#' \code{rownames(x)} to correspond to the ordering in \code{spatial.wt}.
#'
#' @importFrom expm sqrtm
#'
#' @param x data.frame.  This should be a complete data set with no missing
#'   values.  Columns MUST be monitoring sites while rows MUST be time points.
#'   See Details for further information.
#' @param center logical. TRUE (default) will subtract column means from each
#'   column of the \code{x}.
#' @param scale logical.  TRUE will divide values in \code{x} by the column
#'   standard deviation.  Default is FALSE.  See Details for further
#'   information.
#' @param pca.mode character string.  MUST be specified.  Options are "Smode" or
#'   "T-mode".  See Details for further information.
#' @param pca.wt character vector.  "unweighted" (default) performs standard
#'   PCA, "spatial" adjusts for spatial autocorrelation, "temporal" adjusts for
#'   temporal autocorrelation, "spatiotemporal" adjusts for both spatial and
#'   temporal autocorrelation.  See Details for further information.
#' @param spatial.wt character string. Specifies the name of the matrix
#'   containing spatial weights.  Only used if \code{pca.wt} contains
#'   "spatial" or "spatiotemporal".   Can be created by the user or using
#'   \code{\link{createWeightS}}.
#' @param temporal.wt character string. Specifies the name of the matrix
#'   containing temporal weights.  Only used if \code{pca.wt} contains
#'   "temporal" or "spatiotemporal".  Can be created by the user or using
#'   \code{\link{createWeightT}}.
#' @param scree logical. TRUE (default) means a scree plot will be produced for
#'   the first k principal components.
#' @param k integer.  The number of principal components to use when producing
#'   scree plot.  Only used if \code{scree = TRUE}.  Default is 10.
#'
#' @return The function returns a list containing the following: Lists called
#'   code{unweighted}, \code{spatial}, \code{temporal}, \code{spatiotemporal}
#'   are produced depending on what was included in \code{pca.wt}.  Each of
#'   these contains \code{scores} = principal components, \code{loads} =
#'   loadings, \code{var} = variance of each principal component,
#'   \code{var.prop} = proportion of total variance explained by each principal
#'   component, \code{cumvar} = cumulative proportion of total variance
#'   explained.
#'
#' \code{pca.mode} contains the option used for \code{pca.mode}.
#'
#' \code{pca.wt} contains the option used for \code{pca.wt}.
#'
#' \code{call} contains the function call details.
#'
#' \code{row.names} contains the row names of the data.frame on which PCA was
#' performed.  This will be row names of \code{x} if \code{pca.mode = "Smode"}
#' or column names of \code{x} if \code{pca.mode = "Tmode"}.
#'
#' \code{col.names} contains the column names of the data.frame on which PCA was
#' performed.  This will be column names of \code{x} if \code{pca.mode =
#' "Smode"} or row names of \code{x} if \code{pca.mode = "Tmode"}.
#'
#' \code{var.value} contains \code{var} from each of \code{unweighted},
#' \code{spatial}, \code{temporal} and \code{spatiotemporal} combined into a
#' single matrix.
#'
#' \code{mean.timeseries} contains the average time series accross all sites.
#' This is stored for use when \code{plots="meanPM"} in \code{\link{plot_stpca}}
#' and is calculated as \code{apply(x, 1, mean)}.
#'
#' A scree plot of the first \code{k} principal components is produced if
#' \code{scree = TRUE}, with the variance of each component on the y-axis.  The
#' text labels on the plot show the percentage variance explained by each
#' principal component.
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' library(stpca)
#'
#' ## get data.frame with no missing values
#' data(demoY)
#'
#' ## unweighted T-mode PCA
#' Tmode.pca.uw <- stpca(x=demoY, pca.mode="Tmode", pca.wt=c("unweighted"))
#'
#' ## spatial, temporal, and spatiotemporal weighted S-mode PCA
#'
#' # create spatial weights
#' ssndata <- system.file("demoSSN/demoNet.ssn", package="stpca")
#' weightS <- createWeightS(ssndata=ssndata, afvcol="addfunccol")
#'
#' # create temporal weights
#' weightT <- createWeightT(n=nrow(demoY), rho=0.5)
#'
#' # check that the order of monitoring sites in x is the same as the
#' # columns in spatial.wt.  They are the same if the following results
#' # in 1.
#' mean(colnames(demoY) == colnames(weightS))
#'
#' # weighted S-mode PCA
#' Smode.pca.all <- stpca(x = demoY, pca.mode = "Smode",
#'                        spatial.wt = weightS,
#'                        temporal.wt = weightT,
#'                        pca.wt = c("unweighted", "spatial",
#'                                   "temporal", "spatiotemporal"))
#'
#' @export
stpca <- function(x, center=TRUE, scale=FALSE,
                  pca.mode=NULL, pca.wt=c("unweighted"),
                  spatial.wt, temporal.wt,
                  scree=TRUE, k=10) {

  if(missing(pca.mode)) {stop("you must specify pca.mode")}
  if(sum(is.na(x)) > 1) {stop("x must be complete, with no missing values")}

  # read in data and scale as appropriate
  data <- x

  variances <- matrix(nrow=min(ncol(data), nrow(data)), ncol=length(pca.wt))
  variances.prop <- matrix(nrow=min(ncol(data), nrow(data)), ncol=length(pca.wt))
  counter <- 1

  results <- list()

  if(pca.mode=="Smode") {
    data.c <- scale(data, center=center, scale=scale)

    # create spatial weights
    if(sum(c("spatial", "spatiotemporal") %in% pca.wt) > 0) {
      sp.wt <- solve(expm::sqrtm(t(spatial.wt)))
      rownames(sp.wt) <- rownames(spatial.wt)
      colnames(sp.wt) <- colnames(spatial.wt)
    }

    # create temporal weights
    if(sum(c("temporal", "spatiotemporal") %in% pca.wt) > 0) {
      t.wt <- expm::sqrtm(solve(temporal.wt))
    }


    # unweighted PCA
    if("unweighted" %in% pca.wt) {

      uw.svd <- svd(data.c)
      uw.scores <- data.c%*%uw.svd$v
      uw.loads <- uw.svd$v
      uw.var <- uw.svd$d^2
      uw.var.prop <- uw.var/sum(uw.svd$d^2)
      uw.cumvar <- cumsum(uw.var/sum(uw.svd$d^2))

      variances[,counter] <- uw.var
      variances.prop[,counter] <- uw.var.prop
      counter <- counter+1

      uw.results <- list()
      uw.results[["scores"]] <- uw.scores
      uw.results[["loads"]] <- uw.loads
      uw.results[["var"]] <- uw.var
      uw.results[["var.prop"]] <- uw.var.prop
      uw.results[["cumvar"]] <- uw.cumvar

      results[["unweighted"]] <- uw.results
    }

    # spatially weighted PCA
    if("spatial" %in% pca.wt) {

      if(mean(colnames(sp.wt)==colnames(data.c))!=1) {stop("columns and rows in spatial.wt are not in the same order as columns in x")}

      sp.svd <- svd(data.c %*% sp.wt)
      sp.scores <- data.c%*%sp.wt%*%sp.svd$v
      sp.loads <- t(solve(sp.wt))%*%sp.svd$v
      sp.var <- sp.svd$d^2
      sp.var.prop <- sp.var/sum(sp.svd$d^2)
      sp.cumvar <- cumsum(sp.var/sum(sp.svd$d^2))

      variances[,counter] <- sp.var
      variances.prop[,counter] <- sp.var.prop
      counter <- counter+1

      sp.results <- list()
      sp.results[["scores"]] <- sp.scores
      sp.results[["loads"]] <- sp.loads
      sp.results[["var"]] <- sp.var
      sp.results[["var.prop"]] <- sp.var.prop
      sp.results[["cumvar"]] <- sp.cumvar

      results[["spatial"]] <- sp.results
    }

    # temporally weighted PCA
    if("temporal" %in% pca.wt) {

      t.svd <- svd(t.wt%*%data.c)
      t.scores <- data.c%*%t.svd$v
      t.loads <- t.svd$v
      t.var <- t.svd$d^2
      t.var.prop <- t.var/sum(t.svd$d^2)
      t.cumvar <- cumsum(t.var/sum(t.svd$d^2))

      variances[,counter] <- t.var
      variances.prop[,counter] <- t.var.prop
      counter <- counter+1

      t.results <- list()
      t.results[["scores"]] <- t.scores
      t.results[["loads"]] <- t.loads
      t.results[["var"]] <- t.var
      t.results[["var.prop"]] <- t.var.prop
      t.results[["cumvar"]] <- t.cumvar

      results[["temporal"]] <- t.results
    }

    # spatiotemporally weighted PCA
    if("spatiotemporal" %in% pca.wt) {

      if(mean(colnames(sp.wt)==colnames(data.c))!=1) {stop("columns and rows in spatial.wt are not in the same order as columns in x")}

      spt.svd <- svd(t.wt%*%data.c%*%sp.wt)
      spt.scores <- data.c%*%sp.wt%*%spt.svd$v
      spt.loads <- t(solve(sp.wt))%*%spt.svd$v
      spt.var <- spt.svd$d^2
      spt.var.prop <- spt.var/sum(spt.svd$d^2)
      spt.cumvar <- cumsum(spt.var/sum(spt.svd$d^2))

      variances[,counter] <- spt.var
      variances.prop[,counter] <- spt.var.prop
      counter <- counter+1

      spt.results <- list()
      spt.results[["scores"]] <- spt.scores
      spt.results[["loads"]] <- spt.loads
      spt.results[["var"]] <- spt.var
      spt.results[["var.prop"]] <- spt.var.prop
      spt.results[["cumvar"]] <- spt.cumvar

      results[["spatiotemporal"]] <- spt.results
    }
  }


  if(pca.mode=="Tmode") {
    data <- t(data)
    data.c <- scale(data, center=center, scale=scale)

    # create spatial weights
    if(sum(c("spatial", "spatiotemporal") %in% pca.wt) > 0) {
      sp.wt <- expm::sqrtm(solve(spatial.wt))
      rownames(sp.wt) <- rownames(spatial.wt)
      colnames(sp.wt) <- colnames(spatial.wt)
    }

    # create temporal weights
    if(sum(c("temporal", "spatiotemporal") %in% pca.wt) > 0) {
      t.wt <- solve(expm::sqrtm(t(temporal.wt)))
    }

    # unweighted PCA
    if("unweighted" %in% pca.wt) {

      uw.svd <- svd(data.c)
      uw.scores <- data.c%*%uw.svd$v
      uw.loads <- uw.svd$v
      uw.var <- uw.svd$d^2
      uw.var.prop <- uw.var/sum(uw.svd$d^2)
      uw.cumvar <- cumsum(uw.var/sum(uw.svd$d^2))

      variances[,counter] <- uw.var
      variances.prop[,counter] <- uw.var.prop
      counter <- counter+1

      uw.results <- list()
      uw.results[["scores"]] <- uw.scores
      uw.results[["loads"]] <- uw.loads
      uw.results[["var"]] <- uw.var
      uw.results[["var.prop"]] <- uw.var.prop
      uw.results[["cumvar"]] <- uw.cumvar

      results[["unweighted"]] <- uw.results
    }
    # spatially weighted PCA
    if("spatial" %in% pca.wt) {

      if(mean(colnames(sp.wt)==rownames(data.c))!=1) {stop("columns and rows in spatial.wt are not in the same order as rows in x")}

      sp.svd <- svd(sp.wt%*%data.c)
      sp.scores <- data.c%*%sp.svd$v
      sp.loads <- sp.svd$v
      sp.var <- sp.svd$d^2
      sp.var.prop <- sp.var/sum(sp.svd$d^2)
      sp.cumvar <- cumsum(sp.var/sum(sp.svd$d^2))

      variances[,counter] <- sp.var
      variances.prop[,counter] <- sp.var.prop
      counter <- counter+1

      sp.results <- list()
      sp.results[["scores"]] <- sp.scores
      sp.results[["loads"]] <- sp.loads
      sp.results[["var"]] <- sp.var
      sp.results[["var.prop"]] <- sp.var.prop
      sp.results[["cumvar"]] <- sp.cumvar

      results[["spatial"]] <- sp.results
    }

    # temporally weighted PCA
    if("temporal" %in% pca.wt) {

      t.svd <- svd(data.c %*% t.wt)
      t.scores <- data.c%*%t.wt%*%t.svd$v
      t.loads <- t(solve(t.wt))%*%t.svd$v
      t.var <- t.svd$d^2
      t.var.prop <- t.var/sum(t.svd$d^2)
      t.cumvar <- cumsum(t.var/sum(t.svd$d^2))

      variances[,counter] <- t.var
      variances.prop[,counter] <- t.var.prop
      counter <- counter+1

      t.results <- list()
      t.results[["scores"]] <- t.scores
      t.results[["loads"]] <- t.loads
      t.results[["var"]] <- t.var
      t.results[["var.prop"]] <- t.var.prop
      t.results[["cumvar"]] <- t.cumvar

      results[["temporal"]] <- t.results
    }

    # spatiotemporally weighted PCA
    if("spatiotemporal" %in% pca.wt) {

      if(mean(colnames(sp.wt)==rownames(data.c))!=1) {stop("columns and rows in spatial.wt are not in the same order as columns in x")}

      spt.svd <- svd(sp.wt%*%data.c%*%t.wt)
      spt.scores <- data.c%*%t.wt%*%spt.svd$v
      spt.loads <- t(solve(t.wt))%*%spt.svd$v
      spt.var <- spt.svd$d^2
      spt.var.prop <- spt.var/sum(spt.svd$d^2)
      spt.cumvar <- cumsum(spt.var/sum(spt.svd$d^2))

      variances[,counter] <- spt.var
      variances.prop[,counter] <- spt.var.prop
      counter <- counter+1

      spt.results <- list()
      spt.results[["scores"]] <- spt.scores
      spt.results[["loads"]] <- spt.loads
      spt.results[["var"]] <- spt.var
      spt.results[["var.prop"]] <- spt.var.prop
      spt.results[["cumvar"]] <- spt.cumvar

      results[["spatiotemporal"]] <- spt.results
    }
  }

  names(variances) <- pca.wt

  if(scree==TRUE) {
    nplot <- length(pca.wt)
    if(nplot==4) {
      par(mfrow=c(2,2))
    }
    else {
      par(mfrow=c(1,nplot))
    }

    pca.type <- pca.wt
    ymin = min(variances[k,])
    ymax = 1.05*max(variances[1,])

    for(i in 1:nplot) {
      plot(1:k, variances[1:k, i], type="b",
           main=paste(pca.mode, pca.type[i], "PCA"),
           xlab="k", ylab="variance",
           ylim=c(ymin, ymax))
      text(1:k, variances[1:k, i], labels=paste(100*round(variances.prop[1:k, i],2), "%"),
           pos=3)
    }
  }

  mean.timeseries <- apply(x, 1, mean)

  results[["pca.mode"]] <- pca.mode
  results[["pca.wt"]] <- pca.wt
  results[["call"]] <- match.call(stpca)
  results[["row.names"]] <- rownames(data)
  results[["col.names"]] <- colnames(data)
  results[["var.value"]] <- variances
  results[["mean.timeseries"]] <- mean.timeseries


  return(results)

}



