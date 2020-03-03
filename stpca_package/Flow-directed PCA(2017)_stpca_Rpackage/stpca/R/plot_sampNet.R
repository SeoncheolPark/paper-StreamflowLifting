#' Plot the results from \code{\link{sampNet}}
#'
#' This function will plot the results from creating several sampled networks
#' using \code{\link{sampNet}}.  This will allow the user to make comparisons
#' between sampling schemes and proportion of monitoring sites retained in
#' sampled networks.
#'
#' This function can be used to produce plots of several parameters/values
#' estimated for each \code{prop.sites} and \code{samp.scheme} combination.
#'
#' If \code{cov.params.plots = TRUE} and \code{class(x)} is
#' \code{\link[SSN]{SpatialStreamNetwork-class}}, then five plots will be
#' produced showing the covariance function parameters estimated for \code{nsim}
#' sampled networks: Tail-up range, Euclidean range, Tail-up partial sill,
#' Euclidean partial sill, and nugget.  If \code{class(x)} is
#' \code{"data.frame"} then thre plots will be produced: range, partial sill,
#' and nugget.
#'
#' If \code{pred.plots = TRUE}, then three plots will be produced: root mean
#' square predicted error (RMSPE) for observed locations, prediction error for
#' prediction locations, and average kriging standard error (AKSE) ratio.
#'
#' RMSPE is calculated for observed locations as follows:
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
#' Prediction error is calculated for prediction locations as follows:
#'
#' \deqn{\sqrt{\frac{1}{N_{pred}}\sum^{N_{pred}}_{i=1}\left(\hat{Y}_{full}(u_i)-\hat{Y}_{sampled}(u_i)\right)^2}}
#'
#' where \eqn{N_{pred}} is the number of prediction locations,
#' \eqn{\hat{Y}_{full}(u_i)} is the value predicted at prediction (unobserved)
#' location \eqn{i} where the prediction is made using a model fitted to the
#' full monitoring network, and \eqn{\hat{Y}_{sampled}(u_i)} is the value
#' predicted at prediction (unobserved) location \eqn{i} where the prediction is
#' made using a model fitted to a sampled (reduced) network.  This is a measure
#' of how well the model based on fewer monitoring sites predicts at unobserved
#' locations compared to the model based on the full monitoring network.  As
#' with RMSPE, low values are desired.
#'
#' AKSE ratio can be calculated as the ratio of AKSE from a model based on a
#' sampled (reduced) network to AKSE from a model based on the full monitoring
#' network: \eqn{\mathrm{AKSE}_{sampled}/\mathrm{AKSE}_{full}} where
#'
#' \deqn{\mathrm{AKSE} = \sqrt{\frac{1}{N_{pred}}\sum^{N_{pred}}_{i=1}\sigma^2(u_i)}}
#'
#' and \eqn{N_{pred}} is the number of prediction locations, and
#' \eqn{\sigma^2(u_i)} is the kriging variance at prediction location \eqn{i}.
#'
#' If \code{log.TUrange = TRUE} and \code{class(x)} is
#' \code{\link[SSN]{SpatialStreamNetwork-class}}, then values for the Tail-up
#' range parameter will have a natural log transformation applied.  This is
#' useful if the \code{nsim} estimated values vary by orders of magnitude and it
#' is left to the user to decide which display (transformed values or
#' untransformed values) is most useful for their own application.  The
#' transformation is also available for the Euclidean range parameter for
#' \code{class(x)} is \code{\link[SSN]{SpatialStreamNetwork-class}} or
#' \code{"data.frame"}.
#'
#' @importFrom graphics segments
#' @import SSN
#'
#' @param samp.results	character.  Contains the name of the object in which the
#'   results from \code{\link{sampNet}} are stored.
#' @param cov.params.plots	logical.  If \code{TRUE}, then plots of the
#'   covariance parameters estimated for \code{nsim} sampled networks will be
#'   produced.  See Details for further information.
#' @param pred.plots	logical.  If \code{TRUE}, then plots of predicted values and
#' standard errors estimated for \code{nsim} sampled networks will be produced.  See
#' Details for further information.
#' @param log.TUrange	logical.  If \code{TRUE}, then the y-axis and values
#'   plotted in the Tail-up range parameter plot will be transformed using a
#'   natural logarithmic transformation.  See Details.
#' @param log.Eucrange	logical.  If \code{TRUE}, then the y-axis and values
#'   plotted in the Euclidean range parameter plot will be transformed using a
#'   natural logarithmic transformation.  See Details.
#'
#' @return Two sets of plots will be produced depending on argument values
#'   specified.
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' data(sampDemo)
#' plot_sampNet(sampDemo,
#'              cov.params.plots = TRUE,
#'              pred.plots = TRUE)
#'
#' @export
plot_sampNet <- function(samp.results,
                         cov.params.plots=TRUE, pred.plots=TRUE,
                         log.TUrange=FALSE, log.Eucrange=FALSE) {

  n.results <- length(samp.results$results.names)
  n.prop.sites <- length(samp.results$prop.sites)
  n.samp.scheme <- length(samp.results$samp.scheme)

  # get some values for plots
  if(length(samp.results$samp.scheme) == 1) {
    xval <- 0.5
  }

  if(length(samp.results$samp.scheme) == 2) {
    xval <- c(0.33, 0.66)
  }

  if(length(samp.results$samp.scheme) == 3) {
    xval <- c(0.25, 0.5, 0.75)
  }

  if(samp.results$all$class == "SpatialStreamNetwork") {

    if(cov.params.plots == TRUE) {

      # TU range
      TUrange <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(TUrange) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        TUrange[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        TUrange[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        TUrange[i,3:5] <- quantile(subdata$cov.params[,"TUrange"], probs=c(0.25, 0.5, 0.75))
      }
      TUrange[,6] <- as.numeric(as.character(factor(unique(TUrange$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      TUrange[,7] <- as.numeric(as.character(factor(TUrange$samp.scheme, labels=xval)))
      TUrange[,8] <- TUrange$scheme.val/min(TUrange$scheme.val)

      ylab = "Tail-Up range"
      ymax <- max(TUrange$q75, samp.results$all$all.covparams[2,3])
      ymin <- min(TUrange$q25, samp.results$all$all.covparams[2,3])

      if(log.TUrange == TRUE) {
        TUrange[,3:5] <- log(TUrange[,3:5])
        ylab = "ln(Tail-Up range)"
        ymax <- max(TUrange$q75, log(samp.results$all$all.covparams[2,3]))
        ymin <- min(TUrange$q25, log(samp.results$all$all.covparams[2,3]))
      }

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"TUrange"]),
               n.prop.sites+2),
           type="n",
           ylim=c(ymin, ymax),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Tail-up range",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0=TUrange$prop.val+TUrange$scheme.val
      graphics::segments(x0=x0,
                         y0=TUrange$q25,
                         x1=x0,
                         y1=TUrange$q75,
                         lty=TUrange$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(TUrange$prop.val+TUrange$scheme.val, TUrange$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      if(log.TUrange == TRUE) {
        abline(h=log(samp.results$all$all.covparams[2,3]), col="red")
      } else {
        abline(h=samp.results$all$all.covparams[2,3], col="red")
      }

      # add a legend
      legend("topleft", legend=c(unique(TUrange$samp.scheme), "all"),
             lty=c(unique(TUrange$line.val),1), lwd=2,
             col=c(rep(1, length(unique(TUrange$line.val))),2),
             bty="n")


      # Euc range
      Eucrange <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(Eucrange) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        Eucrange[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        Eucrange[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        Eucrange[i,3:5] <- quantile(subdata$cov.params[,"EucRange"], probs=c(0.25, 0.5, 0.75))
      }
      Eucrange[,6] <- as.numeric(as.character(factor(unique(Eucrange$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      Eucrange[,7] <- as.numeric(as.character(factor(Eucrange$samp.scheme, labels=xval)))
      Eucrange[,8] <- Eucrange$scheme.val/min(Eucrange$scheme.val)

      ylab = "Euclidean range"
      ymax <- max(Eucrange$q75, samp.results$all$all.covparams[4,3])
      ymin <- min(Eucrange$q25, samp.results$all$all.covparams[4,3])

      if(log.Eucrange == TRUE) {
        Eucrange[,3:5] <- log(Eucrange[,3:5])
        ylab = "ln(Euclidean range)"
        ymax <- max(Eucrange$q75, log(samp.results$all$all.covparams[4,3]))
        ymin <- min(Eucrange$q25, log(samp.results$all$all.covparams[4,3]))
      }

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"EucRange"]),
               n.prop.sites+2),
           type="n",
           ylim=c(ymin, ymax),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Euclidean range",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- Eucrange$prop.val+Eucrange$scheme.val
      graphics::segments(x0=x0,
                         y0=Eucrange$q25,
                         x1=x0,
                         y1=Eucrange$q75,
                         lty=Eucrange$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(Eucrange$prop.val+Eucrange$scheme.val, Eucrange$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      if(log.Eucrange == TRUE) {
        abline(h=log(samp.results$all$all.covparams[4,3]), col="red")
      } else {
        abline(h=samp.results$all$all.covparams[4,3], col="red")
      }

      # add a legend
      legend("topleft", legend=c(unique(Eucrange$samp.scheme), "all"),
             lty=c(unique(Eucrange$line.val),1), lwd=2,
             col=c(rep(1, length(unique(Eucrange$line.val))),2),
             bty="n")


      # TU partial sill
      TUparsill <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(TUparsill) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        TUparsill[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        TUparsill[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        TUparsill[i,3:5] <- quantile(subdata$cov.params[,"TUparsill"], probs=c(0.25, 0.5, 0.75))
      }
      TUparsill[,6] <- as.numeric(as.character(factor(unique(TUparsill$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      TUparsill[,7] <- as.numeric(as.character(factor(TUparsill$samp.scheme, labels=xval)))
      TUparsill[,8] <- TUparsill$scheme.val/min(TUparsill$scheme.val)

      ylab = "Tail-Up partial sill"

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"TUparsill"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(TUparsill$q25, samp.results$all$all.covparams[1,3]),
                  max(TUparsill$q75, samp.results$all$all.covparams[1,3])),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Tail-up partial sill",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- TUparsill$prop.val+TUparsill$scheme.val
      graphics::segments(x0=x0,
                         y0=TUparsill$q25,
                         x1=x0,
                         y1=TUparsill$q75,
                         lty=TUparsill$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(TUparsill$prop.val+TUparsill$scheme.val, TUparsill$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=samp.results$all$all.covparams[1,3], col="red")
      # add a legend
      legend("topleft", legend=c(unique(TUparsill$samp.scheme), "all"),
             lty=c(unique(TUparsill$line.val),1), lwd=2,
             col=c(rep(1, length(unique(TUparsill$line.val))),2),
             bty="n")


      # Euclidean partial sill
      EucParsill <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(EucParsill) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        EucParsill[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        EucParsill[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        EucParsill[i,3:5] <- quantile(subdata$cov.params[,"EucParsill"], probs=c(0.25, 0.5, 0.75))
      }
      EucParsill[,6] <- as.numeric(as.character(factor(unique(EucParsill$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      EucParsill[,7] <- as.numeric(as.character(factor(EucParsill$samp.scheme, labels=xval)))
      EucParsill[,8] <- EucParsill$scheme.val/min(EucParsill$scheme.val)

      ylab = "Euclidean partial sill"

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"EucParsill"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(EucParsill$q25, samp.results$all$all.covparams[3,3]),
                  max(EucParsill$q75, samp.results$all$all.covparams[3,3])),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Euclidean partial sill",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- EucParsill$prop.val+EucParsill$scheme.val
      graphics::segments(x0=x0,
                         y0=EucParsill$q25,
                         x1=x0,
                         y1=EucParsill$q75,
                         lty=EucParsill$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(EucParsill$prop.val+EucParsill$scheme.val, EucParsill$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=samp.results$all$all.covparams[3,3], col="red")
      # add a legend
      legend("topleft", legend=c(unique(EucParsill$samp.scheme), "all"),
             lty=c(unique(EucParsill$line.val),1), lwd=2,
             col=c(rep(1, length(unique(EucParsill$line.val))),2),
             bty="n")


      # Nugget
      Nugget <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(Nugget) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        Nugget[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        Nugget[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        Nugget[i,3:5] <- quantile(subdata$cov.params[,"Nugget"], probs=c(0.25, 0.5, 0.75))
      }
      Nugget[,6] <- as.numeric(as.character(factor(unique(Nugget$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      Nugget[,7] <- as.numeric(as.character(factor(Nugget$samp.scheme, labels=xval)))
      Nugget[,8] <- Nugget$scheme.val/min(Nugget$scheme.val)

      ylab = "Nugget"

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"Nugget"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(Nugget$q25, samp.results$all$all.covparams[5,3]),
                  max(Nugget$q75, samp.results$all$all.covparams[5,3])),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Nugget",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- Nugget$prop.val+Nugget$scheme.val
      graphics::segments(x0=x0,
                         y0=Nugget$q25,
                         x1=x0,
                         y1=Nugget$q75,
                         lty=Nugget$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(Nugget$prop.val+Nugget$scheme.val, Nugget$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=samp.results$all$all.covparams[5,3], col="red")
      # add a legend
      legend("topleft", legend=c(unique(Nugget$samp.scheme), "all"),
             lty=c(unique(Nugget$line.val),1), lwd=2,
             col=c(rep(1, length(unique(Nugget$line.val))),2),
             bty="n")
    }



    if(pred.plots == TRUE) {

      # RMSPE from LOOCV
      RMSPE <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(RMSPE) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        RMSPE[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        RMSPE[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        RMSPE[i,3:5] <- quantile(subdata$cross.stats[,"RMSPE"], probs=c(0.25, 0.5, 0.75))
      }
      RMSPE[,6] <- as.numeric(as.character(factor(unique(RMSPE$prop.sites), labels=c(length(unique(RMSPE$prop.sites)):1))))
      RMSPE[,7] <- as.numeric(as.character(factor(RMSPE$samp.scheme, labels=xval)))
      RMSPE[,8] <- RMSPE$scheme.val/min(RMSPE$scheme.val)

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cross.stats[,"RMSPE"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(RMSPE$q25, samp.results$all$all.cross.stats$RMSPE),
                  max(RMSPE$q75, samp.results$all$all.cross.stats$RMSPE)),
           xlab="Percentage of sites retained", ylab=expression("RMSPE"[LOOCV]),
           main="RMSPE: Observed sites",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- RMSPE$prop.val+RMSPE$scheme.val
      graphics::segments(x0=x0,
                         y0=RMSPE$q25,
                         x1=x0,
                         y1=RMSPE$q75,
                         lty=RMSPE$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(RMSPE$prop.val+RMSPE$scheme.val, RMSPE$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=samp.results$all$all.cross.stats$RMSPE, col="red")
      # add a legend
      legend("topleft", legend=c(unique(RMSPE$samp.scheme), "all"),
             lty=c(unique(RMSPE$line.val),1), lwd=2,
             col=c(rep(1, length(unique(RMSPE$line.val))),2),
             bty="n")


      # AKSE and prediction error

      # get predicted values
      preds2 <- SSN::getSSNdata.frame(samp.results$all$all.preds$ssn.object, samp.results$call$predpts)
      N <- dim(preds2)[1]
      predsSEcol <- paste(samp.results$call$response.col, "predSE", sep=".")
      AKSE_all <- sqrt((sum(preds2[,predsSEcol]^2))/N)

      # Prediction error with prediction from full network as "true" and prediction
      # from sampled networks as "predicted"
      pred.error <- list()

      for(i in 1:length(samp.results$results.names)) {
        error.mat <- matrix(nrow=dim(preds2)[1],
                            ncol=samp.results$call$nsim)
        pred.data <- samp.results[[samp.results$results.names[i]]]
        for(j in 1:ncol(error.mat)) {
          error.mat[,j] <- preds2[, predsSEcol] - pred.data$preds.value[,j]
        }

        pred.error[[i]] <- sqrt((apply(error.mat^2, 2, sum)/nrow(preds2)))

      }
      names(pred.error) <- samp.results$results.names

      # plot prediction error
      pred.plot <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(pred.plot) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        pred.plot[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        pred.plot[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        pred.plot[i,3:5] <- quantile(pred.error[[i]], probs=c(0.25, 0.5, 0.75))
      }
      pred.plot[,6] <- as.numeric(as.character(factor(unique(pred.plot$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      pred.plot[,7] <- as.numeric(as.character(factor(pred.plot$samp.scheme, labels=xval)))
      pred.plot[,8] <- pred.plot$scheme.val/min(pred.plot$scheme.val)

      plot(0:(n.prop.sites+1),
           rep(1, n.prop.sites+2),
           type="n",
           ylim=c(min(pred.plot$q25), max(pred.plot$q75)),
           xlab="Percentage of sites retained", ylab="prediction error",
           main="Prediction error at prediction sites",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- pred.plot$prop.val+pred.plot$scheme.val
      graphics::segments(x0=x0,
                         y0=pred.plot$q25,
                         x1=x0,
                         y1=pred.plot$q75,
                         lty=pred.plot$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(pred.plot$prop.val+pred.plot$scheme.val, pred.plot$q50, pch=19, col="blue", cex=1.5)
      # add a legend
      legend("topleft", legend=c(unique(pred.plot$samp.scheme)),
             lty=c(unique(pred.plot$line.val)), lwd=2,
             col=c(rep(1, length(unique(pred.plot$line.val)))),
             bty="n")


      # AKSE: average kriging standard error ratio
      akse <- list()

      for(i in 1:length(samp.results$results.names)) {
        akse.data <- samp.results[[samp.results$results.names[i]]]
        akse[[i]] <- (sqrt(apply(akse.data$preds.se^2, 2, sum)/N))/AKSE_all
      }
      names(akse) <- samp.results$results.names

      # plot akse ratio
      akse.plot <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(akse.plot) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        akse.plot[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        akse.plot[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        akse.plot[i,3:5] <- quantile(akse[[i]], probs=c(0.25, 0.5, 0.75))
      }
      akse.plot[,6] <- as.numeric(as.character(factor(unique(akse.plot$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      akse.plot[,7] <- as.numeric(as.character(factor(akse.plot$samp.scheme, labels=xval)))
      akse.plot[,8] <- akse.plot$scheme.val/min(akse.plot$scheme.val)

      plot(0:(n.prop.sites+1),
           rep(1, n.prop.sites+2),
           type="n",
           ylim=c(min(akse.plot$q25), max(akse.plot$q75)),
           xlab="Percentage of sites retained", ylab="ratio",
           main=expression("AKSE"[subset]/"AKSE"[all]),
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- akse.plot$prop.val+akse.plot$scheme.val
      graphics::segments(x0=x0,
                         y0=akse.plot$q25,
                         x1=x0,
                         y1=akse.plot$q75,
                         lty=akse.plot$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(akse.plot$prop.val+akse.plot$scheme.val, akse.plot$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=1, col="red")
      # add a legend
      legend("topleft", legend=c(unique(akse.plot$samp.scheme)),
             lty=c(unique(akse.plot$line.val)), lwd=2,
             col=c(rep(1, length(unique(akse.plot$line.val)))),
             bty="n")

    }

  }

  if(samp.results$all$class == "data.frame") {

    if(cov.params.plots == TRUE) {

      # range
      range <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(range) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        range[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        range[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        range[i,3:5] <- quantile(subdata$cov.params[,"range"], probs=c(0.25, 0.5, 0.75))
      }
      range[,6] <- as.numeric(as.character(factor(unique(range$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      range[,7] <- as.numeric(as.character(factor(range$samp.scheme, labels=xval)))
      range[,8] <- range$scheme.val/min(range$scheme.val)

      ylab = "Range"
      ymax <- max(range$q75, samp.results$all$all.covparams[1,2])
      ymin <- min(range$q25, samp.results$all$all.covparams[1,2])

      if(log.Eucrange == TRUE) {
        range[,3:5] <- log(range[,3:5])
        ylab = "ln(Range)"
        ymax <- max(range$q75, log(samp.results$all$all.covparams[1,2]))
        ymin <- min(range$q25, log(samp.results$all$all.covparams[1,2]))
      }

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"range"]),
               n.prop.sites+2),
           type="n",
           ylim=c(ymin, ymax),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Range",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0=range$prop.val+range$scheme.val
      graphics::segments(x0=x0,
                         y0=range$q25,
                         x1=x0,
                         y1=range$q75,
                         lty=range$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(range$prop.val+range$scheme.val, range$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      if(log.Eucrange == TRUE) {
        abline(h=log(samp.results$all$all.covparams[1,2]), col="red")
      } else {
        abline(h=samp.results$all$all.covparams[1,2], col="red")
      }

      # add a legend
      legend("topleft", legend=c(unique(range$samp.scheme), "all"),
             lty=c(unique(range$line.val),1), lwd=2,
             col=c(rep(1, length(unique(range$line.val))),2),
             bty="n")



      # partial sill
      parsill <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(parsill) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        parsill[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        parsill[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        parsill[i,3:5] <- quantile(subdata$cov.params[,"parsill"], probs=c(0.25, 0.5, 0.75))
      }
      parsill[,6] <- as.numeric(as.character(factor(unique(parsill$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      parsill[,7] <- as.numeric(as.character(factor(parsill$samp.scheme, labels=xval)))
      parsill[,8] <- parsill$scheme.val/min(parsill$scheme.val)

      ylab = "Partial sill"

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"parsill"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(parsill$q25, samp.results$all$all.covparams[1,1]),
                  max(parsill$q75, samp.results$all$all.covparams[1,1])),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: Partial sill",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- parsill$prop.val+parsill$scheme.val
      graphics::segments(x0=x0,
                         y0=parsill$q25,
                         x1=x0,
                         y1=parsill$q75,
                         lty=parsill$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(parsill$prop.val+parsill$scheme.val, parsill$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=samp.results$all$all.covparams[1,1], col="red")
      # add a legend
      legend("topleft", legend=c(unique(parsill$samp.scheme), "all"),
             lty=c(unique(parsill$line.val),1), lwd=2,
             col=c(rep(1, length(unique(parsill$line.val))),2),
             bty="n")



      # Nugget
      nugget <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(nugget) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        nugget[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        nugget[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        nugget[i,3:5] <- quantile(subdata$cov.params[,"Nugget"], probs=c(0.25, 0.5, 0.75))
      }
      nugget[,6] <- as.numeric(as.character(factor(unique(nugget$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      nugget[,7] <- as.numeric(as.character(factor(nugget$samp.scheme, labels=xval)))
      nugget[,8] <- nugget$scheme.val/min(nugget$scheme.val)

      ylab = "Nugget"

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cov.params[,"Nugget"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(nugget$q25, samp.results$all$all.covparams[1,3]),
                  max(nugget$q75, samp.results$all$all.covparams[1,3])),
           xlab="Percentage of sites retained", ylab=ylab,
           main="Covariance parameters: nugget",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- nugget$prop.val+nugget$scheme.val
      graphics::segments(x0=x0,
                         y0=nugget$q25,
                         x1=x0,
                         y1=nugget$q75,
                         lty=nugget$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(nugget$prop.val+nugget$scheme.val, nugget$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=samp.results$all$all.covparams[1,3], col="red")
      # add a legend
      legend("topleft", legend=c(unique(nugget$samp.scheme), "all"),
             lty=c(unique(nugget$line.val),1), lwd=2,
             col=c(rep(1, length(unique(nugget$line.val))),2),
             bty="n")
    }


    if(pred.plots == TRUE) {

      # RMSPE from LOOCV
      all.rmspe <- sqrt(sum((samp.results$all$all.cross.stats$error^2)/length(samp.results$all$all.cross.stats$error)))

      RMSPE <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(RMSPE) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        RMSPE[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        RMSPE[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        RMSPE[i,3:5] <- quantile(subdata$cross.stats[,"RMSPE"], probs=c(0.25, 0.5, 0.75))
      }
      RMSPE[,6] <- as.numeric(as.character(factor(unique(RMSPE$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      RMSPE[,7] <- as.numeric(as.character(factor(RMSPE$samp.scheme, labels=xval)))
      RMSPE[,8] <- RMSPE$scheme.val/min(RMSPE$scheme.val)

      plot(0:(n.prop.sites+1),
           rep(max(samp.results[[samp.results$results.names[1]]]$cross.stats[,"RMSPE"]),
               n.prop.sites+2),
           type="n",
           ylim=c(min(RMSPE$q25, all.rmspe),
                  max(RMSPE$q75, all.rmspe)),
           xlab="Percentage of sites retained", ylab=expression("RMSPE"[LOOCV]),
           main="RMSPE: Observed sites",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- RMSPE$prop.val+RMSPE$scheme.val
      graphics::segments(x0=x0,
                         y0=RMSPE$q25,
                         x1=x0,
                         y1=RMSPE$q75,
                         lty=RMSPE$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(RMSPE$prop.val+RMSPE$scheme.val, RMSPE$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=all.rmspe, col="red")
      # add a legend
      legend("topleft", legend=c(unique(RMSPE$samp.scheme), "all"),
             lty=c(unique(RMSPE$line.val),1), lwd=2,
             col=c(rep(1, length(unique(RMSPE$line.val))),2),
             bty="n")


      # AKSE and prediction error

      # get predicted values
      preds2 <- samp.results$all$all.preds$predict
      N <- length(preds2)
      predsSE <- sqrt(samp.results$all$all.preds$krige.var)
      AKSE_all <- sqrt((sum(predsSE^2))/N)

      # Prediction error with prediction from full network as "true" and prediction
      # from sampled networks as "predicted"
      pred.error <- list()

      for(i in 1:length(samp.results$results.names)) {
        error.mat <- matrix(nrow=length(preds2),
                            ncol=samp.results$call$nsim)
        pred.data <- samp.results[[samp.results$results.names[i]]]
        for(j in 1:ncol(error.mat)) {
          error.mat[,j] <- preds2 - pred.data$preds.value[,j]
        }

        pred.error[[i]] <- sqrt((apply(error.mat^2, 2, sum)/length(preds2)))

      }
      names(pred.error) <- samp.results$results.names

      # plot prediction error
      pred.plot <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(pred.plot) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        pred.plot[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        pred.plot[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        pred.plot[i,3:5] <- quantile(pred.error[[i]], probs=c(0.25, 0.5, 0.75))
      }
      pred.plot[,6] <- as.numeric(as.character(factor(unique(pred.plot$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      pred.plot[,7] <- as.numeric(as.character(factor(pred.plot$samp.scheme, labels=xval)))
      pred.plot[,8] <- pred.plot$scheme.val/min(pred.plot$scheme.val)

      plot(0:(n.prop.sites+1),
           rep(1, n.prop.sites+2),
           type="n",
           ylim=c(min(pred.plot$q25), max(pred.plot$q75)),
           xlab="Percentage of sites retained", ylab="prediction error",
           main="Prediction error at prediction sites",
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- pred.plot$prop.val+pred.plot$scheme.val
      graphics::segments(x0=x0,
                         y0=pred.plot$q25,
                         x1=x0,
                         y1=pred.plot$q75,
                         lty=pred.plot$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(pred.plot$prop.val+pred.plot$scheme.val, pred.plot$q50, pch=19, col="blue", cex=1.5)
      # add a legend
      legend("topleft", legend=c(unique(pred.plot$samp.scheme)),
             lty=c(unique(pred.plot$line.val)), lwd=2,
             col=c(rep(1, length(unique(pred.plot$line.val)))),
             bty="n")


      # AKSE: average kriging standard error ratio
      akse <- list()

      for(i in 1:length(samp.results$results.names)) {
        akse.data <- samp.results[[samp.results$results.names[i]]]
        akse[[i]] <- (sqrt(apply(akse.data$preds.se^2, 2, sum)/N))/AKSE_all
      }
      names(akse) <- samp.results$results.names

      # plot akse ratio
      akse.plot <- data.frame(matrix(nrow=n.samp.scheme*n.prop.sites, ncol=8))
      colnames(akse.plot) <- c("samp.scheme", "prop.sites", "q25", "q50", "q75", "prop.val", "scheme.val", "line.val")
      for(i in 1:length(samp.results$results.names)){
        subdata <- samp.results[[samp.results$results.names[i]]]
        akse.plot[i,1] <- strsplit(samp.results$results.names[i], "[.]")[[1]][1]
        akse.plot[i,2] <- as.numeric(strsplit(samp.results$results.names[i], "[.]")[[1]][2])/100
        akse.plot[i,3:5] <- quantile(akse[[i]], probs=c(0.25, 0.5, 0.75))
      }
      akse.plot[,6] <- as.numeric(as.character(factor(unique(akse.plot$prop.sites), labels=c(length(unique(samp.results$prop.sites)):1))))
      akse.plot[,7] <- as.numeric(as.character(factor(akse.plot$samp.scheme, labels=xval)))
      akse.plot[,8] <- akse.plot$scheme.val/min(akse.plot$scheme.val)

      plot(0:(n.prop.sites+1),
           rep(1, n.prop.sites+2),
           type="n",
           ylim=c(min(akse.plot$q25), max(akse.plot$q75)),
           xlab="Percentage of sites retained", ylab="ratio",
           main=expression("AKSE"[subset]/"AKSE"[all]),
           xaxt="n")
      abline(v=c(1:(n.prop.sites+1)), col="grey")

      axis.vals <- (1:n.prop.sites) + 0.5
      axis(side=1, at=axis.vals, labels=100*samp.results$prop.sites)

      # plot bars for upper and lower quartiles of plotting values
      x0 <- akse.plot$prop.val+akse.plot$scheme.val
      graphics::segments(x0=x0,
                         y0=akse.plot$q25,
                         x1=x0,
                         y1=akse.plot$q75,
                         lty=akse.plot$line.val,
                         lwd=2)
      # plot points for median of plotting values
      points(akse.plot$prop.val+akse.plot$scheme.val, akse.plot$q50, pch=19, col="blue", cex=1.5)
      # add horizontal line for parameter value from all data
      abline(h=1, col="red")
      # add a legend
      legend("topleft", legend=c(unique(akse.plot$samp.scheme)),
             lty=c(unique(akse.plot$line.val)), lwd=2,
             col=c(rep(1, length(unique(akse.plot$line.val)))),
             bty="n")

    }


  }

}
