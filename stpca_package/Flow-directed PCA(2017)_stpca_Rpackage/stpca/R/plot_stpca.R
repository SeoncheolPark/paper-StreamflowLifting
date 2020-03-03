#' Visualize (weighted) PCA results.
#'
#' Visualize the results from applying PCA, adjusted for spatial and/or temporal
#' weights if appropriate.  This will plot the results from the output of
#' \code{\link{stpca}}.
#'
#' The type of plot to use will depend on how many principal components are of
#' interest.  If only 1 PC is of interest then \code{plots = "map"} should be
#' used, or \code{plots = "biplot"} if a pair of PC's  should be plotted.  For 3
#' or more PC's then \code{plots = "glyph"} is recommended.
#'
#' \code{"map"} will display loadings (S-mode PCA) or scores (T-mode PCA) for a
#' single PC.  Color breaks are specified using the \code{probs} argument.
#'
#' \code{"biplot"} displays the standard biplot, commonly used to investigate
#' the results of PCA and displays both scores and loadings on the same plot,
#' for any 2 PC's specified using \code{pc.biplot}.  Arrows represent loadings
#' and points represent principal component scores.
#'
#' \code{"glyph"} is based on \code{\link[GWmodel]{glyph.plot}} in the GWmodel
#' package, and displays a map of the river network (if \code{river} is
#' specified) and glyphs showing the loadings (S-mode PCA) or scores (T-mode
#' PCA) for the principal components specified using \code{pc.glyph}.  The
#' lengths of the spokes are relative to each other and the spoke at the 12
#' o'clock position represents the first PC specified in \code{pc.glyph}, and
#' other PC's are represented moving round clockwise.  Blue indicates positive
#' values and red indicates negative values.  A message will appear when glyph
#' plots are produced prompting the user to check that the matrix of coordinates
#' has monitoring sites in the same order as in the rows (\code{pca.mode =
#' "Tmode"}) or columns (\code{pca.mode = "Smode"}) of the dataframe used for
#' \code{\link{stpca}}.  Incorrect ordering will result in glyphs being
#' displayed at the wrong locations.  If a glyph plot is produced for fewer than
#' three PC's then the plot will be produced but a message suggesting other
#' plots might be more appropriate will also appear.  Glyph plots and code are
#' based on Isabella et. al. (2015).
#'
#' \code{"ts"} displays a time series plot of the principal components when
#' \code{pca.mode = "Smode"} in \code{\link{stpca}}.  This plot is not produced
#' for \code{pca.mode = "Tmode"} since a time series of loadings has no
#' meaningful interpretation.  Rembember that since loadings can be positive or
#' negative, the temporal pattern described by each PC might the negative of
#' what is displayed on the plot.
#'
#' \code{"meanPM"} displays the mean time series for all sites, +/- a multiple
#' (meanPM.factor) of the principal components specified using \code{pc.meanPM}.
#' Red symbols are the mean time series + (meanPM.factor*scores), blue symbols
#' are the mean time series - (meanPM.factor*scores).  This plot shows the
#' variation around the mean captured by each principal component.  The pattern
#' might indicate that a principal component describes a shift in the mean,
#' where some monitoring sites have generally higher/lower values than the mean.
#' The principal component might describe a dampening of the mean signal, where
#' the seasonal pattern oscillates at a lower frequency than the mean pattern.
#' Another possibility is that the prinicpal compnent captures a time lag at
#' some sites compared to the mean.  This list is not exhaustive but has been
#' included to provide some suggestions for interpretation.
#'
#' @importFrom GWmodel glyph.plot
#' @importFrom grDevices colorRampPalette
#' @importFrom sp bpy.colors
#' @importFrom graphics layout
#'
#' @param x character string.  This should specify the name of the object in
#'   which the results of \code{\link{stpca}} are stored.
#' @param plots character string.  Specifies the type of plot to be produced.
#'   Can include any combination of c("map", "biplot", "glyph", "ts", "meanPM").  See
#'   Details for further information.
#' @param river character string.  Name of object containing river network
#'   shapefile.  Used for plotting results on a map of the river network where
#'   \code{plots = c("map", "glyph")}.
#' @param coords matrix of 2 columns.  The columns should contain the x- and
#'   y-axis coordinates for the monitoring site locations respectively and MUST
#'   be in the same order as monitoring sites are specified in rows
#'   (\code{pca.mode = "Tmode"}) or columns (\code{pca.mode = "Smode"}) in the
#'   dataframe specified in \code{\link{stpca}}.
#' @param pc.map vector (of length one).  Specifies which principal component to
#'   plot on a map.  Default is PC1.
#' @param probs vector.  A vector of probabilities used to calculate color
#'   breaks when \code{plots = "map"}.  Values should be in the range [0,1].
#' @param pc.biplot vector.  Specifies the principal components to plot when
#'   \code{plots = "biplot"}.  Default is PC1 on x-axis and PC2 on y-axis.
#' @param biplot.factor scalar.  A scaling value used to produce the biplot.
#'   Higher values make arrows shorter.
#' @param pc.glyph vector.  Specifies the principal components to plot when
#'   \code{plots = "glyph"}.  It is recommended that the glyph plot is used only
#'   to display 3 or more principal components.
#' @param r1 scalar.  Used to control relative lengths of spokes on glyphs in
#'   glyph plot.  Larger values decrease the size of the glyphs.  Default is 10.
#' @param pc.ts vector.  Specifies which principal components should have
#'   principal components/scores plotted as time series.  Should only be used
#'   when \code{pca.mode = "Smode"} in \code{\link{stpca}}.  Default is
#'   c(1,2,3).
#' @param pc.meanPM vector.  Specifies which principal components should be
#'   plotted.  Should only be used when \code{pca.mode = "Smode"} in
#'   \code{\link{stpca}}.  Default is c(1,2,3).
#' @param meanPM.factor scalar.  Used to control the contrast between the mean
#'   time series over all sites and +/- a multiple of the principal
#'   components/scores.  See Details.
#'
#' @return Produces the specified plot(s) for the principal components specified
#'   and for each weighting scheme used in \code{\link{stpca}}.
#'
#' @references
#' Biplot code is based on:
#'
#' \url{http://steviep42.bitbucket.org/YOUTUBE.DIR/BB_phys_stats_ex1.R}
#' and
#' \url{https://www.youtube.com/watch?v=I5GxNzKLIoU}.
#'
#' Glyph plot based on:
#'
#' Isabella Gollini, Binbin Lu, Martin Charlton, Christopher Brunsdon, Paul
#' Harris (2015).  GWmodel: An R Package for Exploring Spatial Heterogeneity
#' Using Geographically Weighted Models.   Journal of Statistical Software,
#' 63(17), 1-50. \url{http://www.jstatsoft.org/v63/i17/}.
#'
#' @author Kelly Gallacher, \email{kelly_gallacher@@hotmail.com}
#'
#' @examples
#' library(stpca)
#' library(maptools) #for reading in shapefile
#' library(SSN) #to import SSN object
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
#'                        spatial.wt = weightS, temporal.wt = weightT,
#'                        pca.wt = c("unweighted", "spatial", "temporal",
#'                                   "spatiotemporal"))
#'
#' ## get the river network shapefile, in this case the edges information in
#' ## demoNet
#' demoNet.path <- system.file("demoSSN/demoNet.ssn", package="stpca")
#' river.path <- paste(demoNet.path, "edges.shp", sep="/")
#' river <- readShapeLines(river.path)
#' demoNet <- importSSN(demoNet.path)
#' data1 <- getSSNdata.frame(demoNet, "Obs")
#' coords <- c("NEAR_X", "NEAR_Y")
#' coords <- data1[,coords]
#'
#' ## plot the results
#'
#' # unweighted T-mode PCA
#' plot_stpca(x = Tmode.pca.uw, plots = "map", coords = coords,
#'            river = river, pc.map = 1, probs = seq(0, 1, by=0.5))
#' plot_stpca(x = Tmode.pca.uw, plots = "biplot")
#' plot_stpca(x = Tmode.pca.uw, plots = "glyph", pc.glyph=1:2,
#'            river = river, coords = coords, r1=5)
#'
#'
#' # spatial, temporal, and spatiotemporal weighted S-mode PCA
#' plot_stpca(x = Smode.pca.all, plots = "map", coords = coords,
#'            river = river)
#' plot_stpca(x = Smode.pca.all, plots = "biplot")
#' plot_stpca(x = Smode.pca.all, plots = "glyph", river = river,
#'            coords = coords)
#' plot_stpca(x = Smode.pca.all, plots = "ts", pc.ts = 1:3)
#' plot_stpca(x = Smode.pca.all, plots = "meanPM")
#'
#' @export
plot_stpca <- function(x,
                       plots = c("map", "biplot", "glyph", "ts", "meanPM"),
                       pc.map = 1, probs = c(0, 0.25, 0.75, 1),
                       pc.biplot = c(1,2), biplot.factor = 100,
                       pc.glyph = c(1,2,3), river, r1 = 10, coords,
                       pc.ts = c(1,2,3), pc.meanPM = c(1,2,3),
                       meanPM.factor = 0.05) {

  pca.mode <- x$pca.mode

  ## biplot
  if("biplot" %in% plots) {

    nplot <- length(x$pca.wt)

    if(nplot==4) {
      par(mfrow=c(2,2))
    }
    else {
      par(mfrow=c(1,nplot))
    }

    pca.type <- x$pca.wt
    xlims <- matrix(nrow=2, ncol=length(pca.type))
    ylims <- matrix(nrow=2, ncol=length(pca.type))

    for(i in 1:nplot) {
      pca.type1 <- pca.type[i]
      scores <- x[[pca.type1]]$scores
      rownames(scores) <- x[["row.names"]]
      loadings <- x[[pca.type1]]$loads
      rownames(loadings) <- x[["col.names"]]
      var.val <- x[[pca.type1]]$var
      sd.val = sqrt(var.val)

      xlims[,i] <- range(scores[,pc.biplot[1]]/sd.val[pc.biplot[1]], loadings[,pc.biplot[1]]*sd.val[pc.biplot[1]]/biplot.factor*1.2)
      ylims[,i] <- range(scores[,pc.biplot[2]]/sd.val[pc.biplot[2]], loadings[,pc.biplot[2]]*sd.val[pc.biplot[2]]/biplot.factor*1.2)
    }

    for(i in 1:nplot) {

      pca.type1 <- pca.type[i]
      scores <- x[[pca.type1]]$scores
      rownames(scores) <- x[["row.names"]]
      loadings <- x[[pca.type1]]$loads
      rownames(loadings) <- x[["col.names"]]

      var.val <- x[[pca.type1]]$var
      var.prop <- x[[pca.type1]]$var.prop

      pc1 = 100*round(var.prop[pc.biplot[1]],digits=2)
      pc2 = 100*round(var.prop[pc.biplot[2]],digits=2)
      xlab=paste("PC", pc.biplot[1], ": ", pc1, "%",sep="")
      ylab=paste("PC", pc.biplot[2], ": ", pc2, "%",sep="")

      # Multiply the scaled data by the scores (principal components)
      sd.val = sqrt(var.val)

      scores.min = min(scores[,pc.biplot])
      scores.max = max(scores[,pc.biplot])

#       xlim <- range(scores[,pc.biplot[1]]/sd.val[pc.biplot[1]], loadings[,pc.biplot[1]]*sd.val[pc.biplot[1]]/biplot.factor*1.2)
#       ylim <- range(scores[,pc.biplot[2]]/sd.val[pc.biplot[2]], loadings[,pc.biplot[2]]*sd.val[pc.biplot[2]]/biplot.factor*1.2)

      xlim <- c(min(xlims[1,]), max(xlims[2,]))
      ylim <- c(min(ylims[1,]), max(ylims[2,]))

      plot(scores[,pc.biplot[1]]/sd.val[pc.biplot[1]],scores[,pc.biplot[2]]/sd.val[pc.biplot[2]],
           main=paste("Weighting =", pca.type1),
           xlab=xlab,ylab=ylab,type="n",
           xlim=xlim, ylim=ylim,
           col.axis="blue",
           cex.main=1.5,
           cex.lab=1.5,
           cex.axis=1.5)
      # First plot the variables as vectors
      arrows(0,0,loadings[,pc.biplot[1]]*sd.val[pc.biplot[1]]/biplot.factor,
             loadings[,pc.biplot[2]]*sd.val[pc.biplot[2]]/biplot.factor,
             length=0.1, lwd=2,angle=20, col="red")
      text(loadings[,pc.biplot[1]]*sd.val[pc.biplot[1]]/biplot.factor*1.2,
           loadings[,pc.biplot[2]]*sd.val[pc.biplot[2]]/biplot.factor*1.2,
           labels=rownames(loadings), col="red", cex=1.2)
      # Second plot the scores as points
      text(scores[,pc.biplot[1]]/sd.val[pc.biplot[1]],scores[,pc.biplot[2]]/sd.val[pc.biplot[2]],
           labels=rownames(scores),col="blue", cex=0.7)
      abline(0,0,col="black")
      abline(v=0,col="green")
    }
  }

  ## glyph plot
  if("glyph" %in% plots) {

    if(length(pc.glyph) < 3) {
      message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
    }

    nplot <- length(x$pca.wt)

    if(nplot==4) {
      par(mfrow=c(2,2))
    }
    else {
      par(mfrow=c(nplot,1))
    }

    pca.type <- x$pca.wt

    for(i in 1:nplot) {

      if(x$pca.mode=="Smode") {
        glyph.data <- x[[pca.type[i]]]
        loadings <- glyph.data$loads[,pc.glyph]
        plot(coords[,1], coords[,2], type="n",
             xlab=names(coords)[1], ylab=names(coords)[2],
             main=paste("Loadings: weight =", pca.type[i]),
             cex.main=1.5, cex.lab=1.5, cex.axis=1.5,
             xlim=c(river@bbox[1,]),
             ylim=c(river@bbox[2,]))
        lines(river, col="lightgrey")
        GWmodel::glyph.plot(as.matrix(loadings), coords,
                            r1=r1, add=TRUE)
      }

      if(x$pca.mode=="Tmode") {
        glyph.data <- x[[pca.type[i]]]
        scores <- glyph.data$scores[,pc.glyph]
        plot(coords[,1], coords[,2], type="n",
             xlab=names(coords)[1], ylab=names(coords)[2],
             main=paste("Scores: weight =", pca.type[i]),
             cex.main=1.5, cex.lab=1.5, cex.axis=1.5,
             xlim=c(river@bbox[1,]),
             ylim=c(river@bbox[2,]))
        lines(river, col="lightgrey")
        GWmodel::glyph.plot(as.matrix(scores), coords,
                            r1=r1, add=TRUE)
      }
    }

    message("For plots = glyph: did you check that coordinates for monitoring site ID's are in \nthe same order as data used for stpca()?  \nIf not then glyphs might not correspond to the correct monitoring site.")
  }

  ## ts plot (time series of scores for Smode)
  if("ts" %in% plots) {

    if(x$pca.mode=="Tmode") {
      message("For plots = ts: ts plot should only be used when pca.mode = Smode")
    }

    if(x$pca.mode=="Smode") {

      nplot <- length(pc.ts)
      pca.type <- x$pca.wt

      par(mfcol=c(nplot, length(pca.type)))

      min.y <- numeric(length(pca.type))
      max.y <- numeric(length(pca.type))
      for(i in 1:length(pca.type)) {
        min.y[i] <- min(x[[pca.type[i]]]$scores[,pc.ts])
        max.y[i] <- max(x[[pca.type[i]]]$scores[,pc.ts])
      }

      for(i in 1:length(pca.type)) {

        ts.data <- x[[pca.type[i]]]

        for(j in pc.ts) {
          plot(1:nrow(ts.data$scores), ts.data$scores[,j], type="l", xaxt="n",
               ylim=c(min(min.y), max(max.y)),
               xlab="Date", ylab="score",
               main=paste("weight = ", pca.type[i], ", PC = ", j, sep=""))
          axis(side=1, at=round(quantile(1:nrow(ts.data$scores), probs=c(0,0.33,0.66,1))),
               labels=rownames(ts.data$scores)[round(quantile(1:nrow(ts.data$scores), probs=c(0,0.33,0.66,1)))])
        }
      }
    }
  }


  ## mean time series plus/minus principal component
  if("meanPM" %in% plots) {

    if(x$pca.mode=="Tmode") {
      message("For plots = meanPM: meanPM plot should only be used when pca.mode = Smode")
    }

    if(x$pca.mode=="Smode") {


      nplot <- length(pc.meanPM)
      pca.type <- x$pca.wt

      par(mfcol=c(nplot, length(pca.type)))

      min.y <- numeric(length(pca.type))
      max.y <- numeric(length(pca.type))
      for(i in 1:length(pca.type)) {
        meanPM.data <- x[[pca.type[i]]]
        min.y[i] <- min((x[["mean.timeseries"]] + (meanPM.factor*meanPM.data$scores[,i])),
                        (x[["mean.timeseries"]] - (meanPM.factor*meanPM.data$scores[,i])))
        max.y[i] <- max((x[["mean.timeseries"]] + (meanPM.factor*meanPM.data$scores[,i])),
                        (x[["mean.timeseries"]] - (meanPM.factor*meanPM.data$scores[,i])))
      }

      for(i in 1:length(pca.type)) {

        meanPM.data <- x[[pca.type[i]]]

        for(j in pc.meanPM) {
          plot(1:nrow(meanPM.data$scores), x[["mean.timeseries"]], type="l", xaxt="n",
               ylim=c(min(min.y), max(max.y)),
               xlab="Date", ylab="score",
               main=paste("weight = ", pca.type[i], ", PC = ", j, sep=""), lwd=2)
          lines(1:nrow(meanPM.data$scores), x[["mean.timeseries"]])
          points(1:nrow(meanPM.data$scores), x[["mean.timeseries"]] + (meanPM.factor*meanPM.data$scores[,j]),
                pch="+", col="red")
          points(1:nrow(meanPM.data$scores), x[["mean.timeseries"]] - (meanPM.factor*meanPM.data$scores[,j]),
                 pch="-", col="blue")
          axis(side=1, at=round(quantile(1:nrow(meanPM.data$scores), probs=c(0,0.33,0.66,1))),
               labels=rownames(meanPM.data$scores)[round(quantile(1:nrow(meanPM.data$scores), probs=c(0,0.33,0.66,1)))])
        }
      }
    }
  }



  ## map plot
  if("map" %in% plots) {

    nplot <- length(x$pca.wt)

    pca.type <- x$pca.wt

    for(i in 1:nplot) {

      if(x$pca.mode=="Smode") {

        title <- paste("Loadings PC", pc.map, "\n weight =", pca.type[i])

        map.data <- x[[pca.type[i]]]
        loadings <- map.data$loads[,pc.map]

        n.colors <- length(probs)-1
        data.vec <- sort(loadings)
        quants <- quantile(data.vec, probs=probs)
        getPalette <- grDevices::colorRampPalette(sp::bpy.colors(n.colors))(n.colors)

        par(mar=c(5, 4, 4, 2) + 0.1)
        graphics::layout(matrix(1:2,ncol=2), widths = c(0.85,0.15), heights = c(1,1))
        plot(coords[,1], coords[,2],
             type="n",
             xlab="Easting", ylab="Northing", main=title,
             cex.main=1.5,
             cex.axis=1.5,
             cex.lab=1.5,
             xlim=c(river@bbox[1,]),
             ylim=c(river@bbox[2,]))
        lines(river)
        points(coords[,1], coords[,2],
               pch=19, cex=2,
               col=getPalette[cut(loadings, breaks=quants)])
        par(mar=c(5.1,0.4,4.1,0.5))
        getPalette <- grDevices::colorRampPalette(sp::bpy.colors(100))(100)
        plot(NA,type="n",ann=FALSE,xlim=c(0,2.1),ylim=c(0,101),xaxt="n",yaxt="n",bty="n")
        rect(
          1,
          c(0:99),
          1.5,
          c(1:100),
          col=c(getPalette), border=NA
        )
        text(x=0.5, y=((0:3)*33), labels=round(quants,2), cex=1.5)
      }

      if(x$pca.mode=="Tmode") {

        title <- paste("Scores PC", pc.map, "\n weight =", pca.type[i])

        map.data <- x[[pca.type[i]]]
        scores <- map.data$scores[,pc.map]

        n.colors <- length(probs)-1
        data.vec <- sort(scores)
        quants <- quantile(data.vec, probs=probs)
        getPalette <- grDevices::colorRampPalette(sp::bpy.colors(n.colors))(n.colors)

        par(mar=c(5, 4, 4, 2) + 0.1)
        graphics::layout(matrix(1:2,ncol=2), widths = c(0.85,0.15), heights = c(1,1))
        plot(coords[,1], coords[,2],
             type="n",
             xlab="Easting", ylab="Northing", main=title,
             cex.main=1.5,
             cex.axis=1.5,
             cex.lab=1.5,
             xlim=c(river@bbox[1,]),
             ylim=c(river@bbox[2,]))
        lines(river)
        points(coords[,1], coords[,2],
               pch=19, cex=2,
               col=getPalette[cut(scores, breaks=quants)])
        par(mar=c(5.1,0.4,4.1,0.5))
        getPalette <- grDevices::colorRampPalette(sp::bpy.colors(100))(100)
        plot(NA,type="n",ann=FALSE,xlim=c(0,2.1),ylim=c(0,101),xaxt="n",yaxt="n",bty="n")
        rect(
          1,
          c(0:99),
          1.5,
          c(1:100),
          col=c(getPalette), border=NA
        )
        a <- n.colors
        b <- 100/n.colors
        text(x=0.5, y=((0:a)*b), labels=round(quants,2), cex=1.5)
      }
    }
  }
}



