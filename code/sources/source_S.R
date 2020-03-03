library(adlift)
library(xts)
library(zoo)
library(riverdist)
#밑의 함수들을 돌릴 때 example_network에 upDist, rid가 반드시 있어야 함

denoise_S <- function (data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE, plot.fig = FALSE, plot.individual = FALSE, pollutant=NULL, polluyear=NULL, plot.thesis=FALSE) {
  #denoise_S <- function (x, f, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", returnall = FALSE) 
  #data로 x, f 해결
  newcoeff <- NULL
  ndetlist <- list()
  tclist <- NULL
  
  #오류 처리
  if(is.null(names(data))){
    data <- matrix(data, nrow=1, ncol=length(data))
    colnames(data) <- as.character(c(1:length(data)))
    #names(data) <- as.character(c(1:length(data)))
  }
  if(is.factor(example_network@obspoints@SSNPoints[[1]]@point.data$rid)){
    example_network@obspoints@SSNPoints[[1]]@point.data$rid <- as.numeric(as.character(example_network@obspoints@SSNPoints[[1]]@point.data$rid))
  }
  if(plot.fig==TRUE & plot.individual==TRUE){
    print("we will only print the denoising result")
    plot.individual <- FALSE
  }
  #오류 처리 끝
  
  adjacency_old <- adjacency
  adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))
  #out <- fwtnp(x, f, LocalPred = pred, neighbours = neigh, intercept = int, closest = clo, nkeep = keep, varonly = TRUE)
  #nkeep=2; int=TRUE; neigh=1; clo=FALSE; pred=streamPred_S; rule = "median"; do.W=TRUE; varonly=FALSE;
  if(plot.thesis==FALSE){
    out <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = keep, intercept = int, neighbours = neigh, closest = clo, LocalPred = pred, varonly = TRUE, plot.individual = plot.individual)
  }else{
    plot.fig=FALSE
    out <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = keep, intercept = int, neighbours = neigh, closest = clo, LocalPred = pred, varonly = TRUE, plot.individual = FALSE, plot.thesis=c(2,14,22,26))
  }
  indsd <- sqrt(out$v)
  norcoeff <- out$coeff/indsd
  lr <- out$lengthsremove
  rem <- out$removelist
  al <- artlev(lr, rem)
  levno <- length(al)
  for (i in 1:levno) {
    ndetlist[[i]] <- norcoeff[al[[i]]]
  }
  sd <- mad(ndetlist[[1]])
  sd1 <- mad(ndetlist[[1]], center = 0)
  sd2 <- mad(norcoeff[rem])
  for (i in 1:levno) {
    tclist <- ebayesthresh(ndetlist[[i]], prior = "cauchy", a = NA, sdev = sd.scale*sd, threshrule = rule)
    newcoeff[al[[i]]] <- tclist * indsd[al[[i]]]
  }
  newcoeff[out$pointsin] <- out$coeff[out$pointsin]
  #fhat <- invtnp(X=as.vector(as.numeric(names(data))), newcoeff, out$lengths, lr, out$pointsin,  rem, out$neighbrs, out$schemehist, out$interhist, length(out$x) - keep, int, neigh, clo, LocalPred=streamPred_S)
  fhat <- invtnp_stream_S(X=as.vector(as.numeric(names(data))), newcoeff, out$lengths, lr, out$pointsin,  rem, out$neighbrs, out$schemehist, out$interhist, length(out$x) - keep, int, neighbours=1, clo, LocalPred=streamPred_S,  data, example_network, adjacency=adjacency_old, realweights)
  
  if(plot.fig==TRUE){
    data_predicted <- out$initresult$weight_matrix%*%as.column(data)
    brks <- range(data_predicted)+c(-0.5,0.5); ngrid <- 121
    brks<-seq(brks[1], brks[2], length = ngrid + 1)
    col.nums<-cut(data_predicted, breaks = brks, labels = FALSE)
    palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black"))
    n.cols <- 121
    main.cols<-palette(n.cols)
    
    data_predicted_denoising <- out$initresult$weight_matrix%*%as.column(fhat$coeff)
    brks.denoising <- range(data_predicted_denoising)+c(-0.5,0.5)
    brks.denoising <-seq(brks.denoising[1], brks.denoising[2], length = ngrid + 1)
    col.nums<-cut(data_predicted_denoising, breaks = brks.denoising, labels = FALSE)
    main.cols.denoising<-palette(n.cols)
    
    par(family = 'sans') 
    par(mar=c(3.1,2.1,3.1,1.1))
    par(mfrow=c(1,2))
    
    if(is.null(pollutant) | is.null(polluyear)){
      #(a) plot Raw data
      #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(a) Raw")
      plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(a) Raw")
      quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], data, main="(a) Raw", cex=2, add=T, zlim=range(data), xlab="", ylab="")
      #(b) plot streamflow lifing scheme denoising result
      plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols.denoising[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(b) Streamflow lifting", yaxt='n')
      quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], fhat$coeff, main="(b) Lifting Denoising", cex=2, add=T, zlim=range(fhat$coeff), xlab="", ylab="", add.legend = T, border="black")
    }else{
      #(a) plot Raw data
      #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(a) Raw")
      plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste("(a) Raw, ", pollutant, ", ", polluyear, sep="")  )
      quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], data, main=paste("(a) Raw, ", pollutant, ", ", polluyear, sep=""), cex=2, add=T, zlim=range(data), xlab="", ylab="")
      #(b) plot streamflow lifing scheme denoising result
      plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols.denoising[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste("(b) Streamflow lifting, ", pollutant, ", ", polluyear, sep=""), yaxt='n')
      quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], fhat$coeff, main=paste("(b) Streamflow lifting, ", pollutant, ", ", polluyear, sep=""), cex=2, add=T, zlim=range(fhat$coeff), xlab="", ylab="", add.legend = T, border="black")
    }
  }
  
  if (returnall) {
    return(list(fhat = fhat, w = out$W, indsd = indsd, al = al, sd = sd))
  }
  else {
    return(fhat$coeff)
  }
}

#adjacency_new = adjacency; adjacency = adjacency_old
#X=as.vector(as.numeric(names(data))); coeff=newcoeff; lengths=out$lengths; lengthsremove=lr; pointsin=out$pointsin; removelist=rem; neighbrs=out$neighbrs; schemehist=out$schemehist; interhist=out$interhist; nadd=length(out$x) - keep; intercept=int; neighbours=1; closest=clo; LocalPred=streamPred_S; adjacency=adjacency_old
invtnp_stream_S <- function (X, coeff, lengths, lengthsremove, pointsin, removelist, neighbrs, schemehist, interhist, nadd = length(X) - 2, intercept = TRUE, neighbours = 1, closest = FALSE, LocalPred = LinearPred,  data, example_network, adjacency, realweights){
  #X: 관찰값이 속해있는 segment number
  adjacency_old <- adjacency
  adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))
  
  #오류 처리
  if(is.factor(example_network@obspoints@SSNPoints[[1]]@point.data$rid)){
    example_network@obspoints@SSNPoints[[1]]@point.data$rid <- as.numeric(as.character(example_network@obspoints@SSNPoints[[1]]@point.data$rid))
  }
  
  if (is.list(X)) {
    coeff <- X$coeff
    lengths <- X$lengths
    lengthsremove <- X$lengthsremove
    removelist <- X$removelist
    neighbrs <- X$neighbrs
    pointsin <- X$pointsin
    schemehist <- X$schemehist
    interhist <- X$interhist
    X <- X$x
  }
  X <- as.row(X)
  coeff <- as.row(coeff)
  n <- length(X)
  N <- length(pointsin)
  m <- length(removelist)
  d <- neighbours
  if (nadd > 0) {
    for (j in 1:nadd) {
      cat(paste("Level " ,j, "\n"))
      N <- length(pointsin)
      remove <- removelist[m - j + 1]
      cat("Remove pts: ", remove, "\n",sep=" ")
      lengthrem <- lengthsremove[m - j + 1]
      nbrs <- neighbrs[[m - j + 1]]
      index <- NULL
      index <- match(nbrs, pointsin)
      #decide r (문제: X[remove], X[nbrs]는 전혀 상관 없음)
      #B <- (X[remove] > X[nbrs])
      B <- (remove > pointsin)
      nt <- sum(B)
      r <- nt +1
      #if (nt == 0) {
      #  r <- which(pointsin == nbrs[1])
      #}
      #if (nt == length(nbrs)) {
      #  r <- which(pointsin == nbrs[length(nbrs)]) + 1
      #}
      #if ((nt > 0) & (nt < length(nbrs))) {
      #  r <- which(pointsin == nbrs[nt + 1])
      #}
      if (is.null(schemehist) == FALSE) {
        if (schemehist[m - j + 1] == "Linear") {
          res <- LinearPred(pointsin, X, coeff, nbrs, remove, intercept = interhist[m - j + 1],  neighbours)
        }
        if (schemehist[m - j + 1] == "Quad") {
          res <- QuadPred(pointsin, X, coeff, nbrs, remove, intercept = interhist[m - j + 1], neighbours)
        }
        if (schemehist[m - j + 1] == "Cubic") {
          res <- CubicPred(pointsin, X, coeff, nbrs, remove, intercept = interhist[m - j + 1], neighbours)
        }
      } else {
        #res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours) #바꿔야 함
        #init_result <- initS_stream(X, data, example_network, adjacency_old, realweights, pointsin=pointsin)
        #I <- init_result$I
        #weight_matrix <- init_result$weight_matrix
        
        #계산 속도 향상을 위해 이 부분을 고쳐야
        #after_result <- initS_stream(X, data, example_network, adjacency_old, realweights, pointsin=sort(c(remove, pointsin)))
        #out2 <- getnbrs_S(X, remove, pointsin=sort(c(remove, pointsin)), neighbours, closest, data, example_network, adjacency_old, realweights) #(1) 작성해야
        #getnbrs_S를 대체할 새로운 함수 생각해야
        out2 <- getnbrs_S_removept_inv(X, remove, pointsin, neighbours, closest, data, example_network, adjacency=adjacency_old, realweights, lengths, weight_matrix_obj=NULL)
        #lengthsnew=lengths_new, nbrs=getnbrs_obj$nbrs, weight=weight_imsi
        cat("Neighbors: ", out2$nbrs, "\n",sep=" ")
        #print(out2$nbrs)
        cat("Weights: ", out2$weight, "\n",sep=" ")
        #print(out2$weight)
        #print(lengthrem)
        #print(out2$lengthsnew[match(remove, names(out2$lengthsnew))])
        if(round(lengthrem,4)!=round(out2$lengthsnew[match(remove, names(out2$lengthsnew))],4) ){
          sdata <- list()
          sdata[["X"]]=X
          sdata[["coeff"]]=coeff
          sdata[["lengths"]]=lengths
          sdata[["lengthsremove"]]=lengthsremove
          sdata[["pointsin"]]=pointsin 
          sdata[["removelist"]]=removelist
          sdata[["neighbrs"]]=neighbrs
          sdata[["schemehist"]]=schemehist
          sdata[["interhist"]]=interhist
          sdata[["nadd"]] = nadd
          saveRDS(sdata, "invmod.RDS")
          print(paste("removelist is ", removelist))
          stop(paste("Lengthsnew is wrong: Original ", round(lengthrem,4), ", New: ", round(out2$lengthsnew[match(remove, names(out2$lengthsnew))],4)) )
        }
        if(sum(out2$nbrs != nbrs)!=0){
          stop("Neighbours do not match")
        }
        #res <- LocalPred(pointsin=sort(c(remove, pointsin)), X, coeff, nbrs, remove, intercept, neighbours, weights=out2$weight)
        #만약 새로운 버전으로 한다면
        res <- LocalPred(pointsin=sort(c(remove, pointsin)), X, coeff, nbrs, remove, intercept, neighbours, weights=out2$weight)
      }
      if (length(res) == 2) {
        l <- res[[1]]
      } else {
        l <- res
      }
      gamweights <- l$weights
      #l1 <- PointsUpdate_S(X, coeff, nbrs, index, remove, pointsin,  weights, lengths, weight_mat=out$weight_matrix) #작성해야
      #l1 <- UndoPointsUpdate_S(X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, lengthnew=after_result$I) #바꿔야 함
      #만약 새로운 버전으로 한다면
      l1 <- UndoPointsUpdate_S(X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, lengthnew=out2$lengthsnew)
      coeff <- l1$coeff
      lengths <- l1$lengths
      pred <- sum(as.column(gamweights) * coeff[nbrs])
      coeff[remove] <- coeff[remove] + pred
      removelist <- setdiff(removelist, remove)
      #if (r == 1) {
      #  #lengths <- c(lengthrem, lengths)
      #  pointsin <- c(remove, pointsin)
      #}
      #if (r == (N + 1)) {
      #  #lengths <- c(lengths, lengthrem)
      #  pointsin <- c(pointsin, remove)
      #}
      #if ((r > 1) & (r < (N + 1))) {
      #  #lengths <- c(lengths[1:(r - 1)], lengthrem, lengths[r:N])
      #  pointsin <- c(pointsin[1:(r - 1)], remove, pointsin[r:N])
      #}
      pointsin <- sort(c(remove, pointsin))
    }
  }
  return(list(coeff = coeff, lengths = lengths, lengthsremove = lengthsremove, 
              pointsin = pointsin, removelist = removelist))
}

#UndoPointsUpdate_S: not needed
UndoPointsUpdate_S <- function(X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, lengthnew){
  alpha <- matrix(0, 1, length(nbrs))
  #if ((r > 1) & (r <= N)) {
  if(length(nbrs)>=2){
    alpha <- lengths[index] * lengthrem/(sum(lengths[index]^2))
    coeff[nbrs] <- coeff[nbrs] - alpha * coeff[remove]
    #lengths[index] <- as.row(lengths[index])
    #prod <- gamweights * lengthrem
    #prod <- as.row(prod)
    #lengths[index] <- lengths[index] - prod
  }else{
    #pointsin2 <- sort(c(remove, pointsin))
    q <- which(pointsin == nbrs)
    alpha <- lengthrem/lengths[q]
    coeff[pointsin[q]] <- coeff[pointsin[q]] - alpha * coeff[remove]
  }
    #alpha <- lengths[index] * lengthrem/(sum(lengths[index]^2))
    #coeff[nbrs] <- coeff[nbrs] - alpha * coeff[remove]
    #lengths[index] <- as.row(lengths[index])
    #prod <- gamweights * lengthrem
    #prod <- as.row(prod)
    #lengths[index] <- lengths[index] - prod
  #}
  #if ((r == 1) | (r == (N + 1))) {
  #  q <- which(pointsin == nbrs)
  #  alpha <- lengthrem/lengths[q]
  #  coeff[pointsin[q]] <- coeff[pointsin[q]] - alpha * coeff[remove]
  #  #lengths[q] <- lengths[q] - lengthrem
  #}
  return(list(coeff = coeff, lengths = lengthnew, alpha = alpha))
}

fwtnp_stream_S <- function(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE, plot.individual = FALSE, plot.thesis=NULL){
  #nkeep = 2; intercept = TRUE; initboundhandl = "reflect";  neighbours = 1; closest = FALSE; LocalPred = streamPred_S; do.W = FALSE; varonly = FALSE
  #TweedData, TweedPredPoints는 그림 출력용으로 받는다
  #MyRivernetwork는 없어도될듯
  
  if(!is.null(plot.thesis)){
    par(family = 'sans') 
    par(mar=c(3.1,2.1,2.1,0.275))
    par(mfrow=c(2,4))
    main.character <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)")
    main.character.num <- 0
  }
  
  adjacency_old <- adjacency
  adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))
  #data: spatial vector
  #I <- initint2_stream(example_network, adjacency=adjacency_old, linear=FALSE, MyRivernetwork)$I
  #I_lines <- initint2_stream(example_network, adjacency=adjacency_old, linear=FALSE, MyRivernetwork)$I_lines
  
  #오류 처리
  if(is.factor(example_network@obspoints@SSNPoints[[1]]@point.data$rid)){
    example_network@obspoints@SSNPoints[[1]]@point.data$rid <- as.numeric(as.character(example_network@obspoints@SSNPoints[[1]]@point.data$rid))
  }
  
  input <- data #data의 segment number를 넣는다
  if(!is.null(names(data))){
    X <- as.numeric(names(data)) #X 따로 만들 필요가 있다
  }else{
    #이름이 없는 경우도 대비합시다
    #X: 관찰값이 속해있는 segment number
    X <- example_network@obspoints@SSNPoints[[1]]@point.data$rid
  }
  
  f <- input
  
  X <- as.row(X)
  f <- as.row(f)
  nkeep <- max(nkeep, 1)
  n <- length(X)
  removelist <- NULL
  lengthsremove <- NULL
  neighbrs <- list()
  gamlist <- list()
  alphalist <- list()
  schemehist <- NULL
  interhist <- NULL
  clolist <- NULL
  pointsin <- matrix(1:n, 1, n)
  #pointsin <- pointsin[order(X)]
  coeff <- f
  matno <- n - nkeep
  W <- v <- NULL
  if ((do.W == 1) & (varonly == 1)) {
    varonly <- FALSE
  }
  ex <- do.W + varonly
  if (ex == 1) {
    W <- diag(n)
  }
  if (varonly) {
    v <- rep(1, times = n)
  }
  #multipleplaces
  
  init_result <- initS_stream(X, data, example_network, adjacency_old, realweights, pointsin=pointsin)
  I <- init_result$I
  #print(paste("lengths at inital level is ", I))
  cat("Lengths at inital level is ", I, "\n",sep=" ")
  #print(sum(I))
  weight_matrix <- init_result$weight_matrix
  
  lengths <- init_result$I
  origlengths <- lengths
  
  ##여기서부터 고쳐야
  #LocalPred <- getnbrs_S
  for (j in 1:matno) {
    cat(paste("Level ", (matno-j+1)), "\n")
    remove <- order(lengths)[1]
    if(adjacency_old$rid_bid[example_network@obspoints@SSNPoints[[1]]@point.data$rid[pointsin[remove]],2]=="1"){
      remove <- order(lengths)[2]
    }
    #if(pointsin[remove]==68){
    #  remove <- order(lengths)[2]
    #}
    remove <- pointsin[remove]
    cat(paste("Remove point is ", remove), "\n")
    removelist[j] <- remove
    #out <- getnbrs(X, remove, pointsin, neighbours, closest)
    #out <- getnbrs_S(X, remove, pointsin, neighbours, closest, data, example_network, adjacency_old, realweights) #(1) 작성해야
  
    #만약 계산속도를 빠르게 하고 싶다면
    start_time <- Sys.time()
    #simulated dataset으로 했을 때 getnbrs_S_removept함수가 이웃을 잘못 찾고 있다
    out <- getnbrs_S_removept(X, remove, pointsin, neighbours, closest, data, example_network, adjacency_old, realweights, weight_matrix_obj=weight_matrix)
    end_time <- Sys.time()
    cat(paste("Computation time at level ", (matno-j+1), " is ", round((end_time - start_time),5), " sec "), "\n")
    weight_matrix_old <- weight_matrix
    weight_matrix <- out$weight_matrix$weight_matrix
    if(sum(is.nan(rowSums(weight_matrix)))!=0){
      stop("Weight matrix computation is wrong")
    }
    cat("Neighbors: ", out$n, "\n",sep=" ")
    cat("Upstream Neighbors", out$upstream_nbrs, "\n",sep=" ")
    cat("Downstream Neighbors", out$downstream_nbrs, "\n",sep=" ")
    cat("Weights: ", out$weight, "\n",sep=" ")
    #print(paste("neighbors: ", out$n))
    nbrs <- out$n
    index <- out$index
    #res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours)
    #start_time <- Sys.time()
    res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours, weights=out$weight) #(2) 작성해야
    #end_time <- Sys.time()
    #print(paste("LocalPred at level ", j, " is ", end_time - start_time))
    if (length(res) == 2) {
      l <- res[[1]]
      clolist[j] <- res[[2]][[1]]
      nbrs <- res[[2]][[2]]
      index <- res[[2]][[3]]
    }else {
      l <- res
    }
    neighbrs[[j]] <- nbrs
    weights <- l[[1]]
    pred <- l[[2]]
    if (length(l) == 3) {
      scheme <- NULL
      int <- NULL
      details <- NULL
    } else {
      scheme <- l[[5]]
      int <- l[[4]]
      details <- l[[6]]
    }
    #print(l[[1]])
    #print(l[[2]])
    coeff[remove] <- coeff[remove] - pred
    #l1 <- PointsUpdate(X, coeff, nbrs, index, remove, pointsin, weights, lengths)
    #start_time <- Sys.time()
    #l1에서 out$weight_matrix 활용
    l1 <- PointsUpdate_S(X, coeff, nbrs, index, remove, pointsin, weights, lengths, weight_mat=out$weight_matrix) #작성해야
    #end_time <- Sys.time()
    #print(paste("PointsUpdate_S ", j, " is ", end_time - start_time))
    coeff <- l1$coeff
    lengths <- l1$lengths
    r <- l1$r
    weights <- l1$weights
    N <- l1$N
    alpha <- l1$alpha
    if (ex) {
      if (varonly) {
        W[r, ] <- W[r, ] - colSums(as.vector(weights) * matrix(W[index, ], nrow = length(nbrs)))
        W[index, ] <- W[index, ] + matrix(alpha) %*% W[r, ]
        v[remove] <- sum(W[r, ]^2)
        np <- setdiff(1:length(pointsin), r)
        W <- W[np, ]
      } else {
        W[remove, ] <- W[remove, ] - colSums(as.vector(weights) * matrix(W[nbrs, ], nrow = length(nbrs)))
        W[nbrs, ] <- W[nbrs, ] + matrix(alpha) %*% W[remove, ]
      }
    }
    lengthsremove[j] <- lengths[r]
    gamlist[[j]] <- weights
    alphalist[[j]] <- alpha
    schemehist[j] <- scheme
    interhist[j] <- int
    lengths <- lengths[setdiff(1:length(pointsin), r)]
    pointsin <- setdiff(pointsin, remove)
    
    #print(paste("lengths at level ", (matno-j+1), " is ", lengths))
    cat(paste("Sum of lengths at level ", (matno-j+1), " is ", sum(lengths)), "\n")
    
    #plotting
    if(plot.individual==TRUE){
      data_predicted <- weight_matrix%*%as.column(coeff[pointsin])
      brks <- range(data_predicted)+c(-0.5,0.5); ngrid <- 121
      brks<-seq(brks[1], brks[2], length = ngrid + 1)
      col.nums<-cut(data_predicted, breaks = brks, labels = FALSE)
      palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black"))
      n.cols <- 121
      main.cols<-palette(n.cols)
      
      #data_predicted_denoising <- out$initresult$weight_matrix%*%as.column(fhat$coeff)
      #brks.denoising <- range(data_predicted_denoising)+c(-0.5,0.5)
      #brks.denoising <-seq(brks.denoising[1], brks.denoising[2], length = ngrid + 1)
      #col.nums<-cut(data_predicted_denoising, breaks = brks.denoising, labels = FALSE)
      #main.cols.denoising<-palette(n.cols)
      
      par(family = 'sans') 
      par(mar=c(2.1,2.1,3.1,1.1))
      par(mfrow=c(1,2))
      #(a) coarser level field
      #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(a) Raw")
      plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste("(a) Coarser, level ", (matno-j+1)))
      quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[pointsin,1], example_network@obspoints@SSNPoints[[1]]@point.coords[pointsin,2], coeff[pointsin], main=paste("(a) Coarser, level ", (matno-j+1)), cex=2, add=T, zlim=range(coeff[pointsin]), xlab="", ylab="")
      #(b) detail coeffs
      if(j!=1){
        plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste("(b) Detail, level ", (matno-j+1)), yaxt='n')
        quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),1], example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),2], coeff[sort(removelist)], main=paste("(b) Detail, level ", (matno-j+1)), cex=2, add=T, zlim=range(coeff[sort(removelist)]), xlab="", ylab="", add.legend = T)
        #quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),1], example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),2], coeff[sort(removelist)], main=paste("(b) Detail, level ", j), cex=2, add=T, zlim=c(min(-1, min(coeff[removelist])), max(1, max(coeff[removelist]))), xlab="", ylab="", add.legend = T)
      }
    }
    if(!is.null(plot.thesis)){
      if(j%in%plot.thesis){
        main.character.num <- main.character.num + 1
        
        data_predicted <- -weight_matrix%*%as.column(coeff[pointsin])
        #brks <- range(data_predicted)+c(-0.5,0.5); ngrid <- 121
        brks <- c(-1,1.2); ngrid <- 121
        brks<-seq(brks[1], brks[2], length = ngrid + 1)
        col.nums<-cut(data_predicted, breaks = brks, labels = FALSE)
        palette<-colorRampPalette(c("cyan", "green", "yellow2", "red", "black"))
        n.cols <- 121
        main.cols<-palette(n.cols)
        
        if(j!=1){
          if(main.character.num==1 | main.character.num==5){
            plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste(main.character[main.character.num]," Detail, level ", (matno-j+1)))
            if(main.character.num==1){
              quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),1], example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),2], coeff[sort(removelist)], main=paste(main.character[main.character.num], " Detail, level ", (matno-j+1)), cex=2, add=T, xlab="", ylab="", add.legend = T, nx=8, ny=8, zlim=c(-1.2,1), legend.mar=1.1, legend.width=1.2)
            }else{
              quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),1], example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),2], coeff[sort(removelist)], main=paste(main.character[main.character.num], " Detail, level ", (matno-j+1)), cex=2, add=T, xlab="", ylab="", add.legend = T, nx=16, ny=16, zlim=c(-1.2,1), legend.mar=1.1, legend.width=1.2)
            }
          }else{
            plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste(main.character[main.character.num]," Detail, level ", (matno-j+1)), yaxt='n')
            quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),1], example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),2], coeff[sort(removelist)], main=paste(main.character[main.character.num], " Detail, level ", (matno-j+1)), cex=2, add=T, xlab="", ylab="", add.legend = T, nx=16, ny=16, zlim=c(-1.2,1), legend.mar=1.1, legend.width=1.2)
          }
          #quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),1], example_network@obspoints@SSNPoints[[1]]@point.coords[sort(removelist),2], coeff[sort(removelist)], main=paste(main.character[main.character.num], " Detail, level ", (matno-j+1)), cex=2, add=T, zlim=range(coeff[sort(removelist)]), xlab="", ylab="", add.legend = T, nx=40, ny=40, zlim=c(-1.2,1))
        }
        main.character.num <- main.character.num + 1
        plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=main.cols[col.nums[TweedPredPoints$StreamUnit]], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main=paste(main.character[main.character.num], " Coarser, level ", (matno-j+1)), yaxt='n')
        #quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[pointsin,1], example_network@obspoints@SSNPoints[[1]]@point.coords[pointsin,2], coeff[pointsin], main=paste(main.character[main.character.num], " Coarser, level ", (matno-j+1)), cex=2, add=T, zlim=range(coeff[pointsin]), xlab="", ylab="")
        #(b) detail coeffs
        if(main.character.num%in%c(2,4,6,8)){
          lut <- palette(n.cols)
          scale = (length(lut)-1)/(0.15)
          for (i in 1:(length(lut)-1)) {
            y = (i-1)/scale + 36.875
            rect(127.225,y,127.275,y+1/scale, col=lut[i], border=NA)
          }
          rect(127.225,36.875,127.275,36.875+0.15, col=NA, border=NULL)
          text(127.31, 36.875, "-1")
          text(127.31, 36.875+0.075, "0.1")
          text(127.31, 36.875+0.15, "1.2")
        }
      }
    }
  }
  if (varonly) {
    v[pointsin] <- rowSums(W^2)
    W <- NULL
  }
  N <- length(pointsin)
  return(list(x = input, coeff = coeff, lengths = lengths, 
              lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist, 
              neighbrs = neighbrs, schemehist = schemehist, interhist = interhist, 
              clolist = clolist, gamlist = gamlist, alphalist = alphalist, 
              W = W, v = v, initresult=init_result))
}

PointsUpdate_S <- function(X, coeff, nbrs, index, remove, pointsin,  weights, lengths, weight_mat){
  r <- which(pointsin == remove) #index랑 같은 역할?
  N <- length(pointsin)
  pointsin2 <- setdiff(pointsin, remove)
  pointsin2 <- as.row(pointsin2)
  
  #lengths update
  if ((r >= 2) & (r <= (N - 1))) {
    lengths <- c(weight_mat$I[c(1:(r-1))], lengths[r] , weight_mat$I[c(r:length(weight_mat$I))])
  } else {
    if (r == 1) {
      lengths <- c(lengths[r], weight_mat$I)
    }
    if (r == N) {
      lengths <- c(weight_mat$I, lengths[r])
    }
  }
  lengths <- as.row(lengths)
  
  #coeff update
  alpha <- matrix(0, 1, length(nbrs))
  if (length(nbrs) >= 2) {
    alpha <- lengths[r] * lengths[index]/(sum(lengths[index]^2))
    coeff[pointsin[index]] <- coeff[pointsin[index]] + alpha *  coeff[remove]
  } else {
    q <- which(pointsin == nbrs)
    alpha <- lengths[r]/lengths[q]
    coeff[pointsin[q]] <- coeff[pointsin[q]] + alpha * coeff[remove]
  }

  return(list(coeff = coeff, lengths = lengths, r = r, N = N, weights = weights, alpha = alpha))
}

streamPred_S <- function(pointsin, X, coeff, nbrs, remove, intercept, neighbours, weights){
  #return
  #return(list(weights = weights, pred = pred, coeff = coeff))
  #coeff: 그냥 반환
  #weights: matrix에서 받으면 됨
  #pred: weight_matrix%*%coeff
  weights <- as.row(weights)
  pred <- weights%*%(as.column(coeff[nbrs]))
  return(list(weights = weights, pred = pred, coeff = as.row(coeff[nbrs])))
}

#index 붙이기
getnbrs_S <- function(X, remove, pointsin, neighbours, closest, data, example_network, adjacency, realweights){
  #이것 역시 시간이 많이 걸린다 (2~3분)
  #remove point가 없어진 river network를 생각
  #이 함수의 리턴값 nbrs=nbrs_new, index=index, weight=weight_new, weight_matrix=weight_matrix_reduced
  pointsin2 <- setdiff(pointsin, remove)
  pointsin2 <- as.row(pointsin2)
  
  #remove point가 없을 때의 weight_matrix 작성
  #이 부분을 바꾸는 것이 좋을 것 같다 (왜냐면 apply 함수를 써도 오래걸린다)
  weight_matrix_reduced <- initS_stream(X, data, example_network, adjacency, realweights, pointsin2)
  #remove가 있는 segment: X[remove]
  #remove가 없을 때 nonnegative 값을 갖는 이웃 찾기(nbrs_new에 저장)]
  #index: the indices into pointsin of the neighbours
  
  nbrs_new <- pointsin2[which(weight_matrix_reduced$weight_matrix[X[remove],]!=0)]
  index <- match(nbrs_new, pointsin)
  weight_new <- weight_matrix_reduced$weight_matrix[X[remove],weight_matrix_reduced$weight_matrix[X[remove],]!=0]
  return(list(nbrs=nbrs_new, index=index, weight=weight_new, weight_matrix=weight_matrix_reduced))
}

#out2 <- getnbrs_S_removept_inv(X, remove, pointsin, neighbours, closest, data, example_network, adjacency=adjacency_old, realweights, lengths, weight_matrix_obj=NULL)
getnbrs_S_removept_inv <- function(X, remove, pointsin, neighbours, closest, data, example_network, adjacency=adjacency_old, realweights, lengths, weight_matrix_obj=NULL, use.internal=TRUE){
  if(!is.null(weight_matrix_obj)){
    weight_matrix=weight_matrix_obj
  }
  
  #(2019년 6월 12일 수정)
  #getnbrs_obj <- getnbrs_fast(segment.num=NULL, pointsin=pointsin, example_network=example_network, adjacency=adjacency, upperExtra=FALSE, remove=remove)
  getnbrs_obj <- getnbrs_fast(segment.num=example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove], pointsin=pointsin, example_network=example_network, adjacency=adjacency, upperExtra=FALSE, remove=NULL)
  #getnbrs_obj$nbrs
  #getnbrs_obj$segs
  
  segsnew <- getnbrs_obj$segs
  if(is.null(getnbrs_obj$direct)){
    #(그대로)
  }else if(length(getnbrs_obj$direct)==1){
    if(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[getnbrs_obj$direct] > example_network@obspoints@SSNPoints[[1]]@point.data$upDist[remove] ){
      #이것을 만족하면 이웃은 상류에
      #이웃의 하류를 segsnew에 추가해 준다
      getnbrs_obj2 <- getnbrs_fast(segment.num=example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove], pointsin=setdiff(pointsin,getnbrs_obj$direct), example_network=example_network, adjacency=adjacency, upperExtra=FALSE, remove=NULL)
      segsnew <- c(segsnew, getnbrs_obj2$segs_downstream)
    }else{
      #이것을 만족하면 이웃은 하류에
      #이웃의 상류를 segsnew에 추가해 준다
      getnbrs_obj2 <- getnbrs_fast(segment.num=example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove], pointsin=setdiff(pointsin,getnbrs_obj$direct), example_network=example_network, adjacency=adjacency, upperExtra=FALSE, remove=NULL)
      segsnew <- c(segsnew, getnbrs_obj2$segs_upstream)
    }
  }else if(length(getnbrs_obj$direct)>=2){
    getnbrs_obj2 <- getnbrs_fast(segment.num=example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove], pointsin=setdiff(pointsin,getnbrs_obj$direct), example_network=example_network, adjacency=adjacency, upperExtra=FALSE, remove=NULL)
    #remove_upstream_order <- order(c(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[remove], example_network@obspoints@SSNPoints[[1]]@point.data$upDist[getnbrs_obj$direct]))[1]
    #동점 고려
    remove_upstream_order <- rank(c(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[remove], example_network@obspoints@SSNPoints[[1]]@point.data$upDist[getnbrs_obj$direct]))[1]
    if(remove_upstream_order==1){
      # 낮은 순서대로이므로 remove point가 해당 세그먼트의 가장 아랫쪽에 있다
      # 하류를 추가해주면 됨
      segsnew <- c(segsnew, getnbrs_obj2$segs_downstream)
    }else if(remove_upstream_order==(length(getnbrs_obj$direct)+1)){
      # 낮은 순서대로이므로 remove point가 해당 세그먼트의 가장 윗쪽에 있다
      # 상류를 추가해 주면 됨
      segsnew <- c(segsnew, getnbrs_obj2$segs_upstream)
    }else{
      #가운데에 있으므로 아무것도 안하면 됨
    }
  }
  
  
  #getnbrs_segment(segment.num=getnbrs_obj$segs[1], pointsin, example_network, adjacency=adjacency_old, upperExtra=FALSE)
  
  weight_vec_cand3 <- compute_shreve_weighted(adjacency=adjacency, realweights)
  
  weight_matrix <- matrix(0, nrow=(length(segsnew)+1), ncol=length(pointsin))
  pointsin2 <- c(pointsin, remove)
  #(6월 14일 수정: 중복관찰이 있는 경우에는 removept를 제거하더라도 그 segment에 계속 관찰값이 남아있을 수 있음)
  #weight_matrix_include_removept <- matrix(0, nrow=length(segsnew), ncol=length(pointsin2))
  weight_matrix_include_removept <- matrix(0, nrow=(length(segsnew)+1), ncol=length(pointsin2))
  
  #print(length(segsnew))
  #print(dim(weight_matrix))
  #print(dim(weight_matrix_include_removept))
  
  if(!is.null(segsnew)){
    if(use.internal==FALSE){
      for(i in 1:length(segsnew)){
        #print(i)
        #(2019년 6월 12일 수정)
        #(1) weight_matrix쪽: pointsin 사용했을 시
        #nbrs_imsi <- getnbrs_segment(segment.num=segsnew[i], pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
        nbrs_imsi <- getnbrs_fast(segment.num=segsnew[i], pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
        if(!is.null(nbrs_imsi$direct)){
          weight_matrix[i, match(nbrs_imsi$direct, pointsin)] <- 1
        }else{
          if(!is.null(nbrs_imsi$nbrs_upstream)){
            weight_denom <- weight_vec_cand3$distweight[segsnew[i]]-realweights[segsnew[i]]
            for(j in 1:length(nbrs_imsi$nbrs_upstream)){
              weight_matrix[i,match(nbrs_imsi$nbrs_upstream[j], pointsin) ] <- weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_upstream[j]]]/weight_denom
            }
          }
          if(!is.null(nbrs_imsi$nbrs_downstream)){
            weight_denom <- weight_vec_cand3$distweight[segsnew[i]]
            for(j in 1:length(nbrs_imsi$nbrs_downstream)){
              weight_matrix[i,match(nbrs_imsi$nbrs_downstream[j], pointsin) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_downstream[j]]]
            }
          }
        }
        #normalizing
        weight_matrix[i,] <- weight_matrix[i,]/sum(weight_matrix[i,])
        
        #2019년 6월 12일 수정
        #(1) weight_matrix_include_removept쪽: pointsin2 사용했을 시
        #nbrs_imsi_removept <- getnbrs_segment(segment.num=segsnew[i], pointsin=pointsin2, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
        nbrs_imsi_removept <- getnbrs_fast(segment.num=segsnew[i], pointsin=pointsin2, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
        if(!is.null(nbrs_imsi_removept$direct)){
          weight_matrix_include_removept[i, match(nbrs_imsi_removept$direct, pointsin2)] <- 1
        }else{
          if(!is.null(nbrs_imsi_removept$nbrs_upstream)){
            weight_denom <- weight_vec_cand3$distweight[segsnew[i]]-realweights[segsnew[i]]
            for(j in 1:length(nbrs_imsi_removept$nbrs_upstream)){
              weight_matrix_include_removept[i,match(nbrs_imsi_removept$nbrs_upstream[j], pointsin2) ] <- weight_vec_cand3$distweight[X[nbrs_imsi_removept$nbrs_upstream[j]]]/weight_denom
            }
          }
          if(!is.null(nbrs_imsi_removept$nbrs_downstream)){
            weight_denom <- weight_vec_cand3$distweight[segsnew[i]]
            for(j in 1:length(nbrs_imsi_removept$nbrs_downstream)){
              weight_matrix_include_removept[i,match(nbrs_imsi_removept$nbrs_downstream[j], pointsin2) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi_removept$nbrs_downstream[j]]]
            }
          }
        }
        #normalizing
        weight_matrix_include_removept[i,] <- weight_matrix_include_removept[i,]/sum(weight_matrix_include_removept[i,])
      }
      #print(dim(weight_matrix))
      #print(dim(weight_matrix_include_removept))
    }else if(use.internal==TRUE){
      weight_matrix[c(1:length(segsnew)),] <- t(unlist(mapply(getweights_internal, segsnew=segsnew[c(1:length(segsnew))], MoreArgs = list(X=X, pointsin=pointsin, example_network=example_network, adjacency=adjacency, weight_vec_cand3=weight_vec_cand3, realweights=realweights))))
      weight_matrix_include_removept[c(1:length(segsnew)),] <- t(unlist(mapply(getweights_internal, segsnew=segsnew[c(1:length(segsnew))], MoreArgs = list(X=X, pointsin=pointsin2, example_network=example_network, adjacency=adjacency, weight_vec_cand3=weight_vec_cand3, realweights=realweights))))
      #weight_matrix <- unlist(mapply(getweights_internal, segsnew=segsnew, MoreArgs = list(X=X, pointsin=pointsin, example_network=example_network, adjacency=adjacency, weight_vec_cand3=weight_vec_cand3, realweights=realweights)))
      #weight_matrix_include_removept <- unlist(mapply(getweights_internal, segsnew=segsnew, MoreArgs = list(X=X, pointsin=pointsin2, example_network=example_network, adjacency=adjacency, weight_vec_cand3=weight_vec_cand3, realweights=realweights)))
      
      #print(dim(weight_matrix))
      #print(dim(weight_matrix_include_removept))
      #getweights_internal(X, segsnew=segs, pointsin, example_network, adjacency, weight_vec_cand3, nbrs_imsi, realweights)
    }
    
  }
  
  removeseg <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove]
  #(2019년 6월 12일 수정)
  #nbrs_imsi <- getnbrs_segment(segment.num=removeseg, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
  nbrs_imsi <- getnbrs_fast(segment.num=removeseg, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE, remove=NULL) #해당 segment의 인접 point 계산
  if(is.null(nbrs_imsi$direct)){
    if(!is.null(nbrs_imsi$nbrs_upstream)){
      weight_denom <- weight_vec_cand3$distweight[removeseg]-realweights[removeseg]
      for(j in 1:length(nbrs_imsi$nbrs_upstream)){
        weight_matrix[nrow(weight_matrix),match(nbrs_imsi$nbrs_upstream[j], pointsin) ] <- weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_upstream[j]]]/weight_denom
      }
    }
    if(!is.null(nbrs_imsi$nbrs_downstream)){
      weight_denom <- weight_vec_cand3$distweight[removeseg]
      for(j in 1:length(nbrs_imsi$nbrs_downstream)){
        weight_matrix[nrow(weight_matrix),match(nbrs_imsi$nbrs_downstream[j], pointsin) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_downstream[j]]]
      }
    }
    #normalizing
    weight_matrix[nrow(weight_matrix),] <- weight_matrix[nrow(weight_matrix),]/sum(weight_matrix[nrow(weight_matrix),])
  }else{
    weight_matrix[nrow(weight_matrix), match(nbrs_imsi$direct, pointsin)] <- 1/length(nbrs_imsi$direct)
  }
  #6월 14일 추가
  nbrs_imsi_removept <- getnbrs_fast(segment.num=removeseg, pointsin=pointsin2, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
  if(is.null(nbrs_imsi_removept$direct)){
    stop("Neighborhood computaion is wrong")
    #if(!is.null(nbrs_imsi$nbrs_upstream)){
    #  weight_denom <- weight_vec_cand3$distweight[removeseg]-realweights[removeseg]
    #  for(j in 1:length(nbrs_imsi$nbrs_upstream)){
    #    weight_matrix[nrow(weight_matrix),match(nbrs_imsi$nbrs_upstream[j], pointsin) ] <- weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_upstream[j]]]/weight_denom
    #  }
    #}
    #if(!is.null(nbrs_imsi$nbrs_downstream)){
    #  weight_denom <- weight_vec_cand3$distweight[removeseg]
    #  for(j in 1:length(nbrs_imsi$nbrs_downstream)){
    #    weight_matrix[nrow(weight_matrix),match(nbrs_imsi$nbrs_downstream[j], pointsin) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_downstream[j]]]
    #  }
    #}
    ##normalizing
    #weight_matrix[nrow(weight_matrix),] <- weight_matrix[nrow(weight_matrix),]/sum(weight_matrix[nrow(weight_matrix),])
  }else{
    #weight_matrix[nrow(weight_matrix), match(nbrs_imsi$direct, pointsin)] <- 1/length(nbrs_imsi$direct)
    weight_matrix_include_removept[nrow(weight_matrix_include_removept), match(nbrs_imsi_removept$direct, pointsin2)] <- 1/length(nbrs_imsi_removept$direct)
  }
  
  #example_network@data$Length
  #변경해야
  #I_old <- t(realweights[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])])%*%weight_matrix
  #I_new <- t(realweights[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])])%*%weight_matrix_include_removept
  
  print(dim(realweights[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])]*example_network@data$Length[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])]))
  
  #I 부분이 (한 segment에 관찰점이 3개 있는 경우) 문제가 있는 듯
  I_old <- t(realweights[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])]*example_network@data$Length[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])])%*%weight_matrix
  I_new <- t(realweights[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])]*example_network@data$Length[c(segsnew,example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove])])%*%weight_matrix_include_removept
  
  #if(!is.null(segsnew)){
  #  I_new <- t(realweights[segsnew])%*%weight_matrix_include_removept
  #}else{
  #  I_new <- rep(0, length(pointsin2))
  #}
  
  
  #I_new[length(I_new)] + realweights[example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove]] #이러면 lr[] 값이 나옴
  #위의 값을 sum(I_old-I_new[-length(I_new)])로도 동일하게 얻을 수 있어야
  
  lengths_new <- rep(0, length(pointsin2))
  names(lengths_new) <- sort(pointsin2)
  lengths_new[-which(remove==names(lengths_new))] <- lengths - (I_old-I_new[-length(I_new)])
  lengths_new[which(remove==names(lengths_new))] <- I_new[length(I_new)] #+ realweights[example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove]]
  
  weight_imsi <- weight_matrix[nrow(weight_matrix),match(getnbrs_obj$nbrs, pointsin)]
  if(abs(sum(weight_imsi)-1)>0.001){
    stop("Weight computation is wrong")
  }
  #output: lengthsnew (after_result$I 대체), nbrs (out2$nbrs 대체), weight (out2$weight 대체)
  return(list(lengthsnew=lengths_new, nbrs=getnbrs_obj$nbrs, weight=weight_imsi))
}

#(5월 24일 목표: getnbrs_S_removept, initS_Stream 등의 함수로 바꾸는 것)

getnbrs_S_removept <- function(X, remove, pointsin, neighbours, closest, data, example_network, adjacency, realweights, weight_matrix_obj=weight_matrix){
  weight_matrix=weight_matrix_obj
  #length쪽은 어차피 나중에 업데이트를 하기 때문에 여기서 weight_matrix 부분만 가져가면 될 것 같다
  #이것 역시 시간이 많이 걸린다 (2~3분)
  #remove point가 없어진 river network를 생각
  pointsin2 <- setdiff(pointsin, remove)
  pointsin2 <- as.row(pointsin2)
  
  #remove point가 없을 때의 weight_matrix 작성
  #이 부분을 바꾸는 것이 좋을 것 같다 (왜냐면 apply 함수를 써도 오래걸린다)
  #weight_matrix_reduced <- initS_stream(X, data, example_network, adjacency, realweights, pointsin2)
  weight_matrix_reduced <- initS_stream_removept(X, data, example_network, adjacency, realweights, pointsin, weight_matrix_obj=weight_matrix, remove)
  #weight_matrix, I가 return값으로 오게 된다
  #remove가 있는 segment: X[remove]
  #remove가 없을 때 nonnegative 값을 갖는 이웃 찾기(nbrs_new에 저장)]: 이것이 이웃
  #index: the indices into pointsin of the neighbours
  
  nbrs_new <- pointsin2[which(weight_matrix_reduced$weight_matrix[X[remove],]!=0)] #중복관찰 처리 가능할 것으로 일단 생각됨
  #추가정보 위해 nbrs_upstream, nbrs_downstream 구분
  upstream_remove <- example_network@obspoints@SSNPoints[[1]]@point.data$upDist[remove]
  nbrs_upstream <- nbrs_new[which(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[nbrs_new]>upstream_remove)]
  nbrs_downstream <- nbrs_new[which(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[nbrs_new]<upstream_remove)]
  
  index <- match(nbrs_new, pointsin)
  weight_new <- weight_matrix_reduced$weight_matrix[X[remove],weight_matrix_reduced$weight_matrix[X[remove],]!=0]
  return(list(nbrs=nbrs_new, index=index, weight=weight_new, weight_matrix=weight_matrix_reduced, weight_matrix_original=weight_matrix_obj, upstream_nbrs=nbrs_upstream, downstream_nbrs=nbrs_downstream))
}


initS_stream_removept <- function(X, data, example_network, adjacency, realweights, pointsin2, weight_matrix_obj=weight_matrix, remove, use.internal=TRUE){
  weight_matrix=weight_matrix_obj
  #remove point가 없을 때의 initS_stream 함수 작성: 기본 목적은 계산속도 향상
  #자료를 키웠을 때 time-consuming (1~2분)
  #return: 어떤 matrix
  pointsin <- setdiff(pointsin2, remove)
  remove_ind <- match(remove, pointsin2)
  n.row <- length(realweights)
  n.col <- length(data[pointsin])

  #이 부분이 핵심: remove_ind가 있는 열이 0이 아닌 부분만 걸러내면 된다
  n.row.candidate <- which(weight_matrix[,remove_ind]!=0)
  
  #weight_matrix <- matrix(0, nrow=n.row, ncol=n.col)
  weight_vec_cand3 <- compute_shreve_weighted(adjacency=adjacency, realweights)
  
  initS_stream_rem <- weight_matrix[,-remove_ind]
  weight_matrix <- initS_stream_rem
  
  if(use.internal==FALSE){
    for(i in n.row.candidate){
      #(2019년 6월 12일 수정)
      #nbrs_imsi <- getnbrs_segment(segment.num=i, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
      nbrs_imsi <- getnbrs_fast(segment.num=i, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
      ##compute shreve_similar weight
      if(!is.null(nbrs_imsi$direct)){
        if(length(nbrs_imsi$direct)!=1){
          N_total <- length(nbrs_imsi$direct)
          weight_matrix[i,match(nbrs_imsi$direct, pointsin)] <- 1/N_total
        }else{
          weight_matrix[i,match(nbrs_imsi$direct, pointsin)] <- 1
        }
      }else{
        #같은 segment에 여러 장소가 있음
        if(!is.null(nbrs_imsi$nbrs_upstream)){
          weight_denom <- weight_vec_cand3$distweight[i]-realweights[i]
          for(j in 1:length(nbrs_imsi$nbrs_upstream)){
            weight_matrix[i,match(nbrs_imsi$nbrs_upstream[j], pointsin) ] <- weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_upstream[j]]]/weight_denom
          }
        }
        if(!is.null(nbrs_imsi$nbrs_downstream)){
          weight_denom <- weight_vec_cand3$distweight[i]
          for(j in 1:length(nbrs_imsi$nbrs_downstream)){
            weight_matrix[i,match(nbrs_imsi$nbrs_downstream[j], pointsin) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_downstream[j]]]
          }
        }
      }
      #normalizing
      weight_matrix[i,] <- weight_matrix[i,]/sum(weight_matrix[i,])
    }
  }else if(use.internal==TRUE){
    weight_matrix[n.row.candidate,] <- t(unlist(mapply(getweights_internal, segsnew=n.row.candidate, MoreArgs = list(X=X, pointsin=pointsin, example_network=example_network, adjacency=adjacency, weight_vec_cand3=weight_vec_cand3, realweights=realweights))))
    #print("initS_Stream_removept")
    #print(dim(weight_matrix))
    #print(n.row)
    #print(dim(t(realweights*example_network@data$Length)))
  }
  
  #I <- t(realweights)%*%weight_matrix
  I <- t(realweights*example_network@data$Length)%*%weight_matrix
  return(list(weight_matrix=weight_matrix, I=I))
}


initS_stream <- function(X, data, example_network, adjacency, realweights, pointsin, use.internal=TRUE){
  #자료를 키웠을 때 time-consuming (1~2분)
  #return: 어떤 matrix
  n.row <- length(realweights)
  n.col <- length(data[pointsin])
  weight_matrix <- matrix(0, nrow=n.row, ncol=n.col)
  weight_vec_cand3 <- compute_shreve_weighted(adjacency=adjacency, realweights)
  
  if(use.internal==FALSE){
    for(i in 1:n.row){
      #(2019년 6월 12일 수정)
      #nbrs_imsi <- getnbrs_segment(segment.num=i, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
      nbrs_imsi <- getnbrs_fast(segment.num=i, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
      ##compute shreve_similar weight
      if(!is.null(nbrs_imsi$direct)){
        if(length(nbrs_imsi$direct)!=1){
          N_total <- length(nbrs_imsi$direct)
          weight_matrix[i,match(nbrs_imsi$direct, pointsin)] <- 1/N_total
        }else{
          weight_matrix[i,match(nbrs_imsi$direct, pointsin)] <- 1
        }
      }else{
        if(!is.null(nbrs_imsi$nbrs_upstream)){
          weight_denom <- weight_vec_cand3$distweight[i]-realweights[i]
          for(j in 1:length(nbrs_imsi$nbrs_upstream)){
            weight_matrix[i,match(nbrs_imsi$nbrs_upstream[j], pointsin) ] <- weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_upstream[j]]]/weight_denom
          }
        }
        if(!is.null(nbrs_imsi$nbrs_downstream)){
          weight_denom <- weight_vec_cand3$distweight[i]
          for(j in 1:length(nbrs_imsi$nbrs_downstream)){
            weight_matrix[i,match(nbrs_imsi$nbrs_downstream[j], pointsin) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_downstream[j]]]
          }
        }
      }
      #normalizing
      weight_matrix[i,] <- weight_matrix[i,]/sum(weight_matrix[i,])
    }
  }else if(use.internal==TRUE){
    weight_matrix <- t(unlist(mapply(getweights_internal, segsnew=c(1:n.row), MoreArgs = list(X=X, pointsin=pointsin, example_network=example_network, adjacency=adjacency, weight_vec_cand3=weight_vec_cand3, realweights=realweights))))
    #print("initS_Stream")
    #print(dim(weight_matrix))  
    #print(n.row)
    #print(dim(t(realweights*example_network@data$Length)))
  }
  #I <- t(realweights)%*%weight_matrix
  I <- t(realweights*example_network@data$Length)%*%weight_matrix
  return(list(weight_matrix=weight_matrix, I=I))
}

getnbrs_segment <- function(segment.num, pointsin, example_network, adjacency, upperExtra=FALSE){
  #initS_stream 함수에 사용됨
  #이 함수는 remove point의 고려 없이 해당 segment의 이웃을 매칭시켜주는 함수이다
  #처음에 이걸 모든 segment들에 다 접합시키는 식으로 작성
  target.segment <- adjacency$rid_bid[segment.num,2]
  points.segmentID <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[pointsin]
  points.segmentID <- as.numeric(points.segmentID) #(20190402 수정)
  points.segment <- adjacency$rid_bid[points.segmentID,2] 
  nbrs_downstream <- c(); nbrs_upstream <- c()
  if(segment.num%in%points.segmentID){
    return(list(nbrs_upstream=NULL, nbrs_downstream=NULL, direct=pointsin[which(segment.num==points.segmentID)]))
  }else{
    points.segment.depth <- sapply(points.segment, nchar)
    target.segment.depth <- nchar(target.segment)
    #(1) upward
    #만약 upstream 이웃이 하나도 없으면 upstream 이웃을 찾는다??
    cand_length <- nchar(adjacency$rid_bid[points.segmentID,2])
    nbrs_cand <- which(nchar(adjacency$rid_bid[points.segmentID,2]) > nchar(adjacency$rid_bid[segment.num,2])) #char가 긴 애들
    #nbrs_cand <- setdiff(nbrs_cand,remove)
    #nbrs_cand <- intersect(nbrs_cand, upstream_index)
    
    nbrs_cand <- nbrs_cand[which(sapply(points.segment[nbrs_cand], function(x) substr(x, start=1, stop=nchar(target.segment))==target.segment)==TRUE)]
    if(length(nbrs_cand)!=0){
      while(TRUE){
        #if(length(upstream_candidate)>1){
        #  upstream_candidate <- upstream_candidate[which.min(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[match(upstream_candidate, pointsin)])]
        #}
        diff_vec <- nchar(adjacency$rid_bid[points.segmentID[nbrs_cand],2])- nchar(adjacency$rid_bid[segment.num,2])
        upstream_candidate <- nbrs_cand[which(diff_vec==min(diff_vec))]
        for(ii in 1:length(upstream_candidate)){
          nbrs_upstream <- c(nbrs_upstream, pointsin[upstream_candidate[ii]])
          nbrs_cand <- setdiff(nbrs_cand, upstream_candidate[ii])
          if(length(nbrs_cand)!=0){
            up_upstream <- nbrs_cand[which(sapply(points.segment[nbrs_cand], function(x) substr(x, start=1, stop=nchar(adjacency$rid_bid[points.segmentID[upstream_candidate[ii]],2]))==adjacency$rid_bid[points.segmentID[upstream_candidate[ii]],2]))]
            if(length(up_upstream)!=0){
              nbrs_cand <- setdiff(nbrs_cand, up_upstream)
            }
          }
        }
        if(length(nbrs_cand)==0){
          break
        }
      }
    }
    
    #(2) downward
    ind_imsi <- segment.num
    while(TRUE){
      mouth_index_imsi <- sum(as.matrix(adjacency$adjacency)[ind_imsi,]==1) #하류 지류의 갯수
      if(mouth_index_imsi==0){
        #이경우는 강물을 따라 내려가는 동안 아무러 관측지점을 못만난다는 얘기인데, 고려해봐야 할 듯
        break
      }else if(mouth_index_imsi==1){
        #하류 지류의 갯수가 2개는 아니라고 가정
        mouth_seg_imsi <- which(as.matrix(adjacency$adjacency)[ind_imsi,]==1)
        if(sum(mouth_seg_imsi==points.segmentID)==0){
          #이 경우에는 다음 segment로 넘어가야 한다
          ind_imsi <- mouth_seg_imsi
          #pt_downstream_dist
        }else{
          #이 경우는 해당 segment에 관측장소가 여러 개 있을 경우이다
          #이제 해당 segment안에 가장 가까운 인덱스를 찾는다
          #if(sum(mouth_seg_imsi==npoints.segmentID)==1){
            #해당 segment에 관측장소가 단 한개만 있을 때
            ind <- which(points.segmentID == mouth_seg_imsi)
            nbrs_downstream <- c(nbrs_downstream, pointsin[ind])
          #}else{
          #  #해당 segment에 관측장소가 여러 개 있을 때
          #  ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi) #해당 segment들에 들어있는 obs들의 좌표들
          #  #그러면 해당 segment와 가장 가까운 장소가 어디인지 판단해야 한다
          #  #example_network@lines에서 들어있는 Line들은 상류 -> 하류 순으로 되어있다
          #  pt_upstream <- example_network@lines[[mouth_seg_imsi]]@Lines[[1]]@coords[1,]
          #  nearest_pt_index <- ind[which.min(as.matrix(dist(rbind(pt_upstream,example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]) ))[-1,1])]
          #  nbrs_downstream <- c(nbrs_downstream, pointsin[nearest_pt_index])
          #}
          break
        }
      }
    }
    return(list(nbrs_upstream=nbrs_upstream, nbrs_downstream=nbrs_downstream, direct=NULL))
  }
}




getnbrs_fast <- function(segment.num=NULL, pointsin, example_network=NULL, adjacency, upperExtra=FALSE, remove=NULL){
  #getnbrs_segment 함수 대용으로 이 함수보다 더 빨리 계산하도록 설계된 함수
  #segment.num (제거하려는 remove observation이 위치한 segment number)
  #pointsin (특정 레벨에 들어있는 obs들의 index) (이것들의 segment index가 아님)
  #example_network: initint2_stream에서 작성된 SpatialStreamNetwork object
  #adjacency: adjacency matrix
  #upperExtra: 추후에 neighborhood size가 더 커질 때 활용할 수 있을듯
  #remove: segment.num 인자를 안받고 remove observation의 index를 받아도 처리 가능케 설계
  
  #output에 대한 설명 필요
  #nbrs
  #segs
  #nbrs_upstream
  #nbrs_downstream
  #segs_upstream
  #segs_downstream
  #direct
  
  #segment.num == remove
  #pointsin == remining segements
  #adjacency: adjacency matrices
  if(is.null(remove) & is.null(segment.num)){
    stop("we need segment.num or remove")
  }
  if(!is.null(example_network)){
    target.segment <- adjacency$rid_bid[segment.num,2]
    points.segmentID <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[pointsin]
    points.segmentID <- as.numeric(points.segmentID) #(20190402 수정)
    points.segment <- adjacency$rid_bid[points.segmentID,2] 
    if(!is.null(remove)){
      segment.num=example_network@obspoints@SSNPoints[[1]]@point.data$rid[remove]
    }
  }
  adjacency_full_mat <- as.matrix(adjacency$adjacency)
  rowSums_result <- rowSums(adjacency_full_mat)
  colSums_result <- colSums(adjacency_full_mat)
  nbrs_candidate <- c()
  segs_candidate <- c()
  
  nbrs_lower_candidate <- c()
  nbrs_upper_candidate <- c()
  segs_lower_candidate <- c()
  segs_upper_candidate <- c()
  
  if(segment.num%in%points.segmentID){
    return(list(nbrs=pointsin[which(segment.num==points.segmentID)], segs=NULL, nbrs_upstream=nbrs_upper_candidate, nbrs_downstream=nbrs_lower_candidate, segs_upstream=segs_upper_candidate, segs_downstream=segs_lower_candidate, direct=pointsin[which(segment.num==points.segmentID)]))
  }else{
    
    #(1) 하류
    #print("compute downstream neighbors")
    if(rowSums_result[segment.num]!=0){
      lower_index_old <- segment.num
      while(TRUE){
        lower_index_imsi <- which(adjacency_full_mat[lower_index_old,]!=0)
        if(lower_index_imsi %in% points.segmentID){
          #nbrs_lower_candidate <- c(nbrs_lower_candidate, lower_index_imsi)
          if(length(which(lower_index_imsi==points.segmentID))==1){
            #(1-1) 해당 대상 segment에 관찰값이 1인 경우
            nbrs_lower_candidate <- c(nbrs_lower_candidate, pointsin[which(lower_index_imsi==points.segmentID)])
            #segs_lower_candidate <- c(segs_lower_candidate, lower_index_imsi)
          }else{
            #(1-2) 해당 대상 segment에 관찰값이 2 또는 3인 경우
            pointsin_candidate <- pointsin[which(lower_index_imsi==points.segmentID)]
            #하류니까 example_network@obspoints@SSNPoints[[1]]@point.data$upDist가 큰 지점
            nbrs_lower_candidate <- c(nbrs_lower_candidate, pointsin_candidate[which.max(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[pointsin_candidate])])
          }
          #그럼 관찰값이 있으므로 종료
          break
        }else{
          segs_lower_candidate <- c(segs_lower_candidate, lower_index_imsi)
          lower_index_old <- lower_index_imsi
        }
        if(rowSums_result[lower_index_old]==0){
          break
        }
      }
    }
    
    #(2) 상류
    #print("compute upstream neighbors")
    adjacency_num_remove <- adjacency$rid_bid[segment.num,2]
    #nchar(adjacency_num_remove)
    segnum_depth <- nchar(adjacency_num_remove)
    #which(sapply(adjacency$rid_bid[,2], nchar) >= nchar(adjacency_num_remove))
    #주의사항 : upper_candidate_index는 segment number임
    upper_candidate_index <- which(sapply(adjacency$rid_bid[,2], function(x) substr(x, start=1, stop=nchar(adjacency_num_remove)))==adjacency_num_remove)
    upper_candidate_index <- setdiff(upper_candidate_index, segment.num)
    if(length(upper_candidate_index)!=0){
      while(TRUE){
        #print(upper_candidate_index)
        #upper_index_imsi <- upper_candidate_index[which.min(adjacency$rid_bid[upper_candidate_index,2])]
        upper_index_imsi <- upper_candidate_index[which(sapply(adjacency_old$rid_bid[upper_candidate_index,2], nchar)==min(sapply(adjacency_old$rid_bid[upper_candidate_index,2], nchar)))]
        remove_index_imsi <- c()
        
        
        #위의 for loop 대신
        #upper_imsi_result <- lapply(upper_index_imsi, getnbrs_fast_internal, pointsin=pointsin, points.segmentID=points.segmentID, example_network=example_network, adjacency=adjacency, upper_candidate_index=upper_candidate_index)
        upper_imsi_result <- lapply(upper_index_imsi, getnbrs_fast_internal, pointsin=pointsin, points.segmentID=points.segmentID, example_network=example_network, adjacency=adjacency, upper_candidate_index=upper_candidate_index, nbrs_upper_candidate=nbrs_upper_candidate, segs_upper_candidate=segs_upper_candidate)
        keys <- unique(unlist(lapply(upper_imsi_result, names)))
        upper_imsi_result2 <- do.call(mapply, c(FUN=c, lapply(upper_imsi_result, `[`, keys)))
        #print(upper_imsi_result2)
        #print(class(upper_imsi_result2))
        #print("---")
        #print(upper_imsi_result2[[1]])
        #print(upper_imsi_result2[[2]])
        #print(upper_imsi_result2[[3]])
        
        nbrs_upper_candidate2 <- c()
        segs_upper_candidate2 <- c()
        remove_index_imsi2 <- c()
        
        if(class(upper_imsi_result2)=="matrix"){
          nbrs_upper_candidate2 <- sort(unique(c(upper_imsi_result2[,1])))
          segs_upper_candidate2 <- sort(unique(c(upper_imsi_result2[,2])))
          remove_index_imsi2 <- sort(upper_imsi_result2[,3])
        }else if(class(upper_imsi_result2)=="list"){
          nbrs_upper_candidate2 <- unique(sort(upper_imsi_result2[[1]]))
          segs_upper_candidate2 <- unique(sort(upper_imsi_result2[[2]]))
          remove_index_imsi2 <- sort(upper_imsi_result2[[3]])
        }else if(class(upper_imsi_result2)=="integer"){
          nbrs_upper_candidate2 <- as.integer(as.numeric(sort(unique(c(upper_imsi_result2[1])))))
          segs_upper_candidate2 <- as.integer(as.numeric(sort(unique(c(upper_imsi_result2[2])))))
          remove_index_imsi2 <- as.integer(as.numeric(sort(upper_imsi_result2[3])))
        }
        
        nbrs_upper_candidate <- nbrs_upper_candidate2
        segs_upper_candidate <- segs_upper_candidate2
        remove_index_imsi <- remove_index_imsi2
        
        
        
        #for(j in 1:length(upper_index_imsi)){
        #  if(upper_index_imsi[j] %in% points.segmentID){
        #    #해당 segment에 관찰값이 있으므로 이웃으로 추가
        #    if(length(which(upper_index_imsi[j]==points.segmentID))==1){
        #      #(2-1) 해당 대상 segment에 관찰값이 1인 경우
        #      #nbrs_upper_candidate <- c(nbrs_upper_candidate, upper_index_imsi[j])
        #      nbrs_upper_candidate <- c(nbrs_upper_candidate, pointsin[which(upper_index_imsi[j]==points.segmentID)])
        #      #segs_upper_candidate <- c(segs_upper_candidate, upper_index_imsi[j])
        #      
        #      #print( pointsin[which(upper_index_imsi[j]==points.segmentID)])
        #    }else{
        #      #(2-2) 해당 대상 segment에 관찰값이 2 또는 3인 경우
        #      pointsin_candidate <- pointsin[which(upper_index_imsi[j]==points.segmentID)]
        #      #상류니까 example_network@obspoints@SSNPoints[[1]]@point.data$upDist가 작은 지점
        #      nbrs_upper_candidate <- c(nbrs_upper_candidate, pointsin_candidate[which.min(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[pointsin_candidate])])
        #    }
        #    
        #    #해당 segment 위의 지점들을 모두 제거해준다
        #    upper_index_upper <- intersect(which(sapply(adjacency$rid_bid[upper_candidate_index,2], function(x) substr(x, start=1, stop=nchar(adjacency$rid_bid[upper_index_imsi[j],2])))==adjacency$rid_bid[upper_index_imsi[j],2] ),
        #                                   which(sapply(adjacency$rid_bid[upper_candidate_index,2], nchar) > nchar(adjacency$rid_bid[upper_index_imsi[j],2]) ) )
        #    remove_index_imsi <- c(remove_index_imsi, sort(upper_candidate_index[upper_index_upper]))                               
        #  }else{
        #    #해당 segment에 관찰값이 없으므로 다음으로 넘어간다
        #    segs_upper_candidate <- c(segs_upper_candidate, upper_index_imsi[j])
        #  }
        #}
        
        
        
        #if(!identical(sort(nbrs_upper_candidate), nbrs_upper_candidate2)){
        #  print("---")
        #  print(sort(nbrs_upper_candidate))
        #  print(nbrs_upper_candidate2)
        #  print(class(sort(nbrs_upper_candidate)))
        #  print(class(nbrs_upper_candidate2))
        #  stop("nbrs_upper_candidates are different")
        #}
        #if(!identical(sort(segs_upper_candidate), segs_upper_candidate2)){
        #  print(sort(segs_upper_candidate))
        #  print(segs_upper_candidate2)
        #  stop("segs_upper_candidates are different")
        #}
        #if(!identical(sort(remove_index_imsi), remove_index_imsi2)){
        #  print(sort(remove_index_imsi))
        #  print(remove_index_imsi2)
        #  stop("remove_index_imsis are different")
        #}
        
        remove_index <- sort(unique(c(remove_index_imsi, upper_index_imsi)))
        upper_candidate_index <- setdiff(upper_candidate_index, remove_index) #segment 숫자를 줄인다?
        if(length(upper_candidate_index)==0){
          break
        }
      }
    }
    
    #정리
    nbrs_upper_candidate <- sort(nbrs_upper_candidate)
    nbrs_lower_candidate <- sort(nbrs_lower_candidate)
    segs_upper_candidate <- sort(segs_upper_candidate)
    segs_lower_candidate <- sort(segs_lower_candidate)
    
    #print(nbrs_upper_candidate)
    #print(nbrs_lower_candidate)
    #print(segs_upper_candidate)
    #print(segs_lower_candidate)
    
    nbrs_final <- sort(c(nbrs_lower_candidate, nbrs_upper_candidate))
    segs_final <- sort(c(segs_lower_candidate, segs_upper_candidate))
    return(list(nbrs=nbrs_final, segs=segs_final, nbrs_upstream=nbrs_upper_candidate, nbrs_downstream=nbrs_lower_candidate, segs_upstream=segs_upper_candidate, segs_downstream=segs_lower_candidate, direct=NULL))
  }
}

########################################
#Define some internal functions
#for the fast computations
########################################
getweights_internal <- function(X, segsnew, pointsin, example_network, adjacency, weight_vec_cand3, realweights){
  #this function is made to use apply function within
  #getnbrs_S_removept_inv for the fast computation
  weight_matrix <- rep(0, length(pointsin))
  
  nbrs_imsi <- getnbrs_fast(segment.num=segsnew, pointsin=pointsin, example_network, adjacency=adjacency, upperExtra=FALSE) #해당 segment의 인접 point 계산
  if(!is.null(nbrs_imsi$direct)){
    weight_matrix[match(nbrs_imsi$direct, pointsin)] <- 1
    ##Maybe we have to add...
    
  }else{
    if(!is.null(nbrs_imsi$nbrs_upstream)){
      weight_denom <- weight_vec_cand3$distweight[segsnew]-realweights[segsnew]
      for(j in 1:length(nbrs_imsi$nbrs_upstream)){
        weight_matrix[match(nbrs_imsi$nbrs_upstream[j], pointsin) ] <- weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_upstream[j]]]/weight_denom
      }
    }
    if(!is.null(nbrs_imsi$nbrs_downstream)){
      weight_denom <- weight_vec_cand3$distweight[segsnew]
      for(j in 1:length(nbrs_imsi$nbrs_downstream)){
        weight_matrix[match(nbrs_imsi$nbrs_downstream[j], pointsin) ] <- weight_denom/weight_vec_cand3$distweight[X[nbrs_imsi$nbrs_downstream[j]]]
      }
    }
  }
  #normalized result
  return(weight_matrix/sum(weight_matrix))
}

getnbrs_fast_internal <- function(upper_index_imsi, pointsin, points.segmentID, example_network, adjacency, upper_candidate_index, nbrs_upper_candidate, segs_upper_candidate){
  #nbrs_upper_candidate <- c()
  remove_index_imsi <- c()
  #segs_upper_candidate <- c()
  if(upper_index_imsi %in% points.segmentID){
    #해당 segment에 관찰값이 있으므로 이웃으로 추가
    if(length(which(upper_index_imsi==points.segmentID))==1){
      #(2-1) 해당 대상 segment에 관찰값이 1인 경우
      #nbrs_upper_candidate <- c(nbrs_upper_candidate, upper_index_imsi[j])
      nbrs_upper_candidate <- c(nbrs_upper_candidate, pointsin[which(upper_index_imsi==points.segmentID)])
      #segs_upper_candidate <- c(segs_upper_candidate, upper_index_imsi[j])
    }else{
      #(2-2) 해당 대상 segment에 관찰값이 2 또는 3인 경우
      pointsin_candidate <- pointsin[which(upper_index_imsi==points.segmentID)]
      #상류니까 example_network@obspoints@SSNPoints[[1]]@point.data$upDist가 작은 지점
      nbrs_upper_candidate <- c(nbrs_upper_candidate, pointsin_candidate[which.min(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[pointsin_candidate])])
    }
    #해당 segment 위의 지점들을 모두 제거해준다
    upper_index_upper <- intersect(which(sapply(adjacency$rid_bid[upper_candidate_index,2], function(x) substr(x, start=1, stop=nchar(adjacency$rid_bid[upper_index_imsi,2])))==adjacency$rid_bid[upper_index_imsi,2] ),
                                   which(sapply(adjacency$rid_bid[upper_candidate_index,2], nchar) > nchar(adjacency$rid_bid[upper_index_imsi,2]) ) )
    remove_index_imsi <- c(remove_index_imsi, sort(upper_candidate_index[upper_index_upper]))                               
  }else{
    #해당 segment에 관찰값이 없으므로 다음으로 넘어간다
    segs_upper_candidate <- c(segs_upper_candidate, upper_index_imsi)
  }
  return(list(nbrs_upper_candidate=nbrs_upper_candidate, segs_upper_candidate=segs_upper_candidate, remove_index_imsi=remove_index_imsi))
}
