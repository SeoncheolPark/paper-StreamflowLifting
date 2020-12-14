library(adlift)
library(xts)
library(zoo)
library(riverdist)
library(nlt)
library(foreach)
library(parallel) #for detect cores
library(doParallel)
library(sjmisc)
#function (input, f, LocalPred = LinearPred, neighbours = 1, intercept = TRUE, 
#         closest = FALSE, nkeep = 2, initboundhandl = "reflect", mod = sample(1:length(input), 
#                                                                               (length(input) - nkeep), FALSE), do.W = FALSE, varonly = FALSE) 
  

#data=data[index_sub]; example_network=example_network2; J=10; endpt=68
#adjacency=adjacency_old
#pred = streamPred_S; neigh = 1; int = TRUE; clo = FALSE; keep = 2; rule = "median"
#sd.scale=1; returnall = FALSE; plot.fig = FALSE; plot.individual = FALSE
#pollutant=NULL; polluyear=NULL; plot.thesis=FALSE; do.parallel=TRUE
#verbose=TRUE; ga=FALSE; index_sub_data=index_sub_data_new
#oremovelist=result_forward$removelist; max.subnum=5

nlt_Stream_S <- function(data, example_network, J=10, endpt=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE, plot.fig = FALSE, plot.individual = FALSE, pollutant=NULL, polluyear=NULL, plot.thesis=FALSE, do.parallel=TRUE, verbose=TRUE, ga=FALSE, index_sub_data=NULL, oremovelist=NULL, max.subnum=5){
  #max.subnum: ???
  n <- length(data)
  vec <- matrix(0, J, n - keep)
  deni <- df <- NULL
  aveghat <- matrix(0, 1, n)
  ghatnat <- NULL
  adjacency2 <- as.data.frame(as.matrix(adjacency_old$adjacency))
  if(is.null(endpt)){
    #stop("You should put endpoint (or endsegment) index")
    endpt <- which(example_network@obspoints@SSNPoints[[1]]@point.data$rid==which(adjacency$rid_bid[,2]=="1"))
  }
  if(do.parallel==FALSE){
    for(ii in 1:J){
      if(verbose){
        cat(ii, "...\n")
      }
      if(ga==FALSE){
        v <- sample(setdiff(c(1:n), endpt), (n - keep), FALSE)
        vec[ii, ] <- as.row(v)
      }else if(ga==TRUE & !is.null(index_sub_data) & !is.null(oremovelist)){
        subnum <- 0
        vc <- oremovelist
        index.candidate <- sort(setdiff(oremovelist, endpt))
        while(TRUE){
          sc <- sample(index.candidate, 2, FALSE)
          #(추가 20 Nov, 2020)
          #if( (index_sub_data[sc[1],2]==index_sub_data[sc[2],2]) & abs(data[sc[1]]-data[sc[2]])<0.25  ){
          #if( (index_sub_data[sc[1],2]==index_sub_data[sc[2],2]) &  (adjacency2[as.numeric(paste(names(data)))[sc[1]], as.numeric(paste(names(data)))[sc[2]]]==1 | adjacency2[as.numeric(paste(names(data)))[sc[2]], as.numeric(paste(names(data)))[sc[1]]]==1) ){
          if( (index_sub_data[sc[1],2]==index_sub_data[sc[2],2])){
            #exchange
            subnum <- subnum + 1
            change.index <- c(which(vc==sc[1]), which(vc==sc[2]))
            vc <- replace(vc, change.index, vc[c(change.index[2], change.index[1])])
            index.candidate <- setdiff(index.candidate, sc)
          }
          if(subnum>=max.subnum){
            break
          }
        }
        v <- vc
        vec[ii,] <- as.row(v)
      }
      
      deni <- denoise_Stream_S_perm(data=data, example_network=example_network, endpt=endpt, per=v, adjacency=adjacency, realweights=realweights, TweedData=TweedData, TweedPredPoints=TweedPredPoints, pred = pred, neigh = neigh, int = int, clo = clo, keep = keep, rule = rule, sd.scale=sd.scale, returnall = FALSE, plot.fig =plot.fig, plot.individual = plot.individual, pollutant=pollutant, polluyear=polluyear, plot.thesis=plot.thesis)
      aveghat <- aveghat + as.row(deni)
    }
    aveghat <- aveghat/J
    

  }else{
    #make clusters
    cl <- parallel::makeCluster(detectCores()-1, setup_strategy = "sequential")
    doParallel::registerDoParallel(cl)
    #parallelized
    aveghat <- c()    
    aveghat <- foreach(ii=1:J, .combine=rbind, .export=c('denoise_Stream_S_perm', 'fwtnp_Stream_S_perm', 'as.row', 'initS_stream', 'compute_shreve_weighted', 'getnbrs_fast', 'adjacency_old', 'getnbrs_S_removept', 'initS_stream_removept', 'as.column', 'PointsUpdate_S', 'artlev', 'ebayesthresh', 'invtnp_stream_S', 'getnbrs_S_removept_inv', 'streamPred_S', 'UndoPointsUpdate_S', 'getweights_internal', 'getnbrs_fast_internal', 'str_contains')) %dopar% {
      if(verbose){
        cat(ii, "...\n")
      }
      if(ga==FALSE){
        v <- sample(setdiff(c(1:n), endpt), (n - keep), FALSE)
        vec[ii, ] <- as.row(v)
        
        #print(v)
      }else if(ga==TRUE & !is.null(index_sub_data) & !is.null(oremovelist)){
        subnum <- 0
        vc <- oremovelist
        index.candidate <- sort(setdiff(oremovelist, endpt))
        while(TRUE){
          if(subnum>=max.subnum){
            break
          }
          sc <- sample(index.candidate, 2, FALSE)
          #(추가 20 Nov, 2020)
          #if( (index_sub_data[sc[1],2]==index_sub_data[sc[2],2]) & abs(data[sc[1]]-data[sc[2]])<0.25){
          if( (index_sub_data[sc[1],2]==index_sub_data[sc[2],2]) ){
          #if( (index_sub_data[sc[1],2]==index_sub_data[sc[2],2]) & (adjacency2[as.numeric(paste(names(data)))[sc[1]], as.numeric(paste(names(data)))[sc[2]]]==1 | adjacency2[as.numeric(paste(names(data)))[sc[2]], as.numeric(paste(names(data)))[sc[1]]]==1)  ){
            #exchange
            subnum <- subnum + 1
            change.index <- c(which(vc==sc[1]), which(vc==sc[2]))
            vc <- replace(vc, change.index, vc[c(change.index[2], change.index[1])])
            index.candidate <- setdiff(index.candidate, sc)
          }
        }
        v <- vc
        vec[ii,] <- as.row(v)
      }
      
      denoise_Stream_S_perm(data=data, example_network=example_network, endpt=endpt, per=vc, adjacency=adjacency, realweights=realweights, TweedData=TweedData, TweedPredPoints=TweedPredPoints, pred = pred, neigh = neigh, int = int, clo = clo, keep = keep, rule = rule, sd.scale=sd.scale, returnall = FALSE, plot.fig =plot.fig, plot.individual = plot.individual, pollutant=pollutant, polluyear=polluyear, plot.thesis=plot.thesis)
      #as.row(deni)
      
    }
    #stopCluster
    parallel::stopCluster(cl)
    aveghat <- colSums(aveghat)/J
  }
  #if(do.orig){
  #  
  #}
  if(returnall) {
    return(list(vec = vec, ghatnat = ghatnat, aveghat = aveghat))
  } else {
    return(aveghat)
  }

}


#adjacency_new = adjacency
#neigh = 1; int = TRUE; clo = FALSE; keep = 125; rule = "median"; sd.scale=1; returnall = FALSE; plot.fig = FALSE; plot.individual = FALSE; pollutant=NULL; polluyear=NULL; plot.thesis=FALSE
#adjacency <- adjacency_old

denoise_Stream_S_perm <- function(data, example_network, endpt=NULL, per=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE, plot.fig = FALSE, plot.individual = FALSE, pollutant=NULL, polluyear=NULL, plot.thesis=FALSE){
  if(is.null(endpt)){
    #stop("You should put endpoint (or endsegment) index")
    endpt <- which(example_network@obspoints@SSNPoints[[1]]@point.data$rid==which(adjacency$rid_bid[,2]=="1"))
  }
  if(is.null(per)){
    per <- sample(setdiff(c(1:length(data)), endpt) , (length(data)-keep), FALSE)
    #per <- as.numeric(names(data)[per])
    print(per)
  }
  
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
    #denoise_S 함수와 비교하여 이 부분만 바꾸면 됨
    #out <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = keep, intercept = int, neighbours = neigh, closest = clo, LocalPred = pred, varonly = TRUE, plot.individual = plot.individual)
    out <- fwtnp_Stream_S_perm(data, example_network, endpt=endpt, mod=per, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = keep, intercept = int, neighbours = neigh, closest = clo, LocalPred = pred, varonly = TRUE, plot.individual = plot.individual)
  }else{
    plot.fig=FALSE
    #denoise_S 함수와 비교하여 이 부분만 바꾸면 됨
    #out <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = keep, intercept = int, neighbours = neigh, closest = clo, LocalPred = pred, varonly = TRUE, plot.individual = FALSE, plot.thesis=c(2,14,22,26))
    out <- fwtnp_Stream_S_perm(data, example_network, endpt=endpt, mod=per, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = keep, intercept = int, neighbours = neigh, closest = clo, LocalPred = pred, varonly = TRUE, plot.individual = FALSE, plot.thesis=c(2,14,22,26))
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

fwtnp_Stream_S_perm <- function(data, example_network, endpt=NULL, mod=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE, plot.individual = FALSE, plot.thesis=NULL){
  if(is.null(endpt)){
    #stop("You should put endpoint (or endsegment) index")
    endpt <- which(example_network@obspoints@SSNPoints[[1]]@point.data$rid==which(adjacency$rid_bid[,2]=="1"))
    print(endpt)
  }
  if(is.null(mod)){
    mod <- sample(setdiff(c(1:length(data)), endpt) , (length(data)-nkeep), FALSE)
    #mod <- as.numeric(names(data)[mod])
    print(mod)
  }
  
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
    remove <- mod[j] #adlift와 비교해서 바뀐 점
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
      print(weight_matrix)
      stop("Weight matrix computation is wrong ")
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
  #saveRDS(mod, "mod.RDS")
  return(list(x = input, coeff = coeff, lengths = lengths, 
              lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist, 
              neighbrs = neighbrs, schemehist = schemehist, interhist = interhist, 
              clolist = clolist, gamlist = gamlist, alphalist = alphalist, 
              W = W, v = v, initresult=init_result, mod=mod))
}
  