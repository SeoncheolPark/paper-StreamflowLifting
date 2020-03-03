invtnp_stream_ST <- function (X=NULL, data, example_network, stream_index, adjacency, st_weight=c(0.995,0.005), MyRivernetwork, coeff, lengths, lengthsremove, pointsin, pointsin_new, pointsin_new_index, removelist, 
            neighbrs, schemehist, interhist, nbrs_total_list, nbrs_spatio_temporal_return_list, nadd = length(X) - 2, intercept = TRUE, 
            neighbours = 1, closest = FALSE, LocalPred = streamPred_ST) 
{
  #추가: data, example_network, stream_index, adjacency, st_weight=c(0.995,0.005), MyRivernetwork
  #추가: nbrs_total_list, nbrs_spatio_temporal_return_list
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
  #my addition
  n.ind <- stream_index[1] #stream_index의 첫 번째 원소를 제거하려는 시계열로
  n.subind <- stream_index[-1] #stream_index의 첫 번째 원소를 제외한 나머지 원소들의 시계열을 이웃들의 시계열로 정한다
  X_subind <- list()
  for(i in 1:length(n.subind)){
    X_subind[[i]] <- which(!is.na(data$normaldata[,n.subind[i]]))
  }
  #end of my addition
  X <- which(!is.na(data$normaldata[,n.ind]))
  X <- as.row(X)
  n <- length(X) #n=70
  if(is.null(nadd)){
    nadd <- n- length(pointsin)
  }
  X_temporal <- X
  #X를 extended vertion으로 다시 정의
  coeff <- as.row(coeff)
  N <- length(pointsin)
  m <- length(removelist)
  d <- neighbours
  #몇가지를 새로 생성하자
  #pointsin_new <- c(pointsin)
  length_X_temporal <- length(which(!is.na(data$normaldata[,n.ind])))
  length_pointsin_new_right <- length(pointsin_new)-length(pointsin)
  
  if (nadd > 0) {
    for (j in 1:nadd) {
      N <- length(pointsin) #N=2,3,...
      remove <- removelist[m - j + 1] #(m-j+1)=68,67,66,...
      lengthrem <- lengthsremove[m - j + 1] #제거된 점의 time series length
      #nbrs <- neighbrs[[m - j + 1]]
      nbrs_original <- neighbrs[[m - j + 1]]
      nbrs <- nbrs_original[nbrs_original<=n] #n=70이라고 가정
      #나의 추가 내용
      nbrs_total <- nbrs_total_list[[m-j+1]]
      nbrs_spatio_temporal <- nbrs_spatio_temporal_return_list[[m-j+1]]
      #나의 추가 내용 끝
      index <- NULL
      index <- match(nbrs, pointsin) #바꿔야하나?
      #나의 추가
      index_extended <- match(nbrs_original, pointsin_new)
      lengths_original <- lengths
      lengths <- lengths_original[c(1:length(pointsin))]
      #나의 추가 끝
      
      B <- (X[remove] > X[nbrs])
      nt <- sum(B)
      if (nt == 0) {
        r <- which(pointsin == nbrs[1])
      }
      if (nt == length(nbrs)) {
        r <- which(pointsin == nbrs[length(nbrs)]) +  1
      }
      if ((nt > 0) & (nt < length(nbrs))) {
        r <- which(pointsin == nbrs[nt + 1])
      }
      if (is.null(schemehist) == FALSE) {
        if (schemehist[m - j + 1] == "Linear") {
          res <- LinearPred(pointsin, X, coeff, nbrs, remove, intercept = interhist[m - j + 1], neighbours)
        }
        if (schemehist[m - j + 1] == "Quad") {
          res <- QuadPred(pointsin, X, coeff, nbrs, remove, intercept = interhist[m - j + 1], neighbours)
        }
        if (schemehist[m - j + 1] == "Cubic") {
          res <- CubicPred(pointsin, X, coeff, nbrs, remove, intercept = interhist[m - j + 1], neighbours)
        }
      }else {
        #res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours)
        #res <- streamPred_ST(pointsin, pointsin_sptatiotemporal=pointsin_new, X_temporal=X_temporal, X_spatial=X_subind,  stream_network=example_network, coeff=coeff, pointsin_index=pointsin_new_index, nbrs, nbrs_total, nbrs_spatio_temporal, remove, adj=adjacency, stream_index, st_weight, method=c("IDW"), intercept, neighbours=NULL)
        #변경해야
        res <- streamPred_ST(pointsin, pointsin_sptatiotemporal=pointsin_new, X_temporal=X_temporal, X_spatial=X_subind,  
                             stream_network=example_network, coeff=coeff, pointsin_index=pointsin_new_index, nbrs, nbrs_total, 
                             nbrs_spatio_temporal, remove, adj=adjacency, stream_index, st_weight, method=c("IDW"), intercept, neighbours=NULL)
        
        #streamPred_ST <- function(pointsin, pointsin_sptatiotemporal=pointsin_new, X_temporal=X_temporal, X_spatial=X_subind,  stream_network=example_network, coeff, 
        #                          pointsin_index=pointsin_new_index, nbrs, nbrs_total, nbrs_spatio_temporal, remove, adj=adjacency, stream_index=c(5,3,4,8), 
        #                          st_weight=c(0.5,0.5), method=c("IDW","LC"), intercept=FALSE, neighbours=NULL)
          
      }
      if (length(res) == 2) {
        l <- res[[1]]
      }
      else {
        l <- res
      }
      #gamweights <- l$weights
      gamweights <- l$weight_extended
      #참고용
      #l1 <- PointsUpdate_ST(X, coeff, nbrs=nbrs_total, index=index_total, remove, pointsin=pointsin_new, weights, lengths, 
      #                      N=length(pointsin), weight_direct=res$weight_direct, nbrs_direct=res$nbrs_direct)
      #l1 <- UndoPointsUpdate(X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem)
      l1 <- UndoPointsUpdate_ST(X, coeff, nbrs=nbrs_total, index=index_total, remove, r, N, pointsin=pointsin_new, gamweights, lengths, lengthrem, weight_direct=l$weight_direct, nbrs_direct=l$nbrs_direct)
      coeff <- l1$coeff
      lengths <- l1$lengths
      pred <- sum(as.column(gamweights) * coeff[nbrs])
      coeff[remove] <- coeff[remove] + pred
      removelist <- setdiff(removelist, remove)
      if (r == 1) {
        lengths <- c(lengthrem, lengths)
        pointsin <- c(remove, pointsin)
        #추가
        pointsin_new <- sort(union(pointsin, pointsin_new))
      }
      if (r == (N + 1)) {
        lengths <- c(lengths, lengthrem)
        pointsin <- c(pointsin, remove)
        #추가
        pointsin_new <- sort(union(pointsin, pointsin_new))
      }
      if ((r > 1) & (r < (N + 1))) {
        lengths <- c(lengths[1:(r - 1)], lengthrem, lengths[r:N])
        pointsin <- c(pointsin[1:(r - 1)], remove, pointsin[r:N])
        #추가
        pointsin_new <- sort(union(pointsin, pointsin_new))
      }
    }
  }
  return(list(coeff = coeff, lengths = lengths, lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist))
}