#spatio-temporal version에 필요한 코드만을 모아놓는다
getnbrs_neighbor_temporal <- function(X_subind, X_remove, nbrs_temporal_index, type="long"){
  #이 함수는 아래 getnbrs_stream3으로 선택된 stream 이웃에 대해 각 장소별로 이웃을 찾아주는 함수다
  #X_remove: remove time 시점 (= X[remove])
  #type은 short과 long 두 가지
  #만약 long을 선택할 경우 특정 시간 간격에 있는 모든 측정시점을 반환
  #만약 short을 선택할 경우 관찰된 시점에 가장 가까운 장소를 1 또는 2개 반환
  nbrs_spatio_temporal_list <- list()
  length_stream3 <- length(X_subind)
  for(i in 1:length_stream3){
    if(type=="long"){
      nbrs_spatio_temporal_list[[i]] <- which(X_subind[[i]]>=nbrs_temporal_index[1] & X_subind[[i]]<=nbrs_temporal_index[2]) 
      #option 1: 특정 간격의 모든 관측시간 반환
    }else if(type=="short"){
      if(sum(X_subind[[1]]==X_remove)==1){
        #대응되는 점 하나만 넣는다
        nbrs_spatio_temporal_list[[i]] <- which(X_subind[[i]]==X_remove)
      }else{
        #이웃점 두개 (또는 한 개) 넣는다
        imsi_larger <- which(X_subind[[i]]>X_remove)[1]
        imsi_smaller <- which(X_subind[[i]]<X_remove)[length(which(X_subind[[i]]<X_remove))]
        nbrs_spatio_temporal_list[[i]] <- c(imsi_smaller, imsi_larger)
      }
    }
    nbrs_spatio_temporal_list[[i]] <-  nbrs_spatio_temporal_list[[i]][complete.cases(nbrs_spatio_temporal_list[[i]])]
  }
  return(nbrs_spatio_temporal_list)
}

#각 장소에 대해 이웃의 수를 출력하는 함수 작성
nbrs_spatial <- function(I_interval, X, X_subind){
  result_matrix <- matrix(0, nrow=length(X_subind), ncol=(length(I_interval)-1))
  result_dist_matrix <- matrix(0, nrow=length(X_subind), ncol=(length(I_interval)-1))
  for(i in 1:(length(I_interval)-1)){
    current_I <- I_interval[(i):(i+1)]
    result_matrix[,i] <- as.numeric(unlist(lapply(X_subind, function(x) length(which(current_I[1]<=x & current_I[2]>=x))  )))
    result_dist_matrix[,i] <- as.numeric(unlist(lapply(X_subind, function(x) sum(exp( -abs( X[i]- x[which(current_I[1]<=x & current_I[2]>=x)])  ))) ))
  }
  return(list(result=result_matrix, result_dist=result_dist_matrix))
}


LinearPred_stream <- function (pointsin, Xnbrs, Xremove, coeff, nbrs, remove, intercept, neighbours){
  #pointsin, neighbours: LinearPred에서는 안쓰임
  #coeff: 이웃 스트림쪽을 받는다
  #input값을 바로 받도록 변경 (우리는 여러 개의 time series를 다뤄야 하므로)
  #Xneighbours <- X[nbrs]
  
  #my addition:
  distvec <- mean(abs(Xnbrs-Xremove))
  
  Xneighbours <- Xnbrs
  Xneighbours <- as.column(Xneighbours)
  #Xremove <- X[remove]
  if (intercept) {
    Xneighbours <- cbind(1, Xneighbours)
    Xremove <- as.row(c(1, Xremove))
  }
  if (length(nbrs) >= 2) {
    temp <- crossprod(Xneighbours)
    mm <- Rmatsolve(temp) %*% t(Xneighbours)
    bhat <- mm %*% matrix(coeff[nbrs], ncol = 1)
    pred <- Xremove %*% bhat
    weights <- matrix(Xremove, nrow = 1) %*% mm
  }
  else {
    mm <- 0
    bhat <- 1
    weights <- 1
    pred <- coeff[nbrs]
  }
  
  #추가: distance를 고려하도록 한다
  
  return(list(weights = weights, pred = pred, coeff = coeff, dist=distvec))
}




#새로운 함수: streamPred_ST를 만든다
streamPred_ST <- function(pointsin, pointsin_sptatiotemporal=pointsin_new, X_temporal=X_temporal, X_spatial=X_subind,  stream_network=example_network, coeff, pointsin_index=pointsin_new_index, nbrs, nbrs_total, nbrs_spatio_temporal, remove, adj=adjacency, stream_index=c(5,3,4,8), st_weight=c(0.5,0.5), method=c("IDW","LC"), intercept=FALSE, neighbours=NULL){
  # pointsin_sptatiotemporal: 필요없음
  #pointsin: 리스트 형태로 받자(elt: spatial, temporal, spatio-temporal)
  #X_temporal
  #X_spatial
  #coeff: 해당 시공간 자료?
  #nbrs_temporal: temporal nbds
  #nbrs_spatial: spatial nbds
  #remove: remove인덱스(Q: spatial로 할 것인가, temporal로 할 것인가? 아니면 spatio-temporal로 할 것인가)
  #adj: stream flow reach에 대한 adjacency 정보다
  #st_weight: spatio-temporal weight를 나타내며 c(0.5 (spatial), 0.5 (tempral))을 default로 한다.
  #method: "IDW"와 "LC"가 있는데 크게 신경쓰지 않아도 되지 않을까?
  #intercept: TRUE & FALSE (어떻게 해야할지 고민해야 한다)
  
  #교수님께서 제안하신 방법: 주변 자료를 smoothing한다(또는 kriging을 한다)
  
  ##step 0: 기초작업
  #X_temporal <- X
  #X_spatial <- X_subind
  pointsin_temporal <- pointsin
  pointsin_spatial <- X_spatial
  #pointsin_sptatiotemporal <- pointsin_new
  stream_network <- example_network
  #coeff <- data_imsi #coeff <- data$normaldata[,c(5,3,4,8)]
  #st_weight <- c(0.999, 0.001); 
  s_weight <- st_weight[1]; t_weight <- st_weight[2]
  nbrs_temporal=nbrs_spatio_temporal$temporal
  nbrs_spatial=nbrs_spatio_temporal$spatial
  nbrs_spatioTemporal <- nbrs_spatio_temporal$spatioTemporal #이게 진짜 spatial 이웃?
  #순서대로 3번지역의 temporal 이웃, 4번지역의 temporal 이웃, 8번지역의 temporal 이웃이다
  #remove
  adj=adjacency #X_temporal 기반
  #method=c("IDW","LC")
  #intercept=FALSE
  #neighbours=NULL
  
  nbrs_combined <- c()
  X_nbrs_combined <- c()
  
  X_remove <- X_temporal[remove]
  
  ##step 1: compute prediction filter: that station #자기 자신의 prediction 필터로 계산
  #resultobj_direct <- LinearPred(pointsin=NULL, X=X, coeff=coeff[pointsin_index[1]:(pointsin_index[2]-1)], nbrs=nbrs_temporal, remove=remove, intercept=1, neighbours=NULL)
  #resultobj_direct <- LinearPred(pointsin=pointsin, X=X_temporal, coeff=coeff[pointsin_index[1]:(pointsin_index[2]-1)], nbrs=nbrs_total[c(1:length(nbrs_temporal))], remove=remove, intercept=TRUE, neighbours=NULL)
  resultobj_direct <- LinearPred(pointsin=pointsin, X=X_temporal, coeff=coeff[pointsin_index[1]:(pointsin_index[2]-1)], nbrs=nbrs_total[c(1:length(nbrs_temporal))], remove=remove, intercept=TRUE, neighbours=NULL)
  weight_direct <- resultobj_direct$weights
  pred_vec_direct <- resultobj_direct$pred
  
  iindex <- length(nbrs_temporal)
  
  ##step 2: compute prediction filter: other stream stations
  weight_imsi_list <- list()
  pred_vec <- c()
  dist_vec <- c()
  for(j in 1:length(nbrs_spatioTemporal)){
    #resultobj_nbrs <- LinearPred_stream(pointsin=NULL, Xnbrs=X_spatial[[j]][nbrs_spatioTemporal[[j]]], Xremove=X_remove, coeff=coeff[X_spatial[[j]],(j+1)], nbrs=nbrs_spatioTemporal[[j]], remove=remove, intercept=TRUE, neighbours=NULL)
    if(j==length(nbrs_spatio_temporal)){
      #resultobj_nbrs <- LinearPred_stream(pointsin=NULL, Xnbrs=X_spatial[[j]][nbrs_spatioTemporal[[j]]], Xremove=X_remove, coeff=coeff[pointsin_index[(j+1)]:length(coeff)], nbrs=nbrs_spatioTemporal[[j]], remove=remove, intercept=TRUE, neighbours=NULL)
      resultobj_nbrs <- LinearPred_stream(pointsin=NULL, Xnbrs=X_spatial[[j]][nbrs_spatioTemporal[[j]]], Xremove=X_remove, coeff=coeff, nbrs=nbrs_total[c((iindex+1):(iindex+length(nbrs_spatioTemporal[[j]])))], remove=remove, intercept=TRUE, neighbours=NULL)
    }else{
      #resultobj_nbrs <- LinearPred_stream(pointsin=NULL, Xnbrs=X_spatial[[j]][nbrs_spatioTemporal[[j]]], Xremove=X_remove, coeff=coeff[pointsin_index[(j+1)]:(pointsin_index[(j+2)]-1)], nbrs=nbrs_spatioTemporal[[j]], remove=remove, intercept=TRUE, neighbours=NULL)
      resultobj_nbrs <- LinearPred_stream(pointsin=NULL, Xnbrs=X_spatial[[j]][nbrs_spatioTemporal[[j]]], Xremove=X_remove, coeff=coeff, nbrs=nbrs_total[c((iindex+1):(iindex+length(nbrs_spatioTemporal[[j]])))], remove=remove, intercept=TRUE, neighbours=NULL)
    }
    iindex <- iindex + length(nbrs_spatioTemporal[[j]])
    weight_imsi_list[[j]] <- resultobj_nbrs$weights
    pred_vec <- c(pred_vec, resultobj_nbrs$pred)
    dist_vec <- c(dist_vec, resultobj_nbrs$dist)
  }
  
  
  ##step 3: compute shreve dist
  #compute_shreve: shreve 거리를 재도록 직접 제작한 함수
  shreve_order <- compute_shreve(adjacency)
  #위의 코드로 작성된 shreve_order를 바로 이용
  example_network_sub_pointdata <- example_network@obspoints@SSNPoints[[1]]@point.data[stream_index,]
  #shreve_order[example_network_sub_pointdata$rid] #8 5 1 9
  shreve_order[as.numeric(example_network_sub_pointdata$rid)]
  #spatial  자기자신 8/8, 나머지 이웃들 5/8, 1/8, 8/9
  spatial_weight <- rep(0, length(stream_index))
  spatial_weight <- apply(rbind(shreve_order[as.numeric(example_network_sub_pointdata$rid)] / shreve_order[as.numeric(example_network_sub_pointdata$rid[1])], shreve_order[as.numeric(example_network_sub_pointdata$rid[1])] / shreve_order[as.numeric(example_network_sub_pointdata$rid)] ), 2, min)
  #temporal 자기자신, 나머지 이웃들 0, 0, 0
  #X랑 X_temporal이랑 같다
  temporal_dist <- c(mean(abs(X_temporal[nbrs_spatio_temporal$temporal]-X_temporal[remove])), rep(0,(length(stream_index)-1)))
  
  ##step 4: compute spatio-temporal dist
  spatiotemporal_dist <- s_weight*(1/spatial_weight) + t_weight*temporal_dist
  final_weight <- (1/spatiotemporal_dist)/sum(1/spatiotemporal_dist)
  
  ##step 5: generate spatio-temporal prediction
  #각각의 point prediction을 우선 구해야
  pred_points <- c(pred_vec_direct, pred_vec)
  final_pred <- final_weight%*%pred_points
  
  ##extra: check the result
  #n.ind <- 5
  #plot(which(!is.na(data$normaldata[,n.ind])), as.numeric(data$normaldata[complete.cases(data$normaldata[,n.ind]),n.ind])-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,n.ind]),n.ind])), xlab="Time", ylim=c(-4,6), ylab="Nitrogen", main=mf04@obspoints@SSNPoints[[1]]@point.data$X[n.ind])
  #points(X[remove], coeff[X_temporal[remove],1], pch=16, col="red") 
  #points(X[remove],final_pred, pch=16, col="red") #-0.821199
  #points(X[remove],pred_points[1], pch=16, col="blue") ##-0.8329162
  
  distobj <- list()
  distobj$spatial <- (1/spatial_weight)
  distobj$temporal <- temporal_dist
  distobj$spatiotemporal <- spatiotemporal_dist
  
  final_weight_extended <- c(weight_direct)*final_weight[1]
  ind_imsi <- 2
  for(kk in 1:length(weight_imsi_list)){
    final_weight_extended <- c(final_weight_extended, weight_imsi_list[[kk]]*final_weight[ind_imsi])
    ind_imsi <- ind_imsi + 1
  }
  
  return(list(weights = final_weight, pred = final_pred, coeff = coeff, distobj=distobj, stream_index=stream_index, st_weight=st_weight, weight_extended=final_weight_extended, weight_direct=weight_direct, nbrs_direct=nbrs_temporal))
}




fwtnp_stream_ST <- function (data, example_network, stream_index, adjacency, st_weight=c(0.995,0.005), MyRivernetwork, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_ST, do.W = FALSE, varonly = FALSE){
  #data=data; example_network=example_network; stream_index=c(5,3,4,8); adjacency; st_weight=c(0.995,0.005); MyRivernetwork=MyRivernetwork; nkeep = 2; intercept = TRUE; initboundhandl = "reflect";  neighbours = 1; closest = FALSE; LocalPred = streamPred_ST; do.W = TRUE; varonly = FALSE
  n.ind <- stream_index[1] #stream_index의 첫 번째 원소를 제거하려는 시계열로
  n.subind <- stream_index[-1] #stream_index의 첫 번째 원소를 제외한 나머지 원소들의 시계열을 이웃들의 시계열로 정한다
  X_subind <- list()
  for(i in 1:length(n.subind)){
    X_subind[[i]] <- which(!is.na(data$normaldata[,n.subind[i]]))
  }
  #X_subind[[1]] <- which(!is.na(data$normaldata[,3])); X_subind[[2]] <- which(!is.na(data$normaldata[,4])); X_subind[[3]] <- which(!is.na(data$normaldata[,8]))
  #만약 그냥 temporal version만 고려한다면...
  I <- intervals(X=which(!is.na(data$normaldata[,n.ind])), initboundhandl = "reflect")
  lengths <- lengthintervals(X=which(!is.na(data$normaldata[,n.ind])), I, type = "midpoints", neighbours=2,  closest=TRUE)
  for(i in 1:length(n.subind)){
    I_imsi <- intervals(X=which(!is.na(data$normaldata[,n.subind[i]])), initboundhandl = "reflect")
    lengths <- c(lengths, lengthintervals(X=which(!is.na(data$normaldata[,n.subind[i]])), I_imsi, type = "midpoints", neighbours=2,  closest=TRUE)) #compute extended lengths
  }
  
  #nbrs_spatial_result: 만들어놓은 함수 이용 (굳이 필요한지는 잘 모르겠음)
  nbrs_spatial_result <- nbrs_spatial(I_interval = I, X=which(!is.na(data$normaldata[,n.ind])), X_subind=X_subind )
  colSums(nbrs_spatial_result$result_dist)
  
  
  #I <- initint2_stream(example_network, adjacency, linear=FALSE, MyRivernetwork)$I
  I_lines <- initint2_stream(example_network, adjacency, linear=FALSE, MyRivernetwork)$I_lines
  
  X <- which(!is.na(data$normaldata[,n.ind])) #X: 70개 (짧다)
  #f <- data$normaldata[,n.ind]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,n.ind]),n.ind])) #f: 70개 (짧다)
  f <- data$normaldata[complete.cases(data$normaldata[,n.ind]),n.ind]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,n.ind]),n.ind])) #점검 필요
  #f: 나중에 쓰지 않기 때문에 큰 상관 없다
  n <- length(X) #70
  
  #nkeep <- 2
  
  origlengths <- lengths
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
  #추가로 반환할 리스트 작성
  nbrs_total_list <- list()
  nbrs_spatio_temporal_return_list <- list()
  #추가로 반환할 리스트 작성 (끝)
  pointsin <- matrix(1:n, 1, n)
  
  #새로운 pointsin 작성할 필요가 있음
  n_new <- 0
  pointsin_new_index <- c(1)
  for(i in 1:length(stream_index)){
    length_imsi <- length(which(!is.na(data$normaldata[,stream_index[i]])))
    n_new <- n_new + length_imsi 
    pointsin_new_index <- c(pointsin_new_index, length_imsi+pointsin_new_index[length(pointsin_new_index)]) #이어붙일 시계열들의 시작점 반환
  }
  pointsin_new_index <- pointsin_new_index[-length(pointsin_new_index)]
  #n_new <- length(which(!is.na(data$normaldata[,5]))) + length(which(!is.na(data$normaldata[,3]))) + length(which(!is.na(data$normaldata[,4]))) + length(which(!is.na(data$normaldata[,8])))
  pointsin_new <- matrix(1:n_new, 1, n_new)
  
  #pointsin_new_index <- c(1, (length(which(!is.na(data$normaldata[,5])))  + 1), (length(which(!is.na(data$normaldata[,5]))) + length(which(!is.na(data$normaldata[,3]))) +1 ), (length(which(!is.na(data$normaldata[,5]))) + length(which(!is.na(data$normaldata[,3]))) + length(which(!is.na(data$normaldata[,4]))) + 1))
  
  #pointsin <- pointsin[order(X)]
  X_temporal <- X
  
  X <- c()
  coeff <- c()
  for(i in 1:length(stream_index)){
    X <- c(X, which(!is.na(data$normaldata[,stream_index[i]]) ))
    coeff <- c(coeff, as.numeric(data$normaldata[complete.cases(data$normaldata[,stream_index[i]]),stream_index[i]]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,stream_index[i]]),stream_index[i]]))))
  }
  
  #X <- c(which(!is.na(data$normaldata[,5])), which(!is.na(data$normaldata[,3])), which(!is.na(data$normaldata[,4])), which(!is.na(data$normaldata[,8])))
  X <- as.row(X)
  #coeff <- f #변경해야
  #coeff <- c(as.numeric(data$normaldata[complete.cases(data$normaldata[,5]),5]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,5]),5]))), 
  #           as.numeric(data$normaldata[complete.cases(data$normaldata[,3]),3]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,3]),3]))), 
  #           as.numeric(data$normaldata[complete.cases(data$normaldata[,4]),4]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,4]),4]))), 
  #           as.numeric(data$normaldata[complete.cases(data$normaldata[,8]),8]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,8]),8]))))
  #(나중에 추가할 것) detrending 시 linear trend 또한 고려해보자
  coeff <- as.row(coeff)
  matno <- n - nkeep #n-nkeep
  
  #변경작업: W matrix의 행과 열을 temporal / spatial 이웃의 수를 모두 고려해 지정한다
  #do.W=T; varonly=F
  W <- v <- NULL
  if ((do.W == 1) & (varonly == 1)) {
    varonly <- FALSE
  }
  ex <- do.W + varonly
  if (ex == 1) {
    W <- diag(n_new)
  }
  if (varonly) {
    v <- rep(1, times = n_new)
  }
  lengths_sub <- lengths[c(1:length(X_temporal))] #70개짜리 저장
  
  #(a) spatial 이웃 찾기 (spatial version에서 만들어놨던 함수 이용)
  nbrs_spatial <- getnbrs_stream3(remove=n.ind, pointsin=c(1:nrow(example_network@obspoints@SSNPoints[[1]]@point.coords)), example_network, I_lines, adjacency, upperExtra=FALSE)
  #remove=n.ind; pointsin=c(1:nrow(example_network@obspoints@SSNPoints[[1]]@point.coords)); example_network; I_lines; adjacency; upperExtra=FALSE
  
  for (j in 1:matno) {
    #j <- 1
    #remove_candidate <- which(lengths_sub==min(lengths_sub))
    #
    #if(length(remove_candidate)!=1){
    #  #select one among candidates
    #  temp_dist_compared <- rep(0, length(remove_candidate))
    #  for(k in 1:length(remove_candidate)){
    #    nbrs_spatio_temporal_list <- getnbrs_neighbor_temporal(X_subind, X_remove=X_temporal[remove_candidate[k]], nbrs_temporal_index, type="short")
    #    for(l in 1:length(nbrs_spatio_temporal_list)){
    #      temp_dist_compared[k] <- sum(abs(X_subind[[l]][nbrs_spatio_temporal_list[[l]]]-X_temporal[remove]))
    #    }
    #  }
    #  remove_candidate <- remove_candidate[which.min(temp_dist_compared)]
    #  remove <- pointsin[remove_candidate]
    #  removelist[j] <- remove
    #}else{
      #한개일 경우에는 그냥 진행
      remove_candidate <- order(lengths_sub)[1]
      remove <- pointsin[remove_candidate]
      removelist[j] <- remove
    #}
    #spatio-temporal 이웃 찾는 함수 만들기
    #(a) spatial 이웃 찾기 (spatial version에서 만들어놨던 함수 이용) -> for loop 밖으로 빼냄
    #(b) spatial 이웃들의 모든 temporal 이웃 반환하기 (X_subind에 time series의 위치가 들어 있다)
    nbrs_temporal_old <- getnbrs(X_temporal, remove=remove, pointsin=pointsin, neighbours, closest)
    nbrs_temporal_index <- X_temporal[nbrs_temporal_old$nbrs]
    nbrs_spatio_temporal_list <- list()
    #nbrs_spatio_temporal_list[[1]] <- X_subind[[1]][which(X_subind[[1]]>=nbrs_temporal_index[1] & X_subind[[1]]<=nbrs_temporal_index[2])]
    #nbrs_spatio_temporal_list[[2]] <- X_subind[[2]][which(X_subind[[2]]>=nbrs_temporal_index[1] & X_subind[[2]]<=nbrs_temporal_index[2])]
    #nbrs_spatio_temporal_list[[3]] <- X_subind[[3]][which(X_subind[[3]]>=nbrs_temporal_index[1] & X_subind[[3]]<=nbrs_temporal_index[2])]
    
    #for(k in 1:length(nbrs_spatial$nbrs)){
    #  if(length(nbrs_temporal_index)>=2){
    #    nbrs_spatio_temporal_list[[k]] <- which(X_subind[[k]]>=nbrs_temporal_index[1] & X_subind[[k]]<=nbrs_temporal_index[2]) #option 1: 특정 간격의 모든 관측시간 반환
    #  }else if(length(nbrs_temporal_index)==1){
    #    if(nbrs_temporal_index < X_temporal[remove]){
    #      nbrs_spatio_temporal_list[[k]] <- which(X_subind[[k]]>=nbrs_temporal_index[1] & X_subind[[k]]<=X_temporal[remove]) #option 1: 특정 간격의 모든 관측시간 반환
    #    }else{
    #      nbrs_spatio_temporal_list[[k]] <- which(X_subind[[k]]<=nbrs_temporal_index[1] & X_subind[[k]]>=X_temporal[remove]) #option 1: 특정 간격의 모든 관측시간 반환
    #    }
    #  }
    #}
    #nbrs_spatio_temporal_list[[1]] <- which(X_subind[[1]]>=nbrs_temporal_index[1] & X_subind[[1]]<=nbrs_temporal_index[2]) #option 1: 특정 간격의 모든 관측시간 반환
    #nbrs_spatio_temporal_list[[2]] <- which(X_subind[[2]]>=nbrs_temporal_index[1] & X_subind[[2]]<=nbrs_temporal_index[2]) #option 1: 특정 간격의 모든 관측시간 반환
    #nbrs_spatio_temporal_list[[3]] <- which(X_subind[[3]]>=nbrs_temporal_index[1] & X_subind[[3]]<=nbrs_temporal_index[2]) #option 1: 특정 간격의 모든 관측시간 반환
    
    nbrs_spatio_temporal_list <- getnbrs_neighbor_temporal(X_subind, X_remove=X_temporal[remove], nbrs_temporal_index, type="short")
    
    nbrs_spatio_temporal <- list()
    nbrs_spatio_temporal$temporal <- nbrs_temporal_old$nbrs
    nbrs_spatio_temporal$spatial <- nbrs_spatial
    nbrs_spatio_temporal$spatioTemporal <- nbrs_spatio_temporal_list
    
    index <- nbrs_temporal_old$index #index_total 함수도 만들어야
    #pointsin_imsi_index: 이 함수 나중에 수정해야
    pointsin_imsi_index <- c(1, (length(pointsin) + 1))
    for(k in 1:length(X_subind)){
      if(k!=length(X_subind)){
        pointsin_imsi_index <- c(pointsin_imsi_index, pointsin_imsi_index[length(pointsin_imsi_index)] + length(which(!is.na(data$normaldata[,n.subind[k]]))))
      }
    }
    #pointsin_imsi_index <- c(1, (length(pointsin) + 1), (length(pointsin) + length(which(!is.na(data$normaldata[,3]))) +1 ), (length(pointsin) + length(which(!is.na(data$normaldata[,3]))) + length(which(!is.na(data$normaldata[,4]))) + 1))
    
    #nbrs_total: 위의 nbrs 함수들은 1번부터 시작하는데 이것을 연속해서 연결해 줄 필요 있음
    nbrs_total <- c(nbrs_spatio_temporal$temporal)
    for(k in 1:length(nbrs_spatio_temporal_list)){
      nbrs_total <- c(nbrs_total, pointsin_new_index[(k+1)]+nbrs_spatio_temporal_list[[k]]-1)
      index <- c(index, pointsin_imsi_index[(k+1)]+nbrs_spatio_temporal_list[[k]]-1 )
    }
    #nbrs_total <- c(nbrs_spatio_temporal$temporal, pointsin_new_index[2]+nbrs_spatio_temporal_list[[1]]-1, pointsin_new_index[3]+nbrs_spatio_temporal_list[[2]]-1, pointsin_new_index[4]+nbrs_spatio_temporal_list[[3]]-1)
    
    #참고: LocaPred 함수의 구성
    #res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours)
    
    #data_imsi2 <- cbind(data$normaldata[,5]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,5]),5])), data$normaldata[,3]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,3]),3])), data$normaldata[,4]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,4]),4])), data$normaldata[,8]-mean(as.numeric(data$normaldata[complete.cases(data$normaldata[,8]),8])))
    
    #res <- streamPred_ST(pointsin, pointsin_sptatiotemporal=pointsin_new, X_temporal=X_temporal, X_spatial=X_subind,  stream_network=example_network, coeff=as.matrix(data_imsi2), pointsin_index=pointsin_new_index, nbrs, nbrs_total, remove, adj=adjacency, stream_index=c(5,3,4,8), st_weight=c(0.995,0.005), method=c("IDW"), intercept=TRUE, neighbours=NULL)
    res <- streamPred_ST(pointsin, pointsin_sptatiotemporal=pointsin_new, X_temporal=X_temporal, X_spatial=X_subind,  stream_network=example_network, coeff=coeff, pointsin_index=pointsin_new_index, nbrs=nbrs_temporal_old, nbrs_total, nbrs_spatio_temporal, remove, adj=adjacency, stream_index, st_weight, method=c("IDW"), intercept, neighbours=NULL)
    
    #streamPred_ST: source_ST에 작성함
    #pointsin;pointsin_sptatiotemporal=pointsin_new; X_temporal=X; X_spatial=X_subind;  stream_network=example_network; coeff=coeff; pointsin_index=pointsin_new_index; nbrs=nbrs_temporal_old; adj=adjacency; stream_index=c(3,4,1,2); st_weight=c(0.995,0.005); method="IDW"; intercept=TRUE; neighbours=NULL
    
    #작성해야
    if (length(res) == 2) {
      l <- res[[1]]
      clolist[j] <- res[[2]][[1]]
      nbrs <- res[[2]][[2]]
      index <- res[[2]][[3]]
    }else {
      l <- res
    }
    #neighbrs[[j]] <- nbrs #고쳐야
    #weights <- l[[1]] #고쳐야
    neighbrs[[j]] <- nbrs_total #고쳐야
    weights <- res$weight_extended #고쳐야
    pred <- l[[2]]
    #if (length(l) == 3) {
    scheme <- NULL
    int <- NULL
    details <- NULL
    #}else {
    #  scheme <- l[[5]]
    #  int <- l[[4]]
    #  details <- l[[6]]
    #}
    
    #coeff를 새로 바꾸자
    coeff[remove] <- coeff[remove] - pred
    
    index_total <- match(nbrs_total, pointsin_new)
    #lengths_total 함수도 만들어야
    print(j)
    #print(remove)
    #print(nbrs_total)
    #print(index_total)
    #print(weights)
    print(removelist)
    #print(res$weight_direct)
    
    #print(length(lengths))
    #print(lengths)
    #l1 <- PointsUpdate(X, coeff, nbrs=nbrs_total, index=index_total, remove, pointsin=pointsin_new, weights, lengths)
    l1 <- PointsUpdate_ST(X, coeff, nbrs=nbrs_total, index=index_total, remove, pointsin=pointsin_new, weights, lengths, N=length(pointsin), weight_direct=res$weight_direct, nbrs_direct=res$nbrs_direct)
    coeff <- l1$coeff
    lengths <- l1$lengths
    r <- l1$r
    weights <- l1$weights
    print(weights)
    print(nbrs_total)
    N <- l1$N
    alpha <- l1$alpha
    if (ex) {
      if (varonly) {
        #W[r, ] <- W[r, ] - colSums(as.vector(weights) * matrix(W[index, ], nrow = length(nbrs_total)))
        W[r, ] <- W[r, ] - colSums(as.vector(weights) * matrix(W[index_total, ], nrow = length(nbrs_total)))
        #W[index, ] <- W[index, ] + matrix(alpha) %*%  W[r, ]
        W[index_total, ] <- W[index_total, ] + matrix(alpha) %*%  W[r, ]
        v[remove] <- sum(W[r, ]^2)
        #np <- setdiff(1:length(pointsin), r)
        np <- setdiff(1:length(pointsin_new), r)
        W <- W[np, ]
      } else {
        print(nbrs_total)
        print(weights)
        W[remove, ] <- W[remove, ] - colSums(as.vector(weights) * matrix(W[nbrs_total, ], nrow = length(nbrs_total)))
        W[nbrs_total, ] <- W[nbrs_total, ] + matrix(alpha) %*% W[remove, ]
      }
    }
    lengthsremove[j] <- lengths[r]
    gamlist[[j]] <- weights
    alphalist[[j]] <- alpha
    schemehist[j] <- scheme
    interhist[j] <- int
    #추가로 반환할 리스트 작성
    nbrs_total_list[[j]] <- nbrs_total
    nbrs_spatio_temporal_return_list[[j]] <- append(list(), nbrs_spatio_temporal)
    #추가로 반환할 리스트 작성(끝)
    lengths_sub <- lengths[setdiff(1:length(pointsin), r)]
    lengths <- lengths[setdiff(1:length(pointsin_new), r)]
    #lengths_sub <- lengths_sub[setdiff(1:length(pointsin), r)]
    pointsin <- setdiff(pointsin, remove)
    #pointsin_new 또한 업데이트 시켜야
    pointsin_new <- setdiff(pointsin_new, remove)
  }#end of for loop
  if (varonly) {
    #v[pointsin] <- rowSums(W^2)
    v[pointsin_new] <- rowSums(W^2)
    W <- NULL
  }
  N <- length(pointsin)
  return(list(x=NULL, coeff = coeff, lengths = lengths, 
              lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist, 
              neighbrs = neighbrs, schemehist = schemehist, interhist = interhist, 
              clolist = clolist, gamlist = gamlist, alphalist = alphalist, nbrs_total_list=nbrs_total_list, nbrs_spatio_temporal_return_list=nbrs_spatio_temporal_return_list,
              W = W, v = v, data=data, example_network=example_network, stream_index=stream_index, adjacency=adjacency, st_weight=st_weight, MyRivernetwork=MyRivernetwork, pointsin_new=pointsin_new, pointsin_new_index=pointsin_new_index))
}

PointsUpdate_ST <- function (X, coeff, nbrs, index, remove, pointsin, weights, lengths, N, weight_direct=NULL, nbrs_direct=NULL){
  r <- which(pointsin == remove)
  
  index_direct <- match(nbrs_direct, pointsin)
  #N <- length(pointsin)
  if ((r >= 2) & (r <= (N - 1))) {
    #lengths[index] <- as.row(lengths[index])
    lengths[index_direct] <- as.row(lengths[index_direct])
    weight_direct <- as.row(weight_direct)
    lengths[index_direct] <- lengths[index_direct] + lengths[r] * weight_direct
  }else {
    if (r == 1) {
      lengths[2] <- lengths[2] + lengths[1]
    }
    if (r == N) {
      lengths[N - 1] <- lengths[N - 1] + lengths[N]
    }
  }
  alpha <- matrix(0, 1, length(nbrs))
  if (length(nbrs) >= 2) {
    alpha <- lengths[r] * lengths[index]/(sum(lengths[index]^2))
    coeff[pointsin[index]] <- coeff[pointsin[index]] + alpha * coeff[remove]
  } else {
    q <- which(pointsin == nbrs)
    alpha <- lengths[r]/lengths[q]
    coeff[pointsin[q]] <- coeff[pointsin[q]] + alpha * coeff[remove]
  }
  return(list(coeff = coeff, lengths = lengths, r = r, N = N, weights = as.row(weights), alpha = alpha))
}

UndoPointsUpdate_ST <- function (X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, weight_direct=NULL, nbrs_direct=NULL) {
  #weight_direct, nbrs_direct를 넣어야 할 지 좀 더 고민해봐야 함
  index_direct <- match(nbrs_direct, pointsin)
  
  alpha <- matrix(0, 1, length(nbrs))
  if ((r > 1) & (r <= N)) {
    #alpha <- lengths[index] * lengthrem/(sum(lengths[index]^2))
    alpha <- lengths[index_direct] * lengthrem/(sum(lengths[index_direct]^2))
    #coeff[nbrs] <- coeff[nbrs] - alpha * coeff[remove]
    coeff[nbrs_direct] <- coeff[nbrs_direct] - alpha * coeff[remove]
    #lengths[index] <- as.row(lengths[index])
    lengths[index_direct] <- as.row(lengths[index_direct])
    #prod <- gamweights * lengthrem
    prod <- weight_direct * lengthrem
    prod <- as.row(prod)
    #lengths[index] <- lengths[index] - prod
    lengths[index_direct] <- lengths[index_direct] - prod
  }
  if ((r == 1) | (r == (N + 1))) {
    #q <- which(pointsin == nbrs)
    q <- which(pointsin == nbrs_direct)
    alpha <- lengthrem/lengths[q]
    coeff[pointsin[q]] <- coeff[pointsin[q]] - alpha * coeff[remove]
    lengths[q] <- lengths[q] - lengthrem
  }
  return(list(coeff = coeff, lengths = lengths, alpha = alpha))
}

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
  X <- which(!is.na(data$normaldata[,n.ind]))
  #my addition
  n.ind <- stream_index[1] #stream_index의 첫 번째 원소를 제거하려는 시계열로
  n.subind <- stream_index[-1] #stream_index의 첫 번째 원소를 제외한 나머지 원소들의 시계열을 이웃들의 시계열로 정한다
  X_subind <- list()
  for(i in 1:length(n.subind)){
    X_subind[[i]] <- which(!is.na(data$normaldata[,n.subind[i]]))
  }
  #end of my addition
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
      #pred <- sum(as.column(gamweights) * coeff[nbrs])
      pred <- sum(as.column(gamweights) * coeff[nbrs_total])
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