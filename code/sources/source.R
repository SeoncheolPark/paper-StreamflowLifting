library(stringr)
library(SSN)
library(rgdal) #for readOGR function
library(riverdist) #for computing upstream distance
library(rivernet)
library(xts)
library(shp2graph) #get_binaryIDs_stream 함수에 사용
library(spam)
library(scales)
########################################
##CODES related to the SSN and smnet package
########################################
initint2_stream <- function(example_network, adjacency, linear=TRUE, MyRivernetwork=NULL){
  #example_network: SSN package의 "SpatialStreamNetwork"를 따라야 하는듯
  if(linear==FALSE & is.null(MyRivernetwork)==TRUE){
    stop("For nonlinear stream network, riverdist object is needed.")
  }
  #example_network: "SpatialStreamNetwork" obj
  #adjacency: list 형태
  #linear: 주어진 stream network가 간단한 직선 형태인지 아닌지 나타내는 indicator
  setClass("ILines", slots=list(Lines="list", ID="character", weight="numeric", range="numeric"))
  
  #추가: 각 segment들이 어떤 node에 배속되었는지 체크하여 결과값을 주면 좋다
  I_lines <- list()
  I <- matrix(rep(0, length(unique(example_network@obspoints@SSNPoints[[1]]@point.data$pid))), nrow=1)
  colnames(I) <- unique(example_network@obspoints@SSNPoints[[1]]@point.data$pid)
  for(j in 1:length(unique(example_network@obspoints@SSNPoints[[1]]@point.data$pid))){
    I_lines[[j]] <- list()
  }
  
  network_SegmentID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID #observation station의 Segment ID
  for(i in 1:length(unique(example_network@network.line.coords$SegmentID))){
    #function matching each segment to the nearest point according to stream network structure
    n_inthe_Segment <- sum(i==network_SegmentID) #i번제 segment에 관찰값이 있는지 확인해보기 위함
    if(n_inthe_Segment==0){
      #segment 안에 point 들이 하나도 없는 경우
      #downstream을 따라 가장 가까운  stream들 중 obs.point가 있는 스트림을 찾는다
      ind_imsi <- i
      #pt_downstream_dist <- 0
      while(TRUE){
        #해당 지류의 하류를 찾는다
        mouth_index_imsi <- sum(as.matrix(adjacency$adjacency)[ind_imsi,]==1) #하류 지류의 갯수
        if(mouth_index_imsi==0){
          #이경우는 강물을 따라 내려가는 동안 아무러 관측지점을 못만난다는 얘기인데, 고려해봐야 할 듯
          break
        }else if(mouth_index_imsi==1){
          #하류 지류의 갯수가 2개는 아니라고 가정
          mouth_seg_imsi <- which(as.matrix(adjacency$adjacency)[ind_imsi,]==1)
          if(sum(mouth_seg_imsi==network_SegmentID)==0){
            #이 경우에는 다음 segment로 넘어가야 한다
            ind_imsi <- mouth_seg_imsi
            #pt_downstream_dist
          }else{
            #이 경우는 해당 segment에 관측장소가 여러 개 있을 경우이다
            #이제 해당 segment안에 가장 가까운 인덱스를 찾는다
            if(sum(mouth_seg_imsi==network_SegmentID)==1){
              #해당 segment에 관측장소가 단 한개만 있을 때
              ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi)
              I[ind] <- I[ind] + example_network@data$Length[i]
              I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range=c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
              #I_lines[[ind]] <- c(I_lines[[ind]], example_network@lines[[i]])
              I_lines[[ind]] <- c(I_lines[[ind]], I_lines_new)
            }else{
              #해당 segment에 관측장소가 여러 개 있을 때
              ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi) #해당 segment들에 들어있는 obs들의 좌표들
              #그러면 해당 segment와 가장 가까운 장소가 어디인지 판단해야 한다
              #example_network@lines에서 들어있는 Line들은 상류 -> 하류 순으로 되어있다
              pt_upstream <- example_network@lines[[mouth_seg_imsi]]@Lines[[1]]@coords[1,]
              nearest_pt_index <- ind[which.min(as.matrix(dist(rbind(pt_upstream,example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]) ))[-1,1])]
              I[nearest_pt_index] <- I[nearest_pt_index] + example_network@data$Length[i]
              I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range=c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
              #I_lines[[nearest_pt_index]] <- c(I_lines[[nearest_pt_index]], example_network@lines[[i]])
              I_lines[[nearest_pt_index]] <- c(I_lines[[nearest_pt_index]], I_lines_new)
            }
            break
          }
        }
      }
    }else if(n_inthe_Segment==1){
      #segment 안에 point 들이 한 개 있는 경우
      #그냥 더해주면 됨
      ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == i)
      I[ind] <- I[ind] + example_network@data$Length[i]
      I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range= c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
      #I_lines[[ind]] <- c(I_lines[[ind]], example_network@lines[[i]])
      I_lines[[ind]] <- c(I_lines[[ind]], I_lines_new)
    }else{
      #segment 안에 point 들이 여러 개 있는 경우
      #point의 좌표를 가지고 등분할 하면 된다
      ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == i)
      
      #example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]
      node_order_index <- order(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[ind], decreasing = TRUE)
      
      if(linear==TRUE){
        #선형인 경우 (비선형인 경우는 따로 만들어야), upDist 이용
        ordered_coords <- example_network@obspoints@SSNPoints[[1]]@point.coords[ind[node_order_index],]
        
        middle_pts <- c() #중점을 저장할 벡터
        for(j in 1: (length(ind)-1)){
          #중점을 계산한다
          middle_pts_imsi <- c((ordered_coords[j,1]+ordered_coords[(j+1),1])/2, (ordered_coords[j,2]+ordered_coords[(j+1),2])/2)
          middle_pts <- rbind(middle_pts, middle_pts_imsi)
        }
        #end_pts <- example_network@lines[[i]]@Lines[[1]]@coords[order(example_network@lines[[i]]@Lines[[1]]@coords[,1]),]
        end_pts <- example_network@lines[[i]]@Lines[[1]]@coords[,]
        middle_pts <- rbind(end_pts[1,], middle_pts, end_pts[2,])
        
        #양 끝점과 중점으로 길이를 계산해 그만큼씩 저장해둔다
        for(k in 1:length(ind)){
          I[ind[node_order_index[k]]] <- I[ind[node_order_index[k]]] + sqrt(sum((middle_pts[k,]-middle_pts[(k+1),])^2))
          example_network_lines_imsi <- example_network@lines[[i]]
          example_network_lines_imsi@ID <- paste(example_network_lines_imsi@ID, "-", k, sep="")
          example_network_lines_imsi@Lines[[1]]@coords <- rbind(middle_pts[k,], middle_pts[(k+1),])
          I_lines_new <- new("ILines", Lines=example_network_lines_imsi@Lines, ID=example_network_lines_imsi@ID, weight=1, range=c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
          #I_lines[[ind[node_order_index[k]]]] <- c(I_lines[[ind[node_order_index[k]]]], example_network_lines_imsi)
          I_lines[[ind[node_order_index[k]]]] <- c(I_lines[[ind[node_order_index[k]]]], I_lines_new)
        }
      }else if(linear==FALSE){
        #주어진 stream set이 linear하지 않은 경우
        #각 segment 안에서 한 개의 토막 당 거리가 같다고 가정을 하고 토막 기준으로 등분할 한다
        
        
        #example_network@obspoints@SSNPoints[[1]]@point.data$rid[ind]
        #example_network@obspoints@SSNPoints[[1]]@point.data$segid[ind]
        
        dist_combined <- rep(0, nrow(example_network@lines[[i]]@Lines[[1]]@coords))
        for(j in 1:nrow(example_network@lines[[i]]@Lines[[1]]@coords)){
          dist_imsi_vec <- rep(0, length(ind))
          for(k in 1:length(ind)){
            dist_imsi_vec[k] <- riverdistance(startseg = i, endseg = i, startvert = min(example_network@obspoints@SSNPoints[[1]]@point.data$segid[ind[k]], j), endvert = max(example_network@obspoints@SSNPoints[[1]]@point.data$segid[ind[k]], j), rivers= MyRivernetwork)
          }
          dist_combined[j] <- ind[which.min(dist_imsi_vec)]
        }
        
        
        
        for(l in 1:length(ind)){
          range_imsi <- range(which(dist_combined==ind[l]))
          I[ind[l]] <- I[ind[l]] + riverdistance(startseg = i, endseg = i, startvert = range_imsi[1], endvert = range_imsi[2], rivers= MyRivernetwork)
          if(range_imsi[1]!=1){
            I[ind[l]] <- I[ind[l]] + riverdistance(startseg = i, endseg = i, startvert = (range_imsi[1]-1), endvert = range_imsi[1], rivers= MyRivernetwork)/2
          }
          if(range_imsi[2]!=nrow(example_network@lines[[i]]@Lines[[1]]@coords)){
            I[ind[l]] <- I[ind[l]] + riverdistance(startseg = i, endseg = i, startvert = range_imsi[2], endvert = (range_imsi[2]+1), rivers= MyRivernetwork)/2
          }
          I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range=range(which(dist_combined==ind[l])))
          #I_lines[[ind[node_order_index[k]]]] <- c(I_lines[[ind[node_order_index[k]]]], example_network_lines_imsi)
          I_lines[[ind[l]]] <- c(I_lines[[ind[l]]], I_lines_new)
        }
      }
      
      
    }
  }
  return(list(I=I, I_lines=I_lines))
  #I: area of node i
  #I_lines: shows that each node i has which segement
}

getnbrs_stream3 <- function(remove, pointsin, example_network, I_lines, adjacency, upperExtra=FALSE){
  #프로토타입 3
  r = which(remove==pointsin)

  nbrs_upstream <- c()
  nbrs_downstream <- c()
  
  #나중에 pointsin 추가해야 할지도
  network_SegmentID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[r]
  n_inthe_Segment <- sum(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID==network_SegmentID)
  r_length <- nchar(adjacency$rid_bid[network_SegmentID,2])
  
  #index <- setdiff()
  
  network_indexID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[-r]
  
  #samestream_index: upstream인지 downstream인지 upDist를 보고 결정해야
  samestream_index <- setdiff(which(adjacency$rid_bid[network_indexID,2]==adjacency$rid_bid[network_SegmentID,2]), remove)
  #upstream_index: 상류
  upstream_index <- setdiff(str_which( substr(adjacency$rid_bid[network_indexID,2], start=1, stop=r_length), adjacency$rid_bid[network_SegmentID,2]), remove)
  upstream_index <- setdiff(upstream_index, samestream_index)
  
  if(length(samestream_index)!=0){
    samestream_index_imsi <- match(samestream_index, pointsin)
    for(i in 1:length(samestream_index_imsi)){
      if(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[r] < example_network@obspoints@SSNPoints[[1]]@point.data$upDist[samestream_index_imsi[i]]){
        #상류에 이웃 추가
        nbrs_upstream <- c(nbrs_upstream, samestream_index[i])
      }else if(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[r] > example_network@obspoints@SSNPoints[[1]]@point.data$upDist[samestream_index_imsi[i]]){
        #하류에 이웃 추가
        nbrs_downstream <- c(nbrs_downstream, samestream_index[i])
      }
    }
  }
  
  if(length(nbrs_upstream)>=2){
    
  }
  if(length(nbrs_downstream)>=2){
    
  }
  
  if(length(nbrs_upstream)==0 & length(upstream_index)!=0){
    #만약 upstream 이웃이 하나도 없으면 upstream 이웃을 찾는다
    cand_length <- nchar(adjacency$rid_bid[network_indexID,2])[-remove]
    nbrs_cand <- pointsin[nchar(adjacency$rid_bid[network_indexID,2]) > nchar(adjacency$rid_bid[network_SegmentID,2])]
    nbrs_cand <- setdiff(nbrs_cand,remove)
    nbrs_cand <- intersect(nbrs_cand, upstream_index)
    while(TRUE){
      diff_vec <- nchar(adjacency$rid_bid[network_indexID[nbrs_cand],2])- nchar(adjacency$rid_bid[network_SegmentID,2])
      upstream_candidate <- nbrs_cand[which(diff_vec==min(diff_vec))]
      if(length(upstream_candidate)>1){
        upstream_candidate <- upstream_candidate[which.min(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[match(upstream_candidate, pointsin)])]
      }
      
      nbrs_upstream <- c(nbrs_upstream, upstream_candidate)
      nchar_upstream_cand <- nchar(adjacency$rid_bid[network_indexID[upstream_candidate],2])
      up_upstream <- str_which( substr(adjacency$rid_bid[network_indexID[nbrs_cand],2], start=1, stop=nchar_upstream_cand ), adjacency$rid_bid[network_indexID[upstream_candidate],2])
      if(length(up_upstream)!=0){
        nbrs_cand <- nbrs_cand[-up_upstream]
      }
      if(length(nbrs_cand)==0){
        break
      }
    }
  }
  if(length(nbrs_downstream)==0){
    #finding downstream neighbour
    #segment 안에 point 들이 하나도 없는 경우
    #downstream을 따라 가장 가까운  stream들 중 obs.point가 있는 스트림을 찾는다
    ind_imsi <- network_SegmentID 
    #pt_downstream_dist <- 0
    while(TRUE){
      #해당 지류의 하류를 찾는다
      mouth_index_imsi <- sum(as.matrix(adjacency$adjacency)[ind_imsi,]==1) #하류 지류의 갯수
      if(mouth_index_imsi==0){
        #이경우는 강물을 따라 내려가는 동안 아무러 관측지점을 못만난다는 얘기인데, 고려해봐야 할 듯
        break
      }else if(mouth_index_imsi==1){
        #하류 지류의 갯수가 2개는 아니라고 가정
        mouth_seg_imsi <- which(as.matrix(adjacency$adjacency)[ind_imsi,]==1)
        if(sum(mouth_seg_imsi==network_indexID)==0){
          #이 경우에는 다음 segment로 넘어가야 한다
          ind_imsi <- mouth_seg_imsi
          #pt_downstream_dist
        }else{
          #이 경우는 해당 segment에 관측장소가 여러 개 있을 경우이다
          #이제 해당 segment안에 가장 가까운 인덱스를 찾는다
          if(sum(mouth_seg_imsi==network_indexID)==1){
            #해당 segment에 관측장소가 단 한개만 있을 때
            ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi)
            nbrs_downstream <- c(nbrs_downstream, pointsin[ind])
          }else{
            #해당 segment에 관측장소가 여러 개 있을 때
            ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi) #해당 segment들에 들어있는 obs들의 좌표들
            #그러면 해당 segment와 가장 가까운 장소가 어디인지 판단해야 한다
            #example_network@lines에서 들어있는 Line들은 상류 -> 하류 순으로 되어있다
            pt_upstream <- example_network@lines[[mouth_seg_imsi]]@Lines[[1]]@coords[1,]
            nearest_pt_index <- ind[which.min(as.matrix(dist(rbind(pt_upstream,example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]) ))[-1,1])]
            nbrs_downstream <- c(nbrs_downstream, pointsin[nearest_pt_index])
          }
          break
        }
      }
    }
  }
  return(list(nbrs=sort(c(nbrs_upstream, nbrs_downstream)), nbrs_upstream=sort(nbrs_upstream), nbrs_downstream=sort(nbrs_downstream)))
}


#고쳐야
getnbrs_stream <- function(remove, pointsin, example_network, I_lines, adjacency, upperExtra=FALSE){
  #upperExtra: Y자형 형태의 네트워크에서 upperExtra부분을 추가할것이냐라는 의미로 이것의 tail-down approach에 해당한다
  r = which(remove==pointsin)
  #나중에 pointsin 추가해야 할지도
  network_SegmentID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[r]
  n_inthe_Segment <- sum(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID==network_SegmentID)
  
  coords_stream <- c()
  for(ii in 1:length(I_lines[[remove]])){
    #ii: r 안에 있는 모든 stream index
    range_up <- I_lines[[remove]][[ii]]@range[1]
    range_down <- I_lines[[remove]][[ii]]@range[2]
    if(range_up!=1){
      range_up <- range_up - 1
    }
    if(range_down!=nrow(I_lines[[remove]][[ii]]@Lines[[1]]@coords)){
      range_down <- range_down + 1
    }
    coords_stream <- rbind(coords_stream, I_lines[[remove]][[ii]]@Lines[[1]]@coords[range_up,], I_lines[[remove]][[ii]]@Lines[[1]]@coords[range_down,])
  }
  coords_stream <- coords_stream[!duplicated(coords_stream),]
  
  n_test <- setdiff(example_network@obspoints@SSNPoints[[1]]@point.data$pid, remove)
  #r 체크하고 나머지들
  
  stream_candidate <- c()
  for(jj in 1:length(n_test)){
    #jj: 노드들 중 r 제거하고 나머지
    lengths_segment <- length(I_lines[[n_test[jj]]])
    match_result <- c()
    for(kk in 1:length(I_lines[[n_test[[jj]]]])){
      #kk: 노드 jj안에 있는 모든 stream index
      match_result_imsi <- apply(I_lines[[n_test[[jj]]]][[kk]]@Lines[[1]]@coords, 1, function(a) sum(which(a%in%coords_stream)))
      match_result <- c(match_result, match_result_imsi)
    }
    if(sum(match_result)!=0){
      if(upperExtra==TRUE){
        #upperExtra==TRUE: Y자의 가지에 해당하는 점도 이웃으로 본다
        stream_candidate <- c(stream_candidate, n_test[jj])
      }else{
        #upperExtra==FALSE: Y자의 가지에 해당하는 점은 이웃으로 안본다
        segment_index <- as.numeric(example_network@obspoints@SSNPoints[[1]]@point.data$rid[[n_test[jj]]])
        #segment_index <- as.numeric(example_network@obspoints@SSNPoints[[1]]@point.data$rid[jj])
        segment_index_r <- as.numeric(example_network@obspoints@SSNPoints[[1]]@point.data$rid[r])
        if(segment_index==segment_index_r){
          #두 점이 같은 segment 안에 있으면 이웃으로 간주
          stream_candidate <- c(stream_candidate, n_test[jj])
        }else{
          distUpstream_index <- example_network@network.line.coords$DistanceUpstream[n_test[jj]]
          #distUpstream_index <- example_network@network.line.coords$DistanceUpstream[segment_index]
          distUpstream_index_r <- example_network@network.line.coords$DistanceUpstream[segment_index_r]
          
          bid_index <- adjacency$rid_bid[segment_index,2]
          bid_index_r <- adjacency$rid_bid[segment_index_r,2]
          if(distUpstream_index_r>distUpstream_index){
            #그러면 제거할 포인트가 이웃보다 상류에 있는 것이다
            if(str_detect(bid_index_r, bid_index) | str_detect(bid_index, bid_index_r)  ){
              adj_imsi <- as.matrix(adjacency$adjacency)
              if(sum(as.matrix(adjacency$adjacency)[,segment_index_r])>=2 | sum(adj_imsi[,as.numeric(as.character(network_SegmentID))])==0){
                #만약 제거지점의 상류가 2개 이상일 경우만 아래지점 이웃 추가
                #또는 제거지점이 최상류에 있을 경우도 이웃 추가
                stream_candidate <- c(stream_candidate, n_test[jj])
              }
            }
          }else{
            #그러면 제거할 포인트가 이웃보다 하류에 있는 것이다
            if(str_detect(bid_index, bid_index_r) | str_detect(bid_index_r, bid_index)){
              stream_candidate <- c(stream_candidate, n_test[jj])
            }
          }
        }
      }
    }
  }
  #최종적으로 같은 노드 안에 이웃이 한개라도 존재하는 경우 그것들만을 subset으로 취한다
  rid_stream_candidate <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[match(stream_candidate,pointsin)]
  rid_stream_r <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[r]
  if(!all(rid_stream_candidate!=rid_stream_r)){
    stream_candidate <- stream_candidate[which(rid_stream_candidate==rid_stream_r)]
  }
  
  return(stream_candidate)
  
  
  #if(length(I_lines[[remove]])==1){
  #  #segment가 1개
  #  coords_upstream <- I_lines[[remove]][[1]]@Lines[[1]]@coords[1,]
  #  coords_downstream <- I_lines[[remove]][[1]]@Lines[[1]]@coords[2,]
  #}else{
  #  #segment가 2개 이상 (여러 개)
  #  for(ii in 1:length(I_lines[[remove]])){
  #    
  #  }
  #}
  
  
  
  #nbd_up <- c(); nbd_down <- c(); nbd_up_extra <- c(); nbd_down_extra <-c()
  #if(n_inthe_Segment!=1){
  #  #(0) 동일 stream 안에서 이웃 탐색
  #  #이 경우에는 같은 segment 내에 관찰점이 적어도 한 개 이상 더 존재하는 것이다
  #  ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID==network_SegmentID)
  #  nbd_ind <- setdiff(ind, remove)
  #  #두 개 이상 존재하더라도 이 이웃들이 전부 remove point 위 또는 아래에만 위치할 수도 있다
  #  node_order_index <- order(example_network@obspoints@SSNPoints[[1]]@point.coords[ind,1])
  #  #node_order_index의 앞쪽에 위치한 것은 상류, 아래쪽에 위치한 것은 하류다
  #  node_order_ind_index <- which(node_order_index==1)
  #  nbd_up <- c(nbd_up, ind[which(node_order_index<node_order_ind_index)])
  #  nbd_down <- c(nbd_down, ind[which(node_order_index>node_order_ind_index)])
  #}
  #if(length(nbd_up)==0 | length(nbd_down)==0){
  #  if(length(nbd_down)==0){
  #    #(우선 같은 segment 내에 하류에 이웃이 있는 경우)
  #    #(1) 상류의 이웃 탐색(없을 수도 있다)
  #    upstream_index <- which(as.matrix(adjacency$adjacency)[,network_SegmentID]!=0) #상류 정보
  #    for(ii in 1:length(upstream_index)){
  #      ind_ii <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID==upstream_index[ii])
  #      #이 ind_ii는 없을 수도 있고, 여러 개 있을 수도 있다
  #      if(length)
  #    }
  #    
  #    as.matrix(adjacency$adjacency)[network_SegmentID,]
  #    as.matrix(adjacency$adjacency)[,network_SegmentID]
  #    
  #    
  #  }else if(length(nbd_up)==0){
  #    #(2) 하류의 이웃 탐색(없을 수도 있다)
  #    #(우선 같은 segment 내에 하류에 이웃이 있는지부터 체크해봐야 한다)
  #  }
  #  
  #  #(3) 바로 밑 segment의 상류 이웃들 탐색 (Y자형 생각해보자)
  #}
  #return(list(nbd_up=nbd_up, nbd_down=nbd_down, nbd_up_extra=nbd_up_extra, nbd_down_extra=nbd_down_extra))
}


#함수: streamPred 함수 만들기
#기본적인 룰: segment 안에서는 inverse distance weight로, 나머지 경우에는 stream 따라서 regression한다
streamPred <- function(pointsin, X=example_network, coeff=sims@obspoints@SSNPoints[[1]]@point.data$Sim_Values, nbrs, remove, adj=adjacency, method=c("IDW","LC"), intercept=FALSE, neighbours=NULL, adaptive=FALSE){
  #prediction weight를 반환하는 함수다
  #method: IDW: inverse distance weight, LC: local constant
  
  #X=example_network; coeff=sims@obspoints@SSNPoints[[1]]@point.data$Sim_Values; adj=adjacency
  
  Xcoords <- X@obspoints@SSNPoints[[1]]@point.coords
  index <- match(nbrs, pointsin)
  r <- which(pointsin == remove)
  Xneighbours <- Xcoords[index,]
  Xremove <- Xcoords[r,]
  #Xneighbours <- Xcoords[nbrs,]
  #Xremove <- Xcoords[remove,]
  
  if (intercept) {
    #Xneighbours <- cbind(1, Xneighbours)
    #Xremove <- as.row(c(1, Xremove))
    #어떤 조치를 취한다
  }
  
  
  if(sum(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)!=0){
    updist_remove <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
    reach_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
    #updist_remove <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==r)]
    #reach_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[which(X@obspoints@SSNPoints[[1]]@point.data$pid==r)]
    addfunccol_remove <- X@obspoints@SSNPoints[[1]]@point.data$shreve[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
  }else{
    updist_remove <- X@predpoints@SSNPoints[[1]]@point.data$upDist[length(X@predpoints@SSNPoints[[1]]@point.data$upDist)]
    reach_remove <- X@predpoints@SSNPoints[[1]]@point.data$rid[length(X@predpoints@SSNPoints[[1]]@point.data$rid)]
    addfunccol_remove <- X@predpoints@SSNPoints[[1]]@point.data$shreve[length(X@predpoints@SSNPoints[[1]]@point.data$addfunccol)]
  }
  rid_remove <- reach_remove
  rid_nbrs <- X@obspoints@SSNPoints[[1]]@point.data$rid[index]
  
  #length_vec <- rep(0, length(nbrs)) #length vec의 inverse weight?
  weight_vec <- rep(0, length(nbrs)) #처음부터 weight를 계산 -> normalize
  
  if(length(nbrs)>=2){
    
    for(i in 1:length(nbrs)){
      nb_index <- nbrs[i]
      #nb_index <- index[i]
      if(method=="IDW" | sum(rid_nbrs==rid_remove)==2){
        #각 nbr point들과 remove point 사이의 거리를 구해야 한다
        
        updist_nb_index <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==nb_index)]
        #updist_remove <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
        
        weight_vec[i] <- 1/abs(updist_nb_index-updist_remove)
        
      }else if(method=="LC"){
        #segment 기반 weight?
        
        addfunccol_nb_index <- X@obspoints@SSNPoints[[1]]@point.data$shreve[which(X@obspoints@SSNPoints[[1]]@point.data$pid==nb_index)]
        #addfunccol_remove <- X@obspoints@SSNPoints[[1]]@point.data$addfunccol[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
        
        weight_vec[i] <- min(addfunccol_nb_index, addfunccol_remove)/max(addfunccol_nb_index, addfunccol_remove)
      }
    }
    #normalizing weight vector
    weights <- weight_vec/sum(weight_vec)
  }else{
    #length가 1인 경우
    #따로 손볼 필요 없음
    mm <- 0
    bhat <- 1
    weights <- 1
    pred <- coeff[nbrs]
  }
  
  #추가: stream segment index도 받는다
  #같은 stream 내에 있을 경우, 그리고 다른 index 내에 있을 경우에 따라 어떻게 prediction을 해야 할 지 결정
  #같은 stream 내에 있을 경우: inverse distance weight
  #다른 stream에 있을 경우: segment를 따라 weight를 준다
  #reach_neighbours <- X@obspoints@SSNPoints[[1]]@point.data$rid[nbrs]
  #reach_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[remove]
  #reach_match <- all(reach_neighbours==reach_remove)
  #
  #if (length(nbrs) >= 2 & reach_match==TRUE & sum(as.numeric(reach_remove))==length(nbrs) ) {
  #  #이 경우 위 아래로 모든 이웃이 다 같은 stream에 있는 것이다
  #  #IDW 적용
  #  
  #}else if(length(nbrs)>=2 & reach_match==FALSE){
  #  temp <- crossprod(Xneighbours)
  #  mm <- Rmatsolve(temp) %*% t(Xneighbours)
  #  bhat <- mm %*% matrix(coeff[nbrs], ncol = 1)
  #  pred <- Xremove %*% bhat
  #  weights <- matrix(Xremove, nrow = 1) %*% mm
  #}else{
  #  #length가 1인 경우
  #  #따로 손볼 필요 없음
  #  mm <- 0
  #  bhat <- 1
  #  weights <- 1
  #  pred <- coeff[nbrs]
  #}
  #return(list(weights = weights, pred = pred, coeff = coeff))
  pred = sum(weights*coeff[nbrs])
  coeff = coeff[remove]
  return(list(weights = weights, pred = pred, coeff = coeff))
}





PointsUpdate_stream <- function (X=example_network, coeff, nbrs, index=NULL, remove, pointsin, weights, lengths=I, lengths_lines=I_lines, method="LC", stream=TRUE, linear=TRUE) 
{
  #X=example_network; lengths=I; lengths_lines=I_lines; method="LC"; index=match(nbrs, pointsin)
  r <- which(pointsin == remove)
  N <- length(pointsin)
  weights_imsi <- (1/weights)/sum(1/weights)
  if(stream==FALSE){
    #stream network가 아닐 경우
    if ((r >= 2) & (r <= (N - 1))) {
      lengths[index] <- as.row(lengths[index])
      weights <- as.row(weights)
      lengths[index] <- lengths[index] + lengths[r] * weights
    }else{
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
    }
    else {
      #q <- which(pointsin == nbrs)
      q <- match(nbrs, pointsin)
      alpha <- lengths[r]/lengths[q]
      coeff[pointsin[q]] <- coeff[pointsin[q]] + alpha * coeff[remove]
    }
    return(list(coeff = coeff, lengths = lengths, r = r, N = N, 
                weights = weights, alpha = alpha))
  }else{
    #제거하려는 장소가 가진 stream이 1개일 때
    #else if(lengths_lines[[remove]][[1]]@Lines==1)
    #stream==TRUE
    rid_nbrs <- X@obspoints@SSNPoints[[1]]@point.data$rid[index]
    rid_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[r]
    
    #I update
    #I_new <- I
    #I_new[nbrs] <- I_new[nbrs] + weights*I[remove]
    lengths[index] <- as.row(lengths[index])
    weights <- as.row(weights)
    lengths[index] <- lengths[index] + lengths[r] * weights
    
    if(linear==TRUE){
      remove_upper <- lengths_lines[[remove]][[1]]@Lines[[1]]@coords[1,]
      remove_lower <- lengths_lines[[remove]][[1]]@Lines[[1]]@coords[2,]
      remove_id <- strsplit(lengths_lines[[remove]][[1]]@ID, "-")[[1]][1]
    }
  
    
    #stream network일 경우
    if(length(nbrs)>=2){
      if(linear==TRUE){
        streamcoord_new <- t(weights_imsi[order(X@network.line.coords$DistanceUpstream[rid_nbrs], decreasing=T)])%*%rbind(remove_upper, remove_lower)
      }
      if(sum(rid_nbrs==rid_remove)==2){
        #이 경우에는 제거점과 그 이웃점들 모두 같은 segment 안에 있는 것이다
        if(method=="IDW"){
          
        }else if(method=="LC"){
          #(1) update lengths_lines
          
          #어떤 이웃이 upper와 닿고 어떤 이웃이 lower와 닿아있는지 분류해야 함
          for(i in 1:length(nbrs)){
            nbr_index <- nbrs[i] #lengths_lines는 remove시키지 않았음에 주의
            for(j in 1:length(lengths_lines[[nbr_index]])){
              if(strsplit(lengths_lines[[nbr_index]][[j]]@ID,"-")[[1]][1] == remove_id){
                if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_upper)){
                  lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,] <- streamcoord_new
                }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_lower)){
                  lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,] <- streamcoord_new
                }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_lower)){
                  #이 경우도 고려 안해도 될듯
                }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_upper)){
                  #upper와 upper가 같을 일은 아마 없을 것이다
                }
              }
            }
          }
          
          #weight만큼의 내분점 찾기
          #lengths_lines[[nbrs[1]]]
          #lengths_lines[[nbrs[2]]]
        }
      }else if(length(unique(c(rid_nbrs, rid_remove)))>=3){
        #Y자형에 제거하려는 점이 밑부분에 있다
        if(length(lengths_lines[[remove]])==1){
          #제거하려는 점이 가진 segement가 1개밖에 없을 경우
          
          for(i in 1:length(nbrs)){
            nbr_index <- nbrs[i]
            org_index <- length(lengths_lines[[nbr_index]])
            #그러면 두 개의 reach는 다른 stream에 있다
            #바로 합쳐주면 됨
            for(j in 1:length(lengths_lines[[remove]])){
              lengths_lines[[nbr_index]][[org_index+1]] <- lengths_lines[[remove]][[j]]
              lengths_lines[[nbr_index]][[org_index+1]]@weight <- weights[i]
            }
          }
          
        }
      }else{
        #이 경우는 이웃한 점들의 reach index와 remove pt의 reach index가 완전히 다른 경우이다
      }
      
    }else{
      if(length(nbrs)==1){
        #special trtment
        #이 경우는 reach index가 서로 같은 경우(같은 물줄기에 있는 경우), 다른 경우로 나누어 구분해야
        nbr_index <- nbrs
        #streamcoord_new 처리할 필요 없음
        if(rid_nbrs==rid_remove){
          #그러면 두 개의 reach는 같은 stream 내에 있다
          for(j in 1:length(lengths_lines[[nbr_index]])){
            if(linear==TRUE){
              if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_upper)){
                lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,] <- remove_lower
              }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_lower)){
                lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,] <- remove_upper
              }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_lower)){
                #이 경우도 고려 안해도 될듯
              }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_upper)){
                #upper와 upper가 같을 일은 아마 없을 것이다
              }
            }else{
              if(lengths_lines[[nbr_index]][[j]]@ID==I_lines[[remove]][[1]]@ID){
                lengths_lines[[nbr_index]][[j]]@range <- range(c(lengths_lines[[nbr_index]][[j]]@range, I_lines[[remove]][[1]]@range))
              }
            }
          }
        }else{
          org_index <- length(lengths_lines[[nbr_index]])
          #그러면 두 개의 reach는 다른 stream에 있다
          #바로 합쳐주면 됨
          for(j in 1:length(lengths_lines[[remove]])){
            lengths_lines[[nbr_index]][[(org_index+j)]] <- lengths_lines[[remove]][[j]]
          }
        }
      }
    }
    
    #update step을 어떻게 할 것인지도 좀 더 고려해봐야
    alpha <- matrix(0, 1, length(nbrs))
    if (length(nbrs) >= 2) {
      alpha <- lengths[r] * lengths[index]/(sum(lengths[index]^2))
      coeff[pointsin[index]] <- coeff[pointsin[index]] + alpha * coeff[remove]
    }
    else {
      q <- which(pointsin == nbrs)
      alpha <- lengths[r]/lengths[q]
      coeff[pointsin[q]] <- coeff[pointsin[q]] + alpha * coeff[remove]
    }
    
    return(list(coeff = coeff, lengths = lengths, lengths_lines = lengths_lines, r = r, N = N, weights = weights, alpha = alpha))
  }
  
}


#lifting forward step
fwtnp_stream <- function(example_network, I, I_lines, adjacency, nkeep=2, linear=FALSE, do.W = TRUE, varonly = FALSE){
  #do.W = TRUE; varonly = FALSE; linear=FALSE; nkeep=2
  
  removelist <- NULL
  lengthsremove <- NULL
  neighbrs <- list()
  gamlist <- list()
  alphalist <- list()
  schemehist <- NULL
  interhist <- NULL
  clolist <- NULL
  W <- v <- NULL
  
  n <- length(I_lines)
  
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
  
  
  pointsin <- matrix(1:n, 1, n)
  
  #nkeep=2
  j<- 1
  
  while(TRUE){
    print(paste("Level ", j))
    if(j==1){
      remove <- order(abs(I))[1] #초기 I만 쓰도록 한다
    }else{
      remove <- order(abs(lengths))[1]
    }
    remove <- pointsin[remove] #원래 인덱스대로 반환됨
    removelist[j] <- remove
    print(paste("remove point is ", remove))
    
    #points(example_network@obspoints@SSNPoints[[1]]@point.coords[which(pointsin==remove),1] , example_network@obspoints@SSNPoints[[1]]@point.coords[which(pointsin==remove),2], cex=2, pch=5, col="red" )
    
    #X_dataframe <- as.data.frame(t(X))
    
    #have to make getnbrs function
    nbrs_obj <- getnbrs_stream3(remove, pointsin, example_network, I_lines, adjacency)
    nbrs <- nbrs_obj$nbrs
    nbrs_upstream <- nbrs_obj$nbrs_upstream
    nbrs_downstream <- nbrs_obj$nbrs_downstream
    print(paste("neighbours of ", remove, " are ", nbrs))
    
    #res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours)
    if(j==1){
      coeff=example_network@obspoints@SSNPoints[[1]]@point.data$Sim_Values
      res <- streamPred(pointsin,  X=example_network, coeff=example_network@obspoints@SSNPoints[[1]]@point.data$Sim_Values, nbrs, remove, adj=adjacency, method="LC", intercept=FALSE, neighbours=NULL)
    }else{
      res <- streamPred(pointsin,  X=example_network, coeff=coeff, nbrs, remove, adj=adjacency, method="LC", intercept=FALSE, neighbours=NULL)
    }
    
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
    }else {
      scheme <- l[[5]]
      int <- l[[4]]
      details <- l[[6]]
    }
    coeff[remove] <- coeff[remove] - pred
    
    #PointsUpdate 함수 만들기: 
    #l1 <- PointsUpdate(X, coeff, nbrs, index, remove, pointsin, weights, lengths)
    #res <- streamPred(pointsin,  X=example_network, coeff=sims@obspoints@SSNPoints[[1]]@point.data$Sim_Values, nbrs, remove, adj=adjacency, method="LC", intercept=FALSE, neighbours=NULL)
    
    l1 <- PointsUpdate_stream(X=example_network, coeff, nbrs, index=match(nbrs, pointsin), remove, pointsin, weights, lengths=I, lengths_lines=I_lines, method="LC", stream=TRUE, linear=linear)
    
    
    coeff <- l1$coeff
    lengths <- l1$lengths
    r <- l1$r
    weights <- l1$weights
    N <- l1$N
    alpha <- l1$alpha
    if (ex){
      if (varonly) {
        W[r, ] <- W[r, ] - colSums(as.vector(weights) *  matrix(W[index, ], nrow = length(nbrs)))
        W[index, ] <- W[index, ] + matrix(alpha) %*% W[r, ]
        v[remove] <- sum(W[r, ]^2)
        np <- setdiff(1:length(pointsin), r)
        W <- W[np, ]
      }else {
        W[remove, ] <- W[remove, ] - colSums(as.vector(weights) * matrix(W[nbrs, ], nrow = length(nbrs)))
        W[nbrs, ] <- W[nbrs, ] + matrix(alpha) %*% W[remove,  ]
      }
    }
    
    lengthsremove[j] <- lengths[r]
    gamlist[[j]] <- weights
    alphalist[[j]] <- alpha
    schemehist[j] <- scheme
    interhist[j] <- int
    lengths <- lengths[setdiff(1:length(pointsin), r)]
    pointsin <- setdiff(pointsin, remove)
    
    ##example network 자체도 업데이트 해야함 : subsetSSN 함수 참고할 것
    #if(j==1){
    #  #example_network@predpoints <- example_network@obspoints
    #  #example_network@predpoints@SSNPoints[[1]]@network.point.coords <- example_network@predpoints@SSNPoints[[1]]@network.point.coords[remove,]
    #  #example_network@predpoints@SSNPoints[[1]]@point.coords <- example_network@predpoints@SSNPoints[[1]]@point.coords[remove,]
    #  #example_network@predpoints@SSNPoints[[1]]@point.data <- example_network@predpoints@SSNPoints[[1]]@point.data[remove,]
    #}else{
    #  
    #}
    example_network@predpoints@SSNPoints[[1]]@network.point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@network.point.coords, example_network@obspoints@SSNPoints[[1]]@network.point.coords[r,])
    example_network@predpoints@SSNPoints[[1]]@point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@point.coords, matrix(example_network@obspoints@SSNPoints[[1]]@point.coords[r,], nrow=1, ncol=2))
    example_network@predpoints@SSNPoints[[1]]@point.data <- rbind(example_network@predpoints@SSNPoints[[1]]@point.data, example_network@obspoints@SSNPoints[[1]]@point.data[r,c(1:ncol(example_network@predpoints@SSNPoints[[1]]@point.data))])
    #if(nrow(example_network@predpoints@SSNPoints[[1]]@network.point.coords)!=1){
      #example_network@predpoints@SSNPoints[[1]]@network.point.coords <- example_network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(example_network@predpoints@SSNPoints[[1]]@network.point.coords),]
      #example_network@predpoints@SSNPoints[[1]]@point.coords <- matrix(example_network@predpoints@SSNPoints[[1]]@point.coords[nrow(example_network@predpoints@SSNPoints[[1]]@point.coords),], nrow=ncol(example_network@predpoints@SSNPoints[[1]]@point.coords), ncol=1)
      #example_network@predpoints@SSNPoints[[1]]@point.data <- example_network@predpoints@SSNPoints[[1]]@point.data[nrow(example_network@predpoints@SSNPoints[[1]]@point.data),]
    #}
    
    example_network@obspoints@SSNPoints[[1]]@network.point.coords <- example_network@obspoints@SSNPoints[[1]]@network.point.coords[-r,]
    example_network@obspoints@SSNPoints[[1]]@point.coords <- example_network@obspoints@SSNPoints[[1]]@point.coords[-r,]
    example_network@obspoints@SSNPoints[[1]]@point.data <- example_network@obspoints@SSNPoints[[1]]@point.data[-r,]
    
    #I update
    I <- l1$lengths[-r]
    #I_lines update (list의 경우에는 subset 하는 것이 불가능)
    I_lines <- l1$lengths_lines
    
    
    ##plotting new data object
    SSN::plot.SpatialStreamNetwork(x=example_network, main=paste("Level", j), VariableName = "Sim_Values", color.palette=rainbow(12), nclasses=10, breaktype = "quantile", brks = NULL, PredPointsID = NULL, add = FALSE, addWithLegend = FALSE, lwdLineCol = "addfunccol", lwdLineEx = 3, lineCol = "black")
    #need to be corrected
    
    #if(length(pointsin)<=nkeep | sum(as.matrix(adjacency$adjacency)[pointsin, pointsin])==0){
    if(length(pointsin)<=nkeep){
      break
    }else{
      #next level
      
      j = j+1
    }
  }
  if (varonly) {
    v[pointsin] <- rowSums(W^2)
    W <- NULL
  }
  return(list(network=example_network, I=I, I_lines=I_lines, adjacency=adjacency, coeff=coeff, lengths=lengths, lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist, 
              neighbrs = neighbrs, schemehist = schemehist, interhist = interhist, 
              clolist = clolist, gamlist = gamlist, alphalist = alphalist, 
              W = W, v=v))
}

#lifting backward step


UndoPointsUpdate_stream <- function (X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, I_lines, stream=TRUE) 
{
  #I_lines=X$I_lines; X=X$network;  stream=TRUE
  
  
  #stream network를 다시 되돌리는 함수를 만들어야 한다
  gamweights_imsi <- (1/gamweights)/sum(1/gamweights)
  if(stream==TRUE){
    #X$network@obspoints를 복원해야 함
    #locID(=자기자신의숫자), upDist, pid(=자기자신의숫자), rid, ratio, shreve, addfunccol, X, X2, Sim_Values
    #point.data_reconstruct <- c()
    
    rid_remove <- X@predpoints@SSNPoints[[1]]@point.data$rid[length(X@predpoints@SSNPoints[[1]]@point.data$rid)]
    rid_nbrs <- X@obspoints@SSNPoints[[1]]@point.data$rid[index]
    
    remove_upper <- I_lines[[remove]][[1]]@Lines[[1]]@coords[1,]
    remove_lower <- I_lines[[remove]][[1]]@Lines[[1]]@coords[2,]
    remove_id <- strsplit(I_lines[[remove]][[1]]@ID, "-")[[1]][1]
    
    remove_coords <- X@predpoints@SSNPoints[[1]]@point.coords[nrow(X@predpoints@SSNPoints[[1]]@point.coords),]
    
    if(length(nbrs)>=2){
      if(sum(rid_nbrs==rid_remove)==2){
        #이 경우에는 제거점과 그 이웃점들 모두 같은 segment 안에 있는 것이다
        coords_nbrs <- X@obspoints@SSNPoints[[1]]@point.coords[index,]
        #midpt_nbrs <- c(mean(coords_nbrs[,1]), mean(coords_nbrs[,2]))
        
        nbrs_order_index <- order(X@obspoints@SSNPoints[[1]]@point.data$upDist[index], decreasing=T)
        ordered_nbrs_coords <- X@obspoints@SSNPoints[[1]]@point.coords[index[nbrs_order_index],]
        
        #middle_pts <- c() #중점들을 저장할 벡터
        middle_pts_imsi <- c((ordered_nbrs_coords[1,1]+remove_coords[1])/2,  (ordered_nbrs_coords[1,2]+remove_coords[2])/2)
        middle_pts_imsi2 <- c((ordered_nbrs_coords[2,1]+remove_coords[1])/2,  (ordered_nbrs_coords[2,2]+remove_coords[2])/2)
        middle_pts <- rbind(middle_pts_imsi, middle_pts_imsi2)
        
        for(ii in 1:length(nbrs)){
          ID_remove <- c()
          for(kk in 1:length(I_lines[[remove]])){
            ID_remove <- c(ID_remove, strsplit(I_lines[[remove]][[1]]@ID, "-")[[kk]][1])
          }
          ll <- length(I_lines[[nbrs[ii]]])
          while(TRUE){
            if(strsplit(I_lines[[nbrs[ii]]][[ll]]@ID, "-")[[1]][1] %in% ID_remove){
              if(length(strsplit(I_lines[[nbrs[ii]]][[ll]]@ID, "-")[[1]])==1){
                I_lines[[nbrs[ii]]] <- I_lines[[nbrs[ii]]][-ll]
              }else{
                if(X@obspoints@SSNPoints[[1]]@point.data$upDist[index[ii]]>X@predpoints@SSNPoints[[1]]@point.data$upDist[length(X@predpoints@SSNPoints[[1]]@point.data$upDist)]){
                  I_lines[[nbrs[ii]]][[ll]]@Lines[[1]]@coords[2,] <- middle_pts[1,]
                }else{
                  I_lines[[nbrs[ii]]][[ll]]@Lines[[1]]@coords[1,] <- middle_pts[2,]
                }
              }
              ll <- ll - 1
            }else{
              ll <- ll - 1
            }
            if(ll==0){
              break
            }
          }
        }
        
      }else if(length(unique(c(rid_nbrs, rid_remove)))>=3){
        #Y자형에 제거하려는 점이 밑부분에 있다
        for(ii in 1:length(nbrs)){
          ID_remove <- c()
          for(kk in 1:length(I_lines[[remove]])){
            ID_remove <- c(ID_remove, strsplit(I_lines[[remove]][[1]]@ID, "-")[[kk]][1])
          }
          ll <- length(I_lines[[nbrs[ii]]])
          while(TRUE){
            if(strsplit(I_lines[[nbrs[ii]]][[ll]]@ID, "-")[[1]][1] %in% ID_remove){
              I_lines[[nbrs[ii]]] <- I_lines[[nbrs[ii]]][-ll]
              ll <- ll - 1
            }else{
              break
            }
          }
        }
      }
      
    }else if(length(nbrs)==1){
      if(rid_remove==rid_nbrs){
        #두 개의 stream edge가 같은 경우, 조금 복잡함
      }else{
        #두 개의 stream edge가 다른 경우, 그냥 nbrs쪽의 I_lines에서 remove가 속한 노드의 것을 삭제해 주면 될듯
        ID_remove <- c()
        for(kk in 1:length(I_lines[[remove]])){
          ID_remove <- c(ID_remove, strsplit(I_lines[[remove]][[1]]@ID, "-")[[kk]][1])
        }
        ll <- length(I_lines[[nbrs]])
        while(TRUE){
          if(strsplit(I_lines[[nbrs]][[ll]]@ID, "-")[[1]][1] %in% ID_remove){
            I_lines[[nbrs]] <- I_lines[[nbrs]][-ll]
            ll <- ll - 1
          }else{
            break
          }
        }
      }
    }
  }
  
  alpha <- matrix(0, 1, length(nbrs))
  #if ((r > 1) & (r <= N)) {
    alpha <- lengths[index] * lengthrem/(sum(lengths[index]^2))
    coeff[nbrs] <- coeff[nbrs] - alpha * coeff[remove]
    lengths[index] <- as.row(lengths[index])
    prod <- gamweights * lengthrem
    prod <- as.row(prod)
    lengths[index] <- lengths[index] - prod
  #}
  #if ((r == 1) | (r == (N + 1))) {
  #  q <- which(pointsin == nbrs)
  #  alpha <- lengthrem/lengths[q]
  #  coeff[pointsin[q]] <- coeff[pointsin[q]] - alpha * coeff[remove]
  #  lengths[q] <- lengths[q] - lengthrem
  #}
  
  return(list(coeff = coeff, lengths = lengths, alpha = alpha, I_lines=I_lines))
}

invtnp_stream <- function (X, coeff, lengths, lengthsremove, pointsin, removelist, 
            neighbrs, schemehist, interhist, nadd = length(X) - 2, intercept = TRUE, 
            neighbours = 1,  closest = FALSE, LocalPred = LinearPred, length_diff) {
  
  # inverse lifting step
  if (is.list(X)) {
    coeff <- X$coeff
    lengths <- X$lengths
    lengthsremove <- X$lengthsremove
    removelist <- X$removelist
    neighbrs <- X$neighbrs
    pointsin <- X$pointsin
    schemehist <- X$schemehist
    interhist <- X$interhist
    #X <- X$network
  }
  
  #X <- as.row(X)
  coeff <- as.row(coeff)
  #n <- length(X)
  n <- length(coeff)
  N <- length(pointsin)
  m <- length(removelist) #n = N+m
  #d <- neighbours
  if (nadd > 0) {
    for (j in 1:nadd) {
      N <- length(pointsin)
      remove <- removelist[m - j + 1]
      lengthrem <- lengthsremove[m - j + 1]
      nbrs <- neighbrs[[m - j + 1]]
      index <- NULL
      index <- match(nbrs, pointsin)
      
      #decide r index
      #B <- (X[remove] > X[nbrs])
      B <- (remove > nbrs)
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
          res <- LinearPred(pointsin, X, coeff, nbrs, 
                            remove, intercept = interhist[m - j + 1], 
                            neighbours)
        }
        if (schemehist[m - j + 1] == "Quad") {
          res <- QuadPred(pointsin, X, coeff, nbrs, remove, 
                          intercept = interhist[m - j + 1], neighbours)
        }
        if (schemehist[m - j + 1] == "Cubic") {
          res <- CubicPred(pointsin, X, coeff, nbrs, 
                           remove, intercept = interhist[m - j + 1], 
                           neighbours)
        }
      }else {
        #schemehist가 null인 경우
        
        res <- streamPred(pointsin,  X=X$network, coeff=coeff, nbrs, remove, adj=X$adjacency, method="LC", intercept=FALSE, neighbours=NULL)
        
        #res <- LocalPred(pointsin, X, coeff, nbrs, remove,intercept, neighbours)
      }
      if (length(res) == 2) {
        l <- res[[1]]
      }else {
        l <- res
      }
      gamweights <- l$weights
      #UndoPointsUpdate함수의 stream version을 만들어야
      #l1 <- UndoPointsUpdate(X, coeff, nbrs, index, remove,  r, N, pointsin, gamweights, lengths, lengthrem)
      l1 <- UndoPointsUpdate_stream(X=X$network, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, I_lines=X$I_lines, stream=TRUE)
      X$I_lines <- l1$I_lines
      coeff <- l1$coeff
      lengths <- l1$lengths
      pred <- sum(as.column(gamweights) * coeff[nbrs])
      coeff[remove] <- coeff[remove] + pred
      removelist <- setdiff(removelist, remove)
      
      if (r == 1) {
        lengths <- c(lengthrem, lengths)
        pointsin <- c(remove, pointsin)
        #update example network
        X$network@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(X$network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),], X$network@obspoints@SSNPoints[[1]]@network.point.coords)
        X$network@obspoints@SSNPoints[[1]]@point.coords <- rbind(X$network@predpoints@SSNPoints[[1]]@point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),], X$network@obspoints@SSNPoints[[1]]@point.coords)
        if(length_diff<=0){
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X$network@obspoints@SSNPoints[[1]]@point.data)
        }else{
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(cbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X=NA, X2=NA, Sim_Values=NA) , X$network@obspoints@SSNPoints[[1]]@point.data)
        }
      }
      if (r == (N + 1)) {
        lengths <- c(lengths, lengthrem)
        pointsin <- c(pointsin, remove)
        #update example network
        X$network@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@network.point.coords, X$network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),])
        X$network@obspoints@SSNPoints[[1]]@point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@point.coords, X$network@predpoints@SSNPoints[[1]]@point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),])
        if(length_diff<=0){
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data, X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),])
        }else{
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data, cbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X=NA, X2=NA, Sim_Values=NA))
        }
      }
      if ((r > 1) & (r < (N + 1))) {
        lengths <- c(lengths[1:(r - 1)], lengthrem, lengths[r:N])
        pointsin <- c(pointsin[1:(r - 1)], remove, pointsin[r:N])
        #update example network
        X$network@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@network.point.coords[c(1:(r - 1)),], X$network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),], X$network@obspoints@SSNPoints[[1]]@network.point.coords[c(r:N),])
        X$network@obspoints@SSNPoints[[1]]@point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@point.coords[c(1:(r - 1)),], X$network@predpoints@SSNPoints[[1]]@point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),], X$network@obspoints@SSNPoints[[1]]@point.coords[c(r:N),])
        if(length_diff<=0){
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data[c(1:(r - 1)),], X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X$network@obspoints@SSNPoints[[1]]@point.data[c(r:N),])
        }else{
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data[c(1:(r - 1)),], cbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X=NA, X2=NA, Sim_Values=NA), X$network@obspoints@SSNPoints[[1]]@point.data[c(r:N),])
        }
      }
      #delete points from predpoint
      X$network@predpoints@SSNPoints[[1]]@network.point.coords <- X$network@predpoints@SSNPoints[[1]]@network.point.coords[-nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),]
      X$network@predpoints@SSNPoints[[1]]@point.coords <- X$network@predpoints@SSNPoints[[1]]@point.coords[-nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),]
      X$network@predpoints@SSNPoints[[1]]@point.data <- X$network@predpoints@SSNPoints[[1]]@point.data[-nrow(X$network@predpoints@SSNPoints[[1]]@point.data),]
    }
  }
  
  
  
  return(list(coeff = coeff, lengths = lengths, lengthsremove = lengthsremove, 
              pointsin = pointsin, removelist = removelist))
}
    
    
########################################
##Data-adaptive stream flow weight selection
########################################  
 









########################################
##Network poltting function도 나중에 만들어보자
########################################








########################################
##SSN에서 SpatialStreamNetwork object 만드는 함수
########################################
importSSN_stream <- function(shreve_obj, location="Miho", multipleplaces=FALSE, RID="original"){
  #location: "Miho", "Full" 있음
  #multipleplaces=FALSE 이면 임의로 고치고 합을 만든다
  #multipleplaces=TRUE이면 그냥 한다
  
  #netID = 1 (왜냐면 네트워크가 1개이므로)
  
  #(a) stream network 자료
  edges  <- readOGR("~/Dropbox/Maps/RiverFlow/KRF_3.0_Geumgang/KRF_ver3_LINE_금강수계.shp")
  #rownames(edges@data) <- edges@data[, "RCH_ID"]
  #reach ID (rid 출력)
  
  #(b) 관측 장소 자료
  sites <- read.csv("~/Dropbox/Data/Riverflow/금강장소(엑셀편집)UTF8.csv",header=T)
  sites <- sites[which(sites$총량==1),]
  sites <- data.frame(name= sites$이름,x=sites$경도.Degree., y=sites$위도.Degree.)
  
  #rownames(sites@data) <- sites@data[,"pid"]
  #rownames(sites@coords) <- sites@data[, "pid"]
  
  #(c) 수질 데이터 자료
  alldata1 <- read.table("~/Dropbox/Data/Riverflow/금강/수질총량/수질20082013.txt",sep=',', header=T)
  alldata2 <- read.table("~/Dropbox/Data/Riverflow/금강/수질총량/수질20132018.txt",sep=',', header=T)
  alldata <- rbind(alldata1, alldata2)
  
  #(d) 예전에 편집한 자료
  if(location=="Miho"){
    data <- readRDS("~/Dropbox/Data/RiverFlow/금강/Geum(miho).RDS")
  }else if(location=="Full"){
    data <- readRDS("~/Dropbox/Data/RiverFlow/금강/Geum.RDS")
  }else{
    stop("Wrong location")
  }
  shape <- data$shape
  #2012년부터 자를까
  #data$eventdata <- data$eventdata["2012/2017",] #이거 없는데 체크해봐야 함 (여기 말고 밑에서는 쓰지 않는다)
  data$streamdata <- data$streamdata["2012/2017",]
  #권역별 출력
  region_index <- unique(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)))[c(1)]
  #ID의 처음 4자리가 권역별 코드
  #3011: 미호천
  #3012: 금강공주
  #3008: 대청댐
  #3007: 보청천
  #3005: 초강
  #3009: 갑천
  #3013: 논산천
  #3006: 대청댐상류
  #3014: 금강하구언
  #3004: 영동천
  #3003: 무주남대천
  #3001: 용담댐
  #3002: 용담댐하류
  #3310: 대청댐하류
  #3301: 만경강(x)
  #3303: 직소천(x)
  #3302: 동진강(x)
  #3101: 삽교천(x)
  #2005: 금강상류부분(필요없음)
  #3202: 부남방조제(x)
  #1101, 3201: 대호방조제(x)
  #3203: 금강서해(x)
  #5301, 5302: 주진천(x)
  
  #for(i in 1:length(region_index)){
  #  shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  #  plot(shape, main=region_index[i])
  #  plot(shape_sub, col="red", add=T, lwd=3)
  #}
  
  name_index <- c("미호천", "금강공주", "대청댐", "보청천", "초강", "갑천", "논산천", "대청댐상류", "금강하구언", "영동천", "무주남대천", "용담댐", "용담댐하류", "대청댐하류")
  name_subindex <- list()
  name_subindex[[1]] <- c("미호댐", "원남댐", "병천천상류", "병천천하류", "보강천", "미호천중류", "미호천상류", "작천보", "한천", "백곡천", "조천", "석화수위표", "미호천하류", "무심천", "백곡댐" ) #15개
  name_subindex[[2]] <- c("세종보", "지천상류", "지천합류후", "어천합류후", "용수천", "백제보", "석성천", "규암수위표", "금천", "논산천합류전", "유구천", "공주보", "대교천", "공주수위표") #14개
  name_subindex[[3]] <- c("대청댐", "대청댐조정지", "소옥천상류", "소옥천하류", "대청댐상류") #5개
  name_subindex[[4]] <- c("항건천", "보청천중류", "삼가천합류후", "삼가천", "보청천상류", "보청천하류") #6개
  name_subindex[[5]] <- c("석천", "초강상류", "초강하류") #3개
  name_subindex[[6]] <- c("갑천하류", "유성수위표", "갑천상류", "대전천", "유등천상류", "유등천하류") #6개
  name_subindex[[7]] <- c("노성천", "논산천하류", "강경천", "논산천상류", "탑정댐") #5개
  name_subindex[[8]] <- c("보청천합류전") #1개
  name_subindex[[9]] <- c("입포수위표", "길산천", "금강하구언") #3개
  name_subindex[[10]] <- c("영동천", "봉황천하류", "초강합류후", "봉황천상류", "봉황천합류후", "호탄수위표") #6개
  name_subindex[[11]] <- c("무주남대천상류", "무주남대천중류", "무주남대천하류") #3개
  name_subindex[[12]] <- c("주자천", "용담댐", "정자천", "구량천", "진안천", "장계천합류후", "장계천", "진안천합류후") #8개
  name_subindex[[13]] <- c("용담댐하류") #1개
  name_subindex[[14]] <- c("미호천합류전", "매포수위표") #2개
  region_subindex <- list()
  
  Geum_RCH_ID <- c()
  for(i in 1:length(region_index)){
    shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
    #plot(shape_sub)
    region_subindex[[i]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6)))
    for(j in 1:length(region_subindex[[i]])){
      shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
      #plot(shape_sub_sub, col="red", add=T)
      Geum_RCH_ID_Imsi <- cbind(shape_sub_sub@data, data.frame(중권역번호=region_index[i], 중권역이름=name_index[i], 소권역번호=region_subindex[[i]][j], 소권역이름=name_subindex[[i]][j]))
      Geum_RCH_ID <- rbind(Geum_RCH_ID, Geum_RCH_ID_Imsi)
    }
  }
  
  #subregion 출력시
  #i=8
  #shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  #plot(shape_sub)
  #region_subindex[[i]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6)))
  #region_subsubindex[[i]] <- list()
  #
  #j=1
  #shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
  #plot(shape_sub_sub, col="red", add=T)
  #region_subsubindex[[i]][[j]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=8)))
  #
  #for(k in 1:length(region_subsubindex[[i]][[j]])){
  #  shape_sub_sub_sub <- subset(shape_sub_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=8))==region_subsubindex[[i]][[j]][k])
  #  plot(shape_sub_sub_sub, col=k, add=T)
  #}
  
  ########################################
  ##Adjacency matrix
  ########################################
  #reference: https://cran.r-project.org/web/packages/riverdist/vignettes/riverdist_vignette.html
  
  #connections: Object of class "matrix", with "numeric" elements. Defined as a square matrix, with elements describing the type of connection detected between line segments.
  #• A value of 1 in element [i,j] indicates that the beginning of segment i is connected to the beginning of segment j.
  #• A value of 2 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j.
  #• A value of 3 in element [i,j] indicates that the end of segment i is connected to the beginning of segment j.
  #• A value of 4 in element [i,j] indicates that the end of segment i is connected to the end of segment j.
  #• A value of 5 in element [i,j] indicates that segments i and j are connected at both beginning and end.
  #• A value of 6 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j, and the end of segment i is connected to the beginning of segment j.
  #• A value of NA in element [i,j] indicates that segments i and j are not connected.
  if(location=="Miho"){
    flow_network_original <- line2network(path = "~/Dropbox/R files/StreamFLow/River/", layer="GeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    flow_network <- line2network(path = "~/Dropbox/R files/StreamFLow/River/", layer="thinGeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    #plot(flow_network)
  }else if(location=="Full"){
    flow_network_original <- line2network(path = "~/Dropbox/R files/StreamFLow/River/", layer="Geum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    flow_network <- line2network(path = "~/Dropbox/R files/StreamFLow/River/", layer="thinGeum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    #plot(flow_network)
  }else{
    stop("Wrong location")
  }
  
  adj_matrix <- flow_network$connections #네트워크의 adjacecy를 표현하였다
  
  #topologydots(rivers=flow_network) #topologydots: 연결되어있으면 초록, 아니면 빨강
  
  #Basic distance calculation in a non-messy river network
  #detectroute(start=4, end=2, rivers=flow_network)
  #riverdistance(startseg=4, startvert=1, endseg=2, endvert=1, rivers=flow_network, map=TRUE)
  
  #subnetwork 출력 예시
  #i=1
  #flow_subnetwork <- trimriver(trimto = c(which(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=6)==region_subindex[[1]][i]) )), rivers=flow_network)
  #plot(flow_subnetwork)
  
  #하구 설정
  if(location=="Miho"){
    flow_network <- setmouth(seg=68,vert=2, rivers=flow_network)  #68번노드가 하구에 해당
    flow_network_original <- setmouth(seg=68,vert=20, rivers=flow_network_original)  #68번노드가 하구에 해당
  }else if(location=="Full"){
    flow_network <- setmouth(seg=614,vert=2, rivers=flow_network)  #68번노드가 하구에 해당
    flow_network_original <- setmouth(seg=614,vert=8, rivers=flow_network_original)  #68번노드가 하구에 해당
  }else{
    stop("Invalid location")
  }
  
  #showends(seg=80, rivers=flow_network) #강물줄기 체크 가능
  
  #내가 추가한 물줄기 길이들
  #reach_length <- shape@data$Shape_Leng[which(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==rregion_subindex[[1]][i]) )]
  
  #edges, sites, alldata를 서로 연동시키는 작업이 필요
  
  #ind1 <- colnames(edges@data) == c("netID")
  #ind2 <- colnames(edges@data) == c("rid")
  #ind3 <- colnames(edges@data) == c("upDist")
  #if (is.factor(edges@data$netID)) {
  #  edges@data$netID <- as.character(edges@data$netID)
  #}
  
  if(location=="Miho"){
    #RID 잘 맞게 변경
    if(multipleplaces==FALSE){
      RCH_ID_manual <- c(38, 30, 40, 44, 45, 63, 50, 46, 1, 88, 22, 87, 27, 84, 62,
                         105, 109, 83, 54, 81, 18, 16, 19, 80, 78, 90, 77, 56, 68)
    }else{
      RCH_ID_manual <- c(38, 30, 40, 44, 45, 63, 50, 46, 1, 88, 22, 87, 27, 84, 62,
                        105, 109, 83, 54, 81, 18, 16, 19, 81, 78, 90, 90, 56, 68)
    }
    
    #data$eventplace$RCH_ID <- shape@data$RCH_ID[c(38, 30, 40, 44, 45, 63, 50, 46, 1, 88, 22, 87, 27, 84, 62,
    #                                              105, 109, 83, 54, 81, 18, 16, 19, 81, 78, 90, 90, 56, 68)]
  
    data$eventplace$RCH_ID <- shape@data$RCH_ID[RCH_ID_manual]
    #data$eventplace$RCH_ID : 원래 있음
    
    #location마다 가장 가까운 lines의 upstream에 매칭을 시켜준다
  }else if(location=="Full" & multipleplaces==FALSE){
    
      #which(data$eventplace$X=="미호천4" | data$eventplace$X=="유등천A" | data$eventplace$X=="삼가천" | data$eventplace$X=="석천")
      #포함이므로 겹치는지점에 메모해 놓은 것과 장소 번호가 다르다.
      #17: 추풍령천, 18: 초강1
      data$eventplace$RCH_ID[17] <- as.integer(30050108)
      #23: 옥천, 24: 우산 (옥천으로 평균)
      #36: 대청댐, 38: 대청
      data$eventplace$RCH_ID[36] <- as.integer(30080409); data$eventplace$RCH_ID[38] <- as.integer(30080419)
      #43: 갑천1, 44: 갑천2
      data$eventplace$RCH_ID[44] <- as.integer(30090208)
      #47: 침산교, 48: 유등천1
      data$eventplace$RCH_ID[47] <- as.integer(30090502)
      #49: 옥계교, 50: 대전천1 (두개 평균을내어야 (대전천1로 평균))
      #55: 갑천5, 56: 갑천5-1 (갑천5로 평균)
      #58: 청원-1 59: 백천
      data$eventplace$RCH_ID[58] <- as.integer(30100203)
      #79: 미호천8, 83: 미호천5
      data$eventplace$RCH_ID[79] <- as.integer(30111301)
      #85: 조천, 86: 조천1
      data$eventplace$RCH_ID[86] <- as.integer(30111505)
      #98: 곰나루, 99: 금강
      data$eventplace$RCH_ID[98] <- as.integer(30120532)
      
      
      #[1] 구량천       가막         진안천       정자천       용담         용포         무주남대천   무주남대천-1 무주남대천-2 부리        
      #[11] 제원         봉황천       제원A        영동         영동천1      영동천2      추풍령천     초강1        금계천       초강2       
      #[21] 이원         옥천         우산         보청천1      보청천2      항건천       보청천3      보청천4      상곡천       추풍천      
      #[31] 옥천천       회인천       주원천       대청댐       품곡천       대청         용호천       현도         두계천1      봉곡2교     
      #[41] 갑천1        갑천2        갑천3        침산교       유등천1      옥계교       대전천1      대전천2      대전천3      유등천5     
      #[51] 갑천4        갑천5        갑천5-1      세종1        청원-1       백천         미호천1      칠장천       미호천1-1    한천        
      #[61] 미호천2      백곡천1      백곡천2      미호천7      초평천       미호천3      보강천1      보강천       성암천       무심천      
      #[71] 무심천1      무심천2      무심천3      석남천       미호천8      병천천A      용두천       병천천       미호천5      미호천5A    
      #[81] 조천         조천1        월하천       미호천6-1    연기         금남         용수천-1     용수천       대교천       대교천2     
      #[91] 세종2        공주1        정안천       곰나루       금강         유구천       목면         공주2        부여         지천        
      #[101] 지천-1       정동         은산천       부여1        금천         부여2        석성천-1     석성천       석성천2      성동        
      #[111] 논산천-1     논산천1      연산천       노성천-1     노성천       논산천2      방축천       논산천4      수철천       마산천      
      #[121] 강경천       강경         산북천       양화-1       길산천       길산천2      금강갑문    
    #}else{
    #  
    #}
    
  }
  dist.pts <- rep(0, nrow(data$eventplace))
  dist.imsi <- c()
  for(jj in 1:nrow(data$eventplace)){
    pts.coords <- c(data$eventplace$경도.Degree.[jj], data$eventplace$위도.Degree.[jj])
    
    kk <- match(data$eventplace$RCH_ID[jj], shape@data$RCH_ID) #이게 부정확한 것 같아서 수정하자
    #실제 RID랑 조금 잘 안맞는 것 같다
    #dist_upper <- apply(shape@data[,c(3,4)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
    #dist_lower <- apply(shape@data[,c(5,6)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
    #
    #kk <- sort(unique(union(which(dist_upper==min(c(dist_upper, dist_lower))), which(dist_lower==min(c(dist_upper, dist_lower))))))
    #
    #dist.imsiimsi <- c()
    #for(pp in 1:length(kk)){
    #  for(ll in 1:nrow(flow_network_original$lines[[kk[pp]]])){
    #    dist.compute <- sqrt(sum((flow_network_original$lines[[kk[pp]]][ll,]-pts.coords)^2))
    #    dist.imsiimsivec <- matrix(c(dist.compute, kk[pp], ll), nrow=1)
    #    dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
    #  }
    #}
    #dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
    #dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
    
    dist.imsiimsi <- c()
    for(ll in 1:nrow(flow_network_original$lines[[kk]])){
      dist.compute <- sqrt(sum((flow_network_original$lines[[kk]][ll,]-pts.coords)^2))
      dist.imsiimsivec <- matrix(c(dist.compute, kk, ll), nrow=1)
      dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
    }
    dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
    dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
  }
  colnames(dist.imsi) <- c("mindist", "rid", "segid")
  rownames(dist.imsi) <- data$eventplace$X
  
  if(location=="Miho"){
    #upDist 또한 구해야 한다
    DistanceUpstream_seg <- rep(0, nrow(shape@data))
    for(ii in 1:length(DistanceUpstream_seg)){
      DistanceUpstream_seg[ii] <- upstream(startseg = 68, endseg = ii, startvert = 20, endvert = 1, rivers=flow_network_original, net=TRUE)
    }
    #관측장소의 DistanceUpstream: dist.imsi활용
    DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
    for(ii in 1:length(DistanceUpstream_obs)){
      DistanceUpstream_obs[ii] <- upstream(startseg = 68, endseg = dist.imsi[ii,2], startvert = 20, endvert = dist.imsi[ii,3], rivers=flow_network_original, net=TRUE)
    }
    shape_tinned <- readOGR("~/Dropbox/R files/StreamFlow/River/thinGeumMiho.shp")
  }else if(location=="Full"){
    #upDist 또한 구해야 한다
    DistanceUpstream_seg <- rep(0, nrow(shape@data))
    for(ii in 1:length(DistanceUpstream_seg)){
      DistanceUpstream_seg[ii] <- upstream(startseg = 614, endseg = ii, startvert = 8, endvert = 1, rivers=flow_network_original, net=TRUE)
    }
    #관측장소의 DistanceUpstream: dist.imsi활용
    DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
    for(ii in 1:length(DistanceUpstream_obs)){
      DistanceUpstream_obs[ii] <- upstream(startseg = 614, endseg = dist.imsi[ii,2], startvert = 8, endvert = dist.imsi[ii,3], rivers=flow_network_original, net=TRUE)
    }
    shape_tinned <- readOGR("~/Dropbox/R files/StreamFlow/River/thinGeum.shp")
  }
  
  #(1)network.line.coords: network.line.coords에 해당
  network.line.coords <- data.frame(NetworkID=as.factor(rep(1,nrow(shape@data))), SegmentID=c(0:(nrow(shape@data)-1)), DistanceUpstream=DistanceUpstream_seg)
  #observeddata <- SSNPoints로 되어있음
  if(location=="Miho"){
    miho <- which(data$eventplace$X=="미호천4")
  }else if(location=="Full"){
    #삼가천, 석천: 시계열에 관찰값이 하나도 없는 장소들
    #유등천A, 미호천4: data$eventplace$X에는 있는데 colnames(data$normaldata_TN)는 없는 장소들
    miho <- which(data$eventplace$X=="미호천4" | data$eventplace$X=="유등천A" | data$eventplace$X=="삼가천" | data$eventplace$X=="석천")
    if(multipleplaces==FALSE){
      miho <- sort(c(miho, which(data$eventplace$X=="무주남대천" | data$eventplace$X=="무주남대천-1" | data$eventplace$X=="우산" | data$eventplace$X=="옥계교" | data$eventplace$X=="갑천5-1")))
    }else{
      miho <- sort(miho)
    }
  }
  #(2)obspoint: obspoints에 해당
  network.line.coords.obsdata <- data.frame(NetworkID=rep(1,(nrow(data$eventplace)-length(miho))), SegmentID=dist.imsi[-miho,2], DistanceUpstream=DistanceUpstream_obs[-miho])
  pointcoords <- cbind(data$eventplace$경도.Degree.[-miho], data$eventplace$위도.Degree.[-miho]); colnames(pointcoords) <- c("coords.x1", "coords.x2")
  #pointdata에 ratio 추가해야
  if(multipleplaces==TRUE & location=="Full"){
    if(RID=="original"){
      RCH_ID_manual <- c(403, 445, 792, 418, 807, 812, 869, 869, 869, 821,
                         826, 827, 830, 833, 788, 789, 228, 228, 861, 867,
                         662, 669, 669, 539, 513, 14, 516, 877, 222, 199,
                         507, 577, 921, 153, 530, 750, 90, 718, 174, 175,
                         579, 579, 584, 569, 569, 182, 182, 856, 855, 591,
                         589, 590, 590, 739, 508, 508, 235, 142, 237, 241,
                         242, 509, 247, 243, 1, 612, 93, 611, 128, 503, 
                         769, 773, 607, 325, 605, 76, 74, 77, 605, 602,
                         754, 754, 448, 592, 738, 737, 62, 49, 481, 476,
                         700, 720, 366, 752, 752, 901, 687, 677, 674, 890,
                         888, 928, 163, 663, 853, 660, 162, 841, 836, 647,
                         775, 554, 190, 560, 555, 551, 252, 646, 785, 919,
                         576, 645, 270, 630, 884, 883, 614) #0번부터 시작하므로 1을 붙여줘야 맞는 듯(2019 08 25)
    }else if(RID=="corrected"){
      RCH_ID_manual <- c(403, 445, 792, 418, 807, 813, 869, 815, 815, 821,
                         826, 294, 830, 833, 789, 789, 858, 228, 861, 867,
                         662, 669, 669, 539, 520, 14, 516, 877, 222, 506,
                         507, 577, 921, 153, 530, 736, 90, 718, 174, 216,
                         579, 581, 584, 569, 569, 182, 182, 856, 572, 573,
                         589, 590, 590, 739, 739, 508, 235, 142, 245, 241,
                         242, 509, 247, 243, 1, 612, 611, 611, 128, 503, 
                         769, 773, 607, 325, 605, 76, 74, 77, 604, 602,
                         754, 754, 448, 592, 737, 737, 56, 728, 482, 476,
                         700, 723, 366, 733, 752, 904, 687, 677, 753, 37,
                         888, 672, 163, 663, 853, 658, 162, 839, 837, 647,
                         775, 554, 190, 560, 555, 551, 252, 646, 296, 575,
                         576, 645, 270, 628, 883, 883, 614)
    }
     ##(738, 737) -> (737, 751)로 바꿈
    dist.imsi <- c()
    lll <- 0
    for(jj in 1:nrow(data$eventplace)){
      pts.coords <- c(data$eventplace$경도.Degree.[jj], data$eventplace$위도.Degree.[jj])
      
      lll <- lll+1
      if((jj%in%miho)){
        lll <-lll - 1
      }
      kk <- RCH_ID_manual[lll]
      #kk <- match(RCH_ID_manual[lll], shape@data$RCH_ID) #이게 부정확한 것 같아서 수정하자
      #kk <- match(data$eventplace$RCH_ID[jj], shape@data$RCH_ID) #이게 부정확한 것 같아서 수정하자
      #실제 RID랑 조금 잘 안맞는 것 같다
      #dist_upper <- apply(shape@data[,c(3,4)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
      #dist_lower <- apply(shape@data[,c(5,6)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
      #
      #kk <- sort(unique(union(which(dist_upper==min(c(dist_upper, dist_lower))), which(dist_lower==min(c(dist_upper, dist_lower))))))
      #
      #dist.imsiimsi <- c()
      #for(pp in 1:length(kk)){
      #  for(ll in 1:nrow(flow_network_original$lines[[kk[pp]]])){
      #    dist.compute <- sqrt(sum((flow_network_original$lines[[kk[pp]]][ll,]-pts.coords)^2))
      #    dist.imsiimsivec <- matrix(c(dist.compute, kk[pp], ll), nrow=1)
      #    dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
      #  }
      #}
      #dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
      #dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
      
      dist.imsiimsi <- c()
      for(ll in 1:nrow(flow_network_original$lines[[kk]])){
        dist.compute <- sqrt(sum((flow_network_original$lines[[kk]][ll,]-pts.coords)^2))
        if(lll==9){
          dist.compute <- sqrt(sum((flow_network_original$lines[[kk]][ll,]-pts.coords-c(0.01,0))^2))
        }
        dist.imsiimsivec <- matrix(c(dist.compute, kk, ll), nrow=1)
        dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
      }
      dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
      dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
    }
    colnames(dist.imsi) <- c("mindist", "rid", "segid")
    rownames(dist.imsi) <- data$eventplace$X
    
    if(RID=="original"){
      dist.imsi[9,3] <- 126
    }else if(RID=="corrected"){
      dist.imsi[9,3] <- 41
    }
    
    DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
    for(ii in 1:length(DistanceUpstream_obs)){
      DistanceUpstream_obs[ii] <- upstream(startseg = 614, endseg = dist.imsi[ii,2], startvert = 8, endvert = dist.imsi[ii,3], rivers=flow_network_original, net=TRUE)
    }
    
    #data$eventplace$RCH_ID <- shape@data$RCH_ID[RCH_ID_manual]
    pointdata <- cbind(data$eventplace[-miho,], rid_old=dist.imsi[-miho,2], upDist=DistanceUpstream_obs[-miho], locID=c(1:(nrow(data$eventplace)-length(miho))), netID=rep(1,(nrow(data$eventplace)-length(miho))), pid=setdiff(c(1:nrow(data$eventplace)), miho), shreve=shreve_obj[dist.imsi[-miho,2]], mindist=dist.imsi[-miho,1], rid=RCH_ID_manual, segid=dist.imsi[-miho,3])
  }else{
    pointdata <- cbind(data$eventplace[-miho,], rid_old=dist.imsi[-miho,2], upDist=DistanceUpstream_obs[-miho], locID=c(1:(nrow(data$eventplace)-length(miho))), netID=rep(1,(nrow(data$eventplace)-length(miho))), pid=setdiff(c(1:nrow(data$eventplace)), miho), shreve=shreve_obj[dist.imsi[-miho,2]], mindist=dist.imsi[-miho,1], rid=dist.imsi[-miho,2], segid=dist.imsi[-miho,3])
  }
  pointbbox <- rbind(range(pointcoords[-miho,1]), range(pointcoords[-miho,2])); colnames(pointbbox) <- c("min", "max"); rownames(pointbbox) <- c("coords.x1", "coords.x2")
  network.SSNPoints <- list()
  network.SSNPoints[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.obsdata, point.coords = pointcoords, point.data=pointdata, points.bbox=pointbbox, proj4string=shape@proj4string)
  obspoints <- new("SSNPoints", SSNPoints=network.SSNPoints, ID="Obs")
  #(3)predpoint
  network.line.coords.preddata <- data.frame(NetworkID=rep(1,length(miho)), SegmentID=dist.imsi[miho,2], DistanceUpstream=DistanceUpstream_obs[miho])
  rownames(network.line.coords.preddata) <- rownames(dist.imsi)[miho]
  pointcoords.preddata <- cbind(data$eventplace$경도.Degree.[miho], data$eventplace$위도.Degree.[miho]); colnames(pointcoords.preddata) <- c("coords.x1", "coords.x2")
  #pointdata에 ratio 추가해야
  pointdata.preddata <- cbind(data$eventplace[miho,], rid_old=dist.imsi[miho,2], upDist=DistanceUpstream_obs[miho], locID=c(1:length(miho)), netID=rep(1,length(miho)), pid=miho, shreve=shreve_obj[dist.imsi[miho,2]], mindist=dist.imsi[miho,1], rid=dist.imsi[miho,2], segid=dist.imsi[miho,3])
  pointbbox.preddata <- rbind(range(pointcoords[miho,1]), range(pointcoords[miho,2])); colnames(pointbbox.preddata) <- c("min", "max"); rownames(pointbbox.preddata) <- c("coords.x1", "coords.x2")
  network.SSNPoints.preddata <- list()
  network.SSNPoints.preddata[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.preddata, point.coords = pointcoords.preddata, point.data=pointdata.preddata, points.bbox=pointbbox.preddata, proj4string=shape@proj4string)
  if(length(miho)==0){
    #(3-1) predpt가 아무것도 없을 경우
    predpoints <- new("SSNPoints", SSNPoints=list(), ID=character(0))
  }else{
    predpoints <- new("SSNPoints", SSNPoints=network.SSNPoints.preddata, ID="Preds")
  }
  #(4)SSNpath: path에 해당
  SSNpath <- "NULL"
  #(5)data: data.frame의 형태이다 ( rid,upDist,    Length, netID가 들어가 있다)
  #data:Object of class "data.frame". 
  #The number of rows in data should equal the number of lines in the lines object. Row names correspond to SegmentID values
  shapedata <- shape@data
  shapedata <- cbind(shapedata, rid=c(1:nrow(shape@data)), upDist=DistanceUpstream_seg, Length=shape@data$Shape_Leng, netID=rep(1,nrow(shape@data)), shreve=shreve_obj )
  #(6)lines: segment polygon들을 담는다
  shapelines <- shape@lines
  #(7)bbox
  shapebbox <- shape@bbox; colnames(shapebbox) = c("min", "max"); rownames(shapebbox) <- c("x", "y")
  #(8)roj4string
  shapeproj4string <- shape@proj4string
  
  Geum_original <- new("SpatialStreamNetwork", network.line.coords=network.line.coords, obspoints=obspoints, predpoints=predpoints, path=SSNpath, data=shapedata, lines=shapelines, bbox=shapebbox, proj4string=shapeproj4string)
  
  #Object of class "SSNPoints" with 2 slots
  #
  #@ SSNPoints: List of SSNPoint objects with 5 slots
  #@ network.point.coords: object of class "data.frame". Row names
  #represent point identifiers (pid) stored in the point.data
  #data.frame.
  #$ NetworkID: factor identifying the NetworkID of that point
  #$ SegmentID: factor identifying the unique stream segment of that point
  #$ DistanceUpstream: numeric value representing the cumulative
  #distance from the network outlet, the most downstream point
  #on a network, to that point
  #@ point.coords: numeric matrix or "data.frame" with x- and y-
  #  coordinates (each row is a point); row names represent point
  #identifiers (pid) stored in the point.data data.frame.
  #@ point.data: object of class "data.frame"; the number of rows in
  #data should equal the number of points in the
  #network.point.coords object; row names are set to the pid
  #attribute.
  #@ points.bbox: Object of class "matrix"; see Spatial-class
  #@ proj4string: Object of class "CRS"; see CRS-class
  #@ ID: character string representing the name of the observation points
  
  
  #network.line <- shape_tinned@lines
  #bbox, proj4string
  
  return(Geum_original)
}


get_binaryIDs_usual <- function(mouth_node, shape){
  #함수설명
  #원래 금강 자료를 분석하기 위해 get_binaryIDs_stream 함수를 만들었는데
  #이 함수는 shape 파일에 RCH_ID, binaryID라는 변수가 없으면 사용할 수 없다
  #그래서 새로운 버전을 만들었다(일반화용?)
  #기준: stpca 패키지에 시뮬레이션 자료에 있는 shape들
  rtNEL1 <-readshpnw(shape, ELComputed=TRUE, longlat=TRUE) #길이 7
  igr1 <-nel2igraph(rtNEL1[[2]], rtNEL1[[3]]) #shape파일을 graph로 변환하였다
  
  edgelist <- rtNEL1[[3]]
  mouth_node_index <- edgelist[mouth_node,c(2,3)]
  mouth_node_final <- c()
  if(sum(edgelist[,c(2,3)]== mouth_node_index[1])==1){
    mouth_node_final <- c(mouth_node_final,mouth_node_index[1])
    key_ind <- mouth_node_index[2]
  }else if(sum(edgelist[,c(2,3)]== mouth_node_index[2])==1){
    mouth_node_final <- c(mouth_node_final,mouth_node_index[2])
    key_ind <- mouth_node_index[1]
  }
  
}

#mouth_node=68; shape=shape_miho
#shape 파일에 RCH_ID,binaryID 팔요함
get_binaryIDs_stream <- function(mouth_node, shape, RCH_ID=NULL){
  if(is.null(RCH_ID)){
    RCH_ID <- shape$RCH_ID
  }
  
  rtNEL1 <-readshpnw(shape, ELComputed=TRUE, longlat=TRUE)
  igr1 <-nel2igraph(rtNEL1[[2]], rtNEL1[[3]])
  
  edgelist <- rtNEL1[[3]]
  mouth_node_index <- edgelist[mouth_node,c(2,3)]
  mouth_node_final <- c()
  if(sum(edgelist[,c(2,3)]== mouth_node_index[1])==1){
    mouth_node_final <- c(mouth_node_final,mouth_node_index[1])
    key_ind <- mouth_node_index[2]
  }else if(sum(edgelist[,c(2,3)]== mouth_node_index[2])==1){
    mouth_node_final <- c(mouth_node_final,mouth_node_index[2])
    key_ind <- mouth_node_index[1]
  }
  
  
  #result_binaryIDs <- data.frame(rid=c(1:length(shape$RCH_ID)), binaryID=rep("", length(shape$RCH_ID)))
  result_binaryIDs <- data.frame(rid=c(1:length(RCH_ID)), binaryID=rep("", length(RCH_ID)))
  result_binaryIDs$binaryID <- as.character(result_binaryIDs$binaryID)
  idx <- 1
  result_binaryIDs[mouth_node,2] <- as.character("1")
  key_RID_ind <- mouth_node
  pointsin <- c(1:nrow(result_binaryIDs))
  pointsin <- setdiff(pointsin, key_RID_ind)
  while(length(key_RID_ind)!=0){
    RID_candidate <- which(edgelist[,3]==key_ind[1])
    if(length(RID_candidate)==0){
      pointsin <- setdiff(pointsin, key_RID_ind[1])
    }else{
      if(result_binaryIDs[RID_candidate[1],2]==""){
        result_binaryIDs[RID_candidate[1],2] <- as.character(paste(result_binaryIDs[key_RID_ind[1],2],"0", sep=""))
        key_RID_ind <- c(key_RID_ind, RID_candidate[1])
        pointsin <- setdiff(pointsin, RID_candidate[1])
        if(length(RID_candidate)==2 & result_binaryIDs[RID_candidate[2],2]==""){
          result_binaryIDs[RID_candidate[2],2] <- as.character(paste(result_binaryIDs[key_RID_ind[1],2],"1", sep=""))
          key_RID_ind <- c(key_RID_ind, RID_candidate[2])
          pointsin <- setdiff(pointsin, RID_candidate[2])
        }
      }else if(length(RID_candidate)==2 & result_binaryIDs[RID_candidate[2],2]==""){
        result_binaryIDs[RID_candidate[2],2] <- as.character(paste(result_binaryIDs[key_RID_ind[1],2],"0", sep=""))
        key_RID_ind <- c(key_RID_ind, RID_candidate[2])
        pointsin <- setdiff(pointsin, RID_candidate[2])
      }
    }
    if(length(key_RID_ind)!=0){
      key_RID_ind <- key_RID_ind[-1]
      if(length(key_RID_ind)!=0){
        key_ind <- edgelist[key_RID_ind[1],2]
      }
    }
    #print(key_RID_ind)
  }
  return(result_binaryIDs)
}

#get_adjacency_stream 돌리기 위해 필요한 함수
re_map_rid <- function(rid_vector, all_rid){
  # in cases where multiple networks are present, gaps could occur in the rid sequence
  # this function maps the rid vector in a simple way so that it agrees with the adjacency matrix
  mapped_rid <- vector("numeric", length = length(rid_vector))
  rng_rid    <- range(all_rid)
  new_key    <- 1:length(unique(all_rid))
  old_rid    <- sort(all_rid)
  for(i in 1:length(rid_vector)){
    mapped_rid[i] <- new_key[which(old_rid == rid_vector[i])]
  }
  mapped_rid
}

#get_adjacency와 정확하게 같으나 binaryIDs_obj를 get_binaryIDs_stream으로부터 얻어낸다
get_adjacency_stream <- function (binaryIDs_obj, netID = 1) 
{
  binaryIDs <- binaryIDs_obj
  #binaryIDs <- get_binaryIDs(ssn_directory, net = netID)
  bid <- binaryIDs[, 2]
  rid <- binaryIDs[, 1]
  nch <- nchar(bid)
  bid.list <- split(bid, nch)
  rid.list <- split(rid, nch)
  n.segments <- nrow(binaryIDs)
  xy <- matrix(ncol = 2, nrow = (n.segments - 1))
  counter <- 0
  for (j in 1:(length(bid.list) - 1)) {
    bid.dn <- bid.list[[j]]
    bid.up <- bid.list[[j + 1]]
    rid.dn <- rid.list[[j]]
    rid.up <- rid.list[[j + 1]]
    bid.up.sub.vec <- bid.up
    rid.up.sub.vec <- rid.up
    for (i in 1:length(bid.dn)) {
      current.dn.bid <- bid.dn[i]
      current.dn.rid <- rid.dn[i]
      inner.count <- 1
      number.upstream <- 0
      n.bid.up <- length(bid.up.sub.vec)
      crit <- FALSE
      while (!crit) {
        if (n.bid.up > 0) {
          current.up.bid <- bid.up.sub.vec[inner.count]
          connected <- substr(current.up.bid, 1, nchar(current.dn.bid)) == 
            current.dn.bid
          if (connected) {
            counter <- counter + 1
            number.upstream <- number.upstream + 1
            xy[counter, ] <- c(current.dn.rid, rid.up.sub.vec[inner.count])
            rid.up.sub.vec <- rid.up.sub.vec[-inner.count]
            bid.up.sub.vec <- bid.up.sub.vec[-inner.count]
          }
          if (!connected) 
            inner.count <- inner.count + 1
          crit <- (number.upstream == 2) | ((number.upstream +  inner.count - 1) == n.bid.up)
        }
        if (n.bid.up == 0) 
          crit <- TRUE
      }
    }
  }
  xy[, 1] <- re_map_rid(rid_vector = xy[, 1], all_rid = rid)
  xy[, 2] <- re_map_rid(rid_vector = xy[, 2], all_rid = rid)
  add.one <- min(xy) == 0
  list(adjacency = spam(list(j = (xy[, 1] + add.one), i = (xy[, 2] + add.one), rep(1, nrow(xy))), nrow = n.segments, 
                        ncol = n.segments), rid_bid = cbind(rid, bid))
}

#compute shreve stream order
compute_shreve <- function(adjacency){
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y자
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
    }
    #case 2: 1자
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  return(shreve)
}


#compute shreve stream order
compute_shreve_and_dist <- function(adjacency, example_network, type="equal", scalevec=c(0.2, 1.5), logdata=FALSE){
  #type="lengthprop": shreve_and_dist[first_stream] is defined by proportional to their lengths
  #type="equal": shreve_and_dist[first_stream] is defined by equal weight (== shreve dist style)
  #scalevec: 최종 weight를 scaling (default는 (0.2,1.5) 이다 (O'Donnell 방법 따라함))
  
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  shreve_and_dist <- rep(0, nrow(adj_mat))
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  if(type=="lengthprop"){
    shreve_and_dist[first_stream] <- example_network@data$Length[first_stream]
  }else if(type=="equal"){
    shreve_and_dist[first_stream] <- 1
  }else{
    stop("wrong types")
  }
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y자
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
      shreve_and_dist[next_stream1[i]] <- sum(shreve_and_dist[which(adj_mat[,next_stream1[i]]!=0)]) + example_network@data$Length[next_stream1[i]]
    }
    #case 2: 1자
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
      shreve_and_dist[next_stream2[j]] <- sum(shreve_and_dist[which(adj_mat[,next_stream2[j]]!=0)]) + example_network@data$Length[next_stream2[j]]
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  #scaling
  #weight_vec_candidate2 <- scales::rescale(log(weight_vec_candidate2$distweight), to=c(0.2, 1.5))
  if(logdata==TRUE){
    shreve_and_dist <- scales::rescale(log(sqrt(shreve_and_dist)), to=c(scalevec[1], scalevec[2]))
    #shreve_and_dist <- scales::rescale(log(shreve_and_dist), to=c(scalevec[1], scalevec[2]))
  }else{
    shreve_and_dist <- scales::rescale((sqrt(shreve_and_dist)), to=c(scalevec[1], scalevec[2]))
  }
  return(list(shreve=shreve, distweight=shreve_and_dist))
}


compute_shreve_weighted <- function(adjacency, weight){
  #compute_shreve_and_dist함수는 example_network의 Length를 인자로 받는데
  #그것 대신 weight object를 따로 받는 식으로 일반화하였다 (실제 자료분석 코드에 쓰지는 않았다)
  #weighted version으로 확장
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  shreve_and_dist <- rep(0, nrow(adj_mat))
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  shreve_and_dist[first_stream] <- weight[first_stream]
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y자
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
      shreve_and_dist[next_stream1[i]] <- sum(shreve_and_dist[which(adj_mat[,next_stream1[i]]!=0)]) + weight[next_stream1[i]]
    }
    #case 2: 1자
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
      shreve_and_dist[next_stream2[j]] <- sum(shreve_and_dist[which(adj_mat[,next_stream2[j]]!=0)]) + weight[next_stream2[j]]
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  return(list(shreve=shreve, distweight=shreve_and_dist))
}


########################################
##Graphics
########################################

source('~/Dropbox/R files/StreamFlow/sources/iwanthue.R', chdir = TRUE)

plot.SpatialStreamNetwork.SC <- 
  function (x, VariableName = NULL, color.palette = NULL, nclasses = NULL, 
            breaktype = "quantile", brks = NULL, PredPointsID = NULL, 
            add = FALSE, addWithLegend = FALSE, lwdLineCol = NULL, lwdLineEx = 1, 
            lineCol = "black", I_lines=NULL, pointsin=NULL, ...) 
  {
    #my addition
    if(is.null(pointsin)){
      pointsin <- c(1:length(I_lines))
    }
    if (missing(lwdLineEx)) 
      lwdLineEx <- 1
    if (missing(lwdLineCol)) {
      x@data$lineWidth <- rep(1, nrow(x@data))
      lwdLineCol <- "lineWidth"
    }
    if (is.null(as.list(match.call()[-1])$pch)) {
      plch = 19
    }
    else plch <- as.list(match.call()[-1])$pch
    if (is.null(as.list(match.call()[-1])$cex)) {
      chex = 1
    }
    else chex <- as.list(match.call()[-1])$cex
    if (is.null(as.list(match.call()[-1])$col)) {
      colr = "black"
    }
    else colr <- as.list(match.call()[-1])$col
    par.orig <- par(no.readonly = TRUE)
    if (!is.null(PredPointsID)) {
      #prediction Point ID를 입력했을 때
      for (i in 1:length(x@predpoints@ID)) {
        if (x@predpoints@ID[i] == PredPointsID) {
          if (add == FALSE & addWithLegend == FALSE) {
            plot(x@bbox[1, ], x@bbox[2, ], type = "n",   ...)
            for (j in 1:length(x@lines)) for (k in 1:length(x@lines[[j]])) if (is.null(lwdLineCol)) 
              lines((x@lines[[j]]@Lines[[k]]@coords), col = lineCol,  ...)
            else lines(x@lines[[j]]@Lines[[k]]@coords, 
                       lwd = lwdLineEx * x@data[i, lwdLineCol], 
                       col = lineCol, ...)
          }
          if (add == TRUE) {
            par(new = TRUE)
            plot(x@bbox[1, ], x@bbox[2, ], type = "n", 
                 bty = "n", xlab = "", ylab = "", ...)
          }
          if (addWithLegend == TRUE) {
            par(new = TRUE)
            layout(matrix(1:2, nrow = 1), widths = c(4, 
                                                     1))
            par(mar = c(5, 5, 3, 0))
            par(mfg = c(1, 1))
            plot(x@bbox[1, ], x@bbox[2, ], type = "n", 
                 bty = "n", xlab = "", ylab = "", ...)
          }
          points(x@predpoints@SSNPoints[[i]]@point.coords, 
                 pch = plch, cex = chex, col = colr)
        }
      }
      par(par.orig)
    }else if (is.null(VariableName)) {
      #variable name이 null인 경우
      plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
      for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
        lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol, 
              ...)
      else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * 
                   x@data[i, lwdLineCol], col = lineCol, ...)
      points(x@obspoints@SSNPoints[[1]]@point.coords, pch = plch, 
             cex = chex, col = colr)
      par(par.orig)
    }else {
      #보통의 경우는 여기에 해당하는 듯
      layout(matrix(1:2, nrow = 1), widths = c(4, 1))
      par(mar = c(5, 5, 3, 0))
      plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
      #for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
      #  lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol,  ...)
      #else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * 
      #             x@data[i, lwdLineCol], col = lineCol, ...)
      
      #define color palette
      #color.palette = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      #set.seed(1)
      #col.stream <- sample(color.palette, length(I_lines))
      #pie(rep(1,length(I_lines)), col=col.stream)
      col.stream <- iwanthue(n=length(I_lines), hmin=0, hmax=360, cmin=30, cmax=80, lmin=35, lmax=80, plot=FALSE, random=1) 
      
      for (i in 1:length(pointsin)){
        for (j in 1:length(I_lines[[pointsin[i]]])){
          lines(I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[c(I_lines[[pointsin[i]]][[j]]@range[1]:I_lines[[pointsin[i]]][[j]]@range[2]),], col=col.stream[i], lwd=lwdLineEx*x@obspoints@SSNPoints[[1]]@point.data[i, lwdLineCol])
          if(I_lines[[pointsin[i]]][[j]]@range[1]!=1){
            mid_pts <- (I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[i]][[j]]@range[1],] + I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[(I_lines[[pointsin[i]]][[j]]@range[1]-1),])/2
            lines(rbind(mid_pts, I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[pointsin[i]]][[j]]@range[1],]), col=col.stream[i], lwd=lwdLineEx*x@obspoints@SSNPoints[[1]]@point.data[i, lwdLineCol])
          }
          if(I_lines[[pointsin[i]]][[j]]@range[2]!=nrow(I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords)){
            mid_pts <- (I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[pointsin[i]]][[j]]@range[2],] + I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[(I_lines[[pointsin[i]]][[j]]@range[2]+1),])/2
            lines(rbind(mid_pts, I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[pointsin[i]]][[j]]@range[2],]), col=col.stream[i], lwd=lwdLineEx*x@obspoints@SSNPoints[[1]]@point.data[i, lwdLineCol])
          }
        }
      }
      
      data <- x@obspoints@SSNPoints[[1]]@point.data
      if (is.null(nclasses)) 
        nclasses <- 10
      lower.breaks <- matrix(0, nrow = nclasses, ncol = 1)
      upper.breaks <- matrix(0, nrow = nclasses, ncol = 1)
      if (breaktype == "quantile") {
        brks <- quantile(data[, VariableName], probs = (1:(nclasses -    1))/nclasses, na.rm = T)
        lower.breaks <- c(min(data[, VariableName], na.rm = T),   brks)
        upper.breaks <- c(brks, max(data[, VariableName],     na.rm = T))
      }
      if (breaktype == "even") {
        brks <- min(data[, VariableName]) + (max(data[, VariableName]) - min(data[, VariableName])) * (1:(nclasses - 1))/nclasses
        lower.breaks <- c(min(data[, VariableName], na.rm = T),   brks)
        upper.breaks <- c(brks, max(data[, VariableName],  na.rm = T))
      }
      if (breaktype == "user") {
        if (is.null(brks)) 
          return("Must specify brks if breaktype = user")
        minD <- min(data[, VariableName], na.rm = TRUE)
        maxD <- max(data[, VariableName], na.rm = TRUE)
        brks <- as.vector(unlist(brks))
        if (minD < min(brks)) 
          brks <- c(brks, minD)
        if (maxD > max(brks)) 
          brks <- c(brks, maxD)
        brks <- sort(unique(unlist(brks)))
        nclasses <- length(brks) - 1
        lower.breaks <- brks[1:nclasses]
        upper.breaks <- brks[2:(nclasses + 1)]
      }
      if (length(color.palette) == 0) 
        color.palette <- rainbow(nclasses, start = 0.66, 
                                 end = 0.99)
      for (j in 1:nclasses) {
        jmax <- upper.breaks[j]
        jmin <- lower.breaks[j]
        indj <- data[, VariableName] >= jmin & data[, VariableName] <=   jmax
        points(x@obspoints@SSNPoints[[1]]@point.coords[indj, , drop = F], col = color.palette[j], pch = plch,  cex = chex)
        #my addition
        text(x@obspoints@SSNPoints[[1]]@point.coords[indj, 1 , drop = F]+0.005, x@obspoints@SSNPoints[[1]]@point.coords[indj, 2 , drop = F]+0.005, pointsin[which(indj)])
      }
      dec.dig <- 2
      left <- as.character(as.numeric(as.integer(lower.breaks * 
                                                   10^dec.dig))/10^dec.dig)
      rght <- as.character(as.numeric(as.integer(upper.breaks * 
                                                   10^dec.dig))/10^dec.dig)
      leglabs <- paste(left, "to", rght)
      par(mar = c(0, 0, 0, 0))
      plot(c(0, 0), c(1, 1), type = "n", xaxt = "n", yaxt = "n", 
           xlab = "", ylab = "", bty = "n")
      legend(x = -1, y = 1.1, legend = leglabs, bty = "n", 
             pch = rep(plch, times = length(leglabs)), col = color.palette, 
             cex = 0.8)
      par(par.orig)
      return(invisible(data.frame(lower.breaks = lower.breaks, 
                                  upper.breaks = upper.breaks)))
    }
  }
















########################################
##CODES not directly related to the SSN and smnet package
########################################

find_orthogonal <- function(x,y){
  #x: 정사영당하는 벡터
  #y: 정사영하려는 점
  cc = crossprod(x,y)/crossprod(x)
  return(c(x[1]*cc,x[2]*cc))
}


forward_line <- function(line_remove, line_nbrs, lengths_remove, lengths_nbrs, adj_matrix, remove, nbrs, type="b2"){
  if(type=="a1" | type=="b1"){
    #case a1
    #이 경우에는 상류가 line_remove, 하류가 line_nbrs가 됨
    newpt <- (line_remove[1,] + find_orthogonal(x=(line_nbrs[2,]- line_remove[1,]), y=(line_remove[2,] - line_remove[1,]) ) )
    
    coarse_line <- rbind(line_remove[1,], line_nbrs[2,])
    detail_line <- rbind(line_remove[2,], newpt)
    
    coarse_lengths <- lengths_nbrs + lengths_remove
    detail_lengths <- lengths_remove - coarse_lengths
    
    #adjacency matrix change
    adj_matrix_imsi <- adj_matrix
    
    if(type=="b1"){
      #nbrs의 하류와 연결된 노드, 그리고 detail의 상류와 연결된 노드를 연결시켜줘야 함
      
      adj_matrix_imsi[nbrs,which(adj_matrix[remove,]%in%c(2))] <- adj_matrix[remove,which(adj_matrix[remove,]%in%c(2))]
      adj_matrix_imsi[which(adj_matrix[,remove]%in%c(2)),nbrs] <- adj_matrix[remove,which(adj_matrix[,remove]%in%c(2))]
    }
    #detail인 부분들은 모두 음수로 바꾸자
    adj_matrix_imsi[remove,] <- -abs(adj_matrix_imsi[remove,])
    adj_matrix_imsi[,remove] <- -abs(adj_matrix_imsi[,remove])
     
  }else if(type=="a2" | type=="b2"){
    #case b2
    #이 경우에는 상류가 line_nbrs, 하류가 line_remove가 됨
    #test_a <- rbind(line_nbrs[1,], line_remove[2,])
    #test_b <- coeff_lines[[remove]][1,]
    
    newpt <- (line_nbrs[1,] + find_orthogonal(x=(line_remove[2,] - line_nbrs[1,]) , y= (line_nbrs[2,] - line_nbrs[1,]) ) )
    
    coarse_line <- rbind(line_nbrs[1,], line_remove[2,])
    detail_line <- rbind(line_remove[1,], newpt)
    
    coarse_lengths <- lengths_nbrs + lengths_remove
    detail_lengths <- lengths_remove - coarse_lengths
    
    
    #adjacency matrix change
    adj_matrix_imsi <- adj_matrix
    
    if(type=="b2"){
      #nbrs의 상류와 연결된 노드, 그리고 detail의 하류와 연결된 노드를 연결시켜줘야 함
      adj_matrix_imsi[nbrs,which(adj_matrix[remove,]%in%c(3,4))] <- adj_matrix[remove,which(adj_matrix[remove,]%in%c(3,4))]
      adj_matrix_imsi[which(adj_matrix[,remove]%in%c(2,4)),nbrs] <- adj_matrix[remove,which(adj_matrix[,remove]%in%c(2,4))]
    }
    #detail인 부분들은 모두 음수로 바꾸자
    adj_matrix_imsi[remove,] <- -abs(adj_matrix_imsi[remove,])
    adj_matrix_imsi[,remove] <- -abs(adj_matrix_imsi[,remove])
  }
  return(list(coarse_line=coarse_line, detail_line=detail_line, coarse_lengths=coarse_lengths, detail_lengths=detail_lengths, adj_matrix=adj_matrix_imsi, type=type))
}




forward_line_Y <- function(line_remove, line_coremove, line_nbrs, lengths_remove, lengths_coremove, lengths_nbrs, adj_matrix, remove, coremove, nbrs){
  #Y자 형태 제거위해 특별히 제작함
  
  
  #plot(rbind(line_remove, line_coremove, line_nbrs), ylim=c(36.57288-(0.04621/2), 36.60227+(0.04621/2)))
  #segments(x0=line_remove[1,1], y0=line_remove[1,2], x1=line_remove[2,1], y1=line_remove[2,2], col="red")
  #segments(x0=line_coremove[1,1], y0=line_coremove[1,2], x1=line_coremove[2,1], y1=line_coremove[2,2], col="green")
  #segments(x0=line_nbrs[1,1], y0=line_nbrs[1,2], x1=line_nbrs[2,1], y1=line_nbrs[2,2], col="blue")
  
  #1차 제거
  newpt_imsi <- -(-(line_remove[1,]-line_remove[2,])-(line_coremove[1,]-line_coremove[2,]))+line_remove[2,]
  
  result_line_remove <- -line_remove
  result_line_coremove <- rbind(newpt_imsi, line_coremove[2,])
  
  #adjacency matrix update
  adj_matrix[remove, ] <- -abs(adj_matrix[remove,])
  adj_matrix[,remove] <- -abs(adj_matrix[,remove])
  
  #update length
  lengths_remove <- lengths_remove - lengths_coremove
  lengths_coremove <- lengths_remove + lengths_coremove
 
  #2차 제거
  result_new <- forward_line(line_remove=result_line_coremove, line_nbrs=line_nbrs, lengths_remove=lengths_coremove, lengths_nbrs=lengths_nbrs, adj_matrix, remove=coremove, nbrs=nbrs, type="b1")
  
  type <- c("c1")
  
  #points(newpt_imsi[1], newpt_imsi[2], col=5)
  
  #test_a <- rbind(coeff_lines[[nbrs]][1,], coeff_lines[[remove]][2,])
  #test_b <- coeff_lines[[remove]][1,]
  #project(test_b , (test_a[2,] -test_a[1,]))
  #
  #coeff_lines[[nbrs]][1,] + find_orthogonal(x=(test_a[2,] -test_a[1,]) , y= (coeff_lines[[nbrs]][2,]-coeff_lines[[nbrs]][1,]) )
  
  #plot(rbind(test_a,test_b), ylim=c(36.581, 36.609), xlim=c(127.292, 127.310))
  #segments(x0=test_a[1,1], y0=test_a[1,2], x1=test_a[2,1], y1=test_a[2,2])
  #segments(x0=coeff_lines[[nbrs]][1,1], y0=coeff_lines[[nbrs]][1,2], x1=coeff_lines[[nbrs]][2,1], y1=coeff_lines[[nbrs]][2,2], col="blue")
  #segments(x0=coeff_lines[[remove]][1,1], y0=coeff_lines[[remove]][1,2], x1=coeff_lines[[remove]][2,1], y1=coeff_lines[[remove]][2,2], col="red")
  #
  #points(127.30154,  36.58159, col="green", pch=16)
  #points(127.3016, 36.60801, col="green", pch=16)
  

  return(list(coarse_line=result_new$coarse_line, detail_line=result_new$detail_line, coarse_lengths=result_new$coarse_lengths, detail_lengths=result_new$detail_lengths, adj_matrix=result_new$adj_matrix, type=type))
}


########################################
##CODES related to graphics
########################################
#reference: https://stackoverflow.com/questions/20127282/r-color-scatterplot-points-by-z-value-with-legend
scatter_fill <- function (x, y, z,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),zlim=c(min(z),max(z)),
                          nlevels = 121, plot.title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", 
                          yaxs = "i", las = 1, y.axes.label=TRUE,
                          axes = TRUE, y.axes=TRUE, frame.plot = axes, smallplot=c(0.85,0.9,0.25,0.75), coldefault=colorRampPalette(c("cyan", "green", "yellow", "red", "black")), plot.legend=TRUE, ...) 
{
  old.par <- par(no.readonly = TRUE)
  bigpar <- old.par$plt
  
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #on.exit(par(par.orig))
  #w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  w <- (3 + mar.orig[2L]) * par("csi") * 2.24
  #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  
  # choose colors to interpolate
  levels <- seq(zlim[1],zlim[2],length.out = nlevels)
  #col <- colorRampPalette(c("red","yellow","dark green"))(nlevels) 
  #following ODonnell's Approach
  col<-coldefault(nlevels)
  colz <- col[cut(z,nlevels)]  
  
  par(bigpar)
  mar <- mar.orig
  #mar[4L] <- 1
  #mar[4L] <- 0.1
  par(mar = mar)
  
  par(las = 0)
  # points
  plot(x,y,type = "n",xaxt='n',yaxt='n',xlab="",ylab="",xlim=xlim,ylim=ylim,bty="n")
  points(x,y,col = colz,xaxt='n',yaxt='n',xlab="",ylab="",bty="n",...)
  
  ## options to make mapping more customizable
  
  if (missing(plot.axes)) {
    if (axes & y.axes==TRUE) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1, labels = FALSE)
      Axis(y, side = 2, labels = FALSE)
    }else if (axes & y.axes==FALSE & y.axes.label==TRUE) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      #Axis(y, side = 2)
    }else if(axes & y.axes==FALSE & y.axes.label==FALSE){
      title(main = "", xlab = "", ylab = "")
      Axis(x, side=1, labels=FALSE, tick=FALSE)
      Axis(y, side=2, labels=FALSE, tick=FALSE)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  
  if(plot.legend==TRUE){
    #par(las = las)
    par(las = 0)
    mar <- mar.orig
    #mar[4L] <- mar[2L]
    #mar[2L] <- 1
    #mar[2L] <- 0.1
    par(mar = mar)
    #   
    #plot.new()
    if(!is.null(smallplot)){
      par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    }
    
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    
    par(las = 2)
    rect(0, levels[-length(levels)], 1, levels[-1L],col=col,border=col) 
    if (missing(key.axes)) {if (axes){axis(4)}}
    else key.axes
    box()
    if (!missing(key.title)) 
      key.title
  }
  
  par(new = TRUE, pty = "m", plt = bigpar, err = -1)
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  par(bigpar)
  mar <- mar.orig
  #mar[4L] <- 1
  #mar[4L] <- 0.1
  par(mar = mar)
  
  par(las = 0)
  plot(x,y,type = "n",xaxt='n',yaxt='n',xlab="",ylab="",xlim=xlim,ylim=ylim,bty="n")
  
  invisible()
}
