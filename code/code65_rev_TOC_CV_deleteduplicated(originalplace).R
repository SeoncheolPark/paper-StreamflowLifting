
deleted_ind_cand <- c(96:114)
for(jkl in 1:length(deleted_ind_cand)){
  library(ggmap)
  register_google(key='AIzaSyAtGCNC-GzkyMQ7ocNYC_G_eW7jzbBJ8Lg')
  library(SpatioTemporal)
  library(plotrix)
  library(maps)
  library(smnet)
  library(SSN)
  library(adlift); library(nlt); library(liftLRD); library(CNLTreg); library(CNLTtsa); library(CliftLRD); library(stringr)
  library(mgcv); library(spam); library(RgoogleMaps)
  library(lubridate)
  library(scales)
  
  ########################################
  #(1) Print Geum-River data
  ########################################
  source('~/Dropbox/R files/StreamFlow/sources/source.R', chdir = TRUE)
  source('~/Dropbox/R files/StreamFlow/sources/source_S.R', chdir = TRUE)
  source('~/Dropbox/R files/StreamFlow/sources/source_S_nlt.R', chdir = TRUE)
  source('~/Dropbox/(논문들)/(논문)River Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/paper_functions.r', chdir = TRUE)
  mouth_node <- 614 #set the mouth of the given network
  
  #길이 112개의 shape
  shape_miho <- readOGR("~/Dropbox/R files/StreamFlow/River/Geum.shp")
  
  #basic plotting
  ##MyRivernetwork: class "rivernetwork"
  MyRivernetwork <- line2network(path="~/Dropbox/R files/StreamFlow/River/", layer="Geum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
  MyRivernetwork <- setmouth(seg=614, vert=8, rivers=MyRivernetwork)
  #If you want to check river network,
  #par(mfrow=c(1,1))
  #plot network
  #plot(MyRivernetwork, segmentnum = FALSE)
  #plot(MyRivernetwork, xlim=c(126.4, 127.0), ylim=c(36.0, 36.2))
  
  #get_binaryIDs_stream: source에 만들어놓은 함수: 뒤의 get_adjacency_stream 실행 위함
  binaryIDs_obj <- get_binaryIDs_stream(mouth_node=614 , shape_miho)
  #수동으로 adjacency를 얻도록 만들어진 함수
  adjacency <- get_adjacency_stream(binaryIDs_obj)
  
  #compute_shreve: shreve 거리를 재도록 직접 제작한 함수
  shreve_order <- compute_shreve(adjacency)
  
  #generate SSN obj
  mf04 <- importSSN_stream(shreve_obj = shreve_order, location="Full", multipleplaces=TRUE) #데이터의 크기: 123개
  #check basics
  dim(mf04@obspoints@SSNPoints[[1]]@point.data) #127 38
  length(unique(mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID)) #114
  which(duplicated(mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID)) #13
  mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID[which(duplicated(mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID))] #3개짜리 하나 있음
  show_weights(mf04, adjacency)
  
  ##(Dec 9, 2020): duplicated pt 제거
  mf04_old <- mf04
  mf04@obspoints@SSNPoints[[1]]@network.point.coords
  which(duplicated(mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID))
  mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID
  mf04@obspoints@SSNPoints[[1]]@point.data$rid
  #"보청천4"   877(rid)
  #"회인천"    577(rid)
  #"주원천"    573(잘못된듯?)
  #"대청댐"    153(rid)
  #"품곡천"    530(rid)
  #"대청"      750(rid,  일부러 바꾼듯)
  #"월하천"   448(rid)
  #"미호천6-1" 592(rid)
  #"연기"      738(rid, 737으로 봐야할수도, 일부러 바꾼듯) -> 737
  #"금남"      737(rid) -> 751
  #"논산천2"   551(rid)
  #"방축천"    252(rid)
  
  mf04@obspoints@SSNPoints[[1]]@point.data$rid== mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID
  which(mf04@obspoints@SSNPoints[[1]]@point.data$rid != mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID)
  rownames(mf04@obspoints@SSNPoints[[1]]@network.point.coords)[which(mf04@obspoints@SSNPoints[[1]]@point.data$rid != mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID)]
  ########################################
  #Load pollutant dataset
  ########################################
  #data <- readRDS("~/Dropbox/Data/RIverFlow/금강/Geum(miho).RDS")
  #data <- readRDS("~/Dropbox/Data/RIverFlow/금강/Geum(miho)(extended).RDS")
  data <- readRDS("~/Dropbox/Data/RiverFlow/금강/Geum(extended).RDS")
  library(extrafont)
  #font_import()
  par(family="NanumGothic")
  
  #case (1): whole
  data$normaldata_TOC <- data$normaldata_TOC["2011-12/2017-11",]
  #case (2): winter
  #data$normaldata_TOC <- data$normaldata_TOC[which(month(data$normaldata_TOC)%in%c(12,1,2)),]
  #case (3): summer
  #data$normaldata_TOC <- data$normaldata_TOC[which(month(data$normaldata_TOC)%in%c(6,7,8)),]
  
  #delete unused point
  #안쓰는 점들 제거 (삼가천, 석천은 관찰값이 없는 장소)
  remove.cand <- which(colnames(data$normaldata_TOC)%in%c("삼가천", "석천"))
  data$normaldata_TOC <- data$normaldata_TOC[,-remove.cand ]
  
  #(추가) TP 자료에 대해서는 data$normaldata_TP==0인 자료 1개를 제거해야
  #나중에 로그변환을 취했을 때 에러가 나지 않는다
  for(ii in 1:ncol(data$normaldata_TOC)){
    if(length(which(data$normaldata_TOC[,ii]==0))!=0){
      data$normaldata_TOC[which(data$normaldata_TOC[,ii]==0),ii] <- NA
    }
  }
  
  normaldata_TOC <- data$normaldata_TOC
  
  #check basics
  data$eventplace$X[which(is.na(match(data$eventplace$X, colnames(data$normaldata_TOC))))]
  dim(data$eventplace)
  dim(data$normaldata_TOC)
  
  #If you want plot time series,
  #시계열 플롯
  #for(i in 1:ncol(data$normaldata_TN)){
  #  if(sum(is.na(data$normaldata_TN[,i]))!=nrow(data$normaldata_TN)){
  #    plot(log(as.numeric(data$normaldata_TN[,i])), main=colnames(data$normaldata_TN)[i], ylab="log(Nitrogen)")
  #  }
  #}
  
  #check matching indexes
  #(SSN obj의 장소 순서, 시계열 장소 순서 매치)
  match(mf04@obspoints@SSNPoints[[1]]@point.data$X, colnames(data$normaldata_TOC))
  #(시계열 장소 순서, SSN obj의 장소 순서 매치)
  match(colnames(data$normaldata_TOC), mf04@obspoints@SSNPoints[[1]]@point.data$X)
  #obsdata에 맞춰 정렬하였다
  #data$normaldata <- data$normaldata[,match(mf04@obspoints@SSNPoints[[1]]@point.data$X, colnames(data$normaldata))]
  
  #내 스스로 리버 network와 관찰 장소 위치를 출력하는 함수들
  plot(mf04@obspoints@SSNPoints[[1]]@point.coords[,1], mf04@obspoints@SSNPoints[[1]]@point.coords[,2], type="n")
  plot(MyRivernetwork, segmentnum = FALSE)
  text(x=mf04@obspoints@SSNPoints[[1]]@point.coords[,1], y=mf04@obspoints@SSNPoints[[1]]@point.coords[,2], labels=c(1:nrow(mf04@obspoints@SSNPoints[[1]]@point.coords)))
  
  #TN: 총질소
  #temp: 수온
  #ph: 수소이온농도
  #elec: 전기전도도
  #O2: 용존산소
  #BOD: BOD
  #COD: COD
  #susp: 부유물질
  #TP: 총인
  #TOC: 총유기탄소
  #flow: 유량
  
  eventplace2 <- data$eventplace[which(!is.na(match(data$eventplace$X, colnames(data$normaldata_TOC)))),] #측정값이 하나도 없는 station은 제거된다
  eventplace2 <- eventplace2[match(eventplace2$X, colnames(data$normaldata_TOC)),] #data$normaldata_TN의 순서에 맞춰 측정소 자료를 정렬한다
  row.names(eventplace2) <- eventplace2$X
  time_index <- as.character(time(data$normaldata_TOC["2011-12/2017-11",]))
  normaldata2 <- as.matrix(data$normaldata_TOC["2011-12/2017-11",])
  row.names(normaldata2) <- time_index
  stlistdata <- as.matrix(data$normaldata_temp["2011-12/2017-11"])
  row.names(stlistdata) <- time_index
  
  na.mismatch <- c()
  for(i in 1:ncol(normaldata2)){
    na.mismatch <- c(na.mismatch, sum(is.na(match(which(!is.na(normaldata2[,i])), which(!is.na(stlistdata[,i]))))) ) #??지금 당장 필요한 정보는 아님, 나중에 회귀분석 시 필요
  }
  #조천 2012년 9월 4일
  eventplace2 <- eventplace2[match(rownames(eventplace2), colnames(normaldata2)),]
  colnames(eventplace2)[1] <- "ID"
  colnames(eventplace2)[16] <- "long"
  colnames(eventplace2)[15] <- "lat"
  eventplace2 <- cbind(eventplace2, 
                       #TN=apply(as.matrix(data$normaldata_TN["2012/2017"]), 2, mean, na.rm=T)#, #TN: 총질소
                       #temp=apply(as.matrix(data$normaldata_TN["2012/2017"]), 2, mean, na.rm=T), #temp: 수온
                       #ph=apply(as.matrix(data$normaldata_ph["2012/2017"]), 2, mean, na.rm=T), #ph: 수소이온농도
                       #elec=apply(as.matrix(data$normaldata_elec["2012/2017"]), 2, mean, na.rm=T), #elec: 전기전도도
                       #O2=apply(as.matrix(data$normaldata_O2["2012/2017"]), 2, mean, na.rm=T), #O2: 용존산소
                       #BOD=apply(as.matrix(data$normaldata_BOD["2012/2017"]), 2, mean, na.rm=T), #BOD: BOD
                       #COD=apply(as.matrix(data$normaldata_COD["2012/2017"]), 2, mean, na.rm=T), #COD: COD
                       #susp=apply(as.matrix(data$normaldata_susp["2012/2017"]), 2, mean, na.rm=T), #susp: 부유물질
                       #TP=apply(as.matrix(data$normaldata_TP["2012/2017"]), 2, mean, na.rm=T)#, TP: 총인
                       TOC=log(apply(as.matrix(data$normaldata_TOC["2011-12/2017-11"]), 2, mean, na.rm=T)) #, TOC: 총유기탄소
                       #apply(as.matrix(data$normaldata_flow["2012/2017"]), 2, mean, na.rm=T) #flow: 유량 (NaN도 있다)
  ) 
  eventplace2 <- cbind(eventplace2, x=eventplace2$long, y=eventplace2$lat, type=ifelse(apply(normaldata2, 2, function(x) sum(!is.na(x))) > 100, "large", "small"))
  ########################################
  ##Start data analysis
  ########################################
  example_network <- mf04
  
  #RID를 0부터 끝까지 정렬하는 업데이트
  for(i in 1:length(example_network@lines)){
    example_network@lines[[i]]@ID <- as.character(i-1)
  }
  #Obs pt의 coord를 rid index에 맞춰 정렬하는 업데이트
  for(i in 1:nrow(example_network@obspoints@SSNPoints[[1]]@point.coords)){
    lineID <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[i]
    segID <- example_network@obspoints@SSNPoints[[1]]@point.data$segid[i]
    example_network@obspoints@SSNPoints[[1]]@point.coords[i,] <- example_network@lines[[lineID]]@Lines[[1]]@coords[segID,]
  }
  
  data$normaldata_TOC <- data$normaldata_TOC[,match(example_network@obspoints@SSNPoints[[1]]@point.data$X, colnames(data$normaldata_TOC))]
  eventplace2 <- eventplace2[match(example_network@obspoints@SSNPoints[[1]]@point.data$X, eventplace2$ID), ]
  
    ########################################
  #flexible smoothing
  #(O'Donnell et al., 2014)
  ########################################
  # vector of proportional flow weights 
  # associated with each river segment 
  #realweights<-read.table("realweights.txt") #없음
  #realweights.RData : vector of flow weights describing proportional flow contribution to nearest downstream neighbour on the River Tweed.
  
  # load adjacency matrix associated with each 
  # of the 298 stream segments of the river tweed
  #(관측 장소가 아닌 298개 stream segment들의 adjacency를 나타낸 것이다)
  adjacency_old <- adjacency 
  #adjacency<-read.table("adjacency.txt")
  #adjacency.RData:  298 by 298 sparse matrix defining neighbourhood structure of river network - if adjacency[i, j] == 1 then segment i flows into segment j.
  ##########나의 해결책
  adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))
  
  # load river tween data: 3 col dataframe:
  # dates, nitrate concs (mg/l) and 
  # location index (corresponding row number of adjacency matrix)
  #TweedData<-read.table("TweedData.txt") #없음
  #TweedData.RData: 3 columns and 12628 rows dataframe containing observed data with variables:
  #  location:: integer describing which stream segment data lie on, corresponds to row and column numbers in adjacency matrix (maybe stream number?)
  #date:: yyyy-mm-dd style date of observation point
  #nitrate:: value in mg/l of TON 
  #Long:: decimal longitude of sampling location
  #Lat:: decimal latitude of sampling location
  #data$normaldata_TN[,1]
  #data qqplot
  par(mfrow=c(1,2))
  for(ii in 1:ncol(data$normaldata_TOC)){
    qqnorm(as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii]))
    qqline(as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii]), col=2)
    qqnorm(log(as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii])))
    qqline(log(as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii])), col=2)
  }
  par(mfrow=c(1,1))
  #data 편집 (detrending)
  Tweed_global_mean <- mean(as.numeric(data$normaldata_TOC), na.rm=T) #3.991815
  Tweed_linear_scale <- matrix(nrow=ncol(data$normaldata_TOC),ncol=3)
  TweedData_old <- data.frame(date=c(), location=c(), nitrate=c(), Long=c(), Lat=c())
  #(OCT 17, 2020) CV를 위해 변경
  data_normaldata_old <- data$normaldata_TOC
  #data$normaldata_TOC <- data$normaldata_TOC[,-deleted_ind]
  #for(ii in setdiff(c(1: ncol(data$normaldata_TOC)), deleted_ind)){
  for(ii in 1: ncol(data$normaldata_TOC)){
    Tweed_linear_scale[ii,1] <- lm(as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii])~which(!is.na(data$normaldata_TOC[,ii])))$coeff[2]
    nitrate_scale <- scale(as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii]), center=T)
    Tweed_linear_scale[ii,2] <- attr(nitrate_scale, "scaled:center"); Tweed_linear_scale[ii,3] <- attr(nitrate_scale, "scaled:scale")
    TweedData_old_imsi <- data.frame(date=date(data$normaldata_TOC)[which(!is.na(data$normaldata_TOC[,ii]))], 
                                     location = rep(example_network@obspoints@SSNPoints[[1]]@point.data$rid[ii], length(which(!is.na(data$normaldata_TOC[,ii])))),
                                     nitrate = as.numeric(data$normaldata_TOC[which(!is.na(data$normaldata_TOC[,ii])), ii]),
                                     #nitrate = scale(as.numeric(data$normaldata_TN[which(!is.na(data$normaldata_TN[,ii])), ii]), center=F),
                                     #nitrate = scale(as.numeric(data$normaldata_TN[which(!is.na(data$normaldata_TN[,ii])), ii]), center=T)+Tweed_global_mean, #(***) 이것 사용
                                     #nitrate = scale(as.numeric(data$normaldata_TN[which(!is.na(data$normaldata_TN[,ii])), ii])-Tweed_linear_scale[ii,1]*which(!is.na(data$normaldata_TN[,ii])), center=T)+Tweed_global_mean,
                                     #Long=rep(eventplace2$long[ii], length(which(!is.na(data$normaldata_TN[,ii])))),
                                     Long=rep(example_network@obspoints@SSNPoints[[1]]@point.coords[ii,1], length(which(!is.na(data$normaldata_TOC[,ii])))),
                                     #Lat=rep(eventplace2$lat[ii], length(which(!is.na(data$normaldata_TN[,ii]))))
                                     Lat=rep(example_network@obspoints@SSNPoints[[1]]@point.coords[ii,2], length(which(!is.na(data$normaldata_TOC[,ii])))),
                                     index=ii
    )
    TweedData_old <- rbind(TweedData_old, TweedData_old_imsi)
  }
  
  data_plot <- apply(as.matrix(data$normaldata_TOC["2011-12/2017-11"]), 2, function(x) log(mean(x, na.rm=T)))
  data_Geum <- readRDS("~/Dropbox/R files/StreamFLow/data/TOCdata.RDS")
  data_Geum_duplicated <- readRDS("~/Dropbox/R files/StreamFLow/data/TOCdata(deleteduplicated).RDS")
  which(duplicated(names(data_Geum)))
  
  
  location_imsi <- example_network@obspoints@SSNPoints[[1]]@point.coords[-which(duplicated(names(data_Geum))),]
  #내가 새로 추가
  #TweedData generation (실제 자료의 duplication 반영)
  TweedData <- c(); n <- 1;   x <- 1/5 #임의의 시간대
  #TweedData generation (실제 자료의 duplication 반영)
  TweedData <- data.frame(date=as.Date("2012-12-31")+x*(365*3),
                          location=as.numeric(paste(names(data_Geum_duplicated)[1])),
                          nitrate=rep(exp(data_Geum_duplicated[1]),n),
                          #nitrate=result_signal[,example_network@obspoints@SSNPoints[[1]]@point.data$rid[1]] + rnorm(1, mean=0, sd.val),
                          Long=rep(location_imsi[1,1],n),
                          #Long=rep(location_imsi_true[1,1],n),
                          Lat=rep(location_imsi[1,2],n),
                          index=1
                          #Lat=rep(location_imsi_true[1,2],n)
  )
  for(i in 2:nrow(location_imsi)){
    dataframe.new <- data.frame(date=as.Date("2012-12-31")+x*(365*3),
                                location=as.numeric(paste(names(data_Geum_duplicated)[i])),
                                nitrate=rep(exp(data_Geum_duplicated[i]),n),
                                #nitrate=result_signal[,example_network@obspoints@SSNPoints[[1]]@point.data$rid[i]] + rnorm(1, mean=0, sd.val),
                                Long=rep(location_imsi[i,1],n),
                                #Long=rep(location_imsi_true[i,1],n),
                                Lat=rep(location_imsi[i,2],n),
                                index=i
                                #Lat=rep(location_imsi_true[i,2],n)
    )
    TweedData <- rbind(TweedData, dataframe.new)
  }
  rownames(TweedData) <- c(1:nrow(TweedData))
  
  #TweedData <- TweedData_old
  #2012년 2월 29일 제거
  #TweedData$date[which(TweedData$date=="2012-02-29")] <- "2012-02-28"
  #TweedData$date[which(TweedData$date=="2016-02-29")] <- "2016-02-28"
  
  # load high resolution matrix of points on the network
  # useful for plotting.
  #TweedPredPoints<-read.table("TweedPredPoints.txt")
  #역시 Lat, Long, SteamUnit, Weights로 구성됨
  #length(unique(TweedPredPoints$Latitude)) = 14114
  #length(unique(TweedPredPoints$Weights)) : 243
  #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col =TweedPredPoints$StreamUnit, pch=16, cex=TweedPredPoints$Weights)
  #points(TweedPredPoints$Longitude[1], TweedPredPoints$Latitude[2], pch=5, cex=2)
  #range(TweedPredPoints$Weights)
  
  #0.2부터 1.5사이로 weight 조정
  weight_vec_candidate <- 1.5*exp(-(example_network@network.line.coords$DistanceUpstream-min(example_network@network.line.coords$DistanceUpstream))/0.4013373)
  weight_vec_candidate2 <- compute_shreve_and_dist(adjacency_old, example_network, scalevec=c(0.2, 1.5))
  
  #weight_vec_candidate2 <- scales::rescale(log(weight_vec_candidate2$distweight), to=c(0.2, 1.5))
  weight_vec_candidate <- weight_vec_candidate2$distweight
  
  TweedPredPoints_old <- data.frame(Latitude=c(), Longitude=c(), StreamUnit=c(), Weights=c())
  for(jj in 1:length(example_network@lines)){
    #if(nrow(example_network@lines[[jj]]@Lines[[1]]@coords)>2){
    #  seg_cand <- example_network@lines[[jj]]@Lines[[1]]@coords[-nrow(example_network@lines[[jj]]@Lines[[1]]@coords),]
    #}else{
    seg_cand <- example_network@lines[[jj]]@Lines[[1]]@coords
    #}
    TweedPredPoints_old_imsi <- data.frame(Latitude=seg_cand[,2], Longitude=seg_cand[,1], 
                                           StreamUnit=rep(example_network@data$rid[jj], nrow(seg_cand)), 
                                           Weights=rep(weight_vec_candidate[jj], nrow(seg_cand)))
    TweedPredPoints_old <- rbind(TweedPredPoints_old, TweedPredPoints_old_imsi)
  }
  TweedPredPoints <- TweedPredPoints_old
  
  #(Oct 16, 2020): CV를 위해 추가

  #example_network@obspoints@SSNPoints[[1]]@network.point.coords
  #example_network@obspoints@SSNPoints[[1]]@point.coords
  #example_network@obspoints@SSNPoints[[1]]@point.data
  #example_network@obspoints@SSNPoints[[1]]@points.bbox
  #example_network@obspoints@SSNPoints[[1]]@proj4string
  example_network@predpoints@SSNPoints[[1]]@network.point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@network.point.coords, example_network@obspoints@SSNPoints[[1]]@network.point.coords[which(duplicated(names(data_Geum))),])
  example_network@obspoints@SSNPoints[[1]]@network.point.coords <- example_network@obspoints@SSNPoints[[1]]@network.point.coords[-which(duplicated(names(data_Geum))),]
  example_network@predpoints@SSNPoints[[1]]@point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@point.coords, example_network@obspoints@SSNPoints[[1]]@point.coords[which(duplicated(names(data_Geum))),])
  example_network@obspoints@SSNPoints[[1]]@point.coords <- example_network@obspoints@SSNPoints[[1]]@point.coords[-which(duplicated(names(data_Geum))),]
  example_network@predpoints@SSNPoints[[1]]@point.data <- rbind(example_network@predpoints@SSNPoints[[1]]@point.data, example_network@obspoints@SSNPoints[[1]]@point.data[which(duplicated(names(data_Geum))),])
  example_network@obspoints@SSNPoints[[1]]@point.data <- example_network@obspoints@SSNPoints[[1]]@point.data[-which(duplicated(names(data_Geum))),]
  
  example_network_old <- example_network
  
  deleted_ind <- deleted_ind_cand[jkl]
  
  TweedData <- TweedData[-deleted_ind,]
  
  example_network@predpoints@SSNPoints[[1]]@network.point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@network.point.coords, example_network@obspoints@SSNPoints[[1]]@network.point.coords[deleted_ind,])
  example_network@obspoints@SSNPoints[[1]]@network.point.coords <- example_network@obspoints@SSNPoints[[1]]@network.point.coords[-deleted_ind,]
  example_network@predpoints@SSNPoints[[1]]@point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@point.coords, example_network@obspoints@SSNPoints[[1]]@point.coords[deleted_ind,])
  example_network@obspoints@SSNPoints[[1]]@point.coords <- example_network@obspoints@SSNPoints[[1]]@point.coords[-deleted_ind,]
  example_network@predpoints@SSNPoints[[1]]@point.data <- rbind(example_network@predpoints@SSNPoints[[1]]@point.data, example_network@obspoints@SSNPoints[[1]]@point.data[deleted_ind,])
  example_network@obspoints@SSNPoints[[1]]@point.data <- example_network@obspoints@SSNPoints[[1]]@point.data[-deleted_ind,]
  
  
  example_network@predpoints@SSNPoints[[1]]@points.bbox <- rbind(range(example_network@predpoints@SSNPoints[[1]]@point.coords[,1]),
                                                                 range(example_network@predpoints@SSNPoints[[1]]@point.coords[,2]))
  rownames(example_network@predpoints@SSNPoints[[1]]@points.bbox) <- c("coords.x1", "coords.x2")
  colnames(example_network@predpoints@SSNPoints[[1]]@points.bbox) <- c("min", "max")
  example_network@obspoints@SSNPoints[[1]]@points.bbox <- rbind(range(example_network@obspoints@SSNPoints[[1]]@point.coords[,1]),
                                                                range(example_network@obspoints@SSNPoints[[1]]@point.coords[,2]))
  rownames(example_network@obspoints@SSNPoints[[1]]@points.bbox) <- c("coords.x1", "coords.x2")
  colnames(example_network@obspoints@SSNPoints[[1]]@points.bbox) <- c("min", "max")
  
  
  
  #역시 Lat, Long, SteamUnit, Weights로 구성됨
  #length(unique(TweedPredPoints$Latitude)) = 14114
  #length(unique(TweedPredPoints$Weights)) : 243
  par(mar=c(1.1,1.1,1.1,1.1))
  par(mfrow=c(1,1))
  plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col =TweedPredPoints$StreamUnit, pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", xaxt="n", yaxt="n")
  text(x=mf04@obspoints@SSNPoints[[1]]@point.coords[,1]+c(rep(0.005,10),0.02,rep(0.005,17)), y=mf04@obspoints@SSNPoints[[1]]@point.coords[,2]+c(rep(0.005,17),0.02,rep(0.005,10)), labels=mf04@obspoints@SSNPoints[[1]]@point.data$rid, col="brown", cex=1.3)
  
  points(TweedPredPoints$Longitude[1], TweedPredPoints$Latitude[2], pch=5, cex=2)
  
  #realweights <- as.matrix(weight_vec_candidate[example_network@obspoints@SSNPoints[[1]]@point.data$rid], ncol=1)
  realweights <- as.matrix(weight_vec_candidate, ncol=1)
  
  #######################################
  #source_Flexible.R code 돌리기
  ########################################
  #(추가: 만약 연도별로 할 것이면...)
  #TweedData <- TweedData[which(year(TweedData$date)=="2012"),]
  #source('~/Dropbox/(논문들)/(논문)RIver Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/Me/source_Flexible.R', chdir = TRUE)
  source('~/Dropbox/R files/StreamFLow/sources/source_Flexible.R', chdir = TRUE)
  penalties_default <- c(50, 50, c(25, 5), 50, c(25, 5), c(25, 25))
  division_factor <- c(30, 25, 20, 15, 10, 5, 1, 0.5, 0.25, 0.125)
  AICc_vec <- rep(0, length(division_factor))
  for(jjj in 1:length(division_factor)){
    result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/division_factor[jjj], plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=TRUE, model.type = c("c", "m"))
    AICc_vec[jjj] <- result_a$AICc
  }
  AICc_vec_min_index <- which.min(AICc_vec)
  
  if(AICc_vec_min_index==1){
    opt_result <- optim(par=penalties_default/division_factor[AICc_vec_min_index], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index)], upper=penalties_default/division_factor[(AICc_vec_min_index+1)], log.y=TRUE, model.type = c("c", "m") )
  }else if(AICc_vec_min_index==length(AICc_vec)){
    opt_result <- optim(par=penalties_default/division_factor[AICc_vec_min_index], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index)], upper=penalties_default/division_factor[(AICc_vec_min_index)], log.y=TRUE, model.type = c("c", "m") )
  }else{
    opt_result <- optim(par=penalties_default/division_factor[AICc_vec_min_index], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index-1)], upper=penalties_default/division_factor[(AICc_vec_min_index+1)], log.y=TRUE, model.type = c("c", "m") )
  }
  opt_val <- opt_result$par
  
  location_unique <- unique(TweedData$location)
  result_7 <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=opt_val, plot.fig=TRUE, station=location_unique, use.optim=FALSE, log.y=TRUE,model.type = c("c", "m") )
  result_7$AICc
  
  result_7_add_residual <- cbind(result_7$TweedData, fit=result_7$fit, resids=result_7$resids) #resids<- fit - log(response) #fit에서 log 뺀 값 (뒤에서 쓰임)
  result_7_split_by_loc <- split(result_7_add_residual, result_7_add_residual$Long)
  newnames <- result_7_add_residual$location[which(!duplicated(result_7_add_residual$Long))]
  ########################################
  #Stream network lifting scheme (출력할 가치 있음, original)
  ########################################
  example_network@obspoints@SSNPoints[[1]]@point.data$rid
  example_network@obspoints@SSNPoints[[1]]@point.data$segid
  
  #split(TweedData, f = TweedData$location)
  result_7_split_by_ind <- split(result_7_add_residual, result_7_add_residual$index)
  data <- sapply(result_7_split_by_ind, function(x) log(mean(x$nitrate)))
  names(data) <- result_7_add_residual$location[which(!duplicated(result_7_add_residual$Long))][order(result_7_add_residual$Long[which(!duplicated(result_7_add_residual$Long))])]
  data <- data[match(example_network@obspoints@SSNPoints[[1]]@point.data$rid, names(data))]
  
  #수정 (Oct 18, 2020)
  data <- readRDS("~/Dropbox/R files/StreamFLow/data/TOCdata(deleteduplicated).RDS")
  data_Tweed <- sapply(result_7_split_by_ind, function(x) log(mean(x$nitrate)))
  if(sum(duplicated(data))!=0){
    stop("data generation error")
  }else if(sum(data_Tweed!=data[-deleted_ind])!=0){
    stop("data does not match")
  }
  data <- data[-deleted_ind]
  
  
  result_forward <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
  
  #data <- readRDS("~/Dropbox/R files/StreamFLow/data/TOCdata(deleteduplicated).RDS")
  #time_init <- Sys.time()
  result_denoise <- denoise_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE)
  #time_end <- Sys.time()
  #time_end - time_init
  
  colnames(result_denoise) <- colnames(data)
  
  cores=detectCores()
  cl <- makeCluster(cores[1]-3) #not to overload your computer
  ##subset groupping
  #이걸 8로 하면 .rds를 (2)로, 6으로 하면 .rds를 (3)으로 저장
  index_sub_data <- data.frame(stations=c(1:length(data)), groups=as.factor(substr(example_network@obspoints@SSNPoints[[1]]@point.data$RCH_ID, start=1, stop=6)))
  registerDoParallel(cl)
  #result_denoise_nlt <- denoise_Stream_S_perm(data, example_network, endpt=NULL, per=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = keep, rule = "median", sd.scale=1, returnall = FALSE)
  
  #time_init <- Sys.time()
  result_nlt <- nlt_Stream_S(data, example_network, J=5, endpt=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = TRUE, index_sub_data=index_sub_data, max.subnum=5, ga=TRUE, oremovelist=result_forward$removelist)
  #time_end <- Sys.time()
  #time_end - time_init
  #J=1, max.subnum=0으로 하면 S-Lifting(M) 결과와 정확하게 같음
  stopCluster(cl)
  
  #data_old <- data
  #data <- readRDS("~/Dropbox/R files/StreamFLow/data/TOCdata(deleteduplicated).RDS")
  #data_Tweed <- sapply(result_7_split_by_ind, function(x) log(mean(x$nitrate)))
  #if(sum(duplicated(data))!=0){
  #  stop("data generation error")
  #}else if(sum(data_Tweed!=data[-deleted_ind])!=0){
  #  stop("data does not match")
  #}
  #data <- data[-deleted_ind]
  
  #initS_obj: needed to prediction
  initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
  data_predicted <- initS_obj$weight_matrix%*%data
  data_predicted0 <- data_predicted
  
  initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
  data_predicted <- initS_obj$weight_matrix%*%t(result_denoise) #result_denoise: proposed method
  
  initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
  data_predicted2 <- initS_obj$weight_matrix%*%result_nlt$aveghat
  
  denoising_ODonell <- (result_7$fit[which(!duplicated(result_7$TweedData$Lat))])
  initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
  data_predicted1 <- initS_obj$weight_matrix%*%denoising_ODonell
  
  zlims <- range(c(data, data_predicted, data_predicted0, data_predicted1, data_predicted2))+c(-0.05, 0.05)
  
  ########################################
  #Prediction error check
  ########################################
  pred_info <- example_network@predpoints@SSNPoints[[1]]@point.data[nrow(example_network@predpoints@SSNPoints[[1]]@point.data),]
  
  TOCdata <- readRDS("~/Dropbox/R files/StreamFLow/data/TOCdata(deleteduplicated).RDS")
  
  data_predicted0[pred_info$rid] #orig_data
  data_predicted[pred_info$rid] #denoise_S
  data_predicted1[pred_info$rid] #denoising_ODonell
  data_predicted2[pred_info$rid] #result_nlt
  
  (TOCdata[deleted_ind] - data_predicted[pred_info$rid])^2
  (TOCdata[deleted_ind] - data_predicted1[pred_info$rid])^2
  (TOCdata[deleted_ind] - data_predicted2[pred_info$rid])^2
  
  result_data_frame <- data.frame(rid=pred_info$rid,
                                  TOCdata=TOCdata[deleted_ind],
                                  Raw=data_predicted0[pred_info$rid],
                                  ODonnell=data_predicted1[pred_info$rid],
                                  SLifting=data_predicted[pred_info$rid],
                                  NLifting=data_predicted2[pred_info$rid],
                                  RMSERaw=sqrt((TOCdata[deleted_ind] - data_predicted0[pred_info$rid])^2),
                                  RMSEODonnell=sqrt((TOCdata[deleted_ind] - data_predicted1[pred_info$rid])^2),
                                  RMSESLifting=sqrt((TOCdata[deleted_ind] - data_predicted[pred_info$rid])^2),
                                  RMSENLifting=sqrt((TOCdata[deleted_ind] - data_predicted2[pred_info$rid])^2)
  )
  
  saveRDS(result_data_frame, paste("~/Dropbox/R files/StreamFLow/CVResult/CV", deleted_ind , "(dd1).RDS", sep=""))
  
}

