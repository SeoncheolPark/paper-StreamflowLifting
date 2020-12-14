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
#binaryIDs_obj <- get_binaryIDs_stream(mouth_node=614 , shape_miho, RID="corrected")
#수동으로 adjacency를 얻도록 만들어진 함수
adjacency <- get_adjacency_stream(binaryIDs_obj)

#compute_shreve: shreve 거리를 재도록 직접 제작한 함수
shreve_order <- compute_shreve(adjacency)

#generate SSN obj
mf04 <- importSSN_stream(shreve_obj = shreve_order, location="Full", multipleplaces=TRUE) #데이터의 크기: 123개
#mf04 <- importSSN_stream(shreve_obj = shreve_order, location="Full", multipleplaces=TRUE, RID="corrected") #데이터의 크기: 123개

which(mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID!=mf04@obspoints@SSNPoints[[1]]@point.data$rid)
mf04@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID
mf04@obspoints@SSNPoints[[1]]@point.data$rid
#check basics
dim(mf04@obspoints@SSNPoints[[1]]@point.data) #127 38
length(unique(mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID)) #114
which(duplicated(mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID)) #13
mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID[which(duplicated(mf04@obspoints@SSNPoints[[1]]@point.data$RCH_ID))] #3개짜리 하나 있음
show_weights(mf04, adjacency)
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
realweights<-read.table("realweights.txt") #없음
#realweights.RData : vector of flow weights describing proportional flow contribution to nearest downstream neighbour on the River Tweed.

# load adjacency matrix associated with each 
# of the 298 stream segments of the river tweed
#(관측 장소가 아닌 298개 stream segment들의 adjacency를 나타낸 것이다)
adjacency_old <- adjacency 
adjacency<-read.table("adjacency.txt")
#adjacency.RData:  298 by 298 sparse matrix defining neighbourhood structure of river network - if adjacency[i, j] == 1 then segment i flows into segment j.
##########나의 해결책
adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))

# load river tween data: 3 col dataframe:
# dates, nitrate concs (mg/l) and 
# location index (corresponding row number of adjacency matrix)
TweedData<-read.table("TweedData.txt") #없음
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
                                   Lat=rep(example_network@obspoints@SSNPoints[[1]]@point.coords[ii,2], length(which(!is.na(data$normaldata_TOC[,ii]))))
  )
  TweedData_old <- rbind(TweedData_old, TweedData_old_imsi)
}
TweedData <- TweedData_old
#2012년 2월 29일 제거
TweedData$date[which(TweedData$date=="2012-02-29")] <- "2012-02-28"
TweedData$date[which(TweedData$date=="2016-02-29")] <- "2016-02-28"

# load high resolution matrix of points on the network
# useful for plotting.
TweedPredPoints<-read.table("TweedPredPoints.txt")
#역시 Lat, Long, SteamUnit, Weights로 구성됨
#length(unique(TweedPredPoints$Latitude)) = 14114
#length(unique(TweedPredPoints$Weights)) : 243
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col =TweedPredPoints$StreamUnit, pch=16, cex=TweedPredPoints$Weights)
points(TweedPredPoints$Longitude[1], TweedPredPoints$Latitude[2], pch=5, cex=2)
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

#Plot Figure 1
par(family = 'sans') 
par(mar=c(1.1,1.1,1.1,4.1))
par(mfrow=c(1,1))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="", xaxt="n", yaxt='n')
data_plot <- apply(as.matrix(data$normaldata_TOC["2011-12/2017-11"]), 2, function(x) log(mean(x, na.rm=T)))
quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], data_plot, main="", cex=2, add=T, zlim=range(data_plot), xlab="", ylab="", add.legend = T, legend.lab="log(mean(TOC))")
#save plot with 500*450
########################################
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

data <- sapply(result_7_split_by_loc, function(x) log(mean(x$nitrate)))
names(data) <- result_7_add_residual$location[which(!duplicated(result_7_add_residual$Long))][order(result_7_add_residual$Long[which(!duplicated(result_7_add_residual$Long))])]
data <- data[match(example_network@obspoints@SSNPoints[[1]]@point.data$rid, names(data))]

#initS_obj: needed to prediction
initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
data_predicted <- initS_obj$weight_matrix%*%data
data_predicted0 <- data_predicted

result_forward <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)

result_denoise <- denoise_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE)
colnames(result_denoise) <- colnames(data)
initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
data_predicted <- initS_obj$weight_matrix%*%t(result_denoise) #result_denoise: proposed method

cores=detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
##subset groupping
index_sub_data <- data.frame(stations=c(1:length(data)), groups=as.factor(substr(example_network@obspoints@SSNPoints[[1]]@point.data$RCH_ID, start=1, stop=8)))
registerDoParallel(cl)
#result_denoise_nlt <- denoise_Stream_S_perm(data, example_network, endpt=NULL, per=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = keep, rule = "median", sd.scale=1, returnall = FALSE)
result_nlt <- nlt_Stream_S(data, example_network, J=10, endpt=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = TRUE, index_sub_data=index_sub_data, max.subnum=5, ga=TRUE, oremovelist=result_forward$removelist)
#J=1, max.subnum=0으로 하면 S-Lifting(M) 결과와 정확하게 같음
stopCluster(cl)
initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
data_predicted2 <- initS_obj$weight_matrix%*%result_nlt$aveghat

denoising_ODonell <- (result_7$fit[which(!duplicated(result_7$TweedData$Lat))])
initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
data_predicted1 <- initS_obj$weight_matrix%*%denoising_ODonell

zlims <- range(c(data, data_predicted, data_predicted0, data_predicted1, data_predicted2))+c(-0.05, 0.05)

#pdf("result_TOC3.pdf", 7, 7)
#png("result_TOC3.png", 700, 700)
par(family = 'sans') 
par(mar=c(3.1,2.1,3.1,1.1))
par(mfrow=c(2,2))
scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted0[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(a) Raw, log(TOC)", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825))
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], data, main="(a) Raw, TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3), bigplot=c(0,0,1,1), xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted1[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(b) ODonnell, log(TOC)", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], denoising_ODonell, main="(b) ODonnell, TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(c) S-Lifting (M), log(TOC)", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825))
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], result_denoise, main="(c) S-Lifting (M), TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(d) S-Lifting (N), log(TOC)", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], result_nlt$aveghat, main="(d) S-Lifting (N), TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))

#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted0[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(a) Raw, TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=data, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))
#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted1[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(b) ODonnell, TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE, quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=denoising_ODonell, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))
#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(c) S-Lifting (M), TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=result_denoise, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))
#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(d) S-Lifting (N), TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE, quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=result_nlt$aveghat, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))

dev.off()

#Table 1과 연계한 plot
#pdf("result_TOC3_Table1.pdf", 7, 7)
#png("result_TOC3_Table1.png", 700, 700)
# choose colors to interpolate
data_range <- c(data_predicted, data_predicted0, data_predicted1, data_predicted2)
#levels <- c(data_range, 2, 3, 4, 5, 6, 8)
#col <- colorRampPalette(c("red","yellow","dark green"))(nlevels) 
#following ODonnell's Approach
#col<-rainbowPalette(7)
#col <- colorRampPalette(c("red","yellow","dark green"))(7) 
col <- colorRampPalette(c("blue", "cyan", "green", "gray", "yellow", "orange", "red"))(7) 
#colz <- col[cut(z, breaks=c(min(data_range),log(c(2,3,4,5,6,8)), max(data_range)))]  

par(family = 'sans') 
par(mar=c(3.1,2.1,3.1,1.1))
par(mfrow=c(2,2))
colz <- col[cut(data_predicted0, breaks=c(min(data_range),log(c(2,3,4,5,6,8)), max(data_range)))]  
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=colz[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(a) Raw, log(TOC)", cex.main=1.5)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], data, main="(a) Raw, TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3), bigplot=c(0,0,1,1), xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude))
legend("bottomleft", c("Ia", "Ib", "II", "III", "IV", "V", "VI"), col=c("blue", "cyan", "green", "gray", "yellow", "orange", "red"), pch=19, ncol=4)

colz <- col[cut(data_predicted1, breaks=c(min(data_range),log(c(2,3,4,5,6,8)), max(data_range)))]  
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=colz[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(b) ODonnell, log(TOC)", cex.main=1.5, plot.legend=TRUE, axes=TRUE, y.axes=FALSE)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], denoising_ODonell, main="(b) ODonnell, TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))
legend("bottomleft", c("Ia", "Ib", "II", "III", "IV", "V", "VI"), col=c("blue", "cyan", "green", "gray", "yellow", "orange", "red"), pch=19, ncol=4)

colz <- col[cut(data_predicted, breaks=c(min(data_range),log(c(2,3,4,5,6,8)), max(data_range)))]  
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=colz[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(c) S-Lifting (M), log(TOC)", cex.main=1.5)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], result_denoise, main="(c) S-Lifting (M), TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))
legend("bottomleft", c("Ia", "Ib", "II", "III", "IV", "V", "VI"), col=c("blue", "cyan", "green", "gray", "yellow", "orange", "red"), pch=19, ncol=4)

colz <- col[cut(data_predicted2, breaks=c(min(data_range),log(c(2,3,4,5,6,8)), max(data_range)))]  
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=colz[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", main="(d) S-Lifting (N), log(TOC)", cex.main=1.5, plot.legend=TRUE, axes=TRUE, y.axes=FALSE)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], result_nlt$aveghat, main="(d) S-Lifting (N), TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))
legend("bottomleft", c("Ia", "Ib", "II", "III", "IV", "V", "VI"), col=c("blue", "cyan", "green", "gray", "yellow", "orange", "red"), pch=19, ncol=4)

#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted0[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(a) Raw, TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=data, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))
#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted1[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(b) ODonnell, TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE, quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=denoising_ODonell, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))
#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(c) S-Lifting (M), TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=result_denoise, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))
#scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(d) S-Lifting (N), TOC", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE, quilt.type="quilt", quilt.x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1], quilt.y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], quilt.z=result_nlt$aveghat, quilt.zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)))

dev.off()

##compute
result_new_mat <- matrix(0, nrow=3, ncol=5)
rownames(result_new_mat) <- c("S-Lifting(M)", "ODonnell", "S-Lifting(N)")
colnames(result_new_mat) <- c("Corr", "RMSE(whole)","WRMSE(whole)", "RMSE(actual)", "WRMSE(actual)")

result_new_mat[1,1] <- cor(data_predicted0, data_predicted) #S-Lifting (M)
result_new_mat[2,1] <- cor(data_predicted0, data_predicted1) #ODonnell
result_new_mat[3,1] <- cor(data_predicted0, data_predicted2) #S-Lifting (N)

#whole
result_new_mat[1,2] <- sqrt(sum((data_predicted0-data_predicted)^2)/length(data_predicted0)) #S-Lifting (M)
result_new_mat[2,2] <- sqrt(sum((data_predicted0-data_predicted1)^2)/length(data_predicted0)) #ODonnell
result_new_mat[3,2] <- sqrt(sum((data_predicted0-data_predicted2)^2)/length(data_predicted0)) #S-Lifting (N)

#whole, weighted ver.
result_new_mat[1,3] <- sqrt(sum(realweights*(data_predicted0-data_predicted)^2)/length(data_predicted0)) #S-Lifting (M)
result_new_mat[2,3] <- sqrt(sum(realweights*(data_predicted0-data_predicted1)^2)/length(data_predicted0)) #ODonnell
result_new_mat[3,3] <- sqrt(sum(realweights*(data_predicted0-data_predicted2)^2)/length(data_predicted0)) #S-Lifting (N)

index_realdata <- sort(unique(as.row(as.numeric(names(data)))))

#actual
result_new_mat[1,4] <- sqrt(sum((data_predicted0[index_realdata,1]-data_predicted[index_realdata,1])^2)/length(index_realdata)) #S-Lifting (M)
result_new_mat[2,4] <- sqrt(sum((data_predicted0[index_realdata,1]-data_predicted1[index_realdata,1])^2)/length(index_realdata)) #ODonnell
result_new_mat[3,4] <- sqrt(sum((data_predicted0[index_realdata,1]-data_predicted2[index_realdata,1])^2)/length(index_realdata)) #S-Lifting (N)

#actual, weighted ver.
result_new_mat[1,5] <- sqrt(sum(realweights[index_realdata,1]*(data_predicted0[index_realdata,1]-data_predicted[index_realdata,1])^2)/length(index_realdata)) #S-Lifting (M)
result_new_mat[2,5] <- sqrt(sum(realweights[index_realdata,1]*(data_predicted0[index_realdata,1]-data_predicted1[index_realdata,1])^2)/length(index_realdata)) #ODonnell
result_new_mat[3,5] <- sqrt(sum(realweights[index_realdata,1]*(data_predicted0[index_realdata,1]-data_predicted2[index_realdata,1])^2)/length(index_realdata)) #S-Lifting (N)

#saveRDS(result_new_mat, "GeumRiver.RDS")

########################################
#교수님의 코멘트 적용
########################################
#1. 가장 fine한 level의 representation
#S-lifting
#이 중에서 한 가지만 선택
#가장 왼쪽 (가장 fine한 level)
#Final fitting 하지 말고 
#Forward lifting scheme의 
#반 (forward lifting scheme 그림 표현)
#반 (forward lifting scheme 그림 표현)
#네 개의 인덱스 중에 하나 가장 fine, 가장 coaser, 중간에 두 개
#raw plot
#result_forward_finerplot_1 <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)

result_forward_finerplot_2 <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 32, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin=result_forward_finerplot_2$pointsin)
data_predicted_finerplot_2 <- initS_obj$weight_matrix%*%result_forward_finerplot_2$coeff[result_forward_finerplot_2$pointsin]

#check
lr <- result_forward_finerplot_2$lengthsremove
rem <- result_forward_finerplot_2$removelist
newcoeff <- result_forward_finerplot_2$coeff
keep=32
int = TRUE; clo = FALSE
fhat <- invtnp_stream_S(X=as.vector(as.numeric(names(data))), newcoeff, result_forward_finerplot_2$lengths, lr, result_forward_finerplot_2$pointsin,  rem, result_forward_finerplot_2$neighbrs, result_forward_finerplot_2$schemehist, result_forward_finerplot_2$interhist, length(result_forward_finerplot_2$x) - keep, int, neighbours=1, clo, LocalPred=streamPred_S,  data, example_network, adjacency=adjacency_old, realweights)
round(data-fhat$coeff,5) #should be zero vectors

initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin=c(setdiff(c(1:length(data)),result_forward_finerplot_2$pointsin),127))
#data_predicted_finerplot_2_detail <- initS_obj$weight_matrix%*%result_forward_finerplot_2$coeff[-setdiff(result_forward_finerplot_2$pointsin,127)]
data_predicted_finerplot_2_detail <- initS_obj$weight_matrix%*%c(result_forward_finerplot_2$coeff[-result_forward_finerplot_2$pointsin], 0)
initS_obj_2 = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin=which(!is.na(data_predicted_finerplot_2_detail)))
data_predicted_finerplot_2_detail_new_2 <- initS_obj_2$weight_matrix%*%data_predicted_finerplot_2_detail[complete.cases(data_predicted_finerplot_2_detail)]


#initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(-data), example_network, adjacency=adjacency_old, realweights, pointsin=setdiff(c(1:length(data)),result_forward_finerplot_2$pointsin))
#data_predicted_finerplot_2_detail_new <- initS_obj$weight_matrix%*%(result_forward_finerplot_2$coeff[-result_forward_finerplot_2$pointsin])
#initS_obj_2 = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin=which(!is.na(data_predicted_finerplot_2_detail_new)))
#data_predicted_finerplot_2_detail_new_2 <- initS_obj_2$weight_matrix%*%data_predicted_finerplot_2_detail_new[complete.cases(data_predicted_finerplot_2_detail_new)]


#next level
datanew <- as.numeric(result_forward_finerplot_2$coeff[result_forward_finerplot_2$pointsin])
names(datanew) <- names(data)[result_forward_finerplot_2$pointsin]
example_network_new <- example_network
example_network_new@obspoints@SSNPoints[[1]]@network.point.coords <- example_network_new@obspoints@SSNPoints[[1]]@network.point.coords[result_forward_finerplot_2$pointsin,]
example_network_new@obspoints@SSNPoints[[1]]@point.coords <- example_network_new@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_2$pointsin,]
example_network_new@obspoints@SSNPoints[[1]]@point.data <- example_network_new@obspoints@SSNPoints[[1]]@point.data[result_forward_finerplot_2$pointsin,]
example_network_new@obspoints@SSNPoints[[1]]@points.bbox[1,] <- range(example_network_new@obspoints@SSNPoints[[1]]@point.data$경도.Degree.)
example_network_new@obspoints@SSNPoints[[1]]@points.bbox[2,] <- range(example_network_new@obspoints@SSNPoints[[1]]@point.data$위도.Degree.)

result_forward_finerplot_3 <- fwtnp_stream_S(datanew, example_network_new, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 8, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
#result_forward_finerplot_3 <- fwtnp_stream_S(data, example_network, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 8, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
initS_obj = initS_stream(X=as.row(as.numeric(names(datanew))), data=as.vector(datanew), example_network_new, adjacency=adjacency_old, realweights, pointsin=result_forward_finerplot_3$pointsin)
data_predicted_finerplot_3 <- initS_obj$weight_matrix%*%result_forward_finerplot_3$coeff[result_forward_finerplot_3$pointsin]

initS_obj = initS_stream(X=as.row(as.numeric(names(datanew))), data=as.vector(datanew), example_network_new, adjacency=adjacency_old, realweights, pointsin=c(setdiff(c(1:length(datanew)),result_forward_finerplot_3$pointsin),length(datanew)))
#data_predicted_finerplot_3_detail <- initS_obj$weight_matrix%*%result_forward_finerplot_3$coeff[-setdiff(result_forward_finerplot_3$pointsin,32)]
data_predicted_finerplot_3_detail <- initS_obj$weight_matrix%*%c(result_forward_finerplot_3$coeff[-result_forward_finerplot_3$pointsin], 0)
initS_obj_2 = initS_stream(X=as.row(as.numeric(names(datanew))), data=as.vector(datanew), example_network_new, adjacency=adjacency_old, realweights, pointsin=which(!is.na(data_predicted_finerplot_3_detail)))
data_predicted_finerplot_3_detail_new_2 <- initS_obj_2$weight_matrix%*%data_predicted_finerplot_3_detail[complete.cases(data_predicted_finerplot_3_detail)]


#initS_obj = initS_stream(X=as.row(as.numeric(names(datanew))), data=as.vector(datanew), example_network_new, adjacency=adjacency_old, realweights, pointsin=setdiff(c(1:length(datanew)),result_forward_finerplot_3$pointsin))
#data_predicted_finerplot_3_detail_new <- initS_obj$weight_matrix%*%(-result_forward_finerplot_3$coeff[-result_forward_finerplot_3$pointsin])
#initS_obj_2 = initS_stream(X=as.row(as.numeric(names(datanew))), data=as.vector(datanew), example_network_new, adjacency=adjacency_old, realweights, pointsin=which(!is.na(data_predicted_finerplot_3_detail_new)))
#data_predicted_finerplot_3_detail_new_2 <- initS_obj_2$weight_matrix%*%data_predicted_finerplot_3_detail_new[complete.cases(data_predicted_finerplot_3_detail_new)]


#next level
datanewnew <- as.numeric(result_forward_finerplot_3$coeff[result_forward_finerplot_3$pointsin])
names(datanewnew) <- names(datanew)[result_forward_finerplot_3$pointsin]
example_network_newnew <- example_network_new
example_network_newnew@obspoints@SSNPoints[[1]]@network.point.coords <- example_network_newnew@obspoints@SSNPoints[[1]]@network.point.coords[result_forward_finerplot_3$pointsin,]
example_network_newnew@obspoints@SSNPoints[[1]]@point.coords <- example_network_newnew@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_3$pointsin,]
example_network_newnew@obspoints@SSNPoints[[1]]@point.data <- example_network_newnew@obspoints@SSNPoints[[1]]@point.data[result_forward_finerplot_3$pointsin,]
example_network_newnew@obspoints@SSNPoints[[1]]@points.bbox[1,] <- range(example_network_newnew@obspoints@SSNPoints[[1]]@point.data$경도.Degree.)
example_network_newnew@obspoints@SSNPoints[[1]]@points.bbox[2,] <- range(example_network_newnew@obspoints@SSNPoints[[1]]@point.data$위도.Degree.)

result_forward_finerplot_4 <- fwtnp_stream_S(datanewnew, example_network_newnew, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
initS_obj = initS_stream(X=as.row(as.numeric(names(datanewnew))), data=as.vector(datanewnew), example_network_newnew, adjacency=adjacency_old, realweights, pointsin=result_forward_finerplot_4$pointsin)
data_predicted_finerplot_4 <- initS_obj$weight_matrix%*%result_forward_finerplot_4$coeff[result_forward_finerplot_4$pointsin]


#initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network, adjacency=adjacency_old, realweights, pointsin=c(setdiff(c(1:length(data)),result_forward_finerplot_4$pointsin),127))
initS_obj = initS_stream(X=as.row(as.numeric(names(datanewnew))), data=as.vector(datanewnew), example_network_newnew, adjacency=adjacency_old, realweights, pointsin=c(setdiff(c(1:length(datanewnew)),result_forward_finerplot_4$pointsin),length(datanewnew)))
data_predicted_finerplot_4_detail <- initS_obj$weight_matrix%*%result_forward_finerplot_4$coeff[-setdiff(result_forward_finerplot_4$pointsin,8)]
data_predicted_finerplot_4_detail <- initS_obj$weight_matrix%*%c(result_forward_finerplot_4$coeff[-result_forward_finerplot_4$pointsin], 0)
initS_obj_2 = initS_stream(X=as.row(as.numeric(names(datanewnew))), data=as.vector(datanewnew), example_network_newnew, adjacency=adjacency_old, realweights, pointsin=which(!is.na(data_predicted_finerplot_4_detail)))
data_predicted_finerplot_4_detail_new_2 <- initS_obj_2$weight_matrix%*%data_predicted_finerplot_4_detail[complete.cases(data_predicted_finerplot_4_detail)]


#initS_obj = initS_stream(X=as.row(as.numeric(names(datanewnew))), data=as.vector(datanewnew), example_network_newnew, adjacency=adjacency_old, realweights, pointsin=setdiff(c(1:length(datanewnew)),result_forward_finerplot_4$pointsin))
#data_predicted_finerplot_4_detail_new <- initS_obj$weight_matrix%*%(-result_forward_finerplot_4$coeff[-result_forward_finerplot_4$pointsin])
#initS_obj_2 = initS_stream(X=as.row(as.numeric(names(datanewnew))), data=as.vector(datanewnew), example_network_newnew, adjacency=adjacency_old, realweights, pointsin=which(!is.na(data_predicted_finerplot_4_detail_new)))
#data_predicted_finerplot_4_detail_new_2 <- initS_obj_2$weight_matrix%*%data_predicted_finerplot_4_detail_new[complete.cases(data_predicted_finerplot_4_detail_new)]


########################################
#교수님의 코멘트 적용
########################################
#2. 4*2 행렬 (이 그림을 그려야 O'Donnell의 논문과 차별이 가능하다)
#(1,1): 가장 finest한 피팅 결과 
#(2,1): 2번째로 finish (2,2): 2번째 결과의 detail
#(3,1): 3번째로 finish (3,2): 3번째 결과의 detail
#(4,1): 2번째로 finish (4,2): 2번째 결과의 detail

#pdf("result_TOC_detail3(rev).pdf", 6, 9)
#png("result_TOC_detail3(rev).png", 600, 900)

par(family = 'sans') 
par(mar=c(3.1,2.1,3.1,1.1))
par(mfrow=c(2,3))

range.val <- max(abs(c(data, data_predicted0, data_predicted_finerplot_2, data_predicted_finerplot_3, data_predicted_finerplot_4, data_predicted_finerplot_2_detail_new_2, data_predicted_finerplot_3_detail_new_2, data_predicted_finerplot_4_detail_new_2)))
range.val.c <- range(c(data, data_predicted0, data_predicted_finerplot_2, data_predicted_finerplot_3, data_predicted_finerplot_4))
range.val.d <- range(c(data_predicted_finerplot_2_detail_new_2, data_predicted_finerplot_3_detail_new_2, data_predicted_finerplot_4_detail_new_2))
palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black"))
palette.c<-colorRampPalette(c("cyan", "green", "yellow", "red", "black"))
palette.d<-colorRampPalette(rev(rainbow(5)))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_finerplot_2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=range.val.c, main="(a) Level 3, Coarser", cex.main=1.5, smallplot=c(0.75,0.8,0.6,0.85), coldefault=palette.c)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_2$pointsin,1], example_network@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_2$pointsin,2], pch=22)

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_finerplot_2_detail[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=range.val.d, main="(b) Level 3, Detail", cex.main=1.5, smallplot=c(0.75,0.8,0.6,0.85), coldefault=palette.d)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_2$pointsin,1], example_network@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_2$pointsin,2], pch=23)

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_finerplot_3[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=range.val.c, main="(c) Level 2, Coarser", cex.main=1.5, smallplot=c(0.75,0.8,0.6,0.85), coldefault=palette.c)
points(example_network_new@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_3$pointsin,1], example_network_new@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_3$pointsin,2], pch=22)

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_finerplot_3_detail[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=range.val.d, main="(d) Level 2, Detail", cex.main=1.5, smallplot=c(0.75,0.8,0.6,0.85), coldefault=palette.d)
#points(example_network@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_3$pointsin,1], example_network@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_3$pointsin,2], pch=23)
#points(example_network@obspoints@SSNPoints[[1]]@point.coords[setdiff(result_forward_finerplot_2$pointsin,result_forward_finerplot_3$pointsin),1], example_network@obspoints@SSNPoints[[1]]@point.coords[setdiff(result_forward_finerplot_2$pointsin,result_forward_finerplot_3$pointsin),2], pch=23)
points(example_network_new@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_3$pointsin,1], example_network_new@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_3$pointsin,2], pch=23)

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_finerplot_4[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=range.val.c, main="(e) Level 1, Coarser", cex.main=1.5, smallplot=c(0.75,0.8,0.6,0.85), coldefault=palette.c)
points(example_network_newnew@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_4$pointsin,1], example_network_newnew@obspoints@SSNPoints[[1]]@point.coords[result_forward_finerplot_4$pointsin,2], pch=22)

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_finerplot_4_detail[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=range.val.d, main="(f) Level 1, Detail", cex.main=1.5, smallplot=c(0.75,0.8,0.6,0.85), coldefault=palette.d)
#points(example_network@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_4$pointsin,1], example_network@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_4$pointsin,2], pch=23)
#points(example_network@obspoints@SSNPoints[[1]]@point.coords[setdiff(result_forward_finerplot_3$pointsin,result_forward_finerplot_4$pointsin),1], example_network@obspoints@SSNPoints[[1]]@point.coords[setdiff(result_forward_finerplot_3$pointsin,result_forward_finerplot_4$pointsin),2], pch=23)
points(example_network_newnew@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_4$pointsin,1], example_network_newnew@obspoints@SSNPoints[[1]]@point.coords[-result_forward_finerplot_4$pointsin,2], pch=23)

dev.off()
