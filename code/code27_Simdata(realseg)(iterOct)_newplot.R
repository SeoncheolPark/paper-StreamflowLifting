library(ggmap)
#thinning folder의 test_SpatioTemporal 함수 보고 작성
library(SpatioTemporal)
library(plotrix)
library(maps)
library(splitstackshape) #stratified sampling
library(riverdist)
library(rivernet)
library(xts)
library(zoo)
library(shape2graph)

library(expm)
library(geoR)
library(GWmodel)
library(matrixcalc)
library(missMDA)
library(sampling)

library(fields)
library(EbayesThresh)
library(waveslim)
library(adlift)
library(wavethresh)
source('~/Dropbox/R files/StreamFlow/sources/source_S.R', chdir = TRUE)
source('~/Dropbox/(논문들)/(논문)RIver Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/Me/source_Flexible.R', chdir = TRUE)

library(foreach)
library(parallel)
library(doParallel)
########################################
#금강 자료 출력
########################################
source('~/Dropbox/R files/StreamFlow/sources/source.R', chdir = TRUE)
#source('~/Dropbox/R files/StreamFlow/sources/source_ST.R', chdir = TRUE)
source('~/Dropbox/R files/StreamFlow/sources/source_S_nlt.R', chdir = TRUE)

library(smnet)
mouth_node <- 68

#길이 112개의 shape
shape_miho <- readOGR("~/Dropbox/R files/StreamFlow/River/GeumMiho.shp")

#기본 plotting
##MyRivernetwork: class "rivernetwork"
MyRivernetwork <- line2network(path="~/Dropbox/R files/StreamFlow/River/", layer="GeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
MyRivernetwork <- setmouth(seg=68, vert=20, rivers=MyRivernetwork)
par(mfrow=c(1,1))
plot(MyRivernetwork)

#get_binaryIDs_stream: source에 만들어놓은 함수: 뒤의 get_adjacency_stream 실행 위함
binaryIDs_obj <- get_binaryIDs_stream(mouth_node=68 , shape_miho)
#수동으로 adjacency를 얻도록 만들어진 함수
adjacency <- get_adjacency_stream(binaryIDs_obj)

#compute_shreve: shreve 거리를 재도록 직접 제작한 함수
shreve_order <- compute_shreve(adjacency)

library(SSN)
#for examples, copy MiddleFork04.ssn directory to R's temporary directory
copyLSN2temp()
mf04p <- importSSN(paste0(tempdir(),'/MiddleFork04.ssn'), predpts = "pred1km",o.write = TRUE)
mf04 <- importSSN_stream(shreve_obj = shreve_order, multipleplaces=TRUE)
#pid 변경
#mf04@obspoints@SSNPoints[[1]]@point.data$pid <- c(1:length(example_network@obspoints@SSNPoints[[1]]@point.data$pid))
#mf04@predpoints@SSNPoints[[1]]@point.data$pid <- length(example_network@obspoints@SSNPoints[[1]]@point.data$pid)+1
plot(mf04)

show_weights(mf04, adjacency)
#show_weights를 작동시키려면
#ims@obspoints@SSNPoints[[1]]@point.data 안에 shreve, addfunccol이 있어야 함
########################################
#data: 시계열 자료 업로드
########################################
#data <- readRDS("~/Dropbox/Data/RIverFlow/금강/Geum(miho).RDS")
data <- readRDS("~/Dropbox/Data/RIverFlow/금강/Geum(miho)(extended).RDS")
#data <- readRDS("~/Dropbox/Data/RIverFlow/금강/Geum(extended).RDS")
library(extrafont)
#font_import()
par(family="NanumGothic")

data$normaldata_TN <- data$normaldata_TN["2012/2017",]
normaldata_TN <- data$normaldata_TN

#data$normaldata_TN["2012-07/2012-08",]
sum(!is.na(data$normaldata_TN)) #총 4372개
sum(!is.na(data$normaldata_TN)) / (sum(!is.na(data$normaldata_TN))  + sum(is.na(data$normaldata_TN)) ) #전체 시계열의 약 7%뿐
apply(as.matrix(data$normaldata_TN), 2, function(x) sum(sum(!is.na(x)))) #각 장소별로 time series value들이 얼마나 있는지 살펴보도록 하자

#시계열 플롯
for(i in 1:ncol(data$normaldata_TN)){
  plot(log(as.numeric(data$normaldata_TN[,i])), main=colnames(data$normaldata_TN)[i], ylab="Nitrogen")
}

#(SSN obj의 장소 순서, 시계열 장소 순서 매치)
match(mf04@obspoints@SSNPoints[[1]]@point.data$X, colnames(data$normaldata_TN))
#(시계열 장소 순서, SSN obj의 장소 순서 매치)
match(colnames(data$normaldata_TN), mf04@obspoints@SSNPoints[[1]]@point.data$X)
#obsdata에 맞춰 정렬하였다
#data$normaldata <- data$normaldata[,match(mf04@obspoints@SSNPoints[[1]]@point.data$X, colnames(data$normaldata))]

#내 스스로 리버 network와 관찰 장소 위치를 출력하는 함수들
plot(mf04@obspoints@SSNPoints[[1]]@point.coords[,1], mf04@obspoints@SSNPoints[[1]]@point.coords[,2], type="n")
plot(MyRivernetwork)
text(x=mf04@obspoints@SSNPoints[[1]]@point.coords[,1], y=mf04@obspoints@SSNPoints[[1]]@point.coords[,2], labels=c(1:nrow(mf04@obspoints@SSNPoints[[1]]@point.coords)))
########################################
#data: 시계열 자료 업로드
########################################
data(mesa.data.raw, package="SpatioTemporal") #참고용
within_ID <- sort(mf04@obspoints@SSNPoints[[1]]@point.data$rid)
#Format
#The structure contains observations, temporal trends, locations, geographic covariates, and spatio-temporal covariates. The data is stored as a list with elements:
#X
#A data.frame containing names, locations, and (geographic) covariates for all the (observation) locations.
#obs
#A time-by-location matrix for the observed data, missing data marked as NA
#lax.conc.1500
#A time-by-location matrix of a spatio-temporal covariate based on output from Caline3QHC.

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
eventplace2 <- data$eventplace[which(!is.na(match(data$eventplace$X, colnames(data$normaldata_TN)))),] #28개 측정소 정보 출력
eventplace2 <- eventplace2[match(eventplace2$X, colnames(data$normaldata_TN)),]
row.names(eventplace2) <- eventplace2$X
time_index <- as.character(time(data$normaldata_TN["2012/2017",]))
normaldata2 <- as.matrix(data$normaldata_TN["2012/2017",])
row.names(normaldata2) <- time_index
stlistdata <- as.matrix(data$normaldata_temp["2012/2017"])
row.names(stlistdata) <- time_index

na.mismatch <- c()
for(i in 1:ncol(normaldata2)){
  na.mismatch <- c(na.mismatch, sum(is.na(match(which(!is.na(normaldata2[,i])), which(!is.na(stlistdata[,i]))))) )
}
#조천 2012년 9월 4일

eventplace2 <- eventplace2[match(rownames(eventplace2), colnames(normaldata2)),]
colnames(eventplace2)[1] <- "ID"
colnames(eventplace2)[16] <- "long"
colnames(eventplace2)[15] <- "lat"
#colnames(eventplace2)[4] <- "type"

#mesa.data$covars에 대한 애들 새로 만들어야 할듯
eventplace2 <- cbind(eventplace2, 
                     TN=apply(as.matrix(data$normaldata_TN["2012/2017"]), 2, mean, na.rm=T), #TN: 총질소
                     temp=apply(as.matrix(data$normaldata_TN["2012/2017"]), 2, mean, na.rm=T), #temp: 수온
                     ph=apply(as.matrix(data$normaldata_ph["2012/2017"]), 2, mean, na.rm=T), #ph: 수소이온농도
                     elec=apply(as.matrix(data$normaldata_elec["2012/2017"]), 2, mean, na.rm=T), #elec: 전기전도도
                     O2=apply(as.matrix(data$normaldata_O2["2012/2017"]), 2, mean, na.rm=T), #O2: 용존산소
                     BOD=apply(as.matrix(data$normaldata_BOD["2012/2017"]), 2, mean, na.rm=T), #BOD: BOD
                     COD=apply(as.matrix(data$normaldata_COD["2012/2017"]), 2, mean, na.rm=T), #COD: COD
                     susp=apply(as.matrix(data$normaldata_susp["2012/2017"]), 2, mean, na.rm=T), #susp: 부유물질
                     TP=apply(as.matrix(data$normaldata_TP["2012/2017"]), 2, mean, na.rm=T), #TP: 총인
                     TOC=apply(as.matrix(data$normaldata_TOC["2012/2017"]), 2, mean, na.rm=T) #, TOC: 총유기탄소
                     #apply(as.matrix(data$normaldata_flow["2012/2017"]), 2, mean, na.rm=T) #flow: 유량 (NaN도 있다)
) 
eventplace2 <- cbind(eventplace2, x=eventplace2$long, y=eventplace2$lat, type=ifelse(apply(normaldata2, 2, function(x) sum(!is.na(x))) > 100, "large", "small"))
########################################
##Start data analysis
########################################
example_network <- mf04
n.duplicated <- sum(duplicated(example_network@obspoints@SSNPoints[[1]]@point.data$rid))
library(adlift); library(nlt); library(liftLRD); library(CNLTreg); library(CNLTtsa); library(CliftLRD); library(stringr)
library(ggmap)
#내가 만든 코드
source('~/Dropbox/R files/StreamFlow/sources/source.R', chdir = TRUE)

#RID를 0부터 112까지 정렬하는 업데이트
for(i in 1:length(example_network@lines)){
  example_network@lines[[i]]@ID <- as.character(i-1)
}
#Obs pt의 coord를 rid index에 맞춰 정렬하는 업데이트
for(i in 1:nrow(example_network@obspoints@SSNPoints[[1]]@point.coords)){
  lineID <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[i]
  segID <- example_network@obspoints@SSNPoints[[1]]@point.data$segid[i]
  example_network@obspoints@SSNPoints[[1]]@point.coords[i,] <- example_network@lines[[lineID]]@Lines[[1]]@coords[segID,]
}


#원래 코드(틀린듯)
#colnames(data$normaldata_TN) <- as.character(example_network@obspoints@SSNPoints[[1]]@point.data$X)
#eventplace2 <- eventplace2[match(example_network@obspoints@SSNPoints[[1]]@point.data$X, eventplace2$ID), ]
#code_Flexible_me.R의 내용 가져옴
data$normaldata_TN <- data$normaldata_TN[,match(example_network@obspoints@SSNPoints[[1]]@point.data$X, colnames(data$normaldata_TN))]
eventplace2 <- eventplace2[match(example_network@obspoints@SSNPoints[[1]]@point.data$X, eventplace2$ID), ]
########################################
#이제부터는 flexible smoothing
########################################
# Load up all the necessary packages and helper functions
#setwd("~/Dropbox/(논문들)/(논문)SpatExtreme관련/Extremes on River Network/Flexible regression models over river networks(2014)_code/")
library(mgcv); library(spam); library(RgoogleMaps)
library(lubridate)
# library(Matrix)
library(scales)

source('~/Dropbox/(논문들)/(논문)River Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/paper_functions.r', chdir = TRUE)

########################################
## Data generation
########################################
#realweights #length 113 vector
adjacency_old <- adjacency

upper_seg <- which(colSums(adjacency_old$adjacency)==0) #이점들에서만 생성하면 됨 (총 55개)
upper_seg_list <- list()
upper_seg_label <- rep("A", 113)
upper_seg_list[[1]] <- c(4, 5, 6)
upper_seg_label[upper_seg_list[[1]]] <- "B"
upper_seg_list[[2]] <- c(7, 8, 14, 15, 16, 25, 26, 31, 32)
upper_seg_label[upper_seg_list[[2]]] <- "C"
upper_seg_list[[3]] <- c(21, 23, 24, 111)
upper_seg_label[upper_seg_list[[3]]] <- "D"
upper_seg_list[[4]] <- c(28, 29, 37, 39, 41, 42, 44, 47)
upper_seg_label[upper_seg_list[[4]]] <- "E"
upper_seg_list[[5]] <- c(33, 34, 35)
upper_seg_label[upper_seg_list[[5]]] <- "F"
upper_seg_list[[6]] <- c(51, 91, 93, 95, 97, 99, 101)
upper_seg_label[upper_seg_list[[6]]] <- "G"
upper_seg_list[[7]] <- c(52, 53)
upper_seg_label[upper_seg_list[[7]]] <- "H"
upper_seg_list[[8]] <- c(58, 61, 102, 104, 106, 108, 110)
upper_seg_label[upper_seg_list[[8]]] <- "I"
upper_seg_list[[9]] <- c(59, 60, 66, 67)
upper_seg_label[upper_seg_list[[9]]] <- "J"
upper_seg_list[[10]] <- c(63); upper_seg_list[[11]] <- c(27)
upper_seg_label[upper_seg_list[[10]]] <- "H"; upper_seg_label[upper_seg_list[[11]]] <- "H"
#upper_seg_list[[12]] <- c(64, 65, 79); upper_seg_list[[13]] <- c(69, 70, 75)
#upper_seg_label[upper_seg_list[[12]]] <- "K"; upper_seg_label[upper_seg_list[[13]]] <- "L"
upper_seg_list[[12]] <- c(64); upper_seg_list[[13]] <- c(65); upper_seg_list[[14]] <- c(79)
upper_seg_label[upper_seg_list[[12]]] <- "H"; upper_seg_label[upper_seg_list[[13]]] <- "H"; upper_seg_label[upper_seg_list[[14]]] <- "H"
upper_seg_list[[15]] <- c(70); upper_seg_list[[16]] <- c(75); upper_seg_list[[17]] <- c(69)
upper_seg_label[upper_seg_list[[15]]] <- "H"; upper_seg_label[upper_seg_list[[16]]] <- "H"; upper_seg_label[upper_seg_list[[17]]] <- "H"

#label 다시 매기기
upper_seg_label <- rep("A", 113)
upper_seg_label[c(1,2,3,4,5,6,112)] <- "B"
upper_seg_label[c(28,29,30,37,38,39,40,41,42,43,44,45,49)] <- "C"
upper_seg_label[c(20,21,22,23,24,87,111)] <- "D"
upper_seg_label[c(33,34,35,36,85)] <- "E"
upper_seg_label[c(58,61,62,83,102,103,104,105,106,107,108,109,110)] <- "F"
upper_seg_label[c(7,8,9,10,11,12,13,14,15,16,17,18,19,25,26,31,32)] <- "G"
upper_seg_label[c(51,90,91,92,93,94,95,96,97,98,99,100,101)] <- "H"
upper_seg_label[c(55,56,57,59,60,66,67)] <- "I"
upper_seg_label[c(52,53,54)] <- "J"
upper_seg_label[c(50,63,27,69,70,75,64,65,79)] <- "K"




weight_vec_candidate <- 1.5*exp(-(example_network@network.line.coords$DistanceUpstream-min(example_network@network.line.coords$DistanceUpstream))/0.4013373)
weight_vec_candidate2 <- compute_shreve_and_dist(adjacency_old, example_network,type="lengthprop")
weight_vec_candidate2 <- scales::rescale(log(sqrt(weight_vec_candidate2$distweight)), to=c(0.2, 1.5))
#weight_vec_candidate2 <- compute_shreve_and_dist(adjacency_old, example_network, scalevec=c(0.2, 1.5),logdata=TRUE,type="lengthprop")
weight_vec_candidate <- weight_vec_candidate2


#plotting (plotting을 위해서는 꼭 TweedPredPoints를 generate해야 한다)
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

#location_imsi
location_imsi <- c()
distupstream_imsi <- c()
for(i in 1:length(example_network@lines)){
  if(nrow(example_network@lines[[i]]@Lines[[1]]@coords)==2){
    location_imsi <- rbind(location_imsi, colMeans(example_network@lines[[i]]@Lines[[1]]@coords))
    distupstream_imsi <- rbind(distupstream_imsi, riverdistance(startseg=i, endseg=68, startvert=1, endvert=20, rivers=MyRivernetwork))
  }else{
    index.imsi2 <- (nrow(example_network@lines[[i]]@Lines[[1]]@coords)/2)-1/2
    location_imsi <- rbind(location_imsi, colMeans(example_network@lines[[i]]@Lines[[1]]@coords[c(index.imsi2 , index.imsi2+1),]))
    distupstream_imsi <- rbind(distupstream_imsi, riverdistance(startseg=i, endseg=68, startvert=index.imsi2, endvert=20, rivers=MyRivernetwork))
  }
}
location_imsi_true <- c()
for(i in 1:nrow(example_network@obspoints@SSNPoints[[1]]@point.data)){
  location_imsi_true <- rbind(location_imsi_true, example_network@lines[[example_network@obspoints@SSNPoints[[1]]@point.data$rid[[i]]]]@Lines[[1]]@coords[example_network@obspoints@SSNPoints[[1]]@point.data$segid[i],])
}

par(family = 'sans') 
par(mar=c(1.1,1.1,1.1,1.1))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="n", ylab="n", xaxt="n", yaxt="n")
points(location_imsi[-upper_seg,1], location_imsi[-upper_seg,2], pch=16, col="orange")
points(location_imsi[upper_seg,1], location_imsi[upper_seg,2], pch=15, col="purple")
#quilt.plot(location_imsi[,1], location_imsi[,2], colMeans(result_signal), cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
library(mixtools)
for(iii in c(1,2,3,4,5,6,8,9)){
  ellipse(mu=c(mean(location_imsi[upper_seg_list[[iii]],1]), mean(location_imsi[upper_seg_list[[iii]],2])), sigma=var(location_imsi[upper_seg_list[[iii]],])/1.5, npoints=200, newplot=FALSE, col="red")
}
library(plotrix)
for(iii in c(7, 10:17)){
  if(iii==7){
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.02, border="red" )
  }else{
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.005, border="red" )
  }
}
legend("topleft", pch=c(15,16), col=c("purple", "orange"), c("Upper-most segments", "Non-upper-most segments"), bty="n")

par(family = 'sans') 
par(mar=c(1.1,1.1,2.1,1.1))
par(mfrow=c(1,2))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="n", ylab="n", xaxt="n", yaxt="n", main="(a)")
points(location_imsi[-upper_seg,1], location_imsi[-upper_seg,2], pch=21, bg="orange", col="black")
points(location_imsi[upper_seg,1], location_imsi[upper_seg,2], pch=22, bg="purple", col="black")
#quilt.plot(location_imsi[,1], location_imsi[,2], colMeans(result_signal), cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
library(mixtools)
for(iii in c(1,2,3,4,5,6,8,9)){
  ellipse(mu=c(mean(location_imsi[upper_seg_list[[iii]],1]), mean(location_imsi[upper_seg_list[[iii]],2])), sigma=var(location_imsi[upper_seg_list[[iii]],])/1.5, npoints=200, newplot=FALSE, col="red")
}
library(plotrix)
for(iii in c(7, 10:17)){
  if(iii==7){
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.02, border="red" )
  }else{
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.01, border="red" )
  }
}
legend("topleft", pch=c(15,16), col=c("purple", "orange"), c("Upper-most", "Non-upper-most"), bty="n")

unique_upper_seg_label <- unique(upper_seg_label)
rainbow_palette <- rainbow(length(unique_upper_seg_label))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=rainbow_palette[match(upper_seg_label[TweedPredPoints$StreamUnit], unique_upper_seg_label)], pch=16, cex=TweedPredPoints$Weights, xlab="n", ylab="n", xaxt="n", yaxt="n", main="(b)")
points(location_imsi[-upper_seg,1], location_imsi[-upper_seg,2], pch=21, bg="orange", col="black")
points(location_imsi[upper_seg,1], location_imsi[upper_seg,2], pch=22, bg="purple", col="black")
legend("topleft", pch=c(15,16), col=c("purple", "orange"), c("Upper-most", "Non-upper-most"), bty="n")


#4 5 6
#7 8 14 15 16 25 26 31 32
#21 23 24 111
#27
#28 29 37 39 41 42 44 47
#33 34 35 
#51 91 93 95 97 99 101
#52 53
#58 61 102 104 106 108 110
#59 60 66 67
#63

#64
#65
#69
#70
#75
#79

obs_pt_mat <- c()
for(i in 1:length(example_network@lines)){
  if(nrow(example_network@lines[[i]]@Lines[[1]]@coords)==2){
    obs_pt_mat <- rbind(obs_pt_mat, colMeans(example_network@lines[[i]]@Lines[[1]]@coords))
  }else{
    median_row <- round(median(nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
    obs_pt_mat <- rbind(obs_pt_mat, example_network@lines[[i]]@Lines[[1]]@coords[median_row,])
  }
}
#plot(obs_pt_mat )

example_network3 <- example_network
example_network3@obspoints@SSNPoints[[1]]@point.coords <- location_imsi
for(i in 1:nrow(location_imsi)){
  example_network3@obspoints@SSNPoints[[1]]@point.data <- rbind(example_network3@obspoints@SSNPoints[[1]]@point.data,
                                                                c(NA,NA,NA,NA,NA,NA,NA,NA,rep(NA,6), location_imsi[i,2], location_imsi[i,1],
                                                                  NA,NA,NA,rep(NA,9), NA, i, distupstream_imsi[i], i, 1, i, shreve_order[i], NA, i, NA    ))
  example_network3@obspoints@SSNPoints[[1]]@network.point.coords <- rbind( example_network3@obspoints@SSNPoints[[1]]@network.point.coords, c(1,i, distupstream_imsi[i]))
}
example_network3@obspoints@SSNPoints[[1]]@point.data <- example_network3@obspoints@SSNPoints[[1]]@point.data[-c(1:28),]
example_network3@obspoints@SSNPoints[[1]]@network.point.coords <-  example_network3@obspoints@SSNPoints[[1]]@network.point.coords[-c(1:28),]



#TweedPredPoints_old <- data.frame(Latitude=c(), Longitude=c(), StreamUnit=c(), Weights=c())
#for(jj in 1:length(example_network@lines)){
#  #if(nrow(example_network@lines[[jj]]@Lines[[1]]@coords)>2){
#  #  seg_cand <- example_network@lines[[jj]]@Lines[[1]]@coords[-nrow(example_network@lines[[jj]]@Lines[[1]]@coords),]
#  #}else{
#  seg_cand <- example_network@lines[[jj]]@Lines[[1]]@coords
#  #}
#  TweedPredPoints_old_imsi <- data.frame(Latitude=seg_cand[,2], Longitude=seg_cand[,1], 
#                                         StreamUnit=rep(example_network@data$rid[jj], nrow(seg_cand)), 
#                                         Weights=rep(weight_vec_candidate[jj], nrow(seg_cand)))
#  TweedPredPoints_old <- rbind(TweedPredPoints_old, TweedPredPoints_old_imsi)
#}
#TweedPredPoints <- TweedPredPoints_old
n.seg <- 113
#n.sub.length <- 59
n.sub.length <- 40
n.sub.length.true <- 28
sd.val <- 1


#adjacency_old 작성
adjacency_old <- get_adjacency_stream(binaryIDs_obj)
adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))


realweights <- as.matrix(weight_vec_candidate, ncol=1)

n.iter <- 1#00
result_mat <- matrix(0, nrow=n.iter, ncol = 4)

result_list <- list()

for(ijk in 1:n.iter){
  ########여기서부터 실행
  #data generation
  print(paste("Iteration ", ijk))
  
  n <- 1; 
  #x <- sort((sample(1:(365*3), n)/(365*3)))
  x <- 1/5 #임의의 시간대
  #generate result_fct.1
  result_fct.1 <- matrix(0, nrow=n, ncol=n.seg)
  for(i in 1:n.seg){
    #if(i%in%upper_seg){
    result_fct.1[,i] <- sin(-x*(2.5*pi))
    #}
  }
  
  #generate result_fct.2
  
  
  result_fct.2 <- matrix(0, nrow=n, ncol=n.seg)
  spat.signal.imsi.orig <- c(3,6,9)
  while(TRUE){
    #spat.signal.imsi <- c(3,4,5)
    spat.signal.imsi <- spat.signal.imsi.orig
    set_length <- c(1:length(upper_seg_list))
    spat.signal.imsi.index <- sample(set_length, length(spat.signal.imsi))
    for(i in 1:length(spat.signal.imsi.index)){
      #result_fct.2[,upper_seg_list[[spat.signal.imsi.index[i]]]] <- rnorm(length(spat.signal.imsi.index[i]), spat.signal.imsi[i], 0.5)
      result_fct.2[,upper_seg_list[[spat.signal.imsi.index[i]]]] <- spat.signal.imsi[i]
    }
    if(sum(result_fct.2!=0)>=30){
      break
    }else{
      set_length <- setdiff(set_length, spat.signal.imsi.index)
      spat.signal.imsi <- spat.signal.imsi.orig[which.min(sapply(spat.signal.imsi.orig, function(x) sum(x==result_fct.2)))]
    }
  }
  
  
  #add noise
  noise_mat <- matrix(rnorm(n*n.seg, 0, sd.val),n,n.seg)
  
  
  result_signal_partial <- result_fct.2 #+ result_fct.3
  lower_seg <- setdiff(1:n.seg,upper_seg)
  
  
  lower_seg <- which(colSums(adjacency_old$adjacency)!=0)
  while(length(lower_seg)>0){
    lower_seg_delete <- c()
    for(i in 1:length(lower_seg)){
      if(sum(!(which(adjacency_old$adjacency[,lower_seg[i]]!=0) %in% lower_seg))==length(which(adjacency_old$adjacency[,lower_seg[i]]!=0) %in% lower_seg)){
        if(length(which(adjacency_old$adjacency[,lower_seg[i]]!=0) %in% lower_seg)==1){
          result_signal_partial[,lower_seg[i]] <- result_signal_partial[,which(adjacency_old$adjacency[,lower_seg[i]]!=0)]
        }else{
          #result_signal_partial[,lower_seg[i]] <- rowSums(result_signal[,which(adjacency_old[,lower_seg[i]]!=0)]) 
          result_signal_partial[,lower_seg[i]] <- (result_signal_partial[,which(adjacency_old$adjacency[,lower_seg[i]]!=0)])%*%(weight_vec_candidate[which(adjacency_old$adjacency[,lower_seg[i]]!=0)]/sum(weight_vec_candidate[which(adjacency_old$adjacency[,lower_seg[i]]!=0)]))
        }
        lower_seg_delete <- c(lower_seg_delete,lower_seg[i])
      }
    }
    lower_seg <- setdiff(lower_seg, lower_seg_delete )
  }
  result_signal <- result_signal_partial + result_fct.1 + 10 #10: 로그변환 했을 시 음수가 되는 것을 방지
  while(TRUE){
    result_data_noisy <- result_signal + noise_mat
    if(min(result_data_noisy)>0){
      break
    }else
      noise_mat <- matrix(rnorm(n*n.seg, 0, sd.val),n,n.seg)
  }
  
  TweedData <- c()
  #TweedData generation (실제 자료의 duplication 반영)
  TweedData <- data.frame(date=as.Date("2012-12-31")+x*(365*3),
                          location=1,
                          nitrate=result_data_noisy[,1],
                          #nitrate=result_signal[,example_network@obspoints@SSNPoints[[1]]@point.data$rid[1]] + rnorm(1, mean=0, sd.val),
                          Long=rep(location_imsi[1,1],n),
                          #Long=rep(location_imsi_true[1,1],n),
                          Lat=rep(location_imsi[1,2],n)
                          #Lat=rep(location_imsi_true[1,2],n)
  )
  for(i in 2:n.seg){
    dataframe.new <- data.frame(date=as.Date("2012-12-31")+x*(365*3),
                                location=i,
                                nitrate=result_data_noisy[,i],
                                #nitrate=result_signal[,example_network@obspoints@SSNPoints[[1]]@point.data$rid[i]] + rnorm(1, mean=0, sd.val),
                                Long=rep(location_imsi[i,1],n),
                                #Long=rep(location_imsi_true[i,1],n),
                                Lat=rep(location_imsi[i,2],n)
                                #Lat=rep(location_imsi_true[i,2],n)
    )
    TweedData <- rbind(TweedData, dataframe.new)
  }
  
  
  reduce.data <- FALSE
  if(reduce.data==TRUE){
    TweedData <- TweedData[which(TweedData$location%in%within_ID),]
    example_network@obspoints@SSNPoints[[1]]@point.data <- example_network@obspoints@SSNPoints[[1]]@point.data[within_ID,]
    example_network@obspoints@SSNPoints[[1]]@network.point.coords <-  example_network@obspoints@SSNPoints[[1]]@network.point.coords[within_ID,]
  }
  
  
  data = sapply(split(TweedData, f = TweedData$location), function(x) mean(x$nitrate)) #꼭 실행시켜야
  
  #names(data) <- example_network@obspoints@SSNPoints[[1]]@point.data$rid
  #TweedData$location <- example_network@obspoints@SSNPoints[[1]]@point.data$rid
  #data <- data[match(example_network@obspoints@SSNPoints[[1]]@point.data$rid, names(data))]
  #par(mfrow=c(1,2))
  #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights)
  #quilt.plot(location_imsi[,1], location_imsi[,2], data, main="(a) Raw", cex=2, add=T, zlim=range(data), xlab="x", ylab="y")
  
  #sapply(split(TweedData$nitrate, f=TweedData$location), function(x) median(x))
  
  
  #추가: 공간자료만으로?
  #TweedData.new = TweedData[seq(1, nrow(TweedData), 5), ]
  #TweedData.new$nitrate <- data
  
  #만약 subindex를 쓰려면
  index_sub <- sort(c(sample(setdiff(c(1:length(example_network3@obspoints@SSNPoints[[1]]@point.data$rid)),68), n.sub.length -1, replace = FALSE),68))
  #strified sampling 활용
  #index_sub_data <- data.frame(stations=setdiff(c(1:length(example_network@obspoints@SSNPoints[[1]]@point.data$rid)),68), groups=as.factor(result_fct.2[-68]))
  index_sub_data_new <- data.frame(stations=c(1:length(example_network3@obspoints@SSNPoints[[1]]@point.data$rid)), groups=as.factor(upper_seg_label))
  index_sub_data <- data.frame(stations=setdiff(c(1:length(example_network3@obspoints@SSNPoints[[1]]@point.data$rid)),68), groups=as.factor(upper_seg_label[-68]))
  index_sub <- sort(c(stratified(indt=index_sub_data, group="groups", size=(n.sub.length -1)/112)$stations, 68))
  
  #그냥 하고싶으면:
  if(n.sub.length==113){
    index_sub <- c(1:length(example_network3@obspoints@SSNPoints[[1]]@point.data$rid))
  }
  
  index_choose <- index_sub
  
  TweedData <- TweedData[index_sub,]
  #그냥 하고싶으면:
  #index_sub <- c(1:length(example_network@obspoints@SSNPoints[[1]]@point.data$rid))
  
  #(1) smnet_ST function
  penalties_default <- c(50, 50, c(25, 5), 50, c(25, 5), c(25, 25))
  division_factor <- c(30, 25, 20, 15, 10, 5, 1, 0.5, 0.25, 0.125)
  AICc_vec <- rep(0, length(division_factor))
  for(jjj in 1:length(division_factor)){
    result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/division_factor[jjj], plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=FALSE, model.type = c("c", "m"))
    AICc_vec[jjj] <- result_a$AICc
  }
  #result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=TRUE, model.type = c("c", "m"))
  #result_a$AICc
  #result_b <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/5, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=TRUE, model.type = c("c", "m"))
  #result_b$AICc
  #result_c <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/25, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=TRUE, model.type = c("c", "m"))
  #result_c$AICc
  #result_d <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/30, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=TRUE, model.type = c("c", "m"))
  #result_d$AICc
  
  #opt_result <- optim(par=penalties_default/55, fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/10, upper=penalties_default/1, log.y=TRUE, model.type = c("c", "m") )
  #opt_val <- c(2.5, 2.0, 1.0, 0.2, 2.0, 1.0, 0.2, 1.0, 1.0)
  #opt_val <- c(5.00002, 5.00000, 2.50000, 0.50000, 5.00000, 2.50000, 0.50000, 2.50000, 2.50000)
  
  AICc_vec_min_index <- which.min(AICc_vec)
  
  if(AICc_vec_min_index==1){
    opt_result <- optim(par=penalties_default/division_factor[AICc_vec_min_index], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index)], upper=penalties_default/division_factor[(AICc_vec_min_index+1)], log.y=FALSE, model.type = c("c", "m") )
  }else if(AICc_vec_min_index==length(AICc_vec)){
    opt_result <- optim(par=penalties_default/division_factor[AICc_vec_min_index], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index)-1], upper=penalties_default/division_factor[(AICc_vec_min_index)], log.y=FALSE, model.type = c("c", "m") )
  }else{
    opt_result <- optim(par=penalties_default/division_factor[AICc_vec_min_index], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index-1)], upper=penalties_default/division_factor[(AICc_vec_min_index+1)], log.y=FALSE, model.type = c("c", "m") )
  }
  opt_val <- opt_result$par
  
  #opt_result <- optim(par=penalties_default/25, fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/30, upper=penalties_default/20, log.y=FALSE,model.type = c("c", "m"))
  #result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties= opt_result$par, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=FALSE,model.type = c("c", "m") )
  result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties= opt_result$par, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=FALSE,model.type = c("c", "m") )
  unique(result_a$fit)
  sapply(split(result_a$fit, f=result_a$TweedData$location), function(x) median(x)) #mean으로 바꿔도 결과 같다
  
  
  example_network2 <- example_network3
  #example_network2@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(example_network@obspoints@SSNPoints[[1]]@network.point.coords, example_network3@obspoints@SSNPoints[[1]]@network.point.coords[index_choose,])
  #example_network2@obspoints@SSNPoints[[1]]@point.coords <- rbind(example_network@obspoints@SSNPoints[[1]]@point.coords, example_network3@obspoints@SSNPoints[[1]]@point.coords[index_choose,])
  #example_network2@obspoints@SSNPoints[[1]]@point.data <- rbind(example_network@obspoints@SSNPoints[[1]]@point.data, example_network3@obspoints@SSNPoints[[1]]@point.data[index_choose,])
  example_network2@obspoints@SSNPoints[[1]]@network.point.coords <- example_network3@obspoints@SSNPoints[[1]]@network.point.coords[index_sub,]
  example_network2@obspoints@SSNPoints[[1]]@point.coords <- example_network3@obspoints@SSNPoints[[1]]@point.coords[index_sub,]
  example_network2@obspoints@SSNPoints[[1]]@point.data <- example_network3@obspoints@SSNPoints[[1]]@point.data[index_sub,]
  
  #(2) stream-flow lifting
  result_denoise <- denoise_S(data[index_sub], example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", returnall = FALSE)
  result_denoise2 <- denoise_S(data[index_sub], example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "hard", returnall = FALSE)
  
  result_forward <- fwtnp_stream_S(data[index_sub], example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
  
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  #denoise_Stream_S_perm(data, example_network2, endpt=1, per=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE, plot.fig = FALSE, plot.individual = FALSE, pollutant=NULL, polluyear=NULL, plot.thesis=FALSE)
  #result_denoise3 <- nlt_Stream_S(data, example_network2, J=10, endpt=1, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = TRUE)$aveghat
  result_denoise3 <- nlt_Stream_S(data[index_sub], example_network2, J=10, endpt=68, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = TRUE, ga=TRUE, index_sub_data=index_sub_data_new, oremovelist=result_forward$removelist)$aveghat
  
  stopCluster(cl)
  
  #adjacency=adjacency_old; pred = streamPred_S; neigh = 1; int = TRUE; clo = FALSE; keep = 35; rule = "hard"; returnall = FALSE
  #plot(unique(result_a$fit), result_denoise)
  
  #location_imsi_true_new <- rbind(location_imsi_true, location_imsi[index_choose,])
  
  par(family = 'sans') 
  par(mar=c(1.1,1.1,2.1,1.1))
  par(mfrow=c(2,3))
  
  zlims <-range(c(data, TweedData$nitrate, result_a$fit, result_denoise, result_denoise2, result_denoise3, colMeans(result_signal)))+c(-0.5,0.5)
  scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,colMeans(result_signal)[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(a) True", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
  
  scatter_fill(TweedPredPoints$Longitude[which(TweedPredPoints$StreamUnit%in%index_choose)], TweedPredPoints$Latitude[which(TweedPredPoints$StreamUnit%in%index_choose)] ,TweedData$nitrate[match(TweedPredPoints$StreamUnit[which(TweedPredPoints$StreamUnit%in%index_choose)], index_choose)], pch=16, cex=TweedPredPoints$Weights[which(TweedPredPoints$StreamUnit%in%index_choose)], xlab="", ylab="", xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude), zlim=zlims, main="(b) Observed", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
  points(TweedPredPoints$Longitude[-which(TweedPredPoints$StreamUnit%in%index_choose)], TweedPredPoints$Latitude[-which(TweedPredPoints$StreamUnit%in%index_choose)] ,col="gray", pch=16, cex=TweedPredPoints$Weights[-which(TweedPredPoints$StreamUnit%in%index_choose)])
  
  #old version
  #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(a) True", xlab="n", ylab="n", xaxt="n", yaxt="n")
  #quilt.plot(location_imsi[,1], location_imsi[,2], colMeans(result_signal), cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
  #
  #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(b) Observed", xlab="n", ylab="n", xaxt="n", yaxt="n")
  #quilt.plot(location_imsi[index_sub,1], location_imsi[index_sub,2], TweedData$nitrate, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")

  result_new_mat <- matrix(0, nrow=4, ncol=5)
  rownames(result_new_mat) <- c( "ODonnell", "S-Lifting(M)",  "S-Lifting(H)", "S-Lifting(N)")
  colnames(result_new_mat) <- c("Corr", "RMSE(whole)","WRMSE(whole)", "RMSE(actual)", "WRMSE(actual)")
  
  if(length(index_sub)==113){
    #예측 단계 불필요
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_a$fit[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(c) O'Donnell", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(d) Proposed (M)", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(e) Proposed (H)", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise3[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(f) Proposed (N)", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(c) O'Donnell", xlab="n", ylab="n", xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[index_sub,1], location_imsi[index_sub,2], result_a$fit, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(d) Proposed (M)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[index_sub,1], location_imsi[index_sub,2], result_denoise, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(e) Proposed (H)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[index_sub,1], location_imsi[index_sub,2], result_denoise2, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(f) Proposed (N)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[index_sub,1], location_imsi[index_sub,2], result_denoise3, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    
    result_mat[ijk,1] <-  sqrt(sum((colMeans(result_signal)-sapply(split(result_a$fit, f=result_a$TweedData$location), function(x) median(x)) )^2)/113)
    print(paste("O'Donnell ", round(sqrt(sum((colMeans(result_signal)-sapply(split(result_a$fit, f=result_a$TweedData$location), function(x) median(x)) )^2)/113),4)))
    result_mat[ijk,2] <- sqrt(sum((colMeans(result_signal)-result_denoise)^2)/113)
    print(paste("Proposed ", round(sqrt(sum((colMeans(result_signal)-result_denoise)^2)/113),4)))
    result_mat[ijk,3] <- sqrt(sum((colMeans(result_signal)-result_denoise2)^2)/113)
    print(paste("Proposed (hard) ", round(sqrt(sum((colMeans(result_signal)-result_denoise2)^2)/113),4)))
    result_mat[ijk,4] <- sqrt(sum((colMeans(result_signal)-result_denoise3)^2)/113)
    print(paste("Proposed (nlt, median) ", round(sqrt(sum((colMeans(result_signal)-result_denoise3)^2)/113),4)))
    
    result_new_mat[1,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_a$fit)) #O'Donnell
    result_new_mat[2,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_denoise)) #S-Lifting(M)
    result_new_mat[3,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_denoise2)) #S-Lifting(H)
    result_new_mat[4,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_denoise3)) #S-Lifting(N)
    
    result_new_mat[1,2] = result_new_mat[1,4] <- result_mat[ijk,1]
    result_new_mat[2,2] = result_new_mat[2,4] <- result_mat[ijk,2]
    result_new_mat[3,2] = result_new_mat[3,4] <- result_mat[ijk,3]
    result_new_mat[4,2] = result_new_mat[4,4] <- result_mat[ijk,4]
    
    result_new_mat[1,3] = result_new_mat[1,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_a$fit))^2)/n.seg)
    result_new_mat[2,3] = result_new_mat[2,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise))^2)/n.seg)
    result_new_mat[3,3] = result_new_mat[3,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise2))^2)/n.seg)
    result_new_mat[4,3] = result_new_mat[4,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise3))^2)/n.seg)
  }else{
    #예측 단계 필요
    initS_obj = initS_stream(X=as.row(as.numeric(names(data)[index_sub])), data=as.vector(data[index_sub]), example_network2, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data[index_sub]), 1, length(data[index_sub])))
    result_a_new <- initS_obj$weight_matrix%*%as.column(result_a$fit)
    result_denoise_new <- initS_obj$weight_matrix%*%as.column(result_denoise)
    result_denoise_new2 <- initS_obj$weight_matrix%*%as.column(result_denoise2)
    result_denoise_new3 <- initS_obj$weight_matrix%*%as.column(result_denoise3)
    
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_a_new[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(c) O'Donnell", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise_new[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(d) Proposed (M)", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise_new2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(e) Proposed (H)", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise_new3[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(f) Proposed (N)", cex.main=1.5, smallplot=c(0.15,0.2,0.65,0.85), y.axes = FALSE, y.axes.label=FALSE)
    
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(c) O'Donnell", xlab="n", ylab="n", xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[,1], location_imsi[,2], result_a_new , cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(d) Proposed (M)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[,1], location_imsi[,2], result_denoise_new, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(e) Proposed (H)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[,1], location_imsi[,2], result_denoise_new2, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(f) Proposed (N)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_imsi[,1], location_imsi[,2], result_denoise_new3, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    
    result_mat[ijk,1] <-  sqrt(sum((colMeans(result_signal)-result_a_new)^2)/113)
    print(paste("O'Donnell ", round(sqrt(sum((colMeans(result_signal)-result_a_new)^2)/113),4)))
    result_mat[ijk,2] <- sqrt(sum((colMeans(result_signal)-result_denoise_new)^2)/113)
    print(paste("Proposed ", round(sqrt(sum((colMeans(result_signal)-result_denoise_new)^2)/113),4)))
    result_mat[ijk,3] <- sqrt(sum((colMeans(result_signal)-result_denoise_new2)^2)/113)
    print(paste("Proposed (hard) ", round(sqrt(sum((colMeans(result_signal)-result_denoise_new2)^2)/113),4)))
    result_mat[ijk,4] <- sqrt(sum((colMeans(result_signal)-result_denoise_new3)^2)/113)
    print(paste("Proposed (nlt, median) ", round(sqrt(sum((colMeans(result_signal)-result_denoise_new3)^2)/113),4)))
    
    result_new_mat[1,1] <- cor(colMeans(result_signal), result_a_new) #O'Donnell
    result_new_mat[2,1] <- cor(colMeans(result_signal), result_denoise_new) #S-Lifting(M)
    result_new_mat[3,1] <- cor(colMeans(result_signal), result_denoise_new2) #S-Lifting(H)
    result_new_mat[4,1] <- cor(colMeans(result_signal), result_denoise_new3) #S-Lifting(N)
    
    result_new_mat[1,2] <- result_mat[ijk,1]
    result_new_mat[2,2] <- result_mat[ijk,2]
    result_new_mat[3,2] <- result_mat[ijk,3]
    result_new_mat[4,2] <- result_mat[ijk,4]
    
    result_new_mat[1,3] <- sqrt(sum(realweights*(colMeans(result_signal)-result_a_new)^2)/80)
    result_new_mat[2,3] <- sqrt(sum(realweights*(colMeans(result_signal)-result_denoise_new)^2)/80)
    result_new_mat[3,3] <- sqrt(sum(realweights*(colMeans(result_signal)-result_denoise_new2)^2)/80)
    result_new_mat[4,3] <- sqrt(sum(realweights*(colMeans(result_signal)-result_denoise_new3)^2)/80)
    
    result_new_mat[1,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_a$fit)^2)/length(index_choose))
    result_new_mat[2,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise)^2)/length(index_choose))
    result_new_mat[3,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise2)^2)/length(index_choose))
    result_new_mat[4,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise3)^2)/length(index_choose))
    
    result_new_mat[1,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-as.numeric(result_a$fit))^2)/length(index_choose))
    result_new_mat[2,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise))^2)/length(index_choose))
    result_new_mat[3,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise2))^2)/length(index_choose))
    result_new_mat[4,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise3))^2)/length(index_choose))
  }
  
  ##If we want to analyze it on the quantile level,
  result_list[[ijk]] <- result_new_mat
}
#saveRDS(result_mat, "StreamSim115(sd1).RDS")
#saveRDS(result_mat, "StreamSim80(sd2)nlt.RDS")

saveRDS(result_list, "ListStreamSim40(sd1)nlt.RDS")

#evaluation
aaaa <- readRDS("~/Dropbox/R files/ListStreamSim40(sd1)nlt.RDS")
bbbb <- array(as.numeric(unlist(aaaa)), dim=c(4,5,100))
mean(bbbb[1,2,])
mean(bbbb[2,2,])
mean(bbbb[3,2,])
mean(bbbb[4,2,])
sd(bbbb[1,2,])
sd(bbbb[2,2,])
sd(bbbb[3,2,])
sd(bbbb[4,2,])

