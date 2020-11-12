#install.packages("~/Dropbox/(논문들)/(논문)River Network/Flow-directed PCA(2017)_stpca_Rpackage/stpca_0.1.tar.gz", repos=NULL, type="source")
library(stpca)
library(SSN)
library(adlift); library(nlt); library(liftLRD); library(CNLTreg); library(CNLTtsa); library(CliftLRD); library(stringr)
library(ggmap)
library(SpatioTemporal)
library(plotrix)
library(maps)
library(smnet)
library(rgdal) #readOGR function
library(riverdist)
library(splitstackshape) #stratified sampling
library(foreach)
library(parallel)
library(doParallel)

#내가 만든 코드
source('~/Dropbox/R files/StreamFlow/sources/source.R', chdir = TRUE)
source('~/Dropbox/R files/StreamFlow/sources/source_S.R', chdir = TRUE)
source('~/Dropbox/R files/StreamFlow/sources/source_S_nlt.R', chdir = TRUE)

#source('~/Dropbox/(논문들)/(논문)RIver Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/Me/source_Flexible.R', chdir = TRUE)
source('~/Dropbox/R files/StreamFLow/sources/source_Flexible.R', chdir = TRUE)
########################################
#stpca package vignette
########################################
#결과를 모두 저장해 두었다
#vignette("Common_Patterns", "stpca")
#vignette("Creating_demo_data", "stpca")
#vignette("Reducing_Networks", "stpca")
## ====================================================================
##
## Creating demo data
##
## ====================================================================
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
data(demoNet)
plot(demoNet, VariableName	= "Sim_Values", 
     xlab="x-coordinate", ylab="y-coordinate",
     cex=1.5) #demoNet: 공간자료
#demoNet@obspoints@SSNPoints의 Sim.value가 출력
text(demoNet@obspoints@SSNPoints[[1]]@point.coords[,1], demoNet@obspoints@SSNPoints[[1]]@point.coords[,2], c(1:30))

#논문에 나오는 additive function value
demoNet@obspoints@SSNPoints[[1]]@point.data$addfunccol
round(demoNet@obspoints@SSNPoints[[1]]@point.data$shreve/33,5)==round(demoNet@obspoints@SSNPoints[[1]]@point.data$addfunccol,5)
#SSN 패키지에 additive.function이라는 함수도 있음
## ---- echo=FALSE---------------------------------------------------------
data(demoY)
##Time series plot (우리는 필요 없음)
#plot(1:25, demoY[,1], type="l", ylim=c(min(demoY), max(demoY)),
#     xlab="Time", ylab="Sim_Values", main="simulated time series")
#for(i in 1:30) {
#  lines(1:25, demoY[,i], col=i)
#} #demoY: 시공간자료


#data(demoYmiss) #demoY에 15%(112개)의 자료에 missing을 만들었다

##길이 112개의 shape
#shape_miho <- readOGR("~/Dropbox/R files/StreamFlow/River/Geum.shp")

#shape_miho: data (OBJECTID   RCH_ID    BRU_X    BRU_Y    BLL_X    BLL_Y   Shape_Leng)
#lines: list obj
#bbox
#shape_miho@proj4string

shape_demoNet <- readOGR("~/Dropbox/R files/StreamFlow/River/GeumMiho.shp")
shape_demoNet@data <-  demoNet@data
shape_demoNet@lines <- demoNet@lines
shape_demoNet@bbox <- demoNet@bbox
shape_demoNet@proj4string <-  demoNet@proj4string

#class(shape_miho)
#[1] "SpatialLinesDataFrame"
#attr(,"package")
#[1] "sp"

#shape_miho$RCH_ID

##현 시점의 문제점: get_binaryIDs_stream 함수를 고쳐야
##get_adjacency_stream 함수의 경우, 

#smnet 패키지의 get_binaryIDs 함수부터 먼저 살펴보자
#get_binaryIDs_stream: source에 만들어놓은 함수: 뒤의 get_adjacency_stream 실행 위함
binaryIDs_obj <- get_binaryIDs_stream(mouth_node=1 , shape_demoNet, RCH_ID=shape_demoNet$rid)
#수동으로 adjacency를 얻도록 만들어진 함수 binaryIDs_obj를 obj로 받기 때문에 위의 것만 해결되면 될 것으로 예상
adjacency <- get_adjacency_stream(binaryIDs_obj)
adjacency_old <- adjacency

########################################
#Data generaation
########################################
upper_seg <- which(colSums(adjacency_old$adjacency)==0) #이점들에서만 생성하면 됨 (총 55개)
upper_seg_list <- list()
upper_seg_label <- rep("A", 80)
upper_seg_list[[1]] <- c(71,77)
upper_seg_list[[2]] <- c(57,78,79)
upper_seg_list[[3]] <- c(27,55,75,76,80)
upper_seg_list[[4]] <- c(25,37,50,51,60,61,69,70)
upper_seg_list[[5]] <- c(52,53,66,73,74)
upper_seg_list[[6]] <- c(39,58)
upper_seg_list[[7]] <- c(29,40,67)
upper_seg_list[[8]] <- c(15)
upper_seg_list[[9]] <- c(13)
upper_seg_list[[10]] <- c(7)
upper_seg_list[[11]] <- c(4,5)

#upper_seg_label[c(15,71,77)] <- "B"
#upper_seg_label[c(13,27,55,57,75,76,78,79,80)] <- "C"
#upper_seg_label[c(22,37,50,51,60,61,69,70)] <- "D"
#upper_seg_label[c(52,53,66,73,74)] <- "E"
#upper_seg_label[c(22,39,40,58,67)] <- "F"
#upper_seg_label[c(4,5,7)] <- "G"

#label 다시 매기기
upper_seg_label[c(9,14,15,47,62,63,71,72,77)] <- "B"
upper_seg_label[c(22,24,25,36,37,42,43,44,45,50,51,60,61,68,69,70)] <- "C"
upper_seg_label[c(23,30,31,32,35,52,53,65,66,73,74)] <- "D"
upper_seg_label[c(16,18,19,38,39,58)] <- "E"
upper_seg_label[c(17,28,29,40,41,46,67)] <- "F"
upper_seg_label[c(1,2,3,4,5,6,7,8,11)] <- "G"


#interlude: plotting
#plot(demoNet)
#plot(demoNet@bbox, type="n", xlim=c(-5,2), ylim=c(0,14))
#for(i in 1:80){
#  coord.imsi <- colMeans(demoNet@lines[[i]]@Lines[[1]]@coords)
#  if(i %in% upper_seg){
#    text(coord.imsi[1], coord.imsi[2], i, col="red")
#  }else{
#    text(coord.imsi[1], coord.imsi[2], i)
#  }
#}

#compute_shreve: shreve 거리를 재도록 직접 제작한 함수
shreve_order <- compute_shreve(adjacency)

########################################
#start data analysis
########################################
data <- demoNet@obspoints@SSNPoints[[1]]@point.data$Sim_Values
eventplace2 <- cbind(demoNet@obspoints@SSNPoints[[1]]@point.data$locID,
                     demoNet@obspoints@SSNPoints[[1]]@point.coords[,1],
                     demoNet@obspoints@SSNPoints[[1]]@point.coords[,2])
colnames(eventplace2) <- c("ID", "x", "y")

example_network <- demoNet


#RID를 0부터 끝까지 정렬하는 업데이트
for(i in 1:length(example_network@lines)){
  example_network@lines[[i]]@ID <- as.character(i-1)
}

example_network@obspoints@SSNPoints[[1]]@point.data$rid
#example_network@obspoints@SSNPoints[[1]]@point.data$rid의 class가 factor이면 오류 발생
example_network@obspoints@SSNPoints[[1]]@point.data$rid <- as.numeric(as.character(example_network@obspoints@SSNPoints[[1]]@point.data$rid))
example_network@obspoints@SSNPoints[[1]]@point.data$rid

mf04 <- example_network
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

weight_vec_candidate2 <- compute_shreve_and_dist(adjacency_old, example_network,type="lengthprop")
weight_vec_candidate2 <- scales::rescale(log(sqrt(weight_vec_candidate2$distweight)), to=c(0.2, 1.5))
weight_vec_candidate <- weight_vec_candidate2
#weight_vec_candidate2 <- compute_shreve_and_dist(adjacency_old, example_network, scalevec=c(0.2, 1.5))
#weight_vec_candidate <- weight_vec_candidate2$distweight


#출력을 위해 직선을 매우 잘게 나누어야 함
TweedPredPoints_old <- data.frame(Latitude=c(), Longitude=c(), StreamUnit=c(), Weights=c())
for(jj in 1:length(example_network@lines)){
  #if(nrow(example_network@lines[[jj]]@Lines[[1]]@coords)>2){
  #  seg_cand <- example_network@lines[[jj]]@Lines[[1]]@coords[-nrow(example_network@lines[[jj]]@Lines[[1]]@coords),]
  #}else{
  #여기서는 직선이가 때문에 30개로 짜른다
  seg_cand <- cbind(seq(from=example_network@lines[[jj]]@Lines[[1]]@coords[1,1], to=example_network@lines[[jj]]@Lines[[1]]@coords[2,1], length.out = 30),
                    seq(from=example_network@lines[[jj]]@Lines[[1]]@coords[1,2], to=example_network@lines[[jj]]@Lines[[1]]@coords[2,2], length.out = 30))
  #}
  TweedPredPoints_old_imsi <- data.frame(Latitude=seg_cand[,2], Longitude=seg_cand[,1], 
                                         StreamUnit=rep(example_network@data$rid[jj], nrow(seg_cand)), 
                                         Weights=rep(weight_vec_candidate[jj], nrow(seg_cand)))
  TweedPredPoints_old <- rbind(TweedPredPoints_old, TweedPredPoints_old_imsi)
}
TweedPredPoints <- TweedPredPoints_old


#plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col =TweedPredPoints$StreamUnit, pch=16, cex=TweedPredPoints$Weights)
#points(TweedPredPoints$Longitude[1], TweedPredPoints$Latitude[2], pch=5, cex=2)


#출력: 나뭇가지 그림과 관찰측정소의 segment number
par(mar=c(1.1,1.1,1.1,1.1))
par(mfrow=c(1,1))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col = TweedPredPoints$StreamUnit, pch=19, cex=TweedPredPoints$Weights*2, xlab="", ylab="", xaxt="n", yaxt="n")
text(x=mf04@obspoints@SSNPoints[[1]]@point.coords[,1], y=mf04@obspoints@SSNPoints[[1]]@point.coords[,2], labels=mf04@obspoints@SSNPoints[[1]]@point.data$rid, col="brown", cex=1.3)
points(TweedPredPoints$Longitude[1], TweedPredPoints$Latitude[2], pch=5, cex=2)

#출력: 나뭇가지 그림과 관찰측정소의 observation numer
par(mar=c(1.1,1.1,1.1,1.1))
par(mfrow=c(1,1))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col = TweedPredPoints$StreamUnit, pch=19, cex=TweedPredPoints$Weights*2, xlab="", ylab="", xaxt="n", yaxt="n")
text(x=mf04@obspoints@SSNPoints[[1]]@point.coords[,1], y=mf04@obspoints@SSNPoints[[1]]@point.coords[,2], labels=mf04@obspoints@SSNPoints[[1]]@point.data$pid, col="brown", cex=1.3)
points(TweedPredPoints$Longitude[1], TweedPredPoints$Latitude[2], pch=5, cex=2)

#출력: 나뭇가지 그림과 stream segment number
par(mar=c(1.1,1.1,1.1,1.1))
par(mfrow=c(1,1))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col = TweedPredPoints$StreamUnit, pch=19, cex=TweedPredPoints$Weights*2, xlab="", ylab="", xaxt="n", yaxt="n")
for(ii in 1:length(mf04@lines)){
  text(x=mean(mf04@lines[[ii]]@Lines[[1]]@coords[,1]), y=mean(mf04@lines[[ii]]@Lines[[1]]@coords[,2]), labels=ii, col="brown", cex=1.3)
}



#plot(example_network@obspoints@SSNPoints[[1]]@point.coords, typ="n")
#text(x=example_network@obspoints@SSNPoints[[1]]@point.coords[,1],
#     y=example_network@obspoints@SSNPoints[[1]]@point.coords[,2], c(1:30))

#location_imsi
location_imsi <- c()
distupstream_imsi <- c()
for(i in 1:length(example_network@lines)){
  if(nrow(example_network@lines[[i]]@Lines[[1]]@coords)==2){
    location_imsi <- rbind(location_imsi, colMeans(example_network@lines[[i]]@Lines[[1]]@coords))
    #distupstream_imsi <- rbind(distupstream_imsi, riverdistance(startseg=i, endseg=68, startvert=1, endvert=20, rivers=MyRivernetwork))
  }else{
    index.imsi2 <- (nrow(example_network@lines[[i]]@Lines[[1]]@coords)/2)-1/2
    location_imsi <- rbind(location_imsi, colMeans(example_network@lines[[i]]@Lines[[1]]@coords[c(index.imsi2 , index.imsi2+1),]))
    #distupstream_imsi <- rbind(distupstream_imsi, riverdistance(startseg=i, endseg=68, startvert=index.imsi2, endvert=20, rivers=MyRivernetwork))
  }
}

demoNet2 <- demoNet
for(j in 1:length(demoNet2@lines)){
  location_sub <- demoNet2@lines[[j]]@Lines[[1]]@coords[2,]-demoNet2@lines[[j]]@Lines[[1]]@coords[1,]
  lines_new <- rbind(c( (demoNet2@lines[[j]]@Lines[[1]]@coords[1,]+ location_sub*c(1/4)) ), c( (demoNet2@lines[[j]]@Lines[[1]]@coords[1,]+ location_sub*c(1/3)) ), c( (demoNet2@lines[[j]]@Lines[[1]]@coords[1,]+ location_sub*c(2/4)) ), c( (demoNet2@lines[[j]]@Lines[[1]]@coords[1,]+ location_sub*c(2/3)) ), c( (demoNet2@lines[[j]]@Lines[[1]]@coords[1,]+ location_sub*c(3/4)) ) )
  
  demoNet2@lines[[j]]@Lines[[1]]@coords <- rbind(demoNet2@lines[[j]]@Lines[[1]]@coords[1,], lines_new, demoNet2@lines[[j]]@Lines[[1]]@coords[2,])
}

Sl <- SpatialLines(demoNet2@lines)
splndf <- SpatialLinesDataFrame(sl = Sl, data = as.data.frame(cbind(rep(0, 80), rep(1,80))), match.ID = FALSE)

proj4string(splndf) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

writeOGR(splndf, dsn="~/Dropbox/R files/StreamFlow/Simulation/" ,layer="Sl",driver="ESRI Shapefile", overwrite_layer = TRUE)
splndf_new <- readOGR("~/Dropbox/R files/StreamFlow/Simulation/Sl.shp")
plot(splndf_new)

proj4string(splndf_new) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

MyRivernetwork <- line2network(path="~/Dropbox/R files/StreamFlow/Simulation/Sl.shp", layer="Sl", reproject="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", tolerance = 0.001)
MyRivernetwork <- setmouth(seg=1, vert=7, rivers=MyRivernetwork)

detectroute(start = 10, end = 2, rivers=MyRivernetwork)
riverdistance(startseg = 1, startvert = 7, endseg = 2, endvert = 1, rivers=MyRivernetwork)


#plotting networks (나중에 손보기)
par(family = 'sans') 
par(mar=c(1.1,1.1,2.1,1.1))
par(mfrow=c(1,2))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, xlab="n", ylab="n", xaxt="n", yaxt="n", main="(c)")
points(location_imsi[-upper_seg,1], location_imsi[-upper_seg,2], pch=21, bg="orange", col="black")
points(location_imsi[upper_seg,1], location_imsi[upper_seg,2], pch=22, bg="purple", col="black")
library(mixtools)
for(iii in c(2,3,4,5,7)){
  ellipse(mu=c(mean(location_imsi[upper_seg_list[[iii]],1]), mean(location_imsi[upper_seg_list[[iii]],2])), sigma=var(location_imsi[upper_seg_list[[iii]],])/1.5, npoints=200, newplot=FALSE, col="red")
}
library(plotrix)
for(iii in c(1,6,8,9,10,11)){
  if(iii==1 | iii==6){
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.65, border="red" )
  }else if(iii==11){
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.45, border="red" )
  }else{
    draw.circle(x=mean(location_imsi[upper_seg_list[[iii]],1]), y=mean(location_imsi[upper_seg_list[[iii]],2]), radius=0.25, border="red" )
  }
}
legend("topleft", pch=c(15,16), col=c("purple", "orange"), c("Most upstream", "Non-most upstream"), bty="n")

unique_upper_seg_label <- unique(upper_seg_label)
rainbow_palette <- rainbow(length(unique_upper_seg_label))
plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col=rainbow_palette[match(upper_seg_label[TweedPredPoints$StreamUnit], unique_upper_seg_label)], pch=16, cex=TweedPredPoints$Weights, xlab="n", ylab="n", xaxt="n", yaxt="n", main="(d)")
points(location_imsi[-upper_seg,1], location_imsi[-upper_seg,2], pch=21, bg="orange", col="black")
points(location_imsi[upper_seg,1], location_imsi[upper_seg,2], pch=22, bg="purple", col="black")
legend("topleft", pch=c(15,16), col=c("purple", "orange"), c("Most upstream", "Non-most upstream"), bty="n")



########################################
#Start iteration
########################################
n.seg <- 80
#n.sub.length <- 80
n.sub.length <- 30
sd.val <- 2


#adjacency_old 작성
adjacency_old <- get_adjacency_stream(binaryIDs_obj)
adjacency <- as.data.frame(as.matrix(adjacency_old$adjacency))


realweights <- as.matrix(weight_vec_candidate, ncol=1)

n.iter <- 100
#n.iter <- 5
result_mat <- matrix(0, nrow=n.iter, ncol = 4)

result_list <- list()

for(ijk in 1:n.iter){
  ########여기서부터 실행
  #data generation
  print(paste("Iteration ", ijk))
  
  #generate binomial: number of the samples for each segment
  num.obs <- rep(1, 80)
  
  #location_binom: matrix
  location_binom <- cbind(demoNet@obspoints@SSNPoints[[1]]@point.coords,
                          as.numeric(paste(demoNet@obspoints@SSNPoints[[1]]@point.data$rid)),
                          rep(1, length(demoNet@obspoints@SSNPoints[[1]]@point.data$rid)),
                          demoNet@obspoints@SSNPoints[[1]]@point.data$upDist
                          ) #segid가 따로 없음(별로 안중요할듯)

  n.seg <- nrow(location_binom)
  colnames(location_binom) <- c("lon", "lat", "rid", "segid", "upDist")
  
  n <- 1; 
  #x <- sort((sample(1:(365*3), n)/(365*3)))
  x <- 1/5 #임의의 시간대
  #add noise
  noise_mat <- matrix(rnorm(n*n.seg, 0, sd.val),n,n.seg)

  #result_signal <- result_signal_partial + result_fct.1 + 5
  result_signal <- demoNet@obspoints@SSNPoints[[1]]@point.data$Sim_Values
  result_signal <- as.row(result_signal)
  colnames(result_signal) <- demoNet@obspoints@SSNPoints[[1]]@point.data$rid
  while(TRUE){
    result_data_noisy <- result_signal + noise_mat
    if(min(result_data_noisy)>0){
      break
    }else
      noise_mat <- matrix(rnorm(n*n.seg, 0, sd.val),n,n.seg)
  }
  
  #(여기서부터 편집해)
  #result_data_noisy2 <- matrix(result_data_noisy[order(as.numeric(colnames(result_data_noisy)))], nrow=1)
  #colnames(result_data_noisy2) <- colnames(result_data_noisy)[order(as.numeric(colnames(result_data_noisy)))]
  #TweedData generation (실제 자료의 duplication 반영)
  TweedData <- data.frame(date=as.Date("2012-12-31")+x*(365*3),
                          #location=1,
                          location=location_binom[1,3],
                          nitrate=result_data_noisy[1,1],
                          #nitrate=result_signal[,example_network@obspoints@SSNPoints[[1]]@point.data$rid[1]] + rnorm(1, mean=0, sd.val),
                          #Long=rep(location_imsi[1,1],n),
                          Long=location_binom[1,1],
                          #Lat=rep(location_imsi[1,2],n)
                          Lat=location_binom[1,2]
  )
  for(i in 2:nrow(location_binom)){
    dataframe.new <- data.frame(date=as.Date("2012-12-31")+x*(365*3),
                                location=location_binom[i,3],
                                nitrate=result_data_noisy[1,i],
                                Long=location_binom[i,1],
                                Lat=location_binom[i,2]
    )
    TweedData <- rbind(TweedData, dataframe.new)
  }
  TweedData <- TweedData[order(TweedData$location),]
  #(Oct 31, 2020): TweedData 추가작업
  TweedData2 <- TweedData
  TweedData2[1,c(3:5)] <- colMeans(TweedData2[c(1:2),c(3:5)])
  TweedData2[11,c(3:5)] <- colMeans(TweedData2[c(11:12),c(3:5)])
  TweedData2[14,c(3:5)] <- colMeans(TweedData2[c(14:15),c(3:5)])
  TweedData2[16,c(3:5)] <- colMeans(TweedData2[c(16:17),c(3:5)])
  TweedData2[19,c(3:5)] <- colMeans(TweedData2[c(19:20),c(3:5)])
  TweedData2[25,c(3:5)] <- colMeans(TweedData2[c(25:26),c(3:5)])
  TweedData2 <- TweedData2[-c(2,12,15,17,20,26),]
  
 index_choose <- sort(as.numeric(paste(sort(unique(demoNet@obspoints@SSNPoints[[1]]@point.data$rid)))))
 index_choose2 <- as.numeric(paste(demoNet@obspoints@SSNPoints[[1]]@point.data$rid))
 
  #sampling을 할 경우 
  
  #(1-5) create a new streamflow network
  example_network2 <- example_network

  data <- as.vector(result_data_noisy)
  names(data) <- index_choose2
  #data <- as.vector(result_data_noisy[order(index_choose2)])
  #names(data) <- index_choose2[order(index_choose2)]
  #TweedData <- TweedData[index_choose,]
  
  
  #(1) smnet_ST function
  penalties_default <- c(50, 50, c(25, 5), 50, c(25, 5), c(25, 25))
  division_factor <- c(30, 25, 20, 15, 10, 5, 1, 0.5, 0.25, 0.125)
  AICc_vec <- rep(0, length(division_factor))
  #AICc_vec_log <- rep(0, length(division_factor))
  for(jjj in 1:length(division_factor)){
    result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/division_factor[jjj], plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=FALSE, model.type = c("c", "m"))
    AICc_vec[jjj] <- result_a$AICc
    #result_b <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties=penalties_default/division_factor[jjj], plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=FALSE, model.type = c("c", "m"))
    #AICc_vec_log[jjj] <- result_b$AICc
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
  
  
  #AICc_vec_min_index_log <- which.min(AICc_vec_log)
  #
  #if(AICc_vec_min_index_log==1){
  #  opt_result_log <- optim(par=penalties_default/division_factor[AICc_vec_min_index_log], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index_log)], upper=penalties_default/division_factor[(AICc_vec_min_index_log+1)], log.y=TRUE, model.type = c("c", "m") )
  #}else if(AICc_vec_min_index_log==length(AICc_vec_log)){
  #  opt_result_log <- optim(par=penalties_default/division_factor[AICc_vec_min_index_log], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index_log)-1], upper=penalties_default/division_factor[(AICc_vec_min_index_log)], log.y=TRUE, model.type = c("c", "m") )
  #}else{
  #  opt_result_log <- optim(par=penalties_default/division_factor[AICc_vec_min_index_log], fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/division_factor[(AICc_vec_min_index_log-1)], upper=penalties_default/division_factor[(AICc_vec_min_index_log+1)], log.y=TRUE, model.type = c("c", "m") )
  #}
  #opt_val_log <- opt_result_log$par
  
  #opt_result <- optim(par=penalties_default, fn=smnet_ST, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints, method="L-BFGS-B", lower=penalties_default/10, upper=penalties_default*2, log.y=FALSE,model.type = c("c", "m"))
  result_a <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties= opt_result$par, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=FALSE,model.type = c("c", "m") )
  unique(result_a$fit)
  #result_a$fit <- result_a$fit
  index_vec_imsi <- cbind(1:30, order(as.numeric(colnames(result_signal))))
  index_vec_imsi <- index_vec_imsi[order(index_vec_imsi[,2]),]
  result_a$fit <- result_a$fit[index_vec_imsi[,1]]
  colnames(result_signal)
  TweedData$location[index_vec_imsi[,1]]
  if(sum(as.numeric(colnames(result_signal))!=TweedData$location[index_vec_imsi[,1]])!=0){
    stop("matching error")
  }
  sapply(split(result_a$fit, f=result_a$TweedData$location), function(x) median(x))
  
  #result_b <- smnet_ST(realweights, adjacency, TweedData, TweedPredPoints, penalties= opt_result_log$par, plot.fig=FALSE, station=NULL, use.optim=FALSE, log.y=TRUE,model.type = c("c", "m") )
  
  #start_time <- Sys.time()
  
  result_denoise <- denoise_S(data, example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", returnall = FALSE)
  result_denoise2 <- denoise_S(data, example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "hard", returnall = FALSE)
  
  #end_time <- Sys.time()
  
  result_forward <- fwtnp_stream_S(data, example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
  
  cores=detectCores()
  cl <- makeCluster(cores[1]-1,  setup_strategy = "sequential") #not to overload your computer
  
  #(Oct 24, 2020): index_sub_data 수정해야
  index_sub_data <- data.frame(stations=c(1:80), groups=as.factor(upper_seg_label))
  registerDoParallel(cl)
  
  #denoise_Stream_S_perm(data, example_network2, endpt=1, per=NULL, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = FALSE, plot.fig = FALSE, plot.individual = FALSE, pollutant=NULL, polluyear=NULL, plot.thesis=FALSE)
  #result_denoise3 <- nlt_Stream_S(data, example_network2, J=10, endpt=1, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = TRUE)$aveghat
  result_denoise3 <- nlt_Stream_S(data, example_network2, J=10, endpt=1, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, pred = streamPred_S, neigh = 1, int = TRUE, clo = FALSE, keep = 2, rule = "median", sd.scale=1, returnall = TRUE, ga=TRUE, index_sub_data=index_sub_data, oremovelist=result_forward$removelist)$aveghat
  
  stopCluster(cl)
  
  #result_forward <- fwtnp_stream_S(data, example_network2, adjacency=adjacency_old, realweights, TweedData, TweedPredPoints, nkeep = 2, intercept = TRUE, initboundhandl = "reflect",  neighbours = 1, closest = FALSE, LocalPred = streamPred_S, do.W = FALSE, varonly = FALSE)
  #lr <- result_forward$lengthsremove
  #rem <- result_forward$removelist
  #newcoeff <- result_forward$coeff
  #keep=2  
  #fhat <- invtnp_stream_S(X=as.vector(as.numeric(names(data))), newcoeff, result_forward$lengths, lr, result_forward$pointsin,  rem, result_forward$neighbrs, result_forward$schemehist, result_forward$interhist, length(result_forward$x) - keep, int, neighbours=1, clo, LocalPred=streamPred_S,  data, example_network, adjacency=adjacency_old, realweights)
  n.observed <- length(as.numeric(unique(names(data))))
  #result_mat[ijk,1] <-  sqrt(sum((colMeans(result_signal)[which(!duplicated(names(data)))]-sapply(split(TweedData$nitrate, f=TweedData$location), function(x) mean(x)))^2)/n.observed)
  
  #result_mat[ijk,1] <-  sqrt(sum((colMeans(result_signal)[index_choose]-result_a$fit)^2)/n.observed)
  #print(paste("O'Donnell ", round(result_mat[ijk,1],4)))
  #result_mat[ijk,2] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise)^2)/n.observed)
  #print(paste("Proposed ", round(result_mat[ijk,2],4)))
  #result_mat[ijk,3] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise2)^2)/n.observed)
  #print(paste("Proposed(hard) ", round(result_mat[ijk,3],4)))
  
  par(family = 'sans') 
  par(mar=c(1.1,1.1,2.1,1.1))
  par(mfrow=c(2,3))
  
  zlims <-range(c(data, TweedData$nitrate, result_a$fit, result_denoise, result_denoise2, result_denoise3, colMeans(result_signal)))+c(-0.5,0.5)
  result_signal2 <- rep(NA, length(sort(unique(TweedPredPoints$StreamUnit))))
  for(iii in 1:length(sort(unique(TweedPredPoints$StreamUnit)))){
    signal_cand <- which(as.numeric(colnames(result_signal))==iii)
    if(length(signal_cand)!=0){
      result_signal2[iii] <- mean(result_signal[signal_cand])
    }
  }
  result_signal2_name <- which(!is.na(result_signal2))
  result_signal2 <- as.vector(result_signal2[complete.cases(result_signal2)])
  result_signal2 <- matrix(result_signal2, nrow=1)
  colnames(result_signal2) <- result_signal2_name
  
  scatter_fill(TweedPredPoints$Longitude[which(TweedPredPoints$StreamUnit %in% as.numeric(colnames(result_signal2)))], TweedPredPoints$Latitude[which(TweedPredPoints$StreamUnit %in% as.numeric(colnames(result_signal2)))], rep(colMeans(result_signal2), each=30), pch=16, cex=TweedPredPoints$Weights[which(TweedPredPoints$StreamUnit%in%index_choose)], xlab="", ylab="", xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude), zlim=zlims, main="(a) True", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
  points(TweedPredPoints$Longitude[-which(TweedPredPoints$StreamUnit%in%index_choose)], TweedPredPoints$Latitude[-which(TweedPredPoints$StreamUnit%in%index_choose)] ,col="gray", pch=16, cex=TweedPredPoints$Weights[-which(TweedPredPoints$StreamUnit%in%index_choose)])
  #points(location_binom[index_choose,1], location_imsi[index_choose,2], pch=0, cex=1)
  
  scatter_fill(TweedPredPoints$Longitude[which(TweedPredPoints$StreamUnit%in%index_choose)], TweedPredPoints$Latitude[which(TweedPredPoints$StreamUnit%in%index_choose)] , rep(TweedData2$nitrate, each=30), pch=16, cex=TweedPredPoints$Weights[which(TweedPredPoints$StreamUnit%in%index_choose)], xlab="", ylab="", xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude), zlim=zlims, main="(b) Observed", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
  #scatter_fill(TweedPredPoints$Longitude[which(TweedPredPoints$StreamUnit%in%index_choose)], TweedPredPoints$Latitude[which(TweedPredPoints$StreamUnit%in%index_choose)] ,TweedData$nitrate[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights[which(TweedPredPoints$StreamUnit%in%index_choose)], xlab="", ylab="", xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude), zlim=zlims, main="(b) Observed", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
  points(TweedPredPoints$Longitude[-which(TweedPredPoints$StreamUnit%in%index_choose)], TweedPredPoints$Latitude[-which(TweedPredPoints$StreamUnit%in%index_choose)] ,col="gray", pch=16, cex=TweedPredPoints$Weights[-which(TweedPredPoints$StreamUnit%in%index_choose)])
  #points(location_binom[index_choose,1], location_imsi[index_choose,2], pch=0, cex=1)
  
  
  #old version
  #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(a) True", xlab="n", ylab="n", xaxt="n", yaxt="n")
  #quilt.plot(location_binom[,1], location_binom[,2], colMeans(result_signal), cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
  
  #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(b) Observed", xlab="n", ylab="n", xaxt="n", yaxt="n")
  #quilt.plot(location_binom[index_choose,1], location_imsi[index_choose,2], TweedData$nitrate, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
  #points(location_imsi[index_sub, 1], location_imsi[index_sub, 2], pch=5, cex=2)
  
  result_new_mat <- matrix(0, nrow=4, ncol=5)
  rownames(result_new_mat) <- c( "ODonnell", "S-Lifting(M)",  "S-Lifting(H)", "S-Lifting(N)")
  colnames(result_new_mat) <- c("Corr", "RMSE(whole)","WRMSE(whole)", "RMSE(actual)", "WRMSE(actual)")
  
  if(length(index_choose)==80){
    #예측 단계 불필요
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_a$fit[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(c) O'Donnell", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(d) Proposed (M)", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(e) Proposed (H)", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise3[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(f) Proposed (N)", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(c) O'Donnell", xlab="n", ylab="n", xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_a$fit, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(d) Proposed (M)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_denoise, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(e) Proposed (H)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_denoise2, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(f) Proposed (N)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_denoise3, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    
    result_mat[ijk,1] <-  sqrt(sum((colMeans(result_signal)[index_choose]-result_a$fit)^2)/n.observed)
    print(paste("O'Donnell ", round(result_mat[ijk,1],4)))
    result_mat[ijk,2] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise)^2)/n.observed)
    print(paste("Proposed ", round(result_mat[ijk,2],4)))
    result_mat[ijk,3] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise2)^2)/n.observed)
    print(paste("Proposed(hard) ", round(result_mat[ijk,3],4)))
    result_mat[ijk,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise3)^2)/n.observed)
    print(paste("Proposed(nlt, median) ", round(result_mat[ijk,4],4)))
    #result_mat[ijk,4] <-  sqrt(sum((colMeans(result_signal)[index_choose]-exp(result_b$fit))^2)/n.observed)
    #print(paste("O'Donnell(log) ", round(result_mat[ijk,4],4)))
    
    result_new_mat[1,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_a$fit)) #O'Donnell
    result_new_mat[2,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_denoise)) #S-Lifting(M)
    result_new_mat[3,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_denoise2)) #S-Lifting(H)
    result_new_mat[4,1] <- cor(colMeans(result_signal)[index_choose], as.numeric(result_denoise3)) #S-Lifting(N)
    
    result_new_mat[1,2] = result_new_mat[1,4] <- result_mat[ijk,1]
    result_new_mat[2,2] = result_new_mat[2,4] <- result_mat[ijk,2]
    result_new_mat[3,2] = result_new_mat[3,4] <- result_mat[ijk,3]
    result_new_mat[4,2] = result_new_mat[4,4] <- result_mat[ijk,4]
    
    result_new_mat[1,3] = result_new_mat[1,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_a$fit))^2)/n.observed)
    result_new_mat[2,3] = result_new_mat[2,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise))^2)/n.observed)
    result_new_mat[3,3] = result_new_mat[3,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise2))^2)/n.observed)
    result_new_mat[4,3] = result_new_mat[4,5] <- sqrt(sum(realweights*(colMeans(result_signal)[index_choose]-as.numeric(result_denoise3))^2)/n.observed)
    
  }else{
    #예측 단계 필요
    time_init <- Sys.time()
    initS_obj = initS_stream(X=as.row(as.numeric(names(data))), data=as.vector(data), example_network2, adjacency=adjacency_old, realweights, pointsin= matrix(1:length(data), 1, length(data)))
    time_end <- Sys.time()
    time_end - time_init
    result_a_new <- initS_obj$weight_matrix%*%as.column(result_a$fit)
    #result_b_new <- initS_obj$weight_matrix%*%as.column(result_b$fit)
    result_denoise_new <- initS_obj$weight_matrix%*%as.column(result_denoise)
    result_denoise_new2 <- initS_obj$weight_matrix%*%as.column(result_denoise2)
    result_denoise_new3 <- initS_obj$weight_matrix%*%as.column(result_denoise3)
    
    rbind(t(result_a$fit), result_denoise, result_denoise2, result_denoise3)
    
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_a_new[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(c) O'Donnell", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise_new[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(d) Proposed (M)", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise_new2[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(e) Proposed (H)", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,result_denoise_new3[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xaxt="n", yaxt="n", zlim=zlims, main="(f) Proposed (N)", cex.main=1.5, smallplot=c(0.2,0.25,0.1,0.3), y.axes = FALSE, y.axes.label=FALSE)
    
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(c) O'Donnell", xlab="n", ylab="n", xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_a_new, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(d) Proposed (M)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_denoise_new, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(e) Proposed (H)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_denoise_new2, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    #
    #plot(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,col="gray", pch=16, cex=TweedPredPoints$Weights, main="(f) Proposed (N)", xlab="n", ylab="n",  xaxt="n", yaxt="n")
    #quilt.plot(location_binom[,1], location_binom[,2], result_denoise_new3, cex=2, add=T, zlim=range(c(data, colMeans(result_signal)))+c(-0.5,0.5), xlab="x", ylab="y")
    
    result_mat[ijk,1] <-  sqrt(sum((colMeans(result_signal)-result_a_new[as.numeric(colnames(result_signal))])^2)/80)
    print(paste("O'Donnell ", round(result_mat[ijk,1],4)))
    result_mat[ijk,2] <- sqrt(sum((colMeans(result_signal)-result_denoise_new[as.numeric(colnames(result_signal))])^2)/80)
    print(paste("Proposed ", round(result_mat[ijk,2],4)))
    result_mat[ijk,3] <- sqrt(sum((colMeans(result_signal)-result_denoise_new2[as.numeric(colnames(result_signal))])^2)/80)
    print(paste("Proposed(hard) ", round(result_mat[ijk,3],4)))
    result_mat[ijk,4] <- sqrt(sum((colMeans(result_signal)-result_denoise_new3[as.numeric(colnames(result_signal))])^2)/80)
    print(paste("Proposed(nlt, median) ", round(result_mat[ijk,4],4)))
    #result_mat[ijk,4] <-  sqrt(sum((colMeans(result_signal)-exp(result_b_new))^2)/80)
    #print(paste("O'Donnell(log) ", round(result_mat[ijk,4],4)))
    
    result_new_mat[1,1] <- cor(colMeans(result_signal), result_a_new[as.numeric(colnames(result_signal))]) #O'Donnell
    result_new_mat[2,1] <- cor(colMeans(result_signal), result_denoise_new[as.numeric(colnames(result_signal))]) #S-Lifting(M)
    result_new_mat[3,1] <- cor(colMeans(result_signal), result_denoise_new2[as.numeric(colnames(result_signal))]) #S-Lifting(H)
    result_new_mat[4,1] <- cor(colMeans(result_signal), result_denoise_new3[as.numeric(colnames(result_signal))]) #S-Lifting(N)
    
    result_new_mat[1,2] <- result_mat[ijk,1]
    result_new_mat[2,2] <- result_mat[ijk,2]
    result_new_mat[3,2] <- result_mat[ijk,3]
    result_new_mat[4,2] <- result_mat[ijk,4]
    
    result_new_mat[1,3] <- sqrt(sum(realweights[as.numeric(colnames(result_signal))]*(colMeans(result_signal)-result_a_new[as.numeric(colnames(result_signal))])^2)/80)
    result_new_mat[2,3] <- sqrt(sum(realweights[as.numeric(colnames(result_signal))]*(colMeans(result_signal)-result_denoise_new[as.numeric(colnames(result_signal))])^2)/80)
    result_new_mat[3,3] <- sqrt(sum(realweights[as.numeric(colnames(result_signal))]*(colMeans(result_signal)-result_denoise_new2[as.numeric(colnames(result_signal))])^2)/80)
    result_new_mat[4,3] <- sqrt(sum(realweights[as.numeric(colnames(result_signal))]*(colMeans(result_signal)-result_denoise_new3[as.numeric(colnames(result_signal))])^2)/80)
    
    result_new_mat[1,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_a$fit)^2)/length(index_choose))
    result_new_mat[2,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise)^2)/length(index_choose))
    result_new_mat[3,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise2)^2)/length(index_choose))
    result_new_mat[4,4] <- sqrt(sum((colMeans(result_signal)[index_choose]-result_denoise3)^2)/length(index_choose))
    
    result_new_mat[1,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-result_a$fit)^2)/length(index_choose))
    result_new_mat[2,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-result_denoise)^2)/length(index_choose))
    result_new_mat[3,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-result_denoise2)^2)/length(index_choose))
    result_new_mat[4,5] <- sqrt(sum(realweights[index_choose]*(colMeans(result_signal)[index_choose]-result_denoise3)^2)/length(index_choose))
  }
  result_list[[ijk]] <- result_new_mat
}
#saveRDS(result_mat, "StreamSTPCA80(sd3).RDS")
#saveRDS(result_mat, "StreamSTPCA40(sd05).RDS")

#saveRDS(result_mat, "StreamSTPCA40(sd1)nlt.RDS")
saveRDS(result_list, "RealListStreamSTPCA40(sd2)nlt.RDS")

#evaluation
aaaa <- readRDS("~/Dropbox/Github/paper-StreamflowLifting/result_RDS/RealListStreamSTPCA40(sd1)nlt.RDS")
#aaaa <- readRDS("~/Dropbox/Github/paper-StreamflowLifting/result_RDS/RealListStreamSTPCA40(sd15)nlt.RDS")
#aaaa <- readRDS("~/Dropbox/Github/paper-StreamflowLifting/result_RDS/RealListStreamSTPCA40(sd2)nlt.RDS")
bbbb <- array(as.numeric(unlist(aaaa)), dim=c(4,5,100))
mean(bbbb[1,2,]*(sqrt(80/30)))
mean(bbbb[2,2,]*(sqrt(80/30)))
mean(bbbb[3,2,]*(sqrt(80/30)))
mean(bbbb[4,2,]*(sqrt(80/30)))
sd(bbbb[1,2,]*(sqrt(80/30)))
sd(bbbb[2,2,]*(sqrt(80/30)))
sd(bbbb[3,2,]*(sqrt(80/30)))
sd(bbbb[4,2,]*(sqrt(80/30)))





