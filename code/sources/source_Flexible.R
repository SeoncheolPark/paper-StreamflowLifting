source('~/Dropbox/R files/StreamFlow/sources/source.R', chdir = TRUE)
source('~/Dropbox/R files/StreamFlow/sources/source_ST.R', chdir = TRUE)

setwd("~/Dropbox/(논문들)/(논문)River Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/")
library(mgcv)
library(spam)
library(RgoogleMaps)
library(ggmap)
library(lubridate)
# library(Matrix)
library(scales)

source('~/Dropbox/(논문들)/(논문)River Network/Extremes on River Network/Flexible regression models over river networks(2014)_code/paper_functions.r', chdir = TRUE)

smnet_ST <- function(realweights, adjacency, TweedData, TweedPredPoints, penalties=c(50, 50, c(25, 5), 50, c(25, 5), c(25, 25)), plot.fig=FALSE, station=NULL, use.optim=TRUE, log.y=FALSE, model.type=NULL){
  #plot.fig: whether to decide plot the figure
  #station.num: when plot.fig is TRUE, we have to choose which one to show
  
  if(log.y==FALSE){
    TweedData$nitrate <- exp(TweedData$nitrate)
  }
  
  # Rgooglemaps plot of the data for February 2004
  if(plot.fig==TRUE){
    # requires an internet connection to retrieve googlemaps tile
    #snap.shot<-c(2004, 02)
    snap.shot <- c(2017, 02)
    snap.shot <- c(year(as.character(TweedData[nrow(TweedData),1])), month(as.character(TweedData[nrow(TweedData),1])))
    x<-TweedData[c("Long", "Lat")]
    dates = as.character(TweedData[,1])
    all.sites<-aggregate(cbind(log(TweedData$nitrate), x[,1], x[,2]), by = list(TweedData$location), mean)
    ind<-which(as.numeric(substr(dates, 1, 4)) == snap.shot[1] & as.numeric(substr(dates, 6, 7)) ==  snap.shot[2])
    station.obs<-aggregate(cbind(log(TweedData$nitrate[ind]), x[ind,1], x[ind,2]), by = list(TweedData$location[ind]), mean) #log (TON)
    colnames(station.obs)<-c("location", "nitrate", "Long", "Lat")
    #n.cols<-n.grid<-121
    n.cols<-n.grid<-131
    brks<-range(log(TweedData$nitrate))
    brks<-seq(brks[1] - 0.1, brks[2] + 0.1, length = n.grid + 1)
    palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black"))
    main.cols<-palette(n.cols)
    col.nums.stats<-cut(station.obs$nitrate, breaks = brks, labels = FALSE)
    #MyMap<-GetMap(center = c(lat = mean(TweedPredPoints$Latitude), lon = mean(TweedPredPoints$Longitude)),  size = c(640, 640),zoom=9,maptype = "terrain")
    #PlotOnStaticMap(MyMap,lat = TweedPredPoints$Latitude, lon = TweedPredPoints$Longitude, pch = 21, add = F,  cex = 0.1, col = "blue")              
    #PlotOnStaticMap(MyMap, lat = all.sites[,4], lon = all.sites[,3],  pch = 21, add = TRUE,cex=1.7, col=1) 
    #PlotOnStaticMap(MyMap,lat = station.obs[,4], lon = station.obs[,3], pch = 21, add = TRUE,cex=1.7, bg=main.cols[col.nums.stats], col=1) 
    #par(mar = c(1.2, 8, 2.7, 8))
    #colour.key(main.cols, brks = range(log(TweedData$nitrate)))
    
    place.for.plot <- unique(TweedData[,c(4,5)])
    get_map(location=c(lat = sum(range(TweedPredPoints$Latitude))/2, lon = sum(range(TweedPredPoints$Longitude))/2 ), zoom=10, maptype="terrain", source="google") %>%
      ggmap(extent = "device") + geom_point(data=TweedPredPoints, aes(x=Longitude, y=Latitude)) + geom_point(data=place.for.plot, aes(x=Long, y=Lat), color="red", size=3)
    #my addition: 장소 중 한 곳을 더 출력해본다 
    #par(mar=c(5.1,4.1,4.1,2.1))
    #plot(which(!is.na(data$normaldata_TN[,8])),as.numeric(log(data$normaldata_TN[which(!is.na(data$normaldata_TN[,8])),8])), xlab="Time", ylab="Log(TON)", main="Miho-cheon5A")
  }
  
  
  
  
  # Fit a model to the river network
  #{
    # set the penalty parameters
    if(is.null(model.type)){
      #model.type = c("c", "m", "s", "si") 
      #model.type = c("c", "m", "s", "t") 
      model.type = c("c", "m", "s", "si", "t", "ti", "ts") #c: constant, m:network(spatial), s:seosonal, si:seosonal*spatial, t:temporal, ti:temporal*spatial, ts:temporal*seasonal
    }
    if(is.null(penalties)){
      penalties <- c(50, 50, c(25, 5), 50, c(25, 5), c(25, 25))  #optimal penalty selection (각 penalty의 의미 또한 조사)
    }
    #model type, penalties adjustment(MY ADDITION)
    #match.result <- match(c("c", "m", "s", "si", "t", "ti", "ts"), model.type)
    #penalties_new <- c()
    #if("c"%in%model.type){
    #  penalties_new <- c(penalties_new, penalties[1])
    #}
    
    response <- TweedData$nitrate
    response.locs <- TweedData$location
    knts<-c(10,10) #아마 the number of knots인듯(순서대로 s모형의 knot 갯수, t모형의 knot 갯수)
    dates = as.character(TweedData[,1]) #이렇게 코딩한 이유는 TweedData[,1]에 date를 넣었기 때문이다
    if(is.null(station)){
      station=38 #station number? 이것은 unique(TweedData$location)에서 선택
      #plotting에서만 쓰이는 듯?
    }
    s
    n.segments<-nrow(adjacency) #강 물줄기 갯수(113개)
    p.dims<-n.segments #p.dims: 강 물줄기 갯수로 정의
    dims<-get.dimension(model.type = model.type, n.segments = n.segments, knts = knts) #[1]    1  113   10 1130   10 1130  100 (c, m, s, si, t, ti, ts)
    c.dims<-cumsum(dims) #dims를 차례대로 더한 인덱스
    n.par<-sum(dims) #총 dims (2494)
    X.list<-vector("list")
    
    if("c"%in%model.type){
      # construct intercept column of model matrix ("c")
      X.list<-c(X.list, list(ones = as.spam(matrix(1, nrow = length(response),  ncol = 1)))) #length(y)*1 1로 구성된 벡터 반환
      model.mat<-X.list$ones
    }
    
    # construct network model matrix and network penalty ("m")
    n<-2 #n: 여기서 처음 등장한다
    p.n<-0
    if("m"%in%model.type){
      # construct network model matrix and network penalty ("m")
      #n<-2 #n: 여기서 처음 등장한다
      #p.n<-0
      X<-matrix(0, nrow = length(response), ncol = n.segments) #(데이터갯수)*(segment개수) 0행렬 생성
      for(i in 1:length(response)){X[i, response.locs[i]] <- 1} #(ex) 1번데이터는 38번 seg에서 관찰되었으므로 X[1,38] <- 1로 둔다
      X.list<-c(X.list, list(spatial = as.spam(X))) #X.list는 ones, spatial로 구성된다
      model.mat<-cbind.spam(model.mat, X.list$spatial) #4372*114 matrix로 확장된다 (ex. 처음 1~70개 자료가 38번 segment에서 나왔으면 1~70행까지 38열은 1, 나머지는 0)
      Z.list<-vector("list")
      rm(list = c("X")) #rm:remove
      D1<-matrix(0, n.segments, n.segments) #113*113 matrix
      D2<-matrix(0, n.segments, n.segments) #113*113 matrix
      # calculate difference matrix D on all pairs of adjacenct stream segments
      for(i in 1:n.segments){
        a <- which(adjacency[,i] == 1) #a가 i의 상류지점에 해당된다
        if(length(a) > 0)
        {
          D1[i, i]<--realweights[a[1],1] #diagonal element에 대해서는 음수를 취함
          D1[i, a[1]]<-realweights[a[1],1] #D1: 113*113 matrix 예를 들어서 1번 segment의 상류지점이 3,6이다 그러면 D1[1,3]=realweights[3], D1[1,1]=-realweights[3]
        }	
        if(length(a) > 1) #D2는 i의 상류지점 a가 2군데일 때만 작동
        {
          D2[i, i]<--realweights[a[2],1] #diagonal element에 대해서는 음수를 취함
          D2[i, a[2]]<-realweights[a[2],1] #D2: 113*113 matrix 예를 들어서 1번 segment의 상류지점이 3,6이다 그러면 D2[1,6]=realweights[6], D2[1,1]=-realweights[6]
        }
      }
      D<-rbind(D1, D2) #((113*2)*113) matrix
      rm.d<-which(rowSums(abs(D))  == 0) #일단 which(apply(adjacency, 2, sum)==0)인 포인트가 들어감(즉 상류가 존재하지 않는 곳, 최상류지점)
      #D2에서는 which(apply(adjacency, 2, sum)<=1)+113인 지점이 반환, 즉 which(apply(adjacency, 2, sum)==c(0,1)) 상류가 없거나 1개인 지점
      P_spat<-as.spam(D[-rm.d,]) #전체 (226*113) matrix들 중 rm.d에서 걸러진 114개 지점 제외한 나머지 112*113 matrix 반환
      PP1<-get.PTP(P = P_spat, dim.list = dims, i = match("m", model.type)) #dims: [1]    1  113   10 1130   10 1130  100 #2494*2494 행렬
      #PP1<-get.PTP(P = P_spat, dim.list = dims, i = 2) #dims: [1]    1  113   10 1130   10 1130  100 #2494*2494 행렬 (오리지널)
      #(i=2인 이유는 아마도 i: spatial이기 때문으로 분석된다)
      
      if(length(dims) == 1) PP1<-crossprod(P_spat) 
      P<-PP1*penalties[1] #penalties[1] 즉 50
      PTP<-t(P) %*% P #PTP: 2494*2494 행렬
    }else{
      PTP <- as.spam(matrix(0, nrow=sum(dims), ncol=sum(dims))) #나의 추가(m이 없을 경우)
      #(CAUTION: PTP spam matrix로 만들지 않을 경우 뒤에서 문제가 생긴다)
    }
    
    # construct seasonal component and seasonal penalty ("s")
    n<-n+1 #n=3이 됨
    #임시 편집
    dates2 <- ymd(dates); year(dates2) <- 2013 #(추가)
    decimal.day<-yday(dates2)/365 #(0,1)까지로 날짜 변환 #(추가)
    if("s"%in%model.type){
      # construct seasonal component and seasonal penalty ("s")
      #n<-n+1 #n=3이 됨
      #decimal.day<-yday(dates)/365 #(0,1)까지로 날짜 변환
      Bas<-cSplineDes(decimal.day, knots = (0:knts[1])/knts[1], ord=4) #윤달 그냥 놔두면 오류 발생 ##여기서 오류발생 #cSplineDes: Uses splineDesign to set up the model matrix for a cyclic B-spline basis. #seasonal이 10개이므로 knot 10
      #dim(Bas): 4372 10
      X.list<-c(X.list, list(seasonal = as.spam(Bas))) #ones, spatial 있는 list에 seasonal 추가함
      model.mat<-cbind.spam(model.mat, X.list$seasonal) #4372*124 matrix 
      # circular penalty matrix
      S<-diff(diag(dim(Bas)[2]), diff = 2) #8*10 matrix (1 -2 1)이 반복되는 구조 (diff가 2이므로 행이 8이 된 듯)
      last.row<-S[nrow(S), ] #last.row: S의 마지막 row 저장
      row1<-c(last.row[ncol(S)], last.row[-ncol(S)]) #row1: last.row에서 마지막 열을 맨 앞으로 붙인 벡터
      S<-as.matrix(rbind(S, row1)) #마지막 행만 바꾸어 다시 저장하였다
      last.row<-S[dim(S)[1],]
      row1<-c(last.row[dim(S)[2]], last.row[-dim(S)[2]]) #같은 작업을 1번 더 함
      P_seasonal<-as.spam(as.matrix(rbind(S, row1))) #(이해함): seasonality니까 diff2인 것을 반복시키게 하기 위함 (10*10 행렬)
      P.seasonal<-get.PTP(P = P_seasonal, dim.list = dims, i = match("s", model.type)) #P.seasonal: 2494*2494 행렬 (penalty matrix의 dimension을 모두 더하면 2494)
      #P.seasonal<-get.PTP(P = P_seasonal, dim.list = dims, i = n) #P.seasonal: 2494*2494 행렬 (penalty matrix의 dimension을 모두 더하면 2494) (오리지널)
      PTP<-PTP+P.seasonal*penalties[n-1] #penalties[2] 즉 50
    }
    
    # create Seasonal-network interaction component ("si")
    n<-n+1 #n <- 4
    if("si"%in%model.type){
      # create Seasonal-network interaction component ("si")
      #n<-n+1 #n <- 4
      # create Seasonal interaction design 
      seas_spat<-box3.prod(X.list$spatial, X.list$seasonal) #4372*1130 matrix) #box3.prod: 만들어놓은 함수 #X.list$spatial은 4372*113 matrix, #X.list$seasonal은 4372*10 matrix
      X.list<-c(X.list, list(seas_spat)) #X.list에 Seasonal-network interaction 추가
      model.mat<-cbind.spam(model.mat, seas_spat) #4372*1254 matrix
      P_si1<-kronecker.spam(P_spat, diag.spam(ncol(P_seasonal))) #kronecker.spam: Computes the generalised kronecker product of two arrays, X and Y. #1120*1130 matrix (앞의 행렬 112*113, 뒤의 행렬 10*10)
      P.seas1<-get.PTP(P = P_si1, dim.list = dims, i = match("si", model.type))*penalties[n-1 + p.n] #2494*2494 matrix (여기서 n-1+p.n = 3이므로, penalties[n-1 + p.n]=25가 된다)
      #P.seas1<-get.PTP(P = P_si1, dim.list = dims, i = n)*penalties[n-1 + p.n] #2494*2494 matrix (여기서 n-1+p.n = 3이므로, penalties[n-1 + p.n]=25가 된다) (오리지널)
      #p.n<-p.n+1 #p.n은 construct network model matrix and network penalty에서 0으로 설정했다, 즉 여기에서 1이 됨
      #P_si2<-kronecker.spam(diag.spam(ncol(P_spat)), P_seasonal) #1130*1130 matrix (앞의 행렬이 113*113, 뒤의 행렬이 10*10)
      #P.seas2<-get.PTP(P = P_si2, dim.list = dims, i = n)*penalties[n-1+p.n] #지금은 n-1+p.n = 4가 됨, 따라서 penalties[n-1+p.n] = 5, P.seas2 또한 2494*2494 matrix
      #PTP<-PTP + P.seas1 + P.seas2
    }
    p.n<-p.n+1 #p.n은 construct network model matrix and network penalty에서 0으로 설정했다, 즉 여기에서 1이 됨
    if("si"%in%model.type){
      P_si2<-kronecker.spam(diag.spam(ncol(P_spat)), P_seasonal) #1130*1130 matrix (앞의 행렬이 113*113, 뒤의 행렬이 10*10)
      P.seas2<-get.PTP(P = P_si2, dim.list = dims, i = match("si", model.type))*penalties[n-1+p.n] #지금은 n-1+p.n = 4가 됨, 따라서 penalties[n-1+p.n] = 5, P.seas2 또한 2494*2494 matrix
      #P.seas2<-get.PTP(P = P_si2, dim.list = dims, i = n)*penalties[n-1+p.n] #지금은 n-1+p.n = 4가 됨, 따라서 penalties[n-1+p.n] = 5, P.seas2 또한 2494*2494 matrix (오리지널)
      PTP<-PTP + P.seas1 + P.seas2
    }
    
    # Long term trend component and penalty ("t")
    n<-n+1	#n <- 5
    if("t"%in%model.type){
      # Long term trend component and penalty ("t")
      #n<-n+1	#n <- 5
      #years<-decimal_date(as.Date(dates))
      #나의 변경작업
      years <- year(dates) + (decimal_date(as.Date(dates2))-2013)
      X.list<-c(X.list, list(trend = as.spam(bbase(years, nseg = (knts[2] - 3))))) #dim as.spam(bbase(years, nseg = (knts[2] - 3))) = 4372*10, (knts[2] - 3)=7, #bbase: Construct B-spline basis
      model.mat<-cbind.spam(model.mat, X.list$trend) #dim(model.mat) : 4372 1264
      P_trend<-as.spam(diff(diag(ncol(X.list$trend)), diff = 1)) #long term trend의 경우에는 difference를 1로 주었다 #dim(P_trend) : 9 * 10
      P.trend<-get.PTP(P = P_trend, dim.list = dims, i = match("t", model.type))*penalties[n-1+p.n] #2494*2494 matrix
      #P.trend<-get.PTP(P = P_trend, dim.list = dims, i = n)*penalties[n-1+p.n] #2494*2494 matrix (원래)
      PTP<-PTP + P.trend
    }
    
    if("ti"%in%model.type){
      # create network-long term trend interaction and penalty ("ti")
      n<-n+1 #n <- 6
      tren_spat<-box3.prod(X.list$spatial, X.list$trend) #4372*1130 matrix
      X.list<-c(X.list, list(tren_spat)) #tren_spat 추가
      model.mat<-cbind.spam(model.mat, tren_spat) #model.mat: 4372*2394 matrix
      P_ti1<-kronecker.spam(P_spat, diag.spam(ncol(P_trend))) #dim(diag.spam(ncol(P_trend))): 10*10 dim(P_spat): 112*113 따라서 dim(P_ti1): 1120*1130
      P_ti2<-kronecker.spam(diag.spam(ncol(P_spat)),  P_trend) #dim(P_trend): 9*10, dim(diag.spam(ncol(P_spat))): 113*113 따라서 dim(P_ti2): 1017:1130
      P.ti1<-get.PTP(P = P_ti1, dim.list = dims, i = match("ti", model.type))*penalties[n-1+p.n] #n-1+p.n=6 따라서 #penalties[n-1+p.n]은 25 dim(P.ti1)은 2494*2494
      #P.ti1<-get.PTP(P = P_ti1, dim.list = dims, i = n)*penalties[n-1+p.n] #n-1+p.n=6 따라서 #penalties[n-1+p.n]은 25 dim(P.ti1)은 2494*2494 (오리지널)
      #p.n<-p.n+1 #p.n을 2로 업데이트
      #P.ti2<-get.PTP(P = P_ti2, dim.list = dims, i = n)*penalties[n-1+p.n] #n-1+p.n=7 따라서 #penalties[n-1+p.n]은 5 dim(P.ti2)은 2494*2494
      #PTP<-PTP + P.ti1 + P.ti2
    }
    p.n<-p.n+1 #p.n을 2로 업데이트
    if("ti"%in%model.type){
      P.ti2<-get.PTP(P = P_ti2, dim.list = dims, i = match("ti", model.type))*penalties[n-1+p.n] #n-1+p.n=7 따라서 #penalties[n-1+p.n]은 5 dim(P.ti2)은 2494*2494
      #P.ti2<-get.PTP(P = P_ti2, dim.list = dims, i = n)*penalties[n-1+p.n] #n-1+p.n=7 따라서 #penalties[n-1+p.n]은 5 dim(P.ti2)은 2494*2494 (오리지널)
      PTP<-PTP + P.ti1 + P.ti2
    }
    
    if("ts"%in%model.type){
      #  create trend-seasonal interaction component ("ts")
      n<-n+1 #n <- 7
      tren_seas <- box3.prod(X.list$trend, X.list$seasonal) #4372*100 matrix #dim(X.list$trend)=4372*10, dim(X.list$seasonal)=4372*10
      model.mat<-cbind.spam(model.mat, tren_seas) ##?? 4372*2494 matrix
      X.list<-c(X.list, list(tren_seas)) #tren.ses 추가
      P_ts1<-kronecker.spam(P_trend, diag.spam(ncol(P_seasonal))) #90*100 matrix
      P_ts2<-diag(ncol(P_trend)) %x% P_seasonal #P_seasonal은 10*10 matrix diag(ncol(P_trend))또한 10*10 matrix
      P.ts1<-get.PTP(P = as.spam(P_ts1), dim.list = dims, i = match("ts", model.type))*penalties[n-1 + p.n] #n-1 + p.n이 8이므로 penalties[n-1 + p.n]은 25, 2494*2494 matrix
      #P.ts1<-get.PTP(P = as.spam(P_ts1), dim.list = dims, i = n)*penalties[n-1 + p.n] #n-1 + p.n이 8이므로 penalties[n-1 + p.n]은 25, 2494*2494 matrix (오리지널)
      #p.n<-p.n+1 #p.n은 3으로 
      #P.ts2<-get.PTP(P = as.spam(P_ts2), dim.list = dims, i = n)*penalties[n-1 + p.n] #2494*2494 matrix
      #PTP<-PTP + P.ts1 + P.ts2
    }
    p.n<-p.n+1 #p.n은 3으로 
    if("ts"%in%model.type){
      P.ts2<-get.PTP(P = as.spam(P_ts2), dim.list = dims, i = match("ts", model.type))*penalties[n-1 + p.n] #2494*2494 matrix
      #P.ts2<-get.PTP(P = as.spam(P_ts2), dim.list = dims, i = n)*penalties[n-1 + p.n] #2494*2494 matrix (오리지널)
      PTP<-PTP + P.ts1 + P.ts2
    }
    
    gc() #gc: Garbage Collection
    
    # FIT THE MODEL using penalised least squares
    ident<-b.diag(lapply(X.list, crossprodspam)) #b.diag: list element들에 bdiag.spam을 넣는 함수 dim(ident)는 2494*2494
    XTX<-t(model.mat) %*% model.mat #model mat는 4372*2494 matrix
    info<-XTX + PTP  +  0.001*ident + diag.spam(c(0.0001, rep(0.1, nrow(XTX)-1))) #맨 첫번째에만 0.0001을 넣고 나머지에만 0.1을 넣는 diagonal 행렬
    #(CAUTION: 만약 "m"을 넣지 않을 경우 info matrix가 spam matrix가 되지 않는다는 문제 생김)
    U<-chol.spam(info, pivot = TRUE) #chol.spam: Cholesky Factorization for Sparse Matrices
    beta_hat<-backsolve(U, forwardsolve(U, t(model.mat) %*% log(response))) #log(TON)에 대해 관심 있으므로 log(response)를 쓴 듯, beta_hat은 2494 벡터
    fit<-model.mat %*% beta_hat #4372*1 행렬
    resids<- fit - log(response) #fit에서 log 뺀 값 (뒤에서 쓰임)
    vec<-forwardsolve(U, t(model.mat))
    pdof<-rowSums(vec^2)
    #   rm("vec");gc()
    dof<-sum(pdof)
    sigma.sq<-sum((log(response) - fit)^2)/(length(response) - dof)
    AICc<-log(sigma.sq) + 1 + ((2 + 2*dof)/(length(response) - dof -2)) #이것을 나중에 모형 평가시 사용하면 된다
    parameters<-as.vector(beta_hat) #2494 vec (sum(dims)가 2494)
    locs<-cumsum(dims) #dims [1]    1  113   10 1130   10 1130  100
    if(length(model.type) == 1){locs<-c(1, locs)}
    comps<-vector("list")	
  #}
  
  if(plot.fig==TRUE){
    # PLOT THE SPATIAL COMPONENT FOR ALL TIME POINTS, INCLUDING INTERCEPT
    #{
    M<-parameters[1]
    comps<-c(comps, list(mean = rep(M, dims[1])))
    if("m"%in%model.type){
      comps<-c(comps, list(spatial = parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]))
      #comps<-c(comps, list(spatial = parameters[(locs[2-1]+1):locs[2]])) #(오리지널)
      fit.point<-comps$spatial
      main.c<-paste("Spatial component ", " DoF = ", 
                    round(sum(pdof[((locs[match("m", model.type)-1]+1):locs[match("m", model.type)])]), 1), sep="")
      #main.c<-paste("Spatial component ", " DoF = ", 
      #              round(sum(pdof[((locs[2-1]+1):locs[2])]), 1), sep="") #오리지널
      n.cols<-121
      #n.cols<-131
      ngrid<-n.cols
      #brks<-range(log(response)) #response보다 많이 차이날 수 있다
      brks <- range(fit.point[as.numeric(TweedPredPoints$StreamUnit)])
      brks<-seq(brks[1], brks[2], length = ngrid + 1)
      col.nums<-cut(fit.point[as.numeric(TweedPredPoints$StreamUnit)], breaks = brks, labels = FALSE)
      palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black"))
      main.cols<-palette(n.cols)
      par(mar = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
      plot(TweedPredPoints$Longitude,TweedPredPoints$Latitude, 
           pch = 20, col = main.cols[col.nums],
           cex = TweedPredPoints$Weights, bty = "n", xlab = "",
           #ylab = "", xaxt = "n", yaxt = "n", ylim = c(55.2, 55.9))
           ylab = "", xaxt = "n", yaxt = "n", ylim = range(TweedPredPoints$Latitude))
      #}
    }
    
    if("s"%in%model.type){
      # PLOT THE SEASONAL COMPONENT WITH PARTIAL RESIDUALS
      #{
      big.zero<-zero.wrap(model.type=model.type, component="s", dims=dims, knts=knts)
      vec<-forwardsolve(U, t(big.zero))
      se.s<-sqrt(colSums(vec^2)*sigma.sq)
      dys.mns<-c(31,28,31,30,31,30,31,31,30,31,30,31)
      mn.dys<-cumsum(c(1,dys.mns)[1:12]) 
      s.bas<-cSplineDes(mn.dys[1]/365, knots = (0:knts[1])/knts[1], ord=4)
      s.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]]
      #s.par<-parameters[(locs[3 - 1] + 1):locs[3]] #오리지널
      S<-s.par %*% t(s.bas)
      comps<-c(comps, list(seasonal = rep(S, 298)))
      seas<-(1:365)/365
      seas.bas<-cSplineDes(1:365/365, knots = 0:knts[1]/knts[1], ord=3)
      seas.comp<-seas.bas %*% s.par
      y.day<-yday(dates)	
      main.s<-"Seasonal component"
      main.s<-paste(main.s, " DoF = ", round(sum(pdof[((locs[match("s", model.type)-1]+1):locs[match("s", model.type)])]), 1), sep="") 
      #main.s<-paste(main.s, " DoF = ", round(sum(pdof[((locs[3-1]+1):locs[3])]), 1), sep="") #오리지널
      which.day<-function(n){seas.comp[n]}
      #partial.residuals<-resids+apply(as.matrix(yday(dates)), 1, which.day)
      partial.residuals<-resids+apply(as.matrix(yday(dates2)), 1, which.day) #윤달에 맞춰 변경? (이것이 y가 된다)
      #day.resids<-yday(dates)
      day.resids<-yday(dates2) #윤달에 맞춰 변경 (이것이 x가 된다)
      par(mar=c(5,5,2,2))
      plot(day.resids,partial.residuals,  pch = ".", xlab = "Day in year", 
           ylim = c(-1.5, 1.5), ylab = "log Nitrate concentrations (mg/l)", 
           type = "n", main = "", cex.axis = 2, cex.lab=2)
      polygon(c(1:365, 365:1), c( 2*se.s+seas.comp, rev(-2*se.s+seas.comp)),
              col = "grey", border = NA)
      points(day.resids,partial.residuals,  col = "brown",pch = ".", cex = 2)
      lines(1:365, seas.comp, xlab = "Day in year", lwd = 3)
      #}
    }
    
    if("t"%in%model.type){
      # PLOT THE TREND COMPONENT
      #{
      years<-year(dates)
      min.y<-min(years)
      max.y<-max(years)
      y.seq<-seq(min.y, max.y+1, by = 1/12)[-1]
      n.t<-which(model.type == "t")
      #t.bas<-bbase(2004, nseg = (knts[2] - 3), xl = min(years), xr = max(years))
      t.bas<-bbase(2014, nseg = (knts[2] - 3), xl = min(years), xr = max(years))
      t.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]]
      #t.par<-parameters[(locs[5 - 1] + 1):locs[5]] #오리지널
      Tr<-t.bas %*% t.par
      comps<-c(comps, list(trend = rep(Tr, 298)))
      tren<-seq(min(years), max(years))
      tren.bas<-bbase(tren, xl = min(years), xr = max(years), nseg = (knts[2] - 3))
      tren.comp<-tren.bas %*% t.par  
      big.zero<-zero.wrap(yrs=(min(years):max(years)), model.type=model.type, component="t", dims=dims, knts=knts)
      vec<-forwardsolve(U, t(big.zero))
      se.t<-sqrt(colSums(vec^2)*sigma.sq)
      main.t<-"Trend component"
      main.t<-paste(main.t, " DoF = ", round(sum(pdof[((locs[match("t", model.type)-1]+1):locs[match("s", model.type)])]), 1), sep="")
      #main.t<-paste(main.t, " DoF = ", round(sum(pdof[((locs[5-1]+1):locs[5])]), 1), sep="") #오리지널		
      which.yr<-function(n){tren.comp[n]}
      yr.index<-year(dates) %% min(year(dates)) + 1
      partial.residuals<-resids+apply(as.matrix(yr.index), 1, which.yr) #이것이 y가 됨
      par(mar=c(5,5,2,2))
      plot(year(dates), partial.residuals, pch = ".", xlab = "Year of observation", 
           ylab = "log Nitrate concentrations (mg/l)", ylim = c(-2, 2), 
           type = "n", main = "", cex.axis = 2, cex.lab=2)
      polygon(c(min(years):max(years),max(years):min(years)),
              c(2*se.t + tren.comp, rev( -2*se.t + tren.comp)), 
              col = "grey", border = NA)
      points(year(dates), partial.residuals, pch = ".", cex = 0.7)
      lines(tren, tren.comp, lwd = 2)
      #}
    }
    
    if("ts"%in%model.type){
      # PLOT THE INTERACTION BETWEEN TREND AND SEASONALITY
      #{
      min.yr<-min(year(dates))
      max.yr<-max(year(dates))
      yers<-seq(min.yr, max.yr+1, length.out = 50)
      n.ts<-which(model.type == "ts")
      ones<-rep(1, length = 50)
      s.bas<-ones %x% cSplineDes(seq(1,365,length.out=50)/365, knots = 0:knts[1]/knts[1], ord=4)
      t.bas<-bbase(yers, nseg = (knts[2] - 3)) %x% ones
      ts.bas<-box.prod(list(t.bas, s.bas))
      ts.par<-parameters[(locs[n.ts - 1] + 1):locs[n.ts]] #locs[6]+1 부터 locs[7] 까지
      ts.fitted<-matrix(as.vector(ts.bas %*% ts.par), nrow = 50, byrow = T)
      par(mar = c(5, 4, 4, 2))
      filled.contour(ts.fitted, xlab = "Trend", ylab = "Seasonal",axes = T)
      #}
    }
    
    # PLOT THE FITTED VALUES WITH TIME AT A SINGLE MONITORING STATION
    #{
    #res<-model.out$res
    parameters<-as.vector(beta_hat)
    if(length(model.type) == 1){locs<-c(1, locs)}
    comps<-vector("list", length(model.type))
    min.yr<-min(year(dates))
    max.yr<-max(year(dates))
    yrs<-seq(min.yr, max.yr) 
    n.yrs<-length(yrs)
    day.seq<-seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs) 
    
    # local segment specific means
    if("m"%in%model.type){
      seg.means<-parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]] #locs[2-1]+1는 2, locs[2]는 114
    }else{
      seg.means<-rep(parameters[1], length(station)) #m이 없을 경우 segment specific means를 정의한다
    }
    #seg.means<-parameters[(locs[2-1]+1):locs[2]] #locs[2-1]+1는 2, locs[2]는 114 #오리지널
    #station.mean<-seg.means[2] #여기서 2가 seg.number
    #나의 변경
    for(kk in 1: length(station)){
      station.mean<-seg.means[station[kk]]
      
      main=0; inter=0;
      if("c"%in%model.type){
        C<-matrix(station.mean, nrow = n.yrs, ncol = 365/5)
        main=main+C
      }
      
      if("m"%in%model.type){
        # catchment wide constant
        M<-parameters[1]
        M<-matrix(M, nrow = n.yrs, ncol = 365/5)
        main=main+M
      }
      
      if("s"%in%model.type){
        # seasonal component
        s.bas<-cSplineDes(seq(5,365,by=5)/365, knots = 0:knts[1]/knts[1], ord=4)
        s.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]] 
        #s.par<-parameters[(locs[3 - 1] + 1):locs[3]] #오리지널
        S<-s.par %*% t(s.bas)
        S<-matrix(S, nrow = n.yrs, ncol = 365/5, byrow = T)	
        main=main+S
      }
      
      if("si"%in%model.type){
        # seasonal spatial interaction component
        all.si.par<-parameters[(locs[match("si", model.type) - 1] + 1):locs[match("si", model.type)]]
        #all.si.par<-parameters[(locs[4 - 1] + 1):locs[4]] #오리지널
        si.par.stat<-all.si.par[(station*ncol(s.bas) + 1):((station+1)*ncol(s.bas))]
        si.add<-s.bas %*% si.par.stat
        SI<-matrix(si.add, nrow = n.yrs, ncol = 365/5, byrow = T)	
        inter=inter+SI
      }
      
      if("t"%in%model.type){
        # trend additive component
        yers<-seq(min.yr, max.yr+1, by=5/365)[-1]
        t.bas<-bbase(yers, nseg = (knts[2] - 3))
        t.par<-parameters[(locs[match("t", model.type) - 1] + 1):locs[match("t", model.type)]]
        #t.par<-parameters[(locs[5 - 1] + 1):locs[5]] #오리지널
        Tr<-t.bas %*% t.par
        Tr<-matrix(Tr, nrow = n.yrs, ncol = 365/5, byrow = T)
        main=main+Tr
      }
      
      if("ti"%in%model.type){
        # trend spatial interaction component
        all.ti.par<-parameters[(locs[match("ti", model.type) - 1] + 1):locs[match("ti", model.type)]]
        #all.ti.par<-parameters[(locs[6 - 1] + 1):locs[6]] #오리지널
        ti.par.stat<-all.ti.par[(station*ncol(t.bas) + 1):((station+1)*ncol(t.bas))]
        ti.add<-t.bas %*% ti.par.stat
        TI<-matrix(ti.add, nrow = n.yrs, ncol = 365/5, byrow = T)	
        inter=inter+TI
      }
      
      if("ts"%in%model.type){
        # trend seasonal interaction
        all.ts.par<-parameters[(locs[match("ts", model.type) - 1] + 1):locs[match("ts", model.type)]]
        #all.ts.par<-parameters[(locs[7 - 1] + 1):locs[7]] #오리지널
        s.bas.ext<-matrix(1, ncol = 1,  nrow = n.yrs) %x% s.bas
        ts.bas<-box3.prod(as.spam(t.bas), as.spam(s.bas.ext))
        ts.add<-ts.bas %*% all.ts.par
        TS<-matrix(ts.add, nrow = n.yrs, ncol = 365/5, byrow = T)
        inter=inter+TS
      }
      
      #main<-S + Tr + C + M
      #inter<-TI + SI + TS
      #obs.ind<-which(response.locs == 2) #이게  unique(response.locs) 보고 바꾸어야
      obs.ind<-which(response.locs == station[kk])
      dates.full<-dates[obs.ind]
      dates.stat<-decimal_date(as.Date(dates.full))
      vals.stat<-response[obs.ind]
      ind<-which(year(dates[obs.ind]) > 2000)
      vals.stat<-vals.stat[ind]
      dates.stat<-decimal_date(as.Date(dates.full[ind]))
      par(mar = c(4,4,2,0.5))
      plot(dates.stat, log(vals.stat), pch = ".", cex = 1, type = "n",ylim = range(log(vals.stat)),
           xlab = "Time", ylab = "log(Nitrate concentrations)",cex.lab = 1, cex.axis=1, main=paste("Segment ", station[kk]))
      lines(day.seq, c(t(main)), type = "l", col = "green", lwd=3)
      lines(day.seq, c(t(inter))+ c(t(main)), type = "l", col = "red", lwd=3)
      points(dates.stat, log(vals.stat))
      # dev.off()
      #나의 추가: legend
      legend("bottomleft", col=c("green", "red"), lty=c(1,1), c("Main", "With Interaction"))
      #}
    }
    
  } 
  if(use.optim==TRUE){
    return(AICc)
  }else{
    return(list(AICc=AICc, beta_hat=beta_hat, fit=fit, resids=resids, model.mat=model.mat, XTX=XTX, PTP=PTP, penalties=penalties, realweights=realweights, adjacency=adjacency, TweedData=TweedData, TweedPredPoints=TweedPredPoints))
  } 
  
}

