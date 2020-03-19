data_predicted_r <- data_predicted - data_predicted0
data_predicted0_r <- data_predicted0 - data_predicted0
data_predicted1_r <- data_predicted1 - data_predicted0
data_predicted2_r <- data_predicted2 - data_predicted0

zlims <- range(c(data_predicted_r, data_predicted0_r, data_predicted1_r, data_predicted2_r))+c(-0.05, 0.05)

#pdf("result_TOC3_residual.pdf", 7, 7)
#png("result_TOC3_residual.png", 700, 700)
par(family = 'sans') 
par(mar=c(3.1,2.1,3.1,1.1))
par(mfrow=c(2,2))
scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted0_r[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(a) Raw - Raw, log(TOC)", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825))
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], data, main="(a) Raw, TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3), bigplot=c(0,0,1,1), xlim=range(TweedPredPoints$Longitude), ylim=range(TweedPredPoints$Latitude))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted1_r[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(b) ODonnell - Raw", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], denoising_ODonell, main="(b) ODonnell, TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted_r[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(c) S-Lifting (M) - Raw", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825))
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], result_denoise, main="(c) S-Lifting (M), TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))

scatter_fill(TweedPredPoints$Longitude, TweedPredPoints$Latitude ,data_predicted2_r[TweedPredPoints$StreamUnit], pch=16, cex=TweedPredPoints$Weights, xlab="", ylab="", zlim=zlims, main="(d) S-Lifting (N) - Raw", cex.main=1.5, smallplot=c(0.8,0.85,0.65,0.825), plot.legend=TRUE, axes=TRUE, y.axes=FALSE)
points(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], pch=22)
#quilt.plot(example_network@obspoints@SSNPoints[[1]]@point.coords[,1], example_network@obspoints@SSNPoints[[1]]@point.coords[,2], result_nlt$aveghat, main="(d) S-Lifting (N), TOC", cex=2, add=T, zlim=range(c(data, data_predicted, data_predicted0, data_predicted1)), xlab="", ylab="",smallplot=c(0.15,0.2,0.15,0.3))

dev.off()




