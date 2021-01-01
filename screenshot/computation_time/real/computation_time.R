aaaa <- readRDS("~/Dropbox/Github/paper-StreamflowLifting/screenshot/computation_time/real/ListStreamSim113(sd1)nlt_withresidualwithtime.RDS")

result_mat <- c()
for(i in 1:length(aaaa)){
  result_mat <- rbind(result_mat, aaaa[[i]]$result_new_mat[,6])
}
colMeans(result_mat)

bbbb <- readRDS("~/Dropbox/Github/paper-StreamflowLifting/screenshot/computation_time/real/ListStreamSTPCA80(sd1)nlt_withresidualwithtime.RDS")

result_mat2 <- c()
for(j in 1:length(bbbb)){
  result_mat2 <- rbind(result_mat2, bbbb[[j]]$result_new_mat[,6])
}
colMeans(result_mat2)