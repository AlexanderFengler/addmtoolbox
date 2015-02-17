# Make a matrix out of a set of parameter combinations supplied

#   Example of possible parameters
#   drift.rates = c(0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007)
#   variances = c(0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007)  
#   non.decision.times = 0 
#   thetas = c(0.6,0.7,0.8,0.9,1)  


param.combs = function(drift.rates,sds,thetas,non.decision.times){
  
# Initialize Matrix
  parameter.combs.mat = matrix(rep(0,4*length(drift.rates)*length(sds)*length(non.decision.times)*length(thetas)),ncol=4)
  
# Fill with all combinations of parameters provided
  row_cnt = 1

for (drift.rate in drift.rates){
  for (sd in sds){
    for (theta in thetas){
      for (non.decision.time in non.decision.times){
        parameter.combs.mat[row_cnt, 1:4] = c(drift.rate,theta,sd,non.decision.time)
        row_cnt = sum(row_cnt,1)
      }
    }
  }
}

return(parameter.combs.mat)
}