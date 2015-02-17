# Generating parameter combinations
aDDM.optimization.fake.drift.1000 = function(set.sizes,subjects){
  
  #source('Scripts/Dealing_with_Workspace/Refresh_workspace_to_core_objects_main.R')
  
  library(plyr)
  library(data.table)
  library(Rcpp)
  library(doMC)
  library(foreach)
  library(iterators)
  registerDoMC(7)
  
  #load("MAIN_EXPERIMENT.RData")
  
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_full_4.cpp")
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_full_6.cpp")
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_full_8.cpp")
  
  source('Scripts/aDDM_Simulation/aDDM_fake_1000.R')
  
  ##############################################################################################################################
  ################################################# PART 1 #####################################################################
  ##############################################################################################################################
  ############################# INITIALIZE PARAMETER COMBINATIONS ##############################################################
  ##############################################################################################################################
  
  drift.rates = c(0.001,0.0012,0.0014,0.0016,0.0018,0.002,0.0022,0.0024,0.0026,0.0028,0.003)
  variances = 0.0002 #c(0.0001,0.00012,0.00014,0.00016,0.00018,0.0002,0.00022,0.00024,0.00026,0.00028,0.0003) 
    non.decision.times = 0 
  thetas = 0.5 #c(0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7)
    
    subjects.amount = length(subjects)
  set.sizes.amount = length(set.sizes)
  out = list(0)
  
  total_iteration_cnt = 1
  set_size_cnt = 1
  
  parameter.combs.mat = matrix(rep(0,6*length(subjects)*length(set.sizes)*length(drift.rates)*length(variances)*length(non.decision.times)*length(thetas)),ncol=6)
  
  row_cnt = 1
  for (subject in subjects){
    for (set.size in set.sizes){
      for (drift.rate in drift.rates){
        for (variance in variances){
          for (theta in thetas){
            for (non.decision.time in non.decision.times){
              parameter.combs.mat[row_cnt, 1:6] = c(subject,set.size,drift.rate,theta,variance,non.decision.time)
              row_cnt = sum(row_cnt,1)
            }
          }
        }
      }
    }
  }
  print("parameter combs ready")
  
  ita = iter(parameter.combs.mat,by="row")
  
  ###############################################################################################################################
  ####################################### PART 2 ################################################################################
  ###############################################################################################################################
  ####################################################### RUN ALL PARAMETERS IN PARALLEL WITH FOREACH ###########################
  ###############################################################################################################################
  
  out[[1]]= foreach (i = ita,.combine='rbind') %dopar% aDDM_fake_1000(i)
  
  temp = do.call(rbind,out)
  
  out = data.frame(Subject = temp[,1],
                   Set.size = temp[,2],
                   Drift.Rate = temp[,3],
                   Variance = temp[,4],
                   Theta = temp[,5],
                   Non.decision.time = temp[,6],
                   Log.Likelihood = temp[,7])
  
  return (out)
}