# Generating parameter combinations
aDDM.optimization.fake = function(set.sizes,subjects){
  
  #source('Scripts/Dealing_with_Workspace/Refresh_workspace_to_core_objects_main.R')
  
  library(plyr)
  library(data.table)
  library(Rcpp)
  library(doMC)
  library(foreach)
  library(iterators)
  registerDoMC(7)
  
  #load("MAIN_EXPERIMENT.RData")
  
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_4.cpp")
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_6.cpp")
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_8.cpp")
  
  source('Scripts/aDDM_Simulation/aDDM_fake.R')
  
  ##############################################################################################################################
  ################################################# PART 1 #####################################################################
  ##############################################################################################################################
  ############################# INITIALIZE PARAMETER COMBINATIONS ##############################################################
  ##############################################################################################################################
  
  drift.rates = c(0.001,0.0015,0.002,0.0025,0.003)
  variances = c(0.0001,0.00015,0.0002,0.00025,0.0003)  
  non.decision.times = 0 
  thetas = c(0.3,0.4,0.5,0.6,0.7)
  
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
  
  #######
  # Initialize Log.txt
  #######
  
  writeLines(c(""), "Sim_logs/log_aDDM_fake_opti.txt")
  sink("Sim_logs/log_aDDM_fake_opti.txt",append=TRUE)
  
  out[[1]]= foreach (i = ita,.combine='rbind') %dopar% aDDM_fake(i)

  sink(NULL)
  
  
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