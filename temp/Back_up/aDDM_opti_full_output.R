# Generating parameter combinations
aDDM.optimization.full.output = function(set.sizes,subjects){
  
  #source('Scripts/Dealing_with_Workspace/Refresh_workspace_to_core_objects_main.R')
  
  library(plyr)
  library(Rcpp)
  library(doMC)
  library(foreach)
  library(iterators)
  library(data.table)
  registerDoMC(7)
  
  #load("MAIN_EXPERIMENT_1_19.RData")
  
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_full_output_4.cpp")
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_full_output_6.cpp")
  sourceCpp(file="Scripts/aDDM_Simulation/aevacc_full_output_8.cpp")
  
  source('Scripts/aDDM_Simulation/aDDM_full_output.R')
  
  ##############################################################################################################################
  ################################################# PART 1 #####################################################################
  ##############################################################################################################################
  ############################# INITIALIZE PARAMETER COMBINATIONS ##############################################################
  ##############################################################################################################################
  drift.rates = c(0.0006666)
  variances = c(0.000196)  
  non.decision.times = 0 
  thetas = c(0.3)
  
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
  
  ita = iter(parameter.combs.mat,by="row")
  
  ###############################################################################################################################
  ####################################### PART 2 ################################################################################
  ###############################################################################################################################
  ####################################################### RUN ALL PARAMETERS IN PARALLEL WITH FOREACH ###########################
  ###############################################################################################################################
  writeLines(c(""), "Sim_logs/log_aDDM_opti.txt")
  sink("Sim_logs/log_aDDM_full_output.txt",append=TRUE)
  
  out[[total_iteration_cnt]]= foreach (i = ita,.combine='rbind') %dopar% aDDM.full.output(i)
  
  sink(NULL)
  sink(NULL)
  
  temp = do.call(rbind,out)
  
  out = as.data.frame(out)
  #out = data.frame(Subject = temp[,1],
    #               Set.size = temp[,2],
     #              Drift.Rate = temp[,3],
      #             Variance = temp[,4],
       #            Theta = temp[,5],
        #           Non.decision.time = temp[,6],
         #          Log.Likelihood = temp[,7])
  
  return (out)
}