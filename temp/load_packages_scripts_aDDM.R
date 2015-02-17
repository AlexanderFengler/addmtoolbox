# Loading all relevant scripts/packages for aDDM optimization
#-------------------------------------------------------------------------------------------------------------------------------------
# Packages
load.aDDM.packages = function(){
  #source('Analysis_Visualization/Dealing_with_Workspace/Refresh_workspace_to_core_objects_main.R')
  library(dplyr)
  library(data.table)
  library(Rcpp)
  library(RcppZiggurat)
  library(data.table)
  library(doMC)
  library(foreach)
  library(iterators)
  registerDoMC(8)
}

load.aDDM.scripts = function(){
  #load("MAIN_EXPERIMENT_1_27.RData")
  
  # Initialize Evidence accumulation functions for all versions of the model that are possible at the moment
  #---------------------------------------------------------------------------------------------------------
  # Normal Output, Normal Fixations
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_fast.cpp")
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_fast_withmem.cpp")
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_fast_withmem_zeronoise.cpp")
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc.cpp")
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_full_output.cpp")
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_full_output_withmem.cpp")
  sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_full_output_withmem_zeronoise.cpp")
  
  
  #----------------------------------------------------------------------------------------------------------
  
  
  # aDDM wrapper functinos
  source('Analysis_Visualization/aDDM_Simulation/aDDM.R')
  source('Analysis_Visualization/aDDM_Simulation/aDDM_fast.R')
}
#-------------------------------------------------------------------------------------------------------------------------------------
