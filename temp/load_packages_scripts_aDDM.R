# Loading all relevant scripts/packages for aDDM optimization
#-------------------------------------------------------------------------------------------------------------------------------------
# Packages
load.aDDM.packages = function(){
  library(dplyr)
  library(data.table)
  library(Rcpp)
  library(RcppZiggurat)
  library(doMC)
  library(foreach)
  library(iterators)
  cores = detectCores()
  registerDoMC(cores)
}

load.aDDM.scripts = function(){

  # Initialize Evidence accumulation functions for all versions of model
  sourceCpp(file="temp/aevacc_fast.cpp")
  sourceCpp(file="temp/aevacc_fast_withmem.cpp")
  sourceCpp(file="temp/aevacc_fast_withmem_zeronoise.cpp")
  sourceCpp(file="temp/aevacc_hist.cpp")
  sourceCpp(file="temp/aevacc_full_output.cpp")
  sourceCpp(file="temp/aevacc_full_output_withmem.cpp")
  sourceCpp(file="temp/aevacc_full_output_withmem_zeronoise.cpp")

  # aDDM wrapper functinos
  source('temp/aDDM.R')
  source('temp/aDDM_fast.R')
}
#-------------------------------------------------------------------------------------------------------------------------------------
