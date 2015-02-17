# Author: Alexander Fengler
# Date: Jan 19th 2015
# Purpose: Generate full aDDM Outputs and store to file

# Load necessary functions
#source('Analysis_Visualization/aDDM_Simulation/load_packages_scripts_aDDM.R')
#load.aDDM.packages()
#load.aDDM.scripts()

# Load the necessary scripts to make this one running from scratch 
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppZiggurat)
library(data.table)

sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_full_output.cpp")
sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_full_output_withmem.cpp")
sourceCpp(file="Analysis_Visualization/aDDM_Simulation/aevacc_full_output_withmem_zeronoise.cpp")
source('Analysis_Visualization/aDDM_Simulation/aDDM.R')

# Define set size
set.size = 8
timesteps = 10 # in milliseconds
nr.reps = 1000
model.type = 'mem'
fixation.model = "Normal"
output.type = "Full"
generate = 0

simple.ddm = 0 # binary tell us whether to look for the best solution with theta = 1 only (1) or across all parameter-combinations present (0)

# Read in Log.Likelihood file and extract optimal parameters
cur.file = paste("Sim_logs/fits_mem_global_train/loglik_mem_all_setsize_", toString(set.size), ".txt",sep="")
cur.log.liks = data.table(read.table(cur.file,header = TRUE))

# Generate all input parameters to aDDM function
if (simple.ddm == 0){
max.pos = which(cur.log.liks[,Log.Likelihood] == max(cur.log.liks[,Log.Likelihood]))
core.parameters = as.vector(unlist(cur.log.liks[max.pos,list(Drift.Rate,Theta,SD,Non.decision.time)]))
}

if (simple.ddm == 1){
  cur.log.liks = cur.log.liks[Theta == 1,]
  max.pos = which(cur.log.liks[,Log.Likelihood] == max(cur.log.liks[,Log.Likelihood]))
  core.parameters = as.vector(unlist(cur.log.liks[max.pos,list(Drift.Rate,Theta,SD,Non.decision.time)]))
}
# GET THE INPUT DATA READY
# ----------------------------------------------------------------------------------------------
choice = as.data.table(core.clean.choices.liking[core.clean.choices.liking$istrain == 0 & 
                                                   core.clean.choices.liking$Set_size == set.size,])
eye = as.data.table(core.clean.eyetracking.general[core.clean.eyetracking.general$istrain == 0 &
                                                     core.clean.eyetracking.general$Set_size == set.size,])
stretch = as.data.table(core.clean.choices.liking.stretched[core.clean.choices.liking.stretched$istrain == 0 &
                                                              core.clean.choices.liking.stretched$Set_size == set.size,])

trials = unique(eye$Trialid)

setkey(choice,Trialid)
setkey(eye,Trialid)
setkey(stretch,Trialid)


# Transfer the "eligible.post.empirical" column from stretch to choice
temp.elig = rep(0,8)

for (trial in trials){
  temp.elig[1:set.size] = stretch[J(trial),Eligible.post.empirical]
  for (i in seq(8)){
    choice[J(trial),paste("eligible", toString(i), sep=""):=temp.elig[i]]
  }
  print(trial)
}


# We need a fixation duration vector that we sample from for post empirical fixations
# Fixation length is taken from real refixations, however capped at extremes (sub-quantiles from 0.05 until 0.95)
fixdur.vec = eye$Duration[eye$Refix.no.last.fix == 1 & eye$Set_size == set.size]
dur.lims = as.vector(quantile(fixdur.vec,c(0.05,0.5,0.95,1)))
fixdur.vec = fixdur.vec[fixdur.vec > dur.lims[1] & fixdur.vec < dur.lims[3]] 
fixdur.vec = timesteps*round(fixdur.vec/timesteps)

# We store the fixation duration vector in eye data.table in column "fixdur.samples"
len = length(eye[,Duration])
eye[,fixdur.sampled:=sample(fixdur.vec,len,replace=TRUE)]

# ----------------------------------------------------------------------------------------------
# Round all fixation times according chosen timestep 
eye$Duration =  timesteps * round(eye$Duration/timesteps)

rts = eye %>% group_by(Trialid)  %>% summarize(rt = sum(Duration))

setkey(rts,Trialid)

trials = unique(eye$Trialid)

# Update the rt's in choice data.table according to sum of all rounded durations for consistency
for(trial in trials){ 
  choice[J(trial),RT_no_sacc:=rts[J(trial),rt]]
  print(trial)
}

# Finalize the two data.table 
choice = choice[,list(Subject,
                      Set_size,
                      Snack_1,
                      Snack_2,
                      Snack_3,
                      Snack_4,
                      Snack_5,
                      Snack_6,
                      Snack_7,
                      Snack_8,
                      Snack_picked,
                      RT_no_sacc,
                      eligible1,
                      eligible2,
                      eligible3,
                      eligible4,
                      eligible5,
                      eligible6,
                      eligible7,
                      eligible8,
                      Trialid)]

setnames(choice,c("subject",
                  "set_size",
                  "v1",
                  "v2",
                  "v3",
                  "v4",
                  "v5",
                  "v6",
                  "v7",
                  "v8",
                  "decision",
                  "rt",
                  "eligible1",
                  "eligible2",
                  "eligible3",
                  "eligible4",
                  "eligible5",
                  "eligible6",
                  "eligible7",
                  "eligible8",
                  "trialid"))

eye = eye[,list(Subject,
                Set_size,
                Item_attended,
                Fixation_nr,
                Duration,
                fixdur.sampled,
                Trialid)]

setnames(eye,c("subject",
               "set_size",
               "fixloc",
               "fixnr",
               "fixdur",
               "fixdur.sampled",
               "trialid"))

setkey(choice,trialid)
setkey(eye,trialid)
# ----------------------------------------------------------------------------------------------

# Run aDDM
out =  aDDM(choice,
                  eye,
                  core.parameters,
                  nr.reps,
                  model.type,
                  output.type,
                  fixation.model,
                  timesteps,
                  generate)

# Store output as file
if (simple.ddm == 0){
  out.file = paste("Data_RAW/aDDM_model_output/addmmem_opti_simplefix_",toString(set.size),'.txt',sep='')
}

if (simple.ddm == 1){
  out.file = paste("Data_RAW/aDDM_model_output/ddm_opti_simplefix_",toString(set.size),'.txt',sep='')
}


write.table(out,out.file,sep=' ',row.names=FALSE)

