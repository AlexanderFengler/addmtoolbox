# Author: Alexander Fengler
# Date: January 16th 2015
# Purpose: This is the script from which to run the aDDM optimization from now on

# According to some model parameters that we chose we need to adjust a few columns in our main data
# ----------------------------------------------------------------------------------------------
library(data.table)
library(tidyr)
library(dplyr)

# Parameters:
timesteps = 10 #in ms // represents the timestep size with which the aDDM shall propagate
rtbinsize = 100 # in ms // represent the rtbin.size that we bin the data into

choice = as.data.table(core.clean.choices.liking[core.clean.choices.liking$istrain == 1,])
eye = as.data.table(core.clean.eyetracking.general[core.clean.eyetracking.general$istrain == 1,])

setkey(choice,Trialid)
setkey(eye,Trialid)
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

# Get two additions to RT, rtup -- the upper limit of the RT-Bin the current trials belongs to, rtdown -- the lower limit of the RT-Bin the current trial belongs to
choice$rtup = ((choice$RT_no_sacc + rtbinsize) %/% rtbinsize) * rtbinsize
choice$rtdown = choice$rtup - rtbinsize


# To make sure that the simulation has enough data to run to the end of the RT-Bin that a given trial is assigned to (real fixation-duration data stops at a random point in this RT-Bin) we need to adjust the length of last fixationas to accomodate simulation runs until at least the end of the assigned RT-Bins

#Step 1 get difference between current RT and end of RT-Bin
choice$rtdiff = choice$rtup - choice$RT_no_sacc
eye$rt.add.last = 0
for(trial in trials){
  eye[J(trial),rt.add.last:=choice[J(trial),rtdiff]]
  print(trial)
}

# Step 2
# We adjust the last fixation length for each trial such that we have data until the end of the rt.bin that the particular trial is assigned to
# Moreover we add a little more to allow the aDDM to go on after the real rt.bin is passed (to assign the simulation run as a failure)
eye$Duration[eye$Last.fix.binary == 1] = eye$Duration[eye$Last.fix.binary == 1] + eye$rt.add.last[eye$Last.fix.binary == 1] + rtbinsize

# Finalize the two data.table
choice = choice[,list(Subject,Set_size,Snack_1,Snack_2,Snack_3,Snack_4,Snack_5,Snack_6,Snack_7,Snack_8,Snack_picked,RT_no_sacc,rtup,rtdown,Trialid)]
setnames(choice,c("subject","set_size","v1","v2","v3","v4","v5","v6","v7","v8","decision","rt","rtup","rtdown","trialid"))
eye = eye[,list(Subject,Set_size,Item_attended,Fixation_nr,Duration,Trialid)]
setnames(eye,c("subject","set_size","fixloc","fixnr","fixdur","trialid"))

# Create List that contains all necessary information to run the aDDM optimization in one place
# ----------------------------------------------------------------------------------------------
core.clean.train.dat = list(subjects = c(0),
                            set.sizes = c(4),
                            choice.dat = choice,
                            eye.dat = eye,
                            drifts = seq(0.00,0.03,0.005),     #default seq(0.00,0.03,0.005)
                            thetas = seq(0.0,1,0.1),           #default seq(0.0,1,0.1)
                            sds = seq(0,0.2,0.025),            #default seq(0,0.2,0.025)
                            non.decision.times = 0,
                            timesteps.ms = 10,
                            nr.reps = 5000,
                            model.type = 'nomem',
                            output.type = "Opti",
                            fixation.model="Normal",
                            allow.extension = 1,
                            allow.fine.grid = 1,
                            generate = 0)
# ----------------------------------------------------------------------------------------------

# Run the model
# ----------------------------------------------------------------------------------------------
# Source Model function
source('temp/aDDMOptiSimple.R')

# Model Function
aDDMOpti(core.clean.train.dat$subjects,
         core.clean.train.dat$set.sizes,
         core.clean.train.dat$choice.dat,
         core.clean.train.dat$eye.dat,
         core.clean.train.dat$drifts,
         core.clean.train.dat$thetas,
         core.clean.train.dat$sds,
         core.clean.train.dat$non.decision.times,
         core.clean.train.dat$timesteps,
         core.clean.train.dat$nr.reps,
         core.clean.train.dat$model.type,
         core.clean.train.dat$output.type,
         core.clean.train.dat$fixation.model,
         core.clean.train.dat$allow.extension,
         core.clean.train.dat$allow.fine.grid,
         core.clean.train.dat$generate)
# ----------------------------------------------------------------------------------------------

