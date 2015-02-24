# Author: Alexander Fengler
# Date: January 16th 2015
# Purpose: This is the script from which to run the aDDM optimization from now on

# According to some model parameters that we chose we need to adjust a few columns in our main data
# ----------------------------------------------------------------------------------------------
addm_fit_artificial = function(conditions = data.table(v1 = c(1,2,3), v2 = c(3,2,1), id = c(1,2,3)),
                               eye.dat = 0,
                               choice.dat = data.table(v1 = c(1,2,3), v2 = c(3,2,1), id = c(1,2,3), rt = 0, decision = 1),
                               nrsims = 1000,
                               timesteps.ms = 10,
                               addmmodel = 'nomem'){

library(data.table)
library(tidyr)
library(dplyr)

# CREATE MODEL FITTING INPUT DATA--------------------------------------------------------------
train.dat = list(conditions = conditions,
                 eye.dat = eye.dat,
                 choice.dat = choice.dat,
                 drifts = seq(0.002,0.01,0.002),      #default seq(0.00,0.03,0.005)
                 thetas = seq(0.2,1,0.2),             #default seq(0.0,1,0.1)
                 sds = seq(0.02,0.1,0.02),            #default seq(0,0.2,0.025)
                 non.decision.times = 0,
                 timesteps.ms = timesteps.ms,
                 nr.reps = nrsims,
                 model.type = addmmodel,
                 output.type = "condition",
                 fixation.model="fakepath",
                 allow.fine.grid = 0,
                 generate = 0)
# ----------------------------------------------------------------------------------------------

# RUN THE MODEL---------------------------------------------------------------------------------
# Source fitting function
source('temp/addm_fit_grid.R')

# Model Function
addm_fit_grid(train.dat$choice.dat,
              train.dat$eye.dat,
              train.dat$drifts,
              train.dat$thetas,
              train.dat$sds,
              train.dat$non.decision.times,
              train.dat$timesteps,
              train.dat$nr.reps,
              train.dat$model.type,
              train.dat$output.type,
              train.dat$fixation.model,
              train.dat$allow.fine.grid,
              train.dat$generate)
# ------------------------------------------------------------------------------------------------
}

