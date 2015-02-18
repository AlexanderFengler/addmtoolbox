# Author: Alexander Fengler
# Date: January 16th 2015
# Purpose: This is the script from which to run the aDDM optimization from now on

# According to some model parameters that we chose we need to adjust a few columns in our main data
# ----------------------------------------------------------------------------------------------


addm_fit_artificial = function(nr.reps = 2000,timesteps.ms = 10,model.type = 'nomem'){

library(data.table)
library(tidyr)
library(dplyr)


# Create List that contains all necessary information to run the aDDM optimization in one place
# ---------------------------------------------------------------------------------------------
train.dat = list(subjects = c(-1),
                 set.sizes = c(0),
                 choice.dat = addm.choice.dat.fake,
                 eye.dat = addm.eye.dat.fake,
                 drifts = seq(0.002,0.01,0.002),      #default seq(0.00,0.03,0.005)
                 thetas = seq(0.2,1,0.2),             #default seq(0.0,1,0.1)
                 sds = seq(0.02,0.1,0.02),            #default seq(0,0.2,0.025)
                 non.decision.times = 0,
                 timesteps.ms = timesteps,
                 nr.reps = nr.reps,
                 model.type = model.type,
                 output.type = "Fake",
                 fixation.model="FakePath",
                 allow.extension = 1,
                 allow.fine.grid = 1,
                 generate = 0)
# ----------------------------------------------------------------------------------------------

# Run the model
# ----------------------------------------------------------------------------------------------
# Source Model function
source('temp/aDDMOpti.R')

# Model Function
aDDMOpti(train.dat$subjects,
         train.dat$set.sizes,
         train.dat$choice.dat,
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
         train.dat$allow.extension,
         train.dat$allow.fine.grid,
         train.dat$generate)
}
# ----------------------------------------------------------------------------------------------

