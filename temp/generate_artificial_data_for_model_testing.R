# Author: Alexander Fengler
# Project: aDDM_4_6_8
# Date: Jan 14th 2015

# This script generates all necessary data needed as input to generate artificial aDDM data
# We can then use this to test the model fitting procedures which currently reside in file: run_aDDM_optimization

# Generating a list that collects all necessary input data to generate artificial aDDM from the aDDM.R scripts (aDDM() function)

addm_generate_artificial = function(set_size, possible_valuations){

# GENERATE MODEL INPUT DATA --------------------------------------------------------------------------------------------------------------
val.dat = as.data.table(matrix(sample(possible_valuations,25*set_size,replace=TRUE),ncol=set_size))
setnames(val.dat,names(val.dat),tolower(names(val.dat)))

decisions = data.table(decision = 0,
                       trialid = 1:25)

rts = data.table(rt = 0)

addm.choice.dat.fake <<- cbind(val.dat,decisions,rts)
addm.eye.dat.fake <<- data.table(fixloc = 0,fixnr=0,fixdur=0,trialid = 1:25)

core.clean.model.test.dat = list(choice.dat = addm.choice.dat.fake,
                                 eye.dat = addm.eye.dat.fake,
                                 timestep.ms = 10,
                                 nr.reps = 10000,
                                 model.type = "nomem",
                                 fixation.model = "FakePath",
                                 output.type = "Fake",
                                 core.parameters = c(0.006,0.6,0.06,0), # order of core parameter: drift.rate, sd, theta, non.decision.time
                                 generate = 1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test
#------------------------------------------------------------------------------------------------------------------------------------------

# RUN MODEL -------------------------------------------------------------------------------------------------------------------------------
aDDM(core.clean.model.test.dat$choice.dat,
     core.clean.model.test.dat$eye.dat,
     core.clean.model.test.dat$core.parameters,
     core.clean.model.test.dat$nr.reps,
     core.clean.model.test.dat$model.type,
     core.clean.model.test.dat$output.type,
     core.clean.model.test.dat$fixation.model,
     core.clean.model.test.dat$timestep.ms,
     core.clean.model.test.dat$generate)
#------------------------------------------------------------------------------------------------------------------------------------------
}


