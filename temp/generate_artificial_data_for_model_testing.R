# Author: Alexander Fengler
# Project: aDDM_4_6_8
# Date: Jan 14th 2015


# This script generates all necessary data needed as input to generate artificial aDDM data
# We can then use this to test the model fitting procedures which currently reside in file: run_aDDM_optimization

# RUN THIS PART ONLY IF CORRESPONDING FILE IS LOST FROM WORKSPACE -- OTHERWISE IT IS REUSABLE
#------------------------------------------------------------------------------------------------------------------------------------------
# Generating a list that collect all necessary input data to generate artificial aDDM from the aDDM.R scripts (aDDM() function)
  
core.clean.model.test.dat = list(items_4 = data.table(v1=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v2=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v3=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v4=sample(c(0,1,2,3),100,replace=TRUE),
                                                      trialid=1:100),
                                 items_6 = data.table(v1=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v2=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v3=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v4=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v5=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v6=sample(c(0,1,2,3),100,replace=TRUE),
                                                      trialid=1:100),
                                 items_8 = data.table(v1=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v2=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v3=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v4=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v5=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v6=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v7=sample(c(0,1,2,3),100,replace=TRUE),
                                                      v8=sample(c(0,1,2,3),100,replace=TRUE),
                                                      trialid=1:100),
                                 decisions = data.table(decision = 0,
                                                        trialid = 1:100),
                                 trialids= 1:100,
                                 rts = data.table(rt = 0,
                                                  trialid = 1:100),
                                 fixations = data.table(fixloc = 0,fixnr=0,fixdur=0,trialid = 1:100),
                                 timestep.ms = 10,
                                 nr.reps = 2000,
                                 fixation.model = "FakePath",
                                 output.type = "Fake",
                                 core.parameters = c(0.05,0.5,0.1,0))  # order of core parameter: drift.rate,sd,theta,non.decision.time
#------------------------------------------------------------------------------------------------------------------------------------------

# GENERATING ARTIFICIAL OUTPUT
#------------------------------------------------------------------------------------------------------------------------------------------

core.clean.artificial.choice.table.4 = aDDM(core.clean.model.test.dat$items_4,
                                            core.clean.model.test.dat$decisions,
                                            core.clean.model.test.dat$fixations,
                                            core.clean.model.test.dat$rts,
                                            core.clean.model.test.dat$trialids,
                                            core.clean.model.test.dat$core.parameters,
                                            core.clean.model.test.dat$nr.reps,
                                            core.clean.model.test.dat$output.type,
                                            core.clean.model.test.dat$fixation.model,
                                            core.clean.model.test.dat$timestep.ms,
                                            1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test


core.clean.artificial.choice.table.6 = aDDM(core.clean.model.test.dat$items_6,
                                            core.clean.model.test.dat$decisions,
                                            core.clean.model.test.dat$fixations,
                                            core.clean.model.test.dat$rts,
                                            core.clean.model.test.dat$trialids,
                                            core.clean.model.test.dat$core.parameters,
                                            core.clean.model.test.dat$nr.reps,
                                            core.clean.model.test.dat$output.type,
                                            core.clean.model.test.dat$fixation.model,
                                            core.clean.model.test.dat$timestep.ms,
                                            1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test


core.clean.artificial.choice.table.8 = aDDM(core.clean.model.test.dat$items_8,
                                            core.clean.model.test.dat$decisions,
                                            core.clean.model.test.dat$fixations,
                                            core.clean.model.test.dat$rts,
                                            core.clean.model.test.dat$trialids,
                                            core.clean.model.test.dat$core.parameters,
                                            core.clean.model.test.dat$nr.reps,
                                            core.clean.model.test.dat$output.type,
                                            core.clean.model.test.dat$fixation.model,
                                            core.clean.model.test.dat$timestep.ms,
                                            1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test

#------------------------------------------------------------------------------------------------------------------------------------------


