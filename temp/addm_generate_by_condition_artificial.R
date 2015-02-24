# Author: Alexander Fengler
# Project: aDDM_4_6_8
# Date: Jan 14th 2015

# This script generates all necessary data needed as input to generate artificial aDDM data
# We can then use this to test the model fitting procedures which currently reside in file: run_aDDM_optimization

# Generating a list that collects all necessary input data to generate artificial aDDM from the aDDM.R scripts (aDDM() function)

addm_generate_by_condition_artificial = function(set_size = 2,
                                                 possible_valuations = c(0,1,2,3),
                                                 timestep_ms = 10,
                                                 simnum = 2000,
                                                 model_type = "standard",
                                                 fixation_model = "fakepath",
                                                 output_type = "fit",
                                                 core_parameters = c(0.006,0.6,0.06,0)){

# GENERATE MODEL INPUT DATA --------------------------------------------------------------------------------------------------------------
  # Conditions
  val.dat = as.data.table(matrix(sample(possible_valuations,50*set_size,replace=TRUE),ncol=set_size))
  setnames(val.dat,names(val.dat),tolower(names(val.dat)))
  ids = data.table(id = 1:50)
  conditions =  cbind(val.dat, ids)

  # Eyetracking
  cur.eye.dat =  data.table(fixloc = 0, fixnr=0, fixdur=0, id = 1:50)

  test.dat = list(conditions = conditions,
                  eye.dat = cur.eye.dat,
                  timestep.ms = timestep_ms,
                  nr.reps = simnum,
                  model.type = model_type,
                  fixation.model = fixation_model,
                  output.type = output_type,
                  core.parameters = core_parameters, # order of core parameter: drift.rate, theta, sd, non.decision.time
                  generate = 1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test
#------------------------------------------------------------------------------------------------------------------------------------------

# RUN MODEL TO GET CHOICE DATA ------------------------------------------------------------------------------------------------------------
choices = addm_by_condition(test.dat$conditions,
                            test.dat$eye.dat,
                            0,
                            test.dat$core.parameters,
                            test.dat$nr.reps,
                            test.dat$model.type,
                            test.dat$output.type,
                            test.dat$fixation.model,
                            test.dat$timestep.ms,
                            test.dat$generate)
#------------------------------------------------------------------------------------------------------------------------------------------
return(list(cur.choice.dat = choices,
            conditions = conditions,
            cur.eye.dat = cur.eye.dat))
}


