#' Generate artificial data for addm model fits
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Generate artificial dataset
#' @return returns a list of three components. A 'data.table' that stores all unique choice set conditions. A 'data.table' that stores eyetracking.data adjusted to be usable for by trial fits. A 'data.table' that stores by trial choice data. An id variable links all three data.tables.
#' \code{addm_generate_by_condition_artificial} Returns artificial data that can be used for test addm model fits
#' @export
#' @param set.size Indicates the number of items that are allowed.
#' @param possible.valuations Vector storing valuations that single items can hold in a choice set.
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param nr.reps An integer number that tells the function how many simulation runs to use.
#' @param model.type A string that indicates which version of the model to run. 'standard' or 'memnoise' when memory effects shall be allowed.
#' @param fixation.model A string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (implemented) 'fixedpath' for following a predetermined fixation path with fixed durations (implemented).
#' @param core.parameters Vector that provide parameters used to generate artificial data from drift diffusion process. (drift,theta,sd,non deicision time)
addm_generate_by_condition_artificial = function(set.size = 2,
                                                 possible.valuations = c(0,1,2,3),
                                                 timestep= 10,
                                                 nr.reps = 2000,
                                                 model.type = "standard",
                                                 fixation.model = "fixedpath",
                                                 core.parameters = c(0.006,0.6,0.06,0)){

# GENERATE MODEL INPUT DATA --------------------------------------------------------------------------------------------------------------
  # Conditions
  val.dat = as.data.table(matrix(sample(possible.valuations,25*set_size,replace=TRUE),ncol=set_size))
  setnames(val.dat,names(val.dat),tolower(names(val.dat)))
  ids = data.table(id = 1:25)
  conditions =  cbind(val.dat, ids)

  # Eyetracking
  cur.eye.dat =  data.table(fixloc = 0, fixnr=0, fixdur=0, id = 1:25)

  test.dat = list(conditions = conditions,
                  eye.dat = cur.eye.dat,
                  choice.dat = 0,
                  timestep = timestep,
                  nr.reps = nr.reps,
                  model.type = model.type,
                  fixation.model = fixation.model,
                  output.type = 'fit',
                  core.parameters = core.parameters, # order of core parameter: drift.rate, theta, sd, non.decision.time
                  generate = 1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test
#------------------------------------------------------------------------------------------------------------------------------------------

# RUN MODEL TO GET CHOICE DATA ------------------------------------------------------------------------------------------------------------
choices = addm_by_condition(test.dat$conditions,
                            test.dat$eye.dat,
                            test.dat$choice.dat,
                            test.dat$core.parameters,
                            test.dat$nr.reps,
                            test.dat$model.type,
                            test.dat$output.type,
                            test.dat$fixation.model,
                            test.dat$timestep,
                            test.dat$generate)
#------------------------------------------------------------------------------------------------------------------------------------------
return(list(cur.choice.dat = choices,
            conditions = conditions,
            cur.eye.dat = cur.eye.dat))
}


