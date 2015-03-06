#' Generate artificial data for addm model fits
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Generate artificial dataset
#' @return returns a list of three components. A 'data.table' that stores all unique choice set conditions. A 'data.table' that stores eyetracking.data adjusted to be usable for by trial fits. A 'data.table' that stores by trial choice data. An id variable links all three data.tables.
#' \code{addm_generate_artificial} Returns artificial data that can be used for test addm model fits
#' @export
#' @param set.size Indicates the number of items that are allowed.
#' @param possible.valuations Vector storing valuations that single items can hold in a choice set.
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param nr.reps An integer number that tells the function how many simulation runs to use.
#' @param model.type A string that indicates which version of the model to run. 'standard' or 'memnoise' when memory effects shall be allowed.
#' @param fixation.model A string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (implemented) 'fixedpath' for following a predetermined fixation path with fixed durations (implemented).
#' @param core.parameters Vector that provide parameters used to generate artificial data from drift diffusion process. (drift,theta,sd,non deicision time)
#' @param nr.conditions integer that provides the number of unique trial conditions to be generated
#' @param rtbinsize integer providing the binsize that reaction times will be sorted into
addm_generate_artificial = function(set.size = 2,
                                    possible.valuations = c(0,1,2,3),
                                    model.parameters = c(0.002,0.5,0.07,0),
                                    nr.reps = 100,
                                    timestep = 10,
                                    rtbinsize = 100,
                                    nr.conditions = 10,
                                    model.type = "standard"){
# GENERATE MODEL INPUT DATA --------------------------------------------------------------------------------------------------------------
  # Conditions
  val.dat = as.data.table(matrix(sample(possible.valuations,nr.conditions*set.size,replace=TRUE),ncol=set.size))
  setnames(val.dat,names(val.dat),tolower(names(val.dat)))
  ids = data.table(condition_id = 1:nr.conditions)
  conditions =  cbind(val.dat, ids)

  test.dat = list(choice.dat = 0,
                  conditions.dat = conditions,
                  model.parameters = model.parameters,
                  nr.reps = nr.reps,
                  timestep = timestep,
                  model.type = model.type,
                  output.type = 'fit',
                  fixation.model = 'fixedpath', # order of core parameter: drift.rate, theta, sd, non.decision.time
                  generate = 1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test
#------------------------------------------------------------------------------------------------------------------------------------------

# RUN MODEL TO GET CHOICE DATA ------------------------------------------------------------------------------------------------------------
choices = do.call(addm_run_by_condition, args = test.dat)
setkey(choices,condition_id)
setkey(conditions,condition_id)
choices = conditions[choices]
choices$id = 1:length(choices[,rt])
choices$rtup = floor((choices$rt + rtbinsize)/ rtbinsize ) * rtbinsize
choices$rtdown = floor(choices$rt / rtbinsize) * rtbinsize
#------------------------------------------------------------------------------------------------------------------------------------------

# GENERATE FIXATION DATA FROM CHOICES OUTPUT ----------------------------------------------------------------------------------------------
temp = choices
setkey(temp,id)
temp$nr.fix = ceiling(temp$rt / 400)
temp$last.dur = ((temp$rt / 400) - floor(temp$rt / 400)) * 400
temp = temp[last.dur == 0, last.dur := 400]
len.eye = sum(temp$nr.fix)
temp = as.matrix(temp[,list(nr.fix,last.dur,condition_id,id)])

eye.dat = matrix(rep(rep(0,6),len.eye),ncol=6)
fixpath = rep(c(seq(1:set.size)),200)
cnt = 1
eye.cnt = 1
nr.fix = 0
last.dur = 0
condition_id = 0
id = 0

while (eye.cnt < len.eye){
  nr.fix[1] = as.integer(round(as.numeric(temp[cnt, 1])))
  last.dur[1] = as.integer(round(as.numeric(temp[cnt, 2])))
  condition_id[1] = as.integer(round(as.numeric(temp[cnt, 3])))
  id[1] = as.integer(round(as.numeric(temp[cnt, 4])))

  eye.dat[eye.cnt:(eye.cnt + nr.fix - 1), 1] = 1:nr.fix # sequential fixation number column
  if (nr.fix > 1){
    eye.dat[eye.cnt:(eye.cnt + nr.fix - 2), 2] = 400 # fixation durations besides last of trial (along addm2_fixation_model_fixedpath())
  }
  eye.dat[(eye.cnt + nr.fix -1), 2] = last.dur # duration of last fixation
  eye.dat[eye.cnt:(eye.cnt + nr.fix - 1), 3] = fixpath[1:nr.fix] # fixation path (along addm2_fixation_model_fixedpath())
  eye.dat[eye.cnt:(eye.cnt + nr.fix - 1), 4] = condition_id # condition id
  eye.dat[eye.cnt:(eye.cnt + nr.fix - 1), 5] = id # by trial id
  eye.dat[eye.cnt:(eye.cnt + nr.fix - 1), 6] = nr.fix # amount of fixations in trial

  eye.cnt[1] = eye.cnt + nr.fix
  cnt[1] = cnt + 1
}

eye.dat.tab = data.table(id = eye.dat[,5],
                     condition_id = eye.dat[,4],
                     fixloc = eye.dat[,3],
                     fixdur = 400,
                     fixnr = eye.dat[,1],
                     max.fix = eye.dat[,6],
                     fixdur.orig = eye.dat[,2]) # fixdur.orig is needed for fits with dynamic model (in our artificial data case this is simply the same as fixdur)

# eye.dat.tab = eye.dat.tab %>%
#   group_by(id) %>%
#   mutate(rt.eye = sum(fixdur.orig))
#
# setkey(eye.dat.tab,id)
# setkey(choices,id)
#
# viewtest = choices[eye.dat.tab]
# viewtest$rt.diff = viewtest$rt - viewtest$rt.eye
#------------------------------------------------------------------------------------------------------------------------------------------

# Generate Artifical fixation da
#------------------------------------------------------------------------------------------------------------------------------------------
return(list(choice.dat = choices,
            conditions.dat = conditions,
            eye.dat = eye.dat.tab))
}
