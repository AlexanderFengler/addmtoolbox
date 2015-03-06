#' Run model by unique trial for one set of parameter values
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Run model by unique trial
#' @return numeric variable storing the negative log likelihood for given parameter set
#' \code{addm_run_by_trial_dynamic}
#' @export
#' @param eye.dat data.table storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and a unique trial ids (id).
#' @param choice.dat data.table storing the item valuations (v1,v2...), reaction times in ms (rt), decisions (decision) and unique trial ids (id).
#' @param model.parameters vector with the four core addm parameters in order (drift.rate, theta, sd, non.decision.time).
#' @param timestep integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param state.step numeric variable between 0 and 1, that indicate how finegrained the vertial grid is supposed to be

addm_run_by_trial_dynamic = function(choice.dat = data.table(v1 = 0,v2 = 0, id = 0),
                                     eye.dat = data.table(fixloc = 0, fixdur = 0, fixnr = 1, id = 0),
                                     model.parameters = c(0.002,0.5,0.07,0),
                                     timestep = 10,
                                     state.step = 0.1){

  # INITIALIZATION OF PARAMETERS -------------------------------------------------------------------------------------------------------
  drift.rate = model.parameters[1]
  theta = model.parameters[2]
  cur.sd = model.parameters[3]
  non.decision.time = model.parameters[4]
  #-------------------------------------------------------------------------------------------------------------------------------------

  # OTHER INITIALIZATIONS --------------------------------------------------------------------------------------------------------------
  # Get set size of current data set
  cur.set_size = length(grep("^v[0-9]$",names(choice.dat)))

  # Initialize amount of trials supplied
  len.trials = length(choice.dat[,id])

  # Matrix that stores all valuations by trial
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=length(choice.dat[,id]))

  # Define position of vector where valuations start
  start.pos = which(names(choice.dat) == "v1")

  # fill up valuations matrix
  for (i in seq(cur.set_size)){
    valuations[i,] = choice.dat[[start.pos + i - 1]]
  }

  # Decision Vector
  decisions = choice.dat[,decision]
  #-------------------------------------------------------------------------------------------------------------------------------------

  # RUN MODEL --------------------------------------------------------------------------------------------------------------------------
  likelihoods = rep(0.5,len.trials)
  trial.cnt = 1
  eye.row.cnt = 1
  len.eye = length(eye.dat[,fixloc])
  eye.mat = as.matrix(eye.dat %>% select(-id,-condition_id))
  eye.mat[,5] = eye.mat[,5] / timestep
  cur.maxfix = 0

  rts = choice.dat[,rt]
  rts = rts / timestep

  while(eye.row.cnt < len.eye){
    cur.maxfix[1] = eye.mat[eye.row.cnt,4]

    # Run Model
     likelihoods[trial.cnt] = dynamicaddm(cur.sd,
                                          theta,
                                          drift.rate,
                                          non.decision.time,
                                          decisions[trial.cnt],
                                          valuations[,trial.cnt],
                                          eye.mat[eye.row.cnt:(eye.row.cnt + cur.maxfix - 1),1], # fixation locations
                                          eye.mat[eye.row.cnt:(eye.row.cnt + cur.maxfix - 1),5], # fixation durations
                                          rts[trial.cnt],
                                          state.step)

    trial.cnt[1] = trial.cnt + 1
    eye.row.cnt[1] = eye.row.cnt + cur.maxfix
  }
  # CONTINUE BY CALCULATING LOG LIKELIHOOD----------------------------------------------------------------------------------------------

  total.log.lik = c(drift.rate,
                    theta,
                    cur.sd,
                    non.decision.time,
                    0,
                    (-1)*sum(log(likelihoods)))

  print(total.log.lik)
  return(total.log.lik)
  # ------------------------------------------------------------------------------------------------------------------------------------
}
