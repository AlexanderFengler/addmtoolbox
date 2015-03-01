#' Runs the model by trial
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Model run by trial
#' @return Returns a log likelihood value.
#' \code{addm_by_trial} Returns log likelihood
#' @export
#' @param eye.dat A 'data.frame' or 'data.table' storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and an id column (id). Can all be initialized as zero columns when a by condition fit is attempted.
#' @param choice.dat A 'data.frame' or  'data.table' storing the item valuations (v1,v2...) , reaction times in ms (rt), decisions as decision and an id column (id). A by trial form is assumed.
#' @param model.parameters A vector with the four core addm parameters in order (drift.rate,theta,sd,non.decision.time).
#' @param nr.reps An integer number that tells the function how many simulation runs to use.
#' @param model.type A string that indicates which version of the model to run. 'standard' or 'memnoise' when memory effects shall be allowed.
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).

addm_by_trial = function(choice.dat = data.table(v1 = 0,v2 = 0, id = 0),
                         eye.dat = data.table(fixloc = 0, fixdur = 0, fixnr = 1, id = 0),
                         model.parameters = c(0.006,0.6,0.06,0),
                         nr.reps = 2000,
                         model.type = 'standard',
                         timestep = 10){

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
  decisions = choice.dat[,decision] - 1 # Minus one because in C++ vectorsstart at indice zero
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZE EVIDENCE ACCUMULATION FUNCTION ------------------------------------------------------------------------------------------
  if (model.type == 'standard'){
    if (cur.set_size == 2){
      aevacc = aevacc2_by_trial
    } else {
      aevacc = aevacc_by_trial
    }
  } else if (model.type == 'memnoise'){
    aevacc = aevacc_by_trial_memnoise
  }
  # ------------------------------------------------------------------------------------------------------------------------------------

  # RUN MODEL --------------------------------------------------------------------------------------------------------------------------
  success.counts = rep(0,len.trials)
  cnt = 1

  cur.rtbin.up = choice.dat[,rtup]
  cur.rtbin.down = choice.dat[,rtdown]
  ids = choice.dat[,id]

  for (trial in ids){
    # Extract Empirical Fixations and Fixations Durations
    cur.fixations = eye.dat[J(trial),fixloc]
    cur.durations = eye.dat[J(trial),fixdur]

    # Run Model
    success.counts[cnt] = aevacc(nr.reps,
                                 cur.rtbin.up[cnt],
                                 cur.rtbin.down[cnt],
                                 decisions[cnt],
                                 cur.sd,
                                 theta,
                                 drift.rate,
                                 non.decision.time,
                                 timestep,
                                 valuations[,cnt],
                                 cur.fixations,
                                 cur.durations)

    cnt[1] = cnt + 1
  }

  # CONTINUE BY CALCULATING LOG LIKELIHOOD----------------------------------------------------------------------------------------------
  success.counts[success.counts == 0] = nr.reps/(nr.reps + 1)

  total.log.lik = c(drift.rate,
                    theta,
                    cur.sd,
                    non.decision.time,
                    nr.reps,
                    (-1)*sum(log(success.counts/nr.reps)))

  print(total.log.lik)
  return(total.log.lik)
  # ------------------------------------------------------------------------------------------------------------------------------------
}
