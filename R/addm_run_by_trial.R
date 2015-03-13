#' Run model by unique trial for one set of parameter values
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Run model by unique trial
#' @return numeric variable storing the negative log likelihood for given parameter set
#' \code{addm_run_by_trial}
#' @export
#' @param eye.dat data.table storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and a unique trial ids (id).
#' @param choice.dat data.table storing the item valuations (v1,v2...), reaction times in ms (rt), decisions (decision) and unique trial ids (id).
#' @param model.parameters vector with the four core addm parameters in order (non.decision.time, drift, sd, theta, gamma, boundary-parameters).
#' @param nr.attributes integer providing the amount of attributes we consider per item
#' @param nr.reps integer that tells the function how many simulation runs to use.
#' @param model.type string that indicates which version of the model to run. 'standard' for traditional (a)ddm, or 'memnoise' when memory effects shall be allowed.
#' @param timestep integer number that provides the timestep-size that is used in the simulations (in ms).

addm_run_by_trial = function(choice.dat = data.table(v1 = 0,v2 = 0, id = 0),
                             eye.dat = data.table(fixloc = 0, fixdur = 0, fixnr = 1, id = 0),
                             model.parameters = c(0,0.002,0.07,1),
                             nr.attributes = 1,
                             boundaryfun = 1,
                             nr.reps = 1000,
                             timestep = 10,
                             model.type = 'standard'){

  # Compute Boundaries -----------------------------------------------------------------------------------------------------------------
  if (class(boundaryfun) == 'function'){
    # COMPUTE
  }
  # ------------------------------------------------------------------------------------------------------------------------------------

  # Initialize parameter theta ---------------------------------------------------------------------------------------------------------
  # only for readability, because it is used later on
  theta = model.parameters[4]
  # ------------------------------------------------------------------------------------------------------------------------------------

  # OTHER INITIALIZATIONS --------------------------------------------------------------------------------------------------------------
  # Get set size of current data set
  cur.set_size = length(grep("^v[0-9]*",names(choice.dat)))
  #return(cur.set_size)

  # Initialize amount of trials supplied
  len.trials = length(choice.dat[,id])

  # Matrix that stores all valuations by trial
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=length(choice.dat[,id]))

  # Define position of vector where valuations start
  start.pos =  grep("^v1",names(choice.dat))[1]

  # fill up valuations matrix
  for (i in seq_len(cur.set_size)){
    valuations[i,] = choice.dat[[start.pos + i - 1]]
  }

  # Decision Vector
  decisions = choice.dat[,decision] - 1 # Minus one because in C++ vectorsstart at indice zero
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZE EVIDENCE ACCUMULATION FUNCTION ------------------------------------------------------------------------------------------
  if (model.type == 'standard'){
    if (cur.set_size / nr.attributes == 2){
      if (nr.attributes == 1){
        if (theta == 1){
          aevacc = evacc2_by_trial
        } else {
          aevacc = aevacc2_by_trial
        }
      } else {
        aevacc = aevaccma2_by_trial
      }
    } else {
      if (nr.attributes == 1){
        aevacc = aevacc_by_trial
      } else {
        stop('You attempted to run a set of > 2 items which have multiple attributes each (by trial): This is not implemented at the moment!')
      }
    }
  } else if (model.type == 'memnoise'){
    if (cur.set_size / nr.attributes == 2){
      stop('You attempted to run a set of 2 items with memory effects (by trial): This is not implemented at the moment!')
    } else {
    aevacc = aevacc_by_trial_memnoise
    }
  }
  # ------------------------------------------------------------------------------------------------------------------------------------
  #return(aevacc)
  # RUN MODEL --------------------------------------------------------------------------------------------------------------------------
  success.counts = rep(0,len.trials)
  trial.cnt = 1
  eye.row.cnt = 1
  len.eye = length(eye.dat[,fixloc])
  eye.mat = as.matrix(eye.dat %>% select(-id,-condition_id))
  cur.maxfix = 0
  cur.rtbin.up = choice.dat[,rtup]
  cur.rtbin.down = choice.dat[,rtdown]
  while(eye.row.cnt < len.eye){
    cur.maxfix[1] = eye.mat[eye.row.cnt,4]

    # Run Model
    success.counts[trial.cnt] = aevacc(model.parameters,
                                       cur.rtbin.up[trial.cnt],
                                       cur.rtbin.down[trial.cnt],
                                       decisions[trial.cnt],
                                       valuations[,trial.cnt],
                                       nr.attributes,
                                       eye.mat[eye.row.cnt:(eye.row.cnt + cur.maxfix - 1),1], # fixation locations
                                       eye.mat[eye.row.cnt:(eye.row.cnt + cur.maxfix - 1),2], # fixation durations
                                       nr.reps,
                                       timestep)

    trial.cnt[1] = trial.cnt[1] + 1
    eye.row.cnt[1] = eye.row.cnt[1] + cur.maxfix
  }
  # CONTINUE BY CALCULATING LOG LIKELIHOOD ---------------------------------------------------------------------------------------------
  success.counts[success.counts == 0] = 1/(1+nr.reps)
  if (nr.attributes == 1){
    total.log.lik = c(model.parameters,
                      (-1)*sum(log(success.counts/nr.reps)))
  } else {
    total.log.lik = c(model.parameters,
                      (-1)*sum(log(success.counts/nr.reps)))
  }

  print(total.log.lik)
  return(total.log.lik)
  # ------------------------------------------------------------------------------------------------------------------------------------
}
