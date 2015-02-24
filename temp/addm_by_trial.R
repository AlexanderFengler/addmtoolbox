###### I AM GOING TO USE THE ADDM BASIS AND ADJUST IT TOWARDS NO ATTENTIONAL EFFECT

#   inputs needed:
# - item distribution matrix
# - each item distribution is assigned a fixation path matrix
# - fixation duration
# - drift rate
# - theta
# - standard deviation

# - number of repetitions used -- nr.reps
# - model.type -- "standard" or "memnoise"  two model types so far, with memory and without memory effects (memory effects currently theta = 0 for unseen items)
# - fixation.model -- 'Normal' or 'Random'
# - timesteps used -- in ms (the model accounts for it then)

addm_by_trial = function(cur.choice.dat,
                         cur.eye.dat,
                         core.parameters,
                         nr.reps,
                         model.type,
                         fixation.model,
                         timestep.ms){

  # INITIALIZATION OF PARAMETERS -------------------------------------------------------------------------------------------------------
  drift.rate = core.parameters[1]
  theta = core.parameters[2]
  cur.sd = core.parameters[3]
  non.decision.time = core.parameters[4]
  #-------------------------------------------------------------------------------------------------------------------------------------

  # OTHER INITIALIZATIONS --------------------------------------------------------------------------------------------------------------
  # Get set size of current data set
  cur.set_size = length(grep("^v[0-9]$",names(cur.choice.dat))) # this works only if the item columns are supplied according to the correct naming (v1,v2 etc.)

  # valuation is a matrix that already takes into account the drift rate to describe the per-timestep mean shift by item for each trial
  # is then supplied to the evidence accumulation function below
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=length(cur.choice.dat[,trialid]))

  # Define position of vector where valuations start
  start.pos = which(names(cur.choice.dat) == "v1")

  # fill up valuations matrix
  for (i in seq(cur.set_size)){
    valuations[i,] = cur.choice.dat[[start.pos + i - 1]]
  }

  # Denerate a decisions vector and adjust the number to accomodate difference in indexing between R and C++
  decisions = cur.choice.dat[,decision] - 1 # Minus one because in C++ vectorsstart at indice zero
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
  success.counts = rep(0,length(cur.choice.dat[,trialid]))
  cnt = 1

  cur.rtbin.up = cur.choice.dat[,rtup]
  cur.rtbin.down = cur.choice.dat[,rtdown]
  trialids = cur.choice.dat[,trialid]

  for (trial in trialids){
    # Extract Empirical Fixations and Fixations Durations
    cur.fixations = cur.eye.dat[J(trial),fixloc]
    cur.durations = cur.eye.dat[J(trial),fixdur]

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
                                 cur.durations,
                                 cur.set_size)

    cnt[1] = cnt + 1
  }

  # CONTINUE BY CALCULATING LOG LIKELIHOOD----------------------------------------------------------------------------------------------
  success.counts[success.counts == 0] = 0.99

  total.log.lik = c(drift.rate,
                    theta,
                    cur.sd,
                    non.decision.time,
                    nr.reps,
                    sum(log(success.counts/nr.reps)))

  print(total.log.lik)
  return(total.log.lik)
  # ------------------------------------------------------------------------------------------------------------------------------------
}
