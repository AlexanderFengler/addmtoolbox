#' Runs parallel grid search for addm fits using foreach loops
#' \code{addm_gridsearch_foreach} Returns data.table with log likelihoods and corresponding parameter combinations
#' #' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Run grid search parallel
#' @return Returns data.table with log likelihoods and corresponding parameter combinations
#' @export
#' @inheritParams addm_fit_grid
#' @param parameter.matrix Matrix that represents the parameter space which is looped over. (drift,theta,sd,non.decision.time)

addm_support_gridsearch_foreach = function(conditions,dat = data.table(v1 = 0, v2 = 0, id = 0),
                                           eye.dat = data.table(fixloc = 0, fixdur = 0, fixnr = 1, id = 0),
                                           choice.dat = data.table(v1 = 0, v2 = 0, rt = 0, decision = 0, id = 0),
                                           parameter.matrix = c(0.006,0.6,0.06,0),
                                           fit.type = 'condition',
                                           nr.reps = 2000,
                                           model.type = 'standard',
                                           fixation.model = 'fixedpath',
                                           log.file = 'defaultlog.txt',
                                           timestep = 10){

  # Initialize iterator and output list --------------------------------------------------------------------------------------------------
  out = list(0)
  ita = iter(parameter.matrix,by="row")
  # --------------------------------------------------------------------------------------------------------------------------------------

  # RUN FIT ------------------------------------------------------------------------------------------------------------------------------
  sink(cur.log.file,append=TRUE)

  if (fit.type == "trial"){
    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% addm_by_trial(choice.dat,
                                                                        eye.dat,
                                                                        i,
                                                                        nr.reps,
                                                                        model.type,
                                                                        fixation.model,
                                                                        timestep)

  } else if (fit.type == 'condition'){
    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% addm_by_condition(conditions.dat,
                                                                            eye.dat,
                                                                            choice.dat,
                                                                            i,
                                                                            nr.reps,
                                                                            model.type,
                                                                            output.type = 'fit',
                                                                            fixation.model,
                                                                            timestep,
                                                                            generate = 0)

  }
  # --------------------------------------------------------------------------------------------------------------------------------------

  # Reset sinks to restore output to console ---------------------------------------------------------------------------------------------
  cur.sink.num = sink.number()
  for (z in 1:cur.sink.num){
    sink(NULL)
  }
  # --------------------------------------------------------------------------------------------------------------------------------------

  # Format and return output -------------------------------------------------------------------------------------------------------------
  temp = do.call(rbind,out)
  out = data.table(drift = temp[,1],
                   theta = temp[,2],
                   sd = temp[,3],
                   non.decision.time = temp[,4],
                   nr.reps = temp[,5],
                   loglik = temp[,6]);
  return(out)
  # --------------------------------------------------------------------------------------------------------------------------------------
}
