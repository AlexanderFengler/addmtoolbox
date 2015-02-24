addm_gridsearch_foreach = function(conditions,dat,
                                   cur.eye.dat,
                                   choice.dat,
                                   parameter.matrix,
                                   output.type,
                                   nr.reps,
                                   model.type,
                                   fixation.model,
                                   cur.log.file,
                                   timestep.ms,
                                   generate){

  # RUN MODEL WITH PROVIDED PARAMETERS --------------------------------------------------------------------------------------------------
  out = list(0)
  ita = iter(parameter.matrix,by="row")


  # RUN FIT BY TRIAL --------------------------------------------------------------------------------------------------------------------
  if (output.type == "trial"){

    sink(cur.log.file,append=TRUE)

    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% addm_by_trial(cur.choice.dat,
                                                                        cur.eye.dat,
                                                                        i,
                                                                        nr.reps,
                                                                        model.type,
                                                                        fixation.model,
                                                                        timestep.ms)

    cur.sink.num = sink.number()
    for (z in 1:cur.sink.num){
      sink(NULL)
    }

    temp = do.call(rbind,out)

    out = data.table(drift = temp[,1],
                     theta = temp[,2],
                     sd = temp[,3],
                     non.decision.time = temp[,4],
                     nr.reps = temp[,5],
                     loglik = temp[,6]);
    return(out)
  }
  # --------------------------------------------------------------------------------------------------------------------------------------

  # RUN FIT BY CONDITION  ----------------------------------------------------------------------------------------------------------------
  if (output.type == "condition"){

    sink(cur.log.file,append=TRUE)

    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% addm_by_condition(conditions.dat,
                                                                            cur.eye.dat,
                                                                            cur.choice.dat,
                                                                            i,
                                                                            nr.reps,
                                                                            model.type,
                                                                            output.type,
                                                                            fixation.model,
                                                                            timestep.ms,
                                                                            generate)

    cur.sink.num = sink.number()
    for (z in 1:cur.sink.num){
      sink(NULL)
    }

    temp = do.call(rbind,out)

    out = data.table(drift = temp[,1],
                     theta = temp[,2],
                     sd = temp[,3],
                     non.decision.time = temp[,4],
                     nr.reps = temp[,5],
                     loglik = temp[,6]);
    return (out)
  }
  #-------------------------------------------------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------------------------------------------------
}
