addm_greadsearch_foreach = function(cur.choice.dat,
                                    cur.eye.dat,
                                    parameter.matrix,
                                    cur.subject,
                                    cur.set_size,
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
  if (output.type == "Opti"){

    sink(cur.log.file,append=TRUE)

    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% add_trial(cur.choice.dat,
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

    out = data.table(Subject = cur.subject,
                     Set.size = cur.set_size,
                     Drift.Rate = temp[,1],
                     Theta = temp[,2],
                     SD = temp[,3],
                     Non.decision.time = temp[,4],
                     Nr.Reps = temp[,5],
                     Log.Likelihood = temp[,6]);
    return(out)
  }
  # --------------------------------------------------------------------------------------------------------------------------------------

  # RUN FIT BY CONDITION  ----------------------------------------------------------------------------------------------------------------
  if (output.type == "Fake"){

    sink(cur.log.file,append=TRUE)

    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% addm_hist(cur.choice.dat,
                                                                    cur.eye.dat,
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

    out = data.table(Subject = cur.subject,
                     Set.size = cur.set_size,
                     Drift.Rate = temp[,1],
                     Theta = temp[,2],
                     SD = temp[,3],
                     Non.decision.time = temp[,4],
                     Nr.Reps = temp[,5],
                     Log.Likelihood = temp[,6]);

    return (out)
  }
  #-------------------------------------------------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------------------------------------------------
}
