# GENERAL INFORMATION ABOUT HOW TO DEAL WITH THE aDDM simulation functions
#############################################################################################

# There are 4 types of aDDM that can be run here

# Two different output types:
# 1 Loglikelihood for optimization purposes
# 2 Full Output for purposes of detailed observation of a particular parameter setup

# Two differen fixation models
# 1 The current best model
# 2 Random fixations and Durations


# There are four different "aevacc..." functions which produce the output for the four aDDM models respectively
# In addition there is a "fake" version of the aDDM, which is designed to take in simulation data for parameter recovery


# IMPORTANT PARAMETER VALUES
# fixation.model ("Normal","Random")
# output.type ("Real","Fake","Full")
# model.type ("mem","nomem")

# Krajbich parameters

# Drift.Rate = 0.000666
# Variance = 0.000196
# Theta = 0.3
# Non-decision_time = 0 (just not modelled)

# Example for using aDDM function
#core.aDDM.ranfix.output.frame = aDDM.optimization(unique(core.clean.choices.liking$Subject),c(4,6,8),0.000666,0.000196,0,0.3,"Full_Output",250,"Random")
#core.DDM.ranfix.output.frame = aDDM.optimization(unique(core.clean.choices.liking$Subject),c(4,6,8),0.000666,0.000196,0,1,"Full_Output",250,"Random")

aDDM.gridsearch = function(cur.choice.dat,
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

  # The parameter matrix needs to be puts into iterator form with the help of iter()

  ita = iter(parameter.matrix,by="row")

  ###################################### RUN ALL PARAMETERS IN PARALLEL WITH FOREACH ############################################
  ###############################################################################################################################
  #-------------------------------------------------------------------------------------------------------------------------------------

  out = list(0)

  if (output.type == "Opti"){

    sink(cur.log.file,append=TRUE)

    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% aDDM_fast(cur.choice.dat,
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


  if (output.type == "Fake"){

    sink(cur.log.file,append=TRUE)

    out[[1]] = foreach (i = ita,.combine='rbind') %dopar% aDDM(cur.choice.dat,
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
}
