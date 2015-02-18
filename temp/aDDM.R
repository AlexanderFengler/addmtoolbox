#   Inputs needed:
# - subject
# - set size
# - drift.rate
# - theta
# - standard deviation
# - non decision time

# - number of repetitions used -- nr.reps
# - fixation.model -- 'Normal'
#                  -- 'Fakepath' (used for fitting pregenerated datasets)
#                  -- 'Random'
# - output.type -- 'Normal' return log.likelihood against real data
#               -- 'Fake' return log.likelihood against a simulated data set
#               -- 'Full' returns the full data output for a specific model (every single simulation repetition as single row)
# - model.type -- "mem" use addm function with memory effects (whichever the best at the moment)
#               -- "nomem" use addm function without memory (normal Krajbich Type)
# - generate -- boolean variable indicating whether we run the function to generate a data frame (1) or to get a log likelihood (0)

aDDM = function(cur.choice.dat,
                cur.eye.dat,
                core.parameters,
                nr.reps,
                model.type,
                output.type,
                fixation.model,
                timestep.ms,
                generate){


  # INITIALIZATION OF PARAMETERS AND RELEVANT SUBSETS OF DATA FRAMES--------------------------------------------------------------------

  # Initialize model parameters
  drift.rate = core.parameters[1]
  theta = core.parameters[2]
  cur.sd = core.parameters[3]
  non.decision.time = core.parameters[4]
  #-------------------------------------------------------------------------------------------------------------------------------------

  # SOME MISCELLANEOUS VARIABLES THAT ARE UTILIZED LATER--------------------------------------------------------------------------------

  # First we need to define max.RT
  cur.max.RT = 40000 %/% timestep.ms

  # Generating column that provides a cutted (binned) version of reaction times
  increment.distances = 100
  cur.breaks = c(rev(seq(20000,0,-increment.distances)),100000)
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZATION OF ALL DATA FRAMES THAT WE NEED TO STORE addm RESULTS IN-------------------------------------------------------------
  nr_rows = length(cur.choice.dat[,trialid])*nr.reps
  len.trials = length(cur.choice.dat[,trialid])
  trialids = cur.choice.dat[,trialid]

  if (output.type == "Real" | output.type == "Fake"){
  addm.output = matrix(rep(-1,nr_rows*2),nrow=nr_rows,ncol=2)
  addm.output[,1] = -1 #Decision
  addm.output[,2] = -1 #RT

  output = rep(0,2*nr.reps)
  output.cols = c(1,2)
  nr.output.cols = length(output.cols)
  }

  if (output.type == "Full"){
    addm.output = matrix(rep(-1,nr_rows*27),nrow=nr_rows,ncol=27)
    addm.output[,1] = drift.rate
    addm.output[,2] = theta
    addm.output[,3] = cur.sd
    addm.output[,4] = non.decision.time

    # output cols are the columns in which the output (Choices and RT's) will be stored
    output = rep(0,22*nr.reps)
    output.cols = seq(5,26)
    nr.output.cols = length(output.cols)
  }

  # Valuations is an intermediate matrix that is then filled into the ev.update (for evidence update) vector that is passed to the DDM function
  cur.set_size = length(which(cur.choice.dat[1,grep("^v[1-9]*",names(cur.choice.dat)),with=FALSE] > -1))
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=len.trials)
  start.pos = which(names(cur.choice.dat) == "v1") # we define startpos, because there may be some columns in the cur.choice.dat data.table prior to the item valuation columns

  for (i in seq(cur.set_size)){
    valuations[i,] = cur.choice.dat[[start.pos + i - 1]]*drift.rate
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # DEFINE CORRCT EVIDENCE ACCUMULATION FUNCTION GIVEN INPUTS---------------------------------------------------------------------------
  if (model.type == "nomem"){
    if (output.type == "Normal" | output.type == "Fake"){
      aevacc = aevacc_hist
    } else if (output.type == "Full"){
      aevacc = aevacc_full_output
    }
  } else if (model.type == "mem"){
    if (output.type == "Normal" | output.type == "Fake"){
      aevacc = aevacc_hist
    } else if (output.type == "Full"){
      aevacc = aevacc_full_output_withmem
    }
  } else if (model.type == "memzeronoise"){
    if (output.type == "Normal" | output.type == "Fake"){
      aevacc = aevacc_hist
    } else if (output.type == "Full"){
      aevacc = aevacc_full_output_withmem_zeronoise
    }
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # WE LOOP THROUGH ALL TRIALS FOR THE addm RUNS----------------------------------------------------------------------------------------
  output.row.min = 1
  output.row.max = nr.reps
  cur.fix.nr = 0
  cnt = 1

  #fixdur.vec = cur.eye.dat$fixdur.sampled %/% timestep.ms

  if (fixation.model == "FakePath"){
    duration.length = 200
    amount.fixations = ceiling(cur.max.RT / duration.length)
    repetitions.of.sequence = ceiling(amount.fixations / cur.set_size)

    cur.fixations = rep(1:cur.set_size,repetitions.of.sequence)
    cur.durations = (rep(duration.length,amount.fixations) + (timestep.ms/2)) %/% timestep.ms
    fixdur.vec = duration.length
  }

  if (fixation.model == "Random"){
    cur.fixations = as.numeric()
    cur.durations = 0
    fixdur.vec = 0
  }

  for (trial in trialids){
        output[] = aevacc(nr.reps,
                          cur.max.RT,
                          cur.sd,
                          theta,
                          valuations[,cnt],
                          cur.fixations,
                          cur.durations,
                          fixdur.vec,
                          cur.set_size,
                          sample)

      addm.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr.reps,ncol=nr.output.cols,byrow=TRUE)
      output.row.min[1] = sum(output.row.min,nr.reps)
      output.row.max[1] = sum(output.row.max,nr.reps)
      cnt = cnt + 1
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # STORING DATA FRAME THAT COLLECTS ALL RELEVANT INFORMATION CONCERNING MODEL OUTPUT
  #-------------------------------------------------------------------------------------------------------------------------------------
  if (output.type == "Full"){
    addm.output.frame = data.table(Drift.Rate = addm.output[,1],
                                   Theta = addm.output[,2], SD = addm.output[,3],
                                   Non.decision.time = addm.output[,4],
                                   Decision = addm.output[,5], Nr.fixations = addm.output[,6],
                                   RT = addm.output[,7] * timestep.ms, Item.last.attended = addm.output[,8],
                                   Value.last.attended = addm.output[,9] / drift.rate, Chosen.last.attended = addm.output[,10],
                                   Duration.1 = addm.output[,11]*timestep.ms, Duration.2 = addm.output[,12]*timestep.ms,
                                   Duration.3 = addm.output[,13]*timestep.ms, Duration.4 = addm.output[,14]*timestep.ms,
                                   Duration.5 = addm.output[,15]*timestep.ms, Duration.6 = addm.output[,16]*timestep.ms,
                                   Duration.7 = addm.output[,17]*timestep.ms, Duration.8 = addm.output[,18]*timestep.ms,
                                   Nr.fix.1 = addm.output[,19], Nr.fix.2 = addm.output[,20],
                                   Nr.fix.3 = addm.output[,21], Nr.fix.4 = addm.output[,22],
                                   Nr.fix.5 = addm.output[,23], Nr.fix.6 = addm.output[,24],
                                   Nr.fix.7 = addm.output[,25], Nr.fix.8 = addm.output[,26],
                                   trialid = rep(trialids,each=nr.reps))
    return(addm.output.frame)
  } else if (output.type == "Real" | output.type == "Fake"){

    # Generate data.table from outputs
    addm.output.frame = data.table(decision=addm.output[,1],
                                   rt=addm.output[,2]*timestep.ms,
                                   trialid=rep(trialids,each=nr.reps))
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # FOR OUTPUT TYPE "FAKE" WE LOAD A PRECOMPUTED SIMULATED DATAFRAME AND TEST LOGLIKELIHOOD AGAINST ITS OUTCOMES
  # THIS IS IMPORTANT FOR TESTING PURPOSES, SPECIFICALLY WHETHER WE CAN RECOVER CERTAIN PARAMETER SETS
  #-------------------------------------------------------------------------------------------------------------------------------------
  if (output.type == "Fake"){

    addm.output.frame$rtbins = cut(addm.output.frame$rt,cur.breaks,include.lowest = TRUE)
    addm.choice.table = addm.output.frame %>% group_by(trialid,decision,rtbins) %>% summarise(count = n()) %>% mutate(choice.p = count / nr.reps)

    # If we run the function to generate a data frame instead of getting log likelihood we just return the frame here
    if (generate == 1){
      addm.choice.table.fake <<- addm.choice.table
      return()
    }

############################### CALCULATING LOGLIKELIHOOD ###########################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------
    setnames(addm.choice.table.fake,'choice.p','choice.p.fake')
    setkey(addm.choice.table,trialid,decision,rtbins)
    setkey(addm.choice.table.fake,trialid,decision,rtbins)

    temp.table = addm.choice.table.fake[addm.choice.table]
    temp.table[,log.lik:=log(choice.p)*choice.p.fake]

    LogLik = sum(temp.table[,log.lik])

  # Now store loglikelihood -----------------------------------------------------------------------------------------------------------
    total.log.lik = c(drift.rate,theta,cur.sd,non.decision.time,nr.reps,LogLik)
    print(total.log.lik)
    return(total.log.lik)
  }
  #-------------------------------------------------------------------------------------------------------------------------------------
}
