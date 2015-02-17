###### I AM GOING TO USE THE ADDM BASIS AND ADJUST IT TOWARDS NO ATTENTIONAL EFFECT

#   inputs needed:

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

  # INITIALIZATION OF PARAMETERS AND RELEVANT SUBSETS OF DATA FRAMES
  #-------------------------------------------------------------------------------------------------------------------------------------
  # initialize model parameters
  drift.rate = core.parameters[1]
  theta = core.parameters[2]
  cur.sd = core.parameters[3]
  non.decision.time = core.parameters[4]
  #-------------------------------------------------------------------------------------------------------------------------------------

  # SOME MISCELLANEOUS VARIABLES THAT ARE UTILIZED LATER
  #-------------------------------------------------------------------------------------------------------------------------------------
  # First we need to define max. RT
  cur.max.RT = 40000 %/% timestep.ms

  # generating column that provides a cutted (binned) version of reaction times
  increment.distances = 100
  cur.breaks = c(rev(seq(20000,0,-increment.distances)),100000)
  #-------------------------------------------------------------------------------------------------------------------------------------
  nr_rows = length(cur.choice.dat[,trialid])*nr.reps
  len.trials = length(cur.choice.dat[,trialid])
  trialids = cur.choice.dat[,trialid]

  # INITIALIZATION OF ALL DATA FRAMES THAT WE NEED TO STORE ADDM RESULTS IN
  #-------------------------------------------------------------------------------------------------------------------------------------
  if (output.type == "Real" | output.type == "Fake"){
  aDDM.output = matrix(rep(-1,nr_rows*2),nrow=nr_rows,ncol=2)
  aDDM.output[,1] = -1 #Decision
  aDDM.output[,2] = -1 #RT

  output = rep(0,2*nr.reps)
  output.cols = c(1,2)
  nr.output.cols = length(output.cols)
  }

  if (output.type == "Full"){
    aDDM.output = matrix(rep(-1,nr_rows*27),nrow=nr_rows,ncol=27)
    aDDM.output[,1] = drift.rate
    aDDM.output[,2] = theta
    aDDM.output[,3] = cur.sd
    aDDM.output[,4] = non.decision.time

    # output cols are just the columns in which the output (Choices and RT's) will be stored
    output = rep(0,22*nr.reps)
    output.cols = seq(5,26)
    nr.output.cols = length(output.cols)
  }

  # Valuations is an intermediate matrix that is then filled into the ev.update (for evidence update) vector that is passed to the DDM function
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=len.trials)
  cur.set_size = length(which(cur.choice.dat[1,grep("^v[1-9]",names(cur.choice.dat)),with=FALSE] > -1))
  start.pos = which(names(cur.choice.dat) == "v1") # we define startpos, because there may be some columns in the cur.choice.dat data.table prior to the item valuation columns
  for (i in seq(cur.set_size)){
    valuations[i,] = cur.choice.dat[[start.pos + i - 1]]*drift.rate
  }

  # DEFINE CORRCT EVIDENCE ACCUMULATION FUNCTION GIVEN INPUTS
  #-------------------------------------------------------------------------------------------------------------------------------------
  if (model.type == "nomem"){
    if (output.type == "Normal"){
      aevacc = aevacc_hist
    } else if (output.type == "Full"){
      aevacc = aevacc_full_output
    }
  } else if (model.type == "mem"){
    if (output.type == "Normal"){
      aevacc = aevacc_hist
    } else if (output.type == "Full"){
      aevacc = aevacc_full_output_withmem
    }
  } else if (model.type == "memzeronoise"){
    if (output.type == "Normal"){
      aevacc = aevacc_hist
    } else if (output.type == "Full"){
      aevacc = aevacc_full_output_withmem_zeronoise
    }
  }


  #-------------------------------------------------------------------------------------------------------------------------------------

  # WE LOOP THROUGH ALL TRIALS FOR THE ADDM RUNS
  #-------------------------------------------------------------------------------------------------------------------------------------
  output.row.min = 1
  output.row.max = nr.reps

  cur.fix.nr = 0
  cnt = 1
  #last.fix.add = 20000 %/% timestep.ms
  fixdur.vec = eye$fixdur.sampled %/% timestep.ms

  if (fixation.model == "FakePath"){
    duration.length = 200
    amount.fixations = ceiling(cur.max.RT / duration.length)
    repetitions.of.sequence = ceiling(amount.fixation / cur.set_size)

    cur.fixations = rep(1:cur.set_size,repetitions.of.sequence)
    cur.durations = (rep(duration.length,repetitions.of.sequence*cur.set_size) + timestep.ms/2) %/% timestep.ms
  }

  if (fixation.model == "Random"){
    cur.fixations = as.numeric()
    cur.durations = 0
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

      aDDM.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr.reps,ncol=nr.output.cols,byrow=TRUE)
      output.row.min[1] = sum(output.row.min,nr.reps)
      output.row.max[1] = sum(output.row.max,nr.reps)
      cnt = cnt + 1
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # STORING DATA FRAME THAT COLLECTS ALL RELEVANT INFORMATION CONCERNING MODEL OUTPUT
  #-------------------------------------------------------------------------------------------------------------------------------------
  if (output.type == "Full"){
    aDDM.output.frame = data.table(Drift.Rate = aDDM.output[,1],
                                           Theta = aDDM.output[,2], SD = aDDM.output[,3],
                                           Non.decision.time = aDDM.output[,4],
                                           Decision = aDDM.output[,5], Nr.fixations = aDDM.output[,6],
                                           RT = aDDM.output[,7] * timestep.ms, Item.last.attended = aDDM.output[,8],
                                           Value.last.attended = aDDM.output[,9] / drift.rate, Chosen.last.attended = aDDM.output[,10],
                                           Duration.1 = aDDM.output[,11]*timestep.ms, Duration.2 = aDDM.output[,12]*timestep.ms,
                                           Duration.3 = aDDM.output[,13]*timestep.ms, Duration.4 = aDDM.output[,14]*timestep.ms,
                                           Duration.5 = aDDM.output[,15]*timestep.ms, Duration.6 = aDDM.output[,16]*timestep.ms,
                                           Duration.7 = aDDM.output[,17]*timestep.ms, Duration.8 = aDDM.output[,18]*timestep.ms,
                                           Nr.fix.1 = aDDM.output[,19], Nr.fix.2 = aDDM.output[,20],
                                           Nr.fix.3 = aDDM.output[,21], Nr.fix.4 = aDDM.output[,22],
                                           Nr.fix.5 = aDDM.output[,23], Nr.fix.6 = aDDM.output[,24],
                                           Nr.fix.7 = aDDM.output[,25], Nr.fix.8 = aDDM.output[,26],
                                           trialid = rep(trialids,each=nr.reps))
    return(aDDM.output.frame)
  } else if (output.type == "Real" | output.type == "Fake"){

    # Generate data.table from outputs
    aDDM.output.frame = data.table(decision=aDDM.output[,1],
                                   rt=aDDM.output[,2],
                                   trialid=rep(trialids,each=nr.reps))
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # FOR OUTPUT TYPE "FAKE" WE LOAD A PRECOMPUTED SIMULATED DATAFRAME AND TEST LOGLIKELIHOOD AGAINST ITS OUTCOMES
  # THIS IS IMPORTANT FOR TESTING PURPOSES, SPECIFICALLY WHETHER WE CAN RECOVER CERTAIN PARAMETER SETS
  #-------------------------------------------------------------------------------------------------------------------------------------
  if (output.type == "Fake"){

    aDDM.output.frame$rtbins = cut(aDDM.output.frame$rt,cur.breaks,include.lowest = TRUE)
    aDDM.choice.table = aDDM.output.frame %>% group_by(trialid,decision,rtbins) %>% summarise(count = n()) %>% mutate(choice.p = count / nr.reps)

    # If we run the function to generate a data frame instead of getting log likelihood we just return the frame here
    if (generate == 1){
      return(aDDM.choice.table)
    }

    if (cur.set_size == 4){
      aDDM.choice.table.fake = core.clean.artificial.choice.table.4
    }

    if (cur.set_size == 6){
      aDDM.choice.table.fake = core.clean.artificial.choice.table.6
    }

    if (cur.set_size == 8){
      aDDM.choice.table.fake = core.clean.artificial.choice.table.8
    }

############################### CALCULATING LOGLIKELIHOOD ###########################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------
    setnames(aDDM.choice.table.fake,'choice.p','choice.p.fake')
    setkey(aDDM.choice.table.fake,Trialid,rtbins)
    setkey(aDDM.choice.table.fake,Trialid,rtbins)

    temp.table = aDDM.choice.table.fake[aDDM.choice.table]
    temp.table[,log.lik:=log(choice.p)*choice.p.fake]

    LogLik = sum(temp.table[,log.lik])

#     vec.pos.cnt = 1
#     aDDM.choice.likelihood = rep(0,length(aDDM.choice.table.fake[,1]))
#
#     for (trial in trialids){
#
#       choices = unique(aDDM.choice.table.fake$decision[aDDM.choice.table.fake$trialid == trial])
#
#       for (choice in choices){
#
#         cur.dat.fake = aDDM.choice.table.fake[aDDM.choice.table.fake$trialid == trial & aDDM.choice.table.fake$decision == choice,]
#         cur.dat.sim = aDDM.choice.table[aDDM.choice.table$trialid == trial & aDDM.choice.table$decision == choice,]
#
#         if (length(cur.dat.sim[,1]) >= 1){
#
#           cur.rtbins = unique(cur.dat.fake$rtbins)
#
#           for (rt.bin in cur.rtbins){
#
#             fake.p = as.vector(cur.dat.fake$choice.p[cur.dat.fake$rtbins == rt.bin])
#             sim.p = as.vector(cur.dat.sim$choice.p[cur.dat.sim$rtbins == rt.bin])
#             cur.length = length(sim.p)
#
#             if (cur.length >= 1){
#               aDDM.choice.likelihood[vec.pos.cnt] = fake.p[1]*log(sim.p[1])
#             }
#
#             if (cur.length < 1){
#               aDDM.choice.likelihood[vec.pos.cnt] = fake.p[1]*log(0.0000001)
#             }
#
#             vec.pos.cnt = vec.pos.cnt + 1
#           }
#         }
#       }
#     }

  # Now store loglikelihood -----------------------------------------------------------------------------------------------------------
    total.log.lik = c(drift.rate,theta,cur.sd,non.decision.time,nr.reps,LogLik)
    print(total.log.lik)
    return(total.log.lik)
  }
  #-------------------------------------------------------------------------------------------------------------------------------------
}
