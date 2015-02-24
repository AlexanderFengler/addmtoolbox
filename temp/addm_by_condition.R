#   Inputs needed:
# - subject
# - set size
# - drift.rate
# - theta
# - standard deviation
# - non decision time

# - number of repetitions used -- nr.reps
# - fixation.model -- 'normal'
#                  -- 'fakepath' (used for fitting pregenerated datasets)
#                  -- 'random'
# - output.type -- 'full'
#               -- 'fit'
# - model.type -- "standard" use addm function with memory effects (whichever the best at the moment)
#               -- "memnoise" use addm function without memory (normal Krajbich Type)
# - generate -- boolean variable indicating whether we run the function to generate a data frame (1) or to get a log likelihood (0)

addm_by_condition = function(conditions,
                             cur.eye.dat,
                             choice.dat,
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
  nr_rows = length(conditions[,id])*nr.reps
  len.trials = length(conditions[,id])
  ids = conditions[,id]
  cur.set_size = length(which(conditions[1,grep("^v[1-9]*",names(conditions)),with=FALSE] > -1))

  # Define Matrix that stores values adjusted for drift rates // will be fed into evidence accumulaiton function
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=len.trials)

  # Define position of vector where valuations start
  start.pos = which(names(conditions) == "v1")

  for (i in seq(cur.set_size)){
    valuations[i,] = conditions[[start.pos + i - 1]]
  }

  if (output.type == "fit"){
    addm.output = matrix(rep(-1,nr_rows*2),nrow=nr_rows,ncol=2)
    addm.output[,1] = -1 #Decision
    addm.output[,2] = -1 #RT

    output = rep(0,2*nr.reps)
    output.cols = c(1,2)
    nr.output.cols = length(output.cols)
  }

  if (output.type == "full"){
    addm.output = matrix(rep(-1,nr_rows*(5 + 6 +(2*cur.set_size)),nrow=nr_rows,ncol=5 + 6 +(2*cur.set_size)))

    addm.output[,1] = drift.rate
    addm.output[,2] = theta
    addm.output[,3] = cur.sd
    addm.output[,4] = non.decision.time

    # output cols are the columns in which the output (Choices and RT's) will be stored
    output = rep(0,(6+(2*cur.set_size))*nr.reps)
    output.cols = seq(5,5 + 6 +(2*cur.set_size))
    nr.output.cols = length(output.cols)
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # DEFINE CORRCT EVIDENCE ACCUMULATION FUNCTION GIVEN INPUTS---------------------------------------------------------------------------
  if (model.type == "standard"){
    if (output.type == "fit"){
      if (cur.set_size ==2){
        aevacc = aevacc2_by_condition
        } else {
          aevacc = aevacc_by_condition
      }
    } else if (output.type == "full"){
      aevacc = aevacc_full_output
    }
  } else if (model.type == "memnoise"){
    if (output.type == "fit"){
      aevacc = aevacc_by_condition_memnoise
    } else if (output.type == "full"){
      aevacc = aevacc_full_output_memnoise
    }
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZE FIXATION PATHWAYS -------------------------------------------------------------------------------------------------------
  if (fixation.model == "fakepath"){
    duration.length = 400
    amount.fixations = ceiling(cur.max.RT / duration.length)
    repetitions.of.sequence = ceiling(amount.fixations / cur.set_size)

    cur.fixations = rep(1:cur.set_size, repetitions.of.sequence)
    cur.durations = rep(duration.length, amount.fixations)
    fixdur.vec = duration.length
  }

  if (fixation.model == "random"){
    cur.fixations = as.numeric()
    cur.durations = 0
    fixdur.vec = 0
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # WE LOOP THROUGH ALL TRIALS FOR THE addm RUNS----------------------------------------------------------------------------------------
  output.row.min = 1
  output.row.max = nr.reps
  cur.fix.nr = 0
  cnt = 1

  for (id in ids){
    output[] = aevacc(cur.sd,
                      theta,
                      drift.rate,
                      non.decision.time,
                      timestep.ms,
                      nr.reps,
                      cur.max.RT,
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

  # STORING DATA FRAME THAT COLLECTS ALL RELEVANT INFORMATION CONCERNING MODEL OUTPUT---------------------------------------------------
  if (output.type == "full"){
    output.names = c("id",
                     "drift.rate",
                     "theta",
                     "sd",
                     "non.decision.time",
                     "decision",
                     "nr.fixations",
                     "rt",
                     "item.last.attended",
                     "value.last.attended",
                     "chosen.last.attended")

    namelen = length(output.names)
    for (i in seq_len(cur.set_size)){
      output.names[namelen + i] = paste("duration.",toString(i),sep='')
      output.names[namelen + cur.set_size + i] = paste("nr.fix.",toString(i),sep='')
    }

    addm.output.frame = cbind(as.data.table(addm.output),as.data.table(id = rep(ids,each=nr.reps)))
    setnames(addm.output.frame,output.names)

    return(addm.output.frame)
  } else if (output.type == "fit"){

    addm.output.frame = data.table(decision=addm.output[,1],
                                   rt=addm.output[,2]*timestep.ms,
                                   id=rep(ids,each=nr.reps))

    if (generate == 1){
      return(addm.output.frame)
    }

    addm.output.frame$rtbins = cut(addm.output.frame$rt,cur.breaks,include.lowest = TRUE)

    addm.choice.table = addm.output.frame %>%
      group_by(id,decision,rtbins) %>%
      summarise(count = n()) %>%
      mutate(choice.p = count / nr.reps)

    choice.dat$rtbins = cut(choice.dat$rt,cur.breaks,include.lowest = TRUE)
    real.choice.table = choice.dat %>% select(id,decision,rtbins)

    # CALCULATE LOGLIKELIHOOD  ----------------------------------------------------------------------------------------------------------
    setkey(addm.choice.table,id,decision,rtbins)
    setkey(real.choice.table,id,decision,rtbins)

    temp.table = addm.choice.table[real.choice.table]
    temp.table[is.na(choice.p),choice.p:=1/(nr.reps + 1)]
    temp.table[,log.lik:=log(choice.p)]
    LogLik = sum(temp.table[,log.lik])
    # -----------------------------------------------------------------------------------------------------------------------------------

    # STORE AND RETURN ------------------------------------------------------------------------------------------------------------------
    total.log.lik = c(drift.rate,theta,cur.sd,non.decision.time,nr.reps,LogLik)
    print(total.log.lik)
    return(total.log.lik)
    # -----------------------------------------------------------------------------------------------------------------------------------
  }
  #--------------------------------------------------------------------------------------------------------------------------------------
}
