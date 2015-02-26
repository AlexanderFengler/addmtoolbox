#' Runs the model by unique trial conditions.
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Model run by condition
#' @return The function has three potential return values. A log likelihood value, utilized when trying to fit the model. A simple data.table providing simulated rts,decisions by condition id, which useful for generating fake data in testing. A full model output with many details, utilized for running the model with optimal parameters and extracting data for plots.
#' \code{addm_by_condition} Returns either log likelihoods, a simple id,decision,rt data frame or a detailed model output
#' @export
#' @param conditions.dat A 'data.frame' storing the item valuations (v1,v2...) by unique trial conditions. An id column (id) needs to be provided that matches by trial data.
#' @param eye.dat A 'data.frame' or 'data.table' storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and an id column (id). Can all be initialized as zero columns when a by condition fit is attempted.
#' @param choice.dat A 'data.frame' or  'data.table' storing the item valuations (v1,v2...) , reaction times in ms (rt), decisions as decision and an id column (id). A by trial form is assumed.
#' @param model.parameters A vector with the four core addm parameters in order (drift.rate,theta,sd,non.decision.time).
#' @param nr.reps An integer number that tells the function how many simulation runs to use.
#' @param model.type A string that indicates which version of the model to run. 'standard' or 'memnoise' when memory effects shall be allowed.
#' @param output.type A string that indicates what output the model shall produce. 'full' for detailed model output, 'fit' for sparse output (rt,decision) by id variable.
#' @param fixation.model A string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (implemented) 'fixedpath' for following a predetermined fixation path with fixed durations (implemented) or 'user' for a user supplied fixation model (function name: user_fixation_model)
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param generate Binary variable that tells the function to return either log likelihood values (0) or rt, decision (1). Relevant only if model.type variable is 'fit'.

addm_by_condition = function(conditions.dat = data.table(v1 = 0, v2 = 0, id = 0),
                             eye.dat = data.table(fixloc = 0, fixdur = 0, fixnr = 1, id = 0),
                             choice.dat = data.table(decision = 0, rt = 0, id = 0),
                             model.parameters = c(0.006,0.6,0.06,0),
                             nr.reps = 2000,
                             model.type = 'standard',
                             output.type = 'fit',
                             fixation.model = 'fixedpath',
                             timestep = 10,
                             generate = 0){

  # INITIALIZATION OF PARAMETERS AND RELEVANT SUBSETS OF DATA FRAMES--------------------------------------------------------------------
  # Initialize model parameters
  drift.rate = model.parameters[1]
  theta = model.parameters[2]
  cur.sd = model.parameters[3]
  non.decision.time = model.parameters[4]
  #-------------------------------------------------------------------------------------------------------------------------------------

  # SOME MISCELLANEOUS VARIABLES THAT ARE UTILIZED LATER--------------------------------------------------------------------------------

  # First we need to define max.RT
  cur.max.RT = 40000

  # Generating column that provides a cutted (binned) version of reaction times
  increment.distances = 100
  cur.breaks = c(rev(seq(20000,0,-increment.distances)),100000)
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZATION OF ALL DATA FRAMES THAT WE NEED TO STORE addm RESULTS IN-------------------------------------------------------------
  nr_rows = length(conditions.dat[,id])*nr.reps
  len.trials = length(conditions.dat[,id])
  ids = conditions.dat[,id]
  cur.set_size = length(which(conditions.dat[1,grep("^v[1-9]*",names(conditions.dat)),with=FALSE] > -1))

  # Define Matrix that stores values adjusted for drift rates // will be fed into evidence accumulaiton function
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=len.trials)

  # Define position of vector where valuations start
  start.pos = which(names(conditions.dat) == "v1")

  for (i in seq(cur.set_size)){
    valuations[i,] = conditions.dat[[start.pos + i - 1]]
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
      aevacc = aevacc_full_output_memnoise
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
  if (fixation.model == "fixedpath"){
    fixation_model = addm2_fixation_model_fixedpath
  }

  if (fixation.model == "random"){
    fixation_model = addm2_fixation_model_random
  }

  if (fixation.model == "user"){
    fixation_model = user_fixation_model
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
                      timestep,
                      nr.reps,
                      cur.max.RT,
                      valuations[,cnt],
                      fixation_model)

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
                                   rt=addm.output[,2]*timestep,
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
