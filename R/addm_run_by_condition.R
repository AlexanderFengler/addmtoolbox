#' Runs model by unique trial condition for one set of parameter values
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Run model by unique trial condition
#' @return The function has three potential return values. A log likelihood value, utilized when trying to fit the model (output.type = 'fit', generate = 0). A simple data.table providing simulated rts,decisions by condition id, which useful for generating fake data in testing (output.type = 'fit', generate = 1). A full model output with many details, utilized for running the model with optimal parameters and extracting data for plots (output.type = 'full').
#' \code{addm_run_by_condition}
#' @export
#' @param conditions.dat  data.table storing the item valuations (v1,v2...) by unique trial conditions. An id column (conditions_id) needs to be provided
#' @param choice.dat data.table storing the item valuations (v1,v2...), reaction times in ms (rt), decisions as decision and an id column (conditions_id).
#' @param model.parameters vector with the four core addm parameters in order (non.decision.time, drift, sd, theta, gamma, boundary-parameters).
#' @param nr.attributes integer providing the amount of attributes we consider per item
#' @param nr.reps integer that tells the function how many simulation runs to use.
#' @param model.type string that indicates which version of the model to run. 'standard' for standard (a)ddm model or 'memnoise' when memory effects shall be allowed.
#' @param output.type string that indicates what output the model shall produce. 'full' for detailed model output, 'fit' for sparse output (rt,decision) by id variable.
#' @param fixation.model string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (supplied) 'fixedpath' for following a predetermined fixation path with fixed durations (supplied) or 'user' for a user supplied fixation model (function name: user_fixation_model).
#' @param timestep integer that provides the timestep-size that is used in the simulations (in ms).
#' @param generate boolean variable that tells the function to return either log likelihood values (0) or rt, decision (1). Relevant only if model.type variable is 'fit'.

addm_run_by_condition = function(choice.dat = data.table(decision = 0, rt = 0, condition_id = 0),
                                 conditions.dat = data.table(v1 = 1, v2 = 100, condition_id = 0),
                                 model.parameters = c(0,0.002,0.07,0.5),
                                 nr.attributes = 1,
                                 boundaryfun = 1,
                                 nr.reps = 2000,
                                 timestep = 10,
                                 model.type = 'standard',
                                 output.type = 'fit',
                                 fixation.model = 'fixedpath',
                                 generate = 0){

  # Compute Boundaries -----------------------------------------------------------------------------------------------------------------
  if (class(boundaryfun) == 'function'){
    # COMPUTE
  }
  # ------------------------------------------------------------------------------------------------------------------------------------

  # Initialize parameter theta ---------------------------------------------------------------------------------------------------------
  # only for readability, because it is used later on
  theta = model.parameters[4]

  # fill up model parameter vector if supplied to little
  if (length(model.parameters) < 7){
    model.parameters[(length(model.parameters)+1):7] = 0
  }
  #return(model.parameters)
  # ------------------------------------------------------------------------------------------------------------------------------------

  # SOME MISCELLANEOUS VARIABLES THAT ARE UTILIZED LATER--------------------------------------------------------------------------------
  # First we need to define max.RT
  cur.max.RT = 40000

  # Generating column that provides a cutted (binned) version of reaction times
  increment.distances = 100
  cur.breaks = c(rev(seq(20000,0,-increment.distances)),100000)
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZATION OF ALL DATA FRAMES THAT WE NEED TO STORE addm RESULTS IN-------------------------------------------------------------
  nr_rows = length(conditions.dat[,condition_id])*nr.reps
  len.trials = length(conditions.dat[,condition_id])
  ids = conditions.dat[,condition_id]
  cur.set_size = length(which(conditions.dat[1, grep("^v[1-9]*", names(conditions.dat)), with=FALSE] > -1))

  # Define Matrix that stores values by condition
  valuations=matrix(rep(0,len.trials*cur.set_size),nrow=cur.set_size,ncol=len.trials)

  # Define position of vector where valuations start
  start.pos =  grep("^v1",names(conditions.dat))[1]

  for (i in seq_len(cur.set_size)){
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
    len.params = length(model.parameters)
    addm.output = matrix(rep(-1,nr_rows*(len.params + 6 +(2*cur.set_size))),nrow = nr_rows, ncol = len.params + 6 +(2*cur.set_size))

    for (i in 1:len.params){
    addm.output[,i] = model.parameters[i]
    }

    # output cols are the columns in which the output (Choices and RT's) will be stored
    output = rep(0,(6+(2*cur.set_size))*nr.reps)
    output.cols = seq(len.params + 1, len.params + 6 +(2*cur.set_size))
    nr.output.cols = length(output.cols)
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # DEFINE CORRCT EVIDENCE ACCUMULATION FUNCTION GIVEN INPUTS---------------------------------------------------------------------------
  if (model.type == "standard"){
    if (output.type == "fit"){
      if (nr.attributes == 1){
        if (cur.set_size / nr.attributes == 2){
          if (theta == 1){
            aevacc = evacc2_by_condition
          } else {
            aevacc = aevacc2_by_condition
          }
        } else {
          aevacc = aevacc_by_condition
        }
      } else {
        if (cur.set_size / nr.attributes == 2){
          aevacc = aevaccma2_by_condition
        } else {
          stop('You attempted to run a set of > 2 items which have multiple attributes each (by condition): This is not yet implemented!')
        }
      }
    } else if (output.type == "full"){
      if (cur.set_size / nr.attributes == 2){
        if (nr.attributes == 1){
          aevacc = aevacc2_full_output
        } else {
          aevacc = aevaccma2_full_output
        }
      } else {
        if (nr.attributes == 1){
          aevacc = aevacc_full_output_memnoise
        } else {
          stop('You attempted to get detailed model output for a multiattributes version of the model: This is not yet implemented!')
        }
      }
    }
  } else if (model.type == "memnoise"){
    if (output.type == "fit"){
      if (nr.attributes == 1){
      aevacc = aevacc_by_condition_memnoise
      } else {
        stop('You attempted to fit a model with memory effects as well as multiple attributes per item: This is not yet implemented!')
      }
    } else if (output.type == "full"){
      if (nr.attributes == 1){
      aevacc = aevacc_full_output_memnoise
      } else {
        stop('You attempted to fit a model with memory effect as well as mutliple attributes per item: This is not yet implemented!')
      }
    }
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # INITIALIZE FIXATION PATHWAYS -------------------------------------------------------------------------------------------------------
  if (fixation.model == "fixedpath"){
    if (cur.set_size / nr.attributes == 2){
      if (nr.attributes == 1){
        fixation_model = addm2_fixation_model_fixedpath
      } else {
        fixation_model = addm_fixation_model_fixedpath
      }
    } else {
      fixation_model = addm_fixation_model_fixedpath
    }
  }

  if (fixation.model == "random"){
    if (cur.set_size / nr.attributes == 2){
      if (nr.attributes == 1){
        fixation_model = addm2_fixation_model_random
      } else {
        fixation_model = addm_fixation_model_random
      }
    } else {
      fixation_model = addm_fixation_model_random
    }
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
    output[] = aevacc(model.parameters,
                      cur.max.RT,
                      valuations[,cnt],
                      nr.attributes,
                      fixation_model,
                      nr.reps,
                      timestep)
    addm.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr.reps,ncol=nr.output.cols,byrow=TRUE)
    output.row.min[1] = sum(output.row.min,nr.reps)
    output.row.max[1] = sum(output.row.max,nr.reps)
    cnt = cnt + 1
  }
  #-------------------------------------------------------------------------------------------------------------------------------------

  # STORING DATA FRAME THAT COLLECTS ALL RELEVANT INFORMATION CONCERNING MODEL OUTPUT---------------------------------------------------
  if (output.type == "full"){
      output.names = c("condition_id",
                       "non.decision.time",
                       "drift.rate",
                       "sd",
                       "theta",
                       "gamma",
                       "scalar_items_not_seen_drift",
                       "scalar_items_not_seen_noise")

      namelen = length(output.names)
      if (length(model.parameters) > 7){
        for (i in 8:length(model.parameters)){
          output.names[namelen + (i-7)] = paste("boundary.param.",toString(i-7),sep='')
        }
      }

      namelen = length(output.names)
      output.names[(namelen+1):(namelen + 6)] = c('decision', 'nr.fixations','rt','item.last.attended','value.last.attended','chosen.last.attended')

      namelen = length(output.names)
      for (i in seq_len(cur.set_size)){
        output.names[namelen + i] = paste("duration.",toString(i),sep='')
        output.names[namelen + cur.set_size + i] = paste("nr.fix.",toString(i),sep='')
      }

    addm.output.frame = cbind(data.table(condition_id = rep(ids,each=nr.reps)),as.data.table(addm.output))
    setnames(addm.output.frame,output.names)
    return(addm.output.frame)

  } else if (output.type == "fit"){

    addm.output.frame = data.table(decision=addm.output[,1],
                                   rt=addm.output[,2],
                                   condition_id=rep(ids,each=nr.reps))

    if (generate == 1){
      return(addm.output.frame)
    }

    addm.output.frame$rtbins = cut(addm.output.frame$rt,cur.breaks,include.lowest = TRUE)

    addm.choice.table = addm.output.frame %>%
      group_by(condition_id,decision,rtbins) %>%
      summarise(count = n()) %>%
      mutate(choice.p = count / nr.reps)

    choice.dat$rtbins = cut(choice.dat$rt,cur.breaks,include.lowest = TRUE)
    real.choice.table = choice.dat %>% select(condition_id,decision,rtbins)
    #-------------------------------------------------------------------------------------------------------------------------------------

    # CALCULATE LOGLIKELIHOOD  ----------------------------------------------------------------------------------------------------------
    setkey(addm.choice.table,condition_id,decision,rtbins)
    setkey(real.choice.table,condition_id,decision,rtbins)

    temp.table = addm.choice.table[real.choice.table]
    temp.table[is.na(choice.p),choice.p:=1/(nr.reps + 1)] # alternative would be to maybe use 1/(nr.reps + 1)
    temp.table[,log.lik:=log(choice.p)]
    LogLik = (-1)*sum(temp.table[,log.lik])
    # -----------------------------------------------------------------------------------------------------------------------------------

    # STORE AND RETURN ------------------------------------------------------------------------------------------------------------------
    if (nr.attributes == 1){
      total.log.lik = c(model.parameters,LogLik)
    } else {
      total.log.lik  = c(model.parameters,LogLik)
    }
    print(total.log.lik)
    return(total.log.lik)
    # -----------------------------------------------------------------------------------------------------------------------------------
  }
  #--------------------------------------------------------------------------------------------------------------------------------------
}
