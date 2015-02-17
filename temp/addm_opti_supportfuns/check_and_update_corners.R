# This function takes in a matrix that provides loglikelihoods as a result of parameter combinations 
# and checks whether the current optimal solution is a corner solution in the parameter space

# If it is a corner solution we update the parameter space and return a new version
# If it is not a corner solution we do nothing an let the main scripts continue with fine grid search


# Log Likelihood file is ordered by columns as follows
# 1 - Subject
# 2 - Set Size
# 3 - Drift Rate
# 4 - theta
# 5 - SD
# 6 - Non.decision.time
# 7 - Repetitions of Simulation
# 8 - Log.Likelihood

check.and.update.corners = function(log.liks,
                                    drift.rates,
                                    thetas,
                                    sds,
                                    non.decision.times,
                                    drift.shift,
                                    drift.step,
                                    sd.shift,
                                    sd.step){
  
  simple.log.lik = round(as.matrix(log.liks[,c(3,4,5,8),with=FALSE]),5)
  
  unique.drifts = unique(simple.log.lik[,1])
  max.drift = max(unique.drifts)
  min.drift = min(unique.drifts)
  
  unique.thetas = unique(simple.log.lik[,2])
  max.theta = max(unique.thetas)
  min.theta = min(unique.thetas)
  
  unique.sds = unique(simple.log.lik[,3])
  max.sd = max(unique.sds)
  min.sd = min(unique.sds)
  
  # in case there are multiple minima, we take a a random sample, moving further with one of them # THIS MAY BE CHANGED TO SEARCH AROUND MULTIPLE MINIMA !!!
  max.val = max(simple.log.lik[,4])
  max.pos = which(simple.log.lik[,4] == max.val)[sample(length(which(simple.log.lik[,4] == max.val)),1)]
  
  # We initialize the "new" parameters that we will pass in the function to generate a parameter matrix later
  
  new.drifts = drift.rates
  new.thetas = thetas
  new.sds = sds
  new.non.decision.times = non.decision.times
  
  corner.cnt = 0
  
  # Test for each corner (we do not need to test theta)
  # ---------------------------------------------------------------------------------------------
  
  if (simple.log.lik[max.pos,1] == max.drift){
    corner.cnt = corner.cnt + 1
    new.drifts = seq(max.drift,max.drift+drift.shift,drift.step)
    new.drifts = new.drifts[new.drifts != max.drift]
  }
  
  if (simple.log.lik[max.pos,1] == min.drift){
    corner.cnt = corner.cnt + 1
    new.drifts = seq(min.drift,max(c(min.drift-drift.shift, 0)), -drift.step)
    new.drifts = new.drifts[new.drifts != min.drift]
  }
  
  if (simple.log.lik[max.pos,3] == max.sd){
    corner.cnt = corner.cnt + 1
    new.sds = seq(max.sd,max.sd+sd.shift,sd.step)
    new.sds = new.sds[new.sds != max.sd]
  }
  
  if (simple.log.lik[max.pos,3] == min.sd){
    corner.cnt = corner.cnt + 1
    new.sds = seq(min.sd,max(c(min.sd - sd.shift,0)), -sd.step)
    new.sds = new.sds[new.sds != min.sd]
  }
  
  # ---------------------------------------------------------------------------------------------
  ###############################################################################################
  
  
  # Cases in which we return without parameter additions
  # ---------------------------------------------------------------------------------------------
  # We need to cover three special cases in which we will not extend the parameter space

  # If the minimum is at drift-rate 0 AND at sd 0
  if (simple.log.lik[max.pos,1] == 0 & simple.log.lik[max.pos,3] == 0){
    return(x = list(Corner = 0, New.parameters = 0))
  }
  # If the minimum is at drift-rate 0 AND at SD < max.sd
  if (simple.log.lik[max.pos,1] == 0 & simple.log.lik[max.pos,3] < max.sd){
    return(x = list(Corner = 0, New.parameters = 0))
  }
  # If the minimum is at drift-rate < max drift-rate AND SD = 0
  if (simple.log.lik[max.pos,1] < max.drift & simple.log.lik[max.pos,3] == 0){
    return(x = list(Corner = 0, New.parameters = 0))   
  }
  
  # If we had no corner solution we return values that will let us pass through to fine grid search
  if (corner.cnt == 0){
      return(x = list(Corner = 0, New.parameters = 0))
  }
  # ---------------------------------------------------------------------------------------------
  ###############################################################################################
  

  # Cases in which we need to return new parameters
  # ---------------------------------------------------------------------------------------------
  # If we touched corner with optimal loglikelihood so far we update the parameter space 
  if (corner.cnt == 1){  
    new.params = param.combs(new.drifts,new.sds,new.thetas,new.non.decision.times)
    return(x = list(Corner = 1, New.parameters = new.params))
  }
  
  # If we touched two corners we need too specify the new parameters in three parts
  # First, all thetas and sds for new.drifts
  # Second, all thetas and drifts for new sds
  # Third, all thetas and new.drifts and new.sds
  # then concatenate
  
  if (corner.cnt > 1){
    
    param.subset1 = param.combs(new.drifts,sds,thetas,non.decision.times)
    param.subset2 = param.combs(drift.rates,new.sds,thetas,non.decision.times)
    param.subset3 = param.combs(new.drifts,new.sds,thetas,non.decision.times)
    new.params = rbind(param.subset1,param.subset2,param.subset3)
    return(x = list(Corner = 1, New.parameters = new.params))
  }
  # ---------------------------------------------------------------------------------------------
  ############################################################################################### 
}