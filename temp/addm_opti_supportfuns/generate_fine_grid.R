# This function takes in the optimal solution of parameters according to the coarse grid search
# and generates a parameter matrix that will serve as the basis for the fine grid search

# Specifically we are running a space fine grid that is centered around the hitherto optimal parameters, that covers
# in shorter intervall the space between the optimal parameters and the next by upper and lower parameters tested in the 
# coarse grid search (example: we finde drift.rate of 0.002 optimal and before tested 0.001,0.002,0.003....now we test 
# 0.001,0.00125,0.0015,0.00175,0.002,0.00225,0.0025,0.00275,0.003...analogous for all other variables)

generate.fine.grid = function(drift.step.fine,
                              sd.step.fine,
                              thetas,
                              non.decision.time.step.fine,
                              coarse.to.fine.ratio,
                              log.liks){
  
  # find max 
  # If multiple solutions currently sampling only one THIS SHOULD PROBABLY BE CHANGED into searching through multiple minima
  
  max.val = max(log.liks$Log.Likelihood)
  max.pos = which(log.liks$Log.Likelihood == max.val)[sample(length(which(log.liks$Log.Likelihood == max.val)),1)]
  
  
  center.parameters = as.numeric(log.liks[max.pos,])
  
  center.drift = center.parameters[3]  
  #center.theta = center.parameters[4] #we actually do not fit a fine thate -- we keep theta spanning the whole range
  center.sd = center.parameters[5]
  center.non.decision.time = center.parameters[6]
  
  subject = as.numeric(log.liks[1,1,with=FALSE])
  set.size = as.numeric(log.liks[1,2,with=FALSE])
  
  fine.non.decision.times = as.numeric(log.liks[1,6,with=FALSE])
    
    #seq(center.non.decision.time-coarse.to.fine.ratio*non.decision.time.step.fine,
     #                   center.non.decision.time+coarse.to.fine.ratio*non.decision.time.step.fine,
      #                  non.decision.time.step.fine)
  
  fine.drifts = seq(center.drift-(coarse.to.fine.ratio-1)*drift.step.fine,
                    center.drift+(coarse.to.fine.ratio-1)*drift.step.fine,
                    drift.step.fine)
  
  fine.thetas = thetas
  
  fine.sds = seq(center.sd-(coarse.to.fine.ratio-1)*sd.step.fine, 
                 center.sd+(coarse.to.fine.ratio-1)*sd.step.fine*2,
                 sd.step.fine)
  
  # generate parameter matrix (making use of param.combs function defined in generate_parameter_combinatinos.R)
  fine.grid = param.combs(fine.drifts,fine.sds,fine.thetas,fine.non.decision.times)
                          
  return(fine.grid)            
}