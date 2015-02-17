check.gridsearch.completeness = function(sim.log.parameters,original.parameters){
  
  original.parameters = original.parameters
  sim.log.parameters = sim.log.parameters
  
  cur.len.sim = length(sim.log.parameters[,1])
  cur.len.orig = length(original.parameters[,1])
  
  missing.params = 0 
  
  all.same = 1
  
  for (i in 1:cur.len.orig){
    cur.param.set = original.parameters[i,] 
    same = 0
    for (j in 1:cur.len.sim){
      same = identical(cur.param.set,as.vector(sim.log.parameters[j,])) 
      
      if (same == 1){
        break
      }
    }

    if (same == 0){
      if (i == 1){
        missing.params = cur.param.set
        all.same[1] = 0
      }
      if (i > 1){
        missing.params = rbind(missing.params,cur.param.set)
        all.same[1] = 0
      }
    }
    
  }
  return(list("All.complete" = all.same, "Missing.parameters" = missing.params))
}