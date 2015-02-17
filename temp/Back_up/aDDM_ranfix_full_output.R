###### I AM GOING TO USE THE ADDM BASIS AND ADJUST IT TOWARDS NO ATTENTIONAL EFFECT

#   inputs needed: 

# - item distribution matrix
# - each item distribution is assigned a fixation path matrix
# - fixation duration
# - drift rate
# - theta
# - variance

aDDM.ranfix.full.output = function(x){
  
  # initialize set_size and subject
  cur.subject =  x[1]
  cur.set_size = x[2]
  drift.rate = x[3]
  theta = x[4]
  variance = x[5]
  non.decision.time = x[6]
  
  sim.reps = seq.int(1:250)
  max.sim.rep = as.integer(max(sim.reps))
  
  # initialize subset of choices.liking that we need here
  choices = core.clean.choices.liking 
  cur.choices = choices[choices$Subject == cur.subject &  
                          choices$Set_size == cur.set_size & 
                          choices$Snack_picked != 999, ]
  
  # initialize subset of core.clean.eyetracking.general that provides the fixations feeded into the aevacc_function
  eye = core.clean.eyetracking.general
  cur.eye = eye[eye$Subject == cur.subject & 
                  eye$Set_size == cur.set_size,]
  
  # Generate a vector that stores all unique refixation times
  cur.refix.dist = unique(cur.eye$Duration[cur.eye$Refix.binary == 1])
  min.refix.time = min(cur.refix.dist)
  # Generate a column that informs us how many fixations to store for a particular trial to be supplied after we ran out of empirical ones
  
  # First we need to define max. RT
  
  RT.fivenum = fivenum(cur.choices$RT_no_sacc)
  cur.max.RT = RT.fivenum[5]
  
  # max.dur gives maximum max.duration that is considered
  max.dur = as.integer(cur.max.RT);
  
  # generating column that provides a cutted (binned) version of reaction times
  cur.min.RT = RT.fivenum[1]
  
  increment.distances = (cur.max.RT - cur.min.RT) / 200
  cur.breaks = c(0,seq(cur.min.RT,cur.max.RT,increment.distances),cur.max.RT+100000)
  cur.choices['RT.bins'] = cut(cur.choices$RT_no_sacc,cur.breaks,include.lowest = TRUE)
  
  
  cur.stretched = data.table(core.clean.choices.liking.stretched[core.clean.choices.liking.stretched$Chosen.binary != -1,])
  setkey(cur.stretched,Subject,Set_size,Trial)
  temp.stretched = cur.stretched[data.table(cur.subject,cur.set_size)] 
  setkey(temp.stretched,Trial,Eligible.post.empirical)
  # In cur.eligible we store the items in the current trial that are eligible for post empirical fixations
  
  # Transition time and non.decision.time model 
  # THERE MIGHT BE A LOT OF STUFF THAT CAN BE TRIED HERE
  
  cur.choices["RT_diff"] = cur.choices$RT_Matlab - cur.choices$RT_no_sacc
  lm.fit = lm(RT_diff ~ Nr.fixations,cur.choices)
  lm.fit.coeff = as.vector(summary(lm.fit)$coefficients[,1])[]
  
  # initialize the data.frame that will be used to store the aDDM outcomes 
  nr_rows = length(unique(cur.choices$Trial))*max.sim.rep 
  trials = unique(cur.choices[,4]) 
  len.trials = length(trials)
  
  aDDM.general.output = matrix(rep(-1,nr_rows*29),nrow=nr_rows,ncol=29)
  aDDM.general.output[,1] = cur.subject
  aDDM.general.output[,2] = cur.set_size
  aDDM.general.output[,3] = rep(trials,each=max.sim.rep)
  aDDM.general.output[,4] = drift.rate
  aDDM.general.output[,5] = theta
  aDDM.general.output[,6] = variance
  
  # output cols are just the columns in which the output (Choices and RT's) will be stored
  output.cols = seq(7,28)
  
  if (cur.set_size == 4){
    # Valuations is an intermediate matrix that is then filled into the ev.update (for evidence update) vector that is passed to the DDM function
    valuations=matrix(rep(0,len.trials*4),nrow=4,ncol=len.trials)
    valuations[1,] = as.vector(cur.choices[,5])
    valuations[2,] = as.vector(cur.choices[,6])
    valuations[3,] = as.vector(cur.choices[,7])
    valuations[4,] = as.vector(cur.choices[,8])
    
    # cur.eligible will be the storage for the screen positions that are eligible for post.empirical fixations
    cur.eligible = rep(0,4)
  }
  
  if (cur.set_size == 6){
    # Valuations is an intermediate matrix that is then filled into the ev.update (for evidence update) vector that is passed to the DDM function
    valuations=matrix(rep(0,len.trials*6),nrow=6,ncol=len.trials)
    valuations[1,] = as.vector(cur.choices[,5])
    valuations[2,] = as.vector(cur.choices[,6])
    valuations[3,] = as.vector(cur.choices[,7])
    valuations[4,] = as.vector(cur.choices[,8])
    valuations[5,] = as.vector(cur.choices[,9])
    valuations[6,] = as.vector(cur.choices[,10])
    # cur.eligible will be the storage for the screen positions that are eligible for post.empirical fixations
    cur.eligible = rep(0,6)  
  }
  
  if (cur.set_size == 8){
    # Valuations is an intermediate matrix that is then filled into the ev.update (for evidence update) vector that is passed to the DDM function
    valuations=matrix(rep(0,len.trials*8),nrow=8,ncol=len.trials)
    valuations[1,] = as.vector(cur.choices[,5])
    valuations[2,] = as.vector(cur.choices[,6])
    valuations[3,] = as.vector(cur.choices[,7])
    valuations[4,] = as.vector(cur.choices[,8])
    valuations[5,] = as.vector(cur.choices[,9])
    valuations[6,] = as.vector(cur.choices[,10])
    valuations[7,] = as.vector(cur.choices[,11])
    valuations[8,] = as.vector(cur.choices[,12])
    # cur.eligible will be the storage for the screen positions that are eligible for post.empirical fixations
    cur.eligible = rep(0,8)  
  }
  
  ev.update = t(as.vector(valuations[,]*drift.rate))
  
  
  # initialization of standard deviation in "cur.sd" which is passed to aevacc_full_4
  cur.sd = sqrt(variance)
  
  nr.fix.dur.sim = 200*max.sim.rep
  nr.fix.loc.sim = 200*max.sim.rep
  
  output.row.min = 1
  output.row.max = max.sim.rep
  
  ev.update.pos.min = 1
  ev.update.pos.max = cur.set_size
  
  mean.fix.duration = mean(cur.eye$Duration)
  sd.fix.duration = sd(cur.eye$Duration)
  
  if (cur.set_size == 4){
    for (trial in trials){
      
      #################
      ##### Simulate Locations
      
      fix.loc.sim = sample(c(1,2,3,4),nr.fix.loc.sim,replace=TRUE)
      
      #################
      ##### Simulate Durations
      
      fix.dur.sim = rnorm(nr.fix.dur.sim,mean = mean.fix.duration,sd=sd.fix.duration)

      #################
      
      # Main function that spits out a vector with all choices and reaction times
      output = aevacc_ranfix_full_output_4(max.sim.rep,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],fix.loc.sim,fix.dur.sim)
      
      # aDDM.general.output is updated with the output from evacc_full_*()
      aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=max.sim.rep,ncol=22,byrow=TRUE) 
      
      output.row.min = sum(output.row.min,max.sim.rep)
      output.row.max = sum(output.row.max,max.sim.rep)
      
      ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
      ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
    }
  }
  
  if (cur.set_size == 6){
    for (trial in trials){
      
      #################
      ##### Simulate Locations
      
      fix.loc.sim = sample(c(1,2,3,4,5,6),nr.fix.loc.sim,replace=TRUE)
      
      #################
      ##### Simulate Durations
      
      fix.dur.sim = rnorm(nr.fix.dur.sim,mean = mean.fix.duration,sd=sd.fix.duration)
      
      
      #################
      
      # Main function that spits out a vector with all choices and reaction times
      output = aevacc_ranfix_full_output_6(max.sim.rep,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],fix.loc.sim,fix.dur.sim)
      
      # aDDM.general.output is updated with the output from evacc_full_*()
      aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=max.sim.rep,ncol=22,byrow=TRUE) 
      
      output.row.min = sum(output.row.min,max.sim.rep)
      output.row.max = sum(output.row.max,max.sim.rep)
      
      ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
      ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
    }
  }
  
  if (cur.set_size == 8){
    for (trial in trials){
      
      #################
      ##### Simulate Locations
      
      fix.loc.sim = sample(c(1,2,3,4,5,6,7,8),nr.fix.loc.sim,replace=TRUE)
      
      #################
      ##### Simulate Durations
      
      fix.dur.sim = rnorm(nr.fix.dur.sim,mean = mean.fix.duration,sd=sd.fix.duration)
      
      #################
      
      # Main function that spits out a vector with all choices and reaction times
      output = aevacc_ranfix_full_output_8(max.sim.rep,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],fix.loc.sim,fix.dur.sim)
      
      # aDDM.general.output is updated with the output from evacc_full_*()
      aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=max.sim.rep,ncol=22,byrow=TRUE) 
      
      output.row.min = sum(output.row.min,max.sim.rep)
      output.row.max = sum(output.row.max,max.sim.rep)
      
      ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
      ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
    }
  }
  
  ######################## FROM HERE WE ONLY GENERATE LOGLIKELIHOOD ETC.
  aDDM.general.output.frame = data.frame(Subject = aDDM.general.output[,1],Set_size = aDDM.general.output[,2],              #2
                                         Trial = aDDM.general.output[,3],Drift.Rate = aDDM.general.output[,4],
                                         Theta = aDDM.general.output[,5], Variance = aDDM.general.output[,6],
                                         Decision = aDDM.general.output[,7], Nr.fixations = aDDM.general.output[,8],
                                         RT_no_sacc = aDDM.general.output[,9], Item.last.attended = aDDM.general.output[,10],
                                         Value.last.attended = aDDM.general.output[,11] / drift.rate, Chosen.last.attended = aDDM.general.output[,12],
                                         Duration.1 = aDDM.general.output[,13], Duration.2 = aDDM.general.output[,14],
                                         Duration.3 = aDDM.general.output[,15], Duration.4 = aDDM.general.output[,16],
                                         Duration.5 = aDDM.general.output[,17], Duration.6 = aDDM.general.output[,18],
                                         Duration.7 = aDDM.general.output[,19], Duration.8 = aDDM.general.output[,20],
                                         Nr.fix.1 = aDDM.general.output[,21], Nr.fix.2 = aDDM.general.output[,22],
                                         Nr.fix.3 = aDDM.general.output[,23], Nr.fix.4 = aDDM.general.output[,24],
                                         Nr.fix.5 = aDDM.general.output[,25], Nr.fix.6 = aDDM.general.output[,26],
                                         Nr.fix.7 = aDDM.general.output[,27], Nr.fix.8 = aDDM.general.output[,28],
                                         RT_Full=aDDM.general.output[,9] + lm.fit.coeff[1] + lm.fit.coeff[2]*aDDM.general.output[,8]) #18 # IMPORTANT DONT FORGET THAT SOME MODEL OF EXTRA TIME WAS ADDED HERE
  
  #aDDM.general.output.frame['Real.choice'] = NA
  #aDDM.general.output.frame['Real.RT'] = NA
  
  #picked.vec = rep(as.vector(cur.choices$Snack_picked),each=max.sim.rep)
  #RT.bins = rep(as.vector(cur.choices$RT.bins),each=max.sim.rep)
  
  #aDDM.general.output.frame$Real.choice = picked.vec
  #aDDM.general.output.frame$Real.RT = RT.bins
  
  #   aDDM.general.output.frame['RT_no_sacc.bins'] = cut(aDDM.general.output.frame$RT_no_sacc,cur.breaks,include.lowest = TRUE)
  #   
  #   aDDM.choice.table = ddply(aDDM.general.output.frame,.(Subject,Theta,Variance,Driftrate,Trial,Choice,RT_no_sacc.bins,Real.RT,Real.choice),
  #                             summarise,
  #                             Choice.percentage = (length(RT_no_sacc.bins))/1000)
  #   
  #   
  #   #####################################################################################################################################################
  #   ############################### CALCULATING LOGLIKELIHOOD ###########################################################################################
  #   #####################################################################################################################################################
  #   
  #   row.cnt = 1
  #   aDDM.choice.likelihood = NULL
  #   aDDM.choice.likelihood = data.frame(Trial = trials, Choice.percentage = NA)
  #   
  #   for (trial in trials){
  #     cur.vec = as.vector(aDDM.choice.table$Choice.percentage[aDDM.choice.table$Trial == trial & aDDM.choice.table$RT_no_sacc.bins == aDDM.choice.table$Real.RT &
  #                                                               aDDM.choice.table$Choice == aDDM.choice.table$Real.choice])
  #     cur.length = length(cur.vec)
  #     
  #     if (cur.length >= 1){
  #       aDDM.choice.likelihood$Choice.percentage[row.cnt] = cur.vec
  #     }
  #     
  #     if (cur.length < 1){ 
  #       aDDM.choice.likelihood$Choice.percentage[row.cnt] = 0.000001
  #     }
  #     
  #     row.cnt = sum(row.cnt,1)
  #   }
  #   
  #   aDDM.choice.likelihood['Loglikelihood'] = log(aDDM.choice.likelihood$Choice.percentage)
  #   
  #   # Now store loglikelihood 
  #   total.log.lik = c(cur.subject,cur.set_size,drift.rate,variance,theta,non.decision.time,sum(aDDM.choice.likelihood$Loglikelihood))
  print(cur.subject)
  return(aDDM.general.output.frame)
  #return(total.log.lik)
  #return(aDDM.choice.likelihood)
  #return(aDDM.choice.table)
}