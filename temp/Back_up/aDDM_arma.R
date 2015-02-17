###### I AM GOING TO USE THE ADDM BASIS AND ADJUST IT TOWARDS NO ATTENTIONAL EFFECT

#   inputs needed: 

# - item distribution matrix
# - each item distribution is assigned a fixation path matrix
# - fixation duration
# - drift rate
# - theta
# - variance

# - number of repetitions used -- nr_reps
# - fixaion.model -- 'Normal' or 'Random'

aDDM_arma = function(x,nr_reps,fixation.model){
  
  # INITIALIZATION OF PARAMETERS AND RELEVANT SUBSETS OF DATA FRAMES
  #-------------------------------------------------------------------------------------------------------------------------------------
  # initialize set_size and subject
  cur.subject =  x[1]
  cur.set_size = x[2]
  
  # initialize model parameters
  drift.rate = x[3]
  theta = x[4]
  cur.sd = x[5]
  non.decision.time = x[6]
  
  # initialize subset of choices.liking that we need here
  cur.choices = core.clean.choices.liking 
  cur.choices = cur.choices[cur.choices$Subject == cur.subject &  
                              cur.choices$Set_size == cur.set_size & 
                              cur.choices$Snack_picked != 999, ]
  
  # initialize subset of core.clean.eyetracking.general that provides the fixations feeded into the aevacc_function
  eye = core.clean.eyetracking.general
  cur.eye = eye[eye$Subject == cur.subject & 
                  eye$Set_size == cur.set_size,]
  
  # initialize subset of core.clean.choices.liking.stretched the long version of the core.clean.choices.liking data.frame
  cur.stretched = data.table(core.clean.choices.liking.stretched[core.clean.choices.liking.stretched$Chosen.binary != -1 & 
                                                                   core.clean.choices.liking.stretched$Subject == cur.subject &
                                                                   core.clean.choices.liking.stretched$Set_size == cur.set_size,]) 
  setkey(cur.stretched,Trial,Eligible.post.empirical)
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  # SOME MISCELLANEOUS VARIABLES THAT ARE UTILIZED LATER
  #-------------------------------------------------------------------------------------------------------------------------------------
  # First we need to define max. RT
  RT.fivenum = fivenum(cur.choices$RT_no_sacc)
  cur.min.RT = RT.fivenum[1]
  cur.max.RT = RT.fivenum[5]
  max.dur = as.integer(cur.max.RT);
  
  # generating column that provides a cutted (binned) version of reaction times
  increment.distances = (cur.max.RT - cur.min.RT) / 200
  cur.breaks = c(0,seq(cur.min.RT,cur.max.RT,increment.distances),cur.max.RT+100000)
  cur.choices['RT.bins'] = cut(cur.choices$RT_no_sacc,cur.breaks,include.lowest = TRUE)
  
  # We need to extract empirical Refixation-durations for subject and set size
  cur.refix.dist = as.vector(cur.eye$Duration[cur.eye$Refix.binary == 1])
  
  # Amount of fixation fixation durations and locations we will simulate in total
  nr.fix.dur.sim = 100*nr_reps
  nr.fix.loc.sim = 100*nr_reps
  
  # We will use the mean fixation duration and standard fixation duration for simulated fixations in the aDDM
  mean.fix.duration = mean(core.clean.eyetracking.general$Duration)
  sd.fix.duration = sd(core.clean.eyetracking.general$Duration)
  
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  # TRANSITION TIME AND NON-DECISION TIME MODEL
  #-------------------------------------------------------------------------------------------------------------------------------------
  # THERE MIGHT BE A LOT OF STUFF THAT CAN BE TRIED HERE
  
  cur.choices["RT_diff"] = cur.choices$RT_Matlab - cur.choices$RT_no_sacc
  lm.fit = lm(RT_diff ~ Nr.fixations,cur.choices)
  lm.fit.coeff = as.vector(summary(lm.fit)$coefficients[,1])[]
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  
  # INITIALIZATION OF ALL DATA FRAMES THAT WE NEED TO STORE ADDM RESULTS IN
  #-------------------------------------------------------------------------------------------------------------------------------------
  nr_rows = length(unique(cur.choices$Trial))*nr_reps 
  trials = unique(cur.choices[,4]) 
  len.trials = length(trials)
  
  aDDM.general.output = matrix(rep(-1,nr_rows*20),nrow=nr_rows,ncol=20)
  aDDM.general.output[,1] = cur.subject
  aDDM.general.output[,2] = cur.set_size
  aDDM.general.output[,3] = rep(trials,each=nr_reps)
  aDDM.general.output[,4] = rep(seq.int(1:nr_reps),times=len.trials)
  aDDM.general.output[,5] = rep(cur.choices[,5],each=nr_reps)
  aDDM.general.output[,6] = rep(cur.choices[,6],each=nr_reps)
  aDDM.general.output[,7] = rep(cur.choices[,7],each=nr_reps)
  aDDM.general.output[,8] = rep(cur.choices[,8],each=nr_reps)
  aDDM.general.output[,9] = rep(cur.choices[,9],each=nr_reps)
  aDDM.general.output[,10] = rep(cur.choices[,10],each=nr_reps)
  aDDM.general.output[,11] = rep(cur.choices[,11],each=nr_reps)
  aDDM.general.output[,12] = rep(cur.choices[,12],each=nr_reps)
  aDDM.general.output[,13] = -1
  aDDM.general.output[,15:16] = -1
  aDDM.general.output[,18] = drift.rate
  aDDM.general.output[,19] = theta
  aDDM.general.output[,20] = variance
  # output cols are just the columns in which the output (Choices and RT's) will be stored
  output.cols = c(13,14,16)
  
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
  
  # ev.update is the fixed shift in mean RDV per time-step by item or each trial, we compute the here once, to avoid unnecessary recomputation in the simulation step
  ev.update = t(as.vector(valuations[,]*drift.rate))
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  # WE LOOP THROUGH ALL TRIALS FOR THE ADDM RUNS
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  output.row.min = 1
  output.row.max = nr_reps
  
  ev.update.pos.min = 1
  ev.update.pos.max = cur.set_size
  
  noisevec = rnorm(100000,0,cur.sd)
  
  if (cur.set_size == 4){
    if (fixation.model == "Normal"){
      for (trial in trials){
        
        #################
        ##### Extract Empirical Fixations
        
        cur.fixations = as.vector(cur.eye$Item_attended[cur.eye$Trial == trial])
        
        #################
        ##### Extract Empirical Durations
        
        cur.durations = as.vector(cur.eye$Duration[cur.eye$Trial == trial])
        
        
        temp = cur.stretched[data.table(trial)]
        cur.eligible = temp$Item[temp$Eligible.post.empirical == 1]
        
        #################
        #### Here I draw fixation location samples : IN C++ WE LIMIT TO 100 FIXATIONS --> SAMPLE 200 per repetition here for absolute worst case 
        
        fix.loc.sim = sample(cur.eligible,nr.fix.loc.sim,replace=TRUE)
        
        #################
        #### I draw fixation durations samples here --> THIS FOLLOWS THE RESTRICTION THAT THERE CAN BE ONLY 100 FIXATIONS IN ADDITION
        
        fix.dur.sim = sample(cur.refix.dist,nr.fix.dur.sim,replace=TRUE)
        
        #################
        
        # Main function that spits out a vector with all choices and reaction times
        output = aevacc_arma_4(nr_reps,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],cur.fixations,cur.durations,fix.loc.sim,fix.dur.sim,noisevec)
        
        # aDDM.general.output is updated with the output from evacc_full_*()
        aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr_reps,ncol=3,byrow=TRUE) 
        
        output.row.min = sum(output.row.min,nr_reps)
        output.row.max = sum(output.row.max,nr_reps)
        
        ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
        ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
      }
    }
    
    if (fixation.model == "Random"){
      for (trial in trials){
        
        #################
        ##### Simulate Locations
        
        fix.loc.sim = sample(c(1,2,3,4),nr.fix.loc.sim,replace=TRUE)
        
        #################
        ##### Simulate Durations
        
        fix.dur.sim = rnorm(nr.fix.dur.sim,mean = mean.fix.duration,sd=sd.fix.duration)
        
        #################
        
        #################
        
        # Main function that spits out a vector with all choices and reaction times
        output = aevacc_ranfix_4(nr_reps,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],fix.loc.sim,fix.dur.sim)
        
        # aDDM.general.output is updated with the output from evacc_full_*()
        aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr_reps,ncol=3,byrow=TRUE) 
        
        output.row.min = sum(output.row.min,nr_reps)
        output.row.max = sum(output.row.max,nr_reps)
        
        ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
        ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
      }
    }
  }
  
  if (cur.set_size == 6){
    if(fixation.model == "Normal"){
      for (trial in trials){
        
        #################
        ##### Extract Empirical Fixations
        
        cur.fixations = as.vector(cur.eye$Item_attended[cur.eye$Trial == trial])
        
        #################
        ##### Extract Empirical Durations
        
        cur.durations = as.vector(cur.eye$Duration[cur.eye$Trial == trial])
        
        
        temp = cur.stretched[data.table(trial)]
        cur.eligible = temp$Item[temp$Eligible.post.empirical == 1]
        
        #################
        #### Here I draw fixation location samples : IN C++ WE LIMIT TO 100 FIXATIONS --> SAMPLE 200 per repetition here for absolute worst case 
        
        fix.loc.sim = sample(cur.eligible,nr.fix.loc.sim,replace=TRUE)
        
        #################
        #### I draw fixation durations samples here --> THIS FOLLOWS THE RESTRICTION THAT THERE CAN BE ONLY 100 FIXATIONS IN ADDITION
        
        fix.dur.sim = sample(cur.refix.dist,nr.fix.dur.sim,replace=TRUE)
        
        #################
        
        # Main function that spits out a vector with all choices and reaction times
        output = aevacc_arma_6(nr_reps,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],cur.fixations,cur.durations,fix.loc.sim,fix.dur.sim,noisevec)
        
        # aDDM.general.output is updated with the output from evacc_full_*()
        aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr_reps,ncol=3,byrow=TRUE) 
        
        output.row.min = sum(output.row.min,nr_reps)
        output.row.max = sum(output.row.max,nr_reps)
        
        ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
        ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
      }
    }
    
    if (fixation.model == "Random"){
      for (trial in trials){
        
        #################
        ##### Simulate Locations
        
        fix.loc.sim = sample(c(1,2,3,4,5,6),nr.fix.loc.sim,replace=TRUE)
        
        #################
        ##### Simulate Durations
        
        fix.dur.sim = rnorm(nr.fix.dur.sim,mean = mean.fix.duration,sd=sd.fix.duration)
        
        #################
        
        #################
        
        # Main function that spits out a vector with all choices and reaction times
        output = aevacc_ranfix_6(nr_reps,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],fix.loc.sim,fix.dur.sim)
        
        # aDDM.general.output is updated with the output from evacc_full_*()
        aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr_reps,ncol=3,byrow=TRUE) 
        
        output.row.min = sum(output.row.min,nr_reps)
        output.row.max = sum(output.row.max,nr_reps)
        
        ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
        ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
      }
    }
    
  }
  
  if (cur.set_size == 8){
    if (fixation.model == "Normal"){
      for (trial in trials){
        
        #################
        ##### Extract Empirical Fixations
        
        cur.fixations = as.vector(cur.eye$Item_attended[cur.eye$Trial == trial])
        
        #################
        ##### Extract Empirical Durations
        
        cur.durations = as.vector(cur.eye$Duration[cur.eye$Trial == trial])
        
        
        temp = cur.stretched[data.table(trial)]
        cur.eligible = temp$Item[temp$Eligible.post.empirical == 1]
        
        #################
        #### Here I draw fixation location samples : IN C++ WE LIMIT TO 100 FIXATIONS --> SAMPLE 200 per repetition here for absolute worst case 
        
        fix.loc.sim = sample(cur.eligible,nr.fix.loc.sim,replace=TRUE)
        
        #################
        #### I draw fixation durations samples here --> THIS FOLLOWS THE RESTRICTION THAT THERE CAN BE ONLY 100 FIXATIONS IN ADDITION
        
        fix.dur.sim = sample(cur.refix.dist,nr.fix.dur.sim,replace=TRUE)
        
        #################
        
        # Main function that spits out a vector with all choices and reaction times
        output = aevacc_arma_8(nr_reps,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],cur.fixations,cur.durations,fix.loc.sim,fix.dur.sim,noisevec)
        
        # aDDM.general.output is updated with the output from evacc_full_*()
        aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr_reps,ncol=3,byrow=TRUE) 
        
        output.row.min = sum(output.row.min,nr_reps)
        output.row.max = sum(output.row.max,nr_reps)
        
        ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
        ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
      }
    }
    
    if (fixation.model == "Random"){
      for (trial in trials){
        
        #################
        ##### Simulate Locations
        
        fix.loc.sim = sample(c(1,2,3,4,5,6,7,8),nr.fix.loc.sim,replace=TRUE)
        
        #################
        ##### Simulate Durations
        
        fix.dur.sim = rnorm(nr.fix.dur.sim,mean = mean.fix.duration,sd=sd.fix.duration)
        
        #################
        
        #################
        
        # Main function that spits out a vector with all choices and reaction times
        output = aevacc_ranfix_8(nr_reps,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],fix.loc.sim,fix.dur.sim)
        
        # aDDM.general.output is updated with the output from evacc_full_*()
        aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=nr_reps,ncol=3,byrow=TRUE) 
        
        output.row.min = sum(output.row.min,nr_reps)
        output.row.max = sum(output.row.max,nr_reps)
        
        ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
        ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
      }
    }
  }
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  # FROM HERE WE ONLY GENERATE LOGLIKELIHOOD ETC.
  #-------------------------------------------------------------------------------------------------------------------------------------
  
  # Generate data.table from outputs
  
  aDDM.general.output.frame = data.table(Subject=aDDM.general.output[,1],Set_size=aDDM.general.output[,2],              #2
                                         Trial=aDDM.general.output[,3],Repetition=aDDM.general.output[,4],              #4
                                         Snack_1=aDDM.general.output[,5],Snack_2=aDDM.general.output[,6],               #6
                                         Snack_3=aDDM.general.output[,7],Snack_4=aDDM.general.output[,8],               #8
                                         Snack_5=aDDM.general.output[,9],Snack_6=aDDM.general.output[,10],             #10
                                         Snack_7=aDDM.general.output[,11],Snack_8=aDDM.general.output[,12],            #12
                                         Choice=aDDM.general.output[,13],Fixation.nr=aDDM.general.output[,14],         #14
                                         Fixation.pos=aDDM.general.output[,15], RT_no_sacc=aDDM.general.output[,16],   #16
                                         RT_Full=aDDM.general.output[,16] + lm.fit.coeff[1] + lm.fit.coeff[2]*aDDM.general.output[,14],Driftrate=aDDM.general.output[,18],    #18 # IMPORTANT DONT FORGET THAT SOME MODEL OF EXTRA TIME WAS ADDED HERE
                                         Theta=aDDM.general.output[,19],Variance=aDDM.general.output[,20])             #20
  
  # Add real choices and RTs to the output data frame to enable checking of likelihood below
  aDDM.general.output.frame$Real.choice = rep(as.vector(cur.choices$Snack_picked),each=nr_reps)
  aDDM.general.output.frame$Real.RT = rep(as.vector(cur.choices$RT.bins),each=nr_reps)  
  aDDM.general.output.frame$RT_no_sacc.bins = cut(aDDM.general.output.frame$RT_no_sacc,cur.breaks,include.lowest = TRUE)
  
  # We generate a table that gives the percentage of choices which go to each particular item, fall under each Reaction time bin, and are matched with Real.choice
  
  aDDM.choice.table = aDDM.general.output.frame %>% group_by(Trial,Choice,RT_no_sacc.bins,Real.RT,Real.choice) %>% summarise(count = n()) %>% mutate(choice.percentage = count / nr_reps)
  
  
  #####################################################################################################################################################
  ############################### CALCULATING LOGLIKELIHOOD ###########################################################################################
  #####################################################################################################################################################
  
  row.cnt = 1
  aDDM.choice.likelihood = data.table(Trial = trials, choice.percentage = 0.000001)
  
  setkey(aDDM.choice.table,Trial)
  
  for (trial in trials){
    cur.subset = aDDM.choice.table[J(trial),]
    cur.vec = as.vector(cur.subset$choice.percentage[cur.subset$RT_no_sacc.bins == cur.subset$Real.RT &
                                                       cur.subset$Choice == cur.subset$Real.choice])
    cur.length = length(cur.vec)
    
    if (cur.length >= 1){
      aDDM.choice.likelihood$choice.percentage[row.cnt] = cur.vec
    }
    
    row.cnt = sum(row.cnt,1)
  }
  
  aDDM.choice.likelihood$Loglikelihood = log(aDDM.choice.likelihood$choice.percentage)
  
  # Now store loglikelihood 
  
  total.log.lik = c(cur.subject,cur.set_size,drift.rate,variance,theta,non.decision.time,sum(aDDM.choice.likelihood$Loglikelihood))
  print(total.log.lik)
  
  #return(aDDM.choice.likelihood)
  #return(aDDM.choice.table)
  #return(aDDM.general.output.frame)
  
  return(total.log.lik)
}