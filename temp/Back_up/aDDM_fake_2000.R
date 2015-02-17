###### I AM GOING TO USE THE ADDM BASIS AND ADJUST IT TOWARDS NO ATTENTIONAL EFFECT

#   inputs needed: 

# - item distribution matrix
# - each item distribution is assigned a fixation path matrix
# - fixation duration
# - drift rate
# - theta
# - variance

aDDM_fake_2000 = function(x){
  
  # initialize set_size and subject
  cur.subject =  x[1]
  cur.set_size = x[2]
  drift.rate = x[3]
  theta = x[4]
  variance = x[5]
  non.decision.time = x[6]
  
  sim.reps = seq.int(1:2000)
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
  refix.vec = unique(cur.eye$Duration[cur.eye$Refix.binary == 1])
  min.refix.time = min(refix.vec)
  # Generate a column that informs us how many fixations to store for a particular trial to be supplied after we ran out of empirical ones
  
  # First we need to define max. RT
  
  RT.fivenum = fivenum(cur.choices$RT_no_sacc[cur.choices$Set_size == cur.set_size & cur.choices$Snack_picked != 999])
  cur.max.RT = RT.fivenum[5]
  
  # max.dur gives maximum max.duration that is considered
  max.dur = as.integer(cur.max.RT);
  
  # generating column that provides a cutted (binned) version of reaction times
  cur.min.RT = RT.fivenum[1]
  
  increment.distances = (cur.max.RT - cur.min.RT) / 200
  cur.breaks = c(0,seq(cur.min.RT,cur.max.RT,increment.distances),cur.max.RT+100000)
  cur.choices['RT.bins'] = cut(cur.choices$RT_no_sacc,cur.breaks,include.lowest = TRUE)
  
  #################################################
  ########## TO BE DONE HERE ######################
  #################################################
  
  # NOW NEED TO GET SUBSET OF ITEMS THAT ARE INCLUDED IN AFTER EMPIRICAL FIXATIONS PER TRIAL
  # PRECOMPUTE A FIXATION PATHWAYS // VARIANTS STARTING OF ON EITHER ITEM....
  
  cur.stretched = data.table(core.clean.choices.liking.stretched[core.clean.choices.liking.stretched$Chosen.binary != -1,])
  setkey(cur.stretched,Subject,Set_size,Trial)
  temp.stretched = cur.stretched[data.table(cur.subject,cur.set_size)] 
  setkey(temp.stretched,Trial,Eligible.post.empirical)
  # In cur.eligible we store the items in the current trial that are eligible for post empirical fixations
  
  
  # We need to extract empirical Refixation-durations for subject and set size
  
  cur.eye = core.clean.eyetracking.general
  
  cur.refix.dist = as.vector(cur.eye$Duration[cur.eye$Subject == cur.subject & cur.eye$Set_size == cur.set_size & cur.eye$Refix.binary == 1])
  
  # Transition time and non.decision.time model 
  # THERE MIGHT BE A LOT OF STUFF THAT CAN BE TRIED HERE
  
  cur.choices["RT_diff"] = cur.choices$RT_Matlab - cur.choices$RT_no_sacc
  lm.fit = lm(RT_diff ~ Nr.fixations,cur.choices)
  lm.fit.coeff = as.vector(summary(lm.fit)$coefficients[,1])[]
  
  # initialize the data.frame that will be used to store the aDDm outcomes 
  nr_rows = length(unique(cur.choices$Trial))*max.sim.rep 
  trials = unique(cur.choices[,4]) 
  len.trials = length(trials)
  
  aDDM.general.output = matrix(rep(-1,nr_rows*20),nrow=nr_rows,ncol=20)
  aDDM.general.output[,1] = cur.subject
  aDDM.general.output[,2] = cur.set_size
  aDDM.general.output[,3] = rep(trials,each=max.sim.rep)
  aDDM.general.output[,4] = rep(sim.reps,times=len.trials)
  aDDM.general.output[,5] = rep(cur.choices[,5],each=max.sim.rep)
  aDDM.general.output[,6] = rep(cur.choices[,6],each=max.sim.rep)
  aDDM.general.output[,7] = rep(cur.choices[,7],each=max.sim.rep)
  aDDM.general.output[,8] = rep(cur.choices[,8],each=max.sim.rep)
  aDDM.general.output[,9] = rep(cur.choices[,9],each=max.sim.rep)
  aDDM.general.output[,10] = rep(cur.choices[,10],each=max.sim.rep)
  aDDM.general.output[,11] = rep(cur.choices[,11],each=max.sim.rep)
  aDDM.general.output[,12] = rep(cur.choices[,12],each=max.sim.rep)
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
  
  ev.update = t(as.vector(valuations[,]*drift.rate))
  
  
  # initialization of standard deviation in "cur.sd" which is passed to aevacc_full_4
  cur.sd = sqrt(variance)
  
  nr.fix.dur.sim = 100*max.sim.rep
  nr.fix.loc.sim = 200*max.sim.rep
  
  output.row.min = 1
  output.row.max = max.sim.rep
  
  ev.update.pos.min = 1
  ev.update.pos.max = cur.set_size
  
  if (cur.set_size == 4){
    for (trial in trials){
      
      #################
      ##### Extract Empirical Fixations
      
      cur.fixations = as.vector(cur.eye$Item_attended[cur.eye$Trial == trial])
      
      #################
      ##### Extract Empirical Durations
      
      cur.durations = as.vector(cur.eye$Duration[cur.eye$Trial == trial])
      
      
      temp = temp.stretched[data.table(trial)]
      cur.eligible = temp$Item[temp$Eligible.post.empirical == 1]
      
      #################
      #### Here I draw fixation location samples : IN C++ WE LIMIT TO 100 FIXATIONS --> SAMPLE 200 per repetition here for absolute worst case 
      
      fix.loc.sim = sample(cur.eligible,nr.fix.loc.sim,replace=TRUE)
      
      #################
      #### I draw fixation durations samples here --> THIS FOLLOWS THE RESTRICTION THAT THERE CAN BE ONLY 100 FIXATIONS IN ADDITION
      
      fix.dur.sim = sample(cur.refix.dist,nr.fix.dur.sim,replace=TRUE)
      
      #################
      
      # Main function that spits out a vector with all choices and reaction times
      output = aevacc_full_4(max.sim.rep,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],cur.fixations,cur.durations,fix.loc.sim,fix.dur.sim)
      
      # aDDM.general.output is updated with the output from evacc_full_*()
      aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=max.sim.rep,ncol=3,byrow=TRUE) 
      
      output.row.min = sum(output.row.min,max.sim.rep)
      output.row.max = sum(output.row.max,max.sim.rep)
      
      ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
      ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
    }
  }
  
  if (cur.set_size == 6){
    for (trial in trials){
      
      #################
      ##### Extract Empirical Fixations
      
      cur.fixations = as.vector(cur.eye$Item_attended[cur.eye$Trial == trial])
      
      #################
      ##### Extract Empirical Durations
      
      cur.durations = as.vector(cur.eye$Duration[cur.eye$Trial == trial])
      
      
      temp = temp.stretched[data.table(trial)]
      cur.eligible = temp$Item[temp$Eligible.post.empirical == 1]
      
      #################
      #### Here I draw fixation location samples : IN C++ WE LIMIT TO 100 FIXATIONS --> SAMPLE 200 per repetition here for absolute worst case 
      
      fix.loc.sim = sample(cur.eligible,nr.fix.loc.sim,replace=TRUE)
      
      #################
      #### I draw fixation durations samples here --> THIS FOLLOWS THE RESTRICTION THAT THERE CAN BE ONLY 100 FIXATIONS IN ADDITION
      
      fix.dur.sim = sample(cur.refix.dist,nr.fix.dur.sim,replace=TRUE)
      
      #################
      
      # Main function that spits out a vector with all choices and reaction times
      output = aevacc_full_6(max.sim.rep,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],cur.fixations,cur.durations,fix.loc.sim,fix.dur.sim)
      
      # aDDM.general.output is updated with the output from evacc_full_*()
      aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=max.sim.rep,ncol=3,byrow=TRUE) 
      
      output.row.min = sum(output.row.min,max.sim.rep)
      output.row.max = sum(output.row.max,max.sim.rep)
      
      ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
      ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
    }
  }
  
  if (cur.set_size == 8){
    for (trial in trials){
      
      #################
      ##### Extract Empirical Fixations
      
      cur.fixations = as.vector(cur.eye$Item_attended[cur.eye$Trial == trial])
      
      #################
      ##### Extract Empirical Durations
      
      cur.durations = as.vector(cur.eye$Duration[cur.eye$Trial == trial])
      
      
      temp = temp.stretched[data.table(trial)]
      cur.eligible = temp$Item[temp$Eligible.post.empirical == 1]
      
      #################
      #### Here I draw fixation location samples : IN C++ WE LIMIT TO 100 FIXATIONS --> SAMPLE 200 per repetition here for absolute worst case 
      
      fix.loc.sim = sample(cur.eligible,nr.fix.loc.sim,replace=TRUE)
      
      #################
      #### I draw fixation durations samples here --> THIS FOLLOWS THE RESTRICTION THAT THERE CAN BE ONLY 100 FIXATIONS IN ADDITION
      
      fix.dur.sim = sample(cur.refix.dist,nr.fix.dur.sim,replace=TRUE)
      
      #################
      
      # Main function that spits out a vector with all choices and reaction times
      output = aevacc_full_8(max.sim.rep,max.dur,cur.sd,theta,ev.update[ev.update.pos.min:ev.update.pos.max],cur.fixations,cur.durations,fix.loc.sim,fix.dur.sim)
      
      # aDDM.general.output is updated with the output from evacc_full_*()
      aDDM.general.output[output.row.min:output.row.max,output.cols] = matrix(output,nrow=max.sim.rep,ncol=3,byrow=TRUE) 
      
      output.row.min = sum(output.row.min,max.sim.rep)
      output.row.max = sum(output.row.max,max.sim.rep)
      
      ev.update.pos.min = sum(ev.update.pos.min,cur.set_size)
      ev.update.pos.max = sum(ev.update.pos.max,cur.set_size)
    }
  }
  print("finished simulating")
  ######################## FROM HERE WE ONLY GENERATE LOGLIKELIHOOD ETC.
  aDDM.general.output.frame = data.frame(Subject=aDDM.general.output[,1],Set_size=aDDM.general.output[,2],              #2
                                         Trial=aDDM.general.output[,3],Repetition=aDDM.general.output[,4],              #4
                                         Snack_1=aDDM.general.output[,5],Snack_2=aDDM.general.output[,6],               #6
                                         Snack_3=aDDM.general.output[,7],Snack_4=aDDM.general.output[,8],               #8
                                         Snack_5=aDDM.general.output[,9],Snack_6=aDDM.general.output[,10],             #10
                                         Snack_7=aDDM.general.output[,11],Snack_8=aDDM.general.output[,12],            #12
                                         Choice=aDDM.general.output[,13],Fixation.nr=aDDM.general.output[,14],         #14
                                         Fixation.pos=aDDM.general.output[,15], RT_no_sacc=aDDM.general.output[,16],   #16
                                         RT_Full=aDDM.general.output[,16] + lm.fit.coeff[1] + lm.fit.coeff[2]*aDDM.general.output[,14],Driftrate=aDDM.general.output[,18],    #18 # IMPORTANT DONT FORGET THAT SOME MODEL OF EXTRA TIME WAS ADDED HERE
                                         Theta=aDDM.general.output[,19],Variance=aDDM.general.output[,20])            #20
  
  
  aDDM.general.output.frame['RT_no_sacc.bins'] = cut(aDDM.general.output.frame$RT_no_sacc,cur.breaks,include.lowest = TRUE)
  
  
  if (cur.set_size == 4){
    aDDM.general.output.frame.fake = core.aDDM.general.output.frame.fake.4
  }
  
  if (cur.set_size == 6){
    aDDM.general.output.frame.fake = core.aDDM.general.output.frame.fake.6
  }
  
  if (cur.set_size == 8){
    aDDM.general.output.frame.fake = core.aDDM.general.output.frame.fake.8
  }
  
  
  aDDM.choice.table = ddply(aDDM.general.output.frame,.(Subject,Theta,Variance,Driftrate,Trial,Choice,RT_no_sacc.bins),
                            summarise,
                            Choice.percentage = (length(RT_no_sacc.bins))/2000)
  
  aDDM.choice.table.fake = ddply(aDDM.general.output.frame.fake,.(Subject,Theta,Variance,Driftrate,Trial,Choice,RT_no_sacc.bins),
                                 summarise,
                                 Choice.percentage = (length(RT_no_sacc.bins))/5000)
  
  #####################################################################################################################################################
  ############################### CALCULATING LOGLIKELIHOOD ###########################################################################################
  #####################################################################################################################################################
  
  vec.pos.cnt = 1
  aDDM.choice.likelihood = rep(0,length(aDDM.choice.table.fake[,1]))
  
  for (trial in trials){
    
    choices = unique(aDDM.choice.table.fake$Choice[aDDM.choice.table.fake$Trial == trial])
    
    for (choice in choices){
      
      cur.dat.fake = aDDM.choice.table.fake[aDDM.choice.table.fake$Trial == trial & aDDM.choice.table.fake$Choice == choice,]
      cur.dat.sim = aDDM.choice.table[aDDM.choice.table$Trial == trial & aDDM.choice.table$Choice == choice,]
      
      if (length(cur.dat.sim[,1]) >= 1){
        
        cur.RT.bins = unique(cur.dat.fake$RT_no_sacc.bins[cur.dat.fake$Trial == trial & cur.dat.fake$Choice == choice])  
        
        for (rt.bin in cur.RT.bins){
          
          fake.p = as.vector(cur.dat.fake$Choice.percentage[cur.dat.fake$RT_no_sacc.bins == rt.bin])
          sim.p = as.vector(cur.dat.sim$Choice.percentage[cur.dat.sim$RT_no_sacc.bins == rt.bin])
          cur.length = length(sim.p)
          
          if (cur.length >= 1){
            aDDM.choice.likelihood[vec.pos.cnt] = fake.p[1]*log(sim.p[1])
          }
          
          if (cur.length < 1){ 
            aDDM.choice.likelihood[vec.pos.cnt] = fake.p[1]*log(0.0000001)
          }
          
          vec.pos.cnt = vec.pos.cnt + 1
        }
      }
    }
  }
  
  # Now store loglikelihood 
  total.log.lik = c(cur.subject,cur.set_size,drift.rate,variance,theta,non.decision.time,sum(aDDM.choice.likelihood))
  
  #return(aDDM.general.output.frame)
  print(total.log.lik)
  return(total.log.lik)
}