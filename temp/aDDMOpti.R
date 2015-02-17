# Author: Alexander Fengler
# Date: Jan 16th 2015
# Purpose: Main function running aDDM optimization based on grid search algorithm

# Several options: run by subject, run with fake dataset, run over all subjects

aDDMOpti= function(subjects,                    # Subject to be tested: -1 = Model testing with fake data // 0 = Fit over all subjects // > 0 fit by subject
         set.sizes,                   # Set sizes to be tested
         choice.dat,                  # Choice data - Needs: All item values per condition as seperate columns by item, item chosen,subject, set.sizes, trialid, RT
         eye.dat,                     # Fixation data - Needs: Fixation Locations, Fixation Durations, subject, set.sizes, trialid (sorted within trial according to fixation number)
         drifts,                      # Minimum Drift Rate we consider in grid
         thetas,                      # Minimum Theta we consider in grid
         sds,                         # Minimum SD we consider in grid
         non.decision.times,          # Minimum Non Decision Time
         timestep.ms,                 # Timestep-size in milliseconds for aDDM simulation propagation
         nr.reps,                     # Number of repitions used in simulation runs
         model.type,                  # Currently choice between model with memory effects, divided in no drift for items not seen ('mem') and no drift and no sd for items not seen ('memzeronoise') and without ('nomem') 
         output.type,                 # Output type can be: "Opti" - Normal loglikelihood for model optimization, "FakeOpti" - loglikelihoods for optimization for artificially created data
         fixation.model,              # Fixation model used in simulations: "Normal" - Real Fixations used, "Random"
         allow.extension,             # Binary, whether we allow extensions of the grid in case we find corner-solution
         allow.fine.grid,             # Binary, whether we allow the fine-grid step in the grid search procedure
         generate)                    # Generate (binary) tells the aDDM function whether to spit out full data (1) or log likelihoods (0)
{
  # LOAD ALL PACKAGES AND FUNCTIONS NECESSARY 
  # ----------------------------------------------------------------------------------------------
  # Load all high level functions needed 
  source('Analysis_Visualization/aDDM_Simulation/addm_opti_supportfuns/generate_parameter_combinations.R')
  source('Analysis_Visualization/aDDM_Simulation/addm_opti_supportfuns/aDDM_gridsearch.R')
  source('Analysis_Visualization/aDDM_Simulation/addm_opti_supportfuns/check_grid_search_completeness.R')
  source('Analysis_Visualization/aDDM_Simulation/addm_opti_supportfuns/check_and_update_corners.R')
  source('Analysis_Visualization/aDDM_Simulation/addm_opti_supportfuns/generate_fine_grid.R')
  
  
  source('Analysis_Visualization/aDDM_Simulation/load_packages_scripts_aDDM.R')
  # Initialize all packages
  load.aDDM.packages()
  
  # Initialize parameters, subjects, set_size
  load.aDDM.scripts()
  
  # INITIALIZING VARIABLES THAT CAN BE MANIPULATED 
  ###############################################################################################
  # ---------------------------------------------------------------------------------------------
  
  # Initialize the logfile we are going to use for storing results from the gridsearch
  cur.log.file = "Sim_logs/cur_log.txt"
  writeLines(c(""),cur.log.file)
  
  # In case we find corner solutions what parameter shift shall we allow per round (balance between degree of exploration and time spend computing)
  drift.shift = (drifts[length(drifts)] - drifts[length(drifts)-1])*3
  sd.shift = (sds[length(sds)] - sds[length(sds)-1])*3
  non.decision.time.shift = (non.decision.times[length(non.decision.times)] - non.decision.times[length(non.decision.times)-1])*3
  
  # Parameters for fine grid search
  # Distance between two points on the grid by dimension (drift.rate,sd)
  coarse.to.fine.ratio = 5
  
  drift.step.fine = (drifts[length(drifts)] - drifts[length(drifts)-1])/coarse.to.fine.ratio
  theta.step.fine = (thetas[length(thetas)] - thetas[length(thetas)-1])/coarse.to.fine.ratio
  sd.step.fine = (sds[length(sds)] - sds[length(sds)-1])/coarse.to.fine.ratio
  non.decision.time.step.fine = (non.decision.times[length(non.decision.times)] - non.decision.times[length(non.decision.times)-1])/coarse.to.fine.ratio
  # ---------------------------------------------------------------------------------------------
  ###############################################################################################
  
  for (cur.subject in subjects){
    for (cur.set_size in set.sizes){
      # GENERATING SUBSETS OF DATA.FRAMES ACCORDING TO SPECIFICATIONS OF SUBJECTS AND SET SIZES TO INCLUDE 
      # WHICH ARE GIVEN ABOVE
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      
      # EXTRACT RELEVANT VALUATIONS
      if (cur.set_size == 4){
        if (cur.subject == -1){
          cur.choice.dat = core.clean.model.test.dat$items_4
        } else if (cur.subject == 0){
          cur.choice.dat =  choice.dat[choice.dat$set_size == cur.set_size, list(v1,v2,v3,v4,decision,rt,rtup,rtdown,trialid)]
        } else {
          cur.choice.dat =  choice.dat[choice.dat$subject == cur.subject & 
                                        choice.dat$set_size == cur.set_size , list(v1,v2,v3,v4,decision,rt,rtup,rtdown,trialid)]
        }
      }
      
      if (cur.set_size == 6){
        if (cur.subject == -1){
          cur.choice.dat = core.clean.model.test.dat$items_6
        } else if (cur.subject == 0){
          cur.choice.dat =  choice.dat[choice.dat$set_size == cur.set_size,list(v1,v2,v3,v4,v5,v6,decision,rt,rtup,rtdown,trialid)]
        } else{ 
          cur.choice.dat =  choice.dat[choice.dat$subject == cur.subject & 
                                        choice.dat$set_size == cur.set_size,list(v1,v2,v3,v4,v5,v6,decision,rt,rtup,rtdown,trialid)]
        }
      }
      
      if (cur.set_size == 8){
        if (cur.subject == -1){
          cur.choice.dat = core.clean.model.test.dat$items_8
        } else if (cur.subject == 0){
          cur.choice.dat =  choice.dat[choice.dat$set_size == cur.set_size,list(v1,v2,v3,v4,v5,v6,v7,v8,decision,rt,rtup,rtdown,trialid)]
        }else{
          cur.choice.dat =  choice.dat[choice.dat$subject == cur.subject & 
                                        choice.dat$set_size == cur.set_size,list(v1,v2,v3,v4,v5,v6,v7,v8,decision,rt,rtup,rtdown,trialid)]
        }
      }
      
      # EXTRACT RELEVANT FIXATIONS
      if (cur.subject == -1) {
        cur.eye.dat = core.clean.mode.test.dat$fixations
      } else if (cur.subject == 0){
        cur.eye.dat = eye.dat[eye.dat$set_size == cur.set_size, list(fixloc,fixnr,fixdur,trialid)]
      } else{
        cur.eye.dat = eye.dat[eye.dat$subject == cur.subject &
                                eye.dat$set_size == cur.set_size, list(fixloc,fixnr,fixdur,trialid)]
      }
  
      # Set keys for the data.tables passed into function
      setkey(cur.choice.dat,trialid)
      setkey(cur.eye.dat,trialid)
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      
      
      # GENERATE PARAMETER MATRIX (ALL COMBINATIONS OF PARAMETERS SPECIFIED ABOVE)
      ###############################################################################################
      # ---------------------------------------------------------------------------------------------
      parameter.matrix = param.combs(drifts,
                                     sds,
                                     thetas,
                                     non.decision.times)
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      
      
      # RUN GRID SEARCH
      ###############################################################################################
      # ---------------------------------------------------------------------------------------------
      log.liks = aDDM.gridsearch(cur.choice.dat,
                                 cur.eye.dat,
                                 parameter.matrix,
                                 cur.subject,
                                 cur.set_size,
                                 output.type,
                                 nr.reps,
                                 model.type,
                                 fixation.model,
                                 cur.log.file,
                                 timestep.ms,
                                 generate)
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      
      
      # CHECK WHETHER GRID SEARCH WAS COMPLETE AND RUN MISSING PARAMETERS IF NOT
      ###############################################################################################
      # ---------------------------------------------------------------------------------------------
      
      # Run completeness check
      missing.parameters.result = check.gridsearch.completeness(as.matrix(log.liks[,3:6,with=FALSE]),
                                                                parameter.matrix)
      
      # In case some are missing we want to rerun and concatenate the missing values to 
      while (missing.parameters.result$All.complete == 0){
        
        missing.parameter.matrix = missing.parameters.result$Missing.parameters
        
        log.liks.update = aDDM.gridsearch(cur.choice.dat,
                                          cur.eye.dat,
                                          missing.parameter.matrix,
                                          cur.subject,
                                          cur.set_size,
                                          output.type,
                                          nr.reps,
                                          model.type,
                                          fixation.model,
                                          cur.log.file,
                                          timestep.ms,
                                          generate)
        
        log.liks = rbind(log.liks.log.liks.update)
        
        missing.parameters.result = check.gridsearch.completeness(as.matrix(log.liks[,3:6,with=FALSE]),
                                                                  parameter.matrix)
      }
      # ---------------------------------------------------------------------------------------------
      
      # We add two columns to log.liks which will be useful later when graphing results
      # ---------------------------------------------------------------------------------------------
      log.liks$Coarse = 1
      log.liks$Extension = 0
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      
      
      # CHECKING FOR CORNER SOLUTION AND EXTEND GRID SEARCH UNTIL NOT CORNER SOLUTION ANYMORE
      ###############################################################################################
      # ---------------------------------------------------------------------------------------------
      if (allow.extension == 1){
        corner = 1
        
        while (corner == 1){
          corner.results = check.and.update.corners(log.liks,
                                                    drifts,
                                                    thetas,
                                                    sds,
                                                    non.decision.times,
                                                    drift.shift,
                                                    drift.step,
                                                    sd.shift,
                                                    sd.step)
          
          if (corner.results$Corner == 0){
            corner = 0
          }
          
          if (corner.results$Corner == 1){
            
            # Run Gridsearch with new set of parameters
            new.parameter.grid = corner.results$New.parameters
            
            log.lik.extension = aDDM.gridsearch(cur.choice.dat,
                                                cur.eye.dat,
                                                new.parameter.grid,
                                                cur.subject,
                                                cur.set_size,
                                                output.type,
                                                nr.reps,
                                                model.type,
                                                fixation.model,
                                                cur.log.file,
                                                timestep.ms,
                                                generate)
            
            # CHECK WHETHER GRID SEARCH WAS COMPLETE AND RUN MISSING PARAMETER COMBINATIONS IF NOT
            # ---------------------------------------------------------------------------------------------
            
            # Run completeness check
            missing.parameters.result = check.gridsearch.completeness(as.matrix(log.lik.extension[,3:6,with=FALSE]),new.parameter.grid)
            
            # In case some are missing we want to rerun and concatenate the missing values to 
            while (missing.parameters.result$All.complete == 0){
              
              missing.parameter.matrix = missing.parameters.result$Missing.parameters
              
              log.lik.update = aDDM.gridsearch(cur.choice.dat,
                                               cur.eye.dat,
                                               missing.parameter.matrix,
                                               cur.subject,
                                               cur.set_size,
                                               output.type,
                                               nr.reps,
                                               model.type,
                                               fixation.model,
                                               timestep.ms,
                                               generate)
              
              log.lik.extension = rbind(log.lik.update,log.lik.extension)
        
              missing.parameters.result = check.gridsearch.completeness(as.matrix(log.lik.extension[,3:6,with=FALSE]),
                                                                        parameter.matrix)
            }
            # ---------------------------------------------------------------------------------------------
            
            # We add the same two columns to log.lok.extension as done to log.lik above 
            # ---------------------------------------------------------------------------------------------
            log.lik.extension$Coarse = 1
            log.lik.extension$Extension = 1
            # ---------------------------------------------------------------------------------------------
            
            log.liks = rbind(log.liks,log.lik.extension)
            
            # We define "cur.log.liks" because we check for corner solution within the extended parameters for each round of checks after round one
            # In round one we check over all precedingly tested parameter compositions (log.liks)
          }
        }
      }
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      
      # IF WE ALLOW A FINE GRID SEARCH THEN WE DO THE FOLLOWING
      ###############################################################################################
      
      if (allow.fine.grid == 1){
        
        # GENERATE FINE GRID 
        ###############################################################################################
        # ---------------------------------------------------------------------------------------------
        fine.parameter.matrix = generate.fine.grid(drift.step.fine,
                                                   sd.step.fine,
                                                   thetas,
                                                   non.decision.time.step.fine,
                                                   coarse.to.fine.ratio,
                                                   log.liks)
        # ---------------------------------------------------------------------------------------------
        ###############################################################################################
        
        
        # RUN FINE GRID SEARCH
        ###############################################################################################
        # ---------------------------------------------------------------------------------------------
        fine.log.liks = aDDM.gridsearch(cur.choice.dat,
                                        cur.eye.dat,
                                        fine.parameter.matrix,
                                        cur.subject,
                                        cur.set_size,
                                        output.type,
                                        nr.reps,
                                        model.type,
                                        fixation.model,
                                        cur.log.file,
                                        timestep.ms,
                                        generate)
        # ---------------------------------------------------------------------------------------------
        ###############################################################################################
        
        # CHECK WHETHER GRID SEARCH WAS COMPLETE
        ###############################################################################################
        # ---------------------------------------------------------------------------------------------
        
        # Run completeness check
        missing.parameters.result = check.gridsearch.completeness(as.matrix(fine.log.liks[,3:6,with=FALSE]),
                                                                  fine.parameter.matrix)
        
        # In case some are missing we want to rerun and concatenate the missing values to 
        while (missing.parameters.result$All.complete == 0){
          
          missing.parameter.matrix = missing.parameters.result$Missing.parameters
          
          added.log.liks = aDDM.gridsearch(cur.choice.dat,
                                           cur.eye.dat,
                                           missing.parameter.matrix,
                                           cur.subject,
                                           cur.set_size,
                                           output.type,
                                           nr.reps,
                                           model.type,
                                           fixation.model,
                                           cur.log.file,
                                           timestep.ms,
                                           generate)
        
          fine.log.liks = rbind(fine.log.liks,added.log.liks)
          
          missing.parameters.result = check.gridsearch.completeness(as.matrix(fine.log.liks[,3:6,with=FALSE]),
                                                                    parameter.matrix)
        }
        
        # We add two columns to fine.log.liks which will be useful later when graphing results (as done before)
        # ---------------------------------------------------------------------------------------------
        fine.log.liks$Coarse = 0
        fine.log.liks$Extension = 0
        
        # Concatenate fine and coarse loglikelihoods first
        log.liks = rbind(log.liks,fine.log.liks)
        # ---------------------------------------------------------------------------------------------
      }
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
      ###############################################################################################
      
      # SAVE STUFF TO FILE
      ###############################################################################################
      # Now we store all log.likelihood outcomes in file
      # ---------------------------------------------------------------------------------------------
      
      if (cur.subject == -1){
        cur.file = paste('Sim_logs/model_testing/loglik_modeltest_',toString(cur.set_size),sep="")
      }
      
      if (cur.subject == 0){
        cur.file = paste('Sim_logs/loglik_',model.type,'_all_setsize_', toString(cur.set_size), sep='') 
      }
      
      if (cur.subject > 0){
        cur.file = paste('Sim_logs/loglik_',model.type,'_subject_', toString(cur.subject), '_setsize_', toString(cur.set_size), sep='')
      }
      
      cur.file = paste(cur.file,".txt",sep="")
      write.table(log.liks,cur.file,quote=FALSE,sep=" ", col.names=TRUE,row.names=FALSE)
      
      #Status 
      print(cur.file)
      # ---------------------------------------------------------------------------------------------
      ###############################################################################################
    }
  }
}