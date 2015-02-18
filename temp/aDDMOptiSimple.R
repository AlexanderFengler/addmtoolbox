# Author: Alexander Fengler
# Date: Jan 16th 2015
# Purpose: Main function running aDDM optimization based on grid search algorithm

# Several options: run by subject, run with fake dataset, run over all subjects

aDDMOpti = function(subjects,                    # Subject to be tested: -1 = Model testing with fake data // 0 = Fit over all subjects // > 0 fit by subject
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
                   fixation.model,              # Fixation model used in simulations: "Normal" - Real Fixations used, "Random", "FakePath" - prespecified path for fixations
                   allow.extension,             # Binary, whether we allow extensions of the grid in case we find corner-solution
                   allow.fine.grid,             # Binary, whether we allow the fine-grid step in the grid search procedure
                   generate){                   # Generate (binary) tells the aDDM function whether to spit out full data (1) or log likelihoods (0)

  # LOAD ALL PACKAGES AND FUNCTIONS NECESSARY ----------------------------------------------------
  # Load all high level functions needed
  source('temp/addm_opti_supportfuns/generate_parameter_combinations.R')
  source('temp/addm_opti_supportfuns/aDDM_gridsearch.R')
  source('temp/addm_opti_supportfuns/check_grid_search_completeness.R')
  source('temp/addm_opti_supportfuns/check_and_update_corners.R')
  source('temp/addm_opti_supportfuns/generate_fine_grid.R')


  source('temp/load_packages_scripts_aDDM.R')
  # Initialize all packages
  load.aDDM.packages()

  # Initialize parameters, subjects, set_size
  load.aDDM.scripts()

  # INITIALIZING VARIABLES THAT CAN BE MANIPULATED ------------------------------------------------

  # Initialize the logfile we are going to use for storing results from the gridsearch
  cur.log.file = "temp/cur_log.txt"
  writeLines(c(""),cur.log.file)

  # In case we find corner solutions what parameter shift shall we allow per round (balance between degree of exploration and time spend computing)
  drift.shift = (drifts[length(drifts)] - drifts[length(drifts)-1])*3
  sd.shift = (sds[length(sds)] - sds[length(sds)-1])*3
  non.decision.time.shift = (non.decision.times[length(non.decision.times)] -
                               non.decision.times[length(non.decision.times)-1])*3

  # Parameters for fine grid search
  # Distance between two points on the grid by dimension (drift.rate,sd)
  coarse.to.fine.ratio = 5

  drift.step.fine = (drifts[length(drifts)] - drifts[length(drifts)-1])/coarse.to.fine.ratio
  theta.step.fine = (thetas[length(thetas)] - thetas[length(thetas)-1])/coarse.to.fine.ratio
  sd.step.fine = (sds[length(sds)] - sds[length(sds)-1])/coarse.to.fine.ratio
  non.decision.time.step.fine = (non.decision.times[length(non.decision.times)] -
                                   non.decision.times[length(non.decision.times)-1])/coarse.to.fine.ratio
  # -----------------------------------------------------------------------------------------------

  for (cur.subject in subjects){
    for (cur.set_size in set.sizes){
      #
      # GENERATING SUBSET OF DATA.FRAMES ----------------------------------------------------------

      # EXTRACT RELEVANT VALUATIONS
        if (cur.subject == -1){
          cur.choice.dat = core.clean.model.test.dat$choice.dat
        } else if (cur.subject == 0){
          cur.choice.dat =  choice.dat[choice.dat$set_size == cur.set_size, names(choice.dat),with=FALSE]
        } else {
          cur.choice.dat =  choice.dat[choice.dat$subject == cur.subject &
                                         choice.dat$set_size == cur.set_size, names(choice.dat),with=FALSE]
        }
      }

      # EXTRACT RELEVANT FIXATIONS
      if (cur.subject == -1) {
        cur.eye.dat = core.clean.mode.test.dat$fixations
      } else if (cur.subject == 0){
        cur.eye.dat = eye.dat[eye.dat$set_size == cur.set_size, list(fixloc,fixnr,fixdur,trialid)]
      } else {
        cur.eye.dat = eye.dat[eye.dat$subject == cur.subject &
                                eye.dat$set_size == cur.set_size, list(fixloc,fixnr,fixdur,trialid)]
      }

      # Set keys for the data.tables passed into function
      setkey(cur.choice.dat,trialid)
      setkey(cur.eye.dat,trialid)
      # ---------------------------------------------------------------------------------------------

      #
      # GENERATE PARAMETER MATRIX -------------------------------------------------------------------
      parameter.matrix = param.combs(drifts,
                                     sds,
                                     thetas,
                                     non.decision.times)
      # ---------------------------------------------------------------------------------------------

      # RUN COARSE GRID SEARCH ----------------------------------------------------------------------
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

      # UPDATE LOGLIKS ------------------------------------------------------------------------------
      log.liks$Coarse = 1
      log.liks$Extension = 0
      # ---------------------------------------------------------------------------------------------

      # IN CASE WE ALLOW FOR FINE GRID SEARCH -------------------------------------------------------
      if (allow.fine.grid == 1){

        # GENERATE FINE GRID ------------------------------------------------------------------------
        fine.parameter.matrix = generate.fine.grid(drift.step.fine,
                                                   sd.step.fine,
                                                   thetas,
                                                   non.decision.time.step.fine,
                                                   coarse.to.fine.ratio,
                                                   log.liks)
        # ---------------------------------------------------------------------------------------------

        # RUN FINE GRID SEARCH ------------------------------------------------------------------------
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

        #
        # UPDATE LOGLIKS ------------------------------------------------------------------------------
        fine.log.liks$Coarse = 0
        fine.log.liks$Extension = 0
        log.liks = rbind(log.liks,fine.log.liks)
        # ---------------------------------------------------------------------------------------------
      }

      # SAVE LOGLIKS ----------------------------------------------------------------------------------
      if (cur.subject == -1){
        cur.file = paste('temp/model_testing/loglik_modeltest_',toString(cur.set_size), sep="")
      }

      if (cur.subject == 0){
        cur.file = paste('temp/loglik_',model.type,'_all_setsize_', toString(cur.set_size), sep="")
      }

      if (cur.subject > 0){
        cur.file = paste('temp/loglik_',model.type,'_subject_', toString(cur.subject), '_setsize_', toString(cur.set_size), sep="")
      }

      cur.file = paste(cur.file,".txt",sep="")
      write.table(log.liks,cur.file,quote=FALSE,sep=" ", col.names=TRUE, row.names=FALSE)

      #Status
      print(cur.file)
      # ---------------------------------------------------------------------------------------------
    }
  }
}
