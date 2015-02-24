# Author: Alexander Fengler
# Date: Jan 16th 2015
# Purpose: Main function running aDDM optimization based on grid search algorithm

# Several options: run by subject, run with fake dataset, run over all subjects

addm_fit_grid = function(conditions = data.table(v1 = c(1,2,3),v2 = c(3,2,1),id = c(1,2,3)), # Choice data - Needs: All item values per condition as seperate columns by item, item chosen,subject, set.sizes, id, RT
                         eye.dat = data.table(fixloc = 0,fixnr = 0, fixdur= 0, id = c(1,2,3)),                     # Fixation data - Needs: Fixation Locations, Fixation Durations, subject, set.sizes, id (sorted within trial according to fixation number)
                         choice.dat = data.table(v1 = c(1,2,3),v2 = c(3,2,1),id = c(1,2,3), rt = c(0,0,0), decision = c(1,1,1)),
                         drifts = seq(0.002,0.01,0.002),                      # Minimum Drift Rate we consider in grid
                         thetas = seq(0.2,1,0.2),                      # Minimum Theta we consider in grid
                         sds = seq(0.02,0.1,0.02),                         # Minimum SD we consider in grid
                         non.decision.times = 0,          # Minimum Non Decision Time
                         timestep.ms = 10,                 # Timestep-size in milliseconds for aDDM simulation propagation
                         nr.reps = 2000,                     # Number of repitions used in simulation runs
                         model.type = 'nomem',                  # Currently choice between model with memory effects, divided in no drift for items not seen ('mem') and no drift and no sd for items not seen ('memzeronoise') and without ('nomem')
                         output.type = 'condition',                 # Output type can be: "Opti" - Normal loglikelihood for model optimization, "FakeOpti" - loglikelihoods for optimization for artificially created data
                         fixation.model = 'FakePath',              # Fixation model used in simulations: "Normal" - Real Fixations used, "Random", "FakePath" - prespecified path for fixations
                         allow.fine.grid = 0,             # Binary, whether we allow the fine-grid step in the grid search procedure
                         generate = 0,
                         logfile = "defaultlog.txt"){                   # Generate (binary) tells the aDDM function whether to spit out full data (1) or log likelihoods (0)

  # LOAD ALL PACKAGES AND FUNCTIONS NECESSARY ----------------------------------------------------
  # Load all high level functions needed
  source('temp/addm_opti_supportfuns/addm_support_compute_finegrid.R')
  source('temp/addm_opti_supportfuns/addm_support_gridsearch_foreach.R')
  source('temp/addm_support_scripts_packages.R')

  # Initialize all packages
  load.aDDM.packages()

  # Initialize parameters, subjects, set_size
  load.aDDM.scripts()

  # INITIALIZING VARIABLES THAT CAN BE MANIPULATED ------------------------------------------------

  # Initialize the logfile we are going to use for storing results from the gridsearch
  cur.log.file = "temp/cur_addm_log.txt"
  writeLines(c(""),cur.log.file)

  # Parameters for fine grid search
  # Distance between two points on the grid by dimension (drift.rate,sd)
  coarse.to.fine.ratio = 5
  drift.step.fine = (drifts[length(drifts)] - drifts[length(drifts)-1])/coarse.to.fine.ratio
  theta.step.fine = (thetas[length(thetas)] - thetas[length(thetas)-1])/coarse.to.fine.ratio
  sd.step.fine = (sds[length(sds)] - sds[length(sds)-1])/coarse.to.fine.ratio
  non.decision.time.step.fine = (non.decision.times[length(non.decision.times)] -
                                   non.decision.times[length(non.decision.times)-1])/coarse.to.fine.ratio
  # -----------------------------------------------------------------------------------------------

      # Make choice.dat and eye.dat data.tables, in case they are not already
      choice.dat = as.data.table(choice.dat)
      eye.dat = as.data.table(eye.dat)
      # Set keys for the data.tables passed into function
      setkey(choice.dat,id)
      setkey(eye.dat,id)
      # ---------------------------------------------------------------------------------------------

      # GENERATE PARAMETER MATRIX -------------------------------------------------------------------
      parameter.matrix = as.matrix(expand.grid(drifts,thetas,sds,non.decision.times))
      # ---------------------------------------------------------------------------------------------

      # RUN coarse GRID SEARCH ----------------------------------------------------------------------
      log.liks = addm_gridsearch_foreach(conditions,
                                         eye.dat,
                                         choice.dat,
                                         parameter.matrix,
                                         output.type,
                                         nr.reps,
                                         model.type,
                                         fixation.model,
                                         cur.log.file,
                                         timestep.ms,
                                         generate)
      # ---------------------------------------------------------------------------------------------

      # UPDATE LOGLIKS ------------------------------------------------------------------------------
      log.liks$coarse = 1
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
        fine.log.liks = addm_gridsearch_foreach(conditions,
                                                eye.dat,
                                                choice.dat,
                                                fine.parameter.matrix,
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
        fine.log.liks$coarse = 0
        log.liks = rbind(log.liks,fine.log.liks)
        # ---------------------------------------------------------------------------------------------
      }

      # SAVE LOGLIKS ----------------------------------------------------------------------------------
      write.table(log.liks,logfile,quote=FALSE, sep=" ", col.names=TRUE, row.names=FALSE)

      #Status
      #print(paste(logfile,"done....",sep=': '))
      # ---------------------------------------------------------------------------------------------
      return(log.liks)
      }
