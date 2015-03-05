#' Run grid search over supplied parameter space
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Fit addm using grid search
#' @return data.table with log likelihoods by parameter combination
#' \code{addm_fit_grid}
#' @export
#' @param data list of three data.tables of each: choice data, eyetracking data, conditions data (as created by addm_dataprep)
#' @param drifts vector of all driftrate values to be tested.
#' @param thetas vector of all theta to be tested (between 0 and 1).
#' @param sds vector of all standard deviation values to be tested.
#' @param non.decision.times vector of all non decision times to be tested (in ms).
#' @param timestep integer that provides the timestep-size that is used in the simulations (in ms).
#' @param nr.reps integer that tells the function how many simulation runs to use.
#' @param model.type string that indicates which version of the model to run. 'standard' for normal model fits. 'memnoise' to allow for memory effects (see vignette for more for detailed explanation of what this is about).
#' @param fit.type string indicating either 'condition' for fits by unique trial conditions, 'trial' for fits by trial, or 'dyn' where you can use a dynamic programming algorithm for fitting the two items case, bypassing simulations for the fitting procedure
#' @param fixation.model string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (example). 'fixedpath' for following a predetermined fixation path with fixed durations (example). 'user' to provide your own fixation model, defined in a function "user_fixation_model" in the global environment.
#' @param allow.fine.grid variable that indicates whether we allow (1) a fine grid to be created and searched around the coarse grid minimum or not (0).
#' @param log.file path to a file for storing fit-logs
#' @param parallel boolean varible that indicates whether to initialize local cluster on start (1) or not (0).
#' @param coarse.to.fine.ratio integer defining the ratio between parameter steps in the coarse versus the fine grid.
#' @param state.step parameter only relevant when using fit.type = 'dyn', for which case it given the precision of the vertical grid utilized in the dynammic programming algorithm

addm_fit_grid = function(data = list(choice.dat = NULL, eye.dat = NULL, conditions.dat = NULL),
                         drifts = seq(0.001,0.0025,0.0005),
                         thetas = seq(0.0,1,0.25),
                         sds = seq(0.05,0.09,0.01),
                         non.decision.times = 0,
                         nr.reps = 1000,
                         timestep = 10,
                         model.type = 'standard',
                         fixation.model = 'fixedpath',
                         fit.type = 'condition',
                         allow.fine.grid = 0,
                         coarse.to.fine.ratio = 4,
                         log.file = "defaultlog.txt",
                         parallel = 1,
                         state.step = 0.1){

  # INITIALIZIONS --------------------------------------------------------------------------------
  choice.dat = data$choice.dat
  eye.dat = data$eye.dat
  conditions.dat = data$conditions.dat

  if (parallel == 1){
    cores = parallel::detectCores()
    doMC::registerDoMC(cores)
  }

  # Initialize the log.file we are going to use for storing intermediate results from the gridsearch
  cur.log.file = "cur_addm_log.txt"
  writeLines(c(""),cur.log.file)

  # Parameters for fine grid search
  # Distance between two points on the grid by dimension (drift.rate,sd)
  drift.step.fine = (drifts[length(drifts)] - drifts[length(drifts)-1])/coarse.to.fine.ratio
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
      log.liks = addm_support_gridsearch_foreach(choice.dat,
                                                 eye.dat,
                                                 conditions.dat,
                                                 parameter.matrix,
                                                 nr.reps,
                                                 timestep,
                                                 model.type,
                                                 fixation.model,
                                                 fit.type,
                                                 cur.log.file,
                                                 state.step)
      # ---------------------------------------------------------------------------------------------

      # UPDATE LOGLIKS ------------------------------------------------------------------------------
      log.liks$coarse = 1
      # ---------------------------------------------------------------------------------------------

      # IN CASE WE ALLOW FOR FINE GRID SEARCH -------------------------------------------------------
      if (allow.fine.grid == 1){

        # GENERATE FINE GRID ------------------------------------------------------------------------
        fine.parameter.matrix = addm_support_compute_finegrid(drift.step.fine,
                                                              sd.step.fine,
                                                              thetas,
                                                              non.decision.time.step.fine,
                                                              coarse.to.fine.ratio,
                                                              log.liks)
        # ---------------------------------------------------------------------------------------------

        # RUN FINE GRID SEARCH ------------------------------------------------------------------------
        fine.log.liks = addm_support_gridsearch_foreach(choice.dat,
                                                        eye.dat,
                                                        conditions.dat,
                                                        fine.parameter.matrix,
                                                        nr.reps,
                                                        timestep,
                                                        model.type,
                                                        fixation.model,
                                                        fit.type,
                                                        cur.log.file,
                                                        state.step)
        # ---------------------------------------------------------------------------------------------

        #
        # UPDATE LOGLIKS ------------------------------------------------------------------------------
        fine.log.liks$coarse = 0
        log.liks = rbind(log.liks,fine.log.liks)
        # ---------------------------------------------------------------------------------------------
      }

      # SAVE LOGLIKS ----------------------------------------------------------------------------------
      write.table(log.liks,log.file,quote=FALSE, sep=" ", col.names=TRUE, row.names=FALSE)
      # -----------------------------------------------------------------------------------------------

      # PRINT FINAL STATUS AND SET OF BEST PARAMETERS -------------------------------------------------
      writeLines(paste("Model was fit successfully and logs are saved in: ",log.file,sep=''))

      setkey(log.liks,loglik)
      writeLines(paste(' \nOptimal Parameters... \n \n',
                       'Drift Rate: ',
                       toString(log.liks[1,drift]),'\n',
                       'Theta: ',toString(log.liks[1,theta]),'\n',
                       'SD: ', toString(log.liks[1,sd]), '\n',
                       'Non decision time: ', toString(log.liks[1,non.decision.time]),sep=''))
      # -----------------------------------------------------------------------------------------------
      file.remove('cur_addm_log.txt')
      return(log.liks)
      }
