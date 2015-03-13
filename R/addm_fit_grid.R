#' Run grid search over supplied parameter space
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Fit addm using grid search
#' @return data.table with log likelihoods by parameter combination
#' \code{addm_fit_grid}
#' @export
#' @param data list of three data.tables of each: choice data, eyetracking data, conditions data (as created by addm_dataprep)
#' @param drifts vector of all driftrate values to be tested.
#' @param thetas vector of all theta values to be tested [0,1].
#' @param gammas vector of all gamma values to be tested [0,1] (matters only when supplying data with multiple attributes by item)
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
#' @param boundaryfun function that is supplied by user for the decision boundaries (has to have at least two inputs: maxrt, timestep)
#' @param boundary.parameters matrix or vector that provides a parameter-space for all parameter sets that shall be tested on the boundary function

addm_fit_grid = function(data = list(choice.dat = NULL, eye.dat = NULL, conditions.dat = NULL, attributes = NULL),
                         drifts = seq(0.001,0.0025,0.0005),
                         thetas = seq(0.0,1,0.25),
                         gammas = 1,
                         sds = seq(0.05,0.09,0.01),
                         non.decision.times = 0,
                         scalar_item_not_seen_drift = 1,
                         scalar_item_not_seen_noise = 1,
                         boundaryfun = 1,
                         boundary.parameters = 0,
                         nr.reps = 1000,
                         timestep = 10,
                         model.type = 'standard',
                         fixation.model = 'fixedpath',
                         fit.type = 'condition',
                         allow.fine.grid = 0,
                         coarse.to.fine.ratio = 4,
                         log.file = 'defaultlog.txt',
                         parallel = 1,
                         state.step = 0.1){

  # INITIALIZIONS --------------------------------------------------------------------------------
  choice.dat = data$choice.dat
  eye.dat = data$eye.dat
  conditions.dat = data$conditions.dat
  nr.attributes = data$attributes

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
      #return(list(choice = choice.dat, eye = eye.dat))
      # GENERATE PARAMETER MATRIX -------------------------------------------------------------------
      parameter.matrix = as.matrix(expand.grid(non.decision.times, drifts, sds, thetas, gammas, scalar_item_not_seen_drift, scalar_item_not_seen_noise))

      # If we have a function supplied for the boundary then we will complement the parameter.matrix
      if (class(boundaryfun) == 'function'){
        parameter.matrix = as.matrix(do.call('expand.grid', c(as.data.frame(parameter.matrix), as.data.frame(boundary.parameters))))
      }
      # ---------------------------------------------------------------------------------------------

      # RUN coarse GRID SEARCH ----------------------------------------------------------------------
      log.liks = addm_support_gridsearch_foreach(choice.dat,
                                                 eye.dat,
                                                 conditions.dat,
                                                 parameter.matrix,
                                                 nr.attributes,
                                                 boundaryfun,
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
                                                              thetas,
                                                              gammas,
                                                              sd.step.fine,
                                                              non.decision.time.step.fine,
                                                              coarse.to.fine.ratio,
                                                              nr.attributes,
                                                              boundary.parameters,
                                                              log.liks)
        # ---------------------------------------------------------------------------------------------

        # RUN FINE GRID SEARCH ------------------------------------------------------------------------
        fine.log.liks = addm_support_gridsearch_foreach(choice.dat,
                                                        eye.dat,
                                                        conditions.dat,
                                                        fine.parameter.matrix,
                                                        nr.attributes,
                                                        boundaryfun,
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
      if (model.type == 'standard'){
        if (nr.attributes == 1){
          writeLines(paste(' \n Optimal Parameters... \n \n',
                           'Drift Rate: ', toString(log.liks[1,drift]),'\n',
                           'Theta: ',toString(log.liks[1,theta]),'\n',
                           'SD: ', toString(log.liks[1,sd]), '\n',
                           'Non decision time: ', toString(log.liks[1,non.decision.time]),sep=''))
        } else {
          writeLines(paste(' \n Optimal Parameters... \n \n',
                           'Drift Rate: ', toString(log.liks[1,drift]),'\n',
                           'Theta: ',toString(log.liks[1,theta]),'\n',
                           'Gamma: ',toString(log.liks[1,gamma]), '\n',
                           'SD: ', toString(log.liks[1,sd]), '\n',
                           'Non decision time: ', toString(log.liks[1,non.decision.time]),sep=''))
        }
      }

      if (model.type == 'memnoise'){
        writeLines(paste(' \n Optimal Parameters... \n \n',
                         'Drift Rate: ', toString(log.liks[1,drift]),'\n',
                         'Theta: ',toString(log.liks[1,theta]),'\n',
                         'SD: ', toString(log.liks[1,sd]), '\n',
                         'Non decision time: ', toString(log.liks[1,non.decision.time]),
                         'Scalar item not seen - drift: ', toString(log.liks[1,scalar_item_not_seen_drift]),
                         'Scalar item not seen - noise: ', toString(log.liks[1,scalar_item_not_seen_noise]), sep=''))
      }
      # -----------------------------------------------------------------------------------------------
      file.remove('cur_addm_log.txt')
      return(log.liks)
      }
