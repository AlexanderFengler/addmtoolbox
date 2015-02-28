#' Runs grid search over supplied parameter space
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Run Gridsearch over supplied parameter space
#' @return Returns data.table with log likelihoods by parameter combination
#' \code{addm_fit_grid} Returns data.table with log likelihoods by parameter combination
#' @export
#' @param eye.dat A 'data.frame' or 'data.table' storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and an id column (id). Can all be initialized as zero columns when a by condition fit is attempted.
#' @param choice.dat A 'data.frame' or  'data.table' storing the item valuations (v1,v2...) , reaction times in ms (rt), decisions as decision and an id column (id). A by trial form is assumed.
#' @param conditions.dat A 'data.frame' storing the item valuations (v1,v2...) by unique trial conditions. An id column (id) needs to be provided that matches by trial data.
#' @param drifts Vector of all driftrate values to be tested.
#' @param thetas Vector of all theta to be tested (between 0 and 1).
#' @param sds Vector of all standard deviation values to be tested.
#' @param non.decision.times Vector of all non decision times to be tested (in ms).
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param nr.reps An integer number that tells the function how many simulation runs to use.
#' @param model.type A string that indicates which version of the model to run. 'standard' or 'memnoise' when memory effects shall be allowed.
#' @param output.type A string that indicates what output the model shall produce. 'full' for detailed model output, 'fit' for sparse output (rt,decision) by id variable.
#' @param fit.type A string indicating either 'condition' for fits by unique trial conditions or 'trial' for fits by trial.
#' @param fixation.model A string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (implemented) 'fixedpath' for following a predetermined fixation path with fixed durations (implemented).
#' @param allow.fine.grid Binary variable that indicates whether we allow (1) a fine grid to be created an searched around the coarse grid minimum or not (0).
#' @param log.file Filepath for storage of logs
#' @param parallel Binary variable that tells whether to initialize local cluster on start (1) or not (0).

addm_fit_grid = function(conditions.dat = data.table(v1 = 0, v2 = 0, id = 0),
                         eye.dat = data.table(fixloc = 0,fixnr = 1, fixdur= 0, id = 0),
                         choice.dat = data.table(v1 = 0,v2 = 0, rt = 0, decision = 0, id = 0),
                         drifts = seq(0.0,0.003,0.0005),
                         thetas = seq(0.0,1,0.25),
                         sds = seq(0.05,0.15,0.025),
                         non.decision.times = 0,
                         timestep = 10,
                         nr.reps = 2000,
                         model.type = 'standard',
                         output.type = 'fit',
                         fit.type = 'condition',
                         fixation.model = 'fixedpath',
                         allow.fine.grid = 0,
                         log.file = "defaultlog.txt",
                         parallel = 1){

  # INITIALIZIONS --------------------------------------------------------------------------------
  if (parallel == 1){
    cores = detectCores()
    registerDoMC(cores)
  }

  # Initialize the log.file we are going to use for storing intermediate results from the gridsearch
  cur.log.file = "cur_addm_log.txt"
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
      log.liks = addm_support_gridsearch_foreach(conditions.dat,
                                                 eye.dat,
                                                 choice.dat,
                                                 parameter.matrix,
                                                 fit.type,
                                                 nr.reps,
                                                 model.type,
                                                 fixation.model,
                                                 log.file,
                                                 timestep)
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
        fine.log.liks = addm_support_gridsearch_foreach(conditions.dat,
                                                        eye.dat,
                                                        choice.dat,
                                                        fine.parameter.matrix,
                                                        fit.type,
                                                        nr.reps,
                                                        model.type,
                                                        fixation.model,
                                                        log.file = cur.log.file,
                                                        timestep)
        # ---------------------------------------------------------------------------------------------

        #
        # UPDATE LOGLIKS ------------------------------------------------------------------------------
        fine.log.liks$coarse = 0
        log.liks = rbind(log.liks,fine.log.liks)
        # ---------------------------------------------------------------------------------------------
      }

      # SAVE LOGLIKS ----------------------------------------------------------------------------------
      write.table(log.liks,log.file,quote=FALSE, sep=" ", col.names=TRUE, row.names=FALSE)

      #Status
      #print(paste(log.file,"done....",sep=': '))
      # ---------------------------------------------------------------------------------------------
      file.remove('cur_addm_log.txt')
      return(log.liks)
      }
