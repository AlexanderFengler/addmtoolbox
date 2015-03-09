#' Compute parameter space for fine grid given the optimal parameter set of a coarse grid search
#' \code{addm_support_compute_finegrid}
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Compute fine grid best parameters from coarse grid search (internal)
#' @return matrix storing new parameter space
#' @export
#' @param drift.step.fine numeric variable storing precision with which drift rate space is searched in the fine grid parameter matrix.
#' @param sd.step.fine numeric variable storing precision with which standard deviation space is searched in the fine grid parameter matrix.
#' @param thetas numeric variable storing testable thetas. The theta space is not changed in the fine parameter matrix.
#' @param non.decision.time.step.fine numeric variable storing precision with which non.decision.time space is searched in the fine grid parameter matrix.
#' @param coarse.to.fine.ratio numeric Variable that stores the precision ratio between the coarse and fine grid for grid search fits.
#' @param log.liks data.table that stores log likelihoods and corresponding parameter combinations resulting from coarse grid search.


addm_support_compute_finegrid = function(drift.step.fine = 0,
                                         thetas = 0,
                                         gammas = 0,
                                         sd.step.fine = 0,
                                         non.decision.time.step.fine = 0,
                                         coarse.to.fine.ratio = 0,
                                         nr.attributes = 1,
                                         log.liks = data.table(drift = 0, theta = 0, sd = 0, non.decision.time = 0, nr.reps = 0, loglik = 0)){

  # Get parameter values at minimum log likelihood -------------------------------------------------------------------
  max.val = max(log.liks$Log.Likelihood)
  max.pos = which(log.liks$Log.Likelihood == max.val)[sample(length(which(log.liks$Log.Likelihood == max.val)),1)]
  center.parameters = as.numeric(log.liks[max.pos,])
  center.drift = center.parameters[1]
  center.sd = center.parameters[3]
  center.non.decision.time = center.parameters[4]
  # ------------------------------------------------------------------------------------------------------------------

  # Generate new parameter space -------------------------------------------------------------------------------------
  fine.non.decision.times = center.non.decision.time
  len.non.decision.times = unique(log.liks$non.decision.time)

  if (len.non.decision.times > 1){
    fine.non.decision.times =  seq(center.non.decision.time - coarse.to.fine.ratio*non.decision.time.step.fine,
                                   center.non.decision.time + coarse.to.fine.ratio*non.decision.time.step.fine,
                                   non.decision.time.step.fine)
  }

  fine.drifts = seq(center.drift-(coarse.to.fine.ratio-1)*drift.step.fine,
                    center.drift+(coarse.to.fine.ratio-1)*drift.step.fine,
                    drift.step.fine)

  fine.thetas = thetas

  fine.sds = seq(center.sd-(coarse.to.fine.ratio-1)*sd.step.fine,
                 center.sd+(coarse.to.fine.ratio-1)*sd.step.fine*2,
                 sd.step.fine)

  fine.gammas = gammas
  # ------------------------------------------------------------------------------------------------------------------

  # Make and return parameter matrix ---------------------------------------------------------------------------------
  if (nr.attributes == 1){
    fine.grid = as.matrix(expand.grid(fine.drifts,
                                      fine.thetas,
                                      fine.sds,
                                      fine.non.decision.times))
  } else {
    fine.grid = as.matrix(expand.grid(fine.drifts,
                                      fine.thetas,
                                      fine.gammas,
                                      fine.sds,
                                      fine.non.decision.times))
  }
  return(fine.grid)
  # ------------------------------------------------------------------------------------------------------------------
}
