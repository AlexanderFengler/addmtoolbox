#' Function internal to grid search, which computes a fine grid around the coarse grid minimum log likelihood
#' \code{addm_support_compute_finegrid} Gives back a parameter matrix
#' #' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Generate fine grid around best parameter combination
#' @return Matrix storing parameter space
#' @export
#' @param drift.step.fine Numeric variable storing precision with which drift rate space is searched in the fine grid parameter matrix.
#' @param sd.step.fine Numeric variable storing precision with which standard deviation space is searched in the fine grid parameter matrix.
#' @param thetas Numeric variable storing testable thetas. The theta space is not changed in the fine parameter matrix.
#' @param non.decision.time.step.fine Numeric variable storing precision with which non.decision.time space is searched in the fine grid parameter matrix.
#' @param coarse.to.fine.ratio Numeric Variable that stores the precision ratio between the coarse and fine grid for grid search fits.
#' @param log.liks 'data.table' that stores log likelihoods and corresponding parameter combinations resulting from coarse grid search.


addm_support_compute_finegrid = function(drift.step.fine,
                                         sd.step.fine,
                                         thetas,
                                         non.decision.time.step.fine,
                                         coarse.to.fine.ratio,
                                         log.liks){

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
  # ------------------------------------------------------------------------------------------------------------------

  # Make and return parameter matrix ---------------------------------------------------------------------------------
  fine.grid = as.matrix(expand.grid(fine.drifts,fine.sds,fine.thetas,fine.non.decision.times))
  return(fine.grid)
  # ------------------------------------------------------------------------------------------------------------------
}
