#' Compute seqeuence of fixation locations and durations (random locations)
#' \code{addm2_fixation_model_random}
#'  @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Example for user-supplied fixation model (aevacc functions)
#' @return matrix with two rows. One providing fixation location and the second providing fixation duration
#' @export
#' @param cur.valuations vector that provides item valuations in current choice set
#' @param fixdur integer variable providing the duration by fixation simulated (in ms)

addm2_fixation_model_random = function(cur.valuations = c(0,0), fixdur = 400){

  # Get fixation durations and location for simulation run -------------------------------------
  x = matrix(rep(0,500),ncol=2)
  locsamp = sample(c(1,2),2,replace=FALSE)
  x[,1] = rep(locsamp,125)
  x[,2] = rep(fixdur,250)
  # --------------------------------------------------------------------------------------------
  return(x)
}
