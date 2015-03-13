#' Compute seqeuence of fixation locations and durations (fixed locations)
#' \code{addm2_fixation_model_fixedpath}
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Example for user-supplied fixation model (aevacc functions)
#' @return matrix with two rows. One providing fixation location and the second providing fixation duration
#' @param cur.valuations vector of current valuations in choice set
#' @param fixdur integer variable that provides the duration of fixations simulated (in ms)
#' @export

addm2_fixation_model_fixedpath = function(cur.valuations = c(0,0), fixdur = 400){
  # Get fixation durations and location for simulation run -------------------------------------
  x = matrix(rep(0,500),ncol=2)
  x[,1] = rep(c(1,2),125)
  x[,2] = rep(fixdur,250)
  # --------------------------------------------------------------------------------------------
  return(x)
}
