#' Compute seqeuence of fixation locations and durations (fixed locations)
#' \code{addm_fixation_model_fixedpath}
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Example for user-supplied fixation model (aevacc functions)
#' @return matrix with two rows. One providing fixation location and the second providing fixation duration
#' @param cur.valuations vector of valuations in current choice set
#' @param set.size number of items present in current choice set
#' @export

addm_fixation_model_fixedpath = function(cur.valuations = 2, set.size = 4){

  # Get fixation durations and location for simulation run -------------------------------------
  reps = ceiling(250/set.size)
  nr.fix = reps * set.size
  x = matrix(rep(0,nr.fix*2),ncol=2)
  x[,1] = rep(c(1:set.size),reps)
  x[,2] = rep(400,nr.fix)
  # --------------------------------------------------------------------------------------------
  return(x)
}
