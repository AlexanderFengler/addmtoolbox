#' Compute seqeuence of fixation locations and durations (fixed locations)
#' \code{addm_fixation_model_random}
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Example for user-supplied fixation model (aevacc functions)
#' @return matrix with two rows. One providing fixation location and the second providing fixation duration
#' @param cur.valuations vector of valuations in current choice set
#' @param set.size number of items present in current choice set
#' @export

addm_fixation_model_random = function(cur.valuations = c(0,0), set.size = 2){

  # Get fixation durations and location for simulation run -------------------------------------
  reps = ceiling(250/set.size)
  nr.fix = reps * set.size
  x = matrix(rep(0,nr.fix*2),ncol=2)
  x[,1] = sample(set.size,nr.fix,replace=TRUE)
  x[,2] = rep(400,nr.fix)
  # --------------------------------------------------------------------------------------------
  return(x)
}
