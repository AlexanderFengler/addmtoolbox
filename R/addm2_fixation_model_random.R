#' Compute seqeuence of fixation locations and durations (random locations)
#' \code{addm2_fixation_model_random}
#'  @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Example for user-supplied fixation model (aevacc functions)
#' @return matrix with two rows. One providing fixation location and the second providing fixation duration
#' @export

addm2_fixation_model_random = function(){

  # Get fixation durations and location for simulation run -------------------------------------
  x = matrix(rep(0,500),ncol=2)
  locsamp = sample(c(1,2),2,replace=FALSE)
  x[,1] = rep(locsamp,125)
  x[,2] = rep(400,250)
  # --------------------------------------------------------------------------------------------
  return(x)
}
