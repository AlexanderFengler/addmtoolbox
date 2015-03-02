#' Compute seqeuence of fixation locations and durations (fixed locations)
#' \code{addm2_fixation_model_fixedpath}
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Example for user-supplied fixation model (aevacc functions)
#' @return matrix with two rows. One providing fixation location and the second providing fixation duration
#' @export

addm2_fixation_model_fixedpath = function(){

  # Get fixation durations and location for simulation run -------------------------------------
  x = matrix(rep(0,500),ncol=2)
  x[,1] = rep(c(1,2),125)
  x[,2] = rep(400,250)
  # --------------------------------------------------------------------------------------------
  return(x)
}
