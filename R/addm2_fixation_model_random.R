#' Compute seqeuence of fixation locations and durations (random locations)
#' \code{addm_fixation_model_random} Returns matrix with row one providing fixation location and row two providing fixation duration
#' #' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Fixation model example (for usage with aevacc functions)
#' @return Returns matrix with row one providing fixation location and row two providing fixation duration
#' @export
#' @param fixdur.vec Vector of fixation-durations to sample from

addm2_fixation_model_random = function(fixdur.vec = c(400,400)){

  # Get fixation durations and location for simulation run -------------------------------------
  x = matrix(rep(0,500),ncol=2)
  locsamp = sample(c(1,2),2,replace=FALSE)
  x[,1] = rep(locsamp,125)
  x[,2] = sample(fixdur.vec,250,replace=TRUE)
  # --------------------------------------------------------------------------------------------
  return(x)
}
