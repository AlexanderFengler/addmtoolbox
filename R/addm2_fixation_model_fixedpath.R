#' Compute seqeuence of fixation locations and durations (fixed locations)
#' \code{addm_fixation_model_fixedpath} Returns matrix with row one providing fixation location and row two providing fixation duration
#' #' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Fixation model example (for usage with aevacc functions)
#' @return Returns matrix with row one providing fixation location and row two providing fixation duration
#' @export
#' @param fixdur.vec Vector of fixation-durations to sample from

addm2_fixation_model_fixedpath = function(){

  # Get fixation durations and location for simulation run -------------------------------------
  x = matrix(rep(0,500),ncol=2)
  x[,1] = rep(c(1,2),125)
  x[,2] = rep(400,250)
  # --------------------------------------------------------------------------------------------
  return(x)
}
