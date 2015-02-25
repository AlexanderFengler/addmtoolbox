#' Compute new fixation location and duration
#' \code{addm_fixation_model_example} Returns vector with [1] fixation location and [0] fixation duration
#' #' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Fixation model example (for usage with aevacc functions)
#' @return Returns vector with [1] fixation location and [0] fixation duration
#' @export
#' @param last.fixloc Numeric Variable providing last fixation location in simulation
#' @param fixdur.vec Vector of fixation-durations to sample from

addm_fixation_model_example = function(last.fixloc = 0,
                                       fixdur.vec = 0){

  # Get new fixation location ------------------------------------------------------------------
  if (last.fixloc == 1){
    fixloc = 1;
  } else {
    fixloc = 2;
  }
  # --------------------------------------------------------------------------------------------

  # Get new fixation duration ------------------------------------------------------------------
  fixdur = sample(fixdur.vec,1)
  # --------------------------------------------------------------------------------------------

  return(c(fixloc,fixdur))
}
