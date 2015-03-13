#' Generate boundaries with convex decay, and potential decay delay
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Convex boundary
#' @return vector that provides boundary levels by timestep
#' \code{addm_boundary_convex()}
#' @export
#' @param maxrt integer providing the maximum rt that should be considered for the boundary
#' @param timestep integer that provides the timestep-size that the current fit utilized
#' @param time_decay_start time delay for start of decay (boundary parameter)
#' @param decay numeric variable providing value for decay concavity

addm_boundary_convex = function(maxrt = 3000, timestep = 10, time_decay_start = 1000, decay = 0.999){
  out = rep(0,ceiling(maxrt/timestep))
  cnt = 1
  time = 0

  while(time <= maxrt){
    while (time <= time_decay_start){
      out[cnt] = 1
      cnt[1] = cnt[1] + 1
      time[1] = time[1] + timestep
    }
    out[cnt] = 1 - (1 - decay^(time - time_decay_start))
    cnt[1] = cnt[1] + 1
    time[1] = time[1] + timestep
  }
  return(out)
}
