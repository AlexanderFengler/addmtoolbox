#' Visualize a single model run (two items)
#' \code{addm2_plot_modelvis()}
#' @param theta numeric variable indicating the attential bias to be use for the model run ([0,1])
#' @param v1 numeric variable indicating the value of the left item
#' @param v2 numeric variable indicating the value of the right item
#' @param sd numeric variable indicating the standard deviation of the noise term used in the model run
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Visualize a single drift diffusion model run with two items
#' @return ggplot object. Costumizable.
#' @export
addm2_plot_modelvis = function(d = 0.002,
                               theta=0.1,
                               sd = 0.02,
                               v1 = 2,
                               v2 = 4){

  # Run model and store essential output ----------------------------------------------------
  rdv = 0
  rdv.vec = seq(1:1000)

  cnt = 1
  for (fix.pos in rep(c(1,2,1,2),each=250)){

    if (fix.pos == 1){
      rdv[1] = rdv[1] + d*(v1 - theta*v2) + rnorm(1,0,sd)
    }
    if (fix.pos == 2){
      rdv[1] = rdv[1] + d*(theta*v1 - v2) + rnorm(1,0,sd)
    }

    rdv.vec[cnt] = rdv
    cnt = cnt + 1
  }

  out.mat = data.frame(rdv = rdv.vec,
                       time = seq(1:1000),
                       fixation.position = rep(c(1,2,1,2),each=250))
  pos = which.min(which(out.mat[,1] >= 1 | out.mat[,1] <= -1))
  decision.time = which(out.mat[,1] >= 1 | out.mat[,1] <= -1)[pos]
  out.mat$decision = 0
  out.mat$decision[decision.time] = 1
  # ------------------------------------------------------------------------------------------

  # Get a matrix that provides text annotation positions ------------------------------------
  x.max = decision.time
  nr.fix = ceiling(decision.time/250)
  last.fix = ((decision.time/250) - floor(decision.time/250))*250

  x = matrix(rep(4000,8),ncol=2)
  for(i in 1:nr.fix){
    if (i < nr.fix | nr.fix == 1){
      if (i == 1){
        x[i,1] = 1
        if (decision.time >= 250){
          x[i,2] = i*250
        } else {
          x[i,2] = decision.time - 1
        }
      } else{
        x[i,1] = (i-1)*(250)
        x[i,2] = i*250
      }
    } else {
      x[i,2] = 250*(i-1) + last.fix
      x[i,1] = 250*(i-1)
    }
  }
  x = as.data.frame(x)
  # ------------------------------------------------------------------------------------------


  # Draw Plot --------------------------------------------------------------------------------
  addm.example.plot = ggplot(out.mat, aes(x=time,y=rdv)) +
    annotate("rect",xmin=x[1,1], xmax=x[1,2], ymin=-1, ymax=1, alpha=0.3, fill = 'grey20') +
    annotate("rect",xmin=x[2,1], xmax=x[2,2], ymin=-1, ymax=1, alpha=0.3, fill = 'grey60') +
    annotate("rect",xmin=x[3,1], xmax=x[3,2], ymin=-1, ymax=1, alpha=0.3, fill = 'grey20') +
    annotate("rect",xmin=x[4,1], xmax=x[4,2], ymin=-1, ymax=1, alpha=0.3, fill = 'grey60') +
    geom_vline(xintercept = decision.time, size=1, color="red", linetype ='dashed') +
    theme_bw(base_size = 14) +
    geom_line(size=1, color="black") +
    geom_hline(yintercep = 0, size = 0.5, linetype='dashed') +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = -1) +
    ylab("RDV") +
    xlab("Time in ms") +
    scale_x_continuous(expand = c(0, 0),limits = c(1,decision.time)) +
    scale_y_continuous(limits = c(-1.2,1.2), breaks = c(-1,1)) + theme(panel.border = element_blank(),
                                                                       axis.title.x=element_text(face="bold"),
                                                                       axis.title.y=element_text(face="bold"))

  for (i in 1:nr.fix){
    if ((i %/% 2 - i / 2) == 0){
      m = 'left'
    } else {
      m = 'right'
    }

    if ((x[i,1] -  x[i,2])  <  -100){
      addm.example.plot = addm.example.plot + annotate('text', label=m, x = (x[i,1] + x[i,2]) / 2, y = -0.8)
    }
  }
  # ------------------------------------------------------------------------------------------
  return(addm.example.plot)
}
