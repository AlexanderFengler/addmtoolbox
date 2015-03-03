#' Visualize Likelihood fitting
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Visualize Likelihood fitting
#' @return list of two ggplot objects. plot.coarse.grid stores a plot corresponding to coarse grid fits. plot.fine.grid analogously stores a plot for fine grid fits.
#' \code{addm_plot_loglik}
#' @export
#' @param logliks 'data.frame' or 'data.table' storing parameter values and corresponding log likelihoods as a result of a model fit.

addm_plot_loglik = function(logliks = data.table(drift = 0, theta = 0, sd = 0, non.decision.time = 0, nr.reps = 0, loglik = 0)){

logliks = as.data.table(logliks)
logliks$theta = as.factor(logliks$theta)
setkey(logliks,loglik)

writeLines(paste(' \nOptimal Parameters... \n \n',
      'Drift Rate: ',
      toString(logliks[1,drift]),'\n',
      'Theta: ',toString(logliks[1,theta]),'\n',
      'SD: ', toString(logliks[1,sd]), '\n',
      'Non decision time: ', toString(logliks[1,non.decision.time]),sep=''))

optimal.parameters = logliks[1,]

p.coarse = ggplot(data = logliks[logliks$coarse == 1,], aes(y=loglik,x=sd,group = theta,color= theta))  +
  facet_wrap(~ drift) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle=-45,hjust=-0.05)) +
  ggtitle("Log Likelihood (separated by drift rate)") +
  geom_line() + geom_point(data=optimal.parameters, size = 5, shape = 4, color = 'red')

if (any(logliks[,coarse] == 0)){
p.fine = ggplot(data = logliks[logliks$coarse == 0,], aes(y=loglik,x=sd,group = theta,color=as.factor(theta)))  +
  facet_wrap(~ drift) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle=-45,hjust=-0.05)) +
  ggtitle("Log Likelihood (separated by drift rate)") +
  geom_line()

return(list(plot.coarse.grid = p.coarse,
            plot.fine.grid = p.fine))
}

return(list(plot.coarse.grid = p.coarse,
            plot.fine.grid = 'Fine grid not supplied because absent in model fit'))
}


