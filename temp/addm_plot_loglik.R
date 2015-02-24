addm_plot_loglik = function(logliks){

logliks = as.data.table(logliks)
setkey(logliks,loglik)

print(logliks[length(logliks[,drift]),])

p.coarse = ggplot(data = logliks[test$coarse == 1,], aes(y=loglik,x=sd,group = as.factor(theta),color=as.factor(theta)))  +
  facet_wrap(~ drift) + theme_bw() +
  geom_line()

p.fine = ggplot(data = logliks[test$coarse == 0,], aes(y=loglik,x=sd,group = as.factor(theta),color=as.factor(theta)))  +
  facet_wrap(~ drift) +
  theme_bw() + geom_line()

return(list(plot.coarse.grid = p.coarse,
            plot.fine.grid = p.fine))
}


