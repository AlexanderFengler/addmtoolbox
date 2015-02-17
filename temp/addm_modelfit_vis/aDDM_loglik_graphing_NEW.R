test = data.table(read.table("Sim_logs/loglik_nomem_all_setsize_4.txt",header = TRUE))

setkey(test,Drift.Rate)

test$NR = seq(1,length(test$Drift.Rate),1)

p.temp = ggplot(data = test[test$Coarse == 1,], aes(y=Log.Likelihood,x=NR,group = as.factor(Theta),color=as.factor(Theta)))  + theme_bw() + geom_line()

p.c = ggplot(data = test[test$Coarse == 1,], aes(y=Log.Likelihood,x=SD,group = as.factor(Theta),color=as.factor(Theta)))  +  facet_wrap(~ Drift.Rate) + theme_bw() + geom_line()

p.f = ggplot(data = test[test$Coarse == 0,], aes(y=Log.Likelihood,x=SD,group = as.factor(Theta),color=as.factor(Theta)))  + facet_wrap(~ Drift.Rate) + theme_bw() + geom_line()

t = as.data.table(test)
setkey(t,Log.Likelihood)
which.max(t$Log.Likelihood)

t[1980,]

test = data.table(read.table("Sim_logs/aDDM_loglik_subject_4_setsize_4.txt",header = TRUE))

p.coarse = ggplot(data = test[test$Coarse == 1,],aes(y=SD,x=Drift.Rate,fill=Log.Likelihood,color=as.factor(Theta))) + geom_raster() + facet_wrap(~ Theta)

p.fine = ggplot(data = test[test$Coarse == 0,],aes(y=SD,x=Drift.Rate,fill=Log.Likelihood,color=as.factor(Theta))) + geom_raster() + facet_wrap(~ Theta)



