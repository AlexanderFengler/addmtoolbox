# Author: Alexander Fengler
# Date: Jan 15th 2015

# Script to run full output aDDM

# Author: Alexander Fengler
# Project: aDDM_4_6_8
# Date: Jan 14th 2015


# This script generates all necessary data needed as input to generate artificial aDDM data
# We can then use this to test the model fitting procedures which currently reside in file: run_aDDM_optimization

# RUN THIS PART ONLY IF CORRESPONDING FILE IS LOST FROM WORKSPACE -- OTHERWISE IT IS REUSABLE

#------------------------------------------------------------------------------------------------------------------------------------------

# Generating values for item 1 and 2 so that trialids are consistent later (we do not resample thos for each parameter combination)
#------------------------------------------------------------------------------------------------------------------------------------------
v1 = sample(c(0,1,2),300,replace=TRUE)
v2 = sample(c(0,1,2),300,replace=TRUE)
#------------------------------------------------------------------------------------------------------------------------------------------


# Generating a list that collect all necessary input data to generate artificial aDDM from the aDDM.R scripts (aDDM() function)
#------------------------------------------------------------------------------------------------------------------------------------------
test.dat = list(items_2 = data.table(v1=v1,
                                     v2=v2,
                                     trialid=1:300),
                                 decisions = data.table(decision = 0,
                                                        trialid = 1:300),
                                 trialids = 1:300,
                                 rts = data.table(rt = 0,
                                                  trialid = 1:300),
                                 fixations = data.table(fixloc = 0,trialid = 1:250),
                                 timestep.ms = 1,
                                 nr.reps = 1000,
                                 fixation.model = "FakePath",
                                 output.type = "Full",
                                 core.parameters = c(0.0005,0.1,0.014,0))  # order of core parameter: drift.rate,sd,theta,non.decision.time
#------------------------------------------------------------------------------------------------------------------------------------------

# GENERATING ARTIFICIAL OUTPUT
#------------------------------------------------------------------------------------------------------------------------------------------

temp.table = aDDM(test.dat$items_2,
                  test.dat$decisions,
                  test.dat$fixations,
                  test.dat$rts,
                  test.dat$trialids,
                  test.dat$core.parameters,
                  test.dat$nr.reps,
                  test.dat$output.type,
                  test.dat$fixation.model,
                  test.dat$timestep.ms,
                  1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test

temp.table$v1 = rep(v1,each=1000)
temp.table$v2 = rep(v2,each=1000)

temp.table$v.other = 0

temp.table$vother[temp.table$Decision == 1] = temp.table$v1[temp.table$Decision == 1]
temp.table$vother[temp.table$Decision == 2] = temp.table$v2[temp.table$Decision == 2]
temp.table$Value.not.last.attended[temp.table$Item.last.attended == 1] = temp.table$v2[temp.table$Item.last.attended == 1]
temp.table$Value.not.last.attended[temp.table$Item.last.attended == 2] = temp.table$v1[temp.table$Item.last.attended == 2]
temp.table$Drift.to.Noise.ratio = temp.table$SD / temp.table$Drift.Rate

temp.table$Last.vs.not.last = temp.table$Value.last.attended - temp.table$Value.not.last.attended

temp.table$Vdiff = temp.table$v1 = temp.table$v2


# antonio.artificial.choice.table.6 = aDDM(core.clean.model.test.dat$items_6,
#                                             core.clean.model.test.dat$decisions,
#                                             core.clean.model.test.dat$fixations,
#                                             core.clean.model.test.dat$rts,
#                                             core.clean.model.test.dat$trialids,
#                                             core.clean.model.test.dat$core.parameters,
#                                             core.clean.model.test.dat$nr.reps,
#                                             core.clean.model.test.dat$output.type,
#                                             core.clean.model.test.dat$fixation.model,
#                                             core.clean.model.test.dat$timestep.ms,
#                                             1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test
# 
# 
# antonio.artificial.choice.table.8 = aDDM(core.clean.model.test.dat$items_8,
#                                             core.clean.model.test.dat$decisions,
#                                             core.clean.model.test.dat$fixations,
#                                             core.clean.model.test.dat$rts,
#                                             core.clean.model.test.dat$trialids,
#                                             core.clean.model.test.dat$core.parameters,
#                                             core.clean.model.test.dat$nr.reps,
#                                             core.clean.model.test.dat$output.type,
#                                             core.clean.model.test.dat$fixation.model,
#                                             core.clean.model.test.dat$timestep.ms,
#                                             1) # the last model parameter tells the model to generate a data.frame instead of running log.likelihood test

#------------------------------------------------------------------------------------------------------------------------------------------

# Graphing it
#------------------------------------------------------------------------------------------------------------------------------------------

antonio.table$Last.vs.not.last = round(antonio.table$Last.vs.not.last,0)

table.reduced = antonio.table %>% group_by(Theta,Drift.Rate,Drift.to.Noise.ratio,Last.vs.not.last,Chosen.last.attended) %>% summarise(count = n()) 

table.reduced = table.reduced %>% group_by(Theta,Drift.Rate,Drift.to.Noise.ratio,Last.vs.not.last) %>% mutate(total = sum(count))

table.reduced$choice.p = table.reduced$count / table.reduced$total

table.reduced.2 =  antonio.table %>% group_by(Theta,Drift.Rate,Drift.to.Noise.ratio,Vdiff,Chosen.last.attended) %>% summarise(count = n()) 

table.reduced.2 = table.reduced.2 %>% group_by(Theta,Drift.Rate,Drift.to.Noise.ratio,Vdiff) %>% mutate(total = sum(count))

table.reduced.2$choice.p = table.reduced.2$count / table.reduced.2$total


p = ggplot(aes(x=Last.vs.not.last,y=choice.p,color=as.factor(Chosen.last.attended),group=as.factor(Chosen.last.attended)),data=table.reduced) + geom_line() + facet_grid(Drift.to.Noise.ratio~Theta)

p1 = ggplot(aes(x=Vdiff,y=choice.p,color=as.factor(Chosen.last.attended),group=as.factor(Chosen.last.attended)),data=table.reduced.2) + geom_line() + facet_grid(Drift.to.Noise.ratio~Theta)









#------------------------------------------------------------------------------------------------------------------------------------------




