# aDDM Graphing LogLikelihood TEMPLATE

#### Vary one parameter

#### CASE: DRIFT RATE
temp = core.aDDM.fake.opti.drift.1000
p = ggplot(temp,aes(x=Drift.Rate,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
    geom_line() + 
    theme_bw() +
    geom_vline(xintercept=core.aDDM.params.4[3],color="green") +
    scale_color_grey(start=0.7,end=0.1) +
    ggtitle("1000 Simulations")

temp = core.aDDM.fake.opti.drift.2000
p = ggplot(temp,aes(x=Drift.Rate,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
  geom_line() + 
  theme_bw() +
  geom_vline(xintercept=core.aDDM.params.4[3],color="green") +
  scale_color_grey(start=0.7,end=0.1) +
  ggtitle("2000 Simulations")


temp = core.aDDM.fake.opti.drift.3000
p = ggplot(temp,aes(x=Drift.Rate,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
  geom_line() + 
  theme_bw() +
  geom_vline(xintercept=core.aDDM.params.4[3],color="green") +
  scale_color_grey(start=0.7,end=0.1) +
  ggtitle("3000 Simulations")


#### CASE: THETA

temp = core.aDDM.fake.opti.theta.1000
p = ggplot(temp,aes(x=Theta,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
    geom_line() + 
    theme_bw() +
    geom_vline(xintercept=core.aDDM.params.4[4],color="green") +
    scale_color_grey(start=0.7,end=0.1) +
    ggtitle("1000 Simulations")

temp = core.aDDM.fake.opti.theta.2000
p = ggplot(temp,aes(x=Theta,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
    geom_line() + 
    theme_bw() +
    geom_vline(xintercept=core.aDDM.params.4[4],color="green") +
    scale_color_grey(start=0.7,end=0.1) +
    ggtitle("2000 Simulations")


temp = core.aDDM.fake.opti.theta.3000
p = ggplot(temp,aes(x=Theta,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
    geom_line() + 
    theme_bw() +
    geom_vline(xintercept=core.aDDM.params.4[4],color="green") +
    scale_color_grey(start=0.7,end=0.1) +
    ggtitle("3000 Simulations")


#### CASE: Variance

temp = core.aDDM.fake.opti.var.1000
p = ggplot(temp,aes(x=Variance,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
  geom_line() + 
  theme_bw() +
  geom_vline(xintercept=core.aDDM.params.4[5],color="green") +
  scale_color_grey(start=0.7,end=0.1) +
  ggtitle("1000 Simulations")

temp = core.aDDM.fake.opti.var.2000
p = ggplot(temp,aes(x=Variance,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
  geom_line() + 
  theme_bw() +
  geom_vline(xintercept=core.aDDM.params.4[5],color="green") +
  scale_color_grey(start=0.7,end=0.1) +
  ggtitle("2000 Simulations")


temp = core.aDDM.fake.opti.var.3000
p = ggplot(temp,aes(x=Variance,y=Log.Likelihood, group = Set.size, color = as.factor(Set.size))) + 
  geom_line() + 
  theme_bw() +
  geom_vline(xintercept=core.aDDM.params.4[5],color="green") +
  scale_color_grey(start=0.7,end=0.1) +
  ggtitle("3000 Simulations")


#### Vary two parameters

#### SET_SIZE 4

#### CASE: Drift and Variance 4
temp = core.aDDM.fake.opti.drift.var.1000[complete.cases(core.aDDM.fake.opti.drift.var.1000),]
temp.4 = temp[temp$Set.size == 4,c(3,4,7)]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Variance == core.aDDM.params.4[5],]
p = ggplot(temp.4,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
    geom_raster() + 
    theme_bw() + 
    scale_fill_gradient(low="black",high="white") +
    geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Variance),shape=8,color="red",size=5) +
    geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Variance),shape=18,color="green",size=7) + 
    ggtitle("aDDM // Theta Constant // Set Size 4")

#### CASE: Drift and Theta 4
temp = core.aDDM.fake.opti.drift.theta.1000[complete.cases(core.aDDM.fake.opti.drift.theta.1000),]
temp.4 = temp[temp$Set.size == 4,c(3,5,7)]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Theta == core.aDDM.params.4[4],]
p = ggplot(temp.4,aes(x=Drift.Rate,y=Theta,fill=Log.Likelihood)) + 
    geom_raster() + 
    theme_bw() + 
    scale_fill_gradient(low="black",high="white") +
    geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Theta),shape=8,color="red",size=5) +
    geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Theta),shape=18,color="green",size=7) +
    ggtitle("aDDM // Variance Constant // Set Size 4")

#### CASE: Variance and Theta 4
temp = core.aDDM.fake.opti.var.theta.1000[complete.cases(core.aDDM.fake.opti.var.theta.1000),]
temp.4 = temp[temp$Set.size == 4,c(4,5,7)]
temp.4 = temp.4[complete.cases(temp.4),]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
temp.4.generated = temp.4[temp.4$Theta == core.aDDM.params.4[4] & temp.4$Variance == core.aDDM.params.4[5],]
p = ggplot(temp.4,aes(x=Variance,y=Theta,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.4.max,aes(x=temp.4.max$Variance,y=temp.4.max$Theta),shape=8,color="red",size=5) +
  geom_point(data=temp.4.generated,aes(x=temp.4.generated$Variance,y=temp.4.generated$Theta),shape=18,color="green",size=7) + 
  ggtitle("aDDM // Drift-Rate Constant // Set Size 4")

#### SET_SIZE 6

#### CASE: Drift and Variance 6
temp = core.aDDM.fake.opti.drift.var.1000[complete.cases(core.aDDM.fake.opti.drift.var.1000),]
temp.6 = temp[temp$Set.size == 6,c(3,4,7)]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Variance == core.aDDM.params.6[5],]
p = ggplot(temp.6,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Variance),shape=8,color="red",size=5) +
  geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM // Theta Constant // Set Size 6")

#### CASE: Drift and Theta 6
temp = core.aDDM.fake.opti.drift.theta.1000[complete.cases(core.aDDM.fake.opti.drift.theta.1000),]
temp.6 = temp[temp$Set.size == 6,c(3,5,7)]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Theta == core.aDDM.params.6[4],]
p = ggplot(temp.6,aes(x=Drift.Rate,y=Theta,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Theta),shape=8,color="red",size=5) +
  geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Theta),shape=18,color="green",size=7) + 
  ggtitle("aDDM // Variance Constant // Set Size 6")

#### CASE: Variance and Theta 6
temp = core.aDDM.fake.opti.var.theta.1000[complete.cases(core.aDDM.fake.opti.var.theta.1000),]
temp.6 = temp[temp$Set.size == 6,c(4,5,7)]
temp.6 = temp.6[complete.cases(temp.6),]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
temp.6.generated = temp.6[temp.6$Theta == core.aDDM.params.6[4] & temp.6$Variance == core.aDDM.params.6[5],]
p = ggplot(temp.6,aes(x=Variance,y=Theta,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Variance,y=temp.6.max$Theta),shape=8,color="red",size=5) +
  geom_point(data=temp.6.generated,aes(x=temp.6.generated$Variance,y=temp.6.generated$Theta),shape=18,color="green",size=7) +
  ggtitle("aDDM // Drift-Rate Constant // Set Size 6")


#### SET_SIZE 8

#### CASE: Drift and Variance 8
temp = core.aDDM.fake.opti.drift.var.1000[complete.cases(core.aDDM.fake.opti.drift.var.1000),]
temp.8 = temp[temp$Set.size == 8,c(3,4,7)]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Variance == core.aDDM.params.8[5],]
p = ggplot(temp.8,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Variance),shape=8,color="red",size=5) +
  geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM // Theta Constant // Set Size 8")

#### CASE: Drift and Theta 8
temp = core.aDDM.fake.opti.drift.theta.1000[complete.cases(core.aDDM.fake.opti.drift.theta.1000),]
temp.8 = temp[temp$Set.size == 8,c(3,5,7)]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Theta == core.aDDM.params.8[4],]
p = ggplot(temp.8,aes(x=Drift.Rate,y=Theta,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Theta),shape=8,color="red",size=5) +
  geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Theta),shape=18,color="green",size=7) +
  ggtitle("aDDM // Variance Constant // Set Size 8")

#### CASE: Variance and Theta 8
temp = core.aDDM.fake.opti.var.theta.1000[complete.cases(core.aDDM.fake.opti.var.theta.1000),]
temp.8 = temp[temp$Set.size == 8,c(4,5,7)]
temp.8 = temp.8[complete.cases(temp.8),]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
temp.8.generated = temp.8[temp.8$Theta == core.aDDM.params.8[4] & temp.8$Variance == core.aDDM.params.8[5],]
p = ggplot(temp.8,aes(x=Variance,y=Theta,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Variance,y=temp.8.max$Theta),shape=8,color="red",size=5) +
  geom_point(data=temp.8.generated,aes(x=temp.8.generated$Variance,y=temp.8.generated$Theta),shape=18,color="green",size=7) + 
  ggtitle("aDDM // Drift-Rate Constant // Set Size 8")


##### VARY ALL PARAMETERS

cur.dat = core.aDDM.fake.opti 
cur.dat = cur.dat[complete.cases(cur.dat),]
set.sizes = c(4,6,8)

for (set.size in set.sizes){
  temp = cur.dat[cur.dat$Set.size == set.size,]
  
  cur.max = max(temp$Log.Likelihood) 
  temp$pcolor[temp$Log.Likelihood == cur.max] = "red"
  temp$pcolor[temp$Log.Likelihood != cur.max] = "black"
  
  main.str= paste("Max.Log.Lik. // aDDM // Set Size",toString(set.size),sep=" ")
  scatterplot3d(temp$Drift.Rate,temp$Variance,temp$Theta,color=temp$pcolor,pch=16,angle=70,main=main.str)  
}

