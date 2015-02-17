###### aDDM Graph real log likelihood
###### I am going to do a simple line plot 

cur.dat = core.aDDM.opti 
cur.dat = cur.dat[complete.cases(cur.dat),]
cur.dat$pcolor[cur.dat$Theta == 0.5] = "grey80"
#cur.dat$pcolor[cur.dat$Theta == 0.5] = "grey60"
cur.dat$pcolor[cur.dat$Theta == 0.6] = "grey40"
#cur.dat$pcolor[cur.dat$Theta == 0.6] = "grey20"
cur.dat$pcolor[cur.dat$Theta == 0.7] = "black"

set.sizes = c(4,6,8)

for (set.size in set.sizes){
  temp = cur.dat[cur.dat$Set.size == set.size,]
  
  #cur.max = max(temp$Log.Likelihood) 
  
  main.str= paste("aDDM // Theta color coded // Set Size",toString(set.size),sep=" ")
  scatterplot3d(temp$Drift.Rate,temp$Variance,temp$Log.Likelihood,color=temp$pcolor,pch=16,angle=70,main=main.str)  
}




temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.4 = temp[temp$Set.size == 4,]
temp.4 = temp.4[temp.4$Theta == 0.6,]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
#temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Variance == core.aDDM.params.4[5],]
p.1 = ggplot(temp.4,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.6 // Set Size 4")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.4 = temp[temp$Set.size == 4,]
temp.4 = temp.4[temp.4$Theta == 0.7,]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
#temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Variance == core.aDDM.params.4[5],]
p.2 = ggplot(temp.4,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.7 // Set Size 4")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.4 = temp[temp$Set.size == 4,]
temp.4 = temp.4[temp.4$Theta == 0.8,]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
#temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Variance == core.aDDM.params.4[5],]
p.3 = ggplot(temp.4,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.8 // Set Size 4")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.4 = temp[temp$Set.size == 4,]
temp.4 = temp.4[temp.4$Theta == 0.9,]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
#temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Variance == core.aDDM.params.4[5],]
p.4 = ggplot(temp.4,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.9 // Set Size 4")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.4 = temp[temp$Set.size == 4,]
temp.4 = temp.4[temp.4$Theta == 1,]
temp.4.max = temp.4[temp.4$Log.Likelihood == max(temp.4$Log.Likelihood),] 
#temp.4.generated = temp.4[temp.4$Drift.Rate == core.aDDM.params.4[3] & temp.4$Variance == core.aDDM.params.4[5],]
p.5 = ggplot(temp.4,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.4.max,aes(x=temp.4.max$Drift.Rate,y=temp.4.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.4.generated,aes(x=temp.4.generated$Drift.Rate,y=temp.4.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 1 // Set Size 4")


##########################################################################################

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.6 = temp[temp$Set.size == 6,]
temp.6 = temp.6[temp.6$Theta == 0.6,]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
#temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Variance == core.aDDM.params.6[5],]
p.6 = ggplot(temp.6,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.6 // Set Size 6")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.6 = temp[temp$Set.size == 6,]
temp.6 = temp.6[temp.6$Theta == 0.7,]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
#temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Variance == core.aDDM.params.6[5],]
p.7 = ggplot(temp.6,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.7 // Set Size 6")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.6 = temp[temp$Set.size == 6,]
temp.6 = temp.6[temp.6$Theta == 0.8,]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
#temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Variance == core.aDDM.params.6[5],]
p.8 = ggplot(temp.6,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.8 // Set Size 6")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.6 = temp[temp$Set.size == 6,]
temp.6 = temp.6[temp.6$Theta == 0.9,]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
#temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Variance == core.aDDM.params.6[5],]
p.9 = ggplot(temp.6,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.9 // Set Size 6")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.6 = temp[temp$Set.size == 6,]
temp.6 = temp.6[temp.6$Theta == 1,]
temp.6.max = temp.6[temp.6$Log.Likelihood == max(temp.6$Log.Likelihood),] 
#temp.6.generated = temp.6[temp.6$Drift.Rate == core.aDDM.params.6[3] & temp.6$Variance == core.aDDM.params.6[5],]
p.10 = ggplot(temp.6,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.6.max,aes(x=temp.6.max$Drift.Rate,y=temp.6.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.6.generated,aes(x=temp.6.generated$Drift.Rate,y=temp.6.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 1 // Set Size 6")

##########################################################################################

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.8 = temp[temp$Set.size == 8,]
temp.8 = temp.8[temp.8$Theta == 0.6,]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
#temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Variance == core.aDDM.params.8[5],]
p.11 = ggplot(temp.8,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.6// Set Size 8")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.8 = temp[temp$Set.size == 8,]
temp.8 = temp.8[temp.8$Theta == 0.7,]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
#temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Variance == core.aDDM.params.8[5],]
p.12 = ggplot(temp.8,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.7// Set Size 8")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.8 = temp[temp$Set.size == 8,]
temp.8 = temp.8[temp.8$Theta == 0.8,]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
#temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Variance == core.aDDM.params.8[5],]
p.13 = ggplot(temp.8,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.8// Set Size 8")

temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.8 = temp[temp$Set.size == 8,]
temp.8 = temp.8[temp.8$Theta == 0.9,]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
#temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Variance == core.aDDM.params.8[5],]
p.14 = ggplot(temp.8,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 0.9// Set Size 8")


temp = core.aDDM.opti.1[complete.cases(core.aDDM.opti.1),]
temp.8 = temp[temp$Set.size == 8,]
temp.8 = temp.8[temp.8$Theta == 1,]
temp.8.max = temp.8[temp.8$Log.Likelihood == max(temp.8$Log.Likelihood),] 
#temp.8.generated = temp.8[temp.8$Drift.Rate == core.aDDM.params.8[3] & temp.8$Variance == core.aDDM.params.8[5],]
p.15 = ggplot(temp.8,aes(x=Drift.Rate,y=Variance,fill=Log.Likelihood)) + 
  geom_raster() + 
  theme_bw() + 
  scale_fill_gradient(low="black",high="white") +
  geom_point(data=temp.8.max,aes(x=temp.8.max$Drift.Rate,y=temp.8.max$Variance),shape=8,color="red",size=5) +
  #geom_point(data=temp.8.generated,aes(x=temp.8.generated$Drift.Rate,y=temp.8.generated$Variance),shape=18,color="green",size=7) + 
  ggtitle("aDDM Subject 3// Theta Constant 1// Set Size 8")