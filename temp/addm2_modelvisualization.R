library(ggplot2)

rdv = 0

theta = 0.1
v1 = 2
v2 = 4
d = 0.002
noise.sd = 0.02

rdv.vec = seq(1:1000)
cnt = 1
for (fix.pos in rep(c(1,2,1,2),each=250)){

  if (fix.pos == 1){
    rdv[1] = rdv[1] + d*(v1 - theta*v2) + rnorm(1,0,noise.sd)
  }
  if (fix.pos == 2){
    rdv[1] = rdv[1] + d*(theta*v1 - v2) + rnorm(1,0,noise.sd)
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


addm.example.plot = ggplot(out.mat, aes(x=time,y=rdv)) +
  geom_vline(xintercept = decision.time, size=1, color="red", linetype='dashed') + theme_bw(base_size = 14) +
  geom_line(size=1, color="black") +
  geom_hline(yintercep = 0, size = 0.5, linetype='dashed') +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = -1) +
  ylab("RDV") +
  xlab("Time in ms") +
  scale_x_continuous(expand = c(0, 0),limits = c(1,decision.time)) +
  scale_y_continuous(limits = c(-1.2,1.2)) + theme(panel.border = element_blank(),
                                                   axis.title.x=element_text(face="bold"),
                                                   axis.title.y=element_text(face="bold"))
