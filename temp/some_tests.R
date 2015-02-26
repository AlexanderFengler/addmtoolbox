a = data.table(id = c(1,2,3),
               vara = c(1,1,1,1,1,1),
               vard = c(1,1,1,1,1,1))

b = data.table(id = c(1,2,3),
               varb = c(1,2,3))

setkey(a,id)
setkey(b,id)

a[b]

t = function(logfile = "curlogfile.txt",x = a){
write.table(x,logfile,sep= " ",quote=FALSE, col.names=TRUE,row.names=FALSE)
}

x = "myname"

x = data.frame(id = rnorm(3,0,1),
               id2 = rnorm(3,0,1))

#hi = data.frame(eval(parse(text = "x2")) = 4) # = data.frame(paste("x",toString(2),sep="") =  4)

y  = c("id3","id4","id5")

x[,y] = 0


t = function(x = mean(c(1,2,3))){
  return(x)
}
