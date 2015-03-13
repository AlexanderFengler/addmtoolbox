test = function(x){
x = x + 1
return(x)
}


x  = matrix(c(1,1,2,2),ncol = 2)
x  = cbind(a = c(seq(1:1000)), b = rep(1,1000))
y  = t(apply(x,1,function(x){x[2] = 2 - 1.01^x[1]; return(x)}))
y = cbind(a = y[,2], b = (1:1000))

x = expand.grid(1:3,1:3)

y  = expand.grid(x,1:3)

data.matrix


tbfun = function(x,z){
  dx = dim(x)[2]
    for (i in 1:dx){
      y = paste('boundary_v',toString(i),'=', 'x[,',toString(i),']',sep='')
      eval(parse(text = y ))
      ztext = paste('expand.grid(z,boundary_v',toString(i),')',sep='')
      z = eval(parse(text = ztext))
    }
    return(z)
}


tbfun2 = function(inmat, addmat){
  add_row = prod(dim(addmat))
  add_col = dim(addmat)[2]
  n = add_row^add_col




  for (i in 1:add_col){
    tempmat = cbind(inmat,rep(addmat[,i],n)
    for (j in 1:add_row){
      tempmat =
    }


  }
}


do.call("rbind", rep(list(x), n))

y = matrix(c(1,2,3,4),ncol=2)
out = as.matrix(do.call('expand.grid', c(as.list(as.data.frame((x))), as.list(as.data.frame((y))))))
out = as.matrix(do.call('expand.grid', c(as.data.frame(x), as.data.frame(y))))



cbind(x,c(1,2,3))




