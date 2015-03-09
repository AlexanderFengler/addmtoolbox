y  = data.table(v1.b1 = 0,
                v1.b2 = 1,
                v1.b3 = 2,
                v2.b1 = 0,
                v2.b2 = 1,
                v2.b3 = 2,
                v3.b3 = 2)
#cur.str = paste("^v[1-9999]",".a[",toString(1),"]",sep='')
#grep("^v[0-9].a[1]*",names(y))
#grep(cur.str,names(y))

num.attributes = 0
for (i in seq(1:1000)){
  cur.str = paste("^v[0-9]",".a[",toString(i),"]",sep='')

  if (length(grep(cur.str,names(y))) < 1){
    break
  }
  num.attributes = num.attributes + 1
}

y$v4.b1

exists("y$v1.b1")


tl = list(a = 0, b = 0)


choice$condition_id = do.call(paste,c(choice[,grep("^v[1-9]*",names(choice),value = TRUE),with=FALSE],sep='_'))

monkeyfun = function(...){
  if (exists("monkey")){
    return('monkey found')
  } else {
    return('monkey not found')
  }
}
