
x = addm_generate_by_condition_artificial(nr.reps = 10)
ch = x$choice.dat
ey = x$eye.dat


ch = as.data.table(read.table(file = 'temp/fakechoice.txt', header=TRUE))
ey = as.data.table(read.table(file = 'temp/fakeeye.txt',header=TRUE))

setkey(ch,id)
setkey(ey,id)

c = ch[1:200,]
e = ey[id %in% 1:200,]

#View(c)
#View(e)





addm_run_by_trial_dynamic(ch,ey)
