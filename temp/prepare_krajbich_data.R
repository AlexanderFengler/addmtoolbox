# Title: Make Krajbich dataset usable for addmtoolbox

# Load necessary packages -------------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(foreign)
# -------------------------------------------------------------------------------------------------------

# Read in Data / rename / extract important information -------------------------------------------------
dat = as.data.table(read.dta("fixations_final.dta"))
dat$id = paste(dat$subject,dat$trial,sep='_')
dat$choice = dat$choice + 1

dat = dat %>% select(-computed_rt,-temp,-corr_fix_duration,-rev_fix_num,-num_fixations,-subject,-trial)
dat = dat %>% group_by(id) %>% mutate(rt = sum(event_duration))

setnames(dat,c("fixnr","fixdur","v1","v2","rt","decision","fixloc","id"))

eye.dat = dat[,list(fixloc,fixdur,fixnr,id)]

choice.dat = dat[,list(rt,decision,id)] %>% group_by(id) %>% summarize(decision = unique(decision),
                                                                       rt = unique(rt))
# -------------------------------------------------------------------------------------------------------

# Write to file -----------------------------------------------------------------------------------------
write.table(eye.dat,file = 'addm_data_eye.txt',row.names = FALSE, sep = ' ')
write.table(choice.dat, file='addm_data_choice.txt',row.names = FALSE,sep = ' ')
# -------------------------------------------------------------------------------------------------------


