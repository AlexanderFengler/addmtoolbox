# Author: Alexander Fengler
# Date: February 23rd 2015
# Title: Prepare data to be supplied to
# choice.dat - "v1","v2"...."vn","decision","rt","id"
# eye,dat - "fixloc","fixnr","fixdur","trialid"

addm_dataprep = function(choice.dat =  data.table(v1 = c(1,2,3),v2 = c(3,2,1),id = c(1,2,3),rt = c(0,0,0), decision = c(0,0,0)),
                         eye.dat = data.table(fixloc = 0,fixnr = 1, fixdur= 0, id = c(1,2,3)),
                         timestep = 10,
                         rtbinsize = 100){

  # TURN INPUT DATA FROM ORIGNAL FORMAT TO KEYED DATA TABLES -------------------------------------
  choice = as.data.table(choice.dat)
  eye = as.data.table(eye.dat)
  setkey(choice,id)
  setkey(eye,id)

  conditions = choice %>% select(-rt,-decision) %>% distinct()
  # ----------------------------------------------------------------------------------------------

  # ROUND FIXATION TIMES ACCODING TO TIMESTEP-SIZE - ADJUST RTS IN CHOICE DATA FRAME ACCORDINGLY -
  eye$fixdur =  timestep * round(eye$fixdur/timestep)
  rts = eye %>% group_by(id) %>% summarize(rt = sum(fixdur))
  setkey(rts,id)
  choice = choice %>% select(-rt)
  choice = choice[rts]
  # ----------------------------------------------------------------------------------------------

  # DEFINE UPPER AND LOWER BOUNDS FOR RTBINS PER TRIAL -------------------------------------------
  choice$rtup = ((choice$rt + rtbinsize) %/% rtbinsize) * rtbinsize
  choice$rtdown = choice$rtup - rtbinsize
  # ----------------------------------------------------------------------------------------------

  # ADJUST LAST FIXATION TIMES SO THAT WE HAVE DATA UNTIL UPPER BOUND OF RT BIN ------------------
  # Step 1: get difference between current RT and end of RT-Bin
  choice$rt.add.last = choice$rtup - choice$rt
  temp  = choice %>% select(id,rt.add.last)
  setkey(temp,id)
  eye = eye[temp]

  # Step 2: Adjust length of last fixation to allow for length of rtbinsize
  eye = eye %>% group_by(id) %>% mutate(max.fix = n())
  eye$last.fix = 0
  eye[fixnr == max.fix,last.fix:= 1]
  eye[last.fix == 1,fixdur := last.fix + rt.add.last + rtbinsize]
  # ----------------------------------------------------------------------------------------------

  # RETURN LIST WITH PREPARED CHOICE AND EYE DATA ------------------------------------------------
  return(list(conditions.dat = conditions,
              eye.dat = eye,
              choice.dat = choice))
  # ----------------------------------------------------------------------------------------------
}
