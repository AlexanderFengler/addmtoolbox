#' Shapes up a 'data.frame' or 'data.table' to be straight usable with the other addm_*** functions
#' \code{addm_dataprep} List of three components. A 'data.table' that stores all unique choice set conditions. A 'data.table' that stores eyetracking.data adjusted to be usable for by trial fits. A 'data.table' that stores by trial choice data. An id variable links all three data.tables.
#' @param choice.dat A 'data.frame' or  'data.table' storing the item valuations (v1,v2...) , reaction times in ms (rt), decisions as decision and an id column (id). A by trial form is assumed.
#' @param eye.dat A 'data.frame' or 'data.table' storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and an id column (id). Can all be initialized as zero columns when a by condition fit is attempted.
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param rtbinsize An integer number that provides the binsize that is used for binning reaction time data (in ms).
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Data Preparation for easy supply to fit functions
#' @return returns a list of three components. A 'data.table' that stores all unique choice set conditions. A 'data.table' that stores eyetracking.data adjusted to be usable for by trial fits. A 'data.table' that stores by trial choice data. An id variable links all three data.tables.
#' @export
#' @useDynLib addmtools
#' @importFrom Rcpp sourceCpp

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
