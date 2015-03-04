#' Preprocess experimental data for easy usage with addmtoolbox package
#' \code{addm_dataprep}
#' @param choice.dat data.frame or data.table storing the item valuations (v1,v2...) , reaction times in ms (rt), decisions (decision) and a trial id column (id).
#' @param eye.dat A data.frame or data.table storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and a trial id column (id). Can all be initialized as zero columns when a by condition fit is attempted.
#' @param timestep integer that provides the timestep-size that is used in the simulations (in ms).
#' @param rtbinsize integer that provides the binsize that is used for binning reaction time data (in ms).
#' @param fit.type string defined as 'condition' or 'trial' for the respective model fits
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Data Preprocessing for easy usage of addmtoolbox package
#' @return list of three components. data.table that stores all unique choice set conditions. data.table that stores eyetracking.data adjusted to be usable for by trial fits. data.table that stores by trial choice data. trial id and condition id variable link all data.tables
#' @export
addm_dataprep = function(choice.dat =  data.table(v1 = 0,v2 = 0,id = 0,rt = 0, decision = 0, id = 0),
                         eye.dat = data.table(fixloc = 0,fixnr = 1, fixdur= 0, id = 0),
                         timestep = 10,
                         rtbinsize = 100){

  # TURN INPUT DATA FROM ORIGNAL FORMAT TO KEYED DATA TABLES -------------------------------------
  choice = as.data.table(choice.dat)
  setkey(choice,id)
  eye = as.data.table(eye.dat)
  setkey(eye,id,fixnr)
  setkey(eye,id)
  # ----------------------------------------------------------------------------------------------

  # PREPARE DATA TABLES --------------------------------------------------------------------------

  # We get a data.table that provides unique trial conditions
  conditions = choice %>% select(-rt,-decision,-id) %>% distinct()
  conditions = as.data.table(conditions)

  # Assumed is that id variable is provided by trial
  # To allow 'by condition' fits we add a 'condition_id' variable
  conditions$condition_id = do.call(paste,c(conditions[,names(conditions),with=FALSE],sep='_'))
  choice$condition_id = do.call(paste,c(choice[,grep("^v[1-9]*",names(choice),value = TRUE),with=FALSE],sep='_'))

  if(length(eye[,fixloc]) > 1){
    temp = choice %>% select(id, condition_id)
    setkey(temp,id)
    eye = temp[eye]
  }

  # EXTRAS NEEDED FOR BY TRIAL FITS ----------------------------------------------------------------

  # ROUND FIXATION TIMES ACCODING TO TIMESTEP-SIZE - ADJUST RTS IN CHOICE DATA FRAME ACCORDINGLY ---
  eye$fixdur =  timestep * round(eye$fixdur/timestep)
  rts = eye %>% group_by(id) %>% summarize(rt = sum(fixdur))
  setkey(rts,id)
  choice = choice %>% select(-rt)
  choice = choice[rts]
  # ------------------------------------------------------------------------------------------------

  # DEFINE UPPER AND LOWER BOUNDS FOR RTBINS PER TRIAL ---------------------------------------------
  choice$rtup = ((choice$rt + rtbinsize) %/% rtbinsize) * rtbinsize
  choice$rtdown = choice$rtup - rtbinsize
  # ------------------------------------------------------------------------------------------------

  # ADJUST LAST FIXATION TIMES SO THAT WE HAVE DATA UNTIL UPPER BOUND OF RT BIN --------------------
  # We simply add rtbinsize to any last fixation, which ensures data to upper bound
  eye = eye %>% group_by(id) %>% mutate(max.fix = n())
  #eye$last.fix = 0
  eye[fixnr == max.fix,last.fix:= 1]
  eye[last.fix == 1,fixdur := fixdur + rtbinsize]
  eye = eye %>% select(-last.fix)
  # -------------------------------------------------------------------------------------------------

  # RETURN LIST WITH CHOICE, EYE, CONDITIONS DATA ---------------------------------------------------
  return(list(choice.dat = choice,
              eye.dat = eye,
              conditions.dat = conditions))
  # -------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------
}
