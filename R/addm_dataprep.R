#' Shapes up a 'data.frame' or 'data.table' to be straight usable with the other addm_*** functions
#' \code{addm_dataprep} List of three components. A 'data.table' that stores all unique choice set conditions. A 'data.table' that stores eyetracking.data adjusted to be usable for by trial fits. A 'data.table' that stores by trial choice data. An id variable links all three data.tables.
#' @param choice.dat A 'data.frame' or  'data.table' storing the item valuations (v1,v2...) , reaction times in ms (rt), decisions as decision and an id column (id). A by trial form is assumed.
#' @param eye.dat A 'data.frame' or 'data.table' storing eyetracking data by trial. Fixation location (fixloc), Fixation number (fixnr), Fixation duration (fixdur) and an id column (id). Can all be initialized as zero columns when a by condition fit is attempted.
#' @param timestep An integer number that provides the timestep-size that is used in the simulations (in ms).
#' @param rtbinsize An integer number that provides the binsize that is used for binning reaction time data (in ms).
#' @param fit.type String defined as 'condition' or 'trial' for the respective model fits
#' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
#' @title Data Preparation for easy supply to fit functions
#' @return returns a list of three components. A 'data.table' that stores all unique choice set conditions. A 'data.table' that stores eyetracking.data adjusted to be usable for by trial fits. A 'data.table' that stores by trial choice data. An id variable links all three data.tables.
#' @export
#' @useDynLib addmtoolbox
#' @importFrom Rcpp sourceCpp

addm_dataprep = function(choice.dat =  data.table(v1 = 0,v2 = 0,id = 0,rt = 0, decision = 0, id = 0),
                         eye.dat = data.table(fixloc = 0,fixnr = 1, fixdur= 0, id = 0),
                         timestep = 10,
                         rtbinsize = 100,
                         fit.type = 'condition'){

  # TURN INPUT DATA FROM ORIGNAL FORMAT TO KEYED DATA TABLES -------------------------------------
  choice = as.data.table(choice.dat)
  setkey(choice,id)

  # Eyetracking data only needed for fit.type 'trial'
  if (fit.type == 'trial'){
    eye = as.data.table(eye.dat)
    setkey(eye,id,fixnr)
    setkey(eye,id)
  }
  # ----------------------------------------------------------------------------------------------

  # IN CASE FIT.TYPE IS 'condition' WE NEED CHOICES AND CONDITIONS TABLES ONLY -------------------


  # Assumed is that id variable is provided by trial
  # This is adjusted towards a by condition id variable here
  if (fit.type == 'condition'){
    conditions = choice %>% select(-rt,-decision,-id) %>% distinct()
    conditions = as.data.table(conditions)

    conditions$id = do.call(paste,c(conditions[,names(conditions),with=FALSE],sep='_'))
    choice$id = do.call(paste,c(choice[,grep("^v[1-9]*",names(choice),value = TRUE),with=FALSE],sep='_'))

    return(list(conditions.dat = conditions,
           choice.dat = choice))
  }
  # ----------------------------------------------------------------------------------------------

  # IN CASE FIT.TYPE IS 'trial' WE NEED TO DO MORE -----------------------------------------------

  # ROUND FIXATION TIMES ACCODING TO TIMESTEP-SIZE - ADJUST RTS IN CHOICE DATA FRAME ACCORDINGLY -
  if (fit.type == 'trial'){
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
    # We simply add rtbinsize to any last fixation, which ensures data to upper bound
    eye = eye %>% group_by(id) %>% mutate(max.fix = n())
    #eye$last.fix = 0
    eye[fixnr == max.fix,last.fix:= 1]
    eye[last.fix == 1,fixdur := fixdur + rtbinsize]
    eye = eye %>% select(-max.fix,-last.fix)
    # ----------------------------------------------------------------------------------------------

    # RETURN LIST WITH PREPARED CHOICE AND EYE DATA ------------------------------------------------
    return(list(eye.dat = eye,
                choice.dat = choice))
    # ----------------------------------------------------------------------------------------------
  }
  # ----------------------------------------------------------------------------------------------
}
