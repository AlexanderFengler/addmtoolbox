
# two data.frames / data.tables will be provided
# choice.dat - "subject", "setsize","v1","v2"...."vn","decision","rt","id"
# eye,dat - "subject","set_size","fixloc","fixnr","fixdur","trialid"

addm_dataprep = function(choice.dat,eye.dat,timestep.ms,rtbinsize){

  # Parameters:
  timesteps = 10 #in ms // represents the timestep size with which the aDDM shall propagate
  rtbinsize = 100 # in ms // represent the rtbin.size that we bin the data into

  choice = as.data.table(choice.dat)
  eye = as.data.table(eye.dat)

  setkey(choice,id)
  setkey(eye,id)
  # ----------------------------------------------------------------------------------------------
  # Round all fixation times according chosen timestep
  eye$fixdur =  timesteps * round(eye$fixdur/timestep.ms)
  rts = eye %>% group_by(id) %>% summarize(rt = sum(fixdur))
  setkey(rts,id)
  trials = unique(eye$id)

  # Update the rt's in choice data.table according to sum of all rounded durations for consistency
  for(trial in trials){
    choice[J(trial),rtr:=rts[J(trial),rt]]
    print(trial)
  }

  # Get two additions to RT, rtup -- the upper limit of the RT-Bin the current trials belongs to, rtdown -- the lower limit of the RT-Bin the current trial belongs to
  choice$rtup = ((choice$rt + rtbinsize) %/% rtbinsize) * rtbinsize
  choice$rtdown = choice$rtup - rtbinsize

  # To make sure that the simulation has enough data to run to the end of the RT-Bin that a given trial is assigned to (real fixation-duration data stops at a random point in this RT-Bin) we need to adjust the length of last fixationas to accomodate simulation runs until at least the end of the assigned RT-Bins

  #Step 1 get difference between current RT and end of RT-Bin
  choice$rtdiff = choice$rtup - choice$rt
  eye$rt.add.last = 0
  for(trial in trials){
    eye[J(trial),rt.add.last:=choice[J(trial),rtdiff]]
    print(trial)
  }

  # Step 2
  # We adjust the last fixation length for each trial such that we have data until the end of the rt.bin that the particular trial is assigned to
  # Moreover we add a little more to allow the aDDM to go on after the real rt.bin is passed (to assign the simulation run as a failure)


  # ADD LAST FIX BINARY MANUALLY HERE
  #####################
  eye$Duration[eye$Last.fix.binary == 1] = eye$Duration[eye$Last.fix.binary == 1] + eye$rt.add.last[eye$Last.fix.binary == 1] + rtbinsize
  #####################

  # Finalize the two data.tables

  # Create List that contains all necessary information to run the aDDM optimization in one place
  # ----------------------------------------------------------------------------------------------
  addm.ready.frame = list(subjects = c(0),
                          set.sizes = c(4),
                          choice.dat = choice,
                          eye.dat = eye,
                          drifts = seq(0.00,0.03,0.005),     #default seq(0.00,0.03,0.005)
                          thetas = seq(0.0,1,0.1),           #default seq(0.0,1,0.1)
                          sds = seq(0,0.2,0.025),            #default seq(0,0.2,0.025)
                          non.decision.times = 0,
                          timesteps.ms = 10,
                          nr.reps = 5000,
                          model.type = 'nomem',
                          output.type = "Opti",
                          fixation.model="Normal",
                          allow.extension = 1,
                          allow.fine.grid = 1,
                          generate = 0)
}
