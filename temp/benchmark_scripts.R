# Benchmark scripts:


# GENERATE FRAMES

# 1 INCREASING TRIAL NUMBER ---------------------------------------------------------------------------------

# Dynamic aDDM ----------------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,7)
mean.time = rep(0,7)

cur.num.conditions = 10
num.reps = c(10,20,40,80,160,320,640)
cnt = 1
for(cur.num.reps in num.reps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial_dynamic(choice.dat = ch, eye.dat = ey),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame.1  = data.frame(model.type = 'dynamic',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Simulation by Trial ---------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,7)
mean.time = rep(0,7)

cur.num.conditions = 10
num.reps = c(10,20,40,80,160,320,640)
cnt = 1
for(cur.num.reps in num.reps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial(choice.dat = ch, eye.dat = ey),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame.2  = data.frame(model.type = 'sim_by_trial',
                        number.trials = num.trials,
                        number.conditions = cur.num.conditions,
                        mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Simulation by Conditions ----------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,7)
mean.time = rep(0,7)

cur.num.conditions = 10
num.reps = c(10,20,40,80,160,320,640)
cnt = 1
for(cur.num.reps in num.reps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat
  co = x$conditions.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_condition(choice.dat = ch, conditions.dat = co),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame.3  = data.frame(model.type = 'sim_by_condition',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Collect and save output -----------------------------------------------------------------------------------

out.frame = rbind(out.frame.1,
                  out.frame.2,
                  out.frame.3)

write.table(x = out.frame, file = 'temp/benchmark_vary_trial_constant_conditions.txt',row.names = FALSE,  sep=' ')

# -----------------------------------------------------------------------------------------------------------

# Draw plot -------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------

# 2 INCREASING NUMBER OF CONDITIONS HOLDING TRIALS CONSTANT -------------------------------------------------
# Dynamic aDDM ----------------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,5)
mean.time = rep(0,5)

cur.num.conditions = c(5,10,20,50,100)
num.reps = c(200,100,50,20,10)
cnt = 1
for(cur.num.reps in num.reps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions[cnt])
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial_dynamic(choice.dat = ch, eye.dat = ey),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions[cnt]
  cnt = cnt + 1
}

out.frame.1  = data.frame(model.type = 'dynamic',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Simulation by Trial ---------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,5)
mean.time = rep(0,5)


cur.num.conditions = c(5,10,20,50,100)
num.reps = c(200,100,50,20,10)
cnt = 1
for(cur.num.reps in num.reps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions[cnt])
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial(choice.dat = ch, eye.dat = ey),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions[cnt]
  cnt = cnt + 1
}

out.frame.2  = data.frame(model.type = 'sim_by_trial',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Simulation by Conditions ----------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,5)
mean.time = rep(0,5)


cur.num.conditions = c(5,10,20,50,100)
num.reps = c(200,100,50,20,10)
cnt = 1
for(cur.num.reps in num.reps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions[cnt])
  ch = x$choice.dat
  ey = x$eye.dat
  co = x$conditions.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_condition(choice.dat = ch, conditions.dat = co),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions[cnt]
  cnt = cnt + 1
}

out.frame.3  = data.frame(model.type = 'sim_by_condition',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Collect and save output -----------------------------------------------------------------------------------

out.frame = rbind(out.frame.1,
                  out.frame.2,
                  out.frame.3)

write.table(x = out.frame, file = 'temp/benchmark_vary_condition_constant_trials.txt',row.names = FALSE,  sep=' ')

# -----------------------------------------------------------------------------------------------------------

# Draw plot -------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------

# 3 CONSTANT  TRIAL NUMBER AND CONDITIONS AND VARYING NUMBER OF SIMULATIONS ---------------------------------
# Dynamic aDDM ----------------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,10)
mean.time = rep(0,10)
sim.nums = c(10,20,40,80,160,320,640,1280,2560,5120)

cur.num.conditions = 100
num.reps = 10
cnt = 1
for(cur.sim.num in sim.nums){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial_dynamic(choice.dat = ch, eye.dat = ey),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame.1  = data.frame(model.type = 'dynamic',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          number.simulations = sim.nums,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Simulation by Trial ---------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,10)
mean.time = rep(0,10)
sim.nums = c(10,20,40,80,160,320,640,1280,2560,5120)

cur.num.conditions = 100
num.reps = 10
cnt = 1
for(cur.sim.num in sim.nums){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial(choice.dat = ch, eye.dat = ey , nr.reps = cur.sim.num),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame.2  = data.frame(model.type = 'sim_by_trial',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          number.simulations = sim.nums,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Simulation by Conditions ----------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,10)
mean.time = rep(0,10)
sim.nums = c(10,20,40,80,160,320,640,1280,2560,5120)

cur.num.conditions = 100
num.reps = 10
cnt = 1
for(cur.sim.num in sim.nums){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat
  co = x$conditions.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_condition(choice.dat = ch, conditions.dat = co, nr.reps = cur.sim.num),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame.3  = data.frame(model.type = 'sim_by_condition',
                          number.trials = num.trials,
                          number.conditions = cur.num.conditions,
                          number.simulations = sim.nums,
                          mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

# Collect and save output -----------------------------------------------------------------------------------

out.frame = rbind(out.frame.1,
                  out.frame.2,
                  out.frame.3)

write.table(x = out.frame, file = 'temp/benchmark_vary_simnumber_constant_trials_conditions.txt',row.names = FALSE,  sep=' ')

# -----------------------------------------------------------------------------------------------------------

# Draw plot -------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------


# VARY stateSTEP IN DYNAMIC DDM -----------------------------------------------------------------------------
num.trials = rep(0,6)
mean.time = rep(0,6)

state.steps = c(0.2,0.1,0.02,0.01,0.002,0.001)
vertical.grid.steps = c(10,20,100,200,1000,2000)
cur.num.conditions = 10
num.reps = 10
cnt = 1
for(cur.state.steps in state.steps){
  x = addm_generate_artificial(nr.reps = cur.num.reps,nr.conditions = cur.num.conditions)
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial_dynamic(choice.dat = ch, eye.dat = ey, state.step = cur.state.steps),times=5)$time/1e+09)
  num.trials[cnt] = cur.num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame  = data.frame(model.type = 'dynamic',
                        number.trials = num.trials,
                        vertical.grid.steps = vertical.grid.steps,
                        number.conditions = cur.num.conditions,
                        mean.time = mean.time)

write.table(x = out.frame, file = 'temp/benchmark_vary_statesteps.txt',row.names = FALSE,  sep=' ')
# ------------------------------------------------------------------------------------------------------------

# Draw plot -------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------


# Vary set size for by trial simulations aDDM ---------------------------------------------------------------

# Simulation by Trial ---------------------------------------------------------------------------------------
# Generate Data
num.trials = rep(0,8)
mean.time = rep(0,8)
set.sizes = c(2,4,6,8,10,12,14,16)

cur.num.conditions = 10
num.reps = 10
cnt = 1
for(cur.set.size in set.sizes){
  x = addm_generate_artificial(nr.reps = num.reps,nr.conditions = cur.num.conditions, set.size = cur.set.size)
  ch = x$choice.dat
  ey = x$eye.dat

  mean.time[cnt] = mean(microbenchmark(addm_run_by_trial(choice.dat = ch, eye.dat = ey),times=5)$time/1e+09)
  num.trials[cnt] = num.reps*cur.num.conditions
  cnt = cnt + 1
}

out.frame  = data.frame(model.type = 'sim_by_trial',
                        number.trials = num.trials,
                        number.conditions = cur.num.conditions,
                        set.sizes = set.sizes,
                        mean.time = mean.time)
# -----------------------------------------------------------------------------------------------------------

write.table(x = out.frame, file = 'temp/benchmark_vary_statesteps.txt', row.names = FALSE, sep=' ')
# ------------------------------------------------------------------------------------------------------------

