# INSTALL PACKAGE IF NECESSARY ---------------------------------------------------------------
devtools::install_github('AlexanderFengler/addmtoolbox',build_vignette=TRUE)
library("addmtoolbox")
# --------------------------------------------------------------------------------------------

# Understanding the format need for input data -----------------------------------------------
# Choice Data
View(head(addm_data_choice))

# Eyetracking Data
View(head(addm_data_eye))
# --------------------------------------------------------------------------------------------

# Preprocess data for easy usage with addmtoolbox --------------------------------------------

# EITHER GENERATE FAKE DATA
addm_dat = addm_generate_artificial(set.size = 2,
                                    nr.reps = 50,
                                    nr.conditions = 10,
                                    timestep = 10,
                                    rtbinsize = 100,
                                    possible.valuations = c(0,1,2,3),
                                    model.parameters = c(0,0.002,0.015,0.5))

cur_dat = addm_dat$choice.dat
cur_dat$val.diff = cur_dat$v1-cur_dat$v2
cur_dat = cur_dat %>% group_by(val.diff) %>% summarise(rt = mean(rt))
hist(addm_dat$choice.dat$rt)
p = ggplot(data=cur_dat,aes(y=rt, x = val.diff)) + geom_line()

# OR READ IN SUPPLIED DATA FRAME
addm_dat  = addm_preprocess(choice.dat = addm_data_choice,
                            eye.dat = addm_data_eye,
                            rtbinsize = 100, timestep = 10)

# Look at the output
addm_dat

# View it separately
addm_dat$choice.dat

addm_dat$eye.dat

addm_dat$conditions.dat
# --------------------------------------------------------------------------------------------

# Run model by trial -------------------------------------------------------------------------
addm_loglik_trial = addm_fit_grid(addm_dat,
                                  fit.type = 'trial',
                                  drifts = c(0.001,0.0015,0.002,0.0025,0.003),
                                  sd = c(0.05,0.06,0.07,0.08,0.09),
                                  nr.reps = 2000)
# Plot Log likelihood ------------------------------------------------------------------------
addm_plot_loglik(addm_loglik_trial)
# --------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------

# Run model by condition ---------------------------------------------------------------------
addm_loglik_condition = addm_fit_grid(addm_dat,
                                      fit.type = 'condition',
                                      fixation.model = 'fixedpath',
                                      parallel = 1,
                                      nr.reps = 5000,
                                      drifts = rep(0.0015,100),
                                      sds = 0.07,
                                      theta = 1)
# Plot Log likelihood ------------------------------------------------------------------------
addm_plot_loglik(addm_loglik_condition)
# --------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------

# Run model dynamic --------------------------------------------------------------------------
addm_loglik_trial = addm_fit_grid(addm_dat,
                                  fit.type = 'dyn',
                                  drifts = c(0.001,0.0015,0.002,0.0025,0.003),
                                  sds = c(0.05,0.06,0.07,0.08,0.09,0.1))
# --------------------------------------------------------------------------------------------

# Simulate Data ------------------------------------------------------------------------------
addm_data_full_output = addm_run_by_condition(choice.dat = addm_dat$choice.dat,
                                              conditions.dat = addm_dat$conditions.dat,
                                              model.parameters = c(0,0.001,0.07,0.75),
                                              output.type = 'full',
                                              fixation.model = 'fixedpath',
                                              nr.reps = 1000)
# --------------------------------------------------------------------------------------------

# Plot some psychometrics --------------------------------------------------------------------
addm2_plot_family(choice.dat = addm_dat$choice.dat,
                  addm.output = addm_data_full_output)
# --------------------------------------------------------------------------------------------
