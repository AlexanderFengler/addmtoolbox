# INSTALL PACKAGE IF NECESSARY ---------------------------------------------------------------
devtools::install_github('AlexanderFengler/addmtoolbox',build_vignette=TRUE)
library("addmtoolbox")
# --------------------------------------------------------------------------------------------

# HERE CHOICE SETS NOT INCLUDED IN PACKAGE FOR THE MOMENT ------------------------------------
addm_data_choice_ma = read.table('temp/data/addm_data_choice_ma.txt', header = TRUE)
addm_data_eye_ma = read.table('temp/data/addm_data_eye_ma.txt', header = TRUE)
# --------------------------------------------------------------------------------------------

# UNDERSTANDING INPUT DATA FORMAT ------------------------------------------------------------
# Choice Data
View(head(addm_data_choice_ma))
# Eyetracking Data
View(head(addm_data_eye_ma))
# --------------------------------------------------------------------------------------------

# PREPROCESS ---------------------------------------------------------------------------------

addm_dat  = addm_preprocess(choice.dat = addm_data_choice_ma,
                            eye.dat = addm_data_eye_ma,
                            rtbinsize = 100, timestep = 10)

# Look at the output
addm_dat

# View it separately
addm_dat$choice.dat

addm_dat$eye.dat

addm_dat$conditions.dat

addm_dat$attributes
# ----------------------------------------------------------------------------------------------

# Run model by trial -------------------------------------------------------------------------
addm_loglik_trial = addm_fit_grid(addm_dat,
                                  drifts = seq(0.001, 0.002, 0.0005),
                                  thetas = c(0, 0.5, 1),
                                  gammas = c(0, 0.5, 1),
                                  sds = 0.05,
                                  fit.type = 'trial')
# --------------------------------------------------------------------------------------------

# Run model by trial -------------------------------------------------------------------------
addm_loglik_trial = addm_fit_grid(addm_dat,
                                  drifts = seq(0.001, 0.002, 0.0005),
                                  thetas = c(0, 0.5, 1),
                                  gammas = c(0, 0.5, 1),
                                  sds = 0.05,
                                  fit.type = 'condition')
# --------------------------------------------------------------------------------------------

# Run model dynamic --------------------------------------------------------------------------
addm_loglik_trial = addm_fit_grid(addm_dat,
                                  drifts = seq(0.001, 0.002, 0.0005),
                                  thetas = c(0, 0.5, 1),
                                  gammas = c(0, 0.5, 1),
                                  sds = 0.05,
                                  fit.type = 'dyn')
# --------------------------------------------------------------------------------------------

# Get full output ----------------------------------------------------------------------------
addm_data_full_output = addm_run_by_condition(choice.dat = addm_dat$choice.dat,
                                              conditions.dat = addm_dat$conditions.dat,
                                              model.parameters = c(0.0015,0.5,0.5,0.07,0),
                                              nr.attributes = 2,
                                              output.type = 'full',
                                              nr.reps = 100)
# --------------------------------------------------------------------------------------------
