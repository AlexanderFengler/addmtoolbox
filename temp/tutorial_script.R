# Understanding the format need for input data -----------------------------------------------

# Choice Data
View(head(addm_data_choice))

# Eyetracking Data
View(head(addm_data_eye))

# --------------------------------------------------------------------------------------------

# Preprocess data for easy usage with addmtoolbox --------------------------------------------
addm_dat  = addm_dataprep(choice.dat = addm_data_choice,
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
                                       fit.type = 'trial')
# --------------------------------------------------------------------------------------------

# Plot Log likelihood ------------------------------------------------------------------------
addm_plot_loglik(addm_loglik_trial)
# --------------------------------------------------------------------------------------------

# Run model by condition ---------------------------------------------------------------------
addm_loglik_condition = addm_fit_grid(addm_dat,
                                           fit.type = 'condition')
# --------------------------------------------------------------------------------------------

# Plot Log likelihood ------------------------------------------------------------------------
addm_plot_loglik(addm_loglik_condition)
# --------------------------------------------------------------------------------------------

# Simulate Data
addm_data_full_output = addm_run_by_condition(choice.dat = addm_dat$choice.dat,
                                              conditions.dat = addm_dat$conditions.dat,
                                              model.parameters = c(0.0015,1,0.07,0),
                                              output.type = 'full',
                                              nr.reps = 250)
# --------------------------------------------------------------------------------------------

# Plot some psychometrics --------------------------------------------------------------------
addm2_plot_family(choice.dat = addm_dat$choice.dat,
                  addm.output = addm_data_full_output)
# --------------------------------------------------------------------------------------------

