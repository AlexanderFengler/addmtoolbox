# Lets Prepare our data
addm.dat  = addm_dataprep(choice.dat = addm_data_choice, eye.dat = addm_data_eye, rtbinsize = 100, timestep = 10)

# Run model by trial
addm_data_loglik_trial = addm_fit_grid(addm.dat, fit.type = 'trial',sds = 0.07)
