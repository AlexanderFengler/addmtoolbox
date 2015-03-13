iandata <- addm_preprocess(choice.dat = addm_data_choice, eye.dat = addm_data_eye, rtbinsize = 100, timestep = 10)

ianaddmtrial <- addm_fit_grid(iandata,
                              drifts = seq(.0001,.0003,.0001),
                              thetas = seq(0,1,.1),
                              sds = seq(0.005,0.045,0.01),
                              non.decision.times = 350,
                              timestep = 10,
                              fit.type = 'trial')

ianloglik <- addm_plot_loglik(ianaddmtrial)

ianfulloutput <- addm_run_by_condition(choice.dat = iandata$choice.dat,
                                       conditions.dat = iandata$conditions.dat,
                                       nr.reps = 100,
                                       model.parameters = c(0.0003, 0.09, 0.045, 350),
                                       output.type = 'full')

ianplotfamily <- addm2_plot_family(choice.dat = iandata$choice.dat, addm.output = ianfulloutput)
