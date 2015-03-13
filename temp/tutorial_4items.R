# Read in my addm project data

afdat_eye4 = read.table(file = 'temp/data/eye4.txt',header = TRUE, sep = ' ')

afdat_choice4 = read.table(file = 'temp/data/choice4.txt',header = TRUE, sep = ' ')

# Preprocess my data

afdat = addm_preprocess(choice.dat = afdat_choice4, eye.dat = afdat_eye4)

afloglik = addm_fit_grid(afdat, drifts = c(0.004,0.005,0.006), sds = c(0.08, 0.06, 0.04, 0.02), theta = c(0, 0.5, 1))

addm_plot_loglik(afloglik)
