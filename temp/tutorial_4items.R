




# Multiple items
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

# SET SIZE 4
# Read in my addm project data -----------------------------------------------------------------------
afdat_eye4 = as.data.table(read.table(file = 'temp/data/eye4.txt',header = TRUE, sep = ' '))
afdat_choice4 = as.data.table(read.table(file = 'temp/data/choice4.txt',header = TRUE, sep = ' '))
# ----------------------------------------------------------------------------------------------------

rows.usable = createDataPartition(1:length(afdat_choice4$id),times=1, p=0.1)
afdat_choice = afdat_choice4[rows.usable[[1]],]

setkey(afdat_choice, id)
setkey(afdat_eye4, id)

afdat_eye = afdat_choice[afdat_eye4]
afdat_eye = afdat_eye[complete.cases(afdat_eye)] %>% select(id,fixloc,fixdur,fixnr)


# Preprocess data
afdat = addm_preprocess(choice.dat = afdat_choice, eye.dat = afdat_eye)

# Run grid search
afloglik = addm_fit_grid(afdat,
                         drifts = c(0.005,0.006,0.007),
                         sds = c(0.08, 0.06, 0.04),
                         theta = c(0,0.5,1),
                         fit.type = 'trial')
addm_plot_loglik(afloglik)
# ----------------------------------------------------------------------------------------------------

# SET SIZE 8
# Read in my addm project data -----------------------------------------------------------------------
afdat_eye8 = as.data.table(read.table(file = 'temp/data/eye8.txt',header = TRUE, sep = ' '))
afdat_choice8 = as.data.table(read.table(file = 'temp/data/choice8.txt',header = TRUE, sep = ' '))
# ----------------------------------------------------------------------------------------------------

rows.usable = createDataPartition(1:length(afdat_choice8$id),times=1, p=0.05)
afdat_choice = afdat_choice8[rows.usable[[1]],]

setkey(afdat_choice, id)
setkey(afdat_eye8, id)

afdat_eye = afdat_choice[afdat_eye8]
afdat_eye = afdat_eye[complete.cases(afdat_eye)] %>% select(id,fixloc,fixdur,fixnr)


# Preprocess my data
afdat = addm_preprocess(choice.dat = afdat_choice, eye.dat = afdat_eye)

# Run Grid Search
afloglik = addm_fit_grid(afdat,
                         drifts = c(0.005,0.006,0.007),
                         sds = c(0.08, 0.06, 0.04),
                         theta = c(0,0.5,1),
                         fit.type = 'trial')
addm_plot_loglik(afloglik)
# ----------------------------------------------------------------------------------------------------

