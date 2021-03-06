% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/addm_run_by_condition.R
\name{addm_run_by_condition}
\alias{addm_run_by_condition}
\title{Run model by unique trial condition}
\usage{
addm_run_by_condition(choice.dat = data.table(decision = 0, rt = 0,
  condition_id = 0), conditions.dat = data.table(v1 = 1, v2 = 100,
  condition_id = 0), model.parameters = c(0, 0.002, 0.07, 0.5),
  nr.attributes = 1, boundaryfun = 1, nr.reps = 2000, timestep = 10,
  model.type = "standard", output.type = "fit",
  fixation.model = "fixedpath", generate = 0)
}
\arguments{
\item{choice.dat}{data.table storing the item valuations (v1,v2...), reaction times in ms (rt), decisions as decision and an id column (conditions_id).}

\item{conditions.dat}{data.table storing the item valuations (v1,v2...) by unique trial conditions. An id column (conditions_id) needs to be provided}

\item{model.parameters}{vector with the four core addm parameters in order (non.decision.time, drift, sd, theta, gamma, boundary-parameters).}

\item{nr.attributes}{integer providing the amount of attributes we consider per item}

\item{nr.reps}{integer that tells the function how many simulation runs to use.}

\item{timestep}{integer that provides the timestep-size that is used in the simulations (in ms).}

\item{model.type}{string that indicates which version of the model to run. 'standard' for standard (a)ddm model or 'memnoise' when memory effects shall be allowed.}

\item{output.type}{string that indicates what output the model shall produce. 'full' for detailed model output, 'fit' for sparse output (rt,decision) by id variable.}

\item{fixation.model}{string that indicates which fixation model will be utilized for simulations. 'random' for random fixations (supplied) 'fixedpath' for following a predetermined fixation path with fixed durations (supplied) or 'user' for a user supplied fixation model (function name: user_fixation_model).}

\item{generate}{boolean variable that tells the function to return either log likelihood values (0) or rt, decision (1). Relevant only if model.type variable is 'fit'.}
}
\value{
The function has three potential return values. A log likelihood value, utilized when trying to fit the model (output.type = 'fit', generate = 0). A simple data.table providing simulated rts,decisions by condition id, which useful for generating fake data in testing (output.type = 'fit', generate = 1). A full model output with many details, utilized for running the model with optimal parameters and extracting data for plots (output.type = 'full').
\code{addm_run_by_condition}
}
\description{
Runs model by unique trial condition for one set of parameter values
}
\author{
Alexander Fengler, \email{alexanderfengler@gmx.de}
}

