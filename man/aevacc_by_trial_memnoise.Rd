% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{aevacc_by_trial_memnoise}
\alias{aevacc_by_trial_memnoise}
\title{Simulate aDDM process by unique trial (>2 items, allow memory effects)}
\usage{
aevacc_by_trial_memnoise(parameters, maxdur, mindur, cur_decision, update,
  nr_attributes, fixpos, fixdur, nr_reps, timestep)
}
\arguments{
\item{parameters}{vector that stores the parameters used for the simulations (Order: [non.decision.time, drift, sd, theta, gamma, scalar_items_seen_drift, scala_items_seen_noise])}

\item{maxdur}{numeric variable that supplies the maximum reaction time considered a success in simulations}

\item{mindur}{numeric variable that supplies the minimum reaction time considered a succes in simulations}

\item{cur_decision}{numeric variable that provides the empirical decision taken in trial}

\item{update}{Vector that stores the item valuations for the trial conditon simulated}

\item{fixpos}{Vector that stores the locations for a supplied fixed fixation pathway}

\item{fixdur}{Vector that stores the fixation durations for a supplied fixed fixation pathway}

\item{nr_reps}{number of repitions (simulation runs)}

\item{timestep}{timestep in ms associated with each step in the drift diffusion process}

\item{items_seen_bias}{Numeric Variable storing the relative amount of drift that unseen items receive}

\item{items_seen_noise_bias}{Numeric Variable storing the relative noise sd that unseen items receive}
}
\value{
numeric variable that provides a success count (runs that predicted a reaction time in the correct rt-bin and simultaneously the correct decision)
}
\description{
Simulate aDDM process by unique trial (>2 items, allow memory effects)
\code{aevacc_by_trial_memnoise}
}
\author{
Alexander Fengler, \email{alexanderfengler@gmx.de}
}

