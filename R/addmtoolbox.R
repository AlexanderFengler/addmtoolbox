#' addmtoolbox: Fit (attentional) drift diffusion models
#'
#' The addmtoolbox package provides four categories of functions:
#' addm_dataprep(), addm_run_...(), addm_fit_...(), addm_plot_...().
#'
#' @section addm_dataprep():
#' This function gets your data in shape for easy usage with the other functions in the package.
#' It is not strictly necessary to use it, however strongly recommended for predictable behavior.
#'
#' @section addm_run_...():
#' These are wrapper functions to run a (a)ddm simulations for a set of provided parameters.
#'
#' @section addm_fit_...():
#' These functions allow you to do maximum likelihood fits of the model, provided with a parameter space. Currently only grid search is possible,
#' however new fitting methods will be developed in the near future.
#'
#' @section addm_plot_...() functions:
#' This set of functions will enable you to plot various psychometrics (currently two plots) and diagnostics of the modelfits.
#'
#' @docType package
#' @name addmtoolbox
#' @import data.table dplyr Rcpp RcppGSL RcppZiggurat parallel doMC foreach iterators ggplot2 pander rmarkdown
#' @useDynLib addmtoolbox
#' @importFrom Rcpp sourceCpp
NULL
