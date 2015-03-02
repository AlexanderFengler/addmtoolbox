#' addmtoolbox: A package to fit (attentional) drift diffusion models
#'
#' The addmtoolbox package provides four categories of functions:
#' addm_dataprep(), addm_plot(), addm_by..(), addm_fit..().
#'
#' @section addm_dataprep():
#' This function gets your data in shape for usage with the rest of the package
#'
#' @section addm_plot() functions:
#' This set of functions will enable you to plot various psychometrics (in development).
#' Currently there is only one function available that lets you plot the log likelihoods
#' for each parameter combination from your model fit.
#'
#' @section addm_by..():
#' These are wrapper functions to run a (a)ddm simulation provided with some parameters.
#'
#' @section addm_fit():
#' These functions allow you to fit the model to a provided parameter space (currently only grid_search).
#'
#' @docType package
#' @name addmtoolbox
#' @import data.table dplyr Rcpp RcppGSL RcppZiggurat doMC foreach iterators ggplot2 pander
#' @useDynLib addmtoolbox
#' @importFrom Rcpp sourceCpp
NULL
