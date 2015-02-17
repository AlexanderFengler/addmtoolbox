install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
library(devtools)
has_devel()


# Initialize Package with

devtools::create("path/to/package/pkgname")

# ignoring files when building a package from source
devtools::use_build_ignore("notes")

# make binary package
devtools::build(binary = TRUE)


# editing and then reloading code one can use
devtools::load_all()

# reset options to old with
on.exit()


# do something when package is loaded
# both normally saved in a file called zzz.R
.onAttach()
.onLoad()

# line out order of sourcing code when package is loaded with "Collate" field in the DESCRIPTION


# Using external packages (by importing them if not present) we also utilitze the decription file
# Hadley provides nice function here to make things easier

devtools::use_package("dplyr")




# set up package to work with rcpp
devtools::use_rcpp()


# making cpp classes available
Rcpp::loadRcppModules()


# Define Namespace in NAMESPACE FILE

# to require namespaces (packages in searchpath at time of execution) in packages use:
requireNamespace("x", quietly = TRUE)

# To import a function fomr a package (and avoid namespace search to save approxiamtely 5microsecons per call....)
@importFrom pgk fun

# to use data internally (store some tables or so)
devtools::use_data() to create this file with the argument internal = TRUE

