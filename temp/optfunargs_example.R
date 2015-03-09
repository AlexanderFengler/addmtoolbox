monkeyfun = function(...){
  return(match.call())
  params <- as.list(match.call()[-1])
  if ("monkey" %in% names(params)){  ## note how I change the test here
    return(params)
    return('monkey found')
  } else {
    return('monkey not found')
  }
}

monkeyfun <- function(...) {
  if (hasArg("monkey")) {
    return('monkey found')
  } else {
    return('monkey not found')
  }
}
monkeyfun()
# [1] "monkey not found"
monkeyfun(monkey=0)
# [1] "monkey found"
