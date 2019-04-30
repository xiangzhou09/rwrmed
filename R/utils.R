demean <- function(x, w) x - weighted.mean(x, w, na.rm = TRUE)

`%||%` <- function(a, b) if (!is.na(a)) a else b
