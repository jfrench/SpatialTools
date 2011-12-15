.onAttach <- function(...) {
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
    cat("##\n## SpatialTools Package v.0.3-3\n")
    cat("## 2011-12-15\n")
    cat("## Copyright (C) 2011-", yr, ", Joshua P. French\n", sep="")
    cat("## Written by Joshua P. French\n")
    cat("##\n## This research was partially supported under NSF Grant ATM-0534173\n")
}