.onAttach <- function(...) {
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

	greet <- paste("#", "# SpatialTools Package v.0.4.3", "# 2012-09-05", 
	paste("# Copyright (C) 2011-", yr, ", Joshua P. French", sep = ""), 
	"# Written by Joshua P. French", 
	"# This research was partially supported under NSF Grant ATM-0534173", "#", sep = "\n")
	packageStartupMessage(greet)
}

