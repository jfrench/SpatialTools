#Take the list of contours from contourLines() and extracts coordinates 
get.contours=function(contours.list)
{
	x <- y <- NULL
	for(i in 1:length(contours.list))
	{
		x <- c(x, contours.list[[i]]$x)
		y <- c(y, contours.list[[i]]$y)
	}
	cbind(x,y)
}

#Plot contour lines from contourLines
plot.contourLines <- function(x, begin=1, end = length(x), ...)
{
	for(i in begin:end)
	{
		lines(x[[i]]$x, x[[i]]$y, ...)
	}
}
