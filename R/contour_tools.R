#Take the list of contours from contourLines() and extracts coordinates 
get.contours=function(x)
{
	x1 <- x2 <- NULL
	for(i in 1:length(x))
	{
		x1 <- c(x1, x[[i]]$x)
		x2 <- c(x2, x[[i]]$y)
	}
	cbind(x = x1,y = x2)
}

#Plot contour lines from contourLines
plot.contourLines <- function(x, begin=1, end = length(x), add = FALSE, ...)
{
	if(!add)
	{
		contours <- get.contours(x)
		rx <- range(contours[,1])
		ry <- range(contours[,2])
		plot(rx, ry, type = "n", ...)
	}
	for(i in begin:end)
	{
		lines(x[[i]]$x, x[[i]]$y, ...)
	}
}
