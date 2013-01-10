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