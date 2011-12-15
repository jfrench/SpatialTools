\name{dist1} 

\alias{dist1}

\title{
	Calculates Euclidean distance of matrix of coordinates
} 

\description{
	Calculates Euclidean distance of \code{coords} matrix of coordinates.  Each row of the 
	\code{coords}
	matrix is assumed to be a distinct coordinate.
} 

\usage{
	dist1(coords)
} 

\arguments{
	\item{coords}{A matrix of dimensions \eqn{nr\times nc}} 
} 

\value{
  \code{dist1} returns a matrix of size \eqn{nr\times nr} containing the Euclidean 
  distances between each pair of coordinates.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{distance} 
\keyword{Euclidean} 

\seealso{\code{dist2} } 

\examples{
	x <- matrix(rnorm(30), ncol = 3)
	dist1(x)
}