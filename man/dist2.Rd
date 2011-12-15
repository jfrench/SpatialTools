\name{dist2} 

\alias{dist2}

\title{
	Calculates Euclidean distance between coordinates of two matrices
} 

\description{
	Calculates Euclidean distance between the coordinates of \code{coords1} 
	and \code{coords2} matrices. It is assumed that \code{coords1} and 
	\code{coords2} have the same number of columns but not 
	necessarily the same number of rows.
} 

\usage{
	dist2(coords1, coords2)
} 

\arguments{
	\item{coords1}{A matrix of dimensions \eqn{nr1\times nc}} 
	\item{coords2}{A matrix of dimensions \eqn{nr2\times nc}} 
} 

\value{
  \code{dist2} returns a matrix of size \eqn{nr1\times nr2} containing the 
  Eucliean distances between
  each pair of coordinates in \code{coords1} and \code{coords2}.
}


%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{distance} 
\keyword{Euclidean} 

\seealso{ \code{dist1} } 

\examples{
	x1 <- matrix(rnorm(30), ncol = 3)
	x2 <- matrix(rnorm(60), ncol = 3)
	dist2(x1, x2)
}