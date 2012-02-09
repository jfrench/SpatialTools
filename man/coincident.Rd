\name{coincident} 

\alias{coincident}

\title{
	Determine coincident locations
} 

\description{
	Determines the coincident locations for two sets of coordinates.
} 

\usage{
	coincident(coords1, coords2)
} 

\arguments{
	\item{coords1}{A numeric matrix of size \eqn{n1 \times 2} containing a set of coordinates.} 
	\item{coords2}{A numeric matrix of size \eqn{n2 \times 2} containing a set of coordinates.} 	
} 

\value{
	Returns a matrix with the indices of the coincident locations.  Specifically, an \eqn{r \times 2} matrix will be returned, with \eqn{r} being the number of coordinates in \code{coords1} coinciding with coordinates in \eqn{coords2}.  If row i of the matrix was (2, 5), then the ith set of coincident locations was between the 2nd row of \code{coords1} and the 5th row of \code{coords2}.  If there are no coincident locations, then a matrix of size \eqn{0 \times 2} will be returned.
}

\details{
	This function calls a C++ algorithm using the Rcpp package.  This (may) result in a significant speedup over pure R code since the search algorithm involves loops.  An implicit assumptions is that there are no duplicate coordinates in \code{coords1} and \code{coords2}, respectively.  In other words, each row of \code{coords1} is unique and each row of \code{coords2} is unique.  There is at most 1 row of \code{coords1} that will match with a row in \code{coords2}.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{duplicate}
\keyword{coords}

\seealso{ \code{cov.st} } 

\examples{
#Generate two sets of coordinates
loc1 <- as.matrix(expand.grid(seq(0, 1, len = 25), seq(0, 1, len = 25)))
loc2 <- as.matrix(expand.grid(seq(0, 1, len = 101), seq(0, 1, len = 101)))

#Return indices of coincident coordinates
coincident(loc1, loc2)


}