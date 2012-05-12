\name{krige.sk} 

\alias{krige.sk}

\title{
	Performs Simple Kriging
} 

\description{
	Performs Simple Kriging using \code{y}, a vector of length \eqn{n},
	\code{V}, the (positive definite) covariance matrix of the 
	observed responses, \code{Vp}, the 
	\eqn{np \times np}
	covariance matrix of the responses to be predicted, \code{Vop}, 
	the \eqn{n \times np} matrix of covariances between the observed 
	responses and the responses to be predicted, and \code{m}, a numeric vector 
  of length 1 identifying the value of the mean
  for each response.
} 

\usage{
	krige.sk(y, V, Vp, Vop, m = 0, return.w = FALSE, nsim = 0, 
	Ve.diag = NULL, method = "eigen")
} 

\arguments{
	\item{y}{The vector of observed responses.
		Should be a matrix of size \eqn{n \times 1} or a vector of 
		length \eqn{n}.}
	\item{V}{The covariance matrix of the observed responses.  
		The size is \eqn{n \times n}.} 
	\item{Vp}{The covariance matrix of the responses to be predicted.
		The size is \eqn{np \times np}} 
	\item{Vop}{The cross-covariance between the observed responses 
		and the responses to be predicted.  The size is 
		\eqn{n \times np}.} 
  \item{m}{A numeric vector of length 1 giving the mean of each response.}
	\item{return.w}{A boolean value (TRUE or FALSE) indicating whether the matrix of prediction weights should be returned.  By default, this is FALSE (and results in faster computation since there is less to return).}
	\item{nsim}{A non-negative value indicating the number of conditional simulations that should be returned.  If this is less than 1, then no conditional simulation is done.}
	\item{Ve.diag}{A vector of length n specifying the measurement error variances of the observed data.} 
	\item{method}{A character vector specifying the method used to decompose \code{V}.  Options are "eigen", "chol", or "svd" (Eigen decomposition, Cholesky decomposition, or Singular value decomposition, respectively).} 
} 

\value{
The function a list containing the following objects:
	\item{pred}{A vector of length \eqn{np} containing the predicted
		responses.}
	\item{mspe}{A vector of length \eqn{np} containing the 
		mean-square prediction error of the predicted
		responses.}
	\item{w}{A \eqn{np \times n} matrix containing the kriging
		weights used to calculate red.}
	\item{simulations}{An \eqn{n \times nsim} matrix containing the \code{nsim} realizations of the conditional realizations.  Each column of the matrix represents a realization of the conditional normal distribution.}
	\item{mean}{The mean value (m) originally provided to the function}.  
}

\details{
It is assumed that there are \eqn{n} observed data values
and that we wish to make predictions at \eqn{np} locations.
The mean is subtracted from each value of \code{y} before determining the kriging weights, 
and then the mean is added onto the predicted response.
}
\references{ Statistical Methods for Spatial Data Analysis, Schabenberger and Gotway (2003).  See p. 226-228.} 

\author{Joshua French} 

\keyword{kriging}
\keyword{simple kriging}
\keyword{simple} 

\examples{
	data(toydata)
	y <- as.vector(toydata$y)
	V <- toydata$V
	Vp <- toydata$Vp
	Vop <- toydata$Vop
	krige.sk(y, V, Vp, Vop, m = 2)
}