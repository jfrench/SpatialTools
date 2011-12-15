\name{krige.uk} 

\alias{krige.uk}

\title{
	Performs universal kriging
} 

\description{
	Performs Universal Kriging using \code{X}, the \eqn{n \times k} 
	design matrix for the regression coefficients of the observed 
	data, \code{y}, the \eqn{n \times 1} matrix of observed responses,
	\code{V}, the (positive definite) covariance matrix of the 
	observed responses, \code{Xp}, the \eqn{np \times k} design matrix
	of the responses to be predicted, \code{Vp}, the 
	\eqn{np \times np}
	covariance matrix of the responses to be predicted, and \code{Vop}, 
	the \eqn{n \times np} matrix of covariances between the observed 
	responses and the responses to be predicted.  Uses Armadillo C++ template via RcppArmadillo to perform most of the operations.
} 

\usage{
	krige.uk(y, V, Vp, Vop, X, Xp)
} 

\arguments{
	\item{y}{The vector of observed responses.  Should be a matrix of size \eqn{n \times 1} or a vector of length \eqn{n}.}
	\item{V}{The covariance matrix of the observed responses.  
		The size is \eqn{n times n}.} 
	\item{Vp}{The covariance matrix of the responses to be predicted.  The size is \eqn{np \times np}} 
	\item{Vop}{The cross-covariance between the observed responses and the responses to be predicted.  The size is \eqn{n \times np}} 
	\item{X}{The design matrix of the observed data.  The size is \eqn{n \times k}} 
	\item{Xp}{The design matrix of the responses to be predicted.  
		The size is \eqn{np \times k}}.
} 

\value{
The function a list containing the following objects:
	\item{pred}{A vector of length \eqn{np} containing the predicted
		responses.}
	\item{mope}{A vector of length \eqn{np} containing the 
		mean-square prediction error of the predicted
		responses.}
	\item{w}{A \eqn{np \times n} matrix containing the kriging
		weights used to calculate \code{red}.} 
	\item{coeff}{A vector of length \eqn{k} containing the estimated
		regression coefficients.} 
	\item{vcov.coeff}{A \eqn{k \times k} matrix containing the
	(estimated) covariance matrix of estimated the regression coefficients.} 
}

\details{

It is assumed that there are \eqn{n} observed data values
and that we wish to make predictions at \eqn{np} locations.  We assume
that there are \eqn{k} regression coefficients (including the intercept).
Both \code{X} and \code{Xp} should contain a column of 1's if an intercept
is desired.

}
\references{ Statistical Methods for Spatial Data Analysis, Schabenberger and Gotway (2003).  See p. 241-243.} 

\author{Joshua French} 

\keyword{kriging}
\keyword{universal kriging}
\keyword{universal} 

\examples{
	data(toydata)
	X <- toydata$X
	y <- toydata$y
	V <- toydata$V
	Xp <- toydata$Xp
	Vp <- toydata$Vp
	Vop <- toydata$Vop
	krige.uk(y = y, V = V, Vp = Vp, Vop = Vop, X = X, Xp = Xp)
}