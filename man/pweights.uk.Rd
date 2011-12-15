\name{pweights.uk} 

\alias{pweights.uk}

\title{
	Calculates Universal Kriging weights
} 

\description{
	Finds Universal Kriging weights using \code{X}, 
	the \eqn{n \times k} 
	design matrix for the regression coefficients of the observed 
	data, 
	\code{V}, the (positive definite) covariance matrix of the 
	observed responses, \code{Xp}, the \eqn{np \times k} design matrix
	of the responses to be predicted, \code{Vp}, the 
	\eqn{np \times np}
	covariance matrix of the responses to be predicted, and \code{Vop}, 
	the \eqn{n x np} matrix of covariances between the observed 
	responses and the responses to be predicted.  Uses Armadillo C++ template via RcppArmadillo to perform most of the operations.
} 

\usage{
	pweights.uk(X, V, Xp, Vp, Vop)
} 

\arguments{
	\item{X}{The design matrix of the observed data.  
		 The size is \eqn{n \times k}} 
	\item{V}{The covariance matrix of the observed responses.  
		The size is \eqn{n times n}.} 
	\item{Xp}{The design matrix of the responses to be predicted.  
		The size is \eqn{np \times k}}.
	\item{Vp}{The covariance matrix of the responses to be predicted.
		The size is \eqn{np \times np}} 
	\item{Vop}{The cross-covariance between the observed responses 
		and the responses to be predicted.  The size is 
		\eqn{n \times np}} 
} 

\value{
The function returns a list containing \code{w}, the \eqn{np \times n} matrix containing the kriging weights used to calculate the predicted values. 
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{kriging}
\keyword{universal kriging}
\keyword{universal} 
\keyword{weights} 

\seealso{krige.uk}

\examples{
	data(toydata)
	X <- toydata$X
	V <- toydata$V
	Xp <- toydata$Xp
	Vp <- toydata$Vp
	Vop <- toydata$Vop
	pweights.uk(X, V, Xp, Vp, Vop)
}