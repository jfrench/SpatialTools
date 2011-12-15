\name{mspe.uk} 

\alias{mspe.uk}

\title{
	Calculates mean-square prediction error of Universal Kriging predictions
} 

\description{
	Calculates mean-square prediction error of Universal Kriging 
	predictions using \code{w}, 
	the \eqn{n \times np} 
	matrix of kriging weights, 
	\code{V}, the (positive definite) covariance matrix of the 
	observed responses, \code{Vp}, the (positive definite) 
	covariance matrix of the 
	predicted responses,and \code{Vop}, 
	the \eqn{n x np} matrix of covariances between the observed 
	responses and the responses to be predicted.  Uses Armadillo C++ template via RcppArmadillo to perform most of the operations.
} 

\usage{
	mspe.uk(w, V, Vp, Vop)
} 

\arguments{
	\item{w}{The matrix of kriging weights for the predicted
		responses.  The size is \eqn{n \times np}} 
	\item{V}{The covariance matrix of the observed responses.  
		The size is \eqn{n times n}.} 
	\item{Vp}{The covariance matrix of the responses to be
		predicted.  The size is \eqn{np times np}.} 
	\item{Vop}{The cross-covariance between the observed responses 
		and the responses to be predicted.  The size is 
		\eqn{n \times np}} 
} 

\value{
The function returns a list containing \code{mspe}, a matrix of size \eqn{n \times 1} containing the mean-square prediction error of the predicted responses. 
}


%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{kriging}
\keyword{universal kriging}
\keyword{universal} 
\keyword{mean-square prediction error} 

\seealso{krige.uk, weights.uk}

\examples{
	data(toydata)
	X <- toydata$X
	V <- toydata$V
	Xp <- toydata$Xp
	Vp <- toydata$Vp
	Vop <- toydata$Vop
	w <- pweights.uk(X, V, Xp, Vp, Vop)
	mspe.uk(w, V, Vp, Vop)
}