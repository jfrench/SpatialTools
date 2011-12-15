\name{condnorm.par} 

\alias{condnorm.par}

\title{
	Calculates conditional mean and covariance matrix of multivariate normal distribution
} 

\description{
	This function calculates the conditional mean and covariance matrix of a multivariate normal distribution based on observed responses.  The conditional mean and a decomposition of the conditional covariance matrix are returned.
} 

\usage{
	condnorm.par(y, V, Vp, Vop, coeff, X, Xp, method = "eigen")
} 

\arguments{
	\item{y}{A numeric vector of length \eqn{n} containing the observed responses.} 
	\item{V}{A numeric matrix of size \eqn{n \times n} containing the covariance of the observed responses.} 
	\item{Vp}{A numeric matrix of size \eqn{np \times np} containing the covariance of the unobserved responses.} 
	\item{Vop}{A numeric matrix of size \eqn{n \times np} containing the covariance between the observed and unobserved responses.} 
	\item{coeff}{A vector of length \eqn{k} containing the trend parameters.} 
	\item{X}{A numeric matrix of size \eqn{n \times k} containing the design matrix of the observed responses.} 
	\item{Xp}{A numeric matrix of size \eqn{np \times k} containing the design matrix of the unobserved responses.} 
	\item{method}{A character vector indicating the method of decomposition.  Current choices include eigen, chol, and svd.} 

} 

\details{
	It is assumed that the unconditional mean of the responses are given by \code{X} \%*\% \code{coeff} and \code{Xp} \%*\% \code{coeff}.  In other words, it is assumed that the mean structure of the response is linear in the coefficients.

The covariance matrix of the unobserved responses conditional on the observed responses in not returned.  Instead a decomposition of this matrix is returned.  Specifically, if \code{Vc} is the conditional covariance matrix, then \code{Vc} = \code{decomp.Vc} \%*\% \code{t(decomp.Vc)}.
}

\value{
	Returns a list with the following elements:
\item{mc}{ The conditional mean. }
\item{decomp.Vc}{ The decomposition of the conditional covariance matrix.}
}


%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{conditional} 
\keyword{normal} 

\seealso{ \code{decomp.cov} } 

\examples{
data(toydata)

y <- toydata$y
V <- toydata$V
Vp <- toydata$Vp
Vop <- toydata$Vop
X <- toydata$X
Xp <- toydata$Xp

temp1 <- condnorm.par(y, V, Vp, Vop, coeff = c(1, 2, 3), X, Xp, 
    method = "eigen")

#conditional covariance matrix
Vc <- temp1$decomp.Vc \%*\% t(temp1$decomp.Vc)
}