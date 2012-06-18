\name{rcondsim} 

\alias{rcondsim}

\title{
	Generates realizations from a conditional normal distribution
} 

\description{
	Generates realizations from a conditional normal distribution.
} 

\usage{
	rcondsim(nsim = 1, krige.obj, y, V, Vp, Vop, Ve.diag, method = "eigen")
} 

\arguments{
	\item{nsim}{An integer indicating the number of realizations from the distribution.}
	\item{krige.obj}{An object returned by krige.sk, krige.ok, or krige.uk.} 
	\item{y}{A vector of length n containing the observed data values.} 
	\item{V}{The covariance matrix of the observed data.  The matrix should be symmetric and positive definite.  The size must be \eqn{n times n}.} 
	\item{Vp}{The covariance matrix of the responses to be predicted.  The matrix should be symmetric and positive definite.  The size is \eqn{np times np}.} 
	\item{Vop}{The cross-covariance matrix between the observed data and the responses to be predicted.  The size must be \eqn{n times np}.} 
	\item{Ve.diag}{A vector of length n specifying the measurement error variances of the observed data.} 
	\item{method}{The method for performing a decomposition of the covariance matrix.  Possible values are "eigen", "chol", and "svd", Eigen value decomposition, Cholesky decomposition, or Singular Value Decomposoition, respectively.} 
} 

\value{
An \eqn{n \times nsim} matrix containing the \code{nsim} realizations of the conditional realizations.  Each column of the matrix represents a realization of the conditional normal distribution.
}

\details{
If doing conditional simulation, the Cholesky decomposition should not work when there are coincident locations between the observed data locations and the predicted data locations.    Both the Eigen and Singular Value Decompositions should work.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{normal}
\keyword{conditional}

\seealso{rmvnorm}

\examples{
# create observed and predicted coordines
ocoords <- matrix(runif(100), ncol = 2)
pcoords <- matrix(runif(200), ncol = 2)

# include some observed locations in the predicted coordinates
acoords <- rbind(ocoords, pcoords)

# create design matrices
X <- as.matrix(cbind(1, ocoords))
Xa <- as.matrix(cbind(1, acoords))

# create covariance matrix
C3 <- cov.sp(coords = ocoords, sp.type = "matern", sp.par = c(2, 1), smoothness = 1, 
	finescale = 0, error = 0.5, pcoords = acoords)

# set values of regression coefficients
coeff <- matrix(c(1, 2, 3), nrow = 1)

# generate data with error
y <- rmvnorm(nsim = 1, mu = tcrossprod(X, coeff), V = C3$V) + rnorm(50, sd = sqrt(.5))

# use universal kriging to make predictions
krige.obj <- krige.uk(as.vector(y), V = C3$V, Vp = C3$Vp, Vop = C3$Vop, 
	X = X, Xp = Xa, return.w = TRUE)

# generate conditional simulations of observations
mycondsim <- rcondsim(nsim = 100, krige.obj = krige.obj, y = as.vector(y), 
	V = C3$V, Vp = C3$Vp, Vop = C3$Vop, 
	Ve.diag = rep(.5, 50), method = "eigen")
}