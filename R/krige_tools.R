krige.uk <- function(y, V, Vp, Vop, X, Xp, ...)
{
	# Check arguments.  Create unspecified arguments if needed.
	nsim <- 0; Ve.diag <- NULL; method <- "eigen"; level <- NULL; alternative <- NULL
	arglist <- list(...)
	argnames <- names(arglist)
	if("nsim" %in% argnames){ nsim <- arglist$nsim }
	if("Ve.diag" %in% argnames){ Ve.diag <- arglist$Ve.diag }
	if("method" %in% argnames){ method <- arglist$method }
	if("level" %in% argnames){ level <- arglist$level }
	if("alternative" %in% argnames){ alternative <- arglist$alternative }
	ins <- krige_arg_check(y, V, Vp, Vop, X, Xp, m = 0, nsim, Ve.diag, method, 
		level, alternative)

	###compute matrix products for future use
	ViX <- solve(V, X)
	XtViX <- crossprod(ViX, X)
	
	#compute gls estimates of regression coefficients
	coeff <- solve(XtViX, crossprod(ViX, y))

	#compute kriging weights
	w <- solve(V, Vop - X%*%solve(XtViX, crossprod(X, solve(V, Vop)) - t(Xp)))

	#blup for yp
	pred <- crossprod(w, y)
	
	#variance of (yp - pred)
	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)

	out <- list(pred = pred, mspe = mspe, coeff = coeff, 
		vcov.coef = solve(XtViX))

	# generate conditional realizations if nsim > 0
	if(nsim > 0)
	{
		# determine number of observed and predicted data values
		n <- nrow(V); np <- nrow(Vp)

		# modify V to account for error
		Vomod <- V - diag(Ve.diag)
		
		# Combine observed and predicted covariance matrices
		Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))
	
		# Generate mean zero normal random variables with appropriate covariance
		# matrix (does not include error)
		newsim <- decomp.cov(Va, method = method) %*% 
			matrix(rnorm(nrow(Va) * nsim), ncol = nsim)
	
		# Take difference of observed data and simulated observations 
		# (now including error)
		newZ <- matrix(y, nrow = n, ncol = nsim) - newsim[1:n,] + 
			matrix(rnorm(n * nsim, sd = sqrt(Ve.diag)), nrow = n, ncol = nsim)
	
		# Create conditional realizations
		sim <- newsim[-(1:n),] + crossprod(w, newZ)
		
		out$sim <- sim
	}
	return(out)
}

krige.ok <- function(y, V, Vp, Vop, ...)
{
	# Check arguments.  Create unspecified arguments if needed.
	nsim <- 0; Ve.diag <- NULL; method <- "eigen"; level <- NULL; alternative <- NULL
	arglist <- list(...)
	argnames <- names(arglist)
	if("nsim" %in% argnames){ nsim <- arglist$nsim }
	if("Ve.diag" %in% argnames){ Ve.diag <- arglist$Ve.diag }
	if("method" %in% argnames){ method <- arglist$method }
	if("level" %in% argnames){ level <- arglist$level }
	if("alternative" %in% argnames){ alternative <- arglist$alternative }
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = 0, nsim, Ve.diag, method, 
		level, alternative)

	out <- .Call( "krige_ok", ys = y, Vs = V, Vps = Vp, Vops = Vop, 
		nsims = nsim, Vediags = ins$Ve.diag, method = ins$method, 
		PACKAGE = "SpatialTools")

	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	out$coeff <- as.vector(out$coeff)
	return(out)
}

krige.sk <- function(y, V, Vp, Vop, m = 0, ...)
{
	# Check arguments.  Create unspecified arguments if needed.
	nsim <- 0; Ve.diag <- NULL; method <- "eigen"; level <- NULL; alternative <- NULL
	arglist <- list(...)
	argnames <- names(arglist)
	if("nsim" %in% argnames){ nsim <- arglist$nsim }
	if("Ve.diag" %in% argnames){ Ve.diag <- arglist$Ve.diag }
	if("method" %in% argnames){ method <- arglist$method }
	if("level" %in% argnames){ level <- arglist$level }
	if("alternative" %in% argnames){ alternative <- arglist$alternative }
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = m, nsim, Ve.diag, method, 
		level, alternative)

	out <- .Call( "krige_sk", ys = y, Vs = V, Vps = Vp, Vops = Vop, ms = m, 
		nsims = nsim, Vediags = ins$Ve.diag, method = ins$method, 
		PACKAGE = "SpatialTools")
	
	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	return(out)
}

pweights.uk <- function(X, V, Xp, Vp, Vop)
{
	pweights_uk_arg_check(X, V, Xp, Vp, Vop)

	.Call( "pweights_uk", Xs = X, Vs = V, Xps = Xp,
		Vps = Vp, Vops = Vop, PACKAGE = "SpatialTools")
}

mspe.uk <- function(w, V, Vp, Vop)
{
	mspe_uk_arg_check(w, V, Vp, Vop)

	as.vector(.Call( "mspe_uk", ws = w, Vs = V, 
		Vps = Vp, Vops = Vop, PACKAGE = "SpatialTools"))
}