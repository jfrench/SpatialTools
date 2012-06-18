rmvnorm <- function(nsim = 1, mu, V, method = "eigen")
{
	mu <- as.vector(mu)
	
	# check arguments of function
	rmvnorm_arg_check(nsim, mu, V, method)
	
	# decompose covariance matrix
	decomp.V <- decomp.cov(V, method = method)
	
	# return simulated values
	return(mu + decomp.V %*% matrix(rnorm(nrow(V) * nsim), nrow = nrow(V), ncol = nsim))
}

rcondsim <- function(nsim = 1, krige.obj, y, V, Vp, Vop, Ve.diag, method = "eigen")
{
	# check arguments
	rcondsim_arg_check(nsim = nsim, y = y, V = V, Vp = Vp, 
		Vop = Vop, Ve.diag = Ve.diag, method = method, krige.obj = krige.obj)
	
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

	# Determine whether the krige.obj is from simple kriging.  If it is
	# then some slight adjustments need to be made to the simulation algorithm.
	if(is.null(krige.obj$mean))
	{
		return(newsim[-(1:n),] + crossprod(krige.obj$w, newZ))
	}else
	{
		return(newsim[-(1:n),] + krige.obj$mean + 
			crossprod(krige.obj$w, newZ - krige.obj$mean))
	}
}

condnorm.par <- function(y, V, Vp, Vop, coeff, X, Xp, method = "eigen")
{
	ViVop <- solve(V, Vop)
	
	#conditional mean of Yp given observed Yo = y
	mc <- Xp %*% coeff + crossprod(ViVop, y - X %*% coeff)
	Vc <- Vp - crossprod(Vop, ViVop)

	decomp.Vc <- decomp.cov(Vc, method = method)
	
	return(list(mc = mc, decomp.Vc = decomp.Vc))
}

