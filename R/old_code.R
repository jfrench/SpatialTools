#krige.uk2 <- function(X, y, V, Xp, Vp, Vop)
#{
#
#	###compute matrix products for future use
#	#compute Vi*X
#	ViX <- solve(V, X)
#	#compute X'*Vi*X
#	XtViX <- crossprod(ViX, X)
#	
#	#compute gls estimates of regression coefficients
#	coef <- solve(XtViX, crossprod(ViX, y))
#
#	#compute kriging weights
#	w <- solve(V, Vop - X%*%solve(XtViX, crossprod(X, solve(V, Vop)) - t(Xp)))
#
#	#blup for Yp
#	pred <- crossprod(w, y)
#	
#	#variance of (Yp - pred)
#	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)
#
#	out <- list(pred = pred, mspe = mspe, w = w, coef = coef, 
#		vcov.coef = solve(XtViX))
#
#	class(out) <- "uk"
#	class(out) <- "krige"
#	return(out)
#}

# krige.sk2 <- function(y, V, Vp, Vop, m = 0)
# {
# 	y <- as.vector(y)
# 	m <- as.vector(m)
# 	
# 	krige_sk_arg_check(y, V, Vp, Vop, m)
# 	
# 	w <- solve(V, Vop)
# 	pred <- m + crossprod(w, y - m)
# 	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)
# 
# 	out <- list(pred = pred, mspe = mspe, w = w)
# 	return(out)
# }



#Generally slower than R based rmvnorm
#rmvnorm2 <- function(nsim = 1, mu, V, method = "eigen")
#{
#	mu <- as.vector(mu)
#	rmvnorm_arg_check(nsim, mu, V, method)
#	if(method == "eigen")
#	{
#		.Call( "rmvnorm", nsims = nsim, mus = mu, Vs = V, methods = 1,
#			PACKAGE = "SpatialTools")
#	}else if(method == "chol")
#	{
#		.Call( "rmvnorm", nsims = nsim, mus = mu, Vs = V, methods = 2,
#			PACKAGE = "SpatialTools")
#	}else
#	{
#		.Call( "rmvnorm", nsims = nsim, mus = mu, Vs = V, methods = 3,
#			PACKAGE = "SpatialTools")
#	}
#}

#condnorm.par2 <- function(y, V, Vp, Vop, coeff, X, Xp, method = "eigen")
#{
#	if(method == "eigen")
#	{
#		.Call( "condnorm_par", ys = y, Vs = V, Vps = Vp, Vops = Vop, coeffs = coeff,
#			Xs = X, Xps = Xp, methods = 1, PACKAGE = "SpatialTools")
#	}else if(method == "chol")
#	{
#		.Call( "condnorm_par", ys = y, Vs = V, Vps = Vp, Vops = Vop, coeffs = coeff,
#			Xs = X, Xps = Xp, methods = 2, PACKAGE = "SpatialTools")
#	}else
#	{
#		.Call( "condnorm_par", ys = y, Vs = V, Vps = Vp, Vops = Vop, coeffs = coeff,
#			Xs = X, Xps = Xp, methods = 3, PACKAGE = "SpatialTools")
#	}
#}
#
# 
# rcondsim <- function(nsim = 1, krige.obj, y, V, Vp, Vop, Ve.diag, method = "eigen")
# {
# 	# check arguments
# 	rcondsim_arg_check(nsim = nsim, y = y, V = V, Vp = Vp, 
# 		Vop = Vop, Ve.diag = Ve.diag, method = method, krige.obj = krige.obj)
# 	
# 	# determine number of observed and predicted data values
# 	n <- nrow(V); np <- nrow(Vp)
# 
# 	# modify V to account for error
# 	Vomod <- V - diag(Ve.diag)
# 	
# 	# Combine observed and predicted covariance matrices
# 	Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))
# 
# 	# Generate mean zero normal random variables with appropriate covariance
# 	# matrix (does not include error)
# 	newsim <- decomp.cov(Va, method = method) %*% 
# 		matrix(rnorm(nrow(Va) * nsim), ncol = nsim)
# 
# 	# Take difference of observed data and simulated observations 
# 	# (now including error)
# 	newZ <- matrix(y, nrow = n, ncol = nsim) - newsim[1:n,] + 
# 		matrix(rnorm(n * nsim, sd = sqrt(Ve.diag)), nrow = n, ncol = nsim)
# 
# 	# Determine whether the krige.obj is from simple kriging.  If it is
# 	# then some slight adjustments need to be made to the simulation algorithm.
# 	if(is.null(krige.obj$mean))
# 	{
# 		return(newsim[-(1:n),] + crossprod(krige.obj$w, newZ))
# 	}else
# 	{
# 		return(newsim[-(1:n),] + krige.obj$mean + 
# 			crossprod(krige.obj$w, newZ - krige.obj$mean))
# 	}
# }
# 
# condnorm.par <- function(y, V, Vp, Vop, coeff, X, Xp, method = "eigen")
# {
# 	ViVop <- solve(V, Vop)
# 	
# 	#conditional mean of Yp given observed Yo = y
# 	mc <- Xp %*% coeff + crossprod(ViVop, y - X %*% coeff)
# 	Vc <- Vp - crossprod(Vop, ViVop)
# 
# 	decomp.Vc <- decomp.cov(Vc, method = method)
# 	
# 	return(list(mc = mc, decomp.Vc = decomp.Vc))
# }
