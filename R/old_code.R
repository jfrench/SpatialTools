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