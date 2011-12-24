krige.uk <- function(y, V, Vp, Vop, X, Xp)
{
	y <- as.vector(y)
	krige_arg_check(y, V, Vp, Vop, X = X, Xp = Xp, coeff = NULL)

	out <- .Call( "krige_uk", Xs = X, ys = y, Vs = V, Xps = Xp,
		Vps = Vp, Vops = Vop, PACKAGE = "SpatialTools")

	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	out$coeff <- as.vector(out$coeff)
	return(out)
}

krige.ok <- function(y, V, Vp, Vop)
{
	y <- as.vector(y)
	krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, coeff = NULL)

	out <- .Call( "krige_ok", ys = y, Vs = V, 
		Vps = Vp, Vops = Vop, PACKAGE = "SpatialTools")
		
	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	out$coeff <- as.vector(out$coeff)
	return(out)
}

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

krige.sk <- function(y, V, Vp, Vop, m = 0)
{
	y <- as.vector(y)
	m <- as.vector(m)
	
	krige_sk_arg_check(y, V, Vp, Vop, m)
	
	out <- .Call( "krige_sk", ys = y, Vs = V, 
		Vps = Vp, Vops = Vop, ms = m, PACKAGE = "SpatialTools")

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

