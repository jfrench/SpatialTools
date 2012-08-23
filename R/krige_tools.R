krige.uk <- function(y, V, Vp, Vop, X, Xp, nsim = 0, Ve.diag = NULL, method = "eigen")
{
	# check arguments, create appropriate values of rws and method for .Call
	ins <- krige_arg_check(y, V, Vp, Vop, X, Xp, m = 0, nsim, Ve.diag, method)

	out <- .Call( "krige_uk", ys = y, Vs = V, Vps = Vp, Vops = Vop, Xs = X, Xps = Xp,
		nsims = nsim, Vediags = ins$Ve.diag, method = ins$method, 
		PACKAGE = "SpatialTools")

	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	out$coeff <- as.vector(out$coeff)
	return(out)
}

krige.ok <- function(y, V, Vp, Vop, nsim = 0, Ve.diag = NULL, method = "eigen")
{
	# check arguments, create appropriate values of rws and method for .Call
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = 0, nsim, Ve.diag, method)

	out <- .Call( "krige_ok", ys = y, Vs = V, Vps = Vp, Vops = Vop, 
		nsims = nsim, Vediags = ins$Ve.diag, method = ins$method, 
		PACKAGE = "SpatialTools")

	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	out$coeff <- as.vector(out$coeff)
	return(out)
}

krige.sk <- function(y, V, Vp, Vop, m = 0, nsim = 0, Ve.diag = NULL, method = "eigen")
{
	# check arguments, create appropriate values of rws and method for .Call
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = 0, nsim, Ve.diag, method)

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