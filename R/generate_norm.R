rmvnorm <- function(nsim = 1, mu, V, method = "eigen")
{
	mu <- as.vector(mu)
	rmvnorm_arg_check(nsim, mu, V, method)
	decomp.V <- decomp.cov(V, method = method)
	return(mu + decomp.V %*% matrix(rnorm(nrow(V) * nsim), nrow = nrow(V), ncol = nsim))
}

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

rmvnorm_arg_check <- function(nsim, mu, V, method)
{
	if(!is.numeric(nsim) || !is.numeric(mu) || ! is.numeric(V))
	{
		stop("nsim, mu, and V arguments must all be numeric")
	}
	if(!isSymmetric(V) || !is.matrix(V))
	{
		stop("V must be a symmetrix matrix")
	}
	if(length(mu) != nrow(V))
	{
		stop("The length of mu must equal nrows of V")
	}
	if(!(method == "eigen" || method == "chol" || method == "svd"))
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
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

condnorm_par_arg_check <- function(y, V, Vp, Vop, coeff, X, Xp, method)
{
	n <- length(y)
	nk <- length(coeff)

	if(!is.numeric(y))
	{
		stop("y must be a numeric vector")
	}
	if(!is.matrix(V) || !is.numeric(V) || (nrow(V)!= ncol(V)))
	{
		stop("V must be a square numeric matrix")
	}
	if(!is.matrix(Vp) || !is.numeric(Vp) || (nrow(Vp)!= ncol(Vp)))
	{
		stop("Vp must be a square numeric matrix")
	}
	if(!is.matrix(Vop) || !is.numeric(Vop))
	{
		stop("Vop must be a numeric matrix")
	}
	if(length(y) != nrow(V))
	{
		stop("length of y must match nrows of V")
	}
	if(length(y) != nrow(Vop))
	{
		stop("length of y must match nrows of Vop")
	}
	if(ncol(Vp) != ncol(Vop))
	{
		stop("ncols of Vp must match ncols of Vop")
	}

	if(!is.numeric(coeff))
	{
		stop("coeff must be a numeric vector")
	}
	if((is.null(X) && !is.null(Xp)) || (!is.null(X) && is.null(Xp)))
	{
		stop("If X is supplied, Xp must also be supplied (and vice versa)")
	}
	if(!is.null(X))
	{
		if(nrow(X) != n)
		{
			stop("nrows of X must match length of y")
		}
		if(ncol(X) != nk)
		{
			stop("ncols of X must match length of coeff")
		}
		if(nrow(Xp) != nrow(Vp))
		{
			stop("nrows of Xp must match nrows of Vp")
		}
		if(ncol(Xp) != ncol(X))
		{
			stop("ncols of Xp must match ncols of X")
		}
	}
	if(!valid_decomp_type)
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
	}
}




