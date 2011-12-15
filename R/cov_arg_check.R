cov_sp_arg_check <- function(coords, sp.type, sp.par, 
	error.var, smoothness, finescale.var, pcoords, 
	D, Dp, Dop)
{
	#check coords argument
	if(!is.numeric(coords) || !is.matrix(coords)){ stop("coords must be a numeric matrix") }

	#check sp.type arguments
	if(!valid_sp_type(sp.type)){ stop("sp.type is not a valid covariance type") }

	#check sp.par argument
	if(!is.numeric(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(!is.vector(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(length(sp.par) != 2){ stop("sp.par must be a numeric vector of length 2") }
	if(!(min(sp.par) > 0)){ stop("sp.par must have positive elements") }


	#check error.var argument
	if(!is.numeric(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(!is.vector(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(length(error.var) > 1){ stop("error.var must be a numeric vector of length 1") }
	if(min(error.var) < 0){ stop("error.var must be non-negative") }


	#check smoothness argument
	if(!is.numeric(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(!is.vector(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(length(smoothness) > 1){ stop("smoothness must be a numeric vector of length 1") }
	if(min(smoothness) <= 0){ stop("smoothness must be positive") }

	#check finescale.var argument
	if(!is.numeric(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(!is.vector(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(length(finescale.var) > 1){ stop("finescale.var must be a numeric vector of length 1") }
	if(min(finescale.var) < 0){ stop("finescale.var must be non-negative") }

	#check pcoords argument
	if(!is.null(pcoords))
	{
		if(!is.numeric(pcoords) || !is.matrix(pcoords)){ stop("pcoords must be a numeric matrix") }
		if(ncol(coords) != ncol(pcoords))
		{ 
			stop("coords and pcoords must have the same number of columns") 
		}
	}
	
	#check D argument
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D) ){ stop("D must be a numeric matrix if provided") } 
		if(min(D) < 0){ stop("D cannot have negative elements") }
		if(!isSymmetric(D)){ stop("D must be a symmetric matrix") } 
		if(nrow(D) != nrow(coords))
		{
			stop("D must be a matrix with nrows and ncols equal to nrows of coords")
		}
	}

	#check Dp argument
	if(!is.null(Dp))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dp is provided") }
		if(!is.numeric(Dp) || !is.matrix(Dp)){ stop("Dp must be a numeric matrix if provided") } 
		if(min(Dp) < 0){ stop("Dp cannot have negative elements") }
		if(!isSymmetric(Dp)){ stop("Dp must be a symmetric matrix") }
		if(nrow(Dp) != nrow(pcoords))
		{
			stop("Dp must be a matrix with nrows and ncols equal to nrows of pcoords")
		}
	}

	#check Dp argument
	if(!is.null(Dop))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dop is provided") }
		if(!is.numeric(Dop)  || !is.matrix(Dop)){ stop("Dop must be a numeric matrix if provided") } 
		if(min(Dop) < 0){ stop("Dop cannot have negative elements") }
		if(ncol(Dop) != nrow(pcoords))
		{
			stop("Dop must have ncols equal to nrows of pcoords")
		}
		if(nrow(Dop) != nrow(coords))
		{
			stop("Dop must have nrows equal to nrows of coords")
		}
	}
}

valid_sp_type <- function(sp.type)
{
	#returns TRUE if sp.type is equal to the below options
	#otherwise it returns FALSE
	return((sp.type == "exponential" ||
		sp.type == "gaussian" ||
		sp.type == "matern" ||
		sp.type == "spherical"))
}

valid_t_type <- function(t.type)
{
	#returns TRUE if sp.type is equal to the below options
	#otherwise it returns FALSE
	return((t.type == "ar1"))
}

valid_decomp_type <- function(method)
{
	return((method == "eigen" || method == "chol" || method == "svd"))
}

cov_st_arg_check <- function(coords, time, sp.type, sp.par, 
	error.var, smoothness, finescale.var, t.type, t.par, pcoords, 
	ptime, D, Dp, Dop, T, Tp, Top)
{
	#check coords argument
	if(!is.numeric(coords) || !is.matrix(coords)){ stop("coords must be a numeric matrix") }

	#check time argument
	if(!is.numeric(time) || nrow(time) != nrow(coords))
	{ 
		stop("time must be a numeric matrix with nrows equal to nrow(coords) (or a vector with length equal to nrow(coords))")
	} 

	#check sp.type arguments
	if(!valid_sp_type(sp.type)){ stop("specified sp.type is not a valid covariance type") }

	#check sp.par argument
	if(!is.numeric(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(!is.vector(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(length(sp.par) != 2){ stop("sp.par must be a numeric vector of length 2") }
	if(!(min(sp.par) > 0)){ stop("sp.par must have positive elements") }


	#check error.var argument
	if(!is.numeric(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(!is.vector(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(length(error.var) > 1){ stop("error.var must be a numeric vector of length 1") }
	if(min(error.var) < 0){ stop("error.var must be non-negative") }


	#check smoothness argument
	if(!is.numeric(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(!is.vector(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(length(smoothness) > 1){ stop("smoothness must be a numeric vector of length 1") }
	if(min(smoothness) <= 0){ stop("smoothness must be positive") }

	#check finescale.var argument
	if(!is.numeric(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(!is.vector(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(length(finescale.var) > 1){ stop("finescale.var must be a numeric vector of length 1") }
	if(min(finescale.var) < 0){ stop("finescale.var must be non-negative") }

	#check t.type argument
	if(!valid_t_type(t.type)){ stop("specified t.type does not match available options") }

	#check t.par argument
	if(!is.numeric(t.par) || t.par < 0 || t.par >= 1){ stop("t.par must be in range [0, 1)") }
	
	#check pcoords argument
	if(!is.null(pcoords))
	{
		if(!is.numeric(pcoords)){ stop("pcoords must be a numeric matrix") }
		if(ncol(coords) != ncol(pcoords))
		{ 
			stop("coords and pcoords must have the same number of columns") 
		}
		if(is.null(ptime)){ stop("ptime must be supplied when pcoords is supplied.") }
		
		#check ptime argument
		if(!is.numeric(ptime) || nrow(ptime) != nrow(pcoords))
		{ 
			stop("ptime must be a numeric matrix with nrows equal to nrow(pcoords) (or a vector with length equal to nrow(pcoords)")
		} 
	}

	
	#check D argument
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D) ){ stop("D must be a numeric matrix if provided") } 
		if(min(D) < 0){ stop("D cannot have negative elements") }
		if(!isSymmetric(D)){ stop("D must be a symmetric matrix") } 
		if(nrow(D) != nrow(coords))
		{
			stop("D must be a matrix with nrows and ncols equal to nrows of coords")
		}
	}

	#check Dp argument
	if(!is.null(Dp))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dp is provided") }
		if(!is.numeric(Dp) || !is.matrix(Dp)){ stop("Dp must be a numeric matrix if provided") } 
		if(min(Dp) < 0){ stop("Dp cannot have negative elements") }
		if(!isSymmetric(Dp)){ stop("Dp must be a symmetric matrix") }
		if(nrow(Dp) != nrow(pcoords))
		{
			stop("Dp must be a matrix with nrows and ncols equal to nrows of pcoords")
		}
	}

	#check Dp argument
	if(!is.null(Dop))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dop is provided") }
		if(!is.numeric(Dop)  || !is.matrix(Dop)){ stop("Dop must be a numeric matrix if provided") } 
		if(min(Dop) < 0){ stop("Dop cannot have negative elements") }
		if(ncol(Dop) != nrow(pcoords))
		{
			stop("Dop must have ncols equal to nrows of pcoords")
		}
		if(nrow(Dop) != nrow(coords))
		{
			stop("Dop must have nrows equal to nrows of coords")
		}
	}

	#check T argument
	if(!is.null(T))
	{
		if(!is.numeric(T) || !is.matrix(T) ){ stop("T must be a numeric matrix if provided") } 
		if(min(T) < 0){ stop("T cannot have negative elements") }
		if(!isSymmetric(T)){ stop("T must be a symmetric matrix") } 
		if(nrow(T) != nrow(coords))
		{
			stop("T must be a matrix with nrows and ncols equal to nrows of coords")
		}
	}

	#check Tp argument
	if(!is.null(Tp))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Tp is provided") }
		if(!is.numeric(Tp) || !is.matrix(Tp)){ stop("Tp must be a numeric matrix if provided") } 
		if(min(Tp) < 0){ stop("Tp cannot have negative elements") }
		if(!isSymmetric(Tp)){ stop("Tp must be a symmetric matrix") }
		if(nrow(Tp) != nrow(pcoords))
		{
			stop("Tp must be a matrix with nrows and ncols equal to nrows of pcoords")
		}
	}

	#check Top argument
	if(!is.null(Top))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Top is provided") }
		if(!is.numeric(Top)  || !is.matrix(Top)){ stop("Top must be a numeric matrix if provided") } 
		if(min(Top) < 0){ stop("Top cannot have negative elements") }
		if(ncol(Top) != nrow(pcoords))
		{
			stop("Top must have ncols equal to nrows of pcoords")
		}
		if(nrow(Top) != nrow(coords))
		{
			stop("Top must have nrows equal to nrows of coords")
		}
	}
}

decomp_cov_check_arg <- function(V, method, checkSymmetric = TRUE)
{
	if(!is.matrix(V) || !is.numeric(V))
	{
		stop("V must be a numeric matrix")
	}
	#Removed because sometimes a symmetric matrix may not be due to numerical imprecision
	#if(checkSymmetric)
	#{
	#	if(!isSymmetric(V)){ stop("V must be symmetric") }
	#}
	if(!(method == "eigen" || method == "chol" || method == "svd"))
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
	}
}

maxlik_cov_sp_check_arg <- function(X, y, coords, sp.type, 
	range.par, error.ratio, smoothness, D, reml, lower, upper)
{
	if(!is.numeric(X) || !is.matrix(X)){ stop("X must be a numeric matrix") }
	if(!is.numeric(y)){ stop("y must be numeric") }
	if(!is.numeric(coords) || !is.numeric(coords)){ stop("coords must be a numeric matrix") }
	if(!valid_sp_type(sp.type)){ stop("specified sp.type is not a valid covariance type") }
	if(!(range.par > 0)){ stop("range.par must be positive") }
	if(!(error.ratio >= 0)){ stop("range.par must be non-negative") }
	if(!(smoothness > 0)){ stop("smoothness must be positive") }
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D))
		{ 
			stop("If supplied, D must be a numeric matrix") 
		}
	}
	if(!is.logical(reml)){ stop("reml must be a logical value") }
	if(!is.null(lower))
	{
		if(!is.numeric(lower)){ stop("lower must be a numeric vector")} 
	}
	if(!is.null(upper))
	{
		if(!is.numeric(upper)){ stop("upper must be a numeric vector")} 
	}
	if(!is.null(lower) && !is.null(upper))
	{
		if(length(lower) != length(upper)){ stop("lower and upper should have the same length") }
	}
}
