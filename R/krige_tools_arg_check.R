krige_arg_check <- function(y, V, Vp, Vop, X, Xp, coeff)
{
	n <- length(y)

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
		stop("length(y) must equal nrow(V)")
	}
	if(length(y) != nrow(Vop))
	{
		stop("length(y) must equal nrow(Vop)")
	}
	if(ncol(Vp) != ncol(Vop))
	{
		stop("ncol(Vp) must equal ncol(Vop)")
	}
	if((is.null(X) && !is.null(Xp)) || (!is.null(X) && is.null(Xp)))
	{
		stop("If X is supplied, Xp must also be supplied (and vice versa)")
	}
	if(!is.null(X))
	{
		if(nrow(X) != n)
		{
			stop("nrow(X) must equal length(y)")
		}
		if(nrow(Xp) != nrow(Vp))
		{
			stop("nrow(Xp) must equal nrow(Vp)")
		}
		if(ncol(Xp) != ncol(X))
		{
			stop("ncol(Xp) must equal ncol(X)")
		}
	}
	if(!is.null(coeff))
	{
		if(!is.numeric(coeff))
		{
			stop("coeff must be a numeric vector")
		}
		if(length(coeff) > 1)
		{
			if(length(coeff) != ncol(X))
			{
				stop("length(coeff) must equal ncol(X)")
			}
		}
	}
}

krige_sk_arg_check <- function(y, V, Vp, Vop, m)
{
	n <- length(y)

	if(!is.numeric(y))
	{
		stop("y must be a numeric vector")
	}
	if(!is.matrix(V) || !is.numeric(V) || (nrow(V)!= ncol(V)))
	{
		stop("V must be a square numeric matrix")
	}
	if(length(y) != nrow(V))
	{
		stop("length(y) must equal nrow(V)")
	}
	if(!is.matrix(Vp) || !is.numeric(Vp) || (nrow(Vp)!= ncol(Vp)))
	{
		stop("Vp must be a square numeric matrix")
	}
	if(!is.matrix(Vop) || !is.numeric(Vop))
	{
		stop("Vop must be a numeric matrix")
	}
	if(length(y) != nrow(Vop))
	{
		stop("length(y) must equal nrow(Vop)")
	}
	if(ncol(Vp) != ncol(Vop))
	{
		stop("ncol(Vp) must equal ncol(Vop)")
	}
	if(!is.numeric(m))
	{
		stop("m must be a numeric vector")
	}
	if(length(m)!=1)
	{
		stop("m must have length 1")
	}
}

pweights_uk_arg_check <- function(X, V, Xp, Vp, Vop)
{
	if(!is.matrix(X))
	{
		stop("X must be a matrix object")
	}
	if(!is.matrix(V))
	{
		stop("V must be a matrix object")
	}
	if(!is.matrix(Xp))
	{
		stop("Xp must be a matrix object")
	}
	if(!is.matrix(Vp))
	{
		stop("Vp must be a matrix object")
	}
	if(!is.matrix(Vop))
	{
		stop("Vop must be a matrix object")
	}
	if(nrow(X) != nrow(V))
	{
		stop("The nrows in X must match nrows in V")
	}
	if(nrow(X) != nrow(Vop))
	{
		stop("The nrows in X must match nrows in Vop")
	}
	if(nrow(Xp) != nrow(Vp))
	{
		stop("The nrows in Xp must match nrows in Vp")
	}
	if(nrow(Xp) != ncol(Vop))
	{
		stop("The nrows in Xp must match ncols in Vop")
	}
	if(ncol(X) != ncol(Xp))
	{
		stop("The ncols in X must match ncols in Xp")
	}
}

mspe_uk_arg_check <- function(w, V, Vp, Vop)
{
	if(!is.matrix(w))
	{
		stop("w must be a matrix object")
	}
	if(!is.matrix(V))
	{
		stop("V must be a matrix object")
	}
	if(!is.matrix(Vop))
	{
		stop("Vop must be a matrix object")
	}
	if(nrow(w) != nrow(V))
	{
		stop("The nrows in w must match nrows in V")
	}
	if(ncol(w) != ncol(Vop))
	{
		stop("The ncols in w must match ncols in Vop")
	}
	if(nrow(w) != nrow(Vop))
	{
		stop("The nrows in w must match nrows in Vop")
	}
	if(ncol(w) != ncol(Vop))
	{
		stop("The ncols in w must match ncols in Vop")
	}

}