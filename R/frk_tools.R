check_frk_arg1 <- function(y, K, S, Sp, X, Xp)
{
	if(!is.vector(y))
	{ stop("y must be a vector") } 
	if(!is.numeric(y))
	{ stop("y must be numeric") } 
	n <- length(y)
	
	if(!is.matrix(K))
	{ stop("K must be a matrix") } 
	if(!is.numeric(K))
	{ stop("K must be numeric") } 
	if(nrow(K) != ncol(K))
	{ stop("nrow(K) must equal ncol(K)") }
	r <- nrow(K)

	if(!is.matrix(S))
	{ stop("S must be a matrix") } 
	if(!is.numeric(S))
	{ stop("S must be numeric") } 
	if(nrow(S) != n )
	{ stop("nrow(S) must equal length(y)") }
	if(ncol(S) != r)
	{ stop("ncol(S) must equal nrow(K)") }

	if(!is.matrix(Sp))
	{ stop("Sp must be a matrix") } 
	if(!is.numeric(Sp))
	{ stop("Sp must be numeric") } 
	if(ncol(Sp) != r)
	{ stop("ncol(Sp) must equal nrow(K)") }

	if(!is.matrix(X))
	{ stop("X must be a matrix") } 
	if(!is.numeric(X))
	{ stop("X must be numeric") } 
	if(nrow(X) != n)
	{ stop("nrow(X) must equal length(y)") }

	if(!is.matrix(Xp))
	{ stop("Xp must be a matrix") } 
	if(!is.numeric(Xp))
	{ stop("Xp must be numeric") } 
	if(nrow(Xp) != nrow(Sp))
	{ stop("nrow(Xp) must equal nrow(Sp)") }
}

check_frk_arg2 <- function(error.var, finescale.var, 
	Ve.diag, Vf.diag, Vfp.diag, coords, pcoords, n, np)
{
	if(!is.vector(error.var))
	{ stop("error.var must be a vector") } 
	if(!is.numeric(error.var))
	{ stop("error.var must be numeric") } 
	if(length(error.var) != 1)
	{ stop("error.var must have length 1") }
	if(error.var <= 0)
	{ stop("error.var must be a positive value") }

	if(!is.vector(finescale.var))
	{ stop("finescale.var must be a vector") } 
	if(!is.numeric(finescale.var))
	{ stop("finescale.var must be numeric") } 
	if(length(finescale.var) != 1)
	{ stop("finescale.var must have length 1") }
	if(finescale.var < 0)
	{ stop("finescale.var must be a non-negative value") }

	if(!is.null(Ve.diag))
	{
		if(!is.vector(Ve.diag))
		{ stop("Ve.diag must be a vector") } 
		if(!is.numeric(Ve.diag))
		{ stop("Ve.diag must be numeric") } 
		if(length(Ve.diag) != n)
		{ stop("length(Ve.diag) must equal length(y)") }
	}

	if(!is.null(Vf.diag))
	{
		if(!is.vector(Vf.diag))
		{ stop("Vf.diag must be a vector") } 
		if(!is.numeric(Vf.diag))
		{ stop("Vf.diag must be numeric") } 
		if(length(Vf.diag) != n)
		{ stop("length(Vf.diag) must equal length(y)") }
	}

	if(!is.null(Vfp.diag))
	{
		if(!is.vector(Vfp.diag))
		{ stop("Vfp.diag must be a vector") } 
		if(!is.numeric(Vfp.diag))
		{ stop("Vfp.diag must be numeric") } 
		if(length(Vfp.diag) != np)
		{ stop("length(Vfp.diag) must equal nrow(Sp)") }
	}
	
	if(finescale.var > 0)
	{
		if(is.null(coords))
		{ stop("coords must be provided when finescale.var > 0") }
		if(!is.matrix(coords))
		{ stop("coords must be a matrix") }
		if(!is.numeric(coords))
		{ stop("coords must be numeric") }
		if(nrow(coords) != n)
		{ stop("nrow(coords) must equal length(y)") }
	
		if(is.null(pcoords))
		{ stop("pcoords must be provided when finescale.var > 0") }
		if(!is.matrix(pcoords))
		{ stop("pcoords must be a matrix") }
		if(!is.numeric(pcoords))
		{ stop("pcoords must be numeric") }
		if(nrow(pcoords) != np)
		{ stop("nrow(pcoords) must equal nrow(Sp)") }
		if(ncol(pcoords) != ncol(coords))
		{ stop("ncol(pcoords) must equal ncol(coords)") }

	}
}

check_arg_coincident <- function(coords1, coords2)
{
	if(!is.matrix(coords1))
	{ stop("coords1 must be a matrix") }
	if(!is.numeric(coords1))
	{ stop("coords1 must be numeric") }
	if(ncol(coords1) != 2)
	{ stop("coords1 must have two columns") }

	
	if(!is.matrix(coords2))
	{ stop("coords2 must be a matrix") }
	if(!is.numeric(coords2))
	{ stop("coords2 must be numeric") }
	if(ncol(coords2) != 2)
	{ stop("coords2 must have two columns") }
}

coincident <- function(coords1, coords2)
{
	check_arg_coincident(coords1, coords2)

	coin <- matrix(NA, nrow = nrow(coords1), ncol = 2)

	tol <- .Machine$double.eps
	
	for(i in 1:nrow(coords1))
	{
		for(j in 1:nrow(coords2))
		{
			if(abs(coords1[i, 1] - coords2[j, 1]) < tol)
			{
				if(abs(coords1[i, 2] - coords2[j, 2]) < tol)
				{
					coin[i,] <- c(i, j)
				}
			}
		}
	}
	return(coin[which(!is.na(coin[,1]) == 1), ])
}

src <- '
		#include <math.h>

		NumericMatrix x(coords1), y(coords2);
		NumericVector tol(eps);
		NumericMatrix coin(x.nrow(), 2);

		int count = 0;

		for(int i = 1; i <= x.nrow(); i++)
		{
			for(int j = 1; j <= y.nrow(); j++)
			{
				if(fabs(x(i - 1, 0) - y(j - 1, 0)) < tol[0])
				{
					if(fabs(x(i - 1, 1) - y(j - 1, 1)) < tol[0])
					{
						coin(i - 1, 0) = i;
						coin(i - 1, 1) = j;
						count++;
					}
				}
			}
		}

		NumericMatrix return_coin(count, 2);

		for(int i = 0; i < x.nrow(); i++)
		{
			if(coin(i, 0) > 0)
			{
				return_coin(return_coin.nrow() - count, 0) = coin(i, 0);
				return_coin(return_coin.nrow() - count, 1) = coin(i, 1);
				count--;
			}
		}

	 	return return_coin;
'

fun <- cxxfunction(signature(coords1 = "numeric", coords2 = "numeric", 
	eps = "numeric"), 
	src, plugin="RcppArmadillo")

my1 <- as.matrix(expand.grid(seq(0, 1, len = 25), seq(0, 1, len = 25)))
my2 <- as.matrix(expand.grid(seq(0, 1, len = 101), seq(0, 1, len = 101)))

ti <- Sys.time()
for(i in 1:1)
{
	mycoin <- coincident(my1, my2)
}
Sys.time() - ti

ti <- Sys.time()
for(i in 1:1)
{
	mycoin2 <- fun(my1, my2, eps = .Machine$double.eps)
}
Sys.time() - ti

my1b <- as.matrix(expand.grid(seq(0, 1, len = 25), seq(0, 1, len = 25)))
my2b <- as.matrix(expand.grid(seq(2, 3, len = 101), seq(0, 1, len = 101)))

	mycoin <- coincident(my1b, my2b)
	mycoin2 <- fun(my1b, my2b, eps = .Machine$double.eps)



mycoin2 <- fun(my1, my2, eps = .Machine$double.eps)
range(mycoin - mycoin2)


plot(my2)
points(my1, col = "red")
points(my1[mycoin[,1], ], col = "green")

library(inline)


src <- '
		NumericMatrix x(coords1), y(coords2);
		NumericVector tol(eps);
		//NumericMatrix coin(x.nrow(), 2);
		NumericMatrix coin(x.nrow(), y.nrow());
		NumericMatrix coin2(x.nrow(), y.nrow());


		double d1 = 0;
		double d2 = 0;
		
		for(int i = 1; i <= x.nrow(); i++)
		{
			for(int j = 1; j <= y.nrow(); j++)
			{
				d1 = x(i - 1, 0) - y(j - 1, 0);
		
				if(d1 < tol[0])
				{
				
					d2 = x(i - 1, 1) - y(j - 1, 1);
	
					if(d2 < tol[0])
					{
						//coin(i - 1, 0) = i;
						//coin(i - 1, 1) = j;
						
						coin(i - 1, j - 1) = d1;
						coin2(i - 1, j - 1) = d2;
					
						std::cout << "i = " << i << "j = " << j << " "; 
						//	<< "d1 = ", abs(x(i - 1, 0) - y(j - 1, 0)) << " " 
						//	<< "d2 = ", abs(x(i - 1, 1) - y(j - 1, 1)) << " ";
						//std::cout << "i = " << i << "j = " << j << " " 
						//	<< "d1 = ", d1 << " " 
						//	<< "d2 = ", d2 << " ";
					}
				}
	
			}
		}

				return Rcpp::List::create(Rcpp::Named("coin") = coin,
   								  Rcpp::Named("coin2") = coin2
								  );
'

src <- '
		#include <math.h>
	
		NumericMatrix x(coords1), y(coords2);
		NumericVector tol(eps);
		//NumericMatrix coin(x.nrow(), 2);
		NumericMatrix coin(x.nrow(), y.nrow());
		NumericMatrix coin2(x.nrow(), y.nrow());

		for(int i = 1; i <= x.nrow(); i++)
		{
			for(int j = 1; j <= y.nrow(); j++)
			{
				std::cout << i << " " << j;
			
				if(fabs(x(i - 1, 0) - y(j - 1, 0)) < tol[0])
				{
					std::cout << " diff = " << fabs(x(i - 1, 0) - y(j - 1, 0)) << 
						" true newline ";
				}
				else
				{	
					std::cout << " diff = " << fabs(x(i - 1, 0) - y(j - 1, 0)) <<
						" false newline ";
				}
			}
		}
'


matrix(my1[,1], nrow = 4, ncol = nrow(my2)) - 
	matrix(my2[, 1], nrow = 4, ncol = nrow(my2), byrow = TRUE)

#
#		//int tot = x.nrow() * y.nrow();
#		//NumericMatrix coin(, 1);
#
#		//int nx = x.nrow();
# 		//int z = 0;
# 
#		/*for(int i = 1; i <= x.nrow(); i++)
#		{
#			for(int j = 1; j <= y.nrow(); j++)
#			{
#				//if(abs(x(i - 1, 0) - y(j - 1, 0)) < tol)
#				//{
#				//	if(abs(x(i - 1, 1) - y(j - 1, 1)) < tol)
#				//	{
#						coin(z, 0) = abs(x(i - 1, 1) - y(j - 1, 1));
#				//	}
#				//}
#			}
#		}*/

krige.frk3 <- function(y, K, S, Sp, X, Xp, error.var, finescale.var = 0, 
	Ve.diag = NULL, Vf.diag = NULL, Vfp.diag = NULL, 
	coords = NULL, pcoords = NULL)
{
	#check that beginning arguments are valid
	check_frk_arg1(y, K, S, Sp, X, Xp)

	n <- length(y); r <- nrow(K); np <- nrow(Sp)
	
	#check that second set of arguments are valid
	check_frk_arg2(error.var, finescale.var, Ve.diag, Vf.diag, 
		Vfp.diag, coords, pcoords, n, np)

	#create Ve, Vf, and Vfp.diag vectors if needed
	if(is.null(Ve.diag)){ Ve.diag <- rep(1, n) }
	if(is.null(Vf.diag)){ Vf.diag <- rep(1, n) }
	if(is.null(Vfp.diag)){ Vfp.diag <- rep(1, np) }

	D.vec <- error.var * ve + finescale.var * vf
	Di.vec <- 1/D.vec
	
	###compute matrix products for future use
	#compute Di %*% S (in a faster way)
	DiS <- Di.vec * S

	#Computer inverse of complete covariance matrix C
	Ci <- diag(Di.vec) - DiS %*% solve(solve(K) + crossprod(S, DiS), t(DiS))

	#inverse of C
	XtCiX <- crossprod(X, Ci %*% X)

	#compute gls estimates of regression coefficients
	coeff <- solve(XtCiX, crossprod(X, Ci %*% y))

	#crosscovariance between predicted and observed responses
	#V <- tcrossprod(S %*% K, S) + diag(1/Di.vec)

	if(finescale.var > 0)
	{
		coin <- coincident(coords, pcoords)
	}

	Vop <- tcrossprod(S %*% K, Sp)
	for(i in 1:nrow(coin))
	{
		Vop[coin[i, 1], coin[i, 2]] <- Vop[coin[i, 1], coin[i, 2]] + 
			Vf.diag[coin[i, 1]]*finescale.var
	}

	#Fast way of computing diag.Vp <- diag(Vp) with
	#Vp = Sp %*% K %*% t(Sp) + diag(finescale.var * Vfp.diag)

	diag.Vp <- rowSums(Sp * (Sp %*% K)) + finescale.var * Vfp.diag
	
	#compute kriging weights
	w <- Ci %*% (Vop - X %*% solve(XtCiX, crossprod(X, Ci %*% Vop) - t(Xp)))

	#determine predicted value and mean-square prediction error
	pred <- crossprod(w, y)
	#mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag.Vp
	mspe <- colSums(((tcrossprod(S %*% K, S) + diag(D.vec)) %*% w) * w) - 
		2 * colSums(w * Vop) + diag.Vp

	return(list(pred = pred, mspe = mspe, w = w, coeff = coeff))
}

krige.frk <- function(y, K, S, Sp, V.error, X, Xp)
{
	###compute matrix products for future use
	Vei <- solve(V.error) #inverse of V.error
	VeiS <- solve(V.error, S)
	#inverse of V
	Vi <- Vei - VeiS %*% solve(solve(K) + crossprod(S, VeiS), t(VeiS))
	ViX <- Vi %*% X
	XtViX <- crossprod(X, ViX)
	#crosscovariance between predicted and observed responses
	V <- tcrossprod(S %*% K, S) + V.error
	Vop <- tcrossprod(S %*% K, Sp)
	Vp <- Sp %*% K %*% t(Sp)
	
	#compute gls estimates of regression coefficients
	coeff <- solve(XtViX, crossprod(ViX, y))
	
	#compute kriging weights
	w <- Vi %*% (Vop - X %*% solve(XtViX, crossprod(X, Vi %*% Vop) - t(Xp)))

	pred <- crossprod(w, y)
	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)

	return(list(pred = pred, mspe = mspe, w = w))
}

krige.frk2 <- function(y, K, S, Sp, error.var, X, Xp)
{
	#overall covariance matrix is V = S %*% K %*% t(S) + V.error + V.finescale 
	#V = S %*% K %*% t(S)
	#Vop = S %*% K %*% t(Sp) + diag(V.finescale)
	#Vp = Sp %*% K %*% t(Sp) + diag(V.finescale)

	#make vectors arguments are vectors, convert them to vector
#	if(!is.vector(y)){ y <- as.vector(y) }
#	if(!is.vector(error.var)){ y <- as.vector(error.var) }
#	if(!is.vector(finescale.var)){ y <- as.vector(finescale.var) }
#			
#	check_frk_arg(y, K, X, Xp, S, Sp, error.var, finescale.var)
#
#	n <- length(y)
#
#	D <- diag(n)*(error.var + finescale.var)
	
	###compute matrix products for future use
	#compute t(S) %*% solve(diag(error.var))
	StVei <- t(S) * matrix(1/error.var, nrow = r, ncol = n, byrow = TRUE) 

	#Computer inverse of covariance matrix
	Vi <- diag(1/error.var) - crossprod(StVei, solve(solve(K) + StVei %*% S)) %*% StVei

	#inverse of V
	XtViX <- crossprod(X, Vi %*% X)
	#crosscovariance between predicted and observed responses
	V <- tcrossprod(S %*% K, S) + diag(error.var)
	Vop <- tcrossprod(S %*% K, Sp)
	Vp = Sp %*% K %*% t(Sp)
	
	#compute gls estimates of regression coefficients
	coeff <- solve(XtViX, crossprod(X, Vi %*% y))
	
	#compute kriging weights
	w <- Vi %*% (Vop - X %*% solve(XtViX, crossprod(X, Vi %*% Vop) - t(Xp)))

	pred <- crossprod(w, y)
	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)

	return(list(pred = pred, mspe = mspe, w = w))
}






rmvnorm.fr2 <- function(nsim = 100, mu, K, S, error.var, method = "eigen")
{
    mu <- as.vector(mu)
    n <- nrow(S)
    r <- nrow(K)
    
    decomp.K <- cbind(decomp.cov(K, method = method), 
    	matrix(0, nrow = r, ncol = n - r))
    
    return(
    	mu + 
    	S %*% decomp.K %*% matrix(rnorm(n * nsim), nrow = n, ncol = nsim) #+ 
    	#matrix(rnorm(n * nsim, sd = rep(sqrt(error.var), times = nsim)), ncol = nsim)
    )
}

rmvnorm.fr <- function(nsim = 100, mu, K, S, method = "eigen")
{
    mu <- as.vector(mu)
    n <- nrow(S)
    r <- nrow(K)
    
    decomp.K <- cbind(decomp.cov(K, method = method), 
    	matrix(0, nrow = r, ncol = n - r))
    
    return(mu + S %*% decomp.K %*% matrix(rnorm(n * nsim), nrow = n, ncol = nsim))
}

rcondsim.fr <- function(nsim = 100, krige.frk.obj, y, K, S, Sp, error.var, 
	method = "eigen")
{
	n <- nrow(S)
	r <- nrow(K)
	np <- length(krige.frk.obj$pred)

	#decomp.K <- cbind(decomp.cov(K, method = method), 
    #	matrix(0, nrow = r, ncol = n - r + np))

	V <- tcrossprod(S %*% K, S)
	Vop <- tcrossprod(S %*% K, Sp)
	Vp <- Sp %*% K %*% t(Sp)

	Va <- rbind(cbind(V, Vop), cbind(t(Vop), Vp))

	newsim <- decomp.cov(Va, method = method) %*% 
		matrix(rnorm(nrow(Va) * nsim), ncol = nsim)
	
	#newsim <- rbind(S, Sp) %*% 
	#	(decomp.K %*% matrix(rnorm((n + np) * nsim), ncol = nsim))
	
	newZ <- matrix(y, nrow = n, ncol = nsim) - newsim[1:n,] + 
		matrix(rnorm(n * nsim, sd = sqrt(error.var)), nrow = n, ncol = nsim)

	return(newsim[-(1:n), ] + crossprod(krige.frk.obj$w, newZ))
}

rcondsim <- function(nsim = 100, krige.obj, y, V, Vp, Vop, Ve, method = "eigen")
{
	n <- nrow(V)
	np <- nrow(Vp)

	Vomod <- V - Ve
	
	Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))

	newsim <- decomp.cov(Va, method = method) %*% 
		matrix(rnorm(nrow(Va) * nsim), ncol = nsim)

	newZ <- matrix(y, nrow = n, ncol = nsim) - newsim[1:n,] + 
		matrix(rnorm(n * nsim, sd = sqrt(diag(Ve))), nrow = n, ncol = nsim)

	return(newsim[-(1:n),] + crossprod(krige.obj$w, newZ))
}

spatial.basis <- function(refcoords, pcoords, rad)
{
	D <- dist2(pcoords, refcoords)
	B <- (1 - (D/rad)^2)^2 * (D <= rad)
	return(B)
}

