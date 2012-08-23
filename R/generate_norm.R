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
