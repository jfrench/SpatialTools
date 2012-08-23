#include "krige.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;

SEXP krige_uk(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP Xs, SEXP Xps, SEXP nsims, 
	SEXP Vediags, SEXP methods){
	
	//Set up data frame
	NumericVector yr(ys);
	arma::colvec y(yr.begin(), yr.size(), false);

	NumericMatrix Vr(Vs);
	arma::mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);

	NumericMatrix Vpr(Vps);
	arma::mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);

	NumericMatrix Vopr(Vops);
	arma::mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);
	
	NumericMatrix Xr(Xs);
	arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

	NumericMatrix Xpr(Xps);
	arma::mat Xp(Xpr.begin(), Xpr.nrow(), Xpr.ncol(), false);

	int nsim = as<int>(nsims);

	NumericVector Vediagr(Vediags);
	arma::colvec Vediag(Vediagr.begin(), Vediagr.size(), false);

	int method = as<int>(methods);

	//compute useful matrices
	arma::mat ViX = solve(V, X);
	arma::mat XtViX = trans(X) * ViX;
	arma::mat vcov_coef = inv(XtViX);
	
	//compute generalized least squares estimates of regression coefficients
	arma::colvec coeff = solve(sympd(XtViX), trans(ViX) * y);
	
	//compute kriging weights
	arma::mat w = solve(V, Vop - X * solve(sympd(XtViX), trans(X) * solve(sympd(V), Vop) - trans(Xp)));
	
	//best linear unbiased predictor of response at prediction locations
	arma::mat pred = trans(w) * y;
	
	//calculate mean-square prediciton error of predicted responses.
	//sum(A, dim = 0) is the equivalent of colSums(A) where A is a matrix
	//diagvec(Vp) is equivalent to diag(Vp)
	arma::mat mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));

	//If nsim = 0, then don't do conditional simulation
	if(nsim == 0)
	{
		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								Rcpp::Named("mspe") = mspe,
								Rcpp::Named("coeff") = coeff,
								Rcpp::Named("vcov.coeff") = vcov_coef);		
	}
	else // Do conditional simulation
	{
		// Modify V matrix to not include measurement error
		arma::mat Vomod = V - diagmat(Vediag);
		
		// Create combined observed, predicted covariance matrix
		// Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))
		arma::mat Va = join_cols(join_rows(Vomod, Vop), join_rows(trans(Vop), Vp));
	
		int na = Va.n_rows;
	
		arma::mat dV = arma::mat(na, na);
		
		dV = decomp_V(Va, method);

		arma::mat simulations = rcondsim(nsim, y, w, Vediag, dV, method, 0);

		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								Rcpp::Named("mspe") = mspe,
								Rcpp::Named("coeff") = coeff,
								Rcpp::Named("vcov.coeff") = vcov_coef,
								Rcpp::Named("w") = w,
								Rcpp::Named("simulations") = simulations);		
	}
}

SEXP krige_sk(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP ms, SEXP nsims, 
	SEXP Vediags, SEXP methods){
    
    NumericVector yr(ys);
    arma::colvec y(yr.begin(), yr.size(), false);
    
    NumericMatrix Vr(Vs);
    arma::mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);
    
    NumericMatrix Vpr(Vps);
    arma::mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);
    
    NumericMatrix Vopr(Vops);
    arma::mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);
    
    NumericVector mr(ms);
    double m = mr[0]; 
    
	int nsim = as<int>(nsims);

	NumericVector Vediagr(Vediags);
	arma::colvec Vediag(Vediagr.begin(), Vediagr.size(), false);

	int method = as<int>(methods);

    //compute kriging weights
    //R version: w <- solve(V, Vop)
    arma::mat w = arma::solve(V, Vop);
    
    //blup for Yp
    //R version:  pred <- m + crossprod(w, y - m)
    arma::colvec pred = m + trans(w) * (y - m);
    
    //variance of (Yp - pred)
    //R version: mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)
    arma::rowvec mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));
    
	//If nsim = 0, then don't do conditional simulation
	if(nsim == 0)
	{
		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								Rcpp::Named("mspe") = mspe,
                        		Rcpp::Named("mean") = m);		
	}
	else // Do conditional simulation
	{
		// Modify V matrix to not include measurement error
		arma::mat Vomod = V - diagmat(Vediag);
		
		// Create combined observed, predicted covariance matrix
		// Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))
		arma::mat Va = join_cols(join_rows(Vomod, Vop), join_rows(trans(Vop), Vp));
	
		int na = Va.n_rows;
	
		arma::mat dV = arma::mat(na, na);
		
		dV = decomp_V(Va, method);

		arma::mat simulations = rcondsim(nsim, y, w, Vediag, dV, method, m);

		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								Rcpp::Named("mspe") = mspe,
								Rcpp::Named("simulations") = simulations,
                                Rcpp::Named("mean") = m
   		                        );		
	}
}

SEXP krige_ok(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP nsims, 
	SEXP Vediags, SEXP methods){
	
	//Set up data frame
	NumericVector yr(ys);
	arma::colvec y(yr.begin(), yr.size(), false);

	NumericMatrix Vr(Vs);
	arma::mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);

	NumericMatrix Vpr(Vps);
	arma::mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);

	NumericMatrix Vopr(Vops);
	arma::mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);
	
	int nsim = as<int>(nsims);

	NumericVector Vediagr(Vediags);
	arma::colvec Vediag(Vediagr.begin(), Vediagr.size(), false);

	int method = as<int>(methods);

  
    arma::colvec X = arma::ones(yr.size());
    
    //compute useful matrices.  Some expressions simply because X is a vector of 1s
    //sum(A, 0) is the equivalent of colSums(A) where A is a matrix
    arma::mat ViX = arma::solve(V, X);
    double XtViX = arma::as_scalar(sum(ViX, 0));
    double vcov_coef = 1/XtViX;
    
    //compute gls estimates of regression coefficients
    //equivalent to solve(XtViX, crossprod(ViX, y))
    double coeff = arma::as_scalar(sum(ViX % y)/XtViX);
    
    //compute kriging weights
    //R version: w <- solve(V, Vop - tcrossprod(X, (colSums(solve(V, Vop)) - 1)/XtViX))
    arma::mat w = arma::solve(V, Vop - repmat( (sum(arma::solve(V, Vop), 0) - 1)/XtViX, y.n_elem, 1));
    
    //blup for Yp
    //R version:  pred <- crossprod(w, y)
    arma::colvec pred = trans(w) * y;
    
    //variance of (Yp - pred)
    //R version: mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)
    arma::rowvec mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));
    
	//If nsim = 0, then don't do conditional simulation
	if(nsim == 0)
	{
		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								Rcpp::Named("mspe") = mspe,
								Rcpp::Named("coeff") = coeff,
								Rcpp::Named("vcov.coeff") = vcov_coef);		
	}
	else // Do conditional simulation
	{
		// Modify V matrix to not include measurement error
		arma::mat Vomod = V - diagmat(Vediag);
		
		// Create combined observed, predicted covariance matrix
		// Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))
		arma::mat Va = join_cols(join_rows(Vomod, Vop), join_rows(trans(Vop), Vp));
	
		int na = Va.n_rows;
	
		arma::mat dV = arma::mat(na, na);
		
		dV = decomp_V(Va, method);

		arma::mat simulations = rcondsim(nsim, y, w, Vediag, dV, method, 0);

		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								Rcpp::Named("mspe") = mspe,
								Rcpp::Named("coeff") = coeff,
								Rcpp::Named("vcov.coeff") = vcov_coef,
								Rcpp::Named("simulations") = simulations);		
	}
}


SEXP pweights_uk(SEXP Xs, SEXP Vs, SEXP Xps, SEXP Vps, SEXP Vops){

	NumericMatrix Xr(Xs);
	mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

	NumericMatrix Vr(Vs);
	mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);

	NumericMatrix Xpr(Xps);
	mat Xp(Xpr.begin(), Xpr.nrow(), Xpr.ncol(), false);

	NumericMatrix Vpr(Vps);
	mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);

	NumericMatrix Vopr(Vops);
	mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);

	//compute useful matrices
	mat XtViX = trans(X) * solve(V, X);;

	//compute kriging weights
	mat w = solve(V, Vop - X * solve(XtViX, trans(X) * solve(V, Vop) - trans(Xp)));

	return Rcpp::wrap(w);
}

SEXP mspe_uk(SEXP ws, SEXP Vs, SEXP Vps, SEXP Vops){

	NumericMatrix wr(ws);
	mat w(wr.begin(), wr.nrow(), wr.ncol(), false);

	NumericMatrix Vr(Vs);
	mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);

	NumericMatrix Vpr(Vps);
	mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);

	NumericMatrix Vopr(Vops);
	mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);

	//calculate mean-square prediciton error of predicted responses.
	//sum(A, dim = 0) is the equivalent of colSums(A) where A is a matrix
	//diagvec(Vp) is equivalent to diag(Vp)
	mat mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));
	
	return Rcpp::wrap(mspe);
}


