#include "krige_ok_tools.h"

using namespace Rcpp;

SEXP krige_ok(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops){

		NumericVector yr(ys);
		arma::colvec y(yr.begin(), yr.size(), false);
		
		NumericMatrix Vr(Vs);
		arma::mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);
		
		NumericMatrix Vpr(Vps);
		arma::mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);
		
		NumericMatrix Vopr(Vops);
		arma::mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);
		
		arma::colvec X = arma::ones(yr.size());
		
		//compute useful matrices.  Some expressions simply because X is a vector of 1s
		//sum(A, 0) is the equivalent of colSums(A) where A is a matrix
		arma::mat ViX = arma::solve(V, X);
		double XtViX = as_scalar(sum(ViX, 0));
		double vcov_coef = 1/XtViX;

		//compute gls estimates of regression coefficients
		//equivalent to solve(XtViX, crossprod(ViX, y))
		double coeff = as_scalar(sum(ViX % y)/XtViX);
		
		//compute kriging weights
		//R version: w <- solve(V, Vop - tcrossprod(X, (colSums(solve(V, Vop)) - 1)/XtViX))
		arma::mat w = arma::solve(V, Vop - repmat( (sum(arma::solve(V, Vop), 0) - 1)/XtViX, y.n_elem, 1));
		
		//blup for Yp
		//R version:  pred <- crossprod(w, y)
		arma::colvec pred = trans(w) * y;
		
		//variance of (Yp - pred)
		//R version: mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)
		arma::rowvec mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));
		
		return Rcpp::List::create(Rcpp::Named("pred") = pred,
								  Rcpp::Named("mspe") = trans(mspe),
								  Rcpp::Named("w") = w,
								  Rcpp::Named("coeff") = coeff,
								  Rcpp::Named("vcov.coeff") = vcov_coef
								  );
}