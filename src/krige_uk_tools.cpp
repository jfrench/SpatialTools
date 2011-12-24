#include "krige_uk_tools.h"

using namespace Rcpp;
using namespace arma;

SEXP krige_uk(SEXP Xs, SEXP ys, SEXP Vs, SEXP Xps, SEXP Vps, SEXP Vops){

	NumericMatrix Xr(Xs);
	mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

	NumericVector yr(ys);
	colvec y(yr.begin(), yr.size(), false);

	NumericMatrix Vr(Vs);
	mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);

	NumericMatrix Xpr(Xps);
	mat Xp(Xpr.begin(), Xpr.nrow(), Xpr.ncol(), false);

	NumericMatrix Vpr(Vps);
	mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);

	NumericMatrix Vopr(Vops);
	mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);

	//compute useful matrices
	mat ViX = solve(V, X);
	mat XtViX = trans(X) * ViX;
	mat vcov_coef = inv(XtViX);

	//compute generalized least squares estimates of regression coefficients
	colvec coeff = solve(sympd(XtViX), trans(ViX) * y);
	
	//compute kriging weights
	mat w = solve(V, Vop - X * solve(sympd(XtViX), trans(X) * solve(sympd(V), Vop) - trans(Xp)));
	
	//best linear unbiased predictor of response at prediction locations
	mat pred = trans(w) * y;
	
	//calculate mean-square prediciton error of predicted responses.
	//sum(A, dim = 0) is the equivalent of colSums(A) where A is a matrix
	//diagvec(Vp) is equivalent to diag(Vp)
	mat mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));
	
	return Rcpp::List::create(Rcpp::Named("pred")= pred,  							  			  Rcpp::Named("mspe")= mspe,  							  			  Rcpp::Named("w")= w,  
							  Rcpp::Named("coeff")= coeff,
							  Rcpp::Named("vcov.coeff")= vcov_coef
							  );
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
