#include "krige.h"

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

SEXP krige_sk(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP ms){
    
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
    
    //compute kriging weights
    //R version: w <- solve(V, Vop)
    arma::mat w = arma::solve(V, Vop);
    
    //blup for Yp
    //R version:  pred <- m + crossprod(w, y - m)
    arma::colvec pred = m + trans(w) * (y - m);
    
    //variance of (Yp - pred)
    //R version: mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)
    arma::rowvec mspe = sum((V * w) % w, 0) - 2 * sum(w % Vop) + trans(diagvec(Vp));
    
    return Rcpp::List::create(Rcpp::Named("pred") = pred,
                              Rcpp::Named("mspe") = mspe,
                              Rcpp::Named("w") = w,
                              Rcpp::Named("mean") = m
                              );
}

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
    
    return Rcpp::List::create(Rcpp::Named("pred") = pred,
                              Rcpp::Named("mspe") = mspe,
                              Rcpp::Named("w") = w,
                              Rcpp::Named("coeff") = coeff,
                              Rcpp::Named("vcov.coeff") = vcov_coef
                              );
}
