#include "krige_sk_tools.h"

using namespace Rcpp;

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
                              Rcpp::Named("w") = w
                              );
}