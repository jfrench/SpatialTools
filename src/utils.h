#ifndef _SpatialTools_UTILS_H
#define _SpatialTools_UTILS_H

#include <RcppArmadillo.h>

arma::mat decomp_V(arma::mat Va, int method);
arma::mat rcondsim(int nsim, arma::vec y, arma::mat w, arma::mat Vediag, arma::mat dV, int method, double m);

#endif
