#ifndef _SpatialTools_KRIGE_UK_TOOLS_H
#define _SpatialTools_KRIGE_UK_TOOLS_H

#include <RcppArmadillo.h>

RcppExport SEXP krige_uk(SEXP Xs, SEXP ys, SEXP Vs, SEXP Xps, SEXP Vps, SEXP Vops);

RcppExport SEXP pweights_uk(SEXP Xs, SEXP Vs, SEXP Xps, SEXP Vps, SEXP Vops);

RcppExport SEXP mspe_uk(SEXP ws, SEXP Vs, SEXP Vps, SEXP Vops);

#endif
