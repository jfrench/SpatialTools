#include "utils.h"

using namespace Rcpp;
using namespace arma;

arma::mat decomp_V(arma::mat Va, int method)
{
	int na = Va.n_rows;
    
	arma::mat dV = arma::mat(na, na);
	
	if(method == 1)
    {
        arma::vec eigval = arma::vec(na);
        arma::mat eigvec = arma::mat(na, na);
     
        //compute eigen values and vectors of V
        eig_sym(eigval, eigvec, Va);
     
        for(int i = 0; i < eigval.n_rows; i++)
        {
            if(eigval(i) < 0)
            {
                eigval(i) = 0;		
            }
        }
     
        dV = eigvec * diagmat(sqrt(eigval));
     }
     else if(method == 2)
     {
         dV = trans(arma::chol(Va));
     }
     else
     {
         arma::mat U = arma::mat(na, na);
         arma::mat U2 = arma::mat(na, na);
         arma::vec sv = arma::vec(na);
     
         svd(U, sv, U2, Va);
     
         dV = U * diagmat(sqrt(sv)) * trans(U2);
     }

    return dV;    
}

arma::mat rcondsim(int nsim, arma::vec y, arma::mat w, arma::mat Vediag, arma::mat dV, int method, double m)
{
    int n = y.n_elem;
	int na = dV.n_rows;
	arma::mat condsim = arma::mat(na, nsim);

	RNGScope scope;
	
	NumericVector Zr = NumericVector(rnorm((na + n) * nsim, 0, 1));
	arma::mat Z(Zr.begin(), (na + n), nsim);
    
	arma::mat newsim = dV * Z.rows(0, na - 1);
  
    // Conditional simulation algorithm changes slightly for simple kriging
    // when m != 0.  
    if(m != 0)
    {    
        arma::mat newZ = repmat(y - m, 1, nsim) - newsim.rows(0, n - 1) + diagmat(sqrt(Vediag)) * Z.rows(na, na + n - 1);
    
        //accounts for the fact that indexing starts at zero
        condsim = m + newsim.rows(n, na - 1) + trans(w) * newZ;
    }
    else
    {    
        arma::mat newZ = repmat(y, 1, nsim) - newsim.rows(0, n - 1) + diagmat(sqrt(Vediag)) * Z.rows(na, na + n - 1);
        
        //accounts for the fact that indexing starts at zero
        condsim = newsim.rows(n, na - 1) + trans(w) * newZ;
    }

    return condsim;
}