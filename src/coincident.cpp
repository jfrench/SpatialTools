#include "coincident_cpp.h"
#include <cmath>

using namespace Rcpp;

SEXP coincident_cpp(SEXP coords1, SEXP coords2, SEXP eps){

		NumericMatrix x(coords1), y(coords2);
		NumericVector tol(eps);
		unsigned int nr = x.rows();
		NumericMatrix coin(nr, 2);

		int count = 0;

		for(unsigned int i = 1; i <= nr; i++)
		{
			for(unsigned int j = 1; j <= nr; j++)
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

		for(unsigned i = 0; i < nr; i++)
		{
			if(coin(i, 0) > 0)
			{
				return_coin(return_coin.nrow() - count, 0) = coin(i, 0);
				return_coin(return_coin.nrow() - count, 1) = coin(i, 1);
				count--;
			}
		}

	 	return return_coin;
}