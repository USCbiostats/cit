#include <R.h>
#include <Rmath.h>
#include <vector>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <iostream>
#include <random>       // std::default_random_engine
#include "linearfunc.h"
#include <Rcpp.h>
#include "maxElementWithNan.h"

using namespace Rcpp;
using namespace std;

/*
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: matrix of 0/1 variables
Programmer: Joshua Millstein
*/


// [[Rcpp::export]]
void citbinm_linear( Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, int &maxit, int &nrow, int &dfz, int &dfx,
	Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval4 , Rcpp::NumericVector pval3nc, int &rseed)
{
	unsigned seed = rseed;
	int rw, cl, xl, i, rind, df, nobs, ip, npos, nperm, nmiss, stride;
	double rss5, pv, pvp, maxp, tmp, rhs, testval;
	bool aa, bb, cc, dd, converged;
	const int firstloop = 1000;
	const int posno = 20;
	const double alpha = .1;
	vector<double> pvec;
	vector<vector<double> > LL;
	vector<vector<double> > GG;
	vector<vector<double> > gpred;
	vector<vector<double> > gresid;
	vector<vector<double> > Gp;
	
	gsl_matrix *Lm, *Gm, *cov, *D;
	gsl_vector *Tm, *c, *Gmj;

	double *designmat = new double[ nrow * (dfz + dfx + 1) ];
	double *phenovec = new double[ nrow ];

	LL.resize( nrow );	
	GetRNGstate();
	
	for(rw = 0; rw < nrow; rw++) {
		LL[rw].resize( dfz );	
	}

	for(cl = 0; cl < dfz; cl++) {
		for(rw = 0; rw < nrow; rw++) {
			LL[rw][cl] = L[rw + nrow * cl];
		}
	}
	
	GG.resize( nrow );	
	
	for(rw = 0; rw < nrow; rw++) {
	  GG[rw].resize( dfx );	
	}
	
	for(cl = 0; cl < dfx; cl++) {
	  for(rw = 0; rw < nrow; rw++) {
	    GG[rw][cl] = G[rw + nrow * cl];
	  }
	}
	
// create analysis vectors w/no missing data
		nobs = 0;
		for(rw = 0; rw < nrow; rw++) {
		  nmiss = 0;
		  for(cl = 0; cl < dfz; cl++) {
		    if( LL[rw][cl] == -9999 ) {
				  nmiss++;
			  }                                                
		  }
		  for(cl = 0; cl < dfx; cl++) {
		    if( GG[rw][cl] == -9999 ) {
		      nmiss++;
		    }                                                
		  }
			aa = nmiss == 0;
			bb = T[rw] != -9999;
			if(aa && bb) {
				nobs++;
			}
		}             

		gpred.resize( nobs );	
		gresid.resize( nobs );
		Gp.resize( nobs );
		for(rw = 0; rw < nobs; rw++) {
		  gpred[rw].resize( dfx );	
		  gresid[rw].resize( dfx );
		  Gp[rw].resize( dfx );
		}
		
		
		Lm = gsl_matrix_alloc (nobs, dfz);
		Gm = gsl_matrix_alloc (nobs, dfx);
		Tm = gsl_vector_alloc (nobs);
		Gmj = gsl_vector_alloc (nobs);
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
		  nmiss = 0;
		  for(cl = 0; cl < dfz; cl++) {
		    if( LL[rw][cl] == -9999 ) {
				  nmiss++;
			  }                                                
		   }
		  for(cl = 0; cl < dfx; cl++) {
		    if( GG[rw][cl] == -9999 ) {
		      nmiss++;
		    }                                                
		  }
		  aa = nmiss == 0;
			bb = T[rw] != -9999;			
			
			if(aa && bb) {
				for(cl = 0; cl < dfz; cl++) {
          gsl_matrix_set(Lm, rind, cl, LL[rw][cl]);
		    }
				for(xl = 0; xl < dfx; xl++) {
				  gsl_matrix_set(Gm, rind, xl, GG[rw][xl]);
				}
				gsl_vector_set(Tm, rind, T[rw]);
				rind++;
			}
		}  	
		// fit model T ~ L
		ip = 1 + dfz;                               // intercept + multiple L variable
		for(rw = 0; rw < nobs; rw++) {
		   phenovec[ rw ] = gsl_vector_get(Tm, rw );
		   designmat[ rw * ip  ] = 1;      // intercept
			for(cl = 0; cl < dfz; cl++) {
                  designmat[ rw * ip + 1 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		     }
		}
		df = dfz;
		converged = linearRegCompare( pv, phenovec, designmat, nobs, ip, df );
		if(!converged)Rcpp::Rcout<< "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pv = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();
		pvec.push_back( pv );  // pval for T ~ L, 9 if it did not converge, p1
		
		// fit model T ~ L + G
		stride = ip + dfx;
		for(rw = 0; rw < nobs; rw++) {
		  designmat[ rw * stride ] = 1;      // intercept
			for(cl = 0; cl < dfz; cl++) {
          	designmat[ rw * stride + 1 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   }
			for(cl = 0; cl < dfx; cl++) {
			  designmat[ rw * stride + 1 + dfz + cl  ]  = gsl_matrix_get (Gm, rw, cl);
			}
		}
		
		df = dfx;
		converged = linearRegCompare( pv, phenovec, designmat, nobs, stride, df );
		if(!converged)Rcpp::Rcout<< "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pv = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();
		pvec.push_back( pv );  // pval for T ~ G|L, 9 if it did not converge, p2

		// fit model T ~ G + L to test L 
		stride = ip + dfx;
		for(rw = 0; rw < nobs; rw++) {
		  designmat[ rw * stride  ] = 1;      // intercept
		  for(cl = 0; cl < dfx; cl++) {
		    designmat[ rw * stride + 1 + cl  ]  = gsl_matrix_get (Gm, rw, cl);
		  }
			for(cl = 0; cl < dfz; cl++) {
        designmat[ rw * stride + 1 + dfx + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		  }
		}
		df = dfz;
		converged = linearRegCompare( pv, phenovec, designmat, nobs, stride, df );
		if(!converged)Rcpp::Rcout<< "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pv = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G
		pval3nc[0] = pv;  // pvalue to be used for non-centrality parameter
		// fit model G ~ L
		D = gsl_matrix_alloc (nobs, ip );
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(D, rw, 0, 1.0);      // intercept
			for(cl = 0; cl < dfz; cl++) {
        gsl_matrix_set(D, rw, cl + 1, gsl_matrix_get (Lm, rw, cl));
		  }
		}
		
		// get linear predictors and residuals
		for(xl = 0; xl < dfx; xl++){
  		c = gsl_vector_alloc (ip);
  		cov = gsl_matrix_alloc (ip, ip);
  		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (nobs, ip);
  		for(rw = 0; rw < nobs; rw++) {
  		  gsl_vector_set(Gmj, rw, gsl_matrix_get (Gm, rw, xl));
  		}
  		gsl_multifit_linear (D, Gmj, c, cov, &rss5, work);
  		gsl_multifit_linear_free (work);
  		gsl_matrix_free (cov);
  
  		// residuals for G ~ L
  		for(rw = 0; rw < nobs; rw++) {
  			rhs = 0;
  			for(cl = 0; cl < ip; cl++) {
          rhs += gsl_vector_get (c, cl) * gsl_matrix_get (D, rw, cl);
  		  }
  			gpred[rw][xl] = rhs;
  			tmp = gsl_vector_get (Gmj, rw) - rhs;
  			gresid[rw][xl] = tmp;
  		}
  		gsl_vector_free (c);
		} // end xl loop
		
		// Conduct an initial set of permutations
		npos = 0;
		for(i = 0; i < firstloop; i++){
			// randomly permute residuals
            
            shuffle( gresid.begin(), gresid.end(), std::default_random_engine(seed) );

			// compute G* based on marginal L effects and permuted residuals
			for(rw = 0; rw < nobs; rw++) {
			  for(xl = 0; xl < dfx; xl++) {
				  Gp[rw][xl] = gpred[rw][xl] + gresid[rw][xl];
			  }
			}
			
			// Recompute p-value for T ~ L|G based on G*
			// fit model T ~ G + L to test L 
  		stride = ip + dfx;
  		for(rw = 0; rw < nobs; rw++) {
  		  designmat[ rw * stride  ] = 1;      // intercept
  		  for(xl = 0; xl < dfx; xl++) {
  		    designmat[ rw * stride + 1 + xl  ]  = Gp[rw][xl];
  		  }
  			for(cl = 0; cl < dfz; cl++) {
          designmat[ rw * stride + 1 + dfx + cl  ]  = gsl_matrix_get (Lm, rw, cl);
  		  }
  		}
		
			df = dfz;
			converged = linearRegCompare( pvp, phenovec, designmat, nobs, stride, df );
			if(!converged)Rcpp::Rcout<< "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
			pvp = ( converged ) ? pvp : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G*
			if( pvp > pv ) npos++;
			
		} // end initial permutation loop

		// Conduct additional permutations if there is some indication of statistical significance
		maxp = maxElementWithNan(pvec);
		nperm = firstloop;
		aa = npos < posno;
		bb = maxp < alpha;
		cc = nperm < maxit;
		testval = (double) (npos + 1) / nperm ;
		dd = maxp < testval; // check that other component p-values are small

		if(aa && bb && cc && dd){
			while(aa && cc) {

				// randomly permute residuals
				
				shuffle( gresid.begin(), gresid.end(), std::default_random_engine(seed) );
				
				// compute G* based on marginal L effects and permuted residuals
				for(rw = 0; rw < nobs; rw++) {
				  for(xl = 0; xl < dfx; xl++) {
				    Gp[rw][xl] = gpred[rw][xl] + gresid[rw][xl];
				  }
				}
				
				// Recompute p-value for T ~ L|G based on G*
				// fit model T ~ G + L to test L 
				stride = ip + dfx;
				for(rw = 0; rw < nobs; rw++) {
				  designmat[ rw * stride  ] = 1;      // intercept
				  for(xl = 0; xl < dfx; xl++) {
				    designmat[ rw * stride + 1 + xl  ]  = Gp[rw][xl];
				  }
				  for(cl = 0; cl < dfz; cl++) {
				    designmat[ rw * stride + 1 + dfx + cl  ]  = gsl_matrix_get (Lm, rw, cl);
				  }
				}
				
				df = dfz;
				converged = linearRegCompare( pvp, phenovec, designmat, nobs, stride, df );
				if(!converged)Rcpp::Rcout<< "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
				pvp = ( converged ) ? pvp : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G*
				if( pvp > pv ) npos++;
				
				aa = npos < posno;
				cc = nperm < ( maxit - 1 );
				nperm++;
			} // end 'while' permutation loop
		} // End if
		pv = 1.0 * npos / nperm;
		pvec.push_back(pv); // pval for L ind T|G

		pval1[0] = pvec[0]; // pval for T ~ L
		pval2[0] = pvec[1]; // pval for T ~ G|L
		pval4[0] = pvec[2]; // pval for L ind T|G

		pvec.clear();
		gresid.clear();
		gpred.clear();
		LL.clear();
		GG.clear();
		Gp.clear();
		gsl_matrix_free (Lm);
		gsl_matrix_free (Gm);
		gsl_vector_free (Tm);
		
	delete [] designmat;
	delete [] phenovec;
		
	PutRNGstate();

} // End citbinm_linear
