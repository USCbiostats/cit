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

#include "logisticfunc.h"

using namespace std;

/*
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: matrix of 0/1 variables
Programmer: Joshua Millstein
*/

extern "C" {

void citconlog2p( double *, double *, double *, int &, 
	int &, double *, double *, double *, double *, int &, int &, int &);
	
double gsl_stats_tss (const double data[], size_t stride, size_t n);

int randwrapper( int n );

int randwrapper( int n ) 
{
	int x ;
	x = (int)(n * unif_rand() );
	return x;
}

// conduct permutations individually for each test so an intersection-union type approach can be applied to permutation-FDR
void citconlog2p( double *L, double *G, double *T, int &nrow, int &ncol, 
	double *pval1, double *pval2, double *pval3, double *pval4, 
	int &maxit, int &permit, int &boots )
{
	int rw, brw, cl, i, j, rind, df, df1, df2, nobs, ip, npos, nperm, nmiss, stride, perm, firstloop;
	int *bootind, *nposperm;
	double rss2, rss3, rss5, F, pv, pvalind, pvp, tmp, rhs;
	double *designmat, *phenovec, *pindep;
	bool aa, bb, cc, converged, permute;
	const int posno = 20;
	vector<vector<double> > LL;
	vector<double> gpred;
	vector<double> gresid;
	vector<double> permindvec;
	gsl_matrix *Lm, *cov, *X;
	gsl_vector *Gm, *Tm, *Gp, *c;

	bootind = new int[nrow];
	nposperm = new int[boots];
	designmat = new double[nrow * (ncol + 2)];
	phenovec = new double[nrow];
	pindep = new double[boots];
	
	for(i = 0; i < boots; i++){
		nposperm[ i ] = 0;
	}	

	firstloop = permit;
    
	LL.resize( nrow );
	GetRNGstate();
	
	for(rw = 0; rw < nrow; rw++) {
		LL[rw].resize( ncol );	
	}
	
	for(cl = 0; cl < ncol; cl++) {
		for(rw = 0; rw < nrow; rw++) {
			LL[rw][cl] = L[rw + nrow * cl];
		}
	}
	
	// create analysis vectors w/no missing data
	nobs = 0;
	for(rw = 0; rw < nrow; rw++) {
		nmiss = 0;
		for(cl = 0; cl < ncol; cl++) {
			if( LL[rw][cl] == -9999 ) {
				nmiss++;
			 }                                                
		}
		aa = nmiss == 0;
		bb = G[rw] != -9999;
		cc = T[rw] != -9999;
		if(aa && bb && cc) {
			nobs++;
		}
	} // end rw loop           

	Lm = gsl_matrix_alloc (nobs, ncol);
	Gm = gsl_vector_alloc (nobs);
	Gp = gsl_vector_alloc (nobs);
	Tm = gsl_vector_alloc (nobs);
	permindvec.resize( nobs );
	
	rind = 0;
	for(rw = 0; rw < nrow; rw++) {
		nmiss = 0;
		for(cl = 0; cl < ncol; cl++) {
			if( LL[rw][cl] == -9999 ) {
				nmiss++;
			}                                                
		} // end cl loop
		aa = nmiss == 0;
		bb = G[rw] != -9999;
		cc = T[rw] != -9999;			
			
		if(aa && bb && cc) {
			for(cl = 0; cl < ncol; cl++) {
                  gsl_matrix_set(Lm, rind, cl, LL[rw][cl]);
			}
			gsl_vector_set(Gm, rind, G[rw]);
			gsl_vector_set(Tm, rind, T[rw]);
			rind++;
		} // end if
	} // end rw loop
	
	for(rw = 0; rw < nobs; rw++) { 
		permindvec[ rw ] = rw;
	}

	// begin bootstrap loop	
	for(perm = 0; perm < (boots + 1); perm++) {
		
		// bootstrap vector of row indicators to be used when perm > 0
		// sample indicators with replacement
		permute = perm > 0;
		random_shuffle( permindvec.begin(), permindvec.end(), randwrapper );
		for(rw = 0; rw < nobs; rw++) {
			bootind[ rw ]  = ( permute ) ? permindvec[ rw ] : rw;
		}
      
		// fit model T ~ L
		ip = 1 + ncol;                               // intercept + multiple L variable
		for(rw = 0; rw < nobs; rw++) {
		   brw = bootind[ rw ] ;
		   phenovec[ rw ] = gsl_vector_get(Tm, rw );
		   designmat[ rw * ip  ] = 1;      // intercept
			for(cl = 0; cl < ncol; cl++) {
                  designmat[ rw * ip + 1 + cl  ]  = gsl_matrix_get (Lm, brw, cl);
			}
		} // end rw loop
		
		df = ncol;
		converged = logisticReg( pv, phenovec, designmat, nobs, ip, df );
		pv = ( converged ) ? pv : 9;
		pval1[perm] = pv;  // pval for T ~ L, 9 if it did not converge, p1
		
		// fit model T ~ L + G
		stride = ip + 1;
		for(rw = 0; rw < nobs; rw++) {
			brw = bootind[ rw ] ;
			designmat[ rw * stride ] = 1;      // intercept
			for(cl = 0; cl < ncol; cl++) {
				designmat[ rw * stride + 1 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
			}
			designmat[ rw * stride + 1 + ncol  ] = gsl_vector_get(Gm, brw );
		}
		
		df = 1;
		converged = logisticReg( pv, phenovec, designmat, nobs, stride, df );
		pv = ( converged ) ? pv : 9;
		pval2[perm]  = pv;  // pval for T ~ G|L, 9 if it did not converge, p2
		
		// fit model G ~ T
		X = gsl_matrix_alloc (nobs,2);
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(X, rw, 0, 1.0);  // intercept
			gsl_matrix_set(X, rw, 1, gsl_vector_get(Tm, rw));
		}
		c = gsl_vector_alloc(2);
		cov = gsl_matrix_alloc(2, 2);
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(nobs, 2);
		gsl_multifit_linear(X, Gm, c, cov, &rss2, work); 
		gsl_multifit_linear_free (work);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_vector_free(c);

		// fit model G ~ L + T
		X = gsl_matrix_alloc(nobs, ip + 1);
		for(rw = 0; rw < nobs; rw++) {
			brw = bootind[ rw ] ;
			gsl_matrix_set(X, rw, 0, 1.0);      // intercept
			for(cl = 0; cl < ncol; cl++) {
                  gsl_matrix_set(X, rw, cl + 1, gsl_matrix_get(Lm, brw, cl));
		     }
		     gsl_matrix_set(X, rw, ip, gsl_vector_get(Tm, rw)); 
		}
		c = gsl_vector_alloc(ip + 1);
		cov = gsl_matrix_alloc(ip + 1, ip + 1);
		work = gsl_multifit_linear_alloc(nobs, ip + 1);
		gsl_multifit_linear(X, Gm, c, cov, &rss3, work);
		gsl_multifit_linear_free(work);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_vector_free(c);
		df1 = ncol;
		df2 = nobs - ip -1;
		F = df2*(rss2-rss3)/(rss3*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		pval3[perm]  = pv; // pval for G ~ L|T, p3

         if( perm == 0 ){
			// fit model T ~ G + L to test L 
			stride = ip + 1;
			for(rw = 0; rw < nobs; rw++) {
		   		designmat[ rw * stride  ] = 1;      // intercept
		   		designmat[ rw * stride + 1  ] = gsl_vector_get(Gm, rw );
				for(cl = 0; cl < ncol; cl++) {
          			designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   		}
		   		phenovec[ rw ] = gsl_vector_get(Tm, rw );
			}
		
			df = ncol;
			converged = logisticReg( pv, phenovec, designmat, nobs, stride, df );
			pvalind = ( converged ) ? pv : 9;    // p-value for T ~ L|G

			// fit model G ~ L
			X = gsl_matrix_alloc (nobs, ip );
			for(rw = 0; rw < nobs; rw++) {
				gsl_matrix_set(X, rw, 0, 1.0);      // intercept
				for(cl = 0; cl < ncol; cl++) {
                 	gsl_matrix_set(X, rw, cl + 1, gsl_matrix_get (Lm, rw, cl));
		     	}
			}
			c = gsl_vector_alloc (ip);
			cov = gsl_matrix_alloc (ip, ip);
			work = gsl_multifit_linear_alloc(nobs, ip);
			gsl_multifit_linear(X, Gm, c, cov, &rss5, work);
			gsl_multifit_linear_free(work);
			gsl_matrix_free(cov);
			
			// residuals for G ~ L
			for(rw = 0; rw < nobs; rw++) {
				rhs = 0;
				for(cl = 0; cl < ip; cl++) {
                  	rhs += gsl_vector_get(c, cl) * gsl_matrix_get(X, rw, cl);
		     	}
			
				gpred.push_back(rhs);
				tmp = gsl_vector_get(Gm, rw) - rhs;
				gresid.push_back(tmp);
			} // end rw loop
			gsl_vector_free (c);
		} // end if perm == 0
		
		// Conduct an initial set of permutations
		
		if( perm > 0 ){
			// bootstrap the residuals like the other tests, but since the outcome is not bootstrapped, this test is conducted under the null.
			// compute G* based on marginal L effects and permuted residuals
			for(rw = 0; rw < nobs; rw++) {
				brw = bootind[ rw ];
				gsl_vector_set(Gp, rw, gpred[rw] + gresid[brw] );
			}
			
			// Recompute p-value for T ~ L|G based on G*
			// fit model T ~ G* + L to test L 
			stride = ip + 1;
			for(rw = 0; rw < nobs; rw++) {
				designmat[ rw * stride  ] = 1;      // intercept
				designmat[ rw * stride + 1  ] = gsl_vector_get(Gp, rw );
				for(cl = 0; cl < ncol; cl++) {
					designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
				}
			} // end rw loop
		
			df = ncol;
			converged = logisticReg( pvp, phenovec, designmat, nobs, stride, df );
			pindep[ perm - 1 ] = ( converged ) ? pvp : 9;    // p-value for T ~ L|G*

		} // end if perm > 0
	} // End perm loop
		
	npos = 0;
	for(i = 0; i < firstloop; i++){
				// randomly permute residuals
				random_shuffle( gresid.begin(), gresid.end(), randwrapper );
			
				// compute G* based on marginal L effects and permuted residuals
				for(rw = 0; rw < nobs; rw++) {
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
				}
			
				// Recompute p-value for T ~ L|G based on G*
				// fit model T ~ G* + L to test L 
				stride = ip + 1;
				for(rw = 0; rw < nobs; rw++) {
					designmat[ rw * stride  ] = 1;      // intercept
					designmat[ rw * stride + 1  ] = gsl_vector_get(Gp, rw );
					for(cl = 0; cl < ncol; cl++) {
						designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
					}
				} // end rw loop
		
				df = ncol;
				converged = logisticReg( pvp, phenovec, designmat, nobs, stride, df );
				pvp = ( converged ) ? pvp : 9;    // p-value for T ~ L|G*
				if( pvp > pvalind ) npos++;
			
				for(j = 0; j < boots; j++){
					if( pvp > pindep[ j ] ) nposperm[ j ] += 1;
				}
	} // end initial i permutation loop
		
	nperm = firstloop;
			
			// Conduct additional permutations if there is some indication of statistical significance
	if( boots == 0 ){
				aa = npos < posno;
				bb = nperm < maxit;
				
				while(aa && bb) {
	
					// randomly permute residuals
					random_shuffle( gresid.begin(), gresid.end(), randwrapper );
				
					// compute G* based on marginal L effects and permuted residuals
					for(rw = 0; rw < nobs; rw++) {
						gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
					}
				
					// Recompute p-value for T ~ L|G based on G*
					// fit model T ~ G* + L to test L 
					stride = ip + 1;
					for(rw = 0; rw < nobs; rw++) {
		   				designmat[ rw * stride  ] = 1;      // intercept
		   				designmat[ rw * stride + 1  ] = gsl_vector_get(Gp, rw );
						for(cl = 0; cl < ncol; cl++) {
          					designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   				}
					} // end rw loop
		
					df = ncol;
					converged = logisticReg( pvp, phenovec, designmat, nobs, stride, df );
					pvp = ( converged ) ? pvp : 9;    // p-value for T ~ L|G*
					if( pvp > pvalind ) npos++;
				
					for(j = 0; j < boots; j++){
			   			if( pvp > pindep[ j ] ) nposperm[ j ] += 1;
					}
				
					aa = npos < posno;
					cc = nperm < ( maxit - 1 );
					nperm++;
				} // end 'while' permutation loop
	} // end if boots == 0 
	
	for(perm = 0; perm < (boots + 1); perm++) {
		pv = ( perm == 0 ) ? 1.0 * npos / nperm : 1.0 * nposperm[ perm - 1 ] / nperm;
		
		// To avoid a p-value = 0, make 0 p-values = 1/nperm
		pval4[perm] = ( pv > 0 ) ? pv : 1.0 * 1.0 / nperm;   // pval for L ind T|G
		
	} // End perm loop

	gresid.clear();
	gpred.clear();
	LL.clear();
	gsl_matrix_free (Lm);
	gsl_vector_free (Gm);
	gsl_vector_free (Tm);
	gsl_vector_free (Gp);
	PutRNGstate();
	
	delete [] bootind;
	delete [] nposperm;
	delete [] designmat;
	delete [] phenovec;
	delete [] pindep;

} // End citconlog2p
} // End extern c wrapper
