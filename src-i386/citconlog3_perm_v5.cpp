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

void citconlog3p( double *, double *, double *, int &, 
	int &, double *, double *, double *, double *, int &, int  &, int & );
	
double gsl_stats_tss (const double data[], size_t stride, size_t n);

int randwrapper3( int n );

int randwrapper3( int n ) 
{
	int x ;
	x = (int)(n * unif_rand() );
	return x;
}

// conduct permutations individually for each test so an intersection-union type approach can be applied to permutation-FDR
void citconlog3p( double *L, double *G, double *T, int &nrow, int &ncol, 
	double *pval1, double *pval2, double *pval3, double *pval4, 
	int &maxit, int &permit, int &boots )
{
	int rw, brw, cl, i, j, rind, df, df1, df2, npos, nperm, dncol, perm, firstloop;
	int *bootind, *nposperm;
	double rss2, rss3, rss5, F, pv, pvalind, pvp, tmp, rhs;
	double *designmat, *phenovec, *pindep;
	bool aa, bb, cc, converged, permute;
	const int posno = 20;
	vector<vector<double> > LL;
	vector<double> gpred;
	vector<double> gresid;
	vector<double> permindvec;
	gsl_matrix *cov, *X;
	gsl_vector *Gm, *Gp, *c;
	
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
	
	for(rw = 0; rw < nrow; rw++) { 
		permindvec[ rw ] = rw;
	}
	
	// begin permutation loop	
	for(perm = 0; perm < (boots + 1); perm++) {
	
		permute = perm > 0;
		random_shuffle( permindvec.begin(), permindvec.end(), randwrapper3 );
		for(rw = 0; rw < nrow; rw++) {
			bootind[ rw ]  = ( permute ) ? permindvec[ rw ] : rw;
		}
	
		// fit model T ~ L
		// create design matrix with no missing values
		dncol = 1 + ncol;                               // intercept + multiple L variable
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			if( aa ){
				phenovec[ rind ] = T[ rw ];
				designmat[ rind * dncol  ] = 1;      // intercept
				for(cl = 0; cl < ncol; cl++) {
                  designmat[ rind * dncol + 1 + cl  ]  = LL[ brw ][ cl ];
				}
				rind++;
			} // end if aa
		} // end for rw
		
		df = ncol;
		converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
		pv = ( converged ) ? pv : 9;
		pval1[perm] = pv;  // pval for T ~ L, 9 if it did not converge, p1
	
		// fit model T ~ L + G
		dncol = 1 + ncol + 1;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ brw ] != -9999 ) ? aa : 0;
			if( aa ){
				phenovec[ rind ] = T[ rw ];
				designmat[ rind * dncol  ] = 1;      // intercept
				for(cl = 0; cl < ncol; cl++) {
                  designmat[ rind * dncol + 1 + cl  ]  = LL[ rw ][ cl ];
				}
				designmat[ rind * dncol + 1 + ncol  ] = G[ brw ];
				rind++;
			} // end if aa
		} // end for rw
		
		df = 1;
		converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
		pv = ( converged ) ? pv : 9;
		pval2[perm]  = pv;  // pval for T ~ G|L, 9 if it did not converge, p2
		
		// fit model G ~ T
		dncol = 2;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				rind++;
			} // end if aa
		} // end for rw
		
		X = gsl_matrix_alloc(rind, dncol);
		Gm = gsl_vector_alloc (rind);
		
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				gsl_matrix_set(X, rind, 0, 1.0);      // intercept
				gsl_matrix_set(X, rind, 1, T[ rw ]); 
				gsl_vector_set(Gm, rind, G[rw]);
				rind++;
			} // end if aa
		} // end for rw
		
		c = gsl_vector_alloc(dncol);
		cov = gsl_matrix_alloc(dncol, dncol);
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(rind, dncol);
		gsl_multifit_linear(X, Gm, c, cov, &rss2, work);
		gsl_multifit_linear_free(work);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_vector_free(c);
	
		// fit model G ~ L + T
		dncol = 1 + ncol + 1;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				rind++;
			} // end if aa
		} // end for rw
		
		X = gsl_matrix_alloc(rind, dncol);
		Gm = gsl_vector_alloc (rind);
		
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				gsl_matrix_set(X, rind, 0, 1.0);      // intercept
				for(cl = 0; cl < ncol; cl++) {
                  gsl_matrix_set(X, rind, cl + 1, LL[ brw ][ cl ] );
		     	}
				gsl_matrix_set(X, rind, 1 + ncol, T[ rw ]); 
				gsl_vector_set(Gm, rind, G[rw]);
				rind++;
			} // end if aa
		} // end for rw
		
		c = gsl_vector_alloc(dncol);
		cov = gsl_matrix_alloc(dncol, dncol);
		work = gsl_multifit_linear_alloc(rind, dncol);
		gsl_multifit_linear(X, Gm, c, cov, &rss3, work);
		gsl_multifit_linear_free(work);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_vector_free(c);
		df1 = ncol;
		df2 = rind - dncol;
		F = df2*(rss2-rss3)/(rss3*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		pval3[perm]  = pv; // pval for G ~ L|T, p3
	
		// fit model T ~ G + L
		if( perm == 0 ){
			dncol = 1 + 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				brw = bootind[ rw ] ;
				aa = 1;
				aa = ( T[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				aa = ( G[ brw ] != -9999 ) ? aa : 0;
				if( aa ){
					phenovec[ rind ] = T[ rw ];
					designmat[ rind * dncol  ] = 1;      // intercept
					designmat[ rind * dncol + 1  ] = G[ brw ];
					for(cl = 0; cl < ncol; cl++) {
                 	designmat[ rind * dncol + 2 + cl  ]  = LL[ rw ][ cl ];
					}
					rind++;
				} // end if aa
			} // end for rw
		
			df = ncol;
			converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
			pvalind = ( converged ) ? pv : 9;    // p-value for T ~ L|G
		
			// fit model G ~ L
			dncol = 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( G[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				if( aa ){
					rind++;
				} // end if aa
			} // end for rw
		
			X = gsl_matrix_alloc(rind, dncol);
			Gm = gsl_vector_alloc (rind);
		
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( G[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				if( aa ){
					gsl_matrix_set(X, rind, 0, 1.0);      // intercept
					for(cl = 0; cl < ncol; cl++) {
                 	gsl_matrix_set(X, rind, cl + 1, LL[ rw ][ cl ] );
		     		}
					gsl_vector_set(Gm, rind, G[rw]);
					rind++;
				} // end if aa
			} // end for rw
		
			c = gsl_vector_alloc(dncol);
			cov = gsl_matrix_alloc(dncol, dncol);
			work = gsl_multifit_linear_alloc(rind, dncol);
			gsl_multifit_linear(X, Gm, c, cov, &rss5, work);
			gsl_multifit_linear_free(work);
			gsl_matrix_free(cov);
			
			// residuals for G ~ L
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( G[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                 aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				if( aa ){
					rhs = gsl_vector_get(c, 0); // intercept
					for(cl = 0; cl < ncol; cl++) {
                  rhs += gsl_vector_get(c, (cl+1)) * LL[ rw ][ cl ];
		     		}
					gpred.push_back(rhs);
					tmp = gsl_vector_get(Gm, rw) - rhs;
					gresid.push_back(tmp);
				} 
				else {
					gpred.push_back(-9999);
					gresid.push_back(-9999);
				}
			} // end rw loop
			gsl_vector_free (c);
		} // end if perm == 0
			
		// Conduct an initial set of permutations
		if( perm > 0 ){
			// bootstrap the residuals like the other tests, but since the outcome is not bootstrapped, this test is conducted under the null.
			// compute G* based on marginal L effects and permuted residuals
			Gp = gsl_vector_alloc (nrow);
			for(rw = 0; rw < nrow; rw++) {
				brw = bootind[ rw ];
				aa = 1;
				aa = ( gpred[rw] != -9999 ) ? aa : 0;
				aa = ( gresid[brw] != -9999 ) ? aa : 0;
				if( aa ){
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[brw] );
				}
				else {
					gsl_vector_set(Gp, rw, -9999 );
				}
			}
			
			// Recompute p-value for T ~ L|G based on G*
			// fit model T ~ G* + L to test L 
			dncol = 1 + 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( T[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				aa = ( gsl_vector_get(Gp, rw ) != -9999 ) ? aa : 0;
				if( aa ){
					phenovec[ rind ] = T[ rw ];
					designmat[ rind * dncol  ] = 1;      // intercept
					designmat[ rind * dncol + 1  ] = gsl_vector_get(Gp, rw );
					for(cl = 0; cl < ncol; cl++) {
                  		designmat[ rind * dncol + 2 + cl  ]  = LL[ rw ][ cl ];
					}
					rind++;
				} // end if aa
			} // end for rw
		
			df = ncol;
			converged = logisticReg( pvp, phenovec, designmat, rind, dncol, df );
			pindep[ perm - 1 ] = ( converged ) ? pvp : 9;    // p-value for T ~ L|G*

		} // end if perm > 0
	} // End perm loop
		
	npos = 0;
	for(i = 0; i < firstloop; i++){
		// randomly permute residuals
		random_shuffle( gresid.begin(), gresid.end(), randwrapper3 );
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ];
			aa = 1;
			aa = ( gpred[rw] != -9999 ) ? aa : 0;
			aa = ( gresid[brw] != -9999 ) ? aa : 0;
			if( aa ){
				gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
			}
			else {
				gsl_vector_set(Gp, rw, -9999 );
			}
		} // end rw loop
		// fit model T ~ G* + L to test L 
		dncol = 1 + 1 + ncol;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
          	aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( gsl_vector_get(Gp, rw ) != -9999 ) ? aa : 0;
			if( aa ){
				phenovec[ rind ] = T[ rw ];
				designmat[ rind * dncol  ] = 1;      // intercept
				designmat[ rind * dncol + 1  ] = gsl_vector_get(Gp, rw );
				for(cl = 0; cl < ncol; cl++) {
                 designmat[ rind * dncol + 2 + cl  ]  = LL[ rw ][ cl ];
				}
				rind++;
			} // end if aa
		} // end for rw
			
		df = ncol;
		converged = logisticReg( pvp, phenovec, designmat, rind, dncol, df );
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
			random_shuffle( gresid.begin(), gresid.end(), randwrapper3 );
			for(rw = 0; rw < nrow; rw++) {
				brw = bootind[ rw ];
				aa = 1;
				aa = ( gpred[rw] != -9999 ) ? aa : 0;
				aa = ( gresid[brw] != -9999 ) ? aa : 0;
				if( aa ){
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
				}
				else {
					gsl_vector_set(Gp, rw, -9999 );
				}
			} // end rw loop
			// fit model T ~ G* + L to test L 
			dncol = 1 + 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( T[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				aa = ( gsl_vector_get(Gp, rw ) != -9999 ) ? aa : 0;
				if( aa ){
					phenovec[ rind ] = T[ rw ];
					designmat[ rind * dncol  ] = 1;      // intercept
					designmat[ rind * dncol + 1  ] = gsl_vector_get(Gp, rw );
					for(cl = 0; cl < ncol; cl++) {
                  		designmat[ rind * dncol + 2 + cl  ]  = LL[ rw ][ cl ];
					}
					rind++;
				} // end if aa
			} // end for rw
			
			df = ncol;
			converged = logisticReg( pvp, phenovec, designmat, rind, dncol, df );
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
	gsl_vector_free (Gm);
	gsl_vector_free (Gp);
	PutRNGstate();
	
	delete [] bootind;
	delete [] nposperm;
	delete [] designmat;
	delete [] phenovec;
	delete [] pindep;

} // End citconlog3p
} // End extern c wrapper
