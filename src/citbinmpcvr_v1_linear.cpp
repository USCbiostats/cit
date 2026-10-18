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
#include <random> // std::default_random_engine
#include "linearfunc.h"
#include <Rcpp.h>
#include "maxElementWithNan.h"

using namespace Rcpp;
using namespace std;

/*
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: vector of 0/1 variables
Programmer: Joshua Millstein
*/

// conduct permutations individually for each test so an intersection-union type approach can be applied to permutation-FDR

// [[Rcpp::export]]
void citbinmpcvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C,
						int &maxit, int &permit, int &permnum, int &nobs, int &dfz, int &dfx, int &dfc,
						Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, Rcpp::IntegerVector perm_index, int &rseed)
{
	unsigned seed = rseed;
	int rw, brw, cl, cl1, xl, i, j, ip, rind, df, npos, nperm, stride, perm, firstloop;
	// JLM
	// I like initializing pointers to NULL.
	// C++ does not initialize them by default.
	// Check the deallocation statement for a little more info
	int *permind = NULL;
	int *nposperm = NULL;
	double rss, pv, pvalind, pvp, rhs;
	double *designmat = NULL;
	double *phenovec = NULL;
	double *pindep = NULL;
	bool aa, bb, converged, permute;
	const int posno = 20;
	vector<vector<double>> LL;
	vector<vector<double>> XX;
	vector<vector<double>> CC;
	vector<vector<double>> Gp;
	vector<vector<double>> gpred;
	vector<vector<double>> gresid;

	// JLM
	// I don't see allocated anywhere, I assume you can set them to NULL
	// but I decided to create them as 10x10 matrices.
	gsl_matrix *Lm = gsl_matrix_alloc(10, 10);
	gsl_matrix *Gm = gsl_matrix_alloc(10, 10);
	gsl_matrix *cov = gsl_matrix_alloc(10, 10);
	gsl_matrix *D = gsl_matrix_alloc(10, 10);
	// I'm assuming the following would work but the deallocation code
	// below would have to change
	// gsl_matrix *Lm = NULL;
	// gsl_matrix *Gm = NULL;
	// gsl_matrix *cov = NULL;
	// gsl_matrix *D = NULL;
	// JLM
	// Same issue with the vectors. Set to length 10
	gsl_vector *c = gsl_vector_alloc(10);
	gsl_vector *Gmj = gsl_vector_alloc(10);
	// gsl_vector *c = NULL;
	// gsl_vector Gmj = NULL;


	designmat = new double[nobs * (dfc + dfz + dfx + 1)];
	phenovec = new double[nobs];
	permind = new int[nobs];
	nposperm = new int[permnum];
	pindep = new double[permnum];

	for (i = 0; i < permnum; i++)
	{
		nposperm[i] = 0;
	}

	firstloop = permit;

	LL.resize(nobs);
	XX.resize(nobs);
	CC.resize(nobs);
	Gp.resize(nobs);
	gpred.resize(nobs);
	gresid.resize(nobs);

	GetRNGstate();

	for (rw = 0; rw < nobs; rw++)
	{
		CC[rw].resize(dfc);
		LL[rw].resize(dfz);
		XX[rw].resize(dfx);
		Gp[rw].resize(dfx);
		gpred[rw].resize(dfx);
		gresid[rw].resize(dfx);
	}

	for (cl = 0; cl < dfz; cl++)
	{
		for (rw = 0; rw < nobs; rw++)
		{
			LL[rw][cl] = L[rw + nobs * cl];
		}
	}
	for (cl = 0; cl < dfx; cl++)
	{
		for (rw = 0; rw < nobs; rw++)
		{
			XX[rw][cl] = G[rw + nobs * cl];
		}
	}
	for (cl = 0; cl < dfc; cl++)
	{
		for (rw = 0; rw < nobs; rw++)
		{
			CC[rw][cl] = C[rw + nobs * cl];
		}
	}

	// begin permutation loop
	for (perm = 0; perm < (permnum + 1); perm++)
	{

		permute = perm > 0;
		for (rw = 0; rw < nobs; rw++)
		{
			permind[rw] = (permute) ? perm_index[(perm - 1) * nobs + rw] : rw;
		}

		// fit model T ~ C + L
		// create design matrix with no missing values
		stride = 1 + dfc + dfz; // intercept + covariates + multiple L variable
		rind = 0;
		for (rw = 0; rw < nobs; rw++)
		{
			brw = permind[rw];
			aa = 1;
			aa = (T[rw] != -9999) ? aa : 0;
			for (cl = 0; cl < dfz; cl++)
			{
				aa = (LL[brw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfc; cl++)
			{
				aa = (CC[brw][cl] != -9999) ? aa : 0;
			}
			if (aa)
			{
				phenovec[rw] = T[rw];
				designmat[rw * stride] = 1; // intercept
				for (cl = 0; cl < dfc; cl++)
				{
					designmat[rw * stride + 1 + cl] = CC[brw][cl];
				}
				for (cl = 0; cl < dfz; cl++)
				{
					designmat[rw * stride + 1 + dfc + cl] = LL[brw][cl];
				}
				rind++;
			} // end if aa
		} // end for rw

		df = dfz;
		converged = linearRegCompare(pv, phenovec, designmat, rind, stride, df);
		if (!converged)
			Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN();
		pval1[perm] = pv; // pval for T ~ C + L, 9 if it did not converge, p1

		// fit model T ~ C + L + G
		stride = 1 + dfc + dfz + dfx;
		rind = 0;
		for (rw = 0; rw < nobs; rw++)
		{
			brw = permind[rw];
			aa = 1;
			aa = (T[rw] != -9999) ? aa : 0;
			for (cl = 0; cl < dfc; cl++)
			{
				aa = (CC[rw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfz; cl++)
			{
				aa = (LL[rw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfx; cl++)
			{
				aa = (XX[rw][cl] != -9999) ? aa : 0;
			}
			if (aa)
			{
				phenovec[rw] = T[rw];
				designmat[rw * stride] = 1; // intercept
				for (cl = 0; cl < dfc; cl++)
				{
					designmat[rw * stride + 1 + cl] = CC[rw][cl];
				}
				for (cl = 0; cl < dfz; cl++)
				{
					designmat[rw * stride + 1 + dfc + cl] = LL[rw][cl];
				}
				for (cl = 0; cl < dfx; cl++)
				{
					designmat[rw * stride + 1 + dfc + dfz + cl] = XX[brw][cl];
				}
				rind++;
			} // end if aa
		} // end for rw

		df = dfx;
		converged = linearRegCompare(pv, phenovec, designmat, rind, stride, df);
		if (!converged)
			Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN();
		pval2[perm] = pv; // pval for T ~ G|L, 9 if it did not converge, p2

		// fit model T ~ C + G + L to test L
		stride = 1 + dfc + dfx + dfz;
		rind = 0;
		for (rw = 0; rw < nobs; rw++)
		{
			brw = permind[rw];
			aa = 1;
			aa = (T[rw] != -9999) ? aa : 0;
			for (cl = 0; cl < dfc; cl++)
			{
				aa = (CC[brw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfz; cl++)
			{
				aa = (LL[brw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfx; cl++)
			{
				aa = (XX[rw][cl] != -9999) ? aa : 0;
			}
			if (aa)
			{
				phenovec[rw] = T[rw];
				designmat[rw * stride] = 1; // intercept
				for (cl = 0; cl < dfc; cl++)
				{
					designmat[rw * stride + 1 + cl] = CC[rw][cl];
				}
				for (cl = 0; cl < dfx; cl++)
				{
					designmat[rw * stride + 1 + dfc + cl] = XX[rw][cl];
				}
				for (cl = 0; cl < dfz; cl++)
				{
					designmat[rw * stride + 1 + dfc + dfx + cl] = LL[brw][cl];
				}
				rind++;
			} // end if aa
		} // end for rw

		df = dfz;
		converged = linearRegCompare(pv, phenovec, designmat, rind, stride, df);
		if (!converged)
			Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|C,G
		pval3nc[perm] = pv;												  // pvalue to be used for non-centrality parameter

		// estimate residuals and predicted values only using non-permuted data
		if (perm == 0)
		{

			pvalind = (converged) ? pv : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|C,G

			// fit model G ~ L
			rind = 0;
			for (rw = 0; rw < nobs; rw++)
			{
				aa = 1;
				for (cl = 0; cl < dfx; cl++)
				{
					aa = (XX[rw][cl] != -9999) ? aa : 0;
				}
				for (cl = 0; cl < dfz; cl++)
				{
					aa = (LL[rw][cl] != -9999) ? aa : 0;
				}
				if (aa)
				{
					rind++;
				} // end if aa
			} // end for rw

			ip = 1 + dfz;
			D = gsl_matrix_alloc(rind, ip);
			Lm = gsl_matrix_alloc(rind, dfz);
			Gm = gsl_matrix_alloc(rind, dfx);
			Gmj = gsl_vector_alloc(rind);
			rind = 0;
			for (rw = 0; rw < nobs; rw++)
			{
				aa = 1;
				for (cl = 0; cl < dfx; cl++)
				{
					aa = (XX[rw][cl] != -9999) ? aa : 0;
				}
				for (cl = 0; cl < dfz; cl++)
				{
					aa = (LL[rw][cl] != -9999) ? aa : 0;
				}
				if (aa)
				{
					gsl_matrix_set(D, rind, 0, 1.0);
					for (cl = 0; cl < dfz; cl++)
					{
						gsl_matrix_set(D, rind, cl + 1, LL[rw][cl]);
					}
					for (cl = 0; cl < dfx; cl++)
					{
						gsl_matrix_set(Gm, rind, cl, XX[rw][cl]);
					}
					rind++;
				} // end if aa
			} // end for rw

			// get linear predictors and residuals
			for (xl = 0; xl < dfx; xl++)
			{
				c = gsl_vector_alloc(ip);
				cov = gsl_matrix_alloc(ip, ip);
				gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(rind, ip);
				for (rw = 0; rw < rind; rw++)
				{
					gsl_vector_set(Gmj, rw, gsl_matrix_get(Gm, rw, xl));
				}
				gsl_multifit_linear(D, Gmj, c, cov, &rss, work);
				gsl_multifit_linear_free(work);
				gsl_matrix_free(cov);

				// residuals for G ~ L
				for (rw = 0; rw < nobs; rw++)
				{
					rhs = 0;
					aa = 1;
					for (cl = 0; cl < dfz; cl++)
					{
						aa = (LL[rw][cl] != -9999) ? aa : 0;
					}
					if (aa)
					{
						rhs = gsl_vector_get(c, 0); // intercept
						for (cl = 0; cl < dfz; cl++)
						{
							cl1 = cl + 1;
							rhs += gsl_vector_get(c, cl1) * LL[rw][cl];
						}
						gpred[rw][xl] = rhs;
						gresid[rw][xl] = XX[rw][xl] - rhs;
					} // end if aa
					else
					{
						gpred[rw][xl] = -9999;
						gresid[rw][xl] = -9999;
					}
				} // end rw loop
			} // end xl loop
			gsl_vector_free(c);
		} // end if perm == 0

		// Conduct an initial set of permutations
		if (perm > 0)
		{
			// compute G* based on marginal L effects and permuted residuals
			for (rw = 0; rw < nobs; rw++)
			{
				brw = permind[rw];
				aa = 1;
				for (xl = 0; xl < dfx; xl++)
				{
					aa = (gpred[rw][xl] != -9999) ? aa : 0;
					aa = (gresid[brw][xl] != -9999) ? aa : 0;
				}
				if (aa)
				{
					for (xl = 0; xl < dfx; xl++)
					{
						Gp[rw][xl] = gpred[rw][xl] + gresid[brw][xl];
					}
				}
				else
				{
					for (xl = 0; xl < dfx; xl++)
					{
						Gp[rw][xl] = -9999;
					}
				} // end if aa
			} // end rw loop

			// Recompute p-value for T ~ L|G based on G*
			// fit model T ~ G* + L to test L
			stride = 1 + dfc + dfx + dfz;
			rind = 0;
			for (rw = 0; rw < nobs; rw++)
			{
				aa = 1;
				aa = (T[rw] != -9999) ? aa : 0;
				for (cl = 0; cl < dfc; cl++)
				{
					aa = (CC[rw][cl] != -9999) ? aa : 0;
				}
				for (cl = 0; cl < dfz; cl++)
				{
					aa = (LL[rw][cl] != -9999) ? aa : 0;
				}
				for (cl = 0; cl < dfx; cl++)
				{
					aa = (Gp[rw][cl] != -9999) ? aa : 0;
				}
				if (aa)
				{
					phenovec[rind] = T[rw];
					designmat[rind * stride] = 1; // intercept
					for (cl = 0; cl < dfc; cl++)
					{
						designmat[rind * stride + 1 + cl] = CC[rw][cl];
					}
					for (cl = 0; cl < dfx; cl++)
					{
						designmat[rind * stride + 1 + dfc + cl] = Gp[rw][cl];
					}
					for (cl = 0; cl < dfz; cl++)
					{
						designmat[rind * stride + 1 + dfc + dfx + cl] = LL[rw][cl];
					}
					rind++;
				} // end if aa
			} // end for rw

			df = dfz;
			converged = linearRegCompare(pvp, phenovec, designmat, rind, stride, df);
			if (!converged)
				Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
			pindep[perm - 1] = (converged) ? pvp : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|G*
		} // end if perm > 0
	} // End perm loop

	npos = 0;
	for (i = 0; i < firstloop; i++)
	{
		// randomly permute residuals

		shuffle(gresid.begin(), gresid.end(), std::default_random_engine(seed));
		seed += 1;
		for (rw = 0; rw < nobs; rw++)
		{
			aa = 1;
			for (xl = 0; xl < dfx; xl++)
			{
				aa = (gpred[rw][xl] != -9999) ? aa : 0;
				aa = (gresid[brw][xl] != -9999) ? aa : 0;
			}
			if (aa)
			{
				for (xl = 0; xl < dfx; xl++)
				{
					Gp[rw][xl] = gpred[rw][xl] + gresid[rw][xl];
				}
			}
			else
			{
				for (xl = 0; xl < dfx; xl++)
				{
					Gp[rw][xl] = -9999;
				}
			} // end if aa
		} // end rw loop

		// Recompute p-value for T ~ L|C,G*
		// fit model T ~ C + G* + L to test L
		stride = 1 + dfc + dfx + dfz;
		rind = 0;
		for (rw = 0; rw < nobs; rw++)
		{
			aa = 1;
			aa = (T[rw] != -9999) ? aa : 0;
			for (cl = 0; cl < dfc; cl++)
			{
				aa = (CC[rw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfz; cl++)
			{
				aa = (LL[rw][cl] != -9999) ? aa : 0;
			}
			for (cl = 0; cl < dfx; cl++)
			{
				aa = (Gp[rw][cl] != -9999) ? aa : 0;
			}
			if (aa)
			{
				phenovec[rind] = T[rw];
				designmat[rind * stride] = 1; // intercept
				for (cl = 0; cl < dfc; cl++)
				{
					designmat[rind * stride + 1 + cl] = CC[rw][cl];
				}
				for (cl = 0; cl < dfx; cl++)
				{
					designmat[rind * stride + 1 + dfc + cl] = Gp[rw][cl];
				}
				for (cl = 0; cl < dfz; cl++)
				{
					designmat[rind * stride + 1 + dfc + dfx + cl] = LL[rw][cl];
				}
				rind++;
			} // end if aa
		} // end for rw

		df = dfz;
		converged = linearRegCompare(pvp, phenovec, designmat, rind, stride, df);
		if (!converged)
			Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pvp = (converged) ? pvp : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|C,G*
		if (pvp > pvalind)
			npos++;

		for (j = 0; j < permnum; j++)
		{
			if (pvp > pindep[j])
				nposperm[j] += 1;
		}
	} // end initial i permutation loop

	nperm = firstloop;

	// Conduct additional permutations if there is some indication of statistical significance
	if (permnum == 0)
	{
		aa = npos < posno;
		bb = nperm < maxit;

		while (aa && bb)
		{

			// randomly permute residuals

			shuffle(gresid.begin(), gresid.end(), std::default_random_engine(seed));
			seed += 1;
			for (rw = 0; rw < nobs; rw++)
			{
				aa = 1;
				for (xl = 0; xl < dfx; xl++)
				{
					aa = (gpred[rw][xl] != -9999) ? aa : 0;
					aa = (gresid[brw][xl] != -9999) ? aa : 0;
				}
				if (aa)
				{
					for (xl = 0; xl < dfx; xl++)
					{
						Gp[rw][xl] = gpred[rw][xl] + gresid[rw][xl];
					}
				}
				else
				{
					for (xl = 0; xl < dfx; xl++)
					{
						Gp[rw][xl] = -9999;
					}
				} // end if aa
			} // end rw loop

			// Recompute p-value for T ~ L|C,G based on G*
			// fit model T ~ C + G* + L to test L
			stride = 1 + dfc + dfx + dfz;
			rind = 0;
			for (rw = 0; rw < nobs; rw++)
			{
				aa = 1;
				aa = (T[rw] != -9999) ? aa : 0;
				for (cl = 0; cl < dfc; cl++)
				{
					aa = (CC[rw][cl] != -9999) ? aa : 0;
				}
				for (cl = 0; cl < dfz; cl++)
				{
					aa = (LL[rw][cl] != -9999) ? aa : 0;
				}
				for (cl = 0; cl < dfx; cl++)
				{
					aa = (Gp[rw][cl] != -9999) ? aa : 0;
				}
				if (aa)
				{
					phenovec[rind] = T[rw];
					designmat[rind * stride] = 1; // intercept
					for (cl = 0; cl < dfc; cl++)
					{
						designmat[rind * stride + 1 + cl] = CC[rw][cl];
					}
					for (cl = 0; cl < dfx; cl++)
					{
						designmat[rind * stride + 1 + dfc + cl] = Gp[rw][cl];
					}
					for (cl = 0; cl < dfz; cl++)
					{
						designmat[rind * stride + 1 + dfc + dfx + cl] = LL[rw][cl];
					}
					rind++;
				} // end if aa
			} // end for rw

			df = dfz;
			converged = linearRegCompare(pvp, phenovec, designmat, rind, stride, df);
			if (!converged)
				Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
			pvp = (converged) ? pvp : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|C,G*
			if (pvp > pvalind)
				npos++;

			for (j = 0; j < permnum; j++)
			{
				if (pvp > pindep[j])
					nposperm[j] += 1;
			}

			aa = npos < posno;
			nperm++;
		} // end 'while' permutation loop
	} // end if permnum == 0

	for (perm = 0; perm < (permnum + 1); perm++)
	{
		pval4[perm] = (perm == 0) ? 1.0 * npos / nperm : 1.0 * nposperm[perm - 1] / nperm; // pval for L ind T|C,G
	} // End perm loop

	LL.clear();
	XX.clear();
	CC.clear();
	Gp.clear();
	gpred.clear();
	gresid.clear();

	if (Lm)
		gsl_matrix_free(Lm);
	if (Gm)
		gsl_matrix_free(Gm);
	if (D)
		gsl_matrix_free(D);
	if (Gmj)
		gsl_vector_free(Gmj);

	PutRNGstate();

	// JLM
	// I don't like using delete without checking the pointer is set to NULL
	// If an allocation fails the pointer is set to NULL. The pointer needs
	// to be initialized to NULL to avoid an error if the pointer is never
	// set to value.
	if (permind)
		delete[] permind;
	if (nposperm)
		delete[] nposperm;
	if (designmat)
		delete[] designmat;
	if (phenovec)
		delete[] phenovec;
	if (pindep)
		delete[] pindep;

} // End citbinmpcvr_linear
