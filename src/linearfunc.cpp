#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<list>
#include<string.h>
#include<map>
#include<set>
#include <numeric>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multifit.h>
#include "linearfunc.h"



using namespace std;

/*
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: matrix of 0/1 variables
Programmer: Joshua Millstein
*/

bool linearReg(double * phenovec_filtered, double * designmat_filtered, int & samplesize, int & rank, int & df, double * SR){
    gsl_matrix_view X_view = gsl_matrix_view_array(designmat_filtered, samplesize, rank);
    gsl_matrix_view subX_view = gsl_matrix_submatrix(&X_view.matrix, 0, 0, samplesize, rank - df);
    gsl_matrix* X = &(subX_view.matrix);

    gsl_matrix *XT= gsl_matrix_alloc(rank-df, samplesize);
    gsl_matrix_transpose_memcpy(XT, X);

    gsl_matrix *XTX= gsl_matrix_alloc(rank-df, rank-df);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, XT, X, 0.0, XTX);

    gsl_matrix_view Y_view = gsl_matrix_view_array(phenovec_filtered,samplesize,1);
    gsl_matrix* Y = &(Y_view.matrix);

    gsl_matrix *XTY= gsl_matrix_alloc(rank-df, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, XT, Y, 0.0, XTY); // bug here

    gsl_matrix *invert_XTX= gsl_matrix_alloc(rank-df, rank-df);
    gsl_matrix_memcpy(invert_XTX, XTX);
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    int returncode = gsl_linalg_cholesky_decomp1(invert_XTX);
    gsl_set_error_handler (old_handler);
    if (returncode!=0){
        cout<<"cannot find inversed matrix"<<endl;
        gsl_matrix_free(XT);
        gsl_matrix_free(XTX);
        gsl_matrix_free(XTY);
        gsl_matrix_free(invert_XTX);
        return false;
    }
    gsl_linalg_cholesky_invert(invert_XTX);


    gsl_matrix *Beta= gsl_matrix_alloc(rank-df, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invert_XTX, XTY, 0.0, Beta);

    // cout<<"When df is "<<df<<", Beta array is:"<<endl;
    // for(int i=0;i<rank-df;i++){
    //     cout<<gsl_matrix_get(Beta,i,0)<< endl;
    // }

    gsl_matrix *Xb= gsl_matrix_alloc(samplesize, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, Beta, 0.0, Xb);


    for(int i=0;i<samplesize;i++){
        SR[i]=gsl_matrix_get(Y,i,0) - gsl_matrix_get(Xb,i,0);
        SR[i]=SR[i]*SR[i];
    }


    gsl_matrix_free(XT);
    gsl_matrix_free(XTX);
    gsl_matrix_free(XTY);
    gsl_matrix_free(invert_XTX);
    gsl_matrix_free(Beta);
    gsl_matrix_free(Xb);
    return true;




}


bool linearRegCompare( double & pvalue, double * phenovec_filtered, double * designmat_filtered, int & samplesize, int & rank, int & df){
    //gsl_multifit_linear_workspace a=*gsl_multifit_linear_alloc(5,5);
    int a=0;
    double *SR_FULL = new double[samplesize];
    bool B1=linearReg(phenovec_filtered, designmat_filtered, samplesize, rank, a, SR_FULL);
    double *SR_LIMIT = new double[samplesize];
    bool B2=linearReg(phenovec_filtered, designmat_filtered, samplesize, rank, df, SR_LIMIT);

    if(!B1 || !B2) {
        delete[] SR_FULL;
        delete[] SR_LIMIT;
        return false;
    }
    double RSS1=0.0, RSS2=0.0;
    for(int i=0;i<samplesize;i++){
        RSS1=RSS1+SR_LIMIT[i];
        RSS2=RSS2+SR_FULL[i];
    }
    double F_statistics=((RSS1-RSS2)/(df))/(RSS2/(samplesize-rank));
    pvalue=gsl_cdf_fdist_Q(F_statistics,static_cast<double>(df),static_cast<double>(samplesize-rank));
   //cout<<"RSS1-RSS2: "<<RSS1-RSS2<<"    F_statistics: "<<F_statistics<<"    df1: "<<df<<"    df2: "<<samplesize-rank<<"    pvalue: "<<pvalue<<endl;


    delete[] SR_FULL;
    delete[] SR_LIMIT;
    return true;

}