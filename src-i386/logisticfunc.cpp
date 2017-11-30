#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<list>
#include<string.h>
#include<map>
#include<set>

#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_statistics_double.h>
#include "logisticfunc.h"

using namespace std;

typedef short unsigned int suint;
typedef unsigned int uint;

bool chDecomp(const double * inputMatrix, double * decompMat,uint dim){
    for(uint i=0;i<dim*dim;++i) decompMat[i] = inputMatrix[i];
    gsl_matrix_const_view A = gsl_matrix_const_view_array(inputMatrix,dim,dim);
    gsl_matrix_view B = gsl_matrix_view_array(decompMat,dim,dim);
    gsl_matrix_memcpy(&B.matrix,&A.matrix);
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    int returncode = gsl_linalg_cholesky_decomp(&B.matrix);
    gsl_set_error_handler (old_handler);
    if (returncode!=0){
        return false;
    }
    return true;
}

void invert(const double * decompMat,double * invertedMat,uint dim){
    gsl_matrix_const_view decomp = gsl_matrix_const_view_array(decompMat,dim,dim);
    gsl_matrix_view invmat = gsl_matrix_view_array(invertedMat,dim,dim);
    gsl_matrix_set_identity(&invmat.matrix);
    for(uint i=0;i<dim;++i){
        gsl_vector_view x = gsl_matrix_column (&invmat.matrix, i);
        gsl_linalg_cholesky_svx(&decomp.matrix,&x.vector);
    }
}

void scoreTest(double & chiSq,double & logL,double * diseaseStatus,
double * effSizeMat,double * betas,double * invInfoMatrix,uint iObsSampleSize,
suint stride,suint params,int & count,int MAX){
    chiSq = logL = 0.;
    double *U = new double[params];
    double *infoMatrix = new double[params*params]; 
    for(uint j=0;j<params;++j){
        U[ j ] = 0;
    }
    for(uint j=0;j<(params*params);++j){
        infoMatrix[ j ] = 0;
    }
    for(uint i=0;i<iObsSampleSize;++i){
        double betaX = 0.;
        uint col = i*stride;
        //cerr<<diseaseStatus[i];
        for(uint j=0;j<params;++j){
           //cerr<<" "<<effSizeMat[col+j];
            betaX+=betas[j] * effSizeMat[col+j];
        }
        //cerr<<endl;
        double pY = exp(betaX);
        pY/=(1+pY);
        if (diseaseStatus[i]==1) logL += log(pY); else logL += log(1-pY);
        double dev = diseaseStatus[i] - pY;
        double info0 = pY * (1-pY);
        for(uint j=0;j<params;++j){
            U[j]+=dev * effSizeMat[col+j];
            for(uint k=j;k<params;++k){
                infoMatrix[j*params+k]+=info0*effSizeMat[col+j]*effSizeMat[col+k];
                //cerr<<" "<<infoMatrix[j*params+k];
                if (k>j) infoMatrix[k*params+j] = infoMatrix[j*params+k];
            }
            //cerr<<endl;
        }
    }
    double *infoMatDecomp = new double[params*params];
    if (!chDecomp(infoMatrix,infoMatDecomp,params)){
        count = MAX;
        return;
    }
    invert(infoMatDecomp,invInfoMatrix,params);
    for(uint i=0;i<params;++i){
        uint col = i*params;
        for(uint j=0;j<params;++j){
            betas[i]+=U[j]*invInfoMatrix[col+j];
            chiSq+=1.0*U[i]*invInfoMatrix[col+j]*U[j];
        }
    }
}

bool fitModel(double & L1,double * phenovec_filtered,double * designmat,double * betas, double * var, int samplesize,int stride, int rank){
  int count=0;
  const int MAX=100;
  double chiSq;
  for (int i=0;i<rank;++i) betas[i] = 0;
  do{
    scoreTest(chiSq,L1,phenovec_filtered,designmat,betas,var,samplesize,stride,rank,count,MAX);
  }while(chiSq>.0001 && count++<MAX);
  
  if (count>=MAX) {
    return false;
  }else{
    return true;
  }
}


bool logisticReg( double & pvalue, double * phenovec_filtered, double * designmat_filtered, int & samplesize, int & rank, int & df){
     // run logistic regression for full and reduced models and compute p-value
     // input: phenovec_filtered, a double vector of zero/ones with no missing values
     // input: designmat_filtered, a double vector of designmatrix rows, in series, including intercept ones
     // input: rank, int variable, no of columns in design matrix
     // input: samplesize, int variable, no rows of design matrix
     // input: df, int variable, no of variables that differ between the reduced and saturated models. It is assumed that the test variables are the last columns in the design matrix.
        
     double L0,L1,dChiSq;
     double *betas = new double[rank];
     double *var = new double[rank*rank];
     bool fitted0, fitted1;
        
     fitted0 = fitModel(L0,phenovec_filtered,designmat_filtered,betas,var,
     samplesize,rank,rank-df);

     fitted1 = fitModel(L1,phenovec_filtered,designmat_filtered,betas,var,samplesize,rank,rank);
     if (fitted0 & fitted1){
       dChiSq = 2.*(L1-L0);
       pvalue = (1.-gsl_cdf_chisq_P(dChiSq,df));
       return true;
     }else{
      return false;
    }
}




