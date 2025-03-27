#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<list>
#include<vector>
#include<string.h>
#include<map>
#include<set>
#include<numeric>
#include<math.h>

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


// Check for collinearity and remove linearly dependent columns
gsl_matrix* removeCollinearColumns(gsl_matrix *X, int &new_rank, int samplesize) {
    int num_columns = new_rank;  // Initial number of columns

    // Allocate a permutation vector for column pivoting
    gsl_permutation *p = gsl_permutation_alloc(num_columns);

    // Allocate the tau vector for storing the Householder coefficients
    gsl_vector *tau = gsl_vector_alloc(std::min(samplesize, num_columns));

    // Perform QR decomposition with column pivoting
    gsl_matrix *QR = gsl_matrix_alloc(samplesize, num_columns);
    gsl_matrix_memcpy(QR, X);

    int signum;
    gsl_vector *norm = gsl_vector_alloc(num_columns);  // To store column norms

    // gsl_matrix *X_copy0 = gsl_matrix_alloc(samplesize, num_columns);      //   just for test
    // gsl_matrix_memcpy(X_copy0, X);

    gsl_linalg_QRPT_decomp(QR, tau, p, &signum, norm);  // QR decomposition with pivoting

    double tolerance = 1e-10;
    std::vector<int> keepColumns;

    // for (int i = 0; i < num_columns; i++) {
    //     if (fabs(gsl_matrix_get(QR, i, i)) >= tolerance) {
    //         keepColumns.push_back(gsl_permutation_get(p, i));
    //     }
    // }


    for (int i = 0; i < num_columns; i++) {
        // if(num_columns == 18){
        // std::cout<<gsl_permutation_get(p, i)<<"th value is: "<<gsl_matrix_get(QR, i, i)<<std::endl;
        //  }
        if (fabs(gsl_matrix_get(QR, i, i)) >= tolerance) {
            keepColumns.push_back(gsl_permutation_get(p, i));
        }
    }

    // Create a new matrix with only independent columns
    gsl_matrix *X_new = gsl_matrix_alloc(samplesize, keepColumns.size());
    for (size_t i = 0; i < keepColumns.size(); i++) {
        for (size_t j = 0; j < samplesize; j++) {
            gsl_matrix_set(X_new, j, i, gsl_matrix_get(X, j, keepColumns[i]));
        }
    }

    // Free allocated memory
    gsl_matrix_free(QR);
    gsl_vector_free(tau);
    gsl_vector_free(norm);
    gsl_permutation_free(p);


    new_rank = keepColumns.size();

    return X_new;
}

bool linearReg(double *phenovec_filtered, double *designmat_filtered, int &samplesize, int &rank, int &df, double *SR, int &actual_df) {
    gsl_matrix_view X_view = gsl_matrix_view_array(designmat_filtered, samplesize, rank);
    gsl_matrix_view subX_view = gsl_matrix_submatrix(&X_view.matrix, 0, 0, samplesize, rank - df);
    gsl_matrix* X = &(subX_view.matrix);

    int new_rank = rank - df;

    gsl_matrix *X_temp = gsl_matrix_alloc(samplesize, new_rank);
    gsl_matrix_memcpy(X_temp, X);  // Copy the original matrix X

    // Remove collinear columns
    gsl_matrix *X_work = removeCollinearColumns(X_temp, new_rank, samplesize);

    gsl_matrix *XT = gsl_matrix_alloc(new_rank, samplesize);
    gsl_matrix_transpose_memcpy(XT, X_work);

    gsl_matrix *XTX = gsl_matrix_alloc(new_rank, new_rank);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, XT, X_work, 0.0, XTX);

    gsl_matrix_view Y_view = gsl_matrix_view_array(phenovec_filtered, samplesize, 1);
    gsl_matrix* Y = &(Y_view.matrix);

    gsl_matrix *XTY = gsl_matrix_alloc(new_rank, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, XT, Y, 0.0, XTY);

    gsl_matrix *invert_XTX = gsl_matrix_alloc(new_rank, new_rank);
    gsl_matrix_memcpy(invert_XTX, XTX);
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    int returncode = gsl_linalg_cholesky_decomp1(invert_XTX);
    gsl_set_error_handler(old_handler);
    if (returncode != 0) {
        gsl_matrix_free(XT);
        gsl_matrix_free(X_temp);
        gsl_matrix_free(X_work);
        gsl_matrix_free(XTX);
        gsl_matrix_free(XTY);
        gsl_matrix_free(invert_XTX);
        return false;
    }
    gsl_linalg_cholesky_invert(invert_XTX);

    gsl_matrix *Beta = gsl_matrix_alloc(new_rank, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invert_XTX, XTY, 0.0, Beta);

    gsl_matrix *Xb = gsl_matrix_alloc(samplesize, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_work, Beta, 0.0, Xb);

    for (int i = 0; i < samplesize; i++) {
        SR[i] = gsl_matrix_get(Y, i, 0) - gsl_matrix_get(Xb, i, 0);
        SR[i] = SR[i] * SR[i];
    }

    gsl_matrix_free(XT);
    gsl_matrix_free(X_temp);
    gsl_matrix_free(X_work);
    gsl_matrix_free(XTX);
    gsl_matrix_free(XTY);
    gsl_matrix_free(invert_XTX);
    gsl_matrix_free(Beta);
    gsl_matrix_free(Xb);
    actual_df = new_rank;
    return true;
}

bool linearRegCompare(double &pvalue, double *phenovec_filtered, double *designmat_filtered, int &samplesize, int &rank, int &df) {
    double *SR_FULL = new double[samplesize];
    int to_reduce = 0;  // Define df as a variable
    int actual_df1 = 0;
    bool B1 = linearReg(phenovec_filtered, designmat_filtered, samplesize, rank, to_reduce, SR_FULL, actual_df1);

    double *SR_LIMIT = new double[samplesize];
    to_reduce = df;  // Use the appropriate value for df
    int actual_df2 = 0;
    bool B2 = linearReg(phenovec_filtered, designmat_filtered, samplesize, rank, to_reduce, SR_LIMIT, actual_df2);

    if (!B1 || !B2) {
        delete[] SR_FULL;
        delete[] SR_LIMIT;
        return false;
    }

    double RSS1=0.0, RSS2=0.0;
    for(int i=0;i<samplesize;i++){
        RSS1=RSS1+SR_LIMIT[i];
        RSS2=RSS2+SR_FULL[i];
    }


    double F_statistics = ((RSS1 - RSS2) / (actual_df1 - actual_df2)) / (RSS2 / (samplesize - actual_df1));
    pvalue = gsl_cdf_fdist_Q(F_statistics, static_cast<double>(df), static_cast<double>(samplesize - rank));

    delete[] SR_FULL;
    delete[] SR_LIMIT;
    return true;
}
