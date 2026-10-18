#ifndef LINEARFUNC_H
#define LINEARFUNC_H
bool linearRegCompare( double & pvalue, double * phenovec_filtered, double * designmat_filtered, int & samplesize, int & rank, int & df);
gsl_matrix* removeCollinearColumns(gsl_matrix *X, int &new_rank, int samplesize);
#endif

