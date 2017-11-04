typedef short unsigned int suint;
typedef unsigned int uint;

bool chDecomp(const double * inputMatrix, double * decompMat,uint dim);

void invert(const double * decompMat,double * invertedMat,uint dim);

void scoreTest(double & chiSq,double & logL,double * diseaseStatus,
double * effSizeMat,double * betas,double * invInfoMatrix,uint iObsSampleSize,
suint stride,suint params,int & count,int MAX);

bool fitModel(double & L1,double * phenovec_filtered,double * designmat,double * betas, double * var, int samplesize,int stride, int rank);

bool logisticReg( double & pvalue, double * phenovec_filtered, double * designmat_filtered, int & samplesize, int & rank, int & df);



