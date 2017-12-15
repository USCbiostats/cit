#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

extern void citconlog2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void citconlog2cvr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void citconlog3p(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void citconlog3pcvr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"citconlog2",     (DL_FUNC) &citconlog2,     11},
    {"citconlog2cvr",  (DL_FUNC) &citconlog2cvr,  13},
    {"citconlog3p",    (DL_FUNC) &citconlog3p,    13},
    {"citconlog3pcvr", (DL_FUNC) &citconlog3pcvr, 15},
    {NULL, NULL, 0}
};

void R_init_cit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
