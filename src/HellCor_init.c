#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void hellcorC(double *x, int *xlen, double *statistic, int *pvalcomp,
		     double *pvalue, int *Kmaxvalue, int *Lmaxvalue,
		     int *Kset, int *Lset, double *alphavalue,
		     double *conflevel, int *B1, int *B2, double *CIetaleft, double *CIetaright);

extern void hilbertpeano(double *x, int *xlen, int *depth, int *setseed);

static const R_CMethodDef CEntries[] = {
    {"hellcorC", (DL_FUNC) &hellcorC, 15},
    {"hilbertpeano", (DL_FUNC) &hilbertpeano, 4},
    {NULL, NULL, 0}
};

void R_init_HellCor(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
