#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void hellcorC(double *x, int *xlen, double *statistic, int *pvalcomp,
		     double *pvalue, int *KLmaxvalue, double *alphavalue);


static const R_CMethodDef CEntries[] = {
    {"hellcorC", (DL_FUNC) &hellcorC, 8},
    {NULL, NULL, 0}
};

void R_init_HellCor(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
