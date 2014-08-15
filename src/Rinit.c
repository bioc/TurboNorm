#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

void bincount(int *x, int *c, int *mm);
void binsum(int *x, double *y, double *s, int *mm);

static const R_CMethodDef cMethods[] = {
  {"bincount", (DL_FUNC) &bincount, 3},
  {"binsum", (DL_FUNC) &binsum, 4},
  NULL
};

void R_init_TurboNorm(DllInfo *info)
{
  R_registerRoutines(info, cMethods,  NULL, NULL, NULL);
};
