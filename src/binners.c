#include <Rinternals.h>
/* Fast counts and sums in bins */

void bincount(int *x, int *c, int *mm) {
  int i, m;
  m = *mm - 1;
  for (i = 0; i < m; i++)  c[x[i]] += 1;
}

void binsum(int *x, double *y, double *s, int *mm) {
  int i, m;
  m = *mm - 1;
  for (i = 0; i < m; i++)  s[x[i]] += y[i];
}


