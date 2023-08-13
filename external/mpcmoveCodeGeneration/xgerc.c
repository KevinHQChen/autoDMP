/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgerc.c
 *
 * Code generation for function 'xgerc'
 *
 */

/* Include files */
#include "xgerc.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void xgerc(int m, int n, double alpha1, int ix0, const double y[2], double A[4],
           int ia0)
{
  int ijA;
  int j;
  if (!(alpha1 == 0.0)) {
    int i;
    int jA;
    jA = ia0;
    i = (unsigned char)n;
    for (j = 0; j < i; j++) {
      double temp;
      temp = y[j];
      if (temp != 0.0) {
        int i1;
        temp *= alpha1;
        i1 = m + jA;
        for (ijA = jA; ijA < i1; ijA++) {
          A[ijA - 1] += A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }
      jA += 2;
    }
  }
}

/* End of code generation (xgerc.c) */
