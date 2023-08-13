/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * qpkwik.h
 *
 * Code generation for function 'qpkwik'
 *
 */

#ifndef QPKWIK_H
#define QPKWIK_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void qpkwik(const double Linv[4], const double Hinv[4], const double f[2],
            const double Ac[124], const double b[62], boolean_T iA[62],
            double x[2], double lambda[62], int *status);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (qpkwik.h) */
