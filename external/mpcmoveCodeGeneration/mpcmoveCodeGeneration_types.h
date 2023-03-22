/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mpcmoveCodeGeneration_types.h
 *
 * Code generation for function 'mpcmoveCodeGeneration'
 *
 */

#ifndef MPCMOVECODEGENERATION_TYPES_H
#define MPCMOVECODEGENERATION_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct8_T
#define typedef_struct8_T
typedef struct {
  double ymin;
  double ymax;
  double umin;
  double umax;
} struct8_T;
#endif /* typedef_struct8_T */

#ifndef typedef_struct10_T
#define typedef_struct10_T
typedef struct {
  double Uopt[31];
  double Yopt[31];
  double Xopt[93];
  double Topt[31];
  double Slack;
  double Iterations;
  double Cost;
} struct10_T;
#endif /* typedef_struct10_T */

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  double Plant[3];
  double LastMove;
  double Covariance[9];
  boolean_T iA[62];
} struct4_T;
#endif /* typedef_struct4_T */

#ifndef typedef_struct6_T
#define typedef_struct6_T
typedef struct {
  double ym;
  double ref;
} struct6_T;
#endif /* typedef_struct6_T */

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  struct6_T signals;
  struct8_T limits;
} struct5_T;
#endif /* typedef_struct5_T */

#endif
/* End of code generation (mpcmoveCodeGeneration_types.h) */
