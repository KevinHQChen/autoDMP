/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_mpcmoveCodeGeneration_api.h
 *
 * Code generation for function 'mpcmoveCodeGeneration'
 *
 */

#ifndef _CODER_MPCMOVECODEGENERATION_API_H
#define _CODER_MPCMOVECODEGENERATION_API_H

/* Include files */
#include "emlrt.h"
#include "tmwtypes.h"
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct8_T
#define typedef_struct8_T
typedef struct {
  real_T ymin;
  real_T ymax;
  real_T umin;
  real_T umax;
} struct8_T;
#endif /* typedef_struct8_T */

#ifndef typedef_struct10_T
#define typedef_struct10_T
typedef struct {
  real_T Uopt[31];
  real_T Yopt[31];
  real_T Xopt[93];
  real_T Topt[31];
  real_T Slack;
  real_T Iterations;
  real_T Cost;
} struct10_T;
#endif /* typedef_struct10_T */

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  real_T Plant[3];
  real_T LastMove;
  real_T Covariance[9];
  boolean_T iA[62];
} struct4_T;
#endif /* typedef_struct4_T */

#ifndef typedef_struct6_T
#define typedef_struct6_T
typedef struct {
  real_T ym;
  real_T ref;
} struct6_T;
#endif /* typedef_struct6_T */

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  struct6_T signals;
  struct8_T limits;
} struct5_T;
#endif /* typedef_struct5_T */

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void mpcmoveCodeGeneration(struct4_T *statedata, struct5_T *onlinedata,
                           real_T *u, struct10_T *Info);

void mpcmoveCodeGeneration_api(const mxArray *const prhs[3], int32_T nlhs,
                               const mxArray *plhs[3]);

void mpcmoveCodeGeneration_atexit(void);

void mpcmoveCodeGeneration_initialize(void);

void mpcmoveCodeGeneration_terminate(void);

void mpcmoveCodeGeneration_xil_shutdown(void);

void mpcmoveCodeGeneration_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_mpcmoveCodeGeneration_api.h) */
