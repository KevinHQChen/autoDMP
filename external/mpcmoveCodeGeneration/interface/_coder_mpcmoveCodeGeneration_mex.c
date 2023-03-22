/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_mpcmoveCodeGeneration_mex.c
 *
 * Code generation for function 'mpcmoveCodeGeneration'
 *
 */

/* Include files */
#include "_coder_mpcmoveCodeGeneration_mex.h"
#include "_coder_mpcmoveCodeGeneration_api.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&mpcmoveCodeGeneration_atexit);
  /* Module initialization. */
  mpcmoveCodeGeneration_initialize();
  /* Dispatch the entry-point. */
  unsafe_mpcmoveCodeGeneration_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  mpcmoveCodeGeneration_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

void unsafe_mpcmoveCodeGeneration_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                              int32_T nrhs,
                                              const mxArray *prhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[3];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs < 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooFewInputsConstants", 9, 4, 21,
                        "mpcmoveCodeGeneration", 4, 21, "mpcmoveCodeGeneration",
                        4, 21, "mpcmoveCodeGeneration");
  }
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4,
                        21, "mpcmoveCodeGeneration");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 21,
                        "mpcmoveCodeGeneration");
  }
  /* Call the function. */
  mpcmoveCodeGeneration_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_mpcmoveCodeGeneration_mex.c) */
