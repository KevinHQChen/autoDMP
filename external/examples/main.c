/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "mpcmoveCodeGeneration.h"
#include "mpcmoveCodeGeneration_terminate.h"
#include "mpcmoveCodeGeneration_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void argInit_3x1_real_T(double result[3]);

static void argInit_3x3_real_T(double result[9]);

static void argInit_62x1_boolean_T(boolean_T result[62]);

static boolean_T argInit_boolean_T(void);

static double argInit_real_T(void);

static void argInit_struct4_T(struct4_T *result);

static struct5_T argInit_struct5_T(void);

static struct6_T argInit_struct6_T(void);

static struct8_T argInit_struct8_T(void);

/* Function Definitions */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

static void argInit_3x3_real_T(double result[9])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 9; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

static void argInit_62x1_boolean_T(boolean_T result[62])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 62; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_boolean_T();
  }
}

static boolean_T argInit_boolean_T(void)
{
  return false;
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void argInit_struct4_T(struct4_T *result)
{
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result->Plant);
  result->LastMove = argInit_real_T();
  argInit_3x3_real_T(result->Covariance);
  argInit_62x1_boolean_T(result->iA);
}

static struct5_T argInit_struct5_T(void)
{
  struct5_T result;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result.signals = argInit_struct6_T();
  result.limits = argInit_struct8_T();
  return result;
}

static struct6_T argInit_struct6_T(void)
{
  struct6_T result;
  double result_tmp;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result.ref = result_tmp;
  result.ym = result_tmp;
  return result;
}

static struct8_T argInit_struct8_T(void)
{
  struct8_T result;
  double result_tmp;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result.ymax = result_tmp;
  result.umin = result_tmp;
  result.umax = result_tmp;
  result.ymin = result_tmp;
  return result;
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_mpcmoveCodeGeneration();
  /* Terminate the application.
You do not need to do this more than one time. */
  mpcmoveCodeGeneration_terminate();
  return 0;
}

void main_mpcmoveCodeGeneration(void)
{
  struct10_T Info;
  struct4_T statedata;
  struct5_T r;
  double u;
  /* Initialize function 'mpcmoveCodeGeneration' input arguments. */
  /* Initialize function input argument 'statedata'. */
  /* Initialize function input argument 'onlinedata'. */
  /* Call the entry-point 'mpcmoveCodeGeneration'. */
  argInit_struct4_T(&statedata);
  r = argInit_struct5_T();
  mpcmoveCodeGeneration(&statedata, &r, &u, &Info);
}

/* End of code generation (main.c) */
