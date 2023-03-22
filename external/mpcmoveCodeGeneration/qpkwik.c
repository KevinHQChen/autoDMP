/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * qpkwik.c
 *
 * Code generation for function 'qpkwik'
 *
 */

/* Include files */
#include "qpkwik.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include "xgerc.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void DropConstraint(int kDrop, boolean_T iA[62], int *nA, int iC[62]);

static double KWIKfactor(const double Ac[124], const int iC[62], int nA,
                         const double Linv[4], double RLinv[4], double D[4],
                         double H[4]);

static double rt_hypotd_snf(double u0, double u1);

static double rt_roundd_snf(double u);

/* Function Definitions */
static void DropConstraint(int kDrop, boolean_T iA[62], int *nA, int iC[62])
{
  int i;
  if (kDrop > 0) {
    int q0;
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      q0 = *nA;
      if (q0 < -2147483647) {
        q0 = MIN_int32_T;
      } else {
        q0--;
      }
      for (i = kDrop; i <= q0; i++) {
        iC[i - 1] = iC[i];
      }
    }
    iC[*nA - 1] = 0;
    q0 = *nA;
    if (q0 < -2147483647) {
      *nA = MIN_int32_T;
    } else {
      *nA = q0 - 1;
    }
  }
}

static double KWIKfactor(const double Ac[124], const int iC[62], int nA,
                         const double Linv[4], double RLinv[4], double D[4],
                         double H[4])
{
  double R[4];
  double TL[4];
  double tau[2];
  double work[2];
  double Q_idx_2;
  double Q_idx_3;
  double Status;
  double atmp;
  double beta1;
  int b_i;
  int exitg1;
  int i;
  int i1;
  int ia0;
  int iac;
  int ii;
  int knt;
  int lastc;
  int lastv;
  Status = 1.0;
  RLinv[0] = 0.0;
  RLinv[1] = 0.0;
  RLinv[2] = 0.0;
  RLinv[3] = 0.0;
  i = (unsigned char)nA;
  for (b_i = 0; b_i < i; b_i++) {
    i1 = iC[b_i];
    beta1 = Ac[i1 - 1];
    atmp = Ac[i1 + 61];
    i1 = b_i << 1;
    RLinv[i1] = Linv[0] * beta1 + Linv[2] * atmp;
    RLinv[i1 + 1] = Linv[1] * beta1 + Linv[3] * atmp;
  }
  TL[0] = RLinv[0];
  TL[1] = RLinv[1];
  TL[2] = RLinv[2];
  TL[3] = RLinv[3];
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  for (b_i = 0; b_i < 2; b_i++) {
    ii = (b_i << 1) + b_i;
    if (b_i + 1 < 2) {
      atmp = TL[ii];
      lastc = ii + 2;
      tau[0] = 0.0;
      beta1 = fabs(TL[ii + 1]);
      if (beta1 != 0.0) {
        beta1 = rt_hypotd_snf(atmp, beta1);
        if (atmp >= 0.0) {
          beta1 = -beta1;
        }
        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          do {
            knt++;
            for (lastv = lastc; lastv <= lastc; lastv++) {
              TL[lastv - 1] *= 9.9792015476736E+291;
            }
            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt < 20));
          beta1 = rt_hypotd_snf(atmp, fabs(TL[ii + 1]));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }
          tau[0] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          for (lastv = lastc; lastv <= lastc; lastv++) {
            TL[lastv - 1] *= atmp;
          }
          for (lastv = 0; lastv < knt; lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }
          atmp = beta1;
        } else {
          tau[0] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          for (lastv = lastc; lastv <= lastc; lastv++) {
            TL[lastv - 1] *= atmp;
          }
          atmp = beta1;
        }
      }
      TL[ii] = 1.0;
      if (tau[0] != 0.0) {
        lastv = 2;
        lastc = ii + 1;
        while ((lastv > 0) && (TL[lastc] == 0.0)) {
          lastv--;
          lastc--;
        }
        lastc = 1;
        knt = ii + 2;
        do {
          exitg1 = 0;
          if (knt + 1 <= (ii + lastv) + 2) {
            if (TL[knt] != 0.0) {
              exitg1 = 1;
            } else {
              knt++;
            }
          } else {
            lastc = 0;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      } else {
        lastv = 0;
        lastc = 0;
      }
      if (lastv > 0) {
        ia0 = ii + 3;
        if (lastc != 0) {
          work[0] = 0.0;
          for (iac = ia0; iac <= ia0; iac += 2) {
            beta1 = 0.0;
            i1 = (iac + lastv) - 1;
            for (knt = iac; knt <= i1; knt++) {
              beta1 += TL[knt - 1] * TL[(ii + knt) - iac];
            }
            knt = ((iac - ii) - 3) >> 1;
            work[knt] += beta1;
          }
        }
        xgerc(lastv, lastc, -tau[0], ii + 1, work, TL, ii + 3);
      }
      TL[ii] = atmp;
    } else {
      tau[1] = 0.0;
    }
  }
  for (ia0 = 0; ia0 < 2; ia0++) {
    for (b_i = 0; b_i <= ia0; b_i++) {
      lastc = b_i + (ia0 << 1);
      R[lastc] = TL[lastc];
    }
    if (ia0 + 2 <= 2) {
      R[(ia0 << 1) + 1] = 0.0;
    }
    work[ia0] = 0.0;
  }
  for (b_i = 1; b_i >= 0; b_i--) {
    ii = b_i + (b_i << 1);
    if (b_i + 1 < 2) {
      TL[ii] = 1.0;
      if (tau[b_i] != 0.0) {
        lastv = 2;
        lastc = ii;
        while ((lastv > 0) && (TL[lastc + 1] == 0.0)) {
          lastv--;
          lastc--;
        }
        lastc = 1;
        knt = ii + 2;
        do {
          exitg1 = 0;
          if (knt + 1 <= (ii + lastv) + 2) {
            if (TL[knt] != 0.0) {
              exitg1 = 1;
            } else {
              knt++;
            }
          } else {
            lastc = 0;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      } else {
        lastv = 0;
        lastc = 0;
      }
      if (lastv > 0) {
        ia0 = ii + 3;
        if (lastc != 0) {
          work[0] = 0.0;
          for (iac = ia0; iac <= ia0; iac += 2) {
            beta1 = 0.0;
            i1 = (iac + lastv) - 1;
            for (knt = iac; knt <= i1; knt++) {
              beta1 += TL[knt - 1] * TL[(ii + knt) - iac];
            }
            knt = ((iac - ii) - 3) >> 1;
            work[knt] += beta1;
          }
        }
        xgerc(lastv, lastc, -tau[b_i], ii + 1, work, TL, ii + 3);
      }
      lastc = ii + 2;
      for (lastv = lastc; lastv <= lastc; lastv++) {
        TL[lastv - 1] *= -tau[b_i];
      }
    }
    TL[ii] = 1.0 - tau[b_i];
    if (b_i - 1 >= 0) {
      TL[ii - 1] = 0.0;
    }
  }
  beta1 = TL[0];
  atmp = TL[1];
  Q_idx_2 = TL[2];
  Q_idx_3 = TL[3];
  b_i = 0;
  do {
    exitg1 = 0;
    if (b_i <= (unsigned char)nA - 1) {
      if (fabs(R[b_i + (b_i << 1)]) < 1.0E-12) {
        Status = -2.0;
        exitg1 = 1;
      } else {
        b_i++;
      }
    } else {
      for (b_i = 0; b_i < 2; b_i++) {
        double Linv_tmp;
        double b_Linv_tmp;
        lastc = b_i << 1;
        Linv_tmp = Linv[lastc];
        b_Linv_tmp = Linv[lastc + 1];
        TL[b_i] = Linv_tmp * beta1 + b_Linv_tmp * atmp;
        TL[b_i + 2] = Linv_tmp * Q_idx_2 + b_Linv_tmp * Q_idx_3;
      }
      RLinv[0] = 0.0;
      RLinv[1] = 0.0;
      RLinv[2] = 0.0;
      RLinv[3] = 0.0;
      for (ia0 = nA; ia0 >= 1; ia0--) {
        i1 = (ia0 + ((ia0 - 1) << 1)) - 1;
        RLinv[i1] = 1.0;
        for (lastv = ia0; lastv <= nA; lastv++) {
          knt = (ia0 + ((lastv - 1) << 1)) - 1;
          RLinv[knt] /= R[i1];
        }
        if (ia0 > 1) {
          for (lastv = 2; lastv <= nA; lastv++) {
            RLinv[2] -= R[2] * RLinv[3];
          }
        }
      }
      if (nA > 2147483646) {
        knt = MAX_int32_T;
      } else {
        knt = nA + 1;
      }
      for (b_i = 0; b_i < 2; b_i++) {
        for (ia0 = b_i + 1; ia0 < 3; ia0++) {
          i1 = b_i + ((ia0 - 1) << 1);
          H[i1] = 0.0;
          for (lastv = knt; lastv < 3; lastv++) {
            H[i1] -= TL[b_i + 2] * TL[ia0 + 1];
          }
          H[(ia0 + (b_i << 1)) - 1] = H[i1];
        }
      }
      for (ia0 = 0; ia0 < i; ia0++) {
        for (b_i = 0; b_i < 2; b_i++) {
          i1 = b_i + (ia0 << 1);
          D[i1] = 0.0;
          for (lastv = ia0 + 1; lastv <= nA; lastv++) {
            knt = (lastv - 1) << 1;
            D[i1] += TL[b_i + knt] * RLinv[ia0 + knt];
          }
        }
      }
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  return Status;
}

static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }
  return y;
}

void qpkwik(const double Linv[4], const double Hinv[4], const double f[2],
            const double Ac[124], const double b[62], boolean_T iA[62],
            double x[2], double lambda[62], int *status)
{
  double cTol[62];
  double D[4];
  double H[4];
  double RLinv[4];
  double Rhs[4];
  double U[4];
  double r[2];
  double lambdamin;
  double rMin;
  double zTa;
  double z_idx_0;
  double z_idx_1;
  int iC[62];
  int b_i;
  int i;
  int iSave;
  int k;
  int kDrop;
  int kNext;
  int nA;
  boolean_T ColdReset;
  boolean_T cTolComputed;
  boolean_T guard1 = false;
  x[0] = 0.0;
  x[1] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 62; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }
  nA = 0;
  for (i = 0; i < 62; i++) {
    if (iA[i]) {
      nA++;
      iC[nA - 1] = i + 1;
    }
  }
  guard1 = false;
  if (nA > 0) {
    double Opt[4];
    int exitg3;
    int tmp;
    Opt[0] = 0.0;
    Opt[1] = 0.0;
    Opt[2] = 0.0;
    Opt[3] = 0.0;
    Rhs[0] = f[0];
    Rhs[2] = 0.0;
    Rhs[1] = f[1];
    Rhs[3] = 0.0;
    tmp = (int)rt_roundd_snf(0.3 * (double)nA);
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((nA > 0) && (*status <= 256)) {
        lambdamin = KWIKfactor(Ac, iC, nA, Linv, RLinv, D, H);
        if (lambdamin < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 1;
          } else {
            nA = 0;
            memset(&iC[0], 0, 62U * sizeof(int));
            for (i = 0; i < 62; i++) {
              iA[i] = false;
            }
            ColdReset = true;
          }
        } else {
          b_i = (unsigned char)nA;
          for (kNext = 0; kNext < b_i; kNext++) {
            Rhs[kNext + 2] = b[iC[kNext] - 1];
            for (i = kNext + 1; i <= nA; i++) {
              iSave = (i + (kNext << 1)) - 1;
              U[iSave] = 0.0;
              for (k = 0; k < b_i; k++) {
                kDrop = k << 1;
                U[iSave] += RLinv[(i + kDrop) - 1] * RLinv[kNext + kDrop];
              }
              U[kNext + ((i - 1) << 1)] = U[iSave];
            }
          }
          for (i = 0; i < 2; i++) {
            Opt[i] = H[i] * Rhs[0] + H[i + 2] * Rhs[1];
            for (k = 0; k < b_i; k++) {
              Opt[i] += D[i + (k << 1)] * Rhs[k + 2];
            }
          }
          for (i = 0; i < b_i; i++) {
            iSave = i << 1;
            Opt[i + 2] = D[iSave] * Rhs[0] + D[iSave + 1] * Rhs[1];
            for (k = 0; k < b_i; k++) {
              Opt[i + 2] += U[i + (k << 1)] * Rhs[k + 2];
            }
          }
          lambdamin = -1.0E-12;
          kDrop = 0;
          for (i = 0; i < b_i; i++) {
            zTa = Opt[i + 2];
            lambda[iC[i] - 1] = zTa;
            if ((zTa < lambdamin) && (i + 1 <= nA)) {
              kDrop = i + 1;
              lambdamin = zTa;
            }
          }
          if (kDrop <= 0) {
            x[0] = Opt[0];
            x[1] = Opt[1];
            exitg3 = 2;
          } else {
            (*status)++;
            if (tmp <= 5) {
              b_i = 5;
            } else {
              b_i = tmp;
            }
            if (*status > b_i) {
              nA = 0;
              memset(&iC[0], 0, 62U * sizeof(int));
              for (i = 0; i < 62; i++) {
                iA[i] = false;
              }
              ColdReset = true;
            } else {
              lambda[iC[kDrop - 1] - 1] = 0.0;
              DropConstraint(kDrop, iA, &nA, iC);
            }
          }
        }
      } else {
        exitg3 = 2;
      }
    } while (exitg3 == 0);
    if (exitg3 != 1) {
      if (nA <= 0) {
        memset(&lambda[0], 0, 62U * sizeof(double));
        x[0] = -Hinv[0] * f[0] + -Hinv[2] * f[1];
        x[1] = -Hinv[1] * f[0] + -Hinv[3] * f[1];
      }
      guard1 = true;
    }
  } else {
    x[0] = -Hinv[0] * f[0] + -Hinv[2] * f[1];
    x[1] = -Hinv[1] * f[0] + -Hinv[3] * f[1];
    guard1 = true;
  }
  if (guard1) {
    double Xnorm0;
    boolean_T exitg2;
    Xnorm0 = b_norm(x);
    exitg2 = false;
    while ((!exitg2) && (*status <= 256)) {
      double cMin;
      cMin = -1.0E-6;
      kNext = -1;
      for (i = 0; i < 62; i++) {
        if (!cTolComputed) {
          lambdamin = fabs(Ac[i] * x[0]);
          zTa = fabs(Ac[i + 62] * x[1]);
          if ((lambdamin < zTa) || (rtIsNaN(lambdamin) && (!rtIsNaN(zTa)))) {
            lambdamin = zTa;
          }
          cTol[i] = fmax(cTol[i], lambdamin);
        }
        if (!iA[i]) {
          lambdamin = ((Ac[i] * x[0] + Ac[i + 62] * x[1]) - b[i]) / cTol[i];
          if (lambdamin < cMin) {
            cMin = lambdamin;
            kNext = i;
          }
        }
      }
      cTolComputed = true;
      if (kNext + 1 <= 0) {
        exitg2 = true;
      } else if (*status == 256) {
        *status = 0;
        exitg2 = true;
      } else {
        int exitg1;
        do {
          exitg1 = 0;
          if ((kNext + 1 > 0) && (*status <= 256)) {
            boolean_T guard2 = false;
            guard2 = false;
            if (nA == 0) {
              lambdamin = Ac[kNext + 62];
              z_idx_0 = Hinv[0] * Ac[kNext] + Hinv[2] * lambdamin;
              z_idx_1 = Hinv[1] * Ac[kNext] + Hinv[3] * lambdamin;
              guard2 = true;
            } else {
              lambdamin = KWIKfactor(Ac, iC, nA, Linv, RLinv, D, H);
              if (lambdamin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                lambdamin = Ac[kNext + 62];
                z_idx_0 = -H[0] * Ac[kNext] + -H[2] * lambdamin;
                z_idx_1 = -H[1] * Ac[kNext] + -H[3] * lambdamin;
                b_i = (unsigned char)nA;
                for (i = 0; i < b_i; i++) {
                  iSave = i << 1;
                  r[i] = Ac[kNext] * D[iSave] + lambdamin * D[iSave + 1];
                }
                guard2 = true;
              }
            }
            if (guard2) {
              double t1;
              boolean_T exitg4;
              boolean_T isT1Inf;
              kDrop = 0;
              t1 = 0.0;
              isT1Inf = true;
              ColdReset = true;
              if (nA > 0) {
                iSave = 0;
                exitg4 = false;
                while ((!exitg4) && (iSave <= (unsigned char)nA - 1)) {
                  if (r[iSave] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    iSave++;
                  }
                }
              }
              if ((nA != 0) && (!ColdReset)) {
                b_i = (unsigned char)nA;
                for (i = 0; i < b_i; i++) {
                  lambdamin = r[i];
                  if (lambdamin > 1.0E-12) {
                    lambdamin = lambda[iC[i] - 1] / lambdamin;
                    if ((kDrop == 0) || (lambdamin < rMin)) {
                      rMin = lambdamin;
                      kDrop = i + 1;
                    }
                  }
                }
                if (kDrop > 0) {
                  t1 = rMin;
                  isT1Inf = false;
                }
              }
              lambdamin = Ac[kNext + 62];
              zTa = z_idx_0 * Ac[kNext] + z_idx_1 * lambdamin;
              if (zTa <= 0.0) {
                lambdamin = 0.0;
                ColdReset = true;
              } else {
                lambdamin =
                    (b[kNext] - (Ac[kNext] * x[0] + lambdamin * x[1])) / zTa;
                ColdReset = false;
              }
              if (isT1Inf && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  cMin = t1;
                } else if (isT1Inf) {
                  cMin = lambdamin;
                } else if (t1 < lambdamin) {
                  cMin = t1;
                } else {
                  cMin = lambdamin;
                }
                b_i = (unsigned char)nA;
                for (i = 0; i < b_i; i++) {
                  iSave = iC[i];
                  lambda[iSave - 1] -= cMin * r[i];
                  if (lambda[iSave - 1] < 0.0) {
                    lambda[iSave - 1] = 0.0;
                  }
                }
                lambda[kNext] += cMin;
                if (fabs(cMin - t1) < 2.2204460492503131E-16) {
                  DropConstraint(kDrop, iA, &nA, iC);
                }
                if (!ColdReset) {
                  x[0] += cMin * z_idx_0;
                  x[1] += cMin * z_idx_1;
                  if (fabs(cMin - lambdamin) < 2.2204460492503131E-16) {
                    if (nA == 2) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      if (nA > 2147483646) {
                        nA = MAX_int32_T;
                      } else {
                        nA++;
                      }
                      iC[nA - 1] = kNext + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while ((!exitg4) && (i + 1 > 1)) {
                        b_i = iC[i - 1];
                        if (iC[i] > b_i) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = b_i;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }
                      iA[kNext] = true;
                      kNext = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            lambdamin = b_norm(x);
            if (fabs(lambdamin - Xnorm0) > 0.001) {
              Xnorm0 = lambdamin;
              for (k = 0; k < 62; k++) {
                cTol[k] = fmax(fabs(b[k]), 1.0);
              }
              cTolComputed = false;
            }
            exitg1 = 2;
          }
        } while (exitg1 == 0);
        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

/* End of code generation (qpkwik.c) */
