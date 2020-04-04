/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include "rtutil_math.h"
#include "rtutil_scat_int.h"



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_integ_gc_init(int n_coef, int n_quad,
                              double *qx, double *qw, double ****P_, int symmetric) {

     int i;

     int ii2;

     double ***P;

     if (! symmetric) {
          P = alloc_array3_d(4, n_quad, n_coef);

          for (i = 0; i < n_quad; ++i)
               gen_spher_funcs_prt(n_coef, qx[i], P[0][i], P[1][i], P[2][i], P[3][i]);
     }
     else {
          P = alloc_array3_d(4, (n_quad + 1) / 2, n_coef);

          if (n_quad % 2 != 0) {
               i = n_quad / 2;

               ii2 = i + 1;

               gen_spher_funcs_prt(n_coef, qx[i], P[0][i], P[1][i], P[2][i], P[3][i]);
          }
          else
               ii2 = n_quad / 2;

          for (i = 0; i < n_quad / 2; ++i) {
               gen_spher_funcs_prt(n_coef, qx[ii2], P[0][i], P[1][i], P[2][i], P[3][i]);

               ii2++;
          }
     }

     *P_ = P;
}



/*******************************************************************************
 *
 ******************************************************************************/
int phase_func_integ_gc(int n_coef, int n_quad, double **pf,
                        double *qx, double *qw, double ***P_init,
                        double **coefs, double accuracy, int allsix, int symmetric) {

     int i;
     int ii1;
     int ii2;
     int j;
     int k;

     int n_elem;

     double a;

     double b1[6];
     double b2[6];

     double *P[4];

     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;

     if (! P_init) {
          P[0] = alloc_array1_d(n_coef);
          P[1] = alloc_array1_d(n_coef);
          P[2] = alloc_array1_d(n_coef);
          P[3] = alloc_array1_d(n_coef);
     }

     for (i = 0; i < n_elem; ++i) {
          for (j = 0; j < n_coef; ++j) {
               coefs[i][j] = 0.;
          }
     }

     if (! symmetric) {
          for (i = 0; i < n_quad; ++i) {
               if (! P_init)
                    gen_spher_funcs_prt(n_coef, qx[i], P[0], P[1], P[2], P[3]);
               else {
                    P[0] = P_init[0][i];
                    P[1] = P_init[1][i];
                    P[2] = P_init[2][i];
                    P[3] = P_init[3][i];
               }

               for (k = 0; k < n_elem; ++k)
                    b1[k] = pf[k][i] * qw[i];
               a = b1[1];
               b1[1] = a + b1[2];
               b1[2] = a - b1[2];

               for (j = 0; j < n_coef; ++j) {
                    coefs[0][j] += P[0][j] * b1[0];
                    coefs[1][j] += P[2][j] * b1[1];
                    coefs[2][j] += P[3][j] * b1[2];
                    coefs[3][j] += P[0][j] * b1[3];
                    coefs[4][j] += P[1][j] * b1[4];
                    coefs[5][j] += P[1][j] * b1[5];
               }
          }
     }
     else {
          if (n_quad % 2 != 0) {
               i = n_quad / 2;

               ii1 = i - 1;
               ii2 = i + 1;

               if (! P_init)
                    gen_spher_funcs_prt(n_coef, qx[i], P[0], P[1], P[2], P[3]);
               else {
                    P[0] = P_init[0][i];
                    P[1] = P_init[1][i];
                    P[2] = P_init[2][i];
                    P[3] = P_init[3][i];
               }

               for (k = 0; k < n_elem; ++k)
                    b1[k] = pf[k][i] * qw[i];
               a = b1[1];
               b1[1] = a + b1[2];
               b1[2] = a - b1[2];

               for (j = 0; j < n_coef; ++j) {
                    coefs[0][j] += P[0][j] * b1[0];
                    coefs[1][j] += P[2][j] * b1[1];
                    coefs[2][j] += P[3][j] * b1[2];
                    coefs[3][j] += P[0][j] * b1[3];
                    coefs[4][j] += P[1][j] * b1[4];
                    coefs[5][j] += P[1][j] * b1[5];
               }
          }
          else {
               ii1 = n_quad / 2 - 1;
               ii2 = n_quad / 2;
          }

          for (i = 0; i < n_quad / 2; ++i) {
               if (! P_init)
                    gen_spher_funcs_prt(n_coef, qx[ii2], P[0], P[1], P[2], P[3]);
               else {
                    P[0] = P_init[0][i];
                    P[1] = P_init[1][i];
                    P[2] = P_init[2][i];
                    P[3] = P_init[3][i];
               }

               for (k = 0; k < n_elem; ++k)
                    b1[k] = pf[k][ii1] * qw[ii2];
               a = b1[1];
               b1[1] = a + b1[2];
               b1[2] = a - b1[2];

               for (k = 0; k < n_elem; ++k)
                    b2[k] = pf[k][ii2] * qw[ii2];
               a = b2[1];
               b2[1] = a + b2[2];
               b2[2] = a - b2[2];

               a = 1.;
               for (j = 0; j < n_coef; ++j) {
                    coefs[0][j] += P[0][j] * (b1[0] * a + b2[0]);
                    coefs[1][j] += P[3][j] *  b1[1] * a + P[2][j] * b2[1];
                    coefs[2][j] += P[2][j] *  b1[2] * a + P[3][j] * b2[2];
                    coefs[3][j] += P[0][j] * (b1[3] * a + b2[3]);
                    coefs[4][j] += P[1][j] * (b1[4] * a + b2[4]);
                    coefs[5][j] += P[1][j] * (b1[5] * a + b2[5]);

                    a *= -1.;
               }

               ii1--;
               ii2++;
          }
     }

     for (j = 0; j < n_coef; ++j) {
          a = coefs[1][j];
          coefs[1][j] = (a + coefs[2][j]) / 2.;
          coefs[2][j] = (a - coefs[2][j]) / 2.;
     }

     for (i = 0; i < n_coef; ++i) {
          a = (2 * i + 1) / 2.;
          for (j = 0; j < n_elem; ++j) {
               coefs[j][i] *= a;
          }
     }

     for (j = 0; j < n_coef; ++j) {
          if (fabs(coefs[0][j]) <= accuracy) {
               j++;
               break;
          }
     }

     if (! P_init) {
          free_array1_d(P[0]);
          free_array1_d(P[1]);
          free_array1_d(P[2]);
          free_array1_d(P[3]);
     }

     return j;
}



int phase_func_integ_gc2(int n_coef, int n_quad, double **pf,
                         double *qx, double *qw, double ***P,
                         double **coefs, double accuracy, int allsix, int symmetric) {

     int i;
     int j;

     int jj1;
     int jj2;

     int n_elem;

     int n_coef2;

     double a;
     double b;

     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;

     for (i = 0; i < n_elem; ++i) {
          for (j = 0; j < n_coef; ++j) {
               coefs[i][j] = 0.;
          }
     }

     if (! symmetric) {
          for (i = 0; i < n_coef; ++i) {
               for (j = 0; j < n_quad; ++j) {
                    coefs[0][i] += P[0][j][i] *  pf[0][j] * qw[j];
                    coefs[1][i] += P[2][j][i] * (pf[1][j] + pf[2][j]) * qw[j];
                    coefs[2][i] += P[3][j][i] * (pf[1][j] - pf[2][j]) * qw[j];
                    coefs[3][i] += P[0][j][i] *  pf[3][j] * qw[j];
                    coefs[4][i] += P[1][j][i] *  pf[4][j] * qw[j];
                    coefs[5][i] += P[1][j][i] *  pf[5][j] * qw[j];
               }

               b = (2 * i + 1) / 2.;
               for (j = 0; j < n_elem; ++j)
                    coefs[j][i] *= b;

               if (fabs(coefs[0][i]) <= accuracy) {
                    i++; break;
               }
          }
     }
     else {
          if (n_quad % 2 != 0) {
               i = n_quad / 2;

               for (j = 0; j < n_coef; ++j) {

                    coefs[0][j] += P[0][i][j] *  pf[0][i] * qw[i];
                    coefs[1][j] += P[2][i][j] * (pf[1][i] + pf[2][i]) * qw[i];
                    coefs[2][j] += P[3][i][j] * (pf[1][i] - pf[2][i]) * qw[i];
                    coefs[3][j] += P[0][i][j] *  pf[3][i] * qw[i];
                    coefs[4][j] += P[1][i][j] *  pf[4][i] * qw[i];
                    coefs[5][j] += P[1][i][j] *  pf[5][i] * qw[i];
               }
          }

          a = 1.;
          for (i = 0; i < n_coef; ++i) {
               if (n_quad % 2 != 0) {
                    j = n_quad / 2;

                    jj1 = j - 1;
                    jj2 = j + 1;
               }
               else {
                    jj1 = n_quad / 2 - 1;
                    jj2 = n_quad / 2;
               }

               for (j = 0; j < n_quad / 2; ++j) {

                    coefs[0][i] +=  P[0][j][i] * (pf[0][jj1] * a + pf[0][jj2])  * qw[jj2];
                    coefs[1][i] += (P[3][j][i] * (pf[1][jj1]     + pf[2][jj1])  * a +
                                    P[2][j][i] * (pf[1][jj2]     + pf[2][jj2])) * qw[jj2];
                    coefs[2][i] += (P[2][j][i] * (pf[1][jj1]     - pf[2][jj1])  * a +
                                    P[3][j][i] * (pf[1][jj2]     - pf[2][jj2])) * qw[jj2];
                    coefs[3][i] +=  P[0][j][i] * (pf[3][jj1] * a + pf[3][jj2])  * qw[jj2];
                    coefs[4][i] +=  P[1][j][i] * (pf[4][jj1] * a + pf[4][jj2])  * qw[jj2];
                    coefs[5][i] +=  P[1][j][i] * (pf[5][jj1] * a + pf[5][jj2])  * qw[jj2];

                    jj1--;
                    jj2++;
               }

               b = (2 * i + 1) / 2.;
               for (j = 0; j < n_elem; ++j)
                    coefs[j][i] *= b;

               if (fabs(coefs[0][i]) <= accuracy) {
                    i++; break;
               }

               a *= -1.;
          }
     }

     n_coef2 = i;

     for (i = 0; i < n_coef2; ++i) {
          a = coefs[1][i];
          coefs[1][i] = (a + coefs[2][i]) / 2.;
          coefs[2][i] = (a - coefs[2][i]) / 2.;
     }

     return n_coef2;
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_integ_lc_init(int n_coef, int n_quad,
                              double *qx, double *qw, double **P_p) {

     int i;
     int ii2;

     if (n_quad % 2 != 0) {
          i = n_quad / 2;

          ii2 = i + 1;

          leg_poly(n_coef, qx[i], P_p[i]);
     }
     else
          ii2 = n_quad / 2;

     for (i = 0; i < n_quad / 2; ++i) {
          leg_poly(n_coef, qx[ii2], P_p[i]);

          ii2++;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int phase_func_integ_lc(int n_coef, int n_quad, double **pf,
                        double *qx, double *qw, double **P_p,
                        double **coefs, double accuracy, int allsix) {

     int i;
     int ii1;
     int ii2;
     int j;
     int k;

     int n_elem;

     double a;

     double b1[6];
     double b2[6];
/*
     double *P_p;
*/
     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;
/*
     P_p = alloc_array1_d(n_coef);
*/
     for (i = 0; i < n_elem; ++i) {
          for (j = 0; j < n_coef; ++j) {
               coefs[i][j] = 0.;
          }
     }

     if (n_quad % 2 != 0) {
          i = n_quad / 2;

          ii1 = i - 1;
          ii2 = i + 1;
/*
          calc_leg_p(n_coef, qx[i], P_p);
*/
          for (k = 0; k < n_elem; ++k) {
               b1[k] = pf[k][i] * qw[i];
          }

          for (j = 0; j < n_coef; ++j) {
               for (k = 0; k < n_elem; ++k) {
                    coefs[k][j] += P_p[i][j] * b1[k];
               }
          }
     }
     else {
          ii1 = n_quad / 2 - 1;
          ii2 = n_quad / 2;
     }

     for (i = 0; i < n_quad / 2; ++i) {
/*
          calc_leg_p(n_coef, qx[ii2], P_p);
*/
          for (k = 0; k < n_elem; ++k) {
               b1[k] = pf[k][ii1] * qw[ii2];
               b2[k] = pf[k][ii2] * qw[ii2];
          }

          a = 1.;
          for (j = 0; j < n_coef; ++j) {
               for (k = 0; k < n_elem; ++k) {
                    coefs[k][j] += P_p[i][j] * (b1[k] * a + b2[k]);
               }

               a *= -1.;
          }

          ii1--;
          ii2++;
     }

     if (allsix) {
          for (i = 0; i < n_coef; ++i) {
               coefs[5][i] = -coefs[5][i];
          }
     }

     for (i = 0; i < n_coef; ++i) {
          a = (2. * i + 1.) / 2.;
          for (j = 0; j < n_elem; ++j) {
               coefs[j][i] *= a;
          }
     }

     for (j = 0; j < n_coef; ++j) {
          if (fabs(coefs[0][j]) <= accuracy) {
               j++; break;
          }
     }
/*
     free_array1_d(P_p);
*/
     return j;
}



/*******************************************************************************
 *
 ******************************************************************************/
double phase_func_d_theta(int n_theta) {

     return 180. / (n_theta - 1);
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_theta(double *theta, int n_theta, double d_theta) {

     int i;

     for (i = 0; i < n_theta; ++i)
          theta[i] = i * d_theta;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void phase_angle_exp_gc(int n_elem, int n_coef, int i,
                               double **coefs, double *p00, double *p0p2,
                               double *p2p2, double *p2m2, double **pf) {

     int j;

     double a[6];

     for (j = 0; j < n_elem; ++j)
          a[j] = 0.;

     for (j = 0; j < 2; ++j) {
          a[0] += coefs[0][j] * p00[j];
          if (n_elem > 1)
               a[3] += coefs[3][j] * p00[j];
     }

     for ( ; j < n_coef; ++j) {
          a[0] += coefs[0][j] * p00[j];
          if (n_elem > 1) {
               a[3] += coefs[3][j] * p00[j];

               a[4] += coefs[4][j] * p0p2[j];
               a[5] += coefs[5][j] * p0p2[j];

               a[1] += (coefs[1][j] + coefs[2][j]) * p2p2[j];
               a[2] += (coefs[1][j] - coefs[2][j]) * p2m2[j];
          }
     }

     pf[0][i] = a[0];
     if (n_elem > 1) {
          pf[3][i] = a[3];

          pf[4][i] = a[4];
          pf[5][i] = a[5];

          pf[1][i] = (a[1] + a[2]) / 2.;
          pf[2][i] = (a[1] - a[2]) / 2.;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_exp_gc(double **coefs, int n_coef,
                       double *theta, double **pf, int n_pf, int allsix) {

     int i;

     int n_elem;

     double *p00;
     double *p0p2;
     double *p2p2;
     double *p2m2;

     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;

     p00  = alloc_array1_d(n_coef);
     p0p2 = alloc_array1_d(n_coef);
     p2p2 = alloc_array1_d(n_coef);
     p2m2 = alloc_array1_d(n_coef);

     for (i = 0; i < n_pf; ++i) {
          gen_spher_funcs_prt(n_coef, cos(theta[i]*D2R), p00, p0p2, p2p2, p2m2);

          phase_angle_exp_gc(n_elem, n_coef, i, coefs, p00, p0p2, p2p2, p2m2, pf);
     }

     free_array1_d(p00);
     free_array1_d(p0p2);
     free_array1_d(p2p2);
     free_array1_d(p2m2);
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_exp_gc_sym(double **coefs, int n_coef,
                           double *theta, double **pf, int n_pf, int allsix) {

     int i1;
     int i2;
     int j;

     int n_pf2;
     int n_elem;

     double a;

     double *p00;
     double *p0p2;
     double *p2p2;
     double *p2m2;

     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;

     p00  = alloc_array1_d(n_coef);
     p0p2 = alloc_array1_d(n_coef);
     p2p2 = alloc_array1_d(n_coef);
     p2m2 = alloc_array1_d(n_coef);

     n_pf2 = n_pf / 2;

     for (i1 = 0; i1 < n_pf2; ++i1) {
          gen_spher_funcs_prt(n_coef, cos(theta[i1]*D2R), p00, p0p2, p2p2, p2m2);

          phase_angle_exp_gc(n_elem, n_coef, i1, coefs, p00, p0p2, p2p2, p2m2, pf);

          a = 1.;
          for (j = 0; j < n_coef; ++j) {
               p00[j]  *= a;
               p0p2[j] *= a;
               p2p2[j] *= a;
               p2m2[j] *= a;

               a *= -1.;
          }

          i2 = n_pf - i1 - 1;

          phase_angle_exp_gc(n_elem, n_coef, i2, coefs, p00, p0p2, p2m2, p2p2, pf);
     }

     if (n_pf % 2) {
          gen_spher_funcs_prt(n_coef, cos(theta[i1]*D2R), p00, p0p2, p2p2, p2m2);

          phase_angle_exp_gc(n_elem, n_coef, i1, coefs, p00, p0p2, p2p2, p2m2, pf);
     }

     free_array1_d(p00);
     free_array1_d(p0p2);
     free_array1_d(p2p2);
     free_array1_d(p2m2);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void phase_angle_exp_lc(int n_elem, int n_coef, int i,
                               double **coefs, double *p, double **pf) {

     int j;
     int k;

     double a;

     for (j = 0; j < n_elem; ++j) {
          a = 0.;
          for (k = 0; k < n_coef; ++k)
               a += coefs[j][k] * p[k];

          pf[j][i] = a;
     }

     if (n_elem > 1)
          pf[5][i] = -pf[5][i];
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_exp_lc(double **coefs, int n_coef,
                       double *theta, double **pf, int n_pf, int allsix) {

     int i;

     int n_elem;

     double *P_l;

     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;

     P_l = alloc_array1_d(n_coef);

     for (i = 0; i < n_pf; ++i) {
          leg_poly(n_coef, cos(theta[i]*D2R), P_l);

          phase_angle_exp_lc(n_elem, n_coef, i, coefs, P_l, pf);
     }

     free_array1_d(P_l);
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_exp_lc_sym(double **coefs, int n_coef,
                           double *theta, double **pf, int n_pf, int allsix) {

     int i1;
     int i2;
     int j;

     int n_pf2;

     int n_elem;

     double a;

     double *P_l;

     if (! allsix)
          n_elem = 1;
     else
          n_elem = 6;

     P_l = alloc_array1_d(n_coef);

     n_pf2 = n_pf / 2;

     for (i1 = 0; i1 < n_pf2; ++i1) {
          leg_poly(n_coef, cos(theta[i1]*D2R), P_l);

          phase_angle_exp_lc(n_elem, n_coef, i1, coefs, P_l, pf);

          a = 1.;
          for (j = 0; j < n_coef; ++j) {
               P_l[j] *= a;
               a *= -1.;
          }

          i2 = n_pf - i1 - 1;

          phase_angle_exp_lc(n_elem, n_coef, i2, coefs, P_l, pf);
     }

     if (n_pf % 2) {
          leg_poly(n_coef, cos(theta[i1]*D2R), P_l);

          phase_angle_exp_lc(n_elem, n_coef, i1, coefs, P_l, pf);
     }

     free_array1_d(P_l);
}



/*******************************************************************************
 *
 ******************************************************************************/
int create_phase_func(int n_coef, int n_theta, int n_derivs, double **xc,
                      double *theta, double **pf, double ***xc_l, double ***pf_l,
                      int flag) {

     int i;

     double d_theta;

     d_theta = phase_func_d_theta(n_theta);

     phase_func_theta(theta, n_theta, d_theta);

     if (flag) {
          phase_func_exp_gc(xc, n_coef, theta, pf, n_theta, 1);

          for (i = 0; i < n_derivs; ++i)
               phase_func_exp_gc(xc_l[i], n_coef, theta, pf_l[i], n_theta, 1);
     }
     else {
          phase_func_exp_lc(xc, n_coef, theta, pf, n_theta, 1);

          for (i = 0; i < n_derivs; ++i)
               phase_func_exp_lc(xc_l[i], n_coef, theta, pf_l[i], n_theta, 1);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
double asymmetry_parameter(int n_theta, double *theta, double **phase) {

     int i;

     double a;
     double b;

     double g;

     g = 0.;

     a = cos(theta[0]*D2R);
     for (i = 0; i < n_theta - 1; ++i) {
          b = cos(theta[i + 1]*D2R);
          g += (phase[0][i] + phase[0][i+1]) * (a + b) / 4. * (a - b);
          a = b;
     }

     return g / 2.;
}
