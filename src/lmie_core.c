/******************************************************************************%
**
**    Copyright (C) 2008-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include "lmie_core.h"


/*******************************************************************************
 *
 ******************************************************************************/
int calc_n1(double x) {

     return (int) (x + 4.05 * pow(x, 1. / 3.) + 8);
}



int calc_n2(int n1, dcomplex z) {

     return (int) (MAX(n1, cabs(z)) + 6.40 * pow(MAX(n1, cabs(z)), 1. / 3.) + 8);
}



/*******************************************************************************
 *
 ******************************************************************************/
void lmie_core_shared_import(lmie_core_shared_data *d,
                             int n_qsize, int n_derivs1, int n_derivs2,
                             int *index1, int *index2,
                             double lambda, double *lambda_l,
                             dcomplex m, dcomplex *m_l, dcomplex m2,
                             double *qx, double *qw, double *nr,
                             double *fv1, double *fv2, double *fv3, double *fv4,
                             double *qx2, double *qw2,
                             double **qx_l, double **qw_l, double **nr_l,
                             int n1, int n2, int n_qang) {

     d->n_qsize    = n_qsize;
     d->n_derivs1 = n_derivs1;
     d->n_derivs2 = n_derivs2;
     d->index1    = index1;
     d->index2    = index2;
     d->lambda    = lambda;
     d->lambda_l  = lambda_l;
     d->m         = m;
     d->m_l       = m_l;
     d->m2        = m2;
     d->qx        = qx;
     d->qw        = qw;
     d->nr        = nr;
     d->fv1       = fv1;
     d->fv2       = fv2;
     d->fv3       = fv3;
     d->fv4       = fv4;
     d->qx2       = qx2;
     d->qw2       = qw2;
     d->qx_l      = qx_l;
     d->qw_l      = qw_l;
     d->nr_l      = nr_l;
     d->n1        = n1;
     d->n2        = n2;
     d->n_qang    = n_qang;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __GNUC__
__attribute__((noinline))
#endif
void *lmie_core(lmie_core_threads_data *d) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int kk;
     int l;

     int n1;
     int n2;

     double a1;
     double a2;
     double a3;
     double a4;
     double a5;
     double a6;
     double a7;
     double a8;

     double *a1_l;
     double *a2_l;

     double x;
     double *x_l;

     double f;
     double g;
     double *f_l;
     double *f_l2;

     double cosx;
     double sinx;

     double *rn;
     double **rn_l;
     double *psi;
     double **psi_l;

     double chin0;
     double *chin0_l;
     double chin1;
     double *chin1_l;
     double chin2;
     double *chin2_l;

     double pin0;
     double pin1;
     double taun;

     double sr1;
     double sr1_l;
     double si1;
     double si1_l;
     double sr2;
     double sr2_l;
     double si2;
     double si2_l;

     dcomplex c1;
     dcomplex c2;
     dcomplex c3;
     dcomplex c4;
     dcomplex c5;
     dcomplex c6;

     dcomplex *c1_l;
/*
     dcomplex *c6_l;
*/
     dcomplex z;
     dcomplex *z_l;

     dcomplex *zeta;
     dcomplex **zeta_l;

     dcomplex *D;
     dcomplex **D_l;

     dcomplex *a;
     dcomplex **a_l;
     dcomplex *b;
     dcomplex **b_l;

     dcomplex *sum[2];
     dcomplex **sum_l[2];
     dcomplex *dif[2];
     dcomplex **dif_l[2];

     dcomplex sp[2];
     dcomplex **sp_l;
     dcomplex sm[2];
     dcomplex **sm_l;

     dcomplex s1;
     dcomplex s1_l;
     dcomplex s2;
     dcomplex s2_l;

     lmie_core_shared_data *s = d->s;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     rn   = alloc_array1_d(s->n1 + 1);
     psi  = alloc_array1_d(s->n1 + 1);

     zeta = alloc_array1_dc(s->n1 + 1);
     D    = alloc_array1_dc(s->n2 + 1);
/*
     D    = alloc_array1_dc(s->n1 + 1);
*/
     a    = alloc_array1_dc(s->n1 + 1);
     b    = alloc_array1_dc(s->n1 + 1);

     if (s->n_derivs1 > 0) {
          a1_l    = alloc_array1_d(s->n_derivs1);
          a2_l    = alloc_array1_d(s->n_derivs1);

          x_l     = alloc_array1_d(s->n_derivs1);
          f_l     = alloc_array1_d(s->n_derivs1);

          rn_l    = alloc_array2_d(s->n1 + 1, s->n_derivs1);
          psi_l   = alloc_array2_d(s->n1 + 1, s->n_derivs1);

          chin0_l = alloc_array1_d(s->n_derivs1);
          chin1_l = alloc_array1_d(s->n_derivs1);
          chin2_l = alloc_array1_d(s->n_derivs1);

          c1_l    = alloc_array1_dc(s->n_derivs1);
/*
          c6_l    = alloc_array1_dc(s->n_derivs1);
*/
          z_l     = alloc_array1_dc(s->n_derivs1);

          zeta_l  = alloc_array2_dc(s->n1 + 1, s->n_derivs1);
          D_l     = alloc_array2_dc(s->n2 + 1, s->n_derivs1);
/*
          D_l     = alloc_array2_dc(s->n1 + 1, s->n_derivs1);
*/
          a_l     = alloc_array2_dc(s->n1 + 1, s->n_derivs1);
          b_l     = alloc_array2_dc(s->n1 + 1, s->n_derivs1);

          sp_l    = alloc_array2_dc(2, s->n_derivs1);
          sm_l    = alloc_array2_dc(2, s->n_derivs1);
     }

     if (s->n_derivs2 > 0)
          f_l2    = alloc_array1_d(s->n_derivs2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < d->n_qsize; ++i) {

          if (! s->flag)
               ii = i;
          else
               ii = d->index[i];


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          x = 2. * PI * s->qx[ii] / s->lambda;
          z = s->m * x;

          f = s->nr[ii] * s->qw[ii];

          for (j = 0; j < s->n_derivs1; ++j) {
               jj = s->index1[j];

               x_l[j] = 2. * PI * (s->qx_l[ii][jj] -
                        s->qx[ii] * s->lambda_l[jj] / s->lambda) / s->lambda;
               z_l[j] = s->m_l[jj] * x + s->m * x_l[j];

               f_l[j] = s->nr_l[ii][jj] * s->qw[ii] + s->nr[ii] * s->qw_l[ii][jj];
          }

          for (j = 0; j < s->n_derivs2; ++j) {
               jj = s->index2[j];

               f_l2[j] = s->nr_l[ii][jj] * s->qw[ii];
          }

          cosx = cos(x);
          sinx = sin(x);


          n1 = calc_n1(x);
          n2 = calc_n2(n1, z);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          rn[n1] = 0.;

          for (j = 0; j < s->n_derivs1; ++j)
               rn_l[n1][j] = 0.;

          for (j = n1-1; j >= 0; --j) {
               a1 = s->fv1[j] / x;
               rn[j] = 1. / (a1 - rn[j+1]);

               if (s->n_derivs1 > 0) {
                    a2 = -a1 / x;
                    a3 = rn[j] * rn[j];
               }

               for (k = 0; k < s->n_derivs1; ++k)
                    rn_l[j][k] = -(a2 * x_l[k] - rn_l[j+1][k]) * a3;
          }


          psi[0] = sinx;
          psi[1] = rn[1] * psi[0];

          for (j = 0; j < s->n_derivs1; ++j) {
               psi_l[0][j]  = x_l[j] * cosx;
               psi_l[1][j]  = rn_l[1][j] * psi[0] + rn[1] * psi_l[0][j];
          }

          for (j = 1; j < n1; ++j) {
               psi[j+1] = rn[j+1] * psi[j];

               for (k = 0; k < s->n_derivs1; ++k)
                    psi_l[j+1][k]  = rn_l[j+1][k] * psi[j] + rn[j+1] * psi_l[j][k];
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a1 = (2. * 0. + 1.) / x;

          if (s->n_derivs1 > 0)
               a2 = -a1 / x;

          for (j = 0; j < s->n_derivs1; ++j)
               a2_l[j] = a2 * x_l[j];

          chin0   = cosx;
          chin1   = a1 * chin0 + sinx;

          zeta[0] = psi[0] + _Complex_I * chin0;
          zeta[1] = psi[1] + _Complex_I * chin1;

          for (j = 0; j < s->n_derivs1; ++j) {
               chin0_l[j]   = x_l[j] * -sinx;
               chin1_l[j]   = a2_l[j] * chin0 + a1 * chin0_l[j] + x_l[j] * cosx;

               zeta_l[0][j] = psi_l[0][j] + _Complex_I * chin0_l[j];
               zeta_l[1][j] = psi_l[1][j] + _Complex_I * chin1_l[j];
          }

          for (j = 1; j < n1; ++j) {
               a1 = s->fv1[j] / x;
               for (k = 0; k < s->n_derivs1; ++k)
                    a2_l[k] = -a1 * x_l[k] / x;

               chin2 = a1 * chin1 - chin0;

               zeta[j+1] = psi[j+1] + _Complex_I * chin2;

               for (k = 0; k < s->n_derivs1; ++k) {
                    chin2_l[k]     = a2_l[k] * chin1 + a1 * chin1_l[k] - chin0_l[k];

                    zeta_l[j+1][k] = psi_l[j+1][k] + _Complex_I * chin2_l[k];
               }

               chin0 = chin1;
               chin1 = chin2;

               for (k = 0; k < s->n_derivs1; ++k) {
                    chin0_l[k] = chin1_l[k];
                    chin1_l[k] = chin2_l[k];
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          D[n2-0] = 0.;

          for (j = 0; j < s->n_derivs1; ++j)
               D_l[n2-0][j] = 0.;

          for (j = n2-1; j >= 0; --j) {
               c1 = (j + 1.) / z ;
               c2 = (D[j+1] + c1);
               D[j] = c1 - 1. / c2;

               if (s->n_derivs1 > 0) {
                    c3 = -c1 / z;
                    c4 = c2 * c2;
               }

               for (k = 0; k < s->n_derivs1; ++k) {
                    c5 = c3 * z_l[k];
                    D_l[j][k] = c5 + (D_l[j+1][k] + c5) / c4;
               }
          }
/*
          c6 = 0.;

          for (j = 0; j < s->n_derivs1; ++j)
               c6_l[j] = 0.;

          for (j = n2-1; j >= n1; --j) {
               c1 = (j + 1.) / z ;
               c2 = (c6 + c1);
               c6 = c1 - 1. / c2;

               if (s->n_derivs1 > 0) {
                    c3 = -c1 / z;
                    c4 = c2 * c2;
               }

               for (k = 0; k < s->n_derivs1; ++k) {
                    c5 = c3 * z_l[k];
                    c6_l[k] = c5 + (c6_l[k] + c5) / c4;
               }
          }

          D[j+1] = c6;

          for (k = 0; k < s->n_derivs1; ++k)
               D_l[j+1][k] = c6_l[k];

          for (        ; j >= 0; --j) {
               c1 = (j + 1.) / z ;
               c2 = (D[j+1] + c1);
               D[j] = c1 - 1. / c2;

               if (s->n_derivs1 > 0) {
                    c3 = -c1 / z;
                    c4 = c2 * c2;
               }

               for (k = 0; k < s->n_derivs1; ++k) {
                    c5 = c3 * z_l[k];
                    D_l[j][k] = c5 + (D_l[j+1][k] + c5) / c4;
               }
          }
*/

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (j = 1; j <= n1; ++j) {
               a1 = j / x;

               for (k = 0; k < s->n_derivs1; ++k)
                    a1_l[k] = -a1 * x_l[k] / x;

               c1 = D[j] / s->m + a1;

               c2 = c1 * psi [j] - psi [j-1];
               c3 = c1 * zeta[j] - zeta[j-1];

               a[j] = c2 / c3;

               if (s->n_derivs1 > 0) {
                    c4 = D[j] / s->m2;
                    c5 = c3 * c3;
               }

               for (k = 0; k < s->n_derivs1; ++k) {
                    kk = s->index1[k];

                    c6 = D_l[j][k] / s->m - c4 * s->m_l[kk] + a1_l[k];

                    a_l[j][k] = (c6 * psi [j] + c1 * psi_l [j][k] - psi_l [j-1][k]) / c3 - c2 *
                                (c6 * zeta[j] + c1 * zeta_l[j][k] - zeta_l[j-1][k]) / c5;
               }

               c1 = s->m * D[j] + a1;

               c2 = c1 * psi [j] - psi [j-1];
               c3 = c1 * zeta[j] - zeta[j-1];

               b[j] = c2 / c3;

               if (s->n_derivs1 > 0)
                    c5 = c3 * c3;

               for (k = 0; k < s->n_derivs1; ++k) {
                    kk = s->index1[k];

                    c4 = s->m_l[kk] * D[j] + s->m * D_l[j][k] + a1_l[k];

                    b_l[j][k] = (c4 * psi [j] + c1 * psi_l [j][k] - psi_l [j-1][k]) / c3 - c2 *
                                (c4 * zeta[j] + c1 * zeta_l[j][k] - zeta_l[j-1][k]) / c5;
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a1 = 0.;
          a2 = 0.;
          c1 = 0.;

          for (j = 0; j < s->n_derivs1; ++j) {
               a1_l[j] = 0.;
               a2_l[j] = 0.;
               c1_l[j] = 0.;
          }

          a3 = -1.;
          for (j = 1; j <= n1; ++j) {
               a1 += s->fv1[j] * creal(a[j] + b[j]);

               a2 += s->fv1[j] * (creal(a[j] * conj(a[j])) + creal(b[j] * conj(b[j])));

               c1 += s->fv1[j] * (a[j] - b[j]) * a3;

               for (k = 0; k < s->n_derivs1; ++k) {
                    a1_l[k] += s->fv1[j] * creal(a_l[j][k] + b_l[j][k]);

                    a2_l[k] += s->fv1[j] * 2. * (creal(a[j] * conj(a_l[j][k])) + creal(b[j] * conj(b_l[j][k])));

                    c1_l[k] += s->fv1[j] * (a_l[j][k] - b_l[j][k]) * a3;
               }

               a3 = -a3;
          }

          *d->cext += f *  a1;
          *d->csca += f *  a2;
          *d->cbak += f * creal(c1 * conj(c1));

          for (j = 0; j < s->n_derivs1; ++j) {
               jj = s->index1[j];

               d->cext_l[jj] += f_l[j] * a1 + f * a1_l[j];
               d->csca_l[jj] += f_l[j] * a2 + f * a2_l[j];
               d->cbak_l[jj] += f_l[j] * creal(c1 * conj(c1)) + f * 2. * creal(c1 * conj(c1_l[j]));
          }

          for (j = 0; j < s->n_derivs2; ++j) {
               jj = s->index2[j];

               d->cext_l[jj] += f_l2[j] * a1;
               d->csca_l[jj] += f_l2[j] * a2;
               d->cbak_l[jj] += f_l2[j] * creal(c1 * conj(c1));
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a1 = 0.;

          for (j = 0; j < s->n_derivs1; ++j)
               a1_l[j] = 0.;

          for (j = 1; j < n1; ++j) {
               a1 += s->fv3[j] * creal(a[j] * conj(a[j + 1]) + b[j] * conj(b[j + 1])) +
                     s->fv4[j] * creal(a[j] * conj(b[j]));

               for (k = 0; k < s->n_derivs1; ++k) {
                    a1_l[k] += s->fv3[j] * creal(a_l[j][k] * conj(a[j + 1]) + a[j] * conj(a_l[j + 1][k]) +
                                                 b_l[j][k] * conj(b[j + 1]) + b[j] * conj(b_l[j + 1][k])) +
                               s->fv4[j] * creal(a_l[j][k] * conj(b[j]) + a[j] * conj(b_l[j][k]));
               }
          }

          *d->g += f * a1;

          for (j = 0; j < s->n_derivs1; ++j) {
               jj = s->index1[j];

               d->g_l[jj] += f_l[j] * a1 + f * a1_l[j];
          }

          for (j = 0; j < s->n_derivs2; ++j) {
               jj = s->index2[j];

               d->g_l[jj] += f_l2[j] * a1;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          sum[0] = zeta;
          dif[0] = D;
          sum[1] = a;
          dif[1] = b;

          if (s->n_derivs1 > 0) {
               sum_l[0] = zeta_l;
               dif_l[0] = D_l;
               sum_l[1] = a_l;
               dif_l[1] = b_l;
          }

          a1 = 1.;
          for (k = 1; k <= n1; ++k) {
               c1 = a[k] + b[k];
               c2 = a[k] - b[k];

               sum[0][k] = s->fv4[k] * c1;
               dif[0][k] = s->fv4[k] * c2;

               sum[1][k] = sum[0][k] * a1;
               dif[1][k] = dif[0][k] * a1;

               for (l = 0; l < s->n_derivs1; ++l) {
                    c1 = a_l[k][l] + b_l[k][l];
                    c2 = a_l[k][l] - b_l[k][l];

                    sum_l[0][k][l] = s->fv4[k] * c1;
                    dif_l[0][k][l] = s->fv4[k] * c2;

                    sum_l[1][k][l] = sum_l[0][k][l] * a1;
                    dif_l[1][k][l] = dif_l[0][k][l] * a1;
               }

               a1 = -a1;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          g = f * .5;

          for (j = 0; j < s->n_qang / 2; ++j) {
               sp[0] = sp[1] = 0.;
               sm[0] = sm[1] = 0.;

               for (k = 0; k < s->n_derivs1; ++k) {
                    sp_l[0][k] = sp_l[1][k] = 0.;
                    sm_l[0][k] = sm_l[1][k] = 0.;
               }

               pin0 = 0.;
               pin1 = 1.;

               if (s->n_derivs1 == 0) {
                    for (k = 1; k <= n1; ++k) {
                         a1 = s->qx2[j] * pin1;
                         a2 = a1 - pin0;

                         taun = k * a2 - pin0;

                         a3 = pin1 + taun;
                         a4 = pin1 - taun;

                         sp[0] += a3 * sum[0][k];
                         sm[0] += a4 * dif[0][k];

                         sp[1] += a4 * sum[1][k];
                         sm[1] += a3 * dif[1][k];

                         pin0 = pin1;

                         pin1 = a1 + s->fv2[k] * a2;
                    }
               }
               else {
                    for (k = 1; k <= n1; ++k) {
                         a1 = s->qx2[j] * pin1;
                         a2 = a1 - pin0;

                         taun = k * a2 - pin0;

                         a3 = pin1 + taun;
                         a4 = pin1 - taun;

                         sp[0] += a3 * sum[0][k];
                         sm[0] += a4 * dif[0][k];

                         sp[1] += a4 * sum[1][k];
                         sm[1] += a3 * dif[1][k];

                         for (l = 0; l < s->n_derivs1; ++l) {
                              sp_l[0][l] += a3 * sum_l[0][k][l];
                              sm_l[0][l] += a4 * dif_l[0][k][l];

                              sp_l[1][l] += a4 * sum_l[1][k][l];
                              sm_l[1][l] += a3 * dif_l[1][k][l];
                         }

                         pin0 = pin1;

                         pin1 = a1 + s->fv2[k] * a2;
                    }
               }

               for (l = 0; l < 2; ++l) {
                    jj = l == 0 ? j : s->n_qang - j - 1;

                    s1 = .5 * (sp[l] + sm[l]);
                    s2 = .5 * (sp[l] - sm[l]);

                    sr1 = creal(s1);
                    si1 = cimag(s1);
                    sr2 = creal(s2);
                    si2 = cimag(s2);

                    a1 = sr1 * sr1 + si1 * si1;
                    c2 = s1 * conj(s2);
                    a2 = creal(c2);
                    a3 = sr2 * sr2 + si2 * si2;
                    c4 = s2 * conj(s1);
                    a4 = creal(c4);

                    a5 =  (a1 + a3) * .5;
                    a6 =  (a2 + a4) * .5;
                    a7 = -(a1 - a3) * .5;
                    a8 = -creal((c2 - c4) * .5 * _Complex_I);

                    d->pf[0][jj] += f * a5;
                    d->pf[1][jj] += f * a6;
                    d->pf[2][jj] += f * a7;
                    d->pf[3][jj] += f * a8;

                    for (k = 0; k < s->n_derivs1; ++k) {
                         s1_l = .5 * (sp_l[l][k] + sm_l[l][k]);
                         s2_l = .5 * (sp_l[l][k] - sm_l[l][k]);

                         sr1_l = creal(s1_l);
                         si1_l = cimag(s1_l);
                         sr2_l = creal(s2_l);
                         si2_l = cimag(s2_l);

                         a1 = 2. * (sr1_l * sr1 + si1_l * si1);
                         c2 = 2. * s1_l * conj(s2);
                         a2 = creal(c2);
                         a3 = 2. * (sr2_l * sr2 + si2_l * si2);
                         c4 = 2. * s2_l * conj(s1);
                         a4 = creal(c4);

                         kk = s->index1[k];

                         d->pf_l[kk][0][jj] += f_l[k] * a5 + g *        (a1 + a3);
                         d->pf_l[kk][1][jj] += f_l[k] * a6 + g *        (a2 + a4);
                         d->pf_l[kk][2][jj] += f_l[k] * a7 + g * -      (a1 - a3);
                         d->pf_l[kk][3][jj] += f_l[k] * a8 + g * -creal((c2 - c4) * _Complex_I);
                    }

                    for (k = 0; k < s->n_derivs2; ++k) {
                         kk = s->index2[k];

                         d->pf_l[kk][0][jj] += a5 * f_l2[k];
                         d->pf_l[kk][1][jj] += a6 * f_l2[k];
                         d->pf_l[kk][2][jj] += a7 * f_l2[k];
                         d->pf_l[kk][3][jj] += a8 * f_l2[k];
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(rn);
     free_array1_d(psi);

     free_array1_dc(zeta);
     free_array1_dc(D);

     free_array1_dc(a);
     free_array1_dc(b);

     if (s->n_derivs1 > 0) {
          free_array1_d(a1_l);
          free_array1_d(a2_l);
          free_array1_d(x_l);
          free_array1_d(f_l);

          free_array2_d(rn_l);
          free_array2_d(psi_l);

          free_array1_d(chin0_l);
          free_array1_d(chin1_l);
          free_array1_d(chin2_l);

          free_array1_dc(c1_l);
/*
          free_array1_dc(c6_l);
*/
          free_array1_dc(z_l);

          free_array2_dc(zeta_l);

          free_array2_dc(D_l);

          free_array2_dc(a_l);
          free_array2_dc(b_l);

          free_array2_dc(sp_l);
          free_array2_dc(sm_l);
     }

     if (s->n_derivs2 > 0)
          free_array1_d(f_l2);


     return d;
}
