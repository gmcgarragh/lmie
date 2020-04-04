/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

/*******************************************************************************
 *
 ******************************************************************************/
__device__ __host__ int kernel_work_n(int n_derivs, int n1, int n2) {

     int n = 0;

     n += (n1 + 1) * sizeof(FLOAT);
     n += (n1 + 1) * sizeof(FLOAT);

     n += (n1 + 1) * sizeof(COMPLEX);
     n += (n2 + 1) * sizeof(COMPLEX);
     n += (n1 + 1) * sizeof(COMPLEX);
     n += (n1 + 1) * sizeof(COMPLEX);

     if (n_derivs > 0) {
          n += (n1 + 1) * n_derivs * sizeof(FLOAT);
          n += (n1 + 1) * n_derivs * sizeof(FLOAT);

          n += (n1 + 1) * n_derivs * sizeof(COMPLEX);
          n += (n2 + 1) * n_derivs * sizeof(COMPLEX);
          n += (n1 + 1) * n_derivs * sizeof(COMPLEX);
          n += (n1 + 1) * n_derivs * sizeof(COMPLEX);
     }

     return n;
}



/*******************************************************************************
 *
 ******************************************************************************/
__global__ void kernel(int i_qsiz, int n_qsiz, int n_derivs, FLOAT lambda,
                       COMPLEX m, COMPLEX *m_l, COMPLEX m2,
                       FLOAT *qx, FLOAT *qw, FLOAT *nr,
                       FLOAT *fv1, FLOAT *fv2, FLOAT *fv3, FLOAT *fv4,
                       FLOAT *qx2, FLOAT *qw2,
                       FLOAT *cext, FLOAT *csca, FLOAT *cbak,
                       FLOAT *g, FLOAT *pf,
                       FLOAT *qx_l, FLOAT *qw_l, FLOAT *nr_l,
                       FLOAT *cext_l, FLOAT *csca_l, FLOAT *cbak_l,
                       FLOAT *g_l, FLOAT *pf_l,
                       int n1, int n2, int n_qang, unsigned char *work) {

     int i;
     int j;
     int k;
     int l;

     FLOAT a1;
     FLOAT a2;
     FLOAT a3;
     FLOAT a4;
     FLOAT a5;
     FLOAT a6;
     FLOAT a7;
     FLOAT a8;

     FLOAT a1_l[MAX_DERIVS];
     FLOAT a2_l[MAX_DERIVS];

     FLOAT x;
     FLOAT x_l[MAX_DERIVS];

     FLOAT f;
     FLOAT f_l1[MAX_DERIVS];
     FLOAT f_l2[MAX_DERIVS];

     FLOAT cosx;
     FLOAT sinx;

     FLOAT *rn;
     FLOAT *rn_l;
     FLOAT *psi;
     FLOAT *psi_l;

     FLOAT chin0;
     FLOAT chin0_l[MAX_DERIVS];
     FLOAT chin1;
     FLOAT chin1_l[MAX_DERIVS];
     FLOAT chin2;
     FLOAT chin2_l[MAX_DERIVS];

     FLOAT pin0;
     FLOAT pin1;
     FLOAT taun;

     FLOAT sr1;
     FLOAT si1;
     FLOAT sr2;
     FLOAT si2;

     COMPLEX c1;
     COMPLEX c2;
     COMPLEX c3;
     COMPLEX c4;
     COMPLEX c5;
     COMPLEX c6;

     COMPLEX c1_l[MAX_DERIVS];

     COMPLEX c6_l[MAX_DERIVS];

     COMPLEX z;
     COMPLEX z_l[MAX_DERIVS];

     COMPLEX *zeta;
     COMPLEX *zeta_l;
     COMPLEX *D;
     COMPLEX *D_l;

     COMPLEX *a;
     COMPLEX *a_l;
     COMPLEX *b;
     COMPLEX *b_l;

     COMPLEX *sum;
     COMPLEX *sum_l;
     COMPLEX *dif;
     COMPLEX *dif_l;

     COMPLEX sp;
     COMPLEX sp_l[MAX_DERIVS];
     COMPLEX sm;
     COMPLEX sm_l[MAX_DERIVS];

     COMPLEX s1;
     COMPLEX s2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i_qsiz = blockDim.x * blockIdx.x + threadIdx.x;

     if (i_qsiz >= n_qsiz)
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i = i_qsiz * kernel_work_n(n_derivs, n1, n2);

     rn   = (FLOAT *)   (work + i); i += (n1 + 1) * sizeof(FLOAT);
     psi  = (FLOAT *)   (work + i); i += (n1 + 1) * sizeof(FLOAT);

     zeta = (COMPLEX *) (work + i); i += (n1 + 1) * sizeof(COMPLEX);
     D    = (COMPLEX *) (work + i); i += (n2 + 1) * sizeof(COMPLEX);
     a    = (COMPLEX *) (work + i); i += (n1 + 1) * sizeof(COMPLEX);
     b    = (COMPLEX *) (work + i); i += (n1 + 1) * sizeof(COMPLEX);

     if (n_derivs > 0) {
          rn_l   = (FLOAT *)   (work + i); i += (n1 + 1) * n_derivs * sizeof(FLOAT);
          psi_l  = (FLOAT *)   (work + i); i += (n1 + 1) * n_derivs * sizeof(FLOAT);

          zeta_l = (COMPLEX *) (work + i); i += (n1 + 1) * n_derivs * sizeof(COMPLEX);
          D_l    = (COMPLEX *) (work + i); i += (n2 + 1) * n_derivs * sizeof(COMPLEX);
          a_l    = (COMPLEX *) (work + i); i += (n1 + 1) * n_derivs * sizeof(COMPLEX);
          b_l    = (COMPLEX *) (work + i); i += (n1 + 1) * n_derivs * sizeof(COMPLEX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i = i_qsiz;

     x = 2. * PI * qx[i] / lambda;
     z = m * x;

     f = nr[i] * qw[i];

     for (j = 0; j < n_derivs; ++j) {
          x_l[j] = 2. * PI * qx_l[A2(i,j)] / lambda;
          z_l[j] = m_l[j] * x + m * x_l[j];

          f_l1[j] = nr_l[A2(i,j)] * qw[i];
          f_l2[j] = nr  [i] * qw_l[A2(i,j)];
     }


     cosx = cos(x);
     sinx = sin(x);


     n1 =               x   + 4.05f * powf(              x,   1.f/3.f) + 8;
     n2 = MAX(n1, abssc(z)) + 6.40f * powf(MAX(n1, abssc(z)), 1.f/3.f) + 8;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     rn[n1] = 0.;

     for (j = 0; j < n_derivs; ++j)
          rn_l[A2(n1,j)] = 0.;

     for (j = n1-1; j >= 0; --j) {
          a1 = fv1[j] / x;
          rn[j] = 1. / (a1 - rn[j+1]);

          if (n_derivs > 0) {
               a2 = -a1 / x;
               a3 = rn[j] * rn[j];
          }

          for (k = 0; k < n_derivs; ++k)
               rn_l[A2(j,k)] = -(a2 * x_l[k] - rn_l[A2(j+1,k)]) * a3;
     }


     psi[0]  = sinx;
     psi[1]  = rn[1] * psi[0];

     for (j = 0; j < n_derivs; ++j) {
          psi_l[A2(0,j)]  = x_l[j] * cosx;
          psi_l[A2(1,j)]  = rn_l[A2(1,j)] * psi[0] + rn[1] * psi_l[A2(0,j)];
     }

     for (j = 1; j < n1; ++j) {
          psi[j+1] = rn[j+1] * psi[j];

          for (k = 0; k < n_derivs; ++k)
               psi_l[A2(j+1,k)]  = rn_l[A2(j+1,k)] * psi[j] + rn[j+1] * psi_l[A2(j,k)];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a1 = (2. * 0. + 1.) / x;

     if (n_derivs > 0)
          a2 = -a1 / x;

     for (j = 0; j < n_derivs; ++j)
          a2_l[j] = a2 * x_l[j];

     chin0   = cosx;
     chin1   = a1 * chin0 + sinx;

     zeta[0] = psi[0] + MAKE_COMPLEX(0., 1.f) * chin0;
     zeta[1] = psi[1] + MAKE_COMPLEX(0., 1.f) * chin1;

     for (j = 0; j < n_derivs; ++j) {
          chin0_l[j]      = x_l[j] * -sinx;
          chin1_l[j]      = a2_l[j] * chin0 + a1 * chin0_l[j] + x_l[j] * cosx;

          zeta_l[A2(0,j)] = psi_l[A2(0,j)] + MAKE_COMPLEX(0., 1.f) * chin0_l[j];
          zeta_l[A2(1,j)] = psi_l[A2(1,j)] + MAKE_COMPLEX(0., 1.f) * chin1_l[j];
     }

     for (j = 1; j < n1; ++j) {
          a1 = fv1[j] / x;
          for (k = 0; k < n_derivs; ++k)
               a2_l[k] = -a1 * x_l[k] / x;

          chin2 = a1 * chin1 - chin0;

          zeta[j+1] = psi[j+1] + MAKE_COMPLEX(0., 1.f) * chin2;

          for (k = 0; k < n_derivs; ++k) {
               chin2_l[k]        = a2_l[k] * chin1 + a1 * chin1_l[k] - chin0_l[k];

               zeta_l[A2(j+1,k)] = psi_l[A2(j+1,k)] + MAKE_COMPLEX(0., 1.f) * chin2_l[k];
          }

          chin0 = chin1;
          chin1 = chin2;

          for (k = 0; k < n_derivs; ++k) {
               chin0_l[k] = chin1_l[k];
               chin1_l[k] = chin2_l[k];
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     D[n2-0] = 0.;

     for (j = 0; j < n_derivs; ++j)
          D_l[A2(n2-0,j)] = 0.;

     for (j = n2-1; j >= 0; --j) {
          c1 = (j + 1.) / z ;
          c2 = (D[j+1] + c1);
          D[j] = c1 - 1. / c2;

          if (n_derivs > 0) {
               c3 = -c1 / z;
               c4 = c2 * c2;
          }

          for (k = 0; k < n_derivs; ++k) {
               c5 = c3 * z_l[k];
               D_l[A2(j,k)] = c5 + (D_l[A2(j+1,k)] + c5) / c4;
          }
     }
*/
     c6 = 0.;

     for (j = 0; j < n_derivs; ++j)
          c6_l[j] = 0.;

     for (j = n2-1; j >= n1; --j) {
          c1 = (j + 1.) / z ;
          c2 = (c6 + c1);
          c6 = c1 - 1. / c2;

          if (n_derivs > 0) {
               c3 = -c1 / z;
               c4 = c2 * c2;
          }

          for (k = 0; k < n_derivs; ++k) {
               c5 = c3 * z_l[k];
               c6_l[k] = c5 + (c6_l[k] + c5) / c4;
          }
     }

     D[j+1] = c6;

     for (k = 0; k < n_derivs; ++k)
          D_l[A2(j+1,k)] = c6_l[k];

     for (        ; j >= 0; --j) {
          c1 = (j + 1.) / z ;
          c2 = (D[j+1] + c1);
          D[j] = c1 - 1. / c2;

          if (n_derivs > 0) {
               c3 = -c1 / z;
               c4 = c2 * c2;
          }

          for (k = 0; k < n_derivs; ++k) {
               c5 = c3 * z_l[k];
               D_l[A2(j,k)] = c5 + (D_l[A2(j+1,k)] + c5) / c4;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 1; j <= n1; ++j) {
          a1 = j / x;

          for (k = 0; k < n_derivs; ++k)
               a1_l[k] = -a1 * x_l[k] / x;

          c1 = D[j] / m + a1;

          c2 = c1 * psi [j] - psi [j-1];
          c3 = c1 * zeta[j] - zeta[j-1];

          a[j] = c2 / c3;

          if (n_derivs > 0) {
               c4 = D[j] / m2;
               c5 = c3 * c3;
          }

          for (k = 0; k < n_derivs; ++k) {
               c6 = D_l[A2(j,k)] / m - c4 * m_l[k] + a1_l[k];

               a_l[A2(j,k)] = (c6 * psi [j] + c1 * psi_l [A2(j,k)] - psi_l [A2(j-1,k)]) / c3 - c2 *
                              (c6 * zeta[j] + c1 * zeta_l[A2(j,k)] - zeta_l[A2(j-1,k)]) / c5;
          }

          c1 = m * D[j] + a1;

          c2 = c1 * psi [j] - psi [j-1];
          c3 = c1 * zeta[j] - zeta[j-1];

          b[j] = c2 / c3;

          if (n_derivs > 0)
               c5 = c3 * c3;

          for (k = 0; k < n_derivs; ++k) {
               c4 = m_l[k] * D[j] + m * D_l[A2(j,k)] + a1_l[k];

               b_l[A2(j,k)] = (c4 * psi [j] + c1 * psi_l [A2(j,k)] - psi_l [A2(j-1,k)]) / c3 - c2 *
                              (c4 * zeta[j] + c1 * zeta_l[A2(j,k)] - zeta_l[A2(j-1,k)]) / c5;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a1 = 0.;
     a2 = 0.;
     c1 = 0.;

     for (j = 0; j < n_derivs; ++j) {
          a1_l[j] = 0.;
          a2_l[j] = 0.;
          c1_l[j] = 0.;
     }

     a3 = -1.;
     for (j = 1; j <= n1; ++j) {
          a1 += fv1[j] * (realxc(a[j]) + realxc(b[j]));

          a2 += fv1[j] * (realxc(a[j] * ~a[j]) + realxc(b[j] * ~b[j]));

          c1 += fv1[j] * (a[j] - b[j]) * a3;

          for (k = 0; k < n_derivs; ++k) {
               a1_l[k] += fv1[j] * (realxc(a_l[A2(j,k)]) + realxc(b_l[A2(j,k)]));

               a2_l[k] += fv1[j] * (2. * (realxc(a[j] * ~a_l[A2(j,k)]) + realxc(b[j] * ~b_l[A2(j,k)])));

               c1_l[k] += fv1[j] * (a_l[A2(j,k)] - b_l[A2(j,k)]) * a3;
          }

          a3 = -a3;
     }
#ifndef USE_ATOMIC
     *cext += f * a1;
     *csca += f * a2;
     *cbak += f * realxc(c1 * ~c1);
#else
     atomicAdd_x(cext, f * a1);
     atomicAdd_x(csca, f * a2);
     atomicAdd_x(cbak, f * realxc(c1 * ~c1));
#endif
     for (j = 0; j < n_derivs; ++j) {
#ifndef USE_ATOMIC
          cext_l[j] += f_l1[j] * a1 + f * a1_l[j] + f_l2[j] * a1;
          csca_l[j] += f_l1[j] * a2 + f * a2_l[j] + f_l2[j] * a2;
          cbak_l[j] += f_l1[j] * realxc(c1 * ~c1) + f * 2. * realxc(c1 * ~c1_l[j]);
#else
          atomicAdd_x(&cext_l[j], f_l1[j] * a1 + f * a1_l[j] + f_l2[j] * a1);
          atomicAdd_x(&csca_l[j], f_l1[j] * a2 + f * a2_l[j] + f_l2[j] * a2);
          atomicAdd_x(&cbak_l[j], f_l1[j] * realxc(c1 * ~c1) + f * 2. * realxc(c1 * ~c1_l[j]));
#endif
     }


     /*--------------------------------------------------------------------
      *
      *------------------------------------------------------------------*/
     a1 = 0.;

     for (j = 0; j < n_derivs; ++j)
          a1_l[j] = 0.;

     for (j = 1; j < n1; ++j) {
          a1 += fv3[j] * realxc(a[j] * ~a[j + 1] + b[j] * ~b[j + 1]) +
                fv4[j] * realxc(a[j] * ~b[j]);

          for (k = 0; k < n_derivs; ++k) {
               a1_l[k] += fv3[j] * realxc(a_l[A2(j,k)] * ~a[j + 1] + a[j] * ~a_l[A2(j + 1,k)] +
                                          b_l[A2(j,k)] * ~b[j + 1] + b[j] * ~b_l[A2(j + 1,k)]) +
                          fv4[j] * realxc(a_l[A2(j,k)] * ~b[j] + a[j] * ~b_l[A2(j,k)]);
          }
     }
#ifndef USE_ATOMIC
     *g += f * a1;
#else
     atomicAdd_x(g, f * a1);
#endif
     for (j = 0; j < n_derivs; ++j)
#ifndef USE_ATOMIC
          g_l[j] += f_l1[j] * a1 + f * a1_l[j];
#else
          atomicAdd_x(&g_l[j], f_l1[j] * a1 + f * a1_l[j]);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     sum = a;
     dif = b;

     if (n_derivs > 0) {
          sum_l = a_l;
          dif_l = b_l;
     }

     for (j = 1; j <= n1; ++j) {
          c1 = a[j] + b[j];
          c2 = a[j] - b[j];

          sum[j] = fv4[j] * c1;
          dif[j] = fv4[j] * c2;

          for (k = 0; k < n_derivs; ++k) {
               c1 = a_l[A2(j,k)] + b_l[A2(j,k)];
               c2 = a_l[A2(j,k)] - b_l[A2(j,k)];

               sum_l[A2(j,k)] = fv4[j] * c1;
               dif_l[A2(j,k)] = fv4[j] * c2;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < n_qang; ++j) {
          sp = 0.;
          sm = 0.;

          for (k = 0; k < n_derivs; ++k) {
               sp_l[k] = 0.;
               sm_l[k] = 0.;
          }

          a5 = qx2[j];

          pin0 = 0.;
          pin1 = 1.;
if (n_derivs == 0) {
          for (k = 1; k <= n1; ++k) {
               a1 = a5 * pin1;
               a2 = a1 - pin0;

               taun = k * a2 - pin0;

               a3 = pin1 + taun;
               a4 = pin1 - taun;

               sp += a3 * sum[k];
               sm += a4 * dif[k];

               for (l = 0; l < n_derivs; ++l) {
                    sp_l[l] += a3 * sum_l[A2(k,l)];
                    sm_l[l] += a4 * dif_l[A2(k,l)];
               }

               pin0 = pin1;

               pin1 = a1 + fv2[k] * a2;
          }
}
else {
          for (k = 1; k <= n1; ++k) {
               a1 = a5 * pin1;
               a2 = a1 - pin0;

               taun = k * a2 - pin0;

               a3 = pin1 + taun;
               a4 = pin1 - taun;

               sp += a3 * sum[k];
               sm += a4 * dif[k];

               for (l = 0; l < n_derivs; ++l) {
                    sp_l[l] += a3 * sum_l[A2(k,l)];
                    sm_l[l] += a4 * dif_l[A2(k,l)];
               }

               pin0 = pin1;

               pin1 = a1 + fv2[k] * a2;
          }
}
          s1 = .5f * (sp + sm);
          s2 = .5f * (sp - sm);

          sr1 = realxc(s1);
          si1 = imagxc(s1);
          sr2 = realxc(s2);
          si2 = imagxc(s2);

          a1 = sr1 * sr1 + si1 * si1;
          c2 = s1 * ~s2;
          a2 = realxc(c2);
          a3 = sr2 * sr2 + si2 * si2;
          c4 = s2 * ~s1;
          a4 = realxc(c4);

          a5 =  (a1 + a3) * .5;
          a6 =  (a2 + a4) * .5;
          a7 = -(a1 - a3) * .5;
          a8 = -realxc((c2 - c4) * MAKE_COMPLEX(0., .5f));
#ifndef USE_ATOMIC
          pf[0 * n_qang + j] += f * a5;
          pf[1 * n_qang + j] += f * a6;
          pf[2 * n_qang + j] += f * a7;
          pf[3 * n_qang + j] += f * a8;
#else
          atomicAdd_x(&pf[0 * n_qang + j], f * a5);
          atomicAdd_x(&pf[1 * n_qang + j], f * a6);
          atomicAdd_x(&pf[2 * n_qang + j], f * a7);
          atomicAdd_x(&pf[3 * n_qang + j], f * a8);
#endif
          for (k = 0; k < n_derivs; ++k) {

          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
extern "C"
int lmie_core_cuda(int n_qsiz, int n_derivs,
                   double lambda, double *lambda_l,
                   doublecomplex m, doublecomplex *m_l, doublecomplex m2,
                   double *qx, double *qw, double *nr,
                   double *fv1, double *fv2, double *fv3, double *fv4,
                   double *qx2, double *qw2,
                   double *cext, double *csca, double *cbak,
                   double *g, double **pf,
                   double **qx_l, double **qw_l, double **nr_l,
                   double *cext_l, double *csca_l, double *cbak_l,
                   double *g_l, double ***pf_l,
                   int n1, int n2, int n_qang) {

     int i;
     int j;
     int k;

     int threadsPerBlock;
     int blocksPerGrid;

     FLOAT h_cext;
     FLOAT h_cext_l[MAX_DERIVS];
     FLOAT h_csca;
     FLOAT h_csca_l[MAX_DERIVS];
     FLOAT h_cbak;
     FLOAT h_cbak_l[MAX_DERIVS];
     FLOAT h_g;
     FLOAT h_g_l[MAX_DERIVS];

     FLOAT *h_qx;
     FLOAT *h_qx_l;
     FLOAT *h_qw;
     FLOAT *h_qw_l;
     FLOAT *h_nr;
     FLOAT *h_nr_l;
     FLOAT *h_qx2;
     FLOAT *h_qw2;
     FLOAT *h_fv1;
     FLOAT *h_fv2;
     FLOAT *h_fv3;
     FLOAT *h_fv4;
     FLOAT *h_pf;
     FLOAT *h_pf_l;

     FLOAT *d_qx;
     FLOAT *d_qx_l;
     FLOAT *d_qw;
     FLOAT *d_qw_l;
     FLOAT *d_nr;
     FLOAT *d_nr_l;
     FLOAT *d_qx2;
     FLOAT *d_qw2;
     FLOAT *d_fv1;
     FLOAT *d_fv2;
     FLOAT *d_fv3;
     FLOAT *d_fv4;

     FLOAT *d_cext;
     FLOAT *d_cext_l;
     FLOAT *d_csca;
     FLOAT *d_csca_l;
     FLOAT *d_cbak;
     FLOAT *d_cbak_l;
     FLOAT *d_g;
     FLOAT *d_g_l;
     FLOAT *d_pf;
     FLOAT *d_pf_l;

     COMPLEX h_m;
     COMPLEX h_m_l[MAX_DERIVS];
     COMPLEX h_m2;

     COMPLEX *d_m_l;

     unsigned char *d_work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     threadsPerBlock = THREADSPERBLOCK;

     blocksPerGrid   = (n_qsiz + threadsPerBlock - 1) / threadsPerBlock;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     h_qx  = (FLOAT *) malloc(n_qsiz     * sizeof(FLOAT));
     h_qw  = (FLOAT *) malloc(n_qsiz     * sizeof(FLOAT));
     h_nr  = (FLOAT *) malloc(n_qsiz     * sizeof(FLOAT));
     h_qx2 = (FLOAT *) malloc(n_qang     * sizeof(FLOAT));
     h_qw2 = (FLOAT *) malloc(n_qang     * sizeof(FLOAT));
     h_fv1 = (FLOAT *) malloc((n1 + 1)   * sizeof(FLOAT));
     h_fv2 = (FLOAT *) malloc((n1 + 1)   * sizeof(FLOAT));
     h_fv3 = (FLOAT *) malloc((n1 + 1)   * sizeof(FLOAT));
     h_fv4 = (FLOAT *) malloc((n1 + 1)   * sizeof(FLOAT));

     h_pf  = (FLOAT *) malloc(4 * n_qang * sizeof(FLOAT));

     if (n_derivs > 0) {
          h_qx_l = (FLOAT *) malloc(n_qsiz * n_derivs * sizeof(FLOAT));
          h_qw_l = (FLOAT *) malloc(n_qsiz * n_derivs * sizeof(FLOAT));
          h_nr_l = (FLOAT *) malloc(n_qsiz * n_derivs * sizeof(FLOAT));

          h_pf_l = (FLOAT *) malloc(n_derivs * 4 * n_qang * sizeof(FLOAT));
     }


     CUDA_CALL_CHECK(cudaMalloc((void **) &d_qx,   n_qsiz     * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_qw,   n_qsiz     * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_nr,   n_qsiz     * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_qx2,  n_qang     * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_qw2,  n_qang     * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_fv1,  (n1 + 1)   * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_fv2,  (n1 + 1)   * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_fv3,  (n1 + 1)   * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_fv4,  (n1 + 1)   * sizeof(FLOAT)));

     CUDA_CALL_CHECK(cudaMalloc((void **) &d_cext, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_csca, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_cbak, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_g,    1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_pf,   4 * n_qang * sizeof(FLOAT)));

     if (n_derivs > 0) {
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_m_l,           n_derivs * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_qx_l, n_qsiz * n_derivs * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_qw_l, n_qsiz * n_derivs * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_nr_l, n_qsiz * n_derivs * sizeof(FLOAT)));

          CUDA_CALL_CHECK(cudaMalloc((void **) &d_cext_l, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_csca_l, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_cbak_l, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_g_l,    n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_pf_l,   n_derivs * 4 * n_qang * sizeof(FLOAT)));
     }

     CUDA_CALL_CHECK(cudaMalloc((void **) &d_work, blocksPerGrid * threadsPerBlock * kernel_work_n(n_derivs, n1, n2)));


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     h_m = MAKE_COMPLEX(realdc(m), imagdc(m));

     for (i = 0; i < n_derivs; ++i)
          h_m_l[i] = MAKE_COMPLEX(realdc(m_l[i]), imagdc(m_l[i]));

     h_m2 = MAKE_COMPLEX(realdc(m2), imagdc(m2));

     for (i = 0; i < n_qsiz; ++i) {
          h_qx[i] = qx[i];
          h_qw[i] = qw[i];
          h_nr[i] = nr[i];

          for (j = 0; j < n_derivs; ++j) {
               h_qx_l[A2(i,j)] = qx_l[i][j];
               h_qw_l[A2(i,j)] = qw_l[i][j];
               h_nr_l[A2(i,j)] = nr_l[i][j];
          }
     }

     for (i = 0; i < n_qang; ++i) {
          h_qx2[i] = qx2[i];
          h_qw2[i] = qw2[i];
     }

     for (i = 0; i < n1 + 1; ++i) {
          h_fv1[i] = fv1[i];
          h_fv2[i] = fv2[i];
          h_fv3[i] = fv3[i];
          h_fv4[i] = fv4[i];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     CUDA_CALL_CHECK(cudaMemcpy(d_qx,  h_qx,  n_qsiz   * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_qw,  h_qw,  n_qsiz   * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_nr,  h_nr,  n_qsiz   * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_qx2, h_qx2, n_qang   * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_qw2, h_qw2, n_qang   * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_fv1, h_fv1, (n1 + 1) * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_fv2, h_fv2, (n1 + 1) * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_fv3, h_fv3, (n1 + 1) * sizeof(FLOAT), cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpy(d_fv4, h_fv4, (n1 + 1) * sizeof(FLOAT), cudaMemcpyHostToDevice));

     if (n_derivs > 0) {
          CUDA_CALL_CHECK(cudaMemcpy(d_m_l,  h_m_l,           n_derivs * sizeof(FLOAT), cudaMemcpyHostToDevice));
          CUDA_CALL_CHECK(cudaMemcpy(d_qx_l, h_qx_l, n_qsiz * n_derivs * sizeof(FLOAT), cudaMemcpyHostToDevice));
          CUDA_CALL_CHECK(cudaMemcpy(d_qw_l, h_qw_l, n_qsiz * n_derivs * sizeof(FLOAT), cudaMemcpyHostToDevice));
          CUDA_CALL_CHECK(cudaMemcpy(d_nr_l, h_nr_l, n_qsiz * n_derivs * sizeof(FLOAT), cudaMemcpyHostToDevice));
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     CUDA_CALL_CHECK(cudaMemset(d_cext, 0, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMemset(d_csca, 0, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMemset(d_cbak, 0, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMemset(d_g,    0, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMemset(d_pf,   0, 4 * n_qang * sizeof(FLOAT)));

     if (n_derivs > 0) {
          CUDA_CALL_CHECK(cudaMemset(d_cext_l, 0, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMemset(d_csca_l, 0, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMemset(d_cbak_l, 0, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMemset(d_g_l,    0, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMemset(d_pf_l,   0, n_derivs * 4 * n_qang * sizeof(FLOAT)));
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     fprintf(stderr, "%d %d %d %d %d %d\n", n_qsiz, n1, n2, n_qang, blocksPerGrid, threadsPerBlock);

     kernel<<<blocksPerGrid, threadsPerBlock>>>(0, n_qsiz, n_derivs, lambda, h_m, d_m_l, h_m2, d_qx, d_qw, d_nr, d_fv1, d_fv2, d_fv3, d_fv4, d_qx2, d_qw2, d_cext, d_csca, d_cbak, d_g, d_pf, d_qx_l, d_qw_l, d_nr_l, d_cext_l, d_csca_l, d_cbak_l, d_g_l, d_pf_l, n1, n2, n_qang, d_work);
     CUDA_CALL_CHECK(cudaGetLastError());


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     CUDA_CALL_CHECK(cudaMemcpy(&h_cext, d_cext, 1          * sizeof(FLOAT), cudaMemcpyDeviceToHost));
     CUDA_CALL_CHECK(cudaMemcpy(&h_csca, d_csca, 1          * sizeof(FLOAT), cudaMemcpyDeviceToHost));
     CUDA_CALL_CHECK(cudaMemcpy(&h_cbak, d_cbak, 1          * sizeof(FLOAT), cudaMemcpyDeviceToHost));
     CUDA_CALL_CHECK(cudaMemcpy(&h_g,    d_g,    1          * sizeof(FLOAT), cudaMemcpyDeviceToHost));
     CUDA_CALL_CHECK(cudaMemcpy(h_pf,    d_pf,   4 * n_qang * sizeof(FLOAT), cudaMemcpyDeviceToHost));

     if (n_derivs > 0) {
          CUDA_CALL_CHECK(cudaMemcpy(h_cext_l, d_cext_l, n_derivs              * sizeof(FLOAT), cudaMemcpyDeviceToHost));
          CUDA_CALL_CHECK(cudaMemcpy(h_csca_l, d_csca_l, n_derivs              * sizeof(FLOAT), cudaMemcpyDeviceToHost));
          CUDA_CALL_CHECK(cudaMemcpy(h_cbak_l, d_cbak_l, n_derivs              * sizeof(FLOAT), cudaMemcpyDeviceToHost));
          CUDA_CALL_CHECK(cudaMemcpy(h_g_l,    d_g_l,    n_derivs              * sizeof(FLOAT), cudaMemcpyDeviceToHost));
          CUDA_CALL_CHECK(cudaMemcpy(h_pf_l,   d_pf_l,   n_derivs * 4 * n_qang * sizeof(FLOAT), cudaMemcpyDeviceToHost));
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     *cext = h_cext;
     *csca = h_csca;
     *cbak = h_cbak;
     *g    = h_g;

     for (j = 0; j < n_qang; ++j)
          pf[0][j] = h_pf[0 * n_qang + j];
     for (j = 0; j < n_qang; ++j)
          pf[1][j] = h_pf[1 * n_qang + j];
     for (j = 0; j < n_qang; ++j)
          pf[2][j] = h_pf[2 * n_qang + j];
     for (j = 0; j < n_qang; ++j)
          pf[3][j] = h_pf[3 * n_qang + j];

     for (i = 0; i < n_derivs; ++i) {
          cext_l[i] = h_cext_l[i];
          csca_l[i] = h_csca_l[i];
          cbak_l[i] = h_cbak_l[i];
          g_l   [i] = h_g_l   [i];

          for (k = 0; k < n_qang; ++k)
               pf_l[i][0][k] = h_pf_l[i * 4 * n_qang + 0 * n_qang + k];
          for (k = 0; k < n_qang; ++k)
               pf_l[i][1][k] = h_pf_l[i * 4 * n_qang + 1 * n_qang + k];
          for (k = 0; k < n_qang; ++k)
               pf_l[i][2][k] = h_pf_l[i * 4 * n_qang + 2 * n_qang + k];
          for (k = 0; k < n_qang; ++k)
               pf_l[i][3][k] = h_pf_l[i * 4 * n_qang + 3 * n_qang + k];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free(h_qx);
     free(h_qw);
     free(h_nr);
     free(h_qx2);
     free(h_qw2);
     free(h_fv1);
     free(h_fv2);
     free(h_fv3);
     free(h_fv4);
     free(h_pf);

     cudaFree(d_qx);
     cudaFree(d_qw);
     cudaFree(d_nr);
     cudaFree(d_qx2);
     cudaFree(d_qw2);
     cudaFree(d_fv1);
     cudaFree(d_fv2);
     cudaFree(d_fv4);
     cudaFree(d_pf);

     if (n_derivs > 0) {
          free(h_qx_l);
          free(h_qw_l);
          free(h_nr_l);
          free(h_pf_l);

          cudaFree(d_qx_l);
          cudaFree(d_qw_l);
          cudaFree(d_nr_l);
          cudaFree(d_pf_l);
     }

     cudaFree(d_work);


     return 0;
}
