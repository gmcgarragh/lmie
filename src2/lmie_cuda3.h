/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#define NUMBER_OF_CORES 112


__constant__ FLOAT fv1[65536 / 4 / 4];
__constant__ FLOAT fv2[65536 / 4 / 4];
__constant__ FLOAT fv3[65536 / 4 / 4];
__constant__ FLOAT fv4[65536 / 4 / 4];


/*******************************************************************************
 *
 ******************************************************************************/
__device__ __host__ int kernel_work_n(int n_derivs, int n1, int n2) {

     int n   = 0;
     int max = 0;
if (0) {
     n += (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* + fv1    */
     n += (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* + fv2    */
     n += (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* + fv3    */
     n += (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* + fv4    */
}
     n += (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* + psi    */
     if (n_derivs > 0) {
          n += (n1 + 1) * n_derivs * sizeof(FLOAT);   if (n > max) max = n;	/* + psi_l  */
     }

     n += (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* + rn     */
     if (n_derivs > 0) {
          n += (n1 + 1) * n_derivs * sizeof(FLOAT);   if (n > max) max = n;	/* + rn_l   */
     }

     n -= (n1 + 1) * sizeof(FLOAT);                   if (n > max) max = n;	/* - rn     */
     if (n_derivs > 0) {
          n -= (n1 + 1) * n_derivs * sizeof(FLOAT);   if (n > max) max = n;	/* - rn_l   */
     }

     n += (n1 + 1) * sizeof(COMPLEX);                 if (n > max) max = n;	/* + zeta   */
     if (n_derivs > 0) {
          n += (n1 + 1) * n_derivs * sizeof(COMPLEX); if (n > max) max = n;	/* + zeta_l */
     }
if (0) {
     n += (n2 + 1) * sizeof(COMPLEX);                 if (n > max) max = n;	/* + D      */
     if (n_derivs > 0) {
          n += (n2 + 1) * n_derivs * sizeof(COMPLEX); if (n > max) max = n;	/* + D_l    */
     }
}
     n += (n1 + 1) * sizeof(COMPLEX);                 if (n > max) max = n;	/* + D      */
     if (n_derivs > 0) {
          n += (n1 + 1) * n_derivs * sizeof(COMPLEX); if (n > max) max = n;	/* + D_l    */
     }

     n += (n1 + 1) * sizeof(COMPLEX);                 if (n > max) max = n;	/* + a      */
     n += (n1 + 1) * sizeof(COMPLEX);                 if (n > max) max = n;	/* + b      */

     if (n_derivs > 0) {
          n += (n1 + 1) * n_derivs * sizeof(COMPLEX); if (n > max) max = n;	/* + a_l    */
          n += (n1 + 1) * n_derivs * sizeof(COMPLEX); if (n > max) max = n;	/* + b_l    */
     }

     return max;
}



/*******************************************************************************
 *
 ******************************************************************************/
__device__ __host__ void thread_schedule(int n_x, int threadsPerBlock, int threadIdx_x, int *i_x_for_thread, int *n_x_for_thread) {

     int n_x_per_thread_w;
     int n_x_per_thread_r;

     n_x_per_thread_w = n_x / threadsPerBlock;
     n_x_per_thread_r = n_x % threadsPerBlock;

     *n_x_for_thread = n_x_per_thread_w;
     if (threadIdx_x < n_x_per_thread_r)
          (*n_x_for_thread)++;

     *i_x_for_thread = threadIdx_x * n_x_per_thread_w + MIN(threadIdx_x, n_x_per_thread_r);
}



/*******************************************************************************
 *
 ******************************************************************************/
__global__ void kernel(int n_qsiz, int n_derivs, FLOAT lambda,
                       COMPLEX m, COMPLEX *m_l, COMPLEX m2,
                       FLOAT *qx, FLOAT *qw, FLOAT *nr,
                       FLOAT *qx2, FLOAT *qw2,
                       FLOAT *cext, FLOAT *csca, FLOAT *cbak,
                       FLOAT *g, FLOAT *pf,
                       FLOAT *qx_l, FLOAT *qw_l, FLOAT *nr_l,
                       FLOAT *cext_l, FLOAT *csca_l, FLOAT *cbak_l,
                       FLOAT *g_l, FLOAT *pf_l,
                       int n1, int n2, int n_qang,
                       int threadsPerBlock, int blocksPerSize, int blocksPerGrid) {
/*
                       int threadsPerBlock, int blocksPerSize, int blocksPerGrid, unsigned char *work) {
*/
     int i;
     int j;
     int jj;
     int k;
     int l;

     int i_work;

     int i_x;
     int n_x;

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

     __shared__ FLOAT x;
     __shared__ FLOAT x_l[MAX_DERIVS];

     __shared__ FLOAT f;
     __shared__ FLOAT f2;
     __shared__ FLOAT f_l[MAX_DERIVS];

     __shared__ FLOAT cosx;
     __shared__ FLOAT sinx;
/*
     __shared__ FLOAT *fv1;
     __shared__ FLOAT *fv2;
     __shared__ FLOAT *fv3;
     __shared__ FLOAT *fv4;
*/
     __shared__ FLOAT *rn;
     __shared__ FLOAT *rn_l;
     __shared__ FLOAT *psi;
     __shared__ FLOAT *psi_l;

     FLOAT chin0;
     __shared__ FLOAT chin0_l[MAX_DERIVS];
     FLOAT chin1;
     __shared__ FLOAT chin1_l[MAX_DERIVS];
     FLOAT chin2;
     __shared__ FLOAT chin2_l[MAX_DERIVS];
__shared__ FLOAT *sum1;
__shared__ FLOAT *sum1_l;
__shared__ FLOAT *sum2;
__shared__ FLOAT *sum2_l;
     FLOAT pin0;
     FLOAT pin1;
     FLOAT taun;

     FLOAT sr1;
     FLOAT sr1_l;
     FLOAT si1;
     FLOAT si1_l;
     FLOAT sr2;
     FLOAT sr2_l;
     FLOAT si2;
     FLOAT si2_l;

     COMPLEX c1;
     COMPLEX c2;
     COMPLEX c3;
     COMPLEX c4;
     COMPLEX c5;
     COMPLEX c6;

     __shared__ COMPLEX c1_l[MAX_DERIVS];

     __shared__ COMPLEX c6_l[MAX_DERIVS];

     __shared__ COMPLEX z;
     __shared__ COMPLEX z_l[MAX_DERIVS];

     __shared__ COMPLEX *zeta;
     __shared__ COMPLEX *zeta_l;

     __shared__ COMPLEX *D;
     __shared__ COMPLEX *D_l;

     __shared__ COMPLEX *a;
     __shared__ COMPLEX *a_l;
     __shared__ COMPLEX *b;
     __shared__ COMPLEX *b_l;
__shared__ COMPLEX *sum3;
__shared__ COMPLEX *sum3_l;
     __shared__ COMPLEX *sum;
     __shared__ COMPLEX *sum_l;
     __shared__ COMPLEX *dif;
     __shared__ COMPLEX *dif_l;
/*
     COMPLEX sp;
     COMPLEX sp_l[MAX_DERIVS];
     COMPLEX sm;
     COMPLEX sm_l[MAX_DERIVS];
*/
     COMPLEX sp0;
     COMPLEX sp1;
     COMPLEX sp0_l[MAX_DERIVS];
     COMPLEX sp1_l[MAX_DERIVS];
     COMPLEX sm0;
     COMPLEX sm1;
     COMPLEX sm0_l[MAX_DERIVS];
     COMPLEX sm1_l[MAX_DERIVS];

     COMPLEX s1;
     COMPLEX s1_l;
     COMPLEX s2;
     COMPLEX s2_l;

     extern __shared__ __align__(4) unsigned char work[];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {

     i_work = 0;
/*
     i_work = blockIdx.x % NUMBER_OF_CORES * kernel_work_n(n_derivs, n1, n2);

     fv1 = (FLOAT *) (work + i_work); i_work += (n1 + 1) * sizeof(FLOAT);
     fv2 = (FLOAT *) (work + i_work); i_work += (n1 + 1) * sizeof(FLOAT);
     fv3 = (FLOAT *) (work + i_work); i_work += (n1 + 1) * sizeof(FLOAT);
     fv4 = (FLOAT *) (work + i_work); i_work += (n1 + 1) * sizeof(FLOAT);
*/
     psi = (FLOAT *) (work + i_work); i_work += (n1 + 1) * sizeof(FLOAT);
     if (n_derivs > 0) {
          psi_l = (FLOAT *) (work + i_work); i_work += (n1 + 1) * n_derivs * sizeof(FLOAT);
     }

     rn = (FLOAT *) (work + i_work); i_work += (n1 + 1) * sizeof(FLOAT);
     if (n_derivs > 0) {
          rn_l = (FLOAT *) (work + i_work); i_work += (n1 + 1) * n_derivs * sizeof(FLOAT);
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
     i = blockIdx.x / blocksPerSize;

     x = 2. * PI * qx[i] / lambda;
     z = m * x;

     f = nr[i] * qw[i];

     f2 = f * .5;

     for (j = 0; j < n_derivs; ++j) {
          x_l[j] = 2. * PI * qx_l[A2(i,j)] / lambda;
/*
          x_l[j] = 2. * PI * (qx_l[A2(i,j)] -
                   qx[i] * lambda_l[j] / lambda) / lambda;
*/
          z_l[j] = m_l[j] * x + m * x_l[j];

          f_l[j] = nr_l[A2(i,j)] * qw[i] + nr[i] * qw_l[A2(i,j)];
     }


     cosx = cos(x);
     sinx = sin(x);
}

__syncthreads();

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (0) {
/*
if (threadIdx.x == 0) {
     for (i = 1; i <= n1; ++i) {
*/
     thread_schedule(n1, threadsPerBlock, threadIdx.x, &i_x, &n_x);

     for (i = i_x + 1; i < i_x + n_x + 1; ++i) {
          fv1[i] = 2. * i + 1.;
          fv2[i] = (i + 1.) / i;
          fv3[i] = (i + 2.) / fv2[i];
          fv4[i] = fv1[i] / (i * (i + 1.));
     }
/*
}
*/
}
__syncthreads();

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n1 =               x   + 4.05f * __powf(              x,   1.f/3.f) + 8;
     n2 = MAX(n1, absxc(z)) + 6.40f * __powf(MAX(n1, absxc(z)), 1.f/3.f) + 8;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
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

          for (k = 0; k < n_derivs; ++k) {
               psi_l[A2(j+1,k)]  = rn_l[A2(j+1,k)] * psi[j] + rn[j+1] * psi_l[A2(j,k)];
          }
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
     i_work -= (n1 + 1) * sizeof(FLOAT);
     if (n_derivs > 0)
          i_work -= (n1 + 1) * n_derivs * sizeof(FLOAT);

     zeta = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * sizeof(COMPLEX);

     if (n_derivs > 0) {
          zeta_l = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * n_derivs * sizeof(COMPLEX);
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
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
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
/*
     D = (COMPLEX *) (work + i_work); i_work += (n2 + 1) * sizeof(COMPLEX);

     if (n_derivs > 0) {
          D_l = (COMPLEX *) (work + i_work); i_work += (n2 + 1) * n_derivs * sizeof(COMPLEX);
     }
*/
     D = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * sizeof(COMPLEX);

     if (n_derivs > 0) {
          D_l = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * n_derivs * sizeof(COMPLEX);
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
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
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (threadIdx.x == 0) {
     a = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * sizeof(COMPLEX);
     b = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * sizeof(COMPLEX);

     if (n_derivs > 0) {
          a_l = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * n_derivs * sizeof(COMPLEX);
          b_l = (COMPLEX *) (work + i_work); i_work += (n1 + 1) * n_derivs * sizeof(COMPLEX);
     }
}

__syncthreads();

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
if (threadIdx.x == 0) {
     for (j = 1; j <= n1; ++j) {
*/
     thread_schedule(n1, threadsPerBlock, threadIdx.x, &i_x, &n_x);

     for (j = i_x + 1; j < i_x + n_x + 1; ++j) {
          a1 = j / x;

          for (k = 0; k < n_derivs; ++k)
               a1_l[k] = -a1 * x_l[k] / x;

          c1 = D[j] / m + a1;

          c2 = c1 * psi [j] - psi [j-1];
          c3 = c1 * zeta[j] - zeta[j-1];

          a[j] = c2 / c3;

          if (n_derivs > 0) {
               c4 = D[j] / m2;
/*
               c5 = c3 * c3;
*/
          }

          for (k = 0; k < n_derivs; ++k) {
               c6 = D_l[A2(j,k)] / m - c4 * m_l[k] + a1_l[k];
/*
               a_l[A2(j,k)] = (c6 * psi [j] + c1 * psi_l [A2(j,k)] - psi_l [A2(j-1,k)]) / c3 - c2 *
                              (c6 * zeta[j] + c1 * zeta_l[A2(j,k)] - zeta_l[A2(j-1,k)]) / c5;
*/
               a_l[A2(j,k)] = ((c6 * psi [j] + c1 * psi_l [A2(j,k)] - psi_l [A2(j-1,k)]) - a[j] *
                               (c6 * zeta[j] + c1 * zeta_l[A2(j,k)] - zeta_l[A2(j-1,k)])) / c3;
          }

          c1 = m * D[j] + a1;

          c2 = c1 * psi [j] - psi [j-1];
          c3 = c1 * zeta[j] - zeta[j-1];

          b[j] = c2 / c3;
/*
          if (n_derivs > 0)
               c5 = c3 * c3;
*/
          for (k = 0; k < n_derivs; ++k) {
               c4 = m_l[k] * D[j] + m * D_l[A2(j,k)] + a1_l[k];
/*
               b_l[A2(j,k)] = (c4 * psi [j] + c1 * psi_l [A2(j,k)] - psi_l [A2(j-1,k)]) / c3 - c2 *
                              (c4 * zeta[j] + c1 * zeta_l[A2(j,k)] - zeta_l[A2(j-1,k)]) / c5;
*/
               b_l[A2(j,k)] = ((c4 * psi [j] + c1 * psi_l [A2(j,k)] - psi_l [A2(j-1,k)]) - b[j] *
                               (c4 * zeta[j] + c1 * zeta_l[A2(j,k)] - zeta_l[A2(j-1,k)])) / c3;
          }
     }
/*
}
*/
__syncthreads();

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
if (blockIdx.x % blocksPerSize == 0 && threadIdx.x == 0) {
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
}

__syncthreads();
*/
     sum1 = (FLOAT *) psi;
     sum2 = (FLOAT *) zeta;
     sum3 = (COMPLEX *) D;

     if (n_derivs > 0) {
          sum1_l = (FLOAT *) psi_l;
          sum2_l = (FLOAT *) zeta_l;
          sum3_l = (COMPLEX *) D_l;
     }

     thread_schedule(n1, threadsPerBlock, threadIdx.x, &i_x, &n_x);

     a3 = pow(-1., i_x + 1);
     for (j = i_x + 1; j < i_x + n_x + 1; ++j) {
          sum1[j] = fv1[j] * (realxc(a[j]) + realxc(b[j]));

          sum2[j] = fv1[j] * (realxc(a[j] * ~a[j]) + realxc(b[j] * ~b[j]));

          sum3[j] = fv1[j] * (a[j] - b[j]) * a3;

          for (k = 0; k < n_derivs; ++k) {
               sum1_l[A2(j,k)] = fv1[j] * (realxc(a_l[A2(j,k)]) + realxc(b_l[A2(j,k)]));

               sum2_l[A2(j,k)] = fv1[j] * (2. * (realxc(a[j] * ~a_l[A2(j,k)]) + realxc(b[j] * ~b_l[A2(j,k)])));

               sum3_l[A2(j,k)] = fv1[j] * (a_l[A2(j,k)] - b_l[A2(j,k)]) * a3;
          }

          a3 = -a3;
     }

__syncthreads();

if (threadIdx.x == 0) {
     a1 = 0.;
     a2 = 0.;
     c1 = 0.;

     for (k = 0; k < n_derivs; ++k) {
          a1_l[k] = 0.;
          a2_l[k] = 0.;
          c1_l[k] = 0.;
     }
     for (j = 1; j < n1; ++j) {
          a1 += sum1[j];
          a2 += sum2[j];
          c1 += sum3[j];

          for (k = 0; k < n_derivs; ++k) {
               a1_l[k] += sum1_l[A2(j,k)];
               a2_l[k] += sum2_l[A2(j,k)];
               c1_l[k] += sum3_l[A2(j,k)];
          }
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
          cext_l[j] += f_l[j] * a1 + f * a1_l[j];
          csca_l[j] += f_l[j] * a2 + f * a2_l[j];
          cbak_l[j] += f_l[j] * realxc(c1 * ~c1) + f * 2. * realxc(c1 * ~c1_l[j]);
#else
          atomicAdd_x(&cext_l[j], f_l[j] * a1 + f * a1_l[j]);
          atomicAdd_x(&csca_l[j], f_l[j] * a2 + f * a2_l[j]);
          atomicAdd_x(&cbak_l[j], f_l[j] * realxc(c1 * ~c1) + f * 2. * realxc(c1 * ~c1_l[j]));
#endif
     }
}

__syncthreads();

     /*--------------------------------------------------------------------
      *
      *------------------------------------------------------------------*/
/*
if (blockIdx.x % blocksPerSize == 0 && threadIdx.x == 0) {
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
}
*/
     thread_schedule(n1 - 1, threadsPerBlock, threadIdx.x, &i_x, &n_x);

     for (j = i_x + 1; j < i_x + n_x + 1; ++j) {
          sum1[j] = fv3[j] * realxc(a[j] * ~a[j + 1] + b[j] * ~b[j + 1]) +
                    fv4[j] * realxc(a[j] * ~b[j]);

          for (k = 0; k < n_derivs; ++k) {
               sum1_l[k] = fv3[j] * realxc(a_l[A2(j,k)] * ~a[j + 1] + a[j] * ~a_l[A2(j + 1,k)] +
                                           b_l[A2(j,k)] * ~b[j + 1] + b[j] * ~b_l[A2(j + 1,k)]) +
                           fv4[j] * realxc(a_l[A2(j,k)] * ~b[j] + a[j] * ~b_l[A2(j,k)]);
          }
     }

__syncthreads();

if (threadIdx.x == 0) {
     a1 = 0.;

     for (k = 0; k < n_derivs; ++k)
          a1_l[k] = 0.;

     for (j = 1; j < n1; ++j) {
          a1 += sum1[j];

          for (k = 0; k < n_derivs; ++k)
               a1_l[k] += sum1_l[A2(j,k)];
     }
#ifndef USE_ATOMIC
     *g += f * a1;
#else
     atomicAdd_x(g, f * a1);
#endif
     for (j = 0; j < n_derivs; ++j)
#ifndef USE_ATOMIC
          g_l[j] += f_l[j] * a1 + f * a1_l[j];
#else
          atomicAdd_x(&g_l[j], f_l[j] * a1 + f * a1_l[j]);
#endif
}

__syncthreads();

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
if (threadIdx.x == 0) {
*/
     sum = a;
     dif = b;

     if (n_derivs > 0) {
          sum_l = a_l;
          dif_l = b_l;
     }
/*
     for (j = 1; j <= n1; ++j) {
*/
     thread_schedule(n1, threadsPerBlock, threadIdx.x, &i_x, &n_x);

     for (j = i_x + 1; j < i_x + n_x + 1; ++j) {

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
/*
}
*/
__syncthreads();

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     j = blockIdx.x % blocksPerSize * threadsPerBlock + threadIdx.x;

     if (j < n_qang) {
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
*/
/*
     j = blockIdx.x % blocksPerSize * threadsPerBlock + threadIdx.x;

     if (j < n_qang / 2) {
*/
     thread_schedule(n_qang / 2, threadsPerBlock, threadIdx.x, &i_x, &n_x);

     for (j = i_x; j < i_x + n_x; ++j) {
          sp0 = 0.;
          sm0 = 0.;
          sp1 = 0.;
          sm1 = 0.;

          for (k = 0; k < n_derivs; ++k) {
               sp0_l[k] = 0.;
               sm0_l[k] = 0.;
               sp1_l[k] = 0.;
               sm1_l[k] = 0.;
          }

          a5 = qx2[j];

          a6 = 1.;

          pin0 = 0.;
          pin1 = 1.;
if (n_derivs == 0) {
          for (k = 1; k <= n1; ++k) {
               a1 = a5 * pin1;
               a2 = a1 - pin0;

               taun = k * a2 - pin0;

               a3 = pin1 + taun;
               a4 = pin1 - taun;

               sp0 += a3 * sum[k];
               sm0 += a4 * dif[k];

               sp1 += a4 * sum[k] * a6;
               sm1 += a3 * dif[k] * a6;

               pin0 = pin1;

               pin1 = a1 + fv2[k] * a2;

               a6 = -a6;
          }
}
else {
          for (k = 1; k <= n1; ++k) {
               a1 = a5 * pin1;
               a2 = a1 - pin0;

               taun = k * a2 - pin0;

               a3 = pin1 + taun;
               a4 = pin1 - taun;

               sp0 += a3 * sum[k];
               sm0 += a4 * dif[k];

               sp1 += a4 * sum[k] * a6;
               sm1 += a3 * dif[k] * a6;

               for (l = 0; l < n_derivs; ++l) {
                    sp0_l[l] += a3 * sum_l[A2(k,l)];
                    sm0_l[l] += a4 * dif_l[A2(k,l)];

                    sp1_l[l] += a4 * sum_l[A2(k,l)] * a6;
                    sm1_l[l] += a3 * dif_l[A2(k,l)] * a6;
               }

               pin0 = pin1;

               pin1 = a1 + fv2[k] * a2;

               a6 = -a6;
          }
}
     l = 0;
     jj = l == 0 ? j : n_qang - j - 1;
          s1 = .5f * (sp0 + sm0);
          s2 = .5f * (sp0 - sm0);

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
          pf[0 * n_qang + jj] += f * a5;
          pf[1 * n_qang + jj] += f * a6;
          pf[2 * n_qang + jj] += f * a7;
          pf[3 * n_qang + jj] += f * a8;
#else
          atomicAdd_x(&pf[0 * n_qang + jj], f * a5);
          atomicAdd_x(&pf[1 * n_qang + jj], f * a6);
          atomicAdd_x(&pf[2 * n_qang + jj], f * a7);
          atomicAdd_x(&pf[3 * n_qang + jj], f * a8);
#endif
          for (k = 0; k < n_derivs; ++k) {
               s1_l = .5 * (sp0_l[k] + sm0_l[k]);
               s2_l = .5 * (sp0_l[k] - sm0_l[k]);

               sr1_l = realxc(s1_l);
               si1_l = imagxc(s1_l);
               sr2_l = realxc(s2_l);
               si2_l = imagxc(s2_l);

               a1 = 2. * (sr1_l * sr1 + si1_l * si1);
               c2 = 2. * s1_l * ~s2;
               a2 = realxc(c2);
               a3 = 2. * (sr2_l * sr2 + si2_l * si2);
               c4 = 2. * s2_l * ~s1;
               a4 = realxc(c4);
#ifndef USE_ATOMIC

#else
               atomicAdd_x(&pf_l[k * 4 * n_qang + 0 * n_qang + jj], f_l[k] * a5 + f2 *         (a1 + a3));
               atomicAdd_x(&pf_l[k * 4 * n_qang + 1 * n_qang + jj], f_l[k] * a6 + f2 *         (a2 + a4));
               atomicAdd_x(&pf_l[k * 4 * n_qang + 2 * n_qang + jj], f_l[k] * a7 + f2 * -       (a1 - a3));
               atomicAdd_x(&pf_l[k * 4 * n_qang + 3 * n_qang + jj], f_l[k] * a8 + f2 * -realxc((c2 - c4) * MAKE_COMPLEX(0., 1.f)));
#endif
          }

     l = 1;
     jj = l == 0 ? j : n_qang - j - 1;
          s1 = .5f * (sp1 + sm1);
          s2 = .5f * (sp1 - sm1);

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
          pf[0 * n_qang + jj] += f * a5;
          pf[1 * n_qang + jj] += f * a6;
          pf[2 * n_qang + jj] += f * a7;
          pf[3 * n_qang + jj] += f * a8;
#else
          atomicAdd_x(&pf[0 * n_qang + jj], f * a5);
          atomicAdd_x(&pf[1 * n_qang + jj], f * a6);
          atomicAdd_x(&pf[2 * n_qang + jj], f * a7);
          atomicAdd_x(&pf[3 * n_qang + jj], f * a8);
#endif
          for (k = 0; k < n_derivs; ++k) {
               s1_l = .5 * (sp1_l[k] + sm1_l[k]);
               s2_l = .5 * (sp1_l[k] - sm1_l[k]);

               sr1_l = realxc(s1_l);
               si1_l = imagxc(s1_l);
               sr2_l = realxc(s2_l);
               si2_l = imagxc(s2_l);

               a1 = 2. * (sr1_l * sr1 + si1_l * si1);
               c2 = 2. * s1_l * ~s2;
               a2 = realxc(c2);
               a3 = 2. * (sr2_l * sr2 + si2_l * si2);
               c4 = 2. * s2_l * ~s1;
               a4 = realxc(c4);
#ifndef USE_ATOMIC

#else
               atomicAdd_x(&pf_l[k * 4 * n_qang + 0 * n_qang + jj], f_l[k] * a5 + f2 *         (a1 + a3));
               atomicAdd_x(&pf_l[k * 4 * n_qang + 1 * n_qang + jj], f_l[k] * a6 + f2 *         (a2 + a4));
               atomicAdd_x(&pf_l[k * 4 * n_qang + 2 * n_qang + jj], f_l[k] * a7 + f2 * -       (a1 - a3));
               atomicAdd_x(&pf_l[k * 4 * n_qang + 3 * n_qang + jj], f_l[k] * a8 + f2 * -realxc((c2 - c4) * MAKE_COMPLEX(0., 1.f)));
#endif
          }
     }


     __syncthreads();
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
     int blocksPerSize;
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

     size_t work_size;
/*
     unsigned char *d_work;
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     threadsPerBlock = THREADSPERBLOCK;
if (0) {
/*
     blocksPerSize   = n_qang /     threadsPerBlock + 1;

     if (n_qang % blocksPerSize == 0)
          threadsPerBlock = n_qang / blocksPerSize;
     else
          threadsPerBlock = n_qang / blocksPerSize + 1;
*/
     blocksPerSize   = n_qang / 2 / threadsPerBlock + 1;
/*
     if (n_qang / 2 % blocksPerSize == 0)
          threadsPerBlock = n_qang / 2 / blocksPerSize;
     else
          threadsPerBlock = n_qang / 2 / blocksPerSize + 1;
*/
     blocksPerGrid   = n_qsiz * blocksPerSize;
}
else {
     blocksPerSize   = 1;
     blocksPerGrid   = n_qsiz;
}

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

     CUDA_CALL_CHECK(cudaMalloc((void **) &d_cext, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_csca, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_cbak, 1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_g,    1          * sizeof(FLOAT)));
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_pf,   4 * n_qang * sizeof(FLOAT)));

     if (n_derivs > 0) {
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_m_l,           n_derivs * sizeof(COMPLEX)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_qx_l, n_qsiz * n_derivs * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_qw_l, n_qsiz * n_derivs * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_nr_l, n_qsiz * n_derivs * sizeof(FLOAT)));

          CUDA_CALL_CHECK(cudaMalloc((void **) &d_cext_l, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_csca_l, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_cbak_l, n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_g_l,    n_derivs              * sizeof(FLOAT)));
          CUDA_CALL_CHECK(cudaMalloc((void **) &d_pf_l,   n_derivs * 4 * n_qang * sizeof(FLOAT)));
     }
/*
     CUDA_CALL_CHECK(cudaMalloc((void **) &d_work, NUMBER_OF_CORES * kernel_work_n(n_derivs, n1, n2)));
*/

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
/*
     int n;
     int r;
     int ii;
     int jj;
     int i_offset;
     int j_offset;

     n = n_qsiz / 112 + 1;
     r = n_qsiz % 112;

     i_offset = 0;
     for (i = 0; i < n; ++i) {
          j_offset = 0;
          for (j = 0; j < 112; ++j) {
               ii = i_offset + j;
               jj = j_offset + i;

               if (ii < n_qsiz) {
                    h_qx[jj] = qx[ii];
                    h_qw[jj] = qw[ii];
                    h_nr[jj] = nr[ii];

                    for (k = 0; k < n_derivs; ++k) {
                         h_qx_l[A2(jj,k)] = qx_l[ii][k];
                         h_qw_l[A2(jj,k)] = qw_l[ii][k];
                         h_nr_l[A2(jj,k)] = nr_l[ii][k];
                    }
               }

               if (j < r)
                    j_offset += n;
               else
                    j_offset += n - 1;
          }

          i_offset += 112;
     }
*/
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

     CUDA_CALL_CHECK(cudaMemcpyToSymbol("fv1", h_fv1, (n1 + 1) * sizeof(FLOAT), 0, cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpyToSymbol("fv2", h_fv2, (n1 + 1) * sizeof(FLOAT), 0, cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpyToSymbol("fv3", h_fv3, (n1 + 1) * sizeof(FLOAT), 0, cudaMemcpyHostToDevice));
     CUDA_CALL_CHECK(cudaMemcpyToSymbol("fv4", h_fv4, (n1 + 1) * sizeof(FLOAT), 0, cudaMemcpyHostToDevice));

     if (n_derivs > 0) {
          CUDA_CALL_CHECK(cudaMemcpy(d_m_l,  h_m_l,           n_derivs * sizeof(COMPLEX), cudaMemcpyHostToDevice));
          CUDA_CALL_CHECK(cudaMemcpy(d_qx_l, h_qx_l, n_qsiz * n_derivs * sizeof(FLOAT),   cudaMemcpyHostToDevice));
          CUDA_CALL_CHECK(cudaMemcpy(d_qw_l, h_qw_l, n_qsiz * n_derivs * sizeof(FLOAT),   cudaMemcpyHostToDevice));
          CUDA_CALL_CHECK(cudaMemcpy(d_nr_l, h_nr_l, n_qsiz * n_derivs * sizeof(FLOAT),   cudaMemcpyHostToDevice));
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
/*
     work_size = 0;
*/
     work_size = kernel_work_n(n_derivs, n1, n2);

     fprintf(stderr, "%d %d %d %d %d %d %d %d\n", n_qsiz, n1, n2, n_qang, blocksPerGrid, blocksPerSize, threadsPerBlock, work_size);

     kernel<<<blocksPerGrid, threadsPerBlock, work_size>>>(n_qsiz, n_derivs, lambda, h_m, d_m_l, h_m2, d_qx, d_qw, d_nr, d_qx2, d_qw2, d_cext, d_csca, d_cbak, d_g, d_pf, d_qx_l, d_qw_l, d_nr_l, d_cext_l, d_csca_l, d_cbak_l, d_g_l, d_pf_l, n1, n2, n_qang, threadsPerBlock, blocksPerSize, blocksPerGrid);
/*
     kernel<<<blocksPerGrid, threadsPerBlock, work_size>>>(n_qsiz, n_derivs, lambda, h_m, d_m_l, h_m2, d_qx, d_qw, d_nr, d_qx2, d_qw2, d_cext, d_csca, d_cbak, d_g, d_pf, d_qx_l, d_qw_l, d_nr_l, d_cext_l, d_csca_l, d_cbak_l, d_g_l, d_pf_l, n1, n2, n_qang, threadsPerBlock, blocksPerSize, blocksPerGrid, d_work);
*/
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
/*
     cudaFree(d_work);
*/

     return 0;
}
