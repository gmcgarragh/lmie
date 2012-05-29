/******************************************************************************%
**
**    Copyright (C) 2008-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#ifdef USE_MPI
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_PTHREADS
#include <pthread.h>
#endif

#include <rtutil_math.h>
#include <rtutil_scat_int.h>

#include "lmie_core.h"
#include "lmie_schedule.h"


/*******************************************************************************
 *
 ******************************************************************************/
#if defined USE_OPENMP || defined USE_PTHREADS
static int size_dist_loop_index(int n_qsize, int i_proc, int n_procs, int **index) {

     int i;
     int ii;

     int n_qsize2;
     int n_qsize2_w;
     int n_qsize2_r;

     if (i_proc >= n_qsize)
          return 0;

     n_qsize2_w = n_qsize / n_procs;
     n_qsize2_r = n_qsize % n_procs;

     n_qsize2 = n_qsize2_w;
     if (i_proc < n_qsize2_r)
          n_qsize2++;

     *index = alloc_array1_i(n_qsize2);

     ii = i_proc;
     for (i = 0; i < n_qsize2; ++i) {
          (*index)[i] = ii;
          ii += n_procs;
     }

     return n_qsize2;
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
static void init_mie_solution(int n_qang, int n_derivs,
                              double *cext, double *csca, double *cbak,
                              double *g, double **pf,
                              double *cext_l, double *csca_l, double *cbak_l,
                              double *g_l, double ***pf_l) {

     int i;
     int j;

     *cext = 0.;
     *csca = 0.;
     *cbak = 0.;
     *g    = 0.;

     for (i = 0; i < n_qang; ++i) {
          pf[0][i] = 0.;
          pf[1][i] = 0.;
          pf[2][i] = 0.;
          pf[3][i] = 0.;
     }

     for (i = 0; i < n_derivs; ++i) {
          cext_l[i] = 0.;
          csca_l[i] = 0.;
          cbak_l[i] = 0.;
          g_l   [i] = 0.;

          for (j = 0; j < n_qang; ++j) {
               pf_l[i][0][j] = 0.;
               pf_l[i][1][j] = 0.;
               pf_l[i][2][j] = 0.;
               pf_l[i][3][j] = 0.;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#if defined USE_MPI || defined USE_OPENMP || defined USE_PTHREADS
static void add_mie_solution(int n_qang, int n_derivs,
                             double *cext1, double *csca1, double *cbak1,
                             double *g1, double **pf1,
                             double *cext_l1, double *csca_l1, double *cbak_l1,
                             double *g_l1, double ***pf_l1,
                             double cext2, double csca2, double cbak2,
                             double g2, double **pf2,
                             double *cext_l2, double *csca_l2, double *cbak_l2,
                             double *g_l2, double ***pf_l2) {

     int i;
     int j;

     *cext1 += cext2;
     *csca1 += csca2;
     *cbak1 += cbak2;
     *g1    += g2;

     for (i = 0; i < n_qang; ++i) {
          pf1[0][i] += pf2[0][i];
          pf1[1][i] += pf2[1][i];
          pf1[2][i] += pf2[2][i];
          pf1[3][i] += pf2[3][i];
     }

     for (i = 0; i < n_derivs; ++i) {
          cext_l1[i] += cext_l2[i];
          csca_l1[i] += csca_l2[i];
          cbak_l1[i] += cbak_l2[i];
          g_l1   [i] += g_l2   [i];

          for (j = 0; j < n_qang; ++j) {
               pf_l1[i][0][j] += pf_l2[i][0][j];
               pf_l1[i][1][j] += pf_l2[i][1][j];
               pf_l1[i][2][j] += pf_l2[i][2][j];
               pf_l1[i][3][j] += pf_l2[i][3][j];
          }
     }
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
static int lmie_core_sequential(int n_qsize, int n_derivs1, int n_derivs2,
                                int *index1, int *index2,
                                double lambda, double *lambda_l,
                                dcomplex m, dcomplex *m_l, dcomplex m2,
                                double *qx, double *qw, double *nr,
                                double *fv1, double *fv2, double *fv3, double *fv4,
                                double *qx2, double *qw2,
                                double *cext, double *csca, double *cbak,
                                double *g, double **pf,
                                double **qx_l, double **qw_l, double **nr_l,
                                double *cext_l, double *csca_l, double *cbak_l,
                                double *g_l, double ***pf_l,
                                int n1, int n2, int n_qang, int n_threads) {

     lmie_core_shared_data s;

     lmie_core_threads_data d;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     lmie_core_shared_import(&s, n_qsize, n_derivs1, n_derivs2, index1, index2,
                             lambda, lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2,
                             fv3, fv4, qx2, qw2, qx_l, qw_l, nr_l, n1, n2, n_qang);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     s.i    = 0;

     s.flag = 0;

     d.s = &s;

     d.n_qsize = n_qsize;

     d.cext   = cext;
     d.csca   = csca;
     d.cbak   = cbak;
     d.g      = g;
     d.pf     = pf;
     d.cext_l = cext_l;
     d.csca_l = csca_l;
     d.cbak_l = cbak_l;
     d.g_l    = g_l;
     d.pf_l   = pf_l;


     lmie_core(&d);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#if defined USE_OPENMP || defined USE_PTHREADS

static int lmie_core_xthreads(int n_qsize, int n_derivs1, int n_derivs2,
                              int *index1, int *index2,
                              double lambda, double *lambda_l,
                              dcomplex m, dcomplex *m_l, dcomplex m2,
                              double *qx, double *qw, double *nr,
                              double *fv1, double *fv2, double *fv3, double *fv4,
                              double *qx2, double *qw2,
                              double *cext, double *csca, double *cbak,
                              double *g, double **pf,
                              double **qx_l, double **qw_l, double **nr_l,
                              double *cext_l, double *csca_l, double *cbak_l,
                              double *g_l, double ***pf_l,
                              int n1, int n2, int n_qang, int n_threads) {

     int i;

     int n_derivs;

     int n_threads2;
#ifdef USE_PTHREADS
     void *t_result;

     pthread_t *pt;

     pthread_attr_t attr;
#endif
     lmie_core_shared_data s;

     lmie_core_threads_data *d;


     n_derivs = n_derivs1 + n_derivs2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     lmie_core_shared_import(&s, n_qsize, n_derivs1, n_derivs2, index1, index2,
                             lambda, lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2,
                             fv3, fv4, qx2, qw2, qx_l, qw_l, nr_l, n1, n2, n_qang);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     s.i    = 0;

     s.flag = 1;

     d = (lmie_core_threads_data *) malloc(n_threads * sizeof(lmie_core_threads_data));

     n_threads2 = 0;
     for (i = 0; i < n_threads; ++i) {
          d[i].s = &s;

          d[i].n_qsize = size_dist_loop_index(n_qsize, i, n_threads, &d[i].index);

          if (d[i].n_qsize > 0) {
               n_threads2++;

               d[i].cext = alloc_array1_d(1);
               d[i].csca = alloc_array1_d(1);
               d[i].cbak = alloc_array1_d(1);
               d[i].g    = alloc_array1_d(1);

               d[i].pf   = alloc_array2_d(4, n_qang);

               if (n_derivs > 0) {
                    d[i].cext_l = alloc_array1_d(n_derivs);
                    d[i].csca_l = alloc_array1_d(n_derivs);
                    d[i].cbak_l = alloc_array1_d(n_derivs);
                    d[i].g_l    = alloc_array1_d(n_derivs);

                    d[i].pf_l   = alloc_array3_d(n_derivs, 4, n_qang);
               }

               init_mie_solution(n_qang, n_derivs, d[i].cext, d[i].csca,
                                 d[i].cbak, d[i].g, d[i].pf, d[i].cext_l,
                                 d[i].csca_l, d[i].cbak_l, d[i].g_l, d[i].pf_l);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifdef USE_OPENMP
     omp_set_num_threads(n_threads2);
#pragma omp parallel
{
#pragma omp for
     for (i = 0; i < n_threads2; ++i)
          lmie_core(&d[i]);
}
#else
     pt = (pthread_t *) malloc(n_threads * sizeof(pthread_t));

     pthread_attr_init(&attr);
     pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

     for (i = 0; i < n_threads2; ++i)
          pthread_create(&pt[i], &attr, (void *(*)(void *)) lmie_core, &d[i]);

     for (i = 0; i < n_threads2; ++i) {
          pthread_join  ( pt[i], &t_result);
          if (t_result == NULL) {
               printf("ERROR: pthread_join()\n");
               return -1;
          }
     }

     pthread_attr_destroy(&attr);

     free(pt);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_threads2; ++i)
          add_mie_solution(n_qang, n_derivs, cext, csca, cbak, g, pf, cext_l,
                           csca_l, cbak_l, g_l, pf_l, *d[i].cext, *d[i].csca,
                           *d[i].cbak, *d[i].g, d[i].pf, d[i].cext_l,
                           d[i].csca_l, d[i].cbak_l, d[i].g_l, d[i].pf_l);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_threads2; ++i) {
          free_array1_i(d[i].index);

          free_array1_d(d[i].cext);
          free_array1_d(d[i].csca);
          free_array1_d(d[i].cbak);
          free_array1_d(d[i].g);
          free_array2_d(d[i].pf);

          if (n_derivs > 0) {
               free_array1_d(d[i].cext_l);
               free_array1_d(d[i].csca_l);
               free_array1_d(d[i].cbak_l);
               free_array1_d(d[i].g_l);
               free_array3_d(d[i].pf_l);
          }
     }

     free(d);


     return 0;
}

#endif


/*******************************************************************************
 *
 ******************************************************************************/
#ifdef USE_CUDA
#ifdef __cplusplus
extern "C"
#endif
int lmie_core_cuda(int n_qsize, int n_derivs,
                   double lambda, double *lambda_l,
                   dcomplex m, dcomplex *m_l, dcomplex m2,
                   double *qx, double *qw, double *nr,
                   double *fv1, double *fv2, double *fv3, double *fv4,
                   double *qx2, double *qw2,
                   double *cext, double *csca, double *cbak,
                   double *g, double **pf,
                   double **qx_l, double **qw_l, double **nr_l,
                   double *cext_l, double *csca_l, double *cbak_l,
                   double *g_l, double ***pf_l,
                   int n1, int n2, int n_qang);
#endif

static int lmie_core_shared(int n_qsize, int n_derivs1, int n_derivs2,
                            int *index1, int *index2,
                            double lambda, double *lambda_l,
                            dcomplex m, dcomplex *m_l, dcomplex m2,
                            double *qx, double *qw, double *nr,
                            double *fv1, double *fv2, double *fv3, double *fv4,
                            double *qx2, double *qw2,
                            double *cext, double *csca, double *cbak,
                            double *g, double **pf,
                            double **qx_l, double **qw_l, double **nr_l,
                            double *cext_l, double *csca_l, double *cbak_l,
                            double *g_l, double ***pf_l,
                            int n1, int n2, int n_qang, int n_threads) {

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_mie_solution(n_qang, n_derivs1 + n_derivs2, cext, csca, cbak, g, pf,
                       cext_l, csca_l, cbak_l, g_l, pf_l);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifdef USE_CUDA
if (0) {
#endif
#if defined USE_OPENMP || defined USE_PTHREADS
     if (n_threads == 1) {
#endif
          if (lmie_core_sequential
                        (n_qsize, n_derivs1, n_derivs2, index1, index2, lambda,
                         lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2, fv3, fv4,
                         qx2, qw2, cext, csca, cbak, g, pf, qx_l, qw_l, nr_l,
                         cext_l, csca_l, cbak_l, g_l, pf_l, n1, n2, n_qang,
                         n_threads)) {
               eprintf("ERROR: lmie_core_sequential()\n");
               return -1;
          }
#if defined USE_OPENMP || defined USE_PTHREADS
     }
     else {
          if (lmie_core_xthreads
                        (n_qsize, n_derivs1, n_derivs2, index1, index2, lambda,
                         lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2, fv3, fv4,
                         qx2, qw2, cext, csca, cbak, g, pf, qx_l, qw_l, nr_l,
                         cext_l, csca_l, cbak_l, g_l, pf_l, n1, n2, n_qang,
                         n_threads)) {
               eprintf("ERROR: lmie_core_xthreads()\n");
               return -1;
          }
     }
#endif
#ifdef USE_CUDA
}
          if (lmie_core_cuda
                        (n_qsize, n_derivs1, lambda, lambda_l, m, m_l, m2, qx, qw,
                         nr, fv1, fv2, fv3, fv4, qx2, qw2, cext, csca, cbak, g,
                         pf, qx_l, qw_l, nr_l, cext_l, csca_l, cbak_l, g_l, pf_l,
                         n1, n2, n_qang)) {
               eprintf("ERROR: lmie_core_cuda()\n");
               return -1;
          }
#endif
     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef USE_MPI

static int lmie_core_master(int n_qsize, int n_derivs1, int n_derivs2,
                            int *index1, int *index2,
                            double lambda, double *lambda_l,
                            dcomplex m, dcomplex *m_l, dcomplex m2,
                            double *qx, double *qw, double *nr,
                            double *fv1, double *fv2, double *fv3, double *fv4,
                            double *qx3, double *qw3,
                            double *cext, double *csca, double *cbak,
                            double *g, double **pf,
                            double **qx_l, double **qw_l, double **nr_l,
                            double *cext_l, double *csca_l, double *cbak_l,
                            double *g_l, double ***pf_l,
                            int n1, int n2, int n_qang, int verbose, int n_threads) {

     int i;
     int j;
     int jj;
     int k;

     int mpi_size;

     int cur_rank;

     int n_active;

     int n_qsize2;
     int n_qsize2_w;
     int n_qsize2_r;

     int i_qsiz;

     int n_derivs;

     double *qx2;
     double **qx_l2;
     double *qw2;
     double **qw_l2;
     double *nr2;
     double **nr_l2;

     double cext2;
     double *cext_l2;
     double csca2;
     double *csca_l2;
     double cbak2;
     double *cbak_l2;
     double g2;
     double *g_l2;

     double **pf2;
     double ***pf_l2;

     MPI_Status status;


     n_derivs = n_derivs1 + n_derivs2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_qsize2_w = n_qsize / (mpi_size - 1);
     n_qsize2_r = n_qsize % (mpi_size - 1);

     if (verbose) {
          printf("%d %d %d\n", n_qsize, n_qsize2_w, n_qsize2_r);
          printf("\n");
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_mie_solution(n_qang, n_derivs, cext, csca, cbak,
                       g, pf, cext_l, csca_l, cbak_l, g_l, pf_l);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     qx2 = alloc_array1_d(n_qsize2_w + 1);
     qw2 = alloc_array1_d(n_qsize2_w + 1);
     nr2 = alloc_array1_d(n_qsize2_w + 1);

     pf2 = alloc_array2_d(4, n_qang);

     if (n_derivs > 0) {
          qx_l2   = alloc_array2_d(n_qsize2_w + 1, n_derivs);
          qw_l2   = alloc_array2_d(n_qsize2_w + 1, n_derivs);
          nr_l2   = alloc_array2_d(n_qsize2_w + 1, n_derivs);

          cext_l2 = alloc_array1_d(n_derivs);
          csca_l2 = alloc_array1_d(n_derivs);
          cbak_l2 = alloc_array1_d(n_derivs);
          g_l2    = alloc_array1_d(n_derivs);

          pf_l2   = alloc_array3_d(n_derivs, 4, n_qang);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_active = 0;

     i_qsiz   = 0;
     n_qsize2  = 0;
     for (i = 0; i < mpi_size - 1; ++i) {
          i_qsiz += n_qsize2;
          n_qsize2 = n_qsize2_w;
          if (i < n_qsize2_r)
               n_qsize2++;

          if (verbose)
               printf("send, %d: %d %d %d %d\n",
                      i, n_qsize, n_qsize2, i_qsiz, i_qsiz + n_qsize2 - 1);

          jj = i;
          for (j = 0; j < n_qsize2; ++j) {
               qx2[j] = qx[jj];
               qw2[j] = qw[jj];
               nr2[j] = nr[jj];

               for (k = 0; k < n_derivs; ++k) {
                    qx_l2[j][k] = qx_l[jj][k];
                    qw_l2[j][k] = qw_l[jj][k];
                    nr_l2[j][k] = nr_l[jj][k];
               }

               jj += mpi_size - 1;
          }

          MPI_Send(&n_qsize2,   1,       MPI_INT,    i + 1,  0, MPI_COMM_WORLD);

          if (n_qsize2 == 0)
               continue;

          n_active++;

          MPI_Send(&n1,        1,         MPI_INT,    i + 1,  1, MPI_COMM_WORLD);
          MPI_Send(&n2,        1,         MPI_INT,    i + 1,  2, MPI_COMM_WORLD);
          MPI_Send(&n_qang,    1,         MPI_INT,    i + 1,  3, MPI_COMM_WORLD);
          MPI_Send(&n_derivs1, 1,         MPI_INT,    i + 1,  4, MPI_COMM_WORLD);
          MPI_Send(&n_derivs2, 1,         MPI_INT,    i + 1,  5, MPI_COMM_WORLD);
          MPI_Send(index1,     n_derivs1, MPI_INT,    i + 1, 35, MPI_COMM_WORLD);
          MPI_Send(index2,     n_derivs2, MPI_INT,    i + 1, 36, MPI_COMM_WORLD);
          MPI_Send(&lambda,    1,         MPI_DOUBLE, i + 1,  6, MPI_COMM_WORLD);
          MPI_Send(&m,         2,         MPI_DOUBLE, i + 1,  7, MPI_COMM_WORLD);
          MPI_Send(&m2,        2,         MPI_DOUBLE, i + 1,  8, MPI_COMM_WORLD);
          MPI_Send(qx2,        n_qsize2,  MPI_DOUBLE, i + 1,  9, MPI_COMM_WORLD);
          MPI_Send(qw2,        n_qsize2,  MPI_DOUBLE, i + 1, 10, MPI_COMM_WORLD);
          MPI_Send(nr2,        n_qsize2,  MPI_DOUBLE, i + 1, 11, MPI_COMM_WORLD);
          MPI_Send(fv1,        n1 + 1,    MPI_DOUBLE, i + 1, 12, MPI_COMM_WORLD);
          MPI_Send(fv2,        n1 + 1,    MPI_DOUBLE, i + 1, 13, MPI_COMM_WORLD);
          MPI_Send(fv3,        n1 + 1,    MPI_DOUBLE, i + 1, 14, MPI_COMM_WORLD);
          MPI_Send(fv4,        n1 + 1,    MPI_DOUBLE, i + 1, 15, MPI_COMM_WORLD);
          MPI_Send(qx3,        n_qang,    MPI_DOUBLE, i + 1, 16, MPI_COMM_WORLD);
          MPI_Send(qw3,        n_qang,    MPI_DOUBLE, i + 1, 17, MPI_COMM_WORLD);
          MPI_Send(&n_threads, 1,         MPI_INT,    i + 1, 18, MPI_COMM_WORLD);

          if (n_derivs > 0) {
               MPI_Send(lambda_l, n_derivs,            MPI_DOUBLE, i + 1, 19, MPI_COMM_WORLD);
               MPI_Send(m_l,      n_derivs * 2,        MPI_DOUBLE, i + 1, 20, MPI_COMM_WORLD);
               MPI_Send(*qx_l2,   n_qsize2 * n_derivs, MPI_DOUBLE, i + 1, 21, MPI_COMM_WORLD);
               MPI_Send(*qw_l2,   n_qsize2 * n_derivs, MPI_DOUBLE, i + 1, 22, MPI_COMM_WORLD);
               MPI_Send(*nr_l2,   n_qsize2 * n_derivs, MPI_DOUBLE, i + 1, 23, MPI_COMM_WORLD);
          }
     }

     if (verbose)
          printf("\n");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_active; ++i) {
/*
          MPI_Recv(&cur_rank, 1,     MPI_INT,    MPI_ANY_SOURCE, 24, MPI_COMM_WORLD, &status);
*/
          MPI_Recv(&cur_rank, 1,     MPI_INT,    i + 1,          24, MPI_COMM_WORLD, &status);

          MPI_Recv(&cext2, 1,        MPI_DOUBLE, cur_rank,       25, MPI_COMM_WORLD, &status);
          MPI_Recv(&csca2, 1,        MPI_DOUBLE, cur_rank,       26, MPI_COMM_WORLD, &status);
          MPI_Recv(&cbak2, 1,        MPI_DOUBLE, cur_rank,       27, MPI_COMM_WORLD, &status);
          MPI_Recv(&g2,    1,        MPI_DOUBLE, cur_rank,       28, MPI_COMM_WORLD, &status);
          MPI_Recv(*pf2,   4*n_qang, MPI_DOUBLE, cur_rank,       29, MPI_COMM_WORLD, &status);

          if (n_derivs > 0) {
               MPI_Recv(cext_l2, n_derivs,          MPI_DOUBLE, cur_rank, 30, MPI_COMM_WORLD, &status);
               MPI_Recv(csca_l2, n_derivs,          MPI_DOUBLE, cur_rank, 31, MPI_COMM_WORLD, &status);
               MPI_Recv(cbak_l2, n_derivs,          MPI_DOUBLE, cur_rank, 32, MPI_COMM_WORLD, &status);
               MPI_Recv(g_l2,    n_derivs,          MPI_DOUBLE, cur_rank, 33, MPI_COMM_WORLD, &status);
               MPI_Recv(**pf_l2, n_derivs*4*n_qang, MPI_DOUBLE, cur_rank, 34, MPI_COMM_WORLD, &status);
          }

          add_mie_solution(n_qang, n_derivs, cext, csca, cbak, g, pf, cext_l,
                           csca_l, cbak_l, g_l, pf_l, cext2, csca2, cbak2, g2,
                           pf2, cext_l2, csca_l2, cbak_l2, g_l2, pf_l2);

          if (verbose)
               printf("recv, %d\n", i);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(qx2);
     free_array1_d(qw2);
     free_array1_d(nr2);

     free_array2_d(pf2);

     if (n_derivs > 0) {
          free_array2_d(qx_l2);
          free_array2_d(qw_l2);
          free_array2_d(nr_l2);

          free_array1_d(cext_l2);
          free_array1_d(csca_l2);
          free_array1_d(cbak_l2);
          free_array1_d(g_l2);

          free_array3_d(pf_l2);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int lmie_core_slave() {

     int mpi_rank;

     int n_qsize;

     int n1;
     int n2;

     int n_qang;

     int n_derivs;
     int n_derivs1;
     int n_derivs2;

     int n_threads;

     int *index1;
     int *index2;

     double lambda;
     double *lambda_l;

     double cext;
     double *cext_l;
     double csca;
     double *csca_l;
     double cbak;
     double *cbak_l;
     double g;
     double *g_l;

     double *qx;
     double **qx_l;
     double *qw;
     double **qw_l;
     double *nr;
     double **nr_l;

     double *fv1;
     double *fv2;
     double *fv3;
     double *fv4;

     double *qx2;
     double *qw2;

     double **pf;
     double ***pf_l;

     dcomplex m;
     dcomplex m2;
     dcomplex *m_l;


     MPI_Status status;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     MPI_Recv(&n_qsize,    1,     MPI_INT,     0,  0, MPI_COMM_WORLD, &status);

     if (n_qsize == 0)
          return 0;

     MPI_Recv(&n1,        1,         MPI_INT,    0,  1, MPI_COMM_WORLD, &status);
     MPI_Recv(&n2,        1,         MPI_INT,    0,  2, MPI_COMM_WORLD, &status);
     MPI_Recv(&n_qang,    1,         MPI_INT,    0,  3, MPI_COMM_WORLD, &status);
     MPI_Recv(&n_derivs1, 1,         MPI_INT,    0,  4, MPI_COMM_WORLD, &status);
     MPI_Recv(&n_derivs2, 1,         MPI_INT,    0,  5, MPI_COMM_WORLD, &status);

     n_derivs = n_derivs1 + n_derivs2;

     index1 = alloc_array1_i(n_derivs1);
     index2 = alloc_array1_i(n_derivs2);

     MPI_Recv(index1,     n_derivs1, MPI_INT,    0, 35, MPI_COMM_WORLD, &status);
     MPI_Recv(index2,     n_derivs2, MPI_INT,    0, 36, MPI_COMM_WORLD, &status);

     MPI_Recv(&lambda,    1,         MPI_DOUBLE, 0,  6, MPI_COMM_WORLD, &status);
     MPI_Recv(&m,         2,         MPI_DOUBLE, 0,  7, MPI_COMM_WORLD, &status);
     MPI_Recv(&m2,        2,         MPI_DOUBLE, 0,  8, MPI_COMM_WORLD, &status);

     qx  = alloc_array1_d(n_qsize);
     qw  = alloc_array1_d(n_qsize);
     nr  = alloc_array1_d(n_qsize);
     fv1 = alloc_array1_d(n1 + 1);
     fv2 = alloc_array1_d(n1 + 1);
     fv3 = alloc_array1_d(n1 + 1);
     fv4 = alloc_array1_d(n1 + 1);
     qx2 = alloc_array1_d(n_qang);
     qw2 = alloc_array1_d(n_qang);

     MPI_Recv(qx,         n_qsize,   MPI_DOUBLE, 0,  9, MPI_COMM_WORLD, &status);
     MPI_Recv(qw,         n_qsize,   MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &status);
     MPI_Recv(nr,         n_qsize,   MPI_DOUBLE, 0, 11, MPI_COMM_WORLD, &status);
     MPI_Recv(fv1,        n1 + 1,    MPI_DOUBLE, 0, 12, MPI_COMM_WORLD, &status);
     MPI_Recv(fv2,        n1 + 1,    MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &status);
     MPI_Recv(fv3,        n1 + 1,    MPI_DOUBLE, 0, 14, MPI_COMM_WORLD, &status);
     MPI_Recv(fv4,        n1 + 1,    MPI_DOUBLE, 0, 15, MPI_COMM_WORLD, &status);
     MPI_Recv(qx2,        n_qang,    MPI_DOUBLE, 0, 16, MPI_COMM_WORLD, &status);
     MPI_Recv(qw2,        n_qang,    MPI_DOUBLE, 0, 17, MPI_COMM_WORLD, &status);
     MPI_Recv(&n_threads, 1,         MPI_INT,    0, 18, MPI_COMM_WORLD, &status);

     if (n_derivs > 0) {
          lambda_l = alloc_array1_d (n_derivs);
          m_l      = alloc_array1_dc(n_derivs);

          MPI_Recv(lambda_l, n_derivs,        MPI_DOUBLE, 0, 19, MPI_COMM_WORLD, &status);
          MPI_Recv(m_l,      n_derivs * 2,    MPI_DOUBLE, 0, 20, MPI_COMM_WORLD, &status);

          qx_l = alloc_array2_d(n_qsize, n_derivs);
          qw_l = alloc_array2_d(n_qsize, n_derivs);
          nr_l = alloc_array2_d(n_qsize, n_derivs);

          MPI_Recv(*qx_l, n_qsize*n_derivs,    MPI_DOUBLE, 0, 21, MPI_COMM_WORLD, &status);
          MPI_Recv(*qw_l, n_qsize*n_derivs,    MPI_DOUBLE, 0, 22, MPI_COMM_WORLD, &status);
          MPI_Recv(*nr_l, n_qsize*n_derivs,    MPI_DOUBLE, 0, 23, MPI_COMM_WORLD, &status);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     pf = alloc_array2_d(4, n_qang);

     if (n_derivs > 0) {
          cext_l = alloc_array1_d (n_derivs);
          csca_l = alloc_array1_d (n_derivs);
          cbak_l = alloc_array1_d (n_derivs);
          g_l    = alloc_array1_d (n_derivs);

          pf_l   = alloc_array3_d(n_derivs, 4, n_qang);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (lmie_core_shared(n_qsize, n_derivs1, n_derivs2, index1, index2, lambda,
                          lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2, fv3, fv4,
                          qx2, qw2, &cext, &csca, &cbak, &g, pf, qx_l, qw_l,
                          nr_l, cext_l, csca_l, cbak_l, g_l, pf_l, n1, n2, n_qang,
                          n_threads)) {
          eprintf("ERROR: lmie_core_shared()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     MPI_Send(&mpi_rank, 1,    MPI_INT,    0, 24, MPI_COMM_WORLD);

     MPI_Send(&cext, 1,        MPI_DOUBLE, 0, 25, MPI_COMM_WORLD);
     MPI_Send(&csca, 1,        MPI_DOUBLE, 0, 26, MPI_COMM_WORLD);
     MPI_Send(&cbak, 1,        MPI_DOUBLE, 0, 27, MPI_COMM_WORLD);
     MPI_Send(&g,    1,        MPI_DOUBLE, 0, 28, MPI_COMM_WORLD);
     MPI_Send(*pf,   4*n_qang, MPI_DOUBLE, 0, 29, MPI_COMM_WORLD);

     if (n_derivs > 0) {
          MPI_Send(cext_l, n_derivs,          MPI_DOUBLE, 0, 30, MPI_COMM_WORLD);
          MPI_Send(csca_l, n_derivs,          MPI_DOUBLE, 0, 31, MPI_COMM_WORLD);
          MPI_Send(cbak_l, n_derivs,          MPI_DOUBLE, 0, 32, MPI_COMM_WORLD);
          MPI_Send(g_l,    n_derivs,          MPI_DOUBLE, 0, 33, MPI_COMM_WORLD);
          MPI_Send(**pf_l, n_derivs*4*n_qang, MPI_DOUBLE, 0, 34, MPI_COMM_WORLD);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_i(index1);
     free_array1_i(index2);

     free_array1_d(qx);
     free_array1_d(qw);
     free_array1_d(nr);

     free_array1_d(fv1);
     free_array1_d(fv2);
     free_array1_d(fv3);
     free_array1_d(fv4);

     free_array1_d(qx2);
     free_array1_d(qw2);

     if (n_derivs > 0) {
          free_array1_dc(m_l);

          free_array1_d(cext_l);
          free_array1_d(csca_l);
          free_array1_d(cbak_l);
          free_array1_d(g_l);

          free_array2_d(qx_l);
          free_array2_d(qw_l);
          free_array2_d(nr_l);
     }

     free_array2_d(pf);

     if (n_derivs > 0)
          free_array3_d(pf_l);


     return 0;
}

#endif


/*******************************************************************************
 *
 ******************************************************************************/
static void pf_4_to_6(double **pf4, double **pf6) {

     pf6[0] = pf4[0];
     pf6[1] = pf4[0];
     pf6[2] = pf4[1];
     pf6[3] = pf4[1];
     pf6[4] = pf4[2];
     pf6[5] = pf4[3];
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_core_solution(int n_qsize, int n_derivs1, int n_derivs2,
                       int *index1, int *index2,
                       double lambda, double *lambda_l,
                       double mr, double mi, double *mr_l, double *mi_l,
                       double *qx, double *qw, double *nr,
                       double r1, double r2,
                       double *cext, double *csca, double *cbak,
                       double *g, double **gc, double **lc,
                       double **qx_l, double **qw_l, double **nr_l,
                       double *cext_l, double *csca_l, double *cbak_l,
                       double *g_l, double ***gc_l, double ***lc_l,
                       double accuracy, int verbose, int n_threads, int use_mpi) {

     int i;
     int j;

     int n1;
     int n2;

     int n_qang;

     int n_coef;
     int n_coef2;

     int n_derivs;
#ifdef USE_MPI
     int mpi_size;
#endif
     double a;
     double a_l;

     double x;

     double *qx2;
     double *qw2;

     double *fv1;
     double *fv2;
     double *fv3;
     double *fv4;

     double **pf;
     double ***pf_l;

     double *pf2[6];

     double **Y;

     double ***P;

     dcomplex z;
     dcomplex m;
     dcomplex m2;
     dcomplex *m_l;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     x = 2. * PI * r2 / lambda;
     z = (mr - _Complex_I * mi) * x;

     n1 = calc_n1(x);
     n2 = calc_n2(n1, z);
/*
     n_qang = 2 * n1 - 1;
*/
     n_qang = 2 * n1;
     n_coef = 2 * n1;

     n_derivs = n_derivs1 + n_derivs2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     fv1 = alloc_array1_d(n1 + 1);
     fv2 = alloc_array1_d(n1 + 1);
     fv3 = alloc_array1_d(n1 + 1);
     fv4 = alloc_array1_d(n1 + 1);

     qx2 = alloc_array1_d(n_qang);
     qw2 = alloc_array1_d(n_qang);

     pf  = alloc_array2_d(4, n_qang);

     if (n_derivs > 0) {
          m_l = alloc_array1_dc(n_derivs);

          pf_l = alloc_array3_d(n_derivs, 4, n_qang);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     m = mr - _Complex_I * mi;

     if (n_derivs > 0)
          m2  = m * m;

     for (i = 0; i < n_derivs; ++i)
          m_l[i] = mr_l[i] - _Complex_I * mi_l[i];


    /*--------------------------------------------------------------------
     *
     *------------------------------------------------------------------*/
     for (i = 1; i <= n1; ++i) {
          fv1[i] = 2. * i + 1.;
          fv2[i] = (i + 1.) / i;
          fv3[i] = (i + 2.) / fv2[i];
          fv4[i] = fv1[i] / (i * (i + 1.));
     }


    /*--------------------------------------------------------------------
     *
     *------------------------------------------------------------------*/
     gauss_leg_quadx(n_qang, -1., 1., qx2, qw2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifdef USE_MPI
     if (use_mpi)
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

     if (! use_mpi || mpi_size == 1) {
#endif
          if (lmie_core_shared
                        (n_qsize, n_derivs1, n_derivs2, index1, index2, lambda,
                         lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2, fv3, fv4,
                         qx2, qw2, cext, csca, cbak, g, pf, qx_l, qw_l, nr_l,
                         cext_l, csca_l, cbak_l, g_l, pf_l, n1, n2, n_qang,
                         n_threads)) {
               eprintf("ERROR: lmie_core_shared()\n");
               return -1;
          }
#ifdef USE_MPI
     }
     else {
          if (lmie_core_master
                        (n_qsize, n_derivs1, n_derivs2, index1, index2, lambda,
                         lambda_l, m, m_l, m2, qx, qw, nr, fv1, fv2, fv3, fv4,
                         qx2, qw2, cext, csca, cbak, g, pf, qx_l, qw_l, nr_l,
                         cext_l, csca_l, cbak_l, g_l, pf_l, n1, n2, n_qang,
                         verbose, n_threads)) {
               eprintf("ERROR: lmie_core_master()\n");
               return -1;
          }
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = 2. / *csca;

     for (i = 0; i < n_derivs; ++i) {
          a_l = -a * csca_l[i] / *csca;

          g_l[i] = g_l[i] * a + *g * a_l;

          for (j = 0; j < n_qang; ++j) {
               pf_l[i][0][j] = pf_l[i][0][j] * a + pf[0][j] * a_l;
               pf_l[i][1][j] = pf_l[i][1][j] * a + pf[1][j] * a_l;
               pf_l[i][2][j] = pf_l[i][2][j] * a + pf[2][j] * a_l;
               pf_l[i][3][j] = pf_l[i][3][j] * a + pf[3][j] * a_l;
          }
     }

     *g = *g / *csca * 2.;

     for (i = 0; i < n_qang; ++i) {
          pf[0][i] *= a;
          pf[1][i] *= a;
          pf[2][i] *= a;
          pf[3][i] *= a;
     }


     /*-------------------------------------------------------------------------
      * multiply by 2. / x^2 to get efficiency, then by PI * r^2 to get cross
      * section ... 2. / x * PI * r^2 = lambda^2 / (2. * PI)
      *-----------------------------------------------------------------------*/
     a = lambda * lambda / (2. * PI);

     for (i = 0; i < n_derivs; ++i) {
          a_l = 2. * lambda_l[i] * lambda / (2. * PI);

          cext_l[i] =  cext_l[i] * a + *cext * a_l;
          csca_l[i] =  csca_l[i] * a + *csca * a_l;
          cbak_l[i] = (cbak_l[i] * a + *cbak * a_l) / (2. * 4. * PI);
     }

     *cext = *cext * a;
     *csca = *csca * a;
     *cbak = *cbak * a / (2. * 4. * PI);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     int init_gsf  = 1;
     int symmetric = 1;

     if (gc) {
          P = NULL;
          if (init_gsf)
               phase_func_integ_gc_init(n_coef, n_qang, qx2, qw2, &P, symmetric);

          pf_4_to_6(pf, pf2);
          n_coef2 = phase_func_integ_gc(n_coef, n_qang, pf2, qx2, qw2, P, gc, accuracy, 1, symmetric);
          for (i = 0; i < n_derivs; ++i) {
               pf_4_to_6(pf_l[i], pf2);
               phase_func_integ_gc(n_coef, n_qang, pf2, qx2, qw2, P, gc_l[i], accuracy, 1, symmetric);

               gc_l[i][0][0] = 0.;
          }

          if (init_gsf)
               free_array3_d(P);
     }

     if (lc) {
          Y = alloc_array2_d((n_qang + 1) / 2, n_coef);

          phase_func_integ_lc_init(n_coef, n_qang, qx2, qw2, Y);

          pf_4_to_6(pf, pf2);
          n_coef2 = phase_func_integ_lc(n_coef, n_qang, pf2, qx2, qw2, Y, lc, accuracy, 1);
          for (i = 0; i < n_derivs; ++i) {
               pf_4_to_6(pf_l[i], pf2);
               phase_func_integ_lc(n_coef, n_qang, pf2, qx2, qw2, Y, lc_l[i], accuracy, 1);

               lc_l[i][0][0] = 0.;
          }

          free_array2_d(Y);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(fv1);
     free_array1_d(fv2);
     free_array1_d(fv3);
     free_array1_d(fv4);

     free_array1_d(qx2);
     free_array1_d(qw2);

     free_array2_d(pf);

     if (n_derivs > 0) {
          free_array1_dc(m_l);

          free_array3_d(pf_l);
     }


     return n_coef2;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef USE_MPI

int lmie_core_solution_slave() {

     if (lmie_core_slave()) {
          eprintf("ERROR: lmie_core_slave()\n");
          return -1;
     }

     return 0;
}

#endif
