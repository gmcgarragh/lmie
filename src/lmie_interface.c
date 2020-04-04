/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <rtutil_scat_int.h>
#include <rtutil_size_dist.h>
#include <rtutil_support.h>

#include "lmie_core.h"
#include "lmie_interface.h"
#include "lmie_schedule.h"


/*******************************************************************************
 *
 ******************************************************************************/
int lmie_calc_max_coef(double lambda, enum size_dist_type dist_type,
                       double a1, double a2, double a3, double r1, double r2) {

     if (dist_type == SIZE_DIST_POWER_LAW)
          get_power_law_range(a1, a2, &r1, &r2);

     return 2 * calc_n1(2. * PI * r2 / lambda);
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_in_alloc(lmie_in_data *in, int n_derivs) {

     if (n_derivs > 0) {
          in->lambda_l = alloc_array1_d(n_derivs);
          in->mr_l     = alloc_array1_d(n_derivs);
          in->mi_l     = alloc_array1_d(n_derivs);
          in->a1_l     = alloc_array1_d(n_derivs);
          in->a2_l     = alloc_array1_d(n_derivs);
          in->a3_l     = alloc_array1_d(n_derivs);
          in->a4_l     = alloc_array1_d(n_derivs);
          in->a5_l     = alloc_array1_d(n_derivs);
          in->r1_l     = alloc_array1_d(n_derivs);
          in->r2_l     = alloc_array1_d(n_derivs);
     }

     return 0;
}



void lmie_in_free(lmie_in_data *in, int flag) {

     if (flag) {
          free_array1_d(in->lambda_l);
          free_array1_d(in->mr_l);
          free_array1_d(in->mi_l);
          free_array1_d(in->a1_l);
          free_array1_d(in->a2_l);
          free_array1_d(in->a3_l);
          free_array1_d(in->a4_l);
          free_array1_d(in->a5_l);
          free_array1_d(in->r1_l);
          free_array1_d(in->r2_l);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_out_alloc(lmie_in_data *in, lmie_out_data *out, int max_coef) {

     if (in->n_derivs > 0) {
          out->r1_l   = alloc_array1_d(in->n_derivs);
          out->r2_l   = alloc_array1_d(in->n_derivs);
          out->norm_l = alloc_array1_d(in->n_derivs);
          out->reff_l = alloc_array1_d(in->n_derivs);
          out->veff_l = alloc_array1_d(in->n_derivs);
          out->gavg_l = alloc_array1_d(in->n_derivs);
          out->vavg_l = alloc_array1_d(in->n_derivs);
          out->ravg_l = alloc_array1_d(in->n_derivs);
          out->rvw_l  = alloc_array1_d(in->n_derivs);
          out->cext_l = alloc_array1_d(in->n_derivs);
          out->csca_l = alloc_array1_d(in->n_derivs);
          out->cbak_l = alloc_array1_d(in->n_derivs);
          out->g_l    = alloc_array1_d(in->n_derivs);
     }

     out->gc = NULL;
     if (in->calc_gc) {
          out->gc = alloc_array2_d(6, max_coef);
          if (in->n_derivs > 0)
               out->gc_l = alloc_array3_d(in->n_derivs, 6, max_coef);
     }

     out->lc = NULL;
     if (in->calc_lc) {
          out->lc = alloc_array2_d(6, max_coef);
          if (in->n_derivs > 0)
               out->lc_l = alloc_array3_d(in->n_derivs, 6, max_coef);
     }

     out->pf = NULL;
     if (in->calc_pf) {
          out->theta = alloc_array1_d(   in->n_angles);
          out->pf    = alloc_array2_d(6, in->n_angles);
          if (in->n_derivs > 0)
               out->pf_l = alloc_array3_d(in->n_derivs, 6, in->n_angles);
     }

     if (in->save_control != 2) {
          out->save1 = NULL;
          out->save2 = NULL;
     }

     return 0;
}



void lmie_out_free(lmie_out_data *out, int flag) {

     if (flag) {
          free_array1_d(out->r1_l);
          free_array1_d(out->r2_l);
          free_array1_d(out->norm_l);
          free_array1_d(out->reff_l);
          free_array1_d(out->veff_l);
          free_array1_d(out->gavg_l);
          free_array1_d(out->vavg_l);
          free_array1_d(out->ravg_l);
          free_array1_d(out->rvw_l);
          free_array1_d(out->cext_l);
          free_array1_d(out->csca_l);
          free_array1_d(out->cbak_l);
          free_array1_d(out->g_l);
     }

     if (out->gc) {
          free_array2_d(out->gc);
          if (flag)
               free_array3_d(out->gc_l);
     }
     if (out->lc) {
          free_array2_d(out->lc);
          if (flag)
               free_array3_d(out->lc_l);
     }
     if (out->pf) {
          free_array1_d(out->theta);
          free_array2_d(out->pf);
          if (flag)
               free_array3_d(out->pf_l);
     }
/*
     if (out->save1) {
          free_array2_d(out->save1);
          free_array3_d(out->save2);
     }
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void lmie_in_zero_derivs(lmie_in_data *in, int n_derivs) {

     int i;

     for (i = 0; i < n_derivs; ++i) {
          in->lambda_l[i] = 0.;
          in->mr_l[i]     = 0.;
          in->mi_l[i]     = 0.;
          in->a1_l[i]     = 0.;
          in->a2_l[i]     = 0.;
          in->a3_l[i]     = 0.;
          in->a4_l[i]     = 0.;
          in->a5_l[i]     = 0.;
          in->r1_l[i]     = 0.;
          in->r2_l[i]     = 0.;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_solution(lmie_in_data *in, lmie_out_data *out, int alloc_out, int verbose,
                  int n_threads, int use_mpi) {

     int max_coef;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (alloc_out) {
          max_coef = lmie_calc_max_coef(in->lambda, in->dist_type, in->a1, in->a2,
                                        in->a3, in->r1, in->r2);

          if (lmie_out_alloc(in, out, max_coef)) {
               fprintf(stderr, "ERROR: lmie_out_alloc()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (lmie_solution2_l(in->calc_gc, in->calc_lc, in->calc_pf,
                          in->dist_type, in->n_int1, in->n_int2, in->n_quad,
                          in->n_angles, in->n_derivs, in->save_control,
                          in->lambda, in->mr, in->mi,
                          in->a1, in->a2, in->a3, in->a4, in->a5,
                          in->r1, in->r2,
                          in->lambda_l, in->mr_l, in->mi_l,
                          in->a1_l, in->a2_l, in->a3_l, in->a4_l, in->a5_l,
                          in->r1_l, in->r2_l,
                          in->accuracy, &out->n_coef,
                          &out->r1, &out->r2,
                          &out->norm, &out->reff, &out->veff,
                          &out->gavg, &out->vavg, &out->ravg, &out->rvw,
                          &out->cext, &out->csca, &out->cbak, &out->g,
                          out->gc, out->lc, out->theta, out->pf,
                          out->r1_l, out->r2_l,
                          out->norm_l, out->reff_l, out->veff_l,
                          out->gavg_l, out->vavg_l, out->ravg_l, out->rvw_l,
                          out->cext_l, out->csca_l, out->cbak_l, out->g_l,
                          out->gc_l, out->lc_l, out->pf_l,
                          &out->save1, &out->save2,
                          verbose, n_threads, use_mpi)) {
          fprintf(stderr, "ERROR: lmie_solution2_l()\n");
          return -1;
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_solution2(int calc_gc, int calc_lc, int calc_pf,
                   enum size_dist_type dist_type,
                   int n_int1, int n_int2,
                   int n_quad, int n_angles, int save_control,
                   double lambda, double mr, double mi,
                   double a1, double a2, double a3, double a4, double a5,
                   double r1, double r2,
                   double accuracy, int *n_coef,
                   double *r21, double *r22,
                   double *norm, double *reff, double *veff,
                   double *gavg, double *vavg, double *ravg, double *rvw,
                   double *cext, double *csca, double *cbak, double *g,
                   double **gc, double **lc, double *theta, double **pf,
                   double ***save1, double ****save2,
                   int verbose, int n_threads, int use_mpi) {

     return lmie_solution2_l(calc_gc, calc_lc, calc_pf,
                             dist_type,
                             n_int1, n_int2,
                             n_quad, n_angles, 0,
                             save_control,
                             lambda, mr, mi,
                             a1, a2, a3, a4, a5,
                             r1, r2,
                             NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL,
                             accuracy, n_coef,
                             r21, r22,
                             norm, reff, veff,
                             gavg, vavg, ravg, rvw,
                             cext, csca, cbak, g,
                             gc, lc, theta, pf,
                             NULL, NULL,
                             NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL,
                             save1, save2,
                             verbose, n_threads, use_mpi);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int create_derivs_indexes(int n_qsize, int n_derivs,
                                 double *lambda_l,
                                 double *mr_l, double *mi_l,
                                 double *r1_l, double *r2_l,
                                 int *n_derivs1, int *n_derivs2,
                                 int **index1, int **index2) {

     int i;
     int i1;
     int i2;

     if (n_derivs <= 0) {
          *n_derivs1 = 0;
          *n_derivs2 = 0;
          return 0;
     }

     *n_derivs1 = 0;
     *n_derivs2 = 0;
     for (i = 0; i < n_derivs; ++i) {
/*
          if (1)
*/
          if (lambda_l[i] != 0. || mr_l[i] != 0. || mi_l[i] != 0. ||
              r1_l[i] != 0. || r2_l[i] != 0.)
               (*n_derivs1)++;
          else
               (*n_derivs2)++;
     }

     if (*n_derivs1 > 0)
          *index1 = alloc_array1_i(*n_derivs1);
     if (*n_derivs2 > 0)
          *index2 = alloc_array1_i(*n_derivs2);

     i1 = 0;
     i2 = 0;
     for (i = 0; i < n_derivs; ++i) {
/*
          if (1)
*/
          if (lambda_l[i] != 0. || mr_l[i] != 0. || mi_l[i] != 0. ||
              r1_l[i] != 0. || r2_l[i] != 0.)
               (*index1)[i1++] = i;
          else
               (*index2)[i2++] = i;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_solution2_l(int calc_gc, int calc_lc, int calc_pf,
                     enum size_dist_type dist_type,
                     int n_int1, int n_int2, int n_quad,
                     int n_angles, int n_derivs, int save_control,
                     double lambda, double mr, double mi,
                     double a1, double a2, double a3, double a4, double a5,
                     double r1, double r2,
                     double *lambda_l, double *mr_l, double *mi_l,
                     double *a1_l, double *a2_l, double *a3_l, double *a4_l, double *a5_l,
                     double *r1_l, double *r2_l,
                     double accuracy, int *n_coef,
                     double *r21, double *r22,
                     double *norm, double *reff, double *veff,
                     double *gavg, double *vavg, double *ravg, double *rvw,
                     double *cext, double *csca, double *cbak, double *g,
                     double **gc, double **lc, double *theta, double **pf,
                     double *r21_l, double *r22_l,
                     double *norm_l, double *reff_l, double *veff_l,
                     double *gavg_l, double *vavg_l, double *ravg_l, double *rvw_l,
                     double *cext_l, double *csca_l, double *cbak_l, double *g_l,
                     double ***gc_l, double ***lc_l, double ***pf_l,
                     double ***save1, double ****save2,
                     int verbose, int n_threads, int use_mpi) {

     int n_qsize;

     int n_derivs1;
     int n_derivs2;

     int *index1;
     int *index2;

     double *qx;
     double *qw;
     double *nr;

     double **qx_l;
     double **qw_l;
     double **nr_l;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (calc_pf && (! calc_gc && ! calc_lc)) {
          fprintf(stderr, "ERROR: either calc_gc or calc_lc must be set to true when "
                  "calc_pf is true\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_qsize = get_particle_dist_n(dist_type, n_int1, n_int2, n_quad);


     qx = alloc_array1_d(n_qsize);
     qw = alloc_array1_d(n_qsize);
     nr = alloc_array1_d(n_qsize);

     if (n_derivs > 0) {
          qx_l = alloc_array2_d(n_qsize, n_derivs);
          qw_l = alloc_array2_d(n_qsize, n_derivs);
          nr_l = alloc_array2_d(n_qsize, n_derivs);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (get_particle_dist(dist_type, n_int1, n_int2, n_quad, n_derivs, a1, a2,
                           a3, a4, a5, r1, r2, a1_l, a2_l, a3_l, a4_l, a5_l,
                           r1_l, r2_l, qx, qw, nr, r21, r22, qx_l, qw_l, nr_l,
                           r21_l, r22_l, norm, norm_l)) {
          fprintf(stderr, "ERROR: get_particle_dist()\n");
          return -1;
     }

     if (get_dist_parameters(n_qsize,  n_derivs,  qx, qw, nr,  qx_l, qw_l, nr_l,
                             reff, veff, gavg, vavg, ravg, rvw, reff_l, veff_l,
                             gavg_l, vavg_l, ravg_l, rvw_l)) {
          fprintf(stderr, "ERROR: get_dist_parameters()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     create_derivs_indexes(n_qsize, n_derivs, lambda_l, mr_l, mi_l, r21_l, r22_l,
                           &n_derivs1, &n_derivs2, &index1, &index2);

     if ((*n_coef = lmie_core_solution(n_qsize, n_derivs1, n_derivs2, save_control,
                                       index1, index2, lambda, lambda_l, mr, mi,
                                       mr_l, mi_l, qx, qw, nr, *r21, *r22, cext,
                                       csca, cbak, g, gc, lc, qx_l, qw_l, nr_l,
                                       cext_l, csca_l, cbak_l, g_l, gc_l, lc_l,
                                       save1, save2, accuracy, verbose, n_threads,
                                       use_mpi)) < 0) {
          fprintf(stderr, "ERROR: lmie_core_solution()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (calc_pf) {
          if (calc_gc) {
               if (create_phase_func(*n_coef, n_angles, n_derivs, gc, theta, pf,
                                     gc_l, pf_l, 1)) {
                    fprintf(stderr, "ERROR: create_phase_func()\n");
                    return -1;
               }
          }
          else {
               if (create_phase_func(*n_coef, n_angles, n_derivs, lc, theta, pf,
                                     lc_l, pf_l, 0)) {
                    fprintf(stderr, "ERROR: create_phase_func()\n");
                    return -1;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(qx);
     free_array1_d(qw);
     free_array1_d(nr);

     if (n_derivs > 0) {
          free_array2_d(qx_l);
          free_array2_d(qw_l);
          free_array2_d(nr_l);
     }

     if (n_derivs1 > 0)
          free_array1_i(index1);

     if (n_derivs2 > 0)
          free_array1_i(index2);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef USE_MPI

int lmie_solution_slave(void) {

     if (lmie_core_solution_slave()) {
          fprintf(stderr, "ERROR: lmie_core_solution_slave()\n");
          return -1;
     }

     return 0;
}

#endif
