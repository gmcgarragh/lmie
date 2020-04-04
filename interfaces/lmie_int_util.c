/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <rtutil_support.h>

#include <lmie_interface.h>


#ifdef __cplusplus
extern "C" {
#endif

int size_dist_code_(char *name);
int lmie_calc_max_coef_(double *lambda, enum size_dist_type *dist_type,
                        double *a1, double *a2, double *a3, double *r1, double *r2);
int lmie_solution2_(int *calc_gc, int *calc_lc, int *calc_pf,
                    enum size_dist_type *dist_type,
                    int *n_int1, int *n_int2,
                    int *n_quad, int *n_angles, int *save_control,
                    double *lambda, double *mr, double *mi,
                    double *a1, double *a2, double *a3, double *a4, double *a5,
                    double *r1, double *r2,
                    double *accuracy, int *n_coef,
                    double *r21, double *r22,
                    double *norm, double *reff, double *veff,
                    double *gavg, double *vavg, double *ravg, double *rvw,
                    double *cext, double *csca, double *cbak, double *g,
                    double *gc, double *lc, double *theta, double *pf,
                    double ***save1, double ****save2,
                    int *max_coef, int *verbose, int *n_threads, int *use_mpi);
int lmie_solution2_l_(int *calc_gc, int *calc_lc, int *calc_pf,
                      enum size_dist_type *dist_type,
                      int *n_int1, int *n_int2, int *n_quad,
                      int *n_angles, int *n_derivs, int *save_control,
                      double *lambda, double *mr, double *mi,
                      double *a1, double *a2, double *a3, double *a4, double *a5,
                      double *r1, double *r2,
                      double *lambda_l, double *mr_l, double *mi_l,
                      double *a1_l, double *a2_l, double *a3_l, double *a4_l, double *a5_l,
                      double *r1_l, double *r2_l,
                      double *accuracy, int *n_coef,
                      double *r21, double *r22,
                      double *norm, double *reff, double *veff,
                      double *gavg, double *vavg, double *ravg, double *rvw,
                      double *cext, double *csca, double *cbak, double *g,
                      double *gc, double *lc, double *theta, double *pf,
                      double *r21_l, double *r22_l,
                      double *norm_l, double *reff_l, double *veff_l,
                      double *gavg_l, double *vavg_l, double *ravg_l, double *rvw_l,
                      double *cext_l, double *csca_l, double *cbak_l, double *g_l,
                      double *gc_l, double *lc_l, double *pf_l,
                      double ***save1, double ****save2,
                      int *max_coef, int *verbose, int *n_threads, int *use_mpi);

#ifdef __cplusplus
}
#endif



/*******************************************************************************
 *
 ******************************************************************************/
int size_dist_code_(char *name) {

     return size_dist_code(name);
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_calc_max_coef_(double *lambda, enum size_dist_type *dist_type,
                        double *a1, double *a2, double *a3, double *r1, double *r2) {

     return lmie_calc_max_coef(*lambda, *dist_type, *a1, *a2, *a3, *r1, *r2);
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_solution2_(int *calc_gc, int *calc_lc, int *calc_pf,
                    enum size_dist_type *dist_type,
                    int *n_int1, int *n_int2,
                    int *n_quad, int *n_angles, int *save_control,
                    double *lambda, double *mr, double *mi,
                    double *a1, double *a2, double *a3, double *a4, double *a5,
                    double *r1, double *r2,
                    double *accuracy, int *n_coef,
                    double *r21, double *r22,
                    double *norm, double *reff, double *veff,
                    double *gavg, double *vavg, double *ravg, double *rvw,
                    double *cext, double *csca, double *cbak, double *g,
                    double *gc, double *lc, double *theta, double *pf,
                    double ***save1, double ****save2,
                    int *max_coef, int *verbose, int *n_threads, int *use_mpi) {

     return lmie_solution2_l_(calc_gc, calc_lc, calc_pf,
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
                              max_coef, verbose, n_threads, use_mpi);
}



/*******************************************************************************
 *
 ******************************************************************************/
int lmie_solution2_l_(int *calc_gc, int *calc_lc, int *calc_pf,
                      enum size_dist_type *dist_type,
                      int *n_int1, int *n_int2, int *n_quad,
                      int *n_angles, int *n_derivs, int *save_control,
                      double *lambda, double *mr, double *mi,
                      double *a1, double *a2, double *a3, double *a4, double *a5,
                      double *r1, double *r2,
                      double *lambda_l, double *mr_l, double *mi_l,
                      double *a1_l, double *a2_l, double *a3_l, double *a4_l, double *a5_l,
                      double *r1_l, double *r2_l,
                      double *accuracy, int *n_coef,
                      double *r21, double *r22,
                      double *norm, double *reff, double *veff,
                      double *gavg, double *vavg, double *ravg, double *rvw,
                      double *cext, double *csca, double *cbak, double *g,
                      double *gc, double *lc, double *theta, double *pf,
                      double *r21_l, double *r22_l,
                      double *norm_l, double *reff_l, double *veff_l,
                      double *gavg_l, double *vavg_l, double *ravg_l, double *rvw_l,
                      double *cext_l, double *csca_l, double *cbak_l, double *g_l,
                      double *gc_l, double *lc_l, double *pf_l,
                      double ***save1, double ****save2,
                      int *max_coef, int *verbose, int *n_threads, int *use_mpi) {

     double **gc2;
     double **lc2;
     double **pf2;

     double ***gc_l2;
     double ***lc_l2;
     double ***pf_l2;

     gc2   = NULL;
     gc_l2 = NULL;
     if (*calc_gc) {
          gc2 = array_from_mem2_d(gc, 6, *max_coef);
          if (*n_derivs > 0)
               gc_l2 = array_from_mem3_d(gc_l, *n_derivs, 6, *max_coef);
     }

     lc2   = NULL;
     lc_l2 = NULL;
     if (*calc_lc) {
          lc2 = array_from_mem2_d(lc, 6, *max_coef);
          if (*n_derivs > 0)
               lc_l2 = array_from_mem3_d(lc_l, *n_derivs, 6, *max_coef);
     }

     pf2   = NULL;
     pf_l2 = NULL;
     if (*calc_pf) {
          pf2   = array_from_mem2_d(pf, 6, *n_angles);
          if (*n_derivs > 0)
               pf_l2 = array_from_mem3_d(pf_l, *n_derivs, 6, *n_angles);
     }

     if (lmie_solution2_l(*calc_gc, *calc_lc, *calc_pf,
                          *dist_type,
                          *n_int1, *n_int2,
                          *n_quad, *n_angles, *n_derivs,
                          *save_control,
                          *lambda, *mr, *mi,
                          *a1, *a2, *a3, *a4, *a5,
                          *r1, *r2,
                          lambda_l, mr_l, mi_l,
                          a1_l, a2_l, a3_l, a4_l, a5_l,
                          r1_l, r2_l,
                          *accuracy, n_coef,
                          r21, r22,
                          norm, reff, veff,
                          gavg, vavg, ravg, rvw,
                          cext, csca, cbak, g,
                          gc2, lc2, theta, pf2,
                          r21_l, r22_l,
                          norm_l, reff_l, veff_l,
                          gavg_l, vavg_l, ravg_l, rvw_l,
                          cext_l, csca_l, cbak_l, g_l,
                          gc_l2, lc_l2, pf_l2,
                          save1, save2,
                          *verbose, *n_threads, *use_mpi)) {
          fprintf(stderr, "ERROR: lmie_solution2_l()\n");
          return -1;
     }

     if (gc) {
          free_array2_d(gc2);
          if (*n_derivs > 0)
               free_array3_d(gc_l2);
     }
     if (lc) {
          free_array2_d(lc2);
          if (*n_derivs > 0)
               free_array3_d(lc_l2);
     }
     if (pf) {
          free_array2_d(pf2);
          if (*n_derivs > 0)
               free_array3_d(pf_l2);
     }

     return 0;
}
