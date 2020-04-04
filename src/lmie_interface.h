/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef LMIE_INTERFACE_H
#define LMIE_INTERFACE_H

#include <rtutil_support.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     int save_control;

     int calc_gc;
     int calc_lc;
     int calc_pf;

     enum size_dist_type dist_type;

     int n_int1;
     int n_int2;
     int n_quad;
     int n_angles;
     int n_derivs;

     double lambda;
     double mr;
     double mi;
     double a1;
     double a2;
     double a3;
     double a4;
     double a5;
     double r1;
     double r2;

     double *lambda_l;
     double *mr_l;
     double *mi_l;
     double *a1_l;
     double *a2_l;
     double *a3_l;
     double *a4_l;
     double *a5_l;
     double *r1_l;
     double *r2_l;

     double accuracy;
} lmie_in_data;


typedef struct {
     int n_coef;

     double r1;
     double r2;
     double norm;
     double reff;
     double veff;
     double gavg;
     double vavg;
     double ravg;
     double rvw;
     double cext;
     double csca;
     double cbak;
     double g;
     double **gc;
     double **lc;
     double *theta;
     double **pf;

     double *r1_l;
     double *r2_l;
     double *norm_l;
     double *reff_l;
     double *veff_l;
     double *gavg_l;
     double *vavg_l;
     double *ravg_l;
     double *rvw_l;
     double *cext_l;
     double *csca_l;
     double *cbak_l;
     double *g_l;
     double ***gc_l;
     double ***lc_l;
     double ***pf_l;

     double **save1;
     double ***save2;
} lmie_out_data;


int lmie_calc_max_coef(double lambda, enum size_dist_type dist_type,
                       double a1, double a2, double a3, double r1, double r2);

int lmie_in_alloc(lmie_in_data *in, int n_derivs);
void lmie_in_free(lmie_in_data *in, int flag);
int lmie_out_alloc(lmie_in_data *in, lmie_out_data *out, int max_coef);
void lmie_out_free(lmie_out_data *out, int flag);

void lmie_in_zero_derivs(lmie_in_data *in, int n_derivs);

int lmie_solution(lmie_in_data *in, lmie_out_data *out, int alloc_out, int verbose,
                  int n_threads, int use_mpi);

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
                   int verbose, int n_threads, int use_mpi);

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
                     int verbose, int n_threads, int use_mpi);

int lmie_solution_slave(void);



#ifdef __cplusplus
}
#endif

#endif /* LMIE_INTERFACE_H */
