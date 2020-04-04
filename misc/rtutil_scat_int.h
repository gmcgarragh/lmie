/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef RTUTIL_SCAT_INT_H
#define RTUTIL_SCAT_INT_H

#ifdef __cplusplus
extern "C" {
#endif


void phase_func_integ_gc_init(int n_coef, int n_quad,
                              double *qx, double *qw, double ****P, int symmetric);
int phase_func_integ_gc(int n_coef, int n_quad, double **pf,
                        double *qx, double *qw, double ***P,
                        double **coefs, double accuracy, int allsix, int symmetric);
int phase_func_integ_gc2(int n_coef, int n_quad, double **pf,
                         double *qx, double *qw, double ***P,
                         double **coefs, double accuracy, int allsix, int symmetric);
void phase_func_integ_lc_init(int n_coef, int n_quad,
                              double *qx, double *qw, double **P_p);
int phase_func_integ_lc(int n_coef, int n_quad, double **pf,
                        double *qx, double *qw, double **P_p,
                        double **coefs, double accuracy, int allsix);

double phase_func_d_theta(int n_theta);

void phase_func_theta(double *theta, int n_theta, double d_theta);

void phase_func_exp_gc(double **coefs, int n_coef,
                       double *mu, double **pf, int n_pf, int allsix);
void phase_func_exp_gc_sym(double **coefs, int n_coef,
                           double *mu, double **pf, int n_pf, int allsix);
void phase_func_exp_lc(double **coefs, int n_coef,
                       double *mu, double **pf, int n_pf, int allsix);
void phase_func_exp_lc_sym(double **coefs, int n_coef,
                           double *mu, double **pf, int n_pf, int allsix);

int create_phase_func(int n_coef, int n_theta, int n_derivs, double **xc,
                      double *theta, double **pf, double ***xc_l, double ***pf_l,
                      int flag);

double asymmetry_parameter(int n_theta, double *theta, double **phase);


#ifdef __cplusplus
}
#endif

#endif /* RTUTIL_SCAT_INT_H */
