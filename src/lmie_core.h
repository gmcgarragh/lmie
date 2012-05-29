/******************************************************************************%
**
**    Copyright (C) 2008-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef LMIE_CORE_H
#define LMIE_CORE_H

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     int i;
     int flag;
     int n_qsize;
     int n_derivs1;
     int n_derivs2;
     int *index1;
     int *index2;
     double lambda;
     double *lambda_l;
     dcomplex m;
     dcomplex *m_l;
     dcomplex m2;
     double *qx;
     double *qw;
     double *nr;
     double *fv1;
     double *fv2;
     double *fv3;
     double *fv4;
     double *qx2;
     double *qw2;
     double **qx_l;
     double **qw_l;
     double **nr_l;
     int n1;
     int n2;
     int n_qang;
} lmie_core_shared_data;


typedef struct {
     int n_qsize;
     int *index;
     double *cext;
     double *csca;
     double *cbak;
     double *g;
     double **pf;
     double *cext_l;
     double *csca_l;
     double *cbak_l;
     double *g_l;
     double ***pf_l;
     lmie_core_shared_data *s;
} lmie_core_threads_data;


int calc_n1(double x);
int calc_n2(int n1, dcomplex z);

void lmie_core_shared_import(lmie_core_shared_data *d,
                             int n_qsize, int n_derivs1, int n_derivs2,
                             int *index1, int *index2,
                             double lambda, double *lambda_l,
                             dcomplex m, dcomplex *m_l, dcomplex m2,
                             double *qx, double *qw, double *nr,
                             double *fv1, double *fv2, double *fv3, double *fv4,
                             double *qx2, double *qw2,
                             double **qx_l, double **qw_l, double **nr_l,
                             int n1, int n2, int n_qang);

void *lmie_core(lmie_core_threads_data *d);


#ifdef __cplusplus
}
#endif

#endif /* LMIE_CORE_H */
