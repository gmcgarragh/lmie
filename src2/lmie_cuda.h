/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef LMIE_CUDA_H
#define LMIE_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif


int lmie_core_cuda(int n_qsiz, int n_derivs,
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


#ifdef __cplusplus
}
#endif

#endif /* LMIE_CUDA_H */
