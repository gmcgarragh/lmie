/*******************************************************************************
**
**    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef LMIE_SCHEDULE_H
#define LMIE_SCHEDULE_H

#ifdef __cplusplus
extern "C" {
#endif


int lmie_core_solution(int n_qsiz, int n_derivs1, int n_derivs2,
                       int save_control, int *index1, int *index2,
                       double lambda, double *lambda_l,
                       double mr, double mi, double *mr_l, double *mi_l,
                       double *qx, double *qw, double *nr,
                       double r1, double r2,
                       double *cext_, double *csca_, double *pbak_,
                       double *g, double **gc, double **lc,
                       double **qx_l, double **qw_l, double **nr_l,
                       double *cext_l_, double *csca_l_, double *pbak_l_,
                       double *g_l, double ***gc_l, double ***lc_l,
                       double ***save1, double ****save2, double accuracy,
                       int verbose, int n_threads, int use_mpi);

int lmie_core_solution_slave();


#ifdef __cplusplus
}
#endif

#endif /* LMIE_SCHEDULE_H */
