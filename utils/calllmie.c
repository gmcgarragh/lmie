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

#include <rtutil_scat_io.h>
#include <rtutil_support.h>

#include <lmie_interface.h>


#ifdef NAME
#undef NAME
#endif
#define NAME    "calllmie"


typedef struct {
     int accuracy;
     int check_cbak;
     int check_g;
     int derivs;
     int findif;
     int dist;
     int dist_l;
     int gc;
     int lc;
     int pf;
     int help;
     int lambda;
     int lambda_l;
     int m;
     int m_l;
     int n_int1;
     int n_int2;
     int n_quad;
     int n_angles;
     int timing;
     int verbose;
     int version;
} options_data;


void usage();
void version();


int main(int argc, char *argv[]) {

     char *arg;

     char *fn_gc;
     char **fn_gc_l;
     char *fn_lc;
     char **fn_lc_l;
     char *fn_pf;
     char **fn_pf_l;

     int i;
     int n;
#ifdef USE_MPI
     int mpi_rank;
     int mpi_size;
#endif
     int n_args;

     enum size_dist_type dist_type;
     int n_int1;
     int n_int2;
     int n_quad;
     int n_angles;
     int n_derivs;

     int n_timing;

     int n_threads;

     int use_mpi;

     int n_coef;
/*
     double dx;
*/
     double check_tol_cbak;
     double check_tol_g;

     double cbak2;
     double g2;

     double lambda;
     double mr;
     double mi;
     double a1;
     double b1;
     double a2;
     double b2;
     double gamma;
     double r11;
     double r12;

     double *lambda_l;
     double *mr_l;
     double *mi_l;
     double *a1_l;
     double *b1_l;
     double *a2_l;
     double *b2_l;
     double *gamma_l;
     double *r11_l;
     double *r12_l;

     double accuracy;

     double r21;
     double r22;

     double norm;

     double reff;
     double veff;

     double area;
     double vol;
     double rmean;
     double rvw;

     double cext;
     double csca;
     double cbak;
     double g;

     double *r21_l;
     double *r22_l;

     double *norm_l;

     double *reff_l;
     double *veff_l;

     double *area_l;
     double *vol_l;
     double *rmean_l;
     double *rvw_l;

     double *cext_l;
     double *csca_l;
     double *cbak_l;
     double *g_l;

     double *zeros;

     double *theta;

     double **gc;
     double ***gc_l;

     double **lc;
     double ***lc_l;

     double **pf;
     double ***pf_l;

     FILE *fp;

     lmie_in_data  in;
     lmie_out_data out;
     lmie_out_data out2;

     options_data o;

     scat_header_data h;

     clock_t time1;
     clock_t time2;
/*
     struct timespec time1;
     struct timespec time2;
*/
#if PLATFORM == WIN32_MSVC
     _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     use_mpi = 0;

#ifdef USE_MPI
     MPI_Init(&argc, &argv);

     MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
     MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

     if (mpi_size != 1) {
          use_mpi = 1;

          if (mpi_rank != 0) {
               if (lmie_solution_slave()) {
                    eprintf("ERROR: lmie_solution_slave()\n");
                    exit(1);
               }

               MPI_Finalize();

               exit(0);
          }
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     o.check_cbak = 0;
     o.check_g    = 0;
     o.derivs     = 0;
     o.findif     = 0;
     o.dist       = 0;
     o.dist_l     = 0;
     o.gc         = 0;
     o.lc         = 0;
     o.pf         = 0;
     o.help       = 0;
     o.lambda     = 0;
     o.lambda_l   = 0;
     o.m          = 0;
     o.m_l        = 0;
     o.n_int1     = 0;
     o.n_int2     = 0;
     o.n_quad     = 0;
     o.n_angles   = 0;
     o.timing     = 0;
     o.verbose    = 0;
     o.version    = 0;

     check_tol_cbak = 1.e-5;
     check_tol_g    = 1.e-5;

     n_int1         = 0;
     n_derivs       = 0;

     accuracy       = 1.e-7;

     n_angles       = 18001;

     n_timing       = 1;

     n_threads      = 1;

     lambda         = 0.;
     mr             = 0.;
     mr             = 0.;
     a1             = 0.;
     b1             = 0.;
     a2             = 0.;
     b2             = 0.;
     gamma          = 0.;
     r11            = 0.;
     r12            = 0.;

     n = 0;
     for (i = 1; i < argc; ++i) {
          if (argv[i][0] == '-') {
               if (strcmp(argv[i], "-accuracy") == 0) {
                    o.accuracy = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    accuracy = strtod_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-check_cbak") == 0) {
                    o.check_cbak = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    check_tol_cbak = strtod_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-check_g") == 0) {
                    o.check_g = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    check_tol_g = strtod_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-derivs") == 0) {
                    o.derivs = 1;
                    o.findif = 0;
                    check_arg_count(i, argc, 1, "-derivs");
                    n_derivs = strtoi_errmsg_exit(argv[++i], "-derivs");

                    zeros = alloc_array1_d(n_derivs);
                    init_array1_d(zeros, n_derivs, 0.);
                    lambda_l = zeros;
                    mr_l     = zeros;
                    mr_l     = zeros;
                    a1_l     = zeros;
                    b1_l     = zeros;
                    a2_l     = zeros;
                    b2_l     = zeros;
                    gamma_l  = zeros;
                    r11_l    = zeros;
                    r12_l    = zeros;
               }
               else if (strcmp(argv[i], "-findif") == 0) {
                    o.derivs = 0;
                    o.findif = 1;
                    check_arg_count(i, argc, 1, "-findif");
/*
                    dx = strtod_errmsg_exit(argv[++i], "-findif");
*/
               }
               else if (strcmp(argv[i], "-dist") == 0) {
                    o.dist = 1;
                    check_arg_count(i, argc, 1, "-dist");
                    if (parse_dist_values(argc, argv, &i, argv[++i], &dist_type,
                                          &a1, &b1, &a2, &b2, &gamma, &r11, &r12, 1)) {
                         eprintf("ERROR: parse_dist_values()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-dist_l") == 0) {
                    if (! o.dist) {
                         eprintf("ERROR: must specify -dist before -dist_l\n");
                         exit(1);
                    }
                    o.dist_l = 1;
                    if (parse_dist_values_l(argc, argv, &i, dist_type, &a1_l, &b1_l,
                                &a2_l, &b2_l, &gamma_l, &r11_l, &r12_l, n_derivs, 1)) {
                         eprintf("ERROR: parse_dist_values_l()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-gc") == 0) {
                    o.gc = 1;
                    n_args = 1;
                    if (o.derivs || o.findif)
                         n_args++;
                    check_arg_count(i, argc, n_args, "-gc");
                    fn_gc = argv[++i];
                    if (o.derivs || o.findif) {
                         fn_gc_l = get_derivs_list_s(n_derivs, argv[i-1], argv[i+1]);
                         ++i;
                    }
               }
               else if (strcmp(argv[i], "-lc") == 0) {
                    o.lc = 1;
                    n_args = 1;
                    if (o.derivs || o.findif)
                         n_args++;
                    check_arg_count(i, argc, n_args, "-lc");
                    fn_lc = argv[++i];
                    if (o.derivs || o.findif) {
                         fn_lc_l = get_derivs_list_s(n_derivs, argv[i-1], argv[i+1]);
                         ++i;
                    }
               }
               else if (strcmp(argv[i], "-pf") == 0) {
                    o.pf = 1;
                    n_args = 1;
                    if (o.derivs || o.findif)
                         n_args++;
                    check_arg_count(i, argc, n_args, "-pf");
                    fn_pf = argv[++i];
                    if (o.derivs || o.findif) {
                         fn_pf_l = get_derivs_list_s(n_derivs, argv[i-1], argv[i+1]);
                         ++i;
                    }
               }
               else if (strcmp(argv[i], "-help") == 0) {
                    usage();
                    exit(0);
               }
               else if (strcmp(argv[i], "-lambda") == 0) {
                    o.lambda = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    lambda = strtod_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-lambda_l") == 0) {
                    o.lambda_l = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    lambda_l = get_derivs_list_d(n_derivs, arg, argv[i+1]);
                    i++;
               }
               else if (strcmp(argv[i], "-m") == 0) {
                    o.m = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 2, arg);
                    mr = strtod_errmsg_exit(argv[++i], arg);
                    mi = strtod_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-m_l") == 0) {
                    o.m_l = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 2, arg);
                    mr_l = get_derivs_list_d(n_derivs, arg, argv[i+1]);
                    mi_l = get_derivs_list_d(n_derivs, arg, argv[i+2]);
                    i += 2;
               }
               else if (strcmp(argv[i], "-n_int1") == 0) {
                    o.n_int1 = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    n_int1 = strtoi_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-n_int2") == 0) {
                    o.n_int2 = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    n_int2 = strtoi_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-n_quad") == 0) {
                    o.n_quad = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    n_quad = strtoi_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-n_angles") == 0) {
                    o.n_angles = 1;
                    arg = argv[i];
                    check_arg_count(i, argc, 1, arg);
                    n_angles = strtoi_errmsg_exit(argv[++i], arg);
               }
               else if (strcmp(argv[i], "-n_threads") == 0) {
                    check_arg_count(i, argc, 1, "-n_threads");
                    n_threads = strtoi_errmsg_exit(argv[++i], "-n_threads");
               }
               else if (strcmp(argv[i], "-timing") == 0) {
                    o.timing = 1;
                    check_arg_count(i, argc, 1, "-timing");
                    n_timing = strtoi_errmsg_exit(argv[++i], "-timing");
               }
               else if (strcmp(argv[i], "-verbose") == 0)
                    o.verbose = 1;
               else if (strcmp(argv[i], "-version") == 0) {
                    version();
                    exit(0);
               }
               else {
                    printf("Invalid option: %s, use -help for more information\n", argv[i]);
                    exit(1);
               }
          }
          else {
/*
               if (n == 0)
                    filename = argv[i];
               else {
*/
                    printf("Too many arguments, use -help for more information\n");
                    exit(1);
/*
               }
*/
               ++n;
          }
     }
/*
     if (n < 0) {
          printf("Not enough arguments, use -help for more information\n");
          exit(1);
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! o.dist) {
          eprintf("ERROR: must specify a size distribution with -dist\n");
          exit(1);
     }

     if (! o.lambda) {
          eprintf("ERROR: must specify the wavelength with -lambda\n");
          exit(1);
     }

     if (! o.m) {
          eprintf("ERROR: must specify the index of refraction with -m\n");
          exit(1);
     }

     if (dist_type == SIZE_DIST_MODIFIED_POWER_LAW && ! o.n_int1) {
          eprintf("ERROR: must specify the # of integration subintervals from "
                 "r0 to r1 with -np (-n01)\n");
          exit(1);
     }

     if (! o.n_int2) {
          eprintf("ERROR: must specify the # of integration subintervals from "
                 "r1 to r2 with -n (-n12)\n");
          exit(1);
     }

     if (! o.n_quad) {
          eprintf("ERROR: must specify the # of guassian division points on "
                 "each interval with -nk\n");
          exit(1);
     }

     if (! o.gc && ! o.lc && ! o.pf) {
          eprintf("ERROR: at least one of -gc -lc or -pf should be used\n");
          exit(1);
     }

     if (o.check_cbak && ! o.pf) {
          eprintf("ERROR: must use -pf when using -check_cbak\n");
          exit(1);
     }

     if (o.check_g && ! o.pf) {
          eprintf("ERROR: must use -gc or -lc when using -check_g\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     lmie_in_alloc(&in, n_derivs);

     in.calc_gc   = o.gc || o.pf;
     in.calc_lc   = o.lc;
     in.calc_pf   = o.pf;

     in.dist_type = dist_type;

     in.n_int1    = n_int1;
     in.n_int2    = n_int2;
     in.n_quad    = n_quad;
     in.n_angles  = n_angles;
     in.n_derivs  = n_derivs;

     in.lambda    = lambda;
     in.mr        = mr;
     in.mi        = mi;
     in.a1        = a1;
     in.b1        = b1;
     in.a2        = a2;
     in.b2        = b2;
     in.gamma     = gamma;
     in.r1        = r11;
     in.r2        = r12;

     for (i = 0; i < n_derivs; ++i) {
          in.lambda_l[i]  = lambda_l[i];
          in.mr_l[i]      = mr_l[i];
          in.mi_l[i]      = mi_l[i];
          in.a1_l[i]      = a1_l[i];
          in.b1_l[i]      = b1_l[i];
          in.a2_l[i]      = a2_l[i];
          in.b2_l[i]      = b2_l[i];
          in.gamma_l[i]   = gamma_l[i];
          in.r1_l[i]      = r11_l[i];
          in.r2_l[i]      = r12_l[i];
     }

     in.accuracy  = accuracy;

     if (lmie_solution(&in, &out, 1, o.verbose, n_threads, use_mpi)) {
          eprintf("ERROR: lmie_solution()\n");
          exit(1);
     }

     if (o.timing) {

          time1 = clock();
/*
          clock_gettime(CLOCK_MONOTONIC, &time1);
*/
          for (i = 0; i < n_timing; ++i) {
               if (lmie_solution(&in, &out2, 1, o.verbose, n_threads, use_mpi)) {
                    eprintf("ERROR: lmie_solution()\n");
                    exit(1);
               }
          }

          lmie_out_free(&out2, o.derivs);

          time2 = clock();
/*
          clock_gettime(CLOCK_MONOTONIC, &time2);
*/
     }


     n_coef  = out.n_coef;

     r21     = out.r1;
     r22     = out.r2;
     norm    = out.norm;
     reff    = out.reff;
     veff    = out.veff;
     area    = out.gavg;
     vol     = out.vavg;
     rmean   = out.ravg;
     rvw     = out.rvw;
     cext    = out.cext;
     csca    = out.csca;
     cbak    = out.cbak;
     g       = out.g;
     gc      = out.gc;
     lc      = out.lc;
     theta   = out.theta;
     pf      = out.pf;

     r21_l   = out.r1_l;
     r22_l   = out.r2_l;
     norm_l  = out.norm_l;
     reff_l  = out.reff_l;
     veff_l  = out.veff_l;
     area_l  = out.gavg_l;
     vol_l   = out.vavg_l;
     rmean_l = out.ravg_l;
     rvw_l   = out.rvw_l;
     cext_l  = out.cext_l;
     csca_l  = out.csca_l;
     cbak_l  = out.cbak_l;
     g_l     = out.g_l;
     gc_l    = out.gc_l;
     lc_l    = out.lc_l;
     pf_l    = out.pf_l;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.check_cbak) {
          cbak2 = pf[0][n_angles - 1] * csca / (4. * PI);

          if (fabs(cbak - cbak2) > check_tol_cbak) {
               eprintf("ERROR: check_g failed: cbak = %.6e and pf[0][n_angles - 1] * "
                      "csca / (4. * PI) = %.6e\n", cbak, cbak2);
               exit(1);
          }
     }

     if (o.check_g) {
          if (o.gc || o.pf)
               g2 = gc[0][1] / 3.;
          else
               g2 = lc[0][1] / 3.;

          if (fabs(g - g2) > check_tol_g) {
               eprintf("ERROR: check_g failed: g = %.6e and gc[0][1] / 3. = %.6e\n", g, g2);
               exit(1);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     scat_header_import(&h, lambda, mr, mi, dist_type, a1, b1, a2, b2, gamma,
                        r21, r22, n_int1, n_int2, n_quad, accuracy, norm,
                        reff, veff, area, vol, rmean, rvw, cext, csca, cbak, g);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.gc) {
          if (fn_gc[0] == '-')
               fp = stdout;
          else {
               if ((fp = fopen(fn_gc, "w")) == NULL) {
                    eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                            fn_gc, strerror(errno));
                    exit(1);
               }
          }

          if (scat_coefs_write_fp(fp, n_coef, &h, gc, 1)) {
               eprintf("ERROR: scat_coefs_write_fp()\n");
               exit(1);
          }

          if (fn_gc[0] != '-')
               fclose(fp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.lc) {
          if (fn_lc[0] == '-')
               fp = stdout;
          else {
               if ((fp = fopen(fn_lc, "w")) == NULL) {
                    eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                            fn_lc, strerror(errno));
                    exit(1);
               }
          }

          if (scat_coefs_write_fp(fp, n_coef, &h, lc, 1)) {
               eprintf("ERROR: scat_coefs_write_fp()\n");
               exit(1);
          }

          if (fn_lc[0] != '-')
               fclose(fp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.pf) {
          if (fn_pf[0] == '-')
               fp = stdout;
          else {
               if ((fp = fopen(fn_pf, "w")) == NULL) {
                    eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                            fn_pf, strerror(errno));
                    exit(1);
               }
          }

          if (phase_func_write_fp(fp, &h, theta, pf, n_angles, 1)) {
               eprintf("ERROR: phase_func_write_fp()\n");
               exit(1);
          }

          if (fn_pf[0] != '-')
               fclose(fp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_derivs; ++i) {

          scat_header_import(&h, lambda_l[i], mr_l[i], mi_l[i], dist_type,
                             a1_l[i], b1_l[i], a2_l[i], b2_l[i], gamma_l[i],
                             r21_l[i], r22_l[i], n_int1, n_int2, n_quad, accuracy,
                             norm_l[i], reff_l[i], veff_l[i], area_l[i], vol_l[i],
                             rmean_l[i], rvw_l[i], cext_l[i], csca_l[i], cbak_l[i],
                             g_l[i]);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (o.gc) {
               if (fn_gc_l[i][0] == '-')
                    fp = stdout;
               else {
                    if ((fp = fopen(fn_gc_l[i], "w")) == NULL) {
                         eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                                 fn_gc_l[i], strerror(errno));
                         exit(1);
                    }
               }

               if (scat_coefs_write_fp(fp, n_coef, &h, gc_l[i], 1)) {
                    eprintf("ERROR: scat_coefs_write_fp()\n");
                    exit(1);
               }

               if (fn_gc_l[i][0] != '-')
                    fclose(fp);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (o.lc) {
               if (fn_lc_l[i][0] == '-')
                    fp = stdout;
               else {
                    if ((fp = fopen(fn_lc_l[i], "w")) == NULL) {
                         eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                                 fn_lc_l[i], strerror(errno));
                         exit(1);
                    }
               }

               if (scat_coefs_write_fp(fp, n_coef, &h, lc_l[i], 1)) {
                    eprintf("ERROR: scat_coefs_write_fp()\n");
                    exit(1);
               }

               if (fn_lc_l[i][0] != '-')
                    fclose(fp);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (o.pf) {
               if (fn_pf_l[i][0] == '-')
                    fp = stdout;
               else {
                    if ((fp = fopen(fn_pf_l[i], "w")) == NULL) {
                         eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                                 fn_pf_l[i], strerror(errno));
                         exit(1);
                    }
               }

               if (phase_func_write_fp(fp, &h, theta, pf_l[i], n_angles, 1)) {
                    eprintf("ERROR: phase_func_write_fp()\n");
                    exit(1);
               }

               if (fn_pf_l[i][0] != '-')
                    fclose(fp);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.timing)
          printf("%.2f\n", (double) (time2 - time1) / (double) CLOCKS_PER_SEC);
/*
          printf("%.2f\n", (double) (time2.tv_sec - time1.tv_sec) + (double) (time2.tv_nsec - time1.tv_nsec) / 1.e9);
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.derivs) {
          if (o.gc) {
               for (i = 0; i < n_derivs; ++i)
                    free(fn_gc_l[i]);
               free((void *) fn_gc_l);
          }

          if (o.lc) {
               for (i = 0; i < n_derivs; ++i)
                    free(fn_lc_l[i]);
               free((void *) fn_lc_l);
          }

          if (o.pf) {
               for (i = 0; i < n_derivs; ++i)
                    free(fn_pf_l[i]);
               free((void *) fn_pf_l);
          }

          if (lambda_l != zeros)
               free(lambda_l);
          if (mr_l != zeros)
               free(mr_l);
          if (mi_l != zeros)
               free(mi_l);
          if (a1_l != zeros)
               free(a1_l);
          if (b1_l != zeros)
               free(b1_l);
          if (a2_l != zeros)
               free(a2_l);
          if (b2_l != zeros)
               free(b2_l);
          if (gamma_l != zeros)
               free(gamma_l);
          if (r11_l  != zeros)
               free(r11_l);
          if (r12_l  != zeros)
               free(r12_l);

          free_array1_d(zeros);
     }

     lmie_in_free(&in, o.derivs);
     lmie_out_free(&out, o.derivs);

#ifdef USE_MPI
     MPI_Finalize();
#endif

     exit(0);
}



void usage() {

     printf("%s <flags ...>\n", NAME);
     printf("\n");

     printf("See the LMie documentation for details on running %s.\n", NAME);
}


void version() {

     printf("%s version %s\n", NAME, VERSION);
}
