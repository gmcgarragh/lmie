/*******************************************************************************
 * This example program calls LMie to calculate the scattering properties of the
 * accumulation mode mineral dust described by d'Almeida, 1991 at 2.25um.  A log
 * normal size distribution is used with a mean radius of 0.39um and a standard
 * deviation of 2.0um ranging from 0.005 to 50um.  The real and imaginary parts
 * of the index of refraction are 1.53 and 0.0054, respectively.  In this case
 * the output includes size distribution statistics, extinction, scattering, and
 * backscattering cross sections, asymmetry parameter, coefficients for the
 * expansion of the scattering matrix in generalized spherical functions, and
 * and the elements of the normalized scattering matrix at 181 angles.  In
 * addition, derivatives of the these output quantities with respect to mean
 * radius, standard deviation, and the real and imaginary parts of the index of
 * refraction are also generated.
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

/* Include the LMie C interface. */
#include <lmie_interface.h>


void output_results(double reff, double veff, double gavg, double vavg,
                    double ravg, double rvw, double cext, double csca,
                    double cbak, double g, int n_coef, double **gc,
                    int n_angles, double *theta, double **pf);


/*******************************************************************************
 * Main program.
 ******************************************************************************/
int main(int argc, char *argv[]) {

     int i;

     lmie_in_data  in;		/* LMie input structure */
     lmie_out_data out;		/* LMie output structure */


     /*-------------------------------------------------------------------------
      * Fill the LMie input structure.  See the LMie documentation for details
      * of each input quantity.
      *-----------------------------------------------------------------------*/

     /* Allocate memory used in the input structure for 4 derivatives */
     lmie_in_alloc(&in, 4);

     in.calc_gc   = 1;
     in.calc_lc   = 0;
     in.calc_pf   = 1;

     in.dist_type = SIZE_DIST_LOG_NORMAL;

     in.n_int1    = 0;			/* not used for log normal */
     in.n_int2    = 64;
     in.n_quad    = 64;
     in.n_angles  = 181;
     in.n_derivs  = 4;

     in.lambda    = 2.250e+00;
     in.mr        = 1.530e+00;
     in.mi        = 5.400e-03;
     in.a1        = 3.900e-01;
     in.b1        = 4.805e-01;
     in.a2        = 0.;			/* not used for log normal */
     in.b2        = 0.;			/* not used for log normal */
     in.gamma     = 0.;			/* not used for log normal */
     in.r1        = 5.000e-03;
     in.r2        = 5.000e+01;

     in.accuracy  = 1.e-7;


     /*-------------------------------------------------------------------------
      * Set linearized inputs.  The strategy here is to set all linearized
      * inputs to zero then set the appropriate linearized inputs to unity.
      *-----------------------------------------------------------------------*/
     for (i = 0; i < in.n_derivs; ++i) {
          in.lambda_l[i] = 0.;
          in.mr_l[i]     = 0.;
          in.mi_l[i]     = 0.;
          in.a1_l[i]     = 0.;
          in.b1_l[i]     = 0.;
          in.a2_l[i]     = 0.;		/* not used for log normal */
          in.b2_l[i]     = 0.;		/* not used for log normal */
          in.gamma_l[i]  = 0.;		/* not used for log normal */
          in.r1_l[i]     = 0.;
          in.r2_l[i]     = 0.;
     }

     /* Or more conveniently just call the function. */
     lmie_in_zero_derivs(&in, in.n_derivs);

     in.mr_l[0]   = 1.;			/* derivative 0 is with respect to mr */
     in.mi_l[1]   = 1.;			/* derivative 1 is with respect to mi */
     in.a1_l[2]   = 1.;			/* derivative 2 is with respect to a1 */
     in.b1_l[3]   = 1.;			/* derivative 3 is with respect to b1 */


     /*-------------------------------------------------------------------------
      * Call LMie to get the Mie solution.  In this case we indicate that memory
      * used in the output structure should be allocated, verbose output is off,
      * two threads are used if LMie was compiled with multi-threading support,
      * and MPI should not be used.
      *-----------------------------------------------------------------------*/
     if (lmie_solution(&in, &out, 1, 0, 2, 0)) {
          fprintf(stderr, "ERROR: lmie_solution()\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      * Output results.
      *-----------------------------------------------------------------------*/
#if PLATFORM == WIN32_MSVC
     _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
     printf("---------------------------------------------------------------------------\n");
     printf("full quantities (non linearized)\n");
     printf("---------------------------------------------------------------------------\n");

     output_results(out.reff, out.veff, out.gavg, out.vavg,
                    out.ravg, out.rvw, out.cext, out.csca,
                    out.cbak, out.g, out.n_coef, out.gc,
                    in.n_angles, out.theta, out.pf);
     printf("\n");

     for (i = 0; i < in.n_derivs; ++i) {
          printf("---------------------------------------------------------------------------\n");
          if (i == 0)
               printf("derivatives wrt the real part of the index of refraction (mr)\n");
          else if (i == 1)
               printf("derivatives wrt the imaginary part of the index of refraction (mi)\n");
          else if (i == 2)
               printf("derivatives wrt log normal size distribution mean radius (a1)\n");
          else if (i == 3)
               printf("derivatives wrt log normal size distribution standard deviation (b1)\n");
          printf("---------------------------------------------------------------------------\n");

          output_results(out.reff_l[i], out.veff_l[i], out.gavg_l[i], out.vavg_l[i],
                         out.ravg_l[i], out.rvw_l[i], out.cext_l[i], out.csca_l[i],
                         out.cbak_l[i], out.g_l[i], out.n_coef, out.gc_l[i],
                         in.n_angles, out.theta, out.pf_l[i]);

          printf("\n");
     }


     /*-------------------------------------------------------------------------
      * Free memory allocated in the input and output structures.
      *-----------------------------------------------------------------------*/
     lmie_in_free (&in,  1);
     lmie_out_free(&out, 1);
 

     exit(0);
}



/*******************************************************************************
 * Function to output results.
 ******************************************************************************/
void output_results(double reff, double veff, double gavg, double vavg,
                    double ravg, double rvw, double cext, double csca,
                    double cbak, double g, int n_coef, double **gc,
                    int n_angles, double *theta, double **pf) {

     int i;
     int j;

     printf("reff   = %14E\n",   reff);
     printf("veff   = %14E\n",   veff);
     printf("gavg   = %14E\n",   gavg);
     printf("vavg   = %14E\n",   vavg);
     printf("ravg   = %14E\n",   ravg);
     printf("rvw    = %14E\n",   rvw);
     printf("cext   = %14E\n",   cext);
     printf("csca   = %14E\n",   csca);
     printf("cbak   = %14E\n",   cbak);
     printf("g      = %14E\n",   g);
     printf("n_coef = %14d\n", n_coef);
     printf("\n");

     printf("generalized spherical function coefficients:\n");
     printf("  alpha1        alpha2        alpha3        alpha4       -beta1        -beta2\n");
     printf("  beta          alpha         zeta          delta         gamma         epsilon\n");
     printf("------------------------------------------------------------------------------------\n");
     for (i = 0; i < n_coef; ++i) {
          for (j = 0; j < 6; ++j) {
               printf("% 14E", gc[j][i]);
          }
          printf("\n");
     }
     printf("\n");

     printf("elements of the normalized scattering matrix:\n");
     printf("  theta         F11           F22           F33           F44           F12           F34\n");
     printf("                a1            a2            a3            a4            b1            b2\n");
     printf("------------------------------------------------------------------------------------------------\n");
     for (i = 0; i < n_angles; ++i) {
          printf("%14E", theta[i]);

          for (j = 0; j < 6; ++j) {
               printf("% 14E", pf[j][i]);
          }
          printf("\n");
     }
}
