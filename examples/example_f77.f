c***********************************************************************
c This example program calls LMie to calculate the scattering properties
c of the accumulation mode mineral dust described by d'Almeida, 1991 at
c 2.25um.  A log normal size distribution is used with a mean radius of
c 0.39um and a standard deviation of 2.0um ranging from 0.005 to 50um.
c The real and imaginary parts of the index of refraction are 1.53 and
c 0.0054, respectively.  In this case the output includes size
c distribution statistics, extinction, scattering, and backscattering
c cross sections, asymmetry parameter, coefficients for the expansion of
c the scattering matrix in generalized spherical functions, and the
c elements of the normalized scattering matrix at 181 angles.  In
c addition, derivatives of the these output quantities with respect to
c mean radius, standard deviation, and the real and imaginary parts of
c the index of refraction are also generated.
c***********************************************************************
      program example_f77

      implicit none

      integer dummy_i

      integer i
      integer error

      logical calc_gc
      logical calc_lc
      logical calc_pf

      character*64 dist_name

c     integer n_int1
      integer n_int2
      integer n_quad
      integer n_angles
      parameter(n_angles = 181)
      integer n_derivs
      parameter(n_derivs = 4)
      integer max_coef
      parameter(max_coef = 336)
      integer max_coef2

      real*8 dummy_r
      real*8 dummy_r_l(n_derivs)

      real*8 lambda
      real*8 mr
      real*8 mi
      real*8 a1
      real*8 b1
c     real*8 a2
c     real*8 b2
c     real*8 gamma
      real*8 r1
      real*8 r2

      real*8 lambda_l(n_derivs)
      real*8 mr_l(n_derivs)
      real*8 mi_l(n_derivs)
      real*8 a1_l(n_derivs)
      real*8 b1_l(n_derivs)
c     real*8 a2_l(n_derivs)
c     real*8 b2_l(n_derivs)
c     real*8 gamma_l(n_derivs)
      real*8 r1_l(n_derivs)
      real*8 r2_l(n_derivs)

      real*8 accuracy

      integer n_coef

      real*8 r21
      real*8 r22
      real*8 norm
      real*8 reff
      real*8 veff
      real*8 gavg
      real*8 vavg
      real*8 ravg
      real*8 rvw
      real*8 cext
      real*8 csca
      real*8 cbak
      real*8 g
      real*8 gc(max_coef, 6)
      real*8 lc(max_coef, 6)
      real*8 theta(n_angles)
      real*8 pf(n_angles, 6)

      real*8 r21_l(n_derivs)
      real*8 r22_l(n_derivs)
      real*8 norm_l(n_derivs)
      real*8 reff_l(n_derivs)
      real*8 veff_l(n_derivs)
      real*8 gavg_l(n_derivs)
      real*8 vavg_l(n_derivs)
      real*8 ravg_l(n_derivs)
      real*8 rvw_l(n_derivs)
      real*8 cext_l(n_derivs)
      real*8 csca_l(n_derivs)
      real*8 cbak_l(n_derivs)
      real*8 g_l(n_derivs)
      real*8 gc_l(max_coef, 6, n_derivs)
      real*8 lc_l(max_coef, 6, n_derivs)
      real*8 pf_l(n_angles, 6, n_derivs)

      integer lmie_calc_max_coef_f77


c     ******************************************************************
c     * Fill the LMie input structure.  See the LMie documentation for
c     * details of each input quantity.
c     ******************************************************************
      calc_gc   = .true.
      calc_lc   = .false.
      calc_pf   = .true.

      dist_name = 'log_normal'

c     n_int1    =  			! not used for log normal
      n_int2    = 64
      n_quad    = 64
c     n_angles  = 181
c     n_derivs  = 4

      lambda    = 2.250d+00
      mr        = 1.530d+00
      mi        = 5.400d-03
      a1        = 3.900d-01
      b1        = 4.805d-01
c     a2        =  			! not used for log normal
c     b2        =  			! not used for log normal
c     gamma     =  			! not used for log normal
      r1        = 5.000d-03
      r2        = 5.000d+01

      accuracy  = 1.d-7


c     ******************************************************************
c     * Set linearized inputs.  The strategy here is to set all
c     * linearized inputs to zero then set the appropriate linearized
c     * inputs to unity.
c     ******************************************************************
      do i = 1, 4
           lambda_l(i) = 0.d0
           mr_l(i)     = 0.d0
           mi_l(i)     = 0.d0
           a1_l(i)     = 0.d0
           b1_l(i)     = 0.d0
c          a2_l(i)  			! not used for log normal
c          b2_l(i)  			! not used for log normal
c          gamma_l(i)  			! not used for log normal
           r1_l(i)     = 0.d0
           r2_l(i)     = 0.d0
      enddo

      mr_l(1)   = 1.d0                  ! derivative 0 is with respect to mr
      mi_l(2)   = 1.d0                  ! derivative 1 is with respect to mi
      a1_l(3)   = 1.d0                  ! derivative 2 is with respect to a1
      b1_l(4)   = 1.d0                  ! derivative 3 is with respect to b1


c     ******************************************************************
c     * Check to see if our max_coef is big enough.
c     ******************************************************************
      max_coef2 = lmie_calc_max_coef_f77(lambda, dist_name,
     &                                   a1, b1, dummy_r, r1, r2)
      if (max_coef2 .gt. max_coef) then
           write(0, *) 'max_coef too small, must be at least ',
     &                  max_coef2
           stop
      endif


c     ******************************************************************
c     * Call LMie to get the Mie solution.  In this case we indicate
c     * that verbose output is off, two threads are used if LMie was
c     * compiled with multi-threading support, and MPI should not be
c     * used.
c     ******************************************************************
      call lmie_solution_l_f77(calc_gc, calc_lc, calc_pf,
     &                         dist_name, dummy_i, n_int2, n_quad,
     &                         n_angles, n_derivs,
     &                         lambda, mr, mi,
     &                         a1, b1, dummy_r, dummy_r,
     &                         dummy_r, r1, r2,
     &                         lambda_l, mr_l, mi_l,
     &                         a1_l, b1_l, dummy_r_l, dummy_r_l,
     &                         dummy_r_l, r1_l, r2_l,
     &                         accuracy, n_coef,
     &                         r21, r22,
     &                         norm, reff, veff,
     &                         gavg, vavg, ravg, rvw,
     &                         cext, csca, cbak, g,
     &                         gc, lc, theta, pf,
     &                         r21_l, r22_l,
     &                         norm_l, reff_l, veff_l,
     &                         gavg_l, vavg_l, ravg_l, rvw_l,
     &                         cext_l, csca_l, cbak_l, g_l,
     &                         gc_l, lc_l, pf_l,
     &                         max_coef, .false., 2, .false., error)
      if (error /= 0) then
           write(0, *) 'lmie_solution_l_f77()'
           error = 1;
           stop
      endif


c     ******************************************************************
c     * Output results.
c     ******************************************************************
      write(*, '("----------------------------------------------------",
     &           "-----------------------")')
      write(*, '("full quantities (non linearized)")')
      write(*, '("----------------------------------------------------",
     &           "-----------------------")')

      call output_result(reff, veff,
     &                   gavg, vavg, ravg, rvw,
     &                   cext, csca, cbak, g,
     &                   n_coef, max_coef, gc,
     &                   n_angles, theta, pf)

      write(*, '()')

      do i = 1, n_derivs
           write(*, '("-----------------------------------------------",
     &                "----------------------------")')
           if (i .eq. 1) then
                write(*, '("derivatives wrt the real part of the index",
     &                     " of refraction (mr)")')
           else if (i .eq. 2) then
                write(*, '("derivatives wrt the imaginary part of the ",
     &                     "index of refraction (mi)")')
           else if (i .eq. 3) then
                write(*, '("derivatives wrt log normal size distributi",
     &                     "on mean radius (a1)")')
           else if (i .eq. 4) then
                write(*, '("derivatives wrt log normal size distributi",
     &                     "on standard deviation (b1)")')
           endif
           write(*, '("-----------------------------------------------",
     &                "----------------------------")')

           call output_result(reff_l(i), veff_l(i),
     &                        gavg_l(i), vavg_l(i), ravg_l(i), rvw_l(i),
     &                        cext_l(i), csca_l(i), cbak_l(i), g_l(i),
     &                        n_coef, max_coef, gc_l(1,1,i),
     &                        n_angles, theta, pf_l(1,1,i))

           write(*, '()')
      enddo

      end



c***********************************************************************
c Subroutine to output results.
c***********************************************************************
      subroutine output_result(reff, veff,
     &                         gavg, vavg, ravg, rvw,
     &                         cext, csca, cbak, g,
     &                         n_coef, max_coef, gc,
     &                         n_angles, theta, pf)

      implicit none

      real*8  reff
      real*8  veff
      real*8  gavg
      real*8  vavg
      real*8  ravg
      real*8  rvw
      real*8  cext
      real*8  csca
      real*8  cbak
      real*8  g
      integer n_coef
      integer max_coef
      real*8  gc(max_coef, 6)
      integer n_angles
      real*8  theta(n_angles)
      real*8  pf(n_angles, 6)

      integer i

      write(*, '("reff   = ", ES14.6)') reff
      write(*, '("veff   = ", ES14.6)') veff
      write(*, '("gavg   = ", ES14.6)') gavg
      write(*, '("vavg   = ", ES14.6)') vavg
      write(*, '("ravg   = ", ES14.6)') ravg
      write(*, '("rvw    = ", ES14.6)') rvw
      write(*, '("cext   = ", ES14.6)') cext
      write(*, '("csca   = ", ES14.6)') csca
      write(*, '("cbak   = ", ES14.6)') cbak
      write(*, '("g      = ", ES14.6)') g
      write(*, '("n_coef = ", I14)') n_coef
      write(*, '()')

      write(*, '("generalized spherical function coefficients:")')
      write(*, '("  alpha1        alpha2        alpha3        alpha4  ",
     &           "     -beta1        -beta2")')
      write(*, '("  beta          alpha         zeta          delta   ",
     &           "      gamma         epsilon")')
      write(*, '("----------------------------------------------------",
     &           "--------------------------------")')
      do i = 1, n_coef
           write(*, '(6ES14.6)') gc(i, 1), gc(i, 2), gc(i, 3), gc(i, 4),
     &                           gc(i, 5), gc(i, 6)
      enddo
      write(*, '()')

      write(*, '("elements of the normalized scattering matrix:")')
      write(*, '("  theta         F11           F22           F33     ",
     &           "      F44           F12           F34")')
      write(*, '("                a1            a2            a3      ",
     &           "      a4            b1            b2")')
      write(*, '("----------------------------------------------------",
     &           "--------------------------------------------")')
      do i = 1, n_angles
           write(*, '(7ES14.6)') theta(i), pf(i, 1), pf(i, 2), pf(i, 3),
     &                                     pf(i, 4), pf(i, 5), pf(i, 6)
      enddo

      end
