c*******************************************************************************
c
c    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
c
c    This source code is licensed under the GNU General Public License (GPL),
c    Version 3.  See the file COPYING for more details.
c
c*******************************************************************************


c***********************************************************************
c
c***********************************************************************
      subroutine terminate_string(s)

      implicit none

      character*64 s

      integer i

      i = 64
      do while (i .gt. 0 .and. s(i:i) .eq. ' ')
           i = i - 1
      enddo

      i = i + 1

      s(i:i) = achar(0)

      end



c***********************************************************************
c
c***********************************************************************
      integer function lmie_calc_max_coef_f77(lambda, dist_name, a1, a2, a3,
     &                                        r1, r2)

      implicit none

      real*8 lambda

      character*64 dist_name

      real*8 a1
      real*8 a2
      real*8 a3
      real*8 r1
      real*8 r2

      integer dist_type

      integer size_dist_code
      integer lmie_calc_max_coef

      call terminate_string(dist_name)
      dist_type = size_dist_code(dist_name);

      lmie_calc_max_coef_f77 =
     &     lmie_calc_max_coef(lambda, dist_type, a1, a2, a3, r1, r2)

      end



c***********************************************************************
c
c***********************************************************************
      subroutine lmie_in_zero_derivs_f77(n_derivs,
     &                                   lambda_l, mr_l, mi_l,
     &                                   a1_l, a2_l, a3_l, a4_l, a5_l,
     &                                   r1_l, r2_l)

      implicit none

      integer n_derivs

      real*8 lambda_l(n_derivs)
      real*8 mr_l(n_derivs)
      real*8 mi_l(n_derivs)
      real*8 a1_l(n_derivs)
      real*8 a2_l(n_derivs)
      real*8 a3_l(n_derivs)
      real*8 a4_l(n_derivs)
      real*8 a5_l(n_derivs)
      real*8 r1_l(n_derivs)
      real*8 r2_l(n_derivs)

      integer i

      do i = 1, n_derivs
           lambda_l(i) = 0.d0
           mr_l(i)     = 0.d0
           mi_l(i)     = 0.d0
           a1_l(i)     = 0.d0
           a2_l(i)     = 0.d0
           a3_l(i)     = 0.d0
           a4_l(i)     = 0.d0
           a5_l(i)     = 0.d0
           r1_l(i)     = 0.d0
           r2_l(i)     = 0.d0
      enddo

      end



c***********************************************************************
c
c***********************************************************************
      subroutine lmie_solution_f77(calc_gc, calc_lc, calc_pf,
     &                             dist_name, n_int1, n_int2, n_quad,
     &                             n_angles, save_control,
     &                             lambda, mr, mi,
     &                             a1, a2, a3, a4, a5,
     &                             r1, r2,
     &                             accuracy, n_coef,
     &                             r21, r22,
     &                             norm, reff, veff,
     &                             gavg, vavg, ravg, rvw,
     &                             cext, csca, cbak, g,
     &                             gc, lc, theta, pf,
     &                             save1, save2,
     &                             max_coef, verbose, n_threads,
     &                             use_mpi, error)

      implicit none

      integer calc_gc
      integer calc_lc
      integer calc_pf

      character*64 dist_name

      integer n_int1
      integer n_int2
      integer n_quad
      integer n_angles

      integer save_control

      real*8 lambda
      real*8 mr
      real*8 mi
      real*8 a1
      real*8 a2
      real*8 a3
      real*8 a4
      real*8 a5
      real*8 r1
      real*8 r2

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

      integer*8 save1
      integer*8 save2

      integer max_coef
      logical verbose
      integer n_threads
      logical use_mpi

      integer error


      integer ret_val

      integer dist_type
      integer use_mpi2

      integer size_dist_code
      integer lmie_solution2


c     ******************************************************************
c     *
c     ******************************************************************
      call terminate_string(dist_name)
      dist_type = size_dist_code(dist_name);

      if (use_mpi) then
           use_mpi2 = 1
      else
           use_mpi2 = 0
      endif

      ret_val = lmie_solution2(calc_gc, calc_lc, calc_pf,
     &                         dist_type, n_int1, n_int2, n_quad,
     &                         n_angles, save_control,
     &                         lambda, mr, mi,
     &                         a1, a2, a3, a4, a5,
     &                         r1, r2,
     &                         accuracy, n_coef,
     &                         r21, r22,
     &                         norm, reff, veff,
     &                         gavg, vavg, ravg, rvw,
     &                         cext, csca, cbak, g,
     &                         gc, lc, theta, pf,
     &                         save1, save2,
     &                         max_coef, verbose, n_threads, use_mpi)
      if (ret_val .ne. 0) then
           write(0, *) 'lmie_solution2()'
           error = -1;
           return
      endif

      error = 0;

      end


c***********************************************************************
c
c***********************************************************************
      subroutine lmie_solution_l_f77(calc_gc, calc_lc, calc_pf,
     &                               dist_name, n_int1, n_int2, n_quad,
     &                               n_angles, n_derivs, save_control,
     &                               lambda, mr, mi,
     &                               a1, a2, a3, a4, a5,
     &                               r1, r2,
     &                               lambda_l, mr_l, mi_l,
     &                               a1_l, a2_l, a3_l, a4_l, a5_l,
     &                               r1_l, r2_l,
     &                               accuracy, n_coef,
     &                               r21, r22,
     &                               norm, reff, veff,
     &                               gavg, vavg, ravg, rvw,
     &                               cext, csca, cbak, g,
     &                               gc, lc, theta, pf,
     &                               r21_l, r22_l,
     &                               norm_l, reff_l, veff_l,
     &                               gavg_l, vavg_l, ravg_l, rvw_l,
     &                               cext_l, csca_l, cbak_l, g_l,
     &                               gc_l, lc_l, pf_l,
     &                               save1, save2,
     &                               max_coef, verbose, n_threads,
     &                               use_mpi, error)

      implicit none

      logical calc_gc
      logical calc_lc
      logical calc_pf

      character*64 dist_name

      integer n_int1
      integer n_int2
      integer n_quad
      integer n_angles
      integer n_derivs

      integer save_control

      real*8 lambda
      real*8 mr
      real*8 mi
      real*8 a1
      real*8 a2
      real*8 a3
      real*8 a4
      real*8 a5
      real*8 r1
      real*8 r2

      real*8 lambda_l(n_derivs)
      real*8 mr_l(n_derivs)
      real*8 mi_l(n_derivs)
      real*8 a1_l(n_derivs)
      real*8 a2_l(n_derivs)
      real*8 a3_l(n_derivs)
      real*8 a4_l(n_derivs)
      real*8 a5_l(n_derivs)
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

      integer*8 save1
      integer*8 save2

      integer max_coef
      logical verbose
      integer n_threads
      logical use_mpi

      integer error


      integer ret_val

      integer calc_gc2
      integer calc_lc2
      integer calc_pf2

      integer dist_type

      integer use_mpi2

      integer size_dist_code
      integer lmie_solution2_l


c     ******************************************************************
c     *
c     ******************************************************************
      call terminate_string(dist_name)
      dist_type = size_dist_code(dist_name);

      if (calc_gc) then
           calc_gc2 = 1
      else
           calc_gc2 = 0
      endif

      if (calc_lc) then
           calc_lc2 = 1
      else
           calc_lc2 = 0
      endif

      if (calc_pf) then
           calc_pf2 = 1
      else
           calc_pf2 = 0
      endif

      if (use_mpi) then
           use_mpi2 = 1
      else
           use_mpi2 = 0
      endif

      ret_val = lmie_solution2_l(calc_gc2, calc_lc2, calc_pf2,
     &                           dist_type, n_int1, n_int2, n_quad,
     &                           n_angles, n_derivs, save_control,
     &                           lambda, mr, mi,
     &                           a1, a2, a3, a4, a5,
     &                           r1, r2,
     &                           lambda_l, mr_l, mi_l,
     &                           a1_l, a2_l, a3_l, a4_l, a5_l,
     &                           r1_l, r2_l,
     &                           accuracy, n_coef,
     &                           r21, r22,
     &                           norm, reff, veff,
     &                           gavg, vavg, ravg, rvw,
     &                           cext, csca, cbak, g,
     &                           gc, lc, theta, pf,
     &                           r21_l, r22_l,
     &                           norm_l, reff_l, veff_l,
     &                           gavg_l, vavg_l, ravg_l, rvw_l,
     &                           cext_l, csca_l, cbak_l, g_l,
     &                           gc_l, lc_l, pf_l,
     &                           save1, save2,
     &                           max_coef, verbose, n_threads, use_mpi2)
      if (ret_val .ne. 0) then
           write(0, *) 'lmie_solution2()'
           error = -1;
           return
      endif

      error = 0;

      end
