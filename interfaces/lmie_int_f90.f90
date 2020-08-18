!*******************************************************************************
!
!    Copyright (C) 2008-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
!
!    This source code is licensed under the GNU General Public License (GPL),
!    Version 3.  See the file COPYING for more details.
!
!*******************************************************************************

module lmie_int_f90

use iso_c_binding


implicit none


private

public :: lmie_in_type, &
          lmie_out_type, &
          lmie_in_allocate_f90, &
          lmie_in_deallocate_f90, &
          lmie_out_allocate_f90, &
          lmie_out_deallocate_f90, &
          lmie_in_zero_derivs_f90, &
          lmie_solution_f90


real(8), target :: zero1(1)     = reshape((/ 0. /), shape(zero1))
real(8), target :: zero2(1,1)   = reshape((/ 0. /), shape(zero2))
real(8), target :: zero3(1,1,1) = reshape((/ 0. /), shape(zero3))


!*******************************************************************************
!
!*******************************************************************************
integer, parameter, public :: DIST_NONE                      = 0
integer, parameter, public :: DIST_MODIFIED_GAMMA            = 1
integer, parameter, public :: DIST_LOG_NORMAL                = 2
integer, parameter, public :: DIST_POWER_LAW                 = 3
integer, parameter, public :: DIST_GAMMA                     = 4
integer, parameter, public :: DIST_MODIFIED_POWER_LAW        = 5
integer, parameter, public :: DIST_BIMODAL_VOLUME_LOG_NORMAL = 6
integer, parameter, public :: DIST_EXPONENTIAL               = 7


!*******************************************************************************
!
!*******************************************************************************
type lmie_in_type
     integer :: save_control

     logical :: calc_gc
     logical :: calc_lc
     logical :: calc_pf

     integer :: dist_type

     integer :: n_int1
     integer :: n_int2
     integer :: n_quad
     integer :: n_angles
     integer :: n_derivs

     real(8) :: lambda
     real(8) :: mr
     real(8) :: mi
     real(8) :: a1
     real(8) :: a2
     real(8) :: a3
     real(8) :: a4
     real(8) :: a5
     real(8) :: r1
     real(8) :: r2

     real(8), pointer :: lambda_l(:)
     real(8), pointer :: mr_l(:)
     real(8), pointer :: mi_l(:)
     real(8), pointer :: a1_l(:)
     real(8), pointer :: a2_l(:)
     real(8), pointer :: a3_l(:)
     real(8), pointer :: a4_l(:)
     real(8), pointer :: a5_l(:)
     real(8), pointer :: r1_l(:)
     real(8), pointer :: r2_l(:)

     real(8) :: accuracy
end type lmie_in_type


type lmie_out_type
     integer :: n_coef

     real(8) :: r1
     real(8) :: r2
     real(8) :: norm
     real(8) :: reff
     real(8) :: veff
     real(8) :: gavg
     real(8) :: vavg
     real(8) :: ravg
     real(8) :: rvw
     real(8) :: cext
     real(8) :: csca
     real(8) :: cbak
     real(8) :: g
     real(8), pointer :: gc(:,:)
     real(8), pointer :: lc(:,:)
     real(8), pointer :: theta(:)
     real(8), pointer :: pf(:,:)

     real(8), pointer :: r1_l(:)
     real(8), pointer :: r2_l(:)
     real(8), pointer :: norm_l(:)
     real(8), pointer :: reff_l(:)
     real(8), pointer :: veff_l(:)
     real(8), pointer :: gavg_l(:)
     real(8), pointer :: vavg_l(:)
     real(8), pointer :: ravg_l(:)
     real(8), pointer :: rvw_l(:)
     real(8), pointer :: cext_l(:)
     real(8), pointer :: csca_l(:)
     real(8), pointer :: cbak_l(:)
     real(8), pointer :: g_l(:)
     real(8), pointer :: gc_l(:,:,:)
     real(8), pointer :: lc_l(:,:,:)
     real(8), pointer :: pf_l(:,:,:)

     type(c_ptr)      :: save1
     type(c_ptr)      :: save2
end type lmie_out_type


!*******************************************************************************
!
!*******************************************************************************
interface
     integer(c_int) function lmie_calc_max_coef(lambda, dist_type, a1, a2, a3, &
                                                r1, r2) bind(c)

     use iso_c_binding

     implicit none

     real(c_double), intent(in) :: lambda
     integer(c_int), intent(in) :: dist_type
     real(c_double), intent(in) :: a1
     real(c_double), intent(in) :: a2
     real(c_double), intent(in) :: a3
     real(c_double), intent(in) :: r1
     real(c_double), intent(in) :: r2

     end function lmie_calc_max_coef
end interface


!*******************************************************************************
!
!*******************************************************************************
interface
     integer(c_int) function lmie_solution2_l(calc_gc, calc_lc, calc_pf, &
                                       dist_type, n_int1, n_int2, n_quad, &
                                       n_angles, n_derivs, save_control, &
                                       lambda, mr, mi, &
                                       a1, a2, a3, a4, a5, &
                                       r1, r2, &
                                       lambda_l, mr_l, mi_l, &
                                       a1_l, a2_l, a3_l, a4_l, a5_l, &
                                       r1_l, r2_l, &
                                       accuracy, n_coef, &
                                       r21, r22, &
                                       norm, reff, veff, &
                                       gavg, vavg, ravg, rvw, &
                                       cext, csca, cbak, g, &
                                       gc, lc, theta, pf, &
                                       r21_l, r22_l, &
                                       norm_l, reff_l, veff_l, &
                                       gavg_l, vavg_l, ravg_l, rvw_l, &
                                       cext_l, csca_l, cbak_l, g_l, &
                                       gc_l, lc_l, pf_l, &
                                       save1, save2, &
                                       max_coef, verbose, n_threads, use_mpi) bind(c)

     use iso_c_binding

     implicit none

     integer(c_int), intent(in)  :: calc_gc
     integer(c_int), intent(in)  :: calc_lc
     integer(c_int), intent(in)  :: calc_pf
     integer(c_int), intent(in)  :: dist_type
     integer(c_int), intent(in)  :: n_int1
     integer(c_int), intent(in)  :: n_int2
     integer(c_int), intent(in)  :: n_quad
     integer(c_int), intent(in)  :: n_angles
     integer(c_int), intent(in)  :: n_derivs
     integer(c_int), intent(in)  :: save_control
     real(c_double), intent(in)  :: lambda
     real(c_double), intent(in)  :: mr
     real(c_double), intent(in)  :: mi
     real(c_double), intent(in)  :: a1
     real(c_double), intent(in)  :: a2
     real(c_double), intent(in)  :: a3
     real(c_double), intent(in)  :: a4
     real(c_double), intent(in)  :: a5
     real(c_double), intent(in)  :: r1
     real(c_double), intent(in)  :: r2
     real(c_double), intent(in)  :: lambda_l(*)
     real(c_double), intent(in)  :: mr_l(*)
     real(c_double), intent(in)  :: mi_l(*)
     real(c_double), intent(in)  :: a1_l(*)
     real(c_double), intent(in)  :: a2_l(*)
     real(c_double), intent(in)  :: a3_l(*)
     real(c_double), intent(in)  :: a4_l(*)
     real(c_double), intent(in)  :: a5_l(*)
     real(c_double), intent(in)  :: r1_l(*)
     real(c_double), intent(in)  :: r2_l(*)
     real(c_double), intent(in)  :: accuracy

     integer(c_int), intent(out) :: n_coef
     real(c_double), intent(out) :: r21
     real(c_double), intent(out) :: r22
     real(c_double), intent(out) :: norm
     real(c_double), intent(out) :: reff
     real(c_double), intent(out) :: veff
     real(c_double), intent(out) :: gavg
     real(c_double), intent(out) :: vavg
     real(c_double), intent(out) :: ravg
     real(c_double), intent(out) :: rvw
     real(c_double), intent(out) :: cext
     real(c_double), intent(out) :: csca
     real(c_double), intent(out) :: cbak
     real(c_double), intent(out) :: g
     real(c_double), intent(out) :: gc(*)
     real(c_double), intent(out) :: lc(*)
     real(c_double), intent(out) :: theta(*)
     real(c_double), intent(out) :: pf(*)
     real(c_double), intent(out) :: r21_l(*)
     real(c_double), intent(out) :: r22_l(*)
     real(c_double), intent(out) :: norm_l(*)
     real(c_double), intent(out) :: reff_l(*)
     real(c_double), intent(out) :: veff_l(*)
     real(c_double), intent(out) :: gavg_l(*)
     real(c_double), intent(out) :: vavg_l(*)
     real(c_double), intent(out) :: ravg_l(*)
     real(c_double), intent(out) :: rvw_l(*)
     real(c_double), intent(out) :: cext_l(*)
     real(c_double), intent(out) :: csca_l(*)
     real(c_double), intent(out) :: cbak_l(*)
     real(c_double), intent(out) :: g_l(*)
     real(c_double), intent(out) :: gc_l(*)
     real(c_double), intent(out) :: lc_l(*)
     real(c_double), intent(out) :: pf_l(*)

     type(c_ptr), intent(inout)  :: save1
     type(c_ptr), intent(inout)  :: save2

     integer(c_int), intent(in)  :: max_coef
     integer(c_int), intent(in)  :: verbose
     integer(c_int), intent(in)  :: n_threads
     integer(c_int), intent(in)  :: use_mpi

     end function lmie_solution2_l
end interface


contains


!*******************************************************************************
!
!*******************************************************************************
subroutine lmie_in_allocate_f90(in, n_derivs)

     implicit none

     type(lmie_in_type), intent(out) :: in
     integer,            intent(in)  :: n_derivs

     if (n_derivs .eq. 0) then
          in%lambda_l => zero1
          in%mr_l     => zero1
          in%mi_l     => zero1
          in%a1_l     => zero1
          in%a2_l     => zero1
          in%a3_l     => zero1
          in%a4_l     => zero1
          in%a5_l     => zero1
          in%r1_l     => zero1
          in%r2_l     => zero1
     else
          allocate(in%lambda_l(n_derivs))
          allocate(in%mr_l(n_derivs))
          allocate(in%mi_l(n_derivs))
          allocate(in%a1_l(n_derivs))
          allocate(in%a2_l(n_derivs))
          allocate(in%a3_l(n_derivs))
          allocate(in%a4_l(n_derivs))
          allocate(in%a5_l(n_derivs))
          allocate(in%r1_l(n_derivs))
          allocate(in%r2_l(n_derivs))
     endif

end subroutine lmie_in_allocate_f90



subroutine lmie_in_deallocate_f90(in, flag)

     implicit none

     type(lmie_in_type), intent(inout) :: in
     logical,            intent(in)    :: flag

     if (flag) then
          deallocate(in%lambda_l)
          deallocate(in%mr_l)
          deallocate(in%mi_l)
          deallocate(in%a1_l)
          deallocate(in%a2_l)
          deallocate(in%a3_l)
          deallocate(in%a4_l)
          deallocate(in%a5_l)
          deallocate(in%r1_l)
          deallocate(in%r2_l)
    endif

end subroutine lmie_in_deallocate_f90



!*******************************************************************************
!
!*******************************************************************************
subroutine lmie_out_allocate_f90(in, out, max_coef)

     implicit none

     type(lmie_in_type),  intent(in)  :: in
     type(lmie_out_type), intent(out) :: out
     integer,             intent(in)  :: max_coef

     if (in%n_derivs .eq. 0) then
        out%r1_l   => zero1
        out%r2_l   => zero1
        out%norm_l => zero1
        out%reff_l => zero1
        out%veff_l => zero1
        out%gavg_l => zero1
        out%vavg_l => zero1
        out%ravg_l => zero1
        out%rvw_l  => zero1
        out%cext_l => zero1
        out%csca_l => zero1
        out%cbak_l => zero1
        out%g_l    => zero1
     else
        allocate(out%r1_l  (in%n_derivs))
        allocate(out%r2_l  (in%n_derivs))
        allocate(out%norm_l(in%n_derivs))
        allocate(out%reff_l(in%n_derivs))
        allocate(out%veff_l(in%n_derivs))
        allocate(out%gavg_l(in%n_derivs))
        allocate(out%vavg_l(in%n_derivs))
        allocate(out%ravg_l(in%n_derivs))
        allocate(out%rvw_l (in%n_derivs))
        allocate(out%cext_l(in%n_derivs))
        allocate(out%csca_l(in%n_derivs))
        allocate(out%cbak_l(in%n_derivs))
        allocate(out%g_l   (in%n_derivs))
     endif

     out%gc   => zero2
     out%gc_l => zero3
     if (in%calc_gc) then
          allocate(out%gc(max_coef, 6))
          if (in%n_derivs > 0) then
               allocate(out%gc_l(max_coef, 6, in%n_derivs))
          endif
     endif

     out%lc   => zero2
     out%lc_l => zero3
     if (in%calc_lc) then
          allocate(out%lc(max_coef, 6))
          if (in%n_derivs > 0) then
               allocate(out%lc_l(max_coef, 6, in%n_derivs))
          endif
     endif

     out%theta => zero1
     out%pf    => zero2
     out%pf_l  => zero3
     if (in%calc_pf) then
          allocate(out%theta(in%n_angles))
          allocate(out%pf   (in%n_angles, 6))
          if (in%n_derivs > 0) then
               allocate(out%pf_l(in%n_angles, 6, in%n_derivs))
          endif
     endif

end subroutine lmie_out_allocate_f90



subroutine lmie_out_deallocate_f90(out, flag)

     implicit none

     type(lmie_out_type), intent(inout) :: out
     logical,             intent(in)    :: flag

     if (flag) then
        deallocate(out%r1_l)
        deallocate(out%r2_l)
        deallocate(out%norm_l)
        deallocate(out%reff_l)
        deallocate(out%veff_l)
        deallocate(out%gavg_l)
        deallocate(out%vavg_l)
        deallocate(out%ravg_l)
        deallocate(out%rvw_l)
        deallocate(out%cext_l)
        deallocate(out%csca_l)
        deallocate(out%cbak_l)
        deallocate(out%g_l)
     endif

     if (.not. associated(out%gc, zero2)) then
          deallocate(out%gc)
          if (flag) then
               deallocate(out%gc_l)
          endif
     endif

     if (.not. associated(out%lc, zero2)) then
          deallocate(out%lc)
          if (flag) then
               deallocate(out%lc_l)
          endif
     endif

     if (.not. associated(out%pf, zero2)) then
          deallocate(out%theta)
          deallocate(out%pf)
          if (flag) then
               deallocate(out%pf_l)
          endif
     endif

end subroutine lmie_out_deallocate_f90



!*******************************************************************************
!
!*******************************************************************************
subroutine lmie_in_zero_derivs_f90(in, n_derivs)

     implicit none

     type(lmie_in_type), intent(inout) :: in
     integer,            intent(in)    :: n_derivs

     integer :: i

     do i = 1, n_derivs
          in%lambda_l(i) = 0.d0
          in%mr_l(i)     = 0.d0
          in%mi_l(i)     = 0.d0
          in%a1_l(i)     = 0.d0
          in%a2_l(i)     = 0.d0
          in%a3_l(i)     = 0.d0
          in%a4_l(i)     = 0.d0
          in%a5_l(i)     = 0.d0
          in%r1_l(i)     = 0.d0
          in%r2_l(i)     = 0.d0
     enddo

end subroutine lmie_in_zero_derivs_f90



!*******************************************************************************
!
!*******************************************************************************
subroutine lmie_solution_f90(in, out, alloc, verbose, n_threads, use_mpi, error)

     implicit none

     type(lmie_in_type),  intent(in)  :: in
     type(lmie_out_type), intent(out) :: out
     logical,             intent(in)  :: alloc
     logical,             intent(in)  :: verbose
     integer,             intent(in)  :: n_threads
     logical,             intent(in)  :: use_mpi
     integer,             intent(out) :: error

     integer                          :: ret_val
     integer                          :: max_coef
     integer                          :: calc_gc
     integer                          :: calc_lc
     integer                          :: calc_pf
     integer                          :: verbose_
     integer                          :: use_mpi_

     integer                          :: lmie_calc_max_coef
     integer                          :: lmie_solution2_l

     if (alloc) then
          max_coef = lmie_calc_max_coef(in%lambda, in%dist_type, in%a1, in%a2, &
                                        in%a3, in%r1, in%r2)

          call lmie_out_allocate_f90(in, out, max_coef)
     endif

     calc_gc = 0
     calc_lc = 0
     calc_pf = 0
     verbose_ = 0
     use_mpi_ = 0

     if (in%calc_gc) calc_gc = 1
     if (in%calc_lc) calc_lc = 1
     if (in%calc_pf) calc_pf = 1
     if (verbose) verbose_ = 1
     if (use_mpi) use_mpi_ = 1

     ret_val = lmie_solution2_l(calc_gc, calc_lc, calc_pf, &
                                in%dist_type, in%n_int1, in%n_int2, in%n_quad, &
                                in%n_angles, in%n_derivs, in%save_control, &
                                in%lambda, in%mr, in%mi, &
                                in%a1, in%a2, in%a3, in%a4, in%a5, &
                                in%r1, in%r2, &
                                in%lambda_l, in%mr_l, in%mi_l, &
                                in%a1_l, in%a2_l, in%a3_l, in%a4_l, in%a5_l, &
                                in%r1_l, in%r2_l, &
                                in%accuracy, out%n_coef, &
                                out%r1, out%r2, &
                                out%norm, out%reff, out%veff, &
                                out%gavg, out%vavg, out%ravg, out%rvw, &
                                out%cext, out%csca, out%cbak, out%g, &
                                out%gc, out%lc, out%theta, out%pf, &
                                out%r1_l, out%r2_l, &
                                out%norm_l, out%reff_l, out%veff_l, &
                                out%gavg_l, out%vavg_l, out%ravg_l, out%rvw_l, &
                                out%cext_l, out%csca_l, out%cbak_l, out%g_l, &
                                out%gc_l, out%lc_l, out%pf_l, &
                                out%save1, out%save2, &
                                max_coef, verbose_, n_threads, use_mpi_)
     if (ret_val /= 0) then
          write(0, *) 'lmie_solution2_l()'
          error = -1;
          return
     endif

     error = 0;

end subroutine lmie_solution_f90

end module lmie_int_f90
