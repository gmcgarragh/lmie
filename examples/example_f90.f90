!*******************************************************************************
! This example program calls LMie to calculate the scattering properties of the
! accumulation mode mineral dust described by d'Almeida, 1991 at 2.25um.  A log
! normal size distribution is used with a mean radius of 0.39um and a standard
! deviation of 2.0um ranging from 0.005 to 50um.  The real and imaginary parts
! of the index of refraction are 1.53 and 0.0054, respectively.  In this case
! the output includes size distribution statistics, extinction, scattering, and
! backscattering cross sections, asymmetry parameter, coefficients for the
! expansion of the scattering matrix in generalized spherical functions, and the
! elements of the normalized scattering matrix at 181 angles.  In addition,
! derivatives of the these output quantities with respect to mean radius,
! standard deviation, and the real and imaginary parts of the index of
! refraction are also generated.
!*******************************************************************************
program example_f90

     ! Use the LMie Fortran90 interface module
     use lmie_int_f90

     implicit none

     integer             :: i
     integer             :: error

     integer, parameter  :: n_derivs = 4

     type(lmie_in_type)  :: in          ! LMie input type
     type(lmie_out_type) :: out         ! LMie output type


     !**************************************************************************
     ! Fill the LMie input structure.  See the LMie documentation for details of
     ! each input quantity.
     !**************************************************************************

     ! Allocate memory used in the input type for n_derivs derivatives.
     call lmie_in_allocate_f90(in, n_derivs)

     in%save_control = 0

     in%calc_gc   = .true.
     in%calc_lc   = .false.
     in%calc_pf   = .true.

     in%dist_type = DIST_LOG_NORMAL

!    in%n_int1    = 			! not used for log normal
     in%n_int2    = 64
     in%n_quad    = 64
     in%n_angles  = 181
     in%n_derivs  = n_derivs

     in%lambda    = 2.250d+00
     in%mr        = 1.530d+00
     in%mi        = 5.400d-03
     in%a1        = 3.900d-01
     in%a2        = 4.805d-01
!    in%a3        = 			! not used for log normal
!    in%a4        = 			! not used for log normal
!    in%a5        = 			! not used for log normal
     in%r1        = 5.000d-03
     in%r2        = 5.000d+01

     in%accuracy  = 1.d-7


     !**************************************************************************
     ! Set linearized inputs.  The strategy here is to set all linearized inputs
     ! to zero then set the appropriate linearized inputs to unity.
     !**************************************************************************
     do i = 1, n_derivs
          in%lambda_l(i) = 0.d0
          in%mr_l(i)     = 0.d0
          in%mi_l(i)     = 0.d0
          in%a1_l(i)     = 0.d0
          in%a2_l(i)     = 0.d0
!         in%a3_l(i)     = 0.d0		! not used for log normal
!         in%a4_l(i)     = 0.d0		! not used for log normal
!         in%a5_l(i)     = 0.d0		! not used for log normal
          in%r1_l(i)     = 0.d0
          in%r2_l(i)     = 0.d0
     enddo

     ! Or more conveniently just call the function.
     call lmie_in_zero_derivs_f90(in, in%n_derivs)

     in%mr_l(1)   = 1.d0                ! derivative 0 is with respect to mr
     in%mi_l(2)   = 1.d0                ! derivative 1 is with respect to mi
     in%a1_l(3)   = 1.d0                ! derivative 2 is with respect to a1
     in%a2_l(4)   = 1.d0                ! derivative 3 is with respect to a2



     !**************************************************************************
     ! Call LMie to get the Mie solution.  In this case we indicate that memory
     ! used in the output structure should be allocated, verbose output is off,
     ! two threads are used if LMie was compiled with multi-threading support,
     ! and MPI should not be used.
     !**************************************************************************
     call lmie_solution_f90(in, out, .true., .false., 2, .false., error)
     if (error /= 0) then
          write(0, *) 'lmie_solution_f90()'
          error = 1;
          stop
     endif


     !**************************************************************************
     ! Output results.
     !**************************************************************************
     write(*, '("---------------------------------------------------------------------------")')
     write(*, '("full quantities (non linearized)")')
     write(*, '("---------------------------------------------------------------------------")')

     call output_result(out%reff, out%veff, &
                        out%gavg, out%vavg, out%ravg, out%rvw, &
                        out%cext, out%csca, out%cbak, out%g, &
                        out%n_coef, size(out%gc, 1), out%gc, &
                        in%n_angles, out%theta, out%pf)

     write(*, '()')

     do i = 1, n_derivs
          write(*, '("---------------------------------------------------------------------------")')
          if (i .eq. 1) then
               write(*, '("derivatives wrt the real part of the index of refraction (mr)")')
          else if (i .eq. 2) then
               write(*, '("derivatives wrt the imaginary part of the index of refraction (mi)")')
          else if (i .eq. 3) then
               write(*, '("derivatives wrt log normal size distribution mean radius (a1)")')
          else if (i .eq. 4) then
               write(*, '("derivatives wrt log normal size distribution standard deviation (a2)")')
          endif
          write(*, '("---------------------------------------------------------------------------")')

          call output_result(out%reff_l(i), out%veff_l(i), &
                             out%gavg_l(i), out%vavg_l(i), out%ravg_l(i), out%rvw_l(i), &
                             out%cext_l(i), out%csca_l(i), out%cbak_l(i), out%g_l(i), &
                             out%n_coef, size(out%gc_l, 1), out%gc_l(:,:,i), &
                             in%n_angles, out%theta, out%pf_l(:,:,i))

          write(*, '()')
     enddo


     !**************************************************************************
     ! Free memory allocated in the input and output types.
     !**************************************************************************
     call lmie_in_deallocate_f90 (in,  in%n_derivs > 0)
     call lmie_out_deallocate_f90(out, in%n_derivs > 0)


end program example_f90



!*******************************************************************************
! Subroutine to output results.
!*******************************************************************************
subroutine output_result(reff, veff, &
                         gavg, vavg, ravg, rvw, &
                         cext, csca, cbak, g, &
                         n_coef, max_coef, gc, &
                         n_angles, theta, pf)

     implicit none

     real(8), intent(in) :: reff
     real(8), intent(in) :: veff
     real(8), intent(in) :: gavg
     real(8), intent(in) :: vavg
     real(8), intent(in) :: ravg
     real(8), intent(in) :: rvw
     real(8), intent(in) :: cext
     real(8), intent(in) :: csca
     real(8), intent(in) :: cbak
     real(8), intent(in) :: g
     integer, intent(in) :: n_coef
     integer, intent(in) :: max_coef
     real(8), intent(in) :: gc(max_coef, 6)
     integer, intent(in) :: n_angles
     real(8), intent(in) :: theta(n_angles)
     real(8), intent(in) :: pf(n_angles, 6)

     integer :: i

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
     write(*, '("  alpha1        alpha2        alpha3        alpha4       -beta1        -beta2")')
     write(*, '("  beta          alpha         zeta          delta         gamma         epsilon")')
     write(*, '("------------------------------------------------------------------------------------")')
     do i = 1, n_coef
          write(*, '(6ES14.6)') gc(i, 1), gc(i, 2), gc(i, 3), gc(i, 4), gc(i, 5), gc(i, 6)
     enddo
     write(*, '()')

     write(*, '("elements of the normalized scattering matrix:")')
     write(*, '("  theta         F11           F22           F33           F44           F12           F34")')
     write(*, '("                a1            a2            a3            a4            b1            b2")')
     write(*, '("------------------------------------------------------------------------------------------------")')
     do i = 1, n_angles
          write(*, '(7ES14.6)') theta(i), pf(i, 1), pf(i, 2), pf(i, 3), pf(i, 4), pf(i, 5), pf(i, 6)
     enddo

end subroutine output_result
