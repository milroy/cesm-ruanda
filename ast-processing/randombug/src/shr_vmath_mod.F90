!===============================================================================
! SVN $Id: shr_vmath_mod.F90 6752 2007-10-04 21:02:15Z jwolfe $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140626/shr/shr_vmath_mod.F90 $
!===============================================================================
! PURPOSE: 
!   provides a uniform, platform-independent API for vector math functions
!===============================================================================

module shr_vmath_mod

   !----------------------------------------------------------------------------
   ! routines that evaluate various math functions for vector arguments
   ! intended to provide platform independent access to vendor optimized code
   !----------------------------------------------------------------------------

   use shr_kind_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none

   private
   public :: shr_vmath_sqrt, &
      shr_vmath_exp, shr_vmath_log, &
      shr_vmath_sin, shr_vmath_cos, &
      shr_vmath_rsqrt, shr_vmath_div

   contains

!===============================================================================

! Subprogram not used subroutine shr_vmath_sqrt(X, Y, n)
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: sqrt for vector arguments, optimized on different platforms
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    Y = sqrt(X)
! Subprogram not used 
! Subprogram not used end subroutine shr_vmath_sqrt

!===============================================================================

! Subprogram not used subroutine shr_vmath_rsqrt(X, Y, n)
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: sqrt for vector arguments, optimized on different platforms
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    Y = 1.0_SHR_KIND_R8/sqrt(X)
! Subprogram not used 
! Subprogram not used end subroutine shr_vmath_rsqrt

!===============================================================================

! Subprogram not used subroutine shr_vmath_exp(X, Y, n)
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: exp for vector arguments, optimized on different platforms
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    Y = exp(X)
! Subprogram not used 
! Subprogram not used end subroutine shr_vmath_exp

!===============================================================================

! Subprogram not used subroutine shr_vmath_div(X, Y, Z, n)
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: Y(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(out) :: Z(n) ! output vector argument
! Subprogram not used 
! Subprogram not used    integer :: i
! Subprogram not used    do i=1,n
! Subprogram not used       Z(i) = X(i)/Y(i)
! Subprogram not used    enddo
! Subprogram not used    return
! Subprogram not used  end subroutine shr_vmath_div

!===============================================================================

subroutine shr_vmath_log(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: log for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------





   Y = log(X)

end subroutine shr_vmath_log

!===============================================================================

! Subprogram not used subroutine shr_vmath_sin(X, Y, n)
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: sin for vector arguments, optimized on different platforms
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    Y = sin(X)
! Subprogram not used 
! Subprogram not used end subroutine shr_vmath_sin

!===============================================================================

! Subprogram not used subroutine shr_vmath_cos(X, Y, n)
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
! Subprogram not used    real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
! Subprogram not used    real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: cos for vector arguments, optimized on different platforms
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    Y = cos(X)
! Subprogram not used 
! Subprogram not used end subroutine shr_vmath_cos

!===============================================================================

end module shr_vmath_mod
