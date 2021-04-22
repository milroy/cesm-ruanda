
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)

      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,3) = 2.9e-12_r8 * exp( -160._r8 * itemp(:,:) )
      rate(:,:,5) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,7) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )

      end subroutine setrxt


! Subprogram not used       subroutine setrxt_hrates( rate, temp, m, ncol, kbot )
! Subprogram not used 
! Subprogram not used       use ppgrid,       only : pver, pcols
! Subprogram not used       use shr_kind_mod, only : r8 => shr_kind_r8
! Subprogram not used       use chem_mods, only : rxntot
! Subprogram not used       use mo_jpl,    only : jpl
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------
! Subprogram not used !       ... dummy arguments
! Subprogram not used !-------------------------------------------------------
! Subprogram not used       integer, intent(in) :: ncol
! Subprogram not used       integer, intent(in) :: kbot
! Subprogram not used       real(r8), intent(in)    :: temp(pcols,pver)
! Subprogram not used       real(r8), intent(in)    :: m(ncol,pver)
! Subprogram not used       real(r8), intent(inout) :: rate(ncol,pver,rxntot)
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------
! Subprogram not used !       ... local variables
! Subprogram not used !-------------------------------------------------------
! Subprogram not used       integer  ::  n
! Subprogram not used       real(r8)  ::  itemp(ncol,kbot)
! Subprogram not used       real(r8)  ::  exp_fac(ncol,kbot)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       end subroutine setrxt_hrates

      end module mo_setrxt
