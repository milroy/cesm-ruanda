
      module mo_airglow

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use abortutils,    only : endrun

      implicit none

      save

      integer , parameter :: nag      = 2
      real(r8), parameter :: secpday  = 86400._r8
      real(r8), parameter :: daypsec  = 1._r8/secpday
      real(r8), parameter :: hc       = 6.62608e-34_r8*2.9979e8_r8/1.e-9_r8
      real(r8), parameter :: wc_o2_1s = 1._r8/762._r8
      real(r8), parameter :: wc_o2_1d = 1._r8/1270._r8

      integer :: rid_ag1, rid_ag2
      logical :: has_airglow

      private
      public :: airglow, init_airglow

      contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine init_airglow
          use mo_chem_utls, only : get_rxt_ndx
          use cam_history,  only : addfld, phys_decomp
          use ppgrid,       only : pver

          implicit none

          rid_ag1 = get_rxt_ndx( 'ag1' )
          rid_ag2 = get_rxt_ndx( 'ag2' )

          has_airglow = rid_ag1 > 0 .and. rid_ag2 > 0 

          if (.not. has_airglow) return

          call addfld( 'AG1',   'K/s ', pver, 'I', 'O2_1S -> O2 + 762nm airglow loss', phys_decomp )
          call addfld( 'AG2',   'K/s ', pver, 'I', 'O2_1D -> O2 + 1.27 micron airglow loss', phys_decomp )
          call addfld( 'AGTOT', 'K/s ', pver, 'I', 'airglow total loss', phys_decomp )

        endsubroutine init_airglow

! Subprogram not used       subroutine airglow( ag_tot, o2_1s, o2_1d, rxt, cp, &
! Subprogram not used                           ncol, lchnk )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !      	... forms the airglow heating rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use chem_mods,     only : rxntot
! Subprogram not used       use ppgrid,        only : pver
! Subprogram not used       use cam_history,   only : outfld
! Subprogram not used       use mo_constants,  only : avo => avogadro
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   ::  ncol                                ! columns in chunck
! Subprogram not used       integer, intent(in)   ::  lchnk                               ! chunk index
! Subprogram not used       real(r8), intent(in)  ::  rxt(ncol,pver,rxntot)               ! rxt rates (1/cm^3/s)
! Subprogram not used       real(r8), intent(in)  ::  o2_1s(ncol,pver)                    ! concentration (mol/mol)
! Subprogram not used       real(r8), intent(in)  ::  o2_1d(ncol,pver)                    ! concentration (mol/mol)
! Subprogram not used       real(r8), intent(in)  ::  cp(ncol,pver)                       ! specific heat capacity
! Subprogram not used       real(r8), intent(out) ::  ag_tot(ncol,pver)                   ! airglow total heating rate (K/s)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer  ::  k
! Subprogram not used       real(r8) ::  tmp(ncol)
! Subprogram not used       real(r8) ::  ag_rate(ncol,pver,nag)
! Subprogram not used 
! Subprogram not used       if (.not. has_airglow) return
! Subprogram not used 
! Subprogram not used       do k = 1,pver
! Subprogram not used          tmp(:)          = hc * avo / cp(:,k)
! Subprogram not used          ag_rate(:,k,1)  = tmp(:)*rxt(:,k,rid_ag1)*o2_1d(:,k)*wc_o2_1d
! Subprogram not used          ag_rate(:,k,2)  = tmp(:)*rxt(:,k,rid_ag2)*o2_1s(:,k)*wc_o2_1s
! Subprogram not used          ag_tot(:,k)     = ag_rate(:,k,1) + ag_rate(:,k,2)
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... output the rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call outfld( 'AG1', ag_rate(:,:,1), ncol, lchnk )
! Subprogram not used       call outfld( 'AG2', ag_rate(:,:,2), ncol, lchnk )
! Subprogram not used       call outfld( 'AGTOT', ag_tot, ncol, lchnk )
! Subprogram not used 
! Subprogram not used       end subroutine airglow

      end module mo_airglow
