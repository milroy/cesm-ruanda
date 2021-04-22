module ebert_curry

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag, get_sw_spectral_boundaries
use abortutils,       only: endrun
use cam_history,      only: outfld

implicit none
private
save

public :: &
   ec_rad_props_init,        &
   cloud_rad_props_get_sw, & ! return SW optical props of total bulk aerosols
   cloud_rad_props_get_lw, & ! return LW optical props of total bulk aerosols
   ec_ice_optics_sw,       &
   ec_ice_get_rad_props_lw


real, public, parameter:: scalefactor = 1._r8 !500._r8/917._r8

! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
! 
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
!
! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
! 
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)

! 
! indexes into pbuf for optical parameters of MG clouds
! 
   integer :: dei_idx     = 0 
   integer :: mu_idx      = 0
   integer :: lambda_idx  = 0 
   integer :: iciwp_idx   = 0 
   integer :: iclwp_idx   = 0 
   integer :: cld_idx     = 0  
   integer :: rei_idx     = 0  

! indexes into constituents for old optics
   integer :: &
        ixcldice,           & ! cloud ice water index
        ixcldliq              ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine ec_rad_props_init()

   use cam_history, only: addfld, phys_decomp
   use netcdf
   use spmd_utils,     only: masterproc
   use ioFileMod,      only: getfil
   use cam_logfile,    only: iulog
   use error_messages, only: handle_ncerr

   use mpishorthand

   use constituents,   only: cnst_get_ind

   integer :: err

   iciwp_idx  = pbuf_get_index('ICIWP',errcode=err)
   iclwp_idx  = pbuf_get_index('ICLWP',errcode=err)
   cld_idx    = pbuf_get_index('CLD')
   rei_idx    = pbuf_get_index('REI')

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   !call addfld ('CLWPTH_OLD','Kg/m2   ',pver, 'I','old In Cloud Liquid Water Path',phys_decomp, sampling_seq='rad_lwsw')
   !call addfld ('KEXT_OLD','m^2/kg',pver,'I','old extinction',phys_decomp)
   !call addfld ('CLDOD_OLD','1',pver,'I','old liquid OD',phys_decomp)
   !call addfld ('REL_OLD','1',pver,'I','old liquid effective radius (liquid)',phys_decomp)

   !call addfld ('CLWPTH_NEW','Kg/m2   ',pver, 'I','In Cloud Liquid Water Path',phys_decomp, sampling_seq='rad_lwsw')
   !call addfld ('KEXT_NEW','m^2/kg',pver,'I','extinction',phys_decomp)
   !call addfld ('CLDOD_NEW','1',pver,'I','liquid OD',phys_decomp)

   !call addfld('CIWPTH_NEW','Kg/m2   ',pver, 'I','In Cloud Ice Water Path',phys_decomp, sampling_seq='rad_lwsw')
   !call addfld('CIWPTH_OLD','Kg/m2   ',pver, 'I','In Cloud Ice Water Path (old)',phys_decomp, sampling_seq='rad_lwsw')

    return

end subroutine ec_rad_props_init

!==============================================================================

! Subprogram not used subroutine cloud_rad_props_get_sw(state, pbuf, &
! Subprogram not used                                   tau, tau_w, tau_w_g, tau_w_f,&
! Subprogram not used                                   diagnosticindex, oldliq, oldice)
! Subprogram not used 
! Subprogram not used ! return totaled (across all species) layer tau, omega, g, f 
! Subprogram not used ! for all spectral interval for aerosols affecting the climate
! Subprogram not used 
! Subprogram not used    ! Arguments
! Subprogram not used    type(physics_state), intent(in) :: state
! Subprogram not used    type(physics_buffer_desc),  pointer :: pbuf(:)
! Subprogram not used    integer, optional,   intent(in) :: diagnosticindex      ! index (if present) to radiation diagnostic information
! Subprogram not used 
! Subprogram not used    real(r8), intent(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
! Subprogram not used    real(r8), intent(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
! Subprogram not used    real(r8), intent(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
! Subprogram not used    real(r8), intent(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w
! Subprogram not used 
! Subprogram not used    logical, optional, intent(in) :: oldliq,oldice
! Subprogram not used 
! Subprogram not used    ! Local variables
! Subprogram not used 
! Subprogram not used    integer :: ncol
! Subprogram not used    integer :: lchnk
! Subprogram not used    integer :: k, i    ! lev and daycolumn indices
! Subprogram not used    integer :: iswband ! sw band indices
! Subprogram not used 
! Subprogram not used    !-----------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ncol  = state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    ! initialize to conditions that would cause failure
! Subprogram not used    tau     (:,:,:) = -100._r8
! Subprogram not used    tau_w   (:,:,:) = -100._r8
! Subprogram not used    tau_w_g (:,:,:) = -100._r8
! Subprogram not used    tau_w_f (:,:,:) = -100._r8
! Subprogram not used 
! Subprogram not used    ! initialize layers to accumulate od's
! Subprogram not used    tau    (:,1:ncol,:) = 0._r8
! Subprogram not used    tau_w  (:,1:ncol,:) = 0._r8
! Subprogram not used    tau_w_g(:,1:ncol,:) = 0._r8
! Subprogram not used    tau_w_f(:,1:ncol,:) = 0._r8
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call ec_ice_optics_sw   (state, pbuf, tau, tau_w, tau_w_g, tau_w_f, oldicewp=.true.)
! Subprogram not used !  call outfld ('CI_OD_SW_OLD', ice_tau(idx_sw_diag,:,:), pcols, lchnk)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used end subroutine cloud_rad_props_get_sw
!==============================================================================

! Subprogram not used subroutine cloud_rad_props_get_lw(state, pbuf, cld_abs_od, diagnosticindex, oldliq, oldice, oldcloud)
! Subprogram not used 
! Subprogram not used ! Purpose: Compute cloud longwave absorption optical depth
! Subprogram not used !    cloud_rad_props_get_lw() is called by radlw() 
! Subprogram not used 
! Subprogram not used    ! Arguments
! Subprogram not used    type(physics_state), intent(in)  :: state
! Subprogram not used    type(physics_buffer_desc), pointer  :: pbuf(:)
! Subprogram not used    real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
! Subprogram not used    integer, optional,   intent(in)  :: diagnosticindex
! Subprogram not used    logical, optional,   intent(in)  :: oldliq  ! use old liquid optics
! Subprogram not used    logical, optional,   intent(in)  :: oldice  ! use old ice optics
! Subprogram not used    logical, optional,   intent(in)  :: oldcloud  ! use old optics for both (b4b)
! Subprogram not used 
! Subprogram not used    ! Local variables
! Subprogram not used 
! Subprogram not used    integer :: bnd_idx     ! LW band index
! Subprogram not used    integer :: i           ! column index
! Subprogram not used    integer :: k           ! lev index
! Subprogram not used    integer :: ncol        ! number of columns
! Subprogram not used    integer :: lchnk
! Subprogram not used 
! Subprogram not used    !-----------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ncol = state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    ! compute optical depths cld_absod 
! Subprogram not used    cld_abs_od = 0._r8
! Subprogram not used 
! Subprogram not used    call ec_ice_get_rad_props_lw(state, pbuf, cld_abs_od, oldicewp=.true.)
! Subprogram not used    !call outfld('CI_OD_LW_OLD', ice_tau_abs_od(idx_lw_diag ,:,:), pcols, lchnk)
! Subprogram not used       
! Subprogram not used end subroutine cloud_rad_props_get_lw

!==============================================================================
! Private methods
!==============================================================================

! Subprogram not used subroutine ec_ice_optics_sw   (state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp)
! Subprogram not used 
! Subprogram not used    use physconst,      only: gravit
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in) :: state
! Subprogram not used    type(physics_buffer_desc), pointer :: pbuf(:)
! Subprogram not used 
! Subprogram not used    real(r8),intent(out) :: ice_tau    (nswbands,pcols,pver) ! extinction optical depth
! Subprogram not used    real(r8),intent(out) :: ice_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
! Subprogram not used    real(r8),intent(out) :: ice_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
! Subprogram not used    real(r8),intent(out) :: ice_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w
! Subprogram not used    logical, intent(in) :: oldicewp
! Subprogram not used 
! Subprogram not used    real(r8), pointer, dimension(:,:) :: rei
! Subprogram not used    real(r8), pointer, dimension(:,:) :: cldn
! Subprogram not used    real(r8), pointer, dimension(:,:) :: tmpptr
! Subprogram not used    real(r8), dimension(pcols,pver) :: cicewp
! Subprogram not used    real(r8), dimension(nswbands) :: wavmin
! Subprogram not used    real(r8), dimension(nswbands) :: wavmax
! Subprogram not used    !
! Subprogram not used    ! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
! Subprogram not used    real(r8) :: abari(4) = &     ! a coefficient for extinction optical depth
! Subprogram not used       (/ 3.448e-03_r8, 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8/)
! Subprogram not used    real(r8) :: bbari(4) = &     ! b coefficient for extinction optical depth
! Subprogram not used       (/ 2.431_r8    , 2.431_r8    ,2.431_r8    ,2.431_r8    /)
! Subprogram not used    real(r8) :: cbari(4) = &     ! c coefficient for single scat albedo
! Subprogram not used       (/ 1.00e-05_r8 , 1.10e-04_r8 ,1.861e-02_r8,.46658_r8   /)
! Subprogram not used    real(r8) :: dbari(4) = &     ! d coefficient for single scat albedo
! Subprogram not used       (/ 0.0_r8      , 1.405e-05_r8,8.328e-04_r8,2.05e-05_r8 /)
! Subprogram not used    real(r8) :: ebari(4) = &     ! e coefficient for asymmetry parameter
! Subprogram not used       (/ 0.7661_r8   , 0.7730_r8   ,0.794_r8    ,0.9595_r8   /)
! Subprogram not used    real(r8) :: fbari(4) = &     ! f coefficient for asymmetry parameter
! Subprogram not used       (/ 5.851e-04_r8, 5.665e-04_r8,7.267e-04_r8,1.076e-04_r8/)
! Subprogram not used 
! Subprogram not used    real(r8) :: abarii           ! A coefficient for current spectral band
! Subprogram not used    real(r8) :: bbarii           ! B coefficient for current spectral band
! Subprogram not used    real(r8) :: cbarii           ! C coefficient for current spectral band
! Subprogram not used    real(r8) :: dbarii           ! D coefficient for current spectral band
! Subprogram not used    real(r8) :: ebarii           ! E coefficient for current spectral band
! Subprogram not used    real(r8) :: fbarii           ! F coefficient for current spectral band
! Subprogram not used 
! Subprogram not used    ! Minimum cloud amount (as a fraction of the grid-box area) to 
! Subprogram not used    ! distinguish from clear sky
! Subprogram not used    real(r8), parameter :: cldmin = 1.0e-80_r8
! Subprogram not used 
! Subprogram not used    ! Decimal precision of cloud amount (0 -> preserve full resolution;
! Subprogram not used    ! 10^-n -> preserve n digits of cloud amount)
! Subprogram not used    real(r8), parameter :: cldeps = 0.0_r8
! Subprogram not used 
! Subprogram not used    integer :: ns, i, k, indxsl, lchnk, Nday
! Subprogram not used    integer :: itim_old
! Subprogram not used    real(r8) :: tmp1i, tmp2i, tmp3i, g
! Subprogram not used 
! Subprogram not used    Nday = state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    itim_old = pbuf_old_tim_idx()
! Subprogram not used    call pbuf_get_field(pbuf, cld_idx,cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
! Subprogram not used    call pbuf_get_field(pbuf, rei_idx,rei)
! Subprogram not used 
! Subprogram not used    if(oldicewp) then
! Subprogram not used      do k=1,pver
! Subprogram not used         do i = 1,Nday
! Subprogram not used            cicewp(i,k) = 1000.0_r8*state%q(i,k,ixcldice)*state%pdel(i,k) /(gravit* max(0.01_r8,cldn(i,k)))
! Subprogram not used         end do
! Subprogram not used      end do
! Subprogram not used    else
! Subprogram not used      if (iciwp_idx<=0) then 
! Subprogram not used         call endrun('ec_ice_optics_sw: oldicewp must be set to true since ICIWP was not found in pbuf')
! Subprogram not used      endif
! Subprogram not used      call pbuf_get_field(pbuf, iciwp_idx, tmpptr)
! Subprogram not used      cicewp(1:pcols,1:pver) =  1000.0_r8*tmpptr(1:pcols,1:pver)
! Subprogram not used    endif
! Subprogram not used    
! Subprogram not used    call get_sw_spectral_boundaries(wavmin,wavmax,'microns')
! Subprogram not used 
! Subprogram not used    do ns = 1, nswbands
! Subprogram not used 
! Subprogram not used       if(wavmax(ns) <= 0.7_r8) then
! Subprogram not used          indxsl = 1
! Subprogram not used       else if(wavmax(ns) <= 1.25_r8) then
! Subprogram not used          indxsl = 2
! Subprogram not used       else if(wavmax(ns) <= 2.38_r8) then
! Subprogram not used          indxsl = 3
! Subprogram not used       else if(wavmax(ns) > 2.38_r8) then
! Subprogram not used          indxsl = 4
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       abarii = abari(indxsl)
! Subprogram not used       bbarii = bbari(indxsl)
! Subprogram not used       cbarii = cbari(indxsl)
! Subprogram not used       dbarii = dbari(indxsl)
! Subprogram not used       ebarii = ebari(indxsl)
! Subprogram not used       fbarii = fbari(indxsl)
! Subprogram not used 
! Subprogram not used       do k=1,pver
! Subprogram not used          do i=1,Nday
! Subprogram not used 
! Subprogram not used             ! note that optical properties for ice valid only
! Subprogram not used             ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
! Subprogram not used             if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
! Subprogram not used                tmp1i = abarii + bbarii/max(13._r8,min(scalefactor*rei(i,k),130._r8))
! Subprogram not used                ice_tau(ns,i,k) = cicewp(i,k)*tmp1i
! Subprogram not used             else
! Subprogram not used                ice_tau(ns,i,k) = 0.0_r8
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             tmp2i = 1._r8 - cbarii - dbarii*min(max(13._r8,scalefactor*rei(i,k)),130._r8)
! Subprogram not used             tmp3i = fbarii*min(max(13._r8,scalefactor*rei(i,k)),130._r8)
! Subprogram not used             ! Do not let single scatter albedo be 1.  Delta-eddington solution
! Subprogram not used             ! for non-conservative case has different analytic form from solution
! Subprogram not used             ! for conservative case, and raddedmx is written for non-conservative case.
! Subprogram not used             ice_tau_w(ns,i,k) = ice_tau(ns,i,k) * min(tmp2i,.999999_r8)
! Subprogram not used             g = ebarii + tmp3i
! Subprogram not used             ice_tau_w_g(ns,i,k) = ice_tau_w(ns,i,k) * g
! Subprogram not used             ice_tau_w_f(ns,i,k) = ice_tau_w(ns,i,k) * g * g
! Subprogram not used 
! Subprogram not used          end do ! End do i=1,Nday
! Subprogram not used       end do    ! End do k=1,pver
! Subprogram not used    end do ! nswbands
! Subprogram not used 
! Subprogram not used end subroutine ec_ice_optics_sw
!==============================================================================

! Subprogram not used subroutine ec_ice_get_rad_props_lw(state, pbuf, abs_od, oldicewp)
! Subprogram not used    use physconst,      only: gravit
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in)   :: state
! Subprogram not used    type(physics_buffer_desc),pointer :: pbuf(:)
! Subprogram not used     real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)
! Subprogram not used    logical, intent(in) :: oldicewp
! Subprogram not used 
! Subprogram not used    real(r8) :: gicewp(pcols,pver)
! Subprogram not used    real(r8) :: gliqwp(pcols,pver)
! Subprogram not used    real(r8) :: cicewp(pcols,pver)
! Subprogram not used    real(r8) :: cliqwp(pcols,pver)
! Subprogram not used    real(r8) :: ficemr(pcols,pver)
! Subprogram not used    real(r8) :: cwp(pcols,pver)
! Subprogram not used    real(r8) :: cldtau(pcols,pver)
! Subprogram not used 
! Subprogram not used    real(r8), pointer, dimension(:,:) :: cldn
! Subprogram not used    real(r8), pointer, dimension(:,:) :: rei
! Subprogram not used    integer :: ncol, itim_old, lwband, i, k, lchnk
! Subprogram not used 
! Subprogram not used     real(r8) :: kabs, kabsi
! Subprogram not used 
! Subprogram not used     real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
! Subprogram not used     parameter (kabsl = 0.090361_r8)
! Subprogram not used 
! Subprogram not used    real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    ncol = state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    itim_old  =  pbuf_old_tim_idx()
! Subprogram not used    call pbuf_get_field(pbuf, rei_idx,   rei)
! Subprogram not used    call pbuf_get_field(pbuf, cld_idx,   cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    if(oldicewp) then
! Subprogram not used      do k=1,pver
! Subprogram not used          do i = 1,ncol
! Subprogram not used             gicewp(i,k) = state%q(i,k,ixcldice)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
! Subprogram not used             gliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
! Subprogram not used             cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
! Subprogram not used             cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
! Subprogram not used             ficemr(i,k) = state%q(i,k,ixcldice) /                 &
! Subprogram not used                  max(1.e-10_r8,(state%q(i,k,ixcldice)+state%q(i,k,ixcldliq)))
! Subprogram not used          end do
! Subprogram not used      end do
! Subprogram not used      cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
! Subprogram not used    else
! Subprogram not used       if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
! Subprogram not used          call endrun('ec_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
! Subprogram not used       endif
! Subprogram not used       call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
! Subprogram not used       call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
! Subprogram not used       do k=1,pver
! Subprogram not used          do i = 1,ncol
! Subprogram not used             cwp(i,k) = 1000.0_r8 *iciwpth(i,k) + 1000.0_r8 *iclwpth(i,k)
! Subprogram not used             ficemr(i,k) = 1000.0_r8*iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    do k=1,pver
! Subprogram not used        do i=1,ncol
! Subprogram not used 
! Subprogram not used           ! Note from Andrew Conley:
! Subprogram not used           !  Optics for RK no longer supported, This is constructed to get
! Subprogram not used           !  close to bit for bit.  Otherwise we could simply use ice water path
! Subprogram not used           !note that optical properties for ice valid only
! Subprogram not used           !in range of 13 > rei > 130 micron (Ebert and Curry 92)
! Subprogram not used           kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
! Subprogram not used           kabs =  kabsi*ficemr(i,k) ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
! Subprogram not used           !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
! Subprogram not used           cldtau(i,k) = kabs*cwp(i,k)
! Subprogram not used        end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do lwband = 1,nlwbands
! Subprogram not used       abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    !if(oldicewp) then
! Subprogram not used    !  call outfld('CIWPTH_OLD',cicewp(:,:)/1000,pcols,lchnk)
! Subprogram not used    !else
! Subprogram not used    !  call outfld('CIWPTH_OLD',iciwpth(:,:),pcols,lchnk)
! Subprogram not used    !endif
! Subprogram not used    !call outfld('CI_OD_LW_OLD',cldtau(:,:),pcols,lchnk)
! Subprogram not used 
! Subprogram not used end subroutine ec_ice_get_rad_props_lw
!==============================================================================

end module ebert_curry
