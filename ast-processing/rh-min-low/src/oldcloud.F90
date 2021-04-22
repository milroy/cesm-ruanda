module oldcloud

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag, get_sw_spectral_boundaries
use abortutils,       only: endrun
use cam_history,      only: outfld
use rad_constituents, only: iceopticsfile, liqopticsfile
use ebert_curry,      only: scalefactor

implicit none
private
save

public :: &
 oldcloud_init, oldcloud_lw, old_liq_get_rad_props_lw, old_ice_get_rad_props_lw

integer :: nmu, nlambda
real(r8), allocatable :: g_mu(:)           ! mu samples on grid
real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
real(r8), allocatable :: ext_sw_liq(:,:,:)
real(r8), allocatable :: ssa_sw_liq(:,:,:)
real(r8), allocatable :: asm_sw_liq(:,:,:)
real(r8), allocatable :: abs_lw_liq(:,:,:)

integer :: n_g_d
real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
real(r8), allocatable :: ext_sw_ice(:,:)
real(r8), allocatable :: ssa_sw_ice(:,:)
real(r8), allocatable :: asm_sw_ice(:,:)
real(r8), allocatable :: abs_lw_ice(:,:)

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
   integer ::   iciwp_idx   = 0 
   integer ::   iclwp_idx   = 0 
   integer ::   cld_idx     = 0 
   integer ::   rel_idx     = 0 
   integer ::   rei_idx     = 0 

! indexes into constituents for old optics
   integer :: &
        ixcldice,           & ! cloud ice water index
        ixcldliq              ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine oldcloud_init()

   use constituents,   only: cnst_get_ind

   integer :: err

   iciwp_idx  = pbuf_get_index('ICIWP',errcode=err)
   iclwp_idx  = pbuf_get_index('ICLWP',errcode=err)
   cld_idx    = pbuf_get_index('CLD')
   rel_idx    = pbuf_get_index('REL')
   rei_idx    = pbuf_get_index('REI')

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   return

end subroutine oldcloud_init

!==============================================================================
! Private methods
!==============================================================================

! Subprogram not used subroutine old_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp)
! Subprogram not used 
! Subprogram not used    use physconst, only: gravit
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in) :: state
! Subprogram not used    type(physics_buffer_desc),pointer :: pbuf(:)
! Subprogram not used 
! Subprogram not used    real(r8),intent(out) :: liq_tau    (nswbands,pcols,pver) ! extinction optical depth
! Subprogram not used    real(r8),intent(out) :: liq_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
! Subprogram not used    real(r8),intent(out) :: liq_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
! Subprogram not used    real(r8),intent(out) :: liq_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w
! Subprogram not used    logical, intent(in) :: oldliqwp
! Subprogram not used 
! Subprogram not used    real(r8), pointer, dimension(:,:) :: rel
! Subprogram not used    real(r8), pointer, dimension(:,:) :: cldn
! Subprogram not used    real(r8), pointer, dimension(:,:) :: tmpptr
! Subprogram not used    real(r8), dimension(pcols,pver) :: cliqwp
! Subprogram not used    real(r8), dimension(nswbands) :: wavmin
! Subprogram not used    real(r8), dimension(nswbands) :: wavmax
! Subprogram not used 
! Subprogram not used    ! Minimum cloud amount (as a fraction of the grid-box area) to 
! Subprogram not used    ! distinguish from clear sky
! Subprogram not used    real(r8), parameter :: cldmin = 1.0e-80_r8
! Subprogram not used 
! Subprogram not used    ! Decimal precision of cloud amount (0 -> preserve full resolution;
! Subprogram not used    ! 10^-n -> preserve n digits of cloud amount)
! Subprogram not used    real(r8), parameter :: cldeps = 0.0_r8
! Subprogram not used 
! Subprogram not used    ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
! Subprogram not used    ! Parameterization for the Shortwave Properties of Water Clouds' JAS
! Subprogram not used    ! vol. 46 may 1989 pp 1419-1427)
! Subprogram not used    real(r8) :: abarl(4) = &  ! A coefficient for extinction optical depth
! Subprogram not used       (/ 2.817e-02_r8, 2.682e-02_r8,2.264e-02_r8,1.281e-02_r8/)
! Subprogram not used    real(r8) :: bbarl(4) = &  ! B coefficient for extinction optical depth
! Subprogram not used       (/ 1.305_r8    , 1.346_r8    ,1.454_r8    ,1.641_r8    /)
! Subprogram not used    real(r8) :: cbarl(4) = &  ! C coefficient for single scat albedo
! Subprogram not used       (/-5.62e-08_r8 ,-6.94e-06_r8 ,4.64e-04_r8 ,0.201_r8    /)
! Subprogram not used    real(r8) :: dbarl(4) = &  ! D coefficient for single  scat albedo
! Subprogram not used       (/ 1.63e-07_r8 , 2.35e-05_r8 ,1.24e-03_r8 ,7.56e-03_r8 /)
! Subprogram not used    real(r8) :: ebarl(4) = &  ! E coefficient for asymmetry parameter
! Subprogram not used       (/ 0.829_r8    , 0.794_r8    ,0.754_r8    ,0.826_r8    /)
! Subprogram not used    real(r8) :: fbarl(4) = &  ! F coefficient for asymmetry parameter
! Subprogram not used       (/ 2.482e-03_r8, 4.226e-03_r8,6.560e-03_r8,4.353e-03_r8/)
! Subprogram not used 
! Subprogram not used    real(r8) :: abarli        ! A coefficient for current spectral band
! Subprogram not used    real(r8) :: bbarli        ! B coefficient for current spectral band
! Subprogram not used    real(r8) :: cbarli        ! C coefficient for current spectral band
! Subprogram not used    real(r8) :: dbarli        ! D coefficient for current spectral band
! Subprogram not used    real(r8) :: ebarli        ! E coefficient for current spectral band
! Subprogram not used    real(r8) :: fbarli        ! F coefficient for current spectral band
! Subprogram not used 
! Subprogram not used    ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
! Subprogram not used    ! greater than 20 micro-meters
! Subprogram not used 
! Subprogram not used    integer :: ns, i, k, indxsl, Nday
! Subprogram not used    integer :: lchnk, itim_old
! Subprogram not used    real(r8) :: tmp1l, tmp2l, tmp3l, g
! Subprogram not used    real(r8) :: kext(pcols,pver)
! Subprogram not used    real(r8), pointer, dimension(:,:) :: iclwpth
! Subprogram not used 
! Subprogram not used    Nday = state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    itim_old = pbuf_old_tim_idx()
! Subprogram not used    call pbuf_get_field(pbuf, cld_idx,cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
! Subprogram not used    call pbuf_get_field(pbuf, rel_idx,rel)
! Subprogram not used 
! Subprogram not used    if (oldliqwp) then
! Subprogram not used      do k=1,pver
! Subprogram not used         do i = 1,Nday
! Subprogram not used            cliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/(gravit*max(0.01_r8,cldn(i,k)))
! Subprogram not used         end do
! Subprogram not used      end do
! Subprogram not used    else
! Subprogram not used      if (iclwp_idx<0) then 
! Subprogram not used         call endrun('old_liquid_optics_sw: oldliqwp must be set to true since ICLWP was not found in pbuf')
! Subprogram not used      endif
! Subprogram not used      ! The following is the eventual target specification for in cloud liquid water path.
! Subprogram not used      call pbuf_get_field(pbuf, iclwp_idx, tmpptr)
! Subprogram not used      cliqwp = tmpptr
! Subprogram not used    endif
! Subprogram not used   
! Subprogram not used    call get_sw_spectral_boundaries(wavmin,wavmax,'microns')
! Subprogram not used 
! Subprogram not used    do ns = 1, nswbands
! Subprogram not used       ! Set index for cloud particle properties based on the wavelength,
! Subprogram not used       ! according to A. Slingo (1989) equations 1-3:
! Subprogram not used       ! Use index 1 (0.25 to 0.69 micrometers) for visible
! Subprogram not used       ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
! Subprogram not used       ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
! Subprogram not used       ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
! Subprogram not used       if(wavmax(ns) <= 0.7_r8) then
! Subprogram not used          indxsl = 1
! Subprogram not used       else if(wavmax(ns) <= 1.25_r8) then
! Subprogram not used          indxsl = 2
! Subprogram not used       else if(wavmax(ns) <= 2.38_r8) then
! Subprogram not used          indxsl = 3
! Subprogram not used       else if(wavmin(ns) > 2.38_r8) then
! Subprogram not used          indxsl = 4
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       ! Set cloud extinction optical depth, single scatter albedo,
! Subprogram not used       ! asymmetry parameter, and forward scattered fraction:
! Subprogram not used       abarli = abarl(indxsl)
! Subprogram not used       bbarli = bbarl(indxsl)
! Subprogram not used       cbarli = cbarl(indxsl)
! Subprogram not used       dbarli = dbarl(indxsl)
! Subprogram not used       ebarli = ebarl(indxsl)
! Subprogram not used       fbarli = fbarl(indxsl)
! Subprogram not used 
! Subprogram not used       do k=1,pver
! Subprogram not used          do i=1,Nday
! Subprogram not used 
! Subprogram not used             ! note that optical properties for liquid valid only
! Subprogram not used             ! in range of 4.2 > rel > 16 micron (Slingo 89)
! Subprogram not used             if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
! Subprogram not used                tmp1l = abarli + bbarli/min(max(4.2_r8,rel(i,k)),16._r8)
! Subprogram not used                liq_tau(ns,i,k) = 1000._r8*cliqwp(i,k)*tmp1l
! Subprogram not used             else
! Subprogram not used                liq_tau(ns,i,k) = 0.0_r8
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             tmp2l = 1._r8 - cbarli - dbarli*min(max(4.2_r8,rel(i,k)),16._r8)
! Subprogram not used             tmp3l = fbarli*min(max(4.2_r8,rel(i,k)),16._r8)
! Subprogram not used             ! Do not let single scatter albedo be 1.  Delta-eddington solution
! Subprogram not used             ! for non-conservative case has different analytic form from solution
! Subprogram not used             ! for conservative case, and raddedmx is written for non-conservative case.
! Subprogram not used             liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * min(tmp2l,.999999_r8)
! Subprogram not used             g = ebarli + tmp3l
! Subprogram not used             liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g
! Subprogram not used             liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g
! Subprogram not used 
! Subprogram not used          end do ! End do i=1,Nday
! Subprogram not used       end do    ! End do k=1,pver
! Subprogram not used    end do ! nswbands
! Subprogram not used 
! Subprogram not used    !call outfld('CL_OD_SW_OLD',liq_tau(idx_sw_diag,:,:), pcols, lchnk)
! Subprogram not used    !call outfld('REL_OLD',rel(:,:), pcols, lchnk)
! Subprogram not used    !call outfld('CLWPTH_OLD',cliqwp(:,:), pcols, lchnk)
! Subprogram not used    !call outfld('KEXT_OLD',kext(:,:), pcols, lchnk)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used end subroutine old_liquid_optics_sw
!==============================================================================

! Subprogram not used subroutine old_ice_optics_sw   (state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp)
! Subprogram not used 
! Subprogram not used    use physconst, only: gravit
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in) :: state
! Subprogram not used    type(physics_buffer_desc),pointer :: pbuf(:)
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
! Subprogram not used         call endrun('old_ice_optics_sw: oldicewp must be set to true since ICIWP was not found in pbuf')
! Subprogram not used      endif
! Subprogram not used      call pbuf_get_field(pbuf, iciwp_idx, tmpptr)
! Subprogram not used      cicewp(1:pcols,1:pver) =  1000.0_r8*tmpptr
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
! Subprogram not used       else if(wavmin(ns) > 2.38_r8) then
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
! Subprogram not used end subroutine old_ice_optics_sw
!==============================================================================

! Subprogram not used subroutine oldcloud_lw(state,pbuf,cld_abs_od,oldwp)
! Subprogram not used    use physconst, only: gravit
! Subprogram not used    type(physics_state), intent(in) :: state
! Subprogram not used    type(physics_buffer_desc),pointer :: pbuf(:)
! Subprogram not used    real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
! Subprogram not used    logical,intent(in)              :: oldwp ! use old definition of waterpath
! Subprogram not used 
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
! Subprogram not used    real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth
! Subprogram not used 
! Subprogram not used     real(r8) :: kabs, kabsi
! Subprogram not used     real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
! Subprogram not used     parameter (kabsl = 0.090361_r8)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    ncol = state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    itim_old  =  pbuf_old_tim_idx()
! Subprogram not used    call pbuf_get_field(pbuf, rei_idx,   rei)
! Subprogram not used    call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
! Subprogram not used 
! Subprogram not used    if (oldwp) then
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
! Subprogram not used          call endrun('oldcloud_lw: oldwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
! Subprogram not used       endif
! Subprogram not used       call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
! Subprogram not used       call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
! Subprogram not used       do k=1,pver
! Subprogram not used          do i = 1,ncol
! Subprogram not used             cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
! Subprogram not used             ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    do k=1,pver
! Subprogram not used        do i=1,ncol
! Subprogram not used 
! Subprogram not used           !note that optical properties for ice valid only
! Subprogram not used           !in range of 13 > rei > 130 micron (Ebert and Curry 92)
! Subprogram not used           kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
! Subprogram not used           kabs = kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
! Subprogram not used           !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
! Subprogram not used           cldtau(i,k) = kabs*cwp(i,k)
! Subprogram not used        end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do lwband = 1,nlwbands
! Subprogram not used       cld_abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used end subroutine oldcloud_lw

!==============================================================================
! Subprogram not used subroutine old_liq_get_rad_props_lw(state, pbuf, abs_od, oldliqwp)
! Subprogram not used    use physconst, only: gravit
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in)  :: state
! Subprogram not used    type(physics_buffer_desc),pointer  :: pbuf(:)
! Subprogram not used    real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)
! Subprogram not used    logical, intent(in) :: oldliqwp
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
! Subprogram not used    real(r8) :: kabs, kabsi
! Subprogram not used    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
! Subprogram not used    parameter (kabsl = 0.090361_r8)
! Subprogram not used 
! Subprogram not used    real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth
! Subprogram not used 
! Subprogram not used    ncol=state%ncol
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    itim_old  =  pbuf_old_tim_idx()
! Subprogram not used    call pbuf_get_field(pbuf, rei_idx,   rei)
! Subprogram not used    call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
! Subprogram not used 
! Subprogram not used    if (oldliqwp) then
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
! Subprogram not used          call endrun('old_liq_get_rad_props_lw: oldliqwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
! Subprogram not used       endif
! Subprogram not used       call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
! Subprogram not used       call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
! Subprogram not used       do k=1,pver
! Subprogram not used          do i = 1,ncol
! Subprogram not used             cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
! Subprogram not used             ficemr(i,k) = 1000.0 * iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    do k=1,pver
! Subprogram not used        do i=1,ncol
! Subprogram not used 
! Subprogram not used           ! Note from Andrew Conley:
! Subprogram not used           !  Optics for RK no longer supported, This is constructed to get
! Subprogram not used           !  close to bit for bit.  Otherwise we could simply use liquid water path
! Subprogram not used           !note that optical properties for ice valid only
! Subprogram not used           !in range of 13 > rei > 130 micron (Ebert and Curry 92)
! Subprogram not used           kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
! Subprogram not used           kabs = kabsl*(1._r8-ficemr(i,k)) ! + kabsi*ficemr(i,k)
! Subprogram not used           !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
! Subprogram not used           cldtau(i,k) = kabs*cwp(i,k)
! Subprogram not used        end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do lwband = 1,nlwbands
! Subprogram not used       abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used end subroutine old_liq_get_rad_props_lw
!==============================================================================

! Subprogram not used subroutine old_ice_get_rad_props_lw(state, pbuf, abs_od, oldicewp)
! Subprogram not used    use physconst, only: gravit
! Subprogram not used    type(physics_state), intent(in)  :: state
! Subprogram not used    type(physics_buffer_desc),pointer  :: pbuf(:)
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
! Subprogram not used    call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
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
! Subprogram not used          call endrun('old_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
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
! Subprogram not used end subroutine old_ice_get_rad_props_lw
!==============================================================================

! Subprogram not used subroutine cloud_total_vis_diag_out(lchnk, nnite, idxnite, tau, radsuffix)
! Subprogram not used 
! Subprogram not used    ! output total aerosol optical depth for the visible band
! Subprogram not used 
! Subprogram not used    use cam_history, only: outfld
! Subprogram not used    use cam_history_support, only : fillvalue
! Subprogram not used 
! Subprogram not used    integer,          intent(in) :: lchnk
! Subprogram not used    integer,          intent(in) :: nnite          ! number of night columns
! Subprogram not used    integer,          intent(in) :: idxnite(nnite) ! local column indices of night columns
! Subprogram not used    real(r8),         intent(in) :: tau(:,:)
! Subprogram not used    character(len=*), intent(in) :: radsuffix ! identifies whether the radiation call
! Subprogram not used                                              ! is for the climate calc or a diagnostic calc
! Subprogram not used  
! Subprogram not used    ! Local variables
! Subprogram not used    integer  :: i
! Subprogram not used    real(r8) :: tmp(pcols)
! Subprogram not used    !-----------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ! compute total aerosol optical depth output where only daylight columns
! Subprogram not used    tmp(:) = sum(tau(:,:), 2)
! Subprogram not used    do i = 1, nnite
! Subprogram not used       tmp(idxnite(i)) = fillvalue
! Subprogram not used    end do
! Subprogram not used    !call outfld('cloudOD_v'//trim(radsuffix), tmp, pcols, lchnk)
! Subprogram not used 
! Subprogram not used end subroutine cloud_total_vis_diag_out

!==============================================================================

end module oldcloud
