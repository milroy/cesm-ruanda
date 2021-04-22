!===========================================================================
! 1 interface to the UNIFIED CONVECTION SCHEME (UNICON)
!
! The USE_UNICON macro converts this module to a stub interface which allows
! 1 to be built without the unicon and unicon_utils modules.
!
!===========================================================================

module unicon_cam

use shr_kind_mod,     only: r8 => shr_kind_r8, i4 => shr_kind_i4
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: rair, cpair, gravit, latvap, latice, zvir, mwdry

use constituents,     only: pcnst, cnst_add, qmin, cnst_get_type_byind, cnst_get_ind, cnst_name
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use camsrfexch,       only: cam_in_t
use physics_buffer,   only: pbuf_add_field, dtype_r8, dyn_time_lvls, pbuf_old_tim_idx, &
                            physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_set_field
use cam_history,      only: outfld, addfld, phys_decomp
use time_manager,     only: is_first_step
use abortutils,       only: endrun


implicit none
private
save

public :: &
   unicon_cam_register,    &
   unicon_cam_init,        &
   unicon_implements_cnst, &
   unicon_init_cnst,       &
   unicon_out_t,           &
   unicon_cam_tend,        &
   unicon_cam_org_diags

real(r8) :: xlv    ! Latent heat of vaporization
real(r8) :: xlf    ! Latent heat of fusion
real(r8) :: xls    ! Latent heat of sublimation
real(r8) :: cp     ! Specific heat of dry air

integer, parameter :: &
   nseg = 1,      &! Number of updraft segments [ # ]
   mix  = pcols,  &! Maximum number of columns
   mkx  = pver,   &! Number of vertical layers
   ncnst = pcnst   ! Number of advected constituents

! For advecting organization-related variables
integer, parameter :: n_org = 5                      ! Number of constituents
character(len=8), dimension(n_org), parameter :: &   ! Constituent names
   cnst_names = (/'ORGawk  ','ORGthl  ','ORGqto  ','ORGuoo  ','ORGvoo  '/)

integer :: awk_cnst_ind, thl_cnst_ind, qt_cnst_ind, u_cnst_ind, v_cnst_ind 

! fields added to physics buffer by this module
integer :: &
   cushavg_idx, &
   cuorg_idx, &
   awk_PBL_idx, &
   delta_thl_PBL_idx, &
   delta_qt_PBL_idx, &
   delta_u_PBL_idx, &
   delta_v_PBL_idx, &
   delta_tr_PBL_idx, &
   cu_cmfr_idx, &
   cu_thlr_idx, &
   cu_qtr_idx, &
   cu_ur_idx, &
   cu_vr_idx, &
   cu_qlr_idx, &
   cu_qir_idx, &
   cu_trr_idx

! fields expected to be in the physics buffer
integer :: &
   ast_idx   = -1,   &
   tke_idx   = -1,   &
   bprod_idx = -1, &
   kpblh_idx  = -1, &
   pblh_idx   = -1, &
   went_idx   = -1,   &
   cush_idx   = -1, &
   shfrc_idx  = -1, &
   icwmrsh_idx = -1, &
   rprdsh_idx  = -1, &
   prec_sh_idx = -1, &
   snow_sh_idx = -1, &
   nevapr_shcu_idx = -1, &
   am_evp_st_idx = -1, &    !  Evaporation area of stratiform precipitation [fraction]
   evprain_st_idx = -1, &   !  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
   evpsnow_st_idx = -1      !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.


! constituent indices
integer :: ixcldliq, ixcldice, ixnumliq, ixnumice

! unicon output fields
type unicon_out_t
   real(r8) :: cmfmc(mix,mkx+1)   ! Upward     convective mass flux at the interface [ kg / s / m2 ]
   real(r8) :: slflx(mix,mkx+1)   ! Net upward convective flux of liquid static energy [ J / s / m2 ]
   real(r8) :: qtflx(mix,mkx+1)   ! Net upward convective flux of total specific humidity [ kg / s / m2 ]
   real(r8) :: rqc(mix,mkx)       ! Prod of suspended LWC+IWC by expel of excessive in-cumulus condensate [ kg / kg / s ] > 0
   real(r8) :: rliq(mix)          ! Vertical integral of 'rqc' in flux unit [ m / s ]
   real(r8) :: cnt(mix)           ! Cloud top  interface index ( ki = kpen )
   real(r8) :: cnb(mix)           ! Cloud base interface index ( ki = krel-1 )
end type unicon_out_t

! logical array to identify constituents that are mode number concentrations
logical :: cnst_is_mam_num(ncnst)
! logical array to identify constituents that are mode specie mass mixing ratios
logical :: cnst_is_mam_mmr(ncnst)

!==================================================================================================
contains
!==================================================================================================
  
subroutine unicon_cam_register

! Register fields in the constituent array and the physics buffer.


   ! Jun.02.2012. Sungsu for advecting organization-related horizontal heterogeneity
   !              within PBL. 
   !              For the time being, in order to save computation time, advection of aerosol perturbations 
   !              are simply neglected.


end subroutine unicon_cam_register

!==================================================================================================

subroutine unicon_cam_init(pbuf2d)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   character(len=*), parameter :: sub='unicon_cam_init: '
   integer :: i, icnst, j, l, nmodes, nspec

   character(len=8)  :: units
   character(len=30) :: varname
   character(len=60) :: surname
   character(len=2)  :: numcha
   integer           :: msfc
   ! ------------------------------------------------------------------------------------------- !


end subroutine unicon_cam_init

!==================================================================================================

function unicon_implements_cnst(name)

   ! Return true if specified constituent is implemented by this package

   character(len=*), intent(in) :: name  ! constituent name
   logical :: unicon_implements_cnst     ! return value

   integer :: m
   !-----------------------------------------------------------------------

   unicon_implements_cnst = .false.


end function unicon_implements_cnst

!==================================================================================================

subroutine unicon_init_cnst(name, q, gcid)

   ! Initialize constituents if they are not read from the initial file

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------


end subroutine unicon_init_cnst

!==================================================================================================

subroutine unicon_cam_tend(dt, state, cam_in, sgh30, &
                           pbuf, ptend, out)


   ! ---------------------- !
   ! Input-output Arguments !
   ! ---------------------- !

   real(r8),                  intent(in)  :: dt         ! Time step [s]
   type(physics_state),       intent(in)  :: state      ! Physics state variables
   type(cam_in_t),            intent(in)  :: cam_in     ! import state
   real(r8),                  intent(in)  :: sgh30(mix) ! Std dev of subgrid topographic height at 30 s horizontal area [ meter ]
   type(physics_buffer_desc), pointer     :: pbuf(:)    ! physics buffer
   type(physics_ptend),       intent(out) :: ptend      ! parameterization tendencies
   type(unicon_out_t),        intent(out) :: out        ! parameterization outputs

   ! -------------------------------------------------------- !
   ! Internal output and local variables for positive tracers !
   ! -------------------------------------------------------- !
   
   real(r8) :: sten_ori(mix,mkx)         !  Tendency of dry static energy [ J / kg / s ]
   real(r8) :: qvten_ori(mix,mkx)        !  Tendency of water vapor specific humidity [ kg / kg / s ]
   real(r8) :: qlten_ori(mix,mkx)        !  Tendency of liquid water mixing ratio [ kg / kg / s ]
   real(r8) :: qiten_ori(mix,mkx)        !  Tendency of ice mixing ratio [ kg / kg / s ]
   real(r8) :: trten_ori(mix,mkx,ncnst)  !  Tendency of tracers [ # / kg / s, kg / kg / s ]

   real(r8) :: slten_pos_inv(mix,mkx)    ! 
   real(r8) :: qtten_pos_inv(mix,mkx)    ! 
   real(r8) :: uten_pos_inv(mix,mkx)     ! 
   real(r8) :: vten_pos_inv(mix,mkx)     ! 
   real(r8) :: sten_pos_inv(mix,mkx)     ! 
   real(r8) :: qvten_pos_inv(mix,mkx)    ! 
   real(r8) :: qlten_pos_inv(mix,mkx)    ! 
   real(r8) :: qiten_pos_inv(mix,mkx)    ! 
   real(r8) :: trten_pos_inv(mix,mkx,ncnst)

   ! --------------- !
   ! Local variables !
   ! --------------- !

   integer :: iend
   integer :: lchnk     
   integer :: itim

   ! fields in physics buffer
   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx)
      ast0_inv    ! Physical stratus fraction [ fraction ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx+1)
      tke0_inv,    &! TKE at the interface [ m2/s2 ]
      bprod0_inv    ! Buoyancy production at the interface [ m2/s3 ]

   integer(i4), pointer, dimension(:) :: & ! (mix)
      kpblh_inv       ! Layer index with PBL top in it or at the base interface

   real(r8), pointer, dimension(:) :: & ! (mix)
      pblh,          &! PBL top height [m]
      went,          &! Entrainment rate at the PBL top interface directly from UW PBL scheme [ m / s ]. went = 0 for STL.
      cush,          &! Cumulus top height [ m ]
      cushavg,       &! Mean cumulus top height weighted by updraft mass flux at surface [ m ]
      cuorg,         &! Convective organization parameter [ 0-1 ]

      awk_PBL,       &! Wake area within PBL [ 0 - 1 ]
      delta_thl_PBL, &! Diff of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
      delta_qt_PBL,  &! Diff of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
      delta_u_PBL,   &! Diff of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
      delta_v_PBL     ! Diff of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,ncnst)
      delta_tr_PBL ! Diff of tr  between off-wake region and grid-mean value avg over the PBL [ kg/kg, #/kg ]

   real(r8), dimension(mix,mkx) :: & ! (mix,mkx)
      cu_cmfum  ! The mass involved in the updraft buoyancy sorting at the previous time step [ kg/s/m2 ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx)
      cu_cmfr, &! The detrained mass from convective up and downdraft at the previous time step [ kg/s/m2 ]
      cu_thlr, &! Mass-flux wghted mean 'thl' of detrained mass from conv up and downdraft at prev step [ K ]
      cu_qtr,  &! Mass-flux wghted mean 'qt'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      cu_ur,   &! Mass-flux wghted mean 'u'   of detrained mass from conv up and downdraft at prev step [ m/s ]
      cu_vr,   &! Mass-flux wghted mean 'v'   of detrained mass from conv up and downdraft at prev step [ m/s ]
      cu_qlr,  &! Mass-flux wghted mean 'ql'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      cu_qir    ! Mass-flux wghted mean 'qi'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]

   real(r8), pointer, dimension(:,:,:) :: & ! (mix,mkx,ncnst)
      cu_trr    ! Mass-flux wghted mean 'tr'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]

   real(r8), dimension(mix,mkx) :: & ! (mix,mkx)
      cu_cmfrd,&! The amount of detrained mass from convective downdraft at the previous time step [ kg/s/m2 ]
      cu_thlrd,&! Mass-flux wghted mean 'thl' of detrained mass from conv downdraft at prev step [ K ]
      cu_qtrd, &! Mass-flux wghted mean 'qt'  of detrained mass from conv downdraft at prev step [ kg/kg ]
      cu_urd,  &! Mass-flux wghted mean 'u'   of detrained mass from conv downdraft at prev step [ m/s ]
      cu_vrd,  &! Mass-flux wghted mean 'v'   of detrained mass from conv downdraft at prev step [ m/s ]
      cu_qlrd, &! Mass-flux wghted mean 'ql'  of detrained mass from conv downdraft at prev step [ kg/kg ]
      cu_qird   ! Mass-flux wghted mean 'qi'  of detrained mass from conv downdraft at prev step [ kg/kg ]

   real(r8), dimension(mix,mkx,ncnst) :: & ! (mix,mkx,ncnst)
      cu_trrd   ! Mass-flux wghted mean 'tr'  of detrained mass from conv downdraft at prev step [ kg/kg ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx)
      shfrc,     &! Convective updraft fractional area
      icwmrsh,   &! In-cloud LWC+IWC within convective updraft
      rprdsh,    &! Prod of rain+snow by lateral expels of cumulus condensate [ kg / kg / s ]
      evapc_inv, &! Evaporation rate of convective precipitation within environment [ kg/kg/s ]
      am_evp_st_inv,   &!  Evaporation area of stratiform precipitation [fraction]
      evprain_st_inv,  &!  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
      evpsnow_st_inv    !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.

   real(r8), pointer, dimension(:) :: & ! (mix)
      precip,    &!  Precipitation flux at surface in flux unit [ m / s ]
      snow        !  Snow flux at surface in flux unit [ m / s ]

   real(r8) :: ps0(mix,0:mkx)            !  Environmental pressure at full sigma levels
   real(r8) :: zs0(mix,0:mkx)            !  Environmental height at full sigma levels
   real(r8) :: p0(mix,mkx)               !  Environmental pressure at half sigma levels
   real(r8) :: z0(mix,mkx)               !  Environmental height at half sigma levels
   real(r8) :: dp0(mix,mkx)              !  Environmental layer pressure thickness
   real(r8) :: dpdry0(mix,mkx)           !  Environmental layer dry pressure thickness
   real(r8) :: u0(mix,mkx)               !  Environmental zonal wind
   real(r8) :: v0(mix,mkx)               !  Environmental meridional wind
   real(r8) :: qv0(mix,mkx)              !  Environmental specific humidity
   real(r8) :: ql0(mix,mkx)              !  Environmental liquid water mixing ratio
   real(r8) :: qi0(mix,mkx)              !  Environmental ice mixing ratio
   real(r8) :: tr0(mix,mkx,ncnst)        !  Environmental tracers [ #/kg, kg/kg ]
   real(r8) :: t0(mix,mkx)               !  Environmental temperature
   real(r8) :: s0(mix,mkx)               !  Environmental dry static energy
   real(r8) :: ast0(mix,mkx)             !  Physical stratiform cloud fraction [ fraction ]
   real(r8) :: tke0(mix,0:mkx)           !  TKE [ m2/s2 ]
   real(r8) :: bprod0(mix,0:mkx)         !  Buoyancy production [ m2/s3 ]
   real(r8) :: am_evp_st(mix,mkx)        !  Evaporation area of stratiform precipitation [fraction]
   real(r8) :: evprain_st(mix,mkx)       !  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
   real(r8) :: evpsnow_st(mix,mkx)       !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.

   integer(i4) :: kpblh(mix)             !  Layer index with PBL top in it or at the base interface


   real(r8) :: am_u(mix,mkx)             !  Convective updraft fractional area
   real(r8) :: qlm_u(mix,mkx)            !  In-cloud LWC within convective updraft [ kg / kg ]
   real(r8) :: qim_u(mix,mkx)            !  In-cloud IWC within convective updraft [ kg / kg ]

   real(r8) :: am_d(mix,mkx)             !  Convective downdraft fractional area
   real(r8) :: qlm_d(mix,mkx)            !  In-cloud LWC within downdraft updraft [ kg / kg ]
   real(r8) :: qim_d(mix,mkx)            !  In-cloud IWC within downdraft updraft [ kg / kg ]

   real(r8) :: cmf_u(mix,0:mkx)          !  Upward     convective mass flux at the interface [ kg / s / m2 ]
   real(r8) :: slflx(mix,0:mkx)          !  Net upward convective flux of liquid static energy [ J / s / m2 ]
   real(r8) :: qtflx(mix,0:mkx)          !  Net upward convective flux of total specific humidity [ kg / s / m2 ]

   real(r8) :: qvten(mix,mkx)            !  Tendency of water vapor specific humidity [ kg / kg / s ]
   real(r8) :: qlten(mix,mkx)            !  Tendency of liquid water mixing ratio [ kg / kg / s ]
   real(r8) :: qiten(mix,mkx)            !  Tendency of ice mixing ratio [ kg / kg / s ]
   real(r8) :: trten(mix,mkx,ncnst)      !  Tendency of tracers [ # / kg / s, kg / kg / s ]

   real(r8) :: sten(mix,mkx)             !  Tendency of dry static energy [ J / kg / s ]
   real(r8) :: uten(mix,mkx)             !  Tendency of zonal wind [ m / s / s ]
   real(r8) :: vten(mix,mkx)             !  Tendency of meridional wind [ m / s / s ]

   real(r8) :: qrten(mix,mkx)            !  Production rate of rain by lateral expels of cumulus condensate [ kg / kg / s ]
   real(r8) :: qsten(mix,mkx)            !  Production rate of snow by lateral expels of cumulus condensate [ kg / kg / s ]

   real(r8) :: evapc(mix,mkx)            !  Evaporation rate of convective precipitation within environment [ kg/kg/s ]

   real(r8) :: rqc(mix,mkx)              !  Production rate of detrained LWC+IWC [ kg / kg / s ] > 0
   real(r8) :: rqc_l(mix,mkx)            !  Production rate of detrained LWC     [ kg / kg / s ] > 0
   real(r8) :: rqc_i(mix,mkx)            !  Production rate of detrained IWC     [ kg / kg / s ] > 0
   real(r8) :: rnc_l(mix,mkx)            !  Production rate of detrained droplet number of cloud liquid droplets [ # / kg / s ] > 0
   real(r8) :: rnc_i(mix,mkx)            !  Production rate of detrained droplet number of cloud    ice droplets [ # / kg / s ] > 0

   real(r8) :: cnt(mix)                  !  Cloud top  interface index ( ki = kpen )
   real(r8) :: cnb(mix)                  !  Cloud base interface index ( ki = krel-1 )

   real(r8) :: pdel0(mix,mkx)            !  Environmental pressure thickness ( either dry or moist ) [ Pa ]
   real(r8) :: trmin                     !  Minimum concentration of tracer [ # / kg ]

   ! For prognostically updated state variables

   real(r8) :: qv0_c(mix,mkx)            !  Environmental specific humidity
   real(r8) :: ql0_c(mix,mkx)            !  Environmental liquid water mixing ratio
   real(r8) :: qi0_c(mix,mkx)            !  Environmental ice mixing ratio
   real(r8) :: t0_c(mix,mkx)             !  Environmental temperature
   real(r8) :: s0_c(mix,mkx)             !  Environmental dry static energy
   real(r8) :: tr0_c(mix,mkx,ncnst)      !  Environmental tracers [ # / kg, kg / kg ]

   ! Layer index variables
   integer  :: k                         !  Vertical index for local fields 
   integer  :: k_inv                     !  Vertical index for incoming fields
   integer  :: mt                        !  Tracer index [ no ]
   integer  :: m

   ! For aerosol tendency output
   character(len=30) :: varname

   logical :: lq(ncnst)

   ! --------- !
   ! Main body !
   ! --------- !


end subroutine unicon_cam_tend

subroutine unicon_cam_org_diags(state, pbuf)

   ! ------------------------------------------------------------------------
   ! Insert the organization-related heterogeneities computed inside the
   ! UNICON into the tracer arrays here before performing advection.
   ! This is necessary to prevent any modifications of organization-related
   ! heterogeneities by non convection-advection process, such as
   ! dry and wet deposition of aerosols, MAM, etc.
   ! Again, note that only UNICON and advection schemes are allowed to
   ! changes to organization at this stage, although we can include the
   ! effects of other physical processes in future.
   ! ------------------------------------------------------------------------

   ! Arguments
   type(physics_state), intent(inout) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   ! Local variables
   real(r8), pointer, dimension(:  ) :: awk_PBL
   real(r8), pointer, dimension(:  ) :: delta_thl_PBL
   real(r8), pointer, dimension(:  ) :: delta_qt_PBL
   real(r8), pointer, dimension(:  ) :: delta_u_PBL
   real(r8), pointer, dimension(:  ) :: delta_v_PBL

   integer :: i, itim, ncol
   ! ------------------------------------------------------------------------


end subroutine unicon_cam_org_diags

end module unicon_cam
