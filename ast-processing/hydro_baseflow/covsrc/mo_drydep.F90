module mo_drydep

  !---------------------------------------------------------------------
  !       ... Dry deposition velocity input data and code for netcdf input
  !---------------------------------------------------------------------

!LKE (10/11/2010): added HCN, CH3CN, HCOOH

  use shr_kind_mod, only : r8 => shr_kind_r8, shr_kind_cl
  use chem_mods,    only : gas_pcnst
  use pmgrid,       only : plev, plevp
  use spmd_utils,   only : masterproc, iam
  use ppgrid,       only : pcols, begchunk, endchunk
  use mo_tracname,  only : solsym
  use abortutils,   only : endrun
  use ioFileMod,    only : getfil

  use mpishorthand, only : mpicom, mpir8, mpiint, mpilog

  use pio
  use cam_pio_utils,only : cam_pio_openfile
  use cam_logfile,  only : iulog
  use dyn_grid,     only : get_dyn_grid_parm, get_horiz_grid_d
  use scamMod,      only:  single_column

  use seq_drydep_mod, only : nddvels =>  n_drydep, drydep_list, mapping
  use physconst,    only : karman

  implicit none

  save

  interface drydep_inti
     module procedure dvel_inti_table
     module procedure dvel_inti_xactive
     module procedure dvel_inti_fromlnd
  end interface

  interface drydep
     module procedure drydep_table
     module procedure drydep_xactive
     module procedure drydep_fromlnd
  end interface

  private
  public :: drydep_inti, drydep, set_soilw, chk_soilw, has_drydep
  public :: drydep_update
  public :: n_land_type, fraction_landuse, drydep_srf_file

  real(r8)              :: dels
  real(r8), allocatable :: days(:)          ! day of year for soilw
  real(r8), allocatable :: dvel(:,:,:,:)    ! depvel array interpolated to model grid
  real(r8), allocatable :: dvel_interp(:,:,:) ! depvel array interpolated to grid and time
  integer :: last, next                     ! day indicies
  integer :: ndays                          ! # of days in soilw file
  integer :: map(gas_pcnst)                 ! indices for drydep species
  integer :: nspecies                       ! number of depvel species in input file

  integer :: pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
             h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
             ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, eooh_ndx, &
             c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
             no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
             hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
             xooh_ndx, ch3cho_ndx, isopooh_ndx
  integer :: alkooh_ndx, mekooh_ndx, tolooh_ndx, terpooh_ndx, ch3cooh_ndx
  integer :: soa_ndx, so4_ndx, cb1_ndx, cb2_ndx, oc1_ndx, oc2_ndx, nh3_ndx, nh4no3_ndx, &
             sa1_ndx, sa2_ndx, sa3_ndx, sa4_ndx, nh4_ndx
  integer :: soam_ndx, soai_ndx, soat_ndx, soab_ndx, soax_ndx, &
             sogm_ndx, sogi_ndx, sogt_ndx, sogb_ndx, sogx_ndx

  logical :: alkooh_dd, mekooh_dd, tolooh_dd, terpooh_dd, ch3cooh_dd
  logical :: soa_dd, so4_dd, cb1_dd, cb2_dd, oc1_dd, oc2_dd, nh3_dd, nh4no3_dd, &
             sa1_dd, sa2_dd, sa3_dd, sa4_dd, nh4_dd
  logical :: soam_dd, soai_dd, soat_dd, soab_dd, soax_dd, &
             sogm_dd, sogi_dd, sogt_dd, sogb_dd, sogx_dd

  logical :: pan_dd, mpan_dd, no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
             h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
             ch3ooh_dd, pooh_dd, ch3coooh_dd, c2h5ooh_dd, eooh_dd, ch3cho_dd, c2h5oh_dd, &
             c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
             glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd

  integer :: so2_ndx
  integer :: ch3cn_ndx, hcn_ndx, hcooh_ndx
  logical :: ch3cn_dd,  hcn_dd, hcooh_dd

  integer :: o3a_ndx,xpan_ndx,xmpan_ndx,xno2_ndx,xhno3_ndx,xonit_ndx,xonitr_ndx,xno_ndx,xho2no2_ndx,xnh4no3_ndx
  logical :: o3a_dd, xpan_dd, xmpan_dd, xno2_dd, xhno3_dd, xonit_dd, xonitr_dd, xno_dd, xho2no2_dd, xnh4no3_dd

  integer :: cohc_ndx=-1, come_ndx=-1, co01_ndx=-1, co02_ndx=-1, co03_ndx=-1, co04_ndx=-1, co05_ndx=-1
  integer :: co06_ndx=-1, co07_ndx=-1, co08_ndx=-1, co09_ndx=-1, co10_ndx=-1
  integer :: co11_ndx=-1, co12_ndx=-1, co13_ndx=-1, co14_ndx=-1, co15_ndx=-1
  integer :: co16_ndx=-1, co17_ndx=-1, co18_ndx=-1, co19_ndx=-1, co20_ndx=-1
  integer :: co21_ndx=-1, co22_ndx=-1, co23_ndx=-1, co24_ndx=-1, co25_ndx=-1
  integer :: co26_ndx=-1, co27_ndx=-1, co28_ndx=-1, co29_ndx=-1, co30_ndx=-1
  integer :: co31_ndx=-1, co32_ndx=-1, co33_ndx=-1, co34_ndx=-1, co35_ndx=-1
  integer :: co36_ndx=-1, co37_ndx=-1, co38_ndx=-1, co39_ndx=-1, co40_ndx=-1
  integer :: co41_ndx=-1, co42_ndx=-1


  integer :: &
       o3_tab_ndx = -1, &
       h2o2_tab_ndx = -1, &
       ch3ooh_tab_ndx = -1, &
       co_tab_ndx = -1, &
       ch3cho_tab_ndx = -1
  logical :: &
       o3_in_tab = .false., &
       h2o2_in_tab = .false., &
       ch3ooh_in_tab = .false., &
       co_in_tab = .false., &
       ch3cho_in_tab = .false.

  real(r8), parameter    :: small_value = 1.e-36_r8
  real(r8), parameter    :: large_value = 1.e36_r8
  real(r8), parameter    :: diffm       = 1.789e-5_r8
  real(r8), parameter    :: diffk       = 1.461e-5_r8
  real(r8), parameter    :: difft       = 2.060e-5_r8
  real(r8), parameter    :: vonkar      = karman
  real(r8), parameter    :: ric         = 0.2_r8
  real(r8), parameter    :: r           = 287.04_r8
  real(r8), parameter    :: cp          = 1004._r8
  real(r8), parameter    :: grav        = 9.81_r8
  real(r8), parameter    :: p00         = 100000._r8
  real(r8), parameter    :: wh2o        = 18.0153_r8
  real(r8), parameter    :: ph          = 1.e-5_r8
  real(r8), parameter    :: ph_inv      = 1._r8/ph
  real(r8), parameter    :: rovcp = r/cp

  integer, pointer :: index_season_lai(:,:)

  logical, public :: has_dvel(gas_pcnst) = .false.
  integer         :: map_dvel(gas_pcnst) = 0
  real(r8) , allocatable            :: soilw_3d(:,:,:)

  logical, parameter :: dyn_soilw = .false.

  real(r8), allocatable  :: fraction_landuse(:,:,:)
  real(r8), allocatable, dimension(:,:,:) :: dep_ra ! [s/m] aerodynamic resistance
  real(r8), allocatable, dimension(:,:,:) :: dep_rb ! [s/m] resistance across sublayer
  integer, parameter :: n_land_type = 11

  integer, allocatable :: spc_ndx(:) ! nddvels
  real(r8), public :: crb 

  type lnd_dvel_type
     real(r8), pointer :: dvel(:,:)   ! deposition velocity over land (cm/s)
  end type lnd_dvel_type

  type(lnd_dvel_type), allocatable :: lnd(:)
  character(len=SHR_KIND_CL) :: drydep_srf_file

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  subroutine dvel_inti_fromlnd 
    use mo_chem_utls,         only : get_spc_ndx
    use abortutils,           only : endrun
    use chem_mods,            only : adv_mass
    use seq_drydep_mod,       only : dfoxd

    implicit none

    integer :: ispc, l

    allocate(spc_ndx(nddvels))
    allocate( lnd(begchunk:endchunk) )

    do ispc = 1,nddvels

       spc_ndx(ispc) = get_spc_ndx(drydep_list(ispc))
       if (spc_ndx(ispc) < 1) then
          write(*,*) 'drydep_inti: '//trim(drydep_list(ispc))//' is not included in species set'
          call endrun('drydep_init: invalid dry deposition species')
       endif

    enddo

    crb = (difft/diffm)**(2._r8/3._r8) !.666666_r8

  endsubroutine dvel_inti_fromlnd

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_update( state, cam_in )
    use physics_types,   only : physics_state
    use camsrfexch,      only : cam_in_t     
    use seq_drydep_mod,  only : drydep_method, DD_XLND

    type(physics_state), intent(in) :: state           ! Physics state variables
    type(cam_in_t),  intent(in) :: cam_in 

    if (nddvels<1) return
    if (drydep_method /= DD_XLND) return

    lnd(state%lchnk)%dvel => cam_in%depvel

  end subroutine drydep_update

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_fromlnd( ocnfrac, icefrac, ncdate, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvelocity, dflx, mmr, &
                             tv, soilw, rh, ncol, lonndx, latndx, lchnk )
                          
    !-------------------------------------------------------------------------------------
    ! combines the deposition velocities provided by the land model with deposition 
    ! velocities over ocean and sea ice 
    !-------------------------------------------------------------------------------------

    use ppgrid,         only : pcols
    use chem_mods,      only : gas_pcnst
      




    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------

    real(r8), intent(in)    :: icefrac(pcols)            
    real(r8), intent(in)    :: ocnfrac(pcols)            

    integer, intent(in)   :: ncol
    integer, intent(in)   :: ncdate                   ! present date (yyyymmdd)
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: rh(ncol,1)               ! relative humidity
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)              
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: soilw(pcols)             ! soil moisture fraction
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvelocity(ncol,gas_pcnst)    ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    integer, intent(in)     ::   latndx(pcols)           ! chunk latitude indicies
    integer, intent(in)     ::   lonndx(pcols)           ! chunk longitude indicies
    integer, intent(in)     ::   lchnk                   ! chunk number

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: ocnice_dvel(ncol,gas_pcnst)
    real(r8) :: ocnice_dflx(pcols,gas_pcnst)

    real(r8), dimension(ncol) :: term    ! work array
    integer  :: ispec
    real(r8)  :: lndfrac(pcols)            





    lndfrac(:ncol) = 1._r8 - ocnfrac(:ncol) - icefrac(:ncol)

    where( lndfrac(:ncol) < 0._r8 ) 
       lndfrac(:ncol) = 0._r8 
    endwhere





    !-------------------------------------------------------------------------------------
    !   ... initialize
    !-------------------------------------------------------------------------------------
    dvelocity(:,:) = 0._r8
    
    !-------------------------------------------------------------------------------------
    !   ... compute the dep velocities over ocean and sea ice
    !       land type 7 is used for ocean
    !       land type 8 is used for sea ice
    !-------------------------------------------------------------------------------------
    call drydep_xactive( ncdate, sfc_temp, pressure_sfc,  &
                         wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                         snow, solar_flux, ocnice_dvel, ocnice_dflx, mmr, &
                         tv, soilw, rh, ncol, lonndx, latndx, lchnk, &



                         ocnfrc=ocnfrac,icefrc=icefrac, beglandtype=7, endlandtype=8 )

    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... merge the land component with the non-land component
       !            ocn and ice already have fractions factored in
       !-------------------------------------------------------------------------------------
       dvelocity(:ncol,spc_ndx(ispec)) = lnd(lchnk)%dvel(:ncol,ispec)*lndfrac(:ncol) &
                                       + ocnice_dvel(:ncol,spc_ndx(ispec))
    enddo
    
    !-------------------------------------------------------------------------------------
    !        ... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx>0 ) then
       dvelocity(:ncol,mpan_ndx) = dvelocity(:ncol,mpan_ndx)/3._r8
    endif
    if( xmpan_ndx>0 ) then
       dvelocity(:ncol,xmpan_ndx) = dvelocity(:ncol,xmpan_ndx)/3._r8
    endif
    if( hcn_ndx>0 ) then
       dvelocity(:ncol,hcn_ndx) = ocnice_dvel(:ncol,hcn_ndx) ! should be zero over land
    endif
    if( ch3cn_ndx>0 ) then
       dvelocity(:ncol,ch3cn_ndx) = ocnice_dvel(:ncol,ch3cn_ndx) ! should be zero over land
    endif

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0 .and. ch3cooh_ndx > 0 ) then
       if( has_dvel(hcooh_ndx) ) then
          dvelocity(:ncol,hcooh_ndx) = dvelocity(:ncol,ch3cooh_ndx)
       end if
    end if
    
    !-------------------------------------------------------------------------------------
    !        ... assign CO tags to CO
    ! put this kludge in for now ...  
    !  -- should be able to set all these via the table mapping in seq_drydep_mod
    !-------------------------------------------------------------------------------------
    if( cohc_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( come_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,come_ndx) = dvelocity(:ncol,co_ndx)
    endif

    if( co01_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co01_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co02_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co02_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co03_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co03_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co04_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co04_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co05_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co05_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co06_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co06_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co07_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co07_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co08_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co08_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co09_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co09_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co10_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co10_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co11_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co11_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co12_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co12_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co13_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co13_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co14_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co14_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co15_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co15_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co16_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co16_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co17_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co17_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co18_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co18_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co19_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co19_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co20_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co20_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co21_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co21_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co22_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co22_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co23_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co23_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co24_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co24_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co25_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co25_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co26_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co26_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co27_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co27_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co28_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co28_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co29_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co29_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co30_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co30_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co31_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co31_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co32_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co32_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co33_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co33_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co34_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co34_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co35_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co35_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co36_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co36_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co37_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co37_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co38_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co38_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co39_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co39_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co40_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co40_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co41_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co41_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if( co42_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,co42_ndx) = dvelocity(:ncol,co_ndx)
    endif
    if ( o3s_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,o3s_ndx) = dvelocity(:ncol,o3_ndx)
    endif

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... compute the deposition flux
       !-------------------------------------------------------------------------------------
       dflx(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,spc_ndx(ispec)) * term(:ncol) * mmr(:ncol,plev,spc_ndx(ispec))
    end do

  end subroutine drydep_fromlnd

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
! Subprogram not used   subroutine dvel_inti_table( depvel_file )
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Initialize module, depvel arrays, and a few other variables.
! Subprogram not used     !           The depvel fields will be linearly interpolated to the correct time
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use mo_constants,  only : d2r, r2d
! Subprogram not used     use ioFileMod,     only : getfil
! Subprogram not used     use string_utils,  only : to_lower, GLC
! Subprogram not used     use mo_chem_utls,  only : get_spc_ndx
! Subprogram not used     use constituents,  only : pcnst
! Subprogram not used     use interpolate_data, only : lininterp_init, lininterp, lininterp_finish,interp_type
! Subprogram not used     use mo_constants,     only : pi
! Subprogram not used     use phys_grid, only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     character(len=*), intent(in) :: depvel_file
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Local variables
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     integer :: nlat, nlon, nmonth, ndims
! Subprogram not used     integer :: dimid_lat, dimid_lon, dimid_species, dimid_time
! Subprogram not used     integer :: dimid(4), count(4), start(4)
! Subprogram not used     integer :: m, ispecies, nchar, ierr
! Subprogram not used     real(r8)    :: scale_factor
! Subprogram not used 
! Subprogram not used     real(r8), allocatable :: dvel_lats(:), dvel_lons(:)
! Subprogram not used     real(r8), allocatable :: dvel_in(:,:,:,:)                          ! input depvel array
! Subprogram not used     character(len=50) :: units
! Subprogram not used     character(len=20), allocatable :: species_names(:)             ! names of depvel species
! Subprogram not used     logical :: found
! Subprogram not used     type(file_desc_t) :: piofile
! Subprogram not used     type(var_desc_t) :: vid, vid_dvel
! Subprogram not used 
! Subprogram not used     character(len=shr_kind_cl) :: locfn
! Subprogram not used     integer :: mm,n
! Subprogram not used     integer :: plat, plon
! Subprogram not used 
! Subprogram not used     integer :: i, c, ncols
! Subprogram not used     real(r8) :: to_lats(pcols), to_lons(pcols)
! Subprogram not used     type(interp_type) :: lon_wgts, lat_wgts
! Subprogram not used     real(r8), parameter :: zero=0._r8, twopi=2._r8*pi
! Subprogram not used 
! Subprogram not used     mm = 1
! Subprogram not used     do m = 1,pcnst
! Subprogram not used        if ( len_trim(drydep_list(m))==0 ) exit
! Subprogram not used        n = get_spc_ndx(drydep_list(m))
! Subprogram not used        if ( n < 1 ) then
! Subprogram not used           write(iulog,*) 'drydep_inti: '//drydep_list(m)//' is not included in species set'
! Subprogram not used           call endrun('drydep_init: invalid dry deposition species')
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if( masterproc ) then
! Subprogram not used        write(iulog,*) 'drydep_inti: following species have dry deposition'
! Subprogram not used        do i=1,nddvels
! Subprogram not used           if( len_trim(drydep_list(i)) > 0 ) then
! Subprogram not used              write(iulog,*) 'drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used        write(iulog,*) 'drydep_inti:'
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if ( nddvels < 1 ) return
! Subprogram not used 
! Subprogram not used     plat = get_dyn_grid_parm('plat')
! Subprogram not used     plon = get_dyn_grid_parm('plon')
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Setup species maps
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     o3a_ndx   = get_spc_ndx( 'O3A')
! Subprogram not used     xpan_ndx  = get_spc_ndx( 'XPAN')
! Subprogram not used     xmpan_ndx = get_spc_ndx( 'XMPAN')
! Subprogram not used     xno2_ndx  = get_spc_ndx( 'XNO2')
! Subprogram not used     xhno3_ndx = get_spc_ndx( 'XHNO3')
! Subprogram not used     xonit_ndx     = get_spc_ndx( 'XONIT')
! Subprogram not used     xonitr_ndx    = get_spc_ndx( 'XONITR')
! Subprogram not used     xno_ndx       = get_spc_ndx( 'XNO')
! Subprogram not used     xho2no2_ndx   = get_spc_ndx( 'XHO2NO2')
! Subprogram not used     o3a_dd   = has_drydep( 'O3A')
! Subprogram not used     xpan_dd  = has_drydep( 'XPAN')
! Subprogram not used     xmpan_dd = has_drydep( 'XMPAN')
! Subprogram not used     xno2_dd  = has_drydep( 'XNO2')
! Subprogram not used     xhno3_dd = has_drydep( 'XHNO3')
! Subprogram not used     xonit_dd     = has_drydep( 'XONIT')
! Subprogram not used     xonitr_dd    = has_drydep( 'XONITR')
! Subprogram not used     xno_dd       = has_drydep( 'XNO')
! Subprogram not used     xho2no2_dd   = has_drydep( 'XHO2NO2')
! Subprogram not used 
! Subprogram not used     pan_ndx  = get_spc_ndx( 'PAN')
! Subprogram not used     mpan_ndx = get_spc_ndx( 'MPAN')
! Subprogram not used     no2_ndx  = get_spc_ndx( 'NO2')
! Subprogram not used     hno3_ndx = get_spc_ndx( 'HNO3')
! Subprogram not used     co_ndx   = get_spc_ndx( 'CO')
! Subprogram not used     o3_ndx   = get_spc_ndx( 'O3')
! Subprogram not used     if( o3_ndx < 1 ) then
! Subprogram not used        o3_ndx = get_spc_ndx( 'OX')
! Subprogram not used     end if
! Subprogram not used     h2o2_ndx     = get_spc_ndx( 'H2O2')
! Subprogram not used     onit_ndx     = get_spc_ndx( 'ONIT')
! Subprogram not used     onitr_ndx    = get_spc_ndx( 'ONITR')
! Subprogram not used     ch4_ndx      = get_spc_ndx( 'CH4')
! Subprogram not used     ch2o_ndx     = get_spc_ndx( 'CH2O')
! Subprogram not used     ch3ooh_ndx   = get_spc_ndx( 'CH3OOH')
! Subprogram not used     ch3cho_ndx   = get_spc_ndx( 'CH3CHO')
! Subprogram not used     ch3cocho_ndx = get_spc_ndx( 'CH3COCHO')
! Subprogram not used     pooh_ndx     = get_spc_ndx( 'POOH')
! Subprogram not used     ch3coooh_ndx = get_spc_ndx( 'CH3COOOH')
! Subprogram not used     c2h5ooh_ndx  = get_spc_ndx( 'C2H5OOH')
! Subprogram not used     eooh_ndx     = get_spc_ndx( 'EOOH')
! Subprogram not used     c3h7ooh_ndx  = get_spc_ndx( 'C3H7OOH')
! Subprogram not used     rooh_ndx     = get_spc_ndx( 'ROOH')
! Subprogram not used     ch3coch3_ndx = get_spc_ndx( 'CH3COCH3')
! Subprogram not used     no_ndx       = get_spc_ndx( 'NO')
! Subprogram not used     ho2no2_ndx   = get_spc_ndx( 'HO2NO2')
! Subprogram not used     glyald_ndx   = get_spc_ndx( 'GLYALD')
! Subprogram not used     hyac_ndx     = get_spc_ndx( 'HYAC')
! Subprogram not used     ch3oh_ndx    = get_spc_ndx( 'CH3OH')
! Subprogram not used     c2h5oh_ndx   = get_spc_ndx( 'C2H5OH')
! Subprogram not used     macrooh_ndx  = get_spc_ndx( 'MACROOH')
! Subprogram not used     isopooh_ndx  = get_spc_ndx( 'ISOPOOH')
! Subprogram not used     xooh_ndx     = get_spc_ndx( 'XOOH')
! Subprogram not used     hydrald_ndx  = get_spc_ndx( 'HYDRALD')
! Subprogram not used     h2_ndx       = get_spc_ndx( 'H2')
! Subprogram not used     Pb_ndx       = get_spc_ndx( 'Pb')
! Subprogram not used     o3s_ndx      = get_spc_ndx( 'O3S')
! Subprogram not used     o3inert_ndx  = get_spc_ndx( 'O3INERT')
! Subprogram not used     alkooh_ndx  = get_spc_ndx( 'ALKOOH')
! Subprogram not used     mekooh_ndx  = get_spc_ndx( 'MEKOOH')
! Subprogram not used     tolooh_ndx  = get_spc_ndx( 'TOLOOH')
! Subprogram not used     terpooh_ndx = get_spc_ndx( 'TERPOOH')
! Subprogram not used     ch3cooh_ndx = get_spc_ndx( 'CH3COOH')
! Subprogram not used     soam_ndx    = get_spc_ndx( 'SOAM' )
! Subprogram not used     soai_ndx    = get_spc_ndx( 'SOAI' )
! Subprogram not used     soat_ndx    = get_spc_ndx( 'SOAT' )
! Subprogram not used     soab_ndx    = get_spc_ndx( 'SOAB' )
! Subprogram not used     soax_ndx    = get_spc_ndx( 'SOAX' )
! Subprogram not used     sogm_ndx    = get_spc_ndx( 'SOGM' )
! Subprogram not used     sogi_ndx    = get_spc_ndx( 'SOGI' )
! Subprogram not used     sogt_ndx    = get_spc_ndx( 'SOGT' )
! Subprogram not used     sogb_ndx    = get_spc_ndx( 'SOGB' )
! Subprogram not used     sogx_ndx    = get_spc_ndx( 'SOGX' )
! Subprogram not used     soa_ndx     = get_spc_ndx( 'SOA' )
! Subprogram not used     so4_ndx     = get_spc_ndx( 'SO4' )
! Subprogram not used     cb1_ndx     = get_spc_ndx( 'CB1' )
! Subprogram not used     cb2_ndx     = get_spc_ndx( 'CB2' )
! Subprogram not used     oc1_ndx     = get_spc_ndx( 'OC1' )
! Subprogram not used     oc2_ndx     = get_spc_ndx( 'OC2' )
! Subprogram not used     nh3_ndx     = get_spc_ndx( 'NH3' )
! Subprogram not used     nh4no3_ndx  = get_spc_ndx( 'NH4NO3' )
! Subprogram not used     xnh4no3_ndx  = get_spc_ndx( 'XNH4NO3' )
! Subprogram not used     sa1_ndx     = get_spc_ndx( 'SA1' )
! Subprogram not used     sa2_ndx     = get_spc_ndx( 'SA2' )
! Subprogram not used     sa3_ndx     = get_spc_ndx( 'SA3' )
! Subprogram not used     sa4_ndx     = get_spc_ndx( 'SA4' )
! Subprogram not used     nh4_ndx     = get_spc_ndx( 'NH4' )
! Subprogram not used     alkooh_dd  = has_drydep( 'ALKOOH')
! Subprogram not used     mekooh_dd  = has_drydep( 'MEKOOH')
! Subprogram not used     tolooh_dd  = has_drydep( 'TOLOOH')
! Subprogram not used     terpooh_dd = has_drydep( 'TERPOOH')
! Subprogram not used     ch3cooh_dd = has_drydep( 'CH3COOH')
! Subprogram not used     soam_dd    = has_drydep( 'SOAM' )
! Subprogram not used     soai_dd    = has_drydep( 'SOAI' )
! Subprogram not used     soat_dd    = has_drydep( 'SOAT' )
! Subprogram not used     soab_dd    = has_drydep( 'SOAB' )
! Subprogram not used     soax_dd    = has_drydep( 'SOAX' )
! Subprogram not used     sogm_dd    = has_drydep( 'SOGM' )
! Subprogram not used     sogi_dd    = has_drydep( 'SOGI' )
! Subprogram not used     sogt_dd    = has_drydep( 'SOGT' )
! Subprogram not used     sogb_dd    = has_drydep( 'SOGB' )
! Subprogram not used     sogx_dd    = has_drydep( 'SOGX' )
! Subprogram not used     soa_dd     = has_drydep( 'SOA' )
! Subprogram not used     so4_dd     = has_drydep( 'SO4' )
! Subprogram not used     cb1_dd     = has_drydep( 'CB1' )
! Subprogram not used     cb2_dd     = has_drydep( 'CB2' )
! Subprogram not used     oc1_dd     = has_drydep( 'OC1' )
! Subprogram not used     oc2_dd     = has_drydep( 'OC2' )
! Subprogram not used     nh3_dd     = has_drydep( 'NH3' )
! Subprogram not used     nh4no3_dd  = has_drydep( 'NH4NO3' )
! Subprogram not used     xnh4no3_dd = has_drydep( 'XNH4NO3' )
! Subprogram not used     sa1_dd     = has_drydep( 'SA1' ) 
! Subprogram not used     sa2_dd     = has_drydep( 'SA2' )
! Subprogram not used     sa3_dd     = has_drydep( 'SA3' ) 
! Subprogram not used     sa4_dd     = has_drydep( 'SA4' )
! Subprogram not used     nh4_dd     = has_drydep( 'NH4' ) 
! Subprogram not used     pan_dd  = has_drydep( 'PAN')
! Subprogram not used     mpan_dd = has_drydep( 'MPAN')
! Subprogram not used     no2_dd  = has_drydep( 'NO2')
! Subprogram not used     hno3_dd = has_drydep( 'HNO3')
! Subprogram not used     co_dd   = has_drydep( 'CO')
! Subprogram not used     o3_dd   = has_drydep( 'O3')
! Subprogram not used     if( .not. o3_dd ) then
! Subprogram not used        o3_dd = has_drydep( 'OX')
! Subprogram not used     end if
! Subprogram not used     h2o2_dd     = has_drydep( 'H2O2')
! Subprogram not used     onit_dd     = has_drydep( 'ONIT')
! Subprogram not used     onitr_dd    = has_drydep( 'ONITR')
! Subprogram not used     ch4_dd      = has_drydep( 'CH4')
! Subprogram not used     ch2o_dd     = has_drydep( 'CH2O')
! Subprogram not used     ch3ooh_dd   = has_drydep( 'CH3OOH')
! Subprogram not used     ch3cho_dd   = has_drydep( 'CH3CHO')
! Subprogram not used     c2h5oh_dd   = has_drydep( 'C2H5OH')
! Subprogram not used     eooh_dd     = has_drydep( 'EOOH')
! Subprogram not used     ch3cocho_dd = has_drydep( 'CH3COCHO')
! Subprogram not used     pooh_dd     = has_drydep( 'POOH')
! Subprogram not used     ch3coooh_dd = has_drydep( 'CH3COOOH')
! Subprogram not used     c2h5ooh_dd  = has_drydep( 'C2H5OOH')
! Subprogram not used     c3h7ooh_dd  = has_drydep( 'C3H7OOH')
! Subprogram not used     rooh_dd     = has_drydep( 'ROOH')
! Subprogram not used     ch3coch3_dd = has_drydep( 'CH3COCH3')
! Subprogram not used     glyald_dd   = has_drydep( 'GLYALD')
! Subprogram not used     hyac_dd     = has_drydep( 'HYAC')
! Subprogram not used     ch3oh_dd    = has_drydep( 'CH3OH')
! Subprogram not used     macrooh_dd  = has_drydep( 'MACROOH')
! Subprogram not used     isopooh_dd  = has_drydep( 'ISOPOOH')
! Subprogram not used     xooh_dd     = has_drydep( 'XOOH')
! Subprogram not used     hydrald_dd  = has_drydep( 'HYDRALD')
! Subprogram not used     h2_dd       = has_drydep( 'H2')
! Subprogram not used     Pb_dd       = has_drydep( 'Pb')
! Subprogram not used     o3s_dd      = has_drydep( 'O3S')
! Subprogram not used     o3inert_dd  = has_drydep( 'O3INERT')
! Subprogram not used     ch3cn_dd    = has_drydep( 'CH3CN')
! Subprogram not used     hcn_dd      = has_drydep( 'HCN')
! Subprogram not used     hcooh_dd    = has_drydep( 'HCOOH')
! Subprogram not used     ch3cn_ndx   = get_spc_ndx( 'CH3CN')
! Subprogram not used     hcn_ndx     = get_spc_ndx( 'HCN')
! Subprogram not used     hcooh_ndx   = get_spc_ndx( 'HCOOH' )
! Subprogram not used 
! Subprogram not used     if( masterproc ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: diagnostics'
! Subprogram not used        write(iulog,'(10i5)') pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
! Subprogram not used             h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
! Subprogram not used             ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, eooh_ndx, &
! Subprogram not used             c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
! Subprogram not used             no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
! Subprogram not used             hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
! Subprogram not used             xooh_ndx, ch3cho_ndx, isopooh_ndx
! Subprogram not used        write(iulog,*) pan_dd, mpan_dd, no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
! Subprogram not used             h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
! Subprogram not used             ch3ooh_dd, pooh_dd, ch3coooh_dd, c2h5ooh_dd, eooh_dd, ch3cho_dd, c2h5oh_dd, &
! Subprogram not used             c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
! Subprogram not used             glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd
! Subprogram not used     endif
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Open NetCDF file
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     call getfil (depvel_file, locfn, 0)
! Subprogram not used     call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Get variable ID for dep vel array
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_varid( piofile, 'dvel', vid_dvel )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Inquire about dimensions
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'species', dimid_species )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_species, nspecies )
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'time', dimid_time )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_time, nmonth )
! Subprogram not used     if(masterproc) write(iulog,*) 'dvel_inti: dimensions (nlon,nlat,nspecies,nmonth) = ',nlon,nlat,nspecies,nmonth
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Check dimensions of dvel variable. Must be (lon, lat, species, month).
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_varndims( piofile, vid_dvel, ndims )
! Subprogram not used 
! Subprogram not used     if( masterproc .and. ndims /= 4 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: dvel has ',ndims,' dimensions. Expecting 4.'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     ierr = pio_inq_vardimid( piofile, vid_dvel, dimid )
! Subprogram not used 
! Subprogram not used     if( dimid(1) /= dimid_lon .or. dimid(2) /= dimid_lat .or. &
! Subprogram not used          dimid(3) /= dimid_species .or. dimid(4) /= dimid_time ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: Dimensions in wrong order for dvel'
! Subprogram not used        write(iulog,*) '...      Expecting (lon, lat, species, month)'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Allocate depvel lats, lons and read
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     allocate( dvel_lats(nlat), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: Failed to allocate dvel_lats vector'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     allocate( dvel_lons(nlon), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: Failed to allocate dvel_lons vector'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_varid( piofile, 'lat', vid )
! Subprogram not used     ierr = pio_get_var( piofile, vid, dvel_lats )
! Subprogram not used     ierr = pio_inq_varid( piofile, 'lon', vid )
! Subprogram not used     ierr = pio_get_var( piofile, vid, dvel_lons )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !       ... Set the transform from inputs lats to simulation lats
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     dvel_lats(:nlat) = d2r * dvel_lats(:nlat)
! Subprogram not used     dvel_lons(:nlon) = d2r * dvel_lons(:nlon)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !     	... Allocate dvel and read data from file
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     allocate( dvel_in(nlon, nlat ,nspecies, nmonth), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: Failed to allocate dvel_in'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     start = (/ 1, 1, 1, 1 /)
! Subprogram not used     count = (/ nlon, nlat, nspecies, nmonth /)
! Subprogram not used 
! Subprogram not used     ierr = pio_get_var( piofile, vid_dvel, start, count, dvel_in )
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !     	... Check units of deposition velocity. If necessary, convert to cm/s.
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     units(:) = ' '
! Subprogram not used     ierr = pio_get_att( piofile, vid_dvel, 'units', units )
! Subprogram not used     if( to_lower(trim(units(:GLC(units)))) == 'm/s' ) then
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        scale_factor = 100._r8
! Subprogram not used     elseif( to_lower(trim(units(:GLC(units)))) == 'cm/s' ) then
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        scale_factor = 1._r8
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        scale_factor = 1._r8
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     dvel_in(:,:,:,:) = scale_factor*dvel_in(:,:,:,:)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !     	... Regrid deposition velocities
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     allocate( dvel(pcols,begchunk:endchunk,nspecies,nmonth),stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: Failed to allocate dvel'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     do c=begchunk,endchunk
! Subprogram not used        ncols = get_ncols_p(c)
! Subprogram not used        call get_rlat_all_p(c, pcols, to_lats)
! Subprogram not used        call get_rlon_all_p(c, pcols, to_lons)
! Subprogram not used        call lininterp_init(dvel_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
! Subprogram not used        call lininterp_init(dvel_lats, nlat, to_lats, ncols, 1, lat_wgts)
! Subprogram not used 
! Subprogram not used        do ispecies = 1,nspecies
! Subprogram not used           do m = 1,12
! Subprogram not used              call lininterp( dvel_in( :,:,ispecies,m ), nlon, nlat, dvel(:,c,ispecies,m), ncols,lon_wgts,lat_wgts)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        call lininterp_finish(lat_wgts)
! Subprogram not used        call lininterp_finish(lon_wgts)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     deallocate( dvel_in )
! Subprogram not used     deallocate( dvel_lats, dvel_lons )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !     	... Read in species names and determine mapping to tracer numbers
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     allocate( species_names(nspecies), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: species_names allocation error = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     ierr = pio_inq_varid( piofile, 'species_name', vid )
! Subprogram not used     ierr = pio_inq_varndims( piofile, vid, ndims )
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_vardimid( piofile, vid, dimid )
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid(1), nchar )
! Subprogram not used     map(:) = 0
! Subprogram not used     do ispecies = 1,nspecies
! Subprogram not used        start(:2) = (/ 1, ispecies /)
! Subprogram not used        count(:2) = (/ nchar, 1 /)
! Subprogram not used        species_names(ispecies)(:) = ' '
! Subprogram not used        ierr = pio_get_var( piofile, vid, start(1:2), count(1:2), species_names(ispecies:ispecies) )
! Subprogram not used        if( species_names(ispecies) == 'O3' ) then
! Subprogram not used           o3_in_tab  = .true.
! Subprogram not used           o3_tab_ndx = ispecies
! Subprogram not used        else if( species_names(ispecies) == 'H2O2' ) then
! Subprogram not used           h2o2_in_tab  = .true.
! Subprogram not used           h2o2_tab_ndx = ispecies
! Subprogram not used        else if( species_names(ispecies) == 'CH3OOH' ) then
! Subprogram not used           ch3ooh_in_tab  = .true.
! Subprogram not used           ch3ooh_tab_ndx = ispecies
! Subprogram not used        else if( species_names(ispecies) == 'CO' ) then
! Subprogram not used           co_in_tab  = .true.
! Subprogram not used           co_tab_ndx = ispecies
! Subprogram not used        else if( species_names(ispecies) == 'CH3CHO' ) then
! Subprogram not used           ch3cho_in_tab  = .true.
! Subprogram not used           ch3cho_tab_ndx = ispecies
! Subprogram not used        end if
! Subprogram not used        found = .false.
! Subprogram not used        do m = 1,gas_pcnst
! Subprogram not used           if( species_names(ispecies) == solsym(m) .or. &
! Subprogram not used                (species_names(ispecies) == 'O3' .and. solsym(m) == 'OX') .or. &
! Subprogram not used                (species_names(ispecies) == 'HNO4' .and. solsym(m) == 'HO2NO2') ) then
! Subprogram not used              if ( has_drydep( solsym(m) ) ) then
! Subprogram not used                 map(m) = ispecies
! Subprogram not used                 found = .true.
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used                 exit
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if( .not. found ) then
! Subprogram not used           write(iulog,*) 'dvel_inti: Warning! DVEL species ',trim(species_names(ispecies)),' not found'
! Subprogram not used        endif
! Subprogram not used     end do
! Subprogram not used     deallocate( species_names )
! Subprogram not used 
! Subprogram not used     call pio_closefile( piofile )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !     	... Allocate dvel_interp array
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     allocate( dvel_interp(pcols,begchunk:endchunk,nspecies),stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'dvel_inti: Failed to allocate dvel_interp; error = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine dvel_inti_table

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine interpdvel( calday, ncol, lchnk )
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     ! 	... Interpolate the fields whose values are required at the
! Subprogram not used     !           begining of a timestep.
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use time_manager,  only : get_calday
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     ! 	... Dummy arguments
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     real(r8), intent(in) :: calday   ! Interpolate the input data to calday
! Subprogram not used     integer, intent(in) :: ncol, lchnk
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     ! 	... Local variables
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     integer :: m, last, next
! Subprogram not used     integer  ::  dates(12) = (/ 116, 214, 316, 415,  516,  615, &
! Subprogram not used                                 716, 816, 915, 1016, 1115, 1216 /)
! Subprogram not used     real(r8) :: calday_loc, last_days, next_days
! Subprogram not used     real(r8), save ::  dys(12)
! Subprogram not used     logical, save  ::  entered = .false.
! Subprogram not used 
! Subprogram not used     if( .not. entered ) then
! Subprogram not used        do m = 1,12
! Subprogram not used           dys(m) = get_calday( dates(m), 0 )
! Subprogram not used        end do
! Subprogram not used        entered = .true.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if( calday < dys(1) ) then
! Subprogram not used        next = 1
! Subprogram not used        last = 12
! Subprogram not used     else if( calday >= dys(12) ) then
! Subprogram not used        next = 1
! Subprogram not used        last = 12
! Subprogram not used     else
! Subprogram not used        do m = 11,1,-1
! Subprogram not used           if( calday >= dys(m) ) then
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        last = m
! Subprogram not used        next = m + 1
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     last_days  = dys( last )
! Subprogram not used     next_days  = dys( next )
! Subprogram not used     calday_loc = calday
! Subprogram not used 
! Subprogram not used     if( next_days < last_days ) then
! Subprogram not used        next_days = next_days + 365._r8
! Subprogram not used     end if
! Subprogram not used     if( calday_loc < last_days ) then
! Subprogram not used        calday_loc = calday_loc + 365._r8
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     do m = 1,nspecies
! Subprogram not used        call intp2d( last_days, next_days, calday_loc, ncol, lchnk, &
! Subprogram not used                     dvel(:,lchnk,m,last), &
! Subprogram not used                     dvel(:,lchnk,m,next), &
! Subprogram not used                     dvel_interp(:,lchnk,m) )
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine interpdvel

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine intp2d( t1, t2, tint, ncol, lchnk, f1, f2, fint )
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... Linearly interpolate between f1(t1) and f2(t2) to fint(tint).
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... Dummy arguments
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     real(r8), intent(in) :: &
! Subprogram not used          t1, &            ! time level of f1
! Subprogram not used          t2, &            ! time level of f2
! Subprogram not used          tint             ! interpolant time
! Subprogram not used     real(r8), dimension(pcols), intent(in) :: &
! Subprogram not used          f1, &            ! field at time t1
! Subprogram not used          f2               ! field at time t2
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: ncol, lchnk
! Subprogram not used 
! Subprogram not used     real(r8), intent(out) :: &
! Subprogram not used          fint(pcols) ! field at time tint
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... Local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer  :: j, plat
! Subprogram not used     real(r8) :: factor
! Subprogram not used 
! Subprogram not used     factor = (tint - t1)/(t2 - t1)
! Subprogram not used 
! Subprogram not used     fint(:ncol) = f1(:ncol) + (f2(:ncol) - f1(:ncol))*factor
! Subprogram not used 
! Subprogram not used   end subroutine intp2d

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine drydep_table( calday, tsurf, zen_angle, &
! Subprogram not used                            depvel, dflx, q, p, &
! Subprogram not used                            tv, ncol, icefrac, ocnfrac, lchnk )
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Form the deposition velocities for this
! Subprogram not used     !           latitude slice
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use physconst,     only : rair,pi
! Subprogram not used     use dycore,        only : dycore_is
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Dummy arguments
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     integer, intent(in)     ::   ncol                    ! columns in chunk
! Subprogram not used     real(r8), intent(in)    ::  q(pcols,plev,gas_pcnst) ! tracer mmr (kg/kg)
! Subprogram not used     real(r8), intent(in)    ::  p(pcols)            ! midpoint pressure in surface layer (Pa)
! Subprogram not used     real(r8), intent(in)    ::  tv(pcols)           ! virtual temperature in surface layer (K)
! Subprogram not used     real(r8), intent(in)    ::  calday              ! time of year in days
! Subprogram not used     real(r8), intent(in)    ::  tsurf(pcols)        ! surface temperature (K)
! Subprogram not used     real(r8), intent(in)    ::  zen_angle(ncol)    ! zenith angle (radians)
! Subprogram not used     real(r8), intent(inout) ::  dflx(pcols,gas_pcnst)   ! flux due to dry deposition (kg/m^2/sec)
! Subprogram not used     real(r8), intent(out)   ::  depvel(ncol,gas_pcnst) ! deposition vel (cm/s)
! Subprogram not used 
! Subprogram not used     real(r8), intent(in) :: icefrac(pcols)          ! sea-ice areal fraction
! Subprogram not used     real(r8), intent(in) :: ocnfrac(pcols)          ! ocean areal fraction
! Subprogram not used     
! Subprogram not used     integer, intent(in)     :: lchnk
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... Local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer :: m, spc_ndx, tmp_ndx, i
! Subprogram not used     real(r8), dimension(ncol) :: vel, glace, temp_fac, wrk, tmp
! Subprogram not used     real(r8), dimension(ncol) :: o3_tab_dvel
! Subprogram not used     real(r8), dimension(ncol) :: ocean 
! Subprogram not used 
! Subprogram not used     real(r8), parameter :: pid2 = .5_r8 * pi
! Subprogram not used 
! Subprogram not used     if(dycore_is('UNSTRUCTURED')) then
! Subprogram not used        call endrun( 'Option not supported for unstructured atmosphere grids ')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... Note the factor 1.e-2 in the wrk array calculation is
! Subprogram not used     !           to transform the incoming dep vel from cm/s to m/s
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     wrk(:ncol) =  1.e-2_r8 * p(:ncol) / (rair * tv(:ncol))
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Initialize all deposition velocities to zero
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     depvel(:,:) = 0._r8
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Time interpolate primary depvel array
! Subprogram not used     !           (also seaice and npp)
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     call interpdvel( calday, ncol, lchnk )
! Subprogram not used 
! Subprogram not used     if( o3_in_tab ) then
! Subprogram not used        do i=1,ncol
! Subprogram not used           o3_tab_dvel(i) = dvel_interp(i,lchnk,o3_tab_ndx)
! Subprogram not used        enddo
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set deposition velocities
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     do m = 1,gas_pcnst
! Subprogram not used        if( map(m) /= 0 ) then
! Subprogram not used           do i = 1,ncol
! Subprogram not used              depvel(i,m) = dvel_interp(i,lchnk,map(m))
! Subprogram not used              dflx(i,m)   = wrk(i) * depvel(i,m) * q(i,plev,m)
! Subprogram not used           enddo
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set some variables needed for some dvel calculations
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     temp_fac(:ncol)   = min( 1._r8, max( 0._r8, (tsurf(:ncol) - 268._r8) / 5._r8 ) )
! Subprogram not used     ocean(:ncol)  = icefrac(:ncol)+ocnfrac(:ncol)
! Subprogram not used     glace(:ncol)  = icefrac(:ncol) + (1._r8 - ocean(:ncol)) * (1._r8 - temp_fac(:ncol))
! Subprogram not used     glace(:ncol)  = min( 1._r8,glace(:ncol) )
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set pan & mpan
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( o3_in_tab ) then
! Subprogram not used        tmp(:ncol) = o3_tab_dvel(:ncol) / 3._r8
! Subprogram not used     else
! Subprogram not used        tmp(:) = 0._r8
! Subprogram not used     end if
! Subprogram not used     if( pan_dd ) then
! Subprogram not used        if( map(pan_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,pan_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,pan_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,pan_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( mpan_dd ) then
! Subprogram not used        if( map(mpan_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,mpan_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,mpan_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,mpan_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set no2 dvel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( no2_dd ) then
! Subprogram not used        if( map(no2_ndx) == 0 .and. o3_in_tab ) then
! Subprogram not used           depvel(:ncol,no2_ndx) = (.6_r8*o3_tab_dvel(:ncol) + .055_r8*ocean(:ncol)) * .9_r8
! Subprogram not used           dflx(:ncol,no2_ndx)   = wrk(:) * depvel(:ncol,no2_ndx) * q(:ncol,plev,no2_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set hno3 dvel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     tmp(:ncol) = (2._r8 - ocnfrac(:ncol)) * (1._r8 - glace(:ncol)) + .05_r8 * glace(:ncol)
! Subprogram not used     if( hno3_dd ) then
! Subprogram not used        if( map(hno3_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,hno3_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,hno3_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,hno3_ndx)
! Subprogram not used        else
! Subprogram not used           tmp(:ncol) = depvel(:ncol,hno3_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( onitr_dd ) then
! Subprogram not used        if( map(onitr_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,onitr_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,onitr_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,onitr_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( isopooh_dd ) then
! Subprogram not used        if( map(isopooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,isopooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,isopooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,isopooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set h2o2 dvel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( .not. h2o2_in_tab ) then
! Subprogram not used        if( o3_in_tab ) then
! Subprogram not used           tmp(:ncol) = .05_r8*glace(:ncol) + ocean(:ncol) - icefrac(:ncol) &
! Subprogram not used                + (1._r8 - (glace(:) + ocean(:ncol)) + icefrac(:ncol)) &
! Subprogram not used                *max( 1._r8,1._r8/(.5_r8 + 1._r8/(6._r8*o3_tab_dvel(:ncol))) )
! Subprogram not used        else
! Subprogram not used           tmp(:ncol) = 0._r8
! Subprogram not used        end if
! Subprogram not used     else
! Subprogram not used        do i=1,ncol
! Subprogram not used           tmp(i) = dvel_interp(i,lchnk,h2o2_tab_ndx)
! Subprogram not used        enddo
! Subprogram not used     end if
! Subprogram not used     if( h2o2_dd ) then
! Subprogram not used        if( map(h2o2_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,h2o2_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,h2o2_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,h2o2_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set hcn dvel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( hcn_dd ) then
! Subprogram not used        if( map(hcn_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,hcn_ndx) = ocnfrac(:ncol)*0.2_r8
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set ch3cn dvel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( ch3cn_dd ) then
! Subprogram not used        if( map(ch3cn_ndx) == 0 ) then
! Subprogram not used           depvel(:,ch3cn_ndx) = ocnfrac(:ncol)*0.2_r8
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set onit
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( onit_dd ) then
! Subprogram not used        if( map(onit_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,onit_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,onit_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,onit_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( ch3cocho_dd ) then
! Subprogram not used        if( map(ch3cocho_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch3cocho_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,ch3cocho_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3cocho_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( ch3ooh_in_tab ) then
! Subprogram not used        do i=1,ncol
! Subprogram not used           tmp(i) = dvel_interp(i,lchnk,ch3ooh_tab_ndx)
! Subprogram not used        enddo
! Subprogram not used     else
! Subprogram not used        tmp(:ncol) = .5_r8 * tmp(:ncol)
! Subprogram not used     end if
! Subprogram not used     if( ch3ooh_dd ) then
! Subprogram not used        if( map(ch3ooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch3ooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,ch3ooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3ooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( pooh_dd ) then
! Subprogram not used        if( map(pooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,pooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,pooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,pooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( ch3coooh_dd ) then
! Subprogram not used        if( map(ch3coooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch3coooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,ch3coooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3coooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( c2h5ooh_dd ) then
! Subprogram not used        if( map(c2h5ooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,c2h5ooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,c2h5ooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,c2h5ooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( c3h7ooh_dd ) then
! Subprogram not used        if( map(c3h7ooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,c3h7ooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,c3h7ooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,c3h7ooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( rooh_dd ) then
! Subprogram not used        if( map(rooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,rooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,rooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,rooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( macrooh_dd ) then
! Subprogram not used        if( map(macrooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,macrooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,macrooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,macrooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( xooh_dd ) then
! Subprogram not used        if( map(xooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,xooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,xooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( ch3oh_dd ) then
! Subprogram not used        if( map(ch3oh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch3oh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,ch3oh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3oh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( c2h5oh_dd ) then
! Subprogram not used        if( map(c2h5oh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,c2h5oh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,c2h5oh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,c2h5oh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( alkooh_dd ) then
! Subprogram not used        if( map(alkooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,alkooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,alkooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,alkooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( mekooh_dd ) then
! Subprogram not used        if( map(mekooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,mekooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,mekooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,mekooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( tolooh_dd ) then
! Subprogram not used        if( map(tolooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,tolooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,tolooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,tolooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( terpooh_dd ) then
! Subprogram not used        if( map(terpooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,terpooh_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,terpooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,terpooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if( o3_in_tab ) then
! Subprogram not used        tmp(:ncol) = o3_tab_dvel(:ncol)
! Subprogram not used     else
! Subprogram not used        tmp(:ncol) = 0._r8
! Subprogram not used     end if
! Subprogram not used     if( ch2o_dd ) then
! Subprogram not used        if( map(ch2o_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch2o_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,ch2o_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch2o_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if( hydrald_dd ) then
! Subprogram not used        if( map(hydrald_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,hydrald_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,hydrald_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,hydrald_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( ch3cooh_dd  ) then
! Subprogram not used        if( map(ch3cooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch3cooh_ndx) = depvel(:ncol,ch2o_ndx)
! Subprogram not used           dflx(:ncol,ch3cooh_ndx) = wrk(:ncol) * depvel(:ncol,ch3cooh_ndx) * q(:ncol,plev,ch3cooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( eooh_dd ) then
! Subprogram not used        if( map(eooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,eooh_ndx) = depvel(:ncol,ch2o_ndx)
! Subprogram not used           dflx(:ncol,eooh_ndx) = wrk(:ncol) * depvel(:ncol,eooh_ndx) * q(:ncol,plev,eooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     ! HCOOH - set to CH3COOH
! Subprogram not used     if( hcooh_dd  ) then
! Subprogram not used        if( map(hcooh_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,hcooh_ndx) = depvel(:ncol,ch2o_ndx)
! Subprogram not used           dflx(:ncol,hcooh_ndx) = wrk(:ncol) * depvel(:ncol,hcooh_ndx) * q(:ncol,plev,hcooh_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set co and related species dep vel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( co_in_tab ) then
! Subprogram not used        do i=1,ncol
! Subprogram not used           tmp(i) = dvel_interp(i,lchnk,co_tab_ndx)
! Subprogram not used        enddo
! Subprogram not used     else
! Subprogram not used        tmp(:) = 0._r8
! Subprogram not used     end if
! Subprogram not used     if( co_dd ) then
! Subprogram not used        if( map(co_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,co_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,co_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,co_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( ch3coch3_dd ) then
! Subprogram not used        if( map(ch3coch3_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,ch3coch3_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,ch3coch3_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3coch3_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( hyac_dd ) then
! Subprogram not used        if( map(hyac_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,hyac_ndx) = tmp(:ncol)
! Subprogram not used           dflx(:ncol,hyac_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,hyac_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if( h2_dd ) then
! Subprogram not used        if( map(h2_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,h2_ndx) = tmp(:ncol) * 1.5_r8                ! Hough(1991)
! Subprogram not used           dflx(:ncol,h2_ndx)   = wrk(:ncol) * depvel(:ncol,h2_ndx) * q(:ncol,plev,h2_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Set glyald
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( glyald_dd ) then
! Subprogram not used        if( map(glyald_ndx) == 0 ) then
! Subprogram not used           if( ch3cho_dd ) then
! Subprogram not used              depvel(:ncol,glyald_ndx) = depvel(:ncol,ch3cho_ndx)
! Subprogram not used           else if( ch3cho_in_tab ) then
! Subprogram not used              do i=1,ncol
! Subprogram not used                 depvel(i,glyald_ndx) = dvel_interp(i,lchnk,ch3cho_tab_ndx)
! Subprogram not used              enddo
! Subprogram not used           else
! Subprogram not used              depvel(:ncol,glyald_ndx) = 0._r8
! Subprogram not used           end if
! Subprogram not used           dflx(:ncol,glyald_ndx)   = wrk(:ncol) * depvel(:ncol,glyald_ndx) * q(:ncol,plev,glyald_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... Lead deposition
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( Pb_dd ) then
! Subprogram not used        if( map(Pb_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,Pb_ndx) = ocean(:ncol)  * .05_r8 + (1._r8 - ocean(:ncol)) * .2_r8
! Subprogram not used           dflx(:ncol,Pb_ndx)   = wrk(:ncol) * depvel(:ncol,Pb_ndx) * q(:ncol,plev,Pb_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !       ... diurnal dependence for OX dvel
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( o3_dd .or. o3s_dd .or. o3inert_dd ) then
! Subprogram not used        if( o3_dd .or. o3_in_tab ) then
! Subprogram not used           if( o3_dd ) then
! Subprogram not used              tmp(:ncol) = max( 1._r8,sqrt( (depvel(:ncol,o3_ndx) - .2_r8)**3/.27_r8 + 4._r8*depvel(:ncol,o3_ndx) + .67_r8 ) )
! Subprogram not used              vel(:ncol) = depvel(:ncol,o3_ndx)
! Subprogram not used           else if( o3_in_tab ) then
! Subprogram not used              tmp(:ncol) = max( 1._r8,sqrt( (o3_tab_dvel(:ncol) - .2_r8)**3/.27_r8 + 4._r8*o3_tab_dvel(:ncol) + .67_r8 ) )
! Subprogram not used              vel(:ncol) = o3_tab_dvel(:ncol)
! Subprogram not used           end if
! Subprogram not used           where( abs( zen_angle(:) ) > pid2 )
! Subprogram not used              vel(:) = vel(:) / tmp(:)
! Subprogram not used           elsewhere
! Subprogram not used              vel(:) = vel(:) * tmp(:)
! Subprogram not used           endwhere
! Subprogram not used 
! Subprogram not used        else
! Subprogram not used           vel(:ncol) = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( o3_dd ) then
! Subprogram not used           depvel(:ncol,o3_ndx) = vel(:ncol)
! Subprogram not used           dflx(:ncol,o3_ndx)   = wrk(:ncol) * vel(:ncol) * q(:ncol,plev,o3_ndx)
! Subprogram not used        end if
! Subprogram not used        !--------------------------------------------------------
! Subprogram not used        !       ... Set stratospheric O3 deposition
! Subprogram not used        !--------------------------------------------------------
! Subprogram not used        if( o3s_dd ) then
! Subprogram not used           depvel(:ncol,o3s_ndx) = vel(:ncol)
! Subprogram not used           dflx(:ncol,o3s_ndx)   = wrk(:ncol) * vel(:ncol) * q(:ncol,plev,o3s_ndx)
! Subprogram not used        end if
! Subprogram not used        if( o3inert_dd ) then
! Subprogram not used           depvel(:ncol,o3inert_ndx) = vel(:ncol)
! Subprogram not used           dflx(:ncol,o3inert_ndx)   = wrk(:ncol) * vel(:ncol) * q(:ncol,plev,o3inert_ndx)
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if( xno2_dd ) then 
! Subprogram not used        if( map(xno2_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xno2_ndx) = depvel(:ncol,no2_ndx)
! Subprogram not used           dflx(:ncol,xno2_ndx)   = wrk(:ncol) * depvel(:ncol,xno2_ndx) * q(:ncol,plev,xno2_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( o3a_dd ) then 
! Subprogram not used        if( map(o3a_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,o3a_ndx) = depvel(:ncol,o3_ndx)
! Subprogram not used           dflx(:ncol,o3a_ndx)   = wrk(:ncol) * depvel(:ncol,o3a_ndx) * q(:ncol,plev,o3a_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xhno3_dd ) then 
! Subprogram not used        if( map(xhno3_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xhno3_ndx) = depvel(:ncol,hno3_ndx)
! Subprogram not used           dflx(:ncol,xhno3_ndx)   = wrk(:ncol) * depvel(:ncol,xhno3_ndx) * q(:ncol,plev,xhno3_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xnh4no3_dd ) then 
! Subprogram not used        if( map(xnh4no3_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xnh4no3_ndx) = depvel(:ncol,nh4no3_ndx)
! Subprogram not used           dflx(:ncol,xnh4no3_ndx)   = wrk(:ncol) * depvel(:ncol,xnh4no3_ndx) * q(:ncol,plev,xnh4no3_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xpan_dd ) then 
! Subprogram not used        if( map(xpan_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xpan_ndx) = depvel(:ncol,pan_ndx)
! Subprogram not used           dflx(:ncol,xpan_ndx)   = wrk(:ncol) * depvel(:ncol,xpan_ndx) * q(:ncol,plev,xpan_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xmpan_dd ) then 
! Subprogram not used        if( map(xmpan_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xmpan_ndx) = depvel(:ncol,mpan_ndx)
! Subprogram not used           dflx(:ncol,xmpan_ndx)   = wrk(:ncol) * depvel(:ncol,xmpan_ndx) * q(:ncol,plev,xmpan_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xonit_dd ) then 
! Subprogram not used        if( map(xonit_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xonit_ndx) = depvel(:ncol,onit_ndx)
! Subprogram not used           dflx(:ncol,xonit_ndx)   = wrk(:ncol) * depvel(:ncol,xonit_ndx) * q(:ncol,plev,xonit_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xonitr_dd ) then 
! Subprogram not used        if( map(xonitr_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xonitr_ndx) = depvel(:ncol,onitr_ndx)
! Subprogram not used           dflx(:ncol,xonitr_ndx)   = wrk(:ncol) * depvel(:ncol,xonitr_ndx) * q(:ncol,plev,xonitr_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xno_dd ) then 
! Subprogram not used        if( map(xno_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xno_ndx) = depvel(:ncol,no_ndx)
! Subprogram not used           dflx(:ncol,xno_ndx)   = wrk(:ncol) * depvel(:ncol,xno_ndx) * q(:ncol,plev,xno_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     if( xho2no2_dd ) then 
! Subprogram not used        if( map(xho2no2_ndx) == 0 ) then
! Subprogram not used           depvel(:ncol,xho2no2_ndx) = depvel(:ncol,ho2no2_ndx)
! Subprogram not used           dflx(:ncol,xho2no2_ndx)   = wrk(:ncol) * depvel(:ncol,xho2no2_ndx) * q(:ncol,plev,xho2no2_ndx)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine drydep_table

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine dvel_inti_xactive( depvel_lnd_file, clim_soilw_file, season_wes_file )
    !-------------------------------------------------------------------------------------
    ! 	... intialize interactive drydep
    !-------------------------------------------------------------------------------------
    use dycore,        only : dycore_is
    use mo_constants,  only : r2d
    use chem_mods,     only : adv_mass
    use mo_chem_utls,  only : get_spc_ndx
    use seq_drydep_mod,only : drydep_method, DD_XATM, DD_XLND
    use phys_control,  only : phys_getopts

    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    character(len=*), intent(in) :: depvel_lnd_file, clim_soilw_file, season_wes_file 

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer :: i, j, ii, jj, jl, ju
    integer :: nlon_veg, nlat_veg, npft_veg
    integer :: nlat_lai, npft_lai, pos_min, imin
    integer :: dimid
    integer :: m, n, l, id
    integer :: length1, astat
    integer, allocatable :: wk_lai(:,:,:)
    integer, allocatable :: index_season_lai_j(:,:)
    integer :: k, num_max, k_max
    integer :: num_seas(5)
    integer :: plon, plat
    integer :: ierr

    real(r8)              :: spc_mass
    real(r8)              :: diff_min, target_lat
    real(r8), allocatable :: vegetation_map(:,:,:)
    real(r8), pointer     :: soilw_map(:,:,:)
    real(r8), allocatable :: work(:,:)
    real(r8), allocatable :: landmask(:,:)
    real(r8), allocatable :: urban(:,:)
    real(r8), allocatable :: lake(:,:)
    real(r8), allocatable :: wetland(:,:)
    real(r8), allocatable :: lon_veg(:)
    real(r8), allocatable :: lon_veg_edge(:)
    real(r8), allocatable :: lat_veg(:)
    real(r8), allocatable :: lat_veg_edge(:)
    real(r8), allocatable :: lat_lai(:)
    real(r8), allocatable :: clat(:)
    character(len=32) :: test_name
    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid
    logical :: do_soilw

    character(len=shr_kind_cl) :: locfn
    logical :: prog_modal_aero

    ! determine if modal aerosols are active so that fraction_landuse array is initialized for modal aerosal dry dep
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    call dvel_inti_fromlnd()

    if( masterproc ) then
       write(iulog,*) 'drydep_inti: following species have dry deposition'
       do i=1,nddvels
          if( len_trim(drydep_list(i)) > 0 ) then
             write(iulog,*) 'drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
          endif
       enddo
       write(iulog,*) 'drydep_inti:'
    endif

    !-------------------------------------------------------------------------------------
    ! 	... get species indices
    !-------------------------------------------------------------------------------------
    xpan_ndx      = get_spc_ndx( 'XPAN' )
    xmpan_ndx     = get_spc_ndx( 'XMPAN' )
    o3a_ndx       = get_spc_ndx( 'O3A' )

    ch4_ndx      = get_spc_ndx( 'CH4' )
    h2_ndx       = get_spc_ndx( 'H2' )
    co_ndx       = get_spc_ndx( 'CO' )
    Pb_ndx       = get_spc_ndx( 'Pb' )
    pan_ndx      = get_spc_ndx( 'PAN' )
    mpan_ndx     = get_spc_ndx( 'MPAN' )
    o3_ndx       = get_spc_ndx( 'OX' )
    if( o3_ndx < 0 ) then
       o3_ndx  = get_spc_ndx( 'O3' )
    end if
    so2_ndx     = get_spc_ndx( 'SO2' )
    alkooh_ndx  = get_spc_ndx( 'ALKOOH')
    mekooh_ndx  = get_spc_ndx( 'MEKOOH')
    tolooh_ndx  = get_spc_ndx( 'TOLOOH')
    terpooh_ndx = get_spc_ndx( 'TERPOOH')
    ch3cooh_ndx = get_spc_ndx( 'CH3COOH')
    soa_ndx     = get_spc_ndx( 'SOA' )
    so4_ndx     = get_spc_ndx( 'SO4' )
    cb1_ndx     = get_spc_ndx( 'CB1' )
    cb2_ndx     = get_spc_ndx( 'CB2' )
    oc1_ndx     = get_spc_ndx( 'OC1' )
    oc2_ndx     = get_spc_ndx( 'OC2' )
    nh3_ndx     = get_spc_ndx( 'NH3' )
    nh4no3_ndx  = get_spc_ndx( 'NH4NO3' )
    sa1_ndx     = get_spc_ndx( 'SA1' )
    sa2_ndx     = get_spc_ndx( 'SA2' )
    sa3_ndx     = get_spc_ndx( 'SA3' )
    sa4_ndx     = get_spc_ndx( 'SA4' )
    nh4_ndx     = get_spc_ndx( 'NH4' )
    alkooh_dd  = has_drydep( 'ALKOOH')
    mekooh_dd  = has_drydep( 'MEKOOH')
    tolooh_dd  = has_drydep( 'TOLOOH')
    terpooh_dd = has_drydep( 'TERPOOH')
    ch3cooh_dd = has_drydep( 'CH3COOH')
    soa_dd     = has_drydep( 'SOA' )
    so4_dd     = has_drydep( 'SO4' )
    cb1_dd     = has_drydep( 'CB1' )
    cb2_dd     = has_drydep( 'CB2' )
    oc1_dd     = has_drydep( 'OC1' )
    oc2_dd     = has_drydep( 'OC2' )
    nh3_dd     = has_drydep( 'NH3' )
    nh4no3_dd  = has_drydep( 'NH4NO3' )
    sa1_dd     = has_drydep( 'SA1' ) 
    sa2_dd     = has_drydep( 'SA2' )
    sa3_dd     = has_drydep( 'SA3' ) 
    sa4_dd     = has_drydep( 'SA4' )
    nh4_dd     = has_drydep( 'NH4' ) 
!
    soam_ndx   = get_spc_ndx( 'SOAM' )
    soai_ndx   = get_spc_ndx( 'SOAI' )
    soat_ndx   = get_spc_ndx( 'SOAT' )
    soab_ndx   = get_spc_ndx( 'SOAB' )
    soax_ndx   = get_spc_ndx( 'SOAX' )
    sogm_ndx   = get_spc_ndx( 'SOGM' )
    sogi_ndx   = get_spc_ndx( 'SOGI' )
    sogt_ndx   = get_spc_ndx( 'SOGT' )
    sogb_ndx   = get_spc_ndx( 'SOGB' )
    sogx_ndx   = get_spc_ndx( 'SOGX' )
    soam_dd    = has_drydep ( 'SOAM' )
    soai_dd    = has_drydep ( 'SOAI' )
    soat_dd    = has_drydep ( 'SOAT' )
    soab_dd    = has_drydep ( 'SOAB' )
    soax_dd    = has_drydep ( 'SOAX' )
    sogm_dd    = has_drydep ( 'SOGM' )
    sogi_dd    = has_drydep ( 'SOGI' )
    sogt_dd    = has_drydep ( 'SOGT' )
    sogb_dd    = has_drydep ( 'SOGB' )
    sogx_dd    = has_drydep ( 'SOGX' )
!
    hcn_ndx     = get_spc_ndx( 'HCN')
    ch3cn_ndx   = get_spc_ndx( 'CH3CN')

  ! for the cotags kludge..
      cohc_ndx     = get_spc_ndx( 'COhc' )
      come_ndx     = get_spc_ndx( 'COme' )
      co01_ndx     = get_spc_ndx( 'CO01' )
      co02_ndx     = get_spc_ndx( 'CO02' )
      co03_ndx     = get_spc_ndx( 'CO03' )
      co04_ndx     = get_spc_ndx( 'CO04' )
      co05_ndx     = get_spc_ndx( 'CO05' )
      co06_ndx     = get_spc_ndx( 'CO06' )
      co07_ndx     = get_spc_ndx( 'CO07' )
      co08_ndx     = get_spc_ndx( 'CO08' )
      co09_ndx     = get_spc_ndx( 'CO09' )
      co10_ndx     = get_spc_ndx( 'CO10' )
      co11_ndx     = get_spc_ndx( 'CO11' )
      co12_ndx     = get_spc_ndx( 'CO12' )
      co13_ndx     = get_spc_ndx( 'CO13' )
      co14_ndx     = get_spc_ndx( 'CO14' )
      co15_ndx     = get_spc_ndx( 'CO15' )
      co16_ndx     = get_spc_ndx( 'CO16' )
      co17_ndx     = get_spc_ndx( 'CO17' )
      co18_ndx     = get_spc_ndx( 'CO18' )
      co19_ndx     = get_spc_ndx( 'CO19' )
      co20_ndx     = get_spc_ndx( 'CO20' )
      co21_ndx     = get_spc_ndx( 'CO21' )
      co22_ndx     = get_spc_ndx( 'CO22' )
      co23_ndx     = get_spc_ndx( 'CO23' )
      co24_ndx     = get_spc_ndx( 'CO24' )
      co25_ndx     = get_spc_ndx( 'CO25' )
      co26_ndx     = get_spc_ndx( 'CO26' )
      co27_ndx     = get_spc_ndx( 'CO27' )
      co28_ndx     = get_spc_ndx( 'CO28' )
      co29_ndx     = get_spc_ndx( 'CO29' )
      co30_ndx     = get_spc_ndx( 'CO30' )
      co31_ndx     = get_spc_ndx( 'CO31' )
      co32_ndx     = get_spc_ndx( 'CO32' )
      co33_ndx     = get_spc_ndx( 'CO33' )
      co34_ndx     = get_spc_ndx( 'CO34' )
      co35_ndx     = get_spc_ndx( 'CO35' )
      co36_ndx     = get_spc_ndx( 'CO36' )
      co37_ndx     = get_spc_ndx( 'CO37' )
      co38_ndx     = get_spc_ndx( 'CO38' )
      co39_ndx     = get_spc_ndx( 'CO39' )
      co40_ndx     = get_spc_ndx( 'CO40' )
      co41_ndx     = get_spc_ndx( 'CO41' )
      co42_ndx     = get_spc_ndx( 'CO42' )

      o3s_ndx      = get_spc_ndx( 'O3S' )

    do i=1,nddvels
       if ( mapping(i) > 0 ) then
          test_name = drydep_list(i)
          m = get_spc_ndx( test_name )
          has_dvel(m) = .true.
          map_dvel(m) = i
       endif
    enddo

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !---------------------------------------------------------------------------
    ! 	... allocate module variables
    !---------------------------------------------------------------------------
    allocate( dep_ra(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_ra; error = ',astat
       call endrun
    end if
    allocate( dep_rb(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_rb; error = ',astat
       call endrun
    end if

    if (drydep_method == DD_XLND .and. (.not.prog_modal_aero)) then
       return
    endif

    do_soilw = .not. dyn_soilw .and. (has_drydep( 'H2' ) .or. has_drydep( 'CO' ))
    allocate( fraction_landuse(pcols,n_land_type, begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate fraction_landuse; error = ',astat
       call endrun
    end if
    if(do_soilw) then
       allocate(soilw_3d(pcols,12,begchunk:endchunk),stat=astat)
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate soilw_3d error = ',astat
          call endrun
       end if
    end if

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')
    allocate( index_season_lai_j(n_land_type,12),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate index_season_lai_j; error = ',astat
       call endrun
    end if
    if(dycore_is('UNSTRUCTURED') ) then
       call get_landuse_and_soilw_from_file(do_soilw)
       allocate( index_season_lai(plon,12),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
          call endrun
       end if
    else
       allocate( index_season_lai(plat,12),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read landuse map
       !---------------------------------------------------------------------------
       call getfil (depvel_lnd_file, locfn, 0)
       call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       ierr = pio_inq_dimid( piofile, 'lon', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlon_veg )
       ierr = pio_inq_dimid( piofile, 'lat', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlat_veg )
       ierr = pio_inq_dimid( piofile, 'pft', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, npft_veg )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
            landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( lon_veg(nlon_veg), lat_veg(nlat_veg), &
            lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation lon, lat arrays; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read the vegetation map and landmask
       !---------------------------------------------------------------------------
       ierr = pio_inq_varid( piofile, 'PCT_PFT', vid )
       ierr = pio_get_var( piofile, vid, vegetation_map )

       ierr = pio_inq_varid( piofile, 'LANDMASK', vid )
       ierr = pio_get_var( piofile, vid, landmask )

       ierr = pio_inq_varid( piofile, 'PCT_URBAN', vid )
       ierr = pio_get_var( piofile, vid, urban )

       ierr = pio_inq_varid( piofile, 'PCT_LAKE', vid )
       ierr = pio_get_var( piofile, vid, lake )

       ierr = pio_inq_varid( piofile, 'PCT_WETLAND', vid )
       ierr = pio_get_var( piofile, vid, wetland )

       call pio_closefile( piofile )

       !---------------------------------------------------------------------------
       ! scale vegetation, urban, lake, and wetland to fraction
       !---------------------------------------------------------------------------
       vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
       wetland(:,:)          = .01_r8 * wetland(:,:)
       lake(:,:)             = .01_r8 * lake(:,:)
       urban(:,:)            = .01_r8 * urban(:,:)







       !---------------------------------------------------------------------------
       ! 	... define lat-lon of vegetation map (1x1)
       !---------------------------------------------------------------------------
       lat_veg(:)      = (/ (-89.5_r8 + (i-1),i=1,nlat_veg  ) /)
       lon_veg(:)      = (/ (  0.5_r8 + (i-1),i=1,nlon_veg  ) /)
       lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
       lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)
       !---------------------------------------------------------------------------
       ! 	... read soilw table if necessary
       !---------------------------------------------------------------------------

       if( do_soilw ) then
          call soilw_inti( clim_soilw_file, nlon_veg, nlat_veg, soilw_map )
       end if

       !---------------------------------------------------------------------------
       ! 	... regrid to model grid
       !---------------------------------------------------------------------------

       call interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg, lat_veg_edge, &
            lon_veg, lon_veg_edge, landmask, urban, lake, &
            wetland, vegetation_map, soilw_map, do_soilw )

       deallocate( vegetation_map, work, stat=astat )
       deallocate( lon_veg, lat_veg, lon_veg_edge, lat_veg_edge, stat=astat )
       deallocate( landmask, urban, lake, wetland, stat=astat )
       if( do_soilw ) then
          deallocate( soilw_map, stat=astat )
       end if
    endif  ! Unstructured grid

    if (drydep_method == DD_XLND) then
       return
    endif

    !---------------------------------------------------------------------------
    ! 	... read LAI based season indeces
    !---------------------------------------------------------------------------
    call getfil (season_wes_file, locfn, 0)
    call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
    !---------------------------------------------------------------------------
    ! 	... get the dimensions
    !---------------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lat', dimid )
    ierr = pio_inq_dimlen( piofile, dimid, nlat_lai )
    ierr = pio_inq_dimid( piofile, 'pft', dimid )
    ierr = pio_inq_dimlen( piofile, dimid, npft_lai )
    !---------------------------------------------------------------------------
    ! 	... allocate arrays
    !---------------------------------------------------------------------------
    allocate( lat_lai(nlat_lai), wk_lai(nlat_lai,npft_lai,12), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
       call endrun
    end if
    !---------------------------------------------------------------------------
    ! 	... read the latitude and the season indicies
    !---------------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'lat', vid )
    ierr = pio_get_var( piofile, vid, lat_lai )

    ierr = pio_inq_varid( piofile, 'season_wes', vid )
    ierr = pio_get_var( piofile, vid, wk_lai )

    call pio_closefile( piofile )


    if(dycore_is('UNSTRUCTURED') ) then
       ! For unstructured grids plon is the 1d horizontal grid size and plat=1
       allocate(clat(plon))
       call get_horiz_grid_d(plon, clat_d_out=clat)
       jl = 1
       ju = plon
    else
       allocate(clat(plat))
       call get_horiz_grid_d(plat, clat_d_out=clat)
       jl = 1
       ju = plat
    end if
    imin = 1
    do j = 1,ju
       diff_min = 10._r8
       pos_min  = -99
       target_lat = clat(j)*r2d
       do i = imin,nlat_lai
          if( abs(lat_lai(i) - target_lat) < diff_min ) then
             diff_min = abs(lat_lai(i) - target_lat)
             pos_min  = i
          end if
       end do
       if( pos_min < 0 ) then
          write(iulog,*) 'dvel_inti: cannot find ',target_lat,' at j,pos_min,diff_min = ',j,pos_min,diff_min
          write(iulog,*) 'dvel_inti: imin,nlat_lai = ',imin,nlat_lai
          write(iulog,*) 'dvel_inti: lat_lai'
          write(iulog,'(1p,10g12.5)') lat_lai(:)
          call endrun
       end if
       if(dycore_is('UNSTRUCTURED') ) then
          imin=1
       else
          imin = pos_min
       end if
       index_season_lai_j(:,:) = wk_lai(pos_min,:,:)

       !---------------------------------------------------------------------------
       ! specify the season as the most frequent in the 11 vegetation classes
       ! this was done to remove a banding problem in dvel (JFL Oct 04)
       !---------------------------------------------------------------------------
       do m = 1,12
          num_seas = 0
          do l = 1,11
             do k = 1,5
                if( index_season_lai_j(l,m) == k ) then
                   num_seas(k) = num_seas(k) + 1
                   exit
                end if
             end do
          end do

          num_max = -1
          do k = 1,5
             if( num_seas(k) > num_max ) then
                num_max = num_seas(k)
                k_max = k
             endif
          end do

          index_season_lai(j,m) = k_max
       end do
    end do

    deallocate( lat_lai, wk_lai, clat, index_season_lai_j)

  end subroutine dvel_inti_xactive

  !-------------------------------------------------------------------------------------
  subroutine get_landuse_and_soilw_from_file(do_soilw)
    use cam_pio_utils, only : cam_pio_openfile
    use ncdio_atm, only : infld
    logical, intent(in) :: do_soilw
    logical :: readvar
    
    type(file_desc_t) :: piofile
    character(len=shr_kind_cl) :: locfn
    logical :: lexist
    
    call getfil (drydep_srf_file, locfn, 1, lexist)
    if(lexist) then
       call cam_pio_openfile(piofile, locfn, PIO_NOWRITE)

       call infld('fraction_landuse', piofile, 'ncol','class',' ',1,pcols,1,n_land_type, begchunk,endchunk, &
            fraction_landuse, readvar, 'PHYS')

       if(do_soilw) then
          call infld('soilw', piofile, 'ncol','month',' ',1,pcols,1,12, begchunk,endchunk, &
               soilw_3d, readvar, 'PHYS')
       end if

       call pio_closefile(piofile)
    else
       call endrun('Unstructured grids require drydep_srf_file ')
    end if


  end subroutine get_landuse_and_soilw_from_file

  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg, lat_veg_edge, &
! Subprogram not used                          lon_veg, lon_veg_edge, landmask, urban, lake, &
! Subprogram not used                          wetland, vegetation_map, soilw_map, do_soilw )
! Subprogram not used 
! Subprogram not used     use mo_constants, only : r2d
! Subprogram not used     use scamMod, only : latiop,loniop,scmlat,scmlon,use_camiop
! Subprogram not used     use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
! Subprogram not used     use filenames, only: ncdata
! Subprogram not used     use dycore, only : dycore_is
! Subprogram not used     use phys_grid, only : scatter_field_to_chunk
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     ! 	... dummy arguments
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     integer, intent(in)      ::  plon, plat, nlon_veg, nlat_veg, npft_veg
! Subprogram not used     real(r8), pointer            :: soilw_map(:,:,:)
! Subprogram not used     real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
! Subprogram not used     real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
! Subprogram not used     real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
! Subprogram not used     real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
! Subprogram not used     real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
! Subprogram not used     real(r8), intent(in)         :: lon_veg(nlon_veg)
! Subprogram not used     real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
! Subprogram not used     real(r8), intent(in)         :: lat_veg(nlat_veg)
! Subprogram not used     real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)
! Subprogram not used     logical,  intent(in)         :: do_soilw
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     ! 	... local variables
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     real(r8) :: closelat,closelon
! Subprogram not used     integer :: latidx,lonidx
! Subprogram not used 
! Subprogram not used     integer, parameter           :: veg_ext = 20
! Subprogram not used     type(file_desc_t)            :: piofile
! Subprogram not used     integer                      :: i, j, ii, jj, jl, ju, i_ndx, n
! Subprogram not used     integer, dimension(plon+1)   :: ind_lon
! Subprogram not used     integer, dimension(plat+1)  :: ind_lat
! Subprogram not used     real(r8)                         :: total_land
! Subprogram not used     real(r8), dimension(plon+1)      :: lon_edge
! Subprogram not used     real(r8), dimension(plat+1)     :: lat_edge
! Subprogram not used     real(r8)                         :: lat1, lat2, lon1, lon2
! Subprogram not used     real(r8)                         :: x1, x2, y1, y2, dx, dy
! Subprogram not used     real(r8)                         :: area, total_area
! Subprogram not used     real(r8), dimension(npft_veg+3)  :: fraction
! Subprogram not used     real(r8)                         :: total_soilw_area
! Subprogram not used     real(r8)                         :: fraction_soilw
! Subprogram not used     real(r8)                         :: total_soilw(12)
! Subprogram not used 
! Subprogram not used     real(r8),    dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
! Subprogram not used     integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext
! Subprogram not used 
! Subprogram not used     real(r8), allocatable :: lam(:), phi(:), garea(:)
! Subprogram not used 
! Subprogram not used     logical, parameter :: has_npole = .true.
! Subprogram not used     integer :: ploniop,platiop
! Subprogram not used     character(len=shr_kind_cl) :: ncdata_loc
! Subprogram not used     real(r8) :: tmp_frac_lu(plon,n_land_type,plat), tmp_soilw_3d(plon,12,plat)
! Subprogram not used 
! Subprogram not used     if(dycore_is('UNSTRUCTURED') ) then
! Subprogram not used        ! For unstructured grids plon is the 1d horizontal grid size and plat=1
! Subprogram not used        allocate(lam(plon), phi(plon))
! Subprogram not used        call get_horiz_grid_d(plon, clat_d_out=phi)
! Subprogram not used     else
! Subprogram not used        allocate(lam(plon), phi(plat))
! Subprogram not used        call get_horiz_grid_d(plat, clat_d_out=phi)
! Subprogram not used     endif
! Subprogram not used     call get_horiz_grid_d(plon, clon_d_out=lam)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     jl = 1
! Subprogram not used     ju = plon
! Subprogram not used 
! Subprogram not used     if (single_column) then
! Subprogram not used        if (use_camiop) then
! Subprogram not used           call getfil (ncdata, ncdata_loc)
! Subprogram not used           call cam_pio_openfile (piofile, trim(ncdata_loc), PIO_NOWRITE)
! Subprogram not used           call shr_scam_getCloseLatLon(piofile%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
! Subprogram not used           call pio_closefile ( piofile)
! Subprogram not used           ploniop=size(loniop)
! Subprogram not used           platiop=size(latiop)
! Subprogram not used        else 
! Subprogram not used           latidx=1
! Subprogram not used           lonidx=1
! Subprogram not used           ploniop=1
! Subprogram not used           platiop=1
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        lon_edge(1) = loniop(lonidx) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
! Subprogram not used 
! Subprogram not used        if (lonidx.lt.ploniop) then
! Subprogram not used           lon_edge(2) = loniop(lonidx+1) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
! Subprogram not used        else
! Subprogram not used           lon_edge(2) = lon_edge(1) + (loniop(2) - loniop(1)) * r2d
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        lat_edge(1) = latiop(latidx) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
! Subprogram not used 
! Subprogram not used        if (latidx.lt.platiop) then
! Subprogram not used           lat_edge(2) = latiop(latidx+1) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
! Subprogram not used        else
! Subprogram not used           lat_edge(2) = lat_edge(1) + (latiop(2) - latiop(1)) * r2d
! Subprogram not used        end if
! Subprogram not used     else
! Subprogram not used        do i = 1,plon
! Subprogram not used           lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
! Subprogram not used        end do
! Subprogram not used        lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
! Subprogram not used        if( .not. has_npole ) then
! Subprogram not used           do j = 1,plat+1
! Subprogram not used              lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
! Subprogram not used           end do
! Subprogram not used        else
! Subprogram not used           do j = 1,plat
! Subprogram not used              lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
! Subprogram not used           end do
! Subprogram not used           lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     do j = 1,plat+1
! Subprogram not used        lat_edge(j) = min( lat_edge(j), 90._r8 )
! Subprogram not used        lat_edge(j) = max( lat_edge(j),-90._r8 )
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     ! wrap around the longitudes
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     do i = -veg_ext,0
! Subprogram not used        lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
! Subprogram not used        mapping_ext     (i) =              nlon_veg+i
! Subprogram not used     end do
! Subprogram not used     do i = 1,nlon_veg
! Subprogram not used        lon_veg_edge_ext(i) = lon_veg_edge(i)
! Subprogram not used        mapping_ext     (i) =              i
! Subprogram not used     end do
! Subprogram not used     do i = nlon_veg+1,nlon_veg+veg_ext
! Subprogram not used        lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
! Subprogram not used        mapping_ext     (i) =              i-nlon_veg
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do j = 1,plon+1
! Subprogram not used        lon1 = lon_edge(j) 
! Subprogram not used        do i = -veg_ext,nlon_veg+veg_ext
! Subprogram not used           dx = lon_veg_edge_ext(i  ) - lon1
! Subprogram not used           dy = lon_veg_edge_ext(i+1) - lon1
! Subprogram not used           if( dx*dy <= 0._r8 ) then
! Subprogram not used              ind_lon(j) = i
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do j = 1,plat+1
! Subprogram not used        lat1 = lat_edge(j)
! Subprogram not used        do i = 1,nlat_veg
! Subprogram not used           dx = lat_veg_edge(i  ) - lat1
! Subprogram not used           dy = lat_veg_edge(i+1) - lat1
! Subprogram not used           if( dx*dy <= 0._r8 ) then
! Subprogram not used              ind_lat(j) = i
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     lat_loop : do j = 1,plat
! Subprogram not used        lon_loop : do i = 1,plon
! Subprogram not used           total_area       = 0._r8
! Subprogram not used           fraction         = 0._r8
! Subprogram not used           total_soilw(:)   = 0._r8
! Subprogram not used           total_soilw_area = 0._r8
! Subprogram not used           do jj = ind_lat(j),ind_lat(j+1)
! Subprogram not used              y1 = max( lat_edge(j),lat_veg_edge(jj) )
! Subprogram not used              y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) ) 
! Subprogram not used              dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
! Subprogram not used              do ii =ind_lon(i),ind_lon(i+1)
! Subprogram not used                 i_ndx = mapping_ext(ii)
! Subprogram not used                 x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
! Subprogram not used                 x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) ) 
! Subprogram not used                 dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
! Subprogram not used                 area = dx * dy
! Subprogram not used                 total_area = total_area + area
! Subprogram not used                 !-----------------------------------------------------------------
! Subprogram not used                 ! 	... special case for ocean grid point 
! Subprogram not used                 !-----------------------------------------------------------------
! Subprogram not used                 if( nint(landmask(i_ndx,jj)) == 0 ) then
! Subprogram not used                    fraction(npft_veg+1) = fraction(npft_veg+1) + area
! Subprogram not used                 else
! Subprogram not used                    do n = 1,npft_veg
! Subprogram not used                       fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
! Subprogram not used                    end do
! Subprogram not used                    fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
! Subprogram not used                    fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
! Subprogram not used                    fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
! Subprogram not used                    !-----------------------------------------------------------------
! Subprogram not used                    ! 	... check if land accounts for the whole area.
! Subprogram not used                    !           If not, the remaining area is in the ocean
! Subprogram not used                    !-----------------------------------------------------------------
! Subprogram not used                    total_land = sum(vegetation_map(i_ndx,jj,:)) &
! Subprogram not used                               + urban  (i_ndx,jj) &
! Subprogram not used                               + lake   (i_ndx,jj) &
! Subprogram not used                               + wetland(i_ndx,jj)
! Subprogram not used                    if( total_land < 1._r8 ) then
! Subprogram not used                       fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
! Subprogram not used                    end if
! Subprogram not used                    !-------------------------------------------------------------------------------------
! Subprogram not used                    ! 	... compute weighted average of soilw over grid (non-water only)
! Subprogram not used                    !-------------------------------------------------------------------------------------
! Subprogram not used                    if( do_soilw ) then
! Subprogram not used                       fraction_soilw = total_land  - (lake(i_ndx,jj) + wetland(i_ndx,jj))
! Subprogram not used                       total_soilw_area = total_soilw_area + fraction_soilw * area
! Subprogram not used                       total_soilw(:)   = total_soilw(:) + fraction_soilw * area * soilw_map(i_ndx,jj,:)
! Subprogram not used                    end if
! Subprogram not used                 end if
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used           !-------------------------------------------------------------------------------------
! Subprogram not used           ! 	... divide by total area of grid box
! Subprogram not used           !-------------------------------------------------------------------------------------
! Subprogram not used           fraction(:) = fraction(:)/total_area
! Subprogram not used           !-------------------------------------------------------------------------------------
! Subprogram not used           ! 	... make sure we don't have too much or too little
! Subprogram not used           !-------------------------------------------------------------------------------------
! Subprogram not used           if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
! Subprogram not used              fraction(:) = fraction(:)/sum(fraction)
! Subprogram not used           end if
! Subprogram not used           !-------------------------------------------------------------------------------------
! Subprogram not used           ! 	... map to Wesely land classification
! Subprogram not used           !-------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           tmp_frac_lu(i, 1, j) =     fraction(20)
! Subprogram not used           tmp_frac_lu(i, 2, j) = sum(fraction(16:17))
! Subprogram not used           tmp_frac_lu(i, 3, j) = sum(fraction(13:15))
! Subprogram not used           tmp_frac_lu(i, 4, j) = sum(fraction( 5: 9))
! Subprogram not used           tmp_frac_lu(i, 5, j) = sum(fraction( 2: 4))
! Subprogram not used           tmp_frac_lu(i, 6, j) =     fraction(19)
! Subprogram not used           tmp_frac_lu(i, 7, j) =     fraction(18)
! Subprogram not used           tmp_frac_lu(i, 8, j) =     fraction( 1)
! Subprogram not used           tmp_frac_lu(i, 9, j) = 0._r8
! Subprogram not used           tmp_frac_lu(i,10, j) = 0._r8
! Subprogram not used           tmp_frac_lu(i,11, j) = sum(fraction(10:12))
! Subprogram not used           if( do_soilw ) then
! Subprogram not used              if( total_soilw_area > 0._r8 ) then
! Subprogram not used                 tmp_soilw_3d(i,:,j) = total_soilw(:)/total_soilw_area
! Subprogram not used              else
! Subprogram not used                 tmp_soilw_3d(i,:,j) = -99._r8
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do lon_loop
! Subprogram not used     end do lat_loop
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     ! 	... reshape according to lat-lon blocks
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     call scatter_field_to_chunk(1,n_land_type,1,plon,tmp_frac_lu,fraction_landuse)
! Subprogram not used     if(do_soilw) call scatter_field_to_chunk(1,12,1,plon,tmp_soilw_3d,soilw_3d)
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     ! 	... make sure there are no out of range values
! Subprogram not used     !-------------------------------------------------------------------------------------
! Subprogram not used     where (fraction_landuse < 0._r8) fraction_landuse = 0._r8
! Subprogram not used     where (fraction_landuse > 1._r8) fraction_landuse = 1._r8
! Subprogram not used 
! Subprogram not used   end subroutine interp_map
  
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_xactive( ncdate, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvel, dflx, mmr, &
                             tv, soilw, rh, ncol, lonndx, latndx, lchnk, &
                             ocnfrc, icefrc, beglandtype, endlandtype )
    !-------------------------------------------------------------------------------------
    !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
    !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
    !   vol. 20, p. 949-964) for calculation of r_a and r_b
    !
    !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
    !   is kept constant where i represents a subgrid environment and a the
    !   grid average environment. thus the calculation proceeds as follows:
    !   va the grid averaged wind is calculated on dots
    !   z0(i) the grid averaged roughness coefficient is calculated
    !   ri(i) the grid averaged richardson number is calculated
    !   --> the grid averaged (u_a)(u*_a) is calculated
    !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
    !   --> final deposotion velocity is weighted average of subgrid scale velocities
    !
    ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
    ! modified by JFL to be used in MOZART-2 (October 2002)
    !-------------------------------------------------------------------------------------

    use seq_drydep_mod, only: z0, rgso, rgss, h2_a, h2_b, h2_c, ri, rclo, rcls, rlu, rac
    use seq_drydep_mod, only: seq_drydep_setHCoeff, foxd, drat
    use physconst,      only: tmelt
    use seq_drydep_mod, only: drydep_method,  DD_XLND

    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)   :: ncol
    integer, intent(in)   :: ncdate                   ! present date (yyyymmdd)
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: rh(ncol,1)               ! relative humidity
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)              
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: soilw(pcols)             ! soil moisture fraction
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvel(ncol,gas_pcnst)        ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    integer, intent(in)     ::   latndx(pcols)           ! chunk latitude indicies
    integer, intent(in)     ::   lonndx(pcols)           ! chunk longitude indicies
    integer, intent(in)     ::   lchnk                   ! chunk number

    integer, intent(in), optional     ::  beglandtype
    integer, intent(in), optional     ::  endlandtype

    real(r8), intent(in), optional      :: ocnfrc(pcols) 
    real(r8), intent(in), optional      :: icefrc(pcols) 

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8), parameter :: scaling_to_cm_per_s = 100._r8
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    integer :: i, ispec, lt, m
    integer :: sndx
    integer :: month

    real(r8) :: slope = 0._r8
    real(r8) :: z0water ! revised z0 over water
    real(r8) :: p       ! pressure at midpoint first layer
    real(r8) :: pg      ! surface pressure
    real(r8) :: es      ! saturation vapor pressure
    real(r8) :: ws      ! saturation mixing ratio
    real(r8) :: hvar    ! constant to compute xmol
    real(r8) :: h       ! constant to compute xmol
    real(r8) :: psih    ! stability correction factor
    real(r8) :: rs      ! constant for calculating rsmx
    real(r8) :: rmx     ! resistance by vegetation
    real(r8) :: zovl    ! ratio of z to  m-o length
    real(r8) :: cvarb   ! cvar averaged over landtypes
    real(r8) :: bb      ! b averaged over landtypes
    real(r8) :: ustarb  ! ustar averaged over landtypes
    real(r8) :: tc(ncol)  ! temperature in celsius
    real(r8) :: cts(ncol) ! correction to rlu rcl and rgs for frost

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,nddvels) :: heff

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location only
    !-------------------------------------------------------------------------------------
    integer                :: index_season(ncol,n_land_type)
    real(r8), dimension(ncol) :: tha     ! atmospheric virtual potential temperature
    real(r8), dimension(ncol) :: thg     ! ground virtual potential temperature
    real(r8), dimension(ncol) :: z       ! height of lowest level
    real(r8), dimension(ncol) :: va      ! magnitude of v on cross points
    real(r8), dimension(ncol) :: ribn    ! richardson number
    real(r8), dimension(ncol) :: qs      ! saturation specific humidity
    real(r8), dimension(ncol) :: crs     ! multiplier to calculate crs
    real(r8), dimension(ncol) :: rdc     ! part of lower canopy resistance
    real(r8), dimension(ncol) :: uustar  ! u*ustar (assumed constant over grid)
    real(r8), dimension(ncol) :: z0b     ! average roughness length over grid
    real(r8), dimension(ncol) :: wrk     ! work array
    real(r8), dimension(ncol) :: term    ! work array
    real(r8), dimension(ncol) :: resc    ! work array
    real(r8), dimension(ncol) :: lnd_frc ! work array
    logical,  dimension(ncol) :: unstable
    logical,  dimension(ncol) :: has_rain
    logical,  dimension(ncol) :: has_dew

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and landtype
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,n_land_type) :: rds   ! resistance for deposition of sulfate
    real(r8), dimension(ncol,n_land_type) :: b     ! buoyancy parameter for unstable conditions
    real(r8), dimension(ncol,n_land_type) :: cvar  ! height parameter
    real(r8), dimension(ncol,n_land_type) :: ustar ! friction velocity
    real(r8), dimension(ncol,n_land_type) :: xmol  ! monin-obukhov length

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location, landtype and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rsmx  ! vegetative resistance (plant mesophyll)
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rclx  ! lower canopy resistance
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rlux  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,n_land_type) :: rlux_o3  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rgsx  ! ground resistance
    real(r8) :: pmid(ncol,1)                             ! for seasalt aerosols
    real(r8) :: tfld(ncol,1)                             ! for seasalt aerosols
    real(r8) :: fact, vds
    real(r8) :: rc                                    ! combined surface resistance
    real(r8) :: var_soilw, dv_soil_h2, fact_h2        ! h2 dvel wrking variables
    logical  :: fr_lnduse(ncol,n_land_type)           ! wrking array
    real(r8) :: dewm                                  ! multiplier for rs when dew occurs

    real(r8) :: lcl_frc_landuse(ncol,n_land_type) 

    integer :: beglt, endlt

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(r8) :: dv_pan
    real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
    real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)

    if (present( beglandtype)) then
      beglt = beglandtype 
    else
      beglt = 1
    endif
    if (present( endlandtype)) then
      endlt = endlandtype 
    else
      endlt = n_land_type
    endif
  
    !-------------------------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------------------------
    do m = 1,gas_pcnst
       dvel(:,m) = 0._r8
    end do

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !-------------------------------------------------------------------------------------
    ! define species-dependent parameters (temperature dependent)
    !-------------------------------------------------------------------------------------
    call seq_drydep_setHCoeff( ncol, sfc_temp, heff )

    do lt = 1,n_land_type
       dep_ra (:,lt,lchnk)   = 0._r8
       dep_rb (:,lt,lchnk)   = 0._r8
       rds(:,lt)   = 0._r8
    end do

    !-------------------------------------------------------------------------------------
    ! 	... set month
    !-------------------------------------------------------------------------------------
    month = mod( ncdate,10000 )/100

    !-------------------------------------------------------------------------------------
    ! define which season (relative to Northern hemisphere climate)
    !-------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------
    ! define season index based on fixed LAI
    !-------------------------------------------------------------------------------------
    if ( drydep_method == DD_XLND ) then
       index_season = 4
    else
       do i = 1,ncol
          index_season(i,:) = index_season_lai(latndx(i),month)
       end do
    endif
    !-------------------------------------------------------------------------------------
    ! special case for snow covered terrain
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       if( snow(i) > .01_r8 ) then
          index_season(i,:) = 4
       end if
    end do
    !-------------------------------------------------------------------------------------
    ! scale rain and define logical arrays
    !-------------------------------------------------------------------------------------
    has_rain(:ncol) = rain(:ncol) > rain_threshold

    !-------------------------------------------------------------------------------------
    ! loop over longitude points
    !-------------------------------------------------------------------------------------
    col_loop :  do i = 1,ncol
       p   = pressure_10m(i)
       pg  = pressure_sfc(i)
       !-------------------------------------------------------------------------------------
       ! potential temperature
       !-------------------------------------------------------------------------------------
       tha(i) = air_temp(i) * (p00/p )**rovcp * (1._r8 + .61_r8*spec_hum(i))
       thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._r8 + .61_r8*spec_hum(i))
       !-------------------------------------------------------------------------------------
       ! height of 1st level
       !-------------------------------------------------------------------------------------
       z(i) = - r/grav * air_temp(i) * (1._r8 + .61_r8*spec_hum(i)) * log(p/pg)
       !-------------------------------------------------------------------------------------
       ! wind speed
       !-------------------------------------------------------------------------------------
       va(i) = max( .01_r8,wind_speed(i) )
       !-------------------------------------------------------------------------------------
       ! Richardson number
       !-------------------------------------------------------------------------------------
       ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
       ribn(i) = min( ribn(i),ric )
       unstable(i) = ribn(i) < 0._r8
       !-------------------------------------------------------------------------------------
       ! saturation vapor pressure (Pascals)
       ! saturation mixing ratio
       ! saturation specific humidity
       !-------------------------------------------------------------------------------------
       es    = 611._r8*exp( 5414.77_r8*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
       ws    = .622_r8*es/(pg - es)
       qs(i) = ws/(1._r8 + ws)
       has_dew(i) = .false.
       if( qs(i) <= spec_hum(i) ) then
          has_dew(i) = .true.
       end if
       if( sfc_temp(i) < tmelt ) then
          has_dew(i) = .false.
       end if
       !-------------------------------------------------------------------------------------
       ! constant in determining rs
       !-------------------------------------------------------------------------------------
       tc(i) = sfc_temp(i) - tmelt
       if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_r8 ) then
          crs(i) = (1._r8 + (200._r8/(solar_flux(i) + .1_r8))**2) * (400._r8/(tc(i)*(40._r8 - tc(i))))
       else
          crs(i) = large_value
       end if
       !-------------------------------------------------------------------------------------
       ! rdc (lower canopy res)
       !-------------------------------------------------------------------------------------
       rdc(i) = 100._r8*(1._r8 + 1000._r8/(solar_flux(i) + 10._r8))/(1._r8 + 1000._r8*slope)
    end do col_loop

    !-------------------------------------------------------------------------------------
    ! 	... form working arrays
    !-------------------------------------------------------------------------------------
    do lt = 1,n_land_type
       do i=1,ncol
          if ( drydep_method == DD_XLND ) then
             lcl_frc_landuse(i,lt) = 0._r8
          else
             lcl_frc_landuse(i,lt) = fraction_landuse(i,lt,lchnk)
          endif
       enddo
    end do
    if ( present(ocnfrc) .and. present(icefrc) ) then
       do i=1,ncol
          ! land type 7 is used for ocean
          ! land type 8 is used for sea ice
          lcl_frc_landuse(i,7) = ocnfrc(i)
          lcl_frc_landuse(i,8) = icefrc(i)
       enddo
    endif
    do lt = 1,n_land_type
       do i=1,ncol
          fr_lnduse(i,lt) = lcl_frc_landuse(i,lt) > 0._r8
       enddo
    end do

    !-------------------------------------------------------------------------------------
    ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
    ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
    !-------------------------------------------------------------------------------------
    z0b(:) = 0._r8
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             z0b(i) = z0b(i) + lcl_frc_landuse(i,lt) * log( z0(index_season(i,lt),lt) )
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! find the constant velocity uu*=(u_i)(u*_i)
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       z0b(i) = exp( z0b(i) )
       cvarb  = vonkar/log( z(i)/z0b(i) )
       !-------------------------------------------------------------------------------------
       ! unstable and stable cases
       !-------------------------------------------------------------------------------------
       if( unstable(i) ) then
          bb = 9.4_r8*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*bb)) )
       else
          ustarb = cvarb * va(i)/(1._r8 + 4.7_r8*ribn(i))
       end if
       uustar(i) = va(i)*ustarb
    end do

    !-------------------------------------------------------------------------------------
    ! calculate the friction velocity for each land type u_i=uustar/u*_i
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( unstable(i) ) then
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)* sqrt( abs(ribn(i))*z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)*sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*b(i,lt))) ) )
             else
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
             end if
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! revise calculation of friction velocity and z0 over water
    !-------------------------------------------------------------------------------------
    lt = 7    
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          if( unstable(i) ) then
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log( z(i)/z0water ))
             b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)*sqrt( abs(ribn(i))*z(i)/z0water )
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)* sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8+ 7.4_r8*b(i,lt))) ) )
          else
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log(z(i)/z0water))
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
          end if
       end if
    end do

    !-------------------------------------------------------------------------------------
    ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             hvar = (va(i)/0.74_r8) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
             if( unstable(i) ) then                      ! unstable
                h = hvar*(1._r8 - (9.4_r8*ribn(i)/(1._r8 + 5.3_r8*b(i,lt))))
             else
                h = hvar/((1._r8+4.7_r8*ribn(i))**2)
             end if
             xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! psih
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( xmol(i,lt) < 0._r8 ) then
                zovl = z(i)/xmol(i,lt)
                zovl = max( -1._r8,zovl )
                psih = exp( .598_r8 + .39_r8*log( -zovl ) - .09_r8*(log( -zovl ))**2 )
                vds  = 2.e-3_r8*ustar(i,lt) * (1._r8 + (300/(-xmol(i,lt)))**0.666_r8)
             else
                zovl = z(i)/xmol(i,lt)
                zovl = min( 1._r8,zovl )
                psih = -5._r8 * zovl
                vds  = 2.e-3_r8*ustar(i,lt)
             end if
             dep_ra (i,lt,lchnk) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
             dep_rb (i,lt,lchnk) = (2._r8/(vonkar*ustar(i,lt))) * crb
             rds(i,lt) = 1._r8/vds
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! surface resistance : depends on both land type and species
    ! land types are computed seperately, then resistance is computed as average of values
    ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
    !
    ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
    !-------------------------------------------------------------------------------------
    species_loop1 :  do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          m = map_dvel(ispec)
          do lt = beglt,endlt
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   sndx = index_season(i,lt)
                   if( ispec == o3_ndx .or. ispec == o3a_ndx .or. ispec == so2_ndx ) then
                      rmx = 0._r8
                   else
                      rmx = 1._r8/(heff(i,m)/3000._r8 + 100._r8*foxd(m))
                   end if
                   cts(i) = 1000._r8*exp( - tc(i) - 4._r8 )                 ! correction for frost
                   rgsx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rgss(sndx,lt))) + (foxd(m)/rgso(sndx,lt)))
                   !-------------------------------------------------------------------------------------
                   ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                   !-------------------------------------------------------------------------------------
                   if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                      if( ispec == co_ndx ) then
                         fact_h2 = 1.0_r8
                      elseif ( ispec == h2_ndx ) then
                         fact_h2 = 0.5_r8
                      elseif ( ispec == ch4_ndx ) then
                         fact_h2 = 50.0_r8
                      end if
                      !-------------------------------------------------------------------------------------
                      ! no deposition on snow, ice, desert, and water
                      !-------------------------------------------------------------------------------------
                      if( lt == 1 .or. lt == 7 .or. lt == 8 .or. sndx == 4 ) then
                         rgsx(i,lt,ispec) = large_value
                      else
                         var_soilw = max( .1_r8,min( soilw(i),.3_r8 ) )
                         if( lt == 3 ) then
                            var_soilw = log( var_soilw )
                         end if
                         dv_soil_h2 = h2_c(lt) + var_soilw*(h2_b(lt) + var_soilw*h2_a(lt))
                         if( dv_soil_h2 > 0._r8 ) then
                            rgsx(i,lt,ispec) = fact_h2/(dv_soil_h2*1.e-4_r8)
                         end if
                      end if
                   end if
                   if( lt == 7 ) then
                      rclx(i,lt,ispec) = large_value
                      rsmx(i,lt,ispec) = large_value
                      rlux(i,lt,ispec) = large_value
                   else
                      rs = ri(sndx,lt)*crs(i)
                      if ( has_dew(i) .or. has_rain(i) ) then
                         dewm = 3._r8
                      else
                         dewm = 1._r8
                      end if
                      rsmx(i,lt,ispec) = (dewm*rs*drat(m) + rmx)
                      !-------------------------------------------------------------------------------------
                      ! jfl : special case for PAN
                      !-------------------------------------------------------------------------------------
                      if( ispec == pan_ndx .or. ispec == xpan_ndx ) then
                         dv_pan =  c0_pan(lt) * (1._r8 - exp( -k_pan(lt)*(dewm*rs*drat(m))*1.e-2_r8 ))
                         if( dv_pan > 0._r8 .and. sndx /= 4 ) then
                            rsmx(i,lt,ispec) = ( 1._r8/dv_pan )
                         end if
                      end if
                      rclx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rcls(sndx,lt))) + (foxd(m)/rclo(sndx,lt)))
                      rlux(i,lt,ispec) = cts(i) + rlu(sndx,lt)/(1.e-5_r8*heff(i,m) + foxd(m))
                   end if
                end if
             end do
          end do
       end if
    end do species_loop1

    do lt = beglt,endlt
       if( lt /= 7 ) then
          do i = 1,ncol
             if( fr_lnduse(i,lt) ) then
                sndx = index_season(i,lt)
                !-------------------------------------------------------------------------------------
                ! 	... no effect if sfc_temp < O C
                !-------------------------------------------------------------------------------------
                if( sfc_temp(i) > tmelt ) then
                   if( has_dew(i) ) then
                      rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   end if
                   if( has_rain(i) ) then
                      ! rlux(i,lt,o3_ndx) = 1./(1.e-3 + (1./(3.*rlu(sndx,lt))))
                      rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + 3._r8*rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   end if
                end if

                if ( o3_ndx > 0 ) then
                   rclx(i,lt,o3_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3_ndx) = cts(i) + rlux(i,lt,o3_ndx)
                end if
                if ( o3a_ndx > 0 ) then
                   rclx(i,lt,o3a_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3a_ndx) = cts(i) + rlux(i,lt,o3a_ndx)
                end if

             end if
          end do
       end if
    end do

    species_loop2 : do ispec = 1,gas_pcnst
       m = map_dvel(ispec)
       if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= o3a_ndx .and. ispec /= so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( has_dew(i) ) then
                               rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
                                    + 1.e-7_r8*heff(i,m) + foxd(m)/rlux_o3(i,lt))
                            end if
                         end if

                      end if
                   end do
                end if
             end do
          else if( ispec == so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( qs(i) <= spec_hum(i) ) then
                               rlux(i,lt,ispec) = 100._r8
                            end if
                            if( has_rain(i) ) then
                               !                               rlux(i,lt,ispec) = 1./(2.e-4 + (1./(3.*rlu(index_season(i,lt),lt))))
                               rlux(i,lt,ispec) = 15._r8*rlu(index_season(i,lt),lt)/(5._r8 + 3.e-3_r8*rlu(index_season(i,lt),lt))
                            end if
                         end if
                         rclx(i,lt,ispec) = cts(i) + rcls(index_season(i,lt),lt)
                         rlux(i,lt,ispec) = cts(i) + rlux(i,lt,ispec)

                      end if
                   end do
                end if
             end do
             do i = 1,ncol
                if( fr_lnduse(i,1) .and. (has_dew(i) .or. has_rain(i)) ) then
                   rlux(i,1,ispec) = 50._r8
                end if
             end do
          end if
       end if
    end do species_loop2

    !-------------------------------------------------------------------------------------
    ! compute rc
    !-------------------------------------------------------------------------------------
    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))
    species_loop3 : do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          wrk(:) = 0._r8
          lt_loop: do lt = beglt,endlt
             do i = 1,ncol
                if (fr_lnduse(i,lt)) then
                   resc(i) = 1._r8/( 1._r8/rsmx(i,lt,ispec) + 1._r8/rlux(i,lt,ispec) &
                                   + 1._r8/(rdc(i) + rclx(i,lt,ispec)) &
                                   + 1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))

                   resc(i) = max( 10._r8,resc(i) )

                   lnd_frc(i) = lcl_frc_landuse(i,lt)
                endif
             enddo
             !-------------------------------------------------------------------------------------
             ! 	... compute average deposition velocity
             !-------------------------------------------------------------------------------------
             select case( solsym(ispec) )
             case( 'SO2' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      ! assume no surface resistance for SO2 over water`
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk)) 
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:))
                   endwhere
                end if
             case( 'SO4' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + rds(:,lt))
                endwhere
             case( 'NH4', 'NH4NO3', 'XNH4NO3' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + 0.5_r8*rds(:,lt))
                endwhere

             !-------------------------------------------------------------------------------------
             !  ... special case for Pb (for consistency with offline code)
             !-------------------------------------------------------------------------------------
             case( 'Pb' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:) * 0.05e-2_r8
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                   endwhere
                end if

             !-------------------------------------------------------------------------------------
             !  ... special case for carbon aerosols
             !-------------------------------------------------------------------------------------
             case( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB','SOAX' )
                if ( drydep_method == DD_XLND ) then
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.10e-2_r8
                   endwhere
                else
                   wrk(:ncol) = 0.10e-2_r8
                endif

             !-------------------------------------------------------------------------------------
             ! deposition over ocean for HCN, CH3CN
             !    velocity estimated from aircraft measurements (E.Apel, INTEX-B)
             !-------------------------------------------------------------------------------------
             case( 'HCN','CH3CN' )
                if( lt == 7 ) then ! over ocean only
                   where( fr_lnduse(:ncol,lt) .and. snow(:ncol) < 0.01_r8  )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                   endwhere
                end if
             case default
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:ncol))
                endwhere
             end select
          end do lt_loop
          dvel(:ncol,ispec) = wrk(:ncol) * scaling_to_cm_per_s
          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)
       end if

    end do species_loop3

    if ( beglt > 1 ) return

    !-------------------------------------------------------------------------------------
    ! 	... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx > 0 ) then
       if( has_dvel(mpan_ndx) ) then
          dvel(:ncol,mpan_ndx) = dvel(:ncol,mpan_ndx)/3._r8
          dflx(:ncol,mpan_ndx) = term(:ncol) * dvel(:ncol,mpan_ndx) * mmr(:ncol,plev,mpan_ndx)
       end if
    end if
    if( xmpan_ndx > 0 ) then
       if( has_dvel(xmpan_ndx) ) then
          dvel(:ncol,xmpan_ndx) = dvel(:ncol,xmpan_ndx)/3._r8
          dflx(:ncol,xmpan_ndx) = term(:ncol) * dvel(:ncol,xmpan_ndx) * mmr(:ncol,plev,xmpan_ndx)
       end if
    end if

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0) then
       if( has_dvel(hcooh_ndx) ) then
          dvel(:ncol,hcooh_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,hcooh_ndx) = term(:ncol) * dvel(:ncol,hcooh_ndx) * mmr(:ncol,plev,hcooh_ndx)
       end if
    end if
!
! SOG species
!
    if( sogm_ndx > 0) then
       if( has_dvel(sogm_ndx) ) then
          dvel(:ncol,sogm_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogm_ndx) = term(:ncol) * dvel(:ncol,sogm_ndx) * mmr(:ncol,plev,sogm_ndx)
       end if
    end if
    if( sogi_ndx > 0) then
       if( has_dvel(sogi_ndx) ) then
          dvel(:ncol,sogi_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogi_ndx) = term(:ncol) * dvel(:ncol,sogi_ndx) * mmr(:ncol,plev,sogi_ndx)
       end if
    end if
    if( sogt_ndx > 0) then
       if( has_dvel(sogt_ndx) ) then
          dvel(:ncol,sogt_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogt_ndx) = term(:ncol) * dvel(:ncol,sogt_ndx) * mmr(:ncol,plev,sogt_ndx)
       end if
    end if
    if( sogb_ndx > 0) then
       if( has_dvel(sogb_ndx) ) then
          dvel(:ncol,sogb_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogb_ndx) = term(:ncol) * dvel(:ncol,sogb_ndx) * mmr(:ncol,plev,sogb_ndx)
       end if
    end if
    if( sogx_ndx > 0) then
       if( has_dvel(sogx_ndx) ) then
          dvel(:ncol,sogx_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogx_ndx) = term(:ncol) * dvel(:ncol,sogx_ndx) * mmr(:ncol,plev,sogx_ndx)
       end if
    end if
!
  end subroutine drydep_xactive

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine soilw_inti( ncfile, nlon_veg, nlat_veg, soilw_map )
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     !	... read primary soil moisture table
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use time_manager,  only : get_calday
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     !	... dummy args
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     integer, intent(in) :: &
! Subprogram not used          nlon_veg, &
! Subprogram not used          nlat_veg
! Subprogram not used     real(r8), pointer :: soilw_map(:,:,:)
! Subprogram not used     character(len=*), intent(in) :: ncfile ! file name of netcdf file containing data
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     !	... local variables
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     integer :: gndx = 0
! Subprogram not used     integer :: nlat, &             ! # of lats in soilw file
! Subprogram not used                nlon                ! # of lons in soilw file
! Subprogram not used     integer :: i, ip, k, m
! Subprogram not used     integer :: j, jl, ju
! Subprogram not used     integer :: lev, day, ierr
! Subprogram not used     type(file_desc_t) :: piofile
! Subprogram not used     type(var_desc_t) :: vid
! Subprogram not used     
! Subprogram not used     integer :: dimid_lat, dimid_lon, dimid_time
! Subprogram not used     integer :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
! Subprogram not used                               716, 816, 915, 1016, 1115, 1216 /)
! Subprogram not used 
! Subprogram not used     character(len=shr_kind_cl) :: locfn
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... open netcdf file
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     call getfil (ncfile, locfn, 0)
! Subprogram not used     call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... get longitudes
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
! Subprogram not used     if( nlon /= nlon_veg ) then
! Subprogram not used        write(iulog,*) 'soilw_inti: soil and vegetation lons differ; ',nlon, nlon_veg
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... get latitudes
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
! Subprogram not used     if( nlat /= nlat_veg ) then
! Subprogram not used        write(iulog,*) 'soilw_inti: soil and vegetation lats differ; ',nlat, nlat_veg
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... set times (days of year)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_dimid( piofile, 'time', dimid_time )
! Subprogram not used     ierr = pio_inq_dimlen( piofile, dimid_time, ndays )
! Subprogram not used     if( ndays /= 12 ) then
! Subprogram not used        write(iulog,*) 'soilw_inti: dataset not a cyclical year'
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     allocate( days(ndays),stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'soilw_inti: days allocation error = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     do m = 1,min(12,ndays)
! Subprogram not used        days(m) = get_calday( dates(m), 0 )
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     !	... allocate arrays
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     allocate( soilw_map(nlon,nlat,ndays), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'soilw_inti: soilw_map allocation error = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     !	... read in the soil moisture
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     ierr = pio_inq_varid( piofile, 'SOILW', vid )
! Subprogram not used     ierr = pio_get_var( piofile, vid, soilw_map )
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     !	... close file
! Subprogram not used     !------------------------------------------------------------------
! Subprogram not used     call pio_closefile( piofile )
! Subprogram not used 
! Subprogram not used   end subroutine soilw_inti
  
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine chk_soilw( calday )
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     !	... check timing for ub values
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use mo_constants, only : dayspy
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     !	... dummy args
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     real(r8), intent(in)    :: calday
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     !	... local variables
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     integer  ::  m, upper
! Subprogram not used     real(r8)     ::  numer, denom
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !	... setup the time interpolation
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( calday < days(1) ) then
! Subprogram not used        next = 1
! Subprogram not used        last = ndays
! Subprogram not used     else
! Subprogram not used        if( days(ndays) < dayspy ) then
! Subprogram not used           upper = ndays
! Subprogram not used        else
! Subprogram not used           upper = ndays - 1
! Subprogram not used        end if
! Subprogram not used        do m = upper,1,-1
! Subprogram not used           if( calday >= days(m) ) then
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        last = m
! Subprogram not used        next = mod( m,ndays ) + 1
! Subprogram not used     end if
! Subprogram not used     numer = calday - days(last)
! Subprogram not used     denom = days(next) - days(last)
! Subprogram not used     if( numer < 0._r8 ) then
! Subprogram not used        numer = dayspy + numer
! Subprogram not used     end if
! Subprogram not used     if( denom < 0._r8 ) then
! Subprogram not used        denom = dayspy + denom
! Subprogram not used     end if
! Subprogram not used     dels = max( min( 1._r8,numer/denom ),0._r8 )
! Subprogram not used 
! Subprogram not used   end subroutine chk_soilw

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
! Subprogram not used   subroutine set_soilw( soilw, lchnk, calday )
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     !	... set the soil moisture
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     !	... dummy args
! Subprogram not used     !--------------------------------------------------------------------
! Subprogram not used     real(r8), intent(inout) :: soilw(pcols)
! Subprogram not used     integer,  intent(in)    :: lchnk           ! chunk indice
! Subprogram not used     real(r8), intent(in)    :: calday
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     integer :: i, ilon,ilat
! Subprogram not used 
! Subprogram not used     call chk_soilw( calday )
! Subprogram not used 
! Subprogram not used     soilw(:) = soilw_3d(:,last,lchnk) + dels *( soilw_3d(:,next,lchnk) - soilw_3d(:,last,lchnk))
! Subprogram not used 
! Subprogram not used   end subroutine set_soilw

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  function has_drydep( name )

    implicit none

    character(len=*), intent(in) :: name

    logical :: has_drydep
    integer :: i

    has_drydep = .false.

    do i=1,nddvels
       if ( trim(name) == trim(drydep_list(i)) ) then
         has_drydep = .true.
         exit
       endif
    enddo

  endfunction has_drydep

end module mo_drydep
