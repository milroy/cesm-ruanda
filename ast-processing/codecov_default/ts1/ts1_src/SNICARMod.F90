module SNICARMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SNICARMod
!
! !DESCRIPTION:
! Calculate albedo of snow containing impurities 
! and the evolution of snow effective radius
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_flush
  use clm_varctl    , only : iulog
  use shr_const_mod , only : SHR_CONST_RHOICE
  use abortutils    , only : endrun

  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SNICAR_RT        ! Snow albedo and vertically-resolved solar absorption
  public :: SnowAge_grain    ! Snow effective grain size evolution
  public :: SnowAge_init     ! Initial read in of snow-aging file
  public :: SnowOptics_init  ! Initial read in of snow-optics file
!
! !PUBLIC DATA MEMBERS:

  real(r8), public, parameter :: snw_rds_min = 54.526_r8  ! minimum allowed snow effective radius
                                                          ! (also "fresh snow" value) [microns]
  integer,  public, parameter :: sno_nbr_aer =   8        ! number of aerosol species in snowpack
                                                          ! (indices described above) [nbr]
  logical,  public, parameter :: DO_SNO_OC =    .false.   ! parameter to include organic carbon (OC)
                                                          ! in snowpack radiative calculations
  logical,  public, parameter :: DO_SNO_AER =   .true.    ! parameter to include aerosols in snowpack radiative calculations

  real(r8), public, parameter :: scvng_fct_mlt_bcphi = 0.20_r8   ! scavenging factor for hydrophillic BC inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_bcpho = 0.03_r8   ! scavenging factor for hydrophobic BC inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocphi = 0.20_r8   ! scavenging factor for hydrophillic OC inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocpho = 0.03_r8   ! scavenging factor for hydrophobic OC inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst1  = 0.02_r8   ! scavenging factor for dust species 1 inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst2  = 0.02_r8   ! scavenging factor for dust species 2 inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst3  = 0.01_r8   ! scavenging factor for dust species 3 inclusion in meltwater
                                                                 ! [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst4  = 0.01_r8   ! scavenging factor for dust species 4 inclusion in meltwater
                                                                 ! [frc]

! !PRIVATE MEMBER FUNCTIONS:

!
! !PRIVATE DATA MEMBERS:
  ! Aerosol species indices:
  !  1= hydrophillic black carbon 
  !  2= hydrophobic black carbon
  !  3= hydrophilic organic carbon
  !  4= hydrophobic organic carbon
  !  5= dust species 1
  !  6= dust species 2
  !  7= dust species 3
  !  8= dust species 4
  integer,  parameter :: numrad_snw  =   5               ! number of spectral bands used in snow model [nbr]
  integer,  parameter :: nir_bnd_bgn =   2               ! first band index in near-IR spectrum [idx]
  integer,  parameter :: nir_bnd_end =   5               ! ending near-IR band index [idx]

  integer,  parameter :: idx_Mie_snw_mx = 1471           ! number of effective radius indices used in Mie lookup table [idx]
  integer,  parameter :: idx_T_max      = 11             ! maxiumum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_T_min      = 1              ! minimum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_max   = 31             ! maxiumum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_min   = 1              ! minimum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_max   = 8              ! maxiumum snow density index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_min   = 1              ! minimum snow density index used in aging lookup table [idx]

  integer,  parameter :: snw_rds_max_tbl = 1500          ! maximum effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_tbl = 30            ! minimium effective radius defined in Mie lookup table [microns]
  real(r8), parameter :: snw_rds_max     = 1500._r8      ! maximum allowed snow effective radius [microns]
  real(r8), parameter :: snw_rds_refrz   = 1000._r8      ! effective radius of re-frozen snow [microns]

  real(r8), parameter :: min_snw = 1.0E-30_r8            ! minimum snow mass required for SNICAR RT calculation [kg m-2]

  !real(r8), parameter :: C1_liq_Brun89 = 1.28E-17_r8    ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89
  real(r8), parameter :: C1_liq_Brun89 = 0._r8           ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: zeroed to accomodate dry snow aging
  real(r8), parameter :: C2_liq_Brun89 = 4.22E-13_r8     ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: corrected for LWC in units of percent

  real(r8), parameter :: tim_cns_bc_rmv  = 2.2E-8_r8     ! time constant for removal of BC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_oc_rmv  = 2.2E-8_r8     ! time constant for removal of OC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_dst_rmv = 2.2E-8_r8     ! time constant for removal of dust in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)

  ! scaling of the snow aging rate (tuning option):
  logical :: flg_snoage_scl    = .false.                 ! flag for scaling the snow aging rate by some arbitrary factor
  real(r8), parameter :: xdrdt = 1.0_r8                  ! arbitrary factor applied to snow aging rate

  ! snow and aerosol Mie parameters:
  ! (arrays declared here, but are set in iniTimeConst)
  ! (idx_Mie_snw_mx is number of snow radii with defined parameters (i.e. from 30um to 1500um))
  
  ! direct-beam weighted ice optical properties
  real(r8) :: ss_alb_snw_drc(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: asm_prm_snw_drc(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: ext_cff_mss_snw_drc(idx_Mie_snw_mx,numrad_snw)

  ! diffuse radiation weighted ice optical properties
  real(r8) :: ss_alb_snw_dfs(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: asm_prm_snw_dfs(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: ext_cff_mss_snw_dfs(idx_Mie_snw_mx,numrad_snw)

  ! hydrophiliic BC
  real(r8) :: ss_alb_bc1(numrad_snw)
  real(r8) :: asm_prm_bc1(numrad_snw)
  real(r8) :: ext_cff_mss_bc1(numrad_snw)

  ! hydrophobic BC
  real(r8) :: ss_alb_bc2(numrad_snw)
  real(r8) :: asm_prm_bc2(numrad_snw)
  real(r8) :: ext_cff_mss_bc2(numrad_snw)

  ! hydrophobic OC
  real(r8) :: ss_alb_oc1(numrad_snw)
  real(r8) :: asm_prm_oc1(numrad_snw)
  real(r8) :: ext_cff_mss_oc1(numrad_snw)

  ! hydrophilic OC
  real(r8) :: ss_alb_oc2(numrad_snw)
  real(r8) :: asm_prm_oc2(numrad_snw)
  real(r8) :: ext_cff_mss_oc2(numrad_snw)

  ! dust species 1:
  real(r8) :: ss_alb_dst1(numrad_snw)
  real(r8) :: asm_prm_dst1(numrad_snw)
  real(r8) :: ext_cff_mss_dst1(numrad_snw)

  ! dust species 2:
  real(r8) :: ss_alb_dst2(numrad_snw)
  real(r8) :: asm_prm_dst2(numrad_snw)
  real(r8) :: ext_cff_mss_dst2(numrad_snw)

  ! dust species 3:
  real(r8) :: ss_alb_dst3(numrad_snw)
  real(r8) :: asm_prm_dst3(numrad_snw)
  real(r8) :: ext_cff_mss_dst3(numrad_snw)

  ! dust species 4:
  real(r8) :: ss_alb_dst4(numrad_snw)
  real(r8) :: asm_prm_dst4(numrad_snw)
  real(r8) :: ext_cff_mss_dst4(numrad_snw)

  ! best-fit parameters for snow aging defined over:
  !  11 temperatures from 225 to 273 K
  !  31 temperature gradients from 0 to 300 K/m
  !   8 snow densities from 0 to 350 kg/m3
  ! (arrays declared here, but are set in iniTimeConst)
  real(r8), pointer :: snowage_tau(:,:,:) ! (idx_rhos_max,idx_Tgrd_max,idx_T_max)
  real(r8), pointer :: snowage_kappa(:,:,:) ! (idx_rhos_max,idx_Tgrd_max,idx_T_max)
  real(r8), pointer :: snowage_drdt0(:,:,:) ! idx_rhos_max,idx_Tgrd_max,idx_T_max)
 
!
! !REVISION HISTORY:
! Created by Mark Flanner
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SNICAR_RT
!
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in module SurfaceAlbedoMod (CLM)
! subroutine albice (CSIM)
!
! !REVISION HISTORY:
! Author: Mark Flanner
!
! !INTERFACE:
  
! Subprogram not used   subroutine SNICAR_RT (flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,  &
! Subprogram not used                         coszen, flg_slr_in, h2osno_liq, h2osno_ice, snw_rds,   &
! Subprogram not used                         mss_cnc_aer_in, albsfc, albout, flx_abs)
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Determine reflectance of, and vertically-resolved solar absorption in, 
! Subprogram not used     ! snow with impurities.
! Subprogram not used     !
! Subprogram not used     ! Original references on physical models of snow reflectance include: 
! Subprogram not used     ! Wiscombe and Warren [1980] and Warren and Wiscombe [1980],
! Subprogram not used     ! Journal of Atmospheric Sciences, 37,
! Subprogram not used     !
! Subprogram not used     ! The multi-layer solution for multiple-scattering used here is from:
! Subprogram not used     ! Toon et al. [1989], Rapid calculation of radiative heating rates 
! Subprogram not used     ! and photodissociation rates in inhomogeneous multiple scattering atmospheres, 
! Subprogram not used     ! J. Geophys. Res., 94, D13, 16287-16301
! Subprogram not used     !
! Subprogram not used     ! The implementation of the SNICAR model in CLM/CSIM is described in:
! Subprogram not used     ! Flanner, M., C. Zender, J. Randerson, and P. Rasch [2007], 
! Subprogram not used     ! Present-day climate forcing and response from black carbon in snow,
! Subprogram not used     ! J. Geophys. Res., 112, D11202, doi: 10.1029/2006JD008003
! Subprogram not used 
! Subprogram not used     
! Subprogram not used     ! !USES:
! Subprogram not used     use clmtype
! Subprogram not used     use clm_varpar       , only : nlevsno, numrad
! Subprogram not used     use clm_time_manager , only : get_nstep
! Subprogram not used     use shr_const_mod    , only : SHR_CONST_PI
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: flg_snw_ice                          ! flag: =1 when called from CLM, =2 when called from CSIM
! Subprogram not used     integer , intent(in)  :: lbc, ubc                             ! column index bounds [unitless]
! Subprogram not used     integer , intent(in)  :: num_nourbanc                         ! number of columns in non-urban filter
! Subprogram not used     integer , intent(in)  :: filter_nourbanc(ubc-lbc+1)           ! column filter for non-urban points
! Subprogram not used     real(r8), intent(in)  :: coszen(lbc:ubc)                      ! cosine of solar zenith angle for next time step
! Subprogram not used                                                                   ! (col) [unitless]
! Subprogram not used     integer , intent(in)  :: flg_slr_in                           ! flag: =1 for direct-beam incident flux,
! Subprogram not used                                                                   ! =2 for diffuse incident flux
! Subprogram not used     real(r8), intent(in)  :: h2osno_liq(lbc:ubc,-nlevsno+1:0)     ! liquid water content (col,lyr) [kg/m2]
! Subprogram not used     real(r8), intent(in)  :: h2osno_ice(lbc:ubc,-nlevsno+1:0)     ! ice content (col,lyr) [kg/m2]
! Subprogram not used     integer,  intent(in)  :: snw_rds(lbc:ubc,-nlevsno+1:0)        ! snow effective radius (col,lyr) [microns, m^-6]
! Subprogram not used     real(r8), intent(in)  :: mss_cnc_aer_in(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)  ! mass concentration of all aerosol species
! Subprogram not used                                                                                ! (col,lyr,aer) [kg/kg]
! Subprogram not used     real(r8), intent(in)  :: albsfc(lbc:ubc,numrad)               ! albedo of surface underlying snow
! Subprogram not used                                                                   ! (col,bnd) [frc]
! Subprogram not used     real(r8), intent(out) :: albout(lbc:ubc,numrad)               ! snow albedo, averaged into 2 bands
! Subprogram not used                                                                   ! (=0 if no sun or no snow) (col,bnd) [frc]
! Subprogram not used     real(r8), intent(out) :: flx_abs(lbc:ubc,-nlevsno+1:1,numrad) ! absorbed flux in each layer per unit flux incident
! Subprogram not used                                                                   ! on top of snowpack (col,lyr,bnd) [frc]
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     !
! Subprogram not used     ! local pointers to implicit in arguments
! Subprogram not used     !
! Subprogram not used     integer,  pointer :: snl(:)              ! negative number of snow layers (col) [nbr]
! Subprogram not used     real(r8), pointer :: h2osno(:)           ! snow liquid water equivalent (col) [kg/m2]   
! Subprogram not used     integer,  pointer :: clandunit(:)        ! corresponding landunit of column (col) [idx] (debugging only)
! Subprogram not used     integer,  pointer :: cgridcell(:)        ! columns's gridcell index (col) [idx] (debugging only)
! Subprogram not used     integer,  pointer :: ltype(:)            ! landunit type (lnd) (debugging only)
! Subprogram not used     real(r8), pointer :: londeg(:)           ! longitude (degrees) (debugging only)
! Subprogram not used     real(r8), pointer :: latdeg(:)           ! latitude (degrees) (debugging only)
! Subprogram not used !
! Subprogram not used ! !OTHER LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! variables for snow radiative transfer calculations
! Subprogram not used 
! Subprogram not used     ! Local variables representing single-column values of arrays:
! Subprogram not used     integer :: snl_lcl                            ! negative number of snow layers [nbr]
! Subprogram not used     integer :: snw_rds_lcl(-nlevsno+1:0)          ! snow effective radius [m^-6]
! Subprogram not used     real(r8):: flx_slrd_lcl(1:numrad_snw)         ! direct beam incident irradiance [W/m2] (set to 1)
! Subprogram not used     real(r8):: flx_slri_lcl(1:numrad_snw)         ! diffuse incident irradiance [W/m2] (set to 1)
! Subprogram not used     real(r8):: mss_cnc_aer_lcl(-nlevsno+1:0,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
! Subprogram not used     real(r8):: h2osno_lcl                         ! total column snow mass [kg/m2]
! Subprogram not used     real(r8):: h2osno_liq_lcl(-nlevsno+1:0)       ! liquid water mass [kg/m2]
! Subprogram not used     real(r8):: h2osno_ice_lcl(-nlevsno+1:0)       ! ice mass [kg/m2]
! Subprogram not used     real(r8):: albsfc_lcl(1:numrad_snw)           ! albedo of underlying surface [frc]
! Subprogram not used     real(r8):: ss_alb_snw_lcl(-nlevsno+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
! Subprogram not used     real(r8):: asm_prm_snw_lcl(-nlevsno+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
! Subprogram not used     real(r8):: ext_cff_mss_snw_lcl(-nlevsno+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
! Subprogram not used     real(r8):: ss_alb_aer_lcl(sno_nbr_aer)        ! single-scatter albedo of aerosol species (aer_nbr) [frc] 
! Subprogram not used     real(r8):: asm_prm_aer_lcl(sno_nbr_aer)       ! asymmetry parameter of aerosol species (aer_nbr) [frc]
! Subprogram not used     real(r8):: ext_cff_mss_aer_lcl(sno_nbr_aer)   ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! Other local variables
! Subprogram not used     integer :: APRX_TYP                           ! two-stream approximation type
! Subprogram not used                                                   ! (1=Eddington, 2=Quadrature, 3=Hemispheric Mean) [nbr]
! Subprogram not used     integer :: DELTA                              ! flag to use Delta approximation (Joseph, 1976)
! Subprogram not used                                                   ! (1= use, 0= don't use)
! Subprogram not used     real(r8):: flx_wgt(1:numrad_snw)              ! weights applied to spectral bands,
! Subprogram not used                                                   ! specific to direct and diffuse cases (bnd) [frc]
! Subprogram not used    
! Subprogram not used     integer :: flg_nosnl                          ! flag: =1 if there is snow, but zero snow layers,
! Subprogram not used                                                   ! =0 if at least 1 snow layer [flg]   
! Subprogram not used     integer :: trip                               ! flag: =1 to redo RT calculation if result is unrealistic
! Subprogram not used     integer :: flg_dover                          ! defines conditions for RT redo (explained below)
! Subprogram not used 
! Subprogram not used     real(r8):: albedo                             ! temporary snow albedo [frc]
! Subprogram not used     real(r8):: flx_sum                            ! temporary summation variable for NIR weighting
! Subprogram not used     real(r8):: albout_lcl(numrad_snw)             ! snow albedo by band [frc]
! Subprogram not used     real(r8):: flx_abs_lcl(-nlevsno+1:1,numrad_snw)! absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]
! Subprogram not used  
! Subprogram not used     real(r8):: L_snw(-nlevsno+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
! Subprogram not used     real(r8):: tau_snw(-nlevsno+1:0)              ! snow optical depth (lyr) [unitless]
! Subprogram not used     real(r8):: L_aer(-nlevsno+1:0,sno_nbr_aer)    ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2] 
! Subprogram not used     real(r8):: tau_aer(-nlevsno+1:0,sno_nbr_aer)  ! aerosol optical depth (lyr,nbr_aer) [unitless]
! Subprogram not used     real(r8):: tau_sum                            ! cumulative (snow+aerosol) optical depth [unitless]
! Subprogram not used     real(r8):: tau_clm(-nlevsno+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless] 
! Subprogram not used     real(r8):: omega_sum                          ! temporary summation of single-scatter albedo of all aerosols [frc]
! Subprogram not used     real(r8):: g_sum                              ! temporary summation of asymmetry parameter of all aerosols [frc]
! Subprogram not used 
! Subprogram not used     real(r8):: tau(-nlevsno+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
! Subprogram not used     real(r8):: omega(-nlevsno+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
! Subprogram not used     real(r8):: g(-nlevsno+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
! Subprogram not used     real(r8):: tau_star(-nlevsno+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
! Subprogram not used                                                   ! (lyr) [unitless]
! Subprogram not used     real(r8):: omega_star(-nlevsno+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
! Subprogram not used     real(r8):: g_star(-nlevsno+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
! Subprogram not used                                                   ! (lyr) [frc]
! Subprogram not used    
! Subprogram not used     integer :: nstep                              ! current timestep [nbr] (debugging only)
! Subprogram not used     integer :: g_idx, c_idx, l_idx                ! gridcell, column, and landunit indices [idx]
! Subprogram not used     integer :: bnd_idx                            ! spectral band index (1 <= bnd_idx <= numrad_snw) [idx]
! Subprogram not used     integer :: rds_idx                            ! snow effective radius index for retrieving
! Subprogram not used                                                   ! Mie parameters from lookup table [idx]
! Subprogram not used     integer :: snl_btm                            ! index of bottom snow layer (0) [idx]
! Subprogram not used     integer :: snl_top                            ! index of top snow layer (-4 to 0) [idx]
! Subprogram not used     integer :: fc                                 ! column filter index
! Subprogram not used     integer :: i                                  ! layer index [idx]
! Subprogram not used     integer :: j                                  ! aerosol number index [idx]
! Subprogram not used     integer :: n                                  ! tridiagonal matrix index [idx]
! Subprogram not used     integer :: m                                  ! secondary layer index [idx]
! Subprogram not used    
! Subprogram not used     real(r8):: F_direct(-nlevsno+1:0)             ! direct-beam radiation at bottom of layer interface (lyr) [W/m^2]
! Subprogram not used     real(r8):: F_net(-nlevsno+1:0)                ! net radiative flux at bottom of layer interface (lyr) [W/m^2]
! Subprogram not used     real(r8):: F_abs(-nlevsno+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
! Subprogram not used     real(r8):: F_abs_sum                          ! total absorbed energy in column [W/m^2]
! Subprogram not used     real(r8):: F_sfc_pls                          ! upward radiative flux at snowpack top [W/m^2]
! Subprogram not used     real(r8):: F_btm_net                          ! net flux at bottom of snowpack [W/m^2]                    
! Subprogram not used     real(r8):: F_sfc_net                          ! net flux at top of snowpack [W/m^2]
! Subprogram not used     real(r8):: energy_sum                         ! sum of all energy terms; should be 0.0 [W/m^2]
! Subprogram not used     real(r8):: F_direct_btm                       ! direct-beam radiation at bottom of snowpack [W/m^2]
! Subprogram not used     real(r8):: mu_not                             ! cosine of solar zenith angle (used locally) [frc]
! Subprogram not used 
! Subprogram not used     integer :: err_idx                            ! counter for number of times through error loop [nbr]
! Subprogram not used     real(r8):: lat_coord                          ! gridcell latitude (debugging only)
! Subprogram not used     real(r8):: lon_coord                          ! gridcell longitude (debugging only)
! Subprogram not used     integer :: sfctype                            ! underlying surface type (debugging only)
! Subprogram not used     real(r8):: pi                                 ! 3.1415...
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! intermediate variables for radiative transfer approximation:
! Subprogram not used     real(r8):: gamma1(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: gamma2(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: gamma3(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: gamma4(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: lambda(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: GAMMA(-nlevsno+1:0)                ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: mu_one                             ! two-stream coefficient from Toon et al. (lyr) [unitless]
! Subprogram not used     real(r8):: e1(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
! Subprogram not used     real(r8):: e2(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
! Subprogram not used     real(r8):: e3(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
! Subprogram not used     real(r8):: e4(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
! Subprogram not used     real(r8):: C_pls_btm(-nlevsno+1:0)            ! intermediate variable: upward flux at bottom interface (lyr) [W/m2]
! Subprogram not used     real(r8):: C_mns_btm(-nlevsno+1:0)            ! intermediate variable: downward flux at bottom interface (lyr) [W/m2]
! Subprogram not used     real(r8):: C_pls_top(-nlevsno+1:0)            ! intermediate variable: upward flux at top interface (lyr) [W/m2]
! Subprogram not used     real(r8):: C_mns_top(-nlevsno+1:0)            ! intermediate variable: downward flux at top interface (lyr) [W/m2]
! Subprogram not used     real(r8):: A(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: B(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: D(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: E(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: AS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: DS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: X(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used     real(r8):: Y(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! Assign local pointers to derived subtypes components (column-level)
! Subprogram not used     ! (CLM-specific)
! Subprogram not used     if (flg_snw_ice == 1) then
! Subprogram not used        snl            => cps%snl
! Subprogram not used        h2osno         => cws%h2osno
! Subprogram not used        clandunit      => col%landunit  ! (debug only)
! Subprogram not used        cgridcell      => col%gridcell  ! (debug only)
! Subprogram not used        ltype          => lun%itype       ! (debug only)
! Subprogram not used        londeg         => grc%londeg        ! (debug only)
! Subprogram not used        latdeg         => grc%latdeg        ! (debug only)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! Define constants
! Subprogram not used     pi = SHR_CONST_PI
! Subprogram not used 
! Subprogram not used     ! always use Delta approximation for snow
! Subprogram not used     DELTA = 1
! Subprogram not used 
! Subprogram not used     ! Get current timestep
! Subprogram not used     nstep = get_nstep()
! Subprogram not used 
! Subprogram not used     ! Loop over all non-urban columns
! Subprogram not used     ! (when called from CSIM, there is only one column)
! Subprogram not used     do fc = 1,num_nourbanc
! Subprogram not used        c_idx = filter_nourbanc(fc)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! Zero absorbed radiative fluxes:
! Subprogram not used        do i=-nlevsno+1,1,1
! Subprogram not used           flx_abs_lcl(:,:)   = 0._r8
! Subprogram not used           flx_abs(c_idx,i,:) = 0._r8
! Subprogram not used        enddo
! Subprogram not used        
! Subprogram not used        ! set snow/ice mass to be used for RT:
! Subprogram not used        if (flg_snw_ice == 1) then
! Subprogram not used           h2osno_lcl = h2osno(c_idx)
! Subprogram not used        else
! Subprogram not used           h2osno_lcl = h2osno_ice(c_idx,0)
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! Qualifier for computing snow RT: 
! Subprogram not used        !  1) sunlight from atmosphere model 
! Subprogram not used        !  2) minimum amount of snow on ground. 
! Subprogram not used        !     Otherwise, set snow albedo to zero
! Subprogram not used        if ((coszen(c_idx) > 0._r8) .and. (h2osno_lcl > min_snw)) then     
! Subprogram not used 
! Subprogram not used           ! Set variables specific to CLM
! Subprogram not used           if (flg_snw_ice == 1) then
! Subprogram not used              ! Assign local (single-column) variables to global values
! Subprogram not used              ! If there is snow, but zero snow layers, we must create a layer locally.
! Subprogram not used              ! This layer is presumed to have the fresh snow effective radius.
! Subprogram not used              if (snl(c_idx) > -1) then
! Subprogram not used                 flg_nosnl         =  1
! Subprogram not used                 snl_lcl           =  -1
! Subprogram not used                 h2osno_ice_lcl(0) =  h2osno_lcl
! Subprogram not used                 h2osno_liq_lcl(0) =  0._r8
! Subprogram not used                 snw_rds_lcl(0)    =  nint(snw_rds_min)
! Subprogram not used              else
! Subprogram not used                 flg_nosnl         =  0
! Subprogram not used                 snl_lcl           =  snl(c_idx)
! Subprogram not used                 h2osno_liq_lcl(:) =  h2osno_liq(c_idx,:)
! Subprogram not used                 h2osno_ice_lcl(:) =  h2osno_ice(c_idx,:)
! Subprogram not used                 snw_rds_lcl(:)    =  snw_rds(c_idx,:)
! Subprogram not used              endif
! Subprogram not used             
! Subprogram not used              snl_btm   = 0
! Subprogram not used              snl_top   = snl_lcl+1
! Subprogram not used 
! Subprogram not used              ! for debugging only
! Subprogram not used              l_idx     = clandunit(c_idx)
! Subprogram not used              g_idx     = cgridcell(c_idx)
! Subprogram not used              sfctype   = ltype(l_idx)
! Subprogram not used              lat_coord = latdeg(g_idx)
! Subprogram not used              lon_coord = londeg(g_idx)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           ! Set variables specific to CSIM
! Subprogram not used           else
! Subprogram not used              flg_nosnl         = 0
! Subprogram not used              snl_lcl           = -1
! Subprogram not used              h2osno_liq_lcl(:) = h2osno_liq(c_idx,:)
! Subprogram not used              h2osno_ice_lcl(:) = h2osno_ice(c_idx,:)
! Subprogram not used              snw_rds_lcl(:)    = snw_rds(c_idx,:)
! Subprogram not used              snl_btm           = 0
! Subprogram not used              snl_top           = 0
! Subprogram not used              sfctype           = -1
! Subprogram not used              lat_coord         = -90
! Subprogram not used              lon_coord         = 0
! Subprogram not used           endif
! Subprogram not used 
! Subprogram not used           ! Set local aerosol array
! Subprogram not used           do j=1,sno_nbr_aer
! Subprogram not used              mss_cnc_aer_lcl(:,j) = mss_cnc_aer_in(c_idx,:,j)
! Subprogram not used           enddo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           ! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
! Subprogram not used           albsfc_lcl(1)                       = albsfc(c_idx,1)
! Subprogram not used           albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(c_idx,2)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           ! Error check for snow grain size:
! Subprogram not used           do i=snl_top,snl_btm,1
! Subprogram not used              if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
! Subprogram not used                 write (iulog,*) "SNICAR ERROR: snow grain radius of ", snw_rds_lcl(i), " out of bounds."
! Subprogram not used                 write (iulog,*) "NSTEP= ", nstep
! Subprogram not used                 write (iulog,*) "flg_snw_ice= ", flg_snw_ice
! Subprogram not used                 write (iulog,*) "column: ", c_idx, " level: ", i, " snl(c)= ", snl_lcl
! Subprogram not used                 write (iulog,*) "lat= ", lat_coord, " lon= ", lon_coord
! Subprogram not used                 write (iulog,*) "h2osno(c)= ", h2osno_lcl
! Subprogram not used                 call endrun()
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used 
! Subprogram not used           ! Incident flux weighting parameters
! Subprogram not used           !  - sum of all VIS bands must equal 1
! Subprogram not used           !  - sum of all NIR bands must equal 1
! Subprogram not used           !
! Subprogram not used           ! Spectral bands (5-band case)
! Subprogram not used           !  Band 1: 0.3-0.7um (VIS)
! Subprogram not used           !  Band 2: 0.7-1.0um (NIR)
! Subprogram not used           !  Band 3: 1.0-1.2um (NIR)
! Subprogram not used           !  Band 4: 1.2-1.5um (NIR)
! Subprogram not used           !  Band 5: 1.5-5.0um (NIR)
! Subprogram not used           !
! Subprogram not used           ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
! Subprogram not used           !
! Subprogram not used           ! 3-band weights
! Subprogram not used           if (numrad_snw==3) then
! Subprogram not used              ! Direct:
! Subprogram not used              if (flg_slr_in == 1) then
! Subprogram not used                 flx_wgt(1) = 1._r8
! Subprogram not used                 flx_wgt(2) = 0.66628670195247_r8
! Subprogram not used                 flx_wgt(3) = 0.33371329804753_r8
! Subprogram not used              ! Diffuse:
! Subprogram not used              elseif (flg_slr_in == 2) then
! Subprogram not used                 flx_wgt(1) = 1._r8
! Subprogram not used                 flx_wgt(2) = 0.77887652162877_r8
! Subprogram not used                 flx_wgt(3) = 0.22112347837123_r8
! Subprogram not used              endif
! Subprogram not used 
! Subprogram not used           ! 5-band weights
! Subprogram not used           elseif(numrad_snw==5) then
! Subprogram not used              ! Direct:
! Subprogram not used              if (flg_slr_in == 1) then
! Subprogram not used                 flx_wgt(1) = 1._r8
! Subprogram not used                 flx_wgt(2) = 0.49352158521175_r8
! Subprogram not used                 flx_wgt(3) = 0.18099494230665_r8
! Subprogram not used                 flx_wgt(4) = 0.12094898498813_r8
! Subprogram not used                 flx_wgt(5) = 0.20453448749347_r8
! Subprogram not used              ! Diffuse:
! Subprogram not used              elseif (flg_slr_in == 2) then
! Subprogram not used                 flx_wgt(1) = 1._r8
! Subprogram not used                 flx_wgt(2) = 0.58581507618433_r8
! Subprogram not used                 flx_wgt(3) = 0.20156903770812_r8
! Subprogram not used                 flx_wgt(4) = 0.10917889346386_r8
! Subprogram not used                 flx_wgt(5) = 0.10343699264369_r8
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used 
! Subprogram not used           ! Loop over snow spectral bands
! Subprogram not used           do bnd_idx = 1,numrad_snw
! Subprogram not used 
! Subprogram not used              mu_not    = coszen(c_idx)  ! must set here, because of error handling
! Subprogram not used              flg_dover = 1              ! default is to redo
! Subprogram not used              err_idx   = 0              ! number of times through loop
! Subprogram not used 
! Subprogram not used              do while (flg_dover > 0)
! Subprogram not used 
! Subprogram not used                 ! DEFAULT APPROXIMATIONS:
! Subprogram not used                 !  VIS:       Delta-Eddington
! Subprogram not used                 !  NIR (all): Delta-Hemispheric Mean
! Subprogram not used                 !  WARNING:   DO NOT USE DELTA-EDDINGTON FOR NIR DIFFUSE - this sometimes results in negative albedo
! Subprogram not used                 !  
! Subprogram not used                 ! ERROR CONDITIONS:
! Subprogram not used                 !  Conditions which cause "trip", resulting in redo of RT approximation:
! Subprogram not used                 !   1. negative absorbed flux
! Subprogram not used                 !   2. total absorbed flux greater than incident flux
! Subprogram not used                 !   3. negative albedo
! Subprogram not used                 !   NOTE: These errors have only been encountered in spectral bands 4 and 5
! Subprogram not used                 !
! Subprogram not used                 ! ERROR HANDLING
! Subprogram not used                 !  1st error (flg_dover=2): switch approximation (Edd->HM or HM->Edd)
! Subprogram not used                 !  2nd error (flg_dover=3): change zenith angle by 0.02 (this happens about 1 in 10^6 cases)
! Subprogram not used                 !  3rd error (flg_dover=4): switch approximation with new zenith
! Subprogram not used                 !  Subsequent errors: repeatedly change zenith and approximations...
! Subprogram not used               
! Subprogram not used                 if (bnd_idx == 1) then
! Subprogram not used                    if (flg_dover == 2) then
! Subprogram not used                       APRX_TYP = 3
! Subprogram not used                    elseif (flg_dover == 3) then
! Subprogram not used                       APRX_TYP = 1
! Subprogram not used                       if (coszen(c_idx) > 0.5_r8) then
! Subprogram not used                          mu_not = mu_not - 0.02_r8
! Subprogram not used                       else
! Subprogram not used                          mu_not = mu_not + 0.02_r8
! Subprogram not used                       endif
! Subprogram not used                    elseif (flg_dover == 4) then
! Subprogram not used                       APRX_TYP = 3
! Subprogram not used                    else
! Subprogram not used                       APRX_TYP = 1
! Subprogram not used                    endif
! Subprogram not used                    
! Subprogram not used                 else
! Subprogram not used                    if (flg_dover == 2) then
! Subprogram not used                       APRX_TYP = 1
! Subprogram not used                    elseif (flg_dover == 3) then
! Subprogram not used                       APRX_TYP = 3
! Subprogram not used                       if (coszen(c_idx) > 0.5_r8) then
! Subprogram not used                          mu_not = mu_not - 0.02_r8
! Subprogram not used                       else
! Subprogram not used                          mu_not = mu_not + 0.02_r8
! Subprogram not used                       endif
! Subprogram not used                    elseif (flg_dover == 4) then
! Subprogram not used                       APRX_TYP = 1
! Subprogram not used                    else
! Subprogram not used                       APRX_TYP = 3
! Subprogram not used                    endif
! Subprogram not used 
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 ! Set direct or diffuse incident irradiance to 1
! Subprogram not used                 ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
! Subprogram not used                 if (flg_slr_in == 1) then
! Subprogram not used                    flx_slrd_lcl(bnd_idx) = 1._r8/(mu_not*pi) ! this corresponds to incident irradiance of 1.0
! Subprogram not used                    flx_slri_lcl(bnd_idx) = 0._r8
! Subprogram not used                 else
! Subprogram not used                    flx_slrd_lcl(bnd_idx) = 0._r8
! Subprogram not used                    flx_slri_lcl(bnd_idx) = 1._r8
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
! Subprogram not used                 ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
! Subprogram not used                 if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
! Subprogram not used                    mss_cnc_aer_lcl(:,:) = 0._r8
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 if ( (numrad_snw == 3).and.(bnd_idx == 3) ) then
! Subprogram not used                    mss_cnc_aer_lcl(:,:) = 0._r8
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 ! Define local Mie parameters based on snow grain size and aerosol species,
! Subprogram not used                 !  retrieved from a lookup table.
! Subprogram not used                 if (flg_slr_in == 1) then
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
! Subprogram not used                       ! snow optical properties (direct radiation)
! Subprogram not used                       ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
! Subprogram not used                       asm_prm_snw_lcl(i)     = asm_prm_snw_drc(rds_idx,bnd_idx)
! Subprogram not used                       ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
! Subprogram not used                    enddo
! Subprogram not used                 elseif (flg_slr_in == 2) then
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
! Subprogram not used                       ! snow optical properties (diffuse radiation)
! Subprogram not used                       ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
! Subprogram not used                       asm_prm_snw_lcl(i)     = asm_prm_snw_dfs(rds_idx,bnd_idx)
! Subprogram not used                       ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
! Subprogram not used                    enddo
! Subprogram not used                 endif
! Subprogram not used                    
! Subprogram not used                 ! aerosol species 1 optical properties
! Subprogram not used                 ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)
! Subprogram not used                 
! Subprogram not used                 ! aerosol species 2 optical properties
! Subprogram not used                 ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)
! Subprogram not used                 
! Subprogram not used                 ! aerosol species 3 optical properties
! Subprogram not used                 ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)
! Subprogram not used                 
! Subprogram not used                 ! aerosol species 4 optical properties
! Subprogram not used                 ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)
! Subprogram not used 
! Subprogram not used                 ! aerosol species 5 optical properties
! Subprogram not used                 ss_alb_aer_lcl(5)        = ss_alb_dst1(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(5)       = asm_prm_dst1(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(5)   = ext_cff_mss_dst1(bnd_idx)
! Subprogram not used                 
! Subprogram not used                 ! aerosol species 6 optical properties
! Subprogram not used                 ss_alb_aer_lcl(6)        = ss_alb_dst2(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(6)       = asm_prm_dst2(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(6)   = ext_cff_mss_dst2(bnd_idx)
! Subprogram not used                 
! Subprogram not used                 ! aerosol species 7 optical properties
! Subprogram not used                 ss_alb_aer_lcl(7)        = ss_alb_dst3(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(7)       = asm_prm_dst3(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(7)   = ext_cff_mss_dst3(bnd_idx)
! Subprogram not used                 
! Subprogram not used                 ! aerosol species 8 optical properties
! Subprogram not used                 ss_alb_aer_lcl(8)        = ss_alb_dst4(bnd_idx)      
! Subprogram not used                 asm_prm_aer_lcl(8)       = asm_prm_dst4(bnd_idx)
! Subprogram not used                 ext_cff_mss_aer_lcl(8)   = ext_cff_mss_dst4(bnd_idx)
! Subprogram not used                 
! Subprogram not used 
! Subprogram not used                 ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
! Subprogram not used                 ! 2. optical Depths (tau_snw, tau_aer)
! Subprogram not used                 ! 3. weighted Mie properties (tau, omega, g)
! Subprogram not used 
! Subprogram not used                 ! Weighted Mie parameters of each layer
! Subprogram not used                 do i=snl_top,snl_btm,1
! Subprogram not used                    L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
! Subprogram not used                    tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)
! Subprogram not used 
! Subprogram not used                    do j=1,sno_nbr_aer
! Subprogram not used                       L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
! Subprogram not used                       tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
! Subprogram not used                    enddo
! Subprogram not used 
! Subprogram not used                    tau_sum   = 0._r8
! Subprogram not used                    omega_sum = 0._r8
! Subprogram not used                    g_sum     = 0._r8
! Subprogram not used   
! Subprogram not used                    do j=1,sno_nbr_aer
! Subprogram not used                       tau_sum    = tau_sum + tau_aer(i,j) 
! Subprogram not used                       omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
! Subprogram not used                       g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
! Subprogram not used                    enddo
! Subprogram not used 
! Subprogram not used                    tau(i)    = tau_sum + tau_snw(i)
! Subprogram not used                    omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
! Subprogram not used                    g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 ! DELTA transformations, if requested
! Subprogram not used                 if (DELTA == 1) then
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       g_star(i)     = g(i)/(1+g(i))
! Subprogram not used                       omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
! Subprogram not used                       tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
! Subprogram not used                    enddo
! Subprogram not used                 else
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       g_star(i)     = g(i)
! Subprogram not used                       omega_star(i) = omega(i)
! Subprogram not used                       tau_star(i)   = tau(i)
! Subprogram not used                    enddo
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 ! Total column optical depth:
! Subprogram not used                 ! tau_clm(i) = total optical depth above the bottom of layer i
! Subprogram not used                 tau_clm(snl_top) = 0._r8
! Subprogram not used                 do i=snl_top+1,snl_btm,1
! Subprogram not used                    tau_clm(i) = tau_clm(i-1)+tau_star(i-1)
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 ! Direct radiation at bottom of snowpack:
! Subprogram not used                 F_direct_btm = albsfc_lcl(bnd_idx)*mu_not*exp(-(tau_clm(snl_btm)+tau_star(snl_btm))/mu_not)*pi*flx_slrd_lcl(bnd_idx)
! Subprogram not used 
! Subprogram not used                 ! Intermediates
! Subprogram not used                 ! Gamma values are approximation-specific.
! Subprogram not used 
! Subprogram not used                 ! Eddington
! Subprogram not used                 if (APRX_TYP==1) then
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       gamma1(i) = (7-(omega_star(i)*(4+(3*g_star(i)))))/4
! Subprogram not used                       gamma2(i) = -(1-(omega_star(i)*(4-(3*g_star(i)))))/4
! Subprogram not used                       gamma3(i) = (2-(3*g_star(i)*mu_not))/4
! Subprogram not used                       gamma4(i) = 1-gamma3(i)
! Subprogram not used                       mu_one    = 0.5
! Subprogram not used                    enddo
! Subprogram not used                    
! Subprogram not used                 ! Quadrature
! Subprogram not used                 elseif (APRX_TYP==2) then
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       gamma1(i) = (3**0.5)*(2-(omega_star(i)*(1+g_star(i))))/2
! Subprogram not used                       gamma2(i) = omega_star(i)*(3**0.5)*(1-g_star(i))/2
! Subprogram not used                       gamma3(i) = (1-((3**0.5)*g_star(i)*mu_not))/2
! Subprogram not used                       gamma4(i) = 1-gamma3(i)
! Subprogram not used                       mu_one    = 1/(3**0.5)
! Subprogram not used                    enddo
! Subprogram not used 
! Subprogram not used                 ! Hemispheric Mean
! Subprogram not used                 elseif (APRX_TYP==3) then
! Subprogram not used                    do i=snl_top,snl_btm,1
! Subprogram not used                       gamma1(i) = 2 - (omega_star(i)*(1+g_star(i)))
! Subprogram not used                       gamma2(i) = omega_star(i)*(1-g_star(i))
! Subprogram not used                       gamma3(i) = (1-((3**0.5)*g_star(i)*mu_not))/2
! Subprogram not used                       gamma4(i) = 1-gamma3(i)
! Subprogram not used                       mu_one    = 0.5
! Subprogram not used                    enddo
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 ! Intermediates for tri-diagonal solution
! Subprogram not used                 do i=snl_top,snl_btm,1
! Subprogram not used                    lambda(i) = sqrt(abs((gamma1(i)**2) - (gamma2(i)**2)))
! Subprogram not used                    GAMMA(i)  = gamma2(i)/(gamma1(i)+lambda(i))
! Subprogram not used 
! Subprogram not used                    e1(i)     = 1+(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
! Subprogram not used                    e2(i)     = 1-(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
! Subprogram not used                    e3(i)     = GAMMA(i) + exp(-lambda(i)*tau_star(i))
! Subprogram not used                    e4(i)     = GAMMA(i) - exp(-lambda(i)*tau_star(i))
! Subprogram not used                 enddo !enddo over snow layers
! Subprogram not used 
! Subprogram not used 
! Subprogram not used                 ! Intermediates for tri-diagonal solution
! Subprogram not used                 do i=snl_top,snl_btm,1
! Subprogram not used                    if (flg_slr_in == 1) then
! Subprogram not used 
! Subprogram not used                       C_pls_btm(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
! Subprogram not used                                      exp(-(tau_clm(i)+tau_star(i))/mu_not)*   &
! Subprogram not used                                      (((gamma1(i)-(1/mu_not))*gamma3(i))+     &
! Subprogram not used                                      (gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))
! Subprogram not used 
! Subprogram not used                       C_mns_btm(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
! Subprogram not used                                      exp(-(tau_clm(i)+tau_star(i))/mu_not)*   &
! Subprogram not used                                      (((gamma1(i)+(1/mu_not))*gamma4(i))+     &
! Subprogram not used                                      (gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))
! Subprogram not used 
! Subprogram not used                       C_pls_top(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
! Subprogram not used                                      exp(-tau_clm(i)/mu_not)*(((gamma1(i)-(1/mu_not))* &
! Subprogram not used                                      gamma3(i))+(gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))
! Subprogram not used 
! Subprogram not used                       C_mns_top(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
! Subprogram not used                                      exp(-tau_clm(i)/mu_not)*(((gamma1(i)+(1/mu_not))* &
! Subprogram not used                                      gamma4(i))+(gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))
! Subprogram not used 
! Subprogram not used                    else
! Subprogram not used                       C_pls_btm(i) = 0._r8
! Subprogram not used                       C_mns_btm(i) = 0._r8
! Subprogram not used                       C_pls_top(i) = 0._r8
! Subprogram not used                       C_mns_top(i) = 0._r8
! Subprogram not used                    endif
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 ! Coefficients for tridiaganol matrix solution
! Subprogram not used                 do i=2*snl_lcl+1,0,1
! Subprogram not used 
! Subprogram not used                    !Boundary values for i=1 and i=2*snl_lcl, specifics for i=odd and i=even    
! Subprogram not used                    if (i==(2*snl_lcl+1)) then
! Subprogram not used                       A(i) = 0
! Subprogram not used                       B(i) = e1(snl_top)
! Subprogram not used                       D(i) = -e2(snl_top)
! Subprogram not used                       E(i) = flx_slri_lcl(bnd_idx)-C_mns_top(snl_top)
! Subprogram not used 
! Subprogram not used                    elseif(i==0) then
! Subprogram not used                       A(i) = e1(snl_btm)-(albsfc_lcl(bnd_idx)*e3(snl_btm))
! Subprogram not used                       B(i) = e2(snl_btm)-(albsfc_lcl(bnd_idx)*e4(snl_btm))
! Subprogram not used                       D(i) = 0
! Subprogram not used                       E(i) = F_direct_btm-C_pls_btm(snl_btm)+(albsfc_lcl(bnd_idx)*C_mns_btm(snl_btm))
! Subprogram not used 
! Subprogram not used                    elseif(mod(i,2)==-1) then   ! If odd and i>=3 (n=1 for i=3)
! Subprogram not used                       n=floor(i/2.0)
! Subprogram not used                       A(i) = (e2(n)*e3(n))-(e4(n)*e1(n))
! Subprogram not used                       B(i) = (e1(n)*e1(n+1))-(e3(n)*e3(n+1))
! Subprogram not used                       D(i) = (e3(n)*e4(n+1))-(e1(n)*e2(n+1))
! Subprogram not used                       E(i) = (e3(n)*(C_pls_top(n+1)-C_pls_btm(n)))+(e1(n)*(C_mns_btm(n)-C_mns_top(n+1)))
! Subprogram not used 
! Subprogram not used                    elseif(mod(i,2)==0) then    ! If even and i<=2*snl_lcl
! Subprogram not used                       n=(i/2)
! Subprogram not used                       A(i) = (e2(n+1)*e1(n))-(e3(n)*e4(n+1))
! Subprogram not used                       B(i) = (e2(n)*e2(n+1))-(e4(n)*e4(n+1))
! Subprogram not used                       D(i) = (e1(n+1)*e4(n+1))-(e2(n+1)*e3(n+1))
! Subprogram not used                       E(i) = (e2(n+1)*(C_pls_top(n+1)-C_pls_btm(n)))+(e4(n+1)*(C_mns_top(n+1)-C_mns_btm(n))) 
! Subprogram not used                    endif
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 AS(0) = A(0)/B(0)
! Subprogram not used                 DS(0) = E(0)/B(0)
! Subprogram not used 
! Subprogram not used                 do i=-1,(2*snl_lcl+1),-1
! Subprogram not used                    X(i)  = 1/(B(i)-(D(i)*AS(i+1)))
! Subprogram not used                    AS(i) = A(i)*X(i)
! Subprogram not used                    DS(i) = (E(i)-(D(i)*DS(i+1)))*X(i)
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 Y(2*snl_lcl+1) = DS(2*snl_lcl+1)
! Subprogram not used                 do i=(2*snl_lcl+2),0,1
! Subprogram not used                    Y(i) = DS(i)-(AS(i)*Y(i-1))
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 ! Downward direct-beam and net flux (F_net) at the base of each layer:
! Subprogram not used                 do i=snl_top,snl_btm,1
! Subprogram not used                    F_direct(i) = mu_not*pi*flx_slrd_lcl(bnd_idx)*exp(-(tau_clm(i)+tau_star(i))/mu_not)
! Subprogram not used                    F_net(i)    = (Y(2*i-1)*(e1(i)-e3(i))) + (Y(2*i)*(e2(i)-e4(i))) + &
! Subprogram not used                                  C_pls_btm(i) - C_mns_btm(i) - F_direct(i)
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 ! Upward flux at snowpack top:
! Subprogram not used                 F_sfc_pls = (Y(2*snl_lcl+1)*(exp(-lambda(snl_top)*tau_star(snl_top))+ &
! Subprogram not used                             GAMMA(snl_top))) + (Y(2*snl_lcl+2)*(exp(-lambda(snl_top)* &
! Subprogram not used                             tau_star(snl_top))-GAMMA(snl_top))) + C_pls_top(snl_top)
! Subprogram not used 
! Subprogram not used                 ! Net flux at bottom = absorbed radiation by underlying surface:
! Subprogram not used                 F_btm_net = -F_net(snl_btm)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used                 ! Bulk column albedo and surface net flux
! Subprogram not used                 albedo    = F_sfc_pls/((mu_not*pi*flx_slrd_lcl(bnd_idx))+flx_slri_lcl(bnd_idx))
! Subprogram not used                 F_sfc_net = F_sfc_pls - ((mu_not*pi*flx_slrd_lcl(bnd_idx))+flx_slri_lcl(bnd_idx))
! Subprogram not used 
! Subprogram not used                 trip = 0
! Subprogram not used                 ! Absorbed flux in each layer
! Subprogram not used                 do i=snl_top,snl_btm,1
! Subprogram not used                    if(i==snl_top) then
! Subprogram not used                       F_abs(i) = F_net(i)-F_sfc_net
! Subprogram not used                    else
! Subprogram not used                       F_abs(i) = F_net(i)-F_net(i-1)
! Subprogram not used                    endif
! Subprogram not used                    flx_abs_lcl(i,bnd_idx) = F_abs(i)
! Subprogram not used                    
! Subprogram not used 
! Subprogram not used                    ! ERROR check: negative absorption
! Subprogram not used                    if (flx_abs_lcl(i,bnd_idx) < -0.00001) then
! Subprogram not used                       trip = 1
! Subprogram not used                    endif
! Subprogram not used                 enddo
! Subprogram not used                 
! Subprogram not used                 flx_abs_lcl(1,bnd_idx) = F_btm_net
! Subprogram not used              
! Subprogram not used                 if (flg_nosnl == 1) then
! Subprogram not used                    ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
! Subprogram not used                    !flx_abs_lcl(:,bnd_idx) = 0._r8
! Subprogram not used                    !flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net
! Subprogram not used 
! Subprogram not used                    ! changed on 20070408:
! Subprogram not used                    ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
! Subprogram not used                    ! handles the case of no snow layers. Then, if a snow layer is addded between now and
! Subprogram not used                    ! SurfaceRadiation (called in Hydrology1), absorbed energy will be properly distributed.
! Subprogram not used                    flx_abs_lcl(0,bnd_idx) = F_abs(0)
! Subprogram not used                    flx_abs_lcl(1,bnd_idx) = F_btm_net
! Subprogram not used                 endif
! Subprogram not used                 
! Subprogram not used                 !Underflow check (we've already tripped the error condition above)
! Subprogram not used                 do i=snl_top,1,1
! Subprogram not used                    if (flx_abs_lcl(i,bnd_idx) < 0._r8) then
! Subprogram not used                       flx_abs_lcl(i,bnd_idx) = 0._r8
! Subprogram not used                    endif
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used                 F_abs_sum = 0._r8
! Subprogram not used                 do i=snl_top,snl_btm,1
! Subprogram not used                    F_abs_sum = F_abs_sum + F_abs(i)
! Subprogram not used                 enddo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used                 !ERROR check: absorption greater than incident flux
! Subprogram not used                 ! (should make condition more generic than "1._r8")
! Subprogram not used                 if (F_abs_sum > 1._r8) then
! Subprogram not used                    trip = 1
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used                 !ERROR check:
! Subprogram not used                 if ((albedo < 0._r8).and.(trip==0)) then
! Subprogram not used                    trip = 1
! Subprogram not used                 endif
! Subprogram not used                 
! Subprogram not used                 ! Set conditions for redoing RT calculation 
! Subprogram not used                 if ((trip == 1).and.(flg_dover == 1)) then
! Subprogram not used                    flg_dover = 2
! Subprogram not used                 elseif ((trip == 1).and.(flg_dover == 2)) then
! Subprogram not used                    flg_dover = 3
! Subprogram not used                 elseif ((trip == 1).and.(flg_dover == 3)) then
! Subprogram not used                    flg_dover = 4
! Subprogram not used                 elseif((trip == 1).and.(flg_dover == 4).and.(err_idx < 20)) then
! Subprogram not used                    flg_dover = 3
! Subprogram not used                    err_idx = err_idx + 1
! Subprogram not used                 elseif((trip == 1).and.(flg_dover == 4).and.(err_idx >= 20)) then
! Subprogram not used                    flg_dover = 0
! Subprogram not used                    write(iulog,*) "SNICAR ERROR: FOUND A WORMHOLE. STUCK IN INFINITE LOOP! Called from: ", flg_snw_ice
! Subprogram not used                    write(iulog,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
! Subprogram not used                    write(iulog,*) "SNICAR STATS: soot1(0)= ", mss_cnc_aer_lcl(0,1)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: soot2(0)= ", mss_cnc_aer_lcl(0,2)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: dust1(0)= ", mss_cnc_aer_lcl(0,3)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: dust2(0)= ", mss_cnc_aer_lcl(0,4)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: dust3(0)= ", mss_cnc_aer_lcl(0,5)
! Subprogram not used                    write(iulog,*) "SNICAR STATS: dust4(0)= ", mss_cnc_aer_lcl(0,6)
! Subprogram not used                   
! Subprogram not used                    call endrun()
! Subprogram not used                 else
! Subprogram not used                    flg_dover = 0
! Subprogram not used                 endif
! Subprogram not used 
! Subprogram not used              enddo !enddo while (flg_dover > 0)
! Subprogram not used              
! Subprogram not used              ! Energy conservation check:
! Subprogram not used              ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
! Subprogram not used              energy_sum = (mu_not*pi*flx_slrd_lcl(bnd_idx)) + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)
! Subprogram not used              if (abs(energy_sum) > 0.00001_r8) then
! Subprogram not used                 write (iulog,"(a,e12.6,a,i6,a,i6)") "SNICAR ERROR: Energy conservation error of : ", energy_sum, &
! Subprogram not used                              " at timestep: ", nstep, " at column: ", c_idx
! Subprogram not used                 call endrun()
! Subprogram not used              endif
! Subprogram not used 
! Subprogram not used              albout_lcl(bnd_idx) = albedo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used              ! Check that albedo is less than 1
! Subprogram not used              if (albout_lcl(bnd_idx) > 1.0) then
! Subprogram not used 
! Subprogram not used                 write (iulog,*) "SNICAR ERROR: Albedo > 1.0 at c: ", c_idx, " NSTEP= ",nstep
! Subprogram not used                 write (iulog,*) "SNICAR STATS: bnd_idx= ",bnd_idx
! Subprogram not used                 write (iulog,*) "SNICAR STATS: albout_lcl(bnd)= ",albout_lcl(bnd_idx), " albsfc_lcl(bnd_idx)= ",albsfc_lcl(bnd_idx)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: landtype= ", sfctype
! Subprogram not used                 write (iulog,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
! Subprogram not used                 write (iulog,*) "SNICAR STATS: coszen= ", coszen(c_idx), " flg_slr= ", flg_slr_in
! Subprogram not used 
! Subprogram not used                 write (iulog,*) "SNICAR STATS: soot(-4)= ", mss_cnc_aer_lcl(-4,1)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: soot(-3)= ", mss_cnc_aer_lcl(-3,1)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: soot(-2)= ", mss_cnc_aer_lcl(-2,1)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: soot(-1)= ", mss_cnc_aer_lcl(-1,1)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: soot(0)= ", mss_cnc_aer_lcl(0,1)
! Subprogram not used 
! Subprogram not used                 write (iulog,*) "SNICAR STATS: L_snw(-4)= ", L_snw(-4)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: L_snw(-3)= ", L_snw(-3)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)
! Subprogram not used 
! Subprogram not used                 write (iulog,*) "SNICAR STATS: snw_rds(-4)= ", snw_rds(c_idx,-4)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: snw_rds(-3)= ", snw_rds(c_idx,-3)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: snw_rds(-2)= ", snw_rds(c_idx,-2)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: snw_rds(-1)= ", snw_rds(c_idx,-1)
! Subprogram not used                 write (iulog,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)
! Subprogram not used                 
! Subprogram not used                 call endrun()
! Subprogram not used              endif
! Subprogram not used              
! Subprogram not used           enddo   ! loop over wvl bands
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           ! Weight output NIR albedo appropriately
! Subprogram not used           albout(c_idx,1) = albout_lcl(1)
! Subprogram not used           flx_sum         = 0._r8
! Subprogram not used           do bnd_idx= nir_bnd_bgn,nir_bnd_end
! Subprogram not used              flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
! Subprogram not used           end do
! Subprogram not used           albout(c_idx,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
! Subprogram not used 
! Subprogram not used           ! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
! Subprogram not used           flx_abs(c_idx,:,1) = flx_abs_lcl(:,1)
! Subprogram not used           do i=snl_top,1,1
! Subprogram not used              flx_sum = 0._r8
! Subprogram not used              do bnd_idx= nir_bnd_bgn,nir_bnd_end
! Subprogram not used                 flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
! Subprogram not used              enddo
! Subprogram not used              flx_abs(c_idx,i,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))          
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
! Subprogram not used        elseif ( (coszen(c_idx) > 0._r8) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0._r8) ) then
! Subprogram not used           albout(c_idx,1) = albsfc(c_idx,1)
! Subprogram not used           albout(c_idx,2) = albsfc(c_idx,2)
! Subprogram not used 
! Subprogram not used        ! There is either zero snow, or no sun
! Subprogram not used        else
! Subprogram not used           albout(c_idx,1) = 0._r8
! Subprogram not used           albout(c_idx,2) = 0._r8
! Subprogram not used        endif    ! if column has snow and coszen > 0
! Subprogram not used 
! Subprogram not used     enddo    ! loop over all columns
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine SNICAR_RT


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowAge_grain
!
! !INTERFACE:
  subroutine SnowAge_grain(lbc, ubc, num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)
    !
    ! !DESCRIPTION:
    ! Updates the snow effective grain size (radius). 
    ! Contributions to grain size evolution are from:
    !   1. vapor redistribution (dry snow) 
    !   2. liquid water redistribution (wet snow)
    !   3. re-freezing of liquid water
    ! 
    ! Vapor redistribution: Method is to retrieve 3 best-bit parameters that
    ! depend on snow temperature, temperature gradient, and density,
    ! that are derived from the microphysical model described in: 
    ! Flanner and Zender (2006), Linking snowpack microphysics and albedo
    ! evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834. 
    ! The parametric equation has the form: 
    ! dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where:
    !   r is the effective radius,
    !   tau and kappa are best-fit parameters,
    !   drdt_0 is the initial rate of change of effective radius, and
    !   dr_fresh is the difference between the current and fresh snow states 
    !  (r_current - r_fresh).
    !
    ! Liquid water redistribution: Apply the grain growth function from:
    !   Brun, E. (1989), Investigation of wet-snow metamorphism in respect of 
    !   liquid-water content, Annals of Glaciology, 13, 22-26.
    !   There are two parameters that describe the grain growth rate as 
    !   a function of snow liquid water content (LWC). The "LWC=0" parameter
    !   is zeroed here because we are accounting for dry snowing with a 
    !   different representation
    !
    ! Re-freezing of liquid water: Assume that re-frozen liquid water clumps
    !   into an arbitrarily large effective grain size (snw_rds_refrz). 
    !   The phenomenon is observed (Grenfell), but so far unquantified, as far as 
    !   I am aware.
    !
    !
    ! !USES:
    use clmtype
    use clm_time_manager , only : get_step_size, get_nstep
    use clm_varpar       , only : nlevsno
    use clm_varcon       , only : spval
    use abortutils       , only : endrun
    use shr_const_mod    , only : SHR_CONST_RHOICE, SHR_CONST_PI
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                  ! column bounds
    integer, intent(in) :: num_snowc                 ! number of column snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1)   ! column filter for snow points
    integer, intent(in) :: num_nosnowc               ! number of column non-snow points in column filter
    integer, intent(in) :: filter_nosnowc(ubc-lbc+1) ! column filter for non-snow points
    !
    !
    ! !CALLED FROM: clm_driver1
    !

    ! !LOCAL VARIABLES:
    !
    ! local pointers to implicit arguments
    !

    real(r8), pointer :: t_soisno(:,:)         ! soil and snow temperature (col,lyr) [K]
    integer,  pointer :: snl(:)                ! negative number of snow layers (col) [nbr]
    real(r8), pointer :: t_grnd(:)             ! ground temperature (col) [K]
    real(r8), pointer :: dz(:,:)               ! layer thickness (col,lyr) [m]
    real(r8), pointer :: h2osno(:)             ! snow water (col) [mm H2O]
    real(r8), pointer :: snw_rds(:,:)          ! effective grain radius (col,lyr) [microns, m-6]
    real(r8), pointer :: snw_rds_top(:)        ! effective grain radius, top layer (col) [microns, m-6]
    real(r8), pointer :: sno_liq_top(:)        ! liquid water fraction (mass) in top snow layer (col) [frc]
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water content (col,lyr) [kg m-2]
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice content (col,lyr) [kg m-2]
    real(r8), pointer :: snot_top(:)           ! snow temperature in top layer (col) [K]
    real(r8), pointer :: dTdz_top(:)           ! temperature gradient in top layer (col) [K m-1]
    real(r8), pointer :: qflx_snow_grnd_col(:) ! snow on ground after interception (col) [kg m-2 s-1]
    real(r8), pointer :: qflx_snwcp_ice(:)     ! excess precipitation due to snow capping [kg m-2 s-1]
    real(r8), pointer :: qflx_snofrz_lyr(:,:)  ! snow freezing rate (col,lyr) [kg m-2 s-1]
    logical , pointer :: do_capsnow(:)         ! true => do snow capping
 
    !
    ! !OTHER LOCAL VARIABLES:
    !
    integer :: snl_top                      ! top snow layer index [idx]
    integer :: snl_btm                      ! bottom snow layer index [idx]
    integer :: i                            ! layer index [idx]
    integer :: c_idx                        ! column index [idx]
    integer :: fc                           ! snow column filter index [idx]
    integer :: T_idx                        ! snow aging lookup table temperature index [idx]
    integer :: Tgrd_idx                     ! snow aging lookup table temperature gradient index [idx]
    integer :: rhos_idx                     ! snow aging lookup table snow density index [idx]
    real(r8) :: t_snotop                    ! temperature at upper layer boundary [K]
    real(r8) :: t_snobtm                    ! temperature at lower layer boundary [K]
    real(r8) :: dTdz(lbc:ubc,-nlevsno:0)    ! snow temperature gradient (col,lyr) [K m-1]
    real(r8) :: bst_tau                     ! snow aging parameter retrieved from lookup table [hour]
    real(r8) :: bst_kappa                   ! snow aging parameter retrieved from lookup table [unitless]
    real(r8) :: bst_drdt0                   ! snow aging parameter retrieved from lookup table [um hr-1]
    real(r8) :: dr                          ! incremental change in snow effective radius [um]
    real(r8) :: dr_wet                      ! incremental change in snow effective radius from wet growth [um]
    real(r8) :: dr_fresh                    ! difference between fresh snow r_e and current r_e [um]
    real(r8) :: newsnow                     ! fresh snowfall [kg m-2]
    real(r8) :: refrzsnow                   ! re-frozen snow [kg m-2]
    real(r8) :: frc_newsnow                 ! fraction of layer mass that is new snow [frc]
    real(r8) :: frc_oldsnow                 ! fraction of layer mass that is old snow [frc]
    real(r8) :: frc_refrz                   ! fraction of layer mass that is re-frozen snow [frc]
    real(r8) :: frc_liq                     ! fraction of layer mass that is liquid water[frc]    
    real(r8) :: dtime                       ! land model time step [sec]
    real(r8) :: rhos                        ! snow density [kg m-3]
    real(r8) :: h2osno_lyr                  ! liquid + solid H2O in snow layer [kg m-2]


    ! Assign local pointers to derived subtypes components (column-level)
    t_soisno           => ces%t_soisno
    snl                => cps%snl
    t_grnd             => ces%t_grnd
    dz                 => cps%dz
    h2osno             => cws%h2osno
    snw_rds            => cps%snw_rds
    h2osoi_liq         => cws%h2osoi_liq
    h2osoi_ice         => cws%h2osoi_ice
    snot_top           => cps%snot_top
    dTdz_top           => cps%dTdz_top
    snw_rds_top        => cps%snw_rds_top
    sno_liq_top        => cps%sno_liq_top
    qflx_snow_grnd_col => pwf_a%qflx_snow_grnd
    qflx_snwcp_ice     => pwf_a%qflx_snwcp_ice
    qflx_snofrz_lyr    => cwf%qflx_snofrz_lyr
    do_capsnow         => cps%do_capsnow
  

    ! set timestep and step interval
    dtime = get_step_size()

    ! loop over columns that have at least one snow layer
    do fc = 1, num_snowc
       c_idx = filter_snowc(fc)

       snl_btm = 0
       snl_top = snl(c_idx) + 1

       ! loop over snow layers
       do i=snl_top,snl_btm,1
          !
          !**********  1. DRY SNOW AGING  ***********
          !
          h2osno_lyr = h2osoi_liq(c_idx,i) + h2osoi_ice(c_idx,i)

          ! temperature gradient
          if (i == snl_top) then 
             ! top layer
             t_snotop = t_grnd(c_idx)
             t_snobtm = (t_soisno(c_idx,i+1)*dz(c_idx,i) + t_soisno(c_idx,i)*dz(c_idx,i+1)) / (dz(c_idx,i)+dz(c_idx,i+1))
          else
             t_snotop = (t_soisno(c_idx,i-1)*dz(c_idx,i) + t_soisno(c_idx,i)*dz(c_idx,i-1)) / (dz(c_idx,i)+dz(c_idx,i-1))
             t_snobtm = (t_soisno(c_idx,i+1)*dz(c_idx,i) + t_soisno(c_idx,i)*dz(c_idx,i+1)) / (dz(c_idx,i)+dz(c_idx,i+1))
          endif
          
          dTdz(c_idx,i) = abs((t_snotop - t_snobtm) / dz(c_idx,i))
          
          ! snow density
          rhos = (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i)) / dz(c_idx,i)

          ! best-fit table indecies
          T_idx    = nint((t_soisno(c_idx,i)-223) / 5) + 1
          Tgrd_idx = nint(dTdz(c_idx,i) / 10) + 1
          rhos_idx = nint((rhos-50) / 50) + 1

          ! boundary check:
          if (T_idx < idx_T_min) then 
             T_idx = idx_T_min
          endif
          if (T_idx > idx_T_max) then 
             T_idx = idx_T_max
          endif
          if (Tgrd_idx < idx_Tgrd_min) then 
             Tgrd_idx = idx_Tgrd_min
          endif
          if (Tgrd_idx > idx_Tgrd_max) then 
             Tgrd_idx = idx_Tgrd_max
          endif
          if (rhos_idx < idx_rhos_min) then 
             rhos_idx = idx_rhos_min
          endif
          if (rhos_idx > idx_rhos_max) then 
             rhos_idx = idx_rhos_max
          endif
             
          ! best-fit parameters
          bst_tau   = snowage_tau(rhos_idx,Tgrd_idx,T_idx)
          bst_kappa = snowage_kappa(rhos_idx,Tgrd_idx,T_idx)     
          bst_drdt0 = snowage_drdt0(rhos_idx,Tgrd_idx,T_idx)


          ! change in snow effective radius, using best-fit parameters
          dr_fresh = snw_rds(c_idx,i)-snw_rds_min
          dr = (bst_drdt0*(bst_tau/(dr_fresh+bst_tau))**(1/bst_kappa)) * (dtime/3600)
          

          !
          !**********  2. WET SNOW AGING  ***********
          !
          ! We are assuming wet and dry evolution occur simultaneously, and 
          ! the contributions from both can be summed. 
          ! This is justified by setting the linear offset constant C1_liq_Brun89 to zero [Brun, 1989]
          
          ! liquid water faction
          frc_liq = min(0.1_r8, (h2osoi_liq(c_idx,i) / (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i))))

          !dr_wet = 1E6_r8*(dtime*(C1_liq_Brun89 + C2_liq_Brun89*(frc_liq**(3))) / (4*SHR_CONST_PI*(snw_rds(c_idx,i)/1E6)**(2)))
          !simplified, units of microns:
          dr_wet = 1E18_r8*(dtime*(C2_liq_Brun89*(frc_liq**(3))) / (4*SHR_CONST_PI*snw_rds(c_idx,i)**(2)))

          dr = dr + dr_wet
    
          !
          !**********  3. SNOWAGE SCALING (TURNED OFF BY DEFAULT)  *************
          !
          ! Multiply rate of change of effective radius by some constant, xdrdt
          if (flg_snoage_scl) then
             dr = dr*xdrdt
          endif

          
          !
          !**********  4. INCREMENT EFFECTIVE RADIUS, ACCOUNTING FOR:  ***********
          !               DRY AGING
          !               WET AGING
          !               FRESH SNOW
          !               RE-FREEZING
          !
          ! new snowfall [kg/m2]
          if (do_capsnow(c_idx)) then
             newsnow = max(0._r8, (qflx_snwcp_ice(c_idx)*dtime))
          else
             newsnow = max(0._r8, (qflx_snow_grnd_col(c_idx)*dtime))
          endif

          ! snow that has re-frozen [kg/m2]
          refrzsnow = max(0._r8, (qflx_snofrz_lyr(c_idx,i)*dtime))
          
          ! fraction of layer mass that is re-frozen
          frc_refrz = refrzsnow / h2osno_lyr
                  
          ! fraction of layer mass that is new snow
          if (i == snl_top) then
             frc_newsnow = newsnow / h2osno_lyr
          else
             frc_newsnow = 0._r8
          endif

          if ((frc_refrz + frc_newsnow) > 1._r8) then
             frc_refrz = frc_refrz / (frc_refrz + frc_newsnow)
             frc_newsnow = 1._r8 - frc_refrz
             frc_oldsnow = 0._r8
          else
             frc_oldsnow = 1._r8 - frc_refrz - frc_newsnow
          endif

          ! mass-weighted mean of fresh snow, old snow, and re-frozen snow effective radius
          snw_rds(c_idx,i) = (snw_rds(c_idx,i)+dr)*frc_oldsnow + snw_rds_min*frc_newsnow + snw_rds_refrz*frc_refrz

             
          !
          !**********  5. CHECK BOUNDARIES   ***********
          !
          ! boundary check
          if (snw_rds(c_idx,i) < snw_rds_min) then
             snw_rds(c_idx,i) = snw_rds_min
          endif

          if (snw_rds(c_idx,i) > snw_rds_max) then
             snw_rds(c_idx,i) = snw_rds_max
          end if

          ! set top layer variables for history files
          if (i == snl_top) then
             snot_top(c_idx)    = t_soisno(c_idx,i)
             dTdz_top(c_idx)    = dTdz(c_idx,i)
             snw_rds_top(c_idx) = snw_rds(c_idx,i)
             sno_liq_top(c_idx) = h2osoi_liq(c_idx,i) / (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i))
          endif

       enddo
    enddo

    ! Special case: snow on ground, but not enough to have defined a snow layer:
    !   set snw_rds to fresh snow grain size:
    do fc = 1, num_nosnowc
       c_idx = filter_nosnowc(fc)
       if (h2osno(c_idx) > 0._r8) then
          snw_rds(c_idx,0) = snw_rds_min
       endif
    enddo
        
  end subroutine SnowAge_grain

  subroutine SnowOptics_init( )
    use fileutils       , only : getfil
    use CLM_varctl      , only : fsnowoptics
    use spmdMod         , only : masterproc
    use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile
   
    type(file_desc_t)  :: ncid                        ! netCDF file id
    character(len=256) :: locfn                       ! local filename
    character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
    integer            :: ier                         ! error status


    !
    ! Open optics file:
    if(masterproc) write(iulog,*) 'Attempting to read snow optical properties .....'
    call getfil (fsnowoptics, locfn, 0)
    call ncd_pio_openfile(ncid, locfn, 0)
    if(masterproc) write(iulog,*) subname,trim(fsnowoptics)

    ! direct-beam snow Mie parameters:
    call ncd_io('ss_alb_ice_drc', ss_alb_snw_drc,            'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_ice_drc',asm_prm_snw_drc,          'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc, 'read', ncid, posNOTonfile=.true.)

    ! diffuse snow Mie parameters
    call ncd_io( 'ss_alb_ice_dfs', ss_alb_snw_dfs,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_ice_dfs', asm_prm_snw_dfs,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs, 'read', ncid, posNOTonfile=.true.)

    ! BC species 1 Mie parameters
    call ncd_io( 'ss_alb_bcphil', ss_alb_bc1,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_bcphil', asm_prm_bc1,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_bcphil', ext_cff_mss_bc1, 'read', ncid, posNOTonfile=.true.)

    ! BC species 2 Mie parameters
    call ncd_io( 'ss_alb_bcphob', ss_alb_bc2,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_bcphob', asm_prm_bc2,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_bcphob', ext_cff_mss_bc2, 'read', ncid, posNOTonfile=.true.)

    ! OC species 1 Mie parameters
    call ncd_io( 'ss_alb_ocphil', ss_alb_oc1,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_ocphil', asm_prm_oc1,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_ocphil', ext_cff_mss_oc1, 'read', ncid, posNOTonfile=.true.)

    ! OC species 2 Mie parameters
    call ncd_io( 'ss_alb_ocphob', ss_alb_oc2,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_ocphob', asm_prm_oc2,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_ocphob', ext_cff_mss_oc2, 'read', ncid, posNOTonfile=.true.)

    ! dust species 1 Mie parameters
    call ncd_io( 'ss_alb_dust01', ss_alb_dst1,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_dust01', asm_prm_dst1,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_dust01', ext_cff_mss_dst1, 'read', ncid, posNOTonfile=.true.)

    ! dust species 2 Mie parameters
    call ncd_io( 'ss_alb_dust02', ss_alb_dst2,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_dust02', asm_prm_dst2,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_dust02', ext_cff_mss_dst2, 'read', ncid, posNOTonfile=.true.)

    ! dust species 3 Mie parameters
    call ncd_io( 'ss_alb_dust03', ss_alb_dst3,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_dust03', asm_prm_dst3,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_dust03', ext_cff_mss_dst3, 'read', ncid, posNOTonfile=.true.)

    ! dust species 4 Mie parameters
    call ncd_io( 'ss_alb_dust04', ss_alb_dst4,           'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'asm_prm_dust04', asm_prm_dst4,         'read', ncid, posNOTonfile=.true.)
    call ncd_io( 'ext_cff_mss_dust04', ext_cff_mss_dst4, 'read', ncid, posNOTonfile=.true.)


    call ncd_pio_closefile(ncid)
    if (masterproc) then

       write(iulog,*) 'Successfully read snow optical properties'
       ! print some diagnostics:
       write (iulog,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
            ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
            ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
       write (iulog,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
            ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
            ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
       if (DO_SNO_OC) then
          write (iulog,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
       else
          write (iulog,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
       endif
       write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
            ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
       write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
            ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
       if (DO_SNO_OC) then
          write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
               ss_alb_oc1(1), ss_alb_oc1(2), ss_alb_oc1(3), ss_alb_oc1(4), ss_alb_oc1(5)
          write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
               ss_alb_oc2(1), ss_alb_oc2(2), ss_alb_oc2(3), ss_alb_oc2(4), ss_alb_oc2(5)
       endif
       write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
            ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
       write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
            ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
       write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
            ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
       write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
            ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
       write(iulog,*)
    end if

  end subroutine SnowOptics_init

  subroutine SnowAge_init( )
   use CLM_varctl      , only : fsnowaging
   use fileutils       , only : getfil
   use spmdMod         , only : masterproc
   use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

   type(file_desc_t)  :: ncid                        ! netCDF file id
   character(len=256) :: locfn                       ! local filename
   character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
   integer            :: varid                       ! netCDF id's
   integer            :: ier                         ! error status

   ! Open snow aging (effective radius evolution) file:
   allocate(snowage_tau(idx_rhos_max,idx_Tgrd_max,idx_T_max))
   allocate(snowage_kappa(idx_rhos_max,idx_Tgrd_max,idx_T_max))
   allocate(snowage_drdt0(idx_rhos_max,idx_Tgrd_max,idx_T_max))

   if(masterproc)  write(iulog,*) 'Attempting to read snow aging parameters .....'
   call getfil (fsnowaging, locfn, 0)
   call ncd_pio_openfile(ncid, locfn, 0)
   if(masterproc) write(iulog,*) subname,trim(fsnowaging)
   
    ! snow aging parameters
   
   call ncd_io('tau', snowage_tau,       'read', ncid, posNOTonfile=.true.)
   call ncd_io('kappa', snowage_kappa,   'read', ncid, posNOTonfile=.true.)
   call ncd_io('drdsdt0', snowage_drdt0, 'read', ncid, posNOTonfile=.true.)

   call ncd_pio_closefile(ncid)
   if (masterproc) then
      
      write(iulog,*) 'Successfully read snow aging properties'
      
      ! print some diagnostics:
      write (iulog,*) 'SNICAR: snowage tau for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_tau(3,11,9)
      write (iulog,*) 'SNICAR: snowage kappa for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_kappa(3,11,9)
      write (iulog,*) 'SNICAR: snowage dr/dt_0 for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_drdt0(3,11,9)
   endif

  end subroutine SnowAge_init


end module SNICARMod
