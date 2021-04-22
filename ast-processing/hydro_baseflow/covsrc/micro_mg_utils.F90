module micro_mg_utils

!--------------------------------------------------------------------------
!
! This module contains process rates and utility functions used by the MG
! microphysics.
!
! Original MG authors: Andrew Gettelman, Hugh Morrison
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
!
! Separated from MG 1.5 by B. Eaton.
! Separated module switched to MG 2.0 and further changes by S. Santos.
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
!--------------------------------------------------------------------------
!
! List of required external functions that must be supplied:
!   gamma --> standard mathematical gamma function (if gamma is an
!       intrinsic, define HAVE_GAMMA_INTRINSICS)
!
!--------------------------------------------------------------------------
!
! Constants that must be specified in the "init" method (module variables):
!
! kind            kind of reals (to verify correct linkage only) -
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                   J kg-1 K-1
! rh2o            gas constant for water vapor                   J kg-1 K-1
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! tmelt           temperature of melting point for water         K
! latvap          latent heat of vaporization                    J kg-1
! latice          latent heat of fusion                          J kg-1
!
!--------------------------------------------------------------------------


use shr_spfn_mod, only: gamma => shr_spfn_gamma


implicit none
private
save

public :: &
     micro_mg_utils_init, &
     size_dist_param_liq, &
     size_dist_param_basic, &
     avg_diameter, &
     ice_deposition_sublimation, &
     kk2000_liq_autoconversion, &
     ice_autoconversion, &
     immersion_freezing, &
     contact_freezing, &
     snow_self_aggregation, &
     accrete_cloud_water_snow, &
     secondary_ice_production, &
     accrete_rain_snow, &
     heterogeneous_rain_freezing, &
     accrete_cloud_water_rain, &
     self_collection_rain, &
     accrete_cloud_ice_snow, &
     evaporate_sublimate_precip, &
     bergeron_process_snow

! 8 byte real and integer
integer, parameter, public :: r8 = selected_real_kind(12)
integer, parameter, public :: i8 = selected_int_kind(18)

public :: MGHydrometeorProps

type :: MGHydrometeorProps
   ! Density (kg/m^3)
   real(r8) :: rho
   ! Information for size calculations.
   ! Basic calculation of mean size is:
   !     lambda = (shape_coef*nic/qic)^(1/eff_dim)
   ! Then lambda is constrained by bounds.
   real(r8) :: eff_dim
   real(r8) :: shape_coef
   real(r8) :: lambda_bounds(2)
   ! Minimum average particle mass (kg).
   ! Limit is applied at the beginning of the size distribution calculations.
   real(r8) :: min_mean_mass
end type MGHydrometeorProps

interface MGHydrometeorProps
   module procedure NewMGHydrometeorProps
end interface

type(MGHydrometeorProps), public :: mg_liq_props
type(MGHydrometeorProps), public :: mg_ice_props
type(MGHydrometeorProps), public :: mg_rain_props
type(MGHydrometeorProps), public :: mg_snow_props

!=================================================
! Public module parameters (mostly for MG itself)
!=================================================

! Pi to 20 digits; more than enough to reach the limit of double precision.
real(r8), parameter, public :: pi = 3.14159265358979323846_r8

! "One minus small number": number near unity for round-off issues.
real(r8), parameter, public :: omsm   = 1._r8 - 1.e-5_r8

! Smallest mixing ratio considered in microphysics.
real(r8), parameter, public :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter, public :: mincld = 0.0001_r8

real(r8), parameter, public :: rhosn = 250._r8  ! bulk density snow
real(r8), parameter, public :: rhoi = 500._r8   ! bulk density ice
real(r8), parameter, public :: rhow = 1000._r8  ! bulk density liquid
real(r8), parameter, public :: rhows = 917._r8  ! bulk density water solid

! autoconversion size threshold for cloud ice to snow (m)
real(r8), parameter, public :: dcs = 90.e-6_r8

! fall speed parameters, V = aD^b (V is in m/s)
! droplets
real(r8), parameter, public :: ac = 3.e7_r8
real(r8), parameter, public :: bc = 2._r8
! snow
real(r8), parameter, public :: as = 11.72_r8
real(r8), parameter, public :: bs = 0.41_r8
! cloud ice
real(r8), parameter, public :: ai = 700._r8
real(r8), parameter, public :: bi = 1._r8
! rain
real(r8), parameter, public :: ar = 841.99667_r8
real(r8), parameter, public :: br = 0.8_r8

! mass of new crystal due to aerosol freezing and growth (kg)
real(r8), parameter, public :: mi0 = &
     4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

!=================================================
! Private module parameters
!=================================================

! Signaling NaN bit pattern that represents a limiter that's turned off.
integer(i8), parameter :: limiter_off = int(Z'7FF1111111111111', i8)

! alternate threshold used for some in-cloud mmr
real(r8), parameter :: icsmall = 1.e-8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(r8), parameter :: dsph = 3._r8

! Bounds for mean diameter for different constituents.
! (E.g. ice must be at least 10 microns but no more than twice the
! threshold for autoconversion to snow.
real(r8), parameter :: lam_bnd_ice(2) = 1._r8/[2._r8*dcs, 10.e-6_r8]
real(r8), parameter :: lam_bnd_rain(2) = 1._r8/[500.e-6_r8, 20.e-6_r8]
real(r8), parameter :: lam_bnd_snow(2) = 1._r8/[2000.e-6_r8, 10.e-6_r8]

! Minimum average mass of particles.
real(r8), parameter :: min_mean_mass_liq = 1.e-20_r8
real(r8), parameter :: min_mean_mass_ice = 1.e-20_r8

! ventilation parameters
! for snow
real(r8), parameter :: f1s = 0.86_r8
real(r8), parameter :: f2s = 0.28_r8
! for rain
real(r8), parameter :: f1r = 0.78_r8
real(r8), parameter :: f2r = 0.308_r8

! collection efficiencies
! aggregation of cloud ice and snow
real(r8), parameter :: eii = 0.1_r8

! immersion freezing parameters, bigg 1953
real(r8), parameter :: bimm = 100._r8
real(r8), parameter :: aimm = 0.66_r8

! Mass of each raindrop created from autoconversion.
real(r8), parameter :: droplet_mass_25um = 4._r8/3._r8*pi*rhow*(25.e-6_r8)**3

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

! additional constants to help speed up code
real(r8) :: gamma_bs_plus3
real(r8) :: gamma_half_br_plus5
real(r8) :: gamma_half_bs_plus5

!==========================================================================
contains
!==========================================================================

! Initialize module variables.
!
! "kind" serves no purpose here except to check for unlikely linking
! issues; always pass in the kind for a double precision real.
!
! "errstring" is the only output; it is blank if there is no error, or set
! to a message if there is an error.
!
! Check the list at the top of this module for descriptions of all other
! arguments.
subroutine micro_mg_utils_init( kind, rh2o, cpair, tmelt_in, latvap, &
     latice, errstring)

  integer,  intent(in)  :: kind
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice

  character(128), intent(out) :: errstring

  !-----------------------------------------------------------------------

  errstring = ' '

  if( kind .ne. r8 ) then
     errstring = 'micro_mg_init: KIND of reals does not match'
     return
  endif

  ! declarations for MG code (transforms variable names)

  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_bs_plus3=gamma(3._r8+bs)
  gamma_half_br_plus5=gamma(5._r8/2._r8+br/2._r8)
  gamma_half_bs_plus5=gamma(5._r8/2._r8+bs/2._r8)

  ! Don't specify lambda bounds for cloud liquid, as they are determined by
  ! pgam dynamically.
  mg_liq_props = MGHydrometeorProps(rhow, dsph, min_mean_mass=min_mean_mass_liq)
  mg_ice_props = MGHydrometeorProps(rhoi, dsph, lam_bnd_ice, min_mean_mass_ice)
  mg_rain_props = MGHydrometeorProps(rhow, dsph, lam_bnd_rain)
  mg_snow_props = MGHydrometeorProps(rhosn, dsph, lam_bnd_snow)

end subroutine micro_mg_utils_init

! Constructor for a constituent property object.
function NewMGHydrometeorProps(rho, eff_dim, lambda_bounds, min_mean_mass) &
     result(res)
  real(r8), intent(in) :: rho, eff_dim
  real(r8), intent(in), optional :: lambda_bounds(2), min_mean_mass
  type(MGHydrometeorProps) :: res

  res%rho = rho
  res%eff_dim = eff_dim
  if (present(lambda_bounds)) then
     res%lambda_bounds = lambda_bounds
  else
     res%lambda_bounds = no_limiter()
  end if
  if (present(min_mean_mass)) then
     res%min_mean_mass = min_mean_mass
  else
     res%min_mean_mass = no_limiter()
  end if

  res%shape_coef = rho*pi*gamma(eff_dim+1._r8)/6._r8

end function NewMGHydrometeorProps

!========================================================================
!FORMULAS
!========================================================================

! Calculate correction due to latent heat for evaporation/sublimation
! Subprogram not used elemental function calc_ab(t, qv, xxl) result(ab)
! Subprogram not used   real(r8), intent(in) :: t     ! Temperature
! Subprogram not used   real(r8), intent(in) :: qv    ! Saturation vapor pressure
! Subprogram not used   real(r8), intent(in) :: xxl   ! Latent heat
! Subprogram not used 
! Subprogram not used   real(r8) :: ab
! Subprogram not used 
! Subprogram not used   real(r8) :: dqsdt
! Subprogram not used 
! Subprogram not used   dqsdt = xxl*qv / (rv * t**2)
! Subprogram not used   ab = 1._r8 + dqsdt*xxl/cpp
! Subprogram not used 
! Subprogram not used end function calc_ab

! get cloud droplet size distribution parameters
elemental subroutine size_dist_param_liq(props, qcic, ncic, rho, pgam, lamc)
  type(MGHydrometeorProps), intent(in) :: props
  real(r8), intent(in) :: qcic
  real(r8), intent(inout) :: ncic
  real(r8), intent(in) :: rho

  real(r8), intent(out) :: pgam
  real(r8), intent(out) :: lamc

  type(MGHydrometeorProps) :: props_loc

  if (qcic > qsmall) then

     ! Local copy of properties that can be modified.
     ! (Elemental routines that operate on arrays can't modify scalar
     ! arguments.)
     props_loc = props

     ! Get pgam from fit to observations of martin et al. 1994
     pgam = 0.0005714_r8*(ncic/1.e6_r8*rho) + 0.2714_r8
     pgam = 1._r8/(pgam**2) - 1._r8
     pgam = max(pgam, 2._r8)
     pgam = min(pgam, 15._r8)

     ! Set coefficient for use in size_dist_param_basic.
     props_loc%shape_coef = pi * props_loc%rho / 6._r8 * &
          rising_factorial(pgam+1._r8, props_loc%eff_dim)

     ! Limit to between 2 and 50 microns mean size.
     props_loc%lambda_bounds = (pgam+1._r8)*1._r8/[50.e-6_r8, 2.e-6_r8]

     call size_dist_param_basic(props_loc, qcic, ncic, lamc)

  else
     ! pgam not calculated in this case, so set it to a value likely to
     ! cause an error if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._r8
     lamc = 0._r8
  end if

contains

  ! Use gamma function to implement rising factorial extended to the reals.
  pure function rising_factorial(x, n)
    real(r8), intent(in) :: x, n
    real(r8) :: rising_factorial

    rising_factorial = gamma(x+n)/gamma(x)

  end function rising_factorial

end subroutine size_dist_param_liq

! Basic routine for getting size distribution parameters.
elemental subroutine size_dist_param_basic(props, qic, nic, lam, n0)
  type(MGHydrometeorProps), intent(in) :: props
  real(r8), intent(in) :: qic
  real(r8), intent(inout) :: nic

  real(r8), intent(out) :: lam
  real(r8), intent(out), optional :: n0

  if (qic > qsmall) then

     ! add upper limit to in-cloud number concentration to prevent
     ! numerical error
     if (limiter_is_on(props%min_mean_mass)) then
        nic = min(nic, qic / props%min_mean_mass)
     end if

     ! lambda = (c n/q)^(1/d)
     lam = (props%shape_coef * nic/qic)**(1._r8/props%eff_dim)

     ! check for slope
     ! adjust vars
     if (lam < props%lambda_bounds(1)) then
        lam = props%lambda_bounds(1)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     else if (lam > props%lambda_bounds(2)) then
        lam = props%lambda_bounds(2)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     end if

  else
     lam = 0._r8
  end if

  if (present(n0)) n0 = nic * lam

end subroutine size_dist_param_basic

real(r8) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  real(r8), intent(in) :: q         ! mass mixing ratio
  real(r8), intent(in) :: n         ! number concentration
  real(r8), intent(in) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi * rho_sub * n/(q*rho_air))**(-1._r8/3._r8)

end function avg_diameter

! Subprogram not used real(r8) elemental function var_coef(relvar, a)
! Subprogram not used   ! Finds a coefficient for process rates based on the relative variance
! Subprogram not used   ! of cloud water.
! Subprogram not used   real(r8), intent(in) :: relvar
! Subprogram not used   real(r8), intent(in) :: a
! Subprogram not used 
! Subprogram not used   var_coef = gamma(relvar + a) / (gamma(relvar) * relvar**a)
! Subprogram not used 
! Subprogram not used end function var_coef

!========================================================================
!MICROPHYSICAL PROCESS CALCULATIONS
!========================================================================

!========================================================================
! Initial ice deposition and sublimation loop.
! Run before the main loop
! This subroutine written by Peter Caldwell

! Subprogram not used elemental subroutine ice_deposition_sublimation(t, qv, qi, ni, &
! Subprogram not used                                                 icldm, rho, dv,qvl, qvi, &
! Subprogram not used                                                 berg, vap_dep, ice_sublim)
! Subprogram not used 
! Subprogram not used   !INPUT VARS:
! Subprogram not used   !===============================================
! Subprogram not used   real(r8), intent(in) :: t
! Subprogram not used   real(r8), intent(in) :: qv
! Subprogram not used   real(r8), intent(in) :: qi
! Subprogram not used   real(r8), intent(in) :: ni
! Subprogram not used   real(r8), intent(in) :: icldm
! Subprogram not used   real(r8), intent(in) :: rho
! Subprogram not used   real(r8), intent(in) :: dv
! Subprogram not used   real(r8), intent(in) :: qvl
! Subprogram not used   real(r8), intent(in) :: qvi
! Subprogram not used 
! Subprogram not used   !OUTPUT VARS:
! Subprogram not used   !===============================================
! Subprogram not used   real(r8), intent(out) :: vap_dep !ice deposition (cell-ave value)
! Subprogram not used   real(r8), intent(out) :: ice_sublim !ice sublimation (cell-ave value)
! Subprogram not used   real(r8), intent(out) :: berg !bergeron enhancement (cell-ave value)
! Subprogram not used 
! Subprogram not used   !INTERNAL VARS:
! Subprogram not used   !===============================================
! Subprogram not used   real(r8) :: ab
! Subprogram not used   real(r8) :: epsi
! Subprogram not used   real(r8) :: qiic
! Subprogram not used   real(r8) :: niic
! Subprogram not used   real(r8) :: lami
! Subprogram not used   real(r8) :: n0i
! Subprogram not used 
! Subprogram not used   if (qi>=qsmall) then
! Subprogram not used 
! Subprogram not used      !GET IN-CLOUD qi, ni
! Subprogram not used      !===============================================
! Subprogram not used      qiic = qi/icldm
! Subprogram not used      niic = ni/icldm
! Subprogram not used 
! Subprogram not used      !Compute linearized condensational heating correction
! Subprogram not used      ab=calc_ab(t, qvi, xxls)
! Subprogram not used      !Get slope and intercept of gamma distn for ice.
! Subprogram not used      call size_dist_param_basic(mg_ice_props, qiic, niic, lami, n0i)
! Subprogram not used      !Get depletion timescale=1/eps
! Subprogram not used      epsi = 2._r8*pi*n0i*rho*Dv/(lami*lami)
! Subprogram not used 
! Subprogram not used      !Compute deposition/sublimation
! Subprogram not used      vap_dep = epsi/ab*(qv - qvi)
! Subprogram not used 
! Subprogram not used      !Make this a grid-averaged quantity
! Subprogram not used      vap_dep=vap_dep*icldm
! Subprogram not used 
! Subprogram not used      !Split into deposition or sublimation.
! Subprogram not used      if (t < tmelt .and. vap_dep>0._r8) then
! Subprogram not used         ice_sublim=0._r8
! Subprogram not used      else
! Subprogram not used      !hm, make ice_sublim negative for consistency with other evap/sub processes
! Subprogram not used         ice_sublim=min(vap_dep,0._r8)
! Subprogram not used         vap_dep=0._r8
! Subprogram not used      end if
! Subprogram not used 
! Subprogram not used      !sublimation occurs @ any T. Not so for berg.
! Subprogram not used      if (t < tmelt) then
! Subprogram not used 
! Subprogram not used         !Compute bergeron rate assuming cloud for whole step.
! Subprogram not used         berg = max(epsi/ab*(qvl - qvi), 0._r8)
! Subprogram not used      else !T>frz
! Subprogram not used         berg=0._r8
! Subprogram not used      end if !T<frz
! Subprogram not used 
! Subprogram not used   else !where qi<qsmall
! Subprogram not used      berg=0._r8
! Subprogram not used      vap_dep=0._r8
! Subprogram not used      ice_sublim=0._r8
! Subprogram not used   end if !qi>qsmall
! Subprogram not used 
! Subprogram not used end subroutine ice_deposition_sublimation

!========================================================================
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error

! Subprogram not used elemental subroutine kk2000_liq_autoconversion(microp_uniform, qcic, &
! Subprogram not used      ncic, rho, relvar, prc, nprc, nprc1)
! Subprogram not used 
! Subprogram not used   logical, intent(in) :: microp_uniform
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: qcic
! Subprogram not used   real(r8), intent(in) :: ncic
! Subprogram not used   real(r8), intent(in) :: rho
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: relvar
! Subprogram not used 
! Subprogram not used   real(r8), intent(out) :: prc
! Subprogram not used   real(r8), intent(out) :: nprc
! Subprogram not used   real(r8), intent(out) :: nprc1
! Subprogram not used 
! Subprogram not used   real(r8) :: prc_coef
! Subprogram not used 
! Subprogram not used   ! Take variance into account, or use uniform value.
! Subprogram not used   if (.not. microp_uniform) then
! Subprogram not used      prc_coef = var_coef(relvar, 2.47_r8)
! Subprogram not used   else
! Subprogram not used      prc_coef = 1._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   if (qcic >= icsmall) then
! Subprogram not used 
! Subprogram not used      ! nprc is increase in rain number conc due to autoconversion
! Subprogram not used      ! nprc1 is decrease in cloud droplet conc due to autoconversion
! Subprogram not used 
! Subprogram not used      ! assume exponential sub-grid distribution of qc, resulting in additional
! Subprogram not used      ! factor related to qcvar below
! Subprogram not used      ! hm switch for sub-columns, don't include sub-grid qc
! Subprogram not used 
! Subprogram not used      prc = prc_coef * &
! Subprogram not used           1350._r8 * qcic**2.47_r8 * (ncic/1.e6_r8*rho)**(-1.79_r8)
! Subprogram not used      nprc = prc/droplet_mass_25um
! Subprogram not used      nprc1 = prc/(qcic/ncic)
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used      prc=0._r8
! Subprogram not used      nprc=0._r8
! Subprogram not used      nprc1=0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine kk2000_liq_autoconversion

!========================================================================
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)

! Subprogram not used elemental subroutine ice_autoconversion(t, qiic, lami, n0i, dcs, prci, nprci)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t
! Subprogram not used   real(r8), intent(in) :: qiic
! Subprogram not used   real(r8), intent(in) :: lami
! Subprogram not used   real(r8), intent(in) :: n0i
! Subprogram not used   real(r8), intent(in) :: dcs
! Subprogram not used 
! Subprogram not used   real(r8), intent(out) :: prci
! Subprogram not used   real(r8), intent(out) :: nprci
! Subprogram not used 
! Subprogram not used   ! Assume autoconversion timescale of 180 seconds.
! Subprogram not used   real(r8), parameter :: ac_time = 180._r8
! Subprogram not used 
! Subprogram not used   if (t <= tmelt .and. qiic >= qsmall) then
! Subprogram not used 
! Subprogram not used      nprci = n0i/(lami*ac_time)*exp(-lami*dcs)
! Subprogram not used 
! Subprogram not used      prci = pi*rhoi*n0i/(6._r8*ac_time)* &
! Subprogram not used           (dcs**3/lami+3._r8*dcs**2/lami**2+ &
! Subprogram not used           6._r8*dcs/lami**3+6._r8/lami**4)*exp(-lami*dcs)
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used      prci=0._r8
! Subprogram not used      nprci=0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine ice_autoconversion

! immersion freezing (Bigg, 1953)
!===================================

! Subprogram not used elemental subroutine immersion_freezing(microp_uniform, t, pgam, lamc, &
! Subprogram not used      cdist1, qcic, relvar, mnuccc, nnuccc)
! Subprogram not used 
! Subprogram not used   logical, intent(in) :: microp_uniform
! Subprogram not used 
! Subprogram not used   ! Temperature
! Subprogram not used   real(r8), intent(in) :: t
! Subprogram not used 
! Subprogram not used   ! Cloud droplet size distribution parameters
! Subprogram not used   real(r8), intent(in) :: pgam
! Subprogram not used   real(r8), intent(in) :: lamc
! Subprogram not used   real(r8), intent(in) :: cdist1
! Subprogram not used 
! Subprogram not used   ! MMR of in-cloud liquid water
! Subprogram not used   real(r8), intent(in) :: qcic
! Subprogram not used 
! Subprogram not used   ! Relative variance of cloud water
! Subprogram not used   real(r8), intent(in) :: relvar
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: mnuccc ! MMR
! Subprogram not used   real(r8), intent(out) :: nnuccc ! Number
! Subprogram not used 
! Subprogram not used   ! Coefficients that will be omitted for sub-columns
! Subprogram not used   real(r8) :: dum, dum1
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if (.not. microp_uniform) then
! Subprogram not used      dum = var_coef(relvar, 2._r8)
! Subprogram not used      dum1 = var_coef(relvar, 1._r8)
! Subprogram not used   else
! Subprogram not used      dum = 1._r8
! Subprogram not used      dum1 = 1._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   if (qcic >= qsmall .and. t < 269.15_r8) then
! Subprogram not used 
! Subprogram not used      mnuccc = dum * &
! Subprogram not used           pi*pi/36._r8*rhow* &
! Subprogram not used           cdist1*gamma(7._r8+pgam)* &
! Subprogram not used           bimm*(exp(aimm*(tmelt - t))-1._r8)/lamc**3/lamc**3
! Subprogram not used 
! Subprogram not used      nnuccc = dum1 * &
! Subprogram not used           pi/6._r8*cdist1*gamma(pgam+4._r8) &
! Subprogram not used           *bimm*(exp(aimm*(tmelt - t))-1._r8)/lamc**3
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used      mnuccc = 0._r8
! Subprogram not used      nnuccc = 0._r8
! Subprogram not used   end if ! qcic > qsmall and t < 4 deg C
! Subprogram not used 
! Subprogram not used end subroutine immersion_freezing

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!===================================================================
! dust size and number in multiple bins are read in from companion routine

! Subprogram not used pure subroutine contact_freezing (microp_uniform, t, p, rndst, nacon, &
! Subprogram not used      pgam, lamc, cdist1, qcic, relvar, mnucct, nnucct)
! Subprogram not used 
! Subprogram not used   logical, intent(in) :: microp_uniform
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t(:)            ! Temperature
! Subprogram not used   real(r8), intent(in) :: p(:)            ! Pressure
! Subprogram not used   real(r8), intent(in) :: rndst(:,:)      ! Radius (for multiple dust bins)
! Subprogram not used   real(r8), intent(in) :: nacon(:,:)      ! Number (for multiple dust bins)
! Subprogram not used 
! Subprogram not used   ! Size distribution parameters for cloud droplets
! Subprogram not used   real(r8), intent(in) :: pgam(:)
! Subprogram not used   real(r8), intent(in) :: lamc(:)
! Subprogram not used   real(r8), intent(in) :: cdist1(:)
! Subprogram not used 
! Subprogram not used   ! MMR of in-cloud liquid water
! Subprogram not used   real(r8), intent(in) :: qcic(:)
! Subprogram not used 
! Subprogram not used   ! Relative cloud water variance
! Subprogram not used   real(r8), intent(in) :: relvar(:)
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: mnucct(:) ! MMR
! Subprogram not used   real(r8), intent(out) :: nnucct(:) ! Number
! Subprogram not used 
! Subprogram not used   real(r8) :: tcnt                  ! scaled relative temperature
! Subprogram not used   real(r8) :: viscosity             ! temperature-specific viscosity (kg/m/s)
! Subprogram not used   real(r8) :: mfp                   ! temperature-specific mean free path (m)
! Subprogram not used 
! Subprogram not used   ! Dimension these according to number of dust bins, inferred from rndst size
! Subprogram not used   real(r8) :: nslip(size(rndst,2))  ! slip correction factors
! Subprogram not used   real(r8) :: ndfaer(size(rndst,2)) ! aerosol diffusivities (m^2/sec)
! Subprogram not used 
! Subprogram not used   ! Coefficients not used for subcolumns
! Subprogram not used   real(r8) :: dum, dum1
! Subprogram not used 
! Subprogram not used   integer  :: i
! Subprogram not used 
! Subprogram not used   do i = 1,size(t)
! Subprogram not used 
! Subprogram not used      if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then
! Subprogram not used 
! Subprogram not used         if (.not. microp_uniform) then
! Subprogram not used            dum = var_coef(relvar(i), 4._r8/3._r8)
! Subprogram not used            dum1 = var_coef(relvar(i), 1._r8/3._r8)
! Subprogram not used         else
! Subprogram not used            dum = 1._r8
! Subprogram not used            dum1 = 1._r8
! Subprogram not used         endif
! Subprogram not used 
! Subprogram not used         tcnt=(270.16_r8-t(i))**1.3_r8
! Subprogram not used         viscosity = 1.8e-5_r8*(t(i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
! Subprogram not used         mfp = 2.0_r8*viscosity/ &                         ! Mean free path (m)
! Subprogram not used                      (p(i)*sqrt( 8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i)) ))
! Subprogram not used 
! Subprogram not used         ! Note that these two are vectors.
! Subprogram not used         nslip = 1.0_r8+(mfp/rndst(i,:))*(1.257_r8+(0.4_r8*exp(-(1.1_r8*rndst(i,:)/mfp))))! Slip correction factor
! Subprogram not used 
! Subprogram not used         ndfaer = 1.381e-23_r8*t(i)*nslip/(6._r8*pi*viscosity*rndst(i,:))  ! aerosol diffusivity (m2/s)
! Subprogram not used 
! Subprogram not used         mnucct(i) = dum *  &
! Subprogram not used              dot_product(ndfaer,nacon(i,:)*tcnt)*pi*pi/3._r8*rhow* &
! Subprogram not used              cdist1(i)*gamma(pgam(i)+5._r8)/lamc(i)**4
! Subprogram not used 
! Subprogram not used         nnucct(i) =  dum1 *  &
! Subprogram not used              dot_product(ndfaer,nacon(i,:)*tcnt)*2._r8*pi*  &
! Subprogram not used              cdist1(i)*gamma(pgam(i)+2._r8)/lamc(i)
! Subprogram not used 
! Subprogram not used      else
! Subprogram not used 
! Subprogram not used         mnucct(i)=0._r8
! Subprogram not used         nnucct(i)=0._r8
! Subprogram not used 
! Subprogram not used      end if ! qcic > qsmall and t < 4 deg C
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine contact_freezing

! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!===================================================================
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

! Subprogram not used elemental subroutine snow_self_aggregation(t, rho, asn, rhosn, qsic, nsic, nsagg)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t     ! Temperature
! Subprogram not used   real(r8), intent(in) :: rho   ! Density
! Subprogram not used   real(r8), intent(in) :: asn   ! fall speed parameter for snow
! Subprogram not used   real(r8), intent(in) :: rhosn ! density of snow
! Subprogram not used 
! Subprogram not used   ! In-cloud snow
! Subprogram not used   real(r8), intent(in) :: qsic ! MMR
! Subprogram not used   real(r8), intent(in) :: nsic ! Number
! Subprogram not used 
! Subprogram not used   ! Output number tendency
! Subprogram not used   real(r8), intent(out) :: nsagg
! Subprogram not used 
! Subprogram not used   if (qsic >= qsmall .and. t <= tmelt) then
! Subprogram not used      nsagg = -1108._r8*asn*eii* &
! Subprogram not used           pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)* &
! Subprogram not used           rho**((2._r8+bs)/3._r8)*qsic**((2._r8+bs)/3._r8)* &
! Subprogram not used           (nsic*rho)**((4._r8-bs)/3._r8) /(4._r8*720._r8*rho)
! Subprogram not used   else
! Subprogram not used      nsagg=0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine snow_self_aggregation

! accretion of cloud droplets onto snow/graupel
!===================================================================
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

! Subprogram not used elemental subroutine accrete_cloud_water_snow(t, rho, asn, uns, mu, qcic, ncic, qsic, &
! Subprogram not used      pgam, lamc, lams, n0s, psacws, npsacws)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t   ! Temperature
! Subprogram not used   real(r8), intent(in) :: rho ! Density
! Subprogram not used   real(r8), intent(in) :: asn ! Fallspeed parameter (snow)
! Subprogram not used   real(r8), intent(in) :: uns ! Current fallspeed   (snow)
! Subprogram not used   real(r8), intent(in) :: mu  ! Viscosity
! Subprogram not used 
! Subprogram not used   ! In-cloud liquid water
! Subprogram not used   real(r8), intent(in) :: qcic ! MMR
! Subprogram not used   real(r8), intent(in) :: ncic ! Number
! Subprogram not used 
! Subprogram not used   ! In-cloud snow
! Subprogram not used   real(r8), intent(in) :: qsic ! MMR
! Subprogram not used 
! Subprogram not used   ! Cloud droplet size parameters
! Subprogram not used   real(r8), intent(in) :: pgam
! Subprogram not used   real(r8), intent(in) :: lamc
! Subprogram not used 
! Subprogram not used   ! Snow size parameters
! Subprogram not used   real(r8), intent(in) :: lams
! Subprogram not used   real(r8), intent(in) :: n0s
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: psacws  ! Mass mixing ratio
! Subprogram not used   real(r8), intent(out) :: npsacws ! Number concentration
! Subprogram not used 
! Subprogram not used   real(r8) :: dc0 ! Provisional mean droplet size
! Subprogram not used   real(r8) :: dum
! Subprogram not used   real(r8) :: eci ! collection efficiency for riming of snow by droplets
! Subprogram not used 
! Subprogram not used   ! ignore collision of snow with droplets above freezing
! Subprogram not used 
! Subprogram not used   if (qsic >= qsmall .and. t <= tmelt .and. qcic >= qsmall) then
! Subprogram not used 
! Subprogram not used      ! put in size dependent collection efficiency
! Subprogram not used      ! mean diameter of snow is area-weighted, since
! Subprogram not used      ! accretion is function of crystal geometric area
! Subprogram not used      ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
! Subprogram not used 
! Subprogram not used      dc0 = (pgam+1._r8)/lamc
! Subprogram not used      dum = dc0*dc0*uns*rhow/(9._r8*mu*(1._r8/lams))
! Subprogram not used      eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
! Subprogram not used 
! Subprogram not used      eci = max(eci,0._r8)
! Subprogram not used      eci = min(eci,1._r8)
! Subprogram not used 
! Subprogram not used      ! no impact of sub-grid distribution of qc since psacws
! Subprogram not used      ! is linear in qc
! Subprogram not used 
! Subprogram not used      psacws = pi/4._r8*asn*qcic*rho*n0s*eci*gamma_bs_plus3 / lams**(bs+3._r8)
! Subprogram not used      npsacws = pi/4._r8*asn*ncic*rho*n0s*eci*gamma_bs_plus3 / lams**(bs+3._r8)
! Subprogram not used   else
! Subprogram not used      psacws = 0._r8
! Subprogram not used      npsacws = 0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine accrete_cloud_water_snow

! add secondary ice production due to accretion of droplets by snow
!===================================================================
! (Hallet-Mossop process) (from Cotton et al., 1986)

! Subprogram not used elemental subroutine secondary_ice_production(t, psacws, msacwi, nsacwi)
! Subprogram not used   real(r8), intent(in) :: t ! Temperature
! Subprogram not used 
! Subprogram not used   ! Accretion of cloud water to snow tendencies
! Subprogram not used   real(r8), intent(inout) :: psacws ! MMR
! Subprogram not used 
! Subprogram not used   ! Output (ice) tendencies
! Subprogram not used   real(r8), intent(out) :: msacwi ! MMR
! Subprogram not used   real(r8), intent(out) :: nsacwi ! Number
! Subprogram not used 
! Subprogram not used   if((t < 270.16_r8) .and. (t >= 268.16_r8)) then
! Subprogram not used      nsacwi = 3.5e8_r8*(270.16_r8-t)/2.0_r8*psacws
! Subprogram not used      msacwi = min(nsacwi*mi0, psacws)
! Subprogram not used   else if((t < 268.16_r8) .and. (t >= 265.16_r8)) then
! Subprogram not used      nsacwi = 3.5e8_r8*(t-265.16_r8)/3.0_r8*psacws
! Subprogram not used      msacwi = min(nsacwi*mi0, psacws)
! Subprogram not used   else
! Subprogram not used      nsacwi = 0.0_r8
! Subprogram not used      msacwi = 0.0_r8
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   psacws = max(0.0_r8,psacws - nsacwi*mi0)
! Subprogram not used 
! Subprogram not used end subroutine secondary_ice_production

! accretion of rain water by snow
!===================================================================
! formula from ikawa and saito, 1991, used by reisner et al., 1998

! Subprogram not used elemental subroutine accrete_rain_snow(t, rho, umr, ums, unr, uns, qric, qsic, &
! Subprogram not used      lamr, n0r, lams, n0s, pracs, npracs )
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t   ! Temperature
! Subprogram not used   real(r8), intent(in) :: rho ! Density
! Subprogram not used 
! Subprogram not used   ! Fallspeeds
! Subprogram not used   ! mass-weighted
! Subprogram not used   real(r8), intent(in) :: umr ! rain
! Subprogram not used   real(r8), intent(in) :: ums ! snow
! Subprogram not used   ! number-weighted
! Subprogram not used   real(r8), intent(in) :: unr ! rain
! Subprogram not used   real(r8), intent(in) :: uns ! snow
! Subprogram not used 
! Subprogram not used   ! In cloud MMRs
! Subprogram not used   real(r8), intent(in) :: qric ! rain
! Subprogram not used   real(r8), intent(in) :: qsic ! snow
! Subprogram not used 
! Subprogram not used   ! Size distribution parameters
! Subprogram not used   ! rain
! Subprogram not used   real(r8), intent(in) :: lamr
! Subprogram not used   real(r8), intent(in) :: n0r
! Subprogram not used   ! snow
! Subprogram not used   real(r8), intent(in) :: lams
! Subprogram not used   real(r8), intent(in) :: n0s
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: pracs  ! MMR
! Subprogram not used   real(r8), intent(out) :: npracs ! Number
! Subprogram not used 
! Subprogram not used   ! Collection efficiency for accretion of rain by snow
! Subprogram not used   real(r8), parameter :: ecr = 1.0_r8
! Subprogram not used 
! Subprogram not used   if (qric >= icsmall .and. qsic >= icsmall .and. t <= tmelt) then
! Subprogram not used 
! Subprogram not used      pracs = pi*pi*ecr*(((1.2_r8*umr-0.95_r8*ums)**2 + &
! Subprogram not used           0.08_r8*ums*umr)**0.5_r8 *  &
! Subprogram not used           rhow * rho * n0r * n0s * &
! Subprogram not used           (5._r8/(lamr**6 * lams)+ &
! Subprogram not used           2._r8/(lamr**5 * lams**2)+ &
! Subprogram not used           0.5_r8/(lamr**4 * lams**3)))
! Subprogram not used 
! Subprogram not used      npracs = pi/2._r8*rho*ecr* (1.7_r8*(unr-uns)**2 + &
! Subprogram not used           0.3_r8*unr*uns)**0.5_r8 * &
! Subprogram not used           n0r*n0s* &
! Subprogram not used           (1._r8/(lamr**3 * lams)+ &
! Subprogram not used           1._r8/(lamr**2 * lams**2)+ &
! Subprogram not used           1._r8/(lamr * lams**3))
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used      pracs = 0._r8
! Subprogram not used      npracs = 0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine accrete_rain_snow

! heterogeneous freezing of rain drops
!===================================================================
! follows from Bigg (1953)

! Subprogram not used elemental subroutine heterogeneous_rain_freezing(t, qric, nric, lamr, mnuccr, nnuccr)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t    ! Temperature
! Subprogram not used 
! Subprogram not used   ! In-cloud rain
! Subprogram not used   real(r8), intent(in) :: qric ! MMR
! Subprogram not used   real(r8), intent(in) :: nric ! Number
! Subprogram not used   real(r8), intent(in) :: lamr ! size parameter
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: mnuccr ! MMR
! Subprogram not used   real(r8), intent(out) :: nnuccr ! Number
! Subprogram not used 
! Subprogram not used   if (t < 269.15_r8 .and. qric >= qsmall) then
! Subprogram not used 
! Subprogram not used      ! Division by lamr**3 twice is old workaround to avoid overflow.
! Subprogram not used      ! Probably no longer necessary
! Subprogram not used      mnuccr = 20._r8*pi*pi*rhow*nric*bimm* &
! Subprogram not used           (exp(aimm*(tmelt - t))-1._r8)/lamr**3 &
! Subprogram not used           /lamr**3
! Subprogram not used 
! Subprogram not used      nnuccr = pi*nric*bimm* &
! Subprogram not used           (exp(aimm*(tmelt - t))-1._r8)/lamr**3
! Subprogram not used   else
! Subprogram not used      mnuccr = 0._r8
! Subprogram not used      nnuccr = 0._r8
! Subprogram not used   end if
! Subprogram not used end subroutine heterogeneous_rain_freezing

! accretion of cloud liquid water by rain
!===================================================================
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

! Subprogram not used elemental subroutine accrete_cloud_water_rain(microp_uniform, qric, qcic, &
! Subprogram not used      ncic, relvar, accre_enhan, pra, npra)
! Subprogram not used 
! Subprogram not used   logical, intent(in) :: microp_uniform
! Subprogram not used 
! Subprogram not used   ! In-cloud rain
! Subprogram not used   real(r8), intent(in) :: qric ! MMR
! Subprogram not used 
! Subprogram not used   ! Cloud droplets
! Subprogram not used   real(r8), intent(in) :: qcic ! MMR
! Subprogram not used   real(r8), intent(in) :: ncic ! Number
! Subprogram not used 
! Subprogram not used   ! SGS variability
! Subprogram not used   real(r8), intent(in) :: relvar
! Subprogram not used   real(r8), intent(in) :: accre_enhan
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: pra  ! MMR
! Subprogram not used   real(r8), intent(out) :: npra ! Number
! Subprogram not used 
! Subprogram not used   ! Coefficient that varies for subcolumns
! Subprogram not used   real(r8) :: pra_coef
! Subprogram not used 
! Subprogram not used   if (.not. microp_uniform) then
! Subprogram not used      pra_coef = accre_enhan * var_coef(relvar, 1.15_r8)
! Subprogram not used   else
! Subprogram not used      pra_coef = 1._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   if (qric >= qsmall .and. qcic >= qsmall) then
! Subprogram not used 
! Subprogram not used      ! include sub-grid distribution of cloud water
! Subprogram not used      pra = pra_coef * 67._r8*(qcic*qric)**1.15_r8
! Subprogram not used 
! Subprogram not used      npra = pra/(qcic/ncic)
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used      pra = 0._r8
! Subprogram not used      npra = 0._r8
! Subprogram not used   end if
! Subprogram not used end subroutine accrete_cloud_water_rain

! Self-collection of rain drops
!===================================================================
! from Beheng(1994)

! Subprogram not used elemental subroutine self_collection_rain(rho, qric, nric, nragg)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: rho  ! Air density
! Subprogram not used 
! Subprogram not used   ! Rain
! Subprogram not used   real(r8), intent(in) :: qric ! MMR
! Subprogram not used   real(r8), intent(in) :: nric ! Number
! Subprogram not used 
! Subprogram not used   ! Output number tendency
! Subprogram not used   real(r8), intent(out) :: nragg
! Subprogram not used 
! Subprogram not used   if (qric >= qsmall) then
! Subprogram not used      nragg = -8._r8*nric*qric*rho
! Subprogram not used   else
! Subprogram not used      nragg = 0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine self_collection_rain

! Accretion of cloud ice by snow
!===================================================================
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

! Subprogram not used elemental subroutine accrete_cloud_ice_snow(t, rho, asn, qiic, niic, qsic, &
! Subprogram not used      lams, n0s, prai, nprai)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t    ! Temperature
! Subprogram not used   real(r8), intent(in) :: rho  ! Density
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: asn  ! Snow fallspeed parameter
! Subprogram not used 
! Subprogram not used   ! Cloud ice
! Subprogram not used   real(r8), intent(in) :: qiic ! MMR
! Subprogram not used   real(r8), intent(in) :: niic ! Number
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: qsic ! Snow MMR
! Subprogram not used 
! Subprogram not used   ! Snow size parameters
! Subprogram not used   real(r8), intent(in) :: lams
! Subprogram not used   real(r8), intent(in) :: n0s
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: prai  ! MMR
! Subprogram not used   real(r8), intent(out) :: nprai ! Number
! Subprogram not used 
! Subprogram not used   if (qsic >= qsmall .and. qiic >= qsmall .and. t <= tmelt) then
! Subprogram not used 
! Subprogram not used      prai = pi/4._r8 * asn * qiic * rho * n0s * eii * gamma_bs_plus3/ &
! Subprogram not used           lams**(bs+3._r8)
! Subprogram not used 
! Subprogram not used      nprai = pi/4._r8 * asn * niic * rho * n0s * eii * gamma_bs_plus3/ &
! Subprogram not used           lams**(bs+3._r8)
! Subprogram not used   else
! Subprogram not used      prai = 0._r8
! Subprogram not used      nprai = 0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine accrete_cloud_ice_snow

! calculate evaporation/sublimation of rain and snow
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

! Subprogram not used elemental subroutine evaporate_sublimate_precip(t, rho, dv, mu, sc, q, qvl, qvi, &
! Subprogram not used      lcldm, cldmax, arn, asn, qcic, qiic, qric, qsic, lamr, n0r, lams, n0s, &
! Subprogram not used      pre, prds)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t    ! temperature
! Subprogram not used   real(r8), intent(in) :: rho  ! air density
! Subprogram not used   real(r8), intent(in) :: dv   ! water vapor diffusivity
! Subprogram not used   real(r8), intent(in) :: mu   ! viscosity
! Subprogram not used   real(r8), intent(in) :: sc   ! schmidt number
! Subprogram not used   real(r8), intent(in) :: q    ! humidity
! Subprogram not used   real(r8), intent(in) :: qvl  ! saturation humidity (water)
! Subprogram not used   real(r8), intent(in) :: qvi  ! saturation humidity (ice)
! Subprogram not used   real(r8), intent(in) :: lcldm  ! liquid cloud fraction
! Subprogram not used   real(r8), intent(in) :: cldmax ! precipitation fraction (maximum overlap)
! Subprogram not used 
! Subprogram not used   ! fallspeed parameters
! Subprogram not used   real(r8), intent(in) :: arn  ! rain
! Subprogram not used   real(r8), intent(in) :: asn  ! snow
! Subprogram not used 
! Subprogram not used   ! In-cloud MMRs
! Subprogram not used   real(r8), intent(in) :: qcic ! cloud liquid
! Subprogram not used   real(r8), intent(in) :: qiic ! cloud ice
! Subprogram not used   real(r8), intent(in) :: qric ! rain
! Subprogram not used   real(r8), intent(in) :: qsic ! snow
! Subprogram not used 
! Subprogram not used   ! Size parameters
! Subprogram not used   ! rain
! Subprogram not used   real(r8), intent(in) :: lamr
! Subprogram not used   real(r8), intent(in) :: n0r
! Subprogram not used   ! snow
! Subprogram not used   real(r8), intent(in) :: lams
! Subprogram not used   real(r8), intent(in) :: n0s
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: pre
! Subprogram not used   real(r8), intent(out) :: prds
! Subprogram not used 
! Subprogram not used   real(r8) :: qclr   ! water vapor mixing ratio in clear air
! Subprogram not used   real(r8) :: ab     ! correction to account for latent heat
! Subprogram not used   real(r8) :: eps    ! 1/ sat relaxation timescale
! Subprogram not used 
! Subprogram not used   real(r8) :: dum
! Subprogram not used 
! Subprogram not used   ! set temporary cloud fraction to zero if cloud water + ice is very small
! Subprogram not used   ! this will ensure that evaporation/sublimation of precip occurs over
! Subprogram not used   ! entire grid cell, since min cloud fraction is specified otherwise
! Subprogram not used   if (qcic+qiic < 1.e-6_r8) then
! Subprogram not used      dum = 0._r8
! Subprogram not used   else
! Subprogram not used      dum = lcldm
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   ! only calculate if there is some precip fraction > cloud fraction
! Subprogram not used 
! Subprogram not used   if (cldmax > dum) then
! Subprogram not used 
! Subprogram not used      ! calculate q for out-of-cloud region
! Subprogram not used      qclr=(q-dum*qvl)/(1._r8-dum)
! Subprogram not used 
! Subprogram not used      ! evaporation of rain
! Subprogram not used      if (qric >= qsmall) then
! Subprogram not used 
! Subprogram not used         ab = calc_ab(t, qvl, xxlv)
! Subprogram not used         eps = 2._r8*pi*n0r*rho*Dv* &
! Subprogram not used              (f1r/(lamr*lamr)+ &
! Subprogram not used              f2r*(arn*rho/mu)**0.5_r8* &
! Subprogram not used              sc**(1._r8/3._r8)*gamma_half_br_plus5/ &
! Subprogram not used              (lamr**(5._r8/2._r8+br/2._r8)))
! Subprogram not used 
! Subprogram not used         pre = eps*(qclr-qvl)/ab
! Subprogram not used 
! Subprogram not used         ! only evaporate in out-of-cloud region
! Subprogram not used         ! and distribute across cldmax
! Subprogram not used         pre=min(pre*(cldmax-dum),0._r8)
! Subprogram not used         pre=pre/cldmax
! Subprogram not used      else
! Subprogram not used         pre = 0._r8
! Subprogram not used      end if
! Subprogram not used 
! Subprogram not used      ! sublimation of snow
! Subprogram not used      if (qsic >= qsmall) then
! Subprogram not used         ab = calc_ab(t, qvi, xxls)
! Subprogram not used         eps = 2._r8*pi*n0s*rho*Dv* &
! Subprogram not used              (f1s/(lams*lams)+ &
! Subprogram not used              f2s*(asn*rho/mu)**0.5_r8* &
! Subprogram not used              sc**(1._r8/3._r8)*gamma_half_bs_plus5/ &
! Subprogram not used              (lams**(5._r8/2._r8+bs/2._r8)))
! Subprogram not used         prds = eps*(qclr-qvi)/ab
! Subprogram not used 
! Subprogram not used         ! only sublimate in out-of-cloud region and distribute over cldmax
! Subprogram not used         prds=min(prds*(cldmax-dum),0._r8)
! Subprogram not used         prds=prds/cldmax
! Subprogram not used      else
! Subprogram not used         prds = 0._r8
! Subprogram not used      end if
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used      prds = 0._r8
! Subprogram not used      pre = 0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine evaporate_sublimate_precip

! bergeron process - evaporation of droplets and deposition onto snow
!===================================================================

! Subprogram not used elemental subroutine bergeron_process_snow(t, rho, dv, mu, sc, qvl, qvi, asn, &
! Subprogram not used      qcic, qsic, lams, n0s, bergs)
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: t    ! temperature
! Subprogram not used   real(r8), intent(in) :: rho  ! air density
! Subprogram not used   real(r8), intent(in) :: dv   ! water vapor diffusivity
! Subprogram not used   real(r8), intent(in) :: mu   ! viscosity
! Subprogram not used   real(r8), intent(in) :: sc   ! schmidt number
! Subprogram not used   real(r8), intent(in) :: qvl  ! saturation humidity (water)
! Subprogram not used   real(r8), intent(in) :: qvi  ! saturation humidity (ice)
! Subprogram not used 
! Subprogram not used   ! fallspeed parameter for snow
! Subprogram not used   real(r8), intent(in) :: asn
! Subprogram not used 
! Subprogram not used   ! In-cloud MMRs
! Subprogram not used   real(r8), intent(in) :: qcic ! cloud liquid
! Subprogram not used   real(r8), intent(in) :: qsic ! snow
! Subprogram not used 
! Subprogram not used   ! Size parameters for snow
! Subprogram not used   real(r8), intent(in) :: lams
! Subprogram not used   real(r8), intent(in) :: n0s
! Subprogram not used 
! Subprogram not used   ! Output tendencies
! Subprogram not used   real(r8), intent(out) :: bergs
! Subprogram not used 
! Subprogram not used   real(r8) :: ab     ! correction to account for latent heat
! Subprogram not used   real(r8) :: eps    ! 1/ sat relaxation timescale
! Subprogram not used 
! Subprogram not used   if (qsic >= qsmall.and. qcic >= qsmall .and. t < tmelt) then
! Subprogram not used      ab = calc_ab(t, qvi, xxls)
! Subprogram not used      eps = 2._r8*pi*n0s*rho*Dv* &
! Subprogram not used           (f1s/(lams*lams)+ &
! Subprogram not used           f2s*(asn*rho/mu)**0.5_r8* &
! Subprogram not used           sc**(1._r8/3._r8)*gamma_half_bs_plus5/ &
! Subprogram not used           (lams**(5._r8/2._r8+bs/2._r8)))
! Subprogram not used      bergs = eps*(qvl-qvi)/ab
! Subprogram not used   else
! Subprogram not used      bergs = 0._r8
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine bergeron_process_snow

!========================================================================
!UTILITIES
!========================================================================

pure function no_limiter()
  real(r8) :: no_limiter

  no_limiter = transfer(limiter_off, no_limiter)

end function no_limiter

pure function limiter_is_on(lim)
  real(r8), intent(in) :: lim
  logical :: limiter_is_on

  limiter_is_on = transfer(lim, limiter_off) /= limiter_off

end function limiter_is_on

end module micro_mg_utils
