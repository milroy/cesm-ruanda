
module cldwat
!----------------------------------------------------------------------- 
! 
! Purpose: Prognostic cloud water data and methods.
! 
! Public interfaces:
!
! inimc -- Initialize constants
! pcond -- Calculate prognostic condensate
!
! Author: P. Rasch, with Modifications by Minghua Zhang
! January 2010, modified by J. Kay to add precip fluxes for COSP simulator
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use spmd_utils,    only: masterproc
   use ppgrid,        only: pcols, pver, pverp
   use physconst,     only: latvap, latice, cpair
   use abortutils,    only: endrun
   use cam_logfile,   only: iulog
   use ref_pres,      only: top_lev => trop_cloud_top_lev

   implicit none

!-----------------------------------------------------------------------
! PUBLIC: Make default data and interfaces private
!-----------------------------------------------------------------------
   private
   save
   public inimc, pcond          ! Public interfaces
   public cldwat_readnl
   integer, public::  ktop      ! Level above 10 hPa

   real(r8),public ::  icritc               ! threshold for autoconversion of cold ice
   real(r8),public ::  icritw               ! threshold for autoconversion of warm ice
!!$   real(r8),public,parameter::  conke  = 1.e-6    ! tunable constant for evaporation of precip
!!$   real(r8),public,parameter::  conke  =  2.e-6    ! tunable constant for evaporation of precip
   real(r8),public ::  conke                ! tunable constant for evaporation of precip
   real(r8),public ::  r3lcrit              ! critical radius where liq conversion begins

!-----------------------------------------------------------------------
! PRIVATE: Everything else is private to this module
!-----------------------------------------------------------------------
   real(r8), private:: rhonot   ! air density at surface
   real(r8), private:: t0       ! Freezing temperature
   real(r8), private:: cldmin   ! assumed minimum cloud amount
   real(r8), private:: small    ! small number compared to unity
   real(r8), private:: c        ! constant for graupel like snow cm**(1-d)/s
   real(r8), private:: d        ! constant for graupel like snow
   real(r8), private:: esi      ! collection efficient for ice by snow
   real(r8), private:: esw      ! collection efficient for water by snow
   real(r8), private:: nos      ! particles snow / cm**4
   real(r8), private:: pi       ! Mathematical constant
   real(r8), private:: gravit   ! Gravitational acceleration at surface
   real(r8), private:: rh2o
   real(r8), private:: prhonos
   real(r8), private:: thrpd    ! numerical three added to d
   real(r8), private:: gam3pd   ! gamma function on (3+d)
   real(r8), private:: gam4pd   ! gamma function on (4+d)
   real(r8), private:: rhoi     ! ice density
   real(r8), private:: rhos     ! snow density
   real(r8), private:: rhow     ! water density
   real(r8), private:: mcon01   ! constants used in cloud microphysics
   real(r8), private:: mcon02   ! constants used in cloud microphysics
   real(r8), private:: mcon03   ! constants used in cloud microphysics
   real(r8), private:: mcon04   ! constants used in cloud microphysics
   real(r8), private:: mcon05   ! constants used in cloud microphysics
   real(r8), private:: mcon06   ! constants used in cloud microphysics
   real(r8), private:: mcon07   ! constants used in cloud microphysics
   real(r8), private:: mcon08   ! constants used in cloud microphysics


! Parameters used in findmcnew
   real(r8) :: capnsi               ! sea ice cloud particles / cm3
   real(r8) :: capnc                ! cold and oceanic cloud particles / cm3
   real(r8) :: capnw                ! warm continental cloud particles / cm3
   real(r8) :: kconst               ! const for terminal velocity (stokes regime)
   real(r8) :: effc                 ! collection efficiency
   real(r8) :: alpha                ! ratio of 3rd moment radius to 2nd
   real(r8) :: capc                 ! constant for autoconversion
   real(r8) :: convfw               ! constant used for fall velocity calculation
   real(r8) :: cracw                ! constant used for rain accreting water
   real(r8) :: critpr               ! critical precip rate collection efficiency changes
   real(r8) :: ciautb               ! coefficient of autoconversion of ice (1/s)








  ! Private data
  real(r8), parameter :: unset_r8 = huge(1.0_r8)

contains
!===============================================================================
  subroutine cldwat_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input




   ! Namelist variables
   real(r8) :: cldwat_icritw  = unset_r8    !   icritw  = threshold for autoconversion of warm ice  
   real(r8) :: cldwat_icritc  = unset_r8    !   icritc  = threshold for autoconversion of cold ice  
   real(r8) :: cldwat_conke   = unset_r8    !   conke   = tunable constant for evaporation of precip
   real(r8) :: cldwat_r3lcrit = unset_r8    !   r3lcrit = critical radius where liq conversion begins

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cldwat_readnl'

   namelist /cldwat_nl/ cldwat_icritw, cldwat_icritc, cldwat_conke, cldwat_r3lcrit

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldwat_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldwat_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      icritw  = cldwat_icritw 
      icritc  = cldwat_icritc
      conke   = cldwat_conke
      r3lcrit = cldwat_r3lcrit

   end if




   ! Broadcast namelist variables
   call mpibcast(icritw,            1, mpir8,  0, mpicom)
   call mpibcast(icritc,            1, mpir8,  0, mpicom)
   call mpibcast(conke,             1, mpir8,  0, mpicom)
   call mpibcast(r3lcrit,           1, mpir8,  0, mpicom)


end subroutine cldwat_readnl

! Subprogram not used subroutine inimc( tmeltx, rhonotx, gravitx, rh2ox)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! initialize constants for the prognostic condensate
! Subprogram not used ! 
! Subprogram not used ! Author: P. Rasch, April 1997
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use pmgrid,       only: plev, plevp
! Subprogram not used    use dycore,       only: dycore_is, get_resolution
! Subprogram not used    use ref_pres,     only: pref_mid
! Subprogram not used 
! Subprogram not used    integer k
! Subprogram not used    real(r8), intent(in) :: tmeltx
! Subprogram not used    real(r8), intent(in) :: rhonotx
! Subprogram not used    real(r8), intent(in) :: gravitx
! Subprogram not used    real(r8), intent(in) :: rh2ox
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    rhonot = rhonotx             ! air density at surface (gm/cm3)
! Subprogram not used    gravit = gravitx
! Subprogram not used    rh2o   = rh2ox
! Subprogram not used    rhos = .1_r8                 ! assumed snow density (gm/cm3)
! Subprogram not used    rhow = 1._r8                 ! water density
! Subprogram not used    rhoi = 1._r8                 ! ice density
! Subprogram not used    esi = 1.0_r8                 ! collection efficient for ice by snow
! Subprogram not used    esw = 0.1_r8                 ! collection efficient for water by snow
! Subprogram not used    t0 = tmeltx                  ! approximate freezing temp
! Subprogram not used    cldmin = 0.02_r8             ! assumed minimum cloud amount
! Subprogram not used    small = 1.e-22_r8            ! a small number compared to unity
! Subprogram not used    c = 152.93_r8                ! constant for graupel like snow cm**(1-d)/s
! Subprogram not used    d = 0.25_r8                  ! constant for graupel like snow
! Subprogram not used    nos = 3.e-2_r8               ! particles snow / cm**4
! Subprogram not used    pi = 4._r8*atan(1.0_r8)
! Subprogram not used    prhonos = pi*rhos*nos
! Subprogram not used    thrpd = 3._r8 + d
! Subprogram not used    if (d==0.25_r8) then
! Subprogram not used       gam3pd = 2.549256966718531_r8 ! only right for d = 0.25
! Subprogram not used       gam4pd = 8.285085141835282_r8
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call endrun(' can only use d ne 0.25 on a cray ')
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used    mcon01 = pi*nos*c*gam3pd/4._r8
! Subprogram not used    mcon02 = 1._r8/(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4._r8)))
! Subprogram not used    mcon03 = -(0.5_r8+d/4._r8)
! Subprogram not used    mcon04 = 4._r8/(4._r8+d)
! Subprogram not used    mcon05 = (3+d)/(4+d)
! Subprogram not used    mcon06 = (3+d)/4._r8
! Subprogram not used    mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
! Subprogram not used    mcon08 = -0.5_r8/(4._r8+d)
! Subprogram not used 
! Subprogram not used    if( masterproc ) write(iulog,*) 'cloud water initialization by inimc complete '
! Subprogram not used 
! Subprogram not used ! Initialize parameters used by findmcnew
! Subprogram not used    capnw = 400._r8              ! warm continental cloud particles / cm3
! Subprogram not used    capnc = 150._r8              ! cold and oceanic cloud particles / cm3
! Subprogram not used !  capnsi = 40._r8              ! sea ice cloud particles density  / cm3
! Subprogram not used    capnsi = 75._r8              ! sea ice cloud particles density  / cm3
! Subprogram not used 
! Subprogram not used    kconst = 1.18e6_r8           ! const for terminal velocity
! Subprogram not used 
! Subprogram not used !  effc = 1._r8                 ! autoconv collection efficiency following boucher 96
! Subprogram not used !  effc = .55*0.05_r8           ! autoconv collection efficiency following baker 93
! Subprogram not used    effc = 0.55_r8               ! autoconv collection efficiency following tripoli and cotton
! Subprogram not used !  effc = 0._r8    ! turn off warm-cloud autoconv
! Subprogram not used    alpha = 1.1_r8**4
! Subprogram not used    capc = pi**(-.333_r8)*kconst*effc *(0.75_r8)**(1.333_r8)*alpha  ! constant for autoconversion
! Subprogram not used 
! Subprogram not used ! critical precip rate at which we assume the collector drops can change the
! Subprogram not used ! drop size enough to enhance the auto-conversion process (mm/day)
! Subprogram not used    critpr = 0.5_r8
! Subprogram not used    convfw = 1.94_r8*2.13_r8*sqrt(rhow*1000._r8*9.81_r8*2.7e-4_r8)
! Subprogram not used 
! Subprogram not used ! liquid microphysics
! Subprogram not used !  cracw = 6_r8                 ! beheng
! Subprogram not used    cracw = .884_r8*sqrt(9.81_r8/(rhow*1000._r8*2.7e-4_r8)) ! tripoli and cotton
! Subprogram not used 
! Subprogram not used ! ice microphysics
! Subprogram not used    ciautb = 5.e-4_r8
! Subprogram not used 
! Subprogram not used    if ( masterproc ) then
! Subprogram not used       write(iulog,*)'tuning parameters cldwat: icritw',icritw,'icritc',icritc,'conke',conke,'r3lcrit',r3lcrit
! Subprogram not used       write(iulog,*)'tuning parameters cldwat: capnw',capnw,'capnc',capnc,'capnsi',capnsi,'kconst',kconst
! Subprogram not used       write(iulog,*)'tuning parameters cldwat: effc',effc,'alpha',alpha,'capc',capc
! Subprogram not used       write(iulog,*)'tuning parameters cldwat: critpr',critpr,'convfw',convfw,'cracw',cracw,'ciautb',ciautb
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used end subroutine inimc

! Subprogram not used subroutine pcond (lchnk   ,ncol    , &
! Subprogram not used                   tn      ,ttend   ,qn      ,qtend   ,omega   , &
! Subprogram not used                   cwat    ,p       ,pdel    ,cldn    ,fice    , fsnow, &
! Subprogram not used                   cme     ,prodprec,prodsnow,evapprec,evapsnow,evapheat, prfzheat, &     
! Subprogram not used                   meltheat,precip  ,snowab  ,deltat  ,fwaut   , &
! Subprogram not used                   fsaut   ,fracw   ,fsacw   ,fsaci   ,lctend  , &
! Subprogram not used                   rhdfda  ,rhu00   ,seaicef, zi      ,ice2pr, liq2pr, &
! Subprogram not used                   liq2snow, snowh, rkflxprc, rkflxsnw, pracwo, psacwo, psacio)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! The public interface to the cloud water parameterization
! Subprogram not used ! returns tendencies to water vapor, temperature and cloud water variables
! Subprogram not used ! 
! Subprogram not used ! For basic method 
! Subprogram not used !  See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
! Subprogram not used !  model climate using diagnosed and 
! Subprogram not used !  predicted condensate parameterizations, 1998, J. Clim., 11,
! Subprogram not used !  pp1587---1614.
! Subprogram not used ! 
! Subprogram not used ! For important modifications to improve the method of determining
! Subprogram not used ! condensation/evaporation see Zhang et al (2001, in preparation)
! Subprogram not used !
! Subprogram not used ! Authors: M. Zhang, W. Lin, P. Rasch and J.E. Kristjansson
! Subprogram not used !          B. A. Boville (latent heat of fusion)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use wv_saturation, only: qsat, estblf, svp_to_qsat, findsp_vc
! Subprogram not used    use physconst, only: epsilo
! Subprogram not used    use cam_control_mod, only: nlvdry
! Subprogram not used !
! Subprogram not used !---------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Input Arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: lchnk                 ! chunk identifier
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: fice(pcols,pver)     ! fraction of cwat that is ice
! Subprogram not used    real(r8), intent(in) :: fsnow(pcols,pver)    ! fraction of rain that freezes to snow
! Subprogram not used    real(r8), intent(in) :: cldn(pcols,pver)     ! new value of cloud fraction    (fraction)
! Subprogram not used    real(r8), intent(in) :: cwat(pcols,pver)     ! cloud water (kg/kg)
! Subprogram not used    real(r8), intent(in) :: omega(pcols,pver)    ! vert pressure vel (Pa/s)
! Subprogram not used    real(r8), intent(in) :: p(pcols,pver)        ! pressure          (K)
! Subprogram not used    real(r8), intent(in) :: pdel(pcols,pver)     ! pressure thickness (Pa)
! Subprogram not used    real(r8), intent(in) :: qn(pcols,pver)       ! new water vapor    (kg/kg)
! Subprogram not used    real(r8), intent(in) :: qtend(pcols,pver)    ! mixing ratio tend  (kg/kg/s)
! Subprogram not used    real(r8), intent(in) :: tn(pcols,pver)       ! new temperature    (K)
! Subprogram not used    real(r8), intent(in) :: ttend(pcols,pver)    ! temp tendencies    (K/s)
! Subprogram not used    real(r8), intent(in) :: deltat               ! time step to advance solution over
! Subprogram not used    real(r8), intent(in) :: lctend(pcols,pver)   ! cloud liquid water tendencies   ====wlin
! Subprogram not used    real(r8), intent(in) :: rhdfda(pcols,pver)   ! dG(a)/da, rh=G(a), when rh>u00  ====wlin
! Subprogram not used    real(r8), intent(in) :: rhu00 (pcols,pver)   ! Rhlim for cloud                 ====wlin
! Subprogram not used    real(r8), intent(in) :: seaicef(pcols)       ! sea ice fraction  (fraction)
! Subprogram not used    real(r8), intent(in) :: zi(pcols,pverp)      ! layer interfaces (m)
! Subprogram not used     real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
! Subprogram not used !
! Subprogram not used ! Output Arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: cme     (pcols,pver) ! rate of cond-evap of condensate (1/s)
! Subprogram not used    real(r8), intent(out) :: prodprec(pcols,pver) ! rate of conversion of condensate to precip (1/s)
! Subprogram not used    real(r8), intent(out) :: evapprec(pcols,pver) ! rate of evaporation of falling precip (1/s)
! Subprogram not used    real(r8), intent(out) :: evapsnow(pcols,pver) ! rate of evaporation of falling snow (1/s)
! Subprogram not used    real(r8), intent(out) :: evapheat(pcols,pver) ! heating rate due to evaporation of precip (W/kg)
! Subprogram not used    real(r8), intent(out) :: prfzheat(pcols,pver) ! heating rate due to freezing of precip (W/kg)
! Subprogram not used    real(r8), intent(out) :: meltheat(pcols,pver) ! heating rate due to snow melt (W/kg)
! Subprogram not used    real(r8), intent(out) :: precip(pcols)        ! rate of precipitation (kg / (m**2 * s))
! Subprogram not used    real(r8), intent(out) :: snowab(pcols)        ! rate of snow (kg / (m**2 * s))
! Subprogram not used    real(r8), intent(out) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
! Subprogram not used    real(r8), intent(out) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
! Subprogram not used    real(r8), intent(out) :: liq2snow(pcols,pver) ! rate of conversion of liquid to snow
! Subprogram not used    real(r8), intent(out) :: rkflxprc(pcols,pverp)   ! grid-box mean RK flux_large_scale_cloud_rain+snow at interfaces (kg m^-2 s^-1)
! Subprogram not used    real(r8), intent(out) :: rkflxsnw(pcols,pverp)   ! grid-box mean RK flux_large_scale_cloud_snow at interfaces (kg m^-2 s^-1)
! Subprogram not used ! intent(out)s here for pcond to pass to stratiform.F90 to be addflded/outflded
! Subprogram not used    real(r8), intent(out) :: pracwo(pcols,pver)      ! accretion of cloud water by rain (1/s)
! Subprogram not used    real(r8), intent(out) :: psacwo(pcols,pver)      ! accretion of cloud water by snow (1/s)
! Subprogram not used    real(r8), intent(out) :: psacio(pcols,pver)      ! accretion of cloud ice by snow (1/s)
! Subprogram not used 
! Subprogram not used    real(r8) nice2pr     ! rate of conversion of ice to snow
! Subprogram not used    real(r8) nliq2pr     ! rate of conversion of liquid to precip
! Subprogram not used    real(r8) nliq2snow   ! rate of conversion of liquid to snow
! Subprogram not used    real(r8) prodsnow(pcols,pver) ! rate of production of snow
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used    real(r8) :: precab(pcols)        ! rate of precipitation (kg / (m**2 * s))
! Subprogram not used    integer i                 ! work variable
! Subprogram not used    integer iter              ! #iterations for precipitation calculation
! Subprogram not used    integer k                 ! work variable
! Subprogram not used    integer l                 ! work variable
! Subprogram not used 
! Subprogram not used    real(r8) cldm(pcols)          ! mean cloud fraction over the time step
! Subprogram not used    real(r8) cldmax(pcols)        ! max cloud fraction above
! Subprogram not used    real(r8) coef(pcols)          ! conversion time scale for condensate to rain
! Subprogram not used    real(r8) cwm(pcols)           ! cwat mixing ratio at midpoint of time step
! Subprogram not used    real(r8) cwn(pcols)           ! cwat mixing ratio at end
! Subprogram not used    real(r8) denom                ! work variable
! Subprogram not used    real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
! Subprogram not used    real(r8) es(pcols)            ! sat. vapor pressure
! Subprogram not used    real(r8) fracw(pcols,pver)    ! relative importance of collection of liquid by rain
! Subprogram not used    real(r8) fsaci(pcols,pver)    ! relative importance of collection of ice by snow
! Subprogram not used    real(r8) fsacw(pcols,pver)    ! relative importance of collection of liquid by snow
! Subprogram not used    real(r8) fsaut(pcols,pver)    ! relative importance of ice auto conversion
! Subprogram not used    real(r8) fwaut(pcols,pver)    ! relative importance of warm cloud autoconversion
! Subprogram not used    real(r8) gamma(pcols)         ! d qs / dT
! Subprogram not used    real(r8) icwc(pcols)          ! in-cloud water content (kg/kg)
! Subprogram not used    real(r8) mincld               ! a small cloud fraction to avoid / zero
! Subprogram not used    real(r8),parameter ::omsm=0.99999_r8                 ! a number just less than unity (for rounding)
! Subprogram not used    real(r8) prprov(pcols)        ! provisional value of precip at btm of layer
! Subprogram not used    real(r8) prtmp                ! work variable
! Subprogram not used    real(r8) q(pcols,pver)        ! mixing ratio before time step ignoring condensate
! Subprogram not used    real(r8) qs(pcols)            ! spec. hum. of water vapor
! Subprogram not used    real(r8) qsn, esn             ! work variable
! Subprogram not used    real(r8) qsp(pcols,pver)      ! sat pt mixing ratio
! Subprogram not used    real(r8) qtl(pcols)           ! tendency which would saturate the grid box in deltat
! Subprogram not used    real(r8) qtmp, ttmp           ! work variable
! Subprogram not used    real(r8) relhum1(pcols)        ! relative humidity
! Subprogram not used    real(r8) relhum(pcols)        ! relative humidity
! Subprogram not used !!$   real(r8) tc                   ! crit temp of transition to ice
! Subprogram not used    real(r8) t(pcols,pver)        ! temp before time step ignoring condensate
! Subprogram not used    real(r8) tsp(pcols,pver)      ! sat pt temperature
! Subprogram not used    real(r8) pol                  ! work variable
! Subprogram not used    real(r8) cdt                  ! work variable
! Subprogram not used    real(r8) wtthick              ! work variable
! Subprogram not used 
! Subprogram not used ! Extra local work space for cloud scheme modification       
! Subprogram not used 
! Subprogram not used    real(r8) cpohl                !cpair/Latvap
! Subprogram not used    real(r8) hlocp                !Latvap/cpair
! Subprogram not used    real(r8) dto2                 !0.5*deltat (delta=2.0*dt)
! Subprogram not used    real(r8) calpha(pcols)        !alpha of new C - E scheme formulation
! Subprogram not used    real(r8) cbeta (pcols)        !beta  of new C - E scheme formulation
! Subprogram not used    real(r8) cbetah(pcols)        !beta_hat at saturation portion 
! Subprogram not used    real(r8) cgamma(pcols)        !gamma of new C - E scheme formulation
! Subprogram not used    real(r8) cgamah(pcols)        !gamma_hat at saturation portion
! Subprogram not used    real(r8) rcgama(pcols)        !gamma/gamma_hat
! Subprogram not used    real(r8) csigma(pcols)        !sigma of new C - E scheme formulation
! Subprogram not used    real(r8) cmec1 (pcols)        !c1    of new C - E scheme formulation
! Subprogram not used    real(r8) cmec2 (pcols)        !c2    of new C - E scheme formulation
! Subprogram not used    real(r8) cmec3 (pcols)        !c3    of new C - E scheme formulation
! Subprogram not used    real(r8) cmec4 (pcols)        !c4    of new C - E scheme formulation
! Subprogram not used    real(r8) cmeres(pcols)        !residual cond of over-sat after cme and evapprec
! Subprogram not used    real(r8) ctmp                 !a scalar representation of cmeres
! Subprogram not used    real(r8) clrh2o               ! Ratio of latvap to water vapor gas const
! Subprogram not used    real(r8) ice(pcols,pver)    ! ice mixing ratio
! Subprogram not used    real(r8) liq(pcols,pver)    ! liquid mixing ratio
! Subprogram not used    real(r8) rcwn(pcols,2,pver), rliq(pcols,2,pver), rice(pcols,2,pver)
! Subprogram not used    real(r8) cwnsave(pcols,2,pver), cmesave(pcols,2,pver)
! Subprogram not used    real(r8) prodprecsave(pcols,2,pver)
! Subprogram not used    logical error_found
! Subprogram not used !
! Subprogram not used !------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    clrh2o = latvap/rh2o   ! Ratio of latvap to water vapor gas const
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    mincld = 1.e-4_r8
! Subprogram not used    iter = 2
! Subprogram not used 
! Subprogram not used !   omsm = 0.99999
! Subprogram not used    cpohl = cpair/latvap
! Subprogram not used    hlocp = latvap/cpair
! Subprogram not used    dto2=0.5_r8*deltat
! Subprogram not used !
! Subprogram not used ! Constant for computing rate of evaporation of precipitation:
! Subprogram not used !
! Subprogram not used !!$   conke = 1.e-5
! Subprogram not used !!$   conke = 1.e-6
! Subprogram not used !
! Subprogram not used ! initialize a few single level fields
! Subprogram not used !
! Subprogram not used    do i = 1,ncol
! Subprogram not used       precip(i) = 0.0_r8
! Subprogram not used       precab(i) = 0.0_r8
! Subprogram not used       snowab(i) = 0.0_r8
! Subprogram not used       cldmax(i) = 0.0_r8
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! initialize multi-level fields 
! Subprogram not used !
! Subprogram not used    do k = 1,pver
! Subprogram not used       do i = 1,ncol
! Subprogram not used          q(i,k) = qn(i,k) 
! Subprogram not used          t(i,k) = tn(i,k)
! Subprogram not used !         q(i,k)=qn(i,k)-qtend(i,k)*deltat
! Subprogram not used !         t(i,k)=tn(i,k)-ttend(i,k)*deltat
! Subprogram not used     end do
! Subprogram not used    end do
! Subprogram not used    cme     (:ncol,:) = 0._r8
! Subprogram not used    evapprec(:ncol,:) = 0._r8
! Subprogram not used    prodprec(:ncol,:) = 0._r8
! Subprogram not used    evapsnow(:ncol,:) = 0._r8
! Subprogram not used    prodsnow(:ncol,:) = 0._r8
! Subprogram not used    evapheat(:ncol,:) = 0._r8
! Subprogram not used    meltheat(:ncol,:) = 0._r8
! Subprogram not used    prfzheat(:ncol,:) = 0._r8
! Subprogram not used    ice2pr(:ncol,:)   = 0._r8
! Subprogram not used    liq2pr(:ncol,:)   = 0._r8
! Subprogram not used    liq2snow(:ncol,:) = 0._r8
! Subprogram not used    fwaut(:ncol,:) = 0._r8
! Subprogram not used    fsaut(:ncol,:) = 0._r8
! Subprogram not used    fracw(:ncol,:) = 0._r8
! Subprogram not used    fsacw(:ncol,:) = 0._r8
! Subprogram not used    fsaci(:ncol,:) = 0._r8
! Subprogram not used    rkflxprc(:ncol,:) = 0._r8
! Subprogram not used    rkflxsnw(:ncol,:) = 0._r8
! Subprogram not used 
! Subprogram not used    pracwo(:ncol,:) = 0._r8
! Subprogram not used    psacwo(:ncol,:) = 0._r8
! Subprogram not used    psacio(:ncol,:) = 0._r8
! Subprogram not used !
! Subprogram not used ! find the wet bulb temp and saturation value
! Subprogram not used ! for the provisional t and q without condensation
! Subprogram not used !
! Subprogram not used    do 800 k = top_lev,pver
! Subprogram not used 
! Subprogram not used       ! "True" means that ice will be taken into account.
! Subprogram not used       call findsp_vc(qn(:ncol,k), tn(:ncol,k), p(:ncol,k), .true., &
! Subprogram not used            tsp(:ncol,k), qsp(:ncol,k))
! Subprogram not used 
! Subprogram not used       call qsat(t(:ncol,k), p(:ncol,k), &
! Subprogram not used            es(:ncol), qs(:ncol), gam=gamma(:ncol))
! Subprogram not used       do i = 1,ncol
! Subprogram not used          relhum(i) = q(i,k)/qs(i)
! Subprogram not used !
! Subprogram not used          cldm(i) = max(cldn(i,k),mincld)
! Subprogram not used !
! Subprogram not used ! the max cloud fraction above this level
! Subprogram not used !
! Subprogram not used          cldmax(i) = max(cldmax(i), cldm(i))
! Subprogram not used 
! Subprogram not used ! define the coefficients for C - E calculation
! Subprogram not used 
! Subprogram not used          calpha(i) = 1.0_r8/qs(i)
! Subprogram not used          cbeta (i) = q(i,k)/qs(i)**2*gamma(i)*cpohl
! Subprogram not used          cbetah(i) = 1.0_r8/qs(i)*gamma(i)*cpohl
! Subprogram not used          cgamma(i) = calpha(i)+latvap*cbeta(i)/cpair
! Subprogram not used          cgamah(i) = calpha(i)+latvap*cbetah(i)/cpair
! Subprogram not used          rcgama(i) = cgamma(i)/cgamah(i)
! Subprogram not used 
! Subprogram not used          if(cldm(i) > mincld) then
! Subprogram not used             icwc(i) = max(0._r8,cwat(i,k)/cldm(i))
! Subprogram not used          else
! Subprogram not used             icwc(i) = 0.0_r8
! Subprogram not used          endif
! Subprogram not used !PJR the above logic give zero icwc with nonzero cwat, dont like it!
! Subprogram not used !PJR generates problems with csigma
! Subprogram not used !PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
! Subprogram not used !         icwc(i) = max(1.e-8_r8,cwat(i,k)/cldm(i))
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! initial guess of evaporation, will be updated within iteration
! Subprogram not used !
! Subprogram not used          evapprec(i,k) = conke*(1._r8 - cldm(i))*sqrt(precab(i)) &
! Subprogram not used                         *(1._r8 - min(relhum(i),1._r8))
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! zero cmeres before iteration for each level
! Subprogram not used !
! Subprogram not used          cmeres(i)=0.0_r8
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used       do i = 1,ncol
! Subprogram not used !
! Subprogram not used ! fractions of ice at this level
! Subprogram not used !
! Subprogram not used !!$         tc = t(i,k) - t0
! Subprogram not used !!$         fice(i,k) = max(0._r8,min(-tc*0.05,1.0_r8))
! Subprogram not used !
! Subprogram not used ! calculate the cooling due to a phase change of the rainwater
! Subprogram not used ! from above
! Subprogram not used !
! Subprogram not used          if (t(i,k) >= t0) then
! Subprogram not used             meltheat(i,k) =  -latice * snowab(i) * gravit/pdel(i,k)
! Subprogram not used             snowab(i) = 0._r8
! Subprogram not used          else
! Subprogram not used             meltheat(i,k) = 0._r8
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! calculate cme and formation of precip. 
! Subprogram not used !
! Subprogram not used ! The cloud microphysics is highly nonlinear and coupled with cme
! Subprogram not used ! Both rain processes and cme are calculated iteratively.
! Subprogram not used ! 
! Subprogram not used       do 100 l = 1,iter
! Subprogram not used 
! Subprogram not used         do i = 1,ncol
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! calculation of cme has 4 scenarios
! Subprogram not used ! ==================================
! Subprogram not used !
! Subprogram not used            if(relhum(i) > rhu00(i,k)) then
! Subprogram not used     
! Subprogram not used            ! 1. whole grid saturation
! Subprogram not used            ! ========================
! Subprogram not used               if(relhum(i) >= 0.999_r8 .or. cldm(i) >= 0.999_r8 ) then
! Subprogram not used                  cme(i,k)=(calpha(i)*qtend(i,k)-cbetah(i)*ttend(i,k))/cgamah(i)
! Subprogram not used 
! Subprogram not used            ! 2. fractional saturation
! Subprogram not used            ! ========================
! Subprogram not used               else
! Subprogram not used                  if (rhdfda(i,k) .eq. 0._r8 .and. icwc(i).eq.0._r8) then
! Subprogram not used                     write (iulog,*) ' cldwat.F90:  empty rh cloud ', i, k, lchnk
! Subprogram not used                     write (iulog,*) ' relhum, iter ', relhum(i), l, rhu00(i,k), cldm(i), cldn(i,k)
! Subprogram not used                     call endrun ()
! Subprogram not used                  endif
! Subprogram not used                   csigma(i) = 1.0_r8/(rhdfda(i,k)+cgamma(i)*icwc(i))
! Subprogram not used                   cmec1(i) = (1.0_r8-cldm(i))*csigma(i)*rhdfda(i,k)
! Subprogram not used                   cmec2(i) = cldm(i)*calpha(i)/cgamah(i)+(1.0_r8-rcgama(i)*cldm(i))*   &
! Subprogram not used                              csigma(i)*calpha(i)*icwc(i)
! Subprogram not used                   cmec3(i) = cldm(i)*cbetah(i)/cgamah(i) +  &
! Subprogram not used                            (cbeta(i)-rcgama(i)*cldm(i)*cbetah(i))*csigma(i)*icwc(i)
! Subprogram not used                   cmec4(i) = csigma(i)*cgamma(i)*icwc(i)
! Subprogram not used 
! Subprogram not used                   ! Q=C-E=-C1*Al + C2*Aq - C3* At + C4*Er
! Subprogram not used   
! Subprogram not used                   cme(i,k) = -cmec1(i)*lctend(i,k) + cmec2(i)*qtend(i,k)  &
! Subprogram not used                              -cmec3(i)*ttend(i,k) + cmec4(i)*evapprec(i,k)
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used            ! 3. when rh < rhu00, evaporate existing cloud water
! Subprogram not used            ! ================================================== 
! Subprogram not used            else if(cwat(i,k) > 0.0_r8)then
! Subprogram not used               ! liquid water should be evaporated but not to exceed 
! Subprogram not used               ! saturation point. if qn > qsp, not to evaporate cwat
! Subprogram not used               cme(i,k)=-min(max(0._r8,qsp(i,k)-qn(i,k)),cwat(i,k))/deltat 
! Subprogram not used 
! Subprogram not used            ! 4. no condensation nor evaporation
! Subprogram not used            ! ==================================
! Subprogram not used            else
! Subprogram not used               cme(i,k)=0.0_r8
! Subprogram not used            endif
! Subprogram not used 
! Subprogram not used   
! Subprogram not used         end do    !end loop for cme update
! Subprogram not used 
! Subprogram not used ! Because of the finite time step, 
! Subprogram not used ! place a bound here not to exceed wet bulb point
! Subprogram not used ! and not to evaporate more than available water
! Subprogram not used !
! Subprogram not used          do i = 1, ncol
! Subprogram not used             qtmp = qn(i,k) - cme(i,k)*deltat
! Subprogram not used 
! Subprogram not used ! possibilities to have qtmp > qsp
! Subprogram not used !
! Subprogram not used !   1. if qn > qs(tn), it condenses; 
! Subprogram not used !      if after applying cme,  qtmp > qsp,  more condensation is applied. 
! Subprogram not used !      
! Subprogram not used !   2. if qn < qs, evaporation should not exceed qsp,
! Subprogram not used     
! Subprogram not used             if(qtmp > qsp(i,k)) then
! Subprogram not used               cme(i,k) = cme(i,k) + (qtmp-qsp(i,k))/deltat
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! if net evaporation, it should not exceed available cwat
! Subprogram not used !
! Subprogram not used             if(cme(i,k) < -cwat(i,k)/deltat)  &
! Subprogram not used                cme(i,k) = -cwat(i,k)/deltat
! Subprogram not used !
! Subprogram not used ! addition of residual condensation from previous step of iteration
! Subprogram not used !
! Subprogram not used             cme(i,k) = cme(i,k) + cmeres(i)
! Subprogram not used 
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          !      limit cme for roundoff errors
! Subprogram not used          do i = 1, ncol
! Subprogram not used             cme(i,k) = cme(i,k)*omsm
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          do i = 1,ncol
! Subprogram not used !
! Subprogram not used ! as a safe limit, condensation should not reduce grid mean rh below rhu00
! Subprogram not used ! 
! Subprogram not used            if(cme(i,k) > 0.0_r8 .and. relhum(i) > rhu00(i,k) )  &
! Subprogram not used               cme(i,k) = min(cme(i,k), (qn(i,k)-qs(i)*rhu00(i,k))/deltat)
! Subprogram not used !
! Subprogram not used ! initial guess for cwm (mean cloud water over time step) if 1st iteration
! Subprogram not used !
! Subprogram not used            if(l < 2) then
! Subprogram not used              cwm(i) = max(cwat(i,k)+cme(i,k)*dto2,  0._r8)
! Subprogram not used            endif
! Subprogram not used 
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used ! provisional precipitation falling through model layer
! Subprogram not used          do i = 1,ncol
! Subprogram not used !!$            prprov(i) =  precab(i) + prodprec(i,k)*pdel(i,k)/gravit
! Subprogram not used ! rain produced in this layer not too effective in collection process
! Subprogram not used             wtthick = max(0._r8,min(0.5_r8,((zi(i,k)-zi(i,k+1))/1000._r8)**2))
! Subprogram not used             prprov(i) =  precab(i) + wtthick*prodprec(i,k)*pdel(i,k)/gravit
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used ! calculate conversion of condensate to precipitation by cloud microphysics 
! Subprogram not used          call findmcnew (lchnk   ,ncol    , &
! Subprogram not used                          k       ,prprov  ,snowab,  t       ,p        , &
! Subprogram not used                          cwm     ,cldm    ,cldmax  ,fice(1,k),coef    , &
! Subprogram not used                          fwaut(1,k),fsaut(1,k),fracw(1,k),fsacw(1,k),fsaci(1,k), &
! Subprogram not used                          seaicef, snowh, pracwo(1,k), psacwo(1,k), psacio(1,k))
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! calculate the precip rate
! Subprogram not used !
! Subprogram not used          error_found = .false.
! Subprogram not used          do i = 1,ncol
! Subprogram not used             if (cldm(i) > 0) then  
! Subprogram not used !
! Subprogram not used ! first predict the cloud water
! Subprogram not used !
! Subprogram not used                cdt = coef(i)*deltat
! Subprogram not used                if (cdt > 0.01_r8) then
! Subprogram not used                   pol = cme(i,k)/coef(i) ! production over loss
! Subprogram not used                   cwn(i) = max(0._r8,(cwat(i,k)-pol)*exp(-cdt)+ pol)
! Subprogram not used                else
! Subprogram not used                   cwn(i) = max(0._r8,(cwat(i,k) + cme(i,k)*deltat)/(1+cdt))
! Subprogram not used                endif
! Subprogram not used !
! Subprogram not used ! now back out the tendency of net rain production
! Subprogram not used !
! Subprogram not used                prodprec(i,k) =  max(0._r8,cme(i,k)-(cwn(i)-cwat(i,k))/deltat)
! Subprogram not used             else
! Subprogram not used                prodprec(i,k) = 0.0_r8
! Subprogram not used                cwn(i) = 0._r8
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             ! provisional calculation of conversion terms
! Subprogram not used             ice2pr(i,k) = prodprec(i,k)*(fsaut(i,k)+fsaci(i,k))
! Subprogram not used             liq2pr(i,k) = prodprec(i,k)*(fwaut(i,k)+fsacw(i,k)+fracw(i,k))
! Subprogram not used !old        liq2snow(i,k) = prodprec(i,k)*fsacw(i,k)
! Subprogram not used 
! Subprogram not used !           revision suggested by Jim McCaa
! Subprogram not used !           it controls the amount of snow hitting the sfc 
! Subprogram not used !           by forcing a lot of conversion of cloud liquid to snow phase
! Subprogram not used !           it might be better done later by an explicit representation of 
! Subprogram not used !           rain accreting ice (and freezing), or by an explicit freezing of raindrops
! Subprogram not used             liq2snow(i,k) = max(prodprec(i,k)*fsacw(i,k), fsnow(i,k)*liq2pr(i,k))
! Subprogram not used 
! Subprogram not used             ! bounds
! Subprogram not used             nice2pr = min(ice2pr(i,k),(cwat(i,k)+cme(i,k)*deltat)*fice(i,k)/deltat)
! Subprogram not used             nliq2pr = min(liq2pr(i,k),(cwat(i,k)+cme(i,k)*deltat)*(1._r8-fice(i,k))/deltat)
! Subprogram not used !            write(iulog,*) ' prodprec ', i, k, prodprec(i,k)
! Subprogram not used !            write(iulog,*) ' nliq2pr, nice2pr ', nliq2pr, nice2pr
! Subprogram not used             if (liq2pr(i,k).ne.0._r8) then
! Subprogram not used                nliq2snow = liq2snow(i,k)*nliq2pr/liq2pr(i,k)   ! correction
! Subprogram not used             else
! Subprogram not used                nliq2snow = liq2snow(i,k)
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used !           avoid roundoff problems generating negatives
! Subprogram not used             nliq2snow = nliq2snow*omsm
! Subprogram not used             nliq2pr = nliq2pr*omsm
! Subprogram not used             nice2pr = nice2pr*omsm
! Subprogram not used             
! Subprogram not used !           final estimates of conversion to precip and snow
! Subprogram not used             prodprec(i,k) = (nliq2pr + nice2pr)
! Subprogram not used             prodsnow(i,k) = (nice2pr + nliq2snow)
! Subprogram not used 
! Subprogram not used             rcwn(i,l,k) =  cwat(i,k) + (cme(i,k)-   prodprec(i,k))*deltat
! Subprogram not used             rliq(i,l,k) = (cwat(i,k) + cme(i,k)*deltat)*(1._r8-fice(i,k)) - nliq2pr * deltat
! Subprogram not used             rice(i,l,k) = (cwat(i,k) + cme(i,k)*deltat)* fice(i,k)      -    nice2pr                     *deltat
! Subprogram not used 
! Subprogram not used !           Save for sanity check later...  
! Subprogram not used !           Putting sanity checks inside loops 100 and 800 screws up the 
! Subprogram not used !           IBM compiler for reasons as yet unknown.  TBH
! Subprogram not used             cwnsave(i,l,k)      = cwn(i)
! Subprogram not used             cmesave(i,l,k)      = cme(i,k)
! Subprogram not used             prodprecsave(i,l,k) = prodprec(i,k)
! Subprogram not used !           End of save for sanity check later...  
! Subprogram not used 
! Subprogram not used !           final version of condensate to precip terms
! Subprogram not used             liq2pr(i,k) = nliq2pr
! Subprogram not used             liq2snow(i,k) = nliq2snow
! Subprogram not used             ice2pr(i,k) = nice2pr
! Subprogram not used 
! Subprogram not used             cwn(i) = rcwn(i,l,k)
! Subprogram not used !
! Subprogram not used ! update any remaining  provisional values
! Subprogram not used !
! Subprogram not used             cwm(i) = (cwn(i) + cwat(i,k))*0.5_r8
! Subprogram not used !
! Subprogram not used ! update in cloud water
! Subprogram not used !
! Subprogram not used             if(cldm(i) > mincld) then
! Subprogram not used                icwc(i) = cwm(i)/cldm(i)
! Subprogram not used             else
! Subprogram not used                icwc(i) = 0.0_r8
! Subprogram not used             endif
! Subprogram not used !PJR the above logic give zero icwc with nonzero cwat, dont like it!
! Subprogram not used !PJR generates problems with csigma
! Subprogram not used !PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
! Subprogram not used !         icwc(i) = max(1.e-8_r8,cwm(i)/cldm(i))
! Subprogram not used 
! Subprogram not used          end do              ! end of do i = 1,ncol
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! calculate provisional value of cloud water for
! Subprogram not used ! evaporation of precipitate (evapprec) calculation
! Subprogram not used !
! Subprogram not used       do i = 1,ncol
! Subprogram not used          qtmp = qn(i,k) - cme(i,k)*deltat
! Subprogram not used          ttmp = tn(i,k) + deltat/cpair * ( meltheat(i,k)       &
! Subprogram not used               + (latvap + latice*fice(i,k)) * cme(i,k) )
! Subprogram not used          esn = estblf(ttmp)
! Subprogram not used          qsn = svp_to_qsat(esn, p(i,k))
! Subprogram not used          qtl(i) = max((qsn - qtmp)/deltat,0._r8)
! Subprogram not used          relhum1(i) = qtmp/qsn
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used       do i = 1,ncol
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used          evapprec(i,k) = conke*(1._r8 - cldm(i))*sqrt(precab(i)) &
! Subprogram not used                          *(1._r8 - min(relhum1(i),1._r8))
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! limit the evaporation to the amount which is entering the box
! Subprogram not used ! or saturates the box
! Subprogram not used !
! Subprogram not used          prtmp = precab(i)*gravit/pdel(i,k)
! Subprogram not used          evapprec(i,k) = min(evapprec(i,k), prtmp, qtl(i))*omsm
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Partition evaporation of precipitate between rain and snow using
! Subprogram not used ! the fraction of snow falling into the box. Determine the heating
! Subprogram not used ! due to evaporation. Note that evaporation is positive (loss of precip,
! Subprogram not used ! gain of vapor) and that heating is negative.
! Subprogram not used          if (evapprec(i,k) > 0._r8) then
! Subprogram not used             evapsnow(i,k) = evapprec(i,k) * snowab(i) / precab(i)
! Subprogram not used             evapheat(i,k) = -latvap * evapprec(i,k) - latice * evapsnow(i,k)
! Subprogram not used          else 
! Subprogram not used             evapsnow(i,k) = 0._r8
! Subprogram not used             evapheat(i,k) = 0._r8
! Subprogram not used          end if
! Subprogram not used ! Account for the latent heat of fusion for liquid drops collected by falling snow
! Subprogram not used          prfzheat(i,k) = latice * liq2snow(i,k)
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used ! now remove the residual of any over-saturation. Normally,
! Subprogram not used ! the oversaturated water vapor should have been removed by 
! Subprogram not used ! cme formulation plus constraints by wet bulb tsp/qsp
! Subprogram not used ! as computed above. However, because of non-linearity,
! Subprogram not used ! addition of (cme-evapprec) to update t and q may still cause
! Subprogram not used ! a very small amount of over saturation. It is called a
! Subprogram not used ! residual of over-saturation because theoretically, cme
! Subprogram not used ! should have taken care of all of large scale condensation.
! Subprogram not used ! 
! Subprogram not used 
! Subprogram not used        do i = 1,ncol
! Subprogram not used           qtmp = qn(i,k)-(cme(i,k)-evapprec(i,k))*deltat
! Subprogram not used           ttmp = tn(i,k) + deltat/cpair * ( meltheat(i,k) + evapheat(i,k) + prfzheat(i,k)      &
! Subprogram not used               + (latvap + latice*fice(i,k)) * cme(i,k) )
! Subprogram not used 
! Subprogram not used           call qsat(ttmp, p(i,k), esn, qsn, dqsdt=dqsdt)
! Subprogram not used 
! Subprogram not used           if( qtmp > qsn ) then
! Subprogram not used              !
! Subprogram not used              !now extra condensation to bring air to just saturation
! Subprogram not used              !
! Subprogram not used              ctmp = (qtmp-qsn)/(1._r8+hlocp*dqsdt)/deltat
! Subprogram not used              cme(i,k) = cme(i,k)+ctmp
! Subprogram not used !
! Subprogram not used ! save residual on cmeres to addtion to cme on entering next iteration
! Subprogram not used ! cme exit here contain the residual but overrided if back to iteration
! Subprogram not used !
! Subprogram not used              cmeres(i) = ctmp
! Subprogram not used           else
! Subprogram not used              cmeres(i) = 0.0_r8
! Subprogram not used           endif
! Subprogram not used        end do
! Subprogram not used               
! Subprogram not used  100 continue              ! end of do l = 1,iter
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! precipitation
! Subprogram not used !
! Subprogram not used       do i = 1,ncol
! Subprogram not used          precip(i) = precip(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
! Subprogram not used          precab(i) = precab(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
! Subprogram not used          if(precab(i).lt.0._r8) precab(i)=0._r8
! Subprogram not used !         snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodprec(i,k)*fice(i,k) - evapsnow(i,k))
! Subprogram not used          snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodsnow(i,k) - evapsnow(i,k))
! Subprogram not used 
! Subprogram not used          ! If temperature above freezing, all precip is rain flux.  if temperature below freezing, all precip is snow flux.
! Subprogram not used          rkflxprc(i,k+1) = precab(i)   !! making this consistent with other precip fluxes.  prc = rain + snow
! Subprogram not used          !!rkflxprc(i,k+1) = precab(i) - snowab(i)
! Subprogram not used          rkflxsnw(i,k+1) = snowab(i)
! Subprogram not used 
! Subprogram not used !!$         if ((precab(i)) < 1.e-10) then      
! Subprogram not used !!$            precab(i) = 0.
! Subprogram not used !!$            snowab(i) = 0.
! Subprogram not used !!$         endif
! Subprogram not used       end do
! Subprogram not used  800 continue                ! level loop (k=1,pver)
! Subprogram not used 
! Subprogram not used ! begin sanity checks
! Subprogram not used    error_found = .false.
! Subprogram not used    do k = top_lev,pver
! Subprogram not used       do l = 1,iter
! Subprogram not used          do i = 1,ncol
! Subprogram not used 	    if (abs(rcwn(i,l,k)).lt.1.e-300_r8) rcwn(i,l,k) = 0._r8
! Subprogram not used 	    if (abs(rliq(i,l,k)).lt.1.e-300_r8) rliq(i,l,k) = 0._r8
! Subprogram not used 	    if (abs(rice(i,l,k)).lt.1.e-300_r8) rice(i,l,k) = 0._r8
! Subprogram not used             if (rcwn(i,l,k).lt.0._r8) error_found = .true.
! Subprogram not used             if (rliq(i,l,k).lt.0._r8) error_found = .true.
! Subprogram not used             if (rice(i,l,k).lt.0._r8) error_found = .true.
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used    if (error_found) then
! Subprogram not used       do k = top_lev,pver
! Subprogram not used          do l = 1,iter
! Subprogram not used             do i = 1,ncol
! Subprogram not used                if (rcwn(i,l,k).lt.0._r8) then
! Subprogram not used                   write(iulog,*) ' prob with neg rcwn1 ', rcwn(i,l,k),  &
! Subprogram not used                      cwnsave(i,l,k)
! Subprogram not used                   write(iulog,*) ' cwat, cme*deltat, prodprec*deltat ', &
! Subprogram not used                      cwat(i,k), cmesave(i,l,k)*deltat,               &
! Subprogram not used                      prodprecsave(i,l,k)*deltat,                     &
! Subprogram not used                      (cmesave(i,l,k)-prodprecsave(i,l,k))*deltat
! Subprogram not used                   call endrun('PCOND')
! Subprogram not used                endif
! Subprogram not used                if (rliq(i,l,k).lt.0._r8) then
! Subprogram not used                   write(iulog,*) ' prob with neg rliq1 ', rliq(i,l,k)
! Subprogram not used                   call endrun('PCOND')
! Subprogram not used                endif
! Subprogram not used                if (rice(i,l,k).lt.0._r8) then
! Subprogram not used                   write(iulog,*) ' prob with neg rice ', rice(i,l,k)
! Subprogram not used                   call endrun('PCOND')
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used    end if
! Subprogram not used ! end sanity checks
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used end subroutine pcond

!##############################################################################

! Subprogram not used subroutine findmcnew (lchnk   ,ncol    , &
! Subprogram not used                       k       ,precab  ,snowab,  t       ,p       , &
! Subprogram not used                       cwm     ,cldm    ,cldmax  ,fice    ,coef    , &
! Subprogram not used                       fwaut   ,fsaut   ,fracw   ,fsacw   ,fsaci   , &
! Subprogram not used                       seaicef ,snowh,  pracwo, psacwo, psacio )
! Subprogram not used  
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! calculate the conversion of condensate to precipitate
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
! Subprogram not used !  model climate using diagnosed and 
! Subprogram not used !  predicted condensate parameterizations, 1998, J. Clim., 11,
! Subprogram not used !  pp1587---1614.
! Subprogram not used ! 
! Subprogram not used ! Author: P. Rasch
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use phys_grid, only: get_rlat_all_p
! Subprogram not used    use comsrf,        only: landm
! Subprogram not used !
! Subprogram not used ! input args
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: lchnk                 ! chunk identifier
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used    integer, intent(in) :: k                     ! level index
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: precab(pcols)        ! rate of precipitation from above (kg / (m**2 * s))
! Subprogram not used    real(r8), intent(in) :: t(pcols,pver)        ! temperature       (K)
! Subprogram not used    real(r8), intent(in) :: p(pcols,pver)        ! pressure          (Pa)
! Subprogram not used    real(r8), intent(in) :: cldm(pcols)          ! cloud fraction
! Subprogram not used    real(r8), intent(in) :: cldmax(pcols)        ! max cloud fraction above this level
! Subprogram not used    real(r8), intent(in) :: cwm(pcols)           ! condensate mixing ratio (kg/kg)
! Subprogram not used    real(r8), intent(in) :: fice(pcols)          ! fraction of cwat that is ice
! Subprogram not used    real(r8), intent(in) :: seaicef(pcols)       ! sea ice fraction 
! Subprogram not used    real(r8), intent(in) :: snowab(pcols)        ! rate of snow from above (kg / (m**2 * s))
! Subprogram not used     real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
! Subprogram not used 
! Subprogram not used ! output arguments
! Subprogram not used    real(r8), intent(out) :: coef(pcols)          ! conversion rate (1/s)
! Subprogram not used    real(r8), intent(out) :: fwaut(pcols)         ! relative importance of liquid autoconversion (a diagnostic)
! Subprogram not used    real(r8), intent(out) :: fsaut(pcols)         ! relative importance of ice autoconversion    (a diagnostic)
! Subprogram not used    real(r8), intent(out) :: fracw(pcols)         ! relative importance of rain accreting liquid (a diagnostic)
! Subprogram not used    real(r8), intent(out) :: fsacw(pcols)         ! relative importance of snow accreting liquid (a diagnostic)
! Subprogram not used    real(r8), intent(out) :: fsaci(pcols)         ! relative importance of snow accreting ice    (a diagnostic)
! Subprogram not used    real(r8), intent(out) :: pracwo(pcols)        ! accretion of cloud water by rain (1/s)
! Subprogram not used    real(r8), intent(out) :: psacwo(pcols)        ! accretion of cloud water by snow (1/s)
! Subprogram not used    real(r8), intent(out) :: psacio(pcols)        ! accretion of cloud ice by snow (1/s)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! work variables
! Subprogram not used 
! Subprogram not used    integer i
! Subprogram not used    integer ii
! Subprogram not used    integer ind(pcols)
! Subprogram not used    integer ncols
! Subprogram not used 
! Subprogram not used    real(r8), parameter :: degrad = 57.296_r8 ! divide by this to convert degrees to radians
! Subprogram not used    real(r8) capn                 ! local cloud particles / cm3
! Subprogram not used    real(r8) capnoice             ! local cloud particles when not over sea ice / cm3
! Subprogram not used    real(r8) ciaut                ! coefficient of autoconversion of ice (1/s)
! Subprogram not used    real(r8) cldloc(pcols)        ! non-zero amount of cloud
! Subprogram not used    real(r8) cldpr(pcols)         ! assumed cloudy volume occupied by rain and cloud
! Subprogram not used    real(r8) con1                 ! work constant
! Subprogram not used    real(r8) con2                 ! work constant
! Subprogram not used    real(r8) csacx                ! constant used for snow accreting liquid or ice
! Subprogram not used !!$   real(r8) dtice                ! interval for transition from liquid to ice
! Subprogram not used    real(r8) icemr(pcols)         ! in-cloud ice mixing ratio
! Subprogram not used    real(r8) icrit                ! threshold for autoconversion of ice
! Subprogram not used    real(r8) liqmr(pcols)         ! in-cloud liquid water mixing ratio
! Subprogram not used    real(r8) pracw                ! rate of rain accreting water
! Subprogram not used    real(r8) prlloc(pcols)        ! local rain flux in mm/day
! Subprogram not used    real(r8) prscgs(pcols)        ! local snow amount in cgs units
! Subprogram not used    real(r8) psaci                ! rate of collection of ice by snow (lin et al 1983)
! Subprogram not used    real(r8) psacw                ! rate of collection of liquid by snow (lin et al 1983)
! Subprogram not used    real(r8) psaut                ! rate of autoconversion of ice condensate
! Subprogram not used    real(r8) ptot                 ! total rate of conversion
! Subprogram not used    real(r8) pwaut                ! rate of autoconversion of liquid condensate
! Subprogram not used    real(r8) r3l                  ! volume radius
! Subprogram not used    real(r8) rainmr(pcols)        ! in-cloud rain mixing ratio
! Subprogram not used    real(r8) rat1                 ! work constant
! Subprogram not used    real(r8) rat2                 ! work constant
! Subprogram not used !!$   real(r8) rdtice               ! recipricol of dtice
! Subprogram not used    real(r8) rho(pcols)           ! density (mks units)
! Subprogram not used    real(r8) rhocgs               ! density (cgs units)
! Subprogram not used    real(r8) rlat(pcols)          ! latitude (radians)
! Subprogram not used    real(r8) snowfr               ! fraction of precipate existing as snow
! Subprogram not used    real(r8) totmr(pcols)         ! in-cloud total condensate mixing ratio
! Subprogram not used    real(r8) vfallw               ! fall speed of precipitate as liquid
! Subprogram not used    real(r8) wp                   ! weight factor used in calculating pressure dep of autoconversion
! Subprogram not used    real(r8) wsi                  ! weight factor for sea ice
! Subprogram not used    real(r8) wt                   ! fraction of ice
! Subprogram not used    real(r8) wland                ! fraction of land
! Subprogram not used 
! Subprogram not used !      real(r8) csaci
! Subprogram not used !      real(r8) csacw
! Subprogram not used !      real(r8) cwaut
! Subprogram not used !      real(r8) efact
! Subprogram not used !      real(r8) lamdas
! Subprogram not used !      real(r8) lcrit
! Subprogram not used !      real(r8) rcwm
! Subprogram not used !      real(r8) r3lc2
! Subprogram not used !      real(r8) snowmr(pcols)
! Subprogram not used !      real(r8) vfalls
! Subprogram not used 
! Subprogram not used    real(8) ftot
! Subprogram not used 
! Subprogram not used !     inline statement functions
! Subprogram not used    real(r8) heavy, heavym, a1, a2, heavyp, heavymp
! Subprogram not used    heavy(a1,a2) = max(0._r8,sign(1._r8,a1-a2))  ! heavyside function
! Subprogram not used    heavym(a1,a2) = max(0.01_r8,sign(1._r8,a1-a2))  ! modified heavyside function
! Subprogram not used !
! Subprogram not used ! New heavyside functions to perhaps address error growth problems
! Subprogram not used !
! Subprogram not used    heavyp(a1,a2) = a1/(a2+a1+1.e-36_r8)
! Subprogram not used    heavymp(a1,a2) = (a1+0.01_r8*a2)/(a2+a1+1.e-36_r8)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! find all the points where we need to do the microphysics
! Subprogram not used ! and set the output variables to zero
! Subprogram not used !
! Subprogram not used    ncols = 0
! Subprogram not used    do i = 1,ncol
! Subprogram not used       coef(i) = 0._r8
! Subprogram not used       fwaut(i) = 0._r8
! Subprogram not used       fsaut(i) = 0._r8
! Subprogram not used       fracw(i) = 0._r8
! Subprogram not used       fsacw(i) = 0._r8
! Subprogram not used       fsaci(i) = 0._r8
! Subprogram not used       liqmr(i) = 0._r8
! Subprogram not used       rainmr(i) = 0._r8
! Subprogram not used       if (cwm(i) > 1.e-20_r8) then
! Subprogram not used          ncols = ncols + 1
! Subprogram not used          ind(ncols) = i
! Subprogram not used       endif
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used    do ii = 1,ncols
! Subprogram not used       i = ind(ii)
! Subprogram not used !
! Subprogram not used ! the local cloudiness at this level
! Subprogram not used !
! Subprogram not used       cldloc(i) = max(cldmin,cldm(i))
! Subprogram not used !
! Subprogram not used ! a weighted mean between max cloudiness above, and this layer
! Subprogram not used !
! Subprogram not used       cldpr(i) = max(cldmin,(cldmax(i)+cldm(i))*0.5_r8)
! Subprogram not used !
! Subprogram not used ! decompose the suspended condensate into
! Subprogram not used ! an incloud liquid and ice phase component
! Subprogram not used !
! Subprogram not used       totmr(i) = cwm(i)/cldloc(i)
! Subprogram not used       icemr(i) = totmr(i)*fice(i)
! Subprogram not used       liqmr(i) = totmr(i)*(1._r8-fice(i))
! Subprogram not used !
! Subprogram not used ! density
! Subprogram not used !
! Subprogram not used       rho(i) = p(i,k)/(287._r8*t(i,k))
! Subprogram not used       rhocgs = rho(i)*1.e-3_r8     ! density in cgs units
! Subprogram not used !
! Subprogram not used ! decompose the precipitate into a liquid and ice phase
! Subprogram not used !
! Subprogram not used       if (t(i,k) > t0) then
! Subprogram not used          vfallw = convfw/sqrt(rho(i))
! Subprogram not used          rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
! Subprogram not used          snowfr = 0
! Subprogram not used !        snowmr(i)
! Subprogram not used       else
! Subprogram not used          snowfr = 1
! Subprogram not used          rainmr(i) = 0._r8
! Subprogram not used       endif
! Subprogram not used !     rainmr(i) = (precab(i)-snowab(i))/(rho(i)*vfallw*cldpr(i))
! Subprogram not used !
! Subprogram not used ! local snow amount in cgs units
! Subprogram not used !
! Subprogram not used       prscgs(i) = precab(i)/cldpr(i)*0.1_r8*snowfr
! Subprogram not used !     prscgs(i) = snowab(i)/cldpr(i)*0.1
! Subprogram not used !
! Subprogram not used ! local rain amount in mm/day
! Subprogram not used !
! Subprogram not used       prlloc(i) = precab(i)*86400._r8/cldpr(i)
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    con1 = 1._r8/(1.333_r8*pi)**0.333_r8 * 0.01_r8 ! meters
! Subprogram not used !
! Subprogram not used ! calculate the conversion terms
! Subprogram not used !
! Subprogram not used    call get_rlat_all_p(lchnk, ncol, rlat)
! Subprogram not used 
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used    do ii = 1,ncols
! Subprogram not used       i = ind(ii)
! Subprogram not used       rhocgs = rho(i)*1.e-3_r8     ! density in cgs units
! Subprogram not used !
! Subprogram not used ! exponential temperature factor
! Subprogram not used !
! Subprogram not used !        efact = exp(0.025*(t(i,k)-t0))
! Subprogram not used !
! Subprogram not used ! some temperature dependent constants
! Subprogram not used !
! Subprogram not used !!$      wt = min(1._r8,max(0._r8,(t0-t(i,k))*rdtice))
! Subprogram not used       wt = fice(i)
! Subprogram not used       icrit = icritc*wt + icritw*(1-wt)
! Subprogram not used !
! Subprogram not used ! jrm Reworked droplet number concentration algorithm
! Subprogram not used       ! Start with pressure-dependent value appropriate for continental air
! Subprogram not used       ! Note: reltab has a temperature dependence here
! Subprogram not used       capn = capnw + (capnc-capnw) * min(1._r8,max(0._r8,1.0_r8-(p(i,k)-0.8_r8*p(i,pver))/(0.2_r8*p(i,pver))))
! Subprogram not used       ! Modify for snow depth over land
! Subprogram not used       capn = capn + (capnc-capn) * min(1.0_r8,max(0.0_r8,snowh(i)*10._r8))
! Subprogram not used       ! Ramp between polluted value over land to clean value over ocean.
! Subprogram not used       capn = capn + (capnc-capn) * min(1.0_r8,max(0.0_r8,1.0_r8-landm(i,lchnk)))
! Subprogram not used       ! Ramp between the resultant value and a sea ice value in the presence of ice.
! Subprogram not used       capn = capn + (capnsi-capn) * min(1.0_r8,max(0.0_r8,seaicef(i)))
! Subprogram not used ! end jrm
! Subprogram not used !      
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! useful terms in following calculations
! Subprogram not used !
! Subprogram not used       rat1 = rhocgs/rhow
! Subprogram not used       rat2 = liqmr(i)/capn
! Subprogram not used       con2 = (rat1*rat2)**0.333_r8
! Subprogram not used !
! Subprogram not used ! volume radius
! Subprogram not used !
! Subprogram not used !        r3l = (rhocgs*liqmr(i)/(1.333*pi*capn*rhow))**0.333 * 0.01 ! meters
! Subprogram not used       r3l = con1*con2
! Subprogram not used !
! Subprogram not used ! critical threshold for autoconversion if modified for mixed phase
! Subprogram not used ! clouds to mimic a bergeron findeisen process
! Subprogram not used ! r3lc2 = r3lcrit*(1.-0.5*fice(i)*(1-fice(i)))
! Subprogram not used !
! Subprogram not used ! autoconversion of liquid
! Subprogram not used !
! Subprogram not used !        cwaut = 2.e-4
! Subprogram not used !        cwaut = 1.e-3
! Subprogram not used !        lcrit = 2.e-4
! Subprogram not used !        lcrit = 5.e-4
! Subprogram not used !        pwaut = max(0._r8,liqmr(i)-lcrit)*cwaut
! Subprogram not used !
! Subprogram not used ! pwaut is following tripoli and cotton (and many others)
! Subprogram not used ! we reduce the autoconversion below critpr, because these are regions where
! Subprogram not used ! the drop size distribution is likely to imply much smaller collector drops than
! Subprogram not used ! those relevant for a cloud distribution corresponding to the value of effc = 0.55
! Subprogram not used ! suggested by cotton (see austin 1995 JAS, baker 1993)
! Subprogram not used 
! Subprogram not used ! easy to follow form
! Subprogram not used !        pwaut = capc*liqmr(i)**2*rhocgs/rhow
! Subprogram not used !    $           *(liqmr(i)*rhocgs/(rhow*capn))**(.333)
! Subprogram not used !    $           *heavy(r3l,r3lcrit)
! Subprogram not used !    $           *max(0.10_r8,min(1._r8,prlloc(i)/critpr))
! Subprogram not used ! somewhat faster form
! Subprogram not used !#ifdef PERGRO
! Subprogram not used       pwaut = capc*liqmr(i)**2*rat1*con2*heavymp(r3l,r3lcrit) * &
! Subprogram not used               max(0.10_r8,min(1._r8,prlloc(i)/critpr))
! Subprogram not used !
! Subprogram not used ! autoconversion of ice
! Subprogram not used !
! Subprogram not used !        ciaut = ciautb*efact
! Subprogram not used       ciaut = ciautb
! Subprogram not used !        psaut = capc*totmr(i)**2*rhocgs/rhoi
! Subprogram not used !     $           *(totmr(i)*rhocgs/(rhoi*capn))**(.333)
! Subprogram not used !
! Subprogram not used ! autoconversion of ice condensate
! Subprogram not used !
! Subprogram not used       psaut = max(0._r8,icemr(i)-icrit)*ciaut
! Subprogram not used !
! Subprogram not used ! collection of liquid by rain
! Subprogram not used !
! Subprogram not used !        pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
! Subprogram not used       pracw = cracw*rho(i)*sqrt(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)
! Subprogram not used 
! Subprogram not used       pracwo(i)=pracw
! Subprogram not used 
! Subprogram not used !!      pracw = 0.
! Subprogram not used !
! Subprogram not used ! the following lines calculate the slope parameter and snow mixing ratio
! Subprogram not used ! from the precip rate using the equations found in lin et al 83
! Subprogram not used ! in the most natural form, but it is expensive, so after some tedious
! Subprogram not used ! algebraic manipulation you can use the cheaper form found below
! Subprogram not used !            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
! Subprogram not used !     $               *0.01   ! convert from cm/s to m/s
! Subprogram not used !            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
! Subprogram not used !            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
! Subprogram not used !            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25
! Subprogram not used !            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd)
! Subprogram not used !
! Subprogram not used ! coefficient for collection by snow independent of phase
! Subprogram not used !
! Subprogram not used       csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! collection of liquid by snow (lin et al 1983)
! Subprogram not used !
! Subprogram not used       psacw = csacx*liqmr(i)*esw
! Subprogram not used 
! Subprogram not used       psacwo(i)=psacw
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! collection of ice by snow (lin et al 1983)
! Subprogram not used !
! Subprogram not used       psaci = csacx*icemr(i)*esi
! Subprogram not used !
! Subprogram not used       psacio(i)=psaci
! Subprogram not used 
! Subprogram not used ! total conversion of condensate to precipitate
! Subprogram not used !
! Subprogram not used       ptot = pwaut + psaut + pracw + psacw + psaci
! Subprogram not used !
! Subprogram not used ! the recipricol of cloud water amnt (or zero if no cloud water)
! Subprogram not used !
! Subprogram not used !         rcwm =  totmr(i)/(max(totmr(i),small)**2)
! Subprogram not used !
! Subprogram not used ! turn the tendency back into a loss rate (1/seconds)
! Subprogram not used !
! Subprogram not used       if (totmr(i) > 0._r8) then
! Subprogram not used          coef(i) = ptot/totmr(i)
! Subprogram not used       else
! Subprogram not used          coef(i) = 0._r8
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (ptot.gt.0._r8) then
! Subprogram not used          fwaut(i) = pwaut/ptot
! Subprogram not used          fsaut(i) = psaut/ptot
! Subprogram not used          fracw(i) = pracw/ptot
! Subprogram not used          fsacw(i) = psacw/ptot
! Subprogram not used          fsaci(i) = psaci/ptot
! Subprogram not used       else
! Subprogram not used          fwaut(i) = 0._r8
! Subprogram not used          fsaut(i) = 0._r8
! Subprogram not used          fracw(i) = 0._r8
! Subprogram not used          fsacw(i) = 0._r8
! Subprogram not used          fsaci(i) = 0._r8
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       ftot = fwaut(i)+fsaut(i)+fracw(i)+fsacw(i)+fsaci(i)
! Subprogram not used !      if (abs(ftot-1._r8).gt.1.e-14_r8.and.ftot.ne.0._r8) then
! Subprogram not used !         write(iulog,*) ' something is wrong in findmcnew ', ftot, &
! Subprogram not used !              fwaut(i),fsaut(i),fracw(i),fsacw(i),fsaci(i)
! Subprogram not used !         write(iulog,*) ' unscaled ', ptot, &
! Subprogram not used !              pwaut,psaut,pracw,psacw,psaci
! Subprogram not used !         write(iulog,*) ' totmr, liqmr, icemr ', totmr(i), liqmr(i), icemr(i)
! Subprogram not used !         call endrun()
! Subprogram not used !      endif
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used end subroutine findmcnew

end module cldwat
