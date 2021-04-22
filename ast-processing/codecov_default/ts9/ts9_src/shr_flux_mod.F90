!===============================================================================
! SVN $Id: shr_flux_mod.F90 26627 2011-02-01 18:09:17Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140626/shr/shr_flux_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: flux_mod -- CCSM shared flux calculations.
!
! !DESCRIPTION:
!
!     CCSM shared flux calculations.
!     
! !REVISION HISTORY:
!     2006-Nov-07 - B. Kauffman - first version, code taken/migrated from cpl6
!
! !INTERFACE: ------------------------------------------------------------------

module shr_flux_mod

! !USES:

   use shr_kind_mod    ! shared kinds
   use shr_const_mod   ! shared constants
   use shr_sys_mod     ! shared system routines
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none

   private ! default private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_flux_atmOcn      ! computes atm/ocn fluxes
   public :: shr_flux_atmIce      ! computes atm/ice fluxes
   public :: shr_flux_MOstability ! boundary layer stability scales/functions

! !PUBLIC DATA MEMBERS:

  integer(SHR_KIND_IN),parameter,public :: shr_flux_MOwScales   = 1 ! w scales  option
  integer(SHR_KIND_IN),parameter,public :: shr_flux_MOfunctions = 2 ! functions option
  real   (SHR_KIND_R8),parameter,public :: shr_flux_MOgammaM = 3.59_SHR_KIND_R8
  real   (SHR_KIND_R8),parameter,public :: shr_flux_MOgammaS = 7.86_SHR_KIND_R8

!EOP

   !--- rename kinds for local readability only ---
   integer,parameter :: R8 = SHR_KIND_R8  ! 8 byte real
   integer,parameter :: IN = SHR_KIND_IN  ! native/default integer

   integer,parameter :: debug = 0 ! internal debug level

!===============================================================================
contains
!===============================================================================
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_atmOcn -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!     
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - code migrated from cpl5 to cpl6
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!     2006-Nov-07 - B. Kauffman - code migrated from cpl6 to share
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_flux_atmOcn(nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   & 
           &               qbot  ,rbot  ,tbot  ,us    ,vs    ,   &
           &               ts    ,mask  ,sen   ,lat   ,lwup  ,   &
           &               evap  ,taux  ,tauy  ,tref  ,qref  ,   &
           &               duu10n,  ustar_sv   ,re_sv ,ssq_sv,   &
           &               missval    )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in) ::       nMax  ! data vector length
   integer(IN),intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
   real(R8)   ,intent(in) :: zbot (nMax) ! atm level height      (m)
   real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind            (m/s)
   real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind            (m/s)
   real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T       (K)
   real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (nMax) ! atm T                 (K) 
   real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature       (K)

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
   real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
   real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
   real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
   real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
   real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

   real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
   real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
   real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)
 
   real(R8),intent(in) ,optional :: missval        ! masked value

! !EOP

   !--- local constants --------------------------------
   real(R8),parameter :: umin  =  0.5_R8 ! minimum wind speed       (m/s)
   real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
   real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

   !--- local variables --------------------------------
   integer(IN) :: n      ! vector loop index
   real(R8)    :: vmag   ! surface wind magnitude   (m/s)
   real(R8)    :: thvbot ! virtual temperature      (K)
   real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
   real(R8)    :: delt   ! potential T difference   (K)
   real(R8)    :: delq   ! humidity difference      (kg/kg)
   real(R8)    :: stable ! stability factor
   real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum) 
   real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)     
   real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)    
   real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)         
   real(R8)    :: rh     ! sqrt of exchange coefficient (heat)             
   real(R8)    :: re     ! sqrt of exchange coefficient (water)            
   real(R8)    :: ustar  ! ustar             
   real(R8)    :: qstar  ! qstar             
   real(R8)    :: tstar  ! tstar             
   real(R8)    :: hol    ! H (at zbot) over L
   real(R8)    :: xsq    ! ?
   real(R8)    :: xqq    ! ?
   real(R8)    :: psimh  ! stability function at zbot (momentum)
   real(R8)    :: psixh  ! stability function at zbot (heat and water)
   real(R8)    :: psix2  ! stability function at ztref reference height
   real(R8)    :: alz    ! ln(zbot/zref)
   real(R8)    :: al2    ! ln(zref/ztref)
   real(R8)    :: u10n   ! 10m neutral wind 
   real(R8)    :: tau    ! stress at zbot
   real(R8)    :: cp     ! specific heat of moist air
   real(R8)    :: bn     ! exchange coef funct for interpolation
   real(R8)    :: bh     ! exchange coef funct for interpolation
   real(R8)    :: fac    ! vertical interpolation factor
   real(R8)    :: spval  ! local missing value

   !--- local functions --------------------------------
   real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
   real(R8)    :: cdn    ! function: neutral drag coeff at 10m
   real(R8)    :: psimhu ! function: unstable part of psimh
   real(R8)    :: psixhu ! function: unstable part of psimx
   real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
   real(R8)    :: Tk     ! dummy arg ~ temperature (K)
   real(R8)    :: xd     ! dummy arg ~ ?
 
   qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
   cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
   psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
   psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)
 
   !--- formats ----------------------------------------
   character(*),parameter :: subName = '(shr_flux_atmOcn) '
   character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES: 
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + sen + lat
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
! 
! ASSUMPTIONS:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-------------------------------------------------------------------------------

   if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) "enter"

   if (present(missval)) then
      spval = missval
   else
      spval = shr_const_spval
   endif
 
   al2 = log(zref/ztref)

   DO n=1,nMax
     if (mask(n) /= 0) then
    
        !--- compute some needed quantities ---
        vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
        thvbot = thbot(n) * (1.0_R8 + shr_const_zvir * qbot(n)) ! virtual temp (K)
        ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)
        delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
        delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
        alz    = log(zbot(n)/zref) 
        cp     = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*ssq) 
   
        !------------------------------------------------------------
        ! first estimate of Z/L and ustar, tstar and qstar
        !------------------------------------------------------------
   
        !--- neutral coefficients, z/L = 0.0 ---
        stable = 0.5_R8 + sign(0.5_R8 , delt)
        rdn    = sqrt(cdn(vmag))
        rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8 
        ren    = 0.0346_R8 
   
        !--- ustar, tstar, qstar ---
        ustar = rdn * vmag
        tstar = rhn * delt  
        qstar = ren * delq  
   
        !--- compute stability & evaluate all stability functions ---
        hol  = shr_const_karman*shr_const_g*zbot(n)*  &
               (tstar/thbot(n)+qstar/(1.0_R8/shr_const_zvir+qbot(n)))/ustar**2
        hol  = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
   
        !--- shift wind speed using old coefficient ---
        rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
        u10n = vmag * rd / rdn 
   
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346_R8
        rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8 
    
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh)) 
        rh = rhn / (1.0_R8 + rhn/shr_const_karman*(alz-psixh)) 
        re = ren / (1.0_R8 + ren/shr_const_karman*(alz-psixh)) 
   
        !--- update ustar, tstar, qstar using updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! iterate to converge on Z/L, ustar, tstar and qstar
        !------------------------------------------------------------
    
        !--- compute stability & evaluate all stability functions ---
        hol  = shr_const_karman*shr_const_g*zbot(n)* &
               (tstar/thbot(n)+qstar/(1.0_R8/shr_const_zvir+qbot(n)))/ustar**2
        hol  = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
    
        !--- shift wind speed using old coeffs ---
        rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
        u10n = vmag * rd/rdn 
    
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346_R8
        rhn = (1.0_R8 - stable)*0.0327_R8 + stable * 0.018_R8 
   
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh)) 
        rh = rhn / (1.0_R8 + rhn/shr_const_karman*(alz-psixh)) 
        re = ren / (1.0_R8 + ren/shr_const_karman*(alz-psixh)) 
    
        !--- update ustar, tstar, qstar using updated, shifted coeffs ---
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! compute the fluxes
        !------------------------------------------------------------
    
        tau = rbot(n) * ustar * ustar 
       
        !--- momentum flux ---
        taux(n) = tau * (ubot(n)-us(n)) / vmag 
        tauy(n) = tau * (vbot(n)-vs(n)) / vmag 
        
        !--- heat flux ---
        sen (n) =                cp * tau * tstar / ustar 
        lat (n) =  shr_const_latvap * tau * qstar / ustar
        lwup(n) = -shr_const_stebol * ts(n)**4 
      
        !--- water flux ---
        evap(n) = lat(n)/shr_const_latvap 
    
        !------------------------------------------------------------
        ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
        !------------------------------------------------------------
        hol = hol*ztref/zbot(n)
        xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
        xqq = sqrt(xsq)
        psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
        fac     = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
        tref(n) = thbot(n) - delt*fac 
        tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
        fac     = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
        qref(n) =  qbot(n) - delq*fac
    
        duu10n(n) = u10n*u10n ! 10m wind speed squared

        !------------------------------------------------------------
        ! optional diagnostics, needed for water tracer fluxes (dcn)
        !------------------------------------------------------------
        if (present(ustar_sv)) ustar_sv(n) = ustar
        if (present(re_sv   )) re_sv(n)    = re
        if (present(ssq_sv  )) ssq_sv(n)   = ssq

     else
        !------------------------------------------------------------
        ! no valid data here -- out of domain
        !------------------------------------------------------------
        sen   (n) = spval  ! sensible         heat flux  (W/m^2)
        lat   (n) = spval  ! latent           heat flux  (W/m^2)
        lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
        evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
        taux  (n) = spval  ! x surface stress (N)
        tauy  (n) = spval  ! y surface stress (N)
        tref  (n) = spval  !  2m reference height temperature (K)
        qref  (n) = spval  !  2m reference height humidity (kg/kg)
        duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

        if (present(ustar_sv)) ustar_sv(n) = spval
        if (present(re_sv   )) re_sv   (n) = spval
        if (present(ssq_sv  )) ssq_sv  (n) = spval
     endif
   ENDDO 

END subroutine shr_flux_atmOcn

!BOP ===========================================================================
!
! !IROUTINE: shr_flux_atmIce -- computes atm/ice fluxes
!
! !DESCRIPTION:
!    Computes atm/ice fluxes
!
! !REVISION HISTORY:
!    2006-Jun-12 - B. Kauffman, first version, adapted from dice6 code
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine shr_flux_atmIce(mask  ,zbot  ,ubot  ,vbot  ,thbot  &
! Subprogram not used                &          ,qbot  ,rbot  ,tbot  ,ts    ,sen    &
! Subprogram not used                &          ,lat   ,lwup  ,evap  ,taux  ,tauy   &
! Subprogram not used                &          ,tref  ,qref                        )
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    !--- input arguments --------------------------------
! Subprogram not used    integer(IN),intent(in)  :: mask (:)    ! 0 <=> cell NOT in model domain
! Subprogram not used    real(R8)   ,intent(in)  :: zbot (:)    ! atm level height  (m)
! Subprogram not used    real(R8)   ,intent(in)  :: ubot (:)    ! atm u wind     (m/s)
! Subprogram not used    real(R8)   ,intent(in)  :: vbot (:)    ! atm v wind     (m/s)
! Subprogram not used    real(R8)   ,intent(in)  :: thbot(:)    ! atm potential T   (K)
! Subprogram not used    real(R8)   ,intent(in)  :: qbot (:)    ! atm specific humidity (kg/kg)
! Subprogram not used    real(R8)   ,intent(in)  :: rbot (:)    ! atm air density   (kg/m^3)
! Subprogram not used    real(R8)   ,intent(in)  :: tbot (:)    ! atm T       (K) 
! Subprogram not used    real(R8)   ,intent(in)  :: ts   (:)    ! surface temperature
! Subprogram not used 
! Subprogram not used    !--- output arguments -------------------------------
! Subprogram not used    real(R8)   ,intent(out) :: sen  (:)    ! sensible      heat flux  (W/m^2)
! Subprogram not used    real(R8)   ,intent(out) :: lat  (:)    ! latent        heat flux  (W/m^2)
! Subprogram not used    real(R8)   ,intent(out) :: lwup (:)    ! long-wave upward heat flux  (W/m^2)
! Subprogram not used    real(R8)   ,intent(out) :: evap (:)    ! evaporative water flux ((kg/s)/m^2)
! Subprogram not used    real(R8)   ,intent(out) :: taux (:)    ! x surface stress (N)
! Subprogram not used    real(R8)   ,intent(out) :: tauy (:)    ! y surface stress (N)
! Subprogram not used    real(R8)   ,intent(out) :: tref (:)    ! 2m reference height temperature
! Subprogram not used    real(R8)   ,intent(out) :: qref (:)    ! 2m reference height humidity
! Subprogram not used  
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local constants --------------------------------
! Subprogram not used    real(R8),parameter :: umin   =  1.0_R8            ! minimum wind speed (m/s)
! Subprogram not used    real(R8),parameter :: zref   = 10.0_R8            ! ref height           ~ m
! Subprogram not used    real(R8),parameter :: ztref  =  2.0_R8            ! ref height for air T ~ m
! Subprogram not used    real(R8),parameter :: spval  = shr_const_spval    ! special value
! Subprogram not used    real(R8),parameter :: g      = shr_const_g        ! gravity
! Subprogram not used    real(R8),parameter :: cpdair = shr_const_cpdair   ! spec heat of dry air
! Subprogram not used    real(R8),parameter :: cpvir  = shr_const_cpvir    ! cpwv/cpdair - 1.0
! Subprogram not used    real(R8),parameter :: zvir   = shr_const_zvir     ! rh2o/rair   - 1.0
! Subprogram not used    real(R8),parameter :: latvap = shr_const_latvap   ! latent heat of evap
! Subprogram not used    real(R8),parameter :: latice = shr_const_latice   ! latent heat of fusion
! Subprogram not used    real(R8),parameter :: stebol = shr_const_stebol   ! Stefan-Boltzmann
! Subprogram not used    real(R8),parameter :: karman = shr_const_karman   ! Von Karman constant
! Subprogram not used    real(R8),parameter :: zzsice = 0.0005_R8          ! ice surface roughness
! Subprogram not used 
! Subprogram not used    !--- local variables --------------------------------
! Subprogram not used    integer(IN) :: lsize  ! array dimensions
! Subprogram not used    integer(IN) :: n      ! array indicies
! Subprogram not used    real(R8)    :: vmag   ! surface wind magnitude   (m/s)
! Subprogram not used    real(R8)    :: thvbot ! virtual temperature      (K)
! Subprogram not used    real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
! Subprogram not used    real(R8)    :: dssqdt ! derivative of ssq wrt Ts (kg/kg/K)
! Subprogram not used    real(R8)    :: delt   ! potential T difference   (K)
! Subprogram not used    real(R8)    :: delq   ! humidity difference      (kg/kg)
! Subprogram not used    real(R8)    :: stable ! stability factor
! Subprogram not used    real(R8)    :: rdn    ! sqrt of neutral exchange coefficient (momentum)
! Subprogram not used    real(R8)    :: rhn    ! sqrt of neutral exchange coefficient (heat)
! Subprogram not used    real(R8)    :: ren    ! sqrt of neutral exchange coefficient (water)
! Subprogram not used    real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
! Subprogram not used    real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
! Subprogram not used    real(R8)    :: re     ! sqrt of exchange coefficient (water)      
! Subprogram not used    real(R8)    :: ustar  ! ustar
! Subprogram not used    real(R8)    :: qstar  ! qstar
! Subprogram not used    real(R8)    :: tstar  ! tstar
! Subprogram not used    real(R8)    :: hol    ! H (at zbot) over L
! Subprogram not used    real(R8)    :: xsq    ! temporary variable
! Subprogram not used    real(R8)    :: xqq    ! temporary variable
! Subprogram not used    real(R8)    :: psimh  ! stability function at zbot (momentum)
! Subprogram not used    real(R8)    :: psixh  ! stability function at zbot (heat and water)
! Subprogram not used    real(R8)    :: alz    ! ln(zbot/z10)
! Subprogram not used    real(R8)    :: ltheat ! latent heat for surface
! Subprogram not used    real(R8)    :: tau    ! stress at zbot
! Subprogram not used    real(R8)    :: cp     ! specific heat of moist air
! Subprogram not used 
! Subprogram not used    real(R8)    :: bn     ! exchange coef funct for interpolation
! Subprogram not used    real(R8)    :: bh     ! exchange coef funct for interpolation
! Subprogram not used    real(R8)    :: fac    ! interpolation factor
! Subprogram not used    real(R8)    :: ln0    ! log factor for interpolation
! Subprogram not used    real(R8)    :: ln3    ! log factor for interpolation
! Subprogram not used 
! Subprogram not used    !--- local functions --------------------------------
! Subprogram not used    real(R8)   :: Tk      ! temperature (K)
! Subprogram not used    real(R8)   :: qsat    ! the saturation humididty of air (kg/m^3)
! Subprogram not used    real(R8)   :: dqsatdt ! derivivative of qsat wrt surface temperature
! Subprogram not used    real(R8)   :: xd      ! dummy argument  
! Subprogram not used    real(R8)   :: psimhu  ! unstable part of psimh
! Subprogram not used    real(R8)   :: psixhu  ! unstable part of psimx
! Subprogram not used 
! Subprogram not used    qsat(Tk)    = 627572.4_R8 / exp(5107.4_R8/Tk)
! Subprogram not used    dqsatdt(Tk) = (5107.4_R8 / Tk**2) * 627572.4_R8 / exp(5107.4_R8/Tk)
! Subprogram not used    psimhu(xd)  = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
! Subprogram not used    psixhu(xd)  =  2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)
! Subprogram not used 
! Subprogram not used    !--- formats ----------------------------------------
! Subprogram not used    character(*),parameter :: subName =  "(shr_flux_atmIce) "
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE:
! Subprogram not used !   using atm & ice state variables, compute atm/ice fluxes 
! Subprogram not used !   and diagnostic 10m air temperature and humidity
! Subprogram not used !
! Subprogram not used ! NOTE: 
! Subprogram not used !   o all fluxes are positive downward
! Subprogram not used !   o net heat flux = net sw + lw up + lw down + sen + lat
! Subprogram not used !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
! Subprogram not used !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
! Subprogram not used ! 
! Subprogram not used ! ASSUME:
! Subprogram not used !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    lsize = size(tbot)
! Subprogram not used 
! Subprogram not used    do n = 1,lsize
! Subprogram not used 
! Subprogram not used      if (mask(n) == 0) then
! Subprogram not used         sen  (n) = spval
! Subprogram not used         lat  (n) = spval
! Subprogram not used         lwup (n) = spval
! Subprogram not used         evap (n) = spval
! Subprogram not used         taux (n) = spval
! Subprogram not used         tauy (n) = spval
! Subprogram not used         tref (n) = spval
! Subprogram not used         qref (n) = spval
! Subprogram not used      else
! Subprogram not used         !--- define some needed variables ---
! Subprogram not used         vmag   = max(umin, sqrt(ubot(n)**2+vbot(n)**2))
! Subprogram not used         thvbot = thbot(n)*(1.0_R8 + zvir * qbot(n)) ! virtual pot temp (K)
! Subprogram not used          ssq   =  qsat  (ts(n)) / rbot(n)           ! sea surf hum (kg/kg)
! Subprogram not used         dssqdt = dqsatdt(ts(n)) / rbot(n)           ! deriv of ssq wrt Ts 
! Subprogram not used         delt   = thbot(n) - ts(n)                   ! pot temp diff (K)
! Subprogram not used         delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)
! Subprogram not used         alz    = log(zbot(n)/zref) 
! Subprogram not used         cp     = cpdair*(1.0_R8 + cpvir*ssq) 
! Subprogram not used         ltheat = latvap + latice
! Subprogram not used 
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used         ! first estimate of Z/L and ustar, tstar and qstar
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used 
! Subprogram not used         !--- neutral coefficients, z/L = 0.0 ---
! Subprogram not used         rdn = karman/log(zref/zzsice)
! Subprogram not used         rhn = rdn
! Subprogram not used         ren = rdn
! Subprogram not used 
! Subprogram not used         !--- ustar,tstar,qstar ----
! Subprogram not used         ustar = rdn * vmag
! Subprogram not used         tstar = rhn * delt  
! Subprogram not used         qstar = ren * delq  
! Subprogram not used 
! Subprogram not used         !--- compute stability & evaluate all stability functions ---
! Subprogram not used         hol    = karman * g * zbot(n) &
! Subprogram not used         &     * (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n))) / ustar**2
! Subprogram not used         hol    = sign( min(abs(hol),10.0_R8), hol )
! Subprogram not used         stable = 0.5_R8 + sign(0.5_R8 , hol)
! Subprogram not used         xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
! Subprogram not used         xqq    = sqrt(xsq)
! Subprogram not used         psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
! Subprogram not used         psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
! Subprogram not used 
! Subprogram not used         !--- shift all coeffs to measurement height and stability ---
! Subprogram not used         rd = rdn / (1.0_R8+rdn/karman*(alz-psimh))
! Subprogram not used         rh = rhn / (1.0_R8+rhn/karman*(alz-psixh))
! Subprogram not used         re = ren / (1.0_R8+ren/karman*(alz-psixh))
! Subprogram not used 
! Subprogram not used         !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
! Subprogram not used         ustar = rd * vmag 
! Subprogram not used         tstar = rh * delt 
! Subprogram not used         qstar = re * delq 
! Subprogram not used 
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used         ! iterate to converge on Z/L, ustar, tstar and qstar
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used 
! Subprogram not used         !--- compute stability & evaluate all stability functions ---
! Subprogram not used         hol    = karman * g * zbot(n) &
! Subprogram not used         &      * (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n))) / ustar**2
! Subprogram not used         hol    = sign( min(abs(hol),10.0_R8), hol )
! Subprogram not used         stable = 0.5_R8 + sign(0.5_R8 , hol)
! Subprogram not used         xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
! Subprogram not used         xqq    = sqrt(xsq)
! Subprogram not used         psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
! Subprogram not used         psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
! Subprogram not used 
! Subprogram not used         !--- shift all coeffs to measurement height and stability ---
! Subprogram not used         rd = rdn / (1.0_R8+rdn/karman*(alz-psimh)) 
! Subprogram not used         rh = rhn / (1.0_R8+rhn/karman*(alz-psixh)) 
! Subprogram not used         re = ren / (1.0_R8+ren/karman*(alz-psixh)) 
! Subprogram not used 
! Subprogram not used         !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
! Subprogram not used         ustar = rd * vmag 
! Subprogram not used         tstar = rh * delt 
! Subprogram not used         qstar = re * delq 
! Subprogram not used 
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used         ! compute the fluxes
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used 
! Subprogram not used         tau = rbot(n) * ustar * ustar 
! Subprogram not used     
! Subprogram not used         !--- momentum flux ---
! Subprogram not used         taux(n) = tau * ubot(n) / vmag 
! Subprogram not used         tauy(n) = tau * vbot(n) / vmag 
! Subprogram not used      
! Subprogram not used         !--- heat flux ---
! Subprogram not used         sen (n) =   cp * tau * tstar / ustar 
! Subprogram not used         lat (n) =  ltheat * tau * qstar / ustar
! Subprogram not used         lwup(n) = -stebol * ts(n)**4 
! Subprogram not used      
! Subprogram not used         !--- water flux ---
! Subprogram not used         evap(n) = lat(n)/ltheat 
! Subprogram not used 
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used         ! compute diagnostic: 2m reference height temperature
! Subprogram not used         !----------------------------------------------------------
! Subprogram not used 
! Subprogram not used         !--- Compute function of exchange coefficients. Assume that 
! Subprogram not used         !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
! Subprogram not used         !--- 1/sqrt(cn(n))=1/rdn and sqrt(cm(n))/ch(n)=1/rh 
! Subprogram not used         bn = karman/rdn
! Subprogram not used         bh = karman/rh
! Subprogram not used 
! Subprogram not used         !--- Interpolation factor for stable and unstable cases
! Subprogram not used         ln0 = log(1.0_R8 + (ztref/zbot(n))*(exp(bn) - 1.0_R8))
! Subprogram not used         ln3 = log(1.0_R8 + (ztref/zbot(n))*(exp(bn - bh) - 1.0_R8))
! Subprogram not used         fac = (ln0 - ztref/zbot(n)*(bn - bh))/bh * stable &
! Subprogram not used         &   + (ln0 - ln3)/bh * (1.0_R8-stable)
! Subprogram not used         fac = min(max(fac,0.0_R8),1.0_R8)
! Subprogram not used 
! Subprogram not used         !--- actual interpolation
! Subprogram not used         tref(n) = ts(n) + (tbot(n) - ts(n))*fac
! Subprogram not used         qref(n) = qbot(n) - delq*fac
! Subprogram not used 
! Subprogram not used      endif
! Subprogram not used    enddo 
! Subprogram not used 
! Subprogram not used end subroutine shr_flux_atmIce
 
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_MOstability -- Monin-Obukhov BL stability functions
!
! !DESCRIPTION:
!
!    Monin-Obukhov boundary layer stability functions, two options:
!    turbulent velocity scales or gradient and integral functions
!    via option = shr_flux_MOwScales or shr_flux_MOfunctions
!     
! !REVISION HISTORY:
!    2007-Sep-19 - B. Kauffman, Bill Large - first version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine shr_flux_MOstability(option,arg1,arg2,arg3,arg4,arg5)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer(IN),intent(in)           :: option ! shr_flux_MOwScales or MOfunctions 
! Subprogram not used    real(R8)   ,intent(in)           :: arg1   ! scales: uStar (in)  funct: zeta (in)
! Subprogram not used    real(R8)   ,intent(inout)        :: arg2   ! scales: zkB   (in)  funct: phim (out)
! Subprogram not used    real(R8)   ,intent(out)          :: arg3   ! scales: phim  (out) funct: phis (out)
! Subprogram not used    real(R8)   ,intent(out)          :: arg4   ! scales: phis  (out) funct: psim (out)
! Subprogram not used    real(R8)   ,intent(out),optional :: arg5   ! scales:    (unused) funct: psis (out)
! Subprogram not used  
! Subprogram not used ! !EOP
! Subprogram not used 
! Subprogram not used    !----- local variables -----
! Subprogram not used    real(R8)           :: zeta  ! z/L
! Subprogram not used    real(R8)           :: uStar ! friction velocity
! Subprogram not used    real(R8)           :: zkB   ! (height)*(von Karman)*(surface bouyancy flux)
! Subprogram not used    real(R8)           :: phim  ! momentum    gradient function or scale
! Subprogram not used    real(R8)           :: phis  ! temperature gradient function or scale
! Subprogram not used    real(R8)           :: psim  ! momentum    integral function or scale
! Subprogram not used    real(R8)           :: psis  ! temperature integral function or scale
! Subprogram not used    real(R8)           :: temp  ! temporary-variable/partial calculation
! Subprogram not used 
! Subprogram not used    !----- local variables, stable case -----
! Subprogram not used    real(R8),parameter :: uStarMin = 0.001_R8 ! lower bound on uStar
! Subprogram not used    real(R8),parameter :: a = 1.000_R8  ! constant from Holtslag & de Bruin, equation 12
! Subprogram not used    real(R8),parameter :: b = 0.667_R8  ! constant from Holtslag & de Bruin, equation 12
! Subprogram not used    real(R8),parameter :: c = 5.000_R8  ! constant from Holtslag & de Bruin, equation 12
! Subprogram not used    real(R8),parameter :: d = 0.350_R8  ! constant from Holtslag & de Bruin, equation 12
! Subprogram not used 
! Subprogram not used    !----- local variables, unstable case -----
! Subprogram not used    real(R8),parameter :: a2 = 3.0_R8   ! constant from Wilson, equation 10 
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(shr_flux_MOstability) '
! Subprogram not used    character(*),parameter ::   F00 = "('(shr_flux_MOstability) ',4a)"
! Subprogram not used    character(*),parameter ::   F01 = "('(shr_flux_MOstability) ',a,i5)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes::
! Subprogram not used !   o this could be two routines, but are one to help keep them aligned
! Subprogram not used !   o the stable calculation is taken from...
! Subprogram not used !     A.A.M. HoltSlag and H.A.R. de Bruin, 1988:
! Subprogram not used !     "Applied Modeling of the Nighttime Surface Energy Balance over Land",
! Subprogram not used !     Journal of Applied Meteorology, Vol. 27, No. 6, June 1988, 659-704
! Subprogram not used !   o the unstable calculation is taken from...
! Subprogram not used !     D. Keith Wilson, 2001: "An Alternative Function for the Wind and 
! Subprogram not used !     Temperature Gradients in Unstable Surface Layers", 
! Subprogram not used !     Boundary-Layer Meteorology, 99 (2001), 151-158
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    !----- check for consistancy between option and arguments ------------------
! Subprogram not used    if (debug > 1 .and. s_loglev > 0) then
! Subprogram not used       if (debug > 2) write(s_logunit,F01) "enter, option = ",option
! Subprogram not used       if ( option == shr_flux_MOwScales .and. present(arg5) ) then
! Subprogram not used          write(s_logunit,F01) "ERROR: option1 must have four arguments"
! Subprogram not used          call shr_sys_abort(subName//"option inconsistant with arguments")
! Subprogram not used       else if ( option == shr_flux_MOfunctions .and. .not. present(arg5) ) then
! Subprogram not used          write(s_logunit,F01) "ERROR: option2 must have five arguments"
! Subprogram not used          call shr_sys_abort(subName//"option inconsistant with arguments")
! Subprogram not used       else 
! Subprogram not used          write(s_logunit,F01) "invalid option = ",option
! Subprogram not used          call shr_sys_abort(subName//"invalid option")
! Subprogram not used       end if
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    !------ velocity scales option ----------------------------------------------
! Subprogram not used    if (option == shr_flux_MOwScales) then
! Subprogram not used 
! Subprogram not used       !--- input ---
! Subprogram not used       uStar = arg1
! Subprogram not used       zkB   = arg2
! Subprogram not used 
! Subprogram not used       if (zkB >= 0.0_R8) then ! ----- stable -----
! Subprogram not used          zeta = zkB/(max(uStar,uStarMin)**3)
! Subprogram not used          temp = exp(-d*zeta)
! Subprogram not used          phim = uStar/(1.0_R8 + zeta*(a + b*(1.0_R8 + c - d*zeta)*temp))
! Subprogram not used          phis = phim
! Subprogram not used       else                    ! ----- unstable -----
! Subprogram not used          temp = (zkB*zkB)**(1.0_R8/a2)   ! note: zkB < 0, zkB*zkB > 0
! Subprogram not used          phim = sqrt(uStar**2 + shr_flux_MOgammaM*temp)
! Subprogram not used          phis = sqrt(uStar**2 + shr_flux_MOgammaS*temp)
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       !--- output ---
! Subprogram not used       arg3 = phim
! Subprogram not used       arg4 = phis
! Subprogram not used    !  arg5 = <unused>
! Subprogram not used 
! Subprogram not used    !------ stability function option -------------------------------------------
! Subprogram not used    else if (option == shr_flux_MOfunctions) then
! Subprogram not used 
! Subprogram not used       !--- input ---
! Subprogram not used       zeta  = arg1
! Subprogram not used 
! Subprogram not used       if (zeta >= 0.0_R8) then ! ----- stable -----
! Subprogram not used          temp = exp(-d*zeta)
! Subprogram not used          phim =        1.0_R8 + zeta*(a + b*(1.0_R8 + c - d*zeta)*temp)
! Subprogram not used          phis = phim
! Subprogram not used          psim = -a*zeta - b*(zeta - c/d)*temp - b*c/d
! Subprogram not used          psis = psim 
! Subprogram not used       else                    ! ----- unstable ----
! Subprogram not used          temp = (zeta*zeta)**(1.0_R8/a2)   ! note: zeta < 0, zeta*zeta > 0
! Subprogram not used          phim = 1.0_R8/sqrt(1.0_R8 + shr_flux_MOgammaM*temp)
! Subprogram not used          phis = 1.0_R8/sqrt(1.0_R8 + shr_flux_MOgammaS*temp)
! Subprogram not used          psim = a2*log(0.5_R8 + 0.5_R8/phim)
! Subprogram not used          psis = a2*log(0.5_R8 + 0.5_R8/phis)
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       !--- output ---
! Subprogram not used       arg2 = phim
! Subprogram not used       arg3 = phis
! Subprogram not used       arg4 = psim
! Subprogram not used       arg5 = psis
! Subprogram not used    !----------------------------------------------------------------------------
! Subprogram not used    else 
! Subprogram not used       write(s_logunit,F01) "invalid option = ",option
! Subprogram not used       call shr_sys_abort(subName//"invalid option")
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used end subroutine shr_flux_MOstability

!===============================================================================
!===============================================================================

end module shr_flux_mod
