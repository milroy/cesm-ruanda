      module mo_aurora
!-----------------------------------------------------------------------
!
! Auroral oval parameterization. See reference:
! R.G. Roble, E.C. Ridley
! An auroral model for the NCAR thermospheric general circulation model (TGCM)
! Annales Geophysicae,5A, (6), 369-382, 1987. 
!
! The aurora oval is a circle in auroral circle coordinates.  Auroral circle
!  coordinates are offset from magnetic coordinates by offa degrees (radians)
!  towards 0 MLT and by dskofa degrees (radians) towards dusk (18 MLT).
! The aurora assumes a Maxwellian in energy, so that the characteristic
!  energy is half of the mean energy (or mean energy = 2*alfa, where alfa
!  is the characteristic energy).  The Maxwellian is approximated in the
!  aion and bion subroutines.
! The aurora oval is assumed to be a Gaussian in auroral latitude, with
!  peak values on the day (=1) and night (=2) sides that change from one to
!  the other using cosines of the auroral longitude coordinate.
! There is provision for a low energy (~75 eV) aurora at the location of the
!  regular (~1-6 keV) aurora in order to simulate the energy flux found
!  at higher altitudes that is non-Maxwellian, but the flux is usually
!  set to zero (1.e-80).
! There is provision for a proton (MeV) aurora, but the flux is usually
!  set to zero (1.e-20).
! The drizzle is a constant low energy electron flux over the polar cap,
!  which goes to 1/e over twice the half-width of the aurora at the
!  radius of the aurora.
! The cusp is a low energy electron flux centered over the dayside convection
!  entrance at phid at the convection reversal boundary theta0.  The cusp
!  falls off over 5 degrees in latitude and over 20 degrees in longitude
!  to 1/e values of the peak at the center.
! 1.e-20 and 1.e-80 are used to give a near zero answer.
!
! The polar drizzle and cusp electron energies are low, and soft particles
!  have great influence on the over-all thermospheric and ionospheric
!  structure, especially on the electron density profiles at mid-latitudes
!  and in winter since low energy electrons produce ionization at high
!  altitudes where loss rates are very low.  (Comment by Wenbin Wang.)
! The original energies for drizzle and cusp were alfad=0.75, alfac=0.5 keV.
! The original guess at energy fluxes were: ed=0.1+2.0*power/100.,ec=0.1+0.9*power/100.
! The next guess at energy fluxes were: ed=0.01+0.2*power/100., ec=0.01+0.09*power/100.
! The values below reflect higher estimates for the electron energy (lower alt)
!
! Calling sequence (all subs in mo_aurora, mo_aurora.F):
!   1) sub aurora_cons called once per time step from advance.
!   2) sub aurora called from dynamics, inside parallel latitude scan.
!   3) subs aurora_cusp and aurora_heat called from sub aurora.
!   4) sub aurora_ions called from sub aurora. 
!
!-----------------------------------------------------------------------

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use mo_constants,  only: pi, &
                               avo => avogadro, &
                               boltz_cgs, &
                               gask => rgas_cgs
      use cam_logfile,   only: iulog
      use spmd_utils,    only: masterproc

      implicit none

      interface aurora
         module procedure aurora_prod
         module procedure aurora_hrate
      end interface

      save

      integer, parameter  :: isouth = 1
      integer, parameter  :: inorth = 2

      ! g = 8.7 m/s^2? Because this is 400 km up?
      real(r8), parameter :: grav   = 870._r8          ! (cm/s^2)

      integer  :: lev1 = 1
      real(r8) :: twopi
      real(r8) :: rmass_o1
      real(r8) :: rmass_o2
      real(r8) :: rmass_n2
      real(r8) :: rmassinv_o1
      real(r8) :: rmassinv_o2
      real(r8) :: rmassinv_n2
      real(r8) :: dtr

!-----------------------------------------------------------------------
! 	... polar drizzle parameters:
!   alfad: Characteristic Maxwellian energy of drizzle electrons (keV)
!   ed   : Column energy input of drizzle electrons (ergs/cm**2/s)
!   fd   : Electron particle flux of drizzle electrons (particles/cm**2/s)
!-----------------------------------------------------------------------
      real(r8), parameter :: alfad = 2.0_r8, &
                             ed    = 0.5_r8     
      real(r8) :: fd                     ! set in sub aurora_ions

!-----------------------------------------------------------------------
! 	... polar cusp parameters:
!   alfac: Characteristic Maxwellian energy of polar cusp electons (keV)
!   ec   : Column energy input of polar cusp electrons (ergs/cm**2/s)
!   fc   : Electron particle flux of polar cusp electrons (particles/cm**2/s)
!-----------------------------------------------------------------------
      real(r8), parameter :: alfac = 1.0_r8, &
                             ec    = 0.5_r8
      real(r8) :: fc                     ! set in sub aurora_ions

!-----------------------------------------------------------------------
! e1: Peak energy flux in noon sector of the aurora (ergs/cm**2/s)
! e2: Peak energy flux in midnight sector of the aurora (ergs/cm**2/s)
! h1: Gaussian half-width of the noon auroral oval in degrees
! h2: Gaussian half-width of the midnight auroral oval in degrees
!-----------------------------------------------------------------------
      real(r8) :: &
        e1, e2, &                        ! set in sub aurora_cons (function of hem power)
        h1, h2                           ! set in sub aurora_cons (function of hem power)

!-----------------------------------------------------------------------
! 	... solar proton parameters for the polar cap (time-gcm only)
!   alfa_sp: Characteristic Maxwellian energy of solar protons (MeV) (was alfad2)
!   e_sp   : Column energy input of solar protons (ergs/cm**2/s) (was ed2)
!   flx_sp : e_sp/1.602e-6, for input to sub bion                (was fd2)
! Add solar protons to ionization if add_sproton is true (time-gcm only)
!-----------------------------------------------------------------------
      logical :: add_sproton = .false.
      real(r8), parameter :: &
        alfa_sp = 10._r8, &
        e_sp    = 1.e-20_r8
      real(r8) :: flx_sp

!-----------------------------------------------------------------------
! 	... high energy electron parameters in the auroral oval (time-gcm only):
!-----------------------------------------------------------------------
      logical :: add_helectron = .false.
      real(r8), parameter :: alfa30 = 40._r8, &                  ! Characteristic energy of auroral electrons (Mev)
                             e30    = .05_r8                     ! Column energy of auroral electrons (ergs/cm**2/s)

!-----------------------------------------------------------------------
! 	... additional auroral parameters
!-----------------------------------------------------------------------
      real(r8) :: &
        alfa0, &        ! average of noon and midnight characteristic Maxw energies
        ralfa,ralfa2, & ! difference ratios of characteristic energies
        rrote, &        ! clockwise rotation from noon of peak dayside energy flux (e1)
        rroth, &        ! clockwise rotation from noon of dayside h1 Gaussian half-width
        h0, &           ! average of noon and midnight Gaussian half-widths
        rh, &           ! difference ratio of half-widths (rh=(h2-h1)/(h2+h1))
        e0,e20, &       ! e0 = average of noon and midnight electrons
        ree,re2, &      ! difference ratios of peak energy fluxes (ree=(e2-e1)/(e2+e1))
        alfa20          ! average of noon and midnight char energies for high alt aurora
      real(r8) :: &
        theta0(2), &    ! convection reversal boundary in radians
        offa(2), &      ! offset of oval towards 0 MLT relative to magnetic pole (rad)
        dskofa(2), &    ! offset of oval in radians towards 18 MLT (f(By))
        phid(2), &      ! dayside convection entrance in MLT converted to radians (f(By))
        rrad(2)         ! radius of auroral circle in radians
      real(r8) :: ctpoten    ! cross-cap potential (kV)
      real(r8) :: byimf      ! BY component of IMF (nT)


      private
      public :: aurora_inti, aurora_timestep_init, aurora
      public :: aurora_register

      logical :: has_ions = .false.
      integer :: indxAIPRS = -1

      contains

        
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
! Subprogram not used       subroutine aurora_register
! Subprogram not used         use ppgrid,       only : pver,pcols
! Subprogram not used         use physics_buffer, only : pbuf_add_field, dtype_r8
! Subprogram not used 
! Subprogram not used         ! add ionization rates to phys buffer for waccmx ionosphere module
! Subprogram not used 
! Subprogram not used         call pbuf_add_field('AurIPRateSum' , 'physpkg', dtype_r8, (/pcols,pver/), indxAIPRS)     ! Sum of ion auroral production rates for O2
! Subprogram not used 
! Subprogram not used       endsubroutine aurora_register

      subroutine aurora_inti
!-----------------------------------------------------------------------
! 	... initialize aurora module
!-----------------------------------------------------------------------

      use ppgrid,       only : pver
      use pmgrid,       only : plev, plevp
      use constituents, only : cnst_get_ind, cnst_mw
      use chem_mods,    only : adv_mass
      use ref_pres,     only : pref_mid
      use mo_chem_utls, only : get_spc_ndx
      use cam_history,  only : addfld, phys_decomp

      implicit none

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer             :: k, m
      real(r8), parameter :: e       = 1.e-10_r8
      real(r8), parameter :: convert = 3.1211e8_r8

      real(r8) :: plb
      real(r8) :: alfa_1, alfa_2, alfa21, alfa22
      real(r8) :: e21, e22

      integer :: op_ndx,o2p_ndx,np_ndx,n2p_ndx,e_ndx

      op_ndx   = get_spc_ndx( 'Op' )
      o2p_ndx  = get_spc_ndx( 'O2p' )
      np_ndx   = get_spc_ndx( 'Np' )
      n2p_ndx  = get_spc_ndx( 'N2p' )
      e_ndx    = get_spc_ndx( 'e' )

      has_ions = op_ndx > 0 .and. o2p_ndx > 0 .and. np_ndx > 0 .and. n2p_ndx > 0 .and. e_ndx > 0

      if (.not. has_ions) return

!-----------------------------------------------------------------------
!	... initialize module variables
!-----------------------------------------------------------------------
      twopi  = 2._r8*pi
      dtr    = pi/180._r8

!-----------------------------------------------------------------------
!	... set molecular weights
!-----------------------------------------------------------------------
      call cnst_get_ind( 'O2', m )
      rmass_o2    = cnst_mw(m)
      rmassinv_o2 = 1._r8/rmass_o2
      call cnst_get_ind( 'O', m )
      rmass_o1    = cnst_mw(m)
      rmassinv_o1 = 1._r8/rmass_o1
      call cnst_get_ind( 'N', m )
      rmass_n2    = 2._r8*cnst_mw(m)
      rmassinv_n2 = 1._r8/rmass_n2

      offa(isouth)   = 4.3_r8*dtr
      offa(inorth)   = 3.7_r8*dtr
      phid(isouth)   = 0._r8
      phid(inorth)   = 0._r8
      alfa_1         = 2._r8
      alfa_2         = 2._r8

!-----------------------------------------------------------------------
! Values from 10/05/94 HPI estimates (50% or more higher than old estimates):
!     alfa_1 = amin1(1.5,1.25+0.05*plevel)
!     alfa_2 = 1.2 + 0.095*plevel
!-----------------------------------------------------------------------
      alfa0  = 0.5_r8*(alfa_1 + alfa_2)
      ralfa  = (alfa_2 - alfa_1) / (alfa_1 + alfa_2 + e)
      alfa21 = 0.075_r8
      alfa22 = 0.075_r8
      alfa20 = 0.5_r8 * (alfa21 + alfa22)
      ralfa2 = (alfa22 - alfa21) / (alfa21 + alfa22 + e)
      e21    = 1.e-80_r8
      e22    = 1.e-80_r8
      e20    = 0.5_r8 * (e21 + e22)
      re2    = (e22 - e21) / (e21 + e22)

!-----------------------------------------------------------------------
! Set cusp and drizzle parameters:
! (conversion between particle number density and characteristic
!  energy and column energy input)
!-----------------------------------------------------------------------
      fc = convert * ec / alfac
      fd = convert * ed / alfad

!-----------------------------------------------------------------------
! Solar proton flux:
!-----------------------------------------------------------------------
      flx_sp = e_sp/1.602e-6_r8

!-----------------------------------------------------------------------
! 	... set auroral lower bndy index
!-----------------------------------------------------------------------
      plb = 5.e-4_r8*exp( 7._r8 ) * .1_r8             ! Pa
      do k = 1,pver
	 if( pref_mid(k) >= plb ) then
	    lev1 = k-1
	    exit
	 end if
      end do

      if (masterproc) write(iulog,*) ' '
      if (masterproc) write(iulog,*) 'aurora_inti: aurora will go down to lev,p = ',lev1,pref_mid(lev1)
      if (masterproc) write(iulog,*) ' '

!-----------------------------------------------------------------------
! Report to stdout:
!-----------------------------------------------------------------------
        call addfld('ALATM   ','RADIANS ',1,'I',&
             'Magnetic latitude at each geographic coordinate',phys_decomp)
        call addfld('ALONM   ','RADIANS ',1,'I',&
             'Magnetic longitude at each geographic coordinate',phys_decomp)
        call addfld( 'QSUM', '/s ', pver, 'I', &
             'total ion production', phys_decomp )

      end subroutine aurora_inti

! Subprogram not used       subroutine aurora_timestep_init
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... per timestep initialization
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use mag_parms,  only : get_mag_parms
! Subprogram not used       use spmd_utils, only : masterproc
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8) :: power, plevel
! Subprogram not used       real(r8) :: roth, rote, rcp, rhp
! Subprogram not used       real(r8) :: arad
! Subprogram not used 
! Subprogram not used       if (.not. has_ions) return
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !	... get hemispheric power
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call get_mag_parms( by = byimf, hpower = power, ctpoten = ctpoten )
! Subprogram not used 
! Subprogram not used       if( power >= .01_r8 ) then
! Subprogram not used          plevel = 2.09_r8*log( power )
! Subprogram not used       else
! Subprogram not used          plevel = 0._r8
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used ! fvitt -- moved the calc of h1, h2, rh, and h0 from aurora_inti
! Subprogram not used !          This was done to for bit-for-bit restarts.
! Subprogram not used !          These aurora oval dimension quantities power dependent and should be updated.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! h1 = Gaussian half-width of the noon auroral oval in degrees
! Subprogram not used ! h2 = Gaussian half-width of the midnight auroral oval in degrees
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !      h1 = 3._r8
! Subprogram not used !      h2 = 10._r8
! Subprogram not used ! modified by LQIAN, 2007
! Subprogram not used ! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! Subprogram not used ! h1 formula given by Wenbin base on POLARVIS image;
! Subprogram not used ! h2 formula based on Emery et al original auroral parameterization report
! Subprogram not used       h1 = min(2.35_r8, 0.83_r8 + 0.33_r8*plevel)
! Subprogram not used       h2 = 2.5_r8+0.025_r8*max(power,55._r8)+0.01_r8*min(0._r8,power-55._r8)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! Values from corrections to Emery et al Parameterization report:
! Subprogram not used !     h1 = amin1(2.35, 0.83 + 0.33*plevel)
! Subprogram not used !     h2 = 2.87 + 0.15*plevel
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       rh = (h2 - h1) / (h1 + h2)
! Subprogram not used       h0     = 0.5_r8 * (h1 + h2) * dtr
! Subprogram not used 
! Subprogram not used       theta0(isouth) = (-3.80_r8 + 8.48_r8*(ctpoten**.1875_r8))*dtr
! Subprogram not used       theta0(inorth) = theta0(isouth)
! Subprogram not used       dskofa(isouth) = (-1.26_r8 + 0.15_r8 * byimf)*dtr
! Subprogram not used       dskofa(inorth) = dskofa(isouth)
! Subprogram not used       roth           = (12.18_r8 - 0.89_r8 * plevel)
! Subprogram not used       rote           = ( 2.62_r8 - 0.55_r8 * plevel)
! Subprogram not used       rroth          = roth * dtr
! Subprogram not used       rrote          = rote * dtr
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! e1 = energy flux in the noon sector of the aurora (ergs/cm**2/s)
! Subprogram not used ! e2 = energy flux in the midnight sector of the aurora (ergs/cm**2/s)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !      e1 = (0.5_r8 + 0.15_r8 * power)
! Subprogram not used !      e2 = (1.5_r8 + 0.25_r8 * power)
! Subprogram not used ! modified by LQIAN, 2008
! Subprogram not used ! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! Subprogram not used ! e1 formula given by Wenbin base on POLARVIS image;
! Subprogram not used ! e2 formula based on Emery et al original auroral parameterization report
! Subprogram not used       e1 = max(0.50_r8, -2.15_r8 + 0.62_r8*plevel)
! Subprogram not used       e2=1._r8+0.11_r8*power
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! Values from corrections to Emery et al Parameterization report:
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       e0  = 0.5_r8 * (e1 + e2)
! Subprogram not used       ree = (e2 - e1) / (e1 + e2)
! Subprogram not used 
! Subprogram not used       rhp          = 14.20_r8 + 0.96_r8*plevel
! Subprogram not used       rcp          = -0.43_r8 + 9.69_r8 * (ctpoten**.1875_r8)
! Subprogram not used       arad         = max( rcp,rhp )
! Subprogram not used       rrad(isouth) = arad*dtr 
! Subprogram not used       rrad(inorth) = arad*dtr
! Subprogram not used 
! Subprogram not used       end subroutine aurora_timestep_init

! Subprogram not used       subroutine aurora_prod( tn, o2, o1, mbar, rlats, &
! Subprogram not used                               qo2p, qop, qn2p, qnp, pmid, &
! Subprogram not used                               lchnk, calday,  ncol, rlons, pbuf )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral parameterization driver
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use mo_apex,     only : alatm, alonm                      ! magnetic latitude,longitude grid (radians)
! Subprogram not used       use ppgrid,      only : pcols, pver
! Subprogram not used       use cam_history, only : outfld
! Subprogram not used       use physics_buffer,only: physics_buffer_desc
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in) ::  &
! Subprogram not used        ncol, &                           ! column count
! Subprogram not used        lchnk                             ! chunk index
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         calday                           ! calendar day of year
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         tn(pcols,pver), &                ! neutral gas temperature (K)
! Subprogram not used         o2(ncol,pver), &                 ! O2 concentration (kg/kg)
! Subprogram not used         o1(ncol,pver), &                 ! O concentration (kg/kg)
! Subprogram not used         mbar(ncol,pver)                  ! mean molecular weight (g/mole)
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         pmid(pcols,pver)                 ! midpoint pressure (Pa)
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         rlats(ncol), &                   ! column latitudes (radians)
! Subprogram not used         rlons(ncol)
! Subprogram not used       real(r8), intent(out) :: &
! Subprogram not used         qo2p(ncol,pver), &               ! o2+ production
! Subprogram not used         qop(ncol,pver), &                ! o+ production
! Subprogram not used         qn2p(ncol,pver), &               ! n2+ production
! Subprogram not used         qnp(ncol,pver)                   ! n+ production
! Subprogram not used 
! Subprogram not used       type(physics_buffer_desc),pointer :: pbuf(:)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer  :: i, k
! Subprogram not used       integer  :: hemis(ncol)
! Subprogram not used       real(r8) :: r2d
! Subprogram not used       real(r8) :: ofda, cosofa, sinofa, aslona
! Subprogram not used       real(r8) :: sunlons(ncol)                     ! sun's mag longitudes
! Subprogram not used       real(r8) :: dlat_aur(ncol)
! Subprogram not used       real(r8) :: dlon_aur(ncol)
! Subprogram not used       real(r8) :: colat(ncol)
! Subprogram not used       real(r8) :: sinlat(ncol)
! Subprogram not used       real(r8) :: coslat(ncol)
! Subprogram not used       real(r8) :: coslon(ncol)
! Subprogram not used       real(r8) :: sinlon(ncol)
! Subprogram not used       real(r8) :: alon(ncol)
! Subprogram not used       real(r8) :: cusp(ncol)
! Subprogram not used       real(r8) :: alfa(ncol)
! Subprogram not used       real(r8) :: alfa2(ncol)
! Subprogram not used       real(r8) :: alfa3(ncol)
! Subprogram not used       real(r8) :: flux(ncol)
! Subprogram not used       real(r8) :: flux2(ncol)
! Subprogram not used       real(r8) :: flux3(ncol)
! Subprogram not used       real(r8) :: drizl(ncol)
! Subprogram not used       real(r8) :: qteaur(ncol)                         ! for electron temperature
! Subprogram not used       logical  :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used       if (.not. has_ions) return
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... initialize ion production
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       do k = 1,pver
! Subprogram not used         qo2p(:,k) = 0._r8
! Subprogram not used         qop(:,k)  = 0._r8
! Subprogram not used         qn2p(:,k) = 0._r8
! Subprogram not used         qnp(:,k)  = 0._r8
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       r2d = 180._r8/pi
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... output mag lons, lats
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call outfld( 'ALONM', r2d*alonm(:ncol,lchnk), pcols, lchnk )
! Subprogram not used       call outfld( 'ALATM', r2d*alatm(:ncol,lchnk), pcols, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... check latitudes, and return if all below 32.5 deg
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       do_aurora(:) = abs( rlats(:) ) > pi/6._r8
! Subprogram not used       if( all( .not. do_aurora(:) ) ) then
! Subprogram not used          return
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... set sun location
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call sunloc( calday, sunlons, lchnk, ncol )
! Subprogram not used 
! Subprogram not used       do i = 1,ncol
! Subprogram not used         if( do_aurora(i) ) then
! Subprogram not used           dlat_aur(i) = alatm(i,lchnk)
! Subprogram not used           dlon_aur(i) = alonm(i,lchnk) - sunlons(i)
! Subprogram not used           if( dlon_aur(i) > pi ) then
! Subprogram not used              dlon_aur(i) = dlon_aur(i) - twopi
! Subprogram not used           else if( dlon_aur(i) < -pi ) then
! Subprogram not used              dlon_aur(i) = dlon_aur(i) + twopi
! Subprogram not used           end if
! Subprogram not used           if( dlat_aur(i) > 0._r8 ) then
! Subprogram not used             hemis(i) = 2
! Subprogram not used           else
! Subprogram not used             hemis(i) = 1
! Subprogram not used           end if
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... find auroral circle coordinates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used             ofda      = sqrt( offa(hemis(i))**2 + dskofa(hemis(i))**2)
! Subprogram not used             cosofa    = cos( ofda )
! Subprogram not used             sinofa    = sin( ofda )
! Subprogram not used             aslona    = asin( dskofa(hemis(i))/ofda )
! Subprogram not used             sinlat(i) = sin( abs( dlat_aur(i) ) )
! Subprogram not used             coslat(i) = cos( dlat_aur(i) )
! Subprogram not used             sinlon(i) = sin( dlon_aur(i) + aslona )
! Subprogram not used             coslon(i) = cos( dlon_aur(i) + aslona )
! Subprogram not used             colat(i)  = acos( cosofa*sinlat(i) - sinofa*coslat(i)*coslon(i))
! Subprogram not used             alon(i)   = mod( atan2( sinlon(i)*coslat(i),sinlat(i)*sinofa &
! Subprogram not used                                     + cosofa*coslat(i)*coslon(i) ) - aslona + 3._r8*pi,twopi) - pi
! Subprogram not used         end if
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... make cusp
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... make alfa, flux, and drizzle
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call aurora_heat( flux, flux2, flux3, alfa, alfa2, &
! Subprogram not used                         alfa3, qteaur, drizl, do_aurora, hemis, &
! Subprogram not used                         alon, colat, ncol )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral additions to ionization rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call aurora_ions( drizl, cusp, alfa, alfa2, alfa3, &
! Subprogram not used                         flux, flux2, flux3, tn, o2, &
! Subprogram not used                         o1, mbar, qo2p, qop, qn2p, &
! Subprogram not used                         qnp, pmid, do_aurora, ncol, lchnk, pbuf )
! Subprogram not used 
! Subprogram not used       end subroutine aurora_prod

! Subprogram not used       subroutine aurora_hrate( tn, o2, o1, mbar, rlats, &
! Subprogram not used                                aur_hrate, cpair, pmid, lchnk, calday, &
! Subprogram not used                                ncol, rlons )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral parameterization driver
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use mo_apex, only : alatm, alonm                      ! magnetic latitude,longitude grid (radians)
! Subprogram not used       use ppgrid,  only : pcols, pver
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in) ::  &
! Subprogram not used        ncol, &                           ! column count
! Subprogram not used        lchnk                             ! chunk index
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         calday                           ! calendar day of year
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         tn(pcols,pver), &                ! neutral gas temperature (K)
! Subprogram not used         o2(ncol,pver), &                 ! O2 concentration (kg/kg)
! Subprogram not used         o1(ncol,pver), &                 ! O concentration (kg/kg)
! Subprogram not used         mbar(ncol,pver)                  ! mean molecular weight (g/mole)
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         cpair(ncol,pver)                 ! specific heat capacity (J/K/kg)
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         pmid(pcols,pver)                 ! midpoint pressure (Pa)
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used         rlats(ncol), &                   ! column latitudes (radians)
! Subprogram not used         rlons(ncol)
! Subprogram not used       real(r8), intent(out) :: &
! Subprogram not used         aur_hrate(ncol,pver)             ! auroral heating rate
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: aur_therm     = 807._r8
! Subprogram not used       real(r8), parameter :: jkcal         = 4184._r8
! Subprogram not used       real(r8), parameter :: aur_heat_eff  = .05_r8
! Subprogram not used       real(r8), parameter :: aur_hconst    = 1.e3_r8*jkcal*aur_therm*aur_heat_eff
! Subprogram not used 
! Subprogram not used       integer  :: i, k
! Subprogram not used       integer  :: hemis(ncol)
! Subprogram not used       real(r8) :: r2d
! Subprogram not used       real(r8) :: ofda, cosofa, sinofa, aslona
! Subprogram not used       real(r8) :: sunlons(ncol)                     ! sun's mag longitudes
! Subprogram not used       real(r8) :: dlat_aur(ncol)
! Subprogram not used       real(r8) :: dlon_aur(ncol)
! Subprogram not used       real(r8) :: colat(ncol)
! Subprogram not used       real(r8) :: sinlat(ncol)
! Subprogram not used       real(r8) :: coslat(ncol)
! Subprogram not used       real(r8) :: coslon(ncol)
! Subprogram not used       real(r8) :: sinlon(ncol)
! Subprogram not used       real(r8) :: alon(ncol)
! Subprogram not used       real(r8) :: cusp(ncol)
! Subprogram not used       real(r8) :: alfa(ncol)
! Subprogram not used       real(r8) :: alfa2(ncol)
! Subprogram not used       real(r8) :: alfa3(ncol)
! Subprogram not used       real(r8) :: flux(ncol)
! Subprogram not used       real(r8) :: flux2(ncol)
! Subprogram not used       real(r8) :: flux3(ncol)
! Subprogram not used       real(r8) :: drizl(ncol)
! Subprogram not used       real(r8) :: qteaur(ncol)                         ! for electron temperature
! Subprogram not used       real(r8) :: qsum(ncol,pver)                      ! total ion production (1/s)
! Subprogram not used       logical  :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... initialize ion production
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       do k = 1,pver
! Subprogram not used         aur_hrate(:,k) = 0._r8
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       if (.not. has_ions) return
! Subprogram not used 
! Subprogram not used       r2d = 180._r8/pi
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... check latitudes, and return if all below 32.5 deg
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       do_aurora(:) = abs( rlats(:) ) > pi/6._r8
! Subprogram not used       if( all( .not. do_aurora(:) ) ) then
! Subprogram not used          return
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... set sun location
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call sunloc( calday, sunlons, lchnk, ncol )
! Subprogram not used 
! Subprogram not used       do i = 1,ncol
! Subprogram not used         if( do_aurora(i) ) then
! Subprogram not used           dlat_aur(i) = alatm(i,lchnk)
! Subprogram not used           dlon_aur(i) = alonm(i,lchnk) - sunlons(i)
! Subprogram not used           if( dlon_aur(i) > pi ) then
! Subprogram not used              dlon_aur(i) = dlon_aur(i) - twopi
! Subprogram not used           else if( dlon_aur(i) < -pi ) then
! Subprogram not used              dlon_aur(i) = dlon_aur(i) + twopi
! Subprogram not used           end if
! Subprogram not used           if( dlat_aur(i) > 0._r8 ) then
! Subprogram not used             hemis(i) = 2
! Subprogram not used           else
! Subprogram not used             hemis(i) = 1
! Subprogram not used           end if
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... find auroral circle coordinates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used             ofda      = sqrt( offa(hemis(i))**2 + dskofa(hemis(i))**2)
! Subprogram not used             cosofa    = cos( ofda )
! Subprogram not used             sinofa    = sin( ofda )
! Subprogram not used             aslona    = asin( dskofa(hemis(i))/ofda )
! Subprogram not used             sinlat(i) = sin( abs( dlat_aur(i) ) )
! Subprogram not used             coslat(i) = cos( dlat_aur(i) )
! Subprogram not used             sinlon(i) = sin( dlon_aur(i) + aslona )
! Subprogram not used             coslon(i) = cos( dlon_aur(i) + aslona )
! Subprogram not used             colat(i)  = acos( cosofa*sinlat(i) - sinofa*coslat(i)*coslon(i))
! Subprogram not used             alon(i)   = mod( atan2( sinlon(i)*coslat(i),sinlat(i)*sinofa &
! Subprogram not used                                     + cosofa*coslat(i)*coslon(i) ) - aslona + 3._r8*pi,twopi) - pi
! Subprogram not used         end if
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... make cusp
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... make alfa, flux, and drizzle
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call aurora_heat( flux, flux2, flux3, alfa, alfa2, &
! Subprogram not used                         alfa3, qteaur, drizl, do_aurora, hemis, &
! Subprogram not used                         alon, colat, ncol )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral additions to ionization rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       call total_ion_prod( drizl, cusp, alfa, alfa2, alfa3, &
! Subprogram not used                            flux, flux2, flux3, tn, o2, &
! Subprogram not used                            o1, mbar, qsum, pmid, do_aurora, &
! Subprogram not used                            ncol, lchnk )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... form auroral heating rate
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       do k = 1,pver
! Subprogram not used          aur_hrate(:,k) = aur_hconst * qsum(:,k) / (cpair(:,k) * mbar(:,k))
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine aurora_hrate

! Subprogram not used       subroutine aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... calculate horizontal variation of polar cusp heating
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)     :: ncol
! Subprogram not used       integer, intent(in)     :: hemis(ncol)
! Subprogram not used       real(r8), intent(in)    :: colat(ncol)
! Subprogram not used       real(r8), intent(in)    :: alon(ncol)
! Subprogram not used       real(r8), intent(out)   :: cusp(ncol)
! Subprogram not used       logical, intent(in)     :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: s5  =.08726646_r8, &
! Subprogram not used                              s20 =.34906585_r8
! Subprogram not used 
! Subprogram not used       where( do_aurora(:) )
! Subprogram not used          cusp(:) = (exp( -((theta0(hemis(:)) - colat(:))/s5)**2 ) &
! Subprogram not used                       + exp( -((pi - theta0(hemis(:)) - colat(:))/s5)**2) ) &
! Subprogram not used                         *exp( -(atan2( sin(alon(:) - phid(hemis(:))), cos(alon(:) - phid(hemis(:))) )/s20)**2 )
! Subprogram not used       elsewhere
! Subprogram not used          cusp(:) = 0._r8
! Subprogram not used       endwhere
! Subprogram not used 
! Subprogram not used       end subroutine aurora_cusp 

! Subprogram not used       subroutine aurora_heat( flux, flux2, flux3, alfa, alfa2, &
! Subprogram not used                               alfa3, qteaur, drizl, do_aurora, hemis, &
! Subprogram not used                               alon, colat, ncol )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... calculate alfa, flux, and drizzle
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)     :: ncol
! Subprogram not used       integer, intent(in)     :: hemis(ncol)
! Subprogram not used       real(r8), intent(in)    :: colat(ncol)
! Subprogram not used       real(r8), intent(in)    :: alon(ncol)
! Subprogram not used       real(r8), intent(inout) :: flux(ncol)
! Subprogram not used       real(r8), intent(inout) :: flux2(ncol)
! Subprogram not used       real(r8), intent(inout) :: flux3(ncol)
! Subprogram not used       real(r8), intent(inout) :: drizl(ncol)
! Subprogram not used       real(r8), intent(inout) :: qteaur(ncol)
! Subprogram not used       real(r8), intent(inout) :: alfa(ncol)
! Subprogram not used       real(r8), intent(inout) :: alfa2(ncol)
! Subprogram not used       real(r8), intent(inout) :: alfa3(ncol)
! Subprogram not used       logical, intent(in)     :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8), dimension(ncol) :: &
! Subprogram not used         coslamda, &                             ! cos(angle from throat)
! Subprogram not used         halfwidth, &                            ! oval half-width
! Subprogram not used         wrk, &                                  ! temp wrk array
! Subprogram not used         dtheta                                  ! latitudinal variation (Gaussian)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! Low-energy protons:
! Subprogram not used !
! Subprogram not used !     alfap0 = 0.5*(alfap1+alfap2)
! Subprogram not used !     e0p = 0.5*(pe1+pe2)
! Subprogram not used !
! Subprogram not used ! coslamda  = cos(lamda)
! Subprogram not used ! halfwidth = auroral half width
! Subprogram not used ! dtheta    = colat-theta0(ihem)
! Subprogram not used ! alfa      = electron energy
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       where( do_aurora(:) )
! Subprogram not used          coslamda(:) = cos( atan2( sin( alon(:) - rrote ),cos( alon(:) - rrote ) ) )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral oval half-width (equation (1) in Roble,1987):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used          halfwidth(:) = h0*(1._r8 - rh*cos( atan2( sin(alon(:) - rroth),cos( alon(:) - rroth ) ) ) )
! Subprogram not used          dtheta(:)    = colat(:) - rrad(hemis(:))
! Subprogram not used       endwhere
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... characteristic energy (equation (2) in Roble,1987):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       if( alfa0 > .01_r8 ) then
! Subprogram not used          where( do_aurora(:) )
! Subprogram not used             alfa(:) =  alfa0*(1._r8 - ralfa*coslamda(:))
! Subprogram not used          endwhere
! Subprogram not used       else
! Subprogram not used          alfa(:) =  0._r8
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       where( do_aurora(:) )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... flux, drizzle, alfa2, flux2
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... energy flux (equation (3) in Roble,1987):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used          wrk(:)   = exp( -(dtheta(:)/halfwidth(:))**2 )
! Subprogram not used          flux(:)  = e0*(1._r8 - ree*coslamda(:))*wrk(:) / (2._r8*alfa(:)*1.602e-9_r8)
! Subprogram not used          drizl(:) = exp( -((dtheta(:) + abs(dtheta(:)))/(2._r8*h0))**2 )
! Subprogram not used          alfa2(:) = alfa20*(1._r8 - ralfa2*coslamda(:))
! Subprogram not used          flux2(:) = e20*(1._r8 - re2*coslamda(:))*wrk(:) / (2._r8*alfa2(:)*1.602e-9_r8)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... alfa3, flux3 for high energy electrons:
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used          alfa3(:) = alfa30
! Subprogram not used          flux3(:) = e30*wrk(:) / 1.602e-6_r8
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... for electron temperature (used in settei):  
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used          qteaur(:) = -7.e8_r8*wrk(:)
! Subprogram not used       endwhere
! Subprogram not used 
! Subprogram not used       end subroutine aurora_heat

! Subprogram not used       subroutine aurora_ions( drizl, cusp, alfa1, alfa2, alfa3, &
! Subprogram not used                               flux1, flux2, flux3, tn, o2, &
! Subprogram not used                               o1, mbar, qo2p, qop, qn2p, &
! Subprogram not used                               qnp, pmid, do_aurora, ncol, lchnk, pbuf )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... calculate auroral additions to ionization rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use ppgrid,      only : pcols, pver
! Subprogram not used       use cam_history, only : outfld
! Subprogram not used 
! Subprogram not used       use physics_buffer,only: physics_buffer_desc, pbuf_set_field, pbuf_get_field
! Subprogram not used       
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in) :: ncol
! Subprogram not used       integer, intent(in) :: lchnk
! Subprogram not used       real(r8), intent(in), dimension(ncol) :: &
! Subprogram not used                              drizl, &
! Subprogram not used                              cusp, &
! Subprogram not used                              alfa1, &
! Subprogram not used                              alfa2, &
! Subprogram not used                              alfa3, &
! Subprogram not used                              flux1, &
! Subprogram not used                              flux2, &
! Subprogram not used                              flux3
! Subprogram not used       real(r8), dimension(pcols,pver), intent(in) :: &
! Subprogram not used                              tn, &                     ! midpoint neutral temperature (K)
! Subprogram not used                              pmid                      ! midpoint pressure (Pa)
! Subprogram not used       real(r8), dimension(ncol,pver), intent(in) :: &
! Subprogram not used                              o2, &                     ! midpoint o2 concentration (kg/kg)
! Subprogram not used                              o1, &                     ! midpoint o  concentration (kg/kg)
! Subprogram not used                              mbar                      ! mean molecular mass (g/mole)
! Subprogram not used       real(r8), dimension(ncol,pver), intent(inout) :: &
! Subprogram not used                              qo2p, &                   ! o2p prod from aurora (molecules/cm^3/s)
! Subprogram not used                              qop, &                    ! op prod from aurora (molecules/cm^3/s)
! Subprogram not used                              qn2p, &                   ! n2p prod from aurora (molecules/cm^3/s)
! Subprogram not used                              qnp                       ! np prod from aurora (molecules/cm^3/s)
! Subprogram not used       logical, intent(in) :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used       type(physics_buffer_desc),pointer :: pbuf(:)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: const0        = 1.e-20_r8
! Subprogram not used 
! Subprogram not used       integer  :: i, k
! Subprogram not used       real(r8), dimension(ncol) :: &
! Subprogram not used         p0ez, &
! Subprogram not used         press, &                                   ! pressure at interface levels (dyne/cm^2)
! Subprogram not used         tempi, &                                   ! temperature at interface levels (K)
! Subprogram not used         xalfa1, &
! Subprogram not used         xalfa2, &
! Subprogram not used         xcusp, &
! Subprogram not used         xdrizl, &                                  ! input to sub aion
! Subprogram not used         xalfa_sp, &
! Subprogram not used         xalfa3, &
! Subprogram not used         flux1_ion, &
! Subprogram not used         flux2_ion, &
! Subprogram not used         cusp_ion, &
! Subprogram not used         drizl_ion, &                               ! output from sub aion
! Subprogram not used         alfa1_ion, &
! Subprogram not used         alfa2_ion, &
! Subprogram not used         alfa3_ion, &                               ! output from sub aion
! Subprogram not used         alfasp_bion, &                             ! output from sub bion
! Subprogram not used         barm_t, &
! Subprogram not used         qsum, &
! Subprogram not used         denom, &
! Subprogram not used         p0ez_mbar, &
! Subprogram not used         tk_mbar, &
! Subprogram not used         barm, &
! Subprogram not used         falfa1, &
! Subprogram not used         falfa2, &
! Subprogram not used         fcusp, &
! Subprogram not used         fdrizl, &
! Subprogram not used         falfa_sp, &
! Subprogram not used         xn2, &
! Subprogram not used         falfa3
! Subprogram not used       real(r8), dimension(ncol) :: &
! Subprogram not used         qo2p_aur, &
! Subprogram not used         qop_aur, &
! Subprogram not used         qn2p_aur                                   ! auroral ionization for O2+, O+, N2+
! Subprogram not used       real(r8) :: qia(5)                           ! low energy proton source (not in use, 1/02)
! Subprogram not used       real(r8) :: wrk(ncol,pver)
! Subprogram not used 
! Subprogram not used       real(r8), pointer   :: aurIPRateSum(:,:) ! Pointer to pbuf auroral ion production sum for O2+,O+,N2+ (s-1 cm-3)
! Subprogram not used 
! Subprogram not used       qia(:) = 0._r8
! Subprogram not used       wrk(:,:) = 0._r8
! Subprogram not used  
! Subprogram not used       !-----------------------------------------------------------
! Subprogram not used       !  Point to production rates array in physics buffer where 
! Subprogram not used       !  rates will be stored for ionosphere module access.  Also, 
! Subprogram not used       !  initialize rates to zero before column loop since only 
! Subprogram not used       !  daylight values are filled
! Subprogram not used       !-----------------------------------------------------------
! Subprogram not used       if (indxAIPRS>0) then
! Subprogram not used         call pbuf_get_field(pbuf, indxAIPRS, aurIPRateSum)
! Subprogram not used         aurIPRateSum(:,:) = 0._r8
! Subprogram not used       endif
! Subprogram not used  
! Subprogram not used level_loop : &
! Subprogram not used       do k = 1,lev1
! Subprogram not used           where( do_aurora(:) )
! Subprogram not used              press(:ncol) = 10._r8*pmid(:ncol,k)              ! from Pa to dyne/cm^2
! Subprogram not used              tempi(:ncol) = tn(:ncol,k)
! Subprogram not used              barm(:)      = mbar(:,k)
! Subprogram not used              p0ez(:)      = (press(:)/(grav*4.e-6_r8))**.606_r8
! Subprogram not used              xalfa1(:)    = p0ez(:)/alfa1(:)
! Subprogram not used              xalfa2(:)    = p0ez(:)/alfa2(:)
! Subprogram not used              xcusp (:)    = p0ez(:)/alfac
! Subprogram not used              xdrizl(:)    = p0ez(:)/alfad
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... initialize (whole array operations):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used              flux1_ion(:) = const0
! Subprogram not used              flux2_ion(:) = const0
! Subprogram not used              alfa1_ion(:) = const0
! Subprogram not used              alfa2_ion(:) = const0
! Subprogram not used              cusp_ion(:)  = const0
! Subprogram not used              drizl_ion(:) = const0
! Subprogram not used           endwhere
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral electrons
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           call aion( xalfa1, alfa1_ion, do_aurora, ncol )
! Subprogram not used           call aion( xalfa2, alfa2_ion, do_aurora, ncol )
! Subprogram not used           call aion( xcusp , cusp_ion, do_aurora, ncol  )
! Subprogram not used           call aion( xdrizl, drizl_ion, do_aurora, ncol )
! Subprogram not used           where( do_aurora(:) )
! Subprogram not used              falfa1(:) = alfa1(:)*flux1(:)  ! s7
! Subprogram not used              falfa2(:) = alfa2(:)*flux2(:)  ! s8
! Subprogram not used              fcusp (:) = cusp(:)*alfac*fc   ! s9
! Subprogram not used              fdrizl(:) = drizl(:)*alfad*fd  ! s10
! Subprogram not used              qsum(:)   = falfa1(:)*alfa1_ion(:) &    ! s7*s3
! Subprogram not used                        + falfa2(:)*alfa2_ion(:) &    ! s8*s4
! Subprogram not used                        + fcusp(:)*cusp_ion (:) &     ! s9*s5
! Subprogram not used                        + fdrizl(:)*drizl_ion(:)       ! s10*s6
! Subprogram not used           endwhere
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... include solar protons if add_sproton is set, 
! Subprogram not used !           and high energy electrons if add_helectron is set
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... solar protons
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           if( add_sproton ) then
! Subprogram not used              if( flx_sp > 1.e-19_r8 ) then
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    p0ez_mbar(:) = press(:)/(boltz_cgs*tempi(:)*barm(:))/avo
! Subprogram not used                    tk_mbar(:)   = gask*tempi(:)/(grav*barm(:))*p0ez_mbar(:)
! Subprogram not used                    xalfa_sp(:)  = ((tk_mbar(:)/.00271_r8)**.58140_r8)/alfa_sp
! Subprogram not used                 endwhere
! Subprogram not used                 call bion( xalfa_sp, alfasp_bion, do_aurora, ncol )
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    falfa_sp(:) = drizl(:)*p0ez_mbar(:)*flx_sp*1.e6_r8/(tk_mbar(:)*35._r8)
! Subprogram not used                    qsum(:)     = qsum(:) + falfa_sp(:)*alfasp_bion(:)
! Subprogram not used                 endwhere
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... high energy electrons
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           if( add_helectron ) then
! Subprogram not used              if( e30 > 1.e-19_r8 ) then
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    xalfa3(:) = p0ez(:)/alfa3(:)        ! alfa3(:)==alfa30
! Subprogram not used                 endwhere
! Subprogram not used                 call aion( xalfa3, alfa3_ion, do_aurora, ncol )
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    falfa3(:) = alfa3(:)*flux3(:)  ! s13 (high energy electrons)
! Subprogram not used                    qsum(:)   = qsum(:) + falfa3(:)*alfa3_ion(:)
! Subprogram not used                 endwhere
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... form production
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           where( do_aurora(:) )
! Subprogram not used              barm_t(:) = grav*barm(:)/(35.e-3_r8*gask*tempi(:))
! Subprogram not used              qsum(:)   = qsum(:)*barm_t(:)               ! s1 = s1*s11
! Subprogram not used              wrk(:,k)  = qsum(:)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... denominator of equations (13-16) in Roble,1987.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used              xn2(:)   = max( (1._r8 - o2(:,k) - o1(:,k)),1.e-8_r8 )
! Subprogram not used              denom(:) = 0.92_r8*xn2(:)*rmassinv_n2 &
! Subprogram not used                       + 1.5_r8*o2(:,k) *rmassinv_o2 + 0.56_r8*o1(:,k) *rmassinv_o1
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... production of O2+ (equation (15) in Roble,1987):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used              qo2p_aur(:) = qsum(:)*o2(:,k)/(rmass_o2*denom(:)) + qia(2)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... production of O+ (equation (16) in Roble,1987):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used              qop_aur(:) = qsum(:)*(.5_r8 *o2(:,k)*rmassinv_o2 &
! Subprogram not used                                    + .56_r8*o1(:,k)*rmassinv_o1)/denom(:) + qia(3)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... production of N2+ (equation (13) in Roble,1987)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used              qn2p_aur(:) = qsum(:)*.7_r8*xn2(:)/(rmass_n2*denom(:)) + qia(1)
! Subprogram not used              qo2p(:,k)   = qo2p(:,k) + qo2p_aur(:)
! Subprogram not used              qop(:,k)    = qop(:,k) + qop_aur(:)
! Subprogram not used              qn2p(:,k)   = qn2p(:,k) + qn2p_aur(:)
! Subprogram not used              qnp(:,k)    = qnp (:,k) + .22_r8/.7_r8 * qn2p_aur(:)
! Subprogram not used           endwhere
! Subprogram not used       end do level_loop
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       !  Store the sum of the ion production rates in pbuf to be used 
! Subprogram not used       !  in the ionosx module 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       if (indxAIPRS>0) then
! Subprogram not used       
! Subprogram not used         aurIPRateSum(1:ncol,1:pver) = wrk(1:ncol,1:pver) 
! Subprogram not used       
! Subprogram not used       endif
! Subprogram not used   
! Subprogram not used       call outfld( 'QSUM', wrk, ncol, lchnk )
! Subprogram not used 
! Subprogram not used       end subroutine aurora_ions

! Subprogram not used       subroutine total_ion_prod( drizl, cusp, alfa1, alfa2, alfa3, &
! Subprogram not used                                  flux1, flux2, flux3, tn, o2, &
! Subprogram not used                                  o1, mbar, tpions, pmid, do_aurora, &
! Subprogram not used                                  ncol, lchnk )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... calculate auroral additions to ionization rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use ppgrid,      only : pcols, pver
! Subprogram not used       use cam_history, only : outfld
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in) :: ncol
! Subprogram not used       integer, intent(in) :: lchnk
! Subprogram not used       real(r8), intent(in), dimension(ncol) :: &
! Subprogram not used                              drizl, &
! Subprogram not used                              cusp, &
! Subprogram not used                              alfa1, &
! Subprogram not used                              alfa2, &
! Subprogram not used                              alfa3, &
! Subprogram not used                              flux1, &
! Subprogram not used                              flux2, &
! Subprogram not used                              flux3
! Subprogram not used       real(r8), dimension(pcols,pver), intent(in) :: &
! Subprogram not used                              tn, &                     ! midpoint neutral temperature (K)
! Subprogram not used                              pmid                      ! midpoint pressure (Pa)
! Subprogram not used       real(r8), dimension(ncol,pver), intent(in) :: &
! Subprogram not used                              o2, &                     ! midpoint o2 concentration (kg/kg)
! Subprogram not used                              o1, &                     ! midpoint o  concentration (kg/kg)
! Subprogram not used                              mbar                      ! mean molecular mass (g/mole)
! Subprogram not used       real(r8), dimension(ncol,pver), intent(inout) :: &
! Subprogram not used                              tpions                    ! total ion production (1/s)
! Subprogram not used       logical, intent(in) :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: const0        = 1.e-20_r8
! Subprogram not used 
! Subprogram not used       integer  :: i, k
! Subprogram not used       real(r8), dimension(ncol) :: &
! Subprogram not used         p0ez, &
! Subprogram not used         press, &                                   ! pressure at interface levels (dyne/cm^2)
! Subprogram not used         tempi, &                                   ! temperature at interface levels (K)
! Subprogram not used         xalfa1, &
! Subprogram not used         xalfa2, &
! Subprogram not used         xcusp, &
! Subprogram not used         xdrizl, &                                  ! input to sub aion
! Subprogram not used         xalfa_sp, &
! Subprogram not used         xalfa3, &
! Subprogram not used         flux1_ion, &
! Subprogram not used         flux2_ion, &
! Subprogram not used         cusp_ion, &
! Subprogram not used         drizl_ion, &                               ! output from sub aion
! Subprogram not used         alfa1_ion, &
! Subprogram not used         alfa2_ion, &
! Subprogram not used         alfa3_ion, &                               ! output from sub aion
! Subprogram not used         alfasp_bion, &                             ! output from sub bion
! Subprogram not used         barm_t, &
! Subprogram not used         qsum, &
! Subprogram not used         denom, &
! Subprogram not used         p0ez_mbar, &
! Subprogram not used         tk_mbar, &
! Subprogram not used         barm, &
! Subprogram not used         falfa1, &
! Subprogram not used         falfa2, &
! Subprogram not used         fcusp, &
! Subprogram not used         fdrizl, &
! Subprogram not used         falfa_sp, &
! Subprogram not used         xn2, &
! Subprogram not used         falfa3
! Subprogram not used       real(r8), dimension(ncol) :: &
! Subprogram not used         qo2p_aur, &
! Subprogram not used         qop_aur, &
! Subprogram not used         qn2p_aur                                   ! auroral ionization for O2+, O+, N2+
! Subprogram not used       real(r8) :: qia(5)                           ! low energy proton source (not in use, 1/02)
! Subprogram not used 
! Subprogram not used       qia(:)      = 0._r8
! Subprogram not used       tpions(:,:) = 0._r8
! Subprogram not used 
! Subprogram not used level_loop : &
! Subprogram not used       do k = 1,lev1
! Subprogram not used           where( do_aurora(:) )
! Subprogram not used              press(:ncol) = 10._r8*pmid(:ncol,k)              ! from Pa to dyne/cm^2
! Subprogram not used              tempi(:ncol) = tn(:ncol,k)
! Subprogram not used              barm(:)      = mbar(:,k)
! Subprogram not used              p0ez(:)      = (press(:)/(grav*4.e-6_r8))**.606_r8
! Subprogram not used              xalfa1(:)    = p0ez(:)/alfa1(:)
! Subprogram not used              xalfa2(:)    = p0ez(:)/alfa2(:)
! Subprogram not used              xcusp (:)    = p0ez(:)/alfac
! Subprogram not used              xdrizl(:)    = p0ez(:)/alfad
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... initiliaze (whole array operations):
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used              flux1_ion(:) = const0
! Subprogram not used              flux2_ion(:) = const0
! Subprogram not used              alfa1_ion(:) = const0
! Subprogram not used              alfa2_ion(:) = const0
! Subprogram not used              cusp_ion(:)  = const0
! Subprogram not used              drizl_ion(:) = const0
! Subprogram not used           endwhere
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... auroral electrons
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           call aion( xalfa1, alfa1_ion, do_aurora, ncol )
! Subprogram not used           call aion( xalfa2, alfa2_ion, do_aurora, ncol )
! Subprogram not used           call aion( xcusp , cusp_ion, do_aurora, ncol  )
! Subprogram not used           call aion( xdrizl, drizl_ion, do_aurora, ncol )
! Subprogram not used           where( do_aurora(:) )
! Subprogram not used              falfa1(:) = alfa1(:)*flux1(:)  ! s7
! Subprogram not used              falfa2(:) = alfa2(:)*flux2(:)  ! s8
! Subprogram not used              fcusp (:) = cusp(:)*alfac*fc   ! s9
! Subprogram not used              fdrizl(:) = drizl(:)*alfad*fd  ! s10
! Subprogram not used              qsum(:)   = falfa1(:)*alfa1_ion(:) &    ! s7*s3
! Subprogram not used                        + falfa2(:)*alfa2_ion(:) &    ! s8*s4
! Subprogram not used                        + fcusp(:)*cusp_ion (:) &     ! s9*s5
! Subprogram not used                        + drizl(:)*drizl_ion(:)       ! s10*s6
! Subprogram not used           endwhere
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... include solar protons if add_sproton is set, 
! Subprogram not used !           and high energy electrons if add_helectron is set
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... solar protons
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           if( add_sproton ) then
! Subprogram not used              if( flx_sp > 1.e-19_r8 ) then
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    p0ez_mbar(:) = press(:)/(boltz_cgs*tempi(:)*barm(:))/avo
! Subprogram not used                    tk_mbar(:)   = gask*tempi(:)/(grav*barm(:))*p0ez_mbar(:)
! Subprogram not used                    xalfa_sp(:)  = ((tk_mbar(:)/.00271_r8)**.58140_r8)/alfa_sp
! Subprogram not used                 endwhere
! Subprogram not used                 call bion( xalfa_sp, alfasp_bion, do_aurora, ncol )
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    falfa_sp(:) = drizl(:)*p0ez_mbar(:)*flx_sp*1.e6_r8/(tk_mbar(:)*35._r8)
! Subprogram not used                    qsum(:)     = qsum(:) + falfa_sp(:)*alfasp_bion(:)
! Subprogram not used                 endwhere
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... high energy electrons
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           if( add_helectron ) then
! Subprogram not used              if( e30 > 1.e-19_r8 ) then
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    xalfa3(:) = p0ez(:)/alfa3(:)        ! alfa3(:)==alfa30
! Subprogram not used                 endwhere
! Subprogram not used                 call aion( xalfa3, alfa3_ion, do_aurora, ncol )
! Subprogram not used                 where( do_aurora(:) )
! Subprogram not used                    falfa3(:) = alfa3(:)*flux3(:)  ! s13 (high energy electrons)
! Subprogram not used                    qsum(:)   = qsum(:) + falfa3(:)*alfa3_ion(:)
! Subprogram not used                 endwhere
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... form production
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used           where( do_aurora(:) )
! Subprogram not used              barm_t(:)   = grav*barm(:)/(35.e-3_r8*gask*tempi(:))
! Subprogram not used              tpions(:,k) = qsum(:)*barm_t(:)               ! s1 = s1*s11
! Subprogram not used           endwhere
! Subprogram not used       end do level_loop
! Subprogram not used 
! Subprogram not used       end subroutine total_ion_prod

! Subprogram not used       subroutine aion( si, so, do_aurora, ncol )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! Calculates integrated f(x) needed for total auroral ionization.
! Subprogram not used ! See equations (10-12) in Roble,1987.
! Subprogram not used ! Coefficients for equation (12) of Roble,1987 are in variable cc 
! Subprogram not used ! (revised since 1987):
! Subprogram not used ! Uses the identity x**y = exp(y*ln(x)) for performance 
! Subprogram not used ! (fewer (1/2) trancendental functions are required).
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used       integer,  intent(in)  :: ncol
! Subprogram not used       real(r8), intent(in)  :: si(ncol)
! Subprogram not used       real(r8), intent(out) :: so(ncol)
! Subprogram not used       logical,  intent(in)  :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: cc(8) = &
! Subprogram not used        (/ 3.2333134511131_r8 ,  2.5658873458085_r8 ,  2.2540957232641_r8 , &
! Subprogram not used           0.72971983372673_r8,  1.1069072431948_r8 ,  1.7134937681128_r8 , &
! Subprogram not used           1.8835442312993_r8 ,  0.86472135072090_r8 /)
! Subprogram not used 
! Subprogram not used       real(r8) :: xlog(ncol)
! Subprogram not used 
! Subprogram not used       where( do_aurora(:) )
! Subprogram not used          xlog(:) = log( si(:) )
! Subprogram not used          so(:)   = cc(1)*exp( cc(2)*xlog(:) - cc(3)*exp( cc(4)*xlog(:) ) ) &
! Subprogram not used                    + cc(5)*exp( cc(6)*xlog(:) - cc(7)*exp( cc(8)*xlog(:) ) )
! Subprogram not used       elsewhere
! Subprogram not used          so(:) = 0._r8
! Subprogram not used       endwhere
! Subprogram not used 
! Subprogram not used       end subroutine aion

! Subprogram not used       subroutine bion( si, so, do_aurora, ncol )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! Calculates integrated f(x) needed for total auroral ionization.
! Subprogram not used ! See equations (10-12) in Roble,1987.
! Subprogram not used ! Use the identity x**y = exp(y*ln(x)) for performance 
! Subprogram not used ! (fewer (1/2) trancendental functions are required).
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol
! Subprogram not used       real(r8), intent(in)  :: si(ncol)
! Subprogram not used       real(r8), intent(out) :: so(ncol)
! Subprogram not used       logical,  intent(in)  :: do_aurora(ncol)
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: cc(8) = &
! Subprogram not used         (/ 0.12718_r8, 4.9119_r8, 1.8429_r8, 0.99336_r8, 0.52472_r8, &
! Subprogram not used            1.5565_r8,  .85732_r8, 1.4116_r8 /)
! Subprogram not used 
! Subprogram not used       real(r8) :: xlog(ncol)
! Subprogram not used 
! Subprogram not used       where( do_aurora(:) )
! Subprogram not used          xlog(:) = log( si(:) )
! Subprogram not used          so(:)   = cc(1)*exp( cc(2)*xlog(:) - cc(3)*exp( cc(4)*xlog(:) ) ) &
! Subprogram not used                    + cc(5)*exp( cc(6)*xlog(:) - cc(7)*exp( cc(8)*xlog(:) ) )
! Subprogram not used       elsewhere
! Subprogram not used          so(:) = 0._r8
! Subprogram not used       endwhere
! Subprogram not used 
! Subprogram not used       end subroutine bion

! Subprogram not used       subroutine sunloc( calday, sunlons, lchnk, ncol )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... calculate sun's longitude in dipole coordinates, defining sunlon
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       use dyn_grid,   only : get_dyn_grid_parm, get_horiz_grid_d
! Subprogram not used       use spmd_utils, only : masterproc
! Subprogram not used       use  phys_grid, only : get_lat_all_p, get_rlat_all_p, get_rlon_all_p
! Subprogram not used       use  mo_apex,   only : glonm                                              ! magnetic longitude grid (radians)
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer,  intent(in)   :: ncol
! Subprogram not used       integer,  intent(in)   :: lchnk
! Subprogram not used       real(r8), intent(in)   :: calday  ! calendar day of year
! Subprogram not used       real(r8), intent(out)  :: sunlons(ncol)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer  :: col
! Subprogram not used       integer  :: latndx(ncol)
! Subprogram not used       integer  :: sunlon_ndx
! Subprogram not used       integer  :: sunlon_ndxp1
! Subprogram not used       integer :: plon, plat
! Subprogram not used       real(r8), allocatable :: clon(:)
! Subprogram not used       real(r8) :: wght1, wght2, dellon, r2d
! Subprogram not used       real(r8) :: rlats(ncol), rlons(ncol)
! Subprogram not used 
! Subprogram not used       r2d = 180._r8/pi
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... sun's geographic coordinates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       plon = get_dyn_grid_parm('plon')
! Subprogram not used       plat = get_dyn_grid_parm('plat')
! Subprogram not used       allocate(clon(plon))
! Subprogram not used       call get_horiz_grid_d(plon,clon_d_out=clon)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       dellon       = .5_r8 - (calday - int(calday))
! Subprogram not used       sunlon_ndx   = mod( nint( dellon*plon ) - 1,plon ) + 1
! Subprogram not used       if( sunlon_ndx < 1 ) then
! Subprogram not used          sunlon_ndx = plon + sunlon_ndx
! Subprogram not used       end if
! Subprogram not used       sunlon_ndxp1 = mod( sunlon_ndx,plon ) + 1
! Subprogram not used       wght2        = min( 1._r8, max( (dellon*twopi - clon(sunlon_ndx))*plon/twopi,0._r8 ) )
! Subprogram not used       wght1        = 1._r8 - wght2
! Subprogram not used       deallocate(clon)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... get chunck latitudes
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call get_lat_all_p( lchnk, ncol, latndx )
! Subprogram not used       call get_rlat_all_p( lchnk, ncol, rlats )
! Subprogram not used       call get_rlon_all_p( lchnk, ncol, rlons )
! Subprogram not used 
! Subprogram not used       do col = 1,ncol
! Subprogram not used !        sunlons(col) = wght1*glonm(sunlon_ndx,latndx(col)) + wght2*glonm(sunlon_ndxp1,latndx(col))
! Subprogram not used 	 sunlons(col) = wght1*glonm(sunlon_ndx,plat/2) + wght2*glonm(sunlon_ndxp1,plat/2)
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       end subroutine sunloc

      end module mo_aurora
