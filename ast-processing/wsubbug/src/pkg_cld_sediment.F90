
module pkg_cld_sediment

!---------------------------------------------------------------------------------
! Purpose:
!
! Contains routines to compute tendencies from sedimentation of cloud liquid and 
! ice particles
!
! Author: Byron Boville  Sept 19, 2002 from code by P. J. Rasch
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, latvap, latice, rair, rhoh2o
  use cldwat,        only: icritc
  use pkg_cldoptics, only: reitab, reltab
  use abortutils,    only: endrun
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: cld_sediment_readnl, cld_sediment_vel, cld_sediment_tend


  real (r8), parameter :: vland  = 1.5_r8            ! liquid fall velocity over land  (cm/s)
  real (r8), parameter :: vocean = 2.8_r8            ! liquid fall velocity over ocean (cm/s)
  real (r8), parameter :: mxsedfac   = 0.99_r8       ! maximum sedimentation flux factor

  logical,   parameter :: stokes = .true.         ! use Stokes velocity instead of McFarquhar and Heymsfield

! parameter for modified McFarquhar and Heymsfield
  real (r8), parameter :: vice_small = 1._r8         ! ice fall velocity for small concentration (cm/s)

! parameters for Stokes velocity
  real (r8), parameter :: eta =  1.7e-5_r8           ! viscosity of air (kg m / s)
  real (r8), parameter :: r40 =  40._r8              !  40 micron radius
  real (r8), parameter :: r400= 400._r8              ! 400 micron radius
  real (r8), parameter :: v400= 1.00_r8              ! fall velocity of 400 micron sphere (m/s)
  real (r8)            :: v40 ! = (2._r8/9._r8) * rhoh2o * gravit/eta * r40**2 * 1.e-12_r8  
                                                     ! Stokes fall velocity of 40 micron sphere (m/s)
  real (r8)            :: vslope !  = (v400 - v40)/(r400 -r40) ! linear slope for large particles m/s/micron

  ! namelist variables
  real(r8) :: cldsed_ice_stokes_fac = huge(1._r8)    ! factor applied to the ice fall velocity computed from 
                                                     ! stokes terminal velocity

!===============================================================================
contains
!===============================================================================

subroutine cld_sediment_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cld_sediment_readnl'

   namelist /cldsed_nl/ cldsed_ice_stokes_fac
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldsed_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldsed_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      write(iulog,*) subname//': cldsed_ice_stokes_fac = ', cldsed_ice_stokes_fac

   end if


   ! Broadcast namelist variables
   call mpibcast(cldsed_ice_stokes_fac, 1, mpir8, 0, mpicom)


end subroutine cld_sediment_readnl

!===============================================================================

! Subprogram not used   subroutine cld_sediment_vel (ncol,                               &
! Subprogram not used        icefrac , landfrac, ocnfrac , pmid    , pdel    , t       , &
! Subprogram not used        cloud   , cldliq  , cldice  , pvliq   , pvice   , landm, snowh)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used ! Compute gravitational sedimentation velocities for cloud liquid water
! Subprogram not used ! and ice, based on Lawrence and Crutzen (1998).
! Subprogram not used 
! Subprogram not used ! LIQUID
! Subprogram not used 
! Subprogram not used ! The fall velocities assume that droplets have a gamma distribution
! Subprogram not used ! with effective radii for land and ocean as assessed by Han et al.;
! Subprogram not used ! see Lawrence and Crutzen (1998) for a derivation.
! Subprogram not used 
! Subprogram not used ! ICE
! Subprogram not used 
! Subprogram not used ! The fall velocities are based on data from McFarquhar and Heymsfield
! Subprogram not used ! or on Stokes terminal velocity for spheres and the effective radius.
! Subprogram not used 
! Subprogram not used ! NEED TO BE CAREFUL - VELOCITIES SHOULD BE AT THE *LOWER* INTERFACE
! Subprogram not used ! (THAT IS, FOR K+1), FLUXES ARE ALSO AT THE LOWER INTERFACE (K+1), 
! Subprogram not used ! BUT MIXING RATIOS ARE AT THE MIDPOINTS (K)...
! Subprogram not used 
! Subprogram not used ! NOTE THAT PVEL IS ON PVERP (INTERFACES), WHEREAS VFALL IS FOR THE CELL
! Subprogram not used ! AVERAGES (I.E., MIDPOINTS); ASSUME THE FALL VELOCITY APPLICABLE TO THE 
! Subprogram not used ! LOWER INTERFACE (K+1) IS THE SAME AS THAT APPLICABLE FOR THE CELL (V(K))
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     MATCH-MPIC version 2.0, Author: mgl, March 1998
! Subprogram not used ! adapted by P. J. Rasch
! Subprogram not used !            B. A. Boville, September 19, 2002
! Subprogram not used !            P. J. Rasch    May 22, 2003 (added stokes flow calc for liquid
! Subprogram not used !                                         drops based on effect radii)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used     integer, intent(in) :: ncol                     ! number of colums to process
! Subprogram not used 
! Subprogram not used     real(r8), intent(in)  :: icefrac (pcols)        ! sea ice fraction (fraction)
! Subprogram not used     real(r8), intent(in)  :: landfrac(pcols)        ! land fraction (fraction)
! Subprogram not used     real(r8), intent(in)  :: ocnfrac (pcols)        ! ocean fraction (fraction)
! Subprogram not used     real(r8), intent(in)  :: pmid  (pcols,pver)     ! pressure of midpoint levels (Pa)
! Subprogram not used     real(r8), intent(in)  :: pdel  (pcols,pver)     ! pressure diff across layer (Pa)
! Subprogram not used     real(r8), intent(in)  :: cloud (pcols,pver)     ! cloud fraction (fraction)
! Subprogram not used     real(r8), intent(in)  :: t     (pcols,pver)     ! temperature (K)
! Subprogram not used     real(r8), intent(in)  :: cldliq(pcols,pver)     ! cloud water, liquid (kg/kg)
! Subprogram not used     real(r8), intent(in)  :: cldice(pcols,pver)     ! cloud water, ice    (kg/kg)
! Subprogram not used     real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
! Subprogram not used 
! Subprogram not used     real(r8), intent(out) :: pvliq (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)
! Subprogram not used     real(r8), intent(out) :: pvice (pcols,pverp)    ! vertical velocity of cloud ice particles (Pa/s)
! Subprogram not used     real(r8), intent(in) :: landm(pcols)            ! land fraction ramped over water
! Subprogram not used ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used     real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
! Subprogram not used     real (r8) :: vfall                              ! settling velocity of cloud particles (m/s)
! Subprogram not used     real (r8) :: icice                              ! in cloud ice water content (kg/kg)
! Subprogram not used     real (r8) :: iciwc                              ! in cloud ice water content in g/m3
! Subprogram not used     real (r8) :: icefac
! Subprogram not used     real (r8) :: logiwc
! Subprogram not used 
! Subprogram not used     real (r8) :: rei(pcols,pver)                    ! effective radius of ice particles (microns)
! Subprogram not used     real (r8) :: rel(pcols,pver)                    ! effective radius of liq particles (microns)
! Subprogram not used     real(r8)  pvliq2 (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)
! Subprogram not used 
! Subprogram not used     integer i,k
! Subprogram not used 
! Subprogram not used     real (r8) :: lbound, ac, bc, cc
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------- initialize linear ramp variables for fall velocity ------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   v40 = (2._r8/9._r8) * rhoh2o * gravit/eta * r40**2 * 1.e-12_r8  
! Subprogram not used   vslope = (v400 - v40)/(r400 -r40)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !--------------------- liquid fall velocity ----------------------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! compute air density
! Subprogram not used     rho(:ncol,:) = pmid(:ncol,:) / (rair * t(:ncol,:))
! Subprogram not used 
! Subprogram not used     pvliq(:ncol,:) = 0._r8
! Subprogram not used 
! Subprogram not used     ! get effective radius of liquid drop
! Subprogram not used     call reltab(ncol, t, landfrac, landm, icefrac, rel, snowh)
! Subprogram not used 
! Subprogram not used     do k = 1,pver
! Subprogram not used        do i = 1,ncol
! Subprogram not used           if (cloud(i,k) > 0._r8 .and. cldliq(i,k) > 0._r8) then
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! newway
! Subprogram not used              if (rel(i,k) < 40._r8 ) then
! Subprogram not used                 vfall = 2._r8/9._r8 * rhoh2o * gravit * rel(i,k)**2 / eta  * 1.e-12_r8  ! micons^2 -> m^2
! Subprogram not used              else
! Subprogram not used                 vfall = v40 + vslope * (rel(i,k)-r40)      ! linear above 40 microns
! Subprogram not used              end if
! Subprogram not used              ! convert the fall speed to pressure units
! Subprogram not used              ! but do not apply the traditional
! Subprogram not used              ! negative convention for pvel.
! Subprogram not used !             pvliq2(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
! Subprogram not used              pvliq(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !--------------------- ice fall velocity -------------------------------
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     pvice(:ncol,:) = 0._r8
! Subprogram not used 
! Subprogram not used     if (stokes) then
! Subprogram not used 
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        !--------------------- stokes terminal velocity < 40 microns -----------
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used        ! get effective radius
! Subprogram not used        call reitab(ncol, t, rei)
! Subprogram not used 
! Subprogram not used        do k = 1,pver
! Subprogram not used           do i = 1,ncol
! Subprogram not used              if (cloud(i,k) > 0._r8 .and. cldice(i,k) > 0._r8) then
! Subprogram not used                 if (rei(i,k) < 40._r8 ) then
! Subprogram not used                    vfall = 2._r8/9._r8 * rhoh2o * gravit * rei(i,k)**2 / eta  * 1.e-12_r8  ! micons^2 -> m^2
! Subprogram not used                    vfall = vfall * cldsed_ice_stokes_fac
! Subprogram not used                 else
! Subprogram not used                    vfall = v40 + vslope * (rei(i,k)-r40)      ! linear above 40 microns
! Subprogram not used                 end if
! Subprogram not used 
! Subprogram not used                 ! convert the fall speed to pressure units, but do not apply the traditional
! Subprogram not used                 ! negative convention for pvel.
! Subprogram not used                 pvice(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        !--------------------- McFarquhar and Heymsfield > icritc --------------
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used        ! lower bound for iciwc
! Subprogram not used 
! Subprogram not used        cc = 128.64_r8 
! Subprogram not used        bc = 53.242_r8
! Subprogram not used        ac = 5.4795_r8
! Subprogram not used        lbound = (-bc + sqrt(bc*bc-4*ac*cc))/(2*ac)
! Subprogram not used        lbound = 10._r8**lbound
! Subprogram not used 
! Subprogram not used        do k = 1,pver
! Subprogram not used           do i = 1,ncol
! Subprogram not used              if (cloud(i,k) > 0._r8 .and. cldice(i,k) > 0._r8) then
! Subprogram not used 
! Subprogram not used                 ! compute the in-cloud ice concentration (kg/kg)
! Subprogram not used                 icice = cldice(i,k) / cloud(i,k)
! Subprogram not used 
! Subprogram not used                 ! compute the ice water content in g/m3
! Subprogram not used                 iciwc = icice * rho(i,k) * 1.e3_r8
! Subprogram not used 
! Subprogram not used                 ! set the fall velocity (cm/s) to depend on the ice water content in g/m3,
! Subprogram not used                 if (iciwc > lbound) then ! need this because of log10
! Subprogram not used                    logiwc = log10(iciwc)
! Subprogram not used                    !          Median - 
! Subprogram not used                    vfall = 128.64_r8 + 53.242_r8*logiwc + 5.4795_r8*logiwc**2
! Subprogram not used                    !          Average - 
! Subprogram not used                 !!$             vfall = 122.63 + 44.111*logiwc + 4.2144*logiwc**2
! Subprogram not used                 else
! Subprogram not used                    vfall = 0._r8
! Subprogram not used                 end if
! Subprogram not used 
! Subprogram not used                 ! set ice velocity to 1 cm/s if ice mixing ratio < icritc, ramp to value
! Subprogram not used                 ! calculated above at 2*icritc
! Subprogram not used                 if (icice <= icritc) then
! Subprogram not used                    vfall = vice_small
! Subprogram not used                 else if(icice < 2*icritc) then
! Subprogram not used                    icefac = (icice-icritc)/icritc
! Subprogram not used                    vfall = vice_small * (1._r8-icefac) + vfall * icefac
! Subprogram not used                 end if
! Subprogram not used 
! Subprogram not used                 ! bound the terminal velocity of ice particles at high concentration
! Subprogram not used                 vfall = min(100.0_r8, vfall)
! Subprogram not used 
! Subprogram not used                 ! convert the fall speed to pressure units, but do not apply the traditional
! Subprogram not used                 ! negative convention for pvel.
! Subprogram not used                 pvice(i,k+1) = vfall     &
! Subprogram not used                    * 0.01_r8                 & ! cm to meters
! Subprogram not used                    * rho(i,k)*gravit        ! meters/sec to pascals/sec
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used  end subroutine cld_sediment_vel


!===============================================================================
! Subprogram not used   subroutine cld_sediment_tend (ncol, dtime  ,               &
! Subprogram not used        pint   , pmid   , pdel   , t      ,                   &
! Subprogram not used        cloud  , cldliq , cldice , pvliq  , pvice  ,          &
! Subprogram not used        liqtend, icetend, wvtend , htend  , sfliq  , sfice   )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     Apply Cloud Particle Gravitational Sedimentation to Condensate
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used     integer,  intent(in)  :: ncol                      ! number of colums to process
! Subprogram not used 
! Subprogram not used     real(r8), intent(in)  :: dtime                     ! time step
! Subprogram not used     real(r8), intent(in)  :: pint  (pcols,pverp)       ! interfaces pressure (Pa)
! Subprogram not used     real(r8), intent(in)  :: pmid  (pcols,pver)        ! midpoint pressures (Pa)
! Subprogram not used     real(r8), intent(in)  :: pdel  (pcols,pver)        ! pressure diff across layer (Pa)
! Subprogram not used     real(r8), intent(in)  :: cloud (pcols,pver)        ! cloud fraction (fraction)
! Subprogram not used     real(r8), intent(in)  :: t     (pcols,pver)        ! temperature (K)
! Subprogram not used     real(r8), intent(in)  :: cldliq(pcols,pver)        ! cloud liquid water (kg/kg)
! Subprogram not used     real(r8), intent(in)  :: cldice(pcols,pver)        ! cloud ice water    (kg/kg)
! Subprogram not used     real(r8), intent(in)  :: pvliq (pcols,pverp)       ! vertical velocity of liquid drops  (Pa/s)
! Subprogram not used     real(r8), intent(in)  :: pvice (pcols,pverp)       ! vertical velocity of ice particles (Pa/s)
! Subprogram not used ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))
! Subprogram not used 
! Subprogram not used     real(r8), intent(out) :: liqtend(pcols,pver)       ! liquid condensate tend
! Subprogram not used     real(r8), intent(out) :: icetend(pcols,pver)       ! ice condensate tend
! Subprogram not used     real(r8), intent(out) :: wvtend (pcols,pver)       ! water vapor tend
! Subprogram not used     real(r8), intent(out) :: htend  (pcols,pver)       ! heating rate
! Subprogram not used     real(r8), intent(out) :: sfliq  (pcols)            ! surface flux of liquid (rain, kg/m/s)
! Subprogram not used     real(r8), intent(out) :: sfice  (pcols)            ! surface flux of ice    (snow, kg/m/s)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used     real(r8) :: fxliq(pcols,pverp)                     ! fluxes at the interfaces, liquid (positive = down)
! Subprogram not used     real(r8) :: fxice(pcols,pverp)                     ! fluxes at the interfaces, ice    (positive = down)
! Subprogram not used     real(r8) :: cldab(pcols)                           ! cloud in layer above
! Subprogram not used     real(r8) :: evapliq                                ! evaporation of cloud liquid into environment
! Subprogram not used     real(r8) :: evapice                                ! evaporation of cloud ice into environment
! Subprogram not used     real(r8) :: cldovrl                                ! cloud overlap factor
! Subprogram not used 
! Subprogram not used     integer :: i,k
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used ! initialize variables
! Subprogram not used     fxliq  (:ncol,:) = 0._r8 ! flux at interfaces (liquid)
! Subprogram not used     fxice  (:ncol,:) = 0._r8 ! flux at interfaces (ice)
! Subprogram not used     liqtend(:ncol,:) = 0._r8 ! condensate tend (liquid)
! Subprogram not used     icetend(:ncol,:) = 0._r8 ! condensate tend (ice)
! Subprogram not used     wvtend(:ncol,:)  = 0._r8 ! environmental moistening
! Subprogram not used     htend(:ncol,:)   = 0._r8 ! evaporative cooling
! Subprogram not used     sfliq(:ncol)     = 0._r8 ! condensate sedimentation flux out bot of column (liquid)
! Subprogram not used     sfice(:ncol)     = 0._r8 ! condensate sedimentation flux out bot of column (ice)
! Subprogram not used 
! Subprogram not used ! fluxes at interior points
! Subprogram not used     call getflx(ncol, pint, cldliq, pvliq, dtime, fxliq)
! Subprogram not used     call getflx(ncol, pint, cldice, pvice, dtime, fxice)
! Subprogram not used 
! Subprogram not used ! calculate fluxes at boundaries
! Subprogram not used     do i = 1,ncol
! Subprogram not used        fxliq(i,1) = 0._r8
! Subprogram not used        fxice(i,1) = 0._r8
! Subprogram not used ! surface flux by upstream scheme
! Subprogram not used        fxliq(i,pverp) = cldliq(i,pver) * pvliq(i,pverp) * dtime
! Subprogram not used        fxice(i,pverp) = cldice(i,pver) * pvice(i,pverp) * dtime
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used ! filter out any negative fluxes from the getflx routine
! Subprogram not used ! (typical fluxes are of order > 1e-3 when clouds are present)
! Subprogram not used     do k = 2,pver
! Subprogram not used        fxliq(:ncol,k) = max(0._r8, fxliq(:ncol,k))
! Subprogram not used        fxice(:ncol,k) = max(0._r8, fxice(:ncol,k))
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used ! Limit the flux out of the bottom of each cell to the water content in each phase.
! Subprogram not used ! Apply mxsedfac to prevent generating very small negative cloud water/ice
! Subprogram not used ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
! Subprogram not used ! ***Should we include the flux in the top, to allow for thin surface layers?
! Subprogram not used ! ***Requires simple treatment of cloud overlap, already included below.
! Subprogram not used     do k = 1,pver
! Subprogram not used        do i = 1,ncol
! Subprogram not used           fxliq(i,k+1) = min( fxliq(i,k+1), mxsedfac * cldliq(i,k) * pdel(i,k) )
! Subprogram not used           fxice(i,k+1) = min( fxice(i,k+1), mxsedfac * cldice(i,k) * pdel(i,k) )
! Subprogram not used !!$        fxliq(i,k+1) = min( fxliq(i,k+1), cldliq(i,k) * pdel(i,k) + fxliq(i,k))
! Subprogram not used !!$        fxice(i,k+1) = min( fxice(i,k+1), cldice(i,k) * pdel(i,k) + fxice(i,k))
! Subprogram not used !!$        fxliq(i,k+1) = min( fxliq(i,k+1), cloud(i,k) * cldliq(i,k) * pdel(i,k) )
! Subprogram not used !!$        fxice(i,k+1) = min( fxice(i,k+1), cloud(i,k) * cldice(i,k) * pdel(i,k) )
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used ! Now calculate the tendencies assuming that condensate evaporates when
! Subprogram not used ! it falls into environment, and does not when it falls into cloud.
! Subprogram not used ! All flux out of the layer comes from the cloudy part.
! Subprogram not used ! Assume maximum overlap for stratiform clouds
! Subprogram not used !  if cloud above < cloud,  all water falls into cloud below
! Subprogram not used !  if cloud above > cloud,  water split between cloud  and environment
! Subprogram not used     do k = 1,pver
! Subprogram not used        cldab(:ncol) = 0._r8
! Subprogram not used        do i = 1,ncol
! Subprogram not used ! cloud overlap cloud factor
! Subprogram not used           cldovrl  = min( cloud(i,k) / (cldab(i)+.0001_r8), 1._r8 )
! Subprogram not used           cldab(i) = cloud(i,k)
! Subprogram not used ! evaporation into environment cause moistening and cooling
! Subprogram not used           evapliq = fxliq(i,k) * (1._r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s)
! Subprogram not used           evapice = fxice(i,k) * (1._r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s)
! Subprogram not used           wvtend(i,k) = evapliq + evapice                          ! evaporation into environment (kg/kg/s)
! Subprogram not used           htend (i,k) = -latvap*evapliq -(latvap+latice)*evapice   ! evaporation (W/kg)
! Subprogram not used ! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)
! Subprogram not used           liqtend(i,k)  = (fxliq(i,k)*cldovrl - fxliq(i,k+1)) / (dtime * pdel(i,k))
! Subprogram not used           icetend(i,k)  = (fxice(i,k)*cldovrl - fxice(i,k+1)) / (dtime * pdel(i,k))
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used ! convert flux out the bottom to mass units Pa -> kg/m2/s
! Subprogram not used     sfliq(:ncol) = fxliq(:ncol,pverp) / (dtime*gravit)
! Subprogram not used     sfice(:ncol) = fxice(:ncol,pverp) / (dtime*gravit)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine cld_sediment_tend

!===============================================================================
! Subprogram not used   subroutine getflx(ncol, xw, phi, vel, deltat, flux)
! Subprogram not used 
! Subprogram not used !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
! Subprogram not used !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
! Subprogram not used !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
! Subprogram not used !.........phi1......phi2.......phi3.....phi4.......phi5.......
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: ncol                      ! number of colums to process
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer k
! Subprogram not used 
! Subprogram not used     real (r8), intent(in) :: vel(pcols,pverp)
! Subprogram not used     real (r8) flux(pcols,pverp)
! Subprogram not used     real (r8) xw(pcols,pverp)
! Subprogram not used     real (r8) psi(pcols,pverp)
! Subprogram not used     real (r8), intent(in) :: phi(pcols,pverp-1)
! Subprogram not used     real (r8) fdot(pcols,pverp)
! Subprogram not used     real (r8) xx(pcols)
! Subprogram not used     real (r8) fxdot(pcols)
! Subprogram not used     real (r8) fxdd(pcols)
! Subprogram not used 
! Subprogram not used     real (r8) psistar(pcols)
! Subprogram not used     real (r8) deltat
! Subprogram not used 
! Subprogram not used     real (r8) xxk(pcols,pver)
! Subprogram not used 
! Subprogram not used     do i = 1,ncol
! Subprogram not used !        integral of phi
! Subprogram not used        psi(i,1) = 0._r8
! Subprogram not used !        fluxes at boundaries
! Subprogram not used        flux(i,1) = 0._r8
! Subprogram not used        flux(i,pverp) = 0._r8
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used !     integral function
! Subprogram not used     do k = 2,pverp
! Subprogram not used        do i = 1,ncol
! Subprogram not used           psi(i,k) = phi(i,k-1)*(xw(i,k)-xw(i,k-1)) + psi(i,k-1)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     calculate the derivatives for the interpolating polynomial
! Subprogram not used     call cfdotmc_pro (ncol, xw, psi, fdot)
! Subprogram not used 
! Subprogram not used !  NEW WAY
! Subprogram not used !     calculate fluxes at interior pts
! Subprogram not used     do k = 2,pver
! Subprogram not used        do i = 1,ncol
! Subprogram not used           xxk(i,k) = xw(i,k)-vel(i,k)*deltat
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used     do k = 2,pver
! Subprogram not used        call cfint2(ncol, xw, psi, fdot, xxk(1,k), fxdot, fxdd, psistar)
! Subprogram not used        do i = 1,ncol
! Subprogram not used           flux(i,k) = (psi(i,k)-psistar(i))
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine getflx



!##############################################################################

! Subprogram not used   subroutine cfint2 (ncol, x, f, fdot, xin, fxdot, fxdd, psistar)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! input
! Subprogram not used     integer ncol                      ! number of colums to process
! Subprogram not used 
! Subprogram not used     real (r8) x(pcols, pverp)
! Subprogram not used     real (r8) f(pcols, pverp)
! Subprogram not used     real (r8) fdot(pcols, pverp)
! Subprogram not used     real (r8) xin(pcols)
! Subprogram not used 
! Subprogram not used ! output
! Subprogram not used     real (r8) fxdot(pcols)
! Subprogram not used     real (r8) fxdd(pcols)
! Subprogram not used     real (r8) psistar(pcols)
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer k
! Subprogram not used     integer intz(pcols)
! Subprogram not used     real (r8) dx
! Subprogram not used     real (r8) s
! Subprogram not used     real (r8) c2
! Subprogram not used     real (r8) c3
! Subprogram not used     real (r8) xx
! Subprogram not used     real (r8) xinf
! Subprogram not used     real (r8) psi1, psi2, psi3, psim
! Subprogram not used     real (r8) cfint
! Subprogram not used     real (r8) cfnew
! Subprogram not used     real (r8) xins(pcols)
! Subprogram not used 
! Subprogram not used !     the minmod function 
! Subprogram not used     real (r8) a, b, c
! Subprogram not used     real (r8) minmod
! Subprogram not used     real (r8) medan
! Subprogram not used     logical found_error
! Subprogram not used 
! Subprogram not used     minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
! Subprogram not used     medan(a,b,c) = a + minmod(b-a,c-a)
! Subprogram not used 
! Subprogram not used     do i = 1,ncol
! Subprogram not used        xins(i) = medan(x(i,1), xin(i), x(i,pverp))
! Subprogram not used        intz(i) = 0
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used ! first find the interval 
! Subprogram not used     do k =  1,pverp-1
! Subprogram not used        do i = 1,ncol
! Subprogram not used           if ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).ge.0) then
! Subprogram not used              intz(i) = k
! Subprogram not used           endif
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     found_error=.false.
! Subprogram not used     do i = 1,ncol
! Subprogram not used        if (intz(i).eq.0._r8) found_error=.true.
! Subprogram not used     end do
! Subprogram not used     if(found_error) then
! Subprogram not used        do i = 1,ncol
! Subprogram not used           if (intz(i).eq.0._r8) then
! Subprogram not used              write(iulog,*) ' interval was not found for col i ', i
! Subprogram not used              call endrun('CFINT2')
! Subprogram not used           endif
! Subprogram not used        end do
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used ! now interpolate
! Subprogram not used     do i = 1,ncol
! Subprogram not used        k = intz(i)
! Subprogram not used        dx = (x(i,k+1)-x(i,k))
! Subprogram not used        s = (f(i,k+1)-f(i,k))/dx
! Subprogram not used        c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
! Subprogram not used        c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
! Subprogram not used        xx = (xins(i)-x(i,k))
! Subprogram not used        fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
! Subprogram not used        fxdd(i) = 6*c3*xx + 2*c2
! Subprogram not used        cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)
! Subprogram not used 
! Subprogram not used !        limit the interpolant
! Subprogram not used        psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
! Subprogram not used        if (k.eq.1) then
! Subprogram not used           psi2 = f(i,1)
! Subprogram not used        else
! Subprogram not used           psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
! Subprogram not used        endif
! Subprogram not used        if (k+1.eq.pverp) then
! Subprogram not used           psi3 = f(i,pverp)
! Subprogram not used        else
! Subprogram not used           psi3 = f(i,k+1) - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
! Subprogram not used        endif
! Subprogram not used        psim = medan(psi1, psi2, psi3)
! Subprogram not used        cfnew = medan(cfint, psi1, psim)
! Subprogram not used        if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
! Subprogram not used !     CHANGE THIS BACK LATER!!!
! Subprogram not used !     $        .gt..1) then
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     UNCOMMENT THIS LATER!!!
! Subprogram not used !            write(iulog,*) ' cfint2 limiting important ', cfint, cfnew
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        endif
! Subprogram not used        psistar(i) = cfnew
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine cfint2



!##############################################################################

! Subprogram not used   subroutine cfdotmc_pro (ncol, x, f, fdot)
! Subprogram not used 
! Subprogram not used !     prototype version; eventually replace with final SPITFIRE scheme
! Subprogram not used 
! Subprogram not used !     calculate the derivative for the interpolating polynomial
! Subprogram not used !     multi column version
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! input
! Subprogram not used     integer ncol                      ! number of colums to process
! Subprogram not used 
! Subprogram not used     real (r8) x(pcols, pverp)
! Subprogram not used     real (r8) f(pcols, pverp)
! Subprogram not used ! output
! Subprogram not used     real (r8) fdot(pcols, pverp)          ! derivative at nodes
! Subprogram not used 
! Subprogram not used ! assumed variable distribution
! Subprogram not used !     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
! Subprogram not used !     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
! Subprogram not used !     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
! Subprogram not used !     .........d2.......d3.......d4.......d5.........     2,pver points
! Subprogram not used !     .........s2.......s3.......s4.......s5.........     2,pver points
! Subprogram not used !     .............dh2......dh3......dh4.............     2,pver-1 points
! Subprogram not used !     .............eh2......eh3......eh4.............     2,pver-1 points
! Subprogram not used !     ..................e3.......e4..................     3,pver-1 points
! Subprogram not used !     .................ppl3......ppl4................     3,pver-1 points
! Subprogram not used !     .................ppr3......ppr4................     3,pver-1 points
! Subprogram not used !     .................t3........t4..................     3,pver-1 points
! Subprogram not used !     ................fdot3.....fdot4................     3,pver-1 points
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! work variables
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer k
! Subprogram not used 
! Subprogram not used     real (r8) a                    ! work var
! Subprogram not used     real (r8) b                    ! work var
! Subprogram not used     real (r8) c                    ! work var
! Subprogram not used     real (r8) s(pcols,pverp)             ! first divided differences at nodes
! Subprogram not used     real (r8) sh(pcols,pverp)            ! first divided differences between nodes
! Subprogram not used     real (r8) d(pcols,pverp)             ! second divided differences at nodes
! Subprogram not used     real (r8) dh(pcols,pverp)            ! second divided differences between nodes
! Subprogram not used     real (r8) e(pcols,pverp)             ! third divided differences at nodes
! Subprogram not used     real (r8) eh(pcols,pverp)            ! third divided differences between nodes
! Subprogram not used     real (r8) pp                   ! p prime
! Subprogram not used     real (r8) ppl(pcols,pverp)           ! p prime on left
! Subprogram not used     real (r8) ppr(pcols,pverp)           ! p prime on right
! Subprogram not used     real (r8) qpl
! Subprogram not used     real (r8) qpr
! Subprogram not used     real (r8) ttt
! Subprogram not used     real (r8) t
! Subprogram not used     real (r8) tmin
! Subprogram not used     real (r8) tmax
! Subprogram not used     real (r8) delxh(pcols,pverp)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     the minmod function 
! Subprogram not used     real (r8) minmod
! Subprogram not used     real (r8) medan
! Subprogram not used     minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
! Subprogram not used     medan(a,b,c) = a + minmod(b-a,c-a)
! Subprogram not used 
! Subprogram not used     do k = 1,pver
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !        first divided differences between nodes
! Subprogram not used        do i = 1, ncol
! Subprogram not used           delxh(i,k) = (x(i,k+1)-x(i,k))
! Subprogram not used           sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used !        first and second divided differences at nodes
! Subprogram not used        if (k.ge.2) then
! Subprogram not used           do i = 1,ncol
! Subprogram not used              d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
! Subprogram not used              s(i,k) = minmod(sh(i,k),sh(i,k-1))
! Subprogram not used           end do
! Subprogram not used        endif
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used !     second and third divided diffs between nodes
! Subprogram not used     do k = 2,pver-1
! Subprogram not used        do i = 1, ncol
! Subprogram not used           eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
! Subprogram not used           dh(i,k) = minmod(d(i,k),d(i,k+1))
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used !     treat the boundaries
! Subprogram not used     do i = 1,ncol
! Subprogram not used        e(i,2) = eh(i,2)
! Subprogram not used        e(i,pver) = eh(i,pver-1)
! Subprogram not used !        outside level
! Subprogram not used        fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)  &
! Subprogram not used             - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
! Subprogram not used        fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
! Subprogram not used        fdot(i,pverp) = sh(i,pver) + d(i,pver)*delxh(i,pver)  &
! Subprogram not used             + eh(i,pver-1)*delxh(i,pver)*(x(i,pverp)-x(i,pver-1))
! Subprogram not used        fdot(i,pverp) = minmod(fdot(i,pverp),3*sh(i,pver))
! Subprogram not used !        one in from boundary
! Subprogram not used        fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
! Subprogram not used        fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
! Subprogram not used        fdot(i,pver) = sh(i,pver) - d(i,pver)*delxh(i,pver)   &
! Subprogram not used             - eh(i,pver-1)*delxh(i,pver)*delxh(i,pver-1)
! Subprogram not used        fdot(i,pver) = minmod(fdot(i,pver),3*s(i,pver))
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do k = 3,pver-1
! Subprogram not used        do i = 1,ncol
! Subprogram not used           e(i,k) = minmod(eh(i,k),eh(i,k-1))
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do k = 3,pver-1
! Subprogram not used 
! Subprogram not used        do i = 1,ncol
! Subprogram not used 
! Subprogram not used !           p prime at k-0.5
! Subprogram not used           ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  
! Subprogram not used !           p prime at k+0.5
! Subprogram not used           ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)
! Subprogram not used 
! Subprogram not used           t = minmod(ppl(i,k),ppr(i,k))
! Subprogram not used 
! Subprogram not used !           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
! Subprogram not used           pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 
! Subprogram not used 
! Subprogram not used !           quartic estimate of fdot
! Subprogram not used           fdot(i,k) = pp                            &
! Subprogram not used                - delxh(i,k-1)*delxh(i,k)            &
! Subprogram not used                *(  eh(i,k-1)*(x(i,k+2)-x(i,k  ))    &
! Subprogram not used                + eh(i,k  )*(x(i,k  )-x(i,k-2))      &
! Subprogram not used                )/(x(i,k+2)-x(i,k-2))
! Subprogram not used 
! Subprogram not used !           now limit it
! Subprogram not used           qpl = sh(i,k-1)       &
! Subprogram not used                + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)*(x(i,k)-x(i,k-2)), &
! Subprogram not used                d(i,k)  -e(i,k)*delxh(i,k))
! Subprogram not used           qpr = sh(i,k)         &
! Subprogram not used                + delxh(i,k  )*minmod(d(i,k)  +e(i,k)*delxh(i,k-1),        &
! Subprogram not used                d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))
! Subprogram not used 
! Subprogram not used           fdot(i,k) = medan(fdot(i,k), qpl, qpr)
! Subprogram not used 
! Subprogram not used           ttt = minmod(qpl, qpr)
! Subprogram not used           tmin = min(0._r8,3*s(i,k),1.5_r8*t,ttt)
! Subprogram not used           tmax = max(0._r8,3*s(i,k),1.5_r8*t,ttt)
! Subprogram not used 
! Subprogram not used           fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine cfdotmc_pro
end module pkg_cld_sediment
