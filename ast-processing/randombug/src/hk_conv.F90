module hk_conv
!
! Moist convection. Primarily data used by both Zhang-McFarlane convection
! and Hack shallow convective schemes.
!
! $Id$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc
   use abortutils,   only: endrun
   implicit none

   private
   save
!
! Public interfaces
!
   public mfinti   !  Initialization of data for moist convection
   public cmfmca   !  Hack shallow convection
   public hkconv_readnl ! read hkconv_nl namelist

!
! Private data used for Hack shallow convection
!
   real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Namelist variables
   real(r8) :: hkconv_c0 = unset_r8    
   real(r8) :: hkconv_cmftau = unset_r8 

   real(r8) :: hlat        ! latent heat of vaporization
   real(r8) :: c0          ! rain water autoconversion coefficient set from namelist input hkconv_c0
   real(r8) :: betamn      ! minimum overshoot parameter
   real(r8) :: rhlat       ! reciprocal of hlat
   real(r8) :: rcp         ! reciprocal of cp
   real(r8) :: cmftau      ! characteristic adjustment time scale set from namelist input hkconv_cmftau
   real(r8) :: rhoh2o      ! density of liquid water (STP)
   real(r8) :: dzmin       ! minimum convective depth for precipitation
   real(r8) :: tiny        ! arbitrary small num used in transport estimates
   real(r8) :: eps         ! convergence criteria (machine dependent)
   real(r8) :: tpmax       ! maximum acceptable t perturbation (degrees C)
   real(r8) :: shpmax      ! maximum acceptable q perturbation (g/g)           

   integer :: iloc         ! longitude location for diagnostics
   integer :: jloc         ! latitude  location for diagnostics
   integer :: nsloc        ! nstep for which to produce diagnostics
!
   logical :: rlxclm       ! logical to relax column versus cloud triplet

   real(r8) cp          ! specific heat of dry air
   real(r8) grav        ! gravitational constant       
   real(r8) rgrav       ! reciprocal of grav
   real(r8) rgas        ! gas constant for dry air
   integer limcnv          ! top interface level limit for convection




contains
subroutine hkconv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'hkconv_readnl'

   namelist /hkconv_nl/ hkconv_cmftau, hkconv_c0
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'hkconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, hkconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      cmftau = hkconv_cmftau
      c0     = hkconv_c0

   end if

   ! Broadcast namelist variables
   call mpibcast(cmftau,            1, mpir8,  0, mpicom)
   call mpibcast(c0,                1, mpir8,  0, mpicom)

end subroutine hkconv_readnl

!================================================================================================

! Subprogram not used subroutine mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhowtr,limcnv_in )
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Initialize moist convective mass flux procedure common block, cmfmca
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! <Describe the algorithm(s) used in the routine.> 
! Subprogram not used ! <Also include any applicable external references.> 
! Subprogram not used ! 
! Subprogram not used ! Author: J. Hack
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use dycore, only: dycore_is, get_resolution
! Subprogram not used    use spmd_utils, only: masterproc
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: rair              ! gas constant for dry air
! Subprogram not used    real(r8), intent(in) :: cpair             ! specific heat of dry air
! Subprogram not used    real(r8), intent(in) :: gravit            ! acceleration due to gravity
! Subprogram not used    real(r8), intent(in) :: latvap            ! latent heat of vaporization
! Subprogram not used    real(r8), intent(in) :: rhowtr            ! density of liquid water (STP)
! Subprogram not used    integer,  intent(in) :: limcnv_in         ! top interface level limit for convection
! Subprogram not used 
! Subprogram not used    ! local variables
! Subprogram not used    character(len=32)    :: hgrid             ! horizontal grid specifier
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Initialize physical constants for moist convective mass flux procedure
! Subprogram not used !
! Subprogram not used    cp     = cpair         ! specific heat of dry air
! Subprogram not used    hlat   = latvap        ! latent heat of vaporization
! Subprogram not used    grav   = gravit        ! gravitational constant
! Subprogram not used    rgas   = rair          ! gas constant for dry air
! Subprogram not used    rhoh2o = rhowtr        ! density of liquid water (STP)
! Subprogram not used 
! Subprogram not used    limcnv = limcnv_in
! Subprogram not used 
! Subprogram not used    ! Initialize free parameters for moist convective mass flux procedure
! Subprogram not used    ! cmftau - characteristic adjustment time scale
! Subprogram not used    ! c0     - rain water autoconversion coeff (1/m)
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       write(iulog,*) 'tuning parameters hk_conv: cmftau',cmftau
! Subprogram not used       write(iulog,*) 'tuning parameters hk_conv: c0',c0
! Subprogram not used    endif
! Subprogram not used    dzmin  = 0.0_r8           ! minimum cloud depth to precipitate (m)
! Subprogram not used    betamn = 0.10_r8          ! minimum overshoot parameter
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    tpmax  = 1.50_r8          ! maximum acceptable t perturbation (deg C)
! Subprogram not used    shpmax = 1.50e-3_r8       ! maximum acceptable q perturbation (g/g)
! Subprogram not used    rlxclm = .true.        ! logical variable to specify that relaxation
! Subprogram not used !                                time scale should applied to column as
! Subprogram not used !                                opposed to triplets individually
! Subprogram not used !
! Subprogram not used ! Initialize miscellaneous (frequently used) constants
! Subprogram not used !
! Subprogram not used    rhlat  = 1.0_r8/hlat      ! reciprocal latent heat of vaporization
! Subprogram not used    rcp    = 1.0_r8/cp        ! reciprocal specific heat of dry air
! Subprogram not used    rgrav  = 1.0_r8/grav      ! reciprocal gravitational constant
! Subprogram not used !
! Subprogram not used ! Initialize diagnostic location information for moist convection scheme
! Subprogram not used !
! Subprogram not used    iloc   = 1             ! longitude point for diagnostic info
! Subprogram not used    jloc   = 1             ! latitude  point for diagnostic info
! Subprogram not used    nsloc  = 1             ! nstep value at which to begin diagnostics
! Subprogram not used !
! Subprogram not used ! Initialize other miscellaneous parameters
! Subprogram not used !
! Subprogram not used    tiny   = 1.0e-36_r8       ! arbitrary small number (scalar transport)
! Subprogram not used    eps    = 1.0e-13_r8       ! convergence criteria (machine dependent)
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used end subroutine mfinti

! Subprogram not used subroutine cmfmca(lchnk   ,ncol    , &
! Subprogram not used                   nstep   ,ztodt     ,pmid    ,pdel    , &
! Subprogram not used                   rpdel   ,zm      ,tpert   ,qpert   ,phis    , &
! Subprogram not used                   pblh    ,t       ,q       ,cmfdt   ,dq      , &
! Subprogram not used                   cmfmc   ,cmfdqr  ,cmfsl   ,cmflq   ,precc   , &
! Subprogram not used                   qc      ,cnt     ,cnb     ,icwmr   ,rliq    , & 
! Subprogram not used                   pmiddry ,pdeldry ,rpdeldry)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Moist convective mass flux procedure:
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! If stratification is unstable to nonentraining parcel ascent,
! Subprogram not used ! complete an adjustment making successive use of a simple cloud model
! Subprogram not used ! consisting of three layers (sometimes referred to as a triplet)
! Subprogram not used !
! Subprogram not used ! Code generalized to allow specification of parcel ("updraft")
! Subprogram not used ! properties, as well as convective transport of an arbitrary
! Subprogram not used ! number of passive constituents (see q array).  The code
! Subprogram not used ! is written so the water vapor field is passed independently
! Subprogram not used ! in the calling list from the block of other transported
! Subprogram not used ! constituents, even though as currently designed, it is the
! Subprogram not used ! first component in the constituents field.
! Subprogram not used ! 
! Subprogram not used ! Author: J. Hack
! Subprogram not used !
! Subprogram not used ! BAB: changed code to report tendencies in cmfdt and dq, instead of
! Subprogram not used ! updating profiles. Cmfdq contains water only, made it a local variable
! Subprogram not used ! made dq (all constituents) the argument.
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !#######################################################################
! Subprogram not used !#                                                                     #
! Subprogram not used !# Debugging blocks are marked this way for easy identification        #
! Subprogram not used !#                                                                     #
! Subprogram not used !#######################################################################
! Subprogram not used    use constituents,  only: pcnst
! Subprogram not used    use constituents,    only: cnst_get_type_byind
! Subprogram not used    use ppgrid,    only: pcols, pver, pverp
! Subprogram not used    use phys_grid, only: get_lat_all_p, get_lon_all_p
! Subprogram not used    use wv_saturation, only: qsat
! Subprogram not used 
! Subprogram not used    real(r8) ssfac               ! supersaturation bound (detrained air)
! Subprogram not used    parameter (ssfac = 1.001_r8)
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: lchnk                ! chunk identifier
! Subprogram not used    integer, intent(in) :: ncol                 ! number of atmospheric columns
! Subprogram not used    integer, intent(in) :: nstep                ! current time step index
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: ztodt               ! 2 delta-t (seconds)
! Subprogram not used    real(r8), intent(in) :: pmid(pcols,pver)    ! pressure
! Subprogram not used    real(r8), intent(in) :: pdel(pcols,pver)    ! delta-p
! Subprogram not used    real(r8), intent(in) :: pmiddry(pcols,pver)    ! pressure
! Subprogram not used    real(r8), intent(in) :: pdeldry(pcols,pver)    ! delta-p
! Subprogram not used    real(r8), intent(in) :: rpdel(pcols,pver)   ! 1./pdel
! Subprogram not used    real(r8), intent(in) :: rpdeldry(pcols,pver)   ! 1./pdel
! Subprogram not used    real(r8), intent(in) :: zm(pcols,pver)      ! height abv sfc at midpoints
! Subprogram not used    real(r8), intent(in) :: tpert(pcols)        ! PBL perturbation theta
! Subprogram not used    real(r8), intent(in) :: qpert(pcols,pcnst)  ! PBL perturbation specific humidity
! Subprogram not used    real(r8), intent(in) :: phis(pcols)         ! surface geopotential
! Subprogram not used    real(r8), intent(in) :: pblh(pcols)         ! PBL height (provided by PBL routine)
! Subprogram not used    real(r8), intent(in) :: t(pcols,pver)       ! temperature (t bar)
! Subprogram not used    real(r8), intent(in) :: q(pcols,pver,pcnst) ! specific humidity (sh bar)
! Subprogram not used !
! Subprogram not used ! Output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: cmfdt(pcols,pver)   ! dt/dt due to moist convection
! Subprogram not used    real(r8), intent(out) :: cmfmc(pcols,pverp)  ! moist convection cloud mass flux
! Subprogram not used    real(r8), intent(out) :: cmfdqr(pcols,pver)  ! dq/dt due to convective rainout
! Subprogram not used    real(r8), intent(out) :: cmfsl(pcols,pver )  ! convective lw static energy flux
! Subprogram not used    real(r8), intent(out) :: cmflq(pcols,pver )  ! convective total water flux
! Subprogram not used    real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
! Subprogram not used ! JJH mod to explicitly export cloud water
! Subprogram not used    real(r8), intent(out) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
! Subprogram not used    real(r8), intent(out) :: cnt(pcols)          ! top level of convective activity
! Subprogram not used    real(r8), intent(out) :: cnb(pcols)          ! bottom level of convective activity
! Subprogram not used    real(r8), intent(out) :: dq(pcols,pver,pcnst) ! constituent tendencies
! Subprogram not used    real(r8), intent(out) :: icwmr(pcols,pver)
! Subprogram not used    real(r8), intent(out) :: rliq(pcols) 
! Subprogram not used !
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    real(r8) pm(pcols,pver)    ! pressure
! Subprogram not used    real(r8) pd(pcols,pver)    ! delta-p
! Subprogram not used    real(r8) rpd(pcols,pver)   ! 1./pdel
! Subprogram not used 
! Subprogram not used    real(r8) cmfdq(pcols,pver)   ! dq/dt due to moist convection
! Subprogram not used    real(r8) gam(pcols,pver)     ! 1/cp (d(qsat)/dT)
! Subprogram not used    real(r8) sb(pcols,pver)      ! dry static energy (s bar)
! Subprogram not used    real(r8) hb(pcols,pver)      ! moist static energy (h bar)
! Subprogram not used    real(r8) shbs(pcols,pver)    ! sat. specific humidity (sh bar star)
! Subprogram not used    real(r8) hbs(pcols,pver)     ! sat. moist static energy (h bar star)
! Subprogram not used    real(r8) shbh(pcols,pverp)   ! specific humidity on interfaces
! Subprogram not used    real(r8) sbh(pcols,pverp)    ! s bar on interfaces
! Subprogram not used    real(r8) hbh(pcols,pverp)    ! h bar on interfaces
! Subprogram not used    real(r8) cmrh(pcols,pverp)   ! interface constituent mixing ratio
! Subprogram not used    real(r8) prec(pcols)         ! instantaneous total precipitation
! Subprogram not used    real(r8) dzcld(pcols)        ! depth of convective layer (m)
! Subprogram not used    real(r8) beta(pcols)         ! overshoot parameter (fraction)
! Subprogram not used    real(r8) betamx(pcols)       ! local maximum on overshoot
! Subprogram not used    real(r8) eta(pcols)          ! convective mass flux (kg/m^2 s)
! Subprogram not used    real(r8) etagdt(pcols)       ! eta*grav*dt
! Subprogram not used    real(r8) cldwtr(pcols)       ! cloud water (mass)
! Subprogram not used    real(r8) rnwtr(pcols)        ! rain water  (mass)
! Subprogram not used !  JJH extension to facilitate export of cloud liquid water
! Subprogram not used    real(r8) totcond(pcols)	! total condensate; mix of precip and cloud water (mass)
! Subprogram not used    real(r8) sc  (pcols)         ! dry static energy   ("in-cloud")
! Subprogram not used    real(r8) shc (pcols)         ! specific humidity   ("in-cloud")
! Subprogram not used    real(r8) hc  (pcols)         ! moist static energy ("in-cloud")
! Subprogram not used    real(r8) cmrc(pcols)         ! constituent mix rat ("in-cloud")
! Subprogram not used    real(r8) dq1(pcols)          ! shb  convective change (lower lvl)
! Subprogram not used    real(r8) dq2(pcols)          ! shb  convective change (mid level)
! Subprogram not used    real(r8) dq3(pcols)          ! shb  convective change (upper lvl)
! Subprogram not used    real(r8) ds1(pcols)          ! sb   convective change (lower lvl)
! Subprogram not used    real(r8) ds2(pcols)          ! sb   convective change (mid level)
! Subprogram not used    real(r8) ds3(pcols)          ! sb   convective change (upper lvl)
! Subprogram not used    real(r8) dcmr1(pcols)        ! q convective change (lower lvl)
! Subprogram not used    real(r8) dcmr2(pcols)        ! q convective change (mid level)
! Subprogram not used    real(r8) dcmr3(pcols)        ! q convective change (upper lvl)
! Subprogram not used    real(r8) estemp(pcols,pver)  ! saturation vapor pressure (scratch)
! Subprogram not used    real(r8) vtemp1(2*pcols)     ! intermediate scratch vector
! Subprogram not used    real(r8) vtemp2(2*pcols)     ! intermediate scratch vector
! Subprogram not used    real(r8) vtemp3(2*pcols)     ! intermediate scratch vector
! Subprogram not used    real(r8) vtemp4(2*pcols)     ! intermediate scratch vector
! Subprogram not used    real(r8) vtemp5(2*pcols)     ! intermediate scratch vector
! Subprogram not used    integer indx1(pcols)     ! longitude indices for condition true
! Subprogram not used    logical etagt0           ! true if eta > 0.0
! Subprogram not used    real(r8) sh1                 ! dummy arg in qhalf statement func.
! Subprogram not used    real(r8) sh2                 ! dummy arg in qhalf statement func.
! Subprogram not used    real(r8) shbs1               ! dummy arg in qhalf statement func.
! Subprogram not used    real(r8) shbs2               ! dummy arg in qhalf statement func.
! Subprogram not used    real(r8) cats                ! modified characteristic adj. time
! Subprogram not used    real(r8) rtdt                ! 1./ztodt
! Subprogram not used    real(r8) qprime              ! modified specific humidity pert.
! Subprogram not used    real(r8) tprime              ! modified thermal perturbation
! Subprogram not used    real(r8) pblhgt              ! bounded pbl height (max[pblh,1m])
! Subprogram not used    real(r8) fac1                ! intermediate scratch variable
! Subprogram not used    real(r8) shprme              ! intermediate specific humidity pert.
! Subprogram not used    real(r8) qsattp              ! sat mix rat for thermally pert PBL parcels
! Subprogram not used    real(r8) dz                  ! local layer depth
! Subprogram not used    real(r8) temp1               ! intermediate scratch variable
! Subprogram not used    real(r8) b1                  ! bouyancy measure in detrainment lvl
! Subprogram not used    real(r8) b2                  ! bouyancy measure in condensation lvl
! Subprogram not used    real(r8) temp2               ! intermediate scratch variable
! Subprogram not used    real(r8) temp3               ! intermediate scratch variable
! Subprogram not used    real(r8) g                   ! bounded vertical gradient of hb
! Subprogram not used    real(r8) tmass               ! total mass available for convective exch
! Subprogram not used    real(r8) denom               ! intermediate scratch variable
! Subprogram not used    real(r8) qtest1              ! used in negative q test (middle lvl)
! Subprogram not used    real(r8) qtest2              ! used in negative q test (lower lvl)
! Subprogram not used    real(r8) fslkp               ! flux lw static energy (bot interface)
! Subprogram not used    real(r8) fslkm               ! flux lw static energy (top interface)
! Subprogram not used    real(r8) fqlkp               ! flux total water (bottom interface)
! Subprogram not used    real(r8) fqlkm               ! flux total water (top interface)
! Subprogram not used    real(r8) botflx              ! bottom constituent mixing ratio flux
! Subprogram not used    real(r8) topflx              ! top constituent mixing ratio flux
! Subprogram not used    real(r8) efac1               ! ratio q to convectively induced chg (btm lvl)
! Subprogram not used    real(r8) efac2               ! ratio q to convectively induced chg (mid lvl)
! Subprogram not used    real(r8) efac3               ! ratio q to convectively induced chg (top lvl)
! Subprogram not used    real(r8) tb(pcols,pver)      ! working storage for temp (t bar)
! Subprogram not used    real(r8) shb(pcols,pver)     ! working storage for spec hum (sh bar)
! Subprogram not used    real(r8) adjfac              ! adjustment factor (relaxation related)
! Subprogram not used    real(r8) rktp
! Subprogram not used    real(r8) rk
! Subprogram not used    integer i,k              ! longitude, level indices
! Subprogram not used    integer ii               ! index on "gathered" vectors
! Subprogram not used    integer len1             ! vector length of "gathered" vectors
! Subprogram not used    integer m                ! constituent index
! Subprogram not used    integer ktp              ! tmp indx used to track top of convective layer
! Subprogram not used !
! Subprogram not used !---------------------------Statement functions-------------------------
! Subprogram not used !
! Subprogram not used    real(r8) qhalf
! Subprogram not used    qhalf(sh1,sh2,shbs1,shbs2) = min(max(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !** BAB initialize output tendencies here
! Subprogram not used !       copy q to dq; use dq below for passive tracer transport
! Subprogram not used    cmfdt(:ncol,:)  = 0._r8
! Subprogram not used    cmfdq(:ncol,:)  = 0._r8
! Subprogram not used    dq(:ncol,:,2:)  = q(:ncol,:,2:)
! Subprogram not used    cmfmc(:ncol,:)  = 0._r8
! Subprogram not used    cmfdqr(:ncol,:) = 0._r8
! Subprogram not used    cmfsl(:ncol,:)  = 0._r8
! Subprogram not used    cmflq(:ncol,:)  = 0._r8
! Subprogram not used    qc(:ncol,:)     = 0._r8
! Subprogram not used    rliq(:ncol)     = 0._r8
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Ensure that characteristic adjustment time scale (cmftau) assumed
! Subprogram not used ! in estimate of eta isn't smaller than model time scale (ztodt)
! Subprogram not used ! The time over which the convection is assumed to act (the adjustment
! Subprogram not used ! time scale) can be applied with each application of the three-level
! Subprogram not used ! cloud model, or applied to the column tendencies after a "hard"
! Subprogram not used ! adjustment (i.e., on a 2-delta t time scale) is evaluated
! Subprogram not used !
! Subprogram not used    if (rlxclm) then
! Subprogram not used       cats   = ztodt             ! relaxation applied to column
! Subprogram not used       adjfac = ztodt/(max(ztodt,cmftau))
! Subprogram not used    else
! Subprogram not used       cats   = max(ztodt,cmftau) ! relaxation applied to triplet
! Subprogram not used       adjfac = 1.0_r8
! Subprogram not used    endif
! Subprogram not used    rtdt = 1.0_r8/ztodt
! Subprogram not used !
! Subprogram not used ! Move temperature and moisture into working storage
! Subprogram not used !
! Subprogram not used    do k=limcnv,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          tb (i,k) = t(i,k)
! Subprogram not used          shb(i,k) = q(i,k,1)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used    do k=1,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          icwmr(i,k) = 0._r8
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Compute sb,hb,shbs,hbs
! Subprogram not used !
! Subprogram not used    call qsat(tb(:ncol,limcnv:pver), pmid(:ncol,limcnv:pver), &
! Subprogram not used         estemp(:ncol,limcnv:pver), shbs(:ncol,limcnv:pver), &
! Subprogram not used         gam=gam(:ncol,limcnv:pver))
! Subprogram not used !
! Subprogram not used    do k=limcnv,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          sb (i,k) = cp*tb(i,k) + zm(i,k)*grav + phis(i)
! Subprogram not used          hb (i,k) = sb(i,k) + hlat*shb(i,k)
! Subprogram not used          hbs(i,k) = sb(i,k) + hlat*shbs(i,k)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Compute sbh, shbh
! Subprogram not used !
! Subprogram not used    do k=limcnv+1,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          sbh (i,k) = 0.5_r8*(sb(i,k-1) + sb(i,k))
! Subprogram not used          shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
! Subprogram not used          hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Specify properties at top of model (not used, but filling anyway)
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       sbh (i,limcnv) = sb(i,limcnv)
! Subprogram not used       shbh(i,limcnv) = shb(i,limcnv)
! Subprogram not used       hbh (i,limcnv) = hb(i,limcnv)
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Zero vertically independent control, tendency & diagnostic arrays
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       prec(i)  = 0.0_r8
! Subprogram not used       dzcld(i) = 0.0_r8
! Subprogram not used       cnb(i)   = 0.0_r8
! Subprogram not used       cnt(i)   = real(pver+1,r8)
! Subprogram not used    end do
! Subprogram not used !#                                                                     #
! Subprogram not used !#                                                                     #
! Subprogram not used !#######################################################################
! Subprogram not used !
! Subprogram not used ! Begin moist convective mass flux adjustment procedure.
! Subprogram not used ! Formalism ensures that negative cloud liquid water can never occur
! Subprogram not used !
! Subprogram not used    do 70 k=pver-1,limcnv+1,-1
! Subprogram not used       do 10 i=1,ncol
! Subprogram not used          etagdt(i) = 0.0_r8
! Subprogram not used          eta   (i) = 0.0_r8
! Subprogram not used          beta  (i) = 0.0_r8
! Subprogram not used          ds1   (i) = 0.0_r8
! Subprogram not used          ds2   (i) = 0.0_r8
! Subprogram not used          ds3   (i) = 0.0_r8
! Subprogram not used          dq1   (i) = 0.0_r8
! Subprogram not used          dq2   (i) = 0.0_r8
! Subprogram not used          dq3   (i) = 0.0_r8
! Subprogram not used !
! Subprogram not used ! Specification of "cloud base" conditions
! Subprogram not used !
! Subprogram not used          qprime    = 0.0_r8
! Subprogram not used          tprime    = 0.0_r8
! Subprogram not used !
! Subprogram not used ! Assign tprime within the PBL to be proportional to the quantity
! Subprogram not used ! tpert (which will be bounded by tpmax), passed to this routine by
! Subprogram not used ! the PBL routine.  Don't allow perturbation to produce a dry
! Subprogram not used ! adiabatically unstable parcel.  Assign qprime within the PBL to be
! Subprogram not used ! an appropriately modified value of the quantity qpert (which will be
! Subprogram not used ! bounded by shpmax) passed to this routine by the PBL routine.  The
! Subprogram not used ! quantity qprime should be less than the local saturation value
! Subprogram not used ! (qsattp=qsat[t+tprime,p]).  In both cases, tpert and qpert are
! Subprogram not used ! linearly reduced toward zero as the PBL top is approached.
! Subprogram not used !
! Subprogram not used          pblhgt = max(pblh(i),1.0_r8)
! Subprogram not used          if ( (zm(i,k+1) <= pblhgt) .and. dzcld(i) == 0.0_r8 ) then
! Subprogram not used             fac1   = max(0.0_r8,1.0_r8-zm(i,k+1)/pblhgt)
! Subprogram not used             tprime = min(tpert(i),tpmax)*fac1
! Subprogram not used             qsattp = shbs(i,k+1) + cp*rhlat*gam(i,k+1)*tprime
! Subprogram not used             shprme = min(min(qpert(i,1),shpmax)*fac1,max(qsattp-shb(i,k+1),0.0_r8))
! Subprogram not used             qprime = max(qprime,shprme)
! Subprogram not used          else
! Subprogram not used             tprime = 0.0_r8
! Subprogram not used             qprime = 0.0_r8
! Subprogram not used          end if
! Subprogram not used !
! Subprogram not used ! Specify "updraft" (in-cloud) thermodynamic properties
! Subprogram not used !
! Subprogram not used          sc (i)    = sb (i,k+1) + cp*tprime
! Subprogram not used          shc(i)    = shb(i,k+1) + qprime
! Subprogram not used          hc (i)    = sc (i    ) + hlat*shc(i)
! Subprogram not used          vtemp4(i) = hc(i) - hbs(i,k)
! Subprogram not used          dz        = pdel(i,k)*rgas*tb(i,k)*rgrav/pmid(i,k)
! Subprogram not used          if (vtemp4(i) > 0.0_r8) then
! Subprogram not used             dzcld(i) = dzcld(i) + dz
! Subprogram not used          else
! Subprogram not used             dzcld(i) = 0.0_r8
! Subprogram not used          end if
! Subprogram not used 10       continue
! Subprogram not used !
! Subprogram not used ! Check on moist convective instability
! Subprogram not used ! Build index vector of points where instability exists
! Subprogram not used !
! Subprogram not used          len1 = 0
! Subprogram not used          do i=1,ncol
! Subprogram not used             if (vtemp4(i) > 0.0_r8) then
! Subprogram not used                len1 = len1 + 1
! Subprogram not used                indx1(len1) = i
! Subprogram not used             end if
! Subprogram not used          end do
! Subprogram not used          if (len1 <= 0) go to 70
! Subprogram not used !
! Subprogram not used ! Current level just below top level => no overshoot
! Subprogram not used !
! Subprogram not used          if (k <= limcnv+1) then
! Subprogram not used             do ii=1,len1
! Subprogram not used                i = indx1(ii)
! Subprogram not used                temp1     = vtemp4(i)/(1.0_r8 + gam(i,k))
! Subprogram not used                cldwtr(i) = max(0.0_r8,(sb(i,k) - sc(i) + temp1))
! Subprogram not used                beta(i)   = 0.0_r8
! Subprogram not used                vtemp3(i) = (1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k))
! Subprogram not used             end do
! Subprogram not used          else
! Subprogram not used !
! Subprogram not used ! First guess at overshoot parameter using crude buoyancy closure
! Subprogram not used ! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
! Subprogram not used ! If pre-existing supersaturation in detrainment layer, beta=0
! Subprogram not used ! cldwtr is temporarily equal to hlat*l (l=> liquid water)
! Subprogram not used !
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used             do ii=1,len1
! Subprogram not used                i = indx1(ii)
! Subprogram not used                temp1     = vtemp4(i)/(1.0_r8 + gam(i,k))
! Subprogram not used                cldwtr(i) = max(0.0_r8,(sb(i,k)-sc(i)+temp1))
! Subprogram not used                betamx(i) = 1.0_r8 - c0*max(0.0_r8,(dzcld(i)-dzmin))
! Subprogram not used                b1        = (hc(i) - hbs(i,k-1))*pdel(i,k-1)
! Subprogram not used                b2        = (hc(i) - hbs(i,k  ))*pdel(i,k  )
! Subprogram not used                beta(i)   = max(betamn,min(betamx(i), 1.0_r8 + b1/b2))
! Subprogram not used                if (hbs(i,k-1) <= hb(i,k-1)) beta(i) = 0.0_r8
! Subprogram not used !
! Subprogram not used ! Bound maximum beta to ensure physically realistic solutions
! Subprogram not used !
! Subprogram not used ! First check constrains beta so that eta remains positive
! Subprogram not used ! (assuming that eta is already positive for beta equal zero)
! Subprogram not used !
! Subprogram not used                vtemp1(i) = -(hbh(i,k+1) - hc(i))*pdel(i,k)*rpdel(i,k+1)+ &
! Subprogram not used                            (1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i))
! Subprogram not used                vtemp2(i) = (1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k))
! Subprogram not used                vtemp3(i) = vtemp2(i)
! Subprogram not used                if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0._r8) then
! Subprogram not used                   betamx(i) = 0.99_r8*(vtemp1(i)/vtemp2(i))
! Subprogram not used                   beta(i)   = max(0.0_r8,min(betamx(i),beta(i)))
! Subprogram not used                end if
! Subprogram not used             end do
! Subprogram not used !
! Subprogram not used ! Second check involves supersaturation of "detrainment layer"
! Subprogram not used ! small amount of supersaturation acceptable (by ssfac factor)
! Subprogram not used !
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used             do ii=1,len1
! Subprogram not used                i = indx1(ii)
! Subprogram not used                if (hb(i,k-1) < hbs(i,k-1)) then
! Subprogram not used                   vtemp1(i) = vtemp1(i)*rpdel(i,k)
! Subprogram not used                   temp2 = gam(i,k-1)*(sbh(i,k) - sc(i) + cldwtr(i)) -  &
! Subprogram not used                           hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
! Subprogram not used                   temp3 = vtemp3(i)*rpdel(i,k)
! Subprogram not used                   vtemp2(i) = (ztodt/cats)*(hc(i) - hbs(i,k))*temp2/ &
! Subprogram not used                               (pdel(i,k-1)*(hbs(i,k-1) - hb(i,k-1))) + temp3
! Subprogram not used                   if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0._r8) then
! Subprogram not used                      betamx(i) = ssfac*(vtemp1(i)/vtemp2(i))
! Subprogram not used                      beta(i)   = max(0.0_r8,min(betamx(i),beta(i)))
! Subprogram not used                   end if
! Subprogram not used                else
! Subprogram not used                   beta(i) = 0.0_r8
! Subprogram not used                end if
! Subprogram not used             end do
! Subprogram not used !
! Subprogram not used ! Third check to avoid introducing 2 delta x thermodynamic
! Subprogram not used ! noise in the vertical ... constrain adjusted h (or theta e)
! Subprogram not used ! so that the adjustment doesn't contribute to "kinks" in h
! Subprogram not used !
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used             do ii=1,len1
! Subprogram not used                i = indx1(ii)
! Subprogram not used                g = min(0.0_r8,hb(i,k) - hb(i,k-1))
! Subprogram not used                temp1 = (hb(i,k) - hb(i,k-1) - g)*(cats/ztodt)/(hc(i) - hbs(i,k))
! Subprogram not used                vtemp1(i) = temp1*vtemp1(i) + (hc(i) - hbh(i,k+1))*rpdel(i,k)
! Subprogram not used                vtemp2(i) = temp1*vtemp3(i)*rpdel(i,k) + (hc(i) - hbh(i,k) - cldwtr(i))* &
! Subprogram not used                            (rpdel(i,k) + rpdel(i,k+1))
! Subprogram not used                if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0._r8) then
! Subprogram not used                   if (vtemp2(i) /= 0.0_r8) then
! Subprogram not used                      betamx(i) = vtemp1(i)/vtemp2(i)
! Subprogram not used                   else
! Subprogram not used                      betamx(i) = 0.0_r8
! Subprogram not used                   end if
! Subprogram not used                   beta(i) = max(0.0_r8,min(betamx(i),beta(i)))
! Subprogram not used                end if
! Subprogram not used             end do
! Subprogram not used          end if
! Subprogram not used !
! Subprogram not used ! Calculate mass flux required for stabilization.
! Subprogram not used !
! Subprogram not used ! Ensure that the convective mass flux, eta, is positive by
! Subprogram not used ! setting negative values of eta to zero..
! Subprogram not used ! Ensure that estimated mass flux cannot move more than the
! Subprogram not used ! minimum of total mass contained in either layer k or layer k+1.
! Subprogram not used ! Also test for other pathological cases that result in non-
! Subprogram not used ! physical states and adjust eta accordingly.
! Subprogram not used !
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used          do ii=1,len1
! Subprogram not used             i = indx1(ii)
! Subprogram not used             beta(i) = max(0.0_r8,beta(i))
! Subprogram not used             temp1 = hc(i) - hbs(i,k)
! Subprogram not used             temp2 = ((1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i)) - &
! Subprogram not used                       beta(i)*vtemp3(i))*rpdel(i,k) - (hbh(i,k+1) - hc(i))*rpdel(i,k+1)
! Subprogram not used             eta(i) = temp1/(temp2*grav*cats)
! Subprogram not used             tmass = min(pdel(i,k),pdel(i,k+1))*rgrav
! Subprogram not used             if (eta(i) > tmass*rtdt .or. eta(i) <= 0.0_r8) eta(i) = 0.0_r8
! Subprogram not used !
! Subprogram not used ! Check on negative q in top layer (bound beta)
! Subprogram not used !
! Subprogram not used             if (shc(i)-shbh(i,k) < 0.0_r8 .and. beta(i)*eta(i) /= 0.0_r8) then
! Subprogram not used                denom = eta(i)*grav*ztodt*(shc(i) - shbh(i,k))*rpdel(i,k-1)
! Subprogram not used                beta(i) = max(0.0_r8,min(-0.999_r8*shb(i,k-1)/denom,beta(i)))
! Subprogram not used             end if
! Subprogram not used !
! Subprogram not used ! Check on negative q in middle layer (zero eta)
! Subprogram not used !
! Subprogram not used             qtest1 = shb(i,k) + eta(i)*grav*ztodt*((shc(i) - shbh(i,k+1)) - &
! Subprogram not used                      (1.0_r8 - beta(i))*cldwtr(i)*rhlat - beta(i)*(shc(i) - shbh(i,k)))* &
! Subprogram not used 	             rpdel(i,k)
! Subprogram not used             if (qtest1 <= 0.0_r8) eta(i) = 0.0_r8
! Subprogram not used !
! Subprogram not used ! Check on negative q in lower layer (bound eta)
! Subprogram not used !
! Subprogram not used             fac1 = -(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
! Subprogram not used             qtest2 = shb(i,k+1) - eta(i)*grav*ztodt*fac1
! Subprogram not used             if (qtest2 < 0.0_r8) then
! Subprogram not used                eta(i) = 0.99_r8*shb(i,k+1)/(grav*ztodt*fac1)
! Subprogram not used             end if
! Subprogram not used             etagdt(i) = eta(i)*grav*ztodt
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Calculate cloud water, rain water, and thermodynamic changes
! Subprogram not used !
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used          do 30 ii=1,len1
! Subprogram not used             i = indx1(ii)
! Subprogram not used             icwmr(i,k) = cldwtr(i)*rhlat
! Subprogram not used             cldwtr(i) = etagdt(i)*cldwtr(i)*rhlat*rgrav
! Subprogram not used ! JJH changes to facilitate export of cloud liquid water --------------------------------
! Subprogram not used             totcond(i) = (1.0_r8 - beta(i))*cldwtr(i)
! Subprogram not used             rnwtr(i) = min(totcond(i),c0*(dzcld(i)-dzmin)*cldwtr(i))
! Subprogram not used             ds1(i) = etagdt(i)*(sbh(i,k+1) - sc(i))*rpdel(i,k+1)
! Subprogram not used             dq1(i) = etagdt(i)*(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
! Subprogram not used             ds2(i) = (etagdt(i)*(sc(i) - sbh(i,k+1)) +  &
! Subprogram not used                      hlat*grav*cldwtr(i) - beta(i)*etagdt(i)*(sc(i) - sbh(i,k)))*rpdel(i,k)
! Subprogram not used ! JJH change for export of cloud liquid water; must use total condensate 
! Subprogram not used ! since rainwater no longer represents total condensate
! Subprogram not used             dq2(i) = (etagdt(i)*(shc(i) - shbh(i,k+1)) - grav*totcond(i) - beta(i)* &
! Subprogram not used                      etagdt(i)*(shc(i) - shbh(i,k)))*rpdel(i,k)
! Subprogram not used             ds3(i) = beta(i)*(etagdt(i)*(sc(i) - sbh(i,k)) - hlat*grav*cldwtr(i))* &
! Subprogram not used                      rpdel(i,k-1)
! Subprogram not used             dq3(i) = beta(i)*etagdt(i)*(shc(i) - shbh(i,k))*rpdel(i,k-1)
! Subprogram not used !
! Subprogram not used ! Isolate convective fluxes for later diagnostics
! Subprogram not used !
! Subprogram not used             fslkp = eta(i)*(sc(i) - sbh(i,k+1))
! Subprogram not used             fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) - hlat*cldwtr(i)*rtdt)
! Subprogram not used             fqlkp = eta(i)*(shc(i) - shbh(i,k+1))
! Subprogram not used             fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))
! Subprogram not used !
! Subprogram not used ! Update thermodynamic profile (update sb, hb, & hbs later)
! Subprogram not used !
! Subprogram not used             tb (i,k+1) = tb(i,k+1)  + ds1(i)*rcp
! Subprogram not used             tb (i,k  ) = tb(i,k  )  + ds2(i)*rcp
! Subprogram not used             tb (i,k-1) = tb(i,k-1)  + ds3(i)*rcp
! Subprogram not used             shb(i,k+1) = shb(i,k+1) + dq1(i)
! Subprogram not used             shb(i,k  ) = shb(i,k  ) + dq2(i)
! Subprogram not used             shb(i,k-1) = shb(i,k-1) + dq3(i)
! Subprogram not used !
! Subprogram not used ! ** Update diagnostic information for final budget **
! Subprogram not used ! Tracking precipitation, temperature & specific humidity tendencies,
! Subprogram not used ! rainout term, convective mass flux, convective liquid
! Subprogram not used ! water static energy flux, and convective total water flux
! Subprogram not used ! The variable afac makes the necessary adjustment to the
! Subprogram not used ! diagnostic fluxes to account for adjustment time scale based on
! Subprogram not used ! how relaxation time scale is to be applied (column vs. triplet)
! Subprogram not used !
! Subprogram not used             prec(i)    = prec(i) + (rnwtr(i)/rhoh2o)*adjfac
! Subprogram not used !
! Subprogram not used ! The following variables have units of "units"/second
! Subprogram not used !
! Subprogram not used             cmfdt (i,k+1) = cmfdt (i,k+1) + ds1(i)*rtdt*adjfac
! Subprogram not used             cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)*rtdt*adjfac
! Subprogram not used             cmfdt (i,k-1) = cmfdt (i,k-1) + ds3(i)*rtdt*adjfac
! Subprogram not used             cmfdq (i,k+1) = cmfdq (i,k+1) + dq1(i)*rtdt*adjfac
! Subprogram not used             cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rtdt*adjfac
! Subprogram not used             cmfdq (i,k-1) = cmfdq (i,k-1) + dq3(i)*rtdt*adjfac
! Subprogram not used ! JJH changes to export cloud liquid water --------------------------------
! Subprogram not used             qc    (i,k  ) = (grav*(totcond(i)-rnwtr(i))*rpdel(i,k))*rtdt*adjfac
! Subprogram not used             cmfdqr(i,k  ) = cmfdqr(i,k  ) + (grav*rnwtr(i)*rpdel(i,k))*rtdt*adjfac
! Subprogram not used             cmfmc (i,k+1) = cmfmc (i,k+1) + eta(i)*adjfac
! Subprogram not used             cmfmc (i,k  ) = cmfmc (i,k  ) + beta(i)*eta(i)*adjfac
! Subprogram not used !
! Subprogram not used ! The following variables have units of w/m**2
! Subprogram not used !
! Subprogram not used             cmfsl (i,k+1) = cmfsl (i,k+1) + fslkp*adjfac
! Subprogram not used             cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm*adjfac
! Subprogram not used             cmflq (i,k+1) = cmflq (i,k+1) + hlat*fqlkp*adjfac
! Subprogram not used             cmflq (i,k  ) = cmflq (i,k  ) + hlat*fqlkm*adjfac
! Subprogram not used 30          continue
! Subprogram not used !
! Subprogram not used ! Next, convectively modify passive constituents
! Subprogram not used ! For now, when applying relaxation time scale to thermal fields after
! Subprogram not used ! entire column has undergone convective overturning, constituents will
! Subprogram not used ! be mixed using a "relaxed" value of the mass flux determined above
! Subprogram not used ! Although this will be inconsistant with the treatment of the thermal
! Subprogram not used ! fields, it's computationally much cheaper, no more-or-less justifiable,
! Subprogram not used ! and consistent with how the history tape mass fluxes would be used in
! Subprogram not used ! an off-line mode (i.e., using an off-line transport model)
! Subprogram not used !
! Subprogram not used             do 50 m=2,pcnst    ! note: indexing assumes water is first field
! Subprogram not used                if (cnst_get_type_byind(m).eq.'dry') then
! Subprogram not used                   pd(:ncol,:) = pdeldry(:ncol,:)
! Subprogram not used                   rpd(:ncol,:) = rpdeldry(:ncol,:)
! Subprogram not used                   pm(:ncol,:) = pmiddry(:ncol,:)
! Subprogram not used                else
! Subprogram not used                   pd(:ncol,:) = pdel(:ncol,:)
! Subprogram not used                   rpd(:ncol,:) = rpdel(:ncol,:)
! Subprogram not used                   pm(:ncol,:) = pmid(:ncol,:)
! Subprogram not used                endif
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used                do 40 ii=1,len1
! Subprogram not used                   i = indx1(ii)
! Subprogram not used !
! Subprogram not used ! If any of the reported values of the constituent is negative in
! Subprogram not used ! the three adjacent levels, nothing will be done to the profile
! Subprogram not used !
! Subprogram not used                   if ((dq(i,k+1,m) < 0.0_r8) .or. (dq(i,k,m) < 0.0_r8) .or. (dq(i,k-1,m) < 0.0_r8)) go to 40
! Subprogram not used !
! Subprogram not used ! Specify constituent interface values (linear interpolation)
! Subprogram not used !
! Subprogram not used                   cmrh(i,k  ) = 0.5_r8*(dq(i,k-1,m) + dq(i,k  ,m))
! Subprogram not used                   cmrh(i,k+1) = 0.5_r8*(dq(i,k  ,m) + dq(i,k+1,m))
! Subprogram not used !
! Subprogram not used ! Specify perturbation properties of constituents in PBL
! Subprogram not used !
! Subprogram not used                   pblhgt = max(pblh(i),1.0_r8)
! Subprogram not used                   if ( (zm(i,k+1) <= pblhgt) .and. dzcld(i) == 0.0_r8 ) then
! Subprogram not used                      fac1 = max(0.0_r8,1.0_r8-zm(i,k+1)/pblhgt)
! Subprogram not used                      cmrc(i) = dq(i,k+1,m) + qpert(i,m)*fac1
! Subprogram not used                   else
! Subprogram not used                      cmrc(i) = dq(i,k+1,m)
! Subprogram not used                   end if
! Subprogram not used !
! Subprogram not used ! Determine fluxes, flux divergence => changes due to convection
! Subprogram not used ! Logic must be included to avoid producing negative values. A bit
! Subprogram not used ! messy since there are no a priori assumptions about profiles.
! Subprogram not used ! Tendency is modified (reduced) when pending disaster detected.
! Subprogram not used !
! Subprogram not used                   botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
! Subprogram not used                   topflx   = beta(i)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
! Subprogram not used                   dcmr1(i) = -botflx*rpd(i,k+1)
! Subprogram not used                   efac1    = 1.0_r8
! Subprogram not used                   efac2    = 1.0_r8
! Subprogram not used                   efac3    = 1.0_r8
! Subprogram not used !
! Subprogram not used                   if (dq(i,k+1,m)+dcmr1(i) < 0.0_r8) then
! Subprogram not used                      if ( abs(dcmr1(i)) > 1.e-300_r8 ) then
! Subprogram not used                         efac1 = max(tiny,abs(dq(i,k+1,m)/dcmr1(i)) - eps)
! Subprogram not used                      else
! Subprogram not used                         efac1 = tiny
! Subprogram not used                      endif
! Subprogram not used                   end if
! Subprogram not used !
! Subprogram not used                   if (efac1 == tiny .or. efac1 > 1.0_r8) efac1 = 0.0_r8
! Subprogram not used                   dcmr1(i) = -efac1*botflx*rpd(i,k+1)
! Subprogram not used                   dcmr2(i) = (efac1*botflx - topflx)*rpd(i,k)
! Subprogram not used !
! Subprogram not used                   if (dq(i,k,m)+dcmr2(i) < 0.0_r8) then
! Subprogram not used                      if ( abs(dcmr2(i)) > 1.e-300_r8 ) then
! Subprogram not used                         efac2 = max(tiny,abs(dq(i,k  ,m)/dcmr2(i)) - eps)
! Subprogram not used                      else
! Subprogram not used                         efac2 = tiny
! Subprogram not used                      endif
! Subprogram not used                   end if
! Subprogram not used !
! Subprogram not used                   if (efac2 == tiny .or. efac2 > 1.0_r8) efac2 = 0.0_r8
! Subprogram not used                   dcmr2(i) = (efac1*botflx - efac2*topflx)*rpd(i,k)
! Subprogram not used                   dcmr3(i) = efac2*topflx*rpd(i,k-1)
! Subprogram not used !
! Subprogram not used                   if ( (dq(i,k-1,m)+dcmr3(i) < 0.0_r8 ) ) then
! Subprogram not used                      if  ( abs(dcmr3(i)) > 1.e-300_r8 ) then
! Subprogram not used                         efac3 = max(tiny,abs(dq(i,k-1,m)/dcmr3(i)) - eps)
! Subprogram not used                      else
! Subprogram not used                         efac3 = tiny
! Subprogram not used                      endif
! Subprogram not used                   end if
! Subprogram not used !
! Subprogram not used                   if (efac3 == tiny .or. efac3 > 1.0_r8) efac3 = 0.0_r8
! Subprogram not used                   efac3    = min(efac2,efac3)
! Subprogram not used                   dcmr2(i) = (efac1*botflx - efac3*topflx)*rpd(i,k)
! Subprogram not used                   dcmr3(i) = efac3*topflx*rpd(i,k-1)
! Subprogram not used !
! Subprogram not used                   dq(i,k+1,m) = dq(i,k+1,m) + dcmr1(i)
! Subprogram not used                   dq(i,k  ,m) = dq(i,k  ,m) + dcmr2(i)
! Subprogram not used                   dq(i,k-1,m) = dq(i,k-1,m) + dcmr3(i)
! Subprogram not used 40                continue
! Subprogram not used 50                continue                ! end of m=2,pcnst loop
! Subprogram not used !
! Subprogram not used ! Constituent modifications complete
! Subprogram not used !
! Subprogram not used                   if (k == limcnv+1) go to 60
! Subprogram not used !
! Subprogram not used ! Complete update of thermodynamic structure at integer levels
! Subprogram not used ! gather/scatter points that need new values of shbs and gamma
! Subprogram not used !
! Subprogram not used                   do ii=1,len1
! Subprogram not used                      i = indx1(ii)
! Subprogram not used                      vtemp1(ii     ) = tb(i,k)
! Subprogram not used                      vtemp1(ii+len1) = tb(i,k-1)
! Subprogram not used                      vtemp2(ii     ) = pmid(i,k)
! Subprogram not used                      vtemp2(ii+len1) = pmid(i,k-1)
! Subprogram not used                   end do
! Subprogram not used                   call qsat(vtemp1(:2*len1), vtemp2(:2*len1), &
! Subprogram not used                        vtemp5(:2*len1), vtemp3(:2*len1), gam=vtemp4(:2*len1))
! Subprogram not used !cdir nodep
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used                   do ii=1,len1
! Subprogram not used                      i = indx1(ii)
! Subprogram not used                      shbs(i,k  ) = vtemp3(ii     )
! Subprogram not used                      shbs(i,k-1) = vtemp3(ii+len1)
! Subprogram not used                      gam(i,k  ) = vtemp4(ii     )
! Subprogram not used                      gam(i,k-1) = vtemp4(ii+len1)
! Subprogram not used                      sb (i,k  ) = sb(i,k  ) + ds2(i)
! Subprogram not used                      sb (i,k-1) = sb(i,k-1) + ds3(i)
! Subprogram not used                      hb (i,k  ) = sb(i,k  ) + hlat*shb(i,k  )
! Subprogram not used                      hb (i,k-1) = sb(i,k-1) + hlat*shb(i,k-1)
! Subprogram not used                      hbs(i,k  ) = sb(i,k  ) + hlat*shbs(i,k  )
! Subprogram not used                      hbs(i,k-1) = sb(i,k-1) + hlat*shbs(i,k-1)
! Subprogram not used                   end do
! Subprogram not used !
! Subprogram not used ! Update thermodynamic information at half (i.e., interface) levels
! Subprogram not used !
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used                   do ii=1,len1
! Subprogram not used                      i = indx1(ii)
! Subprogram not used                      sbh (i,k) = 0.5_r8*(sb(i,k) + sb(i,k-1))
! Subprogram not used                      shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
! Subprogram not used                      hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
! Subprogram not used                      sbh (i,k-1) = 0.5_r8*(sb(i,k-1) + sb(i,k-2))
! Subprogram not used                      shbh(i,k-1) = qhalf(shb(i,k-2),shb(i,k-1),shbs(i,k-2),shbs(i,k-1))
! Subprogram not used                      hbh (i,k-1) = sbh(i,k-1) + hlat*shbh(i,k-1)
! Subprogram not used                   end do
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Ensure that dzcld is reset if convective mass flux zero
! Subprogram not used ! specify the current vertical extent of the convective activity
! Subprogram not used ! top of convective layer determined by size of overshoot param.
! Subprogram not used !
! Subprogram not used 60                do i=1,ncol
! Subprogram not used                      etagt0 = eta(i).gt.0.0_r8
! Subprogram not used                      if ( .not. etagt0) dzcld(i) = 0.0_r8
! Subprogram not used                      if (etagt0 .and. beta(i) > betamn) then
! Subprogram not used                         ktp = k-1
! Subprogram not used                      else
! Subprogram not used                         ktp = k
! Subprogram not used                      end if
! Subprogram not used                      if (etagt0) then
! Subprogram not used                         rk=k
! Subprogram not used                         rktp=ktp
! Subprogram not used                         cnt(i) = min(cnt(i),rktp)
! Subprogram not used                         cnb(i) = max(cnb(i),rk)
! Subprogram not used                      end if
! Subprogram not used                   end do
! Subprogram not used 70                continue                  ! end of k loop
! Subprogram not used !
! Subprogram not used ! ** apply final thermodynamic tendencies **
! Subprogram not used !
! Subprogram not used !**BAB don't update input profiles
! Subprogram not used !!$                  do k=limcnv,pver
! Subprogram not used !!$                     do i=1,ncol
! Subprogram not used !!$                        t (i,k) = t (i,k) + cmfdt(i,k)*ztodt
! Subprogram not used !!$                        q(i,k,1) = q(i,k,1) + cmfdq(i,k)*ztodt
! Subprogram not used !!$                     end do
! Subprogram not used !!$                  end do
! Subprogram not used ! Set output q tendencies 
! Subprogram not used       dq(:ncol,:,1 ) = cmfdq(:ncol,:)
! Subprogram not used       dq(:ncol,:,2:) = (dq(:ncol,:,2:) - q(:ncol,:,2:))/ztodt
! Subprogram not used !
! Subprogram not used ! Kludge to prevent cnb-cnt from being zero (in the event
! Subprogram not used ! someone decides that they want to divide by this quantity)
! Subprogram not used !
! Subprogram not used                   do i=1,ncol
! Subprogram not used                      if (cnb(i) /= 0.0_r8 .and. cnb(i) == cnt(i)) then
! Subprogram not used                         cnt(i) = cnt(i) - 1.0_r8
! Subprogram not used                      end if
! Subprogram not used                   end do
! Subprogram not used !
! Subprogram not used                   do i=1,ncol
! Subprogram not used                      precc(i) = prec(i)*rtdt
! Subprogram not used                   end do
! Subprogram not used !
! Subprogram not used ! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Subprogram not used ! Treat rliq as flux out bottom, to be added back later.
! Subprogram not used    do k = 1, pver
! Subprogram not used       do i = 1, ncol
! Subprogram not used          rliq(i) = rliq(i) + qc(i,k)*pdel(i,k)/grav
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used    rliq(:ncol) = rliq(:ncol) /1000._r8
! Subprogram not used 
! Subprogram not used                   return                 ! we're all done ... return to calling procedure
! Subprogram not used !
! Subprogram not used end subroutine cmfmca
end module hk_conv
