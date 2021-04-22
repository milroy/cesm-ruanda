module seq_flux_mct
  
  use shr_kind_mod,      only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_sys_mod,       only: shr_sys_abort
  use shr_flux_mod,      only: shr_flux_atmocn
  use shr_orb_mod,       only: shr_orb_params, shr_orb_cosz, shr_orb_decl
  use shr_mct_mod,       only: shr_mct_queryConfigFile, shr_mct_sMatReaddnc

  use mct_mod
  use seq_flds_mod
  use seq_comm_mct
  use seq_infodata_mod

  use component_type_mod

  implicit none
  private 	
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_flux_init_mct
  public seq_flux_initexch_mct

  public seq_flux_ocnalb_mct

  public seq_flux_atmocn_mct
  public seq_flux_atmocnexch_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  real(r8), pointer       :: lats(:)  ! latitudes  (degrees)
  real(r8), pointer       :: lons(:)  ! longitudes (degrees)
  integer(in),allocatable :: mask(:)  ! ocn domain mask: 0 <=> inactive cell
  integer(in),allocatable :: emask(:) ! ocn mask on exchange grid decomp

  real(r8), allocatable ::  uocn (:)  ! ocn velocity, zonal
  real(r8), allocatable ::  vocn (:)  ! ocn velocity, meridional
  real(r8), allocatable ::  tocn (:)  ! ocean temperature
  real(r8), allocatable ::  zbot (:)  ! atm level height
  real(r8), allocatable ::  ubot (:)  ! atm velocity, zonal     
  real(r8), allocatable ::  vbot (:)  ! atm velocity, meridional
  real(r8), allocatable ::  thbot(:)  ! atm potential T
  real(r8), allocatable ::  shum (:)  ! atm specific humidity
  real(r8), allocatable ::  dens (:)  ! atm density
  real(r8), allocatable ::  tbot (:)  ! atm bottom surface T
  real(r8), allocatable ::  sen  (:)  ! heat flux: sensible 
  real(r8), allocatable ::  lat  (:)  ! heat flux: latent   
  real(r8), allocatable ::  lwup (:)  ! lwup over ocean
  real(r8), allocatable ::  evap (:)  ! water flux: evaporation
  real(r8), allocatable ::  taux (:)  ! wind stress, zonal
  real(r8), allocatable ::  tauy (:)  ! wind stress, meridional
  real(r8), allocatable ::  tref (:)  ! diagnostic:  2m ref T
  real(r8), allocatable ::  qref (:)  ! diagnostic:  2m ref Q
  real(r8), allocatable :: duu10n(:)  ! diagnostic: 10m wind speed squared

  real(r8), allocatable ::  ustar(:)  ! saved ustar
  real(r8), allocatable ::  re   (:)  ! saved re
  real(r8), allocatable ::  ssq  (:)  ! saved sq

  ! Conversion from degrees to radians

  real(r8),parameter :: const_pi      = SHR_CONST_PI       ! pi
  real(r8),parameter :: const_deg2rad = const_pi/180.0_R8  ! deg to rads

  ! Coupler field indices

  integer :: index_a2x_Sa_z    
  integer :: index_a2x_Sa_u    
  integer :: index_a2x_Sa_v    
  integer :: index_a2x_Sa_tbot 
  integer :: index_a2x_Sa_ptem 
  integer :: index_a2x_Sa_shum 
  integer :: index_a2x_Sa_dens 
  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_xao_So_tref    
  integer :: index_xao_So_qref    
  integer :: index_xao_So_avsdr   
  integer :: index_xao_So_avsdf   
  integer :: index_xao_So_anidr   
  integer :: index_xao_So_anidf   
  integer :: index_xao_Faox_taux  
  integer :: index_xao_Faox_tauy   
  integer :: index_xao_Faox_lat   
  integer :: index_xao_Faox_sen   
  integer :: index_xao_Faox_evap  
  integer :: index_xao_Faox_lwup
  integer :: index_xao_So_ustar
  integer :: index_xao_So_re   
  integer :: index_xao_So_ssq  
  integer :: index_xao_So_duu10n 
  integer :: index_xao_So_u10

  character(len=16) :: fluxsetting = 'unknown'
  character(len=*),parameter  :: fluxsetting_atmocn = 'atmocn'
  character(len=*),parameter  :: fluxsetting_exchange = 'exchange'

  !--- for exchange grid ---
  type(mct_rearr) :: Re_a2e, Re_e2a, Re_o2e, Re_e2o  ! atm/ocn/exch rearrangers
  type(mct_sMat ) :: sMata2o, sMato2a                ! decomp sMat 
  type(mct_gsMap) :: gsmap_ae, gsmap_oe              ! gsmaps for atm/ocn on exch grid
  integer(in)     :: nloc_a2o,nloc_o2a,nloc_o,nloc_a,nloc_ae,nloc_oe 

!===============================================================================
contains
!===============================================================================

  subroutine seq_flux_init_mct(comp, fractions)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(component_type), intent(in) :: comp
    type(mct_aVect), intent(in)  :: fractions
    !
    ! Local variables
    !
    type(mct_gsMap), pointer :: gsMap
    type(mct_gGrid), pointer :: dom
    integer(in)              :: nloc
    integer                  :: kx,kr     ! fractions indices
    integer                  :: ier
    real(r8), pointer        :: rmask(:)  ! ocn domain mask
    character(*),parameter   :: subName =   '(seq_flux_init_mct) '
    !-----------------------------------------------------------------------

    gsmap => component_get_gsmap_cx(comp) 
    dom   => component_get_dom_cx(comp) 

    nloc = mct_avect_lsize(dom%data)

    ! Input fields atm
    allocate( zbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate zbot',ier)
    allocate( ubot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ubot',ier)
    allocate( vbot(nloc))
    if(ier/=0) call mct_die(subName,'allocate vbot',ier)
    allocate(thbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate thbot',ier)
    allocate(shum(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum',ier)
    allocate(dens(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate dens',ier)
    allocate(tbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tbot',ier)
    allocate(ustar(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ustar',ier)
    allocate(re(nloc), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate re',ier)
    allocate(ssq(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ssq',ier)
    allocate( uocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate uocn',ier)
    allocate( vocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate vocn',ier)
    allocate( tocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tocn',ier)

    ! Output fields 
    allocate(sen (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sen',ier)
    allocate(lat (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lat',ier)
    allocate(evap(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap',ier)
    allocate(lwup(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lwup',ier)
    allocate(taux(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate taux',ier)
    allocate(tauy(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tauy',ier)
    allocate(tref(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tref',ier)
    allocate(qref(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qref',ier)
    allocate(duu10n(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate duu10n',ier)

    ! Grid fields
    allocate( lats(nloc),stat=ier )
    if(ier/=0) call mct_die(subName,'allocate lats',ier)
    allocate( lons(nloc),stat=ier )
    if(ier/=0) call mct_die(subName,'allocate lons',ier)
    allocate( emask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate emask',ier)
    allocate(mask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate mask',ier)
    
    ! Get lat, lon, mask, which is time-invariant
    allocate(rmask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate rmask',ier)
    call mct_gGrid_exportRAttr(dom, 'lat' , lats , nloc) 
    call mct_gGrid_exportRAttr(dom, 'lon' , lons , nloc) 
    call mct_gGrid_exportRAttr(dom, 'mask', rmask, nloc)
    !tcx, want to mask properly, but applying this changes answers to roundoff for some reason
    !    kx = mct_aVect_indexRA(fractions,"ofrac")
    !    mask = 0
    !    where (fractions%rAttr(kx,:) > 0.0_r8) mask = nint(rmask)
    mask = nint(rmask)
    deallocate(rmask)
    emask = mask

    fluxsetting = trim(fluxsetting_atmocn)

  end subroutine seq_flux_init_mct

!===============================================================================

! Subprogram not used   subroutine seq_flux_initexch_mct(atm, ocn, mpicom_cplid, cplid)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     type(component_type), intent(in) :: atm 
! Subprogram not used     type(component_type), intent(in) :: ocn
! Subprogram not used     integer(in)         , intent(in) :: mpicom_cplid
! Subprogram not used     integer(in)         , intent(in) :: cplid
! Subprogram not used     !
! Subprogram not used     ! Local variables
! Subprogram not used     !
! Subprogram not used     type(mct_gsMap), pointer :: gsmap_a
! Subprogram not used     type(mct_gGrid), pointer :: dom_a
! Subprogram not used     type(mct_gsMap), pointer :: gsmap_o
! Subprogram not used     type(mct_gGrid), pointer :: dom_o
! Subprogram not used     integer(in)              :: kw,ka,ko,iw,ia,io,n
! Subprogram not used     character(len=128)       :: strat
! Subprogram not used     integer                  :: ier
! Subprogram not used     integer                  :: mytask
! Subprogram not used     integer(in)              :: kmsk            ! field indices
! Subprogram not used     character(len=128)       :: ConfigFileName  ! config file to read
! Subprogram not used     character(len=128)       :: MapLabel        ! map name
! Subprogram not used     character(len=128)       :: MapTypeLabel    ! map type
! Subprogram not used     character(len=256)       :: fileName
! Subprogram not used     character(len=1)         :: maptype
! Subprogram not used     character(len=3)         :: Smaptype
! Subprogram not used     type(mct_aVect)          :: avdom_oe
! Subprogram not used     type(mct_list)           :: sort_keys
! Subprogram not used     character(*),parameter :: subName =   '(seq_flux_initexch_mct) '
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     gsmap_a => component_get_gsmap_cx(atm) ! gsmap_ax
! Subprogram not used     gsmap_o => component_get_gsmap_cx(ocn) ! gsmap_ox
! Subprogram not used     dom_a   => component_get_dom_cx(atm)   ! dom_ax
! Subprogram not used     dom_o   => component_get_dom_cx(ocn)   ! dom_ox
! Subprogram not used 
! Subprogram not used     call shr_mpi_commrank(mpicom_cplid, mytask)
! Subprogram not used 
! Subprogram not used     !--- Get mapping file info
! Subprogram not used     do n = 1,2
! Subprogram not used        ConfigFileName = "seq_maps.rc"
! Subprogram not used        if (n == 1) then
! Subprogram not used           MapLabel = "atm2ocn_fmapname:"
! Subprogram not used           MapTypeLabel = "atm2ocn_fmaptype:"
! Subprogram not used        elseif (n == 2) then
! Subprogram not used           MapLabel = "ocn2atm_fmapname:"
! Subprogram not used           MapTypeLabel = "ocn2atm_fmaptype:"
! Subprogram not used        else
! Subprogram not used           call shr_sys_abort(trim(subname)//' do error1')
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        call shr_mct_queryConfigFile(mpicom_cplid, ConfigFilename, &
! Subprogram not used           trim(MapLabel),fileName,trim(MapTypeLabel),maptype)
! Subprogram not used 
! Subprogram not used        !--- hardwire decomposition to gsmap_o
! Subprogram not used        if (n == 1) then
! Subprogram not used           Smaptype = "src"
! Subprogram not used           call shr_mct_sMatReaddnc(sMata2o, gsmap_a, gsmap_o, Smaptype, &
! Subprogram not used              filename=fileName, mytask=mytask, mpicom=mpicom_cplid)
! Subprogram not used        elseif (n == 2) then
! Subprogram not used           Smaptype = "dst"
! Subprogram not used           call shr_mct_sMatReaddnc(sMato2a, gsmap_o, gsmap_a, Smaptype, &
! Subprogram not used              filename=fileName, mytask=mytask, mpicom=mpicom_cplid)
! Subprogram not used        else
! Subprogram not used           call shr_sys_abort(trim(subname)//' do error2')
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     !--- the two mapping files must have their local indices in identical order
! Subprogram not used     !--- sort the global indices as a starting point
! Subprogram not used 
! Subprogram not used     call mct_list_init(sort_keys,'grow:gcol')
! Subprogram not used     call mct_sMat_SortPermute(sMata2o,sort_keys)
! Subprogram not used     call mct_list_clean(sort_keys)
! Subprogram not used     call mct_list_init(sort_keys,'gcol:grow')
! Subprogram not used     call mct_sMat_SortPermute(sMato2a,sort_keys)
! Subprogram not used     call mct_list_clean(sort_keys)
! Subprogram not used 
! Subprogram not used     !--- now check that they are sorted properly
! Subprogram not used 
! Subprogram not used     nloc_a2o= mct_sMat_lsize(sMata2o)
! Subprogram not used     nloc_o2a= mct_sMat_lsize(sMato2a)
! Subprogram not used 
! Subprogram not used     if (nloc_a2o /= nloc_o2a) then
! Subprogram not used        write(logunit,*) trim(subname),' ERROR: sMat sizes',nloc_a2o,nloc_o2a
! Subprogram not used        call shr_sys_abort(trim(subname)//' ERROR in sMat sizes')
! Subprogram not used     endif
! Subprogram not used     ko = mct_sMat_indexIA(sMata2o,'grow')    ! local row (dst) index
! Subprogram not used     ka = mct_sMat_indexIA(sMato2a,'gcol')    ! local column (src) index
! Subprogram not used     do n = 1,nloc_a2o
! Subprogram not used        io = sMata2o%data%iAttr(ko,n)
! Subprogram not used        ia = sMato2a%data%iAttr(ka,n)
! Subprogram not used        if (io /= ia) then
! Subprogram not used           write(logunit,*) trim(subname),' ERROR: sMat indices1 ',io,ia
! Subprogram not used           call shr_sys_abort(trim(subname)//' ERROR in sMat indices1')
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     ko = mct_sMat_indexIA(sMata2o,'gcol')    ! local column (src) index
! Subprogram not used     ka = mct_sMat_indexIA(sMato2a,'grow')    ! local row (dst) index
! Subprogram not used     do n = 1,nloc_a2o
! Subprogram not used        io = sMata2o%data%iAttr(ko,n)
! Subprogram not used        ia = sMato2a%data%iAttr(ka,n)
! Subprogram not used        if (io /= ia) then
! Subprogram not used           write(logunit,*) trim(subname),' ERROR: sMat indices2 ',io,ia
! Subprogram not used           call shr_sys_abort(trim(subname)//' ERROR in sMat indices2')
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     !--- instantiate/create/compute various datatypes
! Subprogram not used 
! Subprogram not used     call mct_sMat_2XgsMap(sMata2o , gsmap_ae, 0, mpicom_cplid, cplid)
! Subprogram not used     call mct_sMat_2YgsMap(sMata2o , gsmap_oe, 0, mpicom_cplid, cplid)
! Subprogram not used 
! Subprogram not used     call mct_rearr_init(gsmap_a   , gsmap_ae, mpicom_cplid, Re_a2e)
! Subprogram not used     call mct_rearr_init(gsmap_ae  , gsmap_a,  mpicom_cplid, Re_e2a)
! Subprogram not used     call mct_rearr_init(gsmap_o   , gsmap_oe, mpicom_cplid, Re_o2e)
! Subprogram not used     call mct_rearr_init(gsmap_oe  , gsmap_o,  mpicom_cplid, Re_e2o)
! Subprogram not used 
! Subprogram not used     call mct_sMat_g2lMat(sMata2o  , gsmap_ae, 'column',mpicom_cplid)
! Subprogram not used     call mct_sMat_g2lMat(sMata2o  , gsmap_oe, 'row',   mpicom_cplid)
! Subprogram not used     call mct_sMat_g2lMat(sMato2a  , gsmap_ae, 'row',   mpicom_cplid)
! Subprogram not used     call mct_sMat_g2lMat(sMato2a  , gsmap_oe, 'column',mpicom_cplid)
! Subprogram not used 
! Subprogram not used     nloc_a  = mct_gsmap_lsize(gsmap_a  , mpicom_cplid)
! Subprogram not used     nloc_o  = mct_gsmap_lsize(gsmap_o  , mpicom_cplid)
! Subprogram not used     nloc_ae = mct_gsmap_lsize(gsmap_ae , mpicom_cplid)
! Subprogram not used     nloc_oe = mct_gsmap_lsize(gsmap_oe , mpicom_cplid)
! Subprogram not used 
! Subprogram not used     call mct_gsmap_clean(gsmap_ae)
! Subprogram not used     call mct_gsmap_clean(gsmap_oe)
! Subprogram not used 
! Subprogram not used     ! Input fields atm
! Subprogram not used     allocate( emask(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate emask',ier)
! Subprogram not used     allocate( zbot(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate zbot',ier)
! Subprogram not used     allocate( ubot(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate ubot',ier)
! Subprogram not used     allocate( vbot(nloc_a2o))
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate vbot',ier)
! Subprogram not used     allocate(thbot(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate thbot',ier)
! Subprogram not used     allocate(shum(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate shum',ier)
! Subprogram not used     allocate(dens(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate dens',ier)
! Subprogram not used     allocate(tbot(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate tbot',ier)
! Subprogram not used     allocate(ustar(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate ustar',ier)
! Subprogram not used     allocate(re(nloc_a2o), stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate re',ier)
! Subprogram not used     allocate(ssq(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate ssq',ier)
! Subprogram not used     allocate( uocn(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate uocn',ier)
! Subprogram not used     allocate( vocn(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate vocn',ier)
! Subprogram not used     allocate( tocn(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate tocn',ier)
! Subprogram not used 
! Subprogram not used     ! Output fields 
! Subprogram not used     allocate(sen (nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate sen',ier)
! Subprogram not used     allocate(lat (nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate lat',ier)
! Subprogram not used     allocate(evap(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate evap',ier)
! Subprogram not used     allocate(lwup(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate lwup',ier)
! Subprogram not used     allocate(taux(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate taux',ier)
! Subprogram not used     allocate(tauy(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate tauy',ier)
! Subprogram not used     allocate(tref(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate tref',ier)
! Subprogram not used     allocate(qref(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate qref',ier)
! Subprogram not used     allocate(duu10n(nloc_a2o),stat=ier)
! Subprogram not used     if(ier/=0) call mct_die(subName,'allocate duu10n',ier)
! Subprogram not used 
! Subprogram not used     ! set emask
! Subprogram not used 
! Subprogram not used     call mct_avect_init(avdom_oe,dom_o%data,lsize=nloc_oe)
! Subprogram not used     call mct_rearr_rearrange(dom_o%data, avdom_oe, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
! Subprogram not used     ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local dst index
! Subprogram not used     kmsk = mct_aVect_indexRA(avdom_oe,"mask",dieWith=subName)
! Subprogram not used     do n = 1,nloc_a2o
! Subprogram not used        io = sMata2o%data%iAttr(ko,n)
! Subprogram not used        emask(n) = nint(avdom_oe%rAttr(kmsk,io))
! Subprogram not used        if (emask(n) == 0) then
! Subprogram not used           write(logunit,*) trim(subname),' ERROR: weights use masked ocean value'
! Subprogram not used           call shr_sys_abort(trim(subname)//' ERROR: weights use masked ocean value')
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call mct_aVect_clean(avdom_oe)
! Subprogram not used 
! Subprogram not used     fluxsetting = trim(fluxsetting_exchange)
! Subprogram not used 
! Subprogram not used   end subroutine seq_flux_initexch_mct

!===============================================================================

  subroutine seq_flux_ocnalb_mct( infodata, ocn, fractions_o, xao_o )

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_infodata_type) , intent(in)    :: infodata
    type(component_type)    , intent(in)    :: ocn
    type(mct_aVect)         , intent(inout) :: xao_o
    type(mct_aVect)         , intent(inout) :: fractions_o
    !
    ! Local variables
    !
    type(mct_gGrid), pointer :: dom_o
    logical		:: flux_albav		! flux avg option
    integer(in)		:: n,i			! indices
    real(r8)		:: rlat			! gridcell latitude in radians
    real(r8)		:: rlon			! gridcell longitude in radians
    real(r8)		:: cosz			! Cosine of solar zenith angle
    real(r8)		:: eccen		! Earth orbit eccentricity
    real(r8)		:: mvelpp		! Earth orbit
    real(r8)		:: lambm0		! Earth orbit
    real(r8)		:: obliqr		! Earth orbit
    real(r8)		:: delta		! Solar declination angle  in radians
    real(r8)		:: eccf			! Earth orbit eccentricity factor
    real(r8)		:: calday		! calendar day including fraction, at 0e
    real(r8)		:: nextsw_cday		! calendar day of next atm shortwave
    real(r8)		:: anidr		! albedo: near infrared, direct
    real(r8)		:: avsdr		! albedo: visible      , direct
    real(r8)		:: anidf		! albedo: near infrared, diffuse
    real(r8)		:: avsdf		! albedo: visible      , diffuse
    integer(in)		:: ID			! comm ID
    integer(in)		:: ier			! error code
    integer(in)		:: kx,kr		! fractions indices
    integer(in)		:: klat,klon,kmsk	! field indices
    logical		:: update_alb		! was albedo updated
    logical,save	:: first_call = .true. 
    !
    real(R8),parameter :: albdif = 0.06_R8 ! 60 deg reference albedo, diffuse
    real(R8),parameter :: albdir = 0.07_R8 ! 60 deg reference albedo, direct 
    character(*),parameter :: subName =   '(seq_flux_ocnalb_mct) '
    !
    !-----------------------------------------------------------------------

    dom_o => component_get_dom_cx(ocn) ! dom_ox

    call seq_infodata_getData(infodata , &
         flux_albav=flux_albav)

    ! Determine indices

    update_alb = .false.

    if (first_call) then
       index_xao_So_anidr  = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_anidf  = mct_aVect_indexRA(xao_o,'So_anidf')
       index_xao_So_avsdr  = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_avsdf  = mct_aVect_indexRA(xao_o,'So_avsdf')

       nloc_o  = mct_ggrid_lsize(dom_o)
       klat = mct_gGrid_indexRA(dom_o,"lat" ,dieWith=subName)
       klon = mct_gGrid_indexRA(dom_o,"lon" ,dieWith=subName)
       allocate( lats(nloc_o),stat=ier )
       if(ier/=0) call mct_die(subName,'allocate lats',ier)
       allocate( lons(nloc_o),stat=ier )
       if(ier/=0) call mct_die(subName,'allocate lons',ier)
       do n = 1,nloc_o
          lats(n) = dom_o%data%rAttr(klat,n)
          lons(n) = dom_o%data%rAttr(klon,n)
       enddo
       first_call = .false.
    endif

    if (flux_albav) then

       do n=1,nloc_o   
          anidr = albdir
          avsdr = albdir
          anidf = albdif
          avsdf = albdif

          ! Albedo is now function of latitude (will be new implementation)
          !rlat = const_deg2rad * lats(n)
          !anidr = 0.069_R8 - 0.011_R8 * cos(2._R8 * rlat)
          !avsdr = anidr
          !anidf = anidr
          !avsdf = anidr

          xao_o%rAttr(index_xao_So_avsdr,n) = avsdr
          xao_o%rAttr(index_xao_So_anidr,n) = anidr
          xao_o%rAttr(index_xao_So_avsdf,n) = avsdf
          xao_o%rAttr(index_xao_So_anidf,n) = anidf

       end do
       update_alb = .true.

    else

       ! Solar declination 
       ! Will only do albedo calculation if nextsw_cday is not -1.
       
       call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday,orb_eccen=eccen, &
          orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr)
       if (nextsw_cday >= -0.5_r8) then
          calday = nextsw_cday
          call shr_orb_decl(calday, eccen, mvelpp,lambm0, obliqr, delta, eccf)
          ! Compute albedos 
          do n=1,nloc_o
             rlat = const_deg2rad * lats(n)
             rlon = const_deg2rad * lons(n)
             cosz = shr_orb_cosz( calday, rlat, rlon, delta )
             if (cosz  >  0.0_R8) then !--- sun hit --
                anidr = (.026_R8/(cosz**1.7_R8 + 0.065_R8)) +   &
                        (.150_R8*(cosz         - 0.100_R8 ) *   &
                                 (cosz         - 0.500_R8 ) *   &
                                 (cosz         - 1.000_R8 )  )
                avsdr = anidr
                anidf = albdif
                avsdf = albdif
             else !--- dark side of earth ---
                anidr = 1.0_R8
                avsdr = 1.0_R8
                anidf = 1.0_R8
                avsdf = 1.0_R8
             end if

             xao_o%rAttr(index_xao_So_avsdr,n) = avsdr
             xao_o%rAttr(index_xao_So_anidr,n) = anidr
             xao_o%rAttr(index_xao_So_avsdf,n) = avsdf
             xao_o%rAttr(index_xao_So_anidf,n) = anidf

          end do   ! nloc_o
          update_alb = .true.
       endif    ! nextsw_cday
    end if   ! flux_albav

    !--- update current ifrad/ofrad values if albedo was updated

    if (update_alb) then
       kx = mct_aVect_indexRA(fractions_o,"ifrac")
       kr = mct_aVect_indexRA(fractions_o,"ifrad")
       fractions_o%rAttr(kr,:) = fractions_o%rAttr(kx,:)
       kx = mct_aVect_indexRA(fractions_o,"ofrac")
       kr = mct_aVect_indexRA(fractions_o,"ofrad")
       fractions_o%rAttr(kr,:) = fractions_o%rAttr(kx,:)
    endif
       
    end subroutine seq_flux_ocnalb_mct

!===============================================================================

! Subprogram not used   subroutine seq_flux_atmocnexch_mct( infodata, atm, ocn, fractions_a, fractions_o, &
! Subprogram not used        xao_a, xao_o)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     type(seq_infodata_type) , intent(in)    :: infodata
! Subprogram not used     type(component_type)    , intent(in)    :: atm
! Subprogram not used     type(component_type)    , intent(in)    :: ocn
! Subprogram not used     type(mct_aVect)         , intent(in)    :: fractions_a
! Subprogram not used     type(mct_aVect)         , intent(in)    :: fractions_o
! Subprogram not used     type(mct_aVect)         , intent(inout) :: xao_a
! Subprogram not used     type(mct_aVect)         , intent(inout) :: xao_o
! Subprogram not used     !
! Subprogram not used     ! Local variables
! Subprogram not used     !
! Subprogram not used     type(mct_aVect) , pointer :: a2x
! Subprogram not used     type(mct_aVect) , pointer :: o2x
! Subprogram not used     type(mct_gsmap) , pointer :: gsmap_a
! Subprogram not used     type(mct_gsmap) , pointer :: gsmap_o
! Subprogram not used 
! Subprogram not used     type(mct_aVect) :: a2x_e
! Subprogram not used     type(mct_aVect) :: o2x_e
! Subprogram not used     type(mct_aVect) :: xaop_ae
! Subprogram not used     type(mct_aVect) :: xaop_oe
! Subprogram not used     type(mct_aVect) :: xaop_a
! Subprogram not used     type(mct_aVect) :: xaop_o
! Subprogram not used     type(mct_aVect) :: fractions_oe
! Subprogram not used 
! Subprogram not used     integer(in) :: kw,ka,ko,iw,ia,io,kf
! Subprogram not used     integer(in) :: n,i          ! indices
! Subprogram not used     logical     :: dead_comps   ! .true.  => dead components are used
! Subprogram not used     integer(in) :: index_tref  
! Subprogram not used     integer(in) :: index_qref  
! Subprogram not used     integer(in) :: index_duu10n
! Subprogram not used     integer(in) :: index_ustar 
! Subprogram not used     integer(in) :: index_ssq   
! Subprogram not used     integer(in) :: index_re    
! Subprogram not used     integer(in) :: index_u10   
! Subprogram not used     integer(in) :: index_taux  
! Subprogram not used     integer(in) :: index_tauy  
! Subprogram not used     integer(in) :: index_lat   
! Subprogram not used     integer(in) :: index_sen   
! Subprogram not used     integer(in) :: index_evap  
! Subprogram not used     integer(in) :: index_lwup  
! Subprogram not used     integer(in) :: index_sumwt
! Subprogram not used     integer(in) :: atm_nx,atm_ny,ocn_nx,ocn_ny
! Subprogram not used     real(r8)    :: wt
! Subprogram not used     character(len=256) :: fldlist  ! subset of xao fields
! Subprogram not used     !
! Subprogram not used     character(*),parameter :: subName =   '(seq_flux_atmocnexch_mct) '
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     gsmap_a => component_get_gsmap_cx(atm) 
! Subprogram not used     gsmap_o => component_get_gsmap_cx(ocn) 
! Subprogram not used     a2x     => component_get_c2x_cx(atm)  ! a2x_ax
! Subprogram not used     o2x     => component_get_c2x_cx(ocn)  ! o2x_ox 
! Subprogram not used 
! Subprogram not used     if (trim(fluxsetting) /= trim(fluxsetting_exchange)) then
! Subprogram not used        call shr_sys_abort(trim(subname)//' ERROR with init')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! Update ocean surface fluxes 
! Subprogram not used     ! Must fabricate "reasonable" data (using dead components)
! Subprogram not used 
! Subprogram not used     call seq_infodata_GetData(infodata, &
! Subprogram not used          dead_comps=dead_comps,         &
! Subprogram not used          atm_nx=atm_nx, atm_ny=atm_ny,  &
! Subprogram not used          ocn_nx=ocn_nx, ocn_ny=ocn_ny)
! Subprogram not used 
! Subprogram not used     if (dead_comps) then
! Subprogram not used        do n = 1,nloc_a2o
! Subprogram not used           tocn(n) = 290.0_R8 ! ocn temperature            ~ Kelvin
! Subprogram not used           uocn(n) =   0.0_R8 ! ocn velocity, zonal        ~ m/s
! Subprogram not used           vocn(n) =   0.0_R8 ! ocn velocity, meridional   ~ m/s
! Subprogram not used           zbot(n) =  55.0_R8 ! atm height of bottom layer ~ m
! Subprogram not used           ubot(n) =   0.0_R8 ! atm velocity, zonal        ~ m/s
! Subprogram not used           vbot(n) =   2.0_R8 ! atm velocity, meridional   ~ m/s
! Subprogram not used           thbot(n)= 301.0_R8 ! atm potential temperature  ~ Kelvin
! Subprogram not used           shum(n) = 1.e-2_R8 ! atm specific humidity      ~ kg/kg
! Subprogram not used           dens(n) =   1.0_R8 ! atm density                ~ kg/m^3
! Subprogram not used           tbot(n) = 300.0_R8 ! atm temperature            ~ Kelvin
! Subprogram not used        enddo
! Subprogram not used     else        
! Subprogram not used 
! Subprogram not used        !--- instantiate exchange grid aVects
! Subprogram not used        call mct_AVect_init(a2x_e, a2x, nloc_ae)
! Subprogram not used        call mct_AVect_zero(a2x_e)
! Subprogram not used        call mct_AVect_init(o2x_e, o2x, nloc_oe)
! Subprogram not used        call mct_AVect_zero(o2x_e)
! Subprogram not used 
! Subprogram not used        !--- rearrange a2x and o2x into exchange grid
! Subprogram not used 
! Subprogram not used        call mct_rearr_rearrange(a2x, a2x_e, Re_a2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
! Subprogram not used        call mct_rearr_rearrange(o2x, o2x_e, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
! Subprogram not used 
! Subprogram not used        !--- extract fields from a2x and o2x (_e) into local arrays on exchange grid
! Subprogram not used 
! Subprogram not used        ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local row index
! Subprogram not used        ka = mct_sMat_indexIA(sMata2o,'lcol')    ! local column index
! Subprogram not used 
! Subprogram not used        do n = 1,nloc_a2o
! Subprogram not used           io = sMata2o%data%iAttr(ko,n)
! Subprogram not used           ia = sMata2o%data%iAttr(ka,n)
! Subprogram not used           zbot(n) = a2x_e%rAttr(index_a2x_Sa_z   ,ia)
! Subprogram not used           ubot(n) = a2x_e%rAttr(index_a2x_Sa_u   ,ia)
! Subprogram not used           vbot(n) = a2x_e%rAttr(index_a2x_Sa_v   ,ia)
! Subprogram not used           thbot(n)= a2x_e%rAttr(index_a2x_Sa_ptem,ia)
! Subprogram not used           shum(n) = a2x_e%rAttr(index_a2x_Sa_shum,ia)
! Subprogram not used           dens(n) = a2x_e%rAttr(index_a2x_Sa_dens,ia)
! Subprogram not used           tbot(n) = a2x_e%rAttr(index_a2x_Sa_tbot,ia)
! Subprogram not used           tocn(n) = o2x_e%rAttr(index_o2x_So_t   ,io)   
! Subprogram not used           uocn(n) = o2x_e%rAttr(index_o2x_So_u   ,io)
! Subprogram not used           vocn(n) = o2x_e%rAttr(index_o2x_So_v   ,io)
! Subprogram not used        enddo
! Subprogram not used        call mct_aVect_clean(a2x_e)
! Subprogram not used        call mct_aVect_clean(o2x_e)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call shr_flux_atmocn (nloc_a2o , zbot , ubot, vbot, thbot, &
! Subprogram not used                           shum , dens , tbot, uocn, vocn , &
! Subprogram not used                           tocn , emask, sen , lat , lwup , &
! Subprogram not used                           evap , taux , tauy, tref, qref , &
! Subprogram not used                           duu10n,ustar, re  , ssq , missval = 0.0_r8 )
! Subprogram not used 
! Subprogram not used     !--- create temporary aVects on exchange, atm, or ocn decomp as needed
! Subprogram not used 
! Subprogram not used     fldlist = trim(seq_flds_xao_states)//":"//trim(seq_flds_xao_fluxes)//":sumwt"
! Subprogram not used     call mct_aVect_init(xaop_ae,rList=trim(fldlist),lsize=nloc_ae)
! Subprogram not used     call mct_aVect_zero(xaop_ae)
! Subprogram not used     call mct_aVect_init(xaop_oe,rList=trim(fldlist),lsize=nloc_oe)
! Subprogram not used     call mct_aVect_zero(xaop_oe)
! Subprogram not used     call mct_aVect_init(xaop_a, rList=trim(fldlist),lsize=nloc_a)
! Subprogram not used     call mct_aVect_zero(xaop_a)
! Subprogram not used     call mct_aVect_init(xaop_o, rList=trim(fldlist),lsize=nloc_o)
! Subprogram not used     call mct_aVect_zero(xaop_o)
! Subprogram not used 
! Subprogram not used     index_tref   = mct_aVect_indexRA(xaop_ae,"So_tref")
! Subprogram not used     index_qref   = mct_aVect_indexRA(xaop_ae,"So_qref")
! Subprogram not used     index_duu10n = mct_aVect_indexRA(xaop_ae,"So_duu10n")
! Subprogram not used     index_ustar  = mct_aVect_indexRA(xaop_ae,"So_ustar")
! Subprogram not used     index_ssq    = mct_aVect_indexRA(xaop_ae,"So_ssq")
! Subprogram not used     index_re     = mct_aVect_indexRA(xaop_ae,"So_re")
! Subprogram not used     index_u10    = mct_aVect_indexRA(xaop_ae,"So_u10")
! Subprogram not used     index_taux   = mct_aVect_indexRA(xaop_ae,"Faox_taux")
! Subprogram not used     index_tauy   = mct_aVect_indexRA(xaop_ae,"Faox_tauy")
! Subprogram not used     index_lat    = mct_aVect_indexRA(xaop_ae,"Faox_lat")
! Subprogram not used     index_sen    = mct_aVect_indexRA(xaop_ae,"Faox_sen")
! Subprogram not used     index_evap   = mct_aVect_indexRA(xaop_ae,"Faox_evap")
! Subprogram not used     index_lwup   = mct_aVect_indexRA(xaop_ae,"Faox_lwup")
! Subprogram not used     index_sumwt  = mct_aVect_indexRA(xaop_ae,"sumwt")
! Subprogram not used 
! Subprogram not used     !--- aggregate ocean values locally based on exchange grid decomp
! Subprogram not used 
! Subprogram not used     ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local row index
! Subprogram not used     ka = mct_sMat_indexIA(sMata2o,'lcol')    ! local column index
! Subprogram not used     kw = mct_sMat_indexRA(sMata2o,'weight')  ! weight index
! Subprogram not used 
! Subprogram not used     do n = 1,nloc_a2o
! Subprogram not used        io = sMata2o%data%iAttr(ko,n)
! Subprogram not used        ia = sMata2o%data%iAttr(ka,n)
! Subprogram not used        wt = sMata2o%data%rAttr(kw,n)
! Subprogram not used        xaop_oe%rAttr(index_sen   ,io) = xaop_oe%rAttr(index_sen   ,io) + sen(n) * wt
! Subprogram not used        xaop_oe%rAttr(index_lat   ,io) = xaop_oe%rAttr(index_lat   ,io) + lat(n) * wt
! Subprogram not used        xaop_oe%rAttr(index_taux  ,io) = xaop_oe%rAttr(index_taux  ,io) + taux(n)* wt
! Subprogram not used        xaop_oe%rAttr(index_tauy  ,io) = xaop_oe%rAttr(index_tauy  ,io) + tauy(n)* wt
! Subprogram not used        xaop_oe%rAttr(index_evap  ,io) = xaop_oe%rAttr(index_evap  ,io) + evap(n)* wt
! Subprogram not used        xaop_oe%rAttr(index_tref  ,io) = xaop_oe%rAttr(index_tref  ,io) + tref(n)* wt
! Subprogram not used        xaop_oe%rAttr(index_qref  ,io) = xaop_oe%rAttr(index_qref  ,io) + qref(n)* wt
! Subprogram not used        xaop_oe%rAttr(index_ustar ,io) = xaop_oe%rAttr(index_ustar ,io) + ustar(n)*wt   ! friction velocity
! Subprogram not used        xaop_oe%rAttr(index_re    ,io) = xaop_oe%rAttr(index_re    ,io) + re(n)  * wt   ! reynolds number
! Subprogram not used        xaop_oe%rAttr(index_ssq   ,io) = xaop_oe%rAttr(index_ssq   ,io) + ssq(n) * wt   ! s.hum. saturation at Ts
! Subprogram not used        xaop_oe%rAttr(index_lwup  ,io) = xaop_oe%rAttr(index_lwup  ,io) + lwup(n)* wt   
! Subprogram not used        xaop_oe%rAttr(index_duu10n,io) = xaop_oe%rAttr(index_duu10n,io) + duu10n(n)*wt  
! Subprogram not used        xaop_oe%rAttr(index_u10   ,io) = xaop_oe%rAttr(index_u10   ,io) + sqrt(duu10n(n))*wt
! Subprogram not used        xaop_oe%rAttr(index_sumwt ,io) = xaop_oe%rAttr(index_sumwt ,io) + wt
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     !--- aggregate atm values locally based on exchange grid decomp
! Subprogram not used 
! Subprogram not used     ko = mct_sMat_indexIA(sMato2a,'lcol')    ! local column index
! Subprogram not used     ka = mct_sMat_indexIA(sMato2a,'lrow')    ! local row index
! Subprogram not used     kw = mct_sMat_indexRA(sMato2a,'weight')  ! weight index
! Subprogram not used     kf = mct_aVect_indexRA(fractions_o,"ofrac")
! Subprogram not used 
! Subprogram not used     !--- to apply fraction corrections, the indexing must be correct so rearrange
! Subprogram not used     call mct_avect_init(fractions_oe,fractions_o,lsize=nloc_oe)
! Subprogram not used     call mct_rearr_rearrange(fractions_o, fractions_oe, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
! Subprogram not used     do n = 1,nloc_o2a
! Subprogram not used        io = sMato2a%data%iAttr(ko,n)
! Subprogram not used        ia = sMato2a%data%iAttr(ka,n)
! Subprogram not used !tcx   wt = sMato2a%data%rAttr(kw,n)
! Subprogram not used        wt = sMato2a%data%rAttr(kw,n) * fractions_oe%rAttr(kf,io)
! Subprogram not used        xaop_ae%rAttr(index_sen   ,ia) = xaop_ae%rAttr(index_sen   ,ia) + sen(n) * wt
! Subprogram not used        xaop_ae%rAttr(index_lat   ,ia) = xaop_ae%rAttr(index_lat   ,ia) + lat(n) * wt
! Subprogram not used        xaop_ae%rAttr(index_taux  ,ia) = xaop_ae%rAttr(index_taux  ,ia) + taux(n)* wt
! Subprogram not used        xaop_ae%rAttr(index_tauy  ,ia) = xaop_ae%rAttr(index_tauy  ,ia) + tauy(n)* wt
! Subprogram not used        xaop_ae%rAttr(index_evap  ,ia) = xaop_ae%rAttr(index_evap  ,ia) + evap(n)* wt
! Subprogram not used        xaop_ae%rAttr(index_tref  ,ia) = xaop_ae%rAttr(index_tref  ,ia) + tref(n)* wt
! Subprogram not used        xaop_ae%rAttr(index_qref  ,ia) = xaop_ae%rAttr(index_qref  ,ia) + qref(n)* wt
! Subprogram not used        xaop_ae%rAttr(index_ustar ,ia) = xaop_ae%rAttr(index_ustar ,ia) + ustar(n)*wt   ! friction velocity
! Subprogram not used        xaop_ae%rAttr(index_re    ,ia) = xaop_ae%rAttr(index_re    ,ia) + re(n)  * wt   ! reynolds number
! Subprogram not used        xaop_ae%rAttr(index_ssq   ,ia) = xaop_ae%rAttr(index_ssq   ,ia) + ssq(n) * wt   ! s.hum. saturation at Ts
! Subprogram not used        xaop_ae%rAttr(index_lwup  ,ia) = xaop_ae%rAttr(index_lwup  ,ia) + lwup(n)* wt   
! Subprogram not used        xaop_ae%rAttr(index_duu10n,ia) = xaop_ae%rAttr(index_duu10n,ia) + duu10n(n)*wt  
! Subprogram not used        xaop_ae%rAttr(index_u10   ,ia) = xaop_ae%rAttr(index_u10   ,ia) + sqrt(duu10n(n))*wt
! Subprogram not used        xaop_ae%rAttr(index_sumwt ,ia) = xaop_ae%rAttr(index_sumwt ,ia) + wt
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call mct_aVect_clean(fractions_oe)
! Subprogram not used 
! Subprogram not used     !--- rearrange and sum from exchange grid to gsmap_a and gsmap_o decomps
! Subprogram not used 
! Subprogram not used     call mct_rearr_rearrange(xaop_ae, xaop_a, Re_e2a, sum=.true., &
! Subprogram not used          VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
! Subprogram not used     call mct_rearr_rearrange(xaop_oe, xaop_o, Re_e2o, sum=.true., &
! Subprogram not used          VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
! Subprogram not used 
! Subprogram not used     !--- normalize by sum of wts associated with mapping
! Subprogram not used 
! Subprogram not used     do n = 1,nloc_a
! Subprogram not used        wt = xaop_a%rAttr(index_sumwt,n)
! Subprogram not used        if (wt /= 0.0_r8) then 
! Subprogram not used           wt = 1.0_r8/wt
! Subprogram not used        else
! Subprogram not used           wt = 1.0_r8
! Subprogram not used        endif
! Subprogram not used        xaop_a%rAttr(:,n) = xaop_a%rAttr(:,n) * wt
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     do n = 1,nloc_o
! Subprogram not used        wt = xaop_o%rAttr(index_sumwt,n)
! Subprogram not used        if (wt /= 0.0_r8) then 
! Subprogram not used           wt = 1.0_r8/wt
! Subprogram not used        else
! Subprogram not used           wt = 1.0_r8
! Subprogram not used        endif
! Subprogram not used        xaop_o%rAttr(:,n) = xaop_o%rAttr(:,n) * wt
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     !--- copy subset of fields to xao_a and xao_o and clean up
! Subprogram not used 
! Subprogram not used     call mct_avect_clean(xaop_ae)
! Subprogram not used     call mct_avect_clean(xaop_oe)
! Subprogram not used 
! Subprogram not used     call mct_avect_copy(xaop_a, xao_a)
! Subprogram not used     call mct_avect_copy(xaop_o, xao_o)
! Subprogram not used 
! Subprogram not used     call mct_avect_clean(xaop_a)
! Subprogram not used     call mct_avect_clean(xaop_o)
! Subprogram not used 
! Subprogram not used   end subroutine seq_flux_atmocnexch_mct

!===============================================================================

  subroutine seq_flux_atmocn_mct(infodata, comp, c2x, grid, xao)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_infodata_type) , intent(in)         :: infodata
    type(component_type)    , intent(in)         :: comp ! ocn or atm
    type(mct_aVect)         , intent(in), target :: c2x  ! a2x_ox or o2x_ax
    character(len=3)        , intent(in)         :: grid 
    type(mct_aVect)         , intent(inout)      :: xao
    !
    ! Local variables
    !
    type(mct_aVect), pointer :: o2x  
    type(mct_avect), pointer :: a2x
    logical     :: flux_albav   ! flux avg option
    logical     :: dead_comps   ! .true.  => dead components are used
    integer(in) :: n,i          ! indices
    real(r8)    :: rlat         ! gridcell latitude in radians
    real(r8)    :: rlon         ! gridcell longitude in radians
    real(r8)    :: cosz         ! Cosine of solar zenith angle
    real(r8)    :: eccen        ! Earth orbit eccentricity
    real(r8)    :: mvelpp       ! Earth orbit
    real(r8)    :: lambm0       ! Earth orbit
    real(r8)    :: obliqr       ! Earth orbit
    real(r8)    :: delta        ! Solar declination angle  in radians
    real(r8)    :: eccf         ! Earth orbit eccentricity factor
    real(r8)    :: calday       ! calendar day including fraction, at 0e
    real(r8)    :: nextsw_cday  ! calendar day of next atm shortwave
    real(r8)    :: anidr        ! albedo: near infrared, direct
    real(r8)    :: avsdr        ! albedo: visible      , direct
    real(r8)    :: anidf        ! albedo: near infrared, diffuse
    real(r8)    :: avsdf        ! albedo: visible      , diffuse
    integer(in) :: nloc         ! number of gridcells
    integer(in) :: ID           ! comm ID
    logical     :: first_call = .true.
    !
    real(R8),parameter :: albdif = 0.06_R8 ! 60 deg reference albedo, diffuse
    real(R8),parameter :: albdir = 0.07_R8 ! 60 deg reference albedo, direct 
    character(*),parameter :: subName =   '(seq_flux_atmocn_mct) '
    !
    !-----------------------------------------------------------------------

    if (trim(grid) == 'ocn') then  
       o2x => component_get_c2x_cx(comp)
       a2x => c2x
    end if

    if (trim(grid) == 'atm') then  
       o2x => c2x
       a2x => component_get_c2x_cx(comp)
    end if

    call seq_infodata_getData(infodata , &
         flux_albav=flux_albav, &
         dead_comps=dead_comps) 

    if (first_call) then
       index_xao_So_tref   = mct_aVect_indexRA(xao,'So_tref')
       index_xao_So_qref   = mct_aVect_indexRA(xao,'So_qref')
       index_xao_So_ustar  = mct_aVect_indexRA(xao,'So_ustar')  
       index_xao_So_re     = mct_aVect_indexRA(xao,'So_re')  
       index_xao_So_ssq    = mct_aVect_indexRA(xao,'So_ssq')
       index_xao_So_u10    = mct_aVect_indexRA(xao,'So_u10')
       index_xao_So_duu10n = mct_aVect_indexRA(xao,'So_duu10n')
       index_xao_Faox_taux = mct_aVect_indexRA(xao,'Faox_taux')
       index_xao_Faox_tauy = mct_aVect_indexRA(xao,'Faox_tauy')  
       index_xao_Faox_lat  = mct_aVect_indexRA(xao,'Faox_lat')   
       index_xao_Faox_sen  = mct_aVect_indexRA(xao,'Faox_sen')   
       index_xao_Faox_evap = mct_aVect_indexRA(xao,'Faox_evap')   
       index_xao_Faox_lwup = mct_aVect_indexRA(xao,'Faox_lwup')  
       
       index_a2x_Sa_z      = mct_aVect_indexRA(a2x,'Sa_z')
       index_a2x_Sa_u      = mct_aVect_indexRA(a2x,'Sa_u')
       index_a2x_Sa_v      = mct_aVect_indexRA(a2x,'Sa_v')
       index_a2x_Sa_tbot   = mct_aVect_indexRA(a2x,'Sa_tbot')
       index_a2x_Sa_ptem   = mct_aVect_indexRA(a2x,'Sa_ptem')
       index_a2x_Sa_shum   = mct_aVect_indexRA(a2x,'Sa_shum')
       index_a2x_Sa_dens   = mct_aVect_indexRA(a2x,'Sa_dens')
       
       index_o2x_So_t      = mct_aVect_indexRA(o2x,'So_t')
       index_o2x_So_u      = mct_aVect_indexRA(o2x,'So_u')
       index_o2x_So_v      = mct_aVect_indexRA(o2x,'So_v')
       first_call = .false.
    end if
       
    if (trim(fluxsetting) /= trim(fluxsetting_atmocn)) then
       call shr_sys_abort(trim(subname)//' ERROR with init')
    endif

    nloc = mct_aVect_lsize(xao)

    ! Update ocean surface fluxes 
    ! Must fabricate "reasonable" data (using dead components)

    emask = mask
    if (dead_comps) then
       do n = 1,nloc
          mask(n) =   1      ! ocn domain mask            ~ 0 <=> inactive cell
          tocn(n) = 290.0_R8 ! ocn temperature            ~ Kelvin
          uocn(n) =   0.0_R8 ! ocn velocity, zonal        ~ m/s
          vocn(n) =   0.0_R8 ! ocn velocity, meridional   ~ m/s
          zbot(n) =  55.0_R8 ! atm height of bottom layer ~ m
          ubot(n) =   0.0_R8 ! atm velocity, zonal        ~ m/s
          vbot(n) =   2.0_R8 ! atm velocity, meridional   ~ m/s
          thbot(n)= 301.0_R8 ! atm potential temperature  ~ Kelvin
          shum(n) = 1.e-2_R8 ! atm specific humidity      ~ kg/kg
          dens(n) =   1.0_R8 ! atm density                ~ kg/m^3
          tbot(n) = 300.0_R8 ! atm temperature            ~ Kelvin
       enddo
    else	
       do n = 1,nloc
          if (mask(n) /= 0) then	
             zbot(n) = a2x%rAttr(index_a2x_Sa_z   ,n)
             ubot(n) = a2x%rAttr(index_a2x_Sa_u   ,n)
             vbot(n) = a2x%rAttr(index_a2x_Sa_v   ,n)
             thbot(n)= a2x%rAttr(index_a2x_Sa_ptem,n)
             shum(n) = a2x%rAttr(index_a2x_Sa_shum,n)
             dens(n) = a2x%rAttr(index_a2x_Sa_dens,n)
             tbot(n) = a2x%rAttr(index_a2x_Sa_tbot,n)
             tocn(n) = o2x%rAttr(index_o2x_So_t   ,n)   
             uocn(n) = o2x%rAttr(index_o2x_So_u   ,n)
             vocn(n) = o2x%rAttr(index_o2x_So_v   ,n)
             !--- mask missing atm or ocn data
             if (dens(n) < 1.0e-12 .or. tocn(n) < 1.0) then
                emask(n) = 0
                !write(logunit,*) 'aoflux tcx1',n,dens(n),tocn(n)
             endif
          end if
       enddo
    end if

    call shr_flux_atmocn (nloc , zbot , ubot, vbot, thbot, &
                          shum , dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          evap , taux , tauy, tref, qref , &
                          !missval should not be needed if flux calc 
                          !consistent with mrgx2a fraction
                          !duu10n,ustar, re  , ssq, missval = 0.0_r8 )
                          duu10n,ustar, re  , ssq)

    do n = 1,nloc
       if (mask(n) /= 0) then	
          xao%rAttr(index_xao_Faox_sen ,n) = sen(n)
          xao%rAttr(index_xao_Faox_lat ,n) = lat(n)
          xao%rAttr(index_xao_Faox_taux,n) = taux(n)
          xao%rAttr(index_xao_Faox_tauy,n) = tauy(n)
          xao%rAttr(index_xao_Faox_evap,n) = evap(n)
          xao%rAttr(index_xao_So_tref  ,n) = tref(n)
	  xao%rAttr(index_xao_So_qref  ,n) = qref(n)
          xao%rAttr(index_xao_So_ustar ,n) = ustar(n)  ! friction velocity
          xao%rAttr(index_xao_So_re    ,n) = re(n)     ! reynolds number
          xao%rAttr(index_xao_So_ssq   ,n) = ssq(n)    ! s.hum. saturation at Ts
          xao%rAttr(index_xao_Faox_lwup,n) = lwup(n)   
          xao%rAttr(index_xao_So_duu10n,n) = duu10n(n)  
          xao%rAttr(index_xao_So_u10   ,n) = sqrt(duu10n(n))  
       end if
    enddo

  end subroutine seq_flux_atmocn_mct

!===============================================================================

end module seq_flux_mct
