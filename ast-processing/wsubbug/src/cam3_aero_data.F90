module cam3_aero_data
!----------------------------------------------------------------------- 
! 
! Purposes: 
!       read, store, interpolate, and return fields
!         of aerosols to 1.  The initialization
!         file (mass.nc) is assumed to be a monthly climatology
!         of aerosols from MATCH (on a sigma pressure
!         coordinate system).
!       also provide a "background" aerosol field to correct
!         for any deficiencies in the physical parameterizations
!         This fields is a "tuning" parameter.
!       Public methods:
!       (1) - initialization
!          read aerosol masses from external file
!             also pressure coordinates
!          convert from monthly average values to mid-month values
!       (2) - interpolation (time and vertical)
!          interpolate onto pressure levels of 1
!          interpolate to time step of 1
!          return mass of aerosols 
!
!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_scam_mod,   only: shr_scam_GetCloseLatLon
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use phys_grid,      only: get_ncols_p, scatter_field_to_chunk
  use time_manager,   only: get_curr_calday
  use infnan,         only: nan, assignment(=)
  use abortutils,     only: endrun
  use scamMod,        only: scmlon,scmlat,single_column
  use error_messages, only: handle_ncerr
  use physics_types,  only: physics_state
  use boundarydata,   only: boundarydata_init, boundarydata_type
  use perf_mod,       only: t_startf, t_stopf
  use cam_logfile,    only: iulog
  use netcdf

  implicit none
  private
  save

  public :: &
     cam3_aero_data_readnl,       & ! read namelist
     cam3_aero_data_register,     & ! register these aerosols with pbuf2d
     cam3_aero_data_init,         & ! read from file, interpolate onto horiz grid
     cam3_aero_data_timestep_init   ! update data-aerosols to this timestep

  ! namelist variables
  logical, public :: cam3_aero_data_on = .false.
  character(len=256) :: bndtvaer = 'bndtvaer'   ! full pathname for time-variant aerosol mass climatology dataset

  ! naer is number of species in climatology
  integer, parameter :: naer = 11

  real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode

  ! indices to aerosol array (species portion)
  integer, parameter :: &
      idxSUL   =  1, &
      idxSSLTA =  2, & ! accumulation mode
      idxSSLTC =  3, & ! coarse mode
      idxOCPHO =  8, &
      idxBCPHO =  9, &
      idxOCPHI =  10, &
      idxBCPHI = 11

  ! indices to sections of array that represent 
  ! groups of aerosols
  integer, parameter :: &
      idxSSLTfirst    = 2, numSSLT  = 2, &
      idxDUSTfirst    = 4, &
      numDUST         = 4, &
      idxCARBONfirst = 8, &
      numCARBON      = 4

  ! names of aerosols are they are represented in
  ! the climatology file.
  ! Appended '_V' indicates field has been vertically summed.
  character(len=8), parameter :: aerosol_name(naer) =  &
     (/"MSUL_V  "&
      ,"MSSLTA_V"&
      ,"MSSLTC_V"&
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"/)

  ! number of different "groups" of aerosols
  integer, parameter :: num_aer_groups=4

  ! which group does each bin belong to?
  integer, dimension(naer), parameter ::  &
      group =(/1,2,2,3,3,3,3,4,4,4,4/)

  ! name of each group
  character(len=10), dimension(num_aer_groups), parameter :: &
      aerosol_names = (/'sul  ','sslt ','dust ','car  '/)

  ! this boundarydata_type is used for datasets in the ncols format only.
  type(boundarydata_type) :: aerosol_datan

  integer :: aernid = -1           ! netcdf id for aerosol file (init to invalid)
  integer :: species_id(naer) = -1 ! netcdf_id of each aerosol species (init to invalid)
  integer :: Mpsid                 ! netcdf id for MATCH PS
  integer :: nm = 1                ! index to prv month in array. init to 1 and toggle between 1 and 2
  integer :: np = 2                ! index to nxt month in array. init to 2 and toggle between 1 and 2
  integer :: mo_nxt = huge(1)      ! index to nxt month in file

  real(r8) :: cdaym                ! calendar day of prv month
  real(r8) :: cdayp                ! calendar day of next month

  ! aerosol mass 
  real(r8), allocatable :: aer_mass(:, :, :, :)

  ! Days into year for mid month date
  ! This variable is dumb, the dates are in the dataset to be read in but they are
  ! slightly different than this so getting rid of it causes a change which 
  ! exceeds roundoff.
  real(r8) :: Mid(12) = (/16.5_r8,  46.0_r8,  75.5_r8, 106.0_r8, 136.5_r8, 167.0_r8, &
                         197.5_r8, 228.5_r8, 259.0_r8, 289.5_r8, 320.0_r8, 350.5_r8 /)
  
  !  values read from file and temporary values used for interpolation
  !
  !  aerosolc is:
  !  Cumulative Mass at midpoint of each month
  !    on 1's horizontal grid (col)
  !    on MATCH's levels (lev)
  !  aerosolc
  integer, parameter :: paerlev = 28       ! number of levels for aerosol fields (MUST = naerlev)
  integer :: naerlev                       ! size of level dimension in MATCH data
  integer :: naerlon
  integer :: naerlat
  real(r8), pointer :: M_hybi(:)           ! MATCH hybi
  real(r8), pointer :: M_ps(:,:)           ! surface pressure from MATCH file
  real(r8), pointer :: aerosolc(:,:,:,:,:) ! Aerosol cumulative mass from MATCH
  real(r8), pointer :: M_ps_cam_col(:,:,:) ! PS from MATCH on Cam Columns

  ! indices for fields in the physics buffer
  integer :: cam3_sul_idx, cam3_ssam_idx, cam3_sscm_idx, &
      cam3_dust1_idx, cam3_dust2_idx, cam3_dust3_idx, cam3_dust4_idx,&
      cam3_ocpho_idx, cam3_bcpho_idx, cam3_ocphi_idx, cam3_bcphi_idx

!================================================================================================
contains
!================================================================================================

subroutine cam3_aero_data_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cam3_aero_data_readnl'

   namelist /cam3_aero_data_nl/ cam3_aero_data_on, bndtvaer
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam3_aero_data_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cam3_aero_data_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if


   ! Broadcast namelist variables
   call mpibcast(cam3_aero_data_on, 1, mpilog, 0, mpicom)
   call mpibcast(bndtvaer, len(bndtvaer), mpichar, 0, mpicom)


   ! Prevent using these before they are set.
   cdaym = nan
   cdayp = nan

end subroutine cam3_aero_data_readnl

!================================================================================================

! Subprogram not used subroutine cam3_aero_data_register
! Subprogram not used 
! Subprogram not used    ! register old prescribed aerosols with physics buffer
! Subprogram not used    
! Subprogram not used    use physics_buffer, only: pbuf_add_field, dtype_r8
! Subprogram not used 
! Subprogram not used    call pbuf_add_field('cam3_sul',  'physpkg',dtype_r8,(/pcols,pver/),cam3_sul_idx)
! Subprogram not used    call pbuf_add_field('cam3_ssam', 'physpkg',dtype_r8,(/pcols,pver/),cam3_ssam_idx)
! Subprogram not used    call pbuf_add_field('cam3_sscm', 'physpkg',dtype_r8,(/pcols,pver/),cam3_sscm_idx)
! Subprogram not used    call pbuf_add_field('cam3_dust1','physpkg',dtype_r8,(/pcols,pver/),cam3_dust1_idx)
! Subprogram not used    call pbuf_add_field('cam3_dust2','physpkg',dtype_r8,(/pcols,pver/),cam3_dust2_idx)
! Subprogram not used    call pbuf_add_field('cam3_dust3','physpkg',dtype_r8,(/pcols,pver/),cam3_dust3_idx)
! Subprogram not used    call pbuf_add_field('cam3_dust4','physpkg',dtype_r8,(/pcols,pver/),cam3_dust4_idx)
! Subprogram not used    call pbuf_add_field('cam3_ocpho','physpkg',dtype_r8,(/pcols,pver/),cam3_ocpho_idx)
! Subprogram not used    call pbuf_add_field('cam3_bcpho','physpkg',dtype_r8,(/pcols,pver/),cam3_bcpho_idx)
! Subprogram not used    call pbuf_add_field('cam3_ocphi','physpkg',dtype_r8,(/pcols,pver/),cam3_ocphi_idx)
! Subprogram not used    call pbuf_add_field('cam3_bcphi','physpkg',dtype_r8,(/pcols,pver/),cam3_bcphi_idx)
! Subprogram not used 
! Subprogram not used end subroutine cam3_aero_data_register

!================================================================================================

! Subprogram not used subroutine cam3_aero_data_init(phys_state)
! Subprogram not used !------------------------------------------------------------------
! Subprogram not used !  Reads in:
! Subprogram not used !     file from which to read aerosol Masses on 1 grid. Currently
! Subprogram not used !        assumed to be MATCH ncep runs, averaged by month.
! Subprogram not used !     NOTE (Data have been externally interpolated onto 1 grid 
! Subprogram not used !        and backsolved to provide Mid-month values)
! Subprogram not used !     
! Subprogram not used !  Populates:
! Subprogram not used !     module variables:
! Subprogram not used !       aerosolc(pcols,paerlev+1,begchunk:endchunk,naer,2))
! Subprogram not used !       aerosolc(  column_index
! Subprogram not used !                , level_index (match levels)
! Subprogram not used !                , chunk_index 
! Subprogram not used !                , species_index
! Subprogram not used !                , month = 1:2 )
! Subprogram not used !       M_hybi(level_index = Lev_MATCH) = pressure at mid-level.
! Subprogram not used !       M_ps_cam_col(column,chunk,month) ! PS from MATCH on Cam Columns
! Subprogram not used !
! Subprogram not used !  Method:
! Subprogram not used !    read data from file
! Subprogram not used !    allocate memory for storage of aerosol data on 1 horizontal grid
! Subprogram not used !    distribute data to remote nodes
! Subprogram not used !    populates the module variables
! Subprogram not used !
! Subprogram not used !------------------------------------------------------------------
! Subprogram not used    use ioFileMod,    only: getfil
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    use mpishorthand
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
! Subprogram not used 
! Subprogram not used ! local variables
! Subprogram not used 
! Subprogram not used    integer :: naerlev
! Subprogram not used 
! Subprogram not used    integer dateid                       ! netcdf id for date variable
! Subprogram not used    integer secid                        ! netcdf id for seconds variable
! Subprogram not used    integer londimid                     ! netcdf id for longitude dimension
! Subprogram not used    integer latdimid                     ! netcdf id for latitude dimension
! Subprogram not used    integer levdimid                     ! netcdf id for level dimension
! Subprogram not used 
! Subprogram not used    integer timesiz                      ! number of time samples (=12) in netcdf file
! Subprogram not used    integer latid                        ! netcdf id for latitude variable
! Subprogram not used    integer Mhybiid                      ! netcdf id for MATCH hybi
! Subprogram not used    integer timeid                       ! netcdf id for time variable
! Subprogram not used    integer dimids(nf90_max_var_dims)      ! variable shape
! Subprogram not used    integer :: start(4)                  ! start vector for netcdf calls
! Subprogram not used    integer :: kount(4)                  ! count vector for netcdf calls
! Subprogram not used    integer mo                           ! month index
! Subprogram not used    integer m                            ! constituent index
! Subprogram not used    integer :: n                         ! loop index
! Subprogram not used    integer :: i,j,k                     ! spatial indices
! Subprogram not used    integer :: date_aer(12)              ! Date on aerosol dataset (YYYYMMDD)
! Subprogram not used    integer :: attnum                    ! attribute number
! Subprogram not used    integer :: ierr                      ! netcdf return code
! Subprogram not used    real(r8) ::  coldata(paerlev)    ! aerosol field read in from dataset
! Subprogram not used    integer :: ret
! Subprogram not used    integer mo_prv                       ! index to previous month
! Subprogram not used    integer latidx,lonidx
! Subprogram not used 
! Subprogram not used    character(len=8) :: aname                   ! temporary aerosol name
! Subprogram not used    character(len=8) :: tmp_aero_name(naer) ! name for input to boundary data
! Subprogram not used 
! Subprogram not used    character(len=256) :: locfn          ! netcdf local filename to open
! Subprogram not used !
! Subprogram not used ! aerosol_data will be read in from the aerosol boundary dataset, then scattered to chunks
! Subprogram not used ! after filling in the bottom level with zeros
! Subprogram not used ! 
! Subprogram not used    real(r8), allocatable :: aerosol_data(:,:,:)    ! aerosol field read in from dataset
! Subprogram not used    real(r8), allocatable :: aerosol_field(:,:,:)   ! (plon,paerlev+1,plat)  aerosol field to be scattered
! Subprogram not used    real(r8) :: caldayloc                           ! calendar day of current timestep
! Subprogram not used    real(r8) :: closelat,closelon
! Subprogram not used 
! Subprogram not used    character(len=*), parameter :: subname = 'cam3_aero_data_init'
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call t_startf(subname)
! Subprogram not used 
! Subprogram not used    allocate (aer_mass(pcols, pver, naer, begchunk:endchunk) )
! Subprogram not used 
! Subprogram not used    ! set new aerosol names because input file has 1 seasalt bin
! Subprogram not used    do m = 1, naer
! Subprogram not used       tmp_aero_name(m)=aerosol_name(m)
! Subprogram not used       if (aerosol_name(m)=='MSSLTA_V') tmp_aero_name(m) = 'MSSLT_V'
! Subprogram not used       if (aerosol_name(m)=='MSSLTC_V') tmp_aero_name(m) = 'MSSLT_V'
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    allocate (aerosolc(pcols,paerlev+1,begchunk:endchunk,naer,2))
! Subprogram not used    aerosolc(:,:,:,:,:) = 0._r8
! Subprogram not used 
! Subprogram not used    caldayloc = get_curr_calday ()
! Subprogram not used    
! Subprogram not used    if (caldayloc < Mid(1)) then
! Subprogram not used       mo_prv = 12
! Subprogram not used       mo_nxt =  1
! Subprogram not used    else if (caldayloc >= Mid(12)) then
! Subprogram not used       mo_prv = 12
! Subprogram not used       mo_nxt =  1
! Subprogram not used    else
! Subprogram not used       do i = 2 , 12
! Subprogram not used          if (caldayloc < Mid(i)) then
! Subprogram not used             mo_prv = i-1
! Subprogram not used             mo_nxt = i
! Subprogram not used             exit
! Subprogram not used          end if
! Subprogram not used       end do
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    ! Set initial calendar day values
! Subprogram not used    cdaym = Mid(mo_prv)
! Subprogram not used    cdayp = Mid(mo_nxt)
! Subprogram not used 
! Subprogram not used    if (masterproc) &
! Subprogram not used       write(iulog,*) subname//': CAM3 prescribed aerosol dataset is: ', trim(bndtvaer)
! Subprogram not used 
! Subprogram not used    call getfil (bndtvaer, locfn, 0)
! Subprogram not used 
! Subprogram not used    call handle_ncerr( nf90_open (locfn, 0, aernid),&
! Subprogram not used       subname, 339)
! Subprogram not used 
! Subprogram not used    if (single_column) &
! Subprogram not used       call shr_scam_GetCloseLatLon(aernid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
! Subprogram not used 
! Subprogram not used    ! Check to see if this dataset is in ncol format. 
! Subprogram not used    aerosol_datan%isncol=.false.
! Subprogram not used    ierr = nf90_inq_dimid( aernid,  'ncol', londimid )
! Subprogram not used    if ( ierr==NF90_NOERR ) then
! Subprogram not used 
! Subprogram not used       aerosol_datan%isncol=.true.
! Subprogram not used       call handle_ncerr(nf90_close(aernid),subname, 350)
! Subprogram not used 
! Subprogram not used       call boundarydata_init(bndtvaer, phys_state, tmp_aero_name, naer, &
! Subprogram not used                              aerosol_datan, 3)
! Subprogram not used 
! Subprogram not used       aerosolc(:,1:paerlev,:,:,:)=aerosol_datan%fields
! Subprogram not used 
! Subprogram not used       M_ps_cam_col=>aerosol_datan%ps
! Subprogram not used       M_hybi=>aerosol_datan%hybi
! Subprogram not used 
! Subprogram not used    else 
! Subprogram not used 
! Subprogram not used       ! Allocate memory for dynamic arrays local to this module
! Subprogram not used       allocate (M_ps_cam_col(pcols,begchunk:endchunk,2))
! Subprogram not used       allocate (M_hybi(paerlev+1))
! Subprogram not used       ! TBH:  HACK to avoid use of uninitialized values when ncols < pcols
! Subprogram not used       M_ps_cam_col(:,:,:) = 0._r8
! Subprogram not used 
! Subprogram not used       if (masterproc) then
! Subprogram not used 
! Subprogram not used          ! First ensure dataset is 1-ready
! Subprogram not used 
! Subprogram not used          call handle_ncerr(nf90_inquire_attribute (aernid, nf90_global, 'cam-ready', attnum=attnum),&
! Subprogram not used               subname//': interpaerosols needs to be run to create a cam-ready aerosol dataset')
! Subprogram not used 
! Subprogram not used          ! Get and check dimension info
! Subprogram not used 
! Subprogram not used          call handle_ncerr( nf90_inq_dimid( aernid,  'lon', londimid ),&
! Subprogram not used               subname, 378)
! Subprogram not used          call handle_ncerr( nf90_inq_dimid( aernid,  'lev', levdimid ),&
! Subprogram not used               subname, 380)
! Subprogram not used          call handle_ncerr( nf90_inq_dimid( aernid, 'time',   timeid ),&
! Subprogram not used               subname, 382)
! Subprogram not used          call handle_ncerr( nf90_inq_dimid( aernid,  'lat', latdimid ),&
! Subprogram not used               subname, 384)
! Subprogram not used          call handle_ncerr( nf90_inquire_dimension( aernid, londimid, len=naerlon ),&
! Subprogram not used               subname, 386)
! Subprogram not used          call handle_ncerr( nf90_inquire_dimension( aernid, levdimid, len=naerlev ),&
! Subprogram not used               subname, 388)
! Subprogram not used          call handle_ncerr( nf90_inquire_dimension( aernid, latdimid, len=naerlat ),&
! Subprogram not used               subname, 390)
! Subprogram not used          call handle_ncerr( nf90_inquire_dimension( aernid,   timeid, len=timesiz ),&
! Subprogram not used               subname, 392)
! Subprogram not used 
! Subprogram not used          call handle_ncerr( nf90_inq_varid( aernid, 'date',   dateid ),&
! Subprogram not used               subname, 395)
! Subprogram not used          call handle_ncerr( nf90_inq_varid( aernid, 'datesec', secid ),&
! Subprogram not used               subname, 397)
! Subprogram not used 
! Subprogram not used          do m = 1, naer
! Subprogram not used             aname=aerosol_name(m)
! Subprogram not used             ! rename because file has only one seasalt field
! Subprogram not used             if (aname=='MSSLTA_V') aname = 'MSSLT_V'
! Subprogram not used             if (aname=='MSSLTC_V') aname = 'MSSLT_V'
! Subprogram not used             call handle_ncerr( nf90_inq_varid( aernid, TRIM(aname), species_id(m)), &
! Subprogram not used                subname, 405)
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          call handle_ncerr( nf90_inq_varid( aernid, 'lat', latid   ),&
! Subprogram not used               subname, 409)
! Subprogram not used 
! Subprogram not used          ! quick sanity check on one field
! Subprogram not used          call handle_ncerr( nf90_inquire_variable (aernid, species_id(1), dimids=dimids),&
! Subprogram not used               subname, 413)
! Subprogram not used 
! Subprogram not used          if ( (dimids(4) /= timeid) .or. &
! Subprogram not used               (dimids(3) /= levdimid) .or. &
! Subprogram not used               (dimids(2) /= latdimid) .or. &
! Subprogram not used               (dimids(1) /= londimid) ) then
! Subprogram not used             write(iulog,*) subname//': Data must be ordered time, lev, lat, lon'
! Subprogram not used             write(iulog,*) 'data are       ordered as', dimids(4), dimids(3), dimids(2), dimids(1)
! Subprogram not used             write(iulog,*) 'data should be ordered as', timeid, levdimid, latdimid, londimid
! Subprogram not used             call endrun ()
! Subprogram not used          end if
! Subprogram not used 
! Subprogram not used          ! use hybi,PS from MATCH
! Subprogram not used          call handle_ncerr( nf90_inq_varid( aernid, 'hybi', Mhybiid   ),&
! Subprogram not used               subname, 427)
! Subprogram not used          call handle_ncerr( nf90_inq_varid( aernid, 'PS', Mpsid   ),&
! Subprogram not used               subname, 429)
! Subprogram not used 
! Subprogram not used          ! check dimension order for MATCH's surface pressure
! Subprogram not used          call handle_ncerr( nf90_inquire_variable (aernid, Mpsid, dimids=dimids),&
! Subprogram not used               subname, 433)
! Subprogram not used          if ( (dimids(3) /= timeid) .or. &
! Subprogram not used               (dimids(2) /= latdimid) .or. &
! Subprogram not used               (dimids(1) /= londimid) ) then
! Subprogram not used             write(iulog,*) subname//': Pressure must be ordered time, lat, lon'
! Subprogram not used             write(iulog,*) 'data are       ordered as', dimids(3), dimids(2), dimids(1)
! Subprogram not used             write(iulog,*) 'data should be ordered as', timeid, levdimid, latdimid, londimid
! Subprogram not used             call endrun ()
! Subprogram not used          end if
! Subprogram not used 
! Subprogram not used          ! read in hybi from MATCH
! Subprogram not used          call handle_ncerr( nf90_get_var (aernid, Mhybiid, M_hybi),&
! Subprogram not used               subname, 445)
! Subprogram not used 
! Subprogram not used          ! Retrieve date and sec variables.
! Subprogram not used          call handle_ncerr( nf90_get_var (aernid, dateid, date_aer),&
! Subprogram not used               subname, 449)
! Subprogram not used          if (timesiz < 12) then
! Subprogram not used             write(iulog,*) subname//': When cycling aerosols, dataset must have 12 consecutive ', &
! Subprogram not used                  'months of data starting with Jan'
! Subprogram not used             write(iulog,*) 'Current dataset has only ',timesiz,' months'
! Subprogram not used             call endrun ()
! Subprogram not used          end if
! Subprogram not used          do mo = 1,12
! Subprogram not used             if (mod(date_aer(mo),10000)/100 /= mo) then
! Subprogram not used                write(iulog,*) subname//': When cycling aerosols, dataset must have 12 consecutive ', &
! Subprogram not used                     'months of data starting with Jan'
! Subprogram not used                write(iulog,*)'Month ',mo,' of dataset says date=',date_aer(mo)
! Subprogram not used                call endrun ()
! Subprogram not used             end if
! Subprogram not used          end do
! Subprogram not used          if (single_column) then
! Subprogram not used             naerlat=1
! Subprogram not used             naerlon=1
! Subprogram not used          endif
! Subprogram not used          kount(:) = (/naerlon,naerlat,paerlev,1/)
! Subprogram not used       end if          ! masterproc
! Subprogram not used 
! Subprogram not used       ! broadcast hybi to nodes
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call mpibcast (M_hybi, paerlev+1, mpir8, 0, mpicom)
! Subprogram not used       call mpibcast (kount, 3, mpiint, 0, mpicom)
! Subprogram not used       naerlon = kount(1)
! Subprogram not used       naerlat = kount(2)
! Subprogram not used 
! Subprogram not used       allocate(aerosol_field(kount(1),kount(3)+1,kount(2)))
! Subprogram not used       allocate(M_ps(kount(1),kount(2)))
! Subprogram not used       if (masterproc) allocate(aerosol_data(kount(1),kount(2),kount(3)))
! Subprogram not used 
! Subprogram not used       ! Retrieve Aerosol Masses (kg/m^2 in each layer), transpose to model order (lon,lev,lat),
! Subprogram not used       ! then scatter to slaves.
! Subprogram not used       if (nm /= 1 .or. np /= 2) call endrun (subname//': bad nm or np value')
! Subprogram not used       do n=nm,np
! Subprogram not used          if (n == 1) then
! Subprogram not used             mo = mo_prv
! Subprogram not used          else
! Subprogram not used             mo = mo_nxt
! Subprogram not used          end if
! Subprogram not used          
! Subprogram not used          do m=1,naer
! Subprogram not used             if (masterproc) then
! Subprogram not used                if (single_column) then
! Subprogram not used                   start(:) = (/lonidx,latidx,1,mo/)
! Subprogram not used                else
! Subprogram not used                   start(:) = (/1,1,1,mo/)
! Subprogram not used                endif
! Subprogram not used                kount(:) = (/naerlon,naerlat,paerlev,1/)
! Subprogram not used 
! Subprogram not used                call handle_ncerr( nf90_get_var (aernid, species_id(m),aerosol_data, start, kount),&
! Subprogram not used                     subname, 503)
! Subprogram not used                do j=1,naerlat
! Subprogram not used                   do k=1,paerlev
! Subprogram not used                      aerosol_field(:,k,j) = aerosol_data(:,j,k)
! Subprogram not used                   end do
! Subprogram not used                   aerosol_field(:,paerlev+1,j) = 0._r8   ! value at bottom
! Subprogram not used                end do
! Subprogram not used                
! Subprogram not used             end if
! Subprogram not used             call scatter_field_to_chunk (1, paerlev+1, 1, naerlon, aerosol_field, &
! Subprogram not used                  aerosolc(:,:,:,m,n))
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          ! Retrieve PS from Match
! Subprogram not used 
! Subprogram not used          if (masterproc) then
! Subprogram not used             if (single_column) then
! Subprogram not used                start(:) = (/lonidx,latidx,mo,-1/)
! Subprogram not used             else
! Subprogram not used                start(:) = (/1,1,mo,-1/)
! Subprogram not used             endif
! Subprogram not used             kount(:) = (/naerlon,naerlat,1,-1/)
! Subprogram not used             call handle_ncerr( nf90_get_var(aernid, Mpsid, M_ps,start,kount),&
! Subprogram not used                  subname, 526)
! Subprogram not used          end if
! Subprogram not used          call scatter_field_to_chunk (1, 1, 1, naerlon, M_ps(:,:), M_ps_cam_col(:,:,n))
! Subprogram not used       end do     ! n=nm,np (=1,2)
! Subprogram not used 
! Subprogram not used       if(masterproc) deallocate(aerosol_data)
! Subprogram not used       deallocate(aerosol_field)
! Subprogram not used 
! Subprogram not used    end if   ! Check to see if this dataset is in ncol format. 
! Subprogram not used 
! Subprogram not used    call t_stopf(subname)
! Subprogram not used 
! Subprogram not used end subroutine cam3_aero_data_init

!================================================================================================

! Subprogram not used subroutine cam3_aero_data_timestep_init(pbuf2d,  phys_state)
! Subprogram not used !------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  Input:
! Subprogram not used !     time at which aerosol masses are needed (get_curr_calday())
! Subprogram not used !     chunk index
! Subprogram not used !     1's vertical grid (pint)
! Subprogram not used !
! Subprogram not used !  Output:
! Subprogram not used !     values for Aerosol Mass at time specified by get_curr_calday
! Subprogram not used !     on vertical grid specified by pint (aer_mass) :: aerosol at time t
! Subprogram not used !
! Subprogram not used !  Method:
! Subprogram not used !     first determine which indexs of aerosols are the bounding data sets
! Subprogram not used !     interpolate both onto vertical grid aerm(),aerp().
! Subprogram not used !     from those two, interpolate in time.
! Subprogram not used !
! Subprogram not used !------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    use interpolate_data, only: get_timeinterp_factors
! Subprogram not used    
! Subprogram not used    use physics_buffer, only: physics_buffer_desc, dtype_r8, pbuf_set_field, pbuf_get_chunk
! Subprogram not used    use cam_logfile,     only: iulog
! Subprogram not used    use ppgrid,          only: begchunk,endchunk
! Subprogram not used    use physconst,       only: gravit
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! aerosol fields interpolated to current time step
! Subprogram not used !   on pressure levels of this time step.
! Subprogram not used ! these should be made read-only for other modules
! Subprogram not used ! Is allocation done correctly here?
! Subprogram not used !
! Subprogram not used    
! Subprogram not used    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
! Subprogram not used    type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)
! Subprogram not used    real(r8) :: pint(pcols,pverp)  ! interface pres.
! Subprogram not used    integer :: c                           ! chunk index
! Subprogram not used    real(r8) caldayloc                     ! calendar day of current timestep
! Subprogram not used    real(r8) fact1, fact2                  ! time interpolation factors
! Subprogram not used 
! Subprogram not used    integer i, k, j                        ! spatial indices
! Subprogram not used    integer m                              ! constituent index
! Subprogram not used    integer lats(pcols),lons(pcols)        ! latitude and longitudes of column
! Subprogram not used    integer ncol                           ! number of columns
! Subprogram not used    integer lchnk                          ! chunk index
! Subprogram not used    
! Subprogram not used    real(r8) speciesmin(naer)              ! minimal value for each species
! Subprogram not used !
! Subprogram not used ! values before current time step "the minus month"
! Subprogram not used ! aerosolm(pcols,pver) is value of preceeding month's aerosol masses
! Subprogram not used ! aerosolp(pcols,pver) is value of next month's aerosol masses
! Subprogram not used !  (think minus and plus or values to left and right of point to be interpolated)
! Subprogram not used !
! Subprogram not used    real(r8) aerosolm(pcols,pver,naer,begchunk:endchunk) ! aerosol mass from MATCH in column,level at previous (minus) month
! Subprogram not used !
! Subprogram not used ! values beyond (or at) current time step "the plus month"
! Subprogram not used !
! Subprogram not used    real(r8) aerosolp(pcols,pver,naer,begchunk:endchunk) ! aerosol mass from MATCH in column,level at next (plus) month 
! Subprogram not used    real(r8) :: mass_to_mmr(pcols,pver)
! Subprogram not used 
! Subprogram not used    character(len=*), parameter :: subname = 'cam3_aero_data_timestep_init'
! Subprogram not used 
! Subprogram not used    logical error_found
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call aerint(phys_state)
! Subprogram not used 
! Subprogram not used    caldayloc = get_curr_calday ()
! Subprogram not used 
! Subprogram not used    ! Determine time interpolation factors.  1st arg says we are cycling 1 year of data
! Subprogram not used    call get_timeinterp_factors (.true., mo_nxt, cdaym, cdayp, caldayloc, &
! Subprogram not used                     fact1, fact2, 'GET_AEROSOL:')
! Subprogram not used 
! Subprogram not used    ! interpolate (prv and nxt month) bounding datasets onto cam vertical grid.
! Subprogram not used    ! compute mass mixing ratios on CAMS's pressure coordinate
! Subprogram not used    !  for both the "minus" and "plus" months
! Subprogram not used    !
! Subprogram not used    !  This loop over chunk could probably be removed by working with the whole
! Subprogram not used    !  begchunk:endchunk group at once.  It would require a slight generalization 
! Subprogram not used    !  in vert_interpolate.
! Subprogram not used    do c = begchunk,endchunk  
! Subprogram not used                                 
! Subprogram not used       lchnk = phys_state(c)%lchnk
! Subprogram not used       pint = phys_state(c)%pint
! Subprogram not used       ncol = get_ncols_p(c)
! Subprogram not used 
! Subprogram not used       call vert_interpolate (M_ps_cam_col(:,c,nm), pint, nm, aerosolm(:,:,:,c), ncol, c)
! Subprogram not used       call vert_interpolate (M_ps_cam_col(:,c,np), pint, np, aerosolp(:,:,:,c), ncol, c)
! Subprogram not used 
! Subprogram not used       ! Time interpolate.
! Subprogram not used       do m=1,naer
! Subprogram not used          do k=1,pver
! Subprogram not used             do i=1,ncol
! Subprogram not used                aer_mass(i,k,m,c) = aerosolm(i,k,m,c)*fact1 + aerosolp(i,k,m,c)*fact2
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used          ! Partition seasalt aerosol mass
! Subprogram not used          if (m .eq. idxSSLTA) then
! Subprogram not used             aer_mass(:ncol,:,m,c) = (1._r8-wgt_sscm)*aer_mass(:ncol,:,m,c) ! fraction of seasalt mass in accumulation mode
! Subprogram not used          elseif (m .eq. idxSSLTC) then
! Subprogram not used             aer_mass(:ncol,:,m,c) = wgt_sscm*aer_mass(:ncol,:,m,c)      ! fraction of seasalt mass in coarse mode
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       ! exit if mass is negative (we have previously set
! Subprogram not used       !  cumulative mass to be a decreasing function.)
! Subprogram not used       speciesmin(:) = 0._r8 ! speciesmin(m) = 0 is minimum mass for each species
! Subprogram not used  
! Subprogram not used       error_found = .false.
! Subprogram not used       do m=1,naer
! Subprogram not used          do k=1,pver
! Subprogram not used             do i=1,ncol
! Subprogram not used                if (aer_mass(i, k, m,c) < speciesmin(m)) error_found = .true.
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used       if (error_found) then
! Subprogram not used          do m=1,naer
! Subprogram not used             do k=1,pver
! Subprogram not used                do i=1,ncol
! Subprogram not used                   if (aer_mass(i, k, m,c) < speciesmin(m)) then
! Subprogram not used                      write(iulog,*) subname//': negative mass mixing ratio, exiting'
! Subprogram not used                      write(iulog,*) 'm, column, pver',m, i, k ,aer_mass(i, k, m,c)
! Subprogram not used                      call endrun ()
! Subprogram not used                   end if
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end if
! Subprogram not used       do k = 1, pver
! Subprogram not used          mass_to_mmr(1:ncol,k) = gravit/(pint(1:ncol,k+1)-pint(1:ncol,k))
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
! Subprogram not used 
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_sul_idx,   aer_mass(1:ncol,:,        idxSUL,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_ssam_idx,  aer_mass(1:ncol,:,      idxSSLTA,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_sscm_idx,  aer_mass(1:ncol,:,      idxSSLTC,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_dust1_idx, aer_mass(1:ncol,:,  idxDUSTfirst,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_dust2_idx, aer_mass(1:ncol,:,idxDUSTfirst+1,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_dust3_idx, aer_mass(1:ncol,:,idxDUSTfirst+2,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_dust4_idx, aer_mass(1:ncol,:,idxDUSTfirst+3,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_ocpho_idx, aer_mass(1:ncol,:,      idxOCPHO,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_bcpho_idx, aer_mass(1:ncol,:,      idxBCPHO,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_ocphi_idx, aer_mass(1:ncol,:,      idxOCPHI,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used       call pbuf_set_field(phys_buffer_chunk, cam3_bcphi_idx, aer_mass(1:ncol,:,      idxBCPHI,c)*mass_to_mmr(:ncol,:), &
! Subprogram not used            start=(/1,1/), kount=(/ncol,pver/))
! Subprogram not used 
! Subprogram not used    enddo ! c = begchunk:endchunk
! Subprogram not used 
! Subprogram not used end subroutine cam3_aero_data_timestep_init

!================================================================================================

! Subprogram not used subroutine vert_interpolate (Match_ps, pint, n, aerosol_mass, ncol, c)
! Subprogram not used !--------------------------------------------------------------------
! Subprogram not used ! Input: match surface pressure, cam interface pressure, 
! Subprogram not used !        month index, number of columns, chunk index
! Subprogram not used ! 
! Subprogram not used ! Output: Aerosol mass mixing ratio (aerosol_mass)
! Subprogram not used !
! Subprogram not used ! Method:
! Subprogram not used !         interpolate column mass (cumulative) from match onto
! Subprogram not used !           cam's vertical grid (pressure coordinate)
! Subprogram not used !         convert back to mass mixing ratio
! Subprogram not used !
! Subprogram not used !--------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real(r8), intent(out) :: aerosol_mass(pcols,pver,naer)  ! aerosol mass from MATCH
! Subprogram not used    real(r8), intent(in) :: Match_ps(pcols)                ! surface pressure at a particular month
! Subprogram not used    real(r8), intent(in) :: pint(pcols,pverp)              ! interface pressure from 1
! Subprogram not used 
! Subprogram not used    integer, intent(in) :: ncol,c                          ! chunk index and number of columns
! Subprogram not used    integer, intent(in) :: n                               ! prv or nxt month index
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used    integer m                           ! index to aerosol species
! Subprogram not used    integer kupper(pcols)               ! last upper bound for interpolation
! Subprogram not used    integer i, k, kk, kkstart, kount    ! loop vars for interpolation
! Subprogram not used    integer isv, ksv, msv               ! loop indices to save
! Subprogram not used 
! Subprogram not used    logical bad                         ! indicates a bad point found
! Subprogram not used    logical lev_interp_comp             ! interpolation completed for a level 
! Subprogram not used    logical error_found
! Subprogram not used 
! Subprogram not used    real(r8) aerosol(pcols,pverp,naer)  ! cumulative mass of aerosol in column beneath upper 
! Subprogram not used                                        ! interface of level in column at particular month
! Subprogram not used    real(r8) dpl, dpu                   ! lower and upper intepolation factors
! Subprogram not used    real(r8) v_coord                    ! vertical coordinate
! Subprogram not used    real(r8) AER_diff                   ! temp var for difference between aerosol masses
! Subprogram not used 
! Subprogram not used    character(len=*), parameter :: subname = 'cam3_aero_data.vert_interpolate'
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call t_startf ('vert_interpolate')
! Subprogram not used !
! Subprogram not used ! Initialize index array 
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       kupper(i) = 1
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! assign total mass to topmost level
! Subprogram not used !
! Subprogram not used    aerosol(:,1,:) = aerosolc(:,1,c,:,n)
! Subprogram not used !
! Subprogram not used ! At every pressure level, interpolate onto that pressure level
! Subprogram not used !
! Subprogram not used    do k=2,pver
! Subprogram not used !
! Subprogram not used ! Top level we need to start looking is the top level for the previous k
! Subprogram not used ! for all longitude points
! Subprogram not used !
! Subprogram not used       kkstart = paerlev+1
! Subprogram not used       do i=1,ncol
! Subprogram not used          kkstart = min0(kkstart,kupper(i))
! Subprogram not used       end do
! Subprogram not used       kount = 0
! Subprogram not used !
! Subprogram not used ! Store level indices for interpolation
! Subprogram not used !
! Subprogram not used ! for the pressure interpolation should be comparing
! Subprogram not used ! pint(column,lev) with M_hybi(lev)*M_ps_cam_col(month,column,chunk)
! Subprogram not used !
! Subprogram not used       lev_interp_comp = .false.
! Subprogram not used       do kk=kkstart,paerlev
! Subprogram not used          if(.not.lev_interp_comp) then
! Subprogram not used          do i=1,ncol
! Subprogram not used             v_coord = pint(i,k)
! Subprogram not used             if (M_hybi(kk)*Match_ps(i) .lt. v_coord .and. v_coord .le. M_hybi(kk+1)*Match_ps(i)) then
! Subprogram not used                kupper(i) = kk
! Subprogram not used                kount = kount + 1
! Subprogram not used             end if
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! If all indices for this level have been found, do the interpolation and
! Subprogram not used ! go to the next level
! Subprogram not used !
! Subprogram not used ! Interpolate in pressure.
! Subprogram not used !
! Subprogram not used          if (kount.eq.ncol) then
! Subprogram not used             do m=1,naer
! Subprogram not used                do i=1,ncol
! Subprogram not used                   dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
! Subprogram not used                   dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
! Subprogram not used                   aerosol(i,k,m) = &
! Subprogram not used                      (aerosolc(i,kupper(i)  ,c,m,n)*dpl + &
! Subprogram not used                      aerosolc(i,kupper(i)+1,c,m,n)*dpu)/(dpl + dpu)
! Subprogram not used                enddo !i
! Subprogram not used             end do
! Subprogram not used             lev_interp_comp = .true.
! Subprogram not used          end if
! Subprogram not used          end if
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used ! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! Subprogram not used ! must extrapolate from the bottom or top pressure level for at least some
! Subprogram not used ! of the longitude points.
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       if(.not.lev_interp_comp) then
! Subprogram not used          do m=1,naer
! Subprogram not used             do i=1,ncol
! Subprogram not used                if (pint(i,k) .lt. M_hybi(1)*Match_ps(i)) then
! Subprogram not used                   aerosol(i,k,m) =  aerosolc(i,1,c,m,n)
! Subprogram not used                else if (pint(i,k) .gt. M_hybi(paerlev+1)*Match_ps(i)) then
! Subprogram not used                   aerosol(i,k,m) = 0.0_r8
! Subprogram not used                else
! Subprogram not used                   dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
! Subprogram not used                   dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
! Subprogram not used                   aerosol(i,k,m) = &
! Subprogram not used                      (aerosolc(i,kupper(i)  ,c,m,n)*dpl + &
! Subprogram not used                      aerosolc(i,kupper(i)+1,c,m,n)*dpu)/(dpl + dpu)
! Subprogram not used                end if
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          if (kount.gt.ncol) then
! Subprogram not used             call endrun (subname//': Bad data: non-monotonicity suspected in dependent variable')
! Subprogram not used          end if
! Subprogram not used       end if
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !   call t_startf ('vi_checks')
! Subprogram not used !
! Subprogram not used ! aerosol mass beneath lowest interface (pverp) must be 0
! Subprogram not used !
! Subprogram not used    aerosol(1:ncol,pverp,:) = 0._r8
! Subprogram not used !
! Subprogram not used ! Set mass in layer to zero whenever it is less than 
! Subprogram not used !   1.e-40 kg/m^2 in the layer
! Subprogram not used !
! Subprogram not used    do m = 1, naer
! Subprogram not used       do k = 1, pver
! Subprogram not used          do i = 1, ncol
! Subprogram not used             if (aerosol(i,k,m) < 1.e-40_r8) aerosol(i,k,m) = 0._r8
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Set mass in layer to zero whenever it is less than 
! Subprogram not used !   10^-15 relative to column total mass
! Subprogram not used !
! Subprogram not used    error_found = .false.
! Subprogram not used    do m = 1, naer
! Subprogram not used       do k = 1, pver
! Subprogram not used          do i = 1, ncol
! Subprogram not used             AER_diff = aerosol(i,k,m) - aerosol(i,k+1,m)
! Subprogram not used             if( abs(AER_diff) < 1e-15_r8*aerosol(i,1,m)) then
! Subprogram not used                AER_diff = 0._r8
! Subprogram not used             end if
! Subprogram not used             aerosol_mass(i,k,m)= AER_diff 
! Subprogram not used             if (aerosol_mass(i,k,m) < 0) error_found = .true.
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used    if (error_found) then
! Subprogram not used       do m = 1, naer
! Subprogram not used          do k = 1, pver
! Subprogram not used             do i = 1, ncol
! Subprogram not used                if (aerosol_mass(i,k,m) < 0) then
! Subprogram not used                   write(iulog,*) subname//': mass < 0, m, col, lev, mass',m, i, k, aerosol_mass(i,k,m)
! Subprogram not used                   write(iulog,*) subname//': aerosol(k),(k+1)',aerosol(i,k,m),aerosol(i,k+1,m)
! Subprogram not used                   write(iulog,*) subname//': pint(k+1),(k)',pint(i,k+1),pint(i,k)
! Subprogram not used                   write(iulog,*)'n,c',n,c
! Subprogram not used                   call endrun()
! Subprogram not used                end if
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    call t_stopf ('vert_interpolate')
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used end subroutine vert_interpolate

!================================================================================================

! Subprogram not used subroutine aerint (phys_state)
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
! Subprogram not used 
! Subprogram not used    integer :: ntmp                                ! used in index swapping
! Subprogram not used    integer :: start(4)                            ! start vector for netcdf calls
! Subprogram not used    integer :: kount(4)                            ! count vector for netcdf calls
! Subprogram not used    integer :: i,j,k                               ! spatial indices
! Subprogram not used    integer :: m                                   ! constituent index
! Subprogram not used    integer :: cols, cole
! Subprogram not used    integer :: lchnk, ncol
! Subprogram not used    real(r8) :: caldayloc                          ! calendar day of current timestep
! Subprogram not used    real(r8) :: aerosol_data(naerlon,naerlat,paerlev)    ! aerosol field read in from dataset
! Subprogram not used    real(r8) :: aerosol_field(naerlon,paerlev+1,naerlat) ! aerosol field to be scattered
! Subprogram not used    integer latidx,lonidx
! Subprogram not used    real(r8) closelat,closelon
! Subprogram not used 
! Subprogram not used    character(len=*), parameter :: subname = 'cam3_aero_data.aerint'
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (single_column) &
! Subprogram not used       call shr_scam_GetCloseLatLon(aernid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
! Subprogram not used  
! Subprogram not used !
! Subprogram not used ! determine if need to read in next month data
! Subprogram not used ! also determine time interpolation factors
! Subprogram not used !
! Subprogram not used    caldayloc = get_curr_calday ()  
! Subprogram not used !
! Subprogram not used ! If model time is past current forward timeslice, then
! Subprogram not used ! masterproc reads in the next timeslice for time interpolation.  Messy logic is 
! Subprogram not used ! for interpolation between December and January (mo_nxt == 1).  Just like
! Subprogram not used ! ozone_data_timestep_init, sstint.
! Subprogram not used !
! Subprogram not used    if (caldayloc > cdayp .and. .not. (mo_nxt == 1 .and. caldayloc >= cdaym)) then
! Subprogram not used       mo_nxt = mod(mo_nxt,12) + 1
! Subprogram not used       cdaym = cdayp
! Subprogram not used       cdayp = Mid(mo_nxt)
! Subprogram not used !
! Subprogram not used ! Check for valid date info
! Subprogram not used !
! Subprogram not used       if (.not. (mo_nxt == 1 .or. caldayloc <= cdayp)) then
! Subprogram not used          call endrun (subname//': Non-monotonicity suspected in input aerosol data')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       ntmp = nm
! Subprogram not used       nm = np
! Subprogram not used       np = ntmp
! Subprogram not used 
! Subprogram not used       if(aerosol_datan%isncol) then
! Subprogram not used          do lchnk=begchunk,endchunk
! Subprogram not used             ncol=phys_state(lchnk)%ncol
! Subprogram not used             cols=1
! Subprogram not used             cole=cols+aerosol_datan%count(cols,lchnk)-1
! Subprogram not used             do while(cole<=ncol)
! Subprogram not used                start=(/aerosol_datan%start(cols,lchnk),mo_nxt,1,-1/)
! Subprogram not used                kount=(/aerosol_datan%count(cols,lchnk),1,-1,-1/)
! Subprogram not used                call handle_ncerr( nf90_get_var(aerosol_datan%ncid, aerosol_datan%psid , &
! Subprogram not used                     aerosol_datan%ps(cols:cole,lchnk,np), start(1:2), &
! Subprogram not used                     kount(1:2)),&
! Subprogram not used                     subname, 957)
! Subprogram not used                start(2)=1
! Subprogram not used                start(3)=mo_nxt
! Subprogram not used                kount(2)=paerlev
! Subprogram not used                kount(3)=1
! Subprogram not used                do m=1,naer
! Subprogram not used                   call handle_ncerr( nf90_get_var(aerosol_datan%ncid, aerosol_datan%dataid(m) , &
! Subprogram not used                        aerosol_datan%fields(cols:cole,:,lchnk,m,np),  &
! Subprogram not used                        start(1:3), kount(1:3)),&
! Subprogram not used                        subname, 966)
! Subprogram not used 
! Subprogram not used                end do
! Subprogram not used                if(cols==ncol) exit
! Subprogram not used                cols=cols+aerosol_datan%count(cols,lchnk)
! Subprogram not used                cole=cols+aerosol_datan%count(cols,lchnk)-1
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used          aerosolc(:,1:paerlev,:,:,np)=aerosol_datan%fields(:,:,:,:,np)
! Subprogram not used       else
! Subprogram not used          do m=1,naer
! Subprogram not used             if (masterproc) then
! Subprogram not used                if (single_column) then
! Subprogram not used                   naerlon=1
! Subprogram not used                   naerlat=1
! Subprogram not used                   start(:) = (/lonidx,latidx,1,mo_nxt/)
! Subprogram not used                else
! Subprogram not used                   start(:) = (/1,1,1,mo_nxt/)
! Subprogram not used                endif
! Subprogram not used                kount(:) = (/naerlon,naerlat,paerlev,1/)
! Subprogram not used                call handle_ncerr( nf90_get_var (aernid, species_id(m), aerosol_data, start, kount),&
! Subprogram not used                     subname, 987)
! Subprogram not used 
! Subprogram not used                do j=1,naerlat
! Subprogram not used                   do k=1,paerlev
! Subprogram not used                      aerosol_field(:,k,j) = aerosol_data(:,j,k)
! Subprogram not used                   end do
! Subprogram not used                   aerosol_field(:,paerlev+1,j) = 0._r8   ! value at bottom
! Subprogram not used                end do
! Subprogram not used             end if
! Subprogram not used             call scatter_field_to_chunk (1, paerlev+1, 1, naerlon, aerosol_field, &
! Subprogram not used                  aerosolc(:,:,:,m,np))
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! Retrieve PS from Match
! Subprogram not used !
! Subprogram not used          if (masterproc) then
! Subprogram not used                if (single_column) then
! Subprogram not used                   naerlon=1
! Subprogram not used                   naerlat=1
! Subprogram not used                   start(:) = (/lonidx,latidx,mo_nxt,-1/)
! Subprogram not used                else
! Subprogram not used                   start(:) = (/1,1,mo_nxt,-1/)
! Subprogram not used                endif
! Subprogram not used                kount(:) = (/naerlon,naerlat,1,-1/)
! Subprogram not used                call handle_ncerr( nf90_get_var (aernid, Mpsid, M_ps, start, kount),&
! Subprogram not used                     subname, 1012)
! Subprogram not used                write(iulog,*) subname//': Read aerosols data for julian day', Mid(mo_nxt)
! Subprogram not used             end if
! Subprogram not used             call scatter_field_to_chunk (1, 1, 1, naerlon, M_ps(:,:), M_ps_cam_col(:,:,np))
! Subprogram not used          end if
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used end subroutine aerint

end module cam3_aero_data
