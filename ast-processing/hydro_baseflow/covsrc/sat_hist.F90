!-------------------------------------------------------------------------------
! Outputs history field columns as specified by a satellite track data file
!
! Created by Francis Vitt -- 17 Sep 2010
!-------------------------------------------------------------------------------
module sat_hist

  use perf_mod,      only: t_startf, t_stopf
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use cam_logfile,   only: iulog
  use ppgrid,        only: pcols, pver, begchunk, endchunk
  use cam_history_support, only: fieldname_lenp2, max_string_len, ptapes
  use spmd_utils,    only: masterproc, iam
  use abortutils,    only: endrun

  use pio,           only: file_desc_t, iosystem_desc_t, iosystem_desc_t, var_desc_t, io_desc_t
  use pio,           only: pio_openfile, pio_redef, pio_enddef, pio_inq_dimid, pio_inq_varid, pio_seterrorhandling, pio_def_var
  use pio,           only: pio_inq_dimlen, pio_get_att, pio_put_att, pio_get_var, pio_put_var, pio_write_darray
  use pio,           only: pio_real, pio_int, pio_double, pio_copy_att
  use pio,           only: PIO_WRITE,PIO_NOWRITE, PIO_NOERR, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, PIO_Rearr_box, PIO_GLOBAL
  use spmd_utils,    only: mpicom

  use mpishorthand,  only: mpichar, mpiint

   use physconst, only: pi 
  
  implicit none

  private
  save

  public :: sat_hist_readnl
  public :: sat_hist_init
  public :: sat_hist_write
  public :: sat_hist_define
  public :: is_satfile

  character(len=max_string_len)  :: sathist_track_infile
  type(file_desc_t) :: infile

  integer :: half_step
  logical :: has_sat_hist = .false.

  integer :: sathist_nclosest
  integer :: sathist_ntimestep

  real(r8), allocatable :: obs_lats(:)
  real(r8), allocatable :: obs_lons(:)

  logical  :: doy_format
  real(r8) :: first_datetime
  real(r8) :: last_datetime
  integer  :: last_start_index
  integer  :: time_ndx
  integer  :: t_buffer_size
  integer, allocatable :: date_buffer(:), time_buffer(:)
  integer :: sat_tape_num=ptapes-1

  
  ! input file
  integer :: n_profiles
  integer :: time_vid, date_vid, lat_vid, lon_vid, instr_vid, orbit_vid, prof_vid, zenith_vid

  integer :: in_julian_vid
  integer :: in_localtime_vid
  integer :: in_doy_vid
  integer :: in_occ_type_vid

  integer :: in_start_col


  ! output file
  type(var_desc_t) :: out_latid, out_lonid, out_dstid, out_instrid, out_zenithid, out_orbid, out_profid
  type(var_desc_t) :: out_instr_lat_vid, out_instr_lon_vid
  type(var_desc_t) :: out_obs_date_vid, out_obs_time_vid
  type(var_desc_t) :: out_julian_vid
  type(var_desc_t) :: out_localtime_vid
  type(var_desc_t) :: out_doy_vid
  type(var_desc_t) :: out_occ_type_vid

  logical, parameter :: debug = .false.

  real(r8), parameter :: rad2deg = 180._r8/pi            ! degrees per radian

contains
  
!-------------------------------------------------------------------------------

  logical function is_satfile (file_index)
    integer, intent(in) :: file_index ! index of file in question
    is_satfile = file_index == sat_tape_num
  end function is_satfile

!-------------------------------------------------------------------------------
  subroutine sat_hist_readnl(nlfile, hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape)
    
    use namelist_utils,      only: find_group_name
    use units,               only: getunit, freeunit
    use cam_history_support, only: pflds
    use cam_instance,        only: inst_suffix

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    character(len=*), intent(inout) :: hfilename_spec(:)
    character(len=*), intent(inout) :: fincl(:,:)
    character(len=1), intent(inout) :: avgflag_pertape(:)
    integer,          intent(inout) :: mfilt(:), nhtfrq(:)
    
    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'sat_hist_readnl'
    integer :: f, fcnt

    character(len=fieldname_lenp2) :: sathist_fincl(pflds)
    character(len=max_string_len)  :: sathist_hfilename_spec
    integer :: sathist_mfilt, sat_tape_num

    namelist /satellite_options_nl/ sathist_track_infile, sathist_hfilename_spec, sathist_fincl, &
         sathist_mfilt, sathist_nclosest, sathist_ntimestep

    ! set defaults

    sathist_track_infile = ' '
    sathist_hfilename_spec = '%c.cam' // trim(inst_suffix) // '.hs.%y-%m-%d-%s.nc'
    sathist_fincl(:) = ' '
    sathist_mfilt = 100000
    sathist_nclosest = 1
    sathist_ntimestep = 1

    !read namelist options

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'satellite_options_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, satellite_options_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if


    ! broadcast the options to all MPI tasks
    call mpibcast(sathist_track_infile,   len(sathist_track_infile),   mpichar, 0, mpicom)
    call mpibcast(sathist_hfilename_spec, len(sathist_hfilename_spec), mpichar, 0, mpicom)
    call mpibcast(sathist_fincl,          pflds*len(sathist_fincl(1)), mpichar, 0, mpicom)
    call mpibcast(sathist_mfilt,          1,                           mpiint,  0, mpicom)
    call mpibcast(sathist_nclosest,       1,                           mpiint,  0, mpicom)
    call mpibcast(sathist_ntimestep,      1,                           mpiint,  0, mpicom)


    has_sat_hist = len_trim(sathist_track_infile) > 0

    if (.not.has_sat_hist) return

     sat_tape_num=ptapes-1
     hfilename_spec(sat_tape_num) = sathist_hfilename_spec
     mfilt(sat_tape_num) = sathist_mfilt
     fcnt=0
     do f=1, pflds
        fincl(f,sat_tape_num) = sathist_fincl(f)
        if(len_trim(sathist_fincl(f)) > 0) then
           fcnt=fcnt+1
        end if
     enddo
     
     nhtfrq(sat_tape_num) = 1
     avgflag_pertape(sat_tape_num) = 'I'

     if(masterproc) then
        write(iulog,*) 'sathist_track_infile: ',trim(sathist_track_infile)
        write(iulog,*) 'sathist_hfilename_spec: ',trim(sathist_hfilename_spec)
        write(iulog,*) 'sathist_fincl: ',(trim(sathist_fincl(f))//' ', f=1,fcnt)
        write(iulog,*) 'max columns per file sathist_mfilt: ',sathist_mfilt
        write(iulog,*) 'sathist_nclosest: ',sathist_nclosest
        write(iulog,*) 'sathist_ntimestep: ',sathist_ntimestep
     end if

   end subroutine sat_hist_readnl

  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine sat_hist_init
    use cam_pio_utils, only: cam_pio_openfile
    use ioFileMod,     only: getfil
    use spmd_utils,    only: npes
    use time_manager,  only: get_step_size
    use string_utils,  only: to_lower, GLC

    implicit none

    character(len=max_string_len)  :: locfn       ! Local filename
    integer :: ierr, dimid, i

    character(len=128) :: date_format

    if (.not.has_sat_hist) return

    call getfil (sathist_track_infile, locfn)
    call cam_pio_openfile(infile, locfn, PIO_NOWRITE)

    ierr = pio_inq_dimid(infile,'profs',dimid)
    ierr = pio_inq_dimlen(infile, dimid, n_profiles)

    ierr = pio_inq_varid( infile, 'time', time_vid )
    ierr = pio_inq_varid( infile, 'date', date_vid )

    ierr = pio_get_att( infile, date_vid, 'long_name', date_format)
    date_format = to_lower(trim( date_format(:GLC(date_format))))

    if ( index( date_format, 'yyyymmdd') > 0 ) then
       doy_format = .false.
    else if  ( index( date_format, 'yyyyddd') > 0 ) then
       doy_format = .true.
    else
       call endrun('sat_hist_init: date_format not recognized : '//trim(date_format))
    endif

    ierr = pio_inq_varid( infile, 'lat', lat_vid )
    ierr = pio_inq_varid( infile, 'lon', lon_vid )

    call pio_seterrorhandling(infile, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( infile, 'instr_num', instr_vid )
    if(ierr/=PIO_NOERR) instr_vid=-1

    ierr = pio_inq_varid( infile, 'orbit_num', orbit_vid )
    if(ierr/=PIO_NOERR) orbit_vid=-1

    ierr = pio_inq_varid( infile, 'prof_num',  prof_vid )
    if(ierr/=PIO_NOERR) prof_vid=-1

    ierr = pio_inq_varid( infile, 'instr_sza', zenith_vid )
    if(ierr/=PIO_NOERR) zenith_vid=-1

    ierr = pio_inq_varid( infile, 'julian', in_julian_vid )
    if(ierr/=PIO_NOERR) in_julian_vid=-1

    ierr = pio_inq_varid( infile, 'local_time', in_localtime_vid )
    if(ierr/=PIO_NOERR) in_localtime_vid=-1

    ierr = pio_inq_varid( infile, 'doy', in_doy_vid )
    if(ierr/=PIO_NOERR) in_doy_vid=-1

    ierr = pio_inq_varid( infile, 'occ_type', in_occ_type_vid )
    if(ierr/=PIO_NOERR) in_occ_type_vid=-1

    call pio_seterrorhandling(infile, PIO_INTERNAL_ERROR)

    call read_datetime( first_datetime, 1 )
    call read_datetime( last_datetime, n_profiles )
    last_start_index = -1
    t_buffer_size = min(1000,n_profiles)
    allocate( date_buffer(t_buffer_size), time_buffer(t_buffer_size) )
    if (masterproc) write(iulog,*) "sathist_init:", n_profiles, first_datetime, last_datetime
    if ( last_datetime<first_datetime ) then
       call endrun('sat_hist_init: satellite track file has invalid date time info')
    endif

    time_ndx = 1
    half_step = get_step_size()*0.5_r8

  end subroutine sat_hist_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Subprogram not used   subroutine read_datetime( datetime, index )
! Subprogram not used 
! Subprogram not used     real(r8), intent( out ) :: datetime
! Subprogram not used     integer,  intent( in )  :: index
! Subprogram not used 
! Subprogram not used     integer :: ierr
! Subprogram not used     integer :: cnt(1)
! Subprogram not used     integer :: start(1)
! Subprogram not used     integer :: date(1), time(1)
! Subprogram not used 
! Subprogram not used     cnt = (/ 1 /)
! Subprogram not used     start = (/index/)
! Subprogram not used 
! Subprogram not used     ierr = pio_get_var( infile, time_vid, start, cnt, time )
! Subprogram not used     ierr = pio_get_var( infile, date_vid, start, cnt, date )
! Subprogram not used     
! Subprogram not used     datetime = convert_date_time( date(1),time(1) )
! Subprogram not used 
! Subprogram not used   end subroutine read_datetime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Subprogram not used   subroutine read_buffered_datetime( datetime, index )
! Subprogram not used 
! Subprogram not used     real(r8), intent( out ) :: datetime
! Subprogram not used     integer,  intent( in )  :: index
! Subprogram not used 
! Subprogram not used     integer :: ii
! Subprogram not used 
! Subprogram not used     integer :: ierr
! Subprogram not used     integer :: cnt
! Subprogram not used     integer :: start
! Subprogram not used     integer :: date, time
! Subprogram not used     
! Subprogram not used     ! If the request is outside of the buffer then reload the buffer.
! Subprogram not used     if ((last_start_index == -1) .or. (index < last_start_index) &
! Subprogram not used          .or. (index >= (last_start_index + t_buffer_size))) then
! Subprogram not used 
! Subprogram not used        start = (index - 1) / t_buffer_size * t_buffer_size + 1
! Subprogram not used        if ( start+t_buffer_size-1 <= n_profiles ) then
! Subprogram not used           cnt = t_buffer_size 
! Subprogram not used        else
! Subprogram not used           cnt = n_profiles-start+1
! Subprogram not used        endif
! Subprogram not used        ierr = pio_get_var( infile, time_vid, (/ start /), (/ cnt /), time_buffer(1:cnt) )
! Subprogram not used        ierr = pio_get_var( infile, date_vid, (/ start /), (/ cnt /), date_buffer(1:cnt) )
! Subprogram not used 
! Subprogram not used        last_start_index = start
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ii = mod( index - 1, t_buffer_size ) + 1
! Subprogram not used     time = time_buffer(ii)
! Subprogram not used     date = date_buffer(ii)
! Subprogram not used     datetime = convert_date_time( date,time )
! Subprogram not used 
! Subprogram not used   end subroutine read_buffered_datetime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Subprogram not used   function convert_date_time( date,time )
! Subprogram not used     use time_manager, only: set_time_float_from_date
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: date,time
! Subprogram not used     real(r8) :: convert_date_time
! Subprogram not used 
! Subprogram not used     real(r8) :: datetime
! Subprogram not used     integer :: yr, doy, mon, dom
! Subprogram not used 
! Subprogram not used     if ( doy_format ) then
! Subprogram not used        yr = date/1000
! Subprogram not used        doy = date - yr*1000
! Subprogram not used        call set_time_float_from_date( datetime, yr, 1, doy, time )
! Subprogram not used     else 
! Subprogram not used        yr = date/10000
! Subprogram not used        mon = (date - yr*10000)/100
! Subprogram not used        dom = date - yr*10000 - mon*100
! Subprogram not used        call set_time_float_from_date( datetime, yr, mon, dom, time )
! Subprogram not used     endif
! Subprogram not used     convert_date_time = datetime
! Subprogram not used 
! Subprogram not used   end function convert_date_time
!-------------------------------------------------------------------------------
! Subprogram not used   subroutine sat_hist_define(outfile)
! Subprogram not used     use pio, only : pio_inquire
! Subprogram not used     type(file_desc_t), intent(inout) :: outfile
! Subprogram not used 
! Subprogram not used     integer :: coldim
! Subprogram not used     integer :: ierr
! Subprogram not used     
! Subprogram not used     ierr = pio_inquire(outfile, unlimitedDimId=coldim)
! Subprogram not used 
! Subprogram not used     call pio_seterrorhandling(outfile, PIO_BCAST_ERROR)
! Subprogram not used     ierr = define_var( 'instr_lat', coldim, infile, lat_vid,  outfile, out_instr_lat_vid )
! Subprogram not used     ierr = define_var( 'instr_lon', coldim, infile, lon_vid,  outfile, out_instr_lon_vid )
! Subprogram not used     ierr = define_var( 'obs_time', coldim, infile, time_vid,  outfile, out_obs_time_vid )
! Subprogram not used     ierr = define_var( 'obs_date', coldim, infile, date_vid,  outfile, out_obs_date_vid )
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_varid( outfile, 'distance', out_dstid )
! Subprogram not used     if (ierr /= PIO_NOERR) then
! Subprogram not used        ierr = pio_def_var  ( outfile, 'distance', PIO_REAL, (/coldim/), out_dstid )
! Subprogram not used        ierr = pio_put_att  ( outfile, out_dstid, "long_name", "distance from midpoint to observation")
! Subprogram not used        ierr = pio_put_att  ( outfile, out_dstid, "units", "km")
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (orbit_vid>0) then
! Subprogram not used        ierr = define_var( 'orbit_num', coldim, infile, orbit_vid,  outfile, out_orbid )
! Subprogram not used     endif
! Subprogram not used     if (prof_vid>0) then
! Subprogram not used        ierr = define_var( 'prof_num', coldim, infile, prof_vid,  outfile, out_profid )
! Subprogram not used     endif
! Subprogram not used     if (instr_vid>0) then
! Subprogram not used        ierr = define_var( 'instr_num', coldim, infile, instr_vid,  outfile, out_instrid )
! Subprogram not used     endif
! Subprogram not used     if (zenith_vid>0) then
! Subprogram not used        ierr = define_var( 'instr_sza', coldim, infile, zenith_vid,  outfile, out_zenithid )
! Subprogram not used     endif
! Subprogram not used     if (in_occ_type_vid>0) then
! Subprogram not used        ierr = define_var( 'occ_type', coldim, infile, in_occ_type_vid,  outfile, out_occ_type_vid )
! Subprogram not used     endif
! Subprogram not used     if (in_julian_vid>0) then
! Subprogram not used        ierr = define_var( 'julian', coldim, infile, in_julian_vid,  outfile, out_julian_vid )
! Subprogram not used     endif
! Subprogram not used     if (in_localtime_vid>0) then
! Subprogram not used        ierr = define_var( 'local_time', coldim, infile, in_localtime_vid,  outfile, out_localtime_vid )
! Subprogram not used     endif
! Subprogram not used     if (in_doy_vid>0) then
! Subprogram not used        ierr = define_var( 'doy', coldim, infile, in_doy_vid,  outfile, out_doy_vid )
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call pio_seterrorhandling(outfile, PIO_INTERNAL_ERROR)
! Subprogram not used     ierr=pio_put_att (outfile, PIO_GLOBAL, 'satellite_track_file', sathist_track_infile)
! Subprogram not used   end subroutine sat_hist_define


!-------------------------------------------------------------------------------
! Subprogram not used   subroutine sat_hist_write( tape , nflds, nfils)
! Subprogram not used 
! Subprogram not used     use ppgrid,   only : pcols, begchunk, endchunk
! Subprogram not used     use dyn_grid, only : get_dyn_grid_parm
! Subprogram not used     use cam_pio_utils, only: phys_decomp, dyn_decomp
! Subprogram not used     use cam_history_support, only : active_entry
! Subprogram not used     use pio, only : pio_file_is_open
! Subprogram not used     implicit none
! Subprogram not used     type(active_entry) :: tape
! Subprogram not used     integer, intent(in) :: nflds
! Subprogram not used     integer, intent(inout) :: nfils
! Subprogram not used 
! Subprogram not used     integer :: t, f, i, ncols, nocols    
! Subprogram not used     integer :: ierr
! Subprogram not used 
! Subprogram not used     integer, allocatable :: col_ndxs(:)
! Subprogram not used     integer, allocatable :: chk_ndxs(:)
! Subprogram not used     integer, allocatable :: fdyn_ndxs(:)
! Subprogram not used     integer, allocatable :: ldyn_ndxs(:)
! Subprogram not used     integer, allocatable :: phs_owners(:)
! Subprogram not used     integer, allocatable :: dyn_owners(:)
! Subprogram not used     real(r8),allocatable :: mlats(:)
! Subprogram not used     real(r8),allocatable :: mlons(:)
! Subprogram not used     real(r8),allocatable :: phs_dists(:)
! Subprogram not used 
! Subprogram not used     integer :: coldim
! Subprogram not used 
! Subprogram not used     integer :: io_type
! Subprogram not used     logical :: has_dyn_flds
! Subprogram not used 
! Subprogram not used     if (.not.has_sat_hist) return
! Subprogram not used 
! Subprogram not used     call read_next_position( ncols )
! Subprogram not used 
! Subprogram not used     if ( ncols < 1 ) return
! Subprogram not used 
! Subprogram not used     call t_startf ('sat_hist_write')
! Subprogram not used 
! Subprogram not used     ! The n closest columns to the observation will be output,
! Subprogram not used     ! so increase the size of the columns used for output/
! Subprogram not used     nocols = ncols * sathist_nclosest
! Subprogram not used 
! Subprogram not used     allocate( col_ndxs(nocols) )
! Subprogram not used     allocate( chk_ndxs(nocols) )
! Subprogram not used     allocate( fdyn_ndxs(nocols) )
! Subprogram not used     allocate( ldyn_ndxs(nocols) )
! Subprogram not used     allocate( phs_owners(nocols) )
! Subprogram not used     allocate( dyn_owners(nocols) )
! Subprogram not used     allocate( mlats(nocols) )
! Subprogram not used     allocate( mlons(nocols) )
! Subprogram not used     allocate( phs_dists(nocols) )
! Subprogram not used 
! Subprogram not used     has_dyn_flds = .false.
! Subprogram not used     dyn_flds_loop: do f=1,nflds
! Subprogram not used        if ( tape%hlist(f)%field%decomp_type == dyn_decomp ) then
! Subprogram not used           has_dyn_flds = .true.
! Subprogram not used           exit dyn_flds_loop
! Subprogram not used        endif
! Subprogram not used     enddo dyn_flds_loop
! Subprogram not used 
! Subprogram not used     call get_indices( obs_lats, obs_lons, ncols, nocols, has_dyn_flds, col_ndxs, chk_ndxs, &
! Subprogram not used          fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners, mlats, mlons, phs_dists )
! Subprogram not used 
! Subprogram not used     if ( .not. pio_file_is_open(tape%File) ) then
! Subprogram not used        call endrun('sat file not open')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_dimid(tape%File,'ncol',coldim )
! Subprogram not used     
! Subprogram not used     ierr = pio_inq_varid(tape%File, 'lat', out_latid )
! Subprogram not used     ierr = pio_inq_varid(tape%File, 'lon', out_lonid )
! Subprogram not used     ierr = pio_inq_varid(tape%File, 'distance', out_dstid )
! Subprogram not used 
! Subprogram not used     call write_record_coord( tape, mlats(:), mlons(:), phs_dists(:), ncols, nfils )
! Subprogram not used 
! Subprogram not used     do f=1,nflds
! Subprogram not used 
! Subprogram not used        select case (tape%hlist(f)%field%decomp_type)
! Subprogram not used        case (phys_decomp)
! Subprogram not used           call dump_columns(tape%File, tape%hlist(f), nocols, nfils, col_ndxs(:), chk_ndxs(:), phs_owners(:) )
! Subprogram not used        case (dyn_decomp)
! Subprogram not used           call dump_columns(tape%File, tape%hlist(f), nocols, nfils, fdyn_ndxs(:), ldyn_ndxs(:), dyn_owners(:) )
! Subprogram not used        end select
! Subprogram not used 
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     deallocate( col_ndxs, chk_ndxs, fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners )
! Subprogram not used     deallocate( mlons, mlats, phs_dists )
! Subprogram not used     deallocate( obs_lons, obs_lats )
! Subprogram not used 
! Subprogram not used     nfils = nfils + nocols
! Subprogram not used 
! Subprogram not used     call t_stopf ('sat_hist_write')
! Subprogram not used 
! Subprogram not used   end subroutine sat_hist_write

!-------------------------------------------------------------------------------
! Subprogram not used   subroutine dump_columns( File, hitem, ncols, nfils, fdims, ldims, owners  )
! Subprogram not used     use cam_history_support,  only: field_info, hentry, hist_coords
! Subprogram not used     use pionfwrite_mod, only: write_nf
! Subprogram not used     use cam_pio_utils, only : fillvalue
! Subprogram not used     use pio,            only: pio_initdecomp, pio_freedecomp, pio_setframe, pio_offset, pio_iam_iotask, pio_setdebuglevel
! Subprogram not used 
! Subprogram not used     type(File_desc_t),intent(inout)  :: File
! Subprogram not used     type(hentry),     intent(in), target     :: hitem
! Subprogram not used     integer,          intent(in)     :: ncols
! Subprogram not used     integer,          intent(in)     :: nfils
! Subprogram not used     integer,          intent(in)     :: fdims(:)
! Subprogram not used     integer,          intent(in)     :: ldims(:)
! Subprogram not used     integer,          intent(in)     :: owners(:)
! Subprogram not used 
! Subprogram not used     type(field_info), pointer :: field
! Subprogram not used     type(var_desc_t) :: vardesc
! Subprogram not used     type(iosystem_desc_t), pointer :: sat_iosystem
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used     integer :: t, ierr, ndims
! Subprogram not used     integer, allocatable :: dimlens(:)
! Subprogram not used 
! Subprogram not used     real(r8), allocatable :: buf(:)
! Subprogram not used     integer,  allocatable :: dof(:)
! Subprogram not used     integer :: i,k, cnt
! Subprogram not used 
! Subprogram not used     call t_startf ('sat_hist::dump_columns')
! Subprogram not used 
! Subprogram not used     sat_iosystem => File%iosystem
! Subprogram not used     field => hitem%field
! Subprogram not used     vardesc = hitem%varid(1)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ndims=1
! Subprogram not used     if(associated(field%mdims)) then
! Subprogram not used        ndims = size(field%mdims)+1
! Subprogram not used     else if(field%numlev>1) then
! Subprogram not used        ndims=2
! Subprogram not used     end if
! Subprogram not used     allocate(dimlens(ndims))
! Subprogram not used     dimlens(ndims)=ncols
! Subprogram not used     if(ndims>2) then
! Subprogram not used        do i=1,ndims-1
! Subprogram not used           dimlens(i)=hist_coords(field%mdims(i))%dimsize
! Subprogram not used        enddo
! Subprogram not used     else if(field%numlev>1) then
! Subprogram not used        dimlens(1) = field%numlev
! Subprogram not used     end if
! Subprogram not used    
! Subprogram not used 
! Subprogram not used     allocate( buf( product(dimlens) ) )
! Subprogram not used     allocate( dof( product(dimlens) ) )
! Subprogram not used 
! Subprogram not used     cnt = 0
! Subprogram not used     buf = fillvalue
! Subprogram not used     dof = 0
! Subprogram not used 
! Subprogram not used     do i = 1,ncols
! Subprogram not used        do k = 1,field%numlev
! Subprogram not used           cnt = cnt+1
! Subprogram not used           if ( iam == owners(i) ) then
! Subprogram not used              buf(cnt) = hitem%hbuf( fdims(i), k, ldims(i) )
! Subprogram not used              dof(cnt) = cnt
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call pio_setframe(vardesc, int(-1,kind=PIO_OFFSET))
! Subprogram not used 
! Subprogram not used     call pio_initdecomp(sat_iosystem, pio_double, dimlens, dof, iodesc )
! Subprogram not used     if(pio_iam_iotask(sat_iosystem)) &
! Subprogram not used          iodesc%start(ndims)=iodesc%start(ndims)+nfils-1
! Subprogram not used 
! Subprogram not used     call pio_write_darray(File, vardesc, iodesc, buf, ierr, fillval=fillvalue)
! Subprogram not used 
! Subprogram not used     call pio_freedecomp(sat_iosystem, iodesc)
! Subprogram not used 
! Subprogram not used     deallocate( buf )
! Subprogram not used     deallocate( dof )
! Subprogram not used     deallocate( dimlens )
! Subprogram not used 
! Subprogram not used     call t_stopf ('sat_hist::dump_columns')
! Subprogram not used 
! Subprogram not used   end subroutine dump_columns

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Subprogram not used   subroutine read_next_position( ncols )
! Subprogram not used     use time_manager, only: get_curr_date, get_prev_date
! Subprogram not used     use time_manager, only: set_time_float_from_date
! Subprogram not used    
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer,  intent(out) :: ncols
! Subprogram not used 
! Subprogram not used     integer :: ierr
! Subprogram not used     integer :: yr, mon, day, tod
! Subprogram not used     real(r8) :: begdatetime, enddatetime
! Subprogram not used     integer :: beg_ndx, end_ndx, i
! Subprogram not used 
! Subprogram not used     real(r8) :: datetime
! Subprogram not used 
! Subprogram not used     call get_curr_date(yr, mon, day, tod)
! Subprogram not used     call set_time_float_from_date(begdatetime, yr, mon, day, tod-half_step*sathist_ntimestep)
! Subprogram not used     call set_time_float_from_date(enddatetime, yr, mon, day, tod+half_step*sathist_ntimestep)
! Subprogram not used 
! Subprogram not used     ncols = 0
! Subprogram not used 
! Subprogram not used     if ( first_datetime > enddatetime ) then
! Subprogram not used        if (masterproc) write(iulog,'(a,2f16.6)') &
! Subprogram not used             'sat_hist->read_next_position: all of the satellite date times are after the time window', first_datetime, enddatetime
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used     if ( last_datetime < begdatetime ) then
! Subprogram not used        if (masterproc) write(iulog,'(a,2f16.6)') &
! Subprogram not used             'sat_hist->read_next_position: all of the satellite date times are before the time window', begdatetime, last_datetime
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call t_startf ('sat_hist::read_next_position')
! Subprogram not used 
! Subprogram not used     beg_ndx = -99
! Subprogram not used     end_ndx = -99
! Subprogram not used 
! Subprogram not used     bnds_loop: do i = time_ndx,n_profiles
! Subprogram not used 
! Subprogram not used        call read_buffered_datetime( datetime, i )
! Subprogram not used 
! Subprogram not used        if ( datetime>begdatetime .and. beg_ndx<0 ) beg_ndx = i
! Subprogram not used        if ( datetime>enddatetime ) exit bnds_loop
! Subprogram not used        end_ndx = i
! Subprogram not used 
! Subprogram not used     enddo bnds_loop
! Subprogram not used 
! Subprogram not used     if (beg_ndx == -99 .and. end_ndx== -99) then
! Subprogram not used        if (masterproc) write(iulog,'(a)')  'sat_hist->read_next_position: must be beyond last position -- returning.'
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! Advance the search forward, but because of ntimesteps, it is possible
! Subprogram not used     ! for observations used here to be used again. However, we should not go
! Subprogram not used     ! back before the previous beginning time.
! Subprogram not used     if (beg_ndx>0) time_ndx = beg_ndx
! Subprogram not used 
! Subprogram not used     ncols = end_ndx-beg_ndx+1
! Subprogram not used 
! Subprogram not used     if (ncols > 0) then
! Subprogram not used        allocate( obs_lats(ncols), obs_lons(ncols) )
! Subprogram not used        in_start_col = beg_ndx
! Subprogram not used 
! Subprogram not used        ierr = pio_get_var( infile, lat_vid, (/beg_ndx/), (/ncols/), obs_lats )
! Subprogram not used        ierr = pio_get_var( infile, lon_vid, (/beg_ndx/), (/ncols/), obs_lons )
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call t_stopf ('sat_hist::read_next_position')
! Subprogram not used   end subroutine read_next_position
  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Subprogram not used   subroutine write_record_coord( tape, mod_lats, mod_lons, mod_dists, ncols, nfils )
! Subprogram not used 
! Subprogram not used     use time_manager,  only: get_nstep, get_curr_date, get_curr_time
! Subprogram not used     use cam_history_support, only : active_entry
! Subprogram not used     implicit none
! Subprogram not used     type(active_entry), intent(inout) :: tape
! Subprogram not used 
! Subprogram not used     integer,  intent(in) :: ncols
! Subprogram not used     real(r8), intent(in) :: mod_lats(ncols * sathist_nclosest)
! Subprogram not used     real(r8), intent(in) :: mod_lons(ncols * sathist_nclosest)
! Subprogram not used     real(r8), intent(in) :: mod_dists(ncols * sathist_nclosest)
! Subprogram not used     integer,  intent(in) ::  nfils
! Subprogram not used 
! Subprogram not used     integer :: t, ierr, i
! Subprogram not used     integer :: yr, mon, day      ! year, month, and day components of a date
! Subprogram not used     integer :: nstep             ! current timestep number
! Subprogram not used     integer :: ncdate            ! current date in integer format [yyyymmdd]
! Subprogram not used     integer :: ncsec             ! current time of day [seconds]
! Subprogram not used     integer :: ndcur             ! day component of current time
! Subprogram not used     integer :: nscur             ! seconds component of current time
! Subprogram not used     real(r8) :: time             ! current time
! Subprogram not used     integer, allocatable  :: itmp(:)
! Subprogram not used     real(r8), allocatable :: rtmp(:)
! Subprogram not used     real(r8), allocatable :: out_lats(:)
! Subprogram not used     real(r8), allocatable :: out_lons(:)
! Subprogram not used 
! Subprogram not used     call t_startf ('sat_hist::write_record_coord')
! Subprogram not used 
! Subprogram not used     nstep = get_nstep()
! Subprogram not used     call get_curr_date(yr, mon, day, ncsec)
! Subprogram not used     ncdate = yr*10000 + mon*100 + day
! Subprogram not used     call get_curr_time(ndcur, nscur)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     time = ndcur + nscur/86400._r8
! Subprogram not used 
! Subprogram not used     allocate( itmp(ncols * sathist_nclosest) )
! Subprogram not used     allocate( rtmp(ncols * sathist_nclosest) )
! Subprogram not used     
! Subprogram not used     itmp(:) = ncdate
! Subprogram not used     ierr = pio_put_var(tape%File, tape%dateid,(/nfils/), (/ncols * sathist_nclosest/),itmp)
! Subprogram not used     itmp(:) = ncsec
! Subprogram not used     ierr = pio_put_var(tape%File, tape%datesecid,(/nfils/),(/ncols * sathist_nclosest/),itmp)
! Subprogram not used     rtmp(:) = time
! Subprogram not used     ierr = pio_put_var(tape%File, tape%timeid, (/nfils/),(/ncols * sathist_nclosest/),rtmp)
! Subprogram not used     
! Subprogram not used     deallocate(itmp)
! Subprogram not used     deallocate(rtmp)
! Subprogram not used 
! Subprogram not used     ! output model column coordinates
! Subprogram not used     ierr = pio_put_var(tape%File, out_latid, (/nfils/),(/ncols * sathist_nclosest/), mod_lats)
! Subprogram not used     ierr = pio_put_var(tape%File, out_lonid, (/nfils/),(/ncols * sathist_nclosest/), mod_lons)
! Subprogram not used     ierr = pio_put_var(tape%File, out_dstid, (/nfils/),(/ncols * sathist_nclosest/), mod_dists / 1000._r8)
! Subprogram not used     
! Subprogram not used     ! output instrument location
! Subprogram not used     allocate( out_lats(ncols * sathist_nclosest) )
! Subprogram not used     allocate( out_lons(ncols * sathist_nclosest) )
! Subprogram not used     
! Subprogram not used     do i = 1, ncols
! Subprogram not used       out_lats(((i-1)*sathist_nclosest)+1 : (i*sathist_nclosest)) = obs_lats(i)
! Subprogram not used       out_lons(((i-1)*sathist_nclosest)+1 : (i*sathist_nclosest)) = obs_lons(i)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ierr = pio_put_var(tape%File, out_instr_lat_vid, (/nfils/),(/ncols * sathist_nclosest/), out_lats)
! Subprogram not used     ierr = pio_put_var(tape%File, out_instr_lon_vid, (/nfils/),(/ncols * sathist_nclosest/), out_lons)
! Subprogram not used 
! Subprogram not used     deallocate(out_lats)
! Subprogram not used     deallocate(out_lons)
! Subprogram not used     
! Subprogram not used     
! Subprogram not used     ierr = copy_data( infile, date_vid, tape%File, out_obs_date_vid, in_start_col, nfils, ncols )
! Subprogram not used     ierr = copy_data( infile, time_vid, tape%File, out_obs_time_vid, in_start_col, nfils, ncols )
! Subprogram not used     
! Subprogram not used     ! output observation identifiers
! Subprogram not used     if (instr_vid>0) then
! Subprogram not used        ierr = copy_data( infile, instr_vid, tape%File, out_instrid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (orbit_vid>0) then
! Subprogram not used        ierr = copy_data( infile, orbit_vid, tape%File, out_orbid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (prof_vid>0) then
! Subprogram not used        ierr = copy_data( infile, prof_vid, tape%File, out_profid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (zenith_vid>0) then
! Subprogram not used        ierr = copy_data( infile, zenith_vid, tape%File, out_zenithid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (in_julian_vid>0) then
! Subprogram not used        ierr = copy_data( infile, in_julian_vid, tape%File, out_julian_vid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (in_occ_type_vid>0) then
! Subprogram not used        ierr = copy_data( infile, in_occ_type_vid, tape%File, out_occ_type_vid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (in_localtime_vid>0) then
! Subprogram not used        ierr = copy_data( infile, in_localtime_vid, tape%File, out_localtime_vid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used     if (in_doy_vid>0) then
! Subprogram not used        ierr = copy_data( infile, in_doy_vid, tape%File, out_doy_vid, in_start_col, nfils, ncols )
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call t_stopf ('sat_hist::write_record_coord')
! Subprogram not used   end subroutine write_record_coord

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! Subprogram not used   subroutine get_indices( lats, lons, ncols, nocols, has_dyn_flds, col_ndxs, chk_ndxs, &
! Subprogram not used        fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners, mlats, mlons, phs_dists )
! Subprogram not used 
! Subprogram not used     use dyn_grid, only : dyn_grid_get_colndx
! Subprogram not used     use phys_grid, only: get_rlat_p, get_rlon_p
! Subprogram not used 
! Subprogram not used     integer,  intent(in)  :: ncols
! Subprogram not used     real(r8), intent(in)  :: lats(ncols)
! Subprogram not used     real(r8), intent(in)  :: lons(ncols)
! Subprogram not used     integer,  intent(in)  :: nocols
! Subprogram not used     logical,  intent(in)  :: has_dyn_flds
! Subprogram not used     integer,  intent(out) :: col_ndxs(nocols)
! Subprogram not used     integer,  intent(out) :: chk_ndxs(nocols)
! Subprogram not used     integer,  intent(out) :: fdyn_ndxs(nocols)
! Subprogram not used     integer,  intent(out) :: ldyn_ndxs(nocols)
! Subprogram not used     integer,  intent(out) :: phs_owners(nocols)
! Subprogram not used     integer,  intent(out) :: dyn_owners(nocols)
! Subprogram not used     real(r8), intent(out) :: mlats(nocols)
! Subprogram not used     real(r8), intent(out) :: mlons(nocols)
! Subprogram not used     real(r8), intent(out) :: phs_dists(nocols)
! Subprogram not used 
! Subprogram not used     integer :: i, j, ndx
! Subprogram not used     real(r8) :: lat, lon
! Subprogram not used     
! Subprogram not used     integer,  allocatable :: ichks(:),icols(:),idyn1s(:),idyn2s(:), iphs_owners(:), idyn_owners(:)
! Subprogram not used     real(r8), allocatable :: rlats(:), rlons(:), plats(:), plons(:), iphs_dists(:)
! Subprogram not used 
! Subprogram not used     integer :: gcols(sathist_nclosest)
! Subprogram not used 
! Subprogram not used     call t_startf ('sat_hist::get_indices')
! Subprogram not used 
! Subprogram not used     allocate(ichks(sathist_nclosest),icols(sathist_nclosest),idyn1s(sathist_nclosest), &
! Subprogram not used          idyn2s(sathist_nclosest),iphs_owners(sathist_nclosest),idyn_owners(sathist_nclosest))
! Subprogram not used     allocate(rlats(sathist_nclosest), rlons(sathist_nclosest), plats(sathist_nclosest), &
! Subprogram not used          plons(sathist_nclosest), iphs_dists(sathist_nclosest) )
! Subprogram not used 
! Subprogram not used     col_ndxs = -1
! Subprogram not used     chk_ndxs = -1
! Subprogram not used     fdyn_ndxs = -1
! Subprogram not used     ldyn_ndxs = -1
! Subprogram not used     phs_owners = -1
! Subprogram not used     dyn_owners = -1
! Subprogram not used     phs_dists = -1
! Subprogram not used 
! Subprogram not used     ndx = 0
! Subprogram not used     do i = 1,ncols
! Subprogram not used 
! Subprogram not used        lat = lats(i)
! Subprogram not used        lon = lons(i)
! Subprogram not used 
! Subprogram not used        if ( lon >= 360._r8) then
! Subprogram not used          lon = lon-360._r8
! Subprogram not used        endif
! Subprogram not used        if ( lon < 0._r8) then
! Subprogram not used          lon = lon+360._r8
! Subprogram not used        endif
! Subprogram not used        if (lat<-90._r8 .or. lat>90._r8) then
! Subprogram not used           write(iulog,*) 'sat_hist::get_indices lat = ',lat
! Subprogram not used           call endrun('sat_hist::get_indices : lat must be between -90 and 90 degrees (-90<=lat<=90)')
! Subprogram not used        endif
! Subprogram not used        if (lon<0._r8 .or. lon>=360._r8) then
! Subprogram not used           write(iulog,*) 'sat_hist::get_indices lon = ',lon
! Subprogram not used           call endrun('sat_hist::get_indices : lon must be between 0 and 360 degrees (0<=lon<360)')
! Subprogram not used        endif
! Subprogram not used        
! Subprogram not used        call find_cols( lat, lon, sathist_nclosest, iphs_owners, ichks, icols, &
! Subprogram not used                        gcols, iphs_dists, plats, plons )
! Subprogram not used 
! Subprogram not used        if (has_dyn_flds) then
! Subprogram not used           call dyn_grid_get_colndx( gcols, sathist_nclosest, idyn_owners, idyn1s, idyn2s )
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        do j = 1, sathist_nclosest
! Subprogram not used           
! Subprogram not used           if (debug .and. iam==iphs_owners(j) ) then
! Subprogram not used              if ( abs(plats(j)-rlats(j))>1.e-3_r8 ) then
! Subprogram not used                 write(*,'(a,3f20.12)') ' lat, plat, rlat = ', lat, plats(j), rlats(j)
! Subprogram not used                 write(*,'(a,3f20.12)') ' lon, plon, rlon = ', lon, plons(j), rlons(j)
! Subprogram not used                 call endrun('sat_hist::get_indices: dyn lat is different than phys lat ')
! Subprogram not used              endif
! Subprogram not used              if ( abs(plons(j)-rlons(j))>1.e-3_r8 ) then
! Subprogram not used                 write(*,'(a,3f20.12)') ' lat, plat, rlat = ', lat, plats(j), rlats(j)
! Subprogram not used                 write(*,'(a,3f20.12)') ' lon, plon, rlon = ', lon, plons(j), rlons(j)
! Subprogram not used                 call endrun('sat_hist::get_indices: dyn lon is different than phys lon ')
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used           
! Subprogram not used           ndx = ndx+1
! Subprogram not used           
! Subprogram not used           chk_ndxs(ndx)   = ichks(j)
! Subprogram not used           col_ndxs(ndx)   = icols(j)
! Subprogram not used           fdyn_ndxs(ndx)  = idyn1s(j)
! Subprogram not used           ldyn_ndxs(ndx)  = idyn2s(j)
! Subprogram not used           mlats(ndx)      = plats(j)
! Subprogram not used           mlons(ndx)      = plons(j)
! Subprogram not used           phs_owners(ndx) = iphs_owners(j)
! Subprogram not used           dyn_owners(ndx) = idyn_owners(j)
! Subprogram not used           phs_dists(ndx)  = iphs_dists(j)
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     deallocate(ichks, icols, idyn1s, idyn2s, iphs_owners, idyn_owners)
! Subprogram not used     deallocate(rlats, rlons, plats, plons, iphs_dists )
! Subprogram not used 
! Subprogram not used     call t_stopf ('sat_hist::get_indices')
! Subprogram not used   end subroutine get_indices

!-------------------------------------------------------------------------------
! utility function
!-------------------------------------------------------------------------------
! Subprogram not used   integer function define_var( var_name, coldim, infile, in_vid, outfile, out_id ) result(res)
! Subprogram not used 
! Subprogram not used     use pio, only: pio_inq_vartype
! Subprogram not used 
! Subprogram not used     character(len=*), intent(in) :: var_name
! Subprogram not used     integer,          intent(in) :: coldim
! Subprogram not used     type(File_desc_t),intent(inout) :: infile
! Subprogram not used     type(File_desc_t),intent(inout) :: outfile
! Subprogram not used     integer,          intent(in) :: in_vid
! Subprogram not used     type(var_desc_t), intent(out):: out_id
! Subprogram not used 
! Subprogram not used     integer :: type
! Subprogram not used 
! Subprogram not used     res = pio_inq_varid( outfile, var_name, out_id )
! Subprogram not used     if(res/=PIO_NOERR) then
! Subprogram not used 
! Subprogram not used        res = pio_inq_vartype( infile, in_vid, type )
! Subprogram not used 
! Subprogram not used        res = pio_def_var ( outfile, var_name, type, (/coldim/), out_id )
! Subprogram not used 
! Subprogram not used        res = copy_att( infile, in_vid, 'long_name', outfile, out_id )
! Subprogram not used        res = copy_att( infile, in_vid, 'units',     outfile, out_id )
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end function define_var

!-------------------------------------------------------------------------------
! utility function
!-------------------------------------------------------------------------------
! Subprogram not used   integer function copy_data( infile, in_vid, outfile, out_id, instart, outstart, ncols ) result(res)
! Subprogram not used 
! Subprogram not used     type(File_desc_t),intent(in) :: infile
! Subprogram not used     type(File_desc_t),intent(inout) :: outfile
! Subprogram not used     integer,          intent(in) :: in_vid
! Subprogram not used     type(var_desc_t), intent(in) :: out_id
! Subprogram not used     integer,          intent(in) :: instart, outstart, ncols
! Subprogram not used 
! Subprogram not used     real(r8), allocatable :: data(:)
! Subprogram not used     real(r8), allocatable :: outdata(:)
! Subprogram not used     integer               :: i
! Subprogram not used 
! Subprogram not used     allocate( data(ncols) )
! Subprogram not used 
! Subprogram not used     res = pio_get_var( infile,  in_vid, (/instart/),  (/ncols/), data )
! Subprogram not used 
! Subprogram not used     allocate( outdata(ncols * sathist_nclosest) )
! Subprogram not used     
! Subprogram not used     do i = 1, ncols
! Subprogram not used       outdata(((i-1)*sathist_nclosest)+1 : (i*sathist_nclosest)) = data(i)
! Subprogram not used     enddo
! Subprogram not used   
! Subprogram not used     res = pio_put_var( outfile, out_id, (/outstart/), (/ncols * sathist_nclosest/), outdata )
! Subprogram not used 
! Subprogram not used     deallocate(outdata)
! Subprogram not used     deallocate(data)
! Subprogram not used 
! Subprogram not used   end function copy_data

!-------------------------------------------------------------------------------
! utility function
! -- should be able to use pio_copy_att which does not seem to work
!-------------------------------------------------------------------------------
! Subprogram not used   integer function copy_att( infile, in_vid, att_name, outfile, out_id ) result(res)
! Subprogram not used 
! Subprogram not used     type(File_desc_t),intent(inout) :: infile
! Subprogram not used     type(File_desc_t),intent(inout) :: outfile
! Subprogram not used     character(len=*), intent(in) :: att_name
! Subprogram not used     integer,          intent(in) :: in_vid
! Subprogram not used     type(var_desc_t), intent(in) :: out_id
! Subprogram not used 
! Subprogram not used     character(len=1024) :: att
! Subprogram not used     
! Subprogram not used 
! Subprogram not used     res = pio_get_att( infile, in_vid, trim(att_name), att )
! Subprogram not used     if (res==PIO_NOERR) then
! Subprogram not used        res = pio_put_att ( outfile, out_id, trim(att_name), trim(att))
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end function copy_att
  
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
! Subprogram not used   subroutine find_cols(lat, lon, nclosest, owner, lcid, icol, gcol, distmin, mlats, mlons)
! Subprogram not used     use physconst,  only: rearth
! Subprogram not used     use phys_grid,  only: get_rlon_all_p, get_rlat_all_p, get_gcol_p, get_ncols_p
! Subprogram not used     use spmd_utils, only: iam, npes, mpi_integer, mpi_real8, mpicom
! Subprogram not used 
! Subprogram not used     real(r8),intent(in)  :: lat, lon            ! requested location in degrees
! Subprogram not used     integer, intent(in)  :: nclosest            ! number of closest points to find
! Subprogram not used     integer, intent(out) :: owner(nclosest)     ! rank of chunk owner
! Subprogram not used     integer, intent(out) :: lcid(nclosest)      ! local chunk index
! Subprogram not used     integer, intent(out) :: icol(nclosest)      ! column index within the chunk
! Subprogram not used     integer, intent(out) :: gcol(nclosest)      ! global column index 
! Subprogram not used     real(r8),intent(out) :: distmin(nclosest)   ! the distance (m) of the closest column(s)
! Subprogram not used     real(r8),intent(out) :: mlats(nclosest)     ! the latitude of the closest column(s)
! Subprogram not used     real(r8),intent(out) :: mlons(nclosest)     ! the longitude of the closest column(s)
! Subprogram not used 
! Subprogram not used     real(r8) :: dist
! Subprogram not used     real(r8) :: rlats(pcols), rlons(pcols)
! Subprogram not used     real(r8) :: latr, lonr
! Subprogram not used 
! Subprogram not used     integer :: my_owner(nclosest)
! Subprogram not used     integer :: my_lcid(nclosest)
! Subprogram not used     integer :: my_icol(nclosest)
! Subprogram not used     integer :: my_gcol(nclosest)
! Subprogram not used     real(r8) :: my_distmin(nclosest)
! Subprogram not used     real(r8) :: my_mlats(nclosest)
! Subprogram not used     real(r8) :: my_mlons(nclosest)
! Subprogram not used 
! Subprogram not used     integer  :: c, i, j, k, ierr, ncols, mindx(1)
! Subprogram not used     real(r8) :: sendbufr(3)
! Subprogram not used     real(r8) :: recvbufr(3,npes)
! Subprogram not used     integer  :: sendbufi(4)
! Subprogram not used     integer  :: recvbufi(4,npes)
! Subprogram not used 
! Subprogram not used     call t_startf ('sat_hist::find_cols')
! Subprogram not used 
! Subprogram not used     latr = lat/rad2deg              ! to radians
! Subprogram not used     lonr = lon/rad2deg              ! to radians
! Subprogram not used     
! Subprogram not used     my_owner(:)   = -999
! Subprogram not used     my_lcid(:)    = -999
! Subprogram not used     my_icol(:)    = -999
! Subprogram not used     my_gcol(:)    = -999
! Subprogram not used     my_mlats(:)   = -999
! Subprogram not used     my_mlons(:)   = -999
! Subprogram not used     my_distmin(:) = 1.e10_r8
! Subprogram not used 
! Subprogram not used     chk_loop: do c=begchunk,endchunk
! Subprogram not used        ncols = get_ncols_p(c)
! Subprogram not used        call get_rlat_all_p(c, pcols, rlats)
! Subprogram not used        call get_rlon_all_p(c, pcols, rlons)
! Subprogram not used 
! Subprogram not used        col_loop: do i = 1,ncols
! Subprogram not used           ! Use the Spherical Law of Cosines to find the great-circle distance.
! Subprogram not used           dist = acos(sin(latr) * sin(rlats(i)) + cos(latr) * cos(rlats(i)) * cos(rlons(i) - lonr)) * rearth       
! Subprogram not used 
! Subprogram not used           closest_loop: do j = nclosest, 1, -1
! Subprogram not used              if (dist < my_distmin(j)) then
! Subprogram not used 
! Subprogram not used                 if (j < nclosest) then
! Subprogram not used                    my_distmin(j+1) = my_distmin(j)
! Subprogram not used                    my_owner(j+1)   = my_owner(j)
! Subprogram not used                    my_lcid(j+1)    = my_lcid(j)
! Subprogram not used                    my_icol(j+1)    = my_icol(j)
! Subprogram not used                    my_gcol(j+1)    = my_gcol(j)
! Subprogram not used                    my_mlats(j+1)   = my_mlats(j)
! Subprogram not used                    my_mlons(j+1)   = my_mlons(j)
! Subprogram not used                 end if
! Subprogram not used 
! Subprogram not used                 my_distmin(j) = dist
! Subprogram not used                 my_owner(j)   = iam
! Subprogram not used                 my_lcid(j)    = c
! Subprogram not used                 my_icol(j)    = i
! Subprogram not used                 my_gcol(j)    = get_gcol_p(c,i)
! Subprogram not used                 my_mlats(j)   = rlats(i) * rad2deg
! Subprogram not used                 my_mlons(j)   = rlons(i) * rad2deg
! Subprogram not used              else
! Subprogram not used                 exit 
! Subprogram not used              end if
! Subprogram not used           enddo closest_loop
! Subprogram not used 
! Subprogram not used        enddo col_loop
! Subprogram not used     enddo chk_loop
! Subprogram not used 
! Subprogram not used     k = 1
! Subprogram not used 
! Subprogram not used     do j = 1, nclosest
! Subprogram not used 
! Subprogram not used        sendbufr(1) = my_distmin(k)
! Subprogram not used        sendbufr(2) = my_mlats(k)
! Subprogram not used        sendbufr(3) = my_mlons(k)
! Subprogram not used 
! Subprogram not used        call mpi_allgather( sendbufr, 3, mpi_real8, recvbufr, 3, mpi_real8, mpicom, ierr )
! Subprogram not used 
! Subprogram not used        mindx = minloc(recvbufr(1,:))
! Subprogram not used        distmin(j) = recvbufr(1,mindx(1))
! Subprogram not used        mlats(j)   = recvbufr(2,mindx(1))
! Subprogram not used        mlons(j)   = recvbufr(3,mindx(1))
! Subprogram not used 
! Subprogram not used        sendbufi(1) = my_owner(k)
! Subprogram not used        sendbufi(2) = my_lcid(k)
! Subprogram not used        sendbufi(3) = my_icol(k)
! Subprogram not used        sendbufi(4) = my_gcol(k)
! Subprogram not used 
! Subprogram not used        call mpi_allgather( sendbufi, 4, mpi_integer, recvbufi, 4, mpi_integer, mpicom, ierr )
! Subprogram not used 
! Subprogram not used        owner(j)   = recvbufi(1,mindx(1))
! Subprogram not used        lcid(j)    = recvbufi(2,mindx(1))
! Subprogram not used        icol(j)    = recvbufi(3,mindx(1))
! Subprogram not used        gcol(j)    = recvbufi(4,mindx(1))
! Subprogram not used 
! Subprogram not used        if ( iam == owner(j) ) then
! Subprogram not used           k = k+1
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call t_stopf ('sat_hist::find_cols')
! Subprogram not used 
! Subprogram not used   end subroutine find_cols

end module sat_hist
