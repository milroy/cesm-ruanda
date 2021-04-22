module spedata
!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: spedata
!
! !DESCRIPTION
! Handles reading and interpolating solar proton ionization data.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8,r4 => shr_kind_r4, shr_kind_cl
  use shr_cal_mod,  only: shr_cal_gregorian
  use time_manager, only: get_curr_date, get_curr_calday, get_step_size, timemgr_is_caltype
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver, begchunk, endchunk
  use abortutils,   only: endrun
  use pio,          only: var_desc_t, file_desc_t, pio_get_var, pio_get_att, &
                          pio_setdebuglevel, pio_seterrorhandling, pio_bcast_error, &
                          pio_internal_error, pio_noerr, pio_inq_varid, pio_put_var, &
                          pio_put_att, pio_inq_dimid, pio_char, pio_def_dim, pio_def_var, &
                          pio_inq_dimlen, pio_closefile
  use cam_pio_utils,only: cam_pio_openfile
  use perf_mod,     only: t_startf, t_stopf
  use cam_logfile,  only: iulog


  use mpishorthand, only: mpicom, mpir8, mpiint, mpichar


  implicit none

  private  ! all unless made public
  save 

! !PUBLIC MEMBERS

  public spedata_init            ! subroutine to open files, allocate blocked arrays, etc
  public advance_spedata         ! subroutine to read more data and interpolate
  public get_ionpairs_profile    ! interface to get ionization profile
  public spedata_defaultopts
  public spedata_setopts
  public spe_run

  !------------------------------------------------------------------
  ! Interface to access the meteorology fields.  Possible invocations
  ! are as follows:
  !   call get_met_fields( physics_state, us, vs , tend, dt )
  !   call get_met_fields( u, v )
  !   call get_met_fields( cam_in )
  !------------------------------------------------------------------
  Interface get_ionpairs_profile                        ! overload accessors
     Module Procedure get_ionpairs
  End Interface
  
  logical :: remove_spe_file = .false.  ! delete data file when finished with it

! !REVISION HISTORY:
!   25 Oct 2005  Francis Vitt     Creation
!
! EOP
!----------------------------------------------------------------------- 
! $Id: spedata.F90,v 1.1.2.1 2006/05/03 20:53:09 stacy Exp $
! $Author: stacy $
!----------------------------------------------------------------------- 

  type input_profile
     real(r8), dimension(:), pointer :: data => null()
  endtype input_profile


  real(r8), allocatable :: ionpairs(:)  ! interpolated ionization profile
  real(r8), allocatable :: prod_pressures(:) 

  type(input_profile) :: ionpairs_i(2)

  integer :: prod_id ! var id of the data in the netCDF
  integer :: nprod_press

  integer :: dateid           ! var id of the date in the netCDF
  integer :: secid            ! var id of the sec data 
  real(r8) :: datatimem = -1.e36_r8     ! time of prv. values read in
  real(r8) :: datatimep = -1.e36_r8     ! time of nxt. values read in

  integer, parameter :: nm=1    ! array index for previous (minus) data
  integer, parameter :: np=2    ! array indes for next (plus) data

  real(r8) :: curr_mod_time ! model time - calendar day
  real(r8) :: next_mod_time ! model time - calendar day - next time step

  character(len=shr_kind_cl) :: curr_filename = ' '
  character(len=shr_kind_cl) :: next_filename = ' '
  character(len=shr_kind_cl) :: spedata_file = ' '
  type(file_desc_t) :: curr_fileid, next_fileid     ! the pio id of the NetCDF file
  ! IDs of the filename strings in the restart file.
  type(var_desc_t), pointer :: currfnameid => null()
  type(var_desc_t), pointer :: nextfnameid => null()
  real(r8), pointer, dimension(:) :: curr_data_times => null()
  real(r8), pointer, dimension(:) :: next_data_times => null()
  character(len=shr_kind_cl) :: filenames_list = ''

  character(len=16) :: calendar

  logical, protected :: spe_run = .false.

contains

!--------------------------------------------------------------------------
! Get the default runtime options
!--------------------------------------------------------------------------
  subroutine spedata_defaultopts( spe_data_file_out, &
                                  spe_remove_file_out , &
                                  spe_filenames_list_out )

    implicit none

    character(len=shr_kind_cl), intent(out), optional :: spe_data_file_out
    character(len=shr_kind_cl), intent(out), optional :: spe_filenames_list_out
    logical, intent(out), optional :: spe_remove_file_out

    if ( present( spe_data_file_out ) ) then
       spe_data_file_out = spedata_file
    endif

    if ( present( spe_remove_file_out ) ) then
       spe_remove_file_out = remove_spe_file
    endif

    if ( present( spe_filenames_list_out ) ) then
       spe_filenames_list_out =  filenames_list
    endif

  end subroutine spedata_defaultopts

!--------------------------------------------------------------------------
! Set runtime options
!--------------------------------------------------------------------------
  subroutine spedata_setopts( spe_data_file_in, &
                              spe_remove_file_in, &
                              spe_filenames_list_in ) 

    implicit none

    character(len=shr_kind_cl), intent(in), optional :: spe_data_file_in
    character(len=shr_kind_cl), intent(in), optional :: spe_filenames_list_in
    logical, intent(in), optional :: spe_remove_file_in

    integer :: ierr

    if ( present( spe_data_file_in ) ) then
       spedata_file = spe_data_file_in 
    endif

    if ( present( spe_remove_file_in ) ) then
       remove_spe_file = spe_remove_file_in
    endif

    if ( present( spe_filenames_list_in ) ) then
       filenames_list = spe_filenames_list_in 
    endif

    if (len_trim(spedata_file) > 0) then 
      spe_run = .true.
    else
      spe_run = .false.
    endif
      
    if (masterproc) then
       write(iulog,*) 'This run includes NOx and HOx production from solor proton events: ',spe_run
       if ( spe_run ) then
         write(iulog,*) 'Time-variant solar proton ionization dataset (spedata_file) is: ', trim(spedata_file)
         write(iulog,*) 'Solar proton ionization data file will be removed (remove_spe_file): ', remove_spe_file
         write(iulog,*) 'Solar proton ionization data file names list file: ', trim(filenames_list) 
       endif
    endif

  end subroutine spedata_setopts

  subroutine spedata_init()
!--------------------------------------------------------------------------
! Opens file, allocates arrays
!--------------------------------------------------------------------------

    use cam_control_mod, only : nsrest
    use string_utils,    only : to_upper
    use mo_apex,         only : apexmag

    implicit none

!--------------------------------------------------------------------------
! 	... local variables
!--------------------------------------------------------------------------
    integer :: astat
    integer :: ierr

    if (.not. spe_run) return

    ! initialize geo-magnetic coordinate module ...
    call apexmag()

!--------------------------------------------------------------------------
! 	... initial run ?
!--------------------------------------------------------------------------
!    if (nsrest == 0) then
       curr_filename = spedata_file
       next_filename = ''
!    endif

    call open_spe_datafile( curr_filename, curr_fileid, curr_data_times, read_pressures=.true. )

    if ( len_trim(next_filename) > 0 ) &
         call open_spe_datafile( next_filename, next_fileid, next_data_times )
    
    ierr = pio_inq_varid( curr_fileid, 'Prod', prod_id )
    

!--------------------------------------------------------------------------
! allocate space for data arrays ...
!--------------------------------------------------------------------------
    allocate( ionpairs(nprod_press),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'spedata_init: failed to allocate ionpairs; error = ',astat
       call endrun
    end if
    allocate( ionpairs_i(nm)%data(nprod_press),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'spedata_init: failed to allocate ionpairs_i; error = ',astat
       call endrun
    end if
    allocate( ionpairs_i(np)%data(nprod_press),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'spedata_init: failed to allocate ionpairs_i; error = ',astat
       call endrun
    end if

 end subroutine spedata_init


!-----------------------------------------------------------------------
! Reads more data if needed and interpolates data to current model time 
!-----------------------------------------------------------------------
 subroutine advance_spedata()

    implicit none

    if (.not. spe_run) return

    call t_startf('MET__advance')

    call get_model_time()

    if ( ( curr_mod_time > datatimep ) ) then
       call check_files()
    endif

    if ( curr_mod_time > datatimep ) then
       call read_next_spedata()
    end if

! need to inperpolate the data, regardless !
! each mpi tasks needs to interpolate
    call interpolate_spedata()

    call t_stopf('MET__advance')

  end subroutine advance_spedata

!------------------------------------------------------------------------
! get the ion pair production profile
!------------------------------------------------------------------------
! Subprogram not used   subroutine get_ionpairs( model_pressures, pairs )
! Subprogram not used     
! Subprogram not used     use interpolate_data, only : lininterp
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     real(r8), intent(in)  :: model_pressures(:)
! Subprogram not used     real(r8), intent(out) :: pairs(:)
! Subprogram not used 
! Subprogram not used     integer :: npress, ub,lb
! Subprogram not used 
! Subprogram not used     pairs(:) = 0._r8
! Subprogram not used 
! Subprogram not used     if (.not. spe_run) return
! Subprogram not used 
! Subprogram not used     npress = size(model_pressures)
! Subprogram not used 
! Subprogram not used     ub = ubound( prod_pressures,1 )
! Subprogram not used     lb = lbound( prod_pressures,1 )
! Subprogram not used 
! Subprogram not used     ! interpolate to model levels
! Subprogram not used     call lininterp( ionpairs(ub:lb:-1), prod_pressures(ub:lb:-1)*1.e2_r8, nprod_press, pairs, model_pressures, npress )
! Subprogram not used 
! Subprogram not used   end subroutine get_ionpairs

!------------------------------------------------------------------------------
! internal methods :
!------------------------------------------------------------------------------

! Subprogram not used   subroutine get_model_time()
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer yr, mon, day, ncsec  ! components of a date
! Subprogram not used 
! Subprogram not used     call t_startf('MET__get_model_time')
! Subprogram not used 
! Subprogram not used     call get_curr_date(yr, mon, day, ncsec)
! Subprogram not used 
! Subprogram not used     curr_mod_time = get_time_float( yr, mon, day, ncsec )
! Subprogram not used     next_mod_time = curr_mod_time + get_step_size()/86400._r8
! Subprogram not used 
! Subprogram not used     call t_stopf('MET__get_model_time')
! Subprogram not used 
! Subprogram not used   end subroutine get_model_time

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subprogram not used   subroutine check_files()
! Subprogram not used 
! Subprogram not used     use shr_sys_mod, only: shr_sys_system
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used     character(len=128) :: ctmp
! Subprogram not used     integer ::  istat
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (next_mod_time > curr_data_times(size(curr_data_times))) then
! Subprogram not used        if ( .not. associated(next_data_times) ) then
! Subprogram not used           ! open next file...
! Subprogram not used           next_filename = incr_filename( curr_filename )
! Subprogram not used           call open_spe_datafile( next_filename, next_fileid, next_data_times )
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if ( associated(next_data_times) ) then
! Subprogram not used        if (curr_mod_time >= next_data_times(1)) then
! Subprogram not used 
! Subprogram not used           ! close current file ...
! Subprogram not used           call pio_closefile( curr_fileid )
! Subprogram not used           ! remove if requested
! Subprogram not used           if(masterproc) then
! Subprogram not used              if( remove_spe_file ) then
! Subprogram not used                 write(iulog,*) 'check_files: removing file = ',trim(curr_filename) 
! Subprogram not used                 ctmp = 'rm -f ' // trim(curr_filename) 
! Subprogram not used                 write(iulog,*) 'check_files: fsystem issuing command - '
! Subprogram not used                 write(iulog,*) trim(ctmp)
! Subprogram not used                 call shr_sys_system( ctmp, istat )
! Subprogram not used              end if
! Subprogram not used           endif
! Subprogram not used 
! Subprogram not used           curr_filename = next_filename
! Subprogram not used           curr_fileid = next_fileid
! Subprogram not used 
! Subprogram not used           deallocate( curr_data_times )
! Subprogram not used           allocate( curr_data_times( size( next_data_times ) ) )
! Subprogram not used           curr_data_times(:) = next_data_times(:)
! Subprogram not used 
! Subprogram not used           next_filename = ''
! Subprogram not used 
! Subprogram not used           deallocate( next_data_times )
! Subprogram not used           nullify(  next_data_times )
! Subprogram not used 
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine check_files

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subprogram not used   function incr_filename( filename )
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... Increment or decrement a date string withing a filename
! Subprogram not used     !           the filename date section is assumed to be of the form
! Subprogram not used     !           yyyy-dd-mm
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use string_utils,  only : incstr
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     character(len=*), intent(in) :: filename ! present dynamical dataset filename
! Subprogram not used     character(len=shr_kind_cl) :: incr_filename      ! next filename in the sequence
! Subprogram not used 
! Subprogram not used     ! set new next_filename ...
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !	... local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer :: pos, pos1, istat
! Subprogram not used     character(len=shr_kind_cl) :: fn_new, line
! Subprogram not used     character(len=6)   :: seconds
! Subprogram not used     character(len=5)   :: num
! Subprogram not used 
! Subprogram not used     integer :: ios
! Subprogram not used 
! Subprogram not used     if ( len_trim(filenames_list) .eq. 0) then
! Subprogram not used 
! Subprogram not used        ! increment the number in the filename to the next number
! Subprogram not used        pos = len_trim( filename )
! Subprogram not used        fn_new = filename(:pos)
! Subprogram not used        write(iulog,*) 'incr_flnm: old filename = ',trim(fn_new)
! Subprogram not used        if( fn_new(pos-2:) == '.nc' ) then
! Subprogram not used           pos = pos - 3
! Subprogram not used        end if
! Subprogram not used        istat = incstr( fn_new(:pos), 1 )
! Subprogram not used        if( istat /= 0 ) then
! Subprogram not used           write(iulog,*) 'incr_flnm: incstr returned ', istat
! Subprogram not used           write(iulog,*) '           while trying to decrement ',trim( fn_new )
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        ! open filenames_list
! Subprogram not used        write(iulog,*) 'incr_flnm: old filename = ',trim(filename)
! Subprogram not used        write(iulog,*) 'incr_flnm: open filenames_list : ',filenames_list 
! Subprogram not used        open( unit=9, file=filenames_list, iostat=ios, status="OLD")
! Subprogram not used        if (ios /= 0) then
! Subprogram not used           call endrun('not able to open filenames_list file: '//filenames_list)
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        ! read file names
! Subprogram not used        read( unit=9, fmt='(A)', iostat=ios ) line 
! Subprogram not used        if (ios /= 0) then
! Subprogram not used           call endrun('not able to increment file name from filenames_list file: '//filenames_list)
! Subprogram not used        endif
! Subprogram not used        do while( trim(line) /= trim(filename) )
! Subprogram not used           read( unit=9, fmt='(A)', iostat=ios ) line 
! Subprogram not used           if (ios /= 0) then
! Subprogram not used              call endrun('not able to increment file name from filenames_list file: '//filenames_list)
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        read( unit=9, fmt='(A)', iostat=ios ) line 
! Subprogram not used        if (ios /= 0) then
! Subprogram not used           call endrun('not able to increment file name from filenames_list file: '//filenames_list)
! Subprogram not used        endif
! Subprogram not used        fn_new = trim(line)
! Subprogram not used 
! Subprogram not used        close(unit=9)
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     incr_filename = trim(fn_new)
! Subprogram not used     write(iulog,*) 'incr_flnm: new filename = ',incr_filename
! Subprogram not used 
! Subprogram not used   end function incr_filename

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subprogram not used   subroutine find_times( itms, fids, datatm, datatp, time )
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !------------------------------------------------------------------------------
! Subprogram not used     integer, intent(out) :: itms(2) ! record numbers that bracket time
! Subprogram not used     type(file_desc_t), intent(out) :: fids(2) ! ids of files that contains these recs
! Subprogram not used     real(r8), intent(in) :: time    ! time of interest
! Subprogram not used     real(r8), intent(out):: datatm, datatp
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !------------------------------------------------------------------------------
! Subprogram not used     integer np1        ! current forward time index of dataset
! Subprogram not used     integer n,i      ! 
! Subprogram not used     integer :: curr_tsize, next_tsize, all_tsize
! Subprogram not used     real(r8), allocatable, dimension(:):: all_data_times
! Subprogram not used     logical :: found_time
! Subprogram not used 
! Subprogram not used     curr_tsize = size(curr_data_times)
! Subprogram not used     next_tsize = 0
! Subprogram not used     if ( associated(next_data_times)) next_tsize = size(next_data_times)
! Subprogram not used 
! Subprogram not used     all_tsize = curr_tsize + next_tsize
! Subprogram not used 
! Subprogram not used     allocate( all_data_times( all_tsize ) )
! Subprogram not used 
! Subprogram not used     all_data_times(:curr_tsize) = curr_data_times(:)
! Subprogram not used     if (next_tsize > 0) all_data_times(curr_tsize+1:all_tsize) = next_data_times(:)
! Subprogram not used 
! Subprogram not used     ! find bracketing times 
! Subprogram not used     found_time = .false.
! Subprogram not used     do n = 1,all_tsize-1
! Subprogram not used        np1 = n + 1
! Subprogram not used        datatm = all_data_times(n)
! Subprogram not used        datatp = all_data_times(np1)
! Subprogram not used        if ( (time .ge. datatm) .and. (time .le. datatp) ) then
! Subprogram not used           found_time = .true.
! Subprogram not used           exit
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used   
! Subprogram not used     if( found_time ) then
! Subprogram not used        deallocate( all_data_times )
! Subprogram not used   
! Subprogram not used        itms(1) = n
! Subprogram not used        itms(2) = np1
! Subprogram not used        fids(:) = curr_fileid
! Subprogram not used   
! Subprogram not used        do i=1,2
! Subprogram not used           if ( itms(i) > curr_tsize ) then 
! Subprogram not used              itms(i) = itms(i) - curr_tsize 
! Subprogram not used              fids(i) = next_fileid
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'FIND_TIMES: Failed to find dates bracketing desired time =', time
! Subprogram not used        write(iulog,*)' datatm = ',datatm
! Subprogram not used        write(iulog,*)' datatp = ',datatp
! Subprogram not used        write(iulog,*)' all_data_times = ',all_data_times
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine find_times

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! Subprogram not used   subroutine read_next_spedata()
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used     integer :: i
! Subprogram not used     integer :: recnos(2)
! Subprogram not used     type(file_desc_t) :: fids(2)
! Subprogram not used     integer :: cnt(2)            ! array of counts for each dimension
! Subprogram not used     integer :: strt(2)           ! array of starting indices
! Subprogram not used     integer :: ierr
! Subprogram not used     call t_startf('MET__read_next_spedata')
! Subprogram not used 
! Subprogram not used     call find_times( recnos, fids, datatimem, datatimep, curr_mod_time )
! Subprogram not used 
! Subprogram not used     cnt(1)  = nprod_press
! Subprogram not used     cnt(2)  = 1
! Subprogram not used     strt(1) = 1
! Subprogram not used 
! Subprogram not used     do i = 1,2
! Subprogram not used        strt(2) = recnos(i)
! Subprogram not used        ierr = pio_get_var( fids(i), prod_id, strt, cnt,  ionpairs_i(i)%data )
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if (masterproc) write(iulog,*)'READ_NEXT_SPEDATA: Read soloar proton ionization data '
! Subprogram not used 
! Subprogram not used     call t_stopf('MET__read_next_spedata')
! Subprogram not used 
! Subprogram not used   end subroutine read_next_spedata


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subprogram not used   subroutine interpolate_spedata()
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     real(r4) fact1, fact2
! Subprogram not used     real(r8) deltat 
! Subprogram not used 
! Subprogram not used     call t_startf('MET__interpolate_spedata')
! Subprogram not used 
! Subprogram not used     deltat = datatimep - datatimem
! Subprogram not used 
! Subprogram not used     fact1 = (datatimep - curr_mod_time)/deltat
! Subprogram not used     fact2 = 1._r8-fact1
! Subprogram not used 
! Subprogram not used     ionpairs(:) = fact1*ionpairs_i(nm)%data(:) + fact2*ionpairs_i(np)%data(:)
! Subprogram not used 
! Subprogram not used     call t_stopf('MET__interpolate_spedata')
! Subprogram not used 
! Subprogram not used   end subroutine interpolate_spedata

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subprogram not used   subroutine get_dimension( fid, dname, dsize )
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(in) :: fid
! Subprogram not used     character(*), intent(in) :: dname
! Subprogram not used     integer, intent(out) :: dsize
! Subprogram not used 
! Subprogram not used     integer :: dimid, ierr
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_dimid( fid, dname, dimid )
! Subprogram not used     ierr = pio_inq_dimlen( fid, dimid, dsize )
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine get_dimension

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subprogram not used   subroutine open_spe_datafile( fname, fileid, times, read_pressures )
! Subprogram not used 
! Subprogram not used     use ioFileMod, only : getfil
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     character(*), intent(in) :: fname
! Subprogram not used     type(file_desc_t), intent(inout) :: fileid
! Subprogram not used     real(r8), pointer, intent(inout) :: times(:)
! Subprogram not used     logical, optional, intent(in) :: read_pressures
! Subprogram not used 
! Subprogram not used     character(len=shr_kind_cl) :: filen   
! Subprogram not used     integer :: year, month, day, dsize, i, timesize
! Subprogram not used     integer :: dateid,secid,press_id
! Subprogram not used     integer, allocatable , dimension(:) :: dates, datesecs
! Subprogram not used     integer :: ierr
! Subprogram not used     !
! Subprogram not used     ! open file and get fileid
! Subprogram not used     !
! Subprogram not used     call getfil( fname, filen, 0 )
! Subprogram not used     call cam_pio_openfile( fileid, filen, 0 )
! Subprogram not used     if(masterproc)   write(iulog,*)'open_met_datafile: ',trim(filen)
! Subprogram not used 
! Subprogram not used     call get_dimension( fileid, 'time', timesize )
! Subprogram not used 
! Subprogram not used     if ( associated(times) ) deallocate(times)
! Subprogram not used     allocate( times(timesize) )
! Subprogram not used 
! Subprogram not used     allocate( dates(timesize) )
! Subprogram not used     allocate( datesecs(timesize) )
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_varid( fileid, 'date',    dateid  )
! Subprogram not used     ierr = pio_inq_varid( fileid, 'datesec', secid  )
! Subprogram not used 
! Subprogram not used     ierr = pio_get_var( fileid, dateid, dates )
! Subprogram not used     ierr = pio_get_var( fileid, secid,  datesecs  )
! Subprogram not used 
! Subprogram not used     do i=1,timesize
! Subprogram not used        year = dates(i) / 10000
! Subprogram not used        month = mod(dates(i),10000)/100
! Subprogram not used        day = mod(dates(i),100)
! Subprogram not used        times(i) = get_time_float( year, month, day, datesecs(i) )
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     deallocate( dates )
! Subprogram not used     deallocate( datesecs )   
! Subprogram not used 
! Subprogram not used     if (present( read_pressures ) ) then
! Subprogram not used        call get_dimension( fileid, 'pressure', nprod_press )
! Subprogram not used 
! Subprogram not used        allocate( prod_pressures( nprod_press ) )
! Subprogram not used 
! Subprogram not used        ierr = pio_inq_varid( fileid, 'pressure',    press_id  )
! Subprogram not used        ierr = pio_get_var( fileid, press_id, (/ 1 /), (/ nprod_press /), prod_pressures )
! Subprogram not used     endif    
! Subprogram not used 
! Subprogram not used   end subroutine open_spe_datafile

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subprogram not used   function get_time_float( year, month, day, sec )
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! returns float representation of time -- number of days
! Subprogram not used ! since 1 jan 0001 00:00:00.000
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: year, month, day
! Subprogram not used     integer, intent(in) :: sec
! Subprogram not used     real(r8) :: get_time_float
! Subprogram not used 
! Subprogram not used ! ref date is 1 jan 0001
! Subprogram not used 
! Subprogram not used     integer  :: refyr, refmn, refdy
! Subprogram not used     real(r8) :: refsc, fltdy
! Subprogram not used     integer  :: doy(12)
! Subprogram not used 
! Subprogram not used !              jan feb mar apr may jun jul aug sep oct nov dec
! Subprogram not used !              31  28  31  30  31  30  31  31  31  31  30  31
! Subprogram not used     data doy /  1, 32, 60, 91,121,152,182,213,244,274,305,335 /
! Subprogram not used 
! Subprogram not used     refyr = 1
! Subprogram not used     refmn = 1
! Subprogram not used     refdy = 1
! Subprogram not used     refsc = 0._r8
! Subprogram not used 
! Subprogram not used     if ( timemgr_is_caltype(trim(shr_cal_gregorian))) then
! Subprogram not used        fltdy = greg2jday(year, month, day) - greg2jday(refyr,refmn,refdy)
! Subprogram not used     else ! assume no_leap (all years are 365 days)
! Subprogram not used        fltdy = (year - refyr)*365._r8 + &
! Subprogram not used                (doy(month)-doy(refmn)) + &
! Subprogram not used                (day-refdy) 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     get_time_float = fltdy + ((sec-refsc)/86400._r8)
! Subprogram not used 
! Subprogram not used   endfunction get_time_float
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... Return Julian day number given Gregorian date.
! Subprogram not used !
! Subprogram not used ! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! Subprogram not used ! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used   function greg2jday( year, month, day )
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: year, month, day
! Subprogram not used     integer :: greg2jday
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !	... Local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer :: ap, mp
! Subprogram not used     integer :: y, d, n, g
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !     	... Modify year and month numbers
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ap = year - (12 - month)/10
! Subprogram not used     mp = MOD( month-3,12 )
! Subprogram not used     if( mp < 0 ) then
! Subprogram not used        mp = mp + 12
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !     	... Julian day
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     y = INT( 365.25_r8*( ap + 4712 ) )
! Subprogram not used     d = INT( 30.6_r8*mp + .5_r8 )
! Subprogram not used     n = y + d + day  + 59
! Subprogram not used     g = INT( .75_r8*INT( ap/100 + 49 ) ) - 38
! Subprogram not used     greg2jday = n - g
! Subprogram not used 
! Subprogram not used   end function greg2jday

end module spedata
