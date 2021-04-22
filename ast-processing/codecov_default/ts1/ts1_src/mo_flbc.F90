module mo_flbc
  !---------------------------------------------------------------
  ! 	... lower boundary module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use m_types,      only : time_ramp
  use spmd_utils,   only : masterproc,iam
  use abortutils,   only : endrun
  use ioFileMod,    only : getfil
  use ppgrid,       only : pcols, begchunk, endchunk, pver
  use time_manager, only : get_curr_date, get_curr_calday
  use time_utils,   only : flt_date
  use cam_logfile,  only : iulog
  use constituents,  only : pcnst
  use constituents,  only : tracnam=>cnst_name

  implicit none

  type :: flbc
     integer            :: spc_ndx = -1
     real(r8), pointer  :: vmr(:,:,:)
     character(len=16)  :: species = ' '
     logical            :: has_mean
     real(r8), pointer  :: vmr_mean(:)
  end type flbc

  private
  public  :: flbc_inti, flbc_set, flbc_chk, has_flbc
  public  :: flbc_gmean_vmr

  save

  integer, parameter :: time_span = 1

  integer :: ntimes
  integer :: flbc_cnt
  integer :: gndx
  integer :: tim_ndx(2)
  integer :: jlim(2)
  integer, allocatable  :: dates(:)
  real(r8), allocatable     :: times(:)
  logical :: has_flbc(pcnst)
  character(len=256) :: filename, lpath, mspath

  type(time_ramp) :: flbc_timing
  integer ::  ncdate, ncsec

  integer, parameter :: nghg = 5
  integer, parameter :: max_nflbc = pcnst+nghg

  integer, parameter :: co2_ndx = 1
  integer, parameter :: ch4_ndx = 2
  integer, parameter :: n2o_ndx = 3
  integer, parameter :: f11_ndx = 4
  integer, parameter :: f12_ndx = 5
  character(len=5)  :: ghg_names(nghg) = (/ 'CO2  ','CH4  ','N2O  ','CFC11','CFC12' /)
  integer :: ghg_indices(nghg) = -1

  type(flbc) :: flbcs(max_nflbc)

  logical, parameter :: debug = .false.

contains

! Subprogram not used   subroutine flbc_inti( flbc_file, flbc_list, flbc_timing_in, co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr )
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... initialize the fixed lower bndy cond
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use mo_constants,  only : d2r, pi, rearth
! Subprogram not used     use string_utils,  only : to_upper
! Subprogram not used     use constituents,  only : cnst_get_ind
! Subprogram not used     use cam_pio_utils, only : cam_pio_openfile
! Subprogram not used     use pio,           only : pio_get_var,pio_inq_varid,pio_inq_dimid, pio_inq_dimlen
! Subprogram not used     use pio,           only : file_desc_t, pio_closefile, pio_nowrite
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... dummy arguments
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     character(len=*), intent(in) :: flbc_file
! Subprogram not used     character(len=*), intent(in) :: flbc_list(:)
! Subprogram not used     type(time_ramp),  intent(in) :: flbc_timing_in
! Subprogram not used     real(r8),         intent(in) :: co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer :: astat
! Subprogram not used     integer :: j, l, m, n                     ! Indices
! Subprogram not used     integer :: t1, t2
! Subprogram not used     type(file_desc_t) :: ncid
! Subprogram not used     integer :: dimid
! Subprogram not used     integer :: varid
! Subprogram not used     integer :: yr, mon, day, wrk_date, wrk_sec
! Subprogram not used     real(r8)    :: seq
! Subprogram not used     real(r8)    :: wrk_time
! Subprogram not used     character(len=16)  :: species
! Subprogram not used     character(len=16)  :: spc_name
! Subprogram not used     character(len=8)   :: time_type
! Subprogram not used     integer :: ierr
! Subprogram not used 
! Subprogram not used     if ( len_trim( flbc_file ) == 0 ) return
! Subprogram not used 
! Subprogram not used     call get_curr_date( yr, mon, day, ncsec )
! Subprogram not used     ncdate = yr*10000 + mon*100 + day
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... check timing
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     flbc_timing = flbc_timing_in
! Subprogram not used     time_type = to_upper(flbc_timing%type)
! Subprogram not used     flbc_timing%type = time_type
! Subprogram not used     if( time_type /= 'SERIAL' .and. time_type /= 'CYCLICAL' &
! Subprogram not used          .and. time_type /= 'FIXED' ) then
! Subprogram not used        write(iulog,*) 'flbc_inti: time type ',trim(time_type),' is not SERIAL,CYCLICAL, or FIXED'
! Subprogram not used        call endrun('flbc_inti: invalid time type ')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if ( (flbc_timing%cycle_yr>0) .and. (time_type/='CYCLICAL') ) then
! Subprogram not used        call endrun('flbc_inti: cannot specify  flbc_cycle_yr if flbc_type is not CYCLICAL')
! Subprogram not used     endif
! Subprogram not used     if ( ((flbc_timing%fixed_ymd>0).or.(flbc_timing%fixed_tod>0)).and.(time_type/='FIXED') ) then
! Subprogram not used        call endrun('flbc_inti: cannot specify  flbc_fixed_ymd or flbc_fixed_tod if flbc_type is not FIXED')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     wrk_sec  = ncsec
! Subprogram not used     if( time_type == 'SERIAL' ) then
! Subprogram not used        wrk_date = ncdate 
! Subprogram not used     else if( time_type == 'CYCLICAL' ) then
! Subprogram not used 
! Subprogram not used     	! If this is a leap-day, we have to avoid asking for a non-leap-year
! Subprogram not used     	! on a cyclical dataset. When this happens, just use Feb 28 instead
! Subprogram not used     	if (( mon .eq. 2 ) .and. ( day.eq.29 )) then
! Subprogram not used 	   ncdate = yr*10000 + mon*100 + (day-1)
! Subprogram not used            write(iulog,*)'WARNING: flbc_inti using Feb 28 instead of Feb 29 for cyclical dataset'
! Subprogram not used         endif 	
! Subprogram not used        wrk_date = flbc_timing%cycle_yr*10000 + mod(ncdate,10000)
! Subprogram not used     else
! Subprogram not used        wrk_date = flbc_timing%fixed_ymd
! Subprogram not used        wrk_sec  = flbc_timing%fixed_tod
! Subprogram not used     end if
! Subprogram not used     wrk_time = flt_date( wrk_date, wrk_sec )
! Subprogram not used     if (masterproc) write(iulog,*) 'flbc_inti: wrk_date,wrk_sec,wrk_time = ',wrk_date,wrk_sec,wrk_time
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... species with fixed lbc ?
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     has_flbc(:) = .false.
! Subprogram not used     flbc_cnt = 0
! Subprogram not used     
! Subprogram not used     do m = 1,max_nflbc
! Subprogram not used 
! Subprogram not used        if ( len_trim(flbc_list(m))==0 ) exit
! Subprogram not used 
! Subprogram not used        flbc_cnt = flbc_cnt + 1
! Subprogram not used 
! Subprogram not used        call cnst_get_ind (flbc_list(m), n, abort=.false.)
! Subprogram not used 
! Subprogram not used        if (n > 0) then
! Subprogram not used           has_flbc(n) = .true.
! Subprogram not used           flbcs(flbc_cnt)%spc_ndx = n
! Subprogram not used        else ! must be one of the GHGs which is not prognosted
! Subprogram not used           if( .not. any( ghg_names(:) == flbc_list(m) ) ) then
! Subprogram not used              call endrun('flbc_inti: flbc_list member '// trim(flbc_list(m)) //' is not allowed')
! Subprogram not used           endif
! Subprogram not used           flbcs(flbc_cnt)%spc_ndx = -1
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        flbcs(flbc_cnt)%species = trim( flbc_list(m) )
! Subprogram not used 
! Subprogram not used        where( ghg_names(:) == flbc_list(m) )
! Subprogram not used           ghg_indices = m
! Subprogram not used        endwhere
! Subprogram not used 
! Subprogram not used        if( trim(flbcs(flbc_cnt)%species) == 'CFC11' ) then
! Subprogram not used           flbcs(flbc_cnt)%species = 'CFCL3'
! Subprogram not used        elseif( trim(flbcs(flbc_cnt)%species) == 'CFC12' ) then
! Subprogram not used           flbcs(flbc_cnt)%species = 'CF2CL2'
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! check that user has not set vmr namelist values... 
! Subprogram not used     if ( ghg_indices(co2_ndx) > 0 .and. co2vmr>1.e-6_r8) then
! Subprogram not used        call endrun('flbc_inti: cannot specify both co2vmr and CO2 in flbc_file')
! Subprogram not used     endif
! Subprogram not used     if ( ghg_indices(ch4_ndx) > 0 .and. ch4vmr > 0._r8) then
! Subprogram not used        call endrun('flbc_inti: cannot specify both ch4vmr and CH4 in flbc_file')
! Subprogram not used     endif
! Subprogram not used     if ( ghg_indices(n2o_ndx) > 0 .and. n2ovmr > 0._r8) then
! Subprogram not used        call endrun('flbc_inti: cannot specify both n2ovmr and N2O in flbc_file')
! Subprogram not used     endif
! Subprogram not used     if ( ghg_indices(f11_ndx) > 0 .and. f11vmr > 0._r8) then
! Subprogram not used        call endrun('flbc_inti: cannot specify both f11vmr and CFC11 in flbc_file')
! Subprogram not used     endif
! Subprogram not used     if ( ghg_indices(f12_ndx) > 0 .and. f12vmr > 0._r8) then
! Subprogram not used        call endrun('flbc_inti: cannot specify both f12vmr and CFC12 in flbc_file')
! Subprogram not used     endif
! Subprogram not used     
! Subprogram not used     if( flbc_cnt == 0 ) then
! Subprogram not used        return
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if(masterproc) then
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used        if( flbc_cnt > 0 ) then
! Subprogram not used           write(iulog,*) 'flbc_inti: Species with specified lower boundary values'
! Subprogram not used           do n = 1,flbc_cnt
! Subprogram not used              write(iulog,*) trim(flbcs(n)%species)
! Subprogram not used           enddo
! Subprogram not used        else
! Subprogram not used           write(iulog,*) 'There are no species with specified lower boundary values'
! Subprogram not used        end if
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used 
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        ! 	... diagnostics
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used        write(iulog,*) 'flbc_inti: diagnostics'
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used        write(iulog,*) 'lower bndy timing specs'
! Subprogram not used        write(iulog,*) 'type = ',flbc_timing%type
! Subprogram not used        if( time_type == 'CYCLICAL' ) then
! Subprogram not used           write(iulog,*) 'cycle year = ',flbc_timing%cycle_yr
! Subprogram not used        else
! Subprogram not used           write(iulog,*) 'fixed date = ',flbc_timing%fixed_ymd
! Subprogram not used           write(iulog,*) 'fixed time = ',flbc_timing%fixed_tod
! Subprogram not used        end if
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used        write(iulog,*) 'there are ',flbc_cnt,' species with specified lower bndy values'
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used     end if
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... get timing information, allocate arrays, and read in dates
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     call getfil ( flbc_file, filename, 0)
! Subprogram not used     call cam_pio_openfile (ncid, trim(filename), PIO_NOWRITE)
! Subprogram not used     ierr = pio_inq_dimid( ncid, 'time', dimid )
! Subprogram not used     ierr = pio_inq_dimlen( ncid, dimid, ntimes )
! Subprogram not used 
! Subprogram not used     allocate( dates(ntimes),stat=astat )
! Subprogram not used     if( astat/= 0 ) then
! Subprogram not used        write(iulog,*) 'flbc_inti: failed to allocate dates array; error = ',astat
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     allocate( times(ntimes),stat=astat )
! Subprogram not used     if( astat/= 0 ) then
! Subprogram not used        write(iulog,*) 'flbc_inti: failed to allocate times array; error = ',astat
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ierr = pio_inq_varid( ncid, 'date', varid )
! Subprogram not used     ierr = pio_get_var( ncid, varid, dates )
! Subprogram not used 
! Subprogram not used     do n = 1,ntimes
! Subprogram not used        times(n) = flt_date( dates(n), 0 )
! Subprogram not used     end do
! Subprogram not used     if( time_type /= 'CYCLICAL' ) then
! Subprogram not used        if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
! Subprogram not used           write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        do n = 2,ntimes
! Subprogram not used           if( wrk_time <= times(n) ) then
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        tim_ndx(1) = n - 1
! Subprogram not used     else
! Subprogram not used        yr = flbc_timing%cycle_yr
! Subprogram not used        do n = 1,ntimes
! Subprogram not used           if( yr == dates(n)/10000 ) then
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if( n >= ntimes ) then
! Subprogram not used           write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        tim_ndx(1) = n
! Subprogram not used     end if
! Subprogram not used     select case( time_type )
! Subprogram not used     case( 'FIXED' )
! Subprogram not used        tim_ndx(2) = n
! Subprogram not used     case( 'CYCLICAL' )
! Subprogram not used        do n = tim_ndx(1),ntimes
! Subprogram not used           if( yr /= dates(n)/10000 ) then
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        tim_ndx(2) = n - 1
! Subprogram not used        if( (tim_ndx(2) - tim_ndx(1)) < 2 ) then
! Subprogram not used           write(iulog,*) 'flbc_inti: cyclical lb conds require at least two time points'
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used     case( 'SERIAL' )
! Subprogram not used        tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
! Subprogram not used     end select
! Subprogram not used     t1 = tim_ndx(1)
! Subprogram not used     t2 = tim_ndx(2)
! Subprogram not used 
! Subprogram not used     if( masterproc .and. debug ) then
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used        write(iulog,*) 'flbc time cnt = ',ntimes
! Subprogram not used        write(iulog,*) 'flbc times'
! Subprogram not used        write(iulog,'(10i10)') dates(:)
! Subprogram not used        write(iulog,'(1p,5g15.7)') times(:)
! Subprogram not used        write(iulog,*) 'flbc time indicies = ',tim_ndx(:)
! Subprogram not used        write(iulog,'(10i10)') dates(tim_ndx(1):tim_ndx(2))
! Subprogram not used        write(iulog,*) ' '
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     do m = 1,flbc_cnt
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        ! 	... allocate array
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        allocate( flbcs(m)%vmr(pcols,begchunk:endchunk,t1:t2),stat=astat )
! Subprogram not used        if( astat/= 0 ) then
! Subprogram not used           write(iulog,*) 'flbc_inti: failed to allocate lbc vmr; error = ',astat
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        flbcs(m)%has_mean = file_has_gmean(ncid,flbcs(m)%species)
! Subprogram not used        if ( flbcs(m)%has_mean) then
! Subprogram not used           allocate( flbcs(m)%vmr_mean(t1:t2),stat=astat )
! Subprogram not used           if( astat/= 0 ) then
! Subprogram not used              write(iulog,*) 'flbc_inti: failed to allocate lbc vmr_mean; error = ',astat
! Subprogram not used              call endrun
! Subprogram not used           end if
! Subprogram not used        endif
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        ! 	... readin the flbc vmr
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        call flbc_get( ncid, flbcs(m), .true., read_gmean=flbcs(m)%has_mean )
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! 	... close the file
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     call pio_closefile( ncid )
! Subprogram not used 
! Subprogram not used   end subroutine flbc_inti

  subroutine flbc_chk( )
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : file_desc_t, pio_closefile, pio_nowrite
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: m
    integer                     :: t1, t2, tcnt
    integer                     :: astat
    type(file_desc_t)           :: ncid
    real(r8)                        :: wrk_time
    integer ::  yr, mon, day

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    if( flbc_cnt > 0 .and. flbc_timing%type == 'SERIAL' ) then
       wrk_time = flt_date( ncdate, ncsec )
       if( wrk_time > times(tim_ndx(2)) ) then
          tcnt = tim_ndx(2) - tim_ndx(1)
          tim_ndx(1) = tim_ndx(2)
          tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
          t1 = tim_ndx(1)
          t2 = tim_ndx(2)
!!$          if( tcnt /= (t2 - t1) ) then
          !-----------------------------------------------------------------------
          ! 	... allocate array
          !-----------------------------------------------------------------------
          do m = 1,flbc_cnt
             if( associated( flbcs(m)%vmr ) ) then
                deallocate( flbcs(m)%vmr,stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'flbc_chk: failed to deallocate flbc vmr; error = ',astat
                   call endrun
                end if
             end if
             allocate( flbcs(m)%vmr(pcols,begchunk:endchunk,t1:t2),stat=astat )
             if( astat/= 0 ) then
                write(iulog,*) 'flbc_chk: failed to allocate flbc vmr; error = ',astat
                call endrun
             end if
                
             if (flbcs(m)%has_mean) then
                if( associated( flbcs(m)%vmr_mean ) ) then
                   deallocate( flbcs(m)%vmr_mean,stat=astat )
                   if( astat/= 0 ) then
                      write(iulog,*) 'flbc_chk: failed to deallocate flbc vmr; error = ',astat
                      call endrun
                   end if
                end if
                allocate( flbcs(m)%vmr_mean(t1:t2),stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'flbc_chk: failed to allocate flbc vmr; error = ',astat
                   call endrun
                end if

             endif
          end do
!!$          end if

          call cam_pio_openfile (ncid, trim(filename), PIO_NOWRITE)
          !-----------------------------------------------------------------------
          ! 	... readin the lb concentrations
          !-----------------------------------------------------------------------
          do m = 1,flbc_cnt
             call flbc_get( ncid, flbcs(m), .true., read_gmean=flbcs(m)%has_mean )
          end do

          !-----------------------------------------------------------------------
          ! 	... close the file
          !-----------------------------------------------------------------------
          call pio_closefile( ncid )

       end if
    end if

  end subroutine flbc_chk
  
  ! checks for global mean in input file
! Subprogram not used   function file_has_gmean(ncid,species)
! Subprogram not used     use pio, only : file_desc_t, pio_inq_varid, pio_noerr, pio_seterrorhandling, &
! Subprogram not used          pio_bcast_error, pio_internal_error
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type(file_desc_t),      intent(inout) :: ncid
! Subprogram not used     character(*), intent(in) :: species
! Subprogram not used     logical :: file_has_gmean
! Subprogram not used 
! Subprogram not used     integer :: varid, ierr
! Subprogram not used 
! Subprogram not used     ! Allow pio to return the potential error and handle it locally
! Subprogram not used     call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
! Subprogram not used     ierr = pio_inq_varid( ncid, trim(species)//'_LBC_mean', varid)
! Subprogram not used     call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     file_has_gmean = (ierr==PIO_NOERR)
! Subprogram not used 
! Subprogram not used   endfunction file_has_gmean

! Subprogram not used   subroutine flbc_get( ncid, lbcs, initial, read_gmean )
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... read lower bndy values
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     use mo_constants,  only : d2r, pi
! Subprogram not used     use phys_grid,     only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
! Subprogram not used     use pio,           only: file_desc_t, pio_get_var, pio_inq_varndims, &
! Subprogram not used          pio_max_name, pio_inq_varid, pio_inq_dimlen, pio_inq_dimid
! Subprogram not used     use interpolate_data, only : interp_type, lininterp_init, lininterp_finish, lininterp
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... dummy arguments
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     type(file_desc_t), intent(in)           :: ncid
! Subprogram not used     logical, intent(in)           :: initial
! Subprogram not used     type(flbc), intent(inout) :: lbcs
! Subprogram not used 
! Subprogram not used     logical, intent(in), optional :: read_gmean
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer                     :: j, m               ! Indices
! Subprogram not used     integer                     :: t1, t2, tcnt
! Subprogram not used     integer                     :: ierr
! Subprogram not used     integer                     :: vid, nlat, nlon
! Subprogram not used     integer                     :: dimid_lat, dimid_lon
! Subprogram not used     integer                     :: plon, plat
! Subprogram not used     real(r8), allocatable           :: lat(:)
! Subprogram not used     real(r8), allocatable           :: lon(:)
! Subprogram not used     real(r8), allocatable           :: wrk(:,:,:), wrk_zonal(:,:)
! Subprogram not used     real(r8), allocatable           :: wrk2d(:,:)
! Subprogram not used     character(len=pio_max_name)  :: varname
! Subprogram not used     real(r8), allocatable       :: locl_vmr(:,:,:)
! Subprogram not used     integer :: ndims, t, c, ncols
! Subprogram not used     type(interp_type) :: lon_wgts, lat_wgts
! Subprogram not used     real(r8) :: to_lats(pcols), to_lons(pcols)
! Subprogram not used     real(r8), parameter :: twopi=2._r8*pi, zero=0._r8
! Subprogram not used 
! Subprogram not used     t1 = tim_ndx(1)
! Subprogram not used     t2 = tim_ndx(2)
! Subprogram not used     tcnt = t2 - t1 + 1
! Subprogram not used     allocate( locl_vmr(pcols,begchunk:endchunk,tcnt), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'srf_emis_get: locl_emis allocation error = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     locl_vmr(:,:,:) = 0._r8
! Subprogram not used 
! Subprogram not used     initialization : if( initial ) then
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        !       ... get grid dimensions from file
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        !           latitudes
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        ierr = pio_inq_dimid( ncid, 'lat', dimid_lat )
! Subprogram not used        ierr = pio_inq_dimlen( ncid, dimid_lat, nlat )
! Subprogram not used        allocate( lat(nlat),stat=ierr )
! Subprogram not used        if( ierr /= 0 ) then
! Subprogram not used           write(iulog,*) 'flbc_get: lat allocation error = ',ierr
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        ierr = pio_inq_varid( ncid, 'lat', vid )
! Subprogram not used        ierr = pio_get_var( ncid, vid, lat )
! Subprogram not used        lat(:nlat) = lat(:nlat) * d2r
! Subprogram not used        
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        !           longitudes
! Subprogram not used        !-----------------------------------------------------------------------
! Subprogram not used        ierr = pio_inq_dimid( ncid, 'lon', dimid_lon )
! Subprogram not used        ierr = pio_inq_dimlen( ncid, dimid_lon, nlon )
! Subprogram not used        allocate( lon(nlon),stat=ierr )
! Subprogram not used        if( ierr /= 0 ) then
! Subprogram not used           write(iulog,*) 'flbc_get: lon allocation error = ',ierr
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        ierr = pio_inq_varid( ncid, 'lon', vid )
! Subprogram not used        ierr = pio_get_var( ncid, vid, lon )
! Subprogram not used        lon(:nlon) = lon(:nlon) * d2r
! Subprogram not used     end if initialization
! Subprogram not used         
! Subprogram not used     allocate( wrk(nlon,nlat,tcnt), stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'flbc_get: wrk allocation error = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !       ... read data
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     varname = trim(lbcs%species) // '_LBC'
! Subprogram not used     ierr = pio_inq_varid( ncid, trim(varname), vid )
! Subprogram not used     ierr = pio_inq_varndims (ncid, vid, ndims)
! Subprogram not used     
! Subprogram not used     if (ndims==2) then
! Subprogram not used        allocate( wrk_zonal(nlat,tcnt), stat=ierr )
! Subprogram not used        if( ierr /= 0 ) then
! Subprogram not used           write(iulog,*) 'flbc_get: wrk_zonal allocation error = ',ierr
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (ndims==2) then
! Subprogram not used        ierr = pio_get_var( ncid, vid, (/ 1, t1/), &
! Subprogram not used             (/ nlat, tcnt /), wrk_zonal )
! Subprogram not used        do t = 1,tcnt
! Subprogram not used           do j = 1,nlat
! Subprogram not used              wrk(:nlon,j,t) = wrk_zonal(j,t)
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used     else
! Subprogram not used        ierr = pio_get_var( ncid, vid, (/ 1, 1, t1/), &
! Subprogram not used             (/ nlon, nlat, tcnt /), wrk )
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     do c=begchunk,endchunk
! Subprogram not used        ncols = get_ncols_p(c)
! Subprogram not used        call get_rlat_all_p(c, pcols, to_lats)
! Subprogram not used        call get_rlon_all_p(c, pcols, to_lons)
! Subprogram not used        call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
! Subprogram not used        call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)
! Subprogram not used           
! Subprogram not used        do m = 1,tcnt
! Subprogram not used           call lininterp(wrk(:,:,m), nlon, nlat, locl_vmr(:,c,m), ncols, lon_wgts, lat_wgts) 
! Subprogram not used        end do
! Subprogram not used           
! Subprogram not used 
! Subprogram not used        call lininterp_finish(lon_wgts)
! Subprogram not used        call lininterp_finish(lat_wgts)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     deallocate(wrk, stat=ierr)
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'flbc_get: Failed to deallocate wrk, ierr = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (ndims==2) then
! Subprogram not used        deallocate( wrk_zonal,stat=ierr )
! Subprogram not used        if( ierr /= 0 ) then
! Subprogram not used           write(iulog,*) 'flbc_get: Failed to deallocate wrk_zonal, ierr = ',ierr
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if (read_gmean) then
! Subprogram not used        varname = trim(lbcs%species) // '_LBC_mean'
! Subprogram not used        ierr = pio_inq_varid( ncid, trim(varname), vid )
! Subprogram not used        ierr = pio_get_var( ncid, vid, (/t1/), (/tcnt/), lbcs%vmr_mean(t1:t2) )
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do m = t1,t2
! Subprogram not used        lbcs%vmr(:,:,m) = locl_vmr(:,:,m-t1+1)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     deallocate(locl_vmr, stat=ierr )
! Subprogram not used     if( ierr /= 0 ) then
! Subprogram not used        write(iulog,*) 'flbc_get: Failed to deallocate locl_vmr; ierr = ',ierr
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine flbc_get

  subroutine flbc_set( vmr, ncol, lchnk, map )
    !--------------------------------------------------------
    !	... set the lower bndy values
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    ::   ncol
    integer,  intent(in)    ::   lchnk
    integer,  intent(in)    ::   map(:)
    real(r8), intent(inout) ::   vmr(:,:,:)    ! lower bndy concentrations( mol/mol )

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  :: m, n
    integer  :: last, next
    real(r8) :: dels

    if( flbc_cnt < 1 ) then
       return
    end if

    call get_dels( dels, last, next )

    do m = 1,flbc_cnt
       if ( flbcs(m)%spc_ndx > 0 ) then
          n = map( flbcs(m)%spc_ndx )
          vmr(:ncol,pver,n) = flbcs(m)%vmr(:ncol,lchnk,last) &
               + dels * (flbcs(m)%vmr(:ncol,lchnk,next) - flbcs(m)%vmr(:ncol,lchnk,last))
       endif
    end do

  end subroutine flbc_set

! Subprogram not used   subroutine get_dels( dels, last, next )
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     real(r8), intent(out) :: dels
! Subprogram not used     integer,  intent(out) :: last
! Subprogram not used     integer,  intent(out) :: next
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !	... local variables
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     integer  ::  wrk_date, wrk_sec
! Subprogram not used     integer  ::  tcnt, n
! Subprogram not used     real(r8)     ::  wrk_time
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !	... setup the time interpolation
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     wrk_sec  = ncsec
! Subprogram not used     select case( flbc_timing%type )
! Subprogram not used     case( 'SERIAL' )
! Subprogram not used        wrk_date = ncdate
! Subprogram not used     case( 'CYCLICAL' )
! Subprogram not used        wrk_date = flbc_timing%cycle_yr*10000 + mod( ncdate,10000 )
! Subprogram not used     case( 'FIXED' )
! Subprogram not used        wrk_date = flbc_timing%fixed_ymd
! Subprogram not used        wrk_sec  = flbc_timing%fixed_tod
! Subprogram not used     end select
! Subprogram not used 
! Subprogram not used     wrk_time = flt_date( wrk_date, wrk_sec )
! Subprogram not used 
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     !	... set time interpolation factor
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used     if( flbc_timing%type /= 'CYCLICAL' ) then
! Subprogram not used        do n = tim_ndx(1)+1,tim_ndx(2)
! Subprogram not used           if( wrk_time <= times(n) ) then
! Subprogram not used              last = n - 1
! Subprogram not used              next = n
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if( n > ntimes ) then
! Subprogram not used           write(iulog,*) 'flbc_set: interp time is out of bounds'
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        dels = (wrk_time - times(last))/(times(next) - times(last))
! Subprogram not used        !        write(iulog,*) ' '
! Subprogram not used        !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
! Subprogram not used     else
! Subprogram not used        tcnt = tim_ndx(2) - tim_ndx(1) + 1
! Subprogram not used        call findplb( times(tim_ndx(1)), tcnt, wrk_time, n )
! Subprogram not used        if( n < tcnt ) then
! Subprogram not used           last = tim_ndx(1) + n - 1
! Subprogram not used           next = last + 1
! Subprogram not used           dels = (wrk_time - times(last))/(times(next) - times(last))
! Subprogram not used        else
! Subprogram not used           next = tim_ndx(1)
! Subprogram not used           last = tim_ndx(2)
! Subprogram not used           dels = wrk_time - times(last)
! Subprogram not used           if( dels < 0._r8 ) then
! Subprogram not used              dels = 365._r8 + dels
! Subprogram not used           end if
! Subprogram not used           dels = dels/(365._r8 + times(next) - times(last))
! Subprogram not used        end if
! Subprogram not used        !        write(iulog,*) ' '
! Subprogram not used        !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     dels = max( min( 1._r8,dels ),0._r8 )
! Subprogram not used 
! Subprogram not used   end subroutine get_dels

! Subprogram not used   subroutine flbc_gmean_vmr(co2vmr,ch4vmr,n2ovmr,f11vmr,f12vmr)
! Subprogram not used 
! Subprogram not used      implicit none
! Subprogram not used 
! Subprogram not used      real(r8), intent(inout) :: co2vmr
! Subprogram not used      real(r8), intent(inout) :: ch4vmr
! Subprogram not used      real(r8), intent(inout) :: n2ovmr
! Subprogram not used      real(r8), intent(inout) :: f11vmr
! Subprogram not used      real(r8), intent(inout) :: f12vmr
! Subprogram not used 
! Subprogram not used      integer  :: last, next
! Subprogram not used      real(r8) :: dels
! Subprogram not used 
! Subprogram not used      if( flbc_cnt < 1 ) return
! Subprogram not used 
! Subprogram not used      call get_dels( dels, last, next )
! Subprogram not used 
! Subprogram not used      if (ghg_indices(co2_ndx)>0) &
! Subprogram not used           co2vmr = global_mean_vmr(flbcs(ghg_indices(co2_ndx)), dels, last, next )
! Subprogram not used      if (ghg_indices(ch4_ndx)>0) &
! Subprogram not used           ch4vmr = global_mean_vmr(flbcs(ghg_indices(ch4_ndx)), dels, last, next )
! Subprogram not used      if (ghg_indices(n2o_ndx)>0) &
! Subprogram not used           n2ovmr = global_mean_vmr(flbcs(ghg_indices(n2o_ndx)), dels, last, next )
! Subprogram not used      if (ghg_indices(f11_ndx)>0) &
! Subprogram not used           f11vmr = global_mean_vmr(flbcs(ghg_indices(f11_ndx)), dels, last, next )
! Subprogram not used      if (ghg_indices(f12_ndx)>0) &
! Subprogram not used           f12vmr = global_mean_vmr(flbcs(ghg_indices(f12_ndx)), dels, last, next )
! Subprogram not used 
! Subprogram not used   end subroutine flbc_gmean_vmr

! Subprogram not used   function global_mean_vmr( flbcs, dels, last, next  )
! Subprogram not used     use phys_gmean, only: gmean
! Subprogram not used     use phys_grid,  only: get_ncols_p
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type(flbc), intent(in) :: flbcs
! Subprogram not used     real(r8), intent(in) :: dels
! Subprogram not used     integer, intent(in) :: last
! Subprogram not used     integer, intent(in) :: next
! Subprogram not used     real(r8) :: global_mean_vmr
! Subprogram not used     real(r8) :: vmr_arr(pcols,begchunk:endchunk)
! Subprogram not used 
! Subprogram not used     integer  :: lchnk, ncol !, n
! Subprogram not used 
! Subprogram not used     if (flbcs%has_mean) then
! Subprogram not used        global_mean_vmr = flbcs%vmr_mean(last) &
! Subprogram not used             + dels * (flbcs%vmr_mean(next) - flbcs%vmr_mean(last))
! Subprogram not used     else 
! Subprogram not used        do lchnk = begchunk, endchunk
! Subprogram not used           ncol = get_ncols_p(lchnk)
! Subprogram not used           vmr_arr(:ncol,lchnk) = flbcs%vmr(:ncol,lchnk,last) &
! Subprogram not used                + dels * (flbcs%vmr(:ncol,lchnk,next) - flbcs%vmr(:ncol,lchnk,last))
! Subprogram not used        enddo
! Subprogram not used        call gmean (vmr_arr, global_mean_vmr)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   endfunction global_mean_vmr

end module mo_flbc
