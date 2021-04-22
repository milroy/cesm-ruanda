module time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CS
   use shr_cal_mod,  only: shr_cal_noleap, shr_cal_gregorian
   use spmd_utils,   only: masterproc
   use ESMF
   use string_utils, only: to_upper
   use dycore,       only: dycore_is
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog
   use mpishorthand, only: mpicom, mpiint, mpilog, mpichar
   use pio, only : file_desc_t, var_desc_t, &
                   pio_char, pio_int, pio_global, &
                   pio_put_att, pio_def_var, pio_def_dim, &
                   pio_inq_varid, pio_put_var, pio_get_var

   implicit none
   private

! Public methods

   public ::&
      timemgr_init,             &! time manager initialization
      advance_timestep,         &! increment timestep number
      get_step_size,            &! return step size in seconds
      get_nstep,                &! return timestep number
      get_curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
      get_start_date,           &! return components of the start date
      get_ref_date,             &! return components of the reference date
      get_perp_date,            &! return components of the perpetual date, and current time of day
      get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
      get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
      get_curr_calday,          &! return calendar day at end of current timestep
      get_calday,               &! return calendar day from input date
      is_first_step,            &! return true on first step of initial run
      is_first_restart_step,    &! return true on first step of restart or branch run
      timemgr_is_caltype,       &! return true if incoming calendar type string matches actual calendar type in use
      timemgr_get_calendar_cf,  &! return cf standard for calendar type
      is_end_curr_day,          &! return true on last timestep in current day
      is_end_curr_month,        &! return true on last timestep in current month
      is_last_step,             &! return true on last timestep
      is_perpetual,             &! return true if perpetual calendar is in use
      timemgr_init_restart,     &! initialize the pio(netcdf) restart file for writting
      timemgr_write_restart,    &! write info to file needed to restart the time manager
      timemgr_read_restart,     &! read info from file needed to restart the time manager
      timemgr_restart,          &! restart the time manager
      timemgr_check_restart,    &! check that restart agrees with input clock info
      timemgr_datediff,         &! calculate difference between two time instants
      timemgr_time_ge,          &! check if time2 is later than or equal to time1
      timemgr_time_inc,         &! increment time instant by a given interval
      timemgr_set_date_time,    &! set the current date and time
      set_time_float_from_date, &! returns a float representation of time given  yr, mon, day, sec
      set_date_from_time_float   ! returns yr, mon, day, sec given time float

! interfaces


! Public data for namelist input

   integer, parameter :: uninit_int = -999999999
   integer, public ::&
      dtime         = uninit_int ! timestep in seconds

! Private module data

   ! The save attribute in the following declarations are there to make the xlf90 (v8.1.1.6)
   ! compiler happy.  It shouldn't be there since this module has a global save.

   type(ESMF_Calendar), target, save :: tm_cal        ! calendar
   type(ESMF_Clock),    save :: tm_clock      ! Model clock   
   type(ESMF_Time),     save :: tm_perp_date  ! perpetual date

   integer ::&                      ! Data required to restart time manager:
      rst_nstep     = uninit_int,  &! current step number
      rst_step_days = uninit_int,  &! days component of timestep size
      rst_step_sec  = uninit_int,  &! timestep size seconds
      rst_start_ymd = uninit_int,  &! start date
      rst_start_tod = uninit_int,  &! start time of day
      rst_stop_ymd  = uninit_int,  &! stop date
      rst_stop_tod  = uninit_int,  &! stop time of day
      rst_ref_ymd   = uninit_int,  &! reference date
      rst_ref_tod   = uninit_int,  &! reference time of day
      rst_curr_ymd  = uninit_int,  &! current date
      rst_curr_tod  = uninit_int,  &! current time of day
      rst_perp_ymd  = uninit_int    ! perpetual date
   character(len=32) :: rst_calendar  ! Calendar
   logical ::&
      rst_perp_cal  = .false.         ! true when using perpetual calendar

   character(len=32) :: calendar               ! Calendar type
   logical :: tm_first_restart_step = .false.  ! true for first step of a restart or branch run
   logical :: tm_perp_calendar = .false.       ! true when using perpetual calendar
   integer :: cal_type = uninit_int            ! calendar type

   integer, parameter :: varcnt=14
   type var_id
      type(Var_desc_t) :: vardesc
      character(len=18) :: name
   end type var_id
   type(var_id), target :: timevars(varcnt)




!=========================================================================================
contains
!=========================================================================================

subroutine timemgr_init( calendar_in, start_ymd, start_tod, ref_ymd, &
                         ref_tod, stop_ymd, stop_tod,                &
                         perpetual_run, perpetual_ymd )

! Initialize the ESMF time manager.
!
! NOTE - Assumptions:
!      1) The namelist variables have been set before this routine is called.  (set in control/parse_namelist.F90)
! Arguments
   character(len=*), intent(IN) :: calendar_in    ! Calendar type
   integer,          intent(IN) :: start_ymd      ! Start date (YYYYMMDD)
   integer,          intent(IN) :: start_tod      ! Start time of day (sec)
   integer,          intent(IN) :: ref_ymd        ! Reference date (YYYYMMDD)
   integer,          intent(IN) :: ref_tod        ! Reference time of day (sec)
   integer,          intent(IN) :: stop_ymd       ! Stop date (YYYYMMDD)
   integer,          intent(IN) :: stop_tod       ! Stop time of day (sec)
   logical,          intent(IN) :: perpetual_run  ! If in perpetual mode or not
   integer,          intent(IN) :: perpetual_ymd  ! Perpetual date (YYYYMMDD)

! Local variables
   character(len=*), parameter :: sub = 'timemgr_init'
   integer :: rc                            ! return code
   type(ESMF_Time) :: start_date            ! start date for run
   type(ESMF_Time) :: stop_date             ! stop date for run
   type(ESMF_Time) :: curr_date             ! temporary date used in logic
   type(ESMF_Time) :: ref_date              ! reference date for time coordinate
!----------------------------------------------------------------------------------------

! Initalize calendar type.

   calendar = trim(calendar_in)

   call init_calendar()

! Initalize start date.

   start_date = TimeSetymd( start_ymd, start_tod, "start_date" )

! Initalize stop date.

   stop_date = TimeSetymd( stop_ymd, stop_tod, "stop_date" )

! Initalize reference date for time coordinate.

   ref_date = TimeSetymd( ref_ymd, ref_tod, "ref_date" )

   curr_date = start_date

! Initialize clock and stop date

   call initialize_clock( start_date, ref_date, curr_date, stop_date )

! Initialize date used for perpetual calendar day calculation.

   if ( perpetual_run ) then
      tm_perp_calendar = .true.
      tm_perp_date = TimeSetymd( perpetual_ymd, 0, "tm_perp_date" )
   end if

! Print configuration summary to log file (stdout).

   if (masterproc) then
      call timemgr_print()
   end if

end subroutine timemgr_init

!=========================================================================================

subroutine initialize_clock( start_date, ref_date, curr_date, stop_date )
!
! Purpose: Initialize the clock based on the start_date, ref_date, and curr_date
! as well as the settings from the namelist specifying the time to stop
!
! Input variables
   type(ESMF_Time), intent(inout) :: start_date  ! start date for run
   type(ESMF_Time), intent(inout) :: ref_date    ! reference date for time coordinate
   type(ESMF_Time), intent(inout) :: curr_date   ! current date (equal to start_date)
   type(ESMF_Time), intent(inout) :: stop_date   ! stop date for run

! Local variables
   character(len=*), parameter :: sub = 'initialize_clock'
   type(ESMF_TimeInterval) :: step_size       ! timestep size
   type(ESMF_Time) :: current     ! current date (from clock)
   integer :: yr, mon, day, tod   ! Year, month, day, and second as integers
   integer :: rc                  ! return code

   if ( mod(86400,dtime) /= 0 ) then
!!!!      call endrun (sub//': timestep must divide evenly into 1 day')
   end if

   call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

   if ( stop_date <= start_date ) then
      write(iulog,*)sub, ': stop date must be specified later than start date: '
      call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
      call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
      call endrun
   end if
   if ( curr_date >= stop_date ) then
      write(iulog,*)sub, ': stop date must be specified later than current date: '
      call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
      call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
      call endrun
   end if

! Initialize the clock

   tm_clock = ESMF_ClockCreate(name="CAM Time-manager clock", timeStep=step_size, startTime=start_date, &
                               stopTime=stop_date, refTime=ref_date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockSetup')

! Advance clock to the current time (in case of a restart)

   call ESMF_ClockGet(tm_clock, currTime=current, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   do while( curr_date > current )
      call ESMF_ClockAdvance( tm_clock, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
      call ESMF_ClockGet(tm_clock, currTime=current )
      call chkrc(rc, sub//': error return from ESMF_ClockGet')
   end do
end subroutine initialize_clock

!=========================================================================================

function TimeSetymd( ymd, tod, desc )
!
! Set the time by an integer as YYYYMMDD and integer seconds in the day
!
   integer, intent(in) :: ymd            ! Year, month, day YYYYMMDD
   integer, intent(in) :: tod            ! Time of day in seconds
   character(len=*), intent(in) :: desc  ! Description of time to set

   type(ESMF_Time) :: TimeSetymd    ! Return value

   character(len=*), parameter :: sub = 'TimeSetymd'
   integer :: yr, mon, day          ! Year, month, day as integers
   integer :: rc                    ! return code

   if ( (ymd < 0) .or. (tod < 0) .or. (tod > 24*3600) )then
      write(iulog,*) sub//': error yymmdd is a negative number or time-of-day out of bounds', &
                 ymd, tod
      call endrun
   end if
   yr  = ymd / 10000
   mon = (ymd - yr*10000) / 100
   day =  ymd - yr*10000 - mon*100
   call ESMF_TimeSet( TimeSetymd, yy=yr, mm=mon, dd=day, s=tod, &
                      calendar=tm_cal, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeSet: setting '//trim(desc))
end function TimeSetymd

!=========================================================================================

subroutine set_time_float_from_date( time, year, month, day, sec )
!
! Set the time as a float given year, month, day, sec
!
  implicit none

  real(r8),intent(out):: time
  integer, intent(in) :: year
  integer, intent(in) :: month
  integer, intent(in) :: day
  integer, intent(in) :: sec

  integer :: rc                    ! return code
  character(len=*), parameter :: sub = 'set_time_float_from_date'
  type(ESMF_Time) :: date, ref_date
  type(ESMF_TimeInterval) :: diff
  integer :: useday

  call ESMF_TimeSet( date, yy=year, mm=month, dd=day, s=sec, calendar=tm_cal, rc=rc)
  !
  ! If the subroutine returned error, check if it is Feb 29 of a non-leap year
  ! (legitimately used by the time-interpolation routines in tracer_data.F90)
  ! in which case, substitute Feb 28 for the day
  !
  if ( rc .ne. ESMF_SUCCESS ) then
     if ( ( month .eq. 2 ) .and. ( day .eq. 29 ) ) then ! assume the failure is because it is leap day
        useday = 28
        call ESMF_TimeSet( date, yy=year, mm=month, dd=useday, s=sec, calendar=tm_cal, rc=rc)
     else  ! legitimate error, let the model quit
        call chkrc(rc, sub//': error return from ESMF_TimeSet for set_time_float_from_date')        
     endif 
  endif  

  call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockGet for set_time_float_from_date')

  diff = date - ref_date

  call ESMF_TimeIntervalGet( diff, d_r8=time, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet for set_time_float_from_date')

endsubroutine set_time_float_from_date

!=========================================================================================

! Subprogram not used subroutine set_date_from_time_float( time, year, month, day, sec )
! Subprogram not used !
! Subprogram not used ! Set year, month, day, sec given the time as a float
! Subprogram not used !
! Subprogram not used   implicit none
! Subprogram not used 
! Subprogram not used   real(r8),intent(in)  :: time
! Subprogram not used   integer, intent(out) :: year
! Subprogram not used   integer, intent(out) :: month
! Subprogram not used   integer, intent(out) :: day
! Subprogram not used   integer, intent(out) :: sec
! Subprogram not used 
! Subprogram not used   integer :: rc                    ! return code
! Subprogram not used   character(len=*), parameter :: sub = 'set_date_from_time_float'
! Subprogram not used   type(ESMF_Time) :: date, ref_date
! Subprogram not used   type(ESMF_TimeInterval) :: diff
! Subprogram not used 
! Subprogram not used   call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_ClockGet for rset_date_from_time_float')
! Subprogram not used 
! Subprogram not used   call ESMF_TimeIntervalSet( diff, d_r8=time, rc=rc)
! Subprogram not used 
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet for set_date_from_time_float')
! Subprogram not used 
! Subprogram not used   date = ref_date + diff
! Subprogram not used 
! Subprogram not used   call ESMF_TimeGet( date, yy=year, mm=month, dd=day, s=sec, rc=rc)
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_TimeGet for set_date_from_time_float')
! Subprogram not used 
! Subprogram not used endsubroutine set_date_from_time_float

!=========================================================================================

! Subprogram not used integer function TimeGetymd( date, tod )
! Subprogram not used !
! Subprogram not used ! Get the date and time of day in ymd from ESMF Time.
! Subprogram not used !
! Subprogram not used   type(ESMF_Time), intent(inout) :: date    ! Input date to convert to ymd
! Subprogram not used   integer, intent(out), optional :: tod  ! Time of day in seconds
! Subprogram not used 
! Subprogram not used   character(len=*), parameter :: sub = 'TimeGetymd'
! Subprogram not used   integer :: yr, mon, day
! Subprogram not used   integer :: rc                          ! return code
! Subprogram not used 
! Subprogram not used   call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, rc=rc)
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used   TimeGetymd = yr*10000 + mon*100 + day
! Subprogram not used   if ( present( tod ) )then
! Subprogram not used      call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used      call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used   end if
! Subprogram not used   if ( yr < 0 )then
! Subprogram not used      write(iulog,*) sub//': error year is less than zero', yr
! Subprogram not used      call endrun
! Subprogram not used   end if
! Subprogram not used end function TimeGetymd

!=========================================================================================

! Subprogram not used subroutine timemgr_set_date_time( new_ymd, new_tod )
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer,          intent(IN) :: new_ymd  ! date (YYYYMMDD)
! Subprogram not used    integer,          intent(IN) :: new_tod  ! time of day (sec)
! Subprogram not used 
! Subprogram not used    type(ESMF_Time) :: new_time   ! new time obj
! Subprogram not used    type(ESMF_Time) :: clk_time   ! clock's time obj
! Subprogram not used    integer :: rc                 ! return code
! Subprogram not used 
! Subprogram not used    character(len=*), parameter :: sub = 'timemgr_set_date_time'
! Subprogram not used 
! Subprogram not used    new_time = TimeSetymd( new_ymd, new_tod, "curr_date" )
! Subprogram not used 
! Subprogram not used ! Advance clock to the new time
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet(tm_clock, currTime=clk_time, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used    do while( new_time > clk_time )
! Subprogram not used       call ESMF_ClockAdvance( tm_clock, rc=rc )
! Subprogram not used       call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
! Subprogram not used       call ESMF_ClockGet(tm_clock, currTime=clk_time )
! Subprogram not used       call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used    end do
! Subprogram not used ! Print configuration summary to log file (stdout).
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       call timemgr_print()
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used end subroutine timemgr_set_date_time

!=========================================================================================

! Subprogram not used subroutine timemgr_restart( stop_ymd, stop_tod )
! Subprogram not used 
! Subprogram not used ! Restart the ESMF time manager.
! Subprogram not used !
! Subprogram not used ! NOTE - Assumptions:
! Subprogram not used !      1) Restart data have been read on the master process before this routine is called.
! Subprogram not used !         (timemgr_read_restart called from control/restart.F90::cam_read_restart)
! Subprogram not used !      2) Stopping time has been set and input to this routine.
! Subprogram not used ! Arguments
! Subprogram not used    integer,          intent(IN) :: stop_ymd       ! Stop date (YYYYMMDD)
! Subprogram not used    integer,          intent(IN) :: stop_tod       ! Stop time of day (sec)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'timemgr_restart'
! Subprogram not used    integer :: rc                  ! return code
! Subprogram not used    type(ESMF_Time) :: start_date  ! start date for run
! Subprogram not used    type(ESMF_Time) :: stop_date   ! stop date for run
! Subprogram not used    type(ESMF_Time) :: ref_date    ! reference date
! Subprogram not used    type(ESMF_Time) :: curr_date   ! date of data in restart file
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call mpibcast(rst_calendar,  len(rst_calendar), mpichar, 0, mpicom)
! Subprogram not used    call mpibcast(rst_step_sec,  1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_start_ymd, 1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_start_tod, 1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_stop_ymd,  1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_stop_tod,  1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_ref_ymd,   1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_ref_tod,   1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_curr_ymd,  1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_curr_tod,  1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_perp_ymd,  1, mpiint, 0, mpicom)
! Subprogram not used    call mpibcast(rst_perp_cal,  1, mpilog, 0, mpicom)
! Subprogram not used 
! Subprogram not used   calendar = trim(rst_calendar)
! Subprogram not used 
! Subprogram not used ! Initialize calendar type.
! Subprogram not used 
! Subprogram not used    call init_calendar( )
! Subprogram not used 
! Subprogram not used ! Initialize the timestep.
! Subprogram not used 
! Subprogram not used    dtime = rst_step_sec
! Subprogram not used 
! Subprogram not used ! Initialize start date.
! Subprogram not used 
! Subprogram not used    start_date = TimeSetymd( rst_start_ymd, rst_start_tod, "start_date" )
! Subprogram not used 
! Subprogram not used ! Initialize stop date.
! Subprogram not used 
! Subprogram not used    stop_date = TimeSetymd( stop_ymd, stop_tod, "stop_date" )
! Subprogram not used 
! Subprogram not used ! Initialize current date.
! Subprogram not used 
! Subprogram not used    curr_date = TimeSetymd( rst_curr_ymd, rst_curr_tod, "curr_date" )
! Subprogram not used 
! Subprogram not used ! Initialize ref date.
! Subprogram not used 
! Subprogram not used    ref_date = TimeSetymd( rst_ref_ymd, rst_ref_tod, "ref_date" )
! Subprogram not used 
! Subprogram not used ! Initialize clock and the stop date
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call initialize_clock( start_date, ref_date, curr_date, stop_date )
! Subprogram not used 
! Subprogram not used ! Advance the timestep.  Data from the restart file corresponds to the
! Subprogram not used ! last timestep of the previous run.
! Subprogram not used 
! Subprogram not used    call advance_timestep()
! Subprogram not used 
! Subprogram not used !  Set flag that this is the first timestep of the restart run.
! Subprogram not used 
! Subprogram not used    tm_first_restart_step = .true.
! Subprogram not used 
! Subprogram not used ! Initialize date used for perpetual calendar day calculation.
! Subprogram not used 
! Subprogram not used    if ( rst_perp_cal ) then
! Subprogram not used       tm_perp_date = TimeSetymd( rst_perp_ymd, 0, "tm_perp_date" )
! Subprogram not used       tm_perp_calendar = .true.
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used ! Print configuration summary to log file (stdout).
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       call timemgr_print()
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used end subroutine timemgr_restart

!=========================================================================================

! Subprogram not used subroutine timemgr_check_restart( calendar_in, start_ymd, start_tod, ref_ymd, &
! Subprogram not used                                   ref_tod, dtime_in, perpetual_run, perpetual_ymd )
! Subprogram not used 
! Subprogram not used ! Check that time-manager restart agrees with input clock information primitives.
! Subprogram not used !
! Subprogram not used ! Arguments
! Subprogram not used    character(len=*), intent(IN) :: calendar_in    ! Calendar type
! Subprogram not used    integer,          intent(IN) :: start_ymd      ! Start date (YYYYMMDD)
! Subprogram not used    integer,          intent(IN) :: start_tod      ! Start time of day (sec)
! Subprogram not used    integer,          intent(IN) :: ref_ymd        ! Reference date (YYYYMMDD)
! Subprogram not used    integer,          intent(IN) :: ref_tod        ! Reference time of day (sec)
! Subprogram not used    integer,          intent(IN) :: dtime_in       ! Time-step
! Subprogram not used    logical,          intent(IN) :: perpetual_run  ! If in perpetual mode or not
! Subprogram not used    integer,          intent(IN) :: perpetual_ymd  ! Perpetual date (YYYYMMDD)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'timemgr_check_restart'
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   ! Check that input agrees with data on restart file
! Subprogram not used 
! Subprogram not used   if ( (rst_start_ymd /= start_ymd) .or. (rst_start_tod /= start_tod) )then
! Subprogram not used      write(iulog,*) sub,' start ',rst_start_ymd,start_ymd,rst_start_tod,start_tod
! Subprogram not used      call endrun( sub//': input start date does not agree with restart' )
! Subprogram not used   end if
! Subprogram not used   if ( (rst_ref_ymd /= ref_ymd) .or. (rst_ref_tod /= ref_tod) )then
! Subprogram not used      write(iulog,*) sub,' ref ',rst_ref_ymd,ref_ymd,rst_ref_tod,ref_tod
! Subprogram not used      call endrun( sub//': input ref date does not agree with restart' )
! Subprogram not used   end if
! Subprogram not used   if ( rst_perp_cal .neqv. perpetual_run )then
! Subprogram not used      call endrun( sub//': input perpetual mode does not agree with restart' )
! Subprogram not used   end if
! Subprogram not used   if ( rst_step_sec /= dtime_in )then
! Subprogram not used      call endrun( sub//': input dtime does not agree with restart' )
! Subprogram not used   end if
! Subprogram not used   if ( trim(rst_calendar) /= trim(calendar_in) )then
! Subprogram not used      write(iulog,*) 'Input calendar:   ', trim(calendar_in)
! Subprogram not used      write(iulog,*) 'Restart calendar: ', trim(rst_calendar)
! Subprogram not used      call endrun( sub//': input calendar does not agree with restart' )
! Subprogram not used   end if
! Subprogram not used   if ( perpetual_run )then
! Subprogram not used      if ( (rst_perp_ymd /= perpetual_ymd) )then
! Subprogram not used         call endrun( sub//': input perpetual date does not agree with restart' )
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used   
! Subprogram not used end subroutine timemgr_check_restart

!=========================================================================================

subroutine init_calendar( )
!
! Initialize calendar
!
! Local variables
   character(len=*), parameter :: sub = 'init_calendar'
   type(ESMF_CalKind_Flag) :: cal_type        ! calendar type
   character(len=len(calendar)) :: caltmp
   integer :: rc                              ! return code

   caltmp = to_upper(trim(calendar) )
   if ( trim(caltmp) == trim(shr_cal_noleap) ) then
      cal_type = ESMF_CALKIND_NOLEAP
   else if ( trim(caltmp) == trim(shr_cal_gregorian) ) then
      cal_type = ESMF_CALKIND_GREGORIAN
   else
      write(iulog,*)sub,': unrecognized calendar specified: ',calendar
      call endrun
   end if
   tm_cal = ESMF_CalendarCreate( name=caltmp, calkindflag=cal_type, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_CalendarSet')
end subroutine init_calendar

!=========================================================================================

subroutine timemgr_print()

! Local variables
   character(len=*), parameter :: sub = 'timemgr_print'
   integer :: rc
   integer :: yr, mon, day
   integer ::&                  ! Data required to restart time manager:
      nstep     = uninit_int,  &! current step number
      step_sec  = uninit_int,  &! timestep size seconds
      start_yr  = uninit_int,  &! start year
      start_mon = uninit_int,  &! start month
      start_day = uninit_int,  &! start day of month
      start_tod = uninit_int,  &! start time of day
      stop_yr   = uninit_int,  &! stop year
      stop_mon  = uninit_int,  &! stop month
      stop_day  = uninit_int,  &! stop day of month
      stop_tod  = uninit_int,  &! stop time of day
      ref_yr    = uninit_int,  &! reference year
      ref_mon   = uninit_int,  &! reference month
      ref_day   = uninit_int,  &! reference day of month
      ref_tod   = uninit_int,  &! reference time of day
      curr_yr   = uninit_int,  &! current year
      curr_mon  = uninit_int,  &! current month
      curr_day  = uninit_int,  &! current day of month
      curr_tod  = uninit_int    ! current time of day
   integer(ESMF_KIND_I8) :: step_no
   type(ESMF_Time) :: start_date! start date for run
   type(ESMF_Time) :: stop_date ! stop date for run
   type(ESMF_Time) :: curr_date ! date of data in restart file
   type(ESMF_Time) :: ref_date  ! reference date
   type(ESMF_TimeInterval) :: step ! Time-step
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, &
                       refTime=ref_date, stopTime=stop_date, timeStep=step, &
                       advanceCount=step_no, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   nstep = step_no

   write(iulog,*)' ********** Time Manager Configuration **********'

   call ESMF_TimeIntervalGet( step, s=step_sec, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

   call ESMF_TimeGet( start_date, yy=start_yr, mm=start_mon, dd=start_day, &
                      s=start_tod, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeGet( stop_date, yy=stop_yr, mm=stop_mon, dd=stop_day, &
                      s=stop_tod, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeGet( ref_date, yy=ref_yr, mm=ref_mon, dd=ref_day, s=ref_tod, &
                      rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeGet( curr_date, yy=curr_yr, mm=curr_mon, dd=curr_day, &
                      s=curr_tod, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

   write(iulog,*)'  Calendar type:            ',trim(calendar)
   write(iulog,*)'  Timestep size (seconds):  ', step_sec
   write(iulog,*)'  Start date (yr mon day tod):     ', start_yr, start_mon, &
                start_day, start_tod
   write(iulog,*)'  Stop date (yr mon day tod):      ', stop_yr, stop_mon, &
                stop_day, stop_tod
   write(iulog,*)'  Reference date (yr mon day tod): ', ref_yr, ref_mon, &
                ref_day, ref_tod
   write(iulog,*)'  Current step number:      ', nstep
   write(iulog,*)'  Current date (yr mon day tod):   ', curr_yr, curr_mon, &
                curr_day, curr_tod

   if ( tm_perp_calendar ) then
      call ESMF_TimeGet( tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_TimeGet')
      write(iulog,*)'  Use perpetual diurnal cycle date (yr mon day): ', &
                yr, mon, day
   end if

   write(iulog,*)' ************************************************'

end subroutine timemgr_print
!=========================================================================================

subroutine advance_timestep()

! Increment the timestep number.

! Local variables
   character(len=*), parameter :: sub = 'advance_timestep'
   integer :: rc
!-----------------------------------------------------------------------------------------

   call ESMF_ClockAdvance( tm_clock, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockAdvance')

! Set first step flag off.

   tm_first_restart_step = .false.

end subroutine advance_timestep
!=========================================================================================

integer function get_step_size()

! Return the step size in seconds.

! Local variables
   character(len=*), parameter :: sub = 'get_step_size'
   type(ESMF_TimeInterval) :: step_size       ! timestep size
   integer :: rc
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, timeStep=step_size, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   call ESMF_TimeIntervalGet(step_size, s=get_step_size, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')

end function get_step_size
!=========================================================================================

integer function get_nstep()

! Return the timestep number.

! Local variables
   character(len=*), parameter :: sub = 'get_nstep'
   integer :: rc
   integer(ESMF_KIND_I8) :: step_no
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, advanceCount=step_no, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   get_nstep = step_no

end function get_nstep
!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod, offset)

! Return date components valid at end of current timestep with an optional
! offset (positive or negative) in seconds.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.

! Local variables
   character(len=*), parameter :: sub = 'get_curr_date'
   integer :: rc
   type(ESMF_Time) :: date
   type(ESMF_TimeInterval) :: off
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   if (present(offset)) then
      if (offset > 0) then
         call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date + off
      else if (offset < 0) then
         call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date - off
      end if
   end if

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_curr_date
!=========================================================================================

! Subprogram not used subroutine get_perp_date(yr, mon, day, tod, offset)
! Subprogram not used 
! Subprogram not used ! Return time of day valid at end of current timestep and the components
! Subprogram not used ! of the perpetual date (with an optional offset (positive or negative) in seconds.
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
! Subprogram not used                                             ! Positive for future times, negative 
! Subprogram not used                                             ! for previous times.
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'get_perp_date'
! Subprogram not used    integer :: rc
! Subprogram not used    type(ESMF_Time) :: date
! Subprogram not used    type(ESMF_TimeInterval) :: DelTime
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
! Subprogram not used    ! Get time of day add it to perpetual date
! Subprogram not used    ! Get year, month, day so that seconds are time-of-day rather than since start time
! Subprogram not used    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used    call ESMF_TimeIntervalSet(DelTime, s=tod, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
! Subprogram not used    date = tm_perp_date + DelTime
! Subprogram not used    if ( present(offset) )then
! Subprogram not used       call ESMF_TimeIntervalSet(DelTime, s=offset, rc=rc)
! Subprogram not used       call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
! Subprogram not used       date = date + DelTime
! Subprogram not used    end if
! Subprogram not used    ! Get time of day from the result
! Subprogram not used    ! Get year, month, day so that seconds are time-of-day rather than since start time
! Subprogram not used    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used 
! Subprogram not used    ! Get the date from the fixed perpetual date (in case it overflows to next day)
! Subprogram not used    call ESMF_TimeGet(tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used 
! Subprogram not used end subroutine get_perp_date
!=========================================================================================

! Subprogram not used subroutine get_prev_date(yr, mon, day, tod)
! Subprogram not used 
! Subprogram not used ! Return date components valid at beginning of current timestep.
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'get_prev_date'
! Subprogram not used    integer :: rc
! Subprogram not used    type(ESMF_Time) :: date
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used end subroutine get_prev_date
!=========================================================================================

! Subprogram not used subroutine get_start_date(yr, mon, day, tod)
! Subprogram not used 
! Subprogram not used ! Return date components valid at beginning of initial run.
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'get_start_date'
! Subprogram not used    integer :: rc
! Subprogram not used    type(ESMF_Time) :: date
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet(tm_clock, startTime=date, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used end subroutine get_start_date
!=========================================================================================

subroutine get_ref_date(yr, mon, day, tod)

! Return date components of the reference date.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_ref_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, refTime=date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_ref_date
!=========================================================================================

subroutine get_curr_time(days, seconds)

! Return time components valid at end of current timestep.
! Current time is the time interval between the current date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'get_curr_time'
   integer :: rc
   type(ESMF_Time) :: cdate, rdate
   type(ESMF_TimeInterval) :: diff
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=cdate, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_ClockGet( tm_clock, refTime=rdate, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   diff = cdate - rdate

   call ESMF_TimeIntervalGet(diff, d=days, s=seconds, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

end subroutine get_curr_time
!=========================================================================================

! Subprogram not used subroutine get_prev_time(days, seconds)
! Subprogram not used 
! Subprogram not used ! Return time components valid at beg of current timestep.
! Subprogram not used ! prev time is the time interval between the prev date and the reference date.
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       days,   &! number of whole days in time interval
! Subprogram not used       seconds  ! remaining seconds in time interval
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'get_prev_time'
! Subprogram not used    integer :: rc
! Subprogram not used    type(ESMF_Time) :: date, ref_date
! Subprogram not used    type(ESMF_TimeInterval) :: diff
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet for prevTime')
! Subprogram not used    call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet for refTime')
! Subprogram not used    diff = date - ref_date
! Subprogram not used    call ESMF_TimeIntervalGet( diff, d=days, s=seconds, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeintervalGet')
! Subprogram not used 
! Subprogram not used end subroutine get_prev_time
!=========================================================================================

function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value
   real(r8) :: get_curr_calday

! Local variables
   character(len=*), parameter :: sub = 'get_curr_calday'
   integer :: rc
   type(ESMF_Time) :: date
   type(ESMF_TimeInterval) :: off, diurnal
   integer :: year, month, day, tod
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   if (present(offset)) then
      if (offset > 0) then
         call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date + off
      else if (offset < 0) then
         call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date - off
      end if
   end if

   if ( tm_perp_calendar ) then
!     Get current time-of-day from clock
      call ESMF_TimeGet(date, yy=year, mm=month, dd=day, s=tod, rc=rc)
      call chkrc(rc, sub//': error return from ESMF_TimeGet')
!     Get date from perpetual date add time-of-day to it
      call ESMF_TimeIntervalSet( diurnal, s=tod, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
      date = tm_perp_date + diurnal
!!!!  write(iulog,*) ' tod = ', tod
!!!!  call ESMF_TimePrint( date, "string" )
   end if

   call ESMF_TimeGet( date, dayOfYear_r8=get_curr_calday, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

!
!                  WARNING: Gregorian calendar fakes day 366
!
! The zenith angle calculation is only capable of using a 365-day calendar.
! If a Gregorian calendar is being used, the last day of a leap year (day 366)
! is sent to the model as a repetition of the previous day (day 365). 
! This is done by decrementing calday by 1 immediately below.
! bundy, July 2008
!
   if (( get_curr_calday > 366.0_r8 ) .and. ( get_curr_calday <= 367.0_r8 ) &
        .and. (timemgr_is_caltype(trim(shr_cal_gregorian)))) then
      get_curr_calday = get_curr_calday - 1.0_r8
   endif

   if ( (get_curr_calday < 1.0_r8) .or. (get_curr_calday > 366.0_r8) )then
      write(iulog,*) 'atm '//sub//' calday = ', get_curr_calday
      if ( present(offset) ) write(iulog,*) 'offset = ', offset
      call endrun( sub//': error get_curr_calday out of bounds' )
   end if

end function get_curr_calday
!=========================================================================================

function get_calday(ymd, tod)

! Return calendar day corresponding to specified time instant.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, intent(in) :: &
      ymd,   &! date in yearmmdd format
      tod     ! time of day (seconds past 0Z)

! Return value
   real(r8) :: get_calday

! Local variables
   character(len=*), parameter :: sub = 'get_calday'
   integer :: rc                 ! return code
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   date = TimeSetymd( ymd, tod, "get_calday" )
   call ESMF_TimeGet( date, dayOfYear_r8=get_calday, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

!
!                  WARNING: Gregorian calendar fakes day 366
!
! The zenith angle calculation is only capable of using a 365-day calendar.
! If a Gregorian calendar is being used, the last day of a leap year (day 366)
! is sent to the model as a repetition of the previous day (day 365). 
! This is done by decrementing calday by 1 immediately below.
! bundy, July 2008
!
   if (( get_calday > 366.0_r8 ) .and. ( get_calday <= 367.0_r8 ) &
        .and. (timemgr_is_caltype(trim(shr_cal_gregorian)))) then
      get_calday = get_calday - 1.0_r8
   endif

   if ( (get_calday < 1.0_r8) .or. (get_calday > 366.0_r8) )then
      write(iulog,*) 'atm '//sub//' calday = ', get_calday
      call endrun( sub//': error calday out of range' )
   end if

end function get_calday
!=========================================================================================

character(len=SHR_KIND_CS) function timemgr_get_calendar_cf()

! Return cf standard for calendar type

! Local variables
   character(len=*), parameter  :: sub = 'timemgr_get_calendar_cf'
   character(len=len(calendar)) :: caltmp
!-----------------------------------------------------------------------------------------

   caltmp = to_upper(trim(calendar) )
   if ( trim(caltmp) == trim(shr_cal_noleap) ) then
      timemgr_get_calendar_cf = 'noleap'
   else if ( trim(caltmp) == trim(shr_cal_gregorian) ) then
      timemgr_get_calendar_cf = 'gregorian'
   else
      write(iulog,*)sub,': unrecognized calendar specified: ',calendar
      call endrun
   end if

end function timemgr_get_calendar_cf
!=========================================================================================
 
function timemgr_is_caltype( cal_in )

! Return true if incoming calendar type string matches actual calendar type in use

   character(len=*), intent(in) :: cal_in

! Return value
   logical :: timemgr_is_caltype

!-----------------------------------------------------------------------------------------

   timemgr_is_caltype = ( to_upper(trim(calendar)) == to_upper(trim(cal_in)) )

end function timemgr_is_caltype
!=========================================================================================
 
function is_end_curr_day()

! Return true if current timestep is last timestep in current day.

! Return value
   logical :: is_end_curr_day

! Local variables
   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
!-----------------------------------------------------------------------------------------

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_day = (tod == 0)

end function is_end_curr_day
!=========================================================================================

! Subprogram not used logical function is_end_curr_month()
! Subprogram not used 
! Subprogram not used ! Return true if current timestep is last timestep in current month.
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    integer ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call get_curr_date(yr, mon, day, tod)
! Subprogram not used    is_end_curr_month = (day == 1  .and.  tod == 0)
! Subprogram not used 
! Subprogram not used end function is_end_curr_month
!=========================================================================================

logical function is_first_step()

! Return true on first step of initial run only.

! Local variables
   character(len=*), parameter :: sub = 'is_first_step'
   integer :: rc
   integer :: nstep
   integer(ESMF_KIND_I8) :: step_no
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, advanceCount=step_no, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   nstep = step_no
   is_first_step = (nstep == 0)

end function is_first_step
!=========================================================================================

logical function is_first_restart_step()

! Return true on first step of restart run only.

!-----------------------------------------------------------------------------------------

   is_first_restart_step = tm_first_restart_step

end function is_first_restart_step
!=========================================================================================

! Subprogram not used logical function is_last_step()
! Subprogram not used 
! Subprogram not used ! Return true on last timestep.
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'is_last_step'
! Subprogram not used    type(ESMF_Time) :: stop_date
! Subprogram not used    type(ESMF_Time) :: curr_date
! Subprogram not used    type(ESMF_TimeInterval) :: time_step
! Subprogram not used    integer :: rc
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet( tm_clock, stopTime=stop_date, &
! Subprogram not used                        currTime=curr_date, TimeStep=time_step, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used    if ( curr_date+time_step > stop_date ) then
! Subprogram not used       is_last_step = .true.
! Subprogram not used    else
! Subprogram not used       is_last_step = .false.
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used end function is_last_step
!=========================================================================================

! Subprogram not used logical function is_perpetual()
! Subprogram not used 
! Subprogram not used ! Return true on last timestep.
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    is_perpetual = tm_perp_calendar
! Subprogram not used 
! Subprogram not used end function is_perpetual

!=========================================================================================

! Subprogram not used subroutine timemgr_init_restart(File)
! Subprogram not used 
! Subprogram not used   type(File_desc_t) :: File
! Subprogram not used   integer :: i, ret, dims(1)
! Subprogram not used 
! Subprogram not used   call timevars_set_names()
! Subprogram not used   rst_calendar  = trim(calendar)
! Subprogram not used 
! Subprogram not used   ret = PIO_DEF_dim(File, 'cal_strlen', 32, dims(1))
! Subprogram not used   ret = PIO_def_var(File, timevars(1)%name, PIO_CHAR,dims(1:1),timevars(1)%vardesc)
! Subprogram not used   do i=2,varcnt
! Subprogram not used      ret = PIO_def_var(File, timevars(i)%name, PIO_INT, dims(1:0), timevars(i)%vardesc)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine timemgr_init_restart

! Subprogram not used subroutine timemgr_write_restart(File)
! Subprogram not used 
! Subprogram not used   type(File_desc_t) :: File
! Subprogram not used   type(var_id), pointer :: timevar
! Subprogram not used   integer :: rst_perp_cal_int
! Subprogram not used   integer :: i, ret
! Subprogram not used   character(len=*), parameter :: sub = 'timemgr_write_nfrestart'
! Subprogram not used   type(ESMF_Time) :: start_date            ! start date for run
! Subprogram not used   type(ESMF_Time) :: stop_date             ! stop date for run
! Subprogram not used   type(ESMF_Time) :: curr_date             ! temporary date used in logic
! Subprogram not used   type(ESMF_Time) :: ref_date              ! reference date for time coordinate
! Subprogram not used 
! Subprogram not used   if ( tm_perp_calendar ) then
! Subprogram not used      rst_perp_ymd = TimeGetymd( tm_perp_date )
! Subprogram not used      rst_perp_cal = tm_perp_calendar
! Subprogram not used   else
! Subprogram not used      rst_perp_cal = .false.
! Subprogram not used   end if
! Subprogram not used   
! Subprogram not used   call ESMF_ClockGet( tm_clock, startTime=start_date, stopTime=stop_date, &
! Subprogram not used        currTime=curr_date, refTime=ref_date, rc=ret )
! Subprogram not used   call chkrc(ret, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used   rst_step_sec  = dtime
! Subprogram not used   rst_start_ymd = TimeGetymd( start_date, tod=rst_start_tod )
! Subprogram not used   rst_stop_ymd  = TimeGetymd( stop_date,  tod=rst_stop_tod  )
! Subprogram not used   rst_ref_ymd   = TimeGetymd( ref_date,   tod=rst_ref_tod   )
! Subprogram not used   rst_curr_ymd  = TimeGetymd( curr_date,  tod=rst_curr_tod  )
! Subprogram not used   
! Subprogram not used   if ( rst_perp_cal ) then
! Subprogram not used      rst_perp_cal_int = 1
! Subprogram not used   else
! Subprogram not used      rst_perp_cal_int = 0
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   do i=1,varcnt
! Subprogram not used      timevar => timevars(i)
! Subprogram not used      if(timevar%name .eq. 'rst_calendar') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_calendar)
! Subprogram not used      else if(timevar%name == 'rst_nstep') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_nstep)
! Subprogram not used      else if(timevar%name == 'rst_step_days') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_step_days)
! Subprogram not used      else if(timevar%name == 'rst_step_sec') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_step_sec)
! Subprogram not used      else if(timevar%name == 'rst_start_ymd') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_start_ymd)
! Subprogram not used      else if(timevar%name == 'rst_start_tod') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_start_tod)
! Subprogram not used      else if(timevar%name == 'rst_stop_ymd') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_stop_ymd)
! Subprogram not used      else if(timevar%name == 'rst_stop_tod') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_stop_tod)
! Subprogram not used      else if(timevar%name == 'rst_ref_ymd') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_ref_ymd)
! Subprogram not used      else if(timevar%name == 'rst_ref_tod') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_ref_tod)
! Subprogram not used      else if(timevar%name == 'rst_curr_ymd') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_curr_ymd)
! Subprogram not used      else if(timevar%name == 'rst_curr_tod') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_curr_tod)
! Subprogram not used      else if(timevar%name == 'rst_perp_ymd') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_perp_ymd)
! Subprogram not used      else if(timevar%name == 'rst_perp_cal_int') then
! Subprogram not used         ret = PIO_put_var(File, timevar%vardesc%varID, rst_perp_cal_int)
! Subprogram not used      else
! Subprogram not used         write(iulog,*) "<stdin>",1288,'Could not find a match for ', i, timevar%name
! Subprogram not used         call endrun('time_manager')
! Subprogram not used      end if
! Subprogram not used   end do
! Subprogram not used end subroutine timemgr_write_restart

! Subprogram not used subroutine timemgr_read_restart(File)
! Subprogram not used 
! Subprogram not used   type(File_desc_t) :: File
! Subprogram not used   type(var_id), pointer :: timevar
! Subprogram not used   integer :: ret, i
! Subprogram not used   integer :: rst_perp_cal_int
! Subprogram not used 
! Subprogram not used   call timevars_set_names()
! Subprogram not used   
! Subprogram not used   do i=1,varcnt
! Subprogram not used      ret = PIO_inq_varid(File, timevars(i)%name, timevars(i)%vardesc)
! Subprogram not used   end do
! Subprogram not used   
! Subprogram not used   do i=1,varcnt
! Subprogram not used      timevar => timevars(i)
! Subprogram not used      if(timevar%name == 'rst_calendar') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_calendar)
! Subprogram not used         if(rst_calendar(1:7).eq.trim(shr_cal_noleap)) then
! Subprogram not used            rst_calendar(8:)=' '
! Subprogram not used         else if(rst_calendar(1:9).eq.trim(shr_cal_gregorian)) then
! Subprogram not used            rst_calendar(10:)=' '
! Subprogram not used         end if
! Subprogram not used      else if(timevar%name == 'rst_nstep') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_nstep)
! Subprogram not used      else if(timevar%name == 'rst_step_days') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_step_days)
! Subprogram not used      else if(timevar%name == 'rst_step_sec') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_step_sec)
! Subprogram not used      else if(timevar%name == 'rst_start_ymd') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_start_ymd)
! Subprogram not used      else if(timevar%name == 'rst_start_tod') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_start_tod)
! Subprogram not used      else if(timevar%name == 'rst_stop_ymd') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_stop_ymd)
! Subprogram not used      else if(timevar%name == 'rst_stop_tod') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_stop_tod)
! Subprogram not used      else if(timevar%name == 'rst_ref_ymd') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_ref_ymd)
! Subprogram not used      else if(timevar%name == 'rst_ref_tod') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_ref_tod)
! Subprogram not used      else if(timevar%name == 'rst_curr_ymd') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_curr_ymd)
! Subprogram not used      else if(timevar%name == 'rst_curr_tod') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_curr_tod)
! Subprogram not used      else if(timevar%name == 'rst_perp_ymd') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_perp_ymd)
! Subprogram not used      else if(timevar%name == 'rst_perp_cal_int') then
! Subprogram not used         ret = PIO_get_var(File, timevar%vardesc%varID, rst_perp_cal_int)
! Subprogram not used      else
! Subprogram not used         write(iulog,*) "<stdin>",1343,'Could not find a match for ', i, timevar%name
! Subprogram not used         call endrun('time_manager')
! Subprogram not used      end if
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   if ( rst_perp_cal_int /= 0 ) then
! Subprogram not used      rst_perp_cal = .true.
! Subprogram not used   else
! Subprogram not used      rst_perp_cal = .false.
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used end subroutine timemgr_read_restart



! Subprogram not used subroutine timevars_set_names()
! Subprogram not used   timevars(1)%name = 'rst_calendar'
! Subprogram not used   timevars(2)%name = 'rst_nstep'
! Subprogram not used   timevars(3)%name = 'rst_step_days'
! Subprogram not used   timevars(4)%name = 'rst_step_sec'
! Subprogram not used   timevars(5)%name = 'rst_start_ymd'
! Subprogram not used   timevars(6)%name = 'rst_start_tod'
! Subprogram not used   timevars(7)%name = 'rst_stop_ymd'
! Subprogram not used   timevars(8)%name = 'rst_stop_tod'
! Subprogram not used   timevars(9)%name = 'rst_ref_ymd'
! Subprogram not used   timevars(10)%name = 'rst_ref_tod'
! Subprogram not used   timevars(11)%name = 'rst_curr_ymd'
! Subprogram not used   timevars(12)%name = 'rst_curr_tod'
! Subprogram not used   timevars(13)%name = 'rst_perp_ymd'
! Subprogram not used   timevars(14)%name = 'rst_perp_cal_int'
! Subprogram not used end subroutine timevars_set_names

! Subprogram not used subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)
! Subprogram not used 
! Subprogram not used ! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(in) ::&
! Subprogram not used       ymd1,    &! date1 in yyyymmdd format
! Subprogram not used       tod1,    &! time of day relative to date1 (seconds past 0Z)
! Subprogram not used       ymd2,    &! date2 in yyyymmdd format
! Subprogram not used       tod2      ! time of day relative to date2 (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    real(r8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'timemgr_datediff'
! Subprogram not used    integer :: rc   ! return code
! Subprogram not used 
! Subprogram not used    type(ESMF_Time) :: date1
! Subprogram not used    type(ESMF_Time) :: date2
! Subprogram not used    type(ESMF_TimeInterval) :: diff
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    date1 = TimeSetymd( ymd1, tod1, "date1" )
! Subprogram not used    date2 = TimeSetymd( ymd2, tod2, "date2" )
! Subprogram not used    diff = date2 - date1
! Subprogram not used    call ESMF_TimeIntervalGet( diff, d_r8=days, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')
! Subprogram not used 
! Subprogram not used end subroutine timemgr_datediff
!=========================================================================================

! Subprogram not used subroutine timemgr_time_ge(ymd1, tod1, ymd2, tod2, time2_ge_time1)
! Subprogram not used 
! Subprogram not used ! time2_ge_time1 is set to true if (ymd2,tod2) is later than or equal to (ymd1,tod1)
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(in) ::&
! Subprogram not used       ymd1,    &! date1 in yyyymmdd format
! Subprogram not used       tod1,    &! time of day relative to date1 (seconds past 0Z)
! Subprogram not used       ymd2,    &! date2 in yyyymmdd format
! Subprogram not used       tod2      ! time of day relative to date2 (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    logical :: time2_ge_time1
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'timemgr_time_ge'
! Subprogram not used    integer :: rc   ! return code
! Subprogram not used 
! Subprogram not used    type(ESMF_Time) :: time1, time2
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    time1 = TimeSetymd( ymd1, tod1, "date1" )
! Subprogram not used    time2 = TimeSetymd( ymd2, tod2, "date2" )
! Subprogram not used    time2_ge_time1 = (time2 >= time1)
! Subprogram not used 
! Subprogram not used end subroutine timemgr_time_ge

!=========================================================================================

! Subprogram not used subroutine timemgr_time_inc(ymd1, tod1, ymd2, tod2, inc_s, inc_h, inc_d)
! Subprogram not used 
! Subprogram not used ! Increment the time instant (ymd1,tod1) by an interval and return the resulting
! Subprogram not used ! time instant (ymd2,tod2).
! Subprogram not used 
! Subprogram not used    ! Arguments
! Subprogram not used    integer, intent(in) ::&
! Subprogram not used       ymd1,    &! date1 in yyyymmdd format
! Subprogram not used       tod1      ! time of day relative to date1 (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       ymd2,    &! date2 in yyyymmdd format
! Subprogram not used       tod2      ! time of day relative to date2 (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    integer, intent(in), optional ::&
! Subprogram not used       inc_s,   &! number of seconds in interval
! Subprogram not used       inc_h,   &! number of hours in interval
! Subprogram not used       inc_d     ! number of days in interval
! Subprogram not used 
! Subprogram not used    ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'timemgr_time_inc'
! Subprogram not used    integer :: rc   ! return code
! Subprogram not used 
! Subprogram not used    type(ESMF_Time) :: date1
! Subprogram not used    type(ESMF_Time) :: date2
! Subprogram not used    type(ESMF_TimeInterval) :: t_interval
! Subprogram not used    integer :: year, month, day
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ! set esmf time object
! Subprogram not used    date1 = TimeSetymd( ymd1, tod1, "date1" )
! Subprogram not used 
! Subprogram not used    ! set esmf time interval object
! Subprogram not used    if (present(inc_s)) then
! Subprogram not used       call ESMF_TimeIntervalSet(t_interval, s=inc_s, rc=rc)
! Subprogram not used    else if (present(inc_h)) then
! Subprogram not used       call ESMF_TimeIntervalSet(t_interval, h=inc_h, rc=rc)
! Subprogram not used    else if (present(inc_d)) then
! Subprogram not used       call ESMF_TimeIntervalSet(t_interval, d=inc_d, rc=rc)
! Subprogram not used    else
! Subprogram not used       call endrun(sub//': one of the args inc_s, inc_h, or inc_d must be set')
! Subprogram not used    end if
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
! Subprogram not used 
! Subprogram not used    ! increment the time instant
! Subprogram not used    date2 = date1 + t_interval
! Subprogram not used 
! Subprogram not used    ! extract the time components
! Subprogram not used    call ESMF_TimeGet(date2, yy=year, mm=month, dd=day, s=tod2, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used    ymd2 = year*10000 + month*100 + day
! Subprogram not used 
! Subprogram not used end subroutine timemgr_time_inc

!=========================================================================================

subroutine chkrc(rc, mes)
   integer, intent(in)          :: rc   ! return code from time management library
   character(len=*), intent(in) :: mes  ! error message
   if ( rc == ESMF_SUCCESS ) return
   write(iulog,*) mes
   call endrun ('CHKRC')
end subroutine chkrc

end module time_manager
