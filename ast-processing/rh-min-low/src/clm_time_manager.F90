module clm_time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_sys_mod , only: shr_sys_abort
   use spmdMod     , only: masterproc
   use clm_varctl  , only: iulog
   use clm_varcon  , only: isecspday
   use ESMF

   implicit none
   private

   ! Public methods

   public ::&
        get_timemgr_defaults,     &! get startup default values
        set_timemgr_init,         &! setup startup values
        timemgr_init,             &! time manager initialization
        timemgr_restart_io,       &! read/write time manager restart info and restart time manager
        timemgr_restart,          &! restart the time manager using info from timemgr_restart
        timemgr_datediff,         &! calculate difference between two time instants
        advance_timestep,         &! increment timestep number
        get_clock,                &! get the clock from the time-manager
        get_curr_ESMF_Time,       &! get current time in terms of the ESMF_Time
        get_step_size,            &! return step size in seconds
        get_rad_step_size,        &! return radiation step size in seconds
        get_nstep,                &! return timestep number
        get_curr_date,            &! return date components at end of current timestep
        get_prev_date,            &! return date components at beginning of current timestep
        get_start_date,           &! return components of the start date
        get_driver_start_ymd,     &! return year/month/day (as integer in YYYYMMDD format) of driver start date
        get_ref_date,             &! return components of the reference date
        get_perp_date,            &! return components of the perpetual date, and current time of day
        get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
        get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
        get_curr_calday,          &! return calendar day at end of current timestep
        get_calday,               &! return calendar day from input date
        get_calendar,             &! return calendar
        get_days_per_year,        &! return the days per year for current year
        get_curr_yearfrac,        &! return the fractional position in the current year
        get_rest_date,            &! return the date from the restart file
        set_nextsw_cday,          &! set the next radiation calendar day
        is_first_step,            &! return true on first step of initial run
        is_first_restart_step,    &! return true on first step of restart or branch run
        is_beg_curr_day,          &! return true on first timestep in current day
        is_end_curr_day,          &! return true on last timestep in current day
        is_end_curr_month,        &! return true on last timestep in current month
        is_last_step,             &! return true on last timestep
        is_perpetual,             &! return true if perpetual calendar is in use
        is_restart,               &! return true if this is a restart run
        update_rad_dtime           ! track radiation interval via nstep

   ! Public parameter data
   character(len=*), public, parameter :: NO_LEAP_C   = 'NO_LEAP'
   character(len=*), public, parameter :: GREGORIAN_C = 'GREGORIAN'

   ! Private module data

   ! Private data for input

   character(len=ESMF_MAXSTR), save ::&
        calendar   = NO_LEAP_C        ! Calendar to use in date calculations.
   integer,  parameter :: uninit_int = -999999999
   real(r8), parameter :: uninit_r8  = -999999999.0

   ! Input
   integer, save ::&
        dtime          = uninit_int,  &! timestep in seconds
        dtime_rad      = uninit_int,  &! radiation interval in seconds
        nstep_rad_prev = uninit_int    ! radiation interval in seconds

   ! Input from CESM driver
   integer, save ::&
        nelapse       = uninit_int,  &! number of timesteps (or days if negative) to extend a run
        start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
        start_tod     = 0,           &! starting time of day for run in seconds
        stop_ymd      = uninit_int,  &! stopping date for run in yearmmdd format
        stop_tod      = 0,           &! stopping time of day for run in seconds
        ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
        ref_tod       = 0             ! reference time of day for time coordinate in seconds
   type(ESMF_Calendar), target, save   :: tm_cal       ! calendar
   type(ESMF_Clock),    save   :: tm_clock     ! model clock   
   type(ESMF_Time),     save   :: tm_perp_date ! perpetual date

   ! Data required to restart time manager:
   integer, save :: rst_step_sec          = uninit_int ! timestep size seconds
   integer, save :: rst_start_ymd         = uninit_int ! start date
   integer, save :: rst_start_tod         = uninit_int ! start time of day
   integer, save :: rst_ref_ymd           = uninit_int ! reference date
   integer, save :: rst_ref_tod           = uninit_int ! reference time of day
   integer, save :: rst_curr_ymd          = uninit_int ! current date
   integer, save :: rst_curr_tod          = uninit_int ! current time of day

   integer, save :: rst_nstep_rad_prev                 ! nstep of previous radiation call
   integer, save :: perpetual_ymd         = uninit_int ! Perpetual calendar date (YYYYMMDD)
   logical, save :: tm_first_restart_step = .false.    ! true for first step of a restart or branch run
   logical, save :: tm_perp_calendar      = .false.    ! true when using perpetual calendar
   logical, save :: timemgr_set           = .false.    ! true when timemgr initialized
   integer, save :: nestep                = uninit_int ! ending time-step
   !
   ! Next short-wave radiation calendar day
   ! 
   real(r8) :: nextsw_cday = uninit_r8 ! calday from clock of next radiation computation

   ! Private module methods

   private :: timemgr_spmdbcast
   private :: init_calendar
   private :: init_clock
   private :: calc_nestep
   private :: timemgr_print
   private :: TimeGetymd

   !=========================================================================================
contains
  !=========================================================================================

! Subprogram not used   subroutine get_timemgr_defaults( calendar_out,      start_ymd_out,     start_tod_out, ref_ymd_out,        &
! Subprogram not used        ref_tod_out,       stop_ymd_out,      stop_tod_out,  nelapse_out,        &
! Subprogram not used        dtime_out )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! get time manager startup default values
! Subprogram not used     ! 
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), optional, intent(OUT) :: calendar_out       ! Calendar type
! Subprogram not used     integer         , optional, intent(OUT) :: nelapse_out        ! Number of step (or days) to advance
! Subprogram not used     integer         , optional, intent(OUT) :: start_ymd_out      ! Start date       (YYYYMMDD)
! Subprogram not used     integer         , optional, intent(OUT) :: start_tod_out      ! Start time of day (sec)
! Subprogram not used     integer         , optional, intent(OUT) :: ref_ymd_out        ! Reference date   (YYYYMMDD)
! Subprogram not used     integer         , optional, intent(OUT) :: ref_tod_out        ! Reference time of day (sec)
! Subprogram not used     integer         , optional, intent(OUT) :: stop_ymd_out       ! Stop date        (YYYYMMDD)
! Subprogram not used     integer         , optional, intent(OUT) :: stop_tod_out       ! Stop time of day (sec)
! Subprogram not used     integer         , optional, intent(OUT) :: dtime_out          ! Time-step (sec)
! Subprogram not used     !
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_timemgr_defaults'
! Subprogram not used 
! Subprogram not used     if ( timemgr_set ) call shr_sys_abort( sub//":: timemgr_init or timemgr_restart already called" )
! Subprogram not used     if (present(calendar_out)      ) calendar_out       = trim(calendar)
! Subprogram not used     if (present(start_ymd_out)     ) start_ymd_out      = start_ymd
! Subprogram not used     if (present(start_tod_out)     ) start_tod_out      = start_tod
! Subprogram not used     if (present(ref_ymd_out)       ) ref_ymd_out        = ref_ymd
! Subprogram not used     if (present(ref_tod_out)       ) ref_tod_out        = ref_tod
! Subprogram not used     if (present(stop_ymd_out)      ) stop_ymd_out       = stop_ymd
! Subprogram not used     if (present(stop_tod_out)      ) stop_tod_out       = stop_tod
! Subprogram not used     if (present(nelapse_out)       ) nelapse_out        = nelapse
! Subprogram not used     if (present(dtime_out)         ) dtime_out          = dtime
! Subprogram not used 
! Subprogram not used   end subroutine get_timemgr_defaults

  !=========================================================================================

  subroutine set_timemgr_init( calendar_in,      start_ymd_in,     start_tod_in, ref_ymd_in,        &
       ref_tod_in,       stop_ymd_in,      stop_tod_in,  perpetual_run_in,  &
       perpetual_ymd_in, nelapse_in,       dtime_in )

    !---------------------------------------------------------------------------------
    ! set time manager startup values
    ! 
    ! Arguments
    character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
    integer         , optional, intent(IN) :: nelapse_in        ! Number of step (or days) to advance
    integer         , optional, intent(IN) :: start_ymd_in      ! Start date       (YYYYMMDD)
    integer         , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
    integer         , optional, intent(IN) :: ref_ymd_in        ! Reference date   (YYYYMMDD)
    integer         , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
    integer         , optional, intent(IN) :: stop_ymd_in       ! Stop date        (YYYYMMDD)
    integer         , optional, intent(IN) :: stop_tod_in       ! Stop time of day (sec)
    logical         , optional, intent(IN) :: perpetual_run_in  ! If in perpetual mode or not
    integer         , optional, intent(IN) :: perpetual_ymd_in  ! Perpetual date   (YYYYMMDD)
    integer         , optional, intent(IN) :: dtime_in          ! Time-step (sec)
    !
    character(len=*), parameter :: sub = 'clm::set_timemgr_init'

    if ( timemgr_set ) call shr_sys_abort( sub//":: timemgr_init or timemgr_restart already called" )
    if (present(calendar_in)      ) calendar         = trim(calendar_in)
    if (present(start_ymd_in)     ) start_ymd        = start_ymd_in
    if (present(start_tod_in)     ) start_tod        = start_tod_in
    if (present(ref_ymd_in)       ) ref_ymd          = ref_ymd_in
    if (present(ref_tod_in)       ) ref_tod          = ref_tod_in
    if (present(stop_ymd_in)      ) stop_ymd         = stop_ymd_in
    if (present(stop_tod_in)      ) stop_tod         = stop_tod_in
    if (present(perpetual_run_in) )then
       tm_perp_calendar = perpetual_run_in
       if ( tm_perp_calendar ) then
          if ( .not. present(perpetual_ymd_in) .or. perpetual_ymd == uninit_int) &
               call shr_sys_abort( sub//":: perpetual_run set but NOT perpetual_ymd" )
          perpetual_ymd    = perpetual_ymd_in
       end if
    end if
    if (present(nelapse_in)       ) nelapse          = nelapse_in
    if (present(dtime_in)         ) dtime            = dtime_in

  end subroutine set_timemgr_init

  !=========================================================================================

  subroutine timemgr_init( )

    !---------------------------------------------------------------------------------
    ! Initialize the ESMF time manager from the sync clock
    ! 
    ! Arguments
    !
    character(len=*), parameter :: sub = 'clm::timemgr_init'
    integer :: rc                            ! return code
    integer :: yr, mon, day, tod             ! Year, month, day, and second as integers
    type(ESMF_Time) :: start_date            ! start date for run
    type(ESMF_Time) :: stop_date             ! stop date for run
    type(ESMF_Time) :: curr_date             ! temporary date used in logic
    type(ESMF_Time) :: ref_date              ! reference date for time coordinate
    logical :: run_length_specified = .false.
    type(ESMF_Time) :: current               ! current date (from clock)
    type(ESMF_TimeInterval) :: day_step_size ! day step size
    type(ESMF_TimeInterval) :: step_size     ! timestep size
    !---------------------------------------------------------------------------------
    call timemgr_spmdbcast( )

    ! Initalize calendar 

    call init_calendar()

    ! Initalize start date.

    if ( start_ymd == uninit_int ) then
       write(iulog,*)sub,': start_ymd must be specified '
       call shr_sys_abort
    end if
    if ( start_tod == uninit_int ) then
       write(iulog,*)sub,': start_tod must be specified '
       call shr_sys_abort
    end if
    start_date = TimeSetymd( start_ymd, start_tod, "start_date" )

    ! Initialize current date

    curr_date = start_date

    ! Initalize stop date.

    stop_date = TimeSetymd( 99991231, stop_tod, "stop_date" )

    call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

    call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')

    if ( stop_ymd /= uninit_int ) then
       current = TimeSetymd( stop_ymd, stop_tod, "stop_date" )
       if ( current < stop_date ) stop_date = current
       run_length_specified = .true.
    end if
    if ( nelapse /= uninit_int ) then
       if ( nelapse >= 0 ) then
          current = curr_date + step_size*nelapse
       else
          current = curr_date - day_step_size*nelapse
       end if
       if ( current < stop_date ) stop_date = current
       run_length_specified = .true.
    end if
    if ( .not. run_length_specified ) then
       call shr_sys_abort (sub//': Must specify stop_ymd or nelapse')
    end if

    ! Error check 

    if ( stop_date <= start_date ) then
       write(iulog,*)sub, ': stop date must be specified later than start date: '
       call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
       call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
       call shr_sys_abort
    end if
    if ( curr_date >= stop_date ) then
       write(iulog,*)sub, ': stop date must be specified later than current date: '
       call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
       call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
       call shr_sys_abort
    end if

    ! Initalize reference date for time coordinate.

    if ( ref_ymd /= uninit_int ) then
       ref_date = TimeSetymd( ref_ymd, ref_tod, "ref_date" )
    else
       ref_date = start_date
    end if

    ! Initialize clock

    call init_clock( start_date, ref_date, curr_date, stop_date )

    ! Initialize date used for perpetual calendar day calculation.

    if (tm_perp_calendar) then
       tm_perp_date = TimeSetymd( perpetual_ymd, 0, "tm_perp_date" )
    end if

    ! Print configuration summary to log file (stdout).

    if (masterproc) call timemgr_print()

    timemgr_set = .true.

  end subroutine timemgr_init

  !=========================================================================================

  subroutine init_clock( start_date, ref_date, curr_date, stop_date )

    !---------------------------------------------------------------------------------
    ! Purpose: Initialize the clock based on the start_date, ref_date, and curr_date
    ! as well as the settings from the namelist specifying the time to stop
    !
    type(ESMF_Time), intent(in) :: start_date  ! start date for run
    type(ESMF_Time), intent(in) :: ref_date    ! reference date for time coordinate
    type(ESMF_Time), intent(in) :: curr_date   ! current date (equal to start_date)
    type(ESMF_Time), intent(in) :: stop_date   ! stop date for run
    !
    character(len=*), parameter :: sub = 'clm::init_clock'
    type(ESMF_TimeInterval)     :: step_size         ! timestep size
    type(ESMF_Time)             :: current           ! current date (from clock)
    integer                     :: rc                ! return code
    !---------------------------------------------------------------------------------

    call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

    ! Initialize the clock

    tm_clock = ESMF_ClockCreate(name="CLM Time-manager clock", timeStep=step_size, startTime=start_date, &
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
  end subroutine init_clock

  !=========================================================================================

  function TimeSetymd( ymd, tod, desc )
    !---------------------------------------------------------------------------------
    !
    ! Set the time by an integer as YYYYMMDD and integer seconds in the day
    !
    integer, intent(in) :: ymd            ! Year, month, day YYYYMMDD
    integer, intent(in) :: tod            ! Time of day in seconds
    character(len=*), intent(in) :: desc  ! Description of time to set

    type(ESMF_Time) :: TimeSetymd    ! Return value

    character(len=*), parameter :: sub = 'clm::TimeSetymd'
    integer :: yr, mon, day          ! Year, month, day as integers
    integer :: rc                    ! return code
    !---------------------------------------------------------------------------------

    if ( (ymd < 0) .or. (tod < 0) .or. (tod > isecspday) )then
       write(iulog,*) sub//': error yymmdd is a negative number or time-of-day out of bounds', &
            ymd, tod
       call shr_sys_abort
    end if
    yr  = ymd / 10000
    mon = (ymd - yr*10000) / 100
    day =  ymd - yr*10000 - mon*100
    call ESMF_TimeSet( TimeSetymd, yy=yr, mm=mon, dd=day, s=tod, &
         calendar=tm_cal, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeSet: setting '//trim(desc))
  end function TimeSetymd

  !=========================================================================================

! Subprogram not used   integer function TimeGetymd( date, tod )
! Subprogram not used     !
! Subprogram not used     ! Get the date and time of day in ymd from ESMF Time.
! Subprogram not used     !
! Subprogram not used     type(ESMF_Time), intent(inout) :: date ! Input date to convert to ymd
! Subprogram not used     integer, intent(out), optional :: tod  ! Time of day in seconds
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: sub = 'clm::TimeGetymd'
! Subprogram not used     integer :: yr, mon, day
! Subprogram not used     integer :: rc                          ! return code
! Subprogram not used 
! Subprogram not used     call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used     TimeGetymd = yr*10000 + mon*100 + day
! Subprogram not used     if ( present( tod ) )then
! Subprogram not used        call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used        call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used     end if
! Subprogram not used     if ( yr < 0 )then
! Subprogram not used        write(iulog,*) sub//': error year is less than zero', yr
! Subprogram not used        call shr_sys_abort
! Subprogram not used     end if
! Subprogram not used   end function TimeGetymd

  !=========================================================================================

! Subprogram not used   subroutine timemgr_restart_io( ncid, flag )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Read/Write information needed on restart to a netcdf file. 
! Subprogram not used     use ncdio_pio, only: ncd_int
! Subprogram not used     use pio,       only: var_desc_t, file_desc_t
! Subprogram not used     use restUtilMod
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid  ! netcdf id
! Subprogram not used     character(len=*), intent(in) :: flag  ! 'read' or 'write'
! Subprogram not used     !
! Subprogram not used     ! Local variables
! Subprogram not used     character(len=*), parameter :: sub = 'clm::timemgr_restart'
! Subprogram not used     integer :: rc                  ! return code
! Subprogram not used     logical :: readvar             ! determine if variable is on initial file
! Subprogram not used     type(ESMF_Time) :: start_date  ! start date for run
! Subprogram not used     type(ESMF_Time) :: ref_date    ! reference date for run
! Subprogram not used     type(ESMF_Time) :: curr_date   ! date of data in restart file
! Subprogram not used     integer :: rst_caltype         ! calendar type
! Subprogram not used     integer, parameter :: noleap = 1
! Subprogram not used     integer, parameter :: gregorian = 2
! Subprogram not used     character(len=len(calendar)) :: cal
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'write') then
! Subprogram not used        rst_nstep_rad_prev  = nstep_rad_prev
! Subprogram not used     end if
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_nstep_rad_prev', xtype=ncd_int,  &
! Subprogram not used          long_name='previous_radiation_nstep', units='unitless positive integer', &
! Subprogram not used          ifill_value=uninit_int, &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_nstep_rad_prev)
! Subprogram not used     if (flag == 'read') then
! Subprogram not used        nstep_rad_prev = rst_nstep_rad_prev
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag == 'write') then
! Subprogram not used        cal = to_upper(calendar)
! Subprogram not used        if ( trim(cal) == NO_LEAP_C ) then
! Subprogram not used           rst_caltype = noleap
! Subprogram not used        else if ( trim(cal) == GREGORIAN_C ) then
! Subprogram not used           rst_caltype = gregorian
! Subprogram not used        else
! Subprogram not used           call shr_sys_abort(sub//'ERROR: unrecognized calendar specified= '//trim(calendar))
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_type', xtype=ncd_int,  &
! Subprogram not used          long_name='calendar type', units='unitless', flag_meanings=(/ "NO_LEAP_C", "GREGORIAN" /), &
! Subprogram not used          flag_values=(/ noleap, gregorian /), ifill_value=uninit_int, &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_caltype)
! Subprogram not used     if (flag == 'read') then
! Subprogram not used        if ( rst_caltype == noleap ) then
! Subprogram not used           calendar = NO_LEAP_C
! Subprogram not used        else if ( rst_caltype == gregorian ) then
! Subprogram not used           calendar = GREGORIAN_C
! Subprogram not used        else
! Subprogram not used           write(iulog,*)sub,': unrecognized calendar type in restart file: ',rst_caltype
! Subprogram not used           call shr_sys_abort( sub//'ERROR: bad calendar type in restart file')
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag == 'write') then
! Subprogram not used        call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, refTime=ref_date, rc=rc )
! Subprogram not used        call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used        rst_step_sec  = dtime
! Subprogram not used        rst_start_ymd = TimeGetymd( start_date, tod=rst_start_tod )
! Subprogram not used        rst_ref_ymd   = TimeGetymd( ref_date,   tod=rst_ref_tod   )
! Subprogram not used        rst_curr_ymd  = TimeGetymd( curr_date,  tod=rst_curr_tod  )
! Subprogram not used     end if
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_step_sec', xtype=ncd_int, &
! Subprogram not used          long_name='seconds component of timestep size', units='sec',         &
! Subprogram not used          nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_step_sec)
! Subprogram not used     if ((flag == 'read') .and. ( rst_step_sec < 0 .or. rst_step_sec > isecspday )) then
! Subprogram not used        call shr_sys_abort( sub//'ERROR: timemgr_rst_step_sec out of range')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_ymd', xtype=ncd_int, &
! Subprogram not used          long_name='start date', units='YYYYMMDD', ifill_value=uninit_int,     &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_start_ymd)
! Subprogram not used 
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_tod', xtype=ncd_int, &
! Subprogram not used          long_name='start time of day', units='sec',                           &
! Subprogram not used          nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                 &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_start_tod)
! Subprogram not used     if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
! Subprogram not used        call shr_sys_abort( sub//'ERROR: timemgr_rst_strart_tod out of range')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_ymd', xtype=ncd_int,   &
! Subprogram not used          long_name='reference date', units='YYYYMMDD', ifill_value=uninit_int, &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_ref_ymd)
! Subprogram not used 
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_tod', xtype=ncd_int,   &
! Subprogram not used          long_name='reference time of day', units='sec',                       &
! Subprogram not used          nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                 &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_ref_tod)
! Subprogram not used     if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
! Subprogram not used        call shr_sys_abort( sub//'ERROR: timemgr_rst_ref_tod out of range')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_ymd', xtype=ncd_int,  &
! Subprogram not used          long_name='current date', units='YYYYMMDD', ifill_value=uninit_int,   &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_curr_ymd)
! Subprogram not used 
! Subprogram not used     call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_tod', xtype=ncd_int,  &
! Subprogram not used          long_name='current time of day', units='sec',                         &
! Subprogram not used          nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                 &
! Subprogram not used          interpinic_flag='skip', readvar=readvar, data=rst_curr_tod)
! Subprogram not used     if ((flag == 'read') .and. ( rst_curr_tod < 0 .or. rst_curr_tod > isecspday )) then
! Subprogram not used        call shr_sys_abort( sub//'ERROR: timemgr_rst_ref_ymd out of range')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine timemgr_restart_io

  !=========================================================================================

! Subprogram not used   subroutine timemgr_restart( )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Restart the ESMF time manager using the synclock for ending date.
! Subprogram not used     !
! Subprogram not used     character(len=*), parameter :: sub = 'clm::timemgr_restart'
! Subprogram not used     integer :: rc                            ! return code
! Subprogram not used     integer :: yr, mon, day, tod             ! Year, month, day, and second as integers
! Subprogram not used     type(ESMF_Time) :: start_date            ! start date for run
! Subprogram not used     type(ESMF_Time) :: ref_date              ! reference date for run
! Subprogram not used     type(ESMF_Time) :: curr_date             ! date of data in restart file
! Subprogram not used     type(ESMF_Time) :: stop_date             ! stop date for run
! Subprogram not used     type(ESMF_Time) :: current               ! current date (from clock)
! Subprogram not used     type(ESMF_TimeInterval) :: day_step_size ! day step size
! Subprogram not used     type(ESMF_TimeInterval) :: step_size     ! timestep size
! Subprogram not used     logical :: run_length_specified = .false.
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     call timemgr_spmdbcast( )
! Subprogram not used 
! Subprogram not used     ! Initialize calendar from restart info
! Subprogram not used 
! Subprogram not used     call init_calendar()
! Subprogram not used 
! Subprogram not used     ! Initialize the timestep from restart info
! Subprogram not used 
! Subprogram not used     dtime = rst_step_sec
! Subprogram not used 
! Subprogram not used     ! Initialize start date from restart info
! Subprogram not used 
! Subprogram not used     start_date = TimeSetymd( rst_start_ymd, rst_start_tod, "start_date" )
! Subprogram not used 
! Subprogram not used     ! Initialize current date from restart info
! Subprogram not used 
! Subprogram not used     curr_date = TimeSetymd( rst_curr_ymd, rst_curr_tod, "curr_date" )
! Subprogram not used 
! Subprogram not used     ! Initialize stop date from sync clock or namelist input
! Subprogram not used 
! Subprogram not used     stop_date = TimeSetymd( 99991231, stop_tod, "stop_date" )
! Subprogram not used 
! Subprogram not used     call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')
! Subprogram not used 
! Subprogram not used     call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')
! Subprogram not used 
! Subprogram not used     if    ( stop_ymd /= uninit_int ) then
! Subprogram not used        current = TimeSetymd( stop_ymd, stop_tod, "stop_date" )
! Subprogram not used        if ( current < stop_date ) stop_date = current
! Subprogram not used        run_length_specified = .true.
! Subprogram not used     else if ( nelapse /= uninit_int ) then
! Subprogram not used        if ( nelapse >= 0 ) then
! Subprogram not used           current = curr_date + step_size*nelapse
! Subprogram not used        else
! Subprogram not used           current = curr_date - day_step_size*nelapse
! Subprogram not used        end if
! Subprogram not used        if ( current < stop_date ) stop_date = current
! Subprogram not used        run_length_specified = .true.
! Subprogram not used     end if
! Subprogram not used     if ( .not. run_length_specified ) then
! Subprogram not used        call shr_sys_abort (sub//': Must specify stop_ymd or nelapse')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Error check
! Subprogram not used 
! Subprogram not used     if ( stop_date <= start_date ) then
! Subprogram not used        write(iulog,*)sub, ': stop date must be specified later than start date: '
! Subprogram not used        call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used        write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used        call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used        write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used        call shr_sys_abort
! Subprogram not used     end if
! Subprogram not used     if ( curr_date >= stop_date ) then
! Subprogram not used        write(iulog,*)sub, ': stop date must be specified later than current date: '
! Subprogram not used        call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used        write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used        call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used        write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used        call shr_sys_abort
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Initialize nstep_rad_prev from restart info
! Subprogram not used 
! Subprogram not used     nstep_rad_prev = rst_nstep_rad_prev
! Subprogram not used 
! Subprogram not used     ! Initialize ref date from restart info
! Subprogram not used 
! Subprogram not used     ref_date = TimeSetymd( rst_ref_ymd, rst_ref_tod, "ref_date" )
! Subprogram not used 
! Subprogram not used     ! Initialize clock 
! Subprogram not used 
! Subprogram not used     call init_clock( start_date, ref_date, curr_date, stop_date )
! Subprogram not used 
! Subprogram not used     ! Advance the timestep.  
! Subprogram not used     ! Data from the restart file corresponds to the last timestep of the previous run.
! Subprogram not used 
! Subprogram not used     call advance_timestep()
! Subprogram not used 
! Subprogram not used     ! Set flag that this is the first timestep of the restart run.
! Subprogram not used 
! Subprogram not used     tm_first_restart_step = .true.
! Subprogram not used 
! Subprogram not used     ! Calculate ending time step
! Subprogram not used 
! Subprogram not used     call calc_nestep( )
! Subprogram not used 
! Subprogram not used     ! Print configuration summary to log file (stdout).
! Subprogram not used 
! Subprogram not used     if (masterproc) call timemgr_print()
! Subprogram not used 
! Subprogram not used     timemgr_set = .true.
! Subprogram not used 
! Subprogram not used   end subroutine timemgr_restart

  !=========================================================================================

! Subprogram not used   subroutine calc_nestep()
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Calculate ending timestep number
! Subprogram not used     ! Calculation of ending timestep number (nestep) assumes a constant stepsize.
! Subprogram not used     !
! Subprogram not used     character(len=*), parameter :: sub = 'clm::calc_nestep'
! Subprogram not used     integer :: ntspday               ! Number of time-steps per day
! Subprogram not used     type(ESMF_TimeInterval) :: diff  !
! Subprogram not used     type(ESMF_Time) :: start_date    ! start date for run
! Subprogram not used     type(ESMF_Time) :: stop_date     ! stop date for run
! Subprogram not used     integer :: ndays, nsecs          ! Number of days, seconds to ending time
! Subprogram not used     integer :: rc                    ! return code
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call ESMF_ClockGet( tm_clock, stopTime=stop_date, startTime=start_date, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used     ntspday = isecspday/dtime
! Subprogram not used     diff = stop_date - start_date
! Subprogram not used     call ESMF_TimeIntervalGet( diff, d=ndays, s=nsecs, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet calculating nestep')
! Subprogram not used     nestep = ntspday*ndays + nsecs/dtime
! Subprogram not used     if ( mod(nsecs,dtime) /= 0 ) nestep = nestep + 1
! Subprogram not used   end subroutine calc_nestep

  !=========================================================================================

  subroutine init_calendar( )

    !---------------------------------------------------------------------------------
    ! Initialize calendar
    !
    ! Local variables
    !
    character(len=*), parameter :: sub = 'clm::init_calendar'
    type(ESMF_CalKind_Flag) :: cal_type        ! calendar type
    character(len=len(calendar)) :: caltmp
    integer :: rc                              ! return code
    !---------------------------------------------------------------------------------

    caltmp = to_upper(calendar)
    if ( trim(caltmp) == NO_LEAP_C ) then
       cal_type = ESMF_CALKIND_NOLEAP
    else if ( trim(caltmp) == GREGORIAN_C ) then
       cal_type = ESMF_CALKIND_GREGORIAN
    else
       write(iulog,*)sub,': unrecognized calendar specified: ',calendar
       call shr_sys_abort
    end if
    tm_cal = ESMF_CalendarCreate( name=caltmp, calkindflag=cal_type, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_CalendarSet')
  end subroutine init_calendar

  !=========================================================================================

  subroutine timemgr_print()

    !---------------------------------------------------------------------------------
    character(len=*), parameter :: sub = 'clm::timemgr_print'
    integer :: rc
    integer :: yr, mon, day
    integer :: &                   ! Data required to restart time manager:
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
    !---------------------------------------------------------------------------------

    call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, &
         refTime=ref_date, stopTime=stop_date, timeStep=step, &
         advanceCount=step_no, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')
    nstep = step_no

    write(iulog,*)' ******** CLM Time Manager Configuration ********'

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
    write(iulog,*)'  Ending step number:       ', nestep
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

    character(len=*), parameter :: sub = 'clm::advance_timestep'
    integer :: rc

    call ESMF_ClockAdvance( tm_clock, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockAdvance')

    tm_first_restart_step = .false.

  end subroutine advance_timestep

  !=========================================================================================

! Subprogram not used   subroutine get_clock( clock )
! Subprogram not used 
! Subprogram not used     ! Return the ESMF clock
! Subprogram not used 
! Subprogram not used     type(ESMF_Clock), intent(inout) :: clock
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_clock'
! Subprogram not used     type(ESMF_TimeInterval) :: step_size
! Subprogram not used     type(ESMF_Time) :: start_date, stop_date, ref_date
! Subprogram not used     integer :: rc
! Subprogram not used 
! Subprogram not used     call ESMF_ClockGet( tm_clock, timeStep=step_size, startTime=start_date, &
! Subprogram not used          stoptime=stop_date, reftime=ref_date, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used     call ESMF_ClockSet(clock, timeStep=step_size, startTime=start_date, &
! Subprogram not used          stoptime=stop_date, reftime=ref_date, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_ClockSet')
! Subprogram not used 
! Subprogram not used   end subroutine get_clock

  !=========================================================================================

! Subprogram not used   function get_curr_ESMF_Time( )
! Subprogram not used 
! Subprogram not used     ! Return the current time as ESMF_Time
! Subprogram not used 
! Subprogram not used     type(ESMF_Time) :: get_curr_ESMF_Time
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_curr_ESMF_Time'
! Subprogram not used     integer :: rc
! Subprogram not used 
! Subprogram not used     call ESMF_ClockGet( tm_clock, currTime=get_curr_ESMF_Time, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used   end function get_curr_ESMF_Time

  !=========================================================================================

  integer function get_step_size()

    ! Return the step size in seconds.

    character(len=*), parameter :: sub = 'clm::get_step_size'
    type(ESMF_TimeInterval) :: step_size       ! timestep size
    integer :: rc

    call ESMF_ClockGet(tm_clock, timeStep=step_size, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_TimeIntervalGet(step_size, s=get_step_size, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')

  end function get_step_size

  !=========================================================================================

  subroutine update_rad_dtime(doalb)
    !---------------------------------------------------------------------------------
    ! called only on doalb timesteps to save off radiation nsteps
    ! 
    ! Local Arguments
    logical,intent(in) ::  doalb
    integer :: dtime,nstep

    if (doalb) then 

       dtime=get_step_size()
       nstep = get_nstep()

       if (nstep_rad_prev == uninit_int ) then
          dtime_rad = dtime
          nstep_rad_prev = nstep
       else
          dtime_rad = (nstep - nstep_rad_prev) * dtime
          nstep_rad_prev = nstep
       endif
    end if
  end subroutine update_rad_dtime

  !=========================================================================================

! Subprogram not used   integer function get_rad_step_size()
! Subprogram not used 
! Subprogram not used     if (nstep_rad_prev == uninit_int ) then
! Subprogram not used        get_rad_step_size=get_step_size()
! Subprogram not used     else
! Subprogram not used        get_rad_step_size=dtime_rad
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function get_rad_step_size

  !=========================================================================================

  integer function get_nstep()

    ! Return the timestep number.

    character(len=*), parameter :: sub = 'clm::get_nstep'
    integer :: rc
    integer(ESMF_KIND_I8) :: step_no

    call ESMF_ClockGet(tm_clock, advanceCount=step_no, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    get_nstep = step_no

  end function get_nstep

  !=========================================================================================

  subroutine get_curr_date(yr, mon, day, tod, offset)

    !-----------------------------------------------------------------------------------------
    ! Return date components valid at end of current timestep with an optional
    ! offset (positive or negative) in seconds.

    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    character(len=*), parameter :: sub = 'clm::get_curr_date'
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

! Subprogram not used   subroutine get_perp_date(yr, mon, day, tod, offset)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------------------------
! Subprogram not used     ! Return time of day valid at end of current timestep and the components
! Subprogram not used     ! of the perpetual date (with an optional offset (positive or negative) in seconds.
! Subprogram not used 
! Subprogram not used     integer, intent(out) ::&
! Subprogram not used          yr,    &! year
! Subprogram not used          mon,   &! month
! Subprogram not used          day,   &! day of month
! Subprogram not used          tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used     integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
! Subprogram not used     ! Positive for future times, negative 
! Subprogram not used     ! for previous times.
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_perp_date'
! Subprogram not used     integer :: rc
! Subprogram not used     type(ESMF_Time) :: date
! Subprogram not used     type(ESMF_TimeInterval) :: DelTime
! Subprogram not used     !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
! Subprogram not used     ! Get time of day add it to perpetual date
! Subprogram not used     ! Get year, month, day so that seconds are time-of-day rather than since start time
! Subprogram not used     call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used     call ESMF_TimeIntervalSet(DelTime, s=tod, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
! Subprogram not used     date = tm_perp_date + DelTime
! Subprogram not used     if ( present(offset) )then
! Subprogram not used        call ESMF_TimeIntervalSet(DelTime, s=offset, rc=rc)
! Subprogram not used        call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
! Subprogram not used        date = date + DelTime
! Subprogram not used     end if
! Subprogram not used     ! Get time of day from the result
! Subprogram not used     ! Get year, month, day so that seconds are time-of-day rather than since start time
! Subprogram not used     call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used 
! Subprogram not used     ! Get the date from the fixed perpetual date (in case it overflows to next day)
! Subprogram not used     call ESMF_TimeGet(tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used   end subroutine get_perp_date

  !=========================================================================================

  subroutine get_prev_date(yr, mon, day, tod)

    ! Return date components valid at beginning of current timestep.

    ! Arguments
    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_prev_date'
    integer :: rc
    type(ESMF_Time) :: date
    !-----------------------------------------------------------------------------------------

    call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end subroutine get_prev_date

  !=========================================================================================

! Subprogram not used   subroutine get_start_date(yr, mon, day, tod)
! Subprogram not used 
! Subprogram not used     ! Return date components valid at beginning of initial run.
! Subprogram not used 
! Subprogram not used     ! Arguments
! Subprogram not used     integer, intent(out) ::&
! Subprogram not used          yr,    &! year
! Subprogram not used          mon,   &! month
! Subprogram not used          day,   &! day of month
! Subprogram not used          tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_start_date'
! Subprogram not used     integer :: rc
! Subprogram not used     type(ESMF_Time) :: date
! Subprogram not used     !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call ESMF_ClockGet(tm_clock, startTime=date, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used     call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used   end subroutine get_start_date

  !=========================================================================================

! Subprogram not used   integer function get_driver_start_ymd( tod )
! Subprogram not used 
! Subprogram not used     ! Return date of start of simulation from driver (i.e. NOT from restart file)
! Subprogram not used     ! Note: get_start_date gets you the date from the beginning of the simulation
! Subprogram not used     !       on the restart file.
! Subprogram not used 
! Subprogram not used     ! Arguments
! Subprogram not used     integer, optional, intent(out) ::&
! Subprogram not used          tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_driver_start_ymd'
! Subprogram not used     !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if ( start_ymd == uninit_int )then
! Subprogram not used        call shr_sys_abort( sub//': error driver start date is NOT set yet' )
! Subprogram not used     end if
! Subprogram not used     if ( start_ymd < 101 .or. start_ymd > 99991231 )then
! Subprogram not used        call shr_sys_abort( sub//': error driver start date is invalid' )
! Subprogram not used     end if
! Subprogram not used     if ( present(tod) )then
! Subprogram not used        tod = start_tod
! Subprogram not used        if ( (tod < 0) .or. (tod > isecspday) )then
! Subprogram not used           call shr_sys_abort( sub//': error driver start tod is invalid' )
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     get_driver_start_ymd = start_ymd
! Subprogram not used 
! Subprogram not used   end function get_driver_start_ymd

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
    character(len=*), parameter :: sub = 'clm::get_ref_date'
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
    character(len=*), parameter :: sub = 'clm::get_curr_time'
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

  subroutine get_prev_time(days, seconds)

    ! Return time components valid at beg of current timestep.
    ! prev time is the time interval between the prev date and the reference date.

    ! Arguments
    integer, intent(out) ::&
         days,   &! number of whole days in time interval
         seconds  ! remaining seconds in time interval

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_prev_time'
    integer :: rc
    type(ESMF_Time) :: date, ref_date
    type(ESMF_TimeInterval) :: diff
    !-----------------------------------------------------------------------------------------

    call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet for prevTime')
    call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet for refTime')
    diff = date - ref_date
    call ESMF_TimeIntervalGet( diff, d=days, s=seconds, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeintervalGet')

  end subroutine get_prev_time

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
    character(len=*), parameter :: sub = 'clm::get_curr_calday'
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
       call ESMF_TimeGet(date, yy=year, mm=month, dd=day, s=tod, rc=rc)
       call chkrc(rc, sub//': error return from ESMF_TimeGet')
       call ESMF_TimeIntervalSet( diurnal, s=tod, rc=rc )
       call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
       date = tm_perp_date + diurnal
    end if

    call ESMF_TimeGet( date, dayOfYear_r8=get_curr_calday, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    !----------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!! WARNING HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!
    !!!! The following hack fakes day 366 by reusing day 365. This is just because the  !!!!!!
    !!!! current shr_orb_decl calculation can't handle days > 366.                      !!!!!!
    !!!!       Dani Bundy-Coleman and Erik Kluzek Aug/2008                              !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( (get_curr_calday > 366.0) .and. (get_curr_calday <= 367.0) .and. &
         (trim(calendar) == GREGORIAN_C) )then
       get_curr_calday = get_curr_calday - 1.0_r8
    end if
    !!!!!!!!!!!!!! END HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!!!!!
    !----------------------------------------------------------------------------------------!
    if ( (get_curr_calday < 1.0) .or. (get_curr_calday > 366.0) )then
       write(iulog,*) sub, ' = ', get_curr_calday
       if ( present(offset) ) write(iulog,*) 'offset = ', offset
       call shr_sys_abort( sub//': error get_curr_calday out of bounds' )
    end if

  end function get_curr_calday

  !=========================================================================================

! Subprogram not used   function get_calday(ymd, tod)
! Subprogram not used 
! Subprogram not used     ! Return calendar day corresponding to specified time instant.
! Subprogram not used     ! Calendar day 1.0 = 0Z on Jan 1.
! Subprogram not used 
! Subprogram not used     ! Arguments
! Subprogram not used     integer, intent(in) :: &
! Subprogram not used          ymd,   &! date in yearmmdd format
! Subprogram not used          tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used     ! Return value
! Subprogram not used     real(r8) :: get_calday
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_calday'
! Subprogram not used     integer :: rc                 ! return code
! Subprogram not used     type(ESMF_Time) :: date
! Subprogram not used     !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     date = TimeSetymd( ymd, tod, "get_calday" )
! Subprogram not used     call ESMF_TimeGet( date, dayOfYear_r8=get_calday, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used     !----------------------------------------------------------------------------------------!
! Subprogram not used !!!!!!!!!!!!!! WARNING HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!
! Subprogram not used !!!! The following hack fakes day 366 by reusing day 365. This is just because the  !!!!!!
! Subprogram not used !!!! current shr_orb_decl calculation can't handle days > 366.                      !!!!!!
! Subprogram not used !!!!       Dani Bundy-Coleman and Erik Kluzek Aug/2008                              !!!!!!
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used     if ( (get_calday > 366.0) .and. (get_calday <= 367.0) .and. &
! Subprogram not used          (trim(calendar) == GREGORIAN_C) )then
! Subprogram not used        get_calday = get_calday - 1.0_r8
! Subprogram not used     end if
! Subprogram not used !!!!!!!!!!!!!! END HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used     !----------------------------------------------------------------------------------------!
! Subprogram not used     if ( (get_calday < 1.0) .or. (get_calday > 366.0) )then
! Subprogram not used        write(iulog,*) sub, ' = ', get_calday
! Subprogram not used        call shr_sys_abort( sub//': error calday out of range' )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function get_calday

  !=========================================================================================

  function get_calendar()

    ! Return calendar

    ! Return value
    character(len=ESMF_MAXSTR) :: get_calendar

    get_calendar = calendar

  end function get_calendar

  !=========================================================================================

! Subprogram not used   integer function get_days_per_year( offset )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Get the number of days per year for currrent year
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
! Subprogram not used     ! Positive for future times, negative 
! Subprogram not used     ! for previous times.
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_days_per_year'
! Subprogram not used     integer         :: yr, mon, day, tod ! current date year, month, day and time-of-day
! Subprogram not used     type(ESMF_Time) :: eDate             ! ESMF date
! Subprogram not used     integer         :: rc                ! ESMF return code
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if ( present(offset) )then
! Subprogram not used        call get_curr_date(yr, mon, day, tod, offset )
! Subprogram not used     else
! Subprogram not used        call get_curr_date(yr, mon, day, tod )
! Subprogram not used     end if
! Subprogram not used     eDate = TimeSetymd( ymd=yr*10000+1231, tod=0, desc="end of year" )
! Subprogram not used     call ESMF_TimeGet( eDate, dayOfYear=get_days_per_year, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used   end function get_days_per_year

  !=========================================================================================

! Subprogram not used   function get_curr_yearfrac( offset )
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Get the fractional position in the current year. This is 0 at midnight on Jan 1,
! Subprogram not used     ! and 1 at the end of Dec 31.
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     real(r8) :: get_curr_yearfrac  ! function result
! Subprogram not used     
! Subprogram not used     integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
! Subprogram not used     ! Positive for future times, negative 
! Subprogram not used     ! for previous times.
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: sub = 'clm::get_curr_yearfrac'
! Subprogram not used     real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
! Subprogram not used     real(r8) :: days_per_year      ! days per year
! Subprogram not used 
! Subprogram not used     cday          = get_curr_calday(offset=offset)
! Subprogram not used     days_per_year = get_days_per_year()
! Subprogram not used 
! Subprogram not used     get_curr_yearfrac = (cday - 1._r8)/days_per_year
! Subprogram not used 
! Subprogram not used   end function get_curr_yearfrac

  !=========================================================================================

! Subprogram not used   subroutine get_rest_date(ncid, yr)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Get the date from the restart file.
! Subprogram not used     !
! Subprogram not used     ! Currently just returns the year (because the month & day are harder to extract, and
! Subprogram not used     ! currently aren't needed).
! Subprogram not used     use pio,       only: file_desc_t
! Subprogram not used     use ncdio_pio, only: ncd_io
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     type(file_desc_t) , intent(inout) :: ncid ! netcdf id for the restart file
! Subprogram not used     integer           , intent(out)   :: yr   ! year from restart file
! Subprogram not used 
! Subprogram not used     integer :: ymd     ! yyyymmdd from the restart file
! Subprogram not used     logical :: readvar ! whether the variable was read from the file
! Subprogram not used     
! Subprogram not used     integer, parameter :: year_mask = 10000  ! divide by this to get year from ymd
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: subname = 'get_rest_date'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     
! Subprogram not used     ! Get the date (yyyymmdd) from restart file.
! Subprogram not used     ! Note that we cannot simply use the rst_curr_ymd module variable, because that isn't
! Subprogram not used     ! set under some circumstances
! Subprogram not used     call ncd_io(varname='timemgr_rst_curr_ymd', data=ymd, &
! Subprogram not used          ncid=ncid, flag='read', readvar=readvar)
! Subprogram not used     if (.not. readvar) then
! Subprogram not used        call shr_sys_abort(subname//' ERROR: timemgr_rst_curr_ymd not found on restart file')
! Subprogram not used     end if
! Subprogram not used     
! Subprogram not used     ! Extract the year
! Subprogram not used     yr = ymd / year_mask
! Subprogram not used   end subroutine get_rest_date

  !=========================================================================================

  subroutine set_nextsw_cday( nextsw_cday_in )

    ! Set the next radiation calendar day, so that radiation step can be calculated
    !
    ! Arguments
    real(r8), intent(IN) :: nextsw_cday_in ! input calday of next radiation computation

    character(len=*), parameter :: sub = 'clm::set_nextsw_cday'

    nextsw_cday = nextsw_cday_in

  end subroutine set_nextsw_cday

  !=========================================================================================
 
! Subprogram not used   function is_beg_curr_day()
! Subprogram not used  
! Subprogram not used      ! Return true if current timestep is first timestep in current day.
! Subprogram not used      
! Subprogram not used      ! Return value
! Subprogram not used      logical :: is_beg_curr_day
! Subprogram not used   
! Subprogram not used      ! Local variables
! Subprogram not used      integer ::&
! Subprogram not used         yr,    &! year
! Subprogram not used         mon,   &! month
! Subprogram not used         day,   &! day of month
! Subprogram not used         tod     ! time of day (seconds past 0Z)
! Subprogram not used  
! Subprogram not used      call get_curr_date(yr, mon, day, tod)
! Subprogram not used      is_beg_curr_day = ( tod == dtime )
! Subprogram not used  
! Subprogram not used   end function is_beg_curr_day

  !=========================================================================================

  function is_end_curr_day()

    !---------------------------------------------------------------------------------
    ! Return true if current timestep is last timestep in current day.

    ! Return value
    logical :: is_end_curr_day

    ! Local variables
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)
    !---------------------------------------------------------------------------------

    call get_curr_date(yr, mon, day, tod)
    is_end_curr_day = (tod == 0)

  end function is_end_curr_day

  !=========================================================================================

! Subprogram not used   logical function is_end_curr_month()
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Return true if current timestep is last timestep in current month.
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     integer ::&
! Subprogram not used          yr,    &! year
! Subprogram not used          mon,   &! month
! Subprogram not used          day,   &! day of month
! Subprogram not used          tod     ! time of day (seconds past 0Z)
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call get_curr_date(yr, mon, day, tod)
! Subprogram not used     is_end_curr_month = (day == 1  .and.  tod == 0)
! Subprogram not used 
! Subprogram not used   end function is_end_curr_month

  !=========================================================================================

  logical function is_first_step()

    !---------------------------------------------------------------------------------
    ! Return true on first step of initial run only.

    ! Local variables
    character(len=*), parameter :: sub = 'clm::is_first_step'
    integer :: rc
    integer :: nstep
    integer(ESMF_KIND_I8) :: step_no
    !---------------------------------------------------------------------------------

    call ESMF_ClockGet( tm_clock, advanceCount=step_no, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')
    nstep = step_no
    is_first_step = (nstep == 0)

  end function is_first_step
  !=========================================================================================

! Subprogram not used   logical function is_first_restart_step()
! Subprogram not used 
! Subprogram not used     ! Return true on first step of restart run only.
! Subprogram not used 
! Subprogram not used     is_first_restart_step = tm_first_restart_step
! Subprogram not used 
! Subprogram not used   end function is_first_restart_step

  !=========================================================================================

! Subprogram not used   logical function is_last_step()
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used     ! Return true on last timestep.
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     character(len=*), parameter :: sub = 'clm::is_last_step'
! Subprogram not used     type(ESMF_Time) :: stop_date
! Subprogram not used     type(ESMF_Time) :: curr_date
! Subprogram not used     type(ESMF_TimeInterval) :: time_step
! Subprogram not used     integer :: rc
! Subprogram not used     !---------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call ESMF_ClockGet( tm_clock, stopTime=stop_date, &
! Subprogram not used          currTime=curr_date, TimeStep=time_step, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used     if ( curr_date+time_step > stop_date ) then
! Subprogram not used        is_last_step = .true.
! Subprogram not used     else
! Subprogram not used        is_last_step = .false.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function is_last_step

  !=========================================================================================

  logical function is_perpetual()

    ! Return true on last timestep.

    is_perpetual = tm_perp_calendar

  end function is_perpetual

  !=========================================================================================

! Subprogram not used   subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)
! Subprogram not used 
! Subprogram not used     ! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.
! Subprogram not used     ! Arguments
! Subprogram not used     integer, intent(in) ::&
! Subprogram not used          ymd1,    &! date1 in yyyymmdd format
! Subprogram not used          tod1,    &! time of day relative to date1 (seconds past 0Z)
! Subprogram not used          ymd2,    &! date2 in yyyymmdd format
! Subprogram not used          tod2      ! time of day relative to date2 (seconds past 0Z)
! Subprogram not used 
! Subprogram not used     real(r8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     character(len=*), parameter :: sub = 'clm::timemgr_datediff'
! Subprogram not used     integer :: rc   ! return code
! Subprogram not used 
! Subprogram not used     type(ESMF_Time) :: date1
! Subprogram not used     type(ESMF_Time) :: date2
! Subprogram not used     type(ESMF_TimeInterval) :: diff
! Subprogram not used     !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     date1 = TimeSetymd( ymd1, tod1, "date1" )
! Subprogram not used     date2 = TimeSetymd( ymd2, tod2, "date2" )
! Subprogram not used     diff = date2 - date1
! Subprogram not used     call ESMF_TimeIntervalGet( diff, d_r8=days, rc=rc )
! Subprogram not used     call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')
! Subprogram not used     days = days + 1.0_r8
! Subprogram not used 
! Subprogram not used   end subroutine timemgr_datediff

  !=========================================================================================

  subroutine chkrc(rc, mes)
    integer, intent(in)          :: rc   ! return code from time management library
    character(len=*), intent(in) :: mes  ! error message
    if ( rc == ESMF_SUCCESS ) return
    write(iulog,*) mes
    call shr_sys_abort ('CHKRC')
  end subroutine chkrc

  !=========================================================================================

  function to_upper(str)

    !---------------------------------------------------------------------------------
    ! Convert character string to upper case. Use achar and iachar intrinsics
    ! to ensure use of ascii collating sequence.
    !
    ! !INPUT PARAMETERS:
    character(len=*), intent(in) :: str ! String to convert to upper case
    ! !RETURN VALUE:
    character(len=len(str))      :: to_upper
    ! !LOCAL VARIABLES:
    integer :: i                ! Index
    integer :: aseq             ! ascii collating sequence
    character(len=1) :: ctmp    ! Character temporary
    !---------------------------------------------------------------------------------

    do i = 1, len(str)
       ctmp = str(i:i)
       aseq = iachar(ctmp)
       if ( aseq >= 97  .and.  aseq <= 122 ) ctmp = achar(aseq - 32)
       to_upper(i:i) = ctmp
    end do

  end function to_upper

  !=========================================================================================

  logical function is_restart( )
    ! Determine if restart run
    use clm_varctl, only : nsrest, nsrContinue
    if (nsrest == nsrContinue) then
       is_restart = .true.
    else
       is_restart = .false.
    end if
  end function is_restart

  !=========================================================================================

  subroutine timemgr_spmdbcast( )

    use spmdMod, only : mpicom, MPI_INTEGER

    integer :: ier

    call mpi_bcast (dtime    , 1, MPI_INTEGER  , 0, mpicom, ier)

  end subroutine timemgr_spmdbcast

end module clm_time_manager
