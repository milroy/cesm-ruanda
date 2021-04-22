module RtmTimeManager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_sys_mod , only: shr_sys_abort
   use RtmSpmd     , only: masterproc, iam, mpicom_rof, MPI_INTEGER, MPI_CHARACTER
   use RtmVar      , only: isecspday, iulog, nsrest, nsrContinue
   use RtmIO
   use ESMF

   implicit none
   private

! Public methods

   public ::&
      timemgr_setup,            &! setup startup values
      timemgr_init,             &! time manager initialization
      timemgr_restart,          &! read/write time manager restart info and restart time manager
      timemgr_finalize,         &! calls ESMF_ClockDestroy to clean up ESMF clock memory
      advance_timestep,         &! increment timestep number
      get_clock,                &! get the clock from the time-manager
      get_step_size,            &! return step size in seconds
      get_nstep,                &! return timestep number
      get_curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
      get_start_date,           &! return components of the start date
      get_ref_date,             &! return components of the reference date
      get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
      get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
      get_calendar,             &! return calendar
      is_first_step,            &! return true on first step of initial run
      is_first_restart_step,    &! return true on first step of restart or branch run
      is_end_curr_day,          &! return true on last timestep in current day
      is_end_curr_month,        &! return true on last timestep in current month
      is_last_step,             &! return true on last timestep
      is_restart                 ! return true if this is a restart run

! Private module data

! Private data for input

   character(len=ESMF_MAXSTR), save :: calendar   = 'NO_LEAP'     ! Calendar to use in date calculations ('NO_LEAP' or 'GREGORIAN')
   integer,  parameter :: uninit_int = -999999999
   real(r8), parameter :: uninit_r8  = -999999999.0

! Input
   integer, save ::&
      dtime          = uninit_int   ! timestep in seconds

! Input from CESM driver
   integer, save ::&
      nelapse       = uninit_int,  &! number of timesteps (or days if negative) to extend a run
      start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
      start_tod     = 0,           &! starting time of day for run in seconds
      stop_ymd      = uninit_int,  &! stopping date for run in yearmmdd format
      stop_tod      = 0,           &! stopping time of day for run in seconds
      ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
      ref_tod       = 0             ! reference time of day for time coordinate in seconds
   type(ESMF_Calendar), target, save  :: &
        tm_cal       ! calendar
   type(ESMF_Clock), save :: &
        tm_clock     ! model clock   
   integer, save ::&                ! Data required to restart time manager:
      rst_nstep     = uninit_int,  &! current step number
      rst_step_days = uninit_int,  &! days component of timestep size
      rst_step_sec  = uninit_int,  &! timestep size seconds
      rst_start_ymd = uninit_int,  &! start date
      rst_start_tod = uninit_int,  &! start time of day
      rst_ref_ymd   = uninit_int,  &! reference date
      rst_ref_tod   = uninit_int,  &! reference time of day
      rst_curr_ymd  = uninit_int,  &! current date
      rst_curr_tod  = uninit_int    ! current time of day
   character(len=ESMF_MAXSTR), save :: &
   rst_calendar    ! Calendar

   logical, save :: tm_first_restart_step = .false.    ! true for first step of a restart or branch run
   integer, save :: cal_type              = uninit_int ! calendar type
   logical, save :: timemgr_set           = .false.    ! true when timemgr initialized

! Private module methods
   private :: timemgr_spmdbcast
   private :: init_calendar
   private :: init_clock
   private :: timemgr_print
   private :: TimeGetymd

contains

!=========================================================================================

subroutine timemgr_setup( calendar_in, start_ymd_in, start_tod_in, ref_ymd_in, &
                          ref_tod_in,  stop_ymd_in,  stop_tod_in,  nelapse_in)

  ! set time manager startup values
  character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
  integer         , optional, intent(IN) :: nelapse_in        ! Number of step (or days) to advance
  integer         , optional, intent(IN) :: start_ymd_in      ! Start date       (YYYYMMDD)
  integer         , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
  integer         , optional, intent(IN) :: ref_ymd_in        ! Reference date   (YYYYMMDD)
  integer         , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
  integer         , optional, intent(IN) :: stop_ymd_in       ! Stop date        (YYYYMMDD)
  integer         , optional, intent(IN) :: stop_tod_in       ! Stop time of day (sec)
  character(len=*), parameter :: sub = 'rtm::set_timemgr_init'

  ! timemgr_set is called in timemgr_init and timemgr_restart
  if ( timemgr_set ) then
     call shr_sys_abort( sub//":: timemgr_init or timemgr_restart already called" )
  end if
  if (present(calendar_in) ) calendar  = trim(calendar_in)
  if (present(start_ymd_in)) start_ymd = start_ymd_in
  if (present(start_tod_in)) start_tod = start_tod_in
  if (present(ref_ymd_in)  ) ref_ymd   = ref_ymd_in
  if (present(ref_tod_in)  ) ref_tod   = ref_tod_in
  if (present(stop_ymd_in) ) stop_ymd  = stop_ymd_in
  if (present(stop_tod_in) ) stop_tod  = stop_tod_in
  if (present(nelapse_in)  ) nelapse   = nelapse_in

end subroutine timemgr_setup

!=========================================================================================

subroutine timemgr_init( dtime_in )

  ! Initialize the ESMF time manager from the sync clock
  !
  integer, intent(in) :: dtime_in         ! Time-step (sec)
  !
  integer         :: rc                    ! return code
  integer         :: yr, mon, day, tod     ! Year, month, day, and second as integers
  type(ESMF_Time) :: start_date            ! start date for run
  type(ESMF_Time) :: stop_date             ! stop date for run
  type(ESMF_Time) :: curr_date             ! temporary date used in logic
  type(ESMF_Time) :: ref_date              ! reference date for time coordinate
  type(ESMF_Time) :: current               ! current date (from clock)
  type(ESMF_TimeInterval) :: day_step_size ! day step size
  type(ESMF_TimeInterval) :: step_size     ! timestep size
  logical         :: run_length_specified = .false.
  character(len=*), parameter :: sub = 'rtm::timemgr_init'

  !
  dtime = real(dtime_in)
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

  ! Print configuration summary to log file (stdout).
  if (masterproc) call timemgr_print()

  timemgr_set = .true.

end subroutine timemgr_init

!=========================================================================================

subroutine init_clock( start_date, ref_date, curr_date, stop_date )

  ! Initialize the clock based on the start_date, ref_date, and curr_date
  ! as well as the settings from the namelist specifying the time to stop
  !
  type(ESMF_Time), intent(in) :: start_date  ! start date for run
  type(ESMF_Time), intent(in) :: ref_date    ! reference date for time coordinate
  type(ESMF_Time), intent(in) :: curr_date   ! current date (equal to start_date)
  type(ESMF_Time), intent(in) :: stop_date   ! stop date for run
  !
  character(len=*), parameter :: sub = 'rtm::init_clock'
  type(ESMF_TimeInterval) :: step_size       ! timestep size
  type(ESMF_Time) :: current     ! current date (from clock)
  integer :: yr, mon, day, tod   ! Year, month, day, and second as integers
  integer :: rc                  ! return code
  !
  call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

  ! Initialize the clock

  tm_clock = ESMF_ClockCreate(name="RTM Time-manager clock", timeStep=step_size, startTime=start_date, &
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


  ! Set the time by an integer as YYYYMMDD and integer seconds in the day
  !
  integer, intent(in) :: ymd            ! Year, month, day YYYYMMDD
  integer, intent(in) :: tod            ! Time of day in seconds
  character(len=*), intent(in) :: desc  ! Description of time to set
  !
  type(ESMF_Time) :: TimeSetymd         ! Return value
  !
  character(len=*), parameter :: sub = 'rtm::TimeSetymd'
  integer :: yr, mon, day          ! Year, month, day as integers
  integer :: rc                    ! return code
  !
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

! Subprogram not used integer function TimeGetymd( date, tod )
! Subprogram not used 
! Subprogram not used   ! Get the date and time of day in ymd from ESMF Time.
! Subprogram not used   !
! Subprogram not used   type(ESMF_Time), intent(inout) :: date ! Input date to convert to ymd
! Subprogram not used   integer, intent(out), optional :: tod  ! Time of day in seconds
! Subprogram not used   !
! Subprogram not used   character(len=*), parameter :: sub = 'rtm::TimeGetymd'
! Subprogram not used   integer :: yr, mon, day
! Subprogram not used   integer :: rc                          ! return code
! Subprogram not used   !
! Subprogram not used   call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, rc=rc)
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used   TimeGetymd = yr*10000 + mon*100 + day
! Subprogram not used   if ( present( tod ) )then
! Subprogram not used      call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used      call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used   end if
! Subprogram not used   if ( yr < 0 )then
! Subprogram not used      write(iulog,*) sub//': error year is less than zero', yr
! Subprogram not used      call shr_sys_abort
! Subprogram not used   end if
! Subprogram not used end function TimeGetymd

!=========================================================================================

! Subprogram not used subroutine timemgr_restart(ncid, flag)
! Subprogram not used 
! Subprogram not used   ! Read/Write information needed on restart to a netcdf file. 
! Subprogram not used   !
! Subprogram not used   type(file_desc_t), intent(inout) :: ncid  ! netcdf id
! Subprogram not used   character(len=*) , intent(in) :: flag     ! 'read' or 'write'
! Subprogram not used   !
! Subprogram not used   logical :: run_length_specified = .false.
! Subprogram not used   integer :: rc                  ! return code
! Subprogram not used   integer :: yr, mon, day, tod   ! Year, month, day, and second as integers
! Subprogram not used   logical :: readvar             ! determine if variable is on initial file
! Subprogram not used   integer :: rst_caltype         ! calendar type
! Subprogram not used   type(ESMF_Time) :: start_date  ! start date for run
! Subprogram not used   type(ESMF_Time) :: stop_date   ! stop date for run
! Subprogram not used   type(ESMF_Time) :: ref_date    ! reference date for run
! Subprogram not used   type(ESMF_Time) :: curr_date   ! date of data in restart file
! Subprogram not used   type(ESMF_Time) :: current     ! current date (from clock)
! Subprogram not used   type(ESMF_TimeInterval) :: day_step_size ! day step size
! Subprogram not used   type(ESMF_TimeInterval) :: step_size     ! timestep size
! Subprogram not used   integer, parameter :: noleap = 1
! Subprogram not used   integer, parameter :: gregorian = 2
! Subprogram not used   character(len=135) :: varname
! Subprogram not used   character(len=len(calendar)) :: cal
! Subprogram not used   character(len=*), parameter :: sub = 'timemgr_restart'
! Subprogram not used   !
! Subprogram not used   if (flag == 'write') then
! Subprogram not used      rst_calendar  = calendar
! Subprogram not used   else if (flag == 'read') then
! Subprogram not used      calendar = rst_calendar
! Subprogram not used   end if
! Subprogram not used   varname = 'timemgr_rst_type'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='calendar type', units='unitless', flag_meanings=(/ "NO_LEAP_C", "GREGORIAN" /), &
! Subprogram not used           flag_values=(/ noleap, gregorian /), ifill_value=uninit_int )
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      if (flag== 'write') then
! Subprogram not used         cal = to_upper(calendar)
! Subprogram not used         if ( trim(cal) == 'NO_LEAP' ) then
! Subprogram not used            rst_caltype = noleap
! Subprogram not used         else if ( trim(cal) == 'GREGORIAN' ) then
! Subprogram not used            rst_caltype = gregorian
! Subprogram not used         else
! Subprogram not used            call shr_sys_abort(sub//'ERROR: unrecognized calendar specified= '//trim(calendar))
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used      call ncd_io(varname=varname, data=rst_caltype, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used      if (flag == 'read') then
! Subprogram not used         if ( rst_caltype == noleap ) then
! Subprogram not used            calendar = 'NO_LEAP'
! Subprogram not used         else if ( rst_caltype == gregorian ) then
! Subprogram not used            calendar = 'GREGORIAN'
! Subprogram not used         else
! Subprogram not used            write(iulog,*)sub,': unrecognized calendar type in restart file: ',rst_caltype
! Subprogram not used            call shr_sys_abort( sub//'ERROR: bad calendar type in restart file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   if (flag == 'write') then
! Subprogram not used      call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, refTime=ref_date, rc=rc )
! Subprogram not used      call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used      rst_step_sec  = dtime
! Subprogram not used      rst_start_ymd = TimeGetymd( start_date, tod=rst_start_tod )
! Subprogram not used      rst_ref_ymd   = TimeGetymd( ref_date,   tod=rst_ref_tod   )
! Subprogram not used      rst_curr_ymd  = TimeGetymd( curr_date,  tod=rst_curr_tod  )
! Subprogram not used   end if
! Subprogram not used   
! Subprogram not used   varname = 'timemgr_rst_step_sec'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='seconds component of timestep size', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int)
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_step_sec, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used      if ( rst_step_sec < 0 .or. rst_step_sec > isecspday ) then
! Subprogram not used         call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   varname = 'timemgr_rst_start_ymd'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='start date', units='YYYYMMDD', ifill_value=uninit_int)
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_start_ymd, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   varname = 'timemgr_rst_start_tod'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='start time of day', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int)
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_start_tod, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used      if ( rst_start_tod < 0 .or. rst_start_tod > isecspday ) then
! Subprogram not used         call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   varname = 'timemgr_rst_ref_ymd'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='reference date', units='YYYYMMDD', ifill_value=uninit_int)
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_ref_ymd, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   varname = 'timemgr_rst_ref_tod'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='reference time of day', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int)
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_ref_tod, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used      if ( rst_start_tod < 0 .or. rst_start_tod > isecspday ) then
! Subprogram not used         call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   varname = 'timemgr_rst_curr_ymd'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='current date', units='YYYYMMDD', ifill_value=uninit_int)
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_curr_ymd, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   varname = 'timemgr_rst_curr_tod'
! Subprogram not used   if (flag == 'define') then
! Subprogram not used      call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
! Subprogram not used           long_name='current time of day', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int )
! Subprogram not used   else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used      call ncd_io(varname=varname, data=rst_curr_tod, &
! Subprogram not used           ncid=ncid, flag=flag, readvar=readvar)
! Subprogram not used      if (flag=='read' .and. .not. readvar) then
! Subprogram not used         if (is_restart()) then
! Subprogram not used            call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
! Subprogram not used         end if
! Subprogram not used      end if
! Subprogram not used      if ( rst_curr_tod < 0 .or. rst_curr_tod > isecspday ) then
! Subprogram not used         call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if (flag == 'read') then
! Subprogram not used 
! Subprogram not used      ! Restart the ESMF time manager using the synclock for ending date.
! Subprogram not used      call timemgr_spmdbcast( )
! Subprogram not used      
! Subprogram not used      ! Initialize calendar from restart info
! Subprogram not used      call init_calendar()
! Subprogram not used      
! Subprogram not used      ! Initialize the timestep from restart info
! Subprogram not used      dtime = rst_step_sec
! Subprogram not used      
! Subprogram not used      ! Initialize start date from restart info
! Subprogram not used      start_date = TimeSetymd( rst_start_ymd, rst_start_tod, "start_date" )
! Subprogram not used      
! Subprogram not used      ! Initialize current date from restart info
! Subprogram not used      curr_date = TimeSetymd( rst_curr_ymd, rst_curr_tod, "curr_date" )
! Subprogram not used      
! Subprogram not used      ! Initialize stop date from sync clock or namelist input
! Subprogram not used      stop_date = TimeSetymd( 99991231, stop_tod, "stop_date" )
! Subprogram not used      
! Subprogram not used      call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
! Subprogram not used      call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')
! Subprogram not used      
! Subprogram not used      call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
! Subprogram not used      call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')
! Subprogram not used      
! Subprogram not used      if    ( stop_ymd /= uninit_int ) then
! Subprogram not used         current = TimeSetymd( stop_ymd, stop_tod, "stop_date" )
! Subprogram not used         if ( current < stop_date ) stop_date = current
! Subprogram not used         run_length_specified = .true.
! Subprogram not used      else if ( nelapse /= uninit_int ) then
! Subprogram not used         if ( nelapse >= 0 ) then
! Subprogram not used            current = curr_date + step_size*nelapse
! Subprogram not used         else
! Subprogram not used            current = curr_date - day_step_size*nelapse
! Subprogram not used         end if
! Subprogram not used         if ( current < stop_date ) stop_date = current
! Subprogram not used         run_length_specified = .true.
! Subprogram not used      end if
! Subprogram not used      if ( .not. run_length_specified ) then
! Subprogram not used         call shr_sys_abort (sub//': Must specify stop_ymd or nelapse')
! Subprogram not used      end if
! Subprogram not used      
! Subprogram not used      ! Error check
! Subprogram not used      if ( stop_date <= start_date ) then
! Subprogram not used         write(iulog,*)sub, ': stop date must be specified later than start date: '
! Subprogram not used         call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used         write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used         call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used         write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used         call shr_sys_abort
! Subprogram not used      end if
! Subprogram not used      if ( curr_date >= stop_date ) then
! Subprogram not used         write(iulog,*)sub, ': stop date must be specified later than current date: '
! Subprogram not used         call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used         write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used         call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
! Subprogram not used         write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
! Subprogram not used         call shr_sys_abort
! Subprogram not used      end if
! Subprogram not used      
! Subprogram not used      ! Initialize ref date from restart info
! Subprogram not used      ref_date = TimeSetymd( rst_ref_ymd, rst_ref_tod, "ref_date" )
! Subprogram not used      
! Subprogram not used      ! Initialize clock 
! Subprogram not used      call init_clock( start_date, ref_date, curr_date, stop_date )
! Subprogram not used      
! Subprogram not used      ! Set flag that this is the first timestep of the restart run.
! Subprogram not used      tm_first_restart_step = .true.
! Subprogram not used      
! Subprogram not used      ! Print configuration summary to log file (stdout).
! Subprogram not used      if (masterproc) call timemgr_print()
! Subprogram not used 
! Subprogram not used      timemgr_set = .true.
! Subprogram not used 
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used end subroutine timemgr_restart

!=========================================================================================

subroutine init_calendar( )

  !---------------------------------------------------------------------------------
  ! Initialize calendar
  !
  ! Local variables
  !
  character(len=*), parameter :: sub = 'rtm::init_calendar'
  type(ESMF_CalKind_Flag) :: cal_type        ! calendar type
  character(len=len(calendar)) :: caltmp
  integer :: rc                              ! return code
  !---------------------------------------------------------------------------------

  caltmp = to_upper(calendar)
  if ( trim(caltmp) == 'NO_LEAP' ) then
     cal_type = ESMF_CALKIND_NOLEAP
  else if ( trim(caltmp) == 'GREGORIAN' ) then
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
  character(len=*), parameter :: sub = 'rtm::timemgr_print'
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

  write(iulog,*)' ******** RTM Time Manager Configuration ********'

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

  write(iulog,*)' ************************************************'

end subroutine timemgr_print

!=========================================================================================

subroutine advance_timestep()

  ! Increment the timestep number.

  character(len=*), parameter :: sub = 'rtm::advance_timestep'
  integer :: rc
  
  call ESMF_ClockAdvance( tm_clock, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockAdvance')

  tm_first_restart_step = .false.
  
end subroutine advance_timestep

!=========================================================================================

! Subprogram not used subroutine get_clock( clock )
! Subprogram not used 
! Subprogram not used   ! Return the ESMF clock
! Subprogram not used 
! Subprogram not used   type(ESMF_Clock), intent(inout) :: clock
! Subprogram not used 
! Subprogram not used   character(len=*), parameter :: sub = 'rtm::get_clock'
! Subprogram not used   type(ESMF_TimeInterval) :: step_size
! Subprogram not used   type(ESMF_Time) :: start_date, stop_date, ref_date
! Subprogram not used   integer :: rc
! Subprogram not used 
! Subprogram not used   call ESMF_ClockGet( tm_clock, timeStep=step_size, startTime=start_date, &
! Subprogram not used                       stoptime=stop_date, reftime=ref_date, rc=rc )
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used   call ESMF_ClockSet(clock, timeStep=step_size, startTime=start_date, &
! Subprogram not used                             stoptime=stop_date, reftime=ref_date, rc=rc)
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_ClockSet')
! Subprogram not used 
! Subprogram not used end subroutine get_clock

!=========================================================================================

! Subprogram not used integer function get_step_size()
! Subprogram not used 
! Subprogram not used   ! Return the step size in seconds.
! Subprogram not used   
! Subprogram not used   character(len=*), parameter :: sub = 'rtm::get_step_size'
! Subprogram not used   type(ESMF_TimeInterval) :: step_size       ! timestep size
! Subprogram not used   integer :: rc
! Subprogram not used   
! Subprogram not used   call ESMF_ClockGet(tm_clock, timeStep=step_size, rc=rc)
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used   call ESMF_TimeIntervalGet(step_size, s=get_step_size, rc=rc)
! Subprogram not used   call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')
! Subprogram not used   
! Subprogram not used end function get_step_size

!=========================================================================================

integer function get_nstep()

  ! Return the timestep number.

   character(len=*), parameter :: sub = 'rtm::get_nstep'
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

   character(len=*), parameter :: sub = 'rtm::get_curr_date'
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

subroutine get_prev_date(yr, mon, day, tod)

! Return date components valid at beginning of current timestep.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'rtm::get_prev_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_prev_date

!=========================================================================================

! Subprogram not used subroutine get_start_date(yr, mon, day, tod)
! Subprogram not used 
! Subprogram not used    ! Return date components valid at beginning of initial run.
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    character(len=*), parameter :: sub = 'rtm::get_start_date'
! Subprogram not used    integer :: rc
! Subprogram not used    type(ESMF_Time) :: date
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet(tm_clock, startTime=date, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used end subroutine get_start_date

!=========================================================================================

! Subprogram not used subroutine get_ref_date(yr, mon, day, tod)
! Subprogram not used 
! Subprogram not used ! Return date components of the reference date.
! Subprogram not used 
! Subprogram not used ! Arguments
! Subprogram not used    integer, intent(out) ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used    character(len=*), parameter :: sub = 'rtm::get_ref_date'
! Subprogram not used    integer :: rc
! Subprogram not used    type(ESMF_Time) :: date
! Subprogram not used !-----------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet(tm_clock, refTime=date, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used 
! Subprogram not used    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_TimeGet')
! Subprogram not used 
! Subprogram not used end subroutine get_ref_date

!=========================================================================================

subroutine get_curr_time(days, seconds)

! Return time components valid at end of current timestep.
! Current time is the time interval between the current date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'rtm::get_curr_time'
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
   character(len=*), parameter :: sub = 'rtm::get_prev_time'
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

! Subprogram not used function get_calendar()
! Subprogram not used 
! Subprogram not used    ! Return calendar
! Subprogram not used 
! Subprogram not used    character(len=ESMF_MAXSTR) :: get_calendar
! Subprogram not used 
! Subprogram not used    get_calendar = calendar
! Subprogram not used 
! Subprogram not used end function get_calendar

!=========================================================================================
 
! Subprogram not used function is_end_curr_day()
! Subprogram not used 
! Subprogram not used    ! Return true if current timestep is last timestep in current day.
! Subprogram not used    logical :: is_end_curr_day
! Subprogram not used 
! Subprogram not used    integer ::&
! Subprogram not used       yr,    &! year
! Subprogram not used       mon,   &! month
! Subprogram not used       day,   &! day of month
! Subprogram not used       tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    call get_curr_date(yr, mon, day, tod)
! Subprogram not used    is_end_curr_day = (tod == 0)
! Subprogram not used 
! Subprogram not used end function is_end_curr_day

!=========================================================================================

! Subprogram not used logical function is_end_curr_month()
! Subprogram not used 
! Subprogram not used   ! Return true if current timestep is last timestep in current month.
! Subprogram not used    integer ::  yr, mon, day, tod     ! time of day (seconds past 0Z)
! Subprogram not used 
! Subprogram not used    call get_curr_date(yr, mon, day, tod)
! Subprogram not used    is_end_curr_month = (day == 1  .and.  tod == 0)
! Subprogram not used 
! Subprogram not used end function is_end_curr_month

!=========================================================================================

! Subprogram not used logical function is_first_step()
! Subprogram not used 
! Subprogram not used   ! Return true on first step of initial run only.
! Subprogram not used    character(len=*), parameter :: sub = 'rtm::is_first_step'
! Subprogram not used    integer :: rc
! Subprogram not used    integer :: nstep
! Subprogram not used    integer(ESMF_KIND_I8) :: step_no
! Subprogram not used 
! Subprogram not used    call ESMF_ClockGet( tm_clock, advanceCount=step_no, rc=rc )
! Subprogram not used    call chkrc(rc, sub//': error return from ESMF_ClockGet')
! Subprogram not used    nstep = step_no
! Subprogram not used    is_first_step = (nstep == 0)
! Subprogram not used 
! Subprogram not used end function is_first_step

!=========================================================================================

! Subprogram not used logical function is_first_restart_step()
! Subprogram not used 
! Subprogram not used    ! Return true on first step of restart run only.
! Subprogram not used    is_first_restart_step = tm_first_restart_step
! Subprogram not used 
! Subprogram not used end function is_first_restart_step

!=========================================================================================

! Subprogram not used logical function is_last_step()
! Subprogram not used 
! Subprogram not used   ! Return true on last timestep.
! Subprogram not used    character(len=*), parameter :: sub = 'rtm::is_last_step'
! Subprogram not used    type(ESMF_Time) :: stop_date
! Subprogram not used    type(ESMF_Time) :: curr_date
! Subprogram not used    type(ESMF_TimeInterval) :: time_step
! Subprogram not used    integer :: rc
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

subroutine chkrc(rc, mes)
   integer, intent(in)          :: rc   ! return code from time management library
   character(len=*), intent(in) :: mes  ! error message
   if ( rc == ESMF_SUCCESS ) return
   write(iulog,*) mes
   call shr_sys_abort ('CHKRC')
end subroutine chkrc

!=========================================================================================

function to_upper(str)

  ! Convert character string to upper case. Use achar and iachar intrinsics
  ! to ensure use of ascii collating sequence.
  character(len=*), intent(in) :: str ! String to convert to upper case
  character(len=len(str))      :: to_upper

  integer :: i                ! Index
  integer :: aseq             ! ascii collating sequence
  character(len=1) :: ctmp    ! Character temporary
  
  do i = 1, len(str)
     ctmp = str(i:i)
     aseq = iachar(ctmp)
     if ( aseq >= 97  .and.  aseq <= 122 ) ctmp = achar(aseq - 32)
     to_upper(i:i) = ctmp
  end do
  
end function to_upper

!=========================================================================================

! Subprogram not used logical function is_restart( )
! Subprogram not used   ! Determine if restart run
! Subprogram not used   if (nsrest == nsrContinue) then
! Subprogram not used      is_restart = .true.
! Subprogram not used   else
! Subprogram not used      is_restart = .false.
! Subprogram not used   end if
! Subprogram not used end function is_restart

!=========================================================================================

subroutine timemgr_finalize( )
   !
   ! call ESMF_ClockDestroy to clean up ESMF clock memory
   !
   implicit none
   
   integer                     :: rc    ! return code
   character(len=*), parameter :: sub = 'rtm::timemgr_finalize'
   !
   ! tm_clock is a module variable
   !

   ! 
   ! FIX(SPM, 05222014) if you try to compile this with the ESMF interfaces, 
   ! you will get an intel error "This is not a field name that is defined
   ! in the encompassing structure.".  Furthermore, with CME tests,
   ! when building the MCT version we have defined both -DUSE_ESMF_LIB -DMCT_INTERFACE 
   ! during compile.
   ! 


   call ESMF_ClockDestroy( tm_clock, rc )
   call chkrc(rc, sub//': error return from ESMF_ClockDestory')



end subroutine timemgr_finalize

!=========================================================================================

subroutine timemgr_spmdbcast( )

  integer :: ier

  call mpi_bcast (dtime, 1, MPI_INTEGER, 0, mpicom_rof, ier)

end subroutine timemgr_spmdbcast

end module RtmTimeManager
