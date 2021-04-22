module chem_surfvals

!-----------------------------------------------------------------------------------
! Purpose: Provides greenhouse gas (ghg) values at the Earth's surface.
!          These values may be time dependent.
!
! Author: Brian Eaton (assembled module from existing scattered code pieces)
!-----------------------------------------------------------------------------------

   use shr_kind_mod,   only: r8=>shr_kind_r8
   use spmd_utils,     only: masterproc
   use time_manager,   only: get_curr_date, get_start_date, is_end_curr_day, &
                             timemgr_datediff, get_curr_calday
   use abortutils,     only: endrun
   use netcdf
   use error_messages, only: handle_ncerr  
   use cam_logfile,    only: iulog
   use m_types,        only: time_ramp
   use constituents,   only: pcnst

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make default access private
   save

! Public methods
   public ::&
      chem_surfvals_readnl,  &! read namelist input
      chem_surfvals_init,    &! initialize options that depend on namelist input
      chem_surfvals_set,     &! set ghg surface values when scenario_ghg is 'RAMPED' or 'CHEM_LBC_FILE'
      chem_surfvals_get,     &! return surface values for: CO2VMR, CO2MMR, CH4VMR
                              ! N2OVMR, F11VMR, and F12VMR
      chem_surfvals_co2_rad   ! return co2 for radiation

   public :: flbc_list

! Private module data

   ! Default values for namelist variables -- now set by build-namelist
   real(r8) :: o2mmr = .23143_r8               ! o2 mass mixing ratio
   real(r8) :: co2vmr_rad = -1.0_r8            ! co2 vmr override for radiation
   real(r8) :: co2vmr = -1.0_r8                ! co2   volume mixing ratio 
   real(r8) :: n2ovmr = -1.0_r8                ! n2o   volume mixing ratio 
   real(r8) :: ch4vmr = -1.0_r8                ! ch4   volume mixing ratio 
   real(r8) :: f11vmr = -1.0_r8                ! cfc11 volume mixing ratio 
   real(r8) :: f12vmr = -1.0_r8                ! cfc12 volume mixing ratio 
   character(len=16) :: scenario_ghg = 'FIXED' ! 'FIXED','RAMPED' or 'RAMP_CO2_ONLY'
   integer  :: rampYear_ghg = 0                ! ramped gases fixed at this year (if > 0)
   character(len=256) :: bndtvghg = ' '        ! filename for ramped data
   integer  :: ramp_co2_start_ymd = 0          ! start date for co2 ramping (yyyymmdd)
   real(r8) :: ramp_co2_annual_rate = 1.0_r8      ! % amount of co2 ramping per yr; default is 1% 
   real(r8) :: ramp_co2_cap = -9999.0_r8          ! co2 ramp cap if rate>0, floor otherwise 
                                               ! as multiple or fraction of inital value
                                               ! ex. 4.0 => cap at 4x initial co2 setting 
   integer  :: ghg_yearStart_model = 0         ! model start year
   integer  :: ghg_yearStart_data  = 0         ! data  start year   

   logical  :: ghg_use_calendar                ! true => data year = model year
   logical :: doRamp_ghg    ! true => turn on ramping for ghg
   logical :: ramp_just_co2 ! true => ramping to be done just for co2 and not other ghg's
   integer :: fixYear_ghg   ! year at which Ramped gases are fixed
   integer :: co2_start     ! date at which co2 begins ramping
   real(r8) :: co2_daily_factor    ! daily multiplier to achieve annual rate of co2 ramp
   real(r8) :: co2_limit    ! value of co2vmr where ramping ends
   real(r8) :: co2_base     ! initial co2 volume mixing ratio, before any ramping
   integer :: ntim = -1               ! number of yearly data values
   integer,  allocatable :: yrdata(:) ! yearly data values
   real(r8), allocatable :: co2(:)    ! co2 mixing ratios in ppmv 
   real(r8), allocatable :: ch4(:)    ! ppbv
   real(r8), allocatable :: n2o(:)    ! ppbv
   real(r8), allocatable :: f11(:)    ! pptv
   real(r8), allocatable :: f12(:)    ! pptv
   real(r8), allocatable :: adj(:)    ! unitless adjustment factor for f11 & f12
   
   ! fixed lower boundary 
   
   character(len=256) :: flbc_file = ' '
   character(len=16)  :: flbc_list(pcnst) = ''
   type(time_ramp)    :: flbc_timing     != time_ramp( "CYCLICAL",  19970101, 0 )

!=========================================================================================
contains
!=========================================================================================

subroutine chem_surfvals_readnl(nlfile)

   ! Read chem_surfvals_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=*), parameter :: subname = 'chem_surfvals_readnl'
   
   character(len=8)   :: flbc_type = 'CYCLICAL'     ! 'CYCLICAL' | 'SERIAL' | 'FIXED'
   integer            :: flbc_cycle_yr = 0
   integer            :: flbc_fixed_ymd = 0
   integer            :: flbc_fixed_tod = 0

   namelist /chem_surfvals_nl/ co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, &
                               co2vmr_rad, scenario_ghg, rampyear_ghg, bndtvghg, &
                               ramp_co2_start_ymd, ramp_co2_annual_rate, ramp_co2_cap, &
                               ghg_yearStart_model, ghg_yearStart_data
   ! waccm/cam-chem naemlist
   namelist /chem_surfvals_nl/ flbc_type, flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, flbc_list, flbc_file

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'chem_surfvals_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, chem_surfvals_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpibcast (co2vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (n2ovmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (ch4vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (f11vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (f12vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (co2vmr_rad,                      1,   mpir8, 0, mpicom)
   call mpibcast (scenario_ghg,    len(scenario_ghg), mpichar, 0, mpicom)
   call mpibcast (rampyear_ghg,                    1,  mpiint, 0, mpicom)
   call mpibcast (bndtvghg,            len(bndtvghg), mpichar, 0, mpicom)
   call mpibcast (ramp_co2_start_ymd,              1,  mpiint, 0, mpicom)
   call mpibcast (ramp_co2_annual_rate,            1,   mpir8, 0, mpicom)
   call mpibcast (ramp_co2_cap,                    1,   mpir8, 0, mpicom)
   call mpibcast (ghg_yearstart_model,             1,  mpiint, 0, mpicom)
   call mpibcast (ghg_yearstart_data,              1,  mpiint, 0, mpicom)
   
   ! waccm/cam-chem fixed lower boundary 
   
   call mpibcast (flbc_type,         len(flbc_type),                  mpichar, 0, mpicom)
   call mpibcast (flbc_cycle_yr,     1,                               mpiint,  0, mpicom)
   call mpibcast (flbc_fixed_ymd,    1,                               mpiint,  0, mpicom)
   call mpibcast (flbc_fixed_tod,    1,                               mpiint,  0, mpicom)
   call mpibcast (flbc_list,         len(flbc_list(1))*pcnst,         mpichar, 0, mpicom)
   call mpibcast (flbc_file,         len(flbc_file),                  mpichar, 0, mpicom)


   flbc_timing%type      = flbc_type
   flbc_timing%cycle_yr  = flbc_cycle_yr
   flbc_timing%fixed_ymd = flbc_fixed_ymd
   flbc_timing%fixed_tod = flbc_fixed_tod

   if ( len_trim(bndtvghg) > 0 .and. len_trim(flbc_file) > 0 ) then
      call endrun('chem_surfvals_readnl: Cannot specify both bndtvghg and flbc_file ')
   endif

   if (co2vmr_rad > 0._r8) then
      if (masterproc) &
         write(iulog,*) trim(subname)//': co2vmr_rad override is set to ', co2vmr_rad
   end if

end subroutine chem_surfvals_readnl

!================================================================================================

subroutine chem_surfvals_init()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! Set surface values at initial time.
! N.B. This routine must be called after the time manager has been initialized
!      since chem_surfvals_set calls time manager methods.
! 
! Author: B. Eaton - merged code from parse_namelist and rampnl_ghg.
! 
!-----------------------------------------------------------------------

   use infnan,  only : posinf, assignment(=)
   use mo_flbc, only : flbc_inti

   !---------------------------Local variables-----------------------------
   integer :: yr, mon, day, ncsec
   !-----------------------------------------------------------------------

   if (scenario_ghg == 'FIXED') then
      doRamp_ghg = .false.
      ramp_just_co2 = .false.
      if (masterproc) &
         write(iulog,*)'chem_surfvals_init: ghg surface values are fixed as follows'

   else if (scenario_ghg == 'RAMPED') then
      doRamp_ghg = .true.
      ramp_just_co2 = .false.
      call ghg_ramp_read

      fixYear_ghg = rampYear_ghg     ! set private member to namelist var
      if (masterproc) then
         if ( fixYear_ghg > 0 ) then
            write(iulog,*) '  FIXED values from year ',fixYear_ghg
         else
            write(iulog,*) '  RAMPED values initialized to'
         end if
      end if
      call chem_surfvals_set()

   else if (scenario_ghg == 'RAMP_CO2_ONLY') then
      if(ramp_co2_start_ymd == 0) then
         ! by default start the ramp at the initial run time
         call get_start_date(yr, mon, day, ncsec)
         ramp_co2_start_ymd = yr*10000 + mon*100 + day
      end if
      co2_start = ramp_co2_start_ymd

      if(ramp_co2_annual_rate <= -100.0_r8) then
         write(iulog,*) 'RAMP_CO2:  invalid ramp_co2_annual_rate= ',ramp_co2_annual_rate
         call endrun ('chem_surfvals_init: RAMP_CO2_ANNUAL_RATE must be greater than -100.0')
      end if

      doRamp_ghg = .true.
      ramp_just_co2 = .true.
      co2_base = co2vmr        ! save initial setting 
      if (masterproc) &
           write(iulog,*) '  RAMPED values initialized to'

      co2_daily_factor = (ramp_co2_annual_rate*0.01_r8+1.0_r8)**(1.0_r8/365.0_r8)

      if(ramp_co2_cap > 0.0_r8) then  
         co2_limit = ramp_co2_cap * co2_base
      else                                  ! if no cap/floor specified, provide default
         if(ramp_co2_annual_rate < 0.0_r8) then
            co2_limit = 0.0_r8
         else
            co2_limit = posinf
         end if
      end if
      if((ramp_co2_annual_rate<0.0_r8 .and. co2_limit>co2_base) .or. &
         (ramp_co2_annual_rate>0.0_r8 .and. co2_limit<co2_base)) then
         write(iulog,*) 'RAMP_CO2: ramp_co2_cap is unreachable'
         write(iulog,*) 'RAMP_CO2: ramp_co2_annual_rate= ',ramp_co2_annual_rate,' ramp_co2_cap= ',ramp_co2_cap
         call endrun('chem_surfvals_init:  ramp_co2_annual_rate and ramp_co2_cap incompatible')
      end if

      call chem_surfvals_set()
   else if (scenario_ghg == 'CHEM_LBC_FILE') then
      ! set by lower boundary conditions file
      call flbc_inti( flbc_file, flbc_list, flbc_timing, co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr )
      call chem_surfvals_set()
   else
      call endrun ('chem_surfvals_init: input namelist SCENARIO_GHG must be set to either FIXED, RAMPED, RAMP_CO2_ONLY, &
                    or CHEM_LBC_FILE')
   endif

   if (masterproc) then
      write(iulog,*) '  co2 volume mixing ratio = ',co2vmr
      write(iulog,*) '  ch4 volume mixing ratio = ',ch4vmr
      write(iulog,*) '  n2o volume mixing ratio = ',n2ovmr
      write(iulog,*) '  f11 volume mixing ratio = ',f11vmr
      write(iulog,*) '  f12 volume mixing ratio = ',f12vmr
   end if

end subroutine chem_surfvals_init

!=========================================================================================

! Subprogram not used subroutine ghg_ramp_read()
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Read ramped greenhouse gas surface data.  
! Subprogram not used ! 
! Subprogram not used ! Author: T. Henderson
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    use ioFileMod, only: getfil
! Subprogram not used    use mpishorthand, only: mpicom, mpiint, mpir8
! Subprogram not used    character(len=*), parameter :: subname = 'ghg_ramp_read'
! Subprogram not used 
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used    integer :: ncid
! Subprogram not used    integer :: co2_id
! Subprogram not used    integer :: ch4_id
! Subprogram not used    integer :: n2o_id
! Subprogram not used    integer :: f11_id
! Subprogram not used    integer :: f12_id
! Subprogram not used    integer :: adj_id
! Subprogram not used    integer :: date_id
! Subprogram not used    integer :: time_id
! Subprogram not used    integer :: ierror
! Subprogram not used    character(len=256) :: locfn          ! netcdf local filename to open
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used      call getfil (bndtvghg, locfn, 0)
! Subprogram not used      call handle_ncerr( nf90_open (trim(locfn), NF90_NOWRITE, ncid),subname,309)
! Subprogram not used 
! Subprogram not used      write(iulog,*)'GHG_RAMP_READ:  reading ramped greenhouse gas surface data from file ',trim(locfn)
! Subprogram not used 
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'date', date_id ),subname,313)
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'CO2', co2_id ),subname,314)
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'CH4', ch4_id ),subname,315)
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'N2O', n2o_id ),subname,316)
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'f11', f11_id ),subname,317)
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'f12', f12_id ),subname,318)
! Subprogram not used      call handle_ncerr( nf90_inq_varid( ncid, 'adj', adj_id ),subname,319)
! Subprogram not used      call handle_ncerr( nf90_inq_dimid( ncid, 'time', time_id ),subname,320)
! Subprogram not used      call handle_ncerr( nf90_inquire_dimension( ncid, time_id, len=ntim ),subname,321)
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used    call mpibcast (ntim, 1, mpiint, 0, mpicom)
! Subprogram not used    ! these arrays are never deallocated
! Subprogram not used    allocate ( yrdata(ntim), co2(ntim), ch4(ntim), n2o(ntim),    &
! Subprogram not used                  f11(ntim), f12(ntim), adj(ntim), stat=ierror )
! Subprogram not used    if (ierror /= 0) then
! Subprogram not used      write(iulog,*)'GHG_RAMP_READ:  ERROR, allocate() failed!'
! Subprogram not used      call endrun
! Subprogram not used    endif
! Subprogram not used    if (masterproc) then
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, date_id, yrdata ),subname,335)
! Subprogram not used      yrdata = yrdata / 10000
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, co2_id, co2 ),subname,337)
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, ch4_id, ch4 ),subname,338)
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, n2o_id, n2o ),subname,339)
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, f11_id, f11 ),subname,340)
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, f12_id, f12 ),subname,341)
! Subprogram not used      call handle_ncerr( nf90_get_var (ncid, adj_id, adj ),subname,342)
! Subprogram not used      call handle_ncerr( nf90_close (ncid),subname,343)
! Subprogram not used      write(iulog,*)'GHG_RAMP_READ:  successfully read ramped greenhouse gas surface data from years ',&
! Subprogram not used 	yrdata(1),' through ',yrdata(ntim)
! Subprogram not used    endif
! Subprogram not used    call mpibcast (co2, ntim, mpir8, 0, mpicom)
! Subprogram not used    call mpibcast (ch4, ntim, mpir8, 0, mpicom)
! Subprogram not used    call mpibcast (n2o, ntim, mpir8, 0, mpicom)
! Subprogram not used    call mpibcast (f11, ntim, mpir8, 0, mpicom)
! Subprogram not used    call mpibcast (f12, ntim, mpir8, 0, mpicom)
! Subprogram not used    call mpibcast (adj, ntim, mpir8, 0, mpicom)
! Subprogram not used    call mpibcast (yrdata, ntim, mpiint, 0, mpicom)
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used 
! Subprogram not used end subroutine ghg_ramp_read

!=========================================================================================

function chem_surfvals_get(name)
  use physconst,    only: mwdry, mwco2

  character(len=*), intent(in) :: name

  real(r8) :: rmwco2 
  real(r8) :: chem_surfvals_get

  rmwco2 = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
  select case (name)
  case ('CO2VMR')
     chem_surfvals_get = co2vmr
  case ('CO2MMR')
     chem_surfvals_get = rmwco2 * co2vmr
  case ('N2OVMR')
     chem_surfvals_get = n2ovmr
  case ('CH4VMR')
     chem_surfvals_get = ch4vmr
  case ('F11VMR')
     chem_surfvals_get = f11vmr
  case ('F12VMR')
     chem_surfvals_get = f12vmr
  case ('O2MMR')
     chem_surfvals_get = o2mmr
  case default
     call endrun('chem_surfvals_get does not know name')
  end select

end function chem_surfvals_get


!=========================================================================================

function chem_surfvals_co2_rad(vmr_in)
 
   ! Return the value of CO2 (as mmr) that is radiatively active.

   ! This method is used by ghg_data to set the prescribed value of CO2 in
   ! the physics buffer.  If the user has set the co2vmr_rad namelist
   ! variable then that value will override either the value set by the
   ! co2vmr namelist variable, or the values time interpolated from a
   ! dataset.
   
   ! This method is also used by cam_history to write the radiatively active
   ! CO2 to the history file.  The optional argument allows returning the
   ! value as vmr.

   use physconst,    only: mwdry, mwco2

   ! Arguments
   logical, intent(in), optional :: vmr_in  ! return CO2 as vmr

   ! Return value
   real(r8) :: chem_surfvals_co2_rad

   ! Local variables
   real(r8) :: convert_vmr      ! convert vmr to desired output
   !-----------------------------------------------------------------------

   ! by default convert vmr to mmr
   convert_vmr = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
   if (present(vmr_in)) then
      ! if request return vmr
      if (vmr_in) convert_vmr = 1.0_r8
   end if

   if (co2vmr_rad > 0._r8) then
      chem_surfvals_co2_rad = convert_vmr * co2vmr_rad
   else                           
      chem_surfvals_co2_rad = convert_vmr * co2vmr     
   end if

end function chem_surfvals_co2_rad

!=========================================================================================

subroutine chem_surfvals_set()

   use ppgrid,         only: begchunk, endchunk
   use mo_flbc,        only: flbc_gmean_vmr, flbc_chk

!---------------------------Local variables-----------------------------

   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]
   
   if ( doRamp_ghg ) then
      if(ramp_just_co2) then
         call chem_surfvals_set_co2()
      else
         call chem_surfvals_set_all()
      end if
   elseif (scenario_ghg == 'CHEM_LBC_FILE') then
      ! set mixing ratios from cam-chem/waccm lbc file 
      call flbc_chk()
      call flbc_gmean_vmr(co2vmr,ch4vmr,n2ovmr,f11vmr,f12vmr)
   endif

   if (masterproc .and. is_end_curr_day()) then
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      write(iulog,*) 'chem_surfvals_set: ncdate= ',ncdate,' co2vmr=',co2vmr

      if (.not. ramp_just_co2 .and. mon==1 .and. day==1) then
         write(iulog,*) 'chem_surfvals_set: ch4vmr=', ch4vmr, ' n2ovmr=', n2ovmr, &
                        ' f11vmr=', f11vmr, ' f12vmr=', f12vmr
      end if

   end if

   return
end subroutine chem_surfvals_set

!=========================================================================================

! Subprogram not used subroutine chem_surfvals_set_all()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Computes greenhouse gas volume mixing ratios via interpolation of
! Subprogram not used ! yearly input data.
! Subprogram not used ! 
! Subprogram not used ! Author: B. Eaton - updated ramp_ghg for use in chem_surfvals module
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use interpolate_data, only: get_timeinterp_factors
! Subprogram not used 
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used 
! Subprogram not used    integer yrmodel           ! model year
! Subprogram not used    integer nyrm              ! year index
! Subprogram not used    integer nyrp              ! year index
! Subprogram not used    integer :: yr, mon, day   ! components of a date
! Subprogram not used    integer :: ncdate         ! current date in integer format [yyyymmdd]
! Subprogram not used    integer :: ncsec          ! current time of day [seconds]
! Subprogram not used 
! Subprogram not used    real(r8) :: calday            ! current calendar day
! Subprogram not used    real(r8) doymodel             ! model day of year
! Subprogram not used    real(r8) doydatam             ! day of year for input data yrdata(nyrm)
! Subprogram not used    real(r8) doydatap             ! day or year for input data yrdata(nyrp)
! Subprogram not used    real(r8) deltat               ! delta time
! Subprogram not used    real(r8) fact1, fact2         ! time interpolation factors
! Subprogram not used    real(r8) cfcscl               ! cfc scale factor for f11
! Subprogram not used 
! Subprogram not used    integer yearRan_model         ! model ran year
! Subprogram not used !
! Subprogram not used ! ---------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    calday = get_curr_calday()
! Subprogram not used    call get_curr_date(yr, mon, day, ncsec)
! Subprogram not used    ncdate = yr*10000 + mon*100 + day
! Subprogram not used !
! Subprogram not used ! determine ghg_use_calendar      
! Subprogram not used !
! Subprogram not used    if ( ghg_yearStart_model > 0 .and. ghg_yearStart_data > 0 ) then
! Subprogram not used       ghg_use_calendar = .false.
! Subprogram not used    else
! Subprogram not used       ghg_use_calendar = .true.
! Subprogram not used    end if
! Subprogram not used !
! Subprogram not used ! determine index into input data
! Subprogram not used !
! Subprogram not used    if ( fixYear_ghg > 0) then
! Subprogram not used       yrmodel  = fixYear_ghg
! Subprogram not used       nyrm = fixYear_ghg - yrdata(1) + 1
! Subprogram not used    else
! Subprogram not used       if ( ghg_use_calendar) then
! Subprogram not used          yrmodel  = yr          
! Subprogram not used          nyrm = yr - yrdata(1) + 1
! Subprogram not used       else 
! Subprogram not used          yearRan_model = yr - ghg_yearStart_model
! Subprogram not used          if ( yearRan_model < 0 ) then
! Subprogram not used             call endrun('chem_surfvals_set_all: incorrect ghg_yearStart_model')
! Subprogram not used          endif
! Subprogram not used          yrmodel  = yearRan_model + ghg_yearStart_data
! Subprogram not used  
! Subprogram not used          nyrm = ghg_yearStart_data + yearRan_model - yrdata(1) + 1
! Subprogram not used       end if
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    nyrp       = nyrm + 1
! Subprogram not used !
! Subprogram not used ! if current date is before yrdata(1), quit
! Subprogram not used !
! Subprogram not used    if (nyrm < 1) then
! Subprogram not used       write(iulog,*)'chem_surfvals_set_all: data time index is out of bounds'
! Subprogram not used       write(iulog,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
! Subprogram not used       call endrun
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used ! if current date later than yrdata(ntim), call endrun.
! Subprogram not used ! if want to use ntim values - uncomment the following lines
! Subprogram not used ! below and comment the call to endrun and previous write
! Subprogram not used !
! Subprogram not used    if (nyrp > ntim) then
! Subprogram not used       call endrun ('chem_surfvals_set_all: error - current date is past the end of valid data')
! Subprogram not used !         write(iulog,*)'chem_surfvals_set_all: using ghg data for ',yrdata(ntim)
! Subprogram not used !         co2vmr = co2(ntim)*1.e-06
! Subprogram not used !         ch4vmr = ch4(ntim)*1.e-09
! Subprogram not used !         n2ovmr = n2o(ntim)*1.e-09
! Subprogram not used !         f11vmr = f11(ntim)*1.e-12*(1.+cfcscl)
! Subprogram not used !         f12vmr = f12(ntim)*1.e-12
! Subprogram not used !         co2mmr = rmwco2 * co2vmr
! Subprogram not used !         return
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used ! determine time interpolation factors, check sanity
! Subprogram not used ! of interpolation factors to within 32-bit roundoff
! Subprogram not used ! assume that day of year is 1 for all input data
! Subprogram not used !
! Subprogram not used    doymodel = yrmodel*365._r8    + calday
! Subprogram not used    doydatam = yrdata(nyrm)*365._r8 + 1._r8
! Subprogram not used    doydatap = yrdata(nyrp)*365._r8 + 1._r8
! Subprogram not used 
! Subprogram not used    call get_timeinterp_factors(.false.,2,doydatam,doydatap, doymodel, &
! Subprogram not used         fact1, fact2,'chem_surfvals')
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! do time interpolation:
! Subprogram not used !   co2     in ppmv
! Subprogram not used !   n2o,ch4 in ppbv
! Subprogram not used !   f11,f12 in pptv
! Subprogram not used !
! Subprogram not used    co2vmr = (co2(nyrm)*fact1 + co2(nyrp)*fact2)*1.e-06_r8
! Subprogram not used    ch4vmr = (ch4(nyrm)*fact1 + ch4(nyrp)*fact2)*1.e-09_r8
! Subprogram not used    n2ovmr = (n2o(nyrm)*fact1 + n2o(nyrp)*fact2)*1.e-09_r8
! Subprogram not used 
! Subprogram not used    cfcscl = (adj(nyrm)*fact1 + adj(nyrp)*fact2)
! Subprogram not used    f11vmr = (f11(nyrm)*fact1 + f11(nyrp)*fact2)*1.e-12_r8*(1._r8+cfcscl)
! Subprogram not used    f12vmr = (f12(nyrm)*fact1 + f12(nyrp)*fact2)*1.e-12_r8
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used end subroutine chem_surfvals_set_all

!=========================================================================================

! Subprogram not used subroutine chem_surfvals_set_co2()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Computes co2 greenhouse gas volume mixing ratio via ramping info 
! Subprogram not used ! provided in namelist var's
! Subprogram not used ! 
! Subprogram not used ! Author: B. Eaton - updated ramp_ghg for use in chem_surfvals module
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used 
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used 
! Subprogram not used    real(r8) :: daydiff             ! number of days of co2 ramping
! Subprogram not used    integer  :: yr, mon, day, ncsec ! components of a date
! Subprogram not used    integer  :: ncdate              ! current date in integer format [yyyymmdd]
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call get_curr_date(yr, mon, day, ncsec)
! Subprogram not used    ncdate = yr*10000 + mon*100 + day
! Subprogram not used 
! Subprogram not used    call timemgr_datediff(co2_start, 0, ncdate, ncsec, daydiff)
! Subprogram not used 
! Subprogram not used    if (daydiff > 0.0_r8) then
! Subprogram not used 
! Subprogram not used       co2vmr = co2_base*(co2_daily_factor)**daydiff
! Subprogram not used 
! Subprogram not used       if(co2_daily_factor < 1.0_r8) then
! Subprogram not used          co2vmr = max(co2vmr,co2_limit)
! Subprogram not used       else
! Subprogram not used          co2vmr = min(co2vmr,co2_limit)
! Subprogram not used       end if
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used end subroutine chem_surfvals_set_co2


!=========================================================================================

end module chem_surfvals
