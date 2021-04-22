!===============================================================================
! SVN $Id: shr_sys_mod.F90 60325 2014-05-16 20:33:31Z santos@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140626/shr/shr_sys_mod.F90 $
!===============================================================================

! Currently supported by all compilers



! Except this combination?








MODULE shr_sys_mod

   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

   use shr_kind_mod  ! defines real & integer kinds
   use shr_mpi_mod   ! wraps MPI layer
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit








   implicit none

! PUBLIC: Public interfaces

   private

   public :: shr_sys_system  ! make a system call
   public :: shr_sys_chdir   ! change current working dir
   public :: shr_sys_getenv  ! get an environment variable
   public :: shr_sys_abort   ! abort a program
   public :: shr_sys_irtc    ! returns real-time clock tick
   public :: shr_sys_sleep   ! have program sleep for a while
   public :: shr_sys_flush   ! flush an i/o buffer
   public :: shr_sys_backtrace   ! print a backtrace, if possible

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_sys_system(str,rcode)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    character(*)        ,intent(in)  :: str    ! system/shell command string
! Subprogram not used    integer(SHR_KIND_IN),intent(out) :: rcode  ! function return error code
! Subprogram not used 
! Subprogram not used    !----- functions -----
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName =   '(shr_sys_system) '
! Subprogram not used    character(*),parameter :: F00     = "('(shr_sys_system) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: an architecture independent system call
! Subprogram not used ! NOTE: 
! Subprogram not used ! - for Catamount (Cray, pheonix at ORNL) there is no system call -- workarounds 
! Subprogram not used !   exist only for simple "rm" and "cp" commands
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used    rcode = 0
! Subprogram not used 
! Subprogram not used    write(s_logunit,F00) 'ERROR: no implementation of system call for this architecture'
! Subprogram not used    call shr_sys_abort(subName//'no implementation of system call for this architecture')
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_sys_system

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_sys_chdir(path, rcode)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    !----- arguments -----
! Subprogram not used    character(*)        ,intent(in)  :: path    ! chdir to this dir
! Subprogram not used    integer(SHR_KIND_IN),intent(out) :: rcode   ! return code
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    integer(SHR_KIND_IN)             :: lenpath ! length of path
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName =   '(shr_sys_chdir) '
! Subprogram not used    character(*),parameter :: F00     = "('(shr_sys_chdir) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: an architecture independent system call
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    lenpath=len_trim(path)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    write(s_logunit,F00) 'ERROR: no implementation of chdir for this architecture'
! Subprogram not used    call shr_sys_abort(subname//'no implementation of chdir for this machine')
! Subprogram not used 
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_sys_chdir

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_getenv(name, val, rcode)

   IMPLICIT none

   !----- arguments -----
   character(*)        ,intent(in)  :: name    ! env var name
   character(*)        ,intent(out) :: val     ! env var value
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenname ! length of env var name
   integer(SHR_KIND_IN)             :: lenval  ! length of env var value
   character(SHR_KIND_CL)           :: tmpval  ! temporary env var value

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_getenv) '
   character(*),parameter :: F00     = "('(shr_sys_getenv) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!-------------------------------------------------------------------------------

!$OMP master


   call get_environment_variable(name=name,value=val,status=rcode)  ! Intrinsic in F2003
!$OMP end master

END SUBROUTINE shr_sys_getenv

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_sys_abort(string,rc)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    character(*)        ,optional :: string  ! error message string
! Subprogram not used    integer(SHR_KIND_IN),optional :: rc      ! error code
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    logical              :: flag
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName =   '(shr_sys_abort) '
! Subprogram not used    character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"
! Subprogram not used 
! Subprogram not used    ! Local version of the string.
! Subprogram not used    ! (Gets a default value if string is not present.)
! Subprogram not used    character(len=shr_kind_cx) :: local_string
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: consistent stopping mechanism
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (present(string)) then
! Subprogram not used       local_string = trim(string)
! Subprogram not used    else
! Subprogram not used       local_string = "Unknown error submitted to shr_sys_abort."
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    call print_error_to_logs("ERROR", local_string)
! Subprogram not used 
! Subprogram not used    call shr_sys_backtrace()
! Subprogram not used 
! Subprogram not used    call shr_mpi_initialized(flag)
! Subprogram not used 
! Subprogram not used    if (flag) then
! Subprogram not used       if (present(rc)) then
! Subprogram not used          call shr_mpi_abort(trim(local_string),rc)
! Subprogram not used       else
! Subprogram not used          call shr_mpi_abort(trim(local_string))
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used   ! A compiler's abort method may print a backtrace or do other nice
! Subprogram not used   ! things, but in fact we can rarely leverage this, because MPI_Abort
! Subprogram not used   ! usually sends SIGTERM to the process, and we don't catch that signal.
! Subprogram not used    call abort()
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_sys_abort

!===============================================================================
!===============================================================================

integer(SHR_KIND_I8) FUNCTION shr_sys_irtc( rate )

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_I8), optional :: rate

   !----- local -----
   integer(SHR_KIND_IN)      :: count
   integer(SHR_KIND_IN)      :: count_rate
   integer(SHR_KIND_IN)      :: count_max
   integer(SHR_KIND_IN),save :: last_count   = -1
   integer(SHR_KIND_I8),save :: count_offset =  0

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_irtc) '
   character(*),parameter :: F00     = "('(shr_sys_irtc) ',4a)"

!-------------------------------------------------------------------------------
! emulates Cray/SGI irtc function (returns clock tick since last reboot)
!-------------------------------------------------------------------------------

   call system_clock(count=count,count_rate=count_rate, count_max=count_max)
   if ( present(rate) ) rate = count_rate
   shr_sys_irtc = count

   !--- adjust for clock wrap-around ---
   if ( last_count /= -1 ) then
     if ( count < last_count ) count_offset = count_offset + count_max
   end if
   shr_sys_irtc = shr_sys_irtc + count_offset
   last_count   = count

END FUNCTION shr_sys_irtc

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_sys_sleep(sec)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    !----- arguments -----
! Subprogram not used    real   (SHR_KIND_R8),intent(in) :: sec  ! number of seconds to sleep
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    integer(SHR_KIND_IN) :: isec   ! integer number of seconds
! Subprogram not used    integer(SHR_KIND_IN) :: rcode  ! return code
! Subprogram not used    character(90)        :: str    ! system call string
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName =   '(shr_sys_sleep) '
! Subprogram not used    character(*),parameter :: F00     = "('(shr_sys_sleep) ',4a)"
! Subprogram not used    character(*),parameter :: F10     = "('sleep ',i8 )"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: Sleep for approximately sec seconds
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    isec = nint(sec)
! Subprogram not used 
! Subprogram not used    if (isec < 0) then
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) 'ERROR: seconds must be > 0, sec=',sec
! Subprogram not used    else if (isec == 0) then
! Subprogram not used       ! Don't consider this an error and don't call system sleep
! Subprogram not used    else
! Subprogram not used       call sleep(isec)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_sys_sleep

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_flush(unit)

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

   !----- local -----
   integer(SHR_KIND_IN) :: ierr  ! error code

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_flush) '
   character(*),parameter :: F00     = "('(shr_sys_flush) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!
! This is probably no longer needed; the "flush" statement is supported by
! all compilers that CESM supports for years now.
!
!-------------------------------------------------------------------------------
   flush(unit)
!
! The following code was originally present, but there's an obvious issue.
! Since shr_sys_flush is usually used to flush output to a log, when it
! returns an error, does it do any good to print that error to the log?
!
!   if (ierr > 0) then
!      write(s_logunit,*) subname,' Flush reports error: ',ierr
!   endif
!

END SUBROUTINE shr_sys_flush

!===============================================================================
!===============================================================================

subroutine shr_sys_backtrace()

  ! This routine uses compiler-specific facilities to print a backtrace to
  ! error_unit (standard error, usually unit 0).


  ! Currently we have no means to request a backtrace from the NAG runtime,
  ! even though it is capable of emitting backtraces itself, if you use the
  ! "-gline" option.

  ! Similarly, PGI has a -traceback option, but no user interface for
  ! requesting a backtrace to be printed.


  flush(error_unit)

end subroutine shr_sys_backtrace

!===============================================================================
!===============================================================================

!
! This routine prints error messages to s_logunit (which is standard output
! for most tasks in CESM) and also to standard error if s_logunit is a
! file.
!
! It also flushes these output units.
!
subroutine print_error_to_logs(error_type, message)
  character(len=*), intent(in) :: error_type, message

  integer, allocatable :: log_units(:)

  integer :: i

  if (s_logunit == output_unit .or. s_logunit == error_unit) then
     ! If the log unit number is standard output or standard error, just
     ! print to that.
     allocate(log_units(1), source=[s_logunit])
  else
     ! Otherwise print the same message to both the log unit and standard
     ! error.
     allocate(log_units(2), source=[error_unit, s_logunit])
  end if

  do i = 1, size(log_units)
     write(log_units(i),*) trim(error_type), ": ", trim(message)
     flush(log_units(i))
  end do

end subroutine print_error_to_logs

!===============================================================================
!===============================================================================

END MODULE shr_sys_mod
