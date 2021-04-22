!===============================================================================
! SVN $Id: shr_file_mod.F90 43960 2013-02-14 00:35:31Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140626/shr/shr_file_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_file_mod.F90 --- Module to handle various file utilily functions.
!
! !DESCRIPTION:
!
! Miscilaneous methods to handle file and directory utilities as well as FORTRAN
! unit control. Also put/get local files into/from archival location
!
! File utilites used with CCSM Message passing:
!
! shr_file_stdio is the main example here, it changes the working directory, 
!                changes stdin and stdout to a given filename.
!
! This is needed because some implementations of MPI with MPMD so that
! each executable can run in a different working directory and redirect
! output to different files.
!
! File name archival convention, eg.
!    call shr_file_put(rcode,"foo","mss:/USER/foo",rtpd=3650)
! is extensible -- the existence of the option file name prefix, eg. "mss:", 
! and optional arguments, eg. rtpd-3650 can be used to access site-specific 
! storage devices.  Based on CCM (atmosphere) getfile & putfile routines, but
! intended to be a more extensible, shared code.
! 
! !REVISION HISTORY:
!   2006-05-08 E. Kluzek, Add in shr_file_mod and getUnit, freeUnif methods.
!   2000-??-?? B. Kauffman, original version circa 2000
!
! !INTERFACE: ------------------------------------------------------------------

MODULE shr_file_mod

! !USES:

   use shr_kind_mod  ! defines kinds
   use shr_sys_mod   ! system calls
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   
   IMPLICIT none

   PRIVATE           ! By default everything is private to this module
   
! !PUBLIC TYPES:               
   
   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_file_put          ! Put a file to an archive location
   public :: shr_file_get          ! Get a file from an archive location
   public :: shr_file_queryPrefix  ! Get prefix type for a filename
   public :: shr_file_getUnit      ! Get a logical unit for reading or writing
   public :: shr_file_freeUnit     ! Free a logical unit
   public :: shr_file_stdio        ! change dir and stdin and stdout
   public :: shr_file_chDir        ! change current working directory
   public :: shr_file_dirio        ! change stdin and stdout
   public :: shr_file_chStdIn      ! change stdin  (attach to a file)
   public :: shr_file_chStdOut     ! change stdout (attach to a file)
   public :: shr_file_setIO        ! open a log file from namelist
   public :: shr_file_setLogUnit   ! Reset the log unit number
   public :: shr_file_setLogLevel  ! Reset the logging debug level
   public :: shr_file_getLogUnit   ! Get the log unit number
   public :: shr_file_getLogLevel  ! Get the logging debug level

! !PUBLIC DATA MEMBERS:

   ! Integer flags for recognized prefixes on file get/put operations
   integer(SHR_KIND_IN), parameter, public :: shr_file_noPrefix      = 0 ! no recognized prefix
   integer(SHR_KIND_IN), parameter, public :: shr_file_nullPrefix    = 1 ! null:
   integer(SHR_KIND_IN), parameter, public :: shr_file_cpPrefix      = 2 ! cp:
   integer(SHR_KIND_IN), parameter, public :: shr_file_mssPrefix     = 3 ! mss:
   integer(SHR_KIND_IN), parameter, public :: shr_file_hpssPrefix    = 4 ! hpss:

!EOP
   !--- unit numbers, users can ask for unit numbers from 0 to min, but getUnit
   !--- won't give a unit below min, users cannot ask for unit number above max
   !--- for backward compatability.  
   !--- eventually, recommend min as hard lower limit (tcraig, 9/2007)
   integer(SHR_KIND_IN),parameter :: shr_file_minUnit = 10      ! Min unit number to give
   integer(SHR_KIND_IN),parameter :: shr_file_maxUnit = 99      ! Max unit number to give
   logical, save :: UnitTag(0:shr_file_maxUnit) = .false. ! Logical units in use

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_put -- Put a file to an archival location.
!
! !DESCRIPTION:
!    a generic, extensible put-local-file-into-archive routine
! USAGE:
!    call shr_file_put(rcode,"foo","/home/user/foo")
!    if ( rcode /= 0 ) call shr_sys_abort( "error copying foo" )
!    call shr_file_put(rcode,"foo","cp:/home/user/foo",remove=.true.)
!    if ( rcode /= 0 ) call shr_sys_abort( "error copying foo" )
!    call shr_file_put(rcode,"foo","mss:/USER/foo",rtpd=3650)
!    if ( rcode /= 0 ) call shr_sys_abort( "error archiving foo to MSS" )
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_put(rcode,loc_fn,rem_fn,passwd,rtpd,async,remove)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer(SHR_KIND_IN),intent(out)         :: rcode  ! return code (non-zero -- error)
! Subprogram not used    character(*),        intent(in)          :: loc_fn ! local filename
! Subprogram not used    character(*),        intent(in)          :: rem_fn ! remote filename
! Subprogram not used    character(*),        intent(in),optional :: passwd ! password
! Subprogram not used    integer(SHR_KIND_IN),intent(in),optional :: rtpd   ! MSS retention period
! Subprogram not used    logical,             intent(in),optional :: async  ! true <=> asynchronous put
! Subprogram not used    logical,             intent(in),optional :: remove ! true <=> rm after put
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----  
! Subprogram not used    integer(SHR_KIND_IN)   :: rtpd2    ! MSS retention period
! Subprogram not used    logical                :: remove2  ! true <=> rm after put
! Subprogram not used    logical                :: async2   ! true <=> asynchronous put
! Subprogram not used    character(SHR_KIND_CL) :: passwd2  ! password
! Subprogram not used    character(SHR_KIND_CL) :: rfn      ! rem_fn without the destination prefix
! Subprogram not used    character(SHR_KIND_CL) :: cmd      ! command sent to system call
! Subprogram not used    integer(SHR_KIND_IN)   :: prefix   ! remote file prefix type
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(shr_file_put) '
! Subprogram not used    character(*),parameter :: F00 =   "('(shr_file_put) ',4a)"
! Subprogram not used    character(*),parameter :: F01 =   "('(shr_file_put) ',a,i3,2a)"
! Subprogram not used    character(*),parameter :: F02 =   "(a,i4)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes: 
! Subprogram not used ! - On some machines the system call will not return a valid error code
! Subprogram not used ! - when things are sent asynchronously, there probably won't be a error code
! Subprogram not used !   returned.
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    remove2 =.false. ; if ( PRESENT(remove )) remove2 = remove
! Subprogram not used    async2  =.true.  ; if ( PRESENT(async  )) async2  = async
! Subprogram not used    passwd2 = " "    ; if ( PRESENT(passwd )) passwd2 = passwd
! Subprogram not used    rtpd2   = 365    ; if ( PRESENT(rtpd   )) rtpd2   = rtpd
! Subprogram not used    rcode   = 0
! Subprogram not used    prefix = shr_file_queryPrefix( rem_fn )
! Subprogram not used 
! Subprogram not used    if ( trim(rem_fn) == trim(loc_fn) ) then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! (remote file name) == (local file name) => do nothing
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       cmd = 'do nothing: remote file = local file = '//trim(loc_fn)
! Subprogram not used       rcode = 0
! Subprogram not used    else if ( prefix == shr_file_cpPrefix .or. prefix == shr_file_noPrefix )then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! put via 1 cp
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       rfn = rem_fn
! Subprogram not used       if ( rem_fn(1:3) == "cp:") rfn = rem_fn(4:len_trim(rem_fn))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       cmd = '/bin/cp -f '//trim(loc_fn)//' '//trim(rfn)
! Subprogram not used       if (remove2) cmd = trim(cmd)//' && /bin/rm -f '//trim(loc_fn)
! Subprogram not used       if (async2 ) cmd = trim(cmd)//' & '
! Subprogram not used       call shr_sys_system(trim(cmd),rcode)
! Subprogram not used 
! Subprogram not used    else if ( prefix == shr_file_mssPrefix )then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! put onto NCAR's MSS
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       if (rtpd2 > 9999) rtpd2 = 9999
! Subprogram not used       write(cmd,F02) '/usr/local/bin/msrcp -period ',rtpd2
! Subprogram not used       if (async2 .and. (.not. remove2) ) cmd = trim(cmd)//' -async '
! Subprogram not used       if (len_trim(passwd2) > 0 ) cmd = trim(cmd)//' -wpwd '//trim(passwd)
! Subprogram not used       cmd = trim(cmd)//' '//trim(loc_fn)//' '//trim(rem_fn)
! Subprogram not used       if (remove2) cmd = trim(cmd)//' && /bin/rm -f '//trim(loc_fn)
! Subprogram not used       if (async2 .and. remove2 ) cmd = trim(cmd)//' & '
! Subprogram not used       call shr_sys_system(trim(cmd),rcode)
! Subprogram not used    else if ( prefix == shr_file_hpssPrefix )then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! put onto LANL's hpss
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       rcode = -1
! Subprogram not used       cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
! Subprogram not used       write(s_logunit,F00) 'ERROR: hpss option not yet implemented'
! Subprogram not used       call shr_sys_abort( subName//'ERROR: hpss option not yet implemented' )
! Subprogram not used    else if ( prefix == shr_file_nullPrefix )then
! Subprogram not used       ! do nothing
! Subprogram not used       cmd = "null prefix => no file archival, do nothing"
! Subprogram not used       rcode = 0
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F01) 'rcode =',rcode,' cmd = ', trim(cmd)
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_file_put

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_get -- Get a file from archival location.
!
! !DESCRIPTION:
!    a generic, extensible get-local-file-from-archive routine
!
! USAGE:
!    call shr_file_get(rcode,"foo","/home/user/foo")
!    if ( rcode /= 0 ) call shr_sys_abort( "error getting file foo" )
!    call shr_file_get(rcode,"foo","cp:/home/user/foo",remove=.true.)
!    if ( rcode /= 0 ) call shr_sys_abort( "error getting file foo" )
!    call shr_file_get(rcode,"foo","mss:/USER/foo",clobber=.true.)
!    if ( rcode /= 0 ) call shr_sys_abort( "error getting file foo from MSS" )
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_get(rcode,loc_fn,rem_fn,passwd,async,clobber)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer(SHR_KIND_IN),intent(out)         :: rcode   ! return code (non-zero means error)
! Subprogram not used    character(*)        ,intent(in)          :: loc_fn  ! local filename
! Subprogram not used    character(*)        ,intent(in)          :: rem_fn  ! remote filename
! Subprogram not used    character(*)        ,intent(in),optional :: passwd  ! password
! Subprogram not used    logical             ,intent(in),optional :: async   ! true <=> asynchronous get
! Subprogram not used    logical             ,intent(in),optional :: clobber ! true <=> clobber existing file
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    logical                :: async2   ! true <=> asynchronous get
! Subprogram not used    logical                :: clobber2 ! true <=> clobber existing file
! Subprogram not used    logical                :: exists   ! true <=> local file a ready exists
! Subprogram not used    character(SHR_KIND_CL) :: passwd2  ! password
! Subprogram not used    character(SHR_KIND_CL) :: rfn      ! rem_fn without the destination prefix
! Subprogram not used    character(SHR_KIND_CL) :: cmd      ! command sent to system call
! Subprogram not used    integer(SHR_KIND_IN)   :: prefix   ! remote file prefix type
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(shr_file_get) '
! Subprogram not used    character(*),parameter :: F00 =   "('(shr_file_get) ',4a)"
! Subprogram not used    character(*),parameter :: F01 =   "('(shr_file_get) ',a,i3,2a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes: 
! Subprogram not used ! - On some machines the system call will not return a valid error code
! Subprogram not used ! - When things are sent asynchronously, there probably won't be a error code
! Subprogram not used !   returned.
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    passwd2  = " "     ; if (PRESENT(passwd )) passwd2  = passwd
! Subprogram not used    async2   = .false. ; if (PRESENT(async  )) async2   = async
! Subprogram not used    clobber2 = .false. ; if (PRESENT(clobber)) clobber2 = clobber
! Subprogram not used    rcode    = 0
! Subprogram not used 
! Subprogram not used    inquire(file=trim(loc_fn),exist=exists)
! Subprogram not used    prefix = shr_file_queryPrefix( rem_fn )
! Subprogram not used 
! Subprogram not used    if ( exists .and. .not. clobber2 ) then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! (file exists) and (don't clobber) => do nothing
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       cmd = 'do nothing: file exists & no-clobber for '//trim(loc_fn)
! Subprogram not used       rcode = 0
! Subprogram not used    else if ( trim(rem_fn) == trim(loc_fn) ) then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! (remote file name) == (local file name) => do nothing
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       cmd = 'do nothing: remote file = local file for '//trim(loc_fn)
! Subprogram not used       rcode = 0
! Subprogram not used    else if ( prefix == shr_file_cpPrefix .or. prefix == shr_file_noPrefix )then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! get via 1 cp
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       rfn = rem_fn  ! remove prefix from this temp file name
! Subprogram not used       if (rem_fn(1:3) == "cp:") rfn = rem_fn(4:len_trim(rem_fn))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       cmd = '/bin/cp -f '//trim(rfn)//' '//trim(loc_fn)
! Subprogram not used       if (async2) cmd = trim(cmd)//' & '
! Subprogram not used       call shr_sys_system(trim(cmd),rcode)
! Subprogram not used 
! Subprogram not used    else if ( prefix == shr_file_mssPrefix )then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! get from NCAR's MSS
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       cmd = '/usr/local/bin/msrcp '
! Subprogram not used       if (async2)   cmd = trim(cmd)//' -async '
! Subprogram not used       cmd = trim(cmd)//' '//trim(rem_fn)//' '//trim(loc_fn)
! Subprogram not used       call shr_sys_system(trim(cmd),rcode)
! Subprogram not used    else if ( prefix == shr_file_hpssPrefix )then
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       ! get from LANL's hpss
! Subprogram not used       !------------------------------------------------------
! Subprogram not used       rcode = -1
! Subprogram not used       cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
! Subprogram not used       write(s_logunit,F00) 'ERROR: hpss option not yet implemented'
! Subprogram not used       call shr_sys_abort( subName//'ERROR: hpss option not yet implemented' )
! Subprogram not used    else if ( prefix == shr_file_nullPrefix )then
! Subprogram not used       ! do nothing
! Subprogram not used       cmd = "null prefix => no file retrieval, do nothing"
! Subprogram not used       rcode = 0
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F01) 'rcode =',rcode,' cmd = ', trim(cmd)
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_file_get

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_queryPrefix -- Get the prefix type from a filepath.
!
! !DESCRIPTION:
!    
! !INTERFACE: ------------------------------------------------------------------  

integer(SHR_KIND_IN) FUNCTION shr_file_queryPrefix( filepath, prefix )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*), intent(in)            :: filepath ! Input filepath
   character(*), intent(out), optional :: prefix   ! Output prefix description

!EOP

   !----- local -----

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

   if (     filepath(1:5) == "null:" )then
      shr_file_queryPrefix = shr_file_nullPrefix
      if ( present(prefix) ) prefix = "null:"
   else if( filepath(1:3) == "cp:" )then
      shr_file_queryPrefix = shr_file_cpPrefix
      if ( present(prefix) ) prefix = "cp:"
   else if( filepath(1:4) == "mss:" )then
      shr_file_queryPrefix = shr_file_mssPrefix
      if ( present(prefix) ) prefix = "mss:"
   else if( filepath(1:5) == "hpss:" )then
      shr_file_queryPrefix = shr_file_hpssPrefix
      if ( present(prefix) ) prefix = "hpss:"
   else
      shr_file_queryPrefix = shr_file_noPrefix
      if ( present(prefix) ) prefix = ""
   end if

END FUNCTION shr_file_queryPrefix

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_getUnit -- Get a free FORTRAN unit number
!
! !DESCRIPTION: Get the next free FORTRAN unit number.
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------  

INTEGER FUNCTION shr_file_getUnit ( unit )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in),optional :: unit ! desired unit number

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n      ! loop index
   logical                :: opened ! If unit opened or not

   !----- formats -----
   character(*),parameter :: subName = '(shr_file_getUnit) '
   character(*),parameter :: F00   = "('(shr_file_getUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (present (unit)) then
      inquire( unit, opened=opened )
      if (unit < 0 .or. unit > shr_file_maxUnit) then
         write(s_logunit,F00) 'invalid unit number request:', unit
         call shr_sys_abort( 'ERROR: bad input unit number' )
      else if (opened .or. UnitTag(unit) .or. unit == 0 .or. unit == 5 &
               .or. unit == 6) then
         write(s_logunit,F00) 'unit number ', unit, ' is already in use'
         call shr_sys_abort( 'ERROR: Input unit number already in use' )
      else
         shr_file_getUnit = unit
         UnitTag (unit)   = .true.
         return
      end if

   else
      ! --- Choose first available unit other than 0, 5, or 6  ------
      do n=shr_file_maxUnit, shr_file_minUnit, -1
         inquire( n, opened=opened )
         if (n == 5 .or. n == 6 .or. opened) then
            cycle
         end if
         if ( .not. UnitTag(n) ) then
            shr_file_getUnit = n
            UnitTag(n)       = .true.
            return
         end if
      end do
   end if

   call shr_sys_abort( subName//': Error: no available units found' )

END FUNCTION shr_file_getUnit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_freeUnit -- Free up a FORTRAN unit number
!
! !DESCRIPTION: Free up the given unit number
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_freeUnit ( unit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: unit  ! unit number to be freed

!EOP

   !----- local -----

   !----- formats -----
   character(*), parameter :: subName = '(shr_file_freeUnit) '
   character(*), parameter :: F00 =   "('(shr_file_freeUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (unit < 0 .or. unit > shr_file_maxUnit) then
      if (s_loglev > 0) write(s_logunit,F00) 'invalid unit number request:', unit
   else if (unit == 0 .or. unit == 5 .or. unit == 6) then
      call shr_sys_abort( subName//': Error: units 0, 5, and 6 must not be freed' )
   else if (UnitTag(unit)) then
      UnitTag (unit) = .false.
   else
      if (s_loglev > 0) write(s_logunit,F00) 'unit ', unit, ' was not in use'
   end if

   return

END SUBROUTINE shr_file_freeUnit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_stdio -- Change working directory, and redirect stdin/stdout
!
! !DESCRIPTION:
!   1) change the cwd (current working directory) and 
!   2) redirect stdin & stdout (units 5 & 6) to named files,
!   where the desired cwd & files are specified by namelist file.
!
!   Normally this is done to work around limitations in the execution syntax
!   of common MPI implementations.  For example, SGI's mpirun syntax is not 
!   flexible enough to allow MPMD models to select different execution
!   directories or to redirect stdin & stdout on the command line.  
!   Such functionality is highly desireable for CCSM purposes.  
!   ie. mpirun can't handle this:
!   1> cd /usr/tmp/jdoe/csm/case01/atm ; atm < atm.parm > atm.log &
!   1> cd /usr/tmp/jdoe/csm/case01/cpl ; cpl < cpl.parm > cpl.log &
!   etc.
!
! ASSUMPTIONS:
! o if the cwd, stdin, or stdout are to be changed, there must be a namelist
!   file in the cwd named <model>_stdio.nml  where <model> is provided via
!   subroutine dummy argument. 
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_stdio(model)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    character(*),intent(in) :: model ! used to construct env varible name
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(shr_file_stdio) '
! Subprogram not used    character(*),parameter :: F00   = "('(shr_file_stdio) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes:
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call shr_file_chdir   (model) ! changes cwd
! Subprogram not used    call shr_file_chStdOut(model) ! open units 5 & 6 to named files
! Subprogram not used    call shr_file_chStdIn (model) ! open units 5 & 6 to named files
! Subprogram not used  
! Subprogram not used END SUBROUTINE shr_file_stdio

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_chdir -- Change working directory.
!
! !DESCRIPTION:
!   change the cwd (current working directory), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_chdir(model, rcodeOut)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    character(*)        ,intent(in)           :: model    ! used to construct env varible name
! Subprogram not used    integer(SHR_KIND_IN),intent(out),optional :: rcodeOut ! Return error code
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    character(SHR_KIND_CL) :: dir      ! directory to cd to
! Subprogram not used    integer  (SHR_KIND_IN) :: rcode    ! Return error code
! Subprogram not used    character(SHR_KIND_CL) :: filename ! namelist file to read
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(shr_file_chdir) '
! Subprogram not used    character(*),parameter :: F00 =   "('(shr_file_chdir) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes:
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call shr_file_stdioReadNL( model, filename, dirOut=dir, rcodeOut=rcode )
! Subprogram not used    if (dir /= "nochange") then
! Subprogram not used       call shr_sys_chdir(dir ,rcode)
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename),", changed cwd to ",trim(dir)
! Subprogram not used    else
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename),", cwd has *not* been changed"
! Subprogram not used       rcode = 1
! Subprogram not used    endif
! Subprogram not used    if ( present(rcodeOut) ) rcodeOut = rcode
! Subprogram not used  
! Subprogram not used END SUBROUTINE shr_file_chdir

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_dirio --- Change stdin and stdout.
!
! !DESCRIPTION:
!   change the stdin & stdout (units 5 & 6), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_dirio(model)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    character(*),intent(in) :: model ! used to construct env varible name
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(shr_file_dirio) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes:
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call shr_file_chStdIn (model)
! Subprogram not used    call shr_file_chStdOut(model)
! Subprogram not used  
! Subprogram not used END SUBROUTINE shr_file_dirio

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_chStdIn -- Change stdin
!
! !DESCRIPTION:
!   change the stdin (unit 5), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_chStdIn( model, NLFilename, rcodeOut )
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    character(*)          ,intent(in)           :: model      ! used to construct env var name
! Subprogram not used    character(SHR_KIND_CL),intent(out),optional :: NLFilename ! open unit 5 to this
! Subprogram not used    integer  (SHR_KIND_IN),intent(out),optional :: rcodeOut   ! return code
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    character(SHR_KIND_CL) :: stdin    ! open unit 5 to this file
! Subprogram not used    character(SHR_KIND_CL) :: nlfile   ! Namelist filename for model to read from
! Subprogram not used    character(SHR_KIND_CL) :: filename ! namelist file to read
! Subprogram not used    integer  (SHR_KIND_IN) :: rcode    ! return code
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(shr_file_chStdIn) '
! Subprogram not used    character(*),parameter :: F00   = "('(shr_file_chStdIn) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes:
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call shr_file_stdioReadNL( model, filename, stdinOut=stdin, &
! Subprogram not used                               nlfileOut=nlfile, rcodeOut=rcode )
! Subprogram not used    if (stdin  /= "nochange") then
! Subprogram not used       open(unit=5,file=stdin ,status='UNKNOWN',iostat=rcode)
! Subprogram not used       if ( rcode /= 0 )then
! Subprogram not used          if (s_loglev > 0) &
! Subprogram not used             write(s_logunit,F00) "read ",trim(filename),': error opening file as unit 5:', &
! Subprogram not used                       trim(nlfile)
! Subprogram not used       else
! Subprogram not used          if (s_loglev > 0) &
! Subprogram not used             write(s_logunit,F00) "read ",trim(filename),': unit 5 connected to ', &
! Subprogram not used                       trim(stdin)
! Subprogram not used       end if
! Subprogram not used    else
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename), &
! Subprogram not used          ': unit 5 has *not* been redirected'
! Subprogram not used    endif
! Subprogram not used    if ( len_trim(nlfile) > 0) then
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename), &
! Subprogram not used          ': read namelist from file:',trim(nlfile)
! Subprogram not used       if ( .not. present(NLFilename) )then
! Subprogram not used          if (s_loglev > 0) write(s_logunit,F00) "error: namelist filename NOT present"
! Subprogram not used          rcode = 7
! Subprogram not used       end if
! Subprogram not used    else
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename),", "
! Subprogram not used       if ( present(NLFilename) )then
! Subprogram not used          if (s_loglev > 0) write(s_logunit,F00) "error: namelist filename present, but null"
! Subprogram not used          rcode = 8
! Subprogram not used       end if
! Subprogram not used    endif
! Subprogram not used    if ( present(NLFilename) ) NLFilename = nlfile
! Subprogram not used    if ( present(rcodeOut)   ) rcodeOut   = rcode
! Subprogram not used  
! Subprogram not used END SUBROUTINE shr_file_chStdIn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_stdout -- Change stdout
!
! !DESCRIPTION:
!   change the stdout (unit 6), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_chStdOut(model,rcodeOut)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used    !--- arguments ---
! Subprogram not used    character(*)        ,intent(in)           :: model    ! used to construct env varible name
! Subprogram not used    integer(SHR_KIND_IN),intent(out),optional :: rcodeOut ! Return error code
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    character(SHR_KIND_CL) :: filename ! namelist file to read
! Subprogram not used    character(SHR_KIND_CL) :: stdout   ! open unit 6 to this file
! Subprogram not used    integer  (SHR_KIND_IN) :: rcode    ! return code
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(shr_file_chStdOut) '
! Subprogram not used    character(*),parameter :: F00   = "('(shr_file_chStdOut) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes:
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call shr_file_stdioReadNL( model, filename, stdoutOut=stdout, &
! Subprogram not used                               rcodeOut=rcode )
! Subprogram not used    if (stdout /= "nochange") then
! Subprogram not used       close(6)
! Subprogram not used       open(unit=6,file=stdout,position='APPEND')
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename),  &
! Subprogram not used          ': unit 6 connected to ',trim(stdout)
! Subprogram not used       call shr_sys_flush(s_logunit)
! Subprogram not used    else
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "read ",trim(filename),  &
! Subprogram not used          ': unit 6 has *not* been redirected'
! Subprogram not used       rcode = 1
! Subprogram not used    endif
! Subprogram not used  
! Subprogram not used    if ( present(rcodeOut)   ) rcodeOut = rcode
! Subprogram not used  
! Subprogram not used END SUBROUTINE shr_file_chStdOut

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_stdioReadNL -- read in stdio namelist
!
! !DESCRIPTION:
!   Read in the stdio namelist for any given model type. Return any of the
! needed input namelist variables as optional arguments. Return "nochange" in
! dir, stdin, or stdout if shouldn't change.
!
! !INTERFACE: ------------------------------------------------------------------  

! Subprogram not used SUBROUTINE shr_file_stdioReadNL( model, filename, dirOut, stdinOut, stdoutOut, &
! Subprogram not used                                  NLFileOut, rcodeOut )
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    character(*)          ,intent(in)           :: model     ! used to construct env varible name
! Subprogram not used    character(SHR_KIND_CL),intent(out)          :: filename  ! nml file to read from unit 5
! Subprogram not used    character(SHR_KIND_CL),intent(out),optional :: NLFileOut ! open unit 6 to this file
! Subprogram not used    character(SHR_KIND_CL),intent(out),optional :: dirOut    ! directory to cd to
! Subprogram not used    character(SHR_KIND_CL),intent(out),optional :: stdinOut  ! open unit 5 to this file
! Subprogram not used    character(SHR_KIND_CL),intent(out),optional :: stdoutOut ! open unit 6 to this file
! Subprogram not used    integer  (SHR_KIND_IN),intent(out),optional :: rcodeOut  ! return code
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    logical                :: exists   ! true iff file exists
! Subprogram not used    character(SHR_KIND_CL) :: dir      ! directory to cd to
! Subprogram not used    character(SHR_KIND_CL) :: stdin    ! open unit 5 to this file
! Subprogram not used    character(SHR_KIND_CL) :: stdout   ! open unit 6 to this file
! Subprogram not used    character(SHR_KIND_CL) :: NLFile   ! namelist file to read seperately
! Subprogram not used    integer  (SHR_KIND_IN) :: rcode    ! return code
! Subprogram not used    integer  (SHR_KIND_IN) :: unit     ! Unit to read from
! Subprogram not used 
! Subprogram not used    namelist / stdio / dir,stdin,stdout,NLFile
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(shr_file_stdioReadNL) '
! Subprogram not used    character(*),parameter :: F00   = "('(shr_file_stdioReadNL) ',4a)"
! Subprogram not used    character(*),parameter :: F01   = "('(shr_file_stdioReadNL) ',2a,i6)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! Notes:
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    rcode  = 0
! Subprogram not used    dir    = "nochange"
! Subprogram not used    stdin  = "nochange"
! Subprogram not used    stdout = "nochange"
! Subprogram not used    NLFile = " "
! Subprogram not used 
! Subprogram not used    filename = trim(model)//"_stdio.nml"  ! eg. file="cpl_stdio.nml"
! Subprogram not used    inquire(file=filename,exist=exists)
! Subprogram not used 
! Subprogram not used    if (.not. exists) then
! Subprogram not used       if (s_loglev > 0) write(s_logunit,F00) "file ",trim(filename),& 
! Subprogram not used       & " doesn't exist, can not read stdio namelist from it"
! Subprogram not used       rcode = 9
! Subprogram not used    else
! Subprogram not used       unit = shr_file_getUnit()
! Subprogram not used       open (unit,file=filename,action="READ")
! Subprogram not used       read (unit,nml=stdio,iostat=rcode)
! Subprogram not used       close(unit)
! Subprogram not used       call shr_file_freeUnit( unit )
! Subprogram not used       if (rcode /= 0) then
! Subprogram not used          write(s_logunit,F01) 'ERROR: reading ',trim(filename),': iostat=',rcode
! Subprogram not used          call shr_sys_abort(subName//" ERROR reading "//trim(filename) )
! Subprogram not used       end if
! Subprogram not used    endif
! Subprogram not used    if ( len_trim(NLFile) > 0 .and. trim(stdin) /= "nochange" )then
! Subprogram not used       write(s_logunit,F00) "Error: input namelist:"
! Subprogram not used       write(s_logunit,nml=stdio)
! Subprogram not used       call shr_sys_abort(subName//" ERROR trying to both redirect AND "// &
! Subprogram not used                          "open namelist filename" )
! Subprogram not used    end if
! Subprogram not used    if ( present(NLFileOut) ) NLFileOut = NLFile
! Subprogram not used    if ( present(dirOut)    ) dirOut    = dir
! Subprogram not used    if ( present(stdinOut)  ) stdinOut  = stdin
! Subprogram not used    if ( present(stdoutOut) ) stdoutOut = stdout
! Subprogram not used    if ( present(rcodeOut)  ) rcodeOut  = rcode
! Subprogram not used  
! Subprogram not used END SUBROUTINE shr_file_stdioReadNL

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_setIO -- read in stdio namelist
!
! !DESCRIPTION:
!   This opens a namelist file specified as an argument and then opens
!   a log file associated with the unit argument.  This may be extended
!   in the future.
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_setIO( nmlfile, funit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=*)    ,intent(in)  :: nmlfile  ! namelist filename
   integer(SHR_KIND_IN),intent(in)  :: funit    ! unit number for log file

!EOP

   !--- local ---
   logical                :: exists   ! true if file exists
   character(SHR_KIND_CL) :: diri     ! directory to cd to
   character(SHR_KIND_CL) :: diro     ! directory to cd to
   character(SHR_KIND_CL) :: logfile  ! open unit 6 to this file
   integer(SHR_KIND_IN)   :: unit     ! unit number
   integer(SHR_KIND_IN)   :: rcode    ! error code

   namelist / modelio / diri,diro,logfile

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_setIO) '
   character(*),parameter :: F00   = "('(shr_file_setIO) ',4a)"
   character(*),parameter :: F01   = "('(shr_file_setIO) ',2a,i6)"

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

   diri = "."
   diro = "."
   logfile = ""

   inquire(file=nmlfile,exist=exists)

   if (.not. exists) then
      if (s_loglev > 0) write(s_logunit,F00) "file ",trim(nmlfile)," nonexistent"
      return
   else
      unit = shr_file_getUnit()
      open (unit,file=nmlfile,action="READ")
      read (unit,nml=modelio,iostat=rcode)
      close(unit)
      call shr_file_freeUnit( unit )
      if (rcode /= 0) then
         write(s_logunit,F01) 'ERROR: reading ',trim(nmlfile),': iostat=',rcode
         call shr_sys_abort(subName//" ERROR reading "//trim(nmlfile) )
      end if
   endif

   if (len_trim(logfile) > 0) then
      open(funit,file=trim(diro)//"/"//trim(logfile))
   else
      if (s_loglev > 0) write(s_logunit,F00) "logfile not opened"
   endif

END SUBROUTINE shr_file_setIO

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_setLogUnit -- Set the Log I/O Unit number
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_setLogUnit(unit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: unit     ! new unit number

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_setLogUnit) '
   character(*),parameter :: F00   = "('(shr_file_setLogUnit) ',4a)"

!-------------------------------------------------------------------------------
! Notes: Caller must be sure it's a valid unit number
!-------------------------------------------------------------------------------

   if (s_loglev > 1 .and. s_logunit-unit /= 0) then
      write(s_logunit,*) subName,': reset log unit number from/to ',s_logunit, unit
      write(     unit,*) subName,': reset log unit number from/to ',s_logunit, unit
   endif

   s_logunit = unit
 
END SUBROUTINE shr_file_setLogUnit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_setLogLevel -- Set the Log I/O Unit number
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_setLogLevel(newlevel)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: newlevel     ! new log level

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_setLogLevel) '
   character(*),parameter :: F00   = "('(shr_file_setLogLevel) ',4a)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

   if (s_loglev+newlevel > 2 .and. s_loglev-newlevel /= 0) &
      write(s_logunit,*) subName,': reset log level from/to ',s_loglev, newlevel

   s_loglev = newlevel
 
END SUBROUTINE shr_file_setLogLevel

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_getLogUnit -- Set the Log I/O Unit number
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_getLogUnit(unit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(out) :: unit     ! new unit number

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_getLogUnit) '
   character(*),parameter :: F00   = "('(shr_file_getLogUnit) ',4a)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

   unit = s_logunit
 
END SUBROUTINE shr_file_getLogUnit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_getLogLevel -- Set the Log I/O Unit number
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_getLogLevel(curlevel)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(out) :: curlevel     ! new log level

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_getLogLevel) '
   character(*),parameter :: F00   = "('(shr_file_getLogLevel) ',4a)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

   curlevel = s_loglev
 
END SUBROUTINE shr_file_getLogLevel

!===============================================================================
!===============================================================================

END MODULE shr_file_mod
