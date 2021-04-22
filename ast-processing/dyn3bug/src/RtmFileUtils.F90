module RtmFileUtils

!-----------------------------------------------------------------------
! Module containing file I/O utilities
!
! !USES:
  use shr_sys_mod , only : shr_sys_abort
  use shr_file_mod, only : shr_file_get, shr_file_getUnit, shr_file_freeUnit
  use RtmSpmd     , only : masterproc
  use RtmVar      , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: get_filename  !Returns filename given full pathname
  public :: opnfil        !Open local unformatted or formatted file
  public :: getfil        !Obtain local copy of file
  public :: relavu        !Close and release Fortran unit no longer in use
  public :: getavu        !Get next available Fortran unit number
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !PRIVATE MEMBER FUNCTIONS: None
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

  character(len=256) function get_filename (fulpath)

    ! !DESCRIPTION:
    ! Returns filename given full pathname
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: fulpath !full pathname
    !
    ! !LOCAL VARIABLES:
    integer i     !loop index
    integer klen  !length of fulpath character string
    !----------------------------------------------------------

    klen = len_trim(fulpath)
    do i = klen, 1, -1
       if (fulpath(i:i) == '/') go to 10
    end do
    i = 0
10  get_filename = fulpath(i+1:klen)

  end function get_filename

!------------------------------------------------------------------------

   subroutine getfil (fulpath, locfn, iflag)

     ! !DESCRIPTION:
     ! Obtain local copy of file. First check current working directory,
     ! Next check full pathname[fulpath] on disk
     ! 
     ! !ARGUMENTS:
     implicit none
     character(len=*), intent(in)  :: fulpath !Archival or permanent disk full pathname
     character(len=*), intent(out) :: locfn   !output local file name
     integer,          intent(in)  :: iflag   !0=>abort if file not found 1=>do not abort

     ! !LOCAL VARIABLES:
     integer i               !loop index
     integer klen            !length of fulpath character string
     logical lexist          !true if local file exists
     !--------------------------------------------------

     ! get local file name from full name
     locfn = get_filename( fulpath )
     if (len_trim(locfn) == 0) then
	if (masterproc) write(iulog,*)'(GETFIL): local filename has zero length'
        call shr_sys_abort()
     else
        if (masterproc) write(iulog,*)'(GETFIL): attempting to find local file ',  &
             trim(locfn)
     endif

     ! first check if file is in current working directory.
     inquire (file=locfn,exist=lexist)
     if (lexist) then
        if (masterproc) write(iulog,*) '(GETFIL): using ',trim(locfn), &
             ' in current working directory'
        RETURN
     endif

     ! second check for full pathname on disk
     locfn = fulpath

     inquire (file=fulpath,exist=lexist)
     if (lexist) then
        if (masterproc) write(iulog,*) '(GETFIL): using ',trim(fulpath)
        RETURN
     else
        if (masterproc) write(iulog,*)'(GETFIL): failed getting file from full path: ', fulpath
        if (iflag==0) then
           call shr_sys_abort ('GETFIL: FAILED to get '//trim(fulpath))
        else
           RETURN
        endif
     endif

   end subroutine getfil

!------------------------------------------------------------------------

! Subprogram not used    subroutine opnfil (locfn, iun, form)
! Subprogram not used 
! Subprogram not used      ! !DESCRIPTION:
! Subprogram not used      ! Open file locfn in unformatted or formatted form on unit iun
! Subprogram not used      !
! Subprogram not used      ! !ARGUMENTS:
! Subprogram not used      implicit none
! Subprogram not used      character(len=*), intent(in):: locfn  !file name
! Subprogram not used      integer, intent(in):: iun             !fortran unit number
! Subprogram not used      character(len=1), intent(in):: form   !file format: u = unformatted,
! Subprogram not used 
! Subprogram not used      ! !LOCAL VARIABLES:
! Subprogram not used      integer ioe             !error return from fortran open
! Subprogram not used      character(len=11) ft    !format type: formatted. unformatted
! Subprogram not used      !-----------------------------------------------------------
! Subprogram not used 
! Subprogram not used      if (len_trim(locfn) == 0) then
! Subprogram not used         write(iulog,*)'(OPNFIL): local filename has zero length'
! Subprogram not used         call shr_sys_abort()
! Subprogram not used      endif
! Subprogram not used      if (form=='u' .or. form=='U') then
! Subprogram not used         ft = 'unformatted'
! Subprogram not used      else
! Subprogram not used         ft = 'formatted  '
! Subprogram not used      end if
! Subprogram not used      open (unit=iun,file=locfn,status='unknown',form=ft,iostat=ioe)
! Subprogram not used      if (ioe /= 0) then
! Subprogram not used         write(iulog,*)'(OPNFIL): failed to open file ',trim(locfn),        &
! Subprogram not used              &     ' on unit ',iun,' ierr=',ioe
! Subprogram not used         call shr_sys_abort()
! Subprogram not used      else if ( masterproc )then
! Subprogram not used         write(iulog,*)'(OPNFIL): Successfully opened file ',trim(locfn),   &
! Subprogram not used              &     ' on unit= ',iun
! Subprogram not used      end if
! Subprogram not used 
! Subprogram not used    end subroutine opnfil

!------------------------------------------------------------------------

  integer function getavu()

    ! !DESCRIPTION:
    ! Get next available Fortran unit number.
    implicit none

    getavu = shr_file_getunit()

  end function getavu

!------------------------------------------------------------------------

  subroutine relavu (iunit)

    ! !DESCRIPTION:
    ! Close and release Fortran unit no longer in use!

    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iunit    !Fortran unit number
    !----------------------------------------------------

    close(iunit)
    call shr_file_freeUnit(iunit)

  end subroutine relavu

end module RtmFileUtils
