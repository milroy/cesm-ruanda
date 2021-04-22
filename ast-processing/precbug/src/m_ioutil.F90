!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_ioutil - a F90 module for several convenient I/O functions
!
! !DESCRIPTION:
!
!	m\_ioutil is a module containing several portable interfaces for
!	some highly system dependent, but frequently used I/O functions.
!
! !INTERFACE:

	module m_ioutil
	implicit none
	private	! except

	public	:: opntext,clstext ! open/close a text file
	public	:: opnieee,clsieee ! open/close a binary sequential file
	public	:: luavail	   ! return a free logical unit
	public	:: luflush	   ! flush the buffer of a given unit
	!public	:: MX_LU

! !REVISION HISTORY:
! 	16Jul96 - J. Guo	- (to do)
! 	02Apr97 - Jing Guo <guo@eramus> - finished the coding
!	11Feb97 - Jing Guo <guo@thunder> - added luflush()
!       08Nov01  - Jace A Mogill <mogill@cray.com>  FORTRAN only defines 
!                 99 units, three units below unit 10 are often used for
!                 stdin, stdout, and stderr.  Be far more conservative
!                 and stay within FORTRAN standard.
!
!EOP
!_______________________________________________________________________

	character(len=*),parameter :: myname="MCT(MPEU)::m_ioutil"
	integer,parameter :: MX_LU=99

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: opnieee - portablly open an IEEE format file
!
! !DESCRIPTION:
!
!	Open a file in IEEE format.
!
!	IEEE format is refered as a FORTRAN "unformatted" file with
!	"sequantial" access and variable record lengths.  Under common
!	Unix, it is only a file with records packed with a leading 4-
!	byte word and a trailing 4-byte word indicating the size of
!	the record in bytes.  However, under UNICOS, it is also assumed
!	to have numerical data representations represented according to
!	the IEEE standard corresponding KIND conversions.  Under a DEC
!	machine, it means that compilations of the source code should
!	have the "-bigendian" option specified.
!
! !INTERFACE:

! Subprogram not used     subroutine opnieee(lu,fname,status,ier,recl)
! Subprogram not used       use m_stdio,only : stderr
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used       integer,         intent(in) :: lu     ! logical unit number
! Subprogram not used       character(len=*),intent(in) :: fname  ! filename to be opended
! Subprogram not used       character(len=*),intent(in) :: status ! the value for STATUS=
! Subprogram not used       integer,         intent(out):: ier    ! the status
! Subprogram not used       integer,optional,intent(in) :: recl   ! record length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !	02Feb95 - Jing G. - First version included in PSAS.  It is not
! Subprogram not used !		used in the libpsas.a calls, since no binary data input/
! Subprogram not used !		output is to be handled.
! Subprogram not used !
! Subprogram not used ! 	09Oct96 - J. Guo  - Check for any previous assign() call under
! Subprogram not used !		UNICOS.
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 		! local parameter
! Subprogram not used 	character(len=*),parameter :: myname_=myname//'::opnieee'
! Subprogram not used 
! Subprogram not used 	integer,parameter :: iA=ichar('a')
! Subprogram not used 	integer,parameter :: mA=ichar('A')
! Subprogram not used 	integer,parameter :: iZ=ichar('z')
! Subprogram not used 
! Subprogram not used 	logical :: direct
! Subprogram not used 	character(len=16) :: clen
! Subprogram not used 	character(len=len(status)) :: Ustat
! Subprogram not used 	integer :: i,ic
! Subprogram not used 
! Subprogram not used ! Work-around for absoft 9.0 f90, which has trouble understanding that
! Subprogram not used ! ier is an output argument from the write() call below.
! Subprogram not used 
! Subprogram not used 	ier = 0
! Subprogram not used 
! Subprogram not used 	direct=.false.
! Subprogram not used 	if(present(recl)) then
! Subprogram not used 	  if(recl<0) then
! Subprogram not used 	    clen='****************'
! Subprogram not used 	    write(clen,'(i16)',iostat=ier) recl
! Subprogram not used 	    write(stderr,'(3a)') myname_,	&
! Subprogram not used 		': invalid recl, ',trim(adjustl(clen))
! Subprogram not used 	    ier=-1
! Subprogram not used 	    return
! Subprogram not used 	  endif
! Subprogram not used 	  direct = recl>0
! Subprogram not used 	endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 	do i=1,len(status)
! Subprogram not used 	  ic=ichar(status(i:i))
! Subprogram not used 	  if(ic >= iA .and. ic <= iZ) ic=ic+(mA-iA)
! Subprogram not used 	  Ustat(i:i)=char(ic)
! Subprogram not used 	end do
! Subprogram not used 
! Subprogram not used 	select case(Ustat)
! Subprogram not used 
! Subprogram not used 	case ('APPEND')
! Subprogram not used 
! Subprogram not used 	  if(direct) then
! Subprogram not used 	    write(stderr,'(2a)') myname_,		&
! Subprogram not used 		': invalid arguments, (status=="APPEND",recl>0)'
! Subprogram not used 	    ier=1
! Subprogram not used 	    return
! Subprogram not used 	  endif
! Subprogram not used 
! Subprogram not used 	  open(				&
! Subprogram not used 	    unit	=lu,		&
! Subprogram not used 	    file	=fname,		&
! Subprogram not used 	    form	='unformatted',	&
! Subprogram not used 	    access	='sequential',	&
! Subprogram not used 	    status	='unknown',	&
! Subprogram not used 	    position	='append',	&
! Subprogram not used 	    iostat	=ier		)
! Subprogram not used 
! Subprogram not used 	case default
! Subprogram not used 
! Subprogram not used 	  if(direct) then
! Subprogram not used 	    open(			&
! Subprogram not used 	      unit	=lu,		&
! Subprogram not used 	      file	=fname,		&
! Subprogram not used 	      form	='unformatted',	&
! Subprogram not used 	      access	='direct',	&
! Subprogram not used 	      status	=status,	&
! Subprogram not used 	      recl	=recl,		&
! Subprogram not used 	      iostat	=ier		)
! Subprogram not used 
! Subprogram not used 	  else
! Subprogram not used 	    open(			&
! Subprogram not used 	      unit	=lu,		&
! Subprogram not used 	      file	=fname,		&
! Subprogram not used 	      form	='unformatted',	&
! Subprogram not used 	      access	='sequential',	&
! Subprogram not used 	      status	=status,	&
! Subprogram not used 	      position	='asis',	&
! Subprogram not used 	      iostat	=ier		)
! Subprogram not used 	  endif
! Subprogram not used 
! Subprogram not used 	end select
! Subprogram not used 
! Subprogram not used 	end subroutine opnieee
!-----------------------------------------------------------------------
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clsieee - Close a logical unit opened by opnieee()
!
! !DESCRIPTION:
!
!	The reason for a paired clsieee() for opnieee() instead of a
!	simple close(), is for the portability reason.  For example,
!	under UNICOS, special system calls may be need to set up the
!	unit right, and the status of the unit should be restored upon
!	close.
!
! !INTERFACE:

! Subprogram not used 	subroutine clsieee(lu,ier)
! Subprogram not used 	  implicit none
! Subprogram not used 	  integer, intent(in)  :: lu	! the unit used by opnieee()
! Subprogram not used 	  integer, intent(out) :: ier	! the status
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	10Oct96 - J. Guo	- (to do)
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used 	  close(lu,iostat=ier)
! Subprogram not used 
! Subprogram not used 	end subroutine clsieee

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: opntext - portablly open a text file
!
! !DESCRIPTION:
!
!	Open a text (ASCII) file.  Under FORTRAN, it is defined as
!	"formatted" with "sequential" access.
!
! !INTERFACE:

    subroutine opntext(lu,fname,status,ier)
      implicit none

      integer,         intent(in) :: lu     ! logical unit number
      character(len=*),intent(in) :: fname  ! filename to be opended
      character(len=*),intent(in) :: status ! the value for STATUS=<>
      integer,         intent(out):: ier    ! the status


! !REVISION HISTORY:
!
!	02Feb95 - Jing G. - First version included in PSAS and libpsas.a
! 	09Oct96 - J. Guo  - modified to allow assign() call under UNICOS
!			  = and now, it is a module in Fortran 90.
!EOP
!_______________________________________________________________________

		! local parameter
	character(len=*),parameter :: myname_=myname//'::opntext'

	integer,parameter :: iA=ichar('a')
	integer,parameter :: mA=ichar('A')
	integer,parameter :: iZ=ichar('z')

	character(len=len(status)) :: Ustat
	integer :: i,ic


	do i=1,len(status)
	  ic=ichar(status(i:i))
	  if(ic >= iA .and. ic <= iZ) ic=ic+(mA-iA)
	  Ustat(i:i)=char(ic)
	end do

	select case(Ustat)

	case ('APPEND')

	  open(				&
	    unit	=lu,		&
	    file	=fname,		&
	    form	='formatted',	&
	    access	='sequential',	&
	    status	='unknown',	&
	    position	='append',	&
	    iostat	=ier		)

	case default

	  open(				&
	    unit	=lu,		&
	    file	=fname,		&
	    form	='formatted',	&
	    access	='sequential',	&
	    status	=status,	&
	    position	='asis',	&
	    iostat	=ier		)

	end select

	end subroutine opntext

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clstext - close a text file opend with an opntext() call
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clstext(lu,ier)
      implicit none

      integer, intent(in)  :: lu  ! a logical unit to close
      integer, intent(out) :: ier ! the status

! !REVISION HISTORY:
! 	09Oct96 - J. Guo	- (to do)
!EOP
!_______________________________________________________________________

	close(lu,iostat=ier)

	end subroutine clstext

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: luavail - locate the next available unit
!
! !DESCRIPTION:
!
!    luavail() Look for an available (not opened and not statically
!    assigned to any I/O attributes to) logical unit.
!
! !INTERFACE:

	function luavail()
	  use m_stdio
	  implicit none
	  integer :: luavail	! result

! !REVISION HISTORY:
! 	23Apr98 - Jing Guo <guo@thunder> - new prototype/prolog/code
!			- with additional unit constraints for SunOS.
!
! 	: Jing Guo, [09-Oct-96]
! 		+ Checking also Cray assign() attributes, with some
! 		  changes to the code.  See also other routines.
!
! 	: Jing Guo, [01-Apr-94]
! 		+ Initial code.
!   2001-11-08  - Jace A Mogill <mogill@cray.com>  clean up
!		  logic for finding lu.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::luavail'

	integer lu,ios
	logical inuse

	lu=10
	ios=0
	inuse=.true.

	do while(ios.eq.0 .and. inuse .and. lu.le.MX_LU)
	  lu=lu+1
	  inquire(unit=lu,opened=inuse,iostat=ios)
	end do

	if(ios.ne.0) lu=-1
	luavail=lu
end function luavail

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: luflush - a uniform interface of system flush()
!
! !DESCRIPTION:
!
!	Flush() calls available on many systems are often implementation
!	dependent.  This subroutine provides a uniform interface.  It
!	also ignores invalid logical unit value.
!
! !INTERFACE:

! Subprogram not used     subroutine luflush(unit)
! Subprogram not used       use m_stdio, only : stdout
! Subprogram not used       implicit none
! Subprogram not used       integer,optional,intent(in) :: unit
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	13Mar98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !       08Jul02 - E. Ong <eong@mcs.anl.gov> - added flush support for nag95
! Subprogram not used !  2001-11-08  Jace A Mogill <mogill@cray.com>  - Flush is not part of
! Subprogram not used !              the F90 standard.  Default is NO unit flush.
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::luflush'
! Subprogram not used 
! Subprogram not used   integer :: ier
! Subprogram not used   integer :: lu
! Subprogram not used 
! Subprogram not used 	! Which logical unit number?
! Subprogram not used 
! Subprogram not used   lu=stdout
! Subprogram not used   if(present(unit)) lu=unit
! Subprogram not used   if(lu < 0) return
! Subprogram not used 
! Subprogram not used 	! The following call may be system dependent.
! Subprogram not used 
! Subprogram not used   call flush(lu)
! Subprogram not used 
! Subprogram not used end subroutine luflush
!-----------------------------------------------------------------------
end module m_ioutil
!.
