!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_die - die with mpout flushed
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_die
      use m_mpif90, only : MP_perr
      implicit none
      private	! except

      public :: die		! signal an exception
      public :: diex		! a special die() supporting macros
      public :: perr,warn	! message(s) to stderr
      public :: perr_die	! to be phased out
      public :: MP_die		! a special die() for MPI errors
      public :: MP_perr		! perr for MPI errors, from m_mpif90
      public :: MP_perr_die	! a special die() for MPI errors
      public :: assert_		! used by ASSERT() macro of assert.H

      interface die; module procedure	&
	die0_,	& ! die(where)
	die1_,	& ! die(where,message)
	die2_,	& ! die(where,proc,ier)
	die4_	  ! die(where,mesg1,ival1,mesg2,ival2)
      end interface

      interface diex; module procedure	&
	diex_	  ! diex(where,filename,lineno)
      end interface

      interface perr; module procedure	&
	perr1_,	& ! perr(where,message)
	perr2_,	& ! perr(where,proc,ier)
	perr4_	  ! perr(where,mesg1,ival1,mesg2,ival2)
      end interface
      interface warn; module procedure	&
	perr1_,	& ! perr(where,message)
	perr2_,	& ! perr(where,proc,ier)
	perr4_	  ! perr(where,mesg1,ival1,mesg2,ival2)
      end interface

      interface perr_die; module procedure	&
	die2_	  ! perr_die(where,proc,ier)
      end interface

      interface MP_die; module procedure	&
	MPdie2_	  ! MP_die(where,proc,ier)
      end interface
      interface MP_perr_die; module procedure	&
	MPdie2_	  ! MP_die(where,proc,ier)
      end interface


! !REVISION HISTORY:
! 	26Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP
!_______________________________________________________________________
  character(len=*),parameter :: myname='MCT(MPEU)::m_die'
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: die0_ - flush(mpout) before die()
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine die0_(where)
! Subprogram not used       use m_mpout, only : mpout,mpout_flush,mpout_close,mpout_ison
! Subprogram not used       use m_flow, only : flow_flush
! Subprogram not used       use m_dropdead, only : ddie => die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	26Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::die0_'
! Subprogram not used 
! Subprogram not used   call mpout_flush()
! Subprogram not used   if(mpout_ison()) call flow_flush(mpout)
! Subprogram not used   call mpout_close()
! Subprogram not used   call ddie(where)
! Subprogram not used 
! Subprogram not used end subroutine die0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: die1_ - flush(mpout) before die()
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine die1_(where,message)
! Subprogram not used       use m_mpout, only : mpout,mpout_flush,mpout_close,mpout_ison
! Subprogram not used       use m_flow, only : flow_flush
! Subprogram not used       use m_dropdead, only : ddie => die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: message
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	26Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::die1_'
! Subprogram not used 
! Subprogram not used   call mpout_flush()
! Subprogram not used   if(mpout_ison()) call flow_flush(mpout)
! Subprogram not used   call mpout_close()
! Subprogram not used 
! Subprogram not used   call perr1_(where,message)
! Subprogram not used   call ddie(where)
! Subprogram not used 
! Subprogram not used end subroutine die1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: die2_ - flush(mpout) before die()
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine die2_(where,proc,ier)
! Subprogram not used       use m_mpout, only : mpout,mpout_flush,mpout_close,mpout_ison
! Subprogram not used       use m_flow, only : flow_flush
! Subprogram not used       use m_dropdead, only : ddie => die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: proc
! Subprogram not used       integer,intent(in) :: ier
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	26Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::die2_'
! Subprogram not used 
! Subprogram not used   call mpout_flush()
! Subprogram not used   if(mpout_ison()) call flow_flush(mpout)
! Subprogram not used   call mpout_close()
! Subprogram not used 
! Subprogram not used   call perr2_(where,proc,ier)
! Subprogram not used   call ddie(where)
! Subprogram not used 
! Subprogram not used end subroutine die2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: die4_ - flush(mpout) before die()
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine die4_(where,mesg1,ival1,mesg2,ival2)
! Subprogram not used       use m_mpout, only : mpout,mpout_flush,mpout_close,mpout_ison
! Subprogram not used       use m_flow, only : flow_flush
! Subprogram not used       use m_dropdead, only : ddie => die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: mesg1
! Subprogram not used       integer,intent(in) :: ival1
! Subprogram not used       character(len=*),intent(in) :: mesg2
! Subprogram not used       integer,intent(in) :: ival2
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	26Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::die4_'
! Subprogram not used 
! Subprogram not used   call mpout_flush()
! Subprogram not used   if(mpout_ison()) call flow_flush(mpout)
! Subprogram not used   call mpout_close()
! Subprogram not used 
! Subprogram not used   call perr4_(where,mesg1,ival1,mesg2,ival2)
! Subprogram not used   call ddie(where)
! Subprogram not used 
! Subprogram not used end subroutine die4_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: diex_ - flush(mpout) before die()
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine diex_(where,filename,line)
! Subprogram not used       use m_mpout, only : mpout,mpout_flush,mpout_close,mpout_ison
! Subprogram not used       use m_flow, only : flow_flush
! Subprogram not used       use m_dropdead, only : ddie => die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: filename
! Subprogram not used       integer,intent(in) :: line
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	26Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::diex_'
! Subprogram not used 
! Subprogram not used   call mpout_flush()
! Subprogram not used   if(mpout_ison()) call flow_flush(mpout)
! Subprogram not used   call mpout_close()
! Subprogram not used   call ddie(where,filename,line)
! Subprogram not used 
! Subprogram not used end subroutine diex_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: perr1_ - send a simple error message to _stderr_
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine perr1_(where,message)
      use m_stdio,only : stderr
      implicit none
      character(len=*),intent(in) :: where
      character(len=*),intent(in) :: message

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::perr1_'

  write(stderr,'(3a)') where,': ',message

end subroutine perr1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: perr2_ - send a simple error message to _stderr_
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine perr2_(where,proc,ier)
! Subprogram not used       use m_stdio,only : stderr
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: proc
! Subprogram not used       integer,intent(in) :: ier
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::perr2_'
! Subprogram not used   character(len=16) :: cer
! Subprogram not used   integer :: ios
! Subprogram not used 
! Subprogram not used   cer='*******'
! Subprogram not used   write(cer,'(i16)',iostat=ios) ier
! Subprogram not used   write(stderr,'(5a)') where,': ',	&
! Subprogram not used 	proc,' error, stat =',trim(adjustl(cer))
! Subprogram not used 
! Subprogram not used end subroutine perr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: perr4_ - send a simple error message to _stderr_
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine perr4_(where,mesg1,ival1,mesg2,ival2)
! Subprogram not used       use m_stdio,only : stderr
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: mesg1
! Subprogram not used       integer,intent(in) :: ival1
! Subprogram not used       character(len=*),intent(in) :: mesg2
! Subprogram not used       integer,intent(in) :: ival2
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::perr4_'
! Subprogram not used   character(len=16) :: cval1,cval2
! Subprogram not used   integer :: ios
! Subprogram not used 
! Subprogram not used   cval1='*******'
! Subprogram not used   cval2='*******'
! Subprogram not used   write(cval1,'(i16)',iostat=ios) ival1
! Subprogram not used   write(cval2,'(i16)',iostat=ios) ival2
! Subprogram not used 
! Subprogram not used   write(stderr,'(10a)') where,': error, ',	&
! Subprogram not used 	mesg1,'=',trim(adjustl(cval1)),', ',	&
! Subprogram not used 	mesg2,'=',trim(adjustl(cval2)),'.'
! Subprogram not used 
! Subprogram not used end subroutine perr4_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MPdie2_ - invoke MP_perr before die_
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine MPdie2_(where,proc,ier)
! Subprogram not used       use m_mpif90, only : MP_perr
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: where
! Subprogram not used       character(len=*),intent(in) :: proc
! Subprogram not used       integer,intent(in) :: ier
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::MPdie2_'
! Subprogram not used 
! Subprogram not used   call MP_perr(where,proc,ier)
! Subprogram not used   call die0_(where)
! Subprogram not used 
! Subprogram not used end subroutine MPdie2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: assert_ - an utility called by ASSERT() macro only
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine assert_(str, file, line)
! Subprogram not used       use m_mpout,only : mpout,mpout_flush,mpout_close,mpout_ison
! Subprogram not used       use m_flow,only : flow_flush
! Subprogram not used       use m_dropdead,only : ddie => die
! Subprogram not used       implicit none
! Subprogram not used       Character(Len=*), Intent(In) :: str	! a message
! Subprogram not used       Character(Len=*), Intent(In) :: file	! a filename
! Subprogram not used       Integer, Intent(In) :: line		! a line number
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- modified
! Subprogram not used !		- included into m_die for easier module management
! Subprogram not used !	before	- Tom Clune
! Subprogram not used !		- Created for MPI PSAS implementation as a separate
! Subprogram not used !		  module
! Subprogram not used ! 	19Jan01	- J. Larson <larson@mcs.anl.gov> - removed nested 
! Subprogram not used !                 single/double/single quotes in the second argument
! Subprogram not used !                 to the call to perr1_().  This was done for the pgf90
! Subprogram not used !                 port.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_='ASSERT_'
! Subprogram not used 
! Subprogram not used   call mpout_flush()
! Subprogram not used   if(mpout_ison()) call flow_flush(mpout)
! Subprogram not used   call mpout_close()
! Subprogram not used 
! Subprogram not used   call perr1_(myname_,'failed: "//str//")')
! Subprogram not used   call ddie(myname_,file,line)
! Subprogram not used 
! Subprogram not used End subroutine assert_
end module m_die
