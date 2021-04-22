!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_realkinds - real KIND definitions
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_realkinds
      implicit none
      private	! except

      public :: kind_r4		! real*4
      public :: kind_r8		! real*8
      public :: kind_r		! default real
      public :: SP		! default REAL
      public :: DP		! default DOUBLE_PRECISION
      public :: FP		! general floating point precision

      real*4,parameter :: mpeuR4=1.
      real*8,parameter :: mpeuR8=1.
      real,  parameter :: mpeuR =1.

      integer,parameter :: SP = kind(1.  )
      integer,parameter :: DP = kind(1.D0)

!     Set the current default floating point precision 
      integer,parameter :: FP = DP

      integer,parameter :: kind_r4=kind(mpeuR4)
      integer,parameter :: kind_r8=kind(mpeuR8)
      integer,parameter :: kind_r =kind(mpeuR )

! !REVISION HISTORY:
! 	19Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	23Jan03 - R. Jacob <jacob@mcs.anl.gov> - add FP
!EOP
!_______________________________________________________________________
  character(len=*),parameter :: myname='MCT(MPEU)::m_realkinds'

end module m_realkinds
