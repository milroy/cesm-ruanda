!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Permuter - permute/unpermute
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Permuter
      implicit none
      private	! except

      public :: permute	
      public :: unpermute

    interface permute; module procedure	&
	permutei_,	&	! integer in place
	permuteio_,	&	! integer with an output
	permutei1_,	&	! integer in place
	permuteio1_,	&	! integer with an output
	permuter_,	&	! real in place
	permutero_,	&	! real with an output
	permuter1_,	&	! real in place
	permutero1_,	&	! real with an output
	permuted_,	&	! dble in place
	permutedo_,	&	! dble with an output
	permuted1_,	&	! dble in place
	permutedo1_,	&	! dble with an output
	permutel_,	&	! logical in place
	permutelo_,	&	! logical with an output
	permutel1_,	&	! logical in place
	permutelo1_		! logical with an output
    end interface

    interface unpermute; module procedure	&
	unpermutei_,	&	! integer in place
	unpermuteio_,	&	! integer with an output
	unpermutei1_,	&	! integer in place
	unpermuteio1_,	&	! integer with an output
	unpermuter_,	&	! real in place
	unpermutero_,	&	! real with an output
	unpermuter1_,	&	! real in place
	unpermutero1_,	&	! real with an output
	unpermuted_,	&	! dble in place
	unpermutedo_,	&	! dble with an output
	unpermuted1_,	&	! dble in place
	unpermutedo1_,	&	! dble with an output
	unpermutel_,	&	! logical in place
	unpermutelo_,	&	! logical with an output
	unpermutel1_,	&	! logical in place
	unpermutelo1_		! logical with an output
    end interface

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT(MPEU)::m_Permuter'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutei_ - permute an integer array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine permutei_(ary,indx,n)
      use m_die
      implicit none
      integer,dimension(:),intent(inout) :: ary
      integer,dimension(:),intent(in)    :: indx
      integer,             intent(in)    :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::permutei_'

  integer,allocatable,dimension(:) :: wk
  integer :: i,ier

  allocate(wk(n),stat=ier)
	if(ier/=0) call perr_die(myname_,'allocate()',ier)

  call permuteio_(wk,ary,indx,n)

  do i=1,n
    ary(i)=wk(i)
  end do

  deallocate(wk,stat=ier)
	if(ier/=0) call perr_die(myname_,'deallocate()',ier)

end subroutine permutei_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permuteio_ - permute an integer array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine permuteio_(aout,ary,indx,n)
      implicit none
      integer,dimension(:),intent(inout) :: aout
      integer,dimension(:),intent(in ) :: ary
      integer,dimension(:),intent(in)  :: indx
      integer,             intent(in)  :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::permuteio_'

  integer :: i,l

  do i=1,n
    l=indx(i)
    aout(i)=ary(l)
  end do

end subroutine permuteio_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutei_ - unpermute a _permuted_ integer array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpermutei_(ary,indx,n)
      use m_die
      implicit none
      integer,dimension(:),intent(inout) :: ary
      integer,dimension(:),intent(in)    :: indx
      integer,             intent(in)    :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpermutei_'

  integer,allocatable,dimension(:) :: wk
  integer :: i,ier

  allocate(wk(n),stat=ier)
	if(ier/=0) call perr_die(myname_,'allocate()',ier)

  call unpermuteio_(wk,ary,indx,n)

  do i=1,n
    ary(i)=wk(i)
  end do

  deallocate(wk,stat=ier)
	if(ier/=0) call perr_die(myname_,'deallocate()',ier)

end subroutine unpermutei_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermuteio_ - unpermute a _permuted_ integer array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpermuteio_(aout,ary,indx,n)
      implicit none
      integer,dimension(:),intent(inout) :: aout
      integer,dimension(:),intent(in)  :: ary
      integer,dimension(:),intent(in)  :: indx
      integer,             intent(in)  :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpermuteio_'

  integer :: i,l

  do i=1,n
    l=indx(i)
    aout(l)=ary(i)
  end do

end subroutine unpermuteio_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permuter_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permuter_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(inout) :: ary
! Subprogram not used       integer ,dimension(:),intent(in)    :: indx
! Subprogram not used       integer ,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permuter_'
! Subprogram not used 
! Subprogram not used   real(kind(ary)),allocatable,dimension(:) :: wk
! Subprogram not used   integer :: i,ier
! Subprogram not used 
! Subprogram not used   allocate(wk(n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call permutero_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(i)=wk(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine permuter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutero_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutero_(aout,ary,indx,n)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(inout) :: aout
! Subprogram not used       real(SP),dimension(:),intent(in)  :: ary
! Subprogram not used       integer ,dimension(:),intent(in)  :: indx
! Subprogram not used       integer ,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutero_'
! Subprogram not used 
! Subprogram not used   integer :: i,l
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(i)=ary(l)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine permutero_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermuter_ - unpermute a _permuted_ real array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermuter_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(inout) :: ary
! Subprogram not used       integer ,dimension(:),intent(in)    :: indx
! Subprogram not used       integer ,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermuter_'
! Subprogram not used 
! Subprogram not used   real(kind(ary)),allocatable,dimension(:) :: wk
! Subprogram not used   integer :: i,ier
! Subprogram not used 
! Subprogram not used   allocate(wk(n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call unpermutero_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(i)=wk(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine unpermuter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutero_ - unpermute a _permuted_ real array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutero_(aout,ary,indx,n)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(inout) :: aout
! Subprogram not used       real(SP),dimension(:),intent(in)  :: ary
! Subprogram not used       integer ,dimension(:),intent(in)  :: indx
! Subprogram not used       integer ,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutero_'
! Subprogram not used 
! Subprogram not used   integer :: i,l
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(l)=ary(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine unpermutero_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permuted_ - permute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine permuted_(ary,indx,n)
      use m_die
      use m_realkinds,only : DP
      implicit none
      real(DP),dimension(:),intent(inout) :: ary
      integer ,dimension(:),intent(in)    :: indx
      integer ,             intent(in)    :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::permuted_'

  real(kind(ary)),allocatable,dimension(:) :: wk
  integer :: i,ier

  allocate(wk(n),stat=ier)
	if(ier/=0) call perr_die(myname_,'allocate()',ier)

  call permutedo_(wk,ary,indx,n)

  do i=1,n
    ary(i)=wk(i)
  end do

  deallocate(wk,stat=ier)
	if(ier/=0) call perr_die(myname_,'deallocate()',ier)

end subroutine permuted_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutedo_ - permute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine permutedo_(aout,ary,indx,n)
      use m_realkinds,only : DP
      implicit none
      real(DP),dimension(:),intent(inout) :: aout
      real(DP),dimension(:),intent(in)  :: ary
      integer ,dimension(:),intent(in)  :: indx
      integer ,             intent(in)  :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::permutedo_'

  integer :: i,l

  do i=1,n
    l=indx(i)
    aout(i)=ary(l)
  end do

end subroutine permutedo_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermuted_ - unpermute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpermuted_(ary,indx,n)
      use m_die
      use m_realkinds,only : DP
      implicit none
      real(DP),dimension(:),intent(inout) :: ary
      integer ,dimension(:),intent(in)    :: indx
      integer ,             intent(in)    :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpermuted_'

  real(kind(ary)),allocatable,dimension(:) :: wk
  integer :: i,ier

  allocate(wk(n),stat=ier)
	if(ier/=0) call perr_die(myname_,'allocate()',ier)

  call unpermutedo_(wk,ary,indx,n)

  do i=1,n
    ary(i)=wk(i)
  end do

  deallocate(wk,stat=ier)
	if(ier/=0) call perr_die(myname_,'deallocate()',ier)

end subroutine unpermuted_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutedo_ - unpermute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpermutedo_(aout,ary,indx,n)
      use m_realkinds,only : DP
      implicit none
      real(DP),dimension(:),intent(inout) :: aout
      real(DP),dimension(:),intent(in)  :: ary
      integer ,dimension(:),intent(in)  :: indx
      integer ,             intent(in)  :: n

! !REVISION HISTORY:
! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpermutedo_'

  integer :: i,l

  do i=1,n
    l=indx(i)
    aout(l)=ary(i)
  end do

end subroutine unpermutedo_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutel_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutel_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(inout) :: ary
! Subprogram not used       integer,dimension(:),intent(in)    :: indx
! Subprogram not used       integer,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutel_'
! Subprogram not used 
! Subprogram not used   logical,allocatable,dimension(:) :: wk
! Subprogram not used   integer :: i,ier
! Subprogram not used 
! Subprogram not used   allocate(wk(n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call permutelo_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(i)=wk(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine permutel_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutelo_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutelo_(aout,ary,indx,n)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(inout) :: aout
! Subprogram not used       logical,dimension(:),intent(in)  :: ary
! Subprogram not used       integer,dimension(:),intent(in)  :: indx
! Subprogram not used       integer,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutelo_'
! Subprogram not used 
! Subprogram not used   integer :: i,l
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(i)=ary(l)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine permutelo_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutel_ - unpermute a _permuted_ logical array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutel_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(inout) :: ary
! Subprogram not used       integer,dimension(:),intent(in)    :: indx
! Subprogram not used       integer,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutel_'
! Subprogram not used 
! Subprogram not used   logical,allocatable,dimension(:) :: wk
! Subprogram not used   integer :: i,ier
! Subprogram not used 
! Subprogram not used   allocate(wk(n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call unpermutelo_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(i)=wk(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine unpermutel_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutelo_ - unpermute a _permuted_ logical array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutelo_(aout,ary,indx,n)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(inout) :: aout
! Subprogram not used       logical,dimension(:),intent(in)  :: ary
! Subprogram not used       integer,dimension(:),intent(in)  :: indx
! Subprogram not used       integer,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutelo_'
! Subprogram not used 
! Subprogram not used   integer :: i,l
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(l)=ary(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine unpermutelo_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutei1_ - permute an integer array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutei1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer,dimension(:),intent(in)    :: indx
! Subprogram not used       integer,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutei1_'
! Subprogram not used 
! Subprogram not used   integer,allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call permuteio1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine permutei1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permuteio1_ - permute an integer array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permuteio1_(aout,ary,indx,n)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(inout) :: aout
! Subprogram not used       integer,dimension(:,:),intent(in ) :: ary
! Subprogram not used       integer,dimension(:),intent(in)  :: indx
! Subprogram not used       integer,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permuteio1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,i)=ary(1:m,l)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine permuteio1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutei1_ - unpermute a _permuted_ integer array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutei1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer,dimension(:),intent(in)    :: indx
! Subprogram not used       integer,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutei1_'
! Subprogram not used 
! Subprogram not used   integer,allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call unpermuteio1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine unpermutei1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermuteio1_ - unpermute a _permuted_ integer array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermuteio1_(aout,ary,indx,n)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(inout) :: aout
! Subprogram not used       integer,dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer,dimension(:),intent(in)  :: indx
! Subprogram not used       integer,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermuteio1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,l)=ary(1:m,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine unpermuteio1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permuter1_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permuter1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer ,dimension(:),intent(in)    :: indx
! Subprogram not used       integer ,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permuter1_'
! Subprogram not used 
! Subprogram not used   real(kind(ary)),allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call permutero1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine permuter1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutero1_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutero1_(aout,ary,indx,n)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(inout) :: aout
! Subprogram not used       real(SP),dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer ,dimension(:),intent(in)  :: indx
! Subprogram not used       integer ,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutero1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,i)=ary(1:m,l)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine permutero1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermuter1_ - unpermute a _permuted_ real array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermuter1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer ,dimension(:),intent(in)    :: indx
! Subprogram not used       integer ,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermuter1_'
! Subprogram not used 
! Subprogram not used   real(kind(ary)),allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call unpermutero1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine unpermuter1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutero1_ - unpermute a _permuted_ real array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutero1_(aout,ary,indx,n)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(inout) :: aout
! Subprogram not used       real(SP),dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer ,dimension(:),intent(in)  :: indx
! Subprogram not used       integer ,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutero1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,l)=ary(1:m,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine unpermutero1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permuted1_ - permute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permuted1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer ,dimension(:),intent(in)    :: indx
! Subprogram not used       integer ,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permuted1_'
! Subprogram not used 
! Subprogram not used   real(kind(ary)),allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call permutedo1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine permuted1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutedo1_ - permute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutedo1_(aout,ary,indx,n)
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(inout) :: aout
! Subprogram not used       real(DP),dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer ,dimension(:),intent(in)  :: indx
! Subprogram not used       integer ,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutedo1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,i)=ary(1:m,l)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine permutedo1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermuted1_ - unpermute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermuted1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer ,dimension(:),intent(in)    :: indx
! Subprogram not used       integer ,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermuted1_'
! Subprogram not used 
! Subprogram not used   real(kind(ary)),allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call unpermutedo1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine unpermuted1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutedo1_ - unpermute a double precision array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutedo1_(aout,ary,indx,n)
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(inout) :: aout
! Subprogram not used       real(DP),dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer ,dimension(:),intent(in)  :: indx
! Subprogram not used       integer ,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutedo1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,l)=ary(1:m,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine unpermutedo1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutel1_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutel1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer,dimension(:),intent(in)    :: indx
! Subprogram not used       integer,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutel1_'
! Subprogram not used 
! Subprogram not used   logical,allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call permutelo1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine permutel1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permutelo1_ - permute a real array according to indx[]
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine permutelo1_(aout,ary,indx,n)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(inout) :: aout
! Subprogram not used       logical,dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer,dimension(:),intent(in)  :: indx
! Subprogram not used       integer,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::permutelo1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,i)=ary(1:m,l)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine permutelo1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutel1_ - unpermute a _permuted_ logical array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutel1_(ary,indx,n)
! Subprogram not used       use m_die
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(inout) :: ary
! Subprogram not used       integer,dimension(:),intent(in)    :: indx
! Subprogram not used       integer,             intent(in)    :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutel1_'
! Subprogram not used 
! Subprogram not used   logical,allocatable,dimension(:,:) :: wk
! Subprogram not used   integer :: i,l,ier
! Subprogram not used 
! Subprogram not used   l=size(ary,1)
! Subprogram not used   allocate(wk(l,n),stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used   call unpermutelo1_(wk,ary,indx,n)
! Subprogram not used 
! Subprogram not used   do i=1,n
! Subprogram not used     ary(:,i)=wk(:,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   deallocate(wk,stat=ier)
! Subprogram not used 	if(ier/=0) call perr_die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used end subroutine unpermutel1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermutelo1_ - unpermute a _permuted_ logical array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine unpermutelo1_(aout,ary,indx,n)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(inout) :: aout
! Subprogram not used       logical,dimension(:,:),intent(in)  :: ary
! Subprogram not used       integer,dimension(:),intent(in)  :: indx
! Subprogram not used       integer,             intent(in)  :: n
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::unpermutelo1_'
! Subprogram not used 
! Subprogram not used   integer :: i,l,m
! Subprogram not used 
! Subprogram not used   m=min(size(aout,1),size(ary,1))
! Subprogram not used   do i=1,n
! Subprogram not used     l=indx(i)
! Subprogram not used     aout(1:m,l)=ary(1:m,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end subroutine unpermutelo1_

end module m_Permuter
