!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS m_MCTWorld.F90,v 1.26 2007/06/01 19:56:25 rloy Exp
! CVS MCT_2_4_0 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_MCTWorld -- MCTWorld Class
!
! !DESCRIPTION:
! MCTWorld is a datatype which acts as a component model registry.
! All models communicating through MCT must participate in initialization
! of MCTWorld.  The single instance of MCTWorld, {\tt ThisMCTWorld} stores
! the component id and local and global processor rank of each component.
! This module contains methods for creating and destroying {\tt ThisMCTWorld}
! as well as inquiry functions.
!
! !INTERFACE:

 module m_MCTWorld
!
! !USES:
      use m_List, only : List   ! Support for List components.

      implicit none

      private   ! except

! !PUBLIC TYPES:

      public :: MCTWorld        ! The MCTWorld  class data structure

    type MCTWorld
      integer :: MCT_comm                          ! MCT communicator
      integer :: ncomps	                           ! Total number of components
      integer :: mygrank                           ! Rank of this processor in 
                                                   ! global communicator.
      integer,dimension(:),pointer :: nprocspid	   ! Number of processes 
                                                   ! each component is on (e.g. rank of its
						   ! local communicator.
      integer,dimension(:,:),pointer :: idGprocid  ! Translate between local component rank
                                                   ! rank in global communicator.
						   ! idGprocid(modelid,localrank)=globalrank
    end type MCTWorld

! !PUBLIC DATA MEMBERS:

    type(MCTWorld) :: ThisMCTWorld   !  declare the MCTWorld

! !PUBLIC MEMBER FUNCTIONS:
      public :: initialized          ! Determine if MCT is initialized
      public :: init                 ! Create a MCTWorld
      public :: clean                ! Destroy a MCTWorld
      public :: NumComponents        ! Number of Components in the MCTWorld
      public :: ComponentNumProcs    ! Number of processes owned by a given
                                     ! component
      public :: ComponentToWorldRank ! Given the rank of a process on a 
                                     ! component, return its rank on the 
                                     ! world communicator
      public :: ComponentRootRank    ! Return the rank on the world 
                                     ! communicator of the root process of 
                                     ! a component
      public :: ThisMCTWorld         ! Instantiation of the MCTWorld

!  

    interface initialized ; module procedure &
      initialized_
    end interface
    interface init ; module procedure &
      initd_, &
      initm_, &
      initr_
    end interface
    interface clean ; module procedure clean_ ; end interface
    interface NumComponents ; module procedure &
       NumComponents_ 
    end interface
    interface ComponentNumProcs ; module procedure &
       ComponentNumProcs_ 
    end interface
    interface ComponentToWorldRank ; module procedure &
       ComponentToWorldRank_ 
    end interface
    interface ComponentRootRank ; module procedure &
       ComponentRootRank_ 
    end interface



! !REVISION HISTORY:
! 19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! 05Feb01 - J. Larson <larson@mcs.anl.gov> - added query and
!           local-to-global mapping services NumComponents, 
!           ComponentNumProcs, ComponentToWorldRank, and ComponentRootRank
! 08Feb01 - R. Jacob <jacob@mcs.anl.gov> - add mylrank and mygrank
!           to datatype
! 20Apr01 - R. Jacob <jacob@mcs.anl.gov> - remove allids from
!           MCTWorld datatype.  Not needed because component
!           ids are always from 1 to number-of-components.
! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - remove myid, mynprocs
!           and mylrank from MCTWorld datatype because they are not 
!           clearly defined in PCM mode.  Add MCT_comm for future use.
! 03Aug01 - E. Ong <eong@mcs.anl.gov> - explicity specify starting
!           address in mpi_irecv
! 27Nov01 - E. Ong <eong@mcs.anl.gov> - added R. Jacob's version of initd_
!           to support PCM mode. 
! 15Feb02 - R. Jacob - elminate use of MP_COMM_WORLD.  Use
!           argument globalcomm instead.  Create MCT_comm from
!           globalcomm
!EOP __________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_MCTWorld'

 contains



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initialized_ - determine if MCTWorld is initialized
!
! !DESCRIPTION:
! This routine may be used to determine whether {\tt MCTWorld::init} 
! has been called.  If not, the user must call {\tt init} before
! performing any other MCT library calls.
!
! !INTERFACE:

! Subprogram not used  logical function initialized_()
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 01June07 - R. Loy <rloy@mcs.anl.gov> - initial version
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used 
! Subprogram not used   initialized_ = associated(ThisMCTWorld%nprocspid)
! Subprogram not used 
! Subprogram not used   end function initialized_




!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initm_ - initialize MCTWorld
!
! !DESCRIPTION:
! Do a distributed init of MCTWorld for the case where a set of processors
! contains more then one model and the models may not span the set of processors.
! {\tt ncomps} is the total number of components in the entire coupled system.
! {\tt globalcomm} encompasses all the models (typically this can be MPI\_COMM\_WORLD).
! {\tt mycomms} is an array of MPI communicators, each sized for the appropriate model
! and {\tt myids} is a corresponding array of integers containing the model ids for 
! the models on this particular set of processors.
!
! This routine is called once for the models covered by the set of processors.
!
! !INTERFACE:

 subroutine initm_(ncomps,globalcomm,mycomms,myids)
!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)	       :: ncomps          ! number of components
      integer, intent(in)	       :: globalcomm      ! global communicator
      integer, dimension(:),pointer    :: mycomms         ! my communicators
      integer, dimension(:),pointer    :: myids           ! component ids

! !REVISION HISTORY:
! 20Sep07 - T. Craig migrated code from initd routine
! 20Sep07 - T. Craig - made mycomms an array
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initm_'
  integer :: ier,myGid,myLid,i,mysize,Gsize,j

! arrays allocated on the root to coordinate gathring of data
! and non-blocking receives by the root
  integer, dimension(:), allocatable :: compids,reqs,nprocs,Gprocids
  integer, dimension(:), allocatable :: root_nprocs
  integer, dimension(:,:),allocatable :: status,root_idGprocid
  integer, dimension(:,:),pointer :: tmparray
  integer,dimension(:),pointer :: apoint
! ------------------------------------------------------------------

! Check that ncomps is a legal value
  if(ncomps < 1) then
     call die(myname_, "argument ncomps can't less than one!",ncomps)
  endif

  if (size(myids) /= size(mycomms)) then
     call die(myname_, "size of myids and mycomms inconsistent")
  endif

! make sure this has not been called already
  if(associated(ThisMCTWorld%nprocspid) ) then
     write(stderr,'(2a)') myname_, &
      'MCTERROR:  MCTWorld has already been initialized...Continuing'
       RETURN
  endif

! determine overall size
  call MP_comm_size(globalcomm,Gsize,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

! determine my rank in comm_world
  call MP_comm_rank(globalcomm,myGid,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

! allocate space on global root to receive info about 
! the other components
  if(myGid == 0) then
     allocate(nprocs(ncomps),compids(ncomps),&
     reqs(ncomps),status(MP_STATUS_SIZE,ncomps),&
     root_nprocs(ncomps),stat=ier)
     if (ier /= 0) then
        call die(myname_, 'allocate(nprocs,...)',ier)
     endif
  endif


!!!!!!!!!!!!!!!!!!
!  Gather the number of procs from the root of each component
!!!!!!!!!!!!!!!!!!
!
!  First on the global root, post a receive for each component
  if(myGid == 0) then
    do i=1,ncomps
       call MPI_IRECV(root_nprocs(i), 1, MP_INTEGER, MP_ANY_SOURCE,i, &
	 globalcomm, reqs(i), ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(root_nprocs)',ier)
    enddo
  endif

!  The local root on each component sends
  do i=1,size(myids)
    if(mycomms(i)/=MP_COMM_NULL) then
      call MP_comm_size(mycomms(i),mysize,ier)
      if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)
      call MP_comm_rank(mycomms(i),myLid,ier)
      if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)
      if(myLid == 0) then
        call MPI_SEND(mysize,1,MP_INTEGER,0,myids(i),globalcomm,ier)
        if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(mysize)',ier)
      endif
    endif
  enddo

!  Global root waits for all sends
  if(myGid == 0) then
    call MPI_WAITALL(size(reqs), reqs, status, ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif
! Global root now knows how many processors each component is using

!!!!!!!!!!!!!!!!!!
! end of nprocs
!!!!!!!!!!!!!!!!!!


! allocate a tmp array for the receive on root.
  if(myGid == 0) then
    allocate(tmparray(0:Gsize-1,ncomps),stat=ier)
    if(ier/=0) call die(myname_,'allocate(tmparray)',ier)

! fill tmparray with a bad rank value for later error checking
    tmparray = -1
  endif

!!!!!!!!!!!!!!!!!!
!  Gather the Gprocids from each local root
!!!!!!!!!!!!!!!!!!
!
!  First on the global root, post a receive for each component
  if(myGid == 0) then
    do i=1,ncomps
       apoint => tmparray(0:Gsize-1,i)
       call MPI_IRECV(apoint(1), root_nprocs(i),MP_INTEGER, &
       MP_ANY_SOURCE,i,globalcomm, reqs(i), ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV()',ier)
    enddo
  endif

!  The root on each component sends
  do i=1,size(myids)
    if(mycomms(i)/=MP_COMM_NULL) then
      call MP_comm_size(mycomms(i),mysize,ier)
      if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)
      call MP_comm_rank(mycomms(i),myLid,ier)
      if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

! make the master list of global proc ids
!
! allocate space to hold global ids
! only needed on root, but allocate everywhere to avoid complaints.
      allocate(Gprocids(mysize),stat=ier)
      if(ier/=0) call die(myname_,'allocate(Gprocids)',ier)
! gather over the LOCAL comm
      call MPI_GATHER(myGid,1,MP_INTEGER,Gprocids,1,MP_INTEGER,0,mycomms(i),ier)
      if(ier/=0) call die(myname_,'MPI_GATHER Gprocids',ier)

      if(myLid == 0) then
        call MPI_SEND(Gprocids,mysize,MP_INTEGER,0,myids(i),globalcomm,ier)
        if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(Gprocids)',ier)
      endif

      deallocate(Gprocids,stat=ier)
      if(ier/=0) call die(myname_,'deallocate(Gprocids)',ier)
    endif
  enddo

!  Global root waits for all sends
  if(myGid == 0) then
    call MPI_WAITALL(size(reqs), reqs, status, ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(Gprocids)',ier)
  endif

!  Now store the Gprocids in the World description and Broadcast

  if(myGid == 0) then
    allocate(root_idGprocid(ncomps,0:Gsize-1),stat=ier)
    if(ier/=0) call die(myname_,'allocate(root_idGprocid)',ier)

    root_idGprocid = transpose(tmparray)
  endif

  if(myGid /= 0) then
     allocate(root_nprocs(1),root_idGprocid(1,1),stat=ier)
     if(ier/=0) call die(myname_,'non-root allocate(root_idGprocid)',ier)
  endif

!!!!!!!!!!!!!!!!!!
! end of Gprocids
!!!!!!!!!!!!!!!!!!

! now call the init from root.
  call initr_(ncomps,globalcomm,root_nprocs,root_idGprocid)

! if(myGid==0 .or. myGid==17) then
!   write(*,*)'MCTA',myGid,ThisMCTWorld%ncomps,ThisMCTWorld%MCT_comm,ThisMCTWorld%nprocspid
!   do i=1,ThisMCTWorld%ncomps
!     write(*,*)'MCTK',myGid,i,ThisMCTWorld%idGprocid(i,0:ThisMCTWorld%nprocspid(i)-1)
!   enddo
! endif

! deallocate temporary arrays
 deallocate(root_nprocs,root_idGprocid,stat=ier)
 if(ier/=0) call die(myname_,'deallocate(root_nprocs,..)',ier)
 if(myGid == 0) then
   deallocate(compids,reqs,status,nprocs,tmparray,stat=ier)
   if(ier/=0) call die(myname_,'deallocate(compids,..)',ier)
 endif

 end subroutine initm_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: initd_ - initialize MCTWorld
!
! !DESCRIPTION:
! Do a distributed init of MCTWorld using the total number of components 
! {\tt ncomps} and either a unique integer component id {\tt myid} or,
! if more than one model is placed on a processor, an array of integer ids 
! specifying the models {\tt myids}.  Also required is
! the local communicator {\tt mycomm} and global communicator {\tt globalcomm}
! which encompasses all the models (typically this can be MPI\_COMM\_WORLD).
! This routine must be called once by each component (using {\em myid}) or
! component group (using {\em myids}).
!
! !INTERFACE:

! Subprogram not used  subroutine initd_(ncomps,globalcomm,mycomm,myid,myids)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_mpif90
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer, intent(in)	       :: ncomps          ! number of components
! Subprogram not used       integer, intent(in)	       :: globalcomm      ! global communicator
! Subprogram not used       integer, intent(in)	       :: mycomm          ! my communicator
! Subprogram not used       integer, intent(in),optional     :: myid            ! my component id
! Subprogram not used       integer, dimension(:),pointer,optional  :: myids    ! component ids
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! Subprogram not used ! 07Feb01 - R. Jacob <jacob@mcs.anl.gov> - non fatal error
! Subprogram not used !           if init is called a second time.
! Subprogram not used ! 08Feb01 - R. Jacob <jacob@mcs.anl.gov> - initialize the new
! Subprogram not used !           mygrank and mylrank
! Subprogram not used ! 20Apr01 - R. Jacob <jacob@mcs.anl.gov> - remove allids from
! Subprogram not used !           MCTWorld datatype.  Not needed because component
! Subprogram not used !           ids are always from 1 to number-of-components.
! Subprogram not used ! 22Jun01 - R. Jacob <jacob@mcs.anl.gov> - move Bcast and init
! Subprogram not used !           of MCTWorld to initr_
! Subprogram not used ! 20Sep07 - T. Craig migrated code to new initm routine
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initd_'
! Subprogram not used   integer :: msize,ier
! Subprogram not used   integer, dimension(:), pointer :: mycomm1d,myids1d
! Subprogram not used 
! Subprogram not used ! ------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! only one of myid and myids should be present
! Subprogram not used   if(present(myid) .and. present(myids)) then
! Subprogram not used     write(stderr,'(2a)') myname_, &
! Subprogram not used       'MCTERROR:  Must define myid or myids in MCTWord init'
! Subprogram not used       call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(.not.present(myid) .and. .not.present(myids)) then
! Subprogram not used     write(stderr,'(2a)') myname_, &
! Subprogram not used       'MCTERROR:  Must define one of myid or myids in MCTWord init'
! Subprogram not used       call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if (present(myids)) then
! Subprogram not used      msize = size(myids)
! Subprogram not used   else
! Subprogram not used      msize = 1
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   allocate(mycomm1d(msize),myids1d(msize),stat=ier)
! Subprogram not used   if(ier/=0) call die(myname_,'non-root allocate(root_idGprocid)',ier)
! Subprogram not used   mycomm1d(:) = mycomm
! Subprogram not used 
! Subprogram not used   if (present(myids)) then
! Subprogram not used      myids1d(:) = myids(:)
! Subprogram not used   else
! Subprogram not used      myids1d(:) = myid
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call initm_(ncomps,globalcomm,mycomm1d,myids1d)
! Subprogram not used 
! Subprogram not used   deallocate(mycomm1d,myids1d)
! Subprogram not used 
! Subprogram not used  end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initr_ - initialize MCTWorld from global root
!
! !DESCRIPTION:
! Initialize MCTWorld using information valid only on the global root.
! This is called by initm\_ but could also be called by the user
! for very complex model--processor geometries.
!
! !INTERFACE:

 subroutine initr_(ncomps,globalcomm,rnprocspid,ridGprocid)
!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in)                :: ncomps     ! total number of components
      integer, intent(in)                :: globalcomm ! the global communicator
      integer, dimension(:),intent(in)   :: rnprocspid ! number of processors for each component
      integer, dimension(:,:),intent(in) :: ridGprocid ! an array of size (1:ncomps) x (0:Gsize-1) 
						       ! which maps local ranks to global ranks
                                                       ! it's actually 1:Gsize here

! !REVISION HISTORY:
! 22Jun01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: ier,Gsize,myGid,MCTcomm,i,j

! Check that ncomps is a legal value
  if(ncomps < 1) then
     call die(myname_, "argument ncomps can't less than one!",ncomps)
  endif

! determine overall size
  call MP_comm_size(globalcomm,Gsize,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

! determine my rank in comm_world
  call MP_comm_rank(globalcomm,myGid,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

! create the MCT comm world
  call MP_comm_dup(globalcomm,MCTcomm,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_dup()',ier)

  allocate(ThisMCTWorld%nprocspid(ncomps),stat=ier)
  if(ier/=0) call die(myname_,'allocate(MCTWorld%nprocspid(:),...',ier)
  allocate(ThisMCTWorld%idGprocid(ncomps,0:Gsize-1),stat=ier)
  if(ier/=0) call die(myname_,'allocate(MCTWorld%nprocspid(:),...',ier)

!  set the MCTWorld
  ThisMCTWorld%ncomps = ncomps
  ThisMCTWorld%MCT_comm = MCTcomm
  ThisMCTWorld%mygrank = myGid

! Now store the component ids in the World description and Broadcast
  if(myGid == 0) then
    ThisMCTWorld%nprocspid(1:ncomps) = rnprocspid(1:ncomps)
    ThisMCTWorld%idGprocid = ridGprocid
  endif

  call MPI_BCAST(ThisMCTWorld%nprocspid, ncomps, MP_INTEGER, 0, MCTcomm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCast nprocspid',ier)

  call MPI_BCAST(ThisMCTWorld%idGprocid, ncomps*Gsize,MP_INTEGER, 0,MCTcomm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCast Gprocids',ier)

! if(myGid==17) then
!      do i=1,ThisMCTWorld%ncomps
!       do j=1,ThisMCTWorld%nprocspid(i)
!     write(*,*)'MCTK',myGid,i,j-1,ThisMCTWorld%idGprocid(i,j-1)
!    enddo
!   enddo
! endif

 end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a MCTWorld
!
! !DESCRIPTION:
! This routine deallocates the arrays of {\tt ThisMCTWorld}
! It also zeros out the integer components.
!
! !INTERFACE:

! Subprogram not used     subroutine clean_()
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! Subprogram not used ! 08Feb01 - R. Jacob <jacob@mcs.anl.gov> - clean the new
! Subprogram not used !           mygrank and mylrank
! Subprogram not used ! 20Apr01 - R. Jacob <jacob@mcs.anl.gov> - remove allids from
! Subprogram not used !           MCTWorld datatype.  Not needed because component
! Subprogram not used !           ids are always from 1 to number-of-components.
! Subprogram not used ! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - remove myid,mynprocs
! Subprogram not used !           and mylrank.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::clean_'
! Subprogram not used   integer :: ier
! Subprogram not used 
! Subprogram not used   deallocate(ThisMCTWorld%nprocspid,ThisMCTWorld%idGprocid,stat=ier)
! Subprogram not used   if(ier /= 0) call warn(myname_,'deallocate(MCTW,...)',ier)
! Subprogram not used 
! Subprogram not used   ThisMCTWorld%ncomps = 0
! Subprogram not used   ThisMCTWorld%mygrank = 0
! Subprogram not used 
! Subprogram not used  end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: NumComponents_ - Determine number of components in World.
!
! !DESCRIPTION:
! The function {\tt NumComponents\_} takes an input {\tt MCTWorld} 
! argument {\tt World}, and returns the number of component models 
! present.
!
! !INTERFACE:

! Subprogram not used  integer function NumComponents_(World)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(MCTWorld), intent(in)      :: World
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::NumComponents_'
! Subprogram not used 
! Subprogram not used   integer :: ncomps
! Subprogram not used 
! Subprogram not used   ncomps = World%ncomps
! Subprogram not used 
! Subprogram not used   if(ncomps <= 0) then
! Subprogram not used      write(stderr,'(2a,1i3)') myname,":: invalid no. of components = ",ncomps
! Subprogram not used      call die(myname_,'ncomps = ',ncomps)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   NumComponents_ = ncomps
! Subprogram not used 
! Subprogram not used  end function NumComponents_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComponentNumProcs_ - Number of processes a component owns.
!
! !DESCRIPTION:
! The function {\tt ComponentNumProcs\_} takes an input {\tt MCTWorld} 
! argument {\tt World}, and a component ID {\tt comp\_id}, and returns 
! the number of processes owned by that component.
!
! !INTERFACE:

! Subprogram not used  integer function ComponentNumProcs_(World, comp_id)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used       type(MCTWorld), intent(in)      :: World
! Subprogram not used       integer,        intent(in)      :: comp_id
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
! Subprogram not used ! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - modify to use
! Subprogram not used !           nprocspid and comp_id instead of World%mynprocs
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ComponentNumPros_'
! Subprogram not used 
! Subprogram not used   integer :: mynprocs
! Subprogram not used 
! Subprogram not used   mynprocs = World%nprocspid(comp_id)
! Subprogram not used 
! Subprogram not used   if(mynprocs <= 0) then
! Subprogram not used      write(stderr,'(2a,1i6)') myname,":: invalid no. of processes = ",mynprocs
! Subprogram not used      call die(myname_,'Number of processes = ',mynprocs)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   ComponentNumProcs_ = mynprocs
! Subprogram not used 
! Subprogram not used  end function ComponentNumProcs_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComponentToWorldRank_ - Determine rank on COMM_WORLD.
!
! !DESCRIPTION:
! The function {\tt ComponentToWorldRank\_} takes an input component ID 
! {\tt comp\_id} and input rank on that component communicator 
! {\tt comp\_rank}, and returns the rank of that process on the world 
! communicator of {\tt MCTWorld}.  
!
! !INTERFACE:

! Subprogram not used  integer function ComponentToWorldRank_(comp_rank, comp_id, World)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used   
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used       integer, intent(in)	     :: comp_rank ! process rank on the communicator
! Subprogram not used                                                   ! associated with comp_id
! Subprogram not used       integer, intent(in)	     :: comp_id   ! component id
! Subprogram not used       type(MCTWorld), intent(in)     :: World     ! World
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !       05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
! Subprogram not used !       14Jul02 - E. Ong <eong@mcs.anl.gov> - made argument checking required
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ComponentToWorldRank_'
! Subprogram not used 
! Subprogram not used   logical :: valid
! Subprogram not used   integer :: n, world_rank
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       ! Do we want the potentially time-consuming argument checks?
! Subprogram not used       ! The first time we use this function during execution on a
! Subprogram not used       ! given set of components and component ranks, we will.  In
! Subprogram not used       ! later invocations, these argument checks are probably not
! Subprogram not used       ! necessary (unless one alters MCTWorld), and impose a cost
! Subprogram not used       ! one may wish to avoid.
! Subprogram not used 
! Subprogram not used       ! These checks are just conditional statements and are 
! Subprogram not used       ! not particularly time-consuming. It's better to be safe
! Subprogram not used       ! than sorry. -EONG
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       ! Check argument comp_id for validity--assume initially it is not...
! Subprogram not used 
! Subprogram not used   valid = .false.
! Subprogram not used   n = 0
! Subprogram not used 
! Subprogram not used   if((comp_id <= World%ncomps) .and. &
! Subprogram not used        (comp_id > 0)) then
! Subprogram not used      valid = .true.
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used   if(.not. valid) then
! Subprogram not used      write(stderr,'(2a,1i7)') myname,":: invalid component id no. = ",&
! Subprogram not used 	  comp_id
! Subprogram not used      call die(myname_,'invalid comp_id = ',comp_id)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used       ! Check argument comp_rank for validity on the communicator associated
! Subprogram not used       ! with comp_id.  Assume initialy it is invalid.
! Subprogram not used 
! Subprogram not used   valid = .false.
! Subprogram not used   
! Subprogram not used   if((0 <= comp_rank) .or. &
! Subprogram not used        (comp_rank < ComponentNumProcs_(World, comp_id))) then
! Subprogram not used      valid = .true.
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(.not. valid) then
! Subprogram not used      write(stderr,'(2a,1i5,1a,1i2)') myname, &
! Subprogram not used 	  ":: invalid process ID. = ", &
! Subprogram not used 	  comp_rank, "on component ",comp_id
! Subprogram not used      call die(myname_,'invalid comp_rank = ',comp_rank)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       ! If we have reached this point, the input data are valid.
! Subprogram not used       ! Return the global rank for comp_rank on component comp_id
! Subprogram not used 
! Subprogram not used   world_rank = World%idGprocid(comp_id, comp_rank)
! Subprogram not used 
! Subprogram not used   if(world_rank < 0) then
! Subprogram not used      write(stderr,'(2a,1i6)') myname,":: negative world rank = ",world_rank
! Subprogram not used      call die(myname_,'negative world rank = ',world_rank)
! Subprogram not used   endif    
! Subprogram not used 
! Subprogram not used   ComponentToWorldRank_ = world_rank
! Subprogram not used 
! Subprogram not used  end function ComponentToWorldRank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComponentRootRank_ - Rank of component root on COMM_WORLD.
!
! !DESCRIPTION:
! The function {\tt ComponentRootRank\_} takes an input component ID 
! {\tt comp\_id} and input {\tt MCTWorld} variable {\tt World}, and
! returns the global rank of the root of this component.  
!
! !INTERFACE:

! Subprogram not used  integer function ComponentRootRank_(comp_id, World)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used       integer, intent(in)	     :: comp_id   ! component id
! Subprogram not used       type(MCTWorld), intent(in)     :: World     ! World
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !       05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
! Subprogram not used !       14Jul02 - E. Ong <eong@mcs.anl.gov> - made argument checking required
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ComponentRootRank_'
! Subprogram not used 
! Subprogram not used   integer :: world_comp_root
! Subprogram not used 
! Subprogram not used       ! Call ComponentToWorldRank_ assuming the root on a remote component
! Subprogram not used       ! has rank zero on the communicator associated with that component.
! Subprogram not used 
! Subprogram not used   world_comp_root = ComponentToWorldRank_(0, comp_id, World)
! Subprogram not used 
! Subprogram not used   if(world_comp_root < 0) then
! Subprogram not used      write(stderr,'(2a,1i6)') myname,":: negative world rank = ",& 
! Subprogram not used 	  world_comp_root
! Subprogram not used      call die(myname_,'invalid root id = ',world_comp_root)
! Subprogram not used   endif    
! Subprogram not used 
! Subprogram not used   ComponentRootRank_ = world_comp_root
! Subprogram not used 
! Subprogram not used  end function ComponentRootRank_

 end module m_MCTWorld

