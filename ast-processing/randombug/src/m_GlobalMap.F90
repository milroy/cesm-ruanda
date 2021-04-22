!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalMap - One-Dimensional Domain Decomposition Descriptor
!
! !DESCRIPTION:
! The {\tt GlobalMap} is a datatype used to store descriptors of a 
! one-dimensional domain decomposition for a vector on an MPI communicator.  
! It is defined with three assumptions:
! \begin{enumerate}
! \item Each process ID owns only one segment;
! \item No two segments in the decomposition overlap; and
! \item The segments are laid out in identical order to the MPI rank of 
! each process participating in the decomposition.
! \end{enumerate}
! per process ID).  It is the simpler of the two domain decomposition 
! descriptors offerd by MCT (the other being the {\tt GlobalSegMap}).  
! It consists of the following components:
! \begin{itemize}
! \item The MCT component identification number (see the module 
! {\tt m\_MCTWorld} for more information about MCT's component model 
! registry);
! \item The {\em global} number of elements in the distributed vector;
! \item The number of elements {\em stored locally};
! \item The number of elements {\em stored on each process} on the 
! communicator over which the vector is distributed; and
! \item The index of the elemnent {\em immediately before} the starting 
! element of each local segment (this choice allows for direct use of 
! this information with MPI's scatter and gather operations).  We refer 
! to this quantity as the {\em displacement} of the segment, a term used 
! both here and in the definition of the MCT {\tt Navigator} datatype.
! \end{itemize}
!
! Both the segment displacement and length data are stored in arrays 
! whose indices run from zero to $N-1$, where $N$ is the number of MPI 
! processes on the communicator on which the {\tt GlobalMap} is defined.
! This is done so this information corresponds directly to the MPI process 
! ID's on whihc the segments reside.
!
! This module contains the definition of the {\tt GlobalMap} datatype, 
! all-processor and an on-root creation methods (both of which can be 
! used to create a {\tt GlobalMap} on the local communicator), a creation 
! method to create/propagate a {\tt GlobalMap} native to a remote 
! communicator, a destruction method, and a variety of query methods.
! 
! !INTERFACE:

 module m_GlobalMap

! !USES
! No external modules are used in the declaration section of this module.

      implicit none

      private	! except

! !PUBLIC TYPES:

      public :: GlobalMap		! The class data structure

    Type GlobalMap
      integer :: comp_id                        ! Component ID number
      integer :: gsize				! the Global size
      integer :: lsize				! my local size
      integer,dimension(:),pointer :: counts	! all local sizes
      integer,dimension(:),pointer :: displs	! PE ordered locations
    End Type GlobalMap

! !PUBLIC MEMBER FUNCTIONS:

      public :: gsize
      public :: lsize
      public :: init
      public :: init_remote
      public :: clean
      public :: rank
      public :: bounds
      public :: comp_id

    interface gsize; module procedure gsize_; end interface
    interface lsize; module procedure lsize_; end interface
    interface init ; module procedure	&
       initd_,	&	! initialize from all PEs
       initr_		! initialize from the root
    end interface
    interface init_remote; module procedure init_remote_; end interface
    interface clean; module procedure clean_; end interface
    interface rank ; module procedure rank_ ; end interface
    interface bounds; module procedure bounds_; end interface
    interface comp_id ; module procedure comp_id_ ; end interface

! !SEE ALSO:
! The MCT module m_MCTWorld for more information regarding component 
! ID numbers.
!
! !REVISION HISTORY:
! 21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!  9Nov00 - J.W. Larson <larson@mcs.anl.gov> - added init_remote
!           interface.
! 26Jan01 - J.W. Larson <larson@mcs.anl.gov> - added storage for
!           component ID number GlobalMap%comp_id, and associated
!           method comp_id_()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_GlobalMap'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initd_ - Collective Creation on the Local Communicator
!
! !DESCRIPTION:
! This routine creates the {\tt GlobalMap} {\tt GMap} from distributed 
! data spread across the MPI communicatior associated with the input 
! {\tt INTEGER} handle {\tt comm}.  The {\tt INTEGER} input argument 
! {\tt comp\_id} is used to define the MCT component ID for {\tt GMap}.
! The input {\tt INTEGER} argument {\tt ln} is the number of elements 
! in the local vector segment.
!
! !INTERFACE:

! Subprogram not used  subroutine initd_(GMap, comp_id, ln, comm)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used       use m_mpif90
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,         intent(in)  :: comp_id ! Component ID
! Subprogram not used       integer,         intent(in)  :: ln      ! the local size
! Subprogram not used       integer,         intent(in)  :: comm    ! f90 MPI communicator 
! Subprogram not used                                               ! handle 
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(GlobalMap), intent(out) :: GMap
! Subprogram not used 
! Subprogram not used ! !SEE ALSO:
! Subprogram not used ! The MCT module m_MCTWorld for more information regarding component 
! Subprogram not used ! ID numbers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initd_'
! Subprogram not used   integer :: nPEs,myID,ier,l,i
! Subprogram not used 
! Subprogram not used   call MP_comm_size(comm,nPEs,ier)
! Subprogram not used   if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)
! Subprogram not used 
! Subprogram not used   call MP_comm_rank(comm,myID,ier)
! Subprogram not used   if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)
! Subprogram not used 
! Subprogram not used   allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
! Subprogram not used   if(ier /= 0) call die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   call MPI_allgather(ln,1,MP_INTEGER,GMap%counts,1,MP_INTEGER,comm,ier)
! Subprogram not used   if(ier/=0) call MP_perr_die(myname_,'MPI_allgather()',ier)
! Subprogram not used 
! Subprogram not used   l=0
! Subprogram not used   do i=0,nPEs-1
! Subprogram not used     GMap%displs(i)=l
! Subprogram not used     l=l+GMap%counts(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   GMap%lsize=GMap%counts(myID)	! the local size
! Subprogram not used   GMap%gsize=l	! the global size
! Subprogram not used   GMap%comp_id = comp_id ! the component ID number
! Subprogram not used 
! Subprogram not used  end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initr_ Create a GlobalMap from the Root Process
!
! !DESCRIPTION:
! This routine creates the {\tt GlobalMap} {\tt GMap}, and propagates 
! it to all processes on the communicator associated with the MPI 
! {\tt INTEGER} handle {\tt comm}.  The input {\tt INTEGER} arguments 
! {\tt comp\_id} (the MCT component ID number) and {\tt lns(:)} need 
! only be valid on the process whose rank is equal to {\tt root} on 
! {\tt comm}.  The array {\tt lns(:)} should have length equal to the 
! number of processes on {\tt comm}, and contains the length of each
! local segment.
!
! !INTERFACE:

 subroutine initr_(GMap, comp_id, lns, root, comm)

! !USES:

      use m_mpif90
      use m_die
      use m_stdio

      implicit none

! !INPUT PARAMETERS:

      integer,               intent(in)  :: comp_id ! component ID number
      integer, dimension(:), intent(in)  :: lns     ! segment lengths
      integer,               intent(in)  :: root    ! root process ID
      integer,               intent(in)  :: comm    ! communicator ID

! !OUTPUT PARAMETERS:

      type(GlobalMap),       intent(out) :: GMap

! !SEE ALSO:
! The MCT module m_MCTWorld for more information regarding component 
! ID numbers.
!
! !REVISION HISTORY:
! 29May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)


  if(myID == root) then
    if(size(lns(:)) /= nPEs) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(lns) =',size(lns),		&
	', nPEs =',nPEs
      call die(myname_)
    endif

    GMap%counts(:)=lns(:)
  endif

  call MPI_bcast(GMap%counts, nPEs, MP_INTEGER, root, comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

  ! on each process, use GMap%counts(:) to compute GMap%displs(:)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%gsize=l	! the global size

  ! finally, set and broadcast the component ID number GMap%comp_id

  if(myID == root) GMap%comp_id = comp_id

  call MPI_bcast(GMap%comp_id,1,MP_INTEGER,root,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

 end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_remote_ Initialize Remote GlobalMap from the Root
!
! !DESCRIPTION:
! This routine creates and propagates across the local communicator a 
! {\tt GlobalMap} associated with a remote component.  The controlling 
! process in this operation has MPI process ID defined by the input 
! {\tt INTEGER} argument {\tt my\_root}, and its MPI communinicator 
! is defined by the input {\tt INTEGER} argument {\tt my\_comm}.  The 
! input {\tt INTEGER} argument {\tt remote\_npes} is the number of MPI 
! processes on the remote component's communicator (which need be valid 
! only on the process {\tt my\_root}).  The input the {\tt INTEGER} 
! array {\tt remote\_lns(:)}, and the {\tt INTEGER} argument 
! {\tt remote\_comp\_id} need only be valid on the process 
! whose rank on the communicator {\tt my\_comm} is {\tt my\_root}.  The 
! argument {\tt remote\_lns(:)} defines the vector segment length on each 
! process of the remote component's communicator, and the argument 
! {\tt remote\_comp\_id} defines the remote component's ID number in 
! the MCT component registry {\tt MCTWorld}.
!
! !INTERFACE:

! Subprogram not used  subroutine init_remote_(GMap, remote_lns, remote_npes, my_root, &
! Subprogram not used                          my_comm, remote_comp_id)
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used       use m_mpif90
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer, dimension(:), intent(in)  :: remote_lns
! Subprogram not used       integer,               intent(in)  :: remote_npes
! Subprogram not used       integer,               intent(in)  :: my_root
! Subprogram not used       integer,               intent(in)  :: my_comm
! Subprogram not used       integer,               intent(in)  :: remote_comp_id 
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(GlobalMap),       intent(out) :: GMap
! Subprogram not used 
! Subprogram not used ! !SEE ALSO:
! Subprogram not used ! The MCT module m_MCTWorld for more information regarding component 
! Subprogram not used ! ID numbers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  8Nov00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! Subprogram not used ! 26Jan01 - J.W. Larson <larson@mcs.anl.gov> - slight change--remote
! Subprogram not used !           communicator is replaced by remote component ID number
! Subprogram not used !           in argument remote_comp_id.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::init_remote_'
! Subprogram not used   integer :: nPEs,myID,ier,l,i
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         ! Which processor am I on communicator my_comm?  Store
! Subprogram not used         ! the answer in myID:
! Subprogram not used 
! Subprogram not used   call MP_comm_rank(my_comm, myID, ier)
! Subprogram not used   if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)
! Subprogram not used 
! Subprogram not used         ! allocate counts and displacements component arrays
! Subprogram not used         ! for the sake of compactness, store the value of remote_npes
! Subprogram not used         ! in the more tersely named variable nPEs.
! Subprogram not used 
! Subprogram not used   if(myID == my_root) nPEs = remote_npes
! Subprogram not used 
! Subprogram not used   call MPI_bcast(nPEs, 1, MP_INTEGER, my_root, my_comm, ier)
! Subprogram not used   if(ier/=0) call MP_perr_die(myname_,'MPI_bcast(nPEs...)',ier)
! Subprogram not used 
! Subprogram not used   allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
! Subprogram not used   if(ier /= 0) call die(myname_,'allocate()',ier)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         ! On the Root processor, check the size of remote_lns(:)
! Subprogram not used         ! to see it is equal to nPEs, the number of remote processes,
! Subprogram not used         ! then store it as GMap%counts and broadcast it.
! Subprogram not used 
! Subprogram not used   if(myID == my_root) then
! Subprogram not used     if(size(remote_lns(:)) /= nPEs) then
! Subprogram not used       write(stderr,'(2a,2(a,i4))') myname_,	 &
! Subprogram not used 	': _root_ argument error',		 &
! Subprogram not used 	', size(remote_lns) =',size(remote_lns), &
! Subprogram not used 	', nPEs =',nPEs
! Subprogram not used       call die(myname_)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     GMap%counts(:)=remote_lns(:)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call MPI_bcast(GMap%counts, nPEs, MP_INTEGER, my_root, my_comm, ier)
! Subprogram not used   if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)
! Subprogram not used 
! Subprogram not used         ! Now, on each processor of my_comm, compute from 
! Subprogram not used         ! GMap%counts(:) the entries of GMap%displs(:)
! Subprogram not used 
! Subprogram not used   l=0
! Subprogram not used   do i=0,nPEs-1
! Subprogram not used     GMap%displs(i)=l
! Subprogram not used     l=l+GMap%counts(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   GMap%lsize = -1                ! In this case, the local size is invalid!!!
! Subprogram not used   GMap%gsize = l      	         ! the global size
! Subprogram not used 
! Subprogram not used         ! Finally, set GMap's component ID (recall only the value on
! Subprogram not used         ! process my_root is valid).
! Subprogram not used 
! Subprogram not used   if(myID == my_root)  GMap%comp_id = remote_comp_id
! Subprogram not used   call MPI_bcast(GMap%comp_id, 1, MP_INTEGER, my_root, my_comm,ier)
! Subprogram not used   if(ier/=0) call MP_perr_die(myname_,'MPI_bcast(GMap%comp_id...)',ier)
! Subprogram not used 
! Subprogram not used  end subroutine init_remote_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a GlobalMap
!
! !DESCRIPTION:
! This routine deallocates all allocated memory associated with the 
! input/output {\tt GlobalMap} argument {\tt GMap}, and sets to zero 
! all of its statically defined components.  The success (failure) of 
! this operation is signified by the zero (non-zero) value of the 
! optional output {\tt INTEGER} argument {\tt stat}.
!
! !INTERFACE:

 subroutine clean_(GMap, stat)

! !USES:

      use m_die

      implicit none

! !INPUT/OUTPUT PARAMETERS:

      type(GlobalMap),           intent(inout) :: GMap

! !OUTPUT PARAMETERS:

      integer,         optional, intent(out)   :: stat

! !REVISION HISTORY:
! 21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 26Jan01 - J. Larson <larson@mcs.anl.gov> incorporated comp_id.
!  1Mar02 - E.T. Ong <eong@mcs.anl.gov> removed the die to prevent
!           crashes and added stat argument.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(GMap%counts,GMap%displs,stat=ier)

  if(present(stat)) then
     stat=ier
  else
     if(ier /= 0) call warn(myname_,'deallocate(GMap%...)',ier)
  endif
  
  if(ier == 0) then


  endif

  GMap%lsize = 0
  GMap%gsize = 0
  GMap%comp_id = 0

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - Return Local Segment Length
!
! !DESCRIPTION:
! This {\tt INTEGER} function returns the length of the local vector 
! segment as defined by the input {\tt GlobalMap} argument {\tt GMap}.

! !INTERFACE:

 integer function lsize_(GMap)

! !USES:

      implicit none

! !INPUT PARAMETERS:

      type(GlobalMap), intent(in) :: GMap

! !REVISION HISTORY:
! 21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=GMap%lsize

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - Return Global Vector Length
!
! !DESCRIPTION:
! This {\tt INTEGER} function returns the global length of a vector 
! that is decomposed according to the input {\tt GlobalMap} argument 
! {\tt GMap}.
!
! !INTERFACE:

 integer function gsize_(GMap)

! !USES:

      implicit none

! !INPUT PARAMETERS:

      type(GlobalMap), intent(in) :: GMap


! !REVISION HISTORY:
! 21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=GMap%gsize

 end function gsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank_ - Process ID Location of a Given Vector Element
!
! !DESCRIPTION:
! This routine uses the input {\tt GlobalMap} argument {\tt GMap} to 
! determine the process ID (on the communicator on which {\tt GMap} was
! defined) of the vector element with global index {\tt i\_g}.  This 
! process ID is returned in the output {\tt INTEGER} argument {\tt rank}.
!
! !INTERFACE:

! Subprogram not used  subroutine rank_(GMap, i_g, rank)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(GlobalMap), intent(in)  :: GMap
! Subprogram not used       integer,         intent(in)  :: i_g
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,         intent(out) :: rank
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  5May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::rank_'
! Subprogram not used   integer :: i,ilc,ile
! Subprogram not used 
! Subprogram not used   rank=-1	! if nowhere fits
! Subprogram not used   do i=0,size(GMap%displs)-1
! Subprogram not used     ilc=GMap%displs(i)
! Subprogram not used     ile=ilc+GMap%counts(i)
! Subprogram not used 
! Subprogram not used 		! If i_g in (ilc,ile].  Note that i_g := [1:..]
! Subprogram not used 
! Subprogram not used     if(ilc < i_g .and. i_g <= ile) then
! Subprogram not used       rank=i
! Subprogram not used       return
! Subprogram not used     endif
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used  end subroutine rank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bounds_ - First/Last Global Indicies for a Process' Segment
!
! !DESCRIPTION:
! This routine takes as input a process ID (defined by the input 
! {\tt INTEGER} argument {\tt pe\_no}), examines the input {\tt GlobalMap} 
! argument {\tt GMap}, and returns the global indices for the first and 
! last elements of the segment owned by this process in the output 
! {\tt INTEGER} arguments {\tt lbnd} and {\tt ubnd}, respectively.
!
! !INTERFACE:

! Subprogram not used  subroutine bounds_(GMap, pe_no, lbnd, ubnd)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(GlobalMap), intent(in)  :: GMap
! Subprogram not used       integer,         intent(in)  :: pe_no
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,         intent(out) :: lbnd
! Subprogram not used       integer,         intent(out) :: ubnd
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 30Jan01 - J. Larson <larson@mcs.anl.gov> - initial code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::bounds_'
! Subprogram not used 
! Subprogram not used   lbnd = GMap%displs(pe_no) + 1
! Subprogram not used   ubnd = lbnd + GMap%counts(pe_no) - 1
! Subprogram not used 
! Subprogram not used  end subroutine bounds_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: comp_id_ - Return the Component ID Number
!
! !DESCRIPTION:
! This {\tt INTEGER} query function returns the MCT component ID number 
! stored in the input {\tt GlobalMap} argument {\tt GMap}.
!
! !INTERFACE:

! Subprogram not used  integer function comp_id_(GMap)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(GlobalMap), intent(in) :: GMap
! Subprogram not used 
! Subprogram not used ! !SEE ALSO:
! Subprogram not used ! The MCT module m_MCTWorld for more information regarding component 
! Subprogram not used ! ID numbers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 25Jan02 - J. Larson <larson@mcs.anl.gov> - initial version
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::comp_id_'
! Subprogram not used 
! Subprogram not used   comp_id_ = GMap%comp_id
! Subprogram not used 
! Subprogram not used  end function comp_id_

 end module m_GlobalMap
