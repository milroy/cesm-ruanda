!BOP -------------------------------------------------------------------
!
! !MODULE: m_FcComms - MPI collective communication operators
!                      with explict flow control
!
! !DESCRIPTION:
!
! This module includes implementations of MPI collective operators that
! have proven problematic on certain systems when run at scale. By 
! introducing additonal flow control, these problems (exhausting internal
! system resources) can be avoided. These routines were ported from
! the Community Atmosphere Model's spmd_utils.F90.
!
! !INTERFACE:
!
! Workaround for performance issue with rsend on cray systems with
! gemini interconnect
!

 module m_FcComms

      implicit none

      private	! except

      public :: fc_gather_int  ! flow control version of mpi_gather for integer vectors
      public :: fc_gather_fp   ! flow control version of mpi_gather for FP vectors
      public :: fc_gatherv_int ! flow control version of mpi_gatherv for integer vectors
      public :: fc_gatherv_fp  ! flow control version of mpi_gatherv for integer vectors
      public :: get_fcblocksize ! get current value of max_gather_block_size
      public :: set_fcblocksize ! set current value of max_gather_block_size


! !REVISION HISTORY:
! 30Jan09 - P.H. Worley <worleyph@ornl.gov> - imported routines
!           from CAM's spmd_utils to create this module.

  integer, public :: max_gather_block_size = 64
  character(len=*),parameter :: myname='MCT(MPEU)::m_FcComms'

 contains

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gather_int - Gather an array of type integer
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em integer} 
! to the {\tt root} process. Explicit handshaking messages are used
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gather
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          min(max(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, max_gather_block_size is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gather_int (sendbuf, sendcnt, sendtype, &
                             recvbuf, recvcnt, recvtype, &
                             root, comm, flow_cntl )
!
! !USES:
!
      use m_die
      use m_mpif90
!
! !INPUT PARAMETERS: 
!
      integer, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcnt
      integer, intent(in) :: sendtype
      integer, intent(in) :: recvcnt
      integer, intent(in) :: recvtype
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional, intent(in) :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      integer, intent(out) :: recvbuf(*)

! !REVISION HISTORY:
! 30Jan09 - P.H. Worley - imported from spmd_utils.F90
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::fc_gather_int'

 integer :: signal
 logical fc_gather         ! use explicit flow control?
 integer gather_block_size ! number of preposted receive requests

 integer :: mytid, mysize, mtag, p, i, count, displs
 integer :: preposts, head, tail
 integer :: rcvid(max_gather_block_size)
 integer :: status(MP_STATUS_SIZE)
 integer :: ier ! MPI error code

 signal = 1
 if ( present(flow_cntl) ) then
    if (flow_cntl >= 0) then
       gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
       fc_gather = .true.
    else
       fc_gather = .false.
    endif
 else
    gather_block_size = max(1,max_gather_block_size)
    fc_gather = .true.
 endif

 if (fc_gather) then
 
    call mpi_comm_rank (comm, mytid, ier)
    call mpi_comm_size (comm, mysize, ier)
    mtag = 0
    if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
       preposts = min(mysize-1, gather_block_size)
       head = 0
       count = 0
       do p=0, mysize-1
          if (p .ne. root) then
             if (recvcnt > 0) then
                count = count + 1
                if (count > preposts) then
                   tail = mod(head,preposts) + 1
                   call mpi_wait (rcvid(tail), status, ier)
                end if
                head = mod(head,preposts) + 1
                displs = p*recvcnt
                call mpi_irecv ( recvbuf(displs+1), recvcnt, &
                                 recvtype, p, mtag, comm, rcvid(head), &
                                 ier )
                call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
             end if
          end if
       end do

! copy local data
       displs = mytid*recvcnt
       do i=1,sendcnt
          recvbuf(displs+i) = sendbuf(i)
       enddo

! wait for final data
       do i=1,min(count,preposts)
          call mpi_wait (rcvid(i), status, ier)
       enddo

    else

       if (sendcnt > 0) then
          call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                          status, ier )
          call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                           comm, ier )
       end if

    endif
    if (ier /= 0) then
       call MP_perr_die(myname_,':: (point-to-point implementation)',ier)
    end if

 else
 
    call mpi_gather (sendbuf, sendcnt, sendtype, &
                     recvbuf, recvcnt, recvtype, &
                     root, comm, ier)
    if (ier /= 0) then
       call MP_perr_die(myname_,':: MPI_GATHER',ier)
    end if

 endif

 return
 end subroutine fc_gather_int

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gather_fp - Gather an array of type FP
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em FP} to
! the {\tt root} process. Explicit handshaking messages are used
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gather
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          min(max(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, max_gather_block_size is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
! Subprogram not used    subroutine fc_gather_fp (sendbuf, sendcnt, sendtype, &
! Subprogram not used                             recvbuf, recvcnt, recvtype, &
! Subprogram not used                              root, comm, flow_cntl )
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_realkinds, only : FP
! Subprogram not used       use m_die
! Subprogram not used       use m_mpif90
! Subprogram not used !
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       real (FP), intent(in)  :: sendbuf(*)
! Subprogram not used       integer, intent(in) :: sendcnt
! Subprogram not used       integer, intent(in) :: sendtype
! Subprogram not used       integer, intent(in) :: recvcnt
! Subprogram not used       integer, intent(in) :: recvtype
! Subprogram not used       integer, intent(in) :: root
! Subprogram not used       integer, intent(in) :: comm
! Subprogram not used       integer, optional, intent(in) :: flow_cntl
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       real (FP), intent(out) :: recvbuf(*)
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 30Jan09 - P.H. Worley - imported from spmd_utils.F90
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used  character(len=*),parameter :: myname_=myname//'::fc_gather_fp'
! Subprogram not used 
! Subprogram not used  real (FP) :: signal
! Subprogram not used  logical fc_gather         ! use explicit flow control?
! Subprogram not used  integer gather_block_size ! number of preposted receive requests
! Subprogram not used 
! Subprogram not used  integer :: mytid, mysize, mtag, p, i, count, displs
! Subprogram not used  integer :: preposts, head, tail
! Subprogram not used  integer :: rcvid(max_gather_block_size)
! Subprogram not used  integer :: status(MP_STATUS_SIZE)
! Subprogram not used  integer :: ier ! MPI error code
! Subprogram not used 
! Subprogram not used  signal = 1.0
! Subprogram not used  if ( present(flow_cntl) ) then
! Subprogram not used     if (flow_cntl >= 0) then
! Subprogram not used        gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
! Subprogram not used        fc_gather = .true.
! Subprogram not used     else
! Subprogram not used        fc_gather = .false.
! Subprogram not used     endif
! Subprogram not used  else
! Subprogram not used     gather_block_size = max(1,max_gather_block_size)
! Subprogram not used     fc_gather = .true.
! Subprogram not used  endif
! Subprogram not used 
! Subprogram not used  if (fc_gather) then
! Subprogram not used  
! Subprogram not used     call mpi_comm_rank (comm, mytid, ier)
! Subprogram not used     call mpi_comm_size (comm, mysize, ier)
! Subprogram not used     mtag = 0
! Subprogram not used     if (root .eq. mytid) then
! Subprogram not used 
! Subprogram not used ! prepost gather_block_size irecvs, and start receiving data
! Subprogram not used        preposts = min(mysize-1, gather_block_size)
! Subprogram not used        head = 0
! Subprogram not used        count = 0
! Subprogram not used        do p=0, mysize-1
! Subprogram not used           if (p .ne. root) then
! Subprogram not used              if (recvcnt > 0) then
! Subprogram not used                 count = count + 1
! Subprogram not used                 if (count > preposts) then
! Subprogram not used                    tail = mod(head,preposts) + 1
! Subprogram not used                    call mpi_wait (rcvid(tail), status, ier)
! Subprogram not used                 end if
! Subprogram not used                 head = mod(head,preposts) + 1
! Subprogram not used                 displs = p*recvcnt
! Subprogram not used                 call mpi_irecv ( recvbuf(displs+1), recvcnt, &
! Subprogram not used                                  recvtype, p, mtag, comm, rcvid(head), &
! Subprogram not used                                  ier )
! Subprogram not used                 call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used ! copy local data
! Subprogram not used        displs = mytid*recvcnt
! Subprogram not used        do i=1,sendcnt
! Subprogram not used           recvbuf(displs+i) = sendbuf(i)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used ! wait for final data
! Subprogram not used        do i=1,min(count,preposts)
! Subprogram not used           call mpi_wait (rcvid(i), status, ier)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        if (sendcnt > 0) then
! Subprogram not used           call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
! Subprogram not used                           status, ier )
! Subprogram not used           call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
! Subprogram not used                            comm, ier )
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        call MP_perr_die(myname_,':: (point-to-point implementation)',ier)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used  else
! Subprogram not used  
! Subprogram not used     call mpi_gather (sendbuf, sendcnt, sendtype, &
! Subprogram not used                      recvbuf, recvcnt, recvtype, &
! Subprogram not used                       root, comm, ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        call MP_perr_die(myname_,':: MPI_GATHER',ier)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used  endif
! Subprogram not used 
! Subprogram not used  return
! Subprogram not used  end subroutine fc_gather_fp

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_int - Gather an array of type integer
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em integer} 
! to the {\tt root} process. Explicit handshaking messages are used
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          min(max(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, max_gather_block_size is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_int (sendbuf, sendcnt, sendtype, &
                              recvbuf, recvcnts, displs, recvtype, &
                              root, comm, flow_cntl )
!
! !USES:
!
      use m_die
      use m_mpif90
!
! !INPUT PARAMETERS: 
!
      integer, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcnt
      integer, intent(in) :: sendtype
      integer, intent(in) :: recvcnts(*)
      integer, intent(in) :: displs(*)
      integer, intent(in) :: recvtype
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional, intent(in) :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      integer, intent(out) :: recvbuf(*)

! !REVISION HISTORY:
! 30Jan09 - P.H. Worley - imported from spmd_utils.F90
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::fc_gatherv_int'

 integer :: signal
 logical fc_gather         ! use explicit flow control?
 integer gather_block_size ! number of preposted receive requests

 integer :: mytid, mysize, mtag, p, q, i, count
 integer :: preposts, head, tail
 integer :: rcvid(max_gather_block_size)
 integer :: status(MP_STATUS_SIZE)
 integer :: ier ! MPI error code

 signal = 1
 if ( present(flow_cntl) ) then
    if (flow_cntl >= 0) then
       gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
       fc_gather = .true.
    else
       fc_gather = .false.
    endif
 else
    gather_block_size = max(1,max_gather_block_size)
    fc_gather = .true.
 endif

 if (fc_gather) then
 
    call mpi_comm_rank (comm, mytid, ier)
    call mpi_comm_size (comm, mysize, ier)
    mtag = 0
    if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
       preposts = min(mysize-1, gather_block_size)
       head = 0
       count = 0
       do p=0, mysize-1
          if (p .ne. root) then
             q = p+1
             if (recvcnts(q) > 0) then
                count = count + 1
                if (count > preposts) then
                   tail = mod(head,preposts) + 1
                   call mpi_wait (rcvid(tail), status, ier)
                end if
                head = mod(head,preposts) + 1
                call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                 recvtype, p, mtag, comm, rcvid(head), &
                                 ier )
                call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
             end if
          end if
       end do

! copy local data
       q = mytid+1
       do i=1,sendcnt
          recvbuf(displs(q)+i) = sendbuf(i)
       enddo

! wait for final data
       do i=1,min(count,preposts)
          call mpi_wait (rcvid(i), status, ier)
       enddo

    else

       if (sendcnt > 0) then
          call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                          status, ier )
          call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                           comm, ier )
       end if

    endif
    if (ier /= 0) then
       call MP_perr_die(myname_,':: (point-to-point implementation)',ier)
    end if

 else
 
    call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                      recvbuf, recvcnts, displs, recvtype, &
                      root, comm, ier)
    if (ier /= 0) then
       call MP_perr_die(myname_,':: MPI_GATHERV',ier)
    end if

 endif

 return
 end subroutine fc_gatherv_int

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_fp - Gather an array of type FP
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em FP} to
! the {\tt root} process. Explicit handshaking messages are used
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          min(max(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, max_gather_block_size is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_fp (sendbuf, sendcnt, sendtype, &
                             recvbuf, recvcnts, displs, recvtype, &
                             root, comm, flow_cntl )
!
! !USES:
!
      use m_realkinds, only : FP
      use m_die
      use m_mpif90
!
! !INPUT PARAMETERS: 
!
      real (FP), intent(in)  :: sendbuf(*)
      integer, intent(in) :: sendcnt
      integer, intent(in) :: sendtype
      integer, intent(in) :: recvcnts(*)
      integer, intent(in) :: displs(*)
      integer, intent(in) :: recvtype
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional, intent(in) :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      real (FP), intent(out) :: recvbuf(*)

! !REVISION HISTORY:
! 30Jan09 - P.H. Worley - imported from spmd_utils.F90
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::fc_gatherv_fp'

 real (FP) :: signal
 logical fc_gather         ! use explicit flow control?
 integer gather_block_size ! number of preposted receive requests

 integer :: mytid, mysize, mtag, p, q, i, count
 integer :: preposts, head, tail
 integer :: rcvid(max_gather_block_size)
 integer :: status(MP_STATUS_SIZE)
 integer :: ier ! MPI error code

 signal = 1.0
 if ( present(flow_cntl) ) then
    if (flow_cntl >= 0) then
       gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
       fc_gather = .true.
    else
       fc_gather = .false.
    endif
 else
    gather_block_size = max(1,max_gather_block_size)
    fc_gather = .true.
 endif

 if (fc_gather) then
 
    call mpi_comm_rank (comm, mytid, ier)
    call mpi_comm_size (comm, mysize, ier)
    mtag = 0
    if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
       preposts = min(mysize-1, gather_block_size)
       head = 0
       count = 0
       do p=0, mysize-1
          if (p .ne. root) then
             q = p+1
             if (recvcnts(q) > 0) then
                count = count + 1
                if (count > preposts) then
                   tail = mod(head,preposts) + 1
                   call mpi_wait (rcvid(tail), status, ier)
                end if
                head = mod(head,preposts) + 1
                call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                 recvtype, p, mtag, comm, rcvid(head), &
                                 ier )
                call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
             end if
          end if
       end do

! copy local data
       q = mytid+1
       do i=1,sendcnt
          recvbuf(displs(q)+i) = sendbuf(i)
       enddo

! wait for final data
       do i=1,min(count,preposts)
          call mpi_wait (rcvid(i), status, ier)
       enddo

    else

       if (sendcnt > 0) then
          call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                          status, ier )
          call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                           comm, ier )
       end if

    endif
    if (ier /= 0) then
       call MP_perr_die(myname_,':: (point-to-point implementation)',ier)
    end if

 else
 
    call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                      recvbuf, recvcnts, displs, recvtype, &
                      root, comm, ier)
    if (ier /= 0) then
       call MP_perr_die(myname_,':: MPI_GATHERV',ier)
    end if

 endif

 return
 end subroutine fc_gatherv_fp

!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_fcblocksize - return max_gather_block_size
!
! !DESCRIPTION:
! This function returns the current value of max_gather_block_size
!
! !INTERFACE:

! Subprogram not used  function get_fcblocksize()
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! No external modules are used by this function.
! Subprogram not used 
! Subprogram not used      implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used     integer           :: get_fcblocksize
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !       03Mar09 - R. Jacob (jacob@mcs.anl.gov) -- intial version
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::get_fcblocksize'
! Subprogram not used 
! Subprogram not used   get_fcblocksize = max_gather_block_size
! Subprogram not used 
! Subprogram not used  end function  get_fcblocksize

!BOP -------------------------------------------------------------------
!
! !IROUTINE: set_fcblocksize - set max_gather_block_size
!
! !DESCRIPTION:
! This function sets the current value of max_gather_block_size
!
! !INTERFACE:

! Subprogram not used  subroutine set_fcblocksize(gather_block_size)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! No external modules are used by this function.
! Subprogram not used 
! Subprogram not used      implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used     integer           :: gather_block_size
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !       03Mar09 - R. Jacob (jacob@mcs.anl.gov) -- intial version
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//':: set_fcblocksize'
! Subprogram not used 
! Subprogram not used   max_gather_block_size = gather_block_size
! Subprogram not used 
! Subprogram not used  end subroutine  set_fcblocksize

 end module m_FcComms
