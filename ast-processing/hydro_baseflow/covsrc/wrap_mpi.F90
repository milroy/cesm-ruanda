!----------------------------------------------------------------------- 
! 
! Purpose:
!
! 	Wrapper routines for the MPI (Message Passing) library for the
!	distributed memory (1) version of the code. Also data with
!	"shorthand" names for the MPI data types.
!
! Entry points:
!      mpibarrier             Calls mpi_barrier
!      mpifinalize            Calls mpi_finalize
!      mpipack_size           Calls mpi_pack
!      mpipack                Calls mpi_pack
!      mpiunpack              Calls mpi_unpack
!      mpisendrecv            Calls mpi_sendrecv
!      mpiisend               Calls mpi_isend
!      mpiirsend              Calls mpi_irsend
!      mpiissend              Calls mpi_issend
!      mpiirecv               Calls mpi_irecv
!      mpiwait                Calls mpi_wait
!      mpiwaitall             Calls mpi_waitall
!      mpisend                Calls mpi_send
!      mpirsend               Calls mpi_rsend
!      mpissend               Calls mpi_ssend
!      mpirecv                Calls mpi_recv
!      mpigather              Calls mpi_gather
!      mpigatherv             Calls mpi_gatherv
!      mpigathervr4           Calls mpi_gatherv for real*4 data
!      mpigathervint          Calls mpi_gatherv for integer data
!      mpisum                 Calls mpi_sum
!      mpiscatter             Calls mpi_scatter
!      mpiscatterv            Calls mpi_scatterv
!      mpibcast               Calls mpi_bcast
!      mpiallmaxint           Calls mpi_allreduce on integer vector with mpi_max operator
!      mpialltoallv           Calls mpi_alltoallv
!      mpialltoallint         Calls mpi_alltoall for integer data
!      mpiallgatherv          Calls mpi_allgatherv
!      mpiallgatherint        Calls mpi_allgatherv for integer data
!      mpiwincreate           Calls mpi_win_create and mpi_win_fence
!
! Author: Many
! 
!-----------------------------------------------------------------------
!

!
! Performance bug work around for Gemini interconnect
!





!
! Compile these routines only when 1 is defined
!


!****************************************************************

! Subprogram not used    subroutine mpibarrier (comm)
! Subprogram not used !
! Subprogram not used ! MPI barrier, have threads wait until all threads have reached this point
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used    call mpi_barrier (comm, ier)
! Subprogram not used    if (ier.ne.mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_barrier failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpibarrier
 
!****************************************************************
 
! Subprogram not used    subroutine mpifinalize
! Subprogram not used !
! Subprogram not used ! End of all MPI communication
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used    call mpi_finalize (ier)
! Subprogram not used    if (ier.ne.mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_finalize failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpifinalize
 
!****************************************************************
 
! Subprogram not used    subroutine mpipack_size (incount, datatype, comm, size)
! Subprogram not used !
! Subprogram not used ! Returns the size of the packed data
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in):: incount
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used    integer, intent(out):: size
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used    call mpi_pack_size (incount, datatype, comm, size, ier)
! Subprogram not used    if (ier.ne.mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_pack_size failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpipack_size
 
!****************************************************************
 
! Subprogram not used    subroutine mpipack (inbuf, incount, datatype, outbuf, outsize,    &
! Subprogram not used                        position, comm)
! Subprogram not used !
! Subprogram not used ! Pack the data and send it.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real(r8), intent(in):: inbuf(*)
! Subprogram not used    real(r8), intent(out):: outbuf(*)
! Subprogram not used    integer, intent(in):: incount
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(out):: outsize
! Subprogram not used    integer, intent(inout):: position
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used    call mpi_pack (inbuf, incount, datatype, outbuf, outsize,         &
! Subprogram not used                   position, comm, ier)
! Subprogram not used    if (ier.ne.mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_pack failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpipack
 
!****************************************************************
 
! Subprogram not used    subroutine mpiunpack (inbuf, insize, position, outbuf, outcount,  &
! Subprogram not used                          datatype, comm)
! Subprogram not used !
! Subprogram not used ! Un-packs the data from the packed receive buffer
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real(r8), intent(in):: inbuf(*)
! Subprogram not used    real(r8), intent(out):: outbuf(*)
! Subprogram not used    integer, intent(in):: insize
! Subprogram not used    integer, intent(inout):: position
! Subprogram not used    integer, intent(in):: outcount
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used    call mpi_unpack (inbuf, insize, position, outbuf, outcount,       &
! Subprogram not used                     datatype, comm, ier)
! Subprogram not used    if (ier.ne.mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_unpack failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiunpack
 
!****************************************************************
 
! Subprogram not used    subroutine mpisendrecv (sendbuf, sendcount, sendtype, dest, sendtag,  &
! Subprogram not used                            recvbuf, recvcount, recvtype, source,recvtag, &
! Subprogram not used                            comm)
! Subprogram not used !
! Subprogram not used ! Blocking send and receive.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real(r8), intent(in):: sendbuf(*)
! Subprogram not used    real(r8), intent(out):: recvbuf(*)
! Subprogram not used    integer, intent(in):: sendcount
! Subprogram not used    integer, intent(in):: sendtype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: sendtag
! Subprogram not used    integer, intent(in):: recvcount
! Subprogram not used    integer, intent(in):: recvtype
! Subprogram not used    integer, intent(in):: source
! Subprogram not used    integer, intent(in):: recvtag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer :: status(MPI_STATUS_SIZE)
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_sendrecv (sendbuf, sendcount, sendtype, dest, sendtag,   &
! Subprogram not used                       recvbuf, recvcount, recvtype, source, recvtag, &
! Subprogram not used                       comm, status, ier)
! Subprogram not used    if (ier.ne.mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_sendrecv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used !
! Subprogram not used ! ASSUME nrecv = nsend for stats gathering purposes.  This is not actually
! Subprogram not used ! correct, but its the best we can do since recvcount is a Max number
! Subprogram not used !
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nrecv = nrecv + 1
! Subprogram not used    nwsend = nwsend + sendcount
! Subprogram not used    nwrecv = nwrecv + sendcount
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpisendrecv
 
!****************************************************************
 
! Subprogram not used    subroutine mpiisend (buf, count, datatype, dest, tag, comm, request)
! Subprogram not used !
! Subprogram not used ! Does a non-blocking send.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used    integer, intent(out):: request
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_isend (buf, count, datatype, dest, tag, comm, request, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_isend failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nwsend = nwsend + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiisend
 
!****************************************************************
 
! Subprogram not used    subroutine mpiirsend (buf, count, datatype, dest, tag, comm, request)
! Subprogram not used !
! Subprogram not used ! Does a non-blocking ready send.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used    integer, intent(out):: request
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_irsend (buf, count, datatype, dest, tag, comm, request, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_irsend failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nwsend = nwsend + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiirsend
 
!****************************************************************
 
! Subprogram not used    subroutine mpiissend (buf, count, datatype, dest, tag, comm, request)
! Subprogram not used !
! Subprogram not used ! Does a non-blocking synchronous send.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used    integer, intent(out):: request
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_issend (buf, count, datatype, dest, tag, comm, request, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_issend failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nwsend = nwsend + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiissend
 
!****************************************************************
 
! Subprogram not used    subroutine mpiirecv (buf, count, datatype, source, tag, comm, request)
! Subprogram not used !
! Subprogram not used ! Does a non-blocking receive.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(out):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: source
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used    integer, intent(out):: request
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_irecv (buf, count, datatype, source, tag, comm, request, ier )
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_irecv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nrecv = nrecv + 1
! Subprogram not used    nwrecv = nwrecv + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiirecv
 
!****************************************************************
 
! Subprogram not used    subroutine mpiwait (request, status)
! Subprogram not used !
! Subprogram not used ! Waits for a nonblocking operation to complete.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(inout):: request
! Subprogram not used    integer, intent(out):: status
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_wait (request, status, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_wait failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiwait
 
!****************************************************************
 
! Subprogram not used    subroutine mpiwaitall (count, array_of_requests, array_of_statuses)
! Subprogram not used !
! Subprogram not used ! Waits for a collection of nonblocking operations to complete.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(inout):: array_of_requests(*)
! Subprogram not used    integer, intent(out):: array_of_statuses(*)
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_waitall (count, array_of_requests, array_of_statuses, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_waitall failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiwaitall
 
!****************************************************************
 
! Subprogram not used    subroutine mpisend (buf, count, datatype, dest, tag, comm)
! Subprogram not used !
! Subprogram not used ! Does a blocking send
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_send (buf, count, datatype, dest, tag, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_send failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nwsend = nwsend + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpisend
 
!****************************************************************
 
! Subprogram not used    subroutine mpirsend (buf, count, datatype, dest, tag, comm)
! Subprogram not used !
! Subprogram not used ! Does a blocking ready send
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_rsend (buf, count, datatype, dest, tag, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_rsend failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nwsend = nwsend + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpirsend
 
!****************************************************************
 
! Subprogram not used    subroutine mpissend (buf, count, datatype, dest, tag, comm)
! Subprogram not used !
! Subprogram not used ! Does a blocking synchronous send
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: dest
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_ssend (buf, count, datatype, dest, tag, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_ssend failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nsend = nsend + 1
! Subprogram not used    nwsend = nwsend + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpissend
 
!****************************************************************
 
! Subprogram not used    subroutine mpirecv (buf, count, datatype, source, tag, comm)
! Subprogram not used !
! Subprogram not used ! Does a blocking receive
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(out):: buf(*)
! Subprogram not used    integer, intent(in):: count
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: source
! Subprogram not used    integer, intent(in):: tag
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer status (MPI_STATUS_SIZE) ! Status of message
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_recv (buf, count, datatype, source, tag, comm, status, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_recv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    nrecv = nrecv + 1
! Subprogram not used    nwrecv = nwrecv + count
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpirecv
 
!****************************************************************
 
! Subprogram not used    subroutine mpigather (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
! Subprogram not used                          recvtype, root, comm)
! Subprogram not used !
! Subprogram not used ! Collects different messages from each thread on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: sendbuf(*)
! Subprogram not used    real (r8), intent(out):: recvbuf(*)
! Subprogram not used    integer, intent(in):: sendcnt
! Subprogram not used    integer, intent(in):: sendtype
! Subprogram not used    integer, intent(in):: recvcnt
! Subprogram not used    integer, intent(in):: recvtype
! Subprogram not used    integer, intent(in):: root
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_gather (sendbuf, sendcnt, sendtype,                      &
! Subprogram not used                     recvbuf, recvcnt, recvtype, root, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_gather failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpigather
 
!****************************************************************
 
! Subprogram not used    subroutine mpigatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, &
! Subprogram not used                           displs, recvtype, root, comm)
! Subprogram not used !
! Subprogram not used ! Collects different messages from each thread on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils, only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in)  :: sendbuf(*)
! Subprogram not used    real (r8), intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in) :: displs(*)
! Subprogram not used    integer, intent(in) :: sendcnt
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: recvcnts(*)
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: root
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used  
! Subprogram not used    integer ier   ! MPI error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                      root, comm, ier)
! Subprogram not used    if (ier /= mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_gatherv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpigatherv
 
!****************************************************************
 
! Subprogram not used    subroutine mpigathervr4 (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, &
! Subprogram not used                           displs, recvtype, root, comm)
! Subprogram not used !
! Subprogram not used ! Collects different messages from each thread on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r4 => shr_kind_r4, r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils, only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r4), intent(in)  :: sendbuf(*)
! Subprogram not used    real (r4), intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in) :: displs(*)
! Subprogram not used    integer, intent(in) :: sendcnt
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: recvcnts(*)
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: root
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used  
! Subprogram not used    integer ier   ! MPI error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                      root, comm, ier)
! Subprogram not used    if (ier /= mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_gatherv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpigathervr4
 
!****************************************************************
 
! Subprogram not used    subroutine mpigathervint (sendbuf, sendcnt, sendtype, recvbuf, &
! Subprogram not used                              recvcnts, displs, recvtype, root, comm)
! Subprogram not used !
! Subprogram not used ! Collects different messages from each thread on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils, only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in)  :: sendbuf(*)
! Subprogram not used    integer, intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in) :: displs(*)
! Subprogram not used    integer, intent(in) :: sendcnt
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: recvcnts(*)
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: root
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used  
! Subprogram not used    integer ier   ! MPI error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                      root, comm, ier)
! Subprogram not used    if (ier /= mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_gatherv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpigathervint
 
!****************************************************************
 
! Subprogram not used    subroutine mpisum (sendbuf, recvbuf, cnt, datatype, root, comm)
! Subprogram not used !
! Subprogram not used ! Sums sendbuf across all processors on communicator, returning 
! Subprogram not used ! result to root.
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in):: sendbuf(*)
! Subprogram not used    real (r8), intent(out):: recvbuf(*)
! Subprogram not used    integer, intent(in):: cnt
! Subprogram not used    integer, intent(in):: datatype
! Subprogram not used    integer, intent(in):: root
! Subprogram not used    integer, intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_reduce (sendbuf, recvbuf, cnt, datatype, mpi_sum, &
! Subprogram not used                     root, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_reduce failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpisum
 
!****************************************************************
 
! Subprogram not used    subroutine mpiscatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
! Subprogram not used                           recvtype, root, comm)
! Subprogram not used !
! Subprogram not used ! Sends different messages from masterproc to each thread
! Subprogram not used ! 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8),intent(in):: sendbuf(*)
! Subprogram not used    real (r8), intent(out):: recvbuf(*)
! Subprogram not used    integer,intent(in):: sendcnt
! Subprogram not used    integer,intent(in):: sendtype
! Subprogram not used    integer,intent(in):: recvcnt
! Subprogram not used    integer,intent(in):: recvtype
! Subprogram not used    integer,intent(in):: root
! Subprogram not used    integer,intent(in):: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_scatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
! Subprogram not used                      recvtype, root, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_scatter failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiscatter
 
!****************************************************************
 
! Subprogram not used    subroutine mpiscatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, &
! Subprogram not used                            recvcnt, recvtype, root, comm)
! Subprogram not used !
! Subprogram not used ! Sends different messages from masterproc to each thread
! Subprogram not used ! 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in)  :: sendbuf(*)
! Subprogram not used    real (r8), intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in) :: displs(*)
! Subprogram not used    integer, intent(in) :: sendcnts(*)
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: recvcnt
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: root
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_scatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, &
! Subprogram not used                       recvtype, root, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_scatter failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiscatterv
 
!****************************************************************
 
   subroutine mpibcast (buffer, count, datatype, root, comm )
!
! Broadcasts a message from masterproc to all threads
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog




   implicit none

   real (r8), intent(inout):: buffer(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 



   call mpi_bcast (buffer, count, datatype, root, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_bcast failed ier=',ier
      call endrun
   end if



 
   return
   end subroutine mpibcast
!****************************************************************
 
! Subprogram not used    subroutine mpiallmaxint (sendbuf, recvbuf, count, comm)
! Subprogram not used !
! Subprogram not used ! Allreduce integer vector maximum
! Subprogram not used ! 
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in)  :: sendbuf(*)
! Subprogram not used    integer, intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in)  :: count
! Subprogram not used    integer, intent(in)  :: comm
! Subprogram not used  
! Subprogram not used    integer :: ier              ! MPI error code
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_allreduce (sendbuf, recvbuf, count, mpiint, &
! Subprogram not used                        mpimax, comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_allreduce failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpiallmaxint

!****************************************************************
 
! Subprogram not used    subroutine mpialltoallv (sendbuf, sendcnts, sdispls, sendtype, &
! Subprogram not used                             recvbuf, recvcnts, rdispls, recvtype, &
! Subprogram not used                             comm)
! Subprogram not used !
! Subprogram not used ! All-to-all scatter/gather
! Subprogram not used ! 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in)  :: sendbuf(*)
! Subprogram not used    real (r8), intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in) :: sdispls(*)
! Subprogram not used    integer, intent(in) :: sendcnts(*)
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: recvcnts(*)
! Subprogram not used    integer, intent(in) :: rdispls(*)
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used  
! Subprogram not used    integer :: ier              ! MPI error code
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_alltoallv (sendbuf, sendcnts, sdispls, sendtype, &
! Subprogram not used                        recvbuf, recvcnts, rdispls, recvtype, &
! Subprogram not used                        comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_alltoallv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpialltoallv
!****************************************************************
 
! Subprogram not used    subroutine mpialltoallint (sendbuf, sendcnt, recvbuf, recvcnt, &
! Subprogram not used                               comm)
! Subprogram not used !
! Subprogram not used ! All-to-all scatter/gather
! Subprogram not used ! 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in)  :: sendbuf(*)
! Subprogram not used    integer, intent(in)  :: sendcnt
! Subprogram not used    integer, intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in)  :: recvcnt
! Subprogram not used    integer, intent(in)  :: comm
! Subprogram not used  
! Subprogram not used    integer :: ier              ! MPI error code
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_alltoall (sendbuf, sendcnt, mpiint, &
! Subprogram not used                       recvbuf, recvcnt, mpiint, &
! Subprogram not used                       comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_alltoallint failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpialltoallint

!****************************************************************
 
! Subprogram not used    subroutine mpiallgatherv (sendbuf, sendcnt, sendtype, &
! Subprogram not used                              recvbuf, recvcnts, rdispls, recvtype, &
! Subprogram not used                              comm)
! Subprogram not used !
! Subprogram not used ! Collect data from each task and broadcast resulting
! Subprogram not used ! vector to all tasks
! Subprogram not used ! 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real (r8), intent(in)  :: sendbuf(*)
! Subprogram not used    real (r8), intent(out) :: recvbuf(*)
! Subprogram not used    integer, intent(in) :: sendcnt
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: recvcnts(*)
! Subprogram not used    integer, intent(in) :: rdispls(*)
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used  
! Subprogram not used    integer ier   !MP error code
! Subprogram not used  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call mpi_allgatherv (sendbuf, sendcnt, sendtype, &
! Subprogram not used                         recvbuf, recvcnts, rdispls, recvtype, &
! Subprogram not used                         comm, ier)
! Subprogram not used    if (ier/=mpi_success) then
! Subprogram not used       write(iulog,*)'mpi_allgatherv failed ier=',ier
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used  
! Subprogram not used    return
! Subprogram not used    end subroutine mpiallgatherv
!****************************************************************
 
   subroutine mpiallgatherint (sendbuf, scount, recvbuf, rcount, comm)
!
! Collects integer data from each task and broadcasts resulting
! vector to all tasks
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog




   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: scount
   integer, intent(in)  :: rcount
   integer, intent(in)  :: comm
 
   integer ier   !MP error code




   call mpi_allgather (sendbuf, scount, mpiint, recvbuf, rcount, &
                       mpiint, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_allgather failed ier=',ier
      call endrun
   end if



 
   return
   end subroutine mpiallgatherint

!****************************************************************

! Subprogram not used    subroutine mpiwincreate(base,size,comm,win)
! Subprogram not used !
! Subprogram not used ! Creates window for MPI2 one-sided commands
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real(r8), intent(in)  :: base(*)
! Subprogram not used    integer,  intent(in)  :: size
! Subprogram not used    integer,  intent(in)  :: comm
! Subprogram not used    integer,  intent(out) :: win
! Subprogram not used !
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine mpiwincreate
!****************************************************************
!
! If 1 is not turned on
!

