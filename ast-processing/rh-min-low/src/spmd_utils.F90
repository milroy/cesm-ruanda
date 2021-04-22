module spmd_utils

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for miscellaneous 1 utilities
!          and information that are shared between dynamics and 
!          physics packages.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, December 2003
!   swap routines:      P. Worley
!   fc routines:        P. Worley
!   SMP node id logic:  P. Worley
!
! $Id$
! 
!-----------------------------------------------------------------------

!
! Performance bug work around for Gemini interconnect
!





!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils,   only: endrun

   use mpishorthand, only: mpiint, mpii8, mpichar, mpilog, mpipk,      &
                           mpic16, mpir8, mpir4, mpicom, mpimax

   use cam_logfile,  only: iulog

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   include 'mpif.h'          
   private                   ! Make the default access private
   save
!
! Forward from mpishorthand.F with the idea of phasing out use of and removing that file
!





!
!  Forward these from mpif.h (or mpi.mod), the idea being that this should
!  be the only module that uses mpi directly, the rest of cam should use spmd_utils
!
   public :: mpi_max_processor_name,                     &
             mpi_integer, mpi_integer8, mpi_character,   &
             mpi_logical, mpi_real8, mpi_real4,          &
             mpi_complex16,                              &
             mpi_packed, mpi_max,                        &
             mpi_comm_null, mpi_group_null,              &
             mpi_undefined, mpi_status_size, mpi_success,&
             mpi_status_ignore, mpi_double_precision, mpi_sum, mpir8





!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public pair      ! $$$here...  originally from eul|sld/spmd_dyn
   public ceil2     ! $$$here...  originally from eul|sld/spmd_dyn
   public spmdinit
   public spmd_utils_readnl

   public swapm
   public fc_gatherv
   public fc_gathervr4
   public fc_gathervint
   public fc_gathervc
   public altalltoallv


!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------
! physics-motivated dynamics decomposition request
   logical, parameter :: def_mirror = .false.                 ! default
   logical, public    :: phys_mirror_decomp_req = def_mirror 
                         ! flag indicating whether latitudes and their
                         ! reflections across the equator should be 
                         ! assigned to consecutive processes


   public :: mpicom
   public :: mpichar




   logical, public              :: masterproc
   integer, public              :: masterprocid
   integer, public              :: iam
   integer, public              :: npes
   integer, public              :: nsmps
   integer, allocatable, public :: proc_smp_map(:)
   integer, parameter           :: DEFAULT_MASTERPROC=0 
                                      ! the value of iam which is assigned 
                                      ! the masterproc duties

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------
! Swap communication protocol options (reduced set):
!  3, 5:                  nonblocking send
!  2, 3, 4, 5:            nonblocking receive
!  4, 5:                  ready send
   integer, private, parameter :: min_comm_protocol =  2
   integer, private, parameter :: max_comm_protocol =  5
   integer, private, parameter :: def_comm_protocol =  4        ! default
   integer, public :: swap_comm_protocol = def_comm_protocol

! Swap communication maximum request count:
! = -1,0: do not limit number of outstanding send/receive requests
!    > 0: do not allow more than swap_comm_maxreq outstanding
!         nonblocking send requests or nonblocking receive requests
   integer, private, parameter :: def_comm_maxreq = 128        ! default
   integer, public :: swap_comm_maxreq = def_comm_maxreq

! Flow-controlled gather option:
!   < 0: use MPI_Gather
!  >= 0: use point-to-point with handshaking messages and 
!        preposting receive requests up to 
!         min(max(1,fc_gather_flow_cntl),max_gather_block_size) 
!        ahead
   integer, private, parameter :: max_gather_block_size = 64 ! max and default
   integer, public :: fc_gather_flow_cntl = max_gather_block_size

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!========================================================================

! Subprogram not used    integer function pair(np,p,k)
! Subprogram not used 
! Subprogram not used       integer np,p,k,q
! Subprogram not used       q = ieor(p,k)
! Subprogram not used       if(q.gt.np-1) then
! Subprogram not used          pair = -1
! Subprogram not used       else
! Subprogram not used          pair = q
! Subprogram not used       endif
! Subprogram not used       return
! Subprogram not used 
! Subprogram not used    end function pair

!========================================================================

! Subprogram not used   integer function ceil2(n)
! Subprogram not used      integer n,p
! Subprogram not used      p=1
! Subprogram not used      do while(p.lt.n)
! Subprogram not used         p=p*2
! Subprogram not used      enddo
! Subprogram not used      ceil2=p
! Subprogram not used      return
! Subprogram not used   end function ceil2

!========================================================================
  
  subroutine spmdinit( mpicom_atm )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: MPI initialization routine:  
    ! 
    ! Method: get number of cpus, processes, tids, etc
    !         dynamics and physics decompositions are set up later
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------

    implicit none
    integer, intent(in) :: mpicom_atm


    !
    ! Local workspace
    !
    integer i,j,c             ! indices
    integer npthreads         ! thread status
    integer ier               ! return error status    
    integer length            ! length of name
    integer max_len           ! maximum name length
    integer, allocatable :: lengths(:)! max lengths of names for use in gatherv
    integer, allocatable :: displs(:) ! offsets for use in gatherv
    logical done
    character, allocatable                             :: proc_name(:)  ! processor name, this task
    character, allocatable                             :: proc_names(:) ! processor names, all tasks
    character(len=mpi_max_processor_name)              :: tmp_name      ! temporary storage
    character(len=mpi_max_processor_name), allocatable :: smp_names(:)  ! SMP name
    logical mpi_running       ! returned value indicates if MPI_INIT has been called

    !---------------------------------------------------------------------------
    !
    ! Determine 1 MPI communicator group
    !
    mpicom  = mpicom_atm
    !
    ! Set mpishorthand variables.  Need to set as variables rather than parameters since
    ! some MPI implementations set values for MPI tags at run time
    !
    mpiint  = mpi_integer
    mpii8   = mpi_integer8
    mpichar = mpi_character
    mpilog  = mpi_logical
    mpir4   = mpi_real4
    mpir8   = mpi_real8
    mpic16  = mpi_complex16
    mpipk   = mpi_packed
    mpimax  = mpi_max
    !
    ! Get my id  
    !
    call mpi_comm_rank (mpicom, iam, ier) 
    masterprocid = DEFAULT_MASTERPROC
    if (iam == DEFAULT_MASTERPROC) then 
       masterproc = .true.
    else
       masterproc = .false.
    end if
    !
    ! Get number of processors
    !
    max_len = mpi_max_processor_name
    call mpi_comm_size (mpicom, npes, ier)
    allocate ( displs(npes) )
    allocate ( lengths(npes) )
    allocate ( proc_name(max_len) )
    allocate ( proc_names(max_len*npes) )
 
    !
    ! Get processor names and send to root. 
    !
    call mpi_get_processor_name (tmp_name, length, ier)
    proc_name(:) = ' '
    do i = 1, length
       proc_name(i) = tmp_name(i:i)
    end do

    proc_names(:) = ' '
    lengths(:) = max_len
    do i=1,npes
       displs(i) = (i-1)*max_len
    enddo
    call fc_gathervc (proc_name,  max_len, mpichar, &
                      proc_names, lengths, displs, mpichar, &
                      0, mpicom, flow_cntl=-1)
    if (masterproc) then
       write(iulog,*) npes, 'pes participating in computation'
       write(iulog,*) '-----------------------------------'
       write(iulog,*) 'TASK#  NAME'
       do i=0,min(npes-1,256)  ! dont print too many of these
          do c=1,max_len
             tmp_name(c:c) = proc_names(i*max_len+c)
          enddo
          write(iulog,'(i3,2x,a)') i,trim(tmp_name)
       end do
       if(npes-1>256) then
          write(iulog,*) '... list truncated at 256'
       end if
    end if
    !
    ! Identify SMP nodes and process/SMP mapping.
    ! (Assume that processor names are SMP node names on SMP clusters.)
    !
    allocate ( proc_smp_map(0:npes-1) )
    if (masterproc) then
       allocate ( smp_names(0:npes-1) )
       smp_names(:) = ' '
       proc_smp_map(:) = -1
       !
       nsmps = 1
       do c=1,max_len
          tmp_name(c:c) = proc_names(c)
       enddo
       smp_names(0) = trim(tmp_name)
       proc_smp_map(0) = 0
       !
       do i=1,npes-1
          do c=1,max_len
             tmp_name(c:c) = proc_names(i*max_len+c)
          enddo

          j = 0
          done = .false.
          do while ((.not. done) .and. (j < nsmps))
             if (smp_names(j) .eq. trim(tmp_name)) then
                proc_smp_map(i) = j
                done = .true.
             endif
             j = j + 1
          enddo

          if (.not. done) then
             smp_names(nsmps) = trim(tmp_name)
             proc_smp_map(i) = nsmps
             nsmps = nsmps + 1
          endif

       enddo
       deallocate(smp_names)
    endif
    call mpibcast(nsmps, 1, mpiint, 0, mpicom)
    call mpibcast(proc_smp_map, npes, mpiint, 0, mpicom)
    !
    deallocate(displs)
    deallocate(lengths)
    deallocate(proc_name)
    deallocate(proc_names)


  end subroutine spmdinit

!
!========================================================================
!
! Subprogram not used    subroutine swapm (steps, nprocs, swapids,               &
! Subprogram not used                      sndbuf, sbuf_siz, sndlths, sdispls,   &
! Subprogram not used                      rcvbuf, rbuf_siz, rcvlths, rdispls,   &
! Subprogram not used                      comm, comm_protocol, comm_maxreq      )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used !   Reduced version of original swapm (for swap of multiple messages 
! Subprogram not used !   using MPI point-to-point routines), more efficiently implementing a 
! Subprogram not used !   subset of the swap protocols.
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! comm_protocol:
! Subprogram not used !  = 3 or 5: use nonblocking send
! Subprogram not used !  = 2 or 4: use blocking send
! Subprogram not used !  = 4 or 5: use handshaking protocol
! Subprogram not used ! comm_maxreq:
! Subprogram not used !  =-1,0: do not limit number of outstanding send/receive requests
! Subprogram not used !     >0: do not allow more than min(comm_maxreq, steps) outstanding
! Subprogram not used !         nonblocking send requests or nonblocking receive requests
! Subprogram not used !
! Subprogram not used ! Author of original version:  P. Worley
! Subprogram not used ! Ported to 1: P. Worley, December 2003
! Subprogram not used ! Simplified version: P. Worley, October, 2008
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    implicit none
! Subprogram not used !---------------------------Input arguments--------------------------
! Subprogram not used !
! Subprogram not used    integer, intent(in)   :: steps              ! number of swaps to initiate
! Subprogram not used    integer, intent(in)   :: nprocs             ! size of communicator
! Subprogram not used    integer, intent(in)   :: sbuf_siz           ! size of send buffer
! Subprogram not used    integer, intent(in)   :: rbuf_siz           ! size of receive buffer
! Subprogram not used    integer, intent(in)   :: swapids(steps)     ! MPI process id of swap partners
! Subprogram not used 
! Subprogram not used    integer, intent(in)   :: sndlths(0:nprocs-1)! length of outgoing message
! Subprogram not used    integer, intent(in)   :: sdispls(0:nprocs-1)! offset from beginning of send
! Subprogram not used                                                !  buffer where outgoing messages
! Subprogram not used                                                !  should be sent from
! Subprogram not used    integer, intent(in)   :: rcvlths(0:nprocs-1)! length of incoming messages
! Subprogram not used    integer, intent(in)   :: rdispls(0:nprocs-1)! offset from beginning of receive 
! Subprogram not used                                                !  buffer where incoming messages
! Subprogram not used                                                !  should be placed
! Subprogram not used    real(r8), intent(in)  :: sndbuf(sbuf_siz)   ! outgoing message buffer
! Subprogram not used    real(r8), intent(out) :: rcvbuf(rbuf_siz)   ! incoming message buffer
! Subprogram not used 
! Subprogram not used    integer, intent(in)   :: comm               ! MPI communicator
! Subprogram not used    integer, intent(in)   :: comm_protocol      ! swap_comm protocol
! Subprogram not used    integer, intent(in)   :: comm_maxreq        ! maximum number of outstanding 
! Subprogram not used                                                !  nonblocking requests
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    integer :: p                                ! process index
! Subprogram not used    integer :: istep                            ! loop index
! Subprogram not used    integer :: offset_s                         ! index of message beginning in 
! Subprogram not used                                                !  send buffer
! Subprogram not used    integer :: offset_r                         ! index of message beginning in 
! Subprogram not used                                                !  receive buffer
! Subprogram not used    integer :: sndids(steps)                    ! send request ids
! Subprogram not used    integer :: rcvids(steps)                    ! receive request ids
! Subprogram not used    integer :: hs_rcvids(steps)                 ! handshake receive request ids
! Subprogram not used 
! Subprogram not used    integer :: maxreq, maxreqh                  ! maximum number of outstanding 
! Subprogram not used                                                !  nonblocking requests (and half)
! Subprogram not used    integer :: hs_s, hs_r(steps)                ! handshake variables (send/receive)
! Subprogram not used    integer :: rstep                            ! "receive" step index
! Subprogram not used 
! Subprogram not used    logical :: handshake, sendd                 ! protocol option flags
! Subprogram not used 
! Subprogram not used    integer :: ier                              ! return error status    
! Subprogram not used    integer :: status(MPI_STATUS_SIZE)          ! MPI status 
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    if (steps .eq. 0) return
! Subprogram not used 
! Subprogram not used    ! identify communication protocol
! Subprogram not used    if ((comm_protocol < 2) .or. (comm_protocol > 5)) then
! Subprogram not used       sendd = .true.
! Subprogram not used       handshake = .true.
! Subprogram not used    else
! Subprogram not used       if ((comm_protocol .eq. 4) .or. (comm_protocol .eq. 5)) then
! Subprogram not used          handshake = .true.
! Subprogram not used       else
! Subprogram not used          handshake = .false.
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if ((comm_protocol .eq. 2) .or. (comm_protocol .eq. 4)) then
! Subprogram not used          sendd = .true.
! Subprogram not used       else
! Subprogram not used          sendd = .false.
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    ! identify maximum number of outstanding nonblocking requests to permit
! Subprogram not used    if (steps .eq. 1) then
! Subprogram not used       maxreq  = 1
! Subprogram not used       maxreqh = 1
! Subprogram not used    else
! Subprogram not used       if (comm_maxreq >= -1) then
! Subprogram not used          maxreq = comm_maxreq
! Subprogram not used       else
! Subprogram not used          maxreq = steps
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if ((maxreq .le. steps) .and. (maxreq > 0)) then
! Subprogram not used          if (maxreq > 1) then
! Subprogram not used             maxreqh = maxreq/2
! Subprogram not used          else
! Subprogram not used             maxreq  = 2
! Subprogram not used             maxreqh = 1
! Subprogram not used          endif
! Subprogram not used       else
! Subprogram not used          maxreq  = steps
! Subprogram not used          maxreqh = steps
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used ! Four protocol options:
! Subprogram not used !  (1) handshaking + blocking sends
! Subprogram not used    if ((handshake) .and. (sendd)) then
! Subprogram not used 
! Subprogram not used       ! Initialize handshake variable
! Subprogram not used       hs_s = 1
! Subprogram not used 
! Subprogram not used       ! Post initial handshake receive requests
! Subprogram not used       do istep=1,maxreq
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             call mpi_irecv( hs_r(istep), 1, mpiint, p, iam, comm, &
! Subprogram not used                             hs_rcvids(istep), ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! Post initial receive requests
! Subprogram not used       do istep=1,maxreq
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             offset_r = rdispls(p)+1
! Subprogram not used             call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                             comm, rcvids(istep), ier )
! Subprogram not used             call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       rstep = maxreq
! Subprogram not used 
! Subprogram not used       ! Send (and start receiving) data 
! Subprogram not used       do istep=1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used 
! Subprogram not used          ! Submit new rsend request
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             offset_s = sdispls(p)+1
! Subprogram not used             call mpi_wait  ( hs_rcvids(istep), MPI_STATUS_IGNORE, ier )
! Subprogram not used             call mpi_rsend ( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
! Subprogram not used                              comm, ier )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (istep > maxreqh) then
! Subprogram not used 
! Subprogram not used             ! Wait for oldest irecv request to complete
! Subprogram not used             p = swapids(istep-maxreqh)
! Subprogram not used             if (rcvlths(p) > 0) then
! Subprogram not used                call mpi_wait( rcvids(istep-maxreqh), status, ier )
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             if (rstep < steps) then
! Subprogram not used                rstep = rstep + 1
! Subprogram not used                p = swapids(rstep)
! Subprogram not used 
! Subprogram not used                ! Submit a new handshake irecv request
! Subprogram not used                if (sndlths(p) > 0) then
! Subprogram not used                   call mpi_irecv( hs_r(rstep), 1, mpiint, p, iam, comm, &
! Subprogram not used                                   hs_rcvids(rstep), ier )
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                ! Submit a new irecv request
! Subprogram not used                if (rcvlths(p) > 0) then
! Subprogram not used                   offset_r = rdispls(p)+1
! Subprogram not used                   call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                                   comm, rcvids(rstep), ier )
! Subprogram not used                   call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used !
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! wait for rest of receive requests to complete
! Subprogram not used       do istep=steps-maxreqh+1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             call mpi_wait( rcvids(istep), status, ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used !  (2) handshaking + nonblocking sends
! Subprogram not used    elseif ((handshake) .and. (.not. sendd)) then
! Subprogram not used 
! Subprogram not used       ! Initialize handshake variable
! Subprogram not used       hs_s = 1
! Subprogram not used 
! Subprogram not used       ! Post initial handshake receive requests
! Subprogram not used       do istep=1,maxreq
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             call mpi_irecv( hs_r(istep), 1, mpiint, p, iam, comm, &
! Subprogram not used                             hs_rcvids(istep), ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! Post initial receive requests
! Subprogram not used       do istep=1,maxreq
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             offset_r = rdispls(p)+1
! Subprogram not used             call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                             comm, rcvids(istep), ier )
! Subprogram not used             call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       rstep = maxreq
! Subprogram not used 
! Subprogram not used       ! Send (and start receiving) data 
! Subprogram not used       do istep=1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used 
! Subprogram not used          ! Submit new irsend request
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             offset_s = sdispls(p)+1
! Subprogram not used             call mpi_wait  ( hs_rcvids(istep), MPI_STATUS_IGNORE, ier )
! Subprogram not used             call mpi_irsend( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
! Subprogram not used                              comm, sndids(istep), ier )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (istep > maxreqh) then
! Subprogram not used 
! Subprogram not used             ! Wait for oldest irecv request to complete
! Subprogram not used             p = swapids(istep-maxreqh)
! Subprogram not used             if (rcvlths(p) > 0) then
! Subprogram not used                call mpi_wait( rcvids(istep-maxreqh), status, ier )
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             if (rstep < steps) then
! Subprogram not used                rstep = rstep + 1
! Subprogram not used                p = swapids(rstep)
! Subprogram not used 
! Subprogram not used                ! Submit a new handshake irecv request
! Subprogram not used                if (sndlths(p) > 0) then
! Subprogram not used                   call mpi_irecv( hs_r(rstep), 1, mpiint, p, iam, comm, &
! Subprogram not used                                   hs_rcvids(rstep), ier )
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                ! Submit a new irecv request
! Subprogram not used                if (rcvlths(p) > 0) then
! Subprogram not used                   offset_r = rdispls(p)+1
! Subprogram not used                   call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                                   comm, rcvids(rstep), ier )
! Subprogram not used                   call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             ! Wait for outstanding i(r)send request to complete
! Subprogram not used             p = swapids(istep-maxreqh)
! Subprogram not used             if (sndlths(p) > 0) then
! Subprogram not used                call mpi_wait( sndids(istep-maxreqh), status, ier )
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! wait for rest of send and receive requests to complete
! Subprogram not used       do istep=steps-maxreqh+1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             call mpi_wait( rcvids(istep), status, ier )
! Subprogram not used          endif
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             call mpi_wait( sndids(istep), status, ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used !  (3) no handshaking + blocking sends
! Subprogram not used    elseif ((.not. handshake) .and. (sendd)) then
! Subprogram not used 
! Subprogram not used       ! Post receive requests
! Subprogram not used       do istep=1,maxreq
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             offset_r = rdispls(p)+1
! Subprogram not used             call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                             comm, rcvids(istep), ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       rstep = maxreq
! Subprogram not used 
! Subprogram not used       ! Send (and start receiving) data 
! Subprogram not used       do istep=1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used 
! Subprogram not used          ! Submit new send request
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             offset_s = sdispls(p)+1
! Subprogram not used             call mpi_send( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
! Subprogram not used                            comm, ier )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (istep > maxreqh) then
! Subprogram not used 
! Subprogram not used             ! Wait for oldest irecv request to complete
! Subprogram not used             p = swapids(istep-maxreqh)
! Subprogram not used             if (rcvlths(p) > 0) then
! Subprogram not used                call mpi_wait( rcvids(istep-maxreqh), status, ier )
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             ! Submit a new irecv request
! Subprogram not used             if (rstep < steps) then
! Subprogram not used                rstep = rstep + 1
! Subprogram not used                p = swapids(rstep)
! Subprogram not used                if (rcvlths(p) > 0) then
! Subprogram not used                   offset_r = rdispls(p)+1
! Subprogram not used                   call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                                   comm, rcvids(rstep), ier )
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! wait for rest of send and receive requests to complete
! Subprogram not used       do istep=steps-maxreqh+1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             call mpi_wait( rcvids(istep), status, ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used !  (4) no handshaking + nonblocking sends
! Subprogram not used    elseif ((.not. handshake) .and. (.not. sendd)) then
! Subprogram not used 
! Subprogram not used       ! Post receive requests
! Subprogram not used       do istep=1,maxreq
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             offset_r = rdispls(p)+1
! Subprogram not used             call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                             comm, rcvids(istep), ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       rstep = maxreq
! Subprogram not used 
! Subprogram not used       ! Send (and start receiving) data 
! Subprogram not used       do istep=1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used 
! Subprogram not used          ! Submit new isend request
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             offset_s = sdispls(p)+1
! Subprogram not used             call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
! Subprogram not used                             comm, sndids(istep), ier )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (istep > maxreqh) then
! Subprogram not used 
! Subprogram not used             ! Wait for oldest irecv request to complete
! Subprogram not used             p = swapids(istep-maxreqh)
! Subprogram not used             if (rcvlths(p) > 0) then
! Subprogram not used                call mpi_wait( rcvids(istep-maxreqh), status, ier )
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             ! Submit a new irecv request
! Subprogram not used             if (rstep < steps) then
! Subprogram not used                rstep = rstep + 1
! Subprogram not used                p = swapids(rstep)
! Subprogram not used                if (rcvlths(p) > 0) then
! Subprogram not used                   offset_r = rdispls(p)+1
! Subprogram not used                   call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
! Subprogram not used                                   comm, rcvids(rstep), ier )
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             ! Wait for outstanding i(r)send request to complete
! Subprogram not used             p = swapids(istep-maxreqh)
! Subprogram not used             if (sndlths(p) > 0) then
! Subprogram not used                call mpi_wait( sndids(istep-maxreqh), status, ier )
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! wait for rest of send and receive requests to complete
! Subprogram not used       do istep=steps-maxreqh+1,steps
! Subprogram not used          p = swapids(istep)
! Subprogram not used          if (rcvlths(p) > 0) then
! Subprogram not used             call mpi_wait( rcvids(istep), status, ier )
! Subprogram not used          endif
! Subprogram not used          if (sndlths(p) > 0) then
! Subprogram not used             call mpi_wait( sndids(istep), status, ier )
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used 
! Subprogram not used    end subroutine swapm
!
!========================================================================

!----------------------------------------------------------------------- 
! 
! Purpose: gather collective with additional flow control, so as to 
!          be more robust when used with high process counts. 
!          If flow_cntl optional parameter 
!           < 0: use MPI_Gather
!          >= 0: use point-to-point with handshaking messages and 
!                preposting receive requests up to 
!                 min(max(1,flow_cntl),max_gather_block_size) 
!                ahead if optional flow_cntl parameter is present.
!                Otherwise, fc_gather_flow_cntl is used in its place.
!          Default value is 64.
! 
! Entry points:
!      fc_gatherv       functionally equivalent to mpi_gatherv
!      fc_gathervr4     functionally equivalent to mpi_gatherv for real*4 data
!      fc_gathervint    functionally equivalent to mpi_gatherv for integer data
!      fc_gathervc      functionally equivalent to mpi_gatherv for character data
!
! Author: P. Worley
!-----------------------------------------------------------------------

!
!========================================================================
!
! Subprogram not used    subroutine fc_gatherv (sendbuf, sendcnt, sendtype, &
! Subprogram not used                          recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                          root, comm, flow_cntl )
! Subprogram not used !
! Subprogram not used ! Collects different messages from each process on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils, only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
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
! Subprogram not used    integer, optional, intent(in) :: flow_cntl
! Subprogram not used 
! Subprogram not used    real (r8) :: signal
! Subprogram not used    logical fc_gather         ! use explicit flow control?
! Subprogram not used    integer gather_block_size ! number of preposted receive requests
! Subprogram not used 
! Subprogram not used    integer :: mytid, mysize, mtag, p, q, i, count
! Subprogram not used    integer :: preposts, head, tail
! Subprogram not used    integer :: rcvid(max_gather_block_size)
! Subprogram not used    integer :: status(MPI_STATUS_SIZE)
! Subprogram not used    integer ier               ! MPI error code
! Subprogram not used 
! Subprogram not used    if ( present(flow_cntl) ) then
! Subprogram not used       if (flow_cntl >= 0) then
! Subprogram not used          gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
! Subprogram not used          fc_gather = .true.
! Subprogram not used       else
! Subprogram not used          fc_gather = .false.
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       if (fc_gather_flow_cntl >= 0) then
! Subprogram not used          gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
! Subprogram not used          fc_gather = .true.
! Subprogram not used       else
! Subprogram not used          fc_gather = .false.
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (fc_gather) then
! Subprogram not used  
! Subprogram not used       call mpi_comm_rank (comm, mytid, ier)
! Subprogram not used       call mpi_comm_size (comm, mysize, ier)
! Subprogram not used       mtag = 0
! Subprogram not used       if (root .eq. mytid) then
! Subprogram not used 
! Subprogram not used ! prepost gather_block_size irecvs, and start receiving data
! Subprogram not used          preposts = min(mysize-1, gather_block_size)
! Subprogram not used          head = 0
! Subprogram not used          count = 0
! Subprogram not used          do p=0, mysize-1
! Subprogram not used             if (p .ne. root) then
! Subprogram not used                q = p+1
! Subprogram not used                if (recvcnts(q) > 0) then
! Subprogram not used                   count = count + 1
! Subprogram not used                   if (count > preposts) then
! Subprogram not used                      tail = mod(head,preposts) + 1
! Subprogram not used                      call mpi_wait (rcvid(tail), status, ier)
! Subprogram not used                   end if
! Subprogram not used                   head = mod(head,preposts) + 1
! Subprogram not used                   call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
! Subprogram not used                                    recvtype, p, mtag, comm, rcvid(head), &
! Subprogram not used                                    ier )
! Subprogram not used                   call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
! Subprogram not used                end if
! Subprogram not used             end if
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used ! copy local data
! Subprogram not used          q = mytid+1
! Subprogram not used          do i=1,sendcnt
! Subprogram not used             recvbuf(displs(q)+i) = sendbuf(i)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used ! wait for final data
! Subprogram not used          do i=1,min(count,preposts)
! Subprogram not used             call mpi_wait (rcvid(i), status, ier)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       else
! Subprogram not used 
! Subprogram not used          if (sendcnt > 0) then
! Subprogram not used             call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
! Subprogram not used                             status, ier )
! Subprogram not used             call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
! Subprogram not used                              comm, ier )
! Subprogram not used          end if
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used       if (ier /= mpi_success) then
! Subprogram not used          write(iulog,*)'fc_gatherv_r8 failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used  
! Subprogram not used       call mpi_gatherv (sendbuf, sendcnt, sendtype, &
! Subprogram not used                         recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                         root, comm, ier)
! Subprogram not used       if (ier /= mpi_success) then
! Subprogram not used          write(iulog,*)'mpi_gatherv failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine fc_gatherv
!
!========================================================================
!
! Subprogram not used    subroutine fc_gathervr4 (sendbuf, sendcnt, sendtype, &
! Subprogram not used                            recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                            root, comm, flow_cntl )
! Subprogram not used !
! Subprogram not used ! Collects different messages from each process on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r4 => shr_kind_r4, r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils, only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
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
! Subprogram not used    integer, optional, intent(in) :: flow_cntl
! Subprogram not used 
! Subprogram not used    real (r8) :: signal
! Subprogram not used    logical fc_gather         ! use explicit flow control?
! Subprogram not used    integer gather_block_size ! number of preposted receive requests
! Subprogram not used 
! Subprogram not used    integer :: mytid, mysize, mtag, p, q, i, count
! Subprogram not used    integer :: preposts, head, tail
! Subprogram not used    integer :: rcvid(max_gather_block_size)
! Subprogram not used    integer :: status(MPI_STATUS_SIZE)
! Subprogram not used    integer ier               ! MPI error code
! Subprogram not used 
! Subprogram not used    if ( present(flow_cntl) ) then
! Subprogram not used       if (flow_cntl >= 0) then
! Subprogram not used          gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
! Subprogram not used          fc_gather = .true.
! Subprogram not used       else
! Subprogram not used          fc_gather = .false.
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       if (fc_gather_flow_cntl >= 0) then
! Subprogram not used          gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
! Subprogram not used          fc_gather = .true.
! Subprogram not used       else
! Subprogram not used          fc_gather = .false.
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (fc_gather) then
! Subprogram not used  
! Subprogram not used       call mpi_comm_rank (comm, mytid, ier)
! Subprogram not used       call mpi_comm_size (comm, mysize, ier)
! Subprogram not used       mtag = 0
! Subprogram not used       if (root .eq. mytid) then
! Subprogram not used 
! Subprogram not used ! prepost gather_block_size irecvs, and start receiving data
! Subprogram not used          preposts = min(mysize-1, gather_block_size)
! Subprogram not used          head = 0
! Subprogram not used          count = 0
! Subprogram not used          do p=0, mysize-1
! Subprogram not used             if (p .ne. root) then
! Subprogram not used                q = p+1
! Subprogram not used                if (recvcnts(q) > 0) then
! Subprogram not used                   count = count + 1
! Subprogram not used                   if (count > preposts) then
! Subprogram not used                      tail = mod(head,preposts) + 1
! Subprogram not used                      call mpi_wait (rcvid(tail), status, ier)
! Subprogram not used                   end if
! Subprogram not used                   head = mod(head,preposts) + 1
! Subprogram not used                   call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
! Subprogram not used                                    recvtype, p, mtag, comm, rcvid(head), &
! Subprogram not used                                    ier )
! Subprogram not used                   call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
! Subprogram not used                end if
! Subprogram not used             end if
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used ! copy local data
! Subprogram not used          q = mytid+1
! Subprogram not used          do i=1,sendcnt
! Subprogram not used             recvbuf(displs(q)+i) = sendbuf(i)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used ! wait for final data
! Subprogram not used          do i=1,min(count,preposts)
! Subprogram not used             call mpi_wait (rcvid(i), status, ier)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       else
! Subprogram not used 
! Subprogram not used          if (sendcnt > 0) then
! Subprogram not used             call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
! Subprogram not used                             status, ier )
! Subprogram not used             call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
! Subprogram not used                              comm, ier )
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used       if (ier /= mpi_success) then
! Subprogram not used          write(iulog,*)'fc_gatherv_r4 failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used  
! Subprogram not used       call mpi_gatherv (sendbuf, sendcnt, sendtype, &
! Subprogram not used                         recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                         root, comm, ier)
! Subprogram not used       if (ier /= mpi_success) then
! Subprogram not used          write(iulog,*)'mpi_gatherv failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine fc_gathervr4
!
!========================================================================
!
! Subprogram not used    subroutine fc_gathervint (sendbuf, sendcnt, sendtype, &
! Subprogram not used                             recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                             root, comm, flow_cntl )
! Subprogram not used !
! Subprogram not used ! Collects different messages from each process on masterproc
! Subprogram not used !
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use mpishorthand
! Subprogram not used    use abortutils, only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
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
! Subprogram not used    integer, optional, intent(in) :: flow_cntl
! Subprogram not used 
! Subprogram not used    real (r8) :: signal
! Subprogram not used    logical fc_gather         ! use explicit flow control?
! Subprogram not used    integer gather_block_size ! number of preposted receive requests
! Subprogram not used 
! Subprogram not used    integer :: mytid, mysize, mtag, p, q, i, count
! Subprogram not used    integer :: preposts, head, tail
! Subprogram not used    integer :: rcvid(max_gather_block_size)
! Subprogram not used    integer :: status(MPI_STATUS_SIZE)
! Subprogram not used    integer ier               ! MPI error code
! Subprogram not used 
! Subprogram not used    if ( present(flow_cntl) ) then
! Subprogram not used       if (flow_cntl >= 0) then
! Subprogram not used          gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
! Subprogram not used          fc_gather = .true.
! Subprogram not used       else
! Subprogram not used          fc_gather = .false.
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       if (fc_gather_flow_cntl >= 0) then
! Subprogram not used          gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
! Subprogram not used          fc_gather = .true.
! Subprogram not used       else
! Subprogram not used          fc_gather = .false.
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (fc_gather) then
! Subprogram not used  
! Subprogram not used       call mpi_comm_rank (comm, mytid, ier)
! Subprogram not used       call mpi_comm_size (comm, mysize, ier)
! Subprogram not used       mtag = 0
! Subprogram not used       if (root .eq. mytid) then
! Subprogram not used 
! Subprogram not used ! prepost gather_block_size irecvs, and start receiving data
! Subprogram not used          preposts = min(mysize-1, gather_block_size)
! Subprogram not used          head = 0
! Subprogram not used          count = 0
! Subprogram not used          do p=0, mysize-1
! Subprogram not used             if (p .ne. root) then
! Subprogram not used                q = p+1
! Subprogram not used                if (recvcnts(q) > 0) then
! Subprogram not used                   count = count + 1
! Subprogram not used                   if (count > preposts) then
! Subprogram not used                      tail = mod(head,preposts) + 1
! Subprogram not used                      call mpi_wait (rcvid(tail), status, ier)
! Subprogram not used                   end if
! Subprogram not used                   head = mod(head,preposts) + 1
! Subprogram not used                   call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
! Subprogram not used                                    recvtype, p, mtag, comm, rcvid(head), &
! Subprogram not used                                    ier )
! Subprogram not used                   call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
! Subprogram not used                end if
! Subprogram not used             end if
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used ! copy local data
! Subprogram not used          q = mytid+1
! Subprogram not used          do i=1,sendcnt
! Subprogram not used             recvbuf(displs(q)+i) = sendbuf(i)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used ! wait for final data
! Subprogram not used          do i=1,min(count,preposts)
! Subprogram not used             call mpi_wait (rcvid(i), status, ier)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       else
! Subprogram not used 
! Subprogram not used          if (sendcnt > 0) then
! Subprogram not used             call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
! Subprogram not used                             status, ier )
! Subprogram not used             call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
! Subprogram not used                              comm, ier )
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used       if (ier /= mpi_success) then
! Subprogram not used          write(iulog,*)'fc_gatherv_int failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used  
! Subprogram not used       call mpi_gatherv (sendbuf, sendcnt, sendtype, &
! Subprogram not used                         recvbuf, recvcnts, displs, recvtype, &
! Subprogram not used                         root, comm, ier)
! Subprogram not used       if (ier /= mpi_success) then
! Subprogram not used          write(iulog,*)'mpi_gatherv failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine fc_gathervint
!
!========================================================================
!
   subroutine fc_gathervc (sendbuf, sendcnt, sendtype, &
                           recvbuf, recvcnts, displs, recvtype, &
                           root, comm, flow_cntl )
!
! Collects different messages from each process on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun
   use cam_logfile,  only: iulog

   implicit none

   character, intent(in)  :: sendbuf(*)
   character, intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
   integer, optional, intent(in) :: flow_cntl

   real (r8) :: signal
   logical fc_gather         ! use explicit flow control?
   integer gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer ier               ! MPI error code

   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      if (fc_gather_flow_cntl >= 0) then
         gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
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
                  call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
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
            call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
                            status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
          end if

      endif
      if (ier /= mpi_success) then
         write(iulog,*)'fc_gatherv_char failed ier=',ier
         call endrun
      end if

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)
      if (ier /= mpi_success) then
         write(iulog,*)'mpi_gatherv failed ier=',ier
         call endrun
      end if

   endif

   return
   end subroutine fc_gathervc
!
!========================================================================

!----------------------------------------------------------------------- 
! 
! Purpose: implementations of MPI_Alltoall using different messaging
!          layers and different communication protocols, controlled
!          by option argument:
!  0: use mpi_alltoallv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
!  otherwise use mpi_sendrecv implementation
! 
! Entry points:
!      altalltoallv
!
! Author: P. Worley
!-----------------------------------------------------------------------

!****************************************************************
! Subprogram not used    subroutine altalltoallv (option, mytid, nprocs, steps, dests, &
! Subprogram not used                  sendbuf, sbuf_siz, sendcnts, sdispls, sendtype, &
! Subprogram not used                  recvbuf, rbuf_siz, recvcnts, rdispls, recvtype, &
! Subprogram not used                  msgtag, pdispls, desttype, recvwin, comm)
! Subprogram not used !
! Subprogram not used ! All-to-all scatter/gather implemented using Co-Array
! Subprogram not used ! Fortran one-sided commands, MPI-2 one sided commands,
! Subprogram not used ! SWAP module MPI-1 commands, MPI_ALLTOALLV or MPI_SENDRECV.
! Subprogram not used !
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in) :: option               ! 0: mpi_alltoallv
! Subprogram not used                                                ! 1: swap package
! Subprogram not used                                                ! 2: mpi2 
! Subprogram not used                                                ! 3: co-array fortran
! Subprogram not used                                        ! otherwise: sendrecv
! Subprogram not used    integer, intent(in) :: mytid
! Subprogram not used    integer, intent(in) :: nprocs
! Subprogram not used    integer, intent(in) :: steps
! Subprogram not used    integer, intent(in) :: dests(steps)
! Subprogram not used    integer, intent(in) :: sbuf_siz
! Subprogram not used    integer, intent(in) :: sendcnts(0:nprocs-1)
! Subprogram not used    integer, intent(in) :: sdispls(0:nprocs-1)
! Subprogram not used    integer, intent(in) :: sendtype
! Subprogram not used    integer, intent(in) :: rbuf_siz
! Subprogram not used    integer, intent(in) :: recvcnts(0:nprocs-1)
! Subprogram not used    integer, intent(in) :: rdispls(0:nprocs-1)
! Subprogram not used    integer, intent(in) :: recvtype
! Subprogram not used    integer, intent(in) :: msgtag
! Subprogram not used    integer, intent(in) :: pdispls(0:nprocs-1)   ! displacement at 
! Subprogram not used                                                 !  destination
! Subprogram not used    integer, intent(in) :: desttype
! Subprogram not used    integer, intent(in) :: recvwin
! Subprogram not used    integer, intent(in) :: comm
! Subprogram not used 
! Subprogram not used    real (r8), intent(in)  :: sendbuf(sbuf_siz)
! Subprogram not used    real (r8), intent(out) :: recvbuf(rbuf_siz)
! Subprogram not used 
! Subprogram not used    integer :: loption          ! local copy of option
! Subprogram not used    integer :: dest             ! MPI remote process id
! Subprogram not used    integer :: ier              ! MPI error code
! Subprogram not used    integer :: i                ! loop index
! Subprogram not used    integer :: sndids(steps)    ! nonblocking MPI send request ids
! Subprogram not used    integer :: rcvids(steps)    ! nonblocking MPI recv request ids
! Subprogram not used    integer :: status(MPI_STATUS_SIZE)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    loption = option
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  using MPI library collective MPI_ALLTOALLV
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used    if (loption .eq. 0) then
! Subprogram not used 
! Subprogram not used       call mpi_alltoallv (sendbuf, sendcnts, sdispls, sendtype, &
! Subprogram not used                           recvbuf, recvcnts, rdispls, recvtype, &
! Subprogram not used                           comm, ier)
! Subprogram not used !
! Subprogram not used ! test for error
! Subprogram not used       if (ier/=mpi_success) then
! Subprogram not used          write(iulog,*)'altalltoallv (mpi_alltoallv) failed ier=',ier
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  Co-Array Fortran implementation of alltoallv
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used       if (loption .eq. 3) then
! Subprogram not used 
! Subprogram not used          loption = -1
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  MPI-2 one-sided implementation of alltoallv
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used       elseif (loption .eq. 2) then
! Subprogram not used          loption = -1
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  MPI-1 two-sided implementation of alltoallv
! Subprogram not used !  using SWAP routines
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used       elseif (loption .eq. 1) then
! Subprogram not used 
! Subprogram not used          call swapm(steps, nprocs, dests,                      &
! Subprogram not used                     sendbuf, sbuf_siz, sendcnts, sdispls,      &
! Subprogram not used                     recvbuf, rbuf_siz, recvcnts, rdispls,      &
! Subprogram not used                     comm, swap_comm_protocol, swap_comm_maxreq )
! Subprogram not used !
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  Anything else defined to be MPI_SENDRECV implementation
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used       else
! Subprogram not used !
! Subprogram not used          loption = -1
! Subprogram not used !
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  MPI_SENDRECV implementation of alltoallv
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used       if (loption .eq. -1) then
! Subprogram not used          do i=1, steps
! Subprogram not used             dest = dests(i)
! Subprogram not used             call mpi_sendrecv (sendbuf(sdispls(dest)+1), sendcnts(dest), &
! Subprogram not used                                sendtype, dest, msgtag,                   &
! Subprogram not used                                recvbuf(rdispls(dest)+1), recvcnts(dest), &
! Subprogram not used                                recvtype, dest, msgtag,                   &
! Subprogram not used                                comm, status, ier)
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! test for error
! Subprogram not used          if (ier/=mpi_success) then
! Subprogram not used             write(iulog,*)'altalltoallv (mpi1_alltoallv) failed ier=',ier
! Subprogram not used             call endrun
! Subprogram not used          end if
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  Local copy (if necessary)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used       if (sendcnts(mytid) > 0) then
! Subprogram not used          do i=1,sendcnts(iam)
! Subprogram not used             recvbuf(rdispls(mytid)+i) = sendbuf(sdispls(mytid)+i)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine altalltoallv

   
   subroutine spmd_utils_readnl(nlfile)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Read spmd utils namelist to set swap communication protocol options as
!   well as the flow control gather options
! 
! Method: 
! spmd_utils_readnl:
!
! Author of original version:  J. Truesdale
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
     use namelist_utils,  only: find_group_name
     use units,           only: getunit, freeunit
     use mpishorthand
     
     implicit none
!---------------------------Input arguments--------------------------
!
     character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

!---------------------------Local variables--------------------------
!
     integer :: unitn, ierr
     character(len=*), parameter :: subname = 'spmd_utils_readnl'
     
     namelist /spmd_utils_nl/ swap_comm_protocol,swap_comm_maxreq,fc_gather_flow_cntl

!-----------------------------------------------------------------------------

     if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'spmd_utils_nl', status=ierr)
        if (ierr == 0) then
           read(unitn, spmd_utils_nl, iostat=ierr)
           if (ierr /= 0) then
              call endrun(subname // ':: ERROR reading namelist')
           end if
           write(iulog,*) 'Read in spmd_utils_nl namelist from: ', trim(nlfile)
        end if
        close(unitn)
        call freeunit(unitn)
        
           
        if ((swap_comm_protocol < min_comm_protocol) .or. &
             (swap_comm_protocol > max_comm_protocol)) then
           write(iulog,*)                                        &
                'SPMD_UTILS_READNL:  ERROR:  swap_comm_protocol=', &
                swap_comm_protocol, ' is out of range.'
           write(iulog,*)                                        &
                '  It must be between ', min_comm_protocol,' and ',&
                max_comm_protocol
           write(iulog,*)                                        &
                '  Using default value.'
           swap_comm_protocol = def_comm_protocol
        endif
        
        write(iulog,*) 'SPMD SWAP_COMM OPTIONS: '
        write(iulog,*) '  swap_comm_protocol = ', swap_comm_protocol
        write(iulog,*) '  swap_comm_maxreq   = ', swap_comm_maxreq         
        write(iulog,*) 'SPMD FLOW CONTROL GATHER OPTION: '
        write(iulog,*) '  fc_gather_flow_cntl = ', fc_gather_flow_cntl
     endif
        
     ! Broadcast namelist variables
     call mpibcast (swap_comm_protocol , 1,   mpiint ,  0, mpicom)
     call mpibcast (swap_comm_maxreq   , 1,   mpiint ,  0, mpicom)
     call mpibcast (fc_gather_flow_cntl, 1,   mpiint ,  0, mpicom)
      
   end subroutine spmd_utils_readnl
 
 end module spmd_utils

