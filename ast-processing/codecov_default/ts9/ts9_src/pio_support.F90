!>
!! @file pio_support.F90
!! @brief internal code for compiler workarounds, aborts and debug functions
!!
!! $Revision: 823 $
!! $LastChangedDate: 2013-08-29 09:16:58 -0600 (Thu, 29 Aug 2013) $
!<
!>
!! \def _NO_MPI_RSEND
!! Code added as a work around for poor rsend performance on cray systems with
!! Gemini interconnect
!<

module pio_support
  use pio_kinds
  implicit none
  private
  include 'mpif.h'    ! _EXTERNAL
  public :: piodie
  public :: CheckMPIreturn
  public :: pio_readdof
  public :: pio_writedof
  public :: pio_fc_gather_offset


  logical, public :: Debug=.FALSE.
  logical, public :: DebugIO=.FALSE.
  logical, public :: DebugAsync=.FALSE.
  integer,private,parameter :: versno = 1001

  character(len=*), parameter :: modName='pio_support'

contains

! Subprogram not used   subroutine piodie (file,line, msg, ival1, msg2, ival2, msg3, ival3, mpirank)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! Purpose:
! Subprogram not used     !
! Subprogram not used     ! Abort the model for abnormal termination
! Subprogram not used     !
! Subprogram not used     ! Author: Jim Edwards
! Subprogram not used     !
! Subprogram not used     ! Change History
! Subprogram not used     ! 20070608 R. Loy  added optional args
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! $Id$
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     character(len=*), intent(in) :: file
! Subprogram not used     integer,intent(in) :: line
! Subprogram not used     character(len=*), intent(in), optional :: msg,msg2,msg3
! Subprogram not used     integer,intent(in),optional :: ival1,ival2,ival3, mpirank
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: subName=modName//'::pio_die'
! Subprogram not used     integer :: ierr, myrank=-1
! Subprogram not used     
! Subprogram not used     if(present(mpirank)) myrank=mpirank
! Subprogram not used 
! Subprogram not used     if (present(ival3)) then
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ', &
! Subprogram not used             msg,ival1,msg2,ival2,msg3,ival3
! Subprogram not used     else if (present(msg3)) then
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ', &
! Subprogram not used             msg,ival1,msg2,ival2,msg3
! Subprogram not used     else if (present(ival2)) then
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1,msg2,ival2
! Subprogram not used     else if (present(msg2)) then
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1,msg2
! Subprogram not used     else if (present(ival1)) then
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1
! Subprogram not used     else if (present(msg)) then
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg
! Subprogram not used     else
! Subprogram not used        write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': (no message)'
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! passing an argument of 1 to mpi_abort will lead to a STOPALL output 
! Subprogram not used     ! error code of 257
! Subprogram not used     call mpi_abort (MPI_COMM_WORLD, 1, ierr)  
! Subprogram not used 
! Subprogram not used     call abort
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine piodie

!=============================================
!  CheckMPIreturn:
!
!      Check and prints an error message
!  if an error occured in a MPI subroutine.
!=============================================
  subroutine CheckMPIreturn(locmesg, errcode, file, line)

     character(len=*), intent(in) :: locmesg
     integer(i4), intent(in) :: errcode
     character(len=*),optional :: file
     integer, intent(in),optional :: line
     character(len=MPI_MAX_ERROR_STRING) :: errorstring

     integer(i4) :: errorlen

     integer(i4) :: ierr
     if (errcode .ne. MPI_SUCCESS) then
        call MPI_Error_String(errcode,errorstring,errorlen,ierr)
        write(*,*) TRIM(ADJUSTL(locmesg))//errorstring(1:errorlen)
        if(present(file).and.present(line)) then
           call piodie(file,line)
        endif
     end if
  end subroutine CheckMPIreturn

! Subprogram not used   subroutine pio_writedof (file, DOF, comm, punit)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! Purpose:
! Subprogram not used     !
! Subprogram not used     ! Write a DOF to standard format
! Subprogram not used     !
! Subprogram not used     ! Author: T Craig
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     character(len=*),intent(in) :: file
! Subprogram not used     integer(kind=pio_offset)  ,intent(in) :: dof(:)
! Subprogram not used     integer         ,intent(in) :: comm
! Subprogram not used     integer,optional,intent(in) :: punit
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: subName=modName//'::pio_writedof'
! Subprogram not used     integer ierr, myrank, npes, m, n, unit
! Subprogram not used     integer(kind=pio_offset), pointer :: wdof(:)
! Subprogram not used     integer(kind=pio_offset), pointer :: sdof1d(:)
! Subprogram not used     integer(kind=pio_offset) :: sdof, sdof_tmp(1)
! Subprogram not used     integer          :: status(MPI_STATUS_SIZE)
! Subprogram not used     integer, parameter :: masterproc = 0
! Subprogram not used 
! Subprogram not used     integer :: pio_offset_kind                     ! kind of pio_offset
! Subprogram not used     integer :: &
! Subprogram not used       rcv_request    ,&! request id
! Subprogram not used       hs = 1           ! MPI handshaking variable
! Subprogram not used 
! Subprogram not used     unit = 81
! Subprogram not used     if (present(punit)) then
! Subprogram not used        unit = punit
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call MPI_COMM_SIZE(comm,npes,ierr)
! Subprogram not used     call CheckMPIReturn(subName,ierr)
! Subprogram not used     call MPI_COMM_RANK(comm,myrank,ierr)
! Subprogram not used     call CheckMPIReturn(subName,ierr)
! Subprogram not used     sdof = size(dof)
! Subprogram not used 
! Subprogram not used     allocate(sdof1d(0:npes-1))
! Subprogram not used     sdof1d = -1
! Subprogram not used     sdof_tmp(1) = sdof
! Subprogram not used 
! Subprogram not used     if(kind(sdof_tmp) == kind(comm)) then
! Subprogram not used        pio_offset_kind = MPI_INTEGER
! Subprogram not used     else
! Subprogram not used        pio_offset_kind = MPI_INTEGER8
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     call pio_fc_gather_offset(sdof_tmp, 1, PIO_OFFSET_KIND, &
! Subprogram not used        sdof1d, 1, PIO_OFFSET_KIND,masterproc,comm)
! Subprogram not used 
! Subprogram not used     if (myrank == masterproc) then
! Subprogram not used        write(6,*) subName,': writing file ',trim(file),' unit=',unit
! Subprogram not used        open(unit,file=file)
! Subprogram not used        write(unit,*) versno,npes
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     do n = 0,npes-1
! Subprogram not used        if (myrank == masterproc) then
! Subprogram not used           allocate(wdof(sdof1d(n)))
! Subprogram not used        endif
! Subprogram not used        if (myrank == masterproc .and. n == masterproc) then
! Subprogram not used           wdof = dof
! Subprogram not used        else
! Subprogram not used           if (myrank == n .and. sdof > 0) then
! Subprogram not used              call MPI_RECV(hs,1,MPI_INTEGER,masterproc,n,comm,status,ierr)
! Subprogram not used              if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",228,' pio_writedof mpi_recv')
! Subprogram not used              call MPI_SEND(dof,int(sdof),PIO_OFFSET_KIND,masterproc,n,comm,ierr)
! Subprogram not used              if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",231,' pio_writedof mpi_send')
! Subprogram not used           endif
! Subprogram not used           if (myrank == masterproc .and. sdof1d(n) > 0) then
! Subprogram not used              call MPI_IRECV(wdof,int(sdof1d(n)),PIO_OFFSET_KIND,n,n,comm,rcv_request,ierr)
! Subprogram not used              if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",236,' pio_writedof mpi_irecv')
! Subprogram not used              call MPI_SEND(hs,1,MPI_INTEGER,n,n,comm,ierr)
! Subprogram not used              if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",238,' pio_writedof mpi_send')
! Subprogram not used              call MPI_WAIT(rcv_request,status,ierr)
! Subprogram not used              if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",240,' pio_writedof mpi_wait')
! Subprogram not used           endif
! Subprogram not used        endif
! Subprogram not used        if (myrank == masterproc) then
! Subprogram not used           write(unit,*) n,sdof1d(n)
! Subprogram not used           do m = 1,sdof1d(n)
! Subprogram not used              write(unit,*) wdof(m)
! Subprogram not used           enddo
! Subprogram not used           deallocate(wdof)
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if (myrank == masterproc) then
! Subprogram not used        close(unit)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     deallocate(sdof1d)
! Subprogram not used 
! Subprogram not used   end subroutine pio_writedof

! Subprogram not used   subroutine pio_readdof (file, DOF, comm, punit)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! Purpose:
! Subprogram not used     !
! Subprogram not used     ! Read a DOF to standard format
! Subprogram not used     !
! Subprogram not used     ! Author: T Craig
! Subprogram not used     !
! Subprogram not used     ! Change History
! Subprogram not used     ! 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! $Id$
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     character(len=*),intent(in) :: file
! Subprogram not used     integer(kind=pio_offset),pointer:: dof(:)
! Subprogram not used     integer         ,intent(in) :: comm
! Subprogram not used     integer,optional,intent(in) :: punit
! Subprogram not used 
! Subprogram not used     character(len=*), parameter :: subName=modName//'::pio_readdof'
! Subprogram not used     integer :: ierr, myrank, npes, m, n, unit, rn
! Subprogram not used     integer(kind=pio_offset) :: sdof
! Subprogram not used     integer :: rversno, rnpes
! Subprogram not used     integer(kind=pio_offset), pointer :: wdof(:)
! Subprogram not used     integer, parameter :: masterproc = 0
! Subprogram not used     integer :: status(MPI_STATUS_SIZE)
! Subprogram not used     integer :: pio_offset_kind                        ! kind of pio_offset
! Subprogram not used 
! Subprogram not used     unit = 81
! Subprogram not used     if (present(punit)) then
! Subprogram not used        unit = punit
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call MPI_COMM_SIZE(comm,npes,ierr)
! Subprogram not used     call CheckMPIReturn(subName,ierr)
! Subprogram not used     call MPI_COMM_RANK(comm,myrank,ierr)
! Subprogram not used     call CheckMPIReturn(subName,ierr)
! Subprogram not used 
! Subprogram not used     if(kind(sdof) == kind(comm)) then
! Subprogram not used        pio_offset_kind = MPI_INTEGER
! Subprogram not used     else
! Subprogram not used        pio_offset_kind = MPI_INTEGER8
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     allocate(dof(0))   ! default for pes with no dof
! Subprogram not used 
! Subprogram not used     if (myrank == masterproc) then
! Subprogram not used        write(6,*) subName,': reading file ',trim(file),' unit=',unit
! Subprogram not used        open(unit,file=file,status='old')
! Subprogram not used        read(unit,*) rversno,rnpes
! Subprogram not used        write(6,*) subName,': reading file ',trim(file),' versno=',rversno
! Subprogram not used        if (rnpes /= npes) then
! Subprogram not used           call piodie("pio_support.F90",322,'pio_readdof npes incorrect')
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        do n = 0,npes-1
! Subprogram not used           read(unit,*) rn,sdof
! Subprogram not used           if (rn /= n) then
! Subprogram not used              call piodie("pio_support.F90",328,'pio_readdof rn out of sync')
! Subprogram not used           endif
! Subprogram not used           allocate(wdof(sdof))
! Subprogram not used           do m = 1,sdof
! Subprogram not used              read(unit,*) wdof(m)
! Subprogram not used           enddo
! Subprogram not used           if (n == masterproc) then
! Subprogram not used              deallocate(dof)
! Subprogram not used              allocate(dof(sdof))
! Subprogram not used              dof = wdof
! Subprogram not used           else
! Subprogram not used              call MPI_SEND(sdof,1,PIO_OFFSET_KIND,n,n,comm,ierr)
! Subprogram not used              if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",340,' pio_readdof mpi_send1')
! Subprogram not used              if (sdof > 0) then
! Subprogram not used                 call MPI_SEND(wdof,int(sdof),PIO_OFFSET_KIND,n,npes+n,comm,ierr)
! Subprogram not used                 if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",343,' pio_readdof mpi_send2')
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used           deallocate(wdof)
! Subprogram not used        enddo
! Subprogram not used        close(unit)
! Subprogram not used     else
! Subprogram not used        call MPI_RECV(sdof,1,PIO_OFFSET_KIND,masterproc,myrank,comm,status,ierr)
! Subprogram not used        if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",351,' pio_readdof mpi_recv1')
! Subprogram not used        if (sdof > 0) then
! Subprogram not used           deallocate(dof)
! Subprogram not used           allocate(dof(sdof))
! Subprogram not used           call MPI_RECV(dof,int(sdof),PIO_OFFSET_KIND,masterproc,npes+myrank,comm,status,ierr)
! Subprogram not used           if (ierr /= MPI_SUCCESS) call piodie("pio_support.F90",356,' pio_readdof mpi_recv2')
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine pio_readdof


!
!========================================================================
!

   subroutine pio_fc_gather_offset ( sendbuf, sendcnt, sendtype, &
                                  recvbuf, recvcnt, recvtype, &
                                  root, comm, flow_cntl )

!----------------------------------------------------------------------- 
! 
!> Purpose: 
!!   Gather collective with additional flow control, so as to 
!!   be more robust when used with high process counts. 
!!
!! Method: 
!!   If flow_cntl optional parameter 
!!     < 0: use MPI_Gather
!!     >= 0: use point-to-point with handshaking messages and 
!!           preposting receive requests up to 
!!           max(min(1,flow_cntl),max_gather_block_size) 
!!           ahead if optional flow_cntl parameter is present.
!!           Otherwise, fc_gather_flow_cntl is used in its place.
!!     Default value is 64.
!! 
!! Author of original version:  P. Worley
!! Ported from CAM: P. Worley, Jan 2010
!< 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   implicit none

!---------------------------Parameters ---------------------------------
!
   integer, parameter :: max_gather_block_size = 64

!---------------------------Input arguments--------------------------
!
   integer(kind=pio_offset), intent(in)  :: sendbuf(:)       ! outgoing message buffer
   integer, intent(in)  :: sendcnt          ! size of send buffer
   integer, intent(in)  :: sendtype         ! MPI type of send buffer
   integer, intent(in)  :: recvcnt          ! size of receive buffer
   integer, intent(in)  :: recvtype         ! MPI type of receive buffer
   integer, intent(in)  :: root             ! gather destination
   integer, intent(in)  :: comm             ! MPI communicator
   integer,optional, intent(in):: flow_cntl ! flow control variable

!---------------------------Output arguments--------------------------
!
   integer(kind=pio_offset), intent(out) :: recvbuf(*)       ! incoming message buffer
!
!---------------------------Local workspace---------------------------------
!
   character(len=*), parameter :: subName=modName//'::pio_fc_gather_int'

   logical :: fc_gather                     ! use explicit flow control?
   integer :: hs                            ! handshake variable
   integer :: gather_block_size             ! number of preposted receive requests

   integer :: nprocs                        ! size of communicator
   integer :: mytask                        ! MPI task id with communicator
   integer :: mtag                          ! MPI message tag
   integer :: p, i                          ! loop indices
   integer :: displs                        ! offset into receive buffer
   integer :: count, preposts, head, tail   ! variables controlling recv-ahead logic

   integer :: rcvid(max_gather_block_size)  ! receive request ids

   integer :: ier                           ! return error status    
   integer :: status(MPI_STATUS_SIZE)       ! MPI status 

!
!-------------------------------------------------------------------------------------
!
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
 
      ! Determine task id and size of communicator
      call mpi_comm_rank (comm, mytask, ier)
      call mpi_comm_size (comm, nprocs, ier)

      ! Initialize tag and hs variable
      mtag = 2*nprocs
      hs = 1

      if (root .eq. mytask) then

! prepost gather_block_size irecvs, and start receiving data
         preposts = min(nprocs-1, gather_block_size)
         head = 0
         count = 0
         do p=0, nprocs-1
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
                  call mpi_send ( hs, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

! copy local data
         displs = mytask*recvcnt
         do i=1,sendcnt
            recvbuf(displs+i) = sendbuf(i)
         enddo

! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv  ( hs, 1, sendtype, root, mtag, comm, &
                             status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
         end if

      endif
      call CheckMPIReturn(subName,ier)

   else
 
      call mpi_gather (sendbuf, sendcnt, sendtype, &
                       recvbuf, recvcnt, recvtype, &
                       root, comm, ier)
      call CheckMPIReturn(subName,ier)

   endif

   return

   end subroutine pio_fc_gather_offset

end module pio_support
