



module bndry_mod
  use parallel_mod, only : abortmp
  implicit none
  private
  public :: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation
  public :: ghost_exchangeV
!  public :: ghost_exchangev3d
  public :: sort_neighbor_buffer_mapping

  interface bndry_exchangeV
     module procedure bndry_exchangeV_nonth
     module procedure long_bndry_exchangeV_nonth
     module procedure bndry_exchangeV_thsave 
  end interface

contains 

  subroutine bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel

    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success



    type (parallel_t)              :: par
    type (EdgeBuffer_t)            :: buffer

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    logical(kind=log_kind),parameter              :: Debug=.FALSE.

    integer        :: i


    if(omp_in_parallel()) then 
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif


    ! Setup the pointer to proper Schedule



    pSchedule => Schedule(1)

    nlyr = buffer%nlyr

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIreal_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIreal_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle




  end subroutine bndry_exchangeV_nonth

  subroutine long_bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edge_mod, only : LongEdgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel

    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success



    type (parallel_t)              :: par
    type (LongEdgeBuffer_t)            :: buffer

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    logical(kind=log_kind),parameter              :: Debug=.FALSE.

    integer        :: i


    if(omp_in_parallel()) then 
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif


    ! Setup the pointer to proper Schedule



    pSchedule => Schedule(1)

    nlyr = buffer%nlyr

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIinteger_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIinteger_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle




  end subroutine long_bndry_exchangeV_nonth
  !********************************************************************************
  !
  !********************************************************************************
  subroutine bndry_exchangeV_thsave(hybrid,buffer)
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd, np
    use perf_mod, only: t_startf, t_stopf ! _EXTERNAL

    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success



    implicit none

    type (hybrid_t)                   :: hybrid
    type (EdgeBuffer_t)               :: buffer

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    integer        :: i
    logical(kind=log_kind),parameter      :: Debug = .FALSE.


    call t_startf('bndry_exchange')

    !$OMP BARRIER

    if(hybrid%ithr == 0) then 


       ! Setup the pointer to proper Schedule



       pSchedule => Schedule(1)

       nlyr = buffer%nlyr

       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================

       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * pCycle%lengthP
          tag             = pCycle%tag
          iptr            = pCycle%ptrP
          !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * pCycle%lengthP
          tag             = pCycle%tag
          iptr            = pCycle%ptrP
          !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
          endif
       end do    ! icycle


       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length             = pCycle%lengthP
          iptr            = pCycle%ptrP
          do i=0,length-1
             buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
          enddo
       end do   ! icycle



    endif  ! if (hybrid%ithr == 0)

    !$OMP BARRIER

    call t_stopf('bndry_exchange')

  end subroutine bndry_exchangeV_thsave





! Subprogram not used   subroutine ghost_exchangeVfull(hybrid,buffer)
! Subprogram not used !
! Subprogram not used !   MT 2011:  derived from bndry_exchange, but copies an entire
! Subprogram not used !             element of ghost cell information, including corner
! Subprogram not used !             elements.  Requres cubed-sphere grid
! Subprogram not used !
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     use kinds, only : log_kind
! Subprogram not used     use edge_mod, only : Ghostbuffer3D_t
! Subprogram not used     use schedule_mod, only : schedule_t, cycle_t, schedule
! Subprogram not used     use dimensions_mod, only: nelemd
! Subprogram not used 
! Subprogram not used     use parallel_mod, only : abortmp, status, srequest, rrequest, &
! Subprogram not used          mpireal_t, mpiinteger_t, mpi_success
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (hybrid_t)                   :: hybrid
! Subprogram not used     type (GhostBuffer3D_t)               :: buffer
! Subprogram not used 
! Subprogram not used     type (Schedule_t),pointer                     :: pSchedule
! Subprogram not used     type (Cycle_t),pointer                        :: pCycle
! Subprogram not used     integer                                       :: dest,length,tag
! Subprogram not used     integer                                       :: icycle,ierr
! Subprogram not used     integer                                       :: iptr,source,nlyr
! Subprogram not used     integer                                       :: nSendCycles,nRecvCycles
! Subprogram not used     integer                                       :: errorcode,errorlen
! Subprogram not used     character*(80) errorstring
! Subprogram not used 
! Subprogram not used     integer        :: i,i1,i2
! Subprogram not used     logical(kind=log_kind),parameter      :: Debug = .FALSE.
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !$OMP BARRIER
! Subprogram not used 
! Subprogram not used     if(hybrid%ithr == 0) then 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! Setup the pointer to proper Schedule
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        pSchedule => Schedule(1)
! Subprogram not used 
! Subprogram not used        nlyr = buffer%nlyr
! Subprogram not used 
! Subprogram not used        nSendCycles = pSchedule%nSendCycles
! Subprogram not used        nRecvCycles = pSchedule%nRecvCycles
! Subprogram not used 
! Subprogram not used        !==================================================
! Subprogram not used        !  Fire off the sends
! Subprogram not used        !==================================================
! Subprogram not used        do icycle=1,nSendCycles
! Subprogram not used           pCycle      => pSchedule%SendCycle(icycle)
! Subprogram not used           dest            = pCycle%dest - 1
! Subprogram not used           length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
! Subprogram not used           tag             = pCycle%tag
! Subprogram not used           iptr            = pCycle%ptrP_ghost
! Subprogram not used           !print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
! Subprogram not used           call MPI_Isend(buffer%buf(1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
! Subprogram not used           if(ierr .ne. MPI_SUCCESS) then
! Subprogram not used              errorcode=ierr
! Subprogram not used              call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
! Subprogram not used              print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
! Subprogram not used           endif
! Subprogram not used        end do    ! icycle
! Subprogram not used 
! Subprogram not used        !==================================================
! Subprogram not used        !  Post the Receives 
! Subprogram not used        !==================================================
! Subprogram not used        do icycle=1,nRecvCycles
! Subprogram not used           pCycle         => pSchedule%RecvCycle(icycle)
! Subprogram not used           source          = pCycle%source - 1
! Subprogram not used           length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
! Subprogram not used           tag             = pCycle%tag
! Subprogram not used           iptr            = pCycle%ptrP_ghost
! Subprogram not used           !print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
! Subprogram not used           call MPI_Irecv(buffer%receive(1,1,1,iptr),length,MPIreal_t, &
! Subprogram not used                source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
! Subprogram not used           if(ierr .ne. MPI_SUCCESS) then
! Subprogram not used              errorcode=ierr
! Subprogram not used              call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
! Subprogram not used              print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
! Subprogram not used           endif
! Subprogram not used        end do    ! icycle
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        !==================================================
! Subprogram not used        !  Wait for all the receives to complete
! Subprogram not used        !==================================================
! Subprogram not used 
! Subprogram not used        call MPI_Waitall(nSendCycles,Srequest,status,ierr)
! Subprogram not used        call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
! Subprogram not used 
! Subprogram not used        do icycle=1,nRecvCycles
! Subprogram not used           pCycle         => pSchedule%RecvCycle(icycle)
! Subprogram not used           length             = pCycle%lengthP_ghost
! Subprogram not used           iptr            = pCycle%ptrP_ghost
! Subprogram not used           do i=0,length-1
! Subprogram not used              buffer%buf(:,:,1:nlyr,iptr+i) = buffer%receive(:,:,1:nlyr,iptr+i)
! Subprogram not used           enddo
! Subprogram not used        end do   ! icycle
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     endif  ! if (hybrid%ithr == 0)
! Subprogram not used 
! Subprogram not used     !$OMP BARRIER
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine ghost_exchangeVfull

  ! ===========================================
  !  GHOST_EXCHANGEV:
  !  Author: Christoph Erath
  !  derived from bndry_exchange, but copies an entire
  !             element of ghost cell information, including corner
  !             elements.  Requres cubed-sphere grid
  ! =========================================
! Subprogram not used  subroutine ghost_exchangeV(hybrid,buffer,nhc,npoints,ntrac)
! Subprogram not used !
! Subprogram not used !   2011:  derived from bndry_exchange, but copies an entire
! Subprogram not used !             element of ghost cell information, including corner
! Subprogram not used !             elements.  Requres cubed-sphere grid
! Subprogram not used !
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     use kinds, only : log_kind
! Subprogram not used     use edge_mod, only : Ghostbuffertr_t
! Subprogram not used     use schedule_mod, only : schedule_t, cycle_t, schedule
! Subprogram not used     use dimensions_mod, only: nelemd
! Subprogram not used 
! Subprogram not used     use parallel_mod, only : abortmp, status, srequest, rrequest, &
! Subprogram not used          mpireal_t, mpiinteger_t, mpi_success
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (hybrid_t)                   :: hybrid
! Subprogram not used     type (GhostBuffertr_t)               :: buffer
! Subprogram not used     integer :: nhc,npoints,ntrac
! Subprogram not used 
! Subprogram not used     type (Schedule_t),pointer                     :: pSchedule
! Subprogram not used     type (Cycle_t),pointer                        :: pCycle
! Subprogram not used     integer                                       :: dest,length,tag
! Subprogram not used     integer                                       :: icycle,ierr
! Subprogram not used     integer                                       :: iptr,source,nlyr
! Subprogram not used     integer                                       :: nSendCycles,nRecvCycles
! Subprogram not used     integer                                       :: errorcode,errorlen
! Subprogram not used     character*(80) errorstring
! Subprogram not used 
! Subprogram not used     integer        :: i,i1,i2
! Subprogram not used     logical(kind=log_kind),parameter      :: Debug = .FALSE.
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !$OMP BARRIER
! Subprogram not used 
! Subprogram not used     if(hybrid%ithr == 0) then 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! Setup the pointer to proper Schedule
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        pSchedule => Schedule(1)
! Subprogram not used 
! Subprogram not used        nlyr = buffer%nlyr
! Subprogram not used               
! Subprogram not used        nSendCycles = pSchedule%nSendCycles
! Subprogram not used        nRecvCycles = pSchedule%nRecvCycles
! Subprogram not used 
! Subprogram not used        !==================================================
! Subprogram not used        !  Fire off the sends
! Subprogram not used        !==================================================
! Subprogram not used        do icycle=1,nSendCycles
! Subprogram not used           pCycle      => pSchedule%SendCycle(icycle)
! Subprogram not used           dest            = pCycle%dest - 1
! Subprogram not used           length      = nlyr * ntrac * pCycle%lengthP_ghost*nhc*npoints
! Subprogram not used           tag             = pCycle%tag
! Subprogram not used           iptr            = pCycle%ptrP_ghost
! Subprogram not used           !print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
! Subprogram not used           call MPI_Isend(buffer%buf(1,1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
! Subprogram not used           if(ierr .ne. MPI_SUCCESS) then
! Subprogram not used              errorcode=ierr
! Subprogram not used              call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
! Subprogram not used              print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
! Subprogram not used           endif
! Subprogram not used        end do    ! icycle
! Subprogram not used 
! Subprogram not used        !==================================================
! Subprogram not used        !  Post the Receives 
! Subprogram not used        !==================================================
! Subprogram not used        do icycle=1,nRecvCycles
! Subprogram not used           pCycle         => pSchedule%RecvCycle(icycle)
! Subprogram not used           source          = pCycle%source - 1
! Subprogram not used           length      = nlyr * ntrac * pCycle%lengthP_ghost*nhc*npoints
! Subprogram not used           tag             = pCycle%tag
! Subprogram not used           iptr            = pCycle%ptrP_ghost
! Subprogram not used           !print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
! Subprogram not used           call MPI_Irecv(buffer%receive(1,1,1,1,iptr),length,MPIreal_t, &
! Subprogram not used                source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
! Subprogram not used           if(ierr .ne. MPI_SUCCESS) then
! Subprogram not used              errorcode=ierr
! Subprogram not used              call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
! Subprogram not used              print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
! Subprogram not used           endif
! Subprogram not used        end do    ! icycle
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        !==================================================
! Subprogram not used        !  Wait for all the receives to complete
! Subprogram not used        !==================================================
! Subprogram not used 
! Subprogram not used        call MPI_Waitall(nSendCycles,Srequest,status,ierr)
! Subprogram not used        call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
! Subprogram not used 
! Subprogram not used        do icycle=1,nRecvCycles
! Subprogram not used           pCycle         => pSchedule%RecvCycle(icycle)
! Subprogram not used           length             = pCycle%lengthP_ghost
! Subprogram not used           iptr            = pCycle%ptrP_ghost
! Subprogram not used           do i=0,length-1
! Subprogram not used              do i2=1,nhc
! Subprogram not used                 do i1=1,npoints
! Subprogram not used                    buffer%buf(i1,i2,1:nlyr,1:ntrac,iptr+i) = buffer%receive(i1,i2,1:nlyr,1:ntrac,iptr+i)
! Subprogram not used                 enddo
! Subprogram not used              enddo
! Subprogram not used           enddo
! Subprogram not used        end do   ! icycle
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     endif  ! if (hybrid%ithr == 0)
! Subprogram not used 
! Subprogram not used     !$OMP BARRIER
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine ghost_exchangeV





! Subprogram not used   subroutine compute_ghost_corner_orientation(hybrid,elem,nets,nete)
! Subprogram not used !
! Subprogram not used !  this routine can NOT be called in a threaded region because then each thread
! Subprogram not used !  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
! Subprogram not used !
! Subprogram not used   use kinds, only : real_kind
! Subprogram not used   use dimensions_mod, only: nelemd, np
! Subprogram not used   use parallel_mod, only : syncmp
! Subprogram not used   use hybrid_mod, only : hybrid_t
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use edge_mod, only : ghostbuffer3D_t, ghostvpackfull, ghostvunpackfull, &
! Subprogram not used        initghostbuffer3D,freeghostbuffer3D
! Subprogram not used   use control_mod, only : north,south,east,west,neast, nwest, seast, swest
! Subprogram not used 
! Subprogram not used   implicit none
! Subprogram not used 
! Subprogram not used   type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used   type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used   integer :: nets,nete
! Subprogram not used   type (ghostBuffer3D_t)   :: ghostbuf_cv
! Subprogram not used 
! Subprogram not used   real (kind=real_kind) :: cin(2,2,1,nets:nete)  !CE: fvm tracer
! Subprogram not used   real (kind=real_kind) :: cout(-1:4,-1:4,1,nets:nete)  !CE: fvm tracer
! Subprogram not used   integer :: i,j,ie,kptr,np1,np2,nc,nc1,nc2,k,nlev
! Subprogram not used   logical :: fail,fail1,fail2
! Subprogram not used   real (kind=real_kind) :: tol=.1
! Subprogram not used   call syncmp(hybrid%par)
! Subprogram not used !   if (hybrid%par%masterproc) print *,'computing ghost cell corner orientations'
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   ! first test on the Gauss Grid with same number of ghost cells:
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   nc=2   ! test using GLL interior points
! Subprogram not used   nc1=-1
! Subprogram not used   nc2=4
! Subprogram not used 
! Subprogram not used   nlev=1
! Subprogram not used   call initghostbuffer3D(ghostbuf_cv,nlev,nc)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   do ie=nets,nete
! Subprogram not used      cin(1,1,1,ie)=  elem(ie)%gdofp(1,1)
! Subprogram not used      cin(nc,nc,1,ie)=  elem(ie)%gdofp(np,np)
! Subprogram not used      cin(1,nc,1,ie)=   elem(ie)%gdofp(1,np)
! Subprogram not used      cin(nc,1,1,ie)=  elem(ie)%gdofp(np,1)
! Subprogram not used   enddo
! Subprogram not used   cout=0
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  run ghost exchange on c array to get corner orientation
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   do ie=nets,nete
! Subprogram not used      kptr=0
! Subprogram not used      call ghostVpackfull(ghostbuf_cv, cin(:,:,:,ie),1,nc,nc,nlev,kptr,elem(ie)%desc)
! Subprogram not used   end do
! Subprogram not used   call ghost_exchangeVfull(hybrid,ghostbuf_cv)
! Subprogram not used   do ie=nets,nete
! Subprogram not used      kptr=0
! Subprogram not used      call ghostVunpackfull(ghostbuf_cv, cout(:,:,:,ie), nc1,nc2,nc,nlev, kptr, elem(ie)%desc)
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used !       nc +--------+ 
! Subprogram not used !        ^ | nw  ne |    
! Subprogram not used !     j  | |        |
! Subprogram not used !        1 | sw  se |
! Subprogram not used !          +--------+
! Subprogram not used !           1 --> nc
! Subprogram not used !              i
! Subprogram not used 
! Subprogram not used ! check SW corner
! Subprogram not used   do ie=nets,nete
! Subprogram not used      fail1=.false.
! Subprogram not used      fail2=.false.
! Subprogram not used      if ( elem(ie)%desc%putmapP_ghost(swest) /= -1) then
! Subprogram not used         if (abs(cout(nc1,1,1,ie)-cout(nc1,0,1,ie)) .gt. tol )  fail1=.true.
! Subprogram not used         if (abs(cout(1,nc1,1,ie)-cout(0,nc1,1,ie)).gt.tol) fail2=.true.
! Subprogram not used      endif
! Subprogram not used      if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange SW orientation failure')
! Subprogram not used      if (fail1) then
! Subprogram not used         elem(ie)%desc%reverse(swest)=.true.
! Subprogram not used         !print *,'reversion sw orientation ie',ie
! Subprogram not used         !print *,elem(ie)%desc%reverse(nwest),elem(ie)%desc%reverse(north),elem(ie)%desc%reverse(neast)
! Subprogram not used         !print *,elem(ie)%desc%reverse(west),' ',elem(ie)%desc%reverse(east)
! Subprogram not used         !print *,elem(ie)%desc%reverse(swest),elem(ie)%desc%reverse(south),elem(ie)%desc%reverse(seast)
! Subprogram not used      endif
! Subprogram not used   enddo
! Subprogram not used ! check SE corner
! Subprogram not used   do ie=nets,nete
! Subprogram not used      fail1=.false.
! Subprogram not used      fail2=.false.
! Subprogram not used      if ( elem(ie)%desc%putmapP_ghost(seast) /= -1) then
! Subprogram not used         if (abs(cout(nc2,1,1,ie)-cout(nc2,0,1,ie)) .gt. tol )  fail1=.true.
! Subprogram not used         if (abs(cout(nc+1,nc1,1,ie)-cout(nc,nc1,1,ie)).gt.tol) fail2=.true.
! Subprogram not used      endif
! Subprogram not used      if (fail1 .neqv. fail2 ) call abortmp('ghost exchange SE orientation failure')
! Subprogram not used      if (fail1) then
! Subprogram not used         elem(ie)%desc%reverse(seast)=.true.
! Subprogram not used      endif
! Subprogram not used   enddo
! Subprogram not used ! check NW corner
! Subprogram not used   do ie=nets,nete
! Subprogram not used      fail1=.false.
! Subprogram not used      fail2=.false.
! Subprogram not used      if ( elem(ie)%desc%putmapP_ghost(nwest) /= -1) then
! Subprogram not used         if (abs(cout(nc1,nc+1,1,ie)-cout(nc1,nc,1,ie)) .gt. tol )  fail1=.true.
! Subprogram not used         if (abs(cout(1,nc2,1,ie)-cout(0,nc2,1,ie)).gt.tol) fail2=.true.
! Subprogram not used      endif
! Subprogram not used      if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange NW orientation failure')
! Subprogram not used      if (fail1) then
! Subprogram not used         elem(ie)%desc%reverse(nwest)=.true.
! Subprogram not used      endif
! Subprogram not used   enddo
! Subprogram not used ! check NE corner
! Subprogram not used   do ie=nets,nete
! Subprogram not used      fail1=.false.
! Subprogram not used      fail2=.false.
! Subprogram not used      if ( elem(ie)%desc%putmapP_ghost(neast) /= -1) then
! Subprogram not used         if (abs(cout(nc2,nc+1,1,ie)-cout(nc2,nc,1,ie)) .gt. tol )  fail1=.true.
! Subprogram not used         if (abs(cout(nc+1,nc2,1,ie)-cout(nc,nc2,1,ie)).gt.tol) fail2=.true.
! Subprogram not used      endif
! Subprogram not used      if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange NE orientation failure')
! Subprogram not used      if (fail1) then
! Subprogram not used         elem(ie)%desc%reverse(neast)=.true.
! Subprogram not used      endif
! Subprogram not used   enddo
! Subprogram not used   call freeghostbuffer3D(ghostbuf_cv)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  end ghost exchange corner orientation
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   end subroutine





! Subprogram not used   subroutine sort_neighbor_buffer_mapping(hybrid,elem,nets,nete)
! Subprogram not used !
! Subprogram not used !  gather global ID's of all neighbor elements.  Then create sorted (in global ID numbering)
! Subprogram not used !  mapping between edge buffer for each neighbor and a local map.  
! Subprogram not used !
! Subprogram not used !  this routine can NOT be called in a threaded region because then each thread
! Subprogram not used !  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
! Subprogram not used !
! Subprogram not used !  also return num_neigh(ie) = number of neighbors (including onself) for element ie
! Subprogram not used !  
! Subprogram not used !
! Subprogram not used   use kinds, only : real_kind
! Subprogram not used   use dimensions_mod, only: nelemd, np, max_neigh_edges
! Subprogram not used   use parallel_mod, only : syncmp
! Subprogram not used   use hybrid_mod, only : hybrid_t
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use edge_mod, only : ghostbuffer3D_t, ghostvpack_unoriented, ghostvunpack_unoriented, &
! Subprogram not used        initghostbuffer3D,freeghostbuffer3D
! Subprogram not used   use control_mod, only : north,south,east,west,neast, nwest, seast, swest
! Subprogram not used   use coordinate_systems_mod, only: cartesian3D_t
! Subprogram not used   implicit none
! Subprogram not used 
! Subprogram not used   type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used   type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used   integer :: nets,nete
! Subprogram not used   type (ghostBuffer3D_t)   :: ghostbuf_cv
! Subprogram not used 
! Subprogram not used   real (kind=real_kind) :: cin(2,2,4,nets:nete)                    ! 1x1 element input data
! Subprogram not used   real (kind=real_kind) :: cout(2,2,4,max_neigh_edges,nets:nete)   ! 1x1 element output data
! Subprogram not used   integer :: i,j,ie,kptr,np1,np2,nc,k,nlev,actual_neigh_edges,l,l2,sum1,sum2
! Subprogram not used   logical :: fail,fail1,fail2
! Subprogram not used   real (kind=real_kind) :: tol=.1
! Subprogram not used 
! Subprogram not used   type (cartesian3D_t)     :: neigh_corners(4,max_neigh_edges,nelemd)  
! Subprogram not used 
! Subprogram not used   call syncmp(hybrid%par)
! Subprogram not used   if (hybrid%par%masterproc) print *,'checking ghost cell neighbor buffer sorting...'
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   ! first test on the Gauss Grid with same number of ghost cells:
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   nc=2
! Subprogram not used   nlev=4
! Subprogram not used   call initghostbuffer3D(ghostbuf_cv,nlev,nc)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   do ie=nets,nete
! Subprogram not used      cin(:,:,nlev,ie)=  elem(ie)%GlobalID
! Subprogram not used      k=0
! Subprogram not used      do i=1,nc
! Subprogram not used      do j=1,nc
! Subprogram not used         k=k+1
! Subprogram not used         cin(i,j,1,ie) = elem(ie)%corners3D(k)%x
! Subprogram not used         cin(i,j,2,ie) = elem(ie)%corners3D(k)%y
! Subprogram not used         cin(i,j,3,ie) = elem(ie)%corners3D(k)%z
! Subprogram not used      enddo
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used   cout=-1
! Subprogram not used 
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used !  run ghost exchange to get global ID of all neighbors
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   do ie=nets,nete
! Subprogram not used      kptr=0
! Subprogram not used      call ghostVpack_unoriented(ghostbuf_cv, cin(:,:,:,ie),nc,nlev,kptr,elem(ie)%desc)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   ! check for array out of bouds overwriting 
! Subprogram not used   if (int(maxval(  cout(:,:,:,:,:))) /= -1 ) then
! Subprogram not used      call abortmp('ghost excchange unoriented failure ob1')
! Subprogram not used   endif
! Subprogram not used   call ghost_exchangeVfull(hybrid,ghostbuf_cv)
! Subprogram not used   if (int(maxval(  cout(:,:,:,:,:))) /= -1 ) then
! Subprogram not used      call abortmp('ghost excchange unoriented failure ob2')
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   do ie=nets,nete
! Subprogram not used      kptr=0
! Subprogram not used      call ghostVunpack_unoriented(ghostbuf_cv, cout(:,:,:,:,ie),nc,nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used      ! check that we get the count of real neighbors correct
! Subprogram not used      actual_neigh_edges=0
! Subprogram not used 
! Subprogram not used      do l=1,max_neigh_edges
! Subprogram not used         if (int(cout(1,1,nlev,l,ie)) /= -1 ) then
! Subprogram not used            actual_neigh_edges = actual_neigh_edges + 1
! Subprogram not used         endif
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used      if (elem(ie)%desc%actual_neigh_edges /= actual_neigh_edges) then
! Subprogram not used         print *,'desc  actual_neigh_edges: ',elem(ie)%desc%actual_neigh_edges
! Subprogram not used         print *,'check actual_neigh_edges: ',actual_neigh_edges
! Subprogram not used         call abortmp( 'ghost exchange unoriented failure 1')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      ! check that all non-neighbors stayed -1
! Subprogram not used      if ( actual_neigh_edges < max_neigh_edges ) then
! Subprogram not used         do l=actual_neigh_edges+1,max_neigh_edges
! Subprogram not used         if (int(cout(1,1,nlev,l,ie)) /= -1 ) then
! Subprogram not used            call abortmp( 'ghost exchange unoriented failure 2')
! Subprogram not used         endif
! Subprogram not used         enddo
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      ! i am too lazy to check if all id's are identical since they are in
! Subprogram not used      ! different order.  check if there sum is identical
! Subprogram not used      sum1 = sum(int(cout(1,1,nlev,1:actual_neigh_edges,ie)))
! Subprogram not used      sum2=0
! Subprogram not used      do l=1,max_neigh_edges
! Subprogram not used         if (elem(ie)%desc%globalID(l)>0) sum2 = sum2 + elem(ie)%desc%globalID(l)
! Subprogram not used      enddo
! Subprogram not used      if (sum1 /= sum2 ) then
! Subprogram not used         print *,int(cin(1,1,nlev,ie)),elem(ie)%desc%actual_neigh_edges,actual_neigh_edges
! Subprogram not used         write(*,'(a,99i5)') 'ghost=',int(cout(1,1,nlev,1:actual_neigh_edges,ie))
! Subprogram not used         write(*,'(a,99i5)') 'desc =',elem(ie)%desc%globalID(:)
! Subprogram not used 
! Subprogram not used         print *,'cout sum of all neighbor global ids:',sum1 
! Subprogram not used         print *,'desc sum of all neighbor global ids:',sum2  
! Subprogram not used         call abortmp( 'ghost exchange unoriented failure 3')        
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      ! unpack corner data into array
! Subprogram not used      do l=1,elem(ie)%desc%actual_neigh_edges
! Subprogram not used      k=0
! Subprogram not used      do i=1,nc
! Subprogram not used      do j=1,nc
! Subprogram not used         k=k+1
! Subprogram not used         neigh_corners(k,l,ie)%x=cout(i,j,1,l,ie) 
! Subprogram not used         neigh_corners(k,l,ie)%y=cout(i,j,2,l,ie) 
! Subprogram not used         neigh_corners(k,l,ie)%z=cout(i,j,3,l,ie) 
! Subprogram not used      enddo
! Subprogram not used      enddo
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   call freeghostbuffer3D(ghostbuf_cv)
! Subprogram not used   if (hybrid%par%masterproc) print *,'passed.'
! Subprogram not used   end subroutine




end module bndry_mod
