!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_gather_scatter

 module ice_gather_scatter

! !DESCRIPTION:
!  This module contains routines for gathering data to a single
!  processor from a distributed array, and scattering data from a
!  single processor to a distributed array.
!
!  NOTE: The arrays gathered and scattered are assumed to have
!        horizontal dimensions (nx_block, ny_block).
!
! !REVISION HISTORY:
!  SVN:$Id: ice_gather_scatter.F90 131 2008-05-30 16:53:40Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke replaced old routines with new POP
!              infrastructure, added specialized routine scatter_global_stress

! !USES:

   use ice_kinds_mod
   use ice_communicate
   use ice_constants
   use ice_blocks
   use ice_distribution
   use ice_domain
   use ice_domain_size
   use ice_exit

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: gather_global,      &
             scatter_global,     &
             gatherArray,        & 
             scatter_global_stress

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface gather_global
     module procedure gather_global_dbl,  &
                      gather_global_real, &
                      gather_global_int
   end interface 

   interface scatter_global
     module procedure scatter_global_dbl,  &
                      scatter_global_real, &
                      scatter_global_int
   end interface 

   interface gatherArray
     module procedure gatherArray_dbl
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains


! Subprogram not used  subroutine gatherArray_dbl(array_g,array,length,root)
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used    real(dbl_kind) :: array_g(:)  ! The concatonated array
! Subprogram not used    real(dbl_kind) :: array(:)    ! the local piece of the array
! Subprogram not used    integer(int_kind) :: length   ! number of elements in the array
! Subprogram not used    integer(int_kind) :: root     ! root to which to collect the array
! Subprogram not used 
! Subprogram not used    integer(int_kind) :: ierr
! Subprogram not used 
! Subprogram not used    call MPI_Gather(array,length,MPI_REAL8,array_g, length,MPI_REAL8,root, &
! Subprogram not used                    MPI_COMM_ICE,ierr)
! Subprogram not used 
! Subprogram not used  end subroutine gatherArray_dbl


!***********************************************************************
!BOP
! !IROUTINE: gather_global
! !INTERFACE:

! Subprogram not used  subroutine gather_global_dbl(ARRAY_G, ARRAY, dst_task, src_dist)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This subroutine gathers a distributed array to a global-sized
! Subprogram not used !  array on the processor dst_task.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is the specific inteface for double precision arrays 
! Subprogram not used !  corresponding to the generic interface gather_global.  It is shown
! Subprogram not used !  to provide information on the generic interface (the generic
! Subprogram not used !  interface is identical, but chooses a specific inteface based
! Subprogram not used !  on the data type of the input argument).
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used      dst_task   ! task to which array should be gathered
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used      src_dist   ! distribution of blocks in the source array
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used      ARRAY      ! array containing horizontal slab of distributed field
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:), intent(inout) :: &
! Subprogram not used      ARRAY_G    ! array containing global horizontal field on dst_task
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used      i,j,n          ,&! dummy loop counters
! Subprogram not used      nsends         ,&! number of actual sends
! Subprogram not used      src_block      ,&! block locator for send
! Subprogram not used      src_task       ,&! source of message
! Subprogram not used      ierr             ! MPI error flag
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
! Subprogram not used      status
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used      snd_request
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      snd_status
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), allocatable :: &
! Subprogram not used      msg_buffer
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used      this_block  ! block info for current block
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: ib, ig, itask, nprocs, maxBlocks
! Subprogram not used    integer (int_kind) :: iig,ijg,it
! Subprogram not used    integer (int_kind) :: msgLen, msgTag
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if this task is the dst_task, copy local blocks into the global 
! Subprogram not used !  array and post receives for non-local blocks.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    nprocs = get_num_procs()
! Subprogram not used 
! Subprogram not used    if (my_task == dst_task) then
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used 
! Subprogram not used        !*** copy local blocks
! Subprogram not used 
! Subprogram not used        if (src_dist%blockLocation(n) == my_task+1) then
! Subprogram not used 
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          do j=this_block%jlo,this_block%jhi
! Subprogram not used          do i=this_block%ilo,this_block%ihi
! Subprogram not used            ARRAY_G(this_block%i_glob(i), &
! Subprogram not used                    this_block%j_glob(j)) = &
! Subprogram not used                   ARRAY(i,j,src_dist%blockLocalID(n))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used        !*** fill land blocks with zeroes
! Subprogram not used 
! Subprogram not used        else if (src_dist%blockLocation(n) == 0) then
! Subprogram not used 
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          do j=this_block%jlo,this_block%jhi
! Subprogram not used          do i=this_block%ilo,this_block%ihi
! Subprogram not used            ARRAY_G(this_block%i_glob(i), &
! Subprogram not used                    this_block%j_glob(j)) = c0
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      !*** receive blocks to fill up the rest
! Subprogram not used      allocate (msg_buffer(nx_block,ny_block,max_blocks))
! Subprogram not used      do itask = 0,nprocs-1
! Subprogram not used         if(itask /= dst_task) then 
! Subprogram not used            maxBlocks = src_dist%BlockCnt(itask+1)
! Subprogram not used            msgLen = nx_block*ny_block*maxBlocks
! Subprogram not used            msgTag = 3*mpitag_gs+itask
! Subprogram not used            if(maxBlocks>0) then 
! Subprogram not used 
! Subprogram not used              call MPI_RECV(msg_buffer(:,:,1:maxBlocks),msgLen, &
! Subprogram not used                mpiR8, itask, msgTag, MPI_COMM_ICE, status, ierr)
! Subprogram not used 
! Subprogram not used              do ib=1,maxBlocks
! Subprogram not used                ig = src_dist%blockIndex(itask+1,ib)
! Subprogram not used       	       this_block = get_block(ig,ig)
! Subprogram not used                do j=this_block%jlo,this_block%jhi
! Subprogram not used                   do i=this_block%ilo,this_block%ihi
! Subprogram not used                      iig = this_block%i_glob(i) 
! Subprogram not used                      ijg = this_block%j_glob(j) 
! Subprogram not used                      ARRAY_G(iig,ijg) = msg_buffer(i,j,ib)
! Subprogram not used                   end do
! Subprogram not used                end do
! Subprogram not used               enddo
! Subprogram not used             endif
! Subprogram not used         endif
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used      deallocate(msg_buffer)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise send data to dst_task
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      if(nblocks>0) then 
! Subprogram not used         msgLen = nx_block*ny_block*nblocks
! Subprogram not used         msgTag = 3*mpitag_gs+my_task
! Subprogram not used 	call MPI_SEND(ARRAY(:,:,1:nblocks),msgLen,mpiR8,dst_task, &
! Subprogram not used            msgTag,MPI_COMM_ICE,ierr)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine gather_global_dbl

!***********************************************************************

! Subprogram not used  subroutine gather_global_real(ARRAY_G, ARRAY, dst_task, src_dist)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  This subroutine gathers a distributed array to a global-sized
! Subprogram not used !  array on the processor dst_task.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  input variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used      dst_task       ! task to which array should be gathered
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used      src_dist       ! distribution of blocks in the source array
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used      ARRAY          ! array containing distributed field
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  output variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:), intent(inout) :: &
! Subprogram not used      ARRAY_G        ! array containing global field on dst_task
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used      i,j,n          ,&! dummy loop counters
! Subprogram not used      nsends         ,&! number of actual sends
! Subprogram not used      src_block      ,&! block locator for send
! Subprogram not used      ierr             ! MPI error flag
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
! Subprogram not used      status
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used      snd_request
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      snd_status
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), allocatable :: &
! Subprogram not used      msg_buffer
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used      this_block  ! block info for current block
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: ib, ig, itask, nprocs, maxBlocks
! Subprogram not used    integer (int_kind) :: iig,ijg
! Subprogram not used    integer (int_kind) :: msgLen, msgTag
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if this task is the dst_task, copy local blocks into the global 
! Subprogram not used !  array and post receives for non-local blocks.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    nprocs = get_num_procs()
! Subprogram not used 
! Subprogram not used    if (my_task == dst_task) then
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used 
! Subprogram not used        !*** copy local blocks
! Subprogram not used 
! Subprogram not used        if (src_dist%blockLocation(n) == my_task+1) then
! Subprogram not used 
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          do j=this_block%jlo,this_block%jhi
! Subprogram not used          do i=this_block%ilo,this_block%ihi
! Subprogram not used            ARRAY_G(this_block%i_glob(i), &
! Subprogram not used                    this_block%j_glob(j)) = &
! Subprogram not used                   ARRAY(i,j,src_dist%blockLocalID(n))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used        !*** fill land blocks with special values
! Subprogram not used 
! Subprogram not used        else if (src_dist%blockLocation(n) == 0) then
! Subprogram not used 
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          do j=this_block%jlo,this_block%jhi
! Subprogram not used          do i=this_block%ilo,this_block%ihi
! Subprogram not used            ARRAY_G(this_block%i_glob(i), &
! Subprogram not used                    this_block%j_glob(j)) = spval
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      !*** receive blocks to fill up the rest
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      allocate (msg_buffer(nx_block,ny_block,max_blocks))
! Subprogram not used      do itask = 0,nprocs-1
! Subprogram not used         if( itask /= dst_task) then
! Subprogram not used            maxBlocks = src_dist%BlockCnt(itask+1)
! Subprogram not used            msgLen = nx_block*ny_block*maxBlocks
! Subprogram not used            msgTag = 3*mpitag_gs+itask
! Subprogram not used            if( maxBlocks>0) then 
! Subprogram not used 
! Subprogram not used               call MPI_RECV(msg_buffer(:,:,1:maxBlocks), msgLen, &
! Subprogram not used                    mpiR4, itask, msgTag ,MPI_COMM_ICE, status, ierr)
! Subprogram not used 
! Subprogram not used               do ib=1,maxBlocks
! Subprogram not used                  ig = src_dist%blockIndex(itask+1,ib)
! Subprogram not used                  this_block = get_block(ig,ig)
! Subprogram not used                  do j=this_block%jlo,this_block%jhi
! Subprogram not used                    do i=this_block%ilo,this_block%ihi
! Subprogram not used                      iig = this_block%i_glob(i) 
! Subprogram not used                      ijg = this_block%j_glob(j) 
! Subprogram not used                      ARRAY_G(iig,ijg) = msg_buffer(i,j,ib)
! Subprogram not used                    end do
! Subprogram not used                  end do
! Subprogram not used               enddo
! Subprogram not used            endif
! Subprogram not used         endif
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used      deallocate(msg_buffer)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise send data to dst_task
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used      if(nblocks>0) then
! Subprogram not used         msgLen = nx_block*ny_block*nblocks
! Subprogram not used         msgTag = 3*mpitag_gs+my_task 
! Subprogram not used 	call MPI_SEND(ARRAY(:,:,1:nblocks),msgLen,mpiR4,dst_task, &
! Subprogram not used            msgTag,MPI_COMM_ICE,ierr)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine gather_global_real

!***********************************************************************

! Subprogram not used  subroutine gather_global_int(ARRAY_G, ARRAY, dst_task, src_dist)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  This subroutine gathers a distributed array to a global-sized
! Subprogram not used !  array on the processor dst_task.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  input variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used      dst_task       ! task to which array should be gathered
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used      src_dist       ! distribution of blocks in the source array
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used      ARRAY          ! array containing distributed field
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  output variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), intent(inout) :: &
! Subprogram not used      ARRAY_G        ! array containing global field on dst_task
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used      i,j,n          ,&! dummy loop counters
! Subprogram not used      nsends         ,&! number of actual sends
! Subprogram not used      src_block      ,&! block locator for send
! Subprogram not used      ierr             ! MPI error flag
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
! Subprogram not used      status
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used      snd_request
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      snd_status
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), allocatable :: &
! Subprogram not used      msg_buffer
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used      this_block  ! block info for current block
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: ib, ig, len, itask, nprocs, maxBlocks
! Subprogram not used    integer (int_kind) :: iig,ijg
! Subprogram not used    integer (int_kind) :: msgLen, msgTag
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if this task is the dst_task, copy local blocks into the global 
! Subprogram not used !  array and post receives for non-local blocks.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    nprocs = get_num_procs()
! Subprogram not used 
! Subprogram not used    if (my_task == dst_task) then
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used 
! Subprogram not used        !*** copy local blocks
! Subprogram not used 
! Subprogram not used        if (src_dist%blockLocation(n) == my_task+1) then
! Subprogram not used 
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          do j=this_block%jlo,this_block%jhi
! Subprogram not used          do i=this_block%ilo,this_block%ihi
! Subprogram not used            ARRAY_G(this_block%i_glob(i), &
! Subprogram not used                    this_block%j_glob(j)) = &
! Subprogram not used                   ARRAY(i,j,src_dist%blockLocalID(n))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used        !*** fill land blocks with zeroes
! Subprogram not used 
! Subprogram not used        else if (src_dist%blockLocation(n) == 0) then
! Subprogram not used 
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          do j=this_block%jlo,this_block%jhi
! Subprogram not used          do i=this_block%ilo,this_block%ihi
! Subprogram not used            ARRAY_G(this_block%i_glob(i), &
! Subprogram not used                    this_block%j_glob(j)) = 0
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      !*** receive blocks to fill up the rest
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      allocate (msg_buffer(nx_block,ny_block,max_blocks))
! Subprogram not used      do itask = 0,nprocs-1
! Subprogram not used         if( itask /= dst_task) then
! Subprogram not used            maxBlocks = src_dist%BlockCnt(itask+1)
! Subprogram not used            msgLen = nx_block*ny_block*max_blocks
! Subprogram not used            msgTag = 3*mpitag_gs+itask
! Subprogram not used            if(maxBLocks>0) then 
! Subprogram not used              call MPI_RECV(msg_buffer(:,:,1:maxBlocks),msgLen, &
! Subprogram not used                    mpi_integer, itask, msgTag, MPI_COMM_ICE, status, ierr)
! Subprogram not used              do ib=1,maxBlocks
! Subprogram not used                 ig = src_dist%blockIndex(itask+1,ib)
! Subprogram not used                 this_block = get_block(ig,ig)
! Subprogram not used                 do j=this_block%jlo,this_block%jhi
! Subprogram not used                    do i=this_block%ilo,this_block%ihi
! Subprogram not used                       iig = this_block%i_glob(i) 
! Subprogram not used                       ijg = this_block%j_glob(j) 
! Subprogram not used                       ARRAY_G(iig,ijg) = msg_buffer(i,j,ib)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used               enddo
! Subprogram not used            endif
! Subprogram not used         endif
! Subprogram not used      enddo
! Subprogram not used      deallocate(msg_buffer)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise send data to dst_task
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      if(nblocks>0) then
! Subprogram not used         msgLen = nx_block*ny_block*nblocks
! Subprogram not used         msgTag = 3*mpitag_gs+my_task
! Subprogram not used         call MPI_SEND(ARRAY(:,:,1:nblocks),msgLen, &
! Subprogram not used              mpi_integer, dst_task, msgTag, MPI_COMM_ICE, ierr)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine gather_global_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: scatter_global
! !INTERFACE:

 subroutine scatter_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

! !DESCRIPTION:
!  This subroutine scatters a global-sized array to a distributed array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = -1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = -1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         msg_buffer(i,j-yoffset2) = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         ARRAY(i,j-yoffset2,dst_block) &
                           = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR8, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Ensure unused ghost cell values are 0
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
      do n=1,nblocks_tot
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         if (dst_block > 0) then

         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo

         endif
      enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl

!***********************************************************************

! Subprogram not used  subroutine scatter_global_real(ARRAY, ARRAY_G, src_task, dst_dist, &
! Subprogram not used                                field_loc, field_type)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  This subroutine scatters a global-sized array to a distributed array.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  input variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used      src_task       ! task from which array should be scattered
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used      dst_dist       ! distribution of resulting blocks
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:), intent(in) :: &
! Subprogram not used      ARRAY_G        ! array containing global field on src_task
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_type,               &! id for type of field (scalar, vector, angle)
! Subprogram not used       field_loc                  ! id for location on horizontal grid
! Subprogram not used                                  !  (center, NEcorner, Nface, Eface)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  output variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(inout) :: &
! Subprogram not used      ARRAY          ! array containing distributed field
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used      i,j,n,bid,          &! dummy loop indices
! Subprogram not used      nrecvs,             &! actual number of messages received
! Subprogram not used      isrc, jsrc,         &! source addresses
! Subprogram not used      dst_block,          &! location of block in dst array
! Subprogram not used      xoffset, yoffset,   &! offsets for tripole boundary conditions
! Subprogram not used      yoffset2,           &!
! Subprogram not used      isign,              &! sign factor for tripole boundary conditions
! Subprogram not used      ierr                 ! MPI error flag
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used      this_block  ! block info for current block
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used      rcv_request     ! request array for receives
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      rcv_status      ! status array for receives
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:), allocatable :: &
! Subprogram not used      msg_buffer      ! buffer for sending blocks
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  initialize return array to zero and set up tripole quantities
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ARRAY = 0._real_kind
! Subprogram not used 
! Subprogram not used    this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
! Subprogram not used    if (this_block%tripoleTFlag) then
! Subprogram not used      select case (field_loc)
! Subprogram not used      case (field_loc_center)   ! cell center location
! Subprogram not used         xoffset = 2
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_NEcorner) ! cell corner (velocity) location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_Eface)    ! cell face location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_Nface)    ! cell face location
! Subprogram not used         xoffset = 2
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_noupdate) ! ghost cells never used - use cell center
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      end select
! Subprogram not used    else
! Subprogram not used      select case (field_loc)
! Subprogram not used      case (field_loc_center)   ! cell center location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_NEcorner) ! cell corner (velocity) location
! Subprogram not used         xoffset = 0
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_Eface)    ! cell face location
! Subprogram not used         xoffset = 0
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_Nface)    ! cell face location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_noupdate) ! ghost cells never used - use cell center
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      end select
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    select case (field_type)
! Subprogram not used    case (field_type_scalar)
! Subprogram not used       isign =  1
! Subprogram not used    case (field_type_vector)
! Subprogram not used       isign = -1
! Subprogram not used    case (field_type_angle)
! Subprogram not used       isign = -1
! Subprogram not used    case (field_type_noupdate) ! ghost cells never used - use cell center
! Subprogram not used       isign =  1
! Subprogram not used    case default
! Subprogram not used       call abort_ice('Unknown field type in scatter')
! Subprogram not used    end select
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if this task is the src_task, copy blocks of global array into 
! Subprogram not used !  message buffer and send to other processors. also copy local blocks
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task == src_task) then
! Subprogram not used 
! Subprogram not used      !*** send non-local blocks away
! Subprogram not used 
! Subprogram not used      allocate (msg_buffer(nx_block,ny_block))
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) > 0 .and. &
! Subprogram not used            dst_dist%blockLocation(n)-1 /= my_task) then
! Subprogram not used 
! Subprogram not used          msg_buffer = 0._real_kind
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          !*** if this is an interior block, then there is no
! Subprogram not used          !*** padding or update checking required
! Subprogram not used 
! Subprogram not used          if (this_block%iblock > 1         .and. &
! Subprogram not used              this_block%iblock < nblocks_x .and. &
! Subprogram not used              this_block%jblock > 1         .and. &
! Subprogram not used              this_block%jblock < nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used             do i=1,nx_block
! Subprogram not used                msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                          this_block%j_glob(j))
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is an edge block but not a northern edge
! Subprogram not used          !*** we only need to check for closed boundaries and
! Subprogram not used          !*** padding (global index = 0)
! Subprogram not used 
! Subprogram not used          else if (this_block%jblock /= nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) /= 0) then
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                   this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is a northern edge block, we need to check
! Subprogram not used          !*** for and properly deal with tripole boundaries
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) > 0) then ! normal boundary
! Subprogram not used 
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                   this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                else if (this_block%j_glob(j) < 0) then  ! tripole
! Subprogram not used 
! Subprogram not used                   ! for yoffset=0 or 1, yoffset2=0,0
! Subprogram not used                   ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
! Subprogram not used                   do yoffset2=0,max(yoffset,0)-yoffset
! Subprogram not used                     jsrc = ny_global + yoffset + yoffset2 + &
! Subprogram not used                          (this_block%j_glob(j) + ny_global)
! Subprogram not used                     do i=1,nx_block
! Subprogram not used                       if (this_block%i_glob(i) /= 0) then
! Subprogram not used                          isrc = nx_global + xoffset - this_block%i_glob(i)
! Subprogram not used                          if (isrc < 1) isrc = isrc + nx_global
! Subprogram not used                          if (isrc > nx_global) isrc = isrc - nx_global
! Subprogram not used                          msg_buffer(i,j-yoffset2) = isign * ARRAY_G(isrc,jsrc)
! Subprogram not used                       endif
! Subprogram not used                     end do
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          call MPI_SEND(msg_buffer, nx_block*ny_block, &
! Subprogram not used                        mpiR4, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
! Subprogram not used                        MPI_COMM_ICE, ierr)
! Subprogram not used 
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      deallocate(msg_buffer)
! Subprogram not used 
! Subprogram not used      !*** copy any local blocks
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) == my_task+1) then
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          !*** if this is an interior block, then there is no
! Subprogram not used          !*** padding or update checking required
! Subprogram not used 
! Subprogram not used          if (this_block%iblock > 1         .and. &
! Subprogram not used              this_block%iblock < nblocks_x .and. &
! Subprogram not used              this_block%jblock > 1         .and. &
! Subprogram not used              this_block%jblock < nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used             do i=1,nx_block
! Subprogram not used                ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                               this_block%j_glob(j))
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is an edge block but not a northern edge
! Subprogram not used          !*** we only need to check for closed boundaries and
! Subprogram not used          !*** padding (global index = 0)
! Subprogram not used 
! Subprogram not used          else if (this_block%jblock /= nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) /= 0) then
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                        this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is a northern edge block, we need to check
! Subprogram not used          !*** for and properly deal with tripole boundaries
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) > 0) then ! normal boundary
! Subprogram not used 
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                        this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                else if (this_block%j_glob(j) < 0) then  ! tripole
! Subprogram not used 
! Subprogram not used                   ! for yoffset=0 or 1, yoffset2=0,0
! Subprogram not used                   ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
! Subprogram not used                   do yoffset2=0,max(yoffset,0)-yoffset
! Subprogram not used                     jsrc = ny_global + yoffset + yoffset2 + &
! Subprogram not used                          (this_block%j_glob(j) + ny_global)
! Subprogram not used                     do i=1,nx_block
! Subprogram not used                       if (this_block%i_glob(i) /= 0) then
! Subprogram not used                          isrc = nx_global + xoffset - this_block%i_glob(i)
! Subprogram not used                          if (isrc < 1) isrc = isrc + nx_global
! Subprogram not used                          if (isrc > nx_global) isrc = isrc - nx_global
! Subprogram not used                          ARRAY(i,j-yoffset2,dst_block) &
! Subprogram not used                            = isign * ARRAY_G(isrc,jsrc)
! Subprogram not used                       endif
! Subprogram not used                     end do
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise receive data from src_task
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used      allocate (rcv_request(nblocks_tot), &
! Subprogram not used                rcv_status(MPI_STATUS_SIZE, nblocks_tot))
! Subprogram not used 
! Subprogram not used      rcv_request = 0
! Subprogram not used      rcv_status  = 0
! Subprogram not used 
! Subprogram not used      nrecvs = 0
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) == my_task+1) then
! Subprogram not used          nrecvs = nrecvs + 1
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
! Subprogram not used                        mpiR4, src_task, 3*mpitag_gs+n, &
! Subprogram not used                        MPI_COMM_ICE, rcv_request(nrecvs), ierr)
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      if (nrecvs > 0) &
! Subprogram not used        call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)
! Subprogram not used 
! Subprogram not used      deallocate(rcv_request, rcv_status)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    !-----------------------------------------------------------------
! Subprogram not used    ! Ensure unused ghost cell values are 0
! Subprogram not used    !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (field_loc == field_loc_noupdate) then
! Subprogram not used       do n=1,nblocks_tot
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          if (dst_block > 0) then
! Subprogram not used 
! Subprogram not used          ! north edge
! Subprogram not used          do j = this_block%jhi+1,ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             ARRAY (i,j,dst_block) = 0._real_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! east edge
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = this_block%ihi+1,nx_block
! Subprogram not used             ARRAY (i,j,dst_block) = 0._real_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! south edge
! Subprogram not used          do j = 1, this_block%jlo-1
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             ARRAY (i,j,dst_block) = 0._real_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! west edge
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, this_block%ilo-1
! Subprogram not used             ARRAY (i,j,dst_block) = 0._real_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine scatter_global_real

!***********************************************************************

! Subprogram not used  subroutine scatter_global_int(ARRAY, ARRAY_G, src_task, dst_dist, &
! Subprogram not used                                field_loc, field_type)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  This subroutine scatters a global-sized array to a distributed array.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  input variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used      src_task       ! task from which array should be scattered
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_type,               &! id for type of field (scalar, vector, angle)
! Subprogram not used       field_loc                  ! id for location on horizontal grid
! Subprogram not used                                  !  (center, NEcorner, Nface, Eface)
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used      dst_dist       ! distribution of resulting blocks
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), intent(in) :: &
! Subprogram not used      ARRAY_G        ! array containing global field on src_task
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  output variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(inout) :: &
! Subprogram not used      ARRAY          ! array containing distributed field
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used      i,j,n,bid,          &! dummy loop indices
! Subprogram not used      nrecvs,             &! actual number of messages received
! Subprogram not used      isrc, jsrc,         &! source addresses
! Subprogram not used      dst_block,          &! location of block in dst array
! Subprogram not used      xoffset, yoffset,   &! offsets for tripole boundary conditions
! Subprogram not used      yoffset2,           &!
! Subprogram not used      isign,              &! sign factor for tripole boundary conditions
! Subprogram not used      ierr                 ! MPI error flag
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used      this_block  ! block info for current block
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used      rcv_request     ! request array for receives
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      rcv_status      ! status array for receives
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      msg_buffer      ! buffer for sending blocks
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  initialize return array to zero and set up tripole quantities
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ARRAY = 0
! Subprogram not used 
! Subprogram not used    this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
! Subprogram not used    if (this_block%tripoleTFlag) then
! Subprogram not used      select case (field_loc)
! Subprogram not used      case (field_loc_center)   ! cell center location
! Subprogram not used         xoffset = 2
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_NEcorner) ! cell corner (velocity) location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_Eface)    ! cell face location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_Nface)    ! cell face location
! Subprogram not used         xoffset = 2
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_noupdate) ! ghost cells never used - use cell center
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      end select
! Subprogram not used    else
! Subprogram not used      select case (field_loc)
! Subprogram not used      case (field_loc_center)   ! cell center location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_NEcorner) ! cell corner (velocity) location
! Subprogram not used         xoffset = 0
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_Eface)    ! cell face location
! Subprogram not used         xoffset = 0
! Subprogram not used         yoffset = 1
! Subprogram not used      case (field_loc_Nface)    ! cell face location
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 0
! Subprogram not used      case (field_loc_noupdate) ! ghost cells never used - use cell center
! Subprogram not used         xoffset = 1
! Subprogram not used         yoffset = 1
! Subprogram not used      end select
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    select case (field_type)
! Subprogram not used    case (field_type_scalar)
! Subprogram not used       isign =  1
! Subprogram not used    case (field_type_vector)
! Subprogram not used       isign = -1
! Subprogram not used    case (field_type_angle)
! Subprogram not used       isign = -1
! Subprogram not used    case (field_type_noupdate) ! ghost cells never used - use cell center
! Subprogram not used       isign =  1
! Subprogram not used    case default
! Subprogram not used       call abort_ice('Unknown field type in scatter')
! Subprogram not used    end select
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if this task is the src_task, copy blocks of global array into 
! Subprogram not used !  message buffer and send to other processors. also copy local blocks
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task == src_task) then
! Subprogram not used 
! Subprogram not used      !*** send non-local blocks away
! Subprogram not used 
! Subprogram not used      allocate (msg_buffer(nx_block,ny_block))
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) > 0 .and. &
! Subprogram not used            dst_dist%blockLocation(n)-1 /= my_task) then
! Subprogram not used 
! Subprogram not used          msg_buffer = 0
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          !*** if this is an interior block, then there is no
! Subprogram not used          !*** padding or update checking required
! Subprogram not used 
! Subprogram not used          if (this_block%iblock > 1         .and. &
! Subprogram not used              this_block%iblock < nblocks_x .and. &
! Subprogram not used              this_block%jblock > 1         .and. &
! Subprogram not used              this_block%jblock < nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used             do i=1,nx_block
! Subprogram not used                msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                          this_block%j_glob(j))
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is an edge block but not a northern edge
! Subprogram not used          !*** we only need to check for closed boundaries and
! Subprogram not used          !*** padding (global index = 0)
! Subprogram not used 
! Subprogram not used          else if (this_block%jblock /= nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) /= 0) then
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                   this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is a northern edge block, we need to check
! Subprogram not used          !*** for and properly deal with tripole boundaries
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) > 0) then ! normal boundary
! Subprogram not used 
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                   this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                else if (this_block%j_glob(j) < 0) then  ! tripole
! Subprogram not used 
! Subprogram not used                   ! for yoffset=0 or 1, yoffset2=0,0
! Subprogram not used                   ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
! Subprogram not used                   do yoffset2=0,max(yoffset,0)-yoffset
! Subprogram not used                     jsrc = ny_global + yoffset + yoffset2 + &
! Subprogram not used                          (this_block%j_glob(j) + ny_global)
! Subprogram not used                     do i=1,nx_block
! Subprogram not used                       if (this_block%i_glob(i) /= 0) then
! Subprogram not used                          isrc = nx_global + xoffset - this_block%i_glob(i)
! Subprogram not used                          if (isrc < 1) isrc = isrc + nx_global
! Subprogram not used                          if (isrc > nx_global) isrc = isrc - nx_global
! Subprogram not used                          msg_buffer(i,j-yoffset2) = isign * ARRAY_G(isrc,jsrc)
! Subprogram not used                       endif
! Subprogram not used                     end do
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          call MPI_SEND(msg_buffer, nx_block*ny_block, &
! Subprogram not used                        mpi_integer, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
! Subprogram not used                        MPI_COMM_ICE, ierr)
! Subprogram not used 
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      deallocate(msg_buffer)
! Subprogram not used 
! Subprogram not used      !*** copy any local blocks
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) == my_task+1) then
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          !*** if this is an interior block, then there is no
! Subprogram not used          !*** padding or update checking required
! Subprogram not used 
! Subprogram not used          if (this_block%iblock > 1         .and. &
! Subprogram not used              this_block%iblock < nblocks_x .and. &
! Subprogram not used              this_block%jblock > 1         .and. &
! Subprogram not used              this_block%jblock < nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used             do i=1,nx_block
! Subprogram not used                ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                               this_block%j_glob(j))
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is an edge block but not a northern edge
! Subprogram not used          !*** we only need to check for closed boundaries and
! Subprogram not used          !*** padding (global index = 0)
! Subprogram not used 
! Subprogram not used          else if (this_block%jblock /= nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) /= 0) then
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                        this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is a northern edge block, we need to check
! Subprogram not used          !*** for and properly deal with tripole boundaries
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) > 0) then ! normal boundary
! Subprogram not used 
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
! Subprogram not used                                                        this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                else if (this_block%j_glob(j) < 0) then  ! tripole
! Subprogram not used 
! Subprogram not used                   ! for yoffset=0 or 1, yoffset2=0,0
! Subprogram not used                   ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
! Subprogram not used                   do yoffset2=0,max(yoffset,0)-yoffset
! Subprogram not used                     jsrc = ny_global + yoffset + yoffset2 + &
! Subprogram not used                          (this_block%j_glob(j) + ny_global)
! Subprogram not used                     do i=1,nx_block
! Subprogram not used                       if (this_block%i_glob(i) /= 0) then
! Subprogram not used                          isrc = nx_global + xoffset - this_block%i_glob(i)
! Subprogram not used                          if (isrc < 1) isrc = isrc + nx_global
! Subprogram not used                          if (isrc > nx_global) isrc = isrc - nx_global
! Subprogram not used                          ARRAY(i,j-yoffset2,dst_block) &
! Subprogram not used                            = isign * ARRAY_G(isrc,jsrc)
! Subprogram not used                       endif
! Subprogram not used                     end do
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise receive data from src_task
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used      allocate (rcv_request(nblocks_tot), &
! Subprogram not used                rcv_status(MPI_STATUS_SIZE, nblocks_tot))
! Subprogram not used 
! Subprogram not used      rcv_request = 0
! Subprogram not used      rcv_status  = 0
! Subprogram not used 
! Subprogram not used      nrecvs = 0
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) == my_task+1) then
! Subprogram not used          nrecvs = nrecvs + 1
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
! Subprogram not used                        mpi_integer, src_task, 3*mpitag_gs+n, &
! Subprogram not used                        MPI_COMM_ICE, rcv_request(nrecvs), ierr)
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      if (nrecvs > 0) &
! Subprogram not used        call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)
! Subprogram not used 
! Subprogram not used      deallocate(rcv_request, rcv_status)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    !-----------------------------------------------------------------
! Subprogram not used    ! Ensure unused ghost cell values are 0
! Subprogram not used    !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (field_loc == field_loc_noupdate) then
! Subprogram not used       do n=1,nblocks_tot
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          if (dst_block > 0) then
! Subprogram not used 
! Subprogram not used          ! north edge
! Subprogram not used          do j = this_block%jhi+1,ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             ARRAY (i,j,dst_block) = 0
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! east edge
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = this_block%ihi+1,nx_block
! Subprogram not used             ARRAY (i,j,dst_block) = 0
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! south edge
! Subprogram not used          do j = 1, this_block%jlo-1
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             ARRAY (i,j,dst_block) = 0
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! west edge
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, this_block%ilo-1
! Subprogram not used             ARRAY (i,j,dst_block) = 0
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine scatter_global_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: scatter_global_stress
! !INTERFACE:

! Subprogram not used  subroutine scatter_global_stress(ARRAY, ARRAY_G1, ARRAY_G2, &
! Subprogram not used                                   src_task, dst_dist)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This subroutine scatters global stresses to a distributed array.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  Ghost cells in the stress tensor must be handled separately on tripole
! Subprogram not used !  grids, because matching the corner values requires 2 different arrays.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used    include 'mpif.h'
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used      src_task       ! task from which array should be scattered
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used      dst_dist       ! distribution of resulting blocks
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:), intent(in) :: &
! Subprogram not used      ARRAY_G1,     &! array containing global field on src_task
! Subprogram not used      ARRAY_G2       ! array containing global field on src_task
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), intent(inout) :: &
! Subprogram not used      ARRAY          ! array containing distributed field
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used      i,j,n,bid,          &! dummy loop indices
! Subprogram not used      nrecvs,             &! actual number of messages received
! Subprogram not used      isrc, jsrc,         &! source addresses
! Subprogram not used      dst_block,          &! location of block in dst array
! Subprogram not used      xoffset, yoffset,   &! offsets for tripole boundary conditions
! Subprogram not used      yoffset2,           &!
! Subprogram not used      isign,              &! sign factor for tripole boundary conditions
! Subprogram not used      ierr                 ! MPI error flag
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used      this_block  ! block info for current block
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used      rcv_request     ! request array for receives
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:), allocatable :: &
! Subprogram not used      rcv_status      ! status array for receives
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:), allocatable :: &
! Subprogram not used      msg_buffer      ! buffer for sending blocks
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  initialize return array to zero and set up tripole quantities
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ARRAY = c0
! Subprogram not used 
! Subprogram not used    this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
! Subprogram not used    if (this_block%tripoleTFlag) then
! Subprogram not used      xoffset = 2  ! treat stresses as cell-centered scalars (they are not 
! Subprogram not used      yoffset = 0  ! shared with neighboring grid cells)
! Subprogram not used    else
! Subprogram not used      xoffset = 1  ! treat stresses as cell-centered scalars (they are not 
! Subprogram not used      yoffset = 1  ! shared with neighboring grid cells)
! Subprogram not used    endif
! Subprogram not used    isign   = 1
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if this task is the src_task, copy blocks of global array into 
! Subprogram not used !  message buffer and send to other processors. also copy local blocks
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task == src_task) then
! Subprogram not used 
! Subprogram not used      !*** send non-local blocks away
! Subprogram not used 
! Subprogram not used      allocate (msg_buffer(nx_block,ny_block))
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) > 0 .and. &
! Subprogram not used            dst_dist%blockLocation(n)-1 /= my_task) then
! Subprogram not used 
! Subprogram not used          msg_buffer = c0
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          !*** if this is an interior block, then there is no
! Subprogram not used          !*** padding or update checking required
! Subprogram not used 
! Subprogram not used          if (this_block%iblock > 1         .and. &
! Subprogram not used              this_block%iblock < nblocks_x .and. &
! Subprogram not used              this_block%jblock > 1         .and. &
! Subprogram not used              this_block%jblock < nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used             do i=1,nx_block
! Subprogram not used                msg_buffer(i,j) = ARRAY_G1(this_block%i_glob(i),&
! Subprogram not used                                           this_block%j_glob(j))
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is an edge block but not a northern edge
! Subprogram not used          !*** we only need to check for closed boundaries and
! Subprogram not used          !*** padding (global index = 0)
! Subprogram not used 
! Subprogram not used          else if (this_block%jblock /= nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) /= 0) then
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         msg_buffer(i,j) = ARRAY_G1(this_block%i_glob(i),&
! Subprogram not used                                                   this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is a northern edge block, we need to check
! Subprogram not used          !*** for and properly deal with tripole boundaries
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) > 0) then ! normal boundary
! Subprogram not used 
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         msg_buffer(i,j) = ARRAY_G1(this_block%i_glob(i),&
! Subprogram not used                                                   this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                else if (this_block%j_glob(j) < 0) then  ! tripole
! Subprogram not used 
! Subprogram not used                   jsrc = ny_global + yoffset + &
! Subprogram not used                          (this_block%j_glob(j) + ny_global)
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         isrc = nx_global + xoffset - this_block%i_glob(i)
! Subprogram not used                         if (isrc < 1) isrc = isrc + nx_global
! Subprogram not used                         if (isrc > nx_global) isrc = isrc - nx_global
! Subprogram not used                         msg_buffer(i,j) = isign * ARRAY_G2(isrc,jsrc)
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          call MPI_SEND(msg_buffer, nx_block*ny_block, &
! Subprogram not used                        mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
! Subprogram not used                        MPI_COMM_ICE, ierr)
! Subprogram not used 
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      deallocate(msg_buffer)
! Subprogram not used 
! Subprogram not used      !*** copy any local blocks
! Subprogram not used 
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) == my_task+1) then
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          this_block = get_block(n,n)
! Subprogram not used 
! Subprogram not used          !*** if this is an interior block, then there is no
! Subprogram not used          !*** padding or update checking required
! Subprogram not used 
! Subprogram not used          if (this_block%iblock > 1         .and. &
! Subprogram not used              this_block%iblock < nblocks_x .and. &
! Subprogram not used              this_block%jblock > 1         .and. &
! Subprogram not used              this_block%jblock < nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used             do i=1,nx_block
! Subprogram not used                ARRAY(i,j,dst_block) = ARRAY_G1(this_block%i_glob(i),&
! Subprogram not used                                               this_block%j_glob(j))
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is an edge block but not a northern edge
! Subprogram not used          !*** we only need to check for closed boundaries and
! Subprogram not used          !*** padding (global index = 0)
! Subprogram not used 
! Subprogram not used          else if (this_block%jblock /= nblocks_y) then
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) /= 0) then
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         ARRAY(i,j,dst_block) = ARRAY_G1(this_block%i_glob(i),&
! Subprogram not used                                                        this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          !*** if this is a northern edge block, we need to check
! Subprogram not used          !*** for and properly deal with tripole boundaries
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used 
! Subprogram not used             do j=1,ny_block
! Subprogram not used                if (this_block%j_glob(j) > 0) then ! normal boundary
! Subprogram not used 
! Subprogram not used                   do i=1,nx_block
! Subprogram not used                      if (this_block%i_glob(i) /= 0) then
! Subprogram not used                         ARRAY(i,j,dst_block) = ARRAY_G1(this_block%i_glob(i),&
! Subprogram not used                                                        this_block%j_glob(j))
! Subprogram not used                      endif
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                else if (this_block%j_glob(j) < 0) then  ! tripole
! Subprogram not used 
! Subprogram not used                   ! for yoffset=0 or 1, yoffset2=0,0
! Subprogram not used                   ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
! Subprogram not used                   do yoffset2=0,max(yoffset,0)-yoffset
! Subprogram not used                     jsrc = ny_global + yoffset + yoffset2 + &
! Subprogram not used                          (this_block%j_glob(j) + ny_global)
! Subprogram not used                     do i=1,nx_block
! Subprogram not used                       if (this_block%i_glob(i) /= 0) then
! Subprogram not used                          isrc = nx_global + xoffset - this_block%i_glob(i)
! Subprogram not used                          if (isrc < 1) isrc = isrc + nx_global
! Subprogram not used                          if (isrc > nx_global) isrc = isrc - nx_global
! Subprogram not used                          ARRAY(i,j-yoffset2,dst_block) &
! Subprogram not used                            = isign * ARRAY_G2(isrc,jsrc)
! Subprogram not used                       endif
! Subprogram not used                     end do
! Subprogram not used                   end do
! Subprogram not used 
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise receive data from src_task
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used      allocate (rcv_request(nblocks_tot), &
! Subprogram not used                rcv_status(MPI_STATUS_SIZE, nblocks_tot))
! Subprogram not used 
! Subprogram not used      rcv_request = 0
! Subprogram not used      rcv_status  = 0
! Subprogram not used 
! Subprogram not used      nrecvs = 0
! Subprogram not used      do n=1,nblocks_tot
! Subprogram not used        if (dst_dist%blockLocation(n) == my_task+1) then
! Subprogram not used          nrecvs = nrecvs + 1
! Subprogram not used          dst_block = dst_dist%blockLocalID(n)
! Subprogram not used          call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
! Subprogram not used                        mpiR8, src_task, 3*mpitag_gs+n, &
! Subprogram not used                        MPI_COMM_ICE, rcv_request(nrecvs), ierr)
! Subprogram not used        endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      if (nrecvs > 0) &
! Subprogram not used        call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)
! Subprogram not used 
! Subprogram not used      deallocate(rcv_request, rcv_status)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine scatter_global_stress

!***********************************************************************

 end module ice_gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
