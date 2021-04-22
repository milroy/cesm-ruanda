!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_distribution

 module ice_distribution

!
! !DESCRIPTION:
!  This module provides data types and routines for distributing
!  blocks across processors.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_distribution.F90 118 2008-04-08 20:57:17Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke updated to new POP infrastructure
!
! !USES:

   use ice_kinds_mod
   use ice_domain_size
   use ice_communicate
   use ice_blocks
   use ice_exit
   use ice_fileunits, only: nu_diag, nu_timing, ice_stdout, flush_fileunit
   use ice_spacecurve
   use ice_broadcast
!   use ice_probability, only:  EstimateCost, BuildProbabilityStats, WoriteProbabilityStats
   use ice_probability_tools

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: distrb  ! distribution data type
      integer (int_kind) :: &
         nprocs            ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (int_kind), dimension(:), pointer :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block

      integer (int_kind), dimension(:), pointer ::  blockCnt
      integer (int_kind), dimension(:,:), pointer :: blockIndex

   end type

   integer (int_kind), public, parameter ::  &    ! types of blocks:
         lndType     = 0,                    &    ! 	Land
         icefreeType = 1,                    &    !     ice free (ocean only)
         iceType     = 2                          !     sea ice 

!    integer, parameter :: numCoeff = 5



! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution, &
             ice_distributionGet,         &
             ice_distributionGetBlockLoc, &
             ice_distributionGetBlockID, &
             create_local_block_ids

! !PUBLIC DATA MEMBERS:

   character (char_len), public :: &
       processor_shape       ! 'square-pop' (approx) POP default config
                             ! 'square-ice' like square-pop but better for ice
                             ! 'slenderX1' (NPX x 1)
                             ! 'slenderX2' (NPX x 2)

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 function create_distribution(dist_type, nprocs, minBlock, maxBlock, work_per_block, prob_per_block, blockType, bStats, &
      FixMaxBlock, maxDil)

! !DESCRIPTION:
!  This routine determines the distribution of blocks across processors
!  by call the appropriate subroutine based on distribution type
!  requested.  Currently three distributions are supported:
!  2-d Cartesian distribution (cartesian), a load-balanced
!  distribution using a rake algorithm based on an input amount of work 
!  per block, and a space-filling-curve algorithm.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      dist_type             ! method for distributing blocks
                            !  either cartesian or rake

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), intent(in) :: minBlock ! minimum number of blocks to use
   integer (int_kind), intent(in) :: maxBlock ! maximum number of blocks to use

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

   real (dbl_kind), dimension(:), intent(in) :: &
      prob_per_block        ! probability that block contains sea-ice

   integer (int_kind), dimension(:), intent(in) :: &
      blockType             ! type of block

   real (dbl_kind), dimension(:,:), intent(in)  :: bStats  ! block statistics 

   logical, intent(in) :: FixMaxBlock
   real (real_kind), intent(in) :: maxDil

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distribution   ! resulting structure describing
                            !  distribution of blocks

!EOP

!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, nb, nc        ! dummy counters

!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('roundrobin')

      create_distribution = create_distrb_roundrobin(nprocs, work_per_block)

   case('blkrobin')

      create_distribution = create_distrb_blkrobin(nprocs, work_per_block)

   case('blkcart')

      create_distribution = create_distrb_blkcart(nprocs, work_per_block)

   case('cartesian')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('rake')

      create_distribution = create_distrb_rake(nprocs, work_per_block)

   case('spacecurve')

!DBG      write(nu_diag,*) 'before call to create_distrb_spacecurve'
      create_distribution = create_distrb_spacecurve(nprocs, minBlock, maxBlock, &
                            work_per_block,prob_per_block,blockType,bStats, FixMaxBlock, maxDil )
!DBG    stop 'create_distribution: after call to create_distrb_spacecurve'

   case default

      call abort_ice('ice distribution: unknown distribution type')

   end select

   if(my_task == master_task) then 
      write(nu_diag,*) ' '
      nc = 0
      do n = 1,nprocs
         nb = create_distribution%blockcnt(n)
         if (nb > 0) then
            nc = nc + 1
!            write(nu_diag,*) ' Blocks on proc : ',n,nb,':', &
!               create_distribution%blockindex(n,1:nb)
         else
!            write(nu_diag,*) ' Blocks on proc : ',n,nb
         endif
      enddo
      write(nu_diag,*) ' Active processors: ',nc,MAXVAL(create_distribution%blockLocation)
      write(nu_diag,*) ' '
   endif
!DBG write(nu_diag,*) 'end of create_distribution'
!-----------------------------------------------------------------------
!EOC

 end function create_distribution

!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids
! !INTERFACE:

 subroutine create_local_block_ids(block_ids, distribution)

! !DESCRIPTION:
!  This routine determines which blocks in an input distribution are
!  located on the local processor and creates an array of block ids
!  for all local blocks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which local
                             !  blocks required

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), pointer :: &
      block_ids              ! array of block ids for every block
                             ! that resides on the local processor
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bid, bcount        ! dummy counters

   logical (log_kind) :: dbug

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%blockLocation)
      if (distribution%blockLocation(n) == my_task+1) bcount = bcount + 1
   end do


   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

!   dbug = .true.
   dbug = .false.
   if (bcount > 0) then
      do n=1,size(distribution%blockLocation)
         if (distribution%blockLocation(n) == my_task+1) then
            block_ids(distribution%blockLocalID(n)) = n

            if (dbug) then
            write(nu_diag,*) 'block id, proc, local_block: ', &
                             block_ids(distribution%blockLocalID(n)), &
                             distribution%blockLocation(n), &
                             distribution%blockLocalID(n)
            endif
         endif
      end do
   endif

!EOC

 end subroutine create_local_block_ids

!***********************************************************************
!BOP
! !IROUTINE: proc_decomposition
! !INTERFACE:

! Subprogram not used  subroutine proc_decomposition(nprocs, nprocs_x, nprocs_y)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This subroutine attempts to find an optimal (nearly square)
! Subprogram not used !  2d processor decomposition for a given number of processors.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used    use ice_domain_size
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       nprocs                       ! total number or processors
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(out) :: &
! Subprogram not used       nprocs_x, nprocs_y           ! number of procs in each dimension
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       iguess, jguess               ! guesses for nproc_x,y
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       square                       ! square root of nprocs
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  start with an initial guess
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    square = sqrt(real(nprocs,kind=real_kind))
! Subprogram not used    nprocs_x = 0
! Subprogram not used    nprocs_y = 0
! Subprogram not used 
! Subprogram not used    if (processor_shape == 'square-pop') then ! make as square as possible
! Subprogram not used       iguess = nint(square)
! Subprogram not used       jguess = nprocs/iguess
! Subprogram not used    elseif (processor_shape == 'square-ice') then ! better for bipolar ice
! Subprogram not used       jguess = nint(square)
! Subprogram not used       iguess = nprocs/jguess
! Subprogram not used    elseif (processor_shape == 'slenderX1') then ! 1 proc in y direction
! Subprogram not used       jguess = 1
! Subprogram not used       iguess = nprocs/jguess
! Subprogram not used    elseif (processor_shape == 'slenderX2') then ! 2 proc in y direction
! Subprogram not used       jguess = min(2, nprocs)
! Subprogram not used       iguess = nprocs/jguess
! Subprogram not used    else                                  ! abort
! Subprogram not used       call abort_ice('ice: processor_shape not supported, '//trim(processor_shape))
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  try various decompositions to find the best
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    proc_loop: do
! Subprogram not used    if (processor_shape == 'square-pop') then
! Subprogram not used       jguess = nprocs/iguess
! Subprogram not used    else
! Subprogram not used       iguess = nprocs/jguess
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used       if (iguess*jguess == nprocs) then ! valid decomp
! Subprogram not used 
! Subprogram not used          !*** if the blocks can be evenly distributed, it is a
! Subprogram not used          !*** good decomposition
! Subprogram not used          if (mod(nblocks_x,iguess) == 0 .and. &
! Subprogram not used              mod(nblocks_y,jguess) == 0) then
! Subprogram not used             nprocs_x = iguess
! Subprogram not used             nprocs_y = jguess
! Subprogram not used             exit proc_loop
! Subprogram not used 
! Subprogram not used          !*** if the blocks can be evenly distributed in a
! Subprogram not used          !*** transposed direction, it is a good decomposition
! Subprogram not used          else if (mod(nblocks_x,jguess) == 0 .and. &
! Subprogram not used                 mod(nblocks_y,iguess) == 0) then
! Subprogram not used             nprocs_x = jguess
! Subprogram not used             nprocs_y = iguess
! Subprogram not used             exit proc_loop
! Subprogram not used 
! Subprogram not used          !*** A valid decomposition, but keep searching for
! Subprogram not used          !***  a better one
! Subprogram not used          else
! Subprogram not used             if (nprocs_x == 0) then
! Subprogram not used                nprocs_x = iguess
! Subprogram not used                nprocs_y = jguess
! Subprogram not used             endif
! Subprogram not used             if (processor_shape == 'square-pop') then
! Subprogram not used                iguess = iguess - 1
! Subprogram not used                if (iguess == 0) then
! Subprogram not used                   exit proc_loop
! Subprogram not used                else
! Subprogram not used                   cycle proc_loop
! Subprogram not used                endif
! Subprogram not used             else
! Subprogram not used                jguess = jguess - 1
! Subprogram not used                if (jguess == 0) then
! Subprogram not used                   exit proc_loop
! Subprogram not used                else
! Subprogram not used                   cycle proc_loop
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       else ! invalid decomp - keep trying
! Subprogram not used 
! Subprogram not used          if (processor_shape == 'square-pop') then
! Subprogram not used             iguess = iguess - 1
! Subprogram not used             if (iguess == 0) then
! Subprogram not used                exit proc_loop
! Subprogram not used             else
! Subprogram not used                cycle proc_loop
! Subprogram not used             endif
! Subprogram not used          else
! Subprogram not used             jguess = jguess - 1
! Subprogram not used             if (jguess == 0) then
! Subprogram not used                exit proc_loop
! Subprogram not used             else
! Subprogram not used                cycle proc_loop
! Subprogram not used             endif
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    end do proc_loop
! Subprogram not used 
! Subprogram not used    if (nprocs_x == 0) then
! Subprogram not used       call abort_ice('ice: Unable to find 2d processor config')
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (my_task == master_task) then
! Subprogram not used      write(nu_diag,'(a23,i4,a3,i4)') '  Processors (X x Y) = ', &
! Subprogram not used                                         nprocs_x,' x ',nprocs_y
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end subroutine proc_decomposition

!**********************************************************************
!BOP
! !IROUTINE: ice_distributionDestroy
! !INTERFACE:

! Subprogram not used  subroutine ice_distributionDestroy(distribution)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This routine destroys a defined distribution by deallocating
! Subprogram not used !  all memory associated with the distribution.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type (distrb), intent(inout) :: &
! Subprogram not used       distribution          ! distribution to destroy
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: istat  ! status flag for deallocate
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  reset scalars
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    distribution%nprocs       = 0
! Subprogram not used    distribution%communicator   = 0
! Subprogram not used    distribution%numLocalBlocks = 0
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  deallocate arrays
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    deallocate(distribution%blockLocation, stat=istat)
! Subprogram not used    deallocate(distribution%blockLocalID , stat=istat)
! Subprogram not used    deallocate(distribution%blockGlobalID, stat=istat)
! Subprogram not used    deallocate(distribution%blockCnt, stat=istat)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end subroutine ice_distributionDestroy

!***********************************************************************
!BOP
! !IROUTINE: ice_distributionGet
! !INTERFACE:

 subroutine ice_distributionGet(distribution,&
                            nprocs, communicator, numLocalBlocks, &
                            blockLocation, blockLocalID, blockGlobalID)


! !DESCRIPTION:
!  This routine extracts information from a distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

! !OUTPUT PARAMETERS:

      integer (int_kind), intent(out), optional ::   &
         nprocs          ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (int_kind), dimension(:), pointer, optional :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  depending on which optional arguments are present, extract the
!  requested info
!
!-----------------------------------------------------------------------

   if (present(nprocs))       nprocs       = distribution%nprocs
   if (present(communicator))   communicator   = distribution%communicator
   if (present(numLocalBlocks)) numLocalBlocks = distribution%numLocalBlocks

   if (present(blockLocation)) then
      if (associated(distribution%blockLocation)) then
         blockLocation => distribution%blockLocation
      else
        call abort_ice( &
            'ice_distributionGet: blockLocation not allocated')
         return
      endif
   endif

   if (present(blockLocalID)) then
      if (associated(distribution%blockLocalID)) then
         blockLocalID = distribution%blockLocalID
      else
        call abort_ice( &
            'ice_distributionGet: blockLocalID not allocated')
         return
      endif
   endif

   if (present(blockGlobalID)) then
      if (associated(distribution%blockGlobalID)) then
         blockGlobalID = distribution%blockGlobalID
      else
        call abort_ice( &
            'ice_distributionGet: blockGlobalID not allocated')
         return
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionGet

!***********************************************************************
!BOP
! !IROUTINE: ice_distributionGetBlockLoc
! !INTERFACE:

 subroutine ice_distributionGetBlockLoc(distribution, blockID, &
                                        processor, localID)


! !DESCRIPTION:
!  Given a distribution of blocks and a global block ID, return
!  the processor and local index for the block.  A zero for both
!  is returned in the case that the block has been eliminated from
!  the distribution (i.e. has no active points).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (int_kind), intent(in) :: &
      blockID                ! global block id for which location requested

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) ::  &
      processor,            &! processor on which block resides
      localID                ! local index for this block on this proc

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for valid blockID
!
!-----------------------------------------------------------------------

   if (blockID < 0 .or. blockID > nblocks_tot) then
     call abort_ice( &
         'ice_distributionGetBlockLoc: invalid block id')
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the location from the distribution data structure
!
!-----------------------------------------------------------------------

   processor = distribution%blockLocation(blockID)
   localID   = distribution%blockLocalID (blockID)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionGetBlockLoc

!***********************************************************************
!BOP
! !IROUTINE: ice_distributionGetBlockID
! !INTERFACE:

 subroutine ice_distributionGetBlockID(distribution, localID, &
                                       blockID)


! !DESCRIPTION:
!  Given a distribution of blocks and a local block index, return
!  the global block id for the block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (int_kind), intent(in) ::  &
      localID                ! local index for this block on this proc

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      blockID                ! global block id for this local block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for valid localID
!
!-----------------------------------------------------------------------

   if (localID < 0 .or. localID > distribution%numLocalBlocks) then
     call abort_ice( &
         'ice_distributionGetBlockID: invalid local id')
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the global ID from the distribution data structure
!
!-----------------------------------------------------------------------

   blockID   = distribution%blockGlobalID (localID)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionGetBlockID

!***********************************************************************
!BOP
! !IROUTINE: create_distrb_roundrobin
! !INTERFACE:

 function create_distrb_roundrobin(nprocs, workPerBlock) result(newDistrb)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a simple roundrobin algorithm. Mean for prescribed ice or
!  standalone CAM mode.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j,                  &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID                 ! block location on this processor

   integer (int_kind), dimension(:), allocatable :: &
      proc_tmp           ! temp processor id
   
   integer (int_kind) :: pid,n

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (newDistrb%blockLocation(nblocks_tot), &
             newDistrb%blockLocalID (nblocks_tot), stat=istat)

   allocate (newDistrb%blockCnt(nprocs))
!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   processor = 0
   globalID = 0
   proc_tmp = 0

   allocate(newDistrb%blockIndex(nprocs,max_blocks))
   newDistrb%blockIndex(:,:) = 0

   do j=1,nblocks_y
   do i=1,nblocks_x
      
      globalID = globalID + 1

      if (workPerBlock(globalID) /= 0) then
         processor = mod(processor,nprocs) + 1
         proc_tmp(processor) = proc_tmp(processor) + 1
         localID = proc_tmp(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
         newDistrb%blockIndex(processor,localID) = globalID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif

   end do
   end do

   newDistrb%numLocalBlocks = proc_tmp(my_task+1)
   newDistrb%blockCnt(:) = proc_tmp(:)
   deallocate(proc_tmp)

!   write(nu_diag,*) 'my_task,newDistrb%numLocalBlocks',&
!      my_task,newDistrb%numLocalBlocks

!----------------------------------------------------------------------
!
!  now store the local info
!
!----------------------------------------------------------------------

   globalID = 0

   if (newDistrb%numLocalBlocks > 0) then
      allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
                stat=istat)

      processor = my_task + 1
      do localID = 1,newDistrb%numLocalBlocks
         newDistrb%blockGlobalID (localID) = newDistrb%blockIndex(processor,&
                                             localID)
      enddo
   endif

!----------------------------------------------------------------------
!EOC
 end function create_distrb_roundrobin
 
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_blkcart
! !INTERFACE:

! Subprogram not used  function create_distrb_blkcart(nprocs, workPerBlock) result(newDistrb)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This function creates a distribution of blocks across processors
! Subprogram not used !  using a simple blkcart algorithm. Mean for prescribed ice or
! Subprogram not used !  standalone CAM mode.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       nprocs            ! number of processors in this distribution
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), intent(in) :: &
! Subprogram not used       workPerBlock        ! amount of work per block
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type (distrb) :: &
! Subprogram not used       newDistrb           ! resulting structure describing Cartesian
! Subprogram not used                           !  distribution of blocks
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i, j, i2, j2,          &! dummy loop indices
! Subprogram not used       istat,                 &! status flag for allocation
! Subprogram not used       iblock, jblock,        &!
! Subprogram not used       processor,             &! processor position in cartesian decomp
! Subprogram not used       globalID,              &! global block ID
! Subprogram not used       localID,               &! block location on this processor
! Subprogram not used       blktogether,           &! number of blocks together
! Subprogram not used       cnt                     ! counter
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used       proc_tmp           ! temp processor id
! Subprogram not used    
! Subprogram not used    integer (int_kind) :: pid,n
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  create communicator for this distribution
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call create_communicator(newDistrb%communicator, nprocs)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  try to find best processor arrangement
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    newDistrb%nprocs = nprocs
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  allocate space for decomposition
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate (newDistrb%blockLocation(nblocks_tot), &
! Subprogram not used              newDistrb%blockLocalID (nblocks_tot), stat=istat)
! Subprogram not used 
! Subprogram not used    allocate (newDistrb%blockCnt(nprocs))
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  distribute blocks linearly across processors in each direction
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate(proc_tmp(nprocs))
! Subprogram not used    proc_tmp = 0
! Subprogram not used 
! Subprogram not used    allocate(newDistrb%blockIndex(nprocs,max_blocks))
! Subprogram not used    newDistrb%blockIndex(:,:) = 0
! Subprogram not used 
! Subprogram not used    blktogether = max(1,nint(float(nblocks_x*nblocks_y)/float(4*nprocs)))
! Subprogram not used 
! Subprogram not used    ! --- two phases, resetnt processor and cnt for each phase
! Subprogram not used    ! --- phase 1 is south to north, east to west on the left half of the domain
! Subprogram not used    ! --- phase 2 is north to south, east to west on the right half of the domain
! Subprogram not used 
! Subprogram not used    if (mod(nblocks_x,2) /= 0) then
! Subprogram not used       call abort_ice( &
! Subprogram not used          'create_distrb_blkcart: nblocks_x not divisible by 2')
! Subprogram not used       return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    do n=1,2
! Subprogram not used    processor = 1
! Subprogram not used    cnt = 0
! Subprogram not used    do j2=1,nblocks_y
! Subprogram not used    do i2=1,nblocks_x/2
! Subprogram not used       
! Subprogram not used       if (n == 1) then
! Subprogram not used          i = i2
! Subprogram not used          j = j2
! Subprogram not used       else
! Subprogram not used          i = nblocks_x/2 + i2
! Subprogram not used          j = nblocks_y - j2 + 1
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       globalID = (j-1)*nblocks_x + i
! Subprogram not used       if (cnt >= blktogether) then
! Subprogram not used          processor = mod(processor,nprocs) + 1
! Subprogram not used          cnt = 0
! Subprogram not used       endif
! Subprogram not used       cnt = cnt + 1
! Subprogram not used 
! Subprogram not used       if (workPerBlock(globalID) /= 0) then
! Subprogram not used          proc_tmp(processor) = proc_tmp(processor) + 1
! Subprogram not used          localID = proc_tmp(processor)
! Subprogram not used          newDistrb%blockLocation(globalID) = processor
! Subprogram not used          newDistrb%blockLocalID (globalID) = localID
! Subprogram not used          newDistrb%blockIndex(processor,localID) = globalID
! Subprogram not used       else  ! no work - eliminate block from distribution
! Subprogram not used          newDistrb%blockLocation(globalID) = 0
! Subprogram not used          newDistrb%blockLocalID (globalID) = 0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used    end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    newDistrb%numLocalBlocks = proc_tmp(my_task+1)
! Subprogram not used    newDistrb%blockCnt(:) = proc_tmp(:)
! Subprogram not used    deallocate(proc_tmp)
! Subprogram not used 
! Subprogram not used !   write(nu_diag,*) 'my_task,newDistrb%numLocalBlocks',&
! Subprogram not used !      my_task,newDistrb%numLocalBlocks
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now store the local info
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    globalID = 0
! Subprogram not used 
! Subprogram not used    if (newDistrb%numLocalBlocks > 0) then
! Subprogram not used       allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
! Subprogram not used                 stat=istat)
! Subprogram not used 
! Subprogram not used       processor = my_task + 1
! Subprogram not used       do localID = 1,newDistrb%numLocalBlocks
! Subprogram not used          newDistrb%blockGlobalID (localID) = newDistrb%blockIndex(processor,&
! Subprogram not used                                              localID)
! Subprogram not used       enddo
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used  end function create_distrb_blkcart
 
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_blkrobin
! !INTERFACE:

! Subprogram not used  function create_distrb_blkrobin(nprocs, workPerBlock) result(newDistrb)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This function creates a distribution of blocks across processors
! Subprogram not used !  using a simple blkrobin algorithm. Mean for prescribed ice or
! Subprogram not used !  standalone CAM mode.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       nprocs            ! number of processors in this distribution
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), intent(in) :: &
! Subprogram not used       workPerBlock        ! amount of work per block
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type (distrb) :: &
! Subprogram not used       newDistrb           ! resulting structure describing Cartesian
! Subprogram not used                           !  distribution of blocks
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i, j,                  &! dummy loop indices
! Subprogram not used       istat,                 &! status flag for allocation
! Subprogram not used       iblock, jblock,        &!
! Subprogram not used       mblocks,               &! estimate of max blocks per pe
! Subprogram not used       processor,             &! processor position in cartesian decomp
! Subprogram not used       globalID,              &! global block ID
! Subprogram not used       localID                 ! block location on this processor
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used       proc_tmp           ! temp processor id
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:), allocatable :: &
! Subprogram not used       bfree              ! map of assigned blocks
! Subprogram not used    
! Subprogram not used    integer (int_kind) :: pid,n,cnt, blktogether, i2, j2
! Subprogram not used    integer (int_kind) :: totblocks, nchunks
! Subprogram not used    logical (log_kind) :: keepgoing
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  create communicator for this distribution
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call create_communicator(newDistrb%communicator, nprocs)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  try to find best processor arrangement
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    newDistrb%nprocs = nprocs
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  allocate space for decomposition
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate (newDistrb%blockLocation(nblocks_tot), &
! Subprogram not used              newDistrb%blockLocalID (nblocks_tot), stat=istat)
! Subprogram not used 
! Subprogram not used    allocate (newDistrb%blockCnt(nprocs))
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  distribute blocks linearly across processors in each direction
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate(proc_tmp(nprocs))
! Subprogram not used    processor = 0
! Subprogram not used    globalID = 0
! Subprogram not used    proc_tmp = 0
! Subprogram not used 
! Subprogram not used    allocate(newDistrb%blockIndex(nprocs,max_blocks))
! Subprogram not used    newDistrb%blockIndex(:,:) = 0
! Subprogram not used 
! Subprogram not used    allocate(bfree(nblocks_x*nblocks_y))
! Subprogram not used    bfree=.true.
! Subprogram not used 
! Subprogram not used    totblocks = 0
! Subprogram not used    do j=1,nblocks_y
! Subprogram not used    do i=1,nblocks_x
! Subprogram not used       globalID = (j-1)*nblocks_x + i
! Subprogram not used       if (workPerBlock(globalID) /= 0) then
! Subprogram not used          totblocks=totblocks+1
! Subprogram not used       else  ! no work - eliminate block from distribution
! Subprogram not used          bfree(globalID) = .false.
! Subprogram not used          newDistrb%blockLocation(globalID) = 0
! Subprogram not used          newDistrb%blockLocalID (globalID) = 0
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    mblocks = totblocks/nprocs
! Subprogram not used    if (mod(totblocks,nprocs) > 0) mblocks=mblocks+1
! Subprogram not used 
! Subprogram not used    blktogether = max(1,nint(float(totblocks)/float(6*nprocs)))
! Subprogram not used 
! Subprogram not used !   write(nu_diag,*) 'ice_distrb_blkrobin totblocks = ',totblocks,nblocks_y*nblocks_x
! Subprogram not used  
! Subprogram not used    !------------------------------
! Subprogram not used    ! southern group of blocks
! Subprogram not used    !   weave back and forth in i vs j
! Subprogram not used    !   go south to north, low - high pes
! Subprogram not used    !------------------------------
! Subprogram not used 
! Subprogram not used    processor=1
! Subprogram not used    cnt = 0
! Subprogram not used    keepgoing = .true.
! Subprogram not used    do j=1,nblocks_y
! Subprogram not used    do i=1,nblocks_x
! Subprogram not used       if (mod(j,2) == 0) then
! Subprogram not used          i2 = nblocks_x - i + 1
! Subprogram not used       else
! Subprogram not used          i2 = i
! Subprogram not used       endif
! Subprogram not used       globalID = (j-1)*nblocks_x + i2
! Subprogram not used       if (cnt >= blktogether) then
! Subprogram not used          processor = mod(processor,nprocs) + 1
! Subprogram not used          cnt = 0
! Subprogram not used          if (processor == 1) keepgoing = .false.
! Subprogram not used       endif
! Subprogram not used !      write(nu_diag,'(a,6i7,l2)') 'tcx ',i,j,globalID,cnt,blktogether,processor,keepgoing
! Subprogram not used 
! Subprogram not used       if (keepgoing) then
! Subprogram not used          if (bfree(globalID)) then
! Subprogram not used          if (workPerBlock(globalID) /= 0) then
! Subprogram not used             proc_tmp(processor) = proc_tmp(processor) + 1
! Subprogram not used             localID = proc_tmp(processor)
! Subprogram not used             newDistrb%blockLocation(globalID) = processor
! Subprogram not used             newDistrb%blockLocalID (globalID) = localID
! Subprogram not used             newDistrb%blockIndex(processor,localID) = globalID
! Subprogram not used             cnt = cnt + 1
! Subprogram not used             totblocks = totblocks-1
! Subprogram not used             bfree(globalID) = .false.
! Subprogram not used 
! Subprogram not used          else  ! no work - eliminate block from distribution
! Subprogram not used             bfree(globalID) = .false.
! Subprogram not used             newDistrb%blockLocation(globalID) = 0
! Subprogram not used             newDistrb%blockLocalID (globalID) = 0
! Subprogram not used          endif
! Subprogram not used          endif  ! bfree
! Subprogram not used       endif
! Subprogram not used    end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !   write(nu_diag,*) 'ice_distrb_blkrobin totblocks left after southern = ',totblocks
! Subprogram not used 
! Subprogram not used    !------------------------------
! Subprogram not used    ! northern group of blocks
! Subprogram not used    !   weave back and forth in i vs j
! Subprogram not used    !   go north to south, high - low pes
! Subprogram not used    !------------------------------
! Subprogram not used 
! Subprogram not used    processor=nprocs
! Subprogram not used    cnt = 0
! Subprogram not used    keepgoing = .true.
! Subprogram not used    do j=nblocks_y,1,-1
! Subprogram not used    do i=1,nblocks_x
! Subprogram not used       if (mod(j,2) == 1) then
! Subprogram not used          i2 = nblocks_x - i + 1
! Subprogram not used       else
! Subprogram not used          i2 = i
! Subprogram not used       endif
! Subprogram not used       globalID = (j-1)*nblocks_x + i2
! Subprogram not used       if (cnt >= blktogether) then
! Subprogram not used          processor = mod(processor+nprocs-2,nprocs) + 1
! Subprogram not used          cnt = 0
! Subprogram not used          if (processor == nprocs) keepgoing = .false.
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (keepgoing) then
! Subprogram not used          if (bfree(globalID)) then
! Subprogram not used          if (workPerBlock(globalID) /= 0) then
! Subprogram not used             proc_tmp(processor) = proc_tmp(processor) + 1
! Subprogram not used             localID = proc_tmp(processor)
! Subprogram not used             newDistrb%blockLocation(globalID) = processor
! Subprogram not used             newDistrb%blockLocalID (globalID) = localID
! Subprogram not used             newDistrb%blockIndex(processor,localID) = globalID
! Subprogram not used             cnt = cnt + 1
! Subprogram not used             totblocks = totblocks - 1
! Subprogram not used             bfree(globalID) = .false.
! Subprogram not used 
! Subprogram not used          else  ! no work - eliminate block from distribution
! Subprogram not used             bfree(globalID) = .false.
! Subprogram not used             newDistrb%blockLocation(globalID) = 0
! Subprogram not used             newDistrb%blockLocalID (globalID) = 0
! Subprogram not used          endif
! Subprogram not used          endif  ! bfree
! Subprogram not used       endif
! Subprogram not used    end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !   write(nu_diag,*) 'ice_distrb_blkrobin totblocks left after northern = ',totblocks
! Subprogram not used 
! Subprogram not used    !------------------------------
! Subprogram not used    ! central group of blocks
! Subprogram not used    !   weave back and forth in i vs j
! Subprogram not used    !   go north to south, low - high / low - high pes
! Subprogram not used    !------------------------------
! Subprogram not used 
! Subprogram not used    nchunks = 2*nprocs
! Subprogram not used    blktogether = max(1,nint(float(totblocks)/float(nchunks)))
! Subprogram not used    processor=1
! Subprogram not used    cnt = 0
! Subprogram not used    do j=nblocks_y,1,-1
! Subprogram not used    do i=1,nblocks_x
! Subprogram not used       if (mod(j,2) == 1) then
! Subprogram not used          i2 = nblocks_x - i + 1
! Subprogram not used       else
! Subprogram not used          i2 = i
! Subprogram not used       endif
! Subprogram not used       globalID = (j-1)*nblocks_x + i2
! Subprogram not used       if (totblocks > 0) then
! Subprogram not used       do while (proc_tmp(processor) >= mblocks .or. cnt >= blktogether)
! Subprogram not used          nchunks = nchunks - 1
! Subprogram not used          if (nchunks == 0) then
! Subprogram not used             blktogether = 1
! Subprogram not used          else
! Subprogram not used             blktogether = max(1,nint(float(totblocks)/float(nchunks)))
! Subprogram not used          endif
! Subprogram not used          cnt = 0
! Subprogram not used          processor = mod(processor,nprocs) + 1
! Subprogram not used       enddo
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used !      write(nu_diag,*) 'ice_distrb_blkrobin central ',i,j,totblocks,cnt,nchunks,blktogether,processor
! Subprogram not used 
! Subprogram not used       if (bfree(globalID)) then
! Subprogram not used       if (workPerBlock(globalID) /= 0) then
! Subprogram not used          proc_tmp(processor) = proc_tmp(processor) + 1
! Subprogram not used          localID = proc_tmp(processor)
! Subprogram not used          newDistrb%blockLocation(globalID) = processor
! Subprogram not used          newDistrb%blockLocalID (globalID) = localID
! Subprogram not used          newDistrb%blockIndex(processor,localID) = globalID
! Subprogram not used          cnt = cnt + 1
! Subprogram not used          totblocks = totblocks-1
! Subprogram not used          bfree(globalID) = .false.
! Subprogram not used 
! Subprogram not used       else  ! no work - eliminate block from distribution
! Subprogram not used          bfree(globalID) = .false.
! Subprogram not used          newDistrb%blockLocation(globalID) = 0
! Subprogram not used          newDistrb%blockLocalID (globalID) = 0
! Subprogram not used       endif
! Subprogram not used       endif  ! bfree
! Subprogram not used    end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    newDistrb%numLocalBlocks = proc_tmp(my_task+1)
! Subprogram not used    newDistrb%blockCnt(:) = proc_tmp(:)
! Subprogram not used    deallocate(proc_tmp)
! Subprogram not used    deallocate(bfree)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now store the local info
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    globalID = 0
! Subprogram not used 
! Subprogram not used    if (newDistrb%numLocalBlocks > 0) then
! Subprogram not used       allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
! Subprogram not used                 stat=istat)
! Subprogram not used 
! Subprogram not used       processor = my_task + 1
! Subprogram not used       do localID = 1,newDistrb%numLocalBlocks
! Subprogram not used          newDistrb%blockGlobalID (localID) = newDistrb%blockIndex(processor,&
! Subprogram not used                                              localID)
! Subprogram not used       enddo
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used  end function create_distrb_blkrobin
 
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_cart
! !INTERFACE:

! Subprogram not used  function create_distrb_cart(nprocs, workPerBlock) result(newDistrb)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This function creates a distribution of blocks across processors
! Subprogram not used !  using a 2-d Cartesian distribution.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       nprocs            ! number of processors in this distribution
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), intent(in) :: &
! Subprogram not used       workPerBlock        ! amount of work per block
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type (distrb) :: &
! Subprogram not used       newDistrb           ! resulting structure describing Cartesian
! Subprogram not used                           !  distribution of blocks
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i, j,                  &! dummy loop indices
! Subprogram not used       istat,                 &! status flag for allocation
! Subprogram not used       iblock, jblock,        &!
! Subprogram not used       is, ie, js, je,        &! start, end block indices for each proc
! Subprogram not used       processor,             &! processor position in cartesian decomp
! Subprogram not used       globalID,              &! global block ID
! Subprogram not used       localID,               &! block location on this processor
! Subprogram not used       nprocsX,             &! num of procs in x for global domain
! Subprogram not used       nprocsY,             &! num of procs in y for global domain
! Subprogram not used       numBlocksXPerProc,     &! num of blocks per processor in x
! Subprogram not used       numBlocksYPerProc       ! num of blocks per processor in y
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used       proc_tmp           ! temp processor id
! Subprogram not used    
! Subprogram not used    integer (int_kind) :: pid,n
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  create communicator for this distribution
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call create_communicator(newDistrb%communicator, nprocs)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  try to find best processor arrangement
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    newDistrb%nprocs = nprocs
! Subprogram not used 
! Subprogram not used    ! This assumes that we really do know what we are doing.
! Subprogram not used    if (processor_shape == 'blocks') then
! Subprogram not used       nprocsX = nblocks_x / max_blocks
! Subprogram not used       nprocsY = nblocks_y / max_blocks
! Subprogram not used    else
! Subprogram not used       call proc_decomposition(nprocs, nprocsX, nprocsY)
! Subprogram not used    endif
! Subprogram not used                                   
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  allocate space for decomposition
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate (newDistrb%blockLocation(nblocks_tot), &
! Subprogram not used              newDistrb%blockLocalID (nblocks_tot), stat=istat)
! Subprogram not used 
! Subprogram not used    allocate (newDistrb%blockCnt(nprocs))
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  distribute blocks linearly across processors in each direction
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    numBlocksXPerProc = (nblocks_x-1)/nprocsX + 1
! Subprogram not used    numBlocksYPerProc = (nblocks_y-1)/nprocsY + 1
! Subprogram not used 
! Subprogram not used    do j=1,nprocsY
! Subprogram not used    do i=1,nprocsX
! Subprogram not used       processor = (j-1)*nprocsX + i    ! number the processors 
! Subprogram not used                                          ! left to right, bot to top
! Subprogram not used 
! Subprogram not used       is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
! Subprogram not used       ie =  i   *numBlocksXPerProc       ! ending   block in i
! Subprogram not used       if (ie > nblocks_x) ie = nblocks_x
! Subprogram not used       js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
! Subprogram not used       je =  j   *numBlocksYPerProc       ! ending   block in j
! Subprogram not used       if (je > nblocks_y) je = nblocks_y
! Subprogram not used 
! Subprogram not used       localID        = 0  ! initialize counter for local index
! Subprogram not used       do jblock = js,je
! Subprogram not used       do iblock = is,ie
! Subprogram not used          globalID = (jblock - 1)*nblocks_x + iblock
! Subprogram not used          if (workPerBlock(globalID) /= 0) then
! Subprogram not used             localID = localID + 1
! Subprogram not used             newDistrb%blockLocation(globalID) = processor
! Subprogram not used             newDistrb%blockLocalID (globalID) = localID
! Subprogram not used          else  ! no work - eliminate block from distribution
! Subprogram not used             newDistrb%blockLocation(globalID) = 0
! Subprogram not used             newDistrb%blockLocalID (globalID) = 0
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       ! if this is the local processor, set number of local blocks
! Subprogram not used       if (my_task == processor - 1) then
! Subprogram not used          newDistrb%numLocalBlocks = localID
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now store the local info
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (newDistrb%numLocalBlocks > 0) then
! Subprogram not used       allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
! Subprogram not used                 stat=istat)
! Subprogram not used 
! Subprogram not used       do j=1,nprocsY
! Subprogram not used       do i=1,nprocsX
! Subprogram not used          processor = (j-1)*nprocsX + i
! Subprogram not used 
! Subprogram not used          if (processor == my_task + 1) then
! Subprogram not used             is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
! Subprogram not used             ie =  i   *numBlocksXPerProc       ! ending   block in i
! Subprogram not used             if (ie > nblocks_x) ie = nblocks_x
! Subprogram not used             js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
! Subprogram not used             je =  j   *numBlocksYPerProc       ! ending   block in j
! Subprogram not used             if (je > nblocks_y) je = nblocks_y
! Subprogram not used 
! Subprogram not used             localID        = 0  ! initialize counter for local index
! Subprogram not used             do jblock = js,je
! Subprogram not used             do iblock = is,ie
! Subprogram not used                globalID = (jblock - 1)*nblocks_x + iblock
! Subprogram not used                if (workPerBlock(globalID) /= 0) then
! Subprogram not used                   localID = localID + 1
! Subprogram not used                   newDistrb%blockGlobalID (localID) = globalID
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    allocate(proc_tmp(nprocs))
! Subprogram not used    proc_tmp = 0
! Subprogram not used 
! Subprogram not used    allocate(newDistrb%blockIndex(nprocs,max_blocks))
! Subprogram not used    newDistrb%blockIndex(:,:) = 0
! Subprogram not used 
! Subprogram not used    do n=1,nblocks_tot
! Subprogram not used       pid = newDistrb%blockLocation(n)
! Subprogram not used       if(pid>0) then
! Subprogram not used         proc_tmp(pid) = proc_tmp(pid) + 1
! Subprogram not used         if(proc_tmp(pid) <= max_blocks) then 
! Subprogram not used             newDistrb%blockIndex(pid,proc_tmp(pid)) = n
! Subprogram not used         endif
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    newDistrb%blockCnt(:) = proc_tmp(:)
! Subprogram not used    deallocate(proc_tmp)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function create_distrb_cart

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_rake
! !INTERFACE:

! Subprogram not used  function create_distrb_rake(nprocs, workPerBlock) result(newDistrb)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This  function distributes blocks across processors in a
! Subprogram not used !  load-balanced manner based on the amount of work per block.
! Subprogram not used !  A rake algorithm is used in which the blocks are first distributed
! Subprogram not used !  in a Cartesian distribution and then a rake is applied in each
! Subprogram not used !  Cartesian direction.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       nprocs                ! number of processors in this distribution
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), intent(in) :: &
! Subprogram not used       workPerBlock        ! amount of work per block
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type (distrb) :: &
! Subprogram not used       newDistrb           ! resulting structure describing
! Subprogram not used                           ! load-balanced distribution of blocks
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) ::    &
! Subprogram not used       i,j,n              ,&! dummy loop indices
! Subprogram not used       pid                ,&! dummy for processor id
! Subprogram not used       istat              ,&! status flag for allocates
! Subprogram not used       localBlock         ,&! local block position on processor
! Subprogram not used       numOcnBlocks       ,&! number of ocean blocks
! Subprogram not used       maxWork            ,&! max amount of work in any block
! Subprogram not used       nprocsX          ,&! num of procs in x for global domain
! Subprogram not used       nprocsY            ! num of procs in y for global domain
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used       priority           ,&! priority for moving blocks
! Subprogram not used       workTmp            ,&! work per row or column for rake algrthm
! Subprogram not used       procTmp              ! temp processor id for rake algrthm
! Subprogram not used 
! Subprogram not used    type (distrb) :: dist  ! temp hold distribution
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  first set up as Cartesian distribution
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dist = create_distrb_cart(nprocs, workPerBlock)
! Subprogram not used                                     
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if the number of blocks is close to the number of processors,
! Subprogram not used !  only do a 1-d rake on the entire distribution
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    numOcnBlocks = count(workPerBlock /= 0)
! Subprogram not used 
! Subprogram not used    if (numOcnBlocks <= 2*nprocs) then
! Subprogram not used 
! Subprogram not used       allocate(priority(nblocks_tot), stat=istat)
! Subprogram not used 
! Subprogram not used       !*** initialize priority array
! Subprogram not used 
! Subprogram not used       do j=1,nblocks_y
! Subprogram not used       do i=1,nblocks_x
! Subprogram not used          n=(j-1)*nblocks_x + i
! Subprogram not used          if (workPerBlock(n) > 0) then
! Subprogram not used             priority(n) = maxWork + n - workPerBlock(n)
! Subprogram not used          else
! Subprogram not used             priority(n) = 0
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       allocate(workTmp(nblocks_tot), procTmp(nblocks_tot), stat=istat)
! Subprogram not used 
! Subprogram not used       workTmp(:) = 0
! Subprogram not used       do i=1,nprocs
! Subprogram not used          procTmp(i) = i
! Subprogram not used          do n=1,nblocks_tot
! Subprogram not used             if (dist%blockLocation(n) == i) then
! Subprogram not used                workTmp(i) = workTmp(i) + workPerBlock(n)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       call ice_distributionRake (workTmp, procTmp, workPerBlock, &
! Subprogram not used                                  priority, dist)
! Subprogram not used 
! Subprogram not used       deallocate(workTmp, procTmp, stat=istat)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  otherwise re-distribute blocks using a rake in each direction
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used       maxWork = maxval(workPerBlock)
! Subprogram not used 
! Subprogram not used       call proc_decomposition(dist%nprocs, nprocsX, nprocsY)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !     load-balance using a rake algorithm in the x-direction first
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       allocate(priority(nblocks_tot), stat=istat)
! Subprogram not used 
! Subprogram not used       !*** set highest priority such that eastern-most blocks
! Subprogram not used       !*** and blocks with the least amount of work are
! Subprogram not used       !*** moved first
! Subprogram not used 
! Subprogram not used       do j=1,nblocks_y
! Subprogram not used       do i=1,nblocks_x
! Subprogram not used          n=(j-1)*nblocks_x + i
! Subprogram not used          if (workPerBlock(n) > 0) then
! Subprogram not used             priority(n) = (maxWork + 1)*(nblocks_x + i) - &
! Subprogram not used                           workPerBlock(n)
! Subprogram not used          else
! Subprogram not used             priority(n) = 0
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       allocate(workTmp(nprocsX), procTmp(nprocsX), stat=istat)
! Subprogram not used 
! Subprogram not used       do j=1,nprocsY
! Subprogram not used 
! Subprogram not used          workTmp(:) = 0
! Subprogram not used          do i=1,nprocsX
! Subprogram not used             pid = (j-1)*nprocsX + i
! Subprogram not used             procTmp(i) = pid
! Subprogram not used             do n=1,nblocks_tot
! Subprogram not used                if (dist%blockLocation(n) == pid) then
! Subprogram not used                   workTmp(i) = workTmp(i) + workPerBlock(n)
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          call ice_distributionRake (workTmp, procTmp, workPerBlock, &
! Subprogram not used                                     priority, dist)
! Subprogram not used       end do
! Subprogram not used    
! Subprogram not used       deallocate(workTmp, procTmp, stat=istat)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !     use a rake algorithm in the y-direction now
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !*** set highest priority for northern-most blocks
! Subprogram not used 
! Subprogram not used       do j=1,nblocks_y
! Subprogram not used       do i=1,nblocks_x
! Subprogram not used          n=(j-1)*nblocks_x + i
! Subprogram not used          if (workPerBlock(n) > 0) then
! Subprogram not used             priority(n) = (maxWork + 1)*(nblocks_y + j) - &
! Subprogram not used                           workPerBlock(n)
! Subprogram not used          else
! Subprogram not used             priority(n) = 0
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       allocate(workTmp(nprocsY), procTmp(nprocsY), stat=istat)
! Subprogram not used 
! Subprogram not used       do i=1,nprocsX
! Subprogram not used 
! Subprogram not used          workTmp(:) = 0
! Subprogram not used          do j=1,nprocsY
! Subprogram not used             pid = (j-1)*nprocsX + i
! Subprogram not used             procTmp(j) = pid
! Subprogram not used             do n=1,nblocks_tot
! Subprogram not used                if (dist%blockLocation(n) == pid) then
! Subprogram not used                   workTmp(j) = workTmp(j) + workPerBlock(n)
! Subprogram not used                endif
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          call ice_distributionRake (workTmp, procTmp, workPerBlock, &
! Subprogram not used                                     priority, dist)
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       deallocate(workTmp, procTmp, priority, stat=istat)
! Subprogram not used 
! Subprogram not used    endif  ! 1d or 2d rake
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  create new distribution with info extracted from the temporary
! Subprogram not used !  distribution
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    newDistrb%nprocs     = nprocs
! Subprogram not used    newDistrb%communicator = dist%communicator
! Subprogram not used 
! Subprogram not used    allocate(newDistrb%blockLocation(nblocks_tot), &
! Subprogram not used             newDistrb%blockLocalID(nblocks_tot), stat=istat)
! Subprogram not used    allocate (newDistrb%blockCnt(nprocs))
! Subprogram not used 
! Subprogram not used    allocate(procTmp(nprocs), stat=istat)
! Subprogram not used 
! Subprogram not used    allocate(newDistrb%blockIndex(nprocs,max_blocks))
! Subprogram not used    newDistrb%blockIndex(:,:) = 0
! Subprogram not used 
! Subprogram not used    procTmp = 0
! Subprogram not used    do n=1,nblocks_tot
! Subprogram not used       pid = dist%blockLocation(n)  ! processor id
! Subprogram not used       newDistrb%blockLocation(n) = pid
! Subprogram not used 
! Subprogram not used       if (pid > 0) then
! Subprogram not used          procTmp(pid) = procTmp(pid) + 1
! Subprogram not used          newDistrb%blockLocalID (n) = procTmp(pid)
! Subprogram not used 	 if(procTmp(pid) <= max_blocks) then 
! Subprogram not used             newDistrb%blockIndex(pid,procTmp(pid)) = n
! Subprogram not used          endif
! Subprogram not used       else
! Subprogram not used          newDistrb%blockLocalID (n) = 0
! Subprogram not used       endif
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    newDistrb%numLocalBlocks = procTmp(my_task+1)
! Subprogram not used 
! Subprogram not used    if (minval(procTmp) < 1) then
! Subprogram not used       call abort_ice( &
! Subprogram not used          'create_distrb_rake: processors left with no blocks')
! Subprogram not used       return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    newDistrb%blockCnt(:) = procTmp(:) 
! Subprogram not used 
! Subprogram not used    deallocate(procTmp, stat=istat)
! Subprogram not used 
! Subprogram not used    if (istat > 0) then
! Subprogram not used       call abort_ice( &
! Subprogram not used          'create_distrb_rake: error allocating last procTmp')
! Subprogram not used       return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
! Subprogram not used             stat=istat)
! Subprogram not used 
! Subprogram not used    if (istat > 0) then
! Subprogram not used       call abort_ice( &
! Subprogram not used          'create_distrb_rake: error allocating blockGlobalID')
! Subprogram not used       return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    localBlock = 0
! Subprogram not used    do n=1,nblocks_tot
! Subprogram not used       if (newDistrb%blockLocation(n) == my_task+1) then
! Subprogram not used          localBlock = localBlock + 1
! Subprogram not used          newDistrb%blockGlobalID(localBlock) = n
! Subprogram not used       endif
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionDestroy(dist)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function create_distrb_rake

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_spacecurve
! !INTERFACE:

! Subprogram not used  function create_distrb_spacecurve(nprocs, minBlock, maxBlock, work_per_block,prob_per_block,blockType, bStats, &
! Subprogram not used       FixMaxBlock, maxDil )
! Subprogram not used 
! Subprogram not used ! !Description:
! Subprogram not used !  This function distributes blocks across processors in a
! Subprogram not used !  load-balanced manner using space-filling curves
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  added by J. Dennis 3/10/06
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       nprocs                ! number of processors in this distribution
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: minBlock, maxBlock
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), intent(in) :: &
! Subprogram not used       work_per_block        ! amount of work per block
! Subprogram not used 
! Subprogram not used    real(dbl_kind), dimension(:), intent(in) :: &
! Subprogram not used       prob_per_block        ! probability sea-ice within block 
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), intent(in)  :: &
! Subprogram not used       blockType            ! type of block
! Subprogram not used 
! Subprogram not used     real (dbl_kind),  dimension(:,:), intent(in) :: bStats
! Subprogram not used 
! Subprogram not used     logical, intent(in) :: FixMaxBlock
! Subprogram not used 
! Subprogram not used      real (real_kind) :: maxDil
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type (distrb) :: &
! Subprogram not used       create_distrb_spacecurve  ! resulting structure describing
! Subprogram not used                                 ! load-balanced distribution of blocks
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,k,n              ,&! dummy loop indices
! Subprogram not used       pid                  ,&! dummy for processor id
! Subprogram not used       localID              ,&! local block position on processor
! Subprogram not used       max_work             ,&! max amount of work in any block
! Subprogram not used       nprocs_x             ,&! num of procs in x for global domain
! Subprogram not used       nprocs_y               ! num of procs in y for global domain
! Subprogram not used 
! Subprogram not used    character(char_len) :: fname
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: maxB
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:),allocatable :: &
! Subprogram not used         idxT_i,idxT_j,Lindx,Lindx2       ! Temporary indices for SFC
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:),allocatable :: &
! Subprogram not used         Mesh            ,&!   !arrays to hold Space-filling curve
! Subprogram not used         Mesh2             !
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used         nblocksL,nblocks, &! Number of blocks local and total
! Subprogram not used         ii,extra,tmp1,    &! loop tempories used for
! Subprogram not used         s1,ig              ! partitioning curve
! Subprogram not used 
! Subprogram not used    logical, parameter :: Debug = .FALSE.
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:), allocatable :: &
! Subprogram not used       priority           ,&! priority for moving blocks
! Subprogram not used       work_tmp           ,&! work per row or column for rake algrthm
! Subprogram not used       proc_tmp           ,&! temp processor id for rake algrthm
! Subprogram not used       block_count          ! counter to determine local block indx
! Subprogram not used  
! Subprogram not used    integer (int_kind), allocatable, dimension(:) ::  &
! Subprogram not used           blockLocation,   &! block to processor mapping
! Subprogram not used           distance,        &! location in uncompressed SFC
! Subprogram not used           type_on_curve,   &! type of blocks 
! Subprogram not used           work_on_curve
! Subprogram not used 
! Subprogram not used    type (distrb) :: dist  ! temp hold distribution
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: numIce, minblocks
! Subprogram not used    type (factor_t) :: xdim, ydim
! Subprogram not used    
! Subprogram not used    integer (int_kind) :: it,jj,i2,j2
! Subprogram not used    integer (int_kind) :: curveSize, sb_x, sb_y, itmp,numfac
! Subprogram not used    integer (int_kind) :: subNum,sfcNum
! Subprogram not used    logical (log_kind) :: foundX 
! Subprogram not used    
! Subprogram not used    integer (int_kind) :: ns,ntmp
! Subprogram not used    integer (int_kind) :: numLocalBlocks
! Subprogram not used 
! Subprogram not used    real (dbl_kind), allocatable, dimension(:) :: &
! Subprogram not used        Cost_per_proc , &
! Subprogram not used        Cost_per_block
! Subprogram not used 
! Subprogram not used    real (dbl_kind), allocatable, dimension(:,:) :: cStats  ! block statistics on SFC curve 
! Subprogram not used 
! Subprogram not used    integer (int_kind), allocatable, dimension(:) :: work_per_proc,work_per_block2
! Subprogram not used    integer (int_kind) :: ierr  ! error return code 
! Subprogram not used 
! Subprogram not used    character(len=char_len) :: partitioning_type 
! Subprogram not used    logical, parameter :: verbose = .FALSE.
! Subprogram not used !------------------------------------------------------
! Subprogram not used ! Space filling curves only work if:
! Subprogram not used !
! Subprogram not used !    nblocks_x = nblocks_y
! Subprogram not used !       nblocks_x = 2^m 3^n 5^p where m,n,p are integers
! Subprogram not used !------------------------------------------------------
! Subprogram not used    
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1: nblocks_x,nblocks_y: ',nblocks_x,nblocks_y
! Subprogram not used    if((.not. IsFactorable(nblocks_y)) .or. (.not. IsFactorable(nblocks_x))) then
! Subprogram not used      if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.1'
! Subprogram not used      create_distrb_spacecurve = create_distrb_cart(nprocs, work_per_block)
! Subprogram not used      if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.2'
! Subprogram not used      return
! Subprogram not used    endif
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.3'
! Subprogram not used 
! Subprogram not used    !-----------------------------------------------
! Subprogram not used    ! Factor the numbers of blocks in each dimension
! Subprogram not used    !-----------------------------------------------
! Subprogram not used    xdim = Factor(nblocks_x)
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.4'
! Subprogram not used    ydim = Factor(nblocks_y)
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.5'
! Subprogram not used    numfac = xdim%numfact
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #2'
! Subprogram not used 
! Subprogram not used    !---------------------------------------------
! Subprogram not used    ! Match the common factors to create SFC curve
! Subprogram not used    !---------------------------------------------
! Subprogram not used    curveSize=1
! Subprogram not used    do it=1,numfac
! Subprogram not used       call MatchFactor(xdim,ydim,itmp,foundX)
! Subprogram not used       curveSize = itmp*curveSize
! Subprogram not used    enddo
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #3'
! Subprogram not used 
! Subprogram not used    !--------------------------------------
! Subprogram not used    ! determine the size of the sub-blocks 
! Subprogram not used    ! within the space-filling curve 
! Subprogram not used    !--------------------------------------
! Subprogram not used    sb_x = ProdFactor(xdim)
! Subprogram not used    sb_y = ProdFactor(ydim)
! Subprogram not used 
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #4'
! Subprogram not used    call create_communicator(dist%communicator, nprocs)
! Subprogram not used 
! Subprogram not used    dist%nprocs = nprocs
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  allocate space for decomposition
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #5'
! Subprogram not used    allocate(blockLocation(nblocks_tot),type_on_curve(nblocks_tot))
! Subprogram not used    allocate(distance(nblocks_tot))
! Subprogram not used    allocate(work_on_curve(nblocks_tot))
! Subprogram not used    allocate (dist%blockLocation(nblocks_tot), &
! Subprogram not used              dist%blockLocalID (nblocks_tot))
! Subprogram not used 
! Subprogram not used    allocate (dist%blockCnt(nprocs))
! Subprogram not used 
! Subprogram not used    blockLocation = 0
! Subprogram not used    dist%blockLocation=0
! Subprogram not used    dist%blockLocalID =0
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #6'
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  Create the array to hold the SFC and indices into it
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used    allocate(Mesh(curveSize,curveSize))
! Subprogram not used    allocate(Mesh2(nblocks_x,nblocks_y))
! Subprogram not used    allocate(idxT_i(nblocks_tot),idxT_j(nblocks_tot),Lindx(nblocks_tot))
! Subprogram not used    allocate(Lindx2(nblocks_tot))
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #7'
! Subprogram not used 
! Subprogram not used    Mesh  = 0
! Subprogram not used    Mesh2 = 0
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  Cost function estimation... this is a potential replace for work_per_block 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used    allocate(Cost_per_proc(nblocks_tot))
! Subprogram not used    allocate(Cost_per_block(nblocks_tot))
! Subprogram not used    allocate(work_per_proc(nprocs))
! Subprogram not used    allocate(work_per_block2(nblocks_tot))
! Subprogram not used    call EstimateCost(bStats,nblocks_tot,Cost_per_block)
! Subprogram not used    blockLocation=0
! Subprogram not used 
! Subprogram not used    do i=1,nblocks_tot
! Subprogram not used       work_per_block2(i) = NINT(10.0*ABS(Cost_per_block(i)),kind=int_kind)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #8'
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  Generate the space-filling curve
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used    call GenSpaceCurve(Mesh)
! Subprogram not used    Mesh = Mesh + 1    ! make it 1-based indexing
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #8.1'
! Subprogram not used    if(Debug) then
! Subprogram not used      if(my_task ==0) call PrintCurve(Mesh)
! Subprogram not used    endif
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #8.2'
! Subprogram not used    !-----------------------------------------------
! Subprogram not used    ! Reindex the SFC to address internal sub-blocks  
! Subprogram not used    !-----------------------------------------------
! Subprogram not used    do j=1,curveSize
! Subprogram not used    do i=1,curveSize
! Subprogram not used       sfcNum = (Mesh(i,j) - 1)*(sb_x*sb_y) + 1
! Subprogram not used       do jj=1,sb_y
! Subprogram not used       do ii=1,sb_x
! Subprogram not used          subNum = (jj-1)*sb_x + (ii-1)
! Subprogram not used          i2 = (i-1)*sb_x + ii
! Subprogram not used          j2 = (j-1)*sb_y + jj
! Subprogram not used          Mesh2(i2,j2) = sfcNum + subNum
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used    enddo
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #9'
! Subprogram not used    !------------------------------------------------
! Subprogram not used    ! create a linear array of i,j coordinates of SFC
! Subprogram not used    !------------------------------------------------
! Subprogram not used    idxT_i=0;idxT_j=0;Lindx=0;Lindx2=0
! Subprogram not used    do j=1,nblocks_y
! Subprogram not used      do i=1,nblocks_x
! Subprogram not used         n = (j-1)*nblocks_x + i
! Subprogram not used         ig = Mesh2(i,j)
! Subprogram not used         if(work_per_block(n) /= 0) then
! Subprogram not used             idxT_i(ig)=i;idxT_j(ig)=j
! Subprogram not used         endif
! Subprogram not used         Lindx(n) = ig
! Subprogram not used      enddo
! Subprogram not used    enddo
! Subprogram not used    do i=1,nblocks_tot
! Subprogram not used       type_on_curve(Lindx(i)) = blockType(i)
! Subprogram not used    enddo
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #10'
! Subprogram not used 
! Subprogram not used    ! ------------------------------
! Subprogram not used    ! compress out the land blocks
! Subprogram not used    ! ------------------------------
! Subprogram not used    ii=0
! Subprogram not used    do i=1,nblocks_tot
! Subprogram not used       if(IdxT_i(i) .gt. 0) then
! Subprogram not used          ii=ii+1
! Subprogram not used !         Mesh3(idxT_i(i),idxT_j(i)) = ii
! Subprogram not used          n = (idxT_j(i)-1)*nblocks_x + idxT_i(i)
! Subprogram not used          Lindx2(n) = ii 
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used    nblocks=ii
! Subprogram not used !DBG   if(my_task == 0) then 
! Subprogram not used !DBG     write(nu_diag,*) 'work_per_block2:',work_per_block2
! Subprogram not used !DBG   endif
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #11'
! Subprogram not used    allocate(cStats(numCoeff,nblocks)) 
! Subprogram not used    do i=1,nblocks_tot
! Subprogram not used       if(Lindx2(i)>0) then 
! Subprogram not used          work_on_curve(Lindx2(i)) = work_per_block2(i) 
! Subprogram not used          type_on_curve(Lindx2(i))  = blockType(i)
! Subprogram not used          cStats(:,Lindx2(i)) = bStats(:,i)
! Subprogram not used          distance(Lindx2(i)) = Lindx(i)
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #12'
! Subprogram not used 
! Subprogram not used    if(lprint_stats) then 
! Subprogram not used       fname = 'WorkPerBlock.bin'
! Subprogram not used       call WriteIntegerArray(fname,nblocks,work_on_curve) 
! Subprogram not used !DBG     write(nu_diag,*) 'work_per_block2:',work_per_block2
! Subprogram not used !DBG     write(nu_diag,*) 'work_on_curve:', work_on_curve(1:nblocks)
! Subprogram not used !     open(nu_timing,file='WorkPerBlock.bin',recl=4*nblocks, &
! Subprogram not used !          form = 'unformatted', access='direct',status='unknown')
! Subprogram not used !     write(nu_timing,rec=1) work_on_curve(1:nblocks)
! Subprogram not used !     close(nu_timing)
! Subprogram not used    endif 
! Subprogram not used 
! Subprogram not used    if(lprint_stats) then 
! Subprogram not used       call WriteProbabilityStats(bStats,nblocks_tot)
! Subprogram not used    endif
! Subprogram not used      
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #13'
! Subprogram not used 
! Subprogram not used    maxB=MIN(max_blocks,maxBlock)
! Subprogram not used    partitioning_type = 'weight'
! Subprogram not used    select case(partitioning_type)
! Subprogram not used       case ('type') 
! Subprogram not used          ! KLUDGE this is just for testing need to come up with a general solution
! Subprogram not used          numIce     = COUNT(blockType .eq. iceType)
! Subprogram not used          minblocks  = CEILING(REAL(numIce)/REAL(nprocs),kind=int_kind)
! Subprogram not used          !   write(nu_diag,*) 'before TypePartition: {minblocks,maxblocks}: ',minblocks,maxB 
! Subprogram not used          call  TypePartition(type_on_curve,minblocks,maxB,blockLocation)
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #14'
! Subprogram not used          if(MAXVAL(blockLocation) > nprocs) then 
! Subprogram not used               write(nu_diag,*) 'ERROR: problem with partitioning: insufficient processors'
! Subprogram not used          endif
! Subprogram not used          ! re-index blockLocation from curve to physical ordering
! Subprogram not used          do i=1,nblocks_tot
! Subprogram not used              dist%blockLocation(i) = blockLocation(Lindx(i))
! Subprogram not used          enddo
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #15'
! Subprogram not used       case ('weight')
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #16'
! Subprogram not used          call PartitionCurve(work_on_curve(1:nblocks),work_per_proc, &
! Subprogram not used 		blockLocation(1:nblocks),distance(1:nblocks), nprocs,minBlock, maxB,cStats,FixMaxBlock, maxDil, ierr)
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #17'
! Subprogram not used !DBG         write(nu_diag,*) 'After PartitionCurve:'
! Subprogram not used          if(ierr < 0) then 
! Subprogram not used              call abort_ice('create_distrb_spacecurve: PartitionCurve failed')
! Subprogram not used          endif 
! Subprogram not used !DBG         write(nu_diag,*) 'before broadcast_array:'
! Subprogram not used          call broadcast_array(blockLocation,master_task)
! Subprogram not used !DBG         write(nu_diag,*) 'After broadcast_array'
! Subprogram not used          ! re-index blockLocation from curve to physical ordering
! Subprogram not used          numLocalBlocks=0
! Subprogram not used          do i=1,nblocks_tot
! Subprogram not used              if(Lindx2(i)>0) then  
! Subprogram not used                 dist%blockLocation(i) = blockLocation(Lindx2(i))
! Subprogram not used                 numLocalBlocks=numLocalBlocks+1
! Subprogram not used              endif
! Subprogram not used          enddo
! Subprogram not used     end select
! Subprogram not used !   call qsort(dist%blockLocation(1:numLocalBlocks))
! Subprogram not used !   write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: dist%blockLocation(:)', dist%blockLocation(1:nblocks_tot)
! Subprogram not used !   write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: dist%blockLocation(1:numLocalBlocks)', &
! Subprogram not used !        dist%blockLocation(1:numLocalBlocks)
! Subprogram not used   
! Subprogram not used !   call ConvertStatsBlock2Proc(dist%blockLocation,bStats,cStats)
! Subprogram not used !DBG   write(nu_diag,*) 'Before call to BuildProbabilityStats2'
! Subprogram not used    call BuildProbabilityStats2(dist%blockLocation,cStats)
! Subprogram not used !DBG   write(nu_diag,*) 'After call to BuildProbabilityStats2'
! Subprogram not used 
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #19'
! Subprogram not used    call EstimateCost(cStats,nprocs,Cost_per_proc)
! Subprogram not used !DBG   write(nu_diag,*) 'before WriteProbabilityStats'
! Subprogram not used !DBG   write(nu_diag,*) 'after WriteProbabilityStats'
! Subprogram not used 
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #20'
! Subprogram not used    if(lprint_stats) then 
! Subprogram not used       fname = 'Q.bin'
! Subprogram not used       call WriteIntegerArray(fname,nblocks_tot,dist%blockLocation) 
! Subprogram not used 
! Subprogram not used       fname = 'Cost.bin'
! Subprogram not used       call WriteDblArray(fname,nprocs,Cost_per_proc)
! Subprogram not used 
! Subprogram not used       fname = 'perfm.bin'
! Subprogram not used       call WriteDblArray(fname,numCoeff,perfmodel)
! Subprogram not used    endif
! Subprogram not used    
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #21'
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  Reset the dist data structure
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used    allocate(proc_tmp(nprocs))
! Subprogram not used    proc_tmp = 0
! Subprogram not used 
! Subprogram not used    allocate(dist%blockIndex(nprocs,max_blocks))
! Subprogram not used    dist%blockIndex(:,:) = 0
! Subprogram not used    do n=1,nblocks_tot
! Subprogram not used       pid = dist%blockLocation(n)
! Subprogram not used       if(pid>0) then
! Subprogram not used         proc_tmp(pid) = proc_tmp(pid) + 1
! Subprogram not used         dist%blockLocalID(n) = proc_tmp(pid)
! Subprogram not used         if(proc_tmp(pid) <= max_blocks) then 
! Subprogram not used             dist%blockIndex(pid,proc_tmp(pid)) = n
! Subprogram not used         endif
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used    dist%blockCnt(:) = proc_tmp(:) 
! Subprogram not used 
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #22'
! Subprogram not used    !---------------------------------------
! Subprogram not used    ! Set the number of active local blocks
! Subprogram not used    !---------------------------------------
! Subprogram not used    dist%numLocalBlocks = proc_tmp(my_task+1)
! Subprogram not used    if (dist%numLocalBlocks>0) then 
! Subprogram not used       allocate(dist%blockGlobalID(dist%numLocalBlocks))
! Subprogram not used       localID=1
! Subprogram not used       do n=1,nblocks_tot
! Subprogram not used          pid = dist%blockLocation(n)
! Subprogram not used          if(pid == my_task+1) then 
! Subprogram not used               dist%blockGlobalID(localID) = n
! Subprogram not used               localID=localID+1
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used    endif
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #23'
! Subprogram not used    if(Debug) then
! Subprogram not used       if(my_task==0) write(nu_diag,*) 'dist%blockLocation:= ',dist%blockLocation
! Subprogram not used       write(nu_diag,*) 'IAM: ',my_task,' SpaceCurve: Number of blocks {total,local} :=', &
! Subprogram not used                 nblocks_tot,nblocks,proc_tmp(my_task+1)
! Subprogram not used    endif
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate block'  
! Subprogram not used    !---------------------------------
! Subprogram not used    ! Deallocate temporary arrays
! Subprogram not used    !---------------------------------
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(blockLocation)'  
! Subprogram not used    deallocate(blockLocation)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(type_on_curve)'  
! Subprogram not used    deallocate(type_on_curve)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(work_on_curve)'  
! Subprogram not used    deallocate(work_on_curve)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(proc_tmp)'  
! Subprogram not used    deallocate(proc_tmp)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(Mesh,Mesh2)'  
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #24'
! Subprogram not used    deallocate(Mesh,Mesh2)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(idxT_i,idxT_j,Lindx,Lindx2)'  
! Subprogram not used    deallocate(idxT_i,idxT_j,Lindx,Lindx2)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(Cost_per_proc)'  
! Subprogram not used    deallocate(Cost_per_proc) 
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(Cost_per_block)'  
! Subprogram not used    deallocate(Cost_per_block)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(work_per_proc)'  
! Subprogram not used    deallocate(work_per_proc)
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(work_per_block2)'  
! Subprogram not used    deallocate(work_per_block2)  
! Subprogram not used     
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #25'
! Subprogram not used    if(verbose .and. my_task == 0) then 
! Subprogram not used       write(nu_diag,*) 'create_distrb_spacecurve: blockCnt ',dist%blockCnt
! Subprogram not used       write(nu_diag,*) 'create_distrb_spacecurve: blockIndex ',dist%blockIndex
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !DBG   write(nu_diag,*) 'create_distrb_spacecurve: before assignment of result'  
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used    create_distrb_spacecurve = dist  ! return the result
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (verbose .and. my_task==0) then 
! Subprogram not used       write(nu_diag,*) 'create_distrb_spacecurve%blockGlobalID:  ',create_distrb_spacecurve%blockGlobalID
! Subprogram not used       write(nu_diag,*) 'create_distrb_spacecurve%blockCnt: ',create_distrb_spacecurve%blockCnt
! Subprogram not used    endif
! Subprogram not used !DBG   write(nu_diag,*) 'At the end of create_distrb_spacecurve'  
! Subprogram not used    if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #26'
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function create_distrb_spacecurve

! Subprogram not used  subroutine TypePartition(blockType,minblocks,maxblocks,blockLocation)
! Subprogram not used 
! Subprogram not used    integer(kind=int_kind), intent(in)     :: blockType(nblocks_tot)
! Subprogram not used    integer, intent(in)  :: minblocks, & ! Minimum number of blocks per processor 
! Subprogram not used                            maxblocks    ! Maximum number of blocks per processor 
! Subprogram not used    integer(kind=int_kind), intent(inout)  :: blockLocation(nblocks_tot)
! Subprogram not used 
! Subprogram not used    integer :: i,ip,cur,next 
! Subprogram not used    integer :: cntLnd,tcntLnd, &
! Subprogram not used               cntIcefree,tcntIcefree, &
! Subprogram not used               cntIce,tcntIce
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    tcntIce     = 0
! Subprogram not used    tcntLnd     = 0
! Subprogram not used    tcntIcefree = 0
! Subprogram not used    cntIce      = 0
! Subprogram not used    cntLnd      = 0
! Subprogram not used    cntIcefree  = 0
! Subprogram not used    ip = 1
! Subprogram not used    do i=1,nblocks_tot
! Subprogram not used       cur  = blockType(i)
! Subprogram not used       if(i<nblocks_tot) then 
! Subprogram not used          next = blockType(i+1)
! Subprogram not used       else
! Subprogram not used          next = cur
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !------------ 
! Subprogram not used       ! Land point 
! Subprogram not used       !------------ 
! Subprogram not used       if(cur == lndType) then 
! Subprogram not used           tcntLnd = tcntLnd + 1
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used      !----------------------
! Subprogram not used      ! ice Free  point
! Subprogram not used      !----------------------
! Subprogram not used      if (cur == icefreeType) then 
! Subprogram not used         tcntIcefree = tcntIcefree + 1
! Subprogram not used         cntIcefree = cntIcefree + 1
! Subprogram not used         blockLocation(i) = ip
! Subprogram not used         if ((cntIcefree == maxblocks) .or. (next /= cur)) then
! Subprogram not used            ip = ip + 1
! Subprogram not used            cntIcefree=0
! Subprogram not used         endif
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      !----------------------
! Subprogram not used      ! ice Free  point
! Subprogram not used      !----------------------
! Subprogram not used      if (cur == iceType) then 
! Subprogram not used         tcntIce = tcntIce + 1
! Subprogram not used         cntIce  = cntIce + 1
! Subprogram not used         blockLocation(i) = ip
! Subprogram not used         if ((cntIce == minblocks) .or. (next /= cur)) then
! Subprogram not used            ip = ip + 1
! Subprogram not used            cntIce=0
! Subprogram not used         endif
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    if(my_task == 0) then 
! Subprogram not used       write(*,23) tcntLnd+tcntIce+tcntIcefree, tcntIce, tcntIcefree, tcntLnd
! Subprogram not used !      write(nu_diag,*) 'TypePartition: land blks: ',tcntLnd,' Ice blks: ', &
! Subprogram not used !                tcntIce,' IceFree blks: ',tcntIcefree
! Subprogram not used       write(*,24) MAXVAL(blockLocation) 
! Subprogram not used !      write(nu_diag,*) 'TypePartition: Partitioned across ',MAXVAL(blockLocation),' processors'
! Subprogram not used    endif
! Subprogram not used   
! Subprogram not used 23   format('Total blocks: ',i5,' Ice blocks: ',i5,' IceFree blocks: ',i5,' Land blocks: ',i5)
! Subprogram not used 24   format('Partitioned across ',i5,' processors')
! Subprogram not used   
! Subprogram not used 
! Subprogram not used  end subroutine TypePartition

! Subprogram not used   subroutine PartitionCurve(work_per_block, work_per_proc, blockLocation, distance, &
! Subprogram not used              nproc, min_blocks, max_blocks, Stats, FixMaxBlock, maxDil, ierr)
! Subprogram not used 
! Subprogram not used     integer (int_kind), intent(inout) :: work_per_block(:)
! Subprogram not used     integer (int_kind), intent(inout) :: work_per_proc(:)
! Subprogram not used     integer (int_kind), intent(inout) :: blockLocation(:)
! Subprogram not used     integer (int_kind), intent(inout) :: distance(:)
! Subprogram not used     integer (int_kind), intent(in)    :: nproc
! Subprogram not used     integer (int_kind), intent(in)    :: min_blocks
! Subprogram not used     integer (int_kind), intent(in)    :: max_blocks
! Subprogram not used     real    (dbl_kind), intent(in)    :: Stats(:,:) 
! Subprogram not used     logical,            intent(in)    :: FixMaxBlock
! Subprogram not used     real    (real_kind), intent(in)   :: maxDil
! Subprogram not used     integer (int_kind), intent(inout) :: ierr
! Subprogram not used 
! Subprogram not used     integer :: cnt
! Subprogram not used 
! Subprogram not used     integer :: nb,anProc
! Subprogram not used     integer :: n,ip,i,imax
! Subprogram not used     real :: totalCost, avgCost, maxCost,dtCost
! Subprogram not used     integer, allocatable :: saveblockLocation(:)
! Subprogram not used     integer :: maxBlocks,minBlocks
! Subprogram not used     integer :: ivalue
! Subprogram not used     real (real_kind) :: maxValue
! Subprogram not used     real :: maxValue_save
! Subprogram not used     real :: minValue
! Subprogram not used     logical :: contLoop
! Subprogram not used 
! Subprogram not used     integer :: maxB,minB
! Subprogram not used 
! Subprogram not used     logical, parameter :: verbose = .TRUE.
! Subprogram not used     integer :: it,maxiter
! Subprogram not used     integer :: save_maxB 
! Subprogram not used     real :: save_maxCost
! Subprogram not used     
! Subprogram not used     integer :: aminBlocks, amaxBlocks 
! Subprogram not used 
! Subprogram not used     real :: aWork,minWork,maxWork,maxWorkBlock
! Subprogram not used     real :: minCostBlock, maxCostBlock
! Subprogram not used     real :: maxCost_old
! Subprogram not used     real (real_kind) :: amaxDil
! Subprogram not used     real (dbl_kind), allocatable, dimension(:) :: cost_per_proc
! Subprogram not used     real (dbl_kind), allocatable, dimension(:) :: cost_per_block
! Subprogram not used     real (dbl_kind), allocatable, dimension(:,:) :: pStats 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     allocate(pStats(numCoeff,nblocks_tot))
! Subprogram not used 
! Subprogram not used     maxB=max_blocks
! Subprogram not used !    maxDil = 5.0
! Subprogram not used 
! Subprogram not used     if ( my_task .eq. master_task) then  
! Subprogram not used 
! Subprogram not used        nb = SIZE(blockLocation)  ! number of blocks
! Subprogram not used        allocate(saveblockLocation(nb))
! Subprogram not used        allocate(cost_per_proc(nb))
! Subprogram not used        allocate(cost_per_block(nb))
! Subprogram not used 
! Subprogram not used        ! --------------------------------------------
! Subprogram not used        ! Estimate the computational cost of each block  
! Subprogram not used        ! --------------------------------------------
! Subprogram not used        call EstimateCost(Stats,nb,cost_per_block)
! Subprogram not used 
! Subprogram not used        totalCost = SUM(cost_per_block)
! Subprogram not used        save_maxCost = totalCost
! Subprogram not used        maxCostBlock = MAXVAL(cost_per_block)
! Subprogram not used        minCostBlock = MINVAL(cost_per_block)
! Subprogram not used        save_maxB = maxB
! Subprogram not used        avgCost = (totalCost/nproc)
! Subprogram not used 
! Subprogram not used        write(nu_diag,*) 'PartitionCurve: nblocks,nproc ',nb,nproc
! Subprogram not used        write(nu_diag,213) totalCost, avgCost, minCostBlock, maxCostBlock
! Subprogram not used !DBG       write(nu_diag,*) 'PartitionCurve: totalCost,avgCost, maxCostBlock: ',totalCost,avgCost,maxCostBlock
! Subprogram not used !DBG       write(nu_diag,*) distance
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        minB = CEILING(real(nb)/real(nproc),kind=int_kind)
! Subprogram not used        if(maxB < minB ) then
! Subprogram not used           write(nu_diag,*) 'ERROR: unable to partition ',nb,' blocks across ',nproc,' processors'
! Subprogram not used           write(nu_diag,*) 'ERROR: Either increase max_blocks := ',maxB
! Subprogram not used           write(nu_diag,*) 'ERROR: Either increase number of processors'
! Subprogram not used           ierr = -2
! Subprogram not used           return
! Subprogram not used        endif
! Subprogram not used        if(FixMaxBlock) then 
! Subprogram not used          minB = maxB
! Subprogram not used        endif
! Subprogram not used        do while (maxB  >= minB )  
! Subprogram not used 
! Subprogram not used           dtCost = 1.0
! Subprogram not used           maxValue = maxCostBlock*real(maxB)
! Subprogram not used           minValue = maxCostBlock
! Subprogram not used           maxCost_old = maxValue
! Subprogram not used 
! Subprogram not used           contLoop = .true.
! Subprogram not used           maxiter = 20
! Subprogram not used           it = 1
! Subprogram not used           do while(contLoop )
! Subprogram not used 
! Subprogram not used             cost_per_proc=0.0
! Subprogram not used             call wPartition(cost_per_block,blockLocation, distance, nproc,min_blocks, maxB,maxValue,maxDil,aminBlocks, &
! Subprogram not used                  amaxBlocks, amaxDil)
! Subprogram not used             anProc = MAXVAL(blockLocation)
! Subprogram not used             call ConvertStatsBlock2Proc(blockLocation,Stats,pStats)
! Subprogram not used             call EstimateCost(pStats,anProc,cost_per_proc)
! Subprogram not used             maxCost = MAXVAL(cost_per_proc)
! Subprogram not used    
! Subprogram not used 
! Subprogram not used             if(lprint_stats) then 
! Subprogram not used 		write(nu_diag,211) it,anProc,aminBlocks,amaxBlocks,maxB,minValue, maxValue,maxCost, amaxDil
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             if(maxCost > maxValue) then
! Subprogram not used                minValue =  maxValue
! Subprogram not used                dtCost = (maxCost_old-minValue)/2.0
! Subprogram not used                maxValue = maxCost_old - dtCost
! Subprogram not used             else
! Subprogram not used                dtCost = (maxCost-minValue)/2.0
! Subprogram not used                maxValue = maxCost - dtCost
! Subprogram not used                maxCost_old = maxCost
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             if(maxCost == maxCostBlock) contLoop = .false.
! Subprogram not used             if(dtCost < 1.0e-5) contLoop = .false.
! Subprogram not used             if( anProc == nproc .and. it >= maxiter)  contLoop = .false.
! Subprogram not used             if ((save_maxCost > maxCost .and. amaxBlocks <= maxB) .or. &
! Subprogram not used                 (save_maxCost == maxCost .and. save_maxB > maxB)  ) then
! Subprogram not used                 save_maxCost = maxCost
! Subprogram not used                 save_maxB = maxB  
! Subprogram not used                 saveblockLocation = blockLocation
! Subprogram not used             endif
! Subprogram not used             it=it+1
! Subprogram not used           enddo
! Subprogram not used           maxB = maxB - 1
! Subprogram not used        enddo
! Subprogram not used        blockLocation = saveblockLocation
! Subprogram not used        write(nu_diag,214) perfmodel_name
! Subprogram not used 
! Subprogram not used        write(nu_diag,*) '-------------------------wSFC-----------------------'
! Subprogram not used        call PrintPartitionLB(blockLocation,nproc,Stats) 
! Subprogram not used 
! Subprogram not used  211 format('Partition loop: it: ',i4,' anProc: ',i4,' a{min,max}Blocks ', &
! Subprogram not used               2(i4),' maxB: ',i4,' [min,max]Value: ',2(f14.4),' maxCost: ',f14.4,' maxDilation: ',f14.4 )  
! Subprogram not used  213 format('PartitionCurve: TotalCost: ',f12.4,' avgCost: ',f12.4, &
! Subprogram not used 	    ' minBlockCost: ',f12.4,' maxBlockCost: ',f12.4 ) 
! Subprogram not used  214 format('Using performance model: ',a20)
! Subprogram not used 
! Subprogram not used        deallocate(cost_per_proc)
! Subprogram not used        deallocate(cost_per_block)
! Subprogram not used        deallocate(saveblockLocation)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     deallocate(pStats)
! Subprogram not used 
! Subprogram not used     ierr = 0
! Subprogram not used 
! Subprogram not used   end subroutine PartitionCurve

! Subprogram not used   subroutine wPartition(cost_per_block, blockLocation, distance, nproc, min_blocks, max_blocks, maxValue, maxDil, aminBlocks, &
! Subprogram not used        amaxBlocks,amaxDil)
! Subprogram not used 
! Subprogram not used     real (dbl_kind), intent(in) :: cost_per_block(:)
! Subprogram not used     integer (int_kind), intent(inout) :: blockLocation(:)
! Subprogram not used     integer (int_kind), intent(inout) :: distance(:)
! Subprogram not used     integer (int_kind), intent(in) :: nproc
! Subprogram not used     integer (int_kind), intent(in) :: min_blocks
! Subprogram not used     integer (int_kind), intent(in) :: max_blocks
! Subprogram not used     real (real_kind), intent(in) :: maxvalue
! Subprogram not used     real (real_kind), intent(in)  :: maxDil    ! maximum alloable dilation of domains
! Subprogram not used     integer (int_kind), intent(inout) :: aminBlocks, amaxBlocks 
! Subprogram not used     real (real_kind), intent(out) :: amaxDil
! Subprogram not used 
! Subprogram not used     integer (int_kind)  :: n,ip,i,numB
! Subprogram not used     real :: totalCost,avgCost,sumTMP,sumTMP2
! Subprogram not used     
! Subprogram not used !   integer (int_kind), allocatable :: work_per_proc(:)
! Subprogram not used     
! Subprogram not used     integer :: minDist,Dist
! Subprogram not used     integer :: sloc
! Subprogram not used     real    :: dilation, dilation2
! Subprogram not used 
! Subprogram not used     logical, parameter :: Info = .FALSE.
! Subprogram not used     logical  :: break_loop
! Subprogram not used 
! Subprogram not used     n = SIZE(blockLocation)  ! number of blocks
! Subprogram not used 
! Subprogram not used     totalCost = SUM(cost_per_block)
! Subprogram not used 
! Subprogram not used     aminBlocks = n
! Subprogram not used     amaxBlocks = 0
! Subprogram not used     amaxDil  = 1.0
! Subprogram not used     avgCost = (totalCost/nproc)
! Subprogram not used 
! Subprogram not used !DBG    write(nu_diag,*) 'cost_per_block: ',cost_per_block
! Subprogram not used 
! Subprogram not used     ip = 1
! Subprogram not used     i=1
! Subprogram not used !    cost_per_block=0
! Subprogram not used     do while (i<=n)
! Subprogram not used         sumTMP = 0
! Subprogram not used         numB = 0
! Subprogram not used         break_loop = .FALSE.
! Subprogram not used         sloc = distance(i)
! Subprogram not used         do while ( ((sumTMP<maxValue) .or. (ip ==  nproc))  &
! Subprogram not used 		.and. (i<=n) &
! Subprogram not used 		.and. (.not. break_loop))
! Subprogram not used           sumTMP2  = sumTMP + cost_per_block(i)
! Subprogram not used           minDist  = numB + 1
! Subprogram not used           Dist     = distance(i) - sloc + 1
! Subprogram not used           dilation2 = real(Dist)/real(minDist)
! Subprogram not used           if(((sumTMP2 <= maxValue) .and. (numB < max_blocks) .and. dilation2 <= maxDil) &
! Subprogram not used !          if(((sumTMP2 <= maxValue) .and. (numB < max_blocks)) &
! Subprogram not used                  .or. (ip == nproc) .or. (numB < min_blocks) ) then
! Subprogram not used               blockLocation(i) = ip
! Subprogram not used               i=i+1
! Subprogram not used               sumTMP = sumTMP2
! Subprogram not used               dilation = dilation2
! Subprogram not used               numB = numB+1
! Subprogram not used           else
! Subprogram not used               break_loop = .TRUE.
! Subprogram not used           endif
! Subprogram not used         enddo
! Subprogram not used         if(aminBlocks > numB)  aminBlocks = numB
! Subprogram not used         if(amaxBlocks < numB)  amaxBlocks = numB
! Subprogram not used         if(amaxDil < dilation) amaxDil = dilation
! Subprogram not used         ip = ip+1
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used !    call abort_ice('wPartition: at the end of the subroutine')
! Subprogram not used 
! Subprogram not used   end subroutine wPartition

!**********************************************************************
!BOP
! !IROUTINE: ice_distributionRake
! !INTERFACE:

! Subprogram not used  subroutine ice_distributionRake (procWork, procID, blockWork, &
! Subprogram not used                                   priority, distribution)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This subroutine performs a rake algorithm to distribute the work
! Subprogram not used !  along a vector of processors.  In the rake algorithm, a work
! Subprogram not used !  threshold is first set.  Then, moving from left to right, work
! Subprogram not used !  above that threshold is raked to the next processor in line.
! Subprogram not used !  The process continues until the end of the vector is reached
! Subprogram not used !  and then the threshold is reduced by one for a second rake pass.
! Subprogram not used !  In this implementation, a priority for moving blocks is defined
! Subprogram not used !  such that the rake algorithm chooses the highest priority
! Subprogram not used !  block to be moved to the next processor.  This can be used
! Subprogram not used !  for example to always choose the eastern-most block or to
! Subprogram not used !  ensure a block does not stray too far from its neighbors.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in), dimension(:) :: &
! Subprogram not used       blockWork          ,&! amount of work per block
! Subprogram not used       procID               ! global processor number
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(inout), dimension(:) :: &
! Subprogram not used       procWork           ,&! amount of work per processor
! Subprogram not used       priority             ! priority for moving a given block
! Subprogram not used 
! Subprogram not used    type (distrb), intent(inout) :: &
! Subprogram not used       distribution         ! distribution to change
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i, n,                  &! dummy loop indices
! Subprogram not used       np1,                   &! n+1 corrected for cyclical wrap
! Subprogram not used       iproc, inext,          &! processor ids for current and next 
! Subprogram not used       nprocs, numBlocks,   &! number of blocks, processors
! Subprogram not used       lastPriority,          &! priority for most recent block
! Subprogram not used       minPriority,           &! minimum priority
! Subprogram not used       lastLoc,               &! location for most recent block
! Subprogram not used       meanWork, maxWork,     &! mean,max work per processor
! Subprogram not used       diffWork, residual,    &! work differences and residual work
! Subprogram not used       numTransfers            ! counter for number of block transfers
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  initialization
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    nprocs  = size(procWork)
! Subprogram not used    numBlocks = size(blockWork)
! Subprogram not used 
! Subprogram not used    !*** compute mean,max work per processor
! Subprogram not used 
! Subprogram not used    meanWork = sum(procWork)/nprocs + 1
! Subprogram not used    maxWork  = maxval(procWork)
! Subprogram not used    residual = mod(meanWork,nprocs)
! Subprogram not used 
! Subprogram not used    minPriority = 1000000
! Subprogram not used    do n=1,nprocs
! Subprogram not used       iproc = procID(n)
! Subprogram not used       do i=1,numBlocks
! Subprogram not used          if (distribution%blockLocation(i) == iproc) then
! Subprogram not used             minPriority = min(minPriority,priority(i))
! Subprogram not used          endif
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  do two sets of transfers
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    transferLoop: do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !     do rake across the processors
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       numTransfers = 0
! Subprogram not used       do n=1,nprocs
! Subprogram not used          if (n < nprocs) then
! Subprogram not used             np1   = n+1
! Subprogram not used          else
! Subprogram not used             np1   = 1
! Subprogram not used          endif
! Subprogram not used          iproc = procID(n)
! Subprogram not used          inext = procID(np1)
! Subprogram not used 
! Subprogram not used          if (procWork(n) > meanWork) then !*** pass work to next
! Subprogram not used 
! Subprogram not used             diffWork = procWork(n) - meanWork
! Subprogram not used 
! Subprogram not used             rake1: do while (diffWork > 1)
! Subprogram not used 
! Subprogram not used                !*** attempt to find a block with the required
! Subprogram not used                !*** amount of work and with the highest priority
! Subprogram not used                !*** for moving (eg boundary blocks first)
! Subprogram not used 
! Subprogram not used                lastPriority = 0
! Subprogram not used                lastLoc = 0
! Subprogram not used 
! Subprogram not used                do i=1,numBlocks
! Subprogram not used                   if (distribution%blockLocation(i) == iproc) then
! Subprogram not used                      if (priority(i) > lastPriority ) then
! Subprogram not used                         lastPriority = priority(i)
! Subprogram not used                         lastLoc = i
! Subprogram not used                      endif
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used                if (lastLoc == 0) exit rake1 ! could not shift work
! Subprogram not used 
! Subprogram not used                numTransfers = numTransfers + 1
! Subprogram not used                distribution%blockLocation(lastLoc) = inext
! Subprogram not used                if (np1 == 1) priority(lastLoc) = minPriority
! Subprogram not used                diffWork = diffWork - blockWork(lastLoc)
! Subprogram not used 
! Subprogram not used                procWork(n  ) = procWork(n  )-blockWork(lastLoc)
! Subprogram not used                procWork(np1) = procWork(np1)+blockWork(lastLoc)
! Subprogram not used             end do rake1
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !     increment meanWork by one and repeat
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       meanWork = meanWork + 1
! Subprogram not used       if (numTransfers == 0 .or. meanWork > maxWork) exit transferLoop
! Subprogram not used 
! Subprogram not used    end do transferLoop
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used end subroutine ice_distributionRake

!***********************************************************************
! Subprogram not used     subroutine PrintPartitionLB(Location,n,bStats)
! Subprogram not used 
! Subprogram not used       integer :: Location(:)
! Subprogram not used       integer :: n
! Subprogram not used       real (dbl_kind), dimension(:,:) :: bStats 
! Subprogram not used 
! Subprogram not used       real :: maxCost,minCost
! Subprogram not used       integer :: anProc
! Subprogram not used       integer :: minB, maxB 
! Subprogram not used       real :: aCost
! Subprogram not used       real (dbl_kind), allocatable, dimension(:) :: cost_per_proc(:)
! Subprogram not used       integer :: i,ncnt
! Subprogram not used 
! Subprogram not used       real (dbl_kind), allocatable, dimension(:,:) :: pStats
! Subprogram not used 
! Subprogram not used       allocate(cost_per_proc(n))
! Subprogram not used       allocate(pStats(numCoeff,n))
! Subprogram not used       
! Subprogram not used        
! Subprogram not used       call ConvertStatsBlock2Proc(Location,bStats,pStats)
! Subprogram not used       call EstimateCost(pStats,n,cost_per_proc)
! Subprogram not used 
! Subprogram not used !DBG      write(nu_diag,*) 'PrintPartitinoLB: Location:',Location
! Subprogram not used !DBG      write(nu_diag,*) 'PrintPartitinoLB: cost_per_proc:',cost_per_proc
! Subprogram not used 
! Subprogram not used       maxCost = MAXVAL(cost_per_proc)
! Subprogram not used       minCost = MINVAL(cost_per_proc)
! Subprogram not used       aCost = SUM(cost_per_proc)/real(n)
! Subprogram not used       anProc=MAXVAL(Location)
! Subprogram not used       maxB=0
! Subprogram not used       minB=n
! Subprogram not used       do i=1,n
! Subprogram not used            ncnt = COUNT(Location==i)  
! Subprogram not used            maxB = MAX(ncnt,maxB)
! Subprogram not used            minB = MIN(ncnt,minB)
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used !DBG      write(nu_diag,*) 'maxCost: ',maxCost
! Subprogram not used !DBG      write(nu_diag,*) 'minCost: ',minCost
! Subprogram not used !DBG      write(nu_diag,*) 'aCost: ',aCost
! Subprogram not used 
! Subprogram not used !      write(nu_diag,*) 'PrintPartitionLB: on ',anProc,' processors Avg,Min,Max work/proc, imbalance  ', &
! Subprogram not used !                aWork,minWork,maxWork,ABS(aWork-maxWork)/aWork
! Subprogram not used       write(nu_diag,212)  anProc, minB, maxB,aCost,minCost,maxCost,ABS(aCost-maxCost)/aCost
! Subprogram not used       deallocate(cost_per_proc)
! Subprogram not used       deallocate(pStats)
! Subprogram not used 
! Subprogram not used  212 format('PrintPartitionLB: on ',i4,' procs a{min,max}Blocks: ',(2(i4)),' Avg,Min,Max Cost/proc: ',(3(f10.4)), &
! Subprogram not used           ' imbalance: ',f8.2)
! Subprogram not used 
! Subprogram not used     end subroutine PrintPartitionLB

! Subprogram not used    subroutine EstimateCost(coeffMatrix,n,Cost)
! Subprogram not used 
! Subprogram not used      real (dbl_kind) :: coeffMatrix(:,:)
! Subprogram not used      integer (int_kind) :: n
! Subprogram not used      real (dbl_kind):: Cost(:)
! Subprogram not used 
! Subprogram not used      real (dbl_kind) :: tmp
! Subprogram not used 
! Subprogram not used      integer (int_kind) :: i,j
! Subprogram not used 
! Subprogram not used      Cost=0.0_dbl_kind
! Subprogram not used      do i=1,n
! Subprogram not used         tmp = 0.0_dbl_kind
! Subprogram not used         do j=1,numCoeff
! Subprogram not used            tmp = tmp + coeffMatrix(j,i) *perfmodel(j)
! Subprogram not used         enddo
! Subprogram not used         Cost(i) = tmp
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used    end subroutine EstimateCost

! Subprogram not used    subroutine ConvertStatsBlock2Proc(Location,bStats,pStats)
! Subprogram not used 
! Subprogram not used       integer (int_kind) :: Location(:)
! Subprogram not used       real (dbl_kind), intent(in)  :: bStats(:,:)
! Subprogram not used       real (dbl_kind), intent(out) :: pStats(:,:)
! Subprogram not used 
! Subprogram not used       integer (int_kind) :: i,ip,n
! Subprogram not used 
! Subprogram not used       n = size(Location)
! Subprogram not used       pStats = 0.0d0
! Subprogram not used       do i=1,n
! Subprogram not used          ip = Location(i)
! Subprogram not used          if(ip > 0) then 
! Subprogram not used             pStats(:,ip) = pStats(:,ip) + bStats(:,i) 
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used    end subroutine ConvertStatsBlock2Proc

! Subprogram not used    subroutine WriteProbabilityStats(coeffMatrix,n)
! Subprogram not used 
! Subprogram not used     real(dbl_kind)  :: coeffMatrix(:,:)
! Subprogram not used     integer (int_kind) :: n
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if(my_task == master_task) then
! Subprogram not used        open(nu_timing,file='probStats.bin',recl=8*numCoeff*n, &
! Subprogram not used             form = 'unformatted', access = 'direct', status = 'unknown')
! Subprogram not used        write(nu_timing,rec=1) coeffMatrix(:,1:n)
! Subprogram not used        close(nu_timing)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    end subroutine WriteProbabilityStats

! Subprogram not used    subroutine WriteIntegerArray(fname,n,array)
! Subprogram not used      character(char_len) :: fname
! Subprogram not used      integer (int_kind) :: n
! Subprogram not used      integer (int_kind) :: array(:)
! Subprogram not used 
! Subprogram not used      if(my_task == master_task) then 
! Subprogram not used         open(nu_timing,file=TRIM(fname),recl=4*n, &
! Subprogram not used           form = 'unformatted', access = 'direct', status = 'unknown')
! Subprogram not used         write(nu_timing,rec=1) array(1:n)
! Subprogram not used         close(nu_timing)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used    end subroutine WriteIntegerArray

! Subprogram not used    subroutine WriteDblArray(fname,n,array)
! Subprogram not used      character(char_len) :: fname
! Subprogram not used      real (dbl_kind) :: array(:)
! Subprogram not used      integer (int_kind) :: n
! Subprogram not used 
! Subprogram not used      if(my_task == master_task) then 
! Subprogram not used         open(nu_timing,file=TRIM(fname),recl=8*n, &
! Subprogram not used           form = 'unformatted', access = 'direct', status = 'unknown')
! Subprogram not used         write(nu_timing,rec=1) array(1:n)
! Subprogram not used         close(nu_timing)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used    end subroutine WriteDblArray

end module ice_distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
