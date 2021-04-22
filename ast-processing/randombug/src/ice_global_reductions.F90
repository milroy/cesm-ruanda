!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!BOP
! !MODULE: ice_global_reductions

 module ice_global_reductions

! !DESCRIPTION:
!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_global_reductions.F90 112 2008-03-13 21:06:56Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Feb. 2008: Updated from POP version by Elizabeth C. Hunke, LANL
!
! !USES:

   use ice_kinds_mod
   use ice_communicate
   use ice_constants
   use ice_blocks
   use ice_distribution
   use ice_domain_size

   implicit none

   private
   include 'mpif.h'
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: global_sum,      &
             global_sum_prod, &
             global_maxval,   &
             global_minval,   &
             init_global_reductions

   public :: sum_vector_dbl
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface global_sum
     module procedure global_sum_dbl,              &
                      global_sum_real,             &
                      global_sum_int,              &
                      global_sum_scalar_dbl,       &
                      global_sum_scalar_real,      &
                      global_sum_scalar_int
   end interface

   interface global_sum_prod
     module procedure global_sum_prod_dbl,         &
                      global_sum_prod_real,        &
                      global_sum_prod_int
   end interface

   interface global_maxval
     module procedure global_maxval_dbl,           &
                      global_maxval_real,          &
                      global_maxval_int,           &
                      global_maxval_scalar_dbl,    &
                      global_maxval_scalar_real,   &
                      global_maxval_scalar_int
   end interface

   interface global_minval
     module procedure global_minval_dbl,           &
                      global_minval_real,          &
                      global_minval_int,           &
                      global_minval_scalar_dbl,    &
                      global_minval_scalar_real,   &
                      global_minval_scalar_int
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   logical(log_kind) :: ltripole_grid  ! in lieu of use domain

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_global_reductions
! !INTERFACE:

 subroutine init_global_reductions(tripole_flag)

! !DESCRIPTION:
!  Initializes necessary buffers for global reductions.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:
!
   logical(log_kind), intent(in) :: tripole_flag
!
!EOP
!BOC

! This flag is apparently never used; if it were used, it might need
! a corresponding tripoleTFlag to be defined.
   ltripole_grid = tripole_flag

!EOC

 end subroutine init_global_reductions

! Subprogram not used  subroutine sum_vector_dbl(local_vector,global_vector, dist)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  this function returns the sum of vector value across processors
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! distribution from which this is called
! Subprogram not used 
! Subprogram not used    real (dbl_kind), intent(inout) :: &
! Subprogram not used       local_vector(:)                ! local vector to be compared
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: global_vector(:)   ! resulting global sum
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: ierr ! MPI error flag
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: len
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    len = size(local_vector)
! Subprogram not used    if (dist%nprocs > 1) then
! Subprogram not used       if (my_task < dist%nprocs) then
! Subprogram not used          call MPI_ALLREDUCE(local_vector, global_vector, len, &
! Subprogram not used                             mpiR8, MPI_SUM, dist%communicator, ierr)
! Subprogram not used       else
! Subprogram not used          global_vector = c0
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       global_vector = local_vector
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end subroutine sum_vector_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_dbl(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind), dimension(:), allocatable :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum

   integer (int_kind) :: &
      i,j,iblock,n, &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      ierr,         &! mpi error flag
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      nreduce,      &! mpi count
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

!-----------------------------------------------------------------------

   nreduce = 1
   allocate(blockSum(nreduce), &
            globalSumTmp(nreduce))
   blockSum     = 0.0_dbl_kind
   globalSumTmp = 0.0_dbl_kind
   globalSum    = 0.0_dbl_kind

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      n = 1

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum(n) = &
            blockSum(n) + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum(n) = &
               blockSum(n) + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum(n) = blockSum(n) + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum(n) = &
                     blockSum(n) - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     blockSum(n) = blockSum(n) - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum(n) = blockSum(n) - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif
   end do

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(blockSum, globalSumTmp, nreduce, &
                         mpiR8, MPI_SUM, communicator, ierr)
   endif

   do n=1,nreduce
      globalSum = globalSum + globalSumTmp(n)
   enddo

   deallocate(blockSum, globalSumTmp)

!-----------------------------------------------------------------------
!EOC

 end function global_sum_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

! Subprogram not used  function global_sum_real(array, dist, field_loc, mMask, lMask) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of the physical domain of a 2-d array.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_sum
! Subprogram not used !  function corresponding to real arrays.  The generic
! Subprogram not used !  interface is identical but will handle real and integer 2-d slabs
! Subprogram not used !  and real, integer, and double precision scalars.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_loc            ! location of field on staggered grid
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       mMask                ! optional multiplicative mask
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       blockSum,     &! sum of local block domain
! Subprogram not used       localSum       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,   &! local counters
! Subprogram not used       ib,ie,jb,je,  &! beg,end of physical domain
! Subprogram not used       ierr,         &! mpi error flag
! Subprogram not used       blockID,      &! block location
! Subprogram not used       numProcs,     &! number of processor participating
! Subprogram not used       numBlocks,    &! number of local blocks
! Subprogram not used       communicator, &! communicator for this distribution
! Subprogram not used       maxiglob       ! maximum non-redundant value of i_global
! Subprogram not used  
! Subprogram not used    logical (log_kind) :: &
! Subprogram not used       Nrow           ! this field is on a N row (a velocity row)
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block     ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localSum  = 0.0_real_kind
! Subprogram not used    globalSum = 0.0_real_kind
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist,          &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,       &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockSum = 0.0_real_kind
! Subprogram not used 
! Subprogram not used       if (present(mMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = &
! Subprogram not used             blockSum + array(i,j,iblock)*mMask(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockSum = &
! Subprogram not used                blockSum + array(i,j,iblock)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = blockSum + array(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** if this row along or beyond tripole boundary
! Subprogram not used       !*** must eliminate redundant points from global sum
! Subprogram not used 
! Subprogram not used       if (this_block%tripole) then
! Subprogram not used          Nrow=(field_loc == field_loc_Nface .or. &
! Subprogram not used             field_loc == field_loc_NEcorner)
! Subprogram not used          if (Nrow .and. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = 0 ! entire u-row on T-fold grid
! Subprogram not used          elseif (Nrow .or. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
! Subprogram not used          else
! Subprogram not used             maxiglob = -1 ! nothing to do for T-row on u-fold
! Subprogram not used          endif
! Subprogram not used  
! Subprogram not used          if (maxiglob > 0) then
! Subprogram not used 
! Subprogram not used             j = je
! Subprogram not used 
! Subprogram not used             if (present(mMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = &
! Subprogram not used                      blockSum - array(i,j,iblock)*mMask(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else if (present(lMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      if (lMask(i,j,iblock)) &
! Subprogram not used                      blockSum = blockSum - array(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = blockSum - array(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** now add block sum to global sum
! Subprogram not used 
! Subprogram not used       localSum = localSum + blockSum
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localSum, globalSum, 1, &
! Subprogram not used                          mpiR4, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_real

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

! Subprogram not used  function global_sum_int(array, dist, field_loc, mMask, lMask) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of the physical domain of a 2-d array.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_sum
! Subprogram not used !  function corresponding to integer arrays.  The generic
! Subprogram not used !  interface is identical but will handle real and integer 2-d slabs
! Subprogram not used !  and real, integer, and double precision scalars.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_loc            ! location of field on staggered grid
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       mMask                ! optional multiplicative mask
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
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
! Subprogram not used       blockSum,     &! sum of local block domain
! Subprogram not used       localSum       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,   &! local counters
! Subprogram not used       ib,ie,jb,je,  &! beg,end of physical domain
! Subprogram not used       ierr,         &! mpi error flag
! Subprogram not used       blockID,      &! block location
! Subprogram not used       numProcs,     &! number of processor participating
! Subprogram not used       numBlocks,    &! number of local blocks
! Subprogram not used       communicator, &! communicator for this distribution
! Subprogram not used       maxiglob       ! maximum non-redundant value of i_global
! Subprogram not used 
! Subprogram not used    logical (log_kind) :: &
! Subprogram not used       Nrow           ! this field is on a N row (a velocity row)
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block     ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localSum  = 0_int_kind
! Subprogram not used    globalSum = 0_int_kind
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist,          &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,       &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockSum = 0
! Subprogram not used 
! Subprogram not used       if (present(mMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = &
! Subprogram not used             blockSum + array(i,j,iblock)*mMask(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockSum = &
! Subprogram not used                blockSum + array(i,j,iblock)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = blockSum + array(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** if this row along or beyond tripole boundary
! Subprogram not used       !*** must eliminate redundant points from global sum
! Subprogram not used 
! Subprogram not used       if (this_block%tripole) then
! Subprogram not used          Nrow=(field_loc == field_loc_Nface .or. &
! Subprogram not used             field_loc == field_loc_NEcorner)
! Subprogram not used          if (Nrow .and. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = 0 ! entire u-row on T-fold grid
! Subprogram not used          elseif (Nrow .or. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
! Subprogram not used          else
! Subprogram not used             maxiglob = -1 ! nothing to do for T-row on u-fold
! Subprogram not used          endif
! Subprogram not used  
! Subprogram not used          if (maxiglob > 0) then
! Subprogram not used 
! Subprogram not used             j = je
! Subprogram not used 
! Subprogram not used             if (present(mMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = &
! Subprogram not used                      blockSum - array(i,j,iblock)*mMask(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else if (present(lMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      if (lMask(i,j,iblock)) &
! Subprogram not used                      blockSum = blockSum - array(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = blockSum - array(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** now add block sum to global sum
! Subprogram not used 
! Subprogram not used       localSum = localSum + blockSum
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localSum, globalSum, 1, &
! Subprogram not used                          MPI_INTEGER, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_int

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

! Subprogram not used  function global_sum_scalar_dbl(scalar, dist) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of a set of scalars distributed across
! Subprogram not used !  a parallel machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_sum
! Subprogram not used !  function corresponding to double precision scalars.  The generic
! Subprogram not used !  interface is identical but will handle real and integer 2-d slabs
! Subprogram not used !  and real, integer, and double precision scalars.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
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
! Subprogram not used       ierr,         &! mpi error flag
! Subprogram not used       numProcs,     &! number of processor participating
! Subprogram not used       numBlocks,    &! number of local blocks
! Subprogram not used       communicator   ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !#ifdef REPRODUCIBLE
! Subprogram not used !   real (r16_kind) :: &
! Subprogram not used !      scalarTmp, globalSumTmp  ! higher precision for reproducibility
! Subprogram not used !#endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  get communicator for MPI calls
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !  REPRODUCIBLE option is commented out because MPI does not handle 
! Subprogram not used !  REAL16 correctly.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !#ifdef REPRODUCIBLE
! Subprogram not used !   if (my_task < numProcs) then
! Subprogram not used !      scalarTmp = scalar
! Subprogram not used !      call MPI_ALLREDUCE(scalarTmp, globalSumTmp, 1, &
! Subprogram not used !                         mpiR16, MPI_SUM, communicator, ierr)
! Subprogram not used !      globalSum = globalSumTmp
! Subprogram not used !   endif
! Subprogram not used !#else
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalSum, 1, &
! Subprogram not used                          mpiR8, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used !#endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_scalar_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

! Subprogram not used  function global_sum_scalar_real(scalar, dist) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of a set of scalars distributed across
! Subprogram not used !  a parallel machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_sum
! Subprogram not used !  function corresponding to real scalars.  The generic
! Subprogram not used !  interface is identical but will handle real and integer 2-d slabs
! Subprogram not used !  and real, integer, and double precision scalars.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
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
! Subprogram not used       ierr,         &! mpi error flag
! Subprogram not used       numProcs,     &! number of processor participating
! Subprogram not used       numBlocks,    &! number of local blocks
! Subprogram not used       communicator   ! communicator for this distribution
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  get communicator for MPI calls
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalSum, 1, &
! Subprogram not used                          mpiR4, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_scalar_real

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

! Subprogram not used  function global_sum_scalar_int(scalar, dist) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of a set of scalars distributed across
! Subprogram not used !  a parallel machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_sum
! Subprogram not used !  function corresponding to integer scalars.  The generic
! Subprogram not used !  interface is identical but will handle real and integer 2-d slabs
! Subprogram not used !  and real, integer, and double precision scalars.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
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
! Subprogram not used       ierr,         &! mpi error flag
! Subprogram not used       numProcs,     &! number of processor participating
! Subprogram not used       numBlocks,    &! number of local blocks
! Subprogram not used       communicator   ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  get communicator for MPI calls
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalSum, 1, &
! Subprogram not used                          MPI_INTEGER, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_scalar_int

!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

! Subprogram not used  function global_sum_prod_dbl (array1, array2, dist, field_loc, &
! Subprogram not used                                mMask, lMask) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of the physical domain of a product of
! Subprogram not used !  two 2-d arrays.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic 
! Subprogram not used !  global_sum_prod function corresponding to double precision arrays.
! Subprogram not used !  The generic interface is identical but will handle real and integer 
! Subprogram not used !  2-d slabs.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array1, array2       ! arrays whose product is to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for arrays
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_loc            ! location of field on staggered grid
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       mMask                ! optional multiplicative mask
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:), allocatable :: &
! Subprogram not used       blockSum,     &! sum of local block domain
! Subprogram not used       localSum,     &! sum of all local block domains
! Subprogram not used       globalSumTmp   ! higher precision global sum
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,n,    &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       blockID,         &! block location
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       nreduce,         &! mpi count
! Subprogram not used       maxiglob          ! maximum non-redundant value of i_global
! Subprogram not used 
! Subprogram not used    logical (log_kind) :: &
! Subprogram not used       Nrow           ! this field is on a N row (a velocity row)
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block     ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    nreduce = 1
! Subprogram not used    allocate(blockSum(nreduce), &
! Subprogram not used             globalSumTmp(nreduce))
! Subprogram not used    blockSum     = 0.0_dbl_kind
! Subprogram not used    globalSumTmp = 0.0_dbl_kind
! Subprogram not used    globalSum    = 0.0_dbl_kind
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       n = 1
! Subprogram not used 
! Subprogram not used       if (present(mMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum(n) = &
! Subprogram not used             blockSum(n) + array1(i,j,iblock)*array2(i,j,iblock)* &
! Subprogram not used                        mMask(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockSum(n) = &
! Subprogram not used                blockSum(n) + array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum(n) = blockSum(n) + array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** if this row along or beyond tripole boundary
! Subprogram not used       !*** must eliminate redundant points from global sum
! Subprogram not used 
! Subprogram not used       if (this_block%tripole) then
! Subprogram not used          Nrow=(field_loc == field_loc_Nface .or. &
! Subprogram not used             field_loc == field_loc_NEcorner)
! Subprogram not used          if (Nrow .and. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = 0 ! entire u-row on T-fold grid
! Subprogram not used          elseif (Nrow .or. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
! Subprogram not used          else
! Subprogram not used             maxiglob = -1 ! nothing to do for T-row on u-fold
! Subprogram not used          endif
! Subprogram not used  
! Subprogram not used          if (maxiglob > 0) then
! Subprogram not used 
! Subprogram not used             j = je
! Subprogram not used 
! Subprogram not used             if (present(mMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum(n) = &
! Subprogram not used                      blockSum(n) - array1(i,j,iblock)*array2(i,j,iblock)* &
! Subprogram not used                                 mMask(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else if (present(lMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      if (lMask(i,j,iblock)) &
! Subprogram not used                         blockSum(n) = blockSum(n) - &
! Subprogram not used                                    array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum(n) = blockSum(n) - &
! Subprogram not used                                 array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(blockSum, globalSumTmp, nreduce, &
! Subprogram not used                          mpiR8, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    do n=1,nreduce
! Subprogram not used       globalSum = globalSum + globalSumTmp(n)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    deallocate(blockSum, globalSumTmp)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_prod_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

! Subprogram not used  function global_sum_prod_real (array1, array2, dist, field_loc, &
! Subprogram not used                                 mMask, lMask) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of the physical domain of a product of
! Subprogram not used !  two 2-d arrays.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic 
! Subprogram not used !  global_sum_prod function corresponding to single precision arrays.
! Subprogram not used !  The generic interface is identical but will handle real and integer 
! Subprogram not used !  2-d slabs.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array1, array2       ! arrays whose product is to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for arrays
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_loc            ! location of field on staggered grid
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       mMask                ! optional multiplicative mask
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       blockSum,     &! sum of local block domain
! Subprogram not used       localSum       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       blockID,         &! block location
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       maxiglob          ! maximum non-redundant value of i_global
! Subprogram not used  
! Subprogram not used    logical (log_kind) :: &
! Subprogram not used       Nrow           ! this field is on a N row (a velocity row)
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localSum  = 0.0_real_kind
! Subprogram not used    globalSum = 0.0_real_kind
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockSum = 0.0_real_kind
! Subprogram not used 
! Subprogram not used       if (present(mMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = &
! Subprogram not used             blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
! Subprogram not used                        mMask(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockSum = &
! Subprogram not used                blockSum + array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** if this row along or beyond tripole boundary
! Subprogram not used       !*** must eliminate redundant points from global sum
! Subprogram not used 
! Subprogram not used       if (this_block%tripole) then
! Subprogram not used          Nrow=(field_loc == field_loc_Nface .or. &
! Subprogram not used             field_loc == field_loc_NEcorner)
! Subprogram not used          if (Nrow .and. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = 0 ! entire u-row on T-fold grid
! Subprogram not used          elseif (Nrow .or. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
! Subprogram not used          else
! Subprogram not used             maxiglob = -1 ! nothing to do for T-row on u-fold
! Subprogram not used          endif
! Subprogram not used  
! Subprogram not used          if (maxiglob > 0) then
! Subprogram not used 
! Subprogram not used             j = je
! Subprogram not used 
! Subprogram not used             if (present(mMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = &
! Subprogram not used                      blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
! Subprogram not used                                 mMask(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else if (present(lMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      if (lMask(i,j,iblock)) &
! Subprogram not used                         blockSum = blockSum - &
! Subprogram not used                                    array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = blockSum - &
! Subprogram not used                                 array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** now add block sum to global sum
! Subprogram not used 
! Subprogram not used       localSum = localSum + blockSum
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localSum, globalSum, 1, &
! Subprogram not used                          mpiR4, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_prod_real

!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

! Subprogram not used  function global_sum_prod_int (array1, array2, dist, field_loc, &
! Subprogram not used                                mMask, lMask) &
! Subprogram not used           result(globalSum)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global sum of the physical domain of a product of
! Subprogram not used !  two 2-d arrays.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic 
! Subprogram not used !  global_sum_prod function corresponding to integer arrays.
! Subprogram not used !  The generic interface is identical but will handle real and integer 
! Subprogram not used !  2-d slabs.
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array1, array2       ! arrays whose product is to be summed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for arrays
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       field_loc            ! location of field on staggered grid
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       mMask                ! optional multiplicative mask
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       globalSum            ! resulting global sum
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
! Subprogram not used       blockSum,     &! sum of local block domain
! Subprogram not used       localSum       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       blockID,         &! block location
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       maxiglob          ! maximum non-redundant value of i_global
! Subprogram not used  
! Subprogram not used    logical (log_kind) :: &
! Subprogram not used       Nrow           ! this field is on a N row (a velocity row)
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localSum  = 0_int_kind
! Subprogram not used    globalSum = 0_int_kind
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockSum = 0
! Subprogram not used 
! Subprogram not used       if (present(mMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = &
! Subprogram not used             blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
! Subprogram not used                        mMask(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockSum = &
! Subprogram not used                blockSum + array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** if this row along or beyond tripole boundary
! Subprogram not used       !*** must eliminate redundant points from global sum
! Subprogram not used 
! Subprogram not used       if (this_block%tripole) then
! Subprogram not used          Nrow=(field_loc == field_loc_Nface .or. &
! Subprogram not used             field_loc == field_loc_NEcorner)
! Subprogram not used          if (Nrow .and. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = 0 ! entire u-row on T-fold grid
! Subprogram not used          elseif (Nrow .or. this_block%tripoleTFlag) then
! Subprogram not used             maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
! Subprogram not used          else
! Subprogram not used             maxiglob = -1 ! nothing to do for T-row on u-fold
! Subprogram not used          endif
! Subprogram not used  
! Subprogram not used          if (maxiglob > 0) then
! Subprogram not used 
! Subprogram not used             j = je
! Subprogram not used 
! Subprogram not used             if (present(mMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = &
! Subprogram not used                      blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
! Subprogram not used                                 mMask(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else if (present(lMask)) then
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      if (lMask(i,j,iblock)) &
! Subprogram not used                         blockSum = blockSum - &
! Subprogram not used                                    array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             else
! Subprogram not used                do i=ib,ie
! Subprogram not used                   if (this_block%i_glob(i) > maxiglob) then
! Subprogram not used                      blockSum = blockSum - &
! Subprogram not used                                 array1(i,j,iblock)*array2(i,j,iblock)
! Subprogram not used                   endif
! Subprogram not used                end do
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !*** now add block sum to global sum
! Subprogram not used 
! Subprogram not used       localSum = localSum + blockSum
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local sum to global sum
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localSum, globalSum, 1, &
! Subprogram not used                          MPI_INTEGER, MPI_SUM, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end function global_sum_prod_int

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

! Subprogram not used  function global_maxval_dbl (array, dist, lMask) &
! Subprogram not used           result(globalMaxval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global maximum value of the physical domain of a 2-d field
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_maxval
! Subprogram not used !  function corresponding to double precision arrays.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array for which max value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       globalMaxval         ! resulting maximum value of array
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (dbl_kind) ::    &
! Subprogram not used       blockMaxval,     &! sum of local block domain
! Subprogram not used       localMaxval       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       blockID           ! block location
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localMaxval  = -HUGE(0.0_dbl_kind)
! Subprogram not used    globalMaxval = -HUGE(0.0_dbl_kind)
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockMaxval = -HUGE(0.0_dbl_kind)
! Subprogram not used 
! Subprogram not used       if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockMaxval = max(blockMaxval,array(i,j,iblock))
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockMaxval = max(blockMaxval,array(i,j,iblock))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       localMaxval = max(localMaxval,blockMaxval)
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local maxval to global maxval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
! Subprogram not used                          mpiR8, MPI_MAX, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_maxval_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

! Subprogram not used  function global_maxval_real (array, dist, lMask) &
! Subprogram not used           result(globalMaxval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global maximum value of the physical domain of a 2-d field
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_maxval
! Subprogram not used !  function corresponding to single precision arrays.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array for which max value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalMaxval         ! resulting maximum value of array
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (real_kind) ::    &
! Subprogram not used       blockMaxval,     &! sum of local block domain
! Subprogram not used       localMaxval       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       blockID           ! block location
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localMaxval  = -HUGE(0.0_real_kind)
! Subprogram not used    globalMaxval = -HUGE(0.0_real_kind)
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockMaxval = -HUGE(0.0_real_kind)
! Subprogram not used 
! Subprogram not used       if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockMaxval = max(blockMaxval,array(i,j,iblock))
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockMaxval = max(blockMaxval,array(i,j,iblock))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       localMaxval = max(localMaxval,blockMaxval)
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local maxval to global maxval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
! Subprogram not used                          mpiR4, MPI_MAX, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_maxval_real

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

! Subprogram not used  function global_maxval_int (array, dist, lMask) &
! Subprogram not used           result(globalMaxval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global maximum value of the physical domain of a 2-d field
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_maxval
! Subprogram not used !  function corresponding to integer arrays.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array for which max value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       globalMaxval         ! resulting maximum value of array
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) ::    &
! Subprogram not used       blockMaxval,     &! sum of local block domain
! Subprogram not used       localMaxval       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       blockID           ! block location
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localMaxval  = -HUGE(0_int_kind)
! Subprogram not used    globalMaxval = -HUGE(0_int_kind)
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockMaxval = -HUGE(0_int_kind)
! Subprogram not used 
! Subprogram not used       if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockMaxval = max(blockMaxval,array(i,j,iblock))
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockMaxval = max(blockMaxval,array(i,j,iblock))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       localMaxval = max(localMaxval,blockMaxval)
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local maxval to global maxval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
! Subprogram not used                          MPI_INTEGER, MPI_MAX, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_maxval_int

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

! Subprogram not used  function global_maxval_scalar_dbl (scalar, dist) &
! Subprogram not used           result(globalMaxval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global maximum value of a scalar value across
! Subprogram not used !  a distributed machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_maxval
! Subprogram not used !  function corresponding to double precision scalars.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar for which max value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       globalMaxval         ! resulting maximum value
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
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator      ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local maxval to global maxval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
! Subprogram not used                          mpiR8, MPI_MAX, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_maxval_scalar_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

! Subprogram not used  function global_maxval_scalar_real (scalar, dist) &
! Subprogram not used           result(globalMaxval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global maximum value of a scalar value across
! Subprogram not used !  a distributed machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_maxval
! Subprogram not used !  function corresponding to single precision scalars.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar for which max value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalMaxval         ! resulting maximum value
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
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator      ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local maxval to global maxval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
! Subprogram not used                          mpiR4, MPI_MAX, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_maxval_scalar_real

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_scalar_int (scalar, dist) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision scalars.  

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalMaxval         ! resulting maximum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_int

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

! Subprogram not used  function global_minval_dbl (array, dist, lMask) &
! Subprogram not used           result(globalMinval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global minimum value of the physical domain of a 2-d field
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_minval
! Subprogram not used !  function corresponding to double precision arrays.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array for which min value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       globalMinval         ! resulting minimum value of array
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (dbl_kind) ::    &
! Subprogram not used       blockMinval,     &! sum of local block domain
! Subprogram not used       localMinval       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       blockID           ! block location
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localMinval  = HUGE(0.0_dbl_kind)
! Subprogram not used    globalMinval = HUGE(0.0_dbl_kind)
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockMinval = HUGE(0.0_dbl_kind)
! Subprogram not used 
! Subprogram not used       if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockMinval = min(blockMinval,array(i,j,iblock))
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockMinval = min(blockMinval,array(i,j,iblock))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       localMinval = min(localMinval,blockMinval)
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local minval to global minval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
! Subprogram not used                          mpiR8, MPI_MIN, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_minval_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

! Subprogram not used  function global_minval_real (array, dist, lMask) &
! Subprogram not used           result(globalMinval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global minimum value of the physical domain of a 2-d field
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_minval
! Subprogram not used !  function corresponding to single precision arrays.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array for which min value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalMinval         ! resulting minimum value of array
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real (real_kind) ::    &
! Subprogram not used       blockMinval,     &! sum of local block domain
! Subprogram not used       localMinval       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       blockID           ! block location
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localMinval  = HUGE(0.0_real_kind)
! Subprogram not used    globalMinval = HUGE(0.0_real_kind)
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockMinval = HUGE(0.0_real_kind)
! Subprogram not used 
! Subprogram not used       if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockMinval = min(blockMinval,array(i,j,iblock))
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockMinval = min(blockMinval,array(i,j,iblock))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       localMinval = min(localMinval,blockMinval)
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local minval to global minval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
! Subprogram not used                          mpiR4, MPI_MIN, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_minval_real

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

! Subprogram not used  function global_minval_int (array, dist, lMask) &
! Subprogram not used           result(globalMinval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global minimum value of the physical domain of a 2-d field
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_minval
! Subprogram not used !  function corresponding to integer arrays.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), dimension(:,:,:), intent(in) :: &
! Subprogram not used       array                ! array for which min value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution for array X
! Subprogram not used 
! Subprogram not used    logical (log_kind), dimension(:,:,:), intent(in), optional :: &
! Subprogram not used       lMask                ! optional logical mask
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       globalMinval         ! resulting minimum value of array
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  local variables
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    integer (int_kind) ::    &
! Subprogram not used       blockMinval,     &! sum of local block domain
! Subprogram not used       localMinval       ! sum of all local block domains
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       i,j,iblock,      &! local counters
! Subprogram not used       ib,ie,jb,je,     &! beg,end of physical domain
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numBlocks,       &! number of local blocks
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator,    &! communicator for this distribution
! Subprogram not used       blockID           ! block location
! Subprogram not used 
! Subprogram not used    type (block) :: &
! Subprogram not used       this_block          ! holds local block information
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    localMinval  = HUGE(0_int_kind)
! Subprogram not used    globalMinval = HUGE(0_int_kind)
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             numLocalBlocks = numBlocks, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used    do iblock=1,numBlocks
! Subprogram not used       call ice_distributionGetBlockID(dist, iblock, blockID)
! Subprogram not used 
! Subprogram not used       this_block = get_block(blockID, blockID)
! Subprogram not used 
! Subprogram not used       ib = this_block%ilo
! Subprogram not used       ie = this_block%ihi
! Subprogram not used       jb = this_block%jlo
! Subprogram not used       je = this_block%jhi
! Subprogram not used 
! Subprogram not used       blockMinval = HUGE(0_int_kind)
! Subprogram not used 
! Subprogram not used       if (present(lMask)) then
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             if (lMask(i,j,iblock)) then
! Subprogram not used                blockMinval = min(blockMinval,array(i,j,iblock))
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       else
! Subprogram not used          do j=jb,je
! Subprogram not used          do i=ib,ie
! Subprogram not used             blockMinval = min(blockMinval,array(i,j,iblock))
! Subprogram not used          end do
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       localMinval = min(localMinval,blockMinval)
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local minval to global minval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
! Subprogram not used                          MPI_INTEGER, MPI_MIN, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_minval_int

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

! Subprogram not used  function global_minval_scalar_dbl (scalar, dist) &
! Subprogram not used           result(globalMinval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global minimum value of a scalar value across
! Subprogram not used !  a distributed machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_minval
! Subprogram not used !  function corresponding to double precision scalars.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar for which min value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       globalMinval         ! resulting minimum value
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
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator      ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local minval to global minval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalMinval, 1, &
! Subprogram not used                          mpiR8, MPI_MIN, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_minval_scalar_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

! Subprogram not used  function global_minval_scalar_real (scalar, dist) &
! Subprogram not used           result(globalMinval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global minimum value of a scalar value across
! Subprogram not used !  a distributed machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_minval
! Subprogram not used !  function corresponding to single precision scalars.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar for which min value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    real (real_kind) :: &
! Subprogram not used       globalMinval         ! resulting minimum value
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
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator      ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local minval to global minval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalMinval, 1, &
! Subprogram not used                          mpiR4, MPI_MIN, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_minval_scalar_real

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

! Subprogram not used  function global_minval_scalar_int (scalar, dist) &
! Subprogram not used           result(globalMinval)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Computes the global minimum value of a scalar value across
! Subprogram not used !  a distributed machine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used !
! Subprogram not used ! !REMARKS:
! Subprogram not used !  This is actually the specific interface for the generic global_minval
! Subprogram not used !  function corresponding to single precision scalars.  
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       scalar               ! scalar for which min value needed
! Subprogram not used 
! Subprogram not used    type (distrb), intent(in) :: &
! Subprogram not used       dist                 ! block distribution
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind) :: &
! Subprogram not used       globalMinval         ! resulting minimum value
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
! Subprogram not used       ierr,            &! mpi error flag
! Subprogram not used       numProcs,        &! number of processor participating
! Subprogram not used       communicator      ! communicator for this distribution
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call ice_distributionGet(dist, &
! Subprogram not used                             nprocs = numProcs,        &
! Subprogram not used                             communicator = communicator)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  now use MPI global reduction to reduce local minval to global minval
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (my_task < numProcs) then
! Subprogram not used       call MPI_ALLREDUCE(scalar, globalMinval, 1, &
! Subprogram not used                          MPI_INTEGER, MPI_MIN, communicator, ierr)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used  end function global_minval_scalar_int

!***********************************************************************

 end module ice_global_reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
