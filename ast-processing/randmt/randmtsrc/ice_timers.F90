!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_timers

!BOP
! !MODULE: timers
!
! !DESCRIPTION:
!  This module contains routine for supporting multiple CPU timers
!  and accumulates time for each individual block and node (task).
!
! !REVISION HISTORY:
!  SVN:$Id: ice_timers.F90 144 2008-08-12 21:37:19Z eclare $
!
! 2005: Adapted from POP by William Lipscomb
!       Replaced 'stdout' by 'nu_diag'
! 2006 ECH: Replaced 'system_clock' timing mechanism by 'MPI_WTIME'
!           for MPI runs.  Single-processor runs still use system_clock.
!
! !USES:

   use ice_kinds_mod
   use ice_constants
   use ice_domain
   use ice_global_reductions
   use ice_exit
   use ice_fileunits, only: nu_diag, nu_timing
   use ice_communicate, only: my_task, master_task, get_num_procs, lprint_stats
   use ice_gather_scatter, only: gatherArray

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_ice_timers,     &
             get_ice_timer,       &
             ice_timer_clear,     &
             ice_timer_start,     &
             ice_timer_stop,      &
             ice_timer_print,     &
             ice_timer_print_all, &
             ice_timer_check

!EOP
!BOC

!-----------------------------------------------------------------------
! public timers
!-----------------------------------------------------------------------

   integer (int_kind), public ::      &
      timer_total,            &! total time
      timer_step,             &! time stepping
      timer_dynamics,         &! dynamics
      timer_advect,           &! horizontal advection
      timer_column,           &! column
      timer_thermo,           &! thermodynamics
      timer_sw,               &! radiative transfer
      timer_ridge,            &! ridging
      timer_catconv,          &! category conversions
      timer_couple,           &! coupling
      timer_readwrite,        &! read/write
      timer_diags,            &! diagnostics/history
      timer_hist,             &! diagnostics/history
      timer_cplrecv,          &! receive from coupler
      timer_rcvsnd,           &! time between receive to send
      timer_cplsend,          &! send to 1
      timer_sndrcv,           &! time between send to receive
      timer_bound              ! boundary updates
!      timer_tmp                ! for temporary timings

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      max_timers = 50          ! max number of timers

   type timer_data
      character (char_len) :: &
         name                  ! timer name

      logical (log_kind) ::   &
         in_use,              &! true if timer initialized
         node_started          ! true if any thread has started timer

      integer (int_kind) ::   &
         num_blocks,          &! number of blocks using this timer
         num_nodes,           &! number of nodes  using this timer
         num_starts,          &! number of start requests
         num_stops             ! number of stop requests

      real (dbl_kind) :: &
         node_cycles1,        &! cycle number at start for node timer
         node_cycles2          ! cycle number at stop  for node timer

      real (dbl_kind) ::            &
         node_accum_time       ! accumulated time for node timer

      logical (log_kind), dimension(:), pointer :: &
         block_started         ! true if block timer started

      real (dbl_kind), dimension(:), pointer :: &
         block_cycles1,        &! cycle number at start for block timers
         block_cycles2          ! cycle number at stop  for block timers

      real (dbl_kind), dimension(:), pointer :: &
         block_accum_time       ! accumulated time for block timers

   end type

   type (timer_data), dimension(max_timers) :: &
      all_timers               ! timer data for all timers

   real (dbl_kind) ::               &
      clock_rate               ! clock rate in seconds for each cycle

   !----------------------------------------------
   ! some arrays on which to collect timing info
   !---------------------------------------------
   integer (int_kind), public :: timerRoot ! MPI process ID to collect timing information

   real(dbl_kind), public  :: all_ltime(max_timers) ! local times for each timer
   real(dbl_kind), allocatable :: all_gtime(:)      ! global times for each timer


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_ice_timers
! !INTERFACE:

 subroutine init_ice_timers

! !DESCRIPTION:
!  This routine initializes machine parameters and timer structures
!  for computing cpu time from F90 intrinsic timer functions.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n ! dummy loop index

!-----------------------------------------------------------------------
!
!  initialize timer structures
!
!-----------------------------------------------------------------------

   return

   clock_rate = c1

   timerRoot = min(distrb_info%nprocs-1,2)

   do n=1,max_timers
      all_timers(n)%name = 'unknown_timer_name'

      all_timers(n)%in_use       = .false.
      all_timers(n)%node_started = .false.

      all_timers(n)%num_blocks   = 0
      all_timers(n)%num_nodes    = 0
      all_timers(n)%num_starts   = 0
      all_timers(n)%num_stops    = 0
      all_timers(n)%node_cycles1 = c0
      all_timers(n)%node_cycles2 = c0

      all_timers(n)%node_accum_time = c0

      nullify(all_timers(n)%block_started)
      nullify(all_timers(n)%block_cycles1)
      nullify(all_timers(n)%block_cycles2)
      nullify(all_timers(n)%block_accum_time)
   end do

   call get_ice_timer(timer_total,    'Total',    nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_step,     'TimeLoop', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_dynamics, 'Dynamics', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_advect,   'Advection',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_column,   'Column',   nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_thermo,   'Thermo',   nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_sw,       'Shortwave',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_ridge,    'Ridging',  nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_catconv,  'Cat Conv', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_readwrite,'ReadWrite',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_diags,    'Diags    ',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_hist,     'History  ',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_bound,    'Bound',    nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_cplrecv,  'Cpl-Imp', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_cplsend,  'Cpl-Exp', nblocks,distrb_info%nprocs)
!   call get_ice_timer(timer_tmp,      '         ',nblocks,distrb_info%nprocs)

   !------------------------------------------------------
   ! allocate the array of timer values from all processes
   !------------------------------------------------------
   if(my_task .eq. timerRoot) then
     allocate(all_gtime(max_timers*distrb_info%nprocs))
   else
     allocate(all_gtime(1))
   endif


!-----------------------------------------------------------------------
!EOC

   end subroutine init_ice_timers

!***********************************************************************
!BOP
! !IROUTINE: get_ice_timer
! !INTERFACE:

! Subprogram not used  subroutine get_ice_timer(timer_id, name_choice, num_blocks, num_nodes)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This routine initializes a timer with a given name and returns a 
! Subprogram not used !  timer id.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    character (*), intent(in) :: &
! Subprogram not used       name_choice               ! input name for this timer
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       num_nodes,               &! number of nodes(tasks) using this timer
! Subprogram not used       num_blocks                ! number of blocks using this timer
! Subprogram not used                                 ! (can be =1 if timer called outside
! Subprogram not used                                 !  threaded region)
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(out) :: &
! Subprogram not used       timer_id           ! timer number assigned to this timer 
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
! Subprogram not used       n,                 &! dummy loop index
! Subprogram not used       srch_error          ! error flag for search
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  search for next free timer
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    srch_error = 1
! Subprogram not used 
! Subprogram not used    srch_loop: do n=1,max_timers
! Subprogram not used       if (.not. all_timers(n)%in_use) then
! Subprogram not used          srch_error = 0
! Subprogram not used          timer_id = n
! Subprogram not used 
! Subprogram not used          all_timers(n)%name       = ' '
! Subprogram not used          all_timers(n)%name       = name_choice
! Subprogram not used          all_timers(n)%in_use     = .true.
! Subprogram not used          all_timers(n)%num_blocks = num_blocks
! Subprogram not used          all_timers(n)%num_nodes  = num_nodes 
! Subprogram not used 
! Subprogram not used          allocate(all_timers(n)%block_started   (num_blocks), &
! Subprogram not used                   all_timers(n)%block_cycles1   (num_blocks), &
! Subprogram not used                   all_timers(n)%block_cycles2   (num_blocks), &
! Subprogram not used                   all_timers(n)%block_accum_time(num_blocks))
! Subprogram not used 
! Subprogram not used          all_timers(n)%block_started    = .false.
! Subprogram not used          all_timers(n)%block_cycles1    = c0
! Subprogram not used          all_timers(n)%block_cycles2    = c0
! Subprogram not used          all_timers(n)%block_accum_time = c0
! Subprogram not used 
! Subprogram not used          exit srch_loop
! Subprogram not used       endif
! Subprogram not used    end do srch_loop
! Subprogram not used 
! Subprogram not used    if (srch_error /= 0) &
! Subprogram not used       call abort_ice('get_ice_timer: Exceeded maximum number of timers')
! Subprogram not used                     
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end subroutine get_ice_timer

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_clear
! !INTERFACE:

! Subprogram not used  subroutine ice_timer_clear(timer_id)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This routine resets the time for a timer which has already been
! Subprogram not used !  defined.  NOTE: This routine must be called from outside a threaded
! Subprogram not used !  region to ensure correct reset of block timers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       timer_id                ! timer number
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if the timer has been defined, reset all times to 0
! Subprogram not used !  otherwise exit with an error
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    if (all_timers(timer_id)%in_use) then
! Subprogram not used       all_timers(timer_id)%node_started  = .false.
! Subprogram not used       all_timers(timer_id)%num_starts    = 0
! Subprogram not used       all_timers(timer_id)%num_stops     = 0
! Subprogram not used       all_timers(timer_id)%node_cycles1  = c0
! Subprogram not used       all_timers(timer_id)%node_cycles2  = c0
! Subprogram not used 
! Subprogram not used       all_timers(timer_id)%node_accum_time = c0
! Subprogram not used 
! Subprogram not used       all_timers(timer_id)%block_started(:)    = .false.
! Subprogram not used       all_timers(timer_id)%block_cycles1(:)    = c0
! Subprogram not used       all_timers(timer_id)%block_cycles2(:)    = c0
! Subprogram not used       all_timers(timer_id)%block_accum_time(:) = c0
! Subprogram not used    else
! Subprogram not used       call abort_ice &
! Subprogram not used                  ('ice_timer_clear: attempt to reset undefined timer')
! Subprogram not used                     
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end subroutine ice_timer_clear

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_start
! !INTERFACE:

 subroutine ice_timer_start(timer_id, block_id)
 use perf_mod

! !DESCRIPTION:
!  This routine starts a given node timer if it has not already
!  been started by another thread.  If block information is available,
!  the appropriate block timer is also started.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                   ! timer number

   integer (int_kind), intent(in), optional :: &
      block_id                 ! optional block id for this block
                               ! this must be the actual local address
                               ! of the block in the distribution
                               ! from which it is called
                               ! (if timer called outside of block
                               ! region, no block info required)

   double precision MPI_WTIME
   external MPI_WTIME

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if timer is defined, start it up
!
!-----------------------------------------------------------------------

   return

   if (all_timers(timer_id)%in_use) then

      call t_startf(trim(all_timers(timer_id)%name))

      !***
      !*** if called from within a block loop, start block timers
      !***

      if (present(block_id)) then

         !*** if block timer already started, stop it first

         if (all_timers(timer_id)%block_started(block_id)) &
            call ice_timer_stop(timer_id, block_id)

         !*** start block timer

         all_timers(timer_id)%block_started(block_id) = .true.
         all_timers(timer_id)%block_cycles1(block_id) = MPI_WTIME()

         !*** start node timer if not already started by
         !*** another thread.  if already started, keep track
         !*** of number of start requests in order to match
         !*** start and stop requests
 
         !$OMP CRITICAL

         if (.not. all_timers(timer_id)%node_started) then
            all_timers(timer_id)%node_started = .true.
            all_timers(timer_id)%num_starts   = 1
            all_timers(timer_id)%num_stops    = 0
            all_timers(timer_id)%node_cycles1 = MPI_WTIME()
         else
            all_timers(timer_id)%num_starts = &
            all_timers(timer_id)%num_starts + 1
         endif

         !$OMP END CRITICAL

      !***
      !*** if called from outside a block loop, start node timer
      !***

      else

         !*** stop timer if already started
         if (all_timers(timer_id)%node_started)  &
                                        call ice_timer_stop(timer_id)

         !*** start node timer

         all_timers(timer_id)%node_started = .true.
         all_timers(timer_id)%node_cycles1 = MPI_WTIME()

      endif
   else
      call abort_ice &
                 ('ice_timer_start: attempt to start undefined timer')
                    
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_start
 
!***********************************************************************
!BOP
! !IROUTINE: ice_timer_stop
! !INTERFACE:

 subroutine ice_timer_stop(timer_id, block_id)
 use perf_mod

! !DESCRIPTION:
!  This routine stops a given node timer if appropriate.  If block 
!  information is available the appropriate block timer is also stopped.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                   ! timer number

   integer (int_kind), intent(in), optional :: &
      block_id                 ! optional block id for this block
                               ! this must be the actual local address
                               ! of the block in the distribution
                               ! from which it is called
                               ! (if timer called outside of block
                               ! region, no block info required)

   double precision MPI_WTIME
   external MPI_WTIME

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) :: &
      cycles1, cycles2   ! temps to hold cycle info before correction

!-----------------------------------------------------------------------
!
!  get end cycles
!
!-----------------------------------------------------------------------

   return
   cycles2 = MPI_WTIME()

!-----------------------------------------------------------------------
!
!  if timer is defined, stop it
!
!-----------------------------------------------------------------------

   if (all_timers(timer_id)%in_use) then

      !***
      !*** if called from within a block loop, stop block timer
      !***

      if (present(block_id)) then

         all_timers(timer_id)%block_started(block_id) = .false.

         cycles1 = all_timers(timer_id)%block_cycles1(block_id)
         all_timers(timer_id)%block_accum_time(block_id) = &
         all_timers(timer_id)%block_accum_time(block_id) + &
            clock_rate*(cycles2 - cycles1)

         !*** stop node timer if number of requested stops
         !*** matches the number of starts (to avoid stopping
         !*** a node timer started by multiple threads)
 
         cycles1 = all_timers(timer_id)%node_cycles1

         !$OMP CRITICAL

         all_timers(timer_id)%num_stops = &
         all_timers(timer_id)%num_stops + 1

         if (all_timers(timer_id)%num_starts == &
             all_timers(timer_id)%num_stops) then

            all_timers(timer_id)%node_started = .false.
            all_timers(timer_id)%node_accum_time = &
            all_timers(timer_id)%node_accum_time + &
               clock_rate*(cycles2 - cycles1)

            all_timers(timer_id)%num_starts   = 0
            all_timers(timer_id)%num_stops    = 0

            all_ltime(timer_id) = all_timers(timer_id)%node_accum_time

         endif

         !$OMP END CRITICAL

      !***
      !*** if called from outside a block loop, stop node timer
      !***

      else

         all_timers(timer_id)%node_started = .false.
         cycles1 = all_timers(timer_id)%node_cycles1

         all_timers(timer_id)%node_accum_time = &
         all_timers(timer_id)%node_accum_time + &
            clock_rate*(cycles2 - cycles1)

         all_ltime(timer_id) = all_timers(timer_id)%node_accum_time

      endif

      call t_stopf(trim(all_timers(timer_id)%name))
   else
      call abort_ice &
                 ('ice_timer_stop: attempt to stop undefined timer')
                    
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_stop
 
!***********************************************************************
!BOP
! !IROUTINE: ice_timer_print
! !INTERFACE:

! Subprogram not used  subroutine ice_timer_print(timer_id,stats)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  Prints the accumulated time for a given timer and optional
! Subprogram not used !  statistics for that timer. It is assumed that this routine
! Subprogram not used !  is called outside of a block loop.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       timer_id                ! timer number
! Subprogram not used 
! Subprogram not used    logical (log_kind), intent(in), optional :: &
! Subprogram not used       stats                   ! if true, print statistics for node
! Subprogram not used                               !   and block times for this timer
! Subprogram not used 
! Subprogram not used    double precision MPI_WTIME
! Subprogram not used    external MPI_WTIME
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
! Subprogram not used       n,icount,        & ! dummy loop index and counter
! Subprogram not used       nBlocks            
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       local_time,       &! temp space for holding local timer results
! Subprogram not used       min_time,         &! minimum accumulated time
! Subprogram not used       max_time,         &! maximum accumulated time
! Subprogram not used       mean_time          ! mean    accumulated time
! Subprogram not used 
! Subprogram not used    real (dbl_kind) :: &
! Subprogram not used       cycles1, cycles2   ! temps to hold cycle info before correction
! Subprogram not used 
! Subprogram not used    character (41), parameter :: &
! Subprogram not used       timer_format = "('Timer ',i3,': ',a20,f11.2,' seconds')"
! Subprogram not used 
! Subprogram not used    character (49), parameter :: &
! Subprogram not used       stats_fmt1 = "('  Timer stats (node): min = ',f11.2,' seconds')",&
! Subprogram not used       stats_fmt2 = "('                      max = ',f11.2,' seconds')",&
! Subprogram not used       stats_fmt3 = "('                      mean= ',f11.2,' seconds')",&
! Subprogram not used       stats_fmt4 = "('  Timer stats(block): min = ',f11.2,' seconds')"
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  if timer has been defined, check to see whether it is currently
! Subprogram not used !  running.  If it is, update to current running time and print the info
! Subprogram not used !  (without stopping the timer). 
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    if (all_timers(timer_id)%in_use) then
! Subprogram not used 
! Subprogram not used       !*** Find max node time and print that time as default timer
! Subprogram not used       !*** result
! Subprogram not used 
! Subprogram not used       if (my_task < all_timers(timer_id)%num_nodes) then
! Subprogram not used 
! Subprogram not used          if (all_timers(timer_id)%node_started) then
! Subprogram not used             cycles2 = MPI_WTIME()
! Subprogram not used             cycles1 = all_timers(timer_id)%node_cycles1
! Subprogram not used             local_time = all_timers(timer_id)%node_accum_time + &
! Subprogram not used                          clock_rate*(cycles2 - cycles1)
! Subprogram not used          else
! Subprogram not used             local_time = all_timers(timer_id)%node_accum_time
! Subprogram not used          endif
! Subprogram not used       else
! Subprogram not used          local_time = c0
! Subprogram not used       endif
! Subprogram not used       max_time = global_maxval(local_time,distrb_info)
! Subprogram not used       
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used         write (nu_diag,timer_format) timer_id, &
! Subprogram not used               trim(all_timers(timer_id)%name),max_time
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (present(stats)) then
! Subprogram not used       if (stats) then
! Subprogram not used 
! Subprogram not used          !*** compute and print statistics for node timer
! Subprogram not used 
! Subprogram not used          min_time = global_minval(local_time,distrb_info)
! Subprogram not used          mean_time = global_sum(local_time,distrb_info)/ &
! Subprogram not used                      real(all_timers(timer_id)%num_nodes,kind=dbl_kind)
! Subprogram not used          if (my_task == master_task) then
! Subprogram not used             write (nu_diag,stats_fmt1) min_time
! Subprogram not used             write (nu_diag,stats_fmt2) max_time
! Subprogram not used             write (nu_diag,stats_fmt3) mean_time
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          !*** compute and print statistics for block timers
! Subprogram not used          !*** min block time
! Subprogram not used 
! Subprogram not used          local_time = bignum
! Subprogram not used          do n=1,all_timers(timer_id)%num_blocks
! Subprogram not used             local_time = min(local_time, &
! Subprogram not used                              all_timers(timer_id)%block_accum_time(n))
! Subprogram not used          end do
! Subprogram not used          min_time = global_minval(local_time,distrb_info)
! Subprogram not used          if (min_time == bignum) min_time = c0
! Subprogram not used 
! Subprogram not used          !*** max block time
! Subprogram not used 
! Subprogram not used          local_time = -bignum
! Subprogram not used          do n=1,all_timers(timer_id)%num_blocks
! Subprogram not used             local_time = max(local_time, &
! Subprogram not used                              all_timers(timer_id)%block_accum_time(n))
! Subprogram not used          end do
! Subprogram not used          max_time = global_maxval(local_time,distrb_info)
! Subprogram not used          if (max_time == -bignum) min_time = c0
! Subprogram not used 
! Subprogram not used          !*** mean block time
! Subprogram not used 
! Subprogram not used          local_time = c0
! Subprogram not used          nBlocks = all_timers(timer_id)%num_blocks
! Subprogram not used          do n=1,nBlocks
! Subprogram not used             local_time = local_time + &
! Subprogram not used                          all_timers(timer_id)%block_accum_time(n)
! Subprogram not used          end do
! Subprogram not used          icount = global_sum(nBlocks, distrb_info)
! Subprogram not used          if (icount > 0) mean_time=global_sum(local_time,distrb_info)&
! Subprogram not used                                    /real(icount,kind=dbl_kind)
! Subprogram not used 
! Subprogram not used          if (my_task == master_task) then
! Subprogram not used             write (nu_diag,stats_fmt4) min_time
! Subprogram not used             write (nu_diag,stats_fmt2) max_time
! Subprogram not used             write (nu_diag,stats_fmt3) mean_time
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used       call abort_ice &
! Subprogram not used                  ('ice_timer_print: attempt to print undefined timer')
! Subprogram not used                     
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end subroutine ice_timer_print

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_print_all
! !INTERFACE:

 subroutine ice_timer_print_all(stats)

! !DESCRIPTION:
!  Prints the accumulated time for a all timers and optional
!  statistics for that timer. It is assumed that this routine
!  is called outside of a block loop.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in), optional :: &
      stats                   ! if true, print statistics for node
                              !   and block times for this timer

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: nprocs

   integer (int_kind) :: n ! dummy loop index

!-----------------------------------------------------------------------
!
!  loop through timers anc call timer_print for each defined timer
!
!-----------------------------------------------------------------------

   return
   if (my_task == master_task) then
      write(nu_diag,'(/,a19,/)') 'Timing information:'
   endif

   do n=1,max_timers
      if (all_timers(n)%in_use) then
         if (present(stats)) then
            call ice_timer_print(n,stats)
         else
            call ice_timer_print(n)
         endif
      endif
   end do

   !-----------------------------------------------------
   ! gather all timing values onto the timeRoot processor
   !-----------------------------------------------------
   call gatherArray(all_gtime,all_ltime,max_timers,timerRoot)

   !--------------------------
   ! write out the timing data
   !--------------------------
   if(my_task == timerRoot) then
     if(lprint_stats) then
        nprocs = get_num_procs()
        open(nu_timing,file='timing.bin',recl=8*max_timers*nprocs, &
           form = 'unformatted', access = 'direct', status='unknown')
        write(nu_timing,rec=1) all_gtime
        close(nu_timing)
     endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_print_all

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_check
! !INTERFACE:

! Subprogram not used  subroutine ice_timer_check(timer_id,block_id)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !  This routine checks a given timer by stopping and restarting the
! Subprogram not used !  timer.  This is primarily used to periodically accumulate time in 
! Subprogram not used !  the timer to prevent timer cycles from wrapping around max_cycles.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  same as module
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in) :: &
! Subprogram not used       timer_id                   ! timer number
! Subprogram not used 
! Subprogram not used    integer (int_kind), intent(in), optional :: &
! Subprogram not used       block_id                 ! optional block id for this block
! Subprogram not used                                ! this must be the actual local address
! Subprogram not used                                ! of the block in the distribution
! Subprogram not used                                ! from which it is called
! Subprogram not used                                ! (if timer called outside of block
! Subprogram not used                                ! region, no block info required)
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !BOC
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  stop and restart the requested timer
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    if (present(block_id)) then
! Subprogram not used       call ice_timer_stop (timer_id,block_id)
! Subprogram not used       call ice_timer_start(timer_id,block_id)
! Subprogram not used    else
! Subprogram not used       call ice_timer_stop (timer_id)
! Subprogram not used       call ice_timer_start(timer_id)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !EOC
! Subprogram not used 
! Subprogram not used  end subroutine ice_timer_check

!***********************************************************************

 end module ice_timers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
