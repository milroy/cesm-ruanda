!=======================================================================
!BOP
!
! !MODULE: ice_transport_driver - drivers for ice transport
!
! !DESCRIPTION:
!
! Drivers for remapping and upwind ice transport
!
! !REVISION HISTORY:
!  SVN:$Id: ice_transport_upwind.F 28 2006-11-03 20:32:53Z eclare $
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL 
!
! 2004: Revised by William Lipscomb from ice_transport_mpdata.
!       Stripped out mpdata, retained upwind, and added block structure.
! 2006: Incorporated remap transport driver and renamed from
!       ice_transport_upwind.  
!
! !INTERFACE:

      module ice_transport_driver
!
! !USES:
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use perf_mod,        only: t_startf, t_stopf, t_barrierf
!
!EOP
!
      implicit none
      save

      character (len=char_len) ::     &
         advection   ! type of advection scheme used
                     ! 'upwind' => 1st order donor cell scheme
                     ! 'remap' => remapping scheme

      logical, parameter :: & ! if true, prescribe area flux across each edge  
         l_fixed_area = .false.

! NOTE: For remapping, hice, hsno, qice, and qsno are considered tracers.
!       max_ntrace is not equal to max_ntrcr!
!       ntrace is not equal to ntrcr!

      integer (kind=int_kind), parameter ::                      &
         max_ntrace = 2+max_ntrcr+nilyr+nslyr  ! hice,hsno,qice,qsno,trcr

      integer (kind=int_kind) ::                      &
         ntrace              ! number of tracers in use
                          
      integer (kind=int_kind), dimension (:), allocatable ::     &
         tracer_type       ,&! = 1, 2, or 3 (see comments below)
         depend              ! tracer dependencies (see below)

      logical (kind=log_kind), dimension (:), allocatable ::     &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), parameter ::                      &
         integral_order = 3   ! polynomial order of quadrature integrals
                              ! linear=1, quadratic=2, cubic=3

      logical (kind=log_kind), parameter ::     &
         l_dp_midpt = .true.  ! if true, find departure points using
                              ! corrected midpoint velocity
                          
!=======================================================================

      contains

!=======================================================================

!BOP
!
! !IROUTINE: init_transport - initializations for horizontal transport
!
! !INTERFACE:
!
      subroutine init_transport
!
! !DESCRIPTION:
!
! This subroutine is a wrapper for init_remap, which initializes the
! remapping transport scheme.  If the model is run with upwind
! transport, no initializations are necessary.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !USES:
!
      use ice_state, only: trcr_depend, ntrcr
      use ice_exit
      use ice_timers
      use ice_transport_remap, only: init_remap
!
!EOP
!
      integer (kind=int_kind) ::       &
         k, nt, nt1     ! tracer indices

      call ice_timer_start(timer_advect)  ! advection 

      ntrace = 2+ntrcr+nilyr+nslyr  ! hice,hsno,qice,qsno,trcr

      allocate (tracer_type   (ntrace), &
                depend        (ntrace), &
                has_dependents(ntrace))

      if (trim(advection)=='remap') then

!lipscomb - two branches for now; consolidate later

         ! define tracer dependency arrays
         ! see comments in remapping routine

          depend(1:2)         = 0 ! hice, hsno
          tracer_type(1:2)    = 1 ! no dependency
      
          k = 2

          do nt = 1, ntrcr
             depend(k+nt) = trcr_depend(nt) ! 0 for ice area tracers
                                            ! 1 for ice volume tracers
                                            ! 2 for snow volume tracers
             if (trcr_depend(nt) == 0) then
                tracer_type(k+nt) = 1
             else               ! trcr_depend = 1 or 2
                tracer_type(k+nt) = 2
             endif
          enddo

          k = k + ntrcr
          
          depend(k+1:k+nilyr) = 1 ! qice depends on hice
          tracer_type(k+1:k+nilyr) = 2 

          k = k + nilyr

          depend(k+1:k+nslyr) = 2 ! qsno depends on hsno
          tracer_type(k+1:k+nslyr) = 2 

          has_dependents = .false.
          do nt = 1, ntrace
             if (depend(nt) > 0) then
                nt1 = depend(nt)
                has_dependents(nt1) = .true.
                if (nt1 > nt) then
                   write(nu_diag,*)     &
                      'Tracer nt2 =',nt,' depends on tracer nt1 =',nt1
                   call abort_ice       &
                      ('ice: remap transport: Must have nt2 > nt1')
                endif
             endif
          enddo                 ! ntrace

          call init_remap    ! grid quantities

      else   ! upwind

         continue

      endif

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine init_transport

!=======================================================================
!BOP
!
! !IROUTINE: transport_remap - wrapper for remapping transport scheme
!
! !INTERFACE:
!
! Subprogram not used       subroutine transport_remap (dt)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! This subroutine solves the transport equations for one timestep
! Subprogram not used ! using the conservative remapping scheme developed by John Dukowicz
! Subprogram not used ! and John Baumgardner (DB) and modified for sea ice by William
! Subprogram not used ! Lipscomb and Elizabeth Hunke.
! Subprogram not used !
! Subprogram not used ! This scheme preserves monotonicity of ice area and tracers.  That is,
! Subprogram not used ! it does not produce new extrema.  It is second-order accurate in space,
! Subprogram not used ! except where gradients are limited to preserve monotonicity. 
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_boundary
! Subprogram not used       use ice_global_reductions
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_blocks
! Subprogram not used       use ice_state
! Subprogram not used       use ice_grid, only: tarea, HTE, HTN
! Subprogram not used       use ice_exit
! Subprogram not used       use ice_calendar, only: istep1, istep, diagfreq
! Subprogram not used       use ice_timers
! Subprogram not used       use ice_transport_remap, only: horizontal_remap, make_masks
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) ::     &
! Subprogram not used          dt      ! time step
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       ! local variables
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          i, j           ,&! horizontal indices
! Subprogram not used          iblk           ,&! block index
! Subprogram not used          ilo,ihi,jlo,jhi,&! beginning and end of physical domain
! Subprogram not used          n              ,&! ice category index
! Subprogram not used          nt, nt1, nt2     ! tracer indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind),      &
! Subprogram not used          dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
! Subprogram not used          aim            ,&! mean ice category areas in each grid cell
! Subprogram not used          aimask           ! = 1. if ice is present, = 0. otherwise
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind),      &
! Subprogram not used          dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
! Subprogram not used          trm            ,&! mean tracer values in each grid cell
! Subprogram not used          trmask           ! = 1. if tracer is present, = 0. otherwise
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) ::     &
! Subprogram not used          l_stop           ! if true, abort the model
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          istop, jstop     ! indices of grid cell where model aborts 
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
! Subprogram not used          icellsnc         ! number of cells with ice
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind),      &
! Subprogram not used          dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
! Subprogram not used          indxinc, indxjnc   ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used       
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! If l_fixed_area is true, the area of each departure region is
! Subprogram not used     !  computed in advance (e.g., by taking the divergence of the 
! Subprogram not used     !  velocity field and passed to locate_triangles.  The departure 
! Subprogram not used     !  regions are adjusted to obtain the desired area.
! Subprogram not used     ! If false, edgearea is computed in locate_triangles and passed out.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) ::   &
! Subprogram not used          edgearea_e     ,&! area of departure regions for east edges
! Subprogram not used          edgearea_n       ! area of departure regions for north edges
! Subprogram not used 
! Subprogram not used       ! variables related to optional bug checks
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), parameter ::     &
! Subprogram not used          l_conservation_check = .true. ,&! if true, check conservation
! Subprogram not used          l_monotonicity_check = .true.   ! if true, check monotonicity
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(0:ncat) ::     &
! Subprogram not used          asum_init      ,&! initial global ice area
! Subprogram not used          asum_final       ! final global ice area
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(max_ntrace,ncat) ::     &
! Subprogram not used          atsum_init     ,&! initial global ice area*tracer
! Subprogram not used          atsum_final      ! final global ice area*tracer
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable ::     &
! Subprogram not used          tmin         ,&! local min tracer
! Subprogram not used          tmax           ! local max tracer
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: alloc_error
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
! Subprogram not used          work1
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_advect)  ! advection 
! Subprogram not used 
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used !---! Prepare for remapping.
! Subprogram not used !---! Initialize, update ghost cells, fill tracer arrays.
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute open water area in each grid cell.
! Subprogram not used     ! Note: An aggregate_area call is needed only if the open
! Subprogram not used     !       water area has changed since the previous call.
! Subprogram not used     !       Here we assume that aice0 is up to date.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      do iblk = 1, nblocks
! Subprogram not used !         call aggregate_area (nx_block, ny_block,
! Subprogram not used !                              iblk,     &
! Subprogram not used !                              aicen(:,:,:,iblk),     &
! Subprogram not used !                              aice (:,:,  iblk),     &
! Subprogram not used !                              aice0(:,:,  iblk)) 
! Subprogram not used !      enddo
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Ghost cell updates for state variables.
! Subprogram not used     ! Commented out because ghost cells are updated after cleanup_itd.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used !      call ice_timer_start(timer_bound)
! Subprogram not used 
! Subprogram not used !      call ice_HaloUpdate (aice0,            halo_info,     &
! Subprogram not used !                           field_loc_center, field_type_scalar)
! Subprogram not used 
! Subprogram not used !      call bound_state (aicen, trcrn,     &
! Subprogram not used !                        vicen, vsnon,      &
! Subprogram not used !                        eicen, esnon)
! Subprogram not used 
! Subprogram not used !      call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Ghost cell updates for ice velocity.
! Subprogram not used     ! Commented out because ghost cell velocities are computed
! Subprogram not used     !  in ice_dyn_evp.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      call ice_timer_start(timer_bound)
! Subprogram not used !      call ice_HaloUpdate (uvel,               halo_info,     &
! Subprogram not used !                           field_loc_NEcorner, field_type_vector)
! Subprogram not used !      call ice_HaloUpdate (vvel,               halo_info,     &
! Subprogram not used !                           field_loc_NEcorner, field_type_vector)
! Subprogram not used !      call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call t_barrierf('cice_remap_s2t_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_remap_s2t')
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Fill arrays with fields to be remapped.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call state_to_tracers(nx_block,          ny_block,             &
! Subprogram not used                                ntrcr,             ntrace,               &
! Subprogram not used                                aice0(:,:,  iblk),                       &
! Subprogram not used                                aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
! Subprogram not used                                vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
! Subprogram not used                                eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
! Subprogram not used                                aim  (:,:,:,iblk), trm  (:,:,:,:,iblk))
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf   ('cice_remap_s2t')
! Subprogram not used       call t_barrierf('cice_remap_check1_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_remap_check1')
! Subprogram not used 
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used !---! Optional conservation and monotonicity checks.
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (l_conservation_check .and. mod(istep,diagfreq) == 0) then
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute initial values of globally conserved quantities.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do n = 0, ncat
! Subprogram not used             asum_init(n) = global_sum(aim(:,:,n,:),     distrb_info,       &
! Subprogram not used                                       field_loc_center, tarea)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          do n = 1, ncat
! Subprogram not used             do nt = 1, ntrace
! Subprogram not used                if (tracer_type(nt)==1) then ! does not depend on another tracer
! Subprogram not used                   atsum_init(nt,n) =      &
! Subprogram not used                       global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
! Subprogram not used                                       distrb_info,     field_loc_center,   &
! Subprogram not used                                       tarea)
! Subprogram not used                elseif (tracer_type(nt)==2) then ! depends on another tracer
! Subprogram not used                   nt1 = depend(nt)
! Subprogram not used                   do iblk = 1, nblocks
! Subprogram not used                      do j= 1,ny_block  
! Subprogram not used                         do i = 1,nx_block
! Subprogram not used                            work1(i,j,iblk) = trm(i,j,nt,n,iblk)*trm(i,j,nt1,n,iblk)
! Subprogram not used                         end do
! Subprogram not used                      end do
! Subprogram not used                   end do
! Subprogram not used                   atsum_init(nt,n) =     &
! Subprogram not used                       global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
! Subprogram not used                                       distrb_info,  field_loc_center,      &
! Subprogram not used                                       tarea)
! Subprogram not used                elseif (tracer_type(nt)==3) then ! depends on two tracers
! Subprogram not used                   nt1 = depend(nt)
! Subprogram not used                   nt2 = depend(nt1)
! Subprogram not used                   do iblk = 1, nblocks
! Subprogram not used                      do j= 1,ny_block  
! Subprogram not used                         do i = 1,nx_block
! Subprogram not used                            work1(i,j,iblk) = trm(i,j,nt,n,iblk)*trm(i,j,nt1,n,iblk) &
! Subprogram not used 	                                                       *trm(i,j,nt2,n,iblk)
! Subprogram not used                         end do
! Subprogram not used                      end do
! Subprogram not used                   end do
! Subprogram not used                   atsum_init(nt,n) =     &
! Subprogram not used                       global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
! Subprogram not used                                       distrb_info,  field_loc_center,      &
! Subprogram not used                                       tarea)
! Subprogram not used                endif            ! tracer_type
! Subprogram not used             enddo               ! nt
! Subprogram not used          enddo                  ! n
! Subprogram not used 
! Subprogram not used       endif                     ! l_conservation_check
! Subprogram not used       
! Subprogram not used       if (l_monotonicity_check .and. mod(istep,diagfreq) == 0) then
! Subprogram not used 
! Subprogram not used          allocate(tmin(nx_block,ny_block,ntrace,ncat,max_blocks),     &
! Subprogram not used                   tmax(nx_block,ny_block,ntrace,ncat,max_blocks),     &
! Subprogram not used                   STAT=alloc_error)
! Subprogram not used 
! Subprogram not used          if (alloc_error /= 0)      &
! Subprogram not used               call abort_ice ('ice: allocation error')
! Subprogram not used 
! Subprogram not used          tmin(:,:,:,:,:) = c0
! Subprogram not used          tmax(:,:,:,:,:) = c0
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             ilo = this_block%ilo
! Subprogram not used             ihi = this_block%ihi
! Subprogram not used             jlo = this_block%jlo
! Subprogram not used             jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used     ! Compute masks.
! Subprogram not used     ! Masks are used to prevent tracer values in cells without ice
! Subprogram not used     !  from being used in the monotonicity check.
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used             call make_masks (nx_block,          ny_block,              &
! Subprogram not used                              ilo, ihi,          jlo, jhi,              &
! Subprogram not used                              nghost,            ntrace,                &
! Subprogram not used                              has_dependents,    icellsnc(:,iblk),      &
! Subprogram not used                              indxinc(:,:,iblk), indxjnc(:,:,iblk),     &
! Subprogram not used                              aim(:,:,:,iblk),   aimask(:,:,:,iblk),    &
! Subprogram not used                              trm(:,:,:,:,iblk), trmask(:,:,:,:,iblk))
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute local max and min of tracer fields.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             do n = 1, ncat
! Subprogram not used                call local_max_min                                      &  
! Subprogram not used                             (nx_block,           ny_block,             &
! Subprogram not used                              ilo, ihi,           jlo, jhi,             &
! Subprogram not used                              trm (:,:,:,n,iblk),                       &
! Subprogram not used                              tmin(:,:,:,n,iblk), tmax  (:,:,:,n,iblk), &
! Subprogram not used                              aimask(:,:,n,iblk), trmask(:,:,:,n,iblk))
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          call ice_timer_start(timer_bound)
! Subprogram not used          call ice_HaloUpdate (tmin,             halo_info,     &
! Subprogram not used                               field_loc_center, field_type_scalar)
! Subprogram not used          call ice_HaloUpdate (tmax,             halo_info,     &
! Subprogram not used                               field_loc_center, field_type_scalar)
! Subprogram not used          call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             ilo = this_block%ilo
! Subprogram not used             ihi = this_block%ihi
! Subprogram not used             jlo = this_block%jlo
! Subprogram not used             jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used             do n = 1, ncat
! Subprogram not used                call quasilocal_max_min (nx_block, ny_block,     &
! Subprogram not used                                         ilo, ihi, jlo, jhi,     &
! Subprogram not used                                         tmin(:,:,:,n,iblk),      &
! Subprogram not used                                         tmax(:,:,:,n,iblk))
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       endif                     ! l_monotonicity_check
! Subprogram not used 
! Subprogram not used       call t_stopf   ('cice_remap_check1')
! Subprogram not used       call t_barrierf('cice_remap_horz_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_remap_horz')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Main remapping routine: Step ice area and tracers forward in time.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! If l_fixed_area is true, compute edgearea by taking the divergence
! Subprogram not used     !  of the velocity field.  Otherwise, initialize edgearea.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             edgearea_e(i,j,iblk) = c0
! Subprogram not used             edgearea_n(i,j,iblk) = c0
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          if (l_fixed_area) then
! Subprogram not used 
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk,this_block,i,j,ilo,ihi,jlo,jhi)
! Subprogram not used             do iblk = 1, nblocks
! Subprogram not used                this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used                ilo = this_block%ilo
! Subprogram not used                ihi = this_block%ihi
! Subprogram not used                jlo = this_block%jlo
! Subprogram not used                jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used                do j = jlo, jhi
! Subprogram not used                do i = ilo-1, ihi
! Subprogram not used                   edgearea_e(i,j,iblk) = (uvel(i,j,iblk) + uvel(i,j-1,iblk)) &
! Subprogram not used                                         * p5 * HTE(i,j,iblk) * dt
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used 
! Subprogram not used                do j = jlo-1, jhi
! Subprogram not used                do i = ilo, ihi
! Subprogram not used                   edgearea_n(i,j,iblk) = (vvel(i,j,iblk) + vvel(i-1,j,iblk)) &
! Subprogram not used                                         * p5 * HTN(i,j,iblk) * dt
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used 
! Subprogram not used             enddo  ! iblk
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          call horizontal_remap (dt,                ntrace,             &
! Subprogram not used                                 uvel      (:,:,:), vvel      (:,:,:),  &
! Subprogram not used                                 aim     (:,:,:,:), trm   (:,:,:,:,:),  &
! Subprogram not used                                 l_fixed_area,                          &
! Subprogram not used                                 edgearea_e(:,:,:), edgearea_n(:,:,:),  &
! Subprogram not used                                 tracer_type,       depend,             &
! Subprogram not used                                 has_dependents,    integral_order,     &
! Subprogram not used                                 l_dp_midpt)
! Subprogram not used 
! Subprogram not used       call t_stopf   ('cice_remap_horz')
! Subprogram not used       call t_barrierf('cice_remap_t2s_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_remap_t2s')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Given new fields, recompute state variables.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used          call tracers_to_state (nx_block,          ny_block,            &
! Subprogram not used                                 ntrcr,             ntrace,              &
! Subprogram not used                                 aim  (:,:,:,iblk), trm  (:,:,:,:,iblk), &
! Subprogram not used                                 aice0(:,:,  iblk),                      &
! Subprogram not used                                 aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk), &
! Subprogram not used                                 vicen(:,:,:,iblk), vsnon(:,:,  :,iblk), &
! Subprogram not used                                 eicen(:,:,:,iblk), esnon(:,:,  :,iblk)) 
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf   ('cice_remap_t2s')
! Subprogram not used       call t_barrierf('cice_remap_bound1_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_remap_bound1')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Ghost cell updates for state variables.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used 
! Subprogram not used       call bound_state (aicen, trcrn,     &
! Subprogram not used                         vicen, vsnon,      &
! Subprogram not used                         eicen, esnon)
! Subprogram not used 
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used       call t_stopf   ('cice_remap_bound1')
! Subprogram not used       call t_barrierf('cice_remap_check2_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_remap_check2')
! Subprogram not used 
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used !---! Optional conservation and monotonicity checks
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute final values of globally conserved quantities.
! Subprogram not used     ! Check global conservation of area and area*tracers.  (Optional)
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (l_conservation_check .and. mod(istep,diagfreq) == 0) then
! Subprogram not used 
! Subprogram not used          do n = 0, ncat
! Subprogram not used             asum_final(n) = global_sum(aim(:,:,n,:),     distrb_info,      &
! Subprogram not used                                        field_loc_center, tarea)
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          do n = 1, ncat
! Subprogram not used             do nt = 1, ntrace
! Subprogram not used                if (tracer_type(nt)==1) then ! does not depend on another tracer
! Subprogram not used                   atsum_final(nt,n) =      &
! Subprogram not used                       global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
! Subprogram not used                                       distrb_info,     field_loc_center,   &
! Subprogram not used                                       tarea)
! Subprogram not used                elseif (tracer_type(nt)==2) then ! depends on another tracer
! Subprogram not used                   nt1 = depend(nt)
! Subprogram not used                   do iblk = 1, nblocks
! Subprogram not used                      do j= 1,ny_block  
! Subprogram not used                         do i = 1,nx_block
! Subprogram not used                            work1(i,j,iblk) = trm(i,j,nt,n,iblk)*trm(i,j,nt1,n,iblk)
! Subprogram not used                         end do
! Subprogram not used                      end do
! Subprogram not used                   end do
! Subprogram not used                   atsum_final(nt,n) =     &
! Subprogram not used                       global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
! Subprogram not used                                       distrb_info,  field_loc_center,      &
! Subprogram not used                                       tarea)
! Subprogram not used                elseif (tracer_type(nt)==3) then ! depends on two tracers
! Subprogram not used                   nt1 = depend(nt)
! Subprogram not used                   nt2 = depend(nt1)
! Subprogram not used                   do iblk = 1, nblocks
! Subprogram not used                      do j= 1,ny_block  
! Subprogram not used                         do i = 1,nx_block
! Subprogram not used                            work1(i,j,iblk) = trm(i,j,nt,n,iblk)*trm(i,j,nt1,n,iblk) &
! Subprogram not used 	                                                       *trm(i,j,nt2,n,iblk)
! Subprogram not used                         end do
! Subprogram not used                      end do
! Subprogram not used                   end do
! Subprogram not used                   atsum_final(nt,n) =     &
! Subprogram not used                       global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
! Subprogram not used                                       distrb_info,  field_loc_center,      &
! Subprogram not used                                       tarea)
! Subprogram not used                endif            ! tracer_type
! Subprogram not used             enddo               ! nt
! Subprogram not used          enddo                  ! n
! Subprogram not used 
! Subprogram not used 
! Subprogram not used          if (my_task == master_task) then
! Subprogram not used             call global_conservation (l_stop,     &
! Subprogram not used                                       asum_init(0), asum_final(0))
! Subprogram not used 
! Subprogram not used             if (l_stop) then
! Subprogram not used                write (nu_diag,*) 'istep1 =', istep1
! Subprogram not used                write (nu_diag,*) 'transport: conservation error, cat 0'
! Subprogram not used                l_stop = .false.
! Subprogram not used                call abort_ice('ice remap transport: conservation error')
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             do n = 1, ncat               
! Subprogram not used                call global_conservation                                 &
! Subprogram not used                                      (l_stop,                           &
! Subprogram not used                                       asum_init(n),    asum_final(n),   &
! Subprogram not used                                       atsum_init(:,n), atsum_final(:,n))
! Subprogram not used 
! Subprogram not used                if (l_stop) then
! Subprogram not used                   write (nu_diag,*) 'istep1, cat =',     &
! Subprogram not used                                      istep1, n
! Subprogram not used                   write (nu_diag,*) 'transport: conservation error, cat ',n
! Subprogram not used                   l_stop = .false.
! Subprogram not used                   call abort_ice     &
! Subprogram not used                        ('ice remap transport: conservation error')
! Subprogram not used                endif
! Subprogram not used             enddo               ! n
! Subprogram not used 
! Subprogram not used          endif                  ! my_task = master_task
! Subprogram not used 
! Subprogram not used       endif                     ! l_conservation_check
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Check tracer monotonicity.  (Optional)
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (l_monotonicity_check .and. mod(istep,diagfreq) == 0) then
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             ilo = this_block%ilo
! Subprogram not used             ihi = this_block%ihi
! Subprogram not used             jlo = this_block%jlo
! Subprogram not used             jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used             do n = 1, ncat
! Subprogram not used                call check_monotonicity      &
! Subprogram not used                                (nx_block,           ny_block,     &
! Subprogram not used                                 ilo, ihi, jlo, jhi,     &
! Subprogram not used                                 iblk,     &
! Subprogram not used                                 tmin(:,:,:,n,iblk), tmax(:,:,:,n,iblk),  &
! Subprogram not used                                 aim (:,:,  n,iblk), trm (:,:,:,n,iblk),  &
! Subprogram not used                                 l_stop,     &
! Subprogram not used                                 istop,              jstop)
! Subprogram not used 
! Subprogram not used                if (l_stop) then
! Subprogram not used                   write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
! Subprogram not used                                      istep1, my_task, iblk, n
! Subprogram not used                   write (nu_diag,*) 'i_glob, j_glob',this_block%i_glob(istop), &
! Subprogram not used                                                      this_block%j_glob(jstop)
! Subprogram not used                   call abort_ice('ice remap transport: monotonicity error')
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used             enddo               ! n
! Subprogram not used 
! Subprogram not used          enddo                  ! iblk
! Subprogram not used 
! Subprogram not used          deallocate(tmin, tmax, STAT=alloc_error)
! Subprogram not used          if (alloc_error /= 0) call abort_ice ('deallocation error')
! Subprogram not used 
! Subprogram not used       endif                     ! l_monotonicity_check
! Subprogram not used 
! Subprogram not used       call ice_timer_stop(timer_advect)  ! advection 
! Subprogram not used       call t_stopf   ('cice_remap_check2')
! Subprogram not used 
! Subprogram not used       end subroutine transport_remap

!=======================================================================
!BOP
!
! !IROUTINE: transport_upwind - upwind transport
!
! !INTERFACE:
!
! Subprogram not used       subroutine transport_upwind (dt)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Computes the transport equations for one timestep using upwind. Sets
! Subprogram not used ! several fields into a work array and passes it to upwind routine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! same as module
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_boundary
! Subprogram not used       use ice_blocks
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_state
! Subprogram not used       use ice_grid, only: HTE, HTN, tarea
! Subprogram not used       use ice_timers
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) ::     &
! Subprogram not used          dt      ! time step
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          narr               ! number of state variable arrays
! Subprogram not used                             ! not including eicen, esnon
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          i, j, iblk       ,&! horizontal indices
! Subprogram not used          ilo,ihi,jlo,jhi    ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,nblocks) ::     &
! Subprogram not used          uee, vnn           ! cell edge velocities
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind),     &
! Subprogram not used          dimension (:,:,:,:), allocatable ::      &
! Subprogram not used          works              ! work array
! Subprogram not used 
! Subprogram not used       type (block) ::     &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_advect)  ! advection 
! Subprogram not used 
! Subprogram not used       narr = 1 + ncat*(3+ntrcr) ! max number of state variable arrays
! Subprogram not used                                 ! not including eicen, esnon
! Subprogram not used 
! Subprogram not used       allocate (works(nx_block,ny_block,narr,max_blocks))
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Get ghost cell values of state variables.
! Subprogram not used     ! (Assume velocities are already known for ghost cells, also.)
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used !      call bound_state (aicen, trcrn,     &
! Subprogram not used !                        vicen, vsnon,     &
! Subprogram not used !                        eicen, esnon)
! Subprogram not used 
! Subprogram not used       uee(:,:,:) = c0
! Subprogram not used       vnn(:,:,:) = c0
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Average corner velocities to edges.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used       
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,i,j,ilo,ihi,jlo,jhi)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             uee(i,j,iblk) = p5*(uvel(i,j,iblk) + uvel(i,j-1,iblk))
! Subprogram not used             vnn(i,j,iblk) = p5*(vvel(i,j,iblk) + vvel(i-1,j,iblk))
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used       call ice_HaloUpdate (uee,             halo_info,     &
! Subprogram not used                            field_loc_Eface, field_type_vector)
! Subprogram not used       call ice_HaloUpdate (vnn,             halo_info,     &
! Subprogram not used                            field_loc_Nface, field_type_vector)
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! fill work arrays with fields to be advected
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call state_to_work (nx_block,             ny_block,             &
! Subprogram not used                              ntrcr,                                      &
! Subprogram not used                              narr,                 trcr_depend,          &
! Subprogram not used                              aicen (:,:,  :,iblk), trcrn (:,:,:,:,iblk), &
! Subprogram not used                              vicen (:,:,  :,iblk), vsnon (:,:,  :,iblk), &
! Subprogram not used                              aice0 (:,:,    iblk), works (:,:,  :,iblk))
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! advect
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call upwind_field (nx_block,       ny_block,               &
! Subprogram not used                             ilo, ihi,       jlo, jhi,               &
! Subprogram not used                             dt,                                     &
! Subprogram not used                             narr,           works(:,:,:,iblk),      &
! Subprogram not used                             uee(:,:,iblk),  vnn    (:,:,iblk),      &
! Subprogram not used                             HTE(:,:,iblk),  HTN    (:,:,iblk),      &
! Subprogram not used                             tarea(:,:,iblk))
! Subprogram not used 
! Subprogram not used          call upwind_field (nx_block,       ny_block,               &
! Subprogram not used                             ilo, ihi,       jlo, jhi,               &
! Subprogram not used                             dt,                                     &
! Subprogram not used                             ntilyr,         eicen(:,:,:,iblk),      &
! Subprogram not used                             uee(:,:,iblk),  vnn    (:,:,iblk),      &
! Subprogram not used                             HTE(:,:,iblk),  HTN    (:,:,iblk),      &
! Subprogram not used                             tarea(:,:,iblk))
! Subprogram not used 
! Subprogram not used          call upwind_field (nx_block,       ny_block,               &
! Subprogram not used                             ilo, ihi,       jlo, jhi,               &
! Subprogram not used                             dt,                                     &
! Subprogram not used                             ntslyr,         esnon(:,:,:,iblk),      &
! Subprogram not used                             uee(:,:,iblk),  vnn    (:,:,iblk),      &
! Subprogram not used                             HTE(:,:,iblk),  HTN    (:,:,iblk),      &
! Subprogram not used                             tarea(:,:,iblk))
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! convert work arrays back to state variables
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call work_to_state (nx_block,            ny_block,              &
! Subprogram not used                              ntrcr,                                      &
! Subprogram not used                              narr,                trcr_depend,           &
! Subprogram not used                              aicen(:,:,  :,iblk), trcrn (:,:,:,:,iblk),  &
! Subprogram not used                              vicen(:,:,  :,iblk), vsnon (:,:,  :,iblk),  &
! Subprogram not used                              aice0(:,:,    iblk), works (:,:,  :,iblk)) 
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used  
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Ghost cell updates for state variables.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used 
! Subprogram not used       call bound_state (aicen, trcrn,     &
! Subprogram not used                         vicen, vsnon,      &
! Subprogram not used                         eicen, esnon)
! Subprogram not used 
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       call ice_timer_stop(timer_advect)  ! advection 
! Subprogram not used 
! Subprogram not used       deallocate(works)
! Subprogram not used 
! Subprogram not used       end subroutine transport_upwind

!=======================================================================
! The next few subroutines (through check_monotonicity) are called
! by transport_remap.
!=======================================================================
!
!BOP
!
! !IROUTINE: state_to_tracers -fill ice area and tracer arrays
!
! !INTERFACE:
!
! Subprogram not used       subroutine state_to_tracers (nx_block, ny_block,   &
! Subprogram not used                                    ntrcr,    ntrace,     &
! Subprogram not used                                    aice0,                &
! Subprogram not used                                    aicen,    trcrn,      &
! Subprogram not used                                    vicen,    vsnon,      &
! Subprogram not used                                    eicen,    esnon,      &
! Subprogram not used                                    aim,      trm)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Fill ice area and tracer arrays.
! Subprogram not used ! Assume that the advected tracers are hicen, hsnon, trcrn, 
! Subprogram not used !  qicen(1:nilyr), and qsnon(1:nslyr).
! Subprogram not used ! This subroutine must be modified if a different set of tracers
! Subprogram not used !   is to be transported.  The rule for ordering tracers
! Subprogram not used !   is that a dependent tracer (such as qice) must have a larger
! Subprogram not used !   tracer index than the tracer it depends on (i.e., hice).
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_itd, only: ilyr1, slyr1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used            nx_block, ny_block, &  ! block dimensions
! Subprogram not used            ntrcr             , & ! number of tracers in use
! Subprogram not used            ntrace                ! number of tracers in use incl. hi, hs
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            aice0     ! fractional open water area
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            aicen   ,&! fractional ice area
! Subprogram not used            vicen   ,&! volume per unit area of ice          (m)
! Subprogram not used            vsnon     ! volume per unit area of snow         (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            trcrn     ! ice area tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            eicen     ! energy of melting for each ice layer (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            esnon     ! energy of melting for each snow layer (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
! Subprogram not used             intent(out)::     &
! Subprogram not used            aim       ! mean ice area in each grid cell
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace,ncat),  &
! Subprogram not used            intent(out) ::     &
! Subprogram not used            trm       ! mean tracer values in each grid cell
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
! Subprogram not used          worka, &
! Subprogram not used          workb
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            i, j, k, n   ,&! standard indices
! Subprogram not used            it, kt       ,&! tracer indices
! Subprogram not used            ij             ! combined i/j index
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::     &
! Subprogram not used            w1             ! work variable
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat) ::  &
! Subprogram not used            indxi        ,&! compressed i/j indices
! Subprogram not used            indxj
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(0:ncat) ::     &
! Subprogram not used            icells         ! number of cells with ice
! Subprogram not used 
! Subprogram not used       worka(:,:) = c0
! Subprogram not used       workb(:,:) = c0
! Subprogram not used 
! Subprogram not used       aim(:,:,0) = aice0(:,:)
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used 
! Subprogram not used          trm(:,:,:,n) = c0
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Find grid cells where ice is present and fill area array.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          icells(n) = 0
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             aim(i,j,n) = aicen(i,j,n)
! Subprogram not used             if (aim(i,j,n) > puny) then
! Subprogram not used                icells(n) = icells(n) + 1
! Subprogram not used                ij = icells(n)
! Subprogram not used                indxi(ij,n) = i
! Subprogram not used                indxj(ij,n) = j
! Subprogram not used             endif               ! aim > puny
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Fill tracer array
! Subprogram not used     ! Note: If aice > 0, then hice > 0, but we can have hsno = 0.
! Subprogram not used     ! Alse note: We transport qice*nilyr rather than qice, so as to
! Subprogram not used     !  avoid extra operations here and in tracers_to_state.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do ij = 1, icells(n)
! Subprogram not used             i = indxi(ij,n)
! Subprogram not used             j = indxj(ij,n)
! Subprogram not used             w1 = c1 / aim(i,j,n)
! Subprogram not used             worka(i,j) = c1 / vicen(i,j,n)
! Subprogram not used             trm(i,j,1,n) = vicen(i,j,n) * w1 ! hice
! Subprogram not used             trm(i,j,2,n) = vsnon(i,j,n) * w1 ! hsno
! Subprogram not used             if (trm(i,j,2,n) > puny) then
! Subprogram not used                workb(i,j) = c1 / vsnon(i,j,n)
! Subprogram not used             else
! Subprogram not used                workb(i,j) = c0
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          kt = 2
! Subprogram not used 
! Subprogram not used          do it = 1, ntrcr
! Subprogram not used             do ij = 1, icells(n)
! Subprogram not used                i = indxi(ij,n)
! Subprogram not used                j = indxj(ij,n)
! Subprogram not used                trm(i,j,kt+it,n) = trcrn(i,j,it,n) ! ice area tracers
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          kt = kt + ntrcr
! Subprogram not used 
! Subprogram not used          do k =1, nilyr
! Subprogram not used             do ij = 1, icells(n)
! Subprogram not used                i = indxi(ij,n)
! Subprogram not used                j = indxj(ij,n)
! Subprogram not used                trm(i,j,kt+k,n) = eicen(i,j,ilyr1(n)+k-1)*worka(i,j) ! qice
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! ilyr
! Subprogram not used          kt = kt + nilyr
! Subprogram not used 
! Subprogram not used          do k = 1, nslyr
! Subprogram not used             do ij = 1, icells(n)
! Subprogram not used                i = indxi(ij,n)
! Subprogram not used                j = indxj(ij,n)
! Subprogram not used                if (trm(i,j,2,n) > puny)    &    ! hsno > puny
! Subprogram not used                  trm(i,j,kt+k,n) = esnon(i,j,slyr1(n)+k-1)*workb(i,j) & ! qsno
! Subprogram not used                                  + rhos*Lfresh
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! nslyr
! Subprogram not used 
! Subprogram not used       enddo                     ! ncat
! Subprogram not used  
! Subprogram not used       end subroutine state_to_tracers

!=======================================================================
!BOP
!
! !IROUTINE: tracers_to_state - convert tracer array to state variables
!
! !INTERFACE:
!
! Subprogram not used       subroutine tracers_to_state (nx_block, ny_block,   &
! Subprogram not used                                    ntrcr,    ntrace,     &
! Subprogram not used                                    aim,      trm,        &
! Subprogram not used                                    aice0,                &
! Subprogram not used                                    aicen,    trcrn,      &
! Subprogram not used                                    vicen,    vsnon,      &
! Subprogram not used                                    eicen,    esnon) 
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Convert area and tracer arrays back to state variables.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_itd, only: ilyr1, slyr1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used            nx_block, ny_block, & ! block dimensions
! Subprogram not used            ntrcr             , & ! number of tracers in use
! Subprogram not used            ntrace                ! number of tracers in use incl. hi, hs
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            aim       ! fractional ice area
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace,ncat),  &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            trm       ! mean tracer values in each grid cell
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),     &
! Subprogram not used            intent(inout) ::     &
! Subprogram not used            aice0     ! fractional ice area
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
! Subprogram not used            intent(inout) ::     &
! Subprogram not used            aicen   ,&! fractional ice area
! Subprogram not used            vicen   ,&! volume per unit area of ice          (m)
! Subprogram not used            vsnon     ! volume per unit area of snow         (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat),  &
! Subprogram not used            intent(inout) ::     &
! Subprogram not used            trcrn     ! tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),     &
! Subprogram not used            intent(inout) ::     &
! Subprogram not used            eicen ! energy of melting for each ice layer (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),     &
! Subprogram not used            intent(inout) ::     &
! Subprogram not used            esnon ! energy of melting for each snow layer (J/m^2)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            i, j, k, n      ,&! standard indices
! Subprogram not used            it, kt          ,&! tracer indices
! Subprogram not used            icells          ,&! number of cells with ice
! Subprogram not used            ij
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block) ::     &
! Subprogram not used            indxi, indxj      ! compressed indices
! Subprogram not used 
! Subprogram not used       aice0(:,:) = aim(:,:,0)
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used 
! Subprogram not used       icells = 0
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          if (aim(i,j,n) > c0) then
! Subprogram not used             icells = icells + 1
! Subprogram not used             indxi(icells) = i
! Subprogram not used             indxj(icells) = j
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute state variables.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used             aicen(i,j,n) = aim(i,j,n)
! Subprogram not used             vicen(i,j,n) = aim(i,j,n)*trm(i,j,1,n) ! aice*hice
! Subprogram not used             vsnon(i,j,n) = aim(i,j,n)*trm(i,j,2,n) ! aice*hsno
! Subprogram not used          enddo                  ! ij
! Subprogram not used          kt = 2
! Subprogram not used 
! Subprogram not used          do it = 1, ntrcr
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used                trcrn(i,j,it,n) = trm(i,j,kt+it,n)  ! ice tracers
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! ntrcr
! Subprogram not used          kt = kt + ntrcr
! Subprogram not used 
! Subprogram not used          do k = 1, nilyr
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used                eicen(i,j,ilyr1(n)+k-1) = vicen(i,j,n)*trm(i,j,kt+k,n) 
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! nilyr
! Subprogram not used          kt = kt + nilyr
! Subprogram not used 
! Subprogram not used          do k = 1, nslyr
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used                esnon(i,j,slyr1(n)+k-1) = (trm(i,j,kt+k,n) - rhos*Lfresh) &
! Subprogram not used                                          * vsnon(i,j,n)
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! nslyr
! Subprogram not used 
! Subprogram not used       enddo                     ! ncat
! Subprogram not used 
! Subprogram not used       end subroutine tracers_to_state

!=======================================================================
!
!BOP
!
! !IROUTINE: global_conservation - check for changes in conserved quantities
!
! !INTERFACE:
!
! Subprogram not used       subroutine global_conservation (l_stop,                     &
! Subprogram not used                                       asum_init,  asum_final,     &
! Subprogram not used                                       atsum_init, atsum_final)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Check whether values of conserved quantities have changed.
! Subprogram not used ! An error probably means that ghost cells are treated incorrectly.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) ::     &
! Subprogram not used          asum_init   ,&! initial global ice area
! Subprogram not used          asum_final    ! final global ice area
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(max_ntrace), intent(in), optional :: &
! Subprogram not used          atsum_init  ,&! initial global ice area*tracer
! Subprogram not used          atsum_final   ! final global ice area*tracer
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(inout) ::     &
! Subprogram not used          l_stop    ! if true, abort on return
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            nt            ! tracer index
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::     &
! Subprogram not used            diff          ! difference between initial and final values
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       if (asum_init > puny) then
! Subprogram not used          diff = asum_final - asum_init
! Subprogram not used          if (abs(diff/asum_init) > puny) then
! Subprogram not used             l_stop = .true.
! Subprogram not used             write (nu_diag,*)
! Subprogram not used             write (nu_diag,*) 'Ice area conserv error'
! Subprogram not used             write (nu_diag,*) 'Initial global area =', asum_init
! Subprogram not used             write (nu_diag,*) 'Final global area =', asum_final
! Subprogram not used             write (nu_diag,*) 'Fractional error =', abs(diff)/asum_init
! Subprogram not used             write (nu_diag,*) 'asum_final-asum_init =', diff
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (present(atsum_init)) then
! Subprogram not used        do nt = 1, ntrace
! Subprogram not used          if (abs(atsum_init(nt)) > puny) then
! Subprogram not used             diff = atsum_final(nt) - atsum_init(nt)
! Subprogram not used             if (abs(diff/atsum_init(nt)) > puny) then
! Subprogram not used                l_stop = .true.
! Subprogram not used                write (nu_diag,*)
! Subprogram not used                write (nu_diag,*) 'area*tracer conserv error'
! Subprogram not used                write (nu_diag,*) 'tracer index =', nt
! Subprogram not used                write (nu_diag,*) 'Initial global area*tracer =',   &
! Subprogram not used                                   atsum_init(nt)
! Subprogram not used                write (nu_diag,*) 'Final global area*tracer =',     &
! Subprogram not used                                   atsum_final(nt)
! Subprogram not used                write (nu_diag,*) 'Fractional error =',             &
! Subprogram not used                                   abs(diff)/atsum_init(nt)
! Subprogram not used                write (nu_diag,*) 'atsum_final-atsum_init =', diff
! Subprogram not used             endif
! Subprogram not used          endif
! Subprogram not used        enddo
! Subprogram not used       endif                     ! present(atsum_init)
! Subprogram not used 
! Subprogram not used       end subroutine global_conservation

!=======================================================================
!BOP
!
! !IROUTINE: local_max_min - compute local max and min of a scalar field
!
! !INTERFACE:
!
! Subprogram not used       subroutine local_max_min (nx_block, ny_block,     &
! Subprogram not used                                 ilo, ihi, jlo, jhi,     &
! Subprogram not used                                 trm,                    &
! Subprogram not used                                 tmin,     tmax,         &
! Subprogram not used                                 aimask,   trmask)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! At each grid point, compute the local max and min of a scalar
! Subprogram not used ! field phi: i.e., the max and min values in the nine-cell region
! Subprogram not used ! consisting of the home cell and its eight neighbors.
! Subprogram not used ! 
! Subprogram not used ! To extend to the neighbors of the neighbors (25 cells in all),
! Subprogram not used ! follow this call with a call to quasilocal_max_min.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used            nx_block, ny_block,&! block dimensions
! Subprogram not used            ilo,ihi,jlo,jhi     ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),        &
! Subprogram not used            dimension(nx_block,ny_block) ::     &
! Subprogram not used            aimask         ! ice area mask
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),               &
! Subprogram not used            dimension (nx_block,ny_block,max_ntrace) ::    &
! Subprogram not used            trm          ,&! tracer fields
! Subprogram not used            trmask         ! tracer mask
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(out),              &
! Subprogram not used            dimension (nx_block,ny_block,max_ntrace) ::    &
! Subprogram not used            tmin         ,&! local min tracer
! Subprogram not used            tmax           ! local max tracer
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            i, j         ,&! horizontal indices
! Subprogram not used            nt, nt1        ! tracer indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block) ::     &
! Subprogram not used            phimask        ! aimask or trmask, as appropriate
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::     &
! Subprogram not used            phi_nw, phi_n, phi_ne ,&! field values in 8 neighbor cells
! Subprogram not used            phi_w, phi_e          ,&
! Subprogram not used            phi_sw, phi_s, phi_se
! Subprogram not used 
! Subprogram not used       do nt = 1, ntrace
! Subprogram not used 
! Subprogram not used          if (tracer_type(nt)==1) then  ! does not depend on another tracer
! Subprogram not used 
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                phimask(i,j) = aimask(i,j)
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          else   ! depends on another tracer
! Subprogram not used 
! Subprogram not used             nt1 = depend(nt)
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                phimask(i,j) = trmask(i,j,nt1)
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !  Store values of trm in the 8 neighbor cells.
! Subprogram not used !  If aimask = 1, use the true value; otherwise use the home cell value
! Subprogram not used !  so that non-physical values of phi do not contribute to the gradient.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used             do i = ilo, ihi
! Subprogram not used 
! Subprogram not used                phi_nw = phimask(i-1,j+1) * trm(i-1,j+1,nt)     &
! Subprogram not used                   + (c1-phimask(i-1,j+1))* trm(i,  j,  nt)
! Subprogram not used                phi_n  = phimask(i,  j+1) * trm(i,  j+1,nt)     &
! Subprogram not used                   + (c1-phimask(i,  j+1))* trm(i,  j,  nt)
! Subprogram not used                phi_ne = phimask(i+1,j+1) * trm(i+1,j+1,nt)     &
! Subprogram not used                   + (c1-phimask(i+1,j+1))* trm(i,  j,  nt)
! Subprogram not used                phi_w  = phimask(i-1,j)   * trm(i-1,j,  nt)     &
! Subprogram not used                   + (c1-phimask(i-1,j))  * trm(i,  j,  nt)
! Subprogram not used                phi_e  = phimask(i+1,j)   * trm(i+1,j,  nt)     &
! Subprogram not used                   + (c1-phimask(i+1,j))  * trm(i,  j,  nt)
! Subprogram not used                phi_sw = phimask(i-1,j-1) * trm(i-1,j-1,nt)     &
! Subprogram not used                   + (c1-phimask(i-1,j-1))* trm(i,  j,  nt)
! Subprogram not used                phi_s  = phimask(i,  j-1) * trm(i,  j-1,nt)     &
! Subprogram not used                   + (c1-phimask(i,  j-1))* trm(i,  j,  nt)
! Subprogram not used                phi_se = phimask(i+1,j-1) * trm(i+1,j-1,nt)     &
! Subprogram not used                   + (c1-phimask(i+1,j-1))* trm(i,  j,  nt)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     Compute the minimum and maximum among the nine local cells.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used                tmax(i,j,nt) = max (phi_nw, phi_n,  phi_ne, phi_w,     &
! Subprogram not used                       trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)
! Subprogram not used 
! Subprogram not used                tmin(i,j,nt) = min (phi_nw, phi_n,  phi_ne, phi_w,     &
! Subprogram not used                       trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)
! Subprogram not used 
! Subprogram not used             enddo               ! i
! Subprogram not used          enddo                  ! j
! Subprogram not used 
! Subprogram not used       enddo                     ! nt
! Subprogram not used 
! Subprogram not used       end subroutine local_max_min

!=======================================================================
!BOP
!
! !IROUTINE: quasilocal_max_min - look one grid cell farther away
!
! !INTERFACE:
!
! Subprogram not used       subroutine quasilocal_max_min (nx_block, ny_block,     &
! Subprogram not used                                      ilo, ihi, jlo, jhi,     &
! Subprogram not used                                      tmin,     tmax)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Extend the local max and min by one grid cell in each direction.
! Subprogram not used ! Incremental remapping is monotone for the "quasilocal" max and min,
! Subprogram not used ! but in rare cases may violate monotonicity for the local max and min.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used          nx_block, ny_block,&! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi     ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(inout),     &
! Subprogram not used            dimension (nx_block,ny_block,ntrace) ::     &
! Subprogram not used            tmin         ,&! local min tracer
! Subprogram not used            tmax           ! local max tracer
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            i, j          ,&! horizontal indices
! Subprogram not used            nt              ! tracer index
! Subprogram not used 
! Subprogram not used       do nt = 1, ntrace
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used 
! Subprogram not used             tmax(i,j,nt) =     &
! Subprogram not used               max (tmax(i-1,j+1,nt), tmax(i,j+1,nt), tmax(i+1,j+1,nt),     &
! Subprogram not used                    tmax(i-1,j,  nt), tmax(i,j,  nt), tmax(i+1,j,  nt),     &
! Subprogram not used                    tmax(i-1,j-1,nt), tmax(i,j-1,nt), tmax(i+1,j-1,nt))
! Subprogram not used 
! Subprogram not used             tmin(i,j,nt) =     &
! Subprogram not used               min (tmin(i-1,j+1,nt), tmin(i,j+1,nt), tmin(i+1,j+1,nt),     &
! Subprogram not used                    tmin(i-1,j,  nt), tmin(i,j,  nt), tmin(i+1,j,  nt),     &
! Subprogram not used                    tmin(i-1,j-1,nt), tmin(i,j-1,nt), tmin(i+1,j-1,nt))
! Subprogram not used 
! Subprogram not used          enddo                  ! i
! Subprogram not used          enddo                  ! j
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine quasilocal_max_min

!======================================================================
!
!BOP
!
! !IROUTINE: check_monotonicity - check bounds on new tracer values
!
! !INTERFACE:
!
! Subprogram not used       subroutine check_monotonicity (nx_block, ny_block,     &
! Subprogram not used                                      ilo, ihi, jlo, jhi,     &
! Subprogram not used                                      iblk,                   &
! Subprogram not used                                      tmin,     tmax,         &
! Subprogram not used                                      aim,      trm,          &
! Subprogram not used                                      l_stop,                 &
! Subprogram not used                                      istop,    jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! At each grid point, make sure that the new tracer values
! Subprogram not used ! fall between the local max and min values before transport.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used            nx_block, ny_block,&! block dimensions
! Subprogram not used            ilo,ihi,jlo,jhi   ,&! beginning and end of physical domain
! Subprogram not used            iblk                ! block index (diagnostic only)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),         &
! Subprogram not used            dimension (nx_block,ny_block) ::     &
! Subprogram not used            aim            ! new ice area
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),                &
! Subprogram not used            dimension (nx_block,ny_block,max_ntrace) ::     &
! Subprogram not used            trm            ! new tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),                &
! Subprogram not used            dimension (nx_block,ny_block,ntrace) ::     &
! Subprogram not used            tmin         ,&! local min tracer
! Subprogram not used            tmax           ! local max tracer
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(inout) ::     &
! Subprogram not used          l_stop    ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(inout) ::     &
! Subprogram not used          istop, jstop     ! indices of grid cell where model aborts 
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            i, j           ,&! horizontal indices
! Subprogram not used            nt, nt1, nt2     ! tracer indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::     &
! Subprogram not used            w1, w2         ! work variables
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (nx_block, ny_block) ::   &
! Subprogram not used            l_check        ! if true, check monotonicity
! Subprogram not used 
! Subprogram not used       do nt = 1, ntrace
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Load logical array to identify tracers that need checking.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (tracer_type(nt)==1) then ! does not depend on another tracer
! Subprogram not used 
! Subprogram not used             do j = jlo, jhi
! Subprogram not used             do i = ilo, ihi
! Subprogram not used                if (aim(i,j) > puny) then 
! Subprogram not used                   l_check(i,j) = .true.
! Subprogram not used                else
! Subprogram not used                   l_check(i,j) = .false.
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          elseif (tracer_type(nt)==2) then ! depends on another tracer
! Subprogram not used 
! Subprogram not used             nt1 = depend(nt)
! Subprogram not used             do j = jlo, jhi
! Subprogram not used             do i = ilo, ihi
! Subprogram not used                if (abs(trm(i,j,nt1)) > puny .and. aim(i,j) > puny) then
! Subprogram not used                   l_check(i,j) = .true.
! Subprogram not used                else
! Subprogram not used                   l_check(i,j) = .false.
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          elseif (tracer_type(nt)==3) then ! depends on two tracers
! Subprogram not used 
! Subprogram not used             nt1 = depend(nt)
! Subprogram not used             nt2 = depend(nt1)
! Subprogram not used             do j = jlo, jhi
! Subprogram not used             do i = ilo, ihi
! Subprogram not used                if (abs(trm(i,j,nt1)) > puny .and.     &
! Subprogram not used                    abs(trm(i,j,nt2)) > puny .and. aim(i,j) > puny) then
! Subprogram not used                   l_check(i,j) = .true.
! Subprogram not used                else
! Subprogram not used                   l_check(i,j) = .false.
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Make sure new values lie between tmin and tmax
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used 
! Subprogram not used             if (l_check(i,j)) then
! Subprogram not used                ! w1 and w2 allow for roundoff error when abs(trm) is big
! Subprogram not used                w1 = max(c1, abs(tmin(i,j,nt)))
! Subprogram not used                w2 = max(c1, abs(tmax(i,j,nt)))
! Subprogram not used                if (trm(i,j,nt) < tmin(i,j,nt)-w1*puny) then
! Subprogram not used                   l_stop = .true.
! Subprogram not used                   istop = i
! Subprogram not used                   jstop = j
! Subprogram not used                   write (nu_diag,*) ' '
! Subprogram not used                   write (nu_diag,*) 'new tracer < tmin'
! Subprogram not used                   write (nu_diag,*) 'i, j, nt =', i, j, nt
! Subprogram not used                   write (nu_diag,*) 'new tracer =', trm (i,j,nt)
! Subprogram not used                   write (nu_diag,*) 'tmin ='      , tmin(i,j,nt)
! Subprogram not used                   write (nu_diag,*) 'ice area ='  , aim(i,j)
! Subprogram not used                elseif (trm(i,j,nt) > tmax(i,j,nt)+w2*puny) then
! Subprogram not used                   l_stop = .true.
! Subprogram not used                   istop = i
! Subprogram not used                   jstop = j
! Subprogram not used                   write (nu_diag,*) ' '
! Subprogram not used                   write (nu_diag,*) 'new tracer > tmax'
! Subprogram not used                   write (nu_diag,*) 'i, j, nt =', i, j, nt
! Subprogram not used                   write (nu_diag,*) 'new tracer =', trm (i,j,nt)
! Subprogram not used                   write (nu_diag,*) 'tmax ='      , tmax(i,j,nt)
! Subprogram not used                   write (nu_diag,*) 'ice area ='  , aim(i,j)
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          enddo                  ! i
! Subprogram not used          enddo                  ! j
! Subprogram not used 
! Subprogram not used       enddo                     ! nt
! Subprogram not used 
! Subprogram not used       end subroutine check_monotonicity

!=======================================================================
! The remaining subroutines are called by transport_upwind.
!=======================================================================
!BOP
!
! !IROUTINE: state_to_work - fill work arrays with state variables
!
! !INTERFACE:
!
! Subprogram not used       subroutine state_to_work (nx_block, ny_block,        &
! Subprogram not used                                 ntrcr,                     &
! Subprogram not used                                 narr,     trcr_depend,     &
! Subprogram not used                                 aicen,    trcrn,           &
! Subprogram not used                                 vicen,    vsnon,           &
! Subprogram not used                                 aice0,    works)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Fill work array with state variables in preparation for upwind transport
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! same as module
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_itd, only: ilyr1, slyr1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used          nx_block, ny_block ,&! block dimensions
! Subprogram not used          ntrcr             , & ! number of tracers in use
! Subprogram not used          narr        ! number of 2D state variable arrays in works array
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (max_ntrcr), intent(in) ::     &
! Subprogram not used          trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
! Subprogram not used          intent(in) ::     &
! Subprogram not used          aicen   ,&! concentration of ice
! Subprogram not used          vicen   ,&! volume per unit area of ice          (m)
! Subprogram not used          vsnon     ! volume per unit area of snow         (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat),     &
! Subprogram not used          intent(in) ::     &
! Subprogram not used          trcrn     ! ice tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),         &
! Subprogram not used          intent(in) ::        &
! Subprogram not used          aice0     ! concentration of open water
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,narr),     &
! Subprogram not used          intent (out) ::      &
! Subprogram not used          works     ! work array
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::      &
! Subprogram not used          i, j, k, n, it ,&! counting indices
! Subprogram not used          narrays          ! counter for number of state variable arrays
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! This array is used for performance (balance memory/cache vs
! Subprogram not used       ! number of bound calls);  a different number of arrays may perform
! Subprogram not used       ! better depending on the machine used, number of processors, etc.
! Subprogram not used       ! --tested on SGI R2000, using 4 pes for the ice model under MPI
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          works(i,j,1) = aice0(i,j)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       narrays = 1
! Subprogram not used 
! Subprogram not used       do n=1, ncat
! Subprogram not used 
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             works(i,j,narrays+1) = aicen(i,j,n)
! Subprogram not used             works(i,j,narrays+2) = vicen(i,j,n)
! Subprogram not used             works(i,j,narrays+3) = vsnon(i,j,n)
! Subprogram not used          enddo                  ! i
! Subprogram not used          enddo                  ! j
! Subprogram not used          narrays = narrays + 3
! Subprogram not used 
! Subprogram not used          do it = 1, ntrcr
! Subprogram not used             if (trcr_depend(it) == 0) then
! Subprogram not used                do j = 1, ny_block
! Subprogram not used                do i = 1, nx_block
! Subprogram not used                   works(i,j,narrays+it) = aicen(i,j,n)*trcrn(i,j,it,n)
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used             elseif (trcr_depend(it) == 1) then
! Subprogram not used                do j = 1, ny_block
! Subprogram not used                do i = 1, nx_block
! Subprogram not used                   works(i,j,narrays+it) = vicen(i,j,n)*trcrn(i,j,it,n)
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used             elseif (trcr_depend(it) == 2) then
! Subprogram not used                do j = 1, ny_block
! Subprogram not used                do i = 1, nx_block
! Subprogram not used                   works(i,j,narrays+it) = vsnon(i,j,n)*trcrn(i,j,it,n)
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          narrays = narrays + ntrcr
! Subprogram not used 
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       if (narr /= narrays) write(nu_diag,*)      &
! Subprogram not used            "Wrong number of arrays in transport bound call"
! Subprogram not used 
! Subprogram not used       end subroutine state_to_work

!=======================================================================
!BOP
!
! !IROUTINE: work_to_state - convert work arrays back to state variables
!
! !INTERFACE:
!
! Subprogram not used       subroutine work_to_state (nx_block, ny_block,        &
! Subprogram not used                                 ntrcr,                     &
! Subprogram not used                                 narr,     trcr_depend,     &
! Subprogram not used                                 aicen,    trcrn,           &
! Subprogram not used                                 vicen,    vsnon,           &
! Subprogram not used                                 aice0,    works)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Convert work array back to state variables
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! same as module
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_itd, only: ilyr1, slyr1, compute_tracers
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent (in) ::                       &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          ntrcr             , & ! number of tracers in use
! Subprogram not used          narr        ! number of 2D state variable arrays in works array
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (max_ntrcr), intent(in) ::     &
! Subprogram not used          trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent (in) ::                          &
! Subprogram not used          works (nx_block,ny_block,narr)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
! Subprogram not used          intent(out) ::     &
! Subprogram not used          aicen   ,&! concentration of ice
! Subprogram not used          vicen   ,&! volume per unit area of ice          (m)
! Subprogram not used          vsnon     ! volume per unit area of snow         (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
! Subprogram not used          intent(out) ::     &
! Subprogram not used          trcrn     ! ice tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),          &
! Subprogram not used          intent(out) ::     &
! Subprogram not used          aice0     ! concentration of open water
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::      &
! Subprogram not used          i, j, k, n , it,&! counting indices
! Subprogram not used          narrays        ,&! counter for number of state variable arrays
! Subprogram not used          icells           ! number of ocean/ice cells
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block) ::        &
! Subprogram not used         indxi, indxj
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block*ny_block,narr) ::      &
! Subprogram not used          work 
! Subprogram not used 
! Subprogram not used       ! for call to compute_tracers
! Subprogram not used       icells = 0
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          icells = icells + 1
! Subprogram not used          indxi(icells) = i
! Subprogram not used          indxj(icells) = j
! Subprogram not used          work (icells,:) = works(i,j,:)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do j=1,ny_block
! Subprogram not used       do i=1,nx_block
! Subprogram not used          aice0(i,j) = works(i,j,1)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       narrays = 1               ! aice0 is first array
! Subprogram not used 
! Subprogram not used       do n=1,ncat
! Subprogram not used 
! Subprogram not used          do j=1,ny_block
! Subprogram not used          do i=1,nx_block
! Subprogram not used             aicen(i,j,n) = works(i,j,narrays+1)
! Subprogram not used             vicen(i,j,n) = works(i,j,narrays+2)
! Subprogram not used             vsnon(i,j,n) = works(i,j,narrays+3)
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          narrays = narrays + 3
! Subprogram not used 
! Subprogram not used          call compute_tracers (nx_block,     ny_block,               &
! Subprogram not used                                icells,       indxi,   indxj,         &
! Subprogram not used                                ntrcr,        trcr_depend,            &
! Subprogram not used                                work (:,narrays+1:narrays+ntrcr),     &
! Subprogram not used                                aicen(:,:,n),                         &
! Subprogram not used                                vicen(:,:,n), vsnon(:,:,n),           &
! Subprogram not used                                trcrn(:,:,:,n))
! Subprogram not used 
! Subprogram not used          narrays = narrays + ntrcr
! Subprogram not used 
! Subprogram not used       enddo                     ! ncat
! Subprogram not used 
! Subprogram not used       end subroutine work_to_state

!=======================================================================
!BOP
!
! !IROUTINE: upwind_field - advection according to upwind
!
! !INTERFACE:
!
! Subprogram not used       subroutine upwind_field (nx_block, ny_block,   &
! Subprogram not used                                ilo, ihi, jlo, jhi,   &
! Subprogram not used                                dt,                   &
! Subprogram not used                                narrays,  phi,        &
! Subprogram not used                                uee,      vnn,        &
! Subprogram not used                                HTE,      HTN,        &
! Subprogram not used                                tarea)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! upwind transport algorithm
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! same as module
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent (in) ::     &
! Subprogram not used          nx_block, ny_block ,&! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi    ,&! beginning and end of physical domain
! Subprogram not used          narrays              ! number of 2D arrays to be transported
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in) ::         &
! Subprogram not used          dt                   ! time step
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,narrays), &
! Subprogram not used          intent(inout) ::                                         &
! Subprogram not used          phi                  ! scalar field
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block),         &
! Subprogram not used          intent(in)::     &
! Subprogram not used          uee, vnn             ! cell edge velocities
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
! Subprogram not used          HTE                ,&! length of east cell edge 
! Subprogram not used          HTN                ,&! length of north cell edge
! Subprogram not used          tarea                ! grid cell area
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          i, j, k, n           ! standard indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::        &
! Subprogram not used          upwind, y1, y2, a, h   ! function
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
! Subprogram not used          worka, &
! Subprogram not used          workb
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Define upwind function
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       upwind(y1,y2,a,h) = p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! upwind transport
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       worka(:,:) = c0
! Subprogram not used       workb(:,:) = c0
! Subprogram not used 
! Subprogram not used       do n = 1, narrays
! Subprogram not used 
! Subprogram not used          do j = 1, jhi
! Subprogram not used          do i = 1, ihi
! Subprogram not used             worka(i,j)=     &
! Subprogram not used                upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j),HTE(i,j))
! Subprogram not used             workb(i,j)=     &
! Subprogram not used                upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j),HTN(i,j))
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             phi(i,j,n) = phi(i,j,n) - ( worka(i,j)-worka(i-1,j)      &
! Subprogram not used                                       + workb(i,j)-workb(i,j-1) )    &
! Subprogram not used                                       / tarea(i,j)
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       enddo                     ! narrays
! Subprogram not used 
! Subprogram not used       end subroutine upwind_field

!=======================================================================

      end module ice_transport_driver

!=======================================================================
