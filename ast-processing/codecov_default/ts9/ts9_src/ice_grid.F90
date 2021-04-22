!=======================================================================
!BOP
!
! !MODULE: ice_grid - spatial grids, masks and boundary conditions
!
! !DESCRIPTION:
!
! Spatial grids, masks, and boundary conditions
!
! !REVISION HISTORY:
!  SVN:$Id: ice_grid.F90 152 2008-09-24 20:48:59Z eclare $
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!          Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb
!       init_grid split into two parts as in POP 2.0
!       Boundary update routines replaced by POP versions
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: Option to read from netcdf files (A. Keen, Met Office)
!       Grid reading routines reworked by E. Hunke for boundary values
!
! !INTERFACE:
!
      module ice_grid
!
! !USES:
!
      use ice_kinds_mod
      use ice_boundary
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_blocks
      use ice_domain_size
      use ice_domain
      use ice_fileunits
      use ice_gather_scatter
      use ice_read_write
      use ice_timers
      use ice_probability
      use ice_exit
!
!EOP
!
      implicit none
!echmod      save

      character (len=char_len_long), save :: &
         grid_format  , & ! file format ('bin'=binary or 'nc'=netcdf)
         grid_file    , & !  input file for POP grid info
         kmt_file     , & !  input file for POP grid info
         gridcpl_file , & !  input file for POP coupling grid info
         grid_type        !  current options are rectangular (default),
                          !  displaced_pole, tripole, panarctic

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         dxt    , & ! width of T-cell through the middle (m)
         dyt    , & ! height of T-cell through the middle (m)
         dxu    , & ! width of U-cell through the middle (m)
         dyu    , & ! height of U-cell through the middle (m)
         HTE    , & ! length of eastern edge of T-cell (m)
         HTN    , & ! length of northern edge of T-cell (m)
         tarea  , & ! area of T-cell (m^2)
         uarea  , & ! area of U-cell (m^2)
         tarear , & ! 1/tarea
         uarear , & ! 1/uarea
         tinyarea,& ! puny*tarea
         tarean , & ! area of NH T-cells
         tareas , & ! area of SH T-cells
         ULON   , & ! longitude of velocity pts (radians)
         ULAT   , & ! latitude of velocity pts (radians)
         TLON   , & ! longitude of temp pts (radians)
         TLAT   , & ! latitude of temp pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET     ! ANGLE converted to T-cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         cyp    , & ! 1.5*HTE - 0.5*HTE
         cxp    , & ! 1.5*HTN - 0.5*HTN
         cym    , & ! 0.5*HTE - 1.5*HTE
         cxm    , & ! 0.5*HTN - 1.5*HTN
         dxhy   , & ! 0.5*(HTE - HTE)
         dyhx       ! 0.5*(HTN - HTN)

      ! Corners of grid boxes for history output
      real (kind=dbl_kind), dimension (4,nx_block,ny_block,max_blocks), save :: &
         lont_bounds, & ! longitude of gridbox corners for T point
         latt_bounds, & ! latitude of gridbox corners for T point
         lonu_bounds, & ! longitude of gridbox corners for U point
         latu_bounds    ! latitude of gridbox corners for U point       

      ! geometric quantities used for remapping transport
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         xav  , & ! mean T-cell value of x
         yav  , & ! mean T-cell value of y
         xxav , & ! mean T-cell value of xx
         xyav , & ! mean T-cell value of xy
         yyav , & ! mean T-cell value of yy
         xxxav, & ! mean T-cell value of xxx
         xxyav, & ! mean T-cell value of xxy
         xyyav, & ! mean T-cell value of xyy
         yyyav    ! mean T-cell value of yyy

      real (kind=dbl_kind), &
         dimension (2,2,nx_block,ny_block,max_blocks), save :: &
         mne, & ! matrices used for coordinate transformations in remapping
         mnw, & ! ne = northeast corner, nw = northwest, etc.
         mse, & 
         msw

      ! masks
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         bm     , & ! cice block mask (T-cell)
         uvm        ! land/boundary mask, velocity (U-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         ocn_gridcell_frac   ! only relevant for lat-lon grids
                             ! gridcell value of [1 - (land fraction)] (T-cell)

      logical (kind=log_kind), &
         dimension (nx_block,ny_block,max_blocks), save :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         umask  , & ! land/boundary mask, velocity (U-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      ! grid dimensions for rectangular grid
      real (kind=dbl_kind), parameter ::  &
         dxrect = 30.e5_dbl_kind   ,&! uniform HTN (cm)
         dyrect = 30.e5_dbl_kind     ! uniform HTE (cm)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         rndex_global       ! global index for local subdomain (dbl)

      real (kind=dbl_kind), private, &
         dimension(nx_block,ny_block,max_blocks) :: &
         work1

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_grid1 - - distribute blocks across processors
!
! !INTERFACE:
!
      subroutine init_grid1 
!
! !DESCRIPTION:
!
! Distribute blocks across processors.  The distribution is optimized
! based on latitude and topography, contained in the ULAT and KMT arrays. 
!
! !REVISION HISTORY:
!
! authors: William Lipscomb and Phil Jones, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_work, only: work_g1, work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      integer (kind=int_kind) :: &
         i, j, iblk, &
         fid_grid, &     ! file id for netCDF grid file
         fid_kmt         ! file id for netCDF kmt file

      real (dbl_kind), allocatable, dimension(:,:) :: bStats

      integer(kind=int_kind), allocatable :: WorkPerBlock(:),blockType(:)
      real(kind=dbl_kind), allocatable :: ProbPerBlock(:)

      character (char_len) :: &
         fieldname       ! field name in netCDF file

      !-----------------------------------------------------------------
      ! Get global ULAT and KMT arrays used for block decomposition.
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole'      ) then

         if (trim(grid_format) == 'nc') then

            call ice_open_nc(grid_file,fid_grid)
            call ice_open_nc(kmt_file,fid_kmt)

            fieldname='ulat'
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,.true.)
            fieldname='kmt'
            call ice_read_global_nc(fid_kmt,1,fieldname,work_g2,.true.)

            if (my_task == master_task) then
               call ice_close_nc(fid_grid)
               call ice_close_nc(fid_kmt)
            endif

         else

            call ice_open(nu_grid,grid_file,64) ! ULAT
            call ice_open(nu_kmt, kmt_file, 32) ! KMT

            call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)  ! ULAT
            call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT

            if (my_task == master_task) then
               close (nu_grid)
               close (nu_kmt)
            endif

         endif

      elseif (trim(grid_type) == 'panarctic') then

         call ice_open(nu_grid,grid_file,64) ! ULAT, KMT

         call ice_read_global(nu_grid,1,work_g2,'ida8',.true.)  ! KMT
         call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)  ! ULAT



         if (my_task == master_task) close (nu_grid)

      else   ! rectangular grid

         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      endif

      call broadcast_array(work_g1, master_task)   ! ULAT
      call broadcast_array(work_g2, master_task)   ! KMT

      allocate(WorkPerBlock(nblocks_tot),ProbPerBlock(nblocks_tot),blockType(nblocks_tot))
      allocate(bStats(numCoeff,nblocks_tot))
      call CalcWorkPerBlock(distribution_wght, work_g2,work_g1, &
                WorkPerBlock,ProbPerBlock,blockType,bStats)

!      call abort_ice('init_grid1: after call to CalcWorkPerBlock')
!      stop 'init_grid1: after call to CalcWorkPerBlock'

      call init_domain_distribution(work_g2, work_g1,  &
                WorkPerBlock,ProbPerBlock,blockType,bStats,maxDil)  ! KMT, ULAT
!DBG      print *,'init_grid1: after call to init_domain_distribution'

      deallocate(bStats)
      deallocate(WorkPerBlock,ProbPerBlock,blockType)

      deallocate(work_g1)
      deallocate(work_g2)

      !-----------------------------------------------------------------
      ! write additional domain information
      !-----------------------------------------------------------------

      if (my_task == master_task) then
        write(nu_diag,'(a26,i6)') '  Block size:  nx_block = ',nx_block
        write(nu_diag,'(a26,i6)') '               ny_block = ',ny_block
      endif

      end subroutine init_grid1

!=======================================================================
!BOP
!
! !IROUTINE: init_grid2 - horizontal grid initialization
!
! !INTERFACE:
!
      subroutine init_grid2
!
! !DESCRIPTION:
!
! Horizontal grid initialization:
!
!     U{LAT,LONG} = true {latitude,longitude} of U points
!     HT{N,E} = cell widths on {N,E} sides of T cell
!     ANGLE = angle between local x direction and true east
!     hm = land mask (c1 for ocean points, c0 for land points)
!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
!     T-grid and ghost cell values
!     Various grid quantities needed for dynamics and transport
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: work_g1
      use ice_exit
      use ice_blocks, only: get_block, block
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         angle_0, angle_w, angle_s, angle_sw

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks):: &
         out_of_range

      type (block) :: &
         this_block           ! block information for current block
      
      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole'      ) then
         if (trim(grid_format) == 'nc') then
            call popgrid_nc     ! read POP grid lengths from nc file
         else
            call popgrid        ! read POP grid lengths directly
         endif 
      elseif (trim(grid_type) == 'panarctic') then
         call panarctic_grid    ! pan-Arctic grid
      elseif (trim(grid_type) == 'latlon') then
         call latlongrid        ! lat lon grid for sequential CCSM (CAM mode)
         return
      else
         call rectgrid          ! regular rectangular grid
      endif

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            tarea(i,j,iblk) = dxt(i,j,iblk)*dyt(i,j,iblk)
            uarea(i,j,iblk) = dxu(i,j,iblk)*dyu(i,j,iblk)
            if (tarea(i,j,iblk) > c0) then
               tarear(i,j,iblk) = c1/tarea(i,j,iblk)
            else
               tarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (uarea(i,j,iblk) > c0) then
               uarear(i,j,iblk) = c1/uarea(i,j,iblk)
            else
               uarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            dxhy(i,j,iblk) = p5*(HTE(i,j,iblk) - HTE(i-1,j,iblk))
            dyhx(i,j,iblk) = p5*(HTN(i,j,iblk) - HTN(i,j-1,iblk))
         enddo
         enddo

         do j = jlo, jhi+1
         do i = ilo, ihi+1
            cyp(i,j,iblk) = (c1p5*HTE(i,j,iblk) - p5*HTE(i-1,j,iblk))
            cxp(i,j,iblk) = (c1p5*HTN(i,j,iblk) - p5*HTN(i,j-1,iblk))
            ! match order of operations in cyp, cxp for tripole grids
            cym(i,j,iblk) = -(c1p5*HTE(i-1,j,iblk) - p5*HTE(i,j,iblk))
            cxm(i,j,iblk) = -(c1p5*HTN(i,j-1,iblk) - p5*HTN(i,j,iblk))
         enddo
         enddo

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ghost cell updates
      ! On the tripole grid, one must be careful with updates of
      !  quantities that involve a difference of cell lengths.
      ! For example, dyhx and dxhy are cell-centered vector components.
      ! Also note that on the tripole grid, cxp and cxm would swap places,
      !  as would cyp and cym.  These quantities are computed only
      !  in north and east ghost cells (above), not south and west.
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (tarea,              halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (uarea,              halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (tarear,             halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (uarear,             halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (tinyarea,           halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (dxhy,               halo_info, &
                           field_loc_center,   field_type_vector, &
                           fillValue=c1)
      call ice_HaloUpdate (dyhx,               halo_info, &
                           field_loc_center,   field_type_vector, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Calculate ANGLET to be compatible with POP ocean model
      ! First, ensure that -pi <= ANGLE <= pi
      !-----------------------------------------------------------------

      out_of_range = .false.
      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
      if (count(out_of_range) > 0) then
         call abort_ice ('ice: init_grid: ANGLE out of expected range')
      endif

      !-----------------------------------------------------------------
      ! Compute ANGLE on T-grid
      !-----------------------------------------------------------------
      ANGLET = c0

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j,&
      !$OMP                     angle_0,angle_w,angle_s,angle_sw)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            angle_0  = ANGLE(i  ,j  ,iblk) !   w----0
            angle_w  = ANGLE(i-1,j  ,iblk) !   |    |
            angle_s  = ANGLE(i,  j-1,iblk) !   |    |
            angle_sw = ANGLE(i-1,j-1,iblk) !   sw---s

            if ( angle_0 < c0 ) then
               if ( abs(angle_w - angle_0) > pi) &
                        angle_w = angle_w  - pi2
               if ( abs(angle_s - angle_0) > pi) &
                        angle_s = angle_s  - pi2
               if ( abs(angle_sw - angle_0) > pi) &
                        angle_sw = angle_sw - pi2
            endif

            ANGLET(i,j,iblk) = angle_0 * p25 + angle_w * p25 &
                             + angle_s * p25 + angle_sw* p25
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      
      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (ANGLET,           halo_info, &
                           field_loc_center, field_type_angle, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      call makemask          ! velocity mask, hemisphere masks

      call Tlatlon           ! get lat, lon on the T grid

      !----------------------------------------------------------------
      ! Corner coordinates for CF compliant history files
      !----------------------------------------------------------------

      call gridbox_corners

      !-----------------------------------------------------------------
      ! Compute global index (used for unpacking messages from coupler)
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         do j=1,ny_global
         do i=1,nx_global
            work_g1(i,j) = real((j-1)*nx_global + i,kind=dbl_kind)
         enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(rndex_global, work_g1,  &
                          master_task,  distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine init_grid2

!=======================================================================
!BOP
!
! !IROUTINE: popgrid - read and set POP displaced pole (or tripole)
!                      grid and land mask
!
! !INTERFACE:
!
! Subprogram not used       subroutine popgrid
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! POP displaced pole grid and land mask. 
! Subprogram not used ! Grid record number, field and units are: \! (1) ULAT  (radians)    \! (2) ULON  (radians)    \! (3) HTN   (cm)         \! (4) HTE   (cm)         \! (5) HUS   (cm)         \! (6) HUW   (cm)         \! (7) ANGLE (radians)   
! Subprogram not used !
! Subprogram not used ! Land mask record number and field is (1) KMT.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_work, only: work_g1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, iblk, &
! Subprogram not used          ilo,ihi,jlo,jhi      ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: diag
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       call ice_open(nu_grid,grid_file,64)
! Subprogram not used       call ice_open(nu_kmt,kmt_file,32)
! Subprogram not used 
! Subprogram not used       diag = .true.       ! write diagnostic info
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! topography
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_read(nu_kmt,1,work1,'ida4',diag, &
! Subprogram not used                     field_loc=field_loc_center, & 
! Subprogram not used                     field_type=field_type_scalar)
! Subprogram not used 
! Subprogram not used       hm(:,:,:) = c0
! Subprogram not used       bm(:,:,:) = c0
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             hm(i,j,iblk) = work1(i,j,iblk)
! Subprogram not used             if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
! Subprogram not used             bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! lat, lon, angle
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)   ! ULAT
! Subprogram not used       call gridbox_verts(work_g1,latt_bounds)       
! Subprogram not used       call scatter_global(ULAT, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       call ice_HaloExtrapolate(ULAT, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)   ! ULON
! Subprogram not used       call gridbox_verts(work_g1,lont_bounds)       
! Subprogram not used       call scatter_global(ULON, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       call ice_HaloExtrapolate(ULON, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,7,work_g1,'rda8',.true.)   ! ANGLE
! Subprogram not used       call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_angle)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! cell dimensions
! Subprogram not used       ! calculate derived quantities from global arrays to preserve 
! Subprogram not used       ! information on boundaries
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,3,work_g1,'rda8',.true.)   ! HTN
! Subprogram not used       call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTE
! Subprogram not used       call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          close (nu_grid)
! Subprogram not used          close (nu_kmt)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       end subroutine popgrid

!=======================================================================
!BOP
!
! !IROUTINE: popgrid_nc - read and set POP tripole
!                      grid and land mask from netCDF file 
!
! !INTERFACE:
!
! Subprogram not used       subroutine popgrid_nc
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! POP displaced pole grid and land mask.
! Subprogram not used ! Grid record number, field and units are: \! (1) ULAT  (radians)    \! (2) ULON  (radians)    \! (3) HTN   (cm)         \! (4) HTE   (cm)         \! (5) HUS   (cm)         \! (6) HUW   (cm)         \! (7) ANGLE (radians)
! Subprogram not used !
! Subprogram not used ! Land mask record number and field is (1) KMT.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used ! Revised for netcdf input: Ann Keen, Met Office, May 2007
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_work, only: work_g1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, iblk, &
! Subprogram not used          ilo,ihi,jlo,jhi, &     ! beginning and end of physical domain
! Subprogram not used 	 fid_grid, &		! file id for netCDF grid file
! Subprogram not used 	 fid_kmt		! file id for netCDF kmt file
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: diag
! Subprogram not used 
! Subprogram not used       character (char_len) :: &
! Subprogram not used          fieldname		! field name in netCDF file
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       call ice_open_nc(grid_file,fid_grid)
! Subprogram not used       call ice_open_nc(kmt_file,fid_kmt)
! Subprogram not used 
! Subprogram not used       diag = .true.       ! write diagnostic info
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! topography
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       fieldname='kmt'
! Subprogram not used       call ice_read_nc(fid_kmt,1,fieldname,work1,diag, &
! Subprogram not used                        field_loc=field_loc_center, & 
! Subprogram not used                        field_type=field_type_scalar)
! Subprogram not used 
! Subprogram not used       hm(:,:,:) = c0
! Subprogram not used       bm(:,:,:) = c0
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             hm(i,j,iblk) = work1(i,j,iblk)
! Subprogram not used             if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
! Subprogram not used             bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! lat, lon, angle
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       fieldname='ulat'
! Subprogram not used       call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULAT
! Subprogram not used       call gridbox_verts(work_g1,latt_bounds)       
! Subprogram not used       call scatter_global(ULAT, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       call ice_HaloExtrapolate(ULAT, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used 
! Subprogram not used       fieldname='ulon'
! Subprogram not used       call ice_read_global_nc(fid_grid,2,fieldname,work_g1,diag) ! ULON
! Subprogram not used       call gridbox_verts(work_g1,lont_bounds)       
! Subprogram not used       call scatter_global(ULON, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       call ice_HaloExtrapolate(ULON, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used 
! Subprogram not used       fieldname='angle'
! Subprogram not used       call ice_read_global_nc(fid_grid,7,fieldname,work_g1,diag) ! ANGLE    
! Subprogram not used       call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_angle)
! Subprogram not used 
! Subprogram not used       ! fix ANGLE: roundoff error due to single precision
! Subprogram not used       where (ANGLE >  pi) ANGLE =  pi
! Subprogram not used       where (ANGLE < -pi) ANGLE = -pi
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! cell dimensions
! Subprogram not used       ! calculate derived quantities from global arrays to preserve 
! Subprogram not used       ! information on boundaries
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       fieldname='htn'
! Subprogram not used       call ice_read_global_nc(fid_grid,3,fieldname,work_g1,diag) ! HTN
! Subprogram not used       call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt
! Subprogram not used 
! Subprogram not used       fieldname='hte'
! Subprogram not used       call ice_read_global_nc(fid_grid,4,fieldname,work_g1,diag) ! HTE
! Subprogram not used       call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          call ice_close_nc(fid_grid)
! Subprogram not used          call ice_close_nc(fid_kmt)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       end subroutine popgrid_nc

!=======================================================================
!BOP
!
! !IROUTINE: panarctic_grid - read and set Pan-Arctic grid and land mask
!
! !INTERFACE:
!
! Subprogram not used       subroutine panarctic_grid
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Pan-Arctic grid and mask developed by Wieslaw Maslowski
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors: Wieslaw Maslowki, Naval Postgraduate School (based on popgrid)
! Subprogram not used !          William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_domain_size
! Subprogram not used       use ice_work, only: work_g1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! 
! Subprogram not used       ! PIPS rotated spherical grid and land mask
! Subprogram not used       !      rec no.         field         units
! Subprogram not used       !      -------         -----         -----
! Subprogram not used       !   land mask
! Subprogram not used       !         1             KMT         
! Subprogram not used       !   grid
! Subprogram not used       !         2            ULAT         radians
! Subprogram not used       !         3            ULON         radians
! Subprogram not used       !         4             HTN           cm
! Subprogram not used       !         5             HTE           cm
! Subprogram not used       !         6             HUS           cm
! Subprogram not used       !         7             HUW           cm
! Subprogram not used       !         8            ANGLE        radians
! Subprogram not used       !
! Subprogram not used       ! NOTE: There is no separate kmt file.  Land mask is part of grid file.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, iblk, &
! Subprogram not used          ilo,ihi,jlo,jhi      ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: diag
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       call ice_open(nu_grid,grid_file,64)
! Subprogram not used 
! Subprogram not used       diag = .true.       ! write diagnostic info
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) &
! Subprogram not used            write (nu_diag,*) '** Reading pan-Arctic grid **'
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! topography
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_read(nu_grid,1,work1,'ida8',diag, &
! Subprogram not used                     field_loc=field_loc_center, & 
! Subprogram not used                     field_type=field_type_scalar)
! Subprogram not used 
! Subprogram not used       hm(:,:,:) = c0
! Subprogram not used       bm(:,:,:) = c0
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             hm(i,j,iblk) = work1(i,j,iblk)
! Subprogram not used             if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
! Subprogram not used             bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! lat, lon, angle
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)   ! ULAT
! Subprogram not used       call gridbox_verts(work_g1,latt_bounds)       
! Subprogram not used       call scatter_global(ULAT, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       call ice_HaloExtrapolate(ULAT, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,3,work_g1,'rda8',.true.)   ! ULON
! Subprogram not used       call gridbox_verts(work_g1,lont_bounds)       
! Subprogram not used       call scatter_global(ULON, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       call ice_HaloExtrapolate(ULON, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,8,work_g1,'rda8',.true.)   ! ANGLE
! Subprogram not used       call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_angle)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! cell dimensions
! Subprogram not used       ! calculate derived quantities from global arrays to preserve 
! Subprogram not used       ! information on boundaries
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTN
! Subprogram not used       call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt
! Subprogram not used 
! Subprogram not used       call ice_read_global(nu_grid,5,work_g1,'rda8',.true.)   ! HTE
! Subprogram not used       call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) close (nu_grid)
! Subprogram not used 
! Subprogram not used       end subroutine panarctic_grid
!=======================================================================
!BOP
!
! !IROUTINE: latlongrid- lat and lon grid for coupling to standalone CAM
!
! !INTERFACE:
!
      subroutine latlongrid
!
! !DESCRIPTION:
!
! Read in kmt file that matches CAM lat-lon grid and has single column 
! functionality
!
! !REVISION HISTORY:
!
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90 and cice 1 calls
!
! !USES:
!
!     use ice_boundary
      use ice_domain_size
      use ice_scam, only : scmlat, scmlon, single_column
      use netcdf
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk    
      
      integer (kind=int_kind) :: &
         ni, nj, ncid, dimid, varid, ier

      character (len=char_len) :: &
         subname='latlongrid' ! subroutine name

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           closelat, &        ! Single-column latitude value
           closelon, &        ! Single-column longitude value
           closelatidx, &     ! Single-column latitude index to retrieve
           closelonidx        ! Single-column longitude index to retrieve

      integer (kind=int_kind) :: &
           start(2), &        ! Start index to read in
           count(2)           ! Number of points to read in

      integer (kind=int_kind) :: &
           start3(3), &        ! Start index to read in
           count3(3)           ! Number of points to read in

      integer (kind=int_kind) :: &
        status                ! status flag

      real (kind=dbl_kind), allocatable :: &
           lats(:),lons(:),pos_lons(:), glob_grid(:,:)  ! temporaries 

      real (kind=dbl_kind) :: &
         pos_scmlon,&         ! temporary
         scamdata             ! temporary

      !-----------------------------------------------------------------
      ! - kmt file is actually clm fractional land file
      ! - Determine consistency of dimensions
      ! - Read in lon/lat centers in degrees from kmt file
      ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
      !-----------------------------------------------------------------

      ! Determine dimension of domain file and check for consistency

      if (my_task == master_task) then
         call ice_open_nc(kmt_file, ncid)

         status = nf90_inq_dimid (ncid, 'ni', dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=ni)
         status = nf90_inq_dimid (ncid, 'nj', dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=nj)
      end if

      ! Determine start/count to read in for either single column or global lat-lon grid
      ! If single_column, then assume that only master_task is used since there is only one task

      if (single_column) then
         ! Check for consistency
         if (my_task == master_task) then
            if ((nx_global /= 1).or. (ny_global /= 1)) then
               write(nu_diag,*) 'Because you have selected the column model flag'
               write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
               write(nu_diag,*) 'ice_domain_size.F and recompile'
               call abort_ice ('latlongrid: check nx_global, ny_global')
            endif
         end if

         ! Read in domain file for single column
         allocate(lats(nj))
         allocate(lons(ni))
         allocate(pos_lons(ni))
         allocate(glob_grid(ni,nj))

         start3=(/1,1,1/)
         count3=(/ni,nj,1/)
         call check_ret(nf90_inq_varid(ncid, 'xc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, glob_grid, start3, count3), subname)
         do i = 1,ni
            lons(i) = glob_grid(i,1)
         end do

         call check_ret(nf90_inq_varid(ncid, 'yc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, glob_grid, start3, count3), subname)
         do j = 1,nj
            lats(j) = glob_grid(1,j) 
         end do
         
         ! convert lons array and scmlon to 0,360 and find index of value closest to 0
         ! and obtain single-column longitude/latitude indices to retrieve
         
         pos_lons(:)= mod(lons(:) + 360._dbl_kind,360._dbl_kind)
         pos_scmlon = mod(scmlon  + 360._dbl_kind,360._dbl_kind)
         start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
         start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

         deallocate(lats)
         deallocate(lons)
         deallocate(pos_lons)
         deallocate(glob_grid)

         call check_ret(nf90_inq_varid(ncid, 'xc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         TLON = scamdata
         call check_ret(nf90_inq_varid(ncid, 'yc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         TLAT = scamdata
         call check_ret(nf90_inq_varid(ncid, 'area' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         tarea = scamdata
         call check_ret(nf90_inq_varid(ncid, 'mask' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         hm = scamdata
         call check_ret(nf90_inq_varid(ncid, 'frac' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         ocn_gridcell_frac = scamdata
      else
         ! Check for consistency
         if (my_task == master_task) then
            if (nx_global /= ni .and. ny_global /= nj) then
              call abort_ice ('latlongrid: ni,nj not equal to nx_global,ny_global')
            end if
         end if

         ! Read in domain file for global lat-lon grid
         call ice_read_nc(ncid, 1, 'xc'  , TLON             , diag=.true.)
         call ice_read_nc(ncid, 1, 'yc'  , TLAT             , diag=.true.)
         call ice_read_nc(ncid, 1, 'area', tarea            , diag=.true., &
            field_loc=field_loc_center,field_type=field_type_scalar)
         call ice_read_nc(ncid, 1, 'mask', hm               , diag=.true.)
         call ice_read_nc(ncid, 1, 'frac', ocn_gridcell_frac, diag=.true.)
      end if

      if (my_task == master_task) then
         call ice_close_nc(ncid)
      end if

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1,nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            ! Convert from degrees to radians
            TLON(i,j,iblk) = pi*TLON(i,j,iblk)/180._dbl_kind

            ! Convert from degrees to radians
            TLAT(i,j,iblk) = pi*TLAT(i,j,iblk)/180._dbl_kind

            ! Convert from radians^2 to m^2
            ! (area in domain file is in radians^2 and tarea is in m^2)
            tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
            bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
         end do
         end do
      end do
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      ! The U grid (velocity) is not used when run with sequential CAM
      ! because we only use thermodynamic sea ice.  However, ULAT is used
      ! in the default initialization of CICE so we calculate it here as 
      ! a "dummy" so that CICE will initialize with ice.  If a no ice
      ! initialization is OK (or desired) this can be commented out and
      ! ULAT will remain 0 as specified above.  ULAT is located at the
      ! NE corner of the grid cell, TLAT at the center, so here ULAT is
      ! hacked by adding half the latitudinal spacing (in radians) to
      ! TLAT.
      !-----------------------------------------------------------------

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            if (ny_global == 1) then
               uarea(i,j,iblk)  = tarea(i,j,  iblk)
            else
               uarea(i,j,iblk)  = p25*  &
                                 (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                                + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
            endif
            tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
            uarear(i,j,iblk)   = c1/uarea(i,j,iblk)
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            if (single_column) then
               ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/nj)  
            else
               if (ny_global == 1) then
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)
               else
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)  
               endif
            endif
            ULON  (i,j,iblk) = c0
            ANGLE (i,j,iblk) = c0                             

            ANGLET(i,j,iblk) = c0                             
            HTN   (i,j,iblk) = 1.e36_dbl_kind
            HTE   (i,j,iblk) = 1.e36_dbl_kind
            dxt   (i,j,iblk) = 1.e36_dbl_kind
            dyt   (i,j,iblk) = 1.e36_dbl_kind
            dxu   (i,j,iblk) = 1.e36_dbl_kind
            dyu   (i,j,iblk) = 1.e36_dbl_kind
            dxhy  (i,j,iblk) = 1.e36_dbl_kind
            dyhx  (i,j,iblk) = 1.e36_dbl_kind
            cyp   (i,j,iblk) = 1.e36_dbl_kind
            cxp   (i,j,iblk) = 1.e36_dbl_kind
            cym   (i,j,iblk) = 1.e36_dbl_kind
            cxm   (i,j,iblk) = 1.e36_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call makemask

      end subroutine latlongrid

!=======================================================================
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
! Subprogram not used       subroutine check_ret(ret, calling)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Check return status from netcdf call
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used         use netcdf
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used         implicit none
! Subprogram not used         integer, intent(in) :: ret
! Subprogram not used         character(len=*) :: calling
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! author: Mariana Vertenstein
! Subprogram not used ! 2007: Elizabeth Hunke upgraded to netcdf90
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       if (ret /= NF90_NOERR) then
! Subprogram not used          write(nu_diag,*)'netcdf error from ',trim(calling)
! Subprogram not used          write(nu_diag,*)'netcdf strerror = ',trim(NF90_STRERROR(ret))
! Subprogram not used          call abort_ice('ice ice_grid: netcdf check_ret error')
! Subprogram not used       end if
! Subprogram not used         
! Subprogram not used       end subroutine check_ret

!=======================================================================
!BOP
!
! !IROUTINE: rectgrid - regular rectangular grid and mask
!
! !INTERFACE:
!
! Subprogram not used       subroutine rectgrid
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Regular rectangular grid and mask
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_domain_size
! Subprogram not used       use ice_work, only: work_g1
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, iblk, &
! Subprogram not used          imid, jmid
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: length
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Calculate various geometric 2d arrays
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             ANGLE(i,j,iblk) = c0              ! "square with the world"
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       ! Weddell Sea
! Subprogram not used       ! lower left corner of grid is 55W, 75S
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          work_g1 = c0
! Subprogram not used          length = dxrect*cm_to_m/radius*rad_to_deg
! Subprogram not used          work_g1(1,:) = -55._dbl_kind
! Subprogram not used          do j = 1, ny_global
! Subprogram not used          do i = 2, nx_global
! Subprogram not used             work_g1(i,j) = work_g1(i-1,j) + length   ! ULON
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          work_g1(:,:) = work_g1(:,:) / rad_to_deg
! Subprogram not used       endif
! Subprogram not used       call scatter_global(ULON, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          work_g1 = c0
! Subprogram not used          length = dyrect*cm_to_m/radius*rad_to_deg
! Subprogram not used          work_g1(:,1) = -75._dbl_kind
! Subprogram not used          do i = 1, nx_global
! Subprogram not used          do j = 2, ny_global
! Subprogram not used             work_g1(i,j) = work_g1(i,j-1) + length   ! ULAT
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          work_g1(:,:) = work_g1(:,:) / rad_to_deg
! Subprogram not used       endif
! Subprogram not used       call scatter_global(ULAT, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 1, ny_global
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g1(i,j) = dxrect             ! HTN
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call primary_grid_lengths_HTN(work_g1)  ! dxu, dxt
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 1, ny_global
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g1(i,j) = dyrect             ! HTE
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call primary_grid_lengths_HTE(work_g1)  ! dyu, dyt
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Construct T-cell land mask
! Subprogram not used       ! Keyed on ew_boundary_type; ns_boundary_type should be 'open'.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          work_g1(:,:) = c0      ! initialize hm as land
! Subprogram not used 
! Subprogram not used          if (trim(ew_boundary_type) == 'cyclic') then
! Subprogram not used 
! Subprogram not used             do j = 3,ny_global-2      ! closed top and bottom
! Subprogram not used             do i = 1,nx_global        ! open sides
! Subprogram not used                work_g1(i,j) = c1    ! NOTE nx_global > 5
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          elseif (trim(ew_boundary_type) == 'closed') then
! Subprogram not used 
! Subprogram not used             do j = 3,ny_global-2      ! closed top and bottom
! Subprogram not used             do i = 3,nx_global-2      ! closed sides
! Subprogram not used                work_g1(i,j) = c1    ! NOTE nx_global, ny_global > 5
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          elseif (trim(ew_boundary_type) == 'open') then
! Subprogram not used 
! Subprogram not used             ! land in the upper left and lower right corners,
! Subprogram not used             ! otherwise open boundaries
! Subprogram not used             imid = aint(real(nx_global)/c2,kind=int_kind)
! Subprogram not used             jmid = aint(real(ny_global)/c2,kind=int_kind)
! Subprogram not used 
! Subprogram not used             do j = 3,ny_global-2
! Subprogram not used             do i = 3,nx_global-2
! Subprogram not used                work_g1(i,j) = c1    ! open central domain
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used             do j = 1, jmid+2
! Subprogram not used             do i = 1, imid+2
! Subprogram not used                work_g1(i,j) = c1    ! open lower left corner
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used             do j = jmid-2, ny_global
! Subprogram not used             do i = imid-2, nx_global
! Subprogram not used                work_g1(i,j) = c1    ! open upper right corner
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call scatter_global(hm, work_g1, master_task, distrb_info, &
! Subprogram not used                           field_loc_center, field_type_scalar)
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       end subroutine rectgrid

!=======================================================================
!BOP
!
! !IROUTINE: primary_grid_lengths_HTN
!
! !INTERFACE:
!
! Subprogram not used       subroutine primary_grid_lengths_HTN(work_g)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Calculate dxu and dxt from HTN on the global grid, to preserve
! Subprogram not used ! ghost cell and/or land values that might otherwise be lost. Scatter
! Subprogram not used ! dxu, dxt and HTN to all processors.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_work, only: work_g2
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used ! work_g is the global array holding HTN.
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(:,:) :: work_g
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, &
! Subprogram not used          ip1     ! i+1
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g2(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g2(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used       do j = 1, ny_global
! Subprogram not used       do i = 1, nx_global
! Subprogram not used          work_g(i,j) = work_g(i,j) * cm_to_m                ! HTN
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       do j = 1, ny_global
! Subprogram not used       do i = 1, nx_global
! Subprogram not used          ! assume cyclic; noncyclic will be handled during scatter
! Subprogram not used          ip1 = i+1
! Subprogram not used          if (i == nx_global) ip1 = 1
! Subprogram not used          work_g2(i,j) = p5*(work_g(i,j) + work_g(ip1,j))    ! dxu
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(HTN, work_g, master_task, distrb_info, &
! Subprogram not used                           field_loc_Nface, field_type_scalar)
! Subprogram not used       call scatter_global(dxu, work_g2, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used       do j = 2, ny_global
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j-1)) ! dxt
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       ! extrapolate to obtain dxt along j=1
! Subprogram not used       do i = 1, nx_global
! Subprogram not used          work_g2(i,1) = c2*work_g(i,2) - work_g(i,3) ! dxt
! Subprogram not used       enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(dxt, work_g2, master_task, distrb_info, &
! Subprogram not used                           field_loc_center, field_type_scalar)
! Subprogram not used 
! Subprogram not used       deallocate(work_g2)
! Subprogram not used 
! Subprogram not used       end subroutine primary_grid_lengths_HTN

!=======================================================================
!BOP
!
! !IROUTINE: primary_grid_lengths_HTE
!
! !INTERFACE:
!
! Subprogram not used       subroutine primary_grid_lengths_HTE(work_g)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Calculate dyu and dyt from HTE on the global grid, to preserve
! Subprogram not used ! ghost cell and/or land values that might otherwise be lost. Scatter
! Subprogram not used ! dyu, dyt and HTE to all processors.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_work, only: work_g2
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used ! work_g is the global array holding HTE.
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(:,:) :: work_g
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, nyg, &
! Subprogram not used          im1     ! i-1
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g2(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g2(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used       do j = 1, ny_global
! Subprogram not used       do i = 1, nx_global
! Subprogram not used          work_g(i,j) = work_g(i,j) * cm_to_m                ! HTE
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       do j = 1, ny_global-1
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j+1)) ! dyu
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       ! extrapolate to obtain dyu along j=ny_global
! Subprogram not used       ! workaround for intel compiler
! Subprogram not used       nyg = ny_global
! Subprogram not used       do i = 1, nx_global
! Subprogram not used          work_g2(i,nyg) = c2*work_g(i,nyg-1) &
! Subprogram not used                            - work_g(i,nyg-2) ! dyu
! Subprogram not used       enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(HTE, work_g, master_task, distrb_info, &
! Subprogram not used                           field_loc_Eface, field_type_scalar)
! Subprogram not used       call scatter_global(dyu, work_g2, master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used       do j = 1, ny_global
! Subprogram not used       do i = 1, nx_global
! Subprogram not used          ! assume cyclic; noncyclic will be handled during scatter
! Subprogram not used          im1 = i-1
! Subprogram not used          if (i == 1) im1 = nx_global 
! Subprogram not used          work_g2(i,j) = p5*(work_g(i,j) + work_g(im1,j))    ! dyt
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(dyt, work_g2, master_task, distrb_info, &
! Subprogram not used                           field_loc_center, field_type_scalar)
! Subprogram not used 
! Subprogram not used       deallocate(work_g2)
! Subprogram not used 
! Subprogram not used       end subroutine primary_grid_lengths_HTE

!=======================================================================
!BOP
!
! !IROUTINE: makemask - makes logical land masks (T,U) and hemispheric masks
!
! !INTERFACE:
!
      subroutine makemask
!
! !DESCRIPTION:
!
! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (hm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_HaloUpdate (bm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            if (ny_global == 1) then
               uvm(i,j,iblk) =      hm(i,j,  iblk)
            else
               uvm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,  iblk), &
                                    hm(i,j+1,iblk), hm(i+1,j+1,iblk))
            endif
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uvm,                halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            tmask(i,j,iblk) = .false.
            umask(i,j,iblk) = .false.
            if ( hm(i,j,iblk) > p5) tmask(i,j,iblk) = .true.
            if (uvm(i,j,iblk) > p5) umask(i,j,iblk) = .true.
         enddo
         enddo

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

         lmask_n(:,:,iblk) = .false.
         lmask_s(:,:,iblk) = .false.

         tarean(:,:,iblk) = c0
         tareas(:,:,iblk) = c0

         do j = 1, ny_block
         do i = 1, nx_block

            if (ULAT(i,j,iblk) >= -puny) lmask_n(i,j,iblk) = .true. ! N. Hem.
            if (ULAT(i,j,iblk) <  -puny) lmask_s(i,j,iblk) = .true. ! S. Hem.

            ! N hemisphere area mask (m^2)
            if (lmask_n(i,j,iblk)) tarean(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

            ! S hemisphere area mask (m^2)
            if (lmask_s(i,j,iblk)) tareas(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

         enddo
         enddo

      enddo  ! iblk
      !$OMP END PARALLEL DO

      end subroutine makemask

!=======================================================================
!BOP
!
! !IROUTINE: Tlatlon - initializes latitude and longitudes on T grid
!
! !INTERFACE:
!
! Subprogram not used       subroutine Tlatlon
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Initializes latitude and longitude on T grid
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL; code originally based on POP grid
! Subprogram not used ! generation routine
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_domain_size
! Subprogram not used       use ice_global_reductions, only: global_minval, global_maxval
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       save 
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used            i, j, iblk       , & ! horizontal indices
! Subprogram not used            ig, jg           , & ! global horizontal indices
! Subprogram not used            im1              , & ! ig - 1
! Subprogram not used            ilo,ihi,jlo,jhi      ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used            z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used            this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       TLAT(:,:,:) = c0
! Subprogram not used       TLON(:,:,:) = c0
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j, &
! Subprogram not used       !$OMP                    z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4, &
! Subprogram not used       !$OMP                     tx,ty,tz,da)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used 
! Subprogram not used             z1 = cos(ULAT(i-1,j-1,iblk))
! Subprogram not used             x1 = cos(ULON(i-1,j-1,iblk))*z1
! Subprogram not used             y1 = sin(ULON(i-1,j-1,iblk))*z1
! Subprogram not used             z1 = sin(ULAT(i-1,j-1,iblk))
! Subprogram not used 
! Subprogram not used             z2 = cos(ULAT(i,j-1,iblk))
! Subprogram not used             x2 = cos(ULON(i,j-1,iblk))*z2
! Subprogram not used             y2 = sin(ULON(i,j-1,iblk))*z2
! Subprogram not used             z2 = sin(ULAT(i,j-1,iblk))
! Subprogram not used 
! Subprogram not used             z3 = cos(ULAT(i-1,j,iblk))
! Subprogram not used             x3 = cos(ULON(i-1,j,iblk))*z3
! Subprogram not used             y3 = sin(ULON(i-1,j,iblk))*z3
! Subprogram not used             z3 = sin(ULAT(i-1,j,iblk))
! Subprogram not used 
! Subprogram not used             z4 = cos(ULAT(i,j,iblk))
! Subprogram not used             x4 = cos(ULON(i,j,iblk))*z4
! Subprogram not used             y4 = sin(ULON(i,j,iblk))*z4
! Subprogram not used             z4 = sin(ULAT(i,j,iblk))
! Subprogram not used 
! Subprogram not used             tx = (x1+x2+x3+x4)/c4
! Subprogram not used             ty = (y1+y2+y3+y4)/c4
! Subprogram not used             tz = (z1+z2+z3+z4)/c4
! Subprogram not used             da = sqrt(tx**2+ty**2+tz**2)
! Subprogram not used 
! Subprogram not used             tz = tz/da
! Subprogram not used 
! Subprogram not used             ! TLON in radians East
! Subprogram not used             TLON(i,j,iblk) = c0
! Subprogram not used             if (tx /= c0 .or. ty /= c0) TLON(i,j,iblk) = atan2(ty,tx)
! Subprogram not used 
! Subprogram not used             ! TLAT in radians North
! Subprogram not used             TLAT(i,j,iblk) = asin(tz)
! Subprogram not used             
! Subprogram not used          enddo                  ! i
! Subprogram not used          enddo                  ! j         
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used       call ice_HaloUpdate (TLON,             halo_info, &
! Subprogram not used                            field_loc_center, field_type_scalar, &
! Subprogram not used                            fillValue=c1)
! Subprogram not used       call ice_HaloUpdate (TLAT,             halo_info, &
! Subprogram not used                            field_loc_center, field_type_scalar, &
! Subprogram not used                            fillValue=c1)
! Subprogram not used       call ice_HaloExtrapolate(TLON, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used       call ice_HaloExtrapolate(TLAT, distrb_info, &
! Subprogram not used                                ew_boundary_type, ns_boundary_type)
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       x1 = global_minval(TLON, distrb_info, tmask)
! Subprogram not used       x2 = global_maxval(TLON, distrb_info, tmask)
! Subprogram not used       x3 = global_minval(TLAT, distrb_info, tmask)
! Subprogram not used       x4 = global_maxval(TLAT, distrb_info, tmask)
! Subprogram not used 
! Subprogram not used       y1 = global_minval(ULON, distrb_info, umask)
! Subprogram not used       y2 = global_maxval(ULON, distrb_info, umask)
! Subprogram not used       y3 = global_minval(ULAT, distrb_info, umask)
! Subprogram not used       y4 = global_maxval(ULAT, distrb_info, umask)
! Subprogram not used 
! Subprogram not used       if (my_task==master_task) then
! Subprogram not used          write(nu_diag,*) ' '
! Subprogram not used          write(nu_diag,*) 'min/max ULON:', y1*rad_to_deg, y2*rad_to_deg
! Subprogram not used          write(nu_diag,*) 'min/max TLON:', x1*rad_to_deg, x2*rad_to_deg
! Subprogram not used          write(nu_diag,*) 'min/max ULAT:', y3*rad_to_deg, y4*rad_to_deg
! Subprogram not used          write(nu_diag,*) 'min/max TLAT:', x3*rad_to_deg, x4*rad_to_deg
! Subprogram not used       endif                     ! my_task
! Subprogram not used 
! Subprogram not used       end subroutine Tlatlon

!=======================================================================
!BOP
!
! !IROUTINE: t2ugrid_vector - transfer vector from T-cells to U-cells
!
! !INTERFACE:
!
! Subprogram not used       subroutine t2ugrid_vector (work)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Transfer vector component from T-cell centers to U-cell centers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
! Subprogram not used            intent(inout) :: & 
! Subprogram not used            work
! Subprogram not used       
! Subprogram not used       integer (int_kind) :: iblk
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       work1(:,:,:) = work(:,:,:)
! Subprogram not used !      do iblk = 1,nblocks
! Subprogram not used !         work1(:,:,iblk) = work(:,:,iblk)
! Subprogram not used !      enddo
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used       call ice_HaloUpdate (work1,            halo_info, &
! Subprogram not used                            field_loc_center, field_type_vector)
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       call to_ugrid(work1,work)
! Subprogram not used 
! Subprogram not used       end subroutine t2ugrid_vector

!=======================================================================
!BOP
!
! !IROUTINE: to_ugrid - shift from T-cell to U-cell midpoints
!
! !INTERFACE:
!
! Subprogram not used       subroutine to_ugrid(work1,work2)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Shifts quantities from the T-cell midpoint (work1) to the U-cell
! Subprogram not used ! midpoint (work2)
! Subprogram not used ! NOTE: Input array includes ghost cells that must be updated before
! Subprogram not used !       calling this routine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          work1(nx_block,ny_block,max_blocks)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(out) :: &
! Subprogram not used          work2(nx_block,ny_block,max_blocks)
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, iblk, &
! Subprogram not used          ilo,ihi,jlo,jhi      ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       work2(:,:,:) = c0
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             work2(i,j,iblk) = p25 * &
! Subprogram not used                               (work1(i,  j,  iblk)*tarea(i,  j,  iblk)  &
! Subprogram not used                              + work1(i+1,j,  iblk)*tarea(i+1,j,  iblk)  &
! Subprogram not used                              + work1(i,  j+1,iblk)*tarea(i,  j+1,iblk)  &
! Subprogram not used                              + work1(i+1,j+1,iblk)*tarea(i+1,j+1,iblk)) &
! Subprogram not used                              / uarea(i,  j,  iblk)
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       end subroutine to_ugrid

!=======================================================================
!BOP
!
! !IROUTINE: u2tgrid_vector - transfer vector from U-cells to T-cells
!
! !INTERFACE:
!
! Subprogram not used       subroutine u2tgrid_vector (work)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Transfer from U-cell centers to T-cell centers. Writes work into
! Subprogram not used ! another array that has ghost cells
! Subprogram not used ! NOTE: Input array is dimensioned only over physical cells.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          work
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       work1(:,:,:) = work(:,:,:)
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used       call ice_HaloUpdate (work1,            halo_info, &
! Subprogram not used                            field_loc_NEcorner, field_type_vector)
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       call to_tgrid(work1,work)
! Subprogram not used 
! Subprogram not used       end subroutine u2tgrid_vector

!=======================================================================
!BOP
!
! !IROUTINE: to_tgrid - shifts array from U-cell to T-cell midpoints
!
! !INTERFACE:
!
! Subprogram not used       subroutine to_tgrid(work1, work2)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Shifts quantities from the U-cell midpoint (work1) to the T-cell
! Subprogram not used ! midpoint (work2)
! Subprogram not used ! NOTE: Input array includes ghost cells that must be updated before
! Subprogram not used !       calling this routine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind) :: work1(nx_block,ny_block,max_blocks), &
! Subprogram not used                               work2(nx_block,ny_block,max_blocks)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, iblk, &
! Subprogram not used          ilo,ihi,jlo,jhi      ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used       
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             work2(i,j,iblk) = p25 *  &
! Subprogram not used                              (work1(i,  j  ,iblk) * uarea(i,  j,  iblk)  &
! Subprogram not used                             + work1(i-1,j  ,iblk) * uarea(i-1,j,  iblk)  &
! Subprogram not used                             + work1(i,  j-1,iblk) * uarea(i,  j-1,iblk)  & 
! Subprogram not used                             + work1(i-1,j-1,iblk) * uarea(i-1,j-1,iblk)) &
! Subprogram not used                             / tarea(i,  j,  iblk)
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       end subroutine to_tgrid

!=======================================================================
!BOP
!
! !IROUTINE: sinlat - calculates sin of latitudes
!
! !INTERFACE:
!
! Subprogram not used       subroutine sinlat(a, k)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Calculates the sin of latitudes on the grid based on ny_global and
! Subprogram not used ! using code from CAM (/control/gauaw_mod.F90)  In CAM, gauaw_mod.F90 is
! Subprogram not used ! only used to calculate the sin of latitudes (and thence latitude) if
! Subprogram not used ! the dynamical core is set to eul.  If using one of the other dynamical
! Subprogram not used ! cores and coupling to stand alone CAM the latitudes calculated in this
! Subprogram not used ! way may not match the grid from CAM.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Jacob Sewall
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used       use ice_exit
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(ny_global), intent(out) :: & 
! Subprogram not used            a            ! sin of latitudes
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            k            ! number of latitudes (ny_global)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), dimension(k) :: &
! Subprogram not used            sinlats      ! sine of latitudes
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used            eps      , & ! convergence criterion
! Subprogram not used            c        , & ! constant combination
! Subprogram not used            fk       , & ! real k
! Subprogram not used            xz       , & ! abscissa estimate
! Subprogram not used            pkm1     , & ! |
! Subprogram not used            pkm2     , & ! |-polynomials
! Subprogram not used            pkmrk    , & ! |
! Subprogram not used            pk       , & ! |
! Subprogram not used            sp       , & ! current iteration latitude increment
! Subprogram not used            avsp     , & ! |sp|
! Subprogram not used            fn       , & ! real n
! Subprogram not used            avsp_prev, &
! Subprogram not used            testdiff
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), parameter :: &
! Subprogram not used            eps27 = 1.e-27_dbl_kind
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used            n, l     , &
! Subprogram not used            iter     , & ! iteration counter
! Subprogram not used            is       , & ! latitude index
! Subprogram not used            kk           ! k/2 (number of latitudes in a hemisphere)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       c = (c1-(c2/pi)**2)*p25
! Subprogram not used       fk = k
! Subprogram not used       kk = k/2
! Subprogram not used       call bsslzr(sinlats,kk)
! Subprogram not used       do is=1,kk
! Subprogram not used         xz = cos(sinlats(is)/sqrt((fk+p5)**2+c))
! Subprogram not used !
! Subprogram not used ! This is the first approximation to xz
! Subprogram not used !
! Subprogram not used         iter = 0
! Subprogram not used         avsp = c100     !initialize avsp to a very large number
! Subprogram not used 10      continue
! Subprogram not used         avsp_prev = avsp
! Subprogram not used         pkm2 = c1
! Subprogram not used         pkm1 = xz
! Subprogram not used         iter = iter + 1
! Subprogram not used         if (iter > 100) then  ! Error exit
! Subprogram not used            call abort_ice('SINLAT: no convergence in 100 iterations')
! Subprogram not used         end if
! Subprogram not used !
! Subprogram not used ! Computation of the legendre polynomial
! Subprogram not used !
! Subprogram not used         do n=2,k
! Subprogram not used           fn = n
! Subprogram not used           pk = ((c2*fn-1._dbl_kind)*xz*pkm1-(fn-c1)*pkm2)/fn
! Subprogram not used           pkm2 = pkm1
! Subprogram not used           pkm1 = pk
! Subprogram not used         enddo
! Subprogram not used         pkm1 = pkm2
! Subprogram not used         pkmrk = (fk*(pkm1-xz*pk))/(c1-xz**2)
! Subprogram not used         sp = pk/pkmrk
! Subprogram not used         xz = xz - sp
! Subprogram not used         avsp = abs(sp)
! Subprogram not used         testdiff = avsp_prev - avsp
! Subprogram not used         if (testdiff > eps27) go to 10
! Subprogram not used         sinlats(is) = xz
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used ! Complete the sets of abscissas and weights, using the symmetry.
! Subprogram not used ! Also note truncation from real(r8) to real*8
! Subprogram not used !
! Subprogram not used       do n=1,kk
! Subprogram not used         l = k + 1 - n
! Subprogram not used         a(n) = sinlats(n)
! Subprogram not used         a(l) = -sinlats(n)
! Subprogram not used       end do
! Subprogram not used  
! Subprogram not used       end subroutine sinlat

!=======================================================================
!BOP
!
! !IROUTINE: bsslzr - Return n zeros (if n<50) of the Bessel function
!
! !INTERFACE:
!
! Subprogram not used       subroutine bsslzr(bes, n)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Return n zeros (or if n>50, approximate zeros), of the Bessel function
! Subprogram not used ! j0,in the array bes. The first 50 zeros will be given exactly, and the
! Subprogram not used ! remaining zeros are computed by extrapolation,and therefore not exact.
! Subprogram not used !
! Subprogram not used ! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic
! Subprogram not used ! placed in ice_grid.F 4/26/2006 by Jacob Sewall
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! Original version:  CCM1
! Subprogram not used ! Standardized:      J. Rosinski, June 1992
! Subprogram not used ! Reviewed:          J. Hack, D. Williamson, August 1992
! Subprogram not used ! Reviewed:          J. Hack, D. Williamson, April 1996
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            n          ! number of latitudes in hemisphere (ny_global/2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(n), intent(inout) :: & 
! Subprogram not used            bes        ! sin of latitudes
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used !----------------------------------------
! Subprogram not used ! Local Variables
! Subprogram not used !----------------------------------------
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          nn, j
! Subprogram not used      
! Subprogram not used       real (kind=dbl_kind), dimension(50) :: bz
! Subprogram not used       
! Subprogram not used       save bz
! Subprogram not used !
! Subprogram not used !-----------------------------------------
! Subprogram not used ! Local Workspace
! Subprogram not used !-----------------------------------------
! Subprogram not used 
! Subprogram not used       data bz/ &
! Subprogram not used         2.4048255577_dbl_kind,   5.5200781103_dbl_kind,   8.6537279129_dbl_kind, &
! Subprogram not used        11.7915344391_dbl_kind,  14.9309177086_dbl_kind,  18.0710639679_dbl_kind, &
! Subprogram not used        21.2116366299_dbl_kind,  24.3524715308_dbl_kind,  27.4934791320_dbl_kind, &
! Subprogram not used        30.6346064684_dbl_kind,  33.7758202136_dbl_kind,  36.9170983537_dbl_kind, &
! Subprogram not used        40.0584257646_dbl_kind,  43.1997917132_dbl_kind,  46.3411883717_dbl_kind, &
! Subprogram not used        49.4826098974_dbl_kind,  52.6240518411_dbl_kind,  55.7655107550_dbl_kind, &
! Subprogram not used        58.9069839261_dbl_kind,  62.0484691902_dbl_kind,  65.1899648002_dbl_kind, &
! Subprogram not used        68.3314693299_dbl_kind,  71.4729816036_dbl_kind,  74.6145006437_dbl_kind, &
! Subprogram not used        77.7560256304_dbl_kind,  80.8975558711_dbl_kind,  84.0390907769_dbl_kind, &
! Subprogram not used        87.1806298436_dbl_kind,  90.3221726372_dbl_kind,  93.4637187819_dbl_kind, &
! Subprogram not used        96.6052679510_dbl_kind,  99.7468198587_dbl_kind, 102.8883742542_dbl_kind, &
! Subprogram not used       106.0299309165_dbl_kind, 109.1714896498_dbl_kind, 112.3130502805_dbl_kind, &
! Subprogram not used       115.4546126537_dbl_kind, 118.5961766309_dbl_kind, 121.7377420880_dbl_kind, &
! Subprogram not used       124.8793089132_dbl_kind, 128.0208770059_dbl_kind, 131.1624462752_dbl_kind, &
! Subprogram not used       134.3040166383_dbl_kind, 137.4455880203_dbl_kind, 140.5871603528_dbl_kind, &
! Subprogram not used       143.7287335737_dbl_kind, 146.8703076258_dbl_kind, 150.0118824570_dbl_kind, &
! Subprogram not used       153.1534580192_dbl_kind, 156.2950342685_dbl_kind/  
! Subprogram not used 
! Subprogram not used       nn = n 
! Subprogram not used       if (n > 50) then 
! Subprogram not used          bes(50) = bz(50) 
! Subprogram not used          do j = 51, n 
! Subprogram not used             bes(j) = bes(j-1) + pi 
! Subprogram not used          end do 
! Subprogram not used          nn = 49 
! Subprogram not used       endif 
! Subprogram not used       bes(:nn) = bz(:nn) 
! Subprogram not used 
! Subprogram not used       end subroutine bsslzr 

!=======================================================================
! The following code is used for obtaining the coordinates of the grid
! vertices for CF-compliant netCDF history output. Approximate!
!=======================================================================
!
!BOP
!
! !IROUTINE: gridbox_corners - get coordinates of grid box corners
!
! !INTERFACE:
!
! Subprogram not used       subroutine gridbox_corners
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! These fields are only used for netcdf history output, and the
! Subprogram not used ! ghost cell values are not needed.
! Subprogram not used ! NOTE:  Extrapolations were used: these fields are approximate!
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors:   A. McLaren, Met Office
! Subprogram not used !            E. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used       use ice_work, only: work_g2
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used           i,j,iblk,icorner,& ! index counters
! Subprogram not used           ilo,ihi,jlo,jhi    ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       !-------------------------------------------------------------
! Subprogram not used       ! Get coordinates of grid boxes for each block as follows:
! Subprogram not used       ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
! Subprogram not used       !-------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used 
! Subprogram not used             latu_bounds(1,i,j,iblk)=TLAT(i  ,j  ,iblk)*rad_to_deg
! Subprogram not used             latu_bounds(2,i,j,iblk)=TLAT(i+1,j  ,iblk)*rad_to_deg
! Subprogram not used             latu_bounds(3,i,j,iblk)=TLAT(i+1,j+1,iblk)*rad_to_deg
! Subprogram not used             latu_bounds(4,i,j,iblk)=TLAT(i  ,j+1,iblk)*rad_to_deg         
! Subprogram not used 
! Subprogram not used             lonu_bounds(1,i,j,iblk)=TLON(i  ,j  ,iblk)*rad_to_deg
! Subprogram not used             lonu_bounds(2,i,j,iblk)=TLON(i+1,j  ,iblk)*rad_to_deg
! Subprogram not used             lonu_bounds(3,i,j,iblk)=TLON(i+1,j+1,iblk)*rad_to_deg
! Subprogram not used             lonu_bounds(4,i,j,iblk)=TLON(i  ,j+1,iblk)*rad_to_deg         
! Subprogram not used 
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       ! extrapolate on global grid to get edge values
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g2(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g2(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       work1(:,:,:) = latu_bounds(2,:,:,:)
! Subprogram not used       call gather_global(work_g2, work1, master_task, distrb_info)
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 1, ny_global
! Subprogram not used             work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
! Subprogram not used                                     - work_g2(nx_global-2,j)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       latu_bounds(2,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work1(:,:,:) = latu_bounds(3,:,:,:)
! Subprogram not used       call gather_global(work_g2, work1, master_task, distrb_info)
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
! Subprogram not used                                     - work_g2(i,ny_global-2)
! Subprogram not used          enddo
! Subprogram not used          do j = 1, ny_global
! Subprogram not used             work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
! Subprogram not used                                     - work_g2(nx_global-2,j)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       latu_bounds(3,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work1(:,:,:) = latu_bounds(4,:,:,:)
! Subprogram not used       call gather_global(work_g2, work1, master_task, distrb_info)
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
! Subprogram not used                                     - work_g2(i,ny_global-2)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       latu_bounds(4,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work1(:,:,:) = lonu_bounds(2,:,:,:)
! Subprogram not used       call gather_global(work_g2, work1, master_task, distrb_info)
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 1, ny_global
! Subprogram not used             work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
! Subprogram not used                                     - work_g2(nx_global-2,j)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       lonu_bounds(2,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work1(:,:,:) = lonu_bounds(3,:,:,:)
! Subprogram not used       call gather_global(work_g2, work1, master_task, distrb_info)
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
! Subprogram not used                                     - work_g2(i,ny_global-2)
! Subprogram not used          enddo
! Subprogram not used          do j = 1, ny_global
! Subprogram not used             work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
! Subprogram not used                                     - work_g2(nx_global-2,j)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       lonu_bounds(3,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work1(:,:,:) = lonu_bounds(4,:,:,:)
! Subprogram not used       call gather_global(work_g2, work1, master_task, distrb_info)
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
! Subprogram not used                                     - work_g2(i,ny_global-2)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       lonu_bounds(4,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       deallocate(work_g2)
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       ! Convert longitude to Degrees East >0 for history output
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       allocate(work_g2(nx_block,ny_block))  ! not used as global here
! Subprogram not used       !NOTE - this is commented below due to problems with OpenMP
! Subprogram not used       !reproducibility in this loop  
! Subprogram not used       !!$OMP PARALLEL DO PRIVATE(iblk,icorner,work_g2)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do icorner = 1, 4
! Subprogram not used             work_g2(:,:) = lont_bounds(icorner,:,:,iblk) + c360
! Subprogram not used             where (work_g2 > c360) work_g2 = work_g2 - c360
! Subprogram not used             where (work_g2 < c0 )  work_g2 = work_g2 + c360
! Subprogram not used             lont_bounds(icorner,:,:,iblk) = work_g2(:,:)
! Subprogram not used 
! Subprogram not used             work_g2(:,:) = lonu_bounds(icorner,:,:,iblk) + c360
! Subprogram not used             where (work_g2 > c360) work_g2 = work_g2 - c360
! Subprogram not used             where (work_g2 < c0 )  work_g2 = work_g2 + c360
! Subprogram not used             lonu_bounds(icorner,:,:,iblk) = work_g2(:,:)
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !!$OMP END PARALLEL DO
! Subprogram not used       deallocate(work_g2)
! Subprogram not used 
! Subprogram not used       end subroutine gridbox_corners

!=======================================================================
!
!BOP
!
! !IROUTINE: gridbox_verts - coordinates of grid box vertices
!
! !INTERFACE:
!
! Subprogram not used       subroutine gridbox_verts(work_g,vbounds)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! NOTE:  Boundary conditions for fields on NW, SW, SE corners
! Subprogram not used !        have not been implemented; using NE corner location for all.
! Subprogram not used !        Extrapolations are also used: these fields are approximate!
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors:   A. McLaren, Met Office
! Subprogram not used !            E. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used       use ice_work, only: work_g2
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), dimension(:,:), intent(in) :: work_g
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), &
! Subprogram not used           dimension(4,nx_block,ny_block,max_blocks), &
! Subprogram not used           intent(out) :: vbounds
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used           i,j,             & ! index counters
! Subprogram not used           ilo,ihi,jlo,jhi    ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g2(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g2(1,1))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !-------------------------------------------------------------
! Subprogram not used       ! Get coordinates of grid boxes for each block as follows:
! Subprogram not used       ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
! Subprogram not used       !-------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       work_g2(:,:) = c0
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 2, ny_global
! Subprogram not used          do i = 2, nx_global
! Subprogram not used             work_g2(i,j) = work_g(i-1,j-1) * rad_to_deg
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! extrapolate
! Subprogram not used          do j = 1, ny_global
! Subprogram not used             work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
! Subprogram not used          enddo
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,1) = c2*work_g2(i,2) - work_g2(i,3)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       vbounds(1,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work_g2(:,:) = c0
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 2, ny_global
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,j) = work_g(i,j-1) * rad_to_deg
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! extrapolate
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,1) = (c2*work_g2(i,2) - work_g2(i,3))
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       vbounds(2,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work_g2(:,:) = c0
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 1, ny_global
! Subprogram not used          do i = 1, nx_global
! Subprogram not used             work_g2(i,j) = work_g(i,j) * rad_to_deg
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       vbounds(3,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       work_g2(:,:) = c0
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          do j = 1, ny_global
! Subprogram not used          do i = 2, nx_global
! Subprogram not used             work_g2(i,j) = work_g(i-1,j  ) * rad_to_deg         
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          ! extrapolate
! Subprogram not used          do j = 1, ny_global
! Subprogram not used             work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used       call scatter_global(work1, work_g2, &
! Subprogram not used                           master_task, distrb_info, &
! Subprogram not used                           field_loc_NEcorner, field_type_scalar)
! Subprogram not used       vbounds(4,:,:,:) = work1(:,:,:)
! Subprogram not used 
! Subprogram not used       deallocate (work_g2)
! Subprogram not used 
! Subprogram not used       end subroutine gridbox_verts

!=======================================================================

      end module ice_grid

!=======================================================================
