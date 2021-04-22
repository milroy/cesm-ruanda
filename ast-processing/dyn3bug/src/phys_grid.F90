module phys_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of physics computational horizontal grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code.
! 
! Entry points:
!      phys_grid_init       initialize chunk'ed data structure
!      phys_grid_initialized    get physgrid_set flag
!
!      phys_grid_defaultopts   get default runtime options
!      phys_grid_setopts       set runtime options
!
!      get_chunk_indices_p get local chunk index range
!      get_ncols_p         get number of columns for a given chunk
!      get_xxx_all_p       get global indices, coordinates, or values
!                          for a given chunk
!      get_xxx_vec_p       get global indices, coordinates, or values
!                          for a subset of the columns in a chunk
!      get_xxx_p           get global indices, coordinates, or values
!                          for a single column
!      where xxx is
!       area               for column surface area (in radians squared)
!       gcol               for global column index
!       lat                for global latitude index
!       lon                for global longitude index
!       rlat               for latitude coordinate (in radians)
!       rlon               for longitude coordinate (in radians)
!       wght               for column integration weight
!
!      scatter_field_to_chunk
!                          distribute field
!                          to decomposed chunk data structure
!      gather_chunk_to_field
!                          reconstruct field
!                          from decomposed chunk data structure
!
!      read_chunk_from_field
!                          read and distribute field
!                          to decomposed chunk data structure
!      write_field_from_chunk
!                          write field
!                          from decomposed chunk data structure
!
!      block_to_chunk_send_pters
!                          return pointers into send buffer where data
!                          from decomposed fields should
!                          be copied to
!      block_to_chunk_recv_pters
!                          return pointers into receive buffer where data
!                          for decomposed chunk data structures should
!                          be copied from
!      transpose_block_to_chunk
!                          transpose buffer containing decomposed 
!                          fields to buffer
!                          containing decomposed chunk data structures
!
!      chunk_to_block_send_pters
!                          return pointers into send buffer where data
!                          from decomposed chunk data structures should
!                          be copied to
!      chunk_to_block_recv_pters
!                          return pointers into receive buffer where data
!                          for decomposed fields should
!                          be copied from
!      transpose_chunk_to_block
!                          transpose buffer containing decomposed
!                          chunk data structures to buffer
!                          containing decomposed fields
!
!      chunk_index         identify whether index is for a latitude or
!                          a chunk
!
! FOLLOWING ARE NO LONGER USED, AND ARE CURRENTLY COMMENTED OUT
!      get_gcol_owner_p    get owner of column
!                          for given global physics column index
!
!      buff_to_chunk       Copy from local buffer to local chunk data 
!                          structure. (Needed for cpl6.)
!
!      chunk_to_buff       Copy from local chunk data structure to 
!                          local buffer. (Needed for cpl6.)
!
! Author: Patrick Worley and John Drake
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use physconst,    only: pi
   use ppgrid,       only: pcols, pver, begchunk, endchunk
   use spmd_dyn,     only: block_buf_nrecs, chunk_buf_nrecs, &
                           local_dp_map
   use mpishorthand
   use spmd_utils,   only: iam, masterproc, npes, proc_smp_map, nsmps
   use m_MergeSorts, only: IndexSet, IndexSort
   use abortutils,   only: endrun
   use perf_mod
   use cam_logfile,  only: iulog

   implicit none
   save


! dynamics field grid information
   integer, private :: hdim1_d, hdim2_d
                                       ! dimensions of rectangular horizontal grid
                                       ! data structure, If 1D data structure, then
                                       ! hdim2_d == 1.

! physics field data structures
   integer         :: ngcols           ! global column count in physics grid (all)
   integer, public :: ngcols_p         ! global column count in physics grid 
                                       ! (without holes)

   integer, dimension(:), allocatable, private :: dyn_to_latlon_gcol_map
                                       ! map from unsorted (dynamics) to lat/lon sorted grid indices
   integer, dimension(:), allocatable, private :: latlon_to_dyn_gcol_map
                                       ! map from lat/lon sorted grid to unsorted (dynamics) indices
   integer, dimension(:), allocatable, private :: lonlat_to_dyn_gcol_map
                                       ! map from lon/lat sorted grid to unsorted (dynamics) indices

!   integer, private :: clat_p_tot ! number of unique latitudes
!   integer, private :: clon_p_tot ! number of unique longitudes
! these are public to support mozart chemistry in the short term
   integer, public :: clat_p_tot ! number of unique latitudes
   integer, public :: clon_p_tot ! number of unique longitudes

   integer, dimension(:), allocatable, private :: clat_p_cnt ! number of repeats for each latitude
   integer, dimension(:), allocatable, private :: clat_p_idx ! index in latlon ordering for first occurence
                                                             ! of latitude corresponding to given 
                                                             ! latitude index
   real(r8), dimension(:), allocatable :: clat_p  ! unique latitudes (radians, increasing)


   integer, dimension(:), allocatable, private :: clon_p_cnt ! number of repeats for each longitude
   real(r8), dimension(:), allocatable :: clon_p  ! unique longitudes (radians, increasing)

   integer, dimension(:), allocatable, private :: lat_p      ! index into list of unique column latitudes
   integer, dimension(:), allocatable, private :: lon_p      ! index into list of unique column longitudes

! chunk data structures
   type chunk
     integer  :: ncols                 ! number of vertical columns
     integer  :: gcol(pcols)           ! global physics column indices
     integer  :: lon(pcols)            ! global longitude indices
     integer  :: lat(pcols)            ! global latitude indices
     integer  :: owner                 ! id of process where chunk assigned
     integer  :: lcid                  ! local chunk index
   end type chunk

   integer :: nchunks                  ! global chunk count
   type (chunk), dimension(:), allocatable, public :: chunks  
                                       ! global computational grid

   integer, dimension(:), allocatable, private :: npchunks 
                                       ! number of chunks assigned to each process

   type lchunk
     integer  :: ncols                 ! number of vertical columns
     integer  :: cid                   ! global chunk index
     integer  :: gcol(pcols)           ! global physics column indices
     real(r8) :: area(pcols)           ! column surface area (from dynamics)
     real(r8) :: wght(pcols)           ! column integration weight (from dynamics)
   end type lchunk

   integer, private :: nlchunks        ! local chunk count
   type (lchunk), dimension(:), allocatable, private :: lchunks  
                                       ! local chunks

   type knuhc
     integer  :: chunkid               ! chunk id
     integer  :: col                   ! column index in chunk
   end type knuhc

   type (knuhc), dimension(:), allocatable, private :: knuhcs
                                       ! map from global column indices
                                       ! to chunk'ed grid

! column mapping data structures
   type column_map
     integer  :: chunk                 ! global chunk index
     integer  :: ccol                  ! column ordering in chunk
   end type column_map

   integer, private :: nlcols           ! local column count
   type (column_map), dimension(:), allocatable, private :: pgcols
                                       ! ordered list of columns (for use in gather/scatter)
                                       ! NOTE: consistent with local ordering

! column remap data structures
   integer, dimension(:), allocatable, private :: gs_col_num
                                       ! number of columns scattered to each process in
                                       ! field_to_chunk scatter
   integer, dimension(:), allocatable, private :: gs_col_offset
                                       ! offset of columns (-1) in pgcols scattered to
                                       ! each process in field_to_chunk scatter

   integer, dimension(:), allocatable, private :: btofc_blk_num
                                       ! number of grid points scattered to each process in
                                       ! block_to_chunk alltoallv, and gathered from each
                                       ! process in chunk_to_block alltoallv

   integer, dimension(:), allocatable, private :: btofc_chk_num
                                       ! number of grid points gathered from each process in
                                       ! block_to_chunk alltoallv, and scattered to each
                                       ! process in chunk_to_block alltoallv

   type btofc_pters
     integer :: ncols                  ! number of columns in block
     integer :: nlvls                  ! number of levels in columns
     integer, dimension(:,:), pointer :: pter 
   end type btofc_pters
   type (btofc_pters), dimension(:), allocatable, private :: btofc_blk_offset
                                       ! offset in btoc send array (-1) where 
                                       ! (blockid, bcid, k) column should be packed in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob receive array (-1) from which
                                       ! (blockid, bcid, k) column should be unpacked in
                                       ! chunk_to_block alltoallv

   type (btofc_pters), dimension(:), allocatable, private :: btofc_chk_offset
                                       ! offset in btoc receive array (-1) from which
                                       ! (lcid, i, k) data should be unpacked in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob send array (-1) where
                                       ! (lcid, i, k) data should be packed in
                                       ! chunk_to_block alltoallv

! miscellaneous phys_grid data
   integer, private :: dp_coup_steps   ! number of swaps in transpose algorithm
   integer, dimension(:), private, allocatable :: dp_coup_proc
                                       ! swap partner in each step of 
                                       !  transpose algorithm
   logical :: physgrid_set = .false.   ! flag indicates physics grid has been set
   integer, private :: max_nproc_smpx  ! maximum number of processes assigned to a
                                       !  single virtual SMP used to define physics 
                                       !  load balancing
   integer, private :: nproc_busy_d    ! number of processes active during the dynamics
                                       !  (assigned a dynamics block)

! Physics grid decomposition options:  
! -1: each chunk is a dynamics block
!  0: chunk definitions and assignments do not require interprocess comm.
!  1: chunk definitions and assignments do not require internode comm.
!  2: chunk definitions and assignments may require communication between all processes
!  3: chunk definitions and assignments only require communication with one other process
!  4: concatenated blocks, no load balancing, no interprocess communication
   integer, private, parameter :: min_lbal_opt = -1
   integer, private, parameter :: max_lbal_opt = 5
   integer, private, parameter :: def_lbal_opt = 2               ! default
   integer, private :: lbal_opt = def_lbal_opt

! Physics grid load balancing options:  
!  0: assign columns to chunks as single columns, wrap mapped across chunks
!  1: use (day/night; north/south) twin algorithm to determine load-balanced pairs of 
!       columns and assign columns to chunks in pairs, wrap mapped
   integer, private, parameter :: min_twin_alg = 0
   integer, private, parameter :: max_twin_alg = 1
   integer, private, parameter :: def_twin_alg_lonlat = 1         ! default
   integer, private, parameter :: def_twin_alg_unstructured = 0
   integer, private :: twin_alg = def_twin_alg_lonlat

! target number of chunks per thread
   integer, private, parameter :: min_chunks_per_thread = 1
   integer, private, parameter :: def_chunks_per_thread = &
                                    min_chunks_per_thread         ! default
   integer, private :: chunks_per_thread = def_chunks_per_thread

! Dynamics/physics transpose method for nonlocal load-balance:
! -1: use "0" if max_nproc_smpx and nproc_busy_d are both > npes/2; otherwise use "1"
!  0: use mpi_alltoallv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
!  11-13: use mod_comm, choosing any of several methods internal to mod_comm.
!      The method within mod_comm (denoted mod_method) has possible values 0,1,2 and
!      is set according to mod_method = phys_alltoall - modmin_alltoall, where
!      modmin_alltoall is 11.
   integer, private, parameter :: min_alltoall = -1
   integer, private, parameter :: max_alltoall = 3
   integer, private, parameter :: def_alltoall = -1                ! default
   integer, private :: phys_alltoall = def_alltoall

contains
!========================================================================
  integer function get_nlcols_p()
    get_nlcols_p = nlcols
  end function get_nlcols_p

! Subprogram not used   integer function get_clon_p_tot()
! Subprogram not used     get_clon_p_tot = clon_p_tot
! Subprogram not used   end function get_clon_p_tot
! Subprogram not used   integer function get_clat_p_tot()
! Subprogram not used     get_clat_p_tot = clat_p_tot
! Subprogram not used   end function get_clat_p_tot

  subroutine phys_grid_init( )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Physics mapping initialization routine:  
    ! 
    ! Method: 
    ! 
    ! Author: John Drake and Patrick Worley
    ! 
    !-----------------------------------------------------------------------
    use pmgrid, only: plev
    use dyn_grid, only: get_block_bounds_d, &
         get_block_gcol_d, get_block_gcol_cnt_d, &
         get_block_levels_d, get_block_lvl_cnt_d, &
         get_block_owner_d, &
         get_gcol_block_d, get_gcol_block_cnt_d, &
         get_horiz_grid_dim_d, get_horiz_grid_d
       use spmd_utils, only: pair, ceil2
    !
    !------------------------------Arguments--------------------------------
    !
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: i, j, jb, k, p             ! loop indices
    integer :: pre_i                      ! earlier index in loop iteration
    integer :: clat_p_dex, clon_p_dex     ! indices into unique lat. and lon. arrays
    integer :: maxblksiz                  ! maximum number of columns in a dynamics block
    integer :: beg_dex, end_dex           ! index range
    integer :: cid, lcid                  ! global and local chunk ids
    integer :: max_ncols                  ! upper bound on number of columns in a block
    integer :: ncols                      ! number of columns in current chunk
    integer :: curgcol, curgcol_d         ! current global column index
    integer :: firstblock, lastblock      ! global block indices
    integer :: blksiz                     ! current block size
    integer :: glbcnt, curcnt             ! running grid point counts
    integer :: curp                       ! current process id
    integer :: block_cnt                  ! number of blocks containing data
    ! for a given vertical column
    integer :: numlvl                     ! number of vertical levels in block 
    ! column
    integer :: levels(plev+1)             ! vertical level indices
    integer :: owner_d                    ! process owning given block column
    integer :: owner_p                    ! process owning given chunk column
    integer :: blockids(plev+1)           ! block indices
    integer :: bcids(plev+1)              ! block column indices
    integer :: glon, glat                 ! global (lon,lat) indices
    integer :: ntmp1, ntmp2               ! work variables

    logical :: clon_wrap                  ! flag used in initializing lat_p, lon_p

    ! column surface area (from dynamics)
    real(r8), dimension(:), allocatable :: area_d 

    ! column integration weight (from dynamics)
    real(r8), dimension(:), allocatable :: wght_d 

    ! chunk global ordering
    integer, dimension(:), allocatable :: pchunkid                   

    ! permutation array used in physics column sorting;
    ! reused later as work space in (lbal_opt == -1) logic
    integer, dimension(:), allocatable :: cdex

    ! latitudes and longitudes and column area for dynamics columns
    real(r8), dimension(:), allocatable :: clat_d
    real(r8), dimension(:), allocatable :: clon_d
    real(r8) :: clat_p_tmp
    real(r8) :: clon_p_tmp

    integer lons(2), lats(2)

    call t_adj_detailf(-2)
    call t_startf("phys_grid_init")

    !-----------------------------------------------------------------------
    !
    ! Initialize physics grid, using dynamics grid
    ! a) column coordinates

    call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
    ngcols = hdim1_d*hdim2_d
    allocate( clat_d(1:ngcols) )
    allocate( clon_d(1:ngcols) )
    allocate( cdex(1:ngcols) )
    clat_d = 100000.0_r8
    clon_d = 100000.0_r8
    call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d)

    ! count number of "real" column indices
    ngcols_p = 0
    do i=1,ngcols
       if (clon_d(i) < 100000.0_r8) then
          ngcols_p = ngcols_p + 1
       endif
    enddo

    ! sort over longitude and identify unique longitude coordinates
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clon_d,descend=.false.)
    clon_p_tmp = clon_d(cdex(1))
    clon_p_tot = 1

    do i=2,ngcols_p
       if (clon_d(cdex(i)) > clon_p_tmp) then
          clon_p_tot = clon_p_tot + 1
          clon_p_tmp = clon_d(cdex(i))
       endif
    enddo

    allocate( clon_p(1:clon_p_tot) )
    allocate( clon_p_cnt(1:clon_p_tot) )

    pre_i = 1
    clon_p_tot = 1
    clon_p(1) = clon_d(cdex(1))
    do i=2,ngcols_p
       if (clon_d(cdex(i)) > clon_p(clon_p_tot)) then
          clon_p_cnt(clon_p_tot) = i-pre_i
          pre_i = i
          clon_p_tot = clon_p_tot + 1
          clon_p(clon_p_tot) = clon_d(cdex(i))
       endif
    enddo
    clon_p_cnt(clon_p_tot) = (ngcols_p+1)-pre_i

    ! sort over latitude and identify unique latitude coordinates
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clat_d,descend=.false.)
    clat_p_tmp = clat_d(cdex(1))
    clat_p_tot = 1
    do i=2,ngcols_p
       if (clat_d(cdex(i)) > clat_p_tmp) then
          clat_p_tot = clat_p_tot + 1
          clat_p_tmp = clat_d(cdex(i))
       endif
    enddo

    allocate( clat_p(1:clat_p_tot) )
    allocate( clat_p_cnt(1:clat_p_tot) )
    allocate( clat_p_idx(1:clat_p_tot) )

    pre_i = 1
    clat_p_tot = 1
    clat_p(1) = clat_d(cdex(1))
    do i=2,ngcols_p
       if (clat_d(cdex(i)) > clat_p(clat_p_tot)) then
          clat_p_cnt(clat_p_tot) = i-pre_i
          pre_i = i
          clat_p_tot = clat_p_tot + 1
          clat_p(clat_p_tot) = clat_d(cdex(i))
       endif
    enddo
    clat_p_cnt(clat_p_tot) = (ngcols_p+1)-pre_i

    clat_p_idx(1) = 1
    do j=2,clat_p_tot
       clat_p_idx(j) = clat_p_idx(j-1) + clat_p_cnt(j-1)
    enddo

    ! sort by longitude within latitudes
    end_dex = 0
    do j=1,clat_p_tot
       beg_dex = end_dex + 1
       end_dex = end_dex + clat_p_cnt(j)
       call IndexSort(cdex(beg_dex:end_dex),clon_d,descend=.false.)
    enddo

    ! Early clean-up, to minimize memory high water mark
    ! (not executing find_partner or find_twin)
    if (((twin_alg .ne. 1) .and. (lbal_opt .ne. 3)) .or. &
        (lbal_opt .eq. -1)) deallocate( clat_p_cnt)

    ! save "longitude within latitude" column ordering
    ! and determine mapping from unsorted global column index to 
    ! unique latitude/longitude indices
    allocate( lat_p(1:ngcols) )
    allocate( lon_p(1:ngcols) )
    allocate( dyn_to_latlon_gcol_map(1:ngcols) )
    if (lbal_opt .ne. -1) allocate( latlon_to_dyn_gcol_map(1:ngcols_p) )

    clat_p_dex = 1
    lat_p = -1
    dyn_to_latlon_gcol_map = -1
    do i=1,ngcols_p
       if (lbal_opt .ne. -1) latlon_to_dyn_gcol_map(i) = cdex(i)
       dyn_to_latlon_gcol_map(cdex(i)) = i

       do while ((clat_p(clat_p_dex) < clat_d(cdex(i))) .and. &
                 (clat_p_dex < clat_p_tot))
          clat_p_dex = clat_p_dex + 1
       enddo
       lat_p(cdex(i)) = clat_p_dex
    enddo

    ! sort by latitude within longitudes
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clon_d,descend=.false.)
    end_dex = 0
    do i=1,clon_p_tot
       beg_dex = end_dex + 1
       end_dex = end_dex + clon_p_cnt(i)
       call IndexSort(cdex(beg_dex:end_dex),clat_d,descend=.false.)
    enddo

    ! Early clean-up, to minimize memory high water mark
    ! (not executing find_twin)
    if ((twin_alg .ne. 1) .or. (lbal_opt .eq. -1)) deallocate( clon_p_cnt )

    ! save "latitude within longitude" column ordering
    ! (only need in find_twin)
    if ((twin_alg .eq. 1) .and. (lbal_opt .ne. -1)) &
       allocate( lonlat_to_dyn_gcol_map(1:ngcols_p) )

    clon_p_dex = 1
    lon_p = -1
    do i=1,ngcols_p
       if ((twin_alg .eq. 1) .and. (lbal_opt .ne. -1)) &
         lonlat_to_dyn_gcol_map(i) = cdex(i)
       do while ((clon_p(clon_p_dex) < clon_d(cdex(i))) .and. &
                 (clon_p_dex < clon_p_tot))
          clon_p_dex = clon_p_dex + 1
       enddo
       lon_p(cdex(i)) = clon_p_dex
    enddo

    ! Clean-up
    deallocate( clat_d )
    deallocate( clon_d )
    deallocate( cdex )

    !
    ! Determine block index bounds
    !
    call get_block_bounds_d(firstblock,lastblock)

    ! Allocate storage to save number of chunks and columns assigned to each
    ! process during chunk creation and assignment
    !
    allocate( npchunks(0:npes-1) )
    allocate( gs_col_num(0:npes-1) )
    npchunks(:) = 0
    gs_col_num(:) = 0

    !
    ! Option -1: each dynamics block is a single chunk
    !            
    if (lbal_opt == -1) then
       !
       ! Check that pcols >= maxblksiz
       !
       maxblksiz = 0
       do jb=firstblock,lastblock
          maxblksiz = max(maxblksiz,get_block_gcol_cnt_d(jb))
       enddo
       if (pcols < maxblksiz) then
	  write(iulog,*) 'pcols = ',pcols, ' maxblksiz=',maxblksiz
          call endrun ('PHYS_GRID_INIT error: phys_loadbalance -1 specified but PCOLS < MAXBLKSIZ')
       endif

       !
       ! Determine total number of chunks
       !
       nchunks = (lastblock-firstblock+1)

       !
       ! Set max virtual SMP node size
       !
       max_nproc_smpx = 1

       !
       ! Allocate and initialize chunks data structure
       !
       allocate( cdex(1:maxblksiz) )
       allocate( chunks(1:nchunks) )

       do cid=1,nchunks
          ! get number of global column indices in block
          max_ncols = get_block_gcol_cnt_d(cid+firstblock-1)
          ! fill cdex array with global indices from current block
          call get_block_gcol_d(cid+firstblock-1,max_ncols,cdex)

          ncols = 0
          do i=1,max_ncols
             ! check whether global index is for a column that dynamics
             ! intends to pass to the physics
             curgcol_d = cdex(i)
             if (dyn_to_latlon_gcol_map(curgcol_d) .ne. -1) then
                ! yes - then save the information
                ncols = ncols + 1
                chunks(cid)%gcol(ncols) = curgcol_d
                chunks(cid)%lat(ncols) = lat_p(curgcol_d)
                chunks(cid)%lon(ncols) = lon_p(curgcol_d)
             endif
          enddo
          chunks(cid)%ncols = ncols
       enddo

       ! Clean-up
       deallocate( cdex )
       deallocate( lat_p )
       deallocate( lon_p )

       !
       ! Specify parallel decomposition 
       !
       do cid=1,nchunks
          p = get_block_owner_d(cid+firstblock-1)
          chunks(cid)%owner = p
          npchunks(p)       = npchunks(p) + 1
          gs_col_num(p)     = gs_col_num(p) + chunks(cid)%ncols
       enddo
       !
       ! Set flag indicating columns in physics and dynamics 
       ! decompositions reside on the same processes
       !
       local_dp_map = .true. 
       !
    else
       !
       ! Option == 0: split local blocks into chunks,
       !               while attempting to create load-balanced chunks.
       !               Does not work with vertically decomposed blocks.
       !               (default)
       ! Option == 1: split SMP-local blocks into chunks,
       !               while attempting to create load-balanced chunks.
       !               Does not work with vertically decomposed blocks.
       ! Option == 2: load balance chunks with respect to diurnal and
       !               seaonsal cycles and wth respect to latitude, 
       !               and assign chunks to processes
       !               in a way that attempts to minimize communication costs
       ! Option == 3: divide processes into pairs and split 
       !               blocks assigned to these pairs into 
       !               chunks, attempting to create load-balanced chunks.
       !               The process pairs are chosen to maximize load balancing
       !               opportunities.
       !               Does not work with vertically decomposed blocks.
       ! Option == 4: concatenate local blocks, then
       !               divide into chunks.
       !               Does not work with vertically decomposed blocks.
       ! Option == 5: split indiviudal blocks into chunks,
       !               assigning columns using block ordering
       !
       !
       ! Allocate and initialize chunks data structure, then
       ! assign chunks to processes.
       !
       call create_chunks(lbal_opt, chunks_per_thread)

       ! Early clean-up, to minimize memory high water mark
       deallocate( lat_p )
       deallocate( lon_p )
       deallocate( latlon_to_dyn_gcol_map )
       if  (twin_alg .eq. 1) deallocate( lonlat_to_dyn_gcol_map )
       if  (twin_alg .eq. 1) deallocate( clon_p_cnt )
       if ((twin_alg .eq. 1) .or. (lbal_opt .eq. 3)) deallocate( clat_p_cnt )

       !
       ! Determine whether dynamics and physics decompositions
       ! are colocated, not requiring any interprocess communication
       ! in the coupling.
       local_dp_map = .true.   
       do cid=1,nchunks
          do i=1,chunks(cid)%ncols
             curgcol_d = chunks(cid)%gcol(i)
             block_cnt = get_gcol_block_cnt_d(curgcol_d)
             call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
             do jb=1,block_cnt
                owner_d = get_block_owner_d(blockids(jb)) 
                if (owner_d .ne. chunks(cid)%owner) then
                   local_dp_map = .false.   
                endif
             enddo
          enddo
       enddo
    endif
    !
    ! Allocate and initialize data structures for gather/scatter
    !  
    allocate( pgcols(1:ngcols_p) )
    allocate( gs_col_offset(0:npes) )
    allocate( pchunkid(0:npes) )

    ! Initialize pchunkid and gs_col_offset by summing 
    ! number of chunks and columns per process, respectively
    pchunkid(0) = 0
    gs_col_offset(0) = 0
    do p=1,npes-1
       pchunkid(p)      = pchunkid(p-1)      + npchunks(p-1)
       gs_col_offset(p) = gs_col_offset(p-1) + gs_col_num(p-1)
    enddo
    
    ! Determine local ordering via "process id" bin sort
    do cid=1,nchunks
       p = chunks(cid)%owner
       pchunkid(p) = pchunkid(p) + 1

       chunks(cid)%lcid = pchunkid(p) + lastblock

       curgcol = gs_col_offset(p)
       do i=1,chunks(cid)%ncols
          curgcol = curgcol + 1
          pgcols(curgcol)%chunk = cid
          pgcols(curgcol)%ccol = i
       enddo
       gs_col_offset(p) = curgcol
    enddo

    ! Reinitialize pchunkid and gs_col_offset (for real)
    pchunkid(0) = 1
    gs_col_offset(0) = 1
    do p=1,npes-1
       pchunkid(p)      = pchunkid(p-1)      + npchunks(p-1)
       gs_col_offset(p) = gs_col_offset(p-1) + gs_col_num(p-1)
    enddo
    pchunkid(npes)      = pchunkid(npes-1)      + npchunks(npes-1)
    gs_col_offset(npes) = gs_col_offset(npes-1) + gs_col_num(npes-1)

    ! Save local information
    ! (Local chunk index range chosen so that it does not overlap 
    !  {begblock,...,endblock})
    ! 
    nlcols   = gs_col_num(iam)
    nlchunks = npchunks(iam)
    begchunk = pchunkid(iam)   + lastblock
    endchunk = pchunkid(iam+1) + lastblock - 1
    !
    allocate( lchunks(begchunk:endchunk) )
    do cid=1,nchunks
       if (chunks(cid)%owner == iam) then
          lcid = chunks(cid)%lcid
          lchunks(lcid)%ncols = chunks(cid)%ncols
          lchunks(lcid)%cid   = cid
          do i=1,chunks(cid)%ncols
             lchunks(lcid)%gcol(i) = chunks(cid)%gcol(i)
          enddo
       endif
    enddo

    deallocate( pchunkid )
    deallocate( npchunks )
    !
    !-----------------------------------------------------------------------
    !
    ! Initialize physics grid, using dynamics grid
    ! b) column area and integration weight

    allocate( area_d(1:ngcols) )
    allocate( wght_d(1:ngcols) )
    area_d = 0.0_r8
    wght_d = 0.0_r8

    call get_horiz_grid_d(ngcols, area_d_out=area_d, wght_d_out=wght_d)


    if ( abs(sum(area_d) - 4.0_r8*pi) > 1.e-10_r8 ) then
       write(iulog,*) ' ERROR: sum of areas on globe does not equal 4*pi'
       write(iulog,*) ' sum of areas = ', sum(area_d), sum(area_d)-4.0_r8*pi
       call endrun('phys_grid')
    end if

    if ( abs(sum(wght_d) - 4.0_r8*pi) > 1.e-10_r8 ) then
       write(iulog,*) ' ERROR: sum of integration weights on globe does not equal 4*pi'
       write(iulog,*) ' sum of weights = ', sum(wght_d), sum(wght_d)-4.0_r8*pi
       call endrun('phys_grid')
    end if

    do lcid=begchunk,endchunk
       do i=1,lchunks(lcid)%ncols
          lchunks(lcid)%area(i) = area_d(lchunks(lcid)%gcol(i))
          lchunks(lcid)%wght(i) = wght_d(lchunks(lcid)%gcol(i))
       enddo
    enddo

    deallocate( area_d )
    deallocate( wght_d )

    if (.not. local_dp_map) then
       !
       ! allocate and initialize data structures for transposes
       !  
       allocate( btofc_blk_num(0:npes-1) )
       btofc_blk_num = 0
       allocate( btofc_blk_offset(firstblock:lastblock) )
       do jb = firstblock,lastblock
          nullify( btofc_blk_offset(jb)%pter )
       enddo
       !
       glbcnt = 0
       curcnt = 0
       curp = 0
       do curgcol=1,ngcols_p
          cid = pgcols(curgcol)%chunk
          i   = pgcols(curgcol)%ccol
          owner_p   = chunks(cid)%owner
          do while (curp < owner_p)
             btofc_blk_num(curp) = curcnt
             curcnt = 0
             curp = curp + 1
          enddo
          curgcol_d = chunks(cid)%gcol(i)
          block_cnt = get_gcol_block_cnt_d(curgcol_d)
          call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
          do jb = 1,block_cnt
             owner_d = get_block_owner_d(blockids(jb))
             if (iam == owner_d) then
                if (.not. associated(btofc_blk_offset(blockids(jb))%pter)) then
                   blksiz = get_block_gcol_cnt_d(blockids(jb))
                   numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
                   btofc_blk_offset(blockids(jb))%ncols = blksiz
                   btofc_blk_offset(blockids(jb))%nlvls = numlvl
                   allocate( btofc_blk_offset(blockids(jb))%pter(blksiz,numlvl) )
                endif
                do k=1,btofc_blk_offset(blockids(jb))%nlvls
                   btofc_blk_offset(blockids(jb))%pter(bcids(jb),k) = glbcnt
                   curcnt = curcnt + 1
                   glbcnt = glbcnt + 1
                enddo
             endif
          enddo
       enddo
       btofc_blk_num(curp) = curcnt
       block_buf_nrecs = glbcnt
       !  
       allocate( btofc_chk_num(0:npes-1) )
       btofc_chk_num = 0
       allocate( btofc_chk_offset(begchunk:endchunk) )
       do lcid=begchunk,endchunk
          ncols = lchunks(lcid)%ncols
          btofc_chk_offset(lcid)%ncols = ncols
          btofc_chk_offset(lcid)%nlvls = pver+1
          allocate( btofc_chk_offset(lcid)%pter(ncols,pver+1) )
       enddo
       !
       curcnt = 0
       glbcnt = 0
       do p=0,npes-1
          do curgcol=gs_col_offset(iam),gs_col_offset(iam+1)-1
             cid  = pgcols(curgcol)%chunk
             owner_p  = chunks(cid)%owner
             if (iam == owner_p) then
                i    = pgcols(curgcol)%ccol
                lcid = chunks(cid)%lcid
                curgcol_d = chunks(cid)%gcol(i)
                block_cnt = get_gcol_block_cnt_d(curgcol_d)
                call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
                do jb = 1,block_cnt
                   owner_d = get_block_owner_d(blockids(jb))
                   if (p == owner_d) then
                      numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
                      call get_block_levels_d(blockids(jb),bcids(jb),numlvl,levels)
                      do k=1,numlvl
                         btofc_chk_offset(lcid)%pter(i,levels(k)+1) = glbcnt
                         curcnt = curcnt + 1
                         glbcnt = glbcnt + 1
                      enddo
                   endif
                enddo
             endif
          enddo
          btofc_chk_num(p) = curcnt
          curcnt = 0
       enddo
       chunk_buf_nrecs = glbcnt
       !
       ! Precompute swap partners and number of steps in point-to-point
       ! implementations of alltoall algorithm.
       ! First, determine number of swaps.
       !
       dp_coup_steps = 0
       do i=1,ceil2(npes)-1
          p = pair(npes,i,iam)
          if (p >= 0) then
             if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
                dp_coup_steps = dp_coup_steps + 1
             end if
          end if
       end do
       !
       ! Second, determine swap partners.
       !
       allocate( dp_coup_proc(dp_coup_steps) )
       dp_coup_steps = 0
       do i=1,ceil2(npes)-1
          p = pair(npes,i,iam)
          if (p >= 0) then
             if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
                dp_coup_steps = dp_coup_steps + 1
                dp_coup_proc(dp_coup_steps) = p
             end if
          end if
       end do
       !
    endif

    ! Final clean-up
    deallocate( gs_col_offset )
    ! (if eliminate get_lon_xxx, can also deallocate
    !  clat_p_idx, and grid_latlon?))

    !
    physgrid_set = .true.   ! Set flag indicating physics grid is now set
    !
    if (masterproc) then
       write(iulog,*) 'PHYS_GRID_INIT:  Using PCOLS=',pcols,     &
            '  phys_loadbalance=',lbal_opt,            &
            '  phys_twin_algorithm=',twin_alg,         &
            '  phys_alltoall=',phys_alltoall,          &
            '  chunks_per_thread=',chunks_per_thread
    endif
    !

    call t_stopf("phys_grid_init")
    call t_adj_detailf(+2)
    return
  end subroutine phys_grid_init

!========================================================================

! Subprogram not used subroutine phys_grid_find_col(lat, lon, owner, lcid, icol)
! Subprogram not used 
! Subprogram not used    !----------------------------------------------------------------------- 
! Subprogram not used    ! 
! Subprogram not used    ! Purpose: Find the global column closest to the point specified by lat
! Subprogram not used    !          and lon.  Return indices of owning process, local chunk, and 
! Subprogram not used    !          column.
! Subprogram not used    ! 
! Subprogram not used    ! Authors: Phil Rasch / Patrick Worley / B. Eaton
! Subprogram not used    ! 
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: lat, lon    ! requested location in degrees
! Subprogram not used    integer, intent(out) :: owner       ! rank of chunk owner
! Subprogram not used    integer, intent(out) :: lcid      ! local chunk index
! Subprogram not used    integer, intent(out) :: icol        ! column index within the chunk
! Subprogram not used 
! Subprogram not used    ! local
! Subprogram not used    real(r8) dist2           ! the distance (in radians**2 from lat, lon)
! Subprogram not used    real(r8) distmin         ! the distance (in radians**2 from closest column)
! Subprogram not used    real(r8) latr, lonr      ! lat, lon (in radians) of requested location
! Subprogram not used    real(r8) clat, clon      ! lat, lon (in radians) of column being tested
! Subprogram not used    real(r8) const
! Subprogram not used 
! Subprogram not used    integer i
! Subprogram not used    integer cid
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ! Check that input lat and lon are in valid range
! Subprogram not used    if (lon < 0.0_r8 .or. lon >= 360._r8 .or. &
! Subprogram not used        lat < -90._r8 .or. lat > 90._r8) then
! Subprogram not used       if (masterproc) then
! Subprogram not used          write(iulog,*) &
! Subprogram not used             'phys_grid_find_col: ERROR: lon must satisfy 0.<=lon<360. and lat must satisfy -90<=lat<=90.'
! Subprogram not used          write(iulog,*) &
! Subprogram not used             'input lon=', lon, '  input lat=', lat
! Subprogram not used       endif
! Subprogram not used       call endrun('phys_grid_find_col: input ERROR')
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    const = 180._r8/pi            ! degrees per radian
! Subprogram not used    latr = lat/const              ! to radians
! Subprogram not used    lonr = lon/const              ! to radians
! Subprogram not used 
! Subprogram not used    owner   = -999
! Subprogram not used    lcid  = -999
! Subprogram not used    icol    = -999
! Subprogram not used    distmin = 1.e10_r8
! Subprogram not used 
! Subprogram not used    ! scan all chunks for closest point to lat, lon
! Subprogram not used    do cid = 1, nchunks
! Subprogram not used       do i = 1, chunks(cid)%ncols
! Subprogram not used          clat = clat_p(chunks(cid)%lat(i))
! Subprogram not used          clon = clon_p(chunks(cid)%lon(i))
! Subprogram not used          dist2 = (clat-latr)**2 + (clon-lonr)**2
! Subprogram not used          if (dist2 < distmin ) then
! Subprogram not used             distmin = dist2
! Subprogram not used             owner = chunks(cid)%owner
! Subprogram not used             lcid = chunks(cid)%lcid
! Subprogram not used             icol = i
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine phys_grid_find_col

!========================================================================

! Subprogram not used subroutine phys_grid_find_cols(lat, lon, nclosest, owner, lcid, icol, distmin, mlats, mlons)
! Subprogram not used 
! Subprogram not used    !----------------------------------------------------------------------- 
! Subprogram not used    ! 
! Subprogram not used    ! Purpose: Find the global columns closest to the point specified by lat
! Subprogram not used    !          and lon.  Return indices of owning process, local chunk, and 
! Subprogram not used    !          column.
! Subprogram not used    ! 
! Subprogram not used    ! Authors: Phil Rasch / Patrick Worley / B. Eaton
! Subprogram not used    ! 
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used    use physconst,    only : rearth
! Subprogram not used    
! Subprogram not used    real(r8), intent(in) :: lat, lon            ! requested location in degrees
! Subprogram not used    integer, intent(in)  :: nclosest            ! number of closest points to find
! Subprogram not used    integer, intent(out) :: owner(nclosest)     ! rank of chunk owner
! Subprogram not used    integer, intent(out) :: lcid(nclosest)      ! local chunk index
! Subprogram not used    integer, intent(out) :: icol(nclosest)      ! column index within the chunk
! Subprogram not used    real(r8),intent(out) :: distmin(nclosest)   ! the distance (m) of the closest column(s)
! Subprogram not used    real(r8),intent(out) :: mlats(nclosest)     ! the latitude of the closest column(s)
! Subprogram not used    real(r8),intent(out) :: mlons(nclosest)     ! the longitude of the closest column(s)
! Subprogram not used 
! Subprogram not used    ! local
! Subprogram not used    real(r8) dist2           ! the distance (in radians**2 from lat, lon)
! Subprogram not used    real(r8) latr, lonr      ! lat, lon (in radians) of requested location
! Subprogram not used    real(r8) clat, clon      ! lat, lon (in radians) of column being tested
! Subprogram not used    real(r8) const
! Subprogram not used 
! Subprogram not used    integer i, j
! Subprogram not used    integer cid
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ! Check that input lat and lon are in valid range
! Subprogram not used    if (lon < 0.0_r8 .or. lon >= 360._r8 .or. &
! Subprogram not used        lat < -90._r8 .or. lat > 90._r8) then
! Subprogram not used       if (masterproc) then
! Subprogram not used          write(iulog,*) &
! Subprogram not used             'phys_grid_find_cols: ERROR: lon must satisfy 0.<=lon<360. and lat must satisfy -90<=lat<=90.'
! Subprogram not used          write(iulog,*) &
! Subprogram not used             'input lon=', lon, '  input lat=', lat
! Subprogram not used       endif
! Subprogram not used       call endrun('phys_grid_find_cols: input ERROR')
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    const = 180._r8/pi            ! degrees per radian
! Subprogram not used    latr = lat/const              ! to radians
! Subprogram not used    lonr = lon/const              ! to radians
! Subprogram not used 
! Subprogram not used    owner(:)   = -999
! Subprogram not used    lcid(:)    = -999
! Subprogram not used    icol(:)    = -999
! Subprogram not used    mlats(:)   = -999
! Subprogram not used    mlons(:)   = -999
! Subprogram not used    distmin(:) = 1.e10_r8
! Subprogram not used 
! Subprogram not used    ! scan all chunks for closest point to lat, lon
! Subprogram not used    do cid = 1, nchunks
! Subprogram not used       do i = 1, chunks(cid)%ncols
! Subprogram not used          clat = clat_p(chunks(cid)%lat(i))
! Subprogram not used          clon = clon_p(chunks(cid)%lon(i))
! Subprogram not used          dist2 = acos(sin(latr) * sin(clat) + cos(latr) * cos(clat) * cos(clon - lonr)) * rearth       
! Subprogram not used          
! Subprogram not used          do j = nclosest, 1, -1
! Subprogram not used             if (dist2 < distmin(j)) then
! Subprogram not used             
! Subprogram not used                if (j < nclosest) then
! Subprogram not used                  distmin(j+1) = distmin(j)
! Subprogram not used                  owner(j+1)   = owner(j)
! Subprogram not used                  lcid(j+1)    = lcid(j)
! Subprogram not used                  icol(j+1)    = icol(j)
! Subprogram not used                  mlats(j+1)   = mlats(j)
! Subprogram not used                  mlons(j+1)    = mlons(j)
! Subprogram not used                end if
! Subprogram not used              
! Subprogram not used                distmin(j) = dist2
! Subprogram not used                owner(j)   = chunks(cid)%owner
! Subprogram not used                lcid(j)    = chunks(cid)%lcid
! Subprogram not used                icol(j)    = i
! Subprogram not used                mlats(j)   = clat * const
! Subprogram not used                mlons(j)   = clon * const
! Subprogram not used             else
! Subprogram not used                exit
! Subprogram not used             end if
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used    end do
! Subprogram not used    
! Subprogram not used end subroutine phys_grid_find_cols
!
!========================================================================

logical function phys_grid_initialized ()
!----------------------------------------------------------------------- 
! 
! Purpose: Identify whether phys_grid has been called yet or not
! 
! Method: Return physgrid_set
! 
! Author: Pat Worley
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   phys_grid_initialized = physgrid_set
!
   return
   end function phys_grid_initialized

!
!========================================================================
!
   subroutine phys_grid_defaultopts(phys_loadbalance_out, &
                                    phys_twin_algorithm_out, &
                                    phys_alltoall_out, &
                                    phys_chnk_per_thd_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
   use dycore, only: dycore_is
!------------------------------Arguments--------------------------------
     ! physics load balancing option
     integer, intent(out), optional :: phys_loadbalance_out
     ! algorithm to use when determining column pairs to assign to chunks
     integer, intent(out), optional :: phys_twin_algorithm_out
     ! alltoall option
     integer, intent(out), optional :: phys_alltoall_out
     ! number of chunks per thread
     integer, intent(out), optional :: phys_chnk_per_thd_out
!-----------------------------------------------------------------------
     if ( present(phys_loadbalance_out) ) then
       phys_loadbalance_out = def_lbal_opt
     endif
     if ( present(phys_twin_algorithm_out) ) then
       if (dycore_is('UNSTRUCTURED')) then
          phys_twin_algorithm_out = def_twin_alg_unstructured
       else
          phys_twin_algorithm_out = def_twin_alg_lonlat
       endif
     endif
     if ( present(phys_alltoall_out) ) then
       phys_alltoall_out = def_alltoall
     endif
     if ( present(phys_chnk_per_thd_out) ) then
       phys_chnk_per_thd_out = def_chunks_per_thread
     endif
   end subroutine phys_grid_defaultopts
!
!========================================================================
!
   subroutine phys_grid_setopts(phys_loadbalance_in, &
                                phys_twin_algorithm_in, &
                                phys_alltoall_in,    &
                                phys_chnk_per_thd_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
   use spmd_utils, only: phys_mirror_decomp_req
!------------------------------Arguments--------------------------------
     ! physics load balancing option
     integer, intent(in), optional :: phys_loadbalance_in
     ! option to use load balanced column pairs
     integer, intent(in), optional :: phys_twin_algorithm_in
     ! alltoall option
     integer, intent(in), optional :: phys_alltoall_in
     ! number of chunks per thread
     integer, intent(in), optional :: phys_chnk_per_thd_in
!-----------------------------------------------------------------------
     if ( present(phys_loadbalance_in) ) then
        lbal_opt = phys_loadbalance_in
        if ((lbal_opt < min_lbal_opt).or.(lbal_opt > max_lbal_opt)) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_loadbalance=', &
                 phys_loadbalance_in,                             &
                 '  is out of range.  It must be between ',       &
                 min_lbal_opt,' and ',max_lbal_opt
           endif
           call endrun
        endif
        if (lbal_opt .eq. 3) then
           phys_mirror_decomp_req = .true.
        else
           phys_mirror_decomp_req = .false.
        endif
     endif
!
     if ( present(phys_twin_algorithm_in) ) then
        twin_alg = phys_twin_algorithm_in
        if ((twin_alg < min_twin_alg).or.(twin_alg > max_twin_alg)) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_twin_algorithm=', &
                 phys_twin_algorithm_in,                             &
                 '  is out of range.  It must be between ',       &
                 min_twin_alg,' and ',max_twin_alg
           endif
           call endrun
        endif
     endif
!
     if ( present(phys_alltoall_in) ) then
        phys_alltoall = phys_alltoall_in
        if (((phys_alltoall .lt. min_alltoall) .or.    &
             (phys_alltoall .gt. max_alltoall))        &
           ) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SET_OPTS:  ERROR:  phys_alltoall=',   &
                  phys_alltoall_in,                               &
                  '  is out of range.  It must be between ',      &
                  min_alltoall,' and ',max_alltoall
           endif
           call endrun
        endif
     endif
!
     if ( present(phys_chnk_per_thd_in) ) then
        chunks_per_thread = phys_chnk_per_thd_in
        if (chunks_per_thread < min_chunks_per_thread) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_chnk_per_thd=',&
                 phys_chnk_per_thd_in,                            &
                 ' is too small.  It must not be smaller than ',  &
                 min_chunks_per_thread
           endif
           call endrun
        endif
     endif
   end subroutine phys_grid_setopts
!
!========================================================================
!
! Subprogram not used    subroutine get_chunk_indices_p(index_beg, index_end)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return range of indices for local chunks
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(out) :: index_beg  ! first index used for local chunks
! Subprogram not used    integer, intent(out) :: index_end  ! last index used for local chunks
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    index_beg = begchunk
! Subprogram not used    index_end = endchunk
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine get_chunk_indices_p
!
!========================================================================
!
   subroutine get_gcol_all_p(lcid, latdim, gcols)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all global column indices for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
     integer, intent(in)  :: lcid        ! local chunk id
     integer, intent(in)  :: latdim      ! declared size of output array

     integer, intent(out) :: gcols(:)    ! array of global latitude indices
!---------------------------Local workspace-----------------------------
     integer :: i                        ! loop index
     
!-----------------------------------------------------------------------
     gcols=-1
     do i=1,lchunks(lcid)%ncols
        gcols(i) = lchunks(lcid)%gcol(i)
     enddo
     return
   end subroutine get_gcol_all_p

!
!========================================================================
!
   integer function get_gcol_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global physics column index for chunk column
! 
! Method: 
! 
! Author: Jim Edwards / Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!-----------------------------------------------------------------------
   get_gcol_p = lchunks(lcid)%gcol(col)
   
   return
   end function get_gcol_p

!
!========================================================================

! Subprogram not used    subroutine get_gcol_vec_p(lcid, lth, cols, gcols)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return global physics column indices for set of chunk columns
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: lth           ! number of column indices
! Subprogram not used    integer, intent(in)  :: cols(lth)     ! column indices
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: gcols(lth)    ! array of global physics 
! Subprogram not used                                          !  columns indices
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i                          ! loop index
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    do i=1,lth
! Subprogram not used      gcols(i) = lchunks(lcid)%gcol(cols(i))
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine get_gcol_vec_p

!
!========================================================================
!
   integer function get_ncols_p(lcid)
!----------------------------------------------------------------------- 
! 
! Purpose: Return number of columns in chunk given the local chunk id.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid      ! local chunk id

!---------------------------Local workspace-----------------------------
   integer              :: cid       ! global chunk id

!-----------------------------------------------------------------------
   get_ncols_p = lchunks(lcid)%ncols

   return
   end function get_ncols_p
!
!========================================================================
!
   subroutine get_lat_all_p(lcid, latdim, lats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all global latitude indices for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: latdim        ! declared size of output array

   integer, intent(out) :: lats(latdim)  ! array of global latitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     lats(i) = chunks(cid)%lat(i)
   enddo

   return
   end subroutine get_lat_all_p
!
!========================================================================

! Subprogram not used    subroutine get_lat_vec_p(lcid, lth, cols, lats)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return global latitude indices for set of chunk columns
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: lth           ! number of column indices
! Subprogram not used    integer, intent(in)  :: cols(lth)     ! column indices
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: lats(lth)     ! array of global latitude indices
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i                          ! loop index
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    cid = lchunks(lcid)%cid
! Subprogram not used    do i=1,lth
! Subprogram not used      lats(i) = chunks(cid)%lat(cols(i))
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine get_lat_vec_p
!
!========================================================================

   integer function get_lat_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global latitude index for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   get_lat_p = chunks(cid)%lat(col)

   return
   end function get_lat_p
!
!========================================================================
!
   subroutine get_lon_all_p(lcid, londim, lons)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Was: Return all global longitude indices for chunk
!  Now: Return all longitude offsets (+1) for chunk. These are offsets
!       in ordered list of global columns from first
!       column with given latitude to column with given latitude
!       and longitude. This corresponds to the usual longitude indices
!       for full and reduced lon/lat grids.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: londim        ! declared size of output array

   integer, intent(out) :: lons(londim)  ! array of global longitude 
                                         !  indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: lat                        ! latitude index
   integer :: cid                        ! global chunk id
   integer :: gcol                       ! global column id in latlon 
                                         !  ordering

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     lat  = chunks(cid)%lat(i)
     gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(i))
     lons(i) = (gcol - clat_p_idx(lat)) + 1
   enddo

   return
   end subroutine get_lon_all_p
!
!========================================================================

! Subprogram not used    subroutine get_lon_vec_p(lcid, lth, cols, lons)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used !  Was: Return global longitude indices for set of chunk columns.
! Subprogram not used !  Now: Return longitude offsets (+1) for set of chunk columns. 
! Subprogram not used !       These are offsets in ordered list of global columns from first
! Subprogram not used !       column with given latitude to column with given latitude
! Subprogram not used !       and longitude. This corresponds to the usual longitude indices
! Subprogram not used !       for full and reduced lon/lat grids.
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: lth           ! number of column indices
! Subprogram not used    integer, intent(in)  :: cols(lth)     ! column indices
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: lons(lth)     ! array of global longitude indices
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i                          ! loop index
! Subprogram not used    integer :: lat                        ! latitude index
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: gcol                       ! global column id in latlon 
! Subprogram not used                                          !  ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    cid = lchunks(lcid)%cid
! Subprogram not used    do i=1,lth
! Subprogram not used      lat = chunks(cid)%lat(cols(i))
! Subprogram not used      gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(i))
! Subprogram not used      lons(i) = (gcol - clat_p_idx(lat)) + 1
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine get_lon_vec_p
!
!========================================================================

   integer function get_lon_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Was: Return global longitude index for chunk column.
!  Now: Return longitude offset (+1) for chunk column. This is the 
!       offset in ordered list of global columns from first
!       column with given latitude to column with given latitude
!       and longitude. This corresponds to the usual longitude index
!       for full and reduced lon/lat grids.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: cid                        ! global chunk id
   integer :: lat                        ! latitude index
   integer :: gcol                       ! global column id in latlon 
                                         !  ordering

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   lat = chunks(cid)%lat(col)
   gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(col))
   get_lon_p = (gcol - clat_p_idx(lat)) + 1

   return
   end function get_lon_p
!
!========================================================================
!
   subroutine get_rlat_all_p(lcid, rlatdim, rlats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all latitudes (in radians) for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid           ! local chunk id
   integer, intent(in)  :: rlatdim        ! declared size of output array

   real(r8), intent(out) :: rlats(rlatdim)! array of latitudes

!---------------------------Local workspace-----------------------------
   integer :: i                           ! loop index
   integer :: cid                         ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     rlats(i) = clat_p(chunks(cid)%lat(i))
   enddo

   return
   end subroutine get_rlat_all_p
!
!========================================================================
!
   subroutine get_area_all_p(lcid, rdim, area)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all areas for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: rdim          ! declared size of output array

   real(r8), intent(out) :: area(rdim)   ! array of areas

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index

!-----------------------------------------------------------------------
   do i=1,lchunks(lcid)%ncols
     area(i) = lchunks(lcid)%area(i)
   enddo

   return
   end subroutine get_area_all_p
!
!========================================================================
!
! Subprogram not used    real(r8) function get_area_p(lcid, col)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return area for chunk column
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: col           ! column index
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    get_area_p = lchunks(lcid)%area(col)
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end function get_area_p
!
!========================================================================
!
   subroutine get_wght_all_p(lcid, rdim, wght)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all integration weights for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: rdim          ! declared size of output array

   real(r8), intent(out) :: wght(rdim)   ! array of integration weights

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index

!-----------------------------------------------------------------------
   do i=1,lchunks(lcid)%ncols
     wght(i) = lchunks(lcid)%wght(i)
   enddo

   return
   end subroutine get_wght_all_p
!
!========================================================================
!
! Subprogram not used    real(r8) function get_wght_p(lcid, col)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return integration weight for chunk column
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: col           ! column index
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    get_wght_p = lchunks(lcid)%wght(col)
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end function get_wght_p
!
!========================================================================
!
! Subprogram not used    subroutine get_rlat_vec_p(lcid, lth, cols, rlats)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return latitudes (in radians) for set of chunk columns
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: lth           ! number of column indices
! Subprogram not used    integer, intent(in)  :: cols(lth)     ! column indices
! Subprogram not used 
! Subprogram not used    real(r8), intent(out) :: rlats(lth)   ! array of latitudes
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i                          ! loop index
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    cid = lchunks(lcid)%cid
! Subprogram not used    do i=1,lth
! Subprogram not used      rlats(i) = clat_p(chunks(cid)%lat(cols(i)))
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine get_rlat_vec_p
!
!========================================================================

! Subprogram not used    real(r8) function get_rlat_p(lcid, col)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return latitude (in radians) for chunk column
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: col           ! column index
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    cid = lchunks(lcid)%cid
! Subprogram not used    get_rlat_p = clat_p(chunks(cid)%lat(col))
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end function get_rlat_p
!
!========================================================================
!
   subroutine get_rlon_all_p(lcid, rlondim, rlons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all longitudes (in radians) for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid           ! local chunk id
   integer, intent(in)  :: rlondim        ! declared size of output array

   real(r8), intent(out) :: rlons(rlondim)! array of longitudes

!---------------------------Local workspace-----------------------------
   integer :: i                           ! loop index
   integer :: cid                         ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     rlons(i) = clon_p(chunks(cid)%lon(i))
   enddo

   return
   end subroutine get_rlon_all_p
!
!========================================================================

! Subprogram not used    subroutine get_rlon_vec_p(lcid, lth, cols, rlons)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return longitudes (in radians) for set of chunk columns
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid         ! local chunk id
! Subprogram not used    integer, intent(in)  :: lth           ! number of column indices
! Subprogram not used    integer, intent(in)  :: cols(lth)     ! column indices
! Subprogram not used 
! Subprogram not used    real(r8), intent(out) :: rlons(lth)   ! array of longitudes
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i                          ! loop index
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    cid = lchunks(lcid)%cid
! Subprogram not used    do i=1,lth
! Subprogram not used      rlons(i) = clon_p(chunks(cid)%lon(cols(i)))
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine get_rlon_vec_p
!
!========================================================================

! Subprogram not used    real(r8) function get_rlon_p(lcid, col)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return longitude (in radians) for chunk column
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use ppgrid
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: lcid          ! local chunk id
! Subprogram not used    integer, intent(in)  :: col           ! column index
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    cid = lchunks(lcid)%cid
! Subprogram not used    get_rlon_p = clon_p(chunks(cid)%lon(col))
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end function get_rlon_p
!
!========================================================================
!
!  integer function get_gcol_owner_p(gcol)
!----------------------------------------------------------------------- 
! 
! Purpose: Return owner of physics column with indicate index
! 
! Method: 
! 
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!  integer, intent(in)  :: gcol     ! physics column index
!
!-----------------------------------------------------------------------
!
!  get_gcol_owner_p = chunks(knuhcs(gcol)%chunkid)%owner
!
!  return
!  end function get_gcol_owner_p
!
!========================================================================

!  subroutine buff_to_chunk(fdim,mdim,lbuff,localchunks)
!-----------------------------------------------------------------------
!
! Purpose: Copy from local buffer 
!          to local chunk data structure.
!          Needed for cpl6.
!
! Method:
!
! Author: Pat Worley and Robert Jacob
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!  integer, intent(in) :: fdim      ! declared length of first lbuff dimension
!  integer, intent(in) :: mdim      ! declared length of middle lbuff dimension
!  real(r8), intent(in) :: lbuff(fdim, mdim) ! local lon/lat buffer
!
!  real(r8), intent(out):: localchunks(pcols,mdim,begchunk:endchunk) ! local chunks
!
!
!---------------------------Local workspace-----------------------------
!  integer :: i,j,m,n                      ! loop indices
!
!  integer, save :: numcols = 0
!  integer, allocatable, save :: columnid(:), chunkid(:)
!-----------------------------------------------------------------------
!
!  if (numcols .eq. 0) then
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!           endif
!        endif
!     enddo
!     allocate(columnid(1:n))
!     allocate(chunkid(1:n))
!
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!              columnid(n) = knuhcs(i)%col
!              chunkid(n)  = chunks(knuhcs(i)%chunkid)%lcid
!           endif
!        endif
!     end do
!
!     numcols = n
!  endif
!
!  if (numcols .gt. fdim) call endrun('buff_to_chunk')
!  do m=1,mdim
!dir$ concurrent
!dir$ prefervector, preferstream
!     do n = 1, numcols
!        localchunks(columnid(n),m,chunkid(n)) = lbuff(n,m)
!     end do
!  end do
!
!  return
!  end subroutine buff_to_chunk
!
!========================================================================

! Subprogram not used    subroutine scatter_field_to_chunk(fdim,mdim,ldim, &
! Subprogram not used                                      hdim1d,globalfield,localchunks)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Distribute field
! Subprogram not used !          to decomposed chunk data structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    integer, intent(in) :: hdim1d    ! declared first horizontal index 
! Subprogram not used                                     ! dimension
! Subprogram not used    real(r8), intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
! Subprogram not used                                     ! global field
! Subprogram not used 
! Subprogram not used    real(r8), intent(out):: localchunks(fdim,pcols,mdim, &
! Subprogram not used                                        begchunk:endchunk,ldim) 
! Subprogram not used                                     ! local chunks
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: f,i,m,l,p                  ! loop indices
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: lcid                       ! local chunk id
! Subprogram not used    integer :: lid                        ! local column index
! Subprogram not used    integer :: gcol                       ! global column index
! Subprogram not used    integer :: h1                         ! first horizontal dimension index
! Subprogram not used    integer :: h2                         ! second horizontal dimension index
! Subprogram not used 
! Subprogram not used    real(r8) gfield_p(fdim,mdim,ldim,ngcols) 
! Subprogram not used                                          ! vector to be scattered
! Subprogram not used    real(r8) lfield_p(fdim,mdim,ldim,nlcols) 
! Subprogram not used                                          ! local component of scattered
! Subprogram not used                                          !  vector
! Subprogram not used    integer :: displs(0:npes-1)           ! scatter displacements
! Subprogram not used    integer :: sndcnts(0:npes-1)          ! scatter send counts
! Subprogram not used    integer :: recvcnt                    ! scatter receive count
! Subprogram not used    integer :: beglcol                    ! beginning index for local columns
! Subprogram not used                                          !  in global column ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (hdim1d < hdim1_d) then
! Subprogram not used       write(iulog,*) "<stdin>",2030,hdim1d,hdim1_d
! Subprogram not used       call endrun ('SCATTER_FIELD_TO_CHUNK error: hdim1d < hdim1_d')
! Subprogram not used    endif
! Subprogram not used    localchunks(:,:,:,:,:) = 0
! Subprogram not used    displs(0) = 0
! Subprogram not used    sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
! Subprogram not used    beglcol = 0
! Subprogram not used    do p=1,npes-1
! Subprogram not used      displs(p) = displs(p-1) + sndcnts(p-1)
! Subprogram not used      sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
! Subprogram not used      if (p <= iam) then
! Subprogram not used         beglcol = beglcol + gs_col_num(p-1)
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used    recvcnt = fdim*mdim*ldim*nlcols
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used 
! Subprogram not used ! copy field into global (process-ordered) chunked data structure
! Subprogram not used 
! Subprogram not used       do l=1,ldim
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used !DIR$ PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used          do i=1,ngcols_p
! Subprogram not used             cid  = pgcols(i)%chunk
! Subprogram not used             lid  = pgcols(i)%ccol
! Subprogram not used             gcol = chunks(cid)%gcol(lid)
! Subprogram not used             h2   = (gcol-1)/hdim1_d + 1
! Subprogram not used             h1   = mod((gcol-1),hdim1_d) + 1
! Subprogram not used             do m=1,mdim
! Subprogram not used                do f=1,fdim
! Subprogram not used                   gfield_p(f,m,l,i) = &
! Subprogram not used                      globalfield(f, h1, m, h2, l)
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used ! scatter to other processes
! Subprogram not used ! (pgcols ordering consistent with begchunk:endchunk 
! Subprogram not used ! local ordering)
! Subprogram not used 
! Subprogram not used    call t_barrierf('sync_scat_ftoc', mpicom)
! Subprogram not used    call mpiscatterv(gfield_p, sndcnts, displs, mpir8, &
! Subprogram not used                     lfield_p, recvcnt, mpir8, 0, mpicom)
! Subprogram not used 
! Subprogram not used ! copy into local chunked data structure
! Subprogram not used 
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used !DIR$ PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used    do i=1,nlcols
! Subprogram not used       cid = pgcols(beglcol+i)%chunk
! Subprogram not used       lcid = chunks(cid)%lcid
! Subprogram not used       lid = pgcols(beglcol+i)%ccol
! Subprogram not used       do l=1,ldim
! Subprogram not used          do m=1,mdim
! Subprogram not used             do f=1,fdim
! Subprogram not used                localchunks(f,lid,m,lcid,l) = &
! Subprogram not used                  lfield_p(f, m, l, i)
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine scatter_field_to_chunk
!========================================================================

! Subprogram not used    subroutine scatter_field_to_chunk4(fdim,mdim,ldim, &
! Subprogram not used                                       hdim1d,globalfield,localchunks)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Distribute field
! Subprogram not used !          to decomposed chunk data structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    integer, intent(in) :: hdim1d    ! declared first horizontal index 
! Subprogram not used                                     ! dimension
! Subprogram not used    real(r4), intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
! Subprogram not used                                     ! global field
! Subprogram not used 
! Subprogram not used    real(r4), intent(out):: localchunks(fdim,pcols,mdim, &
! Subprogram not used                                        begchunk:endchunk,ldim) 
! Subprogram not used                                     ! local chunks
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: f,i,m,l,p                  ! loop indices
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: lcid                       ! local chunk id
! Subprogram not used    integer :: lid                        ! local column index
! Subprogram not used    integer :: gcol                       ! global column index
! Subprogram not used    integer :: h1                         ! first horizontal dimension index
! Subprogram not used    integer :: h2                         ! second horizontal dimension index
! Subprogram not used 
! Subprogram not used    real(r4) gfield_p(fdim,mdim,ldim,ngcols) 
! Subprogram not used                                          ! vector to be scattered
! Subprogram not used    real(r4) lfield_p(fdim,mdim,ldim,nlcols) 
! Subprogram not used                                          ! local component of scattered
! Subprogram not used                                          !  vector
! Subprogram not used    integer :: displs(0:npes-1)           ! scatter displacements
! Subprogram not used    integer :: sndcnts(0:npes-1)          ! scatter send counts
! Subprogram not used    integer :: recvcnt                    ! scatter receive count
! Subprogram not used    integer :: beglcol                    ! beginning index for local columns
! Subprogram not used                                          !  in global column ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (hdim1d < hdim1_d) then
! Subprogram not used       call endrun ('SCATTER_FIELD_TO_CHUNK4 error: hdim1d < hdim1_d')
! Subprogram not used    endif
! Subprogram not used    displs(0) = 0
! Subprogram not used    sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
! Subprogram not used    beglcol = 0
! Subprogram not used    do p=1,npes-1
! Subprogram not used      displs(p) = displs(p-1) + sndcnts(p-1)
! Subprogram not used      sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
! Subprogram not used      if (p <= iam) then
! Subprogram not used         beglcol = beglcol + gs_col_num(p-1)
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used    recvcnt = fdim*mdim*ldim*nlcols
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       ! copy field into global (process-ordered) chunked data structure
! Subprogram not used       do l=1,ldim
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used !DIR$ PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used          do i=1,ngcols_p
! Subprogram not used             cid  = pgcols(i)%chunk
! Subprogram not used             lid  = pgcols(i)%ccol
! Subprogram not used             gcol = chunks(cid)%gcol(lid)
! Subprogram not used             h2   = (gcol-1)/hdim1_d + 1
! Subprogram not used             h1   = mod((gcol-1),hdim1_d) + 1
! Subprogram not used             do m=1,mdim
! Subprogram not used                do f=1,fdim
! Subprogram not used                   gfield_p(f,m,l,i) = &
! Subprogram not used                      globalfield(f, h1, m, h2, l)
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used ! scatter to other processes
! Subprogram not used ! (pgcols ordering consistent with begchunk:endchunk 
! Subprogram not used !  local ordering)
! Subprogram not used 
! Subprogram not used    call t_barrierf('sync_scat_ftoc', mpicom)
! Subprogram not used    call mpiscatterv(gfield_p, sndcnts, displs, mpir4, &
! Subprogram not used                     lfield_p, recvcnt, mpir4, 0, mpicom)
! Subprogram not used 
! Subprogram not used ! copy into local chunked data structure
! Subprogram not used 
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used !DIR$ PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used    do i=1,nlcols
! Subprogram not used       cid = pgcols(beglcol+i)%chunk
! Subprogram not used       lcid = chunks(cid)%lcid
! Subprogram not used       lid = pgcols(beglcol+i)%ccol
! Subprogram not used       do l=1,ldim
! Subprogram not used          do m=1,mdim
! Subprogram not used             do f=1,fdim
! Subprogram not used                localchunks(f,lid,m,lcid,l) = &
! Subprogram not used                  lfield_p(f, m, l, i)
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine scatter_field_to_chunk4
!========================================================================

! Subprogram not used    subroutine scatter_field_to_chunk_int(fdim,mdim,ldim, &
! Subprogram not used                                          hdim1d,globalfield,localchunks)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Distribute field
! Subprogram not used !          to decomposed chunk data structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    integer, intent(in) :: hdim1d    ! declared first horizontal index 
! Subprogram not used                                     ! dimension
! Subprogram not used    integer, intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
! Subprogram not used                                     ! global field
! Subprogram not used 
! Subprogram not used    integer, intent(out):: localchunks(fdim,pcols,mdim, &
! Subprogram not used                                        begchunk:endchunk,ldim) 
! Subprogram not used                                     ! local chunks
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: f,i,m,l,p                  ! loop indices
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: lcid                       ! local chunk id
! Subprogram not used    integer :: lid                        ! local column index
! Subprogram not used    integer :: gcol                       ! global column index
! Subprogram not used    integer :: h1                         ! first horizontal dimension index
! Subprogram not used    integer :: h2                         ! second horizontal dimension index
! Subprogram not used 
! Subprogram not used    integer gfield_p(fdim,mdim,ldim,ngcols) 
! Subprogram not used                                          ! vector to be scattered
! Subprogram not used    integer lfield_p(fdim,mdim,ldim,nlcols) 
! Subprogram not used                                          ! local component of scattered
! Subprogram not used                                          !  vector
! Subprogram not used    integer :: displs(0:npes-1)           ! scatter displacements
! Subprogram not used    integer :: sndcnts(0:npes-1)          ! scatter send counts
! Subprogram not used    integer :: recvcnt                    ! scatter receive count
! Subprogram not used    integer :: beglcol                    ! beginning index for local columns
! Subprogram not used                                          !  in global column ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (hdim1d < hdim1_d) then
! Subprogram not used       call endrun ('SCATTER_FIELD_TO_CHUNK_INT error: hdim1d < hdim1_d')
! Subprogram not used    endif
! Subprogram not used    displs(0) = 0
! Subprogram not used    sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
! Subprogram not used    beglcol = 0
! Subprogram not used    do p=1,npes-1
! Subprogram not used      displs(p) = displs(p-1) + sndcnts(p-1)
! Subprogram not used      sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
! Subprogram not used      if (p <= iam) then
! Subprogram not used         beglcol = beglcol + gs_col_num(p-1)
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used    recvcnt = fdim*mdim*ldim*nlcols
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used 
! Subprogram not used ! copy field into global (process-ordered) chunked data structure
! Subprogram not used 
! Subprogram not used       do l=1,ldim
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used !DIR$ PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used          do i=1,ngcols_p
! Subprogram not used             cid = pgcols(i)%chunk
! Subprogram not used             lid = pgcols(i)%ccol
! Subprogram not used             gcol = chunks(cid)%gcol(lid)
! Subprogram not used             h2   = (gcol-1)/hdim1_d + 1
! Subprogram not used             h1   = mod((gcol-1),hdim1_d) + 1
! Subprogram not used             do m=1,mdim
! Subprogram not used                do f=1,fdim
! Subprogram not used                   gfield_p(f,m,l,i) = &
! Subprogram not used                      globalfield(f, h1, m, h2, l)
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used ! scatter to other processes
! Subprogram not used ! (pgcols ordering consistent with begchunk:endchunk 
! Subprogram not used !  local ordering)
! Subprogram not used 
! Subprogram not used    call t_barrierf('sync_scat_ftoc', mpicom)
! Subprogram not used    call mpiscatterv(gfield_p, sndcnts, displs, mpiint, &
! Subprogram not used                     lfield_p, recvcnt, mpiint, 0, mpicom)
! Subprogram not used 
! Subprogram not used ! copy into local chunked data structure
! Subprogram not used 
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used !DIR$ PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used    do i=1,nlcols
! Subprogram not used       cid = pgcols(beglcol+i)%chunk
! Subprogram not used       lcid = chunks(cid)%lcid
! Subprogram not used       lid = pgcols(beglcol+i)%ccol
! Subprogram not used       do l=1,ldim
! Subprogram not used          do m=1,mdim
! Subprogram not used             do f=1,fdim
! Subprogram not used                localchunks(f,lid,m,lcid,l) = &
! Subprogram not used                  lfield_p(f, m, l, i)
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine scatter_field_to_chunk_int
!
!========================================================================
!
!  subroutine chunk_to_buff(fdim,mdim,localchunks,lbuff)
!
!-----------------------------------------------------------------------
!
! Purpose: Copy from local chunk data structure
!          to local buffer.  Needed for cpl6.
!          (local = assigned to same process)
!
! Method:
!
! Author: Pat Worley and Robert Jacob
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!  integer, intent(in) :: fdim      ! declared length of first lbuff dimension
!  integer, intent(in) :: mdim      ! declared length of middle lbuff dimension
!  real(r8), intent(in):: localchunks(pcols,mdim, begchunk:endchunk) ! local chunks
!
!  real(r8), intent(out) :: lbuff(fdim,mdim) ! local buff
!
!---------------------------Local workspace-----------------------------
!  integer :: i,j,m,n                  ! loop indices
!
!  integer, save :: numcols = 0
!  integer, allocatable, save :: columnid(:), chunkid(:)
!-----------------------------------------------------------------------
!
!  if (numcols .eq. 0) then
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!           endif
!        endif
!     enddo
!     allocate(columnid(1:n))
!     allocate(chunkid(1:n))
!
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!              columnid(n) = knuhcs(i)%col
!              chunkid(n)  = chunks(knuhcs(i)%chunkid)%lcid
!           endif
!        endif
!     end do
!
!     numcols = n
!  endif
!
!  if (numcols .gt. fdim) call endrun('chunk_to_buff')
!  do m=1,mdim
!dir$ concurrent
!dir$ prefervector, preferstream
!     do n = 1, numcols
!        lbuff(n,m) = localchunks(columnid(n),m,chunkid(n))
!     end do
!  end do
!
!  return
!  end subroutine chunk_to_buff
!
!
!========================================================================
!
! Subprogram not used    subroutine gather_chunk_to_field(fdim,mdim,ldim, &
! Subprogram not used                                      hdim1d,localchunks,globalfield)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Reconstruct field
! Subprogram not used !          from decomposed chunk data structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use spmd_utils,    only: fc_gatherv
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    integer, intent(in) :: hdim1d    ! declared first horizontal index 
! Subprogram not used                                     ! dimension
! Subprogram not used    real(r8), intent(in):: localchunks(fdim,pcols,mdim, &
! Subprogram not used                                       begchunk:endchunk,ldim) 
! Subprogram not used                                     ! local chunks
! Subprogram not used 
! Subprogram not used    real(r8), intent(out) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
! Subprogram not used                                     ! global field
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: f,i,m,l,p                  ! loop indices
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: lcid                       ! local chunk id
! Subprogram not used    integer :: lid                        ! local column index
! Subprogram not used    integer :: gcol                       ! global column index
! Subprogram not used    integer :: h1                         ! first horizontal dimension index
! Subprogram not used    integer :: h2                         ! second horizontal dimension index
! Subprogram not used 
! Subprogram not used    real(r8) gfield_p(fdim,mdim,ldim,ngcols) 
! Subprogram not used                                          ! vector to be gathered
! Subprogram not used    real(r8) lfield_p(fdim,mdim,ldim,nlcols) 
! Subprogram not used                                          ! local component of gather
! Subprogram not used                                          !  vector
! Subprogram not used    integer :: displs(0:npes-1)           ! gather displacements
! Subprogram not used    integer :: rcvcnts(0:npes-1)          ! gather receive count
! Subprogram not used    integer :: sendcnt                    ! gather send counts
! Subprogram not used    integer :: beglcol                    ! beginning index for local columns
! Subprogram not used                                          !  in global column ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (hdim1d < hdim1_d) then
! Subprogram not used       call endrun ('GATHER_CHUNK_TO_FIELD error: hdim1d < hdim1_d')
! Subprogram not used    endif
! Subprogram not used    displs(0) = 0
! Subprogram not used    rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
! Subprogram not used    beglcol = 0
! Subprogram not used    do p=1,npes-1
! Subprogram not used      displs(p) = displs(p-1) + rcvcnts(p-1)
! Subprogram not used      rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
! Subprogram not used      if (p <= iam) then
! Subprogram not used         beglcol = beglcol + gs_col_num(p-1)
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used    sendcnt = fdim*mdim*ldim*nlcols
! Subprogram not used 
! Subprogram not used ! copy into local gather data structure
! Subprogram not used 
! Subprogram not used    do l=1,ldim
! Subprogram not used !DIR$ PREFERVECTOR, PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used       do i=1,nlcols
! Subprogram not used          cid = pgcols(beglcol+i)%chunk
! Subprogram not used          lcid = chunks(cid)%lcid
! Subprogram not used          lid = pgcols(beglcol+i)%ccol
! Subprogram not used          do m=1,mdim
! Subprogram not used             do f=1,fdim
! Subprogram not used                lfield_p(f, m, l, i) = &
! Subprogram not used                   localchunks(f,lid,m,lcid,l)
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used ! gather from other processes
! Subprogram not used 
! Subprogram not used    call t_barrierf('sync_gath_ctof', mpicom)
! Subprogram not used    call fc_gatherv(lfield_p, sendcnt, mpir8, &
! Subprogram not used                    gfield_p, rcvcnts, displs, mpir8, 0, mpicom)
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used 
! Subprogram not used ! copy gathered columns into lon/lat field
! Subprogram not used 
! Subprogram not used !DIR$ PREFERVECTOR, PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used       do i=1,ngcols_p
! Subprogram not used          cid  = pgcols(i)%chunk
! Subprogram not used          lid  = pgcols(i)%ccol
! Subprogram not used          gcol = chunks(cid)%gcol(lid)
! Subprogram not used          h2   = (gcol-1)/hdim1_d + 1
! Subprogram not used          h1   = mod((gcol-1),hdim1_d) + 1
! Subprogram not used          do l=1,ldim
! Subprogram not used             do m=1,mdim
! Subprogram not used                do f=1,fdim
! Subprogram not used                   globalfield(f, h1, m, h2, l)    &
! Subprogram not used                   = gfield_p(f,m,l,i)
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used    call mpibarrier(mpicom)
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine gather_chunk_to_field

!
!========================================================================
!
! Subprogram not used    subroutine gather_chunk_to_field4 (fdim,mdim,ldim, &
! Subprogram not used                                       hdim1d,localchunks,globalfield)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Reconstruct field
! Subprogram not used !          from decomposed chunk data structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use spmd_utils,    only: fc_gathervr4
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    integer, intent(in) :: hdim1d    ! declared first horizontal index 
! Subprogram not used                                     ! dimension
! Subprogram not used    real(r4), intent(in):: localchunks(fdim,pcols,mdim, &
! Subprogram not used                                       begchunk:endchunk,ldim) 
! Subprogram not used                                     ! local chunks
! Subprogram not used 
! Subprogram not used    real(r4), intent(out) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
! Subprogram not used                                     ! global field
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: f,i,m,l,p                  ! loop indices
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: lcid                       ! local chunk id
! Subprogram not used    integer :: lid                        ! local column index
! Subprogram not used    integer :: gcol                       ! global column index
! Subprogram not used    integer :: h1                         ! first horizontal dimension index
! Subprogram not used    integer :: h2                         ! second horizontal dimension index
! Subprogram not used 
! Subprogram not used    real(r4) gfield_p(fdim,mdim,ldim,ngcols) 
! Subprogram not used                                          ! vector to be gathered
! Subprogram not used    real(r4) lfield_p(fdim,mdim,ldim,nlcols) 
! Subprogram not used                                          ! local component of gather
! Subprogram not used                                          !  vector
! Subprogram not used    integer :: displs(0:npes-1)           ! gather displacements
! Subprogram not used    integer :: rcvcnts(0:npes-1)          ! gather receive count
! Subprogram not used    integer :: sendcnt                    ! gather send counts
! Subprogram not used    integer :: beglcol                    ! beginning index for local columns
! Subprogram not used                                          !  in global column ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (hdim1d < hdim1_d) then
! Subprogram not used       call endrun ('GATHER_CHUNK_TO_FIELD4 error: hdim1d < hdim1_d')
! Subprogram not used    endif
! Subprogram not used    displs(0) = 0
! Subprogram not used    rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
! Subprogram not used    beglcol = 0
! Subprogram not used    do p=1,npes-1
! Subprogram not used      displs(p) = displs(p-1) + rcvcnts(p-1)
! Subprogram not used      rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
! Subprogram not used      if (p <= iam) then
! Subprogram not used         beglcol = beglcol + gs_col_num(p-1)
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used    sendcnt = fdim*mdim*ldim*nlcols
! Subprogram not used 
! Subprogram not used ! copy into local gather data structure
! Subprogram not used 
! Subprogram not used    do l=1,ldim
! Subprogram not used !DIR$ PREFERVECTOR, PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used       do i=1,nlcols
! Subprogram not used          cid = pgcols(beglcol+i)%chunk
! Subprogram not used          lcid = chunks(cid)%lcid
! Subprogram not used          lid = pgcols(beglcol+i)%ccol
! Subprogram not used          do m=1,mdim
! Subprogram not used             do f=1,fdim
! Subprogram not used                lfield_p(f, m, l, i) = &
! Subprogram not used                   localchunks(f,lid,m,lcid,l)
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used ! gather from other processes
! Subprogram not used 
! Subprogram not used    call t_barrierf('sync_gath_ctof', mpicom)
! Subprogram not used    call fc_gathervr4(lfield_p, sendcnt, mpir4, &
! Subprogram not used                      gfield_p, rcvcnts, displs, mpir4, 0, mpicom)
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used 
! Subprogram not used ! copy gathered columns into lon/lat field
! Subprogram not used 
! Subprogram not used !DIR$ PREFERVECTOR, PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used       do i=1,ngcols_p
! Subprogram not used          cid  = pgcols(i)%chunk
! Subprogram not used          lid  = pgcols(i)%ccol
! Subprogram not used          gcol = chunks(cid)%gcol(lid)
! Subprogram not used          h2   = (gcol-1)/hdim1_d + 1
! Subprogram not used          h1   = mod((gcol-1),hdim1_d) + 1
! Subprogram not used          do l=1,ldim
! Subprogram not used             do m=1,mdim
! Subprogram not used                do f=1,fdim
! Subprogram not used                   globalfield(f, h1, m, h2, l)    &
! Subprogram not used                   = gfield_p(f,m,l,i)
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine gather_chunk_to_field4

!
!========================================================================
!
! Subprogram not used    subroutine gather_chunk_to_field_int (fdim,mdim,ldim, &
! Subprogram not used                                          hdim1d,localchunks,globalfield)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Reconstruct field
! Subprogram not used !          from decomposed chunk data structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use spmd_utils,    only: fc_gathervint
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    integer, intent(in) :: hdim1d    ! declared first horizontal index 
! Subprogram not used                                     ! dimension
! Subprogram not used    integer, intent(in):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) ! global field
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used 
! Subprogram not used    integer :: f,i,m,l,p                  ! loop indices
! Subprogram not used    integer :: cid                        ! global chunk id
! Subprogram not used    integer :: lcid                       ! local chunk id
! Subprogram not used    integer :: lid                        ! local column index
! Subprogram not used    integer :: gcol                       ! global column index
! Subprogram not used    integer :: h1                         ! first horizontal dimension index
! Subprogram not used    integer :: h2                         ! second horizontal dimension index
! Subprogram not used 
! Subprogram not used    integer gfield_p(fdim,mdim,ldim,ngcols) 
! Subprogram not used                                          ! vector to be gathered
! Subprogram not used    integer lfield_p(fdim,mdim,ldim,nlcols) 
! Subprogram not used                                          ! local component of gather
! Subprogram not used                                          !  vector
! Subprogram not used    integer :: displs(0:npes-1)           ! gather displacements
! Subprogram not used    integer :: rcvcnts(0:npes-1)          ! gather receive count
! Subprogram not used    integer :: sendcnt                    ! gather send counts
! Subprogram not used    integer :: beglcol                    ! beginning index for local columns
! Subprogram not used                                          !  in global column ordering
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (hdim1d < hdim1_d) then
! Subprogram not used       call endrun ('GATHER_CHUNK_TO_FIELD_INT error: hdim1d < hdim1_d')
! Subprogram not used    endif
! Subprogram not used    displs(0) = 0
! Subprogram not used    rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
! Subprogram not used    beglcol = 0
! Subprogram not used    do p=1,npes-1
! Subprogram not used      displs(p) = displs(p-1) + rcvcnts(p-1)
! Subprogram not used      rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
! Subprogram not used      if (p <= iam) then
! Subprogram not used         beglcol = beglcol + gs_col_num(p-1)
! Subprogram not used      endif
! Subprogram not used    enddo
! Subprogram not used    sendcnt = fdim*mdim*ldim*nlcols
! Subprogram not used 
! Subprogram not used ! copy into local gather data structure
! Subprogram not used 
! Subprogram not used    do l=1,ldim
! Subprogram not used !DIR$ PREFERVECTOR, PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used       do i=1,nlcols
! Subprogram not used          cid = pgcols(beglcol+i)%chunk
! Subprogram not used          lcid = chunks(cid)%lcid
! Subprogram not used          lid = pgcols(beglcol+i)%ccol
! Subprogram not used          do m=1,mdim
! Subprogram not used             do f=1,fdim
! Subprogram not used                lfield_p(f, m, l, i) = &
! Subprogram not used                   localchunks(f,lid,m,lcid,l)
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used ! gather from other processes
! Subprogram not used 
! Subprogram not used    call t_barrierf('sync_gath_ctof', mpicom)
! Subprogram not used    call fc_gathervint(lfield_p, sendcnt, mpiint, &
! Subprogram not used                       gfield_p, rcvcnts, displs, mpiint, 0, mpicom)
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used 
! Subprogram not used ! copy gathered columns into lon/lat field
! Subprogram not used 
! Subprogram not used !DIR$ PREFERVECTOR, PREFERSTREAM
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used       do i=1,ngcols_p
! Subprogram not used          cid  = pgcols(i)%chunk
! Subprogram not used          lid  = pgcols(i)%ccol
! Subprogram not used          gcol = chunks(cid)%gcol(lid)
! Subprogram not used          h2   = (gcol-1)/hdim1_d + 1
! Subprogram not used          h1   = mod((gcol-1),hdim1_d) + 1
! Subprogram not used          do l=1,ldim
! Subprogram not used             do m=1,mdim
! Subprogram not used                do f=1,fdim
! Subprogram not used                   globalfield(f, h1, m, h2, l)    &
! Subprogram not used                   = gfield_p(f,m,l,i)
! Subprogram not used                end do
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine gather_chunk_to_field_int

!
!========================================================================
!
! Subprogram not used    subroutine write_field_from_chunk(iu,fdim,mdim,ldim,localchunks)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used !                          
! Subprogram not used ! Purpose: Write field from decomposed chunk data 
! Subprogram not used !          structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: iu        ! logical unit
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used    real(r8), intent(in):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used 
! Subprogram not used    integer :: ioerr                 ! error return
! Subprogram not used 
! Subprogram not used    real(r8), allocatable :: globalfield(:,:,:,:,:)
! Subprogram not used                                     ! global field
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate(globalfield(fdim,hdim1_d,mdim,hdim2_d,ldim))
! Subprogram not used 
! Subprogram not used    call gather_chunk_to_field (fdim,mdim,ldim,hdim1_d,localchunks,globalfield)
! Subprogram not used                                
! Subprogram not used    if (masterproc) then
! Subprogram not used       write (iu,iostat=ioerr) globalfield
! Subprogram not used       if (ioerr /= 0 ) then
! Subprogram not used          write(iulog,*) 'WRITE_FIELD_FROM_CHUNK ioerror ', ioerr,' on i/o unit = ',iu
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    deallocate(globalfield)
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine write_field_from_chunk

!
!========================================================================
!
! Subprogram not used    subroutine read_chunk_from_field(iu,fdim,mdim,ldim,localchunks)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used !                          
! Subprogram not used ! Purpose: Write field from decomposed chunk data 
! Subprogram not used !          structure
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: iu        ! logical unit
! Subprogram not used    integer, intent(in) :: fdim      ! declared length of first dimension
! Subprogram not used    integer, intent(in) :: mdim      ! declared length of middle dimension
! Subprogram not used    integer, intent(in) :: ldim      ! declared length of last dimension
! Subprogram not used 
! Subprogram not used    real(r8), intent(out):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used 
! Subprogram not used    integer :: ioerr                 ! error return
! Subprogram not used 
! Subprogram not used    real(r8), allocatable :: globalfield(:,:,:,:,:)
! Subprogram not used                                     ! global field
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    allocate(globalfield(fdim,hdim1_d,mdim,hdim2_d,ldim))
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       read (iu,iostat=ioerr) globalfield
! Subprogram not used       if (ioerr /= 0 ) then
! Subprogram not used          write(iulog,*) 'READ_CHUNK_FROM_FIELD ioerror ', ioerr,' on i/o unit = ',iu
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    call scatter_field_to_chunk (fdim,mdim,ldim,hdim1_d,globalfield,localchunks)
! Subprogram not used 
! Subprogram not used    deallocate(globalfield)
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine read_chunk_from_field
!
!========================================================================

! Subprogram not used    subroutine transpose_block_to_chunk(record_size, block_buffer, &
! Subprogram not used                                        chunk_buffer, window)
! Subprogram not used                                        
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Transpose buffer containing decomposed 
! Subprogram not used !          fields to buffer
! Subprogram not used !          containing decomposed chunk data structures
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! Modified: Art Mirin, Jan 04, to add support for mod_comm
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use spmd_utils,    only: altalltoallv
! Subprogram not used !------------------------------Parameters-------------------------------
! Subprogram not used !
! Subprogram not used   integer, parameter :: msgtag  = 6000
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: record_size  ! per column amount of data 
! Subprogram not used    real(r8), intent(in) :: block_buffer(record_size*block_buf_nrecs)
! Subprogram not used                                        ! buffer of block data to be
! Subprogram not used                                        ! transposed
! Subprogram not used    real(r8), intent(out):: chunk_buffer(record_size*chunk_buf_nrecs)
! Subprogram not used                                        ! buffer of chunk data 
! Subprogram not used                                        ! transposed into
! Subprogram not used    integer, intent(in), optional :: window
! Subprogram not used                                        ! MPI-2 window id for
! Subprogram not used                                        ! chunk_buffer
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i, p                     ! loop indices
! Subprogram not used    integer :: bbuf_siz                 ! size of block_buffer
! Subprogram not used    integer :: cbuf_siz                 ! size of chunk_buffer
! Subprogram not used    integer :: lwindow                  ! placeholder for missing window
! Subprogram not used    integer :: lopt                     ! local copy of phys_alltoall
! Subprogram not used !
! Subprogram not used    logical, save :: first = .true.
! Subprogram not used    integer, allocatable, save :: sndcnts(:), sdispls(:)
! Subprogram not used    integer, allocatable, save :: rcvcnts(:), rdispls(:)
! Subprogram not used    integer, allocatable, save :: pdispls(:)
! Subprogram not used    integer, save :: prev_record_size = 0
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (first) then
! Subprogram not used ! Compute send/recv/put counts and displacements
! Subprogram not used       allocate(sndcnts(0:npes-1))
! Subprogram not used       allocate(sdispls(0:npes-1))
! Subprogram not used       allocate(rcvcnts(0:npes-1))
! Subprogram not used       allocate(rdispls(0:npes-1))
! Subprogram not used       allocate(pdispls(0:npes-1))
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       first = .false.
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    if (record_size .ne. prev_record_size) then
! Subprogram not used !
! Subprogram not used ! Compute send/recv/put counts and displacements
! Subprogram not used       sdispls(0) = 0
! Subprogram not used       sndcnts(0) = record_size*btofc_blk_num(0)
! Subprogram not used       do p=1,npes-1
! Subprogram not used         sdispls(p) = sdispls(p-1) + sndcnts(p-1)
! Subprogram not used         sndcnts(p) = record_size*btofc_blk_num(p)
! Subprogram not used       enddo
! Subprogram not used !
! Subprogram not used       rdispls(0) = 0
! Subprogram not used       rcvcnts(0) = record_size*btofc_chk_num(0)
! Subprogram not used       do p=1,npes-1
! Subprogram not used          rdispls(p) = rdispls(p-1) + rcvcnts(p-1)
! Subprogram not used          rcvcnts(p) = record_size*btofc_chk_num(p)
! Subprogram not used       enddo
! Subprogram not used !
! Subprogram not used       call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used       prev_record_size = record_size
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    call t_barrierf('sync_tran_btoc', mpicom)
! Subprogram not used    if (phys_alltoall < 0) then
! Subprogram not used       if ((max_nproc_smpx > npes/2) .and. (nproc_busy_d > npes/2)) then
! Subprogram not used          lopt = 0
! Subprogram not used       else
! Subprogram not used          lopt = 1
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       lopt = phys_alltoall
! Subprogram not used       if ((lopt .eq. 2) .and. ( .not. present(window) )) lopt = 1
! Subprogram not used    endif
! Subprogram not used    if (lopt < 4) then
! Subprogram not used !
! Subprogram not used       bbuf_siz = record_size*block_buf_nrecs
! Subprogram not used       cbuf_siz = record_size*chunk_buf_nrecs
! Subprogram not used       if ( present(window) ) then
! Subprogram not used          call altalltoallv(lopt, iam, npes,    &
! Subprogram not used                            dp_coup_steps, dp_coup_proc, &
! Subprogram not used                            block_buffer, bbuf_siz, sndcnts, sdispls, mpir8, &
! Subprogram not used                            chunk_buffer, cbuf_siz, rcvcnts, rdispls, mpir8, &
! Subprogram not used                            msgtag, pdispls, mpir8, window, mpicom)
! Subprogram not used       else
! Subprogram not used          call altalltoallv(lopt, iam, npes,    &
! Subprogram not used                            dp_coup_steps, dp_coup_proc, &
! Subprogram not used                            block_buffer, bbuf_siz, sndcnts, sdispls, mpir8, &
! Subprogram not used                            chunk_buffer, cbuf_siz, rcvcnts, rdispls, mpir8, &
! Subprogram not used                            msgtag, pdispls, mpir8, lwindow, mpicom)
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used    else
! Subprogram not used !
! Subprogram not used       call mpialltoallv(block_buffer, sndcnts, sdispls, mpir8, &
! Subprogram not used                         chunk_buffer, rcvcnts, rdispls, mpir8, &
! Subprogram not used                         mpicom)
! Subprogram not used !
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine transpose_block_to_chunk
!
!========================================================================

! Subprogram not used    subroutine block_to_chunk_send_pters(blockid, fdim, ldim, &
! Subprogram not used                                         record_size, pter)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return pointers into send buffer where column from decomposed 
! Subprogram not used !          fields should be copied to
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: blockid      ! block index
! Subprogram not used    integer, intent(in) :: fdim         ! first dimension of pter array
! Subprogram not used    integer, intent(in) :: ldim         ! last dimension of pter array
! Subprogram not used    integer, intent(in) :: record_size  ! per coordinate amount of data 
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: pter(fdim,ldim)  ! buffer offsets
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i, k                     ! loop indices
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if ((btofc_blk_offset(blockid)%ncols > fdim) .or. &
! Subprogram not used        (btofc_blk_offset(blockid)%nlvls > ldim)) then
! Subprogram not used       write(iulog,*) "BLOCK_TO_CHUNK_SEND_PTERS: pter array dimensions ", &
! Subprogram not used                  "not large enough: (",fdim,",",ldim,") not >= (", &
! Subprogram not used                   btofc_blk_offset(blockid)%ncols,",", &
! Subprogram not used                   btofc_blk_offset(blockid)%nlvls,")"
! Subprogram not used       call endrun()
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    do k=1,btofc_blk_offset(blockid)%nlvls
! Subprogram not used       do i=1,btofc_blk_offset(blockid)%ncols
! Subprogram not used          pter(i,k) = 1 + record_size* &
! Subprogram not used                      (btofc_blk_offset(blockid)%pter(i,k))
! Subprogram not used       enddo
! Subprogram not used       do i=btofc_blk_offset(blockid)%ncols+1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    do k=btofc_blk_offset(blockid)%nlvls+1,ldim
! Subprogram not used       do i=1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine block_to_chunk_send_pters
!
!========================================================================

! Subprogram not used    subroutine block_to_chunk_recv_pters(lcid, fdim, ldim, &
! Subprogram not used                                         record_size, pter)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return pointers into receive buffer where data for
! Subprogram not used !          decomposed chunk data structures should be copied from
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: lcid         ! local chunk id
! Subprogram not used    integer, intent(in) :: fdim         ! first dimension of pter array
! Subprogram not used    integer, intent(in) :: ldim         ! last dimension of pter array
! Subprogram not used    integer, intent(in) :: record_size  ! per coordinate amount of data 
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: pter(fdim,ldim)  ! buffer offset
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i, k                     ! loop indices
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if ((btofc_chk_offset(lcid)%ncols > fdim) .or. &
! Subprogram not used        (btofc_chk_offset(lcid)%nlvls > ldim)) then
! Subprogram not used       write(iulog,*) "BLOCK_TO_CHUNK_RECV_PTERS: pter array dimensions ", &
! Subprogram not used                  "not large enough: (",fdim,",",ldim,") not >= (", &
! Subprogram not used                   btofc_chk_offset(lcid)%ncols,",", &
! Subprogram not used                   btofc_chk_offset(lcid)%nlvls,")"
! Subprogram not used       call endrun()
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    do k=1,btofc_chk_offset(lcid)%nlvls
! Subprogram not used       do i=1,btofc_chk_offset(lcid)%ncols
! Subprogram not used          pter(i,k) = 1 + record_size* &
! Subprogram not used                      (btofc_chk_offset(lcid)%pter(i,k))
! Subprogram not used       enddo
! Subprogram not used       do i=btofc_chk_offset(lcid)%ncols+1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    do k=btofc_chk_offset(lcid)%nlvls+1,ldim
! Subprogram not used       do i=1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine block_to_chunk_recv_pters
!
!========================================================================

! Subprogram not used    subroutine transpose_chunk_to_block(record_size, chunk_buffer, &
! Subprogram not used                                        block_buffer, window)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Transpose buffer containing decomposed 
! Subprogram not used !          chunk data structures to buffer
! Subprogram not used !          containing decomposed fields 
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use spmd_utils,    only: altalltoallv
! Subprogram not used !------------------------------Parameters-------------------------------
! Subprogram not used !
! Subprogram not used   integer, parameter :: msgtag  = 7000
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: record_size  ! per column amount of data 
! Subprogram not used    real(r8), intent(in):: chunk_buffer(record_size*chunk_buf_nrecs)
! Subprogram not used                                        ! buffer of chunk data to be
! Subprogram not used                                        ! transposed
! Subprogram not used    real(r8), intent(out) :: block_buffer(record_size*block_buf_nrecs)
! Subprogram not used                                        ! buffer of block data to
! Subprogram not used                                        ! transpose into
! Subprogram not used    integer, intent(in), optional :: window
! Subprogram not used                                        ! MPI-2 window id for
! Subprogram not used                                        ! chunk_buffer
! Subprogram not used 
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i, p                     ! loop indices
! Subprogram not used    integer :: bbuf_siz                 ! size of block_buffer
! Subprogram not used    integer :: cbuf_siz                 ! size of chunk_buffer
! Subprogram not used    integer :: lwindow                  ! placeholder for missing window
! Subprogram not used    integer :: lopt                     ! local copy of phys_alltoall
! Subprogram not used !
! Subprogram not used    logical, save :: first = .true.
! Subprogram not used    integer, allocatable, save :: sndcnts(:), sdispls(:)
! Subprogram not used    integer, allocatable, save :: rcvcnts(:), rdispls(:)
! Subprogram not used    integer, allocatable, save :: pdispls(:)
! Subprogram not used    integer, save :: prev_record_size = 0
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if (first) then
! Subprogram not used ! Compute send/recv/put counts and displacements
! Subprogram not used       allocate(sndcnts(0:npes-1))
! Subprogram not used       allocate(sdispls(0:npes-1))
! Subprogram not used       allocate(rcvcnts(0:npes-1))
! Subprogram not used       allocate(rdispls(0:npes-1))
! Subprogram not used       allocate(pdispls(0:npes-1))
! Subprogram not used !
! Subprogram not used !
! Subprogram not used       first = .false.
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    if (record_size .ne. prev_record_size) then
! Subprogram not used !
! Subprogram not used ! Compute send/recv/put counts and displacements
! Subprogram not used       sdispls(0) = 0
! Subprogram not used       sndcnts(0) = record_size*btofc_chk_num(0)
! Subprogram not used       do p=1,npes-1
! Subprogram not used         sdispls(p) = sdispls(p-1) + sndcnts(p-1)
! Subprogram not used         sndcnts(p) = record_size*btofc_chk_num(p)
! Subprogram not used       enddo
! Subprogram not used !
! Subprogram not used       rdispls(0) = 0
! Subprogram not used       rcvcnts(0) = record_size*btofc_blk_num(0)
! Subprogram not used       do p=1,npes-1
! Subprogram not used          rdispls(p) = rdispls(p-1) + rcvcnts(p-1)
! Subprogram not used          rcvcnts(p) = record_size*btofc_blk_num(p)
! Subprogram not used       enddo
! Subprogram not used !
! Subprogram not used       call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used       prev_record_size = record_size
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    call t_barrierf('sync_tran_ctob', mpicom)
! Subprogram not used    if (phys_alltoall < 0) then
! Subprogram not used       if ((max_nproc_smpx > npes/2) .and. (nproc_busy_d > npes/2)) then
! Subprogram not used          lopt = 0
! Subprogram not used       else
! Subprogram not used          lopt = 1
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       lopt = phys_alltoall
! Subprogram not used       if ((lopt .eq. 2) .and. ( .not. present(window) )) lopt = 1
! Subprogram not used    endif
! Subprogram not used    if (lopt < 4) then
! Subprogram not used !
! Subprogram not used       bbuf_siz = record_size*block_buf_nrecs
! Subprogram not used       cbuf_siz = record_size*chunk_buf_nrecs
! Subprogram not used       if ( present(window) ) then
! Subprogram not used          call altalltoallv(lopt, iam, npes,    &
! Subprogram not used                            dp_coup_steps, dp_coup_proc, &
! Subprogram not used                            chunk_buffer, cbuf_siz, sndcnts, sdispls, mpir8, &
! Subprogram not used                            block_buffer, bbuf_siz, rcvcnts, rdispls, mpir8, &
! Subprogram not used                            msgtag, pdispls, mpir8, window, mpicom)
! Subprogram not used       else
! Subprogram not used          call altalltoallv(lopt, iam, npes,    &
! Subprogram not used                            dp_coup_steps, dp_coup_proc, &
! Subprogram not used                            chunk_buffer, cbuf_siz, sndcnts, sdispls, mpir8, &
! Subprogram not used                            block_buffer, bbuf_siz, rcvcnts, rdispls, mpir8, &
! Subprogram not used                            msgtag, pdispls, mpir8, lwindow, mpicom)
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used    else
! Subprogram not used       call mpialltoallv(chunk_buffer, sndcnts, sdispls, mpir8, &
! Subprogram not used                         block_buffer, rcvcnts, rdispls, mpir8, &
! Subprogram not used                         mpicom)
! Subprogram not used !
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine transpose_chunk_to_block
!
!========================================================================

! Subprogram not used    subroutine chunk_to_block_send_pters(lcid, fdim, ldim, &
! Subprogram not used                                         record_size, pter)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return pointers into send buffer where data for
! Subprogram not used !          decomposed chunk data structures should be copied to
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: lcid         ! local chunk id
! Subprogram not used    integer, intent(in) :: fdim         ! first dimension of pter array
! Subprogram not used    integer, intent(in) :: ldim         ! last dimension of pter array
! Subprogram not used    integer, intent(in) :: record_size  ! per coordinate amount of data 
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: pter(fdim,ldim)  ! buffer offset
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i, k                     ! loop indices
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if ((btofc_chk_offset(lcid)%ncols > fdim) .or. &
! Subprogram not used        (btofc_chk_offset(lcid)%nlvls > ldim)) then
! Subprogram not used       write(iulog,*) "CHUNK_TO_BLOCK_SEND_PTERS: pter array dimensions ", &
! Subprogram not used                  "not large enough: (",fdim,",",ldim,") not >= (", &
! Subprogram not used                   btofc_chk_offset(lcid)%ncols,",", &
! Subprogram not used                   btofc_chk_offset(lcid)%nlvls,")"
! Subprogram not used       call endrun()
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    do k=1,btofc_chk_offset(lcid)%nlvls
! Subprogram not used       do i=1,btofc_chk_offset(lcid)%ncols
! Subprogram not used          pter(i,k) = 1 + record_size* &
! Subprogram not used                      (btofc_chk_offset(lcid)%pter(i,k))
! Subprogram not used       enddo
! Subprogram not used       do i=btofc_chk_offset(lcid)%ncols+1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    do k=btofc_chk_offset(lcid)%nlvls+1,ldim
! Subprogram not used       do i=1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine chunk_to_block_send_pters
!
!========================================================================

! Subprogram not used    subroutine chunk_to_block_recv_pters(blockid, fdim, ldim, &
! Subprogram not used                                         record_size, pter)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Return pointers into receive buffer where column from decomposed 
! Subprogram not used !          fields should be copied from
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in) :: blockid      ! block index
! Subprogram not used    integer, intent(in) :: fdim         ! first dimension of pter array
! Subprogram not used    integer, intent(in) :: ldim         ! last dimension of pter array
! Subprogram not used    integer, intent(in) :: record_size  ! per coordinate amount of data 
! Subprogram not used 
! Subprogram not used    integer, intent(out) :: pter(fdim,ldim)  ! buffer offsets
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: i, k                     ! loop indices
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if ((btofc_blk_offset(blockid)%ncols > fdim) .or. &
! Subprogram not used        (btofc_blk_offset(blockid)%nlvls > ldim)) then
! Subprogram not used       write(iulog,*) "CHUNK_TO_BLOCK_RECV_PTERS: pter array dimensions ", &
! Subprogram not used                  "not large enough: (",fdim,",",ldim,") not >= (", &
! Subprogram not used                   btofc_blk_offset(blockid)%ncols,",", &
! Subprogram not used                   btofc_blk_offset(blockid)%nlvls,")"
! Subprogram not used       call endrun()
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    do k=1,btofc_blk_offset(blockid)%nlvls
! Subprogram not used       do i=1,btofc_blk_offset(blockid)%ncols
! Subprogram not used          pter(i,k) = 1 + record_size* &
! Subprogram not used                      (btofc_blk_offset(blockid)%pter(i,k))
! Subprogram not used       enddo
! Subprogram not used       do i=btofc_blk_offset(blockid)%ncols+1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    do k=btofc_blk_offset(blockid)%nlvls+1,ldim
! Subprogram not used       do i=1,fdim
! Subprogram not used          pter(i,k) = -1
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine chunk_to_block_recv_pters
!
!========================================================================

   subroutine create_chunks(opt, chunks_per_thread)
!----------------------------------------------------------------------- 
! 
! Purpose: Decompose physics computational grid into chunks, for
!          improved serial efficiency and parallel load balance.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plev
   use dyn_grid, only: get_block_bounds_d, get_block_gcol_cnt_d, &
                       get_gcol_block_cnt_d, get_gcol_block_d, &
                       get_block_owner_d, get_block_gcol_d
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: opt           ! chunking option
      !  0: chunks may cross block boundaries, but retain same
      !     process mapping as blocks. If possible, columns assigned
      !     as day/night pairs. Columns (or pairs) are wrap-mapped.
      !     May not work with vertically decomposed blocks. (default)
      !  1: chunks may cross block boundaries, but retain same
      !     SMP-node mapping as blocks.  If possible, columns assigned
      !     as day/night pairs.  Columns (or pairs) are wrap-mapped.
      !     May not work with vertically decomposed blocks.
      !  2: 2-column day/night and season column pairs wrap-mapped
      !     to chunks to also balance assignment of polar, mid-latitude, 
      !     and equatorial columns across  chunks.
      !  3: same as 1 except that SMP defined to be pairs of consecutive
      !     processes
      !  4: chunks may cross block boundaries, but retain same
      !     process mapping as blocks. Columns assigned to chunks
      !     in block ordering.
      !     May not work with vertically decomposed blocks.
      !  5: Chunks do not cross latitude boundaries, and are block-mapped.
   integer, intent(in)  :: chunks_per_thread 
                                         ! target number of chunks per
                                         !  thread
!---------------------------Local workspace-----------------------------
   integer :: i, j, p                    ! loop indices
   integer :: nlthreads                  ! number of local OpenMP threads
   integer :: npthreads(0:npes-1)        ! number of OpenMP threads per process
   integer :: proc_smp_mapx(0:npes-1)    ! process/virtual SMP node map
   integer :: firstblock, lastblock      ! global block index bounds
   integer :: maxblksiz                  ! maximum number of columns in a dynamics block
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: nsmpx, nsmpy               ! virtual SMP node counts and indices
   integer :: curgcol, twingcol          ! global physics and dynamics column indices
   integer :: smp                        ! SMP node index
   integer :: cid                        ! chunk id
   integer :: jb, ib                     ! global block and columns indices
   integer :: blksiz                     ! current block size
   integer :: ntmp1, ntmp2, nlchunks     ! work variables
   integer :: cbeg                       ! beginning longitude index for 
                                         !  current chunk
   integer :: max_ncols                  ! upper bound on number of columns in a block
   integer :: ncols                      ! number of columns in current chunk
   logical :: error                      ! error flag 

   ! indices for dynamics columns in given block
   integer, dimension(:), allocatable :: cols

   ! number of MPI processes per virtual SMP node (0:nsmpx-1)
   integer, dimension(:), allocatable :: nsmpprocs      

   ! flag indicating whether a process is busy or idle during the dynamics (0:npes-1)
   logical, dimension(:), allocatable :: proc_busy_d

   ! flag indicating whether any of the processes assigned to an SMP node are busy 
   ! during the dynamics, or whether all of them are idle (0:nsmps-1)
   logical, dimension(:), allocatable :: smp_busy_d

   ! actual SMP node/virtual SMP node map (0:nsmps-1)    
   integer, dimension(:), allocatable :: smp_smp_mapx

   ! column/virtual SMP node map (ngcols)
   integer, dimension(:), allocatable :: col_smp_mapx

   ! number of columns assigned to a given virtual SMP node (0:nsmpx-1)
   integer, dimension(:), allocatable :: nsmpcolumns

   ! number of OpenMP threads per virtual SMP node (0:nsmpx-1)
   integer, dimension(:), allocatable :: nsmpthreads

   ! number of chunks assigned to a given virtual SMP node (0:nsmpx-1)
   integer, dimension(:), allocatable :: nsmpchunks
                                         
   ! maximum number of columns assigned to a chunk in a given virtual SMP node (0:nsmpx-1)
   integer, dimension(:), allocatable :: maxcol_chk
                                         
   ! number of chunks in given virtual SMP node receiving maximum number of columns 
   ! (0:nsmpx-1)
   integer, dimension(:), allocatable :: maxcol_chks

   ! chunk id virtual offset (0:nsmpx-1)
   integer, dimension(:), allocatable :: cid_offset

   ! process-local chunk id (0:nsmpx-1)
   integer, dimension(:), allocatable :: local_cid


!-----------------------------------------------------------------------
!
! Determine number of threads per process
!
   nlthreads = 1
!
   call mpiallgatherint(nlthreads, 1, npthreads, 1, mpicom)

!
! Determine index range for dynamics blocks
!
   call get_block_bounds_d(firstblock,lastblock)

!
! Determine maximum number of columns in a block
!
   maxblksiz = 0
   do jb=firstblock,lastblock
      maxblksiz = max(maxblksiz,get_block_gcol_cnt_d(jb))
   enddo

!
!  determine which (and how many) processes are assigned
!  dynamics blocks
!
   allocate( proc_busy_d(0:npes-1) )
   proc_busy_d = .false.
   nproc_busy_d = 0
   do jb=firstblock,lastblock
      p = get_block_owner_d(jb)
      if (.not. proc_busy_d(p) ) then
         proc_busy_d(p) = .true.
         nproc_busy_d = nproc_busy_d + 1
      endif
   enddo

!
! Determine virtual SMP count and processes/virtual SMP map.
!  If option 0 or >3, pretend that each SMP has only one process. 
!  If option 1, use SMP information.
!  If option 2, pretend that all processes are in one SMP node. 
!  If option 3, pretend that each SMP node is made up of two
!     processes, chosen to maximize load-balancing opportunities.
!
!  For all options < 5, if there are "idle" dynamics processes, 
!     assign them to the virtual SMP nodes in wrap fashion.
!     Communication between the active and idle dynamics 
!     processes is scatter/gather (no communications between 
!     idle dynamics processes) so there is no advantage to 
!     blocking the idle processes in these assignments.
!
   if ((opt <= 0) .or. (opt == 4)) then

!     assign active dynamics processes to virtual SMP nodes
      nsmpx = 0
      do p=0,npes-1
         if (proc_busy_d(p)) then
            proc_smp_mapx(p) = nsmpx
            nsmpx = nsmpx + 1
         endif
      enddo
! 
!     assign idle dynamics processes to virtual SMP nodes (wrap map)
      nsmpy = 0
      do p=0,npes-1
         if (.not. proc_busy_d(p)) then
            proc_smp_mapx(p) = nsmpy
            nsmpy = mod(nsmpy+1,nsmpx)
         endif
      enddo

   elseif (opt == 1) then

      allocate( smp_busy_d(0:nsmps-1) )
      allocate( smp_smp_mapx(0:nsmps-1) )

!
!     determine SMP nodes assigned dynamics blocks
      smp_busy_d = .false.
      do p=0,npes-1
         if ( proc_busy_d(p) ) then
            smp = proc_smp_map(p)
            smp_busy_d(smp) = .true.
         endif
      enddo

!
!     determine number of SMP nodes assigned dynamics blocks
      nsmpx = 0
      do smp=0,nsmps-1
         if (smp_busy_d(smp)) then
            smp_smp_mapx(smp) = nsmpx
            nsmpx = nsmpx + 1
         endif
      enddo
!
!     assign processes in active dynamics SMP nodes to virtual SMP nodes
      do p=0,npes-1
         smp = proc_smp_map(p)
         if (smp_busy_d(smp)) then
            proc_smp_mapx(p) = smp_smp_mapx(smp)
         endif
      enddo
! 
!     assign processes in idle dynamics SMP nodes to virtual SMP nodes (wrap map)
      nsmpy = 0
      do p=0,npes-1
         smp = proc_smp_map(p)
         if (.not. smp_busy_d(smp)) then
            proc_smp_mapx(p) = nsmpy
            nsmpy = mod(nsmpy+1,nsmpx)
         endif
      enddo
!
      deallocate( smp_busy_d )
      deallocate( smp_smp_mapx )

   elseif (opt == 2) then

      nsmpx = 1
      do p=0,npes-1
         proc_smp_mapx(p) = 0
      enddo

   elseif (opt == 3) then

!     find active process partners
      proc_smp_mapx = -1
      call find_partners(opt,proc_busy_d,nsmpx,proc_smp_mapx)
! 
!     assign unassigned (idle dynamics) processes to virtual SMP nodes 
!     (wrap map)
      nsmpy = 0
      do p=0,npes-1
         if (proc_smp_mapx(p) .eq. -1) then
            proc_smp_mapx(p) = nsmpy
            nsmpy = mod(nsmpy+1,nsmpx)
         endif
      enddo

   else

      nsmpx = npes
      do p=0,npes-1
         proc_smp_mapx(p) = p
      enddo

   endif
!
   deallocate( proc_busy_d )

!
! Determine maximum number of processes assigned to a single 
! virtual SMP node
!
   allocate( nsmpprocs(0:nsmpx-1) )
!
   nsmpprocs(:) = 0
   do p=0,npes-1
      smp = proc_smp_mapx(p)
      nsmpprocs(smp) = nsmpprocs(smp) + 1
   enddo
   max_nproc_smpx = maxval(nsmpprocs)
!
   deallocate( nsmpprocs )

!
! Determine number of columns assigned to each
! virtual SMP in block decomposition

   allocate( col_smp_mapx(ngcols) )
!
   col_smp_mapx(:) = -1
   error = .false.
   do i=1,ngcols_p
      curgcol = latlon_to_dyn_gcol_map(i)
      block_cnt = get_gcol_block_cnt_d(curgcol)
      call get_gcol_block_d(curgcol,block_cnt,blockids,bcids)
      do jb=1,block_cnt
         p = get_block_owner_d(blockids(jb)) 
         if (col_smp_mapx(i) .eq. -1) then
            col_smp_mapx(i) = proc_smp_mapx(p)
         elseif (col_smp_mapx(i) .ne. proc_smp_mapx(p)) then
            error = .true.
         endif
      enddo
   end do
   if (error) then
      write(iulog,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
               "but vertical decomposition not limited to virtual SMP"
      call endrun()
   endif
!
   allocate( nsmpcolumns(0:nsmpx-1) )
   nsmpcolumns(:) = 0
   do i=1,ngcols_p
      curgcol = latlon_to_dyn_gcol_map(i)
      smp = col_smp_mapx(curgcol)
      nsmpcolumns(smp) = nsmpcolumns(smp) + 1
   end do
!
   deallocate( col_smp_mapx )

!
!  Allocate other work space
!
   allocate( nsmpthreads(0:nsmpx-1) )
   allocate( nsmpchunks (0:nsmpx-1) )
   allocate( maxcol_chk (0:nsmpx-1) )
   allocate( maxcol_chks(0:nsmpx-1) )
   allocate( cid_offset (0:nsmpx-1) )
   allocate( local_cid  (0:nsmpx-1) )
   allocate( cols(1:maxblksiz) )
!
! Options 0-3: split local dynamics blocks into chunks,
!              using wrap-map assignment of columns and
!              day/night and north/south column pairs
!              to chunks to improve load balance
!  Option 0: local is per process
!  Option 1: local is subset of`processes assigned to same SMP node
!  Option 2: local is global
!  Option 3: local is pair of processes chosen to maximize load-balance
!            wrt restriction that only communicate with one other
!            process.
! Option 4: split local dynamics blocks into chunks,
!           using block-map assignment of columns
!             
   if ((opt >= 0) .and. (opt <= 4)) then
!
! Calculate number of threads available in each SMP node. 
!
      nsmpthreads(:) = 0
      do p=0,npes-1
         smp = proc_smp_mapx(p)
         nsmpthreads(smp) = nsmpthreads(smp) + npthreads(p)
      enddo
!
! Determine number of chunks to keep all threads busy
!
      nchunks = 0
      do smp=0,nsmpx-1
         nsmpchunks(smp) = nsmpcolumns(smp)/pcols
         if (mod(nsmpcolumns(smp), pcols) .ne. 0) then
            nsmpchunks(smp) = nsmpchunks(smp) + 1
         endif
         if (nsmpchunks(smp) < chunks_per_thread*nsmpthreads(smp)) then
            nsmpchunks(smp) = chunks_per_thread*nsmpthreads(smp)
         endif
         do while (mod(nsmpchunks(smp), nsmpthreads(smp)) .ne. 0)
            nsmpchunks(smp) = nsmpchunks(smp) + 1
         enddo
         if (nsmpchunks(smp) > nsmpcolumns(smp)) then
            nsmpchunks(smp) = nsmpcolumns(smp)
         endif
         nchunks = nchunks + nsmpchunks(smp)
      enddo
!
! Determine maximum number of columns to assign to chunks
! in a given SMP
!
      do smp=0,nsmpx-1
         if (nsmpchunks(smp) /= 0) then
            ntmp1 = nsmpcolumns(smp)/nsmpchunks(smp)
            ntmp2 = mod(nsmpcolumns(smp),nsmpchunks(smp))
            if (ntmp2 > 0) then
               maxcol_chk(smp) = ntmp1 + 1
               maxcol_chks(smp) = ntmp2
            else
               maxcol_chk(smp) = ntmp1
               maxcol_chks(smp) = nsmpchunks(smp)
            endif
         else
            maxcol_chk(smp) = 0
            maxcol_chks(smp) = 0
         endif
      enddo
!
! Allocate chunks and knuhcs data structures
!
      allocate( chunks(1:nchunks) )
      allocate( knuhcs(1:ngcols) )
!
! Initialize chunks and knuhcs data structures
!
      chunks(:)%ncols = 0
      knuhcs(:)%chunkid = -1
      knuhcs(:)%col = -1
!
! Determine chunk id ranges for each SMP
!
      cid_offset(0) = 1
      local_cid(0) = 0
      do smp=1,nsmpx-1
         cid_offset(smp) = cid_offset(smp-1) + nsmpchunks(smp-1)
         local_cid(smp) = 0
      enddo
!
! Assign columns to chunks
!
      do jb=firstblock,lastblock
         p = get_block_owner_d(jb)
         smp = proc_smp_mapx(p)
         blksiz = get_block_gcol_cnt_d(jb)
         call get_block_gcol_d(jb,blksiz,cols)
         do ib = 1,blksiz
!
! Assign column to a chunk if not already assigned
            curgcol = cols(ib)
            if ((dyn_to_latlon_gcol_map(curgcol) .ne. -1) .and. &
                (knuhcs(curgcol)%chunkid == -1)) then
!
! Find next chunk with space
! (maxcol_chks > 0 test necessary for opt=4 block map)
               cid = cid_offset(smp) + local_cid(smp)
               if (maxcol_chks(smp) > 0) then
                  do while (chunks(cid)%ncols >=  maxcol_chk(smp))
                     local_cid(smp) = mod(local_cid(smp)+1,nsmpchunks(smp))
                     cid = cid_offset(smp) + local_cid(smp)
                  enddo
               else
                  do while (chunks(cid)%ncols >=  maxcol_chk(smp)-1)
                     local_cid(smp) = mod(local_cid(smp)+1,nsmpchunks(smp))
                     cid = cid_offset(smp) + local_cid(smp)
                  enddo
               endif
               chunks(cid)%ncols = chunks(cid)%ncols + 1
               if (chunks(cid)%ncols .eq. maxcol_chk(smp)) &
                  maxcol_chks(smp) = maxcol_chks(smp) - 1
!
               i = chunks(cid)%ncols
               chunks(cid)%gcol(i) = curgcol
               chunks(cid)%lon(i)  = lon_p(curgcol)
               chunks(cid)%lat(i)  = lat_p(curgcol)
               knuhcs(curgcol)%chunkid = cid
               knuhcs(curgcol)%col = i
!
               if (opt < 4) then
!
! If space available, look to assign a load-balancing "twin" to same chunk
                  if ( (chunks(cid)%ncols <  maxcol_chk(smp)) .and. &
                       (maxcol_chks(smp) > 0) .and. (twin_alg > 0)) then

                     call find_twin(curgcol, smp, &
                                    proc_smp_mapx, twingcol)

                     if (twingcol > 0) then
                        chunks(cid)%ncols = chunks(cid)%ncols + 1
                        if (chunks(cid)%ncols .eq. maxcol_chk(smp)) &
                           maxcol_chks(smp) = maxcol_chks(smp) - 1
!
                        i = chunks(cid)%ncols
                        chunks(cid)%gcol(i) = twingcol
                        chunks(cid)%lon(i) = lon_p(twingcol)
                        chunks(cid)%lat(i) = lat_p(twingcol)
                        knuhcs(twingcol)%chunkid = cid
                        knuhcs(twingcol)%col = i
                     endif
!
                  endif
!
! Move on to next chunk (wrap map)
                  local_cid(smp) = mod(local_cid(smp)+1,nsmpchunks(smp))
!
               endif
!
            endif
         enddo
      enddo
!
   else
!
! Option 5: split individual dynamics blocks into chunks,
!            assigning consecutive columns to the same chunk
!
! Determine total number of chunks and
! number of chunks in each "SMP node"
!  (assuming no vertical decomposition)
      nchunks = 0
      nsmpchunks(:) = 0
      do j=firstblock,lastblock
         blksiz = get_block_gcol_cnt_d(j)
         nlchunks = blksiz/pcols
         if (pcols*(blksiz/pcols) /= blksiz) then
            nlchunks = nlchunks + 1
         endif
         nchunks = nchunks + nlchunks
         p = get_block_owner_d(j) 
         nsmpchunks(p) = nsmpchunks(p) + nlchunks
      enddo
!
! Determine chunk id ranges for each SMP
!
      cid_offset(0) = 1
      local_cid(0) = 0
      do smp=1,nsmpx-1
         cid_offset(smp) = cid_offset(smp-1) + nsmpchunks(smp-1)
         local_cid(smp) = 0
      enddo
!
! Allocate chunks and knuhcs data structures
!
      allocate( chunks(1:nchunks) )
      allocate( knuhcs(1:ngcols) )
!
! Initialize chunks and knuhcs data structures
!
      knuhcs(:)%chunkid = -1
      knuhcs(:)%col = -1
      cid = 0
      do jb=firstblock,lastblock
         p = get_block_owner_d(jb)
         smp = proc_smp_mapx(p)
         blksiz = get_block_gcol_cnt_d(jb)
         call get_block_gcol_d(jb,blksiz,cols)

         ib = 0
         do while (ib < blksiz)

            cid = cid_offset(smp) + local_cid(smp)
            max_ncols = min(pcols,blksiz-ib)

            ncols = 0
            do i=1,max_ncols
               ib = ib + 1
               ! check whether global index is for a column that dynamics
               ! intends to pass to the physics
               curgcol = cols(ib)
               if (dyn_to_latlon_gcol_map(curgcol) .ne. -1) then
                  ! yes - then save the information
                  ncols = ncols + 1
                  chunks(cid)%gcol(ncols) = curgcol
                  chunks(cid)%lon(ncols)  = lon_p(curgcol)
                  chunks(cid)%lat(ncols)  = lat_p(curgcol)
                  knuhcs(curgcol)%chunkid = cid
                  knuhcs(curgcol)%col = ncols
               endif
            enddo
            chunks(cid)%ncols = ncols

            local_cid(smp) = local_cid(smp) + 1
         enddo
      enddo
!
! Set number of threads available in each "SMP node". 
!
      do p=0,npes-1
         nsmpthreads(p) = npthreads(p)
      enddo
!
   endif
!
! Assign chunks to processes.
!
   call assign_chunks(npthreads, nsmpx, proc_smp_mapx, &
                      nsmpthreads, nsmpchunks)
!
! Clean up
!
   deallocate( nsmpcolumns )
   deallocate( nsmpthreads )
   deallocate( nsmpchunks  )
   deallocate( maxcol_chk  )
   deallocate( maxcol_chks )
   deallocate( cid_offset  )
   deallocate( local_cid   )
   deallocate( cols )
   deallocate( knuhcs )

   return
   end subroutine create_chunks
!
!========================================================================

! Subprogram not used    subroutine find_partners(opt, proc_busy_d, nsmpx, proc_smp_mapx)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Divide processes into pairs, attempting to maximize the
! Subprogram not used !          the number of columns in one process whose twins are in the 
! Subprogram not used !          other process.
! Subprogram not used ! 
! Subprogram not used ! Method: The day/night and north/south hemisphere complement is defined
! Subprogram not used !         to be the column twin.
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use dyn_grid, only: get_gcol_block_cnt_d, get_gcol_block_d, &
! Subprogram not used                        get_block_owner_d
! Subprogram not used    use pmgrid, only: plev
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: opt           ! chunking option
! Subprogram not used    logical, intent(in)  :: proc_busy_d(0:npes-1)
! Subprogram not used                                          ! active/idle dynamics process flags
! Subprogram not used    integer, intent(out) :: nsmpx         ! calculated number of virtual 
! Subprogram not used                                          !  SMP nodes
! Subprogram not used    integer, intent(out) :: proc_smp_mapx(0:npes-1)
! Subprogram not used                                          ! process/virtual smp map
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: gcol_latlon                ! physics column index (latlon sorted)
! Subprogram not used    integer :: twingcol_latlon            ! physics column index (latlon sorted)
! Subprogram not used    integer :: gcol, twingcol             ! physics column indices
! Subprogram not used    integer :: lon, lat, twinlat          ! longitude and latitude indices
! Subprogram not used    integer :: twinlon_off                ! estimate as to offset of twinlon
! Subprogram not used                                          ! on a latitude line
! Subprogram not used    integer :: block_cnt                  ! number of blocks containing data
! Subprogram not used                                          ! for a given vertical column
! Subprogram not used    integer :: blockids(plev+1)           ! block indices
! Subprogram not used    integer :: bcids(plev+1)              ! block column indices
! Subprogram not used    integer :: jb                         ! block index
! Subprogram not used    integer :: p, twp                     ! process indices
! Subprogram not used    integer :: col_proc_mapx(ngcols)      ! location of columns in 
! Subprogram not used                                          !  dynamics decomposition
! Subprogram not used    integer :: twin_proc_mapx(ngcols)     ! location of column twins in 
! Subprogram not used                                          !  dynamics decomposition
! Subprogram not used    integer :: twin_cnt(0:npes-1)         ! for each process, number of twins 
! Subprogram not used                                          !  in each of the other processes
! Subprogram not used    logical :: assigned(0:npes-1)         ! flag indicating whether process
! Subprogram not used                                          !  assigned to an SMP node yet
! Subprogram not used    integer :: maxpartner, maxcnt         ! process with maximum number of 
! Subprogram not used                                          !  twins and this count
! Subprogram not used 
! Subprogram not used    logical :: error                      ! error flag 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Determine process location of column and its twin in dynamics decomposition
! Subprogram not used !
! Subprogram not used    col_proc_mapx(:) = -1
! Subprogram not used    twin_proc_mapx(:) = -1
! Subprogram not used 
! Subprogram not used    error = .false.
! Subprogram not used    do gcol_latlon=1,ngcols_p
! Subprogram not used 
! Subprogram not used       ! Assume latitude and longitude symmetries and that index manipulations
! Subprogram not used       ! are sufficient to find partners. (Will be true for lon/lat grids.)
! Subprogram not used       gcol = latlon_to_dyn_gcol_map(gcol_latlon)
! Subprogram not used       lat = lat_p(gcol)
! Subprogram not used       twinlat = clat_p_tot+1-lat
! Subprogram not used       lon = lon_p(gcol)
! Subprogram not used       twinlon_off = mod((lon-1)+(clat_p_cnt(twinlat)/2), clat_p_cnt(twinlat))
! Subprogram not used       twingcol_latlon = clat_p_idx(twinlat) + twinlon_off
! Subprogram not used       twingcol = latlon_to_dyn_gcol_map(twingcol_latlon)
! Subprogram not used 
! Subprogram not used       block_cnt = get_gcol_block_cnt_d(gcol)
! Subprogram not used       call get_gcol_block_d(gcol,block_cnt,blockids,bcids)
! Subprogram not used       do jb=1,block_cnt
! Subprogram not used          p = get_block_owner_d(blockids(jb)) 
! Subprogram not used          if (col_proc_mapx(gcol) .eq. -1) then
! Subprogram not used             col_proc_mapx(gcol) = p
! Subprogram not used          elseif (col_proc_mapx(gcol) .ne. p) then
! Subprogram not used             error = .true.
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       block_cnt = get_gcol_block_cnt_d(twingcol)
! Subprogram not used       call get_gcol_block_d(twingcol,block_cnt,blockids,bcids)
! Subprogram not used       do jb=1,block_cnt
! Subprogram not used          p = get_block_owner_d(blockids(jb)) 
! Subprogram not used          if (twin_proc_mapx(gcol) .eq. -1) then
! Subprogram not used             twin_proc_mapx(gcol) = p
! Subprogram not used          elseif (twin_proc_mapx(gcol) .ne. p) then
! Subprogram not used             error = .true.
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    if (error) then
! Subprogram not used       if (masterproc) then
! Subprogram not used          write(iulog,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
! Subprogram not used             "but vertical decomposition not limited to single process"
! Subprogram not used       endif
! Subprogram not used       call endrun()
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Assign process pairs to SMPs, attempting to maximize the number of column,twin
! Subprogram not used ! pairs in same SMP.
! Subprogram not used !
! Subprogram not used    assigned(:) = .false.
! Subprogram not used    twin_cnt(:) = 0
! Subprogram not used    nsmpx = 0
! Subprogram not used    do p=0,npes-1
! Subprogram not used       if ((.not. assigned(p)) .and. (proc_busy_d(p))) then
! Subprogram not used !
! Subprogram not used ! For each process, determine number of twins in each of the other processes
! Subprogram not used ! (running over all columns multiple times to minimize memory requirements).
! Subprogram not used !
! Subprogram not used          do gcol_latlon=1,ngcols_p
! Subprogram not used             gcol = latlon_to_dyn_gcol_map(gcol_latlon)
! Subprogram not used             if (col_proc_mapx(gcol) .eq. p) then
! Subprogram not used                twin_cnt(twin_proc_mapx(gcol)) = &
! Subprogram not used                   twin_cnt(twin_proc_mapx(gcol)) + 1
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used !
! Subprogram not used ! Find process with maximum number of twins that has not yet been designated
! Subprogram not used ! a partner.
! Subprogram not used !
! Subprogram not used          maxpartner = -1
! Subprogram not used          maxcnt = 0
! Subprogram not used          do twp=0,npes-1
! Subprogram not used             if ((.not. assigned(twp)) .and. (twp .ne. p)) then
! Subprogram not used                if (twin_cnt(twp) >= maxcnt) then
! Subprogram not used                   maxcnt = twin_cnt(twp)
! Subprogram not used                   maxpartner = twp
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used !
! Subprogram not used ! Assign p and twp to the same SMP node
! Subprogram not used !
! Subprogram not used          if (maxpartner .ne. -1) then
! Subprogram not used             assigned(p) = .true.
! Subprogram not used             assigned(maxpartner) = .true.
! Subprogram not used             proc_smp_mapx(p) = nsmpx
! Subprogram not used             proc_smp_mapx(maxpartner) = nsmpx
! Subprogram not used             nsmpx = nsmpx + 1
! Subprogram not used          else
! Subprogram not used             if (masterproc) then
! Subprogram not used                write(iulog,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
! Subprogram not used                   "but could not divide processes into pairs."
! Subprogram not used             endif
! Subprogram not used             call endrun()
! Subprogram not used          endif
! Subprogram not used !
! Subprogram not used       endif
! Subprogram not used !      
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine find_partners
!
!========================================================================

! Subprogram not used    subroutine find_twin(gcol, smp, proc_smp_mapx, twingcol_f)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Find column that when paired with gcol in a chunk
! Subprogram not used !          balances the load. A column is a candidate to be paired with
! Subprogram not used !          gcol if it is in the same SMP node as gcol as defined
! Subprogram not used !          by proc_smp_mapx.
! Subprogram not used ! 
! Subprogram not used ! Method: The day/night and north/south hemisphere complement is
! Subprogram not used !         tried first. If it is not a candidate or if it has already been
! Subprogram not used !         assigned, then the day/night complement is tried next. If that
! Subprogram not used !         also is not available, then nothing is returned.
! Subprogram not used ! 
! Subprogram not used ! Author: Patrick Worley
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    use dyn_grid, only: get_gcol_block_d, get_block_owner_d
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used    integer, intent(in)  :: gcol          ! global column index for column
! Subprogram not used                                          ! seeking a twin for
! Subprogram not used    integer, intent(in)  :: smp           ! index of SMP node 
! Subprogram not used                                          ! currently assigned to
! Subprogram not used    integer, intent(in)  :: proc_smp_mapx(0:npes-1)
! Subprogram not used                                          ! process/virtual smp map
! Subprogram not used    integer, intent(out) :: twingcol_f
! Subprogram not used                                          ! global column index for twin
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used    integer :: lon, lat                   ! global lon/lat indices for column
! Subprogram not used                                          ! seeking a twin for
! Subprogram not used    integer :: twinlon, twinlat           ! lon/lat indices of twin candidate
! Subprogram not used    integer :: twinlon_off                ! estimate as to offset of twinlon
! Subprogram not used                                          ! on a latitude line
! Subprogram not used    logical :: found                      ! found flag
! Subprogram not used    integer :: i                          ! loop index
! Subprogram not used    integer :: upper, lower               ! search temporaries
! Subprogram not used    integer :: twingcol_latlon            ! global physics column index (latlon sorted)
! Subprogram not used    integer :: twingcol_lonlat            ! global physics column index (lonlat sorted)
! Subprogram not used    integer :: twingcol                   ! global physics column indes
! Subprogram not used    integer :: diff, min_diff, min_i      ! search temporaries
! Subprogram not used    integer :: jbtwin(npes)               ! global block indices
! Subprogram not used    integer :: ibtwin(npes)               ! global column indices
! Subprogram not used    integer :: twinproc, twinsmp          ! process and smp ids
! Subprogram not used 
! Subprogram not used    integer :: clon_p_idx(clon_p_tot)     ! index in lonlat ordering for first 
! Subprogram not used                                          !  occurrence of longitude corresponding to 
! Subprogram not used                                          !  given latitude index
! Subprogram not used 
! Subprogram not used    real(r8):: twopi                      ! 2*pi
! Subprogram not used    real(r8):: clat, twinclat             ! latitude and twin
! Subprogram not used    real(r8):: clon, twinclon             ! longitude and twin
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    twingcol_f = -1
! Subprogram not used 
! Subprogram not used    ! precompute clon_p_idx
! Subprogram not used    clon_p_idx(1) = 1
! Subprogram not used    do i=2,clon_p_tot
! Subprogram not used       clon_p_idx(i) = clon_p_idx(i-1) + clon_p_cnt(i-1)
! Subprogram not used    enddo
! Subprogram not used !
! Subprogram not used ! Try day/night and north/south hemisphere complement first
! Subprogram not used !
! Subprogram not used    ! determine twin latitude
! Subprogram not used    lat = lat_p(gcol)
! Subprogram not used    clat = clat_p(lat)
! Subprogram not used    twinclat = -clat
! Subprogram not used    twinlat = clat_p_tot+1-lat
! Subprogram not used    if (clat_p(twinlat) .eq. twinclat) then
! Subprogram not used       found = .true.
! Subprogram not used    else
! Subprogram not used       found = .false.
! Subprogram not used       upper = twinlat
! Subprogram not used       lower = twinlat
! Subprogram not used       if (upper < clat_p_tot) upper = twinlat + 1
! Subprogram not used       if (lower > 1) lower = twinlat - 1
! Subprogram not used    endif
! Subprogram not used    do while (.not. found)
! Subprogram not used       if      ((abs(clat_p(upper)-twinclat) < abs(clat_p(twinlat)-twinclat)) .and. &
! Subprogram not used                (upper .ne. twinlat)) then
! Subprogram not used          twinlat = upper
! Subprogram not used          if (upper < clat_p_tot) then
! Subprogram not used             upper = twinlat + 1
! Subprogram not used          else
! Subprogram not used             found = .true.
! Subprogram not used          endif
! Subprogram not used       else if ((abs(clat_p(lower)-twinclat) < abs(clat_p(twinlat)-twinclat)) .and. &
! Subprogram not used                (lower .ne. twinlat))    then
! Subprogram not used          twinlat = lower
! Subprogram not used          if (lower > 1) then
! Subprogram not used             lower = twinlat - 1
! Subprogram not used          else
! Subprogram not used             found = .true.
! Subprogram not used          endif
! Subprogram not used       else
! Subprogram not used          found = .true.
! Subprogram not used       endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used    ! determine twin longitude
! Subprogram not used    twopi = 2.0_r8*pi
! Subprogram not used    lon = lon_p(gcol)
! Subprogram not used    clon = clon_p(lon)
! Subprogram not used    twinclon = mod(clon+pi,twopi)
! Subprogram not used    twinlon = mod((lon-1)+(clon_p_tot/2), clon_p_tot) + 1
! Subprogram not used    if (clon_p(twinlon) .eq. twinclon) then
! Subprogram not used       found = .true.
! Subprogram not used    else
! Subprogram not used       found = .false.
! Subprogram not used       upper = twinlon
! Subprogram not used       lower = twinlon
! Subprogram not used       if (upper < clon_p_tot) upper = twinlon + 1
! Subprogram not used       if (lower > 1) lower = twinlon - 1
! Subprogram not used    endif
! Subprogram not used    do while (.not. found)
! Subprogram not used       if      ((abs(clon_p(upper)-twinclon) < abs(clon_p(twinlon)-twinclon)) .and. &
! Subprogram not used                (upper .ne. twinlon)) then
! Subprogram not used          twinlon = upper
! Subprogram not used          if (upper < clon_p_tot) then
! Subprogram not used             upper = twinlon + 1
! Subprogram not used          else
! Subprogram not used             found = .true.
! Subprogram not used          endif
! Subprogram not used       else if ((abs(clon_p(lower)-twinclon) < abs(clon_p(twinlon)-twinclon)) .and. &
! Subprogram not used                (lower .ne. twinlon))    then
! Subprogram not used          twinlon = lower
! Subprogram not used          if (lower > 1) then
! Subprogram not used             lower = twinlon - 1
! Subprogram not used          else
! Subprogram not used             found = .true.
! Subprogram not used          endif
! Subprogram not used       else
! Subprogram not used          found = .true.
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    ! first, look for an exact match (assuming latitude and longitude symmetries)
! Subprogram not used    twinlon_off = mod((lon-1)+(clat_p_cnt(twinlat)/2), clat_p_cnt(twinlat))
! Subprogram not used    twingcol_latlon = clat_p_idx(twinlat) + twinlon_off
! Subprogram not used    twingcol = latlon_to_dyn_gcol_map(twingcol_latlon)
! Subprogram not used 
! Subprogram not used    ! otherwise, look around for an approximate match using lonlat sorted indices
! Subprogram not used    if ((lon_p(twingcol) .ne. twinlon) .or. (lat_p(twingcol) .ne. twinlat)) then
! Subprogram not used       twingcol_lonlat = clon_p_idx(twinlon)
! Subprogram not used       twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
! Subprogram not used       min_diff = abs(lat_p(twingcol) - twinlat)
! Subprogram not used       min_i = 0
! Subprogram not used       do i = 1, clon_p_cnt(twinlon)-1
! Subprogram not used          twingcol_lonlat = clon_p_idx(twinlon)+i
! Subprogram not used          twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
! Subprogram not used          diff = abs(lat_p(twingcol) - twinlat)
! Subprogram not used          if (diff < min_diff) then
! Subprogram not used             min_diff = diff
! Subprogram not used             min_i = i
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       twingcol_lonlat = clon_p_idx(twinlon) + min_i
! Subprogram not used       twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    ! Check whether twin and original are in same smp
! Subprogram not used    found = .false.
! Subprogram not used    call get_gcol_block_d(twingcol,npes,jbtwin,ibtwin)
! Subprogram not used    twinproc = get_block_owner_d(jbtwin(1))
! Subprogram not used    twinsmp  = proc_smp_mapx(twinproc)
! Subprogram not used !
! Subprogram not used    if ((twinsmp .eq. smp) .and. &
! Subprogram not used        (knuhcs(twingcol)%chunkid == -1)) then
! Subprogram not used       found = .true.
! Subprogram not used       twingcol_f = twingcol
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used ! Try day/night complement next
! Subprogram not used    if (.not. found) then
! Subprogram not used 
! Subprogram not used       ! first, look for an exact match (assuming longitude symmetries)
! Subprogram not used       twinlon_off = mod((lon-1)+(clat_p_cnt(lat)/2), clat_p_cnt(lat))
! Subprogram not used       twingcol_latlon = clat_p_idx(lat) + twinlon_off
! Subprogram not used       twingcol = latlon_to_dyn_gcol_map(twingcol_latlon)
! Subprogram not used 
! Subprogram not used       ! otherwise, look around for an approximate match using lonlat
! Subprogram not used       ! column ordering
! Subprogram not used       if ((lon_p(twingcol) .ne. twinlon) .or. &
! Subprogram not used           (lat_p(twingcol) .ne. lat)) then
! Subprogram not used          twingcol_lonlat = clon_p_idx(twinlon)
! Subprogram not used          twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
! Subprogram not used          min_diff = abs(lat_p(twingcol) - lat)
! Subprogram not used          min_i = 0
! Subprogram not used          do i = 1, clon_p_cnt(twinlon)-1
! Subprogram not used             twingcol_lonlat = clon_p_idx(twinlon)+i
! Subprogram not used             twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
! Subprogram not used             diff = abs(lat_p(twingcol) - lat)
! Subprogram not used             if (diff < min_diff) then
! Subprogram not used                min_diff = diff
! Subprogram not used                min_i = i
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          twingcol_lonlat = clon_p_idx(twinlon) + min_i
! Subprogram not used          twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used       call get_gcol_block_d(twingcol,npes,jbtwin,ibtwin)
! Subprogram not used       twinproc = get_block_owner_d(jbtwin(1))
! Subprogram not used       twinsmp  = proc_smp_mapx(twinproc)
! Subprogram not used !
! Subprogram not used       if ((twinsmp .eq. smp) .and. &
! Subprogram not used           (knuhcs(twingcol)%chunkid == -1)) then
! Subprogram not used          found = .true.
! Subprogram not used          twingcol_f = twingcol
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine find_twin
!
!========================================================================

   subroutine assign_chunks(npthreads, nsmpx, proc_smp_mapx, &
                            nsmpthreads, nsmpchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Assign chunks to processes, balancing the number of
!          chunks per thread and minimizing the communication costs
!          in dp_coupling subject to the restraint that columns
!          do not migrate outside of the current SMP node.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plev
   use dyn_grid, only: get_gcol_block_cnt_d, get_gcol_block_d,&
                       get_block_owner_d 
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: npthreads(0:npes-1)
                                         ! number of OpenMP threads per process
   integer, intent(in)  :: nsmpx         ! virtual smp count
   integer, intent(in)  :: proc_smp_mapx(0:npes-1)
                                         ! process/virtual smp map
   integer, intent(in)  :: nsmpthreads(0:nsmpx-1)
                                         ! number of OpenMP threads 
                                         ! per virtual SMP
   integer, intent(in)  :: nsmpchunks(0:nsmpx-1)
                                         ! number of chunks assigned 
                                         ! to a given virtual SMP
!---------------------------Local workspace-----------------------------
   integer :: i, jb, p                   ! loop indices
   integer :: cid                        ! chunk id
   integer :: smp                        ! SMP index
   integer :: curgcol                    ! global column index
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: ntsks_smpx(0:nsmpx-1)      ! number of processes per virtual SMP
   integer :: smp_proc_mapx(0:nsmpx-1,max_nproc_smpx)   
                                         ! virtual smp to process id map
   integer :: cid_offset(0:nsmpx)        ! chunk id virtual smp offset
   integer :: ntmp1_smp(0:nsmpx-1)       ! minimum number of chunks per thread
                                         !  in a virtual SMP
   integer :: ntmp2_smp(0:nsmpx-1)       ! number of extra chunks to be assigned
                                         !  in a virtual SMP
   integer :: ntmp3_smp(0:nsmpx-1)       ! number of processes in a virtual
                                         !  SMP that get more extra chunks
                                         !  than the others
   integer :: ntmp4_smp(0:nsmpx-1)       ! number of extra chunks per process
                                         !  in a virtual SMP
   integer :: ntmp1, ntmp2               ! work variables
!  integer :: npchunks(0:npes-1)         ! number of chunks to be assigned to
!                                        !  a given process
   integer :: cur_npchunks(0:npes-1)     ! current number of chunks assigned 
                                         !  to a given process
   integer :: column_count(0:npes-1)     ! number of columns from current chunk
                                         !  assigned to each process in dynamics
                                         !  decomposition
!-----------------------------------------------------------------------
!
! Count number of processes per virtual SMP and determine virtual SMP
! to process id map
!
   ntsks_smpx(:) = 0
   smp_proc_mapx(:,:) = -1
   do p=0,npes-1
      smp = proc_smp_mapx(p)
      ntsks_smpx(smp) = ntsks_smpx(smp) + 1
      smp_proc_mapx(smp,ntsks_smpx(smp)) = p
   enddo
!
! Determine chunk id ranges for each virtual SMP
!
   cid_offset(0) = 1
   do smp=1,nsmpx
      cid_offset(smp) = cid_offset(smp-1) + nsmpchunks(smp-1)
   enddo
!
! Determine number of chunks to assign to each process
!
   do smp=0,nsmpx-1
!
! Minimum number of chunks per thread
      ntmp1_smp(smp) = nsmpchunks(smp)/nsmpthreads(smp)

! Number of extra chunks to be assigned
      ntmp2_smp(smp) = mod(nsmpchunks(smp),nsmpthreads(smp))

! Number of processes that get more extra chunks than the others
      ntmp3_smp(smp) = mod(ntmp2_smp(smp),ntsks_smpx(smp))

! Number of extra chunks per process
      ntmp4_smp(smp) = ntmp2_smp(smp)/ntsks_smpx(smp)
      if (ntmp3_smp(smp) > 0) then
         ntmp4_smp(smp) = ntmp4_smp(smp) + 1
      endif
   enddo

   do p=0,npes-1
      smp = proc_smp_mapx(p)

! Update number of extra chunks
      if (ntmp2_smp(smp) > ntmp4_smp(smp)) then
         ntmp2_smp(smp) = ntmp2_smp(smp) - ntmp4_smp(smp)
      else
         ntmp4_smp(smp) = ntmp2_smp(smp)
         ntmp2_smp(smp) = 0
         ntmp3_smp(smp) = 0
      endif

! Set number of chunks
      npchunks(p) = ntmp1_smp(smp)*npthreads(p) + ntmp4_smp(smp)

! Update extra chunk increment
      if (ntmp3_smp(smp) > 0) then
         ntmp3_smp(smp) = ntmp3_smp(smp) - 1
         if (ntmp3_smp(smp) .eq. 0) then
            ntmp4_smp(smp) = ntmp4_smp(smp) - 1
         endif
      endif
   enddo

!
! Assign chunks to processes: 
!
   cur_npchunks(:) = 0
!
   do smp=0,nsmpx-1
      do cid=cid_offset(smp),cid_offset(smp+1)-1
!
         do i=1,ntsks_smpx(smp)
            p = smp_proc_mapx(smp,i)
            column_count(p) = 0
         enddo
!
!  For each chunk, determine number of columns in each
!  process within the dynamics.
         do i=1,chunks(cid)%ncols
            curgcol = chunks(cid)%gcol(i)
            block_cnt = get_gcol_block_cnt_d(curgcol)
            call get_gcol_block_d(curgcol,block_cnt,blockids,bcids)
            do jb=1,block_cnt
               p = get_block_owner_d(blockids(jb)) 
               column_count(p) = column_count(p) + 1
            enddo
         enddo
!
!  Eliminate processes that already have their quota of chunks
         do i=1,ntsks_smpx(smp)
            p = smp_proc_mapx(smp,i)
            if (cur_npchunks(p) == npchunks(p)) then
               column_count(p) = -1
            endif
         enddo
!
!  Assign chunk to process with most
!  columns from chunk, from among those still available
         ntmp1 = -1
         ntmp2 = -1
         do i=1,ntsks_smpx(smp)
            p = smp_proc_mapx(smp,i)
            if (column_count(p) > ntmp1) then
               ntmp1 = column_count(p)
               ntmp2 = p
            endif
         enddo
         cur_npchunks(ntmp2) = cur_npchunks(ntmp2) + 1
         chunks(cid)%owner   = ntmp2

!  Update total number of columns assigned to this process
         gs_col_num(ntmp2)   = gs_col_num(ntmp2) + chunks(cid)%ncols
!
      enddo
!
   enddo
!
   return
   end subroutine assign_chunks
!
!========================================================================

!#######################################################################

end module phys_grid
