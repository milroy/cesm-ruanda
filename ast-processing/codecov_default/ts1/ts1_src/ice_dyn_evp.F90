!=======================================================================
!BOP
!
! !MODULE: ice_dyn_evp - elastic-viscous-plastic sea ice dynamics model
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
!
! !REVISION HISTORY:
!  SVN:$Id: ice_dyn_evp.F90 100 2008-01-29 00:25:32Z eclare $
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)
! 
! !INTERFACE:
!
      module ice_dyn_evp
!
! !USES:
!
      use ice_kinds_mod
      use ice_fileunits
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size
      use ice_constants
      use ice_exit, only: abort_ice
      use perf_mod,        only: t_startf, t_stopf, t_barrierf
!
!EOP
!
      implicit none
      save

      ! namelist parameters

      integer (kind=int_kind) :: &
         kdyn     , & ! type of dynamics ( 1 = evp )
         ndte         ! number of subcycles:  ndte=dt/dte
      logical (kind=log_kind) :: &
         maskhalo_dyn , &  ! turn on masked halo updates in subcycling
         maskhalo_stress , &  ! turn on masked halo updates in stress update for tripole
         splitcomm_dyn        ! turn on overlapping of halo update and computation in subcycling

      logical (kind=log_kind) :: &
         evp_damping  ! if true, use evp damping procedure

      ! other EVP parameters

      character (len=char_len) :: & 
         yield_curve  ! 'ellipse' ('teardrop' needs further testing)
                                                                      ! 
      real (kind=dbl_kind), parameter :: &
         dragw = dragio * rhow, &
                         ! drag coefficient for water on ice *rhow (kg/m^3)
         eyc = 0.36_dbl_kind, &
                         ! coefficient for calculating the parameter E
         cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
         sinw = c0   , & ! sin(ocean turning angle)  ! turning angle = 0
         a_min = p001, & ! minimum ice area
         m_min = p01     ! minimum ice mass (kg/m^2)

      real (kind=dbl_kind) :: &
         ecci     , & ! 1/e^2
         dtei     , & ! 1/dte, where dte is subcycling timestep (1/s)
         dte2T    , & ! dte/2T
         denom1   , & ! constants for stress equation
         denom2   , & !
         rcon         ! for damping criterion (kg/s)

      real (kind=dbl_kind), allocatable :: & 
         fcor_blk(:,:,:)   ! Coriolis parameter (1/s)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: evp - elastic-viscous-plastic dynamics driver
!
! !INTERFACE:
!
! Subprogram not used       subroutine evp (dt)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Elastic-viscous-plastic dynamics driver
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_boundary
! Subprogram not used       use ice_blocks
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_state
! Subprogram not used       use ice_flux
! Subprogram not used       use ice_grid
! Subprogram not used       use ice_timers
! Subprogram not used       use ice_mechred, only: ice_strength
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          dt      ! time step
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: & 
! Subprogram not used          ksub           , & ! subcycle step
! Subprogram not used          iblk           , & ! block index
! Subprogram not used          ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
! Subprogram not used          i, j
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(max_blocks) :: & 
! Subprogram not used          icellt   , & ! no. of cells where icetmask = 1
! Subprogram not used          icellu       ! no. of cells where iceumask = 1
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
! Subprogram not used          indxti   , & ! compressed index in i-direction
! Subprogram not used          indxtj   , & ! compressed index in j-direction
! Subprogram not used          indxui   , & ! compressed index in i-direction
! Subprogram not used          indxuj       ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
! Subprogram not used          tmass    , & ! total mass of ice and snow (kg/m^2)
! Subprogram not used          waterx   , & ! for ocean stress calculation, x (m/s)
! Subprogram not used          watery   , & ! for ocean stress calculation, y (m/s)
! Subprogram not used          forcex   , & ! work array: combined atm stress and ocn tilt, x
! Subprogram not used          forcey   , & ! work array: combined atm stress and ocn tilt, y
! Subprogram not used          aiu      , & ! ice fraction on u-grid
! Subprogram not used          umass    , & ! total mass of ice and snow (u grid)
! Subprogram not used          umassdtei    ! mass of U-cell/dte (kg/m^2 s)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)
! Subprogram not used 
! Subprogram not used !      real (kind=dbl_kind), dimension(nx_block,ny_block,8,max_blocks):: &
! Subprogram not used !         str8         ! stress combinations for momentum equation
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), allocatable ::  &
! Subprogram not used          str8(:,:,:,:)         ! stress combinations for momentum equation
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
! Subprogram not used          icetmask   ! ice extent mask (T-cell)
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block           ! block information for current block
! Subprogram not used       
! Subprogram not used       integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
! Subprogram not used          halomask     ! mask for masked halo creation
! Subprogram not used       type (ice_halo) :: &
! Subprogram not used          halo_info_mask          !  ghost cell update info
! Subprogram not used 
! Subprogram not used !      call ice_timer_start(timer_dynamics) ! dynamics
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used        ! This call is needed only if dt changes during runtime.
! Subprogram not used !echmod: automate this
! Subprogram not used !      call set_evp_parameters (dt)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! boundary updates
! Subprogram not used       ! commented out because the ghost cells are freshly 
! Subprogram not used       ! updated after cleanup_itd
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      call ice_timer_start(timer_bound)
! Subprogram not used !      call ice_HaloUpdate (aice,              halo_info, &
! Subprogram not used !                           field_loc_center,  field_type_scalar)
! Subprogram not used !      call ice_HaloUpdate (vice,              halo_info, &
! Subprogram not used !                           field_loc_center,  field_type_scalar)
! Subprogram not used !      call ice_HaloUpdate (vsno,              halo_info, &
! Subprogram not used !                           field_loc_center,  field_type_scalar)
! Subprogram not used !      call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used !     call t_barrierf ('cice_dyn_evp_prep_BARRIER',MPI_COMM_ICE)
! Subprogram not used 
! Subprogram not used       call t_barrierf ('cice_evp_prep1_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf   ('cice_evp_prep1')
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used          do j = 1, ny_block 
! Subprogram not used          do i = 1, nx_block 
! Subprogram not used             rdg_conv (i,j,iblk) = c0 
! Subprogram not used             rdg_shear(i,j,iblk) = c0 
! Subprogram not used             divu (i,j,iblk) = c0 
! Subprogram not used             shear(i,j,iblk) = c0 
! Subprogram not used             prs_sig(i,j,iblk) = c0 
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! preparation for dynamics
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          call evp_prep1 (nx_block,           ny_block,           & 
! Subprogram not used                          ilo, ihi,           jlo, jhi,           &
! Subprogram not used                          aice    (:,:,iblk), vice    (:,:,iblk), & 
! Subprogram not used                          vsno    (:,:,iblk), tmask   (:,:,iblk), & 
! Subprogram not used                          strairxT_accum(:,:,iblk), strairyT_accum(:,:,iblk), & 
! Subprogram not used                          strairx (:,:,iblk), strairy (:,:,iblk), & 
! Subprogram not used                          tmass   (:,:,iblk), icetmask(:,:,iblk))
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf   ('cice_evp_prep1')
! Subprogram not used       call t_barrierf('cice_evp_bound1_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_evp_bound1')
! Subprogram not used 
! Subprogram not used !      call ice_timer_start(timer_bound)
! Subprogram not used       call ice_HaloUpdate (icetmask,          halo_info, &
! Subprogram not used                            field_loc_center,  field_type_scalar)
! Subprogram not used !      call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       call t_stopf ('cice_evp_bound1')
! Subprogram not used       call t_barrierf ('cice_evp_convert_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_evp_convert')
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! convert fields from T to U grid
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call to_ugrid(tmass,umass)
! Subprogram not used       call to_ugrid(aice, aiu)
! Subprogram not used 
! Subprogram not used       call t2ugrid_vector(strairx)
! Subprogram not used       call t2ugrid_vector(strairy)
! Subprogram not used 
! Subprogram not used      call t_stopf ('cice_evp_convert')
! Subprogram not used      call t_barrierf ('cice_evp_prep2_BARRIER',MPI_COMM_ICE)
! Subprogram not used      call t_startf ('cice_evp_prep2')
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! more preparation for dynamics
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          call evp_prep2 (nx_block,             ny_block,             & 
! Subprogram not used                          ilo, ihi,             jlo, jhi,             &
! Subprogram not used                          icellt(iblk),         icellu(iblk),         & 
! Subprogram not used                          indxti      (:,iblk), indxtj      (:,iblk), & 
! Subprogram not used                          indxui      (:,iblk), indxuj      (:,iblk), & 
! Subprogram not used                          aiu       (:,:,iblk), umass     (:,:,iblk), & 
! Subprogram not used                          umassdtei (:,:,iblk), fcor_blk  (:,:,iblk), & 
! Subprogram not used                          umask     (:,:,iblk),                       & 
! Subprogram not used                          uocn      (:,:,iblk), vocn      (:,:,iblk), & 
! Subprogram not used                          strairx   (:,:,iblk), strairy   (:,:,iblk), & 
! Subprogram not used                          ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk), &  
! Subprogram not used                          icetmask  (:,:,iblk), iceumask  (:,:,iblk), & 
! Subprogram not used                          fm        (:,:,iblk),                       & 
! Subprogram not used                          strtltx   (:,:,iblk), strtlty   (:,:,iblk), & 
! Subprogram not used                          strocnx   (:,:,iblk), strocny   (:,:,iblk), & 
! Subprogram not used                          strintx   (:,:,iblk), strinty   (:,:,iblk), & 
! Subprogram not used                          waterx    (:,:,iblk), watery    (:,:,iblk), & 
! Subprogram not used                          forcex    (:,:,iblk), forcey    (:,:,iblk), & 
! Subprogram not used                          stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
! Subprogram not used                          stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
! Subprogram not used                          stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
! Subprogram not used                          stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
! Subprogram not used                          stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
! Subprogram not used                          stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
! Subprogram not used                          uvel      (:,:,iblk), vvel      (:,:,iblk))
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! ice strength
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call ice_strength (nx_block, ny_block,   & 
! Subprogram not used                             ilo, ihi, jlo, jhi,   &
! Subprogram not used                             icellt(iblk),         & 
! Subprogram not used                             indxti      (:,iblk), & 
! Subprogram not used                             indxtj      (:,iblk), & 
! Subprogram not used                             aice    (:,:,  iblk), & 
! Subprogram not used                             vice    (:,:,  iblk), & 
! Subprogram not used                             aice0   (:,:,  iblk), & 
! Subprogram not used                             aicen   (:,:,:,iblk), &  
! Subprogram not used                             vicen   (:,:,:,iblk), & 
! Subprogram not used                             strength(:,:,  iblk) )
! Subprogram not used 
! Subprogram not used       enddo  ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf ('cice_evp_prep2')
! Subprogram not used       call t_barrierf ('cice_evp_bound2_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_evp_bound2')
! Subprogram not used 
! Subprogram not used       allocate(fld2(nx_block,ny_block,2,max_blocks))
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1,nblocks
! Subprogram not used          fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used          fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call ice_HaloUpdate (strength,           halo_info, &
! Subprogram not used                            field_loc_center,   field_type_scalar)
! Subprogram not used       ! velocities may have changed in evp_prep2
! Subprogram not used       call ice_HaloUpdate (fld2,               halo_info, &
! Subprogram not used                            field_loc_NEcorner, field_type_vector)
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1,nblocks
! Subprogram not used          uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
! Subprogram not used          vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used      call t_stopf ('cice_evp_bound2')
! Subprogram not used      if (maskhalo_dyn) then
! Subprogram not used         call t_barrierf ('cice_evp_halocreate_BARRIER',MPI_COMM_ICE)
! Subprogram not used         call t_startf('cice_evp_halocreate')
! Subprogram not used         call ice_HaloMask(halo_info_mask, halo_info, icetmask)
! Subprogram not used         call t_stopf ('cice_evp_halocreate')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      allocate(str8(nx_block,ny_block,8,nblocks))
! Subprogram not used 
! Subprogram not used   if (splitcomm_dyn) then
! Subprogram not used 
! Subprogram not used       do ksub = 1,ndte        ! subcycling
! Subprogram not used 
! Subprogram not used          call t_barrierf ('cice_evp_subcycling_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf ('cice_evp_subcycling')
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! send halo update, skip on first subcycle
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (ksub > 1) then
! Subprogram not used             call t_startf ('cice_evp_bound3i')
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used             do iblk = 1,nblocks
! Subprogram not used                fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used                fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used             call t_stopf ('cice_evp_bound3i')
! Subprogram not used             call t_startf ('cice_evp_bound3s')
! Subprogram not used             if (maskhalo_dyn) then
! Subprogram not used                call ice_HaloUpdate (fld2,               halo_info_mask, &
! Subprogram not used                                     field_loc_NEcorner, field_type_vector, mode='send')
! Subprogram not used             else
! Subprogram not used                call ice_HaloUpdate (fld2,               halo_info, &
! Subprogram not used                                     field_loc_NEcorner, field_type_vector, mode='send')
! Subprogram not used             endif
! Subprogram not used             call t_stopf ('cice_evp_bound3s')
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! stress tensor equation, total surface stress, phase 1 = non edge pts
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             ilo = this_block%ilo
! Subprogram not used             ihi = this_block%ihi
! Subprogram not used             jlo = this_block%jlo
! Subprogram not used             jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used             call t_startf ('cice_evp_stress1')
! Subprogram not used             call stress (1, ilo, ihi, jlo, jhi,                      &
! Subprogram not used                          nx_block,             ny_block,             & 
! Subprogram not used                          ksub,                 icellt(iblk),         & 
! Subprogram not used                          indxti      (:,iblk), indxtj      (:,iblk), & 
! Subprogram not used                          uvel      (:,:,iblk), vvel      (:,:,iblk), &     
! Subprogram not used                          dxt       (:,:,iblk), dyt       (:,:,iblk), & 
! Subprogram not used                          dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
! Subprogram not used                          cxp       (:,:,iblk), cyp       (:,:,iblk), & 
! Subprogram not used                          cxm       (:,:,iblk), cym       (:,:,iblk), & 
! Subprogram not used                          tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
! Subprogram not used                          strength  (:,:,iblk),                       & 
! Subprogram not used                          stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
! Subprogram not used                          stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
! Subprogram not used                          stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
! Subprogram not used                          stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
! Subprogram not used                          stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
! Subprogram not used                          stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
! Subprogram not used                          shear     (:,:,iblk), divu      (:,:,iblk), & 
! Subprogram not used                          prs_sig   (:,:,iblk),                       & 
! Subprogram not used                          rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
! Subprogram not used                          str8    (:,:,:,iblk) )
! Subprogram not used             call t_stopf ('cice_evp_stress1')
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! recv halo update, skip on first subcycle
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (ksub > 1) then
! Subprogram not used             call t_startf ('cice_evp_bound3r')
! Subprogram not used             if (maskhalo_dyn) then
! Subprogram not used                call ice_HaloUpdate (fld2,               halo_info_mask, &
! Subprogram not used                                     field_loc_NEcorner, field_type_vector, mode='recv')
! Subprogram not used             else
! Subprogram not used                call ice_HaloUpdate (fld2,               halo_info, &
! Subprogram not used                                     field_loc_NEcorner, field_type_vector, mode='recv')
! Subprogram not used             endif
! Subprogram not used             call t_stopf ('cice_evp_bound3r')
! Subprogram not used 
! Subprogram not used             call t_startf ('cice_evp_bound3c')
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used             do iblk = 1,nblocks
! Subprogram not used                uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
! Subprogram not used                vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used             call t_stopf ('cice_evp_bound3c')
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! stress tensor equation, total surface stress, phase 2 = edge points
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             ilo = this_block%ilo
! Subprogram not used             ihi = this_block%ihi
! Subprogram not used             jlo = this_block%jlo
! Subprogram not used             jhi = this_block%jhi
! Subprogram not used             call t_startf ('cice_evp_stress2')
! Subprogram not used             call stress (2, ilo, ihi, jlo, jhi,                      &
! Subprogram not used                          nx_block,             ny_block,             & 
! Subprogram not used                          ksub,                 icellt(iblk),         & 
! Subprogram not used                          indxti      (:,iblk), indxtj      (:,iblk), & 
! Subprogram not used                          uvel      (:,:,iblk), vvel      (:,:,iblk), &     
! Subprogram not used                          dxt       (:,:,iblk), dyt       (:,:,iblk), & 
! Subprogram not used                          dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
! Subprogram not used                          cxp       (:,:,iblk), cyp       (:,:,iblk), & 
! Subprogram not used                          cxm       (:,:,iblk), cym       (:,:,iblk), & 
! Subprogram not used                          tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
! Subprogram not used                          strength  (:,:,iblk),                       & 
! Subprogram not used                          stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
! Subprogram not used                          stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
! Subprogram not used                          stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
! Subprogram not used                          stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
! Subprogram not used                          stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
! Subprogram not used                          stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
! Subprogram not used                          shear     (:,:,iblk), divu      (:,:,iblk), & 
! Subprogram not used                          prs_sig   (:,:,iblk),                       & 
! Subprogram not used                          rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
! Subprogram not used                          str8    (:,:,:,iblk) )
! Subprogram not used             call t_stopf ('cice_evp_stress2')
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! momentum equation
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             call t_startf ('cice_evp_stepu')
! Subprogram not used             call stepu (nx_block,            ny_block,           & 
! Subprogram not used                         icellu       (iblk),                     & 
! Subprogram not used                         indxui     (:,iblk), indxuj    (:,iblk), & 
! Subprogram not used                         aiu      (:,:,iblk), str8  (:,:,:,iblk), & 
! Subprogram not used                         uocn     (:,:,iblk), vocn    (:,:,iblk), &     
! Subprogram not used                         waterx   (:,:,iblk), watery  (:,:,iblk), & 
! Subprogram not used                         forcex   (:,:,iblk), forcey  (:,:,iblk), & 
! Subprogram not used                         umassdtei(:,:,iblk), fm      (:,:,iblk), & 
! Subprogram not used                         uarear   (:,:,iblk),                     & 
! Subprogram not used                         strocnx  (:,:,iblk), strocny (:,:,iblk), & 
! Subprogram not used                         strintx  (:,:,iblk), strinty (:,:,iblk), & 
! Subprogram not used                         uvel     (:,:,iblk), vvel    (:,:,iblk))
! Subprogram not used 
! Subprogram not used             call t_stopf ('cice_evp_stepu')
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          call t_stopf ('cice_evp_subcycling')
! Subprogram not used         
! Subprogram not used       enddo                     ! subcycling
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! one final halo update on the velocities
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call t_barrierf ('cice_evp_bound4_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_evp_bound4')
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1,nblocks
! Subprogram not used          fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used          fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used       if (maskhalo_dyn) then
! Subprogram not used          call ice_HaloUpdate (fld2,               halo_info_mask, &
! Subprogram not used                               field_loc_NEcorner, field_type_vector)
! Subprogram not used        else
! Subprogram not used          call ice_HaloUpdate (fld2,               halo_info, &
! Subprogram not used                               field_loc_NEcorner, field_type_vector)
! Subprogram not used       endif
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1,nblocks
! Subprogram not used          uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
! Subprogram not used          vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used       call t_stopf ('cice_evp_bound4')
! Subprogram not used 
! Subprogram not used   else
! Subprogram not used 
! Subprogram not used       do ksub = 1,ndte        ! subcycling
! Subprogram not used 
! Subprogram not used          call t_barrierf ('cice_evp_subcycling_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf ('cice_evp_subcycling')
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! stress tensor equation, total surface stress, all gridcells
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             ilo = this_block%ilo
! Subprogram not used             ihi = this_block%ihi
! Subprogram not used             jlo = this_block%jlo
! Subprogram not used             jhi = this_block%jhi
! Subprogram not used             call t_startf ('cice_evp_stress0')
! Subprogram not used             call stress (0, ilo, ihi, jlo, jhi,                      &
! Subprogram not used                          nx_block,             ny_block,             & 
! Subprogram not used                          ksub,                 icellt(iblk),         & 
! Subprogram not used                          indxti      (:,iblk), indxtj      (:,iblk), & 
! Subprogram not used                          uvel      (:,:,iblk), vvel      (:,:,iblk), &     
! Subprogram not used                          dxt       (:,:,iblk), dyt       (:,:,iblk), & 
! Subprogram not used                          dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
! Subprogram not used                          cxp       (:,:,iblk), cyp       (:,:,iblk), & 
! Subprogram not used                          cxm       (:,:,iblk), cym       (:,:,iblk), & 
! Subprogram not used                          tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
! Subprogram not used                          strength  (:,:,iblk),                       & 
! Subprogram not used                          stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
! Subprogram not used                          stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
! Subprogram not used                          stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
! Subprogram not used                          stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
! Subprogram not used                          stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
! Subprogram not used                          stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
! Subprogram not used                          shear     (:,:,iblk), divu      (:,:,iblk), & 
! Subprogram not used                          prs_sig   (:,:,iblk),                       & 
! Subprogram not used                          rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
! Subprogram not used                          str8    (:,:,:,iblk) )
! Subprogram not used             call t_stopf ('cice_evp_stress0')
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! momentum equation
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             call t_startf ('cice_evp_stepu')
! Subprogram not used             call stepu (nx_block,            ny_block,           & 
! Subprogram not used                         icellu       (iblk),                     & 
! Subprogram not used                         indxui     (:,iblk), indxuj    (:,iblk), & 
! Subprogram not used                         aiu      (:,:,iblk), str8  (:,:,:,iblk), & 
! Subprogram not used                         uocn     (:,:,iblk), vocn    (:,:,iblk), &     
! Subprogram not used                         waterx   (:,:,iblk), watery  (:,:,iblk), & 
! Subprogram not used                         forcex   (:,:,iblk), forcey  (:,:,iblk), & 
! Subprogram not used                         umassdtei(:,:,iblk), fm      (:,:,iblk), & 
! Subprogram not used                         uarear   (:,:,iblk),                     & 
! Subprogram not used                         strocnx  (:,:,iblk), strocny (:,:,iblk), & 
! Subprogram not used                         strintx  (:,:,iblk), strinty (:,:,iblk), & 
! Subprogram not used                         uvel     (:,:,iblk), vvel    (:,:,iblk))
! Subprogram not used 
! Subprogram not used             call t_stopf ('cice_evp_stepu')
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used          ! halo update
! Subprogram not used          !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             call t_startf ('cice_evp_bound3i')
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used             do iblk = 1,nblocks
! Subprogram not used                fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used                fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used             call t_stopf ('cice_evp_bound3i')
! Subprogram not used             call t_startf ('cice_evp_bound3sr')
! Subprogram not used             if (maskhalo_dyn) then
! Subprogram not used                call ice_HaloUpdate (fld2,               halo_info_mask, &
! Subprogram not used                                     field_loc_NEcorner, field_type_vector)
! Subprogram not used             else
! Subprogram not used                call ice_HaloUpdate (fld2,               halo_info, &
! Subprogram not used                                     field_loc_NEcorner, field_type_vector)
! Subprogram not used             endif
! Subprogram not used             call t_stopf ('cice_evp_bound3sr')
! Subprogram not used 
! Subprogram not used             call t_startf ('cice_evp_bound3c')
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used             do iblk = 1,nblocks
! Subprogram not used                uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
! Subprogram not used                vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used             call t_stopf ('cice_evp_bound3c')
! Subprogram not used 
! Subprogram not used          call t_stopf ('cice_evp_subcycling')
! Subprogram not used         
! Subprogram not used       enddo                     ! subcycling
! Subprogram not used 
! Subprogram not used   endif   ! splitcomm subcycling
! Subprogram not used 
! Subprogram not used       deallocate(str8)
! Subprogram not used       deallocate(fld2)
! Subprogram not used 
! Subprogram not used       if (maskhalo_dyn) then
! Subprogram not used          call t_barrierf ('cice_evp_halodestroy_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf ('cice_evp_halodestroy')
! Subprogram not used          call ice_HaloDestroy(halo_info_mask)
! Subprogram not used          call t_stopf ('cice_evp_halodestroy')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call t_barrierf ('cice_evp_tpupd_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_evp_tpupd')
! Subprogram not used 
! Subprogram not used       ! Force symmetry across the tripole seam
! Subprogram not used       if (trim(grid_type) == 'tripole') then
! Subprogram not used       if (maskhalo_stress) then
! Subprogram not used          !-------------------------------------------------------
! Subprogram not used          ! set halomask to zero because ice_HaloMask always keeps
! Subprogram not used          ! local copies AND tripole zipper communication
! Subprogram not used          !-------------------------------------------------------
! Subprogram not used          call t_barrierf ('cice_evp_tpupdhc_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf('cice_evp_tpupdhc')
! Subprogram not used          halomask = 0
! Subprogram not used          call ice_HaloMask(halo_info_mask, halo_info, halomask)
! Subprogram not used          call t_stopf ('cice_evp_tpupdhc')
! Subprogram not used 
! Subprogram not used          call t_barrierf ('cice_evp_tpupdhu_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf('cice_evp_tpupdhu')
! Subprogram not used          call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used 
! Subprogram not used          call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used 
! Subprogram not used          call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info_mask, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call t_stopf('cice_evp_tpupdhu')
! Subprogram not used 
! Subprogram not used          call t_barrierf ('cice_evp_tpupdhd_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf ('cice_evp_tpupdhd')
! Subprogram not used          call ice_HaloDestroy(halo_info_mask)
! Subprogram not used          call t_stopf ('cice_evp_tpupdhd')
! Subprogram not used       else
! Subprogram not used 
! Subprogram not used          call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used 
! Subprogram not used          call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used 
! Subprogram not used          call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used          call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
! Subprogram not used                               field_loc_center,  field_type_scalar)
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call t_stopf ('cice_evp_tpupd')
! Subprogram not used       call t_barrierf ('cice_evp_finish_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_evp_finish')
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! ice-ocean stress
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used          call evp_finish                               & 
! Subprogram not used               (nx_block,           ny_block,           & 
! Subprogram not used                icellu      (iblk),                     & 
! Subprogram not used                indxui    (:,iblk), indxuj    (:,iblk), & 
! Subprogram not used                uvel    (:,:,iblk), vvel    (:,:,iblk), & 
! Subprogram not used                uocn    (:,:,iblk), vocn    (:,:,iblk), & 
! Subprogram not used                aiu     (:,:,iblk),                     &
! Subprogram not used                strocnx (:,:,iblk), strocny (:,:,iblk), & 
! Subprogram not used                strocnxT(:,:,iblk), strocnyT(:,:,iblk))
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call u2tgrid_vector(strocnxT)    ! shift
! Subprogram not used       call u2tgrid_vector(strocnyT)
! Subprogram not used 
! Subprogram not used !      call ice_timer_stop(timer_dynamics)    ! dynamics
! Subprogram not used       call t_stopf ('cice_evp_finish')
! Subprogram not used 
! Subprogram not used       end subroutine evp

!=======================================================================
!BOP
!
! !IROUTINE: init_evp - initialize parameters needed for evp dynamics
!
! !INTERFACE:
!
      subroutine init_evp (dt)
!
! !DESCRIPTION:
!
! Initialize parameters and variables needed for the evp dynamics
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_flux
      use ice_grid
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, &
         iblk            ! block index

      real (kind=dbl_kind) :: &
         dte         , & ! subcycling timestep for EVP dynamics, s
         ecc         , & ! (ratio of major to minor ellipse axes)^2
         tdamp2          ! 2(wave damping time scale T)

      call set_evp_parameters (dt)

      if (my_task == master_task) then
         write(nu_diag,*) 'dt_dyn  = ',dt
         write(nu_diag,*) 'dte     = ',dt/real(ndte,kind=dbl_kind)
         write(nu_diag,*) 'tdamp   =', eyc*dt
      endif

      allocate(fcor_blk(nx_block,ny_block,max_blocks))

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block

         ! velocity
         uvel(i,j,iblk) = c0    ! m/s
         vvel(i,j,iblk) = c0    ! m/s

         ! strain rates
         divu (i,j,iblk) = c0
         shear(i,j,iblk) = c0
         rdg_conv (i,j,iblk) = c0
         rdg_shear(i,j,iblk) = c0

         ! Coriolis parameter
!!         fcor_blk(i,j,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
         fcor_blk(i,j,iblk) = c2*omega*sin(ULAT(i,j,iblk)) ! 1/s

         ! stress tensor,  kg/s^2
         stressp_1 (i,j,iblk) = c0
         stressp_2 (i,j,iblk) = c0
         stressp_3 (i,j,iblk) = c0
         stressp_4 (i,j,iblk) = c0
         stressm_1 (i,j,iblk) = c0
         stressm_2 (i,j,iblk) = c0
         stressm_3 (i,j,iblk) = c0
         stressm_4 (i,j,iblk) = c0
         stress12_1(i,j,iblk) = c0
         stress12_2(i,j,iblk) = c0
         stress12_3(i,j,iblk) = c0
         stress12_4(i,j,iblk) = c0

         ! ice extent mask on velocity points
         iceumask(i,j,iblk) = .false.

      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_evp

!=======================================================================
!BOP
!
! !IROUTINE: set_evp_parameters - set parameters for evp dynamics
!
! !INTERFACE:
!
      subroutine set_evp_parameters (dt)
!
! !DESCRIPTION:
!
! Set parameters needed for the evp dynamics.
! Note: This subroutine is currently called only during initialization.
!       If the dynamics time step can vary during runtime, it should
!        be called whenever the time step changes.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      real (kind=dbl_kind) :: &
         dte         , & ! subcycling timestep for EVP dynamics, s
         ecc         , & ! (ratio of major to minor ellipse axes)^2
         tdamp2          ! 2*(wave damping time scale T)

      ! elastic time step
      dte = dt/real(ndte,kind=dbl_kind)        ! s
      dtei = c1/dte              ! 1/s

      ! major/minor axis length ratio, squared
      ecc  = c4
      ecci = p25                  ! 1/ecc

      ! constants for stress equation
      tdamp2 = c2*eyc*dt                    ! s
      dte2T = dte/tdamp2                    ! ellipse (unitless)
      denom1 = c1/(c1+dte2T)
      denom2 = c1/(c1+dte2T*ecc)
      rcon = 1230._dbl_kind*eyc*dt*dtei**2  ! kg/s

      end subroutine set_evp_parameters

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep1 - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
! Subprogram not used       subroutine evp_prep1 (nx_block,  ny_block, & 
! Subprogram not used                             ilo, ihi,  jlo, jhi, &
! Subprogram not used                             aice,      vice,     & 
! Subprogram not used                             vsno,      tmask,    & 
! Subprogram not used                             strairxT,  strairyT, & 
! Subprogram not used                             strairx,   strairy,  & 
! Subprogram not used                             tmass,     icetmask)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Computes quantities needed in the stress tensor (sigma)
! Subprogram not used ! and momentum (u) equations, but which do not change during
! Subprogram not used ! the thermodynamics/transport time step:
! Subprogram not used ! ice mass and ice extent masks
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi       ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(in) :: &
! Subprogram not used          aice    , & ! concentration of ice
! Subprogram not used          vice    , & ! volume per unit area of ice          (m)
! Subprogram not used          vsno    , & ! volume per unit area of snow         (m)
! Subprogram not used          strairxT, & ! stress on ice by air, x-direction
! Subprogram not used          strairyT    ! stress on ice by air, y-direction
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(in) :: &
! Subprogram not used          tmask       ! land/boundary mask, thickness (T-cell)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(out) :: &
! Subprogram not used          strairx , & ! stress on ice by air, x-direction
! Subprogram not used          strairy , & ! stress on ice by air, y-direction
! Subprogram not used          tmass       ! total mass of ice and snow (kg/m^2)
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(out) :: &
! Subprogram not used          icetmask    ! ice extent mask (T-cell)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension(nx_block,ny_block) :: &
! Subprogram not used          tmphm               ! temporary mask
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! total mass of ice and snow, centered in T-cell
! Subprogram not used       ! NOTE: vice and vsno must be up to date in all grid cells,
! Subprogram not used       !       including ghost cells
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          if (tmask(i,j)) then
! Subprogram not used             tmass(i,j) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
! Subprogram not used          else
! Subprogram not used             tmass(i,j) = c0
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! ice extent mask (T-cells)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          tmphm(i,j) = tmask(i,j) .and. (aice (i,j) > a_min) & 
! Subprogram not used                                  .and. (tmass(i,j) > m_min)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! prep to convert to U grid
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          ! these quantities include the factor of aice needed for
! Subprogram not used          ! correct treatment of free drift
! Subprogram not used          strairx(i,j) = strairxT(i,j)
! Subprogram not used          strairy(i,j) = strairyT(i,j)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! augmented mask (land + open ocean)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          icetmask (i,j) = 0
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do j = jlo, jhi
! Subprogram not used       do i = ilo, ihi
! Subprogram not used 
! Subprogram not used          ! extend ice extent mask (T-cells) to points around pack
! Subprogram not used          if (tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or. & 
! Subprogram not used              tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or. & 
! Subprogram not used              tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1) ) then
! Subprogram not used             icetmask(i,j) = 1
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (.not.tmask(i,j)) icetmask(i,j) = 0
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine evp_prep1

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep2 - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
! Subprogram not used       subroutine evp_prep2 (nx_block,   ny_block,   & 
! Subprogram not used                             ilo, ihi,  jlo, jhi,    &
! Subprogram not used                             icellt,     icellu,     & 
! Subprogram not used                             indxti,     indxtj,     & 
! Subprogram not used                             indxui,     indxuj,     & 
! Subprogram not used                             aiu,        umass,      & 
! Subprogram not used                             umassdtei,  fcor,       & 
! Subprogram not used                             umask,                  & 
! Subprogram not used                             uocn,       vocn,       & 
! Subprogram not used                             strairx,    strairy,    & 
! Subprogram not used                             ss_tltx,    ss_tlty,    &  
! Subprogram not used                             icetmask,   iceumask,   & 
! Subprogram not used                             fm,                     & 
! Subprogram not used                             strtltx,    strtlty,    & 
! Subprogram not used                             strocnx,    strocny,    &
! Subprogram not used                             strintx,    strinty,    &
! Subprogram not used                             waterx,     watery,     & 
! Subprogram not used                             forcex,     forcey,     &     
! Subprogram not used                             stressp_1,  stressp_2,  &   
! Subprogram not used                             stressp_3,  stressp_4,  & 
! Subprogram not used                             stressm_1,  stressm_2,  & 
! Subprogram not used                             stressm_3,  stressm_4,  & 
! Subprogram not used                             stress12_1, stress12_2, & 
! Subprogram not used                             stress12_3, stress12_4, & 
! Subprogram not used                             uvel,       vvel)
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Computes quantities needed in the stress tensor (sigma)
! Subprogram not used ! and momentum (u) equations, but which do not change during
! Subprogram not used ! the thermodynamics/transport time step:
! Subprogram not used ! --wind stress shift to U grid,
! Subprogram not used ! --ice mass and ice extent masks,
! Subprogram not used ! initializes ice velocity for new points to ocean sfc current
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi       ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(out) :: &
! Subprogram not used          icellt   , & ! no. of cells where icetmask = 1
! Subprogram not used          icellu       ! no. of cells where iceumask = 1
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), & 
! Subprogram not used          intent(out) :: &
! Subprogram not used          indxti   , & ! compressed index in i-direction
! Subprogram not used          indxtj   , & ! compressed index in j-direction
! Subprogram not used          indxui   , & ! compressed index in i-direction
! Subprogram not used          indxuj       ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(in) :: &
! Subprogram not used          umask       ! land/boundary mask, thickness (U-cell)
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(in) :: &
! Subprogram not used          icetmask    ! ice extent mask (T-cell)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(inout) :: &
! Subprogram not used          iceumask    ! ice extent mask (U-cell)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
! Subprogram not used          aiu     , & ! ice fraction on u-grid
! Subprogram not used          umass   , & ! total mass of ice and snow (u grid)
! Subprogram not used          fcor    , & ! Coriolis parameter (1/s)
! Subprogram not used          strairx , & ! stress on ice by air, x-direction
! Subprogram not used          strairy , & ! stress on ice by air, y-direction
! Subprogram not used          uocn    , & ! ocean current, x-direction (m/s)
! Subprogram not used          vocn    , & ! ocean current, y-direction (m/s)
! Subprogram not used          ss_tltx , & ! sea surface slope, x-direction (m/m)
! Subprogram not used          ss_tlty     ! sea surface slope, y-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(out) :: &
! Subprogram not used          umassdtei,& ! mass of U-cell/dte (kg/m^2 s)
! Subprogram not used          waterx  , & ! for ocean stress calculation, x (m/s)
! Subprogram not used          watery  , & ! for ocean stress calculation, y (m/s)
! Subprogram not used          forcex  , & ! work array: combined atm stress and ocn tilt, x
! Subprogram not used          forcey      ! work array: combined atm stress and ocn tilt, y
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(inout) :: &
! Subprogram not used          fm      , & ! Coriolis param. * mass in U-cell (kg/s)
! Subprogram not used          stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
! Subprogram not used          stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
! Subprogram not used          stress12_1,stress12_2,stress12_3,stress12_4, & ! sigma12
! Subprogram not used          uvel    , & ! x-component of velocity (m/s)
! Subprogram not used          vvel    , & ! y-component of velocity (m/s)
! Subprogram not used          strtltx , & ! stress due to sea surface slope, x-direction
! Subprogram not used          strtlty , & ! stress due to sea surface slope, y-direction
! Subprogram not used          strocnx , & ! ice-ocean stress, x-direction
! Subprogram not used          strocny , & ! ice-ocean stress, y-direction
! Subprogram not used          strintx , & ! divergence of internal ice stress, x (N/m^2)
! Subprogram not used          strinty     ! divergence of internal ice stress, y (N/m^2)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, ij
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension(nx_block,ny_block) :: &
! Subprogram not used          iceumask_old      ! old-time iceumask
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          waterx   (i,j) = c0
! Subprogram not used          watery   (i,j) = c0
! Subprogram not used          forcex   (i,j) = c0
! Subprogram not used          forcey   (i,j) = c0
! Subprogram not used          umassdtei(i,j) = c0
! Subprogram not used 
! Subprogram not used          if (icetmask(i,j)==0) then
! Subprogram not used             stressp_1 (i,j) = c0
! Subprogram not used             stressp_2 (i,j) = c0
! Subprogram not used             stressp_3 (i,j) = c0
! Subprogram not used             stressp_4 (i,j) = c0
! Subprogram not used             stressm_1 (i,j) = c0
! Subprogram not used             stressm_2 (i,j) = c0
! Subprogram not used             stressm_3 (i,j) = c0
! Subprogram not used             stressm_4 (i,j) = c0
! Subprogram not used             stress12_1(i,j) = c0
! Subprogram not used             stress12_2(i,j) = c0
! Subprogram not used             stress12_3(i,j) = c0
! Subprogram not used             stress12_4(i,j) = c0
! Subprogram not used          endif                  ! icetmask
! Subprogram not used       enddo                     ! i
! Subprogram not used       enddo                     ! j
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Identify cells where icetmask = 1
! Subprogram not used       ! Note: The icellt mask includes north and east ghost cells
! Subprogram not used       !       where stresses are needed.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       icellt = 0
! Subprogram not used       do j = jlo, jhi+1
! Subprogram not used       do i = ilo, ihi+1
! Subprogram not used          if (icetmask(i,j) == 1) then
! Subprogram not used             icellt = icellt + 1
! Subprogram not used             indxti(icellt) = i
! Subprogram not used             indxtj(icellt) = j
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Define iceumask
! Subprogram not used       ! Identify cells where iceumask is true
! Subprogram not used       ! Initialize velocity where needed
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       icellu = 0
! Subprogram not used       do j = jlo, jhi
! Subprogram not used       do i = ilo, ihi
! Subprogram not used 
! Subprogram not used          ! ice extent mask (U-cells)
! Subprogram not used          iceumask_old(i,j) = iceumask(i,j) ! save
! Subprogram not used          iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j) > a_min) & 
! Subprogram not used                                       .and. (umass(i,j) > m_min)
! Subprogram not used 
! Subprogram not used          if (iceumask(i,j)) then
! Subprogram not used             icellu = icellu + 1
! Subprogram not used             indxui(icellu) = i
! Subprogram not used             indxuj(icellu) = j
! Subprogram not used 
! Subprogram not used             ! initialize velocity for new ice points to ocean sfc current
! Subprogram not used             if (.not. iceumask_old(i,j)) then
! Subprogram not used                uvel(i,j) = uocn(i,j)
! Subprogram not used                vvel(i,j) = vocn(i,j)
! Subprogram not used             endif
! Subprogram not used          else
! Subprogram not used             ! set velocity and stresses to zero for masked-out points
! Subprogram not used             uvel(i,j)    = c0
! Subprogram not used             vvel(i,j)    = c0
! Subprogram not used             strintx(i,j) = c0
! Subprogram not used             strinty(i,j) = c0
! Subprogram not used             strocnx(i,j) = c0
! Subprogram not used             strocny(i,j) = c0
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Define variables for momentum equation
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do ij = 1, icellu
! Subprogram not used          i = indxui(ij)
! Subprogram not used          j = indxuj(ij)
! Subprogram not used 
! Subprogram not used          umassdtei(i,j) = umass(i,j)*dtei ! m/dte, kg/m^2 s
! Subprogram not used          fm(i,j) = fcor(i,j)*umass(i,j)   ! Coriolis * mass
! Subprogram not used 
! Subprogram not used          ! for ocean stress
! Subprogram not used          waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw
! Subprogram not used          watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw
! Subprogram not used 
! Subprogram not used          ! combine tilt with wind stress
! Subprogram not used          strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
! Subprogram not used          strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
! Subprogram not used          forcex(i,j) = strairx(i,j) + strtltx(i,j)
! Subprogram not used          forcey(i,j) = strairy(i,j) + strtlty(i,j)
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine evp_prep2

!=======================================================================
!BOP
!
! !IROUTINE: stress - computes strain rates and internal stress components
!
! !INTERFACE:
!
! Subprogram not used       subroutine stress (phase, ilo, ihi, jlo, jhi, &
! Subprogram not used                          nx_block,   ny_block,   & 
! Subprogram not used                          ksub,       icellt,     & 
! Subprogram not used                          indxti,     indxtj,     & 
! Subprogram not used                          uvel,       vvel,       & 
! Subprogram not used                          dxt,        dyt,        & 
! Subprogram not used                          dxhy,       dyhx,       & 
! Subprogram not used                          cxp,        cyp,        & 
! Subprogram not used                          cxm,        cym,        & 
! Subprogram not used                          tarear,     tinyarea,   & 
! Subprogram not used                          strength,               & 
! Subprogram not used                          stressp_1,  stressp_2,  & 
! Subprogram not used                          stressp_3,  stressp_4,  & 
! Subprogram not used                          stressm_1,  stressm_2,  & 
! Subprogram not used                          stressm_3,  stressm_4,  & 
! Subprogram not used                          stress12_1, stress12_2, & 
! Subprogram not used                          stress12_3, stress12_4, & 
! Subprogram not used                          shear,      divu,       & 
! Subprogram not used                          prs_sig,                & 
! Subprogram not used                          rdg_conv,   rdg_shear,  & 
! Subprogram not used                          str8 )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Computes the rates of strain and internal stress components for
! Subprogram not used ! each of the four corners on each T-grid cell.
! Subprogram not used ! Computes stress terms for the momentum equation
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: & 
! Subprogram not used          phase,              & ! phase
! Subprogram not used          ilo, ihi, jlo, jhi, & ! block dimensions
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          ksub              , & ! subcycling step
! Subprogram not used          icellt                ! no. of cells where icetmask = 1
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), & 
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxti   , & ! compressed index in i-direction
! Subprogram not used          indxtj       ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
! Subprogram not used          strength , & ! ice strength (N/m)
! Subprogram not used          uvel     , & ! x-component of velocity (m/s)
! Subprogram not used          vvel     , & ! y-component of velocity (m/s)
! Subprogram not used          dxt      , & ! width of T-cell through the middle (m)
! Subprogram not used          dyt      , & ! height of T-cell through the middle (m)
! Subprogram not used          dxhy     , & ! 0.5*(HTE - HTE)
! Subprogram not used          dyhx     , & ! 0.5*(HTN - HTN)
! Subprogram not used          cyp      , & ! 1.5*HTE - 0.5*HTE
! Subprogram not used          cxp      , & ! 1.5*HTN - 0.5*HTN
! Subprogram not used          cym      , & ! 0.5*HTE - 1.5*HTE
! Subprogram not used          cxm      , & ! 0.5*HTN - 1.5*HTN
! Subprogram not used          tarear   , & ! 1/tarea
! Subprogram not used          tinyarea     ! puny*tarea
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(inout) :: &
! Subprogram not used          stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
! Subprogram not used          stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
! Subprogram not used          stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), & 
! Subprogram not used          intent(inout) :: &
! Subprogram not used          prs_sig  , & ! replacement pressure, for stress calc
! Subprogram not used          shear    , & ! strain rate II component (1/s)
! Subprogram not used          divu     , & ! strain rate I component, velocity divergence (1/s)
! Subprogram not used          rdg_conv , & ! convergence term for ridging (1/s)
! Subprogram not used          rdg_shear    ! shear term for ridging (1/s)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
! Subprogram not used          intent(out) :: &
! Subprogram not used          str8         ! stress combinations
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, ij, chunk
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used         divune, divunw, divuse, divusw            , & ! divergence
! Subprogram not used         tensionne, tensionnw, tensionse, tensionsw, & ! tension
! Subprogram not used         shearne, shearnw, shearse, shearsw        , & ! shearing
! Subprogram not used         Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
! Subprogram not used         c0ne, c0nw, c0se, c0sw                    , & ! useful combinations
! Subprogram not used         c1ne, c1nw, c1se, c1sw                    , &
! Subprogram not used         ssigpn, ssigps, ssigpe, ssigpw            , &
! Subprogram not used         ssigmn, ssigms, ssigme, ssigmw            , &
! Subprogram not used         ssig12n, ssig12s, ssig12e, ssig12w        , &
! Subprogram not used         ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
! Subprogram not used         csigpne, csigpnw, csigpse, csigpsw        , &
! Subprogram not used         csigmne, csigmnw, csigmse, csigmsw        , &
! Subprogram not used         csig12ne, csig12nw, csig12se, csig12sw    , &
! Subprogram not used         str12ew, str12we, str12ns, str12sn        , &
! Subprogram not used         strp_tmp, strm_tmp, tmp
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       if (phase == 0 .or. phase == 1 .or. phase == 2) then
! Subprogram not used          ! valid, 0=all cells, 1=interior cells, 2=edge cells
! Subprogram not used       else
! Subprogram not used          call abort_ice ('ice_dyn stress: illegal phase')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (phase == 1) then
! Subprogram not used          str8(:,:,:) = c0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used       do ij = 1, icellt
! Subprogram not used          i = indxti(ij)
! Subprogram not used          j = indxtj(ij)
! Subprogram not used 
! Subprogram not used !     if ((phase == 1 .and. (i >  ilo .and. i <  ihi) .and. (j >  jlo .and. j <  jhi)) .or. &
! Subprogram not used !         (phase == 2 .and. (i <= ilo .or.  i >= ihi) .and. (j <= jlo .or.  j >= jhi))) then
! Subprogram not used 
! Subprogram not used      if (i >  ilo .and. i <  ihi .and. j > jlo .and. j < jhi) then
! Subprogram not used         chunk = 1
! Subprogram not used      else
! Subprogram not used         chunk = 2
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if (phase == 0 .or. phase == chunk) then
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! strain rates
! Subprogram not used       ! NOTE these are actually strain rates * area  (m^2/s)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          ! divergence  =  e_11 + e_22
! Subprogram not used          divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
! Subprogram not used                    + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
! Subprogram not used          divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
! Subprogram not used                    + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
! Subprogram not used          divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
! Subprogram not used                    + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
! Subprogram not used          divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
! Subprogram not used                    + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )
! Subprogram not used 
! Subprogram not used          ! tension strain rate  =  e_11 - e_22
! Subprogram not used          tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
! Subprogram not used                    +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
! Subprogram not used          tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
! Subprogram not used                    +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
! Subprogram not used          tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
! Subprogram not used                    +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
! Subprogram not used          tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
! Subprogram not used                    +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )
! Subprogram not used 
! Subprogram not used          ! shearing strain rate  =  e_12
! Subprogram not used          shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
! Subprogram not used                  -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
! Subprogram not used          shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
! Subprogram not used                  -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
! Subprogram not used          shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
! Subprogram not used                  -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
! Subprogram not used          shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
! Subprogram not used                  -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )
! Subprogram not used 
! Subprogram not used          ! Delta (in the denominator of zeta, eta)
! Subprogram not used          Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
! Subprogram not used          Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
! Subprogram not used          Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
! Subprogram not used          Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! on last subcycle, save quantities for mechanical redistribution
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          if (ksub == ndte) then
! Subprogram not used             divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
! Subprogram not used             tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
! Subprogram not used             rdg_conv(i,j)  = -min(divu(i,j),c0)
! Subprogram not used             rdg_shear(i,j) = p5*(tmp-abs(divu(i,j))) 
! Subprogram not used 
! Subprogram not used             ! diagnostic only
! Subprogram not used             ! shear = sqrt(tension**2 + shearing**2)
! Subprogram not used             shear(i,j) = p25*tarear(i,j)*sqrt( &
! Subprogram not used                  (tensionne + tensionnw + tensionse + tensionsw)**2 &
! Subprogram not used                 +  (shearne +   shearnw +   shearse +   shearsw)**2)
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! replacement pressure/Delta                   ! kg/s
! Subprogram not used       ! save replacement pressure for principal stress calculation
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          if (evp_damping) then
! Subprogram not used             ! enforce damping criterion
! Subprogram not used             c0ne = min(strength(i,j)/max(Deltane,c4*tinyarea(i,j)),rcon)
! Subprogram not used             c0nw = min(strength(i,j)/max(Deltanw,c4*tinyarea(i,j)),rcon)
! Subprogram not used             c0sw = min(strength(i,j)/max(Deltasw,c4*tinyarea(i,j)),rcon)
! Subprogram not used             c0se = min(strength(i,j)/max(Deltase,c4*tinyarea(i,j)),rcon)
! Subprogram not used             prs_sig(i,j) = strength(i,j)* &
! Subprogram not used                            Deltane/max(Deltane,c4*tinyarea(i,j)) ! ne
! Subprogram not used          else
! Subprogram not used             ! original version
! Subprogram not used             c0ne = strength(i,j)/max(Deltane,tinyarea(i,j))
! Subprogram not used             c0nw = strength(i,j)/max(Deltanw,tinyarea(i,j))
! Subprogram not used             c0sw = strength(i,j)/max(Deltasw,tinyarea(i,j))
! Subprogram not used             c0se = strength(i,j)/max(Deltase,tinyarea(i,j))
! Subprogram not used             prs_sig(i,j) = c0ne*Deltane ! northeast
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          c1ne = c0ne*dte2T
! Subprogram not used          c1nw = c0nw*dte2T
! Subprogram not used          c1sw = c0sw*dte2T
! Subprogram not used          c1se = c0se*dte2T
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! the stresses                            ! kg/s^2
! Subprogram not used       ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          stressp_1(i,j) = (stressp_1(i,j) + c1ne*(divune - Deltane)) &
! Subprogram not used                           * denom1
! Subprogram not used          stressp_2(i,j) = (stressp_2(i,j) + c1nw*(divunw - Deltanw)) &
! Subprogram not used                           * denom1
! Subprogram not used          stressp_3(i,j) = (stressp_3(i,j) + c1sw*(divusw - Deltasw)) &
! Subprogram not used                           * denom1
! Subprogram not used          stressp_4(i,j) = (stressp_4(i,j) + c1se*(divuse - Deltase)) &
! Subprogram not used                           * denom1
! Subprogram not used 
! Subprogram not used          stressm_1(i,j) = (stressm_1(i,j) + c1ne*tensionne) * denom2
! Subprogram not used          stressm_2(i,j) = (stressm_2(i,j) + c1nw*tensionnw) * denom2
! Subprogram not used          stressm_3(i,j) = (stressm_3(i,j) + c1sw*tensionsw) * denom2
! Subprogram not used          stressm_4(i,j) = (stressm_4(i,j) + c1se*tensionse) * denom2
! Subprogram not used 
! Subprogram not used          stress12_1(i,j) = (stress12_1(i,j) + c1ne*shearne*p5) * denom2
! Subprogram not used          stress12_2(i,j) = (stress12_2(i,j) + c1nw*shearnw*p5) * denom2
! Subprogram not used          stress12_3(i,j) = (stress12_3(i,j) + c1sw*shearsw*p5) * denom2
! Subprogram not used          stress12_4(i,j) = (stress12_4(i,j) + c1se*shearse*p5) * denom2
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Eliminate underflows.
! Subprogram not used       ! The following code is commented out because it is relatively 
! Subprogram not used       ! expensive and most compilers include a flag that accomplishes
! Subprogram not used       ! the same thing more efficiently.  This code is cheaper than
! Subprogram not used       ! handling underflows if the compiler lacks a flag; uncomment
! Subprogram not used       ! it in that case.  The compiler flag is often described with the 
! Subprogram not used       ! phrase "flush to zero".
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
! Subprogram not used !      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
! Subprogram not used !      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
! Subprogram not used !      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))
! Subprogram not used 
! Subprogram not used !      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
! Subprogram not used !      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
! Subprogram not used !      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
! Subprogram not used !      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))
! Subprogram not used 
! Subprogram not used !      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
! Subprogram not used !      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
! Subprogram not used !      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
! Subprogram not used !      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! combinations of the stresses for the momentum equation ! kg/s^2
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          ssigpn  = stressp_1(i,j) + stressp_2(i,j)
! Subprogram not used          ssigps  = stressp_3(i,j) + stressp_4(i,j)
! Subprogram not used          ssigpe  = stressp_1(i,j) + stressp_4(i,j)
! Subprogram not used          ssigpw  = stressp_2(i,j) + stressp_3(i,j)
! Subprogram not used          ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
! Subprogram not used          ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055
! Subprogram not used 
! Subprogram not used          ssigmn  = stressm_1(i,j) + stressm_2(i,j)
! Subprogram not used          ssigms  = stressm_3(i,j) + stressm_4(i,j)
! Subprogram not used          ssigme  = stressm_1(i,j) + stressm_4(i,j)
! Subprogram not used          ssigmw  = stressm_2(i,j) + stressm_3(i,j)
! Subprogram not used          ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
! Subprogram not used          ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055
! Subprogram not used 
! Subprogram not used          ssig12n = stress12_1(i,j) + stress12_2(i,j)
! Subprogram not used          ssig12s = stress12_3(i,j) + stress12_4(i,j)
! Subprogram not used          ssig12e = stress12_1(i,j) + stress12_4(i,j)
! Subprogram not used          ssig12w = stress12_2(i,j) + stress12_3(i,j)
! Subprogram not used          ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
! Subprogram not used          ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111
! Subprogram not used 
! Subprogram not used          csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
! Subprogram not used          csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
! Subprogram not used          csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
! Subprogram not used          csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
! Subprogram not used          
! Subprogram not used          csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
! Subprogram not used          csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
! Subprogram not used          csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
! Subprogram not used          csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
! Subprogram not used          
! Subprogram not used          csig12ne = p222*stress12_1(i,j) + ssig122 &
! Subprogram not used                   + p055*stress12_3(i,j)
! Subprogram not used          csig12nw = p222*stress12_2(i,j) + ssig121 &
! Subprogram not used                   + p055*stress12_4(i,j)
! Subprogram not used          csig12sw = p222*stress12_3(i,j) + ssig122 &
! Subprogram not used                   + p055*stress12_1(i,j)
! Subprogram not used          csig12se = p222*stress12_4(i,j) + ssig121 &
! Subprogram not used                   + p055*stress12_2(i,j)
! Subprogram not used 
! Subprogram not used          str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
! Subprogram not used          str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
! Subprogram not used          str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
! Subprogram not used          str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! for dF/dx (u momentum)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
! Subprogram not used          strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)
! Subprogram not used 
! Subprogram not used          ! northeast (i,j)
! Subprogram not used          str8(i,j,1) = -strp_tmp - strm_tmp - str12ew &
! Subprogram not used               + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne
! Subprogram not used 
! Subprogram not used          ! northwest (i+1,j)
! Subprogram not used          str8(i,j,2) = strp_tmp + strm_tmp - str12we &
! Subprogram not used               + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw
! Subprogram not used 
! Subprogram not used          strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
! Subprogram not used          strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)
! Subprogram not used 
! Subprogram not used          ! southeast (i,j+1)
! Subprogram not used          str8(i,j,3) = -strp_tmp - strm_tmp + str12ew &
! Subprogram not used               + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se
! Subprogram not used 
! Subprogram not used          ! southwest (i+1,j+1)
! Subprogram not used          str8(i,j,4) = strp_tmp + strm_tmp + str12we &
! Subprogram not used               + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! for dF/dy (v momentum)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
! Subprogram not used          strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)
! Subprogram not used 
! Subprogram not used          ! northeast (i,j)
! Subprogram not used          str8(i,j,5) = -strp_tmp + strm_tmp - str12ns &
! Subprogram not used               - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne
! Subprogram not used 
! Subprogram not used          ! southeast (i,j+1)
! Subprogram not used          str8(i,j,6) = strp_tmp - strm_tmp - str12sn &
! Subprogram not used               - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se
! Subprogram not used 
! Subprogram not used          strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
! Subprogram not used          strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)
! Subprogram not used 
! Subprogram not used          ! northwest (i+1,j)
! Subprogram not used          str8(i,j,7) = -strp_tmp + strm_tmp + str12ns &
! Subprogram not used               - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw
! Subprogram not used 
! Subprogram not used          ! southwest (i+1,j+1)
! Subprogram not used          str8(i,j,8) = strp_tmp - strm_tmp + str12sn &
! Subprogram not used               - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw
! Subprogram not used 
! Subprogram not used       endif                     ! phase
! Subprogram not used       enddo                     ! ij
! Subprogram not used 
! Subprogram not used       end subroutine stress

!=======================================================================
!BOP
!
! !IROUTINE: stepu - integrates mom eqn for u,v
!
! !INTERFACE:
!
! Subprogram not used       subroutine stepu (nx_block,   ny_block, &
! Subprogram not used                         icellu,               &
! Subprogram not used                         indxui,     indxuj,   &
! Subprogram not used                         aiu,        str8,      &
! Subprogram not used                         uocn,       vocn,     &
! Subprogram not used                         waterx,     watery,   &
! Subprogram not used                         forcex,     forcey,   &
! Subprogram not used                         umassdtei,  fm,       &
! Subprogram not used                         uarear,               &
! Subprogram not used                         strocnx,    strocny,  &
! Subprogram not used                         strintx,    strinty,  &
! Subprogram not used                         uvel,       vvel)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Calculation of the surface stresses
! Subprogram not used ! Integration of the momentum equation to find velocity (u,v)
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          icellu                ! total count when iceumask is true
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxui  , & ! compressed index in i-direction
! Subprogram not used          indxuj      ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
! Subprogram not used          aiu     , & ! ice fraction on u-grid
! Subprogram not used          waterx  , & ! for ocean stress calculation, x (m/s)
! Subprogram not used          watery  , & ! for ocean stress calculation, y (m/s)
! Subprogram not used          forcex  , & ! work array: combined atm stress and ocn tilt, x
! Subprogram not used          forcey  , & ! work array: combined atm stress and ocn tilt, y
! Subprogram not used          umassdtei,& ! mass of U-cell/dte (kg/m^2 s)
! Subprogram not used          uocn    , & ! ocean current, x-direction (m/s)
! Subprogram not used          vocn    , & ! ocean current, y-direction (m/s)
! Subprogram not used          fm      , & ! Coriolis param. * mass in U-cell (kg/s)
! Subprogram not used          uarear      ! 1/uarea
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,8), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          str8
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          uvel    , & ! x-component of velocity (m/s)
! Subprogram not used          vvel        ! y-component of velocity (m/s)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          strocnx , & ! ice-ocean stress, x-direction
! Subprogram not used          strocny , & ! ice-ocean stress, y-direction
! Subprogram not used          strintx , & ! divergence of internal ice stress, x (N/m^2)
! Subprogram not used          strinty     ! divergence of internal ice stress, y (N/m^2)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, ij
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          uold, vold        , & ! old-time uvel, vvel
! Subprogram not used          vrel              , & ! relative ice-ocean velocity
! Subprogram not used          cca,ccb,ab2,cc1,cc2,& ! intermediate variables
! Subprogram not used          taux, tauy            ! part of ocean stress term          
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! integrate the momentum equation
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do ij =1, icellu
! Subprogram not used          i = indxui(ij)
! Subprogram not used          j = indxuj(ij)
! Subprogram not used 
! Subprogram not used          uold = uvel(i,j)
! Subprogram not used          vold = vvel(i,j)
! Subprogram not used 
! Subprogram not used          ! (magnitude of relative ocean current)*rhow*drag*aice
! Subprogram not used          vrel = aiu(i,j)*dragw*sqrt((uocn(i,j) - uold)**2 + &
! Subprogram not used                                     (vocn(i,j) - vold)**2)  ! m/s
! Subprogram not used          ! ice/ocean stress
! Subprogram not used          taux = vrel*waterx(i,j) ! NOTE this is not the entire
! Subprogram not used          tauy = vrel*watery(i,j) ! ocn stress term
! Subprogram not used 
! Subprogram not used          ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
! Subprogram not used          cca = umassdtei(i,j) + vrel * cosw      ! alpha, kg/m^2 s
! Subprogram not used          ccb = fm(i,j)        + vrel * sinw      ! beta,  kg/m^2 s
! Subprogram not used          ab2 = cca**2 + ccb**2
! Subprogram not used 
! Subprogram not used          ! divergence of the internal stress tensor
! Subprogram not used          strintx(i,j) = uarear(i,j)* &
! Subprogram not used              (str8(i,j,1) + str8(i+1,j,2) + str8(i,j+1,3) + str8(i+1,j+1,4))
! Subprogram not used          strinty(i,j) = uarear(i,j)* &
! Subprogram not used              (str8(i,j,5) + str8(i,j+1,6) + str8(i+1,j,7) + str8(i+1,j+1,8))
! Subprogram not used 
! Subprogram not used          ! finally, the velocity components
! Subprogram not used          cc1 = strintx(i,j) + forcex(i,j) + taux &
! Subprogram not used              + umassdtei(i,j)*uold
! Subprogram not used          cc2 = strinty(i,j) + forcey(i,j) + tauy &
! Subprogram not used              + umassdtei(i,j)*vold
! Subprogram not used 
! Subprogram not used          uvel(i,j) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
! Subprogram not used          vvel(i,j) = (cca*cc2 - ccb*cc1) / ab2
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! ocean-ice stress for coupling
! Subprogram not used       ! here, strocn includes the factor of aice
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          strocnx(i,j) = taux
! Subprogram not used          strocny(i,j) = tauy
! Subprogram not used 
! Subprogram not used       enddo                     ! ij
! Subprogram not used 
! Subprogram not used       end subroutine stepu

!=======================================================================
!BOP
!
! !IROUTINE: evp_finish - calculates ice-ocean stress
!
! !INTERFACE:
!
! Subprogram not used       subroutine evp_finish (nx_block, ny_block, &
! Subprogram not used                              icellu,             &
! Subprogram not used                              indxui,   indxuj,   &
! Subprogram not used                              uvel,     vvel,     &
! Subprogram not used                              uocn,     vocn,     &
! Subprogram not used                              aiu,                &
! Subprogram not used                              strocnx,  strocny,  &
! Subprogram not used                              strocnxT, strocnyT) 
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Calculation of the ice-ocean stress.
! Subprogram not used ! ...the sign will be reversed later...
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          icellu                ! total count when iceumask is true
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxui  , & ! compressed index in i-direction
! Subprogram not used          indxuj      ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
! Subprogram not used          uvel    , & ! x-component of velocity (m/s)
! Subprogram not used          vvel    , & ! y-component of velocity (m/s)
! Subprogram not used          uocn    , & ! ocean current, x-direction (m/s)
! Subprogram not used          vocn    , &  ! ocean current, y-direction (m/s)
! Subprogram not used          aiu         ! ice fraction on u-grid
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          strocnx , & ! ice-ocean stress, x-direction
! Subprogram not used          strocny , & ! ice-ocean stress, y-direction
! Subprogram not used          strocnxT, & ! ice-ocean stress, x-direction
! Subprogram not used          strocnyT    ! ice-ocean stress, y-direction
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, ij
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: vrel
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          strocnxT(i,j) = c0
! Subprogram not used          strocnyT(i,j) = c0
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       ! ocean-ice stress for coupling
! Subprogram not used       do ij =1, icellu
! Subprogram not used          i = indxui(ij)
! Subprogram not used          j = indxuj(ij)
! Subprogram not used 
! Subprogram not used          vrel = dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + &
! Subprogram not used                            (vocn(i,j) - vvel(i,j))**2)  ! m/s
! Subprogram not used          strocnx(i,j) = strocnx(i,j) &
! Subprogram not used                       - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw) * aiu(i,j)
! Subprogram not used          strocny(i,j) = strocny(i,j) &
! Subprogram not used                       - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw) * aiu(i,j)
! Subprogram not used          ! Prepare to convert to T grid
! Subprogram not used          ! divide by aice for coupling
! Subprogram not used          strocnxT(i,j) = strocnx(i,j) / aiu(i,j)
! Subprogram not used          strocnyT(i,j) = strocny(i,j) / aiu(i,j)
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine evp_finish

!=======================================================================
!BOP
!
! !IROUTINE: principal_stress - computes principal stress for yield curve
!
! !INTERFACE:
!
! Subprogram not used       subroutine principal_stress(nx_block,   ny_block,  &
! Subprogram not used                                   stressp_1,  stressm_1, &
! Subprogram not used                                   stress12_1, prs_sig,   &
! Subprogram not used                                   sig1,       sig2)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Computes principal stresses for comparison with the theoretical
! Subprogram not used ! yield curve; northeast values
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block  ! block dimensions
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
! Subprogram not used          stressp_1 , & ! sigma11 + sigma22
! Subprogram not used          stressm_1 , & ! sigma11 - sigma22
! Subprogram not used          stress12_1, & ! sigma12
! Subprogram not used          prs_sig       ! replacement pressure, for stress calc
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
! Subprogram not used          sig1    , & ! principal stress component
! Subprogram not used          sig2        ! principal stress component
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: i, j
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          if (prs_sig(i,j) > puny) then
! Subprogram not used             sig1(i,j) = (p5*(stressp_1(i,j) &
! Subprogram not used                       + sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))) &
! Subprogram not used                       / prs_sig(i,j)
! Subprogram not used             sig2(i,j) = (p5*(stressp_1(i,j) &
! Subprogram not used                       - sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))) &
! Subprogram not used                       / prs_sig(i,j)
! Subprogram not used          else
! Subprogram not used             sig1(i,j) = spval_dbl
! Subprogram not used             sig2(i,j) = spval_dbl
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine principal_stress

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
