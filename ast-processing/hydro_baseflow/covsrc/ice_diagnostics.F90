!=======================================================================
!BOP
!
! !MODULE: ice_diagnostics - diagnostic information output during run
!
! !DESCRIPTION:
!
! Diagnostic information output during run
!
! !REVISION HISTORY:
!  SVN:$Id: ice_diagnostics.F90 52 2007-01-30 18:04:24Z eclare $
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_diagnostics
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_calendar, only: diagfreq, istep1, istep
      use ice_fileunits
      use ice_domain_size
!
!EOP
!
      implicit none
      save
	
      ! diagnostic output file
      character (len=char_len) :: diag_file

      ! point print data

      logical (kind=log_kind) :: &
         print_points     , & ! if true, print point data
         print_global         ! if true, print global data

      integer (kind=int_kind), parameter :: &
         npnt = 2             ! total number of points to be printed

      ! Set to true to identify unstable fast-moving ice.
      logical (kind=log_kind), parameter ::  &
         check_umax = .false. ! if true, check for speed > umax_stab

      real (kind=dbl_kind), parameter :: &
         umax_stab   = 1.0_dbl_kind , & ! ice speed threshold for instability (m/s)
         aice_extmin = 0.15_dbl_kind    ! min aice value for ice extent calc
 
      real (kind=dbl_kind), dimension(npnt) :: &
         latpnt           , & !  latitude of diagnostic points
         lonpnt               ! longitude of diagnostic points

      integer (kind=int_kind) :: &
         iindx            , & ! i index for points
         jindx            , & ! j index for points
         bindx                ! block index for points

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(npnt) :: &
         pdhi             , & ! change in mean ice thickness (m)
         pdhs             , & ! change in mean snow thickness (m)
         pde              , & ! change in ice and snow energy (J m-2)
         plat, plon           ! latitude, longitude of points

      integer (kind=int_kind), dimension(npnt) :: &
         piloc, pjloc, pbloc, pmloc  ! location of diagnostic points

      ! for hemispheric water and heat budgets
      real (kind=dbl_kind) :: &
         totmn            , & ! total ice/snow water mass (nh)
         totms            , & ! total ice/snow water mass (sh)
         totmin           , & ! total ice water mass (nh)
         totmis           , & ! total ice water mass (sh)
         toten            , & ! total ice/snow energy (J)
         totes                ! total ice/snow energy (J)
      real (kind=dbl_kind), dimension(n_aeromx) :: &
         totaeron         , & ! total aerosol mass
         totaeros             ! total aerosol mass

      ! printing info for routine print_state
      ! iblkp, ip, jp, mtask identify the grid cell to print
      character (char_len) :: plabel
      integer (kind=int_kind), parameter :: &
         check_step = 999999999, & ! begin printing at istep1=check_step
         iblkp = 1, &      ! block number
         ip = 3, &         ! i index
         jp = 5, &         ! j index
         mtask = 0         ! my_task

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: runtime_diags - writes max,min,global sums to standard out
!
! !INTERFACE:
!
! Subprogram not used       subroutine runtime_diags (dt)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Writes diagnostic info (max, min, global sums, etc) to standard out
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors: Elizabeth C. Hunke, LANL
! Subprogram not used !          Bruce P. Briegleb, NCAR
! Subprogram not used !          Cecilia M. Bitz, UW
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_broadcast
! Subprogram not used       use ice_global_reductions
! Subprogram not used       use ice_blocks
! Subprogram not used       use ice_domain
! Subprogram not used !MH      use ice_domain_size
! Subprogram not used       use ice_flux
! Subprogram not used       use ice_state
! Subprogram not used       use ice_grid, only: lmask_n, lmask_s, tarean, tareas, grid_type
! Subprogram not used       use ice_therm_vertical, only: calc_Tsfc
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       use ice_prescribed_mod, only : prescribed_ice
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          dt      ! time step
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, k, n, ii,jj, iblk
! Subprogram not used 
! Subprogram not used       ! hemispheric state quantities
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          umaxn,   hmaxn,   shmaxn,    arean,   snwmxn, extentn, &
! Subprogram not used          umaxs,   hmaxs,   shmaxs,    areas,   snwmxs, extents, &
! Subprogram not used          etotn,   mtotn,   micen,     msnwn,   pmaxn,  ketotn, &
! Subprogram not used          etots,   mtots,   mices,     msnws,   pmaxs,  ketots, &
! Subprogram not used          urmsn,   albtotn, arean_alb, &
! Subprogram not used          urmss,   albtots, areas_alb
! Subprogram not used 
! Subprogram not used       ! hemispheric flux quantities
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          rnn, snn, frzn,  hnetn, fhocnn, fhatmn,  fhfrzn, &
! Subprogram not used          rns, sns, frzs,  hnets, fhocns, fhatms,  fhfrzs, &
! Subprogram not used          sfsaltn, sfreshn, evpn, fluxn , delmxn,  delmin, &
! Subprogram not used          sfsalts, sfreshs, evps, fluxs , delmxs,  delmis, &
! Subprogram not used          delein, werrn, herrn, msltn, delmsltn, serrn, &
! Subprogram not used          deleis, werrs, herrs, mslts, delmslts, serrs, &
! Subprogram not used          ftmp,faeron,faeros,fsootn,fsoots
! Subprogram not used 
! Subprogram not used ! MH for aerosol diagnostics
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used         kaero, naero
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used         aeromx1n, aeromx1s, aeromx2n, aeromx2s, &
! Subprogram not used         aeromx3n, aeromx3s, aoermx4, &
! Subprogram not used         aerototn, aerotots     !MH
! Subprogram not used 
! Subprogram not used       ! fields at diagnostic points
! Subprogram not used       real (kind=dbl_kind), dimension(npnt) :: &
! Subprogram not used          paice, pTair, pQa, pfsnow, pfrain, pfsw, pflw, & 
! Subprogram not used          pTsfc, pevap, pfswabs, pflwout, pflat, pfsens, &
! Subprogram not used          pfsurf, pfcondtop, psst,  pTf, hiavg, hsavg, pfhocn, &
! Subprogram not used          pmeltt, pmeltb, pmeltl, psnoice, pfrazil, pcongel
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
! Subprogram not used          work1, work2
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! state of the ice
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! hemispheric quantities
! Subprogram not used 
! Subprogram not used       ! total ice area
! Subprogram not used       arean = global_sum(aice, distrb_info, field_loc_center, tarean)
! Subprogram not used       areas = global_sum(aice, distrb_info, field_loc_center, tareas)
! Subprogram not used       arean = arean * m2_to_km2
! Subprogram not used       areas = areas * m2_to_km2
! Subprogram not used 
! Subprogram not used       ! ice extent (= area of grid cells with aice > aice_extmin)
! Subprogram not used       work1(:,:,:) = c0
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             if (aice(i,j,iblk) >= aice_extmin) work1(i,j,iblk) = c1
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used       extentn = global_sum(work1, distrb_info, field_loc_center, &
! Subprogram not used                            tarean)
! Subprogram not used       extents = global_sum(work1, distrb_info, field_loc_center, &
! Subprogram not used                            tareas)
! Subprogram not used       extentn = extentn * m2_to_km2
! Subprogram not used       extents = extents * m2_to_km2
! Subprogram not used 
! Subprogram not used       ! total ice volume
! Subprogram not used       shmaxn = global_sum(vice, distrb_info, field_loc_center, tarean)
! Subprogram not used       shmaxs = global_sum(vice, distrb_info, field_loc_center, tareas)
! Subprogram not used 
! Subprogram not used       ! total snow volume
! Subprogram not used       snwmxn = global_sum(vsno, distrb_info, field_loc_center, tarean)
! Subprogram not used       snwmxs = global_sum(vsno, distrb_info, field_loc_center, tareas)
! Subprogram not used 
! Subprogram not used       ! total ice-snow kinetic energy
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             work1(i,j,iblk) = p5 &
! Subprogram not used                            * (rhos*vsno(i,j,iblk) + rhoi*vice(i,j,iblk)) &
! Subprogram not used                            * (uvel(i,j,iblk)**2 + vvel(i,j,iblk)**2)
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used       ketotn = global_sum(work1, distrb_info, field_loc_center, tarean)
! Subprogram not used       ketots = global_sum(work1, distrb_info, field_loc_center, tareas)
! Subprogram not used 
! Subprogram not used       ! rms ice speed
! Subprogram not used       urmsn = c2*ketotn/(rhoi*shmaxn + rhos*snwmxn + puny)
! Subprogram not used       if (urmsn > puny) then
! Subprogram not used          urmsn = sqrt(urmsn)
! Subprogram not used       else
! Subprogram not used          urmsn = c0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       urmss = c2*ketots/(rhoi*shmaxs + rhos*snwmxs + puny)
! Subprogram not used       if (urmss > puny) then
! Subprogram not used          urmss = sqrt(urmss)
! Subprogram not used       else
! Subprogram not used          urmss = c0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       ! average ice albedo
! Subprogram not used       ! mask out cells where sun is below horizon (for delta-Eddington)
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             work1(i,j,iblk) = alvdr(i,j,iblk)*awtvdr &
! Subprogram not used                             + alidr(i,j,iblk)*awtidr &
! Subprogram not used                             + alvdf(i,j,iblk)*awtvdf &
! Subprogram not used                             + alidf(i,j,iblk)*awtidf
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO 
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             if (coszen(i,j,iblk) > puny) then
! Subprogram not used                work2(i,j,iblk) = tarean(i,j,iblk)
! Subprogram not used             else
! Subprogram not used                work2(i,j,iblk) = c0
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO 
! Subprogram not used       
! Subprogram not used       arean_alb = global_sum(aice, distrb_info, field_loc_center, work2)      
! Subprogram not used 
! Subprogram not used       albtotn = global_sum_prod(aice, work1, distrb_info, &
! Subprogram not used                                 field_loc_center, work2)
! Subprogram not used 
! Subprogram not used       if (arean_alb > c0) then
! Subprogram not used          albtotn = albtotn / arean_alb
! Subprogram not used       else
! Subprogram not used          albtotn = c0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             if (coszen(i,j,iblk) > puny) then
! Subprogram not used                work2(i,j,iblk) = tareas(i,j,iblk)
! Subprogram not used             else
! Subprogram not used                work2(i,j,iblk) = c0
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO 
! Subprogram not used 
! Subprogram not used       areas_alb = global_sum(aice, distrb_info, field_loc_center, work2)      
! Subprogram not used 
! Subprogram not used       albtots = global_sum_prod(aice, work1, distrb_info, &
! Subprogram not used                                 field_loc_center, work2)
! Subprogram not used 
! Subprogram not used       if (areas_alb > c0) then
! Subprogram not used          albtots = albtots / areas_alb
! Subprogram not used       else
! Subprogram not used          albtots = c0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       ! maximum ice volume (= mean thickness including open water)
! Subprogram not used       hmaxn = global_maxval(vice, distrb_info, lmask_n)
! Subprogram not used       hmaxs = global_maxval(vice, distrb_info, lmask_s)
! Subprogram not used 
! Subprogram not used ! MH put in aerosol diagnostics
! Subprogram not used       if (tr_aero) then
! Subprogram not used          ! aerosols
! Subprogram not used         do naero=1,n_aero
! Subprogram not used          faeron = global_sum_prod(faero(:,:,naero,:), aice_init, distrb_info, &
! Subprogram not used                                field_loc_center, tarean)
! Subprogram not used          faeros = global_sum_prod(faero(:,:,naero,:), aice_init, distrb_info, &
! Subprogram not used                                field_loc_center, tareas)
! Subprogram not used          faeron = faeron*dt
! Subprogram not used          faeros = faeros*dt
! Subprogram not used 
! Subprogram not used          fsootn = global_sum_prod(fsoot(:,:,naero,:), aice, distrb_info, &
! Subprogram not used                                field_loc_center, tarean)
! Subprogram not used          fsoots = global_sum_prod(fsoot(:,:,naero,:), aice, distrb_info, &
! Subprogram not used                                field_loc_center, tareas)
! Subprogram not used          fsootn = fsootn*dt
! Subprogram not used          fsoots = fsoots*dt
! Subprogram not used 
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used            do j = 1, ny_block
! Subprogram not used            do i = 1, nx_block
! Subprogram not used             work1(i,j,iblk) = trcr(i,j,nt_aero  +4*(naero-1),iblk)  *vsno(i,j,iblk) &
! Subprogram not used                             + trcr(i,j,nt_aero+1+4*(naero-1),iblk)*vsno(i,j,iblk) &
! Subprogram not used                             + trcr(i,j,nt_aero+2+4*(naero-1),iblk)*vice(i,j,iblk) &
! Subprogram not used                             + trcr(i,j,nt_aero+3+4*(naero-1),iblk)*vice(i,j,iblk)
! Subprogram not used            enddo
! Subprogram not used            enddo
! Subprogram not used          enddo
! Subprogram not used          aerototn= global_sum(work1, distrb_info, field_loc_center, tarean)
! Subprogram not used          aerotots= global_sum(work1, distrb_info, field_loc_center, tareas)
! Subprogram not used          aeromx1n = global_maxval(work1, distrb_info, lmask_n)
! Subprogram not used          aeromx1s = global_maxval(work1, distrb_info, lmask_s)
! Subprogram not used          if (my_task == master_task) then
! Subprogram not used           write(nu_diag,*) 'aero: ',naero,' faero         : ',&
! Subprogram not used                 faeron, faeros
! Subprogram not used           write(nu_diag,*) 'aero: ',naero,' fsoot         : ',&
! Subprogram not used                 fsootn, fsoots
! Subprogram not used           write(nu_diag,*) 'aero: ',naero,' faero-fsoot   : ',&
! Subprogram not used                 faeron-fsootn, faeros-fsoots
! Subprogram not used           write(nu_diag,*) 'aero: ',naero,' aerotot       : ',&
! Subprogram not used                 aerototn, aerotots
! Subprogram not used           write(nu_diag,*) 'aero: ',naero,' aerotot change: ',&
! Subprogram not used                 aerototn-totaeron(naero), aerotots-totaeros(naero)
! Subprogram not used           write(nu_diag,*) 'aero: ',naero,' aeromax agg: ',&
! Subprogram not used                 aeromx1n,aeromx1s
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used !         do kaero=1,ncat
! Subprogram not used !          do iblk = 1, nblocks
! Subprogram not used !           do j = 1, ny_block
! Subprogram not used !           do i = 1, nx_block
! Subprogram not used !            work1(i,j,iblk) = trcrn(i,j,nt_aero,kaero,iblk)
! Subprogram not used !           enddo
! Subprogram not used !           enddo
! Subprogram not used !          enddo
! Subprogram not used !          aeromx1n = global_maxval(work1, distrb_info, lmask_n)
! Subprogram not used !          aeromx1s = global_maxval(work1, distrb_info, lmask_s)
! Subprogram not used !          if (my_task == master_task) &
! Subprogram not used !           write(nu_diag,*) 'MH aeromx1s: ',aeromx1n,aeromx1s,kaero
! Subprogram not used !         enddo
! Subprogram not used 
! Subprogram not used !       do iblk = 1, nblocks
! Subprogram not used !         do j = 1, ny_block
! Subprogram not used !         do i = 1, nx_block
! Subprogram not used !            work1(i,j,iblk) = trcrn(i,j,nt_aero+1,1,iblk)
! Subprogram not used !         enddo
! Subprogram not used !         enddo
! Subprogram not used !       enddo
! Subprogram not used !       aeromx2n = global_maxval(work1, distrb_info, lmask_n)
! Subprogram not used !       write(nu_diag,*) 'MH aeromx2n: ',aeromx2n
! Subprogram not used !       aeromx2s = global_maxval(work1, distrb_info, lmask_s)
! Subprogram not used !       write(nu_diag,*) 'MH aeromx2s: ',aeromx2s
! Subprogram not used !
! Subprogram not used !       do iblk = 1, nblocks
! Subprogram not used !         do j = 1, ny_block
! Subprogram not used !         do i = 1, nx_block
! Subprogram not used !            work1(i,j,iblk) = trcrn(i,j,nt_aero+2,1,iblk)
! Subprogram not used !         enddo
! Subprogram not used !         enddo
! Subprogram not used !       enddo
! Subprogram not used !       aeromx3n = global_maxval(work1, distrb_info, lmask_n)
! Subprogram not used !       write(nu_diag,*) 'MH aeromx2n: ',aeromx3n
! Subprogram not used !       aeromx3s = global_maxval(work1, distrb_info, lmask_s)
! Subprogram not used !       write(nu_diag,*) 'MH aeromx2s: ',aeromx3s
! Subprogram not used         enddo ! n_aero
! Subprogram not used       endif  ! tr_aero
! Subprogram not used     
! Subprogram not used       ! maximum ice speed
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             work1(i,j,iblk) = sqrt(uvel(i,j,iblk)**2 &
! Subprogram not used                                  + vvel(i,j,iblk)**2)
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       umaxn = global_maxval(work1, distrb_info, lmask_n)
! Subprogram not used       umaxs = global_maxval(work1, distrb_info, lmask_s)
! Subprogram not used 
! Subprogram not used       ! Write warning message if ice speed is too big
! Subprogram not used       ! (Ice speeds of ~1 m/s or more usually indicate instability)
! Subprogram not used 
! Subprogram not used       if (check_umax) then
! Subprogram not used       	 if (umaxn > umax_stab) then
! Subprogram not used             do iblk = 1, nblocks
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                if (abs(work1(i,j,iblk) - umaxn) < puny) then
! Subprogram not used                   write(nu_diag,*) ' '
! Subprogram not used                   write(nu_diag,*) 'Warning, large ice speed'
! Subprogram not used                   write(nu_diag,*) 'my_task, iblk, i, j, umaxn:', &
! Subprogram not used                                     my_task, iblk, i, j, umaxn
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          elseif (umaxs > umax_stab) then
! Subprogram not used             do iblk = 1, nblocks
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                if (abs(work1(i,j,iblk) - umaxs) < puny) then
! Subprogram not used                   write(nu_diag,*) ' '
! Subprogram not used                   write(nu_diag,*) 'Warning, large ice speed'
! Subprogram not used                   write(nu_diag,*) 'my_task, iblk, i, j, umaxs:', &
! Subprogram not used                                     my_task, iblk, i, j, umaxs
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          endif   ! umax
! Subprogram not used       endif      ! check_umax
! Subprogram not used 
! Subprogram not used       ! maximum ice strength
! Subprogram not used 
! Subprogram not used       pmaxn = global_maxval(strength, distrb_info, lmask_n)
! Subprogram not used       pmaxs = global_maxval(strength, distrb_info, lmask_s)
! Subprogram not used 
! Subprogram not used       pmaxn = pmaxn / c1000   ! convert to kN/m
! Subprogram not used       pmaxs = pmaxs / c1000 
! Subprogram not used 
! Subprogram not used       if (print_global) then
! Subprogram not used 
! Subprogram not used          ! total ice/snow internal energy
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                work1(i,j,iblk) = esno(i,j,iblk) + eice(i,j,iblk)
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          etotn = global_sum(work1, distrb_info, &
! Subprogram not used                             field_loc_center, tarean)
! Subprogram not used          etots = global_sum(work1, distrb_info, &
! Subprogram not used                             field_loc_center, tareas)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! various fluxes
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! evap, fsens, and flwout need to be multiplied by aice because
! Subprogram not used       ! regrettably they have been divided by aice for the coupler
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          ! evaporation
! Subprogram not used 
! Subprogram not used          evpn = global_sum_prod(evap, aice, distrb_info, &
! Subprogram not used                                 field_loc_center, tarean)
! Subprogram not used          evps = global_sum_prod(evap, aice, distrb_info, &
! Subprogram not used                                 field_loc_center, tareas)
! Subprogram not used          evpn = evpn*dt
! Subprogram not used          evps = evps*dt
! Subprogram not used 
! Subprogram not used          ! salt flux
! Subprogram not used          sfsaltn = global_sum(fsalt_gbm, distrb_info, &
! Subprogram not used                               field_loc_center, tarean)
! Subprogram not used          sfsalts = global_sum(fsalt_gbm, distrb_info, &
! Subprogram not used                               field_loc_center, tareas)
! Subprogram not used          sfsaltn = sfsaltn*dt
! Subprogram not used          sfsalts = sfsalts*dt
! Subprogram not used 
! Subprogram not used          ! fresh water flux
! Subprogram not used          sfreshn = global_sum(fresh_gbm, distrb_info, &
! Subprogram not used                               field_loc_center, tarean)
! Subprogram not used          sfreshs = global_sum(fresh_gbm, distrb_info, &
! Subprogram not used                               field_loc_center, tareas)
! Subprogram not used          sfreshn = sfreshn*dt
! Subprogram not used          sfreshs = sfreshs*dt
! Subprogram not used 
! Subprogram not used          ! ocean heat
! Subprogram not used          ! Note: fswthru not included because it does not heat ice
! Subprogram not used          fhocnn = global_sum(fhocn_gbm, distrb_info, &
! Subprogram not used                              field_loc_center, tarean)
! Subprogram not used          fhocns = global_sum(fhocn_gbm, distrb_info, &
! Subprogram not used                              field_loc_center, tareas)
! Subprogram not used 
! Subprogram not used          ! latent heat
! Subprogram not used          ! You may be wondering, where is the latent heat flux?
! Subprogram not used          ! It is not included here because it cancels with
! Subprogram not used          ! the evaporative flux times the enthalpy of the
! Subprogram not used          ! ice/snow that evaporated.
! Subprogram not used 
! Subprogram not used          ! atmo heat flux
! Subprogram not used          ! Note: flwout includes the reflected longwave down, needed by the
! Subprogram not used          !  atmosphere as an upwards radiative boundary condition.
! Subprogram not used          ! Also note: fswabs includes solar radiation absorbed in ocean,
! Subprogram not used          !  which must be subtracted here.
! Subprogram not used 
! Subprogram not used          if (calc_Tsfc) then
! Subprogram not used 
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used             do iblk = 1, nblocks
! Subprogram not used                do j = 1, ny_block
! Subprogram not used                do i = 1, nx_block
! Subprogram not used                   work1(i,j,iblk) = &
! Subprogram not used                                (fswabs(i,j,iblk) - fswthru(i,j,iblk) &
! Subprogram not used                               + flwout(i,j,iblk)                          &
! Subprogram not used                               + fsens (i,j,iblk)) * aice(i,j,iblk)        &
! Subprogram not used                               + flw   (i,j,iblk)  * aice_init(i,j,iblk)
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          else   ! fsurf is computed by atmosphere model
! Subprogram not used 
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used             do iblk = 1, nblocks
! Subprogram not used                do j = 1, ny_block
! Subprogram not used                do i = 1, nx_block
! Subprogram not used                   work1(i,j,iblk) = &
! Subprogram not used                              (fsurf(i,j,iblk) - flat(i,j,iblk)) &
! Subprogram not used                               * aice(i,j,iblk)
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used          endif     ! calc_Tsfc
! Subprogram not used 
! Subprogram not used          fhatmn = global_sum(work1, distrb_info, &
! Subprogram not used                              field_loc_center, tarean)
! Subprogram not used          fhatms = global_sum(work1, distrb_info, &
! Subprogram not used                              field_loc_center, tareas)
! Subprogram not used   
! Subprogram not used          ! freezing potential
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                work1(i,j,iblk) = max(c0,frzmlt(i,j,iblk))
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used          fhfrzn = global_sum(work1, distrb_info, &
! Subprogram not used                              field_loc_center, tarean)
! Subprogram not used          fhfrzs = global_sum(work1, distrb_info, &
! Subprogram not used                              field_loc_center, tareas)
! Subprogram not used 
! Subprogram not used          ! rain
! Subprogram not used          rnn = global_sum_prod(frain, aice_init, distrb_info, &
! Subprogram not used                                field_loc_center, tarean)
! Subprogram not used          rns = global_sum_prod(frain, aice_init, distrb_info, &
! Subprogram not used                                field_loc_center, tareas)
! Subprogram not used          rnn = rnn*dt
! Subprogram not used          rns = rns*dt
! Subprogram not used 
! Subprogram not used          ! snow
! Subprogram not used          snn = global_sum_prod(fsnow, aice_init, distrb_info, &
! Subprogram not used                                field_loc_center, tarean)
! Subprogram not used          sns = global_sum_prod(fsnow, aice_init, distrb_info, &
! Subprogram not used                                field_loc_center, tareas)
! Subprogram not used          snn = snn*dt
! Subprogram not used          sns = sns*dt
! Subprogram not used 
! Subprogram not used          ! frazil ice growth !! should not be multiplied by aice
! Subprogram not used          ! m/step->kg/m^2/s
! Subprogram not used          work1(:,:,:) = frazil(:,:,:)*rhoi/dt
! Subprogram not used          frzn = global_sum(work1, distrb_info, &
! Subprogram not used                            field_loc_center, tarean)
! Subprogram not used          frzs = global_sum(work1, distrb_info, field_loc_center, &
! Subprogram not used                            tareas)
! Subprogram not used          frzn = frzn*dt
! Subprogram not used          frzs = frzs*dt
! Subprogram not used 
! Subprogram not used          ! ice and snow mass
! Subprogram not used          micen = rhoi*shmaxn
! Subprogram not used          msnwn = rhos*snwmxn
! Subprogram not used          mices = rhoi*shmaxs
! Subprogram not used          msnws = rhos*snwmxs
! Subprogram not used 
! Subprogram not used          mtotn = micen + msnwn
! Subprogram not used          mtots = mices + msnws
! Subprogram not used   
! Subprogram not used          ! mass change since beginning of time step
! Subprogram not used          delmin = mtotn - totmn
! Subprogram not used          delmis = mtots - totms
! Subprogram not used 
! Subprogram not used          ! ice mass change including frazil ice formation
! Subprogram not used          delmxn = micen - totmin
! Subprogram not used          delmxs = mices - totmis
! Subprogram not used          if (.not. update_ocn_f) then
! Subprogram not used            ! ice mass change excluding frazil ice formation
! Subprogram not used            delmxn = delmxn - frzn
! Subprogram not used            delmxs = delmxs - frzs
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          ! total water flux
! Subprogram not used          fluxn  = c0
! Subprogram not used          fluxs  = c0
! Subprogram not used          if( arean > c0) then
! Subprogram not used            ! water associated with frazil ice included in fresh
! Subprogram not used            fluxn = rnn + snn + evpn - sfreshn
! Subprogram not used            if (.not. update_ocn_f) then
! Subprogram not used              fluxn = fluxn + frzn
! Subprogram not used            endif
! Subprogram not used          endif
! Subprogram not used          if( areas > c0) then
! Subprogram not used            ! water associated with frazil ice included in fresh
! Subprogram not used            fluxs = rns + sns + evps - sfreshs
! Subprogram not used            if (.not. update_ocn_f) then
! Subprogram not used              fluxs = fluxs + frzs
! Subprogram not used            endif
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          werrn = (fluxn-delmin)/(mtotn+c1)
! Subprogram not used          werrs = (fluxs-delmis)/(mtots+c1)
! Subprogram not used 
! Subprogram not used          ! energy change
! Subprogram not used          delein = etotn - toten
! Subprogram not used          deleis = etots - totes
! Subprogram not used 
! Subprogram not used          fhatmn = fhatmn + ( - snn * Lfresh + evpn * Lvap ) / dt
! Subprogram not used          fhatms = fhatms + ( - sns * Lfresh + evps * Lvap ) / dt
! Subprogram not used 
! Subprogram not used          hnetn = (fhatmn - fhocnn - fhfrzn) * dt
! Subprogram not used          hnets = (fhatms - fhocns - fhfrzs) * dt
! Subprogram not used 
! Subprogram not used          herrn = (hnetn - delein) / (etotn - c1)
! Subprogram not used          herrs = (hnets - deleis) / (etots - c1)
! Subprogram not used 
! Subprogram not used          ! salt mass
! Subprogram not used          msltn = micen*ice_ref_salinity*p001
! Subprogram not used          mslts = mices*ice_ref_salinity*p001
! Subprogram not used 
! Subprogram not used          ! change in salt mass
! Subprogram not used          delmsltn = delmxn*ice_ref_salinity*p001
! Subprogram not used          delmslts = delmxs*ice_ref_salinity*p001
! Subprogram not used 
! Subprogram not used          ! salt error
! Subprogram not used          serrn = (sfsaltn + delmsltn) / (msltn + c1)
! Subprogram not used          serrs = (sfsalts + delmslts) / (mslts + c1)
! Subprogram not used 
! Subprogram not used       endif                     ! print_global
! Subprogram not used 
! Subprogram not used       if (print_points) then
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! state of the ice and associated fluxes for 2 defined points
! Subprogram not used       ! NOTE these are computed for the last timestep only (not avg)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do n = 1, npnt
! Subprogram not used             if (my_task == pmloc(n)) then
! Subprogram not used                i = piloc(n)
! Subprogram not used                j = pjloc(n)
! Subprogram not used                iblk = pbloc(n)
! Subprogram not used 
! Subprogram not used                pTair(n) = Tair(i,j,iblk) - Tffresh ! air temperature
! Subprogram not used                pQa(n) = Qa(i,j,iblk)               ! specific humidity
! Subprogram not used                pfsnow(n) = fsnow(i,j,iblk)*dt/rhos ! snowfall
! Subprogram not used                pfrain(n) = frain(i,j,iblk)*dt/rhow ! rainfall
! Subprogram not used                pfsw(n) = fsw(i,j,iblk)             ! shortwave radiation
! Subprogram not used                pflw(n) = flw(i,j,iblk)             ! longwave radiation
! Subprogram not used                paice(n) = aice(i,j,iblk)           ! ice area
! Subprogram not used                
! Subprogram not used                hiavg(n) = c0                       ! avg snow/ice thickness
! Subprogram not used                hsavg(n) = c0
! Subprogram not used                if (paice(n) /= c0) then
! Subprogram not used                   hiavg(n) = vice(i,j,iblk)/paice(n)
! Subprogram not used                   hsavg(n) = vsno(i,j,iblk)/paice(n)
! Subprogram not used                endif
! Subprogram not used                pTsfc(n) = trcr(i,j,nt_Tsfc,iblk)   ! ice/snow sfc temperature
! Subprogram not used                pevap(n) = evap(i,j,iblk)*dt/rhoi   ! sublimation/condensation
! Subprogram not used                pfswabs(n) = fswabs(i,j,iblk)       ! absorbed solar flux
! Subprogram not used                pflwout(n) = flwout(i,j,iblk)       ! outward longwave flux
! Subprogram not used                pflat(n) = flat(i,j,iblk)           ! latent heat flux
! Subprogram not used                pfsens(n) = fsens(i,j,iblk)         ! sensible heat flux
! Subprogram not used                pfsurf(n) = fsurf(i,j,iblk)         ! total sfc heat flux
! Subprogram not used                pfcondtop(n) = fcondtop(i,j,iblk)   ! top sfc cond flux
! Subprogram not used                pmeltt(n) = meltt(i,j,iblk)         ! top melt
! Subprogram not used                pmeltb(n) = meltb(i,j,iblk)         ! bottom melt
! Subprogram not used                pmeltl(n) = meltl(i,j,iblk)         ! lateral melt
! Subprogram not used                psnoice(n) = snoice(i,j,iblk)       ! snow ice
! Subprogram not used                pfrazil(n) = frazil(i,j,iblk)       ! frazil ice
! Subprogram not used                pcongel(n) = congel(i,j,iblk)       ! congelation ice
! Subprogram not used                pdhi(n) = vice(i,j,iblk) - pdhi(n)  ! ice thickness change
! Subprogram not used                pdhs(n) = vsno(i,j,iblk) - pdhs(n)  ! snow thickness change
! Subprogram not used                pde(n) = -(eice(i,j,iblk) &         ! ice/snow energy change 
! Subprogram not used                         + esno(i,j,iblk) - pde(n)) / dt
! Subprogram not used                psst(n) = sst(i,j,iblk)             ! sea surface temperature
! Subprogram not used                pTf(n) = Tf(i,j,iblk)               ! freezing temperature
! Subprogram not used                pfhocn(n) = -fhocn(i,j,iblk)   ! ocean heat used by ice
! Subprogram not used 
! Subprogram not used             endif  ! my_task = pmloc
! Subprogram not used 
! Subprogram not used             call broadcast_scalar(pTair    (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pQa      (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pfsnow   (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pfrain   (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pfsw     (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pflw     (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(paice    (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(hsavg    (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(hiavg    (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pTsfc    (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pevap    (n), pmloc(n))             
! Subprogram not used             call broadcast_scalar(pfswabs  (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pflwout  (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pflat    (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pfsens   (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pfsurf   (n), pmloc(n))
! Subprogram not used             call broadcast_scalar(pfcondtop(n), pmloc(n))
! Subprogram not used             call broadcast_scalar(pmeltt   (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pmeltb   (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pmeltl   (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(psnoice  (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pfrazil  (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pcongel  (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pdhi     (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pdhs     (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pde      (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(psst     (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pTf      (n), pmloc(n)) 
! Subprogram not used             call broadcast_scalar(pfhocn   (n), pmloc(n))
! Subprogram not used             
! Subprogram not used          enddo                  ! npnt
! Subprogram not used       endif                     ! print_points
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! start spewing
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used        if (grid_type == 'panarctic') then   ! Arctic only
! Subprogram not used         write (nu_diag,799) 'Arctic diagnostics'
! Subprogram not used         write (nu_diag,801) 'total ice area  (km^2) = ',arean
! Subprogram not used         write (nu_diag,801) 'total ice extent(km^2) = ',extentn
! Subprogram not used         write (nu_diag,801) 'total ice volume (m^3) = ',shmaxn
! Subprogram not used         write (nu_diag,801) 'total snw volume (m^3) = ',snwmxn
! Subprogram not used         write (nu_diag,801) 'tot kinetic energy (J) = ',ketotn
! Subprogram not used         write (nu_diag,800) 'rms ice speed    (m/s) = ',urmsn
! Subprogram not used         write (nu_diag,800) 'average albedo         = ',albtotn
! Subprogram not used         write (nu_diag,800) 'max ice volume     (m) = ',hmaxn
! Subprogram not used         write (nu_diag,800) 'max ice speed    (m/s) = ',umaxn
! Subprogram not used         write (nu_diag,900) 'max strength    (kN/m) = ',pmaxn
! Subprogram not used 
! Subprogram not used         if (print_global) then  ! global diags for conservations checks
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         if (prescribed_ice) then
! Subprogram not used           write (nu_diag,*) '----------------------------'
! Subprogram not used           write (nu_diag,*)   'This is the prescribed ice option.'
! Subprogram not used           write (nu_diag,*)   'Heat and water will not be conserved.'   
! Subprogram not used         else
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           write (nu_diag,*) '----------------------------'
! Subprogram not used           write (nu_diag,801) 'arwt rain h2o kg in dt = ',rnn
! Subprogram not used           write (nu_diag,801) 'arwt snow h2o kg in dt = ',snn
! Subprogram not used           write (nu_diag,801) 'arwt evap h2o kg in dt = ',evpn
! Subprogram not used           write (nu_diag,801) 'arwt frzl h2o kg in dt = ',frzn
! Subprogram not used           write (nu_diag,801) 'arwt frsh h2o kg in dt = ',sfreshn
! Subprogram not used          
! Subprogram not used           write (nu_diag,801) 'arwt ice mass (kg)     = ',micen
! Subprogram not used           write (nu_diag,801) 'arwt snw mass (kg)     = ',msnwn
! Subprogram not used 
! Subprogram not used           write (nu_diag,801) 'arwt tot mass (kg)     = ',mtotn
! Subprogram not used           write (nu_diag,801) 'arwt tot mass chng(kg) = ',delmin
! Subprogram not used           write (nu_diag,801) 'arwt water flux        = ',fluxn
! Subprogram not used           if (update_ocn_f) then
! Subprogram not used             write (nu_diag,*) '(=rain+snow+evap-fresh)  '
! Subprogram not used           else
! Subprogram not used             write (nu_diag,*) '(=rain+snow+evap+frzl-fresh)  '
! Subprogram not used           endif
! Subprogram not used           write (nu_diag,801) 'water flux error       = ',werrn
! Subprogram not used 
! Subprogram not used          endif                    ! prescribed_ice
! Subprogram not used 
! Subprogram not used          write (nu_diag,*) '----------------------------'
! Subprogram not used          write (nu_diag,801) 'arwt atm heat flux (W) = ',fhatmn
! Subprogram not used          write (nu_diag,801) 'arwt ocn heat flux (W) = ',fhocnn
! Subprogram not used          write (nu_diag,801) 'arwt frzl heat flux(W) = ',fhfrzn
! Subprogram not used          write (nu_diag,801) 'arwt tot energy    (J) = ',etotn
! Subprogram not used          write (nu_diag,801) 'arwt net heat      (J) = ',hnetn
! Subprogram not used          write (nu_diag,801) 'arwt tot energy chng(J)= ',delein
! Subprogram not used          write (nu_diag,801) 'arwt heat error        = ',herrn
! Subprogram not used        
! Subprogram not used          write (nu_diag,*) '----------------------------'
! Subprogram not used          write (nu_diag,801) 'arwt salt mass (kg)    = ',msltn
! Subprogram not used          write (nu_diag,801) 'arwt salt mass chng(kg)= ',delmsltn
! Subprogram not used          write (nu_diag,801) 'arwt salt flx in dt(kg)= ',sfsaltn
! Subprogram not used          write (nu_diag,801) 'arwt salt flx error    = ',serrn
! Subprogram not used          write (nu_diag,*) '----------------------------'
! Subprogram not used 
! Subprogram not used         endif                     ! print_global
! Subprogram not used 
! Subprogram not used        else  ! global grid
! Subprogram not used 
! Subprogram not used         write(nu_diag,899) 'Arctic','Antarctic'
! Subprogram not used 
! Subprogram not used         write(nu_diag,901) 'total ice area  (km^2) = ',arean,  areas
! Subprogram not used         write(nu_diag,901) 'total ice extent(km^2) = ',extentn,extents
! Subprogram not used         write(nu_diag,901) 'total ice volume (m^3) = ',shmaxn, shmaxs
! Subprogram not used         write(nu_diag,901) 'total snw volume (m^3) = ',snwmxn, snwmxs
! Subprogram not used         write(nu_diag,901) 'tot kinetic energy (J) = ',ketotn, ketots
! Subprogram not used         write(nu_diag,900) 'rms ice speed    (m/s) = ',urmsn,  urmss
! Subprogram not used         write(nu_diag,900) 'average albedo         = ',albtotn,albtots
! Subprogram not used         write(nu_diag,900) 'max ice volume     (m) = ',hmaxn,  hmaxs
! Subprogram not used         write(nu_diag,900) 'max ice speed    (m/s) = ',umaxn,  umaxs
! Subprogram not used         write(nu_diag,900) 'max strength    (kN/m) = ',pmaxn,  pmaxs
! Subprogram not used 
! Subprogram not used         if (print_global) then  ! global diags for conservations checks
! Subprogram not used 
! Subprogram not used          write(nu_diag,*) '----------------------------'
! Subprogram not used          write(nu_diag,901) 'arwt rain h2o kg in dt = ',rnn,rns
! Subprogram not used          write(nu_diag,901) 'arwt snow h2o kg in dt = ',snn,sns
! Subprogram not used          write(nu_diag,901) 'arwt evap h2o kg in dt = ',evpn,evps
! Subprogram not used          write(nu_diag,901) 'arwt frzl h2o kg in dt = ',frzn,frzs
! Subprogram not used          write(nu_diag,901) 'arwt frsh h2o kg in dt = ',sfreshn,sfreshs
! Subprogram not used 
! Subprogram not used          write(nu_diag,901) 'arwt ice mass (kg)     = ',micen,mices
! Subprogram not used          write(nu_diag,901) 'arwt snw mass (kg)     = ',msnwn,msnws
! Subprogram not used  
! Subprogram not used          write(nu_diag,901) 'arwt tot mass (kg)     = ',mtotn,mtots
! Subprogram not used          write(nu_diag,901) 'arwt tot mass chng(kg) = ',delmin,delmis
! Subprogram not used          write(nu_diag,901) 'arwt water flux        = ',fluxn,fluxs
! Subprogram not used          if (update_ocn_f) then
! Subprogram not used            write (nu_diag,*) '(=rain+snow+evap-fresh)  '
! Subprogram not used          else
! Subprogram not used            write (nu_diag,*) '(=rain+snow+evap+frzl-fresh)  '
! Subprogram not used          endif
! Subprogram not used          write(nu_diag,901) 'water flux error       = ',werrn,werrs
! Subprogram not used 
! Subprogram not used          write(nu_diag,*) '----------------------------'
! Subprogram not used          write(nu_diag,901) 'arwt atm heat flux (W) = ',fhatmn,fhatms
! Subprogram not used          write(nu_diag,901) 'arwt ocn heat flux (W) = ',fhocnn,fhocns
! Subprogram not used          write(nu_diag,901) 'arwt frzl heat flux(W) = ',fhfrzn,fhfrzs
! Subprogram not used          write(nu_diag,901) 'arwt tot energy    (J) = ',etotn,etots
! Subprogram not used          write(nu_diag,901) 'arwt net heat      (J) = ',hnetn,hnets
! Subprogram not used          write(nu_diag,901) 'arwt tot energy chng(J)= ',delein,deleis
! Subprogram not used          write(nu_diag,901) 'arwt heat error        = ',herrn,herrs
! Subprogram not used 
! Subprogram not used          write(nu_diag,*) '----------------------------'
! Subprogram not used          write(nu_diag,901) 'arwt salt mass (kg)    = ',msltn,mslts
! Subprogram not used          write(nu_diag,901) 'arwt salt mass chng(kg)= ',delmsltn, &
! Subprogram not used                                                         delmslts
! Subprogram not used          write(nu_diag,901) 'arwt salt flx in dt(kg)= ',sfsaltn, &
! Subprogram not used                                                         sfsalts
! Subprogram not used          write(nu_diag,901) 'arwt salt flx error    = ',serrn,serrs
! Subprogram not used          write(nu_diag,*) '----------------------------'
! Subprogram not used 
! Subprogram not used         endif                    ! print_global
! Subprogram not used        endif                     ! grid_type
! Subprogram not used 
! Subprogram not used        call flush_fileunit(nu_diag)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! diagnostics for Arctic and Antarctic points
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used        if (print_points) then
! Subprogram not used 
! Subprogram not used         write(nu_diag,*) '                         '
! Subprogram not used         write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1), &
! Subprogram not used                                                        plat(2),plon(2)
! Subprogram not used         write(nu_diag,903) '  my_task, iblk, i, j     ', &
! Subprogram not used                               pmloc(1),pbloc(1),piloc(1),pjloc(1), &
! Subprogram not used                               pmloc(2),pbloc(2),piloc(2),pjloc(2)
! Subprogram not used         write(nu_diag,*) '----------atm----------'
! Subprogram not used         write(nu_diag,900) 'air temperature (C)    = ',pTair(1),pTair(2)
! Subprogram not used         write(nu_diag,900) 'specific humidity      = ',pQa(1),pQa(2)
! Subprogram not used         write(nu_diag,900) 'snowfall (m)           = ',pfsnow(1), &
! Subprogram not used                                                        pfsnow(2)
! Subprogram not used         write(nu_diag,900) 'rainfall (m)           = ',pfrain(1), &
! Subprogram not used                                                        pfrain(2)
! Subprogram not used         if (.not.calc_Tsfc) then
! Subprogram not used            write(nu_diag,900) 'total surface heat flux= ',pfsurf(1),pfsurf(2)
! Subprogram not used            write(nu_diag,900) 'top sfc conductive flux= ',pfcondtop(1), &
! Subprogram not used                                                           pfcondtop(2)
! Subprogram not used            write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
! Subprogram not used         else
! Subprogram not used            write(nu_diag,900) 'shortwave radiation sum= ',pfsw(1),pfsw(2)
! Subprogram not used            write(nu_diag,900) 'longwave radiation     = ',pflw(1),pflw(2)
! Subprogram not used         endif
! Subprogram not used         write(nu_diag,*) '----------ice----------'
! Subprogram not used         write(nu_diag,900) 'area fraction          = ',paice(1),paice(2)
! Subprogram not used         write(nu_diag,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
! Subprogram not used         write(nu_diag,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
! Subprogram not used         if (calc_Tsfc) then
! Subprogram not used            write(nu_diag,900) 'surface temperature(C) = ',pTsfc(1),pTsfc(2)
! Subprogram not used            write(nu_diag,900) 'absorbed shortwave flx = ',pfswabs(1), &
! Subprogram not used                                                           pfswabs(2)
! Subprogram not used            write(nu_diag,900) 'outward longwave flx   = ',pflwout(1), &
! Subprogram not used                                                           pflwout(2)
! Subprogram not used            write(nu_diag,900) 'sensible heat flx      = ',pfsens(1), &
! Subprogram not used                                                           pfsens(2)
! Subprogram not used            write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
! Subprogram not used         endif
! Subprogram not used         write(nu_diag,900) 'subl/cond (m ice)      = ',pevap(1),pevap(2)
! Subprogram not used         write(nu_diag,900) 'top melt (m)           = ',pmeltt(1) &
! Subprogram not used                                                       ,pmeltt(2)
! Subprogram not used         write(nu_diag,900) 'bottom melt (m)        = ',pmeltb(1) &
! Subprogram not used                                                       ,pmeltb(2)
! Subprogram not used         write(nu_diag,900) 'lateral melt (m)       = ',pmeltl(1) &
! Subprogram not used                                                       ,pmeltl(2)
! Subprogram not used         write(nu_diag,900) 'new ice (m)            = ',pfrazil(1), &
! Subprogram not used                                                        pfrazil(2)
! Subprogram not used         write(nu_diag,900) 'congelation (m)        = ',pcongel(1), &
! Subprogram not used                                                        pcongel(2)
! Subprogram not used         write(nu_diag,900) 'snow-ice (m)           = ',psnoice(1), &
! Subprogram not used                                                        psnoice(2)
! Subprogram not used         write(nu_diag,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
! Subprogram not used         write(nu_diag,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
! Subprogram not used         write(nu_diag,900) 'intnl enrgy chng(W/m^2)= ',pde (1),pde (2)
! Subprogram not used         write(nu_diag,*) '----------ocn----------'
! Subprogram not used         write(nu_diag,900) 'sst (C)                = ',psst(1),psst(2)
! Subprogram not used         write(nu_diag,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
! Subprogram not used         write(nu_diag,900) 'heat used (W/m^2)      = ',pfhocn(1), &
! Subprogram not used                                                        pfhocn(2)
! Subprogram not used 
! Subprogram not used        endif                    ! print_points
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used   799 format (27x,a24)
! Subprogram not used   800 format (a25,2x,f24.17)
! Subprogram not used   801 format (a25,2x,1pe24.17)
! Subprogram not used   899 format (27x,a24,2x,a24)
! Subprogram not used   900 format (a25,2x,f24.17,2x,f24.17)
! Subprogram not used   901 format (a25,2x,1pe24.17,2x,1pe24.17)
! Subprogram not used   902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
! Subprogram not used   903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)
! Subprogram not used 
! Subprogram not used       end subroutine runtime_diags

!=======================================================================
!BOP
!
! !IROUTINE: init_mass_diags - computes global combined ice and snow mass sum
!
! !INTERFACE:
!
      subroutine init_mass_diags
!
! !DESCRIPTION:
!
! Computes global combined ice and snow mass sum
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_global_reductions
      use ice_grid
      use ice_state
      use ice_broadcast
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: n, k, ii, jj, i, j, iblk
      integer (kind=int_kind) :: naero

      real (kind=dbl_kind) :: &
         shmaxn, snwmxn,  shmaxs, snwmxs

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         work1, work2


      ! total ice volume
      shmaxn = global_sum(vice, distrb_info, field_loc_center, tarean)
      shmaxs = global_sum(vice, distrb_info, field_loc_center, tareas)

      ! total snow volume
      snwmxn = global_sum(vsno, distrb_info, field_loc_center, tarean)
      snwmxs = global_sum(vsno, distrb_info, field_loc_center, tareas)

      ! north/south ice mass
      totmin = rhoi*shmaxn
      totmis = rhoi*shmaxs

      ! north/south ice+snow mass
      totmn = totmin + rhos*snwmxn
      totms = totmis + rhos*snwmxs

      ! north/south ice+snow energy
      ! total ice/snow energy
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
      do j=1,ny_block
      do i=1,nx_block
         work1(i,j,iblk) = esno(i,j,iblk) + eice(i,j,iblk)
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO
      
      toten = global_sum(work1, distrb_info, field_loc_center, tarean)
      totes = global_sum(work1, distrb_info, field_loc_center, tareas)

      if (tr_aero) then
       do naero=1,n_aero
        do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = trcr(i,j,nt_aero  +4*(naero-1),iblk)*vsno(i,j,iblk) &
                            + trcr(i,j,nt_aero+1+4*(naero-1),iblk)*vsno(i,j,iblk) &
                            + trcr(i,j,nt_aero+2+4*(naero-1),iblk)*vice(i,j,iblk) &
                            + trcr(i,j,nt_aero+3+4*(naero-1),iblk)*vice(i,j,iblk)
         enddo
         enddo
        enddo
        totaeron(naero)= global_sum(work1, distrb_info, field_loc_center, tarean)
        totaeros(naero)= global_sum(work1, distrb_info, field_loc_center, tareas)
       enddo
      endif

      if (print_points) then

         do n = 1, npnt

            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               pdhi(n) = vice(i,j,iblk)
               pdhs(n) = vsno(i,j,iblk)
               pde(n)  = esno(i,j,iblk) + eice(i,j,iblk)
            endif

         enddo  ! npnt

      endif                     ! print_points

      end subroutine init_mass_diags

!=======================================================================
!BOP
!
! !IROUTINE: init_diags - find tasks for diagnostic points
!
! !INTERFACE:
!
      subroutine init_diags
!
! !DESCRIPTION:
!
!  Find tasks for diagnostic points.
!
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
      use ice_grid
      use ice_blocks
      use ice_broadcast
      use ice_global_reductions
      use ice_gather_scatter
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind) :: &
         latdis  , & ! latitude distance
         londis  , & ! longitude distance
         totdis  , & ! total distance
         mindis  , & ! minimum distance from desired location
         mindis_g    ! global minimum distance from desired location

      integer (kind=int_kind) :: &
         n           , & ! index for point search
         i,j         , & ! grid indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      character (char_len) :: label(npnt)

      type (block) :: &
         this_block           ! block information for current block

      if (print_points) then

         if (my_task==master_task) then
            write(nu_diag,*) ' '
            write(nu_diag,*) ' Find indices of diagnostic points '
         endif

         ! initialize labels
         label(1)(1:40)  = 'Near North Pole pack ice                '
         label(2)(1:40)  = 'Weddell Sea                             '

         piloc(:) = 0
         pjloc(:) = 0
         pbloc(:) = 0
         pmloc(:) = -999
         plat(:)  = -999._dbl_kind
         plon(:)  = -999._dbl_kind

         ! find minimum distance to diagnostic points on this processor 
         do n = 1, npnt
            if (lonpnt(n) > c180) lonpnt(n) = lonpnt(n) - c360

            iindx = 0
            jindx = 0
            bindx = 0
            mindis = 540.0_dbl_kind !  360. + 180.

            !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,j,i, &
            !$OMP                     latdis,londis,totdis,mindis, &
            !$OMP                     jindx,iindx,bindx) 
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               do j = jlo, jhi
               do i = ilo, ihi
                  if (hm(i,j,iblk) > p5) then
                     latdis = abs(latpnt(n)-TLAT(i,j,iblk)*rad_to_deg)
                     londis = abs(lonpnt(n)-TLON(i,j,iblk)*rad_to_deg) &
                            * cos(TLAT(i,j,iblk))
                     totdis = sqrt(latdis**2 + londis**2)
                     if (totdis < mindis) then
                        mindis = totdis
                        jindx = j
                        iindx = i
                        bindx = iblk
                     endif      ! totdis < mindis
                  endif         ! hm > p5
               enddo            ! i
               enddo            ! j
            enddo               ! iblk
            !$OMP END PARALLEL DO 

            ! find global minimum distance to diagnostic points 
            mindis_g = global_minval(mindis, distrb_info)

            ! save indices of minimum-distance grid cell
            if (abs(mindis_g - mindis) < puny) then
               piloc(n) = iindx
               pjloc(n) = jindx
               pbloc(n) = bindx
               pmloc(n) = my_task
               plat(n) = TLAT(iindx,jindx,bindx)*rad_to_deg
               plon(n) = TLON(iindx,jindx,bindx)*rad_to_deg
            endif

            ! communicate to all processors
            piloc(n) = global_maxval(piloc(n), distrb_info)
            pjloc(n) = global_maxval(pjloc(n), distrb_info)
            pbloc(n) = global_maxval(pbloc(n), distrb_info)
            pmloc(n) = global_maxval(pmloc(n), distrb_info)
            plat(n)  = global_maxval(plat(n), distrb_info)
            plon(n)  = global_maxval(plon(n), distrb_info)

            ! write to log file
            if (my_task==master_task) then
               write(nu_diag,*) ' '
               write(nu_diag,100) n,latpnt(n),lonpnt(n),plat(n),plon(n), &
                    piloc(n), pjloc(n), pbloc(n), pmloc(n)
            endif
 100        format(' found point',i4/ &
               '   lat    lon   TLAT   TLON     i     j   block  task'/ &
                4(f6.1,1x),1x,4(i4,2x) )

         enddo                  ! npnt
      endif                     ! print_points

      end subroutine init_diags

!=======================================================================
!BOP
!
! !IROUTINE: print_state - print ice state for specified grid point
!
! !INTERFACE:
!
! Subprogram not used       subroutine print_state(plabel,i,j,iblk)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! This routine is useful for debugging.
! Subprogram not used ! Calls to it should be inserted in the form (after thermo, for example)
! Subprogram not used !      do iblk = 1, nblocks
! Subprogram not used !      do j=jlo,jhi
! Subprogram not used !      do i=ilo,ihi
! Subprogram not used !         plabel = 'post thermo'
! Subprogram not used !         if (istep1 >= check_step .and. iblk==iblkp .and i==ip &
! Subprogram not used !             .and. j==jp .and. my_task == mtask) &
! Subprogram not used !         call print_state(plabel,i,j,iblk)
! Subprogram not used !      enddo
! Subprogram not used !      enddo
! Subprogram not used !      enddo
! Subprogram not used !
! Subprogram not used ! 'use ice_diagnostics' may need to be inserted also
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !MH      use ice_domain_size
! Subprogram not used       use ice_state
! Subprogram not used       use ice_itd
! Subprogram not used       use ice_flux
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character (len=20), intent(in) :: plabel
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in) :: & 
! Subprogram not used           i, j       , & ! horizontal indices
! Subprogram not used           iblk           ! block index
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used            eidebug, esdebug, &
! Subprogram not used            qi, qs, Tsnow
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: n, k
! Subprogram not used 
! Subprogram not used       write(nu_diag,*) plabel
! Subprogram not used       write(nu_diag,*) 'istep1, my_task, i, j, iblk:', &
! Subprogram not used                         istep1, my_task, i, j, iblk
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used       write(nu_diag,*) 'aice0', aice0(i,j,iblk)
! Subprogram not used       do n = 1, ncat
! Subprogram not used          write(nu_diag,*) ' '
! Subprogram not used          write(nu_diag,*) 'n =',n
! Subprogram not used          write(nu_diag,*) 'aicen', aicen(i,j,n,iblk)
! Subprogram not used          write(nu_diag,*) 'vicen', vicen(i,j,n,iblk)
! Subprogram not used          write(nu_diag,*) 'vsnon', vsnon(i,j,n,iblk)
! Subprogram not used          if (aicen(i,j,n,iblk) > puny) then
! Subprogram not used             write(nu_diag,*) 'hin', vicen(i,j,n,iblk)/aicen(i,j,n,iblk)
! Subprogram not used             write(nu_diag,*) 'hsn', vsnon(i,j,n,iblk)/aicen(i,j,n,iblk)
! Subprogram not used          endif
! Subprogram not used          write(nu_diag,*) 'Tsfcn',trcrn(i,j,nt_Tsfc,n,iblk)
! Subprogram not used          write(nu_diag,*) ' '
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       eidebug = c0
! Subprogram not used       do n = 1,ncat
! Subprogram not used          do k = 1,nilyr
! Subprogram not used             write(nu_diag,*) 'eicen, cat ',n,' layer ',k, &
! Subprogram not used                  eicen(i,j,ilyr1(n)+k-1,iblk)
! Subprogram not used             eidebug = eidebug + eicen(i,j,ilyr1(n)+k-1,iblk)
! Subprogram not used             if (aicen(i,j,n,iblk) > puny) then
! Subprogram not used                qi = eicen(i,j,ilyr1(n)+k-1,iblk) / & ! qi, eicen < 0
! Subprogram not used                    (vicen(i,j,n,iblk)/real(nilyr,kind=dbl_kind))
! Subprogram not used                write(nu_diag,*)  'qi/rhoi', qi/rhoi
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          write(nu_diag,*) ' '
! Subprogram not used       enddo
! Subprogram not used       write(nu_diag,*) 'eice(i,j)',eidebug
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used 
! Subprogram not used       esdebug = c0
! Subprogram not used       do n = 1,ncat
! Subprogram not used          if (vsnon(i,j,n,iblk) > puny) then
! Subprogram not used             do k = 1,nslyr
! Subprogram not used                write(nu_diag,*) 'esnon, cat ',n,' layer ',k, &
! Subprogram not used                   esnon(i,j,slyr1(n)+k-1,iblk)
! Subprogram not used                esdebug = esdebug + esnon(i,j,slyr1(n)+k-1,iblk)
! Subprogram not used                qs = esnon(i,j,slyr1(n)+k-1,iblk) / &  ! qs, esnon < 0
! Subprogram not used                    (vsnon(i,j,n,iblk)/real(nslyr,kind=dbl_kind))
! Subprogram not used                Tsnow = (Lfresh + qs/rhos) / cp_ice
! Subprogram not used                write(nu_diag,*) 'qs/rhos', qs/rhos
! Subprogram not used                write(nu_diag,*) 'Tsnow', Tsnow
! Subprogram not used             enddo
! Subprogram not used             write(nu_diag,*) ' '
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       write(nu_diag,*) 'esno(i,j)',esdebug
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used 
! Subprogram not used       write(nu_diag,*) 'uvel(i,j)',uvel(i,j,iblk)
! Subprogram not used       write(nu_diag,*) 'vvel(i,j)',vvel(i,j,iblk)
! Subprogram not used 
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used       write(nu_diag,*) 'atm states and fluxes'
! Subprogram not used       write(nu_diag,*) '            uatm    = ',uatm (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            vatm    = ',vatm (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            potT    = ',potT (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            Tair    = ',Tair (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            Qa      = ',Qa   (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            rhoa    = ',rhoa (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            swvdr   = ',swvdr(i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            swvdf   = ',swvdf(i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            swidr   = ',swidr(i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            swidf   = ',swidf(i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            flw     = ',flw  (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            frain   = ',frain(i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            fsnow   = ',fsnow(i,j,iblk)
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used       write(nu_diag,*) 'ocn states and fluxes'
! Subprogram not used       write(nu_diag,*) '            frzmlt  = ',frzmlt (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            sst     = ',sst    (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            sss     = ',sss    (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            Tf      = ',Tf     (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            uocn    = ',uocn   (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            vocn    = ',vocn   (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            strtltx = ',strtltx(i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            strtlty = ',strtlty(i,j,iblk)
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used       write(nu_diag,*) 'srf states and fluxes'
! Subprogram not used       write(nu_diag,*) '            Tref    = ',Tref  (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            Qref    = ',Qref  (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            fsens   = ',fsens (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            flat    = ',flat  (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            evap    = ',evap  (i,j,iblk)
! Subprogram not used       write(nu_diag,*) '            flwout  = ',flwout(i,j,iblk)
! Subprogram not used       write(nu_diag,*) ' '
! Subprogram not used 
! Subprogram not used       end subroutine print_state

!=======================================================================

      end module ice_diagnostics

!=======================================================================





