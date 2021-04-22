!=======================================================================
!
!BOP
!
! !MODULE: ice_step_mod
!
! !DESCRIPTION:
!
!  Contains CICE component driver routines common to all drivers.
!
! !REVISION HISTORY:
!  SVN:$Id: $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2008 ECH: created module by moving subroutines from drivers/cice4/
!
! !INTERFACE:
!
      module ice_step_mod
!
! !USES:
!
      use ice_atmo
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_fileunits
      use ice_flux
      use ice_grid
      use ice_history
      use ice_restart
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_state
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use perf_mod, only: t_startf, t_stopf, t_barrierf

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: step_therm2, step_dynamics, &
                prep_radiation, step_radiation
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: prep_radiation - step pre-thermo radiation
!
! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
! authors: Mariana Vertenstein, NCAR
!
! !INTERFACE:

      subroutine prep_radiation(dt)
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
      integer (kind=int_kind) :: &
         i,j,n,iblk    ! block index

      if (calc_Tsfc) then

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call prep_radiation_iblk(dt, iblk)
         end do
         !$OMP END PARALLEL DO

      else    ! .not. calc_Tsfc

         ! Initialize for safety
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            fswsfcn(i,j,n,iblk) = c0
            fswintn(i,j,n,iblk) = c0
            fswthrun(i,j,n,iblk) = c0
         enddo   ! i
         enddo   ! j
         enddo   ! ncat
            Iswabsn(:,:,:,iblk) = c0
            Sswabsn(:,:,:,iblk) = c0
         enddo   ! iblk

      endif    ! calc_Tsfc

      end subroutine prep_radiation

!=======================================================================
!BOP
!
! !ROUTINE: prep_radiation_iblk - step pre-thermo radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine prep_radiation_iblk (dt, iblk)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micro-meters)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      real (kind=dbl_kind) :: netsw, netsw_old, ar

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      l_stop = .false.

      fswfac(:,:,iblk) = c1

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j,iblk) > c0 .and. scale_factor(i,j,iblk) > puny) then
            netsw = swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                  + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                  + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                  + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))
            scale_factor(i,j,iblk) = netsw / scale_factor(i,j,iblk)
         else
            scale_factor(i,j,iblk) = c1
         endif
         fswfac(i,j,iblk) = scale_factor(i,j,iblk) ! for history
      enddo               ! i
      enddo               ! j

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

      !-----------------------------------------------------------------
      ! Scale absorbed solar radiation for change in net shortwave
      !-----------------------------------------------------------------

         il1 = ilyr1(n)
         il2 = ilyrn(n)
         sl1 = slyr1(n)
         sl2 = slyrn(n)

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            fswsfcn(i,j,n,iblk) = scale_factor(i,j,iblk)*fswsfcn (i,j,n,iblk)
            fswintn(i,j,n,iblk) = scale_factor(i,j,iblk)*fswintn (i,j,n,iblk)
            fswthrun(i,j,n,iblk)= scale_factor(i,j,iblk)*fswthrun(i,j,n,iblk)
            Sswabsn(i,j,sl1:sl2,iblk) = &
                    scale_factor(i,j,iblk)*Sswabsn(i,j,sl1:sl2,iblk)
            Iswabsn(i,j,il1:il2,iblk) = &
                    scale_factor(i,j,iblk)*Iswabsn(i,j,il1:il2,iblk)
         enddo
      enddo                  ! ncat

      end subroutine prep_radiation_iblk

!=======================================================================
!BOP
!
! !ROUTINE: step_therm2 - step post-coupler thermodynamics
!
! !DESCRIPTION:
!
!-----------------------------------------------------------------------
! Wrapper for driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting. Needed for 
! introducing OpenMP threading more simply.
!
! !REVISION HISTORY:
!
! author: Mariana Vertenstein, NCAR
!
! !INTERFACE:
!
! Subprogram not used       subroutine step_therm2 (dt)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          dt      ! time step
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          iblk,  &  ! block index
! Subprogram not used          i, j
! Subprogram not used 
! Subprogram not used !      call t_barrierf('cice_step2_therm_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf('cice_step2_therm')
! Subprogram not used       call ice_timer_start(timer_column)  ! column physics
! Subprogram not used       call ice_timer_start(timer_thermo)  ! thermodynamics
! Subprogram not used !      call ice_timer_start(timer_tmp)  ! temporary timer
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          call step_therm2_iblk(dt, iblk)
! Subprogram not used       end do
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       !-------------------------------------------------------------------
! Subprogram not used       ! Ghost cell updates for state variables.
! Subprogram not used       !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used       call bound_state (aicen, trcrn, &
! Subprogram not used                         vicen, vsnon, &
! Subprogram not used                         eicen, esnon)
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Aggregate the updated state variables. 
! Subprogram not used       !----------------------------------------------------------------- 
! Subprogram not used  
! Subprogram not used          call aggregate (nx_block,          ny_block,             &
! Subprogram not used                          aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
! Subprogram not used                          vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
! Subprogram not used                          eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
! Subprogram not used                          aice (:,:,  iblk), trcr (:,:,:,  iblk),  &
! Subprogram not used                          vice (:,:,  iblk), vsno (:,:,    iblk),  &
! Subprogram not used                          eice (:,:,  iblk), esno (:,:,    iblk),  &
! Subprogram not used                          aice0(:,:,  iblk), tmask(:,:,    iblk),  &
! Subprogram not used                          ntrcr,             trcr_depend) 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Compute thermodynamic area and volume tendencies.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             daidtt(i,j,iblk) = (aice(i,j,iblk) - daidtt(i,j,iblk)) / dt
! Subprogram not used             dvidtt(i,j,iblk) = (vice(i,j,iblk) - dvidtt(i,j,iblk)) / dt
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf('cice_step2_therm')
! Subprogram not used !      call ice_timer_stop(timer_tmp)  ! temporary timer
! Subprogram not used       call ice_timer_stop(timer_thermo)  ! column physics
! Subprogram not used       call ice_timer_stop(timer_column)  ! column physics
! Subprogram not used 
! Subprogram not used       end subroutine step_therm2

!=======================================================================
!BOP
!
! !ROUTINE: step_therm2_iblk - step post-coupler thermodynamics
!
! !DESCRIPTION:
!
!-----------------------------------------------------------------------
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! NOTE: Ocean fluxes are initialized here.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:

! Subprogram not used       subroutine step_therm2_iblk (dt, iblk)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          dt      ! time step
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          iblk ! block index
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used !lipscomb - delete hicen later?
! Subprogram not used !      real (kind=dbl_kind), &
! Subprogram not used !         dimension (nx_block,ny_block,ncat,max_blocks) :: &
! Subprogram not used !         hicen           ! ice thickness (m)
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
! Subprogram not used          i, j, n
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          icells          ! number of ice/ocean cells 
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block) :: &
! Subprogram not used          indxi, indxj    ! indirect indices for ice/ocean cells
! Subprogram not used 
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block      ! block information for current block
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used          l_stop          ! if true, abort model
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          istop, jstop    ! indices of grid cell where model aborts
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used 
! Subprogram not used       this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used       ilo = this_block%ilo
! Subprogram not used       ihi = this_block%ihi
! Subprogram not used       jlo = this_block%jlo
! Subprogram not used       jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Let rain drain through to the ocean.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          fresh     (i,j,iblk) = fresh(i,j,iblk)       &
! Subprogram not used               + frain(i,j,iblk)*aice(i,j,iblk)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Given thermodynamic growth rates, transport ice between
! Subprogram not used       ! thickness categories.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_catconv, iblk) ! category conversions
! Subprogram not used 
! Subprogram not used       if (kitd == 1) then
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Compute fractional ice area in each grid cell.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          call aggregate_area (nx_block,          ny_block, &
! Subprogram not used                               aicen(:,:,:,iblk),           &
! Subprogram not used                               aice (:,:,  iblk), aice0(:,:,iblk))
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Identify grid cells with ice.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          icells = 0
! Subprogram not used          do j = jlo,jhi
! Subprogram not used          do i = ilo,ihi
! Subprogram not used             if (aice(i,j,iblk) > puny) then
! Subprogram not used                icells = icells + 1
! Subprogram not used                indxi(icells) = i
! Subprogram not used                indxj(icells) = j
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          if (icells > 0) then
! Subprogram not used 
! Subprogram not used             call linear_itd (nx_block, ny_block,       &
! Subprogram not used                              icells, indxi, indxj,     &
! Subprogram not used                              ntrcr,  trcr_depend,      &
! Subprogram not used                              aicen_init(:,:,:,iblk),   &
! Subprogram not used                              vicen_init(:,:,:,iblk),   &
! Subprogram not used                              aicen     (:,:,:,iblk),   &
! Subprogram not used                              trcrn     (:,:,:,:,iblk), & 
! Subprogram not used                              vicen     (:,:,:,iblk),   &
! Subprogram not used                              vsnon     (:,:,:,iblk),   &
! Subprogram not used                              eicen     (:,:,:,iblk),   &
! Subprogram not used                              esnon     (:,:,:,iblk),   &
! Subprogram not used                              aice      (:,:,  iblk),   &
! Subprogram not used                              aice0     (:,:,  iblk),   &
! Subprogram not used                              l_stop,                   &
! Subprogram not used                              istop,    jstop)
! Subprogram not used 
! Subprogram not used             if (l_stop) then
! Subprogram not used                write (nu_diag,*) 'istep1, my_task, iblk =', &
! Subprogram not used                                   istep1, my_task, iblk
! Subprogram not used                write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used                if (istop > 0 .and. jstop > 0) &
! Subprogram not used                     write(nu_diag,*) 'Global i and j:', &
! Subprogram not used                                       this_block%i_glob(istop), &
! Subprogram not used                                       this_block%j_glob(jstop) 
! Subprogram not used                call abort_ice ('ice: Linear ITD error')
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       endif  ! kitd
! Subprogram not used 
! Subprogram not used       call ice_timer_stop(timer_catconv, iblk)  ! category conversions
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Add frazil ice growing in leads.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       ! identify ice-ocean cells
! Subprogram not used       icells = 0
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          if (tmask(i,j,iblk)) then
! Subprogram not used             icells = icells + 1
! Subprogram not used             indxi(icells) = i
! Subprogram not used             indxj(icells) = j
! Subprogram not used          endif
! Subprogram not used       enddo               ! i
! Subprogram not used       enddo               ! j
! Subprogram not used 
! Subprogram not used       call add_new_ice (nx_block,              ny_block, &
! Subprogram not used                         ntrcr,                 icells,   &
! Subprogram not used                         indxi,                 indxj,    &
! Subprogram not used                         tmask    (:,:,  iblk), dt,       &
! Subprogram not used                         aicen    (:,:,:,iblk),           &
! Subprogram not used                         trcrn    (:,:,:,:,iblk),         &
! Subprogram not used                         vicen    (:,:,:,iblk),           &
! Subprogram not used                         eicen    (:,:,:,iblk),           &
! Subprogram not used                         aice0    (:,:,  iblk),           &
! Subprogram not used                         aice     (:,:,  iblk),           &
! Subprogram not used                         frzmlt   (:,:,  iblk),           &
! Subprogram not used                         frazil   (:,:,  iblk),           &
! Subprogram not used                         frz_onset(:,:,  iblk), yday,     &
! Subprogram not used                         fresh    (:,:,  iblk),           &
! Subprogram not used                         fsalt    (:,:,  iblk),           &
! Subprogram not used                         Tf       (:,:,  iblk), l_stop,   &
! Subprogram not used                         istop, jstop)
! Subprogram not used 
! Subprogram not used       if (l_stop) then
! Subprogram not used          write (nu_diag,*) 'istep1, my_task, iblk =', &
! Subprogram not used                             istep1, my_task, iblk
! Subprogram not used          write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used          if (istop > 0 .and. jstop > 0) &
! Subprogram not used               write(nu_diag,*) 'Global i and j:', &
! Subprogram not used                                 this_block%i_glob(istop), &
! Subprogram not used                                 this_block%j_glob(jstop) 
! Subprogram not used          call abort_ice ('ice: add_new_ice error')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Melt ice laterally.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       call lateral_melt (nx_block, ny_block,     &
! Subprogram not used                          ilo, ihi, jlo, jhi,     &
! Subprogram not used                          dt,                     &
! Subprogram not used                          fresh     (:,:,  iblk), &
! Subprogram not used                          fsalt     (:,:,  iblk), &    
! Subprogram not used                          fhocn     (:,:,  iblk), &
! Subprogram not used                          fsoot     (:,:,:,iblk), &
! Subprogram not used                          rside     (:,:,  iblk), &
! Subprogram not used                          meltl     (:,:,  iblk), &
! Subprogram not used                          aicen     (:,:,:,iblk), &
! Subprogram not used                          vicen     (:,:,:,iblk), &
! Subprogram not used                          vsnon     (:,:,:,iblk), &
! Subprogram not used                          eicen     (:,:,:,iblk), &
! Subprogram not used                          esnon     (:,:,:,iblk), &
! Subprogram not used                          trcrn     (:,:,:,:,iblk) )
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! For the special case of a single category, adjust the area and
! Subprogram not used       ! volume (assuming that half the volume change decreases the
! Subprogram not used       ! thickness, and the other half decreases the area).  
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !NOTE - this does not work - hicen_init is not defined - ECH
! Subprogram not used 
! Subprogram not used !         if (ncat==1) &
! Subprogram not used !              call reduce_area (nx_block, ny_block,     &
! Subprogram not used !                                ilo, ihi, jlo, jhi,     &
! Subprogram not used !                                tmask     (:,:,  iblk), &
! Subprogram not used !                                aicen     (:,:,:,iblk), &
! Subprogram not used !                                vicen     (:,:,:,iblk), &
! Subprogram not used !                                hicen_init(:,:,1,iblk), &
! Subprogram not used !                                hicen     (:,:,1,iblk)) 
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! ITD cleanup: Rebin thickness categories if necessary, and remove
! Subprogram not used       !  categories with very small areas.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call cleanup_itd (nx_block,             ny_block,             &
! Subprogram not used                         ilo, ihi,             jlo, jhi,             &
! Subprogram not used                         dt,                   ntrcr,                &
! Subprogram not used                         aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
! Subprogram not used                         vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
! Subprogram not used                         eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
! Subprogram not used                         aice0   (:,:,  iblk), aice      (:,:,iblk), &
! Subprogram not used                         trcr_depend,                                &
! Subprogram not used                         fresh   (:,:,  iblk), fsalt   (:,:,  iblk), &
! Subprogram not used                         fhocn   (:,:,  iblk), fsoot   (:,:,:,iblk), &
! Subprogram not used                         tr_aero,                                    &
! Subprogram not used                         heat_capacity,        l_stop,               &
! Subprogram not used                         istop,                jstop)
! Subprogram not used 
! Subprogram not used       if (l_stop) then
! Subprogram not used          write (nu_diag,*) 'istep1, my_task, iblk =', &
! Subprogram not used                             istep1, my_task, iblk
! Subprogram not used          write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used          if (istop > 0 .and. jstop > 0) &
! Subprogram not used               write(nu_diag,*) 'Global i and j:', &
! Subprogram not used                                 this_block%i_glob(istop), &
! Subprogram not used                                 this_block%j_glob(jstop) 
! Subprogram not used          call abort_ice ('ice: ITD cleanup error')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       end subroutine step_therm2_iblk

!=======================================================================
!BOP
!
! !ROUTINE: step_dynamics - step ice dynamics, transport, and ridging
!
! !DESCRIPTION:
!
! Run one time step of dynamics, horizontal transport, and ridging.
! NOTE: The evp and transport modules include boundary updates, so
!       they cannot be done inside a single block loop.  Ridging
!       and cleanup, on the other hand, are single-column operations. 
!       They are called with argument lists inside block loops
!       to increase modularity.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !INTERFACE:

! Subprogram not used       subroutine step_dynamics(dt_dyn,dt_thm)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          dt_dyn  , & ! dynamic time step
! Subprogram not used          dt_thm      ! thermodynamic time step for diagnostics
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       type (block) :: &
! Subprogram not used          this_block      ! block information for current block
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: & 
! Subprogram not used          iblk        , & ! block index 
! Subprogram not used          i,j         , & ! horizontal indices
! Subprogram not used          ilo,ihi,jlo,jhi ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          icells          ! number of cells with aicen > puny
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block) :: &
! Subprogram not used          indxi, indxj    ! indirect indices for cells with aicen > puny
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used          l_stop          ! if true, abort model
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          istop, jstop    ! indices of grid cell where model aborts
! Subprogram not used 
! Subprogram not used       call init_history_dyn     ! initialize dynamic history variables
! Subprogram not used 
! Subprogram not used       dynCnt = dynCnt + 1
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Elastic-viscous-plastic ice dynamics
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      call t_barrierf ('cice_step_evp_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_step_evp')
! Subprogram not used       if (kdyn == 1) call evp (dt_dyn)
! Subprogram not used       call t_stopf ('cice_step_evp')
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Horizontal ice transport
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      call t_barrierf ('cice_step_horz_transport_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_step_horz_transport')
! Subprogram not used       if (advection == 'upwind') then
! Subprogram not used          call transport_upwind (dt_dyn)    ! upwind
! Subprogram not used       else
! Subprogram not used          call transport_remap (dt_dyn)     ! incremental remapping
! Subprogram not used       endif
! Subprogram not used       call t_stopf ('cice_step_horz_transport')
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Ridging
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call ice_timer_start(timer_column)
! Subprogram not used       call ice_timer_start(timer_ridge)
! Subprogram not used !      call t_barrierf ('cice_step_ridge_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_step_ridge')
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,&
! Subprogram not used       !$OMP	                icells,indxi,indxj,l_stop,istop,jstop)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk), iblk)
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Identify ice-ocean cells.
! Subprogram not used       ! Note:  We can not define icells here using aice>puny because
! Subprogram not used       !        aice has not yet been updated since the transport (and
! Subprogram not used       !        it may be out of whack, which the ridging helps fix).-ECH
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used            
! Subprogram not used          icells = 0
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used !PERF  There is a performance advantage here to only perform the ridging over a limited set of points
! Subprogram not used !PERF            if (tmask(i,j,iblk) .and. aice(i,j,iblk) > c0) then
! Subprogram not used             if (tmask(i,j,iblk)) then
! Subprogram not used                icells = icells + 1
! Subprogram not used                indxi(icells) = i
! Subprogram not used                indxj(icells) = j
! Subprogram not used             endif
! Subprogram not used          enddo               ! i
! Subprogram not used          enddo               ! j
! Subprogram not used 
! Subprogram not used          if (icells > 0) then
! Subprogram not used 
! Subprogram not used          call ridge_ice (nx_block,             ny_block,                 &
! Subprogram not used                          dt_dyn,               dt_thm,                   &
! Subprogram not used                          ntrcr,                icells,                   &
! Subprogram not used                          indxi,                indxj,                    &
! Subprogram not used !!                         Delt    (:,:,  iblk), divu      (:,:,  iblk), &
! Subprogram not used                          rdg_conv(:,:,  iblk), rdg_shear (:,:,  iblk),   &
! Subprogram not used                          aicen   (:,:,:,iblk), trcrn     (:,:,:,:,iblk), &
! Subprogram not used                          vicen   (:,:,:,iblk), vsnon     (:,:,:,iblk),   &
! Subprogram not used                          eicen   (:,:,:,iblk), esnon     (:,:,:,iblk),   &
! Subprogram not used                          aice0   (:,:,  iblk),                           &
! Subprogram not used                          trcr_depend,          l_stop,                   &
! Subprogram not used                          istop,                jstop,                    &   
! Subprogram not used                          dardg1dt(:,:,iblk),   dardg2dt  (:,:,iblk),     &
! Subprogram not used                          dvirdgdt(:,:,iblk),   opening   (:,:,iblk),     &
! Subprogram not used                          fresh   (:,:,iblk),   fhocn     (:,:,iblk),     &
! Subprogram not used                          fsoot   (:,:,:,iblk))
! Subprogram not used 
! Subprogram not used          if (l_stop) then
! Subprogram not used             write (nu_diag,*) 'istep1, my_task, iblk =', &
! Subprogram not used                                istep1, my_task, iblk
! Subprogram not used             write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used             if (istop > 0 .and. jstop > 0) &
! Subprogram not used                  write(nu_diag,*) 'Global i and j:', &
! Subprogram not used                                   this_block%i_glob(istop), &
! Subprogram not used                                   this_block%j_glob(jstop) 
! Subprogram not used             call abort_ice ('ice: Ridging error')
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call ice_timer_stop(timer_ridge)
! Subprogram not used       call t_stopf ('cice_step_ridge')
! Subprogram not used 
! Subprogram not used !      call t_barrierf ('cice_step_column_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_step_column')
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,&
! Subprogram not used       !$OMP	                icells,indxi,indxj,l_stop,istop,jstop)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used          this_block = get_block(blocks_ice(iblk), iblk)
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! ITD cleanup: Rebin thickness categories if necessary, and remove
! Subprogram not used       !  categories with very small areas.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call cleanup_itd (nx_block,             ny_block,             &
! Subprogram not used                            ilo, ihi,             jlo, jhi,             &
! Subprogram not used                            dt_thm,               ntrcr,                &
! Subprogram not used                            aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
! Subprogram not used                            vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
! Subprogram not used                            eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
! Subprogram not used                            aice0   (:,:,  iblk), aice      (:,:,iblk), &
! Subprogram not used                            trcr_depend,                                &
! Subprogram not used                            fresh   (:,:,  iblk), fsalt   (:,:,  iblk), &
! Subprogram not used                            fhocn   (:,:,  iblk), fsoot   (:,:,:,iblk), &
! Subprogram not used                            tr_aero,                                    &
! Subprogram not used                            heat_capacity,        l_stop,               &
! Subprogram not used                            istop,                jstop)
! Subprogram not used 
! Subprogram not used          if (l_stop) then
! Subprogram not used             write (nu_diag,*) 'istep1, my_task, iblk =', &
! Subprogram not used                                istep1, my_task, iblk
! Subprogram not used             write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used             if (istop > 0 .and. jstop > 0) &
! Subprogram not used                  write(nu_diag,*) 'Global i and j:', &
! Subprogram not used                                   this_block%i_glob(istop), &
! Subprogram not used                                   this_block%j_glob(jstop) 
! Subprogram not used             call abort_ice ('ice: ITD cleanup error')
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo              ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf ('cice_step_column')
! Subprogram not used 
! Subprogram not used       !-------------------------------------------------------------------
! Subprogram not used       ! Ghost cell updates for state variables.
! Subprogram not used       !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      call t_barrierf ('cice_step_bound_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_step_bound')
! Subprogram not used       call ice_timer_start(timer_bound)
! Subprogram not used       call bound_state (aicen, trcrn, &
! Subprogram not used                         vicen, vsnon, &
! Subprogram not used                         eicen, esnon)
! Subprogram not used       call ice_timer_stop(timer_bound)
! Subprogram not used       call t_stopf ('cice_step_bound')
! Subprogram not used 
! Subprogram not used !      call t_barrierf ('cice_step_agg_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_step_agg')
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Aggregate the updated state variables. 
! Subprogram not used       !----------------------------------------------------------------- 
! Subprogram not used  
! Subprogram not used          call aggregate (nx_block,          ny_block,             &
! Subprogram not used                          aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
! Subprogram not used                          vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
! Subprogram not used                          eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
! Subprogram not used                          aice (:,:,  iblk), trcr (:,:,:,  iblk),  &
! Subprogram not used                          vice (:,:,  iblk), vsno (:,:,    iblk),  &
! Subprogram not used                          eice (:,:,  iblk), esno (:,:,    iblk),  &
! Subprogram not used                          aice0(:,:,  iblk), tmask(:,:,    iblk),  &
! Subprogram not used                          ntrcr,             trcr_depend) 
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Compute dynamic area and volume tendencies.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used          do j = jlo,jhi
! Subprogram not used          do i = ilo,ihi
! Subprogram not used             dvidtd(i,j,iblk) = (vice(i,j,iblk) - dvidtd(i,j,iblk)) /dt_dyn
! Subprogram not used             daidtd(i,j,iblk) = (aice(i,j,iblk) - daidtd(i,j,iblk)) /dt_dyn
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       enddo              ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf ('cice_step_agg')
! Subprogram not used       call ice_timer_stop(timer_column)
! Subprogram not used 
! Subprogram not used       end subroutine step_dynamics

!=======================================================================
!BOP
!
! !ROUTINE: step_radiation - step pre-coupler radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: Mariana Vertenstein, NCAR
!
! !INTERFACE:

      subroutine step_radiation(dt)
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
      integer (kind=int_kind) :: &
         i,j,n,iblk    ! block index


      alvdr(:,:,:) = c0
      alvdf(:,:,:) = c0
      alidr(:,:,:) = c0
      alidf(:,:,:) = c0
      Sswabsn(:,:,:,:) = c0
!      do iblk = 1,nblocks
!         alvdr(:,:,iblk) = c0
!         alvdf(:,:,iblk) = c0
!         alidr(:,:,iblk) = c0
!         alidf(:,:,iblk) = c0
!         Sswabsn(:,:,:,iblk) = c0
!      enddo

      if (calc_Tsfc) then

         !$OMP PARALLEL DO PRIVATE(iblk)
!         call t_barrierf('cice_step_radiationib_BARRIER',MPI_COMM_ICE)
         do iblk = 1, nblocks
            call t_startf('cice_step_radiationib')
            call step_radiation_iblk(dt, iblk)
            call t_stopf('cice_step_radiationib')
         end do
         !$OMP END PARALLEL DO

      else    ! .not. calc_Tsfc

         ! Initialize for safety
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdrn(i,j,n,iblk) = c0
            alidrn(i,j,n,iblk) = c0
            alvdfn(i,j,n,iblk) = c0
            alidfn(i,j,n,iblk) = c0
            fswsfcn(i,j,n,iblk) = c0
            fswintn(i,j,n,iblk) = c0
            fswthrun(i,j,n,iblk) = c0
         enddo   ! i
         enddo   ! j
         enddo   ! ncat
            Iswabsn(:,:,:,iblk) = c0
            Sswabsn(:,:,:,iblk) = c0
         enddo   ! iblk

      endif    ! calc_Tsfc

      end subroutine step_radiation

!=======================================================================
!BOP
!
! !ROUTINE: step_radiation_iblk - step pre-coupler radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine step_radiation_iblk (dt, iblk)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micro-meters)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      l_stop = .false.

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute cosine of solar zenith angle.
      ! This is used by the delta-Eddington shortwave module.
      ! Albedos are aggregated in merge_fluxes only for cells w/ coszen > 0.
      ! For basic shortwave, simply set coszen to a constant between 0 and 1.
      !-----------------------------------------------------------------

      if (trim(shortwave) == 'dEdd') then ! delta Eddington

!         call t_barrierf('cice_step_rib_swdedd1_BARRIER',MPI_COMM_ICE)
         call t_startf('cice_step_rib_swdedd1')
         ! identify ice-ocean cells
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         call compute_coszen (nx_block,         ny_block,       &
                              icells,                           &
                              indxi,            indxj,          &
                              tlat  (:,:,iblk), tlon(:,:,iblk), &
                              coszen(:,:,iblk), dt)

         call t_stopf('cice_step_rib_swdedd1')
      else                     ! basic (ccsm3) shortwave
         coszen(:,:,iblk) = p5 ! sun above the horizon
      endif

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
           
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

      !-----------------------------------------------------------------
      ! Solar radiation: albedo and absorbed shortwave
      !-----------------------------------------------------------------

         il1 = ilyr1(n)
         il2 = ilyrn(n)
         sl1 = slyr1(n)
         sl2 = slyrn(n)

         if (trim(shortwave) == 'dEdd') then   ! delta Eddington

      ! note that rhoswn, rsnw, fp, hp and Sswabs ARE NOT dimensioned with ncat
      ! BPB 19 Dec 2006

            ! set snow properties
!            call t_barrierf('cice_step_rib_swdedd2_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd2')
            call shortwave_dEdd_set_snow(nx_block, ny_block,           &
                              icells,                                  &
                              indxi,               indxj,              &
                              aicen(:,:,n,iblk),   vsnon(:,:,n,iblk),  &
                              trcrn(:,:,nt_Tsfc,n,iblk), fsn,          &
                              rhosnwn,             rsnwn)
            call t_stopf('cice_step_rib_swdedd2')


            if (.not. tr_pond) then

               ! set pond properties
!               call t_barrierf('cice_step_rib_swdedd3_BARRIER',MPI_COMM_ICE)
               call t_startf('cice_step_rib_swdedd3')
               call shortwave_dEdd_set_pond(nx_block, ny_block,            &
                                 icells,                                   &
                                 indxi,               indxj,               &
                                 aicen(:,:,n,iblk),                        &
                                 trcrn(:,:,nt_Tsfc,n,iblk),                &
                                 fsn,                 fpn,                 &
                                 hpn)
               call t_stopf('cice_step_rib_swdedd3')
            else

               fpn(:,:) = apondn(:,:,n,iblk)
               hpn(:,:) = hpondn(:,:,n,iblk)

            endif


!            call t_barrierf('cice_step_rib_swdedd7_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd7')
            call shortwave_dEdd(nx_block,        ny_block,            &
                                icells,                                 &
                                indxi,             indxj,               &
                                coszen(:,:, iblk),                      &
                                aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                vsnon(:,:,n,iblk), fsn,                 &
                                rhosnwn,           rsnwn,               &
                                fpn,               hpn,                 &
                                trcrn(:,:,:,n,iblk), tarea(:,:,iblk),   &
                                swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                alvdrn(:,:,n,iblk),alvdfn(:,:,n,iblk),  &
                                alidrn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                                fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                                fswthrun(:,:,n,iblk),                   &
                                Sswabsn(:,:,sl1:sl2,iblk),              &
                                Iswabsn(:,:,il1:il2,iblk),              &
                                albicen(:,:,n,iblk),                    &
                                albsnon(:,:,n,iblk),albpndn(:,:,n,iblk))
            call t_stopf('cice_step_rib_swdedd7')

         else

            Sswabsn(:,:,sl1:sl2,iblk) = c0

!            call t_barrierf('cice_step_rib_ccsm3_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_ccsm3')
            call shortwave_ccsm3(nx_block,     ny_block,           &
                           icells,                                 &
                           indxi,             indxj,               &
                           aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                           vsnon(:,:,n,iblk),                      &
                           trcrn(:,:,nt_Tsfc,n,iblk),              &
                           swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                           swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                           alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),  &
                           alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                           fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                           fswthrun(:,:,n,iblk),                   &
                           Iswabsn(:,:,il1:il2,iblk),              &
                           albicen(:,:,n,iblk),albsnon(:,:,n,iblk))

            call t_stopf('cice_step_rib_ccsm3')
         endif

      enddo                  ! ncat

      end subroutine step_radiation_iblk

!=======================================================================

      end module ice_step_mod

!=======================================================================
