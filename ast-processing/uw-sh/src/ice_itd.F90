!=======================================================================
!BOP
!
! !MODULE: ice_itd - initialize and redistribute ice in the ITD
!
! !DESCRIPTION:
!
! Routines to initialize the ice thickness distribution and
! utilities to redistribute ice among categories. These routines
! are not specific to a particular numerical implementation.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669--15,677.
!
! See Bitz, C.M., M.M. Holland, A.J. Weaver, M. Eby, 2001:
! Simulating the ice-thickness distribution in a climate model,
! J. Geophys. Res., 106, 2441--2464.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_itd.F90 138 2008-07-08 20:39:37Z eclare $
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
!
! 2004 WHL: Added multiple snow layers, block structure, cleanup_itd
! 2006 ECH: Added WMO standard ice thickness categories as kcatbound=2
!           Streamlined for efficiency 
!           Converted to free source form (F90)
!
! !INTERFACE:
!
      module ice_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_exit
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         kitd        , & ! type of itd conversions
                         !   0 = delta function
                         !   1 = linear remap
         kcatbound   , & !   0 = old category boundary formula
                         !   1 = new formula giving round numbers
                         !   2 = WMO standard
         ilyr1 (ncat), & ! array position of top ice layer in each cat
         ilyrn (ncat), & ! array position of bottom ice layer in each cat
         slyr1 (ncat), & ! array position of top snow layer in each cat
         slyrn (ncat)    ! array position of bottom snow layer in each cat

      real (kind=dbl_kind), parameter :: &
         hi_min = p01    ! minimum ice thickness allowed (m)

      real (kind=dbl_kind) :: &
         hin_max(0:ncat) ! category limits (m)

      character (len=35) :: c_hi_range(ncat)

!-------------------------------------------------------------------
! a note regarding hi_min and hin_max(0):
! both represent a minimum ice thickness.  hin_max(0) is
! intended to be used for particular numerical implementations
! of category conversions in the ice thickness distribution.
! hi_min is a more general purpose parameter, but is specifically
! for maintaining stability in the thermodynamics.
! hin_max(0) = 0.1 m for the delta function itd
! hin_max(0) = 0.0 m for linear remapping
!
! Also note that the upper limit on the thickest category
! is only used for the linear remapping scheme
! and it is not a true upper limit on the thickness
!-------------------------------------------------------------------

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_itd - initalize area fraction and thickness boundaries for ITD
!
! !INTERFACE:
!
      subroutine init_itd
!
! !DESCRIPTION:
!
! Initialize area fraction and thickness boundaries for the itd model
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
           n    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2

      real (kind=dbl_kind), dimension(5) :: wmo5 ! data for wmo itd
      real (kind=dbl_kind), dimension(6) :: wmo6 ! data for wmo itd
      real (kind=dbl_kind), dimension(7) :: wmo7 ! data for wmo itd

      character(len=8) :: c_hinmax1,c_hinmax2
      character(len=2) :: c_nc

      rncat = real(ncat, kind=dbl_kind)
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of three options.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001) 
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)] 
      ! 
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are 
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.    
      !
      ! The third option provides support for World Meteorological
      !  Organization classification based on thickness.  The full
      !  WMO thickness distribution is used if ncat = 7;  if ncat=5 
      !  or ncat = 6, some of the thinner categories are combined.
      ! For ncat = 5,  boundaries are         30, 70, 120, 200, >200 cm.
      ! For ncat = 6,  boundaries are     15, 30, 70, 120, 200, >200 cm.
      ! For ncat = 7,  boundaries are 10, 15, 30, 70, 120, 200, >200 cm.
      !-----------------------------------------------------------------

      if (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3/rncat
            cc2 = c15*cc1
            cc3 = c3

            hin_max(0) = c0     ! minimum ice thickness, m
         else
            ! delta function itd category limits
            cc1 = max(1.1_dbl_kind/rncat,c1*hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do n = 1, ncat
            x1 = real(n-1,kind=dbl_kind) / rncat
            hin_max(n) = hin_max(n-1) &
                       + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = rn * (d1 + (rn-c1)*d2)
         enddo

      elseif (kcatbound == 2) then  ! WMO standard

        if (ncat == 5) then
         ! thinnest 3 categories combined
         data wmo5 / 0.30_dbl_kind, 0.70_dbl_kind, &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo5(n)
         enddo
       elseif (ncat == 6) then
         ! thinnest 2 categories combined
         data wmo6 / 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo6(n)
         enddo
       elseif (ncat == 7) then
         ! all thickness categories 
         data wmo7 / 0.10_dbl_kind, 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo7(n)
         enddo
       else
         write (nu_diag,*) 'kcatbound=3 (WMO) must have ncat=5, 6 or 7'
         stop
       endif

      endif ! kcatbound

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         do n = 1, ncat
            write (nu_diag,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
            ! Write integer n to character string
            write (c_nc, '(i2)') n    

            ! Write hin_max to character string
            write (c_hinmax1, '(f5.3)') hin_max(n-1)
            write (c_hinmax2, '(f5.3)') hin_max(n)

            ! Save character string to write to history file
            c_hi_range(n)=c_hinmax1//'m < hi Cat '//c_nc//' < '// &
                          c_hinmax2//'m'
         enddo
         write (nu_diag,*) ' '
      endif

      !-----------------------------------------------------------------
      ! vectors identifying first and last layer in each category
      !-----------------------------------------------------------------
      ilyr1(1) = 1                       ! if nilyr  = 4
      ilyrn(1) = nilyr                   !   ilyr1 = { 1,5,9 }
      do n = 2, ncat                     !   ilyrn = { 4,8,12} etc
         ilyr1(n) = ilyrn(n-1) + 1
         ilyrn(n) = ilyrn(n-1) + nilyr
      enddo

      slyr1(1) = 1
      slyrn(1) = nslyr
      do n = 2, ncat
         slyr1(n) = slyrn(n-1) + 1
         slyrn(n) = slyrn(n-1) + nslyr
      enddo

      end subroutine init_itd

!=======================================================================
!BOP
!
! !IROUTINE: aggregate - aggregate ice state variables
!
! !INTERFACE:
!
      subroutine aggregate (nx_block, ny_block, &
                            aicen,    trcrn,    &
                            vicen,    vsnon,    &
                            eicen,    esnon,    &
                            aice,     trcr,     &
                            vice,     vsno,     &
                            eice,     esno,     &
                            aice0,    tmask,    &
                            ntrcr,    trcr_depend)
!
! !DESCRIPTION:
!
! Aggregate ice state variables over thickness categories.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ntrcr                 ! number of tracers in use

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(in) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(in) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(in) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tmask     ! land/boundary mask, thickness (T-cell)

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block),  &
         intent(out) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         eice  , & ! energy of melt. of ice           (J/m^2)
         esno  , & ! energy of melt. of snow layer    (J/m^2)
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr),  &
         intent(out) :: &
         trcr      ! ice tracers
!
!EOP
!
      integer (kind=int_kind) :: &
        icells                ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi, &              ! compressed indices in i/j directions
        indxj

      integer (kind=int_kind) :: &
        i, j, k, n, it, &
        ij                    ! combined i/j horizontal index

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
        atrcr      ! sum of aicen*trcrn or vicen*trcrn or vsnon*trcrn

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! tmask

         aice0(i,j) = c1
         aice (i,j) = c0
         vice (i,j) = c0
         vsno (i,j) = c0
         eice (i,j) = c0
         esno (i,j) = c0
      enddo
      enddo


      allocate (atrcr(icells,ntrcr))

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      atrcr(:,:) = c0

      do n = 1, ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aice(i,j) = aice(i,j) + aicen(i,j,n)
            vice(i,j) = vice(i,j) + vicen(i,j,n)
            vsno(i,j) = vsno(i,j) + vsnon(i,j,n)
         enddo                  ! ij

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then  ! ice area tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcr(ij,it) = atrcr(ij,it)  &
                                + trcrn(i,j,it,n)*aicen(i,j,n)
               enddo            ! ij

            elseif (trcr_depend(it) == 1) then  ! ice volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcr(ij,it) = atrcr(ij,it)  &
                                + trcrn(i,j,it,n)*vicen(i,j,n)
               enddo            ! ij

            elseif (trcr_depend(it) ==2) then ! snow volume tracer

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcr(ij,it) = atrcr(ij,it)  &
                                + trcrn(i,j,it,n)*vsnon(i,j,n)
               enddo            ! ij

            endif               ! trcr_depend
         enddo                  ! ntrcr

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               eice(i,j) = eice(i,j) + eicen(i,j,ilyr1(n)+k-1)
            enddo
         enddo                  ! nilyr

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               esno(i,j) = esno(i,j) + esnon(i,j,slyr1(n)+k-1)
            enddo
         enddo                  ! nslyr

      enddo                     ! ncat

      ! Open water fraction

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         aice0(i,j) = max (c1 - aice(i,j), c0)
      enddo                     ! ij

      ! Tracers

      call compute_tracers (nx_block,     ny_block,   &
                            icells,   indxi,   indxj, &
                            ntrcr,    trcr_depend,    &
                            atrcr(:,:), aice(:,:),    &
                            vice (:,:),   vsno(:,:),  &
                            trcr(:,:,:))

      deallocate (atrcr)

      end subroutine aggregate

!=======================================================================
!BOP
!
! !IROUTINE: aggregate_area - aggregate ice area
!
! !INTERFACE:
!
! Subprogram not used       subroutine aggregate_area (nx_block, ny_block,        &
! Subprogram not used                                  aicen,    aice,     aice0)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Aggregate ice area (but not other state variables) over thickness 
! Subprogram not used ! categories.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors: William H. Lipscomb, LANL
! Subprogram not used !          modified Jan 2004 by Clifford Chen, Fujitsu
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block  ! block dimensions
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
! Subprogram not used          aicen     ! concentration of ice
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
! Subprogram not used          aice, &   ! concentration of ice
! Subprogram not used          aice0     ! concentration of open water
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: i, j, n
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Aggregate
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       aice(:,:) = c0
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             aice(i,j) = aice(i,j) + aicen(i,j,n)
! Subprogram not used          enddo                  ! i
! Subprogram not used          enddo                  ! j
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used 
! Subprogram not used          ! open water fraction
! Subprogram not used          aice0(i,j) = max (c1 - aice(i,j), c0)
! Subprogram not used 
! Subprogram not used       enddo                     ! i
! Subprogram not used       enddo                     ! j
! Subprogram not used 
! Subprogram not used       end subroutine aggregate_area

!=======================================================================
!BOP
!
! !IROUTINE: rebin - rebins thicknesses into defined categories
!
! !INTERFACE:
!
! Subprogram not used       subroutine rebin (nx_block, ny_block,        &
! Subprogram not used                         icells,   indxi,    indxj, &
! Subprogram not used                         ntrcr,    trcr_depend,     &
! Subprogram not used                         aicen,    trcrn,           &
! Subprogram not used                         vicen,    vsnon,           &
! Subprogram not used                         eicen,    esnon,           &
! Subprogram not used                         l_stop,                    &
! Subprogram not used                         istop,    jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Rebins thicknesses into defined categories
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          icells            , & ! number of grid cells with ice
! Subprogram not used          ntrcr                 ! number of tracers in use
! Subprogram not used 
! Subprogram not used        integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxi, indxj      ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
! Subprogram not used          trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          aicen , & ! concentration of ice
! Subprogram not used          vicen , & ! volume per unit area of ice           (m)
! Subprogram not used          vsnon     ! volume per unit area of snow          (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          trcrn     ! ice tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          eicen     ! energy of melting for each ice layer  (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          esnon     ! energy of melting for each snow layer (J/m^2)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(out) :: &
! Subprogram not used          l_stop    ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(out) :: &
! Subprogram not used          istop, jstop    ! indices of grid cell where model aborts
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i,j          , & ! horizontal indices
! Subprogram not used          n            , & ! category index
! Subprogram not used          ij                ! combined horizontal index
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used          shiftflag          ! = .true. if ice must be shifted
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (icells,ncat) :: &
! Subprogram not used          donor              ! donor category index
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (icells,ncat) :: &
! Subprogram not used          daice          , & ! ice area transferred
! Subprogram not used          dvice          , & ! ice volume transferred
! Subprogram not used          hicen              ! ice thickness for each cat (m)
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used          do ij = 1, icells       ! aice(i,j) > puny
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used             donor(ij,n) = 0
! Subprogram not used             daice(ij,n) = c0
! Subprogram not used             dvice(ij,n) = c0
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Compute ice thickness.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used             if (aicen(i,j,n) > puny) then
! Subprogram not used                hicen(ij,n) = vicen(i,j,n) / aicen(i,j,n)
! Subprogram not used             else
! Subprogram not used                hicen(ij,n) = c0
! Subprogram not used             endif
! Subprogram not used          enddo                  ! ij
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! make sure thickness of cat 1 is at least hin_max(0)
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       do ij = 1, icells       ! aice(i,j) > puny
! Subprogram not used          i = indxi(ij)
! Subprogram not used          j = indxj(ij)
! Subprogram not used 
! Subprogram not used          if (aicen(i,j,1) > puny) then
! Subprogram not used             if (hicen(ij,1) <= hin_max(0) .and. hin_max(0) > c0 ) then
! Subprogram not used                aicen(i,j,1) = vicen(i,j,1) / hin_max(0)
! Subprogram not used                hicen(ij,1) = hin_max(0)
! Subprogram not used             endif
! Subprogram not used          endif
! Subprogram not used       enddo                     ! ij
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! If a category thickness is not in bounds, shift the
! Subprogram not used       ! entire area, volume, and energy to the neighboring category
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Move thin categories up
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat-1          ! loop over category boundaries
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! identify thicknesses that are too big
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          shiftflag = .false.
! Subprogram not used          do ij = 1, icells       ! aice(i,j) > puny
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used             if (aicen(i,j,n) > puny .and. &
! Subprogram not used                 hicen(ij,n) > hin_max(n)) then
! Subprogram not used                shiftflag = .true.
! Subprogram not used                donor(ij,n) = n
! Subprogram not used                daice(ij,n) = aicen(i,j,n)
! Subprogram not used                dvice(ij,n) = vicen(i,j,n)
! Subprogram not used             endif
! Subprogram not used          enddo                  ! ij
! Subprogram not used 
! Subprogram not used          if (shiftflag) then
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! shift ice between categories
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used             call shift_ice (nx_block, ny_block,    &
! Subprogram not used                             indxi,    indxj,       &
! Subprogram not used                             icells,                &
! Subprogram not used                             ntrcr,    trcr_depend, &
! Subprogram not used                             aicen,    trcrn,       &
! Subprogram not used                             vicen,    vsnon,       &
! Subprogram not used                             eicen,    esnon,       &
! Subprogram not used                             hicen,    donor,       &
! Subprogram not used                             daice,    dvice,       &
! Subprogram not used                             l_stop,                &
! Subprogram not used                             istop,    jstop)
! Subprogram not used 
! Subprogram not used             if (l_stop) return
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! reset shift parameters
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do ij = 1, icells       ! aice(i,j) > puny
! Subprogram not used             donor(ij,n) = 0
! Subprogram not used             daice(ij,n) = c0
! Subprogram not used             dvice(ij,n) = c0
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          endif                  ! shiftflag
! Subprogram not used 
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Move thick categories down
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = ncat-1, 1, -1      ! loop over category boundaries
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! identify thicknesses that are too small
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          shiftflag = .false.
! Subprogram not used          do ij = 1, icells       ! aice(i,j) > puny
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used             if (aicen(i,j,n+1) > puny .and. &
! Subprogram not used                 hicen(ij,n+1) <= hin_max(n)) then
! Subprogram not used                shiftflag = .true.
! Subprogram not used                donor(ij,n) = n+1
! Subprogram not used                daice(ij,n) = aicen(i,j,n+1)
! Subprogram not used                dvice(ij,n) = vicen(i,j,n+1)
! Subprogram not used             endif
! Subprogram not used          enddo                  ! ij
! Subprogram not used 
! Subprogram not used          if (shiftflag) then
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! shift ice between categories
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used             call shift_ice (nx_block, ny_block,    &
! Subprogram not used                             indxi,    indxj,       &
! Subprogram not used                             icells,                &
! Subprogram not used                             ntrcr,    trcr_depend, &
! Subprogram not used                             aicen,    trcrn,       &
! Subprogram not used                             vicen,    vsnon,       &
! Subprogram not used                             eicen,    esnon,       &
! Subprogram not used                             hicen,    donor,       &
! Subprogram not used                             daice,    dvice,       &
! Subprogram not used                             l_stop,                &
! Subprogram not used                             istop,    jstop)
! Subprogram not used 
! Subprogram not used             if (l_stop) return
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! reset shift parameters
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do ij = 1, icells       ! aice(i,j) > puny
! Subprogram not used             donor(ij,n) = 0
! Subprogram not used             daice(ij,n) = c0
! Subprogram not used             dvice(ij,n) = c0
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          endif                  ! shiftflag
! Subprogram not used 
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       end subroutine rebin

!=======================================================================
!BOP
!
! !IROUTINE: reduce_area - reduce area when ice melts for special case ncat=1
!
! !INTERFACE:
!
! Subprogram not used       subroutine reduce_area (nx_block, ny_block, &
! Subprogram not used                               ilo, ihi, jlo, jhi, &
! Subprogram not used                               tmask,              &
! Subprogram not used                               aicen,     vicen,   &
! Subprogram not used                               aicen_init,vicen_init)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Reduce area when ice melts for special case of ncat=1
! Subprogram not used !
! Subprogram not used ! Use CSM 1.0-like method of reducing ice area
! Subprogram not used ! when melting occurs: assume only half the ice volume
! Subprogram not used ! change goes to thickness decrease, the other half
! Subprogram not used ! to reduction in ice fraction
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors: C. M. Bitz, UW
! Subprogram not used ! modified by: Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi       ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          tmask     ! land/boundary mask, thickness (T-cell)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          aicen , & ! concentration of ice
! Subprogram not used          vicen     ! volume per unit area of ice          (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
! Subprogram not used          aicen_init, & ! old ice area for category 1 (m)
! Subprogram not used          vicen_init    ! old ice volume for category 1 (m)
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j        ! horizontal indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          hi0     , & ! initial hi
! Subprogram not used          hi1     , & ! current hi
! Subprogram not used          dhi         ! hi1 - hi0
! Subprogram not used 
! Subprogram not used       do j = jlo, jhi
! Subprogram not used       do i = ilo, ihi
! Subprogram not used          if (tmask(i,j)) then
! Subprogram not used 
! Subprogram not used             hi0 = c0
! Subprogram not used             if (aicen_init(i,j) > c0) &
! Subprogram not used                 hi0 = vicen_init(i,j) / aicen_init(i,j)
! Subprogram not used 
! Subprogram not used             hi1 = c0
! Subprogram not used             if (aicen(i,j) > c0) &
! Subprogram not used                 hi1 = vicen(i,j) / aicen(i,j)
! Subprogram not used 
! Subprogram not used             ! make sure thickness of cat 1 is at least hin_max(0)
! Subprogram not used             if (hi1 <= hin_max(0) .and. hin_max(0) > c0 ) then
! Subprogram not used                aicen(i,j) = vicen(i,j) / hin_max(0)
! Subprogram not used                hi1 = hin_max(0)
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             if (aicen(i,j) > c0) then
! Subprogram not used                dhi = hi1 - hi0
! Subprogram not used                if (dhi < c0) then
! Subprogram not used                   hi1  = vicen(i,j) / aicen(i,j)
! Subprogram not used                   aicen(i,j) = c2 * vicen(i,j) / (hi1 + hi0)
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          endif                  ! tmask
! Subprogram not used       enddo                     ! i
! Subprogram not used       enddo                     ! j
! Subprogram not used 
! Subprogram not used       end subroutine reduce_area

!=======================================================================
!BOP
!
! !IROUTINE: shift_ice - shift ice across category boundaries
!
! !INTERFACE:
!
! Subprogram not used       subroutine shift_ice (nx_block, ny_block,    &
! Subprogram not used                             indxi,    indxj,       &
! Subprogram not used                             icells,                &
! Subprogram not used                             ntrcr,    trcr_depend, &
! Subprogram not used                             aicen,    trcrn,       &
! Subprogram not used                             vicen,    vsnon,       &
! Subprogram not used                             eicen,    esnon,       &
! Subprogram not used                             hicen,    donor,       &
! Subprogram not used                             daice,    dvice,       &
! Subprogram not used                             l_stop,                &
! Subprogram not used                             istop,    jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Shift ice across category boundaries, conserving area, volume, and
! Subprogram not used ! energy.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          icells            , & ! number of ocean/ice cells
! Subprogram not used          ntrcr                 ! number of tracers in use
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxi             , & ! compressed indices in i/j directions
! Subprogram not used          indxj
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
! Subprogram not used          trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          aicen , & ! concentration of ice
! Subprogram not used          vicen , & ! volume per unit area of ice          (m)
! Subprogram not used          vsnon     ! volume per unit area of snow         (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          trcrn     ! ice tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          eicen     ! energy of melting for each ice layer (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          esnon     ! energy of melting for each snow layer (J/m^2)
! Subprogram not used 
! Subprogram not used       ! NOTE: Third index of donor, daice, dvice should be ncat-1,
! Subprogram not used       !       except that compilers would have trouble when ncat = 1 
! Subprogram not used       integer (kind=int_kind), dimension(icells,ncat), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          donor             ! donor category index
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(icells,ncat), &
! Subprogram not used            intent(inout) :: &
! Subprogram not used          daice         , & ! ice area transferred across boundary
! Subprogram not used          dvice         , & ! ice volume transferred across boundary
! Subprogram not used          hicen             ! ice thickness for each cat        (m)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(out) :: &
! Subprogram not used          l_stop    ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(out) :: &
! Subprogram not used          istop, jstop    ! indices of grid cell where model aborts
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j, m       , & ! horizontal indices
! Subprogram not used          n             , & ! thickness category index
! Subprogram not used          nr            , & ! receiver category
! Subprogram not used          nd            , & ! donor category
! Subprogram not used          k             , & ! ice layer index
! Subprogram not used          it            , & ! tracer index
! Subprogram not used          ilo,ihi,jlo,jhi   ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(icells,max_ntrcr,ncat) :: &
! Subprogram not used          atrcrn            ! aicen*trcrn
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          dvsnow        , & ! snow volume transferred
! Subprogram not used          desnow        , & ! snow energy transferred
! Subprogram not used          deice         , & ! ice energy transferred
! Subprogram not used          datrcr            ! aicen*train transferred
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (icells) :: &
! Subprogram not used         indxii       , & ! compressed indices for i/j directions
! Subprogram not used         indxjj       , &
! Subprogram not used         indxij
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used         ishift      , & ! number of cells with ice to transfer
! Subprogram not used         ij              ! combined i/j horizontal index
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used         daice_negative     , & ! true if daice < -puny
! Subprogram not used         dvice_negative     , & ! true if dvice < -puny
! Subprogram not used         daice_greater_aicen, & ! true if daice > aicen
! Subprogram not used         dvice_greater_vicen    ! true if dvice > vicen
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
! Subprogram not used          worka, &
! Subprogram not used          workb
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used       worka(:,:) = c0
! Subprogram not used       workb(:,:) = c0
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used          do it = 1, ntrcr
! Subprogram not used             if (trcr_depend(it) == 0) then ! ice area tracer
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used                   atrcrn(ij,it,n) = aicen(i,j,n)*trcrn(i,j,it,n)
! Subprogram not used                enddo
! Subprogram not used             elseif (trcr_depend(it) ==1) then  ! ice volume tracer
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used                   atrcrn(ij,it,n) = vicen(i,j,n)*trcrn(i,j,it,n)
! Subprogram not used                enddo
! Subprogram not used             elseif (trcr_depend(it) ==2) then  ! snow volume tracer
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used                   atrcrn(ij,it,n) = vsnon(i,j,n)*trcrn(i,j,it,n)
! Subprogram not used                enddo
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Check for daice or dvice out of range, allowing for roundoff error
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat-1
! Subprogram not used 
! Subprogram not used          daice_negative = .false.
! Subprogram not used          dvice_negative = .false.
! Subprogram not used          daice_greater_aicen = .false.
! Subprogram not used          dvice_greater_vicen = .false.
! Subprogram not used 
! Subprogram not used 
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used             if (donor(ij,n) > 0) then
! Subprogram not used                nd = donor(ij,n)
! Subprogram not used 
! Subprogram not used                if (daice(ij,n) < c0) then
! Subprogram not used                   if (daice(ij,n) > -puny*aicen(i,j,nd)) then
! Subprogram not used                      daice(ij,n) = c0 ! shift no ice
! Subprogram not used                      dvice(ij,n) = c0
! Subprogram not used                   else
! Subprogram not used                      daice_negative = .true.
! Subprogram not used                   endif
! Subprogram not used                endif
! Subprogram not used          
! Subprogram not used                if (dvice(ij,n) < c0) then
! Subprogram not used                   if (dvice(ij,n) > -puny*vicen(i,j,nd)) then   
! Subprogram not used                      daice(ij,n) = c0 ! shift no ice
! Subprogram not used                      dvice(ij,n) = c0
! Subprogram not used                   else
! Subprogram not used                      dvice_negative = .true.
! Subprogram not used                   endif
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                if (daice(ij,n) > aicen(i,j,nd)*(c1-puny)) then
! Subprogram not used                   if (daice(ij,n) < aicen(i,j,nd)*(c1+puny)) then
! Subprogram not used                      daice(ij,n) = aicen(i,j,nd)
! Subprogram not used                      dvice(ij,n) = vicen(i,j,nd)
! Subprogram not used                   else
! Subprogram not used                      daice_greater_aicen = .true.
! Subprogram not used                   endif
! Subprogram not used                endif    
! Subprogram not used 
! Subprogram not used                if (dvice(ij,n) > vicen(i,j,nd)*(c1-puny)) then
! Subprogram not used                   if (dvice(ij,n) < vicen(i,j,nd)*(c1+puny)) then
! Subprogram not used                      daice(ij,n) = aicen(i,j,nd)
! Subprogram not used                      dvice(ij,n) = vicen(i,j,nd)
! Subprogram not used                   else
! Subprogram not used                      dvice_greater_vicen = .true.
! Subprogram not used                   endif
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used             endif               ! donor > 0
! Subprogram not used          enddo                  ! ij
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! error messages
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (daice_negative) then
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used                if (donor(ij,n) > 0 .and.  &
! Subprogram not used                    daice(ij,n) <= -puny*aicen(i,j,nd)) then
! Subprogram not used                   write(nu_diag,*) ' '
! Subprogram not used                   write(nu_diag,*) 'shift_ice: negative daice'
! Subprogram not used                   write(nu_diag,*) 'i, j:', i, j
! Subprogram not used                   write(nu_diag,*) 'boundary, donor cat:', n, nd
! Subprogram not used                   write(nu_diag,*) 'daice =', daice(ij,n)
! Subprogram not used                   write(nu_diag,*) 'dvice =', dvice(ij,n)
! Subprogram not used                   l_stop = .true.
! Subprogram not used                   istop = i
! Subprogram not used                   jstop = j
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used          endif
! Subprogram not used          if (l_stop) return
! Subprogram not used 
! Subprogram not used          if (dvice_negative) then
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used                if (donor(ij,n) > 0 .and.  &
! Subprogram not used                    dvice(ij,n) <= -puny*vicen(i,j,nd)) then
! Subprogram not used                   write(nu_diag,*) ' '
! Subprogram not used                   write(nu_diag,*) 'shift_ice: negative dvice'
! Subprogram not used                   write(nu_diag,*) 'i, j:', i, j
! Subprogram not used                   write(nu_diag,*) 'boundary, donor cat:', n, nd
! Subprogram not used                   write(nu_diag,*) 'daice =', daice(ij,n)
! Subprogram not used                   write(nu_diag,*) 'dvice =', dvice(ij,n)
! Subprogram not used                   l_stop = .true.
! Subprogram not used                   istop = i
! Subprogram not used                   jstop = j
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used          endif
! Subprogram not used          if (l_stop) return
! Subprogram not used 
! Subprogram not used          if (daice_greater_aicen) then
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used                if (donor(ij,n) > 0) then
! Subprogram not used                   nd = donor(ij,n)
! Subprogram not used                   if (daice(ij,n) >= aicen(i,j,nd)*(c1+puny)) then
! Subprogram not used                      write(nu_diag,*) ' '
! Subprogram not used                      write(nu_diag,*) 'shift_ice: daice > aicen'
! Subprogram not used                      write(nu_diag,*) 'i, j:', i, j
! Subprogram not used                      write(nu_diag,*) 'boundary, donor cat:', n, nd
! Subprogram not used                      write(nu_diag,*) 'daice =', daice(ij,n)
! Subprogram not used                      write(nu_diag,*) 'aicen =', aicen(i,j,nd)
! Subprogram not used                      l_stop = .true.
! Subprogram not used                      istop = i
! Subprogram not used                      jstop = j
! Subprogram not used                   endif
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used          endif
! Subprogram not used          if (l_stop) return
! Subprogram not used 
! Subprogram not used          if (dvice_greater_vicen) then
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used                if (donor(ij,n) > 0) then
! Subprogram not used                   nd = donor(ij,n)
! Subprogram not used                   if (dvice(ij,n) >= vicen(i,j,nd)*(c1+puny)) then
! Subprogram not used                      write(nu_diag,*) ' '
! Subprogram not used                      write(nu_diag,*) 'shift_ice: dvice > vicen'
! Subprogram not used                      write(nu_diag,*) 'i, j:', i, j
! Subprogram not used                      write(nu_diag,*) 'boundary, donor cat:', n, nd
! Subprogram not used                      write(nu_diag,*) 'dvice =', dvice(ij,n)
! Subprogram not used                      write(nu_diag,*) 'vicen =', vicen(i,j,nd)
! Subprogram not used                      l_stop = .true.
! Subprogram not used                      istop = i
! Subprogram not used                      jstop = j
! Subprogram not used                   endif
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used          endif
! Subprogram not used          if (l_stop) return
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! transfer volume and energy between categories
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          ishift = 0
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used            if (daice(ij,n) > c0) then ! daice(n) can be < puny
! Subprogram not used              ishift = ishift + 1
! Subprogram not used              indxii(ishift) = i
! Subprogram not used              indxjj(ishift) = j
! Subprogram not used              indxij(ishift) = ij
! Subprogram not used            endif   ! tmask
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used          do ij = 1, ishift
! Subprogram not used             i = indxii(ij)
! Subprogram not used             j = indxjj(ij)
! Subprogram not used             m = indxij(ij)
! Subprogram not used 
! Subprogram not used             nd = donor(m,n)
! Subprogram not used             worka(i,j) = dvice(m,n) / vicen(i,j,nd)
! Subprogram not used             if (nd  ==  n) then
! Subprogram not used                nr = nd+1
! Subprogram not used             else                ! nd = n+1
! Subprogram not used                nr = n
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             aicen(i,j,nd) = aicen(i,j,nd) - daice(m,n)
! Subprogram not used             aicen(i,j,nr) = aicen(i,j,nr) + daice(m,n)
! Subprogram not used 
! Subprogram not used             vicen(i,j,nd) = vicen(i,j,nd) - dvice(m,n)
! Subprogram not used             vicen(i,j,nr) = vicen(i,j,nr) + dvice(m,n)
! Subprogram not used 
! Subprogram not used             dvsnow = vsnon(i,j,nd) * worka(i,j)
! Subprogram not used             vsnon(i,j,nd) = vsnon(i,j,nd) - dvsnow
! Subprogram not used             vsnon(i,j,nr) = vsnon(i,j,nr) + dvsnow
! Subprogram not used             workb(i,j) = dvsnow
! Subprogram not used 
! Subprogram not used          enddo                  ! ij
! Subprogram not used 
! Subprogram not used          do it = 1, ntrcr
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, ishift
! Subprogram not used                i = indxii(ij)
! Subprogram not used                j = indxjj(ij)
! Subprogram not used                m = indxij(ij)
! Subprogram not used 
! Subprogram not used                nd = donor(m,n)
! Subprogram not used                if (nd == n) then
! Subprogram not used                   nr = nd+1
! Subprogram not used                else             ! nd = n+1
! Subprogram not used                   nr = n
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                if (trcr_depend(it) == 0) then
! Subprogram not used                   datrcr = daice(m,n)*trcrn(i,j,it,nd)
! Subprogram not used                elseif (trcr_depend(it) == 1) then
! Subprogram not used                   datrcr = dvice(m,n)*trcrn(i,j,it,nd)
! Subprogram not used                elseif (trcr_depend(it) == 2) then
! Subprogram not used                   datrcr = workb(i,j)  *trcrn(i,j,it,nd)
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                atrcrn(m,it,nd) = atrcrn(m,it,nd) - datrcr
! Subprogram not used                atrcrn(m,it,nr) = atrcrn(m,it,nr) + datrcr
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! ntrcr
! Subprogram not used 
! Subprogram not used          do k = 1, nilyr
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, ishift
! Subprogram not used                i = indxii(ij)
! Subprogram not used                j = indxjj(ij)
! Subprogram not used                m = indxij(ij)
! Subprogram not used 
! Subprogram not used                nd = donor(m,n)
! Subprogram not used                if (nd == n) then
! Subprogram not used                   nr = nd+1
! Subprogram not used                else             ! nd = n+1
! Subprogram not used                   nr = n
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                deice = eicen(i,j,ilyr1(nd)+k-1) * worka(i,j)
! Subprogram not used                eicen(i,j,ilyr1(nd)+k-1) = &
! Subprogram not used                     eicen(i,j,ilyr1(nd)+k-1) - deice
! Subprogram not used                eicen(i,j,ilyr1(nr)+k-1) = &
! Subprogram not used                     eicen(i,j,ilyr1(nr)+k-1) + deice
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! nilyr
! Subprogram not used 
! Subprogram not used          do k = 1, nslyr
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, ishift
! Subprogram not used                i = indxii(ij)
! Subprogram not used                j = indxjj(ij)
! Subprogram not used                m = indxij(ij)
! Subprogram not used 
! Subprogram not used                nd = donor(m,n)
! Subprogram not used                if (nd == n) then
! Subprogram not used                   nr = nd+1
! Subprogram not used                else             ! nd = n+1
! Subprogram not used                   nr = n
! Subprogram not used                endif
! Subprogram not used 
! Subprogram not used                desnow = esnon(i,j,slyr1(nd)+k-1) * worka(i,j)
! Subprogram not used                esnon(i,j,slyr1(nd)+k-1) = &
! Subprogram not used                     esnon(i,j,slyr1(nd)+k-1) - desnow
! Subprogram not used                esnon(i,j,slyr1(nr)+k-1) = &
! Subprogram not used                     esnon(i,j,slyr1(nr)+k-1) + desnow
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! nslyr
! Subprogram not used 
! Subprogram not used       enddo                     ! boundaries, 1 to ncat-1
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Update ice thickness and tracers
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used 
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used             if (aicen(i,j,n) > puny) then
! Subprogram not used                hicen(ij,n)   = vicen (i,j,n)   / aicen(i,j,n)
! Subprogram not used             else
! Subprogram not used                hicen(ij,n)   = c0
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          call compute_tracers (nx_block,        ny_block,       &
! Subprogram not used                                icells,          indxi,   indxj, &
! Subprogram not used                                ntrcr,           trcr_depend,    &
! Subprogram not used                                atrcrn(:,:,n),   aicen(:,:,  n), &
! Subprogram not used                                vicen (:,:,  n), vsnon(:,:,  n), &
! Subprogram not used                                trcrn(:,:,:,n))
! Subprogram not used 
! Subprogram not used       enddo                     ! ncat
! Subprogram not used 
! Subprogram not used       end subroutine shift_ice

!=======================================================================
!BOP
!
! !IROUTINE: column_sum - sum field over all ice categories
!
! !INTERFACE:
!
! Subprogram not used       subroutine column_sum (nx_block, ny_block,       &
! Subprogram not used                              icells,   indxi,   indxj, &
! Subprogram not used                              nsum,                     &
! Subprogram not used                              xin,      xout)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! For each grid cell, sum field over all ice categories.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          nsum              , & ! number of categories/layers
! Subprogram not used          icells                ! number of ice/ocean grid cells
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxi,  indxj          ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,nsum), &
! Subprogram not used            intent(in) :: &
! Subprogram not used            xin              ! input field
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (icells), intent(out) :: &
! Subprogram not used            xout             ! output field
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used            i, j, ij     , & ! horizontal indices
! Subprogram not used            n                ! category/layer index
! Subprogram not used 
! Subprogram not used       do ij = 1, icells
! Subprogram not used          xout(ij) = c0
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do n = 1, nsum
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used             xout(ij) = xout(ij) + xin(i,j,n)
! Subprogram not used          enddo                  ! ij
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       end subroutine column_sum

!=======================================================================
!BOP
!
! !IROUTINE: column_conservation_check
!
! !INTERFACE:
!
! Subprogram not used       subroutine column_conservation_check (nx_block, ny_block,       &
! Subprogram not used                                             icells,   indxi,   indxj, &
! Subprogram not used                                             fieldid,                  &
! Subprogram not used                                             x1,       x2,             &
! Subprogram not used                                             max_err,  l_stop,         &
! Subprogram not used                                             istop,    jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! For each physical grid cell, check that initial and final values
! Subprogram not used ! of a conserved field are equal to within a small value.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          icells                ! number of ice/ocean grid cells
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxi,  indxj     ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(icells), intent(in) :: &
! Subprogram not used          x1            , & ! initial field
! Subprogram not used          x2                ! final field
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          max_err           ! max allowed error
! Subprogram not used 
! Subprogram not used       character (len=char_len), intent(in) :: &
! Subprogram not used          fieldid           ! field identifier
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(inout) :: &
! Subprogram not used          l_stop            ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(inout) :: &
! Subprogram not used          istop, jstop      ! indices of grid cell where model aborts
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          ij                    ! horizontal indices
! Subprogram not used 
! Subprogram not used       do ij = 1, icells
! Subprogram not used          if (abs (x2(ij)-x1(ij)) > max_err) then
! Subprogram not used             l_stop = .true.
! Subprogram not used             istop = indxi(ij)
! Subprogram not used             jstop = indxj(ij)
! Subprogram not used 
! Subprogram not used             write (nu_diag,*) ' '
! Subprogram not used             write (nu_diag,*) 'Conservation error: ', trim(fieldid)
! Subprogram not used             write (nu_diag,*) 'i, j =', istop, jstop
! Subprogram not used             write (nu_diag,*) 'Initial value =', x1(ij)
! Subprogram not used             write (nu_diag,*) 'Final value =',   x2(ij)
! Subprogram not used             write (nu_diag,*) 'Difference =', x2(ij) - x1(ij)
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine column_conservation_check

!=======================================================================
!BOP
!
! !IROUTINE: compute_tracers - compute tracer fields
!
! !INTERFACE:
!
      subroutine compute_tracers (nx_block, ny_block,       &
                                  icells,   indxi,   indxj, &
                                  ntrcr,    trcr_depend,    &
                                  atrcrn,   aicen,          &
                                  vicen,    vsnon,          &
                                  trcrn)
!
! !DESCRIPTION:
!
! Compute tracer fields.
! Given atrcrn = aicen*trcrn (or vicen*trcrn, vsnon*trcrn), compute trcrn.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!          
! !USES:
!
      use ice_state, only: nt_Tsfc
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of ice/ocean grid cells
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj       ! compressed i/j indices

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (icells,ntrcr), &
         intent(in) :: &
         atrcrn    ! aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr), &
         intent(out) :: &
         trcrn     ! ice tracers
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, it, ij       ! counting indices


      trcrn(:,:,:) = c0

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do it = 1, ntrcr
         if (it == nt_Tsfc) then      ! surface temperature
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (aicen(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / aicen(i,j)
               else
                  trcrn(i,j,it) = Tocnfrz
               endif
            enddo

         elseif (trcr_depend(it) == 0) then ! ice area tracers
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (aicen(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / aicen(i,j)
               else
                  trcrn(i,j,it) = c0
               endif
            enddo

         elseif (trcr_depend(it) == 1) then ! ice volume tracers
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (vicen(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / vicen(i,j)
               else
                  trcrn(i,j,it) = c0
               endif
            enddo

         elseif (trcr_depend(it) == 2) then ! snow volume tracers
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (vsnon(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / vsnon(i,j)
               else
                  trcrn(i,j,it) = c0
               endif
            enddo

         endif                  ! trcr_depend
      enddo                     ! ntrcr

      end subroutine compute_tracers

!=======================================================================
!BOP
!
! !IROUTINE: cleanup_itd - rebin if needed, eliminate small ice areas,
!                          and aggregate over categories
!
! !INTERFACE:
!
! Subprogram not used       subroutine cleanup_itd (nx_block,    ny_block,   &
! Subprogram not used                               ilo, ihi,    jlo, jhi,   &
! Subprogram not used                               dt,          ntrcr,      &
! Subprogram not used                               aicen,       trcrn,      &
! Subprogram not used                               vicen,       vsnon,      &
! Subprogram not used                               eicen,       esnon,      &
! Subprogram not used                               aice0,       aice,       &
! Subprogram not used                               trcr_depend, fresh,      &
! Subprogram not used                               fsalt,       fhocn,      &
! Subprogram not used                               fsoot,       tr_aero,    &
! Subprogram not used                               heat_capacity, l_stop,   &
! Subprogram not used                               istop,         jstop,    &
! Subprogram not used                               limit_aice_in)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Cleanup subroutine that rebins thickness categories if necessary,
! Subprogram not used !  eliminates very small ice areas while conserving mass and energy, 
! Subprogram not used !  aggregates state variables, and does a boundary call.  
! Subprogram not used ! It is a good idea to call this subroutine after the thermodynamics
! Subprogram not used !  (thermo_vertical/thermo_itd) and again after the dynamics 
! Subprogram not used !  (evp/transport/ridging).
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: & 
! Subprogram not used          nx_block, ny_block, & ! block dimensions 
! Subprogram not used          ilo,ihi,jlo,jhi   , & ! beginning and end of physical domain
! Subprogram not used          ntrcr                 ! number of tracers in use
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), intent(in) :: & 
! Subprogram not used          dt        ! time step 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),  &
! Subprogram not used          intent(inout) :: & 
! Subprogram not used          aicen , & ! concentration of ice 
! Subprogram not used          vicen , & ! volume per unit area of ice          (m) 
! Subprogram not used          vsnon     ! volume per unit area of snow         (m) 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat),  &
! Subprogram not used          intent(inout) :: & 
! Subprogram not used          trcrn     ! ice tracers 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),  &
! Subprogram not used          intent(inout) :: & 
! Subprogram not used          eicen     ! energy of melting for each ice layer (J/m^2) 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),  &
! Subprogram not used          intent(inout) :: & 
! Subprogram not used          esnon     ! energy of melting for each snow layer (J/m^2) 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),  &
! Subprogram not used          intent(inout) :: & 
! Subprogram not used          aice  , & ! total ice concentration
! Subprogram not used          aice0     ! concentration of open water 
! Subprogram not used      
! Subprogram not used       integer (kind=int_kind), dimension(max_ntrcr), intent(in) :: & 
! Subprogram not used          trcr_depend  ! tracer dependency information
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in) :: &
! Subprogram not used          tr_aero,      &
! Subprogram not used          heat_capacity   ! if false, ice and snow have zero heat capacity
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(out) :: &
! Subprogram not used          l_stop    ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(out) :: &
! Subprogram not used          istop, jstop ! indices of grid cell where model aborts
! Subprogram not used 
! Subprogram not used       ! ice-ocean fluxes (required for strict conservation)
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(inout), optional :: &
! Subprogram not used          fresh    , & ! fresh water flux to ocean (kg/m^2/s)
! Subprogram not used          fsalt    , & ! salt flux to ocean        (kg/m^2/s)
! Subprogram not used          fhocn        ! net heat flux to ocean     (W/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,n_aeromx), &
! Subprogram not used          intent(inout), optional :: &
! Subprogram not used          fsoot        ! soot flux to ocean        (kg/m^2/s)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in), optional ::   &
! Subprogram not used          limit_aice_in      ! if false, allow aice to be out of bounds
! Subprogram not used                             ! may want to allow this for unit tests
! Subprogram not used !    
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j             , & ! horizontal indices
! Subprogram not used          n                , & ! category index
! Subprogram not used          icells               ! number of grid cells with ice
! Subprogram not used 
! Subprogram not used        integer (kind=int_kind), dimension (nx_block*ny_block) :: &
! Subprogram not used          indxi, indxj      ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
! Subprogram not used          dfresh   , & ! zapped fresh water flux (kg/m^2/s)
! Subprogram not used          dfsalt   , & ! zapped salt flux   (kg/m^2/s)
! Subprogram not used          dfhocn       ! zapped energy flux ( W/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,n_aeromx) :: &
! Subprogram not used          dfsoot    ! zapped soot flux   (kg/m^2/s)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) ::   &
! Subprogram not used          limit_aice         ! if true, check for aice out of bounds
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (present(limit_aice_in)) then
! Subprogram not used          limit_aice = limit_aice_in
! Subprogram not used       else
! Subprogram not used          limit_aice = .true.
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Compute total ice area.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call aggregate_area (nx_block, ny_block, &
! Subprogram not used                            aicen(:,:,:), &
! Subprogram not used                            aice,     aice0)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       if (limit_aice) then  ! check for aice out of bounds
! Subprogram not used       
! Subprogram not used          do j = jlo,jhi
! Subprogram not used          do i = ilo,ihi
! Subprogram not used             if (aice(i,j) > c1+puny .or. aice(i,j) < -puny) then
! Subprogram not used                l_stop = .true.
! Subprogram not used                istop = i
! Subprogram not used                jstop = j
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used          if (l_stop) then      ! area out of bounds
! Subprogram not used             i = istop
! Subprogram not used             j = jstop
! Subprogram not used             write(nu_diag,*) ' '
! Subprogram not used             write(nu_diag,*) 'aggregate ice area out of bounds'
! Subprogram not used             write(nu_diag,*) 'my_task, i, j, aice:', &
! Subprogram not used                               my_task, i, j, aice(i,j)
! Subprogram not used             do n = 1, ncat
! Subprogram not used                write(nu_diag,*) 'n, aicen:', n, aicen(i,j,n)
! Subprogram not used             enddo
! Subprogram not used             return
! Subprogram not used          endif                  ! l_stop
! Subprogram not used       endif                     ! limit_aice
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Identify grid cells with ice.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       icells = 0
! Subprogram not used       do j = jlo,jhi
! Subprogram not used       do i = ilo,ihi
! Subprogram not used          if (aice(i,j) > puny) then
! Subprogram not used             icells = icells + 1
! Subprogram not used             indxi(icells) = i
! Subprogram not used             indxj(icells) = j
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Make sure ice in each category is within its thickness bounds.
! Subprogram not used       ! NOTE: The rebin subroutine is needed only in the rare cases
! Subprogram not used       !       when the linear_itd subroutine cannot transfer ice
! Subprogram not used       !       correctly (e.g., very fast ice growth).
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       call rebin (nx_block,     ny_block,       &
! Subprogram not used                   icells,       indxi, indxj,   &
! Subprogram not used                   ntrcr,        trcr_depend,    &
! Subprogram not used                   aicen(:,:,:), trcrn(:,:,:,:), &
! Subprogram not used                   vicen(:,:,:), vsnon(:,:,:),   &
! Subprogram not used                   eicen(:,:,:), esnon(:,:,:),   &
! Subprogram not used                   l_stop,                       &
! Subprogram not used                   istop,      jstop)
! Subprogram not used 
! Subprogram not used       if (l_stop) return
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zero out ice categories with very small areas.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (limit_aice) then
! Subprogram not used          call zap_small_areas (nx_block,     ny_block,       &
! Subprogram not used                                ilo, ihi,     jlo, jhi,       &
! Subprogram not used                                dt,           ntrcr,          &
! Subprogram not used                                aice,         aice0,          &
! Subprogram not used                                aicen(:,:,:), trcrn(:,:,:,:), &
! Subprogram not used                                vicen(:,:,:), vsnon(:,:,:),   &
! Subprogram not used                                eicen(:,:,:), esnon(:,:,:),   &
! Subprogram not used                                dfresh,       dfsalt,         &
! Subprogram not used                                dfhocn,       dfsoot,         &
! Subprogram not used                                tr_aero,                      &
! Subprogram not used                                l_stop,                       &
! Subprogram not used                                istop,        jstop)
! Subprogram not used          if (l_stop) return
! Subprogram not used       endif   ! l_limit_aice
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Update ice-ocean fluxes for strict conservation
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (present(fresh)) &
! Subprogram not used            fresh     (:,:) = fresh(:,:)      + dfresh(:,:) 
! Subprogram not used       if (present(fsalt)) &
! Subprogram not used            fsalt     (:,:) = fsalt(:,:)      + dfsalt(:,:)
! Subprogram not used       if (present(fhocn)) &
! Subprogram not used            fhocn     (:,:) = fhocn(:,:)      + dfhocn(:,:)
! Subprogram not used       if (present(fsoot)) &
! Subprogram not used            fsoot   (:,:,:) = fsoot(:,:,:)    + dfsoot(:,:,:)
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       ! If using zero-layer model (no heat capacity), check that the 
! Subprogram not used       ! energy of snow and ice is correct. 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (.not. heat_capacity) then
! Subprogram not used 
! Subprogram not used          call zerolayer_check(nx_block,    ny_block,   &
! Subprogram not used                               icells,  indxi,   indxj, &
! Subprogram not used                               aicen,                   &
! Subprogram not used                               vicen,       vsnon,      &
! Subprogram not used                               eicen,       esnon,      &
! Subprogram not used                               l_stop,                  &
! Subprogram not used                               istop,       jstop)
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       end subroutine cleanup_itd

!=======================================================================
!BOP
!
! !IROUTINE: zap_small_areas - eliminate very small ice areas
!
! !INTERFACE:
!
! Subprogram not used       subroutine zap_small_areas (nx_block, ny_block, &
! Subprogram not used                                   ilo, ihi, jlo, jhi, &
! Subprogram not used                                   dt,       ntrcr,    &
! Subprogram not used                                   aice,     aice0,    &
! Subprogram not used                                   aicen,    trcrn,    &
! Subprogram not used                                   vicen,    vsnon,    &
! Subprogram not used                                   eicen,    esnon,    &
! Subprogram not used                                   dfresh,   dfsalt,   &
! Subprogram not used                                   dfhocn,   dfsoot,   &
! Subprogram not used                                   tr_aero,            &
! Subprogram not used                                   l_stop,             &
! Subprogram not used                                   istop,    jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! For each ice category in each grid cell, remove ice if the fractional
! Subprogram not used ! area is less than puny.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_state, only: nt_Tsfc, nt_aero
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          nx_block, ny_block, & ! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi   , & ! beginning and end of physical domain
! Subprogram not used          ntrcr                 ! number of tracers in use
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in) :: &
! Subprogram not used          dt                    ! time step
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          aice     , & ! total ice concentration
! Subprogram not used          aice0        ! concentration of open water
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          aicen    , & ! concentration of ice
! Subprogram not used          vicen    , & ! volume per unit area of ice          (m)
! Subprogram not used          vsnon        ! volume per unit area of snow         (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          eicen        ! energy of melting for each ice layer  (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          esnon        ! energy of melting for each snow layer (J/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
! Subprogram not used          intent(inout) :: &
! Subprogram not used          trcrn        ! ice tracers
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), &
! Subprogram not used          intent(out) :: &
! Subprogram not used          dfresh   , & ! zapped fresh water flux (kg/m^2/s)
! Subprogram not used          dfsalt   , & ! zapped salt flux   (kg/m^2/s)
! Subprogram not used          dfhocn       ! zapped energy flux ( W/m^2)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,n_aeromx), &
! Subprogram not used          intent(out) :: &
! Subprogram not used          dfsoot    ! zapped soot flux   (kg/m^2/s)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in) :: &
! Subprogram not used          tr_aero
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(out) :: &
! Subprogram not used          l_stop   ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(out) :: &
! Subprogram not used          istop, jstop ! indices of grid cell where model aborts
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i,j, n, k, it  , & ! counting indices
! Subprogram not used          icells         , & ! number of cells with ice to zap
! Subprogram not used          ij                 ! combined i/j horizontal index
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block) :: &
! Subprogram not used         indxi       , & ! compressed indices for i/j directions
! Subprogram not used         indxj
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: xtmp      ! temporary variable
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used       dfresh(:,:) = c0
! Subprogram not used       dfsalt(:,:) = c0
! Subprogram not used       dfhocn(:,:) = c0
! Subprogram not used       dfsoot(:,:,:) = c0
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zap categories with very small areas.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Count categories to be zapped.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          icells = 0
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             if (aicen(i,j,n) < -puny) then
! Subprogram not used                write (nu_diag,*) 'Zap ice: negative ice area'
! Subprogram not used                write (nu_diag,*) 'i, j, n, aicen =', &
! Subprogram not used                                   i, j, n, aicen(i,j,n)
! Subprogram not used                l_stop = .true.
! Subprogram not used                istop = i
! Subprogram not used                jstop = j
! Subprogram not used                return
! Subprogram not used             elseif ((aicen(i,j,n) >= -puny .and. aicen(i,j,n) < c0) .or. &
! Subprogram not used                     (aicen(i,j,n) > c0 .and. aicen(i,j,n) <= puny)) then
! Subprogram not used                icells = icells + 1
! Subprogram not used                indxi(icells) = i
! Subprogram not used                indxj(icells) = j
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zap ice energy and use ocean heat to melt ice
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do k = 1, nilyr
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, icells
! Subprogram not used                i = indxi(ij)
! Subprogram not used                j = indxj(ij)
! Subprogram not used 
! Subprogram not used                xtmp = eicen(i,j,ilyr1(n)+k-1) / dt ! < 0
! Subprogram not used                dfhocn(i,j) = dfhocn(i,j) + xtmp
! Subprogram not used                eicen(i,j,ilyr1(n)+k-1) = c0
! Subprogram not used 
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! k
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zap snow energy and use ocean heat to melt snow
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do k = 1, nslyr
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, icells
! Subprogram not used                i = indxi(ij)
! Subprogram not used                j = indxj(ij)
! Subprogram not used 
! Subprogram not used                xtmp = esnon(i,j,slyr1(n)+k-1) / dt ! < 0
! Subprogram not used                dfhocn(i,j) = dfhocn(i,j) + xtmp
! Subprogram not used                esnon(i,j,slyr1(n)+k-1) = c0
! Subprogram not used 
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! k
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zap ice and snow volume, add water and salt to ocean
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used          do ij = 1, icells
! Subprogram not used             i = indxi(ij)
! Subprogram not used             j = indxj(ij)
! Subprogram not used 
! Subprogram not used             xtmp = (rhoi*vicen(i,j,n) + rhos*vsnon(i,j,n)) / dt
! Subprogram not used             dfresh(i,j) = dfresh(i,j) + xtmp
! Subprogram not used 
! Subprogram not used             xtmp = rhoi*vicen(i,j,n)*ice_ref_salinity*p001 / dt
! Subprogram not used             dfsalt(i,j) = dfsalt(i,j) + xtmp
! Subprogram not used 
! Subprogram not used             aice0(i,j) = aice0(i,j) + aicen(i,j,n)
! Subprogram not used             aicen(i,j,n) = c0
! Subprogram not used             vicen(i,j,n) = c0
! Subprogram not used             vsnon(i,j,n) = c0
! Subprogram not used             trcrn(i,j,nt_Tsfc,n) = Tocnfrz
! Subprogram not used 
! Subprogram not used          enddo                  ! ij
! Subprogram not used 
! Subprogram not used          if (tr_aero) then
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used           do ij = 1, icells
! Subprogram not used            i = indxi(ij)
! Subprogram not used            j = indxj(ij)
! Subprogram not used            do it=1,n_aero
! Subprogram not used             xtmp &
! Subprogram not used               = (vsnon(i,j,n)*(trcrn(i,j,nt_aero  +4*(it-1),n)   &
! Subprogram not used                               +trcrn(i,j,nt_aero+1+4*(it-1),n))  &
! Subprogram not used               +  vicen(i,j,n)*(trcrn(i,j,nt_aero+2+4*(it-1),n)   &
! Subprogram not used                               +trcrn(i,j,nt_aero+3+4*(it-1),n))) &
! Subprogram not used               / dt
! Subprogram not used             dfsoot(i,j,it) = dfsoot(i,j,it) + xtmp
! Subprogram not used            enddo                 ! n
! Subprogram not used           enddo                  ! ij
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zap tracers
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used          
! Subprogram not used          if (ntrcr >= 2) then
! Subprogram not used             do it = 1, ntrcr   ! this assumes nt_Tsfc = 1
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used                   trcrn(i,j,it,n) = c0
! Subprogram not used                enddo
! Subprogram not used             enddo
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Count cells with excess ice (aice > c1) due to roundoff errors.
! Subprogram not used       ! Zap a little ice in each category so that aice = c1.
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       icells = 0
! Subprogram not used       do j = jlo, jhi
! Subprogram not used       do i = ilo, ihi
! Subprogram not used          if (aice(i,j) > (c1+puny)) then
! Subprogram not used             write (nu_diag,*) 'Zap ice: excess ice area'
! Subprogram not used             write (nu_diag,*) 'i, j, aice =', &
! Subprogram not used                                i, j, aice(i,j)
! Subprogram not used             l_stop = .true.
! Subprogram not used             istop = i
! Subprogram not used             jstop = j
! Subprogram not used             return
! Subprogram not used          elseif (aice(i,j) > c1 .and. aice(i,j) < (c1+puny)) then
! Subprogram not used             icells = icells + 1
! Subprogram not used             indxi(icells) = i
! Subprogram not used             indxj(icells) = j
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------- 
! Subprogram not used       ! Zap ice energy and use ocean heat to melt ice 
! Subprogram not used       !----------------------------------------------------------------- 
! Subprogram not used        
! Subprogram not used          do k = 1, nilyr 
! Subprogram not used !DIR$ CONCURRENT !Cray 
! Subprogram not used !cdir nodep      !NEC 
! Subprogram not used !ocl novrec      !Fujitsu 
! Subprogram not used             do ij = 1, icells 
! Subprogram not used                i = indxi(ij) 
! Subprogram not used                j = indxj(ij) 
! Subprogram not used  
! Subprogram not used                xtmp = eicen(i,j,ilyr1(n)+k-1)  &
! Subprogram not used                     * (aice(i,j)-c1)/aice(i,j) / dt ! < 0 
! Subprogram not used                dfhocn(i,j) = dfhocn(i,j) + xtmp 
! Subprogram not used                eicen(i,j,ilyr1(n)+k-1) = eicen(i,j,ilyr1(n)+k-1) &
! Subprogram not used                                         * (c1/aice(i,j))
! Subprogram not used  
! Subprogram not used             enddo               ! ij 
! Subprogram not used          enddo                  ! k 
! Subprogram not used  
! Subprogram not used       !----------------------------------------------------------------- 
! Subprogram not used       ! Zap snow energy and use ocean heat to melt snow 
! Subprogram not used       !----------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used          do k = 1, nslyr 
! Subprogram not used !DIR$ CONCURRENT !Cray 
! Subprogram not used !cdir nodep      !NEC 
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, icells
! Subprogram not used                i = indxi(ij) 
! Subprogram not used                j = indxj(ij) 
! Subprogram not used  
! Subprogram not used                xtmp = esnon(i,j,slyr1(n)+k-1)  &
! Subprogram not used                     * (aice(i,j)-c1)/aice(i,j) / dt ! < 0 
! Subprogram not used                dfhocn(i,j) = dfhocn(i,j) + xtmp 
! Subprogram not used                esnon(i,j,slyr1(n)+k-1) = esnon(i,j,slyr1(n)+k-1) &
! Subprogram not used                                         *(c1/aice(i,j))
! Subprogram not used  
! Subprogram not used             enddo               ! ij
! Subprogram not used          enddo                  ! k
! Subprogram not used  
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Zap ice and snow volume, add water and salt to ocean
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray 
! Subprogram not used !cdir nodep      !NEC 
! Subprogram not used !ocl novrec      !Fujitsu 
! Subprogram not used          do ij = 1, icells 
! Subprogram not used             i = indxi(ij) 
! Subprogram not used             j = indxj(ij) 
! Subprogram not used  
! Subprogram not used             xtmp = (rhoi*vicen(i,j,n) + rhos*vsnon(i,j,n)) &
! Subprogram not used                  * (aice(i,j)-c1)/aice(i,j) / dt 
! Subprogram not used             dfresh(i,j) = dfresh(i,j) + xtmp 
! Subprogram not used  
! Subprogram not used             xtmp = rhoi*vicen(i,j,n)*ice_ref_salinity*p001 &
! Subprogram not used                  * (aice(i,j)-c1)/aice(i,j) / dt
! Subprogram not used             dfsalt(i,j) = dfsalt(i,j) + xtmp 
! Subprogram not used  
! Subprogram not used             aicen(i,j,n) = aicen(i,j,n) * (c1/aice(i,j)) 
! Subprogram not used             vicen(i,j,n) = vicen(i,j,n) * (c1/aice(i,j)) 
! Subprogram not used             vsnon(i,j,n) = vsnon(i,j,n) * (c1/aice(i,j))
! Subprogram not used  
! Subprogram not used          enddo                  ! ij
! Subprogram not used 
! Subprogram not used       ! Note: Tracers are unchanged.
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used          if (tr_aero) then
! Subprogram not used           do ij = 1, icells
! Subprogram not used            i = indxi(ij)
! Subprogram not used            j = indxj(ij)
! Subprogram not used            do it=1,n_aero
! Subprogram not used             xtmp &
! Subprogram not used               = (vsnon(i,j,n)*(trcrn(i,j,nt_aero  +4*(it-1),n)   &
! Subprogram not used                               +trcrn(i,j,nt_aero+1+4*(it-1),n))  &
! Subprogram not used               +  vicen(i,j,n)*(trcrn(i,j,nt_aero+2+4*(it-1),n)   &
! Subprogram not used                               +trcrn(i,j,nt_aero+3+4*(it-1),n))) &
! Subprogram not used               * (aice(i,j)-c1)/aice(i,j) / dt
! Subprogram not used             dfsoot(i,j,it) = dfsoot(i,j,it) + xtmp
! Subprogram not used            enddo                 ! n
! Subprogram not used           enddo                  ! ij
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo                     ! n
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Correct aice
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray 
! Subprogram not used !cdir nodep      !NEC 
! Subprogram not used !ocl novrec      !Fujitsu 
! Subprogram not used       do ij = 1, icells 
! Subprogram not used          i = indxi(ij) 
! Subprogram not used          j = indxj(ij) 
! Subprogram not used          aice(i,j) = c1
! Subprogram not used          aice0(i,j) = c0
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine zap_small_areas

!=======================================================================
!BOP
!
! !IROUTINE: zerolayer_check - check that snow and ice energy is
!                         correct when using zero layer thermodynamics
!
! !INTERFACE:
!
! Subprogram not used       subroutine zerolayer_check (nx_block,    ny_block,   &
! Subprogram not used                                   icells,  indxi,   indxj, &
! Subprogram not used                                   aicen,                   &
! Subprogram not used                                   vicen,       vsnon,      &
! Subprogram not used                                   eicen,       esnon,      &
! Subprogram not used                                   l_stop,                  &
! Subprogram not used                                   istop,       jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Checks that the snow and ice energy in the zero layer thermodynamics
! Subprogram not used ! model still agrees with the snow and ice volume.
! Subprogram not used ! If there is an error, the model will abort.
! Subprogram not used ! This subroutine is only called if heat_capacity = .false.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Alison McLaren, Met Office
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: & 
! Subprogram not used          nx_block, ny_block, & ! block dimensions 
! Subprogram not used          icells                ! number of grid cells with ice
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block), &
! Subprogram not used          intent(in) :: &
! Subprogram not used          indxi, indxj      ! compressed i/j indices
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),  &
! Subprogram not used          intent(inout) :: & 
! Subprogram not used          aicen , & ! concentration of ice 
! Subprogram not used          vicen , & ! volume per unit area of ice          (m) 
! Subprogram not used          vsnon     ! volume per unit area of snow         (m) 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),  &
! Subprogram not used          intent(in) :: & 
! Subprogram not used          eicen     ! energy of melting for each ice layer (J/m^2) 
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),  &
! Subprogram not used          intent(in) :: & 
! Subprogram not used          esnon     ! energy of melting for each snow layer (J/m^2) 
! Subprogram not used       
! Subprogram not used       logical (kind=log_kind), intent(out) :: &
! Subprogram not used          l_stop    ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(out) :: &
! Subprogram not used          istop, jstop ! indices of grid cell where model aborts
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          i, j             , & ! horizontal indices
! Subprogram not used          n                , & ! category index
! Subprogram not used          ij                   ! combined horizontal index
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), parameter :: &
! Subprogram not used          max_error = puny*Lfresh*rhos ! max error in zero layer energy check
! Subprogram not used                                       ! (so max volume error = puny)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used          ice_energy_correct  , & ! zero layer ice energy check
! Subprogram not used          snow_energy_correct     ! zero layer snow energy check
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
! Subprogram not used          worka, &
! Subprogram not used          workb
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used       ! Initialize
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used       worka(:,:) = c0
! Subprogram not used       workb(:,:) = c0
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       ! Calculate difference between ice and snow energies and the
! Subprogram not used       ! energy values derived from the ice and snow volumes
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       ice_energy_correct  = .true.
! Subprogram not used       snow_energy_correct = .true.
! Subprogram not used 
! Subprogram not used       do n=1,ncat
! Subprogram not used 
! Subprogram not used          do ij=1,icells
! Subprogram not used             i=indxi(ij)
! Subprogram not used             j=indxj(ij)
! Subprogram not used 
! Subprogram not used             worka(i,j) = eicen(i,j,n) + rhoi * Lfresh * vicen(i,j,n)
! Subprogram not used             workb(i,j) = esnon(i,j,n) + rhos * Lfresh * vsnon(i,j,n)
! Subprogram not used 
! Subprogram not used             if(abs(worka(i,j)) > max_error) then
! Subprogram not used                ice_energy_correct = .false.
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             if(abs(workb(i,j)) > max_error) then
! Subprogram not used                snow_energy_correct = .false.
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used       ! If there is a problem, abort with error message
! Subprogram not used       !----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (.not. ice_energy_correct) then
! Subprogram not used 
! Subprogram not used             do ij=1,icells
! Subprogram not used                i=indxi(ij)
! Subprogram not used                j=indxj(ij)
! Subprogram not used 
! Subprogram not used                if(abs(worka(i,j)) > max_error) then
! Subprogram not used                   write(nu_diag,*) ' '
! Subprogram not used                   write(nu_diag,*) &
! Subprogram not used                     'zerolayer check - wrong ice energy'
! Subprogram not used                   write(nu_diag,*) 'i, j, n:', i,j,n
! Subprogram not used                   write(nu_diag,*) 'eicen =', eicen(i,j,n)
! Subprogram not used                   write(nu_diag,*) 'error=',  worka(i,j)
! Subprogram not used                   write(nu_diag,*) 'vicen =', vicen(i,j,n)
! Subprogram not used                   write(nu_diag,*) 'aicen =', aicen(i,j,n)
! Subprogram not used                   l_stop = .true.
! Subprogram not used                   istop = i
! Subprogram not used                   jstop = j
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used          if (l_stop) return
! Subprogram not used 
! Subprogram not used          if (.not. snow_energy_correct) then
! Subprogram not used 
! Subprogram not used             do ij=1,icells
! Subprogram not used                i=indxi(ij)
! Subprogram not used                j=indxj(ij)
! Subprogram not used 
! Subprogram not used                if(abs(workb(i,j)) > max_error) then
! Subprogram not used                   write(nu_diag,*) ' '
! Subprogram not used                   write(nu_diag,*) &
! Subprogram not used                     'zerolayer_check - wrong snow energy'
! Subprogram not used                   write(nu_diag,*) 'i, j, n:', i,j,n
! Subprogram not used                   write(nu_diag,*) 'esnon =', esnon(i,j,n)
! Subprogram not used                   write(nu_diag,*) 'error=',  workb(i,j)
! Subprogram not used                   write(nu_diag,*) 'vsnon =', vsnon(i,j,n)
! Subprogram not used                   write(nu_diag,*) 'aicen =', aicen(i,j,n)
! Subprogram not used                   l_stop = .true.
! Subprogram not used                   istop = i
! Subprogram not used                   jstop = j
! Subprogram not used                   return
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo  ! ncat
! Subprogram not used 
! Subprogram not used       end subroutine zerolayer_check

!=======================================================================

      end module ice_itd

!=======================================================================









