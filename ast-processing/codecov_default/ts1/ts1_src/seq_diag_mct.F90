!===============================================================================
! SVN $Id: seq_diag_mct.F90 61512 2014-06-26 21:59:35Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq5_0_14/driver/seq_diag_mct.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_diag_mod -- computes spatial \& time averages of fluxed quatities
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.
!
! !REMARKS:
!    CESM sign convention for fluxes is positive downward with hierarchy being
!       atm/glc/lnd/rof/ice/ocn
!    Sign convention:
!       positive value <=> the model is gaining water, heat, momentum, etc.
!    Unit convention:
!       heat flux     ~ W/m^2
!       momentum flux ~ N/m^2
!       water flux    ~ (kg/s)/m^2
!       salt  flux    ~ (kg/s)/m^2
!
! !REVISION HISTORY:
!    2012-aug-20 - T. Craig    - add rof component
!    2008-jul-10 - T. Craig    - updated budget implementation
!    2007-may-07 - B. Kauffman - initial port to cpl7. 
!    2002-nov-21 - R. Jacob    - initial port to cpl6. 
!    199x-mmm-dd - B. Kauffman - original version in cpl4.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diag_mct
  
! !USES:

   use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
   use shr_kind_mod, only: i8 => shr_kind_i8,  cl=>shr_kind_cl
   use shr_sys_mod       ! system calls
   use shr_mpi_mod       ! mpi wrappers
   use shr_const_mod     ! shared constants
   use mct_mod           ! mct wrappers
   use esmf

   use seq_comm_mct  ! mpi comm groups & related
   use seq_timemgr_mod
   use component_type_mod

   implicit none
   save
   private

! !PUBLIC TYPES:

   ! none

!PUBLIC MEMBER FUNCTIONS:

   public seq_diag_zero_mct
   public seq_diag_atm_mct
   public seq_diag_lnd_mct
   public seq_diag_rof_mct
   public seq_diag_glc_mct
   public seq_diag_ocn_mct
   public seq_diag_ice_mct
   public seq_diag_accum_mct
   public seq_diag_sum0_mct
   public seq_diag_print_mct
   public seq_diag_avect_mct
   public seq_diag_avdiff_mct

!EOP

   !----------------------------------------------------------------------------
   ! Local data
   !----------------------------------------------------------------------------

   !----- local constants -----
   real(r8),parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
   &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
   &    (shr_const_ocn_ref_sal*shr_const_latice)

   !--- C for component ---
   !--- "r" is recieve in the coupler, "s" is send from the coupler

   integer(in),parameter :: c_size = 22

   integer(in),parameter :: c_atm_as   = 1 ! model index: atm
   integer(in),parameter :: c_atm_ar   = 2 ! model index: atm
   integer(in),parameter :: c_inh_is   = 3 ! model index: ice, northern
   integer(in),parameter :: c_inh_ir   = 4 ! model index: ice, northern
   integer(in),parameter :: c_ish_is   = 5 ! model index: ice, southern
   integer(in),parameter :: c_ish_ir   = 6 ! model index: ice, southern
   integer(in),parameter :: c_lnd_ls   = 7 ! model index: lnd
   integer(in),parameter :: c_lnd_lr   = 8 ! model index: lnd
   integer(in),parameter :: c_ocn_os   = 9 ! model index: ocn
   integer(in),parameter :: c_ocn_or   =10 ! model index: ocn
   integer(in),parameter :: c_rof_rs   =11 ! model index: rof
   integer(in),parameter :: c_rof_rr   =12 ! model index: rof
   integer(in),parameter :: c_glc_gs   =13 ! model index: glc
   integer(in),parameter :: c_glc_gr   =14 ! model index: glc
   ! --- on atm grid ---
   integer(in),parameter :: c_inh_as   =15 ! model index: ice, northern
   integer(in),parameter :: c_inh_ar   =16 ! model index: ice, northern
   integer(in),parameter :: c_ish_as   =17 ! model index: ice, southern
   integer(in),parameter :: c_ish_ar   =18 ! model index: ice, southern
   integer(in),parameter :: c_lnd_as   =19 ! model index: lnd
   integer(in),parameter :: c_lnd_ar   =20 ! model index: lnd
   integer(in),parameter :: c_ocn_as   =21 ! model index: ocn
   integer(in),parameter :: c_ocn_ar   =22 ! model index: ocn

   character(len=8),parameter :: cname(c_size) = &
      (/' c2a_atm',' a2c_atm',' c2i_inh',' i2c_inh',' c2i_ish',' i2c_ish', &
        ' c2l_lnd',' l2c_lnd',' c2o_ocn',' o2c_ocn',' c2r_rof',' r2c_rof', &
        ' c2g_glc',' g2c_glc', &
        ' c2a_inh',' a2c_inh',' c2a_ish',' a2c_ish', &
        ' c2a_lnd',' a2c_lnd',' c2a_ocn',' a2c_ocn' /)

   !--- F for field ---

   integer(in),parameter :: f_size = 17
   integer(in),parameter :: f_a    = 1    ! index for area
   integer(in),parameter :: f_h    = 2    ! 1st index for heat
   integer(in),parameter :: f_w    = 11   ! 1st index for water

   integer(in),parameter :: f_area    = 1 ! area (wrt to unit sphere)
   integer(in),parameter :: f_hfrz    = 2 ! heat : latent, freezing
   integer(in),parameter :: f_hmelt   = 3 ! heat : latent, melting
   integer(in),parameter :: f_hswnet  = 4 ! heat : short wave, net
   integer(in),parameter :: f_hlwdn   = 5 ! heat : longwave down
   integer(in),parameter :: f_hlwup   = 6 ! heat : longwave up
   integer(in),parameter :: f_hlatv   = 7 ! heat : latent, vaporization
   integer(in),parameter :: f_hlatf   = 8 ! heat : latent, fusion, snow       
   integer(in),parameter :: f_hioff   = 9 ! heat : latent, fusion, frozen runoff
   integer(in),parameter :: f_hsen    =10 ! heat : sensible
   integer(in),parameter :: f_wfrz    =11 ! water: freezing
   integer(in),parameter :: f_wmelt   =12 ! water: melting
   integer(in),parameter :: f_wrain   =13 ! water: precip, liquid
   integer(in),parameter :: f_wsnow   =14 ! water: precip, frozen
   integer(in),parameter :: f_wevap   =15 ! water: evaporation
   integer(in),parameter :: f_wroff   =16 ! water: runoff/flood
   integer(in),parameter :: f_wioff   =17 ! water: frozen runoff

   character(len=8),parameter :: fname(f_size) = &
      (/'    area',' hfreeze','   hmelt','  hnetsw','   hlwdn', &
        '   hlwup',' hlatvap',' hlatfus','  hiroff','    hsen', &
        ' wfreeze','   wmelt','   wrain','   wsnow', &
        '   wevap',' wrunoff',' wfrzrof' /)

   !--- P for period ---

   integer(in),parameter :: p_size = 5

   integer(in),parameter :: p_inst = 1
   integer(in),parameter :: p_day  = 2
   integer(in),parameter :: p_mon  = 3
   integer(in),parameter :: p_ann  = 4
   integer(in),parameter :: p_inf  = 5

   character(len=8),parameter :: pname(p_size) = &
      (/'    inst','   daily',' monthly','  annual','all_time' /)

! !PUBLIC DATA MEMBERS

   !--- time-averaged (annual?) global budge diagnostics ---
   !--- note: call sum0 then save budg_dataG and budg_ns on restart from/to root pe ---
   real(r8),public :: budg_dataL(f_size,c_size,p_size) ! local sum, valid on all pes
   real(r8),public :: budg_dataG(f_size,c_size,p_size) ! global sum, valid only on root pe
   real(r8),public :: budg_ns   (f_size,c_size,p_size) ! counter, valid only on root pe

   character(len=*),parameter :: afldname  = 'aream'
   character(len=*),parameter :: latname   = 'lat'
   character(len=*),parameter :: afracname = 'afrac'
   character(len=*),parameter :: lfracname = 'lfrac'
   character(len=*),parameter :: ofracname = 'ofrac'
   character(len=*),parameter :: ifracname = 'ifrac'

   character(*),parameter :: modName = "(seq_diag_mct) "

   integer(in),parameter :: debug = 0 ! internal debug level

! !PRIVATE DATA MEMBERS

   integer :: index_a2x_Faxa_swnet
   integer :: index_a2x_Faxa_lwdn
   integer :: index_a2x_Faxa_rainc
   integer :: index_a2x_Faxa_rainl
   integer :: index_a2x_Faxa_snowc
   integer :: index_a2x_Faxa_snowl

   integer :: index_x2a_Faxx_lwup
   integer :: index_x2a_Faxx_lat
   integer :: index_x2a_Faxx_sen
   integer :: index_x2a_Faxx_evap

   integer :: index_l2x_Fall_swnet
   integer :: index_l2x_Fall_lwup
   integer :: index_l2x_Fall_lat
   integer :: index_l2x_Fall_sen
   integer :: index_l2x_Fall_evap
   integer :: index_l2x_Flrl_rofl
   integer :: index_l2x_Flrl_rofi

   integer :: index_x2l_Faxa_lwdn
   integer :: index_x2l_Faxa_rainc
   integer :: index_x2l_Faxa_rainl
   integer :: index_x2l_Faxa_snowc
   integer :: index_x2l_Faxa_snowl
   integer :: index_x2l_Flrr_flood

   integer :: index_r2x_Forr_rofl
   integer :: index_r2x_Forr_rofi
   integer :: index_r2x_Firr_rofi
   integer :: index_r2x_Flrr_flood

   integer :: index_x2r_Flrl_rofl
   integer :: index_x2r_Flrl_rofi

   integer :: index_o2x_Fioo_q

   integer :: index_xao_Faox_lwup
   integer :: index_xao_Faox_lat
   integer :: index_xao_Faox_sen
   integer :: index_xao_Faox_evap

   integer :: index_x2o_Foxx_lwup
   integer :: index_x2o_Foxx_lat
   integer :: index_x2o_Foxx_sen
   integer :: index_x2o_Foxx_evap
   integer :: index_x2o_Foxx_swnet
   integer :: index_x2o_Foxx_rofl
   integer :: index_x2o_Foxx_rofi
   integer :: index_x2o_Faxa_lwdn
   integer :: index_x2o_Faxa_rain
   integer :: index_x2o_Faxa_snow
   integer :: index_x2o_Fioi_melth
   integer :: index_x2o_Fioi_meltw

   integer :: index_i2x_Fioi_melth
   integer :: index_i2x_Fioi_meltw
   integer :: index_i2x_Faii_swnet
   integer :: index_i2x_Fioi_swpen
   integer :: index_i2x_Faii_lwup
   integer :: index_i2x_Faii_lat
   integer :: index_i2x_Faii_sen
   integer :: index_i2x_Faii_evap

   integer :: index_x2i_Faxa_lwdn
   integer :: index_x2i_Faxa_rain
   integer :: index_x2i_Faxa_snow
   integer :: index_x2i_Fioo_q
   integer :: index_x2i_Fixx_rofi

   integer :: index_g2x_Fogg_rofl
   integer :: index_g2x_Fogg_rofi
   integer :: index_g2x_Figg_rofi

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_zero_mct - zero out global budget diagnostic data.
!
! !DESCRIPTION:
!    Zero out global budget diagnostic data.
!
! !REVISION HISTORY:
!    2008-jul-11 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_zero_mct(EClock,mode)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), intent(in),optional :: EClock
   character(len=*), intent(in),optional :: mode

!EOP

   integer(IN) :: ip,yr,mon,day,sec
   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_zero_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(EClock) .and. .not. present(mode)) then
      call shr_sys_abort(subName//' ERROR EClock or mode should be present')
   endif

   if (present(EClock)) then
      call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
         curr_mon=mon,curr_day=day,curr_tod=sec)

      do ip = 1,p_size
         if (ip == p_inst) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
         if (ip==p_day .and. sec==0) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
         if (ip==p_mon .and. day==1 .and. sec==0) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
         if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
      enddo
   endif

   if (present(mode)) then
      if (trim(mode) == 'inst') then
         budg_dataL(:,:,p_inst) = 0.0_r8
         budg_dataG(:,:,p_inst) = 0.0_r8
         budg_ns(:,:,p_inst) = 0.0_r8
      elseif (trim(mode) == 'day') then
         budg_dataL(:,:,p_day) = 0.0_r8
         budg_dataG(:,:,p_day) = 0.0_r8
         budg_ns(:,:,p_day) = 0.0_r8
      elseif (trim(mode) == 'mon') then
         budg_dataL(:,:,p_mon) = 0.0_r8
         budg_dataG(:,:,p_mon) = 0.0_r8
         budg_ns(:,:,p_mon) = 0.0_r8
      elseif (trim(mode) == 'ann') then
         budg_dataL(:,:,p_ann) = 0.0_r8
         budg_dataG(:,:,p_ann) = 0.0_r8
         budg_ns(:,:,p_ann) = 0.0_r8
      elseif (trim(mode) == 'inf') then
         budg_dataL(:,:,p_inf) = 0.0_r8
         budg_dataG(:,:,p_inf) = 0.0_r8
         budg_ns(:,:,p_inf) = 0.0_r8
      elseif (trim(mode) == 'all') then
         budg_dataL(:,:,:) = 0.0_r8
         budg_dataG(:,:,:) = 0.0_r8
         budg_ns(:,:,:) = 0.0_r8
      else
         call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
      endif
   endif

end subroutine seq_diag_zero_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_accum_mct - accum out global budget diagnostic data.
!
! !DESCRIPTION:
!    Accum out global budget diagnostic data.
!
! !REVISION HISTORY:
!    2008-jul-11 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_accum_mct()
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    integer(in) :: ip
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_accum_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    do ip = p_inst+1,p_size
! Subprogram not used       budg_dataL(:,:,ip) = budg_dataL(:,:,ip) + budg_dataL(:,:,p_inst)
! Subprogram not used    enddo
! Subprogram not used    budg_ns(:,:,:) = budg_ns(:,:,:) + 1.0_r8
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_accum_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_sum0_mct - sum local to global on root
!
! !DESCRIPTION:
!    Sum local values to global on root
!
! !REVISION HISTORY:
!    2008-jul-19 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_sum0_mct()
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    real(r8) :: budg_dataGtmp(f_size,c_size,p_size) ! temporary sum
! Subprogram not used    integer(in)      :: mpicom      ! mpi comm
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_sum0_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call seq_comm_setptrs(CPLID,mpicom=mpicom)
! Subprogram not used    budg_dataGtmp = 0.0_r8
! Subprogram not used    call shr_mpi_sum(budg_dataL,budg_dataGtmp,mpicom,subName)
! Subprogram not used    budg_dataG = budg_dataG + budg_dataGtmp
! Subprogram not used    budg_dataL = 0.0_r8
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_sum0_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_atm_mct - compute global atm input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global atm input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_atm_mct( atm, frac_a, do_a2x, do_x2a )
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(component_type), intent(in) :: atm    ! component type for instance1
! Subprogram not used    type(mct_aVect)     , intent(in) :: frac_a ! frac bundle
! Subprogram not used    logical, optional   , intent(in) :: do_a2x             
! Subprogram not used    logical, optional   , intent(in) :: do_x2a             
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    type(mct_aVect), pointer :: a2x_a        ! model to drv bundle
! Subprogram not used    type(mct_aVect), pointer :: x2a_a        ! drv to model bundle
! Subprogram not used    type(mct_ggrid), pointer :: dom_a
! Subprogram not used    integer(in)              :: k,n,ic,if,ip      ! generic index
! Subprogram not used    integer(in)              :: kArea             ! index of area field in aVect
! Subprogram not used    integer(in)              :: kLat              ! index of lat field in aVect
! Subprogram not used    integer(in)              :: kl,ka,ko,ki       ! fraction indices
! Subprogram not used    integer(in)              :: lSize             ! size of aVect
! Subprogram not used    real(r8)                 :: da,di,do,dl       ! area of a grid cell
! Subprogram not used    logical,save             :: first_time = .true.
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_atm_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dom_a => component_get_dom_cx(atm)
! Subprogram not used    a2x_a => component_get_c2x_cx(atm)  
! Subprogram not used    x2a_a => component_get_x2c_cx(atm)  
! Subprogram not used 
! Subprogram not used    kArea = mct_aVect_indexRA(dom_a%data,afldname)
! Subprogram not used    kLat  = mct_aVect_indexRA(dom_a%data,latname)
! Subprogram not used    ka    = mct_aVect_indexRA(frac_a,afracname)
! Subprogram not used    kl    = mct_aVect_indexRA(frac_a,lfracname)
! Subprogram not used    ko    = mct_aVect_indexRA(frac_a,ofracname)
! Subprogram not used    ki    = mct_aVect_indexRA(frac_a,ifracname)
! Subprogram not used 
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used    ! add values found in this bundle to the budget table
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used 
! Subprogram not used    if (present(do_a2x)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_a2x_Faxa_swnet  = mct_aVect_indexRA(a2x_a,'Faxa_swnet')
! Subprogram not used          index_a2x_Faxa_lwdn   = mct_aVect_indexRA(a2x_a,'Faxa_lwdn')
! Subprogram not used          index_a2x_Faxa_rainc  = mct_aVect_indexRA(a2x_a,'Faxa_rainc')
! Subprogram not used          index_a2x_Faxa_rainl  = mct_aVect_indexRA(a2x_a,'Faxa_rainl')
! Subprogram not used          index_a2x_Faxa_snowc  = mct_aVect_indexRA(a2x_a,'Faxa_snowc')
! Subprogram not used          index_a2x_Faxa_snowl  = mct_aVect_indexRA(a2x_a,'Faxa_snowl')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(a2x_a)
! Subprogram not used       do n=1,lSize
! Subprogram not used       do k=1,4
! Subprogram not used 
! Subprogram not used          if (k == 1) then
! Subprogram not used             ic = c_atm_ar
! Subprogram not used             da = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
! Subprogram not used          elseif (k == 2) then
! Subprogram not used             ic = c_lnd_ar
! Subprogram not used             da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
! Subprogram not used          elseif (k == 3) then
! Subprogram not used             ic = c_ocn_ar
! Subprogram not used             da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
! Subprogram not used          elseif (k == 4) then
! Subprogram not used             if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
! Subprogram not used                ic = c_inh_ar
! Subprogram not used             else
! Subprogram not used                ic = c_ish_ar
! Subprogram not used             endif
! Subprogram not used             da = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da
! Subprogram not used          if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_swnet,n)
! Subprogram not used          if = f_hlwdn ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_lwdn,n)
! Subprogram not used          if = f_wrain ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_rainc,n) &
! Subprogram not used                                                                     + da*a2x_a%rAttr(index_a2x_Faxa_rainl,n)
! Subprogram not used          if = f_wsnow ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_snowc,n) &
! Subprogram not used                                                                     + da*a2x_a%rAttr(index_a2x_Faxa_snowl,n)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used       ! --- heat implied by snow flux ---
! Subprogram not used       ic = c_atm_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       ic = c_lnd_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       ic = c_ocn_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       ic = c_inh_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       ic = c_ish_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (present(do_x2a)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_x2a_Faxx_lwup   = mct_aVect_indexRA(x2a_a,'Faxx_lwup')
! Subprogram not used          index_x2a_Faxx_lat    = mct_aVect_indexRA(x2a_a,'Faxx_lat')
! Subprogram not used          index_x2a_Faxx_sen    = mct_aVect_indexRA(x2a_a,'Faxx_sen')
! Subprogram not used          index_x2a_Faxx_evap   = mct_aVect_indexRA(x2a_a,'Faxx_evap')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(x2a_a)
! Subprogram not used       do n=1,lSize
! Subprogram not used       do k=1,4
! Subprogram not used 
! Subprogram not used          if (k == 1) then
! Subprogram not used             ic = c_atm_as
! Subprogram not used             da = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
! Subprogram not used          elseif (k == 2) then
! Subprogram not used             ic = c_lnd_as
! Subprogram not used             da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
! Subprogram not used          elseif (k == 3) then
! Subprogram not used             ic = c_ocn_as
! Subprogram not used             da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
! Subprogram not used          elseif (k == 4) then
! Subprogram not used             if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
! Subprogram not used                ic = c_inh_as
! Subprogram not used             else
! Subprogram not used                ic = c_ish_as
! Subprogram not used             endif
! Subprogram not used             da = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if = f_area ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da
! Subprogram not used          if = f_hlwup; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_lwup,n)
! Subprogram not used          if = f_hlatv; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_lat,n)
! Subprogram not used          if = f_hsen ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_sen,n)
! Subprogram not used          if = f_wevap; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_evap,n)
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    first_time = .false.
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_atm_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_lnd_mct - compute global lnd input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global lnd input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_lnd_mct( lnd, frac_l, do_l2x, do_x2l)
! Subprogram not used 
! Subprogram not used    type(component_type), intent(in) :: lnd    ! component type for instance1
! Subprogram not used    type(mct_aVect)     , intent(in) :: frac_l ! frac bundle
! Subprogram not used    logical, optional   , intent(in) :: do_l2x             
! Subprogram not used    logical, optional   , intent(in) :: do_x2l             
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    type(mct_aVect), pointer :: l2x_l        ! model to drv bundle
! Subprogram not used    type(mct_aVect), pointer :: x2l_l        ! drv to model bundle
! Subprogram not used    type(mct_ggrid), pointer :: dom_l
! Subprogram not used    integer(in)              :: k,n,ic,if,ip ! generic index
! Subprogram not used    integer(in)              :: kArea        ! index of area field in aVect
! Subprogram not used    integer(in)              :: kLat         ! index of lat field in aVect
! Subprogram not used    integer(in)              :: kl,ka,ko,ki  ! fraction indices
! Subprogram not used    integer(in)              :: lSize        ! size of aVect
! Subprogram not used    real(r8)                 :: da,di,do,dl  ! area of a grid cell
! Subprogram not used    logical,save             :: first_time = .true.
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_lnd_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used    ! add values found in this bundle to the budget table
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dom_l => component_get_dom_cx(lnd)
! Subprogram not used    l2x_l => component_get_c2x_cx(lnd)  
! Subprogram not used    x2l_l => component_get_x2c_cx(lnd)  
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used 
! Subprogram not used    kArea = mct_aVect_indexRA(dom_l%data,afldname)
! Subprogram not used    kl    = mct_aVect_indexRA(frac_l,lfracname)
! Subprogram not used 
! Subprogram not used    if (present(do_l2x)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_l2x_Fall_swnet  = mct_aVect_indexRA(l2x_l,'Fall_swnet')
! Subprogram not used          index_l2x_Fall_lwup   = mct_aVect_indexRA(l2x_l,'Fall_lwup')
! Subprogram not used          index_l2x_Fall_lat    = mct_aVect_indexRA(l2x_l,'Fall_lat')
! Subprogram not used          index_l2x_Fall_sen    = mct_aVect_indexRA(l2x_l,'Fall_sen')
! Subprogram not used          index_l2x_Fall_evap   = mct_aVect_indexRA(l2x_l,'Fall_evap')
! Subprogram not used          index_l2x_Flrl_rofl   = mct_aVect_indexRA(l2x_l,'Flrl_rofl')
! Subprogram not used          index_l2x_Flrl_rofi   = mct_aVect_indexRA(l2x_l,'Flrl_rofi')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(l2x_l)
! Subprogram not used       ic = c_lnd_lr
! Subprogram not used       do n=1,lSize
! Subprogram not used          dl =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
! Subprogram not used          if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl
! Subprogram not used          if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_swnet,n)
! Subprogram not used          if = f_hlwup ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_lwup,n)
! Subprogram not used          if = f_hlatv ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_lat,n)
! Subprogram not used          if = f_hsen  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_sen,n)
! Subprogram not used          if = f_wevap ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_evap,n)
! Subprogram not used          if = f_wroff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dl*l2x_l%rAttr(index_l2x_Flrl_rofl,n)
! Subprogram not used          if = f_wioff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dl*l2x_l%rAttr(index_l2x_Flrl_rofi,n)
! Subprogram not used       end do
! Subprogram not used       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (present(do_x2l)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_x2l_Faxa_lwdn   = mct_aVect_indexRA(x2l_l,'Faxa_lwdn')
! Subprogram not used          index_x2l_Faxa_rainc  = mct_aVect_indexRA(x2l_l,'Faxa_rainc')
! Subprogram not used          index_x2l_Faxa_rainl  = mct_aVect_indexRA(x2l_l,'Faxa_rainl')
! Subprogram not used          index_x2l_Faxa_snowc  = mct_aVect_indexRA(x2l_l,'Faxa_snowc')
! Subprogram not used          index_x2l_Faxa_snowl  = mct_aVect_indexRA(x2l_l,'Faxa_snowl')
! Subprogram not used          index_x2l_Flrr_flood  = mct_aVect_indexRA(x2l_l,'Flrr_flood')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(x2l_l)
! Subprogram not used       ic = c_lnd_ls
! Subprogram not used       do n=1,lSize
! Subprogram not used          dl =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
! Subprogram not used          if = f_area ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl
! Subprogram not used          if = f_hlwdn; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*x2l_l%rAttr(index_x2l_Faxa_lwdn,n)
! Subprogram not used          if = f_wrain; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*x2l_l%rAttr(index_x2l_Faxa_rainc,n) &
! Subprogram not used                                                                    + dl*x2l_l%rAttr(index_x2l_Faxa_rainl,n)
! Subprogram not used          if = f_wsnow; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*x2l_l%rAttr(index_x2l_Faxa_snowc,n) &
! Subprogram not used                                                                    + dl*x2l_l%rAttr(index_x2l_Faxa_snowl,n)
! Subprogram not used          if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dl*x2l_l%rAttr(index_x2l_Flrr_flood,n)
! Subprogram not used       end do
! Subprogram not used       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    first_time = .false.
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_lnd_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_rof_mct - compute global rof input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global rof input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_rof_mct( rof, frac_r)
! Subprogram not used 
! Subprogram not used    type(component_type), intent(in) :: rof    ! component type for instance1
! Subprogram not used    type(mct_aVect)     , intent(in) :: frac_r ! frac bundle
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    type(mct_aVect), pointer :: r2x_r
! Subprogram not used    type(mct_aVect), pointer :: x2r_r
! Subprogram not used    type(mct_ggrid), pointer :: dom_r
! Subprogram not used    integer(in)              :: k,n,ic,if,ip      ! generic index
! Subprogram not used    integer(in)              :: kArea             ! index of area field in aVect
! Subprogram not used    integer(in)              :: kLat              ! index of lat field in aVect
! Subprogram not used    integer(in)              :: kl,ka,ko,ki,kr    ! fraction indices
! Subprogram not used    integer(in)              :: lSize             ! size of aVect
! Subprogram not used    real(r8)                 :: da,di,do,dl,dr    ! area of a grid cell
! Subprogram not used    logical,save             :: first_time = .true.
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_rof_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used    ! add values found in this bundle to the budget table
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dom_r => component_get_dom_cx(rof)
! Subprogram not used    r2x_r => component_get_c2x_cx(rof)  
! Subprogram not used    x2r_r => component_get_x2c_cx(rof)  
! Subprogram not used 
! Subprogram not used    if (first_time) then
! Subprogram not used       index_x2r_Flrl_rofl  = mct_aVect_indexRA(x2r_r,'Flrl_rofl')
! Subprogram not used       index_x2r_Flrl_rofi  = mct_aVect_indexRA(x2r_r,'Flrl_rofi')
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used    ic = c_rof_rr
! Subprogram not used    kArea = mct_aVect_indexRA(dom_r%data,afldname)
! Subprogram not used    lSize = mct_avect_lSize(x2r_r)
! Subprogram not used    do n=1,lSize
! Subprogram not used       dr =  dom_r%data%rAttr(kArea,n)
! Subprogram not used       if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dr*x2r_r%rAttr(index_x2r_Flrl_rofl,n)
! Subprogram not used       if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dr*x2r_r%rAttr(index_x2r_Flrl_rofi,n)
! Subprogram not used    end do
! Subprogram not used    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used 
! Subprogram not used    if (first_time) then
! Subprogram not used       index_r2x_Forr_rofl   = mct_aVect_indexRA(r2x_r,'Forr_rofl')
! Subprogram not used       index_r2x_Forr_rofi   = mct_aVect_indexRA(r2x_r,'Forr_rofi')
! Subprogram not used       index_r2x_Firr_rofi   = mct_aVect_indexRA(r2x_r,'Firr_rofi')
! Subprogram not used       index_r2x_Flrr_flood  = mct_aVect_indexRA(r2x_r,'Flrr_flood')
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used    ic = c_rof_rs
! Subprogram not used    kArea = mct_aVect_indexRA(dom_r%data,afldname)
! Subprogram not used    lSize = mct_avect_lSize(r2x_r)
! Subprogram not used    do n=1,lSize
! Subprogram not used       dr =  dom_r%data%rAttr(kArea,n)
! Subprogram not used       if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dr*r2x_r%rAttr(index_r2x_Forr_rofl,n) &
! Subprogram not used                                                                 + dr*r2x_r%rAttr(index_r2x_Flrr_flood,n)
! Subprogram not used       if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dr*r2x_r%rAttr(index_r2x_Forr_rofi,n) &
! Subprogram not used                                                                 - dr*r2x_r%rAttr(index_r2x_Firr_rofi,n)
! Subprogram not used    end do
! Subprogram not used    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used 
! Subprogram not used    first_time = .false.
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_rof_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_glc_mct - compute global glc input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global glc input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_glc_mct( glc, frac_g)
! Subprogram not used 
! Subprogram not used    type(component_type), intent(in) :: glc    ! component type for instance1
! Subprogram not used    type(mct_aVect)     , intent(in) :: frac_g ! frac bundle
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    type(mct_aVect), pointer :: g2x_g
! Subprogram not used    type(mct_aVect), pointer :: x2g_g
! Subprogram not used    type(mct_ggrid), pointer :: dom_g
! Subprogram not used    integer(in)              :: k,n,ic,if,ip      ! generic index
! Subprogram not used    integer(in)              :: kArea             ! index of area field in aVect
! Subprogram not used    integer(in)              :: kLat              ! index of lat field in aVect
! Subprogram not used    integer(in)              :: kl,ka,ko,ki,kr,kg ! fraction indices
! Subprogram not used    integer(in)              :: lSize             ! size of aVect
! Subprogram not used    real(r8)                 :: da,di,do,dl,dr,dg ! area of a grid cell
! Subprogram not used    logical,save             :: first_time = .true.
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_glc_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used    ! add values found in this bundle to the budget table
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dom_g => component_get_dom_cx(glc)
! Subprogram not used    g2x_g => component_get_c2x_cx(glc)  
! Subprogram not used    x2g_g => component_get_x2c_cx(glc)  
! Subprogram not used 
! Subprogram not used    if (first_time) then
! Subprogram not used       index_g2x_Fogg_rofl   = mct_aVect_indexRA(g2x_g,'Fogg_rofl')
! Subprogram not used       index_g2x_Fogg_rofi   = mct_aVect_indexRA(g2x_g,'Fogg_rofi')
! Subprogram not used       index_g2x_Figg_rofi   = mct_aVect_indexRA(g2x_g,'Figg_rofi')
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used    ic = c_glc_gs
! Subprogram not used    kArea = mct_aVect_indexRA(dom_g%data,afldname)
! Subprogram not used    lSize = mct_avect_lSize(g2x_g)
! Subprogram not used    do n=1,lSize
! Subprogram not used       dg =  dom_g%data%rAttr(kArea,n)
! Subprogram not used       if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dg*g2x_g%rAttr(index_g2x_Fogg_rofl,n)
! Subprogram not used       if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dg*g2x_g%rAttr(index_g2x_Fogg_rofi,n) &
! Subprogram not used                                                                 - dg*g2x_g%rAttr(index_g2x_Figg_rofi,n)
! Subprogram not used    end do
! Subprogram not used    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used 
! Subprogram not used    first_time = .false.
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_glc_mct

!BOP ===========================================================================
!
! !IROUTINE: seq_diag_ocn_mct - compute global ocn input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global ocn input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_ocn_mct( ocn, xao_o, frac_o, do_o2x, do_x2o, do_xao)
! Subprogram not used 
! Subprogram not used    type(component_type) , intent(in)          :: ocn    ! component type for instance1
! Subprogram not used    type(mct_aVect)      , intent(in)          :: frac_o ! frac bundle
! Subprogram not used    type(mct_aVect)      , intent(in)          :: xao_o  
! Subprogram not used    logical              , intent(in),optional :: do_o2x
! Subprogram not used    logical              , intent(in),optional :: do_x2o
! Subprogram not used    logical              , intent(in),optional :: do_xao
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    type(mct_aVect), pointer :: o2x_o        ! model to drv bundle
! Subprogram not used    type(mct_aVect), pointer :: x2o_o        ! drv to model bundle
! Subprogram not used    type(mct_ggrid), pointer :: dom_o
! Subprogram not used    integer(in)              :: k,n,if,ic,ip ! generic index
! Subprogram not used    integer(in)              :: kArea        ! index of area field in aVect
! Subprogram not used    integer(in)              :: kLat         ! index of lat field in aVect
! Subprogram not used    integer(in)              :: kl,ka,ko,ki  ! fraction indices
! Subprogram not used    integer(in)              :: lSize        ! size of aVect
! Subprogram not used    real(r8)                 :: da,di,do,dl  ! area of a grid cell
! Subprogram not used    logical,save             :: first_time = .true.
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_ocn_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (.not. present(do_o2x) .and. &
! Subprogram not used        .not. present(do_x2o) .and. &
! Subprogram not used        .not. present(do_xao)) then
! Subprogram not used       call shr_sys_abort(subName//"ERROR: must input a bundle")
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used    ! add values found in this bundle to the budget table
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dom_o => component_get_dom_cx(ocn)
! Subprogram not used    o2x_o => component_get_c2x_cx(ocn)  
! Subprogram not used    x2o_o => component_get_x2c_cx(ocn)  
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used 
! Subprogram not used    kArea = mct_aVect_indexRA(dom_o%data,afldname)
! Subprogram not used    ko    = mct_aVect_indexRA(frac_o,ofracname)
! Subprogram not used    ki    = mct_aVect_indexRA(frac_o,ifracname)
! Subprogram not used 
! Subprogram not used    if (present(do_o2x)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_o2x_Fioo_q      = mct_aVect_indexRA(o2x_o,'Fioo_q')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(o2x_o)
! Subprogram not used       ic = c_ocn_or
! Subprogram not used       do n=1,lSize
! Subprogram not used          do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
! Subprogram not used          di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
! Subprogram not used          if = f_area; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do
! Subprogram not used          if = f_hfrz; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*max(0.0_r8,o2x_o%rAttr(index_o2x_Fioo_q,n))
! Subprogram not used       end do
! Subprogram not used       budg_dataL(f_wfrz,ic,ip) = budg_dataL(f_hfrz,ic,ip) * HFLXtoWFLX
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (present(do_xao)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_xao_Faox_lwup   = mct_aVect_indexRA(xao_o,'Faox_lwup') 
! Subprogram not used          index_xao_Faox_lat    = mct_aVect_indexRA(xao_o,'Faox_lat')  
! Subprogram not used          index_xao_Faox_sen    = mct_aVect_indexRA(xao_o,'Faox_sen') 
! Subprogram not used          index_xao_Faox_evap   = mct_aVect_indexRA(xao_o,'Faox_evap')  
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(xao_o)
! Subprogram not used       ic = c_ocn_or
! Subprogram not used       do n=1,lSize
! Subprogram not used          do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
! Subprogram not used          if = f_hlwup; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_lwup,n)
! Subprogram not used          if = f_hlatv; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_lat,n)
! Subprogram not used          if = f_hsen ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_sen,n)
! Subprogram not used          if = f_wevap; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_evap,n)
! Subprogram not used       end do
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (present(do_x2o)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_x2o_Fioi_melth  = mct_aVect_indexRA(x2o_o,'Fioi_melth')  
! Subprogram not used          index_x2o_Fioi_meltw  = mct_aVect_indexRA(x2o_o,'Fioi_meltw') 
! Subprogram not used          index_x2o_Foxx_swnet  = mct_aVect_indexRA(x2o_o,'Foxx_swnet')
! Subprogram not used          index_x2o_Faxa_lwdn   = mct_aVect_indexRA(x2o_o,'Faxa_lwdn')
! Subprogram not used          index_x2o_Faxa_rain   = mct_aVect_indexRA(x2o_o,'Faxa_rain') 
! Subprogram not used          index_x2o_Faxa_snow   = mct_aVect_indexRA(x2o_o,'Faxa_snow')  
! Subprogram not used          index_x2o_Foxx_lwup   = mct_aVect_indexRA(x2o_o,'Foxx_lwup') 
! Subprogram not used          index_x2o_Foxx_lat    = mct_aVect_indexRA(x2o_o,'Foxx_lat')  
! Subprogram not used          index_x2o_Foxx_sen    = mct_aVect_indexRA(x2o_o,'Foxx_sen') 
! Subprogram not used          index_x2o_Foxx_evap   = mct_aVect_indexRA(x2o_o,'Foxx_evap')  
! Subprogram not used          index_x2o_Foxx_rofl   = mct_aVect_indexRA(x2o_o,'Foxx_rofl')
! Subprogram not used          index_x2o_Foxx_rofi   = mct_aVect_indexRA(x2o_o,'Foxx_rofi')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       if (.not. present(do_xao)) then
! Subprogram not used          ! these are in x2o but they really are the atm/ocean flux 
! Subprogram not used          ! computed in the coupler and are "like" an o2x
! Subprogram not used          lSize = mct_avect_lSize(x2o_o)
! Subprogram not used          ic = c_ocn_or
! Subprogram not used          do n=1,lSize
! Subprogram not used             do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
! Subprogram not used             di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
! Subprogram not used             if = f_hlwup; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_lwup,n)
! Subprogram not used             if = f_hlatv; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_lat,n)
! Subprogram not used             if = f_hsen ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_sen,n)
! Subprogram not used             if = f_wevap; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_evap,n)
! Subprogram not used          end do
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(x2o_o)
! Subprogram not used       ic = c_ocn_os
! Subprogram not used       do n=1,lSize
! Subprogram not used          do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
! Subprogram not used          di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
! Subprogram not used          if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do
! Subprogram not used          if = f_wmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Fioi_meltw,n)
! Subprogram not used          if = f_hmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Fioi_melth,n)
! Subprogram not used          if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_swnet,n)
! Subprogram not used          if = f_hlwdn ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Faxa_lwdn,n)
! Subprogram not used          if = f_wrain ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Faxa_rain,n)
! Subprogram not used          if = f_wsnow ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Faxa_snow,n)
! Subprogram not used          if = f_wroff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_rofl,n)
! Subprogram not used          if = f_wioff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_rofi,n)
! Subprogram not used       end do
! Subprogram not used       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    first_time = .false.
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_ocn_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_ice_mct - compute global ice input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global ice input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_diag_ice_mct( ice, frac_i, do_i2x, do_x2i)
! Subprogram not used 
! Subprogram not used    type(component_type), intent(in)           :: ice    ! component type for instance1
! Subprogram not used    type(mct_aVect)     , intent(in)           :: frac_i ! frac bundle
! Subprogram not used    logical             , intent(in), optional :: do_i2x
! Subprogram not used    logical             , intent(in), optional :: do_x2i
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    type(mct_aVect), pointer :: i2x_i        ! model to drv bundle
! Subprogram not used    type(mct_aVect), pointer :: x2i_i        ! drv to model bundle
! Subprogram not used    type(mct_ggrid), pointer :: dom_i
! Subprogram not used    integer(in)              :: k,n,ic,if,ip ! generic index
! Subprogram not used    integer(in)              :: kArea        ! index of area field in aVect
! Subprogram not used    integer(in)              :: kLat         ! index of lat field in aVect
! Subprogram not used    integer(in)              :: kl,ka,ko,ki  ! fraction indices
! Subprogram not used    integer(in)              :: lSize        ! size of aVect
! Subprogram not used    real(r8)                 :: da,di,do,dl  ! area of a grid cell
! Subprogram not used    logical,save             :: first_time = .true.
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_ice_mct) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used    ! add values found in this bundle to the budget table
! Subprogram not used    !---------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    dom_i => component_get_dom_cx(ice)
! Subprogram not used    i2x_i => component_get_c2x_cx(ice)  
! Subprogram not used    x2i_i => component_get_x2c_cx(ice)  
! Subprogram not used 
! Subprogram not used    ip = p_inst
! Subprogram not used 
! Subprogram not used    kArea = mct_aVect_indexRA(dom_i%data,afldname)
! Subprogram not used    kLat  = mct_aVect_indexRA(dom_i%data,latname)
! Subprogram not used    ki    = mct_aVect_indexRA(frac_i,ifracname)
! Subprogram not used    ko    = mct_aVect_indexRA(frac_i,ofracname)
! Subprogram not used 
! Subprogram not used    if (present(do_i2x)) then
! Subprogram not used          index_i2x_Fioi_melth  = mct_aVect_indexRA(i2x_i,'Fioi_melth')
! Subprogram not used          index_i2x_Fioi_meltw  = mct_aVect_indexRA(i2x_i,'Fioi_meltw')
! Subprogram not used          index_i2x_Fioi_swpen  = mct_aVect_indexRA(i2x_i,'Fioi_swpen')
! Subprogram not used          index_i2x_Faii_swnet  = mct_aVect_indexRA(i2x_i,'Faii_swnet')
! Subprogram not used          index_i2x_Faii_lwup   = mct_aVect_indexRA(i2x_i,'Faii_lwup')
! Subprogram not used          index_i2x_Faii_lat    = mct_aVect_indexRA(i2x_i,'Faii_lat')
! Subprogram not used          index_i2x_Faii_sen    = mct_aVect_indexRA(i2x_i,'Faii_sen')
! Subprogram not used          index_i2x_Faii_evap   = mct_aVect_indexRA(i2x_i,'Faii_evap')
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(i2x_i)
! Subprogram not used       do n=1,lSize
! Subprogram not used          if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
! Subprogram not used             ic = c_inh_ir
! Subprogram not used          else
! Subprogram not used             ic = c_ish_ir
! Subprogram not used          endif
! Subprogram not used          do =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
! Subprogram not used          di =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
! Subprogram not used          if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di
! Subprogram not used          if = f_hmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - di*i2x_i%rAttr(index_i2x_Fioi_melth,n)
! Subprogram not used          if = f_wmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - di*i2x_i%rAttr(index_i2x_Fioi_meltw,n)
! Subprogram not used          if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_swnet,n) &
! Subprogram not used                                                                     - di*i2x_i%rAttr(index_i2x_Fioi_swpen,n)
! Subprogram not used          if = f_hlwup ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_lwup,n)
! Subprogram not used          if = f_hlatv ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_lat,n)
! Subprogram not used          if = f_hsen  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_sen,n)
! Subprogram not used          if = f_wevap ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_evap,n)
! Subprogram not used       end do
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (present(do_x2i)) then
! Subprogram not used       if (first_time) then
! Subprogram not used          index_x2i_Faxa_lwdn   = mct_aVect_indexRA(x2i_i,'Faxa_lwdn') 
! Subprogram not used          index_x2i_Faxa_rain   = mct_aVect_indexRA(x2i_i,'Faxa_rain')  
! Subprogram not used          index_x2i_Faxa_snow   = mct_aVect_indexRA(x2i_i,'Faxa_snow')  
! Subprogram not used          index_x2i_Fioo_q      = mct_aVect_indexRA(x2i_i,'Fioo_q')  
! Subprogram not used          index_x2i_Fixx_rofi   = mct_aVect_indexRA(x2i_i,'Fixx_rofi')
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       lSize = mct_avect_lSize(x2i_i)
! Subprogram not used       do n=1,lSize
! Subprogram not used          if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
! Subprogram not used             ic = c_inh_is
! Subprogram not used          else
! Subprogram not used             ic = c_ish_is
! Subprogram not used          endif
! Subprogram not used          do =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
! Subprogram not used          di =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
! Subprogram not used          if = f_area ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di
! Subprogram not used          if = f_hlwdn; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Faxa_lwdn,n)
! Subprogram not used          if = f_wrain; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Faxa_rain,n)
! Subprogram not used          if = f_wsnow; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Faxa_snow,n)
! Subprogram not used          if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Fixx_rofi,n)
! Subprogram not used          if = f_hfrz ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - (do+di)*max(0.0_r8,x2i_i%rAttr(index_x2i_Fioo_q,n))
! Subprogram not used       end do
! Subprogram not used       ic = c_inh_is  
! Subprogram not used       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used       budg_dataL(f_wfrz ,ic,ip) =  budg_dataL(f_hfrz ,ic,ip)*HFLXtoWFLX
! Subprogram not used       ic = c_ish_is
! Subprogram not used       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
! Subprogram not used       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
! Subprogram not used       budg_dataL(f_wfrz ,ic,ip) =  budg_dataL(f_hfrz ,ic,ip)*HFLXtoWFLX
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    first_time = .false.
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_ice_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_print_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global budget diagnostics.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used SUBROUTINE seq_diag_print_mct(EClock, stop_alarm, &
! Subprogram not used      budg_print_inst,  budg_print_daily,  budg_print_month,  &
! Subprogram not used      budg_print_ann,  budg_print_ltann,  budg_print_ltend)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(ESMF_Clock) , intent(in) :: EClock
! Subprogram not used    logical          , intent(in) :: stop_alarm
! Subprogram not used    integer          , intent(in) :: budg_print_inst
! Subprogram not used    integer          , intent(in) :: budg_print_daily
! Subprogram not used    integer          , intent(in) :: budg_print_month
! Subprogram not used    integer          , intent(in) :: budg_print_ann
! Subprogram not used    integer          , intent(in) :: budg_print_ltann
! Subprogram not used    integer          , intent(in) :: budg_print_ltend
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    integer(in)      :: ic,if,ip    ! data array indicies
! Subprogram not used    integer(in)      :: ica,icl,icn,ics,ico
! Subprogram not used    integer(in)      :: icar,icxs,icxr,icas
! Subprogram not used    integer(in)      :: n           ! loop counter
! Subprogram not used    integer(in)      :: nday        ! number of days in time avg
! Subprogram not used    integer(in)      :: cdate,sec   ! coded date, seconds
! Subprogram not used    integer(in)      :: yr,mon,day  ! date
! Subprogram not used    integer(in)      :: iam         ! pe number
! Subprogram not used    integer(in)      :: plev        ! print level
! Subprogram not used    logical          :: sumdone     ! has a sum been computed yet
! Subprogram not used    character(len=40):: str         ! string
! Subprogram not used    real(r8) :: dataGpr (f_size,c_size,p_size) ! values to print, scaled and such
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_print_mct) '
! Subprogram not used    character(*),parameter :: F00   = "('(seq_diag_print_mct) ',4a)"
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: FAH="(4a,i9,i6)"
! Subprogram not used    character(*),parameter :: FA0= "('    ',8x,6(6x,a8,1x))"
! Subprogram not used    character(*),parameter :: FA1= "('    ',a8,6f15.8)"
! Subprogram not used    character(*),parameter :: FA0r="('    ',8x,8(6x,a8,1x))"
! Subprogram not used    character(*),parameter :: FA1r="('    ',a8,8f15.8)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! print instantaneous budget data
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    sumdone = .false.
! Subprogram not used    call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used    call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
! Subprogram not used       curr_mon=mon,curr_day=day,curr_tod=sec)
! Subprogram not used    cdate = yr*10000+mon*100+day
! Subprogram not used 
! Subprogram not used    do ip = 1,p_size
! Subprogram not used       plev = 0
! Subprogram not used       if (ip == p_inst) then
! Subprogram not used          plev = max(plev,budg_print_inst)
! Subprogram not used       endif
! Subprogram not used       if (ip==p_day .and. sec==0) then
! Subprogram not used          plev = max(plev,budg_print_daily)
! Subprogram not used       endif
! Subprogram not used       if (ip==p_mon .and. day==1 .and. sec==0) then
! Subprogram not used          plev = max(plev,budg_print_month)
! Subprogram not used       endif
! Subprogram not used       if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
! Subprogram not used          plev = max(plev,budg_print_ann)
! Subprogram not used       endif
! Subprogram not used       if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
! Subprogram not used          plev = max(plev,budg_print_ltann)
! Subprogram not used       endif
! Subprogram not used       if (ip==p_inf .and. stop_alarm) then
! Subprogram not used          plev = max(plev,budg_print_ltend)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    if (plev > 0) then
! Subprogram not used ! ---- doprint ---- doprint ---- doprint ----
! Subprogram not used 
! Subprogram not used    if (.not.sumdone) then
! Subprogram not used       call seq_diag_sum0_mct()
! Subprogram not used       dataGpr = budg_dataG
! Subprogram not used       sumdone = .true.
! Subprogram not used 
! Subprogram not used    !  old budget normalizations (global area and 1e6 for water)
! Subprogram not used       dataGpr = dataGpr/(4.0_r8*shr_const_pi)
! Subprogram not used       dataGpr(f_w:f_size,:,:) = dataGpr(f_w:f_size,:,:) * 1.0e6_r8
! Subprogram not used       dataGpr = dataGpr/budg_ns
! Subprogram not used 
! Subprogram not used       if (iam /= 0) return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    ! ---------------------------------------------------------
! Subprogram not used    ! ---- detail atm budgets and breakdown into components ---
! Subprogram not used    ! ---------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (plev >= 3) then
! Subprogram not used    do ic = 1,2
! Subprogram not used       if (ic == 1) then
! Subprogram not used          ica = c_atm_ar
! Subprogram not used          icl = c_lnd_ar
! Subprogram not used          icn = c_inh_ar
! Subprogram not used          ics = c_ish_ar
! Subprogram not used          ico = c_ocn_ar
! Subprogram not used          str = "ATM_to_CPL"
! Subprogram not used       elseif (ic == 2) then
! Subprogram not used          ica = c_atm_as
! Subprogram not used          icl = c_lnd_as
! Subprogram not used          icn = c_inh_as
! Subprogram not used          ics = c_ish_as
! Subprogram not used          ico = c_ocn_as
! Subprogram not used          str = "CPL_TO_ATM"
! Subprogram not used       else
! Subprogram not used          call shr_sys_abort(subname//' ERROR in ic index code 411')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,trim(str)//' AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
! Subprogram not used       do if = f_a, f_h-1
! Subprogram not used          write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
! Subprogram not used                    dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
! Subprogram not used                                          dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
! Subprogram not used                    dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
! Subprogram not used       do if = f_h, f_w-1
! Subprogram not used          write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
! Subprogram not used                    dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
! Subprogram not used                                          dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
! Subprogram not used                    dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
! Subprogram not used       enddo
! Subprogram not used       write(logunit,FA1)    '   *SUM*',sum(dataGpr(f_h:f_w-1,ica,ip)),sum(dataGpr(f_h:f_w-1,icl,ip)), &
! Subprogram not used          sum(dataGpr(f_h:f_w-1,icn,ip)),sum(dataGpr(f_h:f_w-1,ics,ip)),sum(dataGpr(f_h:f_w-1,ico,ip)), &
! Subprogram not used                                        sum(dataGpr(f_h:f_w-1,ica,ip))+sum(dataGpr(f_h:f_w-1,icl,ip))+ &
! Subprogram not used          sum(dataGpr(f_h:f_w-1,icn,ip))+sum(dataGpr(f_h:f_w-1,ics,ip))+sum(dataGpr(f_h:f_w-1,ico,ip)) 
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
! Subprogram not used       do if = f_w, f_size
! Subprogram not used          write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
! Subprogram not used                    dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
! Subprogram not used                                          dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
! Subprogram not used                    dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
! Subprogram not used       enddo
! Subprogram not used       write(logunit,FA1)    '   *SUM*',sum(dataGpr(f_w:f_size,ica,ip)),sum(dataGpr(f_w:f_size,icl,ip)), &
! Subprogram not used          sum(dataGpr(f_w:f_size,icn,ip)),sum(dataGpr(f_w:f_size,ics,ip)),sum(dataGpr(f_w:f_size,ico,ip)), &
! Subprogram not used                                        sum(dataGpr(f_w:f_size,ica,ip))+sum(dataGpr(f_w:f_size,icl,ip))+ &
! Subprogram not used          sum(dataGpr(f_w:f_size,icn,ip))+sum(dataGpr(f_w:f_size,ics,ip))+sum(dataGpr(f_w:f_size,ico,ip)) 
! Subprogram not used    enddo
! Subprogram not used    endif   ! plev
! Subprogram not used 
! Subprogram not used    ! ---------------------------------------------------------
! Subprogram not used    ! ---- detail lnd/ocn/ice component budgets ----
! Subprogram not used    ! ---------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (plev >= 2) then
! Subprogram not used    do ic = 1,4
! Subprogram not used       if (ic == 1) then
! Subprogram not used          icar = c_lnd_ar
! Subprogram not used          icxs = c_lnd_ls
! Subprogram not used          icxr = c_lnd_lr
! Subprogram not used          icas = c_lnd_as
! Subprogram not used          str = "LND"
! Subprogram not used       elseif (ic == 2) then
! Subprogram not used          icar = c_ocn_ar
! Subprogram not used          icxs = c_ocn_os
! Subprogram not used          icxr = c_ocn_or
! Subprogram not used          icas = c_ocn_as
! Subprogram not used          str = "OCN"
! Subprogram not used       elseif (ic == 3) then
! Subprogram not used          icar = c_inh_ar
! Subprogram not used          icxs = c_inh_is
! Subprogram not used          icxr = c_inh_ir
! Subprogram not used          icas = c_inh_as
! Subprogram not used          str = "ICE_NH"
! Subprogram not used       elseif (ic == 4) then
! Subprogram not used          icar = c_ish_ar
! Subprogram not used          icxs = c_ish_is
! Subprogram not used          icxr = c_ish_ir
! Subprogram not used          icas = c_ish_as
! Subprogram not used          str = "ICE_SH"
! Subprogram not used       else
! Subprogram not used          call shr_sys_abort(subname//' ERROR in ic index code 412')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
! Subprogram not used       do if = f_h, f_w-1
! Subprogram not used          write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
! Subprogram not used                                           dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
! Subprogram not used                                          -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
! Subprogram not used                                           dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
! Subprogram not used       enddo
! Subprogram not used       write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_h:f_w-1,icar,ip)),sum(dataGpr(f_h:f_w-1,icxs,ip)), &
! Subprogram not used                                        sum(dataGpr(f_h:f_w-1,icxr,ip)),-sum(dataGpr(f_h:f_w-1,icas,ip)), &
! Subprogram not used                                        -sum(dataGpr(f_h:f_w-1,icar,ip))+sum(dataGpr(f_h:f_w-1,icxs,ip))+ &
! Subprogram not used                                        sum(dataGpr(f_h:f_w-1,icxr,ip))-sum(dataGpr(f_h:f_w-1,icas,ip))
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
! Subprogram not used       do if = f_w, f_size
! Subprogram not used          write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
! Subprogram not used                                          dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
! Subprogram not used                                          -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
! Subprogram not used                                          dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
! Subprogram not used       enddo
! Subprogram not used       write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_w:f_size,icar,ip)),sum(dataGpr(f_w:f_size,icxs,ip)), &
! Subprogram not used                                        sum(dataGpr(f_w:f_size,icxr,ip)),-sum(dataGpr(f_w:f_size,icas,ip)), &
! Subprogram not used                                        -sum(dataGpr(f_w:f_size,icar,ip))+sum(dataGpr(f_w:f_size,icxs,ip))+ &
! Subprogram not used                                        sum(dataGpr(f_w:f_size,icxr,ip))-sum(dataGpr(f_w:f_size,icas,ip))
! Subprogram not used 
! Subprogram not used    enddo
! Subprogram not used    endif   ! plev
! Subprogram not used 
! Subprogram not used    ! ---------------------------------------------------------
! Subprogram not used    ! ---- net summary budgets ----
! Subprogram not used    ! ---------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (plev >= 1) then
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,'NET AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0) '     atm','     lnd','     ocn','  ice nh','  ice sh',' *SUM*  '
! Subprogram not used       do if = 1,f_h-1
! Subprogram not used          write(logunit,FA1)    fname(if),dataGpr(if,c_atm_ar,ip), &
! Subprogram not used                                          dataGpr(if,c_lnd_lr,ip), &
! Subprogram not used                                          dataGpr(if,c_ocn_or,ip), &
! Subprogram not used                                          dataGpr(if,c_inh_ir,ip), &
! Subprogram not used                                          dataGpr(if,c_ish_ir,ip), &
! Subprogram not used                                          dataGpr(if,c_atm_ar,ip)+ &
! Subprogram not used                                          dataGpr(if,c_lnd_lr,ip)+ &
! Subprogram not used                                          dataGpr(if,c_ocn_or,ip)+ &
! Subprogram not used                                          dataGpr(if,c_inh_ir,ip)+ &
! Subprogram not used                                          dataGpr(if,c_ish_ir,ip)
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,'NET HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
! Subprogram not used       do if = f_h, f_w-1
! Subprogram not used          write(logunit,FA1r)   fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
! Subprogram not used                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
! Subprogram not used                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
! Subprogram not used                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
! Subprogram not used                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
! Subprogram not used                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
! Subprogram not used                                          dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip), &
! Subprogram not used                                          dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
! Subprogram not used                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
! Subprogram not used                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
! Subprogram not used                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
! Subprogram not used                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
! Subprogram not used                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)+ &
! Subprogram not used                                          dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip)
! Subprogram not used       enddo
! Subprogram not used       write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_h:f_w-1,c_atm_ar,ip))+sum(dataGpr(f_h:f_w-1,c_atm_as,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_lnd_lr,ip))+sum(dataGpr(f_h:f_w-1,c_lnd_ls,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_rof_rr,ip))+sum(dataGpr(f_h:f_w-1,c_rof_rs,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_ocn_or,ip))+sum(dataGpr(f_h:f_w-1,c_ocn_os,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_inh_ir,ip))+sum(dataGpr(f_h:f_w-1,c_inh_is,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_ish_ir,ip))+sum(dataGpr(f_h:f_w-1,c_ish_is,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_glc_gr,ip))+sum(dataGpr(f_h:f_w-1,c_glc_gs,ip)), &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_atm_ar,ip))+sum(dataGpr(f_h:f_w-1,c_atm_as,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_lnd_lr,ip))+sum(dataGpr(f_h:f_w-1,c_lnd_ls,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_rof_rr,ip))+sum(dataGpr(f_h:f_w-1,c_rof_rs,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_ocn_or,ip))+sum(dataGpr(f_h:f_w-1,c_ocn_os,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_inh_ir,ip))+sum(dataGpr(f_h:f_w-1,c_inh_is,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_ish_ir,ip))+sum(dataGpr(f_h:f_w-1,c_ish_is,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_h:f_w-1,c_glc_gr,ip))+sum(dataGpr(f_h:f_w-1,c_glc_gs,ip))
! Subprogram not used 
! Subprogram not used       write(logunit,*) ' '
! Subprogram not used       write(logunit,FAH) subname,'NET WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
! Subprogram not used       write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
! Subprogram not used       do if = f_w, f_size
! Subprogram not used          write(logunit,FA1r)   fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
! Subprogram not used                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
! Subprogram not used                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
! Subprogram not used                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
! Subprogram not used                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
! Subprogram not used                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
! Subprogram not used                                          dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip), &
! Subprogram not used                                          dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
! Subprogram not used                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
! Subprogram not used                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
! Subprogram not used                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
! Subprogram not used                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
! Subprogram not used                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)+ &
! Subprogram not used                                          dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip)
! Subprogram not used       enddo
! Subprogram not used       write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_w:f_size,c_atm_ar,ip))+sum(dataGpr(f_w:f_size,c_atm_as,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_lnd_lr,ip))+sum(dataGpr(f_w:f_size,c_lnd_ls,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_rof_rr,ip))+sum(dataGpr(f_w:f_size,c_rof_rs,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_ocn_or,ip))+sum(dataGpr(f_w:f_size,c_ocn_os,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_inh_ir,ip))+sum(dataGpr(f_w:f_size,c_inh_is,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_ish_ir,ip))+sum(dataGpr(f_w:f_size,c_ish_is,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_glc_gr,ip))+sum(dataGpr(f_w:f_size,c_glc_gs,ip)), &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_atm_ar,ip))+sum(dataGpr(f_w:f_size,c_atm_as,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_lnd_lr,ip))+sum(dataGpr(f_w:f_size,c_lnd_ls,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_rof_rr,ip))+sum(dataGpr(f_w:f_size,c_rof_rs,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_ocn_or,ip))+sum(dataGpr(f_w:f_size,c_ocn_os,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_inh_ir,ip))+sum(dataGpr(f_w:f_size,c_inh_is,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_ish_ir,ip))+sum(dataGpr(f_w:f_size,c_ish_is,ip))+ &
! Subprogram not used                                     sum(dataGpr(f_w:f_size,c_glc_gr,ip))+sum(dataGpr(f_w:f_size,c_glc_gs,ip))
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    write(logunit,*) ' '
! Subprogram not used ! ---- doprint ---- doprint ---- doprint ----
! Subprogram not used    endif  ! plev > 0
! Subprogram not used    enddo  ! ip = 1,p_size
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_print_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_avect_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global diagnostics for AV/ID.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used SUBROUTINE seq_diag_avect_mct(infodata, id, av, dom, gsmap, comment)
! Subprogram not used 
! Subprogram not used    use seq_infodata_mod
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(seq_infodata_type) , intent(in)           :: infodata
! Subprogram not used    integer(in)             , intent(in)           :: ID
! Subprogram not used    type(mct_aVect)         , intent(in)           :: av
! Subprogram not used    type(mct_gGrid)         , pointer              :: dom
! Subprogram not used    type(mct_gsMap)         , pointer              :: gsmap
! Subprogram not used    character(len=*)        , intent(in), optional :: comment
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    logical                          :: bfbflag
! Subprogram not used    integer(in)                      :: n,k         ! counters
! Subprogram not used    integer(in)                      :: npts,nptsg  ! number of local/global pts in AV
! Subprogram not used    integer(in)                      :: kflds       ! number of fields in AV
! Subprogram not used    real(r8),                pointer :: sumbuf (:)  ! sum buffer
! Subprogram not used    real(r8),                pointer :: minbuf (:)  ! min buffer
! Subprogram not used    real(r8),                pointer :: maxbuf (:)  ! max buffer
! Subprogram not used    real(r8),                pointer :: sumbufg(:)  ! sum buffer reduced
! Subprogram not used    real(r8),                pointer :: minbufg(:)  ! min buffer reduced
! Subprogram not used    real(r8),                pointer :: maxbufg(:)  ! max buffer reduced
! Subprogram not used    integer(i8),             pointer :: isumbuf (:) ! integer local sum
! Subprogram not used    integer(i8),             pointer :: isumbufg(:) ! integer global sum
! Subprogram not used    integer(i8)                      :: ihuge       ! huge
! Subprogram not used    integer(in)                      :: mpicom      ! mpi comm
! Subprogram not used    integer(in)                      :: iam         ! pe number
! Subprogram not used    integer(in)                      :: km,ka       ! field indices
! Subprogram not used    integer(in)                      :: ns          ! size of local AV
! Subprogram not used    integer(in)                      :: rcode       ! error code
! Subprogram not used    real(r8),                pointer :: weight(:)   ! weight
! Subprogram not used    type(mct_string)                 :: mstring     ! mct char type
! Subprogram not used    character(CL)                    :: lcomment    ! should be long enough
! Subprogram not used    character(CL)                    :: itemc       ! string converted to char
! Subprogram not used 
! Subprogram not used    type(mct_avect)                  :: AV1         ! local avect with one field
! Subprogram not used    type(mct_avect)                  :: AVr1        ! avect on root with one field
! Subprogram not used    type(mct_avect)                  :: AVr2        ! avect on root with one field
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_avect_mct) '
! Subprogram not used    character(*),parameter :: F00   = "('(seq_diag_avect_mct) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! print instantaneous budget data
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call seq_comm_setptrs(ID,&
! Subprogram not used         mpicom=mpicom, iam=iam)
! Subprogram not used 
! Subprogram not used    call seq_infodata_GetData(infodata,&
! Subprogram not used         bfbflag=bfbflag)
! Subprogram not used 
! Subprogram not used    lcomment = ''
! Subprogram not used    if (present(comment)) then
! Subprogram not used       lcomment=trim(comment)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    ns = mct_aVect_lsize(AV)
! Subprogram not used    npts = mct_aVect_lsize(dom%data)
! Subprogram not used    if (ns /= npts) call shr_sys_abort(trim(subname)//' ERROR: size of AV,dom')
! Subprogram not used    km = mct_aVect_indexRA(dom%data,'mask')
! Subprogram not used    ka = mct_aVect_indexRA(dom%data,afldname)
! Subprogram not used    kflds = mct_aVect_nRattr(AV)
! Subprogram not used    allocate(sumbuf(kflds),sumbufg(kflds))
! Subprogram not used 
! Subprogram not used    sumbuf =       0.0_r8
! Subprogram not used 
! Subprogram not used    if (bfbflag) then
! Subprogram not used 
! Subprogram not used       npts = mct_aVect_lsize(AV)
! Subprogram not used       allocate(weight(npts))
! Subprogram not used       weight(:) = 1.0_r8 
! Subprogram not used       do n = 1,npts
! Subprogram not used          if (dom%data%rAttr(km,n) <= 1.0e-06_R8) then
! Subprogram not used             weight(n) = 0.0_r8
! Subprogram not used          else
! Subprogram not used             weight(n) = dom%data%rAttr(ka,n)*shr_const_rearth*shr_const_rearth
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       allocate(maxbuf(kflds),maxbufg(kflds))
! Subprogram not used       maxbuf = 0.0_r8
! Subprogram not used 
! Subprogram not used       do n = 1,npts
! Subprogram not used       do k = 1,kflds
! Subprogram not used          if (.not. shr_const_isspval(AV%rAttr(k,n))) then
! Subprogram not used              maxbuf(k) = max(maxbuf(k),abs(AV%rAttr(k,n)*weight(n)))
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       call shr_mpi_max(maxbuf,maxbufg,mpicom,subname,all=.true.)
! Subprogram not used       call shr_mpi_sum(npts,nptsg,mpicom,subname,all=.true.)
! Subprogram not used 
! Subprogram not used       do k = 1,kflds
! Subprogram not used          if (maxbufg(k) < 1000.0*TINY(maxbufg(k)) .or. &
! Subprogram not used              maxbufg(k) > HUGE(maxbufg(k))/(2.0_r8*nptsg)) then
! Subprogram not used             maxbufg(k) = 0.0_r8
! Subprogram not used          else
! Subprogram not used             maxbufg(k) = (1.1_r8) * maxbufg(k) * nptsg
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       allocate(isumbuf(kflds),isumbufg(kflds))
! Subprogram not used       isumbuf = 0
! Subprogram not used       ihuge = HUGE(isumbuf)
! Subprogram not used 
! Subprogram not used       do n = 1,npts
! Subprogram not used       do k = 1,kflds
! Subprogram not used          if (.not. shr_const_isspval(AV%rAttr(k,n))) then
! Subprogram not used              if (abs(maxbufg(k)) > 1000.0_r8 * TINY(maxbufg)) then
! Subprogram not used                 isumbuf(k) = isumbuf(k) + int((AV%rAttr(k,n)*weight(n)/maxbufg(k))*ihuge,i8)
! Subprogram not used              endif
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       call shr_mpi_sum(isumbuf,isumbufg,mpicom,subname)
! Subprogram not used 
! Subprogram not used       do k = 1,kflds
! Subprogram not used          sumbufg(k) = isumbufg(k)*maxbufg(k)/ihuge
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       deallocate(weight)
! Subprogram not used       deallocate(maxbuf,maxbufg)
! Subprogram not used       deallocate(isumbuf,isumbufg)
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used       npts = mct_aVect_lsize(AV)
! Subprogram not used       allocate(weight(npts))
! Subprogram not used       weight(:) = 1.0_r8 
! Subprogram not used       do n = 1,npts
! Subprogram not used          if (dom%data%rAttr(km,n) <= 1.0e-06_R8) then
! Subprogram not used             weight(n) = 0.0_r8
! Subprogram not used          else
! Subprogram not used             weight(n) = dom%data%rAttr(ka,n)*shr_const_rearth*shr_const_rearth
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do n = 1,npts
! Subprogram not used       do k = 1,kflds
! Subprogram not used          if (.not. shr_const_isspval(AV%rAttr(k,n))) then
! Subprogram not used              sumbuf(k) = sumbuf(k) + AV%rAttr(k,n)*weight(n)
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       !--- global reduction ---
! Subprogram not used       call shr_mpi_sum(sumbuf,sumbufg,mpicom,subname)
! Subprogram not used 
! Subprogram not used       deallocate(weight)
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (iam == 0) then
! Subprogram not used       !      write(logunit,*) 'sdAV: *** writing ',trim(lcomment),': k fld min/max/sum ***'
! Subprogram not used       do k = 1,kflds
! Subprogram not used          call mct_aVect_getRList(mstring,k,AV)
! Subprogram not used          itemc = mct_string_toChar(mstring)
! Subprogram not used          call mct_string_clean(mstring)
! Subprogram not used          if (len_trim(lcomment) > 0) then
! Subprogram not used             write(logunit,100) 'xxx','sorr',k,sumbufg(k),trim(lcomment),trim(itemc)
! Subprogram not used          else
! Subprogram not used             write(logunit,101) 'xxx','sorr',k,sumbufg(k),trim(itemc)
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       call shr_sys_flush(logunit)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    deallocate(sumbuf,sumbufg)
! Subprogram not used 
! Subprogram not used 100  format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a,1x,a)
! Subprogram not used 101  format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a)
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_avect_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_avdiff_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global diagnostics for AV/ID.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used SUBROUTINE seq_diag_avdiff_mct(AV1,AV2,ID,comment)
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(mct_aVect) , intent(in) :: AV1
! Subprogram not used    type(mct_aVect) , intent(in) :: AV2
! Subprogram not used    integer         , intent(in) :: ID
! Subprogram not used    character(len=*), intent(in), optional :: comment
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    integer(in)      :: n,k,n1,k1,n2,k2         ! counters
! Subprogram not used    integer(in)      :: iam         ! pe number
! Subprogram not used    integer(in)      :: cnt         ! counter
! Subprogram not used    real(r8)         :: adiff,rdiff ! diff values
! Subprogram not used    type(mct_string) :: mstring     ! mct char type
! Subprogram not used    character(len=64):: lcomment    ! should be long enough
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName = '(seq_diag_avdiff_mct) '
! Subprogram not used    character(*),parameter :: F00   = "('(seq_diag_avdiff_mct) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! print instantaneous budget data
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call seq_comm_setptrs(ID,iam=iam)
! Subprogram not used 
! Subprogram not used    lcomment = ''
! Subprogram not used    if (present(comment)) then
! Subprogram not used       lcomment=trim(comment)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    n1 = mct_aVect_lsize(AV1)
! Subprogram not used    k1 = mct_aVect_nRattr(AV1)
! Subprogram not used    n2 = mct_aVect_lsize(AV2)
! Subprogram not used    k2 = mct_aVect_nRattr(AV2)
! Subprogram not used 
! Subprogram not used    if (n1 /= n2 .or. k1 /= k2) then
! Subprogram not used       write(s_logunit,*) subname,trim(lcomment),' AV sizes different ',n1,n2,k1,k2
! Subprogram not used       return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    do k = 1,k1
! Subprogram not used       cnt = 0
! Subprogram not used       adiff = 0.
! Subprogram not used       rdiff = 0.
! Subprogram not used       do n = 1,n1
! Subprogram not used          if (AV1%rAttr(k,n) /= AV2%rAttr(k,n)) then
! Subprogram not used             cnt = cnt + 1
! Subprogram not used             adiff = max(adiff, abs(AV1%rAttr(k,n)-AV2%rAttr(k,n)))
! Subprogram not used             rdiff = max(rdiff, abs(AV1%rAttr(k,n)-AV2%rAttr(k,n))/(abs(AV1%rAttr(k,n))+abs(AV2%rAttr(k,n))))
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       if (cnt > 0) then
! Subprogram not used          call mct_aVect_getRList(mstring,k,AV1)
! Subprogram not used          write(s_logunit,*) subname,trim(lcomment),' AVs fld k diff ', &
! Subprogram not used             iam,mct_string_toChar(mstring),cnt,adiff,rdiff, &
! Subprogram not used             minval(AV1%rAttr(k,:)),minval(AV1%rAttr(k,:)), &
! Subprogram not used             maxval(AV1%rAttr(k,:)),maxval(AV2%rAttr(k,:))
! Subprogram not used          call mct_string_clean(mstring)       
! Subprogram not used       endif
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used end subroutine seq_diag_avdiff_mct

!===============================================================================
end module seq_diag_mct
