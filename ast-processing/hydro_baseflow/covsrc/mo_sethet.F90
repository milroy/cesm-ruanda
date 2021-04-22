
module mo_sethet

!
! LKE (10/11/2010): added HCN, CH3CN, HCOOH  to cesm1_0_beta07_offline version
!                   HCN, CH3CN have new Henry's Law coefficients, HCOOH is set to CH3COOH
! LKE (10/18/2010): SO2 washout corrected based on recommendation of R.Easter, PNNL
!
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_logfile,     only: iulog
  use gas_wetdep_opts, only: gas_wetdep_cnt, gas_wetdep_method, gas_wetdep_list
  use phys_control,    only: phys_getopts

  private
  public :: sethet_inti, sethet

  save

  integer :: h2o2_ndx, hno3_ndx, ch2o_ndx, ch3ooh_ndx, ch3coooh_ndx, &
       ho2no2_ndx, ch3cocho_ndx, xooh_ndx, onitr_ndx, glyald_ndx, &
       ch3cho_ndx, mvk_ndx, macr_ndx, pooh_ndx, c2h5ooh_ndx, &
       c3h7ooh_ndx, rooh_ndx, isopno3_ndx, onit_ndx, Pb_ndx, &
       macrooh_ndx, isopooh_ndx, ch3oh_ndx, c2h5oh_ndx, hyac_ndx, hydrald_ndx
  integer :: spc_h2o2_ndx, spc_hno3_ndx
  integer :: spc_so2_ndx
  integer :: spc_sogm_ndx, spc_sogi_ndx, spc_sogt_ndx, spc_sogb_ndx, spc_sogx_ndx

  integer :: alkooh_ndx, mekooh_ndx, tolooh_ndx, terpooh_ndx, ch3cooh_ndx
  integer :: so2_ndx, soa_ndx, so4_ndx, cb2_ndx, oc2_ndx, nh3_ndx, nh4no3_ndx, &
             sa1_ndx, sa2_ndx, sa3_ndx, sa4_ndx, nh4_ndx, h2so4_ndx
  integer :: xisopno3_ndx,xho2no2_ndx,xonitr_ndx,xhno3_ndx,xonit_ndx
  integer :: clono2_ndx, brono2_ndx, hcl_ndx, n2o5_ndx, hocl_ndx, hobr_ndx, hbr_ndx 
  integer :: ch3cn_ndx, hcn_ndx, hcooh_ndx
  integer, allocatable :: wetdep_map(:)
  integer :: sogm_ndx, sogi_ndx, sogt_ndx, sogb_ndx, sogx_ndx
  logical :: do_wetdep

  ! prognostic modal aerosols
  logical :: prog_modal_aero

contains

  subroutine sethet_inti
    !-----------------------------------------------------------------------      
    !       ... intialize the wet removal rate constants routine
    !-----------------------------------------------------------------------      

    use mo_chem_utls, only : get_het_ndx, get_spc_ndx
    use spmd_utils,   only : masterproc
    use abortutils,   only : endrun

    integer :: k, m
    
    do_wetdep = gas_wetdep_cnt>0 .and. gas_wetdep_method=='MOZ'
    if ( .not. do_wetdep) return

    call phys_getopts( prog_modal_aero_out = prog_modal_aero )

    allocate( wetdep_map(gas_wetdep_cnt))

    do k=1,gas_wetdep_cnt
       m = get_het_ndx( trim(gas_wetdep_list(k))) 
       if (m>0) then
          wetdep_map(k) = m
       else
          call endrun('sethet_inti: cannot map '//trim(gas_wetdep_list(k)))
       endif
    enddo

    xisopno3_ndx = get_het_ndx( 'XISOPNO3' )
    xho2no2_ndx  = get_het_ndx( 'XHO2NO2' )
    xonitr_ndx   = get_het_ndx( 'XONITR' )
    xhno3_ndx    = get_het_ndx( 'XHNO3' )
    xonit_ndx    = get_het_ndx( 'XONIT' )

    spc_h2o2_ndx = get_spc_ndx( 'H2O2' )
    spc_hno3_ndx = get_spc_ndx( 'HNO3' )
    spc_so2_ndx  = get_spc_ndx( 'SO2' )

    clono2_ndx = get_het_ndx( 'CLONO2' )
    brono2_ndx = get_het_ndx( 'BRONO2' )
    hcl_ndx    = get_het_ndx( 'HCL' )
    n2o5_ndx   = get_het_ndx( 'N2O5' )
    hocl_ndx   = get_het_ndx( 'HOCL' )
    hobr_ndx   = get_het_ndx( 'HOBR' )
    hbr_ndx    = get_het_ndx( 'HBR' )

    h2o2_ndx   = get_het_ndx( 'H2O2' )
    hno3_ndx   = get_het_ndx( 'HNO3' )
    ch2o_ndx   = get_het_ndx( 'CH2O' )
    ch3ooh_ndx = get_het_ndx( 'CH3OOH' )
    ch3coooh_ndx = get_het_ndx( 'CH3COOOH' )
    ho2no2_ndx  = get_het_ndx( 'HO2NO2' )
    ch3cocho_ndx = get_het_ndx( 'CH3COCHO' )
    xooh_ndx    = get_het_ndx( 'XOOH' )
    onitr_ndx   = get_het_ndx( 'ONITR' )
    glyald_ndx  = get_het_ndx( 'GLYALD' )
    ch3cho_ndx  = get_het_ndx( 'CH3CHO' )
    mvk_ndx     = get_het_ndx( 'MVK' )
    macr_ndx    = get_het_ndx( 'MACR' )
    pooh_ndx    = get_het_ndx( 'POOH' )
    c2h5ooh_ndx = get_het_ndx( 'C2H5OOH' )
    c3h7ooh_ndx = get_het_ndx( 'C3H7OOH' )
    rooh_ndx    = get_het_ndx( 'ROOH' )
    isopno3_ndx = get_het_ndx( 'ISOPNO3' )
    onit_ndx    = get_het_ndx( 'ONIT' )
    Pb_ndx      = get_het_ndx( 'Pb' )
    macrooh_ndx = get_het_ndx( 'MACROOH' )
    isopooh_ndx = get_het_ndx( 'ISOPOOH' )
    ch3oh_ndx   = get_het_ndx( 'CH3OH' )
    c2h5oh_ndx  = get_het_ndx( 'C2H5OH' )
    hyac_ndx    = get_het_ndx( 'HYAC' )
    hydrald_ndx = get_het_ndx( 'HYDRALD' )
    alkooh_ndx  = get_het_ndx( 'ALKOOH' )
    mekooh_ndx  = get_het_ndx( 'MEKOOH' )
    tolooh_ndx  = get_het_ndx( 'TOLOOH' )
    terpooh_ndx = get_het_ndx( 'TERPOOH' )
    ch3cooh_ndx = get_het_ndx( 'CH3COOH' )
    so2_ndx     = get_het_ndx( 'SO2' )
    soa_ndx     = get_het_ndx( 'SOA' )
    sogb_ndx    = get_het_ndx( 'SOGB' )
    sogi_ndx    = get_het_ndx( 'SOGI' )
    sogm_ndx    = get_het_ndx( 'SOGM' )
    sogt_ndx    = get_het_ndx( 'SOGT' )
    sogx_ndx    = get_het_ndx( 'SOGX' )
    so4_ndx     = get_het_ndx( 'SO4' )
    cb2_ndx     = get_het_ndx( 'CB2' )
    oc2_ndx     = get_het_ndx( 'OC2' )
    nh3_ndx     = get_het_ndx( 'NH3' )
    nh4no3_ndx  = get_het_ndx( 'NH4NO3' )
    nh4_ndx     = get_het_ndx( 'NH4' )
    h2so4_ndx   = get_het_ndx( 'H2SO4' )
    sa1_ndx     = get_het_ndx( 'SA1' )
    sa2_ndx     = get_het_ndx( 'SA2' )
    sa3_ndx     = get_het_ndx( 'SA3' )
    sa4_ndx     = get_het_ndx( 'SA4' )
    ch3cn_ndx   = get_het_ndx( 'CH3CN' )
    hcn_ndx     = get_het_ndx( 'HCN' )
    hcooh_ndx   = get_het_ndx( 'HCOOH' )

    if (masterproc) then
       write(iulog,*) 'sethet_inti: new ndx ',so2_ndx,soa_ndx,so4_ndx,cb2_ndx,oc2_ndx, &
            nh3_ndx,nh4no3_ndx,sa1_ndx,sa2_ndx,sa3_ndx,sa4_ndx
       write(iulog,*) ' '
       write(iulog,*) 'sethet_inti: diagnotics '
       write(iulog,'(10i5)') h2o2_ndx, hno3_ndx, ch2o_ndx, ch3ooh_ndx, ch3coooh_ndx, &
            ho2no2_ndx, ch3cocho_ndx, xooh_ndx, onitr_ndx, glyald_ndx, &
            ch3cho_ndx, mvk_ndx, macr_ndx, pooh_ndx, c2h5ooh_ndx, &
            c3h7ooh_ndx, rooh_ndx, isopno3_ndx, onit_ndx, Pb_ndx, &
            macrooh_ndx, isopooh_ndx, ch3oh_ndx, c2h5oh_ndx, hyac_ndx, hydrald_ndx
    endif

  end subroutine sethet_inti

! Subprogram not used   subroutine sethet( het_rates, press, zmid,  phis, tfld, &
! Subprogram not used                      cmfdqr, nrain, nevapr, delt, xhnm, &
! Subprogram not used                      qin, ncol, lchnk )
! Subprogram not used     !-----------------------------------------------------------------------      
! Subprogram not used     !       ... compute rainout loss rates (1/s)
! Subprogram not used     !-----------------------------------------------------------------------      
! Subprogram not used 
! Subprogram not used     use physconst,    only : rga,pi
! Subprogram not used     use chem_mods,    only : gas_pcnst
! Subprogram not used     use ppgrid,       only : pver, pcols
! Subprogram not used     use phys_grid,    only : get_rlat_all_p
! Subprogram not used     use abortutils,   only : endrun
! Subprogram not used     use mo_constants, only : avo => avogadro, boltz_cgs
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------      
! Subprogram not used     !       ... dummy arguments
! Subprogram not used     !-----------------------------------------------------------------------      
! Subprogram not used     integer, intent(in)   ::    ncol                        ! columns in chunk
! Subprogram not used     integer, intent(in)   ::    lchnk                       ! chunk index
! Subprogram not used     real(r8), intent(in)  ::    delt                        ! time step ( s )
! Subprogram not used     real(r8), intent(in)  ::    press(pcols,pver)           ! pressure in pascals
! Subprogram not used     real(r8), intent(in)  ::    cmfdqr(ncol,pver)           ! dq/dt for convection
! Subprogram not used     real(r8), intent(in)  ::    nrain(ncol,pver)            ! stratoform precip
! Subprogram not used     real(r8), intent(in)  ::    nevapr(ncol,pver)           ! evaporation
! Subprogram not used     real(r8), intent(in)  ::    qin(ncol,pver,gas_pcnst)    ! xported species ( vmr )
! Subprogram not used     real(r8), intent(in)  ::    zmid(ncol,pver)             ! midpoint geopot (km)
! Subprogram not used     real(r8), intent(in)  ::    phis(pcols)                 ! surf geopot
! Subprogram not used     real(r8), intent(in)  ::    tfld(pcols,pver)            ! temperature (k)
! Subprogram not used     real(r8), intent(in)  ::    xhnm(ncol,pver)             ! total atms density ( /cm^3)
! Subprogram not used     real(r8), intent(out) ::    het_rates(ncol,pver,gas_pcnst) ! rainout loss rates
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------      
! Subprogram not used     !       ... local variables
! Subprogram not used     !-----------------------------------------------------------------------      
! Subprogram not used     real(r8), parameter ::  xrm   = .189_r8             ! mean diameter of rain drop (cm)
! Subprogram not used     real(r8), parameter ::  xum   = 748._r8             ! mean rain drop terminal velocity (cm/s)
! Subprogram not used     real(r8), parameter ::  xvv   = 6.18e-2_r8          ! kinetic viscosity (cm^2/s)
! Subprogram not used     real(r8), parameter ::  xdg   = .112_r8             ! mass transport coefficient (cm/s)
! Subprogram not used     real(r8), parameter ::  t0    = 298._r8             ! reference temperature (k)
! Subprogram not used     real(r8), parameter ::  xph0  = 1.e-5_r8            ! cloud [h+]
! Subprogram not used     real(r8), parameter ::  satf_hno3  = .016_r8        ! saturation factor for hno3 in clouds 
! Subprogram not used     real(r8), parameter ::  satf_h2o2  = .016_r8        ! saturation factor for h2o2 in clouds 
! Subprogram not used     real(r8), parameter ::  satf_so2   = .016_r8        ! saturation factor for so2 in clouds 
! Subprogram not used     real(r8), parameter ::  satf_ch2o  = .1_r8          ! saturation factor for ch2o in clouds 
! Subprogram not used     real(r8), parameter ::  satf_sog  =  .016_r8        ! saturation factor for sog in clouds
! Subprogram not used     real(r8), parameter ::  const0   = boltz_cgs * 1.e-6_r8 ! (atmospheres/deg k/cm^3)
! Subprogram not used     real(r8), parameter ::  hno3_diss = 15.4_r8         ! hno3 dissociation constant
! Subprogram not used     real(r8), parameter ::  geo_fac  = 6._r8            ! geometry factor (surf area/volume = geo_fac/diameter)
! Subprogram not used     real(r8), parameter ::  mass_air = 29._r8           ! mass of background atmosphere (amu)
! Subprogram not used     real(r8), parameter ::  mass_h2o = 18._r8           ! mass of water vapor (amu)
! Subprogram not used     real(r8), parameter ::  h2o_mol  = 1.e3_r8/mass_h2o ! (gm/mol water)
! Subprogram not used     real(r8), parameter ::  km2cm    = 1.e5_r8          ! convert km to cm
! Subprogram not used     real(r8), parameter ::  m2km     = 1.e-3_r8         ! convert m to km
! Subprogram not used     real(r8), parameter ::  cm3_2_m3 = 1.e-6_r8         ! convert cm^3 to m^3
! Subprogram not used     real(r8), parameter ::  m3_2_cm3 = 1.e6_r8          ! convert m^3 to cm^3
! Subprogram not used     real(r8), parameter ::  liter_per_gram = 1.e-3_r8
! Subprogram not used     real(r8), parameter ::  avo2  = avo * liter_per_gram * cm3_2_m3 ! (liter/gm/mol*(m/cm)^3)
! Subprogram not used 
! Subprogram not used     integer  ::      i, m, k, kk                 ! indicies
! Subprogram not used     real(r8) ::      xkgm                        ! mass flux on rain drop
! Subprogram not used     real(r8) ::      all1, all2                  ! work variables
! Subprogram not used     real(r8) ::      stay                        ! fraction of layer traversed by falling drop in timestep delt
! Subprogram not used     real(r8) ::      xeqca1, xeqca2, xca1, xca2, xdtm
! Subprogram not used     real(r8) ::      xxx1, xxx2, yhno3, yh2o2
! Subprogram not used     real(r8) ::      all3, xeqca3, xca3, xxx3, yso2, so2_diss
! Subprogram not used     real(r8) ::      all4, xeqca4, xca4, xxx4
! Subprogram not used     real(r8) ::      all5, xeqca5, xca5, xxx5
! Subprogram not used     real(r8) ::      all6, xeqca6, xca6, xxx6
! Subprogram not used     real(r8) ::      all7, xeqca7, xca7, xxx7
! Subprogram not used     real(r8) ::      all8, xeqca8, xca8, xxx8
! Subprogram not used     real(r8) ::      ysogm,ysogi,ysogt,ysogb,ysogx
! Subprogram not used 
! Subprogram not used     real(r8), dimension(ncol)  :: &
! Subprogram not used          xk0, work1, work2, work3, zsurf
! Subprogram not used     real(r8), dimension(pver)  :: &
! Subprogram not used          xgas1, xgas2
! Subprogram not used     real(r8), dimension(pver)  :: xgas3, xgas4, xgas5, xgas6, xgas7, xgas8
! Subprogram not used     real(r8), dimension(ncol)  :: &
! Subprogram not used          tmp0_rates, tmp1_rates
! Subprogram not used     real(r8), dimension(ncol,pver)  :: &
! Subprogram not used          delz, &              ! layer depth about interfaces (cm)
! Subprogram not used          xhno3, &             ! hno3 concentration (molecules/cm^3)
! Subprogram not used          xh2o2, &             ! h2o2 concentration (molecules/cm^3)
! Subprogram not used          xso2, &              ! so2 concentration (molecules/cm^3)
! Subprogram not used          xsogm, &             ! sogm concentration (molecules/cm^3)
! Subprogram not used          xsogi, &             ! sogi concentration (molecules/cm^3)
! Subprogram not used          xsogt, &             ! sogt concentration (molecules/cm^3)
! Subprogram not used          xsogb, &             ! sogb concentration (molecules/cm^3)
! Subprogram not used          xsogx, &             ! sogx concentration (molecules/cm^3)
! Subprogram not used          xliq, &              ! liquid rain water content in a grid cell (gm/m^3)
! Subprogram not used          rain                 ! conversion rate of water vapor into rain water (molecules/cm^3/s)
! Subprogram not used     real(r8), dimension(ncol,pver)  :: &
! Subprogram not used          xhen_hno3, xhen_h2o2, xhen_ch2o, xhen_ch3ooh, xhen_ch3co3h, &
! Subprogram not used          xhen_ch3cocho, xhen_xooh, xhen_onitr, xhen_ho2no2, xhen_glyald, &
! Subprogram not used          xhen_ch3cho, xhen_mvk, xhen_macr,xhen_sog
! Subprogram not used     real(r8), dimension(ncol,pver)  :: &
! Subprogram not used          xhen_nh3, xhen_ch3cooh
! Subprogram not used     real(r8), dimension(ncol,pver,8) :: tmp_hetrates
! Subprogram not used     real(r8), dimension(ncol,pver)  :: precip
! Subprogram not used     real(r8), dimension(ncol,pver)  :: xhen_hcn, xhen_ch3cn, xhen_so2
! Subprogram not used 
! Subprogram not used     integer    ::      ktop_all       
! Subprogram not used     integer    ::      ktop(ncol)                  ! 100 mb level
! Subprogram not used 
! Subprogram not used     real(r8) :: rlat(pcols)                       ! latitude in radians for columns
! Subprogram not used     real(r8) :: p_limit
! Subprogram not used     real(r8), parameter :: d2r = pi/180._r8
! Subprogram not used !
! Subprogram not used ! jfl : new variables for rescaling sum of positive values to actual amount
! Subprogram not used !
! Subprogram not used     real(r8) :: total_rain,total_pos
! Subprogram not used     character(len=3) :: hetratestrg
! Subprogram not used     real(r8), parameter :: MISSING = -999999._r8
! Subprogram not used     integer ::  mm
! Subprogram not used 
! Subprogram not used !
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !        note: the press array is in pascals and must be
! Subprogram not used     !              mutiplied by 10 to yield dynes/cm**2.
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !       ... set wet deposition for
! Subprogram not used     !           1. h2o2         2. hno3
! Subprogram not used     !           3. ch2o         4. ch3ooh
! Subprogram not used     !           5. pooh         6. ch3coooh
! Subprogram not used     !           7. ho2no2       8. onit
! Subprogram not used     !           9. mvk         10. macr
! Subprogram not used     !          11. c2h5ooh     12. c3h7ooh
! Subprogram not used     !          13. rooh        14. ch3cocho
! Subprogram not used     !          15. pb          16. macrooh
! Subprogram not used     !          17. xooh        18. onitr
! Subprogram not used     !          19. isopooh     20. ch3oh
! Subprogram not used     !          21. c2h5oh      22. glyald
! Subprogram not used     !          23. hyac        24. hydrald
! Subprogram not used     !          25. ch3cho      26. isopno3
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     het_rates(:,:,:) = 0._r8
! Subprogram not used 
! Subprogram not used     if ( .not. do_wetdep) return
! Subprogram not used 
! Subprogram not used     call get_rlat_all_p(lchnk, ncol, rlat)
! Subprogram not used 
! Subprogram not used     do mm = 1,gas_wetdep_cnt
! Subprogram not used        m = wetdep_map(mm)
! Subprogram not used        if ( m>0 ) then
! Subprogram not used           het_rates(:,:,m) = MISSING
! Subprogram not used        endif
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !	... the 2 and .6 multipliers are from a formula by frossling (1938)
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     xkgm = xdg/xrm * 2._r8 + xdg/xrm * .6_r8 * sqrt( xrm*xum/xvv ) * (xvv/xdg)**(1._r8/3._r8) 
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !	... Find 100 mb level
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     do i = 1,ncol
! Subprogram not used        if ( abs(rlat(i)) > 60._r8*d2r ) then
! Subprogram not used           p_limit = 300.e2_r8
! Subprogram not used        else
! Subprogram not used           p_limit = 100.e2_r8 
! Subprogram not used        endif
! Subprogram not used        k_loop: do k = pver,1,-1
! Subprogram not used           if( press(i,k) < p_limit ) then
! Subprogram not used              ktop(i) = k
! Subprogram not used              exit k_loop
! Subprogram not used           end if
! Subprogram not used        end do k_loop
! Subprogram not used     end do
! Subprogram not used     ktop_all = minval( ktop(:) )
! Subprogram not used !
! Subprogram not used ! jfl
! Subprogram not used !
! Subprogram not used ! this is added to rescale the variable precip (which can only be positive)
! Subprogram not used ! to the actual vertical integral of positive and negative values.  This
! Subprogram not used ! removes point storms
! Subprogram not used !
! Subprogram not used     do i = 1,ncol
! Subprogram not used        total_rain = 0._r8
! Subprogram not used        total_pos  = 0._r8
! Subprogram not used        do k = 1,pver
! Subprogram not used           precip(i,k) = cmfdqr(i,k) + nrain(i,k) - nevapr(i,k)
! Subprogram not used           total_rain = total_rain + precip(i,k)
! Subprogram not used           if ( precip(i,k) < 0._r8 ) precip(i,k) = 0._r8
! Subprogram not used           total_pos  = total_pos  + precip(i,k)
! Subprogram not used        end do
! Subprogram not used        if ( total_rain <= 0._r8 ) then
! Subprogram not used           precip(i,:) = 0._r8
! Subprogram not used        else
! Subprogram not used           do k = 1,pver
! Subprogram not used              precip(i,k) = precip(i,k) * total_rain/total_pos
! Subprogram not used           end do
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do k = 1,pver
! Subprogram not used        !jfl       precip(:ncol,k) = cmfdqr(:ncol,k) + nrain(:ncol,k) - nevapr(:ncol,k)
! Subprogram not used        rain(:ncol,k)   = mass_air*precip(:ncol,k)*xhnm(:ncol,k) / mass_h2o
! Subprogram not used        xliq(:ncol,k)   = precip(:ncol,k) * delt * xhnm(:ncol,k) / avo*mass_air * m3_2_cm3
! Subprogram not used        if( spc_hno3_ndx > 0 ) then
! Subprogram not used           xhno3(:ncol,k)  = qin(:ncol,k,spc_hno3_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xhno3(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_h2o2_ndx > 0 ) then
! Subprogram not used           xh2o2(:ncol,k)  = qin(:ncol,k,spc_h2o2_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xh2o2(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_sogm_ndx > 0 ) then
! Subprogram not used           xsogm(:ncol,k)  = qin(:ncol,k,spc_sogm_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xsogm(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_sogi_ndx > 0 ) then
! Subprogram not used           xsogi(:ncol,k)  = qin(:ncol,k,spc_sogi_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xsogi(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_sogt_ndx > 0 ) then
! Subprogram not used           xsogt(:ncol,k)  = qin(:ncol,k,spc_sogt_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xsogt(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_sogb_ndx > 0 ) then
! Subprogram not used           xsogb(:ncol,k)  = qin(:ncol,k,spc_sogb_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xsogb(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_sogx_ndx > 0 ) then
! Subprogram not used           xsogx(:ncol,k)  = qin(:ncol,k,spc_sogx_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xsogx(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used        if( spc_so2_ndx > 0 ) then
! Subprogram not used           xso2(:ncol,k)  = qin(:ncol,k,spc_so2_ndx) * xhnm(:ncol,k)
! Subprogram not used        else
! Subprogram not used           xso2(:ncol,k)  = 0._r8
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     zsurf(:ncol) = m2km * phis(:ncol) * rga
! Subprogram not used     do k = ktop_all,pver-1
! Subprogram not used        delz(:ncol,k) = abs( (zmid(:ncol,k) - zmid(:ncol,k+1))*km2cm ) 
! Subprogram not used     end do
! Subprogram not used     delz(:ncol,pver) = abs( (zmid(:ncol,pver) - zsurf(:ncol) )*km2cm ) 
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !       ... part 0b,  for temperature dependent of henrys
! Subprogram not used     !                     xxhe1 = henry con for hno3
! Subprogram not used     !                     xxhe2 = henry con for h2o2
! Subprogram not used     !lwh 10/00 -- take henry''s law constants from brasseur et al. [1999],
! Subprogram not used     !             appendix j. for hno3, also consider dissociation to
! Subprogram not used     !             get effective henry''s law constant; equilibrium
! Subprogram not used     !             constant for dissociation from brasseur et al. [1999],
! Subprogram not used     !             appendix k. assume ph=5 (set as xph0 above).
! Subprogram not used     !             heff = h*k/[h+] for hno3 (complete dissociation)
! Subprogram not used     !             heff = h for h2o2 (no dissociation)
! Subprogram not used     !             heff = h * (1 + k/[h+]) (in general)
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     do k = ktop_all,pver
! Subprogram not used        work1(:ncol) = (t0 - tfld(:ncol,k))/(t0*tfld(:ncol,k))
! Subprogram not used        !-----------------------------------------------------------------
! Subprogram not used        ! 	... effective henry''s law constants:
! Subprogram not used        !	hno3, h2o2, ch2o, ch3ooh, ch3coooh (brasseur et al., 1999)
! Subprogram not used        !       xooh, onitr, macrooh (j.-f. muller; brocheton, 1999)
! Subprogram not used        !       isopooh (equal to hno3, as for macrooh)
! Subprogram not used        !       ho2no2 (mozart-1)
! Subprogram not used        !       ch3cocho, hoch2cho (betterton and hoffman, environ. sci. technol., 1988)
! Subprogram not used        !       ch3cho (staudinger and roberts, crit. rev. sci. technol., 1996)
! Subprogram not used        !       mvk, macr (allen et al., environ. toxicol. chem., 1998)
! Subprogram not used        !-----------------------------------------------------------------
! Subprogram not used        xk0(:)             = 2.1e5_r8 *exp( 8700._r8*work1(:) )
! Subprogram not used        xhen_hno3(:,k)     = xk0(:) * ( 1._r8 + hno3_diss / xph0 )
! Subprogram not used        xhen_h2o2(:,k)     = 7.45e4_r8 * exp( 6620._r8 * work1(:) )
! Subprogram not used        xhen_ch2o(:,k)     = 6.3e3_r8 * exp( 6460._r8 * work1(:) )
! Subprogram not used        xhen_ch3ooh(:,k)   = 2.27e2_r8 * exp( 5610._r8 * work1(:) )
! Subprogram not used        xhen_ch3co3h(:,k)  = 4.73e2_r8 * exp( 6170._r8 * work1(:) )
! Subprogram not used        xhen_ch3cocho(:,k) = 3.70e3_r8 * exp( 7275._r8 * work1(:) )
! Subprogram not used        xhen_xooh(:,k)     = 90.5_r8 * exp( 5607._r8 * work1(:) )
! Subprogram not used        xhen_onitr(:,k)    = 7.51e3_r8 * exp( 6485._r8 * work1(:) )
! Subprogram not used        xhen_ho2no2(:,k)   = 2.e4_r8
! Subprogram not used        xhen_glyald(:,k)   = 4.1e4_r8 * exp( 4600._r8 * work1(:) )
! Subprogram not used        xhen_ch3cho(:,k)   = 1.4e1_r8 * exp( 5600._r8 * work1(:) )
! Subprogram not used        xhen_mvk(:,k)      = 21._r8 * exp( 7800._r8 * work1(:) )
! Subprogram not used        xhen_macr(:,k)     = 4.3_r8 * exp( 5300._r8 * work1(:) )
! Subprogram not used        xhen_ch3cooh(:,k)  = 4.1e3_r8 * exp( 6300._r8 * work1(:) )
! Subprogram not used        xhen_sog(:,k)      = 5.e5_r8 * exp (12._r8 * work1(:) )
! Subprogram not used        !
! Subprogram not used        ! calculation for NH3 using the parameters in drydep_tables.F90
! Subprogram not used        !
! Subprogram not used        xhen_nh3 (:,k)     = 1.e6_r8
! Subprogram not used        xhen_ch3cn(:,k)     = 50._r8 * exp( 4000._r8 * work1(:) )
! Subprogram not used        xhen_hcn(:,k)       = 12._r8 * exp( 5000._r8 * work1(:) )
! Subprogram not used        do i = 1, ncol
! Subprogram not used           so2_diss        = 1.23e-2_r8 * exp( 1960._r8 * work1(i) )
! Subprogram not used           xhen_so2(i,k)   = 1.23_r8 * exp( 3120._r8 * work1(i) ) * ( 1._r8 + so2_diss / xph0 )
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        tmp_hetrates(:,k,:) = 0._r8
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !       ... part 1, solve for high henry constant ( hno3, h2o2)
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     col_loop :  do i = 1,ncol
! Subprogram not used        xgas1(:) = xhno3(i,:)                     ! xgas will change during 
! Subprogram not used        xgas2(:) = xh2o2(i,:)                     ! different levels wash 
! Subprogram not used        xgas3(:) = xso2 (i,:)
! Subprogram not used        xgas4(:) = xsogm(i,:)
! Subprogram not used        xgas5(:) = xsogi(i,:)
! Subprogram not used        xgas6(:) = xsogt(i,:)
! Subprogram not used        xgas7(:) = xsogb(i,:)
! Subprogram not used        xgas8(:) = xsogx(i,:)
! Subprogram not used        level_loop1  : do kk = ktop(i),pver
! Subprogram not used           stay = 1._r8
! Subprogram not used           if( rain(i,kk) /= 0._r8 ) then            ! finding rain cloud           
! Subprogram not used              all1 = 0._r8                           ! accumulation to justisfy saturation
! Subprogram not used              all2 = 0._r8 
! Subprogram not used              all3 = 0._r8 
! Subprogram not used              all4 = 0._r8 
! Subprogram not used              all5 = 0._r8 
! Subprogram not used              all6 = 0._r8 
! Subprogram not used              all7 = 0._r8 
! Subprogram not used              all8 = 0._r8 
! Subprogram not used              stay = ((zmid(i,kk) - zsurf(i))*km2cm)/(xum*delt)
! Subprogram not used              stay = min( stay,1._r8 )
! Subprogram not used              !-----------------------------------------------------------------
! Subprogram not used              !       ... calculate the saturation concentration eqca
! Subprogram not used              !-----------------------------------------------------------------
! Subprogram not used              do k = kk,pver                      ! cal washout below cloud
! Subprogram not used                 xeqca1 =  xgas1(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_hno3(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca2 =  xgas2(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_h2o2(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca3 =  xgas3(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_so2( i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca4 =  xgas4(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_sog(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca5 =  xgas5(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_sog(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca6 =  xgas6(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_sog(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca7 =  xgas7(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_sog(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used                 xeqca8 =  xgas8(k) &
! Subprogram not used                      / (xliq(i,kk)*avo2 + 1._r8/(xhen_sog(i,k)*const0*tfld(i,k))) &
! Subprogram not used                      *  xliq(i,kk)*avo2
! Subprogram not used 
! Subprogram not used                 !-----------------------------------------------------------------
! Subprogram not used                 !       ... calculate ca; inside cloud concentration in #/cm3(air)
! Subprogram not used                 !-----------------------------------------------------------------
! Subprogram not used                 xca1 = geo_fac*xkgm*xgas1(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca2 = geo_fac*xkgm*xgas2(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca3 = geo_fac*xkgm*xgas3(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca4 = geo_fac*xkgm*xgas4(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca5 = geo_fac*xkgm*xgas5(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca6 = geo_fac*xkgm*xgas6(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca7 = geo_fac*xkgm*xgas7(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used                 xca8 = geo_fac*xkgm*xgas8(k)/(xrm*xum)*delz(i,k) * xliq(i,kk) * cm3_2_m3
! Subprogram not used 
! Subprogram not used                 !-----------------------------------------------------------------
! Subprogram not used                 !       ... if is not saturated
! Subprogram not used                 !               hno3(gas)_new = hno3(gas)_old - hno3(h2o)
! Subprogram not used                 !           otherwise
! Subprogram not used                 !               hno3(gas)_new = hno3(gas)_old
! Subprogram not used                 !-----------------------------------------------------------------
! Subprogram not used                 all1 = all1 + xca1
! Subprogram not used                 all2 = all2 + xca2
! Subprogram not used                 if( all1 < xeqca1 ) then
! Subprogram not used                    xgas1(k) = max( xgas1(k) - xca1,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 if( all2 < xeqca2 ) then
! Subprogram not used                    xgas2(k) = max( xgas2(k) - xca2,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 all3 = all3 + xca3
! Subprogram not used                 if( all3 < xeqca3 ) then
! Subprogram not used                    xgas3(k) = max( xgas3(k) - xca3,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 all4 = all4 + xca4
! Subprogram not used                 all5 = all5 + xca5
! Subprogram not used                 all6 = all6 + xca6
! Subprogram not used                 all7 = all7 + xca7
! Subprogram not used                 all8 = all8 + xca8
! Subprogram not used                 if( all4 < xeqca4 ) then
! Subprogram not used                    xgas4(k) = max( xgas4(k) - xca4,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 if( all5 < xeqca5 ) then
! Subprogram not used                    xgas5(k) = max( xgas5(k) - xca5,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 if( all6 < xeqca6 ) then
! Subprogram not used                    xgas6(k) = max( xgas6(k) - xca6,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 if( all7 < xeqca7 ) then
! Subprogram not used                    xgas7(k) = max( xgas7(k) - xca7,0._r8 )
! Subprogram not used                 end if
! Subprogram not used                 if( all8 < xeqca8 ) then
! Subprogram not used                    xgas8(k) = max( xgas8(k) - xca8,0._r8 )
! Subprogram not used                 end if
! Subprogram not used              end do
! Subprogram not used           end if
! Subprogram not used           !-----------------------------------------------------------------
! Subprogram not used           !       ... calculate the lifetime of washout (second)
! Subprogram not used           !             after all layers washout 
! Subprogram not used           !             the concentration of hno3 is reduced 
! Subprogram not used           !             then the lifetime xtt is calculated by
! Subprogram not used           !
! Subprogram not used           !                  xtt = (xhno3(ini) - xgas1(new))/(dt*xhno3(ini))
! Subprogram not used           !                  where dt = passing time (s) in vertical
! Subprogram not used           !                             path below the cloud
! Subprogram not used           !                        dt = dz(cm)/um(cm/s)
! Subprogram not used           !-----------------------------------------------------------------
! Subprogram not used           xdtm = delz(i,kk) / xum                     ! the traveling time in each dz
! Subprogram not used           xxx1 = (xhno3(i,kk) - xgas1(kk))
! Subprogram not used           xxx2 = (xh2o2(i,kk) - xgas2(kk))
! Subprogram not used           if( xxx1 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              yhno3  = xhno3(i,kk)/xxx1 * xdtm    
! Subprogram not used           else
! Subprogram not used              yhno3  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           if( xxx2 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              yh2o2  = xh2o2(i,kk)/xxx2 * xdtm     
! Subprogram not used           else
! Subprogram not used              yh2o2  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           tmp_hetrates(i,kk,1) = max( 1._r8 / yh2o2,0._r8 ) * stay
! Subprogram not used           tmp_hetrates(i,kk,2) = max( 1._r8 / yhno3,0._r8 ) * stay
! Subprogram not used           xxx3 = (xso2( i,kk) - xgas3(kk))
! Subprogram not used           if( xxx3 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              yso2  = xso2( i,kk)/xxx3 * xdtm     
! Subprogram not used           else
! Subprogram not used              yso2  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           tmp_hetrates(i,kk,3) = max( 1._r8 / yso2, 0._r8 ) * stay
! Subprogram not used           xxx4 = (xsogm(i,kk) - xgas4(kk))
! Subprogram not used           xxx5 = (xsogi(i,kk) - xgas5(kk))
! Subprogram not used           xxx6 = (xsogt(i,kk) - xgas6(kk))
! Subprogram not used           xxx7 = (xsogb(i,kk) - xgas7(kk))
! Subprogram not used           xxx8 = (xsogx(i,kk) - xgas8(kk))
! Subprogram not used           if( xxx4 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              ysogm  = xsogm(i,kk)/xxx4 * xdtm
! Subprogram not used           else
! Subprogram not used              ysogm  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           if( xxx5 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              ysogi  = xsogi(i,kk)/xxx5 * xdtm
! Subprogram not used           else
! Subprogram not used              ysogi  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           if( xxx6 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              ysogt  = xsogt(i,kk)/xxx6 * xdtm
! Subprogram not used           else
! Subprogram not used              ysogt  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           if( xxx7 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              ysogb  = xsogb(i,kk)/xxx7 * xdtm
! Subprogram not used           else
! Subprogram not used              ysogb  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           if( xxx8 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
! Subprogram not used              ysogx  = xsogx(i,kk)/xxx8 * xdtm
! Subprogram not used           else
! Subprogram not used              ysogx  = 1.e29_r8
! Subprogram not used           end if
! Subprogram not used           tmp_hetrates(i,kk,4) = max( 1._r8 / ysogm,0._r8 ) * stay
! Subprogram not used           tmp_hetrates(i,kk,5) = max( 1._r8 / ysogi,0._r8 ) * stay
! Subprogram not used           tmp_hetrates(i,kk,6) = max( 1._r8 / ysogt,0._r8 ) * stay
! Subprogram not used           tmp_hetrates(i,kk,7) = max( 1._r8 / ysogb,0._r8 ) * stay
! Subprogram not used           tmp_hetrates(i,kk,8) = max( 1._r8 / ysogx,0._r8 ) * stay
! Subprogram not used        end do level_loop1
! Subprogram not used     end do col_loop
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !       ... part 2, in-cloud solve for low henry constant
! Subprogram not used     !                   hno3 and h2o2 have both in and under cloud
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     level_loop2 : do k = ktop_all,pver
! Subprogram not used        Column_loop2 : do i=1,ncol
! Subprogram not used           if ( rain(i,k) <= 0._r8 ) then
! Subprogram not used              het_rates(i,k,:) =  0._r8 
! Subprogram not used              cycle
! Subprogram not used           endif
! Subprogram not used 
! Subprogram not used           work1(i) = avo2 * xliq(i,k)
! Subprogram not used           work2(i) = const0 * tfld(i,k)
! Subprogram not used           work3(i) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch2o(i,k)*work2(i)))),0._r8 ) &
! Subprogram not used                * satf_ch2o
! Subprogram not used           if( ch2o_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch2o_ndx)  = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( isopno3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,isopno3_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( xisopno3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,xisopno3_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( hyac_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,hyac_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( hydrald_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,hydrald_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           work3(i) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3ooh(i,k)*work2(i)))),0._r8 )
! Subprogram not used           if( ch3ooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3ooh_ndx)  = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( pooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,pooh_ndx)  = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( c2h5ooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,c2h5ooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( c3h7ooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,c3h7ooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( rooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,rooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( ch3oh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3oh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( c2h5oh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,c2h5oh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( alkooh_ndx  > 0 ) then
! Subprogram not used              het_rates(i,k,alkooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( mekooh_ndx  > 0 ) then
! Subprogram not used              het_rates(i,k,mekooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( tolooh_ndx  > 0 ) then
! Subprogram not used              het_rates(i,k,tolooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( terpooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,terpooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           if( ch3coooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3coooh_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3co3h(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( ho2no2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ho2no2_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ho2no2(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( xho2no2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,xho2no2_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ho2no2(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( ch3cocho_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3cocho_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3cocho(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( xooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,xooh_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_xooh(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( onitr_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,onitr_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_onitr(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( xonitr_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,xonitr_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_onitr(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( glyald_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,glyald_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_glyald(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( ch3cho_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3cho_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3cho(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( mvk_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,mvk_ndx)  = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_mvk(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( macr_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,macr_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_macr(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( h2o2_ndx > 0 ) then
! Subprogram not used              work3(i) = satf_h2o2 * max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_h2o2(i,k)*work2(i)))),0._r8 )    
! Subprogram not used              het_rates(i,k,h2o2_ndx) =  work3(i) + tmp_hetrates(i,k,1)
! Subprogram not used           end if
! Subprogram not used           if ( prog_modal_aero .and. so2_ndx>0 .and. h2o2_ndx>0 ) then
! Subprogram not used              het_rates(i,k,so2_ndx) = het_rates(i,k,h2o2_ndx)
! Subprogram not used           elseif( so2_ndx > 0 ) then
! Subprogram not used              work3(i) = satf_so2 * max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_so2( i,k)*work2(i)))),0._r8 )    
! Subprogram not used              het_rates(i,k,so2_ndx ) =  work3(i) + tmp_hetrates(i,k,3)
! Subprogram not used           endif
! Subprogram not used !
! Subprogram not used           work3(i) = satf_sog * max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_sog(i,k)*work2(i)))),0._r8 )
! Subprogram not used           if( sogm_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sogm_ndx) =  work3(i) + tmp_hetrates(i,k,4)
! Subprogram not used           end if
! Subprogram not used           if( sogi_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sogi_ndx) =  work3(i) + tmp_hetrates(i,k,5)
! Subprogram not used           end if
! Subprogram not used           if( sogt_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sogt_ndx) =  work3(i) + tmp_hetrates(i,k,6)
! Subprogram not used           end if
! Subprogram not used           if( sogb_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sogb_ndx) =  work3(i) + tmp_hetrates(i,k,7)
! Subprogram not used           end if
! Subprogram not used           if( sogx_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sogx_ndx) =  work3(i) + tmp_hetrates(i,k,8)
! Subprogram not used           end if
! Subprogram not used !
! Subprogram not used           work3(i) = tmp_hetrates(i,k,2) + satf_hno3 * &
! Subprogram not used                max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_hno3(i,k)*work2(i)))),0._r8 )    
! Subprogram not used           tmp0_rates(i)   = work3(i)
! Subprogram not used           tmp1_rates(i)   = .2_r8*work3(i)
! Subprogram not used           if( hno3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,hno3_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( xhno3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,xhno3_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( onit_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,onit_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( xonit_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,xonit_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( Pb_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,Pb_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( macrooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,macrooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( isopooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,isopooh_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           if( clono2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, clono2_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( brono2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, brono2_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( hcl_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, hcl_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( n2o5_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, n2o5_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( hocl_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, hocl_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( hobr_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, hobr_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used           if( hbr_ndx > 0 ) then
! Subprogram not used              het_rates(i,k, hbr_ndx) = work3(i)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           if( soa_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,soa_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( oc2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,oc2_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( cb2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,cb2_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( so4_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,so4_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( sa1_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sa1_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( sa2_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sa2_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( sa3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sa3_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( sa4_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,sa4_ndx) = tmp1_rates(i)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           if( h2so4_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,h2so4_ndx) = tmp0_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( nh4_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,nh4_ndx) = tmp0_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( nh4no3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,nh4no3_ndx ) = tmp0_rates(i)
! Subprogram not used           end if
! Subprogram not used           if( nh3_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,nh3_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_nh3(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           if( ch3cooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3cooh_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3cooh(i,k)*work2(i)))),0._r8 )
! Subprogram not used           end if
! Subprogram not used           if( hcooh_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,hcooh_ndx) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3cooh(i,k)*work2(i)))),0._r8 )
! Subprogram not used           endif
! Subprogram not used           if ( hcn_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,hcn_ndx     ) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_hcn(i,k)*work2(i)))),0._r8 )
! Subprogram not used           endif
! Subprogram not used           if ( ch3cn_ndx > 0 ) then
! Subprogram not used              het_rates(i,k,ch3cn_ndx   ) = max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_ch3cn(i,k)*work2(i)))),0._r8 )
! Subprogram not used           endif
! Subprogram not used 
! Subprogram not used        end do Column_loop2
! Subprogram not used     end do level_loop2
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     !	... Set rates above tropopause = 0.
! Subprogram not used     !-----------------------------------------------------------------
! Subprogram not used     do mm = 1,gas_wetdep_cnt
! Subprogram not used        m = wetdep_map(mm)
! Subprogram not used        do i = 1,ncol
! Subprogram not used           do k = 1,ktop(i)
! Subprogram not used              het_rates(i,k,m) = 0._r8
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        if ( any( het_rates(:ncol,:,m) == MISSING) ) then
! Subprogram not used           write(hetratestrg,'(I3)') m
! Subprogram not used           call endrun('sethet: het_rates (wet dep) not set for het reaction number : '//hetratestrg)
! Subprogram not used        endif
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine sethet

end module mo_sethet
