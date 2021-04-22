
      module mo_strato_rates
!=======================================================================
! ROUTINE
!   ratecon_sfstrat.f
!
!  Date...
!  15 August 2002
!  11 April  2008
!
!  Programmed by...
!   Douglas E. Kinnison
!
! DESCRIPTION
!
! Derivation of the rate constant for reactions on
!   sulfate, NAT, and ICE aerosols.
!
!
! Sulfate Aerosol Reactions               Rxn#   Gamma
!   N2O5   + H2O(l)     =>  2HNO3         (1)    f(wt%)
!   ClONO2 + H2O(l)     =>  HOCl + HNO3   (2)    f(T,P,HCl,H2O,r)
!   BrONO2 + H2O(l)     =>  HOBr + HNO3   (3)    f(T,P,H2O,r)
!   ClONO2 + HCl(l)     =>  Cl2  + HNO3   (4)    f(T,P,HCl,H2O,r)
!   HOCl   + HCl(l)     =>  Cl2  + H2O    (5)    f(T,P,HCl,HOCl,H2O,r)
!   HOBr   + HCl(l)     =>  BrCl + H2O    (6)    f(T,P,HCl,HOBr,H2O,r)
!
! Nitric Acid Di-hydrate Reactions        Rxn#    Gamma   Reference
!   N2O5   + H2O(s)     =>  2HNO3         (7)     4e-4   JPL06-2
!   ClONO2 + H2O(s)     =>  HOCl + HNO3   (8)     4e-3   JPL06-2
!   ClONO2 + HCl(s)     =>  Cl2  + HNO3   (9)     0.2    JPL06-2
!   HOCl   + HCl(s)     =>  Cl2  + H2O    (10)    0.1    JPL06-2
!   BrONO2 + H2O(s)     =>  HOBr + HNO3   (11)    0.3    David Hanson PC
!
! ICE Aersol Reactions                    Rxn#    Gamma
!   N2O5   + H2O(s)     =>  2HNO3         (12)     0.02   JPL06-2
!   ClONO2 + H2O(s)     =>  HOCl + HNO3   (13)     0.3    JPL06-2
!   BrONO2 + H2O(s)     =>  HOBr + HNO3   (14)     0.3    JPL06-2
!   ClONO2 + HCl(s)     =>  Cl2  + HNO3   (15)     0.3    JPL06-2
!   HOCl   + HCl(s)     =>  Cl2  + H2O    (16)     0.2    JPL06-2
!   HOBr   + HCl(s)     =>  BrCl + H2O    (17)     0.3    JPL06-2
!
! NOTE: The rate constants derived from species reacting with H2O are
!       first order (i.e., sec-1 units) - an example is N2O5 + H2O = 2HNO3.
!       Other reactions, e.g., ClONO2 + HCl have rate constants that
!       are second order (i.e., cm+3 molecules-1 sec-1 units). In all
!       of these types of reactions the derived first order rate constant
!       {0.25*(mean Velocity)*SAD*gamma} is divided by the HCl abundance
!       to derive the correct second order units.
!
! NOTE: Liquid Sulfate Aerosols...
!       See coding for references on how the Sulfate Aerosols were handled.
!       Data was used that was more recent than JPL00.
!
!
! INPUT:
!  ad      .    .... air density, molec. cm-3
!  pmid        ..... pressures, hPa
!  temp        ..... temperatures, K
!  rad_sulfate ..... Surface area density, cm2 cm-3
!  sad_sulfate ..... Surface area density, cm2 cm-3
!  sad_nat     ..... Surface area density, cm2 cm-3
!  sad_ice     ..... Surface area density, cm2 cm-3
!  brono2mv    ..... BrONO2 Volume Mixing Ratio
!  clono2mvr   ..... ClONO2 Volume Mixing Ratio
!  h2omvr      ..... H2O Volume Mixing Ratio
!  hclmvr      ..... HCl Volume Mixing Ratio
!  hobrmvr     ..... HOBr Volume Mixing Ratio
!  hoclmvr     ..... HOCl Volume Mixing Ratio
!  n2o5mvr     ..... N2O5 Volume Mixing Ratio
!
! OUTPUT:
!
!  rxt         ..... Rate constant (s-1 and cm3 sec-1 molec-1)
!=======================================================================

      private
      public :: ratecon_sfstrat, init_strato_rates, has_strato_chem

      integer :: id_brono2, id_clono2, id_hcl, id_hocl, &
           id_hobr, id_n2o5
      integer :: rid_het1,  rid_het2,  rid_het3,  rid_het4,  rid_het5, &
           rid_het6,  rid_het7,  rid_het8,  rid_het9,  rid_het10, &
           rid_het11, rid_het12, rid_het13, rid_het14, rid_het15, &
           rid_het16, rid_het17

      logical :: has_strato_chem 

      contains

        subroutine init_strato_rates

          use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
          use mo_aero_settling, only: strat_aer_settl_init
          implicit none

          integer :: ids(23)

          rid_het1  = get_rxt_ndx( 'het1' )
          rid_het2  = get_rxt_ndx( 'het2' )
          rid_het3  = get_rxt_ndx( 'het3' )
          rid_het4  = get_rxt_ndx( 'het4' )
          rid_het5  = get_rxt_ndx( 'het5' )
          rid_het6  = get_rxt_ndx( 'het6' )
          rid_het7  = get_rxt_ndx( 'het7' )
          rid_het8  = get_rxt_ndx( 'het8' )
          rid_het9  = get_rxt_ndx( 'het9' )
          rid_het10 = get_rxt_ndx( 'het10' )
          rid_het11 = get_rxt_ndx( 'het11' )
          rid_het12 = get_rxt_ndx( 'het12' )
          rid_het13 = get_rxt_ndx( 'het13' )
          rid_het14 = get_rxt_ndx( 'het14' )
          rid_het15 = get_rxt_ndx( 'het15' )
          rid_het16 = get_rxt_ndx( 'het16' )
          rid_het17 = get_rxt_ndx( 'het17' )

          id_brono2 = get_spc_ndx( 'BRONO2' )
          id_clono2 = get_spc_ndx( 'CLONO2' )
          id_hcl    = get_spc_ndx( 'HCL' )
          id_hocl   = get_spc_ndx( 'HOCL' )
          id_hobr   = get_spc_ndx( 'HOBR' )
          id_n2o5   = get_spc_ndx( 'N2O5' )

          ids(:) = (/ rid_het1, rid_het2, rid_het3, rid_het4, rid_het5, rid_het6, rid_het7, rid_het8, &
               rid_het9, rid_het10, rid_het11, rid_het12, rid_het13, rid_het14, rid_het15, &
               rid_het16, rid_het17, id_brono2, id_clono2, id_hcl, id_hocl, id_hobr, id_n2o5 /)

          has_strato_chem = all( ids(:) > 0 )

          if (.not. has_strato_chem) return

          call strat_aer_settl_init

        endsubroutine init_strato_rates

! Subprogram not used       subroutine ratecon_sfstrat( ad, pmid, temp, rad_sulfate, sad_sulfate, &
! Subprogram not used                                   sad_nat, sad_ice, h2ovmr, vmr, rxt, ncol )
! Subprogram not used 
! Subprogram not used       use shr_kind_mod, only : r8 => shr_kind_r8
! Subprogram not used       use chem_mods,    only : adv_mass, rxntot, gas_pcnst
! Subprogram not used       use ppgrid,       only : pcols, pver
! Subprogram not used       use mo_sad,       only : sad_top
! Subprogram not used       use cam_logfile,  only : iulog
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer, intent(in) :: ncol                                 ! columns in chunk
! Subprogram not used       real(r8), dimension(ncol,pver,gas_pcnst), intent(in) :: &   ! species concentrations (mol/mol)
! Subprogram not used         vmr
! Subprogram not used       real(r8), dimension(ncol,pver), intent(in) :: &
! Subprogram not used         ad, &                                                     ! Air Density (molec. cm-3)
! Subprogram not used         rad_sulfate, &                                            ! Radius of Sulfate Aerosol (cm)
! Subprogram not used         sad_ice, &                                                ! ICE Surface Area Density (cm-1)
! Subprogram not used         sad_nat, &                                                ! NAT Surface Area Density (cm-1)
! Subprogram not used         sad_sulfate, &                                            ! Sulfate Surface Area Density (cm-1)
! Subprogram not used         h2ovmr                                                    ! water vapor volume mixing ratio( gas phase )
! Subprogram not used       real(r8), dimension(pcols,pver), intent(in) :: &
! Subprogram not used         pmid, &                                                   ! pressure (Pa)
! Subprogram not used         temp                                                      ! temperature (K)
! Subprogram not used 
! Subprogram not used       real(r8), intent(out) :: &
! Subprogram not used         rxt(ncol,pver,rxntot)                                     ! rate constants
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !  	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 	real(r8), parameter :: small_conc = 1.e-30_r8
! Subprogram not used 	real(r8), parameter :: av_const   = 2.117265e4_r8  ! (8*8.31448*1000 / PI)
! Subprogram not used 	real(r8), parameter :: pa2mb      = 1.e-2_r8       ! Pa to mb
! Subprogram not used 	real(r8), parameter :: m2cm       = 100._r8        ! meters to cms
! Subprogram not used 
! Subprogram not used       integer :: &
! Subprogram not used         i, &                      ! altitude loop index
! Subprogram not used         k, &                      ! level loop index
! Subprogram not used         m                         ! species index
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !   	... variables for gamma calculations
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8) :: &
! Subprogram not used         brono2vmr, &                            ! BrONO2 Volume Mixing Ratio
! Subprogram not used         clono2vmr, &                            ! ClONO2 Volume Mixing Ratio
! Subprogram not used         hclvmr, &                               ! HCl Volume Mixing Ratio
! Subprogram not used         hcldeni, &                              ! inverse of HCl density
! Subprogram not used         cntdeni, &                              ! inverse of ClONO2 density
! Subprogram not used         hocldeni, &                             ! inverse of HOCl density
! Subprogram not used         hobrdeni, &                             ! inverse of HOBr density
! Subprogram not used         hoclvmr, &                              ! HOCl Volume Mixing Ratio
! Subprogram not used         hobrvmr, &                              ! HOBr Volume Mixing Ratio
! Subprogram not used         n2o5vmr                                 ! N2O5 Volume Mixing Ratio
! Subprogram not used 
! Subprogram not used         real(r8) :: &
! Subprogram not used         av_n2o5, &                              ! N2O5 Mean Velocity (cm s-1)
! Subprogram not used         av_clono2, &                            ! ClONO2 Mean Velocity (cm s-1)
! Subprogram not used         av_brono2, &                            ! BrONO2Mean Velocity (cm s-1)
! Subprogram not used         av_hocl, &                              ! HOCl Mean Velocity (cm s-1)
! Subprogram not used         av_hobr                                 ! HOBr Mean Velocity (cm s-1)
! Subprogram not used 
! Subprogram not used       real(r8) :: &
! Subprogram not used         pzero_h2o, &                            ! H2O sat vapor press (mbar)
! Subprogram not used         e0, e1, e2, e3, &                       ! coefficients for H2O sat vapor press.
! Subprogram not used         aw, &                                   ! Water activity
! Subprogram not used         m_h2so4, &                              ! H2SO4 molality (mol/kg)
! Subprogram not used         wt, &                                   ! wt % H2SO4
! Subprogram not used         y1, y2, &                               ! used in H2SO4 molality
! Subprogram not used         &  a1, b1, c1, d1, &                    ! used in H2SO4 molality
! Subprogram not used         a2, b2, c2, d2                          ! used in H2SO4 molality
! Subprogram not used 
! Subprogram not used         real(r8) :: &
! Subprogram not used         z1, z2, z3, &                           ! used in H2SO4 soln density
! Subprogram not used         den_h2so4, &                            ! H2SO4 soln density, g/cm3
! Subprogram not used         mol_h2so4, &                            ! Molality of H2SO4, mol / kg
! Subprogram not used         molar_h2so4, &                          ! Molarity of H2SO4, mol / l
! Subprogram not used         x_h2so4, &                              ! H2SO4 mole fraction
! Subprogram not used         aconst, tzero, &                        ! used in viscosity of H2SO4
! Subprogram not used         vis_h2so4, &                            ! H2SO4 viscosity
! Subprogram not used         ah, &                                   ! Acid activity, molarity units
! Subprogram not used         term1,term2,term3,term4, &              ! used in ah
! Subprogram not used         term5,term6,term7,term0, &
! Subprogram not used         T_limit, &                              ! temporary variable for temp (185-260K range)
! Subprogram not used         T_limiti, &                             ! 1./T_limit
! Subprogram not used         T_limitsq, &                            ! sqrt( T_limit )
! Subprogram not used         rad_sulf, &                             ! temporary variable for sulfate radius (cm)
! Subprogram not used         sadsulf, &                              ! temporary variable for sulfate radius (cm)
! Subprogram not used         sadice, &                               ! temporary variable for sulfate radius (cm)
! Subprogram not used         sadnat                                  ! temporary variable for sulfate radius (cm)
! Subprogram not used 
! Subprogram not used       real(r8) :: &
! Subprogram not used         C_cnt, S_cnt, &                         ! used in H_cnt
! Subprogram not used         H_cnt, &                                ! Henry's law coeff. for ClONO2
! Subprogram not used         H_hcl, &                                ! Henry's law coeff. for HCl
! Subprogram not used         D_cnt, &
! Subprogram not used         k_hydr, &
! Subprogram not used         k_h2o, &
! Subprogram not used         k_h, &
! Subprogram not used         k_hcl, &
! Subprogram not used         rdl_cnt, &
! Subprogram not used         f_cnt, &
! Subprogram not used         M_hcl, &
! Subprogram not used         atmos
! Subprogram not used 
! Subprogram not used       real(r8) :: &
! Subprogram not used         Gamma_b_h2o, &
! Subprogram not used         Gamma_cnt_rxn, &
! Subprogram not used         Gamma_b_hcl, &
! Subprogram not used         Gamma_s, &
! Subprogram not used         Fhcl, &
! Subprogram not used         Gamma_s_prime, &
! Subprogram not used         Gamma_b_hcl_prime, &
! Subprogram not used         Gamma_b, &
! Subprogram not used         gprob_n2o5, &
! Subprogram not used         gprob_rxn, &
! Subprogram not used         gprob_tot, &
! Subprogram not used         gprob_cnt, &
! Subprogram not used         gprob_cnt_hcl, &
! Subprogram not used         gprob_cnt_h2o
! Subprogram not used 
! Subprogram not used         real(r8) :: &
! Subprogram not used         D_hocl, &
! Subprogram not used         k_hocl_hcl, &
! Subprogram not used         C_hocl, &
! Subprogram not used         S_hocl, &
! Subprogram not used         H_hocl, &
! Subprogram not used         Gamma_hocl_rxn, &
! Subprogram not used         rdl_hocl, &
! Subprogram not used         f_hocl, &
! Subprogram not used         gprob_hocl_hcl
! Subprogram not used 
! Subprogram not used         real(r8) :: &
! Subprogram not used         h1, h2, h3, &
! Subprogram not used         alpha, &
! Subprogram not used         gprob_bnt_h2o
! Subprogram not used 
! Subprogram not used       real(r8) :: &
! Subprogram not used         C_hobr, &
! Subprogram not used         D_hobr, &
! Subprogram not used         aa, bb, cc, dd, &
! Subprogram not used         k_hobr_hcl, &
! Subprogram not used         k_dl, &
! Subprogram not used         k_wasch, &
! Subprogram not used         H_hobr, &
! Subprogram not used         rdl_hobr, &
! Subprogram not used         Gamma_hobr_rxn, &
! Subprogram not used         f_hobr, &
! Subprogram not used         gprob_hobr_hcl
! Subprogram not used 
! Subprogram not used       real(r8) :: &
! Subprogram not used         pmb,&					! Pressure, mbar (hPa)
! Subprogram not used         pH2O_atm,&				! Partial press. H2O (atm)
! Subprogram not used         pH2O_hPa,&				! Partial press. H2O (hPa)
! Subprogram not used         pHCl_atm,&				! Partial press. HCl (atm)
! Subprogram not used         pCNT_atm                                ! Partial press. ClONO2 (atm)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... Used in pzero h2o calculation
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: wt_e0 = 18.452406985_r8
! Subprogram not used       real(r8), parameter :: wt_e1 = 3505.1578807_r8
! Subprogram not used       real(r8), parameter :: wt_e2 = 330918.55082_r8
! Subprogram not used       real(r8), parameter :: wt_e3 = 12725068.262_r8
! Subprogram not used 
! Subprogram not used       real(r8) :: &
! Subprogram not used         wrk, tmp
! Subprogram not used 
! Subprogram not used       real(r8), parameter :: small = 1.e-16_r8
! Subprogram not used 
! Subprogram not used       if (.not. has_strato_chem) return
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... intialize rate constants
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       do k = 1,pver
! Subprogram not used          rxt(:,k,rid_het1) = 0._r8
! Subprogram not used          rxt(:,k,rid_het2) = 0._r8
! Subprogram not used          rxt(:,k,rid_het3) = 0._r8
! Subprogram not used          rxt(:,k,rid_het4) = 0._r8
! Subprogram not used          rxt(:,k,rid_het5) = 0._r8
! Subprogram not used          rxt(:,k,rid_het6) = 0._r8
! Subprogram not used          rxt(:,k,rid_het7) = 0._r8
! Subprogram not used          rxt(:,k,rid_het8) = 0._r8
! Subprogram not used          rxt(:,k,rid_het9) = 0._r8
! Subprogram not used          rxt(:,k,rid_het10) = 0._r8
! Subprogram not used          rxt(:,k,rid_het11) = 0._r8
! Subprogram not used          rxt(:,k,rid_het12) = 0._r8
! Subprogram not used          rxt(:,k,rid_het13) = 0._r8
! Subprogram not used          rxt(:,k,rid_het14) = 0._r8
! Subprogram not used          rxt(:,k,rid_het15) = 0._r8
! Subprogram not used          rxt(:,k,rid_het16) = 0._r8
! Subprogram not used          rxt(:,k,rid_het17) = 0._r8
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... set rate constants
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used Level_loop : &
! Subprogram not used       do k = sad_top+1,pver
! Subprogram not used column_loop : &
! Subprogram not used          do i = 1,ncol
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !	... set species, pmb, and atmos
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 	    brono2vmr = vmr(i,k,id_brono2)
! Subprogram not used 	    clono2vmr = vmr(i,k,id_clono2)
! Subprogram not used 	    hclvmr    = vmr(i,k,id_hcl)
! Subprogram not used 	    hoclvmr      = vmr(i,k,id_hocl)
! Subprogram not used 	    hobrvmr      = vmr(i,k,id_hobr)
! Subprogram not used 	    if( hclvmr > 0._r8 ) then
! Subprogram not used 	       hcldeni  = 1._r8/(hclvmr*ad(i,k))
! Subprogram not used 	    end if
! Subprogram not used 	    if( clono2vmr > 0._r8 ) then
! Subprogram not used 	       cntdeni  = 1._r8/(clono2vmr*ad(i,k))
! Subprogram not used 	    end if
! Subprogram not used 	    if( hoclvmr > 0._r8 ) then
! Subprogram not used 	       hocldeni  = 1._r8/(hoclvmr*ad(i,k))
! Subprogram not used 	    end if
! Subprogram not used 	    if( hobrvmr > 0._r8 ) then
! Subprogram not used 	       hobrdeni  = 1._r8/(hobrvmr*ad(i,k))
! Subprogram not used 	    end if
! Subprogram not used 	    n2o5vmr      = vmr(i,k,id_n2o5)
! Subprogram not used             sadsulf      = sad_sulfate(i,k)
! Subprogram not used             sadnat       = sad_nat(i,k)
! Subprogram not used 	    sadice       = sad_ice(i,k)
! Subprogram not used             pmb          = pa2mb*pmid(i,k)
! Subprogram not used             atmos        = pmb/1013.25_r8
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !  	... setup for stratospheric aerosols
! Subprogram not used !           data range set: 185K - 240K;    GRL, 24, 1931, 1997
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used             T_limit   = max( temp(i,k),185._r8 )
! Subprogram not used             T_limit   = min( T_limit,240._r8 )
! Subprogram not used             T_limiti  = 1._r8/T_limit
! Subprogram not used             T_limitsq = sqrt( T_limit )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     .... Average velocity (8RT*1000/(PI*MW))**1/2 * 100.(units cm s-1)
! Subprogram not used !     .... or (av_const*T/M2)**1/2
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 	    wrk       = av_const*T_limit
! Subprogram not used             av_n2o5   = sqrt( wrk/adv_mass(id_n2o5) )*m2cm
! Subprogram not used             av_clono2 = sqrt( wrk/adv_mass(id_clono2) )*m2cm
! Subprogram not used             av_brono2 = sqrt( wrk/adv_mass(id_brono2) )*m2cm
! Subprogram not used             av_hocl   = sqrt( wrk/adv_mass(id_hocl) )*m2cm
! Subprogram not used             av_hobr   = sqrt( wrk/adv_mass(id_hobr) )*m2cm
! Subprogram not used has_sadsulf : &
! Subprogram not used             if( sadsulf > 0._r8 ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     .... Partial Pressure of H2O, ClONO2, and HCl in atmospheres
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( hclvmr > 0._r8 ) then
! Subprogram not used                   pHCl_atm  = hclvmr*atmos
! Subprogram not used                else
! Subprogram not used                   pHCl_atm  = 0._r8
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used                if( clono2vmr > 0._r8 ) then
! Subprogram not used                   pCNT_atm  = clono2vmr*atmos
! Subprogram not used                else
! Subprogram not used                   pCNT_atm  = 0._r8
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used                if( h2ovmr(i,k) > 0._r8 ) then
! Subprogram not used                   pH2O_atm  = h2ovmr(i,k)*atmos
! Subprogram not used                else
! Subprogram not used                   pH2O_atm  = 0._r8
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     .... Partial Pressure of H2O in hPa
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                pH2O_hpa = h2ovmr(i,k)*pmb
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the h2so4 Wt% and Activity of H2O - 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... Saturation Water Vapor Pressure (mbar)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                pzero_h2o = exp( wt_e0 - T_limiti*(wt_e1 + T_limiti*(wt_e2 - T_limiti*wt_e3)) )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... H2O activity
! Subprogram not used !     ... if the activity of H2O goes above 1.0, wt% can go negative
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                aw = ph2o_hpa / pzero_h2o
! Subprogram not used                aw = min( aw,1._r8 )
! Subprogram not used                aw = max( aw,.001_r8 )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 Molality (mol/kg)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( aw <= .05_r8 ) then
! Subprogram not used                   a1 = 12.37208932_r8
! Subprogram not used                   b1 = -0.16125516114_r8
! Subprogram not used                   c1 = -30.490657554_r8
! Subprogram not used                   d1 = -2.1133114241_r8
! Subprogram not used                   a2 = 13.455394705_r8
! Subprogram not used                   b2 = -0.1921312255_r8
! Subprogram not used                   c2 = -34.285174607_r8
! Subprogram not used                   d2 = -1.7620073078_r8
! Subprogram not used                else if( aw > .05_r8 .and. aw < .85_r8 ) then
! Subprogram not used                   a1 = 11.820654354_r8
! Subprogram not used                   b1 = -0.20786404244_r8
! Subprogram not used                   c1 = -4.807306373_r8
! Subprogram not used                   d1 = -5.1727540348_r8
! Subprogram not used                   a2 = 12.891938068_r8
! Subprogram not used                   b2 = -0.23233847708_r8
! Subprogram not used                   c2 = -6.4261237757_r8
! Subprogram not used                   d2 = -4.9005471319_r8
! Subprogram not used                else
! Subprogram not used                   a1 = -180.06541028_r8
! Subprogram not used                   b1 = -0.38601102592_r8
! Subprogram not used                   c1 = -93.317846778_r8
! Subprogram not used                   d1 = 273.88132245_r8
! Subprogram not used                   a2 = -176.95814097_r8
! Subprogram not used                   b2 = -0.36257048154_r8
! Subprogram not used                   c2 = -90.469744201_r8
! Subprogram not used                   d2 = 267.45509988_r8
! Subprogram not used                end if
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 mole fraction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                y1       = a1*(aw**b1) + c1*aw + d1
! Subprogram not used                y2       = a2*(aw**b2) + c2*aw + d2
! Subprogram not used                m_h2so4  = y1 + ((T_limit - 190._r8)*(y2 - y1)) / 70._r8
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 Weight Percent
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                wt = 9800._r8*m_h2so4 / (98._r8*m_h2so4  + 1000._r8)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     .... Parameters for h2so4 Solution, JPL-00
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 Solution Density (g/cm3)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 	       wrk = T_limit*T_limit
! Subprogram not used                z1 =  .12364_r8  - 5.6e-7_r8*wrk
! Subprogram not used                z2 = -.02954_r8  + 1.814e-7_r8*wrk
! Subprogram not used                z3 =  2.343e-3_r8 - T_limit*1.487e-6_r8 - 1.324e-8_r8*wrk
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... where mol_h2so4 is molality in mol/kg
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                den_h2so4 = 1._r8 + m_h2so4*(z1 + z2*sqrt(m_h2so4) + z3*m_h2so4)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 Molarity, mol / l
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                molar_h2so4 = den_h2so4*wt/9.8_r8
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 Mole fraction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                x_h2so4   = wt / (wt + ((100._r8 - wt)*98._r8/18._r8))
! Subprogram not used                term1     = .094_r8 - x_h2so4*(.61_r8 - 1.2_r8*x_h2so4)
! Subprogram not used                term2     = (8515._r8 - 10718._r8*(x_h2so4**.7_r8))*T_limiti
! Subprogram not used                H_hcl     = term1 * exp( -8.68_r8 + term2 )
! Subprogram not used                M_hcl     = H_hcl*pHCl_atm
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... h2so4 solution viscosity
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                aconst    = 169.5_r8 + wt*(5.18_r8 - wt*(.0825_r8 - 3.27e-3_r8*wt))
! Subprogram not used                tzero     = 144.11_r8 + wt*(.166_r8 - wt*(.015_r8 - 2.18e-4_r8*wt))
! Subprogram not used                vis_h2so4 = aconst/(T_limit**1.43_r8) * exp( 448._r8/(T_limit - tzero) )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... Acid activity in molarity
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                term1 = 60.51_r8
! Subprogram not used                term2 = .095_r8*wt
! Subprogram not used 	       wrk   = wt*wt
! Subprogram not used                term3 = .0077_r8*wrk
! Subprogram not used                term4 = 1.61e-5_r8*wt*wrk
! Subprogram not used                term5 = (1.76_r8 + 2.52e-4_r8*wrk) * T_limitsq
! Subprogram not used                term6 = -805.89_r8 + (253.05_r8*(wt**.076_r8))
! Subprogram not used                term7 = T_limitsq
! Subprogram not used                ah    = exp( term1 - term2 + term3 - term4 - term5 + term6/term7 )
! Subprogram not used 	       if( ah <= 0._r8 ) then
! Subprogram not used 	          write(iulog,*) 'ratecon: ah <= 0 at i,k, = ',i,k
! Subprogram not used 	          write(iulog,*) 'ratecon: term1,term2,term3,term4,term5,term6,term7,wt,T_limit,ah = ', &
! Subprogram not used 	                               term1,term2,term3,term4,term5,term6,term7,wt,T_limit,ah 
! Subprogram not used 	       end if
! Subprogram not used 
! Subprogram not used 	       wrk      = .25_r8*sadsulf
! Subprogram not used                rad_sulf = max( rad_sulfate(i,k),1.e-7_r8 )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     N2O5 + H2O(liq) =>  2.00*HNO3  Sulfate Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( n2o5vmr > small ) then
! Subprogram not used                   term0 = -25.5265_r8 - wt*(.133188_r8 - wt*(.00930846_r8 - 9.0194e-5_r8*wt))
! Subprogram not used                   term1 = 9283.76_r8 + wt*(115.345_r8 - wt*(5.19258_r8 - .0483464_r8*wt))
! Subprogram not used                   term2 = -851801._r8 - wt*(22191.2_r8 - wt*(766.916_r8 - 6.85427_r8*wt))
! Subprogram not used                   gprob_n2o5 = exp( term0 + T_limiti*(term1 + term2*T_limiti) )
! Subprogram not used                   rxt(i,k,rid_het1) = max( 0._r8,wrk*av_n2o5*gprob_n2o5 )
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ClONO2 + H2O(liq) =  HOCl + HNO3   Sulfate Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used ! 	... NOTE: Aerosol radius in units of cm.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... Radius sulfate set (from sad module)
! Subprogram not used !           Set min radius to 0.001 microns (1e-7 cm)
! Subprogram not used !           Typical radius is 0.1 microns (1e-5 cm)
! Subprogram not used !           f_cnt may go negative under if not set.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   C_cnt         = 1474._r8*T_limitsq
! Subprogram not used                   S_cnt         = .306_r8 + 24._r8*T_limiti
! Subprogram not used                   term1         = exp( -S_cnt*molar_h2so4 )
! Subprogram not used                   H_cnt         = 1.6e-6_r8 * exp( 4710._r8*T_limiti )*term1
! Subprogram not used                   D_cnt         = 5.e-8_r8*T_limit / vis_h2so4
! Subprogram not used                   k_h           = 1.22e12_r8*exp( -6200._r8*T_limiti )
! Subprogram not used                   k_h2o         = 1.95e10_r8*exp( -2800._r8*T_limiti )
! Subprogram not used                   k_hydr        = (k_h2o + k_h*ah)*aw
! Subprogram not used                   k_hcl         = 7.9e11_r8*ah*D_cnt*M_hcl
! Subprogram not used                   rdl_cnt       = sqrt( D_cnt/(k_hydr + k_hcl) )
! Subprogram not used                   term1         = 1._r8/tanh( rad_sulf/rdl_cnt )
! Subprogram not used                   term2         = rdl_cnt/rad_sulf
! Subprogram not used                   f_cnt         = term1 - term2
! Subprogram not used                   if( f_cnt > 0._r8 ) then
! Subprogram not used                      term1         = 4._r8*H_cnt*.082_r8*T_limit
! Subprogram not used                      term2         = sqrt( D_cnt*k_hydr )
! Subprogram not used                      Gamma_b_h2o   = term1*term2/C_cnt
! Subprogram not used                      term1         = sqrt( 1._r8 + k_hcl/k_hydr )
! Subprogram not used                      Gamma_cnt_rxn = f_cnt*Gamma_b_h2o*term1
! Subprogram not used                      Gamma_b_hcl   = Gamma_cnt_rxn*k_hcl/(k_hcl + k_hydr)
! Subprogram not used                      term1         = exp( -1374._r8*T_limiti )
! Subprogram not used                      Gamma_s       = 66.12_r8*H_cnt*M_hcl*term1
! Subprogram not used 		     if( pHCl_atm > 0._r8 ) then
! Subprogram not used                         term1      = .612_r8*(Gamma_s*Gamma_b_hcl)* pCNT_atm/pHCl_atm
! Subprogram not used                         Fhcl       = 1._r8/(1._r8 + term1)
! Subprogram not used 		     else
! Subprogram not used                         Fhcl       = 1._r8
! Subprogram not used 		     end if
! Subprogram not used                      Gamma_s_prime     = Fhcl*Gamma_s
! Subprogram not used                      Gamma_b_hcl_prime = Fhcl*Gamma_b_hcl
! Subprogram not used                      term1         = Gamma_cnt_rxn*k_hydr
! Subprogram not used                      term2         = k_hcl + k_hydr
! Subprogram not used                      Gamma_b       = Gamma_b_hcl_prime + (term1/term2)
! Subprogram not used                      term1         = 1._r8 / (Gamma_s_prime + Gamma_b)
! Subprogram not used                      gprob_cnt     = 1._r8 / (1._r8 + term1)
! Subprogram not used                      term1         = Gamma_s_prime + Gamma_b_hcl_prime
! Subprogram not used                      term2         = Gamma_s_prime + Gamma_b
! Subprogram not used                      gprob_cnt_hcl = gprob_cnt * term1/term2
! Subprogram not used                      gprob_cnt_h2o = gprob_cnt - gprob_cnt_hcl
! Subprogram not used                   else
! Subprogram not used                      gprob_cnt_h2o = 0._r8
! Subprogram not used                      gprob_cnt_hcl = 0._r8
! Subprogram not used                      Fhcl          = 1._r8
! Subprogram not used                   end if
! Subprogram not used                   if( clono2vmr > small ) then
! Subprogram not used                      rxt(i,k,rid_het2) = max( 0._r8,wrk*av_clono2*gprob_cnt_h2o )
! Subprogram not used                   end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !  	... BrONO2 + H2O(liq) =  HOBr + HNO3   Sulfate Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( brono2vmr > small ) then
! Subprogram not used                   h1    = 29.24_r8
! Subprogram not used                   h2    = -.396_r8
! Subprogram not used                   h3    = .114_r8
! Subprogram not used                   alpha = .805_r8
! Subprogram not used                   gprob_rxn = exp( h1 + h2*wt ) + h3
! Subprogram not used                   term1     = 1._r8/alpha
! Subprogram not used                   term2     = 1._r8/gprob_rxn
! Subprogram not used                   gprob_bnt_h2o = 1._r8 / (term1 + term2)
! Subprogram not used                   rxt(i,k,rid_het3) = max( 0._r8,wrk*av_brono2*gprob_bnt_h2o )
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... ClONO2 + HCl(liq) =  Cl2  + HNO3  Sulfate Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( hclvmr > small .and. clono2vmr > small ) then
! Subprogram not used                  if ( hclvmr > clono2vmr ) then
! Subprogram not used                     rxt(i,k,rid_het4) = max( 0._r8,wrk*av_clono2*gprob_cnt_hcl )*hcldeni
! Subprogram not used                  else
! Subprogram not used                     rxt(i,k,rid_het4) = max( 0._r8,wrk*av_clono2*gprob_cnt_hcl )*cntdeni
! Subprogram not used                  end if
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... HOCl + HCl(liq) =  Cl2 + H2O   Sulfate Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( hclvmr > small .and. hoclvmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... Radius sulfate set (from sad module)
! Subprogram not used !           Set min radius to 0.001 microns (1e-7 cm)
! Subprogram not used !           Typical radius is 0.1 microns (1e-5 cm)
! Subprogram not used !           f_hocl may go negative under if not set.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 	          if( pCNT_atm > 0._r8 ) then
! Subprogram not used                      D_hocl          = 6.4e-8_r8*T_limit/vis_h2so4
! Subprogram not used                      k_hocl_hcl      = 1.25e9_r8*ah*D_hocl*M_hcl
! Subprogram not used                      C_hocl          = 2009._r8*T_limitsq
! Subprogram not used                      S_hocl          = .0776_r8 + 59.18_r8*T_limiti
! Subprogram not used                      term1           = exp( -S_hocl*molar_h2so4 )
! Subprogram not used                      H_hocl          = 1.91e-6_r8 * exp( 5862.4_r8*T_limiti )*term1
! Subprogram not used                      term1           = 4._r8*H_hocl*.082_r8*T_limit
! Subprogram not used                      term2           = sqrt( D_hocl*k_hocl_hcl )
! Subprogram not used                      Gamma_hocl_rxn  = term1*term2/C_hocl
! Subprogram not used                      rdl_hocl        = sqrt( D_hocl/k_hocl_hcl )
! Subprogram not used                      term1           = 1._r8/tanh( rad_sulf/rdl_hocl )
! Subprogram not used                      term2           = rdl_hocl/rad_sulf
! Subprogram not used                      f_hocl          = term1 - term2
! Subprogram not used                      if( f_hocl > 0._r8 ) then
! Subprogram not used                         term1           = 1._r8 / (f_hocl*Gamma_hocl_rxn*Fhcl)
! Subprogram not used                         gprob_hocl_hcl  = 1._r8 / (1._r8 + term1)
! Subprogram not used                      else
! Subprogram not used                         gprob_hocl_hcl  = 0._r8
! Subprogram not used                      end if
! Subprogram not used 
! Subprogram not used                      if ( hclvmr > hoclvmr ) then
! Subprogram not used                        rxt(i,k,rid_het5) = max( 0._r8,wrk*av_hocl*gprob_hocl_hcl )*hcldeni
! Subprogram not used                      else
! Subprogram not used                        rxt(i,k,rid_het5) = max( 0._r8,wrk*av_hocl*gprob_hocl_hcl )*hocldeni
! Subprogram not used                      end if
! Subprogram not used 
! Subprogram not used 	          end if
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... HOBr + HCl(liq) =  BrCl + H2O  Sulfate Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( hclvmr > small .and. hobrvmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !   	... Radius sulfate set (from sad module)
! Subprogram not used !           Set min radius to 0.001 microns (1e-7 cm)
! Subprogram not used !           Typical radius is 0.1 microns (1e-5 cm)
! Subprogram not used !           f_hobr may go negative under if not set.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   C_hobr          = 1477._r8*T_limitsq
! Subprogram not used                   D_hobr          = 9.e-9_r8
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	...  Taken from Waschewsky and Abbat
! Subprogram not used !            Dave Hanson (PC) suggested we divide this rc by eight to agee
! Subprogram not used !            with his data (Hanson, in press, 2002).
! Subprogram not used !            k1=k2*Mhcl for gamma(HOBr)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   k_wasch         = .125_r8 * exp( .542_r8*wt - 6440._r8*T_limiti + 10.3_r8)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... Taken from Hanson 2002.
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   H_hobr          = exp( -9.86_r8 + 5427._r8*T_limiti )
! Subprogram not used                   k_dl            = 7.5e14_r8*D_hobr*2._r8                        ! or  7.5e14*D *(2nm)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !  	... If k_wasch is GE than the diffusion limit...
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 		  if( M_hcl > 0._r8 ) then
! Subprogram not used                      if( k_wasch >= k_dl ) then
! Subprogram not used                         k_hobr_hcl   = k_dl * M_hcl
! Subprogram not used                      else
! Subprogram not used                         k_hobr_hcl   = k_wasch * M_hcl
! Subprogram not used                      end if
! Subprogram not used                      term1           = 4._r8*H_hobr*.082_r8*T_limit
! Subprogram not used                      term2           = sqrt( D_hobr*k_hobr_hcl )
! Subprogram not used                      tmp             = rad_sulf/term2
! Subprogram not used                      Gamma_hobr_rxn  = term1*term2/C_hobr
! Subprogram not used                      rdl_hobr        = sqrt( D_hobr/k_hobr_hcl )
! Subprogram not used 		     if( tmp < 1.e2_r8 ) then
! Subprogram not used                         term1           = 1._r8/tanh( rad_sulf/rdl_hobr )
! Subprogram not used 		     else
! Subprogram not used                         term1           = 1._r8
! Subprogram not used 		     end if
! Subprogram not used                      term2           = rdl_hobr/rad_sulf
! Subprogram not used                      f_hobr          = term1 - term2
! Subprogram not used                      if( f_hobr > 0._r8 ) then
! Subprogram not used                         term1            = 1._r8 / (f_hobr*Gamma_hobr_rxn)
! Subprogram not used                         gprob_hobr_hcl   = 1._r8 / (1._r8 + term1)
! Subprogram not used                      else
! Subprogram not used                          gprob_hobr_hcl  = 0._r8
! Subprogram not used                      end if
! Subprogram not used 
! Subprogram not used                      if ( hclvmr > hobrvmr ) then
! Subprogram not used                         rxt(i,k,rid_het6) = max( 0._r8,wrk*av_hobr*gprob_hobr_hcl )*hcldeni
! Subprogram not used                      else
! Subprogram not used                         rxt(i,k,rid_het6) = max( 0._r8,wrk*av_hobr*gprob_hobr_hcl )*hobrdeni    
! Subprogram not used                      end if           
! Subprogram not used 
! Subprogram not used 		  end if
! Subprogram not used                end if
! Subprogram not used             end if has_sadsulf
! Subprogram not used 
! Subprogram not used has_sadnat : &
! Subprogram not used 	    if( sadnat > 0._r8 ) then
! Subprogram not used 	       wrk = .25_r8*sadnat
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... N2O5 + H2O(s) => 2HNO3  NAT Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( n2o5vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... gprob based on JPL06-2 for NAT.
! Subprogram not used !         also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
! Subprogram not used !                 gprob_tot     = 4.e-4
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   rxt(i,k,rid_het7)  = wrk*av_n2o5*4.e-4_r8
! Subprogram not used 	       end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ClONO2 + H2O(s) => HNO3 + HOCl  NAT Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( clono2vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... gprob based on JPL06-2 for NAT.
! Subprogram not used !         also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
! Subprogram not used !                 gprob_tot    = 0.004
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   rxt(i,k,rid_het8) = wrk*av_clono2*4.0e-3_r8
! Subprogram not used    	       end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... ClONO2 + HCl(s) => HNO3 + Cl2, NAT Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( hclvmr > small ) then
! Subprogram not used                   if( clono2vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... gprob based on JPL06-2 for NAT.
! Subprogram not used !         also see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
! Subprogram not used !                 gprob_tot   = 0.2
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                      if ( hclvmr > clono2vmr ) then
! Subprogram not used                         rxt(i,k,rid_het9) = wrk*av_clono2*0.2_r8*hcldeni
! Subprogram not used                      else
! Subprogram not used                         rxt(i,k,rid_het9) = wrk*av_clono2*0.2_r8*cntdeni  
! Subprogram not used                      end if
! Subprogram not used                   end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... HOCl + HCl(s) => H2O + Cl2  NAT Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   if( hoclvmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... gprob based on JPL06-2 for NAT.
! Subprogram not used !         see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
! Subprogram not used !         and      Abbatt and Molina, GRL, 19, 461-464, 1992.
! Subprogram not used !                 gprob_tot   = 0.1
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                      if ( hclvmr > hoclvmr ) then
! Subprogram not used                         rxt(i,k,rid_het10) = wrk*av_hocl*0.1_r8*hcldeni
! Subprogram not used                      else
! Subprogram not used                         rxt(i,k,rid_het10) = wrk*av_hocl*0.1_r8*hocldeni
! Subprogram not used                      end if
! Subprogram not used                   end if
! Subprogram not used                end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... BrONO2 + H2O(s) => HOBr + HNO3  NAT Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( brono2vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !       ... Personel Communication, 11/4/99, David Hanson
! Subprogram not used !                 gprob_tot   = 0.3
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   rxt(i,k,rid_het11) = wrk*av_brono2*0.3_r8
! Subprogram not used                end if
! Subprogram not used             end if has_sadnat
! Subprogram not used 
! Subprogram not used has_sadice : &
! Subprogram not used 	    if( sadice > 0._r8 ) then
! Subprogram not used 	       wrk = .25_r8*sadice
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     N2O5 + H2O(s) => 2HNO3  ICE Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( n2o5vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !       ... gprob based on JPL06-2 for ICE.
! Subprogram not used !           also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
! Subprogram not used !                 gprob_tot    = .02
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   rxt(i,k,rid_het12) = wrk*av_n2o5*0.02_r8
! Subprogram not used  	       end if
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... ClONO2 + H2O(s) => HNO3 + HOCl  ICE Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( clono2vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... gprob based on JPL06-2 for ICE.
! Subprogram not used !     	    also see Hanson and Ravi, JGR, 96, 17307-17314, 1991.
! Subprogram not used !                 gprob_tot    = .3
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   rxt(i,k,rid_het13) = wrk*av_clono2*0.3_r8
! Subprogram not used 	       end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... BrONO2 + H2O(s) => HNO3 + HOBr  ICE Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( brono2vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... gprob based on JPL06-2 for ICE.
! Subprogram not used !           also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
! Subprogram not used !           could be as high as 1.0
! Subprogram not used !                 gprob_tot    = .3
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   rxt(i,k,rid_het14) = wrk*av_brono2*0.3_r8
! Subprogram not used       	       end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ClONO2 + HCl(s) => HNO3 + Cl2, ICE Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                if( hclvmr > small ) then
! Subprogram not used                   if( clono2vmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !       ... gprob based on JPL06-2 for ICE.
! Subprogram not used !           also see Hanson and Ravi, GRL, 15, 17-20, 1988.
! Subprogram not used !           also see Lue et al.,
! Subprogram not used !                 gprob_tot    = .3
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                      if ( hclvmr > clono2vmr ) then
! Subprogram not used                         rxt(i,k,rid_het15) = wrk*av_clono2*0.3_r8*hcldeni
! Subprogram not used                      else
! Subprogram not used                         rxt(i,k,rid_het15) = wrk*av_clono2*0.3_r8*cntdeni
! Subprogram not used                      end if
! Subprogram not used 
! Subprogram not used                   end if
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... HOCl + HCl(s) => H2O + Cl2, ICE Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   if( hoclvmr > small .and. hclvmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !       ... gprob based on JPL06-2 for ICE.
! Subprogram not used !           also see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
! Subprogram not used !           also see Abbatt and Molina, GRL, 19, 461-464, 1992.
! Subprogram not used !                 gprob_tot   = .2
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                      if ( hclvmr > hoclvmr ) then
! Subprogram not used                         rxt(i,k,rid_het16) = wrk*av_hocl*0.2_r8*hcldeni
! Subprogram not used                      else
! Subprogram not used                         rxt(i,k,rid_het16) = wrk*av_hocl*0.2_r8*hocldeni
! Subprogram not used                      end if
! Subprogram not used 
! Subprogram not used                   end if
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     HOBr + HCl(s) => H2O + BrCl, ICE Aerosol Reaction
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                   if( hobrvmr > small .and. hclvmr > small ) then
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !       ... gprob based on JPL06-2 for ICE.
! Subprogram not used !           Abbatt GRL, 21, 665-668, 1994.
! Subprogram not used !                    gprob_tot   = .3
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used                     if ( hclvmr > hobrvmr ) then
! Subprogram not used                        rxt(i,k,rid_het17) = wrk*av_hobr*0.3_r8*hcldeni
! Subprogram not used                     else
! Subprogram not used                        rxt(i,k,rid_het17) = wrk*av_hobr*0.3_r8*hobrdeni
! Subprogram not used                     end if
! Subprogram not used                   end if
! Subprogram not used 	       end if
! Subprogram not used             end if has_sadice
! Subprogram not used          end do column_loop
! Subprogram not used       end do Level_loop
! Subprogram not used 
! Subprogram not used       end subroutine ratecon_sfstrat

      end module mo_strato_rates
