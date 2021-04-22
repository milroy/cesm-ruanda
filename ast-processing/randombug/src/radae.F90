module radae
!------------------------------------------------------------------------------
!
! Description:
!
! Data and subroutines to calculate absorptivities and emissivity needed
! for the LW radiation calculation.
!
! Public interfaces are: 
!
! radae_init ------------ Initialization
! initialize_radbuffer -- Initialize the 3D abs/emis arrays.
! radabs ---------------- Compute absorptivities.
! radems ---------------- Compute emissivity.
! radtpl ---------------- Compute Temperatures and path lengths.
! radoz2 ---------------- Compute ozone path lengths.
! trcpth ---------------- Compute ghg path lengths.
!
! Author:  B. Collins
!
! $Id$
!------------------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use ppgrid,           only: pcols, pverp, begchunk, endchunk, pver
  use infnan,           only: posinf, assignment(=)
  use pmgrid,           only: plev, plevp
  use radconstants,     only: nlwbands, idx_LW_0650_0800, idx_LW_0500_0650, &
        idx_LW_1000_1200, idx_LW_0800_1000,  idx_LW_1200_2000
  use abortutils,       only: endrun
  use cam_logfile,      only: iulog
  use wv_saturation,    only: qsat_water


  implicit none

  save
!-----------------------------------------------------------------------------
! PUBLIC:: By default data and interfaces are private
!-----------------------------------------------------------------------------
  private
  public radabs, radems, radtpl, radae_init, initialize_radbuffer, radoz2, trcpth

  integer, public, parameter :: nbands = 2          ! Number of spectral bands
!
! Following data needed for restarts and in radclwmx
!
  real(r8), public, allocatable, target :: abstot_3d(:,:,:,:) ! Non-adjacent layer absorptivites
  real(r8), public, allocatable, target :: absnxt_3d(:,:,:,:) ! Nearest layer absorptivities
  real(r8), public, allocatable, target :: emstot_3d(:,:,:)   ! Total emissivity
  integer,  public :: ntoplw    ! top level to solve for longwave cooling

!-----------------------------------------------------------------------------
! PRIVATE:: The rest of the data is private to this module.
!-----------------------------------------------------------------------------
  real(r8) :: p0    ! Standard pressure (dynes/cm**2)
  real(r8) :: amd   ! Molecular weight of dry air (g/mol)
  real(r8) :: amco2 ! Molecular weight of co2   (g/mol)
  real(r8) :: mwo3  ! Molecular weight of O3 (g/mol)

  real(r8) :: gravit     ! acceleration due to gravity (m/s**2)
  real(r8) :: gravit_cgs ! acceleration due to gravity (cm/s**2)
  real(r8) :: rga        ! 1./gravit_cgs
  real(r8) :: epsilo     ! Ratio of mol. wght of H2O to dry air
  real(r8) :: omeps      ! 1._r8 - epsilo
  real(r8) :: sslp       ! Standard sea-level pressure (dynes/cm**2)
  real(r8) :: stebol_cgs ! Stefan-Boltzmann's constant (CGS)
  real(r8) :: rgsslp     ! 0.5/(gravit_cgs*sslp)
  real(r8) :: dpfo3      ! Voigt correction factor for O3
  real(r8) :: dpfco2     ! Voigt correction factor for CO2

  integer, parameter :: n_u = 25   ! Number of U in abs/emis tables
  integer, parameter :: n_p = 10   ! Number of P in abs/emis tables
  integer, parameter :: n_tp = 10  ! Number of T_p in abs/emis tables
  integer, parameter :: n_te = 21  ! Number of T_e in abs/emis tables
  integer, parameter :: n_rh = 7   ! Number of RH in abs/emis tables

  real(r8):: ah2onw(n_p, n_tp, n_u, n_te, n_rh)   ! absorptivity (non-window)
  real(r8):: eh2onw(n_p, n_tp, n_u, n_te, n_rh)   ! emissivity   (non-window)
  real(r8):: ah2ow(n_p, n_tp, n_u, n_te, n_rh)    ! absorptivity (window, for adjacent layers)
  real(r8):: cn_ah2ow(n_p, n_tp, n_u, n_te, n_rh)    ! continuum transmission for absorptivity (window)
  real(r8):: cn_eh2ow(n_p, n_tp, n_u, n_te, n_rh)    ! continuum transmission for emissivity   (window)
  real(r8):: ln_ah2ow(n_p, n_tp, n_u, n_te, n_rh)    ! line-only transmission for absorptivity (window)
  real(r8):: ln_eh2ow(n_p, n_tp, n_u, n_te, n_rh)    ! line-only transmission for emissivity   (window)
!
! Constant coefficients for water vapor overlap with trace gases.
! Reference: Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
  real(r8):: coefh(2,4) = reshape(  &
         (/ (/5.46557e+01_r8,-7.30387e-02_r8/), &
            (/1.09311e+02_r8,-1.46077e-01_r8/), &
            (/5.11479e+01_r8,-6.82615e-02_r8/), &
            (/1.02296e+02_r8,-1.36523e-01_r8/) /), (/2,4/) )
!
  real(r8):: coefj(3,2) = reshape( &
            (/ (/2.82096e-02_r8,2.47836e-04_r8,1.16904e-06_r8/), &
               (/9.27379e-02_r8,8.04454e-04_r8,6.88844e-06_r8/) /), (/3,2/) )
!
  real(r8):: coefk(3,2) = reshape( &
            (/ (/2.48852e-01_r8,2.09667e-03_r8,2.60377e-06_r8/) , &
               (/1.03594e+00_r8,6.58620e-03_r8,4.04456e-06_r8/) /), (/3,2/) )
  real(r8):: c16,c17,c26,c27,c28,c29,c30,c31
!
! Farwing correction constants for narrow-band emissivity model,
! introduced to account for the deficiencies in narrow-band model
! used to derive the emissivity; tuned with Arkings line-by-line
! calculations.  Just used for water vapor overlap with trace gases.
!
  real(r8):: fwcoef      ! Farwing correction constant
  real(r8):: fwc1,fwc2   ! Farwing correction constants 
  real(r8):: fc1         ! Farwing correction constant 
!
! Collins/Hackney/Edwards (C/H/E) & Collins/Lee-Taylor/Edwards (C/LT/E) 
!       H2O parameterization
!
! Notation:
! U   = integral (P/P_0 dW)  eq. 15 in Ramanathan/Downey 1986
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
! absorptivity/emissivity in window are fit using an expression:
!
!      a/e = f_a/e * {1.0 - ln_a/e * cn_a/e} 
!
! absorptivity/emissivity in non-window are fit using:
! 
!      a/e = f_a/e * a/e_norm
!
! where
!      a/e = absorptivity/emissivity
! a/e_norm = absorptivity/emissivity normalized to 1
!    f_a/e = value of a/e as U->infinity = f(T_e) only
!   cn_a/e = continuum transmission
!   ln_a/e = line transmission
!
! spectral interval:
!   1 = 0-800 cm^-1 and 1200-2200 cm^-1 (rotation and rotation-vibration)
!   2 = 800-1200 cm^-1                  (window)
!
  real(r8), parameter:: min_tp_h2o = 160.0_r8        ! min T_p for pre-calculated abs/emis 
  real(r8), parameter:: max_tp_h2o = 349.999999_r8   ! max T_p for pre-calculated abs/emis
  integer, parameter :: ntemp = 192 ! Number of temperatures in H2O sat. table for Tp
  integer, parameter :: o_fa = 6   ! Degree+1 of poly of T_e for absorptivity as U->inf.
  integer, parameter :: o_fe = 6   ! Degree+1 of poly of T_e for emissivity as U->inf.
!-----------------------------------------------------------------------------
! Data for f in C/H/E fit -- value of A and E as U->infinity
! New C/LT/E fit (Hitran 2K, CKD 2.4) -- no change
!     These values are determined by integrals of Planck functions or
!     derivatives of Planck functions only.
!-----------------------------------------------------------------------------
!
! fa/fe coefficients for 2 bands (0-800 & 1200-2200, 800-1200 cm^-1)
!
! Coefficients of polynomial for f_a in T_e
!
  real(r8), parameter:: fat(o_fa,nbands) = reshape( (/ &
       (/-1.06665373E-01_r8,  2.90617375E-02_r8, -2.70642049E-04_r8,   &   ! 0-800&1200-2200 cm^-1
          1.07595511E-06_r8, -1.97419681E-09_r8,  1.37763374E-12_r8/), &   !   0-800&1200-2200 cm^-1
       (/ 1.10666537E+00_r8, -2.90617375E-02_r8,  2.70642049E-04_r8,   &   ! 800-1200 cm^-1
         -1.07595511E-06_r8,  1.97419681E-09_r8, -1.37763374E-12_r8/) /) & !   800-1200 cm^-1
       , (/o_fa,nbands/) )
!
! Coefficients of polynomial for f_e in T_e
!
  real(r8), parameter:: fet(o_fe,nbands) = reshape( (/ & 
      (/3.46148163E-01_r8,  1.51240299E-02_r8, -1.21846479E-04_r8,   &   ! 0-800&1200-2200 cm^-1
        4.04970123E-07_r8, -6.15368936E-10_r8,  3.52415071E-13_r8/), &   !   0-800&1200-2200 cm^-1
      (/6.53851837E-01_r8, -1.51240299E-02_r8,  1.21846479E-04_r8,   &   ! 800-1200 cm^-1
       -4.04970123E-07_r8,  6.15368936E-10_r8, -3.52415071E-13_r8/) /) & !   800-1200 cm^-1
      , (/o_fa,nbands/) )
!
! Note: max values should be slightly underestimated to avoid index bound violations
!
  real(r8), parameter:: min_lp_h2o = -3.0_r8         ! min log_10(P) for pre-calculated abs/emis 
  real(r8), parameter:: min_p_h2o = 1.0e-3_r8        ! min log_10(P) for pre-calculated abs/emis 
  real(r8), parameter:: max_lp_h2o = -0.0000001_r8   ! max log_10(P) for pre-calculated abs/emis 
  real(r8), parameter:: dlp_h2o = 0.3333333333333_r8 ! difference in adjacent elements of lp_h2o
 
  real(r8), parameter:: dtp_h2o = 21.111111111111_r8 ! difference in adjacent elements of tp_h2o

  real(r8), parameter:: min_rh_h2o = 0.0_r8          ! min RH for pre-calculated abs/emis 
  real(r8), parameter:: max_rh_h2o = 1.19999999_r8   ! max RH for pre-calculated abs/emis 
  real(r8), parameter:: drh_h2o = 0.2_r8             ! difference in adjacent elements of RH

  real(r8), parameter:: min_te_h2o = -120.0_r8       ! min T_e-T_p for pre-calculated abs/emis 
  real(r8), parameter:: max_te_h2o = 79.999999_r8    ! max T_e-T_p for pre-calculated abs/emis 
  real(r8), parameter:: dte_h2o  = 10.0_r8           ! difference in adjacent elements of te_h2o

  real(r8), parameter:: min_lu_h2o = -8.0_r8         ! min log_10(U) for pre-calculated abs/emis 
  real(r8), parameter:: min_u_h2o  = 1.0e-8_r8       ! min pressure-weighted path-length
  real(r8), parameter:: max_lu_h2o =  3.9999999_r8   ! max log_10(U) for pre-calculated abs/emis 
  real(r8), parameter:: dlu_h2o  = 0.5_r8            ! difference in adjacent elements of lu_h2o

  real(r8), parameter:: g1(6)=(/0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/)
  real(r8), parameter :: g2(6)=(/14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/)
  real(r8), parameter :: g3(6)=(/26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/)
  real(r8), parameter :: g4(6)=(/0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/)
  real(r8), parameter :: ab(6)=(/3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/)
  real(r8), parameter :: bb(6)=(/-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/)
  real(r8), parameter :: abp(6)=(/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/)
  real(r8), parameter :: bbp(6)=(/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/)


! Public Interfaces
!====================================================================================
CONTAINS
!====================================================================================

! Subprogram not used subroutine radabs(lchnk   ,ncol    ,             &
! Subprogram not used    pbr    ,pnm     ,co2em    ,co2eml  ,tplnka  , &
! Subprogram not used    s2c    ,tcg     ,w        ,h2otr   ,plco2   , &
! Subprogram not used    plh2o  ,co2t    ,tint     ,tlayr   ,plol    , &
! Subprogram not used    plos   ,pmln    ,piln     ,ucfc11  ,ucfc12  , &
! Subprogram not used    un2o0  ,un2o1   ,uch4     ,uco211  ,uco212  , &
! Subprogram not used    uco213 ,uco221  ,uco222   ,uco223  ,uptype  , &
! Subprogram not used    bn2o0  ,bn2o1   ,bch4    ,abplnk1  ,abplnk2 , &
! Subprogram not used    abstot ,absnxt  ,plh2ob  ,wb       , &
! Subprogram not used    odap_aer ,aer_trn_ttl, co2mmr)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! h2o  ....  Uses nonisothermal emissivity method for water vapor from
! Subprogram not used !            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
! Subprogram not used !            Emissivity and Absorptivity Formulation for Water Vapor
! Subprogram not used !            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
! Subprogram not used !
! Subprogram not used !            Implementation updated by Collins, Hackney, and Edwards (2001)
! Subprogram not used !               using line-by-line calculations based upon Hitran 1996 and
! Subprogram not used !               CKD 2.1 for absorptivity and emissivity
! Subprogram not used !
! Subprogram not used !            Implementation updated by Collins, Lee-Taylor, and Edwards (2003)
! Subprogram not used !               using line-by-line calculations based upon Hitran 2000 and
! Subprogram not used !               CKD 2.4 for absorptivity and emissivity
! Subprogram not used !
! Subprogram not used ! co2  ....  Uses absorptance parameterization of the 15 micro-meter
! Subprogram not used !            (500 - 800 cm-1) band system of Carbon Dioxide, from
! Subprogram not used !            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
! Subprogram not used !            of the Absorptance Due to the 15 micro-meter Band System
! Subprogram not used !            of Carbon Dioxide Jouranl of Geophysical Research,
! Subprogram not used !            vol. 96., D5, pp 9013-9019.
! Subprogram not used !            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
! Subprogram not used !            are also included.
! Subprogram not used !
! Subprogram not used ! o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
! Subprogram not used !            band system of ozone, from Ramanathan, V. and R.Dickinson,
! Subprogram not used !            1979: The Role of stratospheric ozone in the zonal and
! Subprogram not used !            seasonal radiative energy balance of the earth-troposphere
! Subprogram not used !            system. Journal of the Atmospheric Sciences, Vol. 36,
! Subprogram not used !            pp 1084-1104
! Subprogram not used !
! Subprogram not used ! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
! Subprogram not used !
! Subprogram not used ! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
! Subprogram not used !            bands of nitrous oxide
! Subprogram not used !
! Subprogram not used ! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
! Subprogram not used !            micron bands of CFC11
! Subprogram not used !
! Subprogram not used ! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
! Subprogram not used !            micron bands of CFC12
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Computes individual absorptivities for non-adjacent layers, accounting
! Subprogram not used ! for band overlap, and sums to obtain the total; then, computes the
! Subprogram not used ! nearest layer contribution.
! Subprogram not used ! 
! Subprogram not used ! Author: W. Collins (H2O absorptivity) and J. Kiehl
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: lchnk                       ! chunk identifier
! Subprogram not used    integer, intent(in) :: ncol                        ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: pbr(pcols,pver)            ! Prssr at mid-levels (dynes/cm2)
! Subprogram not used    real(r8), intent(in) :: pnm(pcols,pverp)           ! Prssr at interfaces (dynes/cm2)
! Subprogram not used    real(r8), intent(in) :: co2em(pcols,pverp)         ! Co2 emissivity function
! Subprogram not used    real(r8), intent(in) :: co2eml(pcols,pver)         ! Co2 emissivity function
! Subprogram not used    real(r8), intent(in) :: tplnka(pcols,pverp)        ! Planck fnctn level temperature
! Subprogram not used    real(r8), intent(in) :: s2c(pcols,pverp)           ! H2o continuum path length
! Subprogram not used    real(r8), intent(in) :: tcg(pcols,pverp)           ! H2o-mass-wgted temp. (Curtis-Godson approx.)
! Subprogram not used    real(r8), intent(in) :: w(pcols,pverp)             ! H2o prs wghted path
! Subprogram not used    real(r8), intent(in) :: h2otr(pcols,pverp)         ! H2o trnsmssn fnct for o3 overlap
! Subprogram not used    real(r8), intent(in) :: plco2(pcols,pverp)         ! Co2 prs wghted path length
! Subprogram not used    real(r8), intent(in) :: plh2o(pcols,pverp)         ! H2o prs wfhted path length
! Subprogram not used    real(r8), intent(in) :: co2t(pcols,pverp)          ! Tmp and prs wghted path length
! Subprogram not used    real(r8), intent(in) :: tint(pcols,pverp)          ! Interface temperatures
! Subprogram not used    real(r8), intent(in) :: tlayr(pcols,pverp)         ! K-1 level temperatures
! Subprogram not used    real(r8), intent(in) :: plol(pcols,pverp)          ! Ozone prs wghted path length
! Subprogram not used    real(r8), intent(in) :: plos(pcols,pverp)          ! Ozone path length
! Subprogram not used    real(r8), intent(in) :: pmln(pcols,pver)           ! Ln(pmidm1)
! Subprogram not used    real(r8), intent(in) :: piln(pcols,pverp)          ! Ln(pintm1)
! Subprogram not used    real(r8), intent(in) :: plh2ob(nbands,pcols,pverp) ! Pressure weighted h2o path with 
! Subprogram not used                                                       !    Hulst-Curtis-Godson temp. factor 
! Subprogram not used                                                       !    for H2O bands 
! Subprogram not used    real(r8), intent(in) :: wb(nbands,pcols,pverp)     ! H2o path length with 
! Subprogram not used                                                       !    Hulst-Curtis-Godson temp. factor 
! Subprogram not used                                                       !    for H2O bands 
! Subprogram not used 
! Subprogram not used ! [fraction] absorbtion optical depth, cumulative from top
! Subprogram not used    real(r8), intent(in) :: odap_aer(pcols,pver,nlwbands)
! Subprogram not used 
! Subprogram not used ! [fraction] Total transmission between interfaces k1 and k2
! Subprogram not used    real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Trace gas variables
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: co2mmr(pcols)              ! co2 column mean mass mixing ratio
! Subprogram not used    real(r8), intent(in) :: ucfc11(pcols,pverp)        ! CFC11 path length
! Subprogram not used    real(r8), intent(in) :: ucfc12(pcols,pverp)        ! CFC12 path length
! Subprogram not used    real(r8), intent(in) :: un2o0(pcols,pverp)         ! N2O path length
! Subprogram not used    real(r8), intent(in) :: un2o1(pcols,pverp)         ! N2O path length (hot band)
! Subprogram not used    real(r8), intent(in) :: uch4(pcols,pverp)          ! CH4 path length
! Subprogram not used    real(r8), intent(in) :: uco211(pcols,pverp)        ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco212(pcols,pverp)        ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco213(pcols,pverp)        ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco221(pcols,pverp)        ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco222(pcols,pverp)        ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco223(pcols,pverp)        ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uptype(pcols,pverp)        ! continuum path length
! Subprogram not used    real(r8), intent(in) :: bn2o0(pcols,pverp)         ! pressure factor for n2o
! Subprogram not used    real(r8), intent(in) :: bn2o1(pcols,pverp)         ! pressure factor for n2o
! Subprogram not used    real(r8), intent(in) :: bch4(pcols,pverp)          ! pressure factor for ch4
! Subprogram not used    real(r8), intent(in) :: abplnk1(14,pcols,pverp)    ! non-nearest layer Planck factor
! Subprogram not used    real(r8), intent(in) :: abplnk2(14,pcols,pverp)    ! nearest layer factor
! Subprogram not used !
! Subprogram not used ! Output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: abstot(pcols,ntoplw:pverp,ntoplw:pverp) ! Total absorptivity
! Subprogram not used    real(r8), intent(out) :: absnxt(pcols,pver,4)      ! Total nearest layer absorptivity
! Subprogram not used !
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used !
! Subprogram not used    integer i                   ! Longitude index
! Subprogram not used    integer k                   ! Level index
! Subprogram not used    integer k1                  ! Level index
! Subprogram not used    integer k2                  ! Level index
! Subprogram not used    integer kn                  ! Nearest level index
! Subprogram not used    integer wvl                 ! Wavelength index
! Subprogram not used 
! Subprogram not used    real(r8) abstrc(pcols)              ! total trace gas absorptivity
! Subprogram not used    real(r8) bplnk(14,pcols,4)          ! Planck functions for sub-divided layers
! Subprogram not used    real(r8) pnew(pcols)        ! Effective pressure for H2O vapor linewidth
! Subprogram not used    real(r8) pnewb(nbands)      ! Effective pressure for h2o linewidth w/
! Subprogram not used                                !    Hulst-Curtis-Godson correction for
! Subprogram not used                                !    each band
! Subprogram not used    real(r8) u(pcols)           ! Pressure weighted H2O path length
! Subprogram not used    real(r8) ub(nbands)         ! Pressure weighted H2O path length with
! Subprogram not used                                !    Hulst-Curtis-Godson correction for
! Subprogram not used                                !    each band
! Subprogram not used    real(r8) tbar(pcols,4)      ! Mean layer temperature
! Subprogram not used    real(r8) emm(pcols,4)       ! Mean co2 emissivity
! Subprogram not used    real(r8) o3emm(pcols,4)     ! Mean o3 emissivity
! Subprogram not used    real(r8) o3bndi             ! Ozone band parameter
! Subprogram not used    real(r8) temh2o(pcols,4)    ! Mean layer temperature equivalent to tbar
! Subprogram not used    real(r8) k21                ! Exponential coefficient used to calculate
! Subprogram not used !                              !  rotation band transmissvty in the 650-800
! Subprogram not used !                              !  cm-1 region (tr1)
! Subprogram not used    real(r8) k22                ! Exponential coefficient used to calculate
! Subprogram not used !                              !  rotation band transmissvty in the 500-650
! Subprogram not used !                              !  cm-1 region (tr2)
! Subprogram not used    real(r8) uc1(pcols)         ! H2o continuum pathlength in 500-800 cm-1
! Subprogram not used    real(r8) to3h2o(pcols)      ! H2o trnsmsn for overlap with o3
! Subprogram not used    real(r8) pi                 ! For co2 absorptivity computation
! Subprogram not used    real(r8) sqti(pcols)        ! Used to store sqrt of mean temperature
! Subprogram not used    real(r8) et                 ! Co2 hot band factor
! Subprogram not used    real(r8) et2                ! Co2 hot band factor squared
! Subprogram not used    real(r8) et4                ! Co2 hot band factor to fourth power
! Subprogram not used    real(r8) omet               ! Co2 stimulated emission term
! Subprogram not used    real(r8) f1co2              ! Co2 central band factor
! Subprogram not used    real(r8) f2co2(pcols)       ! Co2 weak band factor
! Subprogram not used    real(r8) f3co2(pcols)       ! Co2 weak band factor
! Subprogram not used    real(r8) t1co2(pcols)       ! Overlap factr weak bands on strong band
! Subprogram not used    real(r8) sqwp               ! Sqrt of co2 pathlength
! Subprogram not used    real(r8) f1sqwp(pcols)      ! Main co2 band factor
! Subprogram not used    real(r8) oneme              ! Co2 stimulated emission term
! Subprogram not used    real(r8) alphat             ! Part of the co2 stimulated emission term
! Subprogram not used    real(r8) co2vmr(pcols)      ! CO2 column mean vmr
! Subprogram not used    real(r8) rmw                ! ratio of molecular weights (air/co2)
! Subprogram not used    real(r8) wco2               ! Constants used to define co2 pathlength
! Subprogram not used    real(r8) posqt              ! Effective pressure for co2 line width
! Subprogram not used    real(r8) u7(pcols)          ! Co2 hot band path length
! Subprogram not used    real(r8) u8                 ! Co2 hot band path length
! Subprogram not used    real(r8) u9                 ! Co2 hot band path length
! Subprogram not used    real(r8) u13                ! Co2 hot band path length
! Subprogram not used    real(r8) rbeta7(pcols)      ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) rbeta8             ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) rbeta9             ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) rbeta13            ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) tpatha             ! For absorptivity computation
! Subprogram not used    real(r8) abso(pcols,4)      ! Absorptivity for various gases/bands
! Subprogram not used    real(r8) dtx(pcols)         ! Planck temperature minus 250 K
! Subprogram not used    real(r8) dty(pcols)         ! Path temperature minus 250 K
! Subprogram not used    real(r8) term7(pcols,2)     ! Kl_inf(i) in eq(r8) of table A3a of R&D
! Subprogram not used    real(r8) term8(pcols,2)     ! Delta kl_inf(i) in eq(r8)
! Subprogram not used    real(r8) tr1                ! Eqn(6) in table A2 of R&D for 650-800
! Subprogram not used    real(r8) tr10(pcols)        ! Eqn (6) times eq(4) in table A2
! Subprogram not used !                              !  of R&D for 500-650 cm-1 region
! Subprogram not used    real(r8) tr2                ! Eqn(6) in table A2 of R&D for 500-650
! Subprogram not used    real(r8) tr5                ! Eqn(4) in table A2 of R&D for 650-800
! Subprogram not used    real(r8) tr6                ! Eqn(4) in table A2 of R&D for 500-650
! Subprogram not used    real(r8) tr9(pcols)         ! Equation (6) times eq(4) in table A2
! Subprogram not used !                              !  of R&D for 650-800 cm-1 region
! Subprogram not used    real(r8) sqrtu(pcols)       ! Sqrt of pressure weighted h20 pathlength
! Subprogram not used    real(r8) fwk(pcols)         ! Equation(33) in R&D far wing correction
! Subprogram not used    real(r8) fwku(pcols)        ! GU term in eqs(1) and (6) in table A2
! Subprogram not used    real(r8) to3co2(pcols)      ! P weighted temp in ozone band model
! Subprogram not used    real(r8) dpnm(pcols)        ! Pressure difference between two levels
! Subprogram not used    real(r8) pnmsq(pcols,pverp) ! Pressure squared
! Subprogram not used    real(r8) dw(pcols)          ! Amount of h2o between two levels
! Subprogram not used    real(r8) uinpl(pcols,4)     ! Nearest layer subdivision factor
! Subprogram not used    real(r8) winpl(pcols,4)     ! Nearest layer subdivision factor
! Subprogram not used    real(r8) zinpl(pcols,4)     ! Nearest layer subdivision factor
! Subprogram not used    real(r8) pinpl(pcols,4)     ! Nearest layer subdivision factor
! Subprogram not used    real(r8) dplh2o(pcols)      ! Difference in press weighted h2o amount
! Subprogram not used    real(r8) r293               ! 1/293
! Subprogram not used    real(r8) r250               ! 1/250
! Subprogram not used    real(r8) r3205              ! Line width factor for o3 (see R&Di)
! Subprogram not used    real(r8) r300               ! 1/300
! Subprogram not used    real(r8) rsslp              ! Reciprocal of sea level pressure
! Subprogram not used    real(r8) r2sslp             ! 1/2 of rsslp
! Subprogram not used    real(r8) ds2c               ! Y in eq(7) in table A2 of R&D
! Subprogram not used    real(r8)  dplos             ! Ozone pathlength eq(A2) in R&Di
! Subprogram not used    real(r8) dplol              ! Presure weighted ozone pathlength
! Subprogram not used    real(r8) tlocal             ! Local interface temperature
! Subprogram not used    real(r8) beta               ! Ozone mean line parameter eq(A3) in R&Di
! Subprogram not used !                               (includes Voigt line correction factor)
! Subprogram not used    real(r8) rphat              ! Effective pressure for ozone beta
! Subprogram not used    real(r8) tcrfac             ! Ozone temperature factor table 1 R&Di
! Subprogram not used    real(r8) tmp1               ! Ozone band factor see eq(A1) in R&Di
! Subprogram not used    real(r8) u1                 ! Effective ozone pathlength eq(A2) in R&Di
! Subprogram not used    real(r8) realnu             ! 1/beta factor in ozone band model eq(A1)
! Subprogram not used    real(r8) tmp2               ! Ozone band factor see eq(A1) in R&Di
! Subprogram not used    real(r8) u2                 ! Effective ozone pathlength eq(A2) in R&Di
! Subprogram not used    real(r8) rsqti              ! Reciprocal of sqrt of path temperature
! Subprogram not used    real(r8) tpath              ! Path temperature used in co2 band model
! Subprogram not used    real(r8) tmp3               ! Weak band factor see K&B
! Subprogram not used    real(r8) rdpnmsq            ! Reciprocal of difference in press^2
! Subprogram not used    real(r8) rdpnm              ! Reciprocal of difference in press
! Subprogram not used    real(r8) p1                 ! Mean pressure factor
! Subprogram not used    real(r8) p2                 ! Mean pressure factor
! Subprogram not used    real(r8) dtym10             ! T - 260 used in eq(9) and (10) table A3a
! Subprogram not used    real(r8) dplco2             ! Co2 path length
! Subprogram not used    real(r8) te                 ! A_0 T factor in ozone model table 1 of R&Di
! Subprogram not used    real(r8) denom              ! Denominator in eq(r8) of table A3a of R&D
! Subprogram not used    real(r8) th2o(pcols)        ! transmission due to H2O
! Subprogram not used    real(r8) tco2(pcols)        ! transmission due to CO2
! Subprogram not used    real(r8) to3(pcols)         ! transmission due to O3
! Subprogram not used !
! Subprogram not used ! Transmission terms for various spectral intervals:
! Subprogram not used !
! Subprogram not used    real(r8) trab2(pcols)       ! H2o   500 -  800 cm-1
! Subprogram not used    real(r8) absbnd             ! Proportional to co2 band absorptance
! Subprogram not used    real(r8) dbvtit(pcols,pverp)! Intrfc drvtv plnck fnctn for o3
! Subprogram not used    real(r8) dbvtly(pcols,pver) ! Level drvtv plnck fnctn for o3
! Subprogram not used !
! Subprogram not used ! Variables for Collins/Hackney/Edwards (C/H/E) & 
! Subprogram not used !       Collins/Lee-Taylor/Edwards (C/LT/E) H2O parameterization
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Notation:
! Subprogram not used ! U   = integral (P/P_0 dW)  eq. 15 in Ramanathan/Downey 1986
! Subprogram not used ! P   = atmospheric pressure
! Subprogram not used ! P_0 = reference atmospheric pressure
! Subprogram not used ! W   = precipitable water path
! Subprogram not used ! T_e = emission temperature
! Subprogram not used ! T_p = path temperature
! Subprogram not used ! RH  = path relative humidity
! Subprogram not used !
! Subprogram not used    real(r8) fa               ! asymptotic value of abs. as U->infinity
! Subprogram not used    real(r8) a_star           ! normalized absorptivity for non-window
! Subprogram not used    real(r8) l_star           ! interpolated line transmission
! Subprogram not used    real(r8) c_star           ! interpolated continuum transmission
! Subprogram not used 
! Subprogram not used    real(r8) te1              ! emission temperature
! Subprogram not used    real(r8) te2              ! te^2
! Subprogram not used    real(r8) te3              ! te^3
! Subprogram not used    real(r8) te4              ! te^4
! Subprogram not used    real(r8) te5              ! te^5
! Subprogram not used 
! Subprogram not used    real(r8) log_u            ! log base 10 of U 
! Subprogram not used    real(r8) log_uc           ! log base 10 of H2O continuum path
! Subprogram not used    real(r8) log_p            ! log base 10 of P
! Subprogram not used    real(r8) t_p              ! T_p
! Subprogram not used    real(r8) t_e              ! T_e (offset by T_p)
! Subprogram not used 
! Subprogram not used    integer iu                ! index for log10(U)
! Subprogram not used    integer iu1               ! iu + 1
! Subprogram not used    integer iuc               ! index for log10(H2O continuum path)
! Subprogram not used    integer iuc1              ! iuc + 1
! Subprogram not used    integer ip                ! index for log10(P)
! Subprogram not used    integer ip1               ! ip + 1
! Subprogram not used    integer itp               ! index for T_p
! Subprogram not used    integer itp1              ! itp + 1
! Subprogram not used    integer ite               ! index for T_e
! Subprogram not used    integer ite1              ! ite + 1
! Subprogram not used    integer irh               ! index for RH
! Subprogram not used    integer irh1              ! irh + 1
! Subprogram not used 
! Subprogram not used    real(r8) dvar             ! normalized variation in T_p/T_e/P/U
! Subprogram not used    real(r8) uvar             ! U * diffusivity factor
! Subprogram not used    real(r8) uscl             ! factor for lineary scaling as U->0
! Subprogram not used 
! Subprogram not used    real(r8) wu               ! weight for U
! Subprogram not used    real(r8) wu1              ! 1 - wu
! Subprogram not used    real(r8) wuc              ! weight for H2O continuum path
! Subprogram not used    real(r8) wuc1             ! 1 - wuc
! Subprogram not used    real(r8) wp               ! weight for P
! Subprogram not used    real(r8) wp1              ! 1 - wp
! Subprogram not used    real(r8) wtp              ! weight for T_p
! Subprogram not used    real(r8) wtp1             ! 1 - wtp
! Subprogram not used    real(r8) wte              ! weight for T_e
! Subprogram not used    real(r8) wte1             ! 1 - wte
! Subprogram not used    real(r8) wrh              ! weight for RH
! Subprogram not used    real(r8) wrh1             ! 1 - wrh
! Subprogram not used 
! Subprogram not used    real(r8) w_0_0_           ! weight for Tp/Te combination
! Subprogram not used    real(r8) w_0_1_           ! weight for Tp/Te combination
! Subprogram not used    real(r8) w_1_0_           ! weight for Tp/Te combination
! Subprogram not used    real(r8) w_1_1_           ! weight for Tp/Te combination
! Subprogram not used 
! Subprogram not used    real(r8) w_0_00           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_0_01           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_0_10           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_0_11           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_00           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_01           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_10           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_11           ! weight for Tp/Te/RH combination
! Subprogram not used 
! Subprogram not used    real(r8) w00_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w00_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w00_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w00_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used 
! Subprogram not used    integer ib                ! spectral interval:
! Subprogram not used                              !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
! Subprogram not used                              !   2 = 800-1200 cm^-1
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    real(r8) pch2o            ! H2O continuum path
! Subprogram not used    real(r8) fch2o            ! temp. factor for continuum
! Subprogram not used    real(r8) uch2o            ! U corresponding to H2O cont. path (window)
! Subprogram not used 
! Subprogram not used    real(r8) fdif             ! secant(zenith angle) for diffusivity approx.
! Subprogram not used 
! Subprogram not used    real(r8) sslp_mks         ! Sea-level pressure in MKS units
! Subprogram not used    real(r8) esx              ! saturation vapor pressure returned by qsat
! Subprogram not used    real(r8) qsx              ! saturation mixing ratio returned by qsat
! Subprogram not used    real(r8) pnew_mks         ! pnew in MKS units
! Subprogram not used    real(r8) q_path           ! effective specific humidity along path
! Subprogram not used    real(r8) rh_path          ! effective relative humidity along path
! Subprogram not used 
! Subprogram not used       integer bnd_idx        ! LW band index
! Subprogram not used       real(r8) aer_pth_dlt   ! [kg m-2] STRAER path between interface levels k1 and k2
! Subprogram not used       real(r8) aer_pth_ngh(pcols)
! Subprogram not used                              ! [kg m-2] STRAER path between neighboring layers
! Subprogram not used       real(r8) odap_aer_ttl  ! [fraction] Total path absorption optical depth
! Subprogram not used       real(r8) aer_trn_ngh(pcols,nlwbands) 
! Subprogram not used                              ! [fraction] Total transmission between 
! Subprogram not used                              !            nearest neighbor sub-levels
! Subprogram not used !
! Subprogram not used !--------------------------Statement function---------------------------
! Subprogram not used !
! Subprogram not used    real(r8) dbvt,t             ! Planck fnctn tmp derivative for o3
! Subprogram not used !
! Subprogram not used    dbvt(t)=(-2.8911366682e-4_r8+(2.3771251896e-6_r8+1.1305188929e-10_r8*t)*t)/ &
! Subprogram not used       (1.0_r8+(-6.1364820707e-3_r8+1.5550319767e-5_r8*t)*t)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Initialize
! Subprogram not used !
! Subprogram not used    do k2=1,4
! Subprogram not used       do k1=1,ntoplw-1
! Subprogram not used          absnxt(:,k1,k2) = posinf    ! set unused portions for lf95 restart write
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    do k=ntoplw,pverp
! Subprogram not used       abstot(:,k,k) = posinf         ! set unused portions for lf95 restart write
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    do k=ntoplw,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          dbvtly(i,k) = dbvt(tlayr(i,k+1))
! Subprogram not used          dbvtit(i,k) = dbvt(tint(i,k))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used    rmw = amd/amco2
! Subprogram not used    do i=1,ncol
! Subprogram not used       dbvtit(i,pverp) = dbvt(tint(i,pverp))
! Subprogram not used       co2vmr(i) = co2mmr(i) * rmw
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    r293    = 1._r8/293._r8
! Subprogram not used    r250    = 1._r8/250._r8
! Subprogram not used    r3205   = 1._r8/.3205_r8
! Subprogram not used    r300    = 1._r8/300._r8
! Subprogram not used    rsslp   = 1._r8/sslp
! Subprogram not used    r2sslp  = 1._r8/(2._r8*sslp)
! Subprogram not used !
! Subprogram not used !Constants for computing U corresponding to H2O cont. path
! Subprogram not used !
! Subprogram not used    fdif       = 1.66_r8
! Subprogram not used    sslp_mks   = sslp / 10.0_r8
! Subprogram not used !
! Subprogram not used ! Non-adjacent layer absorptivity:
! Subprogram not used !
! Subprogram not used ! abso(i,1)     0 -  800 cm-1   h2o rotation band
! Subprogram not used ! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! Subprogram not used ! abso(i,2)   800 - 1200 cm-1   h2o window
! Subprogram not used !
! Subprogram not used ! Separation between rotation and vibration-rotation dropped, so
! Subprogram not used !                only 2 slots needed for H2O absorptivity
! Subprogram not used !
! Subprogram not used ! 500-800 cm^-1 H2o continuum/line overlap already included
! Subprogram not used !                in abso(i,1).  This used to be in abso(i,4)
! Subprogram not used !
! Subprogram not used ! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! Subprogram not used ! abso(i,4)   co2 15  micrometer band system
! Subprogram not used !
! Subprogram not used    do k=ntoplw,pverp
! Subprogram not used       do i=1,ncol
! Subprogram not used          pnmsq(i,k) = pnm(i,k)**2
! Subprogram not used          dtx(i) = tplnka(i,k) - 250._r8
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Non-nearest layer level loops
! Subprogram not used !
! Subprogram not used    do k1=pverp,ntoplw,-1
! Subprogram not used       do k2=pverp,ntoplw,-1
! Subprogram not used          if (k1 == k2) cycle
! Subprogram not used          do i=1,ncol
! Subprogram not used             dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
! Subprogram not used             u(i)      = abs(dplh2o(i))
! Subprogram not used             sqrtu(i)  = sqrt(u(i))
! Subprogram not used             ds2c      = abs(s2c(i,k1) - s2c(i,k2))
! Subprogram not used             dw(i)     = abs(w(i,k1) - w(i,k2))
! Subprogram not used             uc1(i)    = (ds2c + 1.7e-3_r8*u(i))*(1._r8 +  2._r8*ds2c)/(1._r8 + 15._r8*ds2c)
! Subprogram not used             pch2o     = ds2c
! Subprogram not used             pnew(i)   = u(i)/dw(i)
! Subprogram not used             pnew_mks  = pnew(i) * sslp_mks
! Subprogram not used !
! Subprogram not used ! Changed effective path temperature to std. Curtis-Godson form
! Subprogram not used !
! Subprogram not used             tpatha = abs(tcg(i,k1) - tcg(i,k2))/dw(i)
! Subprogram not used             t_p = min(max(tpatha, min_tp_h2o), max_tp_h2o)
! Subprogram not used 
! Subprogram not used             call qsat_water(t_p, pnew_mks, esx, qsx)
! Subprogram not used !
! Subprogram not used ! Compute effective RH along path
! Subprogram not used !
! Subprogram not used             q_path = dw(i) / abs(pnm(i,k1) - pnm(i,k2)) / rga
! Subprogram not used !
! Subprogram not used ! Calculate effective u, pnew for each band using
! Subprogram not used !        Hulst-Curtis-Godson approximation:
! Subprogram not used ! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
! Subprogram not used !           2nd edition, Oxford University Press, 1989.
! Subprogram not used ! Effective H2O path (w)
! Subprogram not used !      eq. 6.24, p. 228
! Subprogram not used ! Effective H2O path pressure (pnew = u/w):
! Subprogram not used !      eq. 6.29, p. 228
! Subprogram not used !
! Subprogram not used             ub(1) = abs(plh2ob(1,i,k1) - plh2ob(1,i,k2)) / psi(t_p,1)
! Subprogram not used             ub(2) = abs(plh2ob(2,i,k1) - plh2ob(2,i,k2)) / psi(t_p,2)
! Subprogram not used             
! Subprogram not used             pnewb(1) = ub(1) / abs(wb(1,i,k1) - wb(1,i,k2)) * phi(t_p,1)
! Subprogram not used             pnewb(2) = ub(2) / abs(wb(2,i,k1) - wb(2,i,k2)) * phi(t_p,2)
! Subprogram not used 
! Subprogram not used             dtx(i)      = tplnka(i,k2) - 250._r8
! Subprogram not used             dty(i)      = tpatha       - 250._r8
! Subprogram not used 
! Subprogram not used             fwk(i)  = fwcoef + fwc1/(1._r8 + fwc2*u(i))
! Subprogram not used             fwku(i) = fwk(i)*u(i)
! Subprogram not used !
! Subprogram not used ! Define variables for C/H/E (now C/LT/E) fit
! Subprogram not used !
! Subprogram not used ! abso(i,1)     0 -  800 cm-1   h2o rotation band
! Subprogram not used ! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! Subprogram not used ! abso(i,2)   800 - 1200 cm-1   h2o window
! Subprogram not used !
! Subprogram not used ! Separation between rotation and vibration-rotation dropped, so
! Subprogram not used !                only 2 slots needed for H2O absorptivity
! Subprogram not used !
! Subprogram not used ! Notation:
! Subprogram not used ! U   = integral (P/P_0 dW)  
! Subprogram not used ! P   = atmospheric pressure
! Subprogram not used ! P_0 = reference atmospheric pressure
! Subprogram not used ! W   = precipitable water path
! Subprogram not used ! T_e = emission temperature
! Subprogram not used ! T_p = path temperature
! Subprogram not used ! RH  = path relative humidity
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Terms for asymptotic value of emissivity
! Subprogram not used !
! Subprogram not used             te1  = tplnka(i,k2)
! Subprogram not used             te2  = te1 * te1
! Subprogram not used             te3  = te2 * te1
! Subprogram not used             te4  = te3 * te1
! Subprogram not used             te5  = te4 * te1
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !  Band-independent indices for lines and continuum tables
! Subprogram not used !
! Subprogram not used             dvar = (t_p - min_tp_h2o) / dtp_h2o
! Subprogram not used             itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
! Subprogram not used             itp1 = itp + 1
! Subprogram not used             wtp = dvar - floor(dvar)
! Subprogram not used             wtp1 = 1.0_r8 - wtp
! Subprogram not used             
! Subprogram not used             t_e = min(max(tplnka(i,k2)-t_p, min_te_h2o), max_te_h2o)
! Subprogram not used             dvar = (t_e - min_te_h2o) / dte_h2o
! Subprogram not used             ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
! Subprogram not used             ite1 = ite + 1
! Subprogram not used             wte = dvar - floor(dvar)
! Subprogram not used             wte1 = 1.0_r8 - wte
! Subprogram not used             
! Subprogram not used             rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
! Subprogram not used             dvar = (rh_path - min_rh_h2o) / drh_h2o
! Subprogram not used             irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
! Subprogram not used             irh1 = irh + 1
! Subprogram not used             wrh = dvar - floor(dvar)
! Subprogram not used             wrh1 = 1.0_r8 - wrh
! Subprogram not used 
! Subprogram not used             w_0_0_ = wtp  * wte
! Subprogram not used             w_0_1_ = wtp  * wte1
! Subprogram not used             w_1_0_ = wtp1 * wte 
! Subprogram not used             w_1_1_ = wtp1 * wte1
! Subprogram not used             
! Subprogram not used             w_0_00 = w_0_0_ * wrh
! Subprogram not used             w_0_01 = w_0_0_ * wrh1
! Subprogram not used             w_0_10 = w_0_1_ * wrh
! Subprogram not used             w_0_11 = w_0_1_ * wrh1
! Subprogram not used             w_1_00 = w_1_0_ * wrh
! Subprogram not used             w_1_01 = w_1_0_ * wrh1
! Subprogram not used             w_1_10 = w_1_1_ * wrh
! Subprogram not used             w_1_11 = w_1_1_ * wrh1
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! H2O Continuum path for 0-800 and 1200-2200 cm^-1
! Subprogram not used !
! Subprogram not used !    Assume foreign continuum dominates total H2O continuum in these bands
! Subprogram not used !    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
! Subprogram not used !    Then the effective H2O path is just 
! Subprogram not used !         U_c = integral[ f(P) dW ]
! Subprogram not used !    where 
! Subprogram not used !           W = water-vapor mass and 
! Subprogram not used !        f(P) = dependence of foreign continuum on pressure 
! Subprogram not used !             = P / sslp
! Subprogram not used !    Then 
! Subprogram not used !         U_c = U (the same effective H2O path as for lines)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Continuum terms for 800-1200 cm^-1
! Subprogram not used !
! Subprogram not used !    Assume self continuum dominates total H2O continuum for this band
! Subprogram not used !    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
! Subprogram not used !    Then the effective H2O self-continuum path is 
! Subprogram not used !         U_c = integral[ h(e,T) dW ]                        (*eq. 1*)
! Subprogram not used !    where 
! Subprogram not used !           W = water-vapor mass and 
! Subprogram not used !           e = partial pressure of H2O along path
! Subprogram not used !           T = temperature along path
! Subprogram not used !      h(e,T) = dependence of foreign continuum on e,T
! Subprogram not used !             = e / sslp * f(T)
! Subprogram not used !
! Subprogram not used !    Replacing
! Subprogram not used !           e =~ q * P / epsilo
! Subprogram not used !           q = mixing ratio of H2O
! Subprogram not used !     epsilo = 0.622
! Subprogram not used !
! Subprogram not used !    and using the definition
! Subprogram not used !           U = integral [ (P / sslp) dW ]
! Subprogram not used !             = (P / sslp) W                                 (homogeneous path)
! Subprogram not used !
! Subprogram not used !    the effective path length for the self continuum is
! Subprogram not used !         U_c = (q / epsilo) f(T) U                         (*eq. 2*)
! Subprogram not used !
! Subprogram not used !    Once values of T, U, and q have been calculated for the inhomogeneous
! Subprogram not used !        path, this sets U_c for the corresponding
! Subprogram not used !        homogeneous atmosphere.  However, this need not equal the
! Subprogram not used !        value of U_c' defined by eq. 1 for the actual inhomogeneous atmosphere
! Subprogram not used !        under consideration.
! Subprogram not used !
! Subprogram not used !    Solution: hold T and q constant, solve for U' that gives U_c' by
! Subprogram not used !        inverting eq. (2):
! Subprogram not used !
! Subprogram not used !        U' = (U_c * epsilo) / (q * f(T))
! Subprogram not used !
! Subprogram not used             fch2o = fh2oself(t_p) 
! Subprogram not used             uch2o = (pch2o * epsilo) / (q_path * fch2o)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Band-dependent indices for non-window
! Subprogram not used !
! Subprogram not used             ib = 1
! Subprogram not used 
! Subprogram not used             uvar = ub(ib) * fdif
! Subprogram not used             log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
! Subprogram not used             dvar = (log_u - min_lu_h2o) / dlu_h2o
! Subprogram not used             iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used             iu1 = iu + 1
! Subprogram not used             wu = dvar - floor(dvar)
! Subprogram not used             wu1 = 1.0_r8 - wu
! Subprogram not used             
! Subprogram not used             log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
! Subprogram not used             dvar = (log_p - min_lp_h2o) / dlp_h2o
! Subprogram not used             ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
! Subprogram not used             ip1 = ip + 1
! Subprogram not used             wp = dvar - floor(dvar)
! Subprogram not used             wp1 = 1.0_r8 - wp
! Subprogram not used          
! Subprogram not used             w00_00 = wp  * w_0_00 
! Subprogram not used             w00_01 = wp  * w_0_01 
! Subprogram not used             w00_10 = wp  * w_0_10 
! Subprogram not used             w00_11 = wp  * w_0_11 
! Subprogram not used             w01_00 = wp  * w_1_00 
! Subprogram not used             w01_01 = wp  * w_1_01 
! Subprogram not used             w01_10 = wp  * w_1_10 
! Subprogram not used             w01_11 = wp  * w_1_11 
! Subprogram not used             w10_00 = wp1 * w_0_00 
! Subprogram not used             w10_01 = wp1 * w_0_01 
! Subprogram not used             w10_10 = wp1 * w_0_10 
! Subprogram not used             w10_11 = wp1 * w_0_11 
! Subprogram not used             w11_00 = wp1 * w_1_00 
! Subprogram not used             w11_01 = wp1 * w_1_01 
! Subprogram not used             w11_10 = wp1 * w_1_10 
! Subprogram not used             w11_11 = wp1 * w_1_11 
! Subprogram not used !
! Subprogram not used ! Asymptotic value of absorptivity as U->infinity
! Subprogram not used !
! Subprogram not used             fa = fat(1,ib) + &
! Subprogram not used                  fat(2,ib) * te1 + &
! Subprogram not used                  fat(3,ib) * te2 + &
! Subprogram not used                  fat(4,ib) * te3 + &
! Subprogram not used                  fat(5,ib) * te4 + &
! Subprogram not used                  fat(6,ib) * te5
! Subprogram not used 
! Subprogram not used             a_star = &
! Subprogram not used                  ah2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
! Subprogram not used             abso(i,ib) = min(max(fa * (1.0_r8 - (1.0_r8 - a_star) * &
! Subprogram not used                                  aer_trn_ttl(i,k1,k2,ib)), &
! Subprogram not used                              0.0_r8), 1.0_r8)
! Subprogram not used !
! Subprogram not used ! Invoke linear limit for scaling wrt u below min_u_h2o
! Subprogram not used !
! Subprogram not used             if (uvar < min_u_h2o) then
! Subprogram not used                uscl = uvar / min_u_h2o
! Subprogram not used                abso(i,ib) = abso(i,ib) * uscl
! Subprogram not used             endif
! Subprogram not used                          
! Subprogram not used !
! Subprogram not used ! Band-dependent indices for window
! Subprogram not used !
! Subprogram not used             ib = 2
! Subprogram not used 
! Subprogram not used             uvar = ub(ib) * fdif
! Subprogram not used             log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
! Subprogram not used             dvar = (log_u - min_lu_h2o) / dlu_h2o
! Subprogram not used             iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used             iu1 = iu + 1
! Subprogram not used             wu = dvar - floor(dvar)
! Subprogram not used             wu1 = 1.0_r8 - wu
! Subprogram not used             
! Subprogram not used             log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
! Subprogram not used             dvar = (log_p - min_lp_h2o) / dlp_h2o
! Subprogram not used             ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
! Subprogram not used             ip1 = ip + 1
! Subprogram not used             wp = dvar - floor(dvar)
! Subprogram not used             wp1 = 1.0_r8 - wp
! Subprogram not used          
! Subprogram not used             w00_00 = wp  * w_0_00 
! Subprogram not used             w00_01 = wp  * w_0_01 
! Subprogram not used             w00_10 = wp  * w_0_10 
! Subprogram not used             w00_11 = wp  * w_0_11 
! Subprogram not used             w01_00 = wp  * w_1_00 
! Subprogram not used             w01_01 = wp  * w_1_01 
! Subprogram not used             w01_10 = wp  * w_1_10 
! Subprogram not used             w01_11 = wp  * w_1_11 
! Subprogram not used             w10_00 = wp1 * w_0_00 
! Subprogram not used             w10_01 = wp1 * w_0_01 
! Subprogram not used             w10_10 = wp1 * w_0_10 
! Subprogram not used             w10_11 = wp1 * w_0_11 
! Subprogram not used             w11_00 = wp1 * w_1_00 
! Subprogram not used             w11_01 = wp1 * w_1_01 
! Subprogram not used             w11_10 = wp1 * w_1_10 
! Subprogram not used             w11_11 = wp1 * w_1_11 
! Subprogram not used 
! Subprogram not used             log_uc  = min(log10(max(uch2o * fdif, min_u_h2o)), max_lu_h2o)
! Subprogram not used             dvar = (log_uc - min_lu_h2o) / dlu_h2o
! Subprogram not used             iuc = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used             iuc1 = iuc + 1
! Subprogram not used             wuc = dvar - floor(dvar)
! Subprogram not used             wuc1 = 1.0_r8 - wuc
! Subprogram not used !
! Subprogram not used ! Asymptotic value of absorptivity as U->infinity
! Subprogram not used !
! Subprogram not used             fa = fat(1,ib) + &
! Subprogram not used                  fat(2,ib) * te1 + &
! Subprogram not used                  fat(3,ib) * te2 + &
! Subprogram not used                  fat(4,ib) * te3 + &
! Subprogram not used                  fat(5,ib) * te4 + &
! Subprogram not used                  fat(6,ib) * te5
! Subprogram not used 
! Subprogram not used             l_star = &
! Subprogram not used                  ln_ah2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
! Subprogram not used                  ln_ah2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
! Subprogram not used                  ln_ah2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
! Subprogram not used 
! Subprogram not used             c_star = &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc , ite , irh ) * w11_11 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc , ite , irh1) * w11_10 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc , ite1, irh ) * w11_01 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc , ite1, irh1) * w11_00 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc1, ite , irh ) * w11_11 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc1, ite , irh1) * w11_10 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc1, ite1, irh ) * w11_01 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp , iuc1, ite1, irh1) * w11_00 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc , ite , irh ) * w10_11 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc , ite , irh1) * w10_10 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc , ite1, irh ) * w10_01 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc , ite1, irh1) * w10_00 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc1, ite , irh ) * w10_11 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc1, ite , irh1) * w10_10 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc1, ite1, irh ) * w10_01 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip , itp1, iuc1, ite1, irh1) * w10_00 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc , ite , irh ) * w01_11 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc , ite , irh1) * w01_10 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc , ite1, irh ) * w01_01 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc , ite1, irh1) * w01_00 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc1, ite , irh ) * w01_11 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc1, ite , irh1) * w01_10 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc1, ite1, irh ) * w01_01 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp , iuc1, ite1, irh1) * w01_00 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc , ite , irh ) * w00_11 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc , ite , irh1) * w00_10 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc , ite1, irh ) * w00_01 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc , ite1, irh1) * w00_00 * wuc1 + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc1, ite , irh ) * w00_11 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc1, ite , irh1) * w00_10 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc1, ite1, irh ) * w00_01 * wuc  + &
! Subprogram not used                  cn_ah2ow(ip1, itp1, iuc1, ite1, irh1) * w00_00 * wuc 
! Subprogram not used             abso(i,ib) = min(max(fa * (1.0_r8 - l_star * c_star * &
! Subprogram not used                                  aer_trn_ttl(i,k1,k2,ib)), &
! Subprogram not used                              0.0_r8), 1.0_r8) 
! Subprogram not used !
! Subprogram not used ! Invoke linear limit for scaling wrt u below min_u_h2o
! Subprogram not used !
! Subprogram not used             if (uvar < min_u_h2o) then
! Subprogram not used                uscl = uvar / min_u_h2o
! Subprogram not used                abso(i,ib) = abso(i,ib) * uscl
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1._r8 + c16*dty(i))
! Subprogram not used             term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1._r8 + c17*dty(i))
! Subprogram not used             term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1._r8 + c26*dty(i))
! Subprogram not used             term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1._r8 + c27*dty(i))
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! 500 -  800 cm-1   h2o rotation band overlap with co2
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             k21    = term7(i,1) + term8(i,1)/ &
! Subprogram not used                (1._r8 + (c30 + c31*(dty(i)-10._r8)*(dty(i)-10._r8))*sqrtu(i))
! Subprogram not used             k22    = term7(i,2) + term8(i,2)/ &
! Subprogram not used                (1._r8 + (c28 + c29*(dty(i)-10._r8))*sqrtu(i))
! Subprogram not used             tr1    = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
! Subprogram not used             tr2    = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
! Subprogram not used             tr1=tr1*aer_trn_ttl(i,k1,k2,idx_LW_0650_0800) 
! Subprogram not used !                                          ! H2O line+STRAER trn 650--800 cm-1
! Subprogram not used             tr2=tr2*aer_trn_ttl(i,k1,k2,idx_LW_0500_0650)
! Subprogram not used !                                          ! H2O line+STRAER trn 500--650 cm-1
! Subprogram not used             tr5    = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
! Subprogram not used             tr6    = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
! Subprogram not used             tr9(i)   = tr1*tr5
! Subprogram not used             tr10(i)  = tr2*tr6
! Subprogram not used             th2o(i) = tr10(i)
! Subprogram not used             trab2(i) = 0.65_r8*tr9(i) + 0.35_r8*tr10(i)
! Subprogram not used          end do
! Subprogram not used          if (k2 < k1) then
! Subprogram not used             do i=1,ncol
! Subprogram not used                to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
! Subprogram not used             end do
! Subprogram not used          else
! Subprogram not used             do i=1,ncol
! Subprogram not used                to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
! Subprogram not used             end do
! Subprogram not used          end if
! Subprogram not used !
! Subprogram not used ! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             dpnm(i)  = pnm(i,k1) - pnm(i,k2)
! Subprogram not used             to3co2(i) = (pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/dpnm(i)
! Subprogram not used             te       = (to3co2(i)*r293)**.7_r8
! Subprogram not used             dplos    = plos(i,k1) - plos(i,k2)
! Subprogram not used             if (dplos == 0._r8) then
! Subprogram not used               abso(i,3) = 0._r8
! Subprogram not used               to3(i) = 1._r8
! Subprogram not used               write(iulog,*) 'radiation ozone error  ',k1,k2,plos(i,k1)
! Subprogram not used             else 
! Subprogram not used               dplol    = plol(i,k1) - plol(i,k2)
! Subprogram not used               u1       = 18.29_r8*abs(dplos)/te
! Subprogram not used               u2       = .5649_r8*abs(dplos)/te
! Subprogram not used               rphat    = dplol/dplos
! Subprogram not used               tlocal   = tint(i,k2)
! Subprogram not used               tcrfac   = sqrt(tlocal*r250)*te
! Subprogram not used               beta     = r3205*(rphat + dpfo3*tcrfac)
! Subprogram not used               realnu   = te/beta
! Subprogram not used               tmp1     = u1/sqrt(4._r8 + u1*(1._r8 + realnu))
! Subprogram not used               tmp2     = u2/sqrt(4._r8 + u2*(1._r8 + realnu))
! Subprogram not used               o3bndi    = 74._r8*te*log(1._r8 + tmp1 + tmp2)
! Subprogram not used               abso(i,3) = o3bndi*to3h2o(i)*dbvtit(i,k2)
! Subprogram not used               to3(i)   = 1.0_r8/(1._r8 + 0.1_r8*tmp1 + 0.1_r8*tmp2)
! Subprogram not used             endif
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! abso(i,4)      co2 15  micrometer band system
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             sqwp      = sqrt(abs(plco2(i,k1) - plco2(i,k2)))
! Subprogram not used             et        = exp(-480._r8/to3co2(i))
! Subprogram not used             sqti(i)   = sqrt(to3co2(i))
! Subprogram not used             rsqti     = 1._r8/sqti(i)
! Subprogram not used             et2       = et*et
! Subprogram not used             et4       = et2*et2
! Subprogram not used             omet      = 1._r8 - 1.5_r8*et2
! Subprogram not used             f1co2     = 899.70_r8*omet*(1._r8 + 1.94774_r8*et + 4.73486_r8*et2)*rsqti
! Subprogram not used             f1sqwp(i) = f1co2*sqwp
! Subprogram not used             t1co2(i)  = 1._r8/(1._r8 + (245.18_r8*omet*sqwp*rsqti))
! Subprogram not used             oneme     = 1._r8 - et2
! Subprogram not used             alphat    = oneme**3*rsqti
! Subprogram not used             pi        = abs(dpnm(i))
! Subprogram not used             wco2      =  2.5221_r8*co2vmr(i)*pi*rga
! Subprogram not used             u7(i)     =  4.9411e4_r8*alphat*et2*wco2
! Subprogram not used             u8        =  3.9744e4_r8*alphat*et4*wco2
! Subprogram not used             u9        =  1.0447e5_r8*alphat*et4*et2*wco2
! Subprogram not used             u13       = 2.8388e3_r8*alphat*et4*wco2
! Subprogram not used             tpath     = to3co2(i)
! Subprogram not used             tlocal    = tint(i,k2)
! Subprogram not used             tcrfac    = sqrt(tlocal*r250*tpath*r300)
! Subprogram not used             posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp + dpfco2*tcrfac)*rsqti
! Subprogram not used             rbeta7(i) = 1._r8/(5.3228_r8*posqt)
! Subprogram not used             rbeta8    = 1._r8/(10.6576_r8*posqt)
! Subprogram not used             rbeta9    = rbeta7(i)
! Subprogram not used             rbeta13   = rbeta9
! Subprogram not used             f2co2(i)  = (u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i)))) + &
! Subprogram not used                (u8   /sqrt(4._r8 + u8*(1._r8 + rbeta8))) + &
! Subprogram not used                (u9   /sqrt(4._r8 + u9*(1._r8 + rbeta9)))
! Subprogram not used             f3co2(i)  = u13/sqrt(4._r8 + u13*(1._r8 + rbeta13))
! Subprogram not used          end do
! Subprogram not used          if (k2 >= k1) then
! Subprogram not used             do i=1,ncol
! Subprogram not used                sqti(i) = sqrt(tlayr(i,k2))
! Subprogram not used             end do
! Subprogram not used          end if
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             tmp1      = log(1._r8 + f1sqwp(i))
! Subprogram not used             tmp2      = log(1._r8 + f2co2(i))
! Subprogram not used             tmp3      = log(1._r8 + f3co2(i))
! Subprogram not used             absbnd    = (tmp1 + 2._r8*t1co2(i)*tmp2 + 2._r8*tmp3)*sqti(i)
! Subprogram not used             abso(i,4) = trab2(i)*co2em(i,k2)*absbnd
! Subprogram not used             tco2(i)   = 1._r8/(1.0_r8+10.0_r8*(u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i)))))
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! Calculate absorptivity due to trace gases, abstrc
! Subprogram not used !
! Subprogram not used          call trcab(ncol     ,                            &
! Subprogram not used             k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
! Subprogram not used             un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
! Subprogram not used             uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
! Subprogram not used             bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
! Subprogram not used             s2c     ,uptype  ,u       ,abplnk1 ,tco2    , &
! Subprogram not used             th2o    ,to3     ,abstrc  , &
! Subprogram not used             aer_trn_ttl)
! Subprogram not used !
! Subprogram not used ! Sum total absorptivity
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             abstot(i,k1,k2) = abso(i,1) + abso(i,2) + &
! Subprogram not used                abso(i,3) + abso(i,4) + abstrc(i)
! Subprogram not used          end do
! Subprogram not used       end do ! do k2 = 
! Subprogram not used    end do ! do k1 = 
! Subprogram not used !
! Subprogram not used ! Adjacent layer absorptivity:
! Subprogram not used !
! Subprogram not used ! abso(i,1)     0 -  800 cm-1   h2o rotation band
! Subprogram not used ! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! Subprogram not used ! abso(i,2)   800 - 1200 cm-1   h2o window
! Subprogram not used !
! Subprogram not used ! Separation between rotation and vibration-rotation dropped, so
! Subprogram not used !                only 2 slots needed for H2O absorptivity
! Subprogram not used !
! Subprogram not used ! 500-800 cm^-1 H2o continuum/line overlap already included
! Subprogram not used !                in abso(i,1).  This used to be in abso(i,4)
! Subprogram not used !
! Subprogram not used ! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! Subprogram not used ! abso(i,4)   co2 15  micrometer band system
! Subprogram not used !
! Subprogram not used ! Nearest layer level loop
! Subprogram not used !
! Subprogram not used    do k2=pver,ntoplw,-1
! Subprogram not used       do i=1,ncol
! Subprogram not used          tbar(i,1)   = 0.5_r8*(tint(i,k2+1) + tlayr(i,k2+1))
! Subprogram not used          emm(i,1)    = 0.5_r8*(co2em(i,k2+1) + co2eml(i,k2))
! Subprogram not used          tbar(i,2)   = 0.5_r8*(tlayr(i,k2+1) + tint(i,k2))
! Subprogram not used          emm(i,2)    = 0.5_r8*(co2em(i,k2) + co2eml(i,k2))
! Subprogram not used          tbar(i,3)   = 0.5_r8*(tbar(i,2) + tbar(i,1))
! Subprogram not used          emm(i,3)    = emm(i,1)
! Subprogram not used          tbar(i,4)   = tbar(i,3)
! Subprogram not used          emm(i,4)    = emm(i,2)
! Subprogram not used          o3emm(i,1)  = 0.5_r8*(dbvtit(i,k2+1) + dbvtly(i,k2))
! Subprogram not used          o3emm(i,2)  = 0.5_r8*(dbvtit(i,k2) + dbvtly(i,k2))
! Subprogram not used          o3emm(i,3)  = o3emm(i,1)
! Subprogram not used          o3emm(i,4)  = o3emm(i,2)
! Subprogram not used          temh2o(i,1) = tbar(i,1)
! Subprogram not used          temh2o(i,2) = tbar(i,2)
! Subprogram not used          temh2o(i,3) = tbar(i,1)
! Subprogram not used          temh2o(i,4) = tbar(i,2)
! Subprogram not used          dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used !  Weighted Planck functions for trace gases
! Subprogram not used !
! Subprogram not used       do wvl = 1,14
! Subprogram not used          do i = 1,ncol
! Subprogram not used             bplnk(wvl,i,1) = 0.5_r8*(abplnk1(wvl,i,k2+1) + abplnk2(wvl,i,k2))
! Subprogram not used             bplnk(wvl,i,2) = 0.5_r8*(abplnk1(wvl,i,k2) + abplnk2(wvl,i,k2))
! Subprogram not used             bplnk(wvl,i,3) = bplnk(wvl,i,1)
! Subprogram not used             bplnk(wvl,i,4) = bplnk(wvl,i,2)
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used       
! Subprogram not used       do i=1,ncol
! Subprogram not used          rdpnmsq    = 1._r8/(pnmsq(i,k2+1) - pnmsq(i,k2))
! Subprogram not used          rdpnm      = 1._r8/dpnm(i)
! Subprogram not used          p1         = .5_r8*(pbr(i,k2) + pnm(i,k2+1))
! Subprogram not used          p2         = .5_r8*(pbr(i,k2) + pnm(i,k2  ))
! Subprogram not used          uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
! Subprogram not used          uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
! Subprogram not used          uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
! Subprogram not used          uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
! Subprogram not used          winpl(i,1) = (.5_r8*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
! Subprogram not used          winpl(i,2) = (.5_r8*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
! Subprogram not used          winpl(i,3) = (.5_r8*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*rdpnm
! Subprogram not used          winpl(i,4) = (.5_r8*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*rdpnm
! Subprogram not used          tmp1       = 1._r8/(piln(i,k2+1) - piln(i,k2))
! Subprogram not used          tmp2       = piln(i,k2+1) - pmln(i,k2)
! Subprogram not used          tmp3       = piln(i,k2  ) - pmln(i,k2)
! Subprogram not used          zinpl(i,1) = (.5_r8*tmp2          )*tmp1
! Subprogram not used          zinpl(i,2) = (        - .5_r8*tmp3)*tmp1
! Subprogram not used          zinpl(i,3) = (.5_r8*tmp2 -    tmp3)*tmp1
! Subprogram not used          zinpl(i,4) = (   tmp2 - .5_r8*tmp3)*tmp1
! Subprogram not used          pinpl(i,1) = 0.5_r8*(p1 + pnm(i,k2+1))
! Subprogram not used          pinpl(i,2) = 0.5_r8*(p2 + pnm(i,k2  ))
! Subprogram not used          pinpl(i,3) = 0.5_r8*(p1 + pnm(i,k2  ))
! Subprogram not used          pinpl(i,4) = 0.5_r8*(p2 + pnm(i,k2+1))
! Subprogram not used       end do
! Subprogram not used       do kn=1,4
! Subprogram not used          do i=1,ncol
! Subprogram not used             u(i)     = uinpl(i,kn)*abs(plh2o(i,k2) - plh2o(i,k2+1))
! Subprogram not used             sqrtu(i) = sqrt(u(i))
! Subprogram not used             dw(i)    = abs(w(i,k2) - w(i,k2+1))
! Subprogram not used             pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
! Subprogram not used             pnew_mks  = pnew(i) * sslp_mks
! Subprogram not used             t_p = min(max(tbar(i,kn), min_tp_h2o), max_tp_h2o)
! Subprogram not used 
! Subprogram not used             call qsat_water(t_p, pnew_mks, esx, qsx)
! Subprogram not used 
! Subprogram not used             q_path = dw(i) / ABS(dpnm(i)) / rga
! Subprogram not used             
! Subprogram not used             ds2c     = abs(s2c(i,k2) - s2c(i,k2+1))
! Subprogram not used             uc1(i)   = uinpl(i,kn)*ds2c
! Subprogram not used             pch2o    = uc1(i)
! Subprogram not used             uc1(i)   = (uc1(i) + 1.7e-3_r8*u(i))*(1._r8 +  2._r8*uc1(i))/(1._r8 + 15._r8*uc1(i))
! Subprogram not used             dtx(i)      = temh2o(i,kn) - 250._r8
! Subprogram not used             dty(i)      = tbar(i,kn) - 250._r8
! Subprogram not used             
! Subprogram not used             fwk(i)    = fwcoef + fwc1/(1._r8 + fwc2*u(i))
! Subprogram not used             fwku(i)   = fwk(i)*u(i)
! Subprogram not used 
! Subprogram not used             aer_trn_ngh(i, 1:nlwbands)= &
! Subprogram not used               exp(-fdif * uinpl(i,kn) * odap_aer(i, k2, 1:nlwbands ) )
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Define variables for C/H/E (now C/LT/E) fit
! Subprogram not used !
! Subprogram not used ! abso(i,1)     0 -  800 cm-1   h2o rotation band
! Subprogram not used ! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! Subprogram not used ! abso(i,2)   800 - 1200 cm-1   h2o window
! Subprogram not used !
! Subprogram not used ! Separation between rotation and vibration-rotation dropped, so
! Subprogram not used !                only 2 slots needed for H2O absorptivity
! Subprogram not used !
! Subprogram not used ! Notation:
! Subprogram not used ! U   = integral (P/P_0 dW)  
! Subprogram not used ! P   = atmospheric pressure
! Subprogram not used ! P_0 = reference atmospheric pressure
! Subprogram not used ! W   = precipitable water path
! Subprogram not used ! T_e = emission temperature
! Subprogram not used ! T_p = path temperature
! Subprogram not used ! RH  = path relative humidity
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Terms for asymptotic value of emissivity
! Subprogram not used !
! Subprogram not used             te1  = temh2o(i,kn)
! Subprogram not used             te2  = te1 * te1
! Subprogram not used             te3  = te2 * te1
! Subprogram not used             te4  = te3 * te1
! Subprogram not used             te5  = te4 * te1
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Indices for lines and continuum tables 
! Subprogram not used ! Note: because we are dealing with the nearest layer,
! Subprogram not used !       the Hulst-Curtis-Godson corrections
! Subprogram not used !       for inhomogeneous paths are not applied.
! Subprogram not used !
! Subprogram not used             uvar = u(i)*fdif
! Subprogram not used             log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
! Subprogram not used             dvar = (log_u - min_lu_h2o) / dlu_h2o
! Subprogram not used             iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used             iu1 = iu + 1
! Subprogram not used             wu = dvar - floor(dvar)
! Subprogram not used             wu1 = 1.0_r8 - wu
! Subprogram not used             
! Subprogram not used             log_p  = min(log10(max(pnew(i), min_p_h2o)), max_lp_h2o)
! Subprogram not used             dvar = (log_p - min_lp_h2o) / dlp_h2o
! Subprogram not used             ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
! Subprogram not used             ip1 = ip + 1
! Subprogram not used             wp = dvar - floor(dvar)
! Subprogram not used             wp1 = 1.0_r8 - wp
! Subprogram not used             
! Subprogram not used             dvar = (t_p - min_tp_h2o) / dtp_h2o
! Subprogram not used             itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
! Subprogram not used             itp1 = itp + 1
! Subprogram not used             wtp = dvar - floor(dvar)
! Subprogram not used             wtp1 = 1.0_r8 - wtp
! Subprogram not used             
! Subprogram not used             t_e = min(max(temh2o(i,kn)-t_p,min_te_h2o),max_te_h2o)
! Subprogram not used             dvar = (t_e - min_te_h2o) / dte_h2o
! Subprogram not used             ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
! Subprogram not used             ite1 = ite + 1
! Subprogram not used             wte = dvar - floor(dvar)
! Subprogram not used             wte1 = 1.0_r8 - wte
! Subprogram not used             
! Subprogram not used             rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
! Subprogram not used             dvar = (rh_path - min_rh_h2o) / drh_h2o
! Subprogram not used             irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
! Subprogram not used             irh1 = irh + 1
! Subprogram not used             wrh = dvar - floor(dvar)
! Subprogram not used             wrh1 = 1.0_r8 - wrh
! Subprogram not used             
! Subprogram not used             w_0_0_ = wtp  * wte
! Subprogram not used             w_0_1_ = wtp  * wte1
! Subprogram not used             w_1_0_ = wtp1 * wte 
! Subprogram not used             w_1_1_ = wtp1 * wte1
! Subprogram not used             
! Subprogram not used             w_0_00 = w_0_0_ * wrh
! Subprogram not used             w_0_01 = w_0_0_ * wrh1
! Subprogram not used             w_0_10 = w_0_1_ * wrh
! Subprogram not used             w_0_11 = w_0_1_ * wrh1
! Subprogram not used             w_1_00 = w_1_0_ * wrh
! Subprogram not used             w_1_01 = w_1_0_ * wrh1
! Subprogram not used             w_1_10 = w_1_1_ * wrh
! Subprogram not used             w_1_11 = w_1_1_ * wrh1
! Subprogram not used             
! Subprogram not used             w00_00 = wp  * w_0_00 
! Subprogram not used             w00_01 = wp  * w_0_01 
! Subprogram not used             w00_10 = wp  * w_0_10 
! Subprogram not used             w00_11 = wp  * w_0_11 
! Subprogram not used             w01_00 = wp  * w_1_00 
! Subprogram not used             w01_01 = wp  * w_1_01 
! Subprogram not used             w01_10 = wp  * w_1_10 
! Subprogram not used             w01_11 = wp  * w_1_11 
! Subprogram not used             w10_00 = wp1 * w_0_00 
! Subprogram not used             w10_01 = wp1 * w_0_01 
! Subprogram not used             w10_10 = wp1 * w_0_10 
! Subprogram not used             w10_11 = wp1 * w_0_11 
! Subprogram not used             w11_00 = wp1 * w_1_00 
! Subprogram not used             w11_01 = wp1 * w_1_01 
! Subprogram not used             w11_10 = wp1 * w_1_10 
! Subprogram not used             w11_11 = wp1 * w_1_11 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Non-window absorptivity
! Subprogram not used !
! Subprogram not used             ib = 1
! Subprogram not used             
! Subprogram not used             fa = fat(1,ib) + &
! Subprogram not used                  fat(2,ib) * te1 + &
! Subprogram not used                  fat(3,ib) * te2 + &
! Subprogram not used                  fat(4,ib) * te3 + &
! Subprogram not used                  fat(5,ib) * te4 + &
! Subprogram not used                  fat(6,ib) * te5
! Subprogram not used             
! Subprogram not used             a_star = &
! Subprogram not used                  ah2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
! Subprogram not used                  ah2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
! Subprogram not used                  ah2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
! Subprogram not used                  ah2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu
! Subprogram not used             
! Subprogram not used             abso(i,ib) = min(max(fa * (1.0_r8 - (1.0_r8 - a_star) * &
! Subprogram not used                                  aer_trn_ngh(i,ib)), &
! Subprogram not used                              0.0_r8), 1.0_r8)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Invoke linear limit for scaling wrt u below min_u_h2o
! Subprogram not used !
! Subprogram not used             if (uvar < min_u_h2o) then
! Subprogram not used                uscl = uvar / min_u_h2o
! Subprogram not used                abso(i,ib) = abso(i,ib) * uscl
! Subprogram not used             endif
! Subprogram not used             
! Subprogram not used !
! Subprogram not used ! Window absorptivity
! Subprogram not used !
! Subprogram not used             ib = 2
! Subprogram not used             
! Subprogram not used             fa = fat(1,ib) + &
! Subprogram not used                  fat(2,ib) * te1 + &
! Subprogram not used                  fat(3,ib) * te2 + &
! Subprogram not used                  fat(4,ib) * te3 + &
! Subprogram not used                  fat(5,ib) * te4 + &
! Subprogram not used                  fat(6,ib) * te5
! Subprogram not used             
! Subprogram not used             a_star = &
! Subprogram not used                  ah2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
! Subprogram not used                  ah2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
! Subprogram not used                  ah2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
! Subprogram not used                  ah2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
! Subprogram not used                  ah2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
! Subprogram not used                  ah2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
! Subprogram not used                  ah2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
! Subprogram not used                  ah2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
! Subprogram not used                  ah2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
! Subprogram not used                  ah2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
! Subprogram not used                  ah2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu
! Subprogram not used             
! Subprogram not used             abso(i,ib) = min(max(fa * (1.0_r8 - (1.0_r8 - a_star) * &
! Subprogram not used                                  aer_trn_ngh(i,ib)), &
! Subprogram not used                              0.0_r8), 1.0_r8)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Invoke linear limit for scaling wrt u below min_u_h2o
! Subprogram not used !
! Subprogram not used             if (uvar < min_u_h2o) then
! Subprogram not used                uscl = uvar / min_u_h2o
! Subprogram not used                abso(i,ib) = abso(i,ib) * uscl
! Subprogram not used             endif
! Subprogram not used             
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1._r8 + c16*dty(i))
! Subprogram not used             term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1._r8 + c17*dty(i))
! Subprogram not used             term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1._r8 + c26*dty(i))
! Subprogram not used             term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1._r8 + c27*dty(i))
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! 500 -  800 cm-1   h2o rotation band overlap with co2
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             dtym10     = dty(i) - 10._r8
! Subprogram not used             denom      = 1._r8 + (c30 + c31*dtym10*dtym10)*sqrtu(i)
! Subprogram not used             k21        = term7(i,1) + term8(i,1)/denom
! Subprogram not used             denom      = 1._r8 + (c28 + c29*dtym10       )*sqrtu(i)
! Subprogram not used             k22        = term7(i,2) + term8(i,2)/denom
! Subprogram not used             tr1     = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
! Subprogram not used             tr2     = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
! Subprogram not used             tr1=tr1*aer_trn_ngh(i,idx_LW_0650_0800) 
! Subprogram not used !                                         ! H2O line+STRAER trn 650--800 cm-1
! Subprogram not used             tr2=tr2*aer_trn_ngh(i,idx_LW_0500_0650) 
! Subprogram not used !                                         ! H2O line+STRAER trn 500--650 cm-1
! Subprogram not used             tr5     = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
! Subprogram not used             tr6     = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
! Subprogram not used             tr9(i)  = tr1*tr5
! Subprogram not used             tr10(i) = tr2*tr6
! Subprogram not used             trab2(i)= 0.65_r8*tr9(i) + 0.35_r8*tr10(i)
! Subprogram not used             th2o(i) = tr10(i)
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! abso(i,3)  o3  9.6 micrometer (nu3 and nu1 bands)
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             te        = (tbar(i,kn)*r293)**.7_r8
! Subprogram not used             dplos     = abs(plos(i,k2+1) - plos(i,k2))
! Subprogram not used             u1        = zinpl(i,kn)*18.29_r8*dplos/te
! Subprogram not used             u2        = zinpl(i,kn)*.5649_r8*dplos/te
! Subprogram not used             tlocal    = tbar(i,kn)
! Subprogram not used             tcrfac    = sqrt(tlocal*r250)*te
! Subprogram not used             beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
! Subprogram not used             realnu    = te/beta
! Subprogram not used             tmp1      = u1/sqrt(4._r8 + u1*(1._r8 + realnu))
! Subprogram not used             tmp2      = u2/sqrt(4._r8 + u2*(1._r8 + realnu))
! Subprogram not used             o3bndi    = 74._r8*te*log(1._r8 + tmp1 + tmp2)
! Subprogram not used             abso(i,3) = o3bndi*o3emm(i,kn)*(h2otr(i,k2+1)/h2otr(i,k2))
! Subprogram not used             to3(i)    = 1.0_r8/(1._r8 + 0.1_r8*tmp1 + 0.1_r8*tmp2)
! Subprogram not used          end do
! Subprogram not used !
! Subprogram not used ! abso(i,4)   co2 15  micrometer band system
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             dplco2   = plco2(i,k2+1) - plco2(i,k2)
! Subprogram not used             sqwp     = sqrt(uinpl(i,kn)*dplco2)
! Subprogram not used             et       = exp(-480._r8/tbar(i,kn))
! Subprogram not used             sqti(i)  = sqrt(tbar(i,kn))
! Subprogram not used             rsqti    = 1._r8/sqti(i)
! Subprogram not used             et2      = et*et
! Subprogram not used             et4      = et2*et2
! Subprogram not used             omet     = (1._r8 - 1.5_r8*et2)
! Subprogram not used             f1co2    = 899.70_r8*omet*(1._r8 + 1.94774_r8*et + 4.73486_r8*et2)*rsqti
! Subprogram not used             f1sqwp(i)= f1co2*sqwp
! Subprogram not used             t1co2(i) = 1._r8/(1._r8 + (245.18_r8*omet*sqwp*rsqti))
! Subprogram not used             oneme    = 1._r8 - et2
! Subprogram not used             alphat   = oneme**3*rsqti
! Subprogram not used             pi       = abs(dpnm(i))*winpl(i,kn)
! Subprogram not used             wco2     = 2.5221_r8*co2vmr(i)*pi*rga
! Subprogram not used             u7(i)    = 4.9411e4_r8*alphat*et2*wco2
! Subprogram not used             u8       = 3.9744e4_r8*alphat*et4*wco2
! Subprogram not used             u9       = 1.0447e5_r8*alphat*et4*et2*wco2
! Subprogram not used             u13      = 2.8388e3_r8*alphat*et4*wco2
! Subprogram not used             tpath    = tbar(i,kn)
! Subprogram not used             tlocal   = tbar(i,kn)
! Subprogram not used             tcrfac   = sqrt((tlocal*r250)*(tpath*r300))
! Subprogram not used             posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
! Subprogram not used             rbeta7(i)= 1._r8/(5.3228_r8*posqt)
! Subprogram not used             rbeta8   = 1._r8/(10.6576_r8*posqt)
! Subprogram not used             rbeta9   = rbeta7(i)
! Subprogram not used             rbeta13  = rbeta9
! Subprogram not used             f2co2(i) = u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i))) + &
! Subprogram not used                  u8   /sqrt(4._r8 + u8*(1._r8 + rbeta8)) + &
! Subprogram not used                  u9   /sqrt(4._r8 + u9*(1._r8 + rbeta9))
! Subprogram not used             f3co2(i) = u13/sqrt(4._r8 + u13*(1._r8 + rbeta13))
! Subprogram not used             tmp1     = log(1._r8 + f1sqwp(i))
! Subprogram not used             tmp2     = log(1._r8 + f2co2(i))
! Subprogram not used             tmp3     = log(1._r8 + f3co2(i))
! Subprogram not used             absbnd   = (tmp1 + 2._r8*t1co2(i)*tmp2 + 2._r8*tmp3)*sqti(i)
! Subprogram not used             abso(i,4)= trab2(i)*emm(i,kn)*absbnd
! Subprogram not used             tco2(i)  = 1.0_r8/(1.0_r8+ 10.0_r8*u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i))))
! Subprogram not used          end do ! do i =
! Subprogram not used !
! Subprogram not used ! Calculate trace gas absorptivity for nearest layer, abstrc
! Subprogram not used !
! Subprogram not used          call trcabn(ncol      ,                            &
! Subprogram not used               k2      ,kn      ,ucfc11  ,ucfc12  ,un2o0   , &
! Subprogram not used               un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
! Subprogram not used               uco221  ,uco222  ,uco223  ,tbar    ,bplnk   , &
! Subprogram not used               winpl   ,pinpl   ,tco2    ,th2o    ,to3     , &
! Subprogram not used               uptype  ,dw      ,s2c     ,u       ,pnew    , &
! Subprogram not used               abstrc  ,uinpl   , &
! Subprogram not used               aer_trn_ngh)
! Subprogram not used !
! Subprogram not used ! Total next layer absorptivity:
! Subprogram not used !
! Subprogram not used          do i=1,ncol
! Subprogram not used             absnxt(i,k2,kn) = abso(i,1) + abso(i,2) + &
! Subprogram not used                  abso(i,3) + abso(i,4) + abstrc(i)
! Subprogram not used          end do
! Subprogram not used       end do ! do kn =
! Subprogram not used    end do ! do k2 =
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used end subroutine radabs

!====================================================================================

! Subprogram not used subroutine radems(lchnk   ,ncol    ,                            &
! Subprogram not used                   s2c     ,tcg     ,w       ,tplnke  ,plh2o   , &
! Subprogram not used                   pnm     ,plco2   ,tint    ,tint4   ,tlayr   , &
! Subprogram not used                   tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , &
! Subprogram not used                   un2o0   ,un2o1   ,uch4    ,uco211 ,uco212   , &
! Subprogram not used                   uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
! Subprogram not used                   bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  , &
! Subprogram not used                   co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  , &
! Subprogram not used                   plh2ob  ,wb      , &
! Subprogram not used                   aer_trn_ttl, co2mmr)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Compute emissivity for H2O, CO2, O3, CH4, N2O, CFC11 and CFC12
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! H2O  ....  Uses nonisothermal emissivity method for water vapor from
! Subprogram not used !            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
! Subprogram not used !            Emissivity and Absorptivity Formulation for Water Vapor
! Subprogram not used !            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
! Subprogram not used !
! Subprogram not used !            Implementation updated by Collins,Hackney, and Edwards 2001
! Subprogram not used !               using line-by-line calculations based upon Hitran 1996 and
! Subprogram not used !               CKD 2.1 for absorptivity and emissivity
! Subprogram not used !
! Subprogram not used !            Implementation updated by Collins, Lee-Taylor, and Edwards (2003)
! Subprogram not used !               using line-by-line calculations based upon Hitran 2000 and
! Subprogram not used !               CKD 2.4 for absorptivity and emissivity
! Subprogram not used !
! Subprogram not used ! CO2  ....  Uses absorptance parameterization of the 15 micro-meter
! Subprogram not used !            (500 - 800 cm-1) band system of Carbon Dioxide, from
! Subprogram not used !            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
! Subprogram not used !            of the Absorptance Due to the 15 micro-meter Band System
! Subprogram not used !            of Carbon Dioxide Jouranl of Geophysical Research,
! Subprogram not used !            vol. 96., D5, pp 9013-9019. Also includes the effects
! Subprogram not used !            of the 9.4 and 10.4 micron bands of CO2.
! Subprogram not used !
! Subprogram not used ! O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
! Subprogram not used !            band system of ozone, from Ramanathan, V. and R. Dickinson,
! Subprogram not used !            1979: The Role of stratospheric ozone in the zonal and
! Subprogram not used !            seasonal radiative energy balance of the earth-troposphere
! Subprogram not used !            system. Journal of the Atmospheric Sciences, Vol. 36,
! Subprogram not used !            pp 1084-1104
! Subprogram not used !
! Subprogram not used ! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
! Subprogram not used !
! Subprogram not used ! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
! Subprogram not used !            bands of nitrous oxide
! Subprogram not used !
! Subprogram not used ! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
! Subprogram not used !            micron bands of CFC11
! Subprogram not used !
! Subprogram not used ! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
! Subprogram not used !            micron bands of CFC12
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Computes individual emissivities, accounting for band overlap, and
! Subprogram not used ! sums to obtain the total.
! Subprogram not used !
! Subprogram not used ! Author: W. Collins (H2O emissivity) and J. Kiehl
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: lchnk                    ! chunk identifier
! Subprogram not used    integer, intent(in) :: ncol                     ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: s2c(pcols,pverp)        ! H2o continuum path length
! Subprogram not used    real(r8), intent(in) :: tcg(pcols,pverp)        ! H2o-mass-wgted temp. (Curtis-Godson approx.)
! Subprogram not used    real(r8), intent(in) :: w(pcols,pverp)          ! H2o path length
! Subprogram not used    real(r8), intent(in) :: tplnke(pcols)           ! Layer planck temperature
! Subprogram not used    real(r8), intent(in) :: plh2o(pcols,pverp)      ! H2o prs wghted path length
! Subprogram not used    real(r8), intent(in) :: pnm(pcols,pverp)        ! Model interface pressure
! Subprogram not used    real(r8), intent(in) :: plco2(pcols,pverp)      ! Prs wghted path of co2
! Subprogram not used    real(r8), intent(in) :: tint(pcols,pverp)       ! Model interface temperatures
! Subprogram not used    real(r8), intent(in) :: tint4(pcols,pverp)      ! Tint to the 4th power
! Subprogram not used    real(r8), intent(in) :: tlayr(pcols,pverp)      ! K-1 model layer temperature
! Subprogram not used    real(r8), intent(in) :: tlayr4(pcols,pverp)     ! Tlayr to the 4th power
! Subprogram not used    real(r8), intent(in) :: plol(pcols,pverp)       ! Pressure wghtd ozone path
! Subprogram not used    real(r8), intent(in) :: plos(pcols,pverp)       ! Ozone path
! Subprogram not used    real(r8), intent(in) :: plh2ob(nbands,pcols,pverp) ! Pressure weighted h2o path with 
! Subprogram not used                                                       !    Hulst-Curtis-Godson temp. factor 
! Subprogram not used                                                       !    for H2O bands 
! Subprogram not used    real(r8), intent(in) :: wb(nbands,pcols,pverp)     ! H2o path length with 
! Subprogram not used                                                       !    Hulst-Curtis-Godson temp. factor 
! Subprogram not used                                                       !    for H2O bands 
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) 
! Subprogram not used !                               ! [fraction] Total strat. aerosol
! Subprogram not used !                               ! transmission between interfaces k1 and k2  
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Trace gas variables
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: co2mmr(pcols)           ! co2 column mean mass mixing ratio
! Subprogram not used    real(r8), intent(in) :: ucfc11(pcols,pverp)     ! CFC11 path length
! Subprogram not used    real(r8), intent(in) :: ucfc12(pcols,pverp)     ! CFC12 path length
! Subprogram not used    real(r8), intent(in) :: un2o0(pcols,pverp)      ! N2O path length
! Subprogram not used    real(r8), intent(in) :: un2o1(pcols,pverp)      ! N2O path length (hot band)
! Subprogram not used    real(r8), intent(in) :: uch4(pcols,pverp)       ! CH4 path length
! Subprogram not used    real(r8), intent(in) :: uco211(pcols,pverp)     ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco212(pcols,pverp)     ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco213(pcols,pverp)     ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco221(pcols,pverp)     ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco222(pcols,pverp)     ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco223(pcols,pverp)     ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
! Subprogram not used    real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
! Subprogram not used    real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
! Subprogram not used    real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type continuum path length
! Subprogram not used !
! Subprogram not used ! Output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: emstot(pcols,pverp)     ! Total emissivity
! Subprogram not used    real(r8), intent(out) :: co2em(pcols,pverp)      ! Layer co2 normalzd plnck funct drvtv
! Subprogram not used    real(r8), intent(out) :: co2eml(pcols,pver)      ! Intrfc co2 normalzd plnck func drvtv
! Subprogram not used    real(r8), intent(out) :: co2t(pcols,pverp)       ! Tmp and prs weighted path length
! Subprogram not used    real(r8), intent(out) :: h2otr(pcols,pverp)      ! H2o transmission over o3 band
! Subprogram not used    real(r8), intent(out) :: abplnk1(14,pcols,pverp) ! non-nearest layer Plack factor
! Subprogram not used    real(r8), intent(out) :: abplnk2(14,pcols,pverp) ! nearest layer factor
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used !
! Subprogram not used    integer i                    ! Longitude index
! Subprogram not used    integer k                    ! Level index]
! Subprogram not used    integer k1                   ! Level index
! Subprogram not used !
! Subprogram not used ! Local variables for H2O:
! Subprogram not used !
! Subprogram not used    real(r8) h2oems(pcols,pverp)     ! H2o emissivity
! Subprogram not used    real(r8) tpathe                  ! Used to compute h2o emissivity
! Subprogram not used    real(r8) dtx(pcols)              ! Planck temperature minus 250 K
! Subprogram not used    real(r8) dty(pcols)              ! Path temperature minus 250 K
! Subprogram not used !
! Subprogram not used ! The 500-800 cm^-1 emission in emis(i,4) has been combined
! Subprogram not used !              into the 0-800 cm^-1 emission in emis(i,1)
! Subprogram not used !
! Subprogram not used    real(r8) emis(pcols,2)           ! H2O emissivity 
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !
! Subprogram not used    real(r8) term7(pcols,2)          ! Kl_inf(i) in eq(r8) of table A3a of R&D
! Subprogram not used    real(r8) term8(pcols,2)          ! Delta kl_inf(i) in eq(r8)
! Subprogram not used    real(r8) tr1(pcols)              ! Equation(6) in table A2 for 650-800
! Subprogram not used    real(r8) tr2(pcols)              ! Equation(6) in table A2 for 500-650
! Subprogram not used    real(r8) tr3(pcols)              ! Equation(4) in table A2 for 650-800
! Subprogram not used    real(r8) tr4(pcols)              ! Equation(4),table A2 of R&D for 500-650
! Subprogram not used    real(r8) tr7(pcols)              ! Equation (6) times eq(4) in table A2
! Subprogram not used !                                      of R&D for 650-800 cm-1 region
! Subprogram not used    real(r8) tr8(pcols)              ! Equation (6) times eq(4) in table A2
! Subprogram not used !                                      of R&D for 500-650 cm-1 region
! Subprogram not used    real(r8) k21(pcols)              ! Exponential coefficient used to calc
! Subprogram not used !                                     rot band transmissivity in the 650-800
! Subprogram not used !                                     cm-1 region (tr1)
! Subprogram not used    real(r8) k22(pcols)              ! Exponential coefficient used to calc
! Subprogram not used !                                     rot band transmissivity in the 500-650
! Subprogram not used !                                     cm-1 region (tr2)
! Subprogram not used    real(r8) u(pcols)                ! Pressure weighted H2O path length
! Subprogram not used    real(r8) ub(nbands)              ! Pressure weighted H2O path length with
! Subprogram not used                                     !  Hulst-Curtis-Godson correction for
! Subprogram not used                                     !  each band
! Subprogram not used    real(r8) pnew                    ! Effective pressure for h2o linewidth
! Subprogram not used    real(r8) pnewb(nbands)           ! Effective pressure for h2o linewidth w/
! Subprogram not used                                     !  Hulst-Curtis-Godson correction for
! Subprogram not used                                     !  each band
! Subprogram not used    real(r8) uc1(pcols)              ! H2o continuum pathlength 500-800 cm-1
! Subprogram not used    real(r8) fwk                     ! Equation(33) in R&D far wing correction
! Subprogram not used    real(r8) troco2(pcols,pverp)     ! H2o overlap factor for co2 absorption
! Subprogram not used    real(r8) emplnk(14,pcols)        ! emissivity Planck factor
! Subprogram not used    real(r8) emstrc(pcols,pverp)     ! total trace gas emissivity
! Subprogram not used !
! Subprogram not used ! Local variables for CO2:
! Subprogram not used !
! Subprogram not used    real(r8) co2vmr(pcols)            ! CO2 column mean vmr
! Subprogram not used    real(r8) rmw                      ! ratio of molecular weights (air/co2)
! Subprogram not used    real(r8) co2ems(pcols,pverp)      ! Co2 emissivity
! Subprogram not used    real(r8) co2plk(pcols)            ! Used to compute co2 emissivity
! Subprogram not used    real(r8) sum(pcols)               ! Used to calculate path temperature
! Subprogram not used    real(r8) t1i                      ! Co2 hot band temperature factor
! Subprogram not used    real(r8) sqti                     ! Sqrt of temperature
! Subprogram not used    real(r8) pi                       ! Pressure used in co2 mean line width
! Subprogram not used    real(r8) et                       ! Co2 hot band factor
! Subprogram not used    real(r8) et2                      ! Co2 hot band factor
! Subprogram not used    real(r8) et4                      ! Co2 hot band factor
! Subprogram not used    real(r8) omet                     ! Co2 stimulated emission term
! Subprogram not used    real(r8) ex                       ! Part of co2 planck function
! Subprogram not used    real(r8) f1co2                    ! Co2 weak band factor
! Subprogram not used    real(r8) f2co2                    ! Co2 weak band factor
! Subprogram not used    real(r8) f3co2                    ! Co2 weak band factor
! Subprogram not used    real(r8) t1co2                    ! Overlap factor weak bands strong band
! Subprogram not used    real(r8) sqwp                     ! Sqrt of co2 pathlength
! Subprogram not used    real(r8) f1sqwp                   ! Main co2 band factor
! Subprogram not used    real(r8) oneme                    ! Co2 stimulated emission term
! Subprogram not used    real(r8) alphat                   ! Part of the co2 stimulated emiss term
! Subprogram not used    real(r8) wco2                     ! Consts used to define co2 pathlength
! Subprogram not used    real(r8) posqt                    ! Effective pressure for co2 line width
! Subprogram not used    real(r8) rbeta7                   ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) rbeta8                   ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) rbeta9                   ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) rbeta13                  ! Inverse of co2 hot band line width par
! Subprogram not used    real(r8) tpath                    ! Path temp used in co2 band model
! Subprogram not used    real(r8) tmp1                     ! Co2 band factor
! Subprogram not used    real(r8) tmp2                     ! Co2 band factor
! Subprogram not used    real(r8) tmp3                     ! Co2 band factor
! Subprogram not used    real(r8) tlayr5                   ! Temperature factor in co2 Planck func
! Subprogram not used    real(r8) rsqti                    ! Reciprocal of sqrt of temperature
! Subprogram not used    real(r8) exm1sq                   ! Part of co2 Planck function
! Subprogram not used    real(r8) u7                       ! Absorber amt for various co2 band systems
! Subprogram not used    real(r8) u8                       ! Absorber amt for various co2 band systems
! Subprogram not used    real(r8) u9                       ! Absorber amt for various co2 band systems
! Subprogram not used    real(r8) u13                      ! Absorber amt for various co2 band systems
! Subprogram not used    real(r8) r250                     ! Inverse 250K
! Subprogram not used    real(r8) r300                     ! Inverse 300K
! Subprogram not used    real(r8) rsslp                    ! Inverse standard sea-level pressure
! Subprogram not used !
! Subprogram not used ! Local variables for O3:
! Subprogram not used !
! Subprogram not used    real(r8) o3ems(pcols,pverp)       ! Ozone emissivity
! Subprogram not used    real(r8) dbvtt(pcols)             ! Tmp drvtv of planck fctn for tplnke
! Subprogram not used    real(r8) dbvt,fo3,t,ux,vx
! Subprogram not used    real(r8) te                       ! Temperature factor
! Subprogram not used    real(r8) u1                       ! Path length factor
! Subprogram not used    real(r8) u2                       ! Path length factor
! Subprogram not used    real(r8) phat                     ! Effecitive path length pressure
! Subprogram not used    real(r8) tlocal                   ! Local planck function temperature
! Subprogram not used    real(r8) tcrfac                   ! Scaled temperature factor
! Subprogram not used    real(r8) beta                     ! Absorption funct factor voigt effect
! Subprogram not used    real(r8) realnu                   ! Absorption function factor
! Subprogram not used    real(r8) o3bndi                   ! Band absorption factor
! Subprogram not used !
! Subprogram not used ! Transmission terms for various spectral intervals:
! Subprogram not used !
! Subprogram not used    real(r8) absbnd                   ! Proportional to co2 band absorptance
! Subprogram not used    real(r8) tco2(pcols)              ! co2 overlap factor
! Subprogram not used    real(r8) th2o(pcols)              ! h2o overlap factor
! Subprogram not used    real(r8) to3(pcols)               ! o3 overlap factor
! Subprogram not used !
! Subprogram not used ! Variables for new H2O parameterization
! Subprogram not used !
! Subprogram not used ! Notation:
! Subprogram not used ! U   = integral (P/P_0 dW)  eq. 15 in Ramanathan/Downey 1986
! Subprogram not used ! P   = atmospheric pressure
! Subprogram not used ! P_0 = reference atmospheric pressure
! Subprogram not used ! W   = precipitable water path
! Subprogram not used ! T_e = emission temperature
! Subprogram not used ! T_p = path temperature
! Subprogram not used ! RH  = path relative humidity
! Subprogram not used !
! Subprogram not used    real(r8) fe               ! asymptotic value of emis. as U->infinity
! Subprogram not used    real(r8) e_star           ! normalized non-window emissivity
! Subprogram not used    real(r8) l_star           ! interpolated line transmission
! Subprogram not used    real(r8) c_star           ! interpolated continuum transmission
! Subprogram not used 
! Subprogram not used    real(r8) te1              ! emission temperature
! Subprogram not used    real(r8) te2              ! te^2
! Subprogram not used    real(r8) te3              ! te^3
! Subprogram not used    real(r8) te4              ! te^4
! Subprogram not used    real(r8) te5              ! te^5
! Subprogram not used 
! Subprogram not used    real(r8) log_u            ! log base 10 of U 
! Subprogram not used    real(r8) log_uc           ! log base 10 of H2O continuum path
! Subprogram not used    real(r8) log_p            ! log base 10 of P
! Subprogram not used    real(r8) t_p              ! T_p
! Subprogram not used    real(r8) t_e              ! T_e (offset by T_p)
! Subprogram not used 
! Subprogram not used    integer iu                ! index for log10(U)
! Subprogram not used    integer iu1               ! iu + 1
! Subprogram not used    integer iuc               ! index for log10(H2O continuum path)
! Subprogram not used    integer iuc1              ! iuc + 1
! Subprogram not used    integer ip                ! index for log10(P)
! Subprogram not used    integer ip1               ! ip + 1
! Subprogram not used    integer itp               ! index for T_p
! Subprogram not used    integer itp1              ! itp + 1
! Subprogram not used    integer ite               ! index for T_e
! Subprogram not used    integer ite1              ! ite + 1
! Subprogram not used    integer irh               ! index for RH
! Subprogram not used    integer irh1              ! irh + 1
! Subprogram not used 
! Subprogram not used    real(r8) dvar             ! normalized variation in T_p/T_e/P/U
! Subprogram not used    real(r8) uvar             ! U * diffusivity factor
! Subprogram not used    real(r8) uscl             ! factor for lineary scaling as U->0
! Subprogram not used 
! Subprogram not used    real(r8) wu               ! weight for U
! Subprogram not used    real(r8) wu1              ! 1 - wu
! Subprogram not used    real(r8) wuc              ! weight for H2O continuum path
! Subprogram not used    real(r8) wuc1             ! 1 - wuc
! Subprogram not used    real(r8) wp               ! weight for P
! Subprogram not used    real(r8) wp1              ! 1 - wp
! Subprogram not used    real(r8) wtp              ! weight for T_p
! Subprogram not used    real(r8) wtp1             ! 1 - wtp
! Subprogram not used    real(r8) wte              ! weight for T_e
! Subprogram not used    real(r8) wte1             ! 1 - wte
! Subprogram not used    real(r8) wrh              ! weight for RH
! Subprogram not used    real(r8) wrh1             ! 1 - wrh
! Subprogram not used 
! Subprogram not used    real(r8) w_0_0_           ! weight for Tp/Te combination
! Subprogram not used    real(r8) w_0_1_           ! weight for Tp/Te combination
! Subprogram not used    real(r8) w_1_0_           ! weight for Tp/Te combination
! Subprogram not used    real(r8) w_1_1_           ! weight for Tp/Te combination
! Subprogram not used 
! Subprogram not used    real(r8) w_0_00           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_0_01           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_0_10           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_0_11           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_00           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_01           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_10           ! weight for Tp/Te/RH combination
! Subprogram not used    real(r8) w_1_11           ! weight for Tp/Te/RH combination
! Subprogram not used 
! Subprogram not used    real(r8) w00_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w00_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w00_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w00_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w01_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w10_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_00           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_01           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_10           ! weight for P/Tp/Te/RH combination
! Subprogram not used    real(r8) w11_11           ! weight for P/Tp/Te/RH combination
! Subprogram not used 
! Subprogram not used    integer ib                ! spectral interval:
! Subprogram not used                              !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
! Subprogram not used                              !   2 = 800-1200 cm^-1
! Subprogram not used 
! Subprogram not used    real(r8) pch2o            ! H2O continuum path
! Subprogram not used    real(r8) fch2o            ! temp. factor for continuum
! Subprogram not used    real(r8) uch2o            ! U corresponding to H2O cont. path (window)
! Subprogram not used 
! Subprogram not used    real(r8) fdif             ! secant(zenith angle) for diffusivity approx.
! Subprogram not used 
! Subprogram not used    real(r8) sslp_mks         ! Sea-level pressure in MKS units
! Subprogram not used    real(r8) esx              ! saturation vapor pressure returned by qsat
! Subprogram not used    real(r8) qsx              ! saturation mixing ratio returned by qsat
! Subprogram not used    real(r8) pnew_mks         ! pnew in MKS units
! Subprogram not used    real(r8) q_path           ! effective specific humidity along path
! Subprogram not used    real(r8) rh_path          ! effective relative humidity along path
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !---------------------------Statement functions-------------------------
! Subprogram not used !
! Subprogram not used ! Derivative of planck function at 9.6 micro-meter wavelength, and
! Subprogram not used ! an absorption function factor:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used    dbvt(t)=(-2.8911366682e-4_r8+(2.3771251896e-6_r8+1.1305188929e-10_r8*t)*t)/ &
! Subprogram not used            (1.0_r8+(-6.1364820707e-3_r8+1.5550319767e-5_r8*t)*t)
! Subprogram not used !
! Subprogram not used    fo3(ux,vx)=ux/sqrt(4._r8+ux*(1._r8+vx))
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Initialize
! Subprogram not used !
! Subprogram not used    r250  = 1._r8/250._r8
! Subprogram not used    r300  = 1._r8/300._r8
! Subprogram not used    rsslp = 1._r8/sslp
! Subprogram not used    rmw   = amd/amco2
! Subprogram not used    do i=1,ncol
! Subprogram not used       co2vmr(i) = co2mmr(i) * rmw
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Constants for computing U corresponding to H2O cont. path
! Subprogram not used !
! Subprogram not used    fdif       = 1.66_r8
! Subprogram not used    sslp_mks   = sslp / 10.0_r8
! Subprogram not used !
! Subprogram not used ! Planck function for co2
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       ex             = exp(960._r8/tplnke(i))
! Subprogram not used       co2plk(i)      = 5.e8_r8/((tplnke(i)**4)*(ex - 1._r8))
! Subprogram not used       co2t(i,ntoplw) = tplnke(i)
! Subprogram not used       sum(i)         = co2t(i,ntoplw)*pnm(i,ntoplw)
! Subprogram not used    end do
! Subprogram not used    k = ntoplw
! Subprogram not used    do k1=pverp,ntoplw+1,-1
! Subprogram not used       k = k + 1
! Subprogram not used       do i=1,ncol
! Subprogram not used          sum(i)         = sum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
! Subprogram not used          ex             = exp(960._r8/tlayr(i,k1))
! Subprogram not used          tlayr5         = tlayr(i,k1)*tlayr4(i,k1)
! Subprogram not used          co2eml(i,k1-1) = 1.2e11_r8*ex/(tlayr5*(ex - 1._r8)**2)
! Subprogram not used          co2t(i,k)      = sum(i)/pnm(i,k)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Initialize planck function derivative for O3
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       dbvtt(i) = dbvt(tplnke(i))
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Calculate trace gas Planck functions
! Subprogram not used !
! Subprogram not used    call trcplk(ncol    ,                                     &
! Subprogram not used                tint    ,tlayr   ,tplnke  ,emplnk  ,abplnk1 , &
! Subprogram not used                abplnk2 )
! Subprogram not used          
! Subprogram not used    if ( ntoplw > 1 )then
! Subprogram not used       emstot(:ncol,:ntoplw-1) = 0._r8
! Subprogram not used    end if
! Subprogram not used           
! Subprogram not used !
! Subprogram not used ! Interface loop
! Subprogram not used !
! Subprogram not used    do k1=ntoplw,pverp
! Subprogram not used !
! Subprogram not used ! H2O emissivity
! Subprogram not used !
! Subprogram not used ! emis(i,1)     0 -  800 cm-1   h2o rotation band
! Subprogram not used ! emis(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! Subprogram not used ! emis(i,2)   800 - 1200 cm-1   h2o window
! Subprogram not used !
! Subprogram not used ! Separation between rotation and vibration-rotation dropped, so
! Subprogram not used !                only 2 slots needed for H2O emissivity
! Subprogram not used !
! Subprogram not used !      emis(i,3)   = 0.0
! Subprogram not used !
! Subprogram not used ! For the p type continuum
! Subprogram not used !
! Subprogram not used       do i=1,ncol
! Subprogram not used          u(i)        = plh2o(i,k1)
! Subprogram not used          pnew        = u(i)/w(i,k1)
! Subprogram not used          pnew_mks    = pnew * sslp_mks
! Subprogram not used !
! Subprogram not used ! Apply scaling factor for 500-800 continuum
! Subprogram not used !
! Subprogram not used          uc1(i)      = (s2c(i,k1) + 1.7e-3_r8*plh2o(i,k1))*(1._r8 + 2._r8*s2c(i,k1))/ &
! Subprogram not used                        (1._r8 + 15._r8*s2c(i,k1))
! Subprogram not used          pch2o       = s2c(i,k1)
! Subprogram not used !
! Subprogram not used ! Changed effective path temperature to std. Curtis-Godson form
! Subprogram not used !
! Subprogram not used          tpathe   = tcg(i,k1)/w(i,k1)
! Subprogram not used          t_p = min(max(tpathe, min_tp_h2o), max_tp_h2o)
! Subprogram not used 
! Subprogram not used          call qsat_water(t_p, pnew_mks, esx, qsx)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Compute effective RH along path
! Subprogram not used !
! Subprogram not used          q_path = w(i,k1) / pnm(i,k1) / rga
! Subprogram not used !
! Subprogram not used ! Calculate effective u, pnew for each band using
! Subprogram not used !        Hulst-Curtis-Godson approximation:
! Subprogram not used ! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
! Subprogram not used !           2nd edition, Oxford University Press, 1989.
! Subprogram not used ! Effective H2O path (w)
! Subprogram not used !      eq. 6.24, p. 228
! Subprogram not used ! Effective H2O path pressure (pnew = u/w):
! Subprogram not used !      eq. 6.29, p. 228
! Subprogram not used !
! Subprogram not used          ub(1) = plh2ob(1,i,k1) / psi(t_p,1)
! Subprogram not used          ub(2) = plh2ob(2,i,k1) / psi(t_p,2)
! Subprogram not used 
! Subprogram not used          pnewb(1) = ub(1) / wb(1,i,k1) * phi(t_p,1)
! Subprogram not used          pnewb(2) = ub(2) / wb(2,i,k1) * phi(t_p,2)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !
! Subprogram not used          dtx(i) = tplnke(i) - 250._r8
! Subprogram not used          dty(i) = tpathe - 250._r8
! Subprogram not used !
! Subprogram not used ! Define variables for C/H/E (now C/LT/E) fit
! Subprogram not used !
! Subprogram not used ! emis(i,1)     0 -  800 cm-1   h2o rotation band
! Subprogram not used ! emis(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! Subprogram not used ! emis(i,2)   800 - 1200 cm-1   h2o window
! Subprogram not used !
! Subprogram not used ! Separation between rotation and vibration-rotation dropped, so
! Subprogram not used !                only 2 slots needed for H2O emissivity
! Subprogram not used !
! Subprogram not used ! emis(i,3)   = 0.0
! Subprogram not used !
! Subprogram not used ! Notation:
! Subprogram not used ! U   = integral (P/P_0 dW)  
! Subprogram not used ! P   = atmospheric pressure
! Subprogram not used ! P_0 = reference atmospheric pressure
! Subprogram not used ! W   = precipitable water path
! Subprogram not used ! T_e = emission temperature
! Subprogram not used ! T_p = path temperature
! Subprogram not used ! RH  = path relative humidity
! Subprogram not used !
! Subprogram not used ! Terms for asymptotic value of emissivity
! Subprogram not used !
! Subprogram not used          te1  = tplnke(i)
! Subprogram not used          te2  = te1 * te1
! Subprogram not used          te3  = te2 * te1
! Subprogram not used          te4  = te3 * te1
! Subprogram not used          te5  = te4 * te1
! Subprogram not used !
! Subprogram not used ! Band-independent indices for lines and continuum tables
! Subprogram not used !
! Subprogram not used          dvar = (t_p - min_tp_h2o) / dtp_h2o
! Subprogram not used          itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
! Subprogram not used          itp1 = itp + 1
! Subprogram not used          wtp = dvar - floor(dvar)
! Subprogram not used          wtp1 = 1.0_r8 - wtp
! Subprogram not used 
! Subprogram not used          t_e = min(max(tplnke(i) - t_p, min_te_h2o), max_te_h2o)
! Subprogram not used          dvar = (t_e - min_te_h2o) / dte_h2o
! Subprogram not used          ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
! Subprogram not used          ite1 = ite + 1
! Subprogram not used          wte = dvar - floor(dvar)
! Subprogram not used          wte1 = 1.0_r8 - wte
! Subprogram not used 
! Subprogram not used          rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
! Subprogram not used          dvar = (rh_path - min_rh_h2o) / drh_h2o
! Subprogram not used          irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
! Subprogram not used          irh1 = irh + 1
! Subprogram not used          wrh = dvar - floor(dvar)
! Subprogram not used          wrh1 = 1.0_r8 - wrh
! Subprogram not used 
! Subprogram not used          w_0_0_ = wtp  * wte
! Subprogram not used          w_0_1_ = wtp  * wte1
! Subprogram not used          w_1_0_ = wtp1 * wte 
! Subprogram not used          w_1_1_ = wtp1 * wte1
! Subprogram not used 
! Subprogram not used          w_0_00 = w_0_0_ * wrh
! Subprogram not used          w_0_01 = w_0_0_ * wrh1
! Subprogram not used          w_0_10 = w_0_1_ * wrh
! Subprogram not used          w_0_11 = w_0_1_ * wrh1
! Subprogram not used          w_1_00 = w_1_0_ * wrh
! Subprogram not used          w_1_01 = w_1_0_ * wrh1
! Subprogram not used          w_1_10 = w_1_1_ * wrh
! Subprogram not used          w_1_11 = w_1_1_ * wrh1
! Subprogram not used !
! Subprogram not used ! H2O Continuum path for 0-800 and 1200-2200 cm^-1
! Subprogram not used !
! Subprogram not used !    Assume foreign continuum dominates total H2O continuum in these bands
! Subprogram not used !    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
! Subprogram not used !    Then the effective H2O path is just 
! Subprogram not used !         U_c = integral[ f(P) dW ]
! Subprogram not used !    where 
! Subprogram not used !           W = water-vapor mass and 
! Subprogram not used !        f(P) = dependence of foreign continuum on pressure 
! Subprogram not used !             = P / sslp
! Subprogram not used !    Then 
! Subprogram not used !         U_c = U (the same effective H2O path as for lines)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! Continuum terms for 800-1200 cm^-1
! Subprogram not used !
! Subprogram not used !    Assume self continuum dominates total H2O continuum for this band
! Subprogram not used !    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
! Subprogram not used !    Then the effective H2O self-continuum path is 
! Subprogram not used !         U_c = integral[ h(e,T) dW ]                        (*eq. 1*)
! Subprogram not used !    where 
! Subprogram not used !           W = water-vapor mass and 
! Subprogram not used !           e = partial pressure of H2O along path
! Subprogram not used !           T = temperature along path
! Subprogram not used !      h(e,T) = dependence of foreign continuum on e,T
! Subprogram not used !             = e / sslp * f(T)
! Subprogram not used !
! Subprogram not used !    Replacing
! Subprogram not used !           e =~ q * P / epsilo
! Subprogram not used !           q = mixing ratio of H2O
! Subprogram not used !     epsilo = 0.622
! Subprogram not used !
! Subprogram not used !    and using the definition
! Subprogram not used !           U = integral [ (P / sslp) dW ]
! Subprogram not used !             = (P / sslp) W                                 (homogeneous path)
! Subprogram not used !
! Subprogram not used !    the effective path length for the self continuum is
! Subprogram not used !         U_c = (q / epsilo) f(T) U                         (*eq. 2*)
! Subprogram not used !
! Subprogram not used !    Once values of T, U, and q have been calculated for the inhomogeneous
! Subprogram not used !        path, this sets U_c for the corresponding
! Subprogram not used !        homogeneous atmosphere.  However, this need not equal the
! Subprogram not used !        value of U_c' defined by eq. 1 for the actual inhomogeneous atmosphere
! Subprogram not used !        under consideration.
! Subprogram not used !
! Subprogram not used !    Solution: hold T and q constant, solve for U' that gives U_c' by
! Subprogram not used !        inverting eq. (2):
! Subprogram not used !
! Subprogram not used !        U' = (U_c * epsilo) / (q * f(T))
! Subprogram not used !
! Subprogram not used          fch2o = fh2oself(t_p)
! Subprogram not used          uch2o = (pch2o * epsilo) / (q_path * fch2o)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Band-dependent indices for non-window
! Subprogram not used !
! Subprogram not used          ib = 1
! Subprogram not used 
! Subprogram not used          uvar = ub(ib) * fdif
! Subprogram not used          log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
! Subprogram not used          dvar = (log_u - min_lu_h2o) / dlu_h2o
! Subprogram not used          iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used          iu1 = iu + 1
! Subprogram not used          wu = dvar - floor(dvar)
! Subprogram not used          wu1 = 1.0_r8 - wu
! Subprogram not used          
! Subprogram not used          log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
! Subprogram not used          dvar = (log_p - min_lp_h2o) / dlp_h2o
! Subprogram not used          ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
! Subprogram not used          ip1 = ip + 1
! Subprogram not used          wp = dvar - floor(dvar)
! Subprogram not used          wp1 = 1.0_r8 - wp
! Subprogram not used 
! Subprogram not used          w00_00 = wp  * w_0_00 
! Subprogram not used          w00_01 = wp  * w_0_01 
! Subprogram not used          w00_10 = wp  * w_0_10 
! Subprogram not used          w00_11 = wp  * w_0_11 
! Subprogram not used          w01_00 = wp  * w_1_00 
! Subprogram not used          w01_01 = wp  * w_1_01 
! Subprogram not used          w01_10 = wp  * w_1_10 
! Subprogram not used          w01_11 = wp  * w_1_11 
! Subprogram not used          w10_00 = wp1 * w_0_00 
! Subprogram not used          w10_01 = wp1 * w_0_01 
! Subprogram not used          w10_10 = wp1 * w_0_10 
! Subprogram not used          w10_11 = wp1 * w_0_11 
! Subprogram not used          w11_00 = wp1 * w_1_00 
! Subprogram not used          w11_01 = wp1 * w_1_01 
! Subprogram not used          w11_10 = wp1 * w_1_10 
! Subprogram not used          w11_11 = wp1 * w_1_11 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Asymptotic value of emissivity as U->infinity
! Subprogram not used !
! Subprogram not used          fe = fet(1,ib) + &
! Subprogram not used               fet(2,ib) * te1 + &
! Subprogram not used               fet(3,ib) * te2 + &
! Subprogram not used               fet(4,ib) * te3 + &
! Subprogram not used               fet(5,ib) * te4 + &
! Subprogram not used               fet(6,ib) * te5
! Subprogram not used 
! Subprogram not used          e_star = &
! Subprogram not used               eh2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
! Subprogram not used               eh2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
! Subprogram not used               eh2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
! Subprogram not used               eh2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
! Subprogram not used               eh2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
! Subprogram not used               eh2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
! Subprogram not used               eh2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
! Subprogram not used               eh2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
! Subprogram not used               eh2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
! Subprogram not used               eh2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
! Subprogram not used               eh2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
! Subprogram not used               eh2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
! Subprogram not used               eh2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
! Subprogram not used               eh2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
! Subprogram not used               eh2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
! Subprogram not used               eh2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
! Subprogram not used               eh2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
! Subprogram not used               eh2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
! Subprogram not used               eh2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
! Subprogram not used               eh2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
! Subprogram not used               eh2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
! Subprogram not used               eh2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
! Subprogram not used               eh2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
! Subprogram not used               eh2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
! Subprogram not used               eh2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
! Subprogram not used          emis(i,ib) = min(max(fe * (1.0_r8 - (1.0_r8 - e_star) * &
! Subprogram not used                               aer_trn_ttl(i,k1,1,ib)), &
! Subprogram not used                           0.0_r8), 1.0_r8)
! Subprogram not used !
! Subprogram not used ! Invoke linear limit for scaling wrt u below min_u_h2o
! Subprogram not used !
! Subprogram not used          if (uvar < min_u_h2o) then
! Subprogram not used             uscl = uvar / min_u_h2o
! Subprogram not used             emis(i,ib) = emis(i,ib) * uscl
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used                       
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Band-dependent indices for window
! Subprogram not used !
! Subprogram not used          ib = 2
! Subprogram not used 
! Subprogram not used          uvar = ub(ib) * fdif
! Subprogram not used          log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
! Subprogram not used          dvar = (log_u - min_lu_h2o) / dlu_h2o
! Subprogram not used          iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used          iu1 = iu + 1
! Subprogram not used          wu = dvar - floor(dvar)
! Subprogram not used          wu1 = 1.0_r8 - wu
! Subprogram not used          
! Subprogram not used          log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
! Subprogram not used          dvar = (log_p - min_lp_h2o) / dlp_h2o
! Subprogram not used          ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
! Subprogram not used          ip1 = ip + 1
! Subprogram not used          wp = dvar - floor(dvar)
! Subprogram not used          wp1 = 1.0_r8 - wp
! Subprogram not used 
! Subprogram not used          w00_00 = wp  * w_0_00 
! Subprogram not used          w00_01 = wp  * w_0_01 
! Subprogram not used          w00_10 = wp  * w_0_10 
! Subprogram not used          w00_11 = wp  * w_0_11 
! Subprogram not used          w01_00 = wp  * w_1_00 
! Subprogram not used          w01_01 = wp  * w_1_01 
! Subprogram not used          w01_10 = wp  * w_1_10 
! Subprogram not used          w01_11 = wp  * w_1_11 
! Subprogram not used          w10_00 = wp1 * w_0_00 
! Subprogram not used          w10_01 = wp1 * w_0_01 
! Subprogram not used          w10_10 = wp1 * w_0_10 
! Subprogram not used          w10_11 = wp1 * w_0_11 
! Subprogram not used          w11_00 = wp1 * w_1_00 
! Subprogram not used          w11_01 = wp1 * w_1_01 
! Subprogram not used          w11_10 = wp1 * w_1_10 
! Subprogram not used          w11_11 = wp1 * w_1_11 
! Subprogram not used 
! Subprogram not used          log_uc  = min(log10(max(uch2o * fdif, min_u_h2o)), max_lu_h2o)
! Subprogram not used          dvar = (log_uc - min_lu_h2o) / dlu_h2o
! Subprogram not used          iuc = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
! Subprogram not used          iuc1 = iuc + 1
! Subprogram not used          wuc = dvar - floor(dvar)
! Subprogram not used          wuc1 = 1.0_r8 - wuc
! Subprogram not used !
! Subprogram not used ! Asymptotic value of emissivity as U->infinity
! Subprogram not used !
! Subprogram not used          fe = fet(1,ib) + &
! Subprogram not used               fet(2,ib) * te1 + &
! Subprogram not used               fet(3,ib) * te2 + &
! Subprogram not used               fet(4,ib) * te3 + &
! Subprogram not used               fet(5,ib) * te4 + &
! Subprogram not used               fet(6,ib) * te5
! Subprogram not used 
! Subprogram not used          l_star = &
! Subprogram not used               ln_eh2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
! Subprogram not used               ln_eh2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
! Subprogram not used               ln_eh2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
! Subprogram not used 
! Subprogram not used          c_star = &
! Subprogram not used               cn_eh2ow(ip , itp , iuc , ite , irh ) * w11_11 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc , ite , irh1) * w11_10 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc , ite1, irh ) * w11_01 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc , ite1, irh1) * w11_00 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc1, ite , irh ) * w11_11 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc1, ite , irh1) * w11_10 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc1, ite1, irh ) * w11_01 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp , iuc1, ite1, irh1) * w11_00 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc , ite , irh ) * w10_11 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc , ite , irh1) * w10_10 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc , ite1, irh ) * w10_01 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc , ite1, irh1) * w10_00 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc1, ite , irh ) * w10_11 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc1, ite , irh1) * w10_10 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc1, ite1, irh ) * w10_01 * wuc  + &
! Subprogram not used               cn_eh2ow(ip , itp1, iuc1, ite1, irh1) * w10_00 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc , ite , irh ) * w01_11 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc , ite , irh1) * w01_10 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc , ite1, irh ) * w01_01 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc , ite1, irh1) * w01_00 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc1, ite , irh ) * w01_11 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc1, ite , irh1) * w01_10 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc1, ite1, irh ) * w01_01 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp , iuc1, ite1, irh1) * w01_00 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc , ite , irh ) * w00_11 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc , ite , irh1) * w00_10 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc , ite1, irh ) * w00_01 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc , ite1, irh1) * w00_00 * wuc1 + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc1, ite , irh ) * w00_11 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc1, ite , irh1) * w00_10 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc1, ite1, irh ) * w00_01 * wuc  + &
! Subprogram not used               cn_eh2ow(ip1, itp1, iuc1, ite1, irh1) * w00_00 * wuc 
! Subprogram not used          emis(i,ib) = min(max(fe * (1.0_r8 - l_star * c_star * &
! Subprogram not used                               aer_trn_ttl(i,k1,1,ib)), &
! Subprogram not used                           0.0_r8), 1.0_r8) 
! Subprogram not used !
! Subprogram not used ! Invoke linear limit for scaling wrt u below min_u_h2o
! Subprogram not used !
! Subprogram not used          if (uvar < min_u_h2o) then
! Subprogram not used             uscl = uvar / min_u_h2o
! Subprogram not used             emis(i,ib) = emis(i,ib) * uscl
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used                       
! Subprogram not used !
! Subprogram not used ! Compute total emissivity for H2O
! Subprogram not used !
! Subprogram not used          h2oems(i,k1) = emis(i,1)+emis(i,2)
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       do i=1,ncol
! Subprogram not used          term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1._r8+c16*dty(i))
! Subprogram not used          term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1._r8+c17*dty(i))
! Subprogram not used          term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1._r8+c26*dty(i))
! Subprogram not used          term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1._r8+c27*dty(i))
! Subprogram not used       end do
! Subprogram not used       do i=1,ncol
! Subprogram not used !
! Subprogram not used ! 500 -  800 cm-1   rotation band overlap with co2
! Subprogram not used !
! Subprogram not used          k21(i) = term7(i,1) + term8(i,1)/ &
! Subprogram not used                  (1._r8 + (c30 + c31*(dty(i)-10._r8)*(dty(i)-10._r8))*sqrt(u(i)))
! Subprogram not used          k22(i) = term7(i,2) + term8(i,2)/ &
! Subprogram not used                  (1._r8 + (c28 + c29*(dty(i)-10._r8))*sqrt(u(i)))
! Subprogram not used          fwk    = fwcoef + fwc1/(1._r8+fwc2*u(i))
! Subprogram not used          tr1(i) = exp(-(k21(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
! Subprogram not used          tr2(i) = exp(-(k22(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
! Subprogram not used          tr1(i)=tr1(i)*aer_trn_ttl(i,k1,1,idx_LW_0650_0800) 
! Subprogram not used !                                            ! H2O line+aer trn 650--800 cm-1
! Subprogram not used          tr2(i)=tr2(i)*aer_trn_ttl(i,k1,1,idx_LW_0500_0650) 
! Subprogram not used !                                            ! H2O line+aer trn 500--650 cm-1
! Subprogram not used          tr3(i) = exp(-((coefh(1,1) + coefh(2,1)*dtx(i))*uc1(i)))
! Subprogram not used          tr4(i) = exp(-((coefh(1,2) + coefh(2,2)*dtx(i))*uc1(i)))
! Subprogram not used          tr7(i) = tr1(i)*tr3(i)
! Subprogram not used          tr8(i) = tr2(i)*tr4(i)
! Subprogram not used          troco2(i,k1) = 0.65_r8*tr7(i) + 0.35_r8*tr8(i)
! Subprogram not used          th2o(i) = tr8(i)
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used ! CO2 emissivity for 15 micron band system
! Subprogram not used !
! Subprogram not used       do i=1,ncol
! Subprogram not used          t1i    = exp(-480._r8/co2t(i,k1))
! Subprogram not used          sqti   = sqrt(co2t(i,k1))
! Subprogram not used          rsqti  = 1._r8/sqti
! Subprogram not used          et     = t1i
! Subprogram not used          et2    = et*et
! Subprogram not used          et4    = et2*et2
! Subprogram not used          omet   = 1._r8 - 1.5_r8*et2
! Subprogram not used          f1co2  = 899.70_r8*omet*(1._r8 + 1.94774_r8*et + 4.73486_r8*et2)*rsqti
! Subprogram not used          sqwp   = sqrt(plco2(i,k1))
! Subprogram not used          f1sqwp = f1co2*sqwp
! Subprogram not used          t1co2  = 1._r8/(1._r8 + 245.18_r8*omet*sqwp*rsqti)
! Subprogram not used          oneme  = 1._r8 - et2
! Subprogram not used          alphat = oneme**3*rsqti
! Subprogram not used          wco2   = 2.5221_r8*co2vmr(i)*pnm(i,k1)*rga
! Subprogram not used          u7     = 4.9411e4_r8*alphat*et2*wco2
! Subprogram not used          u8     = 3.9744e4_r8*alphat*et4*wco2
! Subprogram not used          u9     = 1.0447e5_r8*alphat*et4*et2*wco2
! Subprogram not used          u13    = 2.8388e3_r8*alphat*et4*wco2
! Subprogram not used !
! Subprogram not used          tpath  = co2t(i,k1)
! Subprogram not used          tlocal = tplnke(i)
! Subprogram not used          tcrfac = sqrt((tlocal*r250)*(tpath*r300))
! Subprogram not used          pi     = pnm(i,k1)*rsslp + 2._r8*dpfco2*tcrfac
! Subprogram not used          posqt  = pi/(2._r8*sqti)
! Subprogram not used          rbeta7 =  1._r8/( 5.3288_r8*posqt)
! Subprogram not used          rbeta8 = 1._r8/ (10.6576_r8*posqt)
! Subprogram not used          rbeta9 = rbeta7
! Subprogram not used          rbeta13= rbeta9
! Subprogram not used          f2co2  = (u7/sqrt(4._r8 + u7*(1._r8 + rbeta7))) + &
! Subprogram not used                   (u8/sqrt(4._r8 + u8*(1._r8 + rbeta8))) + &
! Subprogram not used                   (u9/sqrt(4._r8 + u9*(1._r8 + rbeta9)))
! Subprogram not used          f3co2  = u13/sqrt(4._r8 + u13*(1._r8 + rbeta13))
! Subprogram not used          tmp1   = log(1._r8 + f1sqwp)
! Subprogram not used          tmp2   = log(1._r8 +  f2co2)
! Subprogram not used          tmp3   = log(1._r8 +  f3co2)
! Subprogram not used          absbnd = (tmp1 + 2._r8*t1co2*tmp2 + 2._r8*tmp3)*sqti
! Subprogram not used          tco2(i)=1.0_r8/(1.0_r8+10.0_r8*(u7/sqrt(4._r8 + u7*(1._r8 + rbeta7))))
! Subprogram not used          co2ems(i,k1)  = troco2(i,k1)*absbnd*co2plk(i)
! Subprogram not used          ex     = exp(960._r8/tint(i,k1))
! Subprogram not used          exm1sq = (ex - 1._r8)**2
! Subprogram not used          co2em(i,k1) = 1.2e11_r8*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used ! O3 emissivity
! Subprogram not used !
! Subprogram not used       do i=1,ncol
! Subprogram not used          h2otr(i,k1) = exp(-12._r8*s2c(i,k1))
! Subprogram not used           h2otr(i,k1)=h2otr(i,k1)*aer_trn_ttl(i,k1,1,idx_LW_1000_1200)
! Subprogram not used          te          = (co2t(i,k1)/293._r8)**.7_r8
! Subprogram not used          u1          = 18.29_r8*plos(i,k1)/te
! Subprogram not used          u2          = .5649_r8*plos(i,k1)/te
! Subprogram not used          phat        = plos(i,k1)/plol(i,k1)
! Subprogram not used          tlocal      = tplnke(i)
! Subprogram not used          tcrfac      = sqrt(tlocal*r250)*te
! Subprogram not used          beta        = (1._r8/.3205_r8)*((1._r8/phat) + (dpfo3*tcrfac))
! Subprogram not used          realnu      = (1._r8/beta)*te
! Subprogram not used          o3bndi      = 74._r8*te*(tplnke(i)/375._r8)*log(1._r8 + fo3(u1,realnu) + fo3(u2,realnu))
! Subprogram not used          o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
! Subprogram not used          to3(i)=1.0_r8/(1._r8 + 0.1_r8*fo3(u1,realnu) + 0.1_r8*fo3(u2,realnu))
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used !   Calculate trace gas emissivities
! Subprogram not used !
! Subprogram not used       call trcems(ncol    ,                                     &
! Subprogram not used                   k1      ,co2t    ,pnm     ,ucfc11  ,ucfc12  , &
! Subprogram not used                   un2o0   ,un2o1   ,bn2o0   ,bn2o1   ,uch4    , &
! Subprogram not used                   bch4    ,uco211  ,uco212  ,uco213  ,uco221  , &
! Subprogram not used                   uco222  ,uco223  ,uptype  ,w       ,s2c     , &
! Subprogram not used                   u       ,emplnk  ,th2o    ,tco2    ,to3     , &
! Subprogram not used                   emstrc  , &
! Subprogram not used                   aer_trn_ttl)
! Subprogram not used !
! Subprogram not used ! Total emissivity:
! Subprogram not used !
! Subprogram not used       do i=1,ncol
! Subprogram not used          emstot(i,k1) = h2oems(i,k1) + co2ems(i,k1) + o3ems(i,k1)  &
! Subprogram not used                         + emstrc(i,k1)
! Subprogram not used       end do
! Subprogram not used    end do ! End of interface loop
! Subprogram not used 
! Subprogram not used end subroutine radems

!====================================================================================

! Subprogram not used subroutine radtpl(ncol    ,                                     &
! Subprogram not used                   tnm     ,lwupcgs ,qnm     ,pnm     ,plco2   ,plh2o   , &
! Subprogram not used                   tplnka  ,s2c     ,tcg     ,w       ,tplnke  , &
! Subprogram not used                   tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    , &
! Subprogram not used                   piln    ,plh2ob  ,wb      ,co2mmr)
! Subprogram not used !--------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Purpose:
! Subprogram not used ! Compute temperatures and path lengths for longwave radiation
! Subprogram not used !
! Subprogram not used ! Method:
! Subprogram not used ! <Describe the algorithm(s) used in the routine.>
! Subprogram not used ! <Also include any applicable external references.>
! Subprogram not used !
! Subprogram not used ! Author: CCM1
! Subprogram not used !------------------------------Arguments-----------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: tnm(pcols,pver)      ! Model level temperatures
! Subprogram not used    real(r8), intent(in) :: lwupcgs(pcols)       ! Surface longwave up flux
! Subprogram not used    real(r8), intent(in) :: qnm(pcols,pver)      ! Model level specific humidity
! Subprogram not used    real(r8), intent(in) :: pnm(pcols,pverp)     ! Pressure at model interfaces (dynes/cm2)
! Subprogram not used    real(r8), intent(in) :: pmln(pcols,pver)     ! Ln(pmidm1)
! Subprogram not used    real(r8), intent(in) :: piln(pcols,pverp)    ! Ln(pintm1)
! Subprogram not used    real(r8), intent(in) :: co2mmr(pcols)        ! co2 column mean mass mixing ratio
! Subprogram not used !
! Subprogram not used ! Output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: plco2(pcols,pverp)   ! Pressure weighted co2 path
! Subprogram not used    real(r8), intent(out) :: plh2o(pcols,pverp)   ! Pressure weighted h2o path
! Subprogram not used    real(r8), intent(out) :: tplnka(pcols,pverp)  ! Level temperature from interface temperatures
! Subprogram not used    real(r8), intent(out) :: s2c(pcols,pverp)     ! H2o continuum path length
! Subprogram not used    real(r8), intent(out) :: tcg(pcols,pverp)     ! H2o-mass-wgted temp. (Curtis-Godson approx.)
! Subprogram not used    real(r8), intent(out) :: w(pcols,pverp)       ! H2o path length
! Subprogram not used    real(r8), intent(out) :: tplnke(pcols)        ! Equal to tplnka
! Subprogram not used    real(r8), intent(out) :: tint(pcols,pverp)    ! Layer interface temperature
! Subprogram not used    real(r8), intent(out) :: tint4(pcols,pverp)   ! Tint to the 4th power
! Subprogram not used    real(r8), intent(out) :: tlayr(pcols,pverp)   ! K-1 level temperature
! Subprogram not used    real(r8), intent(out) :: tlayr4(pcols,pverp)  ! Tlayr to the 4th power
! Subprogram not used    real(r8), intent(out) :: plh2ob(nbands,pcols,pverp)! Pressure weighted h2o path with 
! Subprogram not used                                                       !    Hulst-Curtis-Godson temp. factor 
! Subprogram not used                                                       !    for H2O bands 
! Subprogram not used    real(r8), intent(out) :: wb(nbands,pcols,pverp)    ! H2o path length with 
! Subprogram not used                                                       !    Hulst-Curtis-Godson temp. factor 
! Subprogram not used                                                       !    for H2O bands 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !---------------------------Local variables--------------------------
! Subprogram not used !
! Subprogram not used    integer i                 ! Longitude index
! Subprogram not used    integer k                 ! Level index
! Subprogram not used    integer kp1               ! Level index + 1
! Subprogram not used 
! Subprogram not used    real(r8) repsil               ! Inver ratio mol weight h2o to dry air
! Subprogram not used    real(r8) dy                   ! Thickness of layer for tmp interp
! Subprogram not used    real(r8) dpnm                 ! Pressure thickness of layer
! Subprogram not used    real(r8) dpnmsq               ! Prs squared difference across layer
! Subprogram not used    real(r8) dw                   ! Increment in H2O path length
! Subprogram not used    real(r8) dplh2o               ! Increment in plh2o
! Subprogram not used    real(r8) cpwpl                ! Const in co2 mix ratio to path length conversn
! Subprogram not used 
! Subprogram not used !--------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    repsil = 1._r8/epsilo
! Subprogram not used !
! Subprogram not used ! Compute co2 and h2o paths
! Subprogram not used !
! Subprogram not used    cpwpl = 0.5_r8/(gravit_cgs*p0)
! Subprogram not used    do i=1,ncol
! Subprogram not used       plh2o(i,ntoplw)  = rgsslp*qnm(i,ntoplw)*pnm(i,ntoplw)*pnm(i,ntoplw)
! Subprogram not used       plco2(i,ntoplw)  = co2mmr(i)*cpwpl*pnm(i,ntoplw)*pnm(i,ntoplw)
! Subprogram not used    end do
! Subprogram not used    do k=ntoplw,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          plh2o(i,k+1)  = plh2o(i,k) + rgsslp* &
! Subprogram not used                          (pnm(i,k+1)**2 - pnm(i,k)**2)*qnm(i,k)
! Subprogram not used          plco2(i,k+1)  = co2mmr(i)*cpwpl*pnm(i,k+1)**2
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Set the top and bottom intermediate level temperatures,
! Subprogram not used ! top level planck temperature and top layer temp**4.
! Subprogram not used !
! Subprogram not used ! Tint is lower interface temperature
! Subprogram not used ! (not available for bottom layer, so use ground temperature)
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       tint4(i,pverp)   = lwupcgs(i)/stebol_cgs
! Subprogram not used       tint(i,pverp)    = sqrt(sqrt(tint4(i,pverp)))
! Subprogram not used       tplnka(i,ntoplw) = tnm(i,ntoplw)
! Subprogram not used       tint(i,ntoplw)   = tplnka(i,ntoplw)
! Subprogram not used       tlayr4(i,ntoplw) = tplnka(i,ntoplw)**4
! Subprogram not used       tint4(i,ntoplw)  = tlayr4(i,ntoplw)
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Intermediate level temperatures are computed using temperature
! Subprogram not used ! at the full level below less dy*delta t,between the full level
! Subprogram not used !
! Subprogram not used    do k=ntoplw+1,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          dy = (piln(i,k) - pmln(i,k))/(pmln(i,k-1) - pmln(i,k))
! Subprogram not used          tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
! Subprogram not used          tint4(i,k) = tint(i,k)**4
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Now set the layer temp=full level temperatures and establish a
! Subprogram not used ! planck temperature for absorption (tplnka) which is the average
! Subprogram not used ! the intermediate level temperatures.  Note that tplnka is not
! Subprogram not used ! equal to the full level temperatures.
! Subprogram not used !
! Subprogram not used    do k=ntoplw+1,pverp
! Subprogram not used       do i=1,ncol
! Subprogram not used          tlayr(i,k)  = tnm(i,k-1)
! Subprogram not used          tlayr4(i,k) = tlayr(i,k)**4
! Subprogram not used          tplnka(i,k) = .5_r8*(tint(i,k) + tint(i,k-1))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Calculate tplank for emissivity calculation.
! Subprogram not used ! Assume isothermal tplnke i.e. all levels=ttop.
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       tplnke(i)       = tplnka(i,ntoplw)
! Subprogram not used       tlayr(i,ntoplw) = tint(i,ntoplw)
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Now compute h2o path fields:
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used !
! Subprogram not used ! Changed effective path temperature to std. Curtis-Godson form
! Subprogram not used !
! Subprogram not used       tcg(i,ntoplw) = rga*qnm(i,ntoplw)*pnm(i,ntoplw)*tnm(i,ntoplw)
! Subprogram not used       w(i,ntoplw)   = sslp * (plh2o(i,ntoplw)*2._r8) / pnm(i,ntoplw)
! Subprogram not used !
! Subprogram not used ! Hulst-Curtis-Godson scaling for H2O path
! Subprogram not used !
! Subprogram not used       wb(1,i,ntoplw) = w(i,ntoplw) * phi(tnm(i,ntoplw),1)
! Subprogram not used       wb(2,i,ntoplw) = w(i,ntoplw) * phi(tnm(i,ntoplw),2)
! Subprogram not used !
! Subprogram not used ! Hulst-Curtis-Godson scaling for effective pressure along H2O path
! Subprogram not used !
! Subprogram not used       plh2ob(1,i,ntoplw) = plh2o(i,ntoplw) * psi(tnm(i,ntoplw),1)
! Subprogram not used       plh2ob(2,i,ntoplw) = plh2o(i,ntoplw) * psi(tnm(i,ntoplw),2)
! Subprogram not used 
! Subprogram not used       s2c(i,ntoplw) = plh2o(i,ntoplw)*fh2oself(tnm(i,ntoplw))*qnm(i,ntoplw)*repsil
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    do k=ntoplw,pver
! Subprogram not used       do i=1,ncol
! Subprogram not used          dpnm       = pnm(i,k+1) - pnm(i,k)
! Subprogram not used          dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
! Subprogram not used          dw         = rga*qnm(i,k)*dpnm
! Subprogram not used          kp1        = k+1
! Subprogram not used          w(i,kp1)   = w(i,k) + dw
! Subprogram not used !
! Subprogram not used ! Hulst-Curtis-Godson scaling for H2O path
! Subprogram not used !
! Subprogram not used          wb(1,i,kp1) = wb(1,i,k) + dw * phi(tnm(i,k),1)
! Subprogram not used          wb(2,i,kp1) = wb(2,i,k) + dw * phi(tnm(i,k),2)
! Subprogram not used !
! Subprogram not used ! Hulst-Curtis-Godson scaling for effective pressure along H2O path
! Subprogram not used !
! Subprogram not used          dplh2o = plh2o(i,kp1) - plh2o(i,k)
! Subprogram not used 
! Subprogram not used          plh2ob(1,i,kp1) = plh2ob(1,i,k) + dplh2o * psi(tnm(i,k),1)
! Subprogram not used          plh2ob(2,i,kp1) = plh2ob(2,i,k) + dplh2o * psi(tnm(i,k),2)
! Subprogram not used !
! Subprogram not used ! Changed effective path temperature to std. Curtis-Godson form
! Subprogram not used !
! Subprogram not used          tcg(i,kp1) = tcg(i,k) + dw*tnm(i,k)
! Subprogram not used          s2c(i,kp1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)* &
! Subprogram not used                       fh2oself(tnm(i,k))*qnm(i,k)*repsil
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine radtpl

!====================================================================================

! Subprogram not used subroutine radae_init(gravx, epsilox, stebol, pstdx, mwdryx, mwco2x, mwo3x)
! Subprogram not used !
! Subprogram not used ! Initialize radae module data
! Subprogram not used !
! Subprogram not used    use pio,          only: file_desc_t, var_desc_t, pio_inq_dimid, pio_inquire_dimension, &
! Subprogram not used                            pio_inquire_variable, pio_inq_varid, pio_get_var, pio_nowrite, &
! Subprogram not used                            pio_max_var_dims, pio_max_name, pio_closefile
! Subprogram not used    use cam_pio_utils,only: cam_pio_openfile
! Subprogram not used    use ioFileMod,    only: getfil
! Subprogram not used    use filenames,    only: absems_data
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Input variables
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: gravx   ! Acceleration due to gravity (m/s**2)
! Subprogram not used    real(r8), intent(in) :: epsilox ! Ratio of mol. wght of H2O to dry air
! Subprogram not used    real(r8), intent(in) :: stebol  ! Stefan-Boltzmann's constant (MKS)
! Subprogram not used    real(r8), intent(in) :: pstdx   ! Standard pressure (pascals)
! Subprogram not used    real(r8), intent(in) :: mwdryx  ! Molecular weight of dry air 
! Subprogram not used    real(r8), intent(in) :: mwco2x  ! Molecular weight of carbon dioxide
! Subprogram not used    real(r8), intent(in) :: mwo3x   ! Molecular weight of ozone
! Subprogram not used !
! Subprogram not used !      Variables for loading absorptivity/emissivity
! Subprogram not used !
! Subprogram not used    type(file_desc_T) :: ncid_ae                ! NetCDF file id for abs/ems file
! Subprogram not used 
! Subprogram not used    integer pdimid                 ! pressure dimension id
! Subprogram not used    integer psize                  ! pressure dimension size
! Subprogram not used 
! Subprogram not used    integer tpdimid                ! path temperature dimension id
! Subprogram not used    integer tpsize                 ! path temperature size
! Subprogram not used 
! Subprogram not used    integer tedimid                ! emission temperature dimension id
! Subprogram not used    integer tesize                 ! emission temperature size
! Subprogram not used 
! Subprogram not used    integer udimid                 ! u (H2O path) dimension id
! Subprogram not used    integer usize                  ! u (H2O path) dimension size
! Subprogram not used 
! Subprogram not used    integer rhdimid                ! relative humidity dimension id
! Subprogram not used    integer rhsize                 ! relative humidity dimension size
! Subprogram not used 
! Subprogram not used    type(var_desc_t) ::    ah2onwid            ! var. id for non-wndw abs.
! Subprogram not used    type(var_desc_t) ::    eh2onwid            ! var. id for non-wndw ems.
! Subprogram not used    type(var_desc_t) ::    ah2owid             ! var. id for wndw abs. (adjacent layers)
! Subprogram not used    type(var_desc_t) :: cn_ah2owid             ! var. id for continuum trans. for wndw abs.
! Subprogram not used    type(var_desc_t) :: cn_eh2owid             ! var. id for continuum trans. for wndw ems.
! Subprogram not used    type(var_desc_t) :: ln_ah2owid             ! var. id for line trans. for wndw abs.
! Subprogram not used    type(var_desc_t) :: ln_eh2owid             ! var. id for line trans. for wndw ems.
! Subprogram not used    
! Subprogram not used    character*(PIO_MAX_NAME) tmpname! dummy variable for var/dim names
! Subprogram not used    character(len=256) locfn       ! local filename
! Subprogram not used    integer tmptype                ! dummy variable for variable type
! Subprogram not used    integer ndims                  ! number of dimensions
! Subprogram not used    integer dims(PIO_MAX_VAR_DIMS)  ! vector of dimension ids
! Subprogram not used    integer natt                   ! number of attributes
! Subprogram not used 
! Subprogram not used    integer ierr                  ! ierr flag returned from pio (pio handles errors internally so it is not checked)
! Subprogram not used !
! Subprogram not used ! Constants to set
! Subprogram not used !
! Subprogram not used    gravit     = gravx
! Subprogram not used    gravit_cgs = 100._r8*gravx
! Subprogram not used    rga        = 1._r8/gravit_cgs
! Subprogram not used    epsilo     = epsilox
! Subprogram not used    omeps      = 1._r8 - epsilo
! Subprogram not used    sslp       = 1.013250e6_r8
! Subprogram not used    stebol_cgs = 1.e3_r8*stebol
! Subprogram not used    rgsslp     = 0.5_r8/(gravit_cgs*sslp)
! Subprogram not used    dpfo3      = 2.5e-3_r8
! Subprogram not used    dpfco2     = 5.0e-3_r8
! Subprogram not used 
! Subprogram not used    p0         = pstdx*10.0_r8
! Subprogram not used    amd        = mwdryx
! Subprogram not used    amco2      = mwco2x
! Subprogram not used    mwo3       = mwo3x
! Subprogram not used !
! Subprogram not used ! Coefficients for h2o emissivity and absorptivity for overlap of H2O 
! Subprogram not used !    and trace gases.
! Subprogram not used !
! Subprogram not used    c16  = coefj(3,1)/coefj(2,1)
! Subprogram not used    c17  = coefk(3,1)/coefk(2,1)
! Subprogram not used    c26  = coefj(3,2)/coefj(2,2)
! Subprogram not used    c27  = coefk(3,2)/coefk(2,2)
! Subprogram not used    c28  = .5_r8
! Subprogram not used    c29  = .002053_r8
! Subprogram not used    c30  = .1_r8
! Subprogram not used    c31  = 3.0e-5_r8
! Subprogram not used !
! Subprogram not used ! Initialize further longwave constants referring to far wing
! Subprogram not used ! correction for overlap of H2O and trace gases; R&D refers to:
! Subprogram not used !
! Subprogram not used !            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
! Subprogram not used !            Emissivity and Absorptivity Formulation for Water Vapor
! Subprogram not used !            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
! Subprogram not used !
! Subprogram not used    fwcoef = .1_r8           ! See eq(33) R&D
! Subprogram not used    fwc1   = .30_r8          ! See eq(33) R&D
! Subprogram not used    fwc2   = 4.5_r8          ! See eq(33) and eq(34) in R&D
! Subprogram not used    fc1    = 2.6_r8          ! See eq(34) R&D
! Subprogram not used 
! Subprogram not used    call getfil(absems_data, locfn)
! Subprogram not used    call cam_pio_openfile(ncid_ae, locfn, PIO_NOWRITE)
! Subprogram not used    
! Subprogram not used    ierr =  pio_inq_dimid(ncid_ae, 'p', pdimid)
! Subprogram not used    ierr =  pio_inquire_dimension(ncid_ae, pdimid, len=psize)
! Subprogram not used    
! Subprogram not used    ierr =  pio_inq_dimid(ncid_ae, 'tp', tpdimid)
! Subprogram not used    ierr =  pio_inquire_dimension(ncid_ae, tpdimid, len=tpsize)
! Subprogram not used    
! Subprogram not used    ierr =  pio_inq_dimid(ncid_ae, 'te', tedimid)
! Subprogram not used    ierr =  pio_inquire_dimension(ncid_ae, tedimid, len=tesize)
! Subprogram not used    
! Subprogram not used    ierr =  pio_inq_dimid(ncid_ae, 'u', udimid)
! Subprogram not used    ierr =  pio_inquire_dimension(ncid_ae, udimid, len=usize)
! Subprogram not used    
! Subprogram not used    ierr =  pio_inq_dimid(ncid_ae, 'rh', rhdimid)
! Subprogram not used    ierr =  pio_inquire_dimension(ncid_ae, rhdimid, len=rhsize)
! Subprogram not used       
! Subprogram not used    if (psize    /= n_p  .or. &
! Subprogram not used         tpsize   /= n_tp .or. &
! Subprogram not used         usize    /= n_u  .or. &
! Subprogram not used         tesize   /= n_te .or. &
! Subprogram not used         rhsize   /= n_rh) then
! Subprogram not used       call endrun ('RADAEINI: dimensions for abs/ems do not match internal def.')
! Subprogram not used    endif
! Subprogram not used    
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'ah2onw',   ah2onwid)
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'eh2onw',   eh2onwid)
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'ah2ow',    ah2owid)
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'cn_ah2ow', cn_ah2owid)
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'cn_eh2ow', cn_eh2owid)
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'ln_ah2ow', ln_ah2owid)
! Subprogram not used    ierr =  pio_inq_varid(ncid_ae, 'ln_eh2ow', ln_eh2owid)
! Subprogram not used       
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, ah2onwid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: non-wndw abs. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, eh2onwid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: non-wndw ems. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, ah2owid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: window abs. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, cn_ah2owid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: cont. trans for abs. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, cn_eh2owid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: cont. trans. for ems. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, ln_ah2owid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: line trans for abs. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    ierr =  pio_inquire_variable(ncid_ae, ln_eh2owid, tmpname, tmptype, ndims, dims, natt)
! Subprogram not used    if (ndims /= 5 .or. &
! Subprogram not used         dims(1) /= pdimid    .or. &
! Subprogram not used         dims(2) /= tpdimid   .or. &
! Subprogram not used         dims(3) /= udimid    .or. &
! Subprogram not used         dims(4) /= tedimid   .or. &
! Subprogram not used         dims(5) /= rhdimid) then
! Subprogram not used       call endrun ('RADAEINI: line trans. for ems. in file /= internal def.')
! Subprogram not used    endif
! Subprogram not used    
! Subprogram not used    ierr =  pio_get_var (ncid_ae, ah2onwid,   ah2onw)
! Subprogram not used    ierr =  pio_get_var (ncid_ae, eh2onwid,   eh2onw)
! Subprogram not used    ierr =  pio_get_var (ncid_ae, ah2owid,    ah2ow)
! Subprogram not used    ierr =  pio_get_var (ncid_ae, cn_ah2owid, cn_ah2ow)
! Subprogram not used    ierr =  pio_get_var (ncid_ae, cn_eh2owid, cn_eh2ow)
! Subprogram not used    ierr =  pio_get_var (ncid_ae, ln_ah2owid, ln_ah2ow)
! Subprogram not used    ierr =  pio_get_var (ncid_ae, ln_eh2owid, ln_eh2ow)
! Subprogram not used       
! Subprogram not used    call pio_closefile(ncid_ae)
! Subprogram not used 
! Subprogram not used end subroutine radae_init

!====================================================================================

subroutine initialize_radbuffer
!
! Initialize radiation buffer data
!

  use ref_pres, only : pref_mid
  use phys_control, only : phys_getopts

  character(len=16) :: radiation_scheme
  integer :: k

! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
! longwave cooling to about 80 km (1 Pa)
   if (pref_mid(1) .lt. 0.1_r8) then
      do k = 1, plev
         if (pref_mid(k) .lt. 1._r8) ntoplw  = k
      end do
   else
      ntoplw  = 1
   end if
   if (masterproc) then
      write(iulog,*) 'INITIALIZE_RADBUFFER: ntoplw =',ntoplw, ' pressure:',pref_mid(ntoplw)
   endif
 
   call phys_getopts(radiation_scheme_out=radiation_scheme)
   
  if(radiation_scheme.eq.'camrt') then
     allocate (abstot_3d(pcols,ntoplw:pverp,ntoplw:pverp,begchunk:endchunk))
     allocate (absnxt_3d(pcols,pver,4,begchunk:endchunk))
     allocate (emstot_3d(pcols,pverp,begchunk:endchunk))
     abstot_3d(:,:,:,:) = posinf
     absnxt_3d(:,:,:,:) = posinf
     emstot_3d(:,:,:) = posinf
  end if
  return
end subroutine initialize_radbuffer

!====================================================================================

! Subprogram not used subroutine radoz2(ncol, o3, pint, plol, plos)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Computes the path length integrals to the model interfaces given the
! Subprogram not used ! ozone volume mixing ratio
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! <Describe the algorithm(s) used in the routine.> 
! Subprogram not used ! <Also include any applicable external references.> 
! Subprogram not used ! 
! Subprogram not used ! Author: CCM1, CMS Contact J. Kiehl
! Subprogram not used ! 
! Subprogram not used !------------------------------Input arguments--------------------------
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                 ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: o3(pcols,pver)      ! ozone mass mixing ratio
! Subprogram not used    real(r8), intent(in) :: pint(pcols,pverp)   ! Model interface pressures
! Subprogram not used !
! Subprogram not used !----------------------------Output arguments---------------------------
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: plol(pcols,pverp)   ! Ozone prs weighted path length (cm)
! Subprogram not used    real(r8), intent(out) :: plos(pcols,pverp)   ! Ozone path length (cm)
! Subprogram not used !
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    integer i                ! longitude index
! Subprogram not used    integer k                ! level index
! Subprogram not used 
! Subprogram not used    real(r8) :: v0    ! Volume of a gas at stp (m**3/kmol)
! Subprogram not used    real(r8) :: p0    ! Standard pressure (pascals)
! Subprogram not used    real(r8) :: cplos ! constant for ozone path length integral
! Subprogram not used    real(r8) :: cplol ! constant for ozone path length integral
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !*******************************************************************
! Subprogram not used ! These hardwired constants need to be replaced with common values.
! Subprogram not used ! They are here for testing infrastructure changes that should not
! Subprogram not used ! change answers.
! Subprogram not used ! Constants for ozone path integrals (multiplication by 100 for unit
! Subprogram not used ! conversion to cgs from mks):
! Subprogram not used !
! Subprogram not used    v0    = 22.4136_r8         ! Volume of a gas at stp (m**3/kmol)
! Subprogram not used    p0    = 0.1_r8*sslp        ! Standard pressure (pascals)
! Subprogram not used    cplos = v0/(mwo3*gravit)       *100.0_r8
! Subprogram not used    cplol = v0/(mwo3*gravit*p0)*0.5_r8*100.0_r8
! Subprogram not used !*******************************************************************
! Subprogram not used !
! Subprogram not used ! Evaluate the ozone path length integrals to interfaces;
! Subprogram not used ! factors of .1 and .01 to convert pressures from cgs to mks:
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       plos(i,ntoplw) = 0.1_r8 *cplos*o3(i,ntoplw)*pint(i,ntoplw)
! Subprogram not used       plol(i,ntoplw) = 0.01_r8*cplol*o3(i,ntoplw)*pint(i,ntoplw)*pint(i,ntoplw)
! Subprogram not used    end do
! Subprogram not used    do k=ntoplw+1,pverp
! Subprogram not used       do i=1,ncol
! Subprogram not used          plos(i,k) = plos(i,k-1) + 0.1_r8*cplos*o3(i,k-1)*(pint(i,k) - pint(i,k-1))
! Subprogram not used          plol(i,k) = plol(i,k-1) + 0.01_r8*cplol*o3(i,k-1)* &
! Subprogram not used                     (pint(i,k)*pint(i,k) - pint(i,k-1)*pint(i,k-1))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine radoz2

!====================================================================================

! Subprogram not used subroutine trcpth(ncol                                        , &
! Subprogram not used                   tnm     ,pnm     ,cfc11   ,cfc12   ,n2o     , &
! Subprogram not used                   ch4     ,qnm     ,ucfc11  ,ucfc12  ,un2o0   , &
! Subprogram not used                   un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
! Subprogram not used                   uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
! Subprogram not used                   bch4    ,uptype  ,co2mmr)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Calculate path lengths and pressure factors for CH4, N2O, CFC11
! Subprogram not used ! and CFC12.
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! See CCM3 description for details
! Subprogram not used ! 
! Subprogram not used ! Author: J. Kiehl
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: tnm(pcols,pver)      ! Model level temperatures
! Subprogram not used    real(r8), intent(in) :: pnm(pcols,pverp)     ! Pres. at model interfaces (dynes/cm2)
! Subprogram not used    real(r8), intent(in) :: qnm(pcols,pver)      ! h2o specific humidity
! Subprogram not used    real(r8), intent(in) :: cfc11(pcols,pver)    ! CFC11 mass mixing ratio
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: cfc12(pcols,pver)    ! CFC12 mass mixing ratio
! Subprogram not used    real(r8), intent(in) :: n2o(pcols,pver)      ! N2O mass mixing ratio
! Subprogram not used    real(r8), intent(in) :: ch4(pcols,pver)      ! CH4 mass mixing ratio
! Subprogram not used    real(r8), intent(in) :: co2mmr(pcols)        ! co2 column mean mass mixing ratio
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: ucfc11(pcols,pverp)  ! CFC11 path length
! Subprogram not used    real(r8), intent(out) :: ucfc12(pcols,pverp)  ! CFC12 path length
! Subprogram not used    real(r8), intent(out) :: un2o0(pcols,pverp)   ! N2O path length
! Subprogram not used    real(r8), intent(out) :: un2o1(pcols,pverp)   ! N2O path length (hot band)
! Subprogram not used    real(r8), intent(out) :: uch4(pcols,pverp)    ! CH4 path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: uco211(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(out) :: uco212(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(out) :: uco213(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(out) :: uco221(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(out) :: uco222(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: uco223(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(out) :: bn2o0(pcols,pverp)   ! pressure factor for n2o
! Subprogram not used    real(r8), intent(out) :: bn2o1(pcols,pverp)   ! pressure factor for n2o
! Subprogram not used    real(r8), intent(out) :: bch4(pcols,pverp)    ! pressure factor for ch4
! Subprogram not used    real(r8), intent(out) :: uptype(pcols,pverp)  ! p-type continuum path length
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used !
! Subprogram not used    integer   i               ! Longitude index
! Subprogram not used    integer   k               ! Level index
! Subprogram not used !
! Subprogram not used    real(r8) co2fac(pcols,1)      ! co2 factor
! Subprogram not used    real(r8) alpha1(pcols)        ! stimulated emission term
! Subprogram not used    real(r8) alpha2(pcols)        ! stimulated emission term
! Subprogram not used    real(r8) rt(pcols)            ! reciprocal of local temperature
! Subprogram not used    real(r8) rsqrt(pcols)         ! reciprocal of sqrt of temp
! Subprogram not used !
! Subprogram not used    real(r8) pbar(pcols)          ! mean pressure
! Subprogram not used    real(r8) dpnm(pcols)          ! difference in pressure
! Subprogram not used    real(r8) diff                 ! diffusivity factor
! Subprogram not used !
! Subprogram not used !--------------------------Data Statements------------------------------
! Subprogram not used !
! Subprogram not used    data diff /1.66_r8/
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !  Calculate path lengths for the trace gases at model top
! Subprogram not used !
! Subprogram not used 
! Subprogram not used    do i = 1,ncol
! Subprogram not used       ucfc11(i,ntoplw) = 1.8_r8 * cfc11(i,ntoplw) * pnm(i,ntoplw) * rga
! Subprogram not used       ucfc12(i,ntoplw) = 1.8_r8 * cfc12(i,ntoplw) * pnm(i,ntoplw) * rga
! Subprogram not used       un2o0(i,ntoplw) = diff * 1.02346e5_r8 * n2o(i,ntoplw) * pnm(i,ntoplw) * rga / sqrt(tnm(i,ntoplw))
! Subprogram not used       un2o1(i,ntoplw) = diff * 2.01909_r8 * un2o0(i,ntoplw) * exp(-847.36_r8/tnm(i,ntoplw))
! Subprogram not used       uch4(i,ntoplw)  = diff * 8.60957e4_r8 * ch4(i,ntoplw) * pnm(i,ntoplw) * rga / sqrt(tnm(i,ntoplw))
! Subprogram not used       co2fac(i,1)     = diff * co2mmr(i) * pnm(i,ntoplw) * rga
! Subprogram not used       alpha1(i) = (1.0_r8 - exp(-1540.0_r8/tnm(i,ntoplw)))**3.0_r8/sqrt(tnm(i,ntoplw))
! Subprogram not used       alpha2(i) = (1.0_r8 - exp(-1360.0_r8/tnm(i,ntoplw)))**3.0_r8/sqrt(tnm(i,ntoplw))
! Subprogram not used       uco211(i,ntoplw) = 3.42217e3_r8 * co2fac(i,1) * alpha1(i) * exp(-1849.7_r8/tnm(i,ntoplw))
! Subprogram not used       uco212(i,ntoplw) = 6.02454e3_r8 * co2fac(i,1) * alpha1(i) * exp(-2782.1_r8/tnm(i,ntoplw))
! Subprogram not used       uco213(i,ntoplw) = 5.53143e3_r8 * co2fac(i,1) * alpha1(i) * exp(-3723.2_r8/tnm(i,ntoplw))
! Subprogram not used       uco221(i,ntoplw) = 3.88984e3_r8 * co2fac(i,1) * alpha2(i) * exp(-1997.6_r8/tnm(i,ntoplw))
! Subprogram not used       uco222(i,ntoplw) = 3.67108e3_r8 * co2fac(i,1) * alpha2(i) * exp(-3843.8_r8/tnm(i,ntoplw))
! Subprogram not used       uco223(i,ntoplw) = 6.50642e3_r8 * co2fac(i,1) * alpha2(i) * exp(-2989.7_r8/tnm(i,ntoplw))
! Subprogram not used       bn2o0(i,ntoplw) = diff * 19.399_r8 * pnm(i,ntoplw)**2.0_r8 * n2o(i,ntoplw) * &
! Subprogram not used                    1.02346e5_r8 * rga / (sslp*tnm(i,ntoplw))
! Subprogram not used       bn2o1(i,ntoplw) = bn2o0(i,ntoplw) * exp(-847.36_r8/tnm(i,ntoplw)) * 2.06646e5_r8
! Subprogram not used       bch4(i,ntoplw) = diff * 2.94449_r8 * ch4(i,ntoplw) * pnm(i,ntoplw)**2.0_r8 * rga * &
! Subprogram not used                   8.60957e4_r8 / (sslp*tnm(i,ntoplw))
! Subprogram not used       uptype(i,ntoplw) = diff * qnm(i,ntoplw) * pnm(i,ntoplw)**2.0_r8 *  &
! Subprogram not used                     exp(1800.0_r8*(1.0_r8/tnm(i,ntoplw) - 1.0_r8/296.0_r8)) * rga / sslp
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Calculate trace gas path lengths through model atmosphere
! Subprogram not used !
! Subprogram not used    do k = ntoplw,pver
! Subprogram not used       do i = 1,ncol
! Subprogram not used          rt(i) = 1._r8/tnm(i,k)
! Subprogram not used          rsqrt(i) = sqrt(rt(i))
! Subprogram not used          pbar(i) = 0.5_r8 * (pnm(i,k+1) + pnm(i,k)) / sslp
! Subprogram not used          dpnm(i) = (pnm(i,k+1) - pnm(i,k)) * rga
! Subprogram not used          alpha1(i) = diff * rsqrt(i) * (1.0_r8 - exp(-1540.0_r8/tnm(i,k)))**3.0_r8
! Subprogram not used          alpha2(i) = diff * rsqrt(i) * (1.0_r8 - exp(-1360.0_r8/tnm(i,k)))**3.0_r8
! Subprogram not used          ucfc11(i,k+1) = ucfc11(i,k) +  1.8_r8 * cfc11(i,k) * dpnm(i)
! Subprogram not used          ucfc12(i,k+1) = ucfc12(i,k) +  1.8_r8 * cfc12(i,k) * dpnm(i)
! Subprogram not used          un2o0(i,k+1) = un2o0(i,k) + diff * 1.02346e5_r8 * n2o(i,k) * rsqrt(i) * dpnm(i)
! Subprogram not used          un2o1(i,k+1) = un2o1(i,k) + diff * 2.06646e5_r8 * n2o(i,k) * &
! Subprogram not used                         rsqrt(i) * exp(-847.36_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          uch4(i,k+1) = uch4(i,k) + diff * 8.60957e4_r8 * ch4(i,k) * rsqrt(i) * dpnm(i)
! Subprogram not used          uco211(i,k+1) = uco211(i,k) + 1.15_r8*3.42217e3_r8 * alpha1(i) * &
! Subprogram not used                          co2mmr(i) * exp(-1849.7_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          uco212(i,k+1) = uco212(i,k) + 1.15_r8*6.02454e3_r8 * alpha1(i) * &
! Subprogram not used                          co2mmr(i) * exp(-2782.1_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          uco213(i,k+1) = uco213(i,k) + 1.15_r8*5.53143e3_r8 * alpha1(i) * &
! Subprogram not used                          co2mmr(i) * exp(-3723.2_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          uco221(i,k+1) = uco221(i,k) + 1.15_r8*3.88984e3_r8 * alpha2(i) * &
! Subprogram not used                          co2mmr(i) * exp(-1997.6_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          uco222(i,k+1) = uco222(i,k) + 1.15_r8*3.67108e3_r8 * alpha2(i) * &
! Subprogram not used                          co2mmr(i) * exp(-3843.8_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          uco223(i,k+1) = uco223(i,k) + 1.15_r8*6.50642e3_r8 * alpha2(i) * &
! Subprogram not used                          co2mmr(i) * exp(-2989.7_r8/tnm(i,k)) * dpnm(i)
! Subprogram not used          bn2o0(i,k+1) = bn2o0(i,k) + diff * 19.399_r8 * pbar(i) * rt(i) &
! Subprogram not used                         * 1.02346e5_r8 * n2o(i,k) * dpnm(i)
! Subprogram not used          bn2o1(i,k+1) = bn2o1(i,k) + diff * 19.399_r8 * pbar(i) * rt(i) &
! Subprogram not used                         * 2.06646e5_r8 * exp(-847.36_r8/tnm(i,k)) * n2o(i,k)*dpnm(i)
! Subprogram not used          bch4(i,k+1) = bch4(i,k) + diff * 2.94449_r8 * rt(i) * pbar(i) &
! Subprogram not used                        * 8.60957e4_r8 * ch4(i,k) * dpnm(i)
! Subprogram not used          uptype(i,k+1) = uptype(i,k) + diff *qnm(i,k) * &
! Subprogram not used                          exp(1800.0_r8*(1.0_r8/tnm(i,k) - 1.0_r8/296.0_r8)) * pbar(i) * dpnm(i)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used end subroutine trcpth



!====================================================================================
! Private Interfaces
!====================================================================================

! Subprogram not used function fh2oself( temp )
! Subprogram not used !
! Subprogram not used ! Short function for H2O self-continuum temperature factor in
! Subprogram not used !   calculation of effective H2O self-continuum path length
! Subprogram not used !
! Subprogram not used ! H2O Continuum: CKD 2.4
! Subprogram not used ! Code for continuum: GENLN3
! Subprogram not used ! Reference: Edwards, D.P., 1992: GENLN2, A General Line-by-Line Atmospheric
! Subprogram not used !                     Transmittance and Radiance Model, Version 3.0 Description
! Subprogram not used !                     and Users Guide, NCAR/TN-367+STR, 147 pp.
! Subprogram not used !
! Subprogram not used ! In GENLN, the temperature scaling of the self-continuum is handled
! Subprogram not used !    by exponential interpolation/extrapolation from observations at
! Subprogram not used !    260K and 296K by:
! Subprogram not used !
! Subprogram not used !         TFAC =  (T(IPATH) - 296.0)/(260.0 - 296.0) 
! Subprogram not used !         CSFFT = CSFF296*(CSFF260/CSFF296)**TFAC
! Subprogram not used !
! Subprogram not used ! For 800-1200 cm^-1, (CSFF260/CSFF296) ranges from ~2.1 to ~1.9
! Subprogram not used !     with increasing wavenumber.  The ratio <CSFF260>/<CSFF296>,
! Subprogram not used !     where <> indicates average over wavenumber, is ~2.07
! Subprogram not used ! 
! Subprogram not used ! fh2oself is (<CSFF260>/<CSFF296>)**TFAC
! Subprogram not used !
! Subprogram not used    real(r8),intent(in) :: temp     ! path temperature
! Subprogram not used    real(r8) fh2oself               ! mean ratio of self-continuum at temp and 296K
! Subprogram not used 
! Subprogram not used    fh2oself = 2.0727484_r8**((296.0_r8 - temp) / 36.0_r8)
! Subprogram not used end function fh2oself

!====================================================================================

! Subprogram not used function phi(tpx,iband)
! Subprogram not used !
! Subprogram not used ! History: First version for Hitran 1996 (C/H/E)
! Subprogram not used !          Current version for Hitran 2000 (C/LT/E)
! Subprogram not used ! Short function for Hulst-Curtis-Godson temperature factors for
! Subprogram not used !   computing effective H2O path
! Subprogram not used ! Line data for H2O: Hitran 2000, plus H2O patches v11.0 for 1341 missing
! Subprogram not used !                    lines between 500 and 2820 cm^-1.
! Subprogram not used !                    See cfa-www.harvard.edu/HITRAN
! Subprogram not used ! Isotopes of H2O: all
! Subprogram not used ! Line widths: air-broadened only (self set to 0)
! Subprogram not used ! Code for line strengths and widths: GENLN3
! Subprogram not used ! Reference: Edwards, D.P., 1992: GENLN2, A General Line-by-Line Atmospheric
! Subprogram not used !                     Transmittance and Radiance Model, Version 3.0 Description
! Subprogram not used !                     and Users Guide, NCAR/TN-367+STR, 147 pp.
! Subprogram not used !
! Subprogram not used ! Note: functions have been normalized by dividing by their values at
! Subprogram not used !       a path temperature of 160K
! Subprogram not used !
! Subprogram not used ! spectral intervals:
! Subprogram not used !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
! Subprogram not used !   2 = 800-1200 cm^-1      
! Subprogram not used !
! Subprogram not used ! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
! Subprogram not used !           2nd edition, Oxford University Press, 1989.
! Subprogram not used ! Phi: function for H2O path
! Subprogram not used !      eq. 6.25, p. 228
! Subprogram not used !
! Subprogram not used    real(r8),intent(in):: tpx      ! path temperature
! Subprogram not used    integer, intent(in):: iband    ! band to process
! Subprogram not used    real(r8) phi                   ! phi for given band
! Subprogram not used    real(r8),parameter ::  phi_r0(nbands) = (/ 9.60917711E-01_r8, -2.21031342E+01_r8/)
! Subprogram not used    real(r8),parameter ::  phi_r1(nbands) = (/ 4.86076751E-04_r8,  4.24062610E-01_r8/)
! Subprogram not used    real(r8),parameter ::  phi_r2(nbands) = (/-1.84806265E-06_r8, -2.95543415E-03_r8/)
! Subprogram not used    real(r8),parameter ::  phi_r3(nbands) = (/ 2.11239959E-09_r8,  7.52470896E-06_r8/)
! Subprogram not used 
! Subprogram not used    phi = (((phi_r3(iband) * tpx) + phi_r2(iband)) * tpx + phi_r1(iband)) &
! Subprogram not used           * tpx + phi_r0(iband)
! Subprogram not used end function phi

!====================================================================================

! Subprogram not used function psi(tpx,iband)
! Subprogram not used !
! Subprogram not used ! History: First version for Hitran 1996 (C/H/E)
! Subprogram not used !          Current version for Hitran 2000 (C/LT/E)
! Subprogram not used ! Short function for Hulst-Curtis-Godson temperature factors for
! Subprogram not used !   computing effective H2O path
! Subprogram not used ! Line data for H2O: Hitran 2000, plus H2O patches v11.0 for 1341 missing
! Subprogram not used !                    lines between 500 and 2820 cm^-1. 
! Subprogram not used !                    See cfa-www.harvard.edu/HITRAN
! Subprogram not used ! Isotopes of H2O: all
! Subprogram not used ! Line widths: air-broadened only (self set to 0)
! Subprogram not used ! Code for line strengths and widths: GENLN3
! Subprogram not used ! Reference: Edwards, D.P., 1992: GENLN2, A General Line-by-Line Atmospheric
! Subprogram not used !                     Transmittance and Radiance Model, Version 3.0 Description
! Subprogram not used !                     and Users Guide, NCAR/TN-367+STR, 147 pp.
! Subprogram not used !
! Subprogram not used ! Note: functions have been normalized by dividing by their values at
! Subprogram not used !       a path temperature of 160K
! Subprogram not used !
! Subprogram not used ! spectral intervals:
! Subprogram not used !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
! Subprogram not used !   2 = 800-1200 cm^-1      
! Subprogram not used !
! Subprogram not used ! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
! Subprogram not used !           2nd edition, Oxford University Press, 1989.
! Subprogram not used ! Psi: function for pressure along path
! Subprogram not used !      eq. 6.30, p. 228
! Subprogram not used !
! Subprogram not used    real(r8),intent(in):: tpx      ! path temperature
! Subprogram not used    integer, intent(in):: iband    ! band to process
! Subprogram not used    real(r8) psi                   ! psi for given band
! Subprogram not used    real(r8),parameter ::  psi_r0(nbands) = (/ 5.65308452E-01_r8, -7.30087891E+01_r8/)
! Subprogram not used    real(r8),parameter ::  psi_r1(nbands) = (/ 4.07519005E-03_r8,  1.22199547E+00_r8/)
! Subprogram not used    real(r8),parameter ::  psi_r2(nbands) = (/-1.04347237E-05_r8, -7.12256227E-03_r8/)
! Subprogram not used    real(r8),parameter ::  psi_r3(nbands) = (/ 1.23765354E-08_r8,  1.47852825E-05_r8/)
! Subprogram not used 
! Subprogram not used    psi = (((psi_r3(iband) * tpx) + psi_r2(iband)) * tpx + psi_r1(iband)) * tpx + psi_r0(iband)
! Subprogram not used 
! Subprogram not used end function psi

!====================================================================================

! Subprogram not used subroutine trcab(ncol    ,                                     &
! Subprogram not used                  k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
! Subprogram not used                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
! Subprogram not used                  uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
! Subprogram not used                  bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
! Subprogram not used                  s2c     ,uptype  ,dplh2o  ,abplnk1 ,tco2    , &
! Subprogram not used                  th2o    ,to3     ,abstrc  , &
! Subprogram not used                  aer_trn_ttl)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
! Subprogram not used ! CFC12.
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! See CCM3 description for equations.
! Subprogram not used ! 
! Subprogram not used ! Author: J. Kiehl
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                     ! number of atmospheric columns
! Subprogram not used    integer, intent(in) :: k1,k2                    ! level indices
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: to3co2(pcols)           ! pressure weighted temperature
! Subprogram not used    real(r8), intent(in) :: pnm(pcols,pverp)        ! interface pressures
! Subprogram not used    real(r8), intent(in) :: ucfc11(pcols,pverp)     ! CFC11 path length
! Subprogram not used    real(r8), intent(in) :: ucfc12(pcols,pverp)     ! CFC12 path length
! Subprogram not used    real(r8), intent(in) :: un2o0(pcols,pverp)      ! N2O path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: un2o1(pcols,pverp)      ! N2O path length (hot band)
! Subprogram not used    real(r8), intent(in) :: uch4(pcols,pverp)       ! CH4 path length
! Subprogram not used    real(r8), intent(in) :: uco211(pcols,pverp)     ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco212(pcols,pverp)     ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco213(pcols,pverp)     ! CO2 9.4 micron band path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: uco221(pcols,pverp)     ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco222(pcols,pverp)     ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco223(pcols,pverp)     ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
! Subprogram not used    real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
! Subprogram not used    real(r8), intent(in) :: dw(pcols)               ! h2o path length
! Subprogram not used    real(r8), intent(in) :: pnew(pcols)             ! pressure
! Subprogram not used    real(r8), intent(in) :: s2c(pcols,pverp)        ! continuum path length
! Subprogram not used    real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type h2o path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: dplh2o(pcols)           ! p squared h2o path length
! Subprogram not used    real(r8), intent(in) :: abplnk1(14,pcols,pverp) ! Planck factor
! Subprogram not used    real(r8), intent(in) :: tco2(pcols)             ! co2 transmission factor
! Subprogram not used    real(r8), intent(in) :: th2o(pcols)             ! h2o transmission factor
! Subprogram not used    real(r8), intent(in) :: to3(pcols)              ! o3 transmission factor
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) ! aer trn.
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !  Output Arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: abstrc(pcols)           ! total trace gas absorptivity
! Subprogram not used !
! Subprogram not used !--------------------------Local Variables------------------------------
! Subprogram not used !
! Subprogram not used    integer  i,l                     ! loop counters
! Subprogram not used 
! Subprogram not used    real(r8) sqti(pcols)             ! square root of mean temp
! Subprogram not used    real(r8) du1                     ! cfc11 path length
! Subprogram not used    real(r8) du2                     ! cfc12 path length
! Subprogram not used    real(r8) acfc1                   ! cfc11 absorptivity 798 cm-1
! Subprogram not used    real(r8) acfc2                   ! cfc11 absorptivity 846 cm-1
! Subprogram not used !
! Subprogram not used    real(r8) acfc3                   ! cfc11 absorptivity 933 cm-1
! Subprogram not used    real(r8) acfc4                   ! cfc11 absorptivity 1085 cm-1
! Subprogram not used    real(r8) acfc5                   ! cfc12 absorptivity 889 cm-1
! Subprogram not used    real(r8) acfc6                   ! cfc12 absorptivity 923 cm-1
! Subprogram not used    real(r8) acfc7                   ! cfc12 absorptivity 1102 cm-1
! Subprogram not used !
! Subprogram not used    real(r8) acfc8                   ! cfc12 absorptivity 1161 cm-1
! Subprogram not used    real(r8) du01                    ! n2o path length
! Subprogram not used    real(r8) dbeta01                 ! n2o pressure factor
! Subprogram not used    real(r8) dbeta11                 !         "
! Subprogram not used    real(r8) an2o1                   ! absorptivity of 1285 cm-1 n2o band
! Subprogram not used !
! Subprogram not used    real(r8) du02                    ! n2o path length
! Subprogram not used    real(r8) dbeta02                 ! n2o pressure factor
! Subprogram not used    real(r8) an2o2                   ! absorptivity of 589 cm-1 n2o band
! Subprogram not used    real(r8) du03                    ! n2o path length
! Subprogram not used    real(r8) dbeta03                 ! n2o pressure factor
! Subprogram not used !
! Subprogram not used    real(r8) an2o3                   ! absorptivity of 1168 cm-1 n2o band
! Subprogram not used    real(r8) duch4                   ! ch4 path length
! Subprogram not used    real(r8) dbetac                  ! ch4 pressure factor
! Subprogram not used    real(r8) ach4                    ! absorptivity of 1306 cm-1 ch4 band
! Subprogram not used    real(r8) du11                    ! co2 path length
! Subprogram not used !
! Subprogram not used    real(r8) du12                    !       "
! Subprogram not used    real(r8) du13                    !       "
! Subprogram not used    real(r8) dbetc1                  ! co2 pressure factor
! Subprogram not used    real(r8) dbetc2                  ! co2 pressure factor
! Subprogram not used    real(r8) aco21                   ! absorptivity of 1064 cm-1 band
! Subprogram not used !
! Subprogram not used    real(r8) du21                    ! co2 path length
! Subprogram not used    real(r8) du22                    !       "
! Subprogram not used    real(r8) du23                    !       "
! Subprogram not used    real(r8) aco22                   ! absorptivity of 961 cm-1 band
! Subprogram not used    real(r8) tt(pcols)               ! temp. factor for h2o overlap factor
! Subprogram not used !
! Subprogram not used    real(r8) psi1                    !                 "
! Subprogram not used    real(r8) phi1                    !                 "
! Subprogram not used    real(r8) p1                      ! h2o overlap factor
! Subprogram not used    real(r8) w1                      !        "
! Subprogram not used    real(r8) ds2c(pcols)             ! continuum path length
! Subprogram not used !
! Subprogram not used    real(r8) duptyp(pcols)           ! p-type path length
! Subprogram not used    real(r8) tw(pcols,6)             ! h2o transmission factor
! Subprogram not used !   real(r8) g1(6)                   !         "
! Subprogram not used !   real(r8) g2(6)                   !         "
! Subprogram not used !   real(r8) g3(6)                   !         "
! Subprogram not used !
! Subprogram not used !   real(r8) g4(6)                   !         "
! Subprogram not used !   real(r8) ab(6)                   ! h2o temp. factor
! Subprogram not used !   real(r8) bb(6)                   !         "
! Subprogram not used !   real(r8) abp(6)                  !         "
! Subprogram not used !   real(r8) bbp(6)                  !         "
! Subprogram not used !
! Subprogram not used    real(r8) tcfc3                   ! transmission for cfc11 band
! Subprogram not used    real(r8) tcfc4                   ! transmission for cfc11 band
! Subprogram not used    real(r8) tcfc6                   ! transmission for cfc12 band
! Subprogram not used    real(r8) tcfc7                   ! transmission for cfc12 band
! Subprogram not used    real(r8) tcfc8                   ! transmission for cfc12 band
! Subprogram not used !
! Subprogram not used    real(r8) tlw                     ! h2o transmission
! Subprogram not used    real(r8) tch4                    ! ch4 transmission
! Subprogram not used !
! Subprogram not used !--------------------------Data Statements------------------------------
! Subprogram not used !
! Subprogram not used !   data g1 /0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/
! Subprogram not used !   data g2 /14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/
! Subprogram not used !   data g3 /26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/
! Subprogram not used !   data g4 /0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/
! Subprogram not used !   data ab /3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/
! Subprogram not used !   data bb /-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/
! Subprogram not used !   data abp/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/
! Subprogram not used !   data bbp/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/
! Subprogram not used !
! Subprogram not used !--------------------------Statement Functions--------------------------
! Subprogram not used !
! Subprogram not used    real(r8) func, u, b
! Subprogram not used    func(u,b) = u/sqrt(4.0_r8 + u*(1.0_r8 + 1.0_r8 / b))
! Subprogram not used !
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    do i = 1,ncol
! Subprogram not used       sqti(i) = sqrt(to3co2(i))
! Subprogram not used !
! Subprogram not used ! h2o transmission
! Subprogram not used !
! Subprogram not used       tt(i) = abs(to3co2(i) - 250.0_r8)
! Subprogram not used       ds2c(i) = abs(s2c(i,k1) - s2c(i,k2))
! Subprogram not used       duptyp(i) = abs(uptype(i,k1) - uptype(i,k2))
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do l = 1,6
! Subprogram not used       do i = 1,ncol
! Subprogram not used          psi1 = exp(abp(l)*tt(i) + bbp(l)*tt(i)*tt(i))
! Subprogram not used          phi1 = exp(ab(l)*tt(i) + bb(l)*tt(i)*tt(i))
! Subprogram not used          p1 = pnew(i)*(psi1/phi1)/sslp
! Subprogram not used          w1 = dw(i)*phi1
! Subprogram not used          tw(i,l) = exp(-g1(l)*p1*(sqrt(1.0_r8 + g2(l)*(w1/p1)) - 1.0_r8) - &
! Subprogram not used                    g3(l)*ds2c(i)-g4(l)*duptyp(i))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       tw(i,1)=tw(i,1)*(0.7_r8*aer_trn_ttl(i,k1,k2,idx_LW_0650_0800)+&! l=1: 0750--0820 cm-1
! Subprogram not used                        0.3_r8*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000)) 
! Subprogram not used       tw(i,2)=tw(i,2)*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000) ! l=2: 0820--0880 cm-1
! Subprogram not used       tw(i,3)=tw(i,3)*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000) ! l=3: 0880--0900 cm-1
! Subprogram not used       tw(i,4)=tw(i,4)*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000) ! l=4: 0900--1000 cm-1
! Subprogram not used       tw(i,5)=tw(i,5)*aer_trn_ttl(i,k1,k2,idx_LW_1000_1200) ! l=5: 1000--1120 cm-1
! Subprogram not used       tw(i,6)=tw(i,6)*aer_trn_ttl(i,k1,k2,idx_LW_1000_1200) ! l=6: 1120--1170 cm-1
! Subprogram not used    end do                    ! end loop over lon
! Subprogram not used    do i = 1,ncol
! Subprogram not used       du1 = abs(ucfc11(i,k1) - ucfc11(i,k2))
! Subprogram not used       du2 = abs(ucfc12(i,k1) - ucfc12(i,k2))
! Subprogram not used !
! Subprogram not used ! cfc transmissions
! Subprogram not used !
! Subprogram not used       tcfc3 = exp(-175.005_r8*du1)
! Subprogram not used       tcfc4 = exp(-1202.18_r8*du1)
! Subprogram not used       tcfc6 = exp(-5786.73_r8*du2)
! Subprogram not used       tcfc7 = exp(-2873.51_r8*du2)
! Subprogram not used       tcfc8 = exp(-2085.59_r8*du2)
! Subprogram not used !
! Subprogram not used ! Absorptivity for CFC11 bands
! Subprogram not used !
! Subprogram not used       acfc1 =  50.0_r8*(1.0_r8 - exp(-54.09_r8*du1))*tw(i,1)*abplnk1(7,i,k2)
! Subprogram not used       acfc2 =  60.0_r8*(1.0_r8 - exp(-5130.03_r8*du1))*tw(i,2)*abplnk1(8,i,k2)
! Subprogram not used       acfc3 =  60.0_r8*(1.0_r8 - tcfc3)*tw(i,4)*tcfc6*abplnk1(9,i,k2)
! Subprogram not used       acfc4 = 100.0_r8*(1.0_r8 - tcfc4)*tw(i,5)*abplnk1(10,i,k2)
! Subprogram not used !
! Subprogram not used ! Absorptivity for CFC12 bands
! Subprogram not used !
! Subprogram not used       acfc5 = 45.0_r8*(1.0_r8 - exp(-1272.35_r8*du2))*tw(i,3)*abplnk1(11,i,k2)
! Subprogram not used       acfc6 = 50.0_r8*(1.0_r8 - tcfc6)* tw(i,4) * abplnk1(12,i,k2)
! Subprogram not used       acfc7 = 80.0_r8*(1.0_r8 - tcfc7)* tw(i,5) * tcfc4*abplnk1(13,i,k2)
! Subprogram not used       acfc8 = 70.0_r8*(1.0_r8 - tcfc8)* tw(i,6) * abplnk1(14,i,k2)
! Subprogram not used !
! Subprogram not used ! Emissivity for CH4 band 1306 cm-1
! Subprogram not used !
! Subprogram not used       tlw = exp(-1.0_r8*sqrt(dplh2o(i)))
! Subprogram not used       tlw=tlw*aer_trn_ttl(i,k1,k2,idx_LW_1200_2000)
! Subprogram not used       duch4 = abs(uch4(i,k1) - uch4(i,k2))
! Subprogram not used       dbetac = abs(bch4(i,k1) - bch4(i,k2))/duch4
! Subprogram not used       ach4 = 6.00444_r8*sqti(i)*log(1.0_r8 + func(duch4,dbetac))*tlw*abplnk1(3,i,k2)
! Subprogram not used       tch4 = 1.0_r8/(1.0_r8 + 0.02_r8*func(duch4,dbetac))
! Subprogram not used !
! Subprogram not used ! Absorptivity for N2O bands
! Subprogram not used !
! Subprogram not used       du01 = abs(un2o0(i,k1) - un2o0(i,k2))
! Subprogram not used       du11 = abs(un2o1(i,k1) - un2o1(i,k2))
! Subprogram not used       dbeta01 = abs(bn2o0(i,k1) - bn2o0(i,k2))/du01
! Subprogram not used       dbeta11 = abs(bn2o1(i,k1) - bn2o1(i,k2))/du11
! Subprogram not used !
! Subprogram not used ! 1285 cm-1 band
! Subprogram not used !
! Subprogram not used       an2o1 = 2.35558_r8*sqti(i)*log(1.0_r8 + func(du01,dbeta01) &
! Subprogram not used               + func(du11,dbeta11))*tlw*tch4*abplnk1(4,i,k2)
! Subprogram not used       du02 = 0.100090_r8*du01
! Subprogram not used       du12 = 0.0992746_r8*du11
! Subprogram not used       dbeta02 = 0.964282_r8*dbeta01
! Subprogram not used !
! Subprogram not used ! 589 cm-1 band
! Subprogram not used !
! Subprogram not used       an2o2 = 2.65581_r8*sqti(i)*log(1.0_r8 + func(du02,dbeta02) + &
! Subprogram not used               func(du12,dbeta02))*th2o(i)*tco2(i)*abplnk1(5,i,k2)
! Subprogram not used       du03 = 0.0333767_r8*du01
! Subprogram not used       dbeta03 = 0.982143_r8*dbeta01
! Subprogram not used !
! Subprogram not used ! 1168 cm-1 band
! Subprogram not used !
! Subprogram not used       an2o3 = 2.54034_r8*sqti(i)*log(1.0_r8 + func(du03,dbeta03))* &
! Subprogram not used               tw(i,6)*tcfc8*abplnk1(6,i,k2)
! Subprogram not used !
! Subprogram not used ! Emissivity for 1064 cm-1 band of CO2
! Subprogram not used !
! Subprogram not used       du11 = abs(uco211(i,k1) - uco211(i,k2))
! Subprogram not used       du12 = abs(uco212(i,k1) - uco212(i,k2))
! Subprogram not used       du13 = abs(uco213(i,k1) - uco213(i,k2))
! Subprogram not used       dbetc1 = 2.97558_r8*abs(pnm(i,k1) + pnm(i,k2))/(2.0_r8*sslp*sqti(i))
! Subprogram not used       dbetc2 = 2.0_r8*dbetc1
! Subprogram not used       aco21 = 3.7571_r8*sqti(i)*log(1.0_r8 + func(du11,dbetc1) &
! Subprogram not used               + func(du12,dbetc2) + func(du13,dbetc2)) &
! Subprogram not used               *to3(i)*tw(i,5)*tcfc4*tcfc7*abplnk1(2,i,k2)
! Subprogram not used !
! Subprogram not used ! Emissivity for 961 cm-1 band
! Subprogram not used !
! Subprogram not used       du21 = abs(uco221(i,k1) - uco221(i,k2))
! Subprogram not used       du22 = abs(uco222(i,k1) - uco222(i,k2))
! Subprogram not used       du23 = abs(uco223(i,k1) - uco223(i,k2))
! Subprogram not used       aco22 = 3.8443_r8*sqti(i)*log(1.0_r8 + func(du21,dbetc1) &
! Subprogram not used               + func(du22,dbetc1) + func(du23,dbetc2)) &
! Subprogram not used               *tw(i,4)*tcfc3*tcfc6*abplnk1(1,i,k2)
! Subprogram not used !
! Subprogram not used ! total trace gas absorptivity
! Subprogram not used !
! Subprogram not used       abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
! Subprogram not used                   acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
! Subprogram not used                   aco21 + aco22
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine trcab

!====================================================================================

! Subprogram not used subroutine trcabn(ncol    ,                                     &
! Subprogram not used                   k2      ,kn      ,ucfc11  ,ucfc12  ,un2o0   , &
! Subprogram not used                   un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
! Subprogram not used                   uco221  ,uco222  ,uco223  ,tbar    ,bplnk   , &
! Subprogram not used                   winpl   ,pinpl   ,tco2    ,th2o    ,to3     , &
! Subprogram not used                   uptype  ,dw      ,s2c     ,up2     ,pnew    , &
! Subprogram not used                   abstrc  ,uinpl   , &
! Subprogram not used                   aer_trn_ngh)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Calculate nearest layer absorptivity due to CH4, N2O, CFC11 and CFC12
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! Equations in CCM3 description
! Subprogram not used ! 
! Subprogram not used ! Author: J. Kiehl
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used    integer, intent(in) :: k2                    ! level index
! Subprogram not used    integer, intent(in) :: kn                    ! level index
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: tbar(pcols,4)        ! pressure weighted temperature
! Subprogram not used    real(r8), intent(in) :: ucfc11(pcols,pverp)  ! CFC11 path length
! Subprogram not used    real(r8), intent(in) :: ucfc12(pcols,pverp)  ! CFC12 path length
! Subprogram not used    real(r8), intent(in) :: un2o0(pcols,pverp)   ! N2O path length
! Subprogram not used    real(r8), intent(in) :: un2o1(pcols,pverp)   ! N2O path length (hot band)
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: uch4(pcols,pverp)    ! CH4 path length
! Subprogram not used    real(r8), intent(in) :: uco211(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco212(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco213(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco221(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: uco222(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco223(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: bplnk(14,pcols,4)    ! weighted Planck fnc. for absorptivity
! Subprogram not used    real(r8), intent(in) :: winpl(pcols,4)       ! fractional path length
! Subprogram not used    real(r8), intent(in) :: pinpl(pcols,4)       ! pressure factor for subdivided layer
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: tco2(pcols)          ! co2 transmission
! Subprogram not used    real(r8), intent(in) :: th2o(pcols)          ! h2o transmission
! Subprogram not used    real(r8), intent(in) :: to3(pcols)           ! o3 transmission
! Subprogram not used    real(r8), intent(in) :: dw(pcols)            ! h2o path length
! Subprogram not used    real(r8), intent(in) :: pnew(pcols)          ! pressure factor
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: s2c(pcols,pverp)     ! h2o continuum factor
! Subprogram not used    real(r8), intent(in) :: uptype(pcols,pverp)  ! p-type path length
! Subprogram not used    real(r8), intent(in) :: up2(pcols)           ! p squared path length
! Subprogram not used    real(r8), intent(in) :: uinpl(pcols,4)       ! Nearest layer subdivision factor
! Subprogram not used    real(r8), intent(in) :: aer_trn_ngh(pcols,nlwbands) 
! Subprogram not used                              ! [fraction] Total transmission between 
! Subprogram not used                              !            nearest neighbor sub-levels
! Subprogram not used !
! Subprogram not used !  Output Arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: abstrc(pcols)        ! total trace gas absorptivity
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !--------------------------Local Variables------------------------------
! Subprogram not used !
! Subprogram not used    integer i,l                   ! loop counters
! Subprogram not used !
! Subprogram not used    real(r8) sqti(pcols)          ! square root of mean temp
! Subprogram not used    real(r8) rsqti(pcols)         ! reciprocal of sqti
! Subprogram not used    real(r8) du1                  ! cfc11 path length
! Subprogram not used    real(r8) du2                  ! cfc12 path length
! Subprogram not used    real(r8) acfc1                ! absorptivity of cfc11 798 cm-1 band
! Subprogram not used !
! Subprogram not used    real(r8) acfc2                ! absorptivity of cfc11 846 cm-1 band
! Subprogram not used    real(r8) acfc3                ! absorptivity of cfc11 933 cm-1 band
! Subprogram not used    real(r8) acfc4                ! absorptivity of cfc11 1085 cm-1 band
! Subprogram not used    real(r8) acfc5                ! absorptivity of cfc11 889 cm-1 band
! Subprogram not used    real(r8) acfc6                ! absorptivity of cfc11 923 cm-1 band
! Subprogram not used !
! Subprogram not used    real(r8) acfc7                ! absorptivity of cfc11 1102 cm-1 band
! Subprogram not used    real(r8) acfc8                ! absorptivity of cfc11 1161 cm-1 band
! Subprogram not used    real(r8) du01                 ! n2o path length
! Subprogram not used    real(r8) dbeta01              ! n2o pressure factors
! Subprogram not used    real(r8) dbeta11              !        "
! Subprogram not used !
! Subprogram not used    real(r8)  an2o1               ! absorptivity of the 1285 cm-1 n2o band
! Subprogram not used    real(r8) du02                 ! n2o path length
! Subprogram not used    real(r8) dbeta02              ! n2o pressure factor
! Subprogram not used    real(r8) an2o2                ! absorptivity of the 589 cm-1 n2o band
! Subprogram not used    real(r8) du03                 ! n2o path length
! Subprogram not used !
! Subprogram not used    real(r8) dbeta03              ! n2o pressure factor
! Subprogram not used    real(r8) an2o3                ! absorptivity of the 1168 cm-1 n2o band
! Subprogram not used    real(r8) duch4                ! ch4 path length
! Subprogram not used    real(r8) dbetac               ! ch4 pressure factor
! Subprogram not used    real(r8) ach4                 ! absorptivity of the 1306 cm-1 ch4 band
! Subprogram not used !
! Subprogram not used    real(r8) du11                 ! co2 path length
! Subprogram not used    real(r8) du12                 !       "
! Subprogram not used    real(r8) du13                 !       "
! Subprogram not used    real(r8) dbetc1               ! co2 pressure factor
! Subprogram not used    real(r8) dbetc2               ! co2 pressure factor
! Subprogram not used !
! Subprogram not used    real(r8) aco21                ! absorptivity of the 1064 cm-1 co2 band
! Subprogram not used    real(r8) du21                 ! co2 path length
! Subprogram not used    real(r8) du22                 !       "
! Subprogram not used    real(r8) du23                 !       "
! Subprogram not used    real(r8) aco22                ! absorptivity of the 961 cm-1 co2 band
! Subprogram not used !
! Subprogram not used    real(r8) tt(pcols)            ! temp. factor for h2o overlap
! Subprogram not used    real(r8) psi1                 !          "
! Subprogram not used    real(r8) phi1                 !          "
! Subprogram not used    real(r8) p1                   ! factor for h2o overlap
! Subprogram not used    real(r8) w1                   !          "
! Subprogram not used !
! Subprogram not used    real(r8) ds2c(pcols)          ! continuum path length
! Subprogram not used    real(r8) duptyp(pcols)        ! p-type path length
! Subprogram not used    real(r8) tw(pcols,6)          ! h2o transmission overlap
! Subprogram not used !   real(r8) g1(6)                ! h2o overlap factor
! Subprogram not used !   real(r8) g2(6)                !         "
! Subprogram not used !
! Subprogram not used !   real(r8) g3(6)                !         "
! Subprogram not used !   real(r8) g4(6)                !         "
! Subprogram not used !   real(r8) ab(6)                ! h2o temp. factor
! Subprogram not used !   real(r8) bb(6)                !         "
! Subprogram not used !   real(r8) abp(6)               !         "
! Subprogram not used !
! Subprogram not used !   real(r8) bbp(6)               !         "
! Subprogram not used    real(r8) tcfc3                ! transmission of cfc11 band
! Subprogram not used    real(r8) tcfc4                ! transmission of cfc11 band
! Subprogram not used    real(r8) tcfc6                ! transmission of cfc12 band
! Subprogram not used    real(r8) tcfc7                !         "
! Subprogram not used !
! Subprogram not used    real(r8) tcfc8                !         "
! Subprogram not used    real(r8) tlw                  ! h2o transmission
! Subprogram not used    real(r8) tch4                 ! ch4 transmission
! Subprogram not used !
! Subprogram not used !--------------------------Data Statements------------------------------
! Subprogram not used !
! Subprogram not used !   data g1 /0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/
! Subprogram not used !   data g2 /14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/
! Subprogram not used !   data g3 /26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/
! Subprogram not used !   data g4 /0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/
! Subprogram not used !   data ab /3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/
! Subprogram not used !   data bb /-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/
! Subprogram not used !   data abp/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/
! Subprogram not used !   data bbp/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/
! Subprogram not used !
! Subprogram not used !--------------------------Statement Functions--------------------------
! Subprogram not used !
! Subprogram not used    real(r8) func, u, b
! Subprogram not used    func(u,b) = u/sqrt(4.0_r8 + u*(1.0_r8 + 1.0_r8 / b))
! Subprogram not used !
! Subprogram not used !------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    do i = 1,ncol
! Subprogram not used       sqti(i) = sqrt(tbar(i,kn))
! Subprogram not used       rsqti(i) = 1._r8 / sqti(i)
! Subprogram not used !
! Subprogram not used ! h2o transmission
! Subprogram not used !
! Subprogram not used       tt(i) = abs(tbar(i,kn) - 250.0_r8)
! Subprogram not used       ds2c(i) = abs(s2c(i,k2+1) - s2c(i,k2))*uinpl(i,kn)
! Subprogram not used       duptyp(i) = abs(uptype(i,k2+1) - uptype(i,k2))*uinpl(i,kn)
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do l = 1,6
! Subprogram not used       do i = 1,ncol
! Subprogram not used          psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
! Subprogram not used          phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
! Subprogram not used          p1 = pnew(i) * (psi1/phi1) / sslp
! Subprogram not used          w1 = dw(i) * winpl(i,kn) * phi1
! Subprogram not used          tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0_r8+g2(l)*(w1/p1))-1.0_r8) &
! Subprogram not used                    - g3(l)*ds2c(i)-g4(l)*duptyp(i))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do i=1,ncol
! Subprogram not used       tw(i,1)=tw(i,1)*(0.7_r8*aer_trn_ngh(i,idx_LW_0650_0800)+&! l=1: 0750--0820 cm-1
! Subprogram not used                        0.3_r8*aer_trn_ngh(i,idx_LW_0800_1000))
! Subprogram not used       tw(i,2)=tw(i,2)*aer_trn_ngh(i,idx_LW_0800_1000) ! l=2: 0820--0880 cm-1
! Subprogram not used       tw(i,3)=tw(i,3)*aer_trn_ngh(i,idx_LW_0800_1000) ! l=3: 0880--0900 cm-1
! Subprogram not used       tw(i,4)=tw(i,4)*aer_trn_ngh(i,idx_LW_0800_1000) ! l=4: 0900--1000 cm-1
! Subprogram not used       tw(i,5)=tw(i,5)*aer_trn_ngh(i,idx_LW_1000_1200) ! l=5: 1000--1120 cm-1
! Subprogram not used       tw(i,6)=tw(i,6)*aer_trn_ngh(i,idx_LW_1000_1200) ! l=6: 1120--1170 cm-1
! Subprogram not used    end do                    ! end loop over lon
! Subprogram not used 
! Subprogram not used    do i = 1,ncol
! Subprogram not used !
! Subprogram not used       du1 = abs(ucfc11(i,k2+1) - ucfc11(i,k2)) * winpl(i,kn)
! Subprogram not used       du2 = abs(ucfc12(i,k2+1) - ucfc12(i,k2)) * winpl(i,kn)
! Subprogram not used !
! Subprogram not used ! cfc transmissions
! Subprogram not used !
! Subprogram not used       tcfc3 = exp(-175.005_r8*du1)
! Subprogram not used       tcfc4 = exp(-1202.18_r8*du1)
! Subprogram not used       tcfc6 = exp(-5786.73_r8*du2)
! Subprogram not used       tcfc7 = exp(-2873.51_r8*du2)
! Subprogram not used       tcfc8 = exp(-2085.59_r8*du2)
! Subprogram not used !
! Subprogram not used ! Absorptivity for CFC11 bands
! Subprogram not used !
! Subprogram not used       acfc1 = 50.0_r8*(1.0_r8 - exp(-54.09_r8*du1)) * tw(i,1)*bplnk(7,i,kn)
! Subprogram not used       acfc2 = 60.0_r8*(1.0_r8 - exp(-5130.03_r8*du1))*tw(i,2)*bplnk(8,i,kn)
! Subprogram not used       acfc3 = 60.0_r8*(1.0_r8 - tcfc3)*tw(i,4)*tcfc6 * bplnk(9,i,kn)
! Subprogram not used       acfc4 = 100.0_r8*(1.0_r8 - tcfc4)* tw(i,5) * bplnk(10,i,kn)
! Subprogram not used !
! Subprogram not used ! Absorptivity for CFC12 bands
! Subprogram not used !
! Subprogram not used       acfc5 = 45.0_r8*(1.0_r8 - exp(-1272.35_r8*du2))*tw(i,3)*bplnk(11,i,kn)
! Subprogram not used       acfc6 = 50.0_r8*(1.0_r8 - tcfc6)*tw(i,4)*bplnk(12,i,kn)
! Subprogram not used       acfc7 = 80.0_r8*(1.0_r8 - tcfc7)* tw(i,5)*tcfc4 *bplnk(13,i,kn)
! Subprogram not used       acfc8 = 70.0_r8*(1.0_r8 - tcfc8)*tw(i,6)*bplnk(14,i,kn)
! Subprogram not used !
! Subprogram not used ! Absorptivity for CH4 band 1306 cm-1
! Subprogram not used !
! Subprogram not used       tlw = exp(-1.0_r8*sqrt(up2(i)))
! Subprogram not used       tlw=tlw*aer_trn_ngh(i,idx_LW_1200_2000)
! Subprogram not used       duch4 = abs(uch4(i,k2+1) - uch4(i,k2)) * winpl(i,kn)
! Subprogram not used       dbetac = 2.94449_r8 * pinpl(i,kn) * rsqti(i) / sslp
! Subprogram not used       ach4 = 6.00444_r8*sqti(i)*log(1.0_r8 + func(duch4,dbetac)) * tlw * bplnk(3,i,kn)
! Subprogram not used       tch4 = 1.0_r8/(1.0_r8 + 0.02_r8*func(duch4,dbetac))
! Subprogram not used !
! Subprogram not used ! Absorptivity for N2O bands
! Subprogram not used !
! Subprogram not used       du01 = abs(un2o0(i,k2+1) - un2o0(i,k2)) * winpl(i,kn)
! Subprogram not used       du11 = abs(un2o1(i,k2+1) - un2o1(i,k2)) * winpl(i,kn)
! Subprogram not used       dbeta01 = 19.399_r8 *  pinpl(i,kn) * rsqti(i) / sslp
! Subprogram not used       dbeta11 = dbeta01
! Subprogram not used !
! Subprogram not used ! 1285 cm-1 band
! Subprogram not used !
! Subprogram not used       an2o1 = 2.35558_r8*sqti(i)*log(1.0_r8 + func(du01,dbeta01) &
! Subprogram not used               + func(du11,dbeta11)) * tlw * tch4 * bplnk(4,i,kn)
! Subprogram not used       du02 = 0.100090_r8*du01
! Subprogram not used       du12 = 0.0992746_r8*du11
! Subprogram not used       dbeta02 = 0.964282_r8*dbeta01
! Subprogram not used !
! Subprogram not used ! 589 cm-1 band
! Subprogram not used !
! Subprogram not used       an2o2 = 2.65581_r8*sqti(i)*log(1.0_r8 + func(du02,dbeta02) &
! Subprogram not used               +  func(du12,dbeta02)) * tco2(i) * th2o(i) * bplnk(5,i,kn)
! Subprogram not used       du03 = 0.0333767_r8*du01
! Subprogram not used       dbeta03 = 0.982143_r8*dbeta01
! Subprogram not used !
! Subprogram not used ! 1168 cm-1 band
! Subprogram not used !
! Subprogram not used       an2o3 = 2.54034_r8*sqti(i)*log(1.0_r8 + func(du03,dbeta03)) * &
! Subprogram not used               tw(i,6) * tcfc8 * bplnk(6,i,kn)
! Subprogram not used !
! Subprogram not used ! Absorptivity for 1064 cm-1 band of CO2
! Subprogram not used !
! Subprogram not used       du11 = abs(uco211(i,k2+1) - uco211(i,k2)) * winpl(i,kn)
! Subprogram not used       du12 = abs(uco212(i,k2+1) - uco212(i,k2)) * winpl(i,kn)
! Subprogram not used       du13 = abs(uco213(i,k2+1) - uco213(i,k2)) * winpl(i,kn)
! Subprogram not used       dbetc1 = 2.97558_r8 * pinpl(i,kn) * rsqti(i) / sslp
! Subprogram not used       dbetc2 = 2.0_r8 * dbetc1
! Subprogram not used       aco21 = 3.7571_r8*sqti(i)*log(1.0_r8 + func(du11,dbetc1) &
! Subprogram not used               + func(du12,dbetc2) + func(du13,dbetc2)) &
! Subprogram not used               * to3(i) * tw(i,5) * tcfc4 * tcfc7 * bplnk(2,i,kn)
! Subprogram not used !
! Subprogram not used ! Absorptivity for 961 cm-1 band of co2
! Subprogram not used !
! Subprogram not used       du21 = abs(uco221(i,k2+1) - uco221(i,k2)) * winpl(i,kn)
! Subprogram not used       du22 = abs(uco222(i,k2+1) - uco222(i,k2)) * winpl(i,kn)
! Subprogram not used       du23 = abs(uco223(i,k2+1) - uco223(i,k2)) * winpl(i,kn)
! Subprogram not used       aco22 = 3.8443_r8*sqti(i)*log(1.0_r8 + func(du21,dbetc1) &
! Subprogram not used               + func(du22,dbetc1) + func(du23,dbetc2)) &
! Subprogram not used               * tw(i,4) * tcfc3 * tcfc6 * bplnk(1,i,kn)
! Subprogram not used !
! Subprogram not used ! total trace gas absorptivity
! Subprogram not used !
! Subprogram not used       abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
! Subprogram not used                   acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
! Subprogram not used                   aco21 + aco22
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine trcabn

!====================================================================================

! Subprogram not used subroutine trcems(ncol    ,                                     &
! Subprogram not used                   k       ,co2t    ,pnm     ,ucfc11  ,ucfc12  , &
! Subprogram not used                   un2o0   ,un2o1   ,bn2o0   ,bn2o1   ,uch4    , &
! Subprogram not used                   bch4    ,uco211  ,uco212  ,uco213  ,uco221  , &
! Subprogram not used                   uco222  ,uco223  ,uptype  ,w       ,s2c     , &
! Subprogram not used                   up2     ,emplnk  ,th2o    ,tco2    ,to3     , &
! Subprogram not used                   emstrc  , &
! Subprogram not used                   aer_trn_ttl)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used !  Calculate emissivity for CH4, N2O, CFC11 and CFC12 bands.
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used !  See CCM3 Description for equations.
! Subprogram not used ! 
! Subprogram not used ! Author: J. Kiehl
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: co2t(pcols,pverp)    ! pressure weighted temperature
! Subprogram not used    real(r8), intent(in) :: pnm(pcols,pverp)     ! interface pressure
! Subprogram not used    real(r8), intent(in) :: ucfc11(pcols,pverp)  ! CFC11 path length
! Subprogram not used    real(r8), intent(in) :: ucfc12(pcols,pverp)  ! CFC12 path length
! Subprogram not used    real(r8), intent(in) :: un2o0(pcols,pverp)   ! N2O path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: un2o1(pcols,pverp)   ! N2O path length (hot band)
! Subprogram not used    real(r8), intent(in) :: uch4(pcols,pverp)    ! CH4 path length
! Subprogram not used    real(r8), intent(in) :: uco211(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco212(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco213(pcols,pverp)  ! CO2 9.4 micron band path length
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: uco221(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco222(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uco223(pcols,pverp)  ! CO2 10.4 micron band path length
! Subprogram not used    real(r8), intent(in) :: uptype(pcols,pverp)  ! continuum path length
! Subprogram not used    real(r8), intent(in) :: bn2o0(pcols,pverp)   ! pressure factor for n2o
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: bn2o1(pcols,pverp)   ! pressure factor for n2o
! Subprogram not used    real(r8), intent(in) :: bch4(pcols,pverp)    ! pressure factor for ch4
! Subprogram not used    real(r8), intent(in) :: emplnk(14,pcols)     ! emissivity Planck factor
! Subprogram not used    real(r8), intent(in) :: th2o(pcols)          ! water vapor overlap factor
! Subprogram not used    real(r8), intent(in) :: tco2(pcols)          ! co2 overlap factor
! Subprogram not used !
! Subprogram not used    real(r8), intent(in) :: to3(pcols)           ! o3 overlap factor
! Subprogram not used    real(r8), intent(in) :: s2c(pcols,pverp)     ! h2o continuum path length
! Subprogram not used    real(r8), intent(in) :: w(pcols,pverp)       ! h2o path length
! Subprogram not used    real(r8), intent(in) :: up2(pcols)           ! pressure squared h2o path length
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: k                 ! level index
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) ! aer trn.
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !  Output Arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: emstrc(pcols,pverp)  ! total trace gas emissivity
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !--------------------------Local Variables------------------------------
! Subprogram not used !
! Subprogram not used    integer i,l               ! loop counters
! Subprogram not used !
! Subprogram not used    real(r8) sqti(pcols)          ! square root of mean temp
! Subprogram not used    real(r8) ecfc1                ! emissivity of cfc11 798 cm-1 band
! Subprogram not used    real(r8) ecfc2                !     "      "    "   846 cm-1 band
! Subprogram not used    real(r8) ecfc3                !     "      "    "   933 cm-1 band
! Subprogram not used    real(r8) ecfc4                !     "      "    "   1085 cm-1 band
! Subprogram not used !
! Subprogram not used    real(r8) ecfc5                !     "      "  cfc12 889 cm-1 band
! Subprogram not used    real(r8) ecfc6                !     "      "    "   923 cm-1 band
! Subprogram not used    real(r8) ecfc7                !     "      "    "   1102 cm-1 band
! Subprogram not used    real(r8) ecfc8                !     "      "    "   1161 cm-1 band
! Subprogram not used    real(r8) u01                  ! n2o path length
! Subprogram not used !
! Subprogram not used    real(r8) u11                  ! n2o path length
! Subprogram not used    real(r8) beta01               ! n2o pressure factor
! Subprogram not used    real(r8) beta11               ! n2o pressure factor
! Subprogram not used    real(r8) en2o1                ! emissivity of the 1285 cm-1 N2O band
! Subprogram not used    real(r8) u02                  ! n2o path length
! Subprogram not used !
! Subprogram not used    real(r8) u12                  ! n2o path length
! Subprogram not used    real(r8) beta02               ! n2o pressure factor
! Subprogram not used    real(r8) en2o2                ! emissivity of the 589 cm-1 N2O band
! Subprogram not used    real(r8) u03                  ! n2o path length
! Subprogram not used    real(r8) beta03               ! n2o pressure factor
! Subprogram not used !
! Subprogram not used    real(r8) en2o3                ! emissivity of the 1168 cm-1 N2O band
! Subprogram not used    real(r8) betac                ! ch4 pressure factor
! Subprogram not used    real(r8) ech4                 ! emissivity of 1306 cm-1 CH4 band
! Subprogram not used    real(r8) betac1               ! co2 pressure factor
! Subprogram not used    real(r8) betac2               ! co2 pressure factor
! Subprogram not used !
! Subprogram not used    real(r8) eco21                ! emissivity of 1064 cm-1 CO2 band
! Subprogram not used    real(r8) eco22                ! emissivity of 961 cm-1 CO2 band
! Subprogram not used    real(r8) tt(pcols)            ! temp. factor for h2o overlap factor
! Subprogram not used    real(r8) psi1                 ! narrow band h2o temp. factor
! Subprogram not used    real(r8) phi1                 !             "
! Subprogram not used !
! Subprogram not used    real(r8) p1                   ! h2o line overlap factor
! Subprogram not used    real(r8) w1                   !          "
! Subprogram not used    real(r8) tw(pcols,6)          ! h2o transmission overlap
! Subprogram not used !   real(r8) g1(6)                ! h2o overlap factor
! Subprogram not used !   real(r8) g2(6)                !          "
! Subprogram not used !
! Subprogram not used !   real(r8) g3(6)                !          "
! Subprogram not used !   real(r8) g4(6)                !          "
! Subprogram not used !   real(r8) ab(6)                !          "
! Subprogram not used !   real(r8) bb(6)                !          "
! Subprogram not used !   real(r8) abp(6)               !          "
! Subprogram not used !
! Subprogram not used !   real(r8) bbp(6)               !          "
! Subprogram not used    real(r8) tcfc3                ! transmission for cfc11 band
! Subprogram not used    real(r8) tcfc4                !          "
! Subprogram not used    real(r8) tcfc6                ! transmission for cfc12 band
! Subprogram not used    real(r8) tcfc7                !          "
! Subprogram not used !
! Subprogram not used    real(r8) tcfc8                !          "
! Subprogram not used    real(r8) tlw                  ! h2o overlap factor
! Subprogram not used    real(r8) tch4                 ! ch4 overlap factor
! Subprogram not used !
! Subprogram not used !--------------------------Data Statements------------------------------
! Subprogram not used !
! Subprogram not used !   data g1 /0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/
! Subprogram not used !   data g2 /14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/
! Subprogram not used !   data g3 /26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/
! Subprogram not used !   data g4 /0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/
! Subprogram not used !   data ab /3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/
! Subprogram not used !   data bb /-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/
! Subprogram not used !   data abp/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/
! Subprogram not used !   data bbp/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/
! Subprogram not used !
! Subprogram not used !--------------------------Statement Functions--------------------------
! Subprogram not used !
! Subprogram not used    real(r8) func, u, b
! Subprogram not used    func(u,b) = u/sqrt(4.0_r8 + u*(1.0_r8 + 1.0_r8 / b))
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    do i = 1,ncol
! Subprogram not used       sqti(i) = sqrt(co2t(i,k))
! Subprogram not used !
! Subprogram not used ! Transmission for h2o
! Subprogram not used !
! Subprogram not used       tt(i) = abs(co2t(i,k) - 250.0_r8)
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used    do l = 1,6
! Subprogram not used       do i = 1,ncol
! Subprogram not used          psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
! Subprogram not used          phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
! Subprogram not used          p1 = pnm(i,k) * (psi1/phi1) / sslp
! Subprogram not used          w1 = w(i,k) * phi1
! Subprogram not used          tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0_r8+g2(l)*(w1/p1))-1.0_r8) &
! Subprogram not used                    - g3(l)*s2c(i,k)-g4(l)*uptype(i,k))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used !     Overlap H2O tranmission with STRAER continuum in 6 trace gas 
! Subprogram not used !                 subbands
! Subprogram not used 
! Subprogram not used       do i=1,ncol
! Subprogram not used          tw(i,1)=tw(i,1)*(0.7_r8*aer_trn_ttl(i,k,1,idx_LW_0650_0800)+&! l=1: 0750--0820 cm-1
! Subprogram not used                           0.3_r8*aer_trn_ttl(i,k,1,idx_LW_0800_1000))
! Subprogram not used          tw(i,2)=tw(i,2)*aer_trn_ttl(i,k,1,idx_LW_0800_1000) ! l=2: 0820--0880 cm-1
! Subprogram not used          tw(i,3)=tw(i,3)*aer_trn_ttl(i,k,1,idx_LW_0800_1000) ! l=3: 0880--0900 cm-1
! Subprogram not used          tw(i,4)=tw(i,4)*aer_trn_ttl(i,k,1,idx_LW_0800_1000) ! l=4: 0900--1000 cm-1
! Subprogram not used          tw(i,5)=tw(i,5)*aer_trn_ttl(i,k,1,idx_LW_1000_1200) ! l=5: 1000--1120 cm-1
! Subprogram not used          tw(i,6)=tw(i,6)*aer_trn_ttl(i,k,1,idx_LW_1000_1200) ! l=6: 1120--1170 cm-1
! Subprogram not used       end do                    ! end loop over lon
! Subprogram not used !
! Subprogram not used    do i = 1,ncol
! Subprogram not used !
! Subprogram not used ! transmission due to cfc bands
! Subprogram not used !
! Subprogram not used       tcfc3 = exp(-175.005_r8*ucfc11(i,k))
! Subprogram not used       tcfc4 = exp(-1202.18_r8*ucfc11(i,k))
! Subprogram not used       tcfc6 = exp(-5786.73_r8*ucfc12(i,k))
! Subprogram not used       tcfc7 = exp(-2873.51_r8*ucfc12(i,k))
! Subprogram not used       tcfc8 = exp(-2085.59_r8*ucfc12(i,k))
! Subprogram not used !
! Subprogram not used ! Emissivity for CFC11 bands
! Subprogram not used !
! Subprogram not used       ecfc1 = 50.0_r8*(1.0_r8 - exp(-54.09_r8*ucfc11(i,k))) * tw(i,1) * emplnk(7,i)
! Subprogram not used       ecfc2 = 60.0_r8*(1.0_r8 - exp(-5130.03_r8*ucfc11(i,k)))* tw(i,2) * emplnk(8,i)
! Subprogram not used       ecfc3 = 60.0_r8*(1.0_r8 - tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
! Subprogram not used       ecfc4 = 100.0_r8*(1.0_r8 - tcfc4)*tw(i,5)*emplnk(10,i)
! Subprogram not used !
! Subprogram not used ! Emissivity for CFC12 bands
! Subprogram not used !
! Subprogram not used       ecfc5 = 45.0_r8*(1.0_r8 - exp(-1272.35_r8*ucfc12(i,k)))*tw(i,3)*emplnk(11,i)
! Subprogram not used       ecfc6 = 50.0_r8*(1.0_r8 - tcfc6)*tw(i,4)*emplnk(12,i)
! Subprogram not used       ecfc7 = 80.0_r8*(1.0_r8 - tcfc7)*tw(i,5)* tcfc4 * emplnk(13,i)
! Subprogram not used       ecfc8 = 70.0_r8*(1.0_r8 - tcfc8)*tw(i,6) * emplnk(14,i)
! Subprogram not used !
! Subprogram not used ! Emissivity for CH4 band 1306 cm-1
! Subprogram not used !
! Subprogram not used       tlw = exp(-1.0_r8*sqrt(up2(i)))
! Subprogram not used 
! Subprogram not used !     Overlap H2O vibration rotation band with STRAER continuum 
! Subprogram not used !             for CH4 1306 cm-1 and N2O 1285 cm-1 bands
! Subprogram not used 
! Subprogram not used             tlw=tlw*aer_trn_ttl(i,k,1,idx_LW_1200_2000)
! Subprogram not used       betac = bch4(i,k)/uch4(i,k)
! Subprogram not used       ech4 = 6.00444_r8*sqti(i)*log(1.0_r8 + func(uch4(i,k),betac)) *tlw * emplnk(3,i)
! Subprogram not used       tch4 = 1.0_r8/(1.0_r8 + 0.02_r8*func(uch4(i,k),betac))
! Subprogram not used !
! Subprogram not used ! Emissivity for N2O bands
! Subprogram not used !
! Subprogram not used       u01 = un2o0(i,k)
! Subprogram not used       u11 = un2o1(i,k)
! Subprogram not used       beta01 = bn2o0(i,k)/un2o0(i,k)
! Subprogram not used       beta11 = bn2o1(i,k)/un2o1(i,k)
! Subprogram not used !
! Subprogram not used ! 1285 cm-1 band
! Subprogram not used !
! Subprogram not used       en2o1 = 2.35558_r8*sqti(i)*log(1.0_r8 + func(u01,beta01) + &
! Subprogram not used               func(u11,beta11))*tlw*tch4*emplnk(4,i)
! Subprogram not used       u02 = 0.100090_r8*u01
! Subprogram not used       u12 = 0.0992746_r8*u11
! Subprogram not used       beta02 = 0.964282_r8*beta01
! Subprogram not used !
! Subprogram not used ! 589 cm-1 band
! Subprogram not used !
! Subprogram not used       en2o2 = 2.65581_r8*sqti(i)*log(1.0_r8 + func(u02,beta02) + &
! Subprogram not used               func(u12,beta02)) * tco2(i) * th2o(i) * emplnk(5,i)
! Subprogram not used       u03 = 0.0333767_r8*u01
! Subprogram not used       beta03 = 0.982143_r8*beta01
! Subprogram not used !
! Subprogram not used ! 1168 cm-1 band
! Subprogram not used !
! Subprogram not used       en2o3 = 2.54034_r8*sqti(i)*log(1.0_r8 + func(u03,beta03)) * &
! Subprogram not used               tw(i,6) * tcfc8 * emplnk(6,i)
! Subprogram not used !
! Subprogram not used ! Emissivity for 1064 cm-1 band of CO2
! Subprogram not used !
! Subprogram not used       betac1 = 2.97558_r8*pnm(i,k) / (sslp*sqti(i))
! Subprogram not used       betac2 = 2.0_r8 * betac1
! Subprogram not used       eco21 = 3.7571_r8*sqti(i)*log(1.0_r8 + func(uco211(i,k),betac1) &
! Subprogram not used               + func(uco212(i,k),betac2) + func(uco213(i,k),betac2)) &
! Subprogram not used               * to3(i) * tw(i,5) * tcfc4 * tcfc7 * emplnk(2,i)
! Subprogram not used !
! Subprogram not used ! Emissivity for 961 cm-1 band
! Subprogram not used !
! Subprogram not used       eco22 = 3.8443_r8*sqti(i)*log(1.0_r8 + func(uco221(i,k),betac1) &
! Subprogram not used               + func(uco222(i,k),betac1) + func(uco223(i,k),betac2)) &
! Subprogram not used               * tw(i,4) * tcfc3 * tcfc6 * emplnk(1,i)
! Subprogram not used !
! Subprogram not used ! total trace gas emissivity
! Subprogram not used !
! Subprogram not used       emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 +ecfc6 + &
! Subprogram not used                     ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 + &
! Subprogram not used                     eco21 + eco22
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine trcems

!====================================================================================

! Subprogram not used subroutine trcplk(ncol    ,                                     &
! Subprogram not used                   tint    ,tlayr   ,tplnke  ,emplnk  ,abplnk1 , &
! Subprogram not used                   abplnk2 )
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used !   Calculate Planck factors for absorptivity and emissivity of
! Subprogram not used !   CH4, N2O, CFC11 and CFC12
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used !   Planck function and derivative evaluated at the band center.
! Subprogram not used ! 
! Subprogram not used ! Author: J. Kiehl
! Subprogram not used ! 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: ncol                 ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: tint(pcols,pverp)   ! interface temperatures
! Subprogram not used    real(r8), intent(in) :: tlayr(pcols,pverp)  ! k-1 level temperatures
! Subprogram not used    real(r8), intent(in) :: tplnke(pcols)       ! Top Layer temperature
! Subprogram not used !
! Subprogram not used ! output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: emplnk(14,pcols)         ! emissivity Planck factor
! Subprogram not used    real(r8), intent(out) :: abplnk1(14,pcols,pverp)  ! non-nearest layer Plack factor
! Subprogram not used    real(r8), intent(out) :: abplnk2(14,pcols,pverp)  ! nearest layer factor
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !--------------------------Local Variables------------------------------
! Subprogram not used !
! Subprogram not used    integer wvl                   ! wavelength index
! Subprogram not used    integer i,k                   ! loop counters
! Subprogram not used !
! Subprogram not used    real(r8) f1(14)                   ! Planck function factor
! Subprogram not used    real(r8) f2(14)                   !        "
! Subprogram not used    real(r8) f3(14)                   !        "
! Subprogram not used !
! Subprogram not used !--------------------------Data Statements------------------------------
! Subprogram not used !
! Subprogram not used    data f1 /5.85713e8_r8,7.94950e8_r8,1.47009e9_r8,1.40031e9_r8,1.34853e8_r8, &
! Subprogram not used             1.05158e9_r8,3.35370e8_r8,3.99601e8_r8,5.35994e8_r8,8.42955e8_r8, &
! Subprogram not used             4.63682e8_r8,5.18944e8_r8,8.83202e8_r8,1.03279e9_r8/
! Subprogram not used    data f2 /2.02493e11_r8,3.04286e11_r8,6.90698e11_r8,6.47333e11_r8, &
! Subprogram not used             2.85744e10_r8,4.41862e11_r8,9.62780e10_r8,1.21618e11_r8, &
! Subprogram not used             1.79905e11_r8,3.29029e11_r8,1.48294e11_r8,1.72315e11_r8, &
! Subprogram not used             3.50140e11_r8,4.31364e11_r8/
! Subprogram not used    data f3 /1383.0_r8,1531.0_r8,1879.0_r8,1849.0_r8,848.0_r8,1681.0_r8, &
! Subprogram not used             1148.0_r8,1217.0_r8,1343.0_r8,1561.0_r8,1279.0_r8,1328.0_r8, &
! Subprogram not used             1586.0_r8,1671.0_r8/
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Calculate emissivity Planck factor
! Subprogram not used !
! Subprogram not used    do wvl = 1,14
! Subprogram not used       do i = 1,ncol
! Subprogram not used          emplnk(wvl,i) = f1(wvl)/(tplnke(i)**4.0_r8*(exp(f3(wvl)/tplnke(i))-1.0_r8))
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used !
! Subprogram not used ! Calculate absorptivity Planck factor for tint and tlayr temperatures
! Subprogram not used !
! Subprogram not used    do wvl = 1,14
! Subprogram not used       do k = ntoplw, pverp
! Subprogram not used          do i = 1, ncol
! Subprogram not used !
! Subprogram not used ! non-nearlest layer function
! Subprogram not used !
! Subprogram not used             abplnk1(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tint(i,k)))  &
! Subprogram not used                                /(tint(i,k)**5.0_r8*(exp(f3(wvl)/tint(i,k))-1.0_r8)**2.0_r8)
! Subprogram not used !
! Subprogram not used ! nearest layer function
! Subprogram not used !
! Subprogram not used             abplnk2(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tlayr(i,k))) &
! Subprogram not used                                /(tlayr(i,k)**5.0_r8*(exp(f3(wvl)/tlayr(i,k))-1.0_r8)**2.0_r8)
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine trcplk


!====================================================================================

end module radae
