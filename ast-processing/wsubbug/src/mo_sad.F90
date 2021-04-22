      module mo_sad

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use physconst,     only : pi
      use ppgrid,        only : pcols, pver
      use m_sad_data,    only : a, b
      use cam_logfile,   only : iulog
      use spmd_utils,    only : masterproc
      

      implicit none

      private
      public  :: sad_inti
      public  :: sad_strat_calc
      public  :: sad_top

      save

      real(r8), parameter :: four_thrd = 4._r8/3._r8
      real(r8), parameter :: one_thrd  = 1._r8/3._r8
      real(r8), parameter :: two_thrd  = 2._r8/3._r8
      real(r8), parameter :: four_pi   = 4._r8*pi

      integer :: sad_top
      integer :: sad_topp

    contains
        

      subroutine sad_inti(pbuf2d)
!----------------------------------------------------------------------
!     ... initialize the sad module
!----------------------------------------------------------------------

      use time_manager, only : is_first_step
      use ref_pres,     only : pref_mid_norm
      use cam_history,  only : addfld, phys_decomp
      use physics_buffer, only : physics_buffer_desc, pbuf_set_field

      implicit none

      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

!---------------------------------------------------------------------- 
!	... Local variables
!---------------------------------------------------------------------- 
      integer  ::  k

!---------------------------------------------------------------------- 
!	... find level where etamids are all > 1 hPa
!---------------------------------------------------------------------- 
      sad_top = 0
      do k = pver,1,-1
	 if( (pref_mid_norm(k)) < .001_r8 ) then
             sad_top = k
             exit
         end if
      end do
      sad_topp = sad_top + 1
      if (masterproc) then
         write(iulog,*) 'sad_inti: sad capped at level ',sad_top
         write(iulog,*) '          whose midpoint is ',pref_mid_norm(sad_topp)*1.e3_r8,' hPa'
      endif

      call addfld( 'H2SO4M_C',  'ug/m3 ',   pver, 'I', 'chemical sulfate aerosol mass',  phys_decomp )

      end subroutine sad_inti
!===============================================================================
! ROUTINE
!   sad_strat_calc
!
!   Date...
!     14 October 1999
!
!   Programmed by...
!     Douglas E. Kinnison
!   Modified by
!     Stacy Walters
!     1 November 1999
!   Modified by 
!     Doug Kinnison
!     1 September 2004; Condensed phase H2O passed in from 1
!     2 November  2004; New treatment of denitrificatoin (NAT)
!    14 November  2004; STS mode of operation.
!    27 March     2008; Using original NAT approach.
!    08 November  2010; STS Approach (HNO3 => STS; then HNO3 => NAT) 
!    24 March     2011; updated mask logic and removed sm NAT logic
!    14 April     2011; updated EQUIL logic
!    19 December  2012; updated using Wegner et al., JGR, 2013a,b.
!    25 April     2013; Removed volcanic heating logic.
!
! DESCRIPTION
!
!     This routine has the logic to derive the surface area density for
!     three types of aerosols: Sulfate (LBS, STS); Nitric Acid Trihydrate (NAT);
!     and Water-ICE. The surface area density is stored in sad_strat(3). The
!     first, second, and third dimensions are SULFATE, NAT, and ICE SAD
!     respectively. The effective radius of each particle is also stored
!     in radius_strat(3).
!
!     NOTE1: For Sulfate and H2O ICE calculations
!     The Surface Area and Volume Densities are calculated from the
!     second and third moment of the LOG Normal Distribution. For an example
!     see Binkowski et al., JGR, 100, 26191-26209, 1995. The Volume Density
!     is substituted into the SAD equation so that the SAD is dependent on
!     the # of particles cm-3, the width of the distribution (sigma), and
!     the volume density of the aerosol. This approach is discussed in 
!     Considine et al., 2000 and Kinnison et al., 2007.
!
!     NOTE2: For the ternary solution calculation
!     the total sulfate mass is derived from the SAGEII SAD data. This approach
!     has been previously used in Considine et al., JGR, 1999. The thermodynamic
!     models used in this routine are from A. Tabazedeh et al, 1994.
!
!     NOTE3: Updates to the PSC scheme is discussed in Wegner et al., 2013a. 
!     80% of the total HNO3 is allowed to see STS, 20% NAT. The number density of 
!     the NAT and ICE particles are is set to 0.01 and 0.1 particle cm-3 respectively. 
!
!     NOTE4: The HCl solubility (in STS) has been added and evalutede in Wegner et al., 2013b. 
!     This solubility is based on Carslaw et al., 1995.
!
!     REFERENCES for this PSC module:
!        Considine, D. B., A. R. Douglass, P. S. Connell, D. E. Kinnison, and D. A., Rotman, 
!          A polar stratospheric cloud parameterization for the three dimensional model of 
!          the global modeling initiative and its response to stratospheric aircraft, 
!          J. Geophys. Res., 105, 3955-3975, 2000.
!
!        Kinnison, D. E.,et al., Sensitivity of chemical tracers to meteorological 
!          parameters in the MOZART-3 chemical transport model, J. Geophys. Res., 
!          112, D20302, doi:10.1029/2006JD007879, 2007.
!
!        Wegner, T, D. E. Kinnison, R. R. Garcia, S. Madronich, S. Solomon, and M. von Hobe, 
!          On the depletion of HCl in the Antarctic polar vortex, 
!          in review J. Geophys. Res., 2013.
!
!        Wegner, T, D. E. Kinnison, R. R. Garcia, S. Madronich, and S. Solomon, 
!          Polar Stratospheric Clouds in SD-WACCM4, in review J. Geophys. Res., 2013.
!
!        Tabazedeh, A., R. P. Turco, K. Drdla, M. Z. Jacobson, and O. B, Toon, A study
!          of the type I polar stratosphere cloud formation, 
!          Geophys. Res. Lett., 21, 1619-1622, 1994.
!
!        Carslaw, K. S., S. L. Clegg, and P. Brimblecombe, A thermodynamic model of the
!          system HCl-HNO3-H2SO4-H2O, including solubilities of HBr, from <200 to 328K, 
!          J. Phys. Chem., 99, 11,557-11,574, doi:1021/100029a039, 1995.
!
!
! ARGUMENTS
!   INPUT:
!     hno3_gas		Nitric Acid gas   phase abundance (mole fraction)
!     hno3_cond(2)	Nitric Acid cond. phase abundance (mole fraction)
!                       (1) in STS; (2) in NAT
!     h2o_cond		Water condensed phase           (mole fraction)
!     h2o_gas           Water gas-phase abundance       (mole fraction)
!
!     hcl_gas           HCL gas-phase abundance         (mole fraction)
!     hcl_cond          HCl condensed phase (STS)       (mole fraction)
!
!     sage_sad		SAGEII surface area density     (cm2-aer cm-3 atm)
!     m                 Airdensity                      (molecules cm-3)
!     press             Pressure                        (hPa)
!
!
!   OUTPUT:
!
!     hno3_gas     = Gas-phase HNO3             Used in chemical solver. 
!     hno3_cond(1) = Condensed HNO3 from STS    Not used in mo_aero_settling.F90
!     hno3_cond(2) = Condensed HNO3 from NAT    Used in mo_aero_settling.F90
!
!     hcl_gas      = Gas-phase HCL              Used in chemical solver. 
!     hcl_cond     = Condensed HCl from STS     
!
!     SAD_strat(1) = Sulfate Aerosol... Used in mo_strato_rates.F90
!     SAD_strat(2) = NAT Aerosol....    Used in mo_strato_rates.F90
!     SAD_strat(3) = Water-Ice......... Used in mo_strato_rates.F90
!
!     RAD_strat(1) = Sulfate Aerosol... Used in mo_strato_rates.F90
!     RAD_strat(2) = NAT large mode.... Used in mo_aero_settling.F90
!     RAD_strat(3) = Water-Ice......... Not used in mo_aero_settling.F90
!
!   NOTE1: The sum of hno3_cond(1-2) will be added to hno3_gas for advection of HNO3 in
!          mo_gas_phase_chemdr.F90.
!
!   NOTE2: The sum of hcl_cond will be added to hcl_gas for advection of HCl in
!          mo_gas_phase_chemdr.F90.
!
!   NOTE3: This routine does not partition H2O.
!
!
! ROUTINES Called (in and below this routine):
!
!       sad_strat_calc
!         nat_sat_temp ...... derives the NAT saturation temp
!
!         sulfate_sad_calc .. Calculates the sulfate SAD;  HNO3, H2O cond. phase
!         calc_radius_lbs ... Calculates the radius for a H2SO4 Binary soln. (T>200K)
!         sad_to_h2so4 ...... Derives H2SO4 abundance (micrograms m-3)
!                             from SAGEII SAD.
!         density............ A. Tabazedeh binary thermo model
!         equil.............. A. Tabazedeh ternary thermo. model
!
!         nat_sad_calc....... Calculates the NAT SAD; HNO3, H2O cond. phase
!         nat_cond........... Derives the NAT HNO3 gas/cond partitioning
!
!         ice_sad_calc....... derives the ICE SAD and H2O gas/cond partitioning
!===============================================================================

! Subprogram not used       subroutine sad_strat_calc( lchnk, m, press, temper, hno3_gas, &
! Subprogram not used                                  hno3_cond, h2o_gas, h2o_cond, hcl_gas, hcl_cond, &
! Subprogram not used                                  sad_sage, radius_strat, sad_strat, ncol, pbuf )
! Subprogram not used 
! Subprogram not used       use cam_history, only : outfld
! Subprogram not used       use physics_buffer, only : physics_buffer_desc
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used       integer, intent(in)     ::  lchnk                      ! chnk id
! Subprogram not used       integer, intent(in)     ::  ncol                       ! columns in chunk
! Subprogram not used       real(r8), intent(in)    ::  m           (ncol,pver)    ! Air density (molecules cm-3)
! Subprogram not used       real(r8), intent(in)    ::  sad_sage    (ncol,pver)    ! SAGEII surface area density (cm2 aer. cm-3 air)
! Subprogram not used       real(r8), intent(in)    ::  press       (ncol,pver)    ! Pressure, hPa
! Subprogram not used       real(r8), intent(in)    ::  temper      (pcols,pver)   ! Temperature (K)
! Subprogram not used       real(r8), intent(inout) ::  h2o_gas     (ncol,pver)    ! H2O gas-phase (mole fraction)
! Subprogram not used       real(r8), intent(inout) ::  h2o_cond    (ncol,pver)    ! H2O condensed phase  (mole fraction)
! Subprogram not used       real(r8), intent(inout) ::  hno3_gas    (ncol,pver)    ! HNO3 condensed phase (mole fraction)
! Subprogram not used       real(r8), intent(inout) ::  hno3_cond   (ncol,pver,2)  ! HNO3 condensed phase (mole fraction)
! Subprogram not used       real(r8), intent(inout) ::  hcl_gas     (ncol,pver)    ! HCL gas-phase        (mole fraction)
! Subprogram not used       real(r8), intent(inout) ::  hcl_cond    (ncol,pver)    ! HCL condensed phase  (mole fraction)
! Subprogram not used       real(r8), intent(out)   ::  radius_strat(ncol,pver,3)  ! Radius of Sulfate, NAT, and ICE (cm)
! Subprogram not used       real(r8), intent(out)   ::  sad_strat   (ncol,pver,3)  ! Surface area density of Sulfate, NAT, ICE (cm2 cm-3)
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: temp_floor = 0._r8
! Subprogram not used 
! Subprogram not used       integer  ::  i, k, n
! Subprogram not used       integer  ::  dims(1)
! Subprogram not used       real(r8) ::  hno3_total    (ncol,pver)     ! HNO3 total = gas-phase + condensed
! Subprogram not used       real(r8) ::  h2o_total     (ncol,pver)     ! H2O total  = gas-phase + condensed
! Subprogram not used       real(r8) ::  radius_lbs    (ncol,pver)     ! Radius of Liquid Binary Sulfate (cm)
! Subprogram not used       real(r8) ::  radius_sulfate(ncol,pver)     ! Radius of Sulfate aerosol (cm)
! Subprogram not used       real(r8) ::  radius_nat    (ncol,pver)     ! Radius of NAT aerosol     (cm)
! Subprogram not used       real(r8) ::  radius_ice    (ncol,pver)     ! Radius of ICE aerosol     (cm)
! Subprogram not used       real(r8) ::  sad_nat       (ncol,pver)     ! SAD of NAT aerosol        (cm2 cm-3)
! Subprogram not used       real(r8) ::  sad_sulfate   (ncol,pver)     ! SAD of Sulfate aerosol    (cm2 cm-3)
! Subprogram not used       real(r8) ::  sad_ice       (ncol,pver)     ! SAD of ICE aerosol        (cm2 cm-3)
! Subprogram not used       real(r8) ::  tsat_nat      (ncol,pver)     ! Temperature for NAT saturation
! Subprogram not used       real(r8) ::  h2o_avail     (ncol,pver)     ! H2O temporary arrays
! Subprogram not used       real(r8) ::  hno3_avail    (ncol,pver)     ! HNO3 temporary array
! Subprogram not used       real(r8) ::  hno3_gas_nat  (ncol,pver)     ! HNO3 after call to NAT routines
! Subprogram not used       real(r8) ::  hno3_gas_sulf (ncol,pver)     ! HNO3 after call to STS routines
! Subprogram not used       real(r8) ::  hno3_cond_nat (ncol,pver)     ! HNO3 condensed after call to NAT
! Subprogram not used       real(r8) ::  hno3_cond_sulf(ncol,pver)     ! HNO3 condensed after call to STS routines
! Subprogram not used       real(r8) ::  hcl_total     (ncol,pver)     ! HCl total  = gas-phase + condensed
! Subprogram not used       real(r8) ::  hcl_avail     (ncol,pver)     ! HCL temporary arrays
! Subprogram not used       real(r8) ::  hcl_gas_sulf  (ncol,pver)     ! HCL after call to STS routines
! Subprogram not used       real(r8) ::  hcl_cond_sulf (ncol,pver)     ! HCL condensed after call to STS routines
! Subprogram not used       real(r8) ::  temp          (pcols,pver)    ! wrk temperature array
! Subprogram not used       real(r8) ::  h2so4m        (ncol,pver)     ! wrk array
! Subprogram not used 
! Subprogram not used       logical  ::  z_val(ncol)
! Subprogram not used 
! Subprogram not used       logical  ::  mask_lbs(ncol,pver)           ! LBS mask T: P>300hPa; P<2hPa2hPa; SAGE<1e-15; T>200K
! Subprogram not used       logical  ::  mask_ice(ncol,pver)           ! ICE mask T: .not. mask_lbs; h2o_cond>0
! Subprogram not used       logical  ::  mask_sts(ncol,pver)           ! STS mask T: .not. mask_sts
! Subprogram not used       logical  ::  mask_nat(ncol,pver)		 ! NAT mask T: mask_sts=T; T<Tsat_nat
! Subprogram not used       type(physics_buffer_desc), pointer :: pbuf(:)
! Subprogram not used 
! Subprogram not used       real(r8), parameter :: eighty_percent = 0.8_r8
! Subprogram not used       real(r8), parameter :: twenty_percent = 0.2_r8
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... initialize to zero
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       do n = 1,3
! Subprogram not used          do k = 1,pver
! Subprogram not used             radius_strat(:,k,n) = 0._r8
! Subprogram not used             sad_strat(:,k,n)    = 0._r8
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used       do n = 1,2
! Subprogram not used          do k = 1,pver
! Subprogram not used             hno3_cond(:,k,n) = 0._r8
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used !
! Subprogram not used       do k = 1,pver
! Subprogram not used             h2o_total     (:,k) = 0._r8
! Subprogram not used !
! Subprogram not used             h2o_avail     (:,k) = 0._r8
! Subprogram not used             hno3_avail    (:,k) = 0._r8
! Subprogram not used             hcl_avail     (:,k) = 0._r8
! Subprogram not used !
! Subprogram not used             hno3_total    (:,k) = 0._r8
! Subprogram not used             hno3_gas_nat  (:,k) = 0._r8
! Subprogram not used             hno3_gas_sulf (:,k) = 0._r8
! Subprogram not used             hno3_cond_nat (:,k) = 0._r8
! Subprogram not used             hno3_cond_sulf(:,k) = 0._r8
! Subprogram not used !
! Subprogram not used             hcl_total     (:,k) = 0._r8
! Subprogram not used             hcl_cond_sulf (:,k) = 0._r8
! Subprogram not used             hcl_gas_sulf  (:,k) = 0._r8
! Subprogram not used       end do       
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... limit temperature
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       do k = 1,pver
! Subprogram not used          h2so4m(:ncol,k)  = 0._r8
! Subprogram not used          temp(:ncol,k)    = max( temp_floor,temper(:ncol,k) )
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... total HNO3 and H2O gas and condensed phases
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used          hno3_total(:,k) = hno3_gas(:,k) + hno3_cond(:,k,1) + hno3_cond(:,k,2) 
! Subprogram not used 	 h2o_total(:,k)  = h2o_gas(:,k)  + h2o_cond(:,k)
! Subprogram not used          hcl_total (:,k) = hcl_gas(:,k)  + hcl_cond(:,k)
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !======================================================================
! Subprogram not used !======================================================================
! Subprogram not used !     ... Set SAD to SAGEII if Temperature is GT 200K and Return
! Subprogram not used !
! Subprogram not used !     ... mask_lbs  = true  .... T > 200K or SAD_SULF <= 1e-15 or 
! Subprogram not used !                                P < 2hPa or P > 300hPa
! Subprogram not used !     ... mask_sts  = false .... .not. mask_lbs
! Subprogram not used !     ... mask_nat  = false .... T <= TSAT_NAT
! Subprogram not used !     ... mask_ice  = false .... H2O_COND > 0.0
! Subprogram not used !======================================================================
! Subprogram not used !======================================================================
! Subprogram not used 
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used          mask_lbs(:,k) = temp(:ncol,k) > 200._r8 .or. sad_sage(:,k) <= 1.e-15_r8 &
! Subprogram not used                          .or. press(:ncol,k) < 2._r8 .or. press(:ncol,k) > 300._r8
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used sage_sad : &
! Subprogram not used       if( any( mask_lbs(:,sad_topp:pver) ) ) then
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_lbs(:,k) )
! Subprogram not used 	       sad_strat(:,k,1)    = sad_sage(:,k)
! Subprogram not used 	       sad_strat(:,k,2)    = 0._r8
! Subprogram not used                sad_strat(:,k,3)    = 0._r8
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... Calculate Liquid Binary Sulfate Radius
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used          call calc_radius_lbs( ncol, mask_lbs, sad_sage, radius_lbs )
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_lbs(:,k) )
! Subprogram not used                radius_strat(:,k,1)     = radius_lbs(:,k)
! Subprogram not used 	       radius_strat(:,k,2)     = 0._r8
! Subprogram not used 	       radius_strat(:,k,3)     = 0._r8
! Subprogram not used 	       hno3_gas    (:,k)       = hno3_total(:,k)
! Subprogram not used                hno3_cond   (:,k,1)     = 0._r8
! Subprogram not used 	       hno3_cond   (:,k,2)     = 0._r8
! Subprogram not used 	       hcl_gas     (:,k)       = hcl_total(:,k)
! Subprogram not used                hcl_cond    (:,k)       = 0._r8
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used          if( all( mask_lbs(:,sad_topp:pver) ) ) then
! Subprogram not used             call outfld( 'H2SO4M_C', h2so4m(:ncol,:), ncol, lchnk )
! Subprogram not used 	    return
! Subprogram not used          end if
! Subprogram not used       end if sage_sad
! Subprogram not used 
! Subprogram not used !======================================================================
! Subprogram not used !======================================================================
! Subprogram not used !     ... Logic for deriving ICE
! Subprogram not used !         Ice formation occurs here if condensed phase H2O exists.
! Subprogram not used !
! Subprogram not used !         mask_lbs  = false.... T > 200K or SAD_SULF < 1e-15 or 
! Subprogram not used !                               P >2hPa or P <300hPa
! Subprogram not used !         mask_ice  = true .... H2O_COND > 0.0
! Subprogram not used !======================================================================
! Subprogram not used !======================================================================
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 do i = 1,ncol
! Subprogram not used 	    if( .not. mask_lbs(i,k) ) then
! Subprogram not used                mask_ice(i,k) = h2o_cond(i,k) > 0._r8
! Subprogram not used 	    else
! Subprogram not used 	       mask_ice(i,k) = .false.
! Subprogram not used 	    end if
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used all_ice : &
! Subprogram not used       if( any( mask_ice(:,sad_topp:pver) ) ) then
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_ice(:,k) )
! Subprogram not used 	       h2o_avail(:,k) = h2o_cond(:,k)
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !        ... ICE 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used          call ice_sad_calc( ncol, press, temp, m, h2o_avail, &
! Subprogram not used 			    sad_ice, radius_ice, mask_ice )
! Subprogram not used 
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_ice(:,k) )
! Subprogram not used                sad_strat   (:,k,3) = sad_ice       (:,k)
! Subprogram not used                radius_strat(:,k,3) = radius_ice    (:,k)
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used       end if all_ice
! Subprogram not used 
! Subprogram not used !======================================================================
! Subprogram not used !======================================================================
! Subprogram not used !  	... LOGIC for STS and NAT 
! Subprogram not used !
! Subprogram not used !           mask_lbs = false .... T > 200K or SAD_SULF <= 1e-15 or 
! Subprogram not used !                                 P < 2hPa or P >300hPa
! Subprogram not used !           mask_sts = true  .... not  mask_lbs
! Subprogram not used !           mask_nat = true  .... T <= TSAT_NAT and  mask_sts = true
! Subprogram not used !======================================================================
! Subprogram not used !======================================================================
! Subprogram not used 
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 do i = 1,ncol
! Subprogram not used 	    if( .not. mask_lbs(i,k) ) then
! Subprogram not used                mask_sts(i,k) = .true.
! Subprogram not used 	    else
! Subprogram not used 	       mask_sts(i,k) = .false.
! Subprogram not used 	    end if
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  	... STS (80% of total HNO3 logic)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !        NOTE: STS only sees 80% of the total HNO3 (Wegner et al., JGR, 2013a)
! Subprogram not used sts_nat_sad : &
! Subprogram not used       if( any( mask_sts(:,sad_topp:pver) ) ) then
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_sts(:,k) )
! Subprogram not used                h2o_avail (:,k)   = h2o_gas   (:,k)
! Subprogram not used                hno3_avail(:,k)   = hno3_total(:,k)*eighty_percent
! Subprogram not used                hcl_avail (:,k)   = hcl_total  (:,k)
! Subprogram not used             endwhere
! Subprogram not used 	    if( any(mask_sts(:,k)) ) then
! Subprogram not used 	       where( mask_sts(:,k) )
! Subprogram not used 	          z_val(:) = hno3_avail(:,k) == 0._r8
! Subprogram not used 	       elsewhere
! Subprogram not used 	          z_val(:) = .false.
! Subprogram not used 	       endwhere
! Subprogram not used 	       if( any( z_val(:) ) ) then
! Subprogram not used 	          write(iulog,*) 'sad_strat_calc: Before CHEM Sulfate_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
! Subprogram not used 	       end if
! Subprogram not used 	    end if
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, hcl_avail, &
! Subprogram not used                                 sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, &
! Subprogram not used                                 hcl_gas_sulf, hcl_cond_sulf, sad_sulfate, &
! Subprogram not used                                 radius_sulfate, mask_sts, lchnk, 1, h2so4m, .true.)
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_sts(:,k) )
! Subprogram not used                sad_strat   (:,k,1) = sad_sulfate   (:,k)
! Subprogram not used                radius_strat(:,k,1) = radius_sulfate(:,k)
! Subprogram not used                hno3_gas    (:,k)   = hno3_gas_sulf (:,k)
! Subprogram not used                hno3_cond   (:,k,1) = hno3_cond_sulf(:,k)
! Subprogram not used                hcl_gas     (:,k)   = hcl_gas_sulf  (:,k)
! Subprogram not used                hcl_cond    (:,k)   = hcl_cond_sulf (:,k)
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... NAT (20% of total HNO3 logic)
! Subprogram not used !     ... using total H2O and gas-phase HNO3 after STS calc
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !        NOTE: NAT only sees 20% of the total HNO3 (Wegner et al., JGR, 2013a)
! Subprogram not used          hno3_avail(:,:) = hno3_total(:,:)*twenty_percent
! Subprogram not used          call nat_sat_temp( ncol, hno3_avail, h2o_avail, press, tsat_nat, mask_sts)
! Subprogram not used 
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used 	   do i = 1,ncol
! Subprogram not used 	     if( .not. mask_lbs(i,k) ) then
! Subprogram not used                mask_nat(i,k) = temp(i,k) <= tsat_nat(i,k)
! Subprogram not used 	     else
! Subprogram not used 	       mask_nat(i,k) = .false.
! Subprogram not used 	     end if
! Subprogram not used            end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_nat(:,k) )
! Subprogram not used                h2o_avail (:,k) = h2o_gas      (:,k)
! Subprogram not used                hno3_avail(:,k) = hno3_total   (:,k)*twenty_percent
! Subprogram not used             endwhere
! Subprogram not used 	    if( any(mask_nat(:,k)) ) then
! Subprogram not used 	       where( mask_nat(:,k) )
! Subprogram not used 	          z_val(:) = hno3_avail(:,k) == 0._r8
! Subprogram not used 	       elsewhere
! Subprogram not used 	          z_val(:) = .false.
! Subprogram not used 	       endwhere
! Subprogram not used 	       if( any( z_val(:) ) ) then
! Subprogram not used 	          write(iulog,*) 'sad_nat_calc: After CHEM Sulf_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
! Subprogram not used 	       end if
! Subprogram not used 	    end if
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used          call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
! Subprogram not used 			    hno3_gas_nat, hno3_cond_nat, &
! Subprogram not used                             sad_nat, radius_nat, mask_nat )
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! NOTE:  Add in gas-phase from STS with gas-phase from NAT
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             where( mask_nat(:,k) )
! Subprogram not used                sad_strat   (:,k,2) = sad_nat       (:,k)
! Subprogram not used                radius_strat(:,k,2) = radius_nat    (:,k)
! Subprogram not used                hno3_gas    (:,k)   = hno3_gas_sulf (:,k) + hno3_gas_nat  (:,k)
! Subprogram not used                hno3_cond   (:,k,2) = hno3_cond_nat (:,k)
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used ! NOTE:  If NAT does not form (in STS region), need to add the 20% of the total HNO3 back to gas-phase
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used 	   do i = 1,ncol
! Subprogram not used              if ( .not. mask_nat(i,k) ) then
! Subprogram not used                if ( mask_sts(i,k) ) then
! Subprogram not used                  hno3_gas   (i,k) = hno3_gas_sulf(i,k) + hno3_total(i,k)*twenty_percent
! Subprogram not used                end if
! Subprogram not used              end if
! Subprogram not used            end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used       end if sts_nat_sad
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call outfld( 'H2SO4M_C', h2so4m(:ncol,:), ncol, lchnk )
! Subprogram not used 
! Subprogram not used       end subroutine sad_strat_calc

! Subprogram not used       subroutine nat_sat_temp( ncol, hno3_total, h2o_avail, press, tsat_nat, mask )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer,  intent(in)   :: ncol
! Subprogram not used       real(r8), intent(in)   :: press(ncol,pver)
! Subprogram not used       real(r8), intent(in)   :: h2o_avail(ncol,pver)
! Subprogram not used       real(r8), intent(in)   :: hno3_total(ncol,pver)
! Subprogram not used       real(r8), intent(out)  :: tsat_nat(ncol,pver)
! Subprogram not used       logical,  intent(in)   :: mask(ncol,pver)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: ssrNAT = 10.0_r8
! Subprogram not used       real(r8), parameter :: ssrNATi = .1_r8
! Subprogram not used       real(r8), parameter :: aa = -2.7836_r8, &
! Subprogram not used                              bb = -0.00088_r8, &
! Subprogram not used                              cc = 38.9855_r8, &
! Subprogram not used                              dd = -11397.0_r8, &
! Subprogram not used                              ee = 0.009179_r8
! Subprogram not used       integer  :: k, i
! Subprogram not used       real(r8) :: bbb                     ! temporary variable
! Subprogram not used       real(r8) :: ccc                     ! temporary variable
! Subprogram not used       real(r8) :: wrk                     ! temporary variable
! Subprogram not used       real(r8) :: tmp                     ! temporary variable
! Subprogram not used       real(r8) :: phno3                   ! hno3 partial pressure
! Subprogram not used       real(r8) :: ph2o                    ! h2o  partial pressure
! Subprogram not used 
! Subprogram not used       tsat_nat(:,:) = 0._r8
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... Derive HNO3 and H2O partial pressure (torr)
! Subprogram not used !          where: 0.7501 = 760/1013.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used          do i = 1,ncol
! Subprogram not used             if( mask(i,k) ) then
! Subprogram not used 
! Subprogram not used                bbb   = press(i,k) * .7501_r8
! Subprogram not used                phno3 = hno3_total(i,k) * bbb
! Subprogram not used                ph2o  = h2o_avail(i,k)  * bbb
! Subprogram not used 
! Subprogram not used                if( phno3 > 0._r8 ) then
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... Calculate NAT Saturation Threshold Temperature
! Subprogram not used !           Hanson and Mauersberger: GRL, Vol.15, 8, p855-858, 1988.
! Subprogram not used !           Substitute m(T) and b(T) into Equation (1). Rearrange and
! Subprogram not used !           solve quadratic eqation. 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                   tmp = log10( ph2o )
! Subprogram not used                   wrk = 1._r8 / (ee + bb*tmp)
! Subprogram not used                   bbb = (aa*tmp - log10( phno3*ssrNATi ) + cc) * wrk
! Subprogram not used                   ccc = dd *wrk
! Subprogram not used                   tsat_nat(i,k) = .5_r8 * (-bbb + sqrt( bbb*bbb - 4._r8*ccc ))
! Subprogram not used                endif
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine nat_sat_temp

! Subprogram not used       subroutine calc_radius_lbs( ncol, mask, sad_sage, radius_lbs )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol
! Subprogram not used       real(r8), intent(in)  :: sad_sage(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: radius_lbs(ncol,pver)
! Subprogram not used       logical, intent(in)   :: mask(ncol,pver)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer  :: k
! Subprogram not used       real(r8) :: lbs_vol_dens(ncol)       ! Vol Density (cm3 aer / cm3 air)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... parameters
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: lbs_part_dens = 10._r8, sigma_lbs = 1.6_r8
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... calculate the volume density (cm3 aerosol / cm3 air)
! Subprogram not used !     	    calculate the mean radius for binary soln
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 where( mask(:,k) )
! Subprogram not used             lbs_vol_dens(:) = ((sad_sage(:,k)**1.5_r8)/3._r8)/sqrt( four_pi*lbs_part_dens ) &
! Subprogram not used                               *exp( 1.5_r8*(log( sigma_lbs ))**2 )
! Subprogram not used             radius_lbs(:,k) = (3._r8*lbs_vol_dens(:)/(four_pi*lbs_part_dens))**one_thrd &
! Subprogram not used                               *exp( -1.5_r8*(log( sigma_lbs ))**2 )
! Subprogram not used 	 endwhere
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine calc_radius_lbs

! Subprogram not used       subroutine ice_sad_calc( ncol, press, temp, m, h2o_avail, &
! Subprogram not used 			       sad_ice, radius_ice, mask )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol
! Subprogram not used       real(r8), intent(in)  :: press     (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: temp      (pcols,pver)
! Subprogram not used       real(r8), intent(in)  :: m         (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: h2o_avail (ncol,pver)
! Subprogram not used       real(r8), intent(out) :: sad_ice   (ncol,pver)
! Subprogram not used       real(r8), intent(out) :: radius_ice(ncol,pver)
! Subprogram not used       logical, intent(in)   :: mask      (ncol,pver)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: &
! Subprogram not used                  avo_num       = 6.02214e23_r8, &
! Subprogram not used                  aconst        = -2663.5_r8, &
! Subprogram not used                  bconst        = 12.537_r8, &
! Subprogram not used                  ice_mass_dens = 1._r8, &
! Subprogram not used                  ice_part_dens = 1.e-1_r8, &
! Subprogram not used                  mwh2o         = 18._r8, &
! Subprogram not used                  sigma_ice     = 1.6_r8, &
! Subprogram not used                  ice_dens_aer  = ice_mass_dens / (mwh2o/avo_num), &
! Subprogram not used                  ice_dens_aeri = 1._r8/ice_dens_aer
! Subprogram not used 
! Subprogram not used       integer  :: k
! Subprogram not used       real(r8) :: h2o_cond_ice(ncol)      ! Condensed phase H2O (from 1)
! Subprogram not used       real(r8) :: voldens_ice (ncol)      ! Volume Density, um3 cm-3
! Subprogram not used 
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 where( mask(:,k) )
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Convert condensed phase to molecules cm-3 units
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 	   h2o_cond_ice(:) = h2o_avail(:,k) * m(:,k)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... ICE volume density .....
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used            voldens_ice(:) = h2o_cond_ice(:)*ice_dens_aeri
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the SAD from log normal distribution .....
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used            sad_ice(:,k) = (four_pi*ice_part_dens)**one_thrd &
! Subprogram not used                          *(3._r8*voldens_ice(:))**two_thrd &
! Subprogram not used                          *exp( -(log( sigma_ice ))**2 )
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !    .... Calculate the radius from log normal distribution .....
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used            radius_ice(:,k) = (3._r8*h2o_cond_ice(:) &
! Subprogram not used                               /(ice_dens_aer*four_pi*ice_part_dens))**one_thrd &
! Subprogram not used                              *exp( -1.5_r8*(log( sigma_ice ))**2 )
! Subprogram not used          endwhere
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine ice_sad_calc

! Subprogram not used       subroutine sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, hcl_avail, &
! Subprogram not used                                    sad_sage, m, hno3_gas, hno3_cond, &
! Subprogram not used                                    hcl_gas, hcl_cond, sad_sulfate, &
! Subprogram not used                                    radius_sulfate, mask, lchnk, flag, h2so4m, is_chem)
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol
! Subprogram not used       integer, intent(in)   :: lchnk, flag
! Subprogram not used       real(r8), intent(in)  :: temp       (pcols,pver)
! Subprogram not used       real(r8), intent(in)  :: press      (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: m          (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: h2o_avail  (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: hno3_avail (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: hcl_avail  (ncol,pver)   
! Subprogram not used       real(r8), intent(in)  :: sad_sage   (ncol,pver)
! Subprogram not used       real(r8), intent(out) :: hno3_gas   (ncol,pver)       ! Gas-phase HNO3, mole fraction
! Subprogram not used       real(r8), intent(out) :: hno3_cond  (ncol,pver)       ! Condensed phase HNO3, mole fraction
! Subprogram not used       real(r8), intent(out) :: hcl_gas    (ncol,pver)       ! Gas-phase HCL, mole fraction
! Subprogram not used       real(r8), intent(out) :: hcl_cond   (ncol,pver)       ! Condensed phase HCL, mole fraction
! Subprogram not used       real(r8), intent(out) :: sad_sulfate(ncol,pver)   
! Subprogram not used       real(r8), intent(out) :: radius_sulfate(ncol,pver)
! Subprogram not used       real(r8), intent(inout) :: h2so4m   (ncol,pver)       ! mass per volume, micro grams m-3
! Subprogram not used       logical, intent(in)   :: is_chem                      ! chemistry calc switch
! Subprogram not used       logical, intent(in)   :: mask       (ncol,pver)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: t_limit           = 200._r8
! Subprogram not used       real(r8), parameter :: avo_num           = 6.02214e23_r8
! Subprogram not used       real(r8), parameter :: mwh2so4           = 98.076_r8
! Subprogram not used       real(r8), parameter :: sigma_sulfate     = 1.6_r8
! Subprogram not used       real(r8), parameter :: sulfate_part_dens = 10._r8
! Subprogram not used 
! Subprogram not used       integer  :: i, k
! Subprogram not used       integer  :: cnt1, cnt2
! Subprogram not used       real(r8) :: ratio
! Subprogram not used       real(r8) :: vals(2)
! Subprogram not used       real(r8) :: h2so4_aer_dens  (ncol,pver)  ! grams cm-3 solution
! Subprogram not used       real(r8) :: h2so4_cond      (ncol,pver)  ! Condensed H2SO4 (moles cm-3 air)
! Subprogram not used       real(r8) :: sulfate_vol_dens(ncol,pver)  ! Volume Density, cm3 aerosol  cm-3 air
! Subprogram not used       real(r8) :: wtf             (ncol,pver)  ! weight fraction of H2SO4 in ternary soln
! Subprogram not used       real(r8) :: wts             (ncol,pver)  ! weight percent of ternary solution
! Subprogram not used 
! Subprogram not used ! Carslaw HCl solubility
! Subprogram not used       real(r8) :: wts0            (ncol,pver)  ! weight percent of H2SO4 is LBS
! Subprogram not used       real(r8) :: wtn             (ncol,pver)  ! weight percent of HNO3 in STS
! Subprogram not used       real(r8) :: ch2so4          (ncol,pver)  ! Total H2SO4 (moles / cm3 of air)
! Subprogram not used       real(r8) :: molh2so4        (ncol,pver)  ! Equil molality of H2SO4 in STS
! Subprogram not used       real(r8) :: molhno3         (ncol,pver)  ! Equil molality of HNO3 in STS
! Subprogram not used       real(r8) :: AD              (ncol,pver)  ! air density (molecules cm-3)
! Subprogram not used       real(r8) :: xmf             (ncol,pver)  ! 
! Subprogram not used       real(r8) :: hhcl            (ncol,pver)  ! henry's solubility of hcl in binary
! Subprogram not used       real(r8) :: phcl0           (ncol,pver)  ! partial pressure of hcl (hPa)
! Subprogram not used       real(r8) :: h2so4vmr        (ncol,pver)  ! atmospheric mole fraction of H2SO4
! Subprogram not used       real(r8) :: nsul            (ncol,pver)  ! moles / m3- H2SO4 pure liquid 
! Subprogram not used       real(r8) :: mcl             (ncol,pver)  ! molality of hcl in ?
! Subprogram not used       real(r8) :: wtcl            (ncol,pver)  ! 
! Subprogram not used       real(r8) :: phcl            (ncol,pver)  ! partial pressure of hcl (over aerosol)
! Subprogram not used       real(r8) :: parthcl         (ncol,pver)  ! fraction of HCl in gas-phase
! Subprogram not used !
! Subprogram not used       real(r8) :: packer          (ncol*pver)
! Subprogram not used       logical  :: do_equil                     ! local mask
! Subprogram not used       logical  :: mask_lt         (ncol,pver)  ! local temperature mask
! Subprogram not used       logical  :: maskx           (ncol,pver)
! Subprogram not used       logical  :: converged       (ncol,pver)  ! EQUIL convergence test
! Subprogram not used 
! Subprogram not used       hcl_gas (:,:) = 0.0_r8
! Subprogram not used       hcl_cond(:,:) = 0.0_r8
! Subprogram not used       parthcl (:,:) = 0.0_r8
! Subprogram not used       phcl0   (:,:) = 0.0_r8
! Subprogram not used 
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used          mask_lt(:,k)  = mask(:,k)
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  	... derive H2SO4 (micro grams / m3) from SAGEII SAD
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       call sad2h2so4( h2o_avail, press, sad_sage, temp, sulfate_vol_dens, & 
! Subprogram not used                       h2so4_aer_dens, h2so4m, mask, ncol )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  	... limit h2so4m
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used             do i = 1,ncol
! Subprogram not used                if( mask(i,k) ) then
! Subprogram not used                   if( h2so4m(i,k) <= 0._r8 ) then
! Subprogram not used                      h2so4m(i,k) = 1.e-2_r8
! Subprogram not used                   end if
! Subprogram not used                end if
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used       if( is_chem ) then
! Subprogram not used       else
! Subprogram not used         do k = sad_topp,pver
! Subprogram not used             where( mask(:ncol,k) )
! Subprogram not used                mask_lt(:ncol,k) = temp(:ncol,k) < t_limit
! Subprogram not used             end where
! Subprogram not used          end do
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       if( is_chem ) then
! Subprogram not used          do_equil = .true.
! Subprogram not used       else
! Subprogram not used          do_equil = any( mask_lt(:,sad_topp:pver) )
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the ternary soln volume density
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       if( do_equil ) then
! Subprogram not used 
! Subprogram not used          call equil( temp, h2so4m, hno3_avail, h2o_avail, press, & 
! Subprogram not used                      hno3_cond, h2so4_cond, wts, wtn, wts0, molh2so4, molhno3, mask_lt, ncol, &
! Subprogram not used                      lchnk, flag, is_chem, converged )       
! Subprogram not used   
! Subprogram not used          do k = sad_topp,pver
! Subprogram not used 
! Subprogram not used 	   where( ( mask_lt(:ncol,k) ) .AND. ( converged(:ncol,k) ) )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... convert h2o, hno3 from moles cm-3 air to molecules cm-3 air
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                hno3_cond(:ncol,k) = min( hno3_cond(:ncol,k),hno3_avail(:ncol,k) )
! Subprogram not used                hno3_gas(:ncol,k)  = hno3_avail(:ncol,k) - hno3_cond(:ncol,k)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Derive ternary volume density (cm3 aer / cm3 air)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                wtf(:ncol,k) = .01_r8* wts(:ncol,k)
! Subprogram not used                sulfate_vol_dens(:ncol,k) = h2so4_cond(:ncol,k)*mwh2so4/(wtf(:ncol,k)*h2so4_aer_dens(:ncol,k))
! Subprogram not used 
! Subprogram not used ! Carslaw solubility
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Partition HCl (gas/condensed) *** Carslaw
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !          THE SOLUBILITY OF HCL 
! Subprogram not used !          HHCl (MOL/KG/ATM) taken form Shi et al., JGR 2001
! Subprogram not used !
! Subprogram not used 
! Subprogram not used !     .... Convert weight % to weight fraction 
! Subprogram not used                 wtn(:ncol,k) = wtn(:ncol,k) * 0.01_r8 
! Subprogram not used                 wts0(:ncol,k) = wts0(:ncol,k) * 0.01_r8
! Subprogram not used 
! Subprogram not used !     .... Derive xmf (mole fraction H2SO4 in LBS )
! Subprogram not used                 xmf(:ncol,k) = (wts0(:ncol,k)*100.0_r8)/((wts0(:ncol,k)*100.0_r8)+ &
! Subprogram not used                            (100.0_r8-(wts0(:ncol,k)*100._r8))*98.0_r8/18.0_r8)
! Subprogram not used 
! Subprogram not used !     .... Derive hhcl (henry's solubility of hcl in binary)
! Subprogram not used                 hhcl(:ncol,k) = (0.094_r8-0.61_r8*xmf(:ncol,k)+1.2_r8*xmf(:ncol,k)**2.0_r8) &
! Subprogram not used                             *exp(-8.68_r8+(8515.0_r8-10718.0_r8*xmf(:ncol,k)**(0.7_r8))/temp(:ncol,k)) 
! Subprogram not used 
! Subprogram not used !     .... Derive phcl0 (partial pressure of hcl( hPa))
! Subprogram not used                 phcl0(:ncol,k) = hcl_avail(:ncol,k)*press(:ncol,k) / 1013.26_r8
! Subprogram not used 
! Subprogram not used !     .... Derive H2SO4 vmr (h2so4_cond = mole / cm-3)
! Subprogram not used                 AD(:ncol,k) = (6.022098e23_r8 * press(:ncol,k) / 1013.26_r8) &
! Subprogram not used                            / (temp(:ncol,k)*8.2058e-2_r8*1000.0_r8)
! Subprogram not used                 h2so4vmr(:ncol,k)  = (h2so4_cond(:ncol,k)*6.022098e23_r8) / AD(:ncol,k)
! Subprogram not used 
! Subprogram not used !     .... Derive nsul (moles / m3 H2SO4 pure liquid )
! Subprogram not used                 nsul(:ncol,k)  = h2so4vmr(:ncol,k) * press(:ncol,k) * 100.0_r8 / 8.314_r8 / temp(:ncol,k)
! Subprogram not used 
! Subprogram not used !     .... Derive  mcl (molality of hcl)
! Subprogram not used                 mcl(:ncol,k) = (1.0_r8/8.314e-5_r8/temp(:ncol,k)*phcl0(:ncol,k))/(nsul(:ncol,k)/molh2so4(:ncol,k) + &
! Subprogram not used                            1.0_r8/(8.314e-5_r8)/temp(:ncol,k)/hhcl(:ncol,k)) 
! Subprogram not used 
! Subprogram not used !     .... Derive wtcl ( )
! Subprogram not used                 wtcl(:ncol,k) = mcl(:ncol,k)*36.5_r8/(1000.0_r8 + 98.12_r8*molh2so4(:ncol,k) + 63.03_r8*molhno3(:ncol,k)) 
! Subprogram not used 
! Subprogram not used !     .... Derive phcl (partial pressure over the aerosol)
! Subprogram not used                 phcl(:ncol,k) = mcl(:ncol,k)/hhcl(:ncol,k) 
! Subprogram not used 
! Subprogram not used !     .... Derive parhcl (fraction of HCl in gas-phase)
! Subprogram not used                 where(phcl0(:ncol,k)>0._r8)
! Subprogram not used                   parthcl(:ncol,k) = 1.0_r8 - (phcl0(:ncol,k) - phcl(:ncol,k)) / phcl0(:ncol,k)
! Subprogram not used                 elsewhere
! Subprogram not used                   parthcl(:ncol,k) = 0._r8
! Subprogram not used                 endwhere
! Subprogram not used 
! Subprogram not used !     .... Partition HCl (gas/condensed)  
! Subprogram not used                 hcl_gas (:ncol,k)  = hcl_avail(:ncol,k) * parthcl(:ncol,k)
! Subprogram not used                 hcl_cond(:ncol,k)  = hcl_avail(:ncol,k) - hcl_gas(:ncol,k)
! Subprogram not used 
! Subprogram not used             elsewhere
! Subprogram not used                hno3_cond(:ncol,k) = 0.0_r8
! Subprogram not used                hno3_gas(:ncol,k)  = hno3_avail(:ncol,k)
! Subprogram not used                hcl_cond (:ncol,k) = 0.0_r8
! Subprogram not used                hcl_gas  (:ncol,k) = hcl_avail(:ncol,k)
! Subprogram not used 
! Subprogram not used             endwhere
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 where( mask(:,k) )
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the SAD (assuming ternary solution)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             sad_sulfate(:,k) = (four_pi*sulfate_part_dens)**one_thrd &
! Subprogram not used                                *(3._r8*sulfate_vol_dens(:,k))**two_thrd &
! Subprogram not used                                *exp( -(log( sigma_sulfate ))**2 )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the radius (assuming ternary solution) (in cm?)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             radius_sulfate(:,k) = (3._r8*sulfate_vol_dens(:,k) &
! Subprogram not used                                    /(four_pi*sulfate_part_dens))**one_thrd &
! Subprogram not used                                   *exp( -1.5_r8*(log( sigma_sulfate ))**2 )
! Subprogram not used 	 endwhere
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine sulfate_sad_calc


! Subprogram not used       subroutine nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
! Subprogram not used 			       hno3_gas, hno3_cond, sad_nat, radius_nat, mask )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol
! Subprogram not used       real(r8), intent(in)  :: press     (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: m         (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: temp      (pcols,pver)
! Subprogram not used       real(r8), intent(in)  :: h2o_avail (ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: hno3_avail(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: hno3_cond (ncol,pver)       ! HNO3 in condensed phase (mole fraction)
! Subprogram not used       real(r8), intent(out) :: hno3_gas  (ncol,pver)       ! HNO3 in gas-phase (mole fraction)
! Subprogram not used       real(r8), intent(out) :: sad_nat   (ncol,pver)   
! Subprogram not used       real(r8), intent(out) :: radius_nat(ncol,pver)
! Subprogram not used       logical, intent(in)   :: mask(ncol,pver)             ! grid mask
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer  :: k, i
! Subprogram not used       real(r8) :: nat_dens_condphase(ncol, pver)      ! Condensed phase NAT, molec cm-3
! Subprogram not used       real(r8) :: voldens_nat       (ncol, pver)      ! Volume Density,  um3 cm-3
! Subprogram not used       real(r8) :: hno3_cond_total   (ncol, pver)      ! Total Condensed phase HNO3 
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... parameters
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: avo_num          = 6.02214e23_r8, &
! Subprogram not used                              nat_mass_dens    = 1.6_r8, &
! Subprogram not used                              nat_part_dens    = 1.0e-2_r8, &
! Subprogram not used                              mwnat            = 117._r8, &
! Subprogram not used                              sigma_nat        = 1.6_r8, &
! Subprogram not used                              nat_dens_aer     = nat_mass_dens / (mwnat/avo_num), &
! Subprogram not used                              nat_dens_aeri    = 1._r8/nat_dens_aer
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     ... Derive HNO3 paritioning (call A. Tabazedeh routine for NAT)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       call nat_cond( ncol, press, temp, h2o_avail, hno3_avail, &
! Subprogram not used 		     hno3_gas, hno3_cond_total, mask )
! Subprogram not used 
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used          do i = 1,ncol
! Subprogram not used masked :   if( mask(i,k) ) then
! Subprogram not used             
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Set Condensed phase for return arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             hno3_cond(i,k) = hno3_cond_total(i,k)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculated Condensed Phase NAT (i.e. HNO3) in
! Subprogram not used !            molecules cm-3 of air units.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             nat_dens_condphase(i,k) = hno3_cond_total(i,k) * m(i,k)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the Volume Density
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             voldens_nat(i,k) = nat_dens_condphase(i,k) * nat_dens_aeri
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the SAD from log normal distribution
! Subprogram not used !     .... Assuming sigma and nad_part_dens (# particles per cm3 of air)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             sad_nat(i,k) = (four_pi*nat_part_dens)**(one_thrd) &
! Subprogram not used                           *(3._r8*voldens_nat(i,k))**(two_thrd) &
! Subprogram not used                           *exp( -(log( sigma_nat )**2 ))
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Calculate the radius of NAT from log normal distribution
! Subprogram not used !     .... Assuming sigma and nat_part_dens (# particles per cm3 
! Subprogram not used !     .... of air)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used             radius_nat(i,k) = (3._r8*nat_dens_condphase(i,k) &
! Subprogram not used                         /(nat_dens_aer*four_pi*nat_part_dens))**(one_thrd) &
! Subprogram not used                         *exp( -1.5_r8*(log( sigma_nat ))**2 )
! Subprogram not used 
! Subprogram not used            end if masked
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine nat_sad_calc

! Subprogram not used       subroutine nat_cond( ncol, press, temp, h2o_avail, hno3_avail, &
! Subprogram not used 			   hno3_gas, hno3_cond, mask )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used ! 	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol
! Subprogram not used       real(r8), intent(in)  :: press(ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: temp(pcols,pver)
! Subprogram not used       real(r8), intent(in)  :: h2o_avail(ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: hno3_avail(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: hno3_gas(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: hno3_cond(ncol,pver)
! Subprogram not used       logical, intent(in)   :: mask(ncol,pver)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used ! 	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter ::  aa = -2.7836_r8,  &
! Subprogram not used                               bb = -0.00088_r8, &
! Subprogram not used                               cc = 38.9855_r8,  &
! Subprogram not used                               dd = -11397.0_r8, &
! Subprogram not used                               ee = 0.009179_r8
! Subprogram not used 
! Subprogram not used       integer  :: i, k
! Subprogram not used       real(r8) :: bt                                 ! temporary variable
! Subprogram not used       real(r8) :: mt                                 ! temporary variable
! Subprogram not used       real(r8) :: t                                  ! temporary variable
! Subprogram not used       real(r8) :: logPhno3                           ! temporary variable
! Subprogram not used       real(r8) :: phno3                              ! hno3 partial pressure
! Subprogram not used       real(r8) :: ph2o                               ! h2o  partial pressure
! Subprogram not used       real(r8) :: phno3_eq                           ! partial pressure above NAT
! Subprogram not used       real(r8) :: wrk      
! Subprogram not used       
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 do i = 1,ncol
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     .... Derive HNO3 and H2O partial pressure (torr)
! Subprogram not used !          where: 0.7501 = 760/1013.
! Subprogram not used !----------------------------------------------------------------------        
! Subprogram not used 	    if( mask(i,k) ) then
! Subprogram not used                wrk   = press(i,k) * .7501_r8
! Subprogram not used                phno3 = hno3_avail(i,k) * wrk
! Subprogram not used                ph2o  = h2o_avail(i,k)  * wrk
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     Calculating the temperature coefficients for the variation of HNO3
! Subprogram not used !     and H2O vapor pressure (torr) over a trihydrate solution of HNO3/H2O
! Subprogram not used !     The coefficients are taken from Hanson and Mauersberger: 
! Subprogram not used !     GRL, Vol.15, 8, p855-858, 1988.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 	       t   = temp(i,k)
! Subprogram not used                bt  = cc + dd/t + ee*t
! Subprogram not used                mt  = aa + bb*t
! Subprogram not used   
! Subprogram not used                logphno3 = mt*log10( ph2o ) + bt
! Subprogram not used                phno3_eq = 10._r8**logphno3
! Subprogram not used 
! Subprogram not used 	       if( phno3 > phno3_eq ) then
! Subprogram not used                   wrk            = 1._r8 / wrk
! Subprogram not used                   hno3_cond(i,k) = (phno3 - phno3_eq) * wrk
! Subprogram not used                   hno3_gas(i,k)  = phno3_eq * wrk
! Subprogram not used 	       else
! Subprogram not used                   hno3_cond(i,k) = 0._r8
! Subprogram not used                   hno3_gas(i,k)  = hno3_avail(i,k)
! Subprogram not used 	       end if
! Subprogram not used 	    end if
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       end subroutine nat_cond

! Subprogram not used       subroutine sad2h2so4( h2o, press, sad_sage, temp, lbs_vol_dens, &
! Subprogram not used                             h2so4_aer_dens, h2so4m, mask, ncol )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: ncol                                    ! columns in chunk
! Subprogram not used       real(r8), intent(in)  :: h2o(ncol,pver)                          ! h2o mole fraction
! Subprogram not used       real(r8), intent(in)  :: sad_sage(ncol,pver)                     ! sad from SAGEII cm2 aer, cm-3 air
! Subprogram not used       real(r8), intent(in)  :: press(ncol,pver)                        ! pressure (hPa)
! Subprogram not used       real(r8), intent(in)  :: temp(pcols,pver)                        ! temperature (K)
! Subprogram not used       real(r8), intent(inout) :: h2so4m(ncol,pver)                       ! microgram/m**3 of air,
! Subprogram not used       real(r8), intent(out) :: h2so4_aer_dens(ncol,pver)               ! units: grams / cm3-aerosol
! Subprogram not used       real(r8), intent(out) :: lbs_vol_dens(ncol,pver)                 ! cm3 aer / cm3 air
! Subprogram not used       logical, intent(in)   :: mask(ncol,pver)                         ! activation mask
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: lbs_part_dens = 10._r8
! Subprogram not used       real(r8), parameter :: sigma_lbs     = 1.6_r8
! Subprogram not used       real(r8), parameter :: t_floor       = 180._r8
! Subprogram not used 
! Subprogram not used       integer  :: i, k, l
! Subprogram not used       real(r8) :: wts0   
! Subprogram not used       real(r8) :: p                         ! pressure, torr
! Subprogram not used       real(r8) :: tr                        ! inverse temperature
! Subprogram not used       real(r8) :: c(6)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... Calculate the volume density (cm3 aerosol / cm3 air)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       do k = sad_topp,pver
! Subprogram not used 	 do i = 1,ncol
! Subprogram not used 	    if( mask(i,k) ) then
! Subprogram not used                lbs_vol_dens(i,k) = ((sad_sage(i,k)**1.5_r8)/3._r8)/sqrt( four_pi*lbs_part_dens ) &
! Subprogram not used                                    *exp( 1.5_r8*(log( sigma_lbs ))**2 )
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... calculate Molality from Tabazadeh EQUIL routine (binary soln)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... DEK, added a minimum to temperature
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                p   = h2o(i,k) * press(i,k) * .7501_r8
! Subprogram not used                tr  = 1._r8 / max( t_floor,temp(i,k) )
! Subprogram not used              
! Subprogram not used 	       do l = 1,6
! Subprogram not used                   c(l) = exp( a(1,l) + tr*(a(2,l) + tr*(a(3,l) + tr*(a(4,l) + tr*a(5,l)))) )
! Subprogram not used 	       end do
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... H2SO4/H2O pure weight percent and molality (moles gram-1)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                wts0 = max( 0._r8,c(1) + p*(-c(2) + p*(c(3) + p*(-c(4) + p*(c(5) - p*c(6))))) )
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... derive weight fraction for density routine
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                wts0 = .01_r8 *wts0
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... calculate the binary solution density, grams / cm3-aerosol
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                h2so4_aer_dens(i,k) = max( 0._r8,density( temp(i,k), wts0 ) )
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !     	... calculate the H2SO4 micrograms m-3 abundance for binary soln
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                h2so4m(i,k) = lbs_vol_dens(i,k)*h2so4_aer_dens(i,k)*wts0*1.e12_r8
! Subprogram not used 	    end if
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used    
! Subprogram not used       end subroutine sad2h2so4

!======================================================================
!
! ROUTINE
!   EQUIL
!
!   Date...
!     7 October 1999
!
!   Programmed by...
!     A. Tabazadeh
!
! DESCRIPTION
!	Ternary solution routine
!
! ARGUMENTS
!
!....  INPUT:
!
!        H2SO4m    = microgram/m**3 of air
!        HNO3r     = mole fraction
!        H2Or      = mole fraction
!        PTOTAL    = hPa
!
!....  Output
!
!        Cwater    = Total moles of liguid water / cm**3 of air
!        hno3_cond = HNO3 Condensed phase (mole fraction)
!        CH2SO4    = Total H2SO4 moles / cm**3 of air
!        WTS       = Weight percent of H2SO4 in the ternary aerosol
!
!======================================================================
! Subprogram not used       subroutine equil( temper, h2so4m, hno3_avail, h2o_avail, press, &
! Subprogram not used                         hno3_cond, ch2so4, wts, wtn, wts0, molh2so4, molhno3, mask, ncol, &
! Subprogram not used                         lchnk, flag, is_chem, converged)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !                       Written by Azadeh Tabazadeh (1993)
! Subprogram not used !           (modified from EQUISOLV -- M. Z. Jacobson -- see below) 
! Subprogram not used !            NASA Ames Research Center , Tel. (415) 604 - 1096
! Subprogram not used !
! Subprogram not used !       This program solves the equilibrium composition for the ternary  
! Subprogram not used !       system of H2SO4/HNO3/H2O under typical stratospheric conditions. 
! Subprogram not used !       The formulation of this work is described by Tabazadeh, A.,      
! Subprogram not used !       Turco, R. P., and Jacobson, M. Z. (1994), "A model for studying  
! Subprogram not used !       the composition and chemical effects of stratospheric aerosols," 
! Subprogram not used !       J. Geophys. Res., 99, 12,897, 1994.        *
! Subprogram not used !
! Subprogram not used !       The solution mechanism for the equilibrium equations is des-     
! Subprogram not used !       cribed by Jacobson, M. Z., Turco, R. P., and Tabazadeh, A.       
! Subprogram not used !       (1994), "Simulating Equilibrium within aerosols and non-equil-   
! Subprogram not used !       ibrium between gases and aerosols," J. Geophys. Res., in review. 
! Subprogram not used !       The mechanism is also codified in the fortran program, EQUISOLV, 
! Subprogram not used !       by M.Z. Jacobson (1991-3). EQUISOLV solves any number of         
! Subprogram not used !       gas / liquid / solid / ionic equilibrium equations simultan-     
! Subprogram not used !       eously and includes treatment of the water equations and act-    
! Subprogram not used !       ivity coefficients. The activity coeffients currently in         
! Subprogram not used !       EQUISOLV are valid for tropospheric temperatures. The acitiv-    
! Subprogram not used !       ities listed here are valid for stratospheric temperatures only. 
! Subprogram not used !
! Subprogram not used !	DEFINING PARAMETERS
! Subprogram not used !
! Subprogram not used !       *NOTE*	  Solver parameters including, F, Z, QN, QD, and deltaX
! Subprogram not used !                 are described in Jacobson et al.
! Subprogram not used !
! Subprogram not used !	PTOTAL	= Total atmospheric pressure in mb
! Subprogram not used !	H2SO4m	= Total mass of H2SO4 (microgram/m**3 of air)
! Subprogram not used !	HNO3r	= HNO3 mixing ratio
! Subprogram not used !	H2Or	= H2O mixing ratio
! Subprogram not used !	P	    = Partial pressure of water in units of torr
! Subprogram not used !	pures	= molality for a pure H2SO4/H2O system
! Subprogram not used !	puren	= molality for a pure HNO3/H2O sytem
! Subprogram not used !	WTS0	= weight percent of H2SO4 in a pure H2SO4/H2O system
! Subprogram not used !	WTN0	= weight percent of HNO3 in a pure HNO3/H2O system
! Subprogram not used !	WTS	    = weight percent of H2SO4 in the ternary aerosol
! Subprogram not used !	WTN 	= weight percent of HNO3 in the ternary aerosol
! Subprogram not used !	PHNO3	= HNO3 vapor pressure over the ternary system in atm
! Subprogram not used !	HNO3	= HNO3 vapor concentration over the ternary system (#/cm3)
! Subprogram not used !	CH2SO4	= Total H2SO4 moles / cm**3 of air
! Subprogram not used !	CHNO3	= Total HNO3 moles / cm**3 of air
! Subprogram not used !	CHplus	= Total H+ moles / cm**3 0f air
! Subprogram not used !	CPHNO3	= Total moles of HNO3 gas / cm**3 of air
! Subprogram not used !	CNO3	= Total moles of NO3- / cm**3 0f air
! Subprogram not used !	Cwater	= Total moles of liguid water / cm**3 of air
! Subprogram not used !	KS  	= Solubility constant for dissolution of HNO3 in
! Subprogram not used !		      water ( HNO3(gas) === H+(aq) + NO3- (aq) )
! Subprogram not used !	nm  	= HNO3 molality at the STREN of the ternary solution
! Subprogram not used !	sm  	= H2SO4 molality at the STREN of the ternary solution
! Subprogram not used !	molHNO3	= Equilibrium molality of HNO3 in the ternary solution
! Subprogram not used !	molH2SO4= Equilibrium molality of H2SO4 in the ternary solution
! Subprogram not used !     STREN   = ionic strenght for the ternary solutin, which in
! Subprogram not used !		      this case is = 3 * molH2SO4 + molHNO3
! Subprogram not used !	acts	= Pure mean binary activity coefficient for the H2SO4/
! Subprogram not used !		      H2O system evaluated at the STREN of the ternary system
! Subprogram not used !	actn	= Pure mean binary activity coefficient for the HNO3/
! Subprogram not used !		      H2O system evaluated at the STREN of the ternary system
! Subprogram not used !     ymix    = Mixed binary activity coefficient for the HNO3/H2O in
! Subprogram not used !		      the ternary solution
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use abortutils, only : endrun
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       integer, intent(in)   :: lchnk
! Subprogram not used       integer, intent(in)   :: flag
! Subprogram not used       integer, intent(in)   :: ncol                         ! columns in chunk
! Subprogram not used       real(r8), intent(in)  :: h2so4m(ncol,pver)    
! Subprogram not used       real(r8), intent(in)  :: hno3_avail(ncol,pver)    
! Subprogram not used       real(r8), intent(in)  :: h2o_avail(ncol,pver)    
! Subprogram not used       real(r8), intent(in)  :: press(ncol,pver)
! Subprogram not used       real(r8), intent(in)  :: temper(pcols,pver)
! Subprogram not used       real(r8), intent(out) :: hno3_cond(ncol,pver)    
! Subprogram not used       real(r8), intent(out) :: ch2so4(ncol,pver)    
! Subprogram not used       real(r8), intent(out) :: wts(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: wtn(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: wts0(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: molh2so4(ncol,pver)
! Subprogram not used       real(r8), intent(out) :: molhno3(ncol,pver)
! Subprogram not used       logical, intent(in)   :: is_chem
! Subprogram not used       logical, intent(in)   :: mask(ncol,pver)              ! activation mask
! Subprogram not used       logical, intent(out)  :: converged(ncol,pver)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !      integer, parameter  :: itermax = 50
! Subprogram not used       integer, parameter  :: itermax = 100
! Subprogram not used       real(r8), parameter :: con_lim  = .00005_r8
! Subprogram not used       real(r8), parameter :: t0       = 298.15_r8
! Subprogram not used       real(r8), parameter :: ks0      = 2.45e6_r8
! Subprogram not used       real(r8), parameter :: lower_delx = 1.e-10_r8
! Subprogram not used       real(r8), parameter :: upper_delx = .98_r8
! Subprogram not used       real(r8), parameter :: con_crit_chem = 5.e-5_r8
! Subprogram not used 
! Subprogram not used       integer  :: i, iter, k, l, nstep
! Subprogram not used       real(r8) :: reduction_factor
! Subprogram not used       real(r8) :: p
! Subprogram not used       real(r8) :: tr
! Subprogram not used       real(r8) :: wtn0
! Subprogram not used       real(r8) :: pures
! Subprogram not used       real(r8) :: puren
! Subprogram not used       real(r8) :: chno3
! Subprogram not used       real(r8) :: chplus
! Subprogram not used       real(r8) :: cno3
! Subprogram not used       real(r8) :: wrk
! Subprogram not used       real(r8) :: z, num, den
! Subprogram not used       real(r8) :: deltax
! Subprogram not used       real(r8) :: chplusnew
! Subprogram not used       real(r8) :: cno3new
! Subprogram not used       real(r8) :: stren
! Subprogram not used       real(r8) :: sm
! Subprogram not used       real(r8) :: actn
! Subprogram not used       real(r8) :: acts
! Subprogram not used       real(r8) :: nm
! Subprogram not used       real(r8) :: ks
! Subprogram not used       real(r8) :: lnks
! Subprogram not used       real(r8) :: lnks0
! Subprogram not used       real(r8) :: mixyln
! Subprogram not used       real(r8) :: wrk_h2so4
! Subprogram not used       real(r8) :: cphno3new
! Subprogram not used       real(r8) :: con_val
! Subprogram not used       real(r8) :: t, t1, t2, f, f1, f2, ymix, hplus, wtotal, ratio 
! Subprogram not used       real(r8) :: con_crit
! Subprogram not used       real(r8) :: h2o_cond(ncol,pver)
! Subprogram not used       real(r8) :: fratio(0:itermax)
! Subprogram not used       real(r8) :: delx(0:itermax)
! Subprogram not used       real(r8) :: delz(0:itermax)
! Subprogram not used       real(r8) :: c(12)
! Subprogram not used       real(r8) :: d(13:22)
! Subprogram not used       logical  :: interval_set
! Subprogram not used       logical  :: positive
! Subprogram not used 
! Subprogram not used       converged(:,:) = .false.
! Subprogram not used 
! Subprogram not used       lnks0 = log( ks0 )
! Subprogram not used       if( is_chem ) then
! Subprogram not used          con_crit = con_crit_chem
! Subprogram not used       else
! Subprogram not used          con_crit = con_crit_chem
! Subprogram not used       end if
! Subprogram not used       Level_loop : do k = sad_topp,pver
! Subprogram not used          Column_loop : do i = 1,ncol
! Subprogram not used             if( mask(i,k) ) then
! Subprogram not used                p = h2o_avail(i,k) * press(i,k) * .7501_r8
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	Calculating the molality for pure binary systems of H2SO4/H2O
! Subprogram not used !	and HNO3/H2O at a given temperature and water vapor pressure
! Subprogram not used !	profile (relative humiditiy). Water activities were used to
! Subprogram not used !	calculate the molalities as described in Tabazadeh et al. (1994).
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                t  = max( 180._r8,temper(i,k) )
! Subprogram not used                tr = 1._r8/t
! Subprogram not used                do l = 1,12
! Subprogram not used                   c(l) = exp( a(1,l) + tr*(a(2,l) + tr*(a(3,l) + tr*(a(4,l) + tr*a(5,l)))) )
! Subprogram not used                end do
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... H2SO4/H2O pure weight percent and molality
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                wts0(i,k)  = max( 0.01_r8,c(1) + p*(-c(2) + p*(c(3) + p*(-c(4) + p*(c(5) - p*c(6))))) )
! Subprogram not used                pures = (wts0(i,k) * 1000._r8)/(100._r8 - wts0(i,k))
! Subprogram not used                pures = pures / 98._r8
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... HNO3/H2O pure weight percent and molality
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                puren = max( 0._r8,c(7) + p*(-c(8) + p*(c(9) + p*(-c(10) + p*(c(11) - p*c(12))))) )
! Subprogram not used !              wtn0 = (puren * 6300._r8) /(puren * 63._r8 + 1000._r8)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	The solving scheme is described both in Jacobson et al. and Tabazadeh
! Subprogram not used !	et al.. Assumptions:
! Subprogram not used !	(1) H2SO4 is present only in the aqueous-phase
! Subprogram not used !	(2) H2SO4 and HNO3 in solution are fully dissocated into H+
! Subprogram not used !	    SO42- and NO3-
! Subprogram not used !	(3) PHNO3 + NO3- = constant
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 	       ch2so4(i,k) = (h2so4m(i,k)*1.e-12_r8) / 98._r8
! Subprogram not used 	       if( pures > 0._r8 ) then
! Subprogram not used 	          wrk_h2so4 = (1000._r8*ch2so4(i,k))/(pures*18._r8)
! Subprogram not used 	       else
! Subprogram not used 	          wrk_h2so4 = 0._r8
! Subprogram not used 	       end if
! Subprogram not used 	       chno3 = 1.2029e-5_r8 * press(i,k) * tr * hno3_avail(i,k)
! Subprogram not used 	       do l = 13,22
! Subprogram not used                   d(l) = b(1,l) + t*(b(2,l) + t*(b(3,l) + t*(b(4,l) + t*b(5,l))))
! Subprogram not used 	       end do
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	Note that KS depends only on the temperature
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                t1	= (t - t0)/(t*t0)
! Subprogram not used                t2	= t0/t - 1._r8 - log( t0/t )
! Subprogram not used                lnks     = lnks0 - 8792.3984_r8 * t1  - 16.8439_r8 * t2
! Subprogram not used                ks	= exp( lnks )
! Subprogram not used 
! Subprogram not used 	       converged(i,k) = .false.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	Setting up initial guesses for the equations above.  Note that
! Subprogram not used !	for the initial choices the mass and the charge must be conserved.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                delx(0)  = .5_r8
! Subprogram not used                z        = .5_r8
! Subprogram not used                delz(0)  = .5_r8
! Subprogram not used                fratio(0) = 0._r8
! Subprogram not used                reduction_factor = .1_r8
! Subprogram not used                interval_set = .false.
! Subprogram not used Iter_loop :    do iter = 1,itermax
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	Cwater is the water equation as described in Tabazadeh et
! Subprogram not used !	al. and Jacobson et al.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                   cno3new   = chno3 * delx(iter-1)
! Subprogram not used                   cphno3new = chno3 * (1._r8 - delx(iter-1))
! Subprogram not used 		  if( puren > 0._r8 ) then
! Subprogram not used 		     t1 = (1000._r8*cno3new)/(puren*18._r8)
! Subprogram not used 		  else
! Subprogram not used 		     t1 = 0._r8
! Subprogram not used 		  end if
! Subprogram not used 	          h2o_cond(i,k) = t1 + wrk_h2so4
! Subprogram not used 		  if( h2o_cond(i,k) > 0._r8 ) then
! Subprogram not used                      wrk      = 1.e3_r8 / (18._r8 * h2o_cond(i,k))
! Subprogram not used                      molhno3(i,k)  = cno3new * wrk
! Subprogram not used                      molh2so4(i,k) = ch2so4(i,k) * wrk
! Subprogram not used 		  else
! Subprogram not used                      molhno3(i,k)  = 0._r8
! Subprogram not used                      molh2so4(i,k) = 0._r8
! Subprogram not used 		  end if
! Subprogram not used                   stren	= molhno3(i,k) + 3._r8 * molh2so4(i,k)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	(1) Calculate the activity of H2SO4 at a given STREN
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                   sm	= stren/3._r8
! Subprogram not used                   acts = d(13) + sm*(d(14) + sm*(d(15) + sm*d(16)))
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	(2) Calculate the activity for HNO3 at a given STREN
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                   nm	= stren
! Subprogram not used                   actn 	= d(17) + nm*(d(18) + nm*(d(19) + nm*(d(20) + nm*(d(21) + nm*d(22)))))
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	(3) Calculate the mixed activity coefficient for HNO3 at STREN
! Subprogram not used !	    as described by Tabazadeh et al.
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                   f1	 = 2._r8 * (molh2so4(i,k) + molhno3(i,k)) * actn
! Subprogram not used                   f2	 = 2.25_r8 * molh2so4(i,k) * acts
! Subprogram not used 
! Subprogram not used                   if (stren > 0._r8) then
! Subprogram not used                     mixyln = (f1 + f2) / (2._r8 * stren)
! Subprogram not used                   else
! Subprogram not used                     mixyln = 0._r8
! Subprogram not used                   end if
! Subprogram not used                   ymix	 = exp( mixyln )
! Subprogram not used                   hplus	 = 2._r8 * molh2so4(i,k) + molhno3(i,k)
! Subprogram not used                   num = ymix**2 * hplus * molhno3(i,k)
! Subprogram not used                   den = 1000._r8 * cphno3new * .0820578_r8 * t * ks
! Subprogram not used 		  if( chno3 == 0._r8 ) then
! Subprogram not used 	             converged(i,k) = .true.
! Subprogram not used 		     exit Iter_loop
! Subprogram not used 		  end if
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !       the denominator
! Subprogram not used !       Calculate the ratio F, check convergence
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !       Calculate the ratio F and reset the deltaX (see Jacobson et al.)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !!DEK
! Subprogram not used !       When the numerator is zero, it can drive the denominator
! Subprogram not used !       to 0, which resulted in a NaN for f and also the fraction
! Subprogram not used !       ratio. Assume that in this case, the limit of f would
! Subprogram not used !       really approach 1, not infinity and thus converge the
! Subprogram not used !       solution.
! Subprogram not used                   if ((num .eq. 0._r8) .and. (den .eq. 0._r8)) then
! Subprogram not used                     f = 1._r8
! Subprogram not used                   else
! Subprogram not used                     f = num / den
! Subprogram not used                   end if
! Subprogram not used                   fratio(iter) = abs( f ) - 1._r8
! Subprogram not used 		  con_val      = abs( f - 1._r8 )
! Subprogram not used 		  if( con_val <= con_lim ) then
! Subprogram not used 		     converged(i,k)  = .true.
! Subprogram not used 		     exit Iter_loop
! Subprogram not used                   end if
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !       non-convergence; setup next iterate
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used                   if( interval_set ) then
! Subprogram not used                      z = reduction_factor * z
! Subprogram not used                      delz(iter) = z
! Subprogram not used 		     if( f > 1._r8 ) then
! Subprogram not used                         deltax = -z
! Subprogram not used 		     else
! Subprogram not used                         deltax = z
! Subprogram not used 		     end if
! Subprogram not used                      delx(iter) = delx(iter-1) + deltax
! Subprogram not used                   else
! Subprogram not used                      if( iter == 1 ) then
! Subprogram not used                         if( fratio(iter) >= 1._r8 ) then
! Subprogram not used                            positive = .false.
! Subprogram not used                         else
! Subprogram not used                            positive = .true.
! Subprogram not used                         end if
! Subprogram not used                      end if
! Subprogram not used                      if( fratio(iter)*fratio(iter-1) < 0._r8 ) then
! Subprogram not used                         interval_set = .true.
! Subprogram not used                         reduction_factor = .5_r8
! Subprogram not used                         delx(iter) = .5_r8*(delx(iter-1) + delx(iter-2))
! Subprogram not used                         z = .5_r8*abs( delx(iter-1) - delx(iter-2) )
! Subprogram not used                      else
! Subprogram not used                         if( .not. positive ) then
! Subprogram not used                            delx(iter) = reduction_factor * delx(iter-1)
! Subprogram not used                         else
! Subprogram not used                            delx(iter) = reduction_factor + delx(iter-1)
! Subprogram not used                            if( delx(iter) > upper_delx ) then
! Subprogram not used                               delx(iter) = .5_r8
! Subprogram not used                               interval_set = .true.
! Subprogram not used                               reduction_factor = .5_r8
! Subprogram not used                            end if
! Subprogram not used                         end if
! Subprogram not used                      end if
! Subprogram not used 		  end if
! Subprogram not used                end do Iter_loop
! Subprogram not used 
! Subprogram not used                wtotal   = molhno3(i,k) * 63._r8 + molh2so4(i,k) * 98._r8 + 1000._r8
! Subprogram not used                wts(i,k) = (molh2so4(i,k) * 9800._r8) / wtotal
! Subprogram not used                wtn(i,k) = (molhno3(i,k) *6300._r8)/ wtotal
! Subprogram not used 	       if( cno3new /= 0._r8 .or. cphno3new /= 0._r8 ) then
! Subprogram not used                   ratio	= max( 0._r8,min( 1._r8,cno3new/(cphno3new + cno3new) ) )
! Subprogram not used                   hno3_cond(i,k) = ratio*hno3_avail(i,k)
! Subprogram not used 	       else
! Subprogram not used                   hno3_cond(i,k) = 0._r8
! Subprogram not used 	       end if
! Subprogram not used                if( .not. converged(i,k) ) then
! Subprogram not used                   write(iulog,*) 'equil: Failed to converge @ is_chem,flag,lchnk,i,k,f = ',is_chem,flag,lchnk,i,k,f
! Subprogram not used                   write(iulog,*) '       wts0,pures,puren,chno3,ch2so4 = ',wts0(i,k),pures,puren,chno3,ch2so4(i,k)
! Subprogram not used                   write(iulog,*) '       stren,mixyln,ymix,hplus,num,den = ',stren,mixyln,ymix,hplus,num,den
! Subprogram not used                   write(iulog,*) '       h2o_avail,hno3_avail,p,t = ',h2o_avail(i,k),hno3_avail(i,k),press(i,k),temper(i,k)
! Subprogram not used                   write(iulog,*) '       molhno3,molh2so4,h2o_cond,hno3_cond = ', &
! Subprogram not used                                  molhno3(i,k),molh2so4(i,k),h2o_cond(i,k),hno3_cond(i,k)
! Subprogram not used                   if( con_val > .05_r8 ) then
! Subprogram not used                      write(iulog,*) ' '
! Subprogram not used                      write(iulog,*) 'equil; diagnostics at lchnk, flag, i, k, iter = ',lchnk,flag,i,k,iter
! Subprogram not used                      write(iulog,*) 'equil; fratio'
! Subprogram not used                      write(iulog,'(5(1pg15.7))') fratio(0:iter-1)
! Subprogram not used                      write(iulog,*) ' '
! Subprogram not used                      write(iulog,*) 'equil; delx'
! Subprogram not used                      write(iulog,'(5(1pg15.7))') delx(0:iter-1)
! Subprogram not used                      write(iulog,*) ' '
! Subprogram not used                      write(iulog,*) 'equil; delz'
! Subprogram not used                      write(iulog,'(5(1pg15.7))') delz(0:iter-1)
! Subprogram not used                      write(iulog,*) ' '
! Subprogram not used                   else if( iter > 50 ) then
! Subprogram not used                      write(iulog,*) 'equil: Iterations are beyond 50, number of iter = ',iter
! Subprogram not used                      write(iulog,*) 'equil: converged @ is_chem,flag,lchnk,i,k = '
! Subprogram not used                      write(iulog,*) is_chem,flag,lchnk,i,k
! Subprogram not used                      write(iulog,*) 'equil: converged @ f, num, den = '
! Subprogram not used                      write(iulog,*) f, num, den
! Subprogram not used                      write(iulog,*) '       h2o_avail,hno3_avail,p,t = '
! Subprogram not used                      write(iulog,*) h2o_avail(i,k),hno3_avail(i,k),press(i,k),temper(i,k)
! Subprogram not used                      write(iulog,*) '       molhno3(i,k),molh2so4(i,k),h2o_cond,hno3_cond = '
! Subprogram not used                      write(iulog,*) molhno3(i,k),molh2so4(i,k),h2o_cond(i,k),hno3_cond(i,k)
! Subprogram not used                   end if
! Subprogram not used 	       end if
! Subprogram not used             end if
! Subprogram not used          end do Column_loop
! Subprogram not used       end do Level_loop
! Subprogram not used 
! Subprogram not used       end subroutine equil

!======================================================================
!
!
! ROUTINE
!   DENSITY
!
!   Date...
!     7 October 1999
!
!   Programmed by...
!     A. Tabazadeh
!
! DESCRIPTION
!     Calculates the density (g cm-3) of a binary sulfate solution.
!
! ARGUMENTS
!   INPUT
!      T           Temperature
!      w           Weight fraction
!
!   OUTPUT
!        den       Density of the Binary Solution (g cm-3)
!
!======================================================================
       
! Subprogram not used       function density( temp, w )
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... Dummy arguments
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), intent(in) :: temp, w
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... Function declarations
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8) :: density
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... Local variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       real(r8), parameter :: a9 = -268.2616e4_r8, a10 = 576.4288e3_r8
! Subprogram not used 
! Subprogram not used       real(r8) :: a0, a1, a2, a3, a4, a5, a6, a7 ,a8
! Subprogram not used       real(r8) :: c1, c2, c3, c4
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... Temperature variables
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       c1 = temp - 273.15_r8
! Subprogram not used       c2 = c1**2
! Subprogram not used       c3 = c1*c2
! Subprogram not used       c4 = c1*c3
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	Polynomial Coefficients
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       a0 = 999.8426_r8 + 334.5402e-4_r8*c1 - 569.1304e-5_r8*c2
! Subprogram not used       a1 = 547.2659_r8 - 530.0445e-2_r8*c1 + 118.7671e-4_r8*c2 + 599.0008e-6_r8*c3
! Subprogram not used       a2 = 526.295e1_r8 + 372.0445e-1_r8*c1 + 120.1909e-3_r8*c2 - 414.8594e-5_r8*c3 + 119.7973e-7_r8*c4
! Subprogram not used       a3 = -621.3958e2_r8 - 287.7670_r8*c1 - 406.4638e-3_r8*c2 + 111.9488e-4_r8*c3 + 360.7768e-7_r8*c4
! Subprogram not used       a4 = 409.0293e3_r8 + 127.0854e1_r8*c1 + 326.9710e-3_r8*c2 - 137.7435e-4_r8*c3 - 263.3585e-7_r8*c4
! Subprogram not used       a5 = -159.6989e4_r8 - 306.2836e1_r8*c1 + 136.6499e-3_r8*c2 + 637.3031e-5_r8*c3
! Subprogram not used       a6 = 385.7411e4_r8 + 408.3717e1_r8*c1 - 192.7785e-3_r8*c2
! Subprogram not used       a7 = -580.8064e4_r8 - 284.4401e1_r8*c1
! Subprogram not used       a8 = 530.1976e4_r8 + 809.1053_r8*c1
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !	... Summation
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used       density = .001_r8*(a0 + w*(a1 + w*(a2 + w*(a3 + w*(a4 + w*(a5 + w*(a6 + w*(a7 + w*(a8 + w*(a9 + w*a10))))))))))
! Subprogram not used 
! Subprogram not used       end function density

      end module mo_sad
