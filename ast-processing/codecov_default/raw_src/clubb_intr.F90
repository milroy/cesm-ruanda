module clubb_intr

  !----------------------------------------------------------------------------------------------------- !
  ! Module to interface 1 with Cloud Layers Unified by Bi-normals (CLUBB), developed                   !
  !    by the University of Wisconsin Milwaukee Group (UWM).                                             !
  !                                                                                                      !
  ! CLUBB replaces the exisiting turbulence, shallow convection, and macrophysics in CAM5                !  
  !                                                                                                      !  
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      ! 
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! Authors:  P. Bogenschutz, C. Craig, A. Gettelman                                                     ! 
  !                                                                                                      ! 
  !----------------------------------------------------------------------------------------------------- !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pver, pverp
  use phys_control,  only: phys_getopts
  use physconst,     only: rair, cpair, gravit, latvap, latice, zvir, rh2o, karman, tms_orocnst, tms_z0fac
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc 
  use constituents,  only: pcnst
  use pbl_utils,     only: calc_ustar, calc_obklen
  use mpishorthand

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: clubb_ini_cam, clubb_register_cam, clubb_tend_cam, &
            stats_end_timestep_clubb, & 
	    clubb_surface



  ! ------------ !
  ! Private data !
  ! ------------ !

  integer, parameter :: &
      grid_type    = 3, &  		! The 2 option specifies stretched thermodynamic levels
      hydromet_dim = 0     		! The hydromet array in SAM-CLUBB is currently 0 elements
   
  real(r8), dimension(0) :: &
      sclr_tol = 1.e-8_r8               ! Total water in kg/kg

  character(len=6), parameter :: &
      saturation_equation = "flatau" 	! Flatau polynomial approximation for SVP

  real(r8), parameter :: &
      theta0   = 300._r8, &		! Reference temperature                     [K]
      ts_nudge = 86400._r8, & 		! Time scale for u/v nudging (not used)     [s]
      p0_clubb = 100000._r8
      
  real(r8), parameter :: &
      host_dx = 100000._r8, &		! Host model deltax [m]			     
      host_dy = 100000._r8		! Host model deltay [m]
      
  integer, parameter :: & 
    sclr_dim = 0   			! Higher-order scalars, set to zero

!  Constant parameters
  logical, parameter, private :: &
    l_uv_nudge       = .false.,       &  ! Use u/v nudging (not used)
    l_implemented    = .true.,        &  ! Implemented in a host model (always true)
    l_host_applies_sfc_fluxes = .false.  ! Whether the host model applies the surface fluxes
    
  logical            :: do_tms
  logical            :: lq(pcnst)
  logical            :: prog_modal_aero

  integer            :: edsclr_dim       ! Number of scalars to transport in CLUBB
 
!  define physics buffer indicies here       
  integer :: &
    wp2_idx, &         			! vertical velocity variances
    wp3_idx, &         			! third moment of vertical velocity
    wpthlp_idx, &      			! turbulent flux of thetal
    wprtp_idx, &       			! turbulent flux of total water
    rtpthlp_idx, &     			! covariance of thetal and rt
    rtp2_idx, &        			! variance of total water
    thlp2_idx, &       			! variance of thetal
    up2_idx, &         			! variance of east-west wind
    vp2_idx, &         			! variance of north-south wind
    upwp_idx, &        			! east-west momentum flux
    vpwp_idx, &        			! north-south momentum flux
    thlm_idx, &        			! mean thetal
    rtm_idx, &         			! mean total water mixing ratio
    um_idx, &         			! mean of east-west wind
    vm_idx, &           		! mean of north-south wind
    cld_idx, &         			! Cloud fraction
    concld_idx, &       		! Convective cloud fraction
    ast_idx, &          		! Stratiform cloud fraction
    alst_idx, &         		! Liquid stratiform cloud fraction
    aist_idx, & 			! Ice stratiform cloud fraction
    qlst_idx, &         		! Physical in-cloud LWC
    qist_idx, &         		! Physical in-cloud IWC
    dp_frac_idx, &      		! deep convection cloud fraction
    sh_frac_idx, &      		! shallow convection cloud fraction
    rel_idx, &             		! Rel
    kvh_idx, &			        ! CLUBB eddy diffusivity on thermo levels
    kvm_idx, &				! CLUBB eddy diffusivity on mom levels
    pblh_idx, &                         ! PBL pbuf
    icwmrdp_idx, &			! In cloud mixing ratio for deep convection
    tke_idx, &                          ! turbulent kinetic energy
    tpert_idx, &                        ! temperature perturbation from PBL
    fice_idx, &                         ! fice_idx index in physics buffer
    cmeliq_idx, &                       ! cmeliq_idx index in physics buffer
    relvar_idx, &                       ! relative cloud water variance
    accre_enhan_idx                     ! optional accretion enhancement factor for MG

  !  Output arrays for CLUBB statistics    
  real(r8), allocatable, dimension(:,:,:) :: out_zt, out_zm, out_radzt, out_radzm, out_sfc

  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90

  contains
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_register_cam( )
!-------------------------------------------------------------------------------
! Description:
!   Register the constituents and fields in the physics buffer
! Author: P. Bogenschutz, C. Craig, A. Gettelman
!
!-------------------------------------------------------------------------------

  end subroutine clubb_register_cam
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_ini_cam(pbuf2d)
!-------------------------------------------------------------------------------
! Description:
!   Initialize UWM CLUBB.
! Author: Cheryl Craig March 2011
! Modifications: Pete Bogenschutz 2011 March and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------




    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, physics_buffer_desc
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    end subroutine clubb_ini_cam
    
    
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    subroutine clubb_tend_cam( &
                              state,   ptend_all,   pbuf,     hdtime, &
                              cmfmc,   cmfmc2,      cam_in,   sgh30,  dlf,   &
			      det_s, det_ice)
			      
!-------------------------------------------------------------------------------
! Description: Provide tendencies of shallow convection, turbulence, and 
!              macrophysics from CLUBB to 1
!   
! Author: Cheryl Craig, March 2011
! Modifications: Pete Bogenschutz, March 2011 and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------

    use physics_types, 	only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init, &
			      physics_ptend_sum, physics_update

    use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              pbuf_set_field, physics_buffer_desc
			      
    use ppgrid, 	only: pver, pverp, pcols
    use constituents, 	only: cnst_get_ind
    use camsrfexch,     only: cam_in_t
      

    implicit none
   
    ! --------------- !
    ! Input Auguments !
    ! --------------- !

    type(physics_state), intent(in)    :: state                    ! Physics state variables  			[vary]
    type(cam_in_t),      intent(in)    :: cam_in
    real(r8), 		 intent(in)    :: hdtime                   ! Host model timestep			[s]
    real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection	[kg/ks/s]
    real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c		[kg/m2/s]
    real(r8),            intent(in)    :: cmfmc2(pcols,pverp)      ! shallow convective mass flux--m subc 	[kg/m2/s]
    real(r8),            intent(in)    :: sgh30(pcols)             ! std deviation of orography			[m]
    
    ! ---------------------- !
    ! Input-Output Auguments !
    ! ---------------------- !
    
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! ---------------------- !
    ! Output Auguments !
    ! ---------------------- !

    type(physics_ptend), intent(out)   :: ptend_all 		           ! package tendencies

    ! These two variables are needed for energy check    
    real(r8),            intent(out)   :: det_s(pcols)               ! Integral of detrained static energy from ice
    real(r8),            intent(out)   :: det_ice(pcols)             ! Integral of detrained ice for energy check

        
    ! --------------- !
    ! Local Variables !
    ! --------------- !

   det_s(:)   = 0.0_r8
   det_ice(:) = 0.0_r8
  end subroutine clubb_tend_cam
    
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
    subroutine clubb_surface ( &
    			     state, ptend, ztodt, cam_in, ustar, obklen)
    
!-------------------------------------------------------------------------------
! Description: Provide the obukov length and the surface friction velocity 
!              for the dry deposition code in routine tphysac.  Since University
!              of Washington Moist Turbulence (UWMT) scheme is not called when 
!              CLUBB is turned on the obukov length and ustar are never initialized
!              nor computed (sometimes never updated from NaN).  In addition, surface
!              fluxes are applied to the constituents.  
!   
! Author: Peter Bogenschutz, August 2011
! Origin: Based heavily on UWMT code (eddy_diff.F90)
! References:
!   None
!-------------------------------------------------------------------------------

    use physics_types,		only: physics_state, physics_ptend, physics_ptend_init
    use physconst, 		only: gravit, zvir, latvap
    use ppgrid, 	        only: pver, pcols
    use constituents, 	        only: pcnst, cnst_get_ind
    use camsrfexch,             only: cam_in_t
    
    implicit none
    
    ! --------------- !
    ! Input Auguments !
    ! --------------- !

    type(physics_state), intent(in)    	:: state		! Physics state variables
    type(cam_in_t),      intent(in)     :: cam_in
    
    real(r8),  		 intent(in)     :: ztodt		! 2 delta-t        [ s ] 

    ! ---------------- !
    ! Output Auguments !
    ! ---------------- !
    
    type(physics_ptend), intent(out)    :: ptend		! Individual parameterization tendencies
    real(r8), 		 intent(out)  	:: obklen(pcols)    	! Obukhov length [ m ]
    real(r8), 		 intent(out)	:: ustar(pcols)		! Surface friction velocity [ m/s ]
    
    obklen(pcols) = 0.0_r8
    ustar(pcols)  = 0.0_r8

    end subroutine clubb_surface


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  
    !-----------------------------------------------------------------------
  subroutine stats_end_timestep_clubb(lchnk,thecol,out_zt,out_zm,out_radzt,out_radzm,out_sfc)

    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------


    integer :: lchnk
    integer :: thecol
    
    real(r8), intent(inout) :: out_zt(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_zm(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_radzt(:,:,:)  ! (pcols,pverp,rad_zt%nn)
    real(r8), intent(inout) :: out_radzm(:,:,:)  ! (pcols,pverp,rad_zm%nn)
    real(r8), intent(inout) :: out_sfc(:,:,:)    ! (pcols,1,sfc%nn)


  end subroutine stats_end_timestep_clubb
  
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  
  
end module clubb_intr
