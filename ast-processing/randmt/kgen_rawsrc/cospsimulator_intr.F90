module cospsimulator_intr
!----------------------------------------------------------------------------------------------------------------------
!Purpose: 1 interface to
!         Name:         CFMIP Observational Simulator Package (COSP)
!         What:         Simulate ISCCP/CloudSat/CALIOP cloud products from GCM inputs
!         Version:      v1.3 released June 2010, updated from v1.1 released May 2009
!         Authors:      Multiple - see http://www.cfmip.net/
!
!Author:  J. Kay (jenkay@ucar.edu) with help from Brian Eaton, John Truesdale, and Y. Zhang/J. Boyle/S. Klein (LLNL)    

! Created: August 2009
! Last modified: Jan 18, 2012
! Status: COSP RUNS WITH BOTH CAM4 and CAM5 physics, CAM5 implementation now includes snow

!! REQUIRED COSP OUTPUTS IF RUNNING COSP OFF-LINE
! If you want to run COSP off-line, the required fields are available and can be added to a history tape via the 1 namelist via fincl calls
! i.e., for CAM4:
! fincl2 = 'Z3:I','Q:I','T:I','PS:I','CLOUD:I','CONCLD:I','CLDICE:I','CLDLIQ:I','LS_FLXPRC:I','LS_FLXSNW:I',
! 'ZMFLXPRC:I','ZMFLXSNW:I','HKFLXPRC:I','HKFLXSNW:I',','REL:I','REI:I','ICLDTWP:I','ICLDIWP:I','EMIS:I'
! i.e., for CAM5: 
! fincl2 = 'Z3:I','Q:I','T:I','PS:I','CLOUD:I','CONCLD:I','CLDICE:I','CLDLIQ:I','LS_FLXPRC:I','LS_FLXSNW:I',
! 'ZMFLXPRC:I','ZMFLXSNW:I','UWFLXPRC:I','UWFLXSNW:I','REL:I','REI:I','ICLDTWP:I','ICLDIWP:I','EMIS:I',
! 'LS_REFFRAIN:I','LS_REFFSNOW:I','CV_REFFLIQ:I','CV_REFFICE:I'
! These can be also set using the namelist variables cosp_histfile_aux (.false.) and cosp_histfile_aux_num (2).

! NOTES from J. Kay on interface design:
! I used ISCCP simulator interface (cloudsimulator_38.F90) as a template.
! Like ISCCP simulator, COSP is called within radiation. F90.
! I have kept the number of changes to COSP core routines to a bare minimum so that it will be easy to add any updates to the code.
! I have also kept this interface as self-contained as possible. e.g., get variables from the physics buffer here
! "Don't pollute the common space."
! I have put "##2check##" next to places in the code where I am unsure what to do or if I have done the right thing.
! Note: These commands start a timer so that you can see how long pieces of the code take to run
!   results in timing file e.g., /ptmp/jenkay/track1_vanilla/ccsm_timing
!   call t_startf ('cospsimulator_intr_run')
!   call t_stopf ('cospsimulator_intr_run')

!----------------------------------------------------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use spmd_utils,      only: masterproc
   use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk 
   use cam_history,     only: addfld, add_default, phys_decomp, outfld
   use cam_history_support,     only: max_fieldname_len 
   use perf_mod,        only: t_startf, t_stopf
   use abortutils,      only: endrun
   use cam_pio_utils,   only: max_chars
   use phys_control,    only: cam_physpkg_is

   implicit none
   private
   save

!Public functions/subroutines

   public :: &
        cospsimulator_intr_readnl,    &
        cospsimulator_intr_init,    &
        cospsimulator_intr_run

   logical, public :: docosp = .false.  ! whether to do COSP calcs and I/O, default is false
                                        ! if docosp is specified in the atm_in namelist,
                                        ! this value is overwritten and cosp is run

   ! frequency at which cosp is called, every cosp_nradsteps radiation timestep
   integer, public :: cosp_nradsteps = 1! 1 namelist variable default, not in COSP namelist

! Private module data

   ! number of dimensions
   integer, parameter ::      &
        nprs_cosp       = 7,  &! number of pressure ranges
        ntau_cosp       = 7,  &! number of optical depth ranges
        ntau_cosp_modis = 6,  &! number of optical depth ranges MODIS
        ndbze_cosp      = 15, &! number of dBZe ranges for COSP radar simulator
        nsr_cosp        = 15, &! number of scattering ranges for COSP lidar simulator
        nhtmisr_cosp    = 16, &! number of heights for misr output (per Marchand)
        nbnds_cosp      = 2,  &! number of bounds of cosp output (for cam_history.F90)
        nsza_cosp       = 5    ! number of solar zenith angle for COSP parasol output
   integer, parameter :: nhtml_cosp = pver  ! number of model levels is pver
   integer    ::  nscol_cosp    ! number of subcolumns for COSP outputs.  
                                ! use namelist input Ncolumns to set
   integer ::  nht_cosp         ! number of height for COSP radar and lidar simulator outputs.  
                                                ! set to 40 if csat_vgrid=.true., else set to Nlr

   ! limits of dimensions, used here to find mid-points, sza_cosp passed to cam_history.F90
   real(r8), parameter :: &
        prslim_cosp_1d(nprs_cosp+1) = (/1000._r8, 800._r8, 680._r8, 560._r8, 440._r8,310._r8,180._r8,0._r8/),  &
        taulim_cosp_1d(ntau_cosp+1) = (/0._r8, 0.3_r8, 1.3_r8, 3.6_r8, 9.4_r8, 23._r8, 60._r8, 379._r8/), &
        taulim_cosp_modis_1d(ntau_cosp_modis+1) = (/0.3_r8, 1.3_r8, 3.6_r8, 9.4_r8, 23._r8, 60._r8, 100000._r8/), &
        dbzelim_cosp_1d(ndbze_cosp+1) = (/-50._r8, -45._r8, -40._r8, -35._r8, -30._r8, -25._r8, -20._r8, -15._r8, &
                -10.0_r8, -5.0_r8, 0.0_r8, 5.0_r8, 10._r8, 15._r8, 20._r8, 25._r8/), &
        srlim_cosp_1d(nsr_cosp+1) = (/0.01_r8, 1.2_r8, 3.0_r8, 5.0_r8, 7.0_r8, 10.0_r8, 15.0_r8, 20.0_r8, 25.0_r8, 30.0_r8, &
                40.0_r8, 50.0_r8, 60.0_r8, 80.0_r8, 999.0_r8, 1009.0_r8/), &
        htmisrlim_cosp_1d(nhtmisr_cosp+1) = (/-99.0_r8, 0.0_r8, 0.5_r8, 1.0_r8, 1.5_r8, 2.0_r8, 2.5_r8, & 
                3.0_r8, 4.0_r8, 5.0_r8, 7.0_r8, 9.0_r8, 11.0_r8, 13.0_r8, 15.0_r8, 17.0_r8, 99.0_r8/)
   ! Removed 'parameter' attribute in order to make it a 'target'
   real(r8), target :: sza_cosp(nsza_cosp) = (/0.0_r8, 15.0_r8, 30.0_r8, 45.0_r8, 60.0_r8/)

   ! limits of dimensions, used by cam_history.F90, 2,nX_cosp  where 1,nX_cosp = min, 2,nX_cosp = max
   real(r8), target :: prslim_cosp(2,nprs_cosp)
   real(r8), target :: taulim_cosp(2,ntau_cosp)
   real(r8), target :: taulim_cosp_modis(2,ntau_cosp_modis)
   real(r8), target :: dbzelim_cosp(2,ndbze_cosp)
   real(r8), target :: srlim_cosp(2,nsr_cosp)
   real(r8), target :: htmisrlim_cosp(2,nhtmisr_cosp)

   ! variable declarations - known sizes

   ! mid points and dimension indices of dimensions
   real(r8), target :: prsmid_cosp(nprs_cosp)            ! pressure midpoints of COSP ISCCP output
   real(r8), target :: taumid_cosp(ntau_cosp)            ! optical depth midpoints of COSP ISCCP output
   real(r8), target :: taumid_cosp_modis(ntau_cosp_modis)! optical depth midpoints of COSP MODIS output
   real(r8), target :: dbzemid_cosp(ndbze_cosp)          ! dbze midpoints of COSP radar output
   real(r8), target :: srmid_cosp(nsr_cosp)              ! sr midpoints of COSP lidar output                                     
   real(r8), target :: htmisrmid_cosp(nhtmisr_cosp)      ! htmisr midpoints of COSP misr simulator output
   real(r8) :: htmlmid_cosp(nhtml_cosp)                  ! model level height midpoints for output

   integer :: prstau_cosp(nprs_cosp*ntau_cosp)           ! ISCCP mixed output dimension index
   integer :: prstau_cosp_modis(nprs_cosp*ntau_cosp_modis)       ! MODIS mixed output dimension index
   integer :: htmisrtau_cosp(nhtmisr_cosp*ntau_cosp)     ! misr mixed output dimension index

   ! real values associated with the collapsed mixed dimensions
   real(r8) :: prstau_prsmid_cosp(nprs_cosp*ntau_cosp)
   real(r8) :: prstau_taumid_cosp(nprs_cosp*ntau_cosp)
   real(r8) :: prstau_prsmid_cosp_modis(nprs_cosp*ntau_cosp_modis)
   real(r8) :: prstau_taumid_cosp_modis(nprs_cosp*ntau_cosp_modis)
   real(r8) :: htmisrtau_htmisrmid_cosp(nhtmisr_cosp*ntau_cosp)
   real(r8) :: htmisrtau_taumid_cosp(nhtmisr_cosp*ntau_cosp)

   ! variable declarations - allocatable sizes
   real(r8),allocatable, target :: htlim_cosp(:,:)       ! height limits for COSP outputs (nht_cosp+1)
   real(r8),allocatable :: htlim_cosp_1d(:)              ! height limits for COSP outputs (nht_cosp+1)
   real(r8),allocatable, target :: htmid_cosp(:)         ! height midpoints of COSP radar/lidar output (nht_cosp)
   integer,allocatable, target :: scol_cosp(:)           ! sub-column number (nscol_cosp)
   integer,allocatable :: htdbze_cosp(:)                 ! radar CFAD mixed output dimension index (nht_cosp*ndbze_cosp)
   integer,allocatable :: htsr_cosp(:)                   ! lidar CFAD mixed output dimension index (nht_cosp*nsr_cosp)
   integer,allocatable :: htmlscol_cosp(:)               ! html-subcolumn mixed output dimension index (nhtml_cosp*nscol_cosp)
   real(r8),allocatable :: htdbze_htmid_cosp(:)          ! (nht_cosp*ndbze_cosp)
   real(r8),public,allocatable :: htdbze_dbzemid_cosp(:)        ! (nht_cosp*ndbze_cosp)
   real(r8),allocatable :: htsr_htmid_cosp(:)            ! (nht_cosp*nsr_cosp)
   real(r8),allocatable :: htsr_srmid_cosp(:)            ! (nht_cosp*nsr_cosp)
   real(r8),allocatable:: htmlscol_htmlmid_cosp(:)       ! (nhtml_cosp*nscol_cosp)
   real(r8),allocatable :: htmlscol_scol_cosp(:)         ! (nhtml_cosp*nscol_cosp)

!! The 1 and COSP namelists defaults are set below.  Some of the COSP namelist 
!! variables are part of the 1 namelist - they all begin with "cosp_" to keep their 
!! names specific to COSP. I set their 1 namelist defaults here, not in namelist_defaults_cam.xml
!!  Variables identified as namelist variables are defined in
!!  ../models/atm/cam/bld/namelist_files/namelist_definition.xml

!! 1 namelist variable defaults
   logical :: cosp_sample_atrain = .false.    ! 1 namelist variable default, not in COSP namelist
   character(len=256) :: cosp_atrainorbitdata   ! 1 namelist variable, no default, need to specify!
   logical :: cosp_amwg = .false.            ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_lite = .false.            ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_passive = .false.         ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_active = .false.          ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_isccp = .false.           ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_3hr = .false.       ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_da = .false.        ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_off = .false.       ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_mon = .false.       ! 1 namelist variable default, not in COSP namelist
   logical :: cosp_lradar_sim = .false.      ! 1 namelist variable default
   logical :: cosp_llidar_sim = .false.      ! 1 namelist variable default
   logical :: cosp_lisccp_sim = .false.      ! 1 namelist variable default
   logical :: cosp_lmisr_sim = .false.       ! 1 namelist variable default
   logical :: cosp_lmodis_sim = .false.      ! 1 namelist variable default
   logical :: cosp_histfile_aux = .false.    ! 1 namelist variable default
   logical :: cosp_lfrac_out = .false.       ! 1 namelist variable default
   logical :: cosp_runall = .false.          ! flag to run all of the cosp simulator package
   integer :: cosp_ncolumns = 50             ! 1 namelist variable default
   integer :: cosp_histfile_num =1           ! 1 namelist variable default, not in COSP namelist 
   integer :: cosp_histfile_aux_num =-1      ! 1 namelist variable default, not in COSP namelist

!! COSP Namelist variables from cosp_output_nl.txt 
   ! Simulator flags
   logical :: lradar_sim = .false.              ! COSP namelist variable, can be changed from default by 1 namelist
   logical :: llidar_sim = .false.              ! ""
   logical :: lisccp_sim = .false.              ! ""
   logical :: lmisr_sim  = .false.              ! ""
   logical :: lmodis_sim = .false.              ! ""
   logical :: lrttov_sim = .false.              ! not running rttov, always set to .false.

   ! Output variables 
   ! All initialized to .false., set to .true. based on 1 namelist in cospsimulator_intr_run
   logical :: lfrac_out = .false.               ! COSP namelist variable, can be changed from default by 1 namelist
   logical :: lalbisccp = .false.
   logical :: latb532 = .false.
   logical :: lboxptopisccp = .false.
   logical :: lboxtauisccp = .false.
   logical :: lcfad_dbze94 = .false.
   logical :: lcfad_lidarsr532 = .false.
   logical :: lclcalipso = .false.
   logical :: lclhcalipso = .false.
   logical :: lclisccp2 = .false.
   logical :: lcllcalipso = .false.
   logical :: lclmcalipso = .false.
   logical :: lcltcalipso = .false.
   logical :: lctpisccp = .false.
   logical :: ldbze94 = .false.
   logical :: lcltradar = .false.
   logical :: lcltradar2 = .false.
   logical :: ltauisccp = .false.
   logical :: ltclisccp = .false.
   logical :: lparasol_refl = .false.
   logical :: lclmisr = .false.
   logical :: lmeantbisccp = .false.
   logical :: lmeantbclrisccp = .false.
   logical :: lclcalipso2 = .false.
   logical :: lcltlidarradar = .false.
   logical :: lbeta_mol532 = .false.
   logical :: Llongitude = .false.
   logical :: Llatitude =.false.
   logical :: lcltmodis = .false.
   logical :: lclwmodis = .false.
   logical :: lclimodis = .false.
   logical :: lclhmodis = .false.
   logical :: lclmmodis = .false.
   logical :: lcllmodis = .false.
   logical :: ltautmodis = .false.
   logical :: ltauwmodis = .false.
   logical :: ltauimodis = .false.
   logical :: ltautlogmodis = .false.
   logical :: ltauwlogmodis = .false.
   logical :: ltauilogmodis = .false.
   logical :: lreffclwmodis = .false.
   logical :: lreffclimodis = .false.
   logical :: lpctmodis = .false.
   logical :: llwpmodis = .false.
   logical :: liwpmodis = .false.
   logical :: lclmodis = .false.
   logical :: ltbrttov = .false.  !! RTTOV OUTPUT (always set to .false.)

! COSP namelist variables from cosp_input_nl.txt
! Set default values, values discussed with Yuying in ()
! Values from cosp_test.f90 case released with cospv1.1 were used as a template
! Note: Unless otherwise specified, these are parameters that cannot be set by the 1 namelist.

   integer, parameter :: Npoints_it = 10000             ! Max # gridpoints to be processed in one iteration (10,000)
   integer :: ncolumns = 50                             ! Number of subcolumns in SCOPS (50), can be changed from default by 1 namelist
   integer :: nlr = 40                                  ! Number of levels in statistical outputs 
                                                        ! (only used if USE_VGRID=.true.)  (40)
   logical :: use_vgrid = .true.                        ! Use fixed vertical grid for outputs? 
                                                        ! (if .true. then define # of levels with nlr)  (.true.)
   logical :: csat_vgrid = .true.                       ! CloudSat vertical grid? 
                                                        ! (if .true. then the CloudSat standard grid is used. 
                                                        ! If set, overides use_vgrid.) (.true.)
   ! namelist variables for COSP input related to radar simulator
   real(r8) :: radar_freq = 94.0_r8                     ! CloudSat radar frequency (GHz) (94.0)
   integer :: surface_radar = 0                         ! surface=1, spaceborne=0 (0)
   integer :: use_mie_tables = 0                        ! use a precomputed lookup table? yes=1,no=0 (0)
   integer :: use_gas_abs = 1                           ! include gaseous absorption? yes=1,no=0 (1)
   integer :: do_ray = 0                                ! calculate/output Rayleigh refl=1, not=0 (0)
   integer :: melt_lay = 0                              ! melting layer model off=0, on=1 (0)
   real(r8) :: k2 = -1                                  ! |K|^2, -1=use frequency dependent default (-1)
   ! namelist variables for COSP input related to lidar simulator
   integer, parameter :: Nprmts_max_hydro = 12          ! Max # params for hydrometeor size distributions (12)
   integer, parameter :: Naero = 1                      ! Number of aerosol species (Not used) (1)
   integer, parameter :: Nprmts_max_aero = 1            ! Max # params for aerosol size distributions (not used) (1)
   integer :: lidar_ice_type = 0                        ! Ice particle shape in lidar calculations 
                                                        ! (0=ice-spheres ; 1=ice-non-spherical) (0)
   integer, parameter :: overlap = 3                    ! overlap type: 1=max, 2=rand, 3=max/rand (3)

   !! namelist variables for COSP input related to ISCCP simulator
   integer :: isccp_topheight = 1                       ! 1 = adjust top height using both a computed infrared 
                                                        ! brightness temperature and the visible
                                                        ! optical depth to adjust cloud top pressure. 
                                                        ! Note that this calculation is most appropriate to compare
                                                        ! to ISCCP data during sunlit hours.
                                                        ! 2 = do not adjust top height, that is cloud top pressure 
                                                        ! is the actual cloud top pressure in the model
                                                        ! 3 = adjust top height using only the computed infrared 
                                                        ! brightness temperature. Note that this calculation is most 
                                                        ! appropriate to compare to ISCCP IR only algortihm (i.e. 
                                                        ! you can compare to nighttime ISCCP data with this option) (1)
   integer :: isccp_topheight_direction = 2             ! direction for finding atmosphere pressure level with 
                                                        ! interpolated temperature equal to the radiance
                                                        ! determined cloud-top temperature
                                                        ! 1 = find the *lowest* altitude (highest pressure) level 
                                                        ! with interpolated temperature 
                                                        ! equal to the radiance determined cloud-top temperature
                                                        ! 2 = find the *highest* altitude (lowest pressure) level 
                                                        ! with interpolated temperature
                                                        ! equal to the radiance determined cloud-top temperature
                                                        ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                                        ! 1 = default setting in COSP v1.1, matches all versions of 
                                                        ! ISCCP simulator with versions numbers 3.5.1 and lower
                                                        ! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator
   ! RTTOV inputs (not used)
   integer, parameter :: Platform = 1                   ! satellite platform (1)
   integer, parameter :: Satellite = 15                 ! satellite (15)
   integer, parameter :: Instrument = 0                 ! instrument (0)
   integer, parameter :: Nchannels = 8                  ! Number of channels to be computed (8)
   integer, parameter :: Channels(Nchannels) =  (/1,3,5,6,8,10,11,13/)  
                                                        ! Channel numbers (match and supply Nchannels) 
                                                        ! (1,3,5,6,8,10,11,13,)
   real(r8), parameter :: Surfem(Nchannels) = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/)               
                                                        ! Surface emissivity (match and supply Nchannels)
                                                        ! (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,)
   real(r8), parameter :: ZenAng = 50._r8               ! Satellite Zenith Angle (50)
   real(r8), parameter :: co =  2.098e-07_r8            ! Mixing ratio CO (2.098e-07), not used in cospv1.1
                                                        ! set to value from cosp_test.F90's cosp_input_nl.txt
   !per Alejandro: cosp rttov gaseous inputs are also mass mixing ratios
   !values from cosp_test.F90, units (kg/kg)
   !Mixing ratio C02 (5.241e-04), Mixing ratio CH4 (9.139e-07), 
   !Mixing ratio N20 (4.665e-07), Mixing ratio CO (2.098e-07)
   !I get CO2, CH4, N20 from cam radiation interface.

!! Other variables
    integer,parameter :: nhydro = 9                     ! number of COSP hydrometeor classes
    logical,allocatable :: first_run_cosp(:)            !.true. if run_cosp has been populated (allocatable->begchunk:endchunk)
    logical,allocatable :: run_cosp(:,:)                !.true. if cosp should be run by column and chunk (allocatable->1:pcols,begchunk:endchunk)

!!! Variables read in from atrain orbit file (private module data)
    integer, parameter :: norbitdata = 9866324
    real(r8), allocatable :: &
        atrainlat(:),         & !  A-train orbit latitudes (float in original data file
        atrainlon(:)            !  A-train orbit longitudes
    integer, allocatable :: &
        atrainday(:),        &  !  A-train calendar day (short in data file)
        atrainhr(:),         &  !  A-train hour (byte in data file)
        atrainmin(:),        &  !  A-train minute (byte in data file)
        atrainsec(:)            !  A-train minute (byte in data file)
    integer :: idxas = 0        ! index to start loop over atrain orbit


! pbuf indices
    integer :: cld_idx, concld_idx, lsreffrain_idx, lsreffsnow_idx, cvreffliq_idx
    integer :: cvreffice_idx, dpcldliq_idx, dpcldice_idx
    integer :: shcldliq_idx, shcldice_idx, shcldliq1_idx, shcldice1_idx, dpflxprc_idx
    integer :: dpflxsnw_idx, shflxprc_idx, shflxsnw_idx, lsflxprc_idx, lsflxsnw_idx
    integer:: rei_idx, rel_idx



CONTAINS

subroutine setcospvalues(Nlr_in,use_vgrid_in,csat_vgrid_in,Ncolumns_in,docosp_in,cosp_nradsteps_in)

   ! input arguments
   integer, intent(in) :: Nlr_in
   logical, intent(in) :: use_vgrid_in
   logical, intent(in) :: csat_vgrid_in
   integer, intent(in) :: Ncolumns_in
   logical, intent(in) :: docosp_in
   integer, intent(in) :: cosp_nradsteps_in

   ! Local variables
   integer :: i,k		! indices
   real(r8) :: zstep

   ! set vertical grid, reference code from cosp_types.F90, line 549
   ! used to set vgrid_bounds in cosp_io.f90, line 844
   if (use_vgrid_in) then		!! using fixed vertical grid
   	if (csat_vgrid_in) then
     	   nht_cosp = 40
     	   zstep = 480.0_r8
   	else
     	   nht_cosp = Nlr_in
     	   zstep = 20000.0_r8/Nlr_in  ! constant vertical spacing, top at 20 km
   	end if
   end if

!  if (use_vgrid_in=.false.) then    !using the model vertical height grid
!	nht_cosp = pver
!	htlim_cosp = (/0._r8/) ##2check##
!  end if

   ! set number of sub-columns using namelist input
   nscol_cosp=Ncolumns_in

   cosp_nradsteps = cosp_nradsteps_in
  
   ! need to allocate memory for these variables
   allocate(htlim_cosp(2,nht_cosp),htlim_cosp_1d(nht_cosp+1),htmid_cosp(nht_cosp),scol_cosp(nscol_cosp),&
	htdbze_cosp(nht_cosp*ndbze_cosp),htsr_cosp(nht_cosp*nsr_cosp),htmlscol_cosp(nhtml_cosp*nscol_cosp),&
	htdbze_htmid_cosp(nht_cosp*ndbze_cosp),htdbze_dbzemid_cosp(nht_cosp*ndbze_cosp),&
	htsr_htmid_cosp(nht_cosp*nsr_cosp),htsr_srmid_cosp(nht_cosp*nsr_cosp),&
	htmlscol_htmlmid_cosp(nhtml_cosp*nscol_cosp),htmlscol_scol_cosp(nhtml_cosp*nscol_cosp))

   if (use_vgrid_in) then		!! using fixed vertical grid
      htlim_cosp_1d(1)= 0.0_r8
      do i=2,nht_cosp+1
         htlim_cosp_1d(i)=(i-1)*zstep      !! based on cosp_types.F90 line 556
      enddo
   end if

   ! calculate mid-points and bounds for cam_history.F90

   do k=1,nprs_cosp
      prsmid_cosp(k) = 0.5_r8*(prslim_cosp_1d(k) + prslim_cosp_1d(k+1))
      prslim_cosp(1,k) = prslim_cosp_1d(k)
      prslim_cosp(2,k) = prslim_cosp_1d(k+1)
   end do

   do k=1,ntau_cosp
      taumid_cosp(k) = 0.5_r8*(taulim_cosp_1d(k) + taulim_cosp_1d(k+1))
      taulim_cosp(1,k) = taulim_cosp_1d(k)
      taulim_cosp(2,k) = taulim_cosp_1d(k+1)
   end do

   do k=1,ntau_cosp_modis
      taumid_cosp_modis(k) = 0.5_r8*(taulim_cosp_modis_1d(k) + taulim_cosp_modis_1d(k+1))
      taulim_cosp_modis(1,k) = taulim_cosp_modis_1d(k)
      taulim_cosp_modis(2,k) = taulim_cosp_modis_1d(k+1)
   end do

   do k=1,ndbze_cosp
      dbzemid_cosp(k) = 0.5_r8*(dbzelim_cosp_1d(k) + dbzelim_cosp_1d(k+1))
      dbzelim_cosp(1,k) = dbzelim_cosp_1d(k)
      dbzelim_cosp(2,k) = dbzelim_cosp_1d(k+1)
   end do

   do k=1,nsr_cosp
      srmid_cosp(k) = 0.5_r8*(srlim_cosp_1d(k) + srlim_cosp_1d(k+1))
      srlim_cosp(1,k) = srlim_cosp_1d(k)
      srlim_cosp(2,k) = srlim_cosp_1d(k+1)
   end do

   htmisrmid_cosp(1) = -99.0_r8
   htmisrlim_cosp(1,1) = htmisrlim_cosp_1d(1)
   htmisrlim_cosp(2,1) = htmisrlim_cosp_1d(2)
   do k=2,nhtmisr_cosp
      htmisrmid_cosp(k) = 0.5_r8*(htmisrlim_cosp_1d(k) + htmisrlim_cosp_1d(k+1))
      htmisrlim_cosp(1,k) = htmisrlim_cosp_1d(k)
      htmisrlim_cosp(2,k) = htmisrlim_cosp_1d(k+1)
   end do

   do k=1,nht_cosp
      htmid_cosp(k) = 0.5_r8*(htlim_cosp_1d(k) + htlim_cosp_1d(k+1))
      htlim_cosp(1,k) = htlim_cosp_1d(k)
      htlim_cosp(2,k) = htlim_cosp_1d(k+1)
   end do

   do k=1,nscol_cosp
      scol_cosp(k) = k
   end do

   !  Just using an index here, model height is a prognostic variable
   do k=1,nhtml_cosp
      htmlmid_cosp(k) = k
   end do

   ! assign mixed dimensions an integer index for cam_history.F90
   do k=1,nprs_cosp*ntau_cosp
      prstau_cosp(k) = k
   end do
   do k=1,nprs_cosp*ntau_cosp_modis
      prstau_cosp_modis(k) = k
   end do
   do k=1,nht_cosp*ndbze_cosp
      htdbze_cosp(k) = k
   end do
   do k=1,nht_cosp*nsr_cosp
      htsr_cosp(k) = k
   end do
   do k=1,nhtml_cosp*nscol_cosp
      htmlscol_cosp(k) = k
   end do
   do k=1,nhtmisr_cosp*ntau_cosp
      htmisrtau_cosp(k) = k
   end do

   ! next, assign collapsed reference vectors for cam_history.F90
   ! convention for saving output = prs1,tau1 ... prs1,tau7 ; prs2,tau1 ... prs2,tau7 etc.
   ! actual output is specified in cospsimulator_intr.F90
   do k=1,nprs_cosp
      prstau_taumid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=taumid_cosp(1:ntau_cosp)
      prstau_prsmid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=prsmid_cosp(k)
      prstau_taumid_cosp_modis(ntau_cosp_modis*(k-1)+1:k*ntau_cosp_modis)=taumid_cosp_modis(1:ntau_cosp_modis)
      prstau_prsmid_cosp_modis(ntau_cosp_modis*(k-1)+1:k*ntau_cosp_modis)=prsmid_cosp(k)
   enddo

   do k=1,nht_cosp
      htdbze_dbzemid_cosp(ndbze_cosp*(k-1)+1:k*ndbze_cosp)=dbzemid_cosp(1:ndbze_cosp)
      htdbze_htmid_cosp(ndbze_cosp*(k-1)+1:k*ndbze_cosp)=htmid_cosp(k)
   enddo

   do k=1,nht_cosp
      htsr_srmid_cosp(nsr_cosp*(k-1)+1:k*nsr_cosp)=srmid_cosp(1:nsr_cosp)
      htsr_htmid_cosp(nsr_cosp*(k-1)+1:k*nsr_cosp)=htmid_cosp(k)
   enddo

   do k=1,nhtml_cosp
      htmlscol_scol_cosp(nscol_cosp*(k-1)+1:k*nscol_cosp)=scol_cosp(1:nscol_cosp)
      htmlscol_htmlmid_cosp(nscol_cosp*(k-1)+1:k*nscol_cosp)=htmlmid_cosp(k)
   enddo

   do k=1,nhtmisr_cosp
      htmisrtau_taumid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=taumid_cosp(1:ntau_cosp)
      htmisrtau_htmisrmid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=htmisrmid_cosp(k)
   enddo

end subroutine setcospvalues

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! subroutine to read namelist variables and run setcospvalues subroutine
! note: cldfrc_readnl is a good template in cloud_fraction.F90
! Make sure that this routine is reading in a namelist.  
! models/atm/cam/bld/build-namelist is the perl script to check

subroutine cospsimulator_intr_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand,    only: mpicom, mpilog, mpiint, mpichar

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input  (nlfile=atm_in)

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cospsimulator_intr_readnl'

    !!! this list should include any variable that you might want to include in the namelist
    !!! philosophy is to not include COSP output flags but just important COSP settings and cfmip controls. 
    namelist /cospsimulator_nl/ docosp, cosp_active, cosp_amwg, cosp_atrainorbitdata, cosp_cfmip_3hr, cosp_cfmip_da, &
        cosp_cfmip_mon, cosp_cfmip_off, cosp_histfile_num, cosp_histfile_aux, cosp_histfile_aux_num, cosp_isccp, cosp_lfrac_out, &
        cosp_lite, cosp_lradar_sim, cosp_llidar_sim, cosp_lisccp_sim,  cosp_lmisr_sim, cosp_lmodis_sim, cosp_ncolumns, &
        cosp_nradsteps, cosp_passive, cosp_sample_atrain, cosp_runall
   !-----------------------------------------------------------------------------

   !! read in the namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )  !! presumably opens the namelist file "nlfile"
      !! position the file to write to the cospsimulator portion of the cam_in namelist
      call find_group_name(unitn, 'cospsimulator_nl', status=ierr)   
      if (ierr == 0) then
           read(unitn, cospsimulator_nl, iostat=ierr)
           if (ierr /= 0) then
              call endrun(subname // ':: ERROR reading namelist')
           end if
      end if
           close(unitn)
           call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpibcast(docosp,               1,  mpilog, 0, mpicom)
   call mpibcast(cosp_atrainorbitdata, len(cosp_atrainorbitdata), mpichar, 0, mpicom)
   call mpibcast(cosp_amwg,            1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lite,            1,  mpilog, 0, mpicom)
   call mpibcast(cosp_passive,         1,  mpilog, 0, mpicom)
   call mpibcast(cosp_active,          1,  mpilog, 0, mpicom)
   call mpibcast(cosp_isccp,           1,  mpilog, 0, mpicom)
   call mpibcast(cosp_runall,          1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_3hr,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_da,        1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_mon,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_off,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lfrac_out,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lradar_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_llidar_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lisccp_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lmisr_sim,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lmodis_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_ncolumns,        1,  mpiint, 0, mpicom)
   call mpibcast(cosp_sample_atrain,   1,  mpilog, 0, mpicom)
   call mpibcast(cosp_histfile_num,    1,  mpiint, 0, mpicom)
   call mpibcast(cosp_histfile_aux_num,1,  mpiint, 0, mpicom)
   call mpibcast(cosp_histfile_aux,    1,  mpilog, 0, mpicom)
   call mpibcast(cosp_nradsteps,       1,  mpiint, 0, mpicom)

   !! reset COSP namelist variables based on input from cam namelist variables
   if (cosp_cfmip_3hr) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
   end if
   if (cosp_cfmip_da) then
      llidar_sim = .true.
      lisccp_sim = .true.
   end if
   if (cosp_cfmip_off) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
   end if
   if (cosp_cfmip_mon) then
      llidar_sim = .true.
      lisccp_sim = .true.
   end if

   if (cosp_lfrac_out) then
      lfrac_out = .true.
   end if
   if (cosp_lradar_sim) then
      lradar_sim = .true.
   end if
   if (cosp_llidar_sim) then
      llidar_sim = .true.
   end if
   if (cosp_lisccp_sim) then
      lisccp_sim = .true.
   end if
   if (cosp_lmisr_sim) then
      lmisr_sim = .true.
   end if
   if (cosp_lmodis_sim) then
      lmodis_sim = .true.
   end if

   if (cosp_histfile_aux .and. cosp_histfile_aux_num == -1) then
      cosp_histfile_aux_num = cosp_histfile_num
   end if

   if (cosp_lite) then
      llidar_sim = .true.
      lisccp_sim = .true.
      lmisr_sim = .true.
      lmodis_sim = .true.
      cosp_ncolumns = 10
      cosp_nradsteps = 3
   end if

   if (cosp_passive) then
      lisccp_sim = .true.
      lmisr_sim = .true.
      lmodis_sim = .true.
      cosp_ncolumns = 10
      cosp_nradsteps = 3
   end if

   if (cosp_active) then
      lradar_sim = .true.
      llidar_sim = .true.
      cosp_ncolumns = 10
      cosp_nradsteps = 3
   end if

   if (cosp_isccp) then
      lisccp_sim = .true.
      cosp_ncolumns = 10
      cosp_nradsteps = 3
   end if

   if (cosp_runall) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
      lmisr_sim = .true.
      lmodis_sim = .true.
      lfrac_out = .true.
   end if

   !! if no simulators are turned on at all and docosp is, set cosp_amwg = .true.
   if((docosp) .and. (.not.lradar_sim) .and. (.not.llidar_sim) .and. (.not.lisccp_sim) .and. &
     (.not.lmisr_sim) .and. (.not.lmodis_sim)) then
      cosp_amwg = .true.
   end if
   if (cosp_amwg) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
      lmisr_sim = .true.
      lmodis_sim = .true.
      cosp_ncolumns = 10
      cosp_nradsteps = 3
   end if


   !! reset COSP namelist variables based on input from cam namelist variables
   if (cosp_ncolumns .ne. ncolumns) then
      ncolumns = cosp_ncolumns
   end if

   !! use the namelist variables to overwrite default COSP output variables above (set to .false.)
   if (lradar_sim) then
      lcfad_dbze94 = .true.
      ldbze94 = .true.
      lcltradar = .true.
      lcltradar2 = .true.
   end if
   if ((lradar_sim) .and. (llidar_sim)) then
      !! turn on the outputs that require both the radar and the lidar simulator
      lclcalipso2 = .true.
      lcltlidarradar = .true.
   end if

   if (llidar_sim) then
      !! turn on the outputs for the lidar simulator
      lcllcalipso = .true.
      lclmcalipso = .true.
      lcltcalipso = .true.
      lclcalipso = .true.
      lclhcalipso = .true.
      lcfad_lidarsr532 = .true.
      latb532 = .true.
      lparasol_refl = .true.
      lbeta_mol532 = .true.
   end if
   if (lisccp_sim) then
      !! turn on the outputs for the isccp simulator
      lalbisccp = .true.
      lboxptopisccp = .true.
      lboxtauisccp = .true.
      lclisccp2 = .true.
      lctpisccp = .true.
      ltauisccp = .true.
      ltclisccp = .true.
      lmeantbisccp = .true.
      lmeantbclrisccp = .true.
   end if
   if (lmisr_sim) then
      !! turn on the outputs for the misr simulator
      lclmisr = .true.
   end if
   if (lmodis_sim) then
      !! turn on the outputs for the modis simulator
      lcltmodis = .true.
      lclwmodis = .true.
      lclimodis = .true.
      lclhmodis = .true.
      lclmmodis = .true.
      lcllmodis = .true.
      ltautmodis = .true.
      ltauwmodis = .true.
      ltauimodis = .true.
      ltautlogmodis = .true.
      ltauwlogmodis = .true.
      ltauilogmodis = .true.
      lreffclwmodis = .true.
      lreffclimodis = .true.
      lpctmodis = .true.
      llwpmodis = .true.
      liwpmodis = .true.
      lclmodis = .true.
   end if

   ! Set vertical coordinate, subcolumn, and calculation frequency cosp options based on namelist inputs
   call setcospvalues(nlr,use_vgrid,csat_vgrid,ncolumns,docosp,cosp_nradsteps)

end subroutine cospsimulator_intr_readnl

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine cospsimulator_intr_init

   use cam_history,         only: addfld, add_default, phys_decomp
   use mpishorthand,        only : mpir8, mpiint, mpicom
   use netcdf,              only : nf90_open, nf90_inq_varid, nf90_get_var, nf90_close, nf90_nowrite
   use error_messages,      only : handle_ncerr, alloc_err
   use cam_history_support, only: add_hist_coord
   
   use physics_buffer,  only: pbuf_get_index
   real(r8),parameter :: R_UNDEF = -1.0E30_r8
   integer ncid,latid,lonid,did,hrid,minid,secid, istat
   !------------------------------------------------------------------------------

if (cosp_sample_atrain) then

!!!! READ IN ATRAIN ORBIT DATA FROM INPUT FILE FOR SUB-SAMPLING
   allocate(atrainlat(norbitdata),atrainlon(norbitdata),atrainday(norbitdata),&
        atrainhr(norbitdata),atrainmin(norbitdata),atrainsec(norbitdata),stat=istat) 
   call alloc_err( istat, 'cospsimulator_intr', 'atrain', norbitdata )
        
   if (masterproc) then
      call handle_ncerr( nf90_open (cosp_atrainorbitdata,NF90_NOWRITE,ncid),&
        'cospsimulator_intr.F90 cospsimulator_intr_init')
      call handle_ncerr( nf90_inq_varid (ncid,'lat',latid),&
        'cospsimulator_intr.F90 latid')
      call handle_ncerr( nf90_inq_varid (ncid,'lon',lonid),&
        'cospsimulator_intr.F90 lonid')
      call handle_ncerr( nf90_inq_varid (ncid,'doy365',did),&
        'cospsimulator_intr.F90 did')
      call handle_ncerr( nf90_inq_varid (ncid,'hour',hrid),&
        'cospsimulator_intr.F90 hrid')
      call handle_ncerr( nf90_inq_varid (ncid,'minute',minid),&
        'cospsimulator_intr.F90 minid')
      call handle_ncerr( nf90_inq_varid (ncid,'second',secid),&
        'cospsimulator_intr.F90 secid')
      call handle_ncerr( nf90_get_var (ncid,latid,atrainlat),&
        'cospsimulator_intr.F90 atrainlat')
      call handle_ncerr( nf90_get_var (ncid,lonid,atrainlon),&
        'cospsimulator_intr.F90 atrainlon')
      call handle_ncerr( nf90_get_var (ncid,did,atrainday),&
        'cospsimulator_intr.F90 atrainday')
      call handle_ncerr( nf90_get_var (ncid,hrid,atrainhr),&
        'cospsimulator_intr.F90 atrainhr')
      call handle_ncerr( nf90_get_var (ncid,minid,atrainmin),&
        'cospsimulator_intr.F90 atrainmin')
      call handle_ncerr( nf90_get_var (ncid,secid,atrainsec),&
        'cospsimulator_intr.F90 atrainsec')
      call handle_ncerr( nf90_close (ncid),&
        'cospsimulator_intr.F90 nf90_close')
   end if
   call mpibcast (atrainlat,norbitdata , mpir8, 0, mpicom)  !! 2 arg = size, see runtime_opts.F90
   call mpibcast (atrainlon,norbitdata , mpir8, 0, mpicom)
   call mpibcast (atrainday,norbitdata , mpiint, 0, mpicom)
   call mpibcast (atrainhr,norbitdata , mpiint, 0, mpicom)
   call mpibcast (atrainmin,norbitdata , mpiint, 0, mpicom)
   call mpibcast (atrainsec,norbitdata , mpiint, 0, mpicom)

endif

! ADDFLD ADD_DEFAULT CALLS FOR COSP OUTPUTS
! notes on addfld/add_default/outfld calls:  
! 1) Dimensions of cosp output should be:  
!       a) ntime=1 (cosp run at every time step)
!       b) nprofile=ncol
! 2) See cam_history.F90, modified for new cosp output dimensions
! 3) Need conditional logic so that addfld,add_default,and outfld calls are only executed when fields are available.
! 4) nhtml_cosp=height_mlev=pver (height of model levels)
! 5) flag_xyfill=.true. is "non-applicable xy points flagged with fillvalue".  
!  per brian, flag_xyfill should be set to true whenever the fields might contain fillvalues.
!  I think this should be .true. for all COSP outputs.
!  Especially because there are cases where there will definitely be fill values (e.g., tau not calculated when no cloud.)
!  For cosp variables with subcolumns (dimension includes nscol_cosp) I have made the outputs instantaneous by default
!  to get around check_accum failing and aborting run in cam_history.90.  Problem is that the vertical dimension
!  can contain a mix of fillvalue and realvalue (e.g., when there is a cloud in one sub-column but not in another).
!  none of these variables are a part of CFMIP.  Also needed to modify cam_history so that the check_accum is
!  not done when output is instantaneous.
! 6) sampling_seq = radiation timestep.  note: cloudsimulator.F90 does not specify anything special.
! 7) Markings for CFMIP output requirements:
!*cfMon* = CFMIP variable for cfMon - monthly-mean cloud diagnostic fields
!*cfOff* = CFMIP variable for cfOff - monthly-mean offline cloud diagnostic fields
!*cfDa* = CFMIP variable for cfDa - daily-mean cloud diagnostic fields
!*cf3hr* = CFMIP variable for cf3hr - 3-hourly cloud diagnostic fields
! 8) Making it the default to produce a separate 1 history file with COSP outputs.  
! "2" is the history tape index.  I specify 2 here to create a separate history file for cosp output. 
! 9) crash fix: add_default was looking through the master list for input field name and not finding it.
! Solution? I removed all of the spaces in the addfld calls, add_default, and outfld calls.

!!! ISCCP OUTPUTS
   if (lisccp_sim) then
      ! register non-standard variable dimensions
      call add_hist_coord('cosp_prs', nprs_cosp, 'COSP Mean ISCCP pressure',  &
           'hPa', prsmid_cosp, bounds_name='cosp_prs_bnds', bounds=prslim_cosp)
      call add_hist_coord('cosp_tau', ntau_cosp,                              &
           'COSP Mean ISCCP optical depth', '1', taumid_cosp,          &
           bounds_name='cosp_tau_bnds', bounds=taulim_cosp)
      call add_hist_coord('cosp_scol', nscol_cosp, 'COSP subcolumn',          &
           values=scol_cosp)

      !! addfld calls for all
      !*cfMon,cfDa* clisccp2 (time,tau,plev,profile), CFMIP wants 7 p bins, 7 tau bins
      call addfld('FISCCP1_COSP','percent   ',nprs_cosp*ntau_cosp,'A', &
                   'Grid-box fraction covered by each ISCCP D level cloud type',phys_decomp,&
                   flag_xyfill=.true., mdimnames=(/'cosp_tau','cosp_prs'/), fill_value=R_UNDEF)

      !*cfMon,cfDa* tclisccp (time,profile), CFMIP wants "gridbox mean cloud cover from ISCCP"
      call addfld('CLDTOT_ISCCP','percent', 1,'A', &
                   'Total Cloud Fraction Calculated by the ISCCP Simulator ',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      !*cfMon,cfDa* albisccp (time,profile)
      !!! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
      call addfld('MEANCLDALB_ISCCP','1',1,'A','Mean cloud albedo*CLDTOT_ISCCP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      !*cfMon,cfDa* ctpisccp (time,profile)
      !!! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)     
      call addfld('MEANPTOP_ISCCP','Pa',1,'A','Mean cloud top pressure*CLDTOT_ISCCP',phys_decomp,flag_xyfill=.true., &
                  fill_value=R_UNDEF)
      ! tauisccp (time,profile)
      !!! For averaging, weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
      call addfld ('MEANTAU_ISCCP','1',1,'A','Mean optical thickness*CLDTOT_ISCCP',phys_decomp,flag_xyfill=.true., &
                   fill_value=R_UNDEF)
      ! meantbisccp (time,profile), at 10.5 um
      call addfld ('MEANTB_ISCCP','K       ',1,'A','Mean Infrared Tb from ISCCP simulator',phys_decomp,flag_xyfill=.true., &
                   fill_value=R_UNDEF)
      ! meantbclrisccp (time,profile)
      call addfld ('MEANTBCLR_ISCCP','K       ',1,'A','Mean Clear-sky Infrared Tb from ISCCP simulator',&
        phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! boxtauisccp (time,column,profile)
      call addfld ('TAU_ISCCP','1',nscol_cosp,'I','Optical Depth in each Subcolumn',&
        phys_decomp,flag_xyfill=.true.,mdimnames=(/'cosp_scol'/), fill_value=R_UNDEF)
      ! boxptopisccp (time,column,profile)
      call addfld ('CLDPTOP_ISCCP','Pa',nscol_cosp,'I','Cloud Top Pressure in each Subcolumn',&
        phys_decomp,flag_xyfill=.true.,mdimnames=(/'cosp_scol'/), fill_value=R_UNDEF)

      !!! add_default calls for CFMIP experiments or else all fields are added to history file except those with sub-column dimension
      if (cosp_cfmip_mon.or.cosp_cfmip_da) then
           !! add cfmip-requested variables to two separate cam history files
         if (cosp_cfmip_da) then
            call add_default ('FISCCP1_COSP',2,' ')
            call add_default ('CLDTOT_ISCCP',2,' ')
            call add_default ('MEANCLDALB_ISCCP',2,' ')
            call add_default ('MEANPTOP_ISCCP',2,' ')
         end if
         if (cosp_cfmip_mon) then
            call add_default ('FISCCP1_COSP',1,' ')
            call add_default ('CLDTOT_ISCCP',1,' ')
            call add_default ('MEANCLDALB_ISCCP',1,' ')
            call add_default ('MEANPTOP_ISCCP',1,' ')
         end if
      else
           !! add all isccp outputs to the history file specified by the 1 namelist variable cosp_histfile_num
            call add_default ('FISCCP1_COSP',cosp_histfile_num,' ')
            call add_default ('CLDTOT_ISCCP',cosp_histfile_num,' ')
            call add_default ('MEANCLDALB_ISCCP',cosp_histfile_num,' ')
            call add_default ('MEANPTOP_ISCCP',cosp_histfile_num,' ')
            call add_default ('MEANTAU_ISCCP',cosp_histfile_num,' ')
            call add_default ('MEANTB_ISCCP',cosp_histfile_num,' ')
            call add_default ('MEANTBCLR_ISCCP',cosp_histfile_num,' ')
      end if
   end if

!!! LIDAR SIMULATOR OUTPUTS
   if (llidar_sim) then
      call add_hist_coord('cosp_ht', nht_cosp,                                &
           'COSP Mean Height for lidar and radar simulator outputs', 'm',     &
           htmid_cosp, bounds_name='cosp_ht_bnds', bounds=htlim_cosp)
      call add_hist_coord('cosp_sr', nsr_cosp,                                &
           'COSP Mean Scattering Ratio for lidar simulator CFAD output', '1', &
           srmid_cosp, bounds_name='cosp_sr_bnds', bounds=srlim_cosp)
      call add_hist_coord('cosp_sza', nsza_cosp, 'COSP Parasol SZA',          &
           'degrees', sza_cosp)
      call add_hist_coord('cosp_scol', nscol_cosp, 'COSP subcolumn',          &
           values=scol_cosp)

      !! addfld calls for all
      !*cfMon,cfOff,cfDa,cf3hr* cllcalipso (time,profile)
      call addfld('CLDLOW_CAL','percent',1,'A','Lidar Low-level Cloud Fraction',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      !*cfMon,cfOff,cfDa,cf3hr* clmcalipso (time,profile)
      call addfld('CLDMED_CAL','percent',1,'A','Lidar Mid-level Cloud Fraction',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      !*cfMon,cfOff,cfDa,cf3hr* clhcalipso (time,profile)
      call addfld('CLDHGH_CAL','percent',1,'A','Lidar High-level Cloud Fraction',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      !*cfMon,cfOff,cfDa,cf3hr* cltcalipso (time,profile)
      call addfld('CLDTOT_CAL','percent',1,'A','Lidar Total Cloud Fraction',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      !*cfMon,cfOff,cfDa,cf3hr* clcalipso (time,height,profile)
      call addfld('CLD_CAL','percent',nht_cosp,'A','Lidar Cloud Fraction (532 nm)',&
                   phys_decomp, flag_xyfill=.true., mdimnames=(/'cosp_ht'/), fill_value=R_UNDEF)
      !*cfMon,cfOff,cfDa,cf3hr* parasol_refl (time,sza,profile)
      call addfld ('RFL_PARASOL','fraction',nsza_cosp,'A','PARASOL-like mono-directional reflectance ',&
                   phys_decomp,flag_xyfill=.true.,mdimnames=(/'cosp_sza'/), fill_value=R_UNDEF)
      !*cfOff,cf3hr* cfad_lidarsr532 (time,height,scat_ratio,profile), %11%, default is 40 vert levs, 15 SR  bins
      call addfld('CFAD_SR532_CAL','fraction',nht_cosp*nsr_cosp,'A',&
                   'Lidar Scattering Ratio CFAD (532 nm)',phys_decomp,&
                   flag_xyfill=.true., mdimnames=(/'cosp_sr','cosp_ht'/), fill_value=R_UNDEF)
      ! beta_mol532 (time,height_mlev,profile)
      call addfld ('MOL532_CAL','m-1sr-1',nhtml_cosp,'A','Lidar Molecular Backscatter (532 nm) ',&
                      phys_decomp,flag_xyfill=.true.,mdimnames=(/'lev'/), fill_value=R_UNDEF)
      ! atb532 (time,height_mlev,column,profile)
      call addfld ('ATB532_CAL','no_unit_log10(x)',nhtml_cosp*nscol_cosp,'I', &
                      'Lidar Attenuated Total Backscatter (532 nm) in each Subcolumn',phys_decomp, &
                      flag_xyfill=.true.,mdimnames=(/'cosp_scol','lev      '/), fill_value=R_UNDEF)

      !!! add_default calls for CFMIP experiments or else all fields are added to history file except those with sub-column dimension/experimental variables
      if (cosp_cfmip_mon .or. cosp_cfmip_off .or. cosp_cfmip_da .or. cosp_cfmip_3hr) then
         if (cosp_cfmip_da) then
            call add_default ('CLDLOW_CAL',2,' ')
            call add_default ('CLDMED_CAL',2,' ')
            call add_default ('CLDHGH_CAL',2,' ')
            call add_default ('CLDTOT_CAL',2,' ')
            call add_default ('CLD_CAL',2,' ')
            call add_default ('RFL_PARASOL',2,' ')
         end if
         if (cosp_cfmip_mon.or.cosp_cfmip_off) then
            call add_default ('CLDLOW_CAL',1,' ')
            call add_default ('CLDMED_CAL',1,' ')
            call add_default ('CLDHGH_CAL',1,' ')
            call add_default ('CLDTOT_CAL',1,' ')
            call add_default ('CLD_CAL',1,' ')
            call add_default ('RFL_PARASOL',1,' ')
         end if
         if (cosp_cfmip_3hr) then
            call add_default ('CFAD_SR532_CAL',3,' ')
            call add_default ('CLDLOW_CAL',3,' ')
            call add_default ('CLDMED_CAL',3,' ')
            call add_default ('CLDHGH_CAL',3,' ')
            call add_default ('CLDTOT_CAL',3,' ')
            call add_default ('CLD_CAL',3,' ')
            call add_default ('RFL_PARASOL',3,' ')
         end if
         if (cosp_cfmip_off) then
            call add_default ('CFAD_SR532_CAL',1,' ')
         end if
      else
         !! add all lidar outputs to the history file specified by the 1 namelist variable cosp_histfile_num
         call add_default ('CLDLOW_CAL',cosp_histfile_num,' ')
         call add_default ('CLDMED_CAL',cosp_histfile_num,' ')
         call add_default ('CLDHGH_CAL',cosp_histfile_num,' ')
         call add_default ('CLDTOT_CAL',cosp_histfile_num,' ')
         call add_default ('CLD_CAL',cosp_histfile_num,' ')
         call add_default ('RFL_PARASOL',cosp_histfile_num,' ')
         call add_default ('CFAD_SR532_CAL',cosp_histfile_num,' ')

         if((.not.cosp_amwg) .and. (.not.cosp_lite) .and. (.not.cosp_passive) .and. (.not.cosp_active) .and. (.not.cosp_isccp)) then
            call add_default ('MOL532_CAL',cosp_histfile_num,' ')
         end if

       end if

   end if

!!! RADAR SIMULATOR OUTPUTS
   if (lradar_sim) then
      call add_hist_coord('cosp_ht', nht_cosp,                                &
           'COSP Mean Height for lidar and radar simulator outputs', 'm',     &
           htmid_cosp, bounds_name='cosp_ht_bnds', bounds=htlim_cosp)
      call add_hist_coord('cosp_dbze', ndbze_cosp,                            &
           'COSP Mean dBZe for radar simulator CFAD output', 'dBZ',           &
           dbzemid_cosp, bounds_name='cosp_dbze_bnds', bounds=dbzelim_cosp)
      call add_hist_coord('cosp_scol', nscol_cosp, 'COSP subcolumn',          &
           values=scol_cosp)

      !!! addfld calls
      !*cfOff,cf3hr* cfad_dbze94 (time,height,dbze,profile), default is 40 vert levs, 15 dBZ bins 
      call addfld('CFAD_DBZE94_CS','fraction',nht_cosp*ndbze_cosp,'A',&
                   'Radar Reflectivity Factor CFAD (94 GHz)',phys_decomp,&
                   flag_xyfill=.true., mdimnames=(/'cosp_dbze','cosp_ht  '/), fill_value=R_UNDEF)
      !*cfOff,cf3hr* clcalipso2 (time,height,profile)
      call addfld ('CLD_CAL_NOTCS','percent',nht_cosp,'A','Cloud occurrence seen by CALIPSO but not CloudSat ',&
                   phys_decomp,flag_xyfill=.true.,mdimnames=(/'cosp_ht'/), fill_value=R_UNDEF)
      ! cltlidarradar (time,profile)
      call addfld ('CLDTOT_CALCS','percent',1,'A',' Lidar and Radar Total Cloud Fraction ',phys_decomp,flag_xyfill=.true., &
                   fill_value=R_UNDEF)
      call addfld ('CLDTOT_CS','percent',1,'A',' Radar total cloud amount ',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('CLDTOT_CS2','percent',1,'A',' Radar total cloud amount without the data for the first kilometer above surface '&
                   ,phys_decomp, flag_xyfill=.true., fill_value=R_UNDEF)
      ! dbze94 (time,height_mlev,column,profile),! height_mlevel = height when vgrid_in = .true. (default)
      call addfld ('DBZE_CS','dBZe    ',nhtml_cosp*nscol_cosp,'I',' Radar dBZe (94 GHz) in each Subcolumn',phys_decomp,&
                      flag_xyfill=.true.,mdimnames=(/'cosp_scol','lev      '/), fill_value=R_UNDEF)

      !!! add_default calls for CFMIP experiments or else all fields are added to history file except those with sub-column dimension
       if (cosp_cfmip_off.or.cosp_cfmip_3hr) then
          if (cosp_cfmip_3hr) then
              call add_default ('CFAD_DBZE94_CS',3,' ')
              call add_default ('CLD_CAL_NOTCS',3,' ')
          end if
          if (cosp_cfmip_off) then
              call add_default ('CFAD_DBZE94_CS',1,' ')
              call add_default ('CLD_CAL_NOTCS',1,' ')
          end if
      else
         !! add all radar outputs to the history file specified by the 1 namelist variable cosp_histfile_num
          call add_default ('CFAD_DBZE94_CS',cosp_histfile_num,' ')
          call add_default ('CLD_CAL_NOTCS',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CALCS',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CS',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CS2',cosp_histfile_num,' ')
      end if
   end if

!!! MISR SIMULATOR OUTPUTS
   if (lmisr_sim) then
      call add_hist_coord('cosp_htmisr', nhtmisr_cosp, 'COSP MISR height',    &
           'km', htmisrmid_cosp,                                              &
           bounds_name='cosp_htmisr_bnds', bounds=htmisrlim_cosp)
      call add_hist_coord('cosp_tau', ntau_cosp,                              &
           'COSP Mean ISCCP optical depth', '1', taumid_cosp,          &
           bounds_name='cosp_tau_bnds', bounds=taulim_cosp)
      call add_hist_coord('cosp_scol', nscol_cosp, 'COSP subcolumn',          &
           values=scol_cosp)

      ! clMISR (time,tau,CTH_height_bin,profile)
      call addfld ('CLD_MISR','percent',nhtmisr_cosp*ntau_cosp,'A','Cloud Fraction from MISR Simulator',&
                   phys_decomp,flag_xyfill=.true.,mdimnames=(/'cosp_tau   ','cosp_htmisr'/), fill_value=R_UNDEF)
      !! add all misr outputs to the history file specified by the 1 namelist variable cosp_histfile_num
      call add_default ('CLD_MISR',cosp_histfile_num,' ')
   end if

!!! MODIS OUTPUT
   if (lmodis_sim) then

      call add_hist_coord('cosp_prs', nprs_cosp, 'COSP Mean ISCCP pressure',  &
           'hPa', prsmid_cosp, bounds_name='cosp_prs_bnds', bounds=prslim_cosp)
      call add_hist_coord('cosp_tau_modis', ntau_cosp_modis,                  &
           'COSP Mean MODIS optical depth', '1', taumid_cosp_modis,    &
           bounds_name='cosp_tau_modis_bnds', bounds=taulim_cosp_modis)

      ! float cltmodis ( time, loc )
      call addfld ('CLTMODIS','%',1,'A','MODIS Total Cloud Fraction',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float clwmodis ( time, loc )
      call addfld ('CLWMODIS','%',1,'A','MODIS Liquid Cloud Fraction',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float climodis ( time, loc )
      call addfld ('CLIMODIS','%',1,'A','MODIS Ice Cloud Fraction',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float clhmodis ( time, loc )
      call addfld ('CLHMODIS','%',1,'A','MODIS High Level Cloud Fraction',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float clmmodis ( time, loc )
      call addfld ('CLMMODIS','%',1,'A','MODIS Mid Level Cloud Fraction',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float cllmodis ( time, loc )
      call addfld ('CLLMODIS','%',1,'A','MODIS Low Level Cloud Fraction',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float tautmodis ( time, loc )
      call addfld ('TAUTMODIS','1',1,'A','MODIS Total Cloud Optical Thickness*CLTMODIS',& 
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float tauwmodis ( time, loc )
      call addfld ('TAUWMODIS','1',1,'A','MODIS Liquid Cloud Optical Thickness*CLWMODIS',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float tauimodis ( time, loc )
      call addfld ('TAUIMODIS','1',1,'A','MODIS Ice Cloud Optical Thickness*CLIMODIS',& 
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float tautlogmodis ( time, loc )
      call addfld ('TAUTLOGMODIS','1',1,'A','MODIS Total Cloud Optical Thickness (Log10 Mean)*CLTMODIS',&  
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float tauwlogmodis ( time, loc )
      call addfld ('TAUWLOGMODIS','1',1,'A','MODIS Liquid Cloud Optical Thickness (Log10 Mean)*CLWMODIS',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float tauilogmodis ( time, loc )
      call addfld ('TAUILOGMODIS','1',1,'A','MODIS Ice Cloud Optical Thickness (Log10 Mean)*CLIMODIS',&  
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float reffclwmodis ( time, loc )
      call addfld ('REFFCLWMODIS','m',1,'A','MODIS Liquid Cloud Particle Size*CLWMODIS',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float reffclimodis ( time, loc )
      call addfld ('REFFCLIMODIS','m',1,'A','MODIS Ice Cloud Particle Size*CLIMODIS',& 
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float pctmodis ( time, loc )
      call addfld ('PCTMODIS','Pa',1,'A','MODIS Cloud Top Pressure*CLTMODIS',&       
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float lwpmodis ( time, loc )
      call addfld ('LWPMODIS','kg m-2',1,'A','MODIS Cloud Liquid Water Path*CLWMODIS',&
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float iwpmodis ( time, loc )
      call addfld ('IWPMODIS','kg m-2',1,'A','MODIS Cloud Ice Water Path*CLIMODIS',&   
                   phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      ! float clmodis ( time, plev, tau, loc )
      call addfld ('CLMODIS','%',nprs_cosp*ntau_cosp_modis,'A','MODIS Cloud Area Fraction',&
                   phys_decomp,flag_xyfill=.true., mdimnames=(/'cosp_tau_modis','cosp_prs      '/), fill_value=R_UNDEF)

      !! add MODIS output to history file specified by the 1 namelist variable cosp_histfile_num
      call add_default ('CLTMODIS',cosp_histfile_num,' ')
      call add_default ('CLWMODIS',cosp_histfile_num,' ')
      call add_default ('CLIMODIS',cosp_histfile_num,' ')
      call add_default ('CLHMODIS',cosp_histfile_num,' ')
      call add_default ('CLMMODIS',cosp_histfile_num,' ')
      call add_default ('CLLMODIS',cosp_histfile_num,' ')
      call add_default ('TAUTMODIS',cosp_histfile_num,' ')
      call add_default ('TAUWMODIS',cosp_histfile_num,' ')
      call add_default ('TAUIMODIS',cosp_histfile_num,' ')
      call add_default ('TAUTLOGMODIS',cosp_histfile_num,' ')
      call add_default ('TAUWLOGMODIS',cosp_histfile_num,' ')
      call add_default ('TAUILOGMODIS',cosp_histfile_num,' ')
      call add_default ('REFFCLWMODIS',cosp_histfile_num,' ')
      call add_default ('REFFCLIMODIS',cosp_histfile_num,' ')
      call add_default ('PCTMODIS',cosp_histfile_num,' ')
      call add_default ('LWPMODIS',cosp_histfile_num,' ')
      call add_default ('IWPMODIS',cosp_histfile_num,' ')
      call add_default ('CLMODIS',cosp_histfile_num,' ')
   end if

!!! SUB-COLUMN OUTPUT
   if (lfrac_out) then
      ! frac_out (time,height_mlev,column,profile)
      call addfld ('SCOPS_OUT','0=nocld,1=strcld,2=cnvcld',nhtml_cosp*nscol_cosp,'I','SCOPS Subcolumn output',&
                   phys_decomp,flag_xyfill=.true.,mdimnames=(/'cosp_scol','lev      '/), fill_value=R_UNDEF)
      !! add scops ouptut to history file specified by the 1 namelist variable cosp_histfile_num
      call add_default ('SCOPS_OUT',cosp_histfile_num,' ')
      ! save sub-column outputs from ISCCP if ISCCP is run
      if (lisccp_sim) then
         call add_default ('TAU_ISCCP',cosp_histfile_num,' ')
         call add_default ('CLDPTOP_ISCCP',cosp_histfile_num,' ')
      end if
      ! save sub-column outputs from lidar if lidar is run
      if (llidar_sim) then
         call add_default ('ATB532_CAL',cosp_histfile_num,' ')
      end if
      ! save sub-column outputs from radar if radar is run
      if (lradar_sim) then
         call add_default ('DBZE_CS',cosp_histfile_num,' ')
      end if
   end if 

!! ADDFLD, ADD_DEFAULT, OUTFLD CALLS FOR COSP OUTPUTS IF RUNNING COSP OFF-LINE
!! Note: A suggestion was to add all of the 1 variables needed to add to make it possible to run COSP off-line
!! These fields are available and can be called from the namelist though.  Here, when the cosp_runall mode is invoked
!! all of the inputs are saved on the cam history file.  This is good de-bugging functionality we should maintain.

   if (cosp_histfile_aux) then
      call addfld ('PS_COSP','Pa',1,'I','PS_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('TS_COSP','K',1,'I','TS_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('P_COSP','Pa',pver,'I','P_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('PH_COSP','Pa',pver,'I','PH_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('ZLEV_COSP','m',pver,'I','ZLEV_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('ZLEV_HALF_COSP','m',pver,'I','ZLEV_HALF_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('T_COSP','K',pver,'I','T_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('RH_COSP','percent',pver,'I','RH_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('Q_COSP','kg/kg',pver,'I','Q_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('CONCLD_COSP','1',pver,'I','CONCLD_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('CLD_COSP','1',pver,'I','CLD_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('O3_COSP','kg/kg',pver,'I','O3_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('U_COSP','m/s',1,'I','U_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)  
      call addfld ('V_COSP','m/s',1,'I','V_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)  
      call addfld ('LSCLIQ_COSP','kg/kg',pver,'I','LSCLIQ_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('LSCICE_COSP','kg/kg',pver,'I','LSCICE_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('CVCLIQ_COSP','kg/kg',pver,'I','CVCLIQ_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('CVCICE_COSP','kg/kg',pver,'I','CVCICE_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('RAIN_LS_COSP','kg/m2/s',pver,'I','RAIN_LS_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('SNOW_LS_COSP','kg/m2/s',pver,'I','SNOW_LS_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('GRPL_LS_COSP','kg/m2/s',pver,'I','GRPL_LS_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('RAIN_CV_COSP','kg/m2/s',pver,'I','RAIN_CV_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('SNOW_CV_COSP','kg/m2/s',pver,'I','SNOW_CV_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_1','m',pver,'I','REFF_COSP_1',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_2','m',pver,'I','REFF_COSP_2',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_3','m',pver,'I','REFF_COSP_3',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_4','m',pver,'I','REFF_COSP_4',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_5','m',pver,'I','REFF_COSP_5',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_6','m',pver,'I','REFF_COSP_6',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_7','m',pver,'I','REFF_COSP_7',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_8','m',pver,'I','REFF_COSP_8',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('REFF_COSP_9','m',pver,'I','REFF_COSP_9',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('DTAU_S_COSP','1',pver,'I','DTAU_S_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('DTAU_C_COSP','1',pver,'I','DTAU_C_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('DEM_S_COSP','1',pver,'I','DEM_S_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('DEM_C_COSP','1',pver,'I','DEM_C_COSP',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('DTAU_S_COSP_SNOW','1',pver,'I','DTAU_S_COSP_SNOW',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)
      call addfld ('DEM_S_COSP_SNOW','1',pver,'I','DEM_S_COSP_SNOW',phys_decomp,flag_xyfill=.true., fill_value=R_UNDEF)

      call add_default ('PS_COSP',cosp_histfile_aux_num,' ')
      call add_default ('TS_COSP',cosp_histfile_aux_num,' ')
      call add_default ('P_COSP',cosp_histfile_aux_num,' ')
      call add_default ('PH_COSP',cosp_histfile_aux_num,' ')
      call add_default ('ZLEV_COSP',cosp_histfile_aux_num,' ')
      call add_default ('ZLEV_HALF_COSP',cosp_histfile_aux_num,' ')
      call add_default ('T_COSP',cosp_histfile_aux_num,' ')
      call add_default ('RH_COSP',cosp_histfile_aux_num,' ')
      call add_default ('Q_COSP',cosp_histfile_aux_num,' ')
      call add_default ('CONCLD_COSP',cosp_histfile_aux_num,' ')
      call add_default ('CLD_COSP',cosp_histfile_aux_num,' ')
      call add_default ('O3_COSP',cosp_histfile_aux_num,' ')
      call add_default ('U_COSP',cosp_histfile_aux_num,' ')
      call add_default ('V_COSP',cosp_histfile_aux_num,' ')
      call add_default ('LSCLIQ_COSP',cosp_histfile_aux_num,' ')
      call add_default ('LSCICE_COSP',cosp_histfile_aux_num,' ')
      call add_default ('CVCLIQ_COSP',cosp_histfile_aux_num,' ')
      call add_default ('CVCICE_COSP',cosp_histfile_aux_num,' ')
      call add_default ('RAIN_LS_COSP',cosp_histfile_aux_num,' ')
      call add_default ('SNOW_LS_COSP',cosp_histfile_aux_num,' ')
      call add_default ('GRPL_LS_COSP',cosp_histfile_aux_num,' ')
      call add_default ('RAIN_CV_COSP',cosp_histfile_aux_num,' ')
      call add_default ('SNOW_CV_COSP',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_1',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_2',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_3',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_4',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_5',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_6',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_7',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_8',cosp_histfile_aux_num,' ')
      call add_default ('REFF_COSP_9',cosp_histfile_aux_num,' ')
      call add_default ('DTAU_S_COSP',cosp_histfile_aux_num,' ')
      call add_default ('DTAU_C_COSP',cosp_histfile_aux_num,' ')
      call add_default ('DEM_S_COSP',cosp_histfile_aux_num,' ')
      call add_default ('DEM_C_COSP',cosp_histfile_aux_num,' ')
      call add_default ('DTAU_S_COSP_SNOW',cosp_histfile_aux_num,' ')
      call add_default ('DEM_S_COSP_SNOW',cosp_histfile_aux_num,' ')
   end if

   rei_idx = pbuf_get_index('REI')
   rel_idx = pbuf_get_index('REL')
   cld_idx = pbuf_get_index('CLD')
   concld_idx = pbuf_get_index('CONCLD')

   if( cam_physpkg_is('cam5') ) then
      lsreffrain_idx = pbuf_get_index('LS_REFFRAIN')
      lsreffsnow_idx = pbuf_get_index('LS_REFFSNOW')
      cvreffliq_idx  = pbuf_get_index('CV_REFFLIQ')
      cvreffice_idx  = pbuf_get_index('CV_REFFICE')
   end if

   dpcldliq_idx  = pbuf_get_index('DP_CLDLIQ')
   dpcldice_idx  = pbuf_get_index('DP_CLDICE')
   shcldliq_idx  = pbuf_get_index('SH_CLDLIQ')
   shcldice_idx  = pbuf_get_index('SH_CLDICE')
   shcldliq1_idx = pbuf_get_index('SH_CLDLIQ1')
   shcldice1_idx = pbuf_get_index('SH_CLDICE1')
   dpflxprc_idx  = pbuf_get_index('DP_FLXPRC')
   dpflxsnw_idx  = pbuf_get_index('DP_FLXSNW')
   shflxprc_idx  = pbuf_get_index('SH_FLXPRC')
   shflxsnw_idx  = pbuf_get_index('SH_FLXSNW')
   lsflxprc_idx  = pbuf_get_index('LS_FLXPRC')
   lsflxsnw_idx  = pbuf_get_index('LS_FLXSNW')

    allocate(first_run_cosp(begchunk:endchunk))
    first_run_cosp(begchunk:endchunk)=.true.
    allocate(run_cosp(1:pcols,begchunk:endchunk))
    run_cosp(1:pcols,begchunk:endchunk)=.false.

end subroutine cospsimulator_intr_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine cospsimulator_intr_run(state,pbuf, cam_in,emis,coszrs,cliqwp_in,cicewp_in,cld_swtau_in,snow_tau_in,snow_emis_in)    
   
   use physics_types,    only: physics_state
   
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
   use camsrfexch,       only: cam_in_t
   use constituents,     only: cnst_get_ind
   use rad_constituents, only: rad_cnst_get_gas
   use wv_saturation,    only: qsat_water
   use phys_control,     only: phys_getopts
   use interpolate_data, only: lininterp_init,lininterp,lininterp_finish,interp_type    
   use physconst,        only: pi, gravit
   use cam_history,      only: outfld,hist_fld_col_active 
   use cmparray_mod,     only: CmpDayNite, ExpDayNite
   use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
   use time_manager,     only: get_curr_calday,get_curr_time,get_ref_date

   real(r8),parameter :: R_UNDEF = -1.0E30_r8
! Arguments
   type(physics_state), intent(in), target :: state
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_in_t),  intent(in) :: cam_in
   !! vars calculated in subroutine param_cldoptics_calc within param_cldoptics.F90
   real(r8), intent(in) :: emis(pcols,pver)             ! cloud longwave emissivity
   real(r8), intent(in) :: coszrs(pcols)                ! cosine solar zenith angle (to tell if day or night)
! make the input arguments optional because they are used differently in CAM4 and CAM5.
   real(r8), intent(in),optional :: cliqwp_in(pcols,pver)   ! in-cloud liquid water path, CAMRT uses this to calculate cld_swtau
   real(r8), intent(in),optional :: cicewp_in(pcols,pver)    ! in-cloud ice water path CAMRT uses this to calculate cld_swtau
   real(r8), intent(in),optional :: cld_swtau_in(pcols,pver) ! RRTM cld_swtau_in, read in using this variable
   real(r8), intent(in),optional :: snow_tau_in(pcols,pver)  ! RRTM grid-box mean SW snow optical depth, used for CAM5 simulations 
   real(r8), intent(in),optional :: snow_emis_in(pcols,pver) ! RRTM grid-box mean LW snow optical depth, used for CAM5 simulations 

! USE_COSP

end subroutine cospsimulator_intr_run

!#######################################################################

end module cospsimulator_intr
