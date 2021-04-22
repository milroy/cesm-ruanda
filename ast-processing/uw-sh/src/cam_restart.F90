module cam_restart
!----------------------------------------------------------------------- 
! 
! module to handle reading and writing of the master restart files.
!
!----------------------------------------------------------------------- 
   use shr_kind_mod,     only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use spmd_utils,       only: masterproc
   use ppgrid,           only: begchunk, endchunk
   use pmgrid,           only: plev, plevp, plat
   use rgrid,            only: nlon, wnummax, fullgrid
   use ioFileMod,        only: getfil, opnfil
   use abortutils,       only: endrun
   use camsrfexch,       only: cam_in_t, cam_out_t     
   use dyn_comp,         only: dyn_import_t, dyn_export_t


   use mpishorthand,     only: mpicom, mpir8, mpiint, mpilog

   use units,            only: getunit
   use shr_kind_mod,     only: shr_kind_cs
   use cam_logfile,      only: iulog
   use pio,              only: file_desc_t, pio_global, pio_noerr, &
                               pio_seterrorhandling, pio_bcast_error, pio_internal_error, &
                               pio_inq_att, pio_def_dim, pio_enddef, &
                               pio_get_att, pio_put_att, pio_closefile

   implicit none
   private
   save

   ! Public interfaces

   public restart_defaultopts    ! initialize namelist variables
   public restart_setopts        ! set values of namelist variables
   public restart_printopts      ! print module options to log
   public cam_write_restart      ! Write the master restart file out
   public cam_read_restart       ! Read the master restart file in
   public get_restcase           ! Get the caseid of the restart file being read in
   public get_restartdir         ! Get the directory name of the restart file being read in

   ! Private data

   integer, parameter :: uninit_int = -999999999

   integer, parameter :: nlen = 256       ! Length of character strings
   character(len=nlen):: pname = ' '      ! Full restart pathname
   character(shr_kind_cs) :: tcase = ' '  ! Read in previous case name

   ! Type of restart run
   logical :: nlres                       ! true => restart or branch run
   logical :: lbrnch                      ! true => branch run

   ! Filename specifiers for master restart filename
   ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = number)
   character(len=nlen) :: rfilename_spec = '%c.cam.r.%y-%m-%d-%s.nc'

   logical :: aeres                ! true => write absorptivities/emissivities to restart file
   integer :: nsds = -1            ! Logical unit number for restart pointer file

   ! Filenames used for restart or branch
   character(len=nlen) :: rest_pfile = './rpointer.atm' ! Restart pointer file contains name of most recently
                                                        ! written restart file
   character(len=nlen) :: cam_branch_file = ' '         ! Filepath of primary restart file for a branch run

!-----------------------------------------------------------------------

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine restart_defaultopts( &
   cam_branch_file_out            )
  !----------------------------------------------------------------------- 
  ! Purpose: Return default runtime options
  !-----------------------------------------------------------------------
  
  character(len=nlen), intent(out), optional :: cam_branch_file_out
  !-----------------------------------------------------------------------

  if ( present(cam_branch_file_out) ) then
     cam_branch_file_out = cam_branch_file
  endif
  
end subroutine restart_defaultopts

!================================================================================================

subroutine restart_setopts( nsrest, cam_branch_file_in )
  !----------------------------------------------------------------------- 
  ! Purpose: Set runtime options
  !-----------------------------------------------------------------------
  use cam_instance, only: inst_suffix
  
  integer,             intent(in)           :: nsrest
  character(len=nlen), intent(in), optional :: cam_branch_file_in
  
  integer :: numset=0
  integer :: rcode
  !-----------------------------------------------------------------------

  ! Set pointer file name based on instance suffix
  rest_pfile = trim(rest_pfile) // trim(inst_suffix)

  rfilename_spec = '%c.cam' // trim(inst_suffix) //'.r.%y-%m-%d-%s.nc'

  ! Set continuation run flags
  if (nsrest==0) then
     nlres  = .false.
     lbrnch = .false.
  else if (nsrest==1) then
     nlres  = .true.
     lbrnch = .false.
  else if (nsrest==3) then
     nlres  = .true.
     lbrnch = .true.
  else
     call endrun ('restart_setopts: not a valid option for nsrest (should be 0, 1 or 3)')
  endif
  
  if ( present(cam_branch_file_in) ) then
     cam_branch_file = cam_branch_file_in
  endif
  
  ! If branch set restart filepath to path given on namelist
  if ( lbrnch ) call set_restart_filepath( cam_branch_file )
  
end subroutine restart_setopts

!=========================================================================================

subroutine restart_printopts

   write(iulog,*)'Summary of restart module options:'
   write(iulog,*)'  Restart pointer file is: ',trim(rest_pfile)
   if (lbrnch) then
      write(iulog,*)'  Branch run will start from: ',trim(cam_branch_file)
   end if
end subroutine restart_printopts

!=========================================================================================

! Subprogram not used    subroutine cam_write_restart( cam_in, cam_out, dyn_out, pbuf2d, &
! Subprogram not used 	                         yr_spec, mon_spec, day_spec, sec_spec )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Write the primary, secondary, and history buffer regeneration files.
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! The cpp 1 definition provides for the funnelling of all program i/o
! Subprogram not used ! through the master processor. Processor 0 either reads restart/history
! Subprogram not used ! data from the disk and distributes it to all processors, or collects
! Subprogram not used ! data from all processors and writes it to disk.
! Subprogram not used ! 
! Subprogram not used ! Author: 
! Subprogram not used ! 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used      use physics_buffer,            only: physics_buffer_desc
! Subprogram not used       use cam_history,      only: write_restart_history, init_restart_history
! Subprogram not used       use radiation,        only: radiation_do
! Subprogram not used       use time_manager,     only: timemgr_write_restart
! Subprogram not used       use filenames,        only: caseid, interpret_filename_spec
! Subprogram not used       use dycore,           only: dycore_is
! Subprogram not used 
! Subprogram not used       use time_manager,     only: timemgr_write_restart, timemgr_init_restart
! Subprogram not used       use restart_dynamics, only: write_restart_dynamics, init_restart_dynamics
! Subprogram not used       use restart_physics,  only: write_restart_physics, init_restart_physics
! Subprogram not used       use cam_pio_utils,    only: cam_pio_createfile
! Subprogram not used       use spmd_utils,       only: iam, mpicom
! Subprogram not used       !
! Subprogram not used       ! Arguments
! Subprogram not used       !
! Subprogram not used       type(cam_in_t),      intent(in) :: cam_in(begchunk:endchunk)
! Subprogram not used       type(cam_out_t),     intent(in) :: cam_out(begchunk:endchunk)
! Subprogram not used       
! Subprogram not used       type(dyn_export_t),  intent(in) :: dyn_out
! Subprogram not used       
! Subprogram not used       type(physics_buffer_desc), pointer  :: pbuf2d(:,:)
! Subprogram not used 
! Subprogram not used       integer            , intent(in), optional :: yr_spec         ! Simulation year
! Subprogram not used       integer            , intent(in), optional :: mon_spec        ! Simulation month
! Subprogram not used       integer            , intent(in), optional :: day_spec        ! Simulation day
! Subprogram not used       integer            , intent(in), optional :: sec_spec        ! Seconds into current simulation day
! Subprogram not used       !
! Subprogram not used       ! Local workspace
! Subprogram not used       !
! Subprogram not used       integer ioerr                 ! write error status
! Subprogram not used       character(len=nlen) :: fname  ! Restart filename
! Subprogram not used       integer :: aeres_int = 0
! Subprogram not used       integer :: ierr
! Subprogram not used       type(file_desc_t) :: File
! Subprogram not used       integer, pointer :: hdimids(:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used       ! Write the primary restart datasets
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       aeres = radiation_do('aeres')
! Subprogram not used       if ( aeres ) aeres_int = 1
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used       ! Write the master restart dataset
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (present(yr_spec).and.present(mon_spec).and.present(day_spec).and.present(sec_spec)) then
! Subprogram not used          fname = interpret_filename_spec( rfilename_spec, &
! Subprogram not used               yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
! Subprogram not used       else
! Subprogram not used          fname = interpret_filename_spec( rfilename_spec )
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       call cam_pio_createfile(File, trim(fname), 0)
! Subprogram not used       call timemgr_init_restart(File)
! Subprogram not used       call init_restart_dynamics(File, hdimids, dyn_out)
! Subprogram not used       call init_restart_physics(File, pbuf2d, hdimids)
! Subprogram not used       call init_restart_history(File)
! Subprogram not used       deallocate(hdimids)
! Subprogram not used 
! Subprogram not used       ierr = PIO_Put_att(File, PIO_GLOBAL, 'caseid', caseid)
! Subprogram not used       ierr = PIO_Put_att(File, PIO_GLOBAL, 'aeres', aeres_int)
! Subprogram not used       if(.not.fullgrid) then
! Subprogram not used          ierr = PIO_Put_att(File, PIO_GLOBAL, 'NLON', nlon)
! Subprogram not used          ierr = PIO_Put_att(File, PIO_GLOBAL, 'WNUMMAX', wnummax)
! Subprogram not used       end if
! Subprogram not used       ierr = pio_enddef(File)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used       ! Dynamics, physics, History
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used       call timemgr_write_restart(File)
! Subprogram not used       call write_restart_dynamics(File, dyn_out)
! Subprogram not used       call write_restart_physics(File, cam_in, cam_out, pbuf2d)
! Subprogram not used 
! Subprogram not used       if (present(yr_spec).and.present(mon_spec).and.&
! Subprogram not used            present(day_spec).and.present(sec_spec)) then
! Subprogram not used          call write_restart_history ( File, &
! Subprogram not used               yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
! Subprogram not used       else
! Subprogram not used          call write_restart_history( File )
! Subprogram not used       end if
! Subprogram not used       call pio_closefile(File)
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used       ! Close the master restart file
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (masterproc) then
! Subprogram not used 	 pname = fname
! Subprogram not used          call write_rest_pfile()
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    end subroutine cam_write_restart

!#######################################################################

! Subprogram not used    subroutine cam_read_restart(cam_in, cam_out, dyn_in, dyn_out, pbuf2d, stop_ymd, stop_tod, NLFileName )
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Acquire and position the restart, master, primary and secondary
! Subprogram not used ! datasets for a continuation run
! Subprogram not used ! 
! Subprogram not used ! Method: 
! Subprogram not used ! 
! Subprogram not used ! Author: 
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used      use physics_buffer, only: physics_buffer_desc
! Subprogram not used       use restart_physics,  only: read_restart_physics
! Subprogram not used       use restart_dynamics, only: read_restart_dynamics
! Subprogram not used       use chem_surfvals,    only: chem_surfvals_init
! Subprogram not used       use phys_grid,        only: phys_grid_init
! Subprogram not used       use camsrfexch,       only: atm2hub_alloc, hub2atm_alloc
! Subprogram not used 
! Subprogram not used       use spmd_dyn,         only: spmdbuf
! Subprogram not used 
! Subprogram not used       use cam_history,      only: read_restart_history
! Subprogram not used       use dycore,           only: dycore_is
! Subprogram not used 
! Subprogram not used       use cam_pio_utils,    only: cam_pio_openfile, clean_iodesc_list
! Subprogram not used       use spmd_utils,       only: iam, mpicom
! Subprogram not used       use time_manager,     only: timemgr_read_restart, timemgr_restart
! Subprogram not used       use filenames,        only: caseid, brnch_retain_casename
! Subprogram not used       use ref_pres,         only: ref_pres_init
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Arguments
! Subprogram not used !
! Subprogram not used    type(cam_in_t),     pointer     :: cam_in(:)
! Subprogram not used    type(cam_out_t),    pointer     :: cam_out(:)
! Subprogram not used    type(dyn_import_t), intent(inout) :: dyn_in
! Subprogram not used    type(dyn_export_t), intent(inout) :: dyn_out
! Subprogram not used    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
! Subprogram not used    character(len=*),   intent(in)  :: NLFileName
! Subprogram not used    integer,            intent(IN)  :: stop_ymd       ! Stop date (YYYYMMDD)
! Subprogram not used    integer,            intent(IN)  :: stop_tod       ! Stop time of day (sec)
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used    character(len=nlen) :: locfn          ! Local filename
! Subprogram not used    character(len=nlen+40) :: errstr
! Subprogram not used    real(r8) :: tmp_rgrid(plat)
! Subprogram not used    integer :: ierr, aeres_int, slen, xtype
! Subprogram not used    type(file_desc_t) :: File
! Subprogram not used    logical :: filefound
! Subprogram not used 
! Subprogram not used       ! lbrnch is false for a restart run (nsrest=1), and true for a 
! Subprogram not used       ! branch run (nsrest=3).  Only read the restart pointer file for
! Subprogram not used       ! a restart run.
! Subprogram not used    aeres = .false.
! Subprogram not used    if (.not.lbrnch) then
! Subprogram not used       call read_rest_pfile
! Subprogram not used    endif
! Subprogram not used   
! Subprogram not used !
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used ! Obtain and read the master restart dataset 
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! Obtain master restart dataset
! Subprogram not used !
! Subprogram not used    call getfil (pname, locfn)
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used    ! Master restart dataset
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used    inquire(FILE=trim(locfn), exist=filefound)
! Subprogram not used 
! Subprogram not used    if(.not.filefound) then
! Subprogram not used       write(errstr,*) 'Could not find restart file ', trim(locfn)
! Subprogram not used       call endrun(errstr)
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    call cam_pio_openfile(File, trim(locfn), 0)
! Subprogram not used    ierr = pio_inq_att(File, pio_global, 'caseid', xtype, slen)
! Subprogram not used    ierr = PIO_Get_att(File, PIO_GLOBAL, 'caseid', tcase)
! Subprogram not used    tcase(slen+1:len(tcase))=''
! Subprogram not used    ierr = PIO_Get_att(File, PIO_GLOBAL, 'aeres', aeres_int)
! Subprogram not used    if(aeres_int==1) aeres=.true.
! Subprogram not used       
! Subprogram not used ! If these variables are not in the file the call will give an error 
! Subprogram not used !       which we can safely ignore
! Subprogram not used    call PIO_SetErrorHandling(File, PIO_BCAST_ERROR)
! Subprogram not used    ierr = pio_get_att(File, PIO_GLOBAL, 'NLON', tmp_rgrid)
! Subprogram not used    if(ierr==PIO_NOERR) nlon=tmp_rgrid
! Subprogram not used    ierr = pio_get_att(File, PIO_GLOBAL, 'WNUMMAX', tmp_rgrid)
! Subprogram not used    if(ierr==PIO_NOERR) wnummax=tmp_rgrid
! Subprogram not used    call PIO_SetErrorHandling(File, PIO_INTERNAL_ERROR)
! Subprogram not used    call timemgr_read_restart(File)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    if(masterproc) then
! Subprogram not used 
! Subprogram not used       if (lbrnch .and. tcase==caseid .and. .not.brnch_retain_casename) then
! Subprogram not used          write(iulog,*) 'READ_RESTART_MASTER: Must change case name on branch run'
! Subprogram not used          write(iulog,*) 'Prev case = ',tcase,' current case = ',caseid
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used       write(iulog,*) 'Files for restart:'
! Subprogram not used 
! Subprogram not used    endif  ! end of if-masterproc
! Subprogram not used 
! Subprogram not used    ! Restart the time manager.
! Subprogram not used    call timemgr_restart( stop_ymd=stop_ymd, stop_tod=stop_tod )
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used       ! Dynamics, physics, History
! Subprogram not used       !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call read_restart_dynamics(File, dyn_in, dyn_out, NLFileName)   
! Subprogram not used 
! Subprogram not used    call initcom ()
! Subprogram not used    call phys_grid_init
! Subprogram not used 
! Subprogram not used    call hub2atm_alloc(cam_in)
! Subprogram not used    call atm2hub_alloc(cam_out)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    ! Initialize physics grid reference pressures (needed by initialize_radbuffer)
! Subprogram not used    call ref_pres_init()
! Subprogram not used 
! Subprogram not used    call read_restart_physics(File, cam_in, cam_out, pbuf2d)
! Subprogram not used 
! Subprogram not used    if (nlres .and. .not.lbrnch) then
! Subprogram not used       call read_restart_history ( File )
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    call pio_closefile(File)
! Subprogram not used    
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used    ! Allocate communication buffers for collective communications
! Subprogram not used    ! between physics and dynamics, if necessary
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    call spmdbuf ()
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    ! Initialize ghg surface values.
! Subprogram not used 
! Subprogram not used    call chem_surfvals_init()
! Subprogram not used    call clean_iodesc_list()
! Subprogram not used 
! Subprogram not used  end subroutine cam_read_restart

!#######################################################################

! Subprogram not used    subroutine write_rest_pfile
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used !
! Subprogram not used ! Write out the restart pointer file
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used    use restart_physics, only: get_abs_restart_filepath
! Subprogram not used    use cam_history,     only: get_ptapes, get_hist_restart_filepath, &
! Subprogram not used                               hstwr, get_hfilepath, nfils, mfilt
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    integer t      ! Tape number
! Subprogram not used    integer mtapes ! Number of tapes that are active
! Subprogram not used 
! Subprogram not used    if ( nsds == -1 ) nsds = getunit()
! Subprogram not used    call opnfil(rest_pfile, nsds, 'f')
! Subprogram not used    rewind nsds
! Subprogram not used    write (nsds,'(a)') trim(pname)
! Subprogram not used    write (nsds,'(//a,a)') '# The following lists the other files needed for restarts', &
! Subprogram not used                         ' (cam only reads the first line of this file).'
! Subprogram not used    write (nsds,'(a,a)') '# The files below refer to the files needed for the master restart file:', &
! Subprogram not used                        trim(pname)
! Subprogram not used    if ( aeres )then
! Subprogram not used       write (nsds,'(a,a)') '# ', trim(get_abs_restart_filepath())
! Subprogram not used    end if
! Subprogram not used !
! Subprogram not used ! History files: Need restart history files when not a time-step to write history info
! Subprogram not used ! Need: history files if they are not full
! Subprogram not used !
! Subprogram not used    mtapes = get_ptapes( )
! Subprogram not used    do t=1,mtapes
! Subprogram not used       if ( .not. hstwr(t) ) then
! Subprogram not used          write (nsds,'(a,a)') '# ', trim(get_hist_restart_filepath( t ))
! Subprogram not used       end if
! Subprogram not used       if ( nfils(t) > 0 .and. nfils(t) < mfilt(t) ) then
! Subprogram not used          write (nsds,'(a,a)') '# ', trim(get_hfilepath( t ))
! Subprogram not used       end if
! Subprogram not used    end do
! Subprogram not used    close (nsds)
! Subprogram not used    write(iulog,*)'(WRITE_REST_PFILE): successfully wrote local restart pointer file ',trim(rest_pfile)
! Subprogram not used    write(iulog,'("---------------------------------------")')
! Subprogram not used    end subroutine write_rest_pfile

!#######################################################################

! Subprogram not used    subroutine read_rest_pfile
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used !
! Subprogram not used ! Read the master restart file from the restart pointer file
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used    character(len=nlen) :: locfn        ! Local pathname for restart pointer file
! Subprogram not used 
! Subprogram not used    nsds = getunit()
! Subprogram not used    call opnfil (rest_pfile, nsds, 'f', status="old")
! Subprogram not used    read (nsds,'(a)') pname
! Subprogram not used    
! Subprogram not used    close(nsds)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    end subroutine read_rest_pfile

!#######################################################################

!-----------------------------------------------------------------------
! BOP
!
! !ROUTINE: set_restart_filepath
!
! !DESCRIPTION: Set the filepath of the specific type of restart file.
!
!-----------------------------------------------------------------------
! !INTERFACE:
! Subprogram not used subroutine set_restart_filepath( rgpath )
! Subprogram not used !
! Subprogram not used ! !PARAMETERS:
! Subprogram not used !
! Subprogram not used   character(len=*), intent(in)  :: rgpath ! Full pathname to restart file
! Subprogram not used !
! Subprogram not used ! EOP
! Subprogram not used !
! Subprogram not used   if ( trim(rgpath) == '' )then
! Subprogram not used      call endrun ('set_restart_filepath: rgpath sent into subroutine is empty')
! Subprogram not used   end if
! Subprogram not used   if ( rgpath(1:1) /= '/' )then
! Subprogram not used      call endrun ('set_restart_filepath: rgpath sent into subroutine is not an absolute pathname')
! Subprogram not used   end if
! Subprogram not used   if ( len_trim(rgpath) > nlen )then
! Subprogram not used      call endrun ('set_restart_filepath: rgpath is too long :'//rgpath)
! Subprogram not used   end if
! Subprogram not used   pname = trim(rgpath)
! Subprogram not used end subroutine set_restart_filepath

!#######################################################################

!-----------------------------------------------------------------------
! BOP
!
! !FUNCTION: get_restcase
!
! !DESCRIPTION: Get the caseid of the case being read in
!
!-----------------------------------------------------------------------
! !INTERFACE:
! Subprogram not used character(len=nlen) function get_restcase()
! Subprogram not used !
! Subprogram not used ! EOP
! Subprogram not used !
! Subprogram not used   if ( trim(tcase) == '' )then
! Subprogram not used      call endrun ('GET_RESTCASE: caseid read in is empty, is this call after cam_read_restart?')
! Subprogram not used   end if
! Subprogram not used   get_restcase = tcase
! Subprogram not used end function get_restcase

!#######################################################################

!-----------------------------------------------------------------------
! BOP
!
! !FUNCTION: get_restartdir
!
! !DESCRIPTION: Get the directory of the restart file being read in
!
!-----------------------------------------------------------------------
! !INTERFACE:
! Subprogram not used character(len=nlen) function get_restartdir()
! Subprogram not used   use filenames,   only: get_dir
! Subprogram not used !
! Subprogram not used ! EOP
! Subprogram not used !
! Subprogram not used   ! Uses pname, so will be updated after a restart file is written out
! Subprogram not used   if ( trim(pname) == '' )then
! Subprogram not used      call endrun ('GET_RESTDIR: restart filename is empty, is this call after cam_read_restart?')
! Subprogram not used   end if
! Subprogram not used   get_restartdir = get_dir(pname)
! Subprogram not used end function get_restartdir

!=========================================================================================

end module cam_restart
