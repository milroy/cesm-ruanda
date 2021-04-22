module restFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restFileMod
!
! !DESCRIPTION:
! Reads from or writes to/ the CLM restart file.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog, use_cn
  use surfrdMod   , only : crop_prog
  use ncdio_pio   , only : file_desc_t, ncd_pio_createfile, ncd_pio_openfile, ncd_global, &
                           ncd_pio_closefile, ncd_defdim, ncd_putatt, ncd_enddef, check_dim
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restFile_read
  public :: restFile_write
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename        ! Sets restart filename
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile    ! Writes restart pointer file
  private :: restFile_closeRestart   ! Close restart file and write restart pointer file
  private :: restFile_dimset
  private :: restFile_dimcheck
  private :: restFile_enddef
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !PRIVATE TYPES: None
  private
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write
!
! !INTERFACE:
! Subprogram not used   subroutine restFile_write( file, nlend, noptr, rdate )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Read/write CLM restart file.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_time_manager , only : timemgr_restart_io, get_nstep
! Subprogram not used     use subgridRestMod   , only : SubgridRest
! Subprogram not used     use BiogeophysRestMod, only : BiogeophysRest
! Subprogram not used     use CNRestMod        , only : CNRest
! Subprogram not used     use CropRestMod      , only : CropRest
! Subprogram not used     use accumulMod       , only : accumulRest
! Subprogram not used     use histFileMod      , only : hist_restart_ncd
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*) , intent(in) :: file            ! output netcdf restart file
! Subprogram not used     logical,           intent(in) :: nlend	     ! if at the end of the simulation
! Subprogram not used     character(len=*) , intent(in) :: rdate           ! restart file time stamp for name
! Subprogram not used     logical,           intent(in), optional :: noptr ! if should NOT write to the restart pointer file
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine clm_driver2
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Mariana Vertenstein
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     type(file_desc_t) :: ncid ! netcdf id
! Subprogram not used     integer :: i       ! index
! Subprogram not used     logical :: ptrfile ! write out the restart pointer file
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if ( present(noptr) )then
! Subprogram not used        ptrfile = .not. noptr
! Subprogram not used     else
! Subprogram not used        ptrfile = .true.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     ! Open restart file
! Subprogram not used     ! --------------------------------------------
! Subprogram not used 
! Subprogram not used     call restFile_open( flag='write', file=file, ncid=ncid )
! Subprogram not used 
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     ! Define dimensions and variables
! Subprogram not used     ! --------------------------------------------
! Subprogram not used 
! Subprogram not used     call restFile_dimset ( ncid )
! Subprogram not used 
! Subprogram not used     ! Define restart file variables
! Subprogram not used 
! Subprogram not used     call timemgr_restart_io( ncid, flag='define' )
! Subprogram not used 
! Subprogram not used     call SubgridRest( ncid, flag='define' )
! Subprogram not used 
! Subprogram not used     call BiogeophysRest( ncid, flag='define' )
! Subprogram not used 
! Subprogram not used     if (use_cn) then
! Subprogram not used        call CNRest( ncid, flag='define' )
! Subprogram not used        if ( crop_prog ) call CropRest( ncid, flag='define' )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call accumulRest( ncid, flag='define' )
! Subprogram not used 
! Subprogram not used     call hist_restart_ncd ( ncid, flag='define', rdate=rdate )
! Subprogram not used 
! Subprogram not used     call restFile_enddef( ncid )
! Subprogram not used 
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     ! Write restart file variables
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     
! Subprogram not used     call timemgr_restart_io( ncid, flag='write' )
! Subprogram not used 
! Subprogram not used     call SubgridRest( ncid, flag='write' )
! Subprogram not used 
! Subprogram not used     call BiogeophysRest( ncid, flag='write' )
! Subprogram not used 
! Subprogram not used     if (use_cn) then
! Subprogram not used        call CNRest( ncid, flag='write' )
! Subprogram not used        if ( crop_prog ) call CropRest( ncid, flag='write' )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call accumulRest( ncid, flag='write' )
! Subprogram not used     
! Subprogram not used     call hist_restart_ncd (ncid, flag='write' )
! Subprogram not used 
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     ! Close restart file and write restart pointer file
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     
! Subprogram not used     call restFile_close( ncid )
! Subprogram not used     call restFile_closeRestart( file, nlend )
! Subprogram not used     
! Subprogram not used     ! Write restart pointer file
! Subprogram not used     
! Subprogram not used     if ( ptrfile ) call restFile_write_pfile( file )
! Subprogram not used     
! Subprogram not used     ! Write out diagnostic info
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
! Subprogram not used        write(iulog,'(72a1)') ("-",i=1,60)
! Subprogram not used     end if
! Subprogram not used     
! Subprogram not used   end subroutine restFile_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read
!
! !INTERFACE:
  subroutine restFile_read( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use BiogeophysRestMod, only : BiogeophysRest
    use CNRestMod        , only : CNRest
    use CropRestMod      , only : CropRest
    use accumulMod       , only : accumulRest
    use histFileMod      , only : hist_restart_ncd
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! output netcdf restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i              ! index
!-----------------------------------------------------------------------

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )

    call BiogeophysRest( ncid, flag='read' )

    if (use_cn) then
       call CNRest( ncid, flag='read' )
       if ( crop_prog ) call CropRest( ncid, flag='read' )
    end if

    call accumulRest( ncid, flag='read' )
    
    call hist_restart_ncd (ncid, flag='read')

    ! Close file 

    call restFile_close( ncid )

    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine restFile_read

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_getfile
!
! !INTERFACE:
  subroutine restFile_getfile( file, path )
!
! !DESCRIPTION:
! Determine and obtain netcdf restart file
!
! !USES:
    use clm_varctl, only : caseid, finidat, nrevsn, nsrest, brnch_retain_casename, &
                           nsrContinue, nsrBranch, nsrStartup
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
!-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 
    
    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if
       
    ! Branch run: 
    ! Restart file pathname is obtained from namelist "nrevsn"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)
    
    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn)
       if (nrevsn(length-2:length) == '.nc') then
          path = trim(nrevsn) 
       else
          path = trim(nrevsn) // '.nc'
       end if
       call getfil( path, file, 0 )
       
       ! tcraig, adding xx. and .clm2 makes this more robust
       ctest = 'xx.'//trim(caseid)//'.clm2'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          write(iulog,*) 'Must change case name on branch run if ',&
               'brnch_retain_casename namelist is not set'
          write(iulog,*) 'previous case filename= ',trim(file),&
               ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
               ' ftest = ',trim(ftest)
          call endrun()
       end if
    end if

    ! Initial run: 
    ! Restart file pathname is obtained from namelist "finidat"
    
    if (nsrest==nsrStartup) then
       call getfil( finidat, file, 0 )
    end if
    
  end subroutine restFile_getfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read_pfile
!
! !INTERFACE:
! Subprogram not used   subroutine restFile_read_pfile( pnamer )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Setup restart file and perform necessary consistency checks
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use fileutils , only : opnfil, getavu, relavu
! Subprogram not used     use clm_varctl, only : rpntfil, rpntdir, inst_suffix
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(out) :: pnamer ! full path of restart file
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine restart in this module
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer :: i                  ! indices
! Subprogram not used     integer :: nio                ! restart unit
! Subprogram not used     integer :: status             ! substring check status
! Subprogram not used     character(len=256) :: locfn   ! Restart pointer file name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Obtain the restart file from the restart pointer file. 
! Subprogram not used     ! For restart runs, the restart pointer file contains the full pathname 
! Subprogram not used     ! of the restart file. For branch runs, the namelist variable 
! Subprogram not used     ! [nrevsn] contains the full pathname of the restart file. 
! Subprogram not used     ! New history files are always created for branch runs.
! Subprogram not used        
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Reading restart pointer file....'
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     nio = getavu()
! Subprogram not used     locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
! Subprogram not used     call opnfil (locfn, nio, 'f')
! Subprogram not used     read (nio,'(a256)') pnamer
! Subprogram not used     call relavu (nio)
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Reading restart data.....'
! Subprogram not used        write(iulog,'(72a1)') ("-",i=1,60)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine restFile_read_pfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_closeRestart
!
! !INTERFACE:
! Subprogram not used   subroutine restFile_closeRestart( file, nlend )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Close restart file and write restart pointer file if
! Subprogram not used ! in write mode, otherwise just close restart file if in read mode
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_time_manager, only : is_last_step
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*) , intent(in) :: file  ! local output filename
! Subprogram not used     logical,           intent(in) :: nlend
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine restart in this module
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Mariana Vertenstein
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer :: i                   !index
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
! Subprogram not used       write(iulog,'(72a1)') ("-",i=1,60)
! Subprogram not used       write(iulog,*)
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used  end subroutine restFile_closeRestart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write_pfile
!
! !INTERFACE:
! Subprogram not used   subroutine restFile_write_pfile( fnamer )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Open restart pointer file. Write names of current netcdf restart file.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varctl, only : rpntdir, rpntfil, inst_suffix
! Subprogram not used     use fileutils , only : relavu
! Subprogram not used     use fileutils , only : getavu, opnfil
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in) :: fnamer
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine restart in this module
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Mariana Vertenstein
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer :: m                    ! index
! Subprogram not used     integer :: nio                  ! restart pointer file
! Subprogram not used     character(len=256) :: filename  ! local file name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        nio = getavu()
! Subprogram not used        filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
! Subprogram not used        call opnfil( filename, nio, 'f' )
! Subprogram not used        
! Subprogram not used        write(nio,'(a)') fnamer
! Subprogram not used        call relavu( nio )
! Subprogram not used        write(iulog,*)'Successfully wrote local restart pointer file'
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine restFile_write_pfile

!-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use clm_time_manager, only : get_nstep
    
    implicit none
    character(len=*),  intent(in) :: flag ! flag to specify read or write
    character(len=*),  intent(in) :: file ! filename
    type(file_desc_t), intent(out):: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance
       
       if (masterproc) then	
          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',get_nstep()
          write(iulog,*)
       end if
       call ncd_pio_createfile(ncid, trim(file))
       
    else if (flag == 'read') then
       
       ! Open netcdf restart file
       
       if (masterproc) then
          write(iulog,*) 'Reading restart dataset'
       end if
       call ncd_pio_openfile (ncid, trim(file), 0)
       
    end if
  
  end subroutine restFile_open

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_filename
!
! !INTERFACE:
! Subprogram not used   character(len=256) function restFile_filename( rdate )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varctl, only : caseid, inst_suffix
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in) :: rdate   ! input date for restart file name 
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine restart in this module
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Mariana Vertenstein
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     restFile_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
! Subprogram not used                         ".r."//trim(rdate)//".nc"
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
! Subprogram not used     end if
! Subprogram not used  
! Subprogram not used   end function restFile_filename

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimset
!
! !INTERFACE:
! Subprogram not used   subroutine restFile_dimset( ncid )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Read/Write initial data from/to netCDF instantaneous initial data file
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use shr_kind_mod, only : r8 => shr_kind_r8
! Subprogram not used     use clm_time_manager, only : get_nstep, get_curr_date
! Subprogram not used     use spmdMod     , only : mpicom, MPI_LOGICAL
! Subprogram not used     use clm_varctl  , only : caseid, ctitle, version, username, hostname, fsurdat, &
! Subprogram not used                              conventions, source
! Subprogram not used     use clm_varpar  , only : numrad, nlevlak, nlevsno, nlevgrnd
! Subprogram not used     use decompMod   , only : get_proc_bounds, get_proc_global
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer :: yr                  ! current year (0 -> ...)
! Subprogram not used     integer :: mon                 ! current month (1 -> 12)
! Subprogram not used     integer :: day                 ! current day (1 -> 31)
! Subprogram not used     integer :: mcsec               ! seconds of current date
! Subprogram not used     integer :: mcdate              ! current date
! Subprogram not used     integer :: dimid               ! netCDF dimension id
! Subprogram not used     integer :: numg                ! total number of gridcells across all processors
! Subprogram not used     integer :: numl                ! total number of landunits across all processors
! Subprogram not used     integer :: numc                ! total number of columns across all processors
! Subprogram not used     integer :: nump                ! total number of pfts across all processors
! Subprogram not used     integer :: ier                 ! error status
! Subprogram not used     integer :: strlen_dimid        ! string dimension id
! Subprogram not used     character(len=  8) :: curdate  ! current date
! Subprogram not used     character(len=  8) :: curtime  ! current time
! Subprogram not used     character(len=256) :: str
! Subprogram not used     character(len= 32) :: subname='restFile_dimset' ! subroutine name
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call get_proc_global(numg, numl, numc, nump)
! Subprogram not used 
! Subprogram not used     ! Define dimensions
! Subprogram not used     
! Subprogram not used     call ncd_defdim(ncid, 'gridcell', numg           , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'landunit', numl           , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'column'  , numc           , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'pft'     , nump           , dimid)
! Subprogram not used     
! Subprogram not used     call ncd_defdim(ncid, 'levgrnd' , nlevgrnd       , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'levlak'  , nlevlak        , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'levsno'  , nlevsno        , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'levsno1'  , nlevsno+1     , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'levtot'  , nlevsno+nlevgrnd, dimid)
! Subprogram not used     call ncd_defdim(ncid, 'numrad'  , numrad         , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'string_length', 64        , dimid)
! Subprogram not used        
! Subprogram not used     ! Define global attributes
! Subprogram not used     
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
! Subprogram not used     call getdatetime(curdate, curtime)
! Subprogram not used     str = 'created on ' // curdate // ' ' // curtime
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
! Subprogram not used     str = '$Id: restFileMod.F90 60694 2014-05-27 22:16:25Z erik $'
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(fsurdat))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'title', &
! Subprogram not used           'CLM Restart information, required to continue a simulation' )
! Subprogram not used 
! Subprogram not used     
! Subprogram not used   end subroutine restFile_dimset
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimcheck
!
! !INTERFACE:
  subroutine restFile_dimcheck( ncid )
!
! !DESCRIPTION:
! Check dimensions of restart file
!
! !USES:
    use decompMod,  only : get_proc_bounds, get_proc_global
    use clm_varpar, only : nlevsno, nlevlak, nlevgrnd
    use clm_varctl, only : single_column, nsrest, nsrStartup
    implicit none
!
! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: numg     ! total number of gridcells across all processors
    integer :: numl     ! total number of landunits across all processors
    integer :: numc     ! total number of columns across all processors
    integer :: nump     ! total number of pfts across all processors
    character(len=32) :: subname='restFile_dimcheck' ! subroutine name
!-----------------------------------------------------------------------

    ! Get relevant sizes

    if ( .not. single_column .or. nsrest /= nsrStartup )then
       call get_proc_global(numg, numl, numc, nump)
       call check_dim(ncid, 'gridcell', numg)
       call check_dim(ncid, 'landunit', numl)
       call check_dim(ncid, 'column'  , numc)
       call check_dim(ncid, 'pft'     , nump)
    end if
    call check_dim(ncid, 'levsno'  , nlevsno)
    call check_dim(ncid, 'levgrnd' , nlevgrnd)
    call check_dim(ncid, 'levlak'  , nlevlak) 

  end subroutine restFile_dimcheck

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_enddef
!
! !INTERFACE:
! Subprogram not used   subroutine restFile_enddef( ncid )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Read a CLM restart file.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Mariana Vertenstein
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call ncd_enddef(ncid)
! Subprogram not used 
! Subprogram not used   end subroutine restFile_enddef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_close
!
! !INTERFACE:
  subroutine restFile_close( ncid )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=32) :: subname='restFile_close' ! subroutine name
!-----------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

  end subroutine restFile_close

end module restFileMod



