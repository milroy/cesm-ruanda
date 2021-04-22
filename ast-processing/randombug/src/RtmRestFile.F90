module RtmRestFile

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restFileMod
!
! !DESCRIPTION:
! Reads from or writes to/ the RTM restart file.
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort
  use RtmSpmd       , only : masterproc 
  use RtmVar        , only : rtmlon, rtmlat, iulog, inst_suffix, rpntfil, &
                             caseid, nsrest, brnch_retain_casename, &
                             finidat_rtm, nrevsn_rtm, &
                             nsrContinue, nsrBranch, nsrStartup, &
                             ctitle, version, username, hostname, conventions, source
  use RtmHistFile   , only : RtmHistRestart
  use RtmFileUtils  , only : relavu, getavu, opnfil, getfil
  use RtmTimeManager, only : timemgr_restart, get_nstep, get_curr_date, is_last_step
  use RunoffMod     , only : runoff
  use RtmIO       
  use RtmDateTime
  use rtm_cpl_indices , only : nt_rtm, rtm_tracers 
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: RtmRestFileName  
  public :: RtmRestFileRead
  public :: RtmRestFileWrite
  public :: RtmRestGetfile
  public :: RtmRestTimeManager
  public :: RtmRestart
  public :: RtmRestFinalize
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile    ! Writes restart pointer file
  private :: restFile_dimset
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
! !PRIVATE TYPES: None
  private

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!=======================================================================

! Subprogram not used   subroutine RtmRestFileWrite( file, rdate )
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Read/write RTM restart file.
! Subprogram not used 
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*) , intent(in) :: file            ! output netcdf restart file
! Subprogram not used     character(len=*) , intent(in) :: rdate           ! restart file time stamp for name
! Subprogram not used 
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     type(file_desc_t) :: ncid ! netcdf id
! Subprogram not used     integer :: i       ! index
! Subprogram not used     logical :: ptrfile ! write out the restart pointer file
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Define dimensions and variables
! Subprogram not used 
! Subprogram not used     if (masterproc) then	
! Subprogram not used        write(iulog,*)
! Subprogram not used        write(iulog,*)'restFile_open: writing RTM restart dataset '
! Subprogram not used        write(iulog,*)
! Subprogram not used     end if
! Subprogram not used     call ncd_pio_createfile(ncid, trim(file))
! Subprogram not used     call restFile_dimset( ncid )
! Subprogram not used     call RtmRestart( ncid, flag='define' )
! Subprogram not used     call RtmHistRestart ( ncid, flag='define', rdate=rdate )
! Subprogram not used     call timemgr_restart( ncid, flag='define' )
! Subprogram not used     call ncd_enddef(ncid)
! Subprogram not used 
! Subprogram not used     ! Write restart file variables
! Subprogram not used     call RtmRestart( ncid, flag='write' )
! Subprogram not used     call RtmHistRestart ( ncid, flag='write' )
! Subprogram not used     call timemgr_restart( ncid, flag='write' )
! Subprogram not used     call ncd_pio_closefile(ncid)
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Successfully wrote local restart file ',trim(file)
! Subprogram not used        write(iulog,'(72a1)') ("-",i=1,60)
! Subprogram not used        write(iulog,*)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Write restart pointer file
! Subprogram not used     call restFile_write_pfile( file )
! Subprogram not used     
! Subprogram not used     ! Write out diagnostic info
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
! Subprogram not used        write(iulog,'(72a1)') ("-",i=1,60)
! Subprogram not used     end if
! Subprogram not used     
! Subprogram not used   end subroutine RtmRestFileWrite

!-----------------------------------------------------------------------

  subroutine RtmRestFileRead( file )

    ! !DESCRIPTION:
    ! Read a RTM restart file.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! output netcdf restart file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i              ! index
    !-------------------------------------

    ! Read file
    if (masterproc) write(iulog,*) 'Reading restart dataset'
    call ncd_pio_openfile (ncid, trim(file), 0)
    call RtmRestart( ncid, flag='read' )
    call RtmHistRestart(ncid, flag='read')
    call ncd_pio_closefile(ncid)

    ! Write out diagnostic info
    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine RtmRestFileRead

!-----------------------------------------------------------------------

! Subprogram not used   subroutine RtmRestTimeManager( file )
! Subprogram not used 
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Read a RTM restart file.
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in) :: file  ! output netcdf restart file
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     type(file_desc_t) :: ncid ! netcdf id
! Subprogram not used     integer :: i              ! index
! Subprogram not used     !-------------------------------------
! Subprogram not used 
! Subprogram not used     ! Read file
! Subprogram not used     if (masterproc) write(iulog,*) 'Reading restart Timemanger'
! Subprogram not used     call ncd_pio_openfile (ncid, trim(file), 0)
! Subprogram not used     call timemgr_restart(ncid, flag='read')
! Subprogram not used     call ncd_pio_closefile(ncid)
! Subprogram not used 
! Subprogram not used     ! Write out diagnostic info
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,'(72a1)') ("-",i=1,60)
! Subprogram not used        write(iulog,*) 'Successfully read restart data for restart run'
! Subprogram not used        write(iulog,*)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine RtmRestTimeManager

!-----------------------------------------------------------------------

  subroutine RtmRestGetfile( file, path )

    !---------------------------------------------------
    ! DESCRIPTION:
    ! Determine and obtain netcdf restart file

    ! ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file

    ! LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !---------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 
    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if
       
    ! Branch run: 
    ! Restart file pathname is obtained from namelist "nrevsn_rtm"
    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn_rtm)
       if (nrevsn_rtm(length-2:length) == '.nc') then
          path = trim(nrevsn_rtm) 
       else
          path = trim(nrevsn_rtm) // '.nc'
       end if
       call getfil( path, file, 0 )
       
       ! Check case name consistency (case name must be different 
       ! for branch run, unless brnch_retain_casename is set)
       ctest = 'xx.'//trim(caseid)//'.rtm'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          write(iulog,*) 'Must change case name on branch run if ',&
               'brnch_retain_casename namelist is not set'
          write(iulog,*) 'previous case filename= ',trim(file),&
               ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
               ' ftest = ',trim(ftest)
          call shr_sys_abort()
       end if
    end if

    ! Initial run 
    if (nsrest==nsrStartup) then
       call getfil( finidat_rtm, file, 0 )
    end if
    
  end subroutine RtmRestGetfile

!-----------------------------------------------------------------------

! Subprogram not used   subroutine restFile_read_pfile( pnamer )
! Subprogram not used 
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Setup restart file and perform necessary consistency checks
! Subprogram not used 
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(out) :: pnamer ! full path of restart file
! Subprogram not used 
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: i                  ! indices
! Subprogram not used     integer :: nio                ! restart unit
! Subprogram not used     integer :: status             ! substring check status
! Subprogram not used     character(len=256) :: locfn   ! Restart pointer file name
! Subprogram not used     !--------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Obtain the restart file from the restart pointer file. 
! Subprogram not used     ! For restart runs, the restart pointer file contains the full pathname 
! Subprogram not used     ! of the restart file. For branch runs, the namelist variable 
! Subprogram not used     ! [nrevsn_rtm] contains the full pathname of the restart file. 
! Subprogram not used     ! New history files are always created for branch runs.
! Subprogram not used        
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Reading restart pointer file....'
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     nio = getavu()
! Subprogram not used     locfn = './'// trim(rpntfil)//trim(inst_suffix)
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

! Subprogram not used   subroutine restFile_write_pfile( fnamer )
! Subprogram not used 
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Open restart pointer file. Write names of current netcdf restart file.
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in) :: fnamer
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: m                    ! index
! Subprogram not used     integer :: nio                  ! restart pointer file
! Subprogram not used     character(len=256) :: filename  ! local file name
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        nio = getavu()
! Subprogram not used        filename= './'// trim(rpntfil)//trim(inst_suffix)
! Subprogram not used        call opnfil( filename, nio, 'f' )
! Subprogram not used        
! Subprogram not used        write(nio,'(a)') fnamer
! Subprogram not used        call relavu( nio )
! Subprogram not used        write(iulog,*)'Successfully wrote local restart pointer file'
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine restFile_write_pfile


!-----------------------------------------------------------------------

! Subprogram not used   character(len=256) function RtmRestFileName( rdate )
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in) :: rdate   ! input date for restart file name 
! Subprogram not used 
! Subprogram not used     RtmRestFileName = "./"//trim(caseid)//".rtm"//trim(inst_suffix)//".r."//trim(rdate)//".nc"
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*)'writing restart file ',trim(RtmRestFileName),' for model date = ',rdate
! Subprogram not used     end if
! Subprogram not used  
! Subprogram not used   end function RtmRestFileName

!------------------------------------------------------------------------

! Subprogram not used   subroutine restFile_dimset( ncid )
! Subprogram not used 
! Subprogram not used     !----------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Read/Write initial data from/to netCDF instantaneous initial data file
! Subprogram not used 
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid
! Subprogram not used 
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: dimid               ! netCDF dimension id
! Subprogram not used     integer :: ier                 ! error status
! Subprogram not used     character(len=  8) :: curdate  ! current date
! Subprogram not used     character(len=  8) :: curtime  ! current time
! Subprogram not used     character(len=256) :: str
! Subprogram not used     character(len= 32) :: subname='restFile_dimset' ! subroutine name
! Subprogram not used     !----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Define dimensions
! Subprogram not used     
! Subprogram not used     call ncd_defdim(ncid, 'rtmlon'  , rtmlon         , dimid)
! Subprogram not used     call ncd_defdim(ncid, 'rtmlat'  , rtmlat         , dimid)
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
! Subprogram not used     str = '$Id: restFileMod.F90 35676 2012-03-22 21:48:04Z tcraig $'
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
! Subprogram not used     call ncd_putatt(ncid, NCD_GLOBAL, 'title', &
! Subprogram not used           'RTM Restart information, required to continue a simulation' )
! Subprogram not used 
! Subprogram not used   end subroutine restFile_dimset
  
!-----------------------------------------------------------------------

  subroutine RtmRestart(ncid, flag)

    !-----------------------------------------------------------------------
    ! DESCRIPTION:
    ! Read/write RTM restart data.
    !
    ! ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout)  :: ncid ! netcdf id
    character(len=*) , intent(in) :: flag   ! 'read' or 'write'
    ! LOCAL VARIABLES:
    logical :: readvar          ! determine if variable is on initial file
    integer :: nt,nv,n          ! indices
    real(r8) , pointer :: dfld(:) ! temporary array
    character(len=32)  :: vname,uname
    character(len=128) :: lname
    !-----------------------------------------------------------------------

    do nv = 1,4
    do nt = 1,nt_rtm

       if (nv == 1) then
          vname = 'RTM_VOLR_'//trim(rtm_tracers(nt))
          lname = 'water volume in cell (volr)'
          uname = 'm3'
          dfld  => runoff%volr(:,nt)
       elseif (nv == 2) then
          vname = 'RTM_FLUXOUT_'//trim(rtm_tracers(nt))
          lname = 'water fluxout in cell (fluxout)'
          uname = 'm3/s'
          dfld  => runoff%fluxout(:,nt)
       elseif (nv == 3) then
          vname = 'RTM_RUNOFF_'//trim(rtm_tracers(nt))
          lname = 'runoff (runoff)'
          uname = 'm3/s'
          dfld  => runoff%runoff(:,nt)
       elseif (nv == 4) then
          vname = 'RTM_DVOLRDT_'//trim(rtm_tracers(nt))
          lname = 'water volume change in cell (dvolrdt)'
          uname = 'mm/s'
          dfld  => runoff%dvolrdt(:,nt)
       else
          write(iulog,*) 'Rtm ERROR: illegal nv value a ',nv
          call shr_sys_abort()
       endif

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(vname), &
               xtype=ncd_double,  dim1name='rtmlon', dim2name='rtmlat', &
               long_name=trim(lname), units=trim(uname))
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(vname), data=dfld, dim1name='allrof', &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (nsrest == nsrContinue) then
                call shr_sys_abort()
             else
                dfld = 0._r8
             end if
          end if
       end if

    enddo
    enddo

    if (flag == 'read') then
       do n = runoff%begr,runoff%endr
          do nt = 1,nt_rtm
             if (abs(runoff%volr(n,nt))    > 1.e30) runoff%volr(n,nt) = 0.
             if (abs(runoff%runoff(n,nt))  > 1.e30) runoff%runoff(n,nt) = 0.
             if (abs(runoff%dvolrdt(n,nt)) > 1.e30) runoff%dvolrdt(n,nt) = 0.
             if (abs(runoff%fluxout(n,nt)) > 1.e30) runoff%fluxout(n,nt) = 0.
          end do
          if (runoff%mask(n) == 1) then
             do nt = 1,nt_rtm
                runoff%runofflnd(n,nt) = runoff%runoff(n,nt)
                runoff%dvolrdtlnd(n,nt)= runoff%dvolrdt(n,nt)
                runoff%volrlnd(n,nt)   = runoff%volr(n,nt)
             end do
          elseif (runoff%mask(n) == 2) then
             do nt = 1,nt_rtm
                runoff%runoffocn(n,nt) = runoff%runoff(n,nt)
                runoff%dvolrdtocn(n,nt)= runoff%dvolrdt(n,nt)
             enddo
          endif
       enddo
    endif

  end subroutine RtmRestart

  subroutine RtmRestFinalize( )

    ! !DESCRIPTION:
    ! clean up memory after a RTM restart handling
    !
    use RtmIO , only : ncd_finalize
    !
    ! !ARGUMENTS:
    implicit none

       call ncd_finalize()

  end subroutine RtmRestFinalize

end module RtmRestFile



