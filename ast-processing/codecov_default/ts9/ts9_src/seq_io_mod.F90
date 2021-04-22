!===============================================================================
! SVN $Id: seq_io_mod.F90 61512 2014-06-26 21:59:35Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq5_0_14/driver/seq_io_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_io_mod -- reads and writes driver files
!
! !DESCRIPTION:
!    Writes attribute vectors to netcdf
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig first version
!    2007-Dec-06 - T. Craig update and improve
!    2008-Feb-16 - J. Edwards convert to PIO
!    2010-Nov    - J. Edwards move PIO init and namelists from components to driver
! Current Problems
!  - the original use of seq_io will now ONLY work with the cpl because
!    of hardwiring cpl_io_type and cpl_io_iosystem.  want the original
!    io capabilities to be usable by any component
!  - the init1 method depends on seq_comm for name consistency but seq_comm_init 
!    wants to be called after init1 so the global_comm can be modified for
!    async IO.  this needs to be reconciled.
!  - this routine stores information for all components but most methods are
!    hardwired to work only for the coupler.  should all the components info
!    be stored here or should this be more a general set of methods that are
!    reusable as it's original intent.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_io_mod

  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in => shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl, cs => shr_kind_cs
  use shr_sys_mod       ! system calls
  use seq_comm_mct
  use seq_flds_mod, only : seq_flds_lookup
  use mct_mod           ! mct wrappers
  use pio
  use component_type_mod

  implicit none
  private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

  public seq_io_wopen
  public seq_io_close
  public seq_io_redef
  public seq_io_enddef
  public seq_io_date2yyyymmdd
  public seq_io_sec2hms
  public seq_io_read
  public seq_io_write
  public seq_io_cpl_init
! !PUBLIC DATA MEMBERS


  ! none

!EOP

  interface seq_io_read
     module procedure seq_io_read_av
     module procedure seq_io_read_avs
     module procedure seq_io_read_avscomp
     module procedure seq_io_read_int
     module procedure seq_io_read_int1d
     module procedure seq_io_read_r8
     module procedure seq_io_read_r81d
     module procedure seq_io_read_char
  end interface
  interface seq_io_write
     module procedure seq_io_write_av
     module procedure seq_io_write_avs
     module procedure seq_io_write_avscomp
     module procedure seq_io_write_int
     module procedure seq_io_write_int1d
     module procedure seq_io_write_r8
     module procedure seq_io_write_r81d
     module procedure seq_io_write_char
     module procedure seq_io_write_time
  end interface

!-------------------------------------------------------------------------------
! Local data
!-------------------------------------------------------------------------------

   character(*),parameter :: prefix = "seq_io_"
   character(CL)          :: wfilename = ''
   real(r8)    ,parameter :: fillvalue = SHR_CONST_SPVAL

   character(*),parameter :: modName = "(seq_io_mod) "
   integer(in) ,parameter :: debug = 1 ! internal debug level


   type(file_desc_t), save :: cpl_io_file
   integer(IN)             :: cpl_pio_iotype
   type(iosystem_desc_t), pointer :: cpl_io_subsystem
   character(*),parameter  :: version ='cpl7v10'
   character(*),parameter  :: version0='cpl7v00'

   character(CL) :: charvar   ! buffer for string read/write
   integer(IN) :: io_comm

!=================================================================================
contains
!=================================================================================

  subroutine seq_io_cpl_init()
    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype

    character(len=seq_comm_namelen) :: cpl_name

    ! For consistency with how io_compname is set, we get the coupler name here
    cpl_name = seq_comm_name(CPLID)
    cpl_io_subsystem=>shr_pio_getiosys(cpl_name)
    cpl_pio_iotype = shr_pio_getiotype(cpl_name)

  end subroutine seq_io_cpl_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_io_wopen - open netcdf file
!
! !DESCRIPTION:
!    open netcdf file
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig - initial version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_io_wopen(filename,clobber,cdf64)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(*),intent(in) :: filename
! Subprogram not used     logical,optional,intent(in):: clobber
! Subprogram not used     logical,optional,intent(in):: cdf64
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lclobber
! Subprogram not used     logical :: lcdf64
! Subprogram not used     integer :: iam,mpicom
! Subprogram not used     integer :: rcode
! Subprogram not used     integer :: nmode
! Subprogram not used     character(CL)  :: lversion
! Subprogram not used     character(*),parameter :: subName = '(seq_io_wopen) '
! Subprogram not used     
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lversion=trim(version0)
! Subprogram not used 
! Subprogram not used     lclobber = .false.
! Subprogram not used     if (present(clobber)) lclobber=clobber
! Subprogram not used 
! Subprogram not used     lcdf64 = .false.
! Subprogram not used     if (present(cdf64)) lcdf64=cdf64
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used 
! Subprogram not used     if (.not. pio_file_is_open(cpl_io_file)) then
! Subprogram not used        ! filename not open
! Subprogram not used        if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used        call shr_mpi_bcast(exists,mpicom,'seq_io_wopen exists')
! Subprogram not used        if (exists) then
! Subprogram not used           if (lclobber) then
! Subprogram not used              nmode = pio_clobber
! Subprogram not used              if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
! Subprogram not used              rcode = pio_createfile(cpl_io_subsystem, cpl_io_file, cpl_pio_iotype, trim(filename), nmode)
! Subprogram not used              if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
! Subprogram not used              rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
! Subprogram not used           else
! Subprogram not used 
! Subprogram not used              rcode = pio_openfile(cpl_io_subsystem, cpl_io_file, cpl_pio_iotype, trim(filename), pio_write)
! Subprogram not used              if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used              call pio_seterrorhandling(cpl_io_file,PIO_BCAST_ERROR)
! Subprogram not used              rcode = pio_get_att(cpl_io_file,pio_global,"file_version",lversion)
! Subprogram not used              call pio_seterrorhandling(cpl_io_file,PIO_INTERNAL_ERROR)
! Subprogram not used              if (trim(lversion) /= trim(version)) then
! Subprogram not used                 rcode = pio_redef(cpl_io_file)
! Subprogram not used                 rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
! Subprogram not used                 rcode = pio_enddef(cpl_io_file)
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used        else
! Subprogram not used           nmode = pio_noclobber
! Subprogram not used           if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
! Subprogram not used           rcode = pio_createfile(cpl_io_subsystem, cpl_io_file, cpl_pio_iotype, trim(filename), nmode)
! Subprogram not used           if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
! Subprogram not used           rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
! Subprogram not used        endif
! Subprogram not used     elseif (trim(wfilename) /= trim(filename)) then
! Subprogram not used        ! filename is open, better match open filename
! Subprogram not used        if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     else
! Subprogram not used        ! filename is already open, just return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used end subroutine seq_io_wopen

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_io_close - close netcdf file
!
! !DESCRIPTION:
!    close netcdf file
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig - initial version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_io_close(filename)
! Subprogram not used 
! Subprogram not used     use pio, only : pio_closefile
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     character(*),intent(in) :: filename
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer :: iam
! Subprogram not used     integer :: rcode
! Subprogram not used     character(*),parameter :: subName = '(seq_io_close) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (.not. pio_file_is_open(cpl_io_file)) then
! Subprogram not used        ! filename not open, just return
! Subprogram not used     elseif (trim(wfilename) /= trim(filename)) then
! Subprogram not used        ! filename matches, close it
! Subprogram not used        call pio_closefile(cpl_io_file)
! Subprogram not used     else
! Subprogram not used        ! different filename is open, abort
! Subprogram not used        if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     wfilename = ''
! Subprogram not used 
! Subprogram not used end subroutine seq_io_close

!===============================================================================

! Subprogram not used subroutine seq_io_redef(filename)
! Subprogram not used     character(len=*), intent(in) :: filename
! Subprogram not used     integer :: rcode
! Subprogram not used 
! Subprogram not used     rcode = pio_redef(cpl_io_file)
! Subprogram not used end subroutine seq_io_redef

!===============================================================================

! Subprogram not used subroutine seq_io_enddef(filename)
! Subprogram not used     character(len=*), intent(in) :: filename
! Subprogram not used     integer :: rcode
! Subprogram not used 
! Subprogram not used     rcode = pio_enddef(cpl_io_file)
! Subprogram not used end subroutine seq_io_enddef

!===============================================================================

! Subprogram not used character(len=10) function seq_io_date2yyyymmdd (date)
! Subprogram not used 
! Subprogram not used ! Input arguments
! Subprogram not used 
! Subprogram not used    integer, intent(in) :: date
! Subprogram not used 
! Subprogram not used ! Local workspace
! Subprogram not used 
! Subprogram not used    integer :: year    ! year of yyyy-mm-dd
! Subprogram not used    integer :: month   ! month of yyyy-mm-dd
! Subprogram not used    integer :: day     ! day of yyyy-mm-dd
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (date < 0) then
! Subprogram not used       call shr_sys_abort ('seq_io_date2yyyymmdd: negative date not allowed')
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    year  = date / 10000
! Subprogram not used    month = (date - year*10000) / 100
! Subprogram not used    day   = date - year*10000 - month*100
! Subprogram not used 
! Subprogram not used    write(seq_io_date2yyyymmdd,80) year, month, day
! Subprogram not used 80 format(i4.4,'-',i2.2,'-',i2.2)
! Subprogram not used 
! Subprogram not used end function seq_io_date2yyyymmdd

!===============================================================================

! Subprogram not used character(len=8) function seq_io_sec2hms (seconds)
! Subprogram not used 
! Subprogram not used ! Input arguments
! Subprogram not used 
! Subprogram not used    integer, intent(in) :: seconds
! Subprogram not used 
! Subprogram not used ! Local workspace
! Subprogram not used 
! Subprogram not used    integer :: hours     ! hours of hh:mm:ss
! Subprogram not used    integer :: minutes   ! minutes of hh:mm:ss
! Subprogram not used    integer :: secs      ! seconds of hh:mm:ss
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if (seconds < 0 .or. seconds > 86400) then
! Subprogram not used       write(logunit,*)'seq_io_sec2hms: bad input seconds:', seconds
! Subprogram not used       call shr_sys_abort()
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    hours   = seconds / 3600
! Subprogram not used    minutes = (seconds - hours*3600) / 60
! Subprogram not used    secs    = (seconds - hours*3600 - minutes*60)
! Subprogram not used 
! Subprogram not used    if (minutes < 0 .or. minutes > 60) then
! Subprogram not used       write(logunit,*)'seq_io_sec2hms: bad minutes = ',minutes
! Subprogram not used       call shr_sys_abort()
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    if (secs < 0 .or. secs > 60) then
! Subprogram not used       write(logunit,*)'seq_io_sec2hms: bad secs = ',secs
! Subprogram not used       call shr_sys_abort()
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    write(seq_io_sec2hms,80) hours, minutes, secs
! Subprogram not used 80 format(i2.2,':',i2.2,':',i2.2)
! Subprogram not used 
! Subprogram not used end function seq_io_sec2hms

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_av - write AV to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_av(filename,gsmap,AV,dname,whead,wdata,nx,ny,nt,fillval,pre,tavg,&
! Subprogram not used 	                     use_float)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     type(mct_gsMap), intent(in) :: gsmap
! Subprogram not used     type(mct_aVect) ,intent(in) :: AV       ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used     integer(in),optional,intent(in) :: nx   ! 2d grid size if available
! Subprogram not used     integer(in),optional,intent(in) :: ny   ! 2d grid size if available
! Subprogram not used     integer(in),optional,intent(in) :: nt   ! time sample
! Subprogram not used     real(r8),optional,intent(in) :: fillval ! fill value
! Subprogram not used     character(len=*),optional,intent(in) :: pre      ! prefix to variable name
! Subprogram not used     logical,optional,intent(in) :: tavg     ! is this a tavg
! Subprogram not used     logical,optional,intent(in) :: use_float ! write output as float rather than double
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: mpicom
! Subprogram not used     integer(in) :: iam
! Subprogram not used     integer(in) :: nf,ns,ng
! Subprogram not used     integer(in) :: i,j,k,n
! Subprogram not used     integer(in),target  :: dimid2(2)
! Subprogram not used     integer(in),target  :: dimid3(3)
! Subprogram not used     integer(in),pointer :: dimid(:)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used     integer(kind=PIO_OffSet) :: frame
! Subprogram not used     type(mct_string) :: mstring     ! mct char type
! Subprogram not used     character(CL)    :: itemc       ! string converted to char
! Subprogram not used     character(CL)    :: name1       ! var name
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     character(CL)    :: lpre        ! local prefix
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     logical :: luse_float
! Subprogram not used     integer(in) :: lnx,lny
! Subprogram not used     real(r8) :: lfillvalue
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_av) '
! Subprogram not used     integer :: lbnum
! Subprogram not used     integer, pointer :: Dof(:)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lfillvalue = fillvalue
! Subprogram not used     if (present(fillval)) then
! Subprogram not used        lfillvalue = fillval
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lpre = trim(dname)
! Subprogram not used     if (present(pre)) then
! Subprogram not used        lpre = trim(pre)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     luse_float = .false.
! Subprogram not used     if (present(use_float)) luse_float = use_float
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)	
! Subprogram not used 
! Subprogram not used     ng = mct_gsmap_gsize(gsmap)
! Subprogram not used     lnx = ng
! Subprogram not used     lny = 1
! Subprogram not used 	
! Subprogram not used     nf = mct_aVect_nRattr(AV)
! Subprogram not used     if (nf < 1) then
! Subprogram not used        write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (present(nx)) then
! Subprogram not used        if (nx /= 0) lnx = nx
! Subprogram not used     endif
! Subprogram not used     if (present(ny)) then
! Subprogram not used        if (ny /= 0) lny = ny
! Subprogram not used     endif
! Subprogram not used     if (lnx*lny /= ng) then
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: grid2d size not consistent ',ng,lnx,lny,trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))
! Subprogram not used 
! Subprogram not used        if (present(nt)) then
! Subprogram not used           dimid3(1:2) = dimid2
! Subprogram not used           rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
! Subprogram not used           dimid => dimid3
! Subprogram not used        else
! Subprogram not used           dimid => dimid2
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        do k = 1,nf
! Subprogram not used           call mct_aVect_getRList(mstring,k,AV)
! Subprogram not used           itemc = mct_string_toChar(mstring)
! Subprogram not used           call mct_string_clean(mstring)
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used           call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
! Subprogram not used 	  if (luse_float) then 
! Subprogram not used              rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
! Subprogram not used           else
! Subprogram not used              rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
! Subprogram not used           end if
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
! Subprogram not used           if (present(tavg)) then
! Subprogram not used              if (tavg) then
! Subprogram not used                 rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
! Subprogram not used        call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
! Subprogram not used        deallocate(dof)
! Subprogram not used 
! Subprogram not used        do k = 1,nf
! Subprogram not used           call mct_aVect_getRList(mstring,k,AV)
! Subprogram not used           itemc = mct_string_toChar(mstring)
! Subprogram not used           call mct_string_clean(mstring)
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used           rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
! Subprogram not used           if (present(nt)) then
! Subprogram not used              frame = nt
! Subprogram not used           else
! Subprogram not used              frame = 1
! Subprogram not used           endif
! Subprogram not used           call pio_setframe(varid,frame)
! Subprogram not used           call pio_write_darray(cpl_io_file, varid, iodesc, av%rattr(k,:), rcode, fillval=lfillvalue)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        call pio_freedecomp(cpl_io_file, iodesc)
! Subprogram not used 
! Subprogram not used     end if
! Subprogram not used   end subroutine seq_io_write_av

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_avs - write AVS to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_avs(filename,gsmap,AVS,dname,whead,wdata,nx,ny,nt,fillval,pre,tavg,&
! Subprogram not used 	                     use_float)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     type(mct_gsMap), intent(in) :: gsmap
! Subprogram not used     type(mct_aVect) ,intent(in) :: AVS(:)   ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used     integer(in),optional,intent(in) :: nx   ! 2d grid size if available
! Subprogram not used     integer(in),optional,intent(in) :: ny   ! 2d grid size if available
! Subprogram not used     integer(in),optional,intent(in) :: nt   ! time sample
! Subprogram not used     real(r8),optional,intent(in) :: fillval ! fill value
! Subprogram not used     character(len=*),optional,intent(in) :: pre      ! prefix to variable name
! Subprogram not used     logical,optional,intent(in) :: tavg     ! is this a tavg
! Subprogram not used     logical,optional,intent(in) :: use_float ! write output as float rather than double
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: mpicom
! Subprogram not used     integer(in) :: iam
! Subprogram not used     integer(in) :: nf,ns,ng,ni
! Subprogram not used     integer(in) :: i,j,k,n,k1,k2
! Subprogram not used     integer(in),target  :: dimid2(2)
! Subprogram not used     integer(in),target  :: dimid3(3)
! Subprogram not used     integer(in),target  :: dimid4(4)
! Subprogram not used     integer(in),pointer :: dimid(:)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used     integer(kind=PIO_OffSet) :: frame
! Subprogram not used     type(mct_string) :: mstring     ! mct char type
! Subprogram not used     character(CL)    :: itemc       ! string converted to char
! Subprogram not used     character(CL)    :: name1       ! var name
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     character(CL)    :: lpre        ! local prefix
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     logical :: luse_float
! Subprogram not used     integer(in) :: lnx,lny
! Subprogram not used     real(r8) :: lfillvalue
! Subprogram not used     real(r8), allocatable :: data(:)
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_avs) '
! Subprogram not used     integer :: lbnum
! Subprogram not used     integer, pointer :: Dof(:)
! Subprogram not used     integer, pointer :: Dofn(:)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lfillvalue = fillvalue
! Subprogram not used     if (present(fillval)) then
! Subprogram not used        lfillvalue = fillval
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lpre = trim(dname)
! Subprogram not used     if (present(pre)) then
! Subprogram not used        lpre = trim(pre)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     luse_float = .false.
! Subprogram not used     if (present(use_float)) luse_float = use_float
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)	
! Subprogram not used 
! Subprogram not used     ni = size(AVS)
! Subprogram not used 
! Subprogram not used     ns = mct_aVect_lsize(AVS(1))
! Subprogram not used     ng = mct_gsmap_gsize(gsmap)
! Subprogram not used     lnx = ng
! Subprogram not used     lny = 1
! Subprogram not used 	
! Subprogram not used     nf = mct_aVect_nRattr(AVS(1))
! Subprogram not used     if (nf < 1) then
! Subprogram not used        write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (present(nx)) then
! Subprogram not used        if (nx /= 0) lnx = nx
! Subprogram not used     endif
! Subprogram not used     if (present(ny)) then
! Subprogram not used        if (ny /= 0) lny = ny
! Subprogram not used     endif
! Subprogram not used     if (lnx*lny /= ng) then
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: grid2d size not consistent ',ng,lnx,lny,trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))
! Subprogram not used 
! Subprogram not used        if (ni > 1) then
! Subprogram not used           rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ni',ni,dimid3(3))
! Subprogram not used           if (present(nt)) then
! Subprogram not used              dimid4(1:2) = dimid2
! Subprogram not used              dimid4(3) = dimid3(3)
! Subprogram not used              rcode = pio_inq_dimid(cpl_io_file,'time',dimid4(4))
! Subprogram not used              dimid => dimid4
! Subprogram not used           else
! Subprogram not used              dimid3(1:2) = dimid2
! Subprogram not used              dimid => dimid3
! Subprogram not used           endif
! Subprogram not used        else
! Subprogram not used           if (present(nt)) then
! Subprogram not used              dimid3(1:2) = dimid2
! Subprogram not used              rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
! Subprogram not used              dimid => dimid3
! Subprogram not used           else
! Subprogram not used              dimid => dimid2
! Subprogram not used           endif
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        do k = 1,nf
! Subprogram not used           call mct_aVect_getRList(mstring,k,AVS(1))
! Subprogram not used           itemc = mct_string_toChar(mstring)
! Subprogram not used           call mct_string_clean(mstring)
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used           call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
! Subprogram not used 	  if (luse_float) then 
! Subprogram not used              rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
! Subprogram not used           else
! Subprogram not used              rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
! Subprogram not used           end if
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
! Subprogram not used           if (present(tavg)) then
! Subprogram not used              if (tavg) then
! Subprogram not used                 rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        allocate(data(ns*ni))
! Subprogram not used        ! note: size of dof is ns
! Subprogram not used        call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
! Subprogram not used        if (ni > 1) then
! Subprogram not used           allocate(dofn(ns*ni))
! Subprogram not used           n = 0
! Subprogram not used           do k1 = 1,ni
! Subprogram not used           do k2 = 1,ns
! Subprogram not used              n = n + 1
! Subprogram not used              dofn(n) = (k1-1)*ng + dof(k2)
! Subprogram not used           enddo
! Subprogram not used           enddo
! Subprogram not used           call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,ni/), dofn, iodesc)
! Subprogram not used           deallocate(dofn)
! Subprogram not used        else
! Subprogram not used           call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
! Subprogram not used        endif
! Subprogram not used        deallocate(dof)
! Subprogram not used 
! Subprogram not used        do k = 1,nf
! Subprogram not used           call mct_aVect_getRList(mstring,k,AVS(1))
! Subprogram not used           itemc = mct_string_toChar(mstring)
! Subprogram not used           call mct_string_clean(mstring)
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used           rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
! Subprogram not used           if (present(nt)) then
! Subprogram not used              frame = nt
! Subprogram not used           else
! Subprogram not used              frame = 1
! Subprogram not used           endif
! Subprogram not used           call pio_setframe(varid,frame)
! Subprogram not used           n = 0
! Subprogram not used           do k1 = 1,ni
! Subprogram not used           do k2 = 1,ns
! Subprogram not used              n = n + 1
! Subprogram not used              data(n) = AVS(k1)%rAttr(k,k2)
! Subprogram not used           enddo
! Subprogram not used           enddo
! Subprogram not used !          call pio_write_darray(cpl_io_file, varid, iodesc, av%rattr(k,:), rcode, fillval=lfillvalue)
! Subprogram not used           call pio_write_darray(cpl_io_file, varid, iodesc, data, rcode, fillval=lfillvalue)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        deallocate(data)
! Subprogram not used        call pio_freedecomp(cpl_io_file, iodesc)
! Subprogram not used 
! Subprogram not used     end if
! Subprogram not used   end subroutine seq_io_write_avs

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_avs - write AVS to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_avscomp(filename, comp, flow, dname, &
! Subprogram not used        whead, wdata, nx, ny, nt, fillval, pre, tavg, use_float)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*) ,intent(in)          :: filename  ! file
! Subprogram not used     type(component_type)  ,intent(in)          :: comp(:)   ! data to be written
! Subprogram not used     character(len=3) ,intent(in)          :: flow      ! 'c2x' or 'x2c'
! Subprogram not used     character(len=*) ,intent(in)          :: dname     ! name of data
! Subprogram not used     logical          ,optional,intent(in) :: whead     ! write header
! Subprogram not used     logical          ,optional,intent(in) :: wdata     ! write data
! Subprogram not used     integer(in)      ,optional,intent(in) :: nx        ! 2d grid size if available
! Subprogram not used     integer(in)      ,optional,intent(in) :: ny        ! 2d grid size if available
! Subprogram not used     integer(in)      ,optional,intent(in) :: nt        ! time sample
! Subprogram not used     real(r8)         ,optional,intent(in) :: fillval   ! fill value
! Subprogram not used     character(len=*) ,optional,intent(in) :: pre       ! prefix to variable name
! Subprogram not used     logical          ,optional,intent(in) :: tavg      ! is this a tavg
! Subprogram not used     logical          ,optional,intent(in) :: use_float ! write output as float rather than double
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     type(mct_gsMap), pointer :: gsmap     ! global seg map on coupler processes
! Subprogram not used     type(mct_avect), pointer :: avcomp1
! Subprogram not used     type(mct_avect), pointer :: avcomp
! Subprogram not used     integer(in)              :: rcode
! Subprogram not used     integer(in)              :: mpicom
! Subprogram not used     integer(in)              :: iam
! Subprogram not used     integer(in)              :: nf,ns,ng,ni
! Subprogram not used     integer(in)              :: i,j,k,n,k1,k2
! Subprogram not used     integer(in),target       :: dimid2(2)
! Subprogram not used     integer(in),target       :: dimid3(3)
! Subprogram not used     integer(in),target       :: dimid4(4)
! Subprogram not used     integer(in),pointer      :: dimid(:)
! Subprogram not used     type(var_desc_t)         :: varid
! Subprogram not used     type(io_desc_t)          :: iodesc
! Subprogram not used     integer(kind=PIO_OffSet) :: frame
! Subprogram not used     type(mct_string)         :: mstring     ! mct char type
! Subprogram not used     character(CL)            :: itemc       ! string converted to char
! Subprogram not used     character(CL)            :: name1       ! var name
! Subprogram not used     character(CL)            :: cunit       ! var units
! Subprogram not used     character(CL)            :: lname       ! long name
! Subprogram not used     character(CL)            :: sname       ! standard name
! Subprogram not used     character(CL)            :: lpre        ! local prefix
! Subprogram not used     logical                  :: exists
! Subprogram not used     logical                  :: lwhead, lwdata
! Subprogram not used     logical                  :: luse_float
! Subprogram not used     integer(in)              :: lnx,lny
! Subprogram not used     real(r8)                 :: lfillvalue
! Subprogram not used     real(r8), allocatable    :: data(:)
! Subprogram not used     character(*),parameter   :: subName = '(seq_io_write_avs) '
! Subprogram not used     integer                  :: lbnum
! Subprogram not used     integer, pointer         :: Dof(:)
! Subprogram not used     integer, pointer         :: Dofn(:)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lfillvalue = fillvalue
! Subprogram not used     if (present(fillval)) then
! Subprogram not used        lfillvalue = fillval
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lpre = trim(dname)
! Subprogram not used     if (present(pre)) then
! Subprogram not used        lpre = trim(pre)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     luse_float = .false.
! Subprogram not used     if (present(use_float)) luse_float = use_float
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)	
! Subprogram not used 
! Subprogram not used     ni = size(comp)
! Subprogram not used     if (trim(flow) == 'x2c') avcomp1 => component_get_x2c_cx(comp(1))  
! Subprogram not used     if (trim(flow) == 'c2x') avcomp1 => component_get_c2x_cx(comp(1))  
! Subprogram not used     gsmap => component_get_gsmap_cx(comp(1))  
! Subprogram not used     ns = mct_aVect_lsize(avcomp1)
! Subprogram not used     ng = mct_gsmap_gsize(gsmap)
! Subprogram not used     lnx = ng
! Subprogram not used     lny = 1
! Subprogram not used 	
! Subprogram not used     nf = mct_aVect_nRattr(avcomp1)
! Subprogram not used     if (nf < 1) then
! Subprogram not used        write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (present(nx)) then
! Subprogram not used        if (nx /= 0) lnx = nx
! Subprogram not used     endif
! Subprogram not used     if (present(ny)) then
! Subprogram not used        if (ny /= 0) lny = ny
! Subprogram not used     endif
! Subprogram not used     if (lnx*lny /= ng) then
! Subprogram not used        if(iam==0) then
! Subprogram not used           write(logunit,*) subname,' ERROR: grid2d size not consistent ',&
! Subprogram not used                ng,lnx,lny,trim(dname)
! Subprogram not used        end if
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))
! Subprogram not used 
! Subprogram not used        if (ni > 1) then
! Subprogram not used           rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ni',ni,dimid3(3))
! Subprogram not used           if (present(nt)) then
! Subprogram not used              dimid4(1:2) = dimid2
! Subprogram not used              dimid4(3) = dimid3(3)
! Subprogram not used              rcode = pio_inq_dimid(cpl_io_file,'time',dimid4(4))
! Subprogram not used              dimid => dimid4
! Subprogram not used           else
! Subprogram not used              dimid3(1:2) = dimid2
! Subprogram not used              dimid => dimid3
! Subprogram not used           endif
! Subprogram not used        else
! Subprogram not used           if (present(nt)) then
! Subprogram not used              dimid3(1:2) = dimid2
! Subprogram not used              rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
! Subprogram not used              dimid => dimid3
! Subprogram not used           else
! Subprogram not used              dimid => dimid2
! Subprogram not used           endif
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        do k = 1,nf
! Subprogram not used           call mct_aVect_getRList(mstring,k,avcomp1)
! Subprogram not used           itemc = mct_string_toChar(mstring)
! Subprogram not used           call mct_string_clean(mstring)
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used           call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
! Subprogram not used 	  if (luse_float) then 
! Subprogram not used              rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
! Subprogram not used           else
! Subprogram not used              rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
! Subprogram not used           end if
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
! Subprogram not used           if (present(tavg)) then
! Subprogram not used              if (tavg) then
! Subprogram not used                 rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        allocate(data(ns*ni))
! Subprogram not used        ! note: size of dof is ns
! Subprogram not used        call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
! Subprogram not used        if (ni > 1) then
! Subprogram not used           allocate(dofn(ns*ni))
! Subprogram not used           n = 0
! Subprogram not used           do k1 = 1,ni
! Subprogram not used           do k2 = 1,ns
! Subprogram not used              n = n + 1
! Subprogram not used              dofn(n) = (k1-1)*ng + dof(k2)
! Subprogram not used           enddo
! Subprogram not used           enddo
! Subprogram not used           call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,ni/), dofn, iodesc)
! Subprogram not used           deallocate(dofn)
! Subprogram not used        else
! Subprogram not used           call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
! Subprogram not used        endif
! Subprogram not used        deallocate(dof)
! Subprogram not used 
! Subprogram not used        do k = 1,nf
! Subprogram not used           call mct_aVect_getRList(mstring,k,avcomp1)
! Subprogram not used           itemc = mct_string_toChar(mstring)
! Subprogram not used           call mct_string_clean(mstring)
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used           rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
! Subprogram not used           if (present(nt)) then
! Subprogram not used              frame = nt
! Subprogram not used           else
! Subprogram not used              frame = 1
! Subprogram not used           endif
! Subprogram not used           call pio_setframe(varid,frame)
! Subprogram not used           n = 0
! Subprogram not used           do k1 = 1,ni
! Subprogram not used              if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(k1))  
! Subprogram not used              if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(k1))  
! Subprogram not used              do k2 = 1,ns
! Subprogram not used                 n = n + 1
! Subprogram not used                 data(n) = avcomp%rAttr(k,k2)
! Subprogram not used              enddo
! Subprogram not used           enddo
! Subprogram not used           call pio_write_darray(cpl_io_file, varid, iodesc, data, rcode, fillval=lfillvalue)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        deallocate(data)
! Subprogram not used        call pio_freedecomp(cpl_io_file, iodesc)
! Subprogram not used 
! Subprogram not used     end if
! Subprogram not used   end subroutine seq_io_write_avscomp

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_int - write scalar integer to netcdf file
  !
  ! !DESCRIPTION:
  !    Write scalar integer to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_int(filename,idata,dname,whead,wdata)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     integer(in)     ,intent(in) :: idata    ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_int) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
! Subprogram not used !       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',1,dimid(1))
! Subprogram not used !       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,dimid,varid)
! Subprogram not used        rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,varid)
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
! Subprogram not used        rcode = pio_put_var(cpl_io_file,varid,idata)
! Subprogram not used 
! Subprogram not used        !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_write_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_int1d - write 1d integer array to netcdf file
  !
  ! !DESCRIPTION:
  !    Write 1d integer array to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_int1d(filename,idata,dname,whead,wdata)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     integer(in)     ,intent(in) :: idata(:) ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam
! Subprogram not used     integer(in) :: dimid(1)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     integer(in) :: lnx
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_int1d) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
! Subprogram not used        lnx = size(idata)
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',lnx,dimid(1))
! Subprogram not used        rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,dimid,varid)
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
! Subprogram not used        rcode = pio_put_var(cpl_io_file,varid,idata)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used        !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_write_int1d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_r8 - write scalar double to netcdf file
  !
  ! !DESCRIPTION:
  !    Write scalar double to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_r8(filename,rdata,dname,whead,wdata)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     real(r8)        ,intent(in) :: rdata    ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_r8) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
! Subprogram not used !       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',1,dimid(1))
! Subprogram not used !       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,dimid,varid)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,varid)
! Subprogram not used        if(rcode==PIO_NOERR) then
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used           rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used           if (lwdata) call seq_io_enddef(filename)
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
! Subprogram not used        rcode = pio_put_var(cpl_io_file,varid,rdata)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_write_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_r81d - write 1d double array to netcdf file
  !
  ! !DESCRIPTION:
  !    Write 1d double array to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_r81d(filename,rdata,dname,whead,wdata)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     real(r8)        ,intent(in) :: rdata(:) ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: mpicom
! Subprogram not used     integer(in) :: iam
! Subprogram not used     integer(in) :: dimid(1)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     integer(in) :: lnx
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_r81d) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
! Subprogram not used        lnx = size(rdata)
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',lnx,dimid(1))
! Subprogram not used        rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,dimid,varid)
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
! Subprogram not used        rcode = pio_put_var(cpl_io_file,varid,rdata)
! Subprogram not used 
! Subprogram not used        !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_write_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_char - write char string to netcdf file
  !
  ! !DESCRIPTION:
  !    Write char string to netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_write_char(filename,rdata,dname,whead,wdata)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     character(len=*),intent(in) :: rdata    ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     logical,optional,intent(in) :: whead    ! write header
! Subprogram not used     logical,optional,intent(in) :: wdata    ! write data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: mpicom
! Subprogram not used     integer(in) :: iam
! Subprogram not used     integer(in) :: dimid(1)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     character(CL)    :: cunit       ! var units
! Subprogram not used     character(CL)    :: lname       ! long name
! Subprogram not used     character(CL)    :: sname       ! standard name
! Subprogram not used     integer(in) :: lnx
! Subprogram not used     logical :: exists
! Subprogram not used     logical :: lwhead, lwdata
! Subprogram not used     character(*),parameter :: subName = '(seq_io_write_char) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lwhead = .true.
! Subprogram not used     lwdata = .true.
! Subprogram not used     if (present(whead)) lwhead = whead
! Subprogram not used     if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used     if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used        ! should we write a warning?
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used     if (lwhead) then
! Subprogram not used        call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
! Subprogram not used        lnx = len(charvar)
! Subprogram not used        rcode = pio_def_dim(cpl_io_file,trim(dname)//'_len',lnx,dimid(1))
! Subprogram not used        rcode = pio_def_var(cpl_io_file,trim(dname),PIO_CHAR,dimid,varid)
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
! Subprogram not used        rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
! Subprogram not used        if (lwdata) call seq_io_enddef(filename)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (lwdata) then
! Subprogram not used        charvar = ''
! Subprogram not used        charvar = trim(rdata)
! Subprogram not used        rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
! Subprogram not used        rcode = pio_put_var(cpl_io_file,varid,charvar)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_write_char

  !===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_io_write_time - write time variable to netcdf file
!
! !DESCRIPTION:
!    Write time variable to netcdf file
!
! !REVISION HISTORY:
!    2009-Feb-11 - M. Vertenstein - initial version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine seq_io_write_time(filename,time_units,time_cal,time_val,nt,whead,wdata,tbnds)
! Subprogram not used 
! Subprogram not used    use shr_cal_mod, only : shr_cal_calMaxLen, shr_cal_calendarName, &
! Subprogram not used                            shr_cal_noleap, shr_cal_gregorian
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used    implicit none
! Subprogram not used    character(len=*),intent(in) :: filename      ! file
! Subprogram not used    character(len=*),intent(in) :: time_units    ! units of time
! Subprogram not used    character(len=*),intent(in) :: time_cal      ! calendar type
! Subprogram not used    real(r8)        ,intent(in) :: time_val      ! data to be written
! Subprogram not used    integer(in),optional,intent(in) :: nt
! Subprogram not used    logical,optional,intent(in) :: whead         ! write header
! Subprogram not used    logical,optional,intent(in) :: wdata         ! write data
! Subprogram not used    real(r8),optional,intent(in) :: tbnds(2)     ! time bounds
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    integer(in) :: rcode
! Subprogram not used    integer(in) :: iam
! Subprogram not used    integer(in) :: dimid(1)
! Subprogram not used    integer(in) :: dimid2(2)
! Subprogram not used    type(var_desc_t) :: varid
! Subprogram not used    integer(in) :: lnx
! Subprogram not used    logical :: exists
! Subprogram not used    logical :: lwhead, lwdata
! Subprogram not used    integer :: start(4),count(4)
! Subprogram not used    character(len=shr_cal_calMaxLen) :: lcalendar
! Subprogram not used    real(r8) :: time_val_1d(1)
! Subprogram not used    character(*),parameter :: subName = '(seq_io_write_time) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    lwhead = .true.
! Subprogram not used    lwdata = .true.
! Subprogram not used    if (present(whead)) lwhead = whead
! Subprogram not used    if (present(wdata)) lwdata = wdata
! Subprogram not used 
! Subprogram not used    if (.not.lwhead .and. .not.lwdata) then
! Subprogram not used       ! should we write a warning?
! Subprogram not used       return
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    call seq_comm_setptrs(CPLID,iam=iam)
! Subprogram not used 
! Subprogram not used    if (lwhead) then
! Subprogram not used       rcode = pio_def_dim(cpl_io_file,'time',PIO_UNLIMITED,dimid(1))
! Subprogram not used       rcode = pio_def_var(cpl_io_file,'time',PIO_DOUBLE,dimid,varid)
! Subprogram not used       rcode = pio_put_att(cpl_io_file,varid,'units',trim(time_units))
! Subprogram not used       lcalendar = shr_cal_calendarName(time_cal,trap=.false.)
! Subprogram not used       if (trim(lcalendar) == trim(shr_cal_noleap)) then
! Subprogram not used          lcalendar = 'noleap'
! Subprogram not used       elseif (trim(lcalendar) == trim(shr_cal_gregorian)) then
! Subprogram not used          lcalendar = 'gregorian'
! Subprogram not used       endif
! Subprogram not used       rcode = pio_put_att(cpl_io_file,varid,'calendar',trim(lcalendar))
! Subprogram not used       if (present(tbnds)) then
! Subprogram not used          rcode = pio_put_att(cpl_io_file,varid,'bounds','time_bnds')
! Subprogram not used          dimid2(2)=dimid(1)
! Subprogram not used          rcode = pio_def_dim(cpl_io_file,'ntb',2,dimid2(1))
! Subprogram not used          rcode = pio_def_var(cpl_io_file,'time_bnds',PIO_DOUBLE,dimid2,varid)
! Subprogram not used       endif
! Subprogram not used       if (lwdata) call seq_io_enddef(filename)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (lwdata) then
! Subprogram not used       start = 1
! Subprogram not used       count = 1
! Subprogram not used       if (present(nt)) then
! Subprogram not used          start(1) = nt
! Subprogram not used       endif
! Subprogram not used       time_val_1d(1) = time_val
! Subprogram not used       rcode = pio_inq_varid(cpl_io_file,'time',varid)
! Subprogram not used       rcode = pio_put_var(cpl_io_file,varid,start,count,time_val_1d)
! Subprogram not used       if (present(tbnds)) then
! Subprogram not used          rcode = pio_inq_varid(cpl_io_file,'time_bnds',varid)
! Subprogram not used          start = 1
! Subprogram not used          count = 1
! Subprogram not used          if (present(nt)) then
! Subprogram not used             start(2) = nt
! Subprogram not used          endif
! Subprogram not used          count(1) = 2
! Subprogram not used          rcode = pio_put_var(cpl_io_file,varid,start,count,tbnds)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       !      write(logunit,*) subname,' wrote time ',lwhead,lwdata
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used  end subroutine seq_io_write_time

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_av - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_av(filename,gsmap,AV,dname,pre)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     type(mct_gsMap), intent(in) :: gsmap
! Subprogram not used     type(mct_aVect) ,intent(inout):: AV     ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     character(len=*),intent(in),optional :: pre      ! prefix name
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam,mpicom
! Subprogram not used     integer(in) :: nf,ns,ng
! Subprogram not used     integer(in) :: i,j,k,n, ndims
! Subprogram not used     type(file_desc_t) :: pioid
! Subprogram not used     integer(in) :: dimid(2)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     integer(in) :: lnx,lny
! Subprogram not used     type(mct_string) :: mstring     ! mct char type
! Subprogram not used     character(CL)    :: itemc       ! string converted to char
! Subprogram not used     logical :: exists
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used     integer(in), pointer :: dof(:)
! Subprogram not used     character(CL)  :: lversion
! Subprogram not used     character(CL)  :: name1
! Subprogram not used     character(CL)  :: lpre
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_av) '
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lversion = trim(version0)
! Subprogram not used 
! Subprogram not used     lpre = trim(dname)
! Subprogram not used     if (present(pre)) then
! Subprogram not used        lpre = trim(pre)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used     call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
! Subprogram not used 
! Subprogram not used     ns = mct_aVect_lsize(AV)
! Subprogram not used     nf = mct_aVect_nRattr(AV)
! Subprogram not used 
! Subprogram not used     if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used     call shr_mpi_bcast(exists,mpicom,'seq_io_read_av exists')
! Subprogram not used     if (exists) then
! Subprogram not used        rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
! Subprogram not used        if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     else
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     do k = 1,nf
! Subprogram not used        call mct_aVect_getRList(mstring,k,AV)
! Subprogram not used        itemc = mct_string_toChar(mstring)
! Subprogram not used        call mct_string_clean(mstring)
! Subprogram not used        if (trim(lversion) == trim(version)) then
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used        else
! Subprogram not used           name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
! Subprogram not used        endif
! Subprogram not used        call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_inq_varid(pioid,trim(name1),varid)
! Subprogram not used        if (rcode == pio_noerr) then
! Subprogram not used           if (k==1) then
! Subprogram not used              rcode = pio_inq_varndims(pioid, varid, ndims)
! Subprogram not used              rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
! Subprogram not used              rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
! Subprogram not used              if (ndims>=2) then
! Subprogram not used                 rcode = pio_inq_dimlen(pioid, dimid(2), lny)
! Subprogram not used              else
! Subprogram not used                 lny = 1
! Subprogram not used              end if
! Subprogram not used              ng = lnx * lny
! Subprogram not used              if (ng /= mct_gsmap_gsize(gsmap)) then
! Subprogram not used                 if (iam==0) write(logunit,*) subname,' ERROR: dimensions do not match',&
! Subprogram not used                      lnx,lny,mct_gsmap_gsize(gsmap)
! Subprogram not used                 call shr_sys_abort()
! Subprogram not used              end if
! Subprogram not used              call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
! Subprogram not used              deallocate(dof)
! Subprogram not used           end if
! Subprogram not used           call pio_read_darray(pioid,varid,iodesc, av%rattr(k,:), rcode)
! Subprogram not used        else
! Subprogram not used           write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
! Subprogram not used           write(logunit,*)'for backwards compatibility will set it to 0'
! Subprogram not used           av%rattr(k,:) = 0.0_r8
! Subprogram not used        end if
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used           
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     !--- zero out fill value, this is somewhat arbitrary
! Subprogram not used     do n = 1,ns
! Subprogram not used     do k = 1,nf
! Subprogram not used        if (AV%rAttr(k,n) == fillvalue) then
! Subprogram not used            AV%rAttr(k,n) = 0.0_r8
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call pio_freedecomp(pioid, iodesc)
! Subprogram not used     call pio_closefile(pioid)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_av

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_avs - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_avs(filename,gsmap,AVS,dname,pre)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     type(mct_gsMap), intent(in) :: gsmap
! Subprogram not used     type(mct_aVect) ,intent(inout):: AVS(:) ! data to be written
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used     character(len=*),intent(in),optional :: pre      ! prefix name
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam,mpicom
! Subprogram not used     integer(in) :: nf,ns,ng,ni
! Subprogram not used     integer(in) :: i,j,k,n,n1,n2,ndims
! Subprogram not used     type(file_desc_t) :: pioid
! Subprogram not used     integer(in) :: dimid(4)
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     integer(in) :: lnx,lny,lni
! Subprogram not used     type(mct_string) :: mstring     ! mct char type
! Subprogram not used     character(CL)    :: itemc       ! string converted to char
! Subprogram not used     logical :: exists
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used     integer(in), pointer :: dof(:)
! Subprogram not used     integer(in), pointer :: dofn(:)
! Subprogram not used     real(r8), allocatable :: data(:)
! Subprogram not used     character(CL)  :: lversion
! Subprogram not used     character(CL)  :: name1
! Subprogram not used     character(CL)  :: lpre
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_avs) '
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lversion = trim(version0)
! Subprogram not used 
! Subprogram not used     lpre = trim(dname)
! Subprogram not used     if (present(pre)) then
! Subprogram not used        lpre = trim(pre)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used     call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
! Subprogram not used 
! Subprogram not used     ni = size(AVS)
! Subprogram not used     ns = mct_aVect_lsize(AVS(1))
! Subprogram not used     nf = mct_aVect_nRattr(AVS(1))
! Subprogram not used     ng = mct_gsmap_gsize(gsmap)
! Subprogram not used 
! Subprogram not used     if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used     call shr_mpi_bcast(exists,mpicom,'seq_io_read_avs exists')
! Subprogram not used     if (exists) then
! Subprogram not used        rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
! Subprogram not used        if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     else
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     allocate(data(ni*ns))
! Subprogram not used 
! Subprogram not used     do k = 1,nf
! Subprogram not used        call mct_aVect_getRList(mstring,k,AVS(1))
! Subprogram not used        itemc = mct_string_toChar(mstring)
! Subprogram not used        call mct_string_clean(mstring)
! Subprogram not used        if (trim(lversion) == trim(version)) then
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used        else
! Subprogram not used           name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
! Subprogram not used        endif
! Subprogram not used        call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_inq_varid(pioid,trim(name1),varid)
! Subprogram not used        if (rcode == pio_noerr) then
! Subprogram not used           if (k==1) then
! Subprogram not used              rcode = pio_inq_varndims(pioid, varid, ndims)
! Subprogram not used              rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
! Subprogram not used              rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
! Subprogram not used              if (ndims>=2) then
! Subprogram not used                 rcode = pio_inq_dimlen(pioid, dimid(2), lny)
! Subprogram not used              else
! Subprogram not used                 lny = 1
! Subprogram not used              end if
! Subprogram not used              if (lnx*lny /= ng) then
! Subprogram not used                 write(logunit,*) subname,' ERROR: dimensions do not match',&
! Subprogram not used                      lnx,lny,mct_gsmap_gsize(gsmap)
! Subprogram not used                 call shr_sys_abort()
! Subprogram not used              end if
! Subprogram not used              if (ndims>=3) then
! Subprogram not used                 rcode = pio_inq_dimlen(pioid, dimid(3), lni)
! Subprogram not used              else
! Subprogram not used                 lni = 1
! Subprogram not used              end if
! Subprogram not used              if (ni /= lni) then
! Subprogram not used                 write(logunit,*) subname,' ERROR: ni dimensions do not match',ni,lni
! Subprogram not used                 call shr_sys_abort()
! Subprogram not used              end if
! Subprogram not used              if (ni > 1) then
! Subprogram not used                 allocate(dofn(ns*ni))
! Subprogram not used                 n = 0
! Subprogram not used                 do n1 = 1,ni
! Subprogram not used                 do n2 = 1,ns
! Subprogram not used                    n = n + 1
! Subprogram not used                    dofn(n) = (n1-1)*ng + dof(n2)
! Subprogram not used                 enddo
! Subprogram not used                 enddo
! Subprogram not used                 call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,lni/), dofn, iodesc)
! Subprogram not used                 deallocate(dofn)
! Subprogram not used              else
! Subprogram not used                 call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
! Subprogram not used              endif
! Subprogram not used              deallocate(dof)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           call pio_read_darray(pioid,varid,iodesc, data, rcode)
! Subprogram not used           n = 0
! Subprogram not used           do n1 = 1,ni
! Subprogram not used           do n2 = 1,ns
! Subprogram not used              n = n + 1
! Subprogram not used              avs(n1)%rAttr(k,n2) = data(n)
! Subprogram not used           enddo
! Subprogram not used           enddo
! Subprogram not used        else
! Subprogram not used           write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
! Subprogram not used           write(logunit,*)'for backwards compatibility will set it to 0'
! Subprogram not used           do n1 = 1,ni
! Subprogram not used              avs(n1)%rattr(k,:) = 0.0_r8
! Subprogram not used           enddo
! Subprogram not used        end if
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     deallocate(data)
! Subprogram not used 
! Subprogram not used     !--- zero out fill value, this is somewhat arbitrary
! Subprogram not used     do n1 = 1,ni
! Subprogram not used     do n2 = 1,ns
! Subprogram not used     do k = 1,nf
! Subprogram not used        if (AVS(n1)%rAttr(k,n2) == fillvalue) then
! Subprogram not used            AVS(n1)%rAttr(k,n2) = 0.0_r8
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call pio_freedecomp(pioid, iodesc)
! Subprogram not used     call pio_closefile(pioid)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_avs

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_avs - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_avscomp(filename, comp, flow, dname, pre)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in)          :: filename ! file
! Subprogram not used     type(component_type),  intent(inout)       :: comp(:)
! Subprogram not used     character(len=3), intent(in)          :: flow     ! 'c2x' or 'x2c'
! Subprogram not used     character(len=*), intent(in)          :: dname    ! name of data
! Subprogram not used     character(len=*), intent(in),optional :: pre      ! prefix name
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     type(mct_gsMap), pointer :: gsmap
! Subprogram not used     type(mct_aVect), pointer :: avcomp  
! Subprogram not used     type(mct_aVect), pointer :: avcomp1 
! Subprogram not used     integer(in)              :: rcode
! Subprogram not used     integer(in)              :: iam,mpicom
! Subprogram not used     integer(in)              :: nf,ns,ng,ni
! Subprogram not used     integer(in)              :: i,j,k,n,n1,n2,ndims
! Subprogram not used     type(file_desc_t)        :: pioid
! Subprogram not used     integer(in)              :: dimid(4)
! Subprogram not used     type(var_desc_t)         :: varid
! Subprogram not used     integer(in)              :: lnx,lny,lni
! Subprogram not used     type(mct_string)         :: mstring ! mct char type
! Subprogram not used     character(CL)            :: itemc   ! string converted to char
! Subprogram not used     logical                  :: exists
! Subprogram not used     type(io_desc_t)          :: iodesc
! Subprogram not used     integer(in), pointer     :: dof(:)
! Subprogram not used     integer(in), pointer     :: dofn(:)
! Subprogram not used     real(r8), allocatable    :: data(:)
! Subprogram not used     character(CL)            :: lversion
! Subprogram not used     character(CL)            :: name1
! Subprogram not used     character(CL)            :: lpre
! Subprogram not used     character(*),parameter   :: subName = '(seq_io_read_avs) '
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     lversion = trim(version0)
! Subprogram not used 
! Subprogram not used     lpre = trim(dname)
! Subprogram not used     if (present(pre)) then
! Subprogram not used        lpre = trim(pre)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     gsmap => component_get_gsmap_cx(comp(1))
! Subprogram not used     if (trim(flow) == 'x2c') avcomp1 => component_get_x2c_cx(comp(1))  
! Subprogram not used     if (trim(flow) == 'c2x') avcomp1 => component_get_c2x_cx(comp(1))  
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used     call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
! Subprogram not used 
! Subprogram not used     ni = size(comp)
! Subprogram not used     ns = mct_aVect_lsize(avcomp1)
! Subprogram not used     nf = mct_aVect_nRattr(avcomp1)
! Subprogram not used     ng = mct_gsmap_gsize(gsmap)
! Subprogram not used 
! Subprogram not used     if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used     call shr_mpi_bcast(exists,mpicom,'seq_io_read_avs exists')
! Subprogram not used     if (exists) then
! Subprogram not used        rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
! Subprogram not used        if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     else
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     allocate(data(ni*ns))
! Subprogram not used 
! Subprogram not used     do k = 1,nf
! Subprogram not used        call mct_aVect_getRList(mstring,k,avcomp1)
! Subprogram not used        itemc = mct_string_toChar(mstring)
! Subprogram not used        call mct_string_clean(mstring)
! Subprogram not used        if (trim(lversion) == trim(version)) then
! Subprogram not used           name1 = trim(lpre)//'_'//trim(itemc)
! Subprogram not used        else
! Subprogram not used           name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
! Subprogram not used        endif
! Subprogram not used        call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_inq_varid(pioid,trim(name1),varid)
! Subprogram not used        if (rcode == pio_noerr) then
! Subprogram not used           if (k==1) then
! Subprogram not used              rcode = pio_inq_varndims(pioid, varid, ndims)
! Subprogram not used              rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
! Subprogram not used              rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
! Subprogram not used              if (ndims>=2) then
! Subprogram not used                 rcode = pio_inq_dimlen(pioid, dimid(2), lny)
! Subprogram not used              else
! Subprogram not used                 lny = 1
! Subprogram not used              end if
! Subprogram not used              if (lnx*lny /= ng) then
! Subprogram not used                 write(logunit,*) subname,' ERROR: dimensions do not match',&
! Subprogram not used                      lnx,lny,mct_gsmap_gsize(gsmap)
! Subprogram not used                 call shr_sys_abort()
! Subprogram not used              end if
! Subprogram not used              if (ndims>=3) then
! Subprogram not used                 rcode = pio_inq_dimlen(pioid, dimid(3), lni)
! Subprogram not used              else
! Subprogram not used                 lni = 1
! Subprogram not used              end if
! Subprogram not used              if (ni /= lni) then
! Subprogram not used                 write(logunit,*) subname,' ERROR: ni dimensions do not match',ni,lni
! Subprogram not used                 call shr_sys_abort()
! Subprogram not used              end if
! Subprogram not used              if (ni > 1) then
! Subprogram not used                 allocate(dofn(ns*ni))
! Subprogram not used                 n = 0
! Subprogram not used                 do n1 = 1,ni
! Subprogram not used                 do n2 = 1,ns
! Subprogram not used                    n = n + 1
! Subprogram not used                    dofn(n) = (n1-1)*ng + dof(n2)
! Subprogram not used                 enddo
! Subprogram not used                 enddo
! Subprogram not used                 call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,lni/), dofn, iodesc)
! Subprogram not used                 deallocate(dofn)
! Subprogram not used              else
! Subprogram not used                 call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
! Subprogram not used              endif
! Subprogram not used              deallocate(dof)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           call pio_read_darray(pioid,varid,iodesc, data, rcode)
! Subprogram not used           n = 0
! Subprogram not used           do n1 = 1,ni
! Subprogram not used              if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(n1))  
! Subprogram not used              if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(n1))  
! Subprogram not used              do n2 = 1,ns
! Subprogram not used                 n = n + 1
! Subprogram not used                 avcomp%rAttr(k,n2) = data(n)
! Subprogram not used              enddo
! Subprogram not used           enddo
! Subprogram not used        else
! Subprogram not used           write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
! Subprogram not used           write(logunit,*)'for backwards compatibility will set it to 0'
! Subprogram not used           do n1 = 1,ni
! Subprogram not used              if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(n1))  
! Subprogram not used              if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(n1))  
! Subprogram not used              avcomp%rattr(k,:) = 0.0_r8
! Subprogram not used           enddo
! Subprogram not used        end if
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     deallocate(data)
! Subprogram not used 
! Subprogram not used     !--- zero out fill value, this is somewhat arbitrary
! Subprogram not used     do n1 = 1,ni
! Subprogram not used        if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(n1))  
! Subprogram not used        if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(n1))  
! Subprogram not used        do n2 = 1,ns
! Subprogram not used           do k = 1,nf
! Subprogram not used              if (avcomp%rAttr(k,n2) == fillvalue) then
! Subprogram not used                 avcomp%rAttr(k,n2) = 0.0_r8
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call pio_freedecomp(pioid, iodesc)
! Subprogram not used     call pio_closefile(pioid)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_avscomp

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int - read scalar integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar integer from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_int(filename,idata,dname)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     integer         ,intent(inout):: idata  ! integer data
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer :: i1d(1)
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_int) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_io_read_int1d(filename,i1d,dname)
! Subprogram not used     idata = i1d(1)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int1d - read 1d integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d integer array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_int1d(filename,idata,dname)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     integer(in)     ,intent(inout):: idata(:)  ! integer data
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam,mpicom
! Subprogram not used     type(file_desc_t) :: pioid 
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     logical :: exists
! Subprogram not used     character(CL)  :: lversion
! Subprogram not used     character(CL)  :: name1
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_int1d) '
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used     lversion=trim(version0)
! Subprogram not used 
! Subprogram not used     if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used     call shr_mpi_bcast(exists,mpicom,'seq_io_read_int1d exists')
! Subprogram not used     if (exists) then
! Subprogram not used        rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
! Subprogram not used        !         write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     else
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (trim(lversion) == trim(version)) then
! Subprogram not used        name1 = trim(dname)
! Subprogram not used     else
! Subprogram not used        name1 = trim(prefix)//trim(dname)
! Subprogram not used     endif
! Subprogram not used     rcode = pio_inq_varid(pioid,trim(name1),varid)
! Subprogram not used     rcode = pio_get_var(pioid,varid,idata)
! Subprogram not used 
! Subprogram not used     call pio_closefile(pioid)
! Subprogram not used 
! Subprogram not used     !      write(logunit,*) subname,' read int ',trim(dname)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_int1d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r8 - read scalar double from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar double from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_r8(filename,rdata,dname)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     real(r8)        ,intent(inout):: rdata  ! real data
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     real(r8) :: r1d(1)
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_r8) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_io_read_r81d(filename,r1d,dname)
! Subprogram not used     rdata = r1d(1)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r81d - read 1d double array from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d double array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_r81d(filename,rdata,dname)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     real(r8)        ,intent(inout):: rdata(:)  ! real data
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam,mpicom
! Subprogram not used     type(file_desc_T) :: pioid 
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     logical :: exists
! Subprogram not used     character(CL)  :: lversion
! Subprogram not used     character(CL)  :: name1
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_r81d) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used 
! Subprogram not used     lversion=trim(version0)
! Subprogram not used 
! Subprogram not used     if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used     call shr_mpi_bcast(exists,mpicom,'seq_io_read_r81d exists')
! Subprogram not used     if (exists) then
! Subprogram not used        rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
! Subprogram not used        !         write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     else
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (trim(lversion) == trim(version)) then
! Subprogram not used        name1 = trim(dname)
! Subprogram not used     else
! Subprogram not used        name1 = trim(prefix)//trim(dname)
! Subprogram not used     endif
! Subprogram not used     rcode = pio_inq_varid(pioid,trim(name1),varid)
! Subprogram not used     rcode = pio_get_var(pioid,varid,rdata)
! Subprogram not used 
! Subprogram not used     call pio_closefile(pioid)
! Subprogram not used 
! Subprogram not used     !      write(logunit,*) subname,' read int ',trim(dname)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_char - read char string from netcdf file
  !
  ! !DESCRIPTION:
  !    Read char string from netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used   subroutine seq_io_read_char(filename,rdata,dname)
! Subprogram not used 
! Subprogram not used     ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used     implicit none
! Subprogram not used     character(len=*),intent(in) :: filename ! file
! Subprogram not used     character(len=*),intent(inout):: rdata  ! character data
! Subprogram not used     character(len=*),intent(in) :: dname    ! name of data
! Subprogram not used 
! Subprogram not used     !EOP
! Subprogram not used 
! Subprogram not used     integer(in) :: rcode
! Subprogram not used     integer(in) :: iam,mpicom
! Subprogram not used     type(file_desc_T) :: pioid 
! Subprogram not used     type(var_desc_t) :: varid
! Subprogram not used     logical :: exists
! Subprogram not used     character(CL)  :: lversion
! Subprogram not used     character(CL)  :: name1
! Subprogram not used     character(*),parameter :: subName = '(seq_io_read_char) '
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
! Subprogram not used 
! Subprogram not used     lversion=trim(version0)
! Subprogram not used 
! Subprogram not used     if (iam==0) inquire(file=trim(filename),exist=exists)
! Subprogram not used     call shr_mpi_bcast(exists,mpicom,'seq_io_read_char exists')
! Subprogram not used     if (exists) then
! Subprogram not used        rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
! Subprogram not used        !         write(logunit,*) subname,' open file ',trim(filename)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
! Subprogram not used        rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
! Subprogram not used        call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
! Subprogram not used     else
! Subprogram not used        if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (trim(lversion) == trim(version)) then
! Subprogram not used        name1 = trim(dname)
! Subprogram not used     else
! Subprogram not used        name1 = trim(prefix)//trim(dname)
! Subprogram not used     endif
! Subprogram not used     rcode = pio_inq_varid(pioid,trim(name1),varid)
! Subprogram not used     rcode = pio_get_var(pioid,varid,charvar)
! Subprogram not used     rdata = trim(charvar)
! Subprogram not used 
! Subprogram not used     call pio_closefile(pioid)
! Subprogram not used 
! Subprogram not used   end subroutine seq_io_read_char

  !===============================================================================
!===============================================================================
end module seq_io_mod
