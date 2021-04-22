module RtmIO

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmIO
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files for RTM
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use shr_sys_mod    , only : shr_sys_flush, shr_sys_abort
  use shr_file_mod   , only : shr_file_getunit, shr_file_freeunit
  use RtmFileUtils   , only : getavu, relavu
  use RtmSpmd        , only : masterproc, mpicom_rof, iam, npes
  use RunoffMod      , only : runoff
  use RtmVar         , only : spval, ispval, iulog, inst_name
  use perf_mod       , only : t_startf, t_stopf
  use mct_mod
  use pio

! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: check_var          ! determine if variable is on netcdf file
  public :: check_dim          ! validity check on dimension
  public :: ncd_pio_openfile   ! open a file
  public :: ncd_pio_createfile ! create a new file
  public :: ncd_pio_closefile  ! close a file
  public :: ncd_pio_init       ! called from rtm_comp
  public :: ncd_enddef         ! end define mode
  public :: ncd_putatt         ! put attribute
  public :: ncd_defdim         ! define dimension
  public :: ncd_inqdid         ! inquire dimension id
  public :: ncd_inqdname       ! inquire dimension name
  public :: ncd_inqdlen        ! inquire dimension length
  public :: ncd_inqfdims       ! inquire file dimnesions 
  public :: ncd_defvar         ! define variables
  public :: ncd_inqvid         ! inquire variable id
  public :: ncd_inqvname       ! inquire variable name
  public :: ncd_inqvdims       ! inquire variable ndims
  public :: ncd_inqvdids       ! inquire variable dimids
  public :: ncd_io             ! write local data
  public :: ncd_finalize       ! clean up any memory

  integer,parameter,public :: ncd_int       = pio_int
  integer,parameter,public :: ncd_log       =-pio_int
  integer,parameter,public :: ncd_float     = pio_real
  integer,parameter,public :: ncd_double    = pio_double
  integer,parameter,public :: ncd_char      = pio_char
  integer,parameter,public :: ncd_global    = pio_global
  integer,parameter,public :: ncd_write     = pio_write
  integer,parameter,public :: ncd_nowrite   = pio_nowrite
  integer,parameter,public :: ncd_clobber   = pio_clobber
  integer,parameter,public :: ncd_noclobber = pio_noclobber
  integer,parameter,public :: ncd_nofill    = pio_nofill
  integer,parameter,public :: ncd_unlimited = pio_unlimited

  ! PIO types needed for ncdio_pio interface calls
  public file_desc_t
  public var_desc_t
  public io_desc_t
!
! !REVISION HISTORY:
!
!
! !PRIVATE MEMBER FUNCTIONS:
!

  interface ncd_putatt
     module procedure ncd_putatt_int
     module procedure ncd_putatt_real
     module procedure ncd_putatt_char
  end interface

  interface ncd_defvar
     module procedure ncd_defvar_bynf
     module procedure ncd_defvar_bygrid
  end interface

  interface ncd_io 
     ! global scalar
     module procedure ncd_io_log_var0_nf
     module procedure ncd_io_int_var0_nf
     module procedure ncd_io_real_var0_nf

     ! global 1d
     module procedure ncd_io_log_var1_nf
     module procedure ncd_io_int_var1_nf
     module procedure ncd_io_real_var1_nf
     module procedure ncd_io_char_var1_nf
     module procedure ncd_io_char_varn_strt_nf

     ! global 2d
     module procedure ncd_io_int_var2_nf
     module procedure ncd_io_real_var2_nf
     module procedure ncd_io_char_var2_nf

     ! local 1d
     module procedure ncd_io_log_var1
     module procedure ncd_io_int_var1
     module procedure ncd_io_real_var1
  end interface

  private :: ncd_getiodesc      ! obtain iodesc

  integer,parameter,private :: debug = 0             ! local debug level

  integer , parameter  , public  :: max_string_len = 256     ! length of strings
  real(r8), parameter  , public  :: fillvalue = 1.e36_r8     ! fill value for netcdf fields

  integer, public :: io_type

  type(iosystem_desc_t), pointer, public  :: pio_subsystem

  type iodesc_plus_type
     character(len=64) :: name
     type(IO_desc_t)   :: iodesc
     integer           :: type
     integer           :: ndims
     integer           :: dims(4)
     integer           :: dimids(4) 
  end type iodesc_plus_type

  integer,parameter      ,private :: max_iodesc = 100
  integer                ,private :: num_iodesc = 0
  type(iodesc_plus_type) ,private, target :: iodesc_list(max_iodesc)

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

  subroutine ncd_pio_init()

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initial PIO
    !
    ! !USES:
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    ! !ARGUMENTS:
    implicit none
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='ncd_pio_init' ! subroutine name
    !-----------------------------------------------------------------------

    PIO_subsystem => shr_pio_getiosys(inst_name)
    io_type       =  shr_pio_getiotype(inst_name)

  end subroutine ncd_pio_init

!-----------------------------------------------------------------------

  subroutine ncd_pio_openfile(file, fname, mode)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Open a NetCDF PIO file
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: file   ! Output PIO file handle
    character(len=*) , intent(in)    :: fname  ! Input filename to open
    integer          , intent(in)    :: mode   ! file mode
    ! !LOCAL VARIABLES:
    integer :: ierr
    character(len=*),parameter :: subname='ncd_pio_openfile' ! subroutine name
    !-----------------------------------------------------------------------

    ierr = pio_openfile(pio_subsystem, file, io_type, fname, mode)

    if(ierr/= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: Failed to open file')
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened existing file ', trim(fname), file%fh
    end if

  end subroutine ncd_pio_openfile

!-----------------------------------------------------------------------

  subroutine ncd_pio_closefile(file)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Close a NetCDF PIO file
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: file   ! PIO file handle to close
    !-----------------------------------------------------------------------

    call pio_closefile(file)

  end subroutine ncd_pio_closefile

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_pio_createfile(file, fname)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Create a new NetCDF file with PIO
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: file    ! PIO file descriptor
! Subprogram not used     character(len=*),  intent(in)    :: fname   ! File name to create
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: ierr
! Subprogram not used     character(len=*),parameter :: subname='ncd_pio_createfile' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ierr = pio_createfile(pio_subsystem, file, io_type, fname, ior(PIO_CLOBBER,PIO_64BIT_OFFSET))
! Subprogram not used 
! Subprogram not used     if(ierr/= PIO_NOERR) then
! Subprogram not used        call shr_sys_abort( subname//' ERROR: Failed to open file to write: '//trim(fname))
! Subprogram not used     else if(pio_subsystem%io_rank==0) then
! Subprogram not used        write(iulog,*) 'Opened file ', trim(fname),  ' to write', file%fh
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine ncd_pio_createfile

!-----------------------------------------------------------------------

! Subprogram not used   subroutine check_var(ncid, varname, vardesc, readvar, print_err )
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Check if variable is on netcdf file
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout)  :: ncid      ! PIO file descriptor
! Subprogram not used     character(len=*) , intent(in)     :: varname   ! Varible name to check
! Subprogram not used     type(Var_desc_t) , intent(out)    :: vardesc   ! Output variable descriptor
! Subprogram not used     logical          , intent(out)    :: readvar   ! If variable exists or not
! Subprogram not used     logical, optional, intent(in)     :: print_err ! If should print about error
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: ret     ! return value
! Subprogram not used     logical :: log_err ! if should log error
! Subprogram not used     character(len=*),parameter :: subname='check_var' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if ( present(print_err) )then
! Subprogram not used        log_err = print_err
! Subprogram not used     else
! Subprogram not used        log_err = .true.
! Subprogram not used     end if
! Subprogram not used     readvar = .true.
! Subprogram not used     call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
! Subprogram not used     ret = PIO_inq_varid (ncid, varname, vardesc)
! Subprogram not used     if (ret /= PIO_noerr) then
! Subprogram not used        readvar = .false.
! Subprogram not used        if (masterproc .and. log_err) &
! Subprogram not used             write(iulog,*) subname//': variable ',trim(varname),' is not on dataset'
! Subprogram not used     end if
! Subprogram not used     call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
! Subprogram not used 
! Subprogram not used   end subroutine check_var

!-----------------------------------------------------------------------

! Subprogram not used   subroutine check_dim(ncid, dimname, value)
! Subprogram not used 
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! Validity check on dimension
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(in) :: ncid      ! PIO file handle
! Subprogram not used     character(len=*), intent(in) :: dimname   ! Dimension name
! Subprogram not used     integer, intent(in)          :: value     ! Expected dimension size
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: dimid, dimlen    ! temporaries
! Subprogram not used     integer :: status           ! error code      
! Subprogram not used     character(len=*),parameter :: subname='check_dim' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     status = pio_inq_dimid (ncid, trim(dimname), dimid)
! Subprogram not used     status = pio_inq_dimlen (ncid, dimid, dimlen)
! Subprogram not used     if (dimlen /= value) then
! Subprogram not used        write(iulog,*) subname//' ERROR: mismatch of input dimension ',dimlen, &
! Subprogram not used             ' with expected value ',value,' for variable ',trim(dimname)
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine check_dim

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_enddef(ncid)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! enddef netcdf file
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status   ! error status
! Subprogram not used     character(len=*),parameter :: subname='ncd_enddef' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     status = PIO_enddef(ncid)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_enddef

  !-----------------------------------------------------------------------

  subroutine ncd_inqdid(ncid,name,dimid,dimexist)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! inquire on a dimension id
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! netcdf file id
    character(len=*), intent(in) :: name      ! dimension name
    integer         , intent(out):: dimid     ! dimension id
    logical,optional, intent(out):: dimexist  ! if this dimension exists or not
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    if ( present(dimexist) )then
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    end if
    status = PIO_inq_dimid(ncid,name,dimid)
    if ( present(dimexist) )then
       if ( status == PIO_NOERR)then
          dimexist = .true.
       else
          dimexist = .false.
       end if
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    end if

  end subroutine ncd_inqdid

!-----------------------------------------------------------------------

  subroutine ncd_inqdlen(ncid,dimid,len,name)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! enddef netcdf file
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid       ! netcdf file id
    integer          , intent(inout) :: dimid      ! dimension id
    integer          , intent(out)   :: len        ! dimension len
    character(len=*), optional, intent(in) :: name ! dimension name
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    if ( present(name) )then
       call ncd_inqdid(ncid,name,dimid)
    end if
    len = -1
    status = PIO_inq_dimlen(ncid,dimid,len)

  end subroutine ncd_inqdlen

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_inqdname(ncid,dimid,dname)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! inquire dim name
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(in) :: ncid      ! netcdf file id
! Subprogram not used     integer          , intent(in) :: dimid     ! dimension id
! Subprogram not used     character(len=*) , intent(out):: dname     ! dimension name
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     status = PIO_inq_dimname(ncid,dimid,dname)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_inqdname

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     type(file_desc_t), intent(inout):: ncid
! Subprogram not used     logical          , intent(out)  :: isgrid2d
! Subprogram not used     integer          , intent(out)  :: ni
! Subprogram not used     integer          , intent(out)  :: nj
! Subprogram not used     integer          , intent(out)  :: ns
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer  :: dimid                                ! netCDF id
! Subprogram not used     integer  :: ier                                  ! error status 
! Subprogram not used     character(len=32) :: subname = 'surfrd_filedims' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ni = 0
! Subprogram not used     nj = 0
! Subprogram not used 
! Subprogram not used     call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
! Subprogram not used     ier = pio_inq_dimid (ncid, 'lon', dimid)
! Subprogram not used     if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, ni)
! Subprogram not used     ier = pio_inq_dimid (ncid, 'lat', dimid)
! Subprogram not used     if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, nj)
! Subprogram not used 
! Subprogram not used     ier = pio_inq_dimid (ncid, 'lsmlon', dimid)
! Subprogram not used     if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, ni)
! Subprogram not used     ier = pio_inq_dimid (ncid, 'lsmlat', dimid)
! Subprogram not used     if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, nj)
! Subprogram not used 
! Subprogram not used     ier = pio_inq_dimid (ncid, 'ni', dimid)
! Subprogram not used     if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, ni)
! Subprogram not used     ier = pio_inq_dimid (ncid, 'nj', dimid)
! Subprogram not used     if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, nj)
! Subprogram not used 
! Subprogram not used     ier = pio_inq_dimid (ncid, 'gridcell', dimid)
! Subprogram not used     if (ier == PIO_NOERR) then
! Subprogram not used        ier = pio_inq_dimlen(ncid, dimid, ni)
! Subprogram not used        if (ier == PIO_NOERR) nj = 1
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
! Subprogram not used 
! Subprogram not used     if (ni == 0 .or. nj == 0) then
! Subprogram not used        write(iulog,*) trim(subname),' ERROR: ni,nj = ',ni,nj,' cannot be zero '
! Subprogram not used        call shr_sys_abort()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (nj == 1) then
! Subprogram not used        isgrid2d = .false.
! Subprogram not used     else
! Subprogram not used        isgrid2d = .true.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ns = ni*nj
! Subprogram not used 
! Subprogram not used   end subroutine ncd_inqfdims

!-----------------------------------------------------------------------

  subroutine ncd_inqvid(ncid,name,varid,vardesc,readvar)
    
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Inquire on a variable ID
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: name      ! variable name
    integer          , intent(out)   :: varid     ! variable id
    type(Var_desc_t) , intent(out)   :: vardesc   ! variable descriptor
    logical, optional, intent(out)   :: readvar   ! does variable exist
    ! !LOCAL VARIABLES:
    integer :: ret               ! return code
    character(len=*),parameter :: subname='ncd_inqvid' ! subroutine name
    !-----------------------------------------------------------------------

    if (present(readvar)) then
       readvar = .false.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ret = PIO_inq_varid(ncid,name,vardesc)
       if (ret /= PIO_noerr) then
          if (masterproc) write(iulog,*) subname//': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    else
       ret = PIO_inq_varid(ncid,name,vardesc)
    endif
    varid = vardesc%varid
 
  end subroutine ncd_inqvid

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_inqvdims(ncid,ndims,vardesc)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! inquire variable dimensions
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
! Subprogram not used     integer          , intent(out)  :: ndims     ! variable ndims
! Subprogram not used     type(Var_desc_t) , intent(inout):: vardesc   ! variable descriptor
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ndims = -1
! Subprogram not used     status = PIO_inq_varndims(ncid,vardesc,ndims)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_inqvdims

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_inqvname(ncid,varid,vname,vardesc)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! inquire variable name
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
! Subprogram not used     integer          , intent(in)   :: varid     ! variable id
! Subprogram not used     character(len=*) , intent(out)  :: vname     ! variable vname
! Subprogram not used     type(Var_desc_t) , intent(inout):: vardesc   ! variable descriptor
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     vname = ''
! Subprogram not used     status = PIO_inq_varname(ncid,vardesc,vname)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_inqvname

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_inqvdids(ncid,dids,vardesc)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! inquire variable dimension ids
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(in)  :: ncid      ! netcdf file id
! Subprogram not used     integer         ,intent(out)  :: dids(:)   ! variable dids
! Subprogram not used     type(Var_desc_t),intent(inout):: vardesc   ! variable descriptor
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     dids = -1
! Subprogram not used     status = PIO_inq_vardimid(ncid,vardesc,dids)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_inqvdids

!-----------------------------------------------------------------------
! Subprogram not used   subroutine ncd_putatt_int(ncid,varid,attrib,value,xtype)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! put integer attributes
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     integer          ,intent(in)    :: varid     ! netcdf var id
! Subprogram not used     character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
! Subprogram not used     integer          ,intent(in)    :: value     ! netcdf attrib value
! Subprogram not used     integer,optional ,intent(in)    :: xtype     ! netcdf data type
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     status = PIO_put_att(ncid,varid,trim(attrib),value)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_putatt_int

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_putatt_char(ncid,varid,attrib,value,xtype)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! put character attributes
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     integer          ,intent(in)    :: varid     ! netcdf var id
! Subprogram not used     character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
! Subprogram not used     character(len=*) ,intent(in)    :: value     ! netcdf attrib value
! Subprogram not used     integer,optional ,intent(in)    :: xtype     ! netcdf data type
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     status = PIO_put_att(ncid,varid,trim(attrib),value)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_putatt_char

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_putatt_real(ncid,varid,attrib,value,xtype)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! put real attributes
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     integer          ,intent(in)    :: varid     ! netcdf var id
! Subprogram not used     character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
! Subprogram not used     real(r8)         ,intent(in)    :: value     ! netcdf attrib value
! Subprogram not used     integer          ,intent(in)    :: xtype     ! netcdf data type
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     real*4  :: value4
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     value4 = value
! Subprogram not used 
! Subprogram not used     if (xtype == pio_double) then
! Subprogram not used        status = PIO_put_att(ncid,varid,trim(attrib),value)
! Subprogram not used     else
! Subprogram not used        status = PIO_put_att(ncid,varid,trim(attrib),value4)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine ncd_putatt_real

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_defdim(ncid,attrib,value,dimid)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! define dimension
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(in) :: ncid      ! netcdf file id
! Subprogram not used     character(len=*) , intent(in) :: attrib    ! netcdf attrib
! Subprogram not used     integer          , intent(in) :: value     ! netcdf attrib value
! Subprogram not used     integer          , intent(out):: dimid     ! netcdf dimension id
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: status
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     status = pio_def_dim(ncid,attrib,value,dimid)
! Subprogram not used 
! Subprogram not used   end subroutine ncd_defdim

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_defvar_bynf(ncid, varname, xtype, ndims, dimid, varid, &
! Subprogram not used        long_name, units, cell_method, missing_value, fill_value, &
! Subprogram not used        imissing_value, ifill_value, comment, flag_meanings, &
! Subprogram not used        flag_values, nvalid_range )
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     !  Define a netcdf variable
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid                  ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)  :: varname                 ! variable name
! Subprogram not used     integer          , intent(in)  :: xtype                   ! external type
! Subprogram not used     integer          , intent(in)  :: ndims                   ! number of dims
! Subprogram not used     integer          , intent(inout) :: varid                 ! returned var id
! Subprogram not used     integer          , intent(in), optional :: dimid(:)       ! dimids
! Subprogram not used     character(len=*) , intent(in), optional :: long_name      ! attribute
! Subprogram not used     character(len=*) , intent(in), optional :: units          ! attribute
! Subprogram not used     character(len=*) , intent(in), optional :: cell_method    ! attribute
! Subprogram not used     character(len=*) , intent(in), optional :: comment        ! attribute
! Subprogram not used     character(len=*) , intent(in), optional :: flag_meanings(:) ! attribute
! Subprogram not used     real(r8)         , intent(in), optional :: missing_value  ! attribute for real
! Subprogram not used     real(r8)         , intent(in), optional :: fill_value     ! attribute for real
! Subprogram not used     integer          , intent(in), optional :: imissing_value ! attribute for int
! Subprogram not used     integer          , intent(in), optional :: ifill_value    ! attribute for int
! Subprogram not used     integer          , intent(in), optional :: flag_values(:)  ! attribute for int
! Subprogram not used     integer          , intent(in), optional :: nvalid_range(2)  ! attribute for int
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: n                   ! indices
! Subprogram not used     integer :: ldimid(4)           ! local dimid
! Subprogram not used     integer :: dimid0(1)           ! local dimid
! Subprogram not used     integer :: status              ! error status 
! Subprogram not used     integer :: lxtype              ! local external type (in case logical variable)
! Subprogram not used     type(var_desc_t)   :: vardesc  ! local vardesc
! Subprogram not used     character(len=128) :: dimname  ! temporary
! Subprogram not used     character(len=256) :: str      ! temporary
! Subprogram not used     character(len=*),parameter :: subname='ncd_defvar_bynf' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     varid = -1
! Subprogram not used 
! Subprogram not used     dimid0 = 0
! Subprogram not used     ldimid = 0
! Subprogram not used     if (present(dimid)) then
! Subprogram not used        ldimid(1:ndims) = dimid(1:ndims)
! Subprogram not used     else   ! ndims must be zero if dimid not present
! Subprogram not used        if (ndims /= 0) then
! Subprogram not used           write(iulog,*) subname//' ERROR: dimid not supplied and ndims ne 0 ',trim(varname),ndims
! Subprogram not used           call shr_sys_abort()
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if ( xtype == ncd_log )then
! Subprogram not used        lxtype = ncd_int
! Subprogram not used     else
! Subprogram not used        lxtype = xtype
! Subprogram not used     end if
! Subprogram not used     if (masterproc .and. debug > 1) then
! Subprogram not used        write(iulog,*) 'Error in defining variable = ', trim(varname)
! Subprogram not used        write(iulog,*) subname//' ',trim(varname),lxtype,ndims,ldimid(1:ndims)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (ndims >  0) then 
! Subprogram not used        status = pio_inq_dimname(ncid,ldimid(ndims),dimname)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Define variable
! Subprogram not used     if (present(dimid)) then
! Subprogram not used        status = PIO_def_var(ncid,trim(varname),lxtype,dimid(1:ndims),vardesc)
! Subprogram not used     else
! Subprogram not used        status = PIO_def_var(ncid,trim(varname),lxtype,dimid0        ,vardesc)
! Subprogram not used     endif
! Subprogram not used     varid = vardesc%varid
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! Add attributes
! Subprogram not used     !
! Subprogram not used     if (present(long_name)) then
! Subprogram not used        call ncd_putatt(ncid, varid, 'long_name', trim(long_name))
! Subprogram not used     end if
! Subprogram not used     if (present(flag_values)) then
! Subprogram not used        status = PIO_put_att(ncid,varid,'flag_values',flag_values)
! Subprogram not used        if ( .not. present(flag_meanings)) then
! Subprogram not used           write(iulog,*) 'Error in defining variable = ', trim(varname)
! Subprogram not used           call shr_sys_abort( subname//" ERROR:: flag_values set -- but not flag_meanings" )
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if (present(flag_meanings)) then
! Subprogram not used        if ( .not. present(flag_values)) then
! Subprogram not used           write(iulog,*) 'Error in defining variable = ', trim(varname)
! Subprogram not used           call shr_sys_abort( subname//" ERROR:: flag_meanings set -- but not flag_values" )
! Subprogram not used        end if
! Subprogram not used        if ( size(flag_values) /= size(flag_meanings) ) then
! Subprogram not used           write(iulog,*) 'Error in defining variable = ', trim(varname)
! Subprogram not used           call shr_sys_abort( subname//" ERROR:: flag_meanings and flag_values dimension different")
! Subprogram not used        end if
! Subprogram not used        str = flag_meanings(1)
! Subprogram not used        do n = 1, size(flag_meanings)
! Subprogram not used           if ( index(flag_meanings(n), ' ') /= 0 )then
! Subprogram not used              write(iulog,*) 'Error in defining variable = ', trim(varname)
! Subprogram not used              call shr_sys_abort( subname//" ERROR:: flag_meanings has an invalid space in it" )
! Subprogram not used           end if
! Subprogram not used           if ( n > 1 ) str = trim(str)//" "//flag_meanings(n)
! Subprogram not used        end do
! Subprogram not used        status = PIO_put_att(ncid,varid,'flag_meanings', trim(str) )
! Subprogram not used     end if
! Subprogram not used     if (present(comment)) then
! Subprogram not used        call ncd_putatt(ncid, varid, 'comment', trim(comment))
! Subprogram not used     end if
! Subprogram not used     if (present(units)) then
! Subprogram not used        call ncd_putatt(ncid, varid, 'units', trim(units))
! Subprogram not used     end if
! Subprogram not used     if (present(cell_method)) then
! Subprogram not used        str = 'time: ' // trim(cell_method)
! Subprogram not used        call ncd_putatt(ncid, varid, 'cell_methods', trim(str))
! Subprogram not used     end if
! Subprogram not used     if (present(fill_value)) then
! Subprogram not used        call ncd_putatt(ncid, varid, '_FillValue', fill_value, lxtype)
! Subprogram not used     end if
! Subprogram not used     if (present(missing_value)) then
! Subprogram not used        call ncd_putatt(ncid, varid, 'missing_value', missing_value, lxtype)
! Subprogram not used     end if
! Subprogram not used     if (present(ifill_value)) then
! Subprogram not used        call ncd_putatt(ncid, varid, '_FillValue', ifill_value, lxtype)
! Subprogram not used     end if
! Subprogram not used     if (present(imissing_value)) then
! Subprogram not used        call ncd_putatt(ncid, varid, 'missing_value', imissing_value, lxtype)
! Subprogram not used     end if
! Subprogram not used     if (present(nvalid_range)) then
! Subprogram not used        status = PIO_put_att(ncid,varid,'valid_range', nvalid_range )
! Subprogram not used     end if
! Subprogram not used     if ( xtype == ncd_log )then
! Subprogram not used        status = PIO_put_att(ncid,varid,'flag_values',     (/0, 1/) )
! Subprogram not used        status = PIO_put_att(ncid,varid,'flag_meanings',  "FALSE TRUE" )
! Subprogram not used        status = PIO_put_att(ncid,varid,'valid_range',    (/0, 1/) )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine ncd_defvar_bynf

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_defvar_bygrid(ncid, varname, xtype, &
! Subprogram not used        dim1name, dim2name, dim3name, dim4name, dim5name, &
! Subprogram not used        long_name, units, cell_method, missing_value, fill_value, &
! Subprogram not used        imissing_value, ifill_value, comment, &
! Subprogram not used        flag_meanings, flag_values, nvalid_range )
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     !  Define a netcdf variable
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid                 ! netcdf file id
! Subprogram not used     character(len=*), intent(in)  :: varname                 ! variable name
! Subprogram not used     integer         , intent(in)  :: xtype                   ! external type
! Subprogram not used     character(len=*), intent(in), optional :: dim1name       ! dimension name
! Subprogram not used     character(len=*), intent(in), optional :: dim2name       ! dimension name
! Subprogram not used     character(len=*), intent(in), optional :: dim3name       ! dimension name
! Subprogram not used     character(len=*), intent(in), optional :: dim4name       ! dimension name
! Subprogram not used     character(len=*), intent(in), optional :: dim5name       ! dimension name
! Subprogram not used     character(len=*), intent(in), optional :: long_name      ! attribute
! Subprogram not used     character(len=*), intent(in), optional :: units          ! attribute
! Subprogram not used     character(len=*), intent(in), optional :: cell_method    ! attribute
! Subprogram not used     character(len=*), intent(in), optional :: comment        ! attribute
! Subprogram not used     character(len=*), intent(in), optional :: flag_meanings(:) ! attribute
! Subprogram not used     real(r8)        , intent(in), optional :: missing_value  ! attribute for real
! Subprogram not used     real(r8)        , intent(in), optional :: fill_value     ! attribute for real
! Subprogram not used     integer         , intent(in), optional :: imissing_value ! attribute for int
! Subprogram not used     integer         , intent(in), optional :: ifill_value    ! attribute for int
! Subprogram not used     integer         , intent(in), optional :: flag_values(:)  ! attribute for int
! Subprogram not used     integer         , intent(in), optional :: nvalid_range(2)  ! attribute for int
! Subprogram not used     !
! Subprogram not used     ! !REVISION HISTORY:
! Subprogram not used     !
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     !EOP
! Subprogram not used     integer :: n              ! indices
! Subprogram not used     integer :: ndims          ! dimension counter
! Subprogram not used     integer :: dimid(5)       ! dimension ids
! Subprogram not used     integer :: varid          ! variable id
! Subprogram not used     integer :: itmp           ! temporary
! Subprogram not used     character(len=256) :: str ! temporary
! Subprogram not used     character(len=*),parameter :: subname='ncd_defvar_bygrid' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     dimid(:) = 0
! Subprogram not used 
! Subprogram not used     ! Determine dimension ids for variable
! Subprogram not used 
! Subprogram not used     if (present(dim1name)) call ncd_inqdid(ncid, dim1name, dimid(1))
! Subprogram not used     if (present(dim2name)) call ncd_inqdid(ncid, dim2name, dimid(2))
! Subprogram not used     if (present(dim3name)) call ncd_inqdid(ncid, dim3name, dimid(3))
! Subprogram not used     if (present(dim4name)) call ncd_inqdid(ncid, dim4name, dimid(4))
! Subprogram not used     if (present(dim5name)) call ncd_inqdid(ncid, dim5name, dimid(5))
! Subprogram not used 
! Subprogram not used     ! Permute dim1 and dim2 if necessary
! Subprogram not used 
! Subprogram not used     ! Define variable
! Subprogram not used 
! Subprogram not used     ndims = 0
! Subprogram not used     if (present(dim1name)) then
! Subprogram not used        do n = 1, size(dimid)
! Subprogram not used           if (dimid(n) /= 0) ndims = ndims + 1
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call ncd_defvar_bynf(ncid,varname,xtype,ndims,dimid,varid, &
! Subprogram not used          long_name=long_name, units=units, cell_method=cell_method, &
! Subprogram not used          missing_value=missing_value, fill_value=fill_value, &
! Subprogram not used          imissing_value=imissing_value, ifill_value=ifill_value, &
! Subprogram not used          comment=comment, flag_meanings=flag_meanings, &
! Subprogram not used          flag_values=flag_values, nvalid_range=nvalid_range )
! Subprogram not used 
! Subprogram not used   end subroutine ncd_defvar_bygrid

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_log_var0_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global integer variable
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname   ! variable name
! Subprogram not used     logical          , intent(inout) :: data      ! raw data
! Subprogram not used     logical, optional, intent(out)   :: readvar   ! was var read?
! Subprogram not used     integer, optional, intent(in)    :: nt        ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(1), count(1)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     integer :: idata                ! raw integer data
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     integer :: temp(1)              ! temporary
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_log_var0_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, varid, idata)
! Subprogram not used           if (      idata == 0 )then
! Subprogram not used              data = .false.
! Subprogram not used           else if ( idata == 1 )then
! Subprogram not used              data = .true.
! Subprogram not used           else
! Subprogram not used              call shr_sys_abort( subname// &
! Subprogram not used                   ' ERROR: bad integer value for logical data' )
! Subprogram not used           end if
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = nt
! Subprogram not used           count(1) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        if ( data )then
! Subprogram not used           temp(1) = 1
! Subprogram not used        else
! Subprogram not used           temp(1) = 0
! Subprogram not used        end if
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, temp)
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_log_var0_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_int_var0_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global integer variable
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname   ! variable name
! Subprogram not used     integer          , intent(inout) :: data      ! raw data
! Subprogram not used     logical, optional, intent(out)   :: readvar   ! was var read?
! Subprogram not used     integer, optional, intent(in)    :: nt        ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(1), count(1)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     integer :: temp(1)              ! temporary
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_int_var0_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, varid, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = nt
! Subprogram not used           count(1) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        temp(1) = data
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, temp)
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_int_var0_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_real_var0_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global real variable
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname   ! variable name
! Subprogram not used     real(r8)         , intent(inout) :: data      ! raw data
! Subprogram not used     logical, optional, intent(out)   :: readvar   ! was var read?
! Subprogram not used     integer, optional, intent(in)    :: nt        ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(1), count(1)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     real(r8):: temp(1)              ! temporary                
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_real_var0_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, vardesc, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     else if (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = nt
! Subprogram not used           count(1) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        temp(1) = data
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, temp)
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_real_var0_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_int_var1_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global integer array
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname   ! variable name
! Subprogram not used     integer          , intent(inout) :: data(:)   ! raw data
! Subprogram not used     logical, optional, intent(out)   :: readvar   ! was var read?
! Subprogram not used     integer, optional, intent(in)    :: nt        ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(2), count(2)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_int_var1_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, varid, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = size(data)
! Subprogram not used           start(2) = nt
! Subprogram not used           count(2) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = size(data)
! Subprogram not used           start(2) = 1
! Subprogram not used           count(2) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, data)
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_int_var1_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_log_var1_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global integer array
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname   ! variable name
! Subprogram not used     logical          , intent(inout) :: data(:)   ! raw data
! Subprogram not used     logical, optional, intent(out)   :: readvar   ! was var read?
! Subprogram not used     integer, optional, intent(in)    :: nt        ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(2), count(2)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     integer, pointer :: idata(:)    ! Temporary integer data to send to file
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_log_var1_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           allocate( idata(size(data)) ) 
! Subprogram not used           status = pio_get_var(ncid, varid, idata)
! Subprogram not used           data = (idata == 1)
! Subprogram not used           if ( any(idata /= 0 .and. idata /= 1) )then
! Subprogram not used              call shr_sys_abort(subname//'ERROR: read in bad integer value(s) for logical data')
! Subprogram not used           end if
! Subprogram not used           deallocate( idata )
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = size(data)
! Subprogram not used           start(2) = nt
! Subprogram not used           count(2) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           count(1) = size(data)
! Subprogram not used           start(2) = 1
! Subprogram not used           count(2) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        allocate( idata(size(data)) ) 
! Subprogram not used        where( data )
! Subprogram not used           idata = 1
! Subprogram not used        elsewhere
! Subprogram not used           idata = 0
! Subprogram not used        end where
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, idata)
! Subprogram not used        deallocate( idata )
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_log_var1_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_real_var1_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global real array
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid                ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag                ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname             ! variable name
! Subprogram not used     real(r8)         , intent(inout) :: data(:)             ! raw data
! Subprogram not used     logical          , optional, intent(out):: readvar      ! was var read?
! Subprogram not used     integer          , optional, intent(in) :: nt           ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(2), count(2)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_real_var1_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, varid, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = 1
! Subprogram not used           start(2) = nt
! Subprogram not used           count(1) = size(data)
! Subprogram not used           count(2) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           start(2) = 1
! Subprogram not used           count(1) = size(data)
! Subprogram not used           count(2) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, data)
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_real_var1_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_char_var1_nf(varname, data, flag, ncid, readvar, nt )
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global char array
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname          ! variable name
! Subprogram not used     character(len=*) , intent(inout) :: data             ! raw data
! Subprogram not used     logical          , optional, intent(out):: readvar   ! was var read?
! Subprogram not used     integer          , optional, intent(in) :: nt        ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                   ! netCDF variable id
! Subprogram not used     integer :: m                       ! indices
! Subprogram not used     integer :: start(2), count(2)      ! output bounds
! Subprogram not used     integer :: status                  ! error code
! Subprogram not used     logical :: varpresent              ! if true, variable is on tape
! Subprogram not used     character(len=32) :: vname         ! variable error checking
! Subprogram not used     character(len=1)  :: tmpString(128)! temp for manipulating output string
! Subprogram not used     type(var_desc_t)  :: vardesc       ! local vardesc pointer
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_char_var1_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, varid, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           count(1) = len_trim(data)
! Subprogram not used           count(2) = 1
! Subprogram not used           do m = 1,count(1)
! Subprogram not used              tmpString(m:m) = data(m:m)
! Subprogram not used           end do
! Subprogram not used           if ( count(1) > size(tmpString) )then
! Subprogram not used              write(iulog,*) subname//' ERROR: input string size is too large:'//trim(data)
! Subprogram not used           end if
! Subprogram not used           start(1) = 1
! Subprogram not used           start(2) = nt
! Subprogram not used           status = pio_put_var(ncid, varid, start=start, count=count, &
! Subprogram not used                ival=tmpString(1:count(1)) )
! Subprogram not used        else
! Subprogram not used           status = pio_put_var(ncid, varid, data )
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used     endif   ! flag
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_char_var1_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_int_var2_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global integer 2D array
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid                 ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag                 ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname              ! variable name
! Subprogram not used     integer          , intent(inout) :: data(:,:)            ! raw data
! Subprogram not used     logical          , optional, intent(out):: readvar       ! was var read?
! Subprogram not used     integer          , optional, intent(in) :: nt            ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(3), count(3)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     logical :: found                ! if true, found lat/lon dims on file
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_int_var2_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_get_var(ncid, varid, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = 1
! Subprogram not used           start(2) = 1
! Subprogram not used           start(3) = nt
! Subprogram not used           count(1) = size(data, dim=1)
! Subprogram not used           count(2) = size(data, dim=2)
! Subprogram not used           count(3) = 1
! Subprogram not used        else
! Subprogram not used           start(1) = 1
! Subprogram not used           start(2) = 1
! Subprogram not used           start(3) = 1
! Subprogram not used           count(1) = size(data, dim=1)
! Subprogram not used           count(2) = size(data, dim=2)
! Subprogram not used           count(3) = 1
! Subprogram not used        end if
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc)
! Subprogram not used        status = pio_put_var(ncid, varid, start, count, data)
! Subprogram not used 
! Subprogram not used     endif   
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_int_var2_nf

!------------------------------------------------------------------------

  subroutine ncd_io_real_var2_nf(varname, data, flag, ncid, readvar, nt)

    !------------------------------------------------------------------------ 
    ! !DESCRIPTION:
    ! netcdf I/O of global real 2D  array
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
    character(len=*), intent(in)    :: flag                ! 'read' or 'write'
    character(len=*), intent(in)    :: varname             ! variable name
    real(r8)        , intent(inout) :: data(:,:)           ! raw data
    logical         , optional, intent(out):: readvar      ! was var read?
    integer         , optional, intent(in) :: nt           ! time sample index
    ! !LOCAL VARIABLES:
    integer :: varid                ! netCDF variable id
    integer :: start(3), count(3)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    logical :: found                ! if true, found lat/lon dims on file
    character(len=*),parameter :: subname='ncd_io_real_var2_nf'
    !-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_get_var(ncid, varid, data)
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          start(2) = 1
          start(3) = nt
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       else
          start(1) = 1
          start(2) = 1
          start(3) = 1
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   

  end subroutine ncd_io_real_var2_nf

!------------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_char_var2_nf(varname, data, flag, ncid, readvar, nt)
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global character array
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
! Subprogram not used     character(len=*), intent(in)    :: flag                ! 'read' or 'write'
! Subprogram not used     character(len=*), intent(in)    :: varname             ! variable name
! Subprogram not used     character(len=*), intent(inout) :: data(:)             ! raw data
! Subprogram not used     logical         , optional, intent(out):: readvar      ! was var read?
! Subprogram not used     integer         , optional, intent(in) :: nt           ! time sample index
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     integer :: varid                ! netCDF variable id
! Subprogram not used     integer :: start(3), count(3)   ! output bounds
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     logical :: varpresent           ! if true, variable is on tape
! Subprogram not used     character(len=32) :: vname      ! variable error checking
! Subprogram not used     type(var_desc_t)  :: vardesc    ! local vardesc pointer
! Subprogram not used     logical :: found                ! if true, found lat/lon dims on file
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_char_var2_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           data   = ' '
! Subprogram not used           status = pio_get_var(ncid, varid, data)
! Subprogram not used        endif
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid  (ncid, varname, varid, vardesc)
! Subprogram not used        if (present(nt))      then
! Subprogram not used           start(1) = 1
! Subprogram not used           start(2) = 1
! Subprogram not used           start(3) = nt
! Subprogram not used           count(1) = size(data)
! Subprogram not used           count(2) = len(data)
! Subprogram not used           count(3) = 1
! Subprogram not used           status = pio_put_var(ncid, varid, start, count, data)
! Subprogram not used        else
! Subprogram not used           status = pio_put_var(ncid, varid, data)
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used     endif   
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_char_var2_nf

  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: ncd_io_char_varn_strt_nf
  !
  ! !INTERFACE:
! Subprogram not used   subroutine ncd_io_char_varn_strt_nf(vardesc, data, flag, ncid, &
! Subprogram not used        start )
! Subprogram not used     !
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! netcdf I/O of global character array with start indices input
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: ncid             ! netcdf file id
! Subprogram not used     character(len=*), intent(in)    :: flag             ! 'read' or 'write'
! Subprogram not used     type(var_desc_t), intent(in)    :: vardesc          ! local vardesc pointer
! Subprogram not used     character(len=*), intent(inout) :: data             ! raw data for this index
! Subprogram not used     integer         , intent(in)    :: start(:)         ! output bounds
! Subprogram not used     !
! Subprogram not used     ! !REVISION HISTORY:
! Subprogram not used     !
! Subprogram not used     !
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     !EOP
! Subprogram not used     integer :: status               ! error code
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_char_varn_strt_nf'
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        status = pio_get_var(ncid, vardesc, start, data )
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        status = pio_put_var(ncid, vardesc, start, data )
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_char_varn_strt_nf

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_int_var1(varname, data, dim1name, flag, ncid, nt, readvar)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! I/O for 1d integer field
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname          ! variable name
! Subprogram not used     integer          , pointer       :: data(:)          ! local decomposition data
! Subprogram not used     character(len=*) , intent(in)    :: dim1name         ! dimension name
! Subprogram not used     integer          , optional, intent(in) :: nt        ! time sample index
! Subprogram not used     logical          , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     character(len=32) :: dimname    ! temporary
! Subprogram not used     integer           :: n          ! index      
! Subprogram not used     integer           :: iodnum     ! iodesc num in list
! Subprogram not used     integer           :: varid      ! varid
! Subprogram not used     integer           :: ndims      ! ndims for var
! Subprogram not used     integer           :: ndims_iod  ! ndims iodesc for var
! Subprogram not used     integer           :: dims(4)    ! dim sizes       
! Subprogram not used     integer           :: dids(4)    ! dim ids
! Subprogram not used     integer           :: start(3)   ! netcdf start index
! Subprogram not used     integer           :: count(3)   ! netcdf count index
! Subprogram not used     integer           :: status     ! error code  
! Subprogram not used     logical           :: varpresent ! if true, variable is on tape
! Subprogram not used     integer           :: xtype      ! netcdf data type
! Subprogram not used     integer                , pointer  :: compDOF(:)
! Subprogram not used     type(iodesc_plus_type) , pointer  :: iodesc_plus
! Subprogram not used     type(var_desc_t)                  :: vardesc
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_int_var1' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (masterproc .and. debug > 1) then
! Subprogram not used        write(iulog,*) subname//' ',trim(flag),' ',trim(varname),' ',trim(dim1name)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           status = pio_inq_varndims(ncid, vardesc, ndims)
! Subprogram not used           status = pio_inq_vardimid(ncid, vardesc, dids)
! Subprogram not used           status = pio_inq_vartype (ncid, vardesc, xtype)
! Subprogram not used           status = pio_inq_dimname(ncid,dids(ndims),dimname)
! Subprogram not used           if ('time' == trim(dimname)) then
! Subprogram not used              ndims_iod = ndims - 1
! Subprogram not used           else
! Subprogram not used              ndims_iod = ndims
! Subprogram not used           end if
! Subprogram not used           do n = 1,ndims_iod
! Subprogram not used              status = pio_inq_dimlen(ncid,dids(n),dims(n))
! Subprogram not used           enddo
! Subprogram not used           call ncd_getiodesc(ncid, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
! Subprogram not used                xtype, iodnum)
! Subprogram not used           iodesc_plus => iodesc_list(iodnum)
! Subprogram not used           if (present(nt)) then
! Subprogram not used              call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
! Subprogram not used           end if
! Subprogram not used           call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
! Subprogram not used        end if
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname ,varid, vardesc)
! Subprogram not used        status = pio_inq_varndims(ncid, vardesc, ndims)
! Subprogram not used        status = pio_inq_vardimid(ncid, vardesc, dids)
! Subprogram not used        status = pio_inq_vartype (ncid, vardesc, xtype)
! Subprogram not used        status = pio_inq_dimname(ncid,dids(ndims),dimname)
! Subprogram not used        if ('time' == trim(dimname)) then
! Subprogram not used           ndims_iod = ndims - 1
! Subprogram not used        else
! Subprogram not used           ndims_iod = ndims
! Subprogram not used        end if
! Subprogram not used        do n = 1,ndims_iod
! Subprogram not used           status = pio_inq_dimlen(ncid,dids(n),dims(n))
! Subprogram not used        enddo
! Subprogram not used        call ncd_getiodesc(ncid, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
! Subprogram not used             xtype, iodnum)
! Subprogram not used        iodesc_plus => iodesc_list(iodnum)
! Subprogram not used        if (present(nt)) then
! Subprogram not used           call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
! Subprogram not used        end if
! Subprogram not used        call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=0)
! Subprogram not used 
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        if (masterproc) then
! Subprogram not used           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
! Subprogram not used           call shr_sys_abort()
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_int_var1

!-----------------------------------------------------------------------

! Subprogram not used   subroutine ncd_io_log_var1(varname, data, dim1name, &
! Subprogram not used        flag, ncid, nt, readvar)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! !DESCRIPTION:
! Subprogram not used     ! I/O for 1d integer field
! Subprogram not used     !
! Subprogram not used     ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
! Subprogram not used     character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
! Subprogram not used     character(len=*) , intent(in)    :: varname          ! variable name
! Subprogram not used     logical          , pointer       :: data(:)          ! local decomposition data
! Subprogram not used     character(len=*) , intent(in)    :: dim1name         ! dimension name
! Subprogram not used     integer          , optional, intent(in) :: nt        ! time sample index
! Subprogram not used     logical          , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! Subprogram not used     ! !LOCAL VARIABLES:
! Subprogram not used     character(len=32) :: dimname    ! temporary
! Subprogram not used     integer           :: n          ! index      
! Subprogram not used     integer           :: iodnum     ! iodesc num in list
! Subprogram not used     integer           :: varid      ! varid
! Subprogram not used     integer           :: ndims      ! ndims for var
! Subprogram not used     integer           :: ndims_iod  ! ndims iodesc for var
! Subprogram not used     integer           :: dims(4)    ! dim sizes       
! Subprogram not used     integer           :: dids(4)    ! dim ids
! Subprogram not used     integer           :: start(3)   ! netcdf start index
! Subprogram not used     integer           :: count(3)   ! netcdf count index
! Subprogram not used     integer           :: status     ! error code  
! Subprogram not used     integer, pointer  :: idata(:)   ! Temporary integer data to send to file
! Subprogram not used     logical           :: varpresent ! if true, variable is on tape
! Subprogram not used     integer           :: xtype      ! netcdf data type
! Subprogram not used     integer                , pointer  :: compDOF(:)
! Subprogram not used     type(iodesc_plus_type) , pointer  :: iodesc_plus
! Subprogram not used     type(var_desc_t)                  :: vardesc
! Subprogram not used     character(len=*),parameter :: subname='ncd_io_log_var1' ! subroutine name
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (masterproc .and. debug > 1) then
! Subprogram not used        write(iulog,*) subname//' ',trim(flag),' ',trim(varname)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag == 'read') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
! Subprogram not used        if (varpresent) then
! Subprogram not used           allocate( idata(size(data)) ) 
! Subprogram not used           status = pio_inq_varndims(ncid, vardesc, ndims)
! Subprogram not used           status = pio_inq_vardimid(ncid, vardesc, dids)
! Subprogram not used           status = pio_inq_vartype (ncid, vardesc, xtype)
! Subprogram not used           status = pio_inq_dimname(ncid,dids(ndims),dimname)
! Subprogram not used           if ('time' == trim(dimname)) then
! Subprogram not used              ndims_iod = ndims - 1
! Subprogram not used           else
! Subprogram not used              ndims_iod = ndims
! Subprogram not used           end if
! Subprogram not used           do n = 1,ndims_iod
! Subprogram not used              status = pio_inq_dimlen(ncid,dids(n),dims(n))
! Subprogram not used           enddo
! Subprogram not used           call ncd_getiodesc(ncid,  ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
! Subprogram not used                xtype, iodnum)
! Subprogram not used           iodesc_plus => iodesc_list(iodnum)
! Subprogram not used           if (present(nt)) then
! Subprogram not used              call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
! Subprogram not used           end if
! Subprogram not used           call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, idata, status)
! Subprogram not used           data = (idata == 1)
! Subprogram not used           if ( any(idata /= 0 .and. idata /= 1) )then
! Subprogram not used              call shr_sys_abort( subname//' ERROR: read in bad integer value(s) for logical data' )
! Subprogram not used           end if
! Subprogram not used           deallocate( idata )
! Subprogram not used        end if
! Subprogram not used        if (present(readvar)) readvar = varpresent
! Subprogram not used 
! Subprogram not used     elseif (flag == 'write') then
! Subprogram not used 
! Subprogram not used        call ncd_inqvid(ncid, varname ,varid, vardesc)
! Subprogram not used        status = pio_inq_varndims(ncid, vardesc, ndims)
! Subprogram not used        status = pio_inq_vardimid(ncid, vardesc, dids)
! Subprogram not used        status = pio_inq_vartype (ncid, vardesc, xtype)
! Subprogram not used        status = pio_inq_dimname(ncid,dids(ndims),dimname)
! Subprogram not used        if ('time' == trim(dimname)) then
! Subprogram not used           ndims_iod = ndims - 1
! Subprogram not used        else
! Subprogram not used           ndims_iod = ndims
! Subprogram not used        end if
! Subprogram not used        do n = 1,ndims_iod
! Subprogram not used           status = pio_inq_dimlen(ncid,dids(n),dims(n))
! Subprogram not used        enddo
! Subprogram not used        call ncd_getiodesc(ncid,  ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
! Subprogram not used             xtype, iodnum)
! Subprogram not used        iodesc_plus => iodesc_list(iodnum)
! Subprogram not used        if (present(nt)) then
! Subprogram not used           call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
! Subprogram not used        end if
! Subprogram not used        allocate( idata(size(data)) ) 
! Subprogram not used        where( data )
! Subprogram not used           idata = 1
! Subprogram not used        elsewhere
! Subprogram not used           idata = 0
! Subprogram not used        end where
! Subprogram not used        call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, idata, status, fillval=0)
! Subprogram not used        deallocate( idata )
! Subprogram not used 
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        if (masterproc) then
! Subprogram not used           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
! Subprogram not used           call shr_sys_abort()
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine ncd_io_log_var1

!-----------------------------------------------------------------------

  subroutine ncd_io_real_var1(varname, data, dim1name, &
                              flag, ncid, nt, readvar)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! I/O for 1d real field
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid                 ! netcdf file id
    character(len=*), intent(in)  :: flag                   ! 'read' or 'write'
    character(len=*), intent(in)  :: varname                ! variable name
    real(r8)        , pointer     :: data(:)                ! local decomposition data
    character(len=*), intent(in)  :: dim1name               ! dimension name
    integer         , optional, intent(in) :: nt            ! time sample index
    logical         , optional, intent(out):: readvar       ! true => variable is on initial dataset (read only)
    ! !LOCAL VARIABLES:
    character(len=32) :: dimname    ! temporary
    integer           :: iodnum     ! iodesc num in list
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: n          ! index      
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    integer           :: status     ! error code  
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    integer                , pointer  :: compDOF(:)
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var1' ! subroutine name
    !-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname)
    endif

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_inq_varndims(ncid, vardesc, ndims)
          status = pio_inq_vardimid(ncid,vardesc, dids)
          status = pio_inq_vartype(ncid, vardesc, xtype)
          status = pio_inq_dimname(ncid,dids(ndims),dimname)
          if ('time' == trim(dimname)) then
             ndims_iod = ndims - 1
          else
             ndims_iod = ndims
          end if
          do n = 1,ndims_iod
             status = pio_inq_dimlen(ncid,dids(n),dims(n))
          enddo
          call ncd_getiodesc(ncid,  ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
               xtype, iodnum)
          iodesc_plus => iodesc_list(iodnum)
          if (present(nt)) then
             call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
          end if
          call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc, dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid,  ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
            xtype, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call shr_sys_abort()
       endif

    endif

  end subroutine ncd_io_real_var1

!------------------------------------------------------------------------

  subroutine ncd_getiodesc(ncid, ndims, dims, dimids, xtype, iodnum)

    !------------------------------------------------------------------------
    ! !DESCRIPTION: 
    ! Returns an index to an io descriptor
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid       ! PIO file descriptor
    integer          , intent(in)    :: ndims      ! ndims for var      
    integer          , intent(in)    :: dims(:)    ! dim sizes
    integer          , intent(in)    :: dimids(:)  ! dim ids
    integer          , intent(in)    :: xtype      ! file external type
    integer          , intent(out)   :: iodnum     ! iodesc num in list
    ! !LOCAL VARIABLES:
    integer :: k,m,n,cnt                     ! indices
    integer :: basetype                      ! pio basetype
    integer :: lsize                         ! local size 
    integer :: gsize                         ! global size
    integer :: status                        ! error status
    logical :: found                         ! true => found created iodescriptor
    integer :: ndims_file                    ! temporary
    character(len=64) dimname_file           ! dimension name on file
    character(len=64) dimname_iodesc         ! dimension name from io descriptor
    integer, pointer  :: compDOF(:)
    character(len=32) :: subname = 'ncd_getiodesc'
    !------------------------------------------------------------------------

    ! Determining if need to create a new io descriptor

    n = 1
    found = .false.
    do while (n <= num_iodesc .and. .not.found)
       if (ndims == iodesc_list(n)%ndims .and. xtype == iodesc_list(n)%type) then
          found = .true.
          ! First found implies that dimension sizes are the same 
          do m = 1,ndims
             if (dims(m) /= iodesc_list(n)%dims(m)) then
                found = .false.
             endif
          enddo
          ! If found - then also check that dimension names are equal - 
          ! dimension ids in iodescriptor are only used to query dimension
          ! names associated with that iodescriptor
          if (found) then
             do m = 1,ndims
                status = PIO_inq_dimname(ncid,dimids(m),dimname_file)
                status = PIO_inquire(ncid, ndimensions=ndims_file)
                if (iodesc_list(n)%dimids(m) > ndims_file) then 
                   found = .false.
                   exit
                else
                   status = PIO_inq_dimname(ncid,iodesc_list(n)%dimids(m),dimname_iodesc)
                   if (trim(dimname_file) .ne. trim(dimname_iodesc)) then
                      found = .false.
                      exit
                   end if
                end if
             end do
          end if
          if (found) then
             iodnum = n
             if (iodnum > num_iodesc) then
                write(iulog,*) trim(subname),' ERROR: iodnum out of range ',iodnum,num_iodesc
                call shr_sys_abort()
             endif
             RETURN
          endif
       endif
       n = n + 1
    enddo

    ! Creating a new io descriptor

    if (ndims > 0) then 
       num_iodesc = num_iodesc + 1
       if (num_iodesc > max_iodesc) then
          write(iulog,*) trim(subname),' ERROR num_iodesc gt max_iodesc ',max_iodesc
          call shr_sys_abort()
       endif
       iodnum = num_iodesc
       if (masterproc .and. debug > 1) then
          write(iulog,*) trim(subname),' creating iodesc at iodnum,ndims,dims(1:ndims),xtype',&
               iodnum,ndims,dims(1:ndims),xtype
       endif
    end if

    ! Initialize the decomposition for PIO

    if (xtype == pio_double ) then
       basetype = PIO_DOUBLE
    else if (xtype == pio_real) then
       basetype  = PIO_DOUBLE
    else if (xtype == pio_int) then
       basetype = PIO_INT
    end if

    gsize = runoff%numr
    lsize = runoff%lnumr

    allocate(compDOF(lsize))

    cnt = 0
    do m = runoff%begr, runoff%endr
       cnt = cnt + 1
       compDOF(cnt) = runoff%gindex(m)
    enddo

    if (debug > 1) then
       do m = 0,npes-1
          if (iam == m) then
             write(iulog,*) trim(subname),' sizes1  = ',iam,gsize,lsize,npes
             write(iulog,*) trim(subname),' compDOF = ',iam,size(compDOF),minval(compDOF),maxval(compDOF)
             call shr_sys_flush(iulog)
          endif
          call mpi_barrier(mpicom_rof,status)
       enddo
    endif

    call pio_initdecomp(pio_subsystem, baseTYPE, dims(1:ndims), compDOF, iodesc_list(iodnum)%iodesc)

    deallocate(compDOF)

    iodesc_list(iodnum)%type  = xtype
    iodesc_list(iodnum)%ndims = ndims
    iodesc_list(iodnum)%dims  = 0
    iodesc_list(iodnum)%dims(1:ndims)   = dims(1:ndims)
    iodesc_list(iodnum)%dimids(1:ndims) = dimids(1:ndims)

  end subroutine ncd_getiodesc

  subroutine ncd_finalize()
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! clean up  PIO
  !
  ! !USES:
  use alloc_mod , only : dealloc_check
  ! !ARGUMENTS:
  implicit none
  ! !LOCAL VARIABLES:
  integer :: n ! index into num_iodesc

  do n = 1,num_iodesc

     if(associated(iodesc_list(n)%iodesc%start)) then
        call dealloc_check(iodesc_list(n)%iodesc%start,'RtmIO::iodesc%start')
        nullify(iodesc_list(n)%iodesc%start)
     end if

     if(associated(iodesc_list(n)%iodesc%count)) then
        call dealloc_check(iodesc_list(n)%iodesc%count,'RtmIO::iodesc%count')
        nullify(iodesc_list(n)%iodesc%count)
     end if

     if(associated(iodesc_list(n)%iodesc%dest_ioproc)) then
        call dealloc_check(iodesc_list(n)%iodesc%dest_ioproc,'RtmIO::iodesc%dest_ioproc')
        nullify(iodesc_list(n)%iodesc%dest_ioproc)
     end if

     if(associated(iodesc_list(n)%iodesc%dest_ioindex)) then
        call dealloc_check(iodesc_list(n)%iodesc%dest_ioindex,'RtmIO::iodesc%ioindex')
        nullify(iodesc_list(n)%iodesc%dest_ioindex)
     end if

  end do
     
  end subroutine ncd_finalize

end module RtmIO
