!>
!! @file
!! @brief NetCDF interface routines
!!
!! $Revision: 894 $
!! $LastChangedDate: 2013-12-13 15:04:58 -0700 (Fri, 13 Dec 2013) $
!!
!<
module nf_mod

  use alloc_mod

  use pio_kinds, only: i4,r4,r8,pio_offset
  use pio_types, only: file_desc_t, iosystem_desc_t, var_desc_t, pio_noerr, pio_iotype_netcdf, &
	pio_iotype_pnetcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c, pio_max_name

  use pio_support, only : Debug, DebugIO, DebugAsync, piodie   
  use pio_utils, only : bad_iotype, check_netcdf

  use netcdf            ! _EXTERNAL
  use pio_support, only : CheckMPIReturn
  use pio_msg_mod
  implicit none
  private
  include 'mpif.h' ! _EXTERNAL
!
! Copyright (C) 2003, Northwestern University and Argonne National Laboratory
! See COPYRIGHT notice in top-level directory.
!
! src/libf/pnetcdf.inc.  Generated from pnetcdf.inc.in by configure.

!
! pnetcdf fortran defines
!

!
! PnetCDF library version numbers
!
      integer PNETCDF_VERSION_MAJOR
      integer PNETCDF_VERSION_MINOR
      integer PNETCDF_VERSION_SUB

      parameter (PNETCDF_VERSION_MAJOR = 1)
      parameter (PNETCDF_VERSION_MINOR = 8)
      parameter (PNETCDF_VERSION_SUB   = 0)

!
! list of PnetCDF options enabled/disabled at configure time
!
      integer PNETCDF_ERANGE_FILL
      integer PNETCDF_SUBFILING
      integer PNETCDF_RELAX_COORD_BOUND
      integer PNETCDF_DEBUG_MODE

      parameter (PNETCDF_ERANGE_FILL       = 1)
      parameter (PNETCDF_SUBFILING         = 0)
      parameter (PNETCDF_RELAX_COORD_BOUND = 0)
      parameter (PNETCDF_DEBUG_MODE        = 0)

!
! external netcdf data types: (must conform with netCDF release)
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
      integer*8         nf_fill_uint
      integer*8         nf_fill_int64
      ! integer*8         nf_fill_uint64    ! no unsigned int*8 in Fortran
      doubleprecision   nf_fill_uint64

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690e+36)
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)


      parameter (nf_fill_uint = 4294967295_8)
      parameter (nf_fill_int64 = -9223372036854775806_8)
      ! parameter (nf_fill_uint64 = 18446744073709551614_8)  ! currently not supported
      parameter (nf_fill_uint64 = 1.8446744073709551614e+19)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_32bit
      integer nf_64bit_data
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_data
      integer nf_format_64bit_offset
      integer nf_format_cdf2
      integer nf_format_cdf5

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_32bit = 16777216)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_cdf2 = 2)
      parameter (nf_format_cdf5 = 5)
      parameter (nf_format_64bit = nf_format_cdf2)
      parameter (nf_format_64bit_offset = nf_format_cdf2)
      parameter (nf_format_64bit_data = nf_format_cdf5)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

      integer*8 nfmpi_unlimited
      parameter (nfmpi_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 512)
      parameter (nf_max_attrs = 4092)
      parameter (nf_max_vars = 4096)
      parameter (nf_max_name = 128)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes: (conform with netCDF release)
!
      integer NF_NOERR
      integer NF2_ERR
      integer NF_EBADID
      integer NF_ENFILE
      integer NF_EEXIST
      integer NF_EINVAL
      integer NF_EPERM
      integer NF_ENOTINDEFINE
      integer NF_EINDEFINE
      integer NF_EINVALCOORDS
      integer NF_EMAXDIMS
      integer NF_ENAMEINUSE
      integer NF_ENOTATT
      integer NF_EMAXATTS
      integer NF_EBADTYPE
      integer NF_EBADDIM
      integer NF_EUNLIMPOS
      integer NF_EMAXVARS
      integer NF_ENOTVAR
      integer NF_EGLOBAL
      integer NF_ENOTNC
      integer NF_ESTS
      integer NF_EMAXNAME
      integer NF_EUNLIMIT
      integer NF_ENORECVARS
      integer NF_ECHAR
      integer NF_EEDGE
      integer NF_ESTRIDE
      integer NF_EBADNAME
      integer NF_ERANGE
      integer NF_ENOMEM
      integer NF_EVARSIZE
      integer NF_EDIMSIZE
      integer NF_ETRUNC
      integer NF_EAXISTYPE
      integer NF_EDAP
      integer NF_ECURL
      integer NF_EIO
      integer NF_ENODATA
      integer NF_EDAPSVC
      integer NF_EDAS
      integer NF_EDDS
      integer NF_EDATADDS
      integer NF_EDAPURL
      integer NF_EDAPCONSTRAINT
      integer NF_ETRANSLATION
      integer NF_EACCESS
      integer NF_EAUTH
      integer NF_ENOTFOUND
      integer NF_ECANTREMOVE

      PARAMETER (NF_NOERR        = 0)     ! No Error
      PARAMETER (NF2_ERR         = -1)    ! Returned for all errors in the v2 API
      PARAMETER (NF_EBADID       = -33)   ! Not a netcdf id
      PARAMETER (NF_ENFILE       = -34)   ! Too many netcdfs open
      PARAMETER (NF_EEXIST       = -35)   ! netcdf file exists and NF_NOCLOBBER
      PARAMETER (NF_EINVAL       = -36)   ! Invalid Argument
      PARAMETER (NF_EPERM        = -37)   ! Write to read only
      PARAMETER (NF_ENOTINDEFINE = -38)   ! Operation not allowed in data mode
      PARAMETER (NF_EINDEFINE    = -39)   ! Operation not allowed in define mode
      PARAMETER (NF_EINVALCOORDS = -40)   ! Index exceeds dimension bound
      PARAMETER (NF_EMAXDIMS     = -41)   ! NF_MAX_DIMS exceeded
      PARAMETER (NF_ENAMEINUSE   = -42)   ! String match to name in use
      PARAMETER (NF_ENOTATT      = -43)   ! Attribute not found
      PARAMETER (NF_EMAXATTS     = -44)   ! NF_MAX_ATTRS exceeded
      PARAMETER (NF_EBADTYPE     = -45)   ! Not a netcdf data type
      PARAMETER (NF_EBADDIM      = -46)   ! Invalid dimension id or name
      PARAMETER (NF_EUNLIMPOS    = -47)   ! NFMPI_UNLIMITED in the wrong index
      PARAMETER (NF_EMAXVARS     = -48)   ! NF_MAX_VARS exceeded
      PARAMETER (NF_ENOTVAR      = -49)   ! Variable not found
      PARAMETER (NF_EGLOBAL      = -50)   ! Action prohibited on NF_GLOBAL varid
      PARAMETER (NF_ENOTNC       = -51)   ! Not a netcdf file
      PARAMETER (NF_ESTS         = -52)   ! In Fortran, string too short
      PARAMETER (NF_EMAXNAME     = -53)   ! NF_MAX_NAME exceeded
      PARAMETER (NF_EUNLIMIT     = -54)   ! NFMPI_UNLIMITED size already in use
      PARAMETER (NF_ENORECVARS   = -55)   ! nc_rec op when there are no record vars
      PARAMETER (NF_ECHAR        = -56)   ! Attempt to convert between text & numbers
      PARAMETER (NF_EEDGE        = -57)   ! Edge+start exceeds dimension bound
      PARAMETER (NF_ESTRIDE      = -58)   ! Illegal stride
      PARAMETER (NF_EBADNAME     = -59)   ! Attribute or variable name contains illegal characters
      PARAMETER (NF_ERANGE       = -60)   ! Math result not representable
      PARAMETER (NF_ENOMEM       = -61)   ! Memory allocation (malloc) failure
      PARAMETER (NF_EVARSIZE     = -62)   ! One or more variable sizes violate format constraints
      PARAMETER (NF_EDIMSIZE     = -63)   ! Invalid dimension size
      PARAMETER (NF_ETRUNC       = -64)   ! File likely truncated or possibly corrupted
      PARAMETER (NF_EAXISTYPE    = -65)   ! Unknown axis type

! Following errors are added for DAP
      PARAMETER (NF_EDAP         = -66)   ! Generic DAP error
      PARAMETER (NF_ECURL        = -67)   ! Generic libcurl error
      PARAMETER (NF_EIO          = -68)   ! Generic IO error
      PARAMETER (NF_ENODATA      = -69)   ! Attempt to access variable with no data
      PARAMETER (NF_EDAPSVC      = -70)   ! DAP server error
      PARAMETER (NF_EDAS         = -71)   ! Malformed or inaccessible DAS
      PARAMETER (NF_EDDS         = -72)   ! Malformed or inaccessible DDS
      PARAMETER (NF_EDATADDS     = -73)   ! Malformed or inaccessible DATADDS
      PARAMETER (NF_EDAPURL      = -74)   ! Malformed DAP URL
      PARAMETER (NF_EDAPCONSTRAINT = -75) ! Malformed DAP Constraint
      PARAMETER (NF_ETRANSLATION = -76)   ! Untranslatable construct
      PARAMETER (NF_EACCESS      = -77)   ! Access Failure
      PARAMETER (NF_EAUTH        = -78)   ! Authorization Failure

! Misc. additional errors
      PARAMETER (NF_ENOTFOUND    = -90)   ! No such file
      PARAMETER (NF_ECANTREMOVE  = -91)   ! Can't remove file

!
! netCDF-4 error codes (copied from netCDF release)
!
      integer NF_EHDFERR
      integer NF_ECANTREAD
      integer NF_ECANTWRITE
      integer NF_ECANTCREATE
      integer NF_EFILEMETA
      integer NF_EDIMMETA
      integer NF_EATTMETA
      integer NF_EVARMETA
      integer NF_ENOCOMPOUND
      integer NF_EATTEXISTS
      integer NF_ENOTNC4
      integer NF_ESTRICTNC3
      integer NF_ENOTNC3
      integer NF_ENOPAR
      integer NF_EPARINIT
      integer NF_EBADGRPID
      integer NF_EBADTYPID
      integer NF_ETYPDEFINED
      integer NF_EBADFIELD
      integer NF_EBADCLASS
      integer NF_EMAPTYPE
      integer NF_ELATEFILL
      integer NF_ELATEDEF
      integer NF_EDIMSCALE
      integer NF_ENOGRP
      integer NF_ESTORAGE
      integer NF_EBADCHUNK
      integer NF_ENOTBUILT
      integer NF_EDISKLESS
      integer NF_ECANTEXTEND
      integer NF_EMPI

      PARAMETER (NF_EHDFERR      = -101)  ! Error at HDF5 layer. 
      PARAMETER (NF_ECANTREAD    = -102)  ! Can't read. 
      PARAMETER (NF_ECANTWRITE   = -103)  ! Can't write. 
      PARAMETER (NF_ECANTCREATE  = -104)  ! Can't create. 
      PARAMETER (NF_EFILEMETA    = -105)  ! Problem with file metadata. 
      PARAMETER (NF_EDIMMETA     = -106)  ! Problem with dimension metadata. 
      PARAMETER (NF_EATTMETA     = -107)  ! Problem with attribute metadata. 
      PARAMETER (NF_EVARMETA     = -108)  ! Problem with variable metadata. 
      PARAMETER (NF_ENOCOMPOUND  = -109)  ! Not a compound type. 
      PARAMETER (NF_EATTEXISTS   = -110)  ! Attribute already exists. 
      PARAMETER (NF_ENOTNC4      = -111)  ! Attempting netcdf-4 operation on netcdf-3 file.   
      PARAMETER (NF_ESTRICTNC3   = -112)  ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      PARAMETER (NF_ENOTNC3      = -113)  ! Attempting netcdf-3 operation on netcdf-4 file.   
      PARAMETER (NF_ENOPAR       = -114)  ! Parallel operation on file opened for non-parallel access.   
      PARAMETER (NF_EPARINIT     = -115)  ! Error initializing for parallel access.   
      PARAMETER (NF_EBADGRPID    = -116)  ! Bad group ID.   
      PARAMETER (NF_EBADTYPID    = -117)  ! Bad type ID.   
      PARAMETER (NF_ETYPDEFINED  = -118)  ! Type has already been defined and may not be edited. 
      PARAMETER (NF_EBADFIELD    = -119)  ! Bad field ID.   
      PARAMETER (NF_EBADCLASS    = -120)  ! Bad class.   
      PARAMETER (NF_EMAPTYPE     = -121)  ! Mapped access for atomic types only.   
      PARAMETER (NF_ELATEFILL    = -122)  ! Attempt to define fill value when data already exists. 
      PARAMETER (NF_ELATEDEF     = -123)  ! Attempt to define var properties, like deflate, after enddef.
      PARAMETER (NF_EDIMSCALE    = -124)  ! Problem with HDF5 dimscales.
      PARAMETER (NF_ENOGRP       = -125)  ! No group found.
      PARAMETER (NF_ESTORAGE     = -126)  ! Can't specify both contiguous and chunking.
      PARAMETER (NF_EBADCHUNK    = -127)  ! Bad chunksize.
      PARAMETER (NF_ENOTBUILT    = -128)  ! Attempt to use feature that was not turned on when netCDF was built.
      PARAMETER (NF_EDISKLESS    = -129)  ! Error in using diskless  access.
      PARAMETER (NF_ECANTEXTEND  = -130)  ! Attempt to extend dataset during ind. I/O operation.
      PARAMETER (NF_EMPI         = -131)  ! MPI operation failed.

!
! PnetCDF error codes start here
!
      integer NF_ESMALL
      integer NF_ENOTINDEP
      integer NF_EINDEP
      integer NF_EFILE
      integer NF_EREAD
      integer NF_EWRITE
      integer NF_EOFILE
      integer NF_EMULTITYPES
      integer NF_EIOMISMATCH
      integer NF_ENEGATIVECNT
      integer NF_EUNSPTETYPE
      integer NF_EINVAL_REQUEST
      integer NF_EAINT_TOO_SMALL
      integer NF_ENOTSUPPORT
      integer NF_ENULLBUF
      integer NF_EPREVATTACHBUF
      integer NF_ENULLABUF
      integer NF_EPENDINGBPUT
      integer NF_EINSUFFBUF
      integer NF_ENOENT
      integer NF_EINTOVERFLOW
      integer NF_ENOTENABLED
      integer NF_EBAD_FILE
      integer NF_ENO_SPACE
      integer NF_EQUOTA
      integer NF_ENULLSTART
      integer NF_ENULLCOUNT
      integer NF_EINVAL_CMODE
      integer NF_ETYPESIZE
      integer NF_ETYPE_MISMATCH
      integer NF_ETYPESIZE_MISMATCH
      integer NF_ESTRICTCDF2
      integer NF_ENOTRECVAR
      integer NF_ENOTFILL
      integer NF_EINVAL_OMODE
      integer NF_EPENDING

      integer NF_EMULTIDEFINE
      integer NF_EMULTIDEFINE_OMODE,      NF_ECMODE
      integer NF_EMULTIDEFINE_DIM_NUM,    NF_EDIMS_NELEMS_MULTIDEFINE
      integer NF_EMULTIDEFINE_DIM_SIZE,   NF_EDIMS_SIZE_MULTIDEFINE
      integer NF_EMULTIDEFINE_DIM_NAME,   NF_EDIMS_NAME_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_NUM,    NF_EVARS_NELEMS_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_NAME,   NF_EVARS_NAME_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_NDIMS,  NF_EVARS_NDIMS_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_DIMIDS, NF_EVARS_DIMIDS_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_TYPE,   NF_EVARS_TYPE_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_LEN,    NF_EVARS_LEN_MULTIDEFINE
      integer NF_EMULTIDEFINE_VAR_BEGIN,  NF_EVARS_BEGIN_MULTIDEFINE
      integer NF_EMULTIDEFINE_NUMRECS,    NF_ENUMRECS_MULTIDEFINE
      integer NF_EMULTIDEFINE_ATTR_NUM
      integer NF_EMULTIDEFINE_ATTR_SIZE
      integer NF_EMULTIDEFINE_ATTR_NAME
      integer NF_EMULTIDEFINE_ATTR_TYPE
      integer NF_EMULTIDEFINE_ATTR_LEN
      integer NF_EMULTIDEFINE_ATTR_VAL
      integer NF_EMULTIDEFINE_FNC_ARGS
      integer NF_EMULTIDEFINE_FILL_MODE
      integer NF_EMULTIDEFINE_VAR_FILL_MODE
      integer NF_EMULTIDEFINE_VAR_FILL_VALUE
      integer NF_EMULTIDEFINE_CMODE

!
! PnetCDF error codes start from -201
!
      PARAMETER (NF_ESMALL                  = -201)   ! size of off_t too small for format
      PARAMETER (NF_ENOTINDEP               = -202)   ! Operation not allowed in collective data mode
      PARAMETER (NF_EINDEP                  = -203)   ! Operation not allowed in independent data mode
      PARAMETER (NF_EFILE                   = -204)   ! Unknown error in file operation
      PARAMETER (NF_EREAD                   = -205)   ! Unknown error in reading file
      PARAMETER (NF_EWRITE                  = -206)   ! Unknown error in writing to file
      PARAMETER (NF_EOFILE                  = -207)   ! file open/creation failed
      PARAMETER (NF_EMULTITYPES             = -208)   ! Multiple types used in memory data
      PARAMETER (NF_EIOMISMATCH             = -209)   ! Input/Output data amount mismatch
      PARAMETER (NF_ENEGATIVECNT            = -210)   ! Negative count is specified
      PARAMETER (NF_EUNSPTETYPE             = -211)   ! Unsupported etype in memory MPI datatype
      PARAMETER (NF_EINVAL_REQUEST          = -212)   ! invalid nonblocking request ID
      PARAMETER (NF_EAINT_TOO_SMALL         = -213)   ! MPI_Aint not large enough to hold requested value
      PARAMETER (NF_ENOTSUPPORT             = -214)   ! feature is not yet supported
      PARAMETER (NF_ENULLBUF                = -215)   ! trying to attach a NULL buffer
      PARAMETER (NF_EPREVATTACHBUF          = -216)   ! previous attached buffer is found
      PARAMETER (NF_ENULLABUF               = -217)   ! no attached buffer is found
      PARAMETER (NF_EPENDINGBPUT            = -218)   ! pending bput is found, cannot detach buffer
      PARAMETER (NF_EINSUFFBUF              = -219)   ! attached buffer is too small
      PARAMETER (NF_ENOENT                  = -220)   ! File does not exist when calling nfmpi_open()
      PARAMETER (NF_EINTOVERFLOW            = -221)   ! Overflow when type cast to 4-byte integer
      PARAMETER (NF_ENOTENABLED             = -222)   ! feature is not enabled
      PARAMETER (NF_EBAD_FILE               = -223)   ! Invalid file name (e.g., path name too long)
      PARAMETER (NF_ENO_SPACE               = -224)   ! Not enough space
      PARAMETER (NF_EQUOTA                  = -225)   ! Quota exceeded
      PARAMETER (NF_ENULLSTART              = -226)   ! argument start is a NULL pointer
      PARAMETER (NF_ENULLCOUNT              = -227)   ! argument count is a NULL pointer
      PARAMETER (NF_EINVAL_CMODE            = -228)   ! Invalid file create mode
      PARAMETER (NF_ETYPESIZE               = -229)   ! MPI derived data type size error (bigger than the variable size)
      PARAMETER (NF_ETYPE_MISMATCH          = -230)   ! element type of the MPI derived data type mismatches the variable type
      PARAMETER (NF_ETYPESIZE_MISMATCH      = -231)   ! file type size mismatches buffer type size
      PARAMETER (NF_ESTRICTCDF2             = -232)   ! Attempting CDF-5 operation on CDF-2 file
      PARAMETER (NF_ENOTRECVAR              = -233)   ! Attempting operation only for record variables
      PARAMETER (NF_ENOTFILL                = -234)   ! Attempting to fill a variable when its fill mode is off
      PARAMETER (NF_EINVAL_OMODE            = -235)   ! Invalid file open mode
      PARAMETER (NF_EPENDING                = -236)   ! Pending nonblocking request is found at file close

!
! PnetCDF header inconsistency errors start from -250
!
      PARAMETER (NF_EMULTIDEFINE            = -250)   ! NC definitions on multiprocesses conflict
      PARAMETER (NF_EMULTIDEFINE_OMODE      = -251)   ! inconsistent file open modes
      PARAMETER (NF_EMULTIDEFINE_DIM_NUM    = -252)   ! inconsistent number of dimensions
      PARAMETER (NF_EMULTIDEFINE_DIM_SIZE   = -253)   ! inconsistent size of dimension
      PARAMETER (NF_EMULTIDEFINE_DIM_NAME   = -254)   ! inconsistent dimension names
      PARAMETER (NF_EMULTIDEFINE_VAR_NUM    = -255)   ! inconsistent number of variables
      PARAMETER (NF_EMULTIDEFINE_VAR_NAME   = -256)   ! inconsistent variable name
      PARAMETER (NF_EMULTIDEFINE_VAR_NDIMS  = -257)   ! inconsistent variable's number of dimensions
      PARAMETER (NF_EMULTIDEFINE_VAR_DIMIDS = -258)   ! inconsistent variable's dimid
      PARAMETER (NF_EMULTIDEFINE_VAR_TYPE   = -259)   ! inconsistent variable's data type
      PARAMETER (NF_EMULTIDEFINE_VAR_LEN    = -260)   ! inconsistent variable's size
      PARAMETER (NF_EMULTIDEFINE_NUMRECS    = -261)   ! inconsistent number of records
      PARAMETER (NF_EMULTIDEFINE_VAR_BEGIN  = -262)   ! inconsistent variable file begin offset (internal use)
      PARAMETER (NF_EMULTIDEFINE_ATTR_NUM   = -263)   ! inconsistent number of attributes
      PARAMETER (NF_EMULTIDEFINE_ATTR_SIZE  = -264)   ! inconsistent memory space used by attribute (internal use)
      PARAMETER (NF_EMULTIDEFINE_ATTR_NAME  = -265)   ! inconsistent attribute name
      PARAMETER (NF_EMULTIDEFINE_ATTR_TYPE  = -266)   ! inconsistent attribute type
      PARAMETER (NF_EMULTIDEFINE_ATTR_LEN   = -267)   ! inconsistent attribute length
      PARAMETER (NF_EMULTIDEFINE_ATTR_VAL   = -268)   ! inconsistent attribute value
      PARAMETER (NF_EMULTIDEFINE_FNC_ARGS   = -269)   ! inconsistent function arguments used in collective API
      PARAMETER (NF_EMULTIDEFINE_FILL_MODE  = -270)   !  inconsistent dataset fill mode
      PARAMETER (NF_EMULTIDEFINE_VAR_FILL_MODE  = -271) ! inconsistent variable fill mode
      PARAMETER (NF_EMULTIDEFINE_VAR_FILL_VALUE = -272) ! inconsistent variable fill value
      PARAMETER (NF_EMULTIDEFINE_CMODE      = -273)   ! inconsistent file create modes

      PARAMETER(NF_ECMODE                  =NF_EMULTIDEFINE_OMODE)
      PARAMETER(NF_EDIMS_NELEMS_MULTIDEFINE=NF_EMULTIDEFINE_DIM_NUM)
      PARAMETER(NF_EDIMS_SIZE_MULTIDEFINE  =NF_EMULTIDEFINE_DIM_SIZE)
      PARAMETER(NF_EDIMS_NAME_MULTIDEFINE  =NF_EMULTIDEFINE_DIM_NAME)
      PARAMETER(NF_EVARS_NELEMS_MULTIDEFINE=NF_EMULTIDEFINE_VAR_NUM)
      PARAMETER(NF_EVARS_NAME_MULTIDEFINE  =NF_EMULTIDEFINE_VAR_NAME)
      PARAMETER(NF_EVARS_NDIMS_MULTIDEFINE =NF_EMULTIDEFINE_VAR_NDIMS)
      PARAMETER(NF_EVARS_DIMIDS_MULTIDEFINE=NF_EMULTIDEFINE_VAR_DIMIDS)
      PARAMETER(NF_EVARS_TYPE_MULTIDEFINE  =NF_EMULTIDEFINE_VAR_TYPE)
      PARAMETER(NF_EVARS_LEN_MULTIDEFINE   =NF_EMULTIDEFINE_VAR_LEN)
      PARAMETER(NF_ENUMRECS_MULTIDEFINE    =NF_EMULTIDEFINE_NUMRECS)
      PARAMETER(NF_EVARS_BEGIN_MULTIDEFINE =NF_EMULTIDEFINE_VAR_BEGIN)

! error handling modes:
!
      integer nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!
! functions in the fortran interface
!

      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!

!     read/write, 0 => readonly
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef
      parameter(nccreat = 2)
!     on create destroy existing file
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!
!     'mode' arguments for nccreate and ncopen
!
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!
!     'size' argument to ncdimdef for an unlimited dimension
!
      integer ncunlim
      parameter(ncunlim = 0)

!
!     attribute id to put/get a global attribute
!
      parameter(ncglobal  = 0)

!
!     advisory maximums:
!
      parameter(maxncop = 32)
      parameter(maxncdim = 100)
      parameter(maxncatt = 2000)
      parameter(maxncvar = 2000)
!     not enforced
      parameter(maxncnam = 128)
      parameter(maxvdims = maxncdim)

!
!     global netcdf error status variable
!     initialized in error.c
!

!     no error
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument
      parameter(nceinval = nf_einval)
!     write to read only
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode
      parameter(ncenotin = nf_enotindefine)
!     operation not allowed in define mode
      parameter(nceindef = nf_eindefine)
!     coordinates out of domain
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use
      parameter(ncename = nf_enameinuse)
!     attribute not found
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid
      parameter(nceglob = nf_eglobal)
!     not a netcdf file
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)

! NULL request for non-blocking I/O APIs
      integer NF_REQ_NULL
      PARAMETER (NF_REQ_NULL = -1)

! indicate to flush all pending non-blocking requests
      integer NF_REQ_ALL
      PARAMETER (NF_REQ_ALL = -1)
      integer NF_GET_REQ_ALL
      PARAMETER (NF_GET_REQ_ALL = -2)
      integer NF_PUT_REQ_ALL
      PARAMETER (NF_PUT_REQ_ALL = -3)

!
! PnetCDF APIs
!

!
! miscellaneous routines:
!
      character*80  nfmpi_inq_libvers
      character*80  nfmpi_strerror
      character*80  nfmpi_strerrno

      external      nfmpi_inq_libvers
      external      nfmpi_strerror
      external      nfmpi_strerrno

      logical       nfmpi_issyserr
      external      nfmpi_issyserr

!
! control routines:
!
      integer  nfmpi_create
      integer  nfmpi_open
      integer  nfmpi_inq_format
      integer  nfmpi_inq_file_format
      integer  nfmpi_inq_file_info
      integer  nfmpi_get_file_info
      integer  nfmpi_delete
      integer  nfmpi_enddef
      integer  nfmpi__enddef
      integer  nfmpi_redef
      integer  nfmpi_set_default_format
      integer  nfmpi_inq_default_format
      integer  nfmpi_sync
      integer  nfmpi_abort
      integer  nfmpi_close
      integer  nfmpi_set_fill
      integer  nfmpi_def_var_fill
      integer  nfmpi_inq_var_fill
      integer  nfmpi_fill_var_rec

      external nfmpi_create
      external nfmpi_open
      external nfmpi_inq_format
      external nfmpi_inq_file_format
      external nfmpi_inq_file_info
      external nfmpi_get_file_info
      external nfmpi_delete
      external nfmpi_enddef
      external nfmpi__enddef
      external nfmpi_redef
      external nfmpi_set_default_format
      external nfmpi_inq_default_format
      external nfmpi_sync
      external nfmpi_abort
      external nfmpi_close
      external nfmpi_set_fill
      external nfmpi_def_var_fill
      external nfmpi_inq_var_fill
      external nfmpi_fill_var_rec

!
! general inquiry routines:
!
      integer  nfmpi_inq
      integer  nfmpi_inq_ndims
      integer  nfmpi_inq_nvars
      integer  nfmpi_inq_num_rec_vars
      integer  nfmpi_inq_num_fix_vars
      integer  nfmpi_inq_natts
      integer  nfmpi_inq_unlimdim
      integer  nfmpi_inq_striping
      integer  nfmpi_inq_malloc_size
      integer  nfmpi_inq_malloc_max_size
      integer  nfmpi_inq_malloc_list
      integer  nfmpi_inq_files_opened
      integer  nfmpi_inq_recsize

      external nfmpi_inq
      external nfmpi_inq_ndims
      external nfmpi_inq_nvars
      external nfmpi_inq_num_rec_vars
      external nfmpi_inq_num_fix_vars
      external nfmpi_inq_natts
      external nfmpi_inq_unlimdim
      external nfmpi_inq_striping
      external nfmpi_inq_malloc_size
      external nfmpi_inq_malloc_max_size
      external nfmpi_inq_malloc_list
      external nfmpi_inq_files_opened
      external nfmpi_inq_recsize
!
! dimension routines:
!
      integer  nfmpi_def_dim
      integer  nfmpi_inq_dimid
      integer  nfmpi_inq_dim
      integer  nfmpi_inq_dimname
      integer  nfmpi_inq_dimlen
      integer  nfmpi_rename_dim

      external nfmpi_def_dim
      external nfmpi_inq_dimid
      external nfmpi_inq_dim
      external nfmpi_inq_dimname
      external nfmpi_inq_dimlen
      external nfmpi_rename_dim
!
! general attribute routines:
!
      integer  nfmpi_inq_att
      integer  nfmpi_inq_attid
      integer  nfmpi_inq_atttype
      integer  nfmpi_inq_attlen
      integer  nfmpi_inq_attname
      integer  nfmpi_copy_att
      integer  nfmpi_rename_att
      integer  nfmpi_del_att

      external nfmpi_inq_att
      external nfmpi_inq_attid
      external nfmpi_inq_atttype
      external nfmpi_inq_attlen
      external nfmpi_inq_attname
      external nfmpi_copy_att
      external nfmpi_rename_att
      external nfmpi_del_att

!
! attribute put/get routines:
!
      integer  nfmpi_put_att,        nfmpi_get_att
      integer  nfmpi_put_att_text,   nfmpi_get_att_text
      integer  nfmpi_put_att_int1,   nfmpi_get_att_int1
      integer  nfmpi_put_att_int2,   nfmpi_get_att_int2
      integer  nfmpi_put_att_int,    nfmpi_get_att_int
      integer  nfmpi_put_att_real,   nfmpi_get_att_real
      integer  nfmpi_put_att_double, nfmpi_get_att_double
      integer  nfmpi_put_att_int8,   nfmpi_get_att_int8

      external nfmpi_put_att,        nfmpi_get_att
      external nfmpi_put_att_text,   nfmpi_get_att_text
      external nfmpi_put_att_int1,   nfmpi_get_att_int1
      external nfmpi_put_att_int2,   nfmpi_get_att_int2
      external nfmpi_put_att_int,    nfmpi_get_att_int
      external nfmpi_put_att_real,   nfmpi_get_att_real
      external nfmpi_put_att_double, nfmpi_get_att_double
      external nfmpi_put_att_int8,   nfmpi_get_att_int8

!
! independent data mode routines:
!
      integer  nfmpi_begin_indep_data
      integer  nfmpi_end_indep_data

      external nfmpi_begin_indep_data
      external nfmpi_end_indep_data

!
! general variable routines:
!
      integer  nfmpi_def_var
      integer  nfmpi_inq_var
      integer  nfmpi_inq_varid
      integer  nfmpi_inq_varname
      integer  nfmpi_inq_vartype
      integer  nfmpi_inq_varndims
      integer  nfmpi_inq_vardimid
      integer  nfmpi_inq_varnatts
      integer  nfmpi_rename_var

      external nfmpi_def_var
      external nfmpi_inq_var
      external nfmpi_inq_varid
      external nfmpi_inq_varname
      external nfmpi_inq_vartype
      external nfmpi_inq_varndims
      external nfmpi_inq_vardimid
      external nfmpi_inq_varnatts
      external nfmpi_rename_var

!
! entire variable put/get routines:
!
      integer  nfmpi_put_var,        nfmpi_put_var_all
      integer  nfmpi_put_var_text,   nfmpi_put_var_text_all
      integer  nfmpi_put_var_int1,   nfmpi_put_var_int1_all
      integer  nfmpi_put_var_int2,   nfmpi_put_var_int2_all
      integer  nfmpi_put_var_int,    nfmpi_put_var_int_all
      integer  nfmpi_put_var_real,   nfmpi_put_var_real_all
      integer  nfmpi_put_var_double, nfmpi_put_var_double_all
      integer  nfmpi_put_var_int8,   nfmpi_put_var_int8_all

      external nfmpi_put_var,        nfmpi_put_var_all
      external nfmpi_put_var_text,   nfmpi_put_var_text_all
      external nfmpi_put_var_int1,   nfmpi_put_var_int1_all
      external nfmpi_put_var_int2,   nfmpi_put_var_int2_all
      external nfmpi_put_var_int,    nfmpi_put_var_int_all
      external nfmpi_put_var_real,   nfmpi_put_var_real_all
      external nfmpi_put_var_double, nfmpi_put_var_double_all
      external nfmpi_put_var_int8,   nfmpi_put_var_int8_all

      integer  nfmpi_get_var,        nfmpi_get_var_all
      integer  nfmpi_get_var_text,   nfmpi_get_var_text_all
      integer  nfmpi_get_var_int1,   nfmpi_get_var_int1_all
      integer  nfmpi_get_var_int2,   nfmpi_get_var_int2_all
      integer  nfmpi_get_var_int,    nfmpi_get_var_int_all
      integer  nfmpi_get_var_real,   nfmpi_get_var_real_all
      integer  nfmpi_get_var_double, nfmpi_get_var_double_all
      integer  nfmpi_get_var_int8,   nfmpi_get_var_int8_all

      external nfmpi_get_var,        nfmpi_get_var_all
      external nfmpi_get_var_text,   nfmpi_get_var_text_all
      external nfmpi_get_var_int1,   nfmpi_get_var_int1_all
      external nfmpi_get_var_int2,   nfmpi_get_var_int2_all
      external nfmpi_get_var_int,    nfmpi_get_var_int_all
      external nfmpi_get_var_real,   nfmpi_get_var_real_all
      external nfmpi_get_var_double, nfmpi_get_var_double_all
      external nfmpi_get_var_int8,   nfmpi_get_var_int8_all

!
! single element variable put/get routines:
!
      integer  nfmpi_put_var1,        nfmpi_put_var1_all
      integer  nfmpi_put_var1_text,   nfmpi_put_var1_text_all
      integer  nfmpi_put_var1_int1,   nfmpi_put_var1_int1_all
      integer  nfmpi_put_var1_int2,   nfmpi_put_var1_int2_all
      integer  nfmpi_put_var1_int,    nfmpi_put_var1_int_all
      integer  nfmpi_put_var1_real,   nfmpi_put_var1_real_all
      integer  nfmpi_put_var1_double, nfmpi_put_var1_double_all
      integer  nfmpi_put_var1_int8,   nfmpi_put_var1_int8_all

      external nfmpi_put_var1,        nfmpi_put_var1_all
      external nfmpi_put_var1_text,   nfmpi_put_var1_text_all
      external nfmpi_put_var1_int1,   nfmpi_put_var1_int1_all
      external nfmpi_put_var1_int2,   nfmpi_put_var1_int2_all
      external nfmpi_put_var1_int,    nfmpi_put_var1_int_all
      external nfmpi_put_var1_real,   nfmpi_put_var1_real_all
      external nfmpi_put_var1_double, nfmpi_put_var1_double_all
      external nfmpi_put_var1_int8,   nfmpi_put_var1_int8_all

      integer  nfmpi_get_var1,        nfmpi_get_var1_all
      integer  nfmpi_get_var1_text,   nfmpi_get_var1_text_all
      integer  nfmpi_get_var1_int1,   nfmpi_get_var1_int1_all
      integer  nfmpi_get_var1_int2,   nfmpi_get_var1_int2_all
      integer  nfmpi_get_var1_int,    nfmpi_get_var1_int_all
      integer  nfmpi_get_var1_real,   nfmpi_get_var1_real_all
      integer  nfmpi_get_var1_double, nfmpi_get_var1_double_all
      integer  nfmpi_get_var1_int8,   nfmpi_get_var1_int8_all

      external nfmpi_get_var1,        nfmpi_get_var1_all
      external nfmpi_get_var1_text,   nfmpi_get_var1_text_all
      external nfmpi_get_var1_int1,   nfmpi_get_var1_int1_all
      external nfmpi_get_var1_int2,   nfmpi_get_var1_int2_all
      external nfmpi_get_var1_int,    nfmpi_get_var1_int_all
      external nfmpi_get_var1_real,   nfmpi_get_var1_real_all
      external nfmpi_get_var1_double, nfmpi_get_var1_double_all
      external nfmpi_get_var1_int8,   nfmpi_get_var1_int8_all

!
! variable sub-array put/get routines:
!
      integer  nfmpi_put_vara,        nfmpi_put_vara_all
      integer  nfmpi_put_vara_text,   nfmpi_put_vara_text_all
      integer  nfmpi_put_vara_int1,   nfmpi_put_vara_int1_all
      integer  nfmpi_put_vara_int2,   nfmpi_put_vara_int2_all
      integer  nfmpi_put_vara_int,    nfmpi_put_vara_int_all
      integer  nfmpi_put_vara_real,   nfmpi_put_vara_real_all
      integer  nfmpi_put_vara_double, nfmpi_put_vara_double_all
      integer  nfmpi_put_vara_int8,   nfmpi_put_vara_int8_all

      external nfmpi_put_vara,        nfmpi_put_vara_all
      external nfmpi_put_vara_text,   nfmpi_put_vara_text_all
      external nfmpi_put_vara_int1,   nfmpi_put_vara_int1_all
      external nfmpi_put_vara_int2,   nfmpi_put_vara_int2_all
      external nfmpi_put_vara_int,    nfmpi_put_vara_int_all
      external nfmpi_put_vara_real,   nfmpi_put_vara_real_all
      external nfmpi_put_vara_double, nfmpi_put_vara_double_all
      external nfmpi_put_vara_int8,   nfmpi_put_vara_int8_all

      integer  nfmpi_get_vara,        nfmpi_get_vara_all
      integer  nfmpi_get_vara_text,   nfmpi_get_vara_text_all
      integer  nfmpi_get_vara_int1,   nfmpi_get_vara_int1_all
      integer  nfmpi_get_vara_int2,   nfmpi_get_vara_int2_all
      integer  nfmpi_get_vara_int,    nfmpi_get_vara_int_all
      integer  nfmpi_get_vara_real,   nfmpi_get_vara_real_all
      integer  nfmpi_get_vara_double, nfmpi_get_vara_double_all
      integer  nfmpi_get_vara_int8,   nfmpi_get_vara_int8_all

      external nfmpi_get_vara,        nfmpi_get_vara_all
      external nfmpi_get_vara_text,   nfmpi_get_vara_text_all
      external nfmpi_get_vara_int1,   nfmpi_get_vara_int1_all
      external nfmpi_get_vara_int2,   nfmpi_get_vara_int2_all
      external nfmpi_get_vara_int,    nfmpi_get_vara_int_all
      external nfmpi_get_vara_real,   nfmpi_get_vara_real_all
      external nfmpi_get_vara_double, nfmpi_get_vara_double_all
      external nfmpi_get_vara_int8,   nfmpi_get_vara_int8_all

!
! strided variable put/get routines:
!
      integer  nfmpi_put_vars,        nfmpi_put_vars_all
      integer  nfmpi_put_vars_text,   nfmpi_put_vars_text_all
      integer  nfmpi_put_vars_int1,   nfmpi_put_vars_int1_all
      integer  nfmpi_put_vars_int2,   nfmpi_put_vars_int2_all
      integer  nfmpi_put_vars_int,    nfmpi_put_vars_int_all
      integer  nfmpi_put_vars_real,   nfmpi_put_vars_real_all
      integer  nfmpi_put_vars_double, nfmpi_put_vars_double_all
      integer  nfmpi_put_vars_int8,   nfmpi_put_vars_int8_all

      external nfmpi_put_vars,        nfmpi_put_vars_all
      external nfmpi_put_vars_text,   nfmpi_put_vars_text_all
      external nfmpi_put_vars_int1,   nfmpi_put_vars_int1_all
      external nfmpi_put_vars_int2,   nfmpi_put_vars_int2_all
      external nfmpi_put_vars_int,    nfmpi_put_vars_int_all
      external nfmpi_put_vars_real,   nfmpi_put_vars_real_all
      external nfmpi_put_vars_double, nfmpi_put_vars_double_all
      external nfmpi_put_vars_int8,   nfmpi_put_vars_int8_all

      integer  nfmpi_get_vars,        nfmpi_get_vars_all
      integer  nfmpi_get_vars_text,   nfmpi_get_vars_text_all
      integer  nfmpi_get_vars_int1,   nfmpi_get_vars_int1_all
      integer  nfmpi_get_vars_int2,   nfmpi_get_vars_int2_all
      integer  nfmpi_get_vars_int,    nfmpi_get_vars_int_all
      integer  nfmpi_get_vars_real,   nfmpi_get_vars_real_all
      integer  nfmpi_get_vars_double, nfmpi_get_vars_double_all
      integer  nfmpi_get_vars_int8,   nfmpi_get_vars_int8_all

      external nfmpi_get_vars,        nfmpi_get_vars_all
      external nfmpi_get_vars_text,   nfmpi_get_vars_text_all
      external nfmpi_get_vars_int1,   nfmpi_get_vars_int1_all
      external nfmpi_get_vars_int2,   nfmpi_get_vars_int2_all
      external nfmpi_get_vars_int,    nfmpi_get_vars_int_all
      external nfmpi_get_vars_real,   nfmpi_get_vars_real_all
      external nfmpi_get_vars_double, nfmpi_get_vars_double_all
      external nfmpi_get_vars_int8,   nfmpi_get_vars_int8_all

!
! mapped variable put/get routines:
!
      integer  nfmpi_put_varm,        nfmpi_put_varm_all
      integer  nfmpi_put_varm_text,   nfmpi_put_varm_text_all
      integer  nfmpi_put_varm_int1,   nfmpi_put_varm_int1_all
      integer  nfmpi_put_varm_int2,   nfmpi_put_varm_int2_all
      integer  nfmpi_put_varm_int,    nfmpi_put_varm_int_all
      integer  nfmpi_put_varm_real,   nfmpi_put_varm_real_all
      integer  nfmpi_put_varm_double, nfmpi_put_varm_double_all
      integer  nfmpi_put_varm_int8,   nfmpi_put_varm_int8_all

      external nfmpi_put_varm,        nfmpi_put_varm_all
      external nfmpi_put_varm_text,   nfmpi_put_varm_text_all
      external nfmpi_put_varm_int1,   nfmpi_put_varm_int1_all
      external nfmpi_put_varm_int2,   nfmpi_put_varm_int2_all
      external nfmpi_put_varm_int,    nfmpi_put_varm_int_all
      external nfmpi_put_varm_real,   nfmpi_put_varm_real_all
      external nfmpi_put_varm_double, nfmpi_put_varm_double_all
      external nfmpi_put_varm_int8,   nfmpi_put_varm_int8_all

      integer  nfmpi_get_varm,        nfmpi_get_varm_all
      integer  nfmpi_get_varm_text,   nfmpi_get_varm_text_all
      integer  nfmpi_get_varm_int1,   nfmpi_get_varm_int1_all
      integer  nfmpi_get_varm_int2,   nfmpi_get_varm_int2_all
      integer  nfmpi_get_varm_int,    nfmpi_get_varm_int_all
      integer  nfmpi_get_varm_real,   nfmpi_get_varm_real_all
      integer  nfmpi_get_varm_double, nfmpi_get_varm_double_all
      integer  nfmpi_get_varm_int8,   nfmpi_get_varm_int8_all

      external nfmpi_get_varm,        nfmpi_get_varm_all
      external nfmpi_get_varm_text,   nfmpi_get_varm_text_all
      external nfmpi_get_varm_int1,   nfmpi_get_varm_int1_all
      external nfmpi_get_varm_int2,   nfmpi_get_varm_int2_all
      external nfmpi_get_varm_int,    nfmpi_get_varm_int_all
      external nfmpi_get_varm_real,   nfmpi_get_varm_real_all
      external nfmpi_get_varm_double, nfmpi_get_varm_double_all
      external nfmpi_get_varm_int8,   nfmpi_get_varm_int8_all

!
! Non-blocking APIs
!
! entire variable iput/iget routines:
!
      integer  nfmpi_iput_var
      integer  nfmpi_iput_var_text
      integer  nfmpi_iput_var_int1
      integer  nfmpi_iput_var_int2
      integer  nfmpi_iput_var_int
      integer  nfmpi_iput_var_real
      integer  nfmpi_iput_var_double
      integer  nfmpi_iput_var_int8

      external nfmpi_iput_var
      external nfmpi_iput_var_text
      external nfmpi_iput_var_int1
      external nfmpi_iput_var_int2
      external nfmpi_iput_var_int
      external nfmpi_iput_var_real
      external nfmpi_iput_var_double
      external nfmpi_iput_var_int8

      integer  nfmpi_iget_var
      integer  nfmpi_iget_var_text
      integer  nfmpi_iget_var_int1
      integer  nfmpi_iget_var_int2
      integer  nfmpi_iget_var_int
      integer  nfmpi_iget_var_real
      integer  nfmpi_iget_var_double
      integer  nfmpi_iget_var_int8

      external nfmpi_iget_var
      external nfmpi_iget_var_text
      external nfmpi_iget_var_int1
      external nfmpi_iget_var_int2
      external nfmpi_iget_var_int
      external nfmpi_iget_var_real
      external nfmpi_iget_var_double
      external nfmpi_iget_var_int8

!
! Nonblocking single-element variable iput/iget routines:
!
      integer  nfmpi_iput_var1
      integer  nfmpi_iput_var1_text
      integer  nfmpi_iput_var1_int1
      integer  nfmpi_iput_var1_int2
      integer  nfmpi_iput_var1_int
      integer  nfmpi_iput_var1_real
      integer  nfmpi_iput_var1_double
      integer  nfmpi_iput_var1_int8

      external nfmpi_iput_var1
      external nfmpi_iput_var1_text
      external nfmpi_iput_var1_int1
      external nfmpi_iput_var1_int2
      external nfmpi_iput_var1_int
      external nfmpi_iput_var1_real
      external nfmpi_iput_var1_double
      external nfmpi_iput_var1_int8

      integer  nfmpi_iget_var1
      integer  nfmpi_iget_var1_text
      integer  nfmpi_iget_var1_int1
      integer  nfmpi_iget_var1_int2
      integer  nfmpi_iget_var1_int
      integer  nfmpi_iget_var1_real
      integer  nfmpi_iget_var1_double
      integer  nfmpi_iget_var1_int8

      external nfmpi_iget_var1
      external nfmpi_iget_var1_text
      external nfmpi_iget_var1_int1
      external nfmpi_iget_var1_int2
      external nfmpi_iget_var1_int
      external nfmpi_iget_var1_real
      external nfmpi_iget_var1_double
      external nfmpi_iget_var1_int8

!
! Nonblocking subarray variable iput/iget routines:
!
      integer  nfmpi_iput_vara
      integer  nfmpi_iput_vara_text
      integer  nfmpi_iput_vara_int1
      integer  nfmpi_iput_vara_int2
      integer  nfmpi_iput_vara_int
      integer  nfmpi_iput_vara_real
      integer  nfmpi_iput_vara_double
      integer  nfmpi_iput_vara_int8

      external nfmpi_iput_vara
      external nfmpi_iput_vara_text
      external nfmpi_iput_vara_int1
      external nfmpi_iput_vara_int2
      external nfmpi_iput_vara_int
      external nfmpi_iput_vara_real
      external nfmpi_iput_vara_double
      external nfmpi_iput_vara_int8

      integer  nfmpi_iget_vara
      integer  nfmpi_iget_vara_text
      integer  nfmpi_iget_vara_int1
      integer  nfmpi_iget_vara_int2
      integer  nfmpi_iget_vara_int
      integer  nfmpi_iget_vara_real
      integer  nfmpi_iget_vara_double
      integer  nfmpi_iget_vara_int8

      external nfmpi_iget_vara
      external nfmpi_iget_vara_text
      external nfmpi_iget_vara_int1
      external nfmpi_iget_vara_int2
      external nfmpi_iget_vara_int
      external nfmpi_iget_vara_real
      external nfmpi_iget_vara_double
      external nfmpi_iget_vara_int8

!
! Nonblocking strided variable iput/iget routines:
!
      integer  nfmpi_iput_vars
      integer  nfmpi_iput_vars_text
      integer  nfmpi_iput_vars_int1
      integer  nfmpi_iput_vars_int2
      integer  nfmpi_iput_vars_int
      integer  nfmpi_iput_vars_real
      integer  nfmpi_iput_vars_double
      integer  nfmpi_iput_vars_int8

      external nfmpi_iput_vars
      external nfmpi_iput_vars_text
      external nfmpi_iput_vars_int1
      external nfmpi_iput_vars_int2
      external nfmpi_iput_vars_int
      external nfmpi_iput_vars_real
      external nfmpi_iput_vars_double
      external nfmpi_iput_vars_int8

      integer  nfmpi_iget_vars
      integer  nfmpi_iget_vars_text
      integer  nfmpi_iget_vars_int1
      integer  nfmpi_iget_vars_int2
      integer  nfmpi_iget_vars_int
      integer  nfmpi_iget_vars_real
      integer  nfmpi_iget_vars_double
      integer  nfmpi_iget_vars_int8

      external nfmpi_iget_vars
      external nfmpi_iget_vars_text
      external nfmpi_iget_vars_int1
      external nfmpi_iget_vars_int2
      external nfmpi_iget_vars_int
      external nfmpi_iget_vars_real
      external nfmpi_iget_vars_double
      external nfmpi_iget_vars_int8

!
! Nonblocking mapped variable iput/iget routines:
!
      integer  nfmpi_iput_varm
      integer  nfmpi_iput_varm_text
      integer  nfmpi_iput_varm_int1
      integer  nfmpi_iput_varm_int2
      integer  nfmpi_iput_varm_int
      integer  nfmpi_iput_varm_real
      integer  nfmpi_iput_varm_double
      integer  nfmpi_iput_varm_int8

      external nfmpi_iput_varm
      external nfmpi_iput_varm_text
      external nfmpi_iput_varm_int1
      external nfmpi_iput_varm_int2
      external nfmpi_iput_varm_int
      external nfmpi_iput_varm_real
      external nfmpi_iput_varm_double
      external nfmpi_iput_varm_int8

      integer  nfmpi_iget_varm
      integer  nfmpi_iget_varm_text
      integer  nfmpi_iget_varm_int1
      integer  nfmpi_iget_varm_int2
      integer  nfmpi_iget_varm_int
      integer  nfmpi_iget_varm_real
      integer  nfmpi_iget_varm_double
      integer  nfmpi_iget_varm_int8

      external nfmpi_iget_varm
      external nfmpi_iget_varm_text
      external nfmpi_iget_varm_int1
      external nfmpi_iget_varm_int2
      external nfmpi_iget_varm_int
      external nfmpi_iget_varm_real
      external nfmpi_iget_varm_double
      external nfmpi_iget_varm_int8

!
! Nonblocking entire variable bput routines:
!
      integer  nfmpi_bput_var
      integer  nfmpi_bput_var_text
      integer  nfmpi_bput_var_int1
      integer  nfmpi_bput_var_int2
      integer  nfmpi_bput_var_int
      integer  nfmpi_bput_var_real
      integer  nfmpi_bput_var_double
      integer  nfmpi_bput_var_int8

      external nfmpi_bput_var
      external nfmpi_bput_var_text
      external nfmpi_bput_var_int1
      external nfmpi_bput_var_int2
      external nfmpi_bput_var_int
      external nfmpi_bput_var_real
      external nfmpi_bput_var_double
      external nfmpi_bput_var_int8

!
! Nonblocking single element variable bput routines:
!
      integer  nfmpi_bput_var1
      integer  nfmpi_bput_var1_text
      integer  nfmpi_bput_var1_int1
      integer  nfmpi_bput_var1_int2
      integer  nfmpi_bput_var1_int
      integer  nfmpi_bput_var1_real
      integer  nfmpi_bput_var1_double
      integer  nfmpi_bput_var1_int8

      external nfmpi_bput_var1
      external nfmpi_bput_var1_text
      external nfmpi_bput_var1_int1
      external nfmpi_bput_var1_int2
      external nfmpi_bput_var1_int
      external nfmpi_bput_var1_real
      external nfmpi_bput_var1_double
      external nfmpi_bput_var1_int8

!
! Nonblocking subarray variable bput routines:
!
      integer  nfmpi_bput_vara
      integer  nfmpi_bput_vara_text
      integer  nfmpi_bput_vara_int1
      integer  nfmpi_bput_vara_int2
      integer  nfmpi_bput_vara_int
      integer  nfmpi_bput_vara_real
      integer  nfmpi_bput_vara_double
      integer  nfmpi_bput_vara_int8

      external nfmpi_bput_vara
      external nfmpi_bput_vara_text
      external nfmpi_bput_vara_int1
      external nfmpi_bput_vara_int2
      external nfmpi_bput_vara_int
      external nfmpi_bput_vara_real
      external nfmpi_bput_vara_double
      external nfmpi_bput_vara_int8

!
! Nonblocking strided variable bput routines:
!
      integer  nfmpi_bput_vars
      integer  nfmpi_bput_vars_text
      integer  nfmpi_bput_vars_int1
      integer  nfmpi_bput_vars_int2
      integer  nfmpi_bput_vars_int
      integer  nfmpi_bput_vars_real
      integer  nfmpi_bput_vars_double
      integer  nfmpi_bput_vars_int8

      external nfmpi_bput_vars
      external nfmpi_bput_vars_text
      external nfmpi_bput_vars_int1
      external nfmpi_bput_vars_int2
      external nfmpi_bput_vars_int
      external nfmpi_bput_vars_real
      external nfmpi_bput_vars_double
      external nfmpi_bput_vars_int8

!
! Nonblocking mapped variable bput routines:
!
      integer  nfmpi_bput_varm
      integer  nfmpi_bput_varm_text
      integer  nfmpi_bput_varm_int1
      integer  nfmpi_bput_varm_int2
      integer  nfmpi_bput_varm_int
      integer  nfmpi_bput_varm_real
      integer  nfmpi_bput_varm_double
      integer  nfmpi_bput_varm_int8

      external nfmpi_bput_varm
      external nfmpi_bput_varm_text
      external nfmpi_bput_varm_int1
      external nfmpi_bput_varm_int2
      external nfmpi_bput_varm_int
      external nfmpi_bput_varm_real
      external nfmpi_bput_varm_double
      external nfmpi_bput_varm_int8

!
! Nonblocking control APIs
!
      integer  nfmpi_wait
      integer  nfmpi_wait_all
      integer  nfmpi_cancel

      external nfmpi_wait
      external nfmpi_wait_all
      external nfmpi_cancel

      integer  nfmpi_buffer_attach
      integer  nfmpi_buffer_detach
      integer  nfmpi_inq_buffer_usage
      integer  nfmpi_inq_buffer_size
      integer  nfmpi_inq_put_size
      integer  nfmpi_inq_get_size
      integer  nfmpi_inq_header_size
      integer  nfmpi_inq_header_extent
      integer  nfmpi_inq_varoffset
      integer  nfmpi_inq_nreqs
      integer  nfmpi_inq_path

      external nfmpi_buffer_attach
      external nfmpi_buffer_detach
      external nfmpi_inq_buffer_usage
      external nfmpi_inq_buffer_size
      external nfmpi_inq_put_size
      external nfmpi_inq_get_size
      external nfmpi_inq_header_size
      external nfmpi_inq_header_extent
      external nfmpi_inq_varoffset
      external nfmpi_inq_nreqs
      external nfmpi_inq_path

!
! varn routines:
!
      integer  nfmpi_put_varn,        nfmpi_put_varn_all
      integer  nfmpi_put_varn_text,   nfmpi_put_varn_text_all
      integer  nfmpi_put_varn_int1,   nfmpi_put_varn_int1_all
      integer  nfmpi_put_varn_int2,   nfmpi_put_varn_int2_all
      integer  nfmpi_put_varn_int,    nfmpi_put_varn_int_all
      integer  nfmpi_put_varn_real,   nfmpi_put_varn_real_all
      integer  nfmpi_put_varn_double, nfmpi_put_varn_double_all
      integer  nfmpi_put_varn_int8,   nfmpi_put_varn_int8_all

      external nfmpi_put_varn,        nfmpi_put_varn_all
      external nfmpi_put_varn_text,   nfmpi_put_varn_text_all
      external nfmpi_put_varn_int1,   nfmpi_put_varn_int1_all
      external nfmpi_put_varn_int2,   nfmpi_put_varn_int2_all
      external nfmpi_put_varn_int,    nfmpi_put_varn_int_all
      external nfmpi_put_varn_real,   nfmpi_put_varn_real_all
      external nfmpi_put_varn_double, nfmpi_put_varn_double_all
      external nfmpi_put_varn_int8,   nfmpi_put_varn_int8_all

      integer  nfmpi_get_varn,        nfmpi_get_varn_all
      integer  nfmpi_get_varn_text,   nfmpi_get_varn_text_all
      integer  nfmpi_get_varn_int1,   nfmpi_get_varn_int1_all
      integer  nfmpi_get_varn_int2,   nfmpi_get_varn_int2_all
      integer  nfmpi_get_varn_int,    nfmpi_get_varn_int_all
      integer  nfmpi_get_varn_real,   nfmpi_get_varn_real_all
      integer  nfmpi_get_varn_double, nfmpi_get_varn_double_all
      integer  nfmpi_get_varn_int8,   nfmpi_get_varn_int8_all

      external nfmpi_get_varn,        nfmpi_get_varn_all
      external nfmpi_get_varn_text,   nfmpi_get_varn_text_all
      external nfmpi_get_varn_int1,   nfmpi_get_varn_int1_all
      external nfmpi_get_varn_int2,   nfmpi_get_varn_int2_all
      external nfmpi_get_varn_int,    nfmpi_get_varn_int_all
      external nfmpi_get_varn_real,   nfmpi_get_varn_real_all
      external nfmpi_get_varn_double, nfmpi_get_varn_double_all
      external nfmpi_get_varn_int8,   nfmpi_get_varn_int8_all

!
! Nonblocking varn routines:
!
      integer  nfmpi_iput_varn
      integer  nfmpi_iput_varn_text
      integer  nfmpi_iput_varn_int1
      integer  nfmpi_iput_varn_int2
      integer  nfmpi_iput_varn_int
      integer  nfmpi_iput_varn_real
      integer  nfmpi_iput_varn_double
      integer  nfmpi_iput_varn_int8

      external nfmpi_iput_varn
      external nfmpi_iput_varn_text
      external nfmpi_iput_varn_int1
      external nfmpi_iput_varn_int2
      external nfmpi_iput_varn_int
      external nfmpi_iput_varn_real
      external nfmpi_iput_varn_double
      external nfmpi_iput_varn_int8

      integer  nfmpi_iget_varn
      integer  nfmpi_iget_varn_text
      integer  nfmpi_iget_varn_int1
      integer  nfmpi_iget_varn_int2
      integer  nfmpi_iget_varn_int
      integer  nfmpi_iget_varn_real
      integer  nfmpi_iget_varn_double
      integer  nfmpi_iget_varn_int8

      external nfmpi_iget_varn
      external nfmpi_iget_varn_text
      external nfmpi_iget_varn_int1
      external nfmpi_iget_varn_int2
      external nfmpi_iget_varn_int
      external nfmpi_iget_varn_real
      external nfmpi_iget_varn_double
      external nfmpi_iget_varn_int8

      integer  nfmpi_bput_varn
      integer  nfmpi_bput_varn_text
      integer  nfmpi_bput_varn_int1
      integer  nfmpi_bput_varn_int2
      integer  nfmpi_bput_varn_int
      integer  nfmpi_bput_varn_real
      integer  nfmpi_bput_varn_double
      integer  nfmpi_bput_varn_int8

      external nfmpi_bput_varn
      external nfmpi_bput_varn_text
      external nfmpi_bput_varn_int1
      external nfmpi_bput_varn_int2
      external nfmpi_bput_varn_int
      external nfmpi_bput_varn_real
      external nfmpi_bput_varn_double
      external nfmpi_bput_varn_int8

!
! vard routines:
!
      integer  nfmpi_put_vard, nfmpi_put_vard_all
      integer  nfmpi_get_vard, nfmpi_get_vard_all

      external nfmpi_put_vard, nfmpi_put_vard_all
      external nfmpi_get_vard, nfmpi_get_vard_all


  !
  !  Attribute functions
  !
  public :: pio_def_var,   &
       pio_inq_attname,    & 
       pio_inq_att,        &
       pio_inq_attlen,     &
       pio_inq_varid,      &
       pio_inq_varname,    &
       pio_inq_vartype,    &
       pio_inq_varndims,   &
       pio_inq_vardimid,   &
       pio_inq_varnatts,   &
       pio_inquire_variable
!>
!! \defgroup PIO_def_var
!<
  interface pio_def_var
     module procedure &
          def_var_0d, &
          def_var_md
  end interface

!>
!! \defgroup PIO_inq_varid
!<
  interface pio_inq_varid
     module procedure inq_varid_vid, &
          inq_varid_vardesc
  end interface
!>
!! \defgroup PIO_inq_att
!<
  interface pio_inq_att
     module procedure inq_att_vid, &
          inq_att_vardesc
  end interface

!>
!! \defgroup PIO_inq_attlen
!<
  interface pio_inq_attlen
     module procedure inq_attlen_vid, &
          inq_attlen_vardesc
  end interface

!>
!! \defgroup PIO_inq_attname
!<
  interface pio_inq_attname
     module procedure inq_attname_vid, &
          inq_attname_vardesc
  end interface

!>
!! \defgroup PIO_inq_varname
!<
  interface pio_inq_varname
     module procedure inq_varname_vid, inq_varname_vdesc
  end interface

!>
!! \defgroup PIO_inq_varndims
!<
  interface pio_inq_varndims
     module procedure inq_varndims_vid, inq_varndims_vdesc
  end interface

!>
!! \defgroup PIO_inq_varnatts
!<
  interface pio_inq_varnatts
     module procedure inq_varnatts_vid, inq_varnatts_vdesc
  end interface

!>
!! \defgroup PIO_inq_vardimid
!<
  interface pio_inq_vardimid
     module procedure inq_vardimid_vid, inq_vardimid_vdesc
  end interface

!>
!! \defgroup PIO_inq_vartype
!<
  interface pio_inq_vartype
     module procedure inq_vartype_vid, inq_vartype_vdesc
  end interface

!>
!! \defgroup PIO_inquire_variable
!<
  interface pio_inquire_variable
     module procedure inquire_variable_vid, inquire_variable_vdesc
  end interface

!>
!! @defgroup PIO_def_dim
!<
   public :: PIO_def_dim

!>
!! @defgroup PIO_enddef
!<
  public :: PIO_enddef

!>
!! \defgroup PIO_redef
!<
  public :: PIO_redef

!> 
!! \defgroup PIO_inquire
!<
   public :: PIO_inquire

!>
!! \defgroup PIO_inq_dimid
!<
   public :: PIO_inq_dimid

!>
!! \defgroup PIO_inq_dimname
!<
   public :: PIO_inq_dimname

!>
!! \defgroup PIO_inq_dimlen
!<
   public :: PIO_inq_dimlen

!>
!! \defgroup PIO_inquire_dimension
!<
  public :: PIO_inquire_dimension

!> 
!! \defgroup PIO_copy_att
!<
  public :: PIO_copy_att
contains 

!>
!! @public 
!! @ingroup PIO_inquire
!! @brief Gets metadata information for netcdf file.
!! @details
!! @param File @copydoc file_desc_t
!! @param nDimensions :  Number of dimensions defined for the netcdf file
!! @param nVariables : Number of variables defined for the netcdf file 
!! @param nAttributes : Number of attributes defined for the netcdf file
!! @param unlimitedDimID : the Unlimited dimension ID
!! @retval ierr @copydoc error_return
!>
  integer function pio_inquire(File,nDimensions,nVariables,nAttributes,unlimitedDimID) result(ierr)
    type (File_desc_t), intent(in) :: File

    integer, optional, intent(out) :: &
         nDimensions,  &! number of dimensions
         nVariables,   &! number of variables
         nAttributes,  & ! number of global attributes
         unlimitedDimID ! ID of unlimited dimension
    integer :: vals(4)

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg
    type(iosystem_desc_t), pointer :: ios

    ierr=PIO_noerr
    vals(:) = -1

    ios => File%iosystem

    if(ios%async_interface .and. .not. ios%ioproc) then
       msg=PIO_MSG_INQUIRE
       if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
    end if

    iotype = File%iotype
    
    if(ios%IOproc) then
       select case(iotype)
          
       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq( File%fh,vals(1),vals(2), &
               vals(3),vals(4))


       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire( File%fh,vals(1),vals(2), &
               vals(3),vals(4))
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire( File%fh,vals(1),vals(2), &
                  vals(3),vals(4))
          endif
          
          if(ios%num_iotasks>1) then
             call MPI_BCAST(vals,4,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if


       case default
          call bad_iotype(iotype,"nf_mod.F90",269)
          
       end select
    endif

    call check_netcdf(File, ierr, "nf_mod.F90",274)

    if(file%iosystem%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(vals,4,MPI_INTEGER,ios%IOMaster, ios%my_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

    if(present(nDimensions)) then	
       ndimensions = vals(1)
    endif
    if(present(nVariables)) then	
       nVariables = vals(2)
    endif
    if(present(nAttributes)) then	
       nAttributes = vals(3)
    endif
    if(present(unlimitedDimID)) then	
       unlimitedDimID = vals(4)
    endif
    
  end function pio_inquire

!>
!! @public 
!! @ingroup PIO_inq_att
!! @brief Gets information about attributes
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable identifier
!! @param name : Name of the attribute 
!! @param xtype : The type of attribute
!! @param len : The length of the attribute 
!! @retval ierr @copydoc error_return
!>
  integer function inq_att_vid(File,varid,name,xtype,len) result(ierr)


    type (File_desc_t), intent(inout) :: File
    integer(i4), intent(in)           :: varid
    character(len=*), intent(in)      :: name
    integer, intent(out)              :: xtype
    integer, intent(out)              :: len !Attribute length

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg, nlen
    integer(kind=PIO_Offset) :: clen
    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem
    iotype = File%iotype
    ierr=PIO_noerr
    nlen = len_trim(name)

    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_ATT
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(nlen,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(name,nlen,MPI_CHARACTER,ios%CompMaster, ios%my_comm , mpierr)
    end if



    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_att(File%fh,varid,name(1:nlen),xtype,clen)

          len = INT(clen,kind=i4)

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_attribute( File%fh,varid,name(1:nlen), &
               xtype=xtype,len=len)          
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)

          if (ios%io_rank==0) then
             ierr=nf90_inquire_attribute( File%fh,varid,name(1:nlen), &
                  xtype=xtype,len=len)
          endif

          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(xtype,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
             call MPI_BCAST(len,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",372)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",376)
    if(ios%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(xtype,1,MPI_INTEGER,ios%IOMaster, ios%my_comm , mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
       call MPI_BCAST(len,1,MPI_INTEGER,ios%IOMaster, ios%my_comm  , mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_att_vid


!>
!! @public 
!! @ingroup PIO_inq_att
!! @brief  Gets information about attributes
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param name : Name of the attribute 
!! @param xtype : The type of attribute
!! @param len : The length of the attribute 
!! @retval ierr @copydoc error_return
!>
  integer function inq_att_vardesc(File,vardesc,name,xtype,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    type(var_desc_t), intent(in)           :: vardesc
    character(len=*), intent(in)      :: name
    integer, intent(out)              :: xtype
    integer, intent(out)              :: len !Attribute length

    ierr = pio_inq_att(file, vardesc%varid, name, xtype, len)

  end function inq_att_vardesc

!>
!! @public 
!! @ingroup PIO_inq_attlen
!! @brief Gets the attribute length 
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : attribute id
!! @param name : name of attribute
!! @param len : Length of attribute
!! @retval ierr @copydoc error_return
!>
  integer function inq_attlen_vid(File,varid,name,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    integer(i4), intent(in)            :: varid
    character(len=*), intent(in)      :: name
    integer, intent(out)              :: len !Attribute length


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg, nlen
    integer(kind=PIO_Offset) :: clen
    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem

    iotype = File%iotype
    ierr=PIO_noerr
    nlen = len_trim(name)

    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_ATTLEN
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(nlen,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(name,nlen,MPI_CHARACTER,ios%CompMaster, ios%my_comm , mpierr)
       
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_attlen(File%fh,varid,name(1:nlen),clen)
          len = INT(clen,kind=i4)

       case(pio_iotype_netcdf4p)
             ierr=nf90_inquire_attribute( File%fh,varid,name(1:nlen), &
                  len=len)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_attribute( File%fh,varid,name(1:nlen), &
                  len=len)
          endif

          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(len,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if


       case default
          call bad_iotype(iotype,"nf_mod.F90",481)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",485)
    if(ios%async_interface.or.ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(len,1,MPI_INTEGER,ios%IOMaster,ios%my_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_attlen_vid

!>
!! @public 
!! @ingroup PIO_inq_attlen
!! @brief  Gets the attribute length 
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param name : name of attribute
!! @param len : Length of attribute
!! @retval ierr @copydoc error_return
!>
  integer function inq_attlen_vardesc(File,vardesc,name,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    type (Var_desc_t), intent(in)            :: vardesc
    character(len=*), intent(in)      :: name
    integer, intent(out),optional     :: len !Attribute length

    ierr = pio_inq_attlen(file, vardesc%varid, name, len)

  end function inq_attlen_vardesc

!> 
!! @public 
!! @ingroup PIO_inq_attname
!! @brief Returns the name of a netcdf attribute 
!! @details
!! @param File @copydoc file_desc_t
!! @param varid :  The variable ID 
!! @param attnum : Attribute number returned from function ????
!! @param name   : Name of the returned attribute
!! @retval ierr @copydoc error_return
!<
  integer function inq_attname_vid(File,varid,attnum,name) result(ierr)

    type (File_desc_t), intent(inout) :: File
    integer(i4), intent(in)           :: varid
    integer, intent(in)              :: attnum !Attribute number
    character(len=*), intent(out)     :: name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg
    type(iosystem_desc_t), pointer :: ios
    character(len=PIO_MAX_NAME) :: tmpname

    ios => File%iosystem

    iotype = File%iotype
    ierr=PIO_noerr
    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_ATTNAME
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(attnum,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_attname(File%fh,varid,attnum,tmpname)


       case(pio_iotype_netcdf4p)
          ierr=nf90_inq_attname(File%fh,varid,attnum,tmpname)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inq_attname(File%fh,varid,attnum,tmpname)
             if(Debug) print *,"nf_mod.F90",570,name
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(tmpname,PIO_MAX_NAME,MPI_CHARACTER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if


       case default
          call bad_iotype(iotype,"nf_mod.F90",580)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",584)
    if(ios%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(tmpname,PIO_MAX_NAME,MPI_CHARACTER,ios%IOMaster,ios%my_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
    name = tmpname(1:len_trim(tmpname))
  end function inq_attname_vid

!> 
!! @public 
!! @ingroup PIO_inq_attname
!! @brief  Returns the name of a netcdf attribute.
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t 
!! @param attnum : Attribute number returned from function ????
!! @param name   : Name of the returned attribute
!! @retval ierr @copydoc error_return
!<
  integer function inq_attname_vardesc(File,vardesc,attnum,name) result(ierr)
    type (File_desc_t), intent(inout) :: File
    type(var_desc_t), intent(in)           :: vardesc
    integer, intent(in)              :: attnum !Attribute number
    character(len=*), intent(out)     :: name

    ierr = pio_inq_attname(file, vardesc%varid, attnum, name)

  end function inq_attname_vardesc

!> 
!! @public 
!! @ingroup PIO_inq_varid
!! @brief  Returns the ID of a netcdf variable given its name 
!! @details
!! @param File @copydoc file_desc_t
!! @param name : Name of the returned attribute
!! @param varid : variable ID
!! @retval ierr @copydoc error_return
!<
  integer function inq_varid_vid(File,name,varid) result(ierr)

    type (File_desc_t), intent(in)   :: File
    character(len=*), intent(in)     :: name
    integer(i4), intent(out)       :: varid
    integer :: ierr2

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg, nlen
    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem

    iotype = File%iotype
    ierr=PIO_noerr
    nlen = len_trim(name)

    if(ios%async_interface) then
       if( .not. ios%ioproc ) then
          msg=PIO_MSG_INQ_VARID
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       
       call MPI_BCAST(nlen,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(name,nlen,MPI_CHARACTER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_varid(File%fh,name(1:nlen),varid)

       case(pio_iotype_netcdf4p)
             ierr=nf90_inq_varid(File%fh,name(1:nlen),varid)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inq_varid(File%fh,name(1:nlen),varid)
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(varid,1,MPI_INTEGER,0,ios%IO_comm,ierr2)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",674)

       end select
    endif

    call check_netcdf(File, ierr,"nf_mod.F90",679)
    if(ios%async_interface.or.ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%IOMaster,ios%my_comm,ierr2)
    end if

  end function inq_varid_vid

!> 
!! @public 
!! @ingroup PIO_inq_varid
!! @brief Returns the ID of a netcdf variable given its name 
!! @details
!! @param File @copydoc file_desc_t
!! @param name   : Name of the returned attribute
!! @param vardesc @copydoc var_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function inq_varid_vardesc(File,name,vardesc) result(ierr)

    type (File_desc_t), intent(in)   :: File
    character(len=*), intent(in)     :: name
    type (Var_desc_t), intent(inout) :: vardesc

    ierr = pio_inq_varid(File, name, vardesc%varid)
    vardesc%rec=-1
    if(ierr==PIO_NOERR) then
       ierr = pio_inq_varndims(File, vardesc%varid, vardesc%ndims) ! needed for nfwrite
    end if
  end function inq_varid_vardesc

!>
!! @public 
!! @ingroup PIO_inq_varname
!! @brief Get the name associated with a variable
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param name : The name of the netcdf variable.
!! @retval ierr @copydoc error_return
!>
  integer function inq_varname_vdesc(File,vardesc,name) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in)    :: vardesc
    character(len=*), intent(out)    :: name
    
    ierr = pio_inq_varname(file,vardesc%varid,name)

  end function inq_varname_vdesc

!>
!! @public 
!! @ingroup PIO_inq_varname
!! @brief Get the name associated with a variable
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable id.
!! @param name : The name of the netcdf variable.
!! @retval ierr @copydoc error_return
!>
  integer function inq_varname_vid(File,varid,name) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer, intent(in) :: varid
    character(len=*), intent(out)    :: name
    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg, nlen

    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem
    iotype = File%iotype
    ierr=PIO_noerr
    nlen = len(name)
    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_VARNAME
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(nlen,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_varname(File%fh,varid,name(1:nlen))


       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_variable(File%fh,varid,name=name(1:nlen))
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,name=name(1:nlen))
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(name,nlen,MPI_CHARACTER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if


       case default
          call bad_iotype(iotype,"nf_mod.F90",789)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",793)
    if(ios%async_interface.or.ios%num_tasks>=ios%num_iotasks) then
       call MPI_BCAST(name,nlen,MPI_CHARACTER,ios%IOMaster,ios%my_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_varname_vid

!>
!! @public 
!! @ingroup PIO_inq_varndims
!! @brief Gets the number of dimension associated with a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The variable identifier
!! @param ndims : The number of dimensions for the variable 
!! @retval ierr @copydoc error_return
!>
  integer function inq_varndims_vid(File,varid,ndims) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer, intent(in) :: varid
    integer(i4), intent(out)    :: ndims


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg

    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem
    iotype = File%iotype
    ierr=PIO_noerr

    if(ios%async_interface) then
       if( .not. ios%ioproc ) then
          msg=PIO_MSG_INQ_VARNDIMS
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if
    
    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_varndims(File%fh,varid,ndims)

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_variable(File%fh,varid,ndims=ndims)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,ndims=ndims)
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(ndims,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",860)

       end select
    endif
    call check_netcdf(File,ierr,"nf_mod.F90",864)



    if(ios%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(ndims,1,MPI_INTEGER,ios%IOMaster,ios%my_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_varndims_vid

!>
!! @public 
!! @ingroup PIO_inq_varndims
!! @brief Gets the number of dimension associated with a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param ndims : The number of dimensions for the variable 
!! @retval ierr @copydoc error_return
!>
  integer function inq_varndims_vdesc(File,vardesc,ndims) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in) :: vardesc
    integer(i4), intent(out)    :: ndims

    ierr = pio_inq_varndims(File, vardesc%varid, ndims)
  end function inq_varndims_vdesc

!>
!! @public 
!! @ingroup PIO_inq_vartype
!! @brief Gets metadata information for netcdf file.
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable id
!! @param type : The type of variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_vartype_vid(File,varid,type) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer, intent(in) :: varid
    integer(i4), intent(out)    :: type


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg

    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem
    iotype = File%iotype
    ierr=PIO_noerr

    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_VARTYPE
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_vartype(File%fh,varid,type)

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_variable(File%fh,varid,xtype=type)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,xtype=type)
          endif

          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(type,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",953)

       end select
    endif
    call check_netcdf(File,ierr,"nf_mod.F90",957)
    if(file%iosystem%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(type,1,MPI_INTEGER,ios%IOMaster,ios%my_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_vartype_vid

!>
!! @public 
!! @ingroup PIO_inq_vartype
!! @brief Gets metadata information for netcdf file.
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param type : The type of variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_vartype_vdesc(File,vardesc,type) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in) :: vardesc
    integer(i4), intent(out)    :: type

    ierr = pio_inq_vartype(File, vardesc%varid, type)
  end function inq_vartype_vdesc

!>
!! @public 
!! @ingroup PIO_inq_vardimid
!! @brief returns the dimids of the variable as an interger array
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The variable id
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @retval ierr @copydoc error_return
!>
  integer function inq_vardimid_vid(File,varid,dimids) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer,            intent(in) :: varid
    integer(i4), intent(out)    :: dimids(:)


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg
    integer :: size_dimids
    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem

    iotype = File%iotype
    ierr=PIO_noerr
    
    size_dimids=size(dimids)

    if(ios%async_interface) then
       if( .not. ios%ioproc ) then
          msg=PIO_MSG_INQ_VARDIMID
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(size_dimids,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_vardimid(File%fh,varid,dimids)

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_variable(File%fh,varid,dimids=dimids)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,dimids=dimids)
          endif

          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(dimids,size(dimids),MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",1047)

       end select
    endif
    call check_netcdf(File,ierr,"nf_mod.F90",1051)
    if(ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(dimids,size_dimids,MPI_INTEGER,ios%IOMaster,ios%My_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_vardimid_vid

!>
!! @public 
!! @ingroup PIO_inq_vardimid
!! @brief returns the dimids of the variable as an interger array
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @retval ierr @copydoc error_return
!>
  integer function inq_vardimid_vdesc(File,vardesc,dimids) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in) :: vardesc
    integer(i4), intent(out)    :: dimids(:)


    ierr = pio_inq_vardimid(File, vardesc%varid, dimids)
  end function inq_vardimid_vdesc

!>
!! @public 
!! @ingroup PIO_inq_varnatts
!! @brief Returns the number of attributes associated with a varaible
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable id
!! @param natts : The number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_varnatts_vid(File,varid,natts) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer           , intent(in) :: varid
    integer(i4), intent(out)         :: natts


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg
    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem

    iotype = File%iotype
    ierr=PIO_noerr

    if(ios%async_interface) then
       if( .not. ios%ioproc ) then
          msg=PIO_MSG_INQ_VARNATTS
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(varid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_varnatts(File%fh,varid,natts)

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_variable(File%fh,varid,nAtts=natts)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,nAtts=natts)
          endif

          call MPI_BCAST(natts,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
          call CheckMPIReturn('nf_mod',mpierr)

       case default
          call bad_iotype(iotype,"nf_mod.F90",1136)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1140)
    if(ios%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(natts,1,MPI_INTEGER,ios%IOMaster,ios%My_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_varnatts_vid

!>
!! @public 
!! @ingroup PIO_inq_varnatts
!! @brief Returns the number of attributes associated with a varaible
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param natts : The number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_varnatts_vdesc(File,vardesc,natts) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in)    :: vardesc
    integer(i4), intent(out)         :: natts


    ierr = pio_inq_varnatts(file, vardesc%varid, natts)
  end function inq_varnatts_vdesc

!>
!! @public 
!! @ingroup PIO_inq_dimid
!! @brief Returns the netcdf dimension id for the name.
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the netcdf dimension.
!! @param dimid : The netcdf dimension id.
!! @retval ierr @copydoc error_return
!!
!! Note that we do not want internal error checking for this funtion.
!>
  integer function pio_inq_dimid(File,name,dimid) result(ierr)

    type (File_desc_t), intent(in) :: File
    character(len=*), intent(in)   :: name
    integer, intent(out)           :: dimid        !dimension ID


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg, nlen
    type(iosystem_desc_t), pointer :: ios

    ios => File%iosystem

    iotype = File%iotype
    ierr=PIO_noerr
    dimid=-1
    nlen = len(name)
    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_DIMID
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(nlen,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(name,nlen,MPI_CHARACTER,ios%CompMaster, ios%my_comm , mpierr)
    end if
    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_dimid(File%fh,name(1:nlen),dimid)

       case (pio_iotype_netcdf4p)
             ierr=nf90_inq_dimid(File%fh,name(1:nlen),dimid)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inq_dimid(File%fh,name(1:nlen),dimid)
          endif
          if(.not. ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(dimid,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",1229)

       end select
    endif

    if(Debug .or. Debugasync) print *,"nf_mod.F90",1234,file%fh, &
      name, dimid, ierr
    call check_netcdf(File, ierr,"nf_mod.F90",1236)

    if(ios%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(dimid,1,MPI_INTEGER,ios%IOMaster,ios%My_comm, mpierr)
       if(Debugasync) print *,"nf_mod.F90",1240,dimid,ierr,mpierr
       call CheckMPIReturn('nf_mod',mpierr)
    end if
 
  end function pio_inq_dimid

!>
!! @public 
!! @ingroup PIO_inq_dimname
!! @brief Gets the name of a dimension given its ID
!! @details
!! @param File @copydoc file_desc_t
!! @param dimid : The netcdf dimension id.
!! @param dimname : The name associated with the netcdf dimension id.
!! @retval ierr @copydoc error_return
!>
  integer function pio_inq_dimname(File,dimid,dimname) result(ierr)

    type (File_desc_t), intent(in) :: File
    integer         , intent(in)   :: dimid
    character(len=*), intent(out)  :: dimname        !dimension name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg, ldn
    type(iosystem_desc_t), pointer :: ios
    
    ios => File%iosystem
    iotype = File%iotype
    ierr=PIO_noerr

    ldn = len(dimname)

    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_DIMNAME
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(dimid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       call MPI_BCAST(ldn,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_dimname(File%fh,dimid,dimname(1:ldn))

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_dimension(File%fh,dimid,name=dimname(1:ldn))
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)

          if (ios%io_rank==0) then
             ierr=nf90_inquire_dimension(File%fh,dimid,name=dimname(1:ldn))
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(dimname,ldn,MPI_CHARACTER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",1308)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1312)
    if(ios%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(dimname,ldn,MPI_CHARACTER,ios%IOMaster,ios%My_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function pio_inq_dimname

!>
!! @public 
!! @ingroup PIO_inq_dimlen
!! @brief Returns the extent of a netCDF dimension 
!! @details
!! @param File @copydoc file_desc_t
!! @param dimid : The netcdf dimension.
!! @param dimlen : The extent of the netcdf dimension.
!! @retval ierr @copydoc error_return
!>
  integer function pio_inq_dimlen(File,dimid,dimlen) result(ierr)

    type (File_desc_t), intent(in) :: File
    integer(i4)     , intent(in)   :: dimid
    integer(i4)     , intent(out)  :: dimlen        !dimension name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg
    integer(kind=PIO_OFFSET) :: clen
    type(iosystem_desc_t), pointer :: ios
    
    ios => File%iosystem
    iotype = File%iotype
    ierr=PIO_noerr

    if(ios%async_interface) then
       if(.not. ios%ioproc ) then
          msg=PIO_MSG_INQ_DIMLEN
          if(debugasync) print *,"nf_mod.F90",1351,msg
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call MPI_BCAST(file%fh,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
       end if
       call MPI_BCAST(dimid,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , mpierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)
          ierr=nfmpi_inq_dimlen(File%fh,dimid,clen)
          dimlen = INT(clen,kind=i4)

       case(pio_iotype_netcdf4p)
          ierr=nf90_inquire_dimension(File%fh,dimid,len=dimlen)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_inquire_dimension(File%fh,dimid,len=dimlen)
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(dimlen,1,MPI_INTEGER,0,ios%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

       case default
          call bad_iotype(iotype,"nf_mod.F90",1381)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1385)
    if(file%iosystem%async_interface .or. ios%num_tasks>ios%num_iotasks) then
       call MPI_BCAST(dimlen,1,MPI_INTEGER,ios%IOMaster,ios%My_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if


  end function pio_inq_dimlen

!> 
!! @public
!! @ingroup PIO_enddef
!! @brief Exits netcdf define mode.
!! @details
!! @param File @copydoc file_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function PIO_enddef(File) result(ierr)
    type (File_desc_t), intent(inout) :: File
    type (iosystem_desc_t), pointer :: ios

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    logical, parameter :: Check = .TRUE.
    integer :: msg = PIO_MSG_ENDDEF
    iotype = File%iotype

    ierr=PIO_noerr

    ios => file%iosystem

    if(ios%async_interface .and. .not. ios%ioproc) then
       if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
    end if
    if(ios%IOproc) then
       select case(iotype)
       case(pio_iotype_pnetcdf)
          ierr=nfmpi_enddef(File%fh)

       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_enddef(File%fh)
          endif
       case(PIO_iotype_netcdf4p)
          ierr=nf90_enddef(File%fh)

       case default
          call bad_iotype(iotype,"nf_mod.F90",1452)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1456)
  end function PIO_enddef

!> 
!! @public
!! @ingroup PIO_redef
!! @brief Re-enters netcdf define mode.   
!! @details 
!! @warning Entering and leaving netcdf define mode causes a file sync operation to 
!!          occur, these operations can be very expensive in parallel systems.   We 
!!          recommend structuring your code to minimize calls to this function.
!! @param File @copydoc file_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function PIO_redef(File) result(ierr)
    type (File_desc_t), intent(inout) :: File

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, msg
    logical, parameter :: Check = .TRUE.
    type(iosystem_desc_t), pointer :: ios
    

    iotype = File%iotype
    ios => file%iosystem
    ierr=PIO_noerr
    if(ios%async_interface .and. .not. ios%ioproc) then
       msg = PIO_MSG_REDEF
       if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
    end if

    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)

          ierr=nfmpi_redef(File%fh)

       case(pio_iotype_netcdf4p)
             ierr=nf90_redef(File%fh)
       case(pio_iotype_netcdf, pio_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_redef(File%fh)
          endif


       case default
          call bad_iotype(iotype,"nf_mod.F90",1510)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1514)
  end function PIO_redef

!> 
!! @public
!! @ingroup PIO_def_dim
!! @brief Defines the netcdf dimension
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the dimension to define
!! @param len :  The size of the dimension
!! @param dimid : The dimension identifier
!<
  integer function PIO_def_dim(File,name,len,dimid) result(ierr)

    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer(i4), intent(in)         :: len
    integer(i4), intent(out)        :: dimid

    !------------------
    ! Local variables
    !------------------
    type(iosystem_desc_t), pointer :: ios
    integer :: iotype, mpierr, nlen
    integer(kind=PIO_Offset)  :: clen
    integer :: msg = PIO_MSG_DEF_DIM

    iotype = File%iotype

    ierr=PIO_noerr
    ios => file%iosystem
    nlen = len_trim(name)
    if(ios%async_interface) then
       if(Debugasync) print *,"nf_mod.F90",1548
       if( .not. ios%ioproc) then
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       end if
       call mpi_bcast(len, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       call mpi_bcast(nlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       call mpi_bcast(name, nlen, mpi_character, ios%compmaster, ios%intercomm, ierr)
       if(Debugasync) print *,"nf_mod.F90",1556,file%fh, name(1:nlen)
    end if
       
    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)

          clen = len
          ierr=nfmpi_def_dim(File%fh,name(1:nlen),clen,dimid)

       case(PIO_iotype_netcdf4p)
          ierr=nf90_def_dim(ncid=File%fh,name=name(1:nlen),len=len,dimid=dimid)
       case(pio_iotype_netcdf,PIO_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr=nf90_def_dim(ncid=File%fh,name=name(1:nlen),len=len,dimid=dimid)
          endif
          if(.not.ios%async_interface .and. ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(dimid, 1, MPI_INTEGER, 0, ios%IO_Comm, ierr)
          end if
       case default
          call bad_iotype(iotype,"nf_mod.F90",1581)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1585)

    if(ios%async_interface .or. ios%num_tasks > ios%num_iotasks) then
       call MPI_BCAST(dimid, 1, MPI_INTEGER, ios%IOMaster, ios%my_Comm, ierr)
    end if
    if(debugasync) print *,"nf_mod.F90",1590,dimid
  end function PIO_def_dim    


!> 
!! @public 
!! @ingroup PIO_def_var
!! @brief Defines a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the variable to define
!! @param type : The type of variable 
!! @param vardesc @copydoc var_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function def_var_0d(File,name,type,vardesc) result(ierr)

    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer, intent(in)             :: type
    type (Var_desc_t), intent(inout) :: vardesc
    integer :: dimids(0)

    ierr = def_var_md(File,name,type,dimids,vardesc)

  end function def_var_0d

!> 
!! @public
!! @ingroup PIO_def_var
!! @brief Defines the a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the variable to define
!! @param type : The type of variable 
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @param vardesc @copydoc var_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function def_var_md(File,name,type,dimids,vardesc) result(ierr)
    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer, intent(in)             :: type
    integer, intent(in)             :: dimids(:)

    type (Var_desc_t), intent(inout) :: vardesc
    type(iosystem_desc_t), pointer :: ios
    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr, nlen
    integer :: msg = PIO_MSG_DEF_VAR


    iotype = File%iotype

    ierr=PIO_noerr
    vardesc%rec=-1
    vardesc%ndims = SIZE(dimids)

    vardesc%type = type

    ios => file%iosystem
    nlen = len_trim(name)

    if(ios%async_interface) then
       if( .not. ios%ioproc) then
          if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
          call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       end if
       call mpi_bcast(type, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       
       call mpi_bcast(nlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       call mpi_bcast(name, nlen, mpi_character, ios%compmaster, ios%intercomm, ierr)
       call mpi_bcast(vardesc%ndims, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       call mpi_bcast(dimids, vardesc%ndims, mpi_integer, ios%compmaster, ios%intercomm, ierr)
    endif
    if(ios%IOproc) then
       select case(iotype)
       case(pio_iotype_pnetcdf)
          if(vardesc%ndims==0) then
             ierr=nfmpi_def_var(File%fh,name(1:nlen),type,vardesc%ndims,dimids,vardesc%varid)
          else
             ierr=nfmpi_def_var(File%fh,name(1:nlen),type,vardesc%ndims,dimids(1:vardesc%ndims),vardesc%varid)
          end if

       case(pio_iotype_netcdf4p)
          if(vardesc%ndims==0) then
             ierr=nf90_def_var( ncid=File%fh,name=name(1:nlen),xtype=type, &
                  varid=vardesc%varid)
          else
             ierr=nf90_def_var( ncid=File%fh,name=name(1:nlen),xtype=type, &
                  dimids=dimids(1:vardesc%ndims),varid=vardesc%varid)
          endif
          ierr = nf90_def_var_fill(File%fh, vardesc%varid, 1, 0)
       case(pio_iotype_netcdf,pio_iotype_netcdf4c)
          ! assuming type valid for both pnetcdf and netcdf
          if (ios%io_rank==0) then
             if(vardesc%ndims==0) then
                ierr=nf90_def_var( ncid=File%fh,name=name(1:nlen),xtype=type, &
                     varid=vardesc%varid)
             else
                ierr=nf90_def_var( ncid=File%fh,name=name(1:nlen),xtype=type, &
                     dimids=dimids(1:vardesc%ndims),varid=vardesc%varid)
             end if
             if (Debug) print *, '0: def_var fh=',File%fh, &
                  'name=',name(1:nlen),' id=',vardesc%varid
             if(iotype==pio_iotype_netcdf4c) then
                if(vardesc%ndims>0 .and. ierr==PIO_NOERR) then
                   ierr = nf90_def_var_deflate(File%fh,vardesc%varid,0,1,1)
                end if
             endif

          endif
          if(.not.ios%async_interface.and.ios%num_tasks==ios%num_iotasks) then
             call MPI_BCAST(vardesc%varid, 1, MPI_INTEGER, 0, ios%IO_Comm, ierr)
          end if
       case default
          call bad_iotype(iotype,"nf_mod.F90",1727)

       end select
    endif
    call check_netcdf(File, ierr,"nf_mod.F90",1731)
    if(ios%async_interface  .or. ios%num_tasks> ios%num_iotasks) then  
       call MPI_BCAST(vardesc%varid, 1, MPI_INTEGER, ios%Iomaster, ios%my_Comm, ierr)
    end if
  end function def_var_md

!>
!! @public
!! @ingroup PIO_copy_att
!! @brief No idea what this function does
!! @details 
!! @param infile @copydoc file_desc_t
!! @param invarid :
!! @param name : 
!! @param outfile :
!! @param outvarid :
!! @retval ierr @copydoc error_return
!<
  integer function pio_copy_att(infile, invarid, name, outfile, outvarid) result(ierr)

    type (File_desc_t), intent(in)  :: infile, outfile
    character(len=*), intent(in)    :: name
    integer, intent(in) :: invarid, outvarid
    integer :: iotype, mpierr, msg
    type(iosystem_desc_t), pointer :: ios


    ios => infile%iosystem
    ierr=PIO_noerr
    iotype = infile%iotype
    if(ios%IOproc) then
       select case(iotype)

       case(pio_iotype_pnetcdf)

          ierr = nfmpi_copy_att(infile%fh, invarid, name, &
               outfile%fh, outvarid)
       case(pio_iotype_netcdf,PIO_iotype_netcdf4c)
          if (ios%io_rank==0) then
             ierr = nf90_copy_att(infile%fh,invarid,name,&
                  outfile%fh,outvarid)     
          end if
       case(PIO_iotype_netcdf4p)
          ierr = nf90_copy_att(infile%fh,invarid,name,&
               outfile%fh,outvarid)     
       end select
    end if
    call check_netcdf(outFile, ierr,"nf_mod.F90",1782)
  end function pio_copy_att


!>
!! @public 
!! @ingroup PIO_inquire_variable
!! @brief Inquires if a NetCDF variable is present and returns its attributes  
!! @details
!! @param ncid : A netcdf file descriptor returned by \ref PIO_openfile or \ref PIO_createfile.
!! @param varid : The netcdf variable ID.
!! @param name : The name of the variable
!! @param xtype : The type of the variable
!! @param ndims : The number of dimensions for the variable.
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @param natts : Number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inquire_variable_vid(ncid, varid, name, xtype, ndims, dimids, natts) result(ierr)
    type(file_desc_t), intent(in) :: ncid
    integer,                         intent( in) :: varid
    character (len = *),   optional, intent(out) :: name
    integer,               optional, intent(out) :: xtype, ndims
    integer, dimension(:), optional, intent(out) :: dimids
    integer,               optional, intent(out) :: natts

    
    if(present(name)) ierr = pio_inq_varname(ncid, varid, name)
    if(present(ndims)) ierr = pio_inq_varndims(ncid, varid, ndims)
    if(present(dimids)) ierr = pio_inq_vardimid(ncid, varid, dimids)
    if(present(natts)) ierr = pio_inq_varnatts(ncid, varid, natts)
    if(present(xtype)) ierr = pio_inq_vartype(ncid, varid, xtype)



  end function inquire_variable_vid

!>
!! @public 
!! @ingroup PIO_inquire_variable
!! @brief Inquires if a NetCDF variable is present and returns its attributes  
!! @details
!! @param ncid : A netcdf file descriptor returned by \ref PIO_openfile or \ref PIO_createfile.
!! @param vardesc @copydoc var_desc_t
!! @param name : The name of the variable
!! @param xtype : The type of the variable
!! @param ndims : The number of dimensions for the variable.
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @param natts : Number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inquire_variable_vdesc(ncid, vardesc, name, xtype, ndims, dimids, natts) result(ierr)
    type(file_desc_t),               intent(in) :: ncid
    type(var_desc_t),                intent( in) :: vardesc
    character (len = *),   optional, intent(out) :: name
    integer,               optional, intent(out) :: xtype, ndims
    integer, dimension(:), optional, intent(out) :: dimids
    integer,               optional, intent(out) :: natts

    if(present(name)) ierr = pio_inq_varname(ncid, vardesc, name)
    if(present(ndims)) ierr = pio_inq_varndims(ncid, vardesc, ndims)
    if(present(dimids)) ierr = pio_inq_vardimid(ncid, vardesc, dimids)
    if(present(natts)) ierr = pio_inq_varnatts(ncid, vardesc, natts)
    if(present(xtype)) ierr = pio_inq_vartype(ncid, vardesc, xtype)

  end function inquire_variable_vdesc

!>
!! @public 
!! @ingroup PIO_inquire_dimension
!! @brief  Get information about a particular dimension in netcdf file 
!! @details
!! @param ncid : A netcdf file descriptor returned by \ref PIO_openfile or \ref PIO_createfile.
!! @param dimid : The netcdf dimension ID.
!! @param name : The name of the dimension.
!! @param len : The length of the dimesions name.
!! @retval ierr @copydoc error_return
!>
  integer function PIO_inquire_dimension(ncid, dimid, name, len) result(ierr)
    type(file_desc_T),             intent(in)  :: ncid
    integer,                       intent( in) :: dimid
    character (len = *), optional, intent(out) :: name
    integer,             optional, intent(out) :: len

    if(present(len)) ierr = pio_inq_dimlen(ncid, dimid, len)
    if(present(name)) ierr = pio_inq_dimname(ncid, dimid,name)

  end function PIO_inquire_dimension


end module nf_mod
