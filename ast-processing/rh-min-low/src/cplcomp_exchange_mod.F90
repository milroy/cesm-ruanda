module cplcomp_exchange_mod

  use shr_kind_mod, only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod, only: CL => SHR_KIND_CL, CX => SHR_KIND_CX, CXX => SHR_KIND_CXX
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod,  only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use mct_mod
  use seq_map_type_mod
  use component_type_mod
  use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other 
  use seq_comm_mct, only: cplid, logunit
  use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs, seq_comm_iamin 
  use seq_diag_mct

  implicit none
  private  ! except
!      
!      
!      (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!      
!      DO NOT EDIT
!      This file created by buildiface 
!      
       INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
       PARAMETER (MPI_SOURCE=3,MPI_TAG=4,MPI_ERROR=5)
       INTEGER MPI_STATUS_SIZE
       PARAMETER (MPI_STATUS_SIZE=5)
       INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
       INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
       INTEGER MPI_ERRCODES_IGNORE(1)
       CHARACTER*1 MPI_ARGVS_NULL(1,1)
       CHARACTER*1 MPI_ARGV_NULL(1)
       INTEGER MPI_SUCCESS
       PARAMETER (MPI_SUCCESS=0)
       INTEGER MPI_ERR_OTHER
       PARAMETER (MPI_ERR_OTHER=15)
       INTEGER MPI_ERR_COUNT
       PARAMETER (MPI_ERR_COUNT=2)
       INTEGER MPI_ERR_SPAWN
       PARAMETER (MPI_ERR_SPAWN=42)
       INTEGER MPI_ERR_LOCKTYPE
       PARAMETER (MPI_ERR_LOCKTYPE=47)
       INTEGER MPI_ERR_OP
       PARAMETER (MPI_ERR_OP=9)
       INTEGER MPI_ERR_DUP_DATAREP
       PARAMETER (MPI_ERR_DUP_DATAREP=24)
       INTEGER MPI_ERR_UNSUPPORTED_DATAREP
       PARAMETER (MPI_ERR_UNSUPPORTED_DATAREP=43)
       INTEGER MPI_ERR_TRUNCATE
       PARAMETER (MPI_ERR_TRUNCATE=14)
       INTEGER MPI_ERR_INFO_NOKEY
       PARAMETER (MPI_ERR_INFO_NOKEY=31)
       INTEGER MPI_ERR_ASSERT
       PARAMETER (MPI_ERR_ASSERT=53)
       INTEGER MPI_ERR_FILE_EXISTS
       PARAMETER (MPI_ERR_FILE_EXISTS=25)
       INTEGER MPI_ERR_PENDING
       PARAMETER (MPI_ERR_PENDING=18)
       INTEGER MPI_ERR_COMM
       PARAMETER (MPI_ERR_COMM=5)
       INTEGER MPI_ERR_KEYVAL
       PARAMETER (MPI_ERR_KEYVAL=48)
       INTEGER MPI_ERR_NAME
       PARAMETER (MPI_ERR_NAME=33)
       INTEGER MPI_ERR_REQUEST
       PARAMETER (MPI_ERR_REQUEST=19)
       INTEGER MPI_ERR_TYPE
       PARAMETER (MPI_ERR_TYPE=3)
       INTEGER MPI_ERR_INFO_VALUE
       PARAMETER (MPI_ERR_INFO_VALUE=30)
       INTEGER MPI_ERR_RMA_SYNC
       PARAMETER (MPI_ERR_RMA_SYNC=50)
       INTEGER MPI_ERR_NO_MEM
       PARAMETER (MPI_ERR_NO_MEM=34)
       INTEGER MPI_ERR_BAD_FILE
       PARAMETER (MPI_ERR_BAD_FILE=22)
       INTEGER MPI_ERR_QUOTA
       PARAMETER (MPI_ERR_QUOTA=39)
       INTEGER MPI_ERR_ROOT
       PARAMETER (MPI_ERR_ROOT=7)
       INTEGER MPI_ERR_SERVICE
       PARAMETER (MPI_ERR_SERVICE=41)
       INTEGER MPI_ERR_IO
       PARAMETER (MPI_ERR_IO=32)
       INTEGER MPI_ERR_ACCESS
       PARAMETER (MPI_ERR_ACCESS=20)
       INTEGER MPI_ERR_NO_SPACE
       PARAMETER (MPI_ERR_NO_SPACE=36)
       INTEGER MPI_ERR_CONVERSION
       PARAMETER (MPI_ERR_CONVERSION=23)
       INTEGER MPI_ERR_WIN
       PARAMETER (MPI_ERR_WIN=45)
       INTEGER MPI_ERR_FILE
       PARAMETER (MPI_ERR_FILE=27)
       INTEGER MPI_ERR_BASE
       PARAMETER (MPI_ERR_BASE=46)
       INTEGER MPI_ERR_IN_STATUS
       PARAMETER (MPI_ERR_IN_STATUS=17)
       INTEGER MPI_ERR_RMA_CONFLICT
       PARAMETER (MPI_ERR_RMA_CONFLICT=49)
       INTEGER MPI_ERR_INFO_KEY
       PARAMETER (MPI_ERR_INFO_KEY=29)
       INTEGER MPI_ERR_ARG
       PARAMETER (MPI_ERR_ARG=12)
       INTEGER MPI_ERR_READ_ONLY
       PARAMETER (MPI_ERR_READ_ONLY=40)
       INTEGER MPI_ERR_SIZE
       PARAMETER (MPI_ERR_SIZE=51)
       INTEGER MPI_ERR_BUFFER
       PARAMETER (MPI_ERR_BUFFER=1)
       INTEGER MPI_ERR_LASTCODE
       PARAMETER (MPI_ERR_LASTCODE=1073741823)
       INTEGER MPI_ERR_DISP
       PARAMETER (MPI_ERR_DISP=52)
       INTEGER MPI_ERR_PORT
       PARAMETER (MPI_ERR_PORT=38)
       INTEGER MPI_ERR_GROUP
       PARAMETER (MPI_ERR_GROUP=8)
       INTEGER MPI_ERR_TOPOLOGY
       PARAMETER (MPI_ERR_TOPOLOGY=10)
       INTEGER MPI_ERR_TAG
       PARAMETER (MPI_ERR_TAG=4)
       INTEGER MPI_ERR_NOT_SAME
       PARAMETER (MPI_ERR_NOT_SAME=35)
       INTEGER MPI_ERR_INFO
       PARAMETER (MPI_ERR_INFO=28)
       INTEGER MPI_ERR_UNKNOWN
       PARAMETER (MPI_ERR_UNKNOWN=13)
       INTEGER MPI_ERR_FILE_IN_USE
       PARAMETER (MPI_ERR_FILE_IN_USE=26)
       INTEGER MPI_ERR_UNSUPPORTED_OPERATION
       PARAMETER (MPI_ERR_UNSUPPORTED_OPERATION=44)
       INTEGER MPI_ERR_AMODE
       PARAMETER (MPI_ERR_AMODE=21)
       INTEGER MPI_ERR_RANK
       PARAMETER (MPI_ERR_RANK=6)
       INTEGER MPI_ERR_DIMS
       PARAMETER (MPI_ERR_DIMS=11)
       INTEGER MPI_ERR_NO_SUCH_FILE
       PARAMETER (MPI_ERR_NO_SUCH_FILE=37)
       INTEGER MPI_ERR_INTERN
       PARAMETER (MPI_ERR_INTERN=16)
       INTEGER MPI_ERRORS_ARE_FATAL
       PARAMETER (MPI_ERRORS_ARE_FATAL=1409286144)
       INTEGER MPI_ERRORS_RETURN
       PARAMETER (MPI_ERRORS_RETURN=1409286145)
       INTEGER MPI_IDENT
       PARAMETER (MPI_IDENT=0)
       INTEGER MPI_CONGRUENT
       PARAMETER (MPI_CONGRUENT=1)
       INTEGER MPI_SIMILAR
       PARAMETER (MPI_SIMILAR=2)
       INTEGER MPI_UNEQUAL
       PARAMETER (MPI_UNEQUAL=3)
       INTEGER MPI_MAX
       PARAMETER (MPI_MAX=1476395009)
       INTEGER MPI_MIN
       PARAMETER (MPI_MIN=1476395010)
       INTEGER MPI_SUM
       PARAMETER (MPI_SUM=1476395011)
       INTEGER MPI_PROD
       PARAMETER (MPI_PROD=1476395012)
       INTEGER MPI_LAND
       PARAMETER (MPI_LAND=1476395013)
       INTEGER MPI_BAND
       PARAMETER (MPI_BAND=1476395014)
       INTEGER MPI_LOR
       PARAMETER (MPI_LOR=1476395015)
       INTEGER MPI_BOR
       PARAMETER (MPI_BOR=1476395016)
       INTEGER MPI_LXOR
       PARAMETER (MPI_LXOR=1476395017)
       INTEGER MPI_BXOR
       PARAMETER (MPI_BXOR=1476395018)
       INTEGER MPI_MINLOC
       PARAMETER (MPI_MINLOC=1476395019)
       INTEGER MPI_MAXLOC
       PARAMETER (MPI_MAXLOC=1476395020)
       INTEGER MPI_REPLACE
       PARAMETER (MPI_REPLACE=1476395021)
       INTEGER MPI_COMM_WORLD
       PARAMETER (MPI_COMM_WORLD=1140850688)
       INTEGER MPI_COMM_SELF
       PARAMETER (MPI_COMM_SELF=1140850689)
       INTEGER MPI_GROUP_EMPTY
       PARAMETER (MPI_GROUP_EMPTY=1207959552)
       INTEGER MPI_COMM_NULL
       PARAMETER (MPI_COMM_NULL=67108864)
       INTEGER MPI_WIN_NULL
       PARAMETER (MPI_WIN_NULL=536870912)
       INTEGER MPI_FILE_NULL
       PARAMETER (MPI_FILE_NULL=0)
       INTEGER MPI_GROUP_NULL
       PARAMETER (MPI_GROUP_NULL=134217728)
       INTEGER MPI_OP_NULL
       PARAMETER (MPI_OP_NULL=402653184)
       INTEGER MPI_DATATYPE_NULL
       PARAMETER (MPI_DATATYPE_NULL=201326592)
       INTEGER MPI_REQUEST_NULL
       PARAMETER (MPI_REQUEST_NULL=738197504)
       INTEGER MPI_ERRHANDLER_NULL
       PARAMETER (MPI_ERRHANDLER_NULL=335544320)
       INTEGER MPI_INFO_NULL
       PARAMETER (MPI_INFO_NULL=469762048)
       INTEGER MPI_TAG_UB
       PARAMETER (MPI_TAG_UB=1681915906)
       INTEGER MPI_HOST
       PARAMETER (MPI_HOST=1681915908)
       INTEGER MPI_IO
       PARAMETER (MPI_IO=1681915910)
       INTEGER MPI_WTIME_IS_GLOBAL
       PARAMETER (MPI_WTIME_IS_GLOBAL=1681915912)
       INTEGER MPI_UNIVERSE_SIZE
       PARAMETER (MPI_UNIVERSE_SIZE=1681915914)
       INTEGER MPI_LASTUSEDCODE
       PARAMETER (MPI_LASTUSEDCODE=1681915916)
       INTEGER MPI_APPNUM
       PARAMETER (MPI_APPNUM=1681915918)
       INTEGER MPI_WIN_BASE
       PARAMETER (MPI_WIN_BASE=1711276034)
       INTEGER MPI_WIN_SIZE
       PARAMETER (MPI_WIN_SIZE=1711276036)
       INTEGER MPI_WIN_DISP_UNIT
       PARAMETER (MPI_WIN_DISP_UNIT=1711276038)
       INTEGER MPI_MAX_ERROR_STRING
       PARAMETER (MPI_MAX_ERROR_STRING=512-1)
       INTEGER MPI_MAX_PORT_NAME
       PARAMETER (MPI_MAX_PORT_NAME=255)
       INTEGER MPI_MAX_OBJECT_NAME
       PARAMETER (MPI_MAX_OBJECT_NAME=127)
       INTEGER MPI_MAX_INFO_KEY
       PARAMETER (MPI_MAX_INFO_KEY=254)
       INTEGER MPI_MAX_INFO_VAL
       PARAMETER (MPI_MAX_INFO_VAL=1023)
       INTEGER MPI_MAX_PROCESSOR_NAME
       PARAMETER (MPI_MAX_PROCESSOR_NAME=128-1)
       INTEGER MPI_MAX_DATAREP_STRING
       PARAMETER (MPI_MAX_DATAREP_STRING=127)
       INTEGER MPI_UNDEFINED
       PARAMETER (MPI_UNDEFINED=(-32766))
       INTEGER MPI_KEYVAL_INVALID
       PARAMETER (MPI_KEYVAL_INVALID=603979776)
       INTEGER MPI_BSEND_OVERHEAD
       PARAMETER (MPI_BSEND_OVERHEAD=88)
       INTEGER MPI_PROC_NULL
       PARAMETER (MPI_PROC_NULL=-1)
       INTEGER MPI_ANY_SOURCE
       PARAMETER (MPI_ANY_SOURCE=-2)
       INTEGER MPI_ANY_TAG
       PARAMETER (MPI_ANY_TAG=-1)
       INTEGER MPI_ROOT
       PARAMETER (MPI_ROOT=-3)
       INTEGER MPI_GRAPH
       PARAMETER (MPI_GRAPH=1)
       INTEGER MPI_CART
       PARAMETER (MPI_CART=2)
       INTEGER MPI_DIST_GRAPH
       PARAMETER (MPI_DIST_GRAPH=3)
       INTEGER MPI_VERSION
       PARAMETER (MPI_VERSION=2)
       INTEGER MPI_SUBVERSION
       PARAMETER (MPI_SUBVERSION=2)
       INTEGER MPI_LOCK_EXCLUSIVE
       PARAMETER (MPI_LOCK_EXCLUSIVE=234)
       INTEGER MPI_LOCK_SHARED
       PARAMETER (MPI_LOCK_SHARED=235)
       INTEGER MPI_COMPLEX
       PARAMETER (MPI_COMPLEX=1275070494)
       INTEGER MPI_DOUBLE_COMPLEX
       PARAMETER (MPI_DOUBLE_COMPLEX=1275072546)
       INTEGER MPI_LOGICAL
       PARAMETER (MPI_LOGICAL=1275069469)
       INTEGER MPI_REAL
       PARAMETER (MPI_REAL=1275069468)
       INTEGER MPI_DOUBLE_PRECISION
       PARAMETER (MPI_DOUBLE_PRECISION=1275070495)
       INTEGER MPI_INTEGER
       PARAMETER (MPI_INTEGER=1275069467)
       INTEGER MPI_2INTEGER
       PARAMETER (MPI_2INTEGER=1275070496)
       INTEGER MPI_2COMPLEX
       PARAMETER (MPI_2COMPLEX=1275072548)
       INTEGER MPI_2DOUBLE_PRECISION
       PARAMETER (MPI_2DOUBLE_PRECISION=1275072547)
       INTEGER MPI_2REAL
       PARAMETER (MPI_2REAL=1275070497)
       INTEGER MPI_2DOUBLE_COMPLEX
       PARAMETER (MPI_2DOUBLE_COMPLEX=1275076645)
       INTEGER MPI_CHARACTER
       PARAMETER (MPI_CHARACTER=1275068698)
       INTEGER MPI_BYTE
       PARAMETER (MPI_BYTE=1275068685)
       INTEGER MPI_UB
       PARAMETER (MPI_UB=1275068433)
       INTEGER MPI_LB
       PARAMETER (MPI_LB=1275068432)
       INTEGER MPI_PACKED
       PARAMETER (MPI_PACKED=1275068687)
       INTEGER MPI_INTEGER1
       PARAMETER (MPI_INTEGER1=1275068717)
       INTEGER MPI_INTEGER2
       PARAMETER (MPI_INTEGER2=1275068975)
       INTEGER MPI_INTEGER4
       PARAMETER (MPI_INTEGER4=1275069488)
       INTEGER MPI_INTEGER8
       PARAMETER (MPI_INTEGER8=1275070513)
       INTEGER MPI_INTEGER16
       PARAMETER (MPI_INTEGER16=MPI_DATATYPE_NULL)
       INTEGER MPI_REAL4
       PARAMETER (MPI_REAL4=1275069479)
       INTEGER MPI_REAL8
       PARAMETER (MPI_REAL8=1275070505)
       INTEGER MPI_REAL16
       PARAMETER (MPI_REAL16=1275072555)
       INTEGER MPI_COMPLEX8
       PARAMETER (MPI_COMPLEX8=1275070504)
       INTEGER MPI_COMPLEX16
       PARAMETER (MPI_COMPLEX16=1275072554)
       INTEGER MPI_COMPLEX32
       PARAMETER (MPI_COMPLEX32=1275076652)
       INTEGER MPI_ADDRESS_KIND, MPI_OFFSET_KIND, MPI_INTEGER_KIND
       PARAMETER (MPI_ADDRESS_KIND=8)
       PARAMETER (MPI_OFFSET_KIND=8)
       PARAMETER (MPI_INTEGER_KIND=8)
       INTEGER MPI_CHAR
       PARAMETER (MPI_CHAR=1275068673)
       INTEGER MPI_SIGNED_CHAR
       PARAMETER (MPI_SIGNED_CHAR=1275068696)
       INTEGER MPI_UNSIGNED_CHAR
       PARAMETER (MPI_UNSIGNED_CHAR=1275068674)
       INTEGER MPI_WCHAR
       PARAMETER (MPI_WCHAR=1275069454)
       INTEGER MPI_SHORT
       PARAMETER (MPI_SHORT=1275068931)
       INTEGER MPI_UNSIGNED_SHORT
       PARAMETER (MPI_UNSIGNED_SHORT=1275068932)
       INTEGER MPI_INT
       PARAMETER (MPI_INT=1275069445)
       INTEGER MPI_UNSIGNED
       PARAMETER (MPI_UNSIGNED=1275069446)
       INTEGER MPI_LONG
       PARAMETER (MPI_LONG=1275070471)
       INTEGER MPI_UNSIGNED_LONG
       PARAMETER (MPI_UNSIGNED_LONG=1275070472)
       INTEGER MPI_FLOAT
       PARAMETER (MPI_FLOAT=1275069450)
       INTEGER MPI_DOUBLE
       PARAMETER (MPI_DOUBLE=1275070475)
       INTEGER MPI_LONG_DOUBLE
       PARAMETER (MPI_LONG_DOUBLE=1275072524)
       INTEGER MPI_LONG_LONG_INT
       PARAMETER (MPI_LONG_LONG_INT=1275070473)
       INTEGER MPI_UNSIGNED_LONG_LONG
       PARAMETER (MPI_UNSIGNED_LONG_LONG=1275070489)
       INTEGER MPI_LONG_LONG
       PARAMETER (MPI_LONG_LONG=1275070473)
       INTEGER MPI_FLOAT_INT
       PARAMETER (MPI_FLOAT_INT=-1946157056)
       INTEGER MPI_DOUBLE_INT
       PARAMETER (MPI_DOUBLE_INT=-1946157055)
       INTEGER MPI_LONG_INT
       PARAMETER (MPI_LONG_INT=-1946157054)
       INTEGER MPI_SHORT_INT
       PARAMETER (MPI_SHORT_INT=-1946157053)
       INTEGER MPI_2INT
       PARAMETER (MPI_2INT=1275070486)
       INTEGER MPI_LONG_DOUBLE_INT
       PARAMETER (MPI_LONG_DOUBLE_INT=-1946157052)
       INTEGER MPI_INT8_T
       PARAMETER (MPI_INT8_T=1275068727)
       INTEGER MPI_INT16_T
       PARAMETER (MPI_INT16_T=1275068984)
       INTEGER MPI_INT32_T
       PARAMETER (MPI_INT32_T=1275069497)
       INTEGER MPI_INT64_T
       PARAMETER (MPI_INT64_T=1275070522)
       INTEGER MPI_UINT8_T
       PARAMETER (MPI_UINT8_T=1275068731)
       INTEGER MPI_UINT16_T
       PARAMETER (MPI_UINT16_T=1275068988)
       INTEGER MPI_UINT32_T
       PARAMETER (MPI_UINT32_T=1275069501)
       INTEGER MPI_UINT64_T
       PARAMETER (MPI_UINT64_T=1275070526)
       INTEGER MPI_C_BOOL
       PARAMETER (MPI_C_BOOL=1275068735)
       INTEGER MPI_C_FLOAT_COMPLEX
       PARAMETER (MPI_C_FLOAT_COMPLEX=1275070528)
       INTEGER MPI_C_COMPLEX
       PARAMETER (MPI_C_COMPLEX=1275070528)
       INTEGER MPI_C_DOUBLE_COMPLEX
       PARAMETER (MPI_C_DOUBLE_COMPLEX=1275072577)
       INTEGER MPI_C_LONG_DOUBLE_COMPLEX
       PARAMETER (MPI_C_LONG_DOUBLE_COMPLEX=1275076674)
       INTEGER MPI_AINT
       PARAMETER (MPI_AINT=1275070531)
       INTEGER MPI_OFFSET
       PARAMETER (MPI_OFFSET=1275070532)
       INTEGER MPI_COMBINER_NAMED
       PARAMETER (MPI_COMBINER_NAMED=1)
       INTEGER MPI_COMBINER_DUP
       PARAMETER (MPI_COMBINER_DUP=2)
       INTEGER MPI_COMBINER_CONTIGUOUS
       PARAMETER (MPI_COMBINER_CONTIGUOUS=3)
       INTEGER MPI_COMBINER_VECTOR
       PARAMETER (MPI_COMBINER_VECTOR=4)
       INTEGER MPI_COMBINER_HVECTOR_INTEGER
       PARAMETER (MPI_COMBINER_HVECTOR_INTEGER=5)
       INTEGER MPI_COMBINER_HVECTOR
       PARAMETER (MPI_COMBINER_HVECTOR=6)
       INTEGER MPI_COMBINER_INDEXED
       PARAMETER (MPI_COMBINER_INDEXED=7)
       INTEGER MPI_COMBINER_HINDEXED_INTEGER
       PARAMETER (MPI_COMBINER_HINDEXED_INTEGER=8)
       INTEGER MPI_COMBINER_HINDEXED
       PARAMETER (MPI_COMBINER_HINDEXED=9)
       INTEGER MPI_COMBINER_INDEXED_BLOCK
       PARAMETER (MPI_COMBINER_INDEXED_BLOCK=10)
       INTEGER MPI_COMBINER_STRUCT_INTEGER
       PARAMETER (MPI_COMBINER_STRUCT_INTEGER=12)
       INTEGER MPI_COMBINER_STRUCT
       PARAMETER (MPI_COMBINER_STRUCT=13)
       INTEGER MPI_COMBINER_SUBARRAY
       PARAMETER (MPI_COMBINER_SUBARRAY=14)
       INTEGER MPI_COMBINER_DARRAY
       PARAMETER (MPI_COMBINER_DARRAY=15)
       INTEGER MPI_COMBINER_F90_REAL
       PARAMETER (MPI_COMBINER_F90_REAL=16)
       INTEGER MPI_COMBINER_F90_COMPLEX
       PARAMETER (MPI_COMBINER_F90_COMPLEX=17)
       INTEGER MPI_COMBINER_F90_INTEGER
       PARAMETER (MPI_COMBINER_F90_INTEGER=18)
       INTEGER MPI_COMBINER_RESIZED
       PARAMETER (MPI_COMBINER_RESIZED=19)
       INTEGER MPIX_COMBINER_HINDEXED_BLOCK
       PARAMETER (MPIX_COMBINER_HINDEXED_BLOCK=11)
       INTEGER MPI_TYPECLASS_REAL
       PARAMETER (MPI_TYPECLASS_REAL=1)
       INTEGER MPI_TYPECLASS_INTEGER
       PARAMETER (MPI_TYPECLASS_INTEGER=2)
       INTEGER MPI_TYPECLASS_COMPLEX
       PARAMETER (MPI_TYPECLASS_COMPLEX=3)
       INTEGER MPI_MODE_NOCHECK
       PARAMETER (MPI_MODE_NOCHECK=1024)
       INTEGER MPI_MODE_NOSTORE
       PARAMETER (MPI_MODE_NOSTORE=2048)
       INTEGER MPI_MODE_NOPUT
       PARAMETER (MPI_MODE_NOPUT=4096)
       INTEGER MPI_MODE_NOPRECEDE
       PARAMETER (MPI_MODE_NOPRECEDE=8192)
       INTEGER MPI_MODE_NOSUCCEED
       PARAMETER (MPI_MODE_NOSUCCEED=16384)
       INTEGER MPIX_COMM_TYPE_SHARED
       PARAMETER (MPIX_COMM_TYPE_SHARED=1)
       INTEGER MPIX_MESSAGE_NULL
       PARAMETER (MPIX_MESSAGE_NULL=MPI_REQUEST_NULL)
       INTEGER MPIX_MESSAGE_NO_PROC
       PARAMETER (MPIX_MESSAGE_NO_PROC=1811939328)
       INTEGER MPI_THREAD_SINGLE
       PARAMETER (MPI_THREAD_SINGLE=0)
       INTEGER MPI_THREAD_FUNNELED
       PARAMETER (MPI_THREAD_FUNNELED=1)
       INTEGER MPI_THREAD_SERIALIZED
       PARAMETER (MPI_THREAD_SERIALIZED=2)
       INTEGER MPI_THREAD_MULTIPLE
       PARAMETER (MPI_THREAD_MULTIPLE=3)
       INTEGER MPI_MODE_RDONLY
       PARAMETER (MPI_MODE_RDONLY=2)
       INTEGER MPI_MODE_RDWR
       PARAMETER (MPI_MODE_RDWR=8)
       INTEGER MPI_MODE_WRONLY
       PARAMETER (MPI_MODE_WRONLY=4)
       INTEGER MPI_MODE_DELETE_ON_CLOSE
       PARAMETER (MPI_MODE_DELETE_ON_CLOSE=16)
       INTEGER MPI_MODE_UNIQUE_OPEN
       PARAMETER (MPI_MODE_UNIQUE_OPEN=32)
       INTEGER MPI_MODE_CREATE
       PARAMETER (MPI_MODE_CREATE=1)
       INTEGER MPI_MODE_EXCL
       PARAMETER (MPI_MODE_EXCL=64)
       INTEGER MPI_MODE_APPEND
       PARAMETER (MPI_MODE_APPEND=128)
       INTEGER MPI_MODE_SEQUENTIAL
       PARAMETER (MPI_MODE_SEQUENTIAL=256)
       INTEGER MPI_SEEK_SET
       PARAMETER (MPI_SEEK_SET=600)
       INTEGER MPI_SEEK_CUR
       PARAMETER (MPI_SEEK_CUR=602)
       INTEGER MPI_SEEK_END
       PARAMETER (MPI_SEEK_END=604)
       INTEGER MPI_ORDER_C
       PARAMETER (MPI_ORDER_C=56)
       INTEGER MPI_ORDER_FORTRAN
       PARAMETER (MPI_ORDER_FORTRAN=57)
       INTEGER MPI_DISTRIBUTE_BLOCK
       PARAMETER (MPI_DISTRIBUTE_BLOCK=121)
       INTEGER MPI_DISTRIBUTE_CYCLIC
       PARAMETER (MPI_DISTRIBUTE_CYCLIC=122)
       INTEGER MPI_DISTRIBUTE_NONE
       PARAMETER (MPI_DISTRIBUTE_NONE=123)
       INTEGER MPI_DISTRIBUTE_DFLT_DARG
       PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
       integer*8 MPI_DISPLACEMENT_CURRENT
       PARAMETER (MPI_DISPLACEMENT_CURRENT=-54278278)
       INTEGER MPI_BOTTOM, MPI_IN_PLACE, MPI_UNWEIGHTED
       EXTERNAL MPI_DUP_FN, MPI_NULL_DELETE_FN, MPI_NULL_COPY_FN
       EXTERNAL MPI_WTIME, MPI_WTICK
       EXTERNAL PMPI_WTIME, PMPI_WTICK
       EXTERNAL MPI_COMM_DUP_FN, MPI_COMM_NULL_DELETE_FN
       EXTERNAL MPI_COMM_NULL_COPY_FN
       EXTERNAL MPI_WIN_DUP_FN, MPI_WIN_NULL_DELETE_FN
       EXTERNAL MPI_WIN_NULL_COPY_FN
       EXTERNAL MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN
       EXTERNAL MPI_TYPE_NULL_COPY_FN
       EXTERNAL MPI_CONVERSION_FN_NULL
       REAL*8 MPI_WTIME, MPI_WTICK
       REAL*8 PMPI_WTIME, PMPI_WTICK


       CHARACTER*1 PADS_A(3), PADS_B(3)
       COMMON /MPIFCMB1/ MPI_STATUS_IGNORE
       COMMON /MPIFCMB2/ MPI_STATUSES_IGNORE
       COMMON /MPIFCMB3/ MPI_BOTTOM
       COMMON /MPIFCMB4/ MPI_IN_PLACE
       COMMON /MPIFCMB5/ MPI_UNWEIGHTED
       COMMON /MPIFCMB6/ MPI_ERRCODES_IGNORE
       COMMON /MPIFCMB7/ MPI_ARGVS_NULL, PADS_A
       COMMON /MPIFCMB8/ MPI_ARGV_NULL, PADS_B
       SAVE /MPIFCMB1/,/MPIFCMB2/
       SAVE /MPIFCMB3/,/MPIFCMB4/,/MPIFCMB5/,/MPIFCMB6/
       SAVE /MPIFCMB7/,/MPIFCMB8/
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_map_init_exchange   ! union of cpl/component pes
  public :: seq_map_map_exchange    ! union of cpl/component pes
  public :: seq_mctext_gsmapInit
  public :: seq_mctext_avInit
  public :: seq_mctext_gGridInit
  public :: seq_mctext_avExtend

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  ! Shared routines for extension and computation of gsmaps, avs, and ggrids
  private :: seq_mctext_gsmapIdentical
  private :: seq_mctext_gsmapExtend
  private :: seq_mctext_gsmapCreate
  private :: seq_mctext_avCreate

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer,public :: seq_mctext_decomp

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  character(*),parameter :: subName = '(seq_mctext_mct)'
  real(r8),parameter :: c1 = 1.0_r8

  !=======================================================================
contains
  !=======================================================================

  subroutine seq_map_init_exchange( comp, mapper, flow, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(component_type), intent(inout)      :: comp 
    type(seq_map)   , intent(inout), pointer :: mapper
    character(len=3), intent(in)             :: flow
    character(len=*), intent(in),optional    :: string
    !
    ! Local Variables
    !
    integer(IN)                :: ID_s
    integer(IN)                :: ID_d
    integer(IN)                :: ID_join
    integer(IN)                :: mapid, mapidmin, mapidmax
    integer(IN)                :: mpicom_s, mpicom_d, mpicom_join
    type(mct_gsmap) , pointer  :: gsmap_s
    type(mct_gsmap) , pointer  :: gsmap_d
    type(mct_gsmap)            :: gsmap_s_join
    type(mct_gsmap)            :: gsmap_d_join
    character(len=*),parameter :: subname = "(seq_map_init_rearrsplit) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    id_join = comp%cplcompid
    call seq_comm_getinfo(ID_join, mpicom=mpicom_join)

    if (flow == 'c2x') then
       gsmap_s => component_get_gsmap_cc(comp)
       gsmap_d => component_get_gsmap_cx(comp)
    end if
    if (flow == 'x2c') then
       gsmap_s => component_get_gsmap_cx(comp)
       gsmap_d => component_get_gsmap_cc(comp)
    end if

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then

       call seq_map_mapmatch(mapid, gsmap_s=gsmap_s, gsmap_d=gsmap_d, strategy="copy")

       if (mapid > 0) then
          call seq_map_mappoint(mapid, mapper)
       else
          call seq_map_mapinit(mapper, mpicom_join)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          if (flow == 'c2x') then
             mapper%gsmap_s => component_get_gsmap_cc(comp)
             mapper%gsmap_d => component_get_gsmap_cx(comp)
          end if
          if (flow == 'x2c') then
             mapper%gsmap_s => component_get_gsmap_cx(comp)
             mapper%gsmap_d => component_get_gsmap_cc(comp)
          end if
       endif

       if (seq_comm_iamroot(ID_join)) then
          write(logunit,'(2A,L2)') subname,' gsmaps ARE IDENTICAL, copyoption = ',mapper%copy_only
       endif

    else

       if (seq_comm_iamroot(ID_join)) write(logunit,'(2A)') subname,' gsmaps are not identical'

       if (flow == 'c2x') then
          id_s = comp%compid 
          id_d = cplid
       end if
       if (flow == 'x2c') then
          id_s = cplid
          id_d = comp%compid 
       end if
       call seq_comm_getinfo(ID_s   , mpicom=mpicom_s)
       call seq_comm_getinfo(ID_d   , mpicom=mpicom_d)
       call seq_comm_getinfo(ID_join, mpicom=mpicom_join)

       ! --- Extend gsmaps to join group of pes

       call seq_mctext_gsmapExtend(gsmap_s, mpicom_s, gsmap_s_join, mpicom_join, ID_join)
       call seq_mctext_gsmapExtend(gsmap_d, mpicom_d, gsmap_d_join, mpicom_join, ID_join)

       ! --- Initialize rearranger based on join gsmaps
       ! --- test for the gsmaps instead of the gsmap joins because the gsmap joins are temporary

       ! -------------------------------
       ! tcx  tcraig mapmatch is a problem here because we're comparing gsmaps that may not be defined
       !      on some pes.  first issue is whether gsmap_identical in underlying routine will abort.
       !      second issue is whether different pes return different values.  use mapidmin, mapidmax to
       !      confirm all mapids returned are the same.  if not, then just set mapid to -1 and compute
       !      a new rearranger.
       ! tcx  not clear this works all the time, so just do not do map matching here for time being
       !      Sept 2013.
       ! -------------------------------
       !       mapid = -1
       !       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")
       !       call shr_mpi_min(mapid,mapidmin,mpicom_join,subname//' min')
       !       call shr_mpi_max(mapid,mapidmax,mpicom_join,subname//' max')
       !       if (mapidmin /= mapidmax) mapid = -1
       ! -------------------------------

       ! --- Initialize rearranger
       ! --- the gsmap joins are temporary so store the regular gsmaps in the mapper
       call seq_map_mapinit(mapper, mpicom_join)
       mapper%rearrange_only = .true.
       mapper%strategy = "rearrange"
       if (flow == 'c2x') then
          mapper%gsmap_s => component_get_gsmap_cc(comp)
          mapper%gsmap_d => component_get_gsmap_cx(comp)
       end if
       if (flow == 'x2c') then
          mapper%gsmap_s => component_get_gsmap_cx(comp)
          mapper%gsmap_d => component_get_gsmap_cc(comp)
       end if
       call seq_map_gsmapcheck(gsmap_s_join, gsmap_d_join)
       call mct_rearr_init(gsmap_s_join, gsmap_d_join, mpicom_join, mapper%rearr)

       ! --- Clean up temporary gsmaps

       call mct_gsMap_clean(gsmap_s_join)
       call mct_gsMap_clean(gsmap_d_join)

    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
          mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
    endif

  end subroutine seq_map_init_exchange

 !===============================================================================

  subroutine seq_map_map_exchange( comp, flow, dom_flag, dom_tmp, string, msgtag )

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(component_type) , intent(inout)               :: comp 
    character(len=3)     , intent(in)                  :: flow
    logical              , intent(in),optional         :: dom_flag
    type(mct_gGrid)      , intent(in),optional, target :: dom_tmp
    character(len=*)     , intent(in),optional         :: string
    integer(IN)          , intent(in),optional         :: msgtag
    !
    ! Local Variables
    !
    type(seq_map)  , pointer :: mapper
    type(mct_aVect), pointer :: av_s
    type(mct_aVect), pointer :: av_d
    type(mct_gGrid), pointer :: dom_s
    type(mct_gGrid), pointer :: dom_d
    integer(IN),save         :: ltag    ! message tag for rearrange
    character(len=*),parameter :: subname = "(seq_map_map) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    if (flow == 'c2x') then
       if (present(dom_flag)) then
          dom_s   => component_get_dom_cc(comp)
          dom_d   => component_get_dom_cx(comp)
          ! Overwrite dom_d pointer if dom_tmp is present
          ! Needed for backwards compatibility with domain checker in component_init_cx
          if (present(dom_tmp)) then
             dom_d => dom_tmp
          end if
       else
          av_s   => component_get_c2x_cc(comp)
          av_d   => component_get_c2x_cx(comp)
       end if
       mapper => component_get_mapper_Cc2x(comp)
    end if
    if (flow == 'x2c') then
       if (present(dom_flag)) then
          dom_s  => component_get_dom_cx(comp)
          dom_d  => component_get_dom_cc(comp)
       else
          av_s   => component_get_x2c_cx(comp)
          av_d   => component_get_x2c_cc(comp)
       end if
       mapper => component_get_mapper_Cx2c(comp)
    end if

    if (present(msgtag)) then
       ltag = msgtag
    else
       ltag = 2000
    endif

    if (mapper%copy_only) then
       !-------------------------------------------
       ! COPY data
       !-------------------------------------------
       if (present(dom_flag)) then
          call mct_aVect_copy(aVin=dom_s%data, aVout=dom_d%data, vector=mct_usevector)
       else
          call mct_aVect_copy(aVin=av_s, aVout=av_d, vector=mct_usevector)
       end if

    else if (mapper%rearrange_only) then
       !-------------------------------------------
       ! REARRANGE data
       !-------------------------------------------
       if (present(dom_flag)) then
          call mct_rearr_rearrange(dom_s%data, dom_d%data, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall)
       else
          call mct_rearr_rearrange(av_s, av_d, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall)
       end if
    end if

  end subroutine seq_map_map_exchange

  !=======================================================================

  subroutine seq_mctext_gsmapInit(comp)

    ! This routine initializes a gsmap based on another gsmap potentially
    ! on other pes.  It addresses non-overlap of pes.

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(component_type), intent(inout) :: comp 
    !
    ! Local Variables
    !
    integer                  :: mpicom_cplid
    integer                  :: mpicom_old
    integer                  :: mpicom_new
    integer                  :: mpicom_join
    integer                  :: ID_old
    integer                  :: ID_new
    integer                  :: ID_join
    type(mct_gsMap), pointer :: gsmap_old
    type(mct_gsMap), pointer :: gsmap_new
    type(mct_gsMap)          :: gsmap_old_join   ! gsmap_old on joined id, temporary
    character(len=*),parameter :: subname = "(seq_mctext_gsmapInit) "
    !-----------------------------------------------------

    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)

    id_new  = cplid
    id_old  = comp%compid 
    id_join = comp%cplcompid

    mpicom_new  = mpicom_cplid
    mpicom_old  = comp%mpicom_compid
    mpicom_join = comp%mpicom_cplcompid

    gsmap_new => component_get_gsmap_cx(comp)
    gsmap_old => component_get_gsmap_cc(comp)

    call seq_comm_getinfo(ID_old ,mpicom=mpicom_old)
    call seq_comm_getinfo(ID_new ,mpicom=mpicom_new)
    call seq_comm_getinfo(ID_join,mpicom=mpicom_join)

    ! --- Set gsmaps
    ! ---   Extend the old one to now span all pes on ID_join
    ! ---   Create a new gsmap on pes associated with ID_new using info from the old one

    call seq_mctext_gsmapExtend(gsmap_old     , mpicom_old  , gsmap_old_join, mpicom_join, ID_join)
    call seq_mctext_gsmapCreate(gsmap_old_join, mpicom_join , gsmap_new     , mpicom_new , ID_new  )

    call mct_gsMap_clean(gsmap_old_join)

  end subroutine seq_mctext_gsmapInit

  !=======================================================================

  subroutine seq_mctext_avInit( comp, flow ) 

    !-----------------------------------------------------
    ! This routine initializes Avs that may need to be extended
    ! 
    ! Arguments
    !
    type(component_type), intent(inout) :: comp 
    character(len=3)    , intent(in)    :: flow 
    !
    ! Local Variables
    !
    integer                  :: lsize
    integer                  :: mpicom_cplid
    integer                  :: mpicom_new
    integer                  :: ID_old
    integer                  :: ID_new
    integer                  :: ID_join
    type(mct_aVect), pointer :: AV1_old
    type(mct_aVect), pointer :: AV1_new
    type(mct_gsmap), pointer :: gsmap_new
    character(len=*),parameter :: subname = "(seq_mctext_avInit) "
    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_comm_getinfo(CPLID ,mpicom=mpicom_CPLID)

    id_new  = cplid
    id_old  = comp%compid 
    id_join = comp%cplcompid

    mpicom_new  = mpicom_cplid

    gsmap_new => component_get_gsmap_cx(comp)

    if (flow == 'c2x') then
       av1_old => component_get_c2x_cc(comp)
       av1_new => component_get_c2x_cx(comp)
    end if
    if (flow == 'x2c') then
       av1_old => component_get_x2c_cc(comp)
       av1_new => component_get_x2c_cx(comp)
    end if

    ! --- Extend old avs and initialize new avs for use in the future

    lsize = 0
    if (seq_comm_iamin(ID_new)) then
       lsize = mct_gsMap_lsize(gsMap_new, mpicom_new)
    endif
    call seq_mctext_avExtend(AV1_old, ID_old, ID_join)
    call seq_mctext_avCreate(AV1_old, ID_old, AV1_new, ID_join, lsize)

  end subroutine seq_mctext_avInit

  !=======================================================================

  subroutine seq_mctext_gGridInit(comp, ggrid_new) 

    !-----------------------------------------------------
    ! This routine initializes gGrids that may need to be extended
    ! 
    ! Arguments
    !
    type(component_type), intent(inout) :: comp 
    type(mct_gGrid), optional, target, intent(inout) :: ggrid_new
    !
    ! Local Variables
    !
    integer                  :: mpicom_cplid
    integer                  :: lsize
    integer                  :: mpicom_new
    integer                  :: ID_old
    integer                  :: ID_new
    integer                  :: ID_join
    type(mct_gGrid), pointer :: GG1_old
    type(mct_gGrid), pointer :: GG1_new
    type(mct_gsmap), pointer :: gsmap_new
    character(len=*),parameter :: subname = "(seq_mctext_gGridInit) "
    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)

    id_new  = cplid
    id_old  = comp%compid 
    id_join = comp%cplcompid

    mpicom_new = mpicom_cplid

    gsmap_new => component_get_gsmap_cx(comp)

    gg1_old => component_get_dom_cc(comp)
    gg1_new => component_get_dom_cx(comp)

    ! --- Extend old ggrids and initialize new ggrids for use in the future

    lsize = 0
    if (seq_comm_iamin(ID_new)) then
       lsize = mct_gsMap_lsize(gsMap_new,mpicom_new)
    endif
    call seq_mctext_avExtend(GG1_old%data, ID_old, ID_join)
    
    if (present(ggrid_new)) then
       call mct_gGrid_init(GGrid=ggrid_new, CoordChars=seq_flds_dom_coord, OtherChars=seq_flds_dom_other, lsize=lsize )
       call mct_avect_zero(ggrid_new%data)
    else
       call mct_gGrid_init(GGrid=GG1_new, CoordChars=seq_flds_dom_coord, OtherChars=seq_flds_dom_other, lsize=lsize )
       call mct_avect_zero(GG1_new%data)
    end if

  end subroutine seq_mctext_gGridInit

  !=======================================================================

  subroutine seq_mctext_gsmapExtend(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !----------------------------------------------------------------
    ! Extend/Convert a gsmap from one mpicom to another mpicom that contains
    ! at least all the pes that gsmap uses, but with different ranks
    !----------------------------------------------------------------

    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_mctext_gsmapExtend) "
    integer :: n
    integer :: ngseg
    integer :: gsize
    integer :: msizei,msizeo
    integer :: mrank,mranko,mrankog   ! sets pe rank of root mpicomi pe in mpicomo
    integer :: mpigrpi,mpigrpo
    integer :: ierr
    integer, pointer :: pei(:),peo(:)
    integer, pointer :: start(:),length(:),peloc(:)

    mranko = -1

    ! --- create the new gsmap on the mpicomi root only

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_rank i')
       if (mrank == 0) then
          call mpi_comm_group(mpicomi,mpigrpi,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group i')
          call mpi_comm_group(mpicomo,mpigrpo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group o')
          call mpi_comm_size(mpicomi,msizei,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size i')
          call mpi_comm_size(mpicomo,msizeo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size o')

          ! --- setup the translation of pe numbers from the old gsmap(mpicom)
          ! --- to the new one, pei -> peo

          allocate(pei(0:msizei-1),peo(0:msizei-1))
          do n = 0,msizei-1
             pei(n) = n
          enddo

          peo = -1
          call mpi_group_translate_ranks(mpigrpi,msizei,pei,mpigrpo,peo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_group_translate_ranks')

          do n = 0,msizei-1
             if (peo(n) < 0 .or. peo(n) > msizeo-1) then
                write(logunit,*) subname,' peo out of bounds ',peo(n),msizeo
                call shr_sys_abort()
             endif
          enddo

          mranko = peo(0)

          ! --- compute the new gsmap which has the same start and length values
          ! --- but peloc is now the mapping of pei to peo

          ngseg = gsmapi%ngseg
          gsize = gsmapi%gsize
          allocate(start(ngseg),length(ngseg),peloc(ngseg))
          do n = 1,ngseg
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
             peloc(n)  = peo(gsmapi%pe_loc(n))
          enddo

          ! --- initialize the gsmap on the root pe

          call mct_gsmap_init(gsmapo,compido,ngseg,gsize,start,length,peloc)

          deallocate(pei,peo,start,length,peloc)
       endif
    endif

    ! --- broadcast via allreduce the mpicomi root pe in mpicomo space
    ! --- mranko is -1 except on the root pe where is it peo of that pe

    call mpi_allreduce(mranko,mrankog,1,MPI_INTEGER,MPI_MAX,mpicomo,ierr)
    call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_allreduce max')

    ! --- broadcast the gsmap to all pes in mpicomo from mrankog

    call mct_gsmap_bcast(gsmapo, mrankog, mpicomo)

    ! tcx summarize decomp info


  end subroutine seq_mctext_gsmapExtend

  !=======================================================================

  subroutine seq_mctext_gsmapCreate(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !---------------------------------------------------------------------
    ! creates a new gsmap on a subset of pes, requires setting a new decomp
    !---------------------------------------------------------------------

    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_mctext_gsmapCreate) "
    integer :: n,m,k
    integer :: ktot            ! number of active cells in gsmap
    integer :: apesi, apeso    ! number of active pes in gsmap
    integer ::        lsizeo   ! local size for lindex
    integer :: ngsegi,ngsego   ! ngseg of mpicomi, mpicomo
    integer :: gsizei,gsizeo   ! gsize of mpicomi, mpicomo
    integer :: msizei,msizeo   ! size of mpicomi, mpicomo
    integer :: mranki,mranko   ! rank in mpicomi, mpicomo
    integer :: ierr
    integer :: decomp_type
    integer, pointer :: start(:),length(:),peloc(:),perm(:),gindex(:),lindex(:)
    real(r8):: rpeloc
    logical :: gsmap_bfbflag = .false. ! normally this should be set to false

    ! --- create a new gsmap on new pes based on the old gsmap
    ! --- gsmapi must be known on all mpicomo pes, compute the same 
    ! --- thing on all pes in parallel

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mranki,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank i')
       call mpi_comm_size(mpicomi,msizei,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size i')
       call mpi_comm_rank(mpicomo,mranko,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank o')
       call mpi_comm_size(mpicomo,msizeo,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size o')

       ngsegi = gsmapi%ngseg
       gsizei = gsmapi%gsize
       gsizeo = gsizei
       call mct_gsMap_activepes(gsmapi,apesi)

       decomp_type = 0

       if (seq_mctext_decomp == 0) then
          if (msizeo == apesi) then      ! preserve segments and decomp
             ! For testing - set decomp_type to 1 - to have gsmapi and gsmapo identical
             if (gsmap_bfbflag) then
                decomp_type = 1     ! better in cpl to have all decomps "same-ish"
             else
                decomp_type = 2
             end if
          elseif (ngsegi >= msizeo) then ! preserve segments, new decomp
             decomp_type = 2
          else                           ! new segments
             decomp_type = 3
          endif
       else
          decomp_type = seq_mctext_decomp
       endif

       !tcx       decomp_type = 3 ! over ride setting above for testing
       !       if (mranko == 0) write(logunit,'(2A,4I)') trim(subname),' decomp_type =',decomp_type,ngsegi,msizeo,apesi

       select case (decomp_type)

       case(1)   ! --- preserve segments and decomp ---------------------

          ! -- copy the gsmap and translate the pes
          call mct_gsMap_copy(gsmapi,gsmapo)
          ngsego = ngsegi
          do n = 1,ngsego
             gsmapo%pe_loc(n) = mod(gsmapo%pe_loc(n),msizeo)    ! translate pes 1:1 from old to new
          enddo

       case(2)   ! --- preserve segments, new decomp --------------------

          ! --- preserve segments, sort the start and length, assign a new pe list
          ngsego = ngsegi
          allocate(start(ngsego),length(ngsego),peloc(ngsego),perm(ngsego))
          do n = 1,ngsego
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
          enddo
          ! --- sort gsmap to minimize permute cost in mct
          call mct_indexset(perm)
          call mct_indexsort(ngsego,perm,start)
          call mct_permute(start,perm,ngsego)
          call mct_permute(length,perm,ngsego)
          ! --- give each pe "equal" number of segments, use reals to avoid integer overflow
          do n = 1,ngsego
             rpeloc = (((msizeo*c1)*((n-1)*c1))/(ngsego*c1))      ! give each pe "equal" number of segments, use reals to avoid integer overflow
             peloc(n) = int(rpeloc)
          enddo
          call mct_gsmap_init(gsmapo,ngsego,start,length,peloc,0,mpicomo,compido,gsizeo)
          deallocate(start,length,peloc,perm)

       case(3)   ! --- new segments, new decomp -------------------------

          ! --- new segments, compute gindex, then parse the gridcells out evenly

          k = 0
          do n = 1,ngsegi
             do m = 1,gsmapi%length(n)
                k = k + 1
                if (k > gsizei) then
                   write(logunit,*) trim(subname),' ERROR in gindex ',k,gsizei
                   call shr_sys_abort()
                endif
             enddo
          enddo
          ktot = k

          allocate(gindex(ktot),perm(ktot))  

          k = 0
          do n = 1,ngsegi
             do m = 1,gsmapi%length(n)
                k = k + 1
                gindex(k) = gsmapi%start(n) + m - 1
             enddo
          enddo
          call mct_indexset(perm)
          call mct_indexsort(ktot,perm,gindex)
          call mct_permute(gindex,perm,ktot)

          k = 0
          do m = 0,msizeo-1
             lsizeo = ktot/msizeo
             if (m < (ktot - lsizeo*msizeo)) lsizeo = lsizeo + 1
             if (mranko == m) then
                allocate(lindex(lsizeo))
                if (k+lsizeo > ktot) then
                   write(logunit,*) trim(subname),' ERROR: decomp out of bounds ',mranko,k,lsizeo,ktot
                   call shr_sys_abort()
                endif
                lindex(1:lsizeo) = gindex(k+1:k+lsizeo)
                !                write(logunit,*) trim(subname),' decomp is ',mranko,lsizeo,k+1,k+lsizeo
             endif
             k = k + lsizeo
          enddo
          if (k /= ktot) then
             write(logunit,*) trim(subname),' ERROR: decomp incomplete ',k,ktot
             call shr_sys_abort()
          endif

          call mct_gsmap_init(gsmapo,lindex,mpicomo,compido,size(lindex),gsizeo)
          deallocate(gindex,perm,lindex)

       case default   ! --- unknown ---
          write(logunit,*) trim(subname),' ERROR decomp_type unknown ',decomp_type
          call shr_sys_abort(trim(subname)//' ERROR decomp_type unknown')

       end select

       if (mranko == 0) then
          write(logunit,102) trim(subname),' created new gsmap decomp_type =',decomp_type
          write(logunit,102) trim(subname),'   ngseg/gsize        = ', &
               mct_gsmap_ngseg(gsmapo),mct_gsmap_gsize(gsmapo)
          call mct_gsmap_activepes(gsmapo,apeso)
          write(logunit,102) trim(subname),'   mpisize/active_pes = ', &
               msizeo,apeso
          write(logunit,102) trim(subname),'   avg seg per pe/ape = ', &
               mct_gsmap_ngseg(gsmapo)/msizeo,mct_gsmap_ngseg(gsmapo)/apeso
          write(logunit,102) trim(subname),'   nlseg/maxnlsegs    = ', &
               mct_gsmap_nlseg(gsmapo,0),mct_gsmap_maxnlseg(gsmapo)
102       format(2A,2I8)
       endif

       !       if (.not. mct_gsmap_increasing(gsmapo) ) then
       !          write(logunit,*) trim(subname),' ERROR: gsmapo not increasing'
       !          call shr_sys_abort()
       !       endif

    endif

  end subroutine seq_mctext_gsmapCreate

  !=======================================================================

  subroutine seq_mctext_avExtend(AVin,IDin,ID)

    !-----------------------------------------------------------------------
    ! Extend an AV to a larger set of pes or
    ! Initialize an AV on another set of pes
    !
    ! Arguments
    !
    type(mct_aVect), intent(INOUT):: AVin
    integer         ,intent(IN)   :: IDin ! ID associated with AVin
    integer        , intent(IN)   :: ID   ! ID to initialize over
    !
    ! Local variables
    !
    character(len=*),parameter :: subname = "(seq_mctext_avExtend) "
    integer :: mpicom
    integer :: rank,rank2
    integer :: lsizei, lsizen
    integer :: srank,srankg
    integer :: ierr
    integer :: nints
    character(len=CXX) :: iList,rList
    !-----------------------------------------------------------------------

    call seq_comm_getinfo(ID,mpicom=mpicom,iam=rank)

    ! --- lsizen is the size of the newly initialized AV, zero is valid
    ! --- lsizei is -1 on any peszero on any pes where AV is not yet initialized

    lsizei = -1  
    if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
    lsizen = 0

    ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
    ! --- set the pe and broadcast it to all other pes using mpi_allreduce

    srank = -1
    srankg = -1
    if (lsizei > 0) srank = rank

    call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

    if (srankg < 0) then
       write(logunit,*) subname,' WARNING AVin empty '
       return
    endif

    ! --- set the iList and rList from the broadcast pe (srankg) and 
    ! --- broadcast the lists

    iList = " "
    rList = " "
    if (rank == srankg) then
       if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
       if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
    endif

    call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
    call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

    ! --- now allocate the AV on any pes where the orig size is zero.  those
    ! --- should be pes that either have no data and may have been allocated
    ! --- before (no harm in doing it again) or have never been allocated

    if (lsizei <= 0) then
       if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
          call mct_aVect_init(AVin,iList=iList,rList=rList,lsize=lsizen)
       elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
          call mct_aVect_init(AVin,iList=iList,lsize=lsizen)
       elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
          call mct_aVect_init(AVin,rList=rList,lsize=lsizen)
       endif
    endif

  end subroutine seq_mctext_avExtend

  !=======================================================================

  subroutine seq_mctext_avCreate(AVin,IDin,AVout,ID,lsize)

    !-----------------------------------------------------------------------
    ! Extend an AV to a larger set of pes or
    ! Initialize an AV on another set of pes
    !-----------------------------------------------------------------------

    implicit none
    type(mct_aVect), intent(INOUT):: AVin
    integer         ,intent(IN)   :: IDin ! ID associated with AVin
    type(mct_aVect), intent(INOUT):: AVout
    integer        , intent(IN)   :: ID   ! ID to initialize over
    integer        , intent(IN)   :: lsize

    ! Local variables

    character(len=*),parameter :: subname = "(seq_mctext_avCreate) "
    integer :: mpicom
    integer :: rank,rank2
    integer :: lsizei, lsizen
    integer :: srank,srankg
    integer :: ierr
    integer :: nints
    character(len=CXX) :: iList,rList

    call seq_comm_getinfo(ID,mpicom=mpicom,iam=rank)

    ! --- lsizen is the size of the newly initialized AV, zero is valid

    lsizei = -1  
    if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
    lsizen = lsize

    ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
    ! --- set the pe and broadcast it to all other pes

    srank = -1
    srankg = -1
    if (lsizei > 0) srank = rank

    call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

    if (srankg < 0) then
       write(logunit,*) subname,' ERROR AVin not initialized '
       call shr_sys_abort()
    endif

    ! --- set the iList and rList from the broadcast pe (srankg) and 
    ! --- broadcast the lists

    iList = " "
    rList = " "
    if (rank == srankg) then
       if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
       if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
    endif

    call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
    call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

    ! --- now allocate the AV on all pes.  the AV should not exist before.
    ! --- If it does, mct should die.

    if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
       call mct_aVect_init(AVout,iList=iList,rList=rList,lsize=lsizen)
    elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
       call mct_aVect_init(AVout,iList=iList,lsize=lsizen)
    elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
       call mct_aVect_init(AVout,rList=rList,lsize=lsizen)
    endif

  end subroutine seq_mctext_avCreate

  !=======================================================================

  logical function seq_mctext_gsmapIdentical(gsmap1,gsmap2)

    implicit none
    type(mct_gsMap), intent(IN):: gsmap1
    type(mct_gsMap), intent(IN):: gsmap2

    ! Local variables

    character(len=*),parameter :: subname = "(seq_mctext_gsmapIdentical) "
    integer :: n
    logical :: identical

    !-----------------------

    identical = .true.

    ! --- continue compare ---
    if (identical) then
       if (mct_gsMap_gsize(gsmap1) /= mct_gsMap_gsize(gsmap2)) identical = .false.
       if (mct_gsMap_ngseg(gsmap1) /= mct_gsMap_ngseg(gsmap2)) identical = .false.
    endif

    ! --- continue compare ---
    if (identical) then
       do n = 1,mct_gsMap_ngseg(gsmap1)
          if (gsmap1%start(n)  /= gsmap2%start(n) ) identical = .false.
          if (gsmap1%length(n) /= gsmap2%length(n)) identical = .false.
          if (gsmap1%pe_loc(n) /= gsmap2%pe_loc(n)) identical = .false.
       enddo
    endif

    seq_mctext_gsmapIdentical = identical

  end function seq_mctext_gsmapIdentical

end module cplcomp_exchange_mod
