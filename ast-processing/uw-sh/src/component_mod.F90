module component_mod

  !----------------------------------------------------------------------------
  ! share code & libs
  !----------------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => SHR_KIND_R8 
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use shr_const_mod,    only: shr_const_cday
  use shr_file_mod,     only: shr_file_setLogLevel, shr_file_setLogUnit
  use shr_file_mod,     only: shr_file_setIO, shr_file_getUnit
  use shr_scam_mod,     only: shr_scam_checkSurface
  use shr_mpi_mod,      only: shr_mpi_min, shr_mpi_max
  use shr_mem_mod,      only: shr_mem_init, shr_mem_getusage
  use shr_cal_mod,      only: shr_cal_date2ymd
  use shr_orb_mod,      only: shr_orb_params
  use shr_reprosum_mod, only: shr_reprosum_setopts
  use seq_comm_mct,     only: GLOID, CPLID, logunit
  use seq_comm_mct,     only: seq_comm_iamin, seq_comm_namelen, num_inst_frc
  use seq_comm_mct,     only: seq_comm_suffix, seq_comm_name, seq_comm_setnthreads
  use seq_comm_mct,     only: seq_comm_getinfo => seq_comm_setptrs
  use seq_comm_mct,     only: seq_comm_petlist 
  use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
  use seq_infodata_mod, only: seq_infodata_exchange, seq_infodata_type
  use seq_diag_mct,     only: seq_diag_avect_mct 
  use seq_map_type_mod  
  use seq_map_mod
  use t_drv_timers_mod
  use component_type_mod
  use seq_cdata_mod,    only : seq_cdata
  use mct_mod   ! mct_ wrappers for mct lib
  use perf_mod
  use ESMF




  implicit none

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
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: component_init_pre
  public :: component_init_cc            ! mct and esmf versions
  public :: component_init_cx
  public :: component_init_aream
  public :: component_init_areacor
  public :: component_run                 ! mct and esmf versions
  public :: component_final               ! mct and esmf versions 
  public :: component_exch
  public :: component_diag


  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

   logical  :: iamroot_GLOID, iamroot_CPLID         ! GLOID, CPLID masterproc
   logical  :: iamin_CPLID                          ! true => pe associated with CPLID
   integer  :: mpicom_GLOID, mpicom_CPLID           ! GLOID, CPLID mpi communicator 
   integer  :: nthreads_GLOID, nthreads_CPLID
   logical  :: drv_threading

   character(*), parameter :: subname = '(component_mod)'
   character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
   character(*), parameter :: F0L = "('"//subname//" : ', A, L6 )"
   character(*), parameter :: F0I = "('"//subname//" : ', A, 2i8 )"
   character(*), parameter :: F0R = "('"//subname//" : ', A, 2g23.15 )"

  !===============================================================================

contains

  !===============================================================================

  subroutine component_init_pre(comp, compid, cplcompid, cplallcompid, &
       infodata, ntype)

    !---------------------------------------------------------------
    ! Initialize driver rearrangers and AVs on driver
    ! Initialize cdata_*x data
    ! Zero out x2*_** in case it never gets used then it'll produce zeros in diags
    ! For ensembles, create only a single dom_*x for the coupler based on the
    !   first ensemble member.  otherwise, just extend the dom_** and dom_*x to
    !   other ensemble members.
    !
    ! Arguments
    type(component_type)     , intent(inout)         :: comp(:)
    integer                  , intent(in)            :: compid(:)
    integer                  , intent(in)            :: cplcompid(:)
    integer                  , intent(in)            :: cplallcompid
    type (seq_infodata_type) , intent(inout), target :: infodata
    character(len=3)         , intent(in)            :: ntype
    !
    ! Local Variables
    integer  :: eci       ! index
    !---------------------------------------------------------------

    ! initialize module variables (this is repetitive here- but does not require a different routine)

    call seq_infodata_getdata(infodata, drv_threading=drv_threading)
    call seq_comm_getinfo(GLOID, mpicom=mpicom_GLOID, iamroot=iamroot_GLOID, nthreads=nthreads_GLOID)
    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID, iamroot=iamroot_CPLID, nthreads=nthreads_CPLID)
    iamin_CPLID = seq_comm_iamin(CPLID)
   
    ! Initialize component type variables
    do eci = 1,size(comp)

       comp(eci)%compid       = compid(eci)
       comp(eci)%cplcompid    = cplcompid(eci)
       comp(eci)%cplallcompid = cplallcompid

       call seq_comm_getinfo(comp(eci)%cplallcompid, mpicom=comp(eci)%mpicom_cplallcompid)
       call seq_comm_getinfo(comp(eci)%cplcompid   , mpicom=comp(eci)%mpicom_cplcompid)
       call seq_comm_getinfo(comp(eci)%compid      , mpicom=comp(eci)%mpicom_compid)
       call seq_comm_getinfo(comp(eci)%compid      , iamroot=comp(eci)%iamroot_compid)
       call seq_comm_getinfo(comp(eci)%compid      , nthreads=comp(eci)%nthreads_compid)

       comp(eci)%iamin_compid       =  seq_comm_iamin (comp(eci)%compid)
       comp(eci)%iamin_cplcompid    =  seq_comm_iamin (comp(eci)%cplcompid)
       comp(eci)%iamin_cplallcompid =  seq_comm_iamin (comp(eci)%cplallcompid)
       comp(eci)%suffix             =  seq_comm_suffix(comp(eci)%compid)
       comp(eci)%name               =  seq_comm_name  (comp(eci)%compid)
       comp(eci)%ntype              =  ntype(1:3)
       comp(eci)%oneletterid        =  ntype(1:1)

       if (eci == 1) then
          allocate(comp(1)%dom_cx)
          allocate(comp(1)%gsmap_cx)
       else
          comp(eci)%dom_cx   => comp(1)%dom_cx 
          comp(eci)%gsmap_cx => comp(1)%gsmap_cx 
       end if

       ! Set cdata_cc - unique for each instance
       allocate(comp(eci)%dom_cc)
       allocate(comp(eci)%gsmap_cc)
       allocate(comp(eci)%cdata_cc)
       comp(eci)%cdata_cc%name     = 'cdata_'//ntype(1:1)//ntype(1:1)
       comp(eci)%cdata_cc%ID       =  comp(eci)%compid
       comp(eci)%cdata_cc%mpicom   =  comp(eci)%mpicom_compid
       comp(eci)%cdata_cc%dom      => comp(eci)%dom_cc 
       comp(eci)%cdata_cc%gsmap    => comp(eci)%gsmap_cc
       comp(eci)%cdata_cc%infodata => infodata

       ! Determine initial value of comp_present in infodata - to do - add this to component 

       if (comp(1)%oneletterid == 'a') call seq_infodata_getData(infodata, atm_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'l') call seq_infodata_getData(infodata, lnd_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'i') call seq_infodata_getData(infodata, ice_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'o') call seq_infodata_getData(infodata, ocn_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'r') call seq_infodata_getData(infodata, rof_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'g') call seq_infodata_getData(infodata, glc_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'w') call seq_infodata_getData(infodata, wav_present=comp(eci)%present)

    end do

  end subroutine component_init_pre

  !===============================================================================

  subroutine component_init_cc(Eclock, comp, comp_init, infodata, NLFilename, &
       seq_flds_x2c_fluxes, seq_flds_c2x_fluxes)

    !---------------------------------------------------------------
    !
    ! Arguments
    type(ESMF_Clock)         , intent(inout) :: EClock
    type(component_type)     , intent(inout) :: comp(:)
    interface 
       subroutine comp_init( Eclock, cdata, x2c, c2x, nlfilename)
         use ESMF         , only: ESMF_Clock
         use seq_cdata_mod, only: seq_cdata
         use mct_mod      , only: mct_avect
         implicit none
         type(ESMF_Clock), intent(inout) :: EClock
         type(seq_cdata) , intent(inout) :: cdata
         type(mct_aVect) , intent(inout) :: x2c
         type(mct_aVect) , intent(inout) :: c2x   
         character(len=*), optional, intent(IN) :: NLFilename ! Namelist filename
       end subroutine comp_init
    end interface 
    type (seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)         , intent(in)           :: NLFilename 
    character(len=*)         , intent(in), optional :: seq_flds_x2c_fluxes
    character(len=*)         , intent(in), optional :: seq_flds_c2x_fluxes
    !
    ! Local Variables
    integer :: k1, k2
    integer :: eci
    !---------------------------------------------------------------

    ! **** Initialize component - this initializes  x2c_cc and c2x_cc ***
    ! the following will call the appropriate comp_init_mct routine

    if (comp(1)%iamin_cplallcompid) then
       call seq_infodata_exchange(infodata, comp(1)%cplallcompid, &
            'cpl2'//comp(1)%ntype(1:3)//'_init')
    end if

    ! The following initializes the component instance cdata_cc (gsmap and dom), 
    ! x2c_cc and c2x_cc

    do eci = 1,size(comp)
       if (iamroot_CPLID .and. comp(eci)%present) then
          write(logunit,F00) 'Initialize component '//trim(comp(eci)%ntype)
          call shr_sys_flush(logunit)
       endif

       if (.not. associated(comp(eci)%x2c_cc)) allocate(comp(eci)%x2c_cc)
       if (.not. associated(comp(eci)%c2x_cc)) allocate(comp(eci)%c2x_cc)

       if (comp(eci)%iamin_compid .and. comp(eci)%present) then
          if (drv_threading) call seq_comm_setnthreads(comp(eci)%nthreads_compid)
          call shr_sys_flush(logunit)
          
          if (present(seq_flds_x2c_fluxes)) then
             call mct_avect_vecmult(comp(eci)%x2c_cc, comp(eci)%drv2mdl, seq_flds_x2c_fluxes, mask_spval=.true.)
          end if

          call comp_init( EClock, comp(eci)%cdata_cc, comp(eci)%x2c_cc, comp(eci)%c2x_cc, &
               NLFilename=NLFilename )
          
          if (present(seq_flds_c2x_fluxes)) then
             call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)
          end if
          
          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       end if
    end do

    if (comp(1)%iamin_cplcompid) then
       call seq_infodata_exchange(infodata, comp(1)%cplcompid, &
            comp(1)%ntype(1:3)//'2cpl_init')
    endif

    ! Determine final value of comp_present in infodata (after component initialization)

    do eci = 1,size(comp) 
       if (comp(1)%oneletterid == 'a') call seq_infodata_getData(infodata, atm_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'l') call seq_infodata_getData(infodata, lnd_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'i') call seq_infodata_getData(infodata, ice_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'o') call seq_infodata_getData(infodata, ocn_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'r') call seq_infodata_getData(infodata, rof_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'g') call seq_infodata_getData(infodata, glc_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'w') call seq_infodata_getData(infodata, wav_present=comp(eci)%present)
    end do


    ! Initialize aream, set it to area for now until maps are read
    !   in some cases, maps are not read at all !!
    ! Entire domain must have reasonable values before calling xxx2xxx init

    do eci = 1,size(comp)
       if (comp(eci)%iamin_compid .and. comp(eci)%present) then
          if (drv_threading) call seq_comm_setnthreads(comp(eci)%nthreads_compid)
          k1 = mct_aVect_indexRa(comp(eci)%cdata_cc%dom%data, "area"  ,perrWith='aa area ')
          k2 = mct_aVect_indexRa(comp(eci)%cdata_cc%dom%data, "aream" ,perrWith='aa aream')

          comp(eci)%cdata_cc%dom%data%rAttr(k2,:) = comp(eci)%cdata_cc%dom%data%rAttr(k1,:)

          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       endif
    end do

  end subroutine component_init_cc

  !===============================================================================


  !===============================================================================


  !===============================================================================
    
  subroutine component_init_cx(comp, infodata) 

    !---------------------------------------------------------------
    ! Uses
    use cplcomp_exchange_mod, only: seq_mctext_gsmapinit, seq_mctext_avInit
    use cplcomp_exchange_mod, only: seq_mctext_avExtend, seq_mctext_gGridInit
    use cplcomp_exchange_mod, only: seq_map_init_exchange, seq_map_map_exchange
    use seq_domain_mct,       only: seq_domain_compare
    use mct_mod,              only: mct_ggrid_clean
    !
    ! Arguments
    type(component_type)     , intent(inout) :: comp(:)
    type (seq_infodata_type) , intent(inout) :: infodata
    !
    ! Local Variables
    integer         :: eci
    integer         :: rc        ! return code
    type(mct_gGrid) :: dom_tmp   ! temporary
    !---------------------------------------------------------------

    ! Initialize driver rearrangers and AVs on driver
    ! Initialize cdata_*x data
    ! Zero out x2*_** in case it never gets used then it'll produce zeros in diags
    ! For ensembles, create only a single dom_*x for the coupler based on the
    !   first ensemble member.  otherwise, just extend the dom_** and dom_*x to
    !   other ensemble members.

    do eci = 1,size(comp)
       if (comp(eci)%present) then

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             call shr_sys_flush(logunit)
          end if

          if (comp(eci)%iamin_cplcompid) then

             ! Create gsmap_cx (note that comp(eci)%gsmap_cx all point to comp(1)%gsmap_cx
             ! This will only be valid on the coupler pes
             if (eci == 1) then
                if (iamroot_CPLID) then
                   write(logunit,F0I) 'creating gsmap_cx for '//comp(eci)%ntype(1:3)
                   call shr_sys_flush(logunit)
                end if
                call seq_mctext_gsmapInit(comp(1))
             endif

             ! Create mapper_Cc2x and mapper_Cx2c
             allocate(comp(eci)%mapper_Cc2x, comp(eci)%mapper_Cx2c)  
             if (iamroot_CPLID) then
                write(logunit,F0I) 'Initializing mapper_C'//comp(eci)%ntype(1:1)//'2x',eci
                call shr_sys_flush(logunit)
             end if
             call seq_map_init_exchange(comp(eci), flow='c2x', mapper=comp(eci)%mapper_Cc2x)    
             if (iamroot_CPLID) then
                write(logunit,F0I) 'Initializing mapper_Cx2'//comp(eci)%ntype(1:1),eci
                call shr_sys_flush(logunit)
             end if
             call seq_map_init_exchange(comp(eci), flow='x2c', mapper=comp(eci)%mapper_Cx2c)  

             ! Create x2c_cx and c2x_cx
             allocate(comp(eci)%x2c_cx, comp(eci)%c2x_cx)
             call seq_mctext_avinit(comp(eci), flow='x2c')
             call seq_mctext_avinit(comp(eci), flow='c2x')

             ! Create dom_cx (note that  comp(eci)%dom_cx all point to  comp(1)%dom_cx
             ! Then verify other ensembles have same domain by comparing to dom_cx
             if (eci == 1) then  ! create dom_cx
                if (iamroot_CPLID) then
                   write(logunit,F0I) 'creating dom_cx'
                   call shr_sys_flush(logunit)
                end if
                call seq_mctext_gGridInit(comp(1))
                call seq_map_map_exchange(comp(1), flow='c2x', dom_flag=.true., msgtag=comp(1)%cplcompid*100+1*10+1)
             else if (eci > 1) then  
                if (iamroot_CPLID) then
                   write(logunit,F0I) 'comparing comp domain ensemble number ',eci
                   call shr_sys_flush(logunit)
                end if
                call seq_mctext_avExtend(comp(eci)%dom_cx%data, cplid, comp(eci)%cplcompid)
                call seq_mctext_gGridInit(comp(eci), dom_tmp)
                call seq_map_map_exchange(comp(eci), flow='c2x', dom_flag=.true., dom_tmp=dom_tmp)
                if (iamin_CPLID) then
                   call seq_domain_compare(comp(eci)%dom_cx, dom_tmp, mpicom_CPLID)
                end if
                call mct_ggrid_clean(dom_tmp,rc)
             endif

             call mct_avect_zero(comp(eci)%x2c_cc)
             call mct_avect_zero(comp(eci)%x2c_cx)

          end if ! if comp(eci)%iamin_cplcompid
       end if  ! if comp(eci)%present
    end do  ! end of eci loop

  end subroutine component_init_cx

  !===============================================================================

  subroutine component_init_aream(infodata, rof_c2_ocn, samegrid_ao, samegrid_al, samegrid_ro)

    !---------------------------------------------------------------
    ! Description
    ! Update (read) aream in domains where appropriate - ON cpl pes
    !
    ! Uses
    use prep_ocn_mod,       only : prep_ocn_get_mapper_Fa2o
    use prep_lnd_mod,       only : prep_lnd_get_mapper_Sa2l
    use prep_ice_mod,       only : prep_ice_get_mapper_SFo2i
    use prep_glc_mod,       only : prep_glc_get_mapper_SFl2g
    use component_type_mod, only : atm, lnd, ice, ocn, rof, glc
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: rof_c2_ocn
    logical                  , intent(in)    :: samegrid_ao
    logical                  , intent(in)    :: samegrid_al
    logical                  , intent(in)    :: samegrid_ro
    !
    ! Local variables
    type(mct_gsmap), pointer :: gsmap_s, gsmap_d
    type(mct_ggrid), pointer :: dom_s, dom_d
    type(seq_map)  , pointer :: mapper_Fa2o
    type(seq_map)  , pointer :: mapper_Sa2l
    type(seq_map)  , pointer :: mapper_SFo2i
    type(seq_map)  , pointer :: mapper_SFl2g
    logical                  :: atm_present ! atm present flag
    logical                  :: lnd_present ! lnd present flag
    logical                  :: ocn_present ! ocn present flag
    logical                  :: ice_present ! ice present flag
    logical                  :: glc_present ! glc present flag
    integer                  :: ka,km
    !---------------------------------------------------------------

    ! Note that the following is assumed to hold - all gsmaps_cx for a given 
    ! instance of a component (e.g. atm(i)) are identical on the coupler processes

    mapper_Fa2o  => prep_ocn_get_mapper_Fa2o()
    mapper_Sa2l  => prep_lnd_get_mapper_Sa2l()
    mapper_SFo2i => prep_ice_get_mapper_SFo2i()
    mapper_SFl2g => prep_glc_get_mapper_SFl2g()

    call seq_infodata_GetData( infodata, &
         atm_present=atm_present,        &
         ocn_present=ocn_present,        &
         ice_present=ice_present,        &
         lnd_present=lnd_present,        &
         glc_present=glc_present)

    if (atm_present .and. ocn_present) then
       if (samegrid_ao) then
          dom_s  => component_get_dom_cx(atm(1))   !dom_ax
          dom_d  => component_get_dom_cx(ocn(1))   !dom_ox
          ka = mct_aVect_indexRa(dom_s%data, "area" )
          km = mct_aVect_indexRa(dom_s%data, "aream" )
          dom_s%data%rAttr(km,:) = dom_s%data%rAttr(ka,:)

          call seq_map_map(mapper_Fa2o, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream')
       else
          gsmap_s => component_get_gsmap_cx(ocn(1)) ! gsmap_ox
          gsmap_d => component_get_gsmap_cx(atm(1)) ! gsmap_ax
          dom_s   => component_get_dom_cx(ocn(1))   ! dom_ox
          dom_d   => component_get_dom_cx(atm(1))   ! dom_ax

          call seq_map_readdata('seq_maps.rc','ocn2atm_fmapname:', mpicom_CPLID, CPLID, &
               gsmap_s=gsmap_s, av_s=dom_s%data, avfld_s='aream', filefld_s='area_a', &
               gsmap_d=gsmap_d, av_d=dom_d%data, avfld_d='aream', filefld_d='area_b', &
               string='ocn2atm aream initialization')
       endif
    end if

    if (ice_present .and. ocn_present) then
       dom_s  => component_get_dom_cx(ocn(1))   !dom_ox
       dom_d  => component_get_dom_cx(ice(1))   !dom_ix

       call seq_map_map(mapper_SFo2i, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream') 
    endif

    if (rof_c2_ocn) then
       if (.not.samegrid_ro) then
          gsmap_s => component_get_gsmap_cx(rof(1)) ! gsmap_rx
          dom_s   => component_get_dom_cx(rof(1))   ! dom_rx

          call seq_map_readdata('seq_maps.rc', 'rof2ocn_rmapname:',mpicom_CPLID, CPLID, &
               gsmap_s=gsmap_s, av_s=dom_s%data, avfld_s='aream', filefld_s='area_a', &
               string='rof2ocn aream initialization')
       endif
    end if

    if (lnd_present .and. atm_present) then
       if (samegrid_al) then
          dom_s  => component_get_dom_cx(atm(1))   !dom_ax
          dom_d  => component_get_dom_cx(lnd(1))   !dom_lx

          call seq_map_map(mapper_Sa2l, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream')
       else
          gsmap_d => component_get_gsmap_cx(lnd(1)) ! gsmap_lx
          dom_d   => component_get_dom_cx(lnd(1))   ! dom_lx

          call seq_map_readdata('seq_maps.rc','atm2lnd_fmapname:',mpicom_CPLID, CPLID, &
               gsmap_d=gsmap_d, av_d=dom_d%data, avfld_d='aream', filefld_d='area_b', &
               string='atm2lnd aream initialization')
       endif
    end if

    if (lnd_present .and. glc_present) then
       dom_s  => component_get_dom_cx(lnd(1))   !dom_lx
       dom_d  => component_get_dom_cx(glc(1))   !dom_gx

       call seq_map_map(mapper_SFl2g, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream')
    endif

  end subroutine component_init_aream

  !===============================================================================

  subroutine component_init_areacor(comp, samegrid, seq_flds_c2x_fluxes)
    !---------------------------------------------------------------
    ! COMPONENT PES and CPL/COMPONENT (for exchange only)
    !
    ! Uses
    use seq_domain_mct, only : seq_domain_areafactinit
    !
    ! Arguments
    type(component_type) , intent(inout) :: comp(:)
    logical              , intent(in)    :: samegrid
    character(len=*)     , intent(in)    :: seq_flds_c2x_fluxes
    !
    ! Local Variables
    integer :: eci, num_inst
    !---------------------------------------------------------------

    num_inst = size(comp) 
    do eci = 1,num_inst

       ! For joint cpl-component pes
       if (comp(eci)%iamin_cplcompid) then

          ! Map component domain from coupler to component processes
          call seq_map_map(comp(eci)%mapper_Cx2c, comp(eci)%dom_cx%data, &
               comp(eci)%dom_cc%data, msgtag=comp(eci)%cplcompid*100+eci*10+5)

          ! For only component pes
          if (comp(eci)%iamin_compid) then

             ! Allocate and initialize area correction factors on component processes
             ! Note that the following call allocates comp(eci)%mld2drv(:) and comp(eci)%drv2mdl(:)
             call seq_domain_areafactinit(comp(eci)%dom_cc,           &
                  comp(eci)%mdl2drv, comp(eci)%drv2mdl, samegrid, &
                  comp(eci)%mpicom_compid, comp(eci)%iamroot_compid,  &
                  'areafact_'//comp(eci)%oneletterid//'_'//trim(comp(eci)%name))

             ! Area correct component initialization output fields
             call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)

          endif

          ! Map corrected initial component AVs from component to coupler pes
          call seq_map_map(comp(eci)%mapper_cc2x, comp(eci)%c2x_cc, &
               comp(eci)%c2x_cx, msgtag=comp(eci)%cplcompid*100+eci*10+7)

       endif
    enddo

  end subroutine component_init_areacor

  !===============================================================================

  subroutine component_run(Eclock, comp, comp_run, infodata,  &
       seq_flds_x2c_fluxes, seq_flds_c2x_fluxes, &
       comp_prognostic, comp_num, timer_barrier, timer_comp_run, &
       run_barriers, ymd, tod, comp_layout)

    !---------------------------------------------------------------
    ! Description
    ! Run component model
    !
    ! Arguments
    type(ESMF_Clock)     , intent(inout)   :: EClock
    type(component_type) , intent(inout)   :: comp(:)
    interface 
       subroutine comp_run( Eclock, cdata, x2c, c2x)
         use ESMF,          only : ESMF_Clock
         use seq_cdata_mod, only : seq_cdata
         use mct_mod,       only : mct_avect
         implicit none
         type(ESMF_Clock), intent(inout) :: EClock
         type(seq_cdata) , intent(inout) :: cdata
         type(mct_aVect) , intent(inout) :: x2c
         type(mct_aVect) , intent(inout) :: c2x   
       end subroutine comp_run
    end interface 
    type (seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)         , intent(in)           :: seq_flds_x2c_fluxes
    character(len=*)         , intent(in)           :: seq_flds_c2x_fluxes
    logical                  , intent(in)           :: comp_prognostic
    integer                  , intent(in), optional :: comp_num
    character(len=*)         , intent(in), optional :: timer_barrier   
    character(len=*)         , intent(in), optional :: timer_comp_run
    logical                  , intent(in), optional :: run_barriers
    integer                  , intent(in), optional :: ymd  ! Current date (YYYYMMDD)
    integer                  , intent(in), optional :: tod  ! Current time of day (seconds) 
    character(len=*)         , intent(in), optional :: comp_layout
    !
    ! Local Variables
    integer  :: eci
    integer  :: ierr
    integer  :: num_inst
    real(r8) :: time_brun         ! Start time
    real(r8) :: time_erun         ! Ending time
    real(r8) :: cktime            ! delta time
    real(r8) :: cktime_acc(10)    ! cktime accumulator array 1 = all, 2 = atm, etc
    integer  :: cktime_cnt(10)    ! cktime counter array
    logical  :: seq_multi_inst    ! a special case of running multiinstances on the same pes. 
    integer  :: phase, phasemin, phasemax  ! phase support
    logical  :: firstloop         ! first time around phase loop
    !---------------------------------------------------------------

    num_inst = size(comp)
    seq_multi_inst = .false.
    phasemin = 1
    phasemax = 1

    if(present(comp_layout)) then
       if(comp_layout .eq. "sequential" .and. num_inst > 1) then
          seq_multi_inst=.true.
          phasemin = 0
       endif
    endif

    do phase = phasemin,phasemax
       if (phase == phasemin) then
          firstloop = .true.
       else
          firstloop = .false.
       endif
       if (comp(1)%oneletterid == 'a') call seq_infodata_putData(infodata, atm_phase=phase)
       if (comp(1)%oneletterid == 'l') call seq_infodata_putData(infodata, lnd_phase=phase)
       if (comp(1)%oneletterid == 'i') call seq_infodata_putData(infodata, ice_phase=phase)
       if (comp(1)%oneletterid == 'o') call seq_infodata_putData(infodata, ocn_phase=phase)
       if (comp(1)%oneletterid == 'r') call seq_infodata_putData(infodata, rof_phase=phase)
       if (comp(1)%oneletterid == 'g') call seq_infodata_putData(infodata, glc_phase=phase)
       if (comp(1)%oneletterid == 'w') call seq_infodata_putData(infodata, wav_phase=phase)

       do eci = 1,num_inst
          if (comp(eci)%iamin_compid) then

             if (present(timer_barrier))  then
                if (present(run_barriers)) then
                   if (run_barriers) then
                      call t_drvstartf (trim(timer_barrier))
                      call mpi_barrier(comp(eci)%mpicom_compid, ierr) 
                      call t_drvstopf (trim(timer_barrier))
                      time_brun = mpi_wtime()
                   endif
                end if
             end if

             if (present(timer_comp_run)) then
                call t_drvstartf (trim(timer_comp_run), barrier=comp(eci)%mpicom_compid)
             end if
             if (drv_threading) call seq_comm_setnthreads(comp(1)%nthreads_compid) 

             if (comp_prognostic .and. firstloop) then
                call mct_avect_vecmult(comp(eci)%x2c_cc, comp(eci)%drv2mdl, seq_flds_x2c_fluxes, mask_spval=.true.)
             end if

             call comp_run(EClock, comp(eci)%cdata_cc, comp(eci)%x2c_cc, comp(eci)%c2x_cc)

             if (phase == 1) then
                call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)
             endif

             if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

             if (present(timer_comp_run)) then
                call t_drvstopf (trim(timer_comp_run))
             end if

             if (present(comp_num)) then
                if (present(run_barriers)) then
                   if (run_barriers) then
                      time_erun = mpi_wtime()
                      cktime = time_erun - time_brun
                      cktime_acc(comp_num) = cktime_acc(comp_num) + cktime
                      cktime_cnt(comp_num) = cktime_cnt(comp_num) + 1
                      if (present(ymd) .and. present(tod)) then
                         write(logunit,107) ' rstamp ',trim(comp(eci)%name),          &
                              '_run_time: model date = ',ymd,tod,                     &
                              ' avg dt = ',cktime_acc(comp_num)/cktime_cnt(comp_num), &
                              ' dt = ',cktime, ' phase = ',phase
                      end if
                   endif
                end if
             end if

          endif
       enddo   ! eci

    enddo   ! phase

107 format( 3A, 2i8, A, f12.4, A, f12.4 )

  end subroutine component_run

  !===============================================================================


  !===============================================================================

  subroutine component_final(Eclock, comp, comp_final)

    !---------------------------------------------------------------
    ! Description
    ! Run component model
    !
    ! Arguments
    type(ESMF_Clock)     , intent(inout) :: EClock
    type(component_type) , intent(inout) :: comp(:)
    interface
       subroutine comp_final( Eclock, cdata, x2c, c2x)
         use ESMF,          only : ESMF_Clock
         use seq_cdata_mod, only : seq_cdata
         use mct_mod,       only : mct_avect
         implicit none
         type(ESMF_Clock), intent(inout) :: EClock
         type(seq_cdata) , intent(inout) :: cdata
         type(mct_aVect) , intent(inout) :: x2c
         type(mct_aVect) , intent(inout) :: c2x   
       end subroutine comp_final
    end interface 
    !
    ! Local Variables
    integer :: eci
    integer :: num_inst
    !---------------------------------------------------------------

    num_inst = size(comp)
    do eci = 1,num_inst
       if (comp(eci)%iamin_compid) then
          if (drv_threading) call seq_comm_setnthreads(comp(1)%nthreads_compid) 
          call comp_final(EClock, comp(eci)%cdata_cc, comp(eci)%x2c_cc, comp(eci)%c2x_cc)
          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       end if
    end do
    
  end subroutine component_final

  !===============================================================================


  !===============================================================================

  subroutine component_exch(comp, flow, infodata, infodata_string, &
       mpicom_barrier, run_barriers, &
       timer_barrier, timer_comp_exch, timer_map_exch, timer_infodata_exch) 

    !---------------------------------------------------------------
    ! Description
    ! Map x2m_mx to x2m_mm (component input av from 
    ! coupler processes to component model processes)
    !
    ! Arguments
    implicit none
    type(component_type)    , intent(inout)        :: comp(:)
    character(len=3)        , intent(in)           :: flow
    type(seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)        , intent(in)           :: infodata_string
    integer                 , intent(in), optional :: mpicom_barrier      ! mpicom for barrier call
    logical                 , intent(in), optional :: run_barriers
    character(len=*)        , intent(in), optional :: timer_barrier       ! timer
    character(len=*)        , intent(in), optional :: timer_comp_exch
    character(len=*)        , intent(in), optional :: timer_map_exch
    character(len=*)        , intent(in), optional :: timer_infodata_exch
    !
    ! Local Variables
    integer :: eci
    integer :: ierr
    !---------------------------------------------------------------

    if (present(timer_barrier))  then
       if (run_barriers) then
          call t_drvstartf (trim(timer_barrier))
          call mpi_barrier(comp(1)%mpicom_cplallcompid,ierr)
          call t_drvstopf (trim(timer_barrier))
       endif
    end if

    if (present(timer_comp_exch)) then
       if (present(mpicom_barrier)) then
          call t_drvstartf (trim(timer_comp_exch), cplcom=.true., barrier=mpicom_barrier)
       end if
    end if

    do eci = 1,size(comp)
       if (comp(eci)%iamin_cplcompid) then
          if (present(timer_map_exch)) then
             call t_drvstartf (trim(timer_map_exch), barrier=comp(eci)%mpicom_cplcompid)
          end if

          if (flow == 'x2c') then ! coupler to component
             call seq_map_map(comp(eci)%mapper_Cx2c, comp(eci)%x2c_cx, comp(eci)%x2c_cc, &
                  msgtag=comp(eci)%cplcompid*100+eci*10+2)
          else if (flow == 'c2x') then ! component to coupler
             call seq_map_map(comp(eci)%mapper_Cc2x, comp(eci)%c2x_cc, comp(eci)%c2x_cx, &
                  msgtag=comp(eci)%cplcompid*100+eci*10+4)
          end if

          if (present(timer_map_exch)) then
             call t_drvstopf (trim(timer_map_exch))
          end if
       endif
    enddo

    if (present(timer_infodata_exch)) then
       call t_drvstartf (trim(timer_infodata_exch), barrier=mpicom_barrier)
    end if
    if (comp(1)%iamin_cplcompid) then
       call seq_infodata_exchange(infodata, comp(1)%cplcompid, trim(infodata_string))
    end if
    if (present(timer_infodata_exch)) then
       call t_drvstopf (trim(timer_infodata_exch))
    end if

    if (present(timer_comp_exch)) then
       if (present(mpicom_barrier)) then
          call t_drvstopf (trim(timer_comp_exch), cplcom=.true.) 
       end if
    end if

  end subroutine component_exch

  !===============================================================================

  subroutine component_diag(infodata, comp, flow, comment, info_debug, timer_diag ) 

    !---------------------------------------------------------------
    ! Description
    ! Component diagnostics for send/recv to coupler
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout)        :: infodata
    type(component_type)     , intent(in)           :: comp(:)
    character(len=3)         , intent(in)           :: flow 
    character(len=*)         , intent(in)           :: comment
    integer                  , intent(in)           :: info_debug
    character(len=*)         , intent(in), optional :: timer_diag
    !
    ! Local Variables
    integer :: eci
    !---------------------------------------------------------------

    if (info_debug > 1) then
       if (present(timer_diag)) then
          call t_drvstartf (trim(timer_diag), barrier=mpicom_CPLID)
       end if

       do eci = 1,size(comp)
          if (flow == 'x2c') then  ! coupler to component
             call seq_diag_avect_mct(infodata, CPLID, comp(eci)%x2c_cx, &
                  comp(eci)%dom_cx, comp(eci)%gsmap_cx, trim(comment)//comp(eci)%suffix)
          end if
          if (flow == 'c2x') then  ! component to coupler
             call seq_diag_avect_mct(infodata, CPLID, comp(eci)%c2x_cx, &
                  comp(eci)%dom_cx, comp(eci)%gsmap_cx, trim(comment)//comp(eci)%suffix)
          end if
       enddo

       if (present(timer_diag)) then
          call t_drvstopf (trim(timer_diag))
       end if
    endif

  end subroutine component_diag

end module component_mod


