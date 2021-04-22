module seq_domain_mct

  use shr_kind_mod, only: R8=>shr_kind_r8, IN=>shr_kind_in
  use shr_kind_mod, only: CL=>shr_kind_cl
  use shr_sys_mod,  only: shr_sys_flush, shr_sys_abort
  use shr_mpi_mod,  only: shr_mpi_min, shr_mpi_max

  use mct_mod
  use seq_comm_mct
  use seq_infodata_mod
  use seq_map_mod     , only: seq_map_map
  use seq_map_type_mod, only: seq_map

  use component_type_mod

  implicit none
  private ! except
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

  public :: seq_domain_check
  public :: seq_domain_compare
  public :: seq_domain_areafactinit

!--------------------------------------------------------------------------
! Public variables
!--------------------------------------------------------------------------

  real(R8), parameter :: eps_tiny   = 1.0e-16_R8 ! roundoff eps
  real(R8), parameter :: eps_big    = 1.0e+02_R8 ! big eps
  real(R8), parameter :: eps_frac_samegrid = 1.0e-14_R8 ! epsilon for fractions for samegrid

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: seq_domain_check_grid

!================================================================================
contains
!================================================================================

!================================================================================

  subroutine seq_domain_check( infodata, &
       atm, ice, lnd, ocn, rof, glc, &
       samegrid_al, samegrid_ao, samegrid_ro)

    !-----------------------------------------------------------
    ! Uses
    !
    use prep_atm_mod, only: prep_atm_get_mapper_Fi2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fl2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fo2a
    use prep_lnd_mod, only: prep_lnd_get_mapper_Fa2l
    use prep_ocn_mod, only: prep_ocn_get_mapper_SFi2o
    use prep_glc_mod, only: prep_glc_get_mapper_SFl2g
    !
    ! Arguments
    !
    type (seq_infodata_type) , intent(inout) :: infodata
    type(component_type)     , intent(in)    :: atm
    type(component_type)     , intent(in)    :: ice
    type(component_type)     , intent(in)    :: lnd
    type(component_type)     , intent(in)    :: ocn
    type(component_type)     , intent(in)    :: rof
    type(component_type)     , intent(in)    :: glc
    logical                  , intent(in)    :: samegrid_al ! atm lnd grid same
    logical                  , intent(in)    :: samegrid_ao ! atm ocn grid same
    logical                  , intent(in)    :: samegrid_ro ! rof ocn grid same
    !
    ! Local variables
    !
    type(seq_map)   , pointer :: mapper_i2a ! inout needed for lower methods
    type(seq_map)   , pointer :: mapper_i2o ! inout needed for lower methods
    type(seq_map)   , pointer :: mapper_o2a !
    type(seq_map)   , pointer :: mapper_l2g !
    type(seq_map)   , pointer :: mapper_a2l !
    type(seq_map)   , pointer :: mapper_l2a !
                                            !
    type(mct_gGrid) , pointer :: atmdom_a   ! atm domain
    type(mct_gGrid) , pointer :: icedom_i   ! ice domain
    type(mct_gGrid) , pointer :: lnddom_l   ! lnd domain
    type(mct_gGrid) , pointer :: ocndom_o   ! ocn domain
    type(mct_gGrid) , pointer :: glcdom_g   ! glc domain
                                            !
    type(mct_gsMap) , pointer :: gsMap_a    ! atm global seg map 
    type(mct_gsMap) , pointer :: gsMap_i    ! ice global seg map 
    type(mct_gsMap) , pointer :: gsMap_l    ! lnd global seg map 
    type(mct_gsMap) , pointer :: gsMap_o    ! ocn global seg map 
    type(mct_gsMap) , pointer :: gsMap_r    ! ocn global seg map 
    type(mct_gsMap) , pointer :: gsMap_g    ! glc global seg map 
    !
    type(mct_gGrid) :: lnddom_a              ! lnd domain info on atm decomp
    type(mct_gGrid) :: lnddom_g              ! lnd domain info on glc decomp
    type(mct_gGrid) :: icedom_a              ! ice domain info on atm decomp (all grids same)
    type(mct_gGrid) :: ocndom_a              ! ocn domain info on atm decomp (all grids same)
    type(mct_gGrid) :: icedom_o              ! ocn domain info on ocn decomp (atm/ocn grid different)
    !
    real(R8), pointer :: fracl(:)            ! land fraction on atm decomp 
    real(R8), pointer :: fraco(:)            ! ocn  fraction on atm decomp 
    real(R8), pointer :: fraci(:)            ! ice  fraction on atm decomp 
    real(R8), pointer :: maskl(:)            ! land mask on atm decomp (all grids same)
    real(R8), pointer :: maski(:)            ! ice  mask on atm decomp (all grids same)
    real(R8), pointer :: masko(:)            ! ocn  mask on atm decomp (all grids same)
    !
    integer(IN) :: n, kl, ko, ki             ! indicies
    integer(IN) :: k1,k2,k3                  ! indicies
    !
    integer(IN) :: mpicom_cplid
    ! 
    logical      :: atm_present              ! atm present flag
    logical      :: lnd_present              ! lnd present flag
    logical      :: ocn_present              ! ocn present flag
    logical      :: ice_present              ! ice present flag
    logical      :: glc_present              ! glc present flag
    logical      :: rof_present              ! rof present flag
    logical      :: ocnrof_prognostic        ! ocn rof prognostic flag
    integer(IN)  :: rcode                    ! error status
    integer(IN)  :: atmsize                  ! local  size of atm  grid
    integer(IN)  :: lndsize                  ! local  size of land grid
    integer(IN)  :: ocnsize                  ! local  size of ocn  grid
    integer(IN)  :: icesize                  ! local  size of ice  grid
    integer(IN)  :: glcsize                  ! local  size of glc  grid
    integer(IN)  :: gatmsize                 ! global size of atm  grid
    integer(IN)  :: glndsize                 ! global size of land grid
    integer(IN)  :: gocnsize                 ! global size of ocn  grid
    integer(IN)  :: grofsize                 ! global size of ocn  grid
    integer(IN)  :: gicesize                 ! global size of ice  grid
    integer(IN)  :: gglcsize                 ! global size of glc  grid
    integer(IN)  :: npts                     ! local size temporary
    integer(IN)  :: ier                      ! error code
    real(R8)     :: diff,dmaxo,dmaxi         ! difference tracker
    logical      :: iamroot                  ! local masterproc
    real(R8)     :: eps_frac                 ! epsilon for fractions
    real(R8)     :: eps_axmask               ! epsilon for masks, atm/lnd
    real(R8)     :: eps_axgrid               ! epsilon for grid coords, atm/lnd
    real(R8)     :: eps_axarea               ! epsilon for areas, atm/lnd
    real(R8)     :: eps_oimask               ! epsilon for masks, ocn/ice
    real(R8)     :: eps_oigrid               ! epsilon for grid coords, ocn/ice
    real(R8)     :: eps_oiarea               ! epsilon for areas, ocn/ice
    real(R8)     :: my_eps_frac              ! local eps_frac value
    real(R8)     :: rmin1,rmax1,rmin,rmax    ! local min max computation
    !
    real(R8),allocatable :: mask (:)         ! temporary real vector, domain mask
    !
    character(*),parameter :: F00 = "('(seq_domain_check) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check) ',a,i6,a)"
    character(*),parameter :: F02 = "('(seq_domain_check) ',a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check) '
    !-----------------------------------------------------------

    mapper_i2a => prep_atm_get_mapper_Fi2a()
    mapper_i2o => prep_ocn_get_mapper_SFi2o()
    mapper_o2a => prep_atm_get_mapper_Fo2a()
    mapper_l2g => prep_glc_get_mapper_SFl2g()
    mapper_a2l => prep_lnd_get_mapper_Fa2l()
    mapper_l2a => prep_atm_get_mapper_Fl2a()

    call seq_comm_setptrs(CPLID,iamroot=iamroot, mpicom=mpicom_cplid)

    call seq_infodata_GetData( infodata,      &
         lnd_present=lnd_present,             &
         ocn_present=ocn_present,             &
         ice_present=ice_present,             &
         glc_present=glc_present,             &
         atm_present=atm_present,             &
         rof_present=rof_present,             &
         ocnrof_prognostic=ocnrof_prognostic, &
         eps_frac=eps_frac,                   &
         eps_amask=eps_axmask,                &
         eps_agrid=eps_axgrid,                &
         eps_aarea=eps_axarea,                &
         eps_omask=eps_oimask,                &
         eps_ogrid=eps_oigrid,                &
         eps_oarea=eps_oiarea )

    ! Get info

    gsmap_a  => component_get_gsmap_cx(atm) ! gsmap_ax
    atmdom_a => component_get_dom_cx(atm)   ! dom_ax
    atmsize  = mct_avect_lsize(atmdom_a%data)
    gatmsize = mct_gsMap_gsize(gsMap_a)

    if (atm_present .and. lnd_present) then
       gsmap_l  => component_get_gsmap_cx(lnd) ! gsmap_lx
       lnddom_l => component_get_dom_cx(lnd)   ! dom_lx
       lndsize  = mct_avect_lsize(lnddom_l%data)
       glndsize = mct_gsMap_gsize(gsMap_l) 

       if (samegrid_al .and. gatmsize /= glndsize) then
          write(logunit,*) subname,' error: global atmsize = ',&
               gatmsize,' global lndsize= ',glndsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' atm and lnd grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking land maskfrac ---'
       call seq_domain_check_fracmask(lnddom_l%data)
       call mct_gGrid_init(oGGrid=lnddom_a, iGGrid=lnddom_l, lsize=atmsize)
       call mct_aVect_zero(lnddom_a%data)
       call seq_map_map(mapper_l2a, lnddom_l%data, lnddom_a%data, norm=.false.)
       allocate(maskl(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate maskl')
       allocate(fracl(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate fracl')
       call mct_aVect_exportRAttr(lnddom_a%data, 'mask', maskl, atmsize)
       call mct_aVect_exportRAttr(lnddom_a%data, 'frac', fracl, atmsize)
    endif

    if (atm_present .and. ocn_present) then
       gsmap_o  => component_get_gsmap_cx(ocn) ! gsmap_ox
       ocndom_o => component_get_dom_cx(ocn)   ! dom_ox
       ocnsize  = mct_avect_lsize(ocndom_o%data)
       gocnsize = mct_gsMap_gsize(gsMap_o)

       if (samegrid_ao .and. gatmsize /= gocnsize) then
          write(logunit,*) subname,' error: global atmsize = ',gatmsize,' global ocnsize= ',gocnsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' atm and ocn grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking ocean maskfrac ---'
       call seq_domain_check_fracmask(ocndom_o%data)
       call mct_gGrid_init(oGGrid=ocndom_a, iGGrid=ocndom_o, lsize=atmsize)
       call mct_aVect_zero(ocndom_a%data)
       call seq_map_map(mapper_o2a, ocndom_o%data, ocndom_a%data, norm=.false.)
       allocate(masko(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate masko')
       allocate(fraco(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate fraco')
       call mct_aVect_exportRAttr(ocndom_a%data, 'mask', masko, atmsize)
       if (samegrid_ao) then
          call mct_aVect_exportRattr(ocndom_a%data, 'frac', fraco, atmsize)
       else
          call mct_aVect_exportRattr(ocndom_a%data, 'mask', fraco, atmsize)
       endif
    endif
   
    if (atm_present .and. ice_present) then
       gsmap_i  => component_get_gsmap_cx(ice) ! gsmap_ix
       icedom_i => component_get_dom_cx(ice)   ! dom_ix
       icesize  = mct_avect_lsize(icedom_i%data)
       gicesize = mct_gsMap_gsize(gsMap_i) 

       if (samegrid_ao .and. gatmsize /= gicesize) then
          write(logunit,*) subname,' error: global atmsize = ',&
               gatmsize,' global icesize= ',gicesize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' atm and ice grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking ice maskfrac ---'
       call seq_domain_check_fracmask(icedom_i%data)
       call mct_gGrid_init(oGGrid=icedom_a, iGGrid=icedom_i, lsize=atmsize)
       call mct_aVect_zero(icedom_a%data)
       call seq_map_map(mapper_i2a, icedom_i%data, icedom_a%data, norm=.false.)
       allocate(maski(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate maski')
       allocate(fraci(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate fraci')
       call mct_aVect_exportRAttr(icedom_a%data, 'mask', maski, atmsize)
       if (samegrid_ao) then
          call mct_aVect_exportRattr(icedom_a%data, 'frac', fraci, atmsize)
       else
          call mct_aVect_exportRattr(icedom_a%data, 'mask', fraci, atmsize)
       endif
    endif

    if (lnd_present .and. glc_present) then
       gsmap_l  => component_get_gsmap_cx(lnd) ! gsmap_lx
       lnddom_l => component_get_dom_cx(lnd)   ! dom_lx
       lndsize  = mct_avect_lsize(lnddom_l%data)
       glndsize = mct_gsMap_gsize(gsMap_l) 

       gsmap_g  => component_get_gsmap_cx(glc) ! gsmap_gx
       glcdom_g => component_get_dom_cx(glc)   ! dom_gx
       glcsize  = mct_avect_lsize(glcdom_g%data)
       gglcsize = mct_gsMap_gsize(gsMap_g) 

       if (gglcsize /= glndsize) then
          write(logunit,*) subname,' error: global glcsize = ',gglcsize,' global lndsize= ',glndsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' glc and lnd grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking glc maskfrac ---'
       call seq_domain_check_fracmask(glcdom_g%data)
       if (iamroot) write(logunit,F00) ' --- checking lnd maskfrac ---'
       call seq_domain_check_fracmask(lnddom_l%data)
       call mct_gGrid_init(oGGrid=lnddom_g, iGGrid=lnddom_l, lsize=glcsize)
       call mct_aVect_zero(lnddom_g%data)
       call seq_map_map(mapper_l2g, lnddom_l%data, lnddom_g%data, norm=.false.)
       if (iamroot) write(logunit,F00) ' --- checking glc/lnd domains ---'
       npts = glcsize
       allocate(mask(npts),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate mask')
       call mct_aVect_getRAttr(lnddom_g%data,"mask",mask,rcode)
       where (mask < eps_axmask) mask = 0.0_R8
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'mask', eps=eps_axmask, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=mask)
       deallocate(mask,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate mask')
    endif

    if (ice_present .and. ocn_present) then
       gsmap_i  => component_get_gsmap_cx(ice) ! gsmap_ix
       icedom_i => component_get_dom_cx(ice)   ! dom_ix
       icesize  = mct_avect_lsize(icedom_i%data)
       gicesize = mct_gsMap_gsize(gsMap_i) 

       gsmap_o  => component_get_gsmap_cx(ocn) ! gsmap_ox
       ocndom_o => component_get_dom_cx(ocn)   ! dom_ox
       ocnsize  = mct_avect_lsize(ocndom_o%data)
       gocnsize = mct_gsMap_gsize(gsMap_o)

       if (gocnsize /= gicesize) then
          write(logunit,*) subname,' error: global ocnsize = ',gocnsize,' global icesize= ',gicesize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' ocean and ice grid must have the same global size')
       endif
       call mct_gGrid_init(oGGrid=icedom_o, iGGrid=icedom_i, lsize=ocnsize)
       call mct_aVect_zero(icedom_o%data)
       call seq_map_map(mapper_i2o, icedom_i%data, icedom_o%data, norm=.false.)
    end if

    if (rof_present .and. ocnrof_prognostic .and. samegrid_ro) then
       gsmap_r  => component_get_gsmap_cx(glc) ! gsmap_gx
       grofsize = mct_gsMap_gsize(gsMap_r)

       if (gocnsize /= grofsize) then
          write(logunit,*) subname,' error: global ocnsize = ',gocnsize,' global rofsize= ',grofsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' ocean and rof grid must have the same global size')
       endif
    end if

    !------------------------------------------------------------------------------
    ! Check ice/ocean grid consistency
    !------------------------------------------------------------------------------

     if (ocn_present .and. ice_present) then
!    if (samegrid_oi) then       ! doesn't yet exist

       npts = ocnsize
       allocate(mask(npts),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate mask')

       if (iamroot) write(logunit,F00) ' --- checking ocn/ice domains ---'
       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'mask', eps=eps_oigrid, mpicom=mpicom_cplid)
       call mct_aVect_getRAttr(ocndom_o%data,"mask",mask,rcode)
       where (mask < eps_oimask) mask = 0.0_R8

       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'lat' , eps=eps_oigrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'lon' , eps=eps_oigrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'area', eps=eps_oiarea, mpicom=mpicom_cplid, mask=mask)

       deallocate(mask,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate mask')

!    endif
     endif

    !------------------------------------------------------------------------------
    ! Check atm/lnd grid consistency
    !------------------------------------------------------------------------------

    if (atm_present .and. lnd_present .and. samegrid_al) then
       if (iamroot) write(logunit,F00) ' --- checking atm/land domains ---'
       call seq_domain_check_grid(atmdom_a%data, lnddom_a%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maskl)
       call seq_domain_check_grid(atmdom_a%data, lnddom_a%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maskl)
       call seq_domain_check_grid(atmdom_a%data, lnddom_a%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=maskl)
    endif

    !------------------------------------------------------------------------------
    ! Check atm/ocn and atm/ice grid consistency (if samegrid)
    !------------------------------------------------------------------------------

    if (atm_present .and. ice_present .and. samegrid_ao) then
       if (iamroot) write(logunit,F00) ' --- checking atm/ice domains ---'
       call seq_domain_check_grid(atmdom_a%data, icedom_a%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maski)
       call seq_domain_check_grid(atmdom_a%data, icedom_a%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maski)
       call seq_domain_check_grid(atmdom_a%data, icedom_a%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=maski)
    endif

    if (atm_present .and. ocn_present .and. samegrid_ao) then
       if (iamroot) write(logunit,F00) ' --- checking atm/ocn domains ---'
       call seq_domain_check_grid(atmdom_a%data, ocndom_a%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=masko)
       call seq_domain_check_grid(atmdom_a%data, ocndom_a%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=masko)
       call seq_domain_check_grid(atmdom_a%data, ocndom_a%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=masko)
    endif

    !------------------------------------------------------------------------------
    ! Check consistency of land fraction with ocean mask on grid
    !------------------------------------------------------------------------------

    my_eps_frac = eps_frac
    if (samegrid_ao) my_eps_frac = eps_frac_samegrid
    if (.not. samegrid_al) my_eps_frac = eps_big

    if (iamroot) write(logunit,F00) ' --- checking fractions in domains ---'
    dmaxi = 0.0_R8
    dmaxo = 0.0_R8
    do n = 1,atmsize
       if (atm_present .and. lnd_present .and. ice_present) then
          diff = abs(1._R8 - fracl(n) - fraci(n))
          dmaxi = max(diff,dmaxi)
          if (diff > my_eps_frac) then
             write(logunit,*)'inconsistency between land fraction and sea ice fraction'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraci= ',fraci(n),' sum= ',fracl(n)+fraci(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//' inconsistency between land fraction and sea ice fraction')
          end if
          if ((1._R8-fraci(n)) > eps_frac .and. fracl(n) < eps_tiny) then
             write(logunit,*)'inconsistency between land mask and sea ice mask'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraci= ',fraci(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'  inconsistency between land mask and sea ice mask')
          end if
       endif
       if (atm_present .and. lnd_present .and. ocn_present) then
          diff = abs(1._R8 - fracl(n) - fraco(n))
          dmaxo = max(diff,dmaxo)
          if (diff > my_eps_frac) then
             write(logunit,*)'inconsistency between land fraction and ocn land fraction'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraco= ',fraco(n),' sum= ',fracl(n)+fraco(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'  inconsistency between land fraction and ocn land fraction')
          end if
          if ((1._R8-fraco(n)) > eps_frac .and. fracl(n) < eps_tiny) then
             write(logunit,*)'inconsistency between land mask and ocn land mask'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraco= ',fraco(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'  inconsistency between land mask and ocn land mask')
          end if
       endif
    end do 
    if (iamroot) then
       write(logunit,F02) ' maximum           difference for ofrac sum ',dmaxo
       write(logunit,F02) ' maximum           difference for ifrac sum ',dmaxi
       write(logunit,F02) ' maximum allowable difference for  frac sum ',my_eps_frac
       write(logunit,F02) ' maximum allowable tolerance for valid frac ',eps_frac
       call shr_sys_flush(logunit)
    endif

    !------------------------------------------------------------------------------
    ! Clean up allocated memory
    !------------------------------------------------------------------------------

    if (atm_present .and. lnd_present) then
       deallocate(fracl,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate fracl')
       deallocate(maskl,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate maskl')
       call mct_gGrid_clean(lnddom_a, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean lnddom_a')
    endif

    if (atm_present .and. ocn_present) then
       deallocate(fraco,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate fraco')
       deallocate(masko,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate masko')
       call mct_gGrid_clean(ocndom_a, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean ocndom_a')
    endif

    if (atm_present .and. ice_present) then
       deallocate(fraci,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate fraci')
       deallocate(maski,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate maski')
       call mct_gGrid_clean(icedom_a, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean icedom_o')
    endif

    if (ocn_present .and. ice_present) then
       call mct_gGrid_clean(icedom_o, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean icedom_o')
    endif

    call shr_sys_flush(logunit)

  end subroutine seq_domain_check

!===============================================================================
  
! Subprogram not used   subroutine seq_domain_compare(dom1, dom2, mpicom, eps)
! Subprogram not used    
! Subprogram not used     !-----------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Arguments
! Subprogram not used 
! Subprogram not used     type(mct_gGrid)  , intent(in) :: dom1
! Subprogram not used     type(mct_gGrid)  , intent(in) :: dom2
! Subprogram not used     integer(IN)      , intent(in) :: mpicom
! Subprogram not used     real(R8),optional, intent(in) :: eps    ! error condition for compare
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     real(R8) :: leps
! Subprogram not used     character(*),parameter :: F00 = "('(seq_domain_compare) ',4a)"
! Subprogram not used     character(*),parameter :: F01 = "('(seq_domain_compare) ',a,i12,a)"
! Subprogram not used     character(*),parameter :: F02 = "('(seq_domain_compare) ',2a,g23.15)"
! Subprogram not used     character(*),parameter :: F0R = "('(seq_domain_compare) ',2A,2g23.15,A )"
! Subprogram not used     character(*),parameter :: subName = '(seq_domain_compare) '
! Subprogram not used 
! Subprogram not used     leps = eps_tiny
! Subprogram not used     if (present(eps)) then
! Subprogram not used        leps = eps
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call seq_domain_check_grid(dom1%data, dom2%data, 'mask', eps=leps, mpicom=mpicom)
! Subprogram not used     call seq_domain_check_grid(dom1%data, dom2%data, 'lat' , eps=leps, mpicom=mpicom)
! Subprogram not used     call seq_domain_check_grid(dom1%data, dom2%data, 'lon' , eps=leps, mpicom=mpicom)
! Subprogram not used     call seq_domain_check_grid(dom1%data, dom2%data, 'area', eps=leps, mpicom=mpicom)
! Subprogram not used 
! Subprogram not used   end subroutine seq_domain_compare

!===============================================================================
  
  subroutine seq_domain_check_fracmask(dom1)
   
    !-----------------------------------------------------------

    ! Arguments

    type(mct_aVect) , intent(in) :: dom1

    ! Local variables
    integer(in) :: n,npts,ndiff
    integer(in) :: rcode
    real(R8), pointer :: dmask(:)           ! temporaries
    real(R8), pointer :: dfrac(:)           ! temporaries

    character(*),parameter :: F00 = "('(seq_domain_check_fracmask) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check_fracmask) ',a,i12,a)"
    character(*),parameter :: F02 = "('(seq_domain_check_fracmask) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check_fracmask) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check_fracmask) '
    !-----------------------------------------------------------

    npts = mct_aVect_lsize(dom1)

    allocate(dmask(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate dmask')
    allocate(dfrac(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate dfrac')

    call mct_aVect_exportRAttr(dom1, 'mask', dmask, npts)
    call mct_aVect_exportRAttr(dom1, 'frac', dfrac, npts)

    ndiff = 0
    do n = 1,npts
       if (abs(dfrac(n)) > eps_tiny .and. abs(dmask(n)) < eps_tiny) then
!debug            write(logunit,*)'n= ',n,' dfrac= ',dfrac(n),' dmask= ',dmask(n)
          ndiff = ndiff + 1
       endif
    enddo

    if (ndiff > 0) then
       write(logunit,*) trim(subname)," ERROR: incompatible domain mask and frac values"
       call shr_sys_flush(logunit)
       call shr_sys_abort(subName//" incompatible domain mask and frac values")
    endif

    deallocate(dmask,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate dmask')
    deallocate(dfrac,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate dfrac')

 end subroutine seq_domain_check_fracmask

!===============================================================================
  
  subroutine seq_domain_check_grid(dom1, dom2, attr, eps, mpicom, mask)
   
    !-----------------------------------------------------------

    ! Arguments

    type(mct_aVect) , intent(in) :: dom1
    type(mct_aVect) , intent(in) :: dom2
    character(len=*), intent(in) :: attr   ! grid attribute to compare
    real(R8)        , intent(in) :: eps    ! error condition for compare
    integer(IN)     , intent(in) :: mpicom
    real(R8)        , intent(in), optional :: mask(:)

    ! Local variables

    integer(in)       :: n,ndiff            ! indices
    integer(in)       :: npts1,npts2,npts   ! counters
    integer(in)       :: rcode              ! error code
    real(R8)          :: diff,max_diff      ! temporaries
    real(R8)          :: tot_diff           ! maximum diff across all pes
    integer(IN)       :: ier                ! error code
    real(R8), pointer :: data1(:)           ! temporaries
    real(R8), pointer :: data2(:)           ! temporaries
    real(R8), pointer :: lmask(:)           ! temporaries
    logical           :: iamroot            ! local masterproc

    character(*),parameter :: F00 = "('(seq_domain_check_grid) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check_grid) ',a,i12,a)"
    character(*),parameter :: F02 = "('(seq_domain_check_grid) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check_grid) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check_grid) '
    !-----------------------------------------------------------

    call seq_comm_setptrs(CPLID,iamroot=iamroot)

    npts1 = mct_aVect_lsize(dom1)
    npts2 = mct_aVect_lsize(dom2)
    npts  = npts1

    if (npts1 == npts2) then
       if (iamroot) write(logunit,F01) " the domain size is = ", npts
    else
       write(logunit,*) trim(subname)," domain size #1 = ", npts1
       write(logunit,*) trim(subname)," domain size #2 = ", npts2
       write(logunit,*) trim(subname)," ERROR: domain size mis-match"
       call shr_sys_abort(subName//" ERROR: domain size mis-match")
    end if

    allocate(data1(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate data1')
    allocate(data2(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate data2')
    allocate(lmask(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate lmask')

    call mct_aVect_exportRAttr(dom1, trim(attr), data1, npts)
    call mct_aVect_exportRAttr(dom2, trim(attr), data2, npts)
    lmask = 1.0_R8
    if (present(mask)) then
       if (size(mask) /= npts) then
	  call shr_sys_abort(subName//" ERROR: mask size mis-match")
       endif
       lmask = mask
    endif

    ! --- adjust lons to address wraparound issues, we're assuming degree here! ---

    if (trim(attr) == "lon") then
       do n = 1,npts
	  if (data2(n) > data1(n)) then
	     do while ( (data1(n)+360.0_R8) < (data2(n)+180.0_R8) ) ! longitude is periodic
		data1(n) = data1(n) + 360.0_R8
	     end do
	  else
	     do while ( (data2(n)+360.0_R8) < (data1(n)+180.0_R8) ) ! longitude is periodic
		data2(n) = data2(n) + 360.0_R8
	     end do
	  endif
       enddo
    endif

    ! Only check consistency where mask is greater than zero, if mask is present

    max_diff = 0.0_R8
    ndiff = 0
    do n=1,npts
       if (lmask(n) > eps_tiny) then
	  diff = abs(data1(n)-data2(n))
	  max_diff = max(max_diff,diff)
	  if (diff > eps) then
      !debug            write(logunit,*)'n= ',n,' data1= ',data1(n),' data2= ',data2(n),' diff= ',diff, ' eps= ',eps
	     ndiff = ndiff + 1
	  endif
       end if
    end do

    call mpi_reduce(max_diff,tot_diff,1,MPI_REAL8,MPI_MAX,0,mpicom,ier)
    if (iamroot) then
       write(logunit,F02) " maximum           difference for ",trim(attr),tot_diff
       write(logunit,F02) " maximum allowable difference for ",trim(attr),eps
       call shr_sys_flush(logunit)
    endif
    call mpi_barrier(mpicom,ier)

    if (ndiff > 0) then
       write(logunit,*) trim(subname)," ERROR: incompatible domain grid coordinates"
       call shr_sys_flush(logunit)
       call shr_sys_abort(subName//" incompatible domain grid coordinates")
    endif

    deallocate(data1,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate data1')
    deallocate(data2,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate data2')
    deallocate(lmask,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate lmask')

  end subroutine seq_domain_check_grid

!===============================================================================

  subroutine seq_domain_areafactinit(domain, mdl2drv, drv2mdl, &
       samegrid, mpicom, iamroot, comment)
    !-----------------------------------------------------------
    !
    ! Arguments
    !
    type(mct_gGrid)  , pointer             :: domain     ! component domain on component pes
    real(R8)         , pointer             :: mdl2drv(:) ! comp->cpl factor on component pes 
    real(R8)         , pointer             :: drv2mdl(:) ! cpl->comp factor on component pes
    logical          , intent(in)          :: samegrid   ! true => two grids are same
    integer          , intent(in)          :: mpicom     ! mpi communicator on component pes  
    logical          , intent(in)          :: iamroot
    character(len=*) , optional,intent(in) :: comment
    !
    ! Local variables
    !
    integer                :: j1,j2,m1,n,rcode
    integer                :: gridsize,m2dsize,d2msize
    real(R8)               :: rmin1,rmax1,rmin,rmax
    real(R8)               :: rmask,rarea,raream
    character(cl)          :: lcomment
    character(len=*),parameter :: subName = '(seq_domain_areafactinit) '
    character(len=*),parameter :: F0R = "(2A,2g23.15,A )"
    !
    !-----------------------------------------------------------

    lcomment = ''
    if (present(comment)) lcomment = comment

    ! get sizes

    gridsize = mct_gGrid_lsize(domain)
    allocate(drv2mdl(gridsize),mdl2drv(gridsize),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate area correction factors')

    j1 = mct_gGrid_indexRA(domain,"area"    ,dieWith=subName)
    j2 = mct_gGrid_indexRA(domain,"aream"   ,dieWith=subName)
    m1 = mct_gGrid_indexRA(domain,"mask"    ,dieWith=subName)

    mdl2drv(:)=1.0_R8
    drv2mdl(:)=1.0_R8

    if (samegrid) then
       ! default 1.0
    else
       do n=1,gridsize
          rmask  = domain%data%rAttr(m1,n)
          rarea  = domain%data%rAttr(j1,n)
          raream = domain%data%rAttr(j2,n)
          if ( abs(rmask) >= 1.0e-06) then
             if (rarea * raream /= 0.0_R8) then
                mdl2drv(n) = rarea/raream
                drv2mdl(n) = 1.0_R8/mdl2drv(n)
                !if (mdl2drv(n) > 10.0 .or. mdl2drv(n) < 0.1) then
                !   write(logunit,*) trim(subname),' WARNING area,aream= ', &
                !      domain%data%rAttr(j1,n),domain%data%rAttr(j2,n),' in ',n,gridsize
                !endif
             else
                write(logunit,*) trim(subname),' ERROR area,aream= ', &
                     rarea,raream,' in ',n,gridsize
                call shr_sys_flush(logunit)
                call shr_sys_abort()
             endif
          endif
       enddo
    end if

    rmin1 = minval(mdl2drv)
    rmax1 = maxval(mdl2drv)
    call shr_mpi_min(rmin1,rmin,mpicom)
    call shr_mpi_max(rmax1,rmax,mpicom)
    if (iamroot) write(logunit,F0R) trim(subname),' : min/max mdl2drv ',rmin,rmax,trim(lcomment)

    rmin1 = minval(drv2mdl)
    rmax1 = maxval(drv2mdl)
    call shr_mpi_min(rmin1,rmin,mpicom)
    call shr_mpi_max(rmax1,rmax,mpicom)
    if (iamroot) write(logunit,F0R) trim(subname),' : min/max drv2mdl ',rmin,rmax,trim(lcomment)
    if (iamroot) call shr_sys_flush(logunit)

  end subroutine seq_domain_areafactinit

!===============================================================================

end module seq_domain_mct



