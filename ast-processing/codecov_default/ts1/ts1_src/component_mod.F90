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
!    Copyright (C) Silicon Graphics International Corp.
!    All rights reserved.
!
!    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
!    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
!    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
!    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
!    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
!    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
!    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
!    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
!    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
!    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
!    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
!    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
!    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
!    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
!
! Copyright Notice
!  + 1993 University of Chicago
!  + 1993 Mississippi State University

!
!    Copyright (C) Silicon Graphics International Corp.
!    All rights reserved.
!
!    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
!    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
!    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
!    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
!    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
!    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
!    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
!    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
!    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
!    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
!    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
!    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
!    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
!    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
!

        integer MPI_VERSION
        integer MPI_SUBVERSION
        logical MPI_SUBARRAYS_SUPPORTED
        logical MPI_ASYNC_PROTECTS_NONBLOCKING

        parameter (MPI_VERSION                    = 3)
        parameter (MPI_SUBVERSION                 = 1)
        parameter (MPI_SUBARRAYS_SUPPORTED        = .FALSE.)
        parameter (MPI_ASYNC_PROTECTS_NONBLOCKING = .FALSE.)

! MPI_Status

        integer MPI_STATUS_SIZE
        parameter (MPI_STATUS_SIZE      = 6)

! Misc Fortran declarations

        integer MPI_BOTTOM
        common /MPI_SGI_PRIVATE/ MPI_BOTTOM

        external MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN, MPI_DUP_FN

        integer MPI_IN_PLACE(1)
        common /MPI_SGI_PRIVATE_INPLACE/ MPI_IN_PLACE

! MPI-2 Section 5.3

        CHARACTER*(1) MPI_ARGV_NULL(1)
        common /MPI_SGI_PRIVATE_CHAR/ MPI_ARGV_NULL

        CHARACTER*(1) MPI_ARGVS_NULL(1)
        EQUIVALENCE(MPI_ARGV_NULL,MPI_ARGVS_NULL(1))

        integer MPI_ERRCODES_IGNORE(1)
        EQUIVALENCE(MPI_BOTTOM,MPI_ERRCODES_IGNORE(1))

! MPI-1 error codes and classes

        integer MPI_SUCCESS
        integer MPI_ERR_BUFFER
        integer MPI_ERR_COUNT
        integer MPI_ERR_TYPE
        integer MPI_ERR_TAG
        integer MPI_ERR_COMM
        integer MPI_ERR_RANK
        integer MPI_ERR_REQUEST
        integer MPI_ERR_ROOT
        integer MPI_ERR_GROUP
        integer MPI_ERR_OP
        integer MPI_ERR_TOPOLOGY
        integer MPI_ERR_DIMS
        integer MPI_ERR_ARG
        integer MPI_ERR_UNKNOWN
        integer MPI_ERR_TRUNCATE
        integer MPI_ERR_OTHER
        integer MPI_ERR_INTERN
        integer MPI_ERR_IN_STATUS
        integer MPI_ERR_PENDING

        parameter (MPI_SUCCESS          = 0)
        parameter (MPI_ERR_BUFFER       = 1)
        parameter (MPI_ERR_COUNT        = 2)
        parameter (MPI_ERR_TYPE         = 3)
        parameter (MPI_ERR_TAG          = 4)
        parameter (MPI_ERR_COMM         = 5)
        parameter (MPI_ERR_RANK         = 6)
        parameter (MPI_ERR_REQUEST      = 7)
        parameter (MPI_ERR_ROOT         = 8)
        parameter (MPI_ERR_GROUP        = 9)
        parameter (MPI_ERR_OP           = 10)
        parameter (MPI_ERR_TOPOLOGY     = 11)
        parameter (MPI_ERR_DIMS         = 12)
        parameter (MPI_ERR_ARG          = 13)
        parameter (MPI_ERR_UNKNOWN      = 14)
        parameter (MPI_ERR_TRUNCATE     = 15)
        parameter (MPI_ERR_OTHER        = 16)
        parameter (MPI_ERR_INTERN       = 17)
        parameter (MPI_ERR_IN_STATUS    = 18)
        parameter (MPI_ERR_PENDING      = 19)

! MPI-2 error codes and classes

        integer MPI_ERR_ACCESS
        integer MPI_ERR_AMODE
        integer MPI_ERR_ASSERT
        integer MPI_ERR_BAD_FILE
        integer MPI_ERR_BASE
        integer MPI_ERR_CONVERSION
        integer MPI_ERR_DISP
        integer MPI_ERR_DUP_DATAREP
        integer MPI_ERR_FILE_EXISTS
        integer MPI_ERR_FILE_IN_USE
        integer MPI_ERR_FILE
        integer MPI_ERR_INFO_KEY
        integer MPI_ERR_INFO_NOKEY
        integer MPI_ERR_INFO_VALUE
        integer MPI_ERR_INFO
        integer MPI_ERR_IO
        integer MPI_ERR_KEYVAL
        integer MPI_ERR_LOCKTYPE
        integer MPI_ERR_NAME
        integer MPI_ERR_NO_MEM
        integer MPI_ERR_NOT_SAME
        integer MPI_ERR_NO_SPACE
        integer MPI_ERR_NO_SUCH_FILE
        integer MPI_ERR_PORT
        integer MPI_ERR_QUOTA
        integer MPI_ERR_READ_ONLY
        integer MPI_ERR_RMA_CONFLICT
        integer MPI_ERR_RMA_SYNC
        integer MPI_ERR_SERVICE
        integer MPI_ERR_SIZE
        integer MPI_ERR_SPAWN
        integer MPI_ERR_UNSUPPORTED_DATAREP
        integer MPI_ERR_UNSUPPORTED_OPERATION
        integer MPI_ERR_WIN
        integer MPI_ERR_RMA_RANGE
        integer MPI_ERR_RMA_ATTACH
        integer MPI_ERR_RMA_SHARED
        integer MPI_ERR_RMA_FLAVOR
        integer MPI_T_ERR_CANNOT_INIT
        integer MPI_T_ERR_NOT_INITIALIZED
        integer MPI_T_ERR_MEMORY
        integer MPI_T_ERR_INVALID_INDEX
        integer MPI_T_ERR_INVALID_ITEM
        integer MPI_T_ERR_INVALID_SESSION
        integer MPI_T_ERR_INVALID_HANDLE
        integer MPI_T_ERR_OUT_OF_HANDLES
        integer MPI_T_ERR_OUT_OF_SESSIONS
        integer MPI_T_ERR_CVAR_SET_NOT_NOW
        integer MPI_T_ERR_CVAR_SET_NEVER
        integer MPI_T_ERR_PVAR_NO_WRITE
        integer MPI_T_ERR_PVAR_NO_STARTSTOP
        integer MPI_T_ERR_PVAR_NO_ATOMIC
        integer MPI_T_ERR_INVALID_NAME
        integer MPI_T_ERR_INVALID

        parameter (MPI_ERR_ACCESS               = 28)
        parameter (MPI_ERR_AMODE                = 29)
        parameter (MPI_ERR_ASSERT               = 30)
        parameter (MPI_ERR_BAD_FILE             = 31)
        parameter (MPI_ERR_BASE                 = 32)
        parameter (MPI_ERR_CONVERSION           = 33)
        parameter (MPI_ERR_DISP                 = 34)
        parameter (MPI_ERR_DUP_DATAREP          = 35)
        parameter (MPI_ERR_FILE_EXISTS          = 36)
        parameter (MPI_ERR_FILE_IN_USE          = 37)
        parameter (MPI_ERR_FILE                 = 38)
        parameter (MPI_ERR_INFO_KEY             = 39)
        parameter (MPI_ERR_INFO_NOKEY           = 40)
        parameter (MPI_ERR_INFO_VALUE           = 41)
        parameter (MPI_ERR_INFO                 = 42)
        parameter (MPI_ERR_IO                   = 43)
        parameter (MPI_ERR_KEYVAL               = 44)
        parameter (MPI_ERR_LOCKTYPE             = 45)
        parameter (MPI_ERR_NAME                 = 46)
        parameter (MPI_ERR_NO_MEM               = 47)
        parameter (MPI_ERR_NOT_SAME             = 48)
        parameter (MPI_ERR_NO_SPACE             = 49)
        parameter (MPI_ERR_NO_SUCH_FILE         = 50)
        parameter (MPI_ERR_PORT                 = 51)
        parameter (MPI_ERR_QUOTA                = 52)
        parameter (MPI_ERR_READ_ONLY            = 53)
        parameter (MPI_ERR_RMA_CONFLICT         = 54)
        parameter (MPI_ERR_RMA_SYNC             = 55)
        parameter (MPI_ERR_SERVICE              = 56)
        parameter (MPI_ERR_SIZE                 = 57)
        parameter (MPI_ERR_SPAWN                = 58)
        parameter (MPI_ERR_UNSUPPORTED_DATAREP  = 59)
        parameter (MPI_ERR_UNSUPPORTED_OPERATION= 60)
        parameter (MPI_ERR_WIN                  = 61)
        parameter (MPI_ERR_RMA_RANGE            = 62)
        parameter (MPI_ERR_RMA_ATTACH           = 63)
        parameter (MPI_ERR_RMA_SHARED           = 64)
        parameter (MPI_ERR_RMA_FLAVOR           = 65)
        parameter (MPI_T_ERR_CANNOT_INIT        = 66)
        parameter (MPI_T_ERR_NOT_INITIALIZED    = 67)
        parameter (MPI_T_ERR_MEMORY             = 68)
        parameter (MPI_T_ERR_INVALID_INDEX      = 69)
        parameter (MPI_T_ERR_INVALID_ITEM       = 70)
        parameter (MPI_T_ERR_INVALID_SESSION    = 71)
        parameter (MPI_T_ERR_INVALID_HANDLE     = 72)
        parameter (MPI_T_ERR_OUT_OF_HANDLES     = 73)
        parameter (MPI_T_ERR_OUT_OF_SESSIONS    = 74)
        parameter (MPI_T_ERR_CVAR_SET_NOT_NOW   = 75)
        parameter (MPI_T_ERR_CVAR_SET_NEVER     = 76)
        parameter (MPI_T_ERR_PVAR_NO_WRITE      = 77)
        parameter (MPI_T_ERR_PVAR_NO_STARTSTOP  = 78)
        parameter (MPI_T_ERR_PVAR_NO_ATOMIC     = 79)
        parameter (MPI_T_ERR_INVALID_NAME       = 80)
        parameter (MPI_T_ERR_INVALID            = 81)

        integer MPI_ERR_LASTCODE
        parameter (MPI_ERR_LASTCODE             = 100)

! Permanent keyvals

        integer MPI_KEYVAL_INVALID
        integer MPI_TAG_UB
        integer MPI_HOST
        integer MPI_IO
        integer MPI_WTIME_IS_GLOBAL
        integer MPI_UNIVERSE_SIZE
        integer MPI_APPNUM
        integer MPI_LASTUSEDCODE

        parameter (MPI_KEYVAL_INVALID   = 0)
        parameter (MPI_TAG_UB           = 5)
        parameter (MPI_HOST             = 6)
        parameter (MPI_IO               = 7)
        parameter (MPI_WTIME_IS_GLOBAL  = 8)
        parameter (MPI_UNIVERSE_SIZE   = 10)
        parameter (MPI_APPNUM          = 12)
        parameter (MPI_LASTUSEDCODE    = 13)

! Results of the compare operations

        integer MPI_IDENT
        integer MPI_CONGRUENT
        integer MPI_SIMILAR
        integer MPI_UNEQUAL

        parameter (MPI_IDENT            = 0)
        parameter (MPI_CONGRUENT        = 1)
        parameter (MPI_SIMILAR          = 2)
        parameter (MPI_UNEQUAL          = 3)

! Topology types

        integer MPI_GRAPH
        integer MPI_CART
        integer MPI_DIST_GRAPH

        parameter (MPI_GRAPH    = 1)
        parameter (MPI_CART     = 2)
        parameter (MPI_DIST_GRAPH = 3)

        integer MPI_UNWEIGHTED
        common /MPI_SGI_PRIVATE_UNWEIGHTED/ MPI_UNWEIGHTED
        integer MPI_WEIGHTS_EMPTY
        common /MPI_SGI_PRIVATE_WEIGHTS_EMPTY/ MPI_WEIGHTS_EMPTY

! Misc constants

        integer MPI_MAX_PROCESSOR_NAME
        parameter (MPI_MAX_PROCESSOR_NAME = 255)

        integer MPI_MAX_ERROR_STRING
        parameter (MPI_MAX_ERROR_STRING = 255)

        integer MPI_MAX_LIBRARY_VERSION_STRING
        parameter (MPI_MAX_LIBRARY_VERSION_STRING = 255)

        integer MPI_MAX_OBJECT_NAME
        parameter (MPI_MAX_OBJECT_NAME = 127)

        integer MPI_BSEND_OVERHEAD
        parameter (MPI_BSEND_OVERHEAD = 32)

        integer MPI_ROOT
        parameter (MPI_ROOT = -4)

        integer MPI_UNDEFINED
        parameter (MPI_UNDEFINED = -3)

        integer MPI_ANY_SOURCE
        parameter (MPI_ANY_SOURCE = -2)

        integer MPI_PROC_NULL
        parameter (MPI_PROC_NULL = -1)

        integer MPI_ANY_TAG
        parameter (MPI_ANY_TAG = -1)

! The following 2 lines are included in the main mpif.h
!       double precision MPI_WTIME, MPI_WTICK
!       external MPI_WTIME, MPI_WTICK

! MPI-2 Section 4.10

        integer MPI_MAX_INFO_KEY
        parameter (MPI_MAX_INFO_KEY = 254)

        integer MPI_MAX_INFO_VAL
        parameter (MPI_MAX_INFO_VAL = 1023)

! MPI-2 Section 5.4

        integer MPI_MAX_PORT_NAME
        parameter (MPI_MAX_PORT_NAME = 255)

! Kind values for MPI-2

        integer MPI_INTEGER_KIND
        parameter (MPI_INTEGER_KIND = 4)

        integer MPI_OFFSET_KIND
        parameter (MPI_OFFSET_KIND = 8)

        integer MPI_ADDRESS_KIND
        parameter (MPI_ADDRESS_KIND = 8)

        integer MPI_COUNT_KIND
        parameter (MPI_COUNT_KIND = 8)

! Section 6.4 bindings for one-sided communication

        integer MPI_MODE_NOCHECK
        integer MPI_MODE_NOSTORE
        integer MPI_MODE_NOPUT
        integer MPI_MODE_NOPRECEDE
        integer MPI_MODE_NOSUCCEED

        parameter (MPI_MODE_NOCHECK             = 1)
        parameter (MPI_MODE_NOSTORE             = 2)
        parameter (MPI_MODE_NOPUT               = 4)
        parameter (MPI_MODE_NOPRECEDE           = 8)
        parameter (MPI_MODE_NOSUCCEED           = 16)

        integer MPI_LOCK_SHARED
        parameter (MPI_LOCK_SHARED              = 1)
        integer MPI_LOCK_EXCLUSIVE
        parameter (MPI_LOCK_EXCLUSIVE           = 2)

! Thread-safety support levels

        integer MPI_THREAD_SINGLE
        integer MPI_THREAD_FUNNELED
        integer MPI_THREAD_SERIALIZED
        integer MPI_THREAD_MULTIPLE

        parameter (MPI_THREAD_SINGLE            = 0)
        parameter (MPI_THREAD_FUNNELED          = 1)
        parameter (MPI_THREAD_SERIALIZED        = 2)
        parameter (MPI_THREAD_MULTIPLE          = 3)

! Datatype Decoding constants

        integer MPI_COMBINER_NAMED
        integer MPI_COMBINER_CONTIGUOUS
        integer MPI_COMBINER_VECTOR
        integer MPI_COMBINER_HVECTOR
        integer MPI_COMBINER_INDEXED
        integer MPI_COMBINER_HINDEXED
        integer MPI_COMBINER_STRUCT
        integer MPI_COMBINER_DARRAY
        integer MPI_COMBINER_DUP
        integer MPI_COMBINER_F90_COMPLEX
        integer MPI_COMBINER_F90_INTEGER
        integer MPI_COMBINER_F90_REAL
        integer MPI_COMBINER_HINDEXED_INTEGER
        integer MPI_COMBINER_HVECTOR_INTEGER
        integer MPI_COMBINER_INDEXED_BLOCK
        integer MPI_COMBINER_RESIZED
        integer MPI_COMBINER_STRUCT_INTEGER
        integer MPI_COMBINER_SUBARRAY
        integer MPI_COMBINER_HINDEXED_BLOCK

        parameter (MPI_COMBINER_NAMED            = (-1))
        parameter (MPI_COMBINER_CONTIGUOUS       = 0)
        parameter (MPI_COMBINER_VECTOR           = 1)
        parameter (MPI_COMBINER_HVECTOR          = 2)
        parameter (MPI_COMBINER_INDEXED          = 3)
        parameter (MPI_COMBINER_HINDEXED         = 4)
        parameter (MPI_COMBINER_STRUCT           = 5)
        parameter (MPI_COMBINER_DARRAY           = 6)
        parameter (MPI_COMBINER_DUP              = 7)
        parameter (MPI_COMBINER_F90_COMPLEX      = 8)
        parameter (MPI_COMBINER_F90_INTEGER      = 9)
        parameter (MPI_COMBINER_F90_REAL         = 10)
        parameter (MPI_COMBINER_HINDEXED_INTEGER = 11)
        parameter (MPI_COMBINER_HVECTOR_INTEGER  = 12)
        parameter (MPI_COMBINER_INDEXED_BLOCK    = 13)
        parameter (MPI_COMBINER_RESIZED          = 14)
        parameter (MPI_COMBINER_STRUCT_INTEGER   = 15)
        parameter (MPI_COMBINER_SUBARRAY         = 16)
        parameter (MPI_COMBINER_HINDEXED_BLOCK   = 17)

! Permanent window keyvals

        integer MPI_WIN_BASE
        integer MPI_WIN_SIZE
        integer MPI_WIN_DISP_UNIT
        integer MPI_WIN_CREATE_FLAVOR
        integer MPI_WIN_MODEL

        parameter (MPI_WIN_BASE         = 4)
        parameter (MPI_WIN_SIZE         = 5)
        parameter (MPI_WIN_DISP_UNIT    = 6)
        parameter (MPI_WIN_CREATE_FLAVOR = 8)
        parameter (MPI_WIN_MODEL        = 10)

! Typeclasses

        integer MPI_TYPECLASS_INTEGER
        integer MPI_TYPECLASS_REAL
        integer MPI_TYPECLASS_COMPLEX

        parameter (MPI_TYPECLASS_INTEGER = 1)
        parameter (MPI_TYPECLASS_REAL    = 2)
        parameter (MPI_TYPECLASS_COMPLEX = 3)

! Communicator types

        integer MPI_COMM_TYPE_SHARED

        parameter (MPI_COMM_TYPE_SHARED = 1)

! Window flavors

        integer MPI_WIN_FLAVOR_CREATE
        integer MPI_WIN_FLAVOR_ALLOCATE
        integer MPI_WIN_FLAVOR_DYNAMIC
        integer MPI_WIN_FLAVOR_SHARED

        parameter (MPI_WIN_FLAVOR_CREATE    = 1)
        parameter (MPI_WIN_FLAVOR_ALLOCATE  = 2)
        parameter (MPI_WIN_FLAVOR_DYNAMIC   = 3)
        parameter (MPI_WIN_FLAVOR_SHARED    = 4)

        integer MPI_WIN_SEPARATE
        integer MPI_WIN_UNIFIED

        parameter (MPI_WIN_SEPARATE = 1)
        parameter (MPI_WIN_UNIFIED  = 2)

! MPI-2 I/O definitions

!
!    Copyright (C) Silicon Graphics International Corp.
!    All rights reserved.
!
!    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
!    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
!    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
!    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
!    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
!    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
!    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
!    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
!    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
!    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
!    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
!    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
!    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
!    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
!
!    Fortran MPI-IO programs 
!    Copyright (C) 1997 University of Chicago. 
!
      INTEGER MPI_MODE_RDONLY, MPI_MODE_RDWR, MPI_MODE_WRONLY
      INTEGER MPI_MODE_DELETE_ON_CLOSE, MPI_MODE_UNIQUE_OPEN
      INTEGER MPI_MODE_CREATE, MPI_MODE_EXCL
      INTEGER MPI_MODE_APPEND, MPI_MODE_SEQUENTIAL
      PARAMETER (MPI_MODE_RDONLY=2, MPI_MODE_RDWR=8, MPI_MODE_WRONLY=4)
      PARAMETER (MPI_MODE_CREATE=1, MPI_MODE_DELETE_ON_CLOSE=16)
      PARAMETER (MPI_MODE_UNIQUE_OPEN=32, MPI_MODE_EXCL=64)
      PARAMETER (MPI_MODE_APPEND=128, MPI_MODE_SEQUENTIAL=256)
!
      INTEGER MPI_FILE_NULL
      PARAMETER (MPI_FILE_NULL=0)
!
      INTEGER MPI_MAX_DATAREP_STRING
      PARAMETER (MPI_MAX_DATAREP_STRING=128)
!
      INTEGER MPI_SEEK_SET, MPI_SEEK_CUR, MPI_SEEK_END
      PARAMETER (MPI_SEEK_SET=600, MPI_SEEK_CUR=602, MPI_SEEK_END=604)
!
      INTEGER MPIO_REQUEST_NULL
      PARAMETER (MPIO_REQUEST_NULL=0)
!
!      INTEGER MPI_OFFSET_KIND
!      PARAMETER (MPI_OFFSET_KIND=8)
!
      integer(kind=8) MPI_DISPLACEMENT_CURRENT
      PARAMETER (MPI_DISPLACEMENT_CURRENT=-54278278)

      INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN
      PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
      INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
      INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
      PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
      PARAMETER (MPI_DISTRIBUTE_NONE=123)
      PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)

      EXTERNAL MPI_CONVERSION_FN_NULL
!
!   End Fortran MPI-IO

!
!    Copyright (C) Silicon Graphics International Corp.
!    All rights reserved.
!
!    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
!    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
!    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
!    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
!    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
!    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
!    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
!    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
!    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
!    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
!    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
!    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
!    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
!    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
!

! MPI_Status

        integer MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
        integer MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
        equivalence (MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE)

        common /MPI_SGI_PRIVATE_STATUS/ MPI_STATUS_IGNORE

! Permanent window keyvals

        external MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN
        external MPI_COMM_DUP_FN
        external MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN
        external MPI_WIN_DUP_FN

!
        integer MPI_FLOAT_INT
        integer MPI_DOUBLE_INT
        integer MPI_LONG_INT
        integer MPI_2INT
        integer MPI_SHORT_INT
        integer MPI_LONG_DOUBLE_INT

        parameter (MPI_FLOAT_INT        = 31)
        parameter (MPI_DOUBLE_INT       = 32)
        parameter (MPI_LONG_INT         = 33)
        parameter (MPI_2INT             = 34)
        parameter (MPI_SHORT_INT        = 35)
        parameter (MPI_LONG_DOUBLE_INT  = 36)

        integer MPI_BYTE
        integer MPI_PACKED
        integer MPI_UB
        integer MPI_LB

        parameter (MPI_BYTE             = 27)
        parameter (MPI_PACKED           = 28)
        parameter (MPI_UB               = 29)
        parameter (MPI_LB               = 30)

        integer MPI_2REAL
        integer MPI_2INTEGER

        parameter (MPI_2REAL            = 37)
        parameter (MPI_2INTEGER         = 39)

        integer MPI_AINT
        integer MPI_OFFSET

        parameter (MPI_AINT             = 55)
        parameter (MPI_OFFSET           = 56)

! MPI_Op

        integer MPI_OP_NULL
        integer MPI_MAX
        integer MPI_MIN
        integer MPI_SUM
        integer MPI_PROD
        integer MPI_LAND
        integer MPI_BAND
        integer MPI_LOR
        integer MPI_BOR
        integer MPI_LXOR
        integer MPI_BXOR
        integer MPI_MAXLOC
        integer MPI_MINLOC
        integer MPI_REPLACE
        integer MPI_NO_OP

        parameter (MPI_OP_NULL  = 0)
        parameter (MPI_MAX      = 1)
        parameter (MPI_MIN      = 2)
        parameter (MPI_SUM      = 3)
        parameter (MPI_PROD     = 4)
        parameter (MPI_LAND     = 5)
        parameter (MPI_BAND     = 6)
        parameter (MPI_LOR      = 7)
        parameter (MPI_BOR      = 8)
        parameter (MPI_LXOR     = 9)
        parameter (MPI_BXOR     = 10)
        parameter (MPI_MAXLOC   = 11)
        parameter (MPI_MINLOC   = 12)
        parameter (MPI_REPLACE  = 13)
        parameter (MPI_NO_OP    = 14)

! MPI_Datatype

        integer MPI_DATATYPE_NULL

        integer MPI_CHAR
        integer MPI_SHORT
        integer MPI_INT
        integer MPI_LONG
        integer MPI_UNSIGNED_CHAR
        integer MPI_UNSIGNED_SHORT
        integer MPI_UNSIGNED
        integer MPI_UNSIGNED_LONG
        integer MPI_FLOAT
        integer MPI_DOUBLE
        integer MPI_LONG_DOUBLE
        integer MPI_LONG_LONG
        integer MPI_LONG_LONG_INT

        integer MPI_INTEGER
        integer MPI_REAL
        integer MPI_DOUBLE_PRECISION
        integer MPI_COMPLEX
        integer MPI_DOUBLE_COMPLEX
        integer MPI_LOGICAL
        integer MPI_CHARACTER
        integer MPI_INTEGER1
        integer MPI_INTEGER2
        integer MPI_INTEGER4
        integer MPI_INTEGER8
        integer MPI_REAL4
        integer MPI_REAL8
        integer MPI_REAL16

        integer MPI_2DOUBLE_PRECISION

        integer MPI_WCHAR
        integer MPI_SIGNED_CHAR
        integer MPI_UNSIGNED_LONG_LONG

        integer MPI_INTEGER16
        integer MPI_COMPLEX8
        integer MPI_COMPLEX16
        integer MPI_COMPLEX32

        integer MPI_INT8_T
        integer MPI_INT16_T
        integer MPI_INT32_T
        integer MPI_INT64_T
        integer MPI_UINT8_T
        integer MPI_UINT16_T
        integer MPI_UINT32_T
        integer MPI_UINT64_T

        integer MPI_C_BOOL
        integer MPI_C_FLOAT_COMPLEX
        integer MPI_C_COMPLEX
        integer MPI_C_DOUBLE_COMPLEX
        integer MPI_C_LONG_DOUBLE_COMPLEX

        integer MPI_COUNT

        integer MPI_CXX_BOOL
        integer MPI_CXX_FLOAT_COMPLEX
        integer MPI_CXX_DOUBLE_COMPLEX
        integer MPI_CXX_LONG_DOUBLE_COMPLEX

        parameter (MPI_DATATYPE_NULL    = 0)

        parameter (MPI_CHAR             = 1)
        parameter (MPI_SHORT            = 2)
        parameter (MPI_INT              = 3)
        parameter (MPI_LONG             = 4)
        parameter (MPI_UNSIGNED_CHAR    = 5)
        parameter (MPI_UNSIGNED_SHORT   = 6)
        parameter (MPI_UNSIGNED         = 7)
        parameter (MPI_UNSIGNED_LONG    = 8)
        parameter (MPI_FLOAT            = 9)
        parameter (MPI_DOUBLE           = 10)
        parameter (MPI_LONG_DOUBLE      = 11)
        parameter (MPI_LONG_LONG        = 12)
        parameter (MPI_LONG_LONG_INT    = 12)

        parameter (MPI_INTEGER          = 13)
        parameter (MPI_REAL             = 14)
        parameter (MPI_DOUBLE_PRECISION = 15)
        parameter (MPI_COMPLEX          = 16)
        parameter (MPI_DOUBLE_COMPLEX   = 17)
        parameter (MPI_LOGICAL          = 18)
        parameter (MPI_CHARACTER        = 19)
        parameter (MPI_INTEGER1         = 20)
        parameter (MPI_INTEGER2         = 21)
        parameter (MPI_INTEGER4         = 22)
        parameter (MPI_INTEGER8         = 23)
        parameter (MPI_REAL4            = 24)
        parameter (MPI_REAL8            = 25)
        parameter (MPI_REAL16           = 26)

        parameter (MPI_2DOUBLE_PRECISION= 38)

        parameter (MPI_WCHAR            = 40)
        parameter (MPI_SIGNED_CHAR      = 41)
        parameter (MPI_UNSIGNED_LONG_LONG = 42)

        parameter (MPI_INTEGER16        = 43)
        parameter (MPI_COMPLEX8         = 44)
        parameter (MPI_COMPLEX16        = 45)
        parameter (MPI_COMPLEX32        = 46)

        parameter (MPI_INT8_T           = 47)
        parameter (MPI_INT16_T          = 48)
        parameter (MPI_INT32_T          = 49)
        parameter (MPI_INT64_T          = 50)
        parameter (MPI_UINT8_T          = 51)
        parameter (MPI_UINT16_T         = 52)
        parameter (MPI_UINT32_T         = 53)
        parameter (MPI_UINT64_T         = 54)

        parameter (MPI_C_BOOL           = 57)
        parameter (MPI_C_FLOAT_COMPLEX  = 58)
        parameter (MPI_C_COMPLEX        = 58)
        parameter (MPI_C_DOUBLE_COMPLEX = 59)
        parameter (MPI_C_LONG_DOUBLE_COMPLEX = 60)

        parameter (MPI_COUNT            = 61)

        parameter (MPI_CXX_BOOL                = 62)
        parameter (MPI_CXX_FLOAT_COMPLEX       = 63)
        parameter (MPI_CXX_DOUBLE_COMPLEX      = 64)
        parameter (MPI_CXX_LONG_DOUBLE_COMPLEX = 65)

! MPI_Comm

        integer MPI_COMM_NULL
        integer MPI_COMM_WORLD
        integer MPI_COMM_SELF

        parameter (MPI_COMM_NULL        = 0)
        parameter (MPI_COMM_WORLD       = 1)
        parameter (MPI_COMM_SELF        = 2)

! MPI_Errhandler

        integer MPI_ERRHANDLER_NULL
        integer MPI_ERRORS_ARE_FATAL
        integer MPI_ERRORS_RETURN

        parameter (MPI_ERRHANDLER_NULL  = 0)
        parameter (MPI_ERRORS_ARE_FATAL = 1)
        parameter (MPI_ERRORS_RETURN    = 2)

        integer MPI_SOURCE
        integer MPI_TAG
        integer MPI_ERROR

        parameter (MPI_SOURCE           = 1)
        parameter (MPI_TAG              = 2)
        parameter (MPI_ERROR            = 3)

! MPI_Group

        integer MPI_GROUP_NULL
        integer MPI_GROUP_EMPTY

        parameter (MPI_GROUP_NULL  = 0)
        parameter (MPI_GROUP_EMPTY = 1)

        integer MPI_INFO_NULL
        parameter (MPI_INFO_NULL = 0)

        integer MPI_INFO_ENV
        parameter (MPI_INFO_ENV = 1)

! MPI_Message

        integer MPI_MESSAGE_NO_PROC
        integer MPI_MESSAGE_NULL

        parameter (MPI_MESSAGE_NO_PROC = -1)
        parameter (MPI_MESSAGE_NULL = 0)

! MPI_Request

        integer MPI_REQUEST_NULL
        parameter (MPI_REQUEST_NULL = 0)

        integer MPI_WIN_NULL
        parameter (MPI_WIN_NULL = 0)

        double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
        external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK

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


