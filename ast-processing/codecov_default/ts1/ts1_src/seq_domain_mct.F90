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



