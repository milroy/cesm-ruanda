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

! Subprogram not used   logical function seq_mctext_gsmapIdentical(gsmap1,gsmap2)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     type(mct_gsMap), intent(IN):: gsmap1
! Subprogram not used     type(mct_gsMap), intent(IN):: gsmap2
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     character(len=*),parameter :: subname = "(seq_mctext_gsmapIdentical) "
! Subprogram not used     integer :: n
! Subprogram not used     logical :: identical
! Subprogram not used 
! Subprogram not used     !-----------------------
! Subprogram not used 
! Subprogram not used     identical = .true.
! Subprogram not used 
! Subprogram not used     ! --- continue compare ---
! Subprogram not used     if (identical) then
! Subprogram not used        if (mct_gsMap_gsize(gsmap1) /= mct_gsMap_gsize(gsmap2)) identical = .false.
! Subprogram not used        if (mct_gsMap_ngseg(gsmap1) /= mct_gsMap_ngseg(gsmap2)) identical = .false.
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! --- continue compare ---
! Subprogram not used     if (identical) then
! Subprogram not used        do n = 1,mct_gsMap_ngseg(gsmap1)
! Subprogram not used           if (gsmap1%start(n)  /= gsmap2%start(n) ) identical = .false.
! Subprogram not used           if (gsmap1%length(n) /= gsmap2%length(n)) identical = .false.
! Subprogram not used           if (gsmap1%pe_loc(n) /= gsmap2%pe_loc(n)) identical = .false.
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     seq_mctext_gsmapIdentical = identical
! Subprogram not used 
! Subprogram not used   end function seq_mctext_gsmapIdentical

end module cplcomp_exchange_mod
