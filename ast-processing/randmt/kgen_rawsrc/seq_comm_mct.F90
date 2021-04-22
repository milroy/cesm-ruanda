module seq_comm_mct

!---------------------------------------------------------------------
!
! Purpose: MCT utitlity functions used in sequential CCSM.
!          Note that if no MPI, will call MCTs fake version
!          (including mpif.h) will be utilized
!
!---------------------------------------------------------------------


!!! NOTE: If all atmospheres are identical in number of processes,
!!! number of threads, and grid layout, we should check that the
!!! user-provided number of processes and threads are consistent
!!! (or else, only accept one entry for these quantities when reading
!!! the namelist).  ARE OTHER PROTECTIONS/CHECKS NEEDED???


  use mct_mod     , only : mct_world_init, mct_die
  use shr_sys_mod , only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod , only : shr_mpi_chkerr, shr_mpi_bcast, shr_mpi_max
  use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit

  implicit none

  private

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

  public seq_comm_init
  public seq_comm_iamin
  public seq_comm_iamroot
  public seq_comm_mpicom
  public seq_comm_iam
  public seq_comm_gloiam
  public seq_comm_gloroot
  public seq_comm_cplpe
  public seq_comm_cmppe
  public seq_comm_name
  public seq_comm_inst
  public seq_comm_suffix
  public seq_comm_petlist
  public seq_comm_setptrs
  public seq_comm_setnthreads
  public seq_comm_getnthreads
  public seq_comm_printcomms

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer, public, parameter :: default_logunit = 6
  integer, public :: logunit  = default_logunit     ! log unit number
  integer, public :: loglevel = 1     ! log level

  integer, public :: global_mype = -1  !! To be initialized

  !!! Note - NUM_COMP_INST_XXX are cpp variables set in buildlib.csm_share

  integer, parameter :: nphysmod = 7  ! number of physical models
  integer, parameter, public :: num_inst_atm = 1
  integer, parameter, public :: num_inst_lnd = 1
  integer, parameter, public :: num_inst_ocn = 1
  integer, parameter, public :: num_inst_ice = 1
  integer, parameter, public :: num_inst_glc = 1
  integer, parameter, public :: num_inst_wav = 1
  integer, parameter, public :: num_inst_rof = 1

  integer, parameter, public :: num_inst_total= num_inst_atm + &
                                                num_inst_lnd + &
                                                num_inst_ocn + &
                                                num_inst_ice + &
                                                num_inst_glc + &
                                                num_inst_wav + &
                                                num_inst_rof + 1

  integer, public :: num_inst_min, num_inst_max
  integer, public :: num_inst_xao    ! for xao flux
  integer, public :: num_inst_frc    ! for fractions

  !!! Each component instance needs two communicators: one internal to the
  !!! instance, and one for communicating with the coupler.
  !!! Additionally, one communicator is needed for the coupler's
  !!! internal communications, and one is needed for the global space.

  integer, parameter, public :: num_inst_phys = num_inst_atm + num_inst_lnd + &
                                                num_inst_ocn + num_inst_ice + &
                                                num_inst_glc + num_inst_rof + &
                                                num_inst_wav
  integer, parameter :: ncomps = (2 + 2*nphysmod + (2 * num_inst_phys))

  integer, public :: GLOID
  integer, public :: CPLID

  integer, public :: ALLATMID
  integer, public :: ALLLNDID
  integer, public :: ALLOCNID
  integer, public :: ALLICEID
  integer, public :: ALLGLCID
  integer, public :: ALLROFID
  integer, public :: ALLWAVID

  integer, public :: CPLALLATMID
  integer, public :: CPLALLLNDID
  integer, public :: CPLALLOCNID
  integer, public :: CPLALLICEID
  integer, public :: CPLALLGLCID
  integer, public :: CPLALLROFID
  integer, public :: CPLALLWAVID

  integer, public :: ATMID(num_inst_atm)
  integer, public :: LNDID(num_inst_lnd)
  integer, public :: OCNID(num_inst_ocn)
  integer, public :: ICEID(num_inst_ice)
  integer, public :: GLCID(num_inst_glc)
  integer, public :: ROFID(num_inst_rof)
  integer, public :: WAVID(num_inst_wav)

  integer, public :: CPLATMID(num_inst_atm)
  integer, public :: CPLLNDID(num_inst_lnd)
  integer, public :: CPLOCNID(num_inst_ocn)
  integer, public :: CPLICEID(num_inst_ice)
  integer, public :: CPLGLCID(num_inst_glc)
  integer, public :: CPLROFID(num_inst_rof)
  integer, public :: CPLWAVID(num_inst_wav)

  integer, parameter, public :: seq_comm_namelen=16
  type seq_comm_type
    character(len=seq_comm_namelen) :: name     ! my name
    character(len=seq_comm_namelen) :: suffix   ! recommended suffix
    integer :: inst            ! my inst index
    integer :: ID              ! my id number
    integer :: mpicom          ! mpicom
    integer :: mpigrp          ! mpigrp
    integer :: npes            ! number of mpi tasks in comm
    integer :: nthreads        ! number of omp threads per task
    integer :: iam             ! my task number in mpicom
    logical :: iamroot         ! am i the root task in mpicom
    integer :: gloiam          ! my task number in mpi_comm_world
    integer :: gloroot         ! the global task number of each comps root on all pes
    integer :: pethreads       ! max number of threads on my task
    integer :: cplpe           ! a common task in mpicom from the cpl group for join mpicoms
    integer :: cmppe           ! a common task in mpicom from the component group for join mpicoms
    logical :: set             ! has this datatype been set
    integer, pointer    :: petlist(:)  ! esmf pet list
  end type seq_comm_type

  type(seq_comm_type) :: seq_comms(ncomps)

  character(*), parameter :: layout_concurrent = 'concurrent'
  character(*), parameter :: layout_sequential = 'sequential'

  character(*), parameter :: F11 = "(a,a,'(',i3,' ',a,')',a,   3i6,' (',a,i6,')',' (',a,i3,')')"
  character(*), parameter :: F12 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')','(',a,2i6,')')"
  character(*), parameter :: F13 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')')"
  character(*), parameter :: F14 = "(a,a,'(',i3,' ',a,')',a,    6x,' (',a,i6,')',' (',a,i3,')')"
  integer :: Global_Comm


  character(len=32), public :: &
       atm_layout, lnd_layout, ice_layout, glc_layout, rof_layout, &
       ocn_layout, wav_layout
  


!=======================================================================
contains
!=======================================================================
  
  subroutine seq_comm_init(Comm_in, nmlfile)
      
    !----------------------------------------------------------
    !
    ! Arguments
    implicit none
    integer, intent(in) :: Comm_in
    character(len=*), intent(IN) :: nmlfile
    !
    ! Local variables
    !
    logical :: error_state
    integer :: ierr, n, count
    character(*), parameter :: subName =   '(seq_comm_init) '
    integer :: mpi_group_world   ! MPI_COMM_WORLD group
    integer :: mype,numpes,myncomps,max_threads,gloroot
    integer :: atm_inst_tasks, lnd_inst_tasks, ocn_inst_tasks, ice_inst_tasks, &
               glc_inst_tasks, rof_inst_tasks, wav_inst_tasks
    integer :: current_task_rootpe, droot
    integer :: amin(num_inst_atm), amax(num_inst_atm), astr(num_inst_atm)
    integer :: lmin(num_inst_lnd), lmax(num_inst_lnd), lstr(num_inst_lnd)
    integer :: imin(num_inst_ice), imax(num_inst_ice), istr(num_inst_ice)
    integer :: omin(num_inst_ocn), omax(num_inst_ocn), ostr(num_inst_ocn)
    integer :: gmin(num_inst_glc), gmax(num_inst_glc), gstr(num_inst_glc)
    integer :: wmin(num_inst_wav), wmax(num_inst_wav), wstr(num_inst_wav)
    integer :: rmin(num_inst_rof), rmax(num_inst_rof), rstr(num_inst_rof)
    integer :: cmin,cmax,cstr
    integer :: pelist(3,1)       ! start, stop, stride for group
    integer, pointer :: comps(:) ! array with component ids
    integer, pointer :: comms(:) ! array with mpicoms
    integer :: nu, i
    logical, save :: first_pass = .true.   ! 

    integer :: &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads
    namelist /ccsm_pes/  &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, atm_layout, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, lnd_layout, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, ice_layout, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, glc_layout, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, wav_layout, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, rof_layout, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, ocn_layout, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads 
    !----------------------------------------------------------

    ! make sure this is first pass and set comms unset
    if (.not. first_pass) then
       write(logunit,*) trim(subname),' ERROR seq_comm_init already called '
       call shr_sys_abort()
    endif
    first_pass = .false.
    Global_Comm = Comm_in

    !! Initialize seq_comms elements

    do n = 1,ncomps
       seq_comms(n)%name = 'unknown'
       seq_comms(n)%suffix = ' '
       seq_comms(n)%inst = 0
       seq_comms(n)%set = .false.
       seq_comms(n)%mpicom = MPI_COMM_NULL    ! do some initialization here 
       seq_comms(n)%iam = -1
       seq_comms(n)%iamroot = .false.
       seq_comms(n)%npes = -1
       seq_comms(n)%nthreads = -1
       seq_comms(n)%gloiam = -1
       seq_comms(n)%gloroot = -1
       seq_comms(n)%pethreads = -1
       seq_comms(n)%cplpe = -1
       seq_comms(n)%cmppe = -1
    enddo


    ! Initialize MPI
    ! Note that if no MPI, will call MCTs fake version

    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')
    call mpi_comm_size(GLOBAL_COMM, numpes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    ! Initialize gloiam on all IDs

    global_mype = mype
 
    do n = 1,ncomps
       seq_comms(n)%gloiam = mype
    enddo

    ! Set ntasks, rootpe, pestride, nthreads for all components

    if (mype == 0) then

       !! Set up default atmosphere process parameters

       atm_ntasks = numpes
       atm_rootpe = 0
       atm_pestride = 1
       atm_nthreads = 1
       atm_layout = trim(layout_concurrent)

       lnd_ntasks = numpes
       lnd_rootpe = 0
       lnd_pestride = 1
       lnd_nthreads = 1
       lnd_layout = trim(layout_concurrent)

       ocn_ntasks = numpes
       ocn_rootpe = 0
       ocn_pestride = 1
       ocn_nthreads = 1
       ocn_layout = trim(layout_concurrent)

       ice_ntasks = numpes
       ice_rootpe = 0
       ice_pestride = 1
       ice_nthreads = 1
       ice_layout = trim(layout_concurrent)

       glc_ntasks = numpes
       glc_rootpe = 0
       glc_pestride = 1
       glc_nthreads = 1
       glc_layout = trim(layout_concurrent)

       rof_ntasks = numpes
       rof_rootpe = 0
       rof_pestride = 1
       rof_nthreads = 1
       rof_layout = trim(layout_concurrent)

       wav_ntasks = numpes
       wav_rootpe = 0
       wav_pestride = 1
       wav_nthreads = 1
       wav_layout = trim(layout_concurrent)

       cpl_ntasks = numpes
       cpl_rootpe = 0
       cpl_pestride = 1
       cpl_nthreads = 1

       ! Read namelist if it exists

       nu = shr_file_getUnit()
       open(nu, file=trim(nmlfile), status='old', iostat=ierr)

       if (ierr == 0) then
          ierr = 1
          do while( ierr > 0 )
             read(nu, nml=ccsm_pes, iostat=ierr)
          end do
          close(nu)
       end if
       call shr_file_freeUnit(nu)

    end if

    !--- compute some other num_inst values
       
    num_inst_xao = max(num_inst_atm,num_inst_ocn)
    num_inst_frc = num_inst_ice

    !--- compute num_inst_min, num_inst_max
    !--- instances must be either 1 or a constant across components
    !--- checks for prognostic/present consistency in the driver

    error_state = .false.
    num_inst_min = num_inst_atm
    num_inst_min = min(num_inst_min, num_inst_lnd)
    num_inst_min = min(num_inst_min, num_inst_ocn)
    num_inst_min = min(num_inst_min, num_inst_ice)
    num_inst_min = min(num_inst_min, num_inst_glc)
    num_inst_min = min(num_inst_min, num_inst_wav)
    num_inst_min = min(num_inst_min, num_inst_rof)
    num_inst_max = num_inst_atm
    num_inst_max = max(num_inst_max, num_inst_lnd)
    num_inst_max = max(num_inst_max, num_inst_ocn)
    num_inst_max = max(num_inst_max, num_inst_ice)
    num_inst_max = max(num_inst_max, num_inst_glc)
    num_inst_max = max(num_inst_max, num_inst_wav)
    num_inst_max = max(num_inst_max, num_inst_rof)

    if (num_inst_min /= num_inst_max .and. num_inst_min /= 1) error_state = .true.
    if (num_inst_atm /= num_inst_min .and. num_inst_atm /= num_inst_max) error_state = .true.
    if (num_inst_lnd /= num_inst_min .and. num_inst_lnd /= num_inst_max) error_state = .true.
    if (num_inst_ocn /= num_inst_min .and. num_inst_ocn /= num_inst_max) error_state = .true.
    if (num_inst_ice /= num_inst_min .and. num_inst_ice /= num_inst_max) error_state = .true.
    if (num_inst_glc /= num_inst_min .and. num_inst_glc /= num_inst_max) error_state = .true.
    if (num_inst_wav /= num_inst_min .and. num_inst_wav /= num_inst_max) error_state = .true.
    if (num_inst_rof /= num_inst_min .and. num_inst_rof /= num_inst_max) error_state = .true.

    if (error_state) then
       write(logunit,*) trim(subname),' ERROR: num_inst inconsistent'
       call shr_sys_abort(trim(subname)//' ERROR: num_inst inconsistent')
    endif

    ! Initialize IDs

    count = 0

    count = count + 1
    GLOID = count
    count = count + 1
    CPLID = count

    count = count + 1
    ALLATMID = count
    count = count + 1
    ALLLNDID = count
    count = count + 1
    ALLOCNID = count
    count = count + 1
    ALLICEID = count
    count = count + 1
    ALLGLCID = count
    count = count + 1
    ALLROFID = count
    count = count + 1
    ALLWAVID = count

    count = count + 1
    CPLALLATMID = count
    count = count + 1
    CPLALLLNDID = count
    count = count + 1
    CPLALLOCNID = count
    count = count + 1
    CPLALLICEID = count
    count = count + 1
    CPLALLGLCID = count
    count = count + 1
    CPLALLROFID = count
    count = count + 1
    CPLALLWAVID = count

    do n = 1, num_inst_atm
       count = count + 1
       ATMID(n) = count
       count = count + 1
       CPLATMID(n) = count
    end do
         
    do n = 1, num_inst_lnd
       count = count + 1
       LNDID(n) = count
       count = count + 1
       CPLLNDID(n) = count
    end do
       
    do n = 1, num_inst_ocn
       count = count + 1
       OCNID(n) = count
       count = count + 1
       CPLOCNID(n) = count
    end do
       
    do n = 1, num_inst_ice
       count = count + 1
       ICEID(n) = count
       count = count + 1
       CPLICEID(n) = count
    end do

    do n = 1, num_inst_glc
       count = count + 1
       GLCID(n) = count
       count = count + 1
       CPLGLCID(n) = count
    end do

    do n = 1, num_inst_rof
       count = count + 1
       ROFID(n) = count
       count = count + 1
       CPLROFID(n) = count
    end do

    do n = 1, num_inst_wav
       count = count + 1
       WAVID(n) = count
       count = count + 1
       CPLWAVID(n) = count
    end do

    if (count /= ncomps) then
       write(logunit,*) trim(subname),' ERROR in ID count ',count,ncomps
       call shr_sys_abort(trim(subname)//' ERROR in ID count')
    endif

    if (mype == 0) then
       !--- validation of inputs ---
       ! rootpes >= 0

       error_state = .false.

       if (atm_rootpe < 0) error_state = .true.
       if (lnd_rootpe < 0) error_state = .true.
       if (ice_rootpe < 0) error_state = .true.
       if (ocn_rootpe < 0) error_state = .true.
       if (glc_rootpe < 0) error_state = .true.
       if (wav_rootpe < 0) error_state = .true.
       if (rof_rootpe < 0) error_state = .true.
       if (cpl_rootpe < 0) error_state = .true.
       
       if (error_state) then
          write(logunit,*) trim(subname),' ERROR: rootpes must be >= 0'
          call shr_sys_abort(trim(subname)//' ERROR: rootpes >= 0')
       endif

!       ! nthreads = 1, temporary
!       if (atm_nthreads /= 1 .or. lnd_nthreads /= 1 .or. ice_nthreads /= 1 .or. &
!           ocn_nthreads /= 1 .or. cpl_nthreads /= 1) then
!          write(logunit,*) trim(subname),' ERROR: nthreads must be 1'
!          call shr_sys_abort()
!       endif

!       ! nthreads should be 1 or something consistent, compute max nthreads
!       amax = max(atm_nthreads,lnd_nthreads)
!       amax = max(amax        ,ice_nthreads)
!       amax = max(amax        ,ocn_nthreads)
!       amax = max(amax        ,cpl_nthreads)

!       ! check that everything is either 1 or max nthreads
!       if ((atm_nthreads /= 1 .and. atm_nthreads /= amax) .or. &
!           (lnd_nthreads /= 1 .and. lnd_nthreads /= amax) .or. &
!           (ice_nthreads /= 1 .and. ice_nthreads /= amax) .or. &
!           (ocn_nthreads /= 1 .and. ocn_nthreads /= amax) .or. &
!           (cpl_nthreads /= 1 .and. cpl_nthreads /= amax)) then
!          write(logunit,*) trim(subname),' ERROR: nthreads must be consistent'
!          call shr_sys_abort()
!       endif

       !! Determine the process layout
       !!
       !! We will assign atm_ntasks / num_inst_atm tasks to each atmosphere
       !! instance.  (This may lead to unallocated tasks if atm_ntasks is
       !! not an integer multiple of num_inst_atm.)

       if (trim(atm_layout) == trim(layout_concurrent)) then
          atm_inst_tasks = atm_ntasks / num_inst_atm
          droot = (atm_inst_tasks * atm_pestride)
       elseif (trim(atm_layout) == trim(layout_sequential)) then
          atm_inst_tasks = atm_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid atm_layout ')
       endif
       current_task_rootpe = atm_rootpe
       do n = 1, num_inst_atm
          amin(n) = current_task_rootpe
          amax(n) = current_task_rootpe &
                    + ((atm_inst_tasks - 1) * atm_pestride)
          astr(n) = atm_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Land instance tasks

       if (trim(lnd_layout) == trim(layout_concurrent)) then
          lnd_inst_tasks = lnd_ntasks / num_inst_lnd
          droot = (lnd_inst_tasks * lnd_pestride)
       elseif (trim(lnd_layout) == trim(layout_sequential)) then
          lnd_inst_tasks = lnd_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid lnd_layout ')
       endif
       current_task_rootpe = lnd_rootpe
       do n = 1, num_inst_lnd
          lmin(n) = current_task_rootpe
          lmax(n) = current_task_rootpe &
                    + ((lnd_inst_tasks - 1) * lnd_pestride)
          lstr(n) = lnd_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Ocean instance tasks
       
       if (trim(ocn_layout) == trim(layout_concurrent)) then
          ocn_inst_tasks = ocn_ntasks / num_inst_ocn
          droot = (ocn_inst_tasks * ocn_pestride)
       elseif (trim(ocn_layout) == trim(layout_sequential)) then
          ocn_inst_tasks = ocn_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid ocn_layout ')
       endif
       current_task_rootpe = ocn_rootpe
       do n = 1, num_inst_ocn
          omin(n) = current_task_rootpe
          omax(n) = current_task_rootpe &
                    + ((ocn_inst_tasks - 1) * ocn_pestride)
          ostr(n) = ocn_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Sea ice instance tasks
       
       if (trim(ice_layout) == trim(layout_concurrent)) then
          ice_inst_tasks = ice_ntasks / num_inst_ice
          droot = (ice_inst_tasks * ice_pestride)
       elseif (trim(ice_layout) == trim(layout_sequential)) then
          ice_inst_tasks = ice_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid ice_layout ')
       endif
       current_task_rootpe = ice_rootpe
       do n = 1, num_inst_ice
          imin(n) = current_task_rootpe
          imax(n) = current_task_rootpe &
                    + ((ice_inst_tasks - 1) * ice_pestride)
          istr(n) = ice_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Glacier instance tasks
       
       if (trim(glc_layout) == trim(layout_concurrent)) then
          glc_inst_tasks = glc_ntasks / num_inst_glc
          droot = (glc_inst_tasks * glc_pestride)
       elseif (trim(glc_layout) == trim(layout_sequential)) then
          glc_inst_tasks = glc_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid glc_layout ')
       endif
       current_task_rootpe = glc_rootpe
       do n = 1, num_inst_glc
          gmin(n) = current_task_rootpe
          gmax(n) = current_task_rootpe &
                    + ((glc_inst_tasks - 1) * glc_pestride)
          gstr(n) = glc_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Runoff instance tasks
       
       if (trim(rof_layout) == trim(layout_concurrent)) then
          rof_inst_tasks = rof_ntasks / num_inst_rof
          droot = (rof_inst_tasks * rof_pestride)
       elseif (trim(rof_layout) == trim(layout_sequential)) then
          rof_inst_tasks = rof_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid rof_layout ')
       endif
       current_task_rootpe = rof_rootpe
       do n = 1, num_inst_rof
          rmin(n) = current_task_rootpe
          rmax(n) = current_task_rootpe &
                    + ((rof_inst_tasks - 1) * rof_pestride)
          rstr(n) = rof_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Wave instance tasks
       
       if (trim(wav_layout) == trim(layout_concurrent)) then
          wav_inst_tasks = wav_ntasks / num_inst_wav
          droot = (wav_inst_tasks * wav_pestride)
       elseif (trim(wav_layout) == trim(layout_sequential)) then
          wav_inst_tasks = wav_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid wav_layout ')
       endif
       current_task_rootpe = wav_rootpe
       do n = 1, num_inst_wav
          wmin(n) = current_task_rootpe
          wmax(n) = current_task_rootpe &
                    + ((wav_inst_tasks - 1) * wav_pestride)
          wstr(n) = wav_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Coupler tasks

       cmin = cpl_rootpe
       cmax = cpl_rootpe + (cpl_ntasks-1)*cpl_pestride
       cstr = cpl_pestride
    end if

       
    call shr_mpi_bcast(atm_nthreads,GLOBAL_COMM,'atm_nthreads')
    call shr_mpi_bcast(lnd_nthreads,GLOBAL_COMM,'lnd_nthreads')
    call shr_mpi_bcast(ocn_nthreads,GLOBAL_COMM,'ocn_nthreads')
    call shr_mpi_bcast(ice_nthreads,GLOBAL_COMM,'ice_nthreads')
    call shr_mpi_bcast(glc_nthreads,GLOBAL_COMM,'glc_nthreads')
    call shr_mpi_bcast(wav_nthreads,GLOBAL_COMM,'wav_nthreads')
    call shr_mpi_bcast(rof_nthreads,GLOBAL_COMM,'rof_nthreads')
    call shr_mpi_bcast(cpl_nthreads,GLOBAL_COMM,'cpl_nthreads')

    ! Create MPI communicator groups

    if (mype == 0) then
       pelist(1,1) = 0
       pelist(2,1) = numpes-1
       pelist(3,1) = 1
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(GLOID, pelist,iname='GLOBAL')

    if (mype == 0) then
       pelist(1,1) = cmin
       pelist(2,1) = cmax
       pelist(3,1) = cstr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(CPLID,pelist,cpl_nthreads,'CPL')

    do n = 1, num_inst_atm
       if (mype == 0) then
          pelist(1,1) = amin(n)
          pelist(2,1) = amax(n)
          pelist(3,1) = astr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ATMID(n), pelist, atm_nthreads, 'ATM', n, num_inst_atm)
       call seq_comm_joincomm(CPLID, ATMID(n), CPLATMID(n), 'CPLATM', n, num_inst_atm)
    end do
    call seq_comm_jcommarr(ATMID,ALLATMID,'ALLATMID',1,1)
    call seq_comm_joincomm(CPLID,ALLATMID,CPLALLATMID,'CPLALLATMID',1,1)

    do n = 1, num_inst_lnd
       if (mype == 0) then
          pelist(1,1) = lmin(n)
          pelist(2,1) = lmax(n)
          pelist(3,1) = lstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(LNDID(n), pelist, lnd_nthreads, 'LND', n, num_inst_lnd)
       call seq_comm_joincomm(CPLID, LNDID(n), CPLLNDID(n), 'CPLLND', n, num_inst_lnd)
    end do
    call seq_comm_jcommarr(LNDID,ALLLNDID,'ALLLNDID',1,1)
    call seq_comm_joincomm(CPLID,ALLLNDID,CPLALLLNDID,'CPLALLLNDID',1,1)

    do n = 1, num_inst_ocn
       if (mype == 0) then
          pelist(1,1) = omin(n)
          pelist(2,1) = omax(n)
          pelist(3,1) = ostr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(OCNID(n), pelist, ocn_nthreads, 'OCN', n, num_inst_ocn)
       call seq_comm_joincomm(CPLID, OCNID(n), CPLOCNID(n), 'CPLOCN', n, num_inst_ocn)
    end do
    call seq_comm_jcommarr(OCNID,ALLOCNID,'ALLOCNID',1,1)
    call seq_comm_joincomm(CPLID,ALLOCNID,CPLALLOCNID,'CPLALLOCNID',1,1)

    do n = 1, num_inst_ice
       if (mype == 0) then
          pelist(1,1) = imin(n)
          pelist(2,1) = imax(n)
          pelist(3,1) = istr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ICEID(n), pelist, ice_nthreads, 'ICE', n, num_inst_ice)
       call seq_comm_joincomm(CPLID, ICEID(n), CPLICEID(n), 'CPLICE', n, num_inst_ice)
    end do
    call seq_comm_jcommarr(ICEID,ALLICEID,'ALLICEID',1,1)
    call seq_comm_joincomm(CPLID,ALLICEID,CPLALLICEID,'CPLALLICEID',1,1)

    do n = 1, num_inst_glc
       if (mype == 0) then
          pelist(1,1) = gmin(n)
          pelist(2,1) = gmax(n)
          pelist(3,1) = gstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(GLCID(n), pelist, glc_nthreads, 'GLC', n, num_inst_glc)
       call seq_comm_joincomm(CPLID, GLCID(n), CPLGLCID(n), 'CPLGLC', n, num_inst_glc)
    end do
    call seq_comm_jcommarr(GLCID,ALLGLCID,'ALLGLCID',1,1)
    call seq_comm_joincomm(CPLID,ALLGLCID,CPLALLGLCID,'CPLALLGLCID',1,1)

    do n = 1, num_inst_rof
       if (mype == 0) then
          pelist(1,1) = rmin(n)
          pelist(2,1) = rmax(n)
          pelist(3,1) = rstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ROFID(n), pelist, rof_nthreads, 'ROF', n, num_inst_rof)
       call seq_comm_joincomm(CPLID, ROFID(n), CPLROFID(n), 'CPLROF', n, num_inst_rof)
    end do
    call seq_comm_jcommarr(ROFID,ALLROFID,'ALLROFID',1,1)
    call seq_comm_joincomm(CPLID,ALLROFID,CPLALLROFID,'CPLALLROFID',1,1)

    do n = 1, num_inst_wav
       if (mype == 0) then
          pelist(1,1) = wmin(n)
          pelist(2,1) = wmax(n)
          pelist(3,1) = wstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(WAVID(n), pelist, wav_nthreads, 'WAV', n, num_inst_wav)
       call seq_comm_joincomm(CPLID, WAVID(n), CPLWAVID(n), 'CPLWAV', n, num_inst_wav)
    end do
    call seq_comm_jcommarr(WAVID,ALLWAVID,'ALLWAVID',1,1)
    call seq_comm_joincomm(CPLID,ALLWAVID,CPLALLWAVID,'CPLALLWAVID',1,1)

    !! Count the total number of threads

    max_threads = -1
    do n = 1,ncomps
       max_threads = max(max_threads,seq_comms(n)%nthreads)
    enddo
    do n = 1,ncomps
       seq_comms(n)%pethreads = max_threads
    enddo

    ! compute each components root pe global id and broadcast so all pes have info

    do n = 1,ncomps
       gloroot = -999
       if (seq_comms(n)%iamroot) gloroot = seq_comms(n)%gloiam
       call shr_mpi_max(gloroot,seq_comms(n)%gloroot,GLOBAL_COMM, &
                        trim(subname)//' gloroot',all=.true.)
    enddo

    ! Initialize MCT

    ! add up valid comps on local pe

    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
       endif
    enddo

    ! set comps and comms

    allocate(comps(myncomps),comms(myncomps),stat=ierr)
    if(ierr/=0) call mct_die(subName,'allocate comps comms',ierr)

    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
          if (myncomps > size(comps)) then
             write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
             call shr_sys_abort()
          endif
          comps(myncomps) = seq_comms(n)%ID
          comms(myncomps) = seq_comms(n)%mpicom
       endif
    enddo

    if (myncomps /= size(comps)) then
       write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
       call shr_sys_abort()
    endif

    call mct_world_init(ncomps, GLOBAL_COMM, comms, comps)

    deallocate(comps,comms)

    call seq_comm_printcomms()

  end subroutine seq_comm_init

!---------------------------------------------------------
  subroutine seq_comm_setcomm(ID,pelist,nthreads,iname,inst,tinst)

    implicit none
    integer,intent(IN) :: ID
    integer,intent(IN) :: pelist(:,:)
    integer,intent(IN),optional :: nthreads
    character(len=*),intent(IN),optional :: iname  ! name of component
    integer,intent(IN),optional :: inst  ! instance of component
    integer,intent(IN),optional :: tinst ! total number of instances for this component

    integer :: mpigrp_world
    integer :: mpigrp
    integer :: mpicom
    integer :: ntask,ntasks,cnt
    integer :: ierr
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    character(*),parameter :: subName =   '(seq_comm_setcomm) '

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif 

    call mpi_comm_group(GLOBAL_COMM, mpigrp_world, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
    call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    ntasks = ((pelist(2,1) - pelist(1,1)) / pelist(3,1)) + 1
    allocate(seq_comms(ID)%petlist(ntasks))
    cnt = 0
    do ntask = pelist(1,1),pelist(2,1),pelist(3,1)
        cnt = cnt + 1
        if (cnt > ntasks) then
           write(logunit,*) subname,' ERROR in petlist init ',ntasks,pelist(1:3,1),ntask,cnt
           call shr_sys_abort(subname//' ERROR in petlist init')
        endif
        seq_comms(ID)%petlist(cnt) = ntask
    enddo

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID

    if (present(inst)) then
       seq_comms(ID)%inst = inst
       set_suffix = .true.
    else
       seq_comms(ID)%inst = 1
       set_suffix = .false.
    endif

    if (present(tinst)) then
       if (tinst == 1) set_suffix = .false.
    endif

    if (present(iname)) then
       seq_comms(ID)%name = trim(iname)
       if (set_suffix) then
          call seq_comm_mkname(cname,iname,seq_comms(ID)%inst)
          seq_comms(ID)%name = trim(cname)
       endif
    endif

    if (set_suffix) then
       call seq_comm_mkname(cname,'_',seq_comms(ID)%inst)
       seq_comms(ID)%suffix = trim(cname)
    else
       seq_comms(ID)%suffix = ' '
    endif

    seq_comms(ID)%mpicom = mpicom
    seq_comms(ID)%mpigrp = mpigrp
    if (present(nthreads)) then
       seq_comms(ID)%nthreads = nthreads
    else
       seq_comms(ID)%nthreads = 1
    endif

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom,seq_comms(ID)%npes,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size')
       call mpi_comm_rank(mpicom,seq_comms(ID)%iam,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank')
       if (seq_comms(ID)%iam == 0) then
          seq_comms(ID)%iamroot = .true.
       else
          seq_comms(ID)%iamroot = .false.
       endif
    else
       seq_comms(ID)%npes = -1
       seq_comms(ID)%iam = -1
       seq_comms(ID)%nthreads = 1
       seq_comms(ID)%iamroot = .false.
    endif

    if (seq_comms(ID)%iamroot) then
       write(logunit,F11) trim(subname),'  initialize ID ',ID,seq_comms(ID)%name, &
         ' pelist   =',pelist,' npes =',seq_comms(ID)%npes,' nthreads =',seq_comms(ID)%nthreads
    endif

  end subroutine seq_comm_setcomm

!---------------------------------------------------------
  subroutine seq_comm_joincomm(ID1,ID2,ID,iname,inst,tinst)

    implicit none
    integer,intent(IN) :: ID1    ! src id
    integer,intent(IN) :: ID2    ! srd id
    integer,intent(IN) :: ID     ! computed id
    character(len=*),intent(IN),optional :: iname  ! comm name
    integer,intent(IN),optional :: inst
    integer,intent(IN),optional :: tinst

    integer :: mpigrp
    integer :: mpicom
    integer :: ierr
    integer :: n,nsize
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    integer,allocatable :: pe_t1(:),pe_t2(:)
    character(*),parameter :: subName =   '(seq_comm_joincomm) '

    ! check that IDs are in valid range, that ID1 and ID2 have
    ! been set, and that ID has not been set

    if (ID1 < 1 .or. ID1 > ncomps) then
       write(logunit,*) subname,' ID1 out of range, abort ',ID1
       call shr_sys_abort()
    endif 
    if (ID2 < 1 .or. ID2 > ncomps) then
       write(logunit,*) subname,' ID2 out of range, abort ',ID2
       call shr_sys_abort()
    endif 
    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif
    if (.not. seq_comms(ID1)%set .or. .not. seq_comms(ID2)%set) then
       write(logunit,*) subname,' ID1 or ID2 not set ',ID1,ID2
       call shr_sys_abort()
    endif
    if (seq_comms(ID)%set) then
       write(logunit,*) subname,' ID already set ',ID
       call shr_sys_abort()
    endif

    call mpi_group_union(seq_comms(ID1)%mpigrp,seq_comms(ID2)%mpigrp,mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_union mpigrp')
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID

    if (present(inst)) then
       seq_comms(ID)%inst = inst
    else
       seq_comms(ID)%inst = 1
    endif

    set_suffix = .true.
    if (present(tinst)) then
       if (tinst == 1) set_suffix = .false.
    endif

    if (present(iname)) then
       seq_comms(ID)%name = trim(iname)
       if (set_suffix) then
          call seq_comm_mkname(cname,iname,seq_comms(ID)%inst)
          seq_comms(ID)%name = trim(cname)
       endif
    endif

    if (set_suffix) then
       call seq_comm_mkname(cname,'_',seq_comms(ID)%inst)
       seq_comms(ID)%suffix = trim(cname)
    else
       seq_comms(ID)%suffix = ' '
    endif

    seq_comms(ID)%mpicom = mpicom
    seq_comms(ID)%mpigrp = mpigrp
    seq_comms(ID)%nthreads = max(seq_comms(ID1)%nthreads,seq_comms(ID2)%nthreads)
    seq_comms(ID)%nthreads = max(seq_comms(ID)%nthreads,1)

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom,seq_comms(ID)%npes,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size')
       call mpi_comm_rank(mpicom,seq_comms(ID)%iam,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank')
       if (seq_comms(ID)%iam == 0) then
          seq_comms(ID)%iamroot = .true.
       else
          seq_comms(ID)%iamroot = .false.
       endif
    else
       seq_comms(ID)%npes = -1
       seq_comms(ID)%iam = -1
       seq_comms(ID)%iamroot = .false.
    endif

! needs to be excluded until mpi_group_size is added to serial mpi in mct

    allocate(pe_t1(1),pe_t2(1))
    pe_t1(1) = 0
    call mpi_group_translate_ranks(seq_comms(ID1)%mpigrp, 1, pe_t1, mpigrp, pe_t2, ierr)
    seq_comms(ID)%cplpe = pe_t2(1)
    pe_t1(1) = 0
    call mpi_group_translate_ranks(seq_comms(ID2)%mpigrp, 1, pe_t1, mpigrp, pe_t2, ierr)
    seq_comms(ID)%cmppe = pe_t2(1)
    deallocate(pe_t1,pe_t2)

    if (seq_comms(ID)%iamroot) then
       if (loglevel > 1) then
          write(logunit,F12) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads, &
          ' cpl/cmp pes =',seq_comms(ID)%cplpe,seq_comms(ID)%cmppe
       else
          write(logunit,F13) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       endif
    endif

  end subroutine seq_comm_joincomm

!---------------------------------------------------------
  subroutine seq_comm_jcommarr(IDs,ID,iname,inst,tinst)

    implicit none
    integer,intent(IN) :: IDs(:) ! src id
    integer,intent(IN) :: ID     ! computed id
    character(len=*),intent(IN),optional :: iname  ! comm name
    integer,intent(IN),optional :: inst
    integer,intent(IN),optional :: tinst

    integer :: mpigrp, mpigrpp
    integer :: mpicom, nids
    integer :: ierr
    integer :: n,nsize
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    integer,allocatable :: pe_t1(:),pe_t2(:)
    character(*),parameter :: subName =   '(seq_comm_jcommarr) '

    ! check that IDs are in valid range, that IDs have
    ! been set, and that ID has not been set

    nids = size(IDs)
    do n = 1,nids
       if (IDs(n) < 1 .or. IDs(n) > ncomps) then
          write(logunit,*) subname,' IDs out of range, abort ',n,IDs(n)
          call shr_sys_abort()
       endif 
       if (.not. seq_comms(IDs(n))%set) then
          write(logunit,*) subname,' IDs not set ',n,IDs(n)
          call shr_sys_abort()
       endif
    enddo

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif
    if (seq_comms(ID)%set) then
       write(logunit,*) subname,' ID already set ',ID
       call shr_sys_abort()
    endif

    mpigrp = seq_comms(IDs(1))%mpigrp
    do n = 1,nids
       mpigrpp = mpigrp
       call mpi_group_union(mpigrpp,seq_comms(IDs(n))%mpigrp,mpigrp,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_union mpigrp')
    enddo
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID

    if (present(inst)) then
       seq_comms(ID)%inst = inst
    else
       seq_comms(ID)%inst = 1
    endif

    set_suffix = .true.
    if (present(tinst)) then
       if (tinst == 1) set_suffix = .false.
    endif

    if (present(iname)) then
       seq_comms(ID)%name = trim(iname)
       if (set_suffix) then
          call seq_comm_mkname(cname,iname,seq_comms(ID)%inst)
          seq_comms(ID)%name = trim(cname)
       endif
    endif

    if (set_suffix) then
       call seq_comm_mkname(cname,'_',seq_comms(ID)%inst)
       seq_comms(ID)%suffix = trim(cname)
    else
       seq_comms(ID)%suffix = ' '
    endif

    seq_comms(ID)%mpicom = mpicom
    seq_comms(ID)%mpigrp = mpigrp
    
    seq_comms(ID)%nthreads = 1
    do n = 1,nids
       seq_comms(ID)%nthreads = max(seq_comms(ID)%nthreads,seq_comms(IDs(n))%nthreads)
    enddo

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom,seq_comms(ID)%npes,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size')
       call mpi_comm_rank(mpicom,seq_comms(ID)%iam,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank')
       if (seq_comms(ID)%iam == 0) then
          seq_comms(ID)%iamroot = .true.
       else
          seq_comms(ID)%iamroot = .false.
       endif
    else
       seq_comms(ID)%npes = -1
       seq_comms(ID)%iam = -1
       seq_comms(ID)%iamroot = .false.
    endif

    seq_comms(ID)%cplpe = -1
    seq_comms(ID)%cmppe = -1

    if (seq_comms(ID)%iamroot) then
       if (loglevel > 1) then
          write(logunit,F14) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join multiple comp IDs',' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       else
          write(logunit,F14) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join multiple comp IDs',' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       endif
    endif

  end subroutine seq_comm_jcommarr

!---------------------------------------------------------
  subroutine seq_comm_printcomms()

    implicit none
    character(*),parameter :: subName =   '(seq_comm_printcomms) '
    integer :: m,n,mype,npes,ierr
    character(len=256) :: iamstring
    character(*),parameter :: F01 = "(4x,a4,4x   ,40(1x,a8))"
    character(*),parameter :: F02 = "(4x,i4,3x,a1,40(2x,i6,1x))"
    character(*),parameter :: F03 = "(4x,i4,3x,a1,a)"

    call mpi_comm_size(GLOBAL_COMM, npes  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')
    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    call shr_sys_flush(logunit)
    call mpi_barrier(GLOBAL_COMM,ierr)
    if (mype == 0) then
       do n = 1,ncomps
          write(logunit,'(a,4i6,2x,3a)') trim(subName),n, &
             seq_comms(n)%gloroot,seq_comms(n)%npes,seq_comms(n)%nthreads, &
             trim(seq_comms(n)%name),':',trim(seq_comms(n)%suffix)
       enddo
!       write(logunit,*) ' '
!       write(logunit,*) trim(subName),' ID layout : global pes vs local pe for each ID'
!       write(logunit,F01) ' gpe',(seq_comms(n)%name,n=1,ncomps),'nthrds'
!       write(logunit,F01) ' ---',(' ------ '       ,n=1,ncomps),'------'
       call shr_sys_flush(logunit)
    endif
!    iamstring = ' '
!   do n = 1,ncomps
!      if (seq_comms(n)%iam >= 0) then
!         write(iamstring((n-1)*9+1:n*9),"(2x,i6,1x)") seq_comms(n)%iam
!      endif
!   enddo
!   n = ncomps + 1
!   write(iamstring((n-1)*9+1:n*9),"(2x,i6,1x)") seq_comms(GLOID)%pethreads

!    call shr_sys_flush(logunit)
!    call mpi_barrier(GLOBAL_COMM,ierr)
!   do m = 0,npes-1
!      if (mype == m) then
!!          write(logunit,F02) mype,':',(seq_comms(n)%iam,n=1,ncomps)
!         write(logunit,F03) mype,':',trim(iamstring)
!         if (m == npes-1) then
!            write(logunit,*) ' '
!         endif
!      endif
!      call shr_sys_flush(logunit)
!      call mpi_barrier(GLOBAL_COMM,ierr)
!   enddo

  end subroutine seq_comm_printcomms

!---------------------------------------------------------
  subroutine seq_comm_setptrs(ID,mpicom,mpigrp,npes,nthreads,iam,iamroot,gloiam,gloroot, &
                                 cplpe,cmppe,pethreads, name)

    implicit none
    integer,intent(in) :: ID
    integer,intent(out),optional :: mpicom
    integer,intent(out),optional :: mpigrp
    integer,intent(out),optional :: npes
    integer,intent(out),optional :: nthreads
    integer,intent(out),optional :: iam
    logical,intent(out),optional :: iamroot
    integer,intent(out),optional :: gloiam
    integer,intent(out),optional :: gloroot
    integer,intent(out),optional :: cplpe
    integer,intent(out),optional :: cmppe
    integer,intent(out),optional :: pethreads
    character(len=seq_comm_namelen)  , intent(out), optional :: name
    character(*),parameter :: subName =   '(seq_comm_setptrs) '

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, return ',ID
       return
    endif 

    if (present(mpicom)) then
       mpicom = seq_comms(ID)%mpicom
    endif

    if (present(mpigrp)) then
       mpigrp = seq_comms(ID)%mpigrp
    endif

    if (present(npes)) then
       npes = seq_comms(ID)%npes
    endif

    if (present(nthreads)) then
       nthreads = seq_comms(ID)%nthreads
    endif

    if (present(iam)) then
       iam = seq_comms(ID)%iam
    endif

    if (present(iamroot)) then
       iamroot = seq_comms(ID)%iamroot
    endif

    if (present(gloiam)) then
       gloiam = seq_comms(ID)%gloiam
    endif

    if (present(gloroot)) then
       gloroot = seq_comms(ID)%gloroot
    endif

    if (present(cplpe)) then
       cplpe = seq_comms(ID)%cplpe
    endif

    if (present(cmppe)) then
       cmppe = seq_comms(ID)%cmppe
    endif

    if (present(pethreads)) then
       pethreads = seq_comms(ID)%pethreads
    endif

    if(present(name)) then
       name = seq_comms(ID)%name
    end if



  end subroutine seq_comm_setptrs
!---------------------------------------------------------
  subroutine seq_comm_setnthreads(nthreads)

    implicit none
    integer,intent(in) :: nthreads
    character(*),parameter :: subName =   '(seq_comm_setnthreads) '


  end subroutine seq_comm_setnthreads
!---------------------------------------------------------
  integer function seq_comm_getnthreads()

    implicit none
    integer :: omp_get_num_threads
    character(*),parameter :: subName =   '(seq_comm_getnthreads) '

    seq_comm_getnthreads = -1

  end function seq_comm_getnthreads
!---------------------------------------------------------
  logical function seq_comm_iamin(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iamin) '

    if (seq_comms(ID)%iam >= 0) then
       seq_comm_iamin = .true.
    else
       seq_comm_iamin = .false.
    endif

  end function seq_comm_iamin
!---------------------------------------------------------
  logical function seq_comm_iamroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iamroot) '

    seq_comm_iamroot = seq_comms(ID)%iamroot

  end function seq_comm_iamroot
!---------------------------------------------------------
  integer function seq_comm_mpicom(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_mpicom) '

    seq_comm_mpicom = seq_comms(ID)%mpicom

  end function seq_comm_mpicom
!---------------------------------------------------------
  integer function seq_comm_iam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iam) '

    seq_comm_iam = seq_comms(ID)%iam

  end function seq_comm_iam
!---------------------------------------------------------
  integer function seq_comm_gloiam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloiam) '

    seq_comm_gloiam = seq_comms(ID)%gloiam

  end function seq_comm_gloiam
!---------------------------------------------------------
  integer function seq_comm_gloroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloroot) '

    seq_comm_gloroot = seq_comms(ID)%gloroot

  end function seq_comm_gloroot
!---------------------------------------------------------
  integer function seq_comm_cplpe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cplpe) '

    seq_comm_cplpe = seq_comms(ID)%cplpe

  end function seq_comm_cplpe
!---------------------------------------------------------
  integer function seq_comm_cmppe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cmppe) '

    seq_comm_cmppe = seq_comms(ID)%cmppe

  end function seq_comm_cmppe
!---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_name(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_name) '

    seq_comm_name = trim(seq_comms(ID)%name)

  end function seq_comm_name
!---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_suffix(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_suffix) '

    seq_comm_suffix = trim(seq_comms(ID)%suffix)

  end function seq_comm_suffix
!---------------------------------------------------------
  subroutine seq_comm_petlist(ID,petlist)

    implicit none
    integer,intent(in) :: ID
    integer,pointer :: petlist(:)
    character(*),parameter :: subName =   '(seq_comm_petlist) '

    petlist => seq_comms(ID)%petlist

  end subroutine seq_comm_petlist
!---------------------------------------------------------
!---------------------------------------------------------
  integer function seq_comm_inst(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_inst) '

    seq_comm_inst = seq_comms(ID)%inst

  end function seq_comm_inst
!---------------------------------------------------------
  subroutine seq_comm_mkname(oname,str1,num)
    implicit none
    character(len=*),intent(out) :: oname
    character(len=*),intent(in)  :: str1
    integer,intent(in)           :: num
    character(*),parameter :: subName =   '(seq_comm_mkname) '

    character(len=8) :: cnum

    write(cnum,'(i4.4)') num
    if (len_trim(str1) + len_trim(cnum) > len(oname)) then
       write(logunit,*) trim(subname),' ERROR in str lens ',len(oname),trim(str1),trim(cnum)
       call shr_sys_abort(trim(subname))
    endif
    oname = trim(str1)//trim(cnum)

  end subroutine seq_comm_mkname
!---------------------------------------------------------
end module seq_comm_mct
