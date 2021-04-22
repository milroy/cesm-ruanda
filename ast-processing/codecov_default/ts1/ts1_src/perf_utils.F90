module perf_utils

!----------------------------------------------------------------------- 
! 
! Purpose: This module supplies the csm_share and CAM utilities
!          needed by perf_mod.F90 (when the csm_share and CAM utilities
!          are not available).
! 
! Author:  P. Worley, October 2007
!
! $Id$
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private

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

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public perfutils_setunit
   public shr_sys_abort
   public shr_mpi_barrier
   public shr_file_getUnit
   public shr_file_freeUnit
   public find_group_name
   public to_lower
   public shr_mpi_bcast

   interface shr_mpi_bcast ; module procedure &
     shr_mpi_bcastl0, &
     shr_mpi_bcasti0
   end interface

!-----------------------------------------------------------------------
! Private interfaces ---------------------------------------------------
!-----------------------------------------------------------------------
   private shr_sys_flush
   private shr_mpi_chkerr
   private shr_mpi_abort

!-----------------------------------------------------------------------
!- include statements --------------------------------------------------
!-----------------------------------------------------------------------
!
! $Id: gptl.inc,v 1.44 2011-03-28 20:55:19 rosinski Exp $
!
! Author: Jim Rosinski
!
! GPTL header file to be included in user code. Values match
! their counterparts in gptl.h. See that file or man pages
! or web-based documenation for descriptions of each value
!
      integer GPTLsync_mpi
      integer GPTLwall
      integer GPTLcpu
      integer GPTLabort_on_error
      integer GPTLoverhead
      integer GPTLdepthlimit
      integer GPTLverbose
      integer GPTLnarrowprint
      integer GPTLpercent
      integer GPTLpersec
      integer GPTLmultiplex
      integer GPTLdopr_preamble
      integer GPTLdopr_threadsort
      integer GPTLdopr_multparent
      integer GPTLdopr_collision
      integer GPTLprint_method
      integer GPTLtablesize

      integer GPTL_IPC
      integer GPTL_CI
      integer GPTL_FPC
      integer GPTL_FPI
      integer GPTL_LSTPI
      integer GPTL_DCMRT
      integer GPTL_LSTPDCM
      integer GPTL_L2MRT
      integer GPTL_LSTPL2M
      integer GPTL_L3MRT

      integer GPTLnanotime
      integer GPTLmpiwtime
      integer GPTLclockgettime
      integer GPTLgettimeofday
      integer GPTLpapitime
      integer GPTLread_real_time

      integer GPTLfirst_parent
      integer GPTLlast_parent
      integer GPTLmost_frequent
      integer GPTLfull_tree

      parameter (GPTLsync_mpi       = 0)
      parameter (GPTLwall           = 1)
      parameter (GPTLcpu            = 2)
      parameter (GPTLabort_on_error = 3)
      parameter (GPTLoverhead       = 4)
      parameter (GPTLdepthlimit     = 5)
      parameter (GPTLverbose        = 6)
      parameter (GPTLnarrowprint    = 7)
      parameter (GPTLpercent        = 9)
      parameter (GPTLpersec         = 10)
      parameter (GPTLmultiplex      = 11)
      parameter (GPTLdopr_preamble  = 12)
      parameter (GPTLdopr_threadsort= 13)
      parameter (GPTLdopr_multparent= 14)
      parameter (GPTLdopr_collision = 15)
      parameter (GPTLprint_method   = 16)
      parameter (GPTLtablesize      = 50)

      parameter (GPTL_IPC           = 17)
      parameter (GPTL_CI            = 18)
      parameter (GPTL_FPC           = 19)
      parameter (GPTL_FPI           = 20)
      parameter (GPTL_LSTPI         = 21)
      parameter (GPTL_DCMRT         = 22)
      parameter (GPTL_LSTPDCM       = 23)
      parameter (GPTL_L2MRT         = 24)
      parameter (GPTL_LSTPL2M       = 25)
      parameter (GPTL_L3MRT         = 26)

      parameter (GPTLgettimeofday   = 1)
      parameter (GPTLnanotime       = 2)
      parameter (GPTLmpiwtime       = 4)
      parameter (GPTLclockgettime   = 5)
      parameter (GPTLpapitime       = 6)
      parameter (GPTLread_real_time = 3)

      parameter (GPTLfirst_parent   = 1)
      parameter (GPTLlast_parent    = 2)
      parameter (GPTLmost_frequent  = 3)
      parameter (GPTLfull_tree      = 4)

! Externals

      integer gptlsetoption
      integer gptlinitialize
      integer gptlstart
      integer gptlstart_handle
      integer gptlstartf
      integer gptlstartf_handle
      integer gptlstop
      integer gptlstop_handle
      integer gptlstopf
      integer gptlstopf_handle
      integer gptlstamp 
      integer gptlpr_set_append
      integer gptlpr_query_append
      integer gptlpr_set_write
      integer gptlpr_query_write
      integer gptlpr
      integer gptlpr_file
      integer gptlpr_summary
      integer gptlpr_summary_file
      integer gptlbarrier
      integer gptlreset 
      integer gptlfinalize
      integer gptlget_memusage
      integer gptlprint_memusage
      integer gptlenable
      integer gptldisable
      integer gptlsetutr
      integer gptlquery
      integer gptlquerycounters
      integer gptlget_wallclock
      integer gptlget_eventvalue
      integer gptlget_nregions
      integer gptlget_regionname
      integer gptl_papilibraryinit
      integer gptlevent_name_to_code
      integer gptlevent_code_to_name

      external gptlsetoption
      external gptlinitialize
      external gptlstart
      external gptlstart_handle
      external gptlstartf
      external gptlstartf_handle
      external gptlstop
      external gptlstop_handle
      external gptlstopf
      external gptlstopf_handle
      external gptlstamp 
      external gptlpr_set_append
      external gptlpr_query_append
      external gptlpr_set_write
      external gptlpr_query_write
      external gptlpr
      external gptlpr_file
      external gptlpr_summary
      external gptlpr_summary_file
      external gptlbarrier
      external gptlreset 
      external gptlfinalize
      external gptlget_memusage
      external gptlprint_memusage
      external gptlenable
      external gptldisable
      external gptlsetutr
      external gptlquery
      external gptlquerycounters
      external gptlget_wallclock
      external gptlget_eventvalue
      external gptlget_nregions
      external gptlget_regionname
      external gptl_papilibraryinit
      external gptlevent_name_to_code
      external gptlevent_code_to_name

!-----------------------------------------------------------------------
! Public data ---------------------------------------------------------
!-----------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! precision/kind constants (from csm_share/shr/shr_kind_mod.F90)
   !----------------------------------------------------------------------------
   integer,parameter,public :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
   integer,parameter,public :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   integer,parameter,public :: SHR_KIND_IN = kind(1)                ! native integer
   integer,parameter,public :: SHR_KIND_CL = 256                    ! long char
   integer,parameter,public :: SHR_KIND_CX = 512                    ! extra-long char

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------

   integer, parameter :: def_pu_logunit = 6                   ! default
   integer, private   :: pu_logunit = def_pu_logunit
                         ! unit number for log output

!=======================================================================
contains
!=======================================================================

!
!========================================================================
!
   subroutine perfutils_setunit(LogUnit)
!----------------------------------------------------------------------- 
! Purpose:  Set log unit number.
! Author:   P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   integer(SHR_KIND_IN), intent(IN) :: LogUnit  ! Unit number for log output
!-----------------------------------------------------------------------
   pu_logunit = LogUnit
!
   return
!
   end subroutine perfutils_setunit

!============== Routines from csm_share/shr/shr_sys_mod.F90 ============
!=======================================================================

! Subprogram not used SUBROUTINE shr_sys_abort(string)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    character(*)        ,optional :: string  ! error message string
! Subprogram not used 
! Subprogram not used    !----- local -----
! Subprogram not used    integer(SHR_KIND_IN) :: ierr
! Subprogram not used    logical              :: flag
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName =   '(shr_sys_abort) '
! Subprogram not used    character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: consistent stopping mechanism
! Subprogram not used ! (dumbed down from original shr_sys_mod.F90 version for use in perf_mod)
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call shr_sys_flush(pu_logunit)
! Subprogram not used 
! Subprogram not used    if ( present(string) ) then
! Subprogram not used       if (len_trim(string) > 0) then
! Subprogram not used          write(pu_logunit,*) trim(subName),' ERROR: ',trim(string)
! Subprogram not used       else
! Subprogram not used          write(pu_logunit,*) trim(subName),' ERROR '
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       write(pu_logunit,*) trim(subName),' ERROR '
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    write(pu_logunit,F00) 'WARNING: calling mpi_abort() and stopping'
! Subprogram not used    call shr_sys_flush(pu_logunit)
! Subprogram not used    call mpi_abort(MPI_COMM_WORLD,0,ierr)
! Subprogram not used    call shr_sys_flush(pu_logunit)
! Subprogram not used    call abort()
! Subprogram not used 
! Subprogram not used    stop
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_sys_abort

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_sys_flush(unit)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    !----- arguments -----
! Subprogram not used    integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit
! Subprogram not used 
! Subprogram not used    !----- formats -----
! Subprogram not used    character(*),parameter :: subName =   '(shr_sys_flush) '
! Subprogram not used    character(*),parameter :: F00     = "('(shr_sys_flush) ',4a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: an architecture independant system call
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call flush(unit)
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_sys_flush

!===============================================================================

!================== Routines from csm_share/shr/shr_mpi_mod.F90 ===============
!===============================================================================

SUBROUTINE shr_mpi_chkerr(rcode,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: rcode  ! input MPI error code
   character(*),         intent(in) :: string ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_chkerr) '
   character(MPI_MAX_ERROR_STRING)  :: lstring
   integer(SHR_KIND_IN)             :: len
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: layer on MPI error checking
!-------------------------------------------------------------------------------

   if (rcode /= MPI_SUCCESS) then
     call MPI_ERROR_STRING(rcode,lstring,len,ierr)
     write(pu_logunit,*) trim(subName),":",lstring(1:len)
     call shr_mpi_abort(string,rcode)
   endif

END SUBROUTINE shr_mpi_chkerr

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_mpi_abort(string,rcode)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    character(*),optional,intent(in)   :: string   ! message
! Subprogram not used    integer,optional,intent(in)        :: rcode    ! optional code
! Subprogram not used 
! Subprogram not used    !----- local ---
! Subprogram not used    character(*),parameter             :: subName = '(shr_mpi_abort) '
! Subprogram not used    integer(SHR_KIND_IN)               :: ierr
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: MPI abort
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if ( present(string) .and. present(rcode) ) then
! Subprogram not used       write(pu_logunit,*) trim(subName),":",trim(string),rcode
! Subprogram not used    endif
! Subprogram not used    call MPI_ABORT(MPI_COMM_WORLD,rcode,ierr)
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_mpi_abort

!===============================================================================
!===============================================================================

! Subprogram not used SUBROUTINE shr_mpi_barrier(comm,string)
! Subprogram not used 
! Subprogram not used    IMPLICIT none
! Subprogram not used 
! Subprogram not used    !----- arguments ---
! Subprogram not used    integer,intent(in)                 :: comm
! Subprogram not used    character(*),optional,intent(in)   :: string   ! message
! Subprogram not used 
! Subprogram not used    !----- local ---
! Subprogram not used    character(*),parameter             :: subName = '(shr_mpi_barrier) '
! Subprogram not used    integer(SHR_KIND_IN)               :: ierr
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! PURPOSE: MPI barrier
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    call MPI_BARRIER(comm,ierr)
! Subprogram not used    if (present(string)) then
! Subprogram not used      call shr_mpi_chkerr(ierr,subName//trim(string))
! Subprogram not used    else
! Subprogram not used      call shr_mpi_chkerr(ierr,subName)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used END SUBROUTINE shr_mpi_barrier

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti0(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcasti0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_BCAST(vec,lsize,MPI_INTEGER,0,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastl0(vec,comm,string)

   IMPLICIT none

   !----- arguments ---
   logical, intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastl0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a logical
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_BCAST(vec,lsize,MPI_LOGICAL,0,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastl0

!===============================================================================

!================== Routines from csm_share/shr/shr_file_mod.F90 ===============
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_getUnit -- Get a free FORTRAN unit number
!
! !DESCRIPTION: Get the next free FORTRAN unit number.
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!     2007-Oct-21 - P. Worley - dumbed down for use in perf_mod
!
! !INTERFACE: ------------------------------------------------------------------  

INTEGER FUNCTION shr_file_getUnit ()

   implicit none

!EOP

   !----- local parameters -----
   integer(SHR_KIND_IN),parameter :: shr_file_minUnit = 10      ! Min unit number to give
   integer(SHR_KIND_IN),parameter :: shr_file_maxUnit = 99      ! Max unit number to give

   !----- local variables -----
   integer(SHR_KIND_IN)   :: n      ! loop index
   logical                :: opened ! If unit opened or not

   !----- formats -----
   character(*),parameter :: subName = '(shr_file_getUnit) '
   character(*),parameter :: F00   = "('(shr_file_getUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   ! --- Choose first available unit other than 0, 5, or 6  ------
   do n=shr_file_minUnit, shr_file_maxUnit
      inquire( n, opened=opened )
      if (n == 5 .or. n == 6 .or. opened) then
         cycle
      end if
      shr_file_getUnit = n
      return
   end do

   call shr_sys_abort( subName//': Error: no available units found' )

END FUNCTION shr_file_getUnit
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_freeUnit -- Free up a FORTRAN unit number
!
! !DESCRIPTION: Free up the given unit number
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!     2007-Oct-21 - P. Worley - dumbed down for use in perf_mod
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_freeUnit ( unit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: unit  ! unit number to be freed

!EOP

   !----- local parameters -----
   integer(SHR_KIND_IN),parameter :: shr_file_minUnit = 10      ! Min unit number to give
   integer(SHR_KIND_IN),parameter :: shr_file_maxUnit = 99      ! Max unit number to give

   !----- formats -----
   character(*), parameter :: subName = '(shr_file_freeUnit) '
   character(*), parameter :: F00 =   "('(shr_file_freeUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (unit < 0 .or. unit > shr_file_maxUnit) then
!pw   if (s_loglev > 0) write(pu_logunit,F00) 'invalid unit number request:', unit
   else if (unit == 0 .or. unit == 5 .or. unit == 6) then
      call shr_sys_abort( subName//': Error: units 0, 5, and 6 must not be freed' )
   end if

   return

END SUBROUTINE shr_file_freeUnit
!===============================================================================

!============= Routines from atm/cam/src/utils/namelist_utils.F90 ==============
!===============================================================================

subroutine find_group_name(unit, group, status)

!---------------------------------------------------------------------------------------
! Purpose: 
! Search a file that contains namelist input for the specified namelist group name.
! Leave the file positioned so that the current record is the first record of the
! input for the specified group.
! 
! Method: 
! Read the file line by line.  Each line is searched for an '&' which may only
! be preceded by blanks, immediately followed by the group name which is case
! insensitive.  If found then backspace the file so the current record is the
! one containing the group name and return success.  Otherwise return -1.
!
! Author:  B. Eaton, August 2007
!---------------------------------------------------------------------------------------

   integer,          intent(in)  :: unit     ! fortran unit attached to file
   character(len=*), intent(in)  :: group    ! namelist group name
   integer,          intent(out) :: status   ! 0 for success, -1 if group name not found

   ! Local variables

   integer           :: len_grp
   integer           :: ios    ! io status
   character(len=80) :: inrec  ! first 80 characters of input record
   character(len=80) :: inrec2 ! left adjusted input record
   character(len=len(group)) :: lc_group

   !---------------------------------------------------------------------------

   len_grp = len_trim(group)
   lc_group = to_lower(group)

   ios = 0
   do while (ios <= 0)

      read(unit, '(a)', iostat=ios, end=102) inrec

      if (ios <= 0) then  ! ios < 0  indicates an end of record condition

         ! look for group name in this record

         ! remove leading blanks
         inrec2 = to_lower(adjustl(inrec))

         ! check for leading '&'
         if (inrec2(1:1) == '&') then

            ! check for case insensitive group name
            if (trim(lc_group) == inrec2(2:len_grp+1)) then

               ! found group name.  backspace to leave file position at this record
               backspace(unit)
               status = 0
               return

            end if
         end if
      end if

   end do

   102 continue  ! end of file processing
   status = -1

end subroutine find_group_name
!===============================================================================

!================ Routines from atm/cam/src/utils/string_utils.F90 =============
!===============================================================================

function to_lower(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to lower case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id$
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to lower case
   character(len=len(str))      :: to_lower

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: upper_to_lower   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
!-----------------------------------------------------------------------
   upper_to_lower = iachar("a") - iachar("A")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
           ctmp = achar(aseq + upper_to_lower)	
      to_lower(i:i) = ctmp
   end do

end function to_lower
!===============================================================================

end module perf_utils
