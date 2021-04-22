module perf_mod

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for controlling the performance
!          timer logic.
! 
! Author:  P. Worley, January 2007
!
! $Id$
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- Uses ----------------------------------------------------------------
!-----------------------------------------------------------------------

   use perf_utils

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
   public t_initf
   public t_setLogUnit
   public t_getLogUnit
   public t_profile_onf
   public t_barrier_onf
   public t_single_filef
   public t_stampf
   public t_startf
   public t_stopf
   public t_enablef
   public t_disablef
   public t_adj_detailf
   public t_barrierf
   public t_prf
   public t_finalizef

!-----------------------------------------------------------------------
! Private interfaces (local) -------------------------------------------
!-----------------------------------------------------------------------
   private perf_defaultopts
   private perf_setopts
   private papi_defaultopts
   private papi_setopts

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
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! perf_mod options
   !----------------------------------------------------------------------------
   integer, parameter :: def_p_logunit = 6                   ! default
   integer, private   :: p_logunit = def_p_logunit
                         ! unit number for log output

   logical, parameter :: def_timing_initialized = .false.      ! default
   logical, private   :: timing_initialized = def_timing_initialized
                         ! flag indicating whether timing library has
                         ! been initialized

   logical, parameter :: def_timing_disable = .false.          ! default
   logical, private   :: timing_disable = def_timing_disable
                         ! flag indicating whether timers are disabled

   logical, parameter :: def_timing_barrier = .false.          ! default
   logical, private   :: timing_barrier = def_timing_barrier
                         ! flag indicating whether the mpi_barrier in
                         ! t_barrierf should be called

   integer, parameter :: def_timer_depth_limit = 99999         ! default
   integer, private   :: timer_depth_limit = def_timer_depth_limit
                         ! integer indicating maximum number of levels of
                         ! timer nesting 

   integer, parameter :: def_timing_detail_limit = 1           ! default
   integer, private   :: timing_detail_limit = def_timing_detail_limit
                         ! integer indicating maximum detail level to
                         ! profile

   integer, parameter :: init_timing_disable_depth = 0         ! init
   integer, private   :: timing_disable_depth = init_timing_disable_depth
                         ! integer indicating depth of t_disablef calls

   integer, parameter :: init_timing_detail = 0                ! init
   integer, private   :: cur_timing_detail = init_timing_detail
                         ! current timing detail level

   logical, parameter :: def_perf_single_file = .false.         ! default
   logical, private   :: perf_single_file = def_perf_single_file
                         ! flag indicating whether the performance timer
                         ! output should be written to a single file 
                         ! (per component communicator) or to a 
                         ! separate file for each process

   integer, parameter :: def_perf_outpe_num = 0                ! default
   integer, private   :: perf_outpe_num = def_perf_outpe_num
                         ! maximum number of processes writing out 
                         ! timing data (for this component communicator)

   integer, parameter :: def_perf_outpe_stride = 1             ! default
   integer, private   :: perf_outpe_stride = def_perf_outpe_stride
                         ! separation between process ids for processes
                         ! that are writing out timing data 
                         ! (for this component communicator)

   logical, parameter :: def_perf_global_stats = .true.        ! default
   logical, private   :: perf_global_stats = def_perf_global_stats
                         ! collect and print out global performance statistics
                         ! (for this component communicator)
   integer, parameter :: def_perf_timer = GPTLmpiwtime         ! default


   integer, private   :: perf_timer = def_perf_timer           ! default
                         ! integer indicating which timer to use
                         ! (as defined in gptl.inc)

   logical, parameter :: def_perf_papi_enable = .false.       ! default
   logical, private   :: perf_papi_enable = def_perf_papi_enable
                         ! flag indicating whether the PAPI namelist
                         ! should be read and HW performance counters
                         ! used in profiling

   ! PAPI counter ids
   integer, parameter :: PAPI_NULL = -1

   integer, parameter :: def_papi_ctr1 = PAPI_NULL           ! default
   integer, private   :: papi_ctr1 = def_papi_ctr1

   integer, parameter :: def_papi_ctr2 = PAPI_NULL           ! default
   integer, private   :: papi_ctr2 = def_papi_ctr2

   integer, parameter :: def_papi_ctr3 = PAPI_NULL           ! default
   integer, private   :: papi_ctr3 = def_papi_ctr3

   integer, parameter :: def_papi_ctr4 = PAPI_NULL           ! default
   integer, private   :: papi_ctr4 = def_papi_ctr4

!=======================================================================
contains
!=======================================================================

!
!========================================================================
!
! Subprogram not used    subroutine t_getLogUnit(LogUnit)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose:  Get log unit number.
! Subprogram not used ! Author:   P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !---------------------------Input arguments-----------------------------
! Subprogram not used !
! Subprogram not used    integer(SHR_KIND_IN), intent(OUT) :: LogUnit  ! Unit number for log output
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    LogUnit = p_logunit
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine t_getLogUnit
!
!========================================================================
!
   subroutine t_setLogUnit(LogUnit)
!----------------------------------------------------------------------- 
! Purpose:  Set log unit number.
! Author:   P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   integer(SHR_KIND_IN), intent(IN) :: LogUnit  ! Unit number for log output
!-----------------------------------------------------------------------

   p_logunit = LogUnit
   call perfutils_setunit(p_logunit)

   return
   end subroutine t_setLogUnit
!
!========================================================================
!
   subroutine perf_defaultopts(timing_disable_out, &
                               perf_timer_out, &
                               timer_depth_limit_out, &
                               timing_detail_limit_out, &
                               timing_barrier_out, &
                               perf_outpe_num_out, &
                               perf_outpe_stride_out, &
                               perf_single_file_out, &
                               perf_global_stats_out, &
                               perf_papi_enable_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
   ! timers disable/enable option
   logical, intent(out), optional :: timing_disable_out
   ! performance timer option
   integer, intent(out), optional :: perf_timer_out
   ! timer depth limit option
   integer, intent(out), optional :: timer_depth_limit_out
   ! timer detail limit option
   integer, intent(out), optional :: timing_detail_limit_out
   ! timing barrier enable/disable option
   logical, intent(out), optional :: timing_barrier_out
   ! number of processes writing out timing data
   integer, intent(out), optional :: perf_outpe_num_out
   ! separation between process ids for processes that are writing out timing data
   integer, intent(out), optional :: perf_outpe_stride_out
   ! timing single / multple output file option
   logical, intent(out), optional :: perf_single_file_out
   ! collect and output global performance statistics option
   logical, intent(out), optional :: perf_global_stats_out
   ! calling PAPI to read HW performance counters option
   logical, intent(out), optional :: perf_papi_enable_out
!-----------------------------------------------------------------------
   if ( present(timing_disable_out) ) then
      timing_disable_out = def_timing_disable
   endif
   if ( present(perf_timer_out) ) then
      perf_timer_out = def_perf_timer
   endif
   if ( present(timer_depth_limit_out) ) then
      timer_depth_limit_out = def_timer_depth_limit
   endif
   if ( present(timing_detail_limit_out) ) then
      timing_detail_limit_out = def_timing_detail_limit
   endif
   if ( present(timing_barrier_out) ) then
      timing_barrier_out = def_timing_barrier
   endif
   if ( present(perf_outpe_num_out) ) then
      perf_outpe_num_out = def_perf_outpe_num
   endif
   if ( present(perf_outpe_stride_out) ) then
      perf_outpe_stride_out = def_perf_outpe_stride
   endif
   if ( present(perf_single_file_out) ) then
      perf_single_file_out = def_perf_single_file
   endif
   if ( present(perf_global_stats_out) ) then
      perf_global_stats_out = def_perf_global_stats
   endif
   if ( present(perf_papi_enable_out) ) then
      perf_papi_enable_out = def_perf_papi_enable
   endif
!
   return
   end subroutine perf_defaultopts
!
!========================================================================
!
   subroutine perf_setopts(mastertask, &
                           LogPrint, &
                           timing_disable_in, &
                           perf_timer_in, &
                           timer_depth_limit_in, &
                           timing_detail_limit_in, &
                           timing_barrier_in, &
                           perf_outpe_num_in, &
                           perf_outpe_stride_in, &
                           perf_single_file_in, &
                           perf_global_stats_in, &
                           perf_papi_enable_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments----------------------------
!
   ! master process?
   logical, intent(in) :: mastertask
   ! Print out to log file?
   logical, intent(IN) :: LogPrint        
   ! timers disable/enable option
   logical, intent(in), optional :: timing_disable_in
   ! performance timer option
   integer, intent(in), optional :: perf_timer_in
   ! timer depth limit option
   integer, intent(in), optional :: timer_depth_limit_in
   ! timer detail limit option
   integer, intent(in), optional :: timing_detail_limit_in
   ! timing barrier enable/disable option
   logical, intent(in), optional :: timing_barrier_in
   ! number of processes writing out timing data
   integer, intent(in), optional :: perf_outpe_num_in
   ! separation between process ids for processes that are writing out timing data
   integer, intent(in), optional :: perf_outpe_stride_in
   ! timing single / multple output file option
   logical, intent(in), optional :: perf_single_file_in
   ! collect and output global performance statistics option
   logical, intent(in), optional :: perf_global_stats_in
   ! calling PAPI to read HW performance counters option
   logical, intent(in), optional :: perf_papi_enable_in
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! error return
!-----------------------------------------------------------------------
   if ( .not. timing_initialized ) then

      if ( present(timing_disable_in) ) then
         timing_disable = timing_disable_in
         if (timing_disable) then
            ierr = GPTLdisable()
         else 
            ierr = GPTLenable()
         endif
      endif
      if ( present(perf_timer_in) ) then
         if ((perf_timer_in .eq. GPTLgettimeofday) .or. &
             (perf_timer_in .eq. GPTLnanotime) .or. &
             (perf_timer_in .eq. GPTLread_real_time) .or. &
             (perf_timer_in .eq. GPTLmpiwtime) .or. &
             (perf_timer_in .eq. GPTLclockgettime) .or. &
             (perf_timer_in .eq. GPTLpapitime)) then
            perf_timer = perf_timer_in
         else
            if (mastertask) then
               write(p_logunit,*) 'PERF_SETOPTS: illegal timer requested=',&
                                  perf_timer_in, '. Request ignored.'
            endif
         endif
      endif
      if ( present(timer_depth_limit_in) ) then
         timer_depth_limit = timer_depth_limit_in
      endif
      if ( present(timing_detail_limit_in) ) then
         timing_detail_limit = timing_detail_limit_in
      endif
      if ( present(timing_barrier_in) ) then
         timing_barrier = timing_barrier_in
      endif
      if ( present(perf_outpe_num_in) ) then
         perf_outpe_num = perf_outpe_num_in
      endif
      if ( present(perf_outpe_stride_in) ) then
         perf_outpe_stride = perf_outpe_stride_in
      endif
      if ( present(perf_single_file_in) ) then
         perf_single_file = perf_single_file_in
      endif
      if ( present(perf_global_stats_in) ) then
         perf_global_stats = perf_global_stats_in
      endif
      if ( present(perf_papi_enable_in) ) then
         if (perf_papi_enable_in) then
            if (mastertask) then
               write(p_logunit,*) 'PERF_SETOPTS: PAPI library not linked in. ',&
                                  'Request to enable PAPI ignored.'
            endif
         endif
         perf_papi_enable = .false.
      endif
!
      if (mastertask .and. LogPrint) then
         write(p_logunit,*) '(t_initf) Using profile_disable=', timing_disable, &             
                            ' profile_timer=', perf_timer
         write(p_logunit,*) '(t_initf)  profile_depth_limit=', timer_depth_limit, &    
                            ' profile_detail_limit=', timing_detail_limit
         write(p_logunit,*) '(t_initf)  profile_barrier=', timing_barrier, &
                            ' profile_outpe_num=', perf_outpe_num
         write(p_logunit,*) '(t_initf)  profile_outpe_stride=', perf_outpe_stride , &
                            ' profile_single_file=', perf_single_file
         write(p_logunit,*) '(t_initf)  profile_global_stats=', perf_global_stats , &
                            ' profile_papi_enable=', perf_papi_enable 
      endif                                                                               
!
   endif
!
   return
   end subroutine perf_setopts

!
!========================================================================
!
! Subprogram not used    subroutine papi_defaultopts(papi_ctr1_out, &
! Subprogram not used                                papi_ctr2_out, &
! Subprogram not used                                papi_ctr3_out, &
! Subprogram not used                                papi_ctr4_out  )
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Return default runtime PAPI counter options
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !---------------------------Input arguments-----------------------------
! Subprogram not used    ! PAPI counter option #1
! Subprogram not used    integer, intent(out), optional :: papi_ctr1_out
! Subprogram not used    ! PAPI counter option #2
! Subprogram not used    integer, intent(out), optional :: papi_ctr2_out
! Subprogram not used    ! PAPI counter option #3
! Subprogram not used    integer, intent(out), optional :: papi_ctr3_out
! Subprogram not used    ! PAPI counter option #4
! Subprogram not used    integer, intent(out), optional :: papi_ctr4_out
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if ( present(papi_ctr1_out) ) then
! Subprogram not used       papi_ctr1_out = def_papi_ctr1
! Subprogram not used    endif
! Subprogram not used    if ( present(papi_ctr2_out) ) then
! Subprogram not used       papi_ctr2_out = def_papi_ctr2
! Subprogram not used    endif
! Subprogram not used    if ( present(papi_ctr3_out) ) then
! Subprogram not used       papi_ctr3_out = def_papi_ctr3
! Subprogram not used    endif
! Subprogram not used    if ( present(papi_ctr4_out) ) then
! Subprogram not used       papi_ctr4_out = def_papi_ctr4
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine papi_defaultopts
!
!========================================================================
!
! Subprogram not used    subroutine papi_setopts(papi_ctr1_in, &
! Subprogram not used                            papi_ctr2_in, &
! Subprogram not used                            papi_ctr3_in, &
! Subprogram not used                            papi_ctr4_in  )
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Set runtime PAPI counter options
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !---------------------------Input arguments----------------------------
! Subprogram not used !
! Subprogram not used    ! performance counter option
! Subprogram not used    integer, intent(in), optional :: papi_ctr1_in
! Subprogram not used    ! performance counter option
! Subprogram not used    integer, intent(in), optional :: papi_ctr2_in
! Subprogram not used    ! performance counter option
! Subprogram not used    integer, intent(in), optional :: papi_ctr3_in
! Subprogram not used    ! performance counter option
! Subprogram not used    integer, intent(in), optional :: papi_ctr4_in
! Subprogram not used !
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    integer  ierr                  ! error return
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used    if ( .not. timing_initialized ) then
! Subprogram not used 
! Subprogram not used       if ( present(papi_ctr1_in) ) then
! Subprogram not used          if (papi_ctr1_in < 0) then
! Subprogram not used             papi_ctr1 = papi_ctr1_in
! Subprogram not used          else
! Subprogram not used             papi_ctr1 = PAPI_NULL
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used       if ( present(papi_ctr2_in) ) then
! Subprogram not used          if (papi_ctr2_in < 0) then
! Subprogram not used             papi_ctr2 = papi_ctr2_in
! Subprogram not used          else
! Subprogram not used             papi_ctr2 = PAPI_NULL
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used       if ( present(papi_ctr3_in) ) then
! Subprogram not used          if (papi_ctr3_in < 0) then
! Subprogram not used             papi_ctr3 = papi_ctr3_in
! Subprogram not used          else
! Subprogram not used             papi_ctr3 = PAPI_NULL
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used       if ( present(papi_ctr4_in) ) then
! Subprogram not used          if (papi_ctr4_in < 0) then
! Subprogram not used             papi_ctr4 = papi_ctr4_in
! Subprogram not used          else
! Subprogram not used             papi_ctr4 = PAPI_NULL
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used    endif
! Subprogram not used !
! Subprogram not used    return
! Subprogram not used    end subroutine papi_setopts
!
!========================================================================
!
! Subprogram not used    logical function t_profile_onf()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Return flag indicating whether profiling is currently active.
! Subprogram not used !          Part of workaround to implement FVbarrierclock before
! Subprogram not used !          communicators exposed in Pilgrim. Does not check level of
! Subprogram not used !          event nesting.
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if ((.not. timing_initialized) .or. &
! Subprogram not used        (timing_disable_depth > 0) .or. &
! Subprogram not used        (cur_timing_detail > timing_detail_limit)) then
! Subprogram not used       t_profile_onf = .false.
! Subprogram not used    else
! Subprogram not used       t_profile_onf = .true.
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    end function t_profile_onf
!
!========================================================================
!
! Subprogram not used    logical function t_barrier_onf()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Return timing_barrier. Part of workaround to implement 
! Subprogram not used !          FVbarrierclock before communicators exposed in Pilgrim. 
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    t_barrier_onf = timing_barrier
! Subprogram not used 
! Subprogram not used    end function t_barrier_onf
!
!========================================================================
!
! Subprogram not used    logical function t_single_filef()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Return perf_single_file. Used to control output of other
! Subprogram not used !          performance data, only spmdstats currently.
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    t_single_filef = perf_single_file
! Subprogram not used 
! Subprogram not used    end function t_single_filef
!
!========================================================================
!
! Subprogram not used    subroutine t_stampf(wall, usr, sys)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Record wallclock, user, and system times (seconds).
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !---------------------------Output arguments-----------------------------
! Subprogram not used !
! Subprogram not used    real(shr_kind_r8), intent(out) :: wall ! wallclock time
! Subprogram not used    real(shr_kind_r8), intent(out) :: usr  ! user time
! Subprogram not used    real(shr_kind_r8), intent(out) :: sys  ! system time
! Subprogram not used !
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    integer  ierr                          ! GPTL error return
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    if ((.not. timing_initialized) .or. &
! Subprogram not used        (timing_disable_depth > 0)) then
! Subprogram not used       wall = 0.0
! Subprogram not used       usr = 0.0
! Subprogram not used       sys = 0.0
! Subprogram not used    else
! Subprogram not used       ierr = GPTLstamp(wall, usr, sys)
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine t_stampf
!
!========================================================================
!
   subroutine t_startf(event, handle)
!----------------------------------------------------------------------- 
! Purpose: Start an event timer
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   ! performance timer event name
   character(len=*), intent(in) :: event  
!
!---------------------------Input/Output arguments----------------------
!
   ! GPTL event handle
   integer(shr_kind_i8), optional :: handle
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if ((timing_initialized) .and. &
       (timing_disable_depth .eq. 0) .and. &
       (cur_timing_detail .le. timing_detail_limit)) then

      if ( present (handle) ) then
         ierr = GPTLstart_handle(event, handle)
      else
         ierr = GPTLstart(event)
      endif

   endif

   return
   end subroutine t_startf
!
!========================================================================
!
   subroutine t_stopf(event, handle)
!----------------------------------------------------------------------- 
! Purpose: Stop an event timer
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   ! performance timer event name
   character(len=*), intent(in) :: event  
!
!---------------------------Input/Output arguments----------------------
!
   ! GPTL event handle
   integer(shr_kind_i8), optional :: handle
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if ((timing_initialized) .and. &
       (timing_disable_depth .eq. 0) .and. &
       (cur_timing_detail .le. timing_detail_limit)) then

      if ( present (handle) ) then
         ierr = GPTLstop_handle(event, handle)
      else
         ierr = GPTLstop(event)
      endif

   endif

   return
   end subroutine t_stopf
!
!========================================================================
!
! Subprogram not used    subroutine t_enablef()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Enable t_startf, t_stopf, t_stampf, and t_barrierf. Ignored
! Subprogram not used !          in threaded regions.
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    integer  ierr                  ! GPTL error return
! Subprogram not used !
! Subprogram not used !---------------------------Externals-----------------------------------
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    if (.not. timing_initialized) return
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    if (timing_disable_depth > 0) then
! Subprogram not used       if (timing_disable_depth .eq. 1) then
! Subprogram not used          ierr = GPTLenable()
! Subprogram not used       endif
! Subprogram not used       timing_disable_depth = timing_disable_depth - 1
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine t_enablef
!
!========================================================================
!
! Subprogram not used    subroutine t_disablef()
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! Purpose: Disable t_startf, t_stopf, t_stampf, and t_barrierf. Ignored
! Subprogram not used !          in threaded regions.
! Subprogram not used ! Author: P. Worley 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !---------------------------Local workspace-----------------------------
! Subprogram not used !
! Subprogram not used    integer  ierr                  ! GPTL error return
! Subprogram not used !
! Subprogram not used !---------------------------Externals-----------------------------------
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used    if (.not. timing_initialized) return
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    if (timing_disable_depth .eq. 0) then
! Subprogram not used       ierr = GPTLdisable()
! Subprogram not used    endif
! Subprogram not used    timing_disable_depth = timing_disable_depth + 1
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used    end subroutine t_disablef
!
!========================================================================
!
   subroutine t_adj_detailf(detail_adjustment)
!----------------------------------------------------------------------- 
! Purpose: Modify current detail level. Ignored in threaded regions.
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   integer, intent(in) :: detail_adjustment ! user defined increase or
                                            ! decrease in detail level
!
!---------------------------Externals-----------------------------------
!
!
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return


   cur_timing_detail = cur_timing_detail + detail_adjustment

   return
   end subroutine t_adj_detailf
!
!========================================================================
!
   subroutine t_barrierf(event, mpicom)
!----------------------------------------------------------------------- 
! Purpose: Call (and time) mpi_barrier. Ignored inside OpenMP
!          threaded regions. Note that barrier executed even if
!          event not recorded because of level of timer event nesting.
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
   ! mpi communicator id
   integer, intent(in), optional :: mpicom
   ! performance timer event name
   character(len=*), intent(in), optional :: event
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! GPTL error return
!
!---------------------------Externals-----------------------------------
!
!
!-----------------------------------------------------------------------
!
   if ((timing_initialized) .and. &
       (timing_disable_depth .eq. 0) .and. &
       (cur_timing_detail .le. timing_detail_limit)) then

      if (timing_barrier) then

         if ( present (event) ) then
            ierr = GPTLstart(event)
         endif

         if ( present (mpicom) ) then
            call shr_mpi_barrier(mpicom, 'T_BARRIERF: bad mpi communicator')
         else
            call shr_mpi_barrier(MPI_COMM_WORLD, 'T_BARRIERF: bad mpi communicator')
         endif

         if ( present (event) ) then
            ierr = GPTLstop(event)
         endif

      endif

   endif

   return
   end subroutine t_barrierf
!
!========================================================================
!
   subroutine t_prf(filename, mpicom, num_outpe, stride_outpe, &
                    single_file, global_stats, output_thispe)
!----------------------------------------------------------------------- 
! Purpose: Write out performance timer data
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   ! performance timer output file name
   character(len=*), intent(in), optional :: filename
   ! mpi communicator id
   integer, intent(in), optional :: mpicom
   ! maximum number of processes writing out timing data
   integer, intent(in), optional :: num_outpe
   ! separation between process ids for processes writing out data 
   integer, intent(in), optional :: stride_outpe
   ! enable/disable the writing of data to a single file
   logical, intent(in), optional :: single_file
   ! enable/disable the collection of global statistics
   logical, intent(in), optional :: global_stats
   ! output timing data for this process
   logical, intent(in), optional :: output_thispe
!
!---------------------------Local workspace-----------------------------
!
   logical  one_file              ! flag indicting whether to write
                                  !  all data to a single file
   logical  glb_stats             ! flag indicting whether to compute
                                  !  global statistics
   logical  pr_write              ! flag indicating whether the current 
                                  !  GPTL output mode is write
   logical  write_data            ! flag indicating whether this process
                                  !  should output its timing data
   integer  i                     ! loop index
   integer  mpicom2               ! local copy of MPI communicator
   integer  me                    ! communicator local process id
   integer  npes                  ! local communicator group size
   integer  gme                   ! global process id
   integer  ierr                  ! MPI error return
   integer  outpe_num             ! max number of processes writing out
                                  !  timing data (excluding output_thispe)
   integer  outpe_stride          ! separation between process ids for
                                  !  processes writing out timing data
   integer  max_outpe             ! max process id for processes
                                  !  writing out timing data
   integer  signal                ! send/recv variable for single
                                  ! output file logic
   integer  str_length            ! string length
   integer  unitn                 ! file unit number
   integer cme_adj                ! length of filename suffix
   integer status (MPI_STATUS_SIZE)    ! Status of message
   character(len=7) cme                ! string representation of process id
   character(len=SHR_KIND_CX+14) fname ! timing output filename
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

   call t_startf("t_prf")
!$OMP MASTER
   call mpi_comm_rank(MPI_COMM_WORLD, gme, ierr)
   if ( present(mpicom) ) then
      mpicom2 = mpicom
      call mpi_comm_size(mpicom2, npes, ierr)
         if (ierr .eq. MPI_ERR_COMM) then
            call shr_sys_abort('T_PRF: bad mpi communicator')
         endif
      call mpi_comm_rank(mpicom2, me, ierr)
   else
      call mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
      mpicom2 = MPI_COMM_WORLD
      me = gme
   endif

   do i=1,SHR_KIND_CX+14
     fname(i:i) = " "
   enddo

   unitn = shr_file_getUnit()

   ! determine what the current output mode is (append or write)
   if (GPTLpr_query_write() == 1) then
     pr_write = .true.
     ierr = GPTLpr_set_append()
   else 
     pr_write=.false.
   endif

   ! Determine whether to write all data to a single fie
   if (present(single_file)) then
      one_file = single_file
   else
      one_file = perf_single_file
   endif

   ! Determine whether to compute global statistics
   if (present(global_stats)) then
      glb_stats = global_stats
   else
      glb_stats = perf_global_stats
   endif

   ! Determine which processes are writing out timing data
   write_data = .false.

   if (present(num_outpe)) then
      if (num_outpe < 0) then
         outpe_num = npes
      else
         outpe_num = num_outpe
      endif
   else
      if (perf_outpe_num < 0) then
         outpe_num = npes
      else
         outpe_num = perf_outpe_num
      endif
   endif

   if (present(stride_outpe)) then
      if (stride_outpe < 1) then
         outpe_stride = 1
      else
         outpe_stride = stride_outpe
      endif
   else
      if (perf_outpe_stride < 1) then
         outpe_stride = 1
      else
         outpe_stride = perf_outpe_stride
      endif
   endif

   max_outpe = min(outpe_num*outpe_stride, npes) - 1

   if ((mod(me, outpe_stride) .eq. 0) .and. (me .le. max_outpe)) &
      write_data = .true.

   if (present(output_thispe)) then
      write_data = output_thispe
   endif

   ! If a single timing output file, take turns writing to it.
   if (one_file) then

      if ( present(filename) ) then
         str_length = min(SHR_KIND_CX,len_trim(filename))
         fname(1:str_length) = filename(1:str_length)
      else
         fname(1:10) = "timing_all"
      endif

      signal = 0
      if (me .eq. 0) then

         if (glb_stats) then
            open( unitn, file=trim(fname), status='UNKNOWN' )
            write( unitn, 100) npes
 100        format(/,"***** GLOBAL STATISTICS (",I6," MPI TASKS) *****",/)
            close( unitn )

            ierr = GPTLpr_summary_file(mpicom2, trim(fname))
         endif

         if (write_data) then
            if (glb_stats) then
               open( unitn, file=trim(fname), status='OLD', position='APPEND' )
            else
               open( unitn, file=trim(fname), status='UNKNOWN' )
            endif

            write( unitn, 101) me, gme
 101        format(/,"************ PROCESS ",I6," (",I6,") ************",/)
            close( unitn )

            ierr = GPTLpr_file(trim(fname))
         endif

      else

         if (glb_stats) then
            ierr = GPTLpr_summary_file(mpicom2, trim(fname))
         endif

         call mpi_recv (signal, 1, mpi_integer, me-1, me-1, mpicom2, status, ierr)
         if (ierr /= mpi_success) then
            write(p_logunit,*) 'T_PRF: mpi_recv failed ierr=',ierr
            call shr_sys_abort()
         end if

         if (write_data) then
            open( unitn, file=trim(fname), status='OLD', position='APPEND' )
            write( unitn, 101) me, gme
            close( unitn )

            ierr = GPTLpr_file(trim(fname))
         endif

      endif

      if (me+1 < npes) &
         call mpi_send (signal, 1, mpi_integer, me+1, me, mpicom2, ierr)

   else

      if (glb_stats) then
         if ( present(filename) ) then
            str_length = min(SHR_KIND_CX-6,len_trim(filename))
            fname(1:str_length) = filename(1:str_length)
         else
            str_length = 6
            fname(1:10) = "timing"
         endif
         fname(str_length+1:str_length+6) = '_stats'

         if (me .eq. 0) then
            open( unitn, file=trim(fname), status='UNKNOWN' )
            write( unitn, 100) npes
            close( unitn )
         endif

         ierr = GPTLpr_summary_file(mpicom2, trim(fname))
         fname(str_length+1:str_length+6) = '      '
      endif

      if (write_data) then
         if (npes .le. 10) then
            write(cme,'(i1.1)') me
            cme_adj = 2
         elseif (npes .le. 100) then
            write(cme,'(i2.2)') me
            cme_adj = 3
         elseif (npes .le. 1000) then
            write(cme,'(i3.3)') me
            cme_adj = 4
         elseif (npes .le. 10000) then
            write(cme,'(i4.4)') me
            cme_adj = 5
         elseif (npes .le. 100000) then
            write(cme,'(i5.5)') me
            cme_adj = 6
         else
            write(cme,'(i6.6)') me
            cme_adj = 7
         endif

         if ( present(filename) ) then
            str_length = min(SHR_KIND_CX-cme_adj,len_trim(filename))
            fname(1:str_length) = filename(1:str_length)
         else
            str_length = 6
            fname(1:10) = "timing"
         endif
         fname(str_length+1:str_length+1) = '.'
         fname(str_length+2:str_length+cme_adj) = cme

         open( unitn, file=trim(fname), status='UNKNOWN' )
         write( unitn, 101) me, gme
         close( unitn )

         ierr = GPTLpr_file(trim(fname))
      endif

   endif

   call shr_file_freeUnit( unitn )

   ! reset GPTL output mode
   if (pr_write) then
     ierr = GPTLpr_set_write()
   endif

!$OMP END MASTER
   call t_stopf("t_prf")

   return
   end subroutine t_prf
!
!========================================================================
!
   subroutine t_initf(NLFilename, LogPrint, LogUnit, mpicom, MasterTask)
!----------------------------------------------------------------------- 
! Purpose:  Set default values of runtime timing options 
!           before namelists prof_inparm and papi_inparm are read,
!           read namelists (and broadcast, if SPMD),
!           then initialize timing library.
! Author:   P. Worley (based on shr_inputinfo_mod and runtime_opts)
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   character(len=*),   intent(IN) :: NLFilename      ! Name-list filename
   logical, optional,  intent(IN) :: LogPrint        ! If print out to log file
   integer, optional,  intent(IN) :: LogUnit         ! Unit number for log output
   integer, optional,  intent(IN) :: mpicom          ! MPI communicator
   logical, optional,  intent(IN) :: MasterTask      ! If MPI master task
!
!---------------------------Local workspace-----------------------------
!
   character(len=*), parameter    :: subname = '(T_INITF) '
   logical                        :: MasterTask2     ! If MPI master task
   logical                        :: LogPrint2       ! If print to log

   integer  me                    ! communicator local process id
   integer  ierr                  ! error return
   integer  unitn                 ! file unit number
   integer  papi_ctr1_id          ! PAPI counter id
   integer  papi_ctr2_id          ! PAPI counter id
   integer  papi_ctr3_id          ! PAPI counter id
   integer  papi_ctr4_id          ! PAPI counter id
!
!---------------------------Namelists ----------------------------------
!
   logical profile_disable
   logical profile_barrier
   logical profile_single_file
   logical profile_global_stats
   integer profile_depth_limit
   integer profile_detail_limit
   integer profile_outpe_num
   integer profile_outpe_stride
   integer profile_timer
   logical profile_papi_enable
   namelist /prof_inparm/ profile_disable, profile_barrier, &
                          profile_single_file, profile_global_stats, &
                          profile_depth_limit, &
                          profile_detail_limit, profile_outpe_num, &
                          profile_outpe_stride, profile_timer, &
                          profile_papi_enable

   character(len=16) papi_ctr1_str
   character(len=16) papi_ctr2_str
   character(len=16) papi_ctr3_str
   character(len=16) papi_ctr4_str
   namelist /papi_inparm/ papi_ctr1_str, papi_ctr2_str,  &
                          papi_ctr3_str, papi_ctr4_str
!-----------------------------------------------------------------------
    if ( timing_initialized ) then
       return
    endif

!$OMP MASTER
    if ( present(LogUnit) ) then
       call t_setLogUnit(LogUnit)
    else
       call t_setLogUnit(def_p_logunit)
    endif

    if ( present(MasterTask) .and. present(mpicom) )then
       call mpi_comm_rank(mpicom, me, ierr)
       if (ierr .eq. MPI_ERR_COMM) then
          call shr_sys_abort('T_INITF: bad mpi communicator')
       endif
       if (me .eq. 0) then
          MasterTask2 = .true.
       else
          MasterTask2 = .false.
       endif
    else
       MasterTask2 = .true.
    end if

    if ( present(LogPrint) ) then
       LogPrint2 = LogPrint
    else
       LogPrint2 = .true.
    endif

    ! Set PERF defaults, then override with user-specified input
    call perf_defaultopts(timing_disable_out=profile_disable, &
                          perf_timer_out=profile_timer, &
                          timer_depth_limit_out=profile_depth_limit, &
                          timing_detail_limit_out=profile_detail_limit, &
                          timing_barrier_out=profile_barrier, &
                          perf_outpe_num_out = profile_outpe_num, &
                          perf_outpe_stride_out = profile_outpe_stride, &
                          perf_single_file_out=profile_single_file, &
                          perf_global_stats_out=profile_global_stats, &
                          perf_papi_enable_out=profile_papi_enable )
    if ( MasterTask2 ) then

       ! Read in the prof_inparm namelist from NLFilename if it exists

       write(p_logunit,*) '(t_initf) Read in prof_inparm namelist from: '//trim(NLFilename)
       unitn = shr_file_getUnit()

       ierr = 1
       open( unitn, file=trim(NLFilename), status='old', iostat=ierr )
       if (ierr .eq. 0) then

          ! Look for prof_inparm group name in the input file.  
          ! If found, leave the file positioned at that namelist group.
          call find_group_name(unitn, 'prof_inparm', status=ierr)

          if (ierr == 0) then  ! found prof_inparm
             read(unitn, nml=prof_inparm, iostat=ierr)  
             if (ierr /= 0) then
                call shr_sys_abort( subname//':: namelist read returns an'// &
                                    ' error condition for prof_inparm' )
             end if
          end if

          close(unitn)

       endif
       call shr_file_freeUnit( unitn )

    endif

    ! This logic assumes that there will be only one MasterTask
    ! per communicator, and that this MasterTask is process 0.
    if ( present(MasterTask) .and. present(mpicom) )then
       call shr_mpi_bcast( profile_disable,      MPICom )
       call shr_mpi_bcast( profile_barrier,      MPICom )
       call shr_mpi_bcast( profile_single_file,  MPICom )
       call shr_mpi_bcast( profile_global_stats, MPICom )
       call shr_mpi_bcast( profile_papi_enable,  MPICom )
       call shr_mpi_bcast( profile_depth_limit,  MPICom )
       call shr_mpi_bcast( profile_detail_limit, MPICom )
       call shr_mpi_bcast( profile_outpe_num,    MPICom )
       call shr_mpi_bcast( profile_outpe_stride, MPICom )
       call shr_mpi_bcast( profile_timer,        MPICom )
    end if
    call perf_setopts    (MasterTask2, LogPrint2, &
                          timing_disable_in=profile_disable, &
                          perf_timer_in=profile_timer, &
                          timer_depth_limit_in=profile_depth_limit, &
                          timing_detail_limit_in=profile_detail_limit, &
                          timing_barrier_in=profile_barrier, &
                          perf_outpe_num_in=profile_outpe_num, &
                          perf_outpe_stride_in=profile_outpe_stride, &
                          perf_single_file_in=profile_single_file, &
                          perf_global_stats_in=profile_global_stats, &
                          perf_papi_enable_in=profile_papi_enable )

    ! Set PAPI defaults, then override with user-specified input
    if (perf_papi_enable) then
       call papi_defaultopts(papi_ctr1_out=papi_ctr1_id, &
                             papi_ctr2_out=papi_ctr2_id, &
                             papi_ctr3_out=papi_ctr3_id, &
                             papi_ctr4_out=papi_ctr4_id )

       if ( MasterTask2 ) then
          papi_ctr1_str = "PAPI_NO_CTR"
          papi_ctr2_str = "PAPI_NO_CTR"
          papi_ctr3_str = "PAPI_NO_CTR"
          papi_ctr4_str = "PAPI_NO_CTR"


          ! Read in the papi_inparm namelist from NLFilename if it exists

          write(p_logunit,*) '(t_initf) Read in papi_inparm namelist from: '//trim(NLFilename)
          unitn = shr_file_getUnit()

          ierr = 1
          open( unitn, file=trim(NLFilename), status='old', iostat=ierr )
          if (ierr .eq. 0) then
             ! Look for papi_inparm group name in the input file.  
             ! If found, leave the file positioned at that namelist group.
             call find_group_name(unitn, 'papi_inparm', status=ierr)

             if (ierr == 0) then  ! found papi_inparm
                read(unitn, nml=papi_inparm, iostat=ierr)  
                if (ierr /= 0) then
                   call shr_sys_abort( subname//':: namelist read returns an'// &
                                      ' error condition for papi_inparm' )
                end if
             end if

             close(unitn)

          endif
          call shr_file_freeUnit( unitn )

          ! if enabled and nothing set, use "defaults"
          if ((papi_ctr1_str(1:11) .eq. "PAPI_NO_CTR") .and. &
              (papi_ctr2_str(1:11) .eq. "PAPI_NO_CTR") .and. &
              (papi_ctr3_str(1:11) .eq. "PAPI_NO_CTR") .and. &
              (papi_ctr4_str(1:11) .eq. "PAPI_NO_CTR")) then
!pw              papi_ctr1_str = "PAPI_TOT_CYC"
!pw              papi_ctr2_str = "PAPI_TOT_INS"
!pw              papi_ctr3_str = "PAPI_FP_OPS"
!pw              papi_ctr4_str = "PAPI_FP_INS"
              papi_ctr1_str = "PAPI_FP_OPS"
          endif

          if (papi_ctr1_str(1:11) /= "PAPI_NO_CTR") then
             ierr = gptlevent_name_to_code(trim(papi_ctr1_str), papi_ctr1_id)
          endif
          if (papi_ctr2_str(1:11) /= "PAPI_NO_CTR") then
             ierr = gptlevent_name_to_code(trim(papi_ctr2_str), papi_ctr2_id)
          endif
          if (papi_ctr3_str(1:11) /= "PAPI_NO_CTR") then
             ierr = gptlevent_name_to_code(trim(papi_ctr3_str), papi_ctr3_id)
          endif
          if (papi_ctr4_str(1:11) /= "PAPI_NO_CTR") then
             ierr = gptlevent_name_to_code(trim(papi_ctr4_str), papi_ctr4_id)
          endif

       endif
       ! This logic assumes that there will be only one MasterTask
       ! per communicator, and that this MasterTask is process 0.
       if ( present(MasterTask) .and. present(mpicom) )then
          call shr_mpi_bcast( papi_ctr1_id,    MPICom )
          call shr_mpi_bcast( papi_ctr2_id,    MPICom )
          call shr_mpi_bcast( papi_ctr3_id,    MPICom )
          call shr_mpi_bcast( papi_ctr4_id,    MPICom )
       end if

       call papi_setopts    (papi_ctr1_in=papi_ctr1_id, &
                             papi_ctr2_in=papi_ctr2_id, &
                             papi_ctr3_in=papi_ctr3_id, &
                             papi_ctr4_in=papi_ctr4_id )
    endif
!$OMP END MASTER
!$OMP BARRIER

   if (timing_disable) return

!$OMP MASTER
   !
   ! Set options and initialize timing library.  
   ! 
   ! Set timer
   if (gptlsetutr (perf_timer) < 0) call shr_sys_abort (subname//':: gptlsetutr')
   !
   ! For logical settings, 2nd arg 0 
   ! to gptlsetoption means disable, non-zero means enable
   !
   ! Turn off CPU timing (expensive)
   !
   if (gptlsetoption (gptlcpu, 0) < 0) call shr_sys_abort (subname//':: gptlsetoption')
   !
   ! Set max timer depth
   !
   if (gptlsetoption (gptldepthlimit, timer_depth_limit) < 0) &
     call shr_sys_abort (subname//':: gptlsetoption')
   !
   ! Next 2 calls only work if PAPI is enabled.  These examples enable counting
   ! of total cycles and floating point ops, respectively
   !
   if (perf_papi_enable) then
      if (papi_ctr1 /= PAPI_NULL) then
         if (gptlsetoption (papi_ctr1, 1) < 0) call shr_sys_abort (subname//':: gptlsetoption')
      endif
      if (papi_ctr2 /= PAPI_NULL) then
         if (gptlsetoption (papi_ctr2, 1) < 0) call shr_sys_abort (subname//':: gptlsetoption')
      endif
      if (papi_ctr3 /= PAPI_NULL) then
         if (gptlsetoption (papi_ctr3, 1) < 0) call shr_sys_abort (subname//':: gptlsetoption')
      endif
      if (papi_ctr4 /= PAPI_NULL) then
         if (gptlsetoption (papi_ctr4, 1) < 0) call shr_sys_abort (subname//':: gptlsetoption')
      endif
   endif
   !
   ! Initialize the timing lib.  This call must occur after all gptlsetoption
   ! calls and before all other timing lib calls.
   !
   if (gptlinitialize () < 0) call shr_sys_abort (subname//':: gptlinitialize')
   timing_initialized = .true.
!$OMP END MASTER
!$OMP BARRIER

   return
   end subroutine t_initf
!
!========================================================================
!
   subroutine t_finalizef()
!----------------------------------------------------------------------- 
! Purpose: shut down timing library
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

!$OMP MASTER
   ierr = GPTLfinalize()
   timing_initialized = .false.
!$OMP END MASTER
!$OMP BARRIER

   return
   end subroutine t_finalizef

!===============================================================================

end module perf_mod
