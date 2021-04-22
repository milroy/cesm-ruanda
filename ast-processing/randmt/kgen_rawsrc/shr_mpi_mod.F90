!===============================================================================
! SVN $Id: shr_mpi_mod.F90 59033 2014-04-11 01:55:15Z santos@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140626/shr/shr_mpi_mod.F90 $
!===============================================================================

Module shr_mpi_mod

!-------------------------------------------------------------------------------
! PURPOSE: general layer on MPI functions
!-------------------------------------------------------------------------------

   use shr_kind_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none
   private

! PUBLIC: Public interfaces

   public :: shr_mpi_chkerr
   public :: shr_mpi_send
   public :: shr_mpi_recv
   public :: shr_mpi_bcast
   public :: shr_mpi_gathScatVInit
   public :: shr_mpi_gatherV
   public :: shr_mpi_scatterV
   public :: shr_mpi_sum
   public :: shr_mpi_min
   public :: shr_mpi_max
   public :: shr_mpi_commsize
   public :: shr_mpi_commrank
   public :: shr_mpi_initialized
   public :: shr_mpi_abort
   public :: shr_mpi_barrier
   public :: shr_mpi_init
   public :: shr_mpi_finalize

   interface shr_mpi_send ; module procedure &
     shr_mpi_sendi0, &
     shr_mpi_sendi1, &
     shr_mpi_sendr0, &
     shr_mpi_sendr1, &
     shr_mpi_sendr3
   end interface
   interface shr_mpi_recv ; module procedure &
     shr_mpi_recvi0, &
     shr_mpi_recvi1, &
     shr_mpi_recvr0, &
     shr_mpi_recvr1, &
     shr_mpi_recvr3
   end interface
   interface shr_mpi_bcast ; module procedure &
     shr_mpi_bcastc0, &
     shr_mpi_bcastc1, &
     shr_mpi_bcastl0, &
     shr_mpi_bcastl1, &
     shr_mpi_bcasti0, &
     shr_mpi_bcasti1, &
     shr_mpi_bcasti2, &
     shr_mpi_bcastr0, &
     shr_mpi_bcastr1, &
     shr_mpi_bcastr2, &
     shr_mpi_bcastr3
   end interface
   interface shr_mpi_gathScatVInit ; module procedure &
     shr_mpi_gathScatVInitr1
   end interface
   interface shr_mpi_gatherv ; module procedure &
     shr_mpi_gatherVr1
   end interface
   interface shr_mpi_scatterv ; module procedure &
     shr_mpi_scatterVr1
   end interface
   interface shr_mpi_sum ; module procedure &
     shr_mpi_sumi0, &
     shr_mpi_sumi1, &
     shr_mpi_sumb0, &
     shr_mpi_sumb1, &
     shr_mpi_sumr0, &
     shr_mpi_sumr1, &
     shr_mpi_sumr2, &
     shr_mpi_sumr3
   end interface
   interface shr_mpi_min ; module procedure &
     shr_mpi_mini0, &
     shr_mpi_mini1, &
     shr_mpi_minr0, &
     shr_mpi_minr1
   end interface
   interface shr_mpi_max ; module procedure &
     shr_mpi_maxi0, &
     shr_mpi_maxi1, &
     shr_mpi_maxr0, &
     shr_mpi_maxr1
   end interface


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

!===============================================================================
CONTAINS
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
     write(s_logunit,*) trim(subName),":",lstring(1:len)
     call shr_mpi_abort(string,rcode)
   endif

END SUBROUTINE shr_mpi_chkerr

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendi0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! send value
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendi0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a single integer
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_SEND(lvec,lsize,MPI_INTEGER,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendi1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendi1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of integers
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_SEND(lvec,lsize,MPI_INTEGER,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendr0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendr0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a real scalar
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_SEND(lvec,lsize,MPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendr1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendr1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_SEND(lvec,lsize,MPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendr3(array,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real   (SHR_KIND_R8), intent(in) :: array(:,:,:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid           ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag           ! tag
   integer(SHR_KIND_IN), intent(in) :: comm          ! mpi communicator
   character(*),optional,intent(in) :: string        ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendr3) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(array)

   call MPI_SEND(array,lsize,MPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvi0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(out):: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvi0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_RECV(lvec,lsize,MPI_INTEGER,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvi1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(out):: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvi1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_RECV(lvec,lsize,MPI_INTEGER,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvr0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(out):: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvr0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_RECV(lvec,lsize,MPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvr1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(out):: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvr1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_RECV(lvec,lsize,MPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvr3(array,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real   (SHR_KIND_R8), intent(out):: array(:,:,:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid           ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag           ! tag
   integer(SHR_KIND_IN), intent(in) :: comm          ! mpi communicator
   character(*),optional,intent(in) :: string        ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvr3) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(array)

   call MPI_RECV(array,lsize,MPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcasti0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
!-------------------------------------------------------------------------------

   lsize = 1
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_INTEGER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastl0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   logical, intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastl0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a logical
!-------------------------------------------------------------------------------

   lsize = 1
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_LOGICAL,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastl0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastc0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   character(len=*), intent(inout)    :: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastc0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a character string
!-------------------------------------------------------------------------------

   lsize = len(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_CHARACTER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastc0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastc1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   character(len=*), intent(inout)    :: vec(:)   ! 1D vector
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastc1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a character string
!-------------------------------------------------------------------------------

   lsize = size(vec)*len(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_CHARACTER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastc1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastr0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a real
!-------------------------------------------------------------------------------

   lsize = 1
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec(:)   ! vector 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcasti1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of integers
!-------------------------------------------------------------------------------

   lsize = size(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_INTEGER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastl1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   logical, intent(inout):: vec(:)      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastl1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a logical
!-------------------------------------------------------------------------------

   lsize = size(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_LOGICAL,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastl1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(inout):: vec(:)   ! vector 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastr1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr2(arr,comm,string,pebcast)

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(inout):: arr(:,:) ! array, 2d 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_bcastr2) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 2d array of reals
!-------------------------------------------------------------------------------

   lsize = size(arr)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(arr,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr2

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti2(arr,comm,string,pebcast)

   IMPLICIT none

   !----- arguments -----
   integer,              intent(inout):: arr(:,:) ! array, 2d 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_bcasti2) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 2d array of integers
!-------------------------------------------------------------------------------

   lsize = size(arr)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(arr,lsize,MPI_INTEGER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti2

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr3(arr,comm,string,pebcast)

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(inout):: arr(:,:,:) ! array, 3d 
   integer(SHR_KIND_IN), intent(in)   :: comm       ! mpi communicator
   character(*),optional,intent(in)   :: string     ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_bcastr3) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 3d array of reals
!-------------------------------------------------------------------------------

   lsize = size(arr)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(arr,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_gathScatvInitr1(comm, rootid, locArr, glob1DArr, globSize, &
                                   displs, string )

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in)   :: comm          ! mpi communicator
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! MPI task to gather/scatter on
   real(SHR_KIND_R8),    intent(in)   :: locArr(:)     ! Local array of distributed data
   real(SHR_KIND_R8),    pointer      :: glob1DArr(:)  ! Global 1D array of gathered data
   integer(SHR_KIND_IN), pointer      :: globSize(:)   ! Size of each distributed piece
   integer(SHR_KIND_IN), pointer      :: displs(:)     ! Displacements for receive
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: npes          ! Number of MPI tasks
   integer(SHR_KIND_IN)               :: locSize       ! Size of local distributed data
   integer(SHR_KIND_IN), pointer      :: sendSize(:)   ! Size to send for initial gather
   integer(SHR_KIND_IN)               :: i             ! Index
   integer(SHR_KIND_IN)               :: rank          ! Rank of this MPI task
   integer(SHR_KIND_IN)               :: nSize         ! Maximum size to send
   integer(SHR_KIND_IN)               :: ierr          ! Error code
   integer(SHR_KIND_IN)               :: nSiz1D        ! Size of 1D global array
   integer(SHR_KIND_IN)               :: maxSize       ! Maximum size

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_gathScatvInitr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Setup arrays for a gatherv/scatterv operation
!-------------------------------------------------------------------------------

   locSize = size(locarr)
   call shr_mpi_commsize( comm, npes )
   call shr_mpi_commrank( comm, rank )
   allocate( globSize(npes) )
   !
   ! --- Gather the send global sizes from each MPI task -----------------------
   !
   allocate( sendSize(npes) )
   sendSize(:) = 1
   globSize(:) = 1
   call MPI_GATHER( locSize, 1, MPI_INTEGER, globSize, sendSize, &
                    MPI_INTEGER, rootid, comm, ierr )
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif
   deallocate( sendSize )
   !
   ! --- Prepare the displacement and allocate arrays -------------------------
   !
   allocate( displs(npes) )
   displs(1) = 0
   if ( rootid /= rank )then
      maxSize = 1
      globSize = 1
   else
      maxSize = maxval(globSize)
   end if
   nsiz1D  = min(maxSize,globSize(1))
   do i = 2, npes
      nSize = min(maxSize,globSize(i-1))
      displs(i) = displs(i-1) + nSize
      nsiz1D = nsiz1D + min(maxSize,globSize(i))
   end do
   allocate( glob1DArr(nsiz1D) )
   !----- Do some error checking for the root task arrays computed ----
   if ( rootid == rank )then
      if ( nsiz1D /= sum(globSize) ) &
         call shr_mpi_abort( subName//" : Error, size of global array not right" )
      if ( any(displs < 0) .or. any(displs >= nsiz1D) ) &
         call shr_mpi_abort( subName//" : Error, displacement array not right" )
      if ( (displs(npes)+globSize(npes)) /= nsiz1D ) &
         call shr_mpi_abort( subName//" : Error, displacement array values too big" )
   end if

END SUBROUTINE shr_mpi_gathScatvInitr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_gathervr1(locarr, locSize, glob1DArr, globSize, displs, rootid, &
                             comm, string )

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(in)   :: locArr(:)     ! Local array
   real(SHR_KIND_R8),    intent(inout):: glob1DArr(:)  ! Global 1D array to receive in on
   integer(SHR_KIND_IN), intent(in)   :: locSize       ! Number to send this PE
   integer(SHR_KIND_IN), intent(in)   :: globSize(:)   ! Number to receive each PE
   integer(SHR_KIND_IN), intent(in)   :: displs(:)     ! Displacements for receive
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! MPI task to gather on
   integer(SHR_KIND_IN), intent(in)   :: comm          ! mpi communicator
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr          ! Error code

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_gathervr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Gather a 1D array of reals
!-------------------------------------------------------------------------------

   call MPI_GATHERV( locarr, locSize, MPI_REAL8, glob1Darr, globSize, displs, &
                     MPI_REAL8, rootid, comm, ierr )
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_gathervr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_scattervr1(locarr, locSize, glob1Darr, globSize, displs, rootid, &
                              comm, string )

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(out)  :: locarr(:)     ! Local array
   real(SHR_KIND_R8),    intent(in)   :: glob1Darr(:)  ! Global 1D array to send from
   integer(SHR_KIND_IN), intent(in)   :: locSize       ! Number to receive this PE
   integer(SHR_KIND_IN), intent(in)   :: globSize(:)   ! Number to send to each PE
   integer(SHR_KIND_IN), intent(in)   :: displs(:)     ! Displacements for send
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! MPI task to scatter on
   integer(SHR_KIND_IN), intent(in)   :: comm          ! mpi communicator
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr          ! Error code

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_scattervr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Scatter a 1D array of reals
!-------------------------------------------------------------------------------


   call MPI_SCATTERV( glob1Darr, globSize, displs, MPI_REAL8, locarr, locSize, &
                      MPI_REAL8, rootid, comm, ierr )
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_scattervr1


!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumi0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumi0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumi1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumi1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumb0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_I8), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_I8), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumb0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumb0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumb1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_I8), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_I8), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumb1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumb1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr2(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:,:)! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:,:)! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr2) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr2

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr3(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:,:,:) ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:,:,:) ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr3) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_mini0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_mini0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_mini0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_mini1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_mini1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_mini1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_minr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_minr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_minr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_minr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_minr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_minr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxi0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxi0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxi1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxi1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_commsize(comm,size,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   integer,intent(out)                :: size
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_commsize) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI commsize
!-------------------------------------------------------------------------------

   call MPI_COMM_SIZE(comm,size,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_commsize

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_commrank(comm,rank,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   integer,intent(out)                :: rank
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_commrank) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI commrank
!-------------------------------------------------------------------------------

   call MPI_COMM_RANK(comm,rank,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_commrank

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_initialized(flag,string)

   IMPLICIT none

   !----- arguments ---
   logical,intent(out)                :: flag
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_initialized) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI initialized
!-------------------------------------------------------------------------------

   call MPI_INITIALIZED(flag,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_initialized

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_abort(string,rcode)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message
   integer,optional,intent(in)        :: rcode    ! optional code

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_abort) '
   integer(SHR_KIND_IN)               :: ierr
   integer                            :: rc       ! return code

!-------------------------------------------------------------------------------
! PURPOSE: MPI abort
!-------------------------------------------------------------------------------

   if ( present(string) .and. present(rcode) ) then
      write(s_logunit,*) trim(subName),":",trim(string),rcode
   endif
   if ( present(rcode) )then
      rc = rcode
   else
      rc = 1001
   end if
   call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)

END SUBROUTINE shr_mpi_abort

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_barrier(comm,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_barrier) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI barrier
!-------------------------------------------------------------------------------

   call MPI_BARRIER(comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_barrier

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_init(string)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_init) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI init
!-------------------------------------------------------------------------------

   call MPI_INIT(ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_init

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_finalize(string)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_finalize) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI finalize
!-------------------------------------------------------------------------------

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_finalize

!===============================================================================
!===============================================================================

END MODULE shr_mpi_mod
